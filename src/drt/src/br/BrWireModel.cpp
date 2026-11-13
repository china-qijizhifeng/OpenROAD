// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrWireModel.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <utility>
#include <vector>

#include "br/BrGeometry.h"

namespace drt::br {

SpacingTable SpacingTable::constant(const Coord spacing)
{
  SpacingTable t;
  t.values_ = {std::max(0, spacing)};
  return t;
}

void SpacingTable::close()
{
  const int rows = numRows();
  const int cols = numColumns();
  auto at = [&](int a, int b, int c) -> Coord& {
    return values_[(static_cast<size_t>(a) * rows + b) * cols + c];
  };
  for (int a = 0; a < rows; a++) {
    for (int b = a + 1; b < rows; b++) {
      for (int c = 0; c < cols; c++) {
        const Coord v = std::max(at(a, b, c), at(b, a, c));
        at(a, b, c) = v;
        at(b, a, c) = v;
      }
    }
  }
  for (int a = 0; a < rows; a++) {
    for (int b = 0; b < rows; b++) {
      for (int c = 0; c < cols; c++) {
        Coord v = at(a, b, c);
        if (a > 0) {
          v = std::max(v, at(a - 1, b, c));
        }
        if (b > 0) {
          v = std::max(v, at(a, b - 1, c));
        }
        if (c > 0) {
          v = std::max(v, at(a, b, c - 1));
        }
        at(a, b, c) = v;
      }
    }
  }
}

int SpacingTable::widthRow(const Coord width) const
{
  const auto it = std::ranges::upper_bound(width_breaks_, width);
  return std::max(0, static_cast<int>(it - width_breaks_.begin()) - 1);
}

int SpacingTable::prlColumn(const Coord prl) const
{
  const auto it = std::ranges::upper_bound(prl_breaks_, prl);
  return std::max(0, static_cast<int>(it - prl_breaks_.begin()) - 1);
}

Coord SpacingTable::maxValue() const
{
  return values_.empty() ? 0 : *std::ranges::max_element(values_);
}

Coord SpacingTable::minValue() const
{
  return values_.empty() ? 0 : *std::ranges::min_element(values_);
}

int LayerRules::cutClassOf(const Coord dx, const Coord dy) const
{
  const Coord w = std::min(dx, dy);
  const Coord l = std::max(dx, dy);
  // drt's frLayer::getCutClassIdx: the last class with these dimensions.
  int found = 0;
  for (size_t k = 0; k < cut_classes.size(); k++) {
    if (cut_classes[k].width == w && cut_classes[k].length == l) {
      found = static_cast<int>(k) + 1;
    }
  }
  return found;
}

Coord LayerRules::cutClassSpacing(const int a, const int b) const
{
  const size_t n = cut_classes.size() + 1;
  if (cut_class_spacing.size() != n * n) {
    return 0;
  }
  return cut_class_spacing[static_cast<size_t>(a) * n + b];
}

Coord LayerRules::sideReach() const
{
  Coord r = corner_spacing;
  for (const SideRule& rule : side_rules) {
    r = std::max({r, rule.forward, rule.backward, rule.side});
  }
  // Regions of stored pieces reach one unit beyond a cut border (see
  // appendSideTests).
  return sideRules() ? r + 1 : 0;
}

Coord LayerRules::reach() const
{
  Coord r = std::max(spacing.maxValue(), projection_spacing);
  for (const Coord s : ndr_spacing) {
    r = std::max(r, s);
  }
  if (!via_layer && line_end.active()) {
    r = std::max({r, line_end.space, line_end.within});
  }
  for (const Coord s : cut_class_spacing) {
    r = std::max(r, s);
  }
  return std::max(r, sideReach());
}

Coord LayerRules::minimumSpacing() const
{
  return spacing.minValue();
}

size_t ShapeClassHash::operator()(const ShapeClass& c) const
{
  uint64_t h = static_cast<uint16_t>(c.layer);
  h = h * 31 + static_cast<uint64_t>(c.kind);
  h = h * 31 + c.width_row;
  h = h * 257 + c.run_column_x;
  h = h * 257 + c.run_column_y;
  h = h * 2 + (c.line_end ? 1 : 0);
  h = h * 257 + c.ndr;
  h = h * 257 + c.cut_class;
  h = h * 17 + c.sides_real;
  h = h * 65537 + c.eol_sides;
  return std::hash<uint64_t>()(h);
}

void ShapeClassTable::setLayers(std::vector<LayerRules> layers)
{
  static std::atomic<uint64_t> generations{0};
  generation_ = ++generations;
  layers_ = std::move(layers);
  const size_t n = layers_.size();
  rule_reach_.clear();
  narrow_rule_reach_.clear();
  wide_row_.clear();
  wide_ndr_.assign(n, {});
  class_reach_ = std::make_unique<std::atomic<Coord>[]>(n);
  narrow_reach_ = std::make_unique<std::atomic<Coord>[]>(n);
  layer_classes_.assign(n, {});
  for (size_t l = 0; l < n; l++) {
    const LayerRules& rules = layers_[l];
    // An overlap rule makes touching shapes interact at distance 1.
    rule_reach_.push_back(std::max<Coord>(rules.reach(), 1));
    class_reach_[l].store(0, std::memory_order_relaxed);
    narrow_reach_[l].store(0, std::memory_order_relaxed);
    // The table is nondecreasing: row r's largest value bounds every pair
    // of rows up to r.
    const SpacingTable& t = rules.spacing;
    const int last = std::max(t.numColumns() - 1, 0);
    int row = std::max(t.numRows(), 1);
    const int64_t base = t.numRows() > 0 && t.numColumns() > 0
                             ? t.at(0, 0, last)
                             : t.maxValue();
    for (const Coord s : rules.ndr_spacing) {
      wide_ndr_[l].push_back(wide_percent_ > 0
                             && int64_t{s} * 100 >= base * wide_percent_);
    }
    if (wide_percent_ > 0 && t.numRows() > 0 && t.numColumns() > 0) {
      for (int r = 1; r < t.numRows(); r++) {
        if (int64_t{t.at(r, r, last)} * 100 >= base * wide_percent_) {
          row = r;
          break;
        }
      }
    }
    wide_row_.push_back(row);
    // Narrow classes: width rows below `row` (the largest value of row
    // row - 1 bounds every pair of them), nondefault rules that are not
    // wide, cuts and cut projections; line-end, cut class and side rules
    // apply to them as to all classes.
    Coord narrow = rules.projection_spacing;
    if (t.numRows() > 0 && t.numColumns() > 0) {
      narrow = std::max(narrow, t.at(row - 1, row - 1, last));
    } else {
      narrow = std::max(narrow, t.maxValue());
    }
    for (size_t k = 0; k < rules.ndr_spacing.size(); k++) {
      if (!wide_ndr_[l][k]) {
        narrow = std::max(narrow, rules.ndr_spacing[k]);
      }
    }
    if (!rules.via_layer && rules.line_end.active()) {
      narrow = std::max({narrow, rules.line_end.space, rules.line_end.within});
    }
    for (const Coord s : rules.cut_class_spacing) {
      narrow = std::max(narrow, s);
    }
    narrow_rule_reach_.push_back(
        std::max<Coord>(std::max(narrow, rules.sideReach()), 1));
  }
  index_.clear();
  classes_.clear();
}

Coord ShapeClassTable::distance(const ClassId a, const ClassId b) const
{
  const Requirement r = requirement(a, b);
  Coord d = std::max({r.facing_x, r.facing_y, r.diagonal});
  if (r.line_end) {
    const LineEndRule& le = layers_[get(a).layer].line_end;
    d = std::max({d, le.space, le.within});
  }
  if (r.overlap) {
    d = std::max<Coord>(d, 1);
  }
  if (r.sides) {
    d = std::max(d, layers_[get(a).layer].sideReach());
  }
  return d;
}

ClassId ShapeClassTable::intern(const ShapeClass& c)
{
  std::lock_guard lock(mutex_);
  const auto it = index_.find(c);
  if (it != index_.end()) {
    return it->second;
  }
  ShapeClass copy = c;
  const size_t id = classes_.append(std::move(copy));
  if (id > 0xffff) {
    throw std::length_error("too many shape classes");
  }
  const auto cls = static_cast<ClassId>(id);
  index_.emplace(c, cls);
  // Grow the interaction distances of the layer by the pairs of the new
  // class with the interned ones (itself included).
  if (c.layer >= 0 && static_cast<size_t>(c.layer) < layers_.size()) {
    std::vector<ClassId>& list = layer_classes_[c.layer];
    list.push_back(cls);
    const bool narrow = !wide(cls);
    Coord all = class_reach_[c.layer].load(std::memory_order_relaxed);
    Coord near = narrow_reach_[c.layer].load(std::memory_order_relaxed);
    for (const ClassId other : list) {
      const Coord d = distance(cls, other);
      all = std::max(all, d);
      if (narrow && !wide(other)) {
        near = std::max(near, d);
      }
    }
    class_reach_[c.layer].store(all, std::memory_order_relaxed);
    narrow_reach_[c.layer].store(near, std::memory_order_relaxed);
  }
  return cls;
}

ShapeClass ShapeClassTable::describe(const int layer,
                                     const Rect& rect,
                                     const ShapeRole role,
                                     const int ndr) const
{
  const LayerRules& rules = layers_[layer];
  ShapeClass c;
  c.layer = static_cast<int16_t>(layer);
  if (role == ShapeRole::kCutProjection) {
    c.kind = ShapeClassKind::kCutProjection;
  } else if (rules.via_layer) {
    c.kind = ShapeClassKind::kCut;
  } else {
    c.kind = ShapeClassKind::kMetal;
  }
  const SpacingTable& table = rules.spacing;
  c.width_row = static_cast<uint8_t>(table.widthRow(rect.width()));
  c.run_column_x = static_cast<uint8_t>(table.prlColumn(rect.dx()));
  c.run_column_y = static_cast<uint8_t>(table.prlColumn(rect.dy()));
  const bool wire
      = role == ShapeRole::kPreferredWire || role == ShapeRole::kRealWire;
  if (c.kind == ShapeClassKind::kMetal && wire) {
    (rules.horizontal ? c.run_column_x : c.run_column_y) = kUnboundedRun;
  }
  if (c.kind == ShapeClassKind::kMetal && rules.line_end.active()
      && (wire || role == ShapeRole::kViaPad)) {
    // The line-end edge runs across the preferred direction.
    const Coord edge = rules.horizontal ? rect.dy() : rect.dx();
    c.line_end = edge < rules.line_end.max_edge;
  }
  if (c.kind == ShapeClassKind::kMetal && rules.ndrSpacing(ndr) > 0) {
    c.ndr = static_cast<uint8_t>(ndr);
  }
  if (c.kind == ShapeClassKind::kCut && !rules.cut_classes.empty()) {
    c.cut_class = static_cast<uint8_t>(rules.cutClassOf(rect.dx(), rect.dy()));
  }
  if (c.kind == ShapeClassKind::kMetal && rules.sideRules()) {
    // Lef58Rules.md: sides along the preferred direction run across it.
    const uint8_t ends = rules.horizontal ? 0x3 : 0xc;
    c.sides_real = role == ShapeRole::kPreferredWire && rules.wire_ends_covered
                       ? static_cast<uint8_t>(0xf & ~ends)
                       : uint8_t{0xf};
    for (size_t k = 0; k < rules.eol_widths.size(); k++) {
      const Coord w = rules.eol_widths[k];
      // x sides (xlo, xhi) are dy long, y sides dx.
      uint16_t bits = (rect.dy() < w ? 0x3 : 0) | (rect.dx() < w ? 0xc : 0);
      if (role == ShapeRole::kPreferredWire) {
        // A routed wire's flanks are as long as the wire (its model is the
        // stick figure's square); kRealWire shapes are actual rectangles.
        bits &= ends;
      }
      bits &= c.sides_real;
      c.eol_sides |= static_cast<uint16_t>(bits << (4 * k));
    }
  }
  return c;
}

Requirement ShapeClassTable::requirement(const ClassId a, const ClassId b) const
{
  const ShapeClass& ca = get(a);
  const ShapeClass& cb = get(b);
  Requirement r;
  if (ca.layer != cb.layer) {
    return r;
  }
  const LayerRules& rules = layers_[ca.layer];
  const SpacingTable& t = rules.spacing;
  const bool metal_a = ca.kind == ShapeClassKind::kMetal;
  const bool metal_b = cb.kind == ShapeClassKind::kMetal;
  if (metal_a && metal_b) {
    const int last = t.numColumns() - 1;
    auto column = [&](uint8_t x, uint8_t y) {
      return std::min<int>(std::min(x, y), last);
    };
    r.facing_y = t.at(
        ca.width_row, cb.width_row, column(ca.run_column_x, cb.run_column_x));
    r.facing_x = t.at(
        ca.width_row, cb.width_row, column(ca.run_column_y, cb.run_column_y));
    r.diagonal = t.at(ca.width_row, cb.width_row, 0);
    // NDR spacing: max with the table value, Euclidean (drt's
    // checkMetalSpacing_prl with checkNDRs).
    const Coord ndr
        = std::max(rules.ndrSpacing(ca.ndr), rules.ndrSpacing(cb.ndr));
    r.facing_x = std::max(r.facing_x, ndr);
    r.facing_y = std::max(r.facing_y, ndr);
    r.diagonal = std::max(r.diagonal, ndr);
    r.overlap = true;
    r.line_end = rules.line_end.active();
    r.sides = rules.sideRules() && (ca.sides_real != 0 || cb.sides_real != 0);
    return r;
  }
  if (metal_a || metal_b) {
    return r;
  }
  const bool cut_a = ca.kind == ShapeClassKind::kCut;
  const bool cut_b = cb.kind == ShapeClassKind::kCut;
  if (cut_a && cut_b) {
    const Coord s = std::max(t.at(ca.width_row, cb.width_row, 0),
                             rules.cutClassSpacing(ca.cut_class, cb.cut_class));
    r.facing_x = r.facing_y = r.diagonal = s;
    r.overlap = true;
    return r;
  }
  if (cut_a || cut_b) {
    const Coord s = rules.projection_spacing;
    r.facing_x = r.facing_y = r.diagonal = s;
    r.overlap = s > 0;
  }
  // Two projections: no rule between them.
  return r;
}

int lineEndExtensions(const ShapeClassTable& classes,
                      const PlacedShape& shape,
                      Rect out[2])
{
  const ShapeClass& c = classes.get(shape.cls);
  if (!c.line_end) {
    return 0;
  }
  const LayerRules& rules = classes.layer(c.layer);
  const LineEndRule& le = rules.line_end;
  if (!le.active()) {
    return 0;
  }
  const Rect& r = shape.rect;
  int n = 0;
  if (rules.horizontal) {
    if (shape.line_end_ends & 1) {
      out[n++]
          = {r.xlo - le.space, r.ylo - le.within, r.xlo, r.yhi + le.within};
    }
    if (shape.line_end_ends & 2) {
      out[n++]
          = {r.xhi, r.ylo - le.within, r.xhi + le.space, r.yhi + le.within};
    }
  } else {
    if (shape.line_end_ends & 1) {
      out[n++]
          = {r.xlo - le.within, r.ylo - le.space, r.xhi + le.within, r.ylo};
    }
    if (shape.line_end_ends & 2) {
      out[n++]
          = {r.xlo - le.within, r.yhi, r.xhi + le.within, r.yhi + le.space};
    }
  }
  return n;
}

namespace {

bool realSide(const ShapeClass& c, const PlacedShape& s, const int side)
{
  return ((c.sides_real >> side) & 1) != 0
         && (s.line_end_ends & cutSideBit(side)) == 0;
}

Rect pointRect(const Coord x, const Coord y)
{
  return {x, y, x, y};
}

// Side `side` of r: position on the normal axis, outward direction and
// lateral range.
struct SideGeometry
{
  bool normal_x;
  Coord pos;
  int dir;
  Coord lo;
  Coord hi;
};

SideGeometry sideGeometry(const Rect& r, const int side)
{
  switch (side) {
    case kSideXlo:
      return {true, r.xlo, -1, r.ylo, r.yhi};
    case kSideXhi:
      return {true, r.xhi, 1, r.ylo, r.yhi};
    case kSideYlo:
      return {false, r.ylo, -1, r.xlo, r.xhi};
    default:
      return {false, r.yhi, 1, r.xlo, r.xhi};
  }
}

Rect fromAxes(const bool normal_x,
              const Coord n_lo,
              const Coord n_hi,
              const Coord l_lo,
              const Coord l_hi)
{
  return normal_x ? Rect{n_lo, l_lo, n_hi, l_hi} : Rect{l_lo, n_lo, l_hi, n_hi};
}

}  // namespace

void appendSideTests(const ShapeClassTable& classes,
                     const PlacedShape& a,
                     const PlacedShape& b,
                     std::vector<SideTest>& out)
{
  const ShapeClass& ca = classes.get(a.cls);
  const ShapeClass& cb = classes.get(b.cls);
  if (ca.layer != cb.layer) {
    return;
  }
  const LayerRules& rules = classes.layer(ca.layer);
  for (int pass = 0; pass < 2; pass++) {
    const PlacedShape& x = pass == 0 ? a : b;
    const PlacedShape& y = pass == 0 ? b : a;
    const ShapeClass& cx = pass == 0 ? ca : cb;
    const ShapeClass& cy = pass == 0 ? cb : ca;
    auto push
        = [&](const Rect& from_x, const Rect& from_y, const uint32_t bit) {
            if (pass == 0) {
              out.push_back({from_x, from_y, bit});
            } else {
              out.push_back({from_y, from_x, bit});
            }
          };
    for (const SideRule& rule : rules.side_rules) {
      const int eol = (cx.eol_sides >> (4 * rule.eol)) & 0xf;
      if (eol == 0) {
        continue;
      }
      for (int s = 0; s < 4; s++) {
        if (((eol >> s) & 1) == 0 || !realSide(cx, x, s)) {
          continue;
        }
        const SideGeometry e = sideGeometry(x.rect, s);
        // The lateral ends of a stored piece's side may be cell borders:
        // the region then reaches one unit across, so that the regions of
        // the pieces cover the open region of the whole side (a corner on
        // the border lies strictly inside one of them).
        const int lat_lo_side = e.normal_x ? kSideYlo : kSideXlo;
        const int lat_hi_side = e.normal_x ? kSideYhi : kSideXhi;
        const Coord lat_lo
            = e.lo - ((x.line_end_ends & cutSideBit(lat_lo_side)) != 0 ? 1 : 0);
        const Coord lat_hi
            = e.hi + ((x.line_end_ends & cutSideBit(lat_hi_side)) != 0 ? 1 : 0);
        switch (rule.kind) {
          case SideRuleKind::kRegionOverlap:
          case SideRuleKind::kRegionCorner: {
            const Coord n_lo
                = e.dir > 0 ? e.pos - rule.backward : e.pos - rule.forward;
            const Coord n_hi
                = e.dir > 0 ? e.pos + rule.forward : e.pos + rule.backward;
            const Rect q = fromAxes(
                e.normal_x, n_lo, n_hi, lat_lo - rule.side, lat_hi + rule.side);
            if (rule.kind == SideRuleKind::kRegionOverlap) {
              push(q, y.rect, kViolRegionOverlap);
              break;
            }
            for (const int sx : {kSideXlo, kSideXhi}) {
              for (const int sy : {kSideYlo, kSideYhi}) {
                if (realSide(cy, y, sx) && realSide(cy, y, sy)) {
                  push(q,
                       pointRect(sx == kSideXlo ? y.rect.xlo : y.rect.xhi,
                                 sy == kSideYlo ? y.rect.ylo : y.rect.yhi),
                       kViolRegionCorner);
                }
              }
            }
            break;
          }
          case SideRuleKind::kFacing:
          case SideRuleKind::kFacingEol: {
            const int opposite = s ^ 1;
            if (!realSide(cy, y, opposite)) {
              break;
            }
            if (rule.kind == SideRuleKind::kFacingEol
                && ((cy.eol_sides >> (4 * rule.eol + opposite)) & 1) == 0) {
              break;
            }
            // Strictly in front of the side, closer than `forward`.
            const Coord n_lo = e.dir > 0 ? e.pos : e.pos - rule.forward;
            const Coord n_hi = e.dir > 0 ? e.pos + rule.forward : e.pos;
            const Rect q = fromAxes(
                e.normal_x, n_lo, n_hi, lat_lo - rule.side, lat_hi + rule.side);
            const SideGeometry f = sideGeometry(y.rect, opposite);
            push(q,
                 fromAxes(f.normal_x, f.pos, f.pos, f.lo, f.hi),
                 rule.kind == SideRuleKind::kFacing ? kViolFacing
                                                    : kViolFacingEol);
            break;
          }
        }
      }
    }
  }
  // LEF58_CORNERSPACING: a convex corner of a and the facing corner of b
  // beyond it (b in its quadrant), both x and y gaps below the spacing.
  const Coord cs = rules.corner_spacing;
  if (cs <= 0) {
    return;
  }
  for (const bool ox : {false, true}) {
    for (const bool oy : {false, true}) {
      const int ax = ox ? kSideXhi : kSideXlo;
      const int ay = oy ? kSideYhi : kSideYlo;
      if (!realSide(ca, a, ax) || !realSide(ca, a, ay)
          || !realSide(cb, b, ax ^ 1) || !realSide(cb, b, ay ^ 1)) {
        continue;
      }
      const Coord bx = ox ? b.rect.xlo : b.rect.xhi;
      const Coord by = oy ? b.rect.ylo : b.rect.yhi;
      // 0 <= gap < cs on both axes, as strict bounds of a's corner.
      const Rect box{ox ? bx - cs : bx - 1,
                     oy ? by - cs : by - 1,
                     ox ? bx + 1 : bx + cs,
                     oy ? by + 1 : by + cs};
      out.push_back({pointRect(ox ? a.rect.xhi : a.rect.xlo,
                               oy ? a.rect.yhi : a.rect.ylo),
                     box,
                     kViolCorner});
    }
  }
}

namespace {

// L-infinity gap of two rectangles.
int64_t gapOf(const Rect& r, const Rect& o)
{
  return std::max<int64_t>({0,
                            static_cast<int64_t>(o.xlo) - r.xhi,
                            static_cast<int64_t>(r.xlo) - o.xhi,
                            static_cast<int64_t>(o.ylo) - r.yhi,
                            static_cast<int64_t>(r.ylo) - o.yhi});
}

uint32_t sideMask(const ShapeClassTable& classes,
                  const PlacedShape& a,
                  const PlacedShape& b,
                  const bool first_only)
{
  const LayerRules& rules = classes.layer(classes.get(a.cls).layer);
  if (gapOf(a.rect, b.rect) >= rules.sideReach()) {
    return 0;
  }
  thread_local std::vector<SideTest> tests;
  tests.clear();
  appendSideTests(classes, a, b, tests);
  uint32_t mask = 0;
  for (const SideTest& t : tests) {
    if (t.from_a.overlaps(t.from_b)) {
      mask |= t.rule;
      if (first_only) {
        break;
      }
    }
  }
  return mask;
}

}  // namespace

bool sideViolation(const ShapeClassTable& classes,
                   const PlacedShape& a,
                   const PlacedShape& b)
{
  return sideMask(classes, a, b, true) != 0;
}

bool shapesViolate(const ShapeClassTable& classes,
                   const PlacedShape& a,
                   const PlacedShape& b)
{
  const Requirement req = classes.requirement(a.cls, b.cls);
  if (req.sides && sideViolation(classes, a, b)) {
    return true;
  }
  const Rect& r = a.rect;
  const Rect& o = b.rect;
  const int64_t ox
      = static_cast<int64_t>(std::min(r.xhi, o.xhi)) - std::max(r.xlo, o.xlo);
  const int64_t oy
      = static_cast<int64_t>(std::min(r.yhi, o.yhi)) - std::max(r.ylo, o.ylo);
  const int64_t gx = std::max<int64_t>({0,
                                        static_cast<int64_t>(o.xlo) - r.xhi,
                                        static_cast<int64_t>(r.xlo) - o.xhi});
  const int64_t gy = std::max<int64_t>({0,
                                        static_cast<int64_t>(o.ylo) - r.yhi,
                                        static_cast<int64_t>(r.ylo) - o.yhi});
  if (ox > 0 && oy > 0) {
    if (req.overlap) {
      return true;
    }
  } else if (ox > 0) {
    if (gy < req.facing_y) {
      return true;
    }
  } else if (oy > 0) {
    if (gx < req.facing_x) {
      return true;
    }
  } else if (gx * gx + gy * gy
             < static_cast<int64_t>(req.diagonal) * req.diagonal) {
    return true;
  }
  if (!req.line_end) {
    return false;
  }
  Rect ext[2];
  for (int k = lineEndExtensions(classes, a, ext) - 1; k >= 0; k--) {
    if (ext[k].overlaps(o)) {
      return true;
    }
  }
  for (int k = lineEndExtensions(classes, b, ext) - 1; k >= 0; k--) {
    if (ext[k].overlaps(r)) {
      return true;
    }
  }
  return false;
}

uint32_t violationMask(const ShapeClassTable& classes,
                       const PlacedShape& a,
                       const PlacedShape& b)
{
  uint32_t mask = 0;
  const Requirement req = classes.requirement(a.cls, b.cls);
  Requirement spacing = req;
  spacing.line_end = false;
  // Spacing alone: the pairwise rule with the line-end rule switched off.
  {
    const Rect& r = a.rect;
    const Rect& o = b.rect;
    const int64_t ox
        = static_cast<int64_t>(std::min(r.xhi, o.xhi)) - std::max(r.xlo, o.xlo);
    const int64_t oy
        = static_cast<int64_t>(std::min(r.yhi, o.yhi)) - std::max(r.ylo, o.ylo);
    const int64_t gx = std::max<int64_t>({0,
                                          static_cast<int64_t>(o.xlo) - r.xhi,
                                          static_cast<int64_t>(r.xlo) - o.xhi});
    const int64_t gy = std::max<int64_t>({0,
                                          static_cast<int64_t>(o.ylo) - r.yhi,
                                          static_cast<int64_t>(r.ylo) - o.yhi});
    bool v = false;
    if (ox > 0 && oy > 0) {
      v = spacing.overlap;
    } else if (ox > 0) {
      v = gy < spacing.facing_y;
    } else if (oy > 0) {
      v = gx < spacing.facing_x;
    } else {
      v = gx * gx + gy * gy
          < static_cast<int64_t>(spacing.diagonal) * spacing.diagonal;
    }
    if (v) {
      mask |= kViolSpacing;
    }
  }
  if (req.sides) {
    mask |= sideMask(classes, a, b, false);
  }
  if (req.line_end) {
    Rect ext[2];
    bool v = false;
    for (int k = lineEndExtensions(classes, a, ext) - 1; k >= 0; k--) {
      v = v || ext[k].overlaps(b.rect);
    }
    for (int k = lineEndExtensions(classes, b, ext) - 1; k >= 0; k--) {
      v = v || ext[k].overlaps(a.rect);
    }
    if (v) {
      mask |= kViolLineEnd;
    }
  }
  return mask;
}

WireType makeWireType(const WireTypeGeometry& geometry,
                      ShapeClassTable& classes)
{
  WireType t;
  t.name = geometry.name;
  const int num_z = static_cast<int>(geometry.width.size());
  t.preferred.resize(num_z);
  t.non_preferred.resize(num_z);
  for (int z = 0; z < num_z; z++) {
    const Coord w = geometry.width[z];
    if (w <= 0) {
      continue;
    }
    const Coord half = w / 2;
    const Rect square{-half, -half, w - half, w - half};
    t.preferred[z] = {true,
                      square,
                      classes.classify(wiringGridLayer(z),
                                       square,
                                       ShapeRole::kPreferredWire,
                                       geometry.ndr)};
    t.non_preferred[z] = {true,
                          square,
                          classes.classify(wiringGridLayer(z),
                                           square,
                                           ShapeRole::kNonPreferredWire,
                                           geometry.ndr)};
  }
  t.via.resize(geometry.via.size());
  for (size_t z = 0; z < geometry.via.size(); z++) {
    const WireTypeGeometry::Via& g = geometry.via[z];
    if (!g.valid || static_cast<int>(z) + 1 >= num_z) {
      continue;
    }
    const int zi = static_cast<int>(z);
    ViaModel& m = t.via[z];
    m.valid = true;
    m.bottom = {
        true,
        g.bottom,
        classes.classify(
            wiringGridLayer(zi), g.bottom, ShapeRole::kViaPad, geometry.ndr)};
    m.cut = {true,
             g.cut,
             classes.classify(viaGridLayer(zi), g.cut, ShapeRole::kCut)};
    m.top = {
        true,
        g.top,
        classes.classify(
            wiringGridLayer(zi + 1), g.top, ShapeRole::kViaPad, geometry.ndr)};
    const int above = viaGridLayer(zi + 1);
    if (above < classes.numLayers()
        && classes.layer(above).projection_spacing > 0) {
      m.projection
          = {true,
             g.cut,
             classes.classify(above, g.cut, ShapeRole::kCutProjection)};
    }
  }
  return t;
}

}  // namespace drt::br
