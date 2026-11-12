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

Coord LayerRules::reach() const
{
  Coord r = std::max(spacing.maxValue(), projection_spacing);
  for (const Coord s : ndr_spacing) {
    r = std::max(r, s);
  }
  if (!via_layer && line_end.active()) {
    r = std::max({r, line_end.space, line_end.within});
  }
  return r;
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
  return std::hash<uint64_t>()(h);
}

void ShapeClassTable::setLayers(std::vector<LayerRules> layers)
{
  static std::atomic<uint64_t> generations{0};
  generation_ = ++generations;
  layers_ = std::move(layers);
  const size_t n = layers_.size();
  rule_reach_.clear();
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
  if (c.kind == ShapeClassKind::kMetal && role == ShapeRole::kPreferredWire) {
    (rules.horizontal ? c.run_column_x : c.run_column_y) = kUnboundedRun;
  }
  if (c.kind == ShapeClassKind::kMetal && rules.line_end.active()
      && (role == ShapeRole::kPreferredWire || role == ShapeRole::kViaPad)) {
    // The line-end edge runs across the preferred direction.
    const Coord edge = rules.horizontal ? rect.dy() : rect.dx();
    c.line_end = edge < rules.line_end.max_edge;
  }
  if (c.kind == ShapeClassKind::kMetal && rules.ndrSpacing(ndr) > 0) {
    c.ndr = static_cast<uint8_t>(ndr);
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
    return r;
  }
  if (metal_a || metal_b) {
    return r;
  }
  const bool cut_a = ca.kind == ShapeClassKind::kCut;
  const bool cut_b = cb.kind == ShapeClassKind::kCut;
  if (cut_a && cut_b) {
    const Coord s = t.at(ca.width_row, cb.width_row, 0);
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

bool shapesViolate(const ShapeClassTable& classes,
                   const PlacedShape& a,
                   const PlacedShape& b)
{
  const Requirement req = classes.requirement(a.cls, b.cls);
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
