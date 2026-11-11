// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrDistanceRules.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrWireModel.h"

namespace drt::br {

namespace {

constexpr int64_t kFar = int64_t{1} << 40;

// Placements u (integers, inclusive) along a line at which a model violates
// a piece; level and owner of the piece.
struct Span
{
  int64_t lo;
  int64_t hi;
  int group;
  RipupLevel level;
  int32_t owner;
  bool removable;
};

// Line coordinates: x = u (the moving axis), y = v (fixed). `m` is the
// model at u = 0 with absolute v; `p` the piece. Appends the spans of u at
// which the spacing rule of Requirement is violated (Sec. 3.1: l2
// distance below the required one, by facing direction).
void appendSpacingSpans(const Rect& m,
                        const Rect& p,
                        const Coord s_u,
                        const Coord s_v,
                        const Coord s_d,
                        const bool overlap,
                        const Span& tmpl,
                        std::vector<Span>& out)
{
  // u-projections overlap positively for u in [a0, a1].
  const int64_t a0 = static_cast<int64_t>(p.xlo) - m.xhi + 1;
  const int64_t a1 = static_cast<int64_t>(p.xhi) - m.xlo - 1;
  if (a0 > a1) {
    return;
  }
  auto push = [&](int64_t lo, int64_t hi) {
    if (lo <= hi) {
      Span s = tmpl;
      s.lo = lo;
      s.hi = hi;
      out.push_back(s);
    }
  };
  const int64_t ov
      = static_cast<int64_t>(std::min(m.yhi, p.yhi)) - std::max(m.ylo, p.ylo);
  if (ov > 0) {
    if (overlap) {
      push(a0 - s_u, a1 + s_u);
    } else if (s_u > 0) {
      push(a0 - s_u, a0 - 1);
      push(a1 + 1, a1 + s_u);
    }
    return;
  }
  const int64_t gv = std::max(static_cast<int64_t>(p.ylo) - m.yhi,
                              static_cast<int64_t>(m.ylo) - p.yhi);
  const bool mid = gv < s_v;
  const int64_t g
      = largestBelowSquare(static_cast<int64_t>(s_d) * s_d - gv * gv);
  if (mid && g >= 0) {
    push(a0 - 1 - g, a1 + 1 + g);
  } else if (mid) {
    push(a0, a1);
  } else if (g >= 0) {
    push(a0 - 1 - g, a0 - 1);
    push(a1 + 1, a1 + 1 + g);
  }
}

// Spans of u at which (moving + u) and `fixed` have intersecting interiors.
void appendOverlapSpan(const Rect& moving,
                       const Rect& fixed,
                       const Span& tmpl,
                       std::vector<Span>& out)
{
  if (!(moving.ylo < fixed.yhi && fixed.ylo < moving.yhi)) {
    return;
  }
  const int64_t lo = static_cast<int64_t>(fixed.xlo) - moving.xhi + 1;
  const int64_t hi = static_cast<int64_t>(fixed.xhi) - moving.xlo - 1;
  if (lo <= hi) {
    Span s = tmpl;
    s.lo = lo;
    s.hi = hi;
    out.push_back(s);
  }
}

// Furthest coordinate from c (in direction dir) up to which the answer
// equals `ref`; dir * kFar if it never changes. Per group with answer M at
// c, the answer changes where a span raising it begins (level > M, or any
// span if M = 0 in allowed/blocked mode) or where the union of the spans
// holding it at c (level M, or any span) ends; O(spans) plus sorting the
// few holding spans.
int64_t sweep(const std::vector<Span>& spans,
              const int64_t c,
              const int dir,
              const std::vector<RipupLevel>& ref,
              const bool levels)
{
  thread_local std::vector<int64_t> rise;
  thread_local std::vector<std::pair<int, std::pair<int64_t, int64_t>>> hold;
  const size_t groups = ref.size();
  rise.assign(groups, kFar);
  hold.clear();
  const int64_t cc = dir * c;
  for (const Span& s : spans) {
    // Mirror the leftward direction onto the rightward one.
    const int64_t lo = dir > 0 ? s.lo : -s.hi;
    const int64_t hi = dir > 0 ? s.hi : -s.lo;
    if (hi < cc) {
      continue;
    }
    const RipupLevel m = ref[s.group];
    if (lo > cc && (levels ? s.level > m : m == 0)) {
      rise[s.group] = std::min(rise[s.group], lo);
    }
    if (m > 0 && (!levels || s.level == m)) {
      hold.push_back({s.group, {lo, hi}});
    }
  }
  int64_t best = kFar;
  for (size_t g = 0; g < groups; g++) {
    best = std::min(best, rise[g]);
  }
  // End of the union of holding spans that contains c, per group (spans
  // chain rarely, so repeated scans beat sorting).
  thread_local std::vector<int64_t> end;
  end.assign(groups, cc - 1);
  for (bool grown = true; grown;) {
    grown = false;
    for (const auto& [g, span] : hold) {
      if (span.first <= end[g] + 1 && span.second > end[g]) {
        end[g] = span.second;
        grown = true;
      }
    }
  }
  for (size_t g = 0; g < groups; g++) {
    if (ref[g] > 0) {
      best = std::min(best, end[g] + 1);
    }
  }
  if (best >= kFar) {
    return dir * kFar;
  }
  return dir * (best - 1);
}

struct PieceExtensions
{
  int count = 0;
  Rect rect[2];
};

struct Scratch
{
  std::vector<int> layers;
  std::vector<Rect> base_windows;
  std::vector<GridWindow> windows;
  std::vector<std::vector<PieceExtensions>> piece_ext;
  std::vector<Span> xs;
  std::vector<Span> ys;
};

}  // namespace

void appendWireTypeModels(const WireType& type,
                          const int z,
                          const int num_wiring_layers,
                          const int base_group,
                          std::vector<LocatedModel>& models)
{
  auto add = [&](int layer, const WireModel& m, int group) {
    if (m.valid) {
      models.push_back({layer, m.rect, m.cls, 3, group});
    }
  };
  if (z < static_cast<int>(type.preferred.size())) {
    add(wiringGridLayer(z), type.preferred[z], base_group);
    add(wiringGridLayer(z), type.non_preferred[z], base_group + 1);
  }
  if (z + 1 < num_wiring_layers && z < static_cast<int>(type.via.size())
      && type.via[z].valid) {
    const ViaModel& v = type.via[z];
    add(wiringGridLayer(z), v.bottom, base_group + 2);
    add(viaGridLayer(z), v.cut, base_group + 2);
    add(wiringGridLayer(z + 1), v.top, base_group + 2);
    add(viaGridLayer(z + 1), v.projection, base_group + 2);
  }
  if (z > 0 && z - 1 < static_cast<int>(type.via.size())
      && type.via[z - 1].valid) {
    const ViaModel& v = type.via[z - 1];
    add(wiringGridLayer(z - 1), v.bottom, base_group + 3);
    add(viaGridLayer(z - 1), v.cut, base_group + 3);
    add(wiringGridLayer(z), v.top, base_group + 3);
    add(viaGridLayer(z), v.projection, base_group + 3);
  }
}

void DistanceRuleChecker::check(const Point l,
                                const std::span<const LocatedModel> models,
                                const CheckOptions& options,
                                CheckResult& result) const
{
  thread_local Scratch scratch;
  const ShapeClassTable& classes = grid_->classes();
  int groups = options.groups;
  for (size_t k = 0; k < models.size(); k++) {
    groups = std::max(
        groups,
        (models[k].group < 0 ? static_cast<int>(k) : models[k].group) + 1);
  }
  result.levels.assign(groups, 0);
  result.nets.clear();
  result.intervals = 0;
  result.x_lo_known_bound = result.x_hi_known_bound = false;
  result.y_lo_known_bound = result.y_hi_known_bound = false;
  // One shape-grid window per layer: all models bloated by the layer's
  // interaction distance.
  scratch.layers.clear();
  for (const LocatedModel& m : models) {
    if (std::ranges::find(scratch.layers, m.layer) == scratch.layers.end()) {
      scratch.layers.push_back(m.layer);
    }
  }
  const size_t num_layers = scratch.layers.size();
  if (scratch.windows.size() < num_layers) {
    scratch.windows.resize(num_layers);
    scratch.piece_ext.resize(num_layers);
  }
  scratch.base_windows.assign(num_layers, Rect{});
  for (size_t k = 0; k < num_layers; k++) {
    const int layer = scratch.layers[k];
    bool first = true;
    Rect hull;
    for (const LocatedModel& m : models) {
      if (m.layer == layer) {
        const Rect placed = m.rect.translated(l);
        hull = first ? placed : hull.hull(placed);
        first = false;
      }
    }
    const Coord reach = classes.reach(layer);
    scratch.base_windows[k] = hull.bloated(reach, reach);
    grid_->collect(
        layer,
        scratch.base_windows[k].bloated(options.margin_x, options.margin_y),
        scratch.windows[k]);
    result.intervals += scratch.windows[k].intervals;
  }
  // Locations whose windows lie in the known area (as shifts of l).
  int64_t kx_lo = -kFar;
  int64_t kx_hi = kFar;
  int64_t ky_lo = -kFar;
  int64_t ky_hi = kFar;
  for (size_t k = 0; k < num_layers; k++) {
    const Rect& w = scratch.base_windows[k];
    const GridWindow& gw = scratch.windows[k];
    if (gw.rows.empty()) {
      continue;
    }
    // Rows run along `along` (x for a horizontal layer).
    const Coord w_along_lo = gw.horizontal ? w.xlo : w.ylo;
    const Coord w_along_hi = gw.horizontal ? w.xhi : w.yhi;
    const Coord w_cross_lo = gw.horizontal ? w.ylo : w.xlo;
    const Coord w_cross_hi = gw.horizontal ? w.yhi : w.xhi;
    int64_t along_lo = -kFar;
    int64_t along_hi = kFar;
    for (const KnownRow& row : gw.rows) {
      if (row.cross_lo > w_cross_hi || row.cross_hi <= w_cross_lo) {
        continue;
      }
      along_lo
          = std::max<int64_t>(along_lo, int64_t{row.along_lo} - w_along_lo);
      along_hi
          = std::min<int64_t>(along_hi, int64_t{row.along_hi} - 1 - w_along_hi);
    }
    const int64_t cross_lo = int64_t{gw.rows.front().cross_lo} - w_cross_lo;
    const int64_t cross_hi = int64_t{gw.rows.back().cross_hi} - 1 - w_cross_hi;
    if (gw.horizontal) {
      kx_lo = std::max(kx_lo, along_lo);
      kx_hi = std::min(kx_hi, along_hi);
      ky_lo = std::max(ky_lo, cross_lo);
      ky_hi = std::min(ky_hi, cross_hi);
    } else {
      ky_lo = std::max(ky_lo, along_lo);
      ky_hi = std::min(ky_hi, along_hi);
      kx_lo = std::max(kx_lo, cross_lo);
      kx_hi = std::min(kx_hi, cross_hi);
    }
  }
  kx_lo = std::min<int64_t>(kx_lo, 0);
  kx_hi = std::max<int64_t>(kx_hi, 0);
  ky_lo = std::min<int64_t>(ky_lo, 0);
  ky_hi = std::max<int64_t>(ky_hi, 0);
  // Violation spans on the horizontal and the vertical line through l;
  // spans outside the known part of a line cannot change its interval.
  const bool want_x = !options.intervals || options.x_interval;
  const bool want_y = options.intervals && options.y_interval;
  const int64_t x_from = l.x + kx_lo - 1;
  const int64_t x_to = l.x + kx_hi + 1;
  const int64_t y_from = l.y + ky_lo - 1;
  const int64_t y_to = l.y + ky_hi + 1;
  auto keep
      = [](std::vector<Span>& spans, size_t from, int64_t lo, int64_t hi) {
          size_t out = from;
          for (size_t k = from; k < spans.size(); k++) {
            if (spans[k].hi >= lo && spans[k].lo <= hi) {
              spans[out++] = spans[k];
            }
          }
          spans.resize(out);
        };
  scratch.xs.clear();
  scratch.ys.clear();
  // Line-end extensions of the stored pieces, once per piece.
  for (size_t k = 0; k < num_layers; k++) {
    const GridWindow& w = scratch.windows[k];
    std::vector<PieceExtensions>& ext = scratch.piece_ext[k];
    ext.resize(w.pieces.size());
    for (size_t p = 0; p < w.pieces.size(); p++) {
      const GridPiece& piece = w.pieces[p];
      ext[p].count
          = piece.line_end_ends == 0
                ? 0
                : lineEndExtensions(
                      classes,
                      PlacedShape{piece.rect, piece.cls, piece.line_end_ends},
                      ext[p].rect);
    }
  }
  for (size_t mi = 0; mi < models.size(); mi++) {
    const LocatedModel& m = models[mi];
    const int group = m.group < 0 ? static_cast<int>(mi) : m.group;
    const size_t layer_index
        = std::ranges::find(scratch.layers, m.layer) - scratch.layers.begin();
    const GridWindow& window = scratch.windows[layer_index];
    const std::vector<PieceExtensions>& piece_ext
        = scratch.piece_ext[layer_index];
    const PlacedShape placed{m.rect.translated(l), m.cls, m.line_end_ends};
    const Rect mx{m.rect.xlo, placed.rect.ylo, m.rect.xhi, placed.rect.yhi};
    const Rect my
        = Rect{placed.rect.xlo, m.rect.ylo, placed.rect.xhi, m.rect.yhi}
              .transposed();
    Rect model_ext[2];
    const int model_exts = lineEndExtensions(classes, placed, model_ext);
    // Requirements between the model's class and the pieces' classes.
    thread_local std::vector<std::pair<uint64_t, Requirement>> reqs;
    thread_local uint64_t req_stamp = 0;
    req_stamp++;
    if (reqs.size() < classes.size()) {
      reqs.resize(classes.size() + 64, {0, Requirement{}});
    }
    for (size_t pi = 0; pi < window.pieces.size(); pi++) {
      const GridPiece& p = window.pieces[pi];
      if (p.owner == options.ignored_owner && p.owner != kNoNet
          && isOwned(p.type)) {
        continue;
      }
      auto& cached = reqs[p.cls];
      if (cached.first != req_stamp) {
        cached = {req_stamp, classes.requirement(m.cls, p.cls)};
      }
      const Requirement& req = cached.second;
      if (!req.overlap && !req.line_end && req.facing_x == 0
          && req.facing_y == 0 && req.diagonal == 0) {
        continue;
      }
      const Span tmpl{0,
                      0,
                      group,
                      grid_->level(p.type, p.owner),
                      p.owner,
                      isRemovable(p.type)};
      const Rect pt = p.rect.transposed();
      const size_t x0 = scratch.xs.size();
      const size_t y0 = scratch.ys.size();
      if (want_x) {
        appendSpacingSpans(mx,
                           p.rect,
                           req.facing_x,
                           req.facing_y,
                           req.diagonal,
                           req.overlap,
                           tmpl,
                           scratch.xs);
      }
      if (want_y) {
        appendSpacingSpans(my,
                           pt,
                           req.facing_y,
                           req.facing_x,
                           req.diagonal,
                           req.overlap,
                           tmpl,
                           scratch.ys);
      }
      if (req.line_end) {
        // Line-end extensions of the model (Sec. 3.1, Fig. 2).
        for (int e = 0; e < model_exts; e++) {
          const Rect& r = model_ext[e];
          if (want_x) {
            appendOverlapSpan({r.xlo - l.x, r.ylo, r.xhi - l.x, r.yhi},
                              p.rect,
                              tmpl,
                              scratch.xs);
          }
          if (want_y) {
            appendOverlapSpan(
                Rect{r.xlo, r.ylo - l.y, r.xhi, r.yhi - l.y}.transposed(),
                pt,
                tmpl,
                scratch.ys);
          }
        }
        // Line-end extensions of the stored shape.
        for (int e = 0; e < piece_ext[pi].count; e++) {
          const Rect& r = piece_ext[pi].rect[e];
          if (want_x) {
            appendOverlapSpan(mx, r, tmpl, scratch.xs);
          }
          if (want_y) {
            appendOverlapSpan(my, r.transposed(), tmpl, scratch.ys);
          }
        }
      }
      keep(scratch.xs, x0, x_from, x_to);
      keep(scratch.ys, y0, y_from, y_to);
    }
  }
  // Answer at l.
  const std::vector<Span>& at_l = want_x ? scratch.xs : scratch.ys;
  const int64_t c = want_x ? l.x : l.y;
  for (const Span& s : at_l) {
    if (s.lo <= c && c <= s.hi) {
      result.levels[s.group] = std::max(result.levels[s.group], s.level);
    }
  }
  if (options.report_nets) {
    for (const Span& s : at_l) {
      if (s.lo > c || c > s.hi || !s.removable
          || result.levels[s.group] >= kNotRemovable) {
        continue;
      }
      auto at = std::ranges::find(result.nets, s.owner, &BlockingNet::net);
      if (at == result.nets.end()) {
        result.nets.push_back({s.owner, 0, 0});
        at = result.nets.end() - 1;
      }
      at->level = std::max(at->level, s.level);
      if (s.group < 64) {
        at->groups |= uint64_t{1} << s.group;
      }
    }
    std::ranges::sort(result.nets, {}, &BlockingNet::net);
  }
  result.x_lo = result.x_hi = l.x;
  result.y_lo = result.y_hi = l.y;
  if (!options.intervals) {
    return;
  }
  auto clampCoord = [](int64_t v) {
    return static_cast<Coord>(
        std::clamp<int64_t>(v, -(int64_t{1} << 30), int64_t{1} << 30));
  };
  auto bound = [&](const std::vector<Span>& spans,
                   int64_t c,
                   int64_t known_lo,
                   int64_t known_hi,
                   Coord& lo,
                   Coord& hi,
                   bool& lo_known,
                   bool& hi_known) {
    const int64_t a_lo = sweep(spans, c, -1, result.levels, options.levels);
    const int64_t a_hi = sweep(spans, c, 1, result.levels, options.levels);
    const int64_t k_lo = c + known_lo;
    const int64_t k_hi = c + known_hi;
    lo_known = k_lo >= a_lo;
    hi_known = k_hi <= a_hi;
    lo = clampCoord(std::max(a_lo, k_lo));
    hi = clampCoord(std::min(a_hi, k_hi));
  };
  if (want_x) {
    bound(scratch.xs,
          l.x,
          kx_lo,
          kx_hi,
          result.x_lo,
          result.x_hi,
          result.x_lo_known_bound,
          result.x_hi_known_bound);
  }
  if (want_y) {
    bound(scratch.ys,
          l.y,
          ky_lo,
          ky_hi,
          result.y_lo,
          result.y_hi,
          result.y_lo_known_bound,
          result.y_hi_known_bound);
  }
}

void DistanceRuleChecker::query(
    const int z,
    const Point location,
    const std::span<const WireType* const> wire_types,
    const CheckOptions& options,
    std::vector<WireTypeAnswer>& answers,
    CheckResult& result) const
{
  thread_local std::vector<LocatedModel> models;
  models.clear();
  const int num_z = (grid_->numLayers() + 1) / 2;
  for (size_t t = 0; t < wire_types.size(); t++) {
    appendWireTypeModels(
        *wire_types[t], z, num_z, 4 * static_cast<int>(t), models);
  }
  // Groups without a model (no via there) stay legal.
  CheckOptions opts = options;
  opts.groups = std::max(opts.groups, 4 * static_cast<int>(wire_types.size()));
  check(location, models, opts, result);
  answers.assign(wire_types.size(), {});
  for (size_t t = 0; t < wire_types.size(); t++) {
    WireTypeAnswer& a = answers[t];
    a.preferred = result.levels[4 * t];
    a.non_preferred = result.levels[4 * t + 1];
    a.via_up = result.levels[4 * t + 2];
    a.via_down = result.levels[4 * t + 3];
  }
}

RipupLevel DistanceRuleChecker::shapeLevel(const int layer,
                                           const Rect& rect,
                                           const ClassId cls,
                                           const uint8_t line_end_ends,
                                           const int32_t ignored_owner,
                                           std::vector<BlockingNet>* nets) const
{
  thread_local GridWindow window;
  const ShapeClassTable& classes = grid_->classes();
  const Coord reach = classes.reach(layer);
  grid_->collect(layer, rect.bloated(reach, reach), window);
  const PlacedShape shape{rect, cls, line_end_ends};
  RipupLevel level = 0;
  if (nets != nullptr) {
    nets->clear();
  }
  for (const GridPiece& p : window.pieces) {
    if (p.owner == ignored_owner && p.owner != kNoNet && isOwned(p.type)) {
      continue;
    }
    if (!shapesViolate(
            classes, shape, PlacedShape{p.rect, p.cls, p.line_end_ends})) {
      continue;
    }
    const RipupLevel l = grid_->level(p.type, p.owner);
    level = std::max(level, l);
    if (nets != nullptr && isRemovable(p.type)) {
      auto at = std::ranges::find(*nets, p.owner, &BlockingNet::net);
      if (at == nets->end()) {
        nets->push_back({p.owner, 0, 1});
        at = nets->end() - 1;
      }
      at->level = std::max(at->level, l);
    }
  }
  if (nets != nullptr) {
    if (level >= kNotRemovable) {
      nets->clear();
    }
    std::ranges::sort(*nets, {}, &BlockingNet::net);
  }
  return level;
}

}  // namespace drt::br
