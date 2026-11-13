// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrOffTrackSearch.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

namespace drt::br {

namespace {

constexpr int kFresh = 4;
// States per grid vertex: four incoming directions and "no incoming
// direction"; one where the incoming direction restricts no arc.
constexpr uint32_t kDirectionStates = 5;
constexpr uint32_t kNoParent = std::numeric_limits<uint32_t>::max();
constexpr int64_t kInf = std::numeric_limits<int64_t>::max() / 4;

int exactIndex(const std::vector<GeoCoord>& v, const GeoCoord c)
{
  auto it = std::lower_bound(v.begin(), v.end(), c);
  return (it != v.end() && *it == c) ? static_cast<int>(it - v.begin()) : -1;
}

// Open intervals of the obstacles whose open interior meets the line
// `along` = c (horizontal line y = c if `rows`), sorted and merged where
// they overlap (touching open intervals stay apart: the common end point is
// free), appended to `out`. `sorted` holds the obstacles in increasing ylo
// (rows) or xlo (columns) order; only those below c can meet the line.
void lineIntervals(const std::vector<Box>& sorted,
                   const GeoCoord c,
                   const bool rows,
                   std::vector<std::pair<GeoCoord, GeoCoord>>& scratch,
                   std::vector<std::pair<GeoCoord, GeoCoord>>& out)
{
  scratch.clear();
  for (const Box& o : sorted) {
    if (rows ? o.ylo >= c : o.xlo >= c) {
      break;
    }
    if (rows ? c < o.yhi : c < o.xhi) {
      scratch.emplace_back(rows ? o.xlo : o.ylo, rows ? o.xhi : o.yhi);
    }
  }
  std::sort(scratch.begin(), scratch.end());
  const size_t first = out.size();
  for (const auto& [lo, hi] : scratch) {
    if (lo >= hi) {
      continue;
    }
    if (out.size() > first && lo < out.back().second) {
      out.back().second = std::max(out.back().second, hi);
    } else {
      out.emplace_back(lo, hi);
    }
  }
}

// True if [a, b] meets one of the open intervals [begin, end).
bool meets(const std::pair<GeoCoord, GeoCoord>* begin,
           const std::pair<GeoCoord, GeoCoord>* end,
           const GeoCoord a,
           const GeoCoord b)
{
  // First interval with hi > a.
  auto it = std::upper_bound(
      begin, end, a, [](GeoCoord v, const auto& p) { return v < p.second; });
  return it != end && it->first < b;
}

// Distance of c to the nearest value of the sorted v (v not empty).
int64_t nearest(const std::vector<GeoCoord>& v, const GeoCoord c)
{
  auto it = std::lower_bound(v.begin(), v.end(), c);
  int64_t d = std::numeric_limits<int64_t>::max() / 8;
  if (it != v.end()) {
    d = int64_t{*it} - c;
  }
  if (it != v.begin()) {
    d = std::min(d, int64_t{c} - *(it - 1));
  }
  return d;
}

// True if consecutive values of the sorted v are at least tau apart.
bool gapsAtLeast(const std::vector<GeoCoord>& v, const GeoCoord tau)
{
  for (size_t i = 1; i < v.size(); i++) {
    if (int64_t{v[i]} - v[i - 1] < tau) {
      return false;
    }
  }
  return true;
}

// Index of each line of `from` in `to` (both sorted), or -1.
void mapLines(const std::vector<GeoCoord>& from,
              const std::vector<GeoCoord>& to,
              std::vector<int>& out)
{
  out.assign(from.size(), -1);
  size_t k = 0;
  for (size_t i = 0; i < from.size(); i++) {
    while (k < to.size() && to[k] < from[i]) {
      k++;
    }
    if (k < to.size() && to[k] == from[i]) {
      out[i] = static_cast<int>(k);
    }
  }
}

// plus[i] = first index with v >= v[i] + tau (n if none); minus[i] = last
// index with v <= v[i] - tau (-1 if none).
void turnLines(const std::vector<GeoCoord>& v,
               const GeoCoord tau,
               std::vector<int>& plus,
               std::vector<int>& minus)
{
  const int n = static_cast<int>(v.size());
  plus.resize(n);
  minus.resize(n);
  for (int i = 0; i < n; i++) {
    plus[i] = static_cast<int>(std::lower_bound(v.begin(), v.end(), v[i] + tau)
                               - v.begin());
    minus[i] = static_cast<int>(std::upper_bound(v.begin(), v.end(), v[i] - tau)
                                - v.begin())
               - 1;
  }
}

}  // namespace

namespace {

// Prefix counts of the unit segments [lines[k], lines[k + 1]] that meet one
// of the sorted, disjoint open intervals [begin, end): out[k] = number of
// blocked segments among the first k.
void blockedPrefix(const std::vector<GeoCoord>& lines,
                   const std::pair<GeoCoord, GeoCoord>* begin,
                   const std::pair<GeoCoord, GeoCoord>* end,
                   int32_t* out)
{
  const int n = static_cast<int>(lines.size());
  int32_t count = 0;
  const std::pair<GeoCoord, GeoCoord>* it = begin;
  for (int k = 0; k + 1 < n; k++) {
    out[k] = count;
    // First interval with hi > lines[k] (the intervals are disjoint and
    // sorted, so their upper ends increase).
    while (it != end && it->second <= lines[k]) {
      ++it;
    }
    count += it != end && it->first < lines[k + 1];
  }
  if (n > 0) {
    out[n - 1] = count;
  }
}

}  // namespace

bool OffTrackSearch::freeH(LayerIndex& l,
                           const int iy,
                           const int ix0,
                           const int ix1)
{
  std::pair<int, int>& span = l.row_span[iy];
  int32_t* blocked = l.row_blocked.data() + static_cast<size_t>(iy) * l.nx;
  if (span.first < 0) {
    thread_local std::vector<std::pair<GeoCoord, GeoCoord>> scratch;
    span.first = static_cast<int>(l.intervals.size());
    lineIntervals(l.by_ylo, l.ys[iy], true, scratch, l.intervals);
    span.second = static_cast<int>(l.intervals.size());
    const auto* data = l.intervals.data();
    blockedPrefix(l.xs, data + span.first, data + span.second, blocked);
  }
  if (ix0 == ix1) {
    const auto* data = l.intervals.data();
    return !meets(data + span.first, data + span.second, l.xs[ix0], l.xs[ix1]);
  }
  return blocked[ix1] == blocked[ix0];
}

bool OffTrackSearch::freeV(LayerIndex& l,
                           const int ix,
                           const int iy0,
                           const int iy1)
{
  std::pair<int, int>& span = l.col_span[ix];
  int32_t* blocked = l.col_blocked.data() + static_cast<size_t>(ix) * l.ny;
  if (span.first < 0) {
    thread_local std::vector<std::pair<GeoCoord, GeoCoord>> scratch;
    span.first = static_cast<int>(l.intervals.size());
    lineIntervals(l.by_xlo, l.xs[ix], false, scratch, l.intervals);
    span.second = static_cast<int>(l.intervals.size());
    const auto* data = l.intervals.data();
    blockedPrefix(l.ys, data + span.first, data + span.second, blocked);
  }
  if (iy0 == iy1) {
    const auto* data = l.intervals.data();
    return !meets(data + span.first, data + span.second, l.ys[iy0], l.ys[iy1]);
  }
  return blocked[iy1] == blocked[iy0];
}

OffTrackResult OffTrackSearch::run(const OffTrackRequest& req)
{
  OffTrackResult res;
  const int num_layers = static_cast<int>(req.layers.size());
  if (static_cast<int>(index_.size()) < num_layers) {
    index_.resize(num_layers);
  }
  uint64_t total = 0;
  const Box& win = req.window;
  for (int z = 0; z < num_layers; z++) {
    const OffTrackLayer& in = req.layers[z];
    LayerIndex& l = index_[z];
    l.xs.clear();
    l.ys.clear();
    for (const GeoCoord x : in.grid.xs) {
      if (x >= win.xlo && x <= win.xhi) {
        l.xs.push_back(x);
      }
    }
    for (const GeoCoord y : in.grid.ys) {
      if (y >= win.ylo && y <= win.yhi) {
        l.ys.push_back(y);
      }
    }
    l.nx = static_cast<int>(l.xs.size());
    l.ny = static_cast<int>(l.ys.size());
    l.base = static_cast<uint32_t>(std::min<uint64_t>(total, kNoParent));
    total += static_cast<uint64_t>(l.nx) * l.ny;
  }
  // One state per grid vertex when every planar arc between neighbouring
  // grid lines is a tau-feasible segment and no arc filter reads the
  // incoming direction: a path of that graph is a tau-feasible path of the
  // grid and a path of G is one of that graph, so both have the same
  // distances.
  bool single_state = req.merge_direction_states && !req.arc_filter;
  for (int z = 0; z < num_layers && single_state; z++) {
    const LayerIndex& l = index_[z];
    single_state = !req.layers[z].wires
                   || (gapsAtLeast(l.xs, req.layers[z].grid.tau)
                       && gapsAtLeast(l.ys, req.layers[z].grid.tau));
  }
  const uint32_t stride = single_state ? 1 : kDirectionStates;
  // The "no incoming direction" state of a grid vertex.
  const uint32_t fresh = single_state ? 0 : kFresh;
  // Tie-break scale: every simple path's off-track preferred length is
  // below it (per layer at most the window length per line), so
  // cost * scale + off-track length orders lexicographically. Off (scale
  // 1) if that would leave less than 2^38 for the cost.
  int64_t scale = 1;
  if (req.on_track_ties) {
    int64_t bound = 1;
    for (int z = 0; z < num_layers; z++) {
      const LayerIndex& l = index_[z];
      if (req.layers[z].tie_tracks == nullptr || l.nx == 0 || l.ny == 0) {
        continue;
      }
      bound += req.layers[z].horizontal
                   ? int64_t{l.ny} * (int64_t{l.xs.back()} - l.xs.front())
                   : int64_t{l.nx} * (int64_t{l.ys.back()} - l.ys.front());
    }
    scale = bound < (int64_t{1} << 24) ? bound : 1;
  }
  res.grid_vertices = static_cast<int64_t>(total);
  res.states = static_cast<int64_t>(total) * stride;
  if (total * stride >= kNoParent) {
    return res;
  }
  for (int z = 0; z < num_layers; z++) {
    const OffTrackLayer& in = req.layers[z];
    LayerIndex& l = index_[z];
    l.obstacles = &in.obstacles;
    l.by_ylo = in.obstacles;
    std::stable_sort(l.by_ylo.begin(),
                     l.by_ylo.end(),
                     [](const Box& a, const Box& b) { return a.ylo < b.ylo; });
    l.by_xlo = in.obstacles;
    std::stable_sort(l.by_xlo.begin(),
                     l.by_xlo.end(),
                     [](const Box& a, const Box& b) { return a.xlo < b.xlo; });
    l.row_span.assign(l.ny, {-1, -1});
    l.col_span.assign(l.nx, {-1, -1});
    l.intervals.clear();
    const size_t cells = static_cast<size_t>(l.nx) * l.ny;
    if (l.row_blocked.size() < cells) {
      l.row_blocked.resize(cells);
      l.col_blocked.resize(cells);
    }
    turnLines(l.xs, in.grid.tau, l.turn_x_plus, l.turn_x_minus);
    turnLines(l.ys, in.grid.tau, l.turn_y_plus, l.turn_y_minus);
    l.row_on_track.clear();
    l.col_on_track.clear();
    l.line_off_track.clear();
    if (req.on_track_ties && in.tie_tracks != nullptr) {
      const std::vector<GeoCoord>& lines = in.horizontal ? l.ys : l.xs;
      l.line_off_track.resize(lines.size());
      for (size_t k = 0; k < lines.size(); k++) {
        l.line_off_track[k] = !std::binary_search(
            in.tie_tracks->begin(), in.tie_tracks->end(), lines[k]);
      }
    }
    if (in.tracks != nullptr) {
      // Planar arcs start on a track line in the preferred direction.
      if (in.horizontal) {
        l.row_on_track.resize(l.ny);
        for (int iy = 0; iy < l.ny; iy++) {
          l.row_on_track[iy] = std::binary_search(
              in.tracks->begin(), in.tracks->end(), l.ys[iy]);
        }
      } else {
        l.col_on_track.resize(l.nx);
        for (int ix = 0; ix < l.nx; ix++) {
          l.col_on_track[ix] = std::binary_search(
              in.tracks->begin(), in.tracks->end(), l.xs[ix]);
        }
      }
    }
  }
  const size_t states = static_cast<size_t>(total) * stride;
  if (labels_.size() < states) {
    labels_.resize(states, Label{0, kNoParent, 0});
  }
  if (via_memo_.size() < total) {
    via_memo_.resize(total);
    via_stamp_.resize(total, 0);
    target_first_.resize(total);
    target_stamp_.resize(total, 0);
  }
  if (++epoch_ == 0) {
    for (Label& label : labels_) {
      label.stamp = 0;
    }
    std::fill(via_stamp_.begin(), via_stamp_.end(), 0);
    std::fill(target_stamp_.begin(), target_stamp_.end(), 0);
    epoch_ = 1;
  }
  for (int z = 0; z < num_layers; z++) {
    LayerIndex& l = index_[z];
    if (z + 1 < num_layers) {
      mapLines(l.xs, index_[z + 1].xs, l.up_x);
      mapLines(l.ys, index_[z + 1].ys, l.up_y);
    }
    if (z > 0) {
      mapLines(l.xs, index_[z - 1].xs, l.down_x);
      mapLines(l.ys, index_[z - 1].ys, l.down_y);
    }
  }
  Label* const labels = labels_.data();
  const uint32_t epoch = epoch_;
  auto dist = [labels, epoch](const uint32_t st) {
    return labels[st].stamp == epoch ? labels[st].dist : kInf;
  };

  auto layer_of = [&](const uint32_t v) {
    int z = num_layers - 1;
    while (z > 0 && index_[z].base > v) {
      z--;
    }
    return z;
  };
  // Future cost (A*, as the on-track search of Sec. 4.1): l1 distance to
  // the lattice spanned by the target coordinates of a layer plus the via
  // costs to reach it. The lattice contains every target, so the bound is
  // admissible; it changes by at most an arc's cost along any arc, so it is
  // consistent, every state is settled with its exact distance, and targets
  // (future cost 0) are settled in increasing distance.
  std::vector<std::vector<GeoCoord>> target_xs(num_layers);
  std::vector<std::vector<GeoCoord>> target_ys(num_layers);
  int64_t unit_cost = std::numeric_limits<int64_t>::max();
  for (const OffTrackLayer& l : req.layers) {
    unit_cost = std::min({unit_cost, l.cost_x, l.cost_y});
  }
  unit_cost = std::max<int64_t>(0, unit_cost);
  for (const OffTrackTarget& t : req.targets) {
    if (t.layer >= 0 && t.layer < num_layers) {
      target_xs[t.layer].push_back(t.p.x);
      target_ys[t.layer].push_back(t.p.y);
    }
  }
  std::vector<int> target_layers;
  for (int z = 0; z < num_layers; z++) {
    std::sort(target_xs[z].begin(), target_xs[z].end());
    target_xs[z].erase(std::unique(target_xs[z].begin(), target_xs[z].end()),
                       target_xs[z].end());
    std::sort(target_ys[z].begin(), target_ys[z].end());
    target_ys[z].erase(std::unique(target_ys[z].begin(), target_ys[z].end()),
                       target_ys[z].end());
    if (!target_xs[z].empty()) {
      target_layers.push_back(z);
    }
  }
  const int num_target_layers = static_cast<int>(target_layers.size());
  if (req.future_cost) {
    for (int z = 0; z < num_layers; z++) {
      LayerIndex& l = index_[z];
      l.near_x.resize(static_cast<size_t>(num_target_layers) * l.nx);
      l.near_y.resize(static_cast<size_t>(num_target_layers) * l.ny);
      for (int k = 0; k < num_target_layers; k++) {
        const int t = target_layers[k];
        for (int ix = 0; ix < l.nx; ix++) {
          l.near_x[static_cast<size_t>(k) * l.nx + ix]
              = nearest(target_xs[t], l.xs[ix]);
        }
        for (int iy = 0; iy < l.ny; iy++) {
          l.near_y[static_cast<size_t>(k) * l.ny + iy]
              = nearest(target_ys[t], l.ys[iy]);
        }
      }
    }
  }
  // Future cost of grid vertex (ix, iy) of layer z.
  auto futureCost = [&](const int z, const int ix, const int iy) -> int64_t {
    const LayerIndex& l = index_[z];
    int64_t h = kInf;
    for (int k = 0; k < num_target_layers; k++) {
      const int t = target_layers[k];
      h = std::min(h,
                   (l.near_x[static_cast<size_t>(k) * l.nx + ix]
                    + l.near_y[static_cast<size_t>(k) * l.ny + iy])
                           * unit_cost
                       + req.via_cost * std::abs(t - z));
    }
    return h == kInf ? 0 : h;
  };
  auto future = [&](const int z, const int ix, const int iy) -> int64_t {
    if (!req.future_cost) {
      return 0;
    }
    return scale * futureCost(z, ix, iy);
  };
  // A binary heap ordered by f only (std::priority_queue's push / pop).
  heap_.clear();
  auto heap_less
      = [](const HeapItem& a, const HeapItem& b) { return a.f > b.f; };
  auto push = [&](const int64_t f, const uint32_t state) {
    heap_.push_back({f, state});
    std::push_heap(heap_.begin(), heap_.end(), heap_less);
  };
  for (const OffTrackSource& s : req.sources) {
    if (s.layer < 0 || s.layer >= num_layers) {
      continue;
    }
    const LayerIndex& l = index_[s.layer];
    const int x0 = static_cast<int>(
        std::lower_bound(l.xs.begin(), l.xs.end(), s.region.xlo)
        - l.xs.begin());
    const int y0 = static_cast<int>(
        std::lower_bound(l.ys.begin(), l.ys.end(), s.region.ylo)
        - l.ys.begin());
    for (int ix = x0; ix < l.nx && l.xs[ix] <= s.region.xhi; ix++) {
      for (int iy = y0; iy < l.ny && l.ys[iy] <= s.region.yhi; iy++) {
        const uint32_t v = l.base + static_cast<uint32_t>(ix) * l.ny + iy;
        const uint32_t st = v * stride + fresh;
        if (dist(st) > s.cost * scale) {
          labels[st] = {s.cost * scale, kNoParent, epoch};
          push(s.cost * scale + future(s.layer, ix, iy), st);
        }
      }
    }
  }
  // Target vertices -> target ids.
  std::vector<std::pair<uint32_t, int>> targets;
  for (const OffTrackTarget& t : req.targets) {
    if (t.layer < 0 || t.layer >= num_layers) {
      continue;
    }
    const LayerIndex& l = index_[t.layer];
    const int ix = exactIndex(l.xs, t.p.x);
    const int iy = exactIndex(l.ys, t.p.y);
    if (ix < 0 || iy < 0) {
      continue;
    }
    targets.emplace_back(l.base + static_cast<uint32_t>(ix) * l.ny + iy, t.id);
  }
  std::sort(targets.begin(), targets.end());
  std::vector<int64_t> reached(targets.size(), -1);
  size_t distinct = 0;
  for (size_t k = 0; k < targets.size(); k++) {
    if (k == 0 || targets[k].first != targets[k - 1].first) {
      distinct++;
      // The entry a settled vertex checks: the first one at or after
      // (vertex, -1).
      const uint32_t v = targets[k].first;
      const auto it = std::lower_bound(
          targets.begin(), targets.end(), std::make_pair(v, -1));
      if (it != targets.end() && it->first == v) {
        target_stamp_[v] = epoch;
        target_first_[v] = static_cast<int32_t>(it - targets.begin());
      }
    }
  }
  size_t remaining = distinct;

  auto via_ok
      = [&](const int lower_z, const uint32_t lower_v, const GridPoint p) {
          if (via_stamp_[lower_v] != epoch) {
            via_stamp_[lower_v] = epoch;
            via_memo_[lower_v]
                = (req.via_allowed && req.via_allowed(lower_z, p)) ? 1 : 0;
          }
          return via_memo_[lower_v] == 1;
        };
  const bool has_filter = static_cast<bool>(req.arc_filter);

  while (!heap_.empty() && remaining > 0) {
    const HeapItem top = heap_.front();
    std::pop_heap(heap_.begin(), heap_.end(), heap_less);
    heap_.pop_back();
    const uint32_t st = top.state;
    const uint32_t v = st / stride;
    const int d = single_state ? kFresh : static_cast<int>(st % stride);
    const int z = layer_of(v);
    LayerIndex& l = index_[z];
    const uint32_t local = v - l.base;
    const int ix = static_cast<int>(local / static_cast<uint32_t>(l.ny));
    const int iy = static_cast<int>(local % static_cast<uint32_t>(l.ny));
    // g = f - future: the future cost of a vertex is a fixed value.
    const int64_t cost = top.f - future(z, ix, iy);
    if (cost > dist(st)) {
      continue;
    }
    if (cost / scale > req.max_cost) {
      break;
    }
    res.settled++;
    const OffTrackLayer& in = req.layers[z];
    const GridPoint p{l.xs[ix], l.ys[iy]};
    if (target_stamp_[v] == epoch) {
      const size_t k = target_first_[v];
      if (reached[k] < 0) {
        for (size_t q = k; q < targets.size() && targets[q].first == v; q++) {
          reached[q] = st;
        }
        remaining--;
        if (req.max_targets > 0
            && distinct - remaining >= static_cast<size_t>(req.max_targets)) {
          remaining = 0;
        }
      }
    }
    auto relax = [&](const uint32_t to,
                     const int64_t to_cost,
                     const int to_layer,
                     const int jx,
                     const int jy) {
      res.arcs++;
      if (to_cost < dist(to)) {
        labels[to] = {to_cost, st, epoch};
        push(to_cost + future(to_layer, jx, jy), to);
      }
    };
    const uint32_t row0 = l.base;
    auto state = [&](const int jx, const int jy, const int dir) {
      return (row0 + static_cast<uint32_t>(jx) * l.ny + jy) * stride
             + (single_state ? 0 : dir);
    };
    auto planar = [&](const int jx, const int jy, const int dir) {
      const GridPoint q{l.xs[jx], l.ys[jy]};
      if (has_filter && !req.arc_filter(z, p, q, d, dir)) {
        return;
      }
      int64_t c = (std::abs(int64_t{q.x} - p.x) * in.cost_x
                   + std::abs(int64_t{q.y} - p.y) * in.cost_y)
                  * scale;
      if (!l.line_off_track.empty() && (dir < 2) == in.horizontal
          && l.line_off_track[in.horizontal ? iy : ix]) {
        // Tie-break: preferred-direction wire off the tracks.
        c += std::abs(int64_t{q.x} - p.x) + std::abs(int64_t{q.y} - p.y);
      }
      relax(state(jx, jy, dir), cost + c, z, jx, jy);
    };
    if (in.wires) {
      // On-track-only layers: arcs along the preferred direction, from a
      // track line.
      bool arcs_x = true;
      bool arcs_y = true;
      if (in.tracks != nullptr) {
        arcs_x = in.horizontal && l.row_on_track[iy];
        arcs_y = !in.horizontal && l.col_on_track[ix];
      }
      // Straight arc to the neighbour in the incoming direction.
      if (d == 0) {
        if (ix + 1 < l.nx && freeH(l, iy, ix, ix + 1) && arcs_x) {
          planar(ix + 1, iy, 0);
        }
      } else if (d == 1) {
        if (ix > 0 && freeH(l, iy, ix - 1, ix) && arcs_x) {
          planar(ix - 1, iy, 1);
        }
      } else if (d == 2) {
        if (iy + 1 < l.ny && freeV(l, ix, iy, iy + 1) && arcs_y) {
          planar(ix, iy + 1, 2);
        }
      } else if (d == 3) {
        if (iy > 0 && freeV(l, ix, iy - 1, iy) && arcs_y) {
          planar(ix, iy - 1, 3);
        }
      }
      // Turn arcs to the nearest vertices at distance >= tau.
      const bool turn_x = !single_state && (d == kFresh || d == 2 || d == 3);
      const bool turn_y = !single_state && (d == kFresh || d == 0 || d == 1);
      if (single_state) {
        // Arcs to the four neighbours.
        if (ix + 1 < l.nx && freeH(l, iy, ix, ix + 1) && arcs_x) {
          planar(ix + 1, iy, 0);
        }
        if (ix > 0 && freeH(l, iy, ix - 1, ix) && arcs_x) {
          planar(ix - 1, iy, 1);
        }
        if (iy + 1 < l.ny && freeV(l, ix, iy, iy + 1) && arcs_y) {
          planar(ix, iy + 1, 2);
        }
        if (iy > 0 && freeV(l, ix, iy - 1, iy) && arcs_y) {
          planar(ix, iy - 1, 3);
        }
      }
      if (turn_x) {
        const int kp = l.turn_x_plus[ix];
        if (kp < l.nx && freeH(l, iy, ix, kp) && arcs_x) {
          planar(kp, iy, 0);
        }
        const int km = l.turn_x_minus[ix];
        if (km >= 0 && freeH(l, iy, km, ix) && arcs_x) {
          planar(km, iy, 1);
        }
      }
      if (turn_y) {
        const int kp = l.turn_y_plus[iy];
        if (kp < l.ny && freeV(l, ix, iy, kp) && arcs_y) {
          planar(ix, kp, 2);
        }
        const int km = l.turn_y_minus[iy];
        if (km >= 0 && freeV(l, ix, km, iy) && arcs_y) {
          planar(ix, km, 3);
        }
      }
    }
    // Via arcs to the "no incoming direction" state of the adjacent layer.
    for (const int nz : {z + 1, z - 1}) {
      if (nz < 0 || nz >= num_layers) {
        continue;
      }
      const LayerIndex& n = index_[nz];
      const int jx = nz > z ? l.up_x[ix] : l.down_x[ix];
      const int jy = nz > z ? l.up_y[iy] : l.down_y[iy];
      if (jx < 0 || jy < 0) {
        continue;
      }
      const uint32_t w = n.base + static_cast<uint32_t>(jx) * n.ny + jy;
      const bool ok = nz > z ? via_ok(z, v, p) : via_ok(nz, w, p);
      if (ok) {
        relax(w * stride + fresh, cost + req.via_cost * scale, nz, jx, jy);
      }
    }
  }

  // Paths to the reached targets.
  for (size_t k = 0; k < targets.size(); k++) {
    if (reached[k] < 0) {
      continue;
    }
    OffTrackPath path;
    path.target_id = targets[k].second;
    path.cost = labels[reached[k]].dist / scale;
    std::vector<RouteNode> full;
    for (uint32_t s = static_cast<uint32_t>(reached[k]); s != kNoParent;
         s = labels[s].parent) {
      const uint32_t v = s / stride;
      const int z = layer_of(v);
      const LayerIndex& l = index_[z];
      const uint32_t local = v - l.base;
      const int ix = static_cast<int>(local / static_cast<uint32_t>(l.ny));
      const int iy = static_cast<int>(local % static_cast<uint32_t>(l.ny));
      full.push_back({z, {l.xs[ix], l.ys[iy]}});
    }
    std::reverse(full.begin(), full.end());
    const size_t n = full.size();
    for (size_t q = 0; q < n; q++) {
      bool keep = q == 0 || q + 1 == n;
      if (!keep) {
        const RouteNode& a = full[q - 1];
        const RouteNode& b = full[q];
        const RouteNode& c = full[q + 1];
        keep = a.layer != b.layer || c.layer != b.layer
               || !((a.p.x == b.p.x && b.p.x == c.p.x)
                    || (a.p.y == b.p.y && b.p.y == c.p.y));
      }
      if (keep) {
        path.nodes.push_back(full[q]);
      }
    }
    res.paths.push_back(std::move(path));
  }
  std::stable_sort(res.paths.begin(),
                   res.paths.end(),
                   [](const OffTrackPath& a, const OffTrackPath& b) {
                     return a.cost < b.cost;
                   });
  return res;
}

namespace {

// Pieces per (layer, horizontal, line) as sorted disjoint intervals.
using WireKey = std::tuple<int, bool, GeoCoord>;
using Intervals = std::vector<std::pair<GeoCoord, GeoCoord>>;
using WireSet = std::map<WireKey, Intervals>;

Intervals mergedIntervals(Intervals in)
{
  std::sort(in.begin(), in.end());
  Intervals out;
  for (const auto& [lo, hi] : in) {
    if (hi <= lo) {
      continue;
    }
    if (!out.empty() && lo <= out.back().second) {
      out.back().second = std::max(out.back().second, hi);
    } else {
      out.emplace_back(lo, hi);
    }
  }
  return out;
}

WireSet wireSet(const std::vector<CatalogueWire>& path)
{
  WireSet raw;
  for (const CatalogueWire& w : path) {
    raw[{w.layer, w.horizontal, w.line}].emplace_back(std::min(w.lo, w.hi),
                                                      std::max(w.lo, w.hi));
  }
  WireSet out;
  for (auto& [key, intervals] : raw) {
    Intervals merged = mergedIntervals(std::move(intervals));
    if (!merged.empty()) {
      out[key] = std::move(merged);
    }
  }
  return out;
}

// Positive-length intersection of two interval lists.
Intervals intersection(const Intervals& a, const Intervals& b)
{
  Intervals out;
  size_t i = 0;
  size_t j = 0;
  while (i < a.size() && j < b.size()) {
    const GeoCoord lo = std::max(a[i].first, b[j].first);
    const GeoCoord hi = std::min(a[i].second, b[j].second);
    if (lo < hi) {
      out.emplace_back(lo, hi);
    }
    if (a[i].second < b[j].second) {
      i++;
    } else {
      j++;
    }
  }
  return out;
}

WireSet intersection(const WireSet& a, const WireSet& b)
{
  WireSet out;
  for (const auto& [key, intervals] : a) {
    const auto it = b.find(key);
    if (it == b.end()) {
      continue;
    }
    Intervals common = intersection(intervals, it->second);
    if (!common.empty()) {
      out[key] = std::move(common);
    }
  }
  return out;
}

int64_t length(const WireSet& s)
{
  int64_t total = 0;
  for (const auto& [key, intervals] : s) {
    for (const auto& [lo, hi] : intervals) {
      total += int64_t{hi} - lo;
    }
  }
  return total;
}

}  // namespace

CatalogueSharing catalogueSharing(
    const std::vector<std::vector<CatalogueWire>>& paths)
{
  CatalogueSharing out;
  std::vector<WireSet> sets;
  sets.reserve(paths.size());
  WireSet all;
  for (const std::vector<CatalogueWire>& path : paths) {
    sets.push_back(wireSet(path));
    out.wire += length(sets.back());
    for (const auto& [key, intervals] : sets.back()) {
      Intervals& u = all[key];
      u.insert(u.end(), intervals.begin(), intervals.end());
    }
  }
  for (auto& [key, intervals] : all) {
    intervals = mergedIntervals(std::move(intervals));
  }
  out.distinct_wire = length(all);
  for (size_t i = 0; i < sets.size(); i++) {
    for (size_t j = i + 1; j < sets.size(); j++) {
      out.pairs++;
      out.sharing_pairs += intersection(sets[i], sets[j]).empty() ? 0 : 1;
    }
  }
  if (sets.size() >= 2) {
    WireSet common = sets[0];
    for (size_t i = 1; i < sets.size() && !common.empty(); i++) {
      common = intersection(common, sets[i]);
    }
    out.common = !common.empty();
  }
  return out;
}

}  // namespace drt::br
