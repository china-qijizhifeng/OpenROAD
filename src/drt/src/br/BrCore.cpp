// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrCore.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <vector>

namespace drt::br {

Cost PathSearch::futureCost(const int z, const Coord x, const Coord y) const
{
  if (req_->future_cost != nullptr) {
    return req_->future_cost->query(z, x, y);
  }
  if (req_->target_boxes.empty()) {
    return 0;
  }
  Cost best = kInfCost;
  for (const TargetBox& b : req_->target_boxes) {
    const Cost dx = std::max({0, b.xlo - x, x - b.xhi});
    const Cost dy = std::max({0, b.ylo - y, y - b.yhi});
    best = std::min(best,
                    (dx + dy) * min_wire_rate_ + viaBound(z, b)
                        + (req_->exact_future_cost ? b.cost : 0));
  }
  return best;
}

Cost PathSearch::viaBound(const int z, const TargetBox& box) const
{
  if (z < box.zlo) {
    return req_->exact_future_cost ? via_prefix_[box.zlo] - via_prefix_[z]
                                   : (box.zlo - z) * min_via_cost_;
  }
  if (z > box.zhi) {
    return req_->exact_future_cost ? via_prefix_[z] - via_prefix_[box.zhi]
                                   : (z - box.zhi) * min_via_cost_;
  }
  return 0;
}

Cost PathSearch::spreadCost(const int z, const int i, const int j) const
{
  if (req_->spread_cost == 0 || req_->spread_tiles == nullptr) {
    return 0;
  }
  const TrackLayer& layer = graph_.layers[z];
  const int t = layer.track_tile[i];
  const int p = layer.pos_tile[j];
  const int tx = layer.horizontal ? p : t;
  const int ty = layer.horizontal ? t : p;
  if ((*req_->spread_tiles)[z][static_cast<size_t>(tx) * graph_.tiles_y + ty]
      == 0) {
    return 0;
  }
  // A neighbouring track is occupied where its wire vertex is not legal
  // (fast grid, not counted as a legality query of the path).
  for (const int di : {-1, 1}) {
    const int i2 = i + di;
    if (i2 >= 0 && i2 < layer.numTracks()
        && reader_.level(z, i2, j, kWire, req_->wire_type) != kLegal) {
      return req_->spread_cost;
    }
  }
  return 0;
}

Cost PathSearch::wireCost(const int z,
                          const int i,
                          const int j,
                          const int j2) const
{
  const TrackLayer& layer = graph_.layers[z];
  const Cost length = std::abs(layer.positions[j2] - layer.positions[j]);
  const Cost spread = spreadCost(z, i, j2);
  if (req_->discounts.empty()) {
    return length * graph_.wire_unit + spread;
  }
  size_t li = 0;
  const bool discounted = local(z, i, std::min(j, j2), li)
                          && (labels_[li].marks & kMarkDiscount) != 0;
  return length * (discounted ? req_->discount_rate : graph_.wire_unit)
         + spread;
}

Cost PathSearch::vertexPenalty(const At a, const Kind k) const
{
  if (level(a, k) == kLegal) {
    return 0;
  }
  return req_->ripup_cost
         * (1 + static_cast<Cost>(grid_.history(graph_.vertex(a.z, a.i, a.j))));
}

bool PathSearch::enter(const At a, const Kind k) const
{
  const RipupLevel l = level(a, k);
  return l == kLegal || (req_->ripup && l <= req_->max_ripup_level);
}

void PathSearch::push(const Item& item)
{
  heap_.push_back(item);
  std::push_heap(heap_.begin(), heap_.end(), std::greater<>());
  res_->heap_pushes++;
}

bool PathSearch::relax(const size_t from_v,
                       const size_t to_v,
                       const Cost g,
                       const Parent parent,
                       const int8_t dir)
{
  (void) from_v;
  int z, i, j;
  graph_.decode(to_v, z, i, j);
  size_t li = 0;
  if (!local(z, i, j, li)) {
    return false;
  }
  if (g >= labels_[li].dist) {
    return false;
  }
  labels_[li].dist = g;
  labels_[li].parent = parent;
  const auto [x, y] = graph_.layers[z].point(i, j);
  push(Item{g + futureCost(z, x, y), g, to_v, dir, 0});
  return true;
}

void PathSearch::settle(const int z,
                        const int i,
                        const int j,
                        const size_t v,
                        const Cost g)
{
  res_->settled++;
  size_t li = 0;
  local(z, i, j, li);
  const bool is_source = (labels_[li].marks & kMarkSource) != 0;
  if ((labels_[li].marks & kMarkTarget) != 0) {
    const auto& t = req_->targets[labels_[li].target_index];
    push(Item{g + t.cost, g, v, 2, labels_[li].target_index});
  }
  const TrackLayer& layer = graph_.layers[z];
  const Region& region = *req_->region;
  const At av{z, i, j};
  // Jogs (non-preferred direction, adjacent tracks).
  if (layer.jogs && (is_source || enter(av, kJog))) {
    for (const int di : {-1, 1}) {
      const int i2 = i + di;
      if (i2 < 0 || i2 >= layer.numTracks()
          || !region.contains(layer, z, i2, j)) {
        continue;
      }
      const size_t w = graph_.vertex(z, i2, j);
      const At aw{z, i2, j};
      const At edge = di > 0 ? av : aw;
      if (!enter(aw, kJog) || !enter(edge, kJogEdge)) {
        continue;
      }
      const Cost c
          = layer.jog_factor * std::abs(layer.tracks[i2] - layer.tracks[i])
            + vertexPenalty(aw, kJog) + vertexPenalty(edge, kJogEdge)
            + (is_source ? 0 : vertexPenalty(av, kJog));
      relax(v, w, g + c, di > 0 ? kFromTrackBelow : kFromTrackAbove, 0);
    }
  }
  // Via up.
  if (z + 1 < graph_.numLayers()) {
    const int i2 = layer.up_track[j];
    const int j2 = layer.up_pos[i];
    if (i2 >= 0 && j2 >= 0
        && region.contains(graph_.layers[z + 1], z + 1, i2, j2)) {
      const size_t w = graph_.vertex(z + 1, i2, j2);
      const At aw{z + 1, i2, j2};
      if (enter(av, kViaUpPad) && enter(av, kCutUp) && enter(aw, kViaDownPad)) {
        const Cost c = graph_.via_cost[z] + vertexPenalty(av, kViaUpPad)
                       + vertexPenalty(av, kCutUp)
                       + vertexPenalty(aw, kViaDownPad);
        relax(v, w, g + c, kFromLayerBelow, 0);
      }
    }
  }
  // Via down.
  if (z > 0) {
    const int i2 = layer.down_track[j];
    const int j2 = layer.down_pos[i];
    if (i2 >= 0 && j2 >= 0
        && region.contains(graph_.layers[z - 1], z - 1, i2, j2)) {
      const size_t w = graph_.vertex(z - 1, i2, j2);
      const At aw{z - 1, i2, j2};
      if (enter(aw, kViaUpPad) && enter(aw, kCutUp) && enter(av, kViaDownPad)) {
        const Cost c = graph_.via_cost[z - 1] + vertexPenalty(aw, kViaUpPad)
                       + vertexPenalty(aw, kCutUp)
                       + vertexPenalty(av, kViaDownPad);
        relax(v, w, g + c, kFromLayerAbove, 0);
      }
    }
  }
}

void PathSearch::setupWindow(const SearchRequest& request)
{
  const Region& region = *request.region;
  int bx_lo = region.bx_lo;
  int by_lo = region.by_lo;
  int bx_hi = region.bx_hi;
  int by_hi = region.by_hi;
  auto include = [&](size_t v) {
    int z, i, j;
    graph_.decode(v, z, i, j);
    const TrackLayer& l = graph_.layers[z];
    const int t = l.track_tile[i];
    const int p = l.pos_tile[j];
    const int tx = l.horizontal ? p : t;
    const int ty = l.horizontal ? t : p;
    bx_lo = std::min(bx_lo, tx);
    bx_hi = std::max(bx_hi, tx);
    by_lo = std::min(by_lo, ty);
    by_hi = std::max(by_hi, ty);
  };
  for (const auto& [v, c] : request.sources) {
    include(v);
  }
  for (const auto& t : request.targets) {
    include(t.vertex);
  }
  window_.assign(graph_.numLayers(), Window{});
  size_t total = 0;
  for (int z = 0; z < graph_.numLayers(); z++) {
    const TrackLayer& l = graph_.layers[z];
    const int tlo = l.horizontal ? by_lo : bx_lo;
    const int thi = l.horizontal ? by_hi : bx_hi;
    const int plo = l.horizontal ? bx_lo : by_lo;
    const int phi = l.horizontal ? bx_hi : by_hi;
    Window& w = window_[z];
    w.i0 = static_cast<int>(
        std::lower_bound(l.track_tile.begin(), l.track_tile.end(), tlo)
        - l.track_tile.begin());
    const int i1 = static_cast<int>(
        std::upper_bound(l.track_tile.begin(), l.track_tile.end(), thi)
        - l.track_tile.begin());
    w.j0 = static_cast<int>(
        std::lower_bound(l.pos_tile.begin(), l.pos_tile.end(), plo)
        - l.pos_tile.begin());
    const int j1 = static_cast<int>(
        std::upper_bound(l.pos_tile.begin(), l.pos_tile.end(), phi)
        - l.pos_tile.begin());
    w.ni = std::max(0, i1 - w.i0);
    w.nj = std::max(0, j1 - w.j0);
    w.wbase = total;
    total += static_cast<size_t>(w.ni) * w.nj;
  }
  epoch_ = labels_.reset(total, request.label_limit);
  for (const DiscountRun& run : request.discounts) {
    for (int j = run.jlo; j < run.jhi; j++) {
      size_t li = 0;
      if (local(run.z, run.i, j, li)) {
        labels_[li].marks |= kMarkDiscount;
      }
    }
  }
}

SearchResult PathSearch::run(const SearchRequest& request)
{
  SearchResult result;
  req_ = &request;
  res_ = &result;
  counters_ = {};
  reader_.reset(&grid_, request.overlay);
  min_via_cost_ = graph_.minViaCost();
  via_prefix_.assign(graph_.numLayers(), 0);
  for (int z = 1; z < graph_.numLayers(); z++) {
    via_prefix_[z] = via_prefix_[z - 1] + graph_.via_cost[z - 1];
  }
  min_wire_rate_ = request.discounts.empty()
                       ? graph_.wire_unit
                       : std::min(graph_.wire_unit, request.discount_rate);
  heap_.clear();
  setupWindow(request);

  for (size_t t = 0; t < request.targets.size(); t++) {
    int z, i, j;
    graph_.decode(request.targets[t].vertex, z, i, j);
    size_t li = 0;
    if (!local(z, i, j, li)) {
      continue;
    }
    if ((labels_[li].marks & kMarkTarget) != 0
        && request.targets[labels_[li].target_index].cost
               <= request.targets[t].cost) {
      continue;
    }
    labels_[li].marks |= kMarkTarget;
    labels_[li].target_index = static_cast<int32_t>(t);
  }
  for (const auto& [v, c] : request.sources) {
    int z, i, j;
    graph_.decode(v, z, i, j);
    size_t li = 0;
    if (!local(z, i, j, li)) {
      continue;
    }
    labels_[li].marks |= kMarkSource;
    if (labels_[li].dist <= c) {
      continue;
    }
    labels_[li].dist = c;
    labels_[li].parent = kSource;
    const auto [x, y] = graph_.layers[z].point(i, j);
    push(Item{c + futureCost(z, x, y), c, v, 0, 0});
  }

  while (!heap_.empty()) {
    if (labels_.overflow()) {
      result.aborted = true;
      break;
    }
    std::pop_heap(heap_.begin(), heap_.end(), std::greater<>());
    const Item item = heap_.back();
    heap_.pop_back();
    result.heap_pops++;
    if (item.dir == 2) {
      result.found = true;
      result.cost = item.key;
      result.target_id = request.targets[item.aux].id;
      result.target_vertex = item.v;
      break;
    }
    int z, i, j;
    graph_.decode(item.v, z, i, j);
    size_t li = 0;
    local(z, i, j, li);
    if (item.g > labels_[li].dist) {
      continue;
    }
    const TrackLayer& layer = graph_.layers[z];
    const Region& region = *request.region;
    settle(z, i, j, item.v, item.g);
    const bool is_source = (labels_[li].marks & kMarkSource) != 0;
    if (item.dir == 0) {
      // Vertex label: start runs in both preferred directions.
      if (!is_source && !enter(At{z, i, j}, kWire)) {
        continue;
      }
      for (const int d : {-1, 1}) {
        const int j2 = j + d;
        if (j2 < 0 || j2 >= layer.numPositions()
            || !region.contains(layer, z, i, j2)) {
          continue;
        }
        const size_t w = graph_.vertex(z, i, j2);
        const At aw{z, i, j2};
        const At edge{z, i, std::min(j, j2)};
        if (!enter(aw, kWire) || !enter(edge, kWireEdge)) {
          continue;
        }
        const Cost c = wireCost(z, i, j, j2) + vertexPenalty(aw, kWire)
                       + vertexPenalty(edge, kWireEdge);
        relax(item.v,
              w,
              item.g + c,
              d > 0 ? kFromLower : kFromHigher,
              static_cast<int8_t>(d));
      }
      continue;
    }
    // Run label: continue along the track. While the reduced cost stays
    // zero, the run is labelled without heap operations.
    const int d = item.dir;
    int cj = j;
    Cost g = item.g;
    while (true) {
      const int j2 = cj + d;
      if (j2 < 0 || j2 >= layer.numPositions()
          || !region.contains(layer, z, i, j2)) {
        break;
      }
      const size_t w = graph_.vertex(z, i, j2);
      const At aw{z, i, j2};
      const At edge{z, i, std::min(cj, j2)};
      if (!enter(aw, kWire) || !enter(edge, kWireEdge)) {
        break;
      }
      const Cost ng = g + wireCost(z, i, cj, j2) + vertexPenalty(aw, kWire)
                      + vertexPenalty(edge, kWireEdge);
      size_t lw = 0;
      if (!local(z, i, j2, lw)) {
        break;
      }
      if (ng >= labels_[lw].dist) {
        break;
      }
      labels_[lw].dist = ng;
      labels_[lw].parent = d > 0 ? kFromLower : kFromHigher;
      const auto [x, y] = layer.point(i, j2);
      const Cost key = ng + futureCost(z, x, y);
      push(Item{key, ng, w, static_cast<int8_t>(d), 0});
      break;
    }
  }
  result.fast_grid_queries = counters_.fast_grid;
  result.checker_queries = counters_.checker;
  if (labels_.overflow()) {
    result.aborted = true;
    result.found = false;
  }
  if (!result.found) {
    return result;
  }
  // Backtrack.
  size_t v = result.target_vertex;
  while (true) {
    result.path.push_back(v);
    int z, i, j;
    graph_.decode(v, z, i, j);
    size_t li = 0;
    local(z, i, j, li);
    const TrackLayer& l = graph_.layers[z];
    switch (static_cast<Parent>(labels_[li].parent)) {
      case kSource:
        return result;
      case kFromLower:
        v = graph_.vertex(z, i, j - 1);
        break;
      case kFromHigher:
        v = graph_.vertex(z, i, j + 1);
        break;
      case kFromTrackBelow:
        v = graph_.vertex(z, i - 1, j);
        break;
      case kFromTrackAbove:
        v = graph_.vertex(z, i + 1, j);
        break;
      case kFromLayerBelow:
        v = graph_.vertex(z - 1, l.down_track[j], l.down_pos[i]);
        break;
      case kFromLayerAbove:
        v = graph_.vertex(z + 1, l.up_track[j], l.up_pos[i]);
        break;
    }
  }
}

}  // namespace drt::br
