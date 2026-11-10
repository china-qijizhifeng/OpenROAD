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
  if (req_->target_boxes.empty()) {
    return 0;
  }
  Cost best = kInfCost;
  for (const TargetBox& b : req_->target_boxes) {
    const Cost dx = std::max({0, b.xlo - x, x - b.xhi});
    const Cost dy = std::max({0, b.ylo - y, y - b.yhi});
    const Cost dz = z < b.zlo ? b.zlo - z : (z > b.zhi ? z - b.zhi : 0);
    best = std::min(best, dx + dy + dz * min_via_cost_);
  }
  return best;
}

Cost PathSearch::vertexPenalty(const size_t v, const Kind k) const
{
  if (grid_.usable(v, k, req_->net)) {
    return 0;
  }
  return req_->ripup_cost * (1 + static_cast<Cost>(grid_.history(v)));
}

bool PathSearch::enter(const size_t v, const Kind k) const
{
  if (grid_.usable(v, k, req_->net)) {
    return true;
  }
  return req_->ripup && grid_.ripupable(v, k, req_->net);
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
  size_t li;
  if (!local(z, i, j, li)) {
    return false;
  }
  if (stamp_[li] == epoch_ && g >= dist_[li]) {
    return false;
  }
  stamp_[li] = epoch_;
  dist_[li] = g;
  parent_[li] = parent;
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
  size_t li;
  local(z, i, j, li);
  const bool is_source = source_stamp_[li] == epoch_;
  if (target_stamp_[li] == epoch_) {
    const auto& t = req_->targets[target_index_[li]];
    push(Item{g + t.cost, g, v, 2, target_index_[li]});
  }
  const TrackLayer& layer = graph_.layers[z];
  const Region& region = *req_->region;
  // Jogs (non-preferred direction, adjacent tracks).
  if (layer.jogs && (is_source || enter(v, kJog))) {
    for (const int di : {-1, 1}) {
      const int i2 = i + di;
      if (i2 < 0 || i2 >= layer.numTracks()
          || !region.contains(layer, z, i2, j)) {
        continue;
      }
      const size_t w = graph_.vertex(z, i2, j);
      const size_t edge = di > 0 ? v : w;
      if (!enter(w, kJog) || !enter(edge, kJogEdge)) {
        continue;
      }
      const Cost c
          = layer.jog_factor * std::abs(layer.tracks[i2] - layer.tracks[i])
            + vertexPenalty(w, kJog) + vertexPenalty(edge, kJogEdge)
            + (is_source ? 0 : vertexPenalty(v, kJog));
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
      if (enter(v, kViaUpPad) && enter(v, kCutUp) && enter(w, kViaDownPad)) {
        const Cost c = graph_.via_cost[z] + vertexPenalty(v, kViaUpPad)
                       + vertexPenalty(v, kCutUp)
                       + vertexPenalty(w, kViaDownPad);
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
      if (enter(w, kViaUpPad) && enter(w, kCutUp) && enter(v, kViaDownPad)) {
        const Cost c = graph_.via_cost[z - 1] + vertexPenalty(w, kViaUpPad)
                       + vertexPenalty(w, kCutUp)
                       + vertexPenalty(v, kViaDownPad);
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
  if (dist_.size() < total) {
    dist_.resize(total);
    parent_.resize(total);
    stamp_.resize(total, 0);
    source_stamp_.resize(total, 0);
    target_stamp_.resize(total, 0);
    target_index_.resize(total);
  }
  if (++epoch_ == 0) {
    std::fill(stamp_.begin(), stamp_.end(), 0);
    std::fill(source_stamp_.begin(), source_stamp_.end(), 0);
    std::fill(target_stamp_.begin(), target_stamp_.end(), 0);
    epoch_ = 1;
  }
}

SearchResult PathSearch::run(const SearchRequest& request)
{
  SearchResult result;
  req_ = &request;
  res_ = &result;
  min_via_cost_ = graph_.minViaCost();
  heap_.clear();
  setupWindow(request);

  for (size_t t = 0; t < request.targets.size(); t++) {
    int z, i, j;
    graph_.decode(request.targets[t].vertex, z, i, j);
    size_t li;
    if (!local(z, i, j, li)) {
      continue;
    }
    if (target_stamp_[li] == epoch_
        && request.targets[target_index_[li]].cost <= request.targets[t].cost) {
      continue;
    }
    target_stamp_[li] = epoch_;
    target_index_[li] = static_cast<int32_t>(t);
  }
  for (const auto& [v, c] : request.sources) {
    int z, i, j;
    graph_.decode(v, z, i, j);
    size_t li;
    if (!local(z, i, j, li)) {
      continue;
    }
    source_stamp_[li] = epoch_;
    if (stamp_[li] == epoch_ && dist_[li] <= c) {
      continue;
    }
    stamp_[li] = epoch_;
    dist_[li] = c;
    parent_[li] = kSource;
    const auto [x, y] = graph_.layers[z].point(i, j);
    push(Item{c + futureCost(z, x, y), c, v, 0, 0});
  }

  while (!heap_.empty()) {
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
    size_t li;
    local(z, i, j, li);
    if (item.g > dist_[li]) {
      continue;
    }
    const TrackLayer& layer = graph_.layers[z];
    const Region& region = *request.region;
    settle(z, i, j, item.v, item.g);
    const bool is_source = source_stamp_[li] == epoch_;
    if (item.dir == 0) {
      // Vertex label: start runs in both preferred directions.
      if (!is_source && !enter(item.v, kWire)) {
        continue;
      }
      for (const int d : {-1, 1}) {
        const int j2 = j + d;
        if (j2 < 0 || j2 >= layer.numPositions()
            || !region.contains(layer, z, i, j2)) {
          continue;
        }
        const size_t w = graph_.vertex(z, i, j2);
        const size_t edge = d > 0 ? item.v : w;
        if (!enter(w, kWire) || !enter(edge, kWireEdge)) {
          continue;
        }
        const Cost c = std::abs(layer.positions[j2] - layer.positions[j])
                       + vertexPenalty(w, kWire)
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
      const size_t edge = d > 0 ? graph_.vertex(z, i, cj) : w;
      if (!enter(w, kWire) || !enter(edge, kWireEdge)) {
        break;
      }
      const Cost ng = g + std::abs(layer.positions[j2] - layer.positions[cj])
                      + vertexPenalty(w, kWire)
                      + vertexPenalty(edge, kWireEdge);
      size_t lw;
      if (!local(z, i, j2, lw)) {
        break;
      }
      if (stamp_[lw] == epoch_ && ng >= dist_[lw]) {
        break;
      }
      stamp_[lw] = epoch_;
      dist_[lw] = ng;
      parent_[lw] = d > 0 ? kFromLower : kFromHigher;
      const auto [x, y] = layer.point(i, j2);
      const Cost key = ng + futureCost(z, x, y);
      push(Item{key, ng, w, static_cast<int8_t>(d), 0});
      break;
    }
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
    size_t li;
    local(z, i, j, li);
    const TrackLayer& l = graph_.layers[z];
    switch (static_cast<Parent>(parent_[li])) {
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
