// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Interval-based shortest path search (Sec. 4.1, Algorithm 4), see BrCore.h.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <vector>

#include "br/BrCore.h"

namespace drt::br {

Cost IntervalSearch::futureCost(const int z, const int i, const int j) const
{
  if (req_->target_boxes.empty()) {
    return 0;
  }
  const auto [x, y] = graph_.layers[z].point(i, j);
  Cost best = kInfCost;
  for (const TargetBox& b : req_->target_boxes) {
    const Cost dx = std::max({0, b.xlo - x, x - b.xhi});
    const Cost dy = std::max({0, b.ylo - y, y - b.yhi});
    const Cost dz = z < b.zlo ? b.zlo - z : (z > b.zhi ? z - b.zhi : 0);
    best = std::min(best, dx + dy + dz * min_via_cost_);
  }
  return best;
}

bool IntervalSearch::usableZero(const size_t v, const Kind k) const
{
  return grid_.usable(v, k, req_->net);
}

bool IntervalSearch::enter(const size_t v, const Kind k) const
{
  return grid_.usable(v, k, req_->net)
         || (req_->ripup && grid_.ripupable(v, k, req_->net));
}

Cost IntervalSearch::penalty(const size_t v, const Kind k) const
{
  if (grid_.usable(v, k, req_->net)) {
    return 0;
  }
  return req_->ripup_cost * (1 + static_cast<Cost>(grid_.history(v)));
}

bool IntervalSearch::inWindow(const int z,
                              const int i,
                              const int j,
                              size_t& li) const
{
  const Window& w = window_[z];
  const int di = i - w.i0;
  const int dj = j - w.j0;
  if (di < 0 || di >= w.ni || dj < 0 || dj >= w.nj) {
    return false;
  }
  li = w.wbase + static_cast<size_t>(di) * w.nj + dj;
  return true;
}

int IntervalSearch::jogSignature(const int z,
                                 const int i,
                                 const int j,
                                 const size_t v) const
{
  const TrackLayer& layer = graph_.layers[z];
  if (!layer.jogs || !usableZero(v, kJog)) {
    return 0;
  }
  int sig = 0;
  for (const int di : {-1, 1}) {
    const int i2 = i + di;
    size_t li;
    if (i2 < 0 || i2 >= layer.numTracks() || !inWindow(z, i2, j, li)
        || !req_->region->contains(layer, z, i2, j)) {
      continue;
    }
    const size_t w = graph_.vertex(z, i2, j);
    const size_t edge = di > 0 ? v : w;
    if (usableZero(w, kJog) && usableZero(edge, kJogEdge)) {
      sig |= di > 0 ? 2 : 1;
    }
  }
  return sig;
}

int IntervalSearch::intervalOf(const int z, const int i, const int j)
{
  size_t li;
  if (!inWindow(z, i, j, li)) {
    return -1;
  }
  if (stamp_[li] == epoch_) {
    return interval_id_[li];
  }
  return buildInterval(z, i, j);
}

int IntervalSearch::buildInterval(const int z, const int i, const int j)
{
  const TrackLayer& layer = graph_.layers[z];
  const Region& region = *req_->region;
  // A vertex may join an interval if it is a plain wire vertex: in the
  // region, unassigned, wire-usable without penalty and not a target.
  auto plain = [&](int jj, size_t& li) {
    if (!inWindow(z, i, jj, li) || stamp_[li] == epoch_
        || target_stamp_[li] == epoch_ || !region.contains(layer, z, i, jj)) {
      return false;
    }
    return usableZero(graph_.vertex(z, i, jj), kWire);
  };
  // Slope class of pi on (a, a + 1): +1 / -1 if it changes by exactly the
  // distance (property (i)), 0 otherwise.
  auto slope = [&](int a) {
    const Cost dpi = futureCost(z, i, a + 1) - futureCost(z, i, a);
    const Cost dist = layer.positions[a + 1] - layer.positions[a];
    if (dpi == dist) {
      return 1;
    }
    if (dpi == -dist) {
      return -1;
    }
    return 0;
  };
  int ja = j;
  int jb = j;
  int cls = 0;
  size_t li0;
  inWindow(z, i, j, li0);
  if (plain(j, li0)) {
    const int sig = jogSignature(z, i, j, graph_.vertex(z, i, j));
    auto joins = [&](int from, int to) {
      size_t li;
      if (!plain(to, li)) {
        return false;
      }
      const int a = std::min(from, to);
      if (!usableZero(graph_.vertex(z, i, a), kWireEdge)) {
        return false;
      }
      if (jogSignature(z, i, to, graph_.vertex(z, i, to)) != sig) {
        return false;
      }
      const int c = slope(a);
      if (c == 0 || (cls != 0 && c != cls)) {
        return false;
      }
      cls = c;
      return true;
    };
    while (jb + 1 < layer.numPositions() && joins(jb, jb + 1)) {
      jb++;
    }
    while (ja - 1 >= 0 && joins(ja, ja - 1)) {
      ja--;
    }
  }
  Interval iv;
  iv.z = z;
  iv.i = i;
  iv.ja = ja;
  iv.jb = jb;
  // pi increases with j (cls +1): the target end is ja.
  iv.dir = cls >= 0 ? 1 : -1;
  const int id = static_cast<int>(intervals_.size());
  intervals_.push_back(std::move(iv));
  for (int jj = ja; jj <= jb; jj++) {
    size_t li;
    inWindow(z, i, jj, li);
    stamp_[li] = epoch_;
    interval_id_[li] = id;
  }
  num_intervals_++;
  return id;
}

Cost IntervalSearch::value(const Interval& iv, const int k) const
{
  const TrackLayer& layer = graph_.layers[iv.z];
  Cost best = kInfCost;
  const Coord pk = layer.positions[iv.jOf(k)];
  for (const Label& l : iv.labels) {
    // Flat toward the target end, slope 2 away from it.
    const Cost d
        = k <= l.k ? l.delta
                   : l.delta + 2 * std::abs(pk - layer.positions[iv.jOf(l.k)]);
    best = std::min(best, d);
  }
  return best;
}

const IntervalSearch::Label* IntervalSearch::argLabel(const Interval& iv,
                                                      const int k) const
{
  const TrackLayer& layer = graph_.layers[iv.z];
  const Coord pk = layer.positions[iv.jOf(k)];
  const Label* best = nullptr;
  Cost best_d = kInfCost;
  for (const Label& l : iv.labels) {
    const Cost d
        = k <= l.k ? l.delta
                   : l.delta + 2 * std::abs(pk - layer.positions[iv.jOf(l.k)]);
    if (d < best_d) {
      best_d = d;
      best = &l;
    }
  }
  return best;
}

void IntervalSearch::pushEvent(const int interval)
{
  const Interval& iv = intervals_[interval];
  if (iv.frontier >= iv.n()) {
    return;
  }
  const Cost key = value(iv, iv.frontier);
  if (key >= kInfCost) {
    return;
  }
  heap_.push_back(Event{key, interval, 0});
  std::push_heap(heap_.begin(), heap_.end(), std::greater<>());
  res_->heap_pushes++;
}

void IntervalSearch::addLabel(const int interval,
                              const int k,
                              const Cost delta,
                              const int pred_interval,
                              const int pred_k)
{
  Interval& iv = intervals_[interval];
  // Settled vertices already have d <= the current level <= delta.
  if (k < iv.frontier || value(iv, k) <= delta) {
    return;
  }
  const TrackLayer& layer = graph_.layers[iv.z];
  const Coord pk = layer.positions[iv.jOf(k)];
  // Drop labels the new one strictly dominates.
  std::erase_if(iv.labels, [&](const Label& l) {
    const Cost d
        = l.k <= k ? delta
                   : delta + 2 * std::abs(layer.positions[iv.jOf(l.k)] - pk);
    return d < l.delta;
  });
  iv.labels.push_back(Label{k, delta, pred_interval, pred_k});
  num_labels_++;
  pushEvent(interval);
}

void IntervalSearch::process(const int id)
{
  // J_I(delta): the settled frontier moves over the flat run at delta.
  const Cost delta = value(intervals_[id], intervals_[id].frontier);
  int g = intervals_[id].frontier;
  const int n = intervals_[id].n();
  while (g < n && value(intervals_[id], g) == delta) {
    g++;
  }
  const int f = intervals_[id].frontier;
  intervals_[id].frontier = g;
  const int z = intervals_[id].z;
  const int i = intervals_[id].i;
  const int ja = intervals_[id].ja;
  const int jb = intervals_[id].jb;
  const int dir = intervals_[id].dir;
  const bool single = n == 1;
  const TrackLayer& layer = graph_.layers[z];
  const Region& region = *req_->region;
  auto j_of = [&](int k) { return dir > 0 ? ja + k : jb - k; };

  // Jog candidates kept per target interval (one label at max pi, line 14).
  struct Pending
  {
    int interval;
    int k;
    Cost delta;
    int from_k;
    Cost pi;
  };
  std::vector<Pending> jog_best;

  for (int k = f; k < g; k++) {
    const int j = j_of(k);
    const size_t v = graph_.vertex(z, i, j);
    size_t li;
    inWindow(z, i, j, li);
    res_->settled++;
    const bool is_source = source_stamp_[li] == epoch_;
    const Cost pi_v = futureCost(z, i, j);
    if (target_stamp_[li] == epoch_) {
      const auto& t = req_->targets[target_index_[li]];
      heap_.push_back(Event{delta + t.cost, -1, target_index_[li]});
      std::push_heap(heap_.begin(), heap_.end(), std::greater<>());
      res_->heap_pushes++;
    }
    auto relax_to = [&](int z2, int i2, int j2, Cost cost) {
      const int iv2 = intervalOf(z2, i2, j2);
      if (iv2 < 0) {
        return;
      }
      const Cost d2 = delta + cost - pi_v + futureCost(z2, i2, j2);
      addLabel(iv2, intervals_[iv2].kOf(j2), d2, id, k);
    };
    // Vias, per vertex.
    if (z + 1 < graph_.numLayers()) {
      const int i2 = layer.up_track[j];
      const int j2 = layer.up_pos[i];
      if (i2 >= 0 && j2 >= 0
          && region.contains(graph_.layers[z + 1], z + 1, i2, j2)) {
        const size_t w = graph_.vertex(z + 1, i2, j2);
        if (enter(v, kViaUpPad) && enter(v, kCutUp) && enter(w, kViaDownPad)) {
          relax_to(z + 1,
                   i2,
                   j2,
                   graph_.via_cost[z] + penalty(v, kViaUpPad)
                       + penalty(v, kCutUp) + penalty(w, kViaDownPad));
        }
      }
    }
    if (z > 0) {
      const int i2 = layer.down_track[j];
      const int j2 = layer.down_pos[i];
      if (i2 >= 0 && j2 >= 0
          && region.contains(graph_.layers[z - 1], z - 1, i2, j2)) {
        const size_t w = graph_.vertex(z - 1, i2, j2);
        if (enter(w, kViaUpPad) && enter(w, kCutUp) && enter(v, kViaDownPad)) {
          relax_to(z - 1,
                   i2,
                   j2,
                   graph_.via_cost[z - 1] + penalty(w, kViaUpPad)
                       + penalty(w, kCutUp) + penalty(v, kViaDownPad));
        }
      }
    }
    // Jogs.
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
        const Cost pen = penalty(w, kJog) + penalty(edge, kJogEdge)
                         + (is_source ? 0 : penalty(v, kJog));
        const int iv2 = intervalOf(z, i2, j);
        if (iv2 < 0) {
          continue;
        }
        const Cost cost
            = layer.jog_factor * std::abs(layer.tracks[i2] - layer.tracks[i])
              + pen;
        const Cost d2 = delta + cost - pi_v + futureCost(z, i2, j);
        const int k2 = intervals_[iv2].kOf(j);
        // Plain interval pairs with the same orientation: the landing of the
        // vertex with maximum pi dominates the others (Algorithm 4, line 14).
        const bool shortcut = !single && pen == 0 && intervals_[iv2].n() > 1
                              && intervals_[iv2].dir == dir;
        if (!shortcut) {
          addLabel(iv2, k2, d2, id, k);
          continue;
        }
        auto it = std::ranges::find_if(
            jog_best, [&](const Pending& p) { return p.interval == iv2; });
        if (it == jog_best.end()) {
          jog_best.push_back(Pending{iv2, k2, d2, k, pi_v});
        } else if (pi_v > it->pi) {
          *it = Pending{iv2, k2, d2, k, pi_v};
        }
      }
    }
    // Preferred-direction edges leaving the interval at its ends.
    for (const int dj : {-1, 1}) {
      if ((dj < 0 && j != ja) || (dj > 0 && j != jb)) {
        continue;
      }
      const int j2 = j + dj;
      if (j2 < 0 || j2 >= layer.numPositions()
          || !region.contains(layer, z, i, j2)) {
        continue;
      }
      if (!is_source && !enter(v, kWire)) {
        continue;
      }
      const size_t w = graph_.vertex(z, i, j2);
      const size_t edge = dj > 0 ? v : w;
      if (!enter(w, kWire) || !enter(edge, kWireEdge)) {
        continue;
      }
      relax_to(z,
               i,
               j2,
               std::abs(layer.positions[j2] - layer.positions[j])
                   + penalty(w, kWire) + penalty(edge, kWireEdge));
    }
  }
  for (const Pending& p : jog_best) {
    addLabel(p.interval, p.k, p.delta, id, p.from_k);
  }
  pushEvent(id);
}

SearchResult IntervalSearch::run(const SearchRequest& request)
{
  SearchResult result;
  req_ = &request;
  res_ = &result;
  min_via_cost_ = graph_.minViaCost();
  heap_.clear();
  intervals_.clear();
  // Search window over the region and the terminals (as in PathSearch).
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
  if (stamp_.size() < total) {
    interval_id_.resize(total);
    stamp_.resize(total, 0);
    source_stamp_.resize(total, 0);
    target_stamp_.resize(total, 0);
    target_index_.resize(total);
  }
  if (++epoch_ == 0) {
    std::ranges::fill(stamp_, 0);
    std::ranges::fill(source_stamp_, 0);
    std::ranges::fill(target_stamp_, 0);
    epoch_ = 1;
  }
  for (size_t t = 0; t < request.targets.size(); t++) {
    int z, i, j;
    graph_.decode(request.targets[t].vertex, z, i, j);
    size_t li;
    if (!inWindow(z, i, j, li)) {
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
    if (inWindow(z, i, j, li)) {
      source_stamp_[li] = epoch_;
    }
  }
  // Algorithm 4, lines 1-4 (with source costs).
  for (const auto& [v, c] : request.sources) {
    int z, i, j;
    graph_.decode(v, z, i, j);
    const int id = intervalOf(z, i, j);
    if (id >= 0) {
      addLabel(id, intervals_[id].kOf(j), c + futureCost(z, i, j), -1, -1);
    }
  }
  while (!heap_.empty()) {
    std::pop_heap(heap_.begin(), heap_.end(), std::greater<>());
    const Event e = heap_.back();
    heap_.pop_back();
    result.heap_pops++;
    if (e.interval < 0) {
      result.found = true;
      result.cost = e.key;
      result.target_id = request.targets[e.aux].id;
      result.target_vertex = request.targets[e.aux].vertex;
      break;
    }
    const Interval& iv = intervals_[e.interval];
    if (iv.frontier >= iv.n() || value(iv, iv.frontier) != e.key) {
      continue;
    }
    process(e.interval);
  }
  if (!result.found) {
    return result;
  }
  // Backtrack through the labels.
  int z, i, j;
  graph_.decode(result.target_vertex, z, i, j);
  size_t li;
  inWindow(z, i, j, li);
  int id = interval_id_[li];
  int k = intervals_[id].kOf(j);
  while (true) {
    const Interval& iv = intervals_[id];
    const Label* l = argLabel(iv, k);
    const int step = l->k >= k ? 1 : -1;
    for (int kk = k;; kk += step) {
      result.path.push_back(graph_.vertex(iv.z, iv.i, iv.jOf(kk)));
      if (kk == l->k) {
        break;
      }
    }
    if (l->pred_interval < 0) {
      break;
    }
    id = l->pred_interval;
    k = l->pred_k;
  }
  return result;
}

}  // namespace drt::br
