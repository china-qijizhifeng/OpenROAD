// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Interval-based shortest path search (Sec. 4.1, Algorithm 4), see BrCore.h.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <utility>
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
    best = std::min(
        best,
        dx + dy + viaBound(z, b) + (req_->exact_future_cost ? b.cost : 0));
  }
  return best;
}

Cost IntervalSearch::viaBound(const int z, const TargetBox& box) const
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

bool IntervalSearch::usableZero(const At v, const Kind k) const
{
  return level(v, k) == kLegal;
}

bool IntervalSearch::enter(const At v, const Kind k) const
{
  const RipupLevel l = level(v, k);
  return l == kLegal || (req_->ripup && l <= req_->max_ripup_level);
}

Cost IntervalSearch::penalty(const At v, const Kind k) const
{
  if (level(v, k) == kLegal) {
    return 0;
  }
  return req_->ripup_cost
         * (1 + static_cast<Cost>(grid_.history(graph_.vertex(v.z, v.i, v.j))));
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

uint8_t IntervalSearch::computeTransitions(const int z,
                                           const int i,
                                           const int j) const
{
  const At v{z, i, j};
  const TrackLayer& layer = graph_.layers[z];
  uint8_t bits = 0;
  if (layer.jogs && enter(v, kJog)) {
    const bool v_zero = usableZero(v, kJog);
    for (const int di : {-1, 1}) {
      const int i2 = i + di;
      size_t li = 0;
      if (i2 < 0 || i2 >= layer.numTracks() || !inWindow(z, i2, j, li)
          || !req_->region->contains(layer, z, i2, j)) {
        continue;
      }
      const At w{z, i2, j};
      const At edge = di > 0 ? v : w;
      if (enter(w, kJog) && enter(edge, kJogEdge)) {
        bits |= di > 0 ? kTransJogAbove : kTransJogBelow;
        if (!v_zero || !usableZero(w, kJog) || !usableZero(edge, kJogEdge)) {
          bits |= kTransPenalty;  // cannot be compressed into an interval
        }
      }
    }
  }
  for (const int dz : {-1, 1}) {
    const int z2 = z + dz;
    if (z2 < 0 || z2 >= graph_.numLayers()) {
      continue;
    }
    const int i2 = dz > 0 ? layer.up_track[j] : layer.down_track[j];
    const int j2 = dz > 0 ? layer.up_pos[i] : layer.down_pos[i];
    size_t li = 0;
    if (i2 < 0 || j2 < 0 || !inWindow(z2, i2, j2, li)
        || !req_->region->contains(graph_.layers[z2], z2, i2, j2)) {
      continue;
    }
    const At w{z2, i2, j2};
    const At lower = dz > 0 ? v : w;
    const At upper = dz > 0 ? w : v;
    const RipupLevel pad = level(lower, kViaUpPad);
    const RipupLevel cut = level(lower, kCutUp);
    const RipupLevel top = level(upper, kViaDownPad);
    const RipupLevel worst = std::max({pad, cut, top});
    if (worst == kLegal || (req_->ripup && worst <= req_->max_ripup_level)) {
      bits |= dz > 0 ? kTransViaAbove : kTransViaBelow;
      if (worst != kLegal) {
        bits |= kTransPenalty;
      }
    }
  }
  return bits;
}

uint8_t IntervalSearch::transitions(const int z,
                                    const int i,
                                    const int j,
                                    VertexLabel& label) const
{
  if (label.transitions_stamp != epoch_) {
    label.transitions = computeTransitions(z, i, j);
    label.transitions_stamp = epoch_;
  }
  return label.transitions;
}

int IntervalSearch::intervalOf(const int z, const int i, const int j)
{
  size_t li = 0;
  if (!inWindow(z, i, j, li)) {
    return -1;
  }
  if (labels_[li].stamp == epoch_) {
    return labels_[li].interval_id;
  }
  return buildInterval(z, i, j);
}

int IntervalSearch::buildInterval(const int z, const int i, const int j)
{
  const TrackLayer& layer = graph_.layers[z];
  const Region& region = *req_->region;
  // A vertex may join an interval if it is a plain wire vertex: in the
  // region, unassigned, wire-usable without penalty and not a target.
  // Its transitions (property (ii)) are computed once and kept.
  const uint8_t property_mask
      = req_->interval_via_property ? uint8_t{0xff}
                                    : uint8_t{kTransJogBelow | kTransJogAbove
                                              | kTransPenalty};
  auto plain = [&](int jj, uint8_t& sig) {
    size_t li = 0;
    if (!inWindow(z, i, jj, li)) {
      return false;
    }
    VertexLabel& label = labels_[li];
    if (label.stamp == epoch_ || label.target_stamp == epoch_
        || !region.contains(layer, z, i, jj)) {
      return false;
    }
    const At v{z, i, jj};
    if (!usableZero(v, kWire)
        || (label.source_stamp == epoch_ && layer.jogs
            && !usableZero(v, kJog))) {
      return false;
    }
    const uint8_t bits = transitions(z, i, jj, label);
    sig = bits & property_mask;
    return (bits & kTransPenalty) == 0;
  };
  int ja = j;
  int jb = j;
  int cls = 0;
  const Cost pi0 = futureCost(z, i, j);
  Cost pi_a = pi0;
  Cost pi_b = pi0;
  uint8_t sig = 0;
  if (!plain(j, sig)) {
    ends_[kEndStartVertex] += 2;
  } else {
    IntervalEnd end = kEndTrack;
    auto joins = [&](int from, int to, Cost& pi_from) {
      uint8_t sig_to = 0;
      if (!plain(to, sig_to)) {
        end = kEndNeighbour;
        return false;
      }
      const int a = std::min(from, to);
      if (!usableZero(At{z, i, a}, kWireEdge)) {
        end = kEndWireEdge;
        return false;
      }
      if (sig_to != sig) {
        end = ((sig_to ^ sig) & (kTransJogBelow | kTransJogAbove)) != 0
                  ? kEndJogProperty
                  : kEndViaProperty;
        return false;
      }
      // Adjacent members must have the same slope, exactly +1 or -1.
      const Cost pi_to = futureCost(z, i, to);
      const Cost step
          = static_cast<Cost>(layer.positions[to]) - layer.positions[from];
      const Cost dpi = pi_to - pi_from;
      const int c = dpi == step ? 1 : (dpi == -step ? -1 : 0);
      if (c == 0 || (cls != 0 && c != cls)) {
        end = kEndFutureCostSlope;
        return false;
      }
      cls = c;
      pi_from = pi_to;
      return true;
    };
    end = kEndTrack;
    while (jb + 1 < layer.numPositions() && joins(jb, jb + 1, pi_b)) {
      jb++;
    }
    ends_[end]++;
    end = kEndTrack;
    while (ja > 0 && joins(ja, ja - 1, pi_a)) {
      ja--;
    }
    ends_[end]++;
  }
  interval_vertices_ += jb - ja + 1;
  Interval iv;
  iv.z = z;
  iv.i = i;
  iv.ja = ja;
  iv.jb = jb;
  // pi increases with j (cls +1): the target end is ja.
  iv.dir = cls >= 0 ? 1 : -1;
  // The checked slope determines pi everywhere in the interval.
  iv.future_at_ja = pi_a;
  const int id = static_cast<int>(intervals_.size());
  intervals_.push_back(std::move(iv));
  for (int jj = ja; jj <= jb; jj++) {
    size_t li = 0;
    inWindow(z, i, jj, li);
    labels_[li].stamp = epoch_;
    labels_[li].interval_id = id;
  }
  num_intervals_++;
  return id;
}

Cost IntervalSearch::value(const Interval& iv, const int k) const
{
  const Label* l = argLabel(iv, k);
  if (l == nullptr) {
    return kInfCost;
  }
  const auto& positions = graph_.layers[iv.z].positions;
  return l->delta
         + (k <= l->k ? 0
                      : 2
                            * std::abs(static_cast<Cost>(positions[iv.jOf(k)])
                                       - positions[iv.jOf(l->k)]));
}

const IntervalSearch::Label* IntervalSearch::argLabel(const Interval& iv,
                                                      const int k) const
{
  if (iv.labels.empty()) {
    return nullptr;
  }
  const auto& positions = graph_.layers[iv.z].positions;
  const int pos = static_cast<int>(
      std::lower_bound(iv.labels.begin(),
                       iv.labels.end(),
                       k,
                       [](const Label& l, int at) { return l.k < at; })
      - iv.labels.begin());
  const Label* best = pos < static_cast<int>(iv.labels.size())
                          ? &iv.labels[iv.suffix_best[pos]]
                          : nullptr;
  if (pos > 0) {
    const Label* left = &iv.labels[iv.prefix_best[pos - 1]];
    const Cost d = left->delta
                   + 2
                         * std::abs(static_cast<Cost>(positions[iv.jOf(k)])
                                    - positions[iv.jOf(left->k)]);
    if (best == nullptr || d < best->delta
        || (d == best->delta && left->order < best->order)) {
      best = left;
    }
  }
  return best;
}

void IntervalSearch::rebuildEnvelope(Interval& iv)
{
  const int n = static_cast<int>(iv.labels.size());
  iv.prefix_best.resize(n);
  iv.suffix_best.resize(n);
  const auto& positions = graph_.layers[iv.z].positions;
  auto intercept = [&](int a) {
    const Label& l = iv.labels[a];
    return l.delta - 2 * iv.dir * static_cast<Cost>(positions[iv.jOf(l.k)]);
  };
  for (int a = 0; a < n; a++) {
    const int prev = a > 0 ? iv.prefix_best[a - 1] : a;
    iv.prefix_best[a]
        = intercept(a) < intercept(prev)
                  || (intercept(a) == intercept(prev)
                      && iv.labels[a].order < iv.labels[prev].order)
              ? a
              : prev;
  }
  for (int a = n - 1; a >= 0; a--) {
    const int next = a + 1 < n ? iv.suffix_best[a + 1] : a;
    iv.suffix_best[a]
        = iv.labels[a].delta < iv.labels[next].delta
                  || (iv.labels[a].delta == iv.labels[next].delta
                      && iv.labels[a].order < iv.labels[next].order)
              ? a
              : next;
  }
}

int IntervalSearch::frontierEnd(const Interval& iv, const Cost delta) const
{
  // d is monotone on an interval; J_I(delta) is a contiguous flat run.
  int lo = iv.frontier;
  int hi = iv.n();
  while (lo < hi) {
    const int mid = lo + (hi - lo) / 2;
    if (value(iv, mid) == delta) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
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
  const int end = frontierEnd(iv, key);
  const Cost pi = iv.futureCost(graph_.layers[iv.z], iv.jOf(end - 1));
  heap_.push_back(Event{key, pi, interval, iv.frontier});
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
  auto at = std::lower_bound(
      iv.labels.begin(), iv.labels.end(), k, [](const Label& l, int pos) {
        return l.k < pos;
      });
  iv.labels.insert(at, Label{k, delta, pred_interval, pred_k, num_labels_});
  rebuildEnvelope(iv);
  num_labels_++;
  pushEvent(interval);
}

void IntervalSearch::process(const int id, const Cost delta)
{
  const int f = intervals_[id].frontier;
  const int g = frontierEnd(intervals_[id], delta);
  intervals_[id].frontier = g;
  const int z = intervals_[id].z;
  const int i = intervals_[id].i;
  const int ja = intervals_[id].ja;
  const int jb = intervals_[id].jb;
  const int dir = intervals_[id].dir;
  const bool single = ja == jb;
  const TrackLayer& layer = graph_.layers[z];
  const Region& region = *req_->region;
  auto j_of = [&](int k) { return dir > 0 ? ja + k : jb - k; };
  auto k_of = [&](int j) { return dir > 0 ? j - ja : jb - j; };
  res_->settled += g - f;
  // Members of a longer interval are plain: their transitions are cached
  // and carry no penalties.
  auto bits_of = [&](int k) {
    size_t li = 0;
    inWindow(z, i, j_of(k), li);
    return labels_[li].transitions;
  };
  const uint8_t signature = single ? 0 : bits_of(f);
  auto relax_to = [&](int from_k, int z2, int i2, int j2, Cost cost) {
    const Cost pi = intervals_[id].futureCost(layer, j_of(from_k));
    const int next_id = intervalOf(z2, i2, j2);
    if (next_id < 0) {
      return;
    }
    const Interval& next = intervals_[next_id];
    const Cost d2 = delta + cost - pi + next.futureCost(graph_.layers[z2], j2);
    addLabel(next_id, next.kOf(j2), d2, id, from_k);
  };

  // Targets are single-vertex intervals, so terminal costs are exact.
  if (single) {
    size_t li = 0;
    inWindow(z, i, ja, li);
    if (labels_[li].target_stamp == epoch_) {
      const auto& t = req_->targets[labels_[li].target_index];
      // delta includes pi(t) <= t.cost (the box cost is the cheapest
      // terminal cost), so the sink key stays >= delta.
      const Cost pi_t = intervals_[id].futureCost(layer, ja);
      heap_.push_back(
          Event{delta - pi_t + t.cost, 0, -1, labels_[li].target_index});
      std::push_heap(heap_.begin(), heap_.end(), std::greater<>());
      res_->heap_pushes++;
    }
  }

  // On perpendicular layers each crossing reaches a different track and
  // hence a different neighboring interval. Enumerate those adjacencies;
  // skip the entire run when the interval has no vias in that direction.
  for (const int dz : {-1, 1}) {
    const int z2 = z + dz;
    const uint8_t via_bit = dz > 0 ? kTransViaAbove : kTransViaBelow;
    if (z2 < 0 || z2 >= graph_.numLayers()
        || (!single && req_->interval_via_property && !(signature & via_bit))) {
      continue;
    }
    for (int k = f; k < g; k++) {
      const int j = j_of(k);
      const int i2 = dz > 0 ? layer.up_track[j] : layer.down_track[j];
      const int j2 = dz > 0 ? layer.up_pos[i] : layer.down_pos[i];
      if (!single) {
        if (bits_of(k) & via_bit) {
          relax_to(k, z2, i2, j2, graph_.via_cost[std::min(z, z2)]);
        }
        continue;
      }
      if (i2 < 0 || j2 < 0 || !region.contains(graph_.layers[z2], z2, i2, j2)) {
        continue;
      }
      const At v{z, i, j};
      const At w{z2, i2, j2};
      const At lower = dz > 0 ? v : w;
      const At upper = dz > 0 ? w : v;
      if (enter(lower, kViaUpPad) && enter(lower, kCutUp)
          && enter(upper, kViaDownPad)) {
        relax_to(k,
                 z2,
                 i2,
                 j2,
                 graph_.via_cost[std::min(z, z2)] + penalty(lower, kViaUpPad)
                     + penalty(lower, kCutUp) + penalty(upper, kViaDownPad));
      }
    }
  }

  // Algorithm 4, lines 13-17: process a neighboring jog interval once,
  // choosing the member of J_I(delta) with maximum pi. Interval membership
  // includes all jog availability; transitions carrying penalties are
  // singletons. No per-vertex jog labeling is needed on a clean run.
  if (layer.jogs) {
    for (const int di : {-1, 1}) {
      const int i2 = i + di;
      const uint8_t jog_bit = di > 0 ? kTransJogAbove : kTransJogBelow;
      if (i2 < 0 || i2 >= layer.numTracks()
          || (!single && !(signature & jog_bit))) {
        continue;
      }
      if (!single) {
        // Every member has the jog (property (ii)); one relaxation per
        // neighbouring interval, at its last member of J_I(delta).
        const Cost cost
            = layer.jog_factor * std::abs(layer.tracks[i2] - layer.tracks[i]);
        for (int k = f; k < g;) {
          const int next_id = intervalOf(z, i2, j_of(k));
          const Interval& next = intervals_[next_id];
          const int last_j = dir > 0 ? std::min(j_of(g - 1), next.jb)
                                     : std::max(j_of(g - 1), next.ja);
          const int end = k_of(last_j) + 1;
          relax_to(end - 1, z, i2, j_of(end - 1), cost);
          k = end;
        }
        continue;
      }
      const int j = ja;
      const At v{z, i, j};
      const At w{z, i2, j};
      const At edge = di > 0 ? v : w;
      size_t li = 0;
      inWindow(z, i, j, li);
      const bool source = labels_[li].source_stamp == epoch_;
      if (!region.contains(layer, z, i2, j) || !(source || enter(v, kJog))
          || !enter(w, kJog) || !enter(edge, kJogEdge)) {
        continue;
      }
      const Cost cost
          = layer.jog_factor * std::abs(layer.tracks[i2] - layer.tracks[i])
            + penalty(w, kJog) + penalty(edge, kJogEdge)
            + (source ? 0 : penalty(v, kJog));
      relax_to(0, z, i2, j, cost);
    }
  }

  // Preferred-direction edges leave an interval only at its endpoints.
  for (const int dj : {-1, 1}) {
    const int j = dj < 0 ? ja : jb;
    const int k = k_of(j);
    const int j2 = j + dj;
    if (k < f || k >= g || j2 < 0 || j2 >= layer.numPositions()
        || !region.contains(layer, z, i, j2)) {
      continue;
    }
    const At v{z, i, j};
    const At w{z, i, j2};
    const At edge = dj > 0 ? v : w;
    size_t li = 0;
    inWindow(z, i, j, li);
    if ((labels_[li].source_stamp != epoch_ && !enter(v, kWire))
        || !enter(w, kWire) || !enter(edge, kWireEdge)) {
      continue;
    }
    relax_to(k,
             z,
             i,
             j2,
             std::abs(layer.positions[j2] - layer.positions[j])
                 + penalty(w, kWire) + penalty(edge, kWireEdge));
  }
  pushEvent(id);
}

SearchResult IntervalSearch::run(const SearchRequest& request)
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
  epoch_ = labels_.reset(total, request.label_limit);
  for (size_t t = 0; t < request.targets.size(); t++) {
    int z, i, j;
    graph_.decode(request.targets[t].vertex, z, i, j);
    size_t li = 0;
    if (!inWindow(z, i, j, li)) {
      continue;
    }
    if (labels_[li].target_stamp == epoch_
        && request.targets[labels_[li].target_index].cost
               <= request.targets[t].cost) {
      continue;
    }
    labels_[li].target_stamp = epoch_;
    labels_[li].target_index = static_cast<int32_t>(t);
  }
  for (const auto& [v, c] : request.sources) {
    int z, i, j;
    graph_.decode(v, z, i, j);
    size_t li = 0;
    if (inWindow(z, i, j, li)) {
      labels_[li].source_stamp = epoch_;
    }
  }
  // Algorithm 4, lines 1-4 (with source costs).
  for (const auto& [v, c] : request.sources) {
    int z, i, j;
    graph_.decode(v, z, i, j);
    const int id = intervalOf(z, i, j);
    if (id >= 0) {
      const Interval& iv = intervals_[id];
      addLabel(id, iv.kOf(j), c + iv.futureCost(graph_.layers[z], j), -1, -1);
    }
  }
  while (!heap_.empty()) {
    if (labels_.overflow()) {
      result.aborted = true;
      break;
    }
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
    if (iv.frontier >= iv.n() || iv.frontier != e.aux
        || value(iv, iv.frontier) != e.key
        || iv.futureCost(graph_.layers[iv.z],
                         iv.jOf(frontierEnd(iv, e.key) - 1))
               != e.max_future) {
      continue;
    }
    process(e.interval, e.key);
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
  // Backtrack through the labels.
  int z, i, j;
  graph_.decode(result.target_vertex, z, i, j);
  size_t li = 0;
  inWindow(z, i, j, li);
  int id = labels_[li].interval_id;
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
