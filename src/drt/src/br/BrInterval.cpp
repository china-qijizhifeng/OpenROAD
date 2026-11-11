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

uint16_t IntervalSearch::computeInfo(const int z,
                                     const int i,
                                     const int j,
                                     const bool vias) const
{
  const TrackLayer& layer = graph_.layers[z];
  const int t = req_->wire_type;
  auto usable = [&](const RipupLevel l) {
    return l == kLegal || (req_->ripup && l <= req_->max_ripup_level);
  };
  uint16_t bits = 0;
  if (t >= grid_.numWireTypes()) {
    // Every level comes from the checker.
    const At v{z, i, j};
    bits |= usableZero(v, kWire) ? kInfoWireLegal : 0;
    const bool jog_usable = layer.jogs && enter(v, kJog);
    bits |= layer.jogs && usableZero(v, kJog) ? kInfoJogLegal : 0;
    if (j + 1 < layer.numPositions() && usableZero(v, kWireEdge)) {
      bits |= kInfoEdgeLegal;
    }
    for (const int di : {-1, 1}) {
      const int i2 = i + di;
      size_t li = 0;
      if (!jog_usable || i2 < 0 || i2 >= layer.numTracks()
          || !inWindow(z, i2, j, li)
          || !req_->region->contains(layer, z, i2, j)) {
        continue;
      }
      const At w{z, i2, j};
      const At edge = di > 0 ? v : w;
      const RipupLevel worst = std::max(level(w, kJog), level(edge, kJogEdge));
      if (usable(worst)) {
        bits |= di > 0 ? kTransJogAbove : kTransJogBelow;
        if (!(bits & kInfoJogLegal) || worst != kLegal) {
          bits |= kTransPenalty;  // cannot be compressed into an interval
        }
      }
    }
    return vias ? bits | computeVias(z, i, j, 0) : bits;
  }
  // Fast wire type: decode the words (detail::levelOf, one word per track).
  const int s_wire = FastGrid::wiringShift(t, 0);
  const int s_jog = FastGrid::wiringShift(t, 1);
  const uint64_t own = reader_.wiringWord(z, i, j);
  const RipupLevel wire = FastGrid::field(own, s_wire);
  const RipupLevel jog = FastGrid::field(own, s_jog);
  int64_t queries = 2;
  bits |= wire == kLegal ? kInfoWireLegal : 0;
  bits |= layer.jogs && jog == kLegal ? kInfoJogLegal : 0;
  if (j + 1 < layer.numPositions()) {
    queries++;
    RipupLevel edge;
    if ((own >> FastGrid::kWireEdgeBit) & 1) {
      edge = level(At{z, i, j}, kWireEdge);
      queries--;
    } else {
      edge = std::max(
          wire, FastGrid::field(reader_.wiringWord(z, i, j + 1), s_wire));
    }
    bits |= edge == kLegal ? kInfoEdgeLegal : 0;
  }
  if (layer.jogs && usable(jog)) {
    for (const int di : {-1, 1}) {
      const int i2 = i + di;
      size_t li = 0;
      if (i2 < 0 || i2 >= layer.numTracks() || !inWindow(z, i2, j, li)
          || !req_->region->contains(layer, z, i2, j)) {
        continue;
      }
      const uint64_t other = reader_.wiringWord(z, i2, j);
      const RipupLevel other_jog = FastGrid::field(other, s_jog);
      const uint64_t lower = di > 0 ? own : other;
      RipupLevel jog_edge;
      queries += 2;
      if ((lower >> FastGrid::kJogEdgeBit) & 1) {
        jog_edge = level(At{z, std::min(i, i2), j}, kJogEdge);
        queries--;
      } else {
        jog_edge = std::max(jog, other_jog);
      }
      const RipupLevel worst = std::max(other_jog, jog_edge);
      if (usable(worst)) {
        bits |= di > 0 ? kTransJogAbove : kTransJogBelow;
        if (jog != kLegal || worst != kLegal) {
          bits |= kTransPenalty;  // cannot be compressed into an interval
        }
      }
    }
  }
  counters_.fast_grid += queries;
  return vias ? bits | computeVias(z, i, j, own) : bits;
}

uint16_t IntervalSearch::computeVias(const int z,
                                     const int i,
                                     const int j,
                                     const uint64_t own) const
{
  const TrackLayer& layer = graph_.layers[z];
  const int t = req_->wire_type;
  const bool fast = t < grid_.numWireTypes();
  uint16_t bits = kInfoVias;
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
    RipupLevel worst;
    if (fast) {
      // Lower pad, cut above it and upper pad (detail::levelOf).
      const int zl = std::min(z, z2);
      const int il = dz > 0 ? i : i2;
      const int jl = dz > 0 ? j : j2;
      const uint64_t lower_word = dz > 0 ? own : reader_.wiringWord(zl, il, jl);
      const uint64_t upper_word = dz > 0 ? reader_.wiringWord(z2, i2, j2) : own;
      const int c = grid_.viaIndex(zl, jl);
      RipupLevel cut = kNotRemovable;
      if (c >= 0) {
        const uint64_t w = reader_.viaWord(zl, il, c);
        cut = std::max(FastGrid::field(w, FastGrid::viaShift(t, 0)),
                       FastGrid::field(w, FastGrid::viaShift(t, 1)));
      }
      worst = std::max(
          {FastGrid::field(lower_word, FastGrid::wiringShift(t, 2)),
           cut,
           FastGrid::field(upper_word, FastGrid::wiringShift(t, 3))});
      counters_.fast_grid += 3;
    } else {
      const At v{z, i, j};
      const At w{z2, i2, j2};
      const At lower = dz > 0 ? v : w;
      const At upper = dz > 0 ? w : v;
      worst = std::max({level(lower, kViaUpPad),
                        level(lower, kCutUp),
                        level(upper, kViaDownPad)});
    }
    if (worst == kLegal || (req_->ripup && worst <= req_->max_ripup_level)) {
      bits |= dz > 0 ? kTransViaAbove : kTransViaBelow;
      if (worst != kLegal) {
        bits |= kTransPenalty;
      }
    }
  }
  return bits;
}

uint16_t IntervalSearch::info(const int z,
                              const int i,
                              const int j,
                              VertexLabel& label) const
{
  if (label.info_stamp != epoch_) {
    label.info = computeInfo(z, i, j, vias_at_build_);
    label.info_stamp = epoch_;
  }
  return label.info;
}

uint16_t IntervalSearch::viaInfo(const int z,
                                 const int i,
                                 const int j,
                                 VertexLabel& label) const
{
  const uint16_t bits = info(z, i, j, label);
  if (bits & kInfoVias) {
    return bits;
  }
  label.info = bits
               | computeVias(z,
                             i,
                             j,
                             req_->wire_type < grid_.numWireTypes()
                                 ? reader_.wiringWord(z, i, j)
                                 : 0);
  return label.info;
}

int IntervalSearch::intervalOf(const int z, const int i, const int j)
{
  size_t li = 0;
  if (!inWindow(z, i, j, li)) {
    return -1;
  }
  if (label(li).stamp == epoch_) {
    return label(li).interval_id;
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
  const uint16_t property_mask
      = kTransJogBelow | kTransJogAbove | kTransPenalty
        | (req_->interval_via_property ? kTransViaBelow | kTransViaAbove : 0);
  auto plain = [&](int jj, uint16_t& sig, uint16_t& bits) {
    size_t li = 0;
    if (!inWindow(z, i, jj, li)) {
      return false;
    }
    VertexLabel& lab = label(li);
    if (lab.stamp == epoch_ || lab.target_stamp == epoch_
        || !region.contains(layer, z, i, jj)) {
      return false;
    }
    bits = info(z, i, jj, lab);
    sig = bits & property_mask;
    return (bits & kInfoWireLegal) != 0
           && !(lab.source_stamp == epoch_ && layer.jogs
                && !(bits & kInfoJogLegal))
           && (bits & kTransPenalty) == 0;
  };
  int ja = j;
  int jb = j;
  int cls = 0;
  // pi_H along track i: per target box, the cost across the track
  // (cross distance, via bound, terminal cost) and its range along it.
  track_future_.clear();
  for (const TargetBox& b : req_->target_boxes) {
    const Coord c = layer.tracks[i];
    const Coord clo = layer.horizontal ? b.ylo : b.xlo;
    const Coord chi = layer.horizontal ? b.yhi : b.xhi;
    track_future_.push_back(
        {std::max({Cost{0}, Cost{clo} - c, Cost{c} - chi}) + viaBound(z, b)
             + (req_->exact_future_cost ? b.cost : 0),
         layer.horizontal ? b.xlo : b.ylo,
         layer.horizontal ? b.xhi : b.yhi});
  }
  auto future = [&](int jj) {
    if (track_future_.empty()) {
      return Cost{0};
    }
    const Cost a = layer.positions[jj];
    Cost best = kInfCost;
    for (const auto& [cross, lo, hi] : track_future_) {
      best = std::min(best, cross + std::max({Cost{0}, lo - a, a - hi}));
    }
    return best;
  };
  const Cost pi0 = future(j);
  Cost pi_a = pi0;
  Cost pi_b = pi0;
  uint16_t sig = 0;
  uint16_t bits_a = 0;
  if (!plain(j, sig, bits_a)) {
    ends_[kEndStartVertex] += 2;
  } else {
    IntervalEnd end = kEndTrack;
    uint16_t bits_b = bits_a;
    auto joins = [&](int from, int to, Cost& pi_from, uint16_t& bits_from) {
      uint16_t sig_to = 0;
      uint16_t bits_to = 0;
      if (!plain(to, sig_to, bits_to)) {
        end = kEndNeighbour;
        return false;
      }
      // The preferred-direction edge is kept at its lower end.
      if (!((to > from ? bits_from : bits_to) & kInfoEdgeLegal)) {
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
      const Cost pi_to = future(to);
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
      bits_from = bits_to;
      return true;
    };
    end = kEndTrack;
    while (jb + 1 < layer.numPositions() && joins(jb, jb + 1, pi_b, bits_b)) {
      jb++;
    }
    ends_[end]++;
    end = kEndTrack;
    while (ja > 0 && joins(ja, ja - 1, pi_a, bits_a)) {
      ja--;
    }
    ends_[end]++;
  }
  interval_vertices_ += jb - ja + 1;
  // Interval objects are reused across searches with their label storage.
  const int id = num_used_;
  if (num_used_ == static_cast<int>(intervals_.size())) {
    intervals_.emplace_back();
  }
  num_used_++;
  Interval& iv = intervals_[id];
  iv.z = z;
  iv.i = i;
  iv.ja = ja;
  iv.jb = jb;
  // pi increases with j (cls +1): the target end is ja.
  iv.dir = cls >= 0 ? 1 : -1;
  // The checked slope determines pi everywhere in the interval.
  iv.future_at_ja = pi_a;
  iv.frontier = 0;
  iv.labels.clear();
  for (int jj = ja; jj <= jb; jj++) {
    size_t li = 0;
    inWindow(z, i, jj, li);
    label(li).stamp = epoch_;
    label(li).interval_id = id;
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
  if (iv.labels.size() <= 1) {
    return iv.labels.empty() ? nullptr : &iv.labels[0];
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
  if (n <= 1) {
    return;  // argLabel needs no envelope
  }
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
  if (iv.labels.size() == 1) {
    // Flat from the label towards the target end, increasing beyond it.
    if (value(iv, iv.frontier) != delta) {
      return iv.frontier;
    }
    return std::max(iv.frontier, iv.labels[0].k) + 1;
  }
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
  res_->heap_pushes++;
  if (direct_) {
    direct_event_ = Event{key, tieKey(pi), interval, iv.frontier};
    return;
  }
  heap_.push_back(Event{key, tieKey(pi), interval, iv.frontier});
  std::push_heap(heap_.begin(), heap_.end(), std::greater<>());
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
  if (iv.ja == iv.jb) {
    // Targets are single-vertex intervals, so terminal costs are exact. The
    // sink is queued with the label: once its key is the minimum, no other
    // path is cheaper (Algorithm 4, line 25, checked as early as possible).
    size_t li = 0;
    inWindow(iv.z, iv.i, iv.ja, li);
    if (label(li).target_stamp == epoch_) {
      const auto& t = req_->targets[label(li).target_index];
      // delta includes pi(t) <= t.cost (the box cost is the cheapest
      // terminal cost), so the sink key stays >= delta.
      const Cost pi_t = iv.futureCost(layer, iv.ja);
      heap_.push_back(Event{
          delta - pi_t + t.cost, kSinkTie, -1, label(li).target_index});
      std::push_heap(heap_.begin(), heap_.end(), std::greater<>());
      res_->heap_pushes++;
    }
  }
  pushEvent(interval);
}

void IntervalSearch::relaxJogs(const JogRun& run)
{
  // Algorithm 4, lines 13-17, for the jogs of J = [f, g) of a plain
  // interval: one relaxation per neighbouring interval, at the member of
  // the overlap with maximum pi (its last one).
  const int id = run.interval;
  const int z = intervals_[id].z;
  const int i = intervals_[id].i;
  const int i2 = i + run.di;
  const TrackLayer& layer = graph_.layers[z];
  const Cost cost
      = layer.jog_factor * std::abs(layer.tracks[i2] - layer.tracks[i]);
  const int last = intervals_[id].jOf(run.g - 1);
  for (int k = run.f; k < run.g;) {
    // intervalOf may add intervals: no references across it.
    const int next_id = intervalOf(z, i2, intervals_[id].jOf(k));
    const Interval& from = intervals_[id];
    const Interval& next = intervals_[next_id];
    const int last_j
        = from.dir > 0 ? std::min(last, next.jb) : std::max(last, next.ja);
    const int end = from.kOf(last_j) + 1;
    const Cost d2 = run.delta + cost - from.futureCost(layer, last_j)
                    + next.futureCost(layer, last_j);
    addLabel(next_id, next.kOf(last_j), d2, id, end - 1);
    k = end;
  }
}

void IntervalSearch::process(const int id, const Cost delta, int g)
{
  const int f = intervals_[id].frontier;
  if (g < 0) {
    g = frontierEnd(intervals_[id], delta);
  }
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
    return label(li).info;
  };
  auto via_bits_of = [&](int k) {
    size_t li = 0;
    inWindow(z, i, j_of(k), li);
    return viaInfo(z, i, j_of(k), label(li));
  };
  const uint16_t signature = single ? 0 : bits_of(f);
  // Relaxation of an edge from member from_k to (z2, i2, j2), in two steps:
  // probe() looks the head up once and says whether a label of cost at
  // least `cost` could improve it (before any legality is read); commit()
  // adds the label with the full cost. A vertex outside every interval
  // gets a deferred label: its interval is built when the label is the
  // cheapest one, so vertices the search only reaches at its frontier
  // never form intervals.
  struct Probe
  {
    VertexLabel* head;
    int z, i, j;
    Cost pi;  // future cost at the head
    Cost d2;  // label value for the probed cost
  };
  auto probe = [&](int from_k, int z2, int i2, int j2, Cost cost, Probe& p) {
    size_t li = 0;
    if (i2 < 0 || j2 < 0 || !inWindow(z2, i2, j2, li)) {
      return false;
    }
    VertexLabel& w = label(li);
    const Cost from = delta + cost - intervals_[id].futureCost(layer, j_of(from_k));
    p = Probe{&w, z2, i2, j2, 0, 0};
    if (w.stamp == epoch_) {
      const Interval& next = intervals_[w.interval_id];
      const int k2 = next.kOf(j2);
      p.pi = next.futureCost(graph_.layers[z2], j2);
      p.d2 = from + p.pi;
      return k2 >= next.frontier && value(next, k2) > p.d2;
    }
    p.pi = futureCost(z2, i2, j2);
    p.d2 = from + p.pi;
    return w.pending_stamp != epoch_ || w.pending_delta > p.d2;
  };
  auto commit = [&](int from_k, const Probe& p, Cost extra) {
    VertexLabel& w = *p.head;
    const Cost d2 = p.d2 + extra;
    if (w.stamp == epoch_) {
      const Interval& next = intervals_[w.interval_id];
      addLabel(w.interval_id, next.kOf(p.j), d2, id, from_k);
      return;
    }
    if (w.pending_stamp == epoch_ && w.pending_delta <= d2) {
      return;
    }
    w.pending_stamp = epoch_;
    w.pending_delta = d2;
    w.pending_pred_interval = id;
    w.pending_pred_k = from_k;
    pending_.push_back({p.z, p.i, p.j});
    heap_.push_back(Event{d2,
                          tieKey(p.pi),
                          kPendingEvent,
                          static_cast<int32_t>(pending_.size() - 1)});
    std::push_heap(heap_.begin(), heap_.end(), std::greater<>());
    res_->heap_pushes++;
  };
  auto relax_to = [&](int from_k, int z2, int i2, int j2, Cost cost) {
    Probe p;
    if (probe(from_k, z2, i2, j2, cost, p)) {
      commit(from_k, p, 0);
    }
  };

  // On perpendicular layers each crossing reaches a different track and
  // hence a different neighboring interval. Enumerate those adjacencies;
  // skip the entire run when the interval has no vias in that direction.
  for (const int dz : {-1, 1}) {
    const int z2 = z + dz;
    const uint16_t via_bit = dz > 0 ? kTransViaAbove : kTransViaBelow;
    if (z2 < 0 || z2 >= graph_.numLayers()
        || (!single && req_->interval_via_property && !(signature & via_bit))) {
      continue;
    }
    for (int k = f; k < g; k++) {
      const int j = j_of(k);
      const int i2 = dz > 0 ? layer.up_track[j] : layer.down_track[j];
      const int j2 = dz > 0 ? layer.up_pos[i] : layer.down_pos[i];
      Probe p;
      if (!probe(k, z2, i2, j2, graph_.via_cost[std::min(z, z2)], p)) {
        continue;
      }
      if (!single) {
        if (via_bits_of(k) & via_bit) {
          commit(k, p, 0);
        }
        continue;
      }
      if (!region.contains(graph_.layers[z2], z2, i2, j2)) {
        continue;
      }
      const At v{z, i, j};
      const At w{z2, i2, j2};
      const At lower = dz > 0 ? v : w;
      const At upper = dz > 0 ? w : v;
      if (enter(lower, kViaUpPad) && enter(lower, kCutUp)
          && enter(upper, kViaDownPad)) {
        commit(k,
               p,
               penalty(lower, kViaUpPad) + penalty(lower, kCutUp)
                   + penalty(upper, kViaDownPad));
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
      const uint16_t jog_bit = di > 0 ? kTransJogAbove : kTransJogBelow;
      if (i2 < 0 || i2 >= layer.numTracks()
          || (!single && !(signature & jog_bit))) {
        continue;
      }
      if (!single) {
        // Every member has the jog (property (ii)). The jogs are relaxed
        // (relaxJogs) when their cheapest possible label is the cheapest
        // event: pi is 1-Lipschitz, so every jog label of the run is at
        // least delta + (beta_z - 1) * |track distance|.
        const Cost distance = std::abs(layer.tracks[i2] - layer.tracks[i]);
        jog_runs_.push_back({id, f, g, di, delta});
        heap_.push_back(Event{
            delta + (layer.jog_factor - 1) * distance,
            tieKey(intervals_[id].futureCost(layer, j_of(g - 1)) - distance),
            kJogEvent,
            static_cast<int32_t>(jog_runs_.size() - 1)});
        std::push_heap(heap_.begin(), heap_.end(), std::greater<>());
        res_->heap_pushes++;
        continue;
      }
      const int j = ja;
      const At v{z, i, j};
      const At w{z, i2, j};
      const At edge = di > 0 ? v : w;
      size_t li = 0;
      inWindow(z, i, j, li);
      const bool source = label(li).source_stamp == epoch_;
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
    Probe p;
    if (!probe(k,
               z,
               i,
               j2,
               std::abs(layer.positions[j2] - layer.positions[j]),
               p)) {
      continue;
    }
    const At v{z, i, j};
    const At w{z, i, j2};
    const At edge = dj > 0 ? v : w;
    size_t li = 0;
    inWindow(z, i, j, li);
    if ((label(li).source_stamp != epoch_ && !enter(v, kWire))
        || !enter(w, kWire) || !enter(edge, kWireEdge)) {
      continue;
    }
    commit(k, p, penalty(w, kWire) + penalty(edge, kWireEdge));
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
  pending_.clear();
  jog_runs_.clear();
  nearest_first_ = request.interval_nearest_first;
  num_used_ = 0;
  if (intervals_.size() > kKeptIntervals) {
    intervals_.resize(kKeptIntervals);
  }
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
  label_page_ = ~size_t{0};
  vias_at_build_ = request.interval_via_property || request.ripup;
  for (size_t t = 0; t < request.targets.size(); t++) {
    int z, i, j;
    graph_.decode(request.targets[t].vertex, z, i, j);
    size_t li = 0;
    if (!inWindow(z, i, j, li)) {
      continue;
    }
    if (label(li).target_stamp == epoch_
        && request.targets[label(li).target_index].cost
               <= request.targets[t].cost) {
      continue;
    }
    label(li).target_stamp = epoch_;
    label(li).target_index = static_cast<int32_t>(t);
  }
  for (const auto& [v, c] : request.sources) {
    int z, i, j;
    graph_.decode(v, z, i, j);
    size_t li = 0;
    if (inWindow(z, i, j, li)) {
      label(li).source_stamp = epoch_;
    }
  }
  // Algorithm 4, lines 1-4 (with source costs): every source gets a label
  // (delta = cost + pi), deferred like the labels of relaxations, so the
  // intervals of sources far from the targets (a routed tree) are built
  // only if the search reaches their level.
  for (const auto& [v, c] : request.sources) {
    int z, i, j;
    graph_.decode(v, z, i, j);
    size_t li = 0;
    if (!inWindow(z, i, j, li)) {
      continue;
    }
    VertexLabel& w = label(li);
    const Cost pi = futureCost(z, i, j);
    const Cost d = c + pi;
    if (w.pending_stamp == epoch_ && w.pending_delta <= d) {
      continue;
    }
    w.pending_stamp = epoch_;
    w.pending_delta = d;
    w.pending_pred_interval = -1;
    w.pending_pred_k = -1;
    pending_.push_back({z, i, j});
    heap_.push_back(Event{
        d, tieKey(pi), kPendingEvent, static_cast<int32_t>(pending_.size() - 1)});
    std::push_heap(heap_.begin(), heap_.end(), std::greater<>());
    result.heap_pushes++;
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
    if (e.interval == kJogEvent) {
      const JogRun run = jog_runs_[e.aux];
      relaxJogs(run);
      continue;
    }
    if (e.interval == kPendingEvent) {
      const Pending p = pending_[e.aux];
      size_t li = 0;
      inWindow(p.z, p.i, p.j, li);
      VertexLabel& w = label(li);
      if (w.pending_stamp != epoch_ || w.pending_delta != e.key) {
        continue;  // superseded by a cheaper deferred label
      }
      w.pending_stamp = 0;
      const int pred_interval = w.pending_pred_interval;
      const int pred_k = w.pending_pred_k;
      const int id = intervalOf(p.z, p.i, p.j);
      const Interval& iv = intervals_[id];
      // The new interval event is usually the next one: process it without
      // a heap round trip when no queued event precedes it.
      direct_ = true;
      addLabel(id, iv.kOf(p.j), e.key, pred_interval, pred_k);
      direct_ = false;
      if (direct_event_.interval >= 0) {
        const Event next = direct_event_;
        direct_event_.interval = -1;
        if (heap_.empty() || !(next > heap_.front())) {
          result.heap_pops++;
          process(next.interval, next.key, -1);
        } else {
          heap_.push_back(next);
          std::push_heap(heap_.begin(), heap_.end(), std::greater<>());
        }
      }
      continue;
    }
    if (e.interval < 0) {
      result.found = true;
      result.cost = e.key;
      result.target_id = request.targets[e.aux].id;
      result.target_vertex = request.targets[e.aux].vertex;
      break;
    }
    const Interval& iv = intervals_[e.interval];
    if (iv.frontier >= iv.n() || iv.frontier != e.aux
        || value(iv, iv.frontier) != e.key) {
      continue;
    }
    const int end = frontierEnd(iv, e.key);
    if (tieKey(iv.futureCost(graph_.layers[iv.z], iv.jOf(end - 1)))
        != e.max_future) {
      continue;
    }
    process(e.interval, e.key, end);
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
  int id = label(li).interval_id;
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
