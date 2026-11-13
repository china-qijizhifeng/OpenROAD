// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrPinAccess.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

namespace drt::br {

void CircuitSignature::addShape(const int layer,
                                const Box& chip_box,
                                const int label)
{
  shapes_.emplace_back(layer, to_local_.apply(chip_box), label);
}

void CircuitSignature::addTrack(const int layer,
                                const bool horizontal,
                                const GeoCoord coord)
{
  // The line y = c (or x = c) through (0, c) (or (c, 0)).
  const GridPoint p
      = to_local_.apply(horizontal ? GridPoint{0, coord} : GridPoint{coord, 0});
  const bool local_horizontal = horizontal != to_local_.swapsAxes();
  tracks_.emplace_back(
      layer, local_horizontal ? 1 : 0, local_horizontal ? p.y : p.x);
}

void CircuitSignature::addPin(const int pin, const int group)
{
  pins_.emplace_back(pin, group);
}

void CircuitSignature::addVia(const int via,
                              const std::vector<std::pair<int, Box>>& shapes)
{
  const Isometry linear = to_local_.linear();
  for (const auto& [layer, box] : shapes) {
    vias_.emplace_back(via, layer, linear.apply(box));
  }
}

std::vector<int64_t> CircuitSignature::finish() const
{
  auto shapes = shapes_;
  auto tracks = tracks_;
  auto pins = pins_;
  auto vias = vias_;
  std::sort(shapes.begin(), shapes.end());
  shapes.erase(std::unique(shapes.begin(), shapes.end()), shapes.end());
  std::sort(tracks.begin(), tracks.end());
  tracks.erase(std::unique(tracks.begin(), tracks.end()), tracks.end());
  std::sort(pins.begin(), pins.end());
  std::sort(vias.begin(), vias.end());
  std::vector<int64_t> out;
  out.push_back(static_cast<int64_t>(values_.size()));
  out.insert(out.end(), values_.begin(), values_.end());
  out.push_back(static_cast<int64_t>(pins.size()));
  for (const auto& [pin, group] : pins) {
    out.push_back(pin);
    out.push_back(group);
  }
  out.push_back(static_cast<int64_t>(vias.size()));
  for (const auto& [via, layer, b] : vias) {
    out.insert(out.end(), {via, layer, b.xlo, b.ylo, b.xhi, b.yhi});
  }
  out.push_back(static_cast<int64_t>(tracks.size()));
  for (const auto& [layer, horizontal, c] : tracks) {
    out.insert(out.end(), {layer, horizontal, c});
  }
  out.push_back(static_cast<int64_t>(shapes.size()));
  for (const auto& [layer, b, label] : shapes) {
    out.insert(out.end(), {layer, b.xlo, b.ylo, b.xhi, b.yhi, label});
  }
  return out;
}

std::pair<int, bool> CircuitClassTable::classify(
    const std::vector<int64_t>& signature)
{
  auto [it, inserted]
      = classes_.emplace(signature, static_cast<int>(classes_.size()));
  return {it->second, inserted};
}

namespace {

int64_t spreadPenalty(const SelectionWeights& w,
                      const AccessCandidate& a,
                      const AccessCandidate& b)
{
  if (w.spread == 0 || a.end_layer != b.end_layer) {
    return 0;
  }
  const int64_t d = std::abs(int64_t{a.end.x} - b.end.x)
                    + std::abs(int64_t{a.end.y} - b.end.y);
  return w.spread * std::max<int64_t>(0, w.spread_distance - d);
}

int64_t individual(const SelectionWeights& w, const AccessCandidate& c)
{
  return c.length + w.direction * c.missing_directions;
}

// Pairwise conflicts of all candidates, flattened.
struct ConflictTable
{
  std::vector<int> offset;  // first flat index of each pin
  int total = 0;
  std::vector<char> bits;
  bool get(int pa, int ca, int pb, int cb) const
  {
    return bits[static_cast<size_t>(offset[pa] + ca) * total + offset[pb] + cb]
           != 0;
  }
};

ConflictTable buildConflicts(const SelectionInput& in)
{
  ConflictTable t;
  const int n = static_cast<int>(in.pins.size());
  t.offset.resize(n);
  for (int p = 0; p < n; p++) {
    t.offset[p] = t.total;
    t.total += static_cast<int>(in.pins[p].size());
  }
  t.bits.assign(static_cast<size_t>(t.total) * t.total, 0);
  for (int pa = 0; pa < n; pa++) {
    for (int pb = pa + 1; pb < n; pb++) {
      if (in.group[pa] == in.group[pb] || !in.conflict) {
        continue;
      }
      for (int ca = 0; ca < static_cast<int>(in.pins[pa].size()); ca++) {
        for (int cb = 0; cb < static_cast<int>(in.pins[pb].size()); cb++) {
          if (in.conflict(pa, ca, pb, cb)) {
            const size_t a = t.offset[pa] + ca;
            const size_t b = t.offset[pb] + cb;
            t.bits[a * t.total + b] = 1;
            t.bits[b * t.total + a] = 1;
          }
        }
      }
    }
  }
  return t;
}

class DestructiveBounding
{
 public:
  DestructiveBounding(const SelectionInput& in, const ConflictTable& conflicts)
      : in_(in), conflicts_(conflicts), n_(static_cast<int>(in.pins.size()))
  {
    // One frame per search depth (a level assigns one pin, so depths
    // 0 .. n), allocated once and reused by every node of that depth.
    frames_.resize(n_ + 1);
    for (Frame& f : frames_) {
      f.reduced.resize(n_);
      f.inc.resize(n_);
    }
    min_individual_.resize(n_);
    for (int p = 0; p < n_; p++) {
      int64_t m = in_.weights.unassigned;
      for (const AccessCandidate& c : in_.pins[p]) {
        m = std::min(m, individual(in_.weights, c));
      }
      min_individual_[p] = m;
    }
  }

  SelectionResult run()
  {
    SelectionResult result;
    int64_t bound = 0;
    for (int p = 0; p < n_; p++) {
      bound += min_individual_[p];
    }
    std::vector<std::vector<int>> domains(n_);
    for (int p = 0; p < n_; p++) {
      for (int c = -1; c < static_cast<int>(in_.pins[p].size()); c++) {
        domains[p].push_back(c);
      }
    }
    while (true) {
      result.iterations++;
      bound_ = bound;
      next_ = std::numeric_limits<int64_t>::max();
      choice_.assign(n_, kOpen);
      track_count_.clear();
      if (search(domains, 0, 0)) {
        result.choice = choice_;
        result.cost = selectionCost(in_, choice_);
        result.optimal = true;
        break;
      }
      if (nodes_ >= in_.node_limit
          || next_ == std::numeric_limits<int64_t>::max()) {
        const int iterations = result.iterations;
        result = selectGreedy(in_);
        result.iterations = iterations;
        result.optimal = false;
        break;
      }
      bound = next_;
    }
    result.nodes = nodes_;
    return result;
  }

 private:
  static constexpr int kOpen = -2;

  // Exact cost increase of adding candidate c (or -1) for pin p.
  int64_t increment(int p, int c) const
  {
    if (c < 0) {
      return in_.weights.unassigned;
    }
    const AccessCandidate& cand = in_.pins[p][c];
    int64_t inc = individual(in_.weights, cand);
    if (in_.weights.track != 0) {
      for (const int64_t t : cand.blocked_tracks) {
        auto it = track_count_.find(t);
        if (it == track_count_.end() || it->second == 0) {
          inc += in_.weights.track;
        }
      }
    }
    for (int q = 0; q < n_; q++) {
      if (choice_[q] >= 0) {
        inc += spreadPenalty(in_.weights, cand, in_.pins[q][choice_[q]]);
      }
    }
    return inc;
  }

  bool compatible(int p, int c) const
  {
    if (c < 0) {
      return true;
    }
    for (int q = 0; q < n_; q++) {
      if (choice_[q] >= 0 && in_.group[q] != in_.group[p]
          && conflicts_.get(p, c, q, choice_[q])) {
        return false;
      }
    }
    return true;
  }

  void assign(int p, int c)
  {
    choice_[p] = c;
    if (c >= 0) {
      for (const int64_t t : in_.pins[p][c].blocked_tracks) {
        track_count_[t]++;
      }
    }
  }

  void unassign(int p)
  {
    const int c = choice_[p];
    if (c >= 0) {
      for (const int64_t t : in_.pins[p][c].blocked_tracks) {
        track_count_[t]--;
      }
    }
    choice_[p] = kOpen;
  }

  // `g`: exact cost of the assigned pins.
  bool search(const std::vector<std::vector<int>>& domains,
              const int64_t g,
              const int depth)
  {
    nodes_++;
    if (nodes_ > in_.node_limit) {
      return false;
    }
    int64_t rest = 0;
    int open = 0;
    for (int p = 0; p < n_; p++) {
      if (choice_[p] == kOpen) {
        rest += min_individual_[p];
        open++;
      }
    }
    if (open == 0) {
      return true;
    }
    // Reduce every open domain under the bound; an empty domain refutes
    // the node.
    Frame& frame = frames_[depth];
    std::vector<std::vector<int>>& reduced = frame.reduced;
    std::vector<std::vector<int64_t>>& inc = frame.inc;
    for (int p = 0; p < n_; p++) {
      reduced[p].clear();
      inc[p].clear();
    }
    int best_pin = -1;
    for (int p = 0; p < n_; p++) {
      if (choice_[p] != kOpen) {
        continue;
      }
      for (const int c : domains[p]) {
        if (!compatible(p, c)) {
          continue;
        }
        const int64_t i = increment(p, c);
        const int64_t lb = g + i + rest - min_individual_[p];
        if (lb > bound_) {
          next_ = std::min(next_, lb);
          continue;
        }
        reduced[p].push_back(c);
        inc[p].push_back(i);
      }
      if (reduced[p].empty()) {
        return false;
      }
      if (best_pin < 0 || reduced[p].size() < reduced[best_pin].size()) {
        best_pin = p;
      }
    }
    const int p = best_pin;
    std::vector<int>& order = frame.order;
    order.resize(reduced[p].size());
    for (size_t k = 0; k < order.size(); k++) {
      order[k] = static_cast<int>(k);
    }
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
      return inc[p][a] < inc[p][b];
    });
    for (const int k : order) {
      assign(p, reduced[p][k]);
      if (search(reduced, g + inc[p][k], depth + 1)) {
        return true;
      }
      unassign(p);
      if (nodes_ > in_.node_limit) {
        return false;
      }
    }
    return false;
  }

  // Scratch of one search depth: reduced domains, their increments and the
  // candidate order of the branching pin.
  struct Frame
  {
    std::vector<std::vector<int>> reduced;
    std::vector<std::vector<int64_t>> inc;
    std::vector<int> order;
  };

  const SelectionInput& in_;
  const ConflictTable& conflicts_;
  const int n_;
  std::vector<Frame> frames_;
  std::vector<int64_t> min_individual_;
  std::vector<int> choice_;
  std::map<int64_t, int> track_count_;
  int64_t bound_ = 0;
  int64_t next_ = 0;
  int64_t nodes_ = 0;
};

}  // namespace

int64_t selectionCost(const SelectionInput& in, const std::vector<int>& choice)
{
  const int n = static_cast<int>(in.pins.size());
  int64_t cost = 0;
  std::vector<int64_t> tracks;
  for (int p = 0; p < n; p++) {
    const int c = choice[p];
    if (c < 0) {
      cost += in.weights.unassigned;
      continue;
    }
    const AccessCandidate& a = in.pins[p][c];
    cost += individual(in.weights, a);
    tracks.insert(
        tracks.end(), a.blocked_tracks.begin(), a.blocked_tracks.end());
    for (int q = p + 1; q < n; q++) {
      if (choice[q] < 0) {
        continue;
      }
      if (in.group[p] != in.group[q] && in.conflict
          && in.conflict(p, c, q, choice[q])) {
        return -1;
      }
      cost += spreadPenalty(in.weights, a, in.pins[q][choice[q]]);
    }
  }
  std::sort(tracks.begin(), tracks.end());
  tracks.erase(std::unique(tracks.begin(), tracks.end()), tracks.end());
  cost += in.weights.track * static_cast<int64_t>(tracks.size());
  return cost;
}

SelectionResult selectConflictFree(const SelectionInput& in)
{
  const ConflictTable conflicts = buildConflicts(in);
  DestructiveBounding bounding(in, conflicts);
  return bounding.run();
}

std::vector<CatalogueVariant> selectionCompletionVariants(
    const bool algorithm3_fallback,
    const int mode,
    const int widenings)
{
  std::vector<CatalogueVariant> variants;
  auto completion = [&variants](const int widening) {
    for (const bool algorithm3 : {false, true}) {
      variants.push_back({.algorithm3 = algorithm3,
                          .all_access_points = true,
                          .extra_radius = widening,
                          .paths_factor = 1 << widening});
    }
  };
  if (mode == 2) {
    completion(0);
  }
  if (algorithm3_fallback) {
    variants.push_back({.algorithm3 = true});
  }
  if (mode == 1) {
    completion(0);
  }
  if (mode >= 1) {
    for (int widening = 1; widening <= widenings; widening++) {
      completion(widening);
    }
  }
  return variants;
}

int bestSelectionAttempt(const int pins,
                         const int assigned,
                         const std::vector<int>& attempts)
{
  int best = -1;
  int best_assigned = assigned;
  for (size_t k = 0; k < attempts.size(); k++) {
    if (attempts[k] >= pins) {
      return static_cast<int>(k);
    }
    if (attempts[k] > best_assigned) {
      best = static_cast<int>(k);
      best_assigned = attempts[k];
    }
  }
  return best;
}

SelectionResult selectGreedy(const SelectionInput& in)
{
  SelectionResult r;
  const int n = static_cast<int>(in.pins.size());
  r.choice.assign(n, -1);
  for (int p = 0; p < n; p++) {
    for (int c = 0; c < static_cast<int>(in.pins[p].size()); c++) {
      bool ok = true;
      for (int q = 0; q < p && ok; q++) {
        ok = r.choice[q] < 0 || in.group[q] == in.group[p] || !in.conflict
             || !in.conflict(p, c, q, r.choice[q]);
      }
      if (ok) {
        r.choice[p] = c;
        break;
      }
    }
  }
  r.cost = selectionCost(in, r.choice);
  return r;
}

int64_t pinContactViolations(const std::vector<RouteNode>& contacts,
                             const std::vector<ContactPin>& pins)
{
  int64_t violations = 0;
  for (const RouteNode& c : contacts) {
    for (const ContactPin& pin : pins) {
      const bool touches = std::any_of(
          pin.shapes.begin(), pin.shapes.end(), [&](const auto& s) {
            return s.first == c.layer && s.second.containsClosed(c.p.x, c.p.y);
          });
      if (touches
          && std::find(pin.access_points.begin(), pin.access_points.end(), c)
                 == pin.access_points.end()) {
        violations++;
      }
    }
  }
  return violations;
}

}  // namespace drt::br
