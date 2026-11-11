// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Off-track pin access (Gester et al., TODAES 18(2), 2013, Sec. 4.3):
// circuit classes (instances in geometrically equal surroundings up to
// translation, mirroring and rotation) and the choice of a conflict-free
// set of primary access paths per circuit by destructive bounding.

#pragma once

#include <cstdint>
#include <functional>
#include <map>
#include <tuple>
#include <vector>

#include "br/BrGeometry.h"

namespace drt::br {

// Canonical description of an instance's surroundings in the instance's
// local frame: every shape within the expanded bounding area (clipped to
// it) with its relation to the instance's pins, the track coordinates and
// the via geometry. Two instances whose signatures are equal see the same
// situation after mapping one local frame onto the other.
class CircuitSignature
{
 public:
  // `to_local` maps chip coordinates into the instance's local frame.
  explicit CircuitSignature(const Isometry& to_local) : to_local_(to_local) {}
  // A shape in chip coordinates; `label` 0 = other net / blockage, k > 0 =
  // shape of the net of the instance's pin group k.
  void addShape(int layer, const Box& chip_box, int label);
  // A track line in chip coordinates: y = coord if `horizontal`, else
  // x = coord.
  void addTrack(int layer, bool horizontal, GeoCoord coord);
  // Pin `pin` of the instance needs access; pins of one net share `group`.
  void addPin(int pin, int group);
  // Via model `via` (index into the caller's via list) with shapes given
  // relative to the via origin (layer, box).
  void addVia(int via, const std::vector<std::pair<int, Box>>& shapes);
  // Free-form parameters that must match (window size, radius, ...).
  void addValue(int64_t value) { values_.push_back(value); }
  std::vector<int64_t> finish() const;

 private:
  Isometry to_local_;
  std::vector<std::tuple<int, Box, int>> shapes_;
  std::vector<std::tuple<int, int, GeoCoord>> tracks_;
  std::vector<std::pair<int, int>> pins_;
  std::vector<std::tuple<int, int, Box>> vias_;
  std::vector<int64_t> values_;
};

// Partition of instances into circuit classes by equal signatures.
class CircuitClassTable
{
 public:
  // Class id of the signature (a new id for an unseen one) and whether it
  // is new.
  std::pair<int, bool> classify(const std::vector<int64_t>& signature);
  int numClasses() const { return static_cast<int>(classes_.size()); }

 private:
  std::map<std::vector<int64_t>, int> classes_;
};

// One catalogue entry of a pin, with the features of Sec. 4.3's
// evaluation: "spreading of endpoints, number of blocked tracks, directions
// of feasible on-track continuation, and length".
struct AccessCandidate
{
  int64_t length = 0;  // path cost (wire length and vias)
  int end_layer = 0;
  GridPoint end;
  // Tracks (caller ids) the path blocks for other nets, sorted and unique.
  std::vector<int64_t> blocked_tracks;
  // Feasible on-track continuation directions missing at the endpoint.
  int missing_directions = 0;
};

struct SelectionWeights
{
  int64_t track = 0;      // per blocked track (union over the solution)
  int64_t direction = 0;  // per missing continuation direction
  int64_t spread = 0;     // per DBU two endpoints are closer than below
  GeoCoord spread_distance = 0;
  int64_t unassigned = int64_t{1} << 40;  // pin without primary path
};

struct SelectionInput
{
  std::vector<std::vector<AccessCandidate>> pins;
  // Pins of equal group belong to one net and never conflict.
  std::vector<int> group;
  // Diff-net conflict of candidate ca of pin pa with candidate cb of pb.
  std::function<bool(int pa, int ca, int pb, int cb)> conflict;
  SelectionWeights weights;
  int64_t node_limit = 100000;
};

struct SelectionResult
{
  std::vector<int> choice;  // candidate per pin, -1 = none
  int64_t cost = 0;
  int64_t nodes = 0;
  int iterations = 0;
  // Proven minimum of the cost model; false when the node limit stopped
  // the enumeration and the greedy solution was taken.
  bool optimal = false;
};

// Exact cost of a (partial) choice under the weights:
//   sum of lengths + direction * missing directions
//   + track * |union of blocked tracks|
//   + spread * sum over endpoint pairs on one layer of
//       max(0, spread_distance - l1 distance)
//   + unassigned * pins without a candidate.
// Returns -1 if two chosen candidates conflict.
int64_t selectionCost(const SelectionInput& in, const std::vector<int>& choice);

// Destructive bounding: for a hypothetical bound F (initially the sum of
// the pins' individual lower bounds), enumerate assignments depth-first,
// removing from every open pin's domain the candidates that conflict with
// the assignment or whose lower bound exceeds F; an empty domain refutes
// the node. If F is refuted, F becomes the smallest lower bound that
// exceeded it. The first complete assignment found has minimum cost.
SelectionResult selectConflictFree(const SelectionInput& in);

// Greedy reference (the previous behaviour): pins in order, each takes its
// first candidate compatible with the earlier choices.
SelectionResult selectGreedy(const SelectionInput& in);

}  // namespace drt::br
