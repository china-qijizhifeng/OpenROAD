// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "GCellGraph.h"

namespace grt::helem {

// Edge of a rectilinear Steiner tree over GCell coordinates. Endpoints may be
// pins or Steiner points; they are connected by pattern routing.
struct SteinerEdge
{
  GCellPoint a;
  GCellPoint b;
};

// Builds a Steiner tree over distinct GCell points. Must be thread safe.
using SteinerBuilder
    = std::function<std::vector<SteinerEdge>(const std::vector<GCellPoint>&)>;

using HelemLogger = std::function<void(const std::string&)>;

// Wire weight w_e used by the 2D LR/LEM stages (paper gap G10).
enum class WireWeightMode
{
  kNormalizedLength,  // edge length / mean edge length
  kScoreLength,       // unit_wire_cost * edge length (ISPD 2024 score units)
};

// Access point choice per pin (paper gap G1).
enum class PinAccessMode
{
  kNetCenter,     // closest to the centre of all access points of the net
  kOtherPinsBox,  // closest to the box of the other pins' first choices
};

// Layer assignment algorithm (paper gap G8).
enum class LayerAssignmentMode
{
  kPath,  // per 2-pin segment DP, segments in BFS order (first version)
  kTree,  // exact DP over the whole Steiner tree, via spans chosen per node
};

// Order of nets inside the layer assignment (paper gap G8).
enum class LayerAssignmentOrder
{
  kShortFirst,
  kLongFirst,
};

struct HelemOptions
{
  // Lagrangian relaxation stage with L-shape routing (paper Section 3.2).
  int lr_iterations = 8;
  // Linearized exponential multiplier stage (paper Algorithm 1).
  int lem_iterations = 3;  // k_max
  int bcd_passes = 1;      // p_max
  double rho_initial = 0.05;
  double rho_growth = 2.0;  // sigma
  double tau = 100.0;
  // y_e^0 (tuned on ISPD24, R8); a negative value selects the paper default
  // (1/(100 L)) sum OvflW.
  double multiplier_initial = 0.2;
  // Segments per BCD block (paper gap G5).
  int batch_size = 2048;
  // 3-bend search box expansion: each side grows by
  // max(bbox_expand_min, ceil(bbox_expand_ratio * side length)) (gap G6).
  double bbox_expand_ratio = 0.25;
  int bbox_expand_min = 1;
  // Upper bound on the expansion so huge nets keep a bounded search box.
  int bbox_expand_max = 32;
  // Cost of one bend in the 3-bend kernel, in units of the via cost
  // expressed as wire weight (gap G13; the paper charges no bend cost).
  double bend_cost_scale = 0.0;
  // Include the segment's current route in g_e when pricing (gap G7).
  bool price_includes_own_route = true;
  WireWeightMode wire_weight_mode = WireWeightMode::kNormalizedLength;
  PinAccessMode pin_access_mode = PinAccessMode::kOtherPinsBox;
  // Among equally close access points prefer the highest layer.
  bool pin_access_prefer_high_layer = false;
  // Scale of the projected 2D capacity (gap G3, 1.0 = sum of layers).
  double capacity_2d_scale = 0.8;
  // Layer assignment.
  int layer_assignment_batch_size = 16;
  LayerAssignmentOrder layer_assignment_order
      = LayerAssignmentOrder::kShortFirst;
  LayerAssignmentMode layer_assignment_mode = LayerAssignmentMode::kTree;
  // Rip-up and reassign passes after the first layer assignment.
  int layer_assignment_passes = 1;
  int threads = 1;
  bool verbose = true;
  // Recompute the demand ledger from scratch after each stage and compare.
  bool verify_ledger = false;
};

struct HelemStageReport
{
  std::string stage;
  double seconds = 0.0;
  int64_t total_overflow_2d = 0;  // sum over 2D edges of max(0, demand - cap)
  int max_overflow_2d = 0;
  int64_t overflowed_edges_2d = 0;
};

struct HelemScore
{
  double wire_cost = 0.0;
  double via_cost = 0.0;
  double overflow_cost = 0.0;
  int64_t via_count = 0;
  double total() const { return wire_cost + via_cost + overflow_cost; }
};

// One access point per pin (nullptr for pins without access points), using
// the pin access rule of gap G1. Shared with the comparison drivers so every
// router sees the same pins.
std::vector<const AccessPoint*> choosePinAccessPoints(const RoutingNet& net,
                                                      PinAccessMode mode,
                                                      bool prefer_high);

class HelemRouter
{
 public:
  HelemRouter(const GCellGraph& graph,
              const std::vector<RoutingNet>& nets,
              const HelemOptions& options,
              SteinerBuilder steiner_builder,
              HelemLogger logger = nullptr);
  ~HelemRouter();

  // Runs preparation, 2D routing and layer assignment. Returns one route per
  // input net, in input order.
  std::vector<NetRoute3D> run();

  const std::vector<HelemStageReport>& stageReports() const;
  // Score of the final 3D solution, computed with the ISPD 2024 evaluator
  // rules (wire, via and overflow cost).
  HelemScore score() const;

 private:
  class Impl;
  Impl* impl_;
};

}  // namespace grt::helem
