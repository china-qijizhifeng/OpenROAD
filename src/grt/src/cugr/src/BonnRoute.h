// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// BonnRoute global routing (Gester et al., "BonnRoute: Algorithms and Data
// Structures for Fast and Good VLSI Routing", TODAES 18(2), 2013, Sec. 2):
// min-max resource sharing (Algorithm 2) over the CUGR 3D grid graph with a
// path-composition Steiner oracle (Algorithm 1), followed by randomized
// rounding and re-choosing among the generated candidates (Sec. 2.4).
//
// Resources are the CUGR wire edges (demand as committed by
// GridGraph::addTreeUsage, including via-pad flank demand) plus one objective
// resource (wire length + vias, CUGR cost weights). Edge capacities follow
// Sec. 2.5 (BonnRouteCapacity.h) or, with capacity_model 0, CUGR's model.
// Pins are vertex sets V_p joined by zero-cost cliques (Sec. 2.1-2.2); the
// trees keep the used clique edges as pin links (GRTreeNode::isPinLink).
// The routing trees left on the nets are committed to the grid graph so that
// every CUGR downstream step (maze clean-up, guides, congestion report)
// applies as is.

#pragma once

#include <array>
#include <atomic>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "BonnRouteCapacity.h"
#include "GRTree.h"
#include "geo.h"

namespace utl {
class Logger;
}

namespace grt {

class GridGraph;
class GRNet;

// The estimator input of Sec. 2.5 (BonnRouteCapacity.h) for CUGR's grid,
// layers and blockages; capacities come back in GridGraph::edgeFlatIndex
// order.
bonnroute::CapacityGrid capacityGridFromCugr(const GridGraph& grid_graph);

struct BonnRouteParams
{
  // Algorithm 2: number of phases t and price update parameter epsilon.
  int phases = 125;
  double epsilon = 1.0;
  // A net keeps its previous solution while its current price is within
  // this factor of a certified lower bound of the current oracle optimum
  // (see the reuse certificate in resourceSharing()).
  double reuse_tolerance = 0.1;
  // Objective resource bound: (1 + slack) * congestion-free objective.
  double objective_slack = 0.05;
  // Objective weights, CUGR units: per M2 pitch of wire / per via.
  double wire_weight = 0.5;
  double via_weight = 4.0;
  // Oracle search region: net bbox expanded by margin + ratio * max(w, h).
  int margin = 3;
  double margin_ratio = 0.25;
  // Capacity floor used to normalise resource usage of blocked edges.
  double min_capacity = 0.25;
  // Nets processed with frozen prices per parallel batch. Fixed size keeps
  // results independent of the thread count.
  int batch_size = 128;
  int threads = 1;
  // Randomized rounding seed and re-choose passes.
  unsigned seed = 1;
  int repair_passes = 8;
  // Pins as vertex sets V_p (Sec. 2.1): the cells of the pin's access
  // points (shape cells for a pin without any); false = one pre-selected
  // cell per pin.
  bool pin_sets = true;
  // Edge capacities: 1 = Sec. 2.5 estimate, 0 = CUGR's capacity model.
  int capacity_model = 1;
  // Sec. 2.5 refinements of capacity_model 1.
  bool capacity_single_tile = true;
  bool capacity_intra_tile = true;
  bool capacity_stacked_vias = true;
  double capacity_blockage_extension = 1.0;
  int capacity_single_tile_layers = 2;
  double capacity_single_tile_margin = 0.25;
  // Oracle on the whole grid graph instead of the expanded net bbox.
  bool full_graph = false;
  // Run CUGR maze rip-up and reroute on nets still on overflowed edges.
  bool maze_cleanup = true;
  // Sec. 2.2: the oracle's Dijkstra uses goal orientation (A*) with a
  // consistent lower bound: l1 distance times the cheapest edge and via
  // costs, and ALT bounds from `landmarks` landmarks [Goldberg and
  // Harrelson 2005] recomputed every `landmark_refresh` phases.
  int goal_orientation = 1;
  int landmarks = 8;
  int landmark_refresh = 25;
  // Debug: rerun every oracle without goal orientation and compare the
  // connection costs (reported as "goal check").
  bool goal_check = false;
  // Sec. 2.3/5.1: threads route nets concurrently and update the shared
  // prices as they go (volatility-tolerant block solver, Mueller et al.
  // 2011), without frozen batches or a barrier per batch; one barrier per
  // phase. Results then depend on thread timing. 0 = deterministic
  // frozen-price batches.
  bool concurrent_prices = false;
  bool verbose = true;

  // Returns false if the name is unknown.
  bool set(const std::string& name, double value);
};

class BonnRouteGlobal
{
 public:
  BonnRouteGlobal(GridGraph* grid_graph,
                  std::vector<GRNet*> nets,
                  const BonnRouteParams& params,
                  utl::Logger* logger);
  ~BonnRouteGlobal();

  // Routes every net; on return each net has a routing tree whose usage is
  // committed to the grid graph.
  void run();
  // CUGR's capacities (GridGraph::edgeFlatIndex order) before capacity
  // model 1 replaced them; empty for capacity model 0.
  const std::vector<double>& cugrCapacity() const { return cugr_capacity_; }

 private:
  // A pin is the set V_p of grid vertices of its shapes (Sec. 2.1); the
  // oracle may connect it at any of them (zero-cost clique K(V_p)).
  struct Terminal
  {
    std::vector<std::array<int, 3>> cells;  // (layer, x, y)
    int pin_index = -1;
  };
  struct NodeRec
  {
    int16_t layer;
    int16_t x;
    int16_t y;
    int32_t parent;  // index into the node list, -1 for the root
    // GRNet pin index if the edge from the parent is a pin link, else -1.
    int32_t pin = -1;
  };
  struct Candidate
  {
    std::vector<NodeRec> nodes;                      // preorder
    std::vector<std::pair<uint32_t, double>> usage;  // resource -> demand
    double objective = 0.0;
    uint64_t hash = 0;
    int count = 0;
    // Representative cell of each terminal (same order as
    // NetState::terminals): a cell of the pin in the tree, preferring one
    // with a wire or via.
    std::vector<std::array<int16_t, 3>> pin_cells;
  };
  struct NetState
  {
    GRNet* net = nullptr;
    std::vector<Terminal> terminals;
    BoxT region;
    std::vector<Candidate> candidates;
    int last = -1;
    double last_log_cost = 0.0;  // log of the oracle value (absolute prices)
    int chosen = -1;
    bool trivial = false;
    // Single-cell net: the pin stack at `trivial_point`.
    PointT trivial_point;
    int trivial_low = 0;
    int trivial_top = 0;
    // ln of the congestion-free objective of the phase-0 oracle; together
    // with the objective price it bounds opt_n(y) from below (reuse
    // certificate, see resourceSharing()).
    double log_base_objective = -1e300;
    int oracle_failures = 0;
  };
  struct Scratch;
  // Cost and path signature of one Algorithm 1 connection (goal_check).
  struct Segment
  {
    double cost = 0.0;
    uint64_t path = 0;
  };
  struct GoalEntry
  {
    double key;  // dist + h
    double dist;
    int64_t node;
    uint32_t version;  // target set the key's h was computed for
  };
  struct GoalTarget
  {
    int l;
    int x;
    int y;
    int64_t vertex;  // global (layer, x, y) index for the landmarks
  };
  // Lower bound h of one oracle call (see goalBound()).
  struct GoalContext
  {
    // Full cost per unit of |dx|, |dy| and per via; objective part alone.
    double cx = 0.0;
    double cy = 0.0;
    double cv = 0.0;
    double ox = 0.0;
    double oy = 0.0;
    double ov = 0.0;
    const uint8_t* via_steps = nullptr;  // [mask][from][to]
    std::vector<GoalTarget> targets;
    bool exact = true;
    int lx = 0;
    int hx = 0;
    int ly = 0;
    int hy = 0;
    int ll = 0;
    int hl = 0;
    // ALT: landmark distances of the targets per landmark, as an interval.
    const float* alt = nullptr;
    double alt_scale = 0.0;
    std::vector<float> alt_lo;
    std::vector<float> alt_hi;
    uint32_t version = 0;
  };
  struct OracleResult
  {
    bool ok = false;
    Candidate candidate;
    double log_cost = 0.0;
    double search_cost = 0.0;  // Dijkstra length of the last connection
    double exact_cost = 0.0;
    bool checked = false;
    bool mismatch = false;
  };

  void prepareNets();
  // Path composition Steiner tree (Algorithm 1) under the current prices.
  OracleResult oracle(const NetState& state,
                      Scratch& scratch,
                      bool objective_only,
                      bool whole_grid,
                      std::vector<Segment>* segments = nullptr) const;
  template <bool kGoal>
  OracleResult oracleImpl(const NetState& state,
                          Scratch& scratch,
                          bool objective_only,
                          bool whole_grid,
                          std::vector<Segment>* segments) const;
  void initGoal(GoalContext& goal,
                const NetState& state,
                bool objective_only,
                int x0,
                int y0,
                int x1,
                int y1,
                int num_layers) const;
  void setGoalTargets(GoalContext& goal,
                      const NetState& state,
                      const std::vector<char>& connected,
                      int num_layers) const;
  double goalBound(const GoalContext& goal, int l, int x, int y) const;
  // Phase-start lower bounds: cheapest wire price per direction and the
  // landmark distances of the congestion part of the costs.
  void updateGoalBounds(bool landmarks);
  void landmarkDijkstra(int64_t source, std::vector<double>& dist) const;
  double viaPriceBound(int l, int x, int y) const;
  void buildViaSteps();
  // Debug comparison of a goal-oriented oracle call with plain Dijkstra.
  void checkGoal(const NetState& state,
                 Scratch& reference,
                 bool objective_only,
                 bool whole_grid,
                 bool ok,
                 const std::vector<Segment>& segments);
  // tree_edges holds vertex pairs; edge_pins[i] is the GRNet pin index of
  // pair i if it is a pin-clique edge, else -1.
  // Without need_usage (the congestion-free pilot reads only the tree and
  // its objective) the demand list and hash stay empty.
  bool buildCandidate(const NetState& state,
                      const std::vector<std::array<int, 3>>& tree_edges,
                      const std::vector<int>& edge_pins,
                      bool need_usage,
                      Candidate& candidate) const;
  // Congestion-free oracle pass: objective bound and reuse certificate;
  // with stack_counts, the via stacks through each (layer, x, y).
  void pilot(std::vector<double>* stack_counts);
  void estimateCapacities();
  void applyCapacities(const std::vector<double>& capacities);
  void reportCapacities() const;
  double padCost(int layer, int x, int y, double pad, const GRNet* net) const;
  double padCostFactor(int layer,
                       int x,
                       int y,
                       double pad,
                       double factor) const;
  std::shared_ptr<GRTreeNode> buildTree(const Candidate& candidate) const;
  double logCost(const Candidate& candidate) const;
  double linearCost(const Candidate& candidate) const;
  // Adds the candidate's usage to the log prices; the weights of the touched
  // resources are refreshed once per batch by refreshWeights().
  void applyPrices(const Candidate& candidate);
  void refreshWeights();
  // Concurrent variant: atomic log-price adds and monotone weight stores.
  void applyPricesConcurrent(const Candidate& candidate);
  // Price reads that tolerate concurrent updates (plain loads on x86).
  double price(uint32_t r) const;
  double objectiveWeight() const;
  void resourceSharingConcurrent(double log_tolerance,
                                 double initial_log_objective);
  void renormalize();
  int insertCandidate(NetState& state, Candidate&& candidate);
  void resourceSharing();
  void roundAndRepair();
  void commit();

  // Resource indexing over the CUGR wire edges.
  uint32_t edgeIndex(int layer, int x, int y) const
  {
    return (static_cast<uint32_t>(layer) * x_size_ + x) * y_size_ + y;
  }
  double wireObjective(int layer, int x, int y) const;
  double wireCost(int layer, int x, int y, double factor) const;

  GridGraph* grid_graph_;
  std::vector<GRNet*> nets_;
  BonnRouteParams params_;
  utl::Logger* logger_;

  int x_size_ = 0;
  int y_size_ = 0;
  int num_layers_ = 0;
  double m2_pitch_ = 1.0;

  std::vector<double> capacity_;       // capacity minus background demand
  std::vector<double> norm_capacity_;  // max(capacity_, min_capacity)
  std::vector<double> log_price_;      // ln y_r (absolute)
  std::vector<double> weight_;         // exp(ln y_r - shift_) / norm_cap
  double log_objective_price_ = 0.0;
  double objective_bound_ = 1.0;
  double objective_weight_ = 0.0;  // exp(ln y_obj - shift_) / bound
  double shift_ = 0.0;
  double max_log_price_ = 0.0;
  // Per-layer via pad lengths of a single via starting at layer l
  // (lower pad on l, upper pad on l + 1).
  std::vector<double> via_pad_lower_;
  std::vector<double> via_pad_upper_;
  // Static grid geometry of the cost functions, per layer: preferred
  // direction, wire objective per edge index, and the flank edges a via pad
  // at index c spreads over (lengths of edges c - 1 and c, and their sum).
  struct PadFlank
  {
    int sum = 0;
    bool lower = false;
    bool higher = false;
  };
  std::vector<int> layer_direction_;
  std::array<std::vector<double>, 2> wire_objective_;  // [direction][index]
  std::vector<std::vector<PadFlank>> pad_flank_;       // [layer][index]
  // Resources whose log price changed in the current batch.
  std::vector<uint32_t> touched_;
  std::vector<uint32_t> touched_stamp_;
  uint32_t touched_epoch_ = 0;
  bool via_demand_ = true;
  // Layers below this hold no routing capacity (pin layers): not priced.
  int min_layer_ = 0;

  std::vector<NetState> states_;
  std::unique_ptr<bonnroute::CapacityEstimator> estimator_;
  std::vector<double> cugr_capacity_;
  double estimate_time_ = 0.0;
  double base_objective_ = 0.0;
  std::vector<int> order_;  // nets taking part in resource sharing

  // Goal orientation (Sec. 2.2).
  std::vector<uint8_t> via_steps_;              // [lo][hi][mask][from][to]
  std::vector<double> ndr_min_;                 // per layer, over the nets
  std::array<double, 2> min_weight_{0.0, 0.0};  // per direction
  std::vector<int64_t> landmark_vertices_;
  std::vector<float> landmark_dist_;  // vertex * landmarks + k
  double landmark_shift_ = 0.0;
  bool landmarks_valid_ = false;
  int landmark_updates_ = 0;
  double landmark_time_ = 0.0;
  int64_t search_pops_ = 0;
  int64_t search_settled_ = 0;
  std::atomic<int64_t> goal_checked_{0};
  std::atomic<int64_t> goal_segments_{0};
  std::atomic<int64_t> goal_ties_{0};
  std::atomic<int64_t> goal_mismatches_{0};

  // Statistics.
  int64_t oracle_calls_ = 0;
  int64_t reuses_ = 0;
  int64_t cost_checks_ = 0;
  int64_t cost_mismatches_ = 0;
};

}  // namespace grt
