// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// BonnRoute global routing (Gester et al., "BonnRoute: Algorithms and Data
// Structures for Fast and Good VLSI Routing", TODAES 18(2), 2013, Sec. 2):
// min-max resource sharing (Algorithm 2) over the CUGR 3D grid graph with a
// path-composition Steiner oracle (Algorithm 1), followed by randomized
// rounding and re-choosing among the generated candidates (Sec. 2.4).
//
// Resources are the CUGR wire edges (capacity in tracks, demand as committed
// by GridGraph::addTreeUsage, including via-pad flank demand) plus one
// objective resource (wire length + vias, CUGR cost weights). The routing
// trees left on the nets are committed to the grid graph so that every CUGR
// downstream step (maze clean-up, guides, congestion report) applies as is.

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "GRTree.h"
#include "geo.h"

namespace utl {
class Logger;
}

namespace grt {

class GridGraph;
class GRNet;

struct BonnRouteParams
{
  // Algorithm 2: number of phases t and price update parameter epsilon.
  int phases = 125;
  double epsilon = 1.0;
  // A net keeps its previous solution while its price grew by at most this
  // factor since the oracle produced it (prices never decrease, so the
  // previous oracle value is a lower bound of the current optimum).
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
  // Run CUGR maze rip-up and reroute on nets still on overflowed edges.
  bool maze_cleanup = true;
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

 private:
  struct Terminal
  {
    PointT point;
    int low_layer;  // pin layer
    int top_layer;  // layer the search connects to
  };
  struct NodeRec
  {
    int16_t layer;
    int16_t x;
    int16_t y;
    int32_t parent;  // index into the node list, -1 for the root
  };
  struct Candidate
  {
    std::vector<NodeRec> nodes;                     // preorder
    std::vector<std::pair<uint32_t, float>> usage;  // resource -> demand
    double objective = 0.0;
    uint64_t hash = 0;
    int count = 0;
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
    int oracle_failures = 0;
  };
  struct Scratch;
  struct OracleResult
  {
    bool ok = false;
    Candidate candidate;
    double log_cost = 0.0;
  };

  void prepareNets();
  // Path composition Steiner tree (Algorithm 1) under the current prices.
  OracleResult oracle(const NetState& state,
                      Scratch& scratch,
                      bool objective_only,
                      bool whole_grid) const;
  bool buildCandidate(const NetState& state,
                      const std::vector<std::array<int, 3>>& tree_edges,
                      Candidate& candidate) const;
  std::shared_ptr<GRTreeNode> buildTree(const Candidate& candidate) const;
  double logCost(const Candidate& candidate) const;
  void applyPrices(const Candidate& candidate);
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
  double viaCost(int layer, int x, int y, const GRNet* net) const;

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
  bool via_demand_ = true;

  std::vector<NetState> states_;
  std::vector<int> order_;  // nets taking part in resource sharing

  // Statistics.
  int64_t oracle_calls_ = 0;
  int64_t reuses_ = 0;
};

}  // namespace grt
