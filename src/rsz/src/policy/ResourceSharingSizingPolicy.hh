// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026-2026, The OpenROAD Authors

#pragma once

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include "OptimizationPolicy.hh"
#include "sta/GraphClass.hh"
#include "sta/LibertyClass.hh"
#include "sta/NetworkClass.hh"

namespace sta {
class Pvt;
class Scene;
class TimingArc;
}  // namespace sta

namespace rsz {

// Settings of the RS_SIZING phase.  Values come from the environment variable
// RSZ_RS_SIZING_OPTIONS ("name=value name=value ..."), read once per run.
struct ResourceSharingSizingConfig
{
  // How arc weights are found.
  enum class Multipliers
  {
    // Resource sharing: multiplicative path weights exp(gamma * usage),
    // summed over all paths through an arc by forward/backward passes.
    kResourceSharing,
    // Lagrangian relaxation: per-arc multipliers updated multiplicatively by
    // the arc's local criticality, then projected onto flow conservation.
    kLagrangian,
  };
  Multipliers multipliers = Multipliers::kResourceSharing;
  // Lagrangian update: lambda *= max(floor, 1 + criticality) ^ exponent,
  // where criticality = -(slack through the arc) / period.
  double lagrangian_exponent = 1.0;
  double lagrangian_floor = 0.05;
  // Lagrangian power multiplier relative to an endpoint multiplier of one.
  double lagrangian_power_weight = 1.0;
  // Resource-sharing iterations (weight computation + sizing oracle).
  int iterations = 6;
  // Oracle rounds per iteration.  Gates that share a net with a gate resized
  // in the same round wait for the next round.
  int rounds_per_iteration = 3;
  // Multiplicative-weights exponent.  A path whose usage (delay / period) is
  // larger by 1/gamma gets e times the weight.
  double gamma = 80.0;
  // Weights are computed from the mean of the per-iteration log weights
  // (multiplicative-weights accumulation) when true, from the current
  // iteration only when false.
  bool accumulate_weights = true;
  // Power budget as a fraction of the initial power proxy of the sized gates.
  double power_budget_ratio = 1.0;
  // Weight of the leakage share in the power proxy; area has weight 1.
  double leakage_weight = 1.0;
  // A cell change is taken only when it lowers the oracle objective by more
  // than this fraction of the gate's weighted delay sum.
  double min_relative_gain = 0.01;
  // Sibling loads considered per input net.
  int max_siblings = 24;
  // Size sequential cells (their clock pins load the clock tree).
  bool size_sequential = false;
  // Iteration acceptance: the iteration is kept when TNS improves and WNS
  // loses at most this much (seconds), or WNS improves and TNS loses at most
  // tns_loss_fraction of |TNS|.
  double wns_loss_limit = 1e-12;
  double tns_loss_fraction = 0.0;
  // Debug: print the weight self-check against STA slacks.
  bool check_weights = false;
  // The phase starts only when the mean violation of the violating endpoints
  // is at least this fraction of the clock period.  Smaller violations are
  // closed by a few local moves.
  double start_violation_ratio = 0.15;
  // Gates whose largest own-arc weight is below this fraction of the largest
  // arc weight in the design are left out of the oracle.
  double critical_weight_ratio = 1e-3;
  // An iteration is undone, and the phase stops, when its TNS gain per power
  // increase is below this fraction of the first kept iteration's.
  double min_efficiency_ratio = 0.5;
};

// RS_SIZING: gate sizing (and Vt choice among swappable cells) by the
// resource-sharing algorithm of Daboul, Hähnle, Held and Schorr ("Provably
// fast and near-optimum gate sizing", TCAD 2018).
//
// Resources are the power budget and every maximal timing path; the usage of a
// path is its delay over the clock period.  Path weights exp(gamma * usage)
// are never enumerated: one forward and one backward log-sum-exp pass over the
// timing graph give, for every timing arc, the summed weight of all paths
// through it (edge weights fulfil flow conservation by construction).
//
// The oracle visits every sizable gate and chooses the cell minimizing
//   w_power * dPower / budget + sum_arcs w_arc * dDelay_arc / period
// over the arcs whose delay the cell choice changes: the gate's own arcs, the
// arcs of the input-net drivers (load change), the arcs of the sibling loads
// (driver slew change) and the arcs of the fanout gates (slew change).  Delay
// changes come from the Liberty table models at the current slews and loads,
// so the oracle reads frozen snapshots and runs on all worker threads.
class ResourceSharingSizingPolicy : public OptimizationPolicy
{
 public:
  ResourceSharingSizingPolicy(Resizer& resizer,
                              MoveCommitter& committer,
                              RepairSetupContext& setup_context,
                              const OptimizerRunConfig& config);
  ~ResourceSharingSizingPolicy() override;

  const char* name() const override { return "ResourceSharingSizingPolicy"; }
  bool start() override;
  void iterate() override;

  // One Liberty arc evaluated with its table model at a frozen operating
  // point.  weight is the summed path weight of the arc divided by the period.
  struct TableArc
  {
    const sta::TimingArc* arc{nullptr};
    const sta::Pvt* pvt{nullptr};
    float input_slew{0.0f};
    float load_cap{0.0f};
    double weight{0.0};
    int input_rf{0};
    int output_rf{0};
    float delay{0.0f};
    float output_slew{0.0f};
  };

  // An input pin of the gate: the arcs of its net driver, which see the input
  // capacitance as load, and the arcs of the other loads of the net, which
  // see the driver's output slew.
  struct InputPinContext
  {
    const sta::LibertyPort* port{nullptr};
    float current_cap{0.0f};
    std::vector<TableArc> driver_arcs;
    std::vector<TableArc> sibling_arcs;
  };

  // An output pin of the gate and the arcs of its fanout gates.
  struct OutputPinContext
  {
    const sta::LibertyPort* port{nullptr};
    float load_cap{0.0f};
    // Smallest max-slew limit on the net; infinity without a limit.
    float slew_limit{0.0f};
    // Largest slew increase from the driver pin to a load pin.
    float wire_slew_degradation{0.0f};
    // Largest current output slew per output transition.
    std::array<float, 2> current_slew{0.0f, 0.0f};
    std::vector<TableArc> fanout_arcs;
  };

  // An arc of the gate itself.
  struct GateArc
  {
    TableArc table_arc;
    int input_index{-1};
    int output_index{-1};
  };

  struct CandidateCell
  {
    sta::LibertyCell* cell{nullptr};
    double power{0.0};
    // Matching arcs aligned with GateSnapshot::gate_arcs.
    std::vector<const sta::TimingArc*> gate_arcs;
    // Input capacitances aligned with GateSnapshot::inputs.
    std::vector<float> input_caps;
    // Max-capacitance limits aligned with GateSnapshot::outputs (infinity
    // without a limit).
    std::vector<float> output_cap_limits;
  };

  struct GateSnapshot
  {
    sta::Instance* instance{nullptr};
    sta::LibertyCell* current_cell{nullptr};
    double current_power{0.0};
    std::vector<InputPinContext> inputs;
    std::vector<OutputPinContext> outputs;
    std::vector<GateArc> gate_arcs;
    std::vector<CandidateCell> candidates;
    // Sum of weight * delay over the gate's own arcs (gain normalization).
    double weighted_delay{0.0};
  };

  struct GateDecision
  {
    sta::Instance* instance{nullptr};
    sta::LibertyCell* cell{nullptr};  // nullptr: keep the current cell
    double gain{0.0};                 // objective decrease (positive)
  };

 private:
  struct IterationResult
  {
    float wns{0.0f};
    float tns{0.0f};
    double power{0.0};
    int resized{0};
  };

  void readConfig();
  void collectSizableGates();
  double clockPeriod() const;
  // Vertices in increasing level order.
  void buildLevelOrder();
  // Forward/backward log-sum-exp passes; fills log_forward_ / log_backward_.
  void computeLogWeights();
  // Multiplicative update and flow projection of the arc multipliers.
  void updateLagrangeMultipliers(int iteration);
  // Arc weight (summed path weight or multiplier) divided by the period.
  double arcWeight(sta::Edge* edge,
                   const sta::TimingArc* arc,
                   sta::Vertex* from,
                   sta::Vertex* to,
                   float delay) const;
  // Power weight on the same scale as arcWeight() * period.
  double powerWeight() const;
  void accumulateLogWeights(int iteration);
  double arcLogWeight(sta::Vertex* from,
                      int from_rf,
                      sta::Vertex* to,
                      int to_rf,
                      float delay) const;
  bool snapshot(sta::Instance* instance, GateSnapshot& snap);
  bool fillTableArc(sta::Edge* edge,
                    const sta::TimingArc* arc,
                    sta::Vertex* from,
                    sta::Vertex* to,
                    float load_cap,
                    TableArc& table_arc) const;
  GateDecision evaluate(const GateSnapshot& snap) const;
  double objective(const GateSnapshot& snap,
                   const CandidateCell& candidate,
                   bool& feasible) const;
  int runRound(std::vector<sta::Instance*>& gates);
  double powerOf(sta::LibertyCell* cell);
  double designPower();
  const sta::Pvt* pvtOf(const sta::Instance* instance) const;
  void checkWeights();
  // Largest weight (times period) of the arcs into the gate's outputs.
  double gateCriticality(sta::Instance* instance) const;
  // Scale the Lagrange multipliers so that the largest arc multiplier is one.
  void normalizeLagrangeMultipliers();
  // Number of endpoints with negative setup slack.
  int violatingEndpointCount();

  ResourceSharingSizingConfig rs_config_;
  const sta::Scene* scene_{nullptr};
  int dcalc_ap_{0};
  double period_{0.0};
  double power_budget_{0.0};
  double log_power_weight_{0.0};
  double log_weight_shift_{0.0};
  double area_scale_{0.0};
  double leakage_scale_{0.0};
  std::vector<sta::Instance*> sizable_gates_;
  // Per vertex id and rise/fall: log of the summed weight of all partial
  // paths from a start to the vertex (forward) and from the vertex to an
  // end (backward).  The *_mean_ arrays hold the iteration means.
  std::vector<double> log_forward_;
  std::vector<double> log_backward_;
  std::vector<double> log_forward_mean_;
  std::vector<double> log_backward_mean_;
  std::vector<sta::Vertex*> level_order_;
  // Lagrangian multipliers per edge id and (from rise/fall, to rise/fall),
  // and per endpoint vertex id and rise/fall.
  std::vector<double> arc_multipliers_;
  std::vector<double> endpoint_multipliers_;
  std::unordered_map<sta::LibertyCell*, double> power_cache_;
  // Wall seconds per step, reported at the end of the phase.
  double seconds_weights_{0.0};
  double seconds_snapshots_{0.0};
  double seconds_oracle_{0.0};
  double seconds_apply_{0.0};
  double seconds_timing_{0.0};
};

}  // namespace rsz
