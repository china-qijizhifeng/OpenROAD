// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026-2026, The OpenROAD Authors

#include "ResourceSharingSizingPolicy.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "OptimizerTypes.hh"
#include "Rebuffer.hh"
#include "RebufferPlannerPool.hh"
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "est/EstimateParasitics.h"
#include "odb/db.h"
#include "rsz/Resizer.hh"
#include "sta/Clock.hh"
#include "sta/Delay.hh"
#include "sta/Graph.hh"
#include "sta/GraphDelayCalc.hh"
#include "sta/Liberty.hh"
#include "sta/MinMax.hh"
#include "sta/Mode.hh"
#include "sta/Network.hh"
#include "sta/PortDirection.hh"
#include "sta/Scene.hh"
#include "sta/Sdc.hh"
#include "sta/Search.hh"
#include "sta/TableModel.hh"
#include "sta/TimingArc.hh"
#include "sta/TimingRole.hh"
#include "sta/Transition.hh"
#include "utl/Logger.h"
#include "utl/ThreadPool.h"

namespace rsz {

using utl::RSZ;

namespace {

constexpr double kNegativeInfinity = -std::numeric_limits<double>::infinity();
constexpr float kFloatInfinity = std::numeric_limits<float>::infinity();

double secondsSince(const std::chrono::steady_clock::time_point start)
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now()
                                       - start)
      .count();
}

double logSumExp(const double a, const double b)
{
  if (a == kNegativeInfinity) {
    return b;
  }
  if (b == kNegativeInfinity) {
    return a;
  }
  const double high = std::max(a, b);
  return high + std::log1p(std::exp(std::min(a, b) - high));
}

bool isPropagatedArcSet(const sta::Edge* edge)
{
  const sta::TimingRole* role = edge->role();
  if (role != nullptr && role->isTimingCheck()) {
    return false;
  }
  if (edge->isDisabledLoop()) {
    return false;
  }
  return role != sta::TimingRole::latchDtoQ();
}

bool tableDelay(const sta::TimingArc* arc,
                const sta::Pvt* pvt,
                const float input_slew,
                const float load_cap,
                float& delay,
                float& output_slew)
{
  const sta::GateTableModel* model
      = dynamic_cast<const sta::GateTableModel*>(arc->model());
  if (model == nullptr) {
    return false;
  }
  output_slew = 0.0f;
  model->gateDelay(pvt, input_slew, load_cap, delay, output_slew);
  return true;
}

int riseFallIndex(const sta::Transition* transition)
{
  const sta::RiseFall* rf = transition->asRiseFall();
  return rf == nullptr ? 0 : static_cast<int>(rf->index());
}

}  // namespace

ResourceSharingSizingPolicy::ResourceSharingSizingPolicy(
    Resizer& resizer,
    MoveCommitter& committer,
    RepairSetupContext& setup_context,
    const OptimizerRunConfig& config)
    : OptimizationPolicy(resizer, committer, setup_context, config)
{
}

ResourceSharingSizingPolicy::~ResourceSharingSizingPolicy() = default;

double ResourceSharingSizingPolicy::netCriticality(
    sta::Vertex* driver_vertex) const
{
  double criticality = 0.0;
  sta::VertexOutEdgeIterator edges(driver_vertex, graph_);
  while (edges.hasNext()) {
    sta::Edge* edge = edges.next();
    if (!edge->isWire()) {
      continue;
    }
    for (const sta::TimingArc* arc : edge->timingArcSet()->arcs()) {
      const float delay = sta::delayAsFloat(
          graph_->arcDelay(edge, const_cast<sta::TimingArc*>(arc), dcalc_ap_));
      criticality = std::max(
          criticality,
          arcWeight(edge, arc, driver_vertex, edge->to(graph_), delay)
              * period_);
    }
  }
  return criticality;
}

int ResourceSharingSizingPolicy::runBufferRound()
{
  const auto start = std::chrono::steady_clock::now();
  sta_->findRequireds();
  struct BufferJob
  {
    const sta::Pin* driver_pin{nullptr};
    const sta::Pin* target_load_pin{nullptr};
    BufferedNetPtr net;
    double criticality{0.0};
  };
  // Candidate nets: violating, critical, small fanout, bufferable.
  std::vector<BufferJob> jobs;
  for (sta::Vertex* vertex : level_order_) {
    const sta::Pin* pin = vertex->pin();
    if (!network_->isDriver(pin) || vertex->isBidirectDriver()) {
      continue;
    }
    const sta::Instance* instance = network_->instance(pin);
    if (network_->isTopInstance(instance)
        || network_->libertyCell(instance) == nullptr
        || sta_->isClock(pin, sta_->cmdMode())) {
      continue;
    }
    if (sta::delayAsFloat(sta_->slack(vertex, max_)) >= 0.0f) {
      continue;
    }
    const double criticality = netCriticality(vertex);
    if (criticality < rs_config_.buffer_weight_ratio) {
      continue;
    }
    int fanout = 0;
    const sta::Pin* worst_load = nullptr;
    float worst_slack = std::numeric_limits<float>::infinity();
    sta::VertexOutEdgeIterator edges(vertex, graph_);
    while (edges.hasNext()) {
      sta::Edge* edge = edges.next();
      if (!edge->isWire()) {
        continue;
      }
      ++fanout;
      sta::Vertex* load = edge->to(graph_);
      const float slack = sta::delayAsFloat(sta_->slack(load, max_));
      if (slack < worst_slack) {
        worst_slack = slack;
        worst_load = load->pin();
      }
    }
    if (fanout < 1 || fanout >= rs_config_.buffer_max_fanout
        || worst_load == nullptr || !resizer_.okToBufferNet(pin)) {
      continue;
    }
    jobs.push_back({.driver_pin = pin,
                    .target_load_pin = worst_load,
                    .net = nullptr,
                    .criticality = criticality});
  }
  std::stable_sort(
      jobs.begin(), jobs.end(), [](const BufferJob& a, const BufferJob& b) {
        return a.criticality > b.criticality;
      });
  if (static_cast<int>(jobs.size()) > rs_config_.max_buffer_nets) {
    jobs.resize(rs_config_.max_buffer_nets);
  }
  Rebuffer& rebuffer = resizer_.rebuffer();
  if (rs_config_.serial_buffering) {
    rebuffer.init();
    rebuffer.initOnCorner(sta_->cmdScene());
    std::unordered_set<const sta::Net*> touched;
    int inserted = 0;
    int improved = 0;
    for (const BufferJob& job : jobs) {
      const sta::Net* net = network_->net(job.driver_pin);
      if (net == nullptr || touched.contains(net)) {
        continue;
      }
      // The driver may have lost its violation to an earlier tree.
      sta::Vertex* vertex = graph_->pinDrvrVertex(job.driver_pin);
      if (vertex == nullptr
          || sta::delayAsFloat(sta_->slack(vertex, max_)) >= 0.0f) {
        continue;
      }
      const int count = rebuffer.rebufferPin(job.driver_pin);
      if (count > 0) {
        touched.insert(net);
        inserted += count;
        ++improved;
      }
    }
    if (inserted > 0) {
      estimate_parasitics_->updateParasitics();
      sta_->findDelays();
    }
    logger_->info(RSZ,
                  3309,
                  "RS_SIZING buffering: {} candidate nets, {} rebuffered, {} "
                  "buffers inserted.",
                  jobs.size(),
                  improved,
                  inserted);
    seconds_buffering_ += secondsSince(start);
    return inserted;
  }
  // Main thread: nets as BufferedNets on the planner corner.
  sta::Scene* scene = sta_->cmdScene();
  rebuffer.initForPlanning(scene);
  if (planners_ == nullptr) {
    planners_ = std::make_unique<RebufferPlannerPool>(resizer_);
  }
  planners_->reset(scene, thread_pool_->threadCount());
  std::vector<BufferJob> ready;
  for (BufferJob& job : jobs) {
    job.net = rebuffer.makeBufferedNet(job.driver_pin);
    if (job.net != nullptr) {
      ready.push_back(std::move(job));
    }
  }
  // Workers: plan from the frozen required times.
  struct BufferChoice
  {
    bool valid{false};
    int plans{0};
    int plans_with_buffers{0};
    float best_gain{-std::numeric_limits<float>::infinity()};
    RebufferTreePlan plan;
  };
  RebufferPlannerPool& pool = *planners_;
  const std::vector<BufferChoice> choices = thread_pool_->parallelMap(
      ready, [this, &pool](const BufferJob& job) {
        BufferChoice choice;
        const RebufferPlannerPool::Lease lease(pool);
        std::vector<RebufferTreePlan> plans
            = lease.planner().makeRebufferTreePlans(
                job.driver_pin, job.net, job.target_load_pin, 4);
        choice.plans = static_cast<int>(plans.size());
        for (RebufferTreePlan& plan : plans) {
          if (plan.buffer_count > 0) {
            ++choice.plans_with_buffers;
            choice.best_gain = std::max(choice.best_gain, plan.driver_slack_gain);
          }
          if (plan.buffer_count == 0 || plan.driver_slack_gain <= 0.0f) {
            continue;
          }
          bool dont_touch_behind_buffer = false;
          for (const BufferTreeLoadPath& load_path : plan.load_paths) {
            for (const BufferTreePathStage& stage : load_path.stages) {
              if (stage.type == BufferTreeStageType::kBuffer
                  && resizer_.dontTouch(load_path.load_pin)) {
                dont_touch_behind_buffer = true;
              }
            }
          }
          if (dont_touch_behind_buffer) {
            continue;
          }
          if (!choice.valid
              || plan.driver_slack_gain > choice.plan.driver_slack_gain) {
            choice.plan = std::move(plan);
            choice.valid = true;
          }
        }
        return choice;
      });
  // Main thread: apply, one tree per net and none on a net that a tree
  // already changed.
  std::unordered_set<const sta::Net*> touched_nets;
  int inserted = 0;
  int plans = 0;
  int plans_with_buffers = 0;
  int valid = 0;
  float best_gain = -std::numeric_limits<float>::infinity();
  for (const BufferChoice& choice : choices) {
    plans += choice.plans;
    plans_with_buffers += choice.plans_with_buffers;
    valid += choice.valid ? 1 : 0;
    best_gain = std::max(best_gain, choice.best_gain);
  }
  for (size_t i = 0; i < choices.size(); ++i) {
    if (!choices[i].valid) {
      continue;
    }
    const sta::Net* net = network_->net(ready[i].driver_pin);
    if (net == nullptr || touched_nets.contains(net)) {
      continue;
    }
    const int count
        = rebuffer.applyRebufferTreePlan(choices[i].plan, "rs_buffer");
    if (count > 0) {
      touched_nets.insert(net);
      inserted += count;
    }
  }
  planners_->reset(scene, 0);
  if (inserted > 0) {
    estimate_parasitics_->updateParasitics();
    sta_->findDelays();
  }
  logger_->info(RSZ,
                3308,
                "RS_SIZING buffering: {} candidate nets, {} planned, {} plans "
                "({} with buffers, best gain {}), {} improving, {} buffers "
                "inserted.",
                jobs.size(),
                ready.size(),
                plans,
                plans_with_buffers,
                std::isfinite(best_gain)
                    ? sta::delayAsString(best_gain, 3, sta_)
                    : std::string("-"),
                valid,
                inserted);
  seconds_buffering_ += secondsSince(start);
  return inserted;
}

void ResourceSharingSizingPolicy::readConfig()
{
  rs_config_ = ResourceSharingSizingConfig{};
  const char* text = std::getenv("RSZ_RS_SIZING_OPTIONS");
  if (text == nullptr) {
    return;
  }
  std::istringstream stream(text);
  std::string item;
  while (stream >> item) {
    const size_t equal = item.find('=');
    if (equal == std::string::npos) {
      logger_->warn(RSZ, 3301, "RS_SIZING: ignoring option '{}'.", item);
      continue;
    }
    const std::string key = item.substr(0, equal);
    const std::string value = item.substr(equal + 1);
    const double number = std::atof(value.c_str());
    if (key == "multipliers") {
      if (value == "rs") {
        rs_config_.multipliers
            = ResourceSharingSizingConfig::Multipliers::kResourceSharing;
      } else if (value == "lr") {
        rs_config_.multipliers
            = ResourceSharingSizingConfig::Multipliers::kLagrangian;
      } else {
        logger_->warn(
            RSZ, 3306, "RS_SIZING: unknown multipliers '{}'.", value);
      }
    } else if (key == "lagrangian_exponent") {
      rs_config_.lagrangian_exponent = number;
    } else if (key == "lagrangian_floor") {
      rs_config_.lagrangian_floor = number;
    } else if (key == "lagrangian_power_weight") {
      rs_config_.lagrangian_power_weight = number;
    } else if (key == "iterations") {
      rs_config_.iterations = static_cast<int>(number);
    } else if (key == "rounds") {
      rs_config_.rounds_per_iteration = static_cast<int>(number);
    } else if (key == "gamma") {
      rs_config_.gamma = number;
    } else if (key == "accumulate") {
      rs_config_.accumulate_weights = number != 0.0;
    } else if (key == "power_budget_ratio") {
      rs_config_.power_budget_ratio = number;
    } else if (key == "leakage_weight") {
      rs_config_.leakage_weight = number;
    } else if (key == "min_relative_gain") {
      rs_config_.min_relative_gain = number;
    } else if (key == "max_siblings") {
      rs_config_.max_siblings = static_cast<int>(number);
    } else if (key == "size_sequential") {
      rs_config_.size_sequential = number != 0.0;
    } else if (key == "wns_loss_limit") {
      rs_config_.wns_loss_limit = number;
    } else if (key == "tns_loss_fraction") {
      rs_config_.tns_loss_fraction = number;
    } else if (key == "start_violation_ratio") {
      rs_config_.start_violation_ratio = number;
    } else if (key == "critical_weight_ratio") {
      rs_config_.critical_weight_ratio = number;
    } else if (key == "min_efficiency_ratio") {
      rs_config_.min_efficiency_ratio = number;
    } else if (key == "buffering") {
      rs_config_.buffering = number != 0.0;
    } else if (key == "buffer_weight_ratio") {
      rs_config_.buffer_weight_ratio = number;
    } else if (key == "max_buffer_nets") {
      rs_config_.max_buffer_nets = static_cast<int>(number);
    } else if (key == "serial_buffering") {
      rs_config_.serial_buffering = number != 0.0;
    } else if (key == "buffer_max_fanout") {
      rs_config_.buffer_max_fanout = static_cast<int>(number);
    } else if (key == "check_weights") {
      rs_config_.check_weights = number != 0.0;
    } else {
      logger_->warn(RSZ, 3302, "RS_SIZING: unknown option '{}'.", key);
    }
  }
}

bool ResourceSharingSizingPolicy::start()
{
  if (!OptimizationPolicy::start()) {
    return false;
  }
  readConfig();
  thread_pool_ = makeWorkerThreadPool();
  return true;
}

double ResourceSharingSizingPolicy::clockPeriod() const
{
  double period = 0.0;
  for (const sta::Clock* clock : sta_->cmdMode()->sdc()->clocks()) {
    const double clock_period = clock->period();
    if (clock_period > 0.0 && (period == 0.0 || clock_period < period)) {
      period = clock_period;
    }
  }
  return period;
}

const sta::Pvt* ResourceSharingSizingPolicy::pvtOf(
    const sta::Instance* instance) const
{
  const sta::Pvt* pvt = scene_->sdc()->pvt(instance, max_);
  if (pvt == nullptr) {
    pvt = scene_->sdc()->operatingConditions(max_);
  }
  return pvt;
}

double ResourceSharingSizingPolicy::powerOf(sta::LibertyCell* cell)
{
  auto found = power_cache_.find(cell);
  if (found != power_cache_.end()) {
    return found->second;
  }
  double area = 0.0;
  sta::dbNetwork* db_network = resizer_.dbNetwork();
  odb::dbMaster* master = db_network->staToDb(db_network->cell(cell));
  if (master != nullptr) {
    area = resizer_.dbuToMeters(master->getWidth())
           * resizer_.dbuToMeters(master->getHeight());
  }
  const double leakage = resizer_.cellLeakage(cell).value_or(0.0f);
  const double power = area * area_scale_
                       + rs_config_.leakage_weight * leakage * leakage_scale_;
  power_cache_[cell] = power;
  return power;
}

double ResourceSharingSizingPolicy::designPower()
{
  double power = 0.0;
  for (sta::Instance* instance : sizable_gates_) {
    power += powerOf(network_->libertyCell(instance));
  }
  return power;
}

void ResourceSharingSizingPolicy::collectSizableGates(const bool rescale)
{
  sizable_gates_.clear();
  sta_->ensureClkNetwork(sta_->cmdMode());
  double area_sum = 0.0;
  double leakage_sum = 0.0;
  sta::dbNetwork* db_network = resizer_.dbNetwork();
  std::unique_ptr<sta::LeafInstanceIterator> instances(
      network_->leafInstanceIterator());
  while (instances->hasNext()) {
    sta::Instance* instance = instances->next();
    sta::LibertyCell* cell = network_->libertyCell(instance);
    if (cell == nullptr || resizer_.dontTouch(instance)) {
      continue;
    }
    if (cell->isSequential() && !rs_config_.size_sequential) {
      continue;
    }
    bool clock_or_unusable = false;
    bool has_output = false;
    std::unique_ptr<sta::InstancePinIterator> pins(
        network_->pinIterator(instance));
    while (pins->hasNext()) {
      const sta::Pin* pin = pins->next();
      if (network_->direction(pin)->isAnyOutput()) {
        has_output = true;
        if (sta_->isClock(pin, sta_->cmdMode())) {
          clock_or_unusable = true;
        }
      }
    }
    if (!has_output || clock_or_unusable) {
      continue;
    }
    if (resizer_.getSwappableCells(cell).size() < 2) {
      continue;
    }
    odb::dbMaster* master = db_network->staToDb(db_network->cell(cell));
    if (master == nullptr) {
      continue;
    }
    area_sum += resizer_.dbuToMeters(master->getWidth())
                * resizer_.dbuToMeters(master->getHeight());
    leakage_sum += resizer_.cellLeakage(cell).value_or(0.0f);
    sizable_gates_.push_back(instance);
  }
  // Area and leakage each contribute one unit for the initial design.
  if (rescale) {
    area_scale_ = area_sum > 0.0 ? 1.0 / area_sum : 0.0;
    leakage_scale_ = leakage_sum > 0.0 ? 1.0 / leakage_sum : 0.0;
    power_cache_.clear();
  }
}

double ResourceSharingSizingPolicy::arcLogWeight(sta::Vertex* from,
                                                 const int from_rf,
                                                 sta::Vertex* to,
                                                 const int to_rf,
                                                 const float delay) const
{
  const std::vector<double>& forward
      = rs_config_.accumulate_weights ? log_forward_mean_ : log_forward_;
  const std::vector<double>& backward
      = rs_config_.accumulate_weights ? log_backward_mean_ : log_backward_;
  const size_t from_index = 2 * graph_->id(from) + from_rf;
  const size_t to_index = 2 * graph_->id(to) + to_rf;
  if (from_index >= forward.size() || to_index >= backward.size()) {
    return kNegativeInfinity;
  }
  const double log_from = forward[from_index];
  const double log_to = backward[to_index];
  if (log_from == kNegativeInfinity || log_to == kNegativeInfinity) {
    return kNegativeInfinity;
  }
  return log_from + rs_config_.gamma * delay / period_ + log_to;
}

void ResourceSharingSizingPolicy::buildLevelOrder()
{
  sta_->ensureLevelized();
  level_order_.clear();
  sta::VertexId max_id = 0;
  {
    sta::VertexIterator vertices(graph_);
    while (vertices.hasNext()) {
      sta::Vertex* vertex = vertices.next();
      level_order_.push_back(vertex);
      max_id = std::max(max_id, graph_->id(vertex));
    }
  }
  std::stable_sort(level_order_.begin(),
                   level_order_.end(),
                   [](const sta::Vertex* a, const sta::Vertex* b) {
                     return a->level() < b->level();
                   });
}

void ResourceSharingSizingPolicy::computeLogWeights()
{
  sta::VertexId max_id = 0;
  for (sta::Vertex* vertex : level_order_) {
    max_id = std::max(max_id, graph_->id(vertex));
  }
  const size_t size = 2 * (static_cast<size_t>(max_id) + 1);
  log_forward_.assign(size, kNegativeInfinity);
  log_backward_.assign(size, kNegativeInfinity);
  const double scale = rs_config_.gamma / period_;
  const sta::SceneSeq scenes{const_cast<sta::Scene*>(scene_)};
  sta::Search* search = sta_->search();

  // Forward: starts take their STA arrival; other vertices sum their fanin.
  for (sta::Vertex* vertex : level_order_) {
    const size_t base = 2 * graph_->id(vertex);
    bool has_fanin = false;
    sta::VertexInEdgeIterator edges(vertex, graph_);
    while (edges.hasNext()) {
      sta::Edge* edge = edges.next();
      if (!isPropagatedArcSet(edge)) {
        continue;
      }
      has_fanin = true;
      sta::Vertex* from = edge->from(graph_);
      const size_t from_base = 2 * graph_->id(from);
      for (const sta::TimingArc* arc : edge->timingArcSet()->arcs()) {
        const int from_rf = riseFallIndex(arc->fromEdge());
        const int to_rf = riseFallIndex(arc->toEdge());
        const double log_from = log_forward_[from_base + from_rf];
        if (log_from == kNegativeInfinity) {
          continue;
        }
        const float delay = sta::delayAsFloat(
            graph_->arcDelay(edge, const_cast<sta::TimingArc*>(arc), dcalc_ap_));
        log_forward_[base + to_rf]
            = logSumExp(log_forward_[base + to_rf], log_from + scale * delay);
      }
    }
    if (!has_fanin) {
      for (const sta::RiseFall* rf : sta::RiseFall::range()) {
        const sta::Arrival arrival = sta_->arrival(
            vertex, rf->asRiseFallBoth(), scenes, max_);
        const float value = sta::delayAsFloat(arrival);
        if (std::abs(value) < 1e10f) {
          log_forward_[base + rf->index()] = scale * value;
        }
      }
    }
  }

  // Backward: ends take minus their STA required time; every vertex also sums
  // its fanout.
  for (auto it = level_order_.rbegin(); it != level_order_.rend(); ++it) {
    sta::Vertex* vertex = *it;
    const size_t base = 2 * graph_->id(vertex);
    if (search->isEndpoint(vertex)) {
      for (const sta::RiseFall* rf : sta::RiseFall::range()) {
        const sta::Required required = sta_->required(
            vertex, rf->asRiseFallBoth(), scenes, max_);
        const float value = sta::delayAsFloat(required);
        if (std::abs(value) < 1e10f) {
          log_backward_[base + rf->index()] = -scale * value;
        }
      }
    }
    sta::VertexOutEdgeIterator edges(vertex, graph_);
    while (edges.hasNext()) {
      sta::Edge* edge = edges.next();
      if (!isPropagatedArcSet(edge)) {
        continue;
      }
      sta::Vertex* to = edge->to(graph_);
      const size_t to_base = 2 * graph_->id(to);
      for (const sta::TimingArc* arc : edge->timingArcSet()->arcs()) {
        const int from_rf = riseFallIndex(arc->fromEdge());
        const int to_rf = riseFallIndex(arc->toEdge());
        const double log_to = log_backward_[to_base + to_rf];
        if (log_to == kNegativeInfinity) {
          continue;
        }
        const float delay = sta::delayAsFloat(
            graph_->arcDelay(edge, const_cast<sta::TimingArc*>(arc), dcalc_ap_));
        log_backward_[base + from_rf]
            = logSumExp(log_backward_[base + from_rf], log_to + scale * delay);
      }
    }
  }
}

void ResourceSharingSizingPolicy::accumulateLogWeights(const int iteration)
{
  // Mean of the per-iteration log weights.  Vertices without a value in one
  // iteration take the current value.
  if (iteration == 0 || log_forward_mean_.size() != log_forward_.size()) {
    log_forward_mean_ = log_forward_;
    log_backward_mean_ = log_backward_;
    return;
  }
  const double old_share = static_cast<double>(iteration) / (iteration + 1);
  const double new_share = 1.0 / (iteration + 1);
  auto mix = [&](std::vector<double>& mean, const std::vector<double>& now) {
    for (size_t i = 0; i < mean.size(); ++i) {
      if (now[i] == kNegativeInfinity) {
        mean[i] = kNegativeInfinity;
      } else if (mean[i] == kNegativeInfinity) {
        mean[i] = now[i];
      } else {
        mean[i] = old_share * mean[i] + new_share * now[i];
      }
    }
  };
  mix(log_forward_mean_, log_forward_);
  mix(log_backward_mean_, log_backward_);
}

void ResourceSharingSizingPolicy::updateLagrangeMultipliers(
    const int iteration)
{
  const sta::SceneSeq scenes{const_cast<sta::Scene*>(scene_)};
  sta::Search* search = sta_->search();
  sta::VertexId max_vertex_id = 0;
  sta::EdgeId max_edge_id = 0;
  for (sta::Vertex* vertex : level_order_) {
    max_vertex_id = std::max(max_vertex_id, graph_->id(vertex));
    sta::VertexOutEdgeIterator edges(vertex, graph_);
    while (edges.hasNext()) {
      max_edge_id = std::max(max_edge_id, graph_->id(edges.next()));
    }
  }
  // Arrival and required time per vertex and rise/fall.
  const size_t vertex_slots = 2 * (static_cast<size_t>(max_vertex_id) + 1);
  std::vector<float> arrival(vertex_slots, -kFloatInfinity);
  std::vector<float> required(vertex_slots, kFloatInfinity);
  for (sta::Vertex* vertex : level_order_) {
    const size_t base = 2 * graph_->id(vertex);
    for (const sta::RiseFall* rf : sta::RiseFall::range()) {
      const float a = sta::delayAsFloat(
          sta_->arrival(vertex, rf->asRiseFallBoth(), scenes, max_));
      const float r = sta::delayAsFloat(
          sta_->required(vertex, rf->asRiseFallBoth(), scenes, max_));
      if (std::abs(a) < 1e10f) {
        arrival[base + rf->index()] = a;
      }
      if (std::abs(r) < 1e10f) {
        required[base + rf->index()] = r;
      }
    }
  }
  const size_t arc_slots = 4 * (static_cast<size_t>(max_edge_id) + 1);
  const bool first = iteration == 0 || arc_multipliers_.size() != arc_slots
                     || endpoint_multipliers_.size() != vertex_slots;
  if (first) {
    arc_multipliers_.assign(arc_slots, 1.0);
    endpoint_multipliers_.assign(vertex_slots, 1.0);
  }
  auto factor = [this](const double criticality) {
    return std::pow(
        std::max(rs_config_.lagrangian_floor, 1.0 + criticality),
        rs_config_.lagrangian_exponent);
  };
  // Multiplicative update by local criticality.  Arcs without a constrained
  // path through them get multiplier zero.
  for (sta::Vertex* vertex : level_order_) {
    const size_t from_base = 2 * graph_->id(vertex);
    sta::VertexOutEdgeIterator edges(vertex, graph_);
    while (edges.hasNext()) {
      sta::Edge* edge = edges.next();
      const size_t edge_base = 4 * graph_->id(edge);
      if (!isPropagatedArcSet(edge)) {
        continue;
      }
      const size_t to_base = 2 * graph_->id(edge->to(graph_));
      for (const sta::TimingArc* arc : edge->timingArcSet()->arcs()) {
        const int from_rf = riseFallIndex(arc->fromEdge());
        const int to_rf = riseFallIndex(arc->toEdge());
        const float a = arrival[from_base + from_rf];
        const float r = required[to_base + to_rf];
        double& multiplier = arc_multipliers_[edge_base + 2 * from_rf + to_rf];
        if (a == -kFloatInfinity || r == kFloatInfinity) {
          multiplier = 0.0;
          continue;
        }
        const float delay = sta::delayAsFloat(
            graph_->arcDelay(edge, const_cast<sta::TimingArc*>(arc), dcalc_ap_));
        const double criticality = (a + delay - r) / period_;
        if (multiplier == 0.0) {
          multiplier = 1.0;
        }
        multiplier *= factor(criticality);
      }
    }
  }
  // Endpoint multipliers by endpoint slack.
  for (sta::Vertex* vertex : level_order_) {
    if (!search->isEndpoint(vertex)) {
      continue;
    }
    const size_t base = 2 * graph_->id(vertex);
    for (int rf = 0; rf < 2; ++rf) {
      if (arrival[base + rf] == -kFloatInfinity
          || required[base + rf] == kFloatInfinity) {
        endpoint_multipliers_[base + rf] = 0.0;
        continue;
      }
      double& multiplier = endpoint_multipliers_[base + rf];
      if (multiplier == 0.0) {
        multiplier = 1.0;
      }
      multiplier
          *= factor((arrival[base + rf] - required[base + rf]) / period_);
    }
  }
  // Flow projection in reverse level order: the multipliers into a vertex
  // transition are scaled to sum to the multipliers out of it (plus its
  // endpoint multiplier).  Out-arcs of a vertex are final when it is visited
  // because their heads have higher levels.
  for (auto it = level_order_.rbegin(); it != level_order_.rend(); ++it) {
    sta::Vertex* vertex = *it;
    const size_t base = 2 * graph_->id(vertex);
    std::array<double, 2> outflow{0.0, 0.0};
    if (search->isEndpoint(vertex)) {
      outflow[0] = endpoint_multipliers_[base];
      outflow[1] = endpoint_multipliers_[base + 1];
    }
    sta::VertexOutEdgeIterator out_edges(vertex, graph_);
    while (out_edges.hasNext()) {
      sta::Edge* edge = out_edges.next();
      if (!isPropagatedArcSet(edge)) {
        continue;
      }
      const size_t edge_base = 4 * graph_->id(edge);
      for (const sta::TimingArc* arc : edge->timingArcSet()->arcs()) {
        const int from_rf = riseFallIndex(arc->fromEdge());
        const int to_rf = riseFallIndex(arc->toEdge());
        outflow[from_rf] += arc_multipliers_[edge_base + 2 * from_rf + to_rf];
      }
    }
    std::array<double, 2> inflow{0.0, 0.0};
    sta::VertexInEdgeIterator in_edges(vertex, graph_);
    while (in_edges.hasNext()) {
      sta::Edge* edge = in_edges.next();
      if (!isPropagatedArcSet(edge)) {
        continue;
      }
      const size_t edge_base = 4 * graph_->id(edge);
      for (const sta::TimingArc* arc : edge->timingArcSet()->arcs()) {
        const int from_rf = riseFallIndex(arc->fromEdge());
        const int to_rf = riseFallIndex(arc->toEdge());
        inflow[to_rf] += arc_multipliers_[edge_base + 2 * from_rf + to_rf];
      }
    }
    std::array<double, 2> scale{0.0, 0.0};
    for (int rf = 0; rf < 2; ++rf) {
      scale[rf] = inflow[rf] > 0.0 ? outflow[rf] / inflow[rf] : 0.0;
    }
    sta::VertexInEdgeIterator scale_edges(vertex, graph_);
    while (scale_edges.hasNext()) {
      sta::Edge* edge = scale_edges.next();
      if (!isPropagatedArcSet(edge)) {
        continue;
      }
      const size_t edge_base = 4 * graph_->id(edge);
      for (const sta::TimingArc* arc : edge->timingArcSet()->arcs()) {
        const int from_rf = riseFallIndex(arc->fromEdge());
        const int to_rf = riseFallIndex(arc->toEdge());
        arc_multipliers_[edge_base + 2 * from_rf + to_rf] *= scale[to_rf];
      }
    }
  }
}

void ResourceSharingSizingPolicy::normalizeLagrangeMultipliers()
{
  // The clock tree carries the flow of every path; normalize by the largest
  // multiplier on data arcs so data-path criticalities spread over [0, 1].
  double largest = 0.0;
  for (sta::Vertex* vertex : level_order_) {
    if (sta_->isClock(vertex->pin(), sta_->cmdMode())) {
      continue;
    }
    sta::VertexOutEdgeIterator edges(vertex, graph_);
    while (edges.hasNext()) {
      const size_t edge_base = 4 * static_cast<size_t>(graph_->id(edges.next()));
      for (size_t slot = 0; slot < 4; ++slot) {
        if (edge_base + slot < arc_multipliers_.size()) {
          largest = std::max(largest, arc_multipliers_[edge_base + slot]);
        }
      }
    }
  }
  if (largest <= 0.0) {
    return;
  }
  for (double& multiplier : arc_multipliers_) {
    multiplier /= largest;
  }
  for (double& multiplier : endpoint_multipliers_) {
    multiplier /= largest;
  }
}

double ResourceSharingSizingPolicy::gateCriticality(
    sta::Instance* instance) const
{
  double criticality = 0.0;
  std::unique_ptr<sta::InstancePinIterator> pins(
      network_->pinIterator(instance));
  while (pins->hasNext()) {
    const sta::Pin* pin = pins->next();
    if (!network_->direction(pin)->isAnyOutput()) {
      continue;
    }
    sta::Vertex* vertex = graph_->pinDrvrVertex(pin);
    if (vertex == nullptr) {
      continue;
    }
    sta::VertexInEdgeIterator edges(vertex, graph_);
    while (edges.hasNext()) {
      sta::Edge* edge = edges.next();
      if (edge->isWire() || !isPropagatedArcSet(edge)) {
        continue;
      }
      for (const sta::TimingArc* arc : edge->timingArcSet()->arcs()) {
        const float delay = sta::delayAsFloat(graph_->arcDelay(
            edge, const_cast<sta::TimingArc*>(arc), dcalc_ap_));
        criticality = std::max(
            criticality,
            arcWeight(edge, arc, edge->from(graph_), vertex, delay)
                * period_);
      }
    }
  }
  return criticality;
}

int ResourceSharingSizingPolicy::violatingEndpointCount()
{
  int count = 0;
  for (sta::Vertex* vertex : sta_->endpoints()) {
    if (sta::delayAsFloat(sta_->slack(vertex, max_)) < 0.0f) {
      ++count;
    }
  }
  return count;
}

double ResourceSharingSizingPolicy::arcWeight(sta::Edge* edge,
                                              const sta::TimingArc* arc,
                                              sta::Vertex* from,
                                              sta::Vertex* to,
                                              const float delay) const
{
  const int from_rf = riseFallIndex(arc->fromEdge());
  const int to_rf = riseFallIndex(arc->toEdge());
  if (rs_config_.multipliers
      == ResourceSharingSizingConfig::Multipliers::kLagrangian) {
    const size_t slot = 4 * static_cast<size_t>(graph_->id(edge))
                        + 2 * from_rf + to_rf;
    return slot < arc_multipliers_.size() ? arc_multipliers_[slot] / period_
                                          : 0.0;
  }
  const double log_weight = arcLogWeight(from, from_rf, to, to_rf, delay);
  if (log_weight == kNegativeInfinity) {
    return 0.0;
  }
  return std::exp(log_weight - log_weight_shift_) / period_;
}

double ResourceSharingSizingPolicy::powerWeight() const
{
  if (rs_config_.multipliers
      == ResourceSharingSizingConfig::Multipliers::kLagrangian) {
    return rs_config_.lagrangian_power_weight;
  }
  return std::exp(log_power_weight_ - log_weight_shift_);
}

void ResourceSharingSizingPolicy::checkWeights()
{
  // With gamma -> infinity, (forward + backward) * period / gamma tends to
  // minus the slack.  Print the comparison on the worst vertices.
  const double scale = rs_config_.gamma / period_;
  int printed = 0;
  for (sta::Vertex* vertex : level_order_) {
    const size_t base = 2 * graph_->id(vertex);
    const double soft
        = -std::max(log_forward_[base] + log_backward_[base],
                    log_forward_[base + 1] + log_backward_[base + 1])
          / scale;
    const float slack = sta::delayAsFloat(sta_->slack(vertex, max_));
    if (slack < 0.0f && printed < 20 && (graph_->id(vertex) % 97) == 0) {
      logger_->report("RS_SIZING check {} sta_slack {:.4g} soft_slack {:.4g}",
                      network_->pathName(vertex->pin()),
                      slack,
                      soft);
      ++printed;
    }
  }
}

bool ResourceSharingSizingPolicy::fillTableArc(sta::Edge* edge,
                                               const sta::TimingArc* arc,
                                               sta::Vertex* from,
                                               sta::Vertex* to,
                                               const float load_cap,
                                               TableArc& table_arc) const
{
  table_arc.arc = arc;
  table_arc.input_rf = riseFallIndex(arc->fromEdge());
  table_arc.output_rf = riseFallIndex(arc->toEdge());
  table_arc.pvt = pvtOf(network_->instance(from->pin()));
  table_arc.input_slew = sta::delayAsFloat(graph_->slew(
      from, sta::RiseFall::find(table_arc.input_rf), dcalc_ap_));
  table_arc.load_cap = load_cap;
  const float delay = sta::delayAsFloat(
      graph_->arcDelay(edge, const_cast<sta::TimingArc*>(arc), dcalc_ap_));
  table_arc.weight = arcWeight(edge, arc, from, to, delay);
  if (!tableDelay(arc,
                  table_arc.pvt,
                  table_arc.input_slew,
                  table_arc.load_cap,
                  table_arc.delay,
                  table_arc.output_slew)) {
    return false;
  }
  return true;
}

bool ResourceSharingSizingPolicy::snapshot(sta::Instance* instance,
                                           GateSnapshot& snap)
{
  sta::LibertyCell* cell = network_->libertyCell(instance);
  if (cell == nullptr) {
    return false;
  }
  snap.instance = instance;
  snap.current_cell = cell;
  snap.current_power = powerOf(cell);
  sta::GraphDelayCalc* graph_delay_calc = sta_->graphDelayCalc();

  // Pins of the gate.
  std::vector<const sta::Pin*> input_pins;
  std::vector<const sta::Pin*> output_pins;
  {
    std::unique_ptr<sta::InstancePinIterator> pins(
        network_->pinIterator(instance));
    while (pins->hasNext()) {
      const sta::Pin* pin = pins->next();
      if (network_->libertyPort(pin) == nullptr) {
        continue;
      }
      const sta::PortDirection* direction = network_->direction(pin);
      if (direction->isAnyOutput()) {
        output_pins.push_back(pin);
      } else if (direction->isAnyInput()) {
        input_pins.push_back(pin);
      }
    }
  }
  if (output_pins.empty()) {
    return false;
  }

  // Inputs: net driver arcs and sibling load arcs.
  for (const sta::Pin* pin : input_pins) {
    InputPinContext input;
    input.port = network_->libertyPort(pin);
    input.current_cap = input.port->capacitance(max_);
    sta::Vertex* load_vertex = graph_->pinLoadVertex(pin);
    if (load_vertex != nullptr) {
      sta::VertexInEdgeIterator wires(load_vertex, graph_);
      while (wires.hasNext()) {
        sta::Edge* wire = wires.next();
        if (!wire->isWire()) {
          continue;
        }
        sta::Vertex* driver_vertex = wire->from(graph_);
        const sta::Pin* driver_pin = driver_vertex->pin();
        const sta::Instance* driver = network_->instance(driver_pin);
        if (network_->isTopInstance(driver)
            || network_->libertyCell(driver) == nullptr
            || sta_->isClock(driver_pin, sta_->cmdMode())) {
          continue;
        }
        const float driver_load
            = graph_delay_calc->loadCap(driver_pin, scene_, max_);
        sta::VertexInEdgeIterator driver_edges(driver_vertex, graph_);
        while (driver_edges.hasNext()) {
          sta::Edge* edge = driver_edges.next();
          if (edge->isWire() || !isPropagatedArcSet(edge)) {
            continue;
          }
          for (const sta::TimingArc* arc : edge->timingArcSet()->arcs()) {
            TableArc table_arc;
            if (fillTableArc(edge,
                             arc,
                             edge->from(graph_),
                             driver_vertex,
                             driver_load,
                             table_arc)
                && table_arc.weight > 0.0) {
              input.driver_arcs.push_back(table_arc);
            }
          }
        }
        // Siblings: the other loads of the driver.
        int siblings = 0;
        sta::VertexOutEdgeIterator loads(driver_vertex, graph_);
        while (loads.hasNext() && siblings < rs_config_.max_siblings) {
          sta::Edge* load_wire = loads.next();
          if (!load_wire->isWire()) {
            continue;
          }
          sta::Vertex* sibling_vertex = load_wire->to(graph_);
          const sta::Pin* sibling_pin = sibling_vertex->pin();
          const sta::Instance* sibling = network_->instance(sibling_pin);
          if (sibling == instance || network_->isTopInstance(sibling)
              || network_->libertyCell(sibling) == nullptr) {
            continue;
          }
          ++siblings;
          sta::VertexOutEdgeIterator sibling_edges(sibling_vertex, graph_);
          while (sibling_edges.hasNext()) {
            sta::Edge* edge = sibling_edges.next();
            if (edge->isWire() || !isPropagatedArcSet(edge)) {
              continue;
            }
            sta::Vertex* sibling_out = edge->to(graph_);
            const float sibling_load = graph_delay_calc->loadCap(
                sibling_out->pin(), scene_, max_);
            for (const sta::TimingArc* arc : edge->timingArcSet()->arcs()) {
              TableArc table_arc;
              if (fillTableArc(edge,
                               arc,
                               sibling_vertex,
                               sibling_out,
                               sibling_load,
                               table_arc)
                  && table_arc.weight > 0.0) {
                input.sibling_arcs.push_back(table_arc);
              }
            }
          }
        }
      }
    }
    snap.inputs.push_back(std::move(input));
  }

  // Outputs: own arcs and fanout arcs.
  for (size_t output_index = 0; output_index < output_pins.size();
       ++output_index) {
    const sta::Pin* pin = output_pins[output_index];
    OutputPinContext output;
    output.port = network_->libertyPort(pin);
    output.load_cap = graph_delay_calc->loadCap(pin, scene_, max_);
    output.slew_limit = kFloatInfinity;
    sta::Vertex* driver_vertex = graph_->pinDrvrVertex(pin);
    if (driver_vertex == nullptr) {
      return false;
    }
    float driver_slew = 0.0f;
    for (const sta::RiseFall* rf : sta::RiseFall::range()) {
      driver_slew = std::max(
          driver_slew,
          sta::delayAsFloat(graph_->slew(driver_vertex, rf, dcalc_ap_)));
    }
    {
      float limit = 0.0f;
      bool exists = false;
      sta_->findSlewLimit(output.port, scene_, max_, limit, exists);
      if (exists && limit > 0.0f) {
        output.slew_limit = std::min(output.slew_limit, limit);
      }
    }
    // Own arcs into this output.
    sta::VertexInEdgeIterator own_edges(driver_vertex, graph_);
    while (own_edges.hasNext()) {
      sta::Edge* edge = own_edges.next();
      if (edge->isWire() || !isPropagatedArcSet(edge)) {
        continue;
      }
      sta::Vertex* from = edge->from(graph_);
      const sta::Pin* from_pin = from->pin();
      int input_index = -1;
      for (size_t i = 0; i < input_pins.size(); ++i) {
        if (input_pins[i] == from_pin) {
          input_index = static_cast<int>(i);
        }
      }
      for (const sta::TimingArc* arc : edge->timingArcSet()->arcs()) {
        GateArc gate_arc;
        gate_arc.input_index = input_index;
        gate_arc.output_index = static_cast<int>(output_index);
        if (!fillTableArc(
                edge, arc, from, driver_vertex, output.load_cap,
                gate_arc.table_arc)) {
          return false;
        }
        output.current_slew[gate_arc.table_arc.output_rf]
            = std::max(output.current_slew[gate_arc.table_arc.output_rf],
                       gate_arc.table_arc.output_slew);
        snap.weighted_delay
            += gate_arc.table_arc.weight * gate_arc.table_arc.delay;
        snap.gate_arcs.push_back(gate_arc);
      }
    }
    // Fanout arcs.
    sta::VertexOutEdgeIterator loads(driver_vertex, graph_);
    while (loads.hasNext()) {
      sta::Edge* wire = loads.next();
      if (!wire->isWire()) {
        continue;
      }
      sta::Vertex* load_vertex = wire->to(graph_);
      const sta::Pin* load_pin = load_vertex->pin();
      for (const sta::RiseFall* rf : sta::RiseFall::range()) {
        output.wire_slew_degradation = std::max(
            output.wire_slew_degradation,
            sta::delayAsFloat(graph_->slew(load_vertex, rf, dcalc_ap_))
                - driver_slew);
      }
      const sta::LibertyPort* load_port = network_->libertyPort(load_pin);
      if (load_port != nullptr) {
        float limit = 0.0f;
        bool exists = false;
        sta_->findSlewLimit(load_port, scene_, max_, limit, exists);
        if (exists && limit > 0.0f) {
          output.slew_limit = std::min(output.slew_limit, limit);
        }
      }
      const sta::Instance* load = network_->instance(load_pin);
      if (network_->isTopInstance(load)
          || network_->libertyCell(load) == nullptr) {
        continue;
      }
      sta::VertexOutEdgeIterator load_edges(load_vertex, graph_);
      while (load_edges.hasNext()) {
        sta::Edge* edge = load_edges.next();
        if (edge->isWire() || !isPropagatedArcSet(edge)) {
          continue;
        }
        sta::Vertex* load_out = edge->to(graph_);
        const float fanout_load
            = graph_delay_calc->loadCap(load_out->pin(), scene_, max_);
        for (const sta::TimingArc* arc : edge->timingArcSet()->arcs()) {
          TableArc table_arc;
          if (fillTableArc(
                  edge, arc, load_vertex, load_out, fanout_load, table_arc)
              && table_arc.weight > 0.0) {
            output.fanout_arcs.push_back(table_arc);
          }
        }
      }
    }
    snap.outputs.push_back(std::move(output));
  }
  if (snap.gate_arcs.empty()) {
    return false;
  }

  // Candidates with matching arcs and pins.
  const std::vector<Resizer::InputCapBudget> cap_budgets
      = resizer_.inputCapBudgets(instance);
  for (sta::LibertyCell* candidate_cell : resizer_.getSwappableCells(cell)) {
    if (candidate_cell == cell) {
      continue;
    }
    if (!resizer_.replacementFitsInputCapBudgets(cap_budgets, candidate_cell)) {
      continue;
    }
    sta::LibertyCell* scene_cell = candidate_cell->sceneCell(scene_, max_);
    if (scene_cell == nullptr) {
      continue;
    }
    CandidateCell candidate;
    candidate.cell = candidate_cell;
    candidate.power = powerOf(candidate_cell);
    bool usable = true;
    for (const InputPinContext& input : snap.inputs) {
      const sta::LibertyPort* port
          = scene_cell->findLibertyPort(input.port->name());
      if (port == nullptr) {
        usable = false;
        break;
      }
      candidate.input_caps.push_back(port->capacitance(max_));
    }
    for (const OutputPinContext& output : snap.outputs) {
      sta::LibertyPort* port = scene_cell->findLibertyPort(output.port->name());
      if (port == nullptr) {
        usable = false;
        break;
      }
      float cap_limit = 0.0f;
      bool exists = false;
      port->capacitanceLimit(max_, cap_limit, exists);
      candidate.output_cap_limits.push_back(
          exists && cap_limit > 0.0f ? cap_limit : kFloatInfinity);
    }
    for (const GateArc& gate_arc : snap.gate_arcs) {
      if (!usable) {
        break;
      }
      const sta::TimingArc* ref = gate_arc.table_arc.arc;
      const sta::LibertyPort* from
          = scene_cell->findLibertyPort(ref->from()->name());
      const sta::LibertyPort* to
          = scene_cell->findLibertyPort(ref->to()->name());
      const sta::TimingArc* match = nullptr;
      if (from != nullptr && to != nullptr) {
        for (const sta::TimingArcSet* arc_set :
             scene_cell->timingArcSets(from, to)) {
          if (arc_set->role()->isTimingCheck()) {
            continue;
          }
          match = findMatchingTimingArc(ref, arc_set);
          if (match != nullptr) {
            break;
          }
        }
      }
      if (match == nullptr) {
        usable = false;
      }
      candidate.gate_arcs.push_back(match);
    }
    if (usable) {
      snap.candidates.push_back(std::move(candidate));
    }
  }
  return !snap.candidates.empty();
}

double ResourceSharingSizingPolicy::objective(const GateSnapshot& snap,
                                              const CandidateCell& candidate,
                                              bool& feasible) const
{
  feasible = true;
  double timing = 0.0;
  // Driver slew change per input and transition.
  std::vector<std::array<float, 2>> driver_slew_change(
      snap.inputs.size(), std::array<float, 2>{0.0f, 0.0f});
  for (size_t i = 0; i < snap.inputs.size(); ++i) {
    const InputPinContext& input = snap.inputs[i];
    const float cap_change = candidate.input_caps[i] - input.current_cap;
    if (cap_change == 0.0f) {
      continue;
    }
    for (const TableArc& arc : input.driver_arcs) {
      float delay = 0.0f;
      float slew = 0.0f;
      if (!tableDelay(arc.arc,
                      arc.pvt,
                      arc.input_slew,
                      std::max(0.0f, arc.load_cap + cap_change),
                      delay,
                      slew)) {
        continue;
      }
      timing += arc.weight * (delay - arc.delay);
      float& change = driver_slew_change[i][arc.output_rf];
      const float slew_change = slew - arc.output_slew;
      if (std::abs(slew_change) > std::abs(change)) {
        change = slew_change;
      }
    }
    for (const TableArc& arc : input.sibling_arcs) {
      const float slew_change = driver_slew_change[i][arc.input_rf];
      if (slew_change == 0.0f) {
        continue;
      }
      float delay = 0.0f;
      float slew = 0.0f;
      if (tableDelay(arc.arc,
                     arc.pvt,
                     std::max(0.0f, arc.input_slew + slew_change),
                     arc.load_cap,
                     delay,
                     slew)) {
        timing += arc.weight * (delay - arc.delay);
      }
    }
  }
  // Own arcs.
  std::vector<std::array<float, 2>> output_slew(
      snap.outputs.size(), std::array<float, 2>{0.0f, 0.0f});
  for (size_t a = 0; a < snap.gate_arcs.size(); ++a) {
    const GateArc& gate_arc = snap.gate_arcs[a];
    const TableArc& arc = gate_arc.table_arc;
    float input_slew = arc.input_slew;
    if (gate_arc.input_index >= 0) {
      input_slew += driver_slew_change[gate_arc.input_index][arc.input_rf];
    }
    float delay = 0.0f;
    float slew = 0.0f;
    if (!tableDelay(candidate.gate_arcs[a],
                    arc.pvt,
                    std::max(0.0f, input_slew),
                    arc.load_cap,
                    delay,
                    slew)) {
      feasible = false;
      return 0.0;
    }
    timing += arc.weight * (delay - arc.delay);
    float& out = output_slew[gate_arc.output_index][arc.output_rf];
    out = std::max(out, slew);
  }
  // Output limits and fanout arcs.
  for (size_t o = 0; o < snap.outputs.size(); ++o) {
    const OutputPinContext& output = snap.outputs[o];
    if (output.load_cap > candidate.output_cap_limits[o]) {
      feasible = false;
      return 0.0;
    }
    const float new_slew = std::max(output_slew[o][0], output_slew[o][1]);
    const float old_slew
        = std::max(output.current_slew[0], output.current_slew[1]);
    if (new_slew + output.wire_slew_degradation > output.slew_limit
        && new_slew > old_slew) {
      feasible = false;
      return 0.0;
    }
    for (const TableArc& arc : output.fanout_arcs) {
      const float slew_change
          = output_slew[o][arc.input_rf] - output.current_slew[arc.input_rf];
      if (slew_change == 0.0f) {
        continue;
      }
      float delay = 0.0f;
      float slew = 0.0f;
      if (tableDelay(arc.arc,
                     arc.pvt,
                     std::max(0.0f, arc.input_slew + slew_change),
                     arc.load_cap,
                     delay,
                     slew)) {
        timing += arc.weight * (delay - arc.delay);
      }
    }
  }
  const double power_term = powerWeight() * (candidate.power - snap.current_power)
                            / power_budget_;
  return timing + power_term;
}

ResourceSharingSizingPolicy::GateDecision
ResourceSharingSizingPolicy::evaluate(const GateSnapshot& snap) const
{
  GateDecision decision;
  decision.instance = snap.instance;
  const double power_weight = powerWeight();
  // Gain threshold: a fraction of the larger of the gate's weighted delay and
  // its weighted power.
  const double scale
      = std::max(snap.weighted_delay,
                 power_weight * snap.current_power / power_budget_);
  const double threshold = rs_config_.min_relative_gain * scale;
  double best = -threshold;
  for (const CandidateCell& candidate : snap.candidates) {
    bool feasible = false;
    const double change = objective(snap, candidate, feasible);
    if (feasible && change < best) {
      best = change;
      decision.cell = candidate.cell;
      decision.gain = -change;
    }
  }
  return decision;
}

int ResourceSharingSizingPolicy::runRound(std::vector<sta::Instance*>& gates)
{
  auto start = std::chrono::steady_clock::now();
  std::vector<GateSnapshot> snapshots;
  snapshots.reserve(gates.size());
  for (sta::Instance* instance : gates) {
    GateSnapshot snap;
    if (snapshot(instance, snap)) {
      snapshots.push_back(std::move(snap));
    }
  }
  seconds_snapshots_ += secondsSince(start);
  start = std::chrono::steady_clock::now();
  const std::vector<GateDecision> decisions = thread_pool_->parallelMap(
      snapshots, [this](const GateSnapshot& snap) { return evaluate(snap); });
  seconds_oracle_ += secondsSince(start);
  start = std::chrono::steady_clock::now();

  std::vector<const GateDecision*> wanted;
  for (const GateDecision& decision : decisions) {
    if (decision.cell != nullptr) {
      wanted.push_back(&decision);
    }
  }
  std::stable_sort(wanted.begin(),
                   wanted.end(),
                   [](const GateDecision* a, const GateDecision* b) {
                     return a->gain > b->gain;
                   });
  // Resize at most one gate per net in a round.
  std::unordered_set<const sta::Net*> touched_nets;
  std::vector<sta::Instance*> deferred;
  int resized = 0;
  for (const GateDecision* decision : wanted) {
    std::vector<const sta::Net*> nets;
    bool conflict = false;
    std::unique_ptr<sta::InstancePinIterator> pins(
        network_->pinIterator(decision->instance));
    while (pins->hasNext()) {
      const sta::Pin* pin = pins->next();
      const sta::Net* net = network_->net(pin);
      if (net == nullptr) {
        continue;
      }
      if (touched_nets.contains(net)) {
        conflict = true;
        break;
      }
      nets.push_back(net);
    }
    if (conflict) {
      deferred.push_back(decision->instance);
      continue;
    }
    if (resizer_.replaceCell(decision->instance, decision->cell)) {
      touched_nets.insert(nets.begin(), nets.end());
      ++resized;
    }
  }
  gates = std::move(deferred);
  if (resized > 0) {
    estimate_parasitics_->updateParasitics();
    sta_->findDelays();
  }
  seconds_apply_ += secondsSince(start);
  return resized;
}

void ResourceSharingSizingPolicy::iterate()
{
  scene_ = sta_->cmdScene();
  dcalc_ap_ = scene_->dcalcAnalysisPtIndex(max_);
  period_ = clockPeriod();
  sta_->findRequireds();
  float wns = sta::delayAsFloat(sta_->worstSlack(max_));
  float tns = sta::delayAsFloat(sta_->totalNegativeSlack(max_));
  if (period_ <= 0.0 || wns >= 0.0f) {
    markRunComplete(true);
    return;
  }
  const int violating = violatingEndpointCount();
  const double mean_violation
      = violating > 0 ? -static_cast<double>(tns) / violating : 0.0;
  if (mean_violation < rs_config_.start_violation_ratio * period_) {
    logger_->info(RSZ,
                  3307,
                  "RS_SIZING: skipped; mean violation {} of {} endpoints is "
                  "below {:.0f}% of the period.",
                  sta::delayAsString(static_cast<float>(mean_violation), 3, sta_),
                  violating,
                  100.0 * rs_config_.start_violation_ratio);
    markRunComplete(true);
    return;
  }
  collectSizableGates(true);
  power_budget_ = rs_config_.power_budget_ratio * designPower();
  logger_->info(RSZ,
                3303,
                "RS_SIZING: {} sizable gates, period {}, WNS {} TNS {}.",
                sizable_gates_.size(),
                sta::delayAsString(static_cast<float>(period_), 3, sta_),
                sta::delayAsString(wns, 3, sta_),
                sta::delayAsString(tns, 1, sta_));

  int total_resized = 0;
  const float initial_tns = tns;
  double first_efficiency = 0.0;
  const bool lagrangian
      = rs_config_.multipliers
        == ResourceSharingSizingConfig::Multipliers::kLagrangian;
  for (int iteration = 0; iteration < rs_config_.iterations; ++iteration) {
    auto weights_start = std::chrono::steady_clock::now();
    if (iteration > 0) {
      collectSizableGates(false);
    }
    buildLevelOrder();
    if (lagrangian) {
      updateLagrangeMultipliers(iteration);
      normalizeLagrangeMultipliers();
    }
    if (!lagrangian || rs_config_.check_weights) {
      computeLogWeights();
      accumulateLogWeights(iteration);
    }
    if (rs_config_.check_weights) {
      checkWeights();
    }
    // Shift so the most critical arc weights stay near one.
    const std::vector<double>& forward = rs_config_.accumulate_weights
                                             ? log_forward_mean_
                                             : log_forward_;
    const std::vector<double>& backward = rs_config_.accumulate_weights
                                              ? log_backward_mean_
                                              : log_backward_;
    // Shift by the most critical data vertex; the clock tree carries the
    // weight of every path.
    log_weight_shift_ = kNegativeInfinity;
    for (sta::Vertex* vertex : level_order_) {
      if (sta_->isClock(vertex->pin(), sta_->cmdMode())) {
        continue;
      }
      const size_t base = 2 * static_cast<size_t>(graph_->id(vertex));
      for (size_t rf = 0; rf < 2 && base + rf < forward.size(); ++rf) {
        if (forward[base + rf] != kNegativeInfinity
            && backward[base + rf] != kNegativeInfinity) {
          log_weight_shift_ = std::max(
              log_weight_shift_, forward[base + rf] + backward[base + rf]);
        }
      }
    }
    if (!lagrangian && log_weight_shift_ == kNegativeInfinity) {
      break;
    }
    // Power usage relative to the budget; a path at zero slack has log
    // weight 0, a power usage of one the same weight.
    const double power = designPower();
    log_power_weight_ = rs_config_.gamma * (power / power_budget_ - 1.0);
    seconds_weights_ += secondsSince(weights_start);

    resizer_.journalBegin();
    int resized = 0;
    std::vector<sta::Instance*> gates;
    for (sta::Instance* instance : sizable_gates_) {
      if (gateCriticality(instance) >= rs_config_.critical_weight_ratio) {
        gates.push_back(instance);
      }
    }
    const size_t critical_gates = gates.size();
    for (int round = 0;
         round < rs_config_.rounds_per_iteration && !gates.empty();
         ++round) {
      resized += runRound(gates);
    }
    if (rs_config_.buffering) {
      const int inserted = runBufferRound();
      buffers_inserted_ += inserted;
      resized += inserted;
    }
    auto timing_start = std::chrono::steady_clock::now();
    sta_->findRequireds();
    const float new_wns = sta::delayAsFloat(sta_->worstSlack(max_));
    const float new_tns = sta::delayAsFloat(sta_->totalNegativeSlack(max_));
    const bool tns_better = new_tns > tns;
    const bool wns_better = new_wns > wns;
    const double new_power = designPower();
    const double power_increase
        = std::max(1e-3, (new_power - power) / power_budget_);
    const double efficiency
        = (static_cast<double>(new_tns) - tns) / std::abs(initial_tns)
          / power_increase;
    bool keep = resized > 0
                && ((tns_better && new_wns >= wns - rs_config_.wns_loss_limit)
                    || (wns_better
                        && new_tns >= tns - rs_config_.tns_loss_fraction
                                                * std::abs(tns)));
    bool stop = false;
    if (keep && first_efficiency > 0.0
        && efficiency < rs_config_.min_efficiency_ratio * first_efficiency) {
      keep = false;
      stop = true;
    }
    if (keep && first_efficiency == 0.0) {
      first_efficiency = efficiency;
    }
    logger_->info(RSZ,
                  3304,
                  "RS_SIZING iteration {}: critical {} resized {} power {:.4f} "
                  "WNS {} TNS {} efficiency {:.3g} {}",
                  iteration + 1,
                  critical_gates,
                  resized,
                  designPower() / power_budget_ * rs_config_.power_budget_ratio,
                  sta::delayAsString(new_wns, 3, sta_),
                  sta::delayAsString(new_tns, 1, sta_),
                  efficiency,
                  keep ? "kept" : "undone");
    seconds_timing_ += secondsSince(timing_start);
    timing_start = std::chrono::steady_clock::now();
    if (keep) {
      resizer_.journalEnd();
      wns = new_wns;
      tns = new_tns;
      total_resized += resized;
    } else {
      resizer_.journalRestore();
      if (resized == 0 || stop) {
        seconds_timing_ += secondsSince(timing_start);
        break;
      }
    }
    seconds_timing_ += secondsSince(timing_start);
    if (wns >= 0.0f) {
      break;
    }
  }
  logger_->info(RSZ,
                3305,
                "RS_SIZING: resized {} gates (buffers {}); seconds: weights "
                "{:.2f} snapshots {:.2f} oracle {:.2f} apply {:.2f} timing "
                "{:.2f} buffering {:.2f}.",
                total_resized,
                buffers_inserted_,
                seconds_weights_,
                seconds_snapshots_,
                seconds_oracle_,
                seconds_apply_,
                seconds_timing_,
                seconds_buffering_);
  markRunComplete(true);
}

}  // namespace rsz
