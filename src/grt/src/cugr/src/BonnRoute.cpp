// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "BonnRoute.h"

#include <omp.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BonnRouteCapacity.h"
#include "BonnRouteResources.h"
#include "Design.h"
#include "GRNet.h"
#include "GRTree.h"
#include "GridGraph.h"
#include "Layers.h"
#include "geo.h"
#include "grt/BonnRouteExtraSpace.h"
#include "odb/db.h"
#include "utl/Logger.h"
#include "utl/timer.h"

namespace grt {

using utl::GRT;

namespace {
// Keep exp(ln y - shift) inside the double range.
constexpr double kRenormalizeAbove = 600.0;
constexpr double kRenormalizeTo = 300.0;
// Destructive interference size of the x86-64/AArch64 targets we build for.
constexpr std::size_t kCacheLine = 64;
// Stage and phase timing of the reports (wall clock).
using SteadyClock = std::chrono::steady_clock;
double secondsBetween(const SteadyClock::time_point from,
                      const SteadyClock::time_point to)
{
  return std::chrono::duration<double>(to - from).count();
}
}  // namespace

void BonnRouteParams::applyProfile(const int which)
{
  profile = (which == 1 || which == 2) ? which : 0;
  // Both: Algorithm 2 with t = 125, epsilon = 1 (Sec. 2.3), pin vertex
  // sets (Sec. 2.1) and the Sec. 2.5 capacity estimate, goal-oriented
  // Dijkstra in the oracle (Sec. 2.2).
  phases = 125;
  epsilon = 1.0;
  pin_sets = true;
  capacity_model = 1;
  goal_orientation = 1;
  gcell_scale = 1;
  if (profile >= 1) {
    // 2013-reimpl and 2021 (the 2021 paper changes DR only): goal orientation
    // with landmarks (Sec. 2.2), volatility-tolerant concurrent price
    // updates (Sec. 2.3) and the extra space assignment with power and
    // yield resources (Sec. 2.1, Fig. 1).
    landmarks = 8;
    concurrent_prices = true;
    extra_space = true;
  } else {
    // 2013-improved: the l1 bound (faster than landmarks here) and
    // deterministic frozen-price batches; extra space stays off until it
    // has been measured.
    landmarks = 0;
    concurrent_prices = false;
    extra_space = false;
  }
}

const char* BonnRouteParams::profileName() const
{
  return profile == 2 ? "2021"
                      : (profile == 1 ? "2013-reimpl" : "2013-improved");
}

bool BonnRouteParams::set(const std::string& name, const double value)
{
  if (name == "profile") {
    applyProfile(static_cast<int>(value));
  } else if (name == "phases") {
    phases = std::max(1, static_cast<int>(value));
  } else if (name == "epsilon") {
    epsilon = value;
  } else if (name == "reuse_tolerance") {
    reuse_tolerance = value;
  } else if (name == "objective_slack") {
    objective_slack = value;
  } else if (name == "wire_weight") {
    wire_weight = value;
  } else if (name == "via_weight") {
    via_weight = value;
  } else if (name == "margin") {
    margin = static_cast<int>(value);
  } else if (name == "margin_ratio") {
    margin_ratio = value;
  } else if (name == "min_capacity") {
    min_capacity = value;
  } else if (name == "batch_size") {
    batch_size = std::max(1, static_cast<int>(value));
  } else if (name == "threads") {
    threads = std::max(1, static_cast<int>(value));
  } else if (name == "seed") {
    seed = static_cast<unsigned>(value);
  } else if (name == "repair_passes") {
    repair_passes = static_cast<int>(value);
  } else if (name == "pin_sets") {
    pin_sets = value != 0.0;
  } else if (name == "full_graph") {
    full_graph = value != 0.0;
  } else if (name == "maze_cleanup") {
    maze_cleanup = value != 0.0;
  } else if (name == "capacity_model") {
    capacity_model = std::clamp(static_cast<int>(value), 0, 1);
  } else if (name == "capacity_single_tile") {
    capacity_single_tile = value != 0.0;
  } else if (name == "capacity_intra_tile") {
    capacity_intra_tile = value != 0.0;
  } else if (name == "capacity_stacked_vias") {
    capacity_stacked_vias = value != 0.0;
  } else if (name == "capacity_blockage_extension") {
    capacity_blockage_extension = std::max(0.0, value);
  } else if (name == "capacity_single_tile_layers") {
    capacity_single_tile_layers = std::max(1, static_cast<int>(value));
  } else if (name == "capacity_single_tile_margin") {
    capacity_single_tile_margin = std::max(0.0, value);
  } else if (name == "goal_orientation") {
    goal_orientation = std::clamp(static_cast<int>(value), 0, 1);
  } else if (name == "landmarks") {
    landmarks = std::clamp(static_cast<int>(value), 0, 64);
  } else if (name == "landmark_refresh") {
    landmark_refresh = std::max(1, static_cast<int>(value));
  } else if (name == "goal_check") {
    goal_check = value != 0.0;
  } else if (name == "concurrent_prices") {
    concurrent_prices = value != 0.0;
  } else if (name == "gcell_scale") {
    gcell_scale = std::clamp(static_cast<int>(value), 1, 64);
  } else if (name == "extra_space") {
    extra_space = value != 0.0;
  } else if (name == "extra_space_step") {
    extra_space_step = std::max(0.0, value);
  } else if (name == "extra_space_max") {
    extra_space_max = std::max(0.0, value);
  } else if (name == "power_target") {
    power_target = value;
  } else if (name == "yield_target") {
    yield_target = value;
  } else if (name == "coupling_fraction") {
    coupling_fraction = std::clamp(value, 0.0, 1.0);
  } else if (name == "resource_adapt") {
    resource_adapt = value != 0.0;
  } else if (name == "verbose") {
    verbose = value != 0.0;
  } else {
    return false;
  }
  return true;
}

// Per-thread Dijkstra state over the search region, reset by epoch stamps.
// Each thread's instance is cache-line aligned: the heap's end pointer is
// written on every push and pop, and must not share a line with another
// thread's vector headers (std::vector<Scratch> stores them contiguously).
struct alignas(kCacheLine) BonnRouteGlobal::Scratch
{
  std::vector<double> dist;
  std::vector<int32_t> parent;
  std::vector<uint32_t> stamp;
  std::vector<uint32_t> tree_stamp;
  std::vector<uint32_t> term_stamp;
  std::vector<int32_t> term;
  std::vector<std::pair<double, int64_t>> heap;
  std::vector<GoalEntry> goal_heap;
  uint32_t epoch = 0;
  int64_t pops = 0;
  int64_t settled = 0;
  // Wall seconds spent in buildCandidate by this thread.
  double build_seconds = 0.0;

  // n vertices; dist / parent / stamp are per (vertex, via state) node.
  void prepare(size_t n)
  {
    if (tree_stamp.size() < n) {
      dist.resize(3 * n);
      parent.resize(3 * n);
      stamp.resize(3 * n, 0);
      tree_stamp.resize(n, 0);
      term_stamp.resize(n, 0);
      term.resize(n);
    }
    if (++epoch == 0) {
      std::ranges::fill(stamp, 0);
      std::ranges::fill(tree_stamp, 0);
      std::ranges::fill(term_stamp, 0);
      epoch = 1;
    }
    heap.clear();
    goal_heap.clear();
  }
};

BonnRouteGlobal::BonnRouteGlobal(GridGraph* grid_graph,
                                 std::vector<GRNet*> nets,
                                 const BonnRouteParams& params,
                                 utl::Logger* logger)
    : grid_graph_(grid_graph),
      nets_(std::move(nets)),
      params_(params),
      logger_(logger)
{
  x_size_ = grid_graph_->getXSize();
  y_size_ = grid_graph_->getYSize();
  num_layers_ = grid_graph_->getNumLayers();
  m2_pitch_ = std::max(1, grid_graph_->getM2Pitch());
  via_demand_ = grid_graph_->isViaDemandEnabled();
  via_pad_lower_.assign(num_layers_, 0.0);
  via_pad_upper_.assign(num_layers_, 0.0);
  for (int l = 0; l + 1 < num_layers_; l++) {
    via_pad_lower_[l] = grid_graph_->getStackPadLength(l, l, l + 1);
    via_pad_upper_[l] = grid_graph_->getStackPadLength(l + 1, l, l + 1);
  }
  // Same expressions as GridGraph's getters, evaluated once.
  for (int direction = 0; direction < 2; direction++) {
    const int size = grid_graph_->getSize(direction);
    wire_objective_[direction].assign(std::max(size, 1), 0.0);
    edge_length_[direction].assign(std::max(size, 1), 0.0);
    for (int c = 0; c + 1 < size; c++) {
      edge_length_[direction][c] = grid_graph_->getEdgeLength(direction, c);
      wire_objective_[direction][c] = params_.wire_weight
                                      * grid_graph_->getEdgeLength(direction, c)
                                      / m2_pitch_;
    }
  }
  layer_direction_.resize(num_layers_);
  pad_flank_.resize(num_layers_);
  for (int l = 0; l < num_layers_; l++) {
    const int direction = grid_graph_->getLayerDirection(l);
    layer_direction_[l] = direction;
    const int size = grid_graph_->getSize(direction);
    pad_flank_[l].resize(size);
    for (int c = 0; c < size; c++) {
      const int lower_len
          = c > 0 ? grid_graph_->getEdgeLength(direction, c - 1) : 0;
      const int higher_len
          = c < size - 1 ? grid_graph_->getEdgeLength(direction, c) : 0;
      pad_flank_[l][c]
          = {lower_len + higher_len, lower_len > 0, higher_len > 0};
    }
  }
  buildViaSteps();
}

BonnRouteGlobal::~BonnRouteGlobal() = default;

double BonnRouteGlobal::wireObjective(const int layer,
                                      const int x,
                                      const int y) const
{
  const int direction = layer_direction_[layer];
  return wire_objective_[direction][direction == MetalLayer::H ? x : y];
}

namespace {
// Relaxed atomic access to a double shared between oracle threads with
// concurrent_prices (GCC/Clang builtins; plain loads and stores on x86).
double loadRelaxed(const double& v)
{
  double out;
  __atomic_load(&v, &out, __ATOMIC_RELAXED);
  return out;
}
// v = max(v, value); returns nothing.
void storeMax(double& v, const double value)
{
  double cur = loadRelaxed(v);
  double desired = value;
  while (cur < desired
         && !__atomic_compare_exchange(
             &v, &cur, &desired, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
  }
}
// v += inc; returns the new value.
double addRelaxed(double& v, const double inc)
{
  double cur = loadRelaxed(v);
  double next = cur + inc;
  while (!__atomic_compare_exchange(
      &v, &cur, &next, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
    next = cur + inc;
  }
  return next;
}
}  // namespace

double BonnRouteGlobal::price(const uint32_t r) const
{
  return loadRelaxed(weight_[r]);
}

double BonnRouteGlobal::objectiveWeight() const
{
  return loadRelaxed(objective_weight_);
}

double BonnRouteGlobal::wireCost(const int layer,
                                 const int x,
                                 const int y,
                                 const double factor) const
{
  return objectiveWeight() * wireObjective(layer, x, y)
         + price(edgeIndex(layer, x, y)) * factor;
}

double BonnRouteGlobal::padCost(const int layer,
                                const int x,
                                const int y,
                                const double pad,
                                const GRNet* net) const
{
  return padCostFactor(layer, x, y, pad, net->getNdrCost(layer));
}

double BonnRouteGlobal::padCostFactor(const int layer,
                                      const int x,
                                      const int y,
                                      const double pad,
                                      const double factor) const
{
  if (!via_demand_ || pad <= 0.0) {
    return 0.0;
  }
  const int direction = layer_direction_[layer];
  const PadFlank& flank = pad_flank_[layer][direction == MetalLayer::H ? x : y];
  if (flank.sum == 0) {
    return 0.0;
  }
  const double demand = pad / flank.sum * factor;
  double cost = 0.0;
  if (flank.lower) {
    const int lx = direction == MetalLayer::H ? x - 1 : x;
    const int ly = direction == MetalLayer::H ? y : y - 1;
    cost += price(edgeIndex(layer, lx, ly)) * demand;
  }
  if (flank.higher) {
    cost += price(edgeIndex(layer, x, y)) * demand;
  }
  return cost;
}

void BonnRouteGlobal::prepareNets()
{
  states_.clear();
  states_.resize(nets_.size());
  for (size_t i = 0; i < nets_.size(); i++) {
    NetState& state = states_[i];
    GRNet* net = nets_[i];
    state.net = net;
    const LayerRange& range = net->getLayerRange();
    const int top_limit = std::min(range.max_layer, num_layers_ - 1);
    // Also records CUGR's preferred access point per pin (the fallback).
    const AccessPointMap selected = grid_graph_->selectAccessPoints(net);
    if (selected.size() <= 1) {
      state.trivial = true;
      if (!selected.empty()) {
        const auto& [point, layers] = *selected.begin();
        state.trivial_point = point;
        state.trivial_top = std::clamp(
            std::max(layers.high(), range.min_layer), 0, top_limit);
        state.trivial_low
            = std::min(std::max(layers.low(), 0), state.trivial_top);
      } else {
        state.trivial_low = state.trivial_top = -1;
      }
      continue;
    }
    const auto& preferred = net->getPreferredAccessPoints();
    for (int p = 0; p < net->getNumPins(); p++) {
      Terminal t;
      t.pin_index = p;
      if (params_.pin_sets) {
        // Sec. 2.1: "Each pin p has shapes in one or more tiles, and is
        // represented by the respective set V_p of vertices." The tiles of
        // the pin's access points are those the detailed router can use;
        // a pin without access points uses the tiles of its shapes.
        for (const GRPoint& g : grid_graph_->pinAccessCells(net, p)) {
          t.cells.push_back(
              {std::clamp(g.getLayerIdx(), 0, top_limit), g.x(), g.y()});
        }
        if (t.cells.empty()) {
          for (const GRPoint& g : net->getPinAccessPoints()[p]) {
            t.cells.push_back(
                {std::clamp(g.getLayerIdx(), 0, top_limit), g.x(), g.y()});
          }
        }
      }
      if (t.cells.empty()) {
        if (auto it = preferred.find(p); it != preferred.end()) {
          t.cells.push_back({std::clamp(it->second.layers.low(), 0, top_limit),
                             it->second.point.x(),
                             it->second.point.y()});
        }
      }
      std::ranges::sort(t.cells);
      t.cells.erase(std::ranges::unique(t.cells).begin(), t.cells.end());
      if (!t.cells.empty()) {
        state.terminals.push_back(std::move(t));
      }
    }
    if (state.terminals.size() <= 1) {
      state.trivial = true;
      const auto& [point, layers] = *selected.begin();
      state.trivial_point = point;
      state.trivial_top
          = std::clamp(std::max(layers.high(), range.min_layer), 0, top_limit);
      state.trivial_low
          = std::min(std::max(layers.low(), 0), state.trivial_top);
      state.terminals.clear();
      continue;
    }
    // Driver first: path composition grows the tree from it.
    if (const auto driver = net->getDriverAccessPoint()) {
      auto it = std::ranges::find_if(state.terminals, [&](const Terminal& t) {
        return std::ranges::any_of(t.cells, [&](const auto& c) {
          return c[1] == driver->x() && c[2] == driver->y();
        });
      });
      if (it != state.terminals.end()) {
        std::rotate(state.terminals.begin(), it, it + 1);
      }
    }
    int lx = std::numeric_limits<int>::max();
    int ly = lx;
    int hx = std::numeric_limits<int>::min();
    int hy = hx;
    for (const Terminal& t : state.terminals) {
      for (const auto& c : t.cells) {
        lx = std::min(lx, c[1]);
        hx = std::max(hx, c[1]);
        ly = std::min(ly, c[2]);
        hy = std::max(hy, c[2]);
      }
    }
    const int extent = std::max(hx - lx, hy - ly);
    const int margin
        = params_.margin
          + static_cast<int>(std::ceil(params_.margin_ratio * extent));
    state.region = BoxT(std::max(0, lx - margin),
                        std::max(0, ly - margin),
                        std::min(x_size_ - 1, hx + margin),
                        std::min(y_size_ - 1, hy + margin));
    order_.push_back(static_cast<int>(i));
  }
}

// Algorithm 1 under the current prices. Search states are (vertex, how a
// via column was entered) so that a stacked via is charged its column pad
// max(upper, lower) on intermediate layers, as GridGraph commits it.
BonnRouteGlobal::OracleResult BonnRouteGlobal::oracle(
    const NetState& state,
    Scratch& scratch,
    const bool objective_only,
    const bool whole_grid,
    std::vector<Segment>* segments) const
{
  if (params_.goal_orientation != 0) {
    return oracleImpl<true>(
        state, scratch, objective_only, whole_grid, segments);
  }
  return oracleImpl<false>(
      state, scratch, objective_only, whole_grid, segments);
}

// kGoal = Sec. 2.2's goal orientation: the same search with the heap keyed
// by dist + h, h a consistent lower bound of the remaining distance to the
// unconnected terminals (goalBound()). Without it the search is plain
// Dijkstra keyed by dist.
template <bool kGoal>
BonnRouteGlobal::OracleResult BonnRouteGlobal::oracleImpl(
    const NetState& state,
    Scratch& scratch,
    const bool objective_only,
    const bool whole_grid,
    std::vector<Segment>* segments) const
{
  OracleResult result;
  const GRNet* net = state.net;
  const BoxT region = (whole_grid || params_.full_graph)
                          ? BoxT(0, 0, x_size_ - 1, y_size_ - 1)
                          : state.region;
  const int x0 = region.lx();
  const int y0 = region.ly();
  const int x1 = region.hx();
  const int y1 = region.hy();
  const int width = x1 - x0 + 1;
  const int height = y1 - y0 + 1;
  const LayerRange& range = net->getLayerRange();
  const int num_layers = std::min(range.max_layer, num_layers_ - 1) + 1;
  const int64_t plane = static_cast<int64_t>(width) * height;
  const size_t n = static_cast<size_t>(plane) * num_layers;
  scratch.prepare(n);
  auto index = [&](int l, int x, int y) -> int64_t {
    return l * plane + static_cast<int64_t>(x - x0) * height + (y - y0);
  };
  auto greater = [](const std::pair<double, int64_t>& a,
                    const std::pair<double, int64_t>& b) { return a > b; };
  // Goal orientation: heap entries carry dist + h and the target-set
  // version their h was computed for.
  [[maybe_unused]] GoalContext goal;
  if constexpr (kGoal) {
    initGoal(goal, state, objective_only, x0, y0, x1, y1, num_layers);
  }
  [[maybe_unused]] auto goal_greater
      = [](const GoalEntry& a, const GoalEntry& b) {
          if (a.key != b.key) {
            return a.key > b.key;
          }
          if (a.dist != b.dist) {
            return a.dist < b.dist;  // deeper first among equal keys
          }
          return a.node > b.node;
        };
  [[maybe_unused]] auto vertex_of = [&](int64_t v) {
    const int l = static_cast<int>(v / plane);
    const int64_t rem = v % plane;
    return std::array<int, 3>{l,
                              x0 + static_cast<int>(rem / height),
                              y0 + static_cast<int>(rem % height)};
  };
  [[maybe_unused]] auto bound = [&](int64_t node) {
    const auto [l, x, y] = vertex_of(node / 3);
    return goalBound(goal, l, x, y);
  };
  auto push = [&](double d, int64_t node) {
    if constexpr (kGoal) {
      scratch.goal_heap.push_back({d + bound(node), d, node, goal.version});
      std::push_heap(
          scratch.goal_heap.begin(), scratch.goal_heap.end(), goal_greater);
    } else {
      scratch.heap.emplace_back(d, node);
      std::push_heap(scratch.heap.begin(), scratch.heap.end(), greater);
    }
  };
  auto relax = [&](int64_t from, int64_t to, double d) {
    if (scratch.stamp[to] != scratch.epoch || d < scratch.dist[to]) {
      scratch.stamp[to] = scratch.epoch;
      scratch.dist[to] = d;
      scratch.parent[to] = static_cast<int32_t>(from);
      push(d, to);
    }
  };

  const int num_terminals = static_cast<int>(state.terminals.size());
  std::vector<std::vector<int>> pins_at;  // terminal lists per marked vertex
  std::vector<int64_t> marked;
  for (int t = 0; t < num_terminals; t++) {
    for (const auto& c : state.terminals[t].cells) {
      if (c[0] >= num_layers) {
        continue;
      }
      const int64_t v = index(c[0], c[1], c[2]);
      if (scratch.term_stamp[v] != scratch.epoch) {
        scratch.term_stamp[v] = scratch.epoch;
        scratch.term[v] = static_cast<int32_t>(marked.size());
        marked.push_back(v);
        pins_at.emplace_back();
      }
      pins_at[scratch.term[v]].push_back(t);
    }
  }
  std::vector<char> connected(num_terminals, 0);
  int remaining = num_terminals;
  std::vector<std::array<int64_t, 2>> edges;  // vertex pairs
  std::vector<int> edge_pins;                 // pin of a clique edge, or -1
  auto add_source = [&](int64_t v) {
    scratch.tree_stamp[v] = scratch.epoch;
    const int64_t node = v * 3;
    scratch.stamp[node] = scratch.epoch;
    scratch.dist[node] = 0.0;
    scratch.parent[node] = -1;
    push(0.0, node);
  };
  // Sec. 2.2: "To find T*, we can add the edges of the clique K(V_p) for
  // each pin p in n with zero cost, and use an approximation algorithm for
  // the Steiner tree problem in weighted graphs." Algorithm 1 runs in G plus
  // these cliques: once a vertex of V_p joins the tree, every vertex of V_p
  // is reached at no cost and becomes a source. The clique edges are kept as
  // pin links; buildCandidate prunes the ones no branch uses.
  std::vector<int64_t> pending;
  auto connect_at = [&](int64_t first) {
    pending.assign(1, first);
    for (size_t k = 0; k < pending.size(); k++) {
      const int64_t v = pending[k];
      if (scratch.term_stamp[v] != scratch.epoch) {
        continue;
      }
      for (const int t : pins_at[scratch.term[v]]) {
        if (connected[t]) {
          continue;
        }
        connected[t] = 1;
        remaining--;
        for (const auto& c : state.terminals[t].cells) {
          if (c[0] >= num_layers) {
            continue;
          }
          const int64_t w = index(c[0], c[1], c[2]);
          if (scratch.tree_stamp[w] != scratch.epoch) {
            edges.push_back({v, w});
            edge_pins.push_back(state.terminals[t].pin_index);
            add_source(w);
            pending.push_back(w);
          }
        }
      }
    }
  };
  auto refresh_targets = [&]() {
    if constexpr (kGoal) {
      setGoalTargets(goal, state, connected, num_layers);
    }
  };
  // Algorithm 1, step 1: the component of the first terminal (driver).
  {
    const auto& c = state.terminals[0].cells.front();
    const int64_t root = index(c[0], c[1], c[2]);
    add_source(root);
    connect_at(root);
    refresh_targets();
  }
  auto heap_empty = [&]() {
    if constexpr (kGoal) {
      return scratch.goal_heap.empty();
    } else {
      return scratch.heap.empty();
    }
  };
  const double obj = objectiveWeight();
  const bool extra = params_.extra_space && !objective_only;
  double last_connection_cost = 0.0;
  while (!heap_empty() && remaining > 0) {
    double d;
    int64_t node;
    if constexpr (kGoal) {
      std::pop_heap(
          scratch.goal_heap.begin(), scratch.goal_heap.end(), goal_greater);
      const GoalEntry entry = scratch.goal_heap.back();
      scratch.goal_heap.pop_back();
      scratch.pops++;
      d = entry.dist;
      node = entry.node;
      if (d > scratch.dist[node]) {
        continue;
      }
      if (entry.version != goal.version) {
        // Fewer targets since the push: h can only have grown. Re-key the
        // entry if it is no longer the minimum.
        const double key = d + bound(node);
        if (key > entry.key && !scratch.goal_heap.empty()
            && goal_greater(GoalEntry{key, d, node, goal.version},
                            scratch.goal_heap.front())) {
          scratch.goal_heap.push_back({key, d, node, goal.version});
          std::push_heap(
              scratch.goal_heap.begin(), scratch.goal_heap.end(), goal_greater);
          continue;
        }
      }
    } else {
      std::pop_heap(scratch.heap.begin(), scratch.heap.end(), greater);
      d = scratch.heap.back().first;
      node = scratch.heap.back().second;
      scratch.heap.pop_back();
      scratch.pops++;
      if (d > scratch.dist[node]) {
        continue;
      }
    }
    scratch.settled++;
    const int64_t u = node / 3;
    const int s = static_cast<int>(node % 3);
    if (scratch.tree_stamp[u] != scratch.epoch
        && scratch.term_stamp[u] == scratch.epoch
        && std::ranges::any_of(pins_at[scratch.term[u]],
                               [&](int t) { return !connected[t]; })) {
      // Algorithm 1, step 5: add the path to the tree, then let the pins
      // it reached join with all their vertices.
      last_connection_cost = d;
      int64_t cur = node;
      std::vector<int64_t> added;
      while (scratch.tree_stamp[cur / 3] != scratch.epoch) {
        const int64_t p = scratch.parent[cur];
        edges.push_back({cur / 3, p / 3});
        edge_pins.push_back(-1);
        scratch.tree_stamp[cur / 3] = scratch.epoch;
        added.push_back(cur / 3);
        cur = p;
      }
      if (segments != nullptr) {
        uint64_t key = 1469598103934665603ULL;
        for (const int64_t v : added) {
          key = (key ^ static_cast<uint64_t>(v)) * 1099511628211ULL;
        }
        segments->push_back({d, key});
      }
      for (const int64_t v : added) {
        add_source(v);
      }
      for (const int64_t v : added) {
        connect_at(v);
      }
      refresh_targets();
      continue;
    }
    const int l = static_cast<int>(u / plane);
    const int64_t rem = u % plane;
    const int x = x0 + static_cast<int>(rem / height);
    const int y = y0 + static_cast<int>(rem % height);
    if (l >= range.min_layer && l <= range.max_layer) {
      const double factor = net->getNdrCost(l);
      auto wire = [&](int ex, int ey) {
        if (objective_only) {
          return wireObjective(l, ex, ey);
        }
        double c = wireCost(l, ex, ey, factor);
        if (extra) {
          // Eq. (1): the edge cost with the optimum extra space.
          c += spaceCost(l, ex, ey, nullptr);
        }
        return c;
      };
      if (layer_direction_[l] == MetalLayer::H) {
        if (x > x0) {
          relax(node, (u - height) * 3, d + wire(x - 1, y));
        }
        if (x < x1) {
          relax(node, (u + height) * 3, d + wire(x, y));
        }
      } else {
        if (y > y0) {
          relax(node, (u - 1) * 3, d + wire(x, y - 1));
        }
        if (y < y1) {
          relax(node, (u + 1) * 3, d + wire(x, y));
        }
      }
    }
    if (l + 1 < num_layers) {
      double c = params_.via_weight;
      if (!objective_only) {
        // Pad on l: lower pad of this via, or the column increment when the
        // column already came up into l.
        const double lower = via_pad_lower_[l];
        const double pad_l = (s == 1 && l > 0)
                                 ? std::max(via_pad_upper_[l - 1], lower)
                                       - via_pad_upper_[l - 1]
                                 : lower;
        c = obj * params_.via_weight + padCost(l, x, y, pad_l, net)
            + padCost(l + 1, x, y, via_pad_upper_[l], net);
      }
      relax(node, (u + plane) * 3 + 1, d + c);
    }
    if (l > 0) {
      double c = params_.via_weight;
      if (!objective_only) {
        const double upper = via_pad_upper_[l - 1];
        const double pad_l
            = (s == 2) ? std::max(upper, via_pad_lower_[l]) - via_pad_lower_[l]
                       : upper;
        c = obj * params_.via_weight + padCost(l, x, y, pad_l, net)
            + padCost(l - 1, x, y, via_pad_lower_[l - 1], net);
      }
      relax(node, (u - plane) * 3 + 2, d + c);
    }
  }
  if (remaining > 0) {
    return result;
  }
  result.search_cost = last_connection_cost;
  auto cell = [&](int64_t v) -> std::array<int, 3> {
    const int l = static_cast<int>(v / plane);
    const int64_t rem = v % plane;
    return {l,
            x0 + static_cast<int>(rem / height),
            y0 + static_cast<int>(rem % height)};
  };
  // A tree the net's oracle returned before maps to the stored candidate;
  // buildCandidate would rebuild it unchanged. With extra space the demand
  // depends on the prices, so those calls always build.
  if (!objective_only && !params_.extra_space) {
    result.keyed = true;
    TreeKey& key = result.key;
    key.edges.clear();
    key.edges.reserve(edges.size());
    for (size_t i = 0; i < edges.size(); i++) {
      const auto a = cell(edges[i][0]);
      const auto b = cell(edges[i][1]);
      uint32_t ka = edgeIndex(a[0], a[1], a[2]);
      uint32_t kb = edgeIndex(b[0], b[1], b[2]);
      if (ka == kb) {
        continue;
      }
      if (ka > kb) {
        std::swap(ka, kb);
      }
      key.edges.push_back({ka, kb, static_cast<uint32_t>(edge_pins[i] + 1)});
    }
    std::ranges::sort(key.edges);
    // One entry per vertex pair, the smallest pin first (as buildCandidate's
    // adjacency keeps it).
    size_t kept = 0;
    for (size_t i = 0; i < key.edges.size(); i++) {
      if (kept > 0 && key.edges[kept - 1][0] == key.edges[i][0]
          && key.edges[kept - 1][1] == key.edges[i][1]) {
        continue;
      }
      key.edges[kept++] = key.edges[i];
    }
    key.edges.resize(kept);
    uint64_t hash = 1469598103934665603ULL;
    for (const auto& e : key.edges) {
      for (const uint32_t v : e) {
        hash = (hash ^ v) * 1099511628211ULL;
      }
    }
    key.hash = hash;
    for (const TreeKey& known : state.tree_keys) {
      if (known.hash == key.hash && known.edges == key.edges) {
        result.known = known.candidate;
        result.candidate.objective = known.objective;
        result.ok = true;
        return result;
      }
    }
  }
  std::vector<std::array<int, 3>> tree_edges;  // (layer, x, y) pairs packed
  tree_edges.reserve(edges.size() * 2);
  for (const auto& [a, b] : edges) {
    tree_edges.push_back(cell(a));
    tree_edges.push_back(cell(b));
  }
  const SteadyClock::time_point build_start = SteadyClock::now();
  const bool built = buildCandidate(
      state, tree_edges, edge_pins, !objective_only, result.candidate);
  scratch.build_seconds += secondsBetween(build_start, SteadyClock::now());
  if (!built) {
    return result;
  }
  result.ok = true;
  return result;
}

// ---------------------------------------------------------------------------
// Goal orientation (Sec. 2.2: "a variant of goal-orientation with landmarks
// [Goldberg and Harrelson 2005]").
//
// Every oracle edge cost is (objective weight) * (wire objective or via
// weight) + (price part: wire weight * NDR factor, or via-pad prices). Prices
// never decrease (applyPrices only adds), so bounds taken from earlier prices
// stay lower bounds. h(v) = max(L1(v), O(v) + A(v)), each term a consistent
// lower bound of the remaining distance to the nearest unconnected target:
//   L1: min over targets of cx|dx| + cy|dy| + cv * vias(dx, dy, l -> l_t),
//       with cx = objective weight * cheapest wire objective in the region
//       + cheapest wire price per direction (phase start) * the net's
//       smallest NDR factor, cv = objective weight * via weight, and
//       vias() the fewest layer changes that reach a wire layer of each
//       needed direction and then the target layer (via_steps_);
//   O:  the objective part of L1 alone;
//   A:  ALT on the price part: max over landmarks of the distance of the
//       landmark distance d_k(v) to the targets' interval [min, max] of
//       d_k(t), with d_k from a Dijkstra over the per-edge minimum price
//       part (smallest NDR factor, cheapest via-entry state).
// The objective and price parts of a path bound O and A separately, so
// O + A is a lower bound; L1, O and A are 1-Lipschitz along every edge with
// respect to its cost, so the max is consistent.

namespace {
constexpr int kExactTargets = 16;
constexpr uint8_t kUnreachable = 255;
}  // namespace

void BonnRouteGlobal::buildViaSteps()
{
  const int L = num_layers_;
  via_steps_.assign(static_cast<size_t>(L) * L * 4 * L * L, kUnreachable);
  for (int lo = 0; lo < L; lo++) {
    for (int hi = lo; hi < L; hi++) {
      // Wire layers lo..hi; vias pass any layer 0..hi.
      auto bits = [&](int l) {
        if (l < lo || l > hi) {
          return 0;
        }
        return layer_direction_[l] == MetalLayer::H ? 1 : 2;
      };
      uint8_t* table
          = via_steps_.data() + (static_cast<size_t>(lo) * L + hi) * 4 * L * L;
      for (int from = 0; from <= hi; from++) {
        // BFS over (layer, directions seen).
        std::vector<int> dist(static_cast<size_t>(L) * 4, -1);
        std::vector<std::pair<int, int>> queue;
        const int start_bits = bits(from);
        dist[from * 4 + start_bits] = 0;
        queue.emplace_back(from, start_bits);
        for (size_t k = 0; k < queue.size(); k++) {
          const auto [l, seen] = queue[k];
          for (const int next : {l - 1, l + 1}) {
            if (next < 0 || next > hi) {
              continue;
            }
            const int nseen = seen | bits(next);
            if (dist[next * 4 + nseen] < 0) {
              dist[next * 4 + nseen] = dist[l * 4 + seen] + 1;
              queue.emplace_back(next, nseen);
            }
          }
        }
        for (int mask = 0; mask < 4; mask++) {
          for (int to = 0; to <= hi; to++) {
            int best = -1;
            for (int seen = 0; seen < 4; seen++) {
              const int dd = dist[to * 4 + seen];
              if ((seen & mask) == mask && dd >= 0 && (best < 0 || dd < best)) {
                best = dd;
              }
            }
            // A direction without a wire layer: no path; fall back to the
            // layer distance (any value is admissible then).
            table[(mask * L + from) * L + to] = static_cast<uint8_t>(
                best >= 0 ? std::min(best, 254) : std::abs(from - to));
          }
        }
      }
    }
  }
}

void BonnRouteGlobal::initGoal(GoalContext& goal,
                               const NetState& state,
                               const bool objective_only,
                               const int x0,
                               const int y0,
                               const int x1,
                               const int y1,
                               const int num_layers) const
{
  const GRNet* net = state.net;
  const LayerRange& range = net->getLayerRange();
  auto cheapest = [&](int direction, int lo, int hi) {
    double best = 0.0;
    bool any = false;
    for (int c = lo; c < hi; c++) {
      const double w = wire_objective_[direction][c];
      best = any ? std::min(best, w) : w;
      any = true;
    }
    return best;
  };
  const double objective_x = cheapest(MetalLayer::H, x0, x1);
  const double objective_y = cheapest(MetalLayer::V, y0, y1);
  const double scale = objective_only ? 1.0 : objectiveWeight();
  goal.ox = scale * objective_x;
  goal.oy = scale * objective_y;
  goal.ov = scale * params_.via_weight;
  goal.cx = goal.ox;
  goal.cy = goal.oy;
  goal.cv = goal.ov;
  if (!objective_only) {
    double factor = std::numeric_limits<double>::infinity();
    for (int l = std::max(0, range.min_layer); l < num_layers; l++) {
      factor = std::min(factor, net->getNdrCost(l));
    }
    if (std::isfinite(factor)) {
      goal.cx += min_weight_[MetalLayer::H] * factor;
      goal.cy += min_weight_[MetalLayer::V] * factor;
    }
  }
  const int lo = std::clamp(range.min_layer, 0, num_layers - 1);
  const int L = num_layers_;
  goal.via_steps
      = via_steps_.data()
        + (static_cast<size_t>(lo) * L + (num_layers - 1)) * 4 * L * L;
  const int K = params_.landmarks;
  goal.alt = nullptr;
  if (!objective_only && K > 0 && landmarks_valid_) {
    goal.alt = landmark_dist_.data();
    goal.alt_scale = std::exp(landmark_shift_ - shift_);
    goal.alt_lo.assign(K, 0.0F);
    goal.alt_hi.assign(K, 0.0F);
  }
  goal.version = 0;
}

void BonnRouteGlobal::setGoalTargets(GoalContext& goal,
                                     const NetState& state,
                                     const std::vector<char>& connected,
                                     const int num_layers) const
{
  goal.targets.clear();
  for (size_t t = 0; t < state.terminals.size(); t++) {
    if (connected[t]) {
      continue;
    }
    for (const auto& c : state.terminals[t].cells) {
      if (c[0] < num_layers) {
        goal.targets.push_back(
            {c[0],
             c[1],
             c[2],
             static_cast<int64_t>(edgeIndex(c[0], c[1], c[2]))});
      }
    }
  }
  goal.exact = goal.targets.size() <= kExactTargets;
  goal.lx = goal.ly = goal.ll = std::numeric_limits<int>::max();
  goal.hx = goal.hy = goal.hl = std::numeric_limits<int>::min();
  for (const GoalTarget& t : goal.targets) {
    goal.lx = std::min(goal.lx, t.x);
    goal.hx = std::max(goal.hx, t.x);
    goal.ly = std::min(goal.ly, t.y);
    goal.hy = std::max(goal.hy, t.y);
    goal.ll = std::min(goal.ll, t.l);
    goal.hl = std::max(goal.hl, t.l);
  }
  if (goal.alt != nullptr) {
    const int K = params_.landmarks;
    for (int k = 0; k < K; k++) {
      float lo = std::numeric_limits<float>::infinity();
      float hi = -std::numeric_limits<float>::infinity();
      for (const GoalTarget& t : goal.targets) {
        const float d = goal.alt[t.vertex * K + k];
        lo = std::min(lo, d);
        hi = std::max(hi, d);
      }
      goal.alt_lo[k] = lo;
      goal.alt_hi[k] = hi;
    }
  }
  goal.version++;
}

double BonnRouteGlobal::goalBound(const GoalContext& goal,
                                  const int l,
                                  const int x,
                                  const int y) const
{
  if (goal.targets.empty()) {
    return 0.0;
  }
  const int L = num_layers_;
  double full;
  double objective;
  if (goal.exact) {
    full = std::numeric_limits<double>::infinity();
    objective = full;
    for (const GoalTarget& t : goal.targets) {
      const int dx = std::abs(x - t.x);
      const int dy = std::abs(y - t.y);
      const int mask = (dx != 0 ? 1 : 0) | (dy != 0 ? 2 : 0);
      const int vias = goal.via_steps[(mask * L + l) * L + t.l];
      full = std::min(full, goal.cx * dx + goal.cy * dy + goal.cv * vias);
      objective
          = std::min(objective, goal.ox * dx + goal.oy * dy + goal.ov * vias);
    }
  } else {
    auto gap = [](int v, int lo, int hi) {
      return v < lo ? lo - v : (v > hi ? v - hi : 0);
    };
    const int dx = gap(x, goal.lx, goal.hx);
    const int dy = gap(y, goal.ly, goal.hy);
    const int dl = gap(l, goal.ll, goal.hl);
    full = goal.cx * dx + goal.cy * dy + goal.cv * dl;
    objective = goal.ox * dx + goal.oy * dy + goal.ov * dl;
  }
  if (goal.alt == nullptr) {
    return full;
  }
  const int K = params_.landmarks;
  const float* d = goal.alt + static_cast<int64_t>(edgeIndex(l, x, y)) * K;
  double alt = 0.0;
  for (int k = 0; k < K; k++) {
    const double dv = d[k];
    const double below = goal.alt_lo[k] - dv;
    const double above = dv - goal.alt_hi[k];
    // Float storage: give back the rounding of both distances.
    const double slack
        = 1.2e-7
          * (std::abs(dv)
             + std::max(std::abs(goal.alt_lo[k]), std::abs(goal.alt_hi[k])));
    alt = std::max(alt, std::max(below, above) - slack);
  }
  return std::max(full, objective + alt * goal.alt_scale);
}

// Lower bound of the price part of a via between l and l + 1 at (x, y):
// the cheaper direction and via-entry state, smallest NDR factors.
double BonnRouteGlobal::viaPriceBound(const int l,
                                      const int x,
                                      const int y) const
{
  const double lower = via_pad_lower_[l];
  const double up_l = l > 0 ? std::min(lower,
                                       std::max(via_pad_upper_[l - 1], lower)
                                           - via_pad_upper_[l - 1])
                            : lower;
  const double up
      = padCostFactor(l, x, y, up_l, ndr_min_[l])
        + padCostFactor(l + 1, x, y, via_pad_upper_[l], ndr_min_[l + 1]);
  const double upper = via_pad_upper_[l];
  const double down_h = std::min(
      upper, std::max(upper, via_pad_lower_[l + 1]) - via_pad_lower_[l + 1]);
  const double down = padCostFactor(l + 1, x, y, down_h, ndr_min_[l + 1])
                      + padCostFactor(l, x, y, lower, ndr_min_[l]);
  return std::min(up, down);
}

void BonnRouteGlobal::landmarkDijkstra(const int64_t source,
                                       std::vector<double>& dist) const
{
  const int64_t X = x_size_;
  const int64_t Y = y_size_;
  const int64_t P = X * Y;
  dist.assign(P * num_layers_, std::numeric_limits<double>::infinity());
  using Entry = std::pair<double, int64_t>;
  std::vector<Entry> heap;
  auto greater = [](const Entry& a, const Entry& b) { return a > b; };
  dist[source] = 0.0;
  heap.emplace_back(0.0, source);
  while (!heap.empty()) {
    std::pop_heap(heap.begin(), heap.end(), greater);
    const double d = heap.back().first;
    const int64_t v = heap.back().second;
    heap.pop_back();
    if (d > dist[v]) {
      continue;
    }
    const int y = static_cast<int>(v % Y);
    const int x = static_cast<int>((v / Y) % X);
    const int l = static_cast<int>(v / P);
    auto relax = [&](int64_t w, double c) {
      if (d + c < dist[w]) {
        dist[w] = d + c;
        heap.emplace_back(d + c, w);
        std::push_heap(heap.begin(), heap.end(), greater);
      }
    };
    if (l >= min_layer_) {
      const double f = ndr_min_[l];
      if (layer_direction_[l] == MetalLayer::H) {
        if (x > 0) {
          relax(v - Y, weight_[edgeIndex(l, x - 1, y)] * f);
        }
        if (x + 1 < X) {
          relax(v + Y, weight_[edgeIndex(l, x, y)] * f);
        }
      } else {
        if (y > 0) {
          relax(v - 1, weight_[edgeIndex(l, x, y - 1)] * f);
        }
        if (y + 1 < Y) {
          relax(v + 1, weight_[edgeIndex(l, x, y)] * f);
        }
      }
    }
    if (l + 1 < num_layers_) {
      relax(v + P, viaPriceBound(l, x, y));
    }
    if (l > 0) {
      relax(v - P, viaPriceBound(l - 1, x, y));
    }
  }
}

void BonnRouteGlobal::updateGoalBounds(const bool landmarks)
{
  // Cheapest current wire price per direction on the priced layers.
  min_weight_ = {std::numeric_limits<double>::infinity(),
                 std::numeric_limits<double>::infinity()};
  for (int l = min_layer_; l < num_layers_; l++) {
    const int direction = layer_direction_[l];
    const int xe = direction == MetalLayer::H ? x_size_ - 1 : x_size_;
    const int ye = direction == MetalLayer::H ? y_size_ : y_size_ - 1;
    double best = min_weight_[direction];
    for (int x = 0; x < xe; x++) {
      for (int y = 0; y < ye; y++) {
        best = std::min(best, weight_[edgeIndex(l, x, y)]);
      }
    }
    min_weight_[direction] = best;
  }
  for (double& w : min_weight_) {
    if (!std::isfinite(w)) {
      w = 0.0;
    }
  }
  const int K = params_.landmarks;
  if (!landmarks || K <= 0) {
    return;
  }
  utl::Timer timer;
  const int64_t R = static_cast<int64_t>(num_layers_) * x_size_ * y_size_;
  std::vector<std::vector<double>> dist(K);
  if (landmark_vertices_.empty()) {
    // Spread over the die boundary (where farthest-point selection ends up
    // on a nearly uniform grid), alternating the lowest priced and the top
    // layer.
    const int64_t perimeter = 2 * static_cast<int64_t>(x_size_ - 1)
                              + 2 * static_cast<int64_t>(y_size_ - 1);
    for (int k = 0; k < K; k++) {
      int64_t p = perimeter > 0 ? perimeter * k / K : 0;
      int x = 0;
      int y = 0;
      if (p < x_size_ - 1) {
        x = static_cast<int>(p);
      } else if ((p -= x_size_ - 1) < y_size_ - 1) {
        x = x_size_ - 1;
        y = static_cast<int>(p);
      } else if ((p -= y_size_ - 1) < x_size_ - 1) {
        x = x_size_ - 1 - static_cast<int>(p);
        y = y_size_ - 1;
      } else {
        p -= x_size_ - 1;
        y = y_size_ - 1 - static_cast<int>(p);
      }
      const int l = k % 2 == 0 ? std::min(min_layer_, num_layers_ - 1)
                               : num_layers_ - 1;
      landmark_vertices_.push_back(edgeIndex(l, x, y));
    }
  }
  {
    const int threads = std::max(1, params_.threads);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (int k = 0; k < K; k++) {
      landmarkDijkstra(landmark_vertices_[k], dist[k]);
    }
  }
  landmark_dist_.resize(R * K);
  for (int64_t v = 0; v < R; v++) {
    for (int k = 0; k < K; k++) {
      // Unreachable vertices give no bound.
      const double d = dist[k][v];
      landmark_dist_[v * K + k]
          = std::isfinite(d) ? static_cast<float>(d) : 0.0F;
    }
  }
  // An unreachable vertex has distance 0 for every landmark; with a finite
  // target interval this could overstate h. Disable ALT if any exist.
  bool all_finite = true;
  for (int k = 0; k < K && all_finite; k++) {
    for (int64_t v = 0; v < R; v++) {
      if (!std::isfinite(dist[k][v])) {
        all_finite = false;
        break;
      }
    }
  }
  landmarks_valid_ = all_finite;
  landmark_shift_ = shift_;
  landmark_updates_++;
  landmark_time_ += timer.elapsed();
}

void BonnRouteGlobal::checkGoal(const NetState& state,
                                Scratch& reference,
                                const bool objective_only,
                                const bool whole_grid,
                                const bool ok,
                                const std::vector<Segment>& segments)
{
  std::vector<Segment> expected;
  const OracleResult r = oracleImpl<false>(
      state, reference, objective_only, whole_grid, &expected);
  goal_checked_++;
  bool mismatch = r.ok != ok;
  const size_t n = std::min(segments.size(), expected.size());
  for (size_t i = 0; i < n && !mismatch; i++) {
    const double a = segments[i].cost;
    const double b = expected[i].cost;
    goal_segments_++;
    if (std::abs(a - b) > 1e-9 * std::max(1.0, std::abs(b))) {
      mismatch = true;
      if (goal_mismatches_ < 5) {
        logger_->report(
            "BonnRoute goal check: net {} connection {} cost {:.12g}, "
            "Dijkstra {:.12g}.",
            state.net->getName(),
            i,
            a,
            b);
      }
      break;
    }
    if (segments[i].path != expected[i].path) {
      goal_ties_++;  // equal cost, other path: later connections differ
      break;
    }
  }
  if (mismatch) {
    goal_mismatches_++;
  }
}

bool BonnRouteGlobal::buildCandidate(
    const NetState& state,
    const std::vector<std::array<int, 3>>& tree_edges,
    const std::vector<int>& edge_pins,
    const bool need_usage,
    Candidate& candidate) const
{
  auto key = [&](int l, int x, int y) -> uint64_t {
    return (static_cast<uint64_t>(l) * x_size_ + x) * y_size_ + y;
  };
  auto decode = [&](uint64_t k) -> std::array<int, 3> {
    const int y = static_cast<int>(k % y_size_);
    k /= y_size_;
    const int x = static_cast<int>(k % x_size_);
    return {static_cast<int>(k / x_size_), x, y};
  };
  // Neighbor and the GRNet pin of the edge (-1 = wire or via).
  using Link = std::pair<uint64_t, int>;
  std::unordered_map<uint64_t, std::vector<Link>> adjacency;
  for (size_t i = 0; i + 1 < tree_edges.size(); i += 2) {
    const auto& a = tree_edges[i];
    const auto& b = tree_edges[i + 1];
    const uint64_t ka = key(a[0], a[1], a[2]);
    const uint64_t kb = key(b[0], b[1], b[2]);
    if (ka == kb) {
      continue;
    }
    const int pin = edge_pins[i / 2];
    adjacency[ka].emplace_back(kb, pin);
    adjacency[kb].emplace_back(ka, pin);
  }
  for (auto& [k, list] : adjacency) {
    // A wire or via (-1) sorts before a pin link between the same cells.
    std::ranges::sort(list);
    list.erase(std::ranges::unique(list, {}, &Link::first).begin(), list.end());
  }
  const int num_terminals = static_cast<int>(state.terminals.size());
  std::vector<std::vector<uint64_t>> terminal_cells(num_terminals);
  std::unordered_map<uint64_t, std::vector<int>> terminals_at;
  for (int t = 0; t < num_terminals; t++) {
    for (const auto& c : state.terminals[t].cells) {
      const uint64_t k = key(c[0], c[1], c[2]);
      terminal_cells[t].push_back(k);
      terminals_at[k].push_back(t);
    }
  }
  auto in_forest = [&](uint64_t k) { return adjacency.contains(k); };
  // Prune clique edges that end in a cell no branch uses, unless that cell
  // is the only one of some pin in the forest.
  std::vector<int> present(num_terminals, 0);
  for (int t = 0; t < num_terminals; t++) {
    for (const uint64_t k : terminal_cells[t]) {
      present[t] += in_forest(k);
    }
  }
  auto prunable = [&](uint64_t v) {
    auto it = adjacency.find(v);
    if (it == adjacency.end() || it->second.size() != 1
        || it->second.front().second < 0) {
      return false;
    }
    auto at = terminals_at.find(v);
    if (at != terminals_at.end()) {
      for (const int t : at->second) {
        if (present[t] <= 1) {
          return false;
        }
      }
    }
    return true;
  };
  // Without any edge left, every pin contains this one cell.
  const auto& start = state.terminals[0].cells.front();
  uint64_t sole = key(start[0], start[1], start[2]);
  std::vector<uint64_t> leaves;
  for (const auto& [v, list] : adjacency) {
    if (prunable(v)) {
      leaves.push_back(v);
    }
  }
  std::ranges::sort(leaves);  // deterministic order
  for (size_t i = 0; i < leaves.size(); i++) {
    const uint64_t v = leaves[i];
    if (!prunable(v)) {
      continue;
    }
    const uint64_t w = adjacency[v].front().first;
    std::erase_if(adjacency[w],
                  [&](const Link& link) { return link.first == v; });
    adjacency.erase(v);
    if (auto at = terminals_at.find(v); at != terminals_at.end()) {
      for (const int t : at->second) {
        present[t]--;
      }
    }
    if (adjacency[w].empty()) {
      adjacency.erase(w);
      sole = w;
    } else if (prunable(w)) {
      leaves.push_back(w);
    }
  }
  // Representative cell of each pin: one in the forest, preferably with a
  // wire or via.
  std::vector<uint64_t> representative(num_terminals);
  const bool single_cell = adjacency.empty();
  for (int t = 0; t < num_terminals; t++) {
    int best = -1;
    for (size_t i = 0; i < terminal_cells[t].size(); i++) {
      const uint64_t k = terminal_cells[t][i];
      if (single_cell ? k != sole : !in_forest(k)) {
        continue;
      }
      const bool physical
          = !single_cell
            && std::ranges::any_of(adjacency.at(k),
                                   [](const Link& l) { return l.second < 0; });
      if (best < 0 || physical) {
        best = static_cast<int>(i);
        if (physical) {
          break;
        }
      }
    }
    if (best < 0) {
      return false;
    }
    representative[t] = terminal_cells[t][best];
  }
  std::unordered_map<uint64_t, char> pinned;
  candidate.pin_cells.clear();
  for (const uint64_t k : representative) {
    pinned[k] = 1;
    const auto c = decode(k);
    candidate.pin_cells.push_back({static_cast<int16_t>(c[0]),
                                   static_cast<int16_t>(c[1]),
                                   static_cast<int16_t>(c[2])});
  }
  // Axis of a unit edge: 0 = x, 1 = y, 2 = via.
  auto axis = [&](uint64_t a, uint64_t b) {
    const auto pa = decode(a);
    const auto pb = decode(b);
    if (pa[0] != pb[0]) {
      return 2;
    }
    return pa[1] != pb[1] ? 0 : 1;
  };

  const uint64_t root = representative[0];
  candidate.nodes.clear();
  candidate.objective = 0.0;
  std::unordered_map<uint64_t, char> visited;
  {
    const auto p = decode(root);
    candidate.nodes.push_back({static_cast<int16_t>(p[0]),
                               static_cast<int16_t>(p[1]),
                               static_cast<int16_t>(p[2]),
                               -1});
  }
  visited[root] = 1;
  // (node index, vertex) whose unexplored neighbours are pending.
  std::vector<std::pair<int32_t, uint64_t>> stack{{0, root}};
  int num_vias = 0;
  while (!stack.empty()) {
    const auto [node_index, vertex] = stack.back();
    stack.pop_back();
    auto it = adjacency.find(vertex);
    if (it == adjacency.end()) {
      continue;
    }
    const std::vector<Link> links = it->second;
    for (const auto& [next, pin] : links) {
      if (visited.contains(next)) {
        continue;
      }
      if (pin >= 0) {
        // Pin link: the pin's own metal, no wire or via.
        visited[next] = 1;
        const auto pc = decode(next);
        candidate.nodes.push_back({static_cast<int16_t>(pc[0]),
                                   static_cast<int16_t>(pc[1]),
                                   static_cast<int16_t>(pc[2]),
                                   node_index,
                                   pin});
        stack.emplace_back(static_cast<int32_t>(candidate.nodes.size() - 1),
                           next);
        continue;
      }
      // Walk a straight run of wire or via unit edges.
      const int run_axis = axis(vertex, next);
      uint64_t prev = vertex;
      uint64_t cur = next;
      while (true) {
        visited[cur] = 1;
        const auto pp = decode(prev);
        if (run_axis == 2) {
          num_vias++;
        } else {
          const auto pc = decode(cur);
          candidate.objective += wireObjective(
              pp[0], std::min(pp[1], pc[1]), std::min(pp[2], pc[2]));
        }
        const auto& list = adjacency[cur];
        if (pinned.contains(cur) || list.size() != 2) {
          break;
        }
        const Link& after = list[0].first == prev ? list[1] : list[0];
        if (after.second >= 0 || visited.contains(after.first)
            || axis(cur, after.first) != run_axis) {
          break;
        }
        prev = cur;
        cur = after.first;
      }
      const auto pc = decode(cur);
      candidate.nodes.push_back({static_cast<int16_t>(pc[0]),
                                 static_cast<int16_t>(pc[1]),
                                 static_cast<int16_t>(pc[2]),
                                 node_index});
      stack.emplace_back(static_cast<int32_t>(candidate.nodes.size() - 1), cur);
    }
  }
  candidate.objective += params_.via_weight * num_vias;
  candidate.usage.clear();
  candidate.hash = 0;
  candidate.count = 0;
  if (!need_usage) {
    return true;
  }

  const std::shared_ptr<GRTreeNode> tree = buildTree(candidate);
  std::vector<std::pair<uint32_t, double>> usage;
  grid_graph_->forEachTreeDemand(
      tree,
      state.net->getNdrCosts(),
      [&](const int layer, const PointT edge, const CapacityT demand) {
        usage.emplace_back(edgeIndex(layer, edge.x(), edge.y()),
                           static_cast<double>(demand));
      });
  candidate.space.clear();
  if (params_.extra_space) {
    // Sec. 2.2: b* = (chi(T*), s) with the extra space s that minimizes
    // each edge's cost; the edge uses w(n, e) + s, power and yield their
    // gamma(s) (Fig. 1).
    std::vector<std::array<int, 3>> wires;
    forEachWireEdge(candidate,
                    [&](int l, int x, int y) { wires.push_back({l, x, y}); });
    std::ranges::sort(wires);
    wires.erase(std::ranges::unique(wires).begin(), wires.end());
    double power = 0.0;
    double yield = 0.0;
    for (const auto& [l, x, y] : wires) {
      double space = 0.0;
      spaceCost(l, x, y, &space);
      const uint32_t r = edgeIndex(l, x, y);
      if (space > 0.0) {
        usage.emplace_back(r, space);
        candidate.space.emplace_back(r, static_cast<float>(space));
      }
      const double length = edgeLength(l, x, y);
      power += length * space_model_[l].powerPerLength(space);
      yield += length * space_model_[l].yieldPerLength(space);
    }
    if (use_power_) {
      usage.emplace_back(power_resource_, power);
    }
    if (use_yield_) {
      usage.emplace_back(yield_resource_, yield);
    }
  }
  std::ranges::sort(usage);
  candidate.usage.clear();
  for (const auto& [r, d] : usage) {
    if (!candidate.usage.empty() && candidate.usage.back().first == r) {
      candidate.usage.back().second += d;
    } else {
      candidate.usage.emplace_back(r, d);
    }
  }
  uint64_t hash = 1469598103934665603ULL;
  for (const auto& [r, d] : candidate.usage) {
    hash = (hash ^ r) * 1099511628211ULL;
    hash = (hash ^ static_cast<uint64_t>(std::llround(d * 1048576.0)))
           * 1099511628211ULL;
  }
  candidate.hash = hash;
  candidate.count = 0;
  return true;
}

std::shared_ptr<GRTreeNode> BonnRouteGlobal::buildTree(
    const Candidate& candidate) const
{
  std::vector<std::shared_ptr<GRTreeNode>> nodes;
  nodes.reserve(candidate.nodes.size());
  for (const NodeRec& rec : candidate.nodes) {
    nodes.push_back(std::make_shared<GRTreeNode>(rec.layer, rec.x, rec.y));
    if (rec.parent >= 0) {
      nodes[rec.parent]->addChild(nodes.back(), rec.pin);
    }
  }
  return nodes.front();
}

double BonnRouteGlobal::edgeLength(const int layer,
                                   const int x,
                                   const int y) const
{
  const int direction = layer_direction_[layer];
  return edge_length_[direction][direction == MetalLayer::H ? x : y];
}

template <typename F>
void BonnRouteGlobal::forEachWireEdge(const Candidate& candidate, F&& fn) const
{
  for (const NodeRec& rec : candidate.nodes) {
    if (rec.parent < 0 || rec.pin >= 0) {
      continue;
    }
    const NodeRec& up = candidate.nodes[rec.parent];
    if (up.layer != rec.layer) {
      continue;
    }
    if (layer_direction_[rec.layer] == MetalLayer::H) {
      for (int c = std::min(rec.x, up.x); c < std::max(rec.x, up.x); c++) {
        fn(rec.layer, c, rec.y);
      }
    } else {
      for (int c = std::min(rec.y, up.y); c < std::max(rec.y, up.y); c++) {
        fn(rec.layer, rec.x, c);
      }
    }
  }
}

void BonnRouteGlobal::zeroSpaceResources(const Candidate& candidate,
                                         double& power,
                                         double& yield) const
{
  power = 0.0;
  yield = 0.0;
  forEachWireEdge(candidate, [&](int l, int x, int y) {
    const double length = edgeLength(l, x, y);
    power += length * space_model_[l].powerPerLength(0.0);
    yield += length * space_model_[l].yieldPerLength(0.0);
  });
}

// Fig. 1 resource models per wiring layer. Interpretation choices (see
// BonnRouteResources.h and src/grt/README.md):
//  - gap g = pitch - width (the gap to a wire on the neighbor track),
//    growth a = pitch / 2 (s tracks of space split over both sides);
//  - capacitance per unit length = odb's area capacitance * width + 2 *
//    edge capacitance (LEF CPERSQDIST / EDGECAPACITANCE, or the values
//    set_layer_rc stores, which put the whole wire capacitance in the area
//    term); coupling_fraction of it couples to the neighbor wires at gap g.
//    Activity 1: power is proportional to the switched capacitance;
//  - minimum defect size x0 = the smallest width or gap of the routing
//    layers, so x0 <= d on every layer.
// NDR nets use w(n, e) = their NDR factor for the space resource and the
// default-wire power / yield models.
void BonnRouteGlobal::initSpaceModels()
{
  num_edges_ = static_cast<uint32_t>(num_layers_) * x_size_ * y_size_;
  power_resource_ = num_edges_;
  yield_resource_ = num_edges_ + 1;
  use_power_ = params_.extra_space && params_.power_target > 0.0;
  use_yield_ = params_.extra_space && params_.yield_target > 0.0;
  space_grid_.step = params_.extra_space_step;
  space_grid_.max = params_.extra_space ? params_.extra_space_max : 0.0;
  const Design* design = grid_graph_->getDesign();
  double defect = std::numeric_limits<double>::max();
  for (int l = std::max(0, min_layer_); l < num_layers_; l++) {
    const MetalLayer& layer = design->getLayer(l);
    const double width = std::max(1, layer.getWidth());
    const double pitch = std::max(1, layer.getPitch());
    const double gap
        = pitch - width > 0.0 ? pitch - width : std::max(1, layer.getSpacing());
    defect = std::min({defect, width, gap});
  }
  if (defect == std::numeric_limits<double>::max()) {
    defect = 1.0;
  }
  space_model_.clear();
  for (int l = 0; l < num_layers_; l++) {
    const MetalLayer& layer = design->getLayer(l);
    odb::dbTechLayer* tech = layer.getTechLayer();
    space_model_.push_back(bonnroute::makeSpaceLayerModel(
        layer.getWidth(),
        layer.getPitch(),
        layer.getSpacing(),
        tech != nullptr ? tech->getCapacitance() : 0.0,
        tech != nullptr ? tech->getEdgeCapacitance() : 0.0,
        grid_graph_->getLibDBU(),
        params_.coupling_fraction,
        defect));
    if (params_.extra_space && params_.verbose && l >= min_layer_) {
      const bonnroute::SpaceLayerModel& m = space_model_.back();
      logger_->report(
          "BonnRoute extra space model {}: gap {:.0f}, growth {:.0f} DBU per "
          "track, ground {:.4g} / coupling {:.4g} fF per DBU, x0 {:.0f}.",
          grid_graph_->getLayerName(l),
          m.gap,
          m.growth,
          m.ground_capacitance,
          m.coupling_capacitance,
          m.defect_size);
    }
  }
}

double BonnRouteGlobal::spaceCost(const int layer,
                                  const int x,
                                  const int y,
                                  double* space) const
{
  const bonnroute::SpaceLayerModel& m = space_model_[layer];
  const double length = edgeLength(layer, x, y);
  const double p = use_power_ ? price(power_resource_) : 0.0;
  const double q = use_yield_ ? price(yield_resource_) : 0.0;
  const double shape = length * (p * m.powerShape() + q * m.yieldShape());
  const double constant
      = length
        * (p * m.ground_capacitance
           + q * m.defect_size * m.defect_size / (2.0 * m.wire_width));
  bonnroute::SpaceGrid grid = space_grid_;
  if (layer < min_layer_) {
    grid.max = 0.0;
  }
  const bonnroute::SpaceChoice choice = bonnroute::minimizeSpace(
      price(edgeIndex(layer, x, y)), shape, m.gap, m.growth, grid);
  if (space != nullptr) {
    *space = choice.space;
  }
  return choice.cost + constant;
}

// Sec. 2.1: "For the objective function, we do not have an upper bound a
// priori, but we can guess a value that we expect to be achievable and
// adapt it if needed." The power and yield bounds start at their targets
// times the congestion-free value without extra space. A bound grows when
// its resource's average load lambda_r exceeds the larger of 1 and the
// largest edge load by more than 5%: then the bound is not achievable
// alongside the edge capacities. The log price is rescaled to the new bound
// (ln y_r = eps * sum of g_n^r); prices of this resource may then fall, so
// the reuse certificates based on earlier oracle values are dropped.
std::string BonnRouteGlobal::resourceLoads(const int done) const
{
  if (!params_.extra_space) {
    return "";
  }
  const double scale = params_.epsilon * done;
  return fmt::format(", power load {:.3f}, yield load {:.3f}",
                     use_power_ ? log_price_[power_resource_] / scale : 0.0,
                     use_yield_ ? log_price_[yield_resource_] / scale : 0.0);
}

void BonnRouteGlobal::adaptResourceBounds(const int done)
{
  const double scale = params_.epsilon * done;
  if (scale <= 0.0) {
    return;
  }
  double edge = 0.0;
  for (uint32_t r = 0; r < num_edges_; r++) {
    edge = std::max(edge, log_price_[r]);
  }
  const double reference = std::max(1.0, edge / scale);
  bool changed = false;
  for (const uint32_t r : {power_resource_, yield_resource_}) {
    if ((r == power_resource_ && !use_power_)
        || (r == yield_resource_ && !use_yield_)) {
      continue;
    }
    const double lambda = log_price_[r] / scale;
    if (lambda <= 1.05 * reference) {
      continue;
    }
    const double factor = lambda / reference;
    norm_capacity_[r] *= factor;
    log_price_[r] /= factor;
    weight_[r] = std::exp(log_price_[r] - shift_) / norm_capacity_[r];
    bound_adaptations_++;
    changed = true;
    if (params_.verbose) {
      logger_->report(
          "BonnRoute phase {:3d}: {} bound raised by {:.3f} (load {:.3f}, "
          "edge load {:.3f}).",
          done,
          r == power_resource_ ? "power" : "yield",
          factor,
          lambda,
          edge / scale);
    }
  }
  if (changed) {
    for (NetState& state : states_) {
      state.last_log_cost = -1e300;
    }
  }
}

// After rounding and re-choosing, s is trimmed so that every edge holds it
// next to all nets' wires: the edge's extra space is scaled down to the
// capacity left after the committed wires (and rounded down to the grid).
// CUGR's demand keeps the wires only, so its maze clean-up may still use
// that room; trees it reroutes lose their extra space (CUGR exports only
// the space of trees left as committed here).
void BonnRouteGlobal::trimAndExportSpace()
{
  extra_space_.clear();
  std::vector<double> total(num_edges_, 0.0);
  for (const int index : order_) {
    const NetState& state = states_[index];
    if (state.chosen < 0) {
      continue;
    }
    for (const auto& [r, sp] : state.candidates[state.chosen].space) {
      total[r] += sp;
    }
  }
  std::vector<float> keep(num_edges_, 1.0f);
  double before = 0.0;
  for (uint32_t r = 0; r < num_edges_; r++) {
    if (total[r] <= 0.0) {
      continue;
    }
    before += total[r];
    const int y = static_cast<int>(r % y_size_);
    const int x = static_cast<int>((r / y_size_) % x_size_);
    const int l = static_cast<int>(r / y_size_ / x_size_);
    const GraphEdge edge = grid_graph_->getEdge(l, x, y);
    const double slack = std::max(0.0, edge.capacity - edge.demand);
    if (total[r] > slack) {
      keep[r] = static_cast<float>(slack / total[r]);
    }
  }
  const double step = space_grid_.step > 0.0 ? space_grid_.step : 1.0;
  double after = 0.0;
  double power = 0.0;
  double power0 = 0.0;
  double yield = 0.0;
  double yield0 = 0.0;
  int64_t edges = 0;
  for (const int index : order_) {
    const NetState& state = states_[index];
    if (state.chosen < 0) {
      continue;
    }
    const Candidate& chosen = state.candidates[state.chosen];
    NetSpace ns;
    ns.net = state.net;
    ns.tree = state.net->getRoutingTree();
    std::vector<std::pair<uint32_t, float>> trimmed;
    for (const auto& [r, sp] : chosen.space) {
      const double v = std::floor(sp * keep[r] / step + 1e-6) * step;
      if (v <= 0.0) {
        continue;
      }
      trimmed.emplace_back(r, static_cast<float>(v));
      const int y = static_cast<int>(r % y_size_);
      const int x = static_cast<int>((r / y_size_) % x_size_);
      const int l = static_cast<int>(r / y_size_ / x_size_);
      ns.edges.push_back({l + 1, x, y, static_cast<float>(v)});
      after += v;
      edges++;
    }
    double p0 = 0.0;
    double y0 = 0.0;
    zeroSpaceResources(chosen, p0, y0);
    power0 += p0;
    yield0 += y0;
    forEachWireEdge(chosen, [&](int l, int x, int y) {
      const uint32_t r = edgeIndex(l, x, y);
      const auto it = std::ranges::lower_bound(
          trimmed, r, {}, &std::pair<uint32_t, float>::first);
      const double sp
          = (it != trimmed.end() && it->first == r) ? it->second : 0.0;
      const double length = edgeLength(l, x, y);
      power += length * space_model_[l].powerPerLength(sp);
      yield += length * space_model_[l].yieldPerLength(sp);
    });
    if (!ns.edges.empty()) {
      extra_space_.push_back(std::move(ns));
    }
  }
  if (params_.verbose) {
    logger_->report(
        "BonnRoute extra space: {} edges of {} nets, {:.1f} tracks ({:.1f} "
        "before trimming to capacity); wire power {:.4g} fF ({:+.2f}% vs. "
        "no extra space), yield critical area {:.4g} ({:+.2f}%); bounds "
        "adapted {} times.",
        edges,
        extra_space_.size(),
        after,
        before,
        power,
        power0 > 0.0 ? 100.0 * (power - power0) / power0 : 0.0,
        yield,
        yield0 > 0.0 ? 100.0 * (yield - yield0) / yield0 : 0.0,
        bound_adaptations_);
  }
}

double BonnRouteGlobal::linearCost(
    const double objective,
    const std::vector<std::pair<uint32_t, double>>& usage) const
{
  double cost = objectiveWeight() * objective;
  for (const auto& [r, d] : usage) {
    cost += price(r) * d;
  }
  return cost;
}

double BonnRouteGlobal::linearCost(const Candidate& candidate) const
{
  return linearCost(candidate.objective, candidate.usage);
}

double BonnRouteGlobal::logCost(const Candidate& candidate) const
{
  return std::log(std::max(linearCost(candidate), 1e-300)) + shift_;
}

const std::vector<std::pair<uint32_t, double>>& BonnRouteGlobal::resultUsage(
    const NetState& state,
    const OracleResult& result) const
{
  return result.known >= 0 ? state.candidates[result.known].usage
                           : result.candidate.usage;
}

int BonnRouteGlobal::storeResult(NetState& state, OracleResult&& result)
{
  if (result.known >= 0) {
    return result.known;
  }
  const double objective = result.candidate.objective;
  const int index = insertCandidate(state, std::move(result.candidate));
  if (result.keyed) {
    result.key.candidate = index;
    result.key.objective = objective;
    state.tree_keys.push_back(std::move(result.key));
  }
  return index;
}

void BonnRouteGlobal::applyPrices(const Candidate& candidate)
{
  const double eps = params_.epsilon;
  for (const auto& [r, d] : candidate.usage) {
    const double lp = (log_price_[r] += eps * d / norm_capacity_[r]);
    max_log_price_ = std::max(max_log_price_, lp);
    if (touched_stamp_[r] != touched_epoch_) {
      touched_stamp_[r] = touched_epoch_;
      touched_.push_back(r);
    }
  }
  log_objective_price_ += eps * candidate.objective / objective_bound_;
  if (std::max(max_log_price_, log_objective_price_) - shift_
      > kRenormalizeAbove) {
    renormalize();
  }
}

// Within a batch nothing reads the weights, so each touched resource gets
// exp(ln y_r - shift_) / u_r once, from its final log price: the same value
// the last per-net update (or a renormalize after it) would have stored.
void BonnRouteGlobal::refreshWeights()
{
  for (const uint32_t r : touched_) {
    weight_[r] = std::exp(log_price_[r] - shift_) / norm_capacity_[r];
  }
  touched_.clear();
  if (++touched_epoch_ == 0) {
    std::ranges::fill(touched_stamp_, 0);
    touched_epoch_ = 1;
  }
  objective_weight_
      = std::exp(log_objective_price_ - shift_) / objective_bound_;
}

void BonnRouteGlobal::renormalize()
{
  shift_ = std::max(max_log_price_, log_objective_price_) - kRenormalizeTo;
  for (size_t r = 0; r < log_price_.size(); r++) {
    weight_[r] = std::exp(log_price_[r] - shift_) / norm_capacity_[r];
  }
  objective_weight_
      = std::exp(log_objective_price_ - shift_) / objective_bound_;
}

int BonnRouteGlobal::insertCandidate(NetState& state, Candidate&& candidate)
{
  for (size_t i = 0; i < state.candidates.size(); i++) {
    const Candidate& c = state.candidates[i];
    if (c.hash == candidate.hash && c.usage == candidate.usage) {
      return static_cast<int>(i);
    }
  }
  state.candidates.push_back(std::move(candidate));
  return static_cast<int>(state.candidates.size() - 1);
}

void BonnRouteGlobal::pilot(std::vector<double>* stack_counts)
{
  const int threads = std::max(1, params_.threads);
  std::vector<Scratch> scratch(threads);
  const bool check = params_.goal_check && params_.goal_orientation != 0;
  std::vector<Scratch> reference(check ? threads : 0);
  const int num_nets = static_cast<int>(order_.size());
  // Congestion-free objective: guess for the objective resource bound.
  std::vector<double> objective(num_nets, 0.0);
  const bool space = params_.extra_space;
  base_power_.assign(space ? num_nets : 0, 0.0);
  base_yield_.assign(space ? num_nets : 0, 0.0);
  // Via columns (x, y, low, high) of each net's congestion-free tree.
  std::vector<std::vector<std::array<int, 4>>> columns(
      stack_counts != nullptr ? num_nets : 0);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (int i = 0; i < num_nets; i++) {
    Scratch& s = scratch[omp_get_thread_num()];
    std::vector<Segment> segments;
    std::vector<Segment>* record = check ? &segments : nullptr;
    OracleResult r = oracle(states_[order_[i]], s, true, false, record);
    bool whole = false;
    if (!r.ok) {
      segments.clear();
      r = oracle(states_[order_[i]], s, true, true, record);
      whole = true;
    }
    if (check) {
      checkGoal(states_[order_[i]],
                reference[omp_get_thread_num()],
                true,
                whole,
                r.ok,
                segments);
    }
    objective[i] = r.ok ? r.candidate.objective : 0.0;
    if (space && r.ok) {
      zeroSpaceResources(r.candidate, base_power_[i], base_yield_[i]);
    }
    states_[order_[i]].log_base_objective
        = r.ok && r.candidate.objective > 0.0 ? std::log(r.candidate.objective)
                                              : -1e300;
    if (stack_counts != nullptr && r.ok) {
      const auto& nodes = r.candidate.nodes;
      for (const NodeRec& rec : nodes) {
        if (rec.parent < 0 || rec.pin >= 0) {
          continue;
        }
        const NodeRec& up = nodes[rec.parent];
        if (up.layer != rec.layer) {
          columns[i].push_back({rec.x,
                                rec.y,
                                std::min<int>(rec.layer, up.layer),
                                std::max<int>(rec.layer, up.layer)});
        }
      }
    }
  }
  base_objective_ = 0.0;
  for (const double o : objective) {
    base_objective_ += o;
  }
  for (const Scratch& s : scratch) {
    search_pops_ += s.pops;
    search_settled_ += s.settled;
    pilot_build_time_ += s.build_seconds;
  }
  if (stack_counts == nullptr) {
    return;
  }
  // Sec. 2.5: "A stacked via from layer l to layer l + 2 of course also
  // consumes space on layer l + 1." Merge each net's columns per cell like
  // GridGraph's via-column model and count the layers strictly inside.
  stack_counts->assign(static_cast<size_t>(num_layers_) * x_size_ * y_size_,
                       0.0);
  for (auto& spans : columns) {
    std::ranges::sort(spans);
    for (size_t i = 0; i < spans.size();) {
      const auto [x, y, low, first_high] = spans[i];
      int high = first_high;
      for (i++; i < spans.size() && spans[i][0] == x && spans[i][1] == y
                && spans[i][2] <= high;
           i++) {
        high = std::max(high, spans[i][3]);
      }
      for (int l = low + 1; l < high; l++) {
        (*stack_counts)[grid_graph_->edgeFlatIndex(l, x, y)] += 1.0;
      }
    }
  }
}

bonnroute::CapacityGrid capacityGridFromCugr(const GridGraph& grid_graph)
{
  const Design* design = grid_graph.getDesign();
  bonnroute::CapacityGrid grid;
  for (int d = 0; d < 2; d++) {
    for (int i = 0; i <= grid_graph.getSize(d); i++) {
      grid.gridlines[d].push_back(grid_graph.getGridline(d, i));
    }
  }
  const BoxT die = design->getDieRegion();
  grid.die = odb::Rect(die.lx(), die.ly(), die.hx(), die.hy());
  grid.linear_via_demand = grid_graph.isViaDemandEnabled();
  std::vector<std::vector<BoxT>> obstacles;
  design->getAllObstacles(obstacles, true);
  const int num_layers = design->getNumLayers();
  for (int l = 0; l < num_layers; l++) {
    const MetalLayer& metal = design->getLayer(l);
    bonnroute::CapacityLayer layer;
    layer.name = metal.getName();
    layer.direction = metal.getDirection() == MetalLayer::H
                          ? bonnroute::kHorizontal
                          : bonnroute::kVertical;
    for (int t = 0; t < metal.getNumTracks(); t++) {
      layer.tracks.push_back(metal.getTrackLocation(t));
    }
    layer.pitch = metal.getPitch();
    layer.width = metal.getWidth();
    layer.spacing_widths = metal.getParallelWidths();
    layer.spacing_lengths = metal.getParallelLengths();
    layer.spacing = metal.getParallelSpacingTable();
    layer.routable = l >= grid_graph.getMinRoutingLayer();
    layer.adjustment = metal.getAdjustment();
    if (l >= 1 && l + 1 < num_layers) {
      // A column through l puts the larger of the pads of the vias below
      // and above on it (GridGraph::stackPadLength).
      const double lower = design->getViaDemandLengthLower(l);
      const double upper = design->getViaDemandLengthUpper(l - 1);
      const bool use_lower = lower >= upper;
      layer.through_pad_length = std::max(lower, upper);
      layer.through_pad_along = use_lower ? design->getViaPadAlongLower(l)
                                          : design->getViaPadAlongUpper(l - 1);
      layer.through_pad_tracks = use_lower
                                     ? design->getViaPadTracksLower(l)
                                     : design->getViaPadTracksUpper(l - 1);
    }
    if (l < static_cast<int>(obstacles.size())) {
      for (const BoxT& box : obstacles[l]) {
        layer.obstacles.emplace_back(box.lx(), box.ly(), box.hx(), box.hy());
      }
    }
    grid.layers.push_back(std::move(layer));
  }
  return grid;
}

void BonnRouteGlobal::estimateCapacities()
{
  utl::Timer timer;
  cugr_capacity_.assign(static_cast<size_t>(num_layers_) * x_size_ * y_size_,
                        0.0);
  for (int l = 0; l < num_layers_; l++) {
    for (int x = 0; x < x_size_; x++) {
      for (int y = 0; y < y_size_; y++) {
        cugr_capacity_[grid_graph_->edgeFlatIndex(l, x, y)]
            = grid_graph_->getEdge(l, x, y).capacity;
      }
    }
  }
  bonnroute::CapacityEstimatorOptions options;
  options.blockage_extension = params_.capacity_blockage_extension;
  options.single_tile_nets = params_.capacity_single_tile;
  options.intra_tile = params_.capacity_intra_tile;
  options.stacked_vias = params_.capacity_stacked_vias;
  options.single_tile_layers = params_.capacity_single_tile_layers;
  options.single_tile_margin = params_.capacity_single_tile_margin;
  estimator_ = std::make_unique<bonnroute::CapacityEstimator>(
      capacityGridFromCugr(*grid_graph_), options);
  estimator_->buildTrackGraph();
  std::vector<bonnroute::EstimatorNet> nets;
  nets.reserve(states_.size());
  for (const NetState& state : states_) {
    if (state.trivial && state.trivial_low < 0) {
      continue;
    }
    const GRNet* net = state.net;
    bonnroute::EstimatorNet en;
    en.single_tile = state.trivial;
    en.tile_x = state.trivial_point.x();
    en.tile_y = state.trivial_point.y();
    en.min_layer = net->getLayerRange().min_layer;
    en.max_layer = std::min(net->getLayerRange().max_layer, num_layers_ - 1);
    const auto& preferred = net->getPreferredAccessPoints();
    for (int p = 0; p < net->getNumPins(); p++) {
      bonnroute::EstimatorPin pin;
      pin.points = grid_graph_->pinAccessLocations(net, p);
      if (auto it = preferred.find(p); it != preferred.end()) {
        pin.cell_x = it->second.point.x();
        pin.cell_y = it->second.point.y();
        if (pin.points.empty()) {
          pin.points.emplace_back(
              odb::Point(grid_graph_->getGridCenter(0, pin.cell_x),
                         grid_graph_->getGridCenter(1, pin.cell_y)),
              it->second.layers.low());
        }
      }
      en.pins.push_back(std::move(pin));
    }
    nets.push_back(std::move(en));
  }
  estimator_->addNets(nets);
  applyCapacities(estimator_->capacities());
  estimate_time_ = timer.elapsed();
}

void BonnRouteGlobal::applyCapacities(const std::vector<double>& capacities)
{
  for (int l = std::max(0, grid_graph_->getMinRoutingLayer()); l < num_layers_;
       l++) {
    for (int x = 0; x < x_size_; x++) {
      for (int y = 0; y < y_size_; y++) {
        grid_graph_->setEdgeCapacity(
            l, x, y, capacities[grid_graph_->edgeFlatIndex(l, x, y)]);
      }
    }
  }
}

void BonnRouteGlobal::reportCapacities() const
{
  if (!params_.verbose || !estimator_) {
    return;
  }
  const bonnroute::CapacityTotals& totals = estimator_->totals();
  std::vector<double> cugr_layer(num_layers_, 0.0);
  for (int l = 0; l < num_layers_; l++) {
    for (int x = 0; x < x_size_; x++) {
      for (int y = 0; y < y_size_; y++) {
        cugr_layer[l] += cugr_capacity_[grid_graph_->edgeFlatIndex(l, x, y)];
      }
    }
  }
  double cugr = 0.0;
  double estimate = 0.0;
  for (int l = std::max(0, grid_graph_->getMinRoutingLayer()); l < num_layers_;
       l++) {
    cugr += cugr_layer[l];
    estimate += totals.estimate[l];
    logger_->report(
        "BonnRoute capacity {}: CUGR {:.1f}, estimate {:.1f} (tracks {:.1f}, "
        "single-tile nets -{:.1f}, intra-tile -{:.1f}, stacked vias -{:.1f}, "
        "adjustment {:.2f}).",
        grid_graph_->getLayerName(l),
        cugr_layer[l],
        totals.estimate[l],
        totals.tracks[l],
        totals.single_tile[l],
        totals.intra_tile[l],
        totals.stacked_vias[l],
        grid_graph_->getDesign()->getLayer(l).getAdjustment());
  }
  logger_->report(
      "BonnRoute capacity total: CUGR {:.1f}, estimate {:.1f}; single-tile "
      "nets routed {}, unroutable {}, vertices {}; estimation {:.2f}s.",
      cugr,
      estimate,
      totals.single_tile_routed,
      totals.single_tile_failed,
      totals.single_tile_vertices,
      estimate_time_);
}

void BonnRouteGlobal::resourceSharing()
{
  const int threads = std::max(1, params_.threads);
  std::vector<Scratch> scratch(threads);
  const bool goal = params_.goal_orientation != 0;
  const bool check = params_.goal_check && goal;
  std::vector<Scratch> reference(check ? threads : 0);
  const int num_nets = static_cast<int>(order_.size());

  const double base_objective = base_objective_;
  objective_bound_
      = std::max(1.0, (1.0 + params_.objective_slack) * base_objective);

  // Start the objective price where one gcell of wire costs about as much as
  // one track of a typical edge, so early phases already weigh length.
  std::vector<double> caps;
  for (const double c : capacity_) {
    if (c > 0.0 && c < std::numeric_limits<double>::max()) {
      caps.push_back(c);
    }
  }
  double typical_capacity = 1.0;
  if (!caps.empty()) {
    std::nth_element(caps.begin(), caps.begin() + caps.size() / 2, caps.end());
    typical_capacity = std::max(1.0, caps[caps.size() / 2]);
  }
  double gcell_objective = 0.0;
  {
    const int mid = std::max(0, x_size_ / 2 - 1);
    gcell_objective
        = params_.wire_weight
          * std::max(1, grid_graph_->getEdgeLength(MetalLayer::H, mid))
          / m2_pitch_;
  }
  log_objective_price_ = std::log(
      objective_bound_ / (typical_capacity * std::max(gcell_objective, 1e-9)));
  const double initial_log_objective = log_objective_price_;
  shift_ = 0.0;
  max_log_price_ = 0.0;
  for (size_t r = 0; r < log_price_.size(); r++) {
    weight_[r] = 1.0 / norm_capacity_[r];
  }
  touched_.clear();
  touched_stamp_.assign(log_price_.size(), 0);
  touched_epoch_ = 1;
  objective_weight_ = std::exp(log_objective_price_) / objective_bound_;
  if (log_objective_price_ > kRenormalizeAbove) {
    renormalize();
  }
  if (params_.verbose) {
    logger_->report(
        "BonnRoute: {} nets, base objective {:.4g}, bound {:.4g}, typical "
        "capacity {:.2f}.",
        num_nets,
        base_objective,
        objective_bound_,
        typical_capacity);
  }

  enum class Kind : uint8_t
  {
    kReuse,
    kNew,
    kFail
  };
  const int batch = params_.batch_size;
  std::vector<Kind> kinds(batch);
  std::vector<OracleResult> results(batch);
  const double log_tolerance
      = std::log1p(std::max(0.0, params_.reuse_tolerance));
  if (params_.concurrent_prices) {
    resourceSharingConcurrent(log_tolerance, initial_log_objective);
    return;
  }
  // Wall seconds per thread: busy time of the current batch, and the phase's
  // split into reuse tests, search, buildCandidate and cost evaluation.
  struct alignas(kCacheLine) ThreadTime
  {
    double batch = 0.0;
    double reuse = 0.0;
    double search = 0.0;
    double build = 0.0;
    double cost = 0.0;
  };
  // Wall seconds of one phase (or the sum over all phases).
  struct PhaseTime
  {
    double wall = 0.0;
    double goal = 0.0;
    double batches = 0.0;
    double longest = 0.0;  // sum over batches of the busiest thread
    double busy = 0.0;     // sum over batches and threads
    double reuse = 0.0;
    double search = 0.0;
    double build = 0.0;
    double cost = 0.0;
    double apply = 0.0;
    double weights = 0.0;
    void add(const PhaseTime& o)
    {
      wall += o.wall;
      goal += o.goal;
      batches += o.batches;
      longest += o.longest;
      busy += o.busy;
      reuse += o.reuse;
      search += o.search;
      build += o.build;
      cost += o.cost;
      apply += o.apply;
      weights += o.weights;
    }
  };
  std::vector<ThreadTime> thread_time(threads);
  PhaseTime total_time;
  int64_t num_batches = 0;
  double report_time = 0.0;
  for (int phase = 0; phase < params_.phases; phase++) {
    const SteadyClock::time_point phase_start = SteadyClock::now();
    PhaseTime phase_time;
    int64_t phase_calls = 0;
    int64_t phase_reuses = 0;
    if (goal) {
      updateGoalBounds(phase % params_.landmark_refresh == 0);
    }
    phase_time.goal = secondsBetween(phase_start, SteadyClock::now());
    for (int begin = 0; begin < num_nets; begin += batch) {
      const int end = std::min(num_nets, begin + batch);
      // Prices are frozen during the parallel part of the batch.
      const double batch_log_objective = log_objective_price_;
      const SteadyClock::time_point batch_start = SteadyClock::now();
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
      for (int i = begin; i < end; i++) {
        const SteadyClock::time_point net_start = SteadyClock::now();
        ThreadTime& own_time = thread_time[omp_get_thread_num()];
        const NetState& state = states_[order_[i]];
        const int slot = i - begin;
        results[slot].ok = false;
        if (state.last >= 0) {
          const double now = logCost(state.candidates[state.last]);
          // Relative to the objective price, which every phase raises.
          // Reuse certificate. Prices never decrease, so for the current
          // prices y: opt_n(y) >= opt_n(y_old) >= C_old / sigma, and
          // opt_n(y) >= y_obj * o* / U >= y_obj * o_base / (sigma U) with
          // o_base the phase-0 congestion-free oracle objective. Reusing b
          // when y.g(b) <= (1 + tol) * max(C_old, y_obj o_base / U) keeps
          // it a (1 + tol) sigma approximation, the oracle's own guarantee.
          const double log_objective_bound = batch_log_objective
                                             - std::log(objective_bound_)
                                             + state.log_base_objective;
          const double bound
              = std::max(state.last_log_cost, log_objective_bound);
          if (now - bound <= log_tolerance) {
            kinds[slot] = Kind::kReuse;
            const double t = secondsBetween(net_start, SteadyClock::now());
            own_time.reuse += t;
            own_time.batch += t;
            continue;
          }
        }
        Scratch& s = scratch[omp_get_thread_num()];
        const SteadyClock::time_point oracle_start = SteadyClock::now();
        const double build_before = s.build_seconds;
        std::vector<Segment> segments;
        std::vector<Segment>* record = check ? &segments : nullptr;
        OracleResult r = oracle(state, s, false, false, record);
        bool whole = false;
        if (!r.ok) {
          segments.clear();
          r = oracle(state, s, false, true, record);
          whole = true;
        }
        const SteadyClock::time_point oracle_end = SteadyClock::now();
        if (check) {
          checkGoal(state,
                    reference[omp_get_thread_num()],
                    false,
                    whole,
                    r.ok,
                    segments);
        }
        if (r.ok) {
          const auto& usage = resultUsage(state, r);
          r.log_cost
              = std::log(std::max(linearCost(r.candidate.objective, usage),
                                  1e-300))
                + shift_;
          // Differential check of the oracle's edge costs against the
          // committed demand: exact for a single path between single cells.
          if (state.terminals.size() == 2
              && state.terminals[0].cells.size() == 1
              && state.terminals[1].cells.size() == 1) {
            const double exact = linearCost(r.candidate.objective, usage);
            r.checked = true;
            r.mismatch = std::abs(exact - r.search_cost)
                         > 1e-9 * std::max(1.0, std::abs(exact));
            r.exact_cost = exact;
          }
          kinds[slot] = Kind::kNew;
        } else {
          kinds[slot] = state.last >= 0 ? Kind::kReuse : Kind::kFail;
        }
        results[slot] = std::move(r);
        const SteadyClock::time_point net_end = SteadyClock::now();
        const double build = s.build_seconds - build_before;
        own_time.reuse += secondsBetween(net_start, oracle_start);
        own_time.search += secondsBetween(oracle_start, oracle_end) - build;
        own_time.build += build;
        own_time.cost += secondsBetween(oracle_end, net_end);
        own_time.batch += secondsBetween(net_start, net_end);
      }
      const SteadyClock::time_point apply_start = SteadyClock::now();
      phase_time.batches += secondsBetween(batch_start, apply_start);
      double longest = 0.0;
      for (ThreadTime& t : thread_time) {
        longest = std::max(longest, t.batch);
        phase_time.busy += t.batch;
        t.batch = 0.0;
      }
      phase_time.longest += longest;
      num_batches++;
      // Apply the batch in net order: deterministic for any thread count.
      for (int i = begin; i < end; i++) {
        NetState& state = states_[order_[i]];
        const int slot = i - begin;
        if (kinds[slot] == Kind::kFail) {
          state.oracle_failures++;
          continue;
        }
        if (kinds[slot] == Kind::kNew) {
          phase_calls++;
          cost_checks_ += results[slot].checked;
          cost_mismatches_ += results[slot].mismatch;
          if (results[slot].mismatch && cost_mismatches_ <= 5) {
            const Candidate& c = results[slot].known >= 0
                                     ? state.candidates[results[slot].known]
                                     : results[slot].candidate;
            std::string nodes;
            for (const NodeRec& n : c.nodes) {
              nodes
                  += fmt::format(" ({},{},{})<{}", n.layer, n.x, n.y, n.parent);
            }
            logger_->report(
                "BonnRoute cost check: net {} search {:.6g} exact {:.6g} "
                "objective {:.6g} objective weight {:.6g}; tree{}",
                state.net->getName(),
                results[slot].search_cost,
                results[slot].exact_cost,
                c.objective,
                objective_weight_,
                nodes);
          }
          const double log_cost = results[slot].log_cost;
          state.last = storeResult(state, std::move(results[slot]));
          state.last_log_cost = log_cost;
        } else {
          reuses_++;
          phase_reuses++;
        }
        Candidate& chosen = state.candidates[state.last];
        chosen.count++;
        applyPrices(chosen);
      }
      const SteadyClock::time_point weights_start = SteadyClock::now();
      phase_time.apply += secondsBetween(apply_start, weights_start);
      refreshWeights();
      phase_time.weights += secondsBetween(weights_start, SteadyClock::now());
    }
    for (ThreadTime& t : thread_time) {
      phase_time.reuse += t.reuse;
      phase_time.search += t.search;
      phase_time.build += t.build;
      phase_time.cost += t.cost;
      t = ThreadTime{};
    }
    phase_time.wall = secondsBetween(phase_start, SteadyClock::now());
    total_time.add(phase_time);
    oracle_calls_ += phase_calls;
    const int done = phase + 1;
    if (done == params_.phases) {
      for (const Scratch& sc : scratch) {
        search_pops_ += sc.pops;
        search_settled_ += sc.settled;
      }
    }
    const SteadyClock::time_point report_start = SteadyClock::now();
    if (params_.verbose) {
      logger_->report(
          "BonnRoute phase timing {:3d}: {:.4f}s, oracle calls {}, reuses {}, "
          "max edge load {:.3f}; goal bounds {:.4f}s, batches {:.4f}s "
          "(longest thread {:.4f}s, busy {:.4f}s: reuse {:.4f}s, search "
          "{:.4f}s, build {:.4f}s, cost {:.4f}s), apply {:.4f}s, weights "
          "{:.4f}s.",
          done,
          phase_time.wall,
          phase_calls,
          phase_reuses,
          maxEdgeLoad(done),
          phase_time.goal,
          phase_time.batches,
          phase_time.longest,
          phase_time.busy,
          phase_time.reuse,
          phase_time.search,
          phase_time.build,
          phase_time.cost,
          phase_time.apply,
          phase_time.weights);
    }
    if (params_.verbose
        && (done == 1 || done % 25 == 0 || done == params_.phases)) {
      const double lambda_max = maxEdgeLoad(done);
      const double lambda_obj = (log_objective_price_ - initial_log_objective)
                                / (params_.epsilon * done);
      logger_->report(
          "BonnRoute phase {:3d}: oracle calls {:7d}, max edge load {:.3f}, "
          "objective load {:.3f}{}.",
          done,
          phase_calls,
          lambda_max,
          lambda_obj,
          resourceLoads(done));
    }
    report_time += secondsBetween(report_start, SteadyClock::now());
    if (params_.extra_space && params_.resource_adapt
        && (done == 1 || done % 25 == 0) && done < params_.phases) {
      adaptResourceBounds(done);
    }
  }
  if (params_.verbose) {
    logger_->report(
        "BonnRoute sharing timing: {} phases {:.3f}s; goal bounds {:.3f}s, "
        "batches {:.3f}s in {} batches (longest thread {:.3f}s, busy {:.3f}s "
        "on {} threads: reuse {:.3f}s, search {:.3f}s, build {:.3f}s, cost "
        "{:.3f}s), apply {:.3f}s, weights {:.3f}s, phase reports {:.3f}s.",
        params_.phases,
        total_time.wall,
        total_time.goal,
        total_time.batches,
        num_batches,
        total_time.longest,
        total_time.busy,
        threads,
        total_time.reuse,
        total_time.search,
        total_time.build,
        total_time.cost,
        total_time.apply,
        total_time.weights,
        report_time);
  }
}

double BonnRouteGlobal::maxEdgeLoad(const int done) const
{
  // ln y_r = eps * (sum of normalized loads); divide by phases for the
  // average fractional load lambda_r.
  double lambda_max = 0.0;
  for (uint32_t r = 0; r < num_edges_; r++) {
    lambda_max = std::max(lambda_max, log_price_[r]);
  }
  return lambda_max / (params_.epsilon * done);
}

void BonnRouteGlobal::applyPricesConcurrent(const Candidate& candidate)
{
  const double eps = params_.epsilon;
  for (const auto& [r, d] : candidate.usage) {
    const double lp = addRelaxed(log_price_[r], eps * d / norm_capacity_[r]);
    // Stores from threads that added earlier may land later: keep the
    // largest, which belongs to the latest log price (weights stay monotone).
    storeMax(weight_[r], std::exp(lp - shift_) / norm_capacity_[r]);
    storeMax(max_log_price_, lp);
  }
  const double lo = addRelaxed(log_objective_price_,
                               eps * candidate.objective / objective_bound_);
  storeMax(objective_weight_, std::exp(lo - shift_) / objective_bound_);
}

// Sec. 2.3: "the algorithm can be parallelized, and a highly scalable
// parallelized version has been implemented"; Sec. 5.1: "volatility-
// tolerant block solvers [Mueller et al. 2011] can be used to guarantee the
// desired approximation ratio even with concurrent access to the same
// resources by multiple threads". Each thread takes the next net, reads the
// current prices, calls the oracle and applies the price update at once.
// Prices only grow while a net is routed, which the reuse certificate and
// the goal-orientation bounds already assume. Nets of a phase are taken
// largest search region first so the phase ends without a long tail.
void BonnRouteGlobal::resourceSharingConcurrent(
    const double log_tolerance,
    const double initial_log_objective)
{
  const int threads = std::max(1, params_.threads);
  std::vector<Scratch> scratch(threads);
  const bool goal = params_.goal_orientation != 0;
  const int num_nets = static_cast<int>(order_.size());
  std::vector<int> order = order_;
  auto work = [&](int index) {
    const BoxT& r = states_[index].region;
    return static_cast<int64_t>(r.hx() - r.lx() + 1) * (r.hy() - r.ly() + 1);
  };
  std::ranges::stable_sort(order,
                           [&](int a, int b) { return work(a) > work(b); });
  std::atomic<int64_t> reuses{0};
  for (int phase = 0; phase < params_.phases; phase++) {
    const SteadyClock::time_point phase_start = SteadyClock::now();
    const int64_t reuses_before = reuses.load();
    if (goal) {
      updateGoalBounds(phase % params_.landmark_refresh == 0);
    }
    std::atomic<int64_t> phase_calls{0};
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (int i = 0; i < num_nets; i++) {
      NetState& state = states_[order[i]];
      bool reuse = false;
      if (state.last >= 0) {
        const double now = logCost(state.candidates[state.last]);
        const double log_objective_bound = loadRelaxed(log_objective_price_)
                                           - std::log(objective_bound_)
                                           + state.log_base_objective;
        const double bound = std::max(state.last_log_cost, log_objective_bound);
        reuse = now - bound <= log_tolerance;
      }
      if (!reuse) {
        Scratch& s = scratch[omp_get_thread_num()];
        OracleResult r = oracle(state, s, false, false);
        if (!r.ok) {
          r = oracle(state, s, false, true);
        }
        if (r.ok) {
          // No exact cost check: other threads change the prices between
          // the search and any recomputation.
          const double log_cost
              = std::log(std::max(
                    linearCost(r.candidate.objective, resultUsage(state, r)),
                    1e-300))
                + shift_;
          state.last = storeResult(state, std::move(r));
          state.last_log_cost = log_cost;
          phase_calls++;
        } else if (state.last < 0) {
          state.oracle_failures++;
          continue;
        } else {
          reuse = true;
        }
      }
      if (reuse) {
        reuses++;
      }
      Candidate& chosen = state.candidates[state.last];
      chosen.count++;
      applyPricesConcurrent(chosen);
    }
    oracle_calls_ += phase_calls;
    // Between phases no thread runs: renormalize here only. Within one
    // phase a log price grows by its normalized load, far below the
    // remaining headroom of exp().
    if (std::max(max_log_price_, log_objective_price_) - shift_
        > kRenormalizeAbove) {
      renormalize();
    }
    const int done = phase + 1;
    if (params_.verbose) {
      logger_->report(
          "BonnRoute phase timing {:3d}: {:.4f}s, oracle calls {}, reuses {}, "
          "max edge load {:.3f} (concurrent prices).",
          done,
          secondsBetween(phase_start, SteadyClock::now()),
          phase_calls.load(),
          reuses.load() - reuses_before,
          maxEdgeLoad(done));
    }
    if (params_.verbose
        && (done == 1 || done % 25 == 0 || done == params_.phases)) {
      const double lambda_max = maxEdgeLoad(done);
      const double lambda_obj = (log_objective_price_ - initial_log_objective)
                                / (params_.epsilon * done);
      logger_->report(
          "BonnRoute phase {:3d}: oracle calls {:7d}, max edge load {:.3f}, "
          "objective load {:.3f}{} (concurrent prices).",
          done,
          phase_calls.load(),
          lambda_max,
          lambda_obj,
          resourceLoads(done));
    }
    if (params_.extra_space && params_.resource_adapt
        && (done == 1 || done % 25 == 0) && done < params_.phases) {
      adaptResourceBounds(done);
    }
  }
  reuses_ += reuses;
  for (const Scratch& sc : scratch) {
    search_pops_ += sc.pops;
    search_settled_ += sc.settled;
  }
}

void BonnRouteGlobal::roundAndRepair()
{
  std::vector<double> demand(capacity_.size(), 0.0);
  std::mt19937_64 rng(params_.seed);
  for (const int index : order_) {
    NetState& state = states_[index];
    int total = 0;
    for (const Candidate& c : state.candidates) {
      total += c.count;
    }
    if (total == 0) {
      state.chosen = -1;
      continue;
    }
    std::uniform_int_distribution<int> pick(0, total - 1);
    int draw = pick(rng);
    for (size_t i = 0; i < state.candidates.size(); i++) {
      draw -= state.candidates[i].count;
      if (draw < 0) {
        state.chosen = static_cast<int>(i);
        break;
      }
    }
    for (const auto& [r, d] : state.candidates[state.chosen].usage) {
      demand[r] += d;
    }
  }
  auto overflow
      = [&](uint32_t r, double d) { return std::max(0.0, d - capacity_[r]); };
  auto report = [&](const char* tag) {
    double total = 0.0;
    int edges = 0;
    for (size_t r = 0; r < demand.size(); r++) {
      const double o = overflow(r, demand[r]);
      if (o > 1e-9) {
        total += o;
        edges++;
      }
    }
    if (params_.verbose) {
      logger_->report(
          "BonnRoute {}: overflow {:.2f} on {} edges.", tag, total, edges);
    }
  };
  report("after rounding");

  int changes = 0;
  for (int pass = 0; pass < params_.repair_passes; pass++) {
    int pass_changes = 0;
    for (const int index : order_) {
      NetState& state = states_[index];
      if (state.chosen < 0 || state.candidates.size() < 2) {
        continue;
      }
      const auto& current = state.candidates[state.chosen].usage;
      bool on_overflow = false;
      for (const auto& [r, d] : current) {
        if (demand[r] > capacity_[r] + 1e-9) {
          on_overflow = true;
          break;
        }
      }
      if (!on_overflow) {
        continue;
      }
      // Overflow change of replacing the current candidate by `alt`.
      auto delta = [&](const std::vector<std::pair<uint32_t, double>>& alt) {
        double change = 0.0;
        size_t i = 0;
        size_t j = 0;
        while (i < current.size() || j < alt.size()) {
          uint32_t r;
          double before_d = 0.0;
          double after_d = 0.0;
          if (j == alt.size()
              || (i < current.size() && current[i].first < alt[j].first)) {
            r = current[i].first;
            before_d = current[i++].second;
          } else if (i == current.size() || alt[j].first < current[i].first) {
            r = alt[j].first;
            after_d = alt[j++].second;
          } else {
            r = current[i].first;
            before_d = current[i++].second;
            after_d = alt[j++].second;
          }
          const double base = demand[r] - before_d;
          change += overflow(r, base + after_d) - overflow(r, demand[r]);
        }
        return change;
      };
      int best = state.chosen;
      double best_delta = 0.0;
      double best_objective = state.candidates[state.chosen].objective;
      for (size_t i = 0; i < state.candidates.size(); i++) {
        const Candidate& c = state.candidates[i];
        if (static_cast<int>(i) == state.chosen || c.count == 0) {
          continue;
        }
        const double d = delta(c.usage);
        if (d < best_delta - 1e-9
            || (d < -1e-9 && std::abs(d - best_delta) <= 1e-9
                && c.objective < best_objective)) {
          best = static_cast<int>(i);
          best_delta = d;
          best_objective = c.objective;
        }
      }
      if (best != state.chosen) {
        for (const auto& [r, d] : current) {
          demand[r] -= d;
        }
        state.chosen = best;
        for (const auto& [r, d] : state.candidates[best].usage) {
          demand[r] += d;
        }
        pass_changes++;
      }
    }
    changes += pass_changes;
    if (pass_changes == 0) {
      break;
    }
  }
  if (params_.verbose) {
    logger_->report("BonnRoute re-choose: {} candidate changes.", changes);
  }
  report("after re-choose");
}

void BonnRouteGlobal::commit()
{
  int failed = 0;
  for (const int index : order_) {
    NetState& state = states_[index];
    if (state.chosen < 0) {
      failed++;
      state.net->clearRoutingTree();
      continue;
    }
    const Candidate& chosen = state.candidates[state.chosen];
    std::shared_ptr<GRTreeNode> tree = buildTree(chosen);
    state.net->setRoutingTree(tree);
    // The oracle chose the access cell of every pin.
    for (size_t t = 0; t < state.terminals.size(); t++) {
      const auto& c = chosen.pin_cells[t];
      state.net->addPreferredAccessPoint(
          state.terminals[t].pin_index,
          AccessPoint{PointT(c[1], c[2]), IntervalT(c[0])});
    }
    grid_graph_->addTreeUsage(tree, state.net->getNdrCosts());
  }
  if (params_.extra_space) {
    trimAndExportSpace();
  }
  if (failed > 0) {
    logger_->warn(GRT,
                  760,
                  "BonnRoute found no route for {} nets; they fall back to "
                  "pattern routing.",
                  failed);
  }
}

void BonnRouteGlobal::run()
{
  utl::Timer timer;
  // Wall seconds of the stages, reported with verbose.
  SteadyClock::time_point mark = SteadyClock::now();
  auto lap = [&mark]() {
    const SteadyClock::time_point now = SteadyClock::now();
    const double seconds = secondsBetween(mark, now);
    mark = now;
    return seconds;
  };
  prepareNets();
  const double prepare_time = lap();
  min_layer_ = num_layers_;
  for (const GRNet* net : nets_) {
    min_layer_ = std::min(min_layer_, net->getLayerRange().min_layer);
  }
  const bool estimate = params_.capacity_model == 1;
  if (estimate) {
    estimateCapacities();
  }
  const double capacity_time = lap();
  // Single-gcell nets keep their pin stack; their usage is background.
  // With the Sec. 2.5 estimate, their routes inside the tile are already
  // part of the capacities, so they keep only the pin's cell.
  const bool single_tile_in_capacity = estimate && params_.capacity_single_tile;
  for (NetState& state : states_) {
    if (!state.trivial || state.trivial_low < 0) {
      continue;
    }
    const PointT& pt = state.trivial_point;
    auto node = std::make_shared<GRTreeNode>(state.trivial_low, pt.x(), pt.y());
    if (state.trivial_top != state.trivial_low && !single_tile_in_capacity) {
      node->addChild(
          std::make_shared<GRTreeNode>(state.trivial_top, pt.x(), pt.y()));
    }
    state.net->setRoutingTree(node);
    grid_graph_->addTreeUsage(node, state.net->getNdrCosts());
  }
  ndr_min_.assign(num_layers_, std::numeric_limits<double>::infinity());
  for (const int index : order_) {
    for (int l = 0; l < num_layers_; l++) {
      ndr_min_[l] = std::min(ndr_min_[l], states_[index].net->getNdrCost(l));
    }
  }
  for (double& f : ndr_min_) {
    if (!std::isfinite(f)) {
      f = 1.0;
    }
  }
  initSpaceModels();
  const double trivial_time = lap();
  std::vector<double> stack_counts;
  const bool count_stacks = estimate && params_.capacity_stacked_vias;
  pilot(count_stacks ? &stack_counts : nullptr);
  const double pilot_time = lap();
  if (estimate) {
    utl::Timer stack_timer;
    if (count_stacks) {
      estimator_->addStackedVias(stack_counts);
      applyCapacities(estimator_->capacities());
    }
    estimate_time_ += stack_timer.elapsed();
    reportCapacities();
  }
  const double stacked_vias_time = lap();
  // Edge resources, then power and yield (used with extra_space only).
  const size_t num_resources = static_cast<size_t>(num_edges_) + 2;
  capacity_.assign(num_resources, 0.0);
  norm_capacity_.assign(num_resources, params_.min_capacity);
  log_price_.assign(num_resources, 0.0);
  weight_.assign(num_resources, 0.0);
  for (int l = 0; l < num_layers_; l++) {
    for (int x = 0; x < x_size_; x++) {
      for (int y = 0; y < y_size_; y++) {
        const GraphEdge edge = grid_graph_->getEdge(l, x, y);
        const uint32_t r = edgeIndex(l, x, y);
        if (l < min_layer_) {
          // Pin layers: via pads there are not a routing resource.
          capacity_[r] = std::numeric_limits<double>::max();
          norm_capacity_[r] = 1e30;
          continue;
        }
        capacity_[r] = edge.capacity - edge.demand;
        norm_capacity_[r] = std::max(capacity_[r], params_.min_capacity);
      }
    }
  }
  // Power and yield never count as overflow in rounding (Sec. 2.4 repairs
  // edge capacities); their bounds u_r normalize the prices.
  for (const uint32_t r : {power_resource_, yield_resource_}) {
    capacity_[r] = std::numeric_limits<double>::max();
    norm_capacity_[r] = 1e30;
  }
  if (params_.extra_space) {
    double power = 0.0;
    double yield = 0.0;
    for (size_t i = 0; i < base_power_.size(); i++) {
      power += base_power_[i];
      yield += base_yield_[i];
    }
    use_power_ = use_power_ && power > 0.0;
    use_yield_ = use_yield_ && yield > 0.0;
    if (use_power_) {
      norm_capacity_[power_resource_] = params_.power_target * power;
    }
    if (use_yield_) {
      norm_capacity_[yield_resource_] = params_.yield_target * yield;
    }
    if (params_.verbose) {
      logger_->report(
          "BonnRoute extra space: step {:.2f}, max {:.2f} tracks; power "
          "bound {:.4g} fF ({:.2f} x {:.4g}), yield bound {:.4g} ({:.2f} x "
          "{:.4g}).",
          space_grid_.step,
          space_grid_.max,
          use_power_ ? norm_capacity_[power_resource_] : 0.0,
          params_.power_target,
          power,
          use_yield_ ? norm_capacity_[yield_resource_] : 0.0,
          params_.yield_target,
          yield);
    }
  }
  const double resources_time = lap();
  resourceSharing();
  const double sharing_time = timer.elapsed();
  const double sharing_proper_time = lap();
  roundAndRepair();
  const double rounding_time = lap();
  commit();
  const double commit_time = lap();
  if (params_.verbose) {
    logger_->report(
        "BonnRoute stage timing: prepare nets {:.3f}s, capacity estimate "
        "{:.3f}s, trivial nets {:.3f}s, pilot {:.3f}s (build {:.3f}s thread "
        "time), stacked vias {:.3f}s, resources {:.3f}s, resource sharing "
        "{:.3f}s, rounding {:.3f}s, commit {:.3f}s.",
        prepare_time,
        capacity_time,
        trivial_time,
        pilot_time,
        pilot_build_time_,
        stacked_vias_time,
        resources_time,
        sharing_proper_time,
        rounding_time,
        commit_time);
    size_t candidates = 0;
    for (const int index : order_) {
      candidates += states_[index].candidates.size();
    }
    logger_->report(
        "BonnRoute: oracle cost check {} of {} two-pin searches differ "
        "from the committed demand cost.",
        cost_mismatches_,
        cost_checks_);
    logger_->report(
        "BonnRoute search: goal orientation {}, {} heap pops, {} settled; "
        "{} landmarks, {} landmark updates {:.2f}s.",
        params_.goal_orientation,
        search_pops_,
        search_settled_,
        params_.goal_orientation != 0 ? params_.landmarks : 0,
        landmark_updates_,
        landmark_time_);
    if (params_.goal_check && params_.goal_orientation != 0) {
      logger_->report(
          "BonnRoute goal check: {} oracle calls, {} connections compared, "
          "{} equal-cost ties, {} mismatches.",
          goal_checked_.load(),
          goal_segments_.load(),
          goal_ties_.load(),
          goal_mismatches_.load());
    }
    logger_->report(
        "BonnRoute: {} oracle calls, {} reuses, {:.2f} candidates/net, "
        "resource sharing {:.2f}s, total {:.2f}s.",
        oracle_calls_,
        reuses_,
        order_.empty() ? 0.0 : static_cast<double>(candidates) / order_.size(),
        sharing_time,
        timer.elapsed());
  }
}

}  // namespace grt
