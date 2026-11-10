// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "BonnRoute.h"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "GRNet.h"
#include "GRTree.h"
#include "GridGraph.h"
#include "Layers.h"
#include "geo.h"
#include "utl/Logger.h"
#include "utl/timer.h"

namespace grt {

using utl::GRT;

namespace {
// Keep exp(ln y - shift) inside the double range.
constexpr double kRenormalizeAbove = 600.0;
constexpr double kRenormalizeTo = 300.0;
}  // namespace

bool BonnRouteParams::set(const std::string& name, const double value)
{
  if (name == "phases") {
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
  } else if (name == "maze_cleanup") {
    maze_cleanup = value != 0.0;
  } else if (name == "verbose") {
    verbose = value != 0.0;
  } else {
    return false;
  }
  return true;
}

// Per-thread Dijkstra state over the search region, reset by epoch stamps.
struct BonnRouteGlobal::Scratch
{
  std::vector<double> dist;
  std::vector<int32_t> parent;
  std::vector<uint32_t> stamp;
  std::vector<uint32_t> tree_stamp;
  std::vector<uint32_t> term_stamp;
  std::vector<int32_t> term;
  std::vector<std::pair<double, int64_t>> heap;
  uint32_t epoch = 0;

  void prepare(size_t n)
  {
    if (dist.size() < n) {
      dist.resize(n);
      parent.resize(n);
      stamp.resize(n, 0);
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
}

BonnRouteGlobal::~BonnRouteGlobal() = default;

double BonnRouteGlobal::wireObjective(const int layer,
                                      const int x,
                                      const int y) const
{
  const int direction = grid_graph_->getLayerDirection(layer);
  const int index = direction == MetalLayer::H ? x : y;
  return params_.wire_weight * grid_graph_->getEdgeLength(direction, index)
         / m2_pitch_;
}

double BonnRouteGlobal::wireCost(const int layer,
                                 const int x,
                                 const int y,
                                 const double factor) const
{
  return objective_weight_ * wireObjective(layer, x, y)
         + weight_[edgeIndex(layer, x, y)] * factor;
}

double BonnRouteGlobal::viaCost(const int layer,
                                const int x,
                                const int y,
                                const GRNet* net) const
{
  double cost = objective_weight_ * params_.via_weight;
  if (!via_demand_) {
    return cost;
  }
  for (int side = 0; side < 2; side++) {
    const int pad_layer = layer + side;
    const double pad
        = side == 0 ? via_pad_lower_[layer] : via_pad_upper_[layer];
    if (pad <= 0.0) {
      continue;
    }
    const int direction = grid_graph_->getLayerDirection(pad_layer);
    const int c = direction == MetalLayer::H ? x : y;
    const int size = grid_graph_->getSize(direction);
    const int lower_len
        = c > 0 ? grid_graph_->getEdgeLength(direction, c - 1) : 0;
    const int higher_len
        = c < size - 1 ? grid_graph_->getEdgeLength(direction, c) : 0;
    const int sum = lower_len + higher_len;
    if (sum == 0) {
      continue;
    }
    const double demand = pad / sum * net->getNdrCost(pad_layer);
    if (lower_len > 0) {
      const int lx = direction == MetalLayer::H ? x - 1 : x;
      const int ly = direction == MetalLayer::H ? y : y - 1;
      cost += weight_[edgeIndex(pad_layer, lx, ly)] * demand;
    }
    if (higher_len > 0) {
      cost += weight_[edgeIndex(pad_layer, x, y)] * demand;
    }
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
    const AccessPointMap selected = grid_graph_->selectAccessPoints(net);
    std::vector<std::pair<PointT, IntervalT>> points;
    points.reserve(selected.size());
    for (const auto& [point, layers] : selected) {
      points.emplace_back(point, layers);
    }
    std::ranges::sort(points, [](const auto& a, const auto& b) {
      return std::make_pair(a.first.x(), a.first.y())
             < std::make_pair(b.first.x(), b.first.y());
    });
    // Driver first: path composition grows the tree from it.
    if (const auto driver = net->getDriverAccessPoint()) {
      auto it = std::ranges::find_if(
          points, [&](const auto& p) { return p.first == *driver; });
      if (it != points.end()) {
        std::rotate(points.begin(), it, it + 1);
      }
    }
    int lx = std::numeric_limits<int>::max();
    int ly = lx;
    int hx = std::numeric_limits<int>::min();
    int hy = hx;
    for (const auto& [point, layers] : points) {
      Terminal t;
      t.point = point;
      t.top_layer
          = std::clamp(std::max(layers.high(), range.min_layer), 0, top_limit);
      t.low_layer = std::min(std::max(layers.low(), 0), t.top_layer);
      state.terminals.push_back(t);
      lx = std::min(lx, point.x());
      hx = std::max(hx, point.x());
      ly = std::min(ly, point.y());
      hy = std::max(hy, point.y());
    }
    if (state.terminals.size() <= 1) {
      state.trivial = true;
      continue;
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

BonnRouteGlobal::OracleResult BonnRouteGlobal::oracle(
    const NetState& state,
    Scratch& scratch,
    const bool objective_only,
    const bool whole_grid) const
{
  OracleResult result;
  const GRNet* net = state.net;
  const BoxT region
      = whole_grid ? BoxT(0, 0, x_size_ - 1, y_size_ - 1) : state.region;
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
  const uint32_t epoch = scratch.epoch;
  auto index = [&](int l, int x, int y) -> int64_t {
    return l * plane + static_cast<int64_t>(x - x0) * height + (y - y0);
  };
  auto greater = [](const std::pair<double, int64_t>& a,
                    const std::pair<double, int64_t>& b) { return a > b; };
  auto push = [&](double d, int64_t v) {
    scratch.heap.emplace_back(d, v);
    std::push_heap(scratch.heap.begin(), scratch.heap.end(), greater);
  };
  auto relax = [&](int64_t from, int64_t to, double d) {
    if (scratch.stamp[to] != epoch || d < scratch.dist[to]) {
      scratch.stamp[to] = epoch;
      scratch.dist[to] = d;
      scratch.parent[to] = static_cast<int32_t>(from);
      push(d, to);
    }
  };

  const int num_terminals = static_cast<int>(state.terminals.size());
  std::vector<char> connected(num_terminals, 0);
  for (int i = 0; i < num_terminals; i++) {
    const Terminal& t = state.terminals[i];
    const int64_t v = index(t.top_layer, t.point.x(), t.point.y());
    scratch.term_stamp[v] = epoch;
    scratch.term[v] = i;
  }
  std::vector<std::array<int64_t, 2>> edges;  // local vertex pairs
  {
    const Terminal& t = state.terminals[0];
    const int64_t v = index(t.top_layer, t.point.x(), t.point.y());
    scratch.stamp[v] = epoch;
    scratch.dist[v] = 0.0;
    scratch.parent[v] = -1;
    scratch.tree_stamp[v] = epoch;
    push(0.0, v);
    connected[0] = 1;
  }
  // Terminals sharing the root vertex are connected already (cannot happen
  // with one terminal per gcell, kept for robustness).
  int remaining = num_terminals - 1;

  while (!scratch.heap.empty() && remaining > 0) {
    std::pop_heap(scratch.heap.begin(), scratch.heap.end(), greater);
    const auto [d, u] = scratch.heap.back();
    scratch.heap.pop_back();
    if (d > scratch.dist[u]) {
      continue;
    }
    if (scratch.term_stamp[u] == epoch && !connected[scratch.term[u]]
        && scratch.tree_stamp[u] != epoch) {
      // Algorithm 1, step 5: add the path to the tree; its vertices become
      // sources at distance 0 for the next connection.
      int64_t w = u;
      while (scratch.tree_stamp[w] != epoch) {
        const int64_t p = scratch.parent[w];
        edges.push_back({w, p});
        scratch.tree_stamp[w] = epoch;
        // The path may pass through further terminals.
        if (scratch.term_stamp[w] == epoch && !connected[scratch.term[w]]) {
          connected[scratch.term[w]] = 1;
          remaining--;
        }
        scratch.dist[w] = 0.0;
        push(0.0, w);
        w = p;
      }
      continue;
    }
    const int l = static_cast<int>(u / plane);
    const int64_t rem = u % plane;
    const int x = x0 + static_cast<int>(rem / height);
    const int y = y0 + static_cast<int>(rem % height);
    if (l >= range.min_layer && l <= range.max_layer) {
      const double factor = net->getNdrCost(l);
      if (grid_graph_->getLayerDirection(l) == MetalLayer::H) {
        if (x > x0) {
          const double c = objective_only ? wireObjective(l, x - 1, y)
                                          : wireCost(l, x - 1, y, factor);
          relax(u, u - height, d + c);
        }
        if (x < x1) {
          const double c = objective_only ? wireObjective(l, x, y)
                                          : wireCost(l, x, y, factor);
          relax(u, u + height, d + c);
        }
      } else {
        if (y > y0) {
          const double c = objective_only ? wireObjective(l, x, y - 1)
                                          : wireCost(l, x, y - 1, factor);
          relax(u, u - 1, d + c);
        }
        if (y < y1) {
          const double c = objective_only ? wireObjective(l, x, y)
                                          : wireCost(l, x, y, factor);
          relax(u, u + 1, d + c);
        }
      }
    }
    if (l + 1 < num_layers) {
      const double c
          = objective_only ? params_.via_weight : viaCost(l, x, y, net);
      relax(u, u + plane, d + c);
    }
    if (l > 0) {
      const double c
          = objective_only ? params_.via_weight : viaCost(l - 1, x, y, net);
      relax(u, u - plane, d + c);
    }
  }
  if (remaining > 0) {
    return result;
  }
  std::vector<std::array<int, 3>> tree_edges;  // (layer, x, y) pairs packed
  tree_edges.reserve(edges.size() * 2);
  for (const auto& [a, b] : edges) {
    for (const int64_t v : {a, b}) {
      const int l = static_cast<int>(v / plane);
      const int64_t rem = v % plane;
      tree_edges.push_back({l,
                            x0 + static_cast<int>(rem / height),
                            y0 + static_cast<int>(rem % height)});
    }
  }
  if (!buildCandidate(state, tree_edges, result.candidate)) {
    return result;
  }
  result.ok = true;
  return result;
}

bool BonnRouteGlobal::buildCandidate(
    const NetState& state,
    const std::vector<std::array<int, 3>>& tree_edges,
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
  std::unordered_map<uint64_t, std::vector<uint64_t>> adjacency;
  auto add_edge = [&](uint64_t a, uint64_t b) {
    adjacency[a].push_back(b);
    adjacency[b].push_back(a);
  };
  for (size_t i = 0; i + 1 < tree_edges.size(); i += 2) {
    const auto& a = tree_edges[i];
    const auto& b = tree_edges[i + 1];
    add_edge(key(a[0], a[1], a[2]), key(b[0], b[1], b[2]));
  }
  std::unordered_map<uint64_t, char> pinned;
  for (const Terminal& t : state.terminals) {
    for (int l = t.low_layer; l < t.top_layer; l++) {
      add_edge(key(l, t.point.x(), t.point.y()),
               key(l + 1, t.point.x(), t.point.y()));
    }
    pinned[key(t.low_layer, t.point.x(), t.point.y())] = 1;
    pinned[key(t.top_layer, t.point.x(), t.point.y())] = 1;
  }
  for (auto& [k, list] : adjacency) {
    std::ranges::sort(list);
    list.erase(std::ranges::unique(list).begin(), list.end());
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

  const Terminal& root_terminal = state.terminals[0];
  const uint64_t root = key(root_terminal.low_layer,
                            root_terminal.point.x(),
                            root_terminal.point.y());
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
    for (const uint64_t next : adjacency[vertex]) {
      if (visited.contains(next)) {
        continue;
      }
      // Walk a straight run of unit edges.
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
        const uint64_t after = list[0] == prev ? list[1] : list[0];
        if (visited.contains(after) || axis(cur, after) != run_axis) {
          break;
        }
        prev = cur;
        cur = after;
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

  const std::shared_ptr<GRTreeNode> tree = buildTree(candidate);
  std::vector<std::pair<uint32_t, float>> usage;
  grid_graph_->forEachTreeDemand(
      tree,
      state.net->getNdrCosts(),
      [&](const int layer, const PointT edge, const CapacityT demand) {
        usage.emplace_back(edgeIndex(layer, edge.x(), edge.y()),
                           static_cast<float>(demand));
      });
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
    hash = (hash ^ static_cast<uint64_t>(std::lround(d * 1024.0f)))
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
      nodes[rec.parent]->addChild(nodes.back());
    }
  }
  return nodes.front();
}

double BonnRouteGlobal::logCost(const Candidate& candidate) const
{
  double cost = objective_weight_ * candidate.objective;
  for (const auto& [r, d] : candidate.usage) {
    cost += weight_[r] * d;
  }
  return std::log(std::max(cost, 1e-300)) + shift_;
}

void BonnRouteGlobal::applyPrices(const Candidate& candidate)
{
  const double eps = params_.epsilon;
  for (const auto& [r, d] : candidate.usage) {
    const double lp = (log_price_[r] += eps * d / norm_capacity_[r]);
    weight_[r] = std::exp(lp - shift_) / norm_capacity_[r];
    max_log_price_ = std::max(max_log_price_, lp);
  }
  log_objective_price_ += eps * candidate.objective / objective_bound_;
  objective_weight_
      = std::exp(log_objective_price_ - shift_) / objective_bound_;
  if (std::max(max_log_price_, log_objective_price_) - shift_
      > kRenormalizeAbove) {
    renormalize();
  }
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

void BonnRouteGlobal::resourceSharing()
{
  const int threads = std::max(1, params_.threads);
  std::vector<Scratch> scratch(threads);
  const int num_nets = static_cast<int>(order_.size());

  // Congestion-free objective: guess for the objective resource bound.
  double base_objective = 0.0;
  {
    std::vector<double> objective(num_nets, 0.0);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
    for (int i = 0; i < num_nets; i++) {
      Scratch& s = scratch[omp_get_thread_num()];
      OracleResult r = oracle(states_[order_[i]], s, true, false);
      if (!r.ok) {
        r = oracle(states_[order_[i]], s, true, true);
      }
      objective[i] = r.ok ? r.candidate.objective : 0.0;
    }
    for (const double o : objective) {
      base_objective += o;
    }
  }
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
  for (int phase = 0; phase < params_.phases; phase++) {
    int64_t phase_calls = 0;
    for (int begin = 0; begin < num_nets; begin += batch) {
      const int end = std::min(num_nets, begin + batch);
      // Prices are frozen during the parallel part of the batch.
      const double batch_log_objective = log_objective_price_;
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
      for (int i = begin; i < end; i++) {
        const NetState& state = states_[order_[i]];
        const int slot = i - begin;
        results[slot].ok = false;
        if (state.last >= 0) {
          const double now = logCost(state.candidates[state.last]);
          // Relative to the objective price, which every phase raises.
          const double growth
              = (now - state.last_log_cost)
                - (batch_log_objective - state.last_log_objective);
          if (growth <= log_tolerance) {
            kinds[slot] = Kind::kReuse;
            continue;
          }
        }
        Scratch& s = scratch[omp_get_thread_num()];
        OracleResult r = oracle(state, s, false, false);
        if (!r.ok) {
          r = oracle(state, s, false, true);
        }
        if (r.ok) {
          r.log_cost = logCost(r.candidate);
          kinds[slot] = Kind::kNew;
        } else {
          kinds[slot] = state.last >= 0 ? Kind::kReuse : Kind::kFail;
        }
        results[slot] = std::move(r);
      }
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
          const double log_cost = results[slot].log_cost;
          state.last
              = insertCandidate(state, std::move(results[slot].candidate));
          state.last_log_cost = log_cost;
          state.last_log_objective = batch_log_objective;
        } else {
          reuses_++;
        }
        Candidate& chosen = state.candidates[state.last];
        chosen.count++;
        applyPrices(chosen);
      }
    }
    oracle_calls_ += phase_calls;
    const int done = phase + 1;
    if (params_.verbose
        && (done == 1 || done % 25 == 0 || done == params_.phases)) {
      // ln y_r = eps * (sum of normalized loads); divide by phases for the
      // average fractional load lambda_r.
      double lambda_max = 0.0;
      for (const double lp : log_price_) {
        lambda_max = std::max(lambda_max, lp);
      }
      lambda_max /= params_.epsilon * done;
      const double lambda_obj = (log_objective_price_ - initial_log_objective)
                                / (params_.epsilon * done);
      logger_->report(
          "BonnRoute phase {:3d}: oracle calls {:7d}, max edge load {:.3f}, "
          "objective load {:.3f}.",
          done,
          phase_calls,
          lambda_max,
          lambda_obj);
    }
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
      auto delta = [&](const std::vector<std::pair<uint32_t, float>>& alt) {
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
    std::shared_ptr<GRTreeNode> tree
        = buildTree(state.candidates[state.chosen]);
    state.net->setRoutingTree(tree);
    grid_graph_->addTreeUsage(tree, state.net->getNdrCosts());
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
  prepareNets();
  min_layer_ = num_layers_;
  for (const GRNet* net : nets_) {
    min_layer_ = std::min(min_layer_, net->getLayerRange().min_layer);
  }
  // Single-gcell nets keep their pin stack; their usage is background.
  for (NetState& state : states_) {
    if (!state.trivial || state.terminals.empty()) {
      continue;
    }
    const Terminal& t = state.terminals[0];
    auto node
        = std::make_shared<GRTreeNode>(t.low_layer, t.point.x(), t.point.y());
    if (t.top_layer != t.low_layer) {
      node->addChild(
          std::make_shared<GRTreeNode>(t.top_layer, t.point.x(), t.point.y()));
    }
    state.net->setRoutingTree(node);
    grid_graph_->addTreeUsage(node, state.net->getNdrCosts());
  }
  const size_t num_resources
      = static_cast<size_t>(num_layers_) * x_size_ * y_size_;
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
  resourceSharing();
  const double sharing_time = timer.elapsed();
  roundAndRepair();
  commit();
  if (params_.verbose) {
    size_t candidates = 0;
    for (const int index : order_) {
      candidates += states_[index].candidates.size();
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
