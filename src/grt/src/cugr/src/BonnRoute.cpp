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

#include "BonnRouteCapacity.h"
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
// Destructive interference size of the x86-64/AArch64 targets we build for.
constexpr std::size_t kCacheLine = 64;
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
  uint32_t epoch = 0;

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

double BonnRouteGlobal::padCost(const int layer,
                                const int x,
                                const int y,
                                const double pad,
                                const GRNet* net) const
{
  if (!via_demand_ || pad <= 0.0) {
    return 0.0;
  }
  const int direction = grid_graph_->getLayerDirection(layer);
  const int c = direction == MetalLayer::H ? x : y;
  const int size = grid_graph_->getSize(direction);
  const int lower_len
      = c > 0 ? grid_graph_->getEdgeLength(direction, c - 1) : 0;
  const int higher_len
      = c < size - 1 ? grid_graph_->getEdgeLength(direction, c) : 0;
  const int sum = lower_len + higher_len;
  if (sum == 0) {
    return 0.0;
  }
  const double demand = pad / sum * net->getNdrCost(layer);
  double cost = 0.0;
  if (lower_len > 0) {
    const int lx = direction == MetalLayer::H ? x - 1 : x;
    const int ly = direction == MetalLayer::H ? y : y - 1;
    cost += weight_[edgeIndex(layer, lx, ly)] * demand;
  }
  if (higher_len > 0) {
    cost += weight_[edgeIndex(layer, x, y)] * demand;
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
    const bool whole_grid) const
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
  auto push = [&](double d, int64_t node) {
    scratch.heap.emplace_back(d, node);
    std::push_heap(scratch.heap.begin(), scratch.heap.end(), greater);
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
  // Algorithm 1, step 1: the component of the first terminal (driver).
  {
    const auto& c = state.terminals[0].cells.front();
    const int64_t root = index(c[0], c[1], c[2]);
    add_source(root);
    connect_at(root);
  }
  const double obj = objective_weight_;
  double last_connection_cost = 0.0;
  while (!scratch.heap.empty() && remaining > 0) {
    std::pop_heap(scratch.heap.begin(), scratch.heap.end(), greater);
    const auto [d, node] = scratch.heap.back();
    scratch.heap.pop_back();
    if (d > scratch.dist[node]) {
      continue;
    }
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
      for (const int64_t v : added) {
        add_source(v);
      }
      for (const int64_t v : added) {
        connect_at(v);
      }
      continue;
    }
    const int l = static_cast<int>(u / plane);
    const int64_t rem = u % plane;
    const int x = x0 + static_cast<int>(rem / height);
    const int y = y0 + static_cast<int>(rem % height);
    if (l >= range.min_layer && l <= range.max_layer) {
      const double factor = net->getNdrCost(l);
      auto wire = [&](int ex, int ey) {
        return objective_only ? wireObjective(l, ex, ey)
                              : wireCost(l, ex, ey, factor);
      };
      if (grid_graph_->getLayerDirection(l) == MetalLayer::H) {
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
  std::vector<std::array<int, 3>> tree_edges;  // (layer, x, y) pairs packed
  tree_edges.reserve(edges.size() * 2);
  for (const auto& [a, b] : edges) {
    tree_edges.push_back(cell(a));
    tree_edges.push_back(cell(b));
  }
  if (!buildCandidate(state, tree_edges, edge_pins, result.candidate)) {
    return result;
  }
  result.ok = true;
  return result;
}

bool BonnRouteGlobal::buildCandidate(
    const NetState& state,
    const std::vector<std::array<int, 3>>& tree_edges,
    const std::vector<int>& edge_pins,
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

  const std::shared_ptr<GRTreeNode> tree = buildTree(candidate);
  std::vector<std::pair<uint32_t, double>> usage;
  grid_graph_->forEachTreeDemand(
      tree,
      state.net->getNdrCosts(),
      [&](const int layer, const PointT edge, const CapacityT demand) {
        usage.emplace_back(edgeIndex(layer, edge.x(), edge.y()),
                           static_cast<double>(demand));
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

double BonnRouteGlobal::linearCost(const Candidate& candidate) const
{
  double cost = objective_weight_ * candidate.objective;
  for (const auto& [r, d] : candidate.usage) {
    cost += weight_[r] * d;
  }
  return cost;
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

void BonnRouteGlobal::pilot(std::vector<double>* stack_counts)
{
  const int threads = std::max(1, params_.threads);
  std::vector<Scratch> scratch(threads);
  const int num_nets = static_cast<int>(order_.size());
  // Congestion-free objective: guess for the objective resource bound.
  std::vector<double> objective(num_nets, 0.0);
  // Via columns (x, y, low, high) of each net's congestion-free tree.
  std::vector<std::vector<std::array<int, 4>>> columns(
      stack_counts != nullptr ? num_nets : 0);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (int i = 0; i < num_nets; i++) {
    Scratch& s = scratch[omp_get_thread_num()];
    OracleResult r = oracle(states_[order_[i]], s, true, false);
    if (!r.ok) {
      r = oracle(states_[order_[i]], s, true, true);
    }
    objective[i] = r.ok ? r.candidate.objective : 0.0;
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
          // Differential check of the oracle's edge costs against the
          // committed demand: exact for a single path between single cells.
          if (state.terminals.size() == 2
              && state.terminals[0].cells.size() == 1
              && state.terminals[1].cells.size() == 1) {
            const double exact = linearCost(r.candidate);
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
          cost_checks_ += results[slot].checked;
          cost_mismatches_ += results[slot].mismatch;
          if (results[slot].mismatch && cost_mismatches_ <= 5) {
            const Candidate& c = results[slot].candidate;
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
          state.last
              = insertCandidate(state, std::move(results[slot].candidate));
          state.last_log_cost = log_cost;
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
  const bool estimate = params_.capacity_model == 1;
  if (estimate) {
    estimateCapacities();
  }
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
  std::vector<double> stack_counts;
  const bool count_stacks = estimate && params_.capacity_stacked_vias;
  pilot(count_stacks ? &stack_counts : nullptr);
  if (estimate) {
    utl::Timer stack_timer;
    if (count_stacks) {
      estimator_->addStackedVias(stack_counts);
      applyCapacities(estimator_->capacities());
    }
    estimate_time_ += stack_timer.elapsed();
    reportCapacities();
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
        "BonnRoute: oracle cost check {} of {} two-pin searches differ "
        "from the committed demand cost.",
        cost_mismatches_,
        cost_checks_);
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
