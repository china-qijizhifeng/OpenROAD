// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "HelemRouter.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numeric>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace grt::helem {

namespace {

constexpr double kInfinity = std::numeric_limits<double>::infinity();
// Exponents are clamped so that prices stay finite on badly congested edges.
constexpr double kMaxExponent = 40.0;

double clampedExp(double exponent)
{
  return std::exp(std::min(exponent, kMaxExponent));
}

uint64_t splitMix64(uint64_t value)
{
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

// Axis-aligned polyline with at most five corners (3-bend pattern).
struct Path2D
{
  uint8_t count = 0;
  std::array<GCellPoint, 5> points;

  void push(const GCellPoint& point)
  {
    if (count > 0 && points[count - 1] == point) {
      return;
    }
    // Drop a middle corner that is collinear with its neighbours.
    if (count >= 2) {
      const GCellPoint& a = points[count - 2];
      const GCellPoint& b = points[count - 1];
      if ((a.x == b.x && b.x == point.x) || (a.y == b.y && b.y == point.y)) {
        points[count - 1] = point;
        if (points[count - 2] == point) {
          --count;
        }
        return;
      }
    }
    points[count++] = point;
  }
};

struct Segment
{
  int net = 0;
  int node_a = 0;
  int node_b = 0;
  GCellPoint a;
  GCellPoint b;
  // 3-bend search box, inclusive.
  int x_low = 0, y_low = 0, x_high = 0, y_high = 0;
  Path2D path;
};

struct NetTree
{
  std::vector<GCellPoint> nodes;  // pins first, then Steiner points
  int pin_node_count = 0;
  // Layer interval of the chosen access points at each pin node.
  std::vector<std::pair<int, int>> pin_layers;
  int first_segment = 0;
  int segment_count = 0;
  int64_t hpwl = 0;
};

struct LayerInterval
{
  int low = std::numeric_limits<int>::max();
  int high = std::numeric_limits<int>::min();
  bool empty() const { return low > high; }
  void add(int layer)
  {
    low = std::min(low, layer);
    high = std::max(high, layer);
  }
  int distance(int layer) const
  {
    if (layer < low) {
      return low - layer;
    }
    if (layer > high) {
      return layer - high;
    }
    return 0;
  }
};

struct ViaSpan
{
  int x;
  int y;
  int low;
  int high;
};

// Layer assignment result of one net before it is committed.
struct NetAssignment
{
  std::vector<uint64_t> wire_edges;  // (layer << 40) | layer storage index
  std::vector<ViaSpan> vias;
};

uint64_t wireKey(int layer, int index)
{
  return (static_cast<uint64_t>(layer) << 40) | static_cast<uint64_t>(index);
}

// Thread-local buffers of the 3-bend kernel and the layer assignment.
struct Scratch
{
  std::vector<double> horizontal_prefix;
  std::vector<double> vertical_prefix;
  std::vector<double> distance;
  std::vector<int16_t> via_from;
  std::vector<GCellPoint> cells;
  std::vector<int> step_layers;
  std::vector<double> via_cost_prefix;  // per (position, layer) cumulative
};

}  // namespace

class HelemRouter::Impl
{
 public:
  Impl(const GCellGraph& graph,
       const std::vector<RoutingNet>& nets,
       const HelemOptions& options,
       SteinerBuilder steiner_builder,
       HelemLogger logger)
      : graph_(graph),
        nets_(nets),
        options_(options),
        steiner_builder_(std::move(steiner_builder)),
        logger_(std::move(logger))
  {
  }

  std::vector<NetRoute3D> run();
  HelemScore score() const { return score_; }
  const std::vector<HelemStageReport>& reports() const { return reports_; }

 private:
  using Clock = std::chrono::steady_clock;

  void log(const std::string& message) const
  {
    if (options_.verbose && logger_) {
      logger_(message);
    }
  }
  int threads() const { return std::max(1, options_.threads); }

  // Preparation
  void buildGrid2D();
  void choosePinsAndBuildTrees();

  // 2D edge helpers. Horizontal 2D edge (x, y)-(x+1, y) has index
  // y * x_size + x; vertical edge (x, y)-(x, y+1) has index x * y_size + y.
  int horizontalIndex(int x, int y) const { return y * graph_.x_size + x; }
  int verticalIndex(int x, int y) const { return x * graph_.y_size + y; }

  template <typename Fn>
  void forEachEdge(const Path2D& path, Fn&& fn) const
  {
    for (int i = 0; i + 1 < path.count; ++i) {
      const GCellPoint& u = path.points[i];
      const GCellPoint& v = path.points[i + 1];
      if (u.y == v.y) {
        const int low = std::min(u.x, v.x);
        const int high = std::max(u.x, v.x);
        for (int x = low; x < high; ++x) {
          fn(true, horizontalIndex(x, u.y));
        }
      } else {
        const int low = std::min(u.y, v.y);
        const int high = std::max(u.y, v.y);
        for (int y = low; y < high; ++y) {
          fn(false, verticalIndex(u.x, y));
        }
      }
    }
  }

  void addPathDemand(const Path2D& path, int delta);
  void recomputeDemand2D();
  HelemStageReport overflowReport2D(const std::string& stage) const;
  void verifyLedger(const std::string& stage) const;

  // LR stage
  Path2D routeLShape(const Segment& segment) const;
  void runLagrangianRelaxation();

  // LEM stage
  double edgePrice(bool horizontal, int index, int demand_offset) const;
  void refreshPrice(bool horizontal, int index);
  void refreshAllPrices();
  Path2D routeThreeBend(const Segment& segment, Scratch& scratch) const;
  void runExponentialMultiplier();

  // Layer assignment
  void runLayerAssignment();
  NetAssignment assignNet(int net, Scratch& scratch) const;
  double wireCost3D(int layer, const GCellPoint& u, const GCellPoint& v) const;
  // Cost of the via cut between `layer` and `layer + 1` at a GCell: unit via
  // cost plus the overflow change of the stacked-via demand on `layer`.
  double viaCost3D(int layer, const GCellPoint& cell) const;
  void commitNet(int net, const NetAssignment& assignment);

  const GCellGraph& graph_;
  const std::vector<RoutingNet>& nets_;
  HelemOptions options_;
  SteinerBuilder steiner_builder_;
  HelemLogger logger_;

  // 2D grid
  std::vector<double> capacity_h_, capacity_v_;
  std::vector<int> demand_h_, demand_v_;
  std::vector<double> weight_x_, weight_y_;  // w_e per column / row edge
  std::vector<double> multiplier_h_, multiplier_v_;
  std::vector<double> price_h_, price_v_;
  double rho_ = 0.0;
  double bend_cost_ = 0.0;  // in 3-bend kernel cost units

  std::vector<NetTree> trees_;
  std::vector<Segment> segments_;

  // 3D state
  std::vector<std::vector<int>> demand_3d_;  // half-track units, per layer
  std::vector<NetRoute3D> routes_;
  int64_t wire_length_dbu_ = 0;
  int64_t via_count_ = 0;

  std::vector<HelemStageReport> reports_;
  HelemScore score_;
};

// ---------------------------------------------------------------------------
// Preparation

void HelemRouter::Impl::buildGrid2D()
{
  const int x_size = graph_.x_size;
  const int y_size = graph_.y_size;
  const int cells = graph_.cellCount();
  capacity_h_.assign(cells, 0.0);
  capacity_v_.assign(cells, 0.0);
  for (int layer = 0; layer < graph_.num_layers; ++layer) {
    if (!graph_.isRoutingLayer(layer)) {
      continue;
    }
    // Storage order of a layer matches the 2D index of its direction.
    std::vector<double>& target
        = graph_.layers[layer].horizontal ? capacity_h_ : capacity_v_;
    const std::vector<float>& source = graph_.capacity[layer];
    for (int i = 0; i < cells; ++i) {
      target[i] += source[i];
    }
  }
  for (int i = 0; i < cells; ++i) {
    capacity_h_[i] *= options_.capacity_2d_scale;
    capacity_v_[i] *= options_.capacity_2d_scale;
  }
  // Edges that leave the grid have no capacity.
  for (int y = 0; y < y_size; ++y) {
    capacity_h_[horizontalIndex(x_size - 1, y)] = 0.0;
  }
  for (int x = 0; x < x_size; ++x) {
    capacity_v_[verticalIndex(x, y_size - 1)] = 0.0;
  }
  demand_h_.assign(cells, 0);
  demand_v_.assign(cells, 0);

  weight_x_.assign(std::max(0, x_size - 1), 1.0);
  weight_y_.assign(std::max(0, y_size - 1), 1.0);
  double length_sum = 0.0;
  int64_t length_count = 0;
  for (int64_t length : graph_.horizontal_edge_length) {
    length_sum += length;
    ++length_count;
  }
  for (int64_t length : graph_.vertical_edge_length) {
    length_sum += length;
    ++length_count;
  }
  const double mean_length = length_count > 0 ? length_sum / length_count : 1;
  for (int x = 0; x + 1 < x_size; ++x) {
    const double length = graph_.horizontal_edge_length[x];
    weight_x_[x]
        = options_.wire_weight_mode == WireWeightMode::kNormalizedLength
              ? length / mean_length
              : graph_.unit_wire_cost * length;
  }
  for (int y = 0; y + 1 < y_size; ++y) {
    const double length = graph_.vertical_edge_length[y];
    weight_y_[y]
        = options_.wire_weight_mode == WireWeightMode::kNormalizedLength
              ? length / mean_length
              : graph_.unit_wire_cost * length;
  }
  // One bend needs at least one via; express the via cost in wire weight
  // units and scale it like the edge prices (2 tau * price).
  const double via_in_weight
      = options_.wire_weight_mode == WireWeightMode::kNormalizedLength
            ? graph_.unit_via_cost / (graph_.unit_wire_cost * mean_length)
            : graph_.unit_via_cost;
  bend_cost_ = options_.bend_cost_scale * 2.0 * options_.tau * via_in_weight;
  multiplier_h_.assign(cells, 0.0);
  multiplier_v_.assign(cells, 0.0);
  price_h_.assign(cells, 0.0);
  price_v_.assign(cells, 0.0);
}

void HelemRouter::Impl::choosePinsAndBuildTrees()
{
  const int net_count = static_cast<int>(nets_.size());
  trees_.assign(net_count, NetTree());
  std::vector<std::vector<SteinerEdge>> tree_edges(net_count);

#pragma omp parallel for num_threads(threads()) schedule(dynamic, 256)
  for (int n = 0; n < net_count; ++n) {
    const RoutingNet& net = nets_[n];
    NetTree& tree = trees_[n];
    // Pin access (gap G1): the access point closest to the net centre, then
    // the lowest layer.
    int x_low = std::numeric_limits<int>::max();
    int y_low = x_low;
    int x_high = std::numeric_limits<int>::min();
    int y_high = x_high;
    for (const RoutingPin& pin : net.pins) {
      for (const AccessPoint& ap : pin.access_points) {
        x_low = std::min(x_low, ap.x);
        x_high = std::max(x_high, ap.x);
        y_low = std::min(y_low, ap.y);
        y_high = std::max(y_high, ap.y);
      }
    }
    const double center_x = 0.5 * (x_low + x_high);
    const double center_y = 0.5 * (y_low + y_high);
    const bool prefer_high = options_.pin_access_prefer_high_layer;
    auto better_layer = [&](const AccessPoint& a, const AccessPoint& b) {
      return prefer_high ? a.layer > b.layer : a.layer < b.layer;
    };
    // First choice: closest to the net centre.
    std::vector<const AccessPoint*> choice(net.pins.size(), nullptr);
    for (size_t p = 0; p < net.pins.size(); ++p) {
      double best_distance = kInfinity;
      for (const AccessPoint& ap : net.pins[p].access_points) {
        const double distance
            = std::abs(ap.x - center_x) + std::abs(ap.y - center_y);
        if (distance < best_distance
            || (distance == best_distance && better_layer(ap, *choice[p]))) {
          choice[p] = &ap;
          best_distance = distance;
        }
      }
    }
    if (options_.pin_access_mode == PinAccessMode::kOtherPinsBox
        && net.pins.size() > 1) {
      // Second choice: closest to the box spanned by the other pins.
      std::vector<const AccessPoint*> second = choice;
      for (size_t p = 0; p < net.pins.size(); ++p) {
        int bx_low = std::numeric_limits<int>::max(), by_low = bx_low;
        int bx_high = std::numeric_limits<int>::min(), by_high = bx_high;
        for (size_t q = 0; q < net.pins.size(); ++q) {
          if (q == p || choice[q] == nullptr) {
            continue;
          }
          bx_low = std::min(bx_low, choice[q]->x);
          bx_high = std::max(bx_high, choice[q]->x);
          by_low = std::min(by_low, choice[q]->y);
          by_high = std::max(by_high, choice[q]->y);
        }
        if (bx_low > bx_high) {
          continue;
        }
        double best_distance = kInfinity;
        for (const AccessPoint& ap : net.pins[p].access_points) {
          const double distance
              = std::max({0, bx_low - ap.x, ap.x - bx_high})
                + std::max({0, by_low - ap.y, ap.y - by_high});
          if (distance < best_distance
              || (distance == best_distance && better_layer(ap, *second[p]))) {
            second[p] = &ap;
            best_distance = distance;
          }
        }
      }
      choice = second;
    }
    std::unordered_map<int64_t, int> node_of_cell;
    for (size_t p = 0; p < net.pins.size(); ++p) {
      const AccessPoint* best = choice[p];
      if (best == nullptr) {
        continue;
      }
      const int64_t key = static_cast<int64_t>(best->x) * graph_.y_size
                          + best->y;
      auto [it, inserted] = node_of_cell.try_emplace(
          key, static_cast<int>(tree.nodes.size()));
      if (inserted) {
        tree.nodes.push_back({best->x, best->y});
        tree.pin_layers.emplace_back(best->layer, best->layer);
      } else {
        auto& interval = tree.pin_layers[it->second];
        interval.first = std::min(interval.first, best->layer);
        interval.second = std::max(interval.second, best->layer);
      }
    }
    tree.pin_node_count = static_cast<int>(tree.nodes.size());
    if (tree.pin_node_count < 2) {
      continue;
    }
    int px_low = std::numeric_limits<int>::max();
    int py_low = px_low;
    int px_high = std::numeric_limits<int>::min();
    int py_high = px_high;
    for (const GCellPoint& p : tree.nodes) {
      px_low = std::min(px_low, p.x);
      px_high = std::max(px_high, p.x);
      py_low = std::min(py_low, p.y);
      py_high = std::max(py_high, p.y);
    }
    tree.hpwl = static_cast<int64_t>(px_high - px_low) + (py_high - py_low);
    if (tree.pin_node_count == 2) {
      tree_edges[n].push_back({tree.nodes[0], tree.nodes[1]});
    } else {
      tree_edges[n] = steiner_builder_(tree.nodes);
    }
  }

  // Turn Steiner edges into segments over unique tree nodes. Coincident
  // points are merged and edges that would close a cycle are dropped.
  segments_.clear();
  for (int n = 0; n < net_count; ++n) {
    NetTree& tree = trees_[n];
    tree.first_segment = static_cast<int>(segments_.size());
    if (tree.pin_node_count < 2) {
      continue;
    }
    std::unordered_map<int64_t, int> node_of_cell;
    for (int i = 0; i < static_cast<int>(tree.nodes.size()); ++i) {
      node_of_cell.emplace(
          static_cast<int64_t>(tree.nodes[i].x) * graph_.y_size
              + tree.nodes[i].y,
          i);
    }
    auto node_id = [&](const GCellPoint& p) {
      const int64_t key = static_cast<int64_t>(p.x) * graph_.y_size + p.y;
      auto [it, inserted]
          = node_of_cell.try_emplace(key, static_cast<int>(tree.nodes.size()));
      if (inserted) {
        tree.nodes.push_back(p);
      }
      return it->second;
    };
    std::vector<int> parent;
    auto find = [&](int v) {
      while (parent[v] != v) {
        parent[v] = parent[parent[v]];
        v = parent[v];
      }
      return v;
    };
    std::vector<std::pair<int, int>> edges;
    for (const SteinerEdge& edge : tree_edges[n]) {
      const int a = node_id(edge.a);
      const int b = node_id(edge.b);
      if (a != b) {
        edges.emplace_back(a, b);
      }
    }
    parent.resize(tree.nodes.size());
    std::iota(parent.begin(), parent.end(), 0);
    auto add_segment = [&](int a, int b) {
      Segment segment;
      segment.net = n;
      segment.node_a = a;
      segment.node_b = b;
      segment.a = tree.nodes[a];
      segment.b = tree.nodes[b];
      const int x_low = std::min(segment.a.x, segment.b.x);
      const int x_high = std::max(segment.a.x, segment.b.x);
      const int y_low = std::min(segment.a.y, segment.b.y);
      const int y_high = std::max(segment.a.y, segment.b.y);
      auto expansion = [&](int side) {
        const int grow = static_cast<int>(
            std::ceil(options_.bbox_expand_ratio * side));
        return std::min(options_.bbox_expand_max,
                        std::max(options_.bbox_expand_min, grow));
      };
      const int expand_x = expansion(x_high - x_low);
      const int expand_y = expansion(y_high - y_low);
      segment.x_low = std::max(0, x_low - expand_x);
      segment.x_high = std::min(graph_.x_size - 1, x_high + expand_x);
      segment.y_low = std::max(0, y_low - expand_y);
      segment.y_high = std::min(graph_.y_size - 1, y_high + expand_y);
      segments_.push_back(segment);
    };
    for (const auto& [a, b] : edges) {
      const int root_a = find(a);
      const int root_b = find(b);
      if (root_a == root_b) {
        continue;
      }
      parent[root_a] = root_b;
      add_segment(a, b);
    }
    // Pins left disconnected by the Steiner builder join node 0 directly.
    for (int pin = 1; pin < tree.pin_node_count; ++pin) {
      const int root_pin = find(pin);
      const int root_zero = find(0);
      if (root_pin != root_zero) {
        parent[root_pin] = root_zero;
        add_segment(0, pin);
      }
    }
    tree.segment_count
        = static_cast<int>(segments_.size()) - tree.first_segment;
  }
}

// ---------------------------------------------------------------------------
// 2D ledger

void HelemRouter::Impl::addPathDemand(const Path2D& path, int delta)
{
  forEachEdge(path, [&](bool horizontal, int index) {
    (horizontal ? demand_h_ : demand_v_)[index] += delta;
  });
}

void HelemRouter::Impl::recomputeDemand2D()
{
  std::fill(demand_h_.begin(), demand_h_.end(), 0);
  std::fill(demand_v_.begin(), demand_v_.end(), 0);
  for (const Segment& segment : segments_) {
    addPathDemand(segment.path, 1);
  }
}

void HelemRouter::Impl::verifyLedger(const std::string& stage) const
{
  std::vector<int> horizontal(demand_h_.size(), 0);
  std::vector<int> vertical(demand_v_.size(), 0);
  for (const Segment& segment : segments_) {
    forEachEdge(segment.path, [&](bool is_horizontal, int index) {
      (is_horizontal ? horizontal : vertical)[index] += 1;
    });
  }
  if (horizontal != demand_h_ || vertical != demand_v_) {
    log("[helem] ledger mismatch after " + stage);
    std::abort();
  }
}

HelemStageReport HelemRouter::Impl::overflowReport2D(
    const std::string& stage) const
{
  HelemStageReport report;
  report.stage = stage;
  auto scan = [&](const std::vector<int>& demand,
                  const std::vector<double>& capacity) {
    for (size_t i = 0; i < demand.size(); ++i) {
      const double over = demand[i] - capacity[i];
      if (over > 0) {
        const int rounded = static_cast<int>(std::ceil(over));
        report.total_overflow_2d += rounded;
        report.max_overflow_2d = std::max(report.max_overflow_2d, rounded);
        ++report.overflowed_edges_2d;
      }
    }
  };
  scan(demand_h_, capacity_h_);
  scan(demand_v_, capacity_v_);
  return report;
}

// ---------------------------------------------------------------------------
// LR stage: L-shape routing with linear multipliers (paper eq. 8)

Path2D HelemRouter::Impl::routeLShape(const Segment& segment) const
{
  const GCellPoint& s = segment.a;
  const GCellPoint& t = segment.b;
  auto cost_of = [&](const Path2D& path) {
    double cost = 0.0;
    forEachEdge(path, [&](bool horizontal, int index) {
      cost += horizontal ? price_h_[index] : price_v_[index];
    });
    return cost;
  };
  Path2D horizontal_first;
  horizontal_first.push(s);
  horizontal_first.push({t.x, s.y});
  horizontal_first.push(t);
  if (s.x == t.x || s.y == t.y) {
    return horizontal_first;
  }
  Path2D vertical_first;
  vertical_first.push(s);
  vertical_first.push({s.x, t.y});
  vertical_first.push(t);
  return cost_of(vertical_first) < cost_of(horizontal_first) ? vertical_first
                                                             : horizontal_first;
}

void HelemRouter::Impl::runLagrangianRelaxation()
{
  const int x_size = graph_.x_size;
  const int y_size = graph_.y_size;
  // Price = w_e + y_e with y starting at zero.
  auto set_prices = [&]() {
    for (int y = 0; y < y_size; ++y) {
      for (int x = 0; x + 1 < x_size; ++x) {
        const int index = horizontalIndex(x, y);
        price_h_[index] = weight_x_[x] + multiplier_h_[index];
      }
    }
    for (int x = 0; x < x_size; ++x) {
      for (int y = 0; y + 1 < y_size; ++y) {
        const int index = verticalIndex(x, y);
        price_v_[index] = weight_y_[y] + multiplier_v_[index];
      }
    }
  };
  std::fill(multiplier_h_.begin(), multiplier_h_.end(), 0.0);
  std::fill(multiplier_v_.begin(), multiplier_v_.end(), 0.0);
  const int segment_count = static_cast<int>(segments_.size());
  const int iterations = std::max(1, options_.lr_iterations);
  for (int k = 1; k <= iterations; ++k) {
    set_prices();
#pragma omp parallel for num_threads(threads()) schedule(dynamic, 1024)
    for (int i = 0; i < segment_count; ++i) {
      segments_[i].path = routeLShape(segments_[i]);
    }
    recomputeDemand2D();
    if (k == iterations) {
      break;
    }
    const double step = 1.0 / (100.0 * k);
    const int cells = graph_.cellCount();
    for (int i = 0; i < cells; ++i) {
      multiplier_h_[i] = std::max(
          0.0, multiplier_h_[i] + step * (demand_h_[i] - capacity_h_[i]));
      multiplier_v_[i] = std::max(
          0.0, multiplier_v_[i] + step * (demand_v_[i] - capacity_v_[i]));
    }
  }
}

// ---------------------------------------------------------------------------
// LEM stage (paper Algorithm 1, eq. 13 and 23)

double HelemRouter::Impl::edgePrice(bool horizontal,
                                    int index,
                                    int demand_offset) const
{
  // pi_e = w_e + y_e * exp(rho * (d_e - Cap_e))
  double weight;
  if (horizontal) {
    weight = weight_x_[index % graph_.x_size];
  } else {
    weight = weight_y_[index % graph_.y_size];
  }
  const int demand
      = (horizontal ? demand_h_[index] : demand_v_[index]) + demand_offset;
  const double capacity = horizontal ? capacity_h_[index] : capacity_v_[index];
  const double multiplier
      = horizontal ? multiplier_h_[index] : multiplier_v_[index];
  return weight + multiplier * clampedExp(rho_ * (demand - capacity));
}

void HelemRouter::Impl::refreshPrice(bool horizontal, int index)
{
  (horizontal ? price_h_ : price_v_)[index] = edgePrice(horizontal, index, 0);
}

void HelemRouter::Impl::refreshAllPrices()
{
  const int x_size = graph_.x_size;
  const int y_size = graph_.y_size;
#pragma omp parallel for num_threads(threads()) schedule(static)
  for (int y = 0; y < y_size; ++y) {
    for (int x = 0; x + 1 < x_size; ++x) {
      refreshPrice(true, horizontalIndex(x, y));
    }
  }
#pragma omp parallel for num_threads(threads()) schedule(static)
  for (int x = 0; x < x_size; ++x) {
    for (int y = 0; y + 1 < y_size; ++y) {
      refreshPrice(false, verticalIndex(x, y));
    }
  }
}

Path2D HelemRouter::Impl::routeThreeBend(const Segment& segment,
                                         Scratch& scratch) const
{
  const int x_low = segment.x_low;
  const int y_low = segment.y_low;
  const int width = segment.x_high - x_low + 1;
  const int height = segment.y_high - y_low + 1;
  const double two_tau = 2.0 * options_.tau;

  // Edge cost of eq. (23): 1 - 2 x_old + 2 tau * price.
  std::vector<double>& h = scratch.horizontal_prefix;
  std::vector<double>& v = scratch.vertical_prefix;
  h.assign(static_cast<size_t>(width) * height, 0.0);
  v.assign(static_cast<size_t>(width) * height, 0.0);
  for (int r = 0; r < height; ++r) {
    const int y = y_low + r;
    double* row = &h[static_cast<size_t>(r) * width];
    for (int i = 0; i + 1 < width; ++i) {
      row[i + 1] = 1.0 + two_tau * price_h_[horizontalIndex(x_low + i, y)];
    }
  }
  for (int c = 0; c < width; ++c) {
    const int x = x_low + c;
    double* column = &v[static_cast<size_t>(c) * height];
    for (int j = 0; j + 1 < height; ++j) {
      column[j + 1] = 1.0 + two_tau * price_v_[verticalIndex(x, y_low + j)];
    }
  }
  // Proximal term for the current route, and optionally remove the route's
  // own demand from the price (gap G7).
  forEachEdge(segment.path, [&](bool horizontal, int index) {
    double adjust = -2.0;
    if (!options_.price_includes_own_route) {
      adjust += two_tau
                * (edgePrice(horizontal, index, -1)
                   - (horizontal ? price_h_[index] : price_v_[index]));
    }
    if (horizontal) {
      const int y = index / graph_.x_size;
      const int x = index % graph_.x_size;
      if (y >= y_low && y < y_low + height && x >= x_low
          && x + 1 < x_low + width) {
        h[static_cast<size_t>(y - y_low) * width + (x - x_low) + 1] += adjust;
      }
    } else {
      const int x = index / graph_.y_size;
      const int y = index % graph_.y_size;
      if (x >= x_low && x < x_low + width && y >= y_low
          && y + 1 < y_low + height) {
        v[static_cast<size_t>(x - x_low) * height + (y - y_low) + 1] += adjust;
      }
    }
  });
  // Prefix sums along rows and columns.
  for (int r = 0; r < height; ++r) {
    double* row = &h[static_cast<size_t>(r) * width];
    for (int i = 1; i < width; ++i) {
      row[i] += row[i - 1];
    }
  }
  for (int c = 0; c < width; ++c) {
    double* column = &v[static_cast<size_t>(c) * height];
    for (int j = 1; j < height; ++j) {
      column[j] += column[j - 1];
    }
  }
  auto h_cost = [&](int y, int x1, int x2) {
    const double* row = &h[static_cast<size_t>(y - y_low) * width];
    return std::abs(row[x2 - x_low] - row[x1 - x_low]);
  };
  auto v_cost = [&](int x, int y1, int y2) {
    const double* column = &v[static_cast<size_t>(x - x_low) * height];
    return std::abs(column[y2 - y_low] - column[y1 - y_low]);
  };

  const GCellPoint& s = segment.a;
  const GCellPoint& t = segment.b;
  double best_cost = kInfinity;
  int best_x = s.x;
  int best_y = s.y;
  bool best_first_horizontal = true;
  bool best_second_horizontal = true;
  // Mid points are scanned in index order; ties keep the first (gap G9).
  for (int j = y_low; j < y_low + height; ++j) {
    for (int i = x_low; i < x_low + width; ++i) {
      const double first_h = h_cost(s.y, s.x, i) + v_cost(i, s.y, j);
      const double first_v = v_cost(s.x, s.y, j) + h_cost(j, s.x, i);
      const double second_h = h_cost(j, i, t.x) + v_cost(t.x, j, t.y);
      const double second_v = v_cost(i, j, t.y) + h_cost(t.y, i, t.x);
      const bool use_first_h = first_h <= first_v;
      const bool use_second_h = second_h <= second_v;
      double cost = std::min(first_h, first_v) + std::min(second_h, second_v);
      if (bend_cost_ > 0) {
        Path2D candidate;
        candidate.push(s);
        candidate.push(use_first_h ? GCellPoint{i, s.y} : GCellPoint{s.x, j});
        candidate.push({i, j});
        candidate.push(use_second_h ? GCellPoint{t.x, j} : GCellPoint{i, t.y});
        candidate.push(t);
        cost += bend_cost_ * std::max(0, candidate.count - 2);
      }
      if (cost < best_cost) {
        best_cost = cost;
        best_x = i;
        best_y = j;
        best_first_horizontal = use_first_h;
        best_second_horizontal = use_second_h;
      }
    }
  }
  Path2D path;
  path.push(s);
  path.push(best_first_horizontal ? GCellPoint{best_x, s.y}
                                  : GCellPoint{s.x, best_y});
  path.push({best_x, best_y});
  path.push(best_second_horizontal ? GCellPoint{t.x, best_y}
                                   : GCellPoint{best_x, t.y});
  path.push(t);
  return path;
}

void HelemRouter::Impl::runExponentialMultiplier()
{
  const int segment_count = static_cast<int>(segments_.size());
  // Batch order: deterministic pseudo-random interleave (gap G5).
  std::vector<int> order(segment_count);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [](int a, int b) {
    const uint64_t ha = splitMix64(a);
    const uint64_t hb = splitMix64(b);
    return ha != hb ? ha < hb : a < b;
  });

  double multiplier_initial = options_.multiplier_initial;
  if (multiplier_initial < 0) {
    double weight_sum = 0.0;
    int layers = 0;
    for (int layer = 0; layer < graph_.num_layers; ++layer) {
      if (graph_.isRoutingLayer(layer)) {
        weight_sum += graph_.layers[layer].overflow_weight;
        ++layers;
      }
    }
    multiplier_initial = layers > 0 ? weight_sum / (100.0 * layers) : 0.0;
  }
  std::fill(multiplier_h_.begin(), multiplier_h_.end(), multiplier_initial);
  std::fill(multiplier_v_.begin(), multiplier_v_.end(), multiplier_initial);
  rho_ = options_.rho_initial;

  const int batch_size = std::max(1, options_.batch_size);
  const int thread_count = threads();
  std::vector<Scratch> scratch(thread_count);
  std::vector<Path2D> new_paths(batch_size);

  for (int k = 0; k < options_.lem_iterations; ++k) {
    const auto start = Clock::now();
    refreshAllPrices();
    for (int pass = 0; pass < options_.bcd_passes; ++pass) {
      for (int begin = 0; begin < segment_count; begin += batch_size) {
        const int end = std::min(segment_count, begin + batch_size);
#pragma omp parallel for num_threads(thread_count) schedule(dynamic, 8)
        for (int i = begin; i < end; ++i) {
#ifdef _OPENMP
          Scratch& local = scratch[omp_get_thread_num()];
#else
          Scratch& local = scratch[0];
#endif
          new_paths[i - begin] = routeThreeBend(segments_[order[i]], local);
        }
        // Commit in batch order, then refresh prices on touched edges.
        for (int i = begin; i < end; ++i) {
          Segment& segment = segments_[order[i]];
          addPathDemand(segment.path, -1);
          addPathDemand(new_paths[i - begin], 1);
        }
        for (int i = begin; i < end; ++i) {
          Segment& segment = segments_[order[i]];
          forEachEdge(segment.path, [&](bool horizontal, int index) {
            refreshPrice(horizontal, index);
          });
          segment.path = new_paths[i - begin];
          forEachEdge(segment.path, [&](bool horizontal, int index) {
            refreshPrice(horizontal, index);
          });
        }
      }
    }
    // y_e <- y_e * exp(rho * g_e), rho <- sigma * rho
    const int cells = graph_.cellCount();
    for (int i = 0; i < cells; ++i) {
      multiplier_h_[i]
          *= clampedExp(rho_ * (demand_h_[i] - capacity_h_[i]));
      multiplier_v_[i]
          *= clampedExp(rho_ * (demand_v_[i] - capacity_v_[i]));
    }
    rho_ *= options_.rho_growth;
    HelemStageReport report
        = overflowReport2D("lem-" + std::to_string(k + 1));
    report.seconds
        = std::chrono::duration<double>(Clock::now() - start).count();
    reports_.push_back(report);
    log("[helem] " + report.stage
        + " overflow_total=" + std::to_string(report.total_overflow_2d)
        + " overflow_max=" + std::to_string(report.max_overflow_2d)
        + " overflow_edges=" + std::to_string(report.overflowed_edges_2d)
        + " seconds=" + std::to_string(report.seconds));
    if (options_.verify_ledger) {
      verifyLedger(report.stage);
    }
  }
}

// ---------------------------------------------------------------------------
// Layer assignment: shortest path over (position along the 2D path) x
// (layer), solved exactly by dynamic programming with up/down via sweeps.

double HelemRouter::Impl::wireCost3D(int layer,
                                     const GCellPoint& u,
                                     const GCellPoint& v) const
{
  const RoutingLayer& routing_layer = graph_.layers[layer];
  const bool horizontal_step = u.y == v.y;
  if (!graph_.isRoutingLayer(layer)
      || routing_layer.horizontal != horizontal_step) {
    return kInfinity;
  }
  const int x = std::min(u.x, v.x);
  const int y = std::min(u.y, v.y);
  const int index = graph_.layerEdgeIndex(layer, x, y);
  const int demand = demand_3d_[layer][index];
  const double capacity = graph_.capacity[layer][index];
  return graph_.unit_wire_cost * graph_.layerEdgeLength(layer, x, y)
         + routing_layer.overflow_weight
               * (overflowCost(capacity, demand + 2)
                  - overflowCost(capacity, demand));
}

double HelemRouter::Impl::viaCost3D(int layer, const GCellPoint& cell) const
{
  double cost = graph_.unit_via_cost;
  if (!graph_.isRoutingLayer(layer)) {
    return cost;
  }
  const bool horizontal = graph_.layers[layer].horizontal;
  const int position = horizontal ? cell.x : cell.y;
  const int limit = (horizontal ? graph_.x_size : graph_.y_size) - 1;
  const std::vector<int>& demand = demand_3d_[layer];
  const std::vector<float>& capacity = graph_.capacity[layer];
  auto delta = [&](int index, int amount) {
    return overflowCost(capacity[index], demand[index] + amount)
           - overflowCost(capacity[index], demand[index]);
  };
  const int here = graph_.layerEdgeIndex(layer, cell.x, cell.y);
  const int before = horizontal
                         ? graph_.layerEdgeIndex(layer, cell.x - 1, cell.y)
                         : graph_.layerEdgeIndex(layer, cell.x, cell.y - 1);
  double overflow = 0.0;
  if (position > 0 && position < limit) {
    overflow = delta(before, 1) + delta(here, 1);
  } else if (position > 0) {
    overflow = delta(before, 2);
  } else if (position < limit) {
    overflow = delta(here, 2);
  }
  return cost + graph_.layers[layer].overflow_weight * overflow;
}

NetAssignment HelemRouter::Impl::assignNet(int net, Scratch& scratch) const
{
  NetAssignment result;
  const NetTree& tree = trees_[net];
  const int layers = graph_.num_layers;
  const int y_size = graph_.y_size;

  std::unordered_map<int64_t, LayerInterval> spans;
  auto cell_key = [&](const GCellPoint& p) {
    return static_cast<int64_t>(p.x) * y_size + p.y;
  };
  for (int pin = 0; pin < tree.pin_node_count; ++pin) {
    LayerInterval& interval = spans[cell_key(tree.nodes[pin])];
    interval.add(tree.pin_layers[pin].first);
    interval.add(tree.pin_layers[pin].second);
  }

  // Orient segments by breadth-first search from node 0.
  std::vector<std::vector<std::pair<int, int>>> adjacency(tree.nodes.size());
  for (int s = 0; s < tree.segment_count; ++s) {
    const Segment& segment = segments_[tree.first_segment + s];
    adjacency[segment.node_a].emplace_back(segment.node_b, s);
    adjacency[segment.node_b].emplace_back(segment.node_a, s);
  }
  std::vector<std::pair<int, bool>> visit_order;  // segment, forward
  std::vector<char> visited(tree.nodes.size(), 0);
  std::vector<int> queue{0};
  visited[0] = 1;
  for (size_t head = 0; head < queue.size(); ++head) {
    const int node = queue[head];
    for (const auto& [next, s] : adjacency[node]) {
      if (visited[next]) {
        continue;
      }
      visited[next] = 1;
      queue.push_back(next);
      visit_order.emplace_back(
          s, segments_[tree.first_segment + s].node_a == node);
    }
  }

  std::vector<GCellPoint>& cells = scratch.cells;
  std::vector<double>& distance = scratch.distance;
  std::vector<int16_t>& via_from = scratch.via_from;
  std::vector<int>& step_layers = scratch.step_layers;
  for (const auto& [s, forward] : visit_order) {
    const Segment& segment = segments_[tree.first_segment + s];
    // Cell sequence from the parent node to the child node.
    cells.clear();
    const Path2D& path = segment.path;
    for (int i = 0; i + 1 < path.count; ++i) {
      const GCellPoint& u = path.points[i];
      const GCellPoint& v = path.points[i + 1];
      const int dx = (v.x > u.x) - (v.x < u.x);
      const int dy = (v.y > u.y) - (v.y < u.y);
      GCellPoint p = u;
      while (p != v) {
        cells.push_back(p);
        p.x += dx;
        p.y += dy;
      }
    }
    cells.push_back(path.points[path.count - 1]);
    if (!forward) {
      std::reverse(cells.begin(), cells.end());
    }
    const int steps = static_cast<int>(cells.size()) - 1;
    if (steps <= 0) {
      continue;
    }
    // via_prefix[p * (layers + 1) + l] = sum of via cut costs below layer l.
    std::vector<double>& via_prefix = scratch.via_cost_prefix;
    via_prefix.assign(static_cast<size_t>(steps + 1) * (layers + 1), 0.0);
    for (int p = 0; p <= steps; ++p) {
      double* prefix = &via_prefix[static_cast<size_t>(p) * (layers + 1)];
      for (int l = 0; l + 1 < layers; ++l) {
        prefix[l + 1] = prefix[l] + viaCost3D(l, cells[p]);
      }
    }
    auto interval_cost = [&](int p, const LayerInterval& interval, int l) {
      const double* prefix = &via_prefix[static_cast<size_t>(p) * (layers + 1)];
      if (l < interval.low) {
        return prefix[interval.low] - prefix[l];
      }
      if (l > interval.high) {
        return prefix[l] - prefix[interval.high];
      }
      return 0.0;
    };
    const LayerInterval start = spans[cell_key(cells.front())];
    auto end_it = spans.find(cell_key(cells.back()));
    const bool has_end = end_it != spans.end();
    const LayerInterval end = has_end ? end_it->second : LayerInterval();

    distance.assign(static_cast<size_t>(steps + 1) * layers, kInfinity);
    via_from.assign(static_cast<size_t>(steps + 1) * layers, -1);
    for (int l = 0; l < layers; ++l) {
      distance[l] = start.empty() ? 0.0 : interval_cost(0, start, l);
      via_from[l] = static_cast<int16_t>(l);
    }
    for (int p = 0; p < steps; ++p) {
      const double* current = &distance[static_cast<size_t>(p) * layers];
      double* next = &distance[static_cast<size_t>(p + 1) * layers];
      int16_t* next_from = &via_from[static_cast<size_t>(p + 1) * layers];
      for (int l = 0; l < layers; ++l) {
        const double wire = wireCost3D(l, cells[p], cells[p + 1]);
        next[l] = current[l] + wire;
        next_from[l] = static_cast<int16_t>(l);
      }
      if (p + 1 == steps) {
        break;
      }
      // Via sweeps at the next cell (up, then down).
      const double* prefix
          = &via_prefix[static_cast<size_t>(p + 1) * (layers + 1)];
      for (int l = 1; l < layers; ++l) {
        const double cut = prefix[l] - prefix[l - 1];
        if (next[l - 1] + cut < next[l]) {
          next[l] = next[l - 1] + cut;
          next_from[l] = next_from[l - 1];
        }
      }
      for (int l = layers - 2; l >= 0; --l) {
        const double cut = prefix[l + 1] - prefix[l];
        if (next[l + 1] + cut < next[l]) {
          next[l] = next[l + 1] + cut;
          next_from[l] = next_from[l + 1];
        }
      }
    }
    const double* last = &distance[static_cast<size_t>(steps) * layers];
    int best_layer = -1;
    double best = kInfinity;
    for (int l = 0; l < layers; ++l) {
      const double cost
          = last[l] + (has_end ? interval_cost(steps, end, l) : 0.0);
      if (cost < best) {
        best = cost;
        best_layer = l;
      }
    }
    if (best_layer < 0) {
      continue;  // no legal layer; leaves the segment open
    }
    // Backtrack: the wire of step p - 1 arrives at cell p on `layer`; the via
    // sweep at cell p - 1 brought the wire from via_from.
    step_layers.assign(steps, 0);
    int layer = best_layer;
    for (int p = steps; p >= 1; --p) {
      step_layers[p - 1] = layer;
      if (p - 1 >= 1) {
        layer = via_from[static_cast<size_t>(p - 1) * layers + layer];
      }
    }
    // Record wires and via spans.
    for (int p = 0; p < steps; ++p) {
      const int l = step_layers[p];
      const int x = std::min(cells[p].x, cells[p + 1].x);
      const int y = std::min(cells[p].y, cells[p + 1].y);
      result.wire_edges.push_back(wireKey(l, graph_.layerEdgeIndex(l, x, y)));
    }
    spans[cell_key(cells.front())].add(step_layers.front());
    for (int p = 1; p < steps; ++p) {
      if (step_layers[p - 1] != step_layers[p]) {
        LayerInterval& interval = spans[cell_key(cells[p])];
        interval.add(step_layers[p - 1]);
        interval.add(step_layers[p]);
      }
    }
    spans[cell_key(cells.back())].add(step_layers.back());
  }

  std::sort(result.wire_edges.begin(), result.wire_edges.end());
  result.wire_edges.erase(
      std::unique(result.wire_edges.begin(), result.wire_edges.end()),
      result.wire_edges.end());
  bool has_via = false;
  for (const auto& [key, interval] : spans) {
    if (interval.low < interval.high) {
      result.vias.push_back({static_cast<int>(key / y_size),
                             static_cast<int>(key % y_size),
                             interval.low,
                             interval.high});
      has_via = true;
    }
  }
  // A multi-pin net whose pins share one GCell and layer still needs a
  // routing object so that the evaluator sees the pins connected.
  if (result.wire_edges.empty() && !has_via && nets_[net].pins.size() > 1
      && tree.pin_node_count == 1) {
    const int pin_layer = tree.pin_layers[0].first;
    const int low = pin_layer + 1 < layers ? pin_layer : pin_layer - 1;
    result.vias.push_back(
        {tree.nodes[0].x, tree.nodes[0].y, low, low + 1});
  }
  std::sort(result.vias.begin(),
            result.vias.end(),
            [](const ViaSpan& a, const ViaSpan& b) {
              return std::tie(a.x, a.y) < std::tie(b.x, b.y);
            });
  return result;
}

void HelemRouter::Impl::commitNet(int net, const NetAssignment& assignment)
{
  NetRoute3D& route = routes_[net];
  route.segments.clear();
  const int x_size = graph_.x_size;
  const int y_size = graph_.y_size;
  // Cells covered by this net's wires, for the stacked-via rule.
  std::vector<uint64_t> covered;
  covered.reserve(assignment.wire_edges.size() * 2);

  size_t i = 0;
  while (i < assignment.wire_edges.size()) {
    const uint64_t key = assignment.wire_edges[i];
    const int layer = static_cast<int>(key >> 40);
    size_t j = i;
    while (j + 1 < assignment.wire_edges.size()
           && assignment.wire_edges[j + 1] == assignment.wire_edges[j] + 1) {
      ++j;
    }
    const int first = static_cast<int>(assignment.wire_edges[i] & 0xffffffffffULL);
    const int last = static_cast<int>(assignment.wire_edges[j] & 0xffffffffffULL);
    std::vector<int>& demand = demand_3d_[layer];
    for (int index = first; index <= last; ++index) {
      demand[index] += 2;
      covered.push_back(wireKey(layer, index));
      covered.push_back(wireKey(layer, index + 1));
    }
    if (graph_.layers[layer].horizontal) {
      const int y = first / x_size;
      const int x_first = first % x_size;
      const int x_last = last % x_size;
      for (int x = x_first; x <= x_last; ++x) {
        wire_length_dbu_ += graph_.horizontal_edge_length[x];
      }
      route.segments.push_back({x_first, y, layer, x_last + 1, y, layer});
    } else {
      const int x = first / y_size;
      const int y_first = first % y_size;
      const int y_last = last % y_size;
      for (int y = y_first; y <= y_last; ++y) {
        wire_length_dbu_ += graph_.vertical_edge_length[y];
      }
      route.segments.push_back({x, y_first, layer, x, y_last + 1, layer});
    }
    i = j + 1;
  }
  std::sort(covered.begin(), covered.end());
  for (const ViaSpan& via : assignment.vias) {
    route.segments.push_back(
        {via.x, via.y, via.low, via.x, via.y, via.high});
    via_count_ += via.high - via.low;
    for (int layer = via.low; layer < via.high; ++layer) {
      const int index = graph_.layerEdgeIndex(layer, via.x, via.y);
      if (std::binary_search(
              covered.begin(), covered.end(), wireKey(layer, index))) {
        continue;
      }
      // Stacked via: half a track on the adjacent edges of the lower layer.
      std::vector<int>& demand = demand_3d_[layer];
      if (graph_.layers[layer].horizontal) {
        const bool has_left = via.x > 0;
        const bool has_right = via.x < x_size - 1;
        if (has_left && has_right) {
          demand[graph_.layerEdgeIndex(layer, via.x - 1, via.y)] += 1;
          demand[index] += 1;
        } else if (has_left) {
          demand[graph_.layerEdgeIndex(layer, via.x - 1, via.y)] += 2;
        } else if (has_right) {
          demand[index] += 2;
        }
      } else {
        const bool has_below = via.y > 0;
        const bool has_above = via.y < y_size - 1;
        if (has_below && has_above) {
          demand[graph_.layerEdgeIndex(layer, via.x, via.y - 1)] += 1;
          demand[index] += 1;
        } else if (has_below) {
          demand[graph_.layerEdgeIndex(layer, via.x, via.y - 1)] += 2;
        } else if (has_above) {
          demand[index] += 2;
        }
      }
    }
  }
}

void HelemRouter::Impl::runLayerAssignment()
{
  const int net_count = static_cast<int>(nets_.size());
  demand_3d_.assign(graph_.num_layers,
                    std::vector<int>(graph_.cellCount(), 0));
  routes_.assign(net_count, NetRoute3D());
  wire_length_dbu_ = 0;
  via_count_ = 0;

  std::vector<int> order(net_count);
  std::iota(order.begin(), order.end(), 0);
  const bool short_first
      = options_.layer_assignment_order == LayerAssignmentOrder::kShortFirst;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
    return short_first ? trees_[a].hpwl < trees_[b].hpwl
                       : trees_[a].hpwl > trees_[b].hpwl;
  });

  const int batch_size = std::max(1, options_.layer_assignment_batch_size);
  const int thread_count = threads();
  std::vector<Scratch> scratch(thread_count);
  std::vector<NetAssignment> batch(batch_size);
  for (int begin = 0; begin < net_count; begin += batch_size) {
    const int end = std::min(net_count, begin + batch_size);
#pragma omp parallel for num_threads(thread_count) schedule(dynamic, 4)
    for (int i = begin; i < end; ++i) {
#ifdef _OPENMP
      Scratch& local = scratch[omp_get_thread_num()];
#else
      Scratch& local = scratch[0];
#endif
      batch[i - begin] = assignNet(order[i], local);
    }
    for (int i = begin; i < end; ++i) {
      commitNet(order[i], batch[i - begin]);
    }
  }

  score_ = HelemScore();
  score_.wire_cost = graph_.unit_wire_cost * wire_length_dbu_;
  score_.via_count = via_count_;
  score_.via_cost = graph_.unit_via_cost * via_count_;
  for (int layer = 0; layer < graph_.num_layers; ++layer) {
    if (!graph_.isRoutingLayer(layer)) {
      continue;
    }
    double layer_cost = 0.0;
    const std::vector<int>& demand = demand_3d_[layer];
    const std::vector<float>& capacity = graph_.capacity[layer];
    for (int i = 0; i < graph_.cellCount(); ++i) {
      layer_cost += overflowCost(capacity[i], demand[i]);
    }
    score_.overflow_cost += graph_.layers[layer].overflow_weight * layer_cost;
  }
}

// ---------------------------------------------------------------------------

std::vector<NetRoute3D> HelemRouter::Impl::run()
{
  auto stage = [&](const std::string& name, auto&& body, bool overflow) {
    const auto start = Clock::now();
    body();
    HelemStageReport report;
    if (overflow) {
      report = overflowReport2D(name);
    }
    report.stage = name;
    report.seconds
        = std::chrono::duration<double>(Clock::now() - start).count();
    reports_.push_back(report);
    std::string message = "[helem] " + name
                          + " seconds=" + std::to_string(report.seconds);
    if (overflow) {
      message += " overflow_total=" + std::to_string(report.total_overflow_2d)
                 + " overflow_max=" + std::to_string(report.max_overflow_2d)
                 + " overflow_edges="
                 + std::to_string(report.overflowed_edges_2d);
    }
    log(message);
  };

  stage("prepare", [&] {
    buildGrid2D();
    choosePinsAndBuildTrees();
    log("[helem] nets=" + std::to_string(nets_.size())
        + " segments=" + std::to_string(segments_.size()));
  }, false);
  stage("lr", [&] { runLagrangianRelaxation(); }, true);
  if (options_.verify_ledger) {
    verifyLedger("lr");
  }
  if (options_.lem_iterations > 0) {
    const auto start = Clock::now();
    runExponentialMultiplier();
    HelemStageReport report = overflowReport2D("lem");
    report.seconds
        = std::chrono::duration<double>(Clock::now() - start).count();
    reports_.push_back(report);
  }
  stage("layer-assignment", [&] { runLayerAssignment(); }, false);
  char buffer[256];
  std::snprintf(buffer,
                sizeof(buffer),
                "[helem] score wire=%.4f via=%.4f overflow=%.4f total=%.4f "
                "vias=%lld",
                score_.wire_cost,
                score_.via_cost,
                score_.overflow_cost,
                score_.total(),
                static_cast<long long>(score_.via_count));
  log(buffer);
  return std::move(routes_);
}

// ---------------------------------------------------------------------------

HelemRouter::HelemRouter(const GCellGraph& graph,
                         const std::vector<RoutingNet>& nets,
                         const HelemOptions& options,
                         SteinerBuilder steiner_builder,
                         HelemLogger logger)
    : impl_(new Impl(graph,
                     nets,
                     options,
                     std::move(steiner_builder),
                     std::move(logger)))
{
}

HelemRouter::~HelemRouter()
{
  delete impl_;
}

std::vector<NetRoute3D> HelemRouter::run()
{
  return impl_->run();
}

const std::vector<HelemStageReport>& HelemRouter::stageReports() const
{
  return impl_->reports();
}

HelemScore HelemRouter::score() const
{
  return impl_->score();
}

}  // namespace grt::helem
