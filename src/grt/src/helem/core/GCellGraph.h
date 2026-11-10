// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace grt::helem {

// One routing layer of the 3D GCell grid.
struct RoutingLayer
{
  std::string name;
  bool horizontal = true;
  double overflow_weight = 1.0;  // OvflW_l in the ISPD 2024 score
  double min_length = 0.0;
};

// 3D GCell grid with per-layer capacities, following the ISPD 2024 resource
// file semantics: the capacity stored for GCell (x, y) on a horizontal layer
// is the capacity of the wire edge (x, y)-(x+1, y); on a vertical layer it is
// the capacity of the edge (x, y)-(x, y+1).
//
// Storage order per layer: horizontal layers use index y * x_size + x,
// vertical layers use x * y_size + y, so that a run of edges along the
// preferred direction is contiguous.
struct GCellGraph
{
  int num_layers = 0;
  int x_size = 0;
  int y_size = 0;
  // Layers below this index carry pins only.
  int min_routing_layer = 1;
  int max_routing_layer = -1;  // inclusive; -1 means the top layer

  double unit_wire_cost = 1.0;  // per DBU of wire length
  double unit_via_cost = 1.0;   // per via cut

  std::vector<RoutingLayer> layers;
  std::vector<int64_t> horizontal_edge_length;  // x_size - 1 entries
  std::vector<int64_t> vertical_edge_length;    // y_size - 1 entries
  std::vector<std::vector<float>> capacity;     // [layer][storage index]

  int topRoutingLayer() const
  {
    return max_routing_layer < 0 ? num_layers - 1 : max_routing_layer;
  }
  bool isRoutingLayer(int layer) const
  {
    return layer >= min_routing_layer && layer <= topRoutingLayer();
  }
  int cellCount() const { return x_size * y_size; }

  // Index of the wire edge starting at GCell (x, y) on a layer.
  int layerEdgeIndex(int layer, int x, int y) const
  {
    return layers[layer].horizontal ? y * x_size + x : x * y_size + y;
  }
  int64_t layerEdgeLength(int layer, int x, int y) const
  {
    return layers[layer].horizontal ? horizontal_edge_length[x]
                                    : vertical_edge_length[y];
  }
};

// ISPD 2024 overflow cost of one GCell edge. `demand` is in half-track units
// (a wire adds 2, a stacked via adds 1 to each adjacent edge), matching the
// official evaluator.
inline double overflowCost(double capacity, int demand)
{
  if (capacity > 0.001) {
    return std::exp(0.5 * (0.5 * demand - capacity));
  }
  if (capacity >= 0 && demand > 0) {
    return std::exp(0.5 * 1.5 * demand);
  }
  return 0.0;
}

struct GCellPoint
{
  int x = 0;
  int y = 0;
  bool operator==(const GCellPoint& other) const
  {
    return x == other.x && y == other.y;
  }
  bool operator!=(const GCellPoint& other) const { return !(*this == other); }
};

struct AccessPoint
{
  int layer = 0;
  int x = 0;
  int y = 0;
};

struct RoutingPin
{
  std::vector<AccessPoint> access_points;
};

struct RoutingNet
{
  std::string name;
  std::vector<RoutingPin> pins;
};

// A routed 3D wire or via in GCell coordinates; low <= high componentwise.
struct RouteSegment3D
{
  int x_low, y_low, layer_low;
  int x_high, y_high, layer_high;
};

struct NetRoute3D
{
  std::vector<RouteSegment3D> segments;
};

}  // namespace grt::helem
