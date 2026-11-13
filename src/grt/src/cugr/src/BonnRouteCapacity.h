// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Global routing edge capacities after Gester et al., "BonnRoute: Algorithms
// and Data Structures for Fast and Good VLSI Routing", TODAES 18(2), 2013,
// Section 2.5 ("Global Routing Edge Capacities"):
//
//  (a) count the usable vertices of the track graph (Sec. 3.5) between the
//      tile centers of an edge, after extending blockages, and divide by the
//      number of vertices on one track of that span;
//  (b) route nets whose pins lie in a single tile inside that tile first and
//      remove the track-graph vertices they use;
//  (c) reduce capacities by the intra-tile connection length of longer nets,
//      estimated by Steiner lengths;
//  (d) reduce intermediate-layer capacities by the expected loss from
//      stacked vias, computed from a lattice counting model.
//
// The estimator depends only on the plain data of CapacityGrid and
// EstimatorNet below (tracks, blockages, gcell grid, pins), so any global
// router on a gcell grid can use it: BonnRoute fills them from CUGR's
// design (capacityGridFromCugr in BonnRoute.h), and a FastRoute adapter can
// fill them from odb and GlobalRouter's grid.
//
// Units: coordinates and lengths are DBU. Capacities are in tracks, i.e.
// the number of standard wires an edge can carry, as in CUGR and FastRoute
// before their conversion to their own edge units. Edge (l, x, y) joins
// gcell (x, y) to (x + 1, y) on a horizontal layer and to (x, y + 1) on a
// vertical layer; the last gcell of a row or column has no edge.
//
// The helpers at the top and the estimator are unit tested without a
// database (src/grt/test/cpp/TestBonnRouteCapacity.cpp).

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "odb/geom.h"

namespace grt::bonnroute {

inline constexpr int kHorizontal = 0;
inline constexpr int kVertical = 1;

// Sec. 2.5, stacked vias: a lattice of `rows` x `columns` vertices holds k
// pairwise disjoint sets of `length` consecutive vertices of one row ("in
// x-direction"), each placement equally likely. Returns E[max over columns
// of the number of chosen vertices in the column] for k = 0..max_count by
// counting, for every bound m, the placements with at most m chosen
// vertices per column. Counts for k beyond the lattice's capacity repeat
// the last feasible value.
std::vector<double> expectedMaxColumnLoad(int rows,
                                          int columns,
                                          int length,
                                          int max_count);

// Linear interpolation of a table indexed by integer counts, clamped to the
// table's range.
double interpolate(const std::vector<double>& table, double count);

// Length of a rectilinear minimum spanning tree split into its horizontal
// and vertical parts; the Steiner length estimate of Sec. 2.5 for the
// intra-tile connections (within a factor 3/2 of the Steiner minimum).
std::pair<int64_t, int64_t> rectilinearSpanningLength(
    const std::vector<odb::Point>& points);

// One routing layer of the track graph (Sec. 3.5): tracks in the preferred
// direction, and vertices where projected tracks of the neighboring wiring
// layers cross them. Stores one usable/unusable bit per vertex.
class TrackLayer
{
 public:
  TrackLayer() = default;
  // `tracks`: track coordinates across the preferred direction; `vertices`:
  // vertex coordinates along it. Both are sorted and unique.
  TrackLayer(int direction, std::vector<int> tracks, std::vector<int> vertices);

  int direction() const { return direction_; }
  const std::vector<int>& tracks() const { return tracks_; }
  const std::vector<int>& vertices() const { return vertices_; }
  int numTracks() const { return static_cast<int>(tracks_.size()); }
  int numVertices() const { return static_cast<int>(vertices_.size()); }

  bool usable(int track, int vertex) const;
  // Marks vertices [first, last] of `track` unusable.
  void block(int track, int first, int last);
  // Number of usable vertices of `track` among indices [first, last).
  int countUsable(int track, int first, int last) const;

  // Indices of the tracks / vertices with coordinate in [lo, hi]; the pair
  // is empty (first > second) if there is none.
  std::pair<int, int> trackRange(int lo, int hi) const;
  std::pair<int, int> vertexRange(int lo, int hi) const;
  // Index of the first vertex with coordinate >= coord (numVertices() if
  // none).
  int lowerVertex(int coord) const;
  // Index of the nearest track / vertex; -1 if there is none.
  int nearestTrack(int coord) const;
  int nearestVertex(int coord) const;

 private:
  int direction_ = 0;
  std::vector<int> tracks_;
  std::vector<int> vertices_;
  std::vector<uint64_t> blocked_;  // track-major bit matrix
};

// One wiring layer of the estimator's input (bottom-up order, index =
// routing level - 1).
struct CapacityLayer
{
  std::string name;
  int direction = kHorizontal;  // preferred direction
  // Track coordinates across the preferred direction (y on a horizontal
  // layer), sorted.
  std::vector<int> tracks;
  int pitch = 0;
  int width = 0;
  // Parallel-run spacing table as in LEF: spacing[w][l] applies from
  // spacing_widths[w] and spacing_lengths[l] on; one entry for a plain rule.
  std::vector<int> spacing_widths{0};
  std::vector<int> spacing_lengths{0};
  std::vector<std::vector<int>> spacing{{0}};
  // Layers below the minimum routing layer get capacity 0.
  bool routable = false;
  // User layer adjustment: fraction of the capacity removed (0 = none).
  double adjustment = 0.0;
  // Pad that a via stack passing through this layer (from the layer below
  // to the layer above) leaves on it: extent along the preferred
  // direction, tracks blocked across, and their product as the router's
  // linear demand model charges it.
  double through_pad_along = 0.0;
  double through_pad_tracks = 1.0;
  double through_pad_length = 0.0;
  // Blockages (other nets' pins, obstructions, special nets).
  std::vector<odb::Rect> obstacles;

  // Spacing for a wire of `width` next to a shape over parallel `length`,
  // looked up like CUGR's MetalLayer::getParallelSpacing.
  int parallelSpacing(int width, int length) const;
};

// The gcell grid and wiring layers.
struct CapacityGrid
{
  // gridlines[0] / [1]: x / y coordinates; gcell i spans [g[i], g[i + 1]].
  std::array<std::vector<int>, 2> gridlines;
  std::vector<CapacityLayer> layers;
  odb::Rect die;
  // The router's demand model already charges each via pad linearly (CUGR
  // via-pad flank demand); the stacked-via estimate then adds only its
  // excess over that charge.
  bool linear_via_demand = false;

  int numLayers() const { return static_cast<int>(layers.size()); }
  // Number of gcells along dimension 0 (x) or 1 (y).
  int size(int dimension) const
  {
    return static_cast<int>(gridlines[dimension].size()) - 1;
  }
  int gridline(int dimension, int index) const
  {
    return gridlines[dimension][index];
  }
  int center(int dimension, int index) const
  {
    return (gridlines[dimension][index] + gridlines[dimension][index + 1]) / 2;
  }
  // Distance between the centers of gcells `edge` and `edge + 1`.
  int edgeLength(int dimension, int edge) const
  {
    return center(dimension, edge + 1) - center(dimension, edge);
  }
  // Flat index of edge (layer, x, y), the order of capacities().
  size_t index(int layer, int x, int y) const
  {
    return (static_cast<size_t>(layer) * size(0) + x) * size(1) + y;
  }
};

// A pin as the estimator sees it: its access locations (DBU, layer index)
// and the grid cell holding the pin's preferred access point.
struct EstimatorPin
{
  std::vector<std::pair<odb::Point, int>> points;
  int cell_x = -1;
  int cell_y = -1;
};

struct EstimatorNet
{
  std::vector<EstimatorPin> pins;
  // All pins in tile (tile_x, tile_y): routed inside it by (b).
  bool single_tile = false;
  int tile_x = -1;
  int tile_y = -1;
  int min_layer = 0;
  int max_layer = 0;
};

struct CapacityEstimatorOptions
{
  // Blockage extension in the preferred direction, in units of the
  // wire's spacing plus half width (the perpendicular keep-out).
  double blockage_extension = 1.0;
  bool single_tile_nets = true;
  bool intra_tile = true;
  bool stacked_vias = true;
  // Wiring layers above the lowest routable one used for single-tile nets.
  int single_tile_layers = 2;
  // Single-tile routing area: tile grown by this fraction of its size.
  double single_tile_margin = 0.25;
};

// Per-layer totals in tracks (sum of edge capacities) for the report.
struct CapacityTotals
{
  std::vector<double> tracks;        // (a)
  std::vector<double> single_tile;   // removed by (b)
  std::vector<double> intra_tile;    // removed by (c)
  std::vector<double> stacked_vias;  // removed by (d)
  std::vector<double> estimate;      // final, after the layer adjustment
  int single_tile_routed = 0;
  int single_tile_failed = 0;
  int64_t single_tile_vertices = 0;
};

class CapacityEstimator
{
 public:
  CapacityEstimator(CapacityGrid grid, const CapacityEstimatorOptions& options);

  // (a): track graph and blockages. Must be called first.
  void buildTrackGraph();
  // (b), then (c) for the other nets.
  void addNets(const std::vector<EstimatorNet>& nets);
  // (d): `stack_counts[grid.index(l, x, y)]` is the expected number of via
  // stacks passing through layer l in gcell (x, y).
  void addStackedVias(const std::vector<double>& stack_counts);
  // Capacities in tracks indexed by CapacityGrid::index, with the user
  // layer adjustment applied as CUGR and FastRoute apply it (at least one
  // track is kept on an edge with any capacity, unless the adjustment is
  // 1); unroutable layers and the last gcell of each row keep 0.
  std::vector<double> capacities();
  const CapacityTotals& totals() const { return totals_; }
  const TrackLayer& layer(int index) const { return layers_[index]; }
  const CapacityGrid& grid() const { return grid_; }

 private:
  size_t flat(int layer, int x, int y) const
  {
    return grid_.index(layer, x, y);
  }
  // Vertex capacity (a)/(b) of every edge, before reductions.
  std::vector<double> vertexCapacities() const;
  bool routeSingleTile(const EstimatorNet& net);
  void addIntraTile(int layer, int tile_x, int tile_y, double length);
  void addIntraTileNet(const EstimatorNet& net, bool include_center);
  int lowestLayer(int direction, int min_layer, int max_layer) const;

  CapacityGrid grid_;
  CapacityEstimatorOptions options_;
  int num_layers_ = 0;
  int x_size_ = 0;
  int y_size_ = 0;
  int min_layer_ = 0;
  std::vector<TrackLayer> layers_;
  std::vector<double> stacked_;  // (d) only, per edge
  std::vector<double> intra_;    // (c) only, per edge
  std::vector<double> before_single_tile_;
  CapacityTotals totals_;
};

}  // namespace grt::bonnroute
