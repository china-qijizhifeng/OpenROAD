// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Off-track path search on the path-preserving digraph G of the blockage
// grid (Gester et al., TODAES 18(2), 2013, Sec. 3.8; Theorem 3.2 from
// Maßberg and Nieberg 2012) and the diff-net legality adapter it is driven
// with.

#pragma once

#include <cstdint>
#include <functional>
#include <limits>
#include <vector>

#include "br/BrBlockageGrid.h"
#include "br/BrGeometry.h"

namespace drt::br {

// Diff-net legality of off-track wiring: the role of the distance rule
// checking module (Sec. 3.4) for the off-track search, the access path
// catalogues and the post-processing. BonnRouteDR implements it with its
// current shape index (shapeConflicts); the shape grid of Sec. 3.3-3.4 can
// replace that implementation without touching the callers.
class DiffNetLegality
{
 public:
  virtual ~DiffNetLegality() = default;
  // True if a shape of `net` covering `box` on drt layer `layer` keeps all
  // diff-net distances to the routing space.
  virtual bool allowed(int layer, const Box& box, int net) const = 0;
  // Appends open boxes of reference points that a shape `model` (relative
  // to its reference point) of `net` on `layer` must avoid inside `window`:
  // other nets' shapes bloated by spacing, Minkowski difference with the
  // model. For a wire model [-h, h]^2 these are the centre-line obstacles;
  // a segment avoiding them keeps `allowed` for its swept box.
  virtual void blockages(int layer,
                         const Box& window,
                         int net,
                         const Box& model,
                         std::vector<Box>& out) const = 0;
  // True if shapes of two different nets at `a` and `b` on `layer` violate
  // a distance rule (conflicts between candidate access paths).
  virtual bool conflict(int layer, const Box& a, const Box& b) const = 0;
};

// One wiring layer of the off-track search.
struct OffTrackLayer
{
  BlockageGridLayer grid;
  // Open interiors the centre line must avoid.
  std::vector<Box> obstacles;
  // Planar arcs allowed on this layer (false: the layer is only passed by
  // vias, e.g. a pin layer outside the track graph).
  bool wires = true;
  // Cost per DBU in x and y (1, 1 = length, as in Theorem 3.2).
  int64_t cost_x = 1;
  int64_t cost_y = 1;
};

struct OffTrackSource
{
  int layer = 0;
  // All grid vertices in the closed box are sources.
  Box region;
};

struct OffTrackTarget
{
  int layer = 0;
  GridPoint p;
  int id = 0;
};

struct OffTrackRequest
{
  std::vector<OffTrackLayer> layers;
  // Grid vertices outside the closed window are not used.
  Box window;
  std::vector<OffTrackSource> sources;
  std::vector<OffTrackTarget> targets;
  int64_t via_cost = 0;
  // Via between layer z and z + 1 at p allowed (Sec. 3.8: "Vertices of two
  // adjacent wiring layers are connected if a via is allowed here").
  std::function<bool(int z, const GridPoint& p)> via_allowed;
  // Optional rule hook (Sec. 3.8: "More sophisticated rules like short
  // edge avoidance can be added ... by removing edges that induce a
  // violation"): return false to remove a planar arc of direction `dir`
  // (0 +x, 1 -x, 2 +y, 3 -y) from `from` to `to` entered with
  // `incoming` (same codes, 4 = no incoming direction).
  std::function<bool(int z,
                     const GridPoint& from,
                     const GridPoint& to,
                     int incoming,
                     int dir)>
      arc_filter;
  int64_t max_cost = std::numeric_limits<int64_t>::max() / 4;
  // Stop after this many distinct target vertices are settled (0 = all):
  // the catalogue keeps the nearest on-track points.
  int max_targets = 0;
  // Use the admissible future cost (A*); false = plain Dijkstra.
  bool future_cost = true;
};

struct OffTrackPath
{
  int target_id = 0;
  int64_t cost = 0;
  // Source to target, one node per bend and two per via.
  std::vector<RouteNode> nodes;
};

struct OffTrackResult
{
  // At most one (shortest) path per reached target, in increasing cost.
  std::vector<OffTrackPath> paths;
  int64_t grid_vertices = 0;  // blockage-grid vertices inside the window
  int64_t states = 0;         // vertices of G (5 per grid vertex)
  int64_t settled = 0;
  int64_t arcs = 0;
};

// The path-preserving digraph G (Sec. 3.8) on the blockage-grid vertices of
// all layers: per grid vertex one state per incoming direction (+x, -x, +y,
// -y) and one state "no incoming direction" for sources and via landings.
//  - straight arcs: to the neighbouring grid vertex in the incoming
//    direction;
//  - turn arcs: to the nearest grid vertex at distance >= tau in each
//    direction perpendicular to the incoming one (all four directions from
//    the "no incoming direction" state), so every segment, the first and
//    the last included, is at least tau long;
//  - via arcs: to the "no incoming direction" state of the same point on an
//    adjacent layer where `via_allowed`.
// An arc exists only if its segment avoids the open obstacle interiors.
// Dijkstra (with an admissible, consistent future cost) from all sources
// settles the targets; the paths are the shortest
// tau-feasible paths among those with grid vertices as segment ends, which
// by Theorem 3.2 are shortest tau-feasible paths.
class OffTrackSearch
{
 public:
  OffTrackResult run(const OffTrackRequest& request);

 private:
  struct LayerIndex
  {
    std::vector<GeoCoord> xs;
    std::vector<GeoCoord> ys;
    size_t base = 0;
    const std::vector<Box>* obstacles = nullptr;
    // Per row (y) / column (x), built on first use: the sorted open
    // intervals (lo, hi) of obstacle interiors crossing the line.
    std::vector<std::vector<std::pair<GeoCoord, GeoCoord>>> rows;
    std::vector<std::vector<std::pair<GeoCoord, GeoCoord>>> cols;
    std::vector<char> row_built;
    std::vector<char> col_built;
    // Index of each line in the layer above / below, or -1.
    std::vector<int> up_x;
    std::vector<int> up_y;
    std::vector<int> down_x;
    std::vector<int> down_y;
    int nx() const { return static_cast<int>(xs.size()); }
    int ny() const { return static_cast<int>(ys.size()); }
  };
  // True if the segment between grid lines ix0 < ix1 on row iy (resp.
  // iy0 < iy1 on column ix) avoids all open obstacle interiors.
  bool freeH(LayerIndex& l, int iy, int ix0, int ix1);
  bool freeV(LayerIndex& l, int ix, int iy0, int iy1);

  std::vector<LayerIndex> index_;
  // State labels, valid where stamp_ == epoch_ (buffers are reused across
  // runs without clearing).
  std::vector<int64_t> dist_;
  std::vector<int64_t> parent_;
  std::vector<uint32_t> stamp_;
  std::vector<int8_t> via_memo_;
  std::vector<uint32_t> via_stamp_;
  // Future cost per grid vertex, valid where future_stamp_ == epoch_.
  std::vector<int64_t> future_memo_;
  std::vector<uint32_t> future_stamp_;
  uint32_t epoch_ = 0;
};

}  // namespace drt::br
