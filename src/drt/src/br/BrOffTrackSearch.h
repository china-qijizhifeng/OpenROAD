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
#include <utility>
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
  // On-track-only layer (Sec. 3.8 "removing edges that induce a
  // violation"): with `tracks` set, planar arcs run only in the preferred
  // direction (`horizontal`) and start on one of the sorted track lines
  // (y = track if horizontal, else x = track). The search evaluates this
  // per grid line once; it equals an arc_filter with the same rule.
  const std::vector<GeoCoord>* tracks = nullptr;
  bool horizontal = true;
};

struct OffTrackSource
{
  int layer = 0;
  // All grid vertices in the closed box are sources.
  Box region;
  // Initial distance (e.g. a penalty for an access point pin access did
  // not choose); the future cost stays admissible.
  int64_t cost = 0;
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
//
// The result (paths, their order and every counter) is a function of the
// request: the heap holds (f, state) pairs ordered by f only, and every
// push and pop happens in the order of the reference formulation above
// (BonnRouteOffTrackTest keeps that formulation and compares results on
// random requests). A window with 2^32 / 5 or more grid vertices is not
// searched (no paths; grid_vertices reports its size).
class OffTrackSearch
{
 public:
  OffTrackResult run(const OffTrackRequest& request);

 private:
  struct LayerIndex
  {
    std::vector<GeoCoord> xs;
    std::vector<GeoCoord> ys;
    uint32_t base = 0;
    int nx = 0;
    int ny = 0;
    const std::vector<Box>* obstacles = nullptr;
    // Obstacles in increasing ylo / xlo order: the candidates of a row y
    // (column x) are the prefix with ylo < y (xlo < x).
    std::vector<Box> by_ylo;
    std::vector<Box> by_xlo;
    // Per row (y) / column (x), built on first use: the sorted, merged
    // open intervals (lo, hi) of obstacle interiors crossing the line, as
    // [begin, end) of `intervals`; begin = -1 until built.
    std::vector<std::pair<int, int>> row_span;
    std::vector<std::pair<int, int>> col_span;
    std::vector<std::pair<GeoCoord, GeoCoord>> intervals;
    // Index of each line in the layer above / below, or -1.
    std::vector<int> up_x;
    std::vector<int> up_y;
    std::vector<int> down_x;
    std::vector<int> down_y;
    // Turn arcs: the nearest line at distance >= tau above (`_plus`, nx /
    // ny if none) and below (`_minus`, -1 if none).
    std::vector<int> turn_x_plus;
    std::vector<int> turn_x_minus;
    std::vector<int> turn_y_plus;
    std::vector<int> turn_y_minus;
    // On-track-only layers: lines on which planar arcs may start.
    std::vector<char> row_on_track;
    std::vector<char> col_on_track;
    // Future cost terms per target layer t: l1 distance of each line to
    // the nearest target coordinate of layer t (index t * nx + ix).
    std::vector<int64_t> near_x;
    std::vector<int64_t> near_y;
  };
  // True if the segment between grid lines ix0 < ix1 on row iy (resp.
  // iy0 < iy1 on column ix) avoids all open obstacle interiors.
  bool freeH(LayerIndex& l, int iy, int ix0, int ix1);
  bool freeV(LayerIndex& l, int ix, int iy0, int iy1);

  // Label of a state, valid where stamp == epoch_ (buffers are reused
  // across runs without clearing).
  struct Label
  {
    int64_t dist;
    uint32_t parent;  // kNoParent for sources
    uint32_t stamp;
  };
  struct HeapItem
  {
    int64_t f;
    uint32_t state;
  };
  std::vector<LayerIndex> index_;
  std::vector<Label> labels_;
  std::vector<HeapItem> heap_;
  // Per grid vertex: via memo (stamp, allowed) and the first index of the
  // vertex in the sorted target list (valid where target_stamp_ == epoch_).
  std::vector<uint32_t> via_stamp_;
  std::vector<int8_t> via_memo_;
  std::vector<uint32_t> target_stamp_;
  std::vector<int32_t> target_first_;
  uint32_t epoch_ = 0;
};

}  // namespace drt::br
