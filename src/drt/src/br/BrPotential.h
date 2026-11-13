// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Goal-oriented potentials for the BonnRoute path search (Ahrens, Henke,
// Rabenstein, Vygen, "Faster Goal-Oriented Shortest Path Search for Bulk and
// Incremental Detailed Routing", arXiv:2111.06169v2).
//
// Every potential is pi(v) = dist_(G,c)(v, T) (or a lower bound of it for
// L1Potential) in the infinite grid G = Z x Z x {0..l-1} of Sec. 1.1, where
// T is a union of target rectangles, each with an additive terminal cost:
// pi(v) = min over targets R of (dist(v, R) + terminal cost of R). All
// arithmetic is exact int64; there is no floating point.
//
// Cost conventions (shared by all classes):
//  * Costs are per unit length (DBU); a rate of kInfCost forbids the
//    direction. Finite rates must be >= 0 and at most kMaxFiniteRate, and
//    coordinates must fit in Coord, so that finite path costs inside the
//    query domain stay far below kInfCost. Results >= kInfCost are reported
//    as kInfCost (unreachable).
//  * Layers are 0-based; via_up[z] is the via between z and z + 1.
//
// The classes are reusable: prepare() may be called many times and keeps the
// capacity of its buffers. query() performs no allocation, except that
// GeneralPotential in kTileFunctions mode appends the lazily built function
// list of a tile to a pool whose capacity is kept across prepare() calls.
// Instances are not thread-safe; use one per worker thread.

#pragma once

#include <cstdint>
#include <vector>

#include "br/BrCore.h"

namespace drt::br {

constexpr Cost kMaxFiniteRate = Cost{1} << 24;
// Upper bound on the number of layers (stack buffers in the queries).
constexpr int kMaxPotentialLayers = 64;

// Per-layer rates of the simple model (Sec. 1.1).
struct LayerCostRates
{
  std::vector<Cost> horizontal;  // cost per DBU along x on layer z
  std::vector<Cost> vertical;    // cost per DBU along y on layer z
  std::vector<Cost> via_up;      // via z -> z + 1, size numLayers() - 1

  int numLayers() const { return static_cast<int>(horizontal.size()); }
};

// Target rectangle {x in [xlo, xhi], y in [ylo, yhi], layer z}.
struct TargetRectangle
{
  Coord xlo = 0;
  Coord ylo = 0;
  Coord xhi = 0;
  Coord yhi = 0;
  int z = 0;
  Cost terminal_cost = 0;  // >= 0, added to the distance to this rectangle
};

// Cumulative counters for reports; reset() clears them.
struct PotentialCounters
{
  int64_t prepares = 0;
  int64_t queries = 0;
  // GeneralPotential preprocessing (Theorem 13).
  int64_t functions_offered = 0;   // affine functions offered to envelopes
  int64_t functions_inserted = 0;  // ... that were not dominated
  int64_t functions_removed = 0;   // removed later as dominated
  int64_t heap_pops = 0;           // functions processed (propagated)
  int64_t stale_heap_entries = 0;  // lazy heap entries skipped or re-keyed
  int64_t boundary_elements = 0;   // line elements (both directions)
  // GeneralPotential kTileFunctions mode (Sec. 5.1 variant).
  int64_t tile_lists_built = 0;
  int64_t tile_candidates = 0;  // 2-D affine functions generated
  int64_t tile_functions = 0;   // ... kept after domination pruning

  void reset() { *this = PotentialCounters{}; }
};

// l1-distance potential of Table 3's caption: distances in x and y scaled by
// the minimum horizontal / vertical rate over all layers, plus the via cost
// difference c(0, z') - c(0, z) of the cheapest via per layer pair. Feasible
// for every edge whose cost is at least the corresponding layer rate.
class L1Potential
{
 public:
  void prepare(const LayerCostRates& rates,
               const std::vector<TargetRectangle>& targets);
  Cost query(int z, Coord x, Coord y) const;
  const PotentialCounters& counters() const { return counters_; }
  void resetCounters() { counters_.reset(); }

 private:
  Cost viaLowerBound(int z0, int z1) const;

  Cost min_horizontal_ = kInfCost;
  Cost min_vertical_ = kInfCost;
  // Prefix sums of the finite via costs and prefix counts of forbidden vias.
  std::vector<Cost> via_prefix_;
  std::vector<int> via_forbidden_prefix_;
  std::vector<TargetRectangle> targets_;
  mutable PotentialCounters counters_;
};

// Exact distance in the simple model, Theorem 3. A path of Lemma 1 (vias,
// a run on layer a, vias, a run on layer b, vias) costs C + h * dx + v * dy
// with C its via cost and (h, v) the rates of its runs, so the distance
// between two layers is the lower envelope of these linear functions over
// all (a, b, run order). prepare() tabulates, per layer pair, the functions
// that no other function bounds from below everywhere (a handful: the rates
// take few values); the table is kept while the rates stay the same, so a
// query costs a few multiply-adds per target.
class SimplePotential
{
 public:
  void prepare(const LayerCostRates& rates,
               const std::vector<TargetRectangle>& targets);
  Cost query(int z, Coord x, Coord y) const;
  const PotentialCounters& counters() const { return counters_; }
  void resetCounters() { counters_.reset(); }
  // dist((x_target, y_target, z_target), (x, y, z)) for dx, dy >= 0: from
  // the table, and by the two O(l) sweeps of Theorem 3 (reference).
  Cost pointDistance(int target_z, int z, Cost dx, Cost dy) const;
  Cost pointDistanceSweep(int target_z, int z, Cost dx, Cost dy) const;
  // Linear functions in the table (all layer pairs).
  size_t tableFunctions() const { return functions_.size(); }

 private:
  Cost viaCost(int z0, int z1) const;
  void buildTable();

  // Lower envelope of layer pair (target_z, z): functions_[first_[k] ..
  // first_[k + 1]) with k = target_z * num_layers_ + z.
  struct LinearCost
  {
    Cost constant;
    Cost horizontal;
    Cost vertical;
  };
  std::vector<LinearCost> functions_;
  std::vector<int> first_;
  bool table_valid_ = false;

  int num_layers_ = 0;
  LayerCostRates rates_;
  // Minimum rates over the layers: the L1 bound used to order and skip
  // targets in query().
  Cost min_horizontal_ = kInfCost;
  Cost min_vertical_ = kInfCost;
  std::vector<Cost> via_prefix_;
  std::vector<int> via_forbidden_prefix_;
  std::vector<TargetRectangle> targets_;
  mutable PotentialCounters counters_;
};

// General model of Sec. 1.1: cuts xi (x_cuts) and upsilon (y_cuts) define
// closed tiles [xi_i, xi_{i+1}] x [upsilon_j, upsilon_{j+1}] per layer, with
// xi_0 = upsilon_0 = -infinity and xi_{p+1} = upsilon_{q+1} = +infinity. An
// edge lying in several closed tiles costs the minimum of their rates. Two
// (never three) consecutive equal cuts form a degenerate tile, which models a
// cheap line (Sec. 5.2); neutral degenerate tiles use kInfCost.
//
// Rates are indexed by index(z, i, j), i in [0, p], j in [0, q]. via_up of a
// tile on layer z is the via z -> z + 1 (unused on the top layer).
//
// The domain is the box in which the potential is exact. The effective
// domain is the bounding box of `domain`, all cuts and all targets; queries
// outside it are projected onto it, which still yields a feasible lower
// bound (the outer tiles are invariant along their unbounded direction).
struct TileCostModel
{
  int num_layers = 0;
  std::vector<Coord> x_cuts;  // xi_1 <= ... <= xi_p
  std::vector<Coord> y_cuts;  // upsilon_1 <= ... <= upsilon_q
  Coord domain_xlo = 0;
  Coord domain_ylo = 0;
  Coord domain_xhi = 0;
  Coord domain_yhi = 0;
  std::vector<Cost> horizontal;
  std::vector<Cost> vertical;
  std::vector<Cost> via_up;

  int numColumns() const { return static_cast<int>(x_cuts.size()) + 1; }
  int numRows() const { return static_cast<int>(y_cuts.size()) + 1; }
  size_t index(int z, int i, int j) const
  {
    return (static_cast<size_t>(z) * numColumns() + i) * numRows() + j;
  }
  // Sizes the rate arrays for the current cuts and fills them with `rate`.
  void resize(int layers, Cost rate);
  // Simple model as a general model without cuts (p = q = 0).
  void assignSimple(const LayerCostRates& rates,
                    Coord xlo,
                    Coord ylo,
                    Coord xhi,
                    Coord yhi);
};

// Per-layer minimum of the finite tile rates (the rates of the simple and
// L1 potentials that remain feasible for this model, e.g. with discounts).
LayerCostRates minimumLayerRates(const TileCostModel& model);

// Global routing cell grid as BonnRouteDR maps coordinates to Region tiles:
// tile tx covers [origin_x + tx * pitch_x, origin_x + (tx + 1) * pitch_x),
// and the first / last tile extend to the die boundary.
struct GCellGrid
{
  Coord origin_x = 0;
  Coord origin_y = 0;
  Coord pitch_x = 1;
  Coord pitch_y = 1;
  Coord die_xlo = 0;
  Coord die_ylo = 0;
  Coord die_xhi = 0;
  Coord die_yhi = 0;
};

// Reservation of the routed net (Sec. 5.2): the edges of an axis-parallel
// segment on layer z cost rate * numerator / denominator (rounded down)
// instead of the rate of the surrounding tiles, in the direction of the
// segment. Exact integer costs need rate * numerator divisible by
// denominator; scale all base rates by the denominator to guarantee it (for
// delta = 3/4: base rates times 4, reserved edges cost 3 * unscaled rate).
struct DiscountedWire
{
  int z = 0;
  Coord x1 = 0;
  Coord y1 = 0;
  Coord x2 = 0;
  Coord y2 = 0;
  int64_t numerator = 3;
  int64_t denominator = 4;
};

// Builds the general model of one search from the routing area: allowed
// Region tiles get the base rates, all other tiles kInfCost (closed tiles,
// so edges on the boundary of an allowed GCell are finite); a via z -> z + 1
// is finite where the GCell is allowed on both layers. Each discounted wire
// becomes a degenerate tile (Sec. 5.2). Adjacent columns / rows with equal
// rates are merged, which does not change any edge cost.
class TileCostModelBuilder
{
 public:
  void build(const Region& region,
             const GCellGrid& grid,
             const LayerCostRates& base,
             const std::vector<DiscountedWire>& wires,
             TileCostModel& model);

 private:
  struct Line
  {
    Coord position;
    bool doubled;  // carries a degenerate tile (reservation line)
  };
  void collectLines(std::vector<Line>& lines, std::vector<Coord>& cuts);
  void mergeEqualColumns(TileCostModel& model);
  void mergeEqualRows(TileCostModel& model);

  std::vector<Line> x_lines_;
  std::vector<Line> y_lines_;
  std::vector<DiscountedWire> horizontal_wires_;  // clipped, x1 < x2
  std::vector<DiscountedWire> vertical_wires_;    // clipped, y1 < y2
  std::vector<Cost> scratch_;
  std::vector<Coord> scratch_cuts_;
  std::vector<int> kept_indices_;
  std::vector<uint8_t> keep_;
};

// Exact distance in the general model, Sec. 4. prepare() refines the grid
// with respect to the targets (Corollary 15) and runs the boundary
// propagation of Theorem 13 twice: once for the horizontal grid lines and
// once for the vertical ones. Each grid line is split into line elements:
// the grid points (intersections with the perpendicular grid lines) and the
// open intervals between them. Every element keeps the lower envelope of
// affine functions (Lemma 9) over its integer points.
//
// Queries: binary search of the tile; on a grid line the element envelope is
// evaluated directly; in the open interior of a tile either
//  * kBoundary: the post-Theorem-13 query over the 4l boundary elements plus
//    via-only paths, or
//  * kTileFunctions (default): the Sec. 5.1 variant, a per (tile, layer)
//    list of 2-D affine functions c + a*x + b*y built lazily from the same
//    boundary envelopes and pruned by provable domination (pairwise >= at
//    the four corners of the tile interior).
class GeneralPotential
{
 public:
  enum class QueryMode : uint8_t
  {
    kBoundary,
    kTileFunctions
  };

  void prepare(const TileCostModel& model,
               const std::vector<TargetRectangle>& targets);
  Cost query(int z, Coord x, Coord y) const;
  void setQueryMode(QueryMode mode) { mode_ = mode; }
  QueryMode queryMode() const { return mode_; }
  const PotentialCounters& counters() const { return counters_; }
  void resetCounters() { counters_.reset(); }

 private:
  enum RateKind : uint8_t
  {
    kHorizontalRate = 0,
    kVerticalRate = 1,
    kViaRate = 2
  };
  // One coordinate axis after refinement. Slot 2k + 1 is the grid line at
  // positions[k]; slot 2k is the open interval between positions[k - 1] and
  // positions[k] (k = 0 and k = P are the outer intervals up to the domain
  // boundary). Each slot records the range of model tiles containing it.
  struct Axis
  {
    std::vector<Coord> positions;
    std::vector<int> first_tile;
    std::vector<int> last_tile;
    Coord domain_lo = 0;
    Coord domain_hi = 0;

    int numPositions() const { return static_cast<int>(positions.size()); }
    int numSlots() const { return 2 * numPositions() + 1; }
    int locate(Coord c) const;
    Coord slotLo(int slot) const;
    Coord slotHi(int slot) const;
    void build(const std::vector<Coord>& cuts,
               const std::vector<Coord>& extra,
               Coord lo,
               Coord hi);
  };
  // Affine function offset + slope * u over the integer points of a line
  // element, u the coordinate along the line. [owned_lo, owned_hi] is where
  // it attains the envelope minimum (Lemma 9).
  struct LineFunction
  {
    Cost offset;
    Cost slope;
    Cost key;
    Coord owned_lo;
    Coord owned_hi;
    int32_t next;
    int32_t element;
    bool alive;
    bool processed;
    Cost at(Coord u) const { return offset + slope * u; }
  };
  struct HeapEntry
  {
    Cost key;
    int32_t function;
    bool operator>(const HeapEntry& o) const
    {
      return key != o.key ? key > o.key : function > o.function;
    }
  };
  // Grid lines of one direction (Theorem 13 runs once per direction). For
  // horizontal lines, the line index m is the y position index and the
  // along slot is an x slot; vertical lines swap the axes.
  struct LineRun
  {
    bool horizontal_lines = true;
    int num_lines = 0;
    int num_slots = 0;
    std::vector<int32_t> head;  // per element: first function, -1 if none
    size_t element(int z, int line, int slot) const
    {
      return (static_cast<size_t>(z) * num_lines + line) * num_slots + slot;
    }
  };
  struct TileFunction
  {
    Cost offset;
    Cost slope_x;
    Cost slope_y;
  };
  struct TileCandidate
  {
    TileFunction function;
    Cost corner[4];
  };
  struct TileListEntry
  {
    uint32_t stamp = 0;
    int32_t begin = 0;
    int32_t count = 0;
  };

  Cost rate(RateKind kind, int x_slot, int y_slot, int z) const;
  Cost runRate(const LineRun& run,
               RateKind kind,
               int along_slot,
               int line_slot,
               int z) const;
  const Axis& alongAxis(const LineRun& run) const
  {
    return run.horizontal_lines ? x_axis_ : y_axis_;
  }
  const Axis& lineAxis(const LineRun& run) const
  {
    return run.horizontal_lines ? y_axis_ : x_axis_;
  }
  void runPropagation(LineRun& run);
  void offer(LineRun& run,
             size_t element,
             Coord lo,
             Coord hi,
             Cost offset,
             Cost slope);
  void process(LineRun& run, int32_t function);
  Cost evaluate(const LineRun& run, int z, int line, int slot, Coord u) const;
  Cost tileInteriorBoundary(int z,
                            int x_slot,
                            int y_slot,
                            Coord x,
                            Coord y) const;
  Cost tileInteriorFunctions(int z,
                             int x_slot,
                             int y_slot,
                             Coord x,
                             Coord y) const;
  void buildTileList(int z, int x_slot, int y_slot, TileListEntry& entry) const;
  void addTileCandidate(const TileFunction& f,
                        Coord xa,
                        Coord ya,
                        Coord xb,
                        Coord yb) const;

  QueryMode mode_ = QueryMode::kTileFunctions;
  int num_layers_ = 0;
  bool has_targets_ = false;
  // Copy of the model rates (queries read them in kBoundary mode and when
  // building tile lists lazily).
  int model_columns_ = 1;
  int model_rows_ = 1;
  std::vector<Cost> rates_[3];
  Axis x_axis_;
  Axis y_axis_;
  std::vector<TargetRectangle> targets_;
  LineRun horizontal_run_;
  LineRun vertical_run_;
  std::vector<LineFunction> functions_;
  std::vector<HeapEntry> heap_;
  std::vector<Coord> scratch_positions_;
  // Lazy tile lists (kTileFunctions).
  mutable std::vector<TileListEntry> tile_lists_;
  mutable std::vector<TileFunction> tile_functions_;
  mutable std::vector<TileCandidate> tile_scratch_;
  mutable uint32_t tile_epoch_ = 0;
  mutable PotentialCounters counters_;
};

}  // namespace drt::br
