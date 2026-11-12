// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// BonnRoute detailed routing core (Gester et al., TODAES 18(2), 2013,
// Sec. 4.1, 4.2): the on-track vertex and interval path searches with
// future cost. They read legality from the fast grid (Sec. 3.6, see
// BrFastGrid.h), which forwards marked edges and other wire types to the
// distance rule checker. Independent of the drt database so it can be
// tested in isolation; BonnRouteDR adapts it to frDesign.

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "br/BrFastGrid.h"
#include "br/BrTrackGraph.h"

namespace drt::br {

// A track-graph vertex by layer, track and position.
struct At
{
  int z;
  int i;
  int j;
};

// Routing area R of one search: allowed global routing tiles per layer.
struct Region
{
  int tiles_x = 1;
  int tiles_y = 1;
  std::vector<std::vector<uint8_t>> allowed;  // [z][tx * tiles_y + ty]

  void init(int num_layers, int tx, int ty)
  {
    tiles_x = tx;
    tiles_y = ty;
    allowed.assign(num_layers, std::vector<uint8_t>(size_t(tx) * ty, 0));
  }
  std::vector<std::pair<int, size_t>> touched;
  // Clears the allowed tiles in O(#allowed).
  void reset()
  {
    for (const auto& [z, index] : touched) {
      allowed[z][index] = 0;
    }
    touched.clear();
    bx_lo = by_lo = std::numeric_limits<int>::max();
    bx_hi = by_hi = -1;
  }
  bool allowedTile(int z, int tx, int ty) const
  {
    return z >= 0 && z < static_cast<int>(allowed.size()) && tx >= 0
           && tx < tiles_x && ty >= 0 && ty < tiles_y
           && allowed[z][size_t(tx) * tiles_y + ty] != 0;
  }
  // Bounding box of all allowed tiles.
  int bx_lo = std::numeric_limits<int>::max();
  int by_lo = std::numeric_limits<int>::max();
  int bx_hi = -1;
  int by_hi = -1;

  void allow(int z, int tx, int ty)
  {
    if (z >= 0 && z < static_cast<int>(allowed.size()) && tx >= 0
        && tx < tiles_x && ty >= 0 && ty < tiles_y) {
      uint8_t& cell = allowed[z][size_t(tx) * tiles_y + ty];
      if (!cell) {
        cell = 1;
        touched.emplace_back(z, size_t(tx) * tiles_y + ty);
      }
      bx_lo = std::min(bx_lo, tx);
      by_lo = std::min(by_lo, ty);
      bx_hi = std::max(bx_hi, tx);
      by_hi = std::max(by_hi, ty);
    }
  }
  bool contains(const TrackLayer& l, int z, int i, int j) const
  {
    const int t = l.track_tile[i];
    const int p = l.pos_tile[j];
    const int tx = l.horizontal ? p : t;
    const int ty = l.horizontal ? t : p;
    return allowed[z][size_t(tx) * tiles_y + ty] != 0;
  }
};

// Target cluster for the future cost: a rectangle on a layer range.
struct TargetBox
{
  Coord xlo, ylo, xhi, yhi;
  int zlo, zhi;
  // Smallest terminal cost of the targets in the box; added to the bound.
  Cost cost = 0;
};

// Lower bound on the cost from a location to the targets of a search
// (future cost of Sec. 4.1; potentials of arXiv:2111.06169). It must be
// feasible for the search graph: pi(v) <= c(v, w) + pi(w) on every edge.
class FutureCostOracle
{
 public:
  virtual ~FutureCostOracle() = default;
  virtual Cost query(int z, Coord x, Coord y) const = 0;
};

// Same-net rules between the vias of the searching net (Lef58Rules.md).
// A via is given by its lower vertex (zv, i, j) on via layer zv (between
// track layers zv and zv + 1). The vertex search refuses a via that
// allowed() rejects (against the net's cuts outside the search) or that
// pairAllowed() rejects together with the via by which the path entered
// the current layer, when both are on via layer zv: the path turns back on
// `turn_layer` (zv or zv + 1), and the pads of the two vias on the other
// layer are not connected by it.
class SameNetVias
{
 public:
  virtual ~SameNetVias() = default;
  virtual bool allowed(int zv, int i, int j) const = 0;
  virtual bool pairAllowed(int zv, int i, int j, int i0, int j0, int turn_layer)
      const = 0;
  // Via layers with rules (the search checks only their vias).
  bool active(const int zv) const
  {
    return zv >= 0 && zv < static_cast<int>(layers.size()) && layers[zv] != 0;
  }
  std::vector<uint8_t> layers;
  // The search looks for the entering via within this distance (DBU, in x
  // and y) of the vertex.
  Coord reach = 0;
};

// Preferred-direction edges (j, j + 1) for jlo <= j < jhi on track i of
// layer z that the searching net may use at a discount (reservations of
// its input wiring, arXiv:2111.06169 Sec. 5.2).
struct DiscountRun
{
  int z, i, jlo, jhi;
};

struct SearchRequest
{
  int32_t net = -1;
  const Region* region = nullptr;
  // Source vertices with their initial cost (0 on the routed tree, stub
  // cost at pin attachments).
  std::vector<std::pair<size_t, Cost>> sources;
  // Target vertices: (vertex, terminal cost, target id).
  struct Target
  {
    size_t vertex;
    Cost cost;
    int id;
  };
  std::vector<Target> targets;
  std::vector<TargetBox> target_boxes;
  // Wire type of the connection (Sec. 4.1: "a wiretype w"); wire types
  // outside the fast grid are answered by the distance rule checker.
  int wire_type = 0;
  // Rip-up mode (Sec. 4.2): vertices and edges whose shapes become legal
  // after ripping up shapes of ripup level <= max_ripup_level are usable at
  // ripup_cost * (1 + history).
  bool ripup = false;
  Cost ripup_cost = 0;
  RipupLevel max_ripup_level = kMaxRemovableLevel;
  // Fast-grid words near the searching net's own shapes with those shapes
  // removed (Sec. 4.4); nullptr: none.
  const FastGridOverlay* overlay = nullptr;
  // Per-vertex label budget of the search (0: unlimited); a search that
  // exceeds it stops without a path and reports `aborted`.
  size_t label_limit = 0;
  bool use_intervals = true;
  // Interval property (ii) of Sec. 4.1 for vias as well as jogs. Without
  // it, intervals ignore via availability and relax vias per vertex.
  bool interval_via_property = true;
  // Among intervals at the same key, process the one whose J_I(delta)
  // reaches the smallest future cost first (the vertex search's deeper-
  // first rule) instead of the largest (Algorithm 4, line 12). Both are
  // exact; this one reaches the targets early in the last key level.
  bool interval_nearest_first = false;
  // pi_H of Sec. 4.1 with lb_via(z) = cost of the cheapest via stack from z
  // to a target layer (exact); false keeps dz * (cheapest via).
  bool exact_future_cost = true;
  // Optional future cost; nullptr selects the l1 bound pi_H on the target
  // boxes. The interval search uses an oracle as pi with property (i)
  // checked between neighbours (it needs wire_unit 1); only the vertex
  // search supports discounts.
  const FutureCostOracle* future_cost = nullptr;
  std::vector<DiscountRun> discounts;
  Cost discount_rate = 0;  // cost per DBU of a discounted edge
  // Wire spreading (Gester et al. 2013, Sec. 4.2): in tiles where global
  // routing left room ([z][tx * tiles_y + ty] != 0), entering a wire vertex
  // next to an occupied track costs spread_cost, which keeps the tracks
  // beside routed wires free where space allows. In the interval search
  // such vertices are singletons and the cost enters the wire relaxations.
  const std::vector<std::vector<uint8_t>>* spread_tiles = nullptr;
  Cost spread_cost = 0;
  // Same-net via rules (vertex search only; nullptr: none).
  const SameNetVias* same_net_vias = nullptr;
  // Vertex search: in tiles where global routing left most tracks of the
  // layer free ([z][tx * tiles_y + ty] != 0), jogs cost the layer's
  // jog_factor_low instead of jog_factor (a wrong-way wire blocks the
  // crossing tracks it uses, which matters where they are in demand).
  const std::vector<std::vector<uint8_t>>* jog_low_tiles = nullptr;
};

struct SearchResult
{
  bool found = false;
  Cost cost = kInfCost;
  int target_id = -1;
  size_t target_vertex = 0;
  std::vector<size_t> path;  // from target back to a source
  // Statistics.
  int64_t heap_pops = 0;
  int64_t heap_pushes = 0;
  int64_t settled = 0;
  // Legality queries answered by the fast grid / the distance rule checker.
  int64_t fast_grid_queries = 0;
  int64_t checker_queries = 0;
  // Run labels popped, and those whose continuation (the run label they
  // pushed one position further along the track) is the very next pop
  // (vertex search).
  int64_t run_label_pops = 0;
  int64_t run_continuations = 0;
  bool aborted = false;  // label budget exceeded
};

// pi_H of Sec. 4.1 over the target boxes of a search, in flat arrays: per
// box its bounds, per (layer, box) the via bound plus terminal cost. The
// value at (z, x, y) is min over boxes of l1 distance * rate + constant.
class TargetBoxBound
{
 public:
  // via_prefix[z] = cost of the vias 0 -> z; min_via_cost bounds a via
  // when exact_future_cost is off.
  void setup(const SearchRequest& request,
             const std::vector<Cost>& via_prefix,
             Cost min_via_cost);
  bool empty() const { return xlo_.empty(); }
  size_t size() const { return xlo_.size(); }
  // Via bound plus terminal cost of box b from layer z.
  Cost constant(const int z, const size_t b) const
  {
    return constant_[static_cast<size_t>(z) * xlo_.size() + b];
  }
  Cost at(const int z, const Coord x, const Coord y, const Cost rate) const
  {
    const size_t n = xlo_.size();
    const Cost* constant = constant_.data() + static_cast<size_t>(z) * n;
    Cost best = kInfCost;
    for (size_t b = 0; b < n; b++) {
      const Cost dx = std::max({0, xlo_[b] - x, x - xhi_[b]});
      const Cost dy = std::max({0, ylo_[b] - y, y - yhi_[b]});
      best = std::min(best, (dx + dy) * rate + constant[b]);
    }
    return best;
  }

 private:
  std::vector<Coord> xlo_;
  std::vector<Coord> xhi_;
  std::vector<Coord> ylo_;
  std::vector<Coord> yhi_;
  std::vector<Cost> constant_;
};

// Per-vertex labels of a search window, stored in pages that a search
// allocates when it touches them: memory follows the explored part of the
// window, not its bounding box (a net spanning the die has a window of
// most of the track graph). Pages are reused by later searches. With
// kClear, a page is reset to Label{} when a search first touches it, so
// labels need no epoch of their own; otherwise stale labels carry older
// epochs. A search may be limited to a number of pages; beyond it, labels
// go to a scratch page and overflow() is set (the search then stops).
// Between searches at most kKeptBytes of pages stay allocated.
template <typename Label, bool kClear = false>
class PagedLabels
{
 public:
  static constexpr int kPageBits = 10;
  static constexpr size_t kPageSize = size_t{1} << kPageBits;
  static constexpr size_t kKeptBytes = size_t{16} << 20;

  // Starts a search over `size` window vertices with at most `limit`
  // labels (0: unlimited); returns its epoch.
  uint32_t reset(const size_t size, const size_t limit = 0)
  {
    const size_t kept = kKeptBytes / (kPageSize * sizeof(Label));
    if (pool_.size() > kept) {
      pool_.resize(kept);
    }
    limit_pages_ = limit == 0 ? 0 : (limit + kPageSize - 1) / kPageSize;
    overflow_ = false;
    const size_t pages = (size >> kPageBits) + 1;
    if (slots_.size() < pages) {
      slots_.resize(pages, Slot{});
    }
    used_ = 0;
    if (++epoch_ == 0) {
      std::fill(slots_.begin(), slots_.end(), Slot{});
      for (auto& page : pool_) {
        std::fill(page.get(), page.get() + kPageSize, Label{});
      }
      epoch_ = 1;
    }
    return epoch_;
  }
  Label& operator[](const size_t index)
  {
    const Slot& slot = slots_[index >> kPageBits];
    if (slot.epoch == epoch_) {
      return slot.page[index & (kPageSize - 1)];
    }
    return touch(index);
  }
  bool overflow() const { return overflow_; }
  size_t bytes() const
  {
    return slots_.capacity() * sizeof(Slot)
           + pool_.size() * kPageSize * sizeof(Label);
  }
  // Largest page count a search has used.
  size_t peakPages() const { return peak_pages_; }

 private:
  // A page of the window: its storage while `epoch` is the current one.
  struct Slot
  {
    Label* page = nullptr;
    uint32_t epoch = 0;
  };
  // First touch of a page in this search.
  [[gnu::noinline]] Label& touch(const size_t index)
  {
    Slot& slot = slots_[index >> kPageBits];
    if (limit_pages_ != 0 && used_ >= limit_pages_) {
      overflow_ = true;
      if (!scratch_) {
        scratch_ = std::make_unique<Label[]>(1);
      }
      scratch_[0] = Label{};
      return scratch_[0];
    }
    if (used_ == pool_.size()) {
      pool_.push_back(std::make_unique<Label[]>(kPageSize));
    }
    slot.page = pool_[used_++].get();
    slot.epoch = epoch_;
    peak_pages_ = std::max(peak_pages_, used_);
    if constexpr (kClear) {
      std::fill(slot.page, slot.page + kPageSize, Label{});
    }
    return slot.page[index & (kPageSize - 1)];
  }

  std::vector<Slot> slots_;
  std::vector<std::unique_ptr<Label[]>> pool_;
  std::unique_ptr<Label[]> scratch_;
  size_t used_ = 0;
  size_t limit_pages_ = 0;
  size_t peak_pages_ = 0;
  bool overflow_ = false;
  uint32_t epoch_ = 0;
};

// Vertex search with future cost on the track graph (use_intervals = 0).
class PathSearch
{
 public:
  PathSearch(const TrackGraph& graph, const FastGrid& grid)
      : graph_(graph), grid_(grid)
  {
  }
  SearchResult run(const SearchRequest& request);
  size_t scratchBytes() const { return labels_.bytes(); }
  size_t peakScratchBytes() const
  {
    return labels_.peakPages() * Labels::kPageSize * sizeof(Label);
  }

 private:
  enum Parent : uint8_t
  {
    kSource = 0,
    kFromLower = 1,       // along: came from j - 1
    kFromHigher = 2,      // along: came from j + 1
    kFromTrackBelow = 3,  // jog: came from i - 1
    kFromTrackAbove = 4,  // jog: came from i + 1
    kFromLayerBelow = 5,  // via from z - 1
    kFromLayerAbove = 6,  // via from z + 1
  };
  // A heap entry at vertex (z, i, j). Ties of key and g are broken by the
  // vertex in (z, i, j) order, which is the order of the vertex ids.
  struct Item
  {
    Cost key;
    Cost g;
    int32_t i;
    int32_t j;
    int32_t aux;
    int8_t z;
    int8_t dir;  // 0 vertex label, +1 / -1 ray, 2 sink
    bool operator>(const Item& o) const
    {
      if (key != o.key) {
        return key > o.key;
      }
      if (g != o.g) {
        return g < o.g;  // deeper first on ties
      }
      if (z != o.z) {
        return z > o.z;
      }
      if (i != o.i) {
        return i > o.i;
      }
      return j > o.j;
    }
  };

  Cost futureCost(int z, Coord x, Coord y) const;
  // FastGridReader::level, inlined: the run cursors answer most queries
  // without a call.
  RipupLevel level(const At a, Kind k) const
  {
    return detail::levelOf(
        grid_,
        a.z,
        a.i,
        a.j,
        k,
        req_->wire_type,
        reader_.ignoredOwner(),
        &counters_,
        [this](int z, int i, int j) { return reader_.wiringWord(z, i, j); },
        [this](int z, int i, int c) { return reader_.viaWord(z, i, c); });
  }
  // A vertex or edge of level l may be entered: legal, or removable in
  // rip-up mode.
  bool usable(const RipupLevel l) const
  {
    return l == kLegal || (req_->ripup && l <= req_->max_ripup_level);
  }
  // Rip-up penalty of entering a of level l.
  Cost penalty(At a, RipupLevel l) const;
  bool enter(const At a, const Kind k) const { return usable(level(a, k)); }
  // 4-ary min-heap on heap_. Item::operator> is a strict total order (no
  // two queued items share key, g and vertex), so the pop sequence is that
  // of any other heap.
  void push(const Item& item);
  Item pop();
  bool relax(int z, int i, int j, Cost g, Parent parent, int8_t dir);
  // Cost of the preferred-direction edge between positions j and j2,
  // including the wire spreading cost of entering j2.
  Cost wireCost(int z, int i, int j, int j2) const;
  Cost spreadCost(int z, int i, int j) const;
  Cost jogFactor(int z, int i, int j) const;
  void settle(int z, int i, int j, Cost g);
  // The via by which the path to the settled vertex (z, i, j) entered layer
  // z, as via layer and lower vertex; false if the path starts on the layer
  // or the via is farther than `reach` in x or y.
  bool enteringVia(int z, int i, int j, Coord reach, At& via) const;

  const TrackGraph& graph_;
  const FastGrid& grid_;
  const SearchRequest* req_ = nullptr;
  SearchResult* res_ = nullptr;
  mutable LegalityCounters counters_;
  mutable FastGridReader reader_;
  Cost min_via_cost_ = 0;
  std::vector<Cost> via_prefix_;  // via_prefix_[z] = cost of vias 0 -> z
  Cost min_wire_rate_ = 1;
  TargetBoxBound boxes_;

  // Search window: per layer track range [i0, i0 + ni) and position range
  // [j0, j0 + nj) covering the region; local index = wbase + di * nj + dj.
  struct Window
  {
    int i0 = 0, ni = 0, j0 = 0, nj = 0;
    size_t wbase = 0;
  };
  bool local(int z, int i, int j, size_t& index) const
  {
    const Window& w = window_[z];
    const int di = i - w.i0;
    const int dj = j - w.j0;
    if (di < 0 || di >= w.ni || dj < 0 || dj >= w.nj) {
      return false;
    }
    index = w.wbase + static_cast<size_t>(di) * w.nj + dj;
    return true;
  }
  void setupWindow(const SearchRequest& request);
  std::vector<Window> window_;

  // Labels of one search (pages reset on first touch, see PagedLabels):
  // dist = kInfCost until the vertex is labelled; marks tell sources,
  // targets (target_index) and discounted edges (j, j + 1) at j.
  enum Mark : uint8_t
  {
    kMarkSource = 1,
    kMarkTarget = 2,
    kMarkDiscount = 4,
  };
  struct Label
  {
    Cost dist = kInfCost;
    int32_t target_index = 0;
    uint8_t parent = 0;
    uint8_t marks = 0;
  };
  using Labels = PagedLabels<Label, /*kClear=*/true>;
  mutable Labels labels_;
  uint32_t epoch_ = 0;
  std::vector<Item> heap_;
};

// Interval-based shortest path search (Sec. 4.1, Algorithm 4; Hetzel 1998,
// Peyer et al. 2009). Vertices of a track are grouped on the fly into
// maximal intervals I on which (i) the future cost changes by exactly the
// distance between vertices and (ii) the jog edges are the same. Each
// interval keeps a set L(I) of non-redundant labels (v, delta) with reduced
// distance d(v') = delta + |v - v'| - pi(v) + pi(v'); intervals are
// processed at their frontier J_I(delta) in increasing delta. Via edges are
// relaxed per vertex of J_I(delta): their positions interleave on layers
// whose neighbours have different pitches, so requiring (ii) for vias would
// reduce most intervals to single vertices. Target vertices and vertices
// with rip-up penalties form single-vertex intervals.
class IntervalSearch
{
 public:
  IntervalSearch(const TrackGraph& graph, const FastGrid& grid)
      : graph_(graph), grid_(grid)
  {
  }
  SearchResult run(const SearchRequest& request);
  size_t scratchBytes() const { return labels_.bytes(); }
  size_t peakScratchBytes() const
  {
    return labels_.peakPages() * Labels::kPageSize * sizeof(VertexLabel);
  }
  int64_t intervals() const { return num_intervals_; }
  int64_t labels() const { return num_labels_; }
  // Why intervals end (Sec. 4.1 properties), counted per interval end.
  enum IntervalEnd : uint8_t
  {
    kEndTrack = 0,        // end of the track or of the search window
    kEndStartVertex,      // the start vertex cannot join any interval
    kEndNeighbour,        // neighbour unusable, target, source or outside R
    kEndWireEdge,         // preferred-direction edge unusable
    kEndJogProperty,      // property (ii) differs for jogs
    kEndViaProperty,      // property (ii) differs for vias only
    kEndFutureCostSlope,  // property (i) fails
    kNumIntervalEnds
  };
  int64_t intervalEnds(IntervalEnd reason) const { return ends_[reason]; }
  int64_t intervalVertices() const { return interval_vertices_; }

 private:
  struct Label
  {
    int k;       // position index within the interval (0 = target end)
    Cost delta;  // reduced distance at k
    int32_t pred_interval;  // -1 for a source label
    int pred_k;
    int64_t order;  // preserve the oldest label when envelope values tie
  };
  struct Interval
  {
    int z, i;
    int ja, jb;  // position range, inclusive
    int dir;     // +1: k = j - ja (target end at ja); -1: k = jb - j
    Cost future_at_ja;
    std::vector<Label> labels;  // ordered by k
    // Prefix minima of delta - 2*x and suffix minima of delta. They query
    // the lower envelope in O(log labels) without scanning every label.
    std::vector<int> prefix_best;
    std::vector<int> suffix_best;
    int frontier = 0;  // vertices with k < frontier are settled
    int n() const { return jb - ja + 1; }
    int jOf(int k) const { return dir > 0 ? ja + k : jb - k; }
    int kOf(int j) const { return dir > 0 ? j - ja : jb - j; }
    Cost futureCost(const TrackLayer& layer, int j) const
    {
      return future_at_ja
             + dir
                   * (static_cast<Cost>(layer.positions[j])
                      - layer.positions[ja]);
    }
  };
  struct Event
  {
    Cost key;
    Cost max_future;
    int32_t interval;  // -1: sink (target reached), kPendingEvent
    int32_t aux;       // frontier at enqueue time, target index, pending index
    bool operator>(const Event& o) const
    {
      if (key != o.key) {
        return key > o.key;
      }
      // Algorithm 4, line 12: among intervals at delta, choose the one
      // whose settled run contains the largest future cost (max_future is
      // its negation with interval_nearest_first). Sinks come first.
      if (max_future != o.max_future) {
        return max_future < o.max_future;
      }
      return interval > o.interval;
    }
  };

  Cost futureCost(int z, int i, int j) const;
  // Wire spreading cost of entering wire vertex (z, i, j) along its track
  // (as PathSearch::spreadCost).
  Cost spreadCost(int z, int i, int j) const;
  // FastGridReader::level, inlined: the run cursors answer most queries
  // of a track scan without a call.
  RipupLevel level(const At a, Kind k) const
  {
    return detail::levelOf(
        grid_,
        a.z,
        a.i,
        a.j,
        k,
        req_->wire_type,
        reader_.ignoredOwner(),
        &counters_,
        [this](int z, int i, int j) { return reader_.wiringWord(z, i, j); },
        [this](int z, int i, int c) { return reader_.viaWord(z, i, c); });
  }
  bool usableZero(At v, Kind k) const;
  bool enter(At v, Kind k) const;
  Cost penalty(At v, Kind k) const;
  bool inWindow(int z, int i, int j, size_t& li) const;
  int intervalOf(int z, int i, int j);
  int buildInterval(int z, int i, int j);
  // What the interval search needs to know about a vertex, computed once
  // per search into its vertex label: its jog and via transitions
  // (property (ii)) and the levels its interval membership depends on.
  enum : uint16_t
  {
    kTransJogBelow = 1,
    kTransJogAbove = 2,
    kTransViaBelow = 4,
    kTransViaAbove = 8,
    kTransPenalty = 16,    // an existing transition carries a rip-up penalty
    kInfoWireLegal = 32,   // kWire
    kInfoJogLegal = 64,    // kJog
    kInfoEdgeLegal = 128,  // kWireEdge to position j + 1
    kInfoVias = 256,       // the via bits are computed
    kInfoSpread = 512,     // entering the vertex costs the spreading cost
  };
  struct VertexLabel;
  uint16_t info(int z, int i, int j, VertexLabel& label) const;
  uint16_t viaInfo(int z, int i, int j, VertexLabel& label) const;
  // The fast-grid words of the vertex and its neighbours give all levels
  // at once; wire types outside the fast grid and marked edges go through
  // levelOf to the checker.
  uint16_t computeInfo(int z, int i, int j, bool vias) const;
  uint16_t computeVias(int z, int i, int j, uint64_t own) const;
  // Via transitions are part of the interval property, or carry rip-up
  // penalties that make vertices singletons: needed to build intervals.
  bool vias_at_build_ = true;
  Cost value(const Interval& iv, int k) const;
  const Label* argLabel(const Interval& iv, int k) const;
  void rebuildEnvelope(Interval& iv);
  int frontierEnd(const Interval& iv, Cost delta) const;
  void addLabel(int interval, int k, Cost delta, int pred_interval, int pred_k);
  // Settles J_I(delta) = [frontier, g) (g < 0: computed) and relaxes it.
  void process(int interval, Cost delta, int g);
  void pushEvent(int interval);

  const TrackGraph& graph_;
  const FastGrid& grid_;
  const SearchRequest* req_ = nullptr;
  SearchResult* res_ = nullptr;
  mutable LegalityCounters counters_;
  mutable FastGridReader reader_;
  Cost min_via_cost_ = 0;
  std::vector<Cost> via_prefix_;  // via_prefix_[z] = cost of vias 0 -> z
  TargetBoxBound boxes_;

  struct Window
  {
    int i0 = 0, ni = 0, j0 = 0, nj = 0;
    size_t wbase = 0;
  };
  std::vector<Window> window_;
  // Per-vertex state of one search (pages reset on first touch, see
  // PagedLabels); marks tell which fields hold data of this search.
  enum VertexMark : uint8_t
  {
    kMarkInterval = 1,  // interval_id: the vertex's interval
    kMarkSource = 2,
    kMarkTarget = 4,    // target_index
    kMarkInfo = 8,      // info
    kMarkPending = 16,  // pending_*: best deferred label
  };
  struct VertexLabel
  {
    // Best deferred label of a vertex outside every interval (see
    // relaxations in process): delta and predecessor.
    Cost pending_delta = 0;
    int32_t interval_id = 0;
    int32_t target_index = 0;
    int32_t pending_pred_interval = 0;
    int32_t pending_pred_k = 0;
    uint16_t info = 0;
    uint8_t marks = 0;
  };
  bool nearest_first_ = false;
  Cost tieKey(const Cost pi) const { return nearest_first_ ? -pi : pi; }
  static constexpr Cost kSinkTie = kInfCost;
  // Event::interval of a deferred label; Event::aux indexes pending_.
  static constexpr int32_t kPendingEvent = -2;
  struct Pending
  {
    int32_t z, i, j;
  };
  std::vector<Pending> pending_;
  // Event::interval of the jogs of a settled run J = [f, g) of a plain
  // interval (relaxed by relaxJogs); Event::aux indexes jog_runs_.
  static constexpr int32_t kJogEvent = -3;
  struct JogRun
  {
    int32_t interval;
    int32_t f, g;
    int32_t di;
    Cost delta;
  };
  std::vector<JogRun> jog_runs_;
  void relaxJogs(const JogRun& run);
  // While direct_, pushEvent stores the event in direct_event_ instead.
  bool direct_ = false;
  Event direct_event_{0, 0, -1, 0};
  using Labels = PagedLabels<VertexLabel, /*kClear=*/true>;
  mutable Labels labels_;
  // labels_[li] through the last page used (most lookups stay on a track).
  VertexLabel& label(const size_t li) const
  {
    if ((li >> kLabelPageBits) == label_page_) {
      return label_base_[li & kLabelPageMask];
    }
    VertexLabel& l = labels_[li];
    if (!labels_.overflow()) {
      label_page_ = li >> kLabelPageBits;
      label_base_ = &l - (li & kLabelPageMask);
    }
    return l;
  }
  static constexpr int kLabelPageBits = Labels::kPageBits;
  static constexpr size_t kLabelPageMask = Labels::kPageSize - 1;
  mutable size_t label_page_ = ~size_t{0};
  mutable VertexLabel* label_base_ = nullptr;
  // Per target box of the track being built: cross cost, along range.
  std::vector<std::array<Cost, 3>> track_future_;
  uint32_t epoch_ = 0;
  // intervals_[0, num_used_) belong to the current search; the others keep
  // their label storage for later searches (at most kKeptIntervals).
  std::vector<Interval> intervals_;
  int num_used_ = 0;
  static constexpr size_t kKeptIntervals = size_t{1} << 16;
  std::vector<Event> heap_;
  int64_t num_intervals_ = 0;
  std::array<int64_t, kNumIntervalEnds> ends_{};
  int64_t interval_vertices_ = 0;
  int64_t num_labels_ = 0;
};

}  // namespace drt::br
