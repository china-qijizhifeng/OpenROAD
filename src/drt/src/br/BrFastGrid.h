// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Fast grid of BonnRoute (Gester et al., TODAES 18(2), 2013, Sec. 3.6).
// No drt dependencies.
//
// "BonnRoute stores continuously updated data computed by the distance rule
// checking module for a restricted set of wire types to be placed at
// on-track locations. [...] the fast grid stores information for vertices,
// but not for edges. This allows to group longer sequences of vertices
// with identical data to intervals. [...] the fast grid also stores
// intervals in via layers that run in the same direction as the next lower
// wiring layer. If legality of a wire on an edge {v, w} cannot be deduced
// from legality at v and w because off-track wires or other shapes are
// present in the vicinity, an extra bit at one of the vertices v and w
// encodes this."
//
// "Per wire type, legality information for the wire models used in
// preferred and non-preferred wiring direction, respectively, has to be
// stored, plus information for via bottom and top pads, that is, 4 types of
// shapes. Hence, 12 bits per wire type and interval suffice to encode
// legality of these shapes in eight different ripup levels. A double word,
// that is, 64 bits, therefore suffices to store information for five wire
// types. In via layers fewer bits are required as there are only two types
// of shapes (via cut shapes and their projections for inter-layer via
// distance checking)."
//
// Wiring-layer word: bits 12 t + 3 k .. + 2 hold the level of shape kind k
// (0 preferred wire, 1 non-preferred wire, 2 via bottom pad, 3 via top pad)
// of fast wire type t < 5; bit 60 marks the preferred-direction edge to
// the next position and bit 61 the jog edge to the next track as edges that
// the distance rule checker must answer. Via-layer word: bits 6 t + 3 k ..
// + 2 for k = 0 cut, 1 projection.
//
// Intervals are kept per track and global routing tile column (chunk), so
// threads routing disjoint tile regions (Sec. 5.1) write disjoint chunks.

#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "br/BrDistanceRules.h"
#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrTrackGraph.h"
#include "br/BrWireModel.h"

namespace drt::br {

// Exact legality from the distance rule checking module, for wire types
// outside the fast grid and for edges whose bit is set.
class LegalityOracle
{
 public:
  virtual ~LegalityOracle() = default;
  // ignored_owner: shapes of this net count as removed (Sec. 4.4).
  virtual RipupLevel vertexLevel(int z,
                                 int i,
                                 int j,
                                 int wire_type,
                                 Kind kind,
                                 int32_t ignored_owner) const = 0;
  // Level of an edge kind (see Kind): includes both endpoints.
  virtual RipupLevel edgeLevel(int z,
                               int i,
                               int j,
                               int wire_type,
                               Kind kind,
                               int32_t ignored_owner) const = 0;
};

#define BR_FAST_GRID_READER 1

// Who answered the path search's legality queries.
struct LegalityCounters
{
  int64_t fast_grid = 0;
  int64_t checker = 0;
};

class FastGrid
{
 public:
  static constexpr int kMaxWireTypes = 5;
  static constexpr int kWiringKinds = 4;
  static constexpr int kViaKinds = 2;
  static constexpr int kWireEdgeBit = 60;
  static constexpr int kJogEdgeBit = 61;
  static constexpr uint64_t kEdgeBits
      = (uint64_t{1} << kWireEdgeBit) | (uint64_t{1} << kJogEdgeBit);

  static int wiringShift(const int wire_type, const int kind)
  {
    return 12 * wire_type + 3 * kind;
  }
  static int viaShift(const int wire_type, const int kind)
  {
    return 6 * wire_type + 3 * kind;
  }
  static RipupLevel field(const uint64_t word, const int shift)
  {
    return static_cast<RipupLevel>((word >> shift) & 7);
  }

  // All vertices legal, no edge bits.
  void init(const TrackGraph* graph, int num_wire_types);
  void setOracle(const LegalityOracle* oracle) { oracle_ = oracle; }
  const LegalityOracle* oracle() const { return oracle_; }
  int numWireTypes() const { return num_wire_types_; }
  const TrackGraph& graph() const { return *graph_; }

  uint64_t wiringWord(int z, int i, int j) const;
  uint64_t viaWord(int z, int i, int c) const;
  // The word at (z, i, j) / via position c and the positions [lo, hi] of
  // its run (within one tile).
  uint64_t wiringRun(int z, int i, int j, int& lo, int& hi) const;
  // The runs of track i in the tile of position j: positions first ..
  // first + count - 1, n runs (n = 0: all zero).
  struct RunsView
  {
    const uint16_t* starts = nullptr;
    const uint64_t* words = nullptr;
    int n = 0;
    int first = 0;
    int count = 0;
  };
  RunsView wiringRuns(int z, int i, int j) const;
  uint64_t viaRun(int z, int i, int c, int& lo, int& hi) const;
  // Level of a vertex or edge kind for wire type `wire_type`; wire types
  // outside the fast grid and marked edges go to the oracle.
  RipupLevel level(int z,
                   int i,
                   int j,
                   Kind kind,
                   int wire_type,
                   LegalityCounters* counters = nullptr) const;
  RipupLevel level(size_t v,
                   Kind kind,
                   int wire_type = 0,
                   LegalityCounters* counters = nullptr) const
  {
    int z, i, j;
    graph_->decode(v, z, i, j);
    return level(z, i, j, kind, wire_type, counters);
  }
  bool edgeQuery(int z, int i, int j, Kind edge_kind) const;

  // Positions [j0, j1] of track i get (word & mask) in the bits of mask.
  void assignWiring(int z, int i, int j0, int j1, uint64_t word, uint64_t mask);
  void assignVia(int z, int i, int c0, int c1, uint64_t word, uint64_t mask);
  // Positions j0 .. j0 + n - 1 of track i get (values[k] & mask).
  void writeWiring(int z,
                   int i,
                   int j0,
                   const uint64_t* values,
                   int n,
                   uint64_t mask);
  void writeVia(int z,
                int i,
                int c0,
                const uint64_t* values,
                int n,
                uint64_t mask);
  void setLevel(size_t v, Kind kind, int wire_type, RipupLevel level);
  void setEdgeQuery(size_t v, Kind edge_kind, bool query);

  // Via positions of layer z: positions with a track of layer z + 1.
  int viaIndex(const int z, const int j) const
  {
    return layers_[z].via_index[j];
  }
  int numViaPositions(const int z) const
  {
    return static_cast<int>(layers_[z].via_position.size());
  }
  int viaPosition(const int z, const int c) const
  {
    return layers_[z].via_position[c];
  }

  // Rip-up history (Sec. 4.2 cost increase); not part of the paper's fast
  // grid, kept sparse beside it. Bumped serially only.
  uint16_t history(const size_t v) const
  {
    const auto it = history_.find(v);
    return it == history_.end() ? 0 : it->second;
  }
  void bumpHistory(size_t v);

  struct Statistics
  {
    int64_t wiring_intervals = 0;
    int64_t via_intervals = 0;
    int64_t vertices = 0;
    int64_t via_vertices = 0;
    int64_t marked_edges = 0;
    int64_t tiles = 0;  // allocated tiles
    size_t bytes = 0;   // run storage, tile tables, per-position maps
  };
  Statistics statistics() const;
  // Releases spare capacity after a bulk build.
  void compact();

 private:
  // Runs of the tracks of one global routing tile: track t's runs are
  // [track_first[t], track_first[t + 1]); a run holds the word from its
  // start (relative position) to the next start. A track without runs is
  // all zero (all legal, no marks). Tiles are allocated on the first
  // non-zero word; a tile is written by one thread at a time (Sec. 5.1).
  struct TileRuns
  {
    std::vector<uint32_t> track_first;
    std::vector<uint16_t> starts;
    std::vector<uint64_t> words;
  };
  struct Layer
  {
    int pos_chunks = 1;
    int track_chunks = 1;
    std::vector<int> pos_chunk;                     // per position
    std::vector<int> chunk_first_pos;               // pos_chunks + 1
    std::vector<int> track_chunk;                   // per track
    std::vector<int> chunk_first_track;             // track_chunks + 1
    std::vector<std::unique_ptr<TileRuns>> wiring;  // [track chunk][chunk]
    std::vector<int> via_index;                     // per position, -1
    std::vector<int> via_position;                  // per via position
    std::vector<int> via_chunk;                     // per via position
    std::vector<int> via_chunk_first;               // pos_chunks + 1
    std::vector<std::unique_ptr<TileRuns>> via;
  };
  static uint64_t lookup(const TileRuns* tile, int track, int index);
  static uint64_t lookupRun(const TileRuns* tile,
                            int track,
                            int index,
                            int count,
                            int& lo,
                            int& hi);
  // Positions [lo, hi] (relative) of `track` in a tile of `tracks` tracks
  // and `count` positions get (word & mask) in the bits of mask.
  static void assign(std::unique_ptr<TileRuns>& tile,
                     int tracks,
                     int track,
                     int count,
                     int lo,
                     int hi,
                     uint64_t word,
                     uint64_t mask);
  static size_t bytes(const std::vector<std::unique_ptr<TileRuns>>& tiles);
  static void assignWords(std::unique_ptr<TileRuns>& tile,
                          int tracks,
                          int track,
                          int count,
                          int lo,
                          const uint64_t* values,
                          int n,
                          uint64_t mask);

  const TrackGraph* graph_ = nullptr;
  const LegalityOracle* oracle_ = nullptr;
  int num_wire_types_ = 1;
  std::vector<Layer> layers_;
  std::unordered_map<size_t, uint16_t> history_;
};

// Exact fast-grid words of the tracks near a net's own shapes with those
// shapes removed (Sec. 4.4: "The shapes of the components that contributed
// to S and T [...] are removed temporarily from routing space"). The
// searches of the net read the fast grid through it instead of the shapes
// being removed and reinserted physically.
struct FastGridOverlay
{
  struct Run
  {
    int32_t lo;
    int32_t hi;
    uint64_t word;
  };
  int32_t owner = kNoNet;
  // Key (z << 32) | track; runs sorted by position (via: via position).
  std::unordered_map<uint64_t, std::vector<Run>> wiring;
  std::unordered_map<uint64_t, std::vector<Run>> via;
  void clear()
  {
    owner = kNoNet;
    wiring.clear();
    via.clear();
  }
  static uint64_t key(const int z, const int i)
  {
    return (static_cast<uint64_t>(z) << 32) | static_cast<uint32_t>(i);
  }
};

// The search-side view of the fast grid: run cursors (the last run seen
// per track, direct mapped) make sequential access along a track and its
// neighbours O(1), and an optional overlay replaces the words near the
// searching net's own shapes. One reader per search object / thread.
class FastGridReader
{
 public:
  void reset(const FastGrid* grid, const FastGridOverlay* overlay);
  RipupLevel level(int z,
                   int i,
                   int j,
                   Kind kind,
                   int wire_type,
                   LegalityCounters* counters = nullptr);
  RipupLevel level(size_t v,
                   Kind kind,
                   int wire_type,
                   LegalityCounters* counters = nullptr)
  {
    int z, i, j;
    grid_->graph().decode(v, z, i, j);
    return level(z, i, j, kind, wire_type, counters);
  }
  uint64_t wiringWord(const int z, const int i, const int j)
  {
    Cursor& c = wiring_[slot(z, i)];
    if (c.z == z && c.i == i && j >= c.lo && j <= c.hi) {
      return c.word;
    }
    return fillWiring(c, z, i, j);
  }
  uint64_t viaWord(const int z, const int i, const int v)
  {
    Cursor& c = via_[slot(z, i)];
    if (c.z == z && c.i == i && v >= c.lo && v <= c.hi) {
      return c.word;
    }
    return fillVia(c, z, i, v);
  }
  const FastGrid& grid() const { return *grid_; }
  int32_t ignoredOwner() const
  {
    return overlay_ != nullptr ? overlay_->owner : kNoNet;
  }

 private:
  struct Cursor
  {
    int32_t z = -1;
    int32_t i = 0;
    int32_t lo = 0;
    int32_t hi = -1;
    uint64_t word = 0;
    // The runs of the cursor's tile (wiring, without overlay).
    FastGrid::RunsView runs;
    int32_t run = -1;
  };
  static constexpr int kSlots = 256;
  static int slot(const int z, const int i)
  {
    return ((i << 2) ^ z) & (kSlots - 1);
  }
  uint64_t fillWiring(Cursor& c, int z, int i, int j);
  uint64_t fillVia(Cursor& c, int z, int i, int v);
  static uint64_t applyOverlay(
      const std::unordered_map<uint64_t, std::vector<FastGridOverlay::Run>>&
          runs,
      int z,
      int i,
      int j,
      uint64_t word,
      int& lo,
      int& hi);

  const FastGrid* grid_ = nullptr;
  const FastGridOverlay* overlay_ = nullptr;
  std::array<Cursor, kSlots> wiring_;
  std::array<Cursor, kSlots> via_;
};

namespace detail {

inline int wiringKind(const Kind kind)
{
  switch (kind) {
    case kWire:
      return 0;
    case kJog:
      return 1;
    case kViaUpPad:
      return 2;
    case kViaDownPad:
      return 3;
    default:
      return -1;
  }
}

// Level of a vertex or edge kind from word sources (FastGrid, reader).
template <typename Wiring, typename Via>
inline RipupLevel levelOf(const FastGrid& grid,
                          const int z,
                          const int i,
                          const int j,
                          const Kind kind,
                          const int wire_type,
                          const int32_t ignored,
                          LegalityCounters* counters,
                          Wiring&& wiring,
                          Via&& via)
{
  const bool edge = kind == kWireEdge || kind == kJogEdge;
  const LegalityOracle* oracle = grid.oracle();
  if (wire_type >= grid.numWireTypes()) {
    if (counters != nullptr) {
      counters->checker++;
    }
    if (oracle == nullptr) {
      return kNotRemovable;
    }
    return edge ? oracle->edgeLevel(z, i, j, wire_type, kind, ignored)
                : oracle->vertexLevel(z, i, j, wire_type, kind, ignored);
  }
  if (edge) {
    const uint64_t w = wiring(z, i, j);
    const bool along = kind == kWireEdge;
    if ((w >> (along ? FastGrid::kWireEdgeBit : FastGrid::kJogEdgeBit)) & 1) {
      if (counters != nullptr) {
        counters->checker++;
      }
      return oracle != nullptr
                 ? oracle->edgeLevel(z, i, j, wire_type, kind, ignored)
                 : kNotRemovable;
    }
    if (counters != nullptr) {
      counters->fast_grid++;
    }
    // Deduced from the vertex data at both endpoints.
    const int shift = FastGrid::wiringShift(wire_type, along ? 0 : 1);
    const uint64_t other = along ? wiring(z, i, j + 1) : wiring(z, i + 1, j);
    return std::max(FastGrid::field(w, shift), FastGrid::field(other, shift));
  }
  if (counters != nullptr) {
    counters->fast_grid++;
  }
  if (kind == kCutUp) {
    const int c = grid.viaIndex(z, j);
    if (c < 0) {
      return kNotRemovable;
    }
    const uint64_t w = via(z, i, c);
    return std::max(FastGrid::field(w, FastGrid::viaShift(wire_type, 0)),
                    FastGrid::field(w, FastGrid::viaShift(wire_type, 1)));
  }
  return FastGrid::field(wiring(z, i, j),
                         FastGrid::wiringShift(wire_type, wiringKind(kind)));
}

}  // namespace detail

// Keeps a fast grid equal to the distance rule checker's answers ("We
// store continuously updated data computed by the distance rule checking
// module", Sec. 3.6) and answers the remaining queries exactly through the
// checker. Wire types 0 .. fast - 1 of `types` are in the fast grid.
//
// Updates are incremental and exact: an inserted shape raises the levels
// it violates to max(old, its level), computed from that shape alone; a
// removed shape can only lower levels equal to its own, and only those
// (with the affected edge bits) are recomputed exactly from the shape
// grid, one window per track.
class FastGridBuilder final : public LegalityOracle
{
 public:
  void init(const TrackGraph* graph,
            const std::vector<WireType>* types,
            int fast,
            const DistanceRuleChecker* checker,
            FastGrid* grid);

  // Initial build in two phases (each parallel over track ranges):
  // vertex data, preferred-edge bits and via-layer data of tracks [i0, i1],
  // then the jog-edge bits of the pairs (i, i + 1), which read both tracks.
  void rebuildVertices(int z, int i0, int i1);
  // Jog bits of the pairs (i, i + 1), i in [i0, i1]: computed reading both
  // tracks, then written; across tile rows, compute all before writing any.
  struct JogBits
  {
    int z;
    int p;
    std::vector<uint8_t> flags;  // per position
  };
  void computeJogBits(int z, int i0, int i1, std::vector<JogBits>& out) const;
  void writeJogBits(const std::vector<JogBits>& bits) const;
  // Serial convenience: compute and write.
  void rebuildJogBits(int z, int i0, int i1);
  // Track ranges (z, i0, i1) of whole tile rows: jobs that may run on
  // parallel threads within one phase (a tile is written by one thread).
  std::vector<std::array<int, 3>> rebuildJobs() const;

  // Incremental updates; call after the shape grid changed. `shape` is the
  // inserted / removed shape (absolute) and `level` its ripup level.
  void shapeInserted(int grid_layer,
                     const PlacedShape& shape,
                     RipupLevel level);
  void shapeRemoved(int grid_layer, const PlacedShape& shape, RipupLevel level);
  // Removals of a group of shapes (a net's rip-up, its reservations): the
  // candidates are collected and recomputed once per layer by flush().
  struct PendingRemovals;
  void shapeRemoved(int grid_layer,
                    const PlacedShape& shape,
                    RipupLevel level,
                    PendingRemovals& pending);
  void flush(PendingRemovals& pending);
  // Exact recomputation of everything a shape meeting `rect` can change.
  void shapeChanged(int grid_layer, const Rect& rect);

  // The overlay of `owner` for its own shapes (grid layer, shape, level).
  struct OwnShape
  {
    int grid_layer;
    PlacedShape shape;
    RipupLevel level;
  };
  void buildOverlay(int32_t owner,
                    const std::vector<OwnShape>& shapes,
                    FastGridOverlay& overlay) const;

  RipupLevel vertexLevel(int z,
                         int i,
                         int j,
                         int wire_type,
                         Kind kind,
                         int32_t ignored_owner) const override;
  RipupLevel edgeLevel(int z,
                       int i,
                       int j,
                       int wire_type,
                       Kind kind,
                       int32_t ignored_owner) const override;
  // Exact level of a vertex or edge kind with the shapes owned by
  // `ignored_owner` treated as removed (Sec. 4.4); kNoNet: all shapes.
  RipupLevel exactLevel(int z,
                        int i,
                        int j,
                        int wire_type,
                        Kind kind,
                        int32_t ignored_owner) const;

  struct UpdateStatistics
  {
    int64_t inserts = 0;
    int64_t removes = 0;
    int64_t recomputed_positions = 0;  // exact per-track recomputation
    int64_t exact_bits = 0;            // edge bits checked exactly
  };

 private:
  struct Range
  {
    int i0, i1;  // tracks
    int j0, j1;  // positions
    int p0, p1;  // jog pairs (p, p + 1)
    int e0, e1;  // edges (j, j + 1)
  };
  // Exact point levels of track i over positions [j0, j1] (level fields
  // only) and, per fast wire type, the largest preferred-wire point level
  // on each edge (j, j + 1), j in [j0, j1 - 1].
  struct TrackLevels
  {
    std::vector<uint64_t> words;
    std::vector<RipupLevel> interior;  // [t * edges + (j - j0)]
    GridWindow pieces;
  };
  // Tracks / positions of layer z whose models can meet `rect`.
  Range wiringRange(int z, const Rect& rect, int grid_layer) const;
  void trackLevels(int z,
                   int i,
                   int j0,
                   int j1,
                   int32_t ignored,
                   TrackLevels& out) const;
  // Jog hull levels of pair (p, p + 1) over [j0, j1]: [t * n + (j - j0)].
  void jogLevels(int z,
                 int p,
                 int j0,
                 int j1,
                 int32_t ignored,
                 std::vector<RipupLevel>& out) const;
  // Exact full words (levels and both edge bits) of tracks [i0, i1] over
  // positions [j0, j1]: [(i - i0) * n + (j - j0)].
  void computeWiring(int z,
                     int i0,
                     int i1,
                     int j0,
                     int j1,
                     int32_t ignored,
                     std::vector<uint64_t>& words) const;
  // Exact via words of tracks [i0, i1] over via positions [c0, c1].
  void computeVia(int z,
                  int i0,
                  int i1,
                  int c0,
                  int c1,
                  int32_t ignored,
                  std::vector<uint64_t>& words) const;
  void writeWiring(int z,
                   int i0,
                   int i1,
                   int j0,
                   int j1,
                   const std::vector<uint64_t>& words,
                   uint64_t mask) const;
  void writeVia(int z,
                int i0,
                int i1,
                int c0,
                int c1,
                const std::vector<uint64_t>& words) const;
  void viaRange(int z,
                int grid_layer,
                const Rect& rect,
                const Rect& box,
                int& i0,
                int& i1,
                int& c0,
                int& c1) const;
  // Candidates of a removal (or of an overlay): (track, position) whose
  // word can change: a level equal to the removed shape's, or a marked
  // edge whose swept wire the shape violates.
  struct Candidates
  {
    std::vector<std::pair<int, int>> at;
    int i0 = 1 << 30, i1 = -1, j0 = 1 << 30, j1 = -1;
    void add(int i, int j)
    {
      at.emplace_back(i, j);
      i0 = std::min(i0, i);
      i1 = std::max(i1, i);
      j0 = std::min(j0, j);
      j1 = std::max(j1, j);
    }
    bool empty() const { return at.empty(); }
  };

 public:
  struct PendingRemovals
  {
    std::vector<Candidates> wiring;  // per wiring layer
    bool empty() const
    {
      for (const Candidates& c : wiring) {
        if (!c.empty()) {
          return false;
        }
      }
      return true;
    }
  };

 private:
  // Exact words at the candidates of wiring layer z and around them (edge
  // bits at both sides, jog bits of the pairs below and above), the shapes
  // of `ignored` treated as removed: (track, position, word) appended to
  // out. Words elsewhere are unchanged.
  void recomputeCandidates(
      int z,
      Candidates& cand,
      int32_t ignored,
      std::vector<std::tuple<int, int, uint64_t>>& out) const;
  void wiringCandidates(int z,
                        int grid_layer,
                        const PlacedShape& shape,
                        RipupLevel level,
                        Candidates& out) const;
  void viaCandidates(int z,
                     int grid_layer,
                     const PlacedShape& shape,
                     RipupLevel level,
                     bool projection,
                     Candidates& out) const;
  void insertWiring(int z, int grid_layer, const PlacedShape& s, RipupLevel l);
  void writeChanged(int z,
                    std::vector<std::tuple<int, int, uint64_t>>& changed) const;
  void insertVia(int z,
                 int grid_layer,
                 const PlacedShape& s,
                 RipupLevel l,
                 bool projection);
  void pointModels(int z, std::vector<LocatedModel>& out) const;
  void viaModels(int z, std::vector<LocatedModel>& out) const;
  LocatedModel jogHull(int z, int p, int t) const;
  Rect prefHull(int z, int i, int j, int t) const;
  Point location(int z, int i, int j) const;
  RipupLevel sweptLevel(int z,
                        int i,
                        int j,
                        int wire_type,
                        Kind kind,
                        int32_t ignored_owner) const;

  const TrackGraph* graph_ = nullptr;
  const std::vector<WireType>* types_ = nullptr;
  int fast_ = 0;
  const DistanceRuleChecker* checker_ = nullptr;
  FastGrid* grid_ = nullptr;
  // Bounding boxes of the models with respect to the location, per layer.
  std::vector<Rect> wiring_model_box_;
  std::vector<Rect> cut_model_box_;
  std::vector<Rect> projection_model_box_;
  std::vector<bool> has_projection_;
};

}  // namespace drt::br
