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
  virtual RipupLevel vertexLevel(int z,
                                 int i,
                                 int j,
                                 int wire_type,
                                 Kind kind) const = 0;
  // Level of an edge kind (see Kind): includes both endpoints.
  virtual RipupLevel edgeLevel(int z,
                               int i,
                               int j,
                               int wire_type,
                               Kind kind) const = 0;
};

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

  const TrackGraph* graph_ = nullptr;
  const LegalityOracle* oracle_ = nullptr;
  int num_wire_types_ = 1;
  std::vector<Layer> layers_;
  std::unordered_map<size_t, uint16_t> history_;
};

// Keeps a fast grid equal to the distance rule checker's answers ("We
// store continuously updated data computed by the distance rule checking
// module", Sec. 3.6) and answers the remaining queries exactly through the
// checker. Wire types 0 .. fast - 1 of `types` are in the fast grid.
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
  // Recomputes everything a shape on `grid_layer` meeting `rect` (absolute)
  // can change; call after the shape grid changed.
  void shapeChanged(int grid_layer, const Rect& rect);

  RipupLevel vertexLevel(int z,
                         int i,
                         int j,
                         int wire_type,
                         Kind kind) const override;
  RipupLevel edgeLevel(int z,
                       int i,
                       int j,
                       int wire_type,
                       Kind kind) const override;
  // Exact level of a vertex or edge kind with the shapes owned by
  // `ignored_owner` treated as removed (Sec. 4.4); kNoNet: all shapes.
  RipupLevel exactLevel(int z,
                        int i,
                        int j,
                        int wire_type,
                        Kind kind,
                        int32_t ignored_owner) const;

 private:
  struct Range
  {
    int i0, i1;  // tracks
    int j0, j1;  // positions
    int p0, p1;  // jog pairs (p, p + 1)
    int e0, e1;  // edges (j, j + 1)
  };
  // Tracks / positions of layer z whose models can meet `rect`.
  Range wiringRange(int z, const Rect& rect, int grid_layer) const;
  void recomputeWiring(int z, const Range& r, bool vertices, bool jogs) const;
  void recomputeVia(int z, int i0, int i1, int c0, int c1) const;
  // Jog-edge bits of pair (p, p + 1) over positions [j0, j1] (reads only).
  void jogFlags(int z,
                int p,
                int j0,
                int j1,
                std::vector<uint8_t>& flags) const;
  void writeFlags(int z,
                  int p,
                  int j0,
                  const std::vector<uint8_t>& flags,
                  int bit) const;
  void pointModels(int z, std::vector<LocatedModel>& out) const;
  void viaModels(int z, std::vector<LocatedModel>& out) const;
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
