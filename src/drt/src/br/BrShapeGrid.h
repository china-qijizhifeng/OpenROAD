// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Shape grid of BonnRoute (Gester et al., TODAES 18(2), 2013, Sec. 3.3).
// No drt dependencies.
//
// "The shape grid partitions the chip area on each wiring layer and on the
// via layers into axis-parallel rectangular cells, such that each of them
// has at most one neighbor cell to the left, right, bottom and top,
// respectively. The size of the cells is small enough such that shapes of
// different nets cannot be legally present in the same cell."
//
// "In each cell, all intersections of shapes with its area are stored with
// coordinates relative to an anchor point [...], plus their respective type
// and minimum distance requirements. [...] we store this data indirectly by
// a cell configuration number which is used as an index into a lookup table
// that stores the actual data."
//
// "The data volume is further reduced by grouping neighboring cells in
// preferred routing direction to intervals which are stored in an AVL-tree
// in each row or column of cells [...]. Intervals of cells with
// configuration 0 are not stored."
//
// "Moreover, for each nonempty interval we store the net that the shapes of
// this interval belong to if they are removable. The type and net
// identifier of a shape are used to assign a ripup level to it."
//
// Concurrency (Sec. 5.1): each row is split at the global routing tile
// borders into chunks with their own AVL tree, node pool and shared-cell
// table. Threads that route in disjoint tile regions therefore never touch
// the same tree. The configuration table is shared, append-only and locked
// on insertion; writers memoise configuration transitions per thread.
//
// Wide shapes (ShapeClassTable::wide(), e.g. macro blockages and power
// straps) interact over larger distances than all other shapes of their
// layer. They are kept whole in a separate index of per-tile-chunk lists
// instead of the cells, and collect() adds those within the extra distance
// (classReach - narrowReach) of the area it reports, so that windows bloated by
// the narrow reach see them exactly.

#pragma once

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <utility>
#include <vector>

#include "br/BrGeometry.h"
#include "br/BrWireModel.h"

namespace drt::br {

enum class ShapeType : uint8_t
{
  kBlockage = 0,     // fixed, no net, never removable
  kPin = 1,          // fixed, owned by its net
  kReservation = 2,  // owned by a net, non-removable for other nets
  kWire = 3,         // routed wire, removable
  kVia = 4,          // routed via shape, removable
};

inline bool isRemovable(const ShapeType t)
{
  return t == ShapeType::kWire || t == ShapeType::kVia;
}
inline bool isOwned(const ShapeType t)
{
  return t != ShapeType::kBlockage;
}

// Ripup levels (Sec. 3.3, 3.6): 0 = legal, 1..6 = legal after ripping up
// shapes of at most that level, 7 = blocked by non-removable shapes.
using RipupLevel = uint8_t;
constexpr RipupLevel kLegal = 0;
constexpr RipupLevel kMaxRemovableLevel = 6;
constexpr RipupLevel kNotRemovable = 7;
constexpr int32_t kNoNet = -1;

// One shape intersection stored in a cell configuration, relative to the
// cell's anchor (its lower-left corner).
struct CellShape
{
  Coord xlo = 0;
  Coord ylo = 0;
  Coord xhi = 0;
  Coord yhi = 0;
  ClassId cls = 0;
  ShapeType type = ShapeType::kBlockage;
  // Line-end ends of the original shape inside this piece (see PlacedShape).
  uint8_t line_end_ends = 0;
  friend auto operator<=>(const CellShape&, const CellShape&) = default;
};

struct CellConfiguration
{
  std::vector<CellShape> shapes;  // sorted multiset
  bool has_owned = false;         // contains shapes that belong to a net
};

// The cell configuration table; configuration 0 is the empty cell.
class CellConfigurationTable
{
 public:
  CellConfigurationTable();
  CellConfigurationTable(const CellConfigurationTable&) = delete;
  CellConfigurationTable& operator=(const CellConfigurationTable&) = delete;
  // Thread-safe.
  uint32_t intern(std::vector<CellShape> shapes);
  const CellConfiguration& operator[](const uint32_t id) const
  {
    return table_[id];
  }
  size_t size() const { return table_.size(); }
  void clear();

 private:
  struct Hash
  {
    size_t operator()(const std::vector<CellShape>& v) const;
  };
  std::mutex mutex_;
  std::unordered_map<std::vector<CellShape>, uint32_t, Hash> index_;
  AppendOnlyTable<CellConfiguration> table_;
};

// Cell geometry of one shape-grid layer. Rows run along the preferred
// direction; cell (row, cell) covers along [along_origin + cell *
// along_cell, + along_cell) and cross [cross_origin + row * cross_cell, +
// cross_cell).
struct ShapeGridLayerGeometry
{
  bool horizontal = true;
  Coord along_origin = 0;
  Coord cross_origin = 0;
  Coord along_cell = 1;
  Coord cross_cell = 1;
  int along_cells = 1;
  int cross_cells = 1;
};

// Global routing tile axis (the partition regions of Sec. 5.1 are unions
// of tiles). Tile of a coordinate as in BonnRouteDR::tileX / tileY.
struct ChunkAxis
{
  Coord origin = 0;
  Coord size = 1 << 30;
  int count = 1;
  int of(const Coord c) const
  {
    return std::clamp(static_cast<int>((c - origin) / size), 0, count - 1);
  }
};

// A stored shape piece at its absolute position, as the distance rule
// checker sees it.
struct GridPiece
{
  Rect rect;
  ClassId cls = 0;
  ShapeType type = ShapeType::kBlockage;
  uint8_t line_end_ends = 0;
  int32_t owner = kNoNet;
};

// Area [along_lo, along_hi) x [cross_lo, cross_hi) whose content is known
// completely after a collect().
struct KnownRow
{
  Coord cross_lo = 0;
  Coord cross_hi = 0;
  Coord along_lo = 0;
  Coord along_hi = 0;
};

struct GridWindow
{
  int layer = 0;
  bool horizontal = true;
  std::vector<GridPiece> pieces;
  std::vector<KnownRow> rows;  // consecutive rows, increasing cross
  int64_t intervals = 0;       // shape-grid intervals visited
  void clear()
  {
    pieces.clear();
    rows.clear();
    intervals = 0;
  }
};

// Per-thread memo of configuration transitions.
class ShapeGridWriteCache
{
 public:
  void clear() { transitions_.clear(); }

 private:
  friend class ShapeGrid;
  struct Key
  {
    uint32_t config;
    CellShape shape;
    bool insert;
    friend bool operator==(const Key&, const Key&) = default;
  };
  struct Hash
  {
    size_t operator()(const Key& k) const;
  };
  std::unordered_map<Key, uint32_t, Hash> transitions_;
};

class ShapeGrid
{
 public:
  // Owner value of an interval whose single cell holds shapes of several
  // nets (only for inputs that already violate rules, e.g. overlapping
  // pins of different nets); its parts are kept per owner.
  static constexpr int32_t kSharedCell = -2;

  ShapeGrid() = default;
  ShapeGrid(const ShapeGrid&) = delete;
  ShapeGrid& operator=(const ShapeGrid&) = delete;
  ~ShapeGrid();

  void init(const ShapeClassTable* classes,
            std::vector<ShapeGridLayerGeometry> layers,
            ChunkAxis x_chunks,
            ChunkAxis y_chunks);
  int numLayers() const { return static_cast<int>(layers_.size()); }
  const ShapeGridLayerGeometry& geometry(const int layer) const
  {
    return layers_[layer].geom;
  }
  const ShapeClassTable& classes() const { return *classes_; }
  const CellConfigurationTable& configurations() const { return configs_; }

  // Adds / removes a shape (rect at absolute coordinates). A removal must
  // repeat the arguments of an insertion. Shapes are clipped to the grid
  // area. Returns false if a removal did not find the shape.
  bool insert(int layer,
              const Rect& rect,
              ClassId cls,
              ShapeType type,
              int32_t owner,
              ShapeGridWriteCache& cache,
              uint8_t line_end_ends = 3);
  bool remove(int layer,
              const Rect& rect,
              ClassId cls,
              ShapeType type,
              int32_t owner,
              ShapeGridWriteCache& cache,
              uint8_t line_end_ends = 3);

  // Enumerates the stored pieces of all intervals of the rows and cells
  // meeting `window` (closed), plus the area whose content is thereby
  // known: per row from the end of the preceding interval to the start of
  // the following one, within the tile chunks the window meets and at most
  // kKnownHorizon cells beyond the window. Wide shapes are reported whole
  // (clipped to the grid area), once each, if they meet the window or the
  // known area bloated by classReach - narrowReach of the layer; a window that
  // extends at least narrowReach beyond a model therefore holds every
  // shape the model can violate, also for its shifts inside the known area.
  static constexpr int kKnownHorizon = 128;
  void collect(int layer, const Rect& window, GridWindow& out) const;

  // Sec. 3.3: ripup level from shape type and net.
  RipupLevel level(const ShapeType type, const int32_t owner) const
  {
    if (!isRemovable(type)) {
      return kNotRemovable;
    }
    return owner >= 0 && owner < static_cast<int32_t>(net_levels_.size())
               ? net_levels_[owner]
               : 1;
  }
  // Removable shapes of `net` get ripup level `level` (1..6); default 1.
  void setNetRipupLevel(int32_t net, RipupLevel level);

  struct Statistics
  {
    int64_t intervals = 0;
    int64_t interval_cells = 0;
    int64_t shared_cells = 0;
    int64_t chunks = 0;
    size_t configurations = 0;
    size_t configuration_shapes = 0;
    int64_t wide_shapes = 0;
    size_t bytes = 0;  // trees, cell maps, configuration table, wide index
  };
  Statistics statistics() const;
  // Releases spare capacity after a bulk build.
  void compact();

  // Tests: visit every stored interval of a layer.
  struct IntervalView
  {
    int row;
    int begin;
    int end;  // inclusive cells
    int chunk;
    uint32_t config;
    int32_t owner;
    // Shared cells: (owner, config) parts.
    std::vector<std::pair<int32_t, uint32_t>> parts;
  };
  void forEachInterval(
      int layer,
      const std::function<void(const IntervalView&)>& fn) const;
  int chunkOfCell(int layer, int cell) const;
  // Tests: visit every stored wide shape of a layer once.
  void forEachWideShape(int layer,
                        const std::function<void(const GridPiece&)>& fn) const;

 private:
  // A wide shape in the lists of the tile chunks it meets.
  struct WideShape
  {
    GridPiece piece;
    int cell_chunk_lo = 0;
    int cell_chunk_hi = 0;
    int row_chunk_lo = 0;
    int row_chunk_hi = 0;
  };
  // Writers are the threads owning the chunk (Sec. 5.1); readers of
  // collect() windows that reach into another region's chunk only take
  // pieces that do not affect their own answers, under the lock.
  struct WideList
  {
    mutable std::mutex mutex;
    std::atomic<uint32_t> size{0};
    std::vector<WideShape> shapes;
  };
  // One stored interval of a row chunk: cells [begin, end] (absolute cell
  // indices) with equal configuration and owner. `value` is the
  // configuration, or the shared-cell index when owner == kSharedCell.
  struct Interval
  {
    int32_t begin;
    int32_t end;
    uint32_t value;
    int32_t owner;
  };
  struct SharedPart
  {
    int32_t owner;
    uint32_t config;
  };
  // Rows and cells of one layer anchored in one global routing tile. The
  // intervals of all its rows lie in one array in row order (compressed
  // rows): row r holds pool[row_offset[r], row_offset[r + 1]), sorted by
  // begin and disjoint. An update shifts the tail of the block's array.
  struct Block
  {
    int row_first = 0;
    int rows = 0;
    int cell_first = 0;
    int cells = 0;
    std::vector<uint32_t> row_offset;  // rows + 1
    std::vector<Interval> pool;
    std::vector<std::vector<SharedPart>> shared;
    std::vector<uint32_t> free_shared;
  };
  struct Layer
  {
    ShapeGridLayerGeometry geom;
    int along_chunks = 1;
    int cross_chunks = 1;
    std::vector<int> chunk_first_cell;  // along_chunks + 1
    std::vector<int> chunk_first_row;   // cross_chunks + 1
    std::vector<int> cell_chunk;        // per cell
    std::vector<int> row_chunk;         // per row
    std::vector<std::unique_ptr<Block>> blocks;
    // Wide shapes per tile chunk (row chunk * along_chunks + cell chunk).
    std::unique_ptr<WideList[]> wide;
  };
  // State of a run of cells (an interval, or a gap with value 0 and no
  // owner).
  using Segment = Interval;

  bool update(int layer,
              const Rect& rect,
              ClassId cls,
              ShapeType type,
              int32_t owner,
              bool insert,
              uint8_t line_end_ends,
              ShapeGridWriteCache& cache);
  bool updateWide(int layer,
                  const Rect& rect,
                  ClassId cls,
                  ShapeType type,
                  int32_t owner,
                  bool insert,
                  uint8_t line_end_ends);
  void collectWide(const Layer& layer, const Rect& area, GridWindow& out) const;
  bool apply(int layer,
             const Rect& rect,
             ClassId cls,
             ShapeType type,
             int32_t owner,
             bool insert,
             uint8_t line_end_ends,
             ShapeGridWriteCache& cache,
             bool dry_run);
  bool containsRun(const Layer& layer,
                   int row,
                   int cell_begin,
                   int cell_end,
                   const CellShape& piece,
                   int32_t owner) const;
  bool applyRun(Layer& layer,
                int row,
                int cell_begin,
                int cell_end,
                const CellShape& piece,
                int32_t owner,
                bool insert,
                ShapeGridWriteCache& cache);
  Block& block(Layer& layer, int row_chunk, int cell_chunk);
  const Block* block(const Layer& layer, int row_chunk, int cell_chunk) const;
  uint32_t transition(uint32_t config,
                      const CellShape& piece,
                      bool insert,
                      ShapeGridWriteCache& cache);
  // Applies a piece to one segment; may expand it into shared cells.
  bool applyToSegment(Block& b,
                      const Segment& seg,
                      const CellShape& piece,
                      int32_t owner,
                      bool insert,
                      ShapeGridWriteCache& cache,
                      std::vector<Segment>& out);
  bool applyToShared(Block& b,
                     uint32_t index,
                     const CellShape& piece,
                     int32_t owner,
                     bool insert,
                     ShapeGridWriteCache& cache,
                     Segment& cell);
  void emitPieces(const Layer& layer,
                  int row,
                  int begin,
                  int end,
                  uint32_t config,
                  int32_t owner,
                  int clip_lo,
                  int clip_hi,
                  GridWindow& out) const;

  // Sorted interval arrays of the rows (Sec. 3.3 interval sets per row).
  static const Interval* rowBegin(const Block& b, int local)
  {
    return b.pool.data() + b.row_offset[local];
  }
  static const Interval* rowEnd(const Block& b, int local)
  {
    return b.pool.data() + b.row_offset[local + 1];
  }
  // First interval of the row with end >= cell.
  static const Interval* firstEndingAtOrAfter(const Block& b,
                                              int local,
                                              int32_t cell);
  // Last interval with begin <= cell / first with begin >= cell; nullptr
  // if none.
  static const Interval* findLast(const Block& b, int local, int32_t cell);
  static const Interval* findFirst(const Block& b, int local, int32_t cell);
  // Replaces the row's intervals [first, last) by `with`.
  static void replaceRow(Block& b,
                         int local,
                         uint32_t first,
                         uint32_t last,
                         const std::vector<Interval>& with);

  const ShapeClassTable* classes_ = nullptr;
  std::vector<Layer> layers_;
  ChunkAxis x_chunks_;
  ChunkAxis y_chunks_;
  CellConfigurationTable configs_;
  std::vector<RipupLevel> net_levels_;
};

// Largest cell edge for a layer such that two pieces of different nets in
// one cell violate a rule: the cell diagonal must not exceed the smallest
// diff-net spacing (Sec. 3.3). Via layers with an inter-layer rule also
// hold projections of the lower via layer, which must keep that layer's cut
// spacing.
Coord maximalCellSize(const ShapeClassTable& classes, int layer);

}  // namespace drt::br
