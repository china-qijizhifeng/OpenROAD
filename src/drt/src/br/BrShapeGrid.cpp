// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrShapeGrid.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

#include "br/BrGeometry.h"
#include "br/BrWireModel.h"

namespace drt::br {

namespace {

inline size_t mix(size_t h, const uint64_t v)
{
  h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  return h;
}

size_t hashShape(size_t h, const CellShape& s)
{
  h = mix(h,
          (static_cast<uint64_t>(static_cast<uint32_t>(s.xlo)) << 32)
              | static_cast<uint32_t>(s.ylo));
  h = mix(h,
          (static_cast<uint64_t>(static_cast<uint32_t>(s.xhi)) << 32)
              | static_cast<uint32_t>(s.yhi));
  h = mix(h,
          (static_cast<uint64_t>(s.cls) << 16)
              | (static_cast<uint64_t>(s.type) << 8) | s.line_end_ends);
  return h;
}

}  // namespace

CellConfigurationTable::CellConfigurationTable()
{
  clear();
}

void CellConfigurationTable::clear()
{
  std::lock_guard lock(mutex_);
  table_.clear();
  index_.clear();
  table_.append(CellConfiguration{});
  index_.emplace(std::vector<CellShape>{}, 0);
}

size_t CellConfigurationTable::Hash::operator()(
    const std::vector<CellShape>& v) const
{
  size_t h = v.size();
  for (const CellShape& s : v) {
    h = hashShape(h, s);
  }
  return h;
}

uint32_t CellConfigurationTable::intern(std::vector<CellShape> shapes)
{
  std::ranges::sort(shapes);
  std::lock_guard lock(mutex_);
  const auto it = index_.find(shapes);
  if (it != index_.end()) {
    return it->second;
  }
  CellConfiguration c;
  c.shapes = shapes;
  c.has_owned = std::ranges::any_of(
      shapes, [](const CellShape& s) { return isOwned(s.type); });
  const uint32_t id = static_cast<uint32_t>(table_.append(std::move(c)));
  index_.emplace(std::move(shapes), id);
  return id;
}

size_t ShapeGridWriteCache::Hash::operator()(const Key& k) const
{
  return hashShape(mix(k.config, k.insert ? 1 : 2), k.shape);
}

ShapeGrid::~ShapeGrid() = default;

void ShapeGrid::init(const ShapeClassTable* classes,
                     std::vector<ShapeGridLayerGeometry> layers,
                     const ChunkAxis x_chunks,
                     const ChunkAxis y_chunks)
{
  classes_ = classes;
  x_chunks_ = x_chunks;
  y_chunks_ = y_chunks;
  configs_.clear();
  layers_.clear();
  layers_.resize(layers.size());
  for (size_t l = 0; l < layers.size(); l++) {
    Layer& layer = layers_[l];
    layer.geom = layers[l];
    const ShapeGridLayerGeometry& g = layer.geom;
    const ChunkAxis& along = g.horizontal ? x_chunks : y_chunks;
    const ChunkAxis& cross = g.horizontal ? y_chunks : x_chunks;
    layer.along_chunks = along.count;
    layer.cross_chunks = cross.count;
    layer.cell_chunk.resize(g.along_cells);
    for (int c = 0; c < g.along_cells; c++) {
      layer.cell_chunk[c] = along.of(static_cast<Coord>(
          g.along_origin + static_cast<int64_t>(c) * g.along_cell));
    }
    layer.row_chunk.resize(g.cross_cells);
    for (int r = 0; r < g.cross_cells; r++) {
      layer.row_chunk[r] = cross.of(static_cast<Coord>(
          g.cross_origin + static_cast<int64_t>(r) * g.cross_cell));
    }
    auto firsts = [](const std::vector<int>& of, int count) {
      std::vector<int> first(count + 1, static_cast<int>(of.size()));
      for (int k = count - 1, c = static_cast<int>(of.size()) - 1; k >= 0;
           k--) {
        while (c >= 0 && of[c] >= k) {
          c--;
        }
        first[k] = c + 1;
      }
      return first;
    };
    layer.chunk_first_cell = firsts(layer.cell_chunk, layer.along_chunks);
    layer.chunk_first_row = firsts(layer.row_chunk, layer.cross_chunks);
    layer.blocks.clear();
    layer.blocks.resize(static_cast<size_t>(layer.along_chunks)
                        * layer.cross_chunks);
    layer.wide = std::make_unique<WideList[]>(
        static_cast<size_t>(layer.along_chunks) * layer.cross_chunks);
  }
}

void ShapeGrid::setNetRipupLevel(const int32_t net, const RipupLevel level)
{
  if (net < 0) {
    return;
  }
  if (net >= static_cast<int32_t>(net_levels_.size())) {
    net_levels_.resize(net + 1, 1);
  }
  net_levels_[net] = std::clamp<RipupLevel>(level, 1, kMaxRemovableLevel);
}

int ShapeGrid::chunkOfCell(const int layer, const int cell) const
{
  return layers_[layer].cell_chunk[cell];
}

ShapeGrid::Block& ShapeGrid::block(Layer& layer,
                                   const int row_chunk,
                                   const int cell_chunk)
{
  std::unique_ptr<Block>& b
      = layer.blocks[static_cast<size_t>(row_chunk) * layer.along_chunks
                     + cell_chunk];
  if (!b) {
    b = std::make_unique<Block>();
    b->row_first = layer.chunk_first_row[row_chunk];
    b->rows = layer.chunk_first_row[row_chunk + 1] - b->row_first;
    b->cell_first = layer.chunk_first_cell[cell_chunk];
    b->cells = layer.chunk_first_cell[cell_chunk + 1] - b->cell_first;
    b->row_offset.assign(b->rows + 1, 0);
  }
  return *b;
}

const ShapeGrid::Block* ShapeGrid::block(const Layer& layer,
                                         const int row_chunk,
                                         const int cell_chunk) const
{
  return layer
      .blocks[static_cast<size_t>(row_chunk) * layer.along_chunks + cell_chunk]
      .get();
}

// Interval sets of the rows (Sec. 3.3: "intervals which are stored in an
// AVL-tree in each row"): per row chunk a sorted array of disjoint
// intervals, the rows of a block concatenated (compressed rows). Row
// chunks hold few intervals, so binary search and shifting the block's
// tail replace the tree's rotations; the enumeration order (increasing
// first cell) is the same.

const ShapeGrid::Interval* ShapeGrid::firstEndingAtOrAfter(const Block& b,
                                                           const int local,
                                                           const int32_t cell)
{
  return std::partition_point(
      rowBegin(b, local), rowEnd(b, local), [cell](const Interval& n) {
        return n.end < cell;
      });
}

const ShapeGrid::Interval* ShapeGrid::findLast(const Block& b,
                                               const int local,
                                               const int32_t cell)
{
  const Interval* first = rowBegin(b, local);
  const Interval* it = std::partition_point(
      first, rowEnd(b, local), [cell](const Interval& n) {
        return n.begin <= cell;
      });
  return it == first ? nullptr : it - 1;
}

const ShapeGrid::Interval* ShapeGrid::findFirst(const Block& b,
                                                const int local,
                                                const int32_t cell)
{
  const Interval* last = rowEnd(b, local);
  const Interval* it = std::partition_point(
      rowBegin(b, local), last, [cell](const Interval& n) {
        return n.begin < cell;
      });
  return it == last ? nullptr : it;
}

void ShapeGrid::replaceRow(Block& b,
                           const int local,
                           const uint32_t first,
                           const uint32_t last,
                           const std::vector<Interval>& with)
{
  const uint32_t base = b.row_offset[local];
  const auto n = static_cast<int64_t>(with.size());
  const auto m = static_cast<int64_t>(last - first);
  const int64_t delta = n - m;
  if (delta > 0) {
    // Grow by an eighth: blocks keep little spare room.
    const size_t need = b.pool.size() + delta;
    if (need > b.pool.capacity()) {
      b.pool.reserve(need + need / 8 + 4);
    }
    b.pool.insert(b.pool.begin() + base + last, delta, Interval{});
  } else if (delta < 0) {
    b.pool.erase(b.pool.begin() + base + first + n,
                 b.pool.begin() + base + last);
  }
  std::ranges::copy(with, b.pool.begin() + base + first);
  if (delta != 0) {
    for (int r = local + 1; r <= b.rows; r++) {
      b.row_offset[r] = static_cast<uint32_t>(b.row_offset[r] + delta);
    }
  }
}

uint32_t ShapeGrid::transition(const uint32_t config,
                               const CellShape& piece,
                               const bool insert,
                               ShapeGridWriteCache& cache)
{
  const ShapeGridWriteCache::Key key{config, piece, insert};
  const auto it = cache.transitions_.find(key);
  if (it != cache.transitions_.end()) {
    return it->second;
  }
  std::vector<CellShape> shapes = configs_[config].shapes;
  if (insert) {
    shapes.insert(std::ranges::upper_bound(shapes, piece), piece);
  } else {
    const auto at = std::ranges::lower_bound(shapes, piece);
    shapes.erase(at);
  }
  const uint32_t id = configs_.intern(std::move(shapes));
  cache.transitions_.emplace(key, id);
  return id;
}

namespace {

bool containsShape(const CellConfiguration& c, const CellShape& piece)
{
  return std::ranges::binary_search(c.shapes, piece);
}

}  // namespace

bool ShapeGrid::applyToShared(Block& b,
                              const uint32_t index,
                              const CellShape& piece,
                              const int32_t owner,
                              const bool insert,
                              ShapeGridWriteCache& cache,
                              Segment& cell)
{
  std::vector<SharedPart> parts = b.shared[index];
  const int32_t target = isOwned(piece.type) ? owner : kNoNet;
  auto at = std::ranges::find(parts, target, &SharedPart::owner);
  if (insert) {
    if (at == parts.end()) {
      parts.push_back({target, configs_.intern({piece})});
    } else {
      at->config = transition(at->config, piece, true, cache);
    }
  } else {
    if (at == parts.end() || !containsShape(configs_[at->config], piece)) {
      return false;
    }
    at->config = transition(at->config, piece, false, cache);
    if (at->config == 0) {
      parts.erase(at);
    }
  }
  std::ranges::sort(parts, {}, &SharedPart::owner);
  const auto owned = std::ranges::count_if(
      parts, [](const SharedPart& p) { return p.owner != kNoNet; });
  if (owned > 1) {
    b.shared[index] = std::move(parts);
    return true;
  }
  // Back to a single owner: one ordinary cell configuration.
  std::vector<CellShape> shapes;
  int32_t single = kNoNet;
  for (const SharedPart& p : parts) {
    const auto& s = configs_[p.config].shapes;
    shapes.insert(shapes.end(), s.begin(), s.end());
    if (p.owner != kNoNet) {
      single = p.owner;
    }
  }
  b.shared[index].clear();
  b.free_shared.push_back(index);
  cell.value = configs_.intern(std::move(shapes));
  cell.owner = cell.value == 0 ? kNoNet : single;
  return true;
}

bool ShapeGrid::applyToSegment(Block& b,
                               const Segment& seg,
                               const CellShape& piece,
                               const int32_t owner,
                               const bool insert,
                               ShapeGridWriteCache& cache,
                               std::vector<Segment>& out)
{
  if (seg.owner == kSharedCell) {
    Segment cell = seg;
    const bool ok
        = applyToShared(b, seg.value, piece, owner, insert, cache, cell);
    out.push_back(cell);
    return ok;
  }
  const bool owned = isOwned(piece.type);
  if (insert) {
    if (owned && seg.owner != kNoNet && seg.owner != owner) {
      // Shapes of two nets in one cell: keep their parts per owner.
      const CellConfiguration& c = configs_[seg.value];
      std::vector<CellShape> fixed;
      std::vector<CellShape> owned_shapes;
      for (const CellShape& s : c.shapes) {
        (isOwned(s.type) ? owned_shapes : fixed).push_back(s);
      }
      std::vector<SharedPart> parts;
      if (!fixed.empty()) {
        parts.push_back({kNoNet, configs_.intern(std::move(fixed))});
      }
      parts.push_back({seg.owner, configs_.intern(std::move(owned_shapes))});
      parts.push_back({owner, configs_.intern({piece})});
      std::ranges::sort(parts, {}, &SharedPart::owner);
      for (int cell = seg.begin; cell <= seg.end; cell++) {
        uint32_t index;
        if (!b.free_shared.empty()) {
          index = b.free_shared.back();
          b.free_shared.pop_back();
        } else {
          index = static_cast<uint32_t>(b.shared.size());
          b.shared.emplace_back();
        }
        b.shared[index] = parts;
        out.push_back({cell, cell, index, kSharedCell});
      }
      return true;
    }
    const uint32_t v = transition(seg.value, piece, true, cache);
    out.push_back({seg.begin, seg.end, v, owned ? owner : seg.owner});
    return true;
  }
  if ((owned && seg.owner != owner)
      || !containsShape(configs_[seg.value], piece)) {
    out.push_back(seg);
    return false;
  }
  const uint32_t v = transition(seg.value, piece, false, cache);
  out.push_back(
      {seg.begin, seg.end, v, configs_[v].has_owned ? seg.owner : kNoNet});
  return true;
}

bool ShapeGrid::containsRun(const Layer& layer,
                            const int row,
                            const int cell_begin,
                            const int cell_end,
                            const CellShape& piece,
                            const int32_t owner) const
{
  const Block* b
      = block(layer, layer.row_chunk[row], layer.cell_chunk[cell_begin]);
  if (b == nullptr) {
    return false;
  }
  const int local = row - b->row_first;
  int next = cell_begin;
  const int32_t target = isOwned(piece.type) ? owner : kNoNet;
  for (const Interval *n = firstEndingAtOrAfter(*b, local, cell_begin),
                      *e = rowEnd(*b, local);
       n != e && n->begin <= cell_end;
       ++n) {
    const Interval& node = *n;
    if (node.begin > next) {
      return false;  // a gap holds nothing
    }
    next = node.end + 1;
    if (node.owner == kSharedCell) {
      const auto& parts = b->shared[node.value];
      const auto at = std::ranges::find(parts, target, &SharedPart::owner);
      if (at == parts.end() || !containsShape(configs_[at->config], piece)) {
        return false;
      }
    } else if ((isOwned(piece.type) && node.owner != owner)
               || !containsShape(configs_[node.value], piece)) {
      return false;
    }
  }
  return next > cell_end;
}

bool ShapeGrid::applyRun(Layer& layer,
                         const int row,
                         const int cell_begin,
                         const int cell_end,
                         const CellShape& piece,
                         const int32_t owner,
                         const bool insert,
                         ShapeGridWriteCache& cache)
{
  const int row_chunk = layer.row_chunk[row];
  const int cell_chunk = layer.cell_chunk[cell_begin];
  if (!insert
      && block(std::as_const(layer), row_chunk, cell_chunk) == nullptr) {
    return false;
  }
  Block& b = block(layer, row_chunk, cell_chunk);
  const int local = row - b.row_first;
  thread_local std::vector<Segment> segments;
  thread_local std::vector<Segment> changed;
  thread_local std::vector<Segment> merged;
  segments.clear();
  changed.clear();
  merged.clear();
  // The intervals meeting [cell_begin - 1, cell_end + 1] are replaced.
  const Interval* row_first = rowBegin(b, local);
  const Interval* row_last = rowEnd(b, local);
  const Interval* from = firstEndingAtOrAfter(b, local, cell_begin - 1);
  const Interval* to = from;
  while (to != row_last && to->begin <= cell_end + 1) {
    ++to;
  }
  const auto first = static_cast<uint32_t>(from - row_first);
  const auto last = static_cast<uint32_t>(to - row_first);
  int lo = cell_begin;
  int hi = cell_end;
  for (const Interval* n = from; n != to; ++n) {
    lo = std::min(lo, n->begin);
    hi = std::max(hi, n->end);
  }
  int pos = lo;
  for (const Interval* n = from; n != to; ++n) {
    if (n->begin > pos) {
      segments.push_back({pos, n->begin - 1, 0, kNoNet});
    }
    segments.push_back(*n);
    pos = n->end + 1;
  }
  if (pos <= hi) {
    segments.push_back({pos, hi, 0, kNoNet});
  }
  bool ok = true;
  for (const Segment& s : segments) {
    if (s.end < cell_begin || s.begin > cell_end) {
      changed.push_back(s);
      continue;
    }
    if (s.begin < cell_begin) {
      changed.push_back({s.begin, cell_begin - 1, s.value, s.owner});
    }
    const Segment middle{std::max(s.begin, cell_begin),
                         std::min(s.end, cell_end),
                         s.value,
                         s.owner};
    ok &= applyToSegment(b, middle, piece, owner, insert, cache, changed);
    if (s.end > cell_end) {
      changed.push_back({cell_end + 1, s.end, s.value, s.owner});
    }
  }
  // Re-merge neighbouring cells with equal configuration and owner.
  Segment open{0, -2, 0, kNoNet};
  auto flush = [&]() {
    if (open.end >= open.begin && (open.value != 0 || open.owner != kNoNet)) {
      merged.push_back(open);
    }
  };
  for (const Segment& s : changed) {
    if (open.end + 1 == s.begin && open.value == s.value
        && open.owner == s.owner && s.owner != kSharedCell) {
      open.end = s.end;
      continue;
    }
    flush();
    open = s;
  }
  flush();
  replaceRow(b, local, first, last, merged);
  return ok;
}

bool ShapeGrid::update(const int layer_index,
                       const Rect& rect,
                       const ClassId cls,
                       const ShapeType type,
                       const int32_t owner,
                       const bool insert,
                       const uint8_t line_end_ends,
                       ShapeGridWriteCache& cache)
{
  const ShapeGridLayerGeometry& g = layers_[layer_index].geom;
  const Coord a_lo = g.horizontal ? rect.xlo : rect.ylo;
  const Coord a_hi = g.horizontal ? rect.xhi : rect.yhi;
  const Coord c_lo = g.horizontal ? rect.ylo : rect.xlo;
  const Coord c_hi = g.horizontal ? rect.yhi : rect.xhi;
  if (a_lo >= a_hi || c_lo >= c_hi) {
    return true;
  }
  const int64_t c0 = floorDivide(a_lo - g.along_origin, g.along_cell);
  const int64_t c1 = floorDivide(a_hi - 1 - g.along_origin, g.along_cell);
  const int64_t r0 = floorDivide(c_lo - g.cross_origin, g.cross_cell);
  const int64_t r1 = floorDivide(c_hi - 1 - g.cross_origin, g.cross_cell);
  if (c1 < 0 || c0 >= g.along_cells || r1 < 0 || r0 >= g.cross_cells) {
    return true;
  }
  // A removal is applied only if every piece is present.
  if (!insert
      && !apply(layer_index,
                rect,
                cls,
                type,
                owner,
                insert,
                line_end_ends,
                cache,
                /*dry_run=*/true)) {
    return false;
  }
  return apply(
      layer_index, rect, cls, type, owner, insert, line_end_ends, cache, false);
}

bool ShapeGrid::apply(const int layer_index,
                      const Rect& rect,
                      const ClassId cls,
                      const ShapeType type,
                      const int32_t owner,
                      const bool insert,
                      const uint8_t line_end_ends,
                      ShapeGridWriteCache& cache,
                      const bool dry_run)
{
  Layer& layer = layers_[layer_index];
  const ShapeGridLayerGeometry& g = layer.geom;
  const Coord a_lo = g.horizontal ? rect.xlo : rect.ylo;
  const Coord a_hi = g.horizontal ? rect.xhi : rect.yhi;
  const Coord c_lo = g.horizontal ? rect.ylo : rect.xlo;
  const Coord c_hi = g.horizontal ? rect.yhi : rect.xhi;
  const int64_t c0_raw = floorDivide(a_lo - g.along_origin, g.along_cell);
  const int64_t c1_raw = floorDivide(a_hi - 1 - g.along_origin, g.along_cell);
  const int64_t r0_raw = floorDivide(c_lo - g.cross_origin, g.cross_cell);
  const int64_t r1_raw = floorDivide(c_hi - 1 - g.cross_origin, g.cross_cell);
  const int c0 = static_cast<int>(std::max<int64_t>(c0_raw, 0));
  const int c1 = static_cast<int>(std::min<int64_t>(c1_raw, g.along_cells - 1));
  const int r0 = static_cast<int>(std::max<int64_t>(r0_raw, 0));
  const int r1 = static_cast<int>(std::min<int64_t>(r1_raw, g.cross_cells - 1));
  const bool line_end = classes_->get(cls).line_end;
  const bool low_end = line_end && (line_end_ends & 1) && c0 == c0_raw;
  const bool high_end = line_end && (line_end_ends & 2) && c1 == c1_raw;
  bool ok = true;
  for (int r = r0; r <= r1 && (ok || !dry_run); r++) {
    const Coord row_lo
        = static_cast<Coord>(g.cross_origin + int64_t(r) * g.cross_cell);
    const Coord pc_lo = std::max(c_lo, row_lo) - row_lo;
    const Coord pc_hi = std::min(c_hi, row_lo + g.cross_cell) - row_lo;
    auto run = [&](int first, int last) {
      if (first > last) {
        return;
      }
      const Coord cell_lo
          = static_cast<Coord>(g.along_origin + int64_t(first) * g.along_cell);
      CellShape piece;
      const Coord pa_lo = std::max(a_lo, cell_lo) - cell_lo;
      const Coord pa_hi = std::min(a_hi, cell_lo + g.along_cell) - cell_lo;
      if (g.horizontal) {
        piece.xlo = pa_lo;
        piece.xhi = pa_hi;
        piece.ylo = pc_lo;
        piece.yhi = pc_hi;
      } else {
        piece.ylo = pa_lo;
        piece.yhi = pa_hi;
        piece.xlo = pc_lo;
        piece.xhi = pc_hi;
      }
      piece.cls = cls;
      piece.type = type;
      piece.line_end_ends = static_cast<uint8_t>(
          (low_end && first == c0 ? 1 : 0) | (high_end && last == c1 ? 2 : 0));
      for (int begin = first; begin <= last;) {
        const int chunk = layer.cell_chunk[begin];
        const int end = std::min(last, layer.chunk_first_cell[chunk + 1] - 1);
        ok &= dry_run
                  ? containsRun(layer, r, begin, end, piece, owner)
                  : applyRun(layer, r, begin, end, piece, owner, insert, cache);
        begin = end + 1;
      }
    };
    run(c0, c0);
    run(c0 + 1, c1 - 1);
    if (c1 > c0) {
      run(c1, c1);
    }
  }
  return ok;
}

bool ShapeGrid::updateWide(const int layer_index,
                           const Rect& rect,
                           const ClassId cls,
                           const ShapeType type,
                           const int32_t owner,
                           const bool insert,
                           const uint8_t line_end_ends)
{
  Layer& layer = layers_[layer_index];
  const ShapeGridLayerGeometry& g = layer.geom;
  const int64_t a_lo = g.horizontal ? rect.xlo : rect.ylo;
  const int64_t a_hi = g.horizontal ? rect.xhi : rect.yhi;
  const int64_t c_lo = g.horizontal ? rect.ylo : rect.xlo;
  const int64_t c_hi = g.horizontal ? rect.yhi : rect.xhi;
  if (a_lo >= a_hi || c_lo >= c_hi) {
    return true;
  }
  // Clipped to the grid area like the cell pieces; a clipped end is no
  // line-end.
  const int64_t ga_lo = g.along_origin;
  const int64_t ga_hi = ga_lo + int64_t{g.along_cells} * g.along_cell;
  const int64_t gc_lo = g.cross_origin;
  const int64_t gc_hi = gc_lo + int64_t{g.cross_cells} * g.cross_cell;
  const auto pa_lo = static_cast<Coord>(std::max(a_lo, ga_lo));
  const auto pa_hi = static_cast<Coord>(std::min(a_hi, ga_hi));
  const auto pc_lo = static_cast<Coord>(std::max(c_lo, gc_lo));
  const auto pc_hi = static_cast<Coord>(std::min(c_hi, gc_hi));
  if (pa_lo >= pa_hi || pc_lo >= pc_hi) {
    return true;
  }
  WideShape w;
  GridPiece& p = w.piece;
  p.rect = g.horizontal ? Rect{pa_lo, pc_lo, pa_hi, pc_hi}
                        : Rect{pc_lo, pa_lo, pc_hi, pa_hi};
  p.cls = cls;
  p.type = type;
  p.owner = owner;
  if (classes_->get(cls).line_end) {
    p.line_end_ends = static_cast<uint8_t>(
        ((line_end_ends & 1) != 0 && a_lo >= ga_lo ? 1 : 0)
        | ((line_end_ends & 2) != 0 && a_hi <= ga_hi ? 2 : 0));
  }
  const auto cell0 = floorDivide(pa_lo - ga_lo, g.along_cell);
  const auto cell1 = floorDivide(pa_hi - 1 - ga_lo, g.along_cell);
  const auto row0 = floorDivide(pc_lo - gc_lo, g.cross_cell);
  const auto row1 = floorDivide(pc_hi - 1 - gc_lo, g.cross_cell);
  w.cell_chunk_lo = layer.cell_chunk[cell0];
  w.cell_chunk_hi = layer.cell_chunk[cell1];
  w.row_chunk_lo = layer.row_chunk[row0];
  w.row_chunk_hi = layer.row_chunk[row1];
  auto list = [&](int q, int k) -> WideList& {
    return layer.wide[static_cast<size_t>(q) * layer.along_chunks + k];
  };
  auto same = [&](const WideShape& o) {
    return o.piece.rect == p.rect && o.piece.cls == p.cls
           && o.piece.type == p.type && o.piece.owner == p.owner
           && o.piece.line_end_ends == p.line_end_ends;
  };
  if (!insert) {
    // A removal is applied only if the shape is present.
    for (int q = w.row_chunk_lo; q <= w.row_chunk_hi; q++) {
      for (int k = w.cell_chunk_lo; k <= w.cell_chunk_hi; k++) {
        WideList& l = list(q, k);
        std::lock_guard lock(l.mutex);
        if (std::ranges::none_of(l.shapes, same)) {
          return false;
        }
      }
    }
  }
  for (int q = w.row_chunk_lo; q <= w.row_chunk_hi; q++) {
    for (int k = w.cell_chunk_lo; k <= w.cell_chunk_hi; k++) {
      WideList& l = list(q, k);
      std::lock_guard lock(l.mutex);
      if (insert) {
        l.shapes.push_back(w);
      } else {
        const auto it = std::ranges::find_if(l.shapes, same);
        *it = l.shapes.back();
        l.shapes.pop_back();
      }
      l.size.store(static_cast<uint32_t>(l.shapes.size()),
                   std::memory_order_release);
    }
  }
  return true;
}

bool ShapeGrid::insert(const int layer,
                       const Rect& rect,
                       const ClassId cls,
                       const ShapeType type,
                       const int32_t owner,
                       ShapeGridWriteCache& cache,
                       const uint8_t line_end_ends)
{
  if (classes_->wide(cls)) {
    return updateWide(layer,
                      rect,
                      cls,
                      type,
                      isOwned(type) ? owner : kNoNet,
                      true,
                      line_end_ends);
  }
  return update(layer,
                rect,
                cls,
                type,
                isOwned(type) ? owner : kNoNet,
                true,
                line_end_ends,
                cache);
}

bool ShapeGrid::remove(const int layer,
                       const Rect& rect,
                       const ClassId cls,
                       const ShapeType type,
                       const int32_t owner,
                       ShapeGridWriteCache& cache,
                       const uint8_t line_end_ends)
{
  if (classes_->wide(cls)) {
    return updateWide(layer,
                      rect,
                      cls,
                      type,
                      isOwned(type) ? owner : kNoNet,
                      false,
                      line_end_ends);
  }
  return update(layer,
                rect,
                cls,
                type,
                isOwned(type) ? owner : kNoNet,
                false,
                line_end_ends,
                cache);
}

void ShapeGrid::emitPieces(const Layer& layer,
                           const int row,
                           const int begin,
                           const int end,
                           const uint32_t config,
                           const int32_t owner,
                           const int clip_lo,
                           const int clip_hi,
                           GridWindow& out) const
{
  const int cb = std::max(begin, clip_lo);
  const int ce = std::min(end, clip_hi);
  if (cb > ce) {
    return;
  }
  const ShapeGridLayerGeometry& g = layer.geom;
  const Coord row_lo
      = static_cast<Coord>(g.cross_origin + int64_t(row) * g.cross_cell);
  auto cell_lo = [&](int c) {
    return static_cast<Coord>(g.along_origin + int64_t(c) * g.along_cell);
  };
  for (const CellShape& s : configs_[config].shapes) {
    const Coord pa_lo = g.horizontal ? s.xlo : s.ylo;
    const Coord pa_hi = g.horizontal ? s.xhi : s.yhi;
    const Coord pc_lo = g.horizontal ? s.ylo : s.xlo;
    const Coord pc_hi = g.horizontal ? s.yhi : s.xhi;
    const int32_t piece_owner = isOwned(s.type) ? owner : kNoNet;
    auto push = [&](Coord along_lo, Coord along_hi) {
      GridPiece p;
      if (g.horizontal) {
        p.rect = {along_lo, row_lo + pc_lo, along_hi, row_lo + pc_hi};
      } else {
        p.rect = {row_lo + pc_lo, along_lo, row_lo + pc_hi, along_hi};
      }
      p.cls = s.cls;
      p.type = s.type;
      p.line_end_ends = s.line_end_ends;
      p.owner = piece_owner;
      out.pieces.push_back(p);
    };
    if (pa_lo == 0 && pa_hi == g.along_cell && s.line_end_ends == 0) {
      // The piece spans its cells: one rectangle for the whole run.
      push(cell_lo(cb), cell_lo(ce) + g.along_cell);
      continue;
    }
    for (int c = cb; c <= ce; c++) {
      push(cell_lo(c) + pa_lo, cell_lo(c) + pa_hi);
    }
  }
}

void ShapeGrid::collect(const int layer_index,
                        const Rect& window,
                        GridWindow& out) const
{
  out.clear();
  out.layer = layer_index;
  const Layer& layer = layers_[layer_index];
  const ShapeGridLayerGeometry& g = layer.geom;
  out.horizontal = g.horizontal;
  const Coord a_lo = g.horizontal ? window.xlo : window.ylo;
  const Coord a_hi = g.horizontal ? window.xhi : window.yhi;
  const Coord c_lo = g.horizontal ? window.ylo : window.xlo;
  const Coord c_hi = g.horizontal ? window.yhi : window.xhi;
  auto clampCell = [&](int64_t c) {
    return static_cast<int>(std::clamp<int64_t>(c, 0, g.along_cells - 1));
  };
  auto clampRow = [&](int64_t r) {
    return static_cast<int>(std::clamp<int64_t>(r, 0, g.cross_cells - 1));
  };
  const int c0 = clampCell(floorDivide(a_lo - g.along_origin, g.along_cell));
  const int c1 = clampCell(floorDivide(a_hi - g.along_origin, g.along_cell));
  const int r0 = clampRow(floorDivide(c_lo - g.cross_origin, g.cross_cell));
  const int r1 = clampRow(floorDivide(c_hi - g.cross_origin, g.cross_cell));
  const int k0 = layer.cell_chunk[c0];
  const int k1 = layer.cell_chunk[c1];
  const int horizon_lo
      = std::max(c0 - kKnownHorizon, layer.chunk_first_cell[k0]);
  const int horizon_hi
      = std::min(c1 + kKnownHorizon, layer.chunk_first_cell[k1 + 1] - 1);
  for (int r = r0; r <= r1; r++) {
    const int row_chunk = layer.row_chunk[r];
    int known_lo = horizon_lo;
    int known_hi = horizon_hi;
    for (int k = k0; k <= k1; k++) {
      const Block* b = block(layer, row_chunk, k);
      if (b == nullptr) {
        continue;
      }
      const int local = r - b->row_first;
      const int q_lo = std::max(c0, layer.chunk_first_cell[k]);
      const int q_hi = std::min(c1, layer.chunk_first_cell[k + 1] - 1);
      // found = [from, to): the intervals meeting [q_lo, q_hi].
      const Interval* from = firstEndingAtOrAfter(*b, local, q_lo);
      const Interval* row_last = rowEnd(*b, local);
      const Interval* to = from;
      for (; to != row_last && to->begin <= q_hi; ++to) {
        const Interval& node = *to;
        out.intervals++;
        if (node.owner == kSharedCell) {
          for (const SharedPart& p : b->shared[node.value]) {
            emitPieces(layer,
                       r,
                       node.begin,
                       node.end,
                       p.config,
                       p.owner,
                       horizon_lo,
                       horizon_hi,
                       out);
          }
        } else {
          emitPieces(layer,
                     r,
                     node.begin,
                     node.end,
                     node.value,
                     node.owner,
                     horizon_lo,
                     horizon_hi,
                     out);
        }
      }
      if (k == k0) {
        const int x = from != to && from->begin <= c0 ? from->begin : c0;
        const Interval* pred = findLast(*b, local, x - 1);
        if (pred != nullptr) {
          known_lo = std::max(known_lo, pred->end + 1);
        }
      }
      if (k == k1) {
        const int y = from != to && (to - 1)->end >= c1 ? (to - 1)->end : c1;
        const Interval* succ = findFirst(*b, local, y + 1);
        if (succ != nullptr) {
          known_hi = std::min(known_hi, succ->begin - 1);
        }
      }
    }
    // Nothing is stored beyond the grid area, so it is known to be empty.
    constexpr Coord kOutside = 1 << 30;
    KnownRow kr;
    kr.cross_lo
        = static_cast<Coord>(g.cross_origin + int64_t(r) * g.cross_cell);
    kr.cross_hi = kr.cross_lo + g.cross_cell;
    kr.along_lo
        = static_cast<Coord>(g.along_origin + int64_t(known_lo) * g.along_cell);
    kr.along_hi = static_cast<Coord>(g.along_origin
                                     + int64_t(known_hi + 1) * g.along_cell);
    if (r == 0) {
      kr.cross_lo = -kOutside;
    }
    if (r == g.cross_cells - 1) {
      kr.cross_hi = kOutside;
    }
    if (known_lo == 0) {
      kr.along_lo = -kOutside;
    }
    if (known_hi == g.along_cells - 1) {
      kr.along_hi = kOutside;
    }
    out.rows.push_back(kr);
  }
  // Wide shapes near the window and the known rows (which lie inside the
  // rows r0..r1 and the cells horizon_lo..horizon_hi).
  const Coord extra = std::max<Coord>(
      0,
      classes_->classReach(layer_index) - classes_->narrowReach(layer_index));
  const int64_t along_lo = std::min<int64_t>(
      a_lo, g.along_origin + int64_t{horizon_lo} * g.along_cell);
  const int64_t along_hi = std::max<int64_t>(
      a_hi, g.along_origin + int64_t{horizon_hi + 1} * g.along_cell);
  const int64_t cross_lo
      = std::min<int64_t>(c_lo, g.cross_origin + int64_t{r0} * g.cross_cell);
  const int64_t cross_hi = std::max<int64_t>(
      c_hi, g.cross_origin + int64_t{r1 + 1} * g.cross_cell);
  auto clampCoord = [](int64_t v) {
    return static_cast<Coord>(
        std::clamp<int64_t>(v,
                            std::numeric_limits<Coord>::min(),
                            std::numeric_limits<Coord>::max()));
  };
  const Rect area{clampCoord(along_lo - extra),
                  clampCoord(cross_lo - extra),
                  clampCoord(along_hi + extra),
                  clampCoord(cross_hi + extra)};
  collectWide(layer, area, out);
}

void ShapeGrid::collectWide(const Layer& layer,
                            const Rect& area,
                            GridWindow& out) const
{
  // `area` holds (along, cross) coordinates as (x, y), closed.
  const ShapeGridLayerGeometry& g = layer.geom;
  auto cellOf = [&](Coord a) {
    return static_cast<int>(std::clamp<int64_t>(
        floorDivide(int64_t{a} - g.along_origin, g.along_cell),
        0,
        g.along_cells - 1));
  };
  auto rowOf = [&](Coord c) {
    return static_cast<int>(std::clamp<int64_t>(
        floorDivide(int64_t{c} - g.cross_origin, g.cross_cell),
        0,
        g.cross_cells - 1));
  };
  const int k0 = layer.cell_chunk[cellOf(area.xlo)];
  const int k1 = layer.cell_chunk[cellOf(area.xhi)];
  const int q0 = layer.row_chunk[rowOf(area.ylo)];
  const int q1 = layer.row_chunk[rowOf(area.yhi)];
  for (int q = q0; q <= q1; q++) {
    for (int k = k0; k <= k1; k++) {
      const WideList& l
          = layer.wide[static_cast<size_t>(q) * layer.along_chunks + k];
      if (l.size.load(std::memory_order_acquire) == 0) {
        continue;
      }
      std::lock_guard lock(l.mutex);
      for (const WideShape& w : l.shapes) {
        // Once per shape: from the first of its chunks inside the area.
        if (k != std::max(w.cell_chunk_lo, k0)
            || q != std::max(w.row_chunk_lo, q0)) {
          continue;
        }
        const Rect& r = w.piece.rect;
        const Coord pa_lo = g.horizontal ? r.xlo : r.ylo;
        const Coord pa_hi = g.horizontal ? r.xhi : r.yhi;
        const Coord pc_lo = g.horizontal ? r.ylo : r.xlo;
        const Coord pc_hi = g.horizontal ? r.yhi : r.xhi;
        if (pa_hi < area.xlo || pa_lo > area.xhi || pc_hi < area.ylo
            || pc_lo > area.yhi) {
          continue;
        }
        out.pieces.push_back(w.piece);
      }
    }
  }
}

void ShapeGrid::forEachWideShape(
    const int layer_index,
    const std::function<void(const GridPiece&)>& fn) const
{
  const Layer& layer = layers_[layer_index];
  for (int q = 0; q < layer.cross_chunks; q++) {
    for (int k = 0; k < layer.along_chunks; k++) {
      const WideList& l
          = layer.wide[static_cast<size_t>(q) * layer.along_chunks + k];
      std::lock_guard lock(l.mutex);
      for (const WideShape& w : l.shapes) {
        if (w.cell_chunk_lo == k && w.row_chunk_lo == q) {
          fn(w.piece);
        }
      }
    }
  }
}

void ShapeGrid::compact()
{
  for (Layer& layer : layers_) {
    const size_t lists
        = static_cast<size_t>(layer.along_chunks) * layer.cross_chunks;
    for (size_t k = 0; k < lists; k++) {
      std::lock_guard lock(layer.wide[k].mutex);
      layer.wide[k].shapes.shrink_to_fit();
    }
    for (auto& b : layer.blocks) {
      if (b) {
        b->pool.shrink_to_fit();
        b->shared.shrink_to_fit();
      }
    }
  }
}

ShapeGrid::Statistics ShapeGrid::statistics() const
{
  Statistics s;
  s.configurations = configs_.size();
  for (size_t c = 0; c < configs_.size(); c++) {
    s.configuration_shapes += configs_[c].shapes.size();
  }
  // Table entries, their shape vectors and the hash index (one copy).
  s.bytes += s.configurations * (sizeof(CellConfiguration) + 48)
             + 2 * s.configuration_shapes * sizeof(CellShape);
  for (const Layer& layer : layers_) {
    s.bytes += (layer.cell_chunk.capacity() + layer.row_chunk.capacity()
                + layer.chunk_first_cell.capacity()
                + layer.chunk_first_row.capacity())
                   * sizeof(int)
               + layer.blocks.capacity() * sizeof(layer.blocks[0]);
    const size_t lists
        = static_cast<size_t>(layer.along_chunks) * layer.cross_chunks;
    s.bytes += lists * sizeof(WideList);
    for (size_t k = 0; k < lists; k++) {
      const WideList& l = layer.wide[k];
      std::lock_guard lock(l.mutex);
      s.bytes += l.shapes.capacity() * sizeof(WideShape);
      for (const WideShape& w : l.shapes) {
        s.wide_shapes
            += static_cast<size_t>(w.row_chunk_lo) * layer.along_chunks
                   + w.cell_chunk_lo
               == k;
      }
    }
    for (const auto& b : layer.blocks) {
      if (!b) {
        continue;
      }
      s.bytes += sizeof(Block) + b->row_offset.capacity() * sizeof(uint32_t)
                 + b->pool.capacity() * sizeof(Interval)
                 + b->shared.capacity() * sizeof(b->shared[0]);
      for (const auto& parts : b->shared) {
        s.bytes += parts.capacity() * sizeof(SharedPart);
      }
      s.chunks++;
      s.shared_cells
          += static_cast<int64_t>(b->shared.size() - b->free_shared.size());
      for (int local = 0; local < b->rows; local++) {
        for (const Interval* n = rowBegin(*b, local); n != rowEnd(*b, local);
             ++n) {
          s.intervals++;
          s.interval_cells += n->end - n->begin + 1;
        }
      }
    }
  }
  return s;
}

void ShapeGrid::forEachInterval(
    const int layer_index,
    const std::function<void(const IntervalView&)>& fn) const
{
  const Layer& layer = layers_[layer_index];
  for (int rc = 0; rc < layer.cross_chunks; rc++) {
    for (int k = 0; k < layer.along_chunks; k++) {
      const Block* b = block(layer, rc, k);
      if (b == nullptr) {
        continue;
      }
      for (int local = 0; local < b->rows; local++) {
        for (const Interval* n = rowBegin(*b, local); n != rowEnd(*b, local);
             ++n) {
          const Interval& node = *n;
          IntervalView v{b->row_first + local,
                         node.begin,
                         node.end,
                         k,
                         node.value,
                         node.owner,
                         {}};
          if (node.owner == kSharedCell) {
            for (const SharedPart& p : b->shared[node.value]) {
              v.parts.emplace_back(p.owner, p.config);
            }
          }
          fn(v);
        }
      }
    }
  }
}

Coord maximalCellSize(const ShapeClassTable& classes, const int layer)
{
  const LayerRules& rules = classes.layer(layer);
  Coord s = rules.minimumSpacing();
  if (rules.via_layer && rules.projection_spacing > 0) {
    s = std::min(s, rules.projection_spacing);
    if (layer >= 2) {
      s = std::min(s, classes.layer(layer - 2).minimumSpacing());
    }
  }
  if (s <= 0) {
    return 0;
  }
  // Largest d with 2 d^2 <= s^2: two points of one half-open d x d cell
  // are closer than s.
  int64_t d = static_cast<int64_t>(s / 1.41421356237);
  while (2 * d * d > static_cast<int64_t>(s) * s) {
    d--;
  }
  while (2 * (d + 1) * (d + 1) <= static_cast<int64_t>(s) * s) {
    d++;
  }
  return static_cast<Coord>(std::max<int64_t>(d, 1));
}

}  // namespace drt::br
