// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrShapeGrid.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
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
    b->roots.assign(b->rows, -1);
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

// AVL tree (Sec. 3.3) keyed by the first cell of an interval.

int32_t ShapeGrid::height(const Block& b, const int32_t n)
{
  return n < 0 ? 0 : b.nodes[n].height;
}

void ShapeGrid::update(Block& b, const int32_t n)
{
  b.nodes[n].height
      = 1 + std::max(height(b, b.nodes[n].left), height(b, b.nodes[n].right));
}

int32_t ShapeGrid::rotateLeft(Block& b, const int32_t n)
{
  const int32_t r = b.nodes[n].right;
  b.nodes[n].right = b.nodes[r].left;
  b.nodes[r].left = n;
  update(b, n);
  update(b, r);
  return r;
}

int32_t ShapeGrid::rotateRight(Block& b, const int32_t n)
{
  const int32_t l = b.nodes[n].left;
  b.nodes[n].left = b.nodes[l].right;
  b.nodes[l].right = n;
  update(b, n);
  update(b, l);
  return l;
}

int32_t ShapeGrid::rebalance(Block& b, const int32_t n)
{
  update(b, n);
  const int32_t l = b.nodes[n].left;
  const int32_t r = b.nodes[n].right;
  const int balance = height(b, l) - height(b, r);
  if (balance > 1) {
    if (height(b, b.nodes[l].left) < height(b, b.nodes[l].right)) {
      b.nodes[n].left = rotateLeft(b, l);
    }
    return rotateRight(b, n);
  }
  if (balance < -1) {
    if (height(b, b.nodes[r].right) < height(b, b.nodes[r].left)) {
      b.nodes[n].right = rotateRight(b, r);
    }
    return rotateLeft(b, n);
  }
  return n;
}

int32_t ShapeGrid::insertNode(Block& b, const int32_t root, const int32_t n)
{
  if (root < 0) {
    return n;
  }
  if (b.nodes[n].begin < b.nodes[root].begin) {
    const int32_t l = insertNode(b, b.nodes[root].left, n);
    b.nodes[root].left = l;
  } else {
    const int32_t r = insertNode(b, b.nodes[root].right, n);
    b.nodes[root].right = r;
  }
  return rebalance(b, root);
}

int32_t ShapeGrid::eraseMin(Block& b, const int32_t root, int32_t& min)
{
  if (b.nodes[root].left < 0) {
    min = root;
    return b.nodes[root].right;
  }
  const int32_t l = eraseMin(b, b.nodes[root].left, min);
  b.nodes[root].left = l;
  return rebalance(b, root);
}

int32_t ShapeGrid::eraseNode(Block& b, const int32_t root, const int32_t begin)
{
  if (root < 0) {
    return -1;
  }
  if (begin < b.nodes[root].begin) {
    const int32_t l = eraseNode(b, b.nodes[root].left, begin);
    b.nodes[root].left = l;
  } else if (begin > b.nodes[root].begin) {
    const int32_t r = eraseNode(b, b.nodes[root].right, begin);
    b.nodes[root].right = r;
  } else {
    const int32_t l = b.nodes[root].left;
    const int32_t r = b.nodes[root].right;
    b.free_nodes.push_back(root);
    if (r < 0) {
      return l;
    }
    int32_t min = -1;
    const int32_t rest = eraseMin(b, r, min);
    b.nodes[min].left = l;
    b.nodes[min].right = rest;
    return rebalance(b, min);
  }
  return rebalance(b, root);
}

int32_t ShapeGrid::newNode(Block& b, const Segment& s)
{
  int32_t n;
  if (!b.free_nodes.empty()) {
    n = b.free_nodes.back();
    b.free_nodes.pop_back();
  } else {
    n = static_cast<int32_t>(b.nodes.size());
    b.nodes.emplace_back();
  }
  b.nodes[n] = Node{s.begin, s.end, s.value, s.owner, -1, -1, 1};
  return n;
}

int32_t ShapeGrid::findLast(const Block& b, int32_t root, const int32_t cell)
{
  int32_t best = -1;
  while (root >= 0) {
    if (b.nodes[root].begin <= cell) {
      best = root;
      root = b.nodes[root].right;
    } else {
      root = b.nodes[root].left;
    }
  }
  return best;
}

int32_t ShapeGrid::findFirst(const Block& b, int32_t root, const int32_t cell)
{
  int32_t best = -1;
  while (root >= 0) {
    if (b.nodes[root].begin >= cell) {
      best = root;
      root = b.nodes[root].left;
    } else {
      root = b.nodes[root].right;
    }
  }
  return best;
}

void ShapeGrid::collectNodes(const Block& b,
                             const int32_t root,
                             const int32_t lo,
                             const int32_t hi,
                             std::vector<int32_t>& out)
{
  if (root < 0) {
    return;
  }
  const Node& n = b.nodes[root];
  if (n.begin > hi) {
    collectNodes(b, n.left, lo, hi, out);
    return;
  }
  if (n.end < lo) {
    collectNodes(b, n.right, lo, hi, out);
    return;
  }
  collectNodes(b, n.left, lo, hi, out);
  out.push_back(root);
  collectNodes(b, n.right, lo, hi, out);
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
  thread_local std::vector<int32_t> found;
  found.clear();
  collectNodes(*b, b->roots[row - b->row_first], cell_begin, cell_end, found);
  int next = cell_begin;
  const int32_t target = isOwned(piece.type) ? owner : kNoNet;
  for (const int32_t n : found) {
    const Node& node = b->nodes[n];
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
  int32_t root = b.roots[local];
  thread_local std::vector<int32_t> found;
  thread_local std::vector<Segment> segments;
  thread_local std::vector<Segment> changed;
  found.clear();
  segments.clear();
  changed.clear();
  collectNodes(b, root, cell_begin - 1, cell_end + 1, found);
  int lo = cell_begin;
  int hi = cell_end;
  for (const int32_t n : found) {
    lo = std::min(lo, b.nodes[n].begin);
    hi = std::max(hi, b.nodes[n].end);
  }
  int pos = lo;
  for (const int32_t n : found) {
    const Node& node = b.nodes[n];
    if (node.begin > pos) {
      segments.push_back({pos, node.begin - 1, 0, kNoNet});
    }
    segments.push_back({node.begin, node.end, node.value, node.owner});
    pos = node.end + 1;
  }
  if (pos <= hi) {
    segments.push_back({pos, hi, 0, kNoNet});
  }
  for (const Segment& s : segments) {
    if (s.value != 0 || s.owner != kNoNet) {
      root = eraseNode(b, root, s.begin);
    }
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
      root = insertNode(b, root, newNode(b, open));
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
  b.roots[local] = root;
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
  const int64_t c0 = floorDiv(a_lo - g.along_origin, g.along_cell);
  const int64_t c1 = floorDiv(a_hi - 1 - g.along_origin, g.along_cell);
  const int64_t r0 = floorDiv(c_lo - g.cross_origin, g.cross_cell);
  const int64_t r1 = floorDiv(c_hi - 1 - g.cross_origin, g.cross_cell);
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
  const int64_t c0_raw = floorDiv(a_lo - g.along_origin, g.along_cell);
  const int64_t c1_raw = floorDiv(a_hi - 1 - g.along_origin, g.along_cell);
  const int64_t r0_raw = floorDiv(c_lo - g.cross_origin, g.cross_cell);
  const int64_t r1_raw = floorDiv(c_hi - 1 - g.cross_origin, g.cross_cell);
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

bool ShapeGrid::insert(const int layer,
                       const Rect& rect,
                       const ClassId cls,
                       const ShapeType type,
                       const int32_t owner,
                       ShapeGridWriteCache& cache,
                       const uint8_t line_end_ends)
{
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
  const int c0 = clampCell(floorDiv(a_lo - g.along_origin, g.along_cell));
  const int c1 = clampCell(floorDiv(a_hi - g.along_origin, g.along_cell));
  const int r0 = clampRow(floorDiv(c_lo - g.cross_origin, g.cross_cell));
  const int r1 = clampRow(floorDiv(c_hi - g.cross_origin, g.cross_cell));
  const int k0 = layer.cell_chunk[c0];
  const int k1 = layer.cell_chunk[c1];
  thread_local std::vector<int32_t> found;
  for (int r = r0; r <= r1; r++) {
    const int row_chunk = layer.row_chunk[r];
    const int horizon_lo
        = std::max(c0 - kKnownHorizon, layer.chunk_first_cell[k0]);
    const int horizon_hi
        = std::min(c1 + kKnownHorizon, layer.chunk_first_cell[k1 + 1] - 1);
    int known_lo = horizon_lo;
    int known_hi = horizon_hi;
    for (int k = k0; k <= k1; k++) {
      const Block* b = block(layer, row_chunk, k);
      if (b == nullptr) {
        continue;
      }
      const int32_t root = b->roots[r - b->row_first];
      const int q_lo = std::max(c0, layer.chunk_first_cell[k]);
      const int q_hi = std::min(c1, layer.chunk_first_cell[k + 1] - 1);
      found.clear();
      collectNodes(*b, root, q_lo, q_hi, found);
      for (const int32_t n : found) {
        const Node& node = b->nodes[n];
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
        const int x = !found.empty() && b->nodes[found.front()].begin <= c0
                          ? b->nodes[found.front()].begin
                          : c0;
        const int32_t pred = findLast(*b, root, x - 1);
        if (pred >= 0) {
          known_lo = std::max(known_lo, b->nodes[pred].end + 1);
        }
      }
      if (k == k1) {
        const int y = !found.empty() && b->nodes[found.back()].end >= c1
                          ? b->nodes[found.back()].end
                          : c1;
        const int32_t succ = findFirst(*b, root, y + 1);
        if (succ >= 0) {
          known_hi = std::min(known_hi, b->nodes[succ].begin - 1);
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
}

void ShapeGrid::compact()
{
  for (Layer& layer : layers_) {
    for (auto& b : layer.blocks) {
      if (b) {
        b->nodes.shrink_to_fit();
        b->free_nodes.shrink_to_fit();
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
    for (const auto& b : layer.blocks) {
      if (!b) {
        continue;
      }
      s.bytes += sizeof(Block) + b->roots.capacity() * sizeof(int32_t)
                 + b->nodes.capacity() * sizeof(Node)
                 + b->free_nodes.capacity() * sizeof(int32_t)
                 + b->shared.capacity() * sizeof(b->shared[0]);
      for (const auto& parts : b->shared) {
        s.bytes += parts.capacity() * sizeof(SharedPart);
      }
      s.chunks++;
      s.intervals
          += static_cast<int64_t>(b->nodes.size() - b->free_nodes.size());
      s.shared_cells
          += static_cast<int64_t>(b->shared.size() - b->free_shared.size());
      for (const int32_t root : b->roots) {
        std::vector<int32_t> all;
        collectNodes(*b, root, b->cell_first, b->cell_first + b->cells, all);
        for (const int32_t n : all) {
          s.interval_cells += b->nodes[n].end - b->nodes[n].begin + 1;
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
        std::vector<int32_t> all;
        collectNodes(
            *b, b->roots[local], b->cell_first, b->cell_first + b->cells, all);
        for (const int32_t n : all) {
          const Node& node = b->nodes[n];
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
