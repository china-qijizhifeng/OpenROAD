// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Shape grid (Sec. 3.3) against a brute-force model: after random insert /
// remove sequences every cell holds exactly the intersections of the live
// shapes with it (relative to its anchor, with type and class), intervals
// are maximal runs of equal configuration and owner inside a tile chunk,
// configuration 0 is never stored, and collect() reports exactly the
// content of the area it declares known.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <map>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include "BonnRouteRoutingSpaceTestUtil.h"
#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrWireModel.h"

using namespace drt::br;
using namespace drt::br::test;

namespace {

int64_t g_wide_shapes = 0;

using CellKey = std::tuple<int, int, int>;  // layer, row, cell
using Content = std::vector<std::pair<int32_t, CellShape>>;

// Cells of a shape and its pieces, as the shape grid should store them.
void addPieces(const ShapeGrid& grid,
               const StoredShape& s,
               std::map<CellKey, Content>& cells)
{
  const ShapeGridLayerGeometry& g = grid.geometry(s.layer);
  const Coord a_lo = g.horizontal ? s.rect.xlo : s.rect.ylo;
  const Coord a_hi = g.horizontal ? s.rect.xhi : s.rect.yhi;
  const Coord c_lo = g.horizontal ? s.rect.ylo : s.rect.xlo;
  const Coord c_hi = g.horizontal ? s.rect.yhi : s.rect.xhi;
  const bool line_end = grid.classes().get(s.cls).line_end;
  for (int r = 0; r < g.cross_cells; r++) {
    const Coord row_lo = g.cross_origin + r * g.cross_cell;
    const Coord pc_lo = std::max(c_lo, row_lo) - row_lo;
    const Coord pc_hi = std::min(c_hi, row_lo + g.cross_cell) - row_lo;
    if (pc_lo >= pc_hi) {
      continue;
    }
    for (int c = 0; c < g.along_cells; c++) {
      const Coord cell_lo = g.along_origin + c * g.along_cell;
      const Coord pa_lo = std::max(a_lo, cell_lo) - cell_lo;
      const Coord pa_hi = std::min(a_hi, cell_lo + g.along_cell) - cell_lo;
      if (pa_lo >= pa_hi) {
        continue;
      }
      CellShape p;
      if (g.horizontal) {
        p = {pa_lo, pc_lo, pa_hi, pc_hi, s.cls, s.type, 0};
      } else {
        p = {pc_lo, pa_lo, pc_hi, pa_hi, s.cls, s.type, 0};
      }
      if (line_end) {
        p.line_end_ends = static_cast<uint8_t>(
            (a_lo >= cell_lo && a_lo < cell_lo + g.along_cell ? 1 : 0)
            | (a_hi > cell_lo && a_hi <= cell_lo + g.along_cell ? 2 : 0));
      }
      cells[{s.layer, r, c}].emplace_back(isOwned(s.type) ? s.owner : kNoNet,
                                          p);
    }
  }
}

Content sorted(Content c)
{
  std::ranges::sort(c);
  return c;
}

using WideKey = std::tuple<Coord, Coord, Coord, Coord, int, int, int, int32_t>;

WideKey wideKey(const GridPiece& p)
{
  return {p.rect.xlo,
          p.rect.ylo,
          p.rect.xhi,
          p.rect.yhi,
          p.cls,
          static_cast<int>(p.type),
          p.line_end_ends,
          p.owner};
}

// A wide shape as the index stores it: clipped to the grid area, owner
// only for owned types, line-ends only at unclipped ends; false if nothing
// is left.
bool widePiece(const ShapeGrid& grid, const StoredShape& s, GridPiece& out)
{
  const ShapeGridLayerGeometry& g = grid.geometry(s.layer);
  const Coord ga_lo = g.along_origin;
  const Coord ga_hi = g.along_origin + g.along_cells * g.along_cell;
  const Coord gc_lo = g.cross_origin;
  const Coord gc_hi = g.cross_origin + g.cross_cells * g.cross_cell;
  const Coord a_lo = g.horizontal ? s.rect.xlo : s.rect.ylo;
  const Coord a_hi = g.horizontal ? s.rect.xhi : s.rect.yhi;
  const Coord c_lo = g.horizontal ? s.rect.ylo : s.rect.xlo;
  const Coord c_hi = g.horizontal ? s.rect.yhi : s.rect.xhi;
  const Coord pa_lo = std::max(a_lo, ga_lo);
  const Coord pa_hi = std::min(a_hi, ga_hi);
  const Coord pc_lo = std::max(c_lo, gc_lo);
  const Coord pc_hi = std::min(c_hi, gc_hi);
  if (pa_lo >= pa_hi || pc_lo >= pc_hi) {
    return false;
  }
  out.rect = g.horizontal ? Rect{pa_lo, pc_lo, pa_hi, pc_hi}
                          : Rect{pc_lo, pa_lo, pc_hi, pa_hi};
  out.cls = s.cls;
  out.type = s.type;
  out.owner = isOwned(s.type) ? s.owner : kNoNet;
  out.line_end_ends = 0;
  if (grid.classes().get(s.cls).line_end) {
    out.line_end_ends = static_cast<uint8_t>((a_lo >= ga_lo ? 1 : 0)
                                             | (a_hi <= ga_hi ? 2 : 0));
  }
  return true;
}

int runCase(std::mt19937& rng, const int case_id)
{
  ShapeClassTable classes;
  const int num_z = randomInt(rng, 1, 3);
  classes.setWideRatioPercent(randomWideRatio(rng));
  classes.setLayers(randomRules(rng, num_z, rng() % 2));
  const Coord extent = randomInt(rng, 120, 420);
  ShapeGrid grid;
  ChunkAxis cx{0, randomInt(rng, 30, 200), 0};
  cx.count = extent / cx.size + 1;
  ChunkAxis cy{0, randomInt(rng, 30, 200), 0};
  cy.count = extent / cy.size + 1;
  grid.init(&classes, gridGeometry(rng, classes, extent, rng() % 2), cx, cy);
  const int nets = randomInt(rng, 1, 5);
  for (int n = 0; n < nets; n++) {
    if (rng() % 3 == 0) {
      grid.setNetRipupLevel(n, static_cast<RipupLevel>(randomInt(rng, 1, 6)));
    }
  }
  ShapeGridWriteCache cache;
  std::vector<StoredShape> live;
  int fails = 0;
  const int steps = randomInt(rng, 10, 120);
  for (int step = 0; step < steps; step++) {
    if (live.empty() || rng() % 3 != 0) {
      const StoredShape s
          = randomShape(rng,
                        classes,
                        randomInt(rng, 0, classes.numLayers() - 1),
                        extent,
                        nets);
      grid.insert(s.layer, s.rect, s.cls, s.type, s.owner, cache);
      live.push_back(s);
    } else {
      const size_t at = rng() % live.size();
      const StoredShape s = live[at];
      live.erase(live.begin() + at);
      if (!grid.remove(s.layer, s.rect, s.cls, s.type, s.owner, cache)) {
        fails++;
        std::printf("case %d: remove failed\n", case_id);
      }
    }
  }
  // A removal that does not match an insertion is rejected.
  if (!live.empty()) {
    StoredShape s = live.front();
    s.rect.xhi += 1;
    if (grid.remove(s.layer, s.rect, s.cls, s.type, s.owner, cache)) {
      std::printf("case %d: bogus removal accepted\n", case_id);
      fails++;
    }
  }
  // Narrow shapes live in the cells, wide ones in the wide index.
  std::map<CellKey, Content> expected;
  std::vector<std::vector<GridPiece>> wide_live(classes.numLayers());
  for (const StoredShape& s : live) {
    if (classes.wide(s.cls)) {
      GridPiece p;
      if (widePiece(grid, s, p)) {
        wide_live[s.layer].push_back(p);
      }
      continue;
    }
    addPieces(grid, s, expected);
  }
  int64_t wide_count = 0;
  for (int layer = 0; layer < classes.numLayers(); layer++) {
    std::vector<WideKey> want;
    std::vector<WideKey> have;
    for (const GridPiece& p : wide_live[layer]) {
      want.push_back(wideKey(p));
    }
    grid.forEachWideShape(
        layer, [&](const GridPiece& p) { have.push_back(wideKey(p)); });
    std::ranges::sort(want);
    std::ranges::sort(have);
    wide_count += static_cast<int64_t>(want.size());
    if (want != have) {
      std::printf("case %d: wide index of layer %d differs (%zu vs %zu)\n",
                  case_id,
                  layer,
                  want.size(),
                  have.size());
      fails++;
    }
  }
  // Stored content per cell and interval invariants.
  std::map<CellKey, Content> stored;
  int64_t shared = 0;
  for (int layer = 0; layer < classes.numLayers(); layer++) {
    std::map<std::pair<int, int>, std::vector<ShapeGrid::IntervalView>> rows;
    grid.forEachInterval(layer, [&](const ShapeGrid::IntervalView& v) {
      rows[{v.row, v.chunk}].push_back(v);
      for (int c = v.begin; c <= v.end; c++) {
        if (grid.chunkOfCell(layer, c) != v.chunk) {
          std::printf("case %d: interval crosses a chunk\n", case_id);
          fails++;
        }
        Content& content = stored[{layer, v.row, c}];
        if (v.owner == ShapeGrid::kSharedCell) {
          shared++;
          if (v.begin != v.end || v.parts.size() < 2) {
            std::printf("case %d: bad shared cell\n", case_id);
            fails++;
          }
          int owned_parts = 0;
          for (const auto& [owner, config] : v.parts) {
            owned_parts += owner != kNoNet;
            for (const CellShape& s : grid.configurations()[config].shapes) {
              content.emplace_back(isOwned(s.type) ? owner : kNoNet, s);
            }
          }
          if (owned_parts < 2) {
            std::printf("case %d: shared cell with one owner\n", case_id);
            fails++;
          }
          continue;
        }
        if (v.config == 0) {
          std::printf("case %d: configuration 0 stored\n", case_id);
          fails++;
        }
        const CellConfiguration& cfg = grid.configurations()[v.config];
        if ((v.owner != kNoNet) != cfg.has_owned) {
          std::printf("case %d: owner without owned shapes\n", case_id);
          fails++;
        }
        for (const CellShape& s : cfg.shapes) {
          content.emplace_back(isOwned(s.type) ? v.owner : kNoNet, s);
        }
      }
    });
    for (auto& [key, list] : rows) {
      for (size_t k = 1; k < list.size(); k++) {
        const auto& a = list[k - 1];
        const auto& b = list[k];
        if (b.begin <= a.end) {
          std::printf("case %d: overlapping intervals\n", case_id);
          fails++;
        }
        if (b.begin == a.end + 1 && a.owner != ShapeGrid::kSharedCell
            && a.owner == b.owner && a.config == b.config) {
          std::printf("case %d: non-maximal intervals\n", case_id);
          fails++;
        }
      }
    }
  }
  for (auto& [key, content] : expected) {
    content = sorted(content);
  }
  for (auto& [key, content] : stored) {
    content = sorted(content);
  }
  if (expected != stored) {
    std::printf("case %d: cell contents differ (%zu vs %zu cells)\n",
                case_id,
                expected.size(),
                stored.size());
    fails++;
    int shown = 0;
    for (const auto& [key, content] : expected) {
      auto it = stored.find(key);
      if ((it == stored.end() || it->second != content) && shown++ < 3) {
        const auto [layer, r, c] = key;
        std::printf("  cell layer %d row %d cell %d: expected %zu stored %zu\n",
                    layer,
                    r,
                    c,
                    content.size(),
                    it == stored.end() ? size_t{0} : it->second.size());
        for (const auto& [owner, p] : content) {
          std::printf(
              "    want owner %d (%d,%d)-(%d,%d) cls %d type %d ends %d\n",
              owner,
              p.xlo,
              p.ylo,
              p.xhi,
              p.yhi,
              p.cls,
              (int) p.type,
              p.line_end_ends);
        }
        if (it != stored.end()) {
          for (const auto& [owner, p] : it->second) {
            std::printf(
                "    have owner %d (%d,%d)-(%d,%d) cls %d type %d ends %d\n",
                owner,
                p.xlo,
                p.ylo,
                p.xhi,
                p.yhi,
                p.cls,
                (int) p.type,
                p.line_end_ends);
          }
        }
      }
    }
  }
  // collect(): pieces cover exactly the content of the known area.
  for (int q = 0; q < 30; q++) {
    const int layer = randomInt(rng, 0, classes.numLayers() - 1);
    const Coord x = randomInt(rng, -20, extent + 20);
    const Coord y = randomInt(rng, -20, extent + 20);
    const Rect window{
        x, y, x + randomInt(rng, 0, 80), y + randomInt(rng, 0, 80)};
    GridWindow out;
    grid.collect(layer, window, out);
    const ShapeGridLayerGeometry& g = grid.geometry(layer);
    // Wide pieces: live ones, each once, all of those within reach -
    // narrowReach of the window or of a known row (closed).
    const Coord extra = classes.classReach(layer) - classes.narrowReach(layer);
    std::vector<WideKey> wide_got;
    for (const GridPiece& p : out.pieces) {
      if (classes.wide(p.cls)) {
        wide_got.push_back(wideKey(p));
      }
    }
    std::ranges::sort(wide_got);
    std::vector<WideKey> wide_all;
    for (const GridPiece& p : wide_live[layer]) {
      wide_all.push_back(wideKey(p));
    }
    std::ranges::sort(wide_all);
    if (!std::ranges::includes(wide_all, wide_got)) {
      std::printf(
          "case %d query %d: collect reports a wide shape that is "
          "not stored (or twice)\n",
          case_id,
          q);
      fails++;
    }
    auto meets = [&](const Rect& r, const Rect& a) {
      return r.xhi >= a.xlo && r.xlo <= a.xhi && r.yhi >= a.ylo
             && r.ylo <= a.yhi;
    };
    std::vector<Rect> needed{window.bloated(extra, extra)};
    for (const KnownRow& row : out.rows) {
      const Coord lo = std::max(row.along_lo, -(1 << 28));
      const Coord hi = std::min(row.along_hi, 1 << 28);
      const Coord clo = std::max(row.cross_lo, -(1 << 28));
      const Coord chi = std::min(row.cross_hi, 1 << 28);
      const Rect known = g.horizontal ? Rect{lo, clo, hi - 1, chi - 1}
                                      : Rect{clo, lo, chi - 1, hi - 1};
      needed.push_back(known.bloated(extra, extra));
    }
    for (const GridPiece& p : wide_live[layer]) {
      const bool need = std::ranges::any_of(
          needed, [&](const Rect& a) { return meets(p.rect, a); });
      if (need && !std::ranges::binary_search(wide_got, wideKey(p))) {
        std::printf(
            "case %d query %d: collect misses a wide shape\n", case_id, q);
        fails++;
        break;
      }
    }
    std::map<CellKey, Content> got;
    for (const GridPiece& p : out.pieces) {
      if (classes.wide(p.cls)) {
        continue;
      }
      StoredShape s{layer, p.rect, p.cls, p.type, p.owner};
      std::map<CellKey, Content> one;
      addPieces(grid, s, one);
      for (auto& [key, content] : one) {
        for (auto [owner, piece] : content) {
          // Run rectangles carry no line-end; per-cell pieces keep theirs.
          piece.line_end_ends = 0;
          got[key].emplace_back(owner, piece);
        }
      }
    }
    std::map<CellKey, Content> want;
    for (const KnownRow& row : out.rows) {
      for (int r = 0; r < g.cross_cells; r++) {
        const Coord row_lo = g.cross_origin + r * g.cross_cell;
        if (row_lo < row.cross_lo || row_lo + g.cross_cell > row.cross_hi) {
          continue;
        }
        for (int c = 0; c < g.along_cells; c++) {
          const Coord cell_lo = g.along_origin + c * g.along_cell;
          if (cell_lo < row.along_lo || cell_lo + g.along_cell > row.along_hi) {
            continue;
          }
          auto it = expected.find({layer, r, c});
          if (it != expected.end()) {
            for (auto [owner, piece] : it->second) {
              piece.line_end_ends = 0;
              want[{layer, r, c}].emplace_back(owner, piece);
            }
          }
        }
      }
    }
    for (auto& [key, content] : got) {
      content = sorted(content);
    }
    for (auto& [key, content] : want) {
      content = sorted(content);
    }
    if (got != want) {
      std::printf("case %d query %d: collect differs from the known area\n",
                  case_id,
                  q);
      fails++;
    }
  }
  const ShapeGrid::Statistics st = grid.statistics();
  if (st.shared_cells != shared) {
    std::printf("case %d: shared cell count\n", case_id);
    fails++;
  }
  if (st.wide_shapes != wide_count) {
    std::printf("case %d: wide shape count\n", case_id);
    fails++;
  }
  g_wide_shapes += wide_count;
  return fails;
}

}  // namespace

int main()
{
  std::mt19937 rng(20261111);
  int fails = 0;
  const int cases = 300;
  for (int c = 0; c < cases && fails < 10; c++) {
    fails += runCase(rng, c);
  }
  std::printf(
      "shape grid: %d cases, %d failures, %lld wide shapes at the end\n",
      cases,
      fails,
      static_cast<long long>(g_wide_shapes));
  return fails != 0;
}
