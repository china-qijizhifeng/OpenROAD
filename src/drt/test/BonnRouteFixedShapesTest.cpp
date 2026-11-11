// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// insertFixedShapes (tile-row bands in parallel) against inserting the same
// shapes one by one: every layer's intervals (row, cells, chunk, owner,
// configuration content, shared-cell parts), its wide shapes, and the
// pieces collect() reports for random windows are equal. Random rules,
// cell geometries, tile sizes, band counts, shape types (blockages, pins,
// wide shapes, shapes reaching beyond the grid) and std::thread workers.

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <random>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include "BonnRouteRoutingSpaceTestUtil.h"
#include "br/BrFixedShapes.h"
#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrWireModel.h"

using namespace drt::br;
using namespace drt::br::test;

namespace {

// Runs body(k) for k < count on `threads` std::threads (dynamic order).
BandParallelFor threadPool(const int threads)
{
  return [threads](const int count, const std::function<void(int)>& body) {
    std::atomic<int> next{0};
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++) {
      pool.emplace_back([&]() {
        for (int k = next++; k < count; k = next++) {
          body(k);
        }
      });
    }
    for (std::thread& t : pool) {
      t.join();
    }
  };
}

using IntervalKey
    = std::tuple<int,
                 int,
                 int,
                 int,
                 int32_t,
                 std::vector<CellShape>,
                 std::vector<std::pair<int32_t, std::vector<CellShape>>>>;

// Intervals of a layer with configuration contents (ids may differ).
std::vector<IntervalKey> intervals(const ShapeGrid& grid, const int layer)
{
  std::vector<IntervalKey> out;
  grid.forEachInterval(layer, [&](const ShapeGrid::IntervalView& v) {
    std::vector<std::pair<int32_t, std::vector<CellShape>>> parts;
    for (const auto& [owner, config] : v.parts) {
      parts.emplace_back(owner, grid.configurations()[config].shapes);
    }
    // A shared cell keeps its content in the parts.
    const bool shared = v.owner == ShapeGrid::kSharedCell;
    out.emplace_back(v.row,
                     v.begin,
                     v.end,
                     v.chunk,
                     v.owner,
                     shared ? std::vector<CellShape>{}
                            : grid.configurations()[v.config].shapes,
                     parts);
  });
  return out;
}

using PieceKey
    = std::tuple<Coord, Coord, Coord, Coord, ClassId, int, int, int32_t>;

std::vector<PieceKey> pieces(const std::vector<GridPiece>& in)
{
  std::vector<PieceKey> out;
  for (const GridPiece& p : in) {
    out.emplace_back(p.rect.xlo,
                     p.rect.ylo,
                     p.rect.xhi,
                     p.rect.yhi,
                     p.cls,
                     static_cast<int>(p.type),
                     p.line_end_ends,
                     p.owner);
  }
  return out;
}

int runCase(std::mt19937& rng, const int case_id, int64_t& shapes_total)
{
  ShapeClassTable classes;
  const int num_z = randomInt(rng, 1, 3);
  classes.setWideRatioPercent(randomWideRatio(rng));
  classes.setLayers(randomRules(rng, num_z, rng() % 2));
  const Coord extent = randomInt(rng, 120, 420);
  ChunkAxis cx{0, randomInt(rng, 20, 200), 0};
  cx.count = extent / cx.size + 1;
  ChunkAxis cy{0, randomInt(rng, 20, 200), 0};
  cy.count = extent / cy.size + 1;
  const std::vector<ShapeGridLayerGeometry> geometry
      = gridGeometry(rng, classes, extent, rng() % 2);
  // Shapes per layer in insertion order; some reach beyond the grid.
  std::vector<std::vector<FixedShapeInsert>> inserts(classes.numLayers());
  const int count = randomInt(rng, 10, 300);
  for (int k = 0; k < count; k++) {
    const int layer = randomInt(rng, 0, classes.numLayers() - 1);
    StoredShape s = randomShape(rng, classes, layer, extent, 4);
    if (s.type != ShapeType::kBlockage && s.type != ShapeType::kPin) {
      s.type = rng() % 2 ? ShapeType::kBlockage : ShapeType::kPin;
      s.owner = s.type == ShapeType::kPin ? randomInt(rng, 0, 3) : kNoNet;
    }
    if (rng() % 10 == 0) {
      s.rect.xlo -= randomInt(rng, 0, 60);
      s.rect.yhi += randomInt(rng, 0, 60);
      s.cls = classes.classify(layer, s.rect, ShapeRole::kFixed);
    }
    inserts[layer].push_back({s.rect, s.cls, s.type, s.owner});
  }
  ShapeGrid serial;
  serial.init(&classes, geometry, cx, cy);
  ShapeGridWriteCache cache;
  for (int g = 0; g < classes.numLayers(); g++) {
    for (const FixedShapeInsert& in : inserts[g]) {
      serial.insert(g, in.rect, in.cls, in.type, in.owner, cache);
    }
  }
  ShapeGrid banded;
  banded.init(&classes, geometry, cx, cy);
  const int bands = randomInt(rng, 1, 12);
  const int threads = randomInt(rng, 1, 6);
  insertFixedShapes(banded, inserts, cx, cy, bands, threadPool(threads));
  int fails = 0;
  for (int g = 0; g < classes.numLayers(); g++) {
    shapes_total += static_cast<int64_t>(inserts[g].size());
    if (intervals(serial, g) != intervals(banded, g)) {
      fails++;
      std::printf(
          "case %d layer %d: intervals differ (bands %d)\n", case_id, g, bands);
    }
    std::vector<PieceKey> wide_a;
    std::vector<PieceKey> wide_b;
    serial.forEachWideShape(
        g, [&](const GridPiece& p) { wide_a.push_back(pieces({p}).front()); });
    banded.forEachWideShape(
        g, [&](const GridPiece& p) { wide_b.push_back(pieces({p}).front()); });
    if (wide_a != wide_b) {
      fails++;
      std::printf("case %d layer %d: wide shapes differ\n", case_id, g);
    }
    for (int q = 0; q < 20; q++) {
      const Coord x = randomInt(rng, -20, extent);
      const Coord y = randomInt(rng, -20, extent);
      const Rect window{
          x, y, x + randomInt(rng, 0, 80), y + randomInt(rng, 0, 80)};
      GridWindow a;
      GridWindow b;
      serial.collect(g, window, a);
      banded.collect(g, window, b);
      if (pieces(a.pieces) != pieces(b.pieces)) {
        fails++;
        std::printf("case %d layer %d: collect differs\n", case_id, g);
        break;
      }
    }
  }
  return fails;
}

}  // namespace

int main()
{
  std::mt19937 rng(20261112);
  int fails = 0;
  int64_t shapes = 0;
  const int cases = 400;
  for (int c = 0; c < cases; c++) {
    fails += runCase(rng, c, shapes);
  }
  std::printf(
      "fixed shapes: %d cases, %lld shapes, banded == serial, %d failures\n",
      cases,
      static_cast<long long>(shapes),
      fails);
  return fails == 0 ? 0 : 1;
}
