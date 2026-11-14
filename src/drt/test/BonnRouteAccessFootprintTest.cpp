// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// AccessFootprint (BrAccessFootprint.h): the global routing tiles recorded
// for shape grid and fast grid accesses equal a brute-force enumeration of
// the cells / vertices each call touches, on random rules, cell geometries,
// shapes and windows; accesses outside the allowed tile box are reported
// with their kind, layer, tiles, unit and step; a footprint records only
// its own thread.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <random>
#include <thread>
#include <utility>
#include <vector>

#include "BonnRouteRoutingSpaceTestUtil.h"
#include "br/BrAccessFootprint.h"
#include "br/BrFastGrid.h"
#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrTrackGraph.h"
#include "br/BrWireModel.h"

using namespace drt::br;
using namespace drt::br::test;

namespace {

// Different tile sizes and counts in x and y expose swapped axes.
constexpr Coord kExtent = 1500;
const ChunkAxis kAxisX{0, 200, 8};
const ChunkAxis kAxisY{0, 260, 6};
const TileBox kEverywhere{0, 0, kAxisX.count - 1, kAxisY.count - 1};

int failures = 0;

void expect(const bool ok, const char* what, const int c, const int k)
{
  if (!ok) {
    failures++;
    if (failures <= 20) {
      std::printf("case %d op %d: %s\n", c, k, what);
    }
  }
}

// Tile range of the cells [origin + c * cell, + cell) meeting [lo, hi).
std::pair<int, int> cellTiles(const ChunkAxis& axis,
                              const Coord origin,
                              const Coord cell,
                              const int cells,
                              const int64_t lo,
                              const int64_t hi)
{
  int first = std::numeric_limits<int>::max();
  int last = std::numeric_limits<int>::min();
  for (int c = 0; c < cells; c++) {
    const int64_t cell_lo = origin + int64_t{c} * cell;
    if (cell_lo < hi && cell_lo + cell > lo) {
      const int tile = axis.of(static_cast<Coord>(cell_lo));
      first = std::min(first, tile);
      last = std::max(last, tile);
    }
  }
  return {first, last};
}

// Tiles of the cells of layer geometry `g` meeting [xlo, xhi) x [ylo, yhi).
TileBox cellBox(const ShapeGridLayerGeometry& g,
                const int64_t xlo,
                const int64_t ylo,
                const int64_t xhi,
                const int64_t yhi)
{
  const bool h = g.horizontal;
  const auto [x0, x1] = cellTiles(kAxisX,
                                  h ? g.along_origin : g.cross_origin,
                                  h ? g.along_cell : g.cross_cell,
                                  h ? g.along_cells : g.cross_cells,
                                  xlo,
                                  xhi);
  const auto [y0, y1] = cellTiles(kAxisY,
                                  h ? g.cross_origin : g.along_origin,
                                  h ? g.cross_cell : g.along_cell,
                                  h ? g.cross_cells : g.along_cells,
                                  ylo,
                                  yhi);
  return {x0, y0, x1, y1};
}

void shapeGridCases(std::mt19937& rng, const int cases)
{
  for (int c = 0; c < cases; c++) {
    const int num_z = randomInt(rng, 2, 3);
    ShapeClassTable classes;
    classes.setWideRatioPercent(randomWideRatio(rng));
    classes.setLayers(randomRules(rng, num_z, rng() % 2));
    const std::vector<ShapeGridLayerGeometry> layers
        = gridGeometry(rng, classes, kExtent, c % 2 == 1);
    ShapeGrid grid;
    grid.init(&classes, layers, kAxisX, kAxisY);
    ShapeGridWriteCache cache;
    AccessFootprint footprint;
    const AccessFootprintScope scope(&footprint);
    std::vector<StoredShape> live;
    GridWindow window;
    for (int k = 0; k < 300; k++) {
      footprint.begin(k, kEverywhere);
      const int op = rng() % 3;
      if (op == 0 || live.empty()) {
        const StoredShape s = randomShape(
            rng, classes, randomInt(rng, 0, classes.numLayers() - 1), kExtent, 5);
        grid.insert(s.layer, s.rect, s.cls, s.type, s.owner, cache);
        live.push_back(s);
        const TileBox want = cellBox(
            layers[s.layer], s.rect.xlo, s.rect.ylo, s.rect.xhi, s.rect.yhi);
        expect(footprint.writes(0) == want, "insert: written tiles", c, k);
        expect(footprint.reads(0).empty(), "insert: no block read", c, k);
      } else if (op == 1) {
        const size_t at = rng() % live.size();
        const StoredShape s = live[at];
        live.erase(live.begin() + at);
        const bool removed
            = grid.remove(s.layer, s.rect, s.cls, s.type, s.owner, cache);
        const TileBox want = cellBox(
            layers[s.layer], s.rect.xlo, s.rect.ylo, s.rect.xhi, s.rect.yhi);
        expect(removed, "remove: shape present", c, k);
        expect(footprint.writes(0) == want, "remove: written tiles", c, k);
        // Cell shapes are tested for presence first; wide shapes are tested
        // inside the write.
        expect(classes.wide(s.cls) ? footprint.reads(0).empty()
                                   : footprint.reads(0) == want,
               "remove: tiles of the presence test",
               c,
               k);
      } else {
        const int layer = randomInt(rng, 0, classes.numLayers() - 1);
        const Coord x = randomInt(rng, 0, kExtent - 2);
        const Coord y = randomInt(rng, 0, kExtent - 2);
        const Rect rect{x,
                        y,
                        randomInt(rng, x, std::min(kExtent - 1, x + 400)),
                        randomInt(rng, y, std::min(kExtent - 1, y + 400))};
        grid.collect(layer, rect, window);
        // A window is closed: the cells of its upper corner are read too.
        const TileBox want = cellBox(layers[layer],
                                     rect.xlo,
                                     rect.ylo,
                                     int64_t{rect.xhi} + 1,
                                     int64_t{rect.yhi} + 1);
        expect(footprint.reads(0) == want, "collect: read tiles", c, k);
        expect(footprint.writes(0).empty(), "collect: no write", c, k);
        expect(footprint.wideReads().contains(want),
               "collect: wide lists of the window's tiles",
               c,
               k);
      }
      expect(footprint.unitViolations() == 0, "no violation", c, k);
    }
    // The halo derivation relies on the order of the three reaches for
    // every set of interned classes.
    for (int layer = 0; layer < classes.numLayers(); layer++) {
      expect(classes.narrowReach(layer) <= classes.narrowRuleReach(layer),
             "narrowReach <= narrowRuleReach",
             c,
             layer);
      expect(classes.narrowRuleReach(layer) <= classes.ruleReach(layer),
             "narrowRuleReach <= ruleReach",
             c,
             layer);
      expect(classes.classReach(layer) <= classes.ruleReach(layer),
             "classReach <= ruleReach",
             c,
             layer);
    }
  }
}

struct Graph
{
  TrackGraph graph;
  FastGrid fast;
};

void buildGraph(std::mt19937& rng, Graph& g, const int num_z)
{
  g.graph.layers.resize(num_z);
  const bool first_horizontal = rng() % 2;
  for (int z = 0; z < num_z; z++) {
    TrackLayer& l = g.graph.layers[z];
    l.horizontal = (z % 2 == 0) == first_horizontal;
    const Coord limit = l.horizontal ? kAxisY.size * kAxisY.count
                                     : kAxisX.size * kAxisX.count;
    for (Coord t = randomInt(rng, 3, 20); t < std::min(limit, kExtent) - 10;
         t += randomInt(rng, 24, 40)) {
      l.tracks.push_back(t);
    }
  }
  for (int z = 0; z < num_z; z++) {
    std::vector<Coord> positions;
    for (const int nz : {z - 1, z + 1}) {
      if (nz >= 0 && nz < num_z) {
        positions.insert(positions.end(),
                         g.graph.layers[nz].tracks.begin(),
                         g.graph.layers[nz].tracks.end());
      }
    }
    std::ranges::sort(positions);
    positions.erase(std::ranges::unique(positions).begin(), positions.end());
    g.graph.layers[z].positions = positions;
  }
  g.graph.via_cost.assign(num_z - 1, 10);
  g.graph.finalize([](Coord x) { return kAxisX.of(x); },
                   [](Coord y) { return kAxisY.of(y); });
  g.fast.init(&g.graph, 1);
}

// Tiles of the vertices (z, i, j0 .. j1).
TileBox vertexBox(const TrackGraph& graph,
                  const int z,
                  const int i,
                  const int j0,
                  const int j1)
{
  TileBox box;
  for (int j = j0; j <= j1; j++) {
    const auto [x, y] = graph.layers[z].point(i, j);
    const int tx = kAxisX.of(x);
    const int ty = kAxisY.of(y);
    box.add({tx, ty, tx, ty});
  }
  return box;
}

void fastGridCases(std::mt19937& rng, const int cases)
{
  for (int c = 0; c < cases; c++) {
    const int num_z = randomInt(rng, 2, 4);
    Graph g;
    buildGraph(rng, g, num_z);
    AccessFootprint footprint;
    const AccessFootprintScope scope(&footprint);
    FastGridReader reader;
    reader.reset(&g.fast, nullptr);
    for (int k = 0; k < 400; k++) {
      footprint.begin(k, kEverywhere);
      const int z = randomInt(rng, 0, num_z - 1);
      const TrackLayer& l = g.graph.layers[z];
      const int i = randomInt(rng, 0, l.numTracks() - 1);
      const int j0 = randomInt(rng, 0, l.numPositions() - 1);
      const int j1 = randomInt(rng, j0, l.numPositions() - 1);
      const int vias = z + 1 < num_z ? g.fast.numViaPositions(z) : 0;
      const int op = rng() % 8;
      if (op == 0) {
        g.fast.wiringWord(z, i, j0);
        expect(footprint.reads(0) == vertexBox(g.graph, z, i, j0, j0),
               "wiringWord",
               c,
               k);
      } else if (op == 1) {
        int lo;
        int hi;
        g.fast.wiringRun(z, i, j0, lo, hi);
        expect(footprint.reads(0) == vertexBox(g.graph, z, i, j0, j0),
               "wiringRun",
               c,
               k);
        // A run stays inside the tile that was recorded.
        expect(vertexBox(g.graph, z, i, lo, hi) == footprint.reads(0),
               "wiringRun: run inside the recorded tile",
               c,
               k);
      } else if (op == 2) {
        const FastGrid::RunsView view = g.fast.wiringRuns(z, i, j0);
        expect(footprint.reads(0) == vertexBox(g.graph, z, i, j0, j0),
               "wiringRuns",
               c,
               k);
        expect(vertexBox(g.graph, z, i, view.first, view.first + view.count - 1)
                   == footprint.reads(0),
               "wiringRuns: view inside the recorded tile",
               c,
               k);
      } else if (op == 3) {
        const uint64_t word = uint64_t{rng() % 8} << FastGrid::wiringShift(0, 0);
        g.fast.assignWiring(z, i, j0, j1, word, uint64_t{7});
        expect(footprint.writes(0) == vertexBox(g.graph, z, i, j0, j1),
               "assignWiring",
               c,
               k);
        expect(footprint.reads(0).empty(), "assignWiring: no read", c, k);
      } else if (op == 4) {
        std::vector<uint64_t> words(j1 - j0 + 1);
        for (uint64_t& w : words) {
          w = rng() % 8;
        }
        g.fast.writeWiring(z, i, j0, words.data(), j1 - j0 + 1, uint64_t{7});
        expect(footprint.writes(0) == vertexBox(g.graph, z, i, j0, j1),
               "writeWiring",
               c,
               k);
      } else if (op == 5 && vias > 0) {
        const int c0 = randomInt(rng, 0, vias - 1);
        const int c1 = randomInt(rng, c0, vias - 1);
        g.fast.assignVia(z, i, c0, c1, rng() % 8, uint64_t{7});
        TileBox want;
        for (int v = c0; v <= c1; v++) {
          const int j = g.fast.viaPosition(z, v);
          want.add(vertexBox(g.graph, z, i, j, j));
        }
        expect(footprint.writes(0) == want, "assignVia", c, k);
      } else if (op == 6 && vias > 0) {
        const int v = randomInt(rng, 0, vias - 1);
        const int j = g.fast.viaPosition(z, v);
        int lo;
        int hi;
        if (rng() % 2 == 0) {
          g.fast.viaWord(z, i, v);
        } else {
          g.fast.viaRun(z, i, v, lo, hi);
        }
        expect(footprint.reads(0) == vertexBox(g.graph, z, i, j, j),
               "viaWord / viaRun",
               c,
               k);
      } else {
        // The search-side reader: a level query records the tiles of the
        // vertices whose words it fetched (the vertex, and for an edge its
        // other endpoint).
        const bool edge = j0 + 1 < l.numPositions() && rng() % 2 == 0;
        reader.reset(&g.fast, nullptr);
        reader.level(z, i, j0, edge ? kWireEdge : kWire, 0);
        TileBox want = vertexBox(g.graph, z, i, j0, j0);
        if (edge) {
          want.add(vertexBox(g.graph, z, i, j0 + 1, j0 + 1));
        }
        expect(want.contains(footprint.reads(0))
                   && footprint.reads(0).contains(
                       vertexBox(g.graph, z, i, j0, j0)),
               "reader level",
               c,
               k);
      }
      expect(footprint.unitViolations() == 0, "no violation", c, k);
    }
  }
}

// A commit journal records writes even when its reads hit warm caches.
// A nested search still installs a complete read footprint, then restores the
// journal's mode. History writes participate in the same invalidation box.
void writeOnlyCase(std::mt19937& rng)
{
  Graph g;
  buildGraph(rng, g, 2);
  FastGridReader reader;
  reader.reset(&g.fast, nullptr);
  const size_t v = g.graph.vertex(0, 0, 0);
  const TileBox tile = vertexBox(g.graph, 0, 0, 0, 0);
  const int via_position = g.fast.viaPosition(0, 0);
  TileBox want = tile;
  want.add(vertexBox(g.graph, 0, 0, via_position, via_position));
  reader.level(0, 0, 0, kWire, 0);
  AccessFootprint journal(AccessFootprint::Mode::kWrites);
  journal.begin(0, kEverywhere);
  {
    const AccessFootprintScope scope(&journal);
    expect(AccessFootprint::currentForRead() == nullptr,
           "journal skips read instrumentation",
           0,
           0);
    g.fast.setLevel(v, kWire, 0, RipupLevel{1});
    g.fast.assignVia(0, 0, 0, 0, 0, uint64_t{7});
    g.fast.bumpHistory(v);
    reader.level(0, 0, 0, kWire, 0);
    g.fast.history(v);
    g.fast.viaWord(0, 0, 0);
    expect(journal.writes(0) == want && journal.reads(0).empty(),
           "journal retains wiring via and history writes",
           0,
           1);
    AccessFootprint search;
    search.begin(1, kEverywhere);
    {
      const AccessFootprintScope nested(&search);
      reader.level(0, 0, 0, kWire, 0);
      g.fast.history(v);
    }
    expect(
        search.reads(0) == tile && AccessFootprint::currentForRead() == nullptr,
        "nested search records warm reads and restores journal mode",
        0,
        2);
    recordAccess(AccessKind::kWideShapeRead, 0, kEverywhere);
    expect(
        journal.wideReads().empty(), "journal ignores wide-list reads", 0, 3);
    recordAccess(AccessKind::kWideShapeWrite, 0, kEverywhere);
    expect(journal.writes(0) == kEverywhere,
           "journal retains wide-list writes",
           0,
           4);
  }
}

void violationCases(std::mt19937& rng)
{
  ShapeClassTable classes;
  classes.setWideRatioPercent(0);
  classes.setLayers(randomRules(rng, 2, true));
  const std::vector<ShapeGridLayerGeometry> layers
      = gridGeometry(rng, classes, kExtent, false);
  ShapeGrid grid;
  grid.init(&classes, layers, kAxisX, kAxisY);
  ShapeGridWriteCache cache;
  AccessFootprint footprint;
  GridWindow window;
  const ClassId cls
      = classes.classify(0, {0, 0, 10, 10}, ShapeRole::kViaPad, 0);
  // Tile (2, 3) is [400, 600) x [780, 1040).
  const TileBox tile{2, 3, 2, 3};
  const Rect inside{450, 800, 470, 820};
  // Into the next tile by more than a cell, so a cell anchored there is met.
  const Rect right{590, 800, 660, 820};    // reaches into tile (3, 3)
  const Rect above{450, 1030, 470, 1100};  // reaches into tile (2, 4)
  {
    const AccessFootprintScope scope(&footprint);
    footprint.begin(7, tile);
    footprint.setStep(2);
    grid.insert(0, inside, cls, ShapeType::kWire, 1, cache);
    grid.collect(0, inside, window);
    expect(footprint.unitViolations() == 0, "inside: no violation", 0, 0);
    expect(footprint.writes(2) == tile && footprint.reads(2) == tile,
           "inside: step boxes",
           0,
           0);
    expect(footprint.accessed() == tile, "inside: accessed", 0, 0);
    grid.insert(0, right, cls, ShapeType::kWire, 1, cache);
    expect(footprint.unitViolations() == 1, "right: one violation", 0, 1);
    footprint.setStep(4);
    grid.collect(0, above, window);
    expect(footprint.unitViolations() == 2, "above: second violation", 0, 2);
    expect(footprint.violations().size() == 2, "two violations kept", 0, 2);
    if (footprint.violations().size() == 2) {
      const AccessFootprint::Violation& a = footprint.violations()[0];
      const AccessFootprint::Violation& b = footprint.violations()[1];
      expect(a.kind == AccessKind::kShapeGridWrite && a.layer == 0
                 && a.tiles == TileBox{2, 3, 3, 3} && a.allowed == tile
                 && a.unit == 7 && a.step == 2,
             "right: violation details",
             0,
             1);
      expect(b.kind == AccessKind::kShapeGridRead && b.layer == 0
                 && b.tiles == TileBox{2, 3, 2, 4} && b.unit == 7
                 && b.step == 4,
             "above: violation details",
             0,
             2);
    }
    expect(tile.reachOf(footprint.accessed()) == 1, "reach 1 tile", 0, 2);
    // A new unit starts with empty boxes; the totals stay.
    footprint.begin(8, kEverywhere);
    expect(footprint.accessed().empty() && footprint.unitViolations() == 0
               && footprint.violationCount() == 2 && footprint.step() == 0,
           "begin resets the unit",
           0,
           3);
    // At most kMaxViolations are kept; all are counted.
    footprint.begin(9, tile);
    for (int k = 0; k < 40; k++) {
      grid.collect(0, right, window);
    }
    expect(footprint.unitViolations() == 40
               && footprint.violationCount() == 42
               && footprint.violations().size()
                      == AccessFootprint::kMaxViolations,
           "violation cap",
           0,
           4);
    // Nested scopes restore the outer footprint.
    AccessFootprint inner;
    {
      const AccessFootprintScope nested(&inner);
      inner.begin(1, kEverywhere);
      grid.collect(0, inside, window);
    }
    expect(inner.records() > 0 && AccessFootprint::current() == &footprint,
           "nested scope",
           0,
           5);
    // Another thread's accesses are its own.
    const int64_t before = footprint.records();
    std::thread other([&]() {
      GridWindow w;
      grid.collect(0, inside, w);
    });
    other.join();
    expect(footprint.records() == before, "thread local", 0, 6);
  }
  expect(AccessFootprint::current() == nullptr, "scope uninstalls", 0, 7);
  const int64_t before = footprint.records();
  grid.collect(0, inside, window);
  expect(footprint.records() == before, "no footprint: no record", 0, 8);
}

// The cut-to-projection rule of a via layer is a rule between narrow
// classes: it enters narrowReach() once both classes exist, and
// narrowRuleReach() from the start.
void projectionReachCase()
{
  std::vector<LayerRules> rules(5);
  for (int l = 0; l < 5; l++) {
    rules[l].via_layer = (l % 2) == 1;
    rules[l].horizontal = (l / 2) % 2 == 0;
    rules[l].spacing = SpacingTable::constant(6);
  }
  rules[3].projection_spacing = 14;
  ShapeClassTable classes;
  classes.setLayers(rules);
  expect(classes.narrowRuleReach(3) == 14 && classes.ruleReach(3) == 14,
         "projection spacing in the rule reaches",
         0,
         0);
  expect(classes.narrowRuleReach(1) == 6, "plain via layer", 0, 1);
  classes.classify(3, {0, 0, 6, 6}, ShapeRole::kCut, 0);
  expect(classes.narrowReach(3) == 6, "cuts alone", 0, 2);
  classes.classify(3, {0, 0, 6, 6}, ShapeRole::kCutProjection, 0);
  expect(classes.narrowReach(3) == 14, "cut and projection", 0, 3);
  expect(classes.narrowReach(3) <= classes.narrowRuleReach(3),
         "narrowReach <= narrowRuleReach",
         0,
         4);
}

}  // namespace

int main()
{
  std::mt19937 rng(20261113);
  projectionReachCase();
  shapeGridCases(rng, 60);
  fastGridCases(rng, 60);
  violationCases(rng);
  writeOnlyCase(rng);
  std::printf("access footprint: %d failures\n", failures);
  return failures != 0;
}
