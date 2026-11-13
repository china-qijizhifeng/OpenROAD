// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Sec. 5.1 partitions: threads change the routing space only inside their
// regions of global routing tiles. Shape grid and fast grid are chunked by
// tile, so threads update them in place without touching each other's
// chunks. This test runs random insert / remove / query sequences in
// disjoint tile regions on parallel threads and checks that the final
// shape grid content, the fast grid words and every query answer equal a
// serial run of the same sequences (also run it with -fsanitize=thread).
// Every region thread records its access footprint (BrAccessFootprint.h):
// all block and tile accesses stay inside the region's own tiles.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <map>
#include <random>
#include <thread>
#include <tuple>
#include <vector>

#include "BonnRouteRoutingSpaceTestUtil.h"
#include "br/BrAccessFootprint.h"
#include "br/BrDistanceRules.h"
#include "br/BrFastGrid.h"
#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrTrackGraph.h"
#include "br/BrWireModel.h"

using namespace drt::br;
using namespace drt::br::test;

namespace {

constexpr Coord kTile = 240;
constexpr int kTiles = 9;
constexpr Coord kExtent = kTile * kTiles;

struct Op
{
  bool insert;
  StoredShape shape;
  Point query;
};

struct Space
{
  ShapeClassTable* classes;
  TrackGraph graph;
  std::vector<WireType> types;
  ShapeGrid grid;
  DistanceRuleChecker checker;
  FastGrid fast;
  FastGridBuilder builder;
};

void build(Space& s,
           ShapeClassTable& classes,
           const std::vector<WireTypeGeometry>& geometry,
           const int num_z,
           const std::vector<StoredShape>& initial)
{
  s.classes = &classes;
  s.graph.layers.resize(num_z);
  for (int z = 0; z < num_z; z++) {
    TrackLayer& l = s.graph.layers[z];
    l.horizontal = classes.layer(wiringGridLayer(z)).horizontal;
    for (Coord x = 10 + 3 * z; x < kExtent - 10; x += 28 + 2 * z) {
      l.tracks.push_back(x);
    }
  }
  for (int z = 0; z < num_z; z++) {
    std::vector<Coord> pos;
    for (const int nz : {z - 1, z + 1}) {
      if (nz >= 0 && nz < num_z) {
        pos.insert(pos.end(),
                   s.graph.layers[nz].tracks.begin(),
                   s.graph.layers[nz].tracks.end());
      }
    }
    std::ranges::sort(pos);
    pos.erase(std::ranges::unique(pos).begin(), pos.end());
    s.graph.layers[z].positions = pos;
  }
  s.graph.via_cost.assign(num_z - 1, 10);
  const ChunkAxis axis{0, kTile, kTiles};
  s.graph.finalize([&](Coord x) { return axis.of(x); },
                   [&](Coord y) { return axis.of(y); });
  std::vector<ShapeGridLayerGeometry> layers;
  for (int l = 0; l < classes.numLayers(); l++) {
    ShapeGridLayerGeometry g;
    g.horizontal = classes.layer(l).horizontal;
    g.along_cell = g.cross_cell
        = std::max<Coord>(1, maximalCellSize(classes, l));
    g.along_cells = g.cross_cells = kExtent / g.along_cell + 1;
    layers.push_back(g);
  }
  s.grid.init(&classes, layers, axis, axis);
  for (const WireTypeGeometry& g : geometry) {
    s.types.push_back(makeWireType(g, classes));
  }
  s.checker.init(&s.grid);
  s.fast.init(&s.graph, 1);
  s.builder.init(&s.graph, &s.types, 1, &s.checker, &s.fast);
  ShapeGridWriteCache cache;
  for (const StoredShape& sh : initial) {
    s.grid.insert(sh.layer, sh.rect, sh.cls, sh.type, sh.owner, cache);
  }
  // The initial build in parallel, as BonnRouteDR::rebuildFastGrid.
  const auto jobs = s.builder.rebuildJobs();
  std::vector<std::vector<FastGridBuilder::JogBits>> bits(jobs.size());
  for (int phase = 0; phase < 3; phase++) {
    std::vector<std::thread> threads;
    for (int t = 0; t < 4; t++) {
      threads.emplace_back([&, t]() {
        for (size_t k = t; k < jobs.size(); k += 4) {
          if (phase == 0) {
            s.builder.rebuildVertices(jobs[k][0], jobs[k][1], jobs[k][2]);
          } else if (phase == 1) {
            s.builder.computeJogBits(
                jobs[k][0], jobs[k][1], jobs[k][2], bits[k]);
          } else {
            s.builder.writeJogBits(bits[k]);
          }
        }
      });
    }
    for (std::thread& t : threads) {
      t.join();
    }
  }
}

// Runs one region's operations; returns the query answers.
std::vector<std::vector<RipupLevel>> run(Space& s,
                                         const std::vector<Op>& ops,
                                         ShapeGridWriteCache& cache)
{
  std::vector<std::vector<RipupLevel>> answers;
  std::vector<LocatedModel> models;
  const int num_z = s.graph.numLayers();
  for (int z = 0; z < num_z; z++) {
    appendWireTypeModels(s.types[0], z, num_z, 4 * z, models);
  }
  for (const Op& op : ops) {
    const StoredShape& sh = op.shape;
    if (op.insert) {
      s.grid.insert(sh.layer, sh.rect, sh.cls, sh.type, sh.owner, cache);
      s.builder.shapeInserted(sh.layer,
                              PlacedShape{sh.rect, sh.cls, 3},
                              s.grid.level(sh.type, sh.owner));
    } else {
      s.grid.remove(sh.layer, sh.rect, sh.cls, sh.type, sh.owner, cache);
      s.builder.shapeRemoved(sh.layer,
                             PlacedShape{sh.rect, sh.cls, 3},
                             s.grid.level(sh.type, sh.owner));
    }
    CheckOptions opts;
    opts.report_nets = true;
    opts.groups = 4 * num_z;
    CheckResult res;
    s.checker.check(op.query, models, opts, res);
    answers.push_back(res.levels);
    for (const BlockingNet& n : res.nets) {
      answers.back().push_back(static_cast<RipupLevel>(n.net));
      answers.back().push_back(n.level);
    }
    answers.back().push_back(static_cast<RipupLevel>(res.x_hi - res.x_lo));
    // The fast grid as the path search reads it, near the query.
    const TrackLayer& l = s.graph.layers[0];
    const int i
        = static_cast<int>(std::ranges::lower_bound(
                               l.tracks, l.horizontal ? op.query.y : op.query.x)
                           - l.tracks.begin());
    const int j = static_cast<int>(
        std::ranges::lower_bound(l.positions,
                                 l.horizontal ? op.query.x : op.query.y)
        - l.positions.begin());
    if (i < l.numTracks() && j + 1 < l.numPositions()) {
      answers.back().push_back(s.fast.level(0, i, j, kWire, 0));
      answers.back().push_back(s.fast.level(0, i, j, kWireEdge, 0));
    }
  }
  return answers;
}

using Content = std::map<std::tuple<int, int, int>, std::vector<CellShape>>;

Content content(const Space& s)
{
  Content c;
  for (int layer = 0; layer < s.grid.numLayers(); layer++) {
    s.grid.forEachInterval(layer, [&](const ShapeGrid::IntervalView& v) {
      for (int cell = v.begin; cell <= v.end; cell++) {
        auto& list = c[{layer, v.row, cell}];
        if (v.owner == ShapeGrid::kSharedCell) {
          for (const auto& [owner, config] : v.parts) {
            const auto& sh = s.grid.configurations()[config].shapes;
            list.insert(list.end(), sh.begin(), sh.end());
          }
        } else {
          const auto& sh = s.grid.configurations()[v.config].shapes;
          list.insert(list.end(), sh.begin(), sh.end());
        }
        std::ranges::sort(list);
      }
    });
  }
  return c;
}

}  // namespace

int main()
{
  std::mt19937 rng(97);
  int fails = 0;
  for (int c = 0; c < 6; c++) {
    const int num_z = randomInt(rng, 2, 3);
    // One class table: classes are interned before the threads start.
    ShapeClassTable classes;
    classes.setWideRatioPercent(randomWideRatio(rng));
    classes.setLayers(randomRules(rng, num_z, rng() % 2));
    std::vector<WireTypeGeometry> geometry(1);
    for (int z = 0; z < num_z; z++) {
      geometry[0].width.push_back(10);
    }
    for (int z = 0; z + 1 < num_z; z++) {
      geometry[0].via.push_back(
          {true, {-7, -5, 7, 5}, {-3, -3, 3, 3}, {-5, -7, 5, 7}});
    }
    // Fixed shapes all over the layout before the (parallel) initial build.
    std::vector<StoredShape> initial;
    for (int k = 0; k < 400; k++) {
      initial.push_back(randomShape(rng,
                                    classes,
                                    randomInt(rng, 0, classes.numLayers() - 1),
                                    kExtent,
                                    5));
    }
    Space parallel;
    Space serial;
    build(parallel, classes, geometry, num_z, initial);
    build(serial, classes, geometry, num_z, initial);
    // Regions: 2 x 2 tiles at (3a, 3b), one free tile between regions.
    // Shapes keep the router's halo from the region border: updates and
    // queries reach at most 2 (reach + model + cell) + pitch beyond them.
    std::vector<std::vector<Op>> ops;
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        const Coord x0 = 3 * a * kTile + 150;
        const Coord y0 = 3 * b * kTile + 150;
        const Coord size = 2 * kTile - 300;
        std::vector<Op> list;
        std::vector<StoredShape> live;
        for (int k = 0; k < 150; k++) {
          Op op;
          op.query
              = {x0 + randomInt(rng, 0, size), y0 + randomInt(rng, 0, size)};
          if (live.empty() || rng() % 3 != 0) {
            StoredShape sh
                = randomShape(rng,
                              classes,
                              randomInt(rng, 0, classes.numLayers() - 1),
                              size,
                              5);
            sh.rect = sh.rect.translated({x0, y0});
            op.insert = true;
            op.shape = sh;
            live.push_back(sh);
          } else {
            const size_t at = rng() % live.size();
            op.insert = false;
            op.shape = live[at];
            live.erase(live.begin() + at);
          }
          list.push_back(op);
        }
        ops.push_back(list);
      }
    }
    std::vector<std::vector<std::vector<RipupLevel>>> answers(ops.size());
    std::vector<ShapeGridWriteCache> caches(ops.size());
    std::vector<AccessFootprint> footprints(ops.size());
    std::vector<std::thread> threads;
    for (size_t r = 0; r < ops.size(); r++) {
      threads.emplace_back([&, r]() {
        // Region r = (a, b): tiles [3a, 3a + 1] x [3b, 3b + 1].
        const int a = static_cast<int>(r) / 3;
        const int b = static_cast<int>(r) % 3;
        const AccessFootprintScope scope(&footprints[r]);
        footprints[r].begin(static_cast<int32_t>(r),
                            {3 * a, 3 * b, 3 * a + 1, 3 * b + 1});
        answers[r] = run(parallel, ops[r], caches[r]);
      });
    }
    for (std::thread& t : threads) {
      t.join();
    }
    for (size_t r = 0; r < ops.size(); r++) {
      const AccessFootprint& f = footprints[r];
      if (f.records() == 0 || f.violationCount() != 0) {
        std::printf(
            "case %d region %zu: %lld accesses recorded, %lld outside the "
            "region\n",
            c,
            r,
            static_cast<long long>(f.records()),
            static_cast<long long>(f.violationCount()));
        for (const AccessFootprint::Violation& v : f.violations()) {
          std::printf("  %s layer %d tiles [%d, %d] - [%d, %d]\n",
                      accessKindName(v.kind),
                      v.layer,
                      v.tiles.xlo,
                      v.tiles.ylo,
                      v.tiles.xhi,
                      v.tiles.yhi);
        }
        fails++;
      }
    }
    ShapeGridWriteCache serial_cache;
    for (size_t r = 0; r < ops.size(); r++) {
      if (run(serial, ops[r], serial_cache) != answers[r]) {
        std::printf("case %d region %zu: answers differ\n", c, r);
        fails++;
      }
    }
    if (content(parallel) != content(serial)) {
      std::printf("case %d: shape grid content differs\n", c);
      fails++;
    }
    for (int z = 0; z < num_z; z++) {
      const TrackLayer& l = parallel.graph.layers[z];
      for (int i = 0; i < l.numTracks(); i++) {
        for (int j = 0; j < l.numPositions(); j++) {
          if (parallel.fast.wiringWord(z, i, j)
              != serial.fast.wiringWord(z, i, j)) {
            fails++;
          }
        }
        for (int v = 0; z + 1 < num_z && v < parallel.fast.numViaPositions(z);
             v++) {
          if (parallel.fast.viaWord(z, i, v) != serial.fast.viaWord(z, i, v)) {
            fails++;
          }
        }
      }
    }
  }
  std::printf("partitions: %d failures\n", fails);
  return fails != 0;
}
