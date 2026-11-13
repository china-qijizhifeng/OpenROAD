// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Dense fast-grid encoding (fast_grid_dense) against the run encoding and
// the brute-force reference: two fast grids over one shape grid, one in
// runs and one with a random set of dense layers (switched on and off in
// the middle of the sequence), receive the same random insert / remove /
// group-remove / exact-recompute sequence. After every step every wiring
// and via word, the runs reported by wiringRun / viaRun, and the
// statistics are compared; at the end the readers (with and without an
// own-shape overlay) are compared along random access walks and the dense
// grid's vertex levels against the brute-force reference.

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <random>
#include <thread>
#include <vector>

#include "BonnRouteRoutingSpaceTestUtil.h"
#include "br/BrDistanceRules.h"
#include "br/BrFastGrid.h"
#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrTrackGraph.h"
#include "br/BrWireModel.h"

using namespace drt::br;
using namespace drt::br::test;

namespace {

struct Totals
{
  int64_t cases = 0;
  int64_t steps = 0;
  int64_t word_checks = 0;
  int64_t run_checks = 0;
  int64_t reader_checks = 0;
  int64_t brute_checks = 0;
  int64_t dense_layers = 0;
  int64_t switches = 0;
  int64_t overlays = 0;
  int64_t nonzero = 0;
  int64_t large = 0;
  int64_t fails = 0;
};

void fail(Totals& t, const char* what, const int c, const int step)
{
  if (t.fails++ < 12) {
    std::printf("case %d step %d: %s\n", c, step, what);
  }
}

// Every word of both grids, the runs of the dense one (a maximal run of
// equal words inside the tile chunk containing the position) and the
// statistics.
void compareGrids(const TrackGraph& graph,
                  const FastGrid& runs,
                  const FastGrid& dense,
                  Totals& t,
                  const int c,
                  const int step)
{
  const int num_z = graph.numLayers();
  for (int z = 0; z < num_z; z++) {
    const TrackLayer& l = graph.layers[z];
    for (int i = 0; i < l.numTracks(); i++) {
      for (int j = 0; j < l.numPositions(); j++) {
        const uint64_t a = runs.wiringWord(z, i, j);
        const uint64_t b = dense.wiringWord(z, i, j);
        t.word_checks++;
        t.nonzero += a != 0;
        if (a != b) {
          fail(t, "wiring word differs", c, step);
        }
        int alo, ahi, blo, bhi;
        const uint64_t ra = runs.wiringRun(z, i, j, alo, ahi);
        const uint64_t rb = dense.wiringRun(z, i, j, blo, bhi);
        t.run_checks++;
        // Both encodings keep maximal runs per tile chunk.
        if (ra != a || rb != b || alo != blo || ahi != bhi || blo > j
            || bhi < j) {
          fail(t, "wiring run differs", c, step);
        }
        const FastGrid::RunsView va = runs.wiringRuns(z, i, j);
        const FastGrid::RunsView vb = dense.wiringRuns(z, i, j);
        if (va.first != vb.first || va.count != vb.count) {
          fail(t, "runs view chunk differs", c, step);
        }
        if (dense.dense(z) && (vb.dense == nullptr || vb.dense[j - vb.first] != b)) {
          fail(t, "dense runs view word differs", c, step);
        }
      }
      for (int cc = 0; z + 1 < num_z && cc < runs.numViaPositions(z); cc++) {
        const uint64_t a = runs.viaWord(z, i, cc);
        const uint64_t b = dense.viaWord(z, i, cc);
        t.word_checks++;
        if (a != b) {
          fail(t, "via word differs", c, step);
        }
        int alo, ahi, blo, bhi;
        runs.viaRun(z, i, cc, alo, ahi);
        dense.viaRun(z, i, cc, blo, bhi);
        t.run_checks++;
        if (alo != blo || ahi != bhi) {
          fail(t, "via run differs", c, step);
        }
      }
    }
  }
  const FastGrid::Statistics sa = runs.statistics();
  const FastGrid::Statistics sb = dense.statistics();
  if (sa.wiring_intervals != sb.wiring_intervals
      || sa.via_intervals != sb.via_intervals
      || sa.marked_edges != sb.marked_edges || sa.vertices != sb.vertices
      || sa.via_vertices != sb.via_vertices) {
    fail(t, "statistics differ", c, step);
  }
}

// Brute-force level of one model (rect relative to l), the shapes of
// `ignored` treated as removed.
RipupLevel brute(const ShapeClassTable& classes,
                 const ShapeGrid& grid,
                 const std::vector<StoredShape>& shapes,
                 const LocatedModel& m,
                 const Point l,
                 const int32_t ignored)
{
  std::vector<LocatedModel> one{m};
  one[0].group = 0;
  return bruteLevels(classes, grid, shapes, one, 1, l, ignored)[0];
}

void runCase(std::mt19937& rng, const int c, Totals& t)
{
  t.cases++;
  const int num_z = randomInt(rng, 1, 3);
  const bool first_horizontal = rng() % 2;
  ShapeClassTable classes;
  classes.setWideRatioPercent(randomWideRatio(rng));
  classes.setLayers(randomRules(rng, num_z, first_horizontal));
  // Large cases have more than 64 tracks per layer, so tracks share
  // reader cursor slots (FastGridReader::slot).
  const bool large = rng() % 4 == 0;
  const Coord extent = large ? randomInt(rng, 900, 1300) : randomInt(rng, 150, 320);
  t.large += large;
  TrackGraph graph;
  graph.layers.resize(num_z);
  std::vector<std::vector<Coord>> tracks(num_z);
  for (int z = 0; z < num_z; z++) {
    const int pitch = large ? randomInt(rng, 7, 12) : randomInt(rng, 16, 40);
    for (Coord x = randomInt(rng, 4, pitch); x < extent - 4;
         x += pitch + (rng() % 6 == 0 ? randomInt(rng, 1, 30) : 0)) {
      tracks[z].push_back(x);
    }
  }
  for (int z = 0; z < num_z; z++) {
    TrackLayer& l = graph.layers[z];
    l.horizontal = classes.layer(wiringGridLayer(z)).horizontal;
    l.tracks = tracks[z];
    std::vector<Coord> pos;
    for (const int nz : {z - 1, z + 1}) {
      if (nz >= 0 && nz < num_z) {
        pos.insert(pos.end(), tracks[nz].begin(), tracks[nz].end());
      }
    }
    if (pos.empty()) {
      for (Coord x = 8; x < extent - 8; x += randomInt(rng, 10, 30)) {
        pos.push_back(x);
      }
    }
    std::ranges::sort(pos);
    pos.erase(std::ranges::unique(pos).begin(), pos.end());
    l.positions = pos;
    l.jogs = rng() % 4 != 0;
  }
  for (int z = 0; z + 1 < num_z; z++) {
    graph.via_cost.push_back(10);
  }
  // Small chunks: many tile borders inside a track.
  ChunkAxis cx{0, randomInt(rng, 30, 300), 0};
  cx.count = extent / cx.size + 1;
  ChunkAxis cy{0, randomInt(rng, 30, 300), 0};
  cy.count = extent / cy.size + 1;
  graph.finalize([&](Coord x) { return cx.of(x); },
                 [&](Coord y) { return cy.of(y); });
  ShapeGrid grid;
  grid.init(&classes, gridGeometry(rng, classes, extent, rng() % 2), cx, cy);
  const int nets = randomInt(rng, 1, 4);
  for (int n = 0; n < nets; n++) {
    grid.setNetRipupLevel(n, static_cast<RipupLevel>(randomInt(rng, 1, 6)));
  }
  std::vector<WireType> types;
  const int fast
      = rng() % 4 == 0 ? FastGrid::kMaxWireTypes : randomInt(rng, 1, 2);
  for (int k = 0; k < fast + 1; k++) {
    WireTypeGeometry g;
    g.ndr = k == 0 ? 0 : randomInt(rng, 0, 2);
    for (int z = 0; z < num_z; z++) {
      g.width.push_back(randomInt(rng, 6, 16));
    }
    for (int z = 0; z + 1 < num_z; z++) {
      WireTypeGeometry::Via v;
      v.valid = rng() % 6 != 0;
      const Coord cut = randomInt(rng, 2, 5);
      v.cut = {-cut, -cut, cut, cut};
      const Coord b = cut + randomInt(rng, 0, 4);
      const Coord bx = b + randomInt(rng, 0, 5);
      v.bottom = {-bx, -b, bx, b};
      v.top = {-b, -bx, b, bx};
      g.via.push_back(v);
    }
    types.push_back(makeWireType(g, classes));
  }
  DistanceRuleChecker checker(&grid);
  FastGrid runs;
  runs.init(&graph, fast);
  FastGridBuilder runs_builder;
  runs_builder.init(&graph, &types, fast, &checker, &runs);
  FastGrid dense;
  dense.init(&graph, fast);
  // Random dense layers before the initial build (as BonnRouteDR does);
  // more are switched later.
  for (int z = 0; z < num_z; z++) {
    if (rng() % 4 != 0) {
      dense.setDense(z, true);
    }
  }
  FastGridBuilder dense_builder;
  dense_builder.init(&graph, &types, fast, &checker, &dense);
  // Half of the cases recompute removal groups of the dense grid on
  // several threads (BonnRouteDR's flush outside parallel regions).
  if (rng() % 2 == 0) {
    dense_builder.setParallelFor(
        [](int count, const std::function<void(int)>& body) {
          std::atomic<int> next{0};
          std::vector<std::thread> threads;
          for (int k = 0; k < 4; k++) {
            threads.emplace_back([&] {
              for (int q = next++; q < count; q = next++) {
                body(q);
              }
            });
          }
          for (std::thread& th : threads) {
            th.join();
          }
        },
        1);
  }
  ShapeGridWriteCache cache;
  std::vector<StoredShape> live;
  const int before = randomInt(rng, 0, 20);
  for (int k = 0; k < before; k++) {
    const StoredShape s = randomShape(
        rng, classes, randomInt(rng, 0, classes.numLayers() - 1), extent, nets);
    grid.insert(s.layer, s.rect, s.cls, s.type, s.owner, cache);
    live.push_back(s);
  }
  for (FastGridBuilder* b : {&runs_builder, &dense_builder}) {
    for (int z = 0; z < num_z; z++) {
      b->rebuildVertices(z, 0, graph.layers[z].numTracks() - 1);
    }
    for (int z = 0; z < num_z; z++) {
      b->rebuildJogBits(z, 0, graph.layers[z].numTracks() - 1);
    }
  }
  compareGrids(graph, runs, dense, t, c, -1);
  const int steps = randomInt(rng, 0, 30);
  for (int k = 0; k < steps; k++) {
    t.steps++;
    if (rng() % 6 == 0) {
      // Switch a layer between the encodings; every word is kept.
      const int z = randomInt(rng, 0, num_z - 1);
      dense.setDense(z, !dense.dense(z));
      t.switches++;
    }
    if (live.empty() || rng() % 3 != 0) {
      StoredShape s = randomShape(rng,
                                  classes,
                                  randomInt(rng, 0, classes.numLayers() - 1),
                                  extent,
                                  nets);
      const int z = s.layer / 2;
      if (s.layer % 2 == 0 && rng() % 2 && graph.layers[z].numTracks() > 0) {
        const TrackLayer& l = graph.layers[z];
        const WireModel& m = types[0].preferred[z];
        const Coord tr = l.tracks[rng() % l.numTracks()];
        const Coord a = l.positions[rng() % l.numPositions()];
        const Coord b = a + randomInt(rng, 0, 80);
        s.rect = l.horizontal ? Rect{a + m.rect.xlo,
                                     tr + m.rect.ylo,
                                     b + m.rect.xhi,
                                     tr + m.rect.yhi}
                              : Rect{tr + m.rect.xlo,
                                     a + m.rect.ylo,
                                     tr + m.rect.xhi,
                                     b + m.rect.yhi};
        s.cls = m.cls;
        s.type = ShapeType::kWire;
        s.owner = randomInt(rng, 0, nets - 1);
      }
      grid.insert(s.layer, s.rect, s.cls, s.type, s.owner, cache);
      if (rng() % 8 == 0) {
        runs_builder.shapeChanged(s.layer, s.rect);
        dense_builder.shapeChanged(s.layer, s.rect);
      } else {
        const PlacedShape placed{s.rect, s.cls, 3};
        const RipupLevel level = grid.level(s.type, s.owner);
        runs_builder.shapeInserted(s.layer, placed, level);
        dense_builder.shapeInserted(s.layer, placed, level);
      }
      live.push_back(s);
    } else if (rng() % 3 == 0) {
      FastGridBuilder::PendingRemovals pending_runs;
      FastGridBuilder::PendingRemovals pending_dense;
      for (int r = randomInt(rng, 1, rng() % 2 ? 4 : 12);
           r > 0 && !live.empty();
           r--) {
        const size_t at = rng() % live.size();
        const StoredShape s = live[at];
        live.erase(live.begin() + at);
        grid.remove(s.layer, s.rect, s.cls, s.type, s.owner, cache);
        const PlacedShape placed{s.rect, s.cls, 3};
        const RipupLevel level = grid.level(s.type, s.owner);
        runs_builder.shapeRemoved(s.layer, placed, level, pending_runs);
        dense_builder.shapeRemoved(s.layer, placed, level, pending_dense);
      }
      runs_builder.flush(pending_runs);
      dense_builder.flush(pending_dense);
    } else {
      const size_t at = rng() % live.size();
      const StoredShape s = live[at];
      live.erase(live.begin() + at);
      grid.remove(s.layer, s.rect, s.cls, s.type, s.owner, cache);
      const PlacedShape placed{s.rect, s.cls, 3};
      const RipupLevel level = grid.level(s.type, s.owner);
      runs_builder.shapeRemoved(s.layer, placed, level);
      dense_builder.shapeRemoved(s.layer, placed, level);
    }
    compareGrids(graph, runs, dense, t, c, k);
  }
  for (int z = 0; z < num_z; z++) {
    t.dense_layers += dense.dense(z);
  }
  // Readers along random walks (track scans in both directions, jogs,
  // vias, jumps), with and without the overlay of one net.
  const int32_t own = randomInt(rng, 0, nets - 1);
  std::vector<FastGridBuilder::OwnShape> own_shapes;
  for (const StoredShape& s : live) {
    if (s.owner == own && isOwned(s.type)) {
      own_shapes.push_back({s.layer,
                            PlacedShape{s.rect, s.cls, 3},
                            grid.level(s.type, s.owner)});
    }
  }
  t.overlays += !own_shapes.empty();
  FastGridOverlay overlay_runs;
  FastGridOverlay overlay_dense;
  runs_builder.buildOverlay(own, own_shapes, overlay_runs);
  dense_builder.buildOverlay(own, own_shapes, overlay_dense);
  if (overlay_runs.wiring.size() != overlay_dense.wiring.size()
      || overlay_runs.via.size() != overlay_dense.via.size()) {
    fail(t, "overlay sizes differ", c, steps);
  }
  for (const bool with_overlay : {false, true}) {
    FastGridReader ra;
    FastGridReader rb;
    ra.reset(&runs, with_overlay ? &overlay_runs : nullptr);
    rb.reset(&dense, with_overlay ? &overlay_dense : nullptr);
    for (int walk = 0; walk < 40; walk++) {
      int z = randomInt(rng, 0, num_z - 1);
      if (graph.layers[z].numTracks() == 0
          || graph.layers[z].numPositions() == 0) {
        continue;
      }
      int i = randomInt(rng, 0, graph.layers[z].numTracks() - 1);
      int j = randomInt(rng, 0, graph.layers[z].numPositions() - 1);
      for (int s = 0; s < 200; s++) {
        const TrackLayer& l = graph.layers[z];
        for (int wt = 0; wt <= fast; wt++) {
          for (const Kind k : {kWire, kJog, kViaUpPad, kViaDownPad, kCutUp}) {
            if (k == kCutUp && (z + 1 >= num_z || l.up_track[j] < 0)) {
              continue;
            }
            t.reader_checks++;
            if (ra.level(z, i, j, k, wt) != rb.level(z, i, j, k, wt)) {
              fail(t,
                   with_overlay ? "overlay reader level differs"
                                : "reader level differs",
                   c,
                   steps);
            }
          }
          for (const Kind e : {kWireEdge, kJogEdge}) {
            const bool along = e == kWireEdge;
            if (along ? j + 1 >= l.numPositions() : i + 1 >= l.numTracks()) {
              continue;
            }
            t.reader_checks++;
            if (ra.level(z, i, j, e, wt) != rb.level(z, i, j, e, wt)) {
              fail(t, "reader edge level differs", c, steps);
            }
          }
        }
        if (ra.wiringWord(z, i, j) != rb.wiringWord(z, i, j)) {
          fail(t, "reader word differs", c, steps);
        }
        if (z + 1 < num_z && l.up_track[j] >= 0) {
          const int cc = runs.viaIndex(z, j);
          if (ra.viaWord(z, i, cc) != rb.viaWord(z, i, cc)) {
            fail(t, "reader via word differs", c, steps);
          }
        }
        // Next step: mostly along the track, sometimes a jog, a via or a
        // jump.
        const int move = randomInt(rng, 0, 15);
        if (move < 10) {
          j += move < 7 ? 1 : -1;
        } else if (move < 12) {
          i += move == 10 ? 1 : -1;
        } else if (move < 14) {
          const int nz = z + (move == 12 ? 1 : -1);
          if (nz >= 0 && nz < num_z) {
            const int i2 = nz > z ? l.up_track[j] : l.down_track[j];
            const int j2 = nz > z ? l.up_pos[i] : l.down_pos[i];
            if (i2 >= 0 && j2 >= 0) {
              z = nz;
              i = i2;
              j = j2;
            }
          }
        } else {
          i = randomInt(rng, 0, graph.layers[z].numTracks() - 1);
          j = randomInt(rng, 0, graph.layers[z].numPositions() - 1);
        }
        i = std::clamp(i, 0, graph.layers[z].numTracks() - 1);
        j = std::clamp(j, 0, graph.layers[z].numPositions() - 1);
      }
    }
  }
  // The dense grid against the brute-force reference: vertex kinds of the
  // first wire type, with and without the overlay's owner removed.
  for (const bool with_overlay : {false, true}) {
    FastGridReader rd;
    rd.reset(&dense, with_overlay ? &overlay_dense : nullptr);
    const int32_t ignored = with_overlay ? own : kNoNet;
    const WireType& w = types[0];
    for (int z = 0; z < num_z; z++) {
      const TrackLayer& l = graph.layers[z];
      for (int i = 0; i < l.numTracks(); i++) {
        for (int j = 0; j < l.numPositions(); j++) {
          const auto [x, y] = l.point(i, j);
          const Point p{x, y};
          const RipupLevel wire = brute(
              classes,
              grid,
              live,
              LocatedModel{2 * z, w.preferred[z].rect, w.preferred[z].cls, 3, 0},
              p,
              ignored);
          const RipupLevel jog
              = brute(classes,
                      grid,
                      live,
                      LocatedModel{2 * z,
                                   w.non_preferred[z].rect,
                                   w.non_preferred[z].cls,
                                   3,
                                   0},
                      p,
                      ignored);
          t.brute_checks += 2;
          if (rd.level(z, i, j, kWire, 0) != wire
              || rd.level(z, i, j, kJog, 0) != jog) {
            fail(t,
                 with_overlay ? "dense overlay level differs from reference"
                              : "dense level differs from reference",
                 c,
                 steps);
          }
        }
      }
    }
  }
}

}  // namespace

int main()
{
  std::mt19937 rng(23);
  Totals t;
  const int cases = 300;
  for (int c = 0; c < cases && t.fails < 12; c++) {
    runCase(rng, c, t);
  }
  std::printf(
      "fast grid dense: %lld cases, %lld steps, %lld word checks (%lld "
      "non-zero), %lld run checks, %lld reader checks, %lld brute-force "
      "checks, %lld dense layers at the end, %lld encoding switches, %lld "
      "overlays, %lld large cases, %lld failures\n",
      static_cast<long long>(t.cases),
      static_cast<long long>(t.steps),
      static_cast<long long>(t.word_checks),
      static_cast<long long>(t.nonzero),
      static_cast<long long>(t.run_checks),
      static_cast<long long>(t.reader_checks),
      static_cast<long long>(t.brute_checks),
      static_cast<long long>(t.dense_layers),
      static_cast<long long>(t.switches),
      static_cast<long long>(t.overlays),
      static_cast<long long>(t.large),
      static_cast<long long>(t.fails));
  return t.fails != 0;
}
