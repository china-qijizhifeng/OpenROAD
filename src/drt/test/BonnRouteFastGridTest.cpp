// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Fast grid (Sec. 3.6) against the brute-force distance-rule reference
// after random insert / remove sequences: every on-track vertex kind and
// every edge of every fast wire type, edge "query" bits set exactly where
// vertex legality does not imply edge legality, incremental updates equal
// to a rebuild from scratch, and the checker path for other wire types.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <random>
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
  int64_t vertex_checks = 0;
  int64_t edge_checks = 0;
  int64_t marked = 0;
  int64_t blocked = 0;
  int64_t fails = 0;
  int64_t intervals = 0;
  int64_t via_intervals = 0;
  int64_t vertices = 0;
  int64_t overlays = 0;
};

void fail(Totals& t, const char* what, int c)
{
  if (t.fails++ < 12) {
    std::printf("case %d: %s\n", c, what);
  }
}

// Brute-force level of one model (rect relative to l), the shapes of
// `ignored` treated as removed.
thread_local int32_t g_ignored = kNoNet;
RipupLevel brute(const ShapeClassTable& classes,
                 const ShapeGrid& grid,
                 const std::vector<StoredShape>& shapes,
                 const LocatedModel& m,
                 const Point l)
{
  std::vector<LocatedModel> one{m};
  one[0].group = 0;
  return bruteLevels(classes, grid, shapes, one, 1, l, g_ignored)[0];
}

void runCase(std::mt19937& rng, const int c, Totals& t)
{
  t.cases++;
  const int num_z = randomInt(rng, 1, 3);
  const bool first_horizontal = rng() % 2;
  ShapeClassTable classes;
  classes.setLayers(randomRules(rng, num_z, first_horizontal));
  const Coord extent = randomInt(rng, 150, 320);
  // Track graph: random pitches and offsets, irregular gaps.
  TrackGraph graph;
  graph.layers.resize(num_z);
  std::vector<std::vector<Coord>> tracks(num_z);
  for (int z = 0; z < num_z; z++) {
    const int pitch = randomInt(rng, 16, 40);
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
  ChunkAxis cx{0, randomInt(rng, 40, 400), 0};
  cx.count = extent / cx.size + 1;
  ChunkAxis cy{0, randomInt(rng, 40, 400), 0};
  cy.count = extent / cy.size + 1;
  graph.finalize([&](Coord x) { return cx.of(x); },
                 [&](Coord y) { return cy.of(y); });
  ShapeGrid grid;
  grid.init(&classes, gridGeometry(rng, classes, extent, rng() % 2), cx, cy);
  const int nets = randomInt(rng, 1, 4);
  for (int n = 0; n < nets; n++) {
    grid.setNetRipupLevel(n, static_cast<RipupLevel>(randomInt(rng, 1, 6)));
  }
  // Wire types: fast ones plus one answered by the checker only.
  std::vector<WireType> types;
  const int fast = randomInt(rng, 1, 2);
  for (int k = 0; k < fast + 1; k++) {
    WireTypeGeometry g;
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
  FastGrid fg;
  fg.init(&graph, fast);
  FastGridBuilder builder;
  builder.init(&graph, &types, fast, &checker, &fg);
  ShapeGridWriteCache cache;
  std::vector<StoredShape> live;
  // Some shapes before the initial build, the rest incrementally.
  const int before = randomInt(rng, 0, 20);
  for (int k = 0; k < before; k++) {
    const StoredShape s = randomShape(
        rng, classes, randomInt(rng, 0, classes.numLayers() - 1), extent, nets);
    grid.insert(s.layer, s.rect, s.cls, s.type, s.owner, cache);
    live.push_back(s);
  }
  for (int z = 0; z < num_z; z++) {
    builder.rebuildVertices(z, 0, graph.layers[z].numTracks() - 1);
  }
  for (int z = 0; z < num_z; z++) {
    builder.rebuildJogBits(z, 0, graph.layers[z].numTracks() - 1);
  }
  const int steps = randomInt(rng, 0, 40);
  for (int k = 0; k < steps; k++) {
    if (live.empty() || rng() % 3 != 0) {
      StoredShape s = randomShape(rng,
                                  classes,
                                  randomInt(rng, 0, classes.numLayers() - 1),
                                  extent,
                                  nets);
      // On-track wires of routed nets, as the router creates them.
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
      // Small off-track shapes between two vertices (Fig. 4 zigzag edge).
      if (s.layer % 2 == 0 && rng() % 3 == 0 && graph.layers[z].numTracks() > 1
          && graph.layers[z].numPositions() > 1) {
        const TrackLayer& l = graph.layers[z];
        const int i = randomInt(rng, 0, l.numTracks() - 2);
        const int j = randomInt(rng, 0, l.numPositions() - 2);
        const auto [x0, y0] = l.point(i, j);
        const auto [x1, y1] = rng() % 2 ? l.point(i, j + 1) : l.point(i + 1, j);
        const Coord h = randomInt(rng, 1, 3);
        const Coord mx = (x0 + x1) / 2;
        const Coord my = (y0 + y1) / 2;
        s.rect = {mx - h, my - h, mx + h, my + h};
        s.cls = classes.classify(s.layer, s.rect, ShapeRole::kFixed);
      }
      grid.insert(s.layer, s.rect, s.cls, s.type, s.owner, cache);
      // Incremental insert (max with the shape's induced levels); now and
      // then the exact recomputation.
      if (rng() % 8 == 0) {
        builder.shapeChanged(s.layer, s.rect);
      } else {
        builder.shapeInserted(s.layer,
                              PlacedShape{s.rect, s.cls, 3},
                              grid.level(s.type, s.owner));
      }
      live.push_back(s);
    } else if (rng() % 3 == 0) {
      // A group of removals, recomputed once (a net's rip-up).
      FastGridBuilder::PendingRemovals pending;
      for (int r = randomInt(rng, 1, 4); r > 0 && !live.empty(); r--) {
        const size_t at = rng() % live.size();
        const StoredShape s = live[at];
        live.erase(live.begin() + at);
        grid.remove(s.layer, s.rect, s.cls, s.type, s.owner, cache);
        builder.shapeRemoved(s.layer,
                             PlacedShape{s.rect, s.cls, 3},
                             grid.level(s.type, s.owner),
                             pending);
      }
      builder.flush(pending);
    } else {
      const size_t at = rng() % live.size();
      const StoredShape s = live[at];
      live.erase(live.begin() + at);
      grid.remove(s.layer, s.rect, s.cls, s.type, s.owner, cache);
      builder.shapeRemoved(
          s.layer, PlacedShape{s.rect, s.cls, 3}, grid.level(s.type, s.owner));
    }
  }
  // Searches read through a reader (run cursors); with an overlay, the
  // shapes of one net count as removed (Sec. 4.4). Both are checked against
  // the brute-force reference, the base also against fg.level().
  FastGridReader rd;
  FastGridOverlay overlay;
  const int32_t own = randomInt(rng, 0, nets - 1);
  std::vector<FastGridBuilder::OwnShape> own_shapes;
  for (const StoredShape& s : live) {
    if (s.owner == own && isOwned(s.type)) {
      own_shapes.push_back({s.layer,
                            PlacedShape{s.rect, s.cls, 3},
                            grid.level(s.type, s.owner)});
    }
  }
  for (const bool with_overlay : {false, true}) {
    if (with_overlay) {
      builder.buildOverlay(own, own_shapes, overlay);
      rd.reset(&fg, &overlay);
      g_ignored = own;
    } else {
      rd.reset(&fg, nullptr);
      g_ignored = kNoNet;
    }
    // Every vertex and edge against the reference.
    auto model = [&](int layer, const WireModel& m) {
      return LocatedModel{layer, m.rect, m.cls, 3, 0};
    };
    for (int z = 0; z < num_z; z++) {
      const TrackLayer& l = graph.layers[z];
      for (int i = 0; i < l.numTracks(); i++) {
        for (int j = 0; j < l.numPositions(); j++) {
          const auto [x, y] = l.point(i, j);
          const Point p{x, y};
          for (int wt = 0; wt <= fast; wt++) {
            const WireType& w = types[wt];
            auto expect_vertex = [&](Kind k) -> RipupLevel {
              switch (k) {
                case kWire:
                  return brute(
                      classes, grid, live, model(2 * z, w.preferred[z]), p);
                case kJog:
                  return brute(
                      classes, grid, live, model(2 * z, w.non_preferred[z]), p);
                case kViaUpPad:
                  return z + 1 < num_z && w.via[z].valid
                             ? brute(classes,
                                     grid,
                                     live,
                                     model(2 * z, w.via[z].bottom),
                                     p)
                             : kLegal;
                case kViaDownPad:
                  return z > 0 && w.via[z - 1].valid
                             ? brute(classes,
                                     grid,
                                     live,
                                     model(2 * z, w.via[z - 1].top),
                                     p)
                             : kLegal;
                case kCutUp: {
                  if (z + 1 >= num_z || l.up_track[j] < 0 || !w.via[z].valid) {
                    return kLegal;
                  }
                  RipupLevel v = brute(
                      classes, grid, live, model(2 * z + 1, w.via[z].cut), p);
                  if (w.via[z].projection.valid) {
                    v = std::max(v,
                                 brute(classes,
                                       grid,
                                       live,
                                       model(2 * z + 3, w.via[z].projection),
                                       p));
                  }
                  return v;
                }
                default:
                  return kLegal;
              }
            };
            for (const Kind k : {kWire, kJog, kViaUpPad, kCutUp, kViaDownPad}) {
              if (k == kCutUp && (z + 1 >= num_z || l.up_track[j] < 0)) {
                continue;
              }
              const bool absent
                  = (k == kViaUpPad && (z + 1 >= num_z || !w.via[z].valid))
                    || (k == kViaDownPad && (z == 0 || !w.via[z - 1].valid))
                    || (k == kCutUp && !w.via[z].valid);
              const RipupLevel want = expect_vertex(k);
              const RipupLevel got = rd.level(z, i, j, k, wt);
              if (!with_overlay && got != fg.level(z, i, j, k, wt)) {
                fail(t, "reader differs from fast grid", c);
              }
              t.vertex_checks++;
              t.blocked += want != kLegal;
              // Absent models: the fast grid stores legal, the checker says
              // not available; the search never asks.
              if (!absent && got != want) {
                fail(t,
                     with_overlay ? "overlay vertex level differs"
                                  : "vertex level differs",
                     c);
              }
            }
            // Edges: level = max(endpoints, swept wire); bit exactness.
            for (const Kind e : {kWireEdge, kJogEdge}) {
              const bool along = e == kWireEdge;
              if (along ? j + 1 >= l.numPositions() : i + 1 >= l.numTracks()) {
                continue;
              }
              const WireModel& m = along ? w.preferred[z] : w.non_preferred[z];
              const Coord d = along ? l.positions[j + 1] - l.positions[j]
                                    : l.tracks[i + 1] - l.tracks[i];
              const bool shift_x = along == l.horizontal;
              LocatedModel swept = model(2 * z, m);
              swept.rect = m.rect.hull(
                  m.rect.translated(shift_x ? Point{d, 0} : Point{0, d}));
              const RipupLevel s = brute(classes, grid, live, swept, p);
              const auto [x2, y2]
                  = along ? l.point(i, j + 1) : l.point(i + 1, j);
              const RipupLevel ea
                  = brute(classes, grid, live, model(2 * z, m), p);
              const RipupLevel eb
                  = brute(classes, grid, live, model(2 * z, m), Point{x2, y2});
              const RipupLevel want = std::max({ea, eb, s});
              t.edge_checks++;
              if (!along && !l.jogs) {
                continue;
              }
              if (rd.level(z, i, j, e, wt) != want) {
                fail(t,
                     with_overlay ? "overlay edge level differs"
                                  : "edge level differs",
                     c);
              }
            }
          }
          // Bits: set iff some fast wire type's swept wire exceeds its
          // endpoints.
          for (const Kind e : {kWireEdge, kJogEdge}) {
            if (with_overlay) {
              continue;
            }
            const bool along = e == kWireEdge;
            if ((along ? j + 1 >= l.numPositions() : i + 1 >= l.numTracks())
                || (!along && !l.jogs)) {
              continue;
            }
            bool need = false;
            for (int wt = 0; wt < fast; wt++) {
              const WireModel& m
                  = along ? types[wt].preferred[z] : types[wt].non_preferred[z];
              const Coord d = along ? l.positions[j + 1] - l.positions[j]
                                    : l.tracks[i + 1] - l.tracks[i];
              const bool shift_x = along == l.horizontal;
              LocatedModel swept = model(2 * z, m);
              swept.rect = m.rect.hull(
                  m.rect.translated(shift_x ? Point{d, 0} : Point{0, d}));
              const auto [x2, y2]
                  = along ? l.point(i, j + 1) : l.point(i + 1, j);
              const RipupLevel ends = std::max(
                  brute(classes, grid, live, model(2 * z, m), p),
                  brute(classes, grid, live, model(2 * z, m), Point{x2, y2}));
              need |= brute(classes, grid, live, swept, p) > ends;
            }
            t.marked += need;
            if (fg.edgeQuery(z, i, j, e) != need) {
              fail(t, "edge query bit not exact", c);
            }
          }
        }
      }
    }
  }
  g_ignored = kNoNet;
  t.overlays += !own_shapes.empty();
  // Incremental updates equal a rebuild from scratch.
  FastGrid fresh;
  fresh.init(&graph, fast);
  FastGridBuilder fresh_builder;
  fresh_builder.init(&graph, &types, fast, &checker, &fresh);
  for (int z = 0; z < num_z; z++) {
    fresh_builder.rebuildVertices(z, 0, graph.layers[z].numTracks() - 1);
  }
  for (int z = 0; z < num_z; z++) {
    fresh_builder.rebuildJogBits(z, 0, graph.layers[z].numTracks() - 1);
  }
  for (int z = 0; z < num_z; z++) {
    const TrackLayer& l = graph.layers[z];
    for (int i = 0; i < l.numTracks(); i++) {
      for (int j = 0; j < l.numPositions(); j++) {
        uint64_t a = fg.wiringWord(z, i, j);
        uint64_t b = fresh.wiringWord(z, i, j);
        if (!l.jogs) {
          a &= ~(uint64_t{1} << FastGrid::kJogEdgeBit);
          b &= ~(uint64_t{1} << FastGrid::kJogEdgeBit);
        }
        if (a != b) {
          fail(t, "incremental word differs from rebuild", c);
        }
      }
      for (int cc = 0; z + 1 < num_z && cc < fg.numViaPositions(z); cc++) {
        if (fg.viaWord(z, i, cc) != fresh.viaWord(z, i, cc)) {
          fail(t, "incremental via word differs from rebuild", c);
        }
      }
    }
  }
  const FastGrid::Statistics st = fg.statistics();
  t.intervals += st.wiring_intervals;
  t.via_intervals += st.via_intervals;
  t.vertices += st.vertices;
}

}  // namespace

int main()
{
  std::mt19937 rng(11);
  Totals t;
  const int cases = 500;
  for (int c = 0; c < cases && t.fails < 12; c++) {
    runCase(rng, c, t);
  }
  std::printf(
      "fast grid: %lld cases, %lld vertex checks (%lld blocked), %lld edge "
      "checks, %lld marked edges, %lld failures | %lld vertices in %lld "
      "wiring intervals, %lld via intervals, %lld overlays checked\n",
      static_cast<long long>(t.cases),
      static_cast<long long>(t.vertex_checks),
      static_cast<long long>(t.blocked),
      static_cast<long long>(t.edge_checks),
      static_cast<long long>(t.marked),
      static_cast<long long>(t.fails),
      static_cast<long long>(t.vertices),
      static_cast<long long>(t.intervals),
      static_cast<long long>(t.via_intervals),
      static_cast<long long>(t.overlays));
  return t.fails != 0;
}
