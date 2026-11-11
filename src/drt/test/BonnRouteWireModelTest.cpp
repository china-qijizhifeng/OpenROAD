// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Wire models, shape classes and distance rules (Sec. 3.1, 3.2): sampled
// spacing tables reproduce nondecreasing rule functions with strict and
// non-strict thresholds exactly and dominate non-monotone ones, shape
// classes follow the line-end and run-length modelling, the pairwise rule
// is symmetric and translation invariant, and wire types carry the models
// of every layer.

#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "BonnRouteRoutingSpaceTestUtil.h"
#include "br/BrGeometry.h"
#include "br/BrWireModel.h"

using namespace drt::br;
using namespace drt::br::test;

namespace {

int spacingTables(std::mt19937& rng)
{
  int fails = 0;
  for (int c = 0; c < 500; c++) {
    std::vector<Coord> wt;
    std::vector<Coord> pt;
    std::vector<bool> ws;
    std::vector<bool> ps;
    for (int k = randomInt(rng, 0, 4); k > 0; k--) {
      wt.push_back(randomInt(rng, 1, 200));
      ws.push_back(rng() % 2);
    }
    for (int k = randomInt(rng, 0, 4); k > 0; k--) {
      pt.push_back(randomInt(rng, 1, 400));
      ps.push_back(rng() % 2);
    }
    const bool monotone = rng() % 4 != 0;
    std::vector<int> inc(64);
    for (int& v : inc) {
      v = monotone ? randomInt(rng, 0, 5) : randomInt(rng, -5, 5);
    }
    auto f = [&](Coord a, Coord b, Coord prl) {
      int s = 20;
      for (size_t k = 0; k < wt.size(); k++) {
        const Coord w = std::max(a, b);
        if (ws[k] ? w > wt[k] : w >= wt[k]) {
          s += inc[k];
        }
        if ((ws[k] ? std::min(a, b) > wt[k] : std::min(a, b) >= wt[k])) {
          s += inc[8 + k];
        }
      }
      for (size_t k = 0; k < pt.size(); k++) {
        if (ps[k] ? prl > pt[k] : prl >= pt[k]) {
          s += inc[16 + k];
        }
      }
      return std::max(0, s);
    };
    std::vector<Coord> wb;
    std::vector<Coord> pb;
    for (const Coord t : wt) {
      wb.push_back(t);
      wb.push_back(t + 1);
    }
    for (const Coord t : pt) {
      pb.push_back(t);
      pb.push_back(t + 1);
    }
    const SpacingTable table = SpacingTable::sample(wb, pb, f);
    for (int q = 0; q < 200; q++) {
      const Coord a = randomInt(rng, 0, 260);
      const Coord b = randomInt(rng, 0, 260);
      const Coord prl = randomInt(rng, 0, 500);
      const Coord got = table.at(
          table.widthRow(a), table.widthRow(b), table.prlColumn(prl));
      if (monotone ? got != f(a, b, prl) : got < f(a, b, prl)) {
        if (fails++ < 5) {
          std::printf("spacing table mismatch: case %d (%d %d %d): %d vs %d\n",
                      c,
                      a,
                      b,
                      prl,
                      got,
                      f(a, b, prl));
        }
      }
      // Nondecreasing in every argument.
      const Coord more = table.at(table.widthRow(a + randomInt(rng, 0, 40)),
                                  table.widthRow(b),
                                  table.prlColumn(prl + randomInt(rng, 0, 80)));
      if (more < got) {
        fails++;
      }
    }
  }
  return fails;
}

int classesAndRules(std::mt19937& rng)
{
  int fails = 0;
  for (int c = 0; c < 300; c++) {
    ShapeClassTable classes;
    const int num_z = randomInt(rng, 1, 3);
    classes.setLayers(randomRules(rng, num_z, rng() % 2));
    const Coord extent = 300;
    std::vector<StoredShape> shapes;
    for (int k = 0; k < 40; k++) {
      shapes.push_back(randomShape(
          rng, classes, randomInt(rng, 0, classes.numLayers() - 1), extent, 4));
    }
    for (const StoredShape& s : shapes) {
      const ShapeClass& sc = classes.get(s.cls);
      const LayerRules& rules = classes.layer(s.layer);
      if (sc.layer != s.layer
          || (sc.kind == ShapeClassKind::kMetal) == rules.via_layer) {
        fails++;
      }
      if (sc.line_end
          && (!rules.line_end.active()
              || (rules.horizontal ? s.rect.dy() : s.rect.dx())
                     >= rules.line_end.max_edge)) {
        fails++;
      }
      if (classes.reach(s.layer) < rules.spacing.maxValue()) {
        fails++;
      }
    }
    for (int q = 0; q < 400; q++) {
      const StoredShape& a = shapes[rng() % shapes.size()];
      const StoredShape& b = shapes[rng() % shapes.size()];
      if (a.layer != b.layer) {
        continue;
      }
      const PlacedShape pa{a.rect, a.cls, 3};
      const PlacedShape pb{b.rect, b.cls, 3};
      const bool v = shapesViolate(classes, pa, pb);
      if (v != shapesViolate(classes, pb, pa)) {
        fails++;
      }
      const Point d{randomInt(rng, -50, 50), randomInt(rng, -50, 50)};
      if (v
          != shapesViolate(classes,
                           PlacedShape{a.rect.translated(d), a.cls, 3},
                           PlacedShape{b.rect.translated(d), b.cls, 3})) {
        fails++;
      }
      // Beyond the reach no rule applies.
      const Coord r = classes.reach(a.layer);
      const Rect far = b.rect.translated(
          {a.rect.xhi - b.rect.xlo + r + randomInt(rng, 0, 10), 0});
      if (shapesViolate(classes, pa, PlacedShape{far, b.cls, 3})) {
        fails++;
      }
    }
    // Wire types.
    WireTypeGeometry g;
    for (int z = 0; z < num_z; z++) {
      g.width.push_back(randomInt(rng, 5, 20));
    }
    for (int z = 0; z + 1 < num_z; z++) {
      g.via.push_back({true, {-8, -6, 8, 6}, {-3, -3, 3, 3}, {-6, -8, 6, 8}});
    }
    const WireType t = makeWireType(g, classes);
    for (int z = 0; z < num_z; z++) {
      const WireModel& p = t.preferred[z];
      const WireModel& n = t.non_preferred[z];
      if (!p.valid || p.rect.dx() != g.width[z] || p.rect.dy() != g.width[z]
          || classes.get(n.cls).line_end) {
        fails++;
      }
      const ShapeClass& pc = classes.get(p.cls);
      const bool horizontal = classes.layer(wiringGridLayer(z)).horizontal;
      if ((horizontal ? pc.run_column_x : pc.run_column_y) != kUnboundedRun) {
        fails++;
      }
    }
    for (int z = 0; z + 1 < num_z; z++) {
      const ViaModel& v = t.via[z];
      const bool projection
          = viaGridLayer(z + 1) < classes.numLayers()
            && classes.layer(viaGridLayer(z + 1)).projection_spacing > 0;
      if (!v.valid || v.projection.valid != projection
          || classes.get(v.cut.cls).kind != ShapeClassKind::kCut
          || (projection
              && classes.get(v.projection.cls).kind
                     != ShapeClassKind::kCutProjection)) {
        fails++;
      }
    }
  }
  if (fails) {
    std::printf("class / rule failures: %d\n", fails);
  }
  return fails;
}

}  // namespace

int main()
{
  std::mt19937 rng(5);
  const int fails = spacingTables(rng) + classesAndRules(rng);
  std::printf("wire models: %d failures\n", fails);
  return fails != 0;
}
