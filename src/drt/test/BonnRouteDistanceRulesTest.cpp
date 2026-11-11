// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Distance rule checking module (Sec. 3.4) against the brute-force
// pairwise reference with the same shape classes: allowed models and ripup
// levels, maximal horizontal and vertical intervals (constant answer inside;
// the answer changes just outside unless the interval ends at the bound of
// the known area), blocking nets with ripup levels, the wire-type query
// interface and the placed-shape check.

#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "BonnRouteRoutingSpaceTestUtil.h"
#include "br/BrDistanceRules.h"
#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrWireModel.h"

using namespace drt::br;
using namespace drt::br::test;

namespace {

struct Totals
{
  int64_t queries = 0;
  int64_t fails = 0;
  int64_t blocked = 0;
  int64_t known_ends = 0;
  int64_t answer_ends = 0;
  int64_t interval_length = 0;
  int64_t nets_reported = 0;
};

void fail(Totals& t, const char* what, int c, int q)
{
  if (t.fails++ < 12) {
    std::printf("case %d query %d: %s\n", c, q, what);
  }
}

void runCase(std::mt19937& rng, const int c, Totals& t)
{
  ShapeClassTable classes;
  const int num_z = randomInt(rng, 1, 3);
  classes.setLayers(randomRules(rng, num_z, rng() % 2));
  const Coord extent = randomInt(rng, 120, 360);
  ShapeGrid grid;
  ChunkAxis cx{0, randomInt(rng, 40, 400), 0};
  cx.count = extent / cx.size + 1;
  ChunkAxis cy{0, randomInt(rng, 40, 400), 0};
  cy.count = extent / cy.size + 1;
  grid.init(&classes, gridGeometry(rng, classes, extent, rng() % 2), cx, cy);
  const int nets = randomInt(rng, 1, 5);
  for (int n = 0; n < nets; n++) {
    grid.setNetRipupLevel(n, static_cast<RipupLevel>(randomInt(rng, 1, 6)));
  }
  ShapeGridWriteCache cache;
  std::vector<StoredShape> live;
  const int count = randomInt(rng, 0, 40);
  for (int k = 0; k < count + 10; k++) {
    if (k < count || live.empty()) {
      const StoredShape s
          = randomShape(rng,
                        classes,
                        randomInt(rng, 0, classes.numLayers() - 1),
                        extent,
                        nets);
      grid.insert(s.layer, s.rect, s.cls, s.type, s.owner, cache);
      live.push_back(s);
    } else if (rng() % 2) {
      const size_t at = rng() % live.size();
      grid.remove(live[at].layer,
                  live[at].rect,
                  live[at].cls,
                  live[at].type,
                  live[at].owner,
                  cache);
      live.erase(live.begin() + at);
    }
  }
  // Wire types of random widths and vias.
  std::vector<WireType> types;
  const int num_types = randomInt(rng, 1, 3);
  for (int k = 0; k < num_types; k++) {
    WireTypeGeometry g;
    for (int z = 0; z < num_z; z++) {
      g.width.push_back(randomInt(rng, 6, 20));
    }
    for (int z = 0; z + 1 < num_z; z++) {
      WireTypeGeometry::Via v;
      v.valid = true;
      const Coord cut = randomInt(rng, 2, 6);
      v.cut = {-cut, -cut, cut, cut};
      const Coord b = cut + randomInt(rng, 0, 6);
      const Coord bx = b + randomInt(rng, 0, 6);
      v.bottom = {-bx, -b, bx, b};
      v.top = {-b, -bx, b, bx};
      g.via.push_back(v);
    }
    types.push_back(makeWireType(g, classes));
  }
  DistanceRuleChecker checker(&grid);
  for (int q = 0; q < 40; q++) {
    t.queries++;
    const Point l{randomInt(rng, -10, extent + 10),
                  randomInt(rng, -10, extent + 10)};
    const int z = randomInt(rng, 0, num_z - 1);
    std::vector<LocatedModel> models;
    for (int k = 0; k < num_types; k++) {
      appendWireTypeModels(types[k], z, num_z, 4 * k, models);
    }
    // Also an arbitrary model such as a swept edge (Sec. 3.6).
    if (rng() % 2) {
      const LocatedModel& base = models[rng() % models.size()];
      LocatedModel stretched = base;
      if (rng() % 2) {
        stretched.rect.xhi += randomInt(rng, 1, 40);
      } else {
        stretched.rect.yhi += randomInt(rng, 1, 40);
      }
      stretched.group = 4 * num_types;
      models.push_back(stretched);
    }
    const int groups = 4 * num_types + 1;
    CheckOptions opts;
    opts.levels = rng() % 3 != 0;
    opts.report_nets = true;
    opts.ignored_owner = rng() % 3 == 0 ? randomInt(rng, 0, nets - 1) : kNoNet;
    opts.margin_x = rng() % 2 ? randomInt(rng, 0, 60) : 0;
    opts.margin_y = rng() % 2 ? randomInt(rng, 0, 60) : 0;
    opts.groups = groups;
    CheckResult res;
    checker.check(l, models, opts, res);
    std::vector<BlockingNet> want_nets;
    const std::vector<RipupLevel> want = bruteLevels(
        classes, grid, live, models, groups, l, opts.ignored_owner, &want_nets);
    if (res.levels != want) {
      fail(t, "levels differ", c, q);
      continue;
    }
    for (const RipupLevel lv : want) {
      t.blocked += lv != 0;
    }
    std::vector<BlockingNet> got_nets = res.nets;
    if (got_nets != want_nets) {
      fail(t, "blocking nets differ", c, q);
    }
    t.nets_reported += static_cast<int64_t>(got_nets.size());
    // Removing every reported net makes each reported group legal.
    for (const BlockingNet& bn : got_nets) {
      if (bn.level < 1 || bn.level > kMaxRemovableLevel) {
        fail(t, "blocking net level out of range", c, q);
      }
    }
    // Intervals: constant answer inside, change just outside unless known
    // bound.
    auto at = [&](Point p) {
      return bruteLevels(
          classes, grid, live, models, groups, p, opts.ignored_owner);
    };
    if (res.x_lo > l.x || res.x_hi < l.x || res.y_lo > l.y || res.y_hi < l.y) {
      fail(t, "interval does not contain the location", c, q);
      continue;
    }
    const Coord lim = extent + 400;
    // Lengths inside the layout area (intervals reaching beyond the grid
    // are unbounded there).
    t.interval_length += int64_t{std::min(res.x_hi, lim)}
                         - std::max(res.x_lo, -lim) + std::min(res.y_hi, lim)
                         - std::max(res.y_lo, -lim);
    for (Coord x = std::max(res.x_lo, -lim); x <= std::min(res.x_hi, lim);
         x++) {
      if (!sameAnswer(at({x, l.y}), want, opts.levels)) {
        fail(t, "x interval not constant", c, q);
        break;
      }
    }
    for (Coord y = std::max(res.y_lo, -lim); y <= std::min(res.y_hi, lim);
         y++) {
      if (!sameAnswer(at({l.x, y}), want, opts.levels)) {
        fail(t, "y interval not constant", c, q);
        break;
      }
    }
    auto check_end = [&](bool known, Point outside, const char* what) {
      if (known) {
        t.known_ends++;
        return;
      }
      t.answer_ends++;
      if (outside.x < -lim || outside.x > lim || outside.y < -lim
          || outside.y > lim) {
        fail(t, "unbounded interval end without known bound", c, q);
        return;
      }
      if (sameAnswer(at(outside), want, opts.levels)) {
        fail(t, what, c, q);
      }
    };
    check_end(res.x_lo_known_bound, {res.x_lo - 1, l.y}, "x_lo not maximal");
    check_end(res.x_hi_known_bound, {res.x_hi + 1, l.y}, "x_hi not maximal");
    check_end(res.y_lo_known_bound, {l.x, res.y_lo - 1}, "y_lo not maximal");
    check_end(res.y_hi_known_bound, {l.x, res.y_hi + 1}, "y_hi not maximal");
    // Paper interface: per wire type and via model.
    std::vector<const WireType*> ptrs;
    for (const WireType& w : types) {
      ptrs.push_back(&w);
    }
    std::vector<WireTypeAnswer> answers;
    CheckResult qres;
    checker.query(z, l, ptrs, opts, answers, qres);
    for (int k = 0; k < num_types; k++) {
      if (answers[k].preferred != res.levels[4 * k]
          || answers[k].non_preferred != res.levels[4 * k + 1]
          || answers[k].via_up != res.levels[4 * k + 2]
          || answers[k].via_down != res.levels[4 * k + 3]) {
        fail(t, "wire type query differs", c, q);
      }
    }
    // Placed shape check.
    const LocatedModel& m = models[rng() % models.size()];
    std::vector<BlockingNet> shape_nets;
    const RipupLevel lv = checker.shapeLevel(m.layer,
                                             m.rect.translated(l),
                                             m.cls,
                                             m.line_end_ends,
                                             opts.ignored_owner,
                                             &shape_nets);
    std::vector<LocatedModel> one{m};
    one[0].group = 0;
    std::vector<BlockingNet> one_nets;
    const auto one_want = bruteLevels(
        classes, grid, live, one, 1, l, opts.ignored_owner, &one_nets);
    for (BlockingNet& n : one_nets) {
      n.groups = 1;
    }
    if (lv != one_want[0] || shape_nets != one_nets) {
      fail(t, "shape level differs", c, q);
    }
  }
}

}  // namespace

int main()
{
  std::mt19937 rng(31);
  Totals t;
  const int cases = 400;
  for (int c = 0; c < cases && t.fails < 12; c++) {
    runCase(rng, c, t);
  }
  std::printf(
      "distance rules: %d cases, %lld queries, %lld failures | %lld blocked "
      "groups, %lld reported nets | interval ends: %lld answer changes, %lld "
      "known-area bounds, mean length %.1f\n",
      cases,
      static_cast<long long>(t.queries),
      static_cast<long long>(t.fails),
      static_cast<long long>(t.blocked),
      static_cast<long long>(t.nets_reported),
      static_cast<long long>(t.answer_ends),
      static_cast<long long>(t.known_ends),
      t.queries ? static_cast<double>(t.interval_length) / (2.0 * t.queries)
                : 0.0);
  return t.fails != 0;
}
