// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Exactness test of the goal-oriented potentials (BrPotential.h) against a
// plain multi-source Dijkstra on the explicit integer grid: SimplePotential
// and GeneralPotential (both query modes) must equal the distance at every
// vertex, L1Potential must be a feasible lower bound.
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <queue>
#include <random>
#include <utility>
#include <vector>

#include "br/BrCore.h"
#include "br/BrPotential.h"

using namespace drt::br;

namespace {

int failures = 0;

#define CHECK(cond, ...)                                            \
  do {                                                              \
    if (!(cond)) {                                                  \
      failures++;                                                   \
      if (failures <= 20) {                                         \
        std::printf("FAIL %s:%d: %s: ", __FILE__, __LINE__, #cond); \
        std::printf(__VA_ARGS__);                                   \
        std::printf("\n");                                          \
      }                                                             \
    }                                                               \
  } while (0)

Cost add(Cost a, Cost b)
{
  return (a >= kInfCost || b >= kInfCost) ? kInfCost
                                          : std::min(a + b, kInfCost);
}

// Plain Dijkstra on the grid [lo, hi]^2 x {0..l-1} with the closed-tile
// minimum rule evaluated directly from the model.
struct Reference
{
  const TileCostModel& model;
  Coord lo, hi;
  int layers;
  int width;
  std::vector<Cost> east;   // edge (x, y) - (x + 1, y)
  std::vector<Cost> north;  // edge (x, y) - (x, y + 1)
  std::vector<Cost> up;     // edge z - z + 1
  std::vector<Cost> dist;

  Reference(const TileCostModel& m, Coord l, Coord h)
      : model(m), lo(l), hi(h), layers(m.num_layers), width(h - l + 1)
  {
    const size_t n = size_t(layers) * width * width;
    east.assign(n, kInfCost);
    north.assign(n, kInfCost);
    up.assign(n, kInfCost);
    std::vector<std::pair<int, int>> col_point(width), col_edge(width);
    std::vector<std::pair<int, int>> row_point(width), row_edge(width);
    for (Coord c = lo; c <= hi; c++) {
      col_point[c - lo] = columns(c, c);
      col_edge[c - lo] = columns(c, c + 1);
      row_point[c - lo] = rows(c, c);
      row_edge[c - lo] = rows(c, c + 1);
    }
    auto minimum = [&](const std::vector<Cost>& rates,
                       int z,
                       std::pair<int, int> cr,
                       std::pair<int, int> rr) {
      Cost best = kInfCost;
      for (int i = cr.first; i <= cr.second; i++) {
        for (int j = rr.first; j <= rr.second; j++) {
          best = std::min(best, rates[model.index(z, i, j)]);
        }
      }
      return best;
    };
    for (int z = 0; z < layers; z++) {
      for (Coord x = lo; x <= hi; x++) {
        for (Coord y = lo; y <= hi; y++) {
          const size_t v = index(z, x, y);
          if (x < hi) {
            east[v] = minimum(
                model.horizontal, z, col_edge[x - lo], row_point[y - lo]);
          }
          if (y < hi) {
            north[v] = minimum(
                model.vertical, z, col_point[x - lo], row_edge[y - lo]);
          }
          if (z + 1 < layers) {
            up[v] = minimum(
                model.via_up, z, col_point[x - lo], row_point[y - lo]);
          }
        }
      }
    }
  }
  size_t index(int z, Coord x, Coord y) const
  {
    return (size_t(z) * width + (x - lo)) * width + (y - lo);
  }
  int64_t columnLo(int i) const
  {
    return i == 0 ? INT64_MIN : model.x_cuts[i - 1];
  }
  int64_t columnHi(int i) const
  {
    return i == int(model.x_cuts.size()) ? INT64_MAX : model.x_cuts[i];
  }
  int64_t rowLo(int j) const
  {
    return j == 0 ? INT64_MIN : model.y_cuts[j - 1];
  }
  int64_t rowHi(int j) const
  {
    return j == int(model.y_cuts.size()) ? INT64_MAX : model.y_cuts[j];
  }
  // Range of closed columns / rows containing [a, b] (contiguous).
  std::pair<int, int> columns(Coord a, Coord b) const
  {
    int first = INT_MAX, last = -1;
    for (int i = 0; i < model.numColumns(); i++) {
      if (columnLo(i) <= a && b <= columnHi(i)) {
        first = std::min(first, i);
        last = std::max(last, i);
      }
    }
    return {first, last};
  }
  std::pair<int, int> rows(Coord a, Coord b) const
  {
    int first = INT_MAX, last = -1;
    for (int j = 0; j < model.numRows(); j++) {
      if (rowLo(j) <= a && b <= rowHi(j)) {
        first = std::min(first, j);
        last = std::max(last, j);
      }
    }
    return {first, last};
  }
  void run(const std::vector<TargetRectangle>& targets)
  {
    dist.assign(east.size(), kInfCost);
    using Item = std::pair<Cost, size_t>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> heap;
    for (const TargetRectangle& t : targets) {
      for (Coord x = std::max(t.xlo, lo); x <= std::min(t.xhi, hi); x++) {
        for (Coord y = std::max(t.ylo, lo); y <= std::min(t.yhi, hi); y++) {
          const size_t v = index(t.z, x, y);
          if (t.terminal_cost < dist[v]) {
            dist[v] = t.terminal_cost;
            heap.push({t.terminal_cost, v});
          }
        }
      }
    }
    while (!heap.empty()) {
      const auto [d, v] = heap.top();
      heap.pop();
      if (d != dist[v]) {
        continue;
      }
      forEachEdge(v, [&](size_t w, Cost c) {
        const Cost nd = add(d, c);
        if (nd < dist[w]) {
          dist[w] = nd;
          heap.push({nd, w});
        }
      });
    }
  }
  template <typename F>
  void forEachEdge(size_t v, F&& fn) const
  {
    const int z = int(v / (size_t(width) * width));
    const int x = int((v / width) % width);
    const int y = int(v % width);
    if (x + 1 < width && east[v] < kInfCost) {
      fn(v + width, east[v]);
    }
    if (x > 0 && east[v - width] < kInfCost) {
      fn(v - width, east[v - width]);
    }
    if (y + 1 < width && north[v] < kInfCost) {
      fn(v + 1, north[v]);
    }
    if (y > 0 && north[v - 1] < kInfCost) {
      fn(v - 1, north[v - 1]);
    }
    const size_t plane = size_t(width) * width;
    if (z + 1 < layers && up[v] < kInfCost) {
      fn(v + plane, up[v]);
    }
    if (z > 0 && up[v - plane] < kInfCost) {
      fn(v - plane, up[v - plane]);
    }
  }
};

struct Potentials
{
  L1Potential l1;
  SimplePotential simple;
  GeneralPotential general_boundary;
  GeneralPotential general_tiles;
  Potentials()
  {
    general_boundary.setQueryMode(GeneralPotential::QueryMode::kBoundary);
    general_tiles.setQueryMode(GeneralPotential::QueryMode::kTileFunctions);
  }
};

// Compares all potentials with the reference on the box [lo, hi]^2; the
// model's domain may be smaller (queries outside it must still be feasible
// lower bounds, and exact inside).
void checkInstance(Potentials& p,
                   const TileCostModel& model,
                   const std::vector<TargetRectangle>& targets,
                   const LayerCostRates* simple_rates,
                   Coord lo,
                   Coord hi,
                   const char* name,
                   int id)
{
  Reference ref(model, lo, hi);
  ref.run(targets);
  const LayerCostRates rates = minimumLayerRates(model);
  p.l1.prepare(rates, targets);
  p.general_boundary.prepare(model, targets);
  p.general_tiles.prepare(model, targets);
  if (simple_rates) {
    p.simple.prepare(*simple_rates, targets);
  }
  auto inside = [&](Coord x, Coord y) {
    return model.domain_xlo <= x && x <= model.domain_xhi
           && model.domain_ylo <= y && y <= model.domain_yhi;
  };
  const size_t n = ref.dist.size();
  std::vector<Cost> l1(n), general(n), tiles(n);
  for (int z = 0; z < model.num_layers; z++) {
    for (Coord x = lo; x <= hi; x++) {
      for (Coord y = lo; y <= hi; y++) {
        const size_t v = ref.index(z, x, y);
        const Cost expected = ref.dist[v];
        l1[v] = p.l1.query(z, x, y);
        general[v] = p.general_boundary.query(z, x, y);
        tiles[v] = p.general_tiles.query(z, x, y);
        CHECK(l1[v] <= expected,
              "%s %d: l1 %lld > dist %lld at (%d,%d,%d)",
              name,
              id,
              (long long) l1[v],
              (long long) expected,
              x,
              y,
              z);
        if (inside(x, y)) {
          CHECK(general[v] == expected,
                "%s %d: general/boundary %lld != dist %lld at (%d,%d,%d)",
                name,
                id,
                (long long) general[v],
                (long long) expected,
                x,
                y,
                z);
          CHECK(tiles[v] == expected,
                "%s %d: general/tiles %lld != dist %lld at (%d,%d,%d)",
                name,
                id,
                (long long) tiles[v],
                (long long) expected,
                x,
                y,
                z);
        } else {
          CHECK(general[v] <= expected && tiles[v] == general[v],
                "%s %d: projected general %lld / %lld vs dist %lld at "
                "(%d,%d,%d)",
                name,
                id,
                (long long) general[v],
                (long long) tiles[v],
                (long long) expected,
                x,
                y,
                z);
        }
        if (simple_rates) {
          const Cost s = p.simple.query(z, x, y);
          CHECK(s == expected,
                "%s %d: simple %lld != dist %lld at (%d,%d,%d)",
                name,
                id,
                (long long) s,
                (long long) expected,
                x,
                y,
                z);
        }
      }
    }
  }
  // Feasibility: pi(v) <= c(v, w) + pi(w) for every edge, both directions.
  for (size_t v = 0; v < n; v++) {
    ref.forEachEdge(v, [&](size_t w, Cost c) {
      CHECK(l1[v] <= add(c, l1[w]),
            "%s %d: l1 infeasible edge %zu-%zu",
            name,
            id,
            v,
            w);
      CHECK(general[v] <= add(c, general[w]),
            "%s %d: general infeasible edge %zu-%zu",
            name,
            id,
            v,
            w);
    });
  }
}

std::vector<Coord> randomCuts(std::mt19937& rng, Coord n, int max_count)
{
  std::vector<Coord> cuts;
  const int count = std::uniform_int_distribution<int>(0, max_count)(rng);
  for (int k = 0; k < count; k++) {
    cuts.push_back(std::uniform_int_distribution<Coord>(0, n)(rng));
  }
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
  std::vector<Coord> out;
  for (const Coord c : cuts) {
    out.push_back(c);
    if (rng() % 4 == 0) {
      out.push_back(c);  // degenerate tile
    }
  }
  return out;
}

Cost randomRate(std::mt19937& rng, int inf_percent, int zero_percent, int max)
{
  const int r = int(rng() % 100);
  if (r < inf_percent) {
    return kInfCost;
  }
  if (r < inf_percent + zero_percent) {
    return 0;
  }
  return std::uniform_int_distribution<Cost>(1, max)(rng);
}

std::vector<TargetRectangle> randomTargets(std::mt19937& rng,
                                           Coord n,
                                           int layers,
                                           int max_count,
                                           bool points)
{
  std::vector<TargetRectangle> targets;
  const int count = std::uniform_int_distribution<int>(1, max_count)(rng);
  for (int k = 0; k < count; k++) {
    TargetRectangle t;
    t.xlo = std::uniform_int_distribution<Coord>(0, n)(rng);
    t.ylo = std::uniform_int_distribution<Coord>(0, n)(rng);
    t.xhi = points
                ? t.xlo
                : std::min<Coord>(n,
                                  t.xlo
                                      + std::uniform_int_distribution<Coord>(
                                          0, std::max<Coord>(n / 3, 0))(rng));
    t.yhi = points
                ? t.ylo
                : std::min<Coord>(n,
                                  t.ylo
                                      + std::uniform_int_distribution<Coord>(
                                          0, std::max<Coord>(n / 3, 0))(rng));
    t.z = std::uniform_int_distribution<int>(0, layers - 1)(rng);
    t.terminal_cost
        = rng() % 2 ? 0 : std::uniform_int_distribution<Cost>(0, 12)(rng);
    targets.push_back(t);
  }
  return targets;
}

LayerCostRates randomLayerRates(std::mt19937& rng, int layers)
{
  LayerCostRates rates;
  for (int z = 0; z < layers; z++) {
    rates.horizontal.push_back(randomRate(rng, 15, 5, 9));
    rates.vertical.push_back(randomRate(rng, 15, 5, 9));
    if (z + 1 < layers) {
      rates.via_up.push_back(randomRate(rng, 15, 5, 12));
    }
  }
  return rates;
}

TileCostModel randomGeneralModel(std::mt19937& rng, Coord n, int layers)
{
  TileCostModel model;
  model.x_cuts = randomCuts(rng, n, 6);
  model.y_cuts = randomCuts(rng, n, 6);
  model.domain_xlo = 0;
  model.domain_ylo = 0;
  model.domain_xhi = n;
  model.domain_yhi = n;
  model.resize(layers, kInfCost);
  auto degenerate = [](const std::vector<Coord>& cuts, int i) {
    return i >= 1 && i < int(cuts.size()) && cuts[i - 1] == cuts[i];
  };
  for (int z = 0; z < layers; z++) {
    for (int i = 0; i < model.numColumns(); i++) {
      for (int j = 0; j < model.numRows(); j++) {
        const size_t k = model.index(z, i, j);
        if (degenerate(model.x_cuts, i) || degenerate(model.y_cuts, j)) {
          // Mostly neutral; sometimes a discounted reservation line.
          model.horizontal[k] = randomRate(rng, 65, 5, 3);
          model.vertical[k] = randomRate(rng, 65, 5, 3);
          model.via_up[k] = randomRate(rng, 80, 5, 4);
        } else {
          model.horizontal[k] = randomRate(rng, 15, 3, 9);
          model.vertical[k] = randomRate(rng, 15, 3, 9);
          model.via_up[k] = randomRate(rng, 20, 3, 12);
        }
      }
    }
  }
  return model;
}

// Paper Figure 7: two horizontally adjacent tiles, target (0, 0).
void figureSeven(Potentials& p)
{
  TileCostModel model;
  model.x_cuts = {0, 4, 7};
  model.y_cuts = {0, 1};
  model.domain_xlo = 0;
  model.domain_ylo = 0;
  model.domain_xhi = 7;
  model.domain_yhi = 1;
  model.resize(1, kInfCost);
  model.horizontal[model.index(0, 1, 1)] = 2;
  model.vertical[model.index(0, 1, 1)] = 20;
  model.horizontal[model.index(0, 2, 1)] = 1;
  model.vertical[model.index(0, 2, 1)] = 10;
  const std::vector<TargetRectangle> targets = {{0, 0, 0, 0, 0, 0}};
  const Cost bottom[8] = {0, 2, 4, 6, 8, 9, 10, 11};
  const Cost top[8] = {20, 22, 22, 20, 18, 19, 20, 21};
  for (GeneralPotential* g : {&p.general_boundary, &p.general_tiles}) {
    g->prepare(model, targets);
    for (Coord x = 0; x <= 7; x++) {
      CHECK(g->query(0, x, 0) == bottom[x], "figure 7 bottom x=%d", x);
      CHECK(g->query(0, x, 1) == top[x], "figure 7 top x=%d", x);
    }
  }
  checkInstance(p, model, targets, nullptr, 0, 7, "figure7", 0);
}

// A pocket of finite tiles that cannot reach the target.
void unreachable(Potentials& p)
{
  TileCostModel model;
  model.x_cuts = {0, 10, 20, 30};
  model.y_cuts = {0, 10};
  model.domain_xlo = 0;
  model.domain_ylo = 0;
  model.domain_xhi = 30;
  model.domain_yhi = 10;
  model.resize(2, kInfCost);
  for (int z = 0; z < 2; z++) {
    // Tiles of columns 1 and 3 are routable, column 2 is blocked.
    for (int i : {1, 3}) {
      model.horizontal[model.index(z, i, 1)] = 1 + z;
      model.vertical[model.index(z, i, 1)] = 2 - z;
    }
  }
  model.via_up[model.index(0, 1, 1)] = 5;
  model.via_up[model.index(0, 3, 1)] = 5;
  // Column 2 is blocked; the cut lines x = 10 and x = 20 still carry the
  // vertical edges of their routable neighbours (closed tiles).
  const std::vector<TargetRectangle> targets = {{2, 3, 4, 3, 1, 7}};
  for (GeneralPotential* g : {&p.general_boundary, &p.general_tiles}) {
    g->prepare(model, targets);
    CHECK(g->query(0, 25, 5) == kInfCost, "pocket must be unreachable");
    CHECK(g->query(1, 15, 5) == kInfCost, "blocked column");
    CHECK(g->query(1, 3, 3) == 7, "inside the target");
    CHECK(g->query(0, 3, 3) == 12, "via below the target");
  }
  checkInstance(p, model, targets, nullptr, 0, 30, "unreachable", 0);
  // Everything forbidden.
  TileCostModel blocked;
  blocked.domain_xhi = blocked.domain_yhi = 5;
  blocked.resize(3, kInfCost);
  checkInstance(p, blocked, targets, nullptr, 0, 5, "blocked", 0);
  // No targets.
  p.general_tiles.prepare(model, {});
  CHECK(p.general_tiles.query(0, 3, 3) == kInfCost, "no targets");
}

int randomGeneral(Potentials& p,
                  int cases,
                  uint32_t seed,
                  int min_layers,
                  int max_layers,
                  bool points,
                  int max_targets,
                  const char* name)
{
  std::mt19937 rng(seed);
  for (int c = 0; c < cases && failures == 0; c++) {
    const Coord n = rng() % 3 == 0
                        ? 40
                        : std::uniform_int_distribution<Coord>(1, 24)(rng);
    const int layers
        = std::uniform_int_distribution<int>(min_layers, max_layers)(rng);
    const TileCostModel model = randomGeneralModel(rng, n, layers);
    const auto targets
        = randomTargets(rng, n, layers, points ? 1 : max_targets, points);
    // Every fourth case also checks the projection outside the domain.
    const Coord margin = c % 4 == 0 ? 3 : 0;
    checkInstance(p, model, targets, nullptr, -margin, n + margin, name, c);
  }
  return failures;
}

int randomSimple(Potentials& p,
                 int cases,
                 uint32_t seed,
                 int min_layers,
                 int max_layers,
                 bool points,
                 int max_targets,
                 const char* name)
{
  std::mt19937 rng(seed);
  for (int c = 0; c < cases && failures == 0; c++) {
    const Coord n = rng() % 3 == 0
                        ? 40
                        : std::uniform_int_distribution<Coord>(1, 24)(rng);
    const int layers
        = std::uniform_int_distribution<int>(min_layers, max_layers)(rng);
    const LayerCostRates rates = randomLayerRates(rng, layers);
    TileCostModel model;
    model.assignSimple(rates, 0, 0, n, n);
    const auto targets
        = randomTargets(rng, n, layers, points ? 1 : max_targets, points);
    checkInstance(p, model, targets, &rates, 0, n, name, c);
  }
  return failures;
}

// Builder: the model's edge costs must follow the Region / reservation
// semantics, and the potential must be exact on the built model.
int builderCases(Potentials& p, int cases, uint32_t seed)
{
  std::mt19937 rng(seed);
  TileCostModelBuilder builder;
  TileCostModel model;
  for (int c = 0; c < cases && failures == 0; c++) {
    const int layers = std::uniform_int_distribution<int>(1, 4)(rng);
    const int tiles = std::uniform_int_distribution<int>(1, 5)(rng);
    GCellGrid grid;
    grid.origin_x = grid.origin_y = 0;
    grid.pitch_x = std::uniform_int_distribution<Coord>(3, 8)(rng);
    grid.pitch_y = std::uniform_int_distribution<Coord>(3, 8)(rng);
    grid.die_xhi = grid.pitch_x * tiles;
    grid.die_yhi = grid.pitch_y * tiles;
    Region region;
    region.init(layers, tiles, tiles);
    for (int z = 0; z < layers; z++) {
      for (int tx = 0; tx < tiles; tx++) {
        for (int ty = 0; ty < tiles; ty++) {
          if (rng() % 3 != 0) {
            region.allow(z, tx, ty);
          }
        }
      }
    }
    LayerCostRates base;
    for (int z = 0; z < layers; z++) {
      // Scaled by the denominator 4 so that 3/4 discounts stay exact.
      const Cost pref = 4 * std::uniform_int_distribution<Cost>(1, 3)(rng);
      base.horizontal.push_back(z % 2 ? 10 * pref : pref);
      base.vertical.push_back(z % 2 ? pref : 10 * pref);
      if (z + 1 < layers) {
        base.via_up.push_back(4
                              * std::uniform_int_distribution<Cost>(1, 9)(rng));
      }
    }
    std::vector<DiscountedWire> wires;
    const int num_wires = std::uniform_int_distribution<int>(0, 6)(rng);
    const Coord nx = grid.die_xhi, ny = grid.die_yhi;
    for (int k = 0; k < num_wires; k++) {
      DiscountedWire w;
      w.z = std::uniform_int_distribution<int>(0, layers - 1)(rng);
      w.numerator = rng() % 2 ? 3 : 1;
      w.denominator = rng() % 2 ? 4 : 2;
      if (rng() % 2) {
        w.y1 = w.y2 = std::uniform_int_distribution<Coord>(0, ny)(rng);
        w.x1 = std::uniform_int_distribution<Coord>(0, nx)(rng);
        w.x2 = std::uniform_int_distribution<Coord>(0, nx)(rng);
      } else {
        w.x1 = w.x2 = std::uniform_int_distribution<Coord>(0, nx)(rng);
        w.y1 = std::uniform_int_distribution<Coord>(0, ny)(rng);
        w.y2 = std::uniform_int_distribution<Coord>(0, ny)(rng);
      }
      wires.push_back(w);
    }
    builder.build(region, grid, base, wires, model);
    if (region.bx_hi < 0) {
      continue;
    }
    // Independent edge cost reference on the window.
    const Coord lo = std::min(model.domain_xlo, model.domain_ylo);
    const Coord hi = std::max(model.domain_xhi, model.domain_yhi);
    Reference ref(model, lo, hi);
    auto gcells
        = [](Coord a, Coord b, Coord pitch, int count, int* first, int* last) {
            // Closed GCells [t * pitch, (t + 1) * pitch] (outer ones unbounded)
            // containing [a, b].
            *first = INT_MAX;
            *last = -1;
            for (int t = 0; t < count; t++) {
              const int64_t glo = t == 0 ? INT64_MIN : int64_t(t) * pitch;
              const int64_t ghi
                  = t == count - 1 ? INT64_MAX : int64_t(t + 1) * pitch;
              if (glo <= a && b <= ghi) {
                *first = std::min(*first, t);
                *last = std::max(*last, t);
              }
            }
          };
    for (int z = 0; z < layers; z++) {
      for (Coord x = lo; x <= hi; x++) {
        for (Coord y = lo; y <= hi; y++) {
          const size_t v = ref.index(z, x, y);
          for (int kind = 0; kind < 3; kind++) {
            if ((kind == 0 && x == hi) || (kind == 1 && y == hi)
                || (kind == 2 && z + 1 == layers)) {
              continue;
            }
            const Coord x1 = kind == 0 ? x + 1 : x;
            const Coord y1 = kind == 1 ? y + 1 : y;
            int fx, lx, fy, ly;
            gcells(x, x1, grid.pitch_x, tiles, &fx, &lx);
            gcells(y, y1, grid.pitch_y, tiles, &fy, &ly);
            Cost expected = kInfCost;
            for (int tx = fx; tx <= lx; tx++) {
              for (int ty = fy; ty <= ly; ty++) {
                if (!region.allowedTile(z, tx, ty)) {
                  continue;
                }
                if (kind == 2) {
                  if (region.allowedTile(z + 1, tx, ty)) {
                    expected = std::min(expected, base.via_up[z]);
                  }
                } else {
                  expected = std::min(
                      expected,
                      kind == 0 ? base.horizontal[z] : base.vertical[z]);
                }
              }
            }
            const Cost undiscounted = expected;
            for (const DiscountedWire& w : wires) {
              if (w.z != z || undiscounted >= kInfCost || kind == 2) {
                continue;
              }
              const bool horizontal_wire = w.y1 == w.y2 && w.x1 != w.x2;
              const bool vertical_wire = w.x1 == w.x2 && w.y1 != w.y2;
              if (kind == 0 && horizontal_wire && w.y1 == y
                  && std::min(w.x1, w.x2) <= x && x1 <= std::max(w.x1, w.x2)) {
                expected = std::min(expected,
                                    undiscounted * w.numerator / w.denominator);
              }
              if (kind == 1 && vertical_wire && w.x1 == x
                  && std::min(w.y1, w.y2) <= y && y1 <= std::max(w.y1, w.y2)) {
                expected = std::min(expected,
                                    undiscounted * w.numerator / w.denominator);
              }
            }
            const Cost actual = kind == 0   ? ref.east[v]
                                : kind == 1 ? ref.north[v]
                                            : ref.up[v];
            const bool in_window
                = x >= model.domain_xlo && x1 <= model.domain_xhi
                  && y >= model.domain_ylo && y1 <= model.domain_yhi;
            if (in_window) {
              CHECK(actual == expected,
                    "builder %d: kind %d at (%d,%d,%d): model %lld expected "
                    "%lld",
                    c,
                    kind,
                    x,
                    y,
                    z,
                    (long long) actual,
                    (long long) expected);
            }
          }
        }
      }
    }
    const auto targets = randomTargets(
        rng, std::min(model.domain_xhi, model.domain_yhi), layers, 3, false);
    checkInstance(p, model, targets, nullptr, lo, hi, "builder", c);
  }
  return failures;
}

}  // namespace

int main()
{
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  Potentials p;
  figureSeven(p);
  unreachable(p);
  if (failures) {
    return 1;
  }
  struct Suite
  {
    const char* name;
    int cases;
    bool simple;
    int min_layers, max_layers;
    bool points;
    int max_targets;
    uint32_t seed;
  };
  const Suite suites[] = {
      {"simple", 1500, true, 2, 5, false, 4, 11},
      {"simple-1layer", 200, true, 1, 1, false, 4, 12},
      {"simple-point", 300, true, 1, 5, true, 1, 13},
      {"simple-many", 200, true, 2, 5, false, 24, 14},
      {"general", 3000, false, 2, 5, false, 4, 21},
      {"general-1layer", 300, false, 1, 1, false, 4, 22},
      {"general-point", 500, false, 1, 5, true, 1, 23},
      {"general-many", 200, false, 2, 5, false, 24, 24},
  };
  for (const Suite& s : suites) {
    const int before = failures;
    if (s.simple) {
      randomSimple(p,
                   s.cases,
                   s.seed,
                   s.min_layers,
                   s.max_layers,
                   s.points,
                   s.max_targets,
                   s.name);
    } else {
      randomGeneral(p,
                    s.cases,
                    s.seed,
                    s.min_layers,
                    s.max_layers,
                    s.points,
                    s.max_targets,
                    s.name);
    }
    std::printf(
        "%-15s cases=%d failures=%d\n", s.name, s.cases, failures - before);
    if (failures) {
      return 1;
    }
  }
  builderCases(p, 400, 31);
  std::printf("%-15s cases=%d failures=%d\n", "builder", 400, failures);
  const PotentialCounters& g = p.general_tiles.counters();
  std::printf(
      "general/tiles counters: prepares=%lld queries=%lld offered=%lld "
      "inserted=%lld removed=%lld pops=%lld tile_lists=%lld "
      "tile_candidates=%lld tile_functions=%lld\n",
      (long long) g.prepares,
      (long long) g.queries,
      (long long) g.functions_offered,
      (long long) g.functions_inserted,
      (long long) g.functions_removed,
      (long long) g.heap_pops,
      (long long) g.tile_lists_built,
      (long long) g.tile_candidates,
      (long long) g.tile_functions);
  return failures ? 1 : 0;
}
