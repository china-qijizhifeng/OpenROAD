// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Differential test: interval search == vertex search == plain Dijkstra, all
// reading legality through the fast grid API (Sec. 3.6) including edges
// marked for the distance rule checker and wire types outside the fast
// grid.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

#include "br/BrCore.h"
#include "br/BrDistanceRules.h"
#include "br/BrPotentialOracle.h"
#include "br/BrShapeGrid.h"
#include "br/BrWireModel.h"

using namespace drt::br;

// Stands in for the distance rule checker: exact levels of marked edges,
// and wire type t >= fast answered like fast wire type t - 1.
struct TestOracle : LegalityOracle
{
  const FastGrid* grid = nullptr;
  std::map<std::pair<size_t, int>, RipupLevel> edges;
  RipupLevel vertexLevel(int z, int i, int j, int t, Kind k, int32_t)
      const override
  {
    return grid->level(z, i, j, k, t - 1);
  }
  RipupLevel edgeLevel(int z, int i, int j, int t, Kind k, int32_t)
      const override
  {
    if (t >= grid->numWireTypes()) {
      return grid->level(z, i, j, k, t - 1);
    }
    const auto it = edges.find({grid->graph().vertex(z, i, j), k});
    return it == edges.end() ? kNotRemovable : it->second;
  }
};

struct Ref
{
  const TrackGraph& g;
  const FastGrid& f;
  const SearchRequest& r;
  bool enter(size_t v, Kind k) const
  {
    const RipupLevel l = f.level(v, k, r.wire_type);
    return l == kLegal || (r.ripup && l <= r.max_ripup_level);
  }
  Cost pen(size_t v, Kind k) const
  {
    return f.level(v, k, r.wire_type) == kLegal
               ? 0
               : r.ripup_cost * (1 + f.history(v));
  }
  // Cost of one edge u->w as the searches price it; kInfCost if absent.
  Cost edge(size_t u, size_t w, bool u_source) const
  {
    int z, i, j, z2, i2, j2;
    g.decode(u, z, i, j);
    g.decode(w, z2, i2, j2);
    const TrackLayer& l = g.layers[z];
    if (!r.region->contains(g.layers[z2], z2, i2, j2)) {
      return kInfCost;
    }
    if (z == z2 && i == i2 && std::abs(j - j2) == 1) {
      size_t e = j2 > j ? u : w;
      if (!(u_source || enter(u, kWire)) || !enter(w, kWire)
          || !enter(e, kWireEdge)) {
        return kInfCost;
      }
      Cost rate = g.wire_unit;
      for (const DiscountRun& d : r.discounts) {
        if (d.z == z && d.i == i && std::min(j, j2) >= d.jlo
            && std::min(j, j2) < d.jhi) {
          rate = r.discount_rate;
        }
      }
      Cost spread = 0;
      if (r.spread_cost != 0 && r.spread_tiles != nullptr) {
        const int t = l.track_tile[i];
        const int p = l.pos_tile[j2];
        const int tx = l.horizontal ? p : t;
        const int ty = l.horizontal ? t : p;
        if ((*r.spread_tiles)[z][static_cast<size_t>(tx) * g.tiles_y + ty]) {
          for (const int di : {-1, 1}) {
            if (i + di >= 0 && i + di < l.numTracks()
                && f.level(g.vertex(z, i + di, j2), kWire, r.wire_type)
                       != kLegal) {
              spread = r.spread_cost;
            }
          }
        }
      }
      return std::abs(l.positions[j2] - l.positions[j]) * rate + pen(w, kWire)
             + pen(e, kWireEdge) + spread;
    }
    if (z == z2 && j == j2 && std::abs(i - i2) == 1) {
      size_t e = i2 > i ? u : w;
      if (!l.jogs || !(u_source || enter(u, kJog)) || !enter(w, kJog)
          || !enter(e, kJogEdge)) {
        return kInfCost;
      }
      return l.jog_factor * std::abs(l.tracks[i2] - l.tracks[i]) + pen(w, kJog)
             + pen(e, kJogEdge) + (u_source ? 0 : pen(u, kJog));
    }
    if (z2 == z + 1 && l.up_track[j] == i2 && l.up_pos[i] == j2) {
      if (!(enter(u, kViaUpPad) && enter(u, kCutUp) && enter(w, kViaDownPad))) {
        return kInfCost;
      }
      return g.via_cost[z] + pen(u, kViaUpPad) + pen(u, kCutUp)
             + pen(w, kViaDownPad);
    }
    if (z2 == z - 1 && l.down_track[j] == i2 && l.down_pos[i] == j2) {
      if (!(enter(w, kViaUpPad) && enter(w, kCutUp) && enter(u, kViaDownPad))) {
        return kInfCost;
      }
      return g.via_cost[z - 1] + pen(w, kViaUpPad) + pen(w, kCutUp)
             + pen(u, kViaDownPad);
    }
    return kInfCost;
  }
  // Recomputes a returned path's cost (source cost + edges + target cost).
  Cost pathCost(const SearchResult& res) const
  {
    if (res.path.empty()) {
      return kInfCost;
    }
    std::set<size_t> src;
    Cost c = kInfCost;
    for (auto [v, sc] : r.sources) {
      src.insert(v);
      if (v == res.path.back()) {
        c = std::min(c, sc);
      }
    }
    if (c == kInfCost) {
      return kInfCost;
    }
    for (size_t k = res.path.size() - 1; k > 0; k--) {
      Cost e = edge(res.path[k], res.path[k - 1], src.count(res.path[k]) > 0);
      if (e >= kInfCost) {
        return kInfCost;
      }
      c += e;
    }
    Cost tc = kInfCost;
    for (auto& t : r.targets) {
      if (t.vertex == res.path.front()) {
        tc = std::min(tc, t.cost);
      }
    }
    return tc >= kInfCost ? kInfCost : c + tc;
  }
  Cost run()
  {
    std::vector<Cost> d(g.numVertices(), kInfCost);
    std::set<size_t> src;
    using P = std::pair<Cost, size_t>;
    std::priority_queue<P, std::vector<P>, std::greater<>> q;
    for (auto [v, c] : r.sources) {
      src.insert(v);
      if (c < d[v]) {
        d[v] = c;
        q.push({c, v});
      }
    }
    Cost best = kInfCost;
    while (!q.empty()) {
      auto [dv, v] = q.top();
      q.pop();
      if (dv > d[v]) {
        continue;
      }
      for (auto& t : r.targets) {
        if (t.vertex == v) {
          best = std::min(best, dv + t.cost);
        }
      }
      int z, i, j;
      g.decode(v, z, i, j);
      const TrackLayer& l = g.layers[z];
      bool s = src.count(v);
      auto go = [&](size_t w, Cost c) {
        if (dv + c < d[w]) {
          d[w] = dv + c;
          q.push({d[w], w});
        }
      };
      if (s || enter(v, kWire)) {
        for (int dd : {-1, 1}) {
          int j2 = j + dd;
          if (j2 < 0 || j2 >= l.numPositions()
              || !r.region->contains(l, z, i, j2)) {
            continue;
          }
          size_t w = g.vertex(z, i, j2);
          size_t e = dd > 0 ? v : w;
          if (!enter(w, kWire) || !enter(e, kWireEdge)) {
            continue;
          }
          go(w, edge(v, w, s));
        }
      }
      if (l.jogs && (s || enter(v, kJog))) {
        for (int di : {-1, 1}) {
          int i2 = i + di;
          if (i2 < 0 || i2 >= l.numTracks()
              || !r.region->contains(l, z, i2, j)) {
            continue;
          }
          size_t w = g.vertex(z, i2, j);
          size_t e = di > 0 ? v : w;
          if (!enter(w, kJog) || !enter(e, kJogEdge)) {
            continue;
          }
          go(w,
             l.jog_factor * std::abs(l.tracks[i2] - l.tracks[i]) + pen(w, kJog)
                 + pen(e, kJogEdge) + (s ? 0 : pen(v, kJog)));
        }
      }
      if (z + 1 < g.numLayers()) {
        int i2 = l.up_track[j], j2 = l.up_pos[i];
        if (i2 >= 0 && j2 >= 0
            && r.region->contains(g.layers[z + 1], z + 1, i2, j2)) {
          size_t w = g.vertex(z + 1, i2, j2);
          if (enter(v, kViaUpPad) && enter(v, kCutUp)
              && enter(w, kViaDownPad)) {
            go(w,
               g.via_cost[z] + pen(v, kViaUpPad) + pen(v, kCutUp)
                   + pen(w, kViaDownPad));
          }
        }
      }
      if (z > 0) {
        int i2 = l.down_track[j], j2 = l.down_pos[i];
        if (i2 >= 0 && j2 >= 0
            && r.region->contains(g.layers[z - 1], z - 1, i2, j2)) {
          size_t w = g.vertex(z - 1, i2, j2);
          if (enter(w, kViaUpPad) && enter(w, kCutUp)
              && enter(v, kViaDownPad)) {
            go(w,
               g.via_cost[z - 1] + pen(w, kViaUpPad) + pen(w, kCutUp)
                   + pen(v, kViaDownPad));
          }
        }
      }
    }
    return best;
  }
};

// Path cost recomputed from the returned vertex sequence.
static bool validPath(const SearchResult& r, const SearchRequest& q)
{
  if (r.path.empty()) {
    return false;
  }
  size_t last = r.path.back();
  bool is_src = false;
  for (auto [v, c] : q.sources) {
    is_src |= v == last;
  }
  return is_src && r.path.front() == r.target_vertex;
}

// modern: wire unit 4 with discounted runs (arXiv:2111.06169 Sec. 5.2);
// only the vertex search supports it. potential: 0 = l1 target boxes,
// 1..3 = PotentialOracle l1 / simple / general on the request's Region.
static int randomCases(const int cases,
                       const bool sparse,
                       const bool modern,
                       const int potential = 0)
{
  std::mt19937 rng(sparse ? 23456 : 12345);
  int64_t pops_i = 0, pops_v = 0, fails = 0, found = 0, labels = 0, ivs = 0,
          settled_i = 0, settled_v = 0, fast_queries = 0, checker_queries = 0;
  for (int c = 0; c < cases; c++) {
    TrackGraph g;
    int Z = 2 + rng() % 3;
    const int span = 200 + rng() % 400;
    std::vector<std::vector<Coord>> tr(Z);
    for (int z = 0; z < Z; z++) {
      int pitch = 20 + rng() % 40;
      int off = rng() % pitch;
      for (int x = off; x < span;
           x += pitch + (rng() % 5 == 0 ? rng() % 30 : 0)) {
        tr[z].push_back(x);
      }
    }
    g.layers.resize(Z);
    for (int z = 0; z < Z; z++) {
      TrackLayer& l = g.layers[z];
      l.horizontal = (z % 2) == 0;
      l.tracks = tr[z];
      std::vector<Coord> pos;
      if (z > 0) {
        pos.insert(pos.end(), tr[z - 1].begin(), tr[z - 1].end());
      }
      if (z + 1 < Z) {
        pos.insert(pos.end(), tr[z + 1].begin(), tr[z + 1].end());
      }
      std::sort(pos.begin(), pos.end());
      pos.erase(std::unique(pos.begin(), pos.end()), pos.end());
      l.positions = pos;
      l.jogs = rng() % 3 != 0;
      l.jog_factor = 1 + rng() % 5;
    }
    if (modern) {
      g.wire_unit = 4;
      for (TrackLayer& l : g.layers) {
        l.jog_factor *= 4;
      }
    }
    for (int z = 0; z + 1 < Z; z++) {
      g.via_cost.push_back((10 + rng() % 200) * g.wire_unit);
    }
    const int tile = 100 + rng() % 300;
    auto tf = [&](Coord c) { return c / tile; };
    g.finalize(tf, tf);
    const int T = span / tile + 1;
    g.tiles_x = g.tiles_y = T;
    FastGrid f;
    f.init(&g, 2);
    TestOracle legality;
    legality.grid = &f;
    f.setOracle(&legality);
    double pfix = sparse ? 0.005 : (rng() % 30) / 100.0;
    double prout = sparse ? 0.005 : (rng() % 20) / 100.0;
    std::uniform_real_distribution<double> u(0, 1);
    for (size_t v = 0; v < g.numVertices(); v++) {
      for (int t = 0; t < 2; t++) {
        for (int k = 0; k <= kViaDownPad; k++) {
          // Non-removable shapes (blockages, pins, other nets'
          // reservations), routed shapes of ripup level 1..6.
          if (u(rng) < pfix) {
            f.setLevel(v, (Kind) k, t, kNotRemovable);
          } else if (u(rng) < prout) {
            f.setLevel(v, (Kind) k, t, 1 + rng() % kMaxRemovableLevel);
          } else if (u(rng) < (sparse ? 0.001 : 0.03)) {
            f.setLevel(v, (Kind) k, t, kNotRemovable);
          }
        }
      }
      // Edges whose swept wire meets shapes between the endpoints.
      for (const Kind e : {kWireEdge, kJogEdge}) {
        if (u(rng) < (sparse ? 0.002 : 0.05)) {
          f.setEdgeQuery(v, e, true);
          legality.edges[{v, e}]
              = u(rng) < 0.5 ? kNotRemovable : 1 + rng() % kMaxRemovableLevel;
        }
      }
      if (u(rng) < (sparse ? 0.0 : 0.05)) {
        f.bumpHistory(v);
      }
    }
    Region reg;
    reg.init(Z, T, T);
    bool full = rng() % 2;
    for (int z = 0; z < Z; z++) {
      for (int x = 0; x < T; x++) {
        for (int y = 0; y < T; y++) {
          if (full || u(rng) < 0.7) {
            reg.allow(z, x, y);
          }
        }
      }
    }
    SearchRequest q;
    q.net = rng() % 3;
    q.region = &reg;
    q.ripup = rng() % 2;
    q.ripup_cost = 50 + rng() % 500;
    q.max_ripup_level = 1 + rng() % kMaxRemovableLevel;
    // Wire type 2 is outside the fast grid: the checker answers.
    q.wire_type = rng() % 3;
    std::vector<std::vector<uint8_t>> spread_tiles;
    if (modern && rng() % 2) {
      spread_tiles.assign(Z, std::vector<uint8_t>(size_t(T) * T, 0));
      for (auto& layer_tiles : spread_tiles) {
        for (auto& tile_flag : layer_tiles) {
          tile_flag = rng() % 2;
        }
      }
      q.spread_tiles = &spread_tiles;
      q.spread_cost = 4 * (1 + rng() % 60);
    }
    if (modern) {
      q.discount_rate = 3;
      const int discounts = rng() % 6;
      for (int d = 0; d < discounts; d++) {
        const int z = rng() % Z;
        const TrackLayer& l = g.layers[z];
        const int i = rng() % l.numTracks();
        const int jlo = rng() % l.numPositions();
        const int jhi
            = std::min(l.numPositions() - 1, jlo + 1 + (int) (rng() % 40));
        q.discounts.push_back({z, i, jlo, jhi});
      }
    }
    int ns = 1 + rng() % 3, nt = 1 + rng() % 4;
    auto rv = [&]() { return (size_t) (rng() % g.numVertices()); };
    for (int s = 0; s < ns; s++) {
      q.sources.push_back({rv(), (Cost) (rng() % 3 == 0 ? rng() % 100 : 0)});
    }
    for (int t = 0; t < nt; t++) {
      size_t v = rv();
      int z, i, j;
      g.decode(v, z, i, j);
      auto [x, y] = g.layers[z].point(i, j);
      q.targets.push_back({v, (Cost) (rng() % 100), t});
      // Boxes may cover several targets; their cost is the cheapest one.
      if (t == 0 || rng() % 2) {
        q.target_boxes.push_back({x, y, x, y, z, z, q.targets.back().cost});
      } else {
        auto& b = q.target_boxes.back();
        b.cost = std::min(b.cost, q.targets.back().cost);
        b.xlo = std::min(b.xlo, x);
        b.xhi = std::max(b.xhi, x);
        b.ylo = std::min(b.ylo, y);
        b.yhi = std::max(b.yhi, y);
        b.zlo = std::min(b.zlo, z);
        b.zhi = std::max(b.zhi, z);
      }
    }
    PotentialOracle oracle;
    if (potential == 3) {
      // Contract of the general model (and of BonnRouteDR): sources and
      // targets lie in the routing area, so every search edge is an edge of
      // the grid graph whose distances the potential holds.
      for (const auto& [v, c] : q.sources) {
        int z, i, j;
        g.decode(v, z, i, j);
        const TrackLayer& l = g.layers[z];
        const int t = l.track_tile[i];
        const int p = l.pos_tile[j];
        reg.allow(z, l.horizontal ? p : t, l.horizontal ? t : p);
      }
    }
    if (potential > 0) {
      const GCellGrid gcells{0, 0, tile, tile, 0, 0, span, span};
      oracle.prepare(static_cast<PotentialOracle::Model>(potential - 1),
                     g,
                     q,
                     gcells,
                     3,
                     4);
      q.future_cost = &oracle;
    }
    PathSearch ps(g, f);
    IntervalSearch is(g, f);
    SearchResult ri = modern ? SearchResult{} : is.run(q);
    SearchResult rvx = ps.run(q);
    Ref ref{g, f, q};
    Cost rc = ref.run();
    Cost ci = ri.found ? ri.cost : kInfCost,
         cv = rvx.found ? rvx.cost : kInfCost;
    Cost pi = ri.found ? ref.pathCost(ri) : kInfCost;
    const Cost pv = rvx.found ? ref.pathCost(rvx) : kInfCost;
    if (modern) {
      ci = pi = rc;
    }
    if (ci != rc || cv != rc || pi != ci || pv != cv
        || (ri.found && !validPath(ri, q))
        || (rvx.found && !validPath(rvx, q))) {
      fails++;
      if (fails < 10) {
        printf(
            "case %d mismatch: interval %lld (path %lld) vertex %lld (path "
            "%lld) ref %lld\n",
            c,
            (long long) ci,
            (long long) pi,
            (long long) cv,
            (long long) pv,
            (long long) rc);
      }
    }
    found += modern ? rvx.found : ri.found;
    fast_queries += rvx.fast_grid_queries + ri.fast_grid_queries;
    checker_queries += rvx.checker_queries + ri.checker_queries;
    pops_i += ri.heap_pops;
    labels += is.labels();
    ivs += is.intervals();
    settled_i += ri.settled;
    settled_v += rvx.settled;
    pops_v += rvx.heap_pops;
  }
  printf(
      "cases %d found %lld fails %lld | heap pops interval %lld vertex %lld "
      "(%.2fx) | intervals %lld labels %lld | settled interval %lld vertex "
      "%lld | legality queries fast grid %lld checker %lld\n",
      cases,
      (long long) found,
      (long long) fails,
      (long long) pops_i,
      (long long) pops_v,
      pops_i ? double(pops_v) / pops_i : 0.0,
      (long long) ivs,
      (long long) labels,
      (long long) settled_i,
      (long long) settled_v,
      (long long) fast_queries,
      (long long) checker_queries);
  return fails != 0;
}

// Enumerate geometry directly, including shapes between legal endpoints,
// touching endpoints, and boxes outside the track extent in both directions.
static int blockedGeometry()
{
  const std::vector<std::pair<int, int>> bounds{
      {-20, -5}, {-20, 0}, {0, 0}, {1, 99}, {100, 100}, {101, 120}};
  for (const bool horizontal : {false, true}) {
    TrackLayer layer;
    layer.horizontal = horizontal;
    layer.tracks = {0, 100};
    layer.positions = {0, 100};
    for (const Kind kind : {kWire, kJog, kWireEdge, kJogEdge}) {
      for (const auto& [xlo, xhi] : bounds) {
        for (const auto& [ylo, yhi] : bounds) {
          std::set<std::pair<int, int>> actual;
          forEachBlocked(layer, kind, xlo, ylo, xhi, yhi, [&](int i, int j) {
            actual.emplace(i, j);
          });
          std::set<std::pair<int, int>> expected;
          for (int i = 0; i < layer.numTracks(); i++) {
            for (int j = 0; j < layer.numPositions(); j++) {
              const int ni = i + (kind == kJogEdge);
              const int nj = j + (kind == kWireEdge);
              if (ni >= layer.numTracks() || nj >= layer.numPositions()) {
                continue;
              }
              const auto [x0, y0] = layer.point(i, j);
              const auto [x1, y1] = layer.point(ni, nj);
              if (std::max(x0, x1) >= xlo && std::min(x0, x1) <= xhi
                  && std::max(y0, y1) >= ylo && std::min(y0, y1) <= yhi) {
                expected.emplace(i, j);
              }
            }
          }
          if (actual != expected) {
            std::printf(
                "blocked geometry mismatch: horizontal=%d kind=%d "
                "box=(%d,%d)-(%d,%d)\n",
                horizontal,
                kind,
                xlo,
                ylo,
                xhi,
                yhi);
            return 1;
          }
        }
      }
    }
  }
  return 0;
}

static int obstacleDetour()
{
  for (const bool horizontal : {false, true}) {
    TrackGraph graph;
    graph.layers.resize(1);
    TrackLayer& layer = graph.layers[0];
    layer.horizontal = horizontal;
    layer.tracks = {0, 100};
    layer.positions = {0, 100, 200, 300};
    graph.finalize([](int) { return 0; }, [](int) { return 0; });
    FastGrid grid;
    grid.init(&graph, 1);
    TestOracle oracle;
    oracle.grid = &grid;
    grid.setOracle(&oracle);
    // The obstacle hits the middle of edge (100, 200), while its endpoints
    // remain legal (Fig. 4 zigzag edge): the edge is marked for the
    // checker. The shortest path must jog to the other track and back.
    grid.setEdgeQuery(graph.vertex(0, 0, 1), kWireEdge, true);
    oracle.edges[{graph.vertex(0, 0, 1), kWireEdge}] = kNotRemovable;
    Region region;
    region.init(1, 1, 1);
    region.allow(0, 0, 0);
    SearchRequest request;
    request.net = 0;
    request.region = &region;
    request.sources = {{graph.vertex(0, 0, 0), 0}};
    request.targets = {{graph.vertex(0, 0, 3), 0, 0}};
    const auto [x, y] = layer.point(0, 3);
    request.target_boxes = {{x, y, x, y, 0, 0}};
    PathSearch vertex(graph, grid);
    IntervalSearch interval(graph, grid);
    Ref reference{graph, grid, request};
    for (const auto& result : {vertex.run(request), interval.run(request)}) {
      if (!result.found || result.cost != 1100
          || reference.pathCost(result) != result.cost) {
        std::printf("failed to detour around blocked edge: horizontal=%d\n",
                    horizontal);
        return 1;
      }
    }
  }
  return 0;
}

// Real workers reuse search objects as targets, windows and occupancy change.
// Exercise long track runs and multiple source labels against plain Dijkstra.
static int reusedSearch()
{
  for (const bool horizontal : {false, true}) {
    TrackGraph graph;
    graph.layers.resize(2);
    graph.via_cost = {90};
    for (int z = 0; z < 2; z++) {
      TrackLayer& layer = graph.layers[z];
      layer.horizontal = horizontal != (z == 1);
      layer.tracks = {0, 100, 200, 300, 400};
      for (int j = 0; j <= 40; j++) {
        layer.positions.push_back(j * 10);
      }
    }
    auto tile = [](Coord c) { return c / 100; };
    graph.finalize(tile, tile);
    FastGrid grid;
    grid.init(&graph, 1);
    TestOracle oracle;
    oracle.grid = &grid;
    grid.setOracle(&oracle);
    Region region;
    region.init(2, 5, 5);
    PathSearch vertex(graph, grid);
    IntervalSearch interval(graph, grid);
    std::mt19937 rng(287);
    for (int query = 0; query < 48; query++) {
      region.reset();
      const int first = query % 3 == 0 ? 1 : 0;
      const int last = query % 3 == 1 ? 3 : 4;
      for (int z = 0; z < 2; z++) {
        for (int x = first; x <= last; x++) {
          for (int y = first; y <= last; y++) {
            region.allow(z, x, y);
          }
        }
      }
      const size_t obstacle = graph.vertex(0, 2, 20);
      const int occupied = query % 2;
      grid.setEdgeQuery(obstacle, kWireEdge, occupied != 0);
      oracle.edges[{obstacle, kWireEdge}] = 1;
      SearchRequest request;
      request.net = 0;
      request.region = &region;
      request.ripup = query % 4 == 1;
      request.ripup_cost = 80;
      for (int s = 0; s < 4; s++) {
        request.sources.emplace_back(
            graph.vertex(0, 2, first * 10 + rng() % ((last - first) * 10 + 1)),
            rng() % 200);
      }
      for (int t = 0; t < 3; t++) {
        const size_t v
            = graph.vertex(t % 2, first + t, (query % 2 ? first : last) * 10);
        request.targets.push_back({v, static_cast<Cost>(rng() % 80), t});
        const auto [x, y] = graph.point(v);
        request.target_boxes.push_back({x, y, x, y, t % 2, t % 2});
      }
      if (query % 5 == 0) {
        request.target_boxes.clear();
      }
      Ref reference{graph, grid, request};
      const Cost expected = reference.run();
      for (const auto& result : {vertex.run(request), interval.run(request)}) {
        const Cost actual = result.found ? result.cost : kInfCost;
        if (actual != expected
            || (result.found && reference.pathCost(result) != actual)) {
          std::printf("reused search mismatch: horizontal=%d query=%d\n",
                      horizontal,
                      query);
          return 1;
        }
      }
      grid.setEdgeQuery(obstacle, kWireEdge, false);
    }
  }
  return 0;
}

// Reservations (pin access paths, input wires) are shape-grid shapes owned
// by their net (Sec. 4.3, 4.4): against a multiset model, a model placed on
// a reservation is blocked for a net iff a live reservation of another
// owner violates it (the querying net's own shapes count as removed).
static int reservationCases()
{
  std::mt19937 rng(7);
  int fails = 0;
  for (int c = 0; c < 300; c++) {
    ShapeClassTable classes;
    LayerRules rules;
    rules.spacing = SpacingTable::constant(4 + rng() % 4);
    classes.setLayers({rules});
    ShapeGrid grid;
    ShapeGridLayerGeometry g;
    g.along_cell = g.cross_cell
        = std::max<Coord>(1, maximalCellSize(classes, 0));
    g.along_cells = g.cross_cells = 120 / g.along_cell + 1;
    grid.init(&classes, {g}, ChunkAxis{}, ChunkAxis{});
    DistanceRuleChecker checker(&grid);
    ShapeGridWriteCache cache;
    const int owners = 1 + rng() % 4;
    // A few reservation footprints, shared by several owners.
    std::vector<Rect> spots;
    for (int k = 0; k < 4; k++) {
      const Coord x = 10 + rng() % 80;
      const Coord y = 10 + rng() % 80;
      spots.push_back({x, y, x + 3 + (Coord) (rng() % 20), y + 6});
    }
    struct Live
    {
      Rect rect;
      ClassId cls;
      int owner;
    };
    std::vector<Live> live;
    for (int step = 0; step < 120; step++) {
      if (live.empty() || rng() % 2) {
        const Rect r = spots[rng() % spots.size()];
        const Live l{r,
                     classes.classify(0, r, ShapeRole::kFixed),
                     (int) (rng() % owners)};
        grid.insert(0, l.rect, l.cls, ShapeType::kReservation, l.owner, cache);
        live.push_back(l);
      } else {
        const size_t at = rng() % live.size();
        const Live l = live[at];
        live.erase(live.begin() + at);
        if (!grid.remove(
                0, l.rect, l.cls, ShapeType::kReservation, l.owner, cache)) {
          fails++;
        }
      }
      const Rect model{-2, -2, 2, 2};
      const LocatedModel m{
          0, model, classes.classify(0, model, ShapeRole::kFixed), 3, 0};
      for (int q = 0; q < 4; q++) {
        const Point at{(Coord) (rng() % 110), (Coord) (rng() % 110)};
        for (int net = -1; net <= owners; net++) {
          bool expect = false;
          for (const Live& l : live) {
            expect
                |= l.owner != net
                   && shapesViolate(classes,
                                    PlacedShape{model.translated(at), m.cls, 3},
                                    PlacedShape{l.rect, l.cls, 3});
          }
          CheckOptions opts;
          opts.ignored_owner = net;
          opts.intervals = false;
          CheckResult res;
          checker.check(at, std::span(&m, 1), opts, res);
          if ((res.levels[0] == kNotRemovable) != expect
              || (res.levels[0] != kLegal && res.levels[0] != kNotRemovable)) {
            if (fails++ < 5) {
              printf("reservation mismatch: case %d step %d\n", c, step);
            }
          }
        }
      }
    }
  }
  return fails != 0;
}

// A search whose label budget is exceeded stops without a path and says
// so; the same search without a budget finds the path.
static int labelBudget()
{
  TrackGraph graph;
  graph.layers.resize(1);
  TrackLayer& layer = graph.layers[0];
  layer.tracks = {0, 100, 200};
  for (int j = 0; j < 4000; j++) {
    layer.positions.push_back(j * 10);
  }
  graph.finalize([](int) { return 0; }, [](int) { return 0; });
  FastGrid grid;
  grid.init(&graph, 1);
  Region region;
  region.init(1, 1, 1);
  region.allow(0, 0, 0);
  SearchRequest request;
  request.region = &region;
  request.sources = {{graph.vertex(0, 0, 0), 0}};
  request.targets = {{graph.vertex(0, 2, 3999), 0, 0}};
  PathSearch search(graph, grid);
  request.label_limit = 1500;
  const SearchResult capped = search.run(request);
  request.label_limit = 0;
  const SearchResult full = search.run(request);
  if (!capped.aborted || capped.found || full.aborted || !full.found) {
    std::printf("label budget: capped aborted=%d found=%d, full found=%d\n",
                capped.aborted,
                capped.found,
                full.found);
    return 1;
  }
  return 0;
}

int main()
{
  if (blockedGeometry() || obstacleDetour() || reusedSearch()
      || reservationCases() || labelBudget()) {
    return 1;
  }
  return randomCases(400, false, false) || randomCases(400, true, false)
         || randomCases(400, false, true) || randomCases(400, true, true)
         || randomCases(300, false, true, 1) || randomCases(300, true, true, 2)
         || randomCases(300, false, true, 3) || randomCases(300, true, true, 3)
         || randomCases(300, false, false, 3);
}
