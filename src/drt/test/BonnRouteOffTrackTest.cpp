// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Randomized tests of the off-track modules (Gester et al., TODAES 18(2),
// 2013, Sec. 3.7, 3.8, 4.3, 4.4):
//  - Algorithm 3 against a literal transcription;
//  - Theorem 3.2: shortest tau-feasible path on the path-preserving digraph
//    == brute-force shortest tau-feasible path on the integer grid (one and
//    two layers), and every returned path is tau-feasible and obstacle-free;
//  - the same-net mapping tau = min_step + w against the boundary polygon;
//  - destructive bounding against exhaustive enumeration;
//  - circuit signatures under the eight isometries;
//  - the same-net post-processing invariants.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "br/BrBlockageGrid.h"
#include "br/BrGeometry.h"
#include "br/BrOffTrackSearch.h"
#include "br/BrPathCleanup.h"
#include "br/BrPinAccess.h"
#include "br/BrSameNetRules.h"

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

constexpr int64_t kInf = std::numeric_limits<int64_t>::max() / 4;

// ---------------------------------------------------------------------------
// Algorithm 3, transcribed line by line (sentinels -inf, inf as int64).
std::vector<GeoCoord> algorithm3Literal(const std::vector<GeoCoord>& borders,
                                        const GeoCoord tau)
{
  const int64_t ninf = std::numeric_limits<int64_t>::min() / 4;
  const int64_t pinf = std::numeric_limits<int64_t>::max() / 4;
  std::vector<int64_t> vbar{ninf, pinf};  // line 1 (s, t are in borders)
  for (const GeoCoord x : borders) {      // line 2
    vbar.push_back(x);
  }
  std::sort(vbar.begin(), vbar.end());  // line 3
  vbar.erase(std::unique(vbar.begin(), vbar.end()), vbar.end());
  std::vector<int64_t> v;
  const int nv = static_cast<int>(vbar.size());
  for (int i = 1; i + 1 < nv; i++) {  // line 4 (finite x_i only)
    int k = 0;                        // line 5: min k in N (0 included)
    while (!(vbar[i - k] - vbar[i - k - 1] >= 4 * int64_t{tau})) {
      k++;
    }
    const int cmin = i - k;
    k = 0;  // line 6
    while (!(vbar[i + k + 1] - vbar[i + k] >= 4 * int64_t{tau})) {
      k++;
    }
    const int cmax = i + k;
    // line 7: lambda in Z with x_i + lambda tau in [x_cmin - 2tau,
    // x_cmax + 2tau]; enumerate by stepping from x_i.
    for (int64_t x = vbar[i]; x >= vbar[cmin] - 2 * tau; x -= tau) {
      v.push_back(x);  // line 8
    }
    for (int64_t x = vbar[i] + tau; x <= vbar[cmax] + 2 * tau; x += tau) {
      v.push_back(x);
    }
  }
  std::sort(v.begin(), v.end());  // line 11
  v.erase(std::unique(v.begin(), v.end()), v.end());
  return {v.begin(), v.end()};
}

void testAlgorithm3(std::mt19937& rng)
{
  int cases = 0;
  double max_ratio = 0;
  int64_t max_lines = 0;
  for (int it = 0; it < 4000; it++) {
    const int n = 1 + static_cast<int>(rng() % 30);
    const GeoCoord tau = 1 + static_cast<GeoCoord>(rng() % 12);
    const int spread = 1 + static_cast<int>(rng() % 400);
    std::vector<GeoCoord> b;
    for (int k = 0; k < n; k++) {
      b.push_back(static_cast<GeoCoord>(rng() % spread) - spread / 2);
    }
    const std::vector<GeoCoord> got = blockageGridCoordinates(b, tau);
    const std::vector<GeoCoord> want = algorithm3Literal(b, tau);
    CHECK(got == want,
          "n=%d tau=%d size %zu vs %zu",
          n,
          tau,
          got.size(),
          want.size());
    for (const GeoCoord x : b) {
      CHECK(std::binary_search(got.begin(), got.end(), x), "hanan line lost");
    }
    std::vector<GeoCoord> u = b;
    std::sort(u.begin(), u.end());
    u.erase(std::unique(u.begin(), u.end()), u.end());
    const int64_t m = static_cast<int64_t>(u.size());
    // Each line adds at most 4 (run length - 1) + 5 coordinates.
    CHECK(static_cast<int64_t>(got.size()) <= 4 * m * m + m, "bound");
    max_ratio = std::max(max_ratio, static_cast<double>(got.size()) / m);
    max_lines = std::max<int64_t>(max_lines, got.size());
    cases++;
  }
  std::printf(
      "algorithm3: %d cases equal to the literal transcription; max |V|/n "
      "= %.2f, max |V| = %lld (bound 4n^2+n)\n",
      cases,
      max_ratio,
      static_cast<long long>(max_lines));
}

// ---------------------------------------------------------------------------
// Theorem 3.2: brute force on the integer grid.
struct Instance
{
  struct Layer
  {
    GeoCoord tau = 1;
    bool wires = true;
    std::vector<Box> obstacles;
    // Cost per unit of x / y wire (axis weights: the jog factor of a
    // non-preferred direction).
    int64_t cost_x = 1;
    int64_t cost_y = 1;
  };
  std::vector<Layer> layers;
  std::vector<Box> via_obstacles;  // between layer 0 and 1
  // If not empty: vias only at these points (pin access on a layer without
  // wires: only at the access points).
  std::vector<GridPoint> via_points;
  int64_t via_cost = 0;
  RouteNode s;
  RouteNode t;
};

bool edgeBlocked(const std::vector<Box>& obs,
                 const GridPoint& a,
                 const GridPoint& b)
{
  for (const Box& o : obs) {
    if (segmentMeetsInterior(a, b, o)) {
      return true;
    }
  }
  return false;
}

bool viaBlocked(const Instance& in, const GridPoint& p)
{
  if (!in.via_points.empty()
      && std::find(in.via_points.begin(), in.via_points.end(), p)
             == in.via_points.end()) {
    return true;
  }
  for (const Box& o : in.via_obstacles) {
    if (o.containsOpen(p.x, p.y)) {
      return true;
    }
  }
  return false;
}

// Dijkstra over (layer, x, y, direction, run length capped at tau).
int64_t bruteForce(const Instance& in, const Box& dom)
{
  const int nl = static_cast<int>(in.layers.size());
  const int w = dom.xhi - dom.xlo + 1;
  const int h = dom.yhi - dom.ylo + 1;
  int tmax = 0;
  for (const auto& l : in.layers) {
    tmax = std::max(tmax, l.tau);
  }
  const int nlen = tmax + 1;
  auto id = [&](int z, int x, int y, int d, int len) {
    return ((((static_cast<int64_t>(z) * w + (x - dom.xlo)) * h + (y - dom.ylo))
                 * 5
             + d)
            * nlen)
           + len;
  };
  std::vector<int64_t> dist(static_cast<size_t>(nl) * w * h * 5 * nlen, kInf);
  using Item = std::pair<int64_t, int64_t>;
  std::priority_queue<Item, std::vector<Item>, std::greater<>> pq;
  const int64_t s0 = id(in.s.layer, in.s.p.x, in.s.p.y, 4, 0);
  dist[s0] = 0;
  pq.emplace(0, s0);
  const int dx[4] = {1, -1, 0, 0};
  const int dy[4] = {0, 0, 1, -1};
  while (!pq.empty()) {
    auto [c, st] = pq.top();
    pq.pop();
    if (c > dist[st]) {
      continue;
    }
    int64_t r = st;
    const int len = static_cast<int>(r % nlen);
    r /= nlen;
    const int d = static_cast<int>(r % 5);
    r /= 5;
    const int y = static_cast<int>(r % h) + dom.ylo;
    r /= h;
    const int x = static_cast<int>(r % w) + dom.xlo;
    const int z = static_cast<int>(r / w);
    const Instance::Layer& l = in.layers[z];
    const bool done = d == 4 || len >= l.tau;
    if (z == in.t.layer && x == in.t.p.x && y == in.t.p.y && done) {
      return c;
    }
    auto push = [&](int64_t to, int64_t nc) {
      if (nc < dist[to]) {
        dist[to] = nc;
        pq.emplace(nc, to);
      }
    };
    if (l.wires) {
      for (int nd = 0; nd < 4; nd++) {
        const bool straight = nd == d;
        const bool perpendicular = d == 4 || (nd / 2 != d / 2);
        if (!straight && !(perpendicular && done)) {
          continue;
        }
        const int nx = x + dx[nd];
        const int ny = y + dy[nd];
        if (nx < dom.xlo || nx > dom.xhi || ny < dom.ylo || ny > dom.yhi) {
          continue;
        }
        if (edgeBlocked(l.obstacles, {x, y}, {nx, ny})) {
          continue;
        }
        const int nlen2 = straight ? std::min(len + 1, tmax) : 1;
        push(id(z, nx, ny, nd, std::min(nlen2, tmax)),
             c + (nd < 2 ? l.cost_x : l.cost_y));
      }
    }
    if (done) {
      for (const int nz : {z - 1, z + 1}) {
        if (nz < 0 || nz >= nl || viaBlocked(in, {x, y})) {
          continue;
        }
        push(id(nz, x, y, 4, 0), c + in.via_cost);
      }
    }
  }
  return kInf;
}

OffTrackRequest makeRequest(const Instance& in,
                            std::vector<BlockageGridLayer>& grid)
{
  BlockageGridInput gin;
  for (const auto& l : in.layers) {
    // points_only is valid where vias into the layer exist only at points.
    gin.layers.push_back(
        {l.tau, l.obstacles, !l.wires && !in.via_points.empty()});
  }
  if (in.layers.size() > 1) {
    gin.via_blockages.push_back(in.via_obstacles);
  }
  gin.points = {in.s.p, in.t.p};
  grid = buildBlockageGrid(gin);
  OffTrackRequest req;
  req.window = {-1000000, -1000000, 1000000, 1000000};
  for (size_t z = 0; z < in.layers.size(); z++) {
    OffTrackLayer l;
    l.grid = grid[z];
    l.obstacles = in.layers[z].obstacles;
    l.wires = in.layers[z].wires;
    l.cost_x = in.layers[z].cost_x;
    l.cost_y = in.layers[z].cost_y;
    req.layers.push_back(std::move(l));
  }
  req.sources.push_back({in.s.layer, {in.s.p.x, in.s.p.y, in.s.p.x, in.s.p.y}});
  req.targets.push_back({in.t.layer, in.t.p, 0});
  req.via_cost = in.via_cost;
  req.via_allowed
      = [&in](int, const GridPoint& p) { return !viaBlocked(in, p); };
  return req;
}

// Checks a returned path: endpoints, axis-parallel segments of length >=
// tau, obstacle-free, legal vias; returns its cost or -1.
int64_t verifyPath(const Instance& in, const OffTrackPath& path)
{
  const auto& n = path.nodes;
  if (n.empty() || !(n.front() == in.s) || !(n.back() == in.t)) {
    return -1;
  }
  int64_t cost = 0;
  for (size_t k = 0; k + 1 < n.size(); k++) {
    const RouteNode& a = n[k];
    const RouteNode& b = n[k + 1];
    if (a.layer != b.layer) {
      if (std::abs(a.layer - b.layer) != 1 || !(a.p == b.p)
          || viaBlocked(in, a.p)) {
        return -1;
      }
      cost += in.via_cost;
      continue;
    }
    if (a.p.x != b.p.x && a.p.y != b.p.y) {
      return -1;
    }
    const Instance::Layer& l = in.layers[a.layer];
    const int64_t len
        = std::abs(int64_t{a.p.x} - b.p.x) + std::abs(int64_t{a.p.y} - b.p.y);
    if (len < l.tau || !l.wires || edgeBlocked(l.obstacles, a.p, b.p)) {
      return -1;
    }
    cost += len * (a.p.y == b.p.y ? l.cost_x : l.cost_y);
  }
  return cost;
}

GridPoint freePoint(std::mt19937& rng,
                    const std::vector<Box>& obs,
                    const int size)
{
  for (int tries = 0;; tries++) {
    const GridPoint p{static_cast<GeoCoord>(rng() % size),
                      static_cast<GeoCoord>(rng() % size)};
    bool inside = false;
    for (const Box& o : obs) {
      inside = inside || o.containsOpen(p.x, p.y);
    }
    if (!inside || tries > 100) {
      return p;
    }
  }
}

std::vector<Box> randomBoxes(std::mt19937& rng,
                             int count,
                             int size,
                             int max_side)
{
  std::vector<Box> out;
  for (int k = 0; k < count; k++) {
    const GeoCoord x = static_cast<GeoCoord>(rng() % size);
    const GeoCoord y = static_cast<GeoCoord>(rng() % size);
    const GeoCoord wx = 1 + static_cast<GeoCoord>(rng() % max_side);
    const GeoCoord wy = 1 + static_cast<GeoCoord>(rng() % max_side);
    out.push_back({x, y, x + wx, y + wy});
  }
  return out;
}

void testTheorem(std::mt19937& rng, const int layers, const int count)
{
  static const GeoCoord kTaus[] = {1, 2, 3, 4, 5, 7};
  int equal = 0;
  int unreachable = 0;
  int bad_paths = 0;
  int64_t max_vertices = 0;
  std::map<int, std::pair<int64_t, int>> vertices_by_n;  // n -> (sum, cases)
  OffTrackSearch search;
  int64_t astar_settled = 0;
  int64_t dijkstra_settled = 0;
  int points_only_cases = 0;
  int weighted_cases = 0;
  for (int it = 0; it < count; it++) {
    Instance in;
    const int size = 20 + static_cast<int>(rng() % 20);
    int n = 0;
    for (int z = 0; z < layers; z++) {
      Instance::Layer l;
      l.tau = kTaus[rng() % 6];
      const int k = static_cast<int>(rng() % 9);
      l.obstacles = randomBoxes(rng, k, size, 10);
      n += k;
      l.wires = layers == 1 || z == 1 || rng() % 3 != 0;
      if (rng() % 2 == 0) {
        // Preferred direction 1, the other the jog factor.
        const int64_t jog = 1 + static_cast<int64_t>(rng() % 4);
        (rng() % 2 == 0 ? l.cost_x : l.cost_y) = jog;
        weighted_cases++;
      }
      in.layers.push_back(l);
    }
    if (layers > 1) {
      const int k = static_cast<int>(rng() % 4);
      in.via_obstacles = randomBoxes(rng, k, size, 8);
      n += k;
      in.via_cost = static_cast<int64_t>(rng() % 6);
    }
    in.s.layer = static_cast<int>(rng() % layers);
    in.t.layer = static_cast<int>(rng() % layers);
    in.s.p = freePoint(rng, in.layers[in.s.layer].obstacles, size);
    in.t.p = freePoint(rng, in.layers[in.t.layer].obstacles, size);
    if (layers > 1 && !in.layers[0].wires && rng() % 2 == 0) {
      in.via_points = {in.s.p, in.t.p};
      points_only_cases++;
    }
    std::vector<BlockageGridLayer> grid;
    OffTrackRequest req = makeRequest(in, grid);
    const OffTrackResult res = search.run(req);
    // The future cost changes no distance: plain Dijkstra agrees, also on
    // the nearest of several targets with max_targets.
    {
      OffTrackRequest multi = req;
      for (int k = 1; k < 6; k++) {
        const int layer = static_cast<int>(rng() % layers);
        multi.targets.push_back(
            {layer, freePoint(rng, in.layers[layer].obstacles, size), k});
      }
      multi.future_cost = false;
      const OffTrackResult all = search.run(multi);
      multi.future_cost = true;
      multi.max_targets = 3;
      const OffTrackResult nearest = search.run(multi);
      std::vector<int64_t> want;
      std::vector<int64_t> got_costs;
      for (const OffTrackPath& p : all.paths) {
        want.push_back(p.cost);
      }
      for (const OffTrackPath& p : nearest.paths) {
        got_costs.push_back(p.cost);
      }
      want.resize(std::min(want.size(), got_costs.size()));
      CHECK(got_costs.size() >= std::min<size_t>(3, all.paths.size()),
            "A* found %zu of %zu",
            got_costs.size(),
            all.paths.size());
      CHECK(std::equal(want.begin(), want.end(), got_costs.begin()),
            "A* nearest targets differ from Dijkstra");
      astar_settled += nearest.settled;
      dijkstra_settled += all.settled;
    }
    int tmax = 1;
    for (const auto& l : in.layers) {
      tmax = std::max(tmax, l.tau);
    }
    const GeoCoord margin = 4 * tmax + 4;
    const Box dom{-margin, -margin, size + 12 + margin, size + 12 + margin};
    const int64_t brute = bruteForce(in, dom);
    const int64_t got = res.paths.empty() ? kInf : res.paths.front().cost;
    if (!res.paths.empty()) {
      const int64_t verified = verifyPath(in, res.paths.front());
      if (verified != got) {
        bad_paths++;
      }
    }
    CHECK(got == brute,
          "layers=%d it=%d G=%lld brute=%lld",
          layers,
          it,
          static_cast<long long>(got == kInf ? -1 : got),
          static_cast<long long>(brute == kInf ? -1 : brute));
    equal += got == brute;
    unreachable += brute == kInf;
    max_vertices = std::max(max_vertices, res.grid_vertices);
    auto& acc = vertices_by_n[n];
    acc.first += res.grid_vertices;
    acc.second++;
  }
  CHECK(bad_paths == 0, "%d returned paths not tau-feasible", bad_paths);
  std::printf(
      "theorem 3.2 (%d layer%s): %d/%d instances G == brute force (%d without "
      "path), %d infeasible returned paths; max grid vertices %lld; %d "
      "with a points_only first layer, %d layers with axis weights\n",
      layers,
      layers > 1 ? "s" : "",
      equal,
      count,
      unreachable,
      bad_paths,
      static_cast<long long>(max_vertices),
      points_only_cases,
      weighted_cases);
  std::printf(
      "  A* (3 nearest of 6 targets) settled %.1f%% of the states plain "
      "Dijkstra settles for all targets, same costs\n",
      100.0 * astar_settled / std::max<int64_t>(1, dijkstra_settled));
  std::printf("  mean grid vertices by obstacle count n:");
  for (const auto& [n, acc] : vertices_by_n) {
    std::printf(" n=%d:%.0f", n, double(acc.first) / acc.second);
  }
  std::printf("\n");
}

// Worst-case growth: n obstacles packed with gaps < 4 tau.
void reportVertexGrowth()
{
  std::printf("vertex growth (tau = 10, single layer, 2 points):\n");
  std::printf("  %6s %14s %14s %14s %14s\n",
              "n",
              "packed |X||Y|",
              "packed/n^2",
              "sparse |X||Y|",
              "sparse/n^2");
  for (const int n : {1, 2, 4, 8, 16, 32, 64}) {
    BlockageGridInput packed;
    BlockageGridInput sparse;
    packed.layers.push_back({10, {}});
    sparse.layers.push_back({10, {}});
    for (int k = 0; k < n; k++) {
      // Packed: borders 7 apart (< 4 tau); sparse: 100 apart (>= 4 tau).
      packed.layers[0].blockages.push_back(
          {14 * k, 14 * k, 14 * k + 7, 14 * k + 7});
      sparse.layers[0].blockages.push_back(
          {200 * k, 200 * k, 200 * k + 100, 200 * k + 100});
    }
    packed.points = {{0, 0}, {1, 1}};
    sparse.points = {{0, 0}, {1, 1}};
    const auto gp = buildBlockageGrid(packed);
    const auto gs = buildBlockageGrid(sparse);
    const double nn = double(n) * n;
    std::printf("  %6d %14zu %14.2f %14zu %14.2f\n",
                n,
                gp[0].numVertices(),
                gp[0].numVertices() / nn,
                gs[0].numVertices(),
                gs[0].numVertices() / nn);
  }
}

// ---------------------------------------------------------------------------
// B1: tau = min_step + w makes every pair of adjacent boundary edges of a
// tau-feasible simple path contain an edge >= min_step.
struct BoundaryCheck
{
  bool simple = true;
  int violations = 0;  // adjacent edge pairs both shorter than min step
  int edges = 0;
};

BoundaryCheck boundaryEdges(const std::vector<GridPoint>& path,
                            const GeoCoord half,
                            const GeoCoord min_step)
{
  BoundaryCheck out;
  std::vector<Box> boxes;
  for (size_t k = 0; k + 1 < path.size(); k++) {
    boxes.push_back(wireBox(path[k], path[k + 1], half));
  }
  Box bb = boxes.front();
  for (const Box& b : boxes) {
    bb = bb.merged(b);
  }
  const int w = bb.xhi - bb.xlo + 2;
  const int h = bb.yhi - bb.ylo + 2;
  auto filled = [&](int cx, int cy) {
    // Cell [cx, cx + 1) x [cy, cy + 1) in absolute coordinates.
    for (const Box& b : boxes) {
      if (cx >= b.xlo && cx + 1 <= b.xhi && cy >= b.ylo && cy + 1 <= b.yhi) {
        return true;
      }
    }
    return false;
  };
  std::vector<char> cell(static_cast<size_t>(w + 2) * (h + 2), 0);
  auto at = [&](int cx, int cy) -> char& {
    return cell[static_cast<size_t>(cx - bb.xlo + 1) * (h + 2)
                + (cy - bb.ylo + 1)];
  };
  for (int cx = bb.xlo - 1; cx <= bb.xhi; cx++) {
    for (int cy = bb.ylo - 1; cy <= bb.yhi; cy++) {
      at(cx, cy) = filled(cx, cy) ? 1 : 0;
    }
  }
  // Directed unit edges with the filled cell on the left.
  std::map<GridPoint, std::vector<GridPoint>> next;
  for (int cx = bb.xlo; cx < bb.xhi; cx++) {
    for (int cy = bb.ylo; cy < bb.yhi; cy++) {
      if (!at(cx, cy)) {
        continue;
      }
      if (!at(cx, cy - 1)) {
        next[{cx, cy}].push_back({cx + 1, cy});
      }
      if (!at(cx, cy + 1)) {
        next[{cx + 1, cy + 1}].push_back({cx, cy + 1});
      }
      if (!at(cx - 1, cy)) {
        next[{cx, cy + 1}].push_back({cx, cy});
      }
      if (!at(cx + 1, cy)) {
        next[{cx + 1, cy}].push_back({cx + 1, cy + 1});
      }
    }
  }
  for (const auto& [p, outs] : next) {
    if (outs.size() != 1) {
      out.simple = false;
      return out;
    }
  }
  // One cycle (the union of a simple path is simply connected).
  std::vector<GridPoint> cycle;
  GridPoint p = next.begin()->first;
  const GridPoint start = p;
  do {
    cycle.push_back(p);
    p = next[p].front();
  } while (!(p == start) && cycle.size() <= next.size());
  if (cycle.size() != next.size()) {
    out.simple = false;  // several boundary components
    return out;
  }
  // Corners and edge lengths.
  std::vector<int64_t> lengths;
  const size_t n = cycle.size();
  std::vector<size_t> corners;
  for (size_t k = 0; k < n; k++) {
    const GridPoint& a = cycle[(k + n - 1) % n];
    const GridPoint& b = cycle[k];
    const GridPoint& c = cycle[(k + 1) % n];
    if (!((a.x == b.x && b.x == c.x) || (a.y == b.y && b.y == c.y))) {
      corners.push_back(k);
    }
  }
  for (size_t k = 0; k < corners.size(); k++) {
    const size_t a = corners[k];
    const size_t b = corners[(k + 1) % corners.size()];
    lengths.push_back(static_cast<int64_t>((b + n - a) % n));
  }
  out.edges = static_cast<int>(lengths.size());
  for (size_t k = 0; k < lengths.size(); k++) {
    const int64_t l1 = lengths[k];
    const int64_t l2 = lengths[(k + 1) % lengths.size()];
    if (l1 < min_step && l2 < min_step) {
      out.violations++;
    }
  }
  return out;
}

// Random rectilinear paths with segment lengths in [tau, 3 tau] whose
// non-adjacent segments keep a gap (so the metal is a simple polygon).
std::vector<GridPoint> randomPath(std::mt19937& rng,
                                  const GeoCoord tau,
                                  const GeoCoord lo_len,
                                  const GeoCoord half,
                                  const int segments)
{
  std::vector<GridPoint> p{{0, 0}};
  bool horizontal = rng() % 2 == 0;
  for (int k = 0; k < segments; k++) {
    for (int tries = 0; tries < 20; tries++) {
      const GeoCoord len
          = lo_len + static_cast<GeoCoord>(rng() % (2 * tau + 1));
      const GeoCoord sign = rng() % 2 == 0 ? 1 : -1;
      const GridPoint q = horizontal
                              ? GridPoint{p.back().x + sign * len, p.back().y}
                              : GridPoint{p.back().x, p.back().y + sign * len};
      const Box nb = wireBox(p.back(), q, half);
      bool ok = true;
      for (size_t j = 0; j + 2 < p.size() && ok; j++) {
        // Gap of at least 1 to non-adjacent segments.
        ok = !wireBox(p[j], p[j + 1], half).bloated(1).touches(nb);
      }
      if (ok) {
        p.push_back(q);
        horizontal = !horizontal;
        break;
      }
    }
  }
  return p;
}

void testSameNetMapping(std::mt19937& rng)
{
  int checked = 0;
  int violations_mapped = 0;
  int violations_width = 0;
  int checked_width = 0;
  for (int it = 0; it < 3000; it++) {
    const GeoCoord half = 1 + static_cast<GeoCoord>(rng() % 3);
    const GeoCoord w = 2 * half;
    const GeoCoord min_step = 1 + static_cast<GeoCoord>(rng() % 8);
    SameNetLayerRules rules;
    rules.width = w;
    rules.min_step_length = min_step;
    const GeoCoord tau = minimumSegmentLength(rules, false, 0, w).tau;
    CHECK(tau == min_step + w || (min_step + w < w && tau == w), "tau");
    const auto path = randomPath(rng, tau, tau, half, 2 + rng() % 6);
    if (path.size() < 2) {
      continue;
    }
    const BoundaryCheck b = boundaryEdges(path, half, min_step);
    if (b.simple) {
      checked++;
      violations_mapped += b.violations;
    }
    // Teeth: segments as short as the width only.
    const auto loose = randomPath(rng, tau, w, half, 2 + rng() % 6);
    if (loose.size() >= 2) {
      const BoundaryCheck c = boundaryEdges(loose, half, min_step);
      if (c.simple) {
        checked_width++;
        violations_width += c.violations > 0;
      }
    }
  }
  CHECK(violations_mapped == 0,
        "%d short-edge pairs with tau = step + w",
        violations_mapped);
  std::printf(
      "same-net mapping: %d simple paths with tau = min_step + w, %d "
      "adjacent short-edge pairs; control (segments >= w): %d/%d paths "
      "violate\n",
      checked,
      violations_mapped,
      violations_width,
      checked_width);
}

// ---------------------------------------------------------------------------
// Destructive bounding vs exhaustive enumeration.
void testSelection(std::mt19937& rng)
{
  int cases = 0;
  int greedy_worse = 0;
  int64_t nodes = 0;
  int64_t max_nodes = 0;
  int iterations = 0;
  double greedy_excess = 0;
  for (int it = 0; it < 3000; it++) {
    SelectionInput in;
    const int pins = 1 + static_cast<int>(rng() % 6);
    in.pins.resize(pins);
    for (int p = 0; p < pins; p++) {
      in.group.push_back(rng() % 5 == 0 && p > 0 ? in.group[p - 1] : p);
      const int cands = static_cast<int>(rng() % 6);
      for (int c = 0; c < cands; c++) {
        AccessCandidate a;
        a.length = static_cast<int64_t>(rng() % 50);
        a.end_layer = static_cast<int>(rng() % 2);
        a.end = {static_cast<GeoCoord>(rng() % 20),
                 static_cast<GeoCoord>(rng() % 20)};
        a.missing_directions = static_cast<int>(rng() % 3);
        const int nt = static_cast<int>(rng() % 3);
        for (int k = 0; k < nt; k++) {
          a.blocked_tracks.push_back(static_cast<int64_t>(rng() % 6));
        }
        std::sort(a.blocked_tracks.begin(), a.blocked_tracks.end());
        a.blocked_tracks.erase(
            std::unique(a.blocked_tracks.begin(), a.blocked_tracks.end()),
            a.blocked_tracks.end());
        in.pins[p].push_back(a);
      }
    }
    std::map<std::tuple<int, int, int, int>, bool> table;
    const int density = static_cast<int>(rng() % 4);
    for (int pa = 0; pa < pins; pa++) {
      for (int pb = pa + 1; pb < pins; pb++) {
        for (size_t ca = 0; ca < in.pins[pa].size(); ca++) {
          for (size_t cb = 0; cb < in.pins[pb].size(); cb++) {
            const bool c = static_cast<int>(rng() % 4) < density;
            table[{pa, (int) ca, pb, (int) cb}] = c;
            table[{pb, (int) cb, pa, (int) ca}] = c;
          }
        }
      }
    }
    in.conflict = [&table](int pa, int ca, int pb, int cb) {
      return table.at({pa, ca, pb, cb});
    };
    in.weights.track = static_cast<int64_t>(rng() % 10);
    in.weights.direction = static_cast<int64_t>(rng() % 10);
    in.weights.spread = static_cast<int64_t>(rng() % 3);
    in.weights.spread_distance = static_cast<GeoCoord>(rng() % 10);
    in.weights.unassigned = 1000;
    // Exhaustive.
    std::vector<int> choice(pins, -1);
    int64_t best = kInf;
    std::function<void(int)> rec = [&](int p) {
      if (p == pins) {
        const int64_t c = selectionCost(in, choice);
        if (c >= 0) {
          best = std::min(best, c);
        }
        return;
      }
      for (int c = -1; c < static_cast<int>(in.pins[p].size()); c++) {
        choice[p] = c;
        rec(p + 1);
      }
    };
    rec(0);
    const SelectionResult r = selectConflictFree(in);
    CHECK(r.optimal, "not optimal at it=%d", it);
    CHECK(r.cost == best,
          "it=%d cost %lld vs exhaustive %lld",
          it,
          static_cast<long long>(r.cost),
          static_cast<long long>(best));
    CHECK(selectionCost(in, r.choice) == r.cost, "cost mismatch");
    const SelectionResult g = selectGreedy(in);
    CHECK(g.cost >= best, "greedy below optimum");
    greedy_worse += g.cost > best;
    if (best > 0) {
      greedy_excess += double(g.cost - best) / best;
    }
    nodes += r.nodes;
    max_nodes = std::max(max_nodes, r.nodes);
    iterations += r.iterations;
    cases++;
  }
  std::printf(
      "destructive bounding: %d cases equal to exhaustive enumeration; mean "
      "%.1f nodes (max %lld), mean %.2f bound iterations; greedy worse in "
      "%d cases (mean excess %.1f%%)\n",
      cases,
      double(nodes) / cases,
      static_cast<long long>(max_nodes),
      double(iterations) / cases,
      greedy_worse,
      100.0 * greedy_excess / cases);
}

// ---------------------------------------------------------------------------
// Circuit signatures under isometries.
void testSignatures(std::mt19937& rng)
{
  int cases = 0;
  for (int it = 0; it < 500; it++) {
    // A local situation.
    std::vector<std::tuple<int, Box, int>> shapes;
    for (int k = 0; k < 12; k++) {
      const Box b = randomBoxes(rng, 1, 40, 8).front();
      shapes.emplace_back(
          static_cast<int>(rng() % 3), b, static_cast<int>(rng() % 3));
    }
    std::vector<std::tuple<int, bool, GeoCoord>> tracks;
    for (int k = 0; k < 6; k++) {
      tracks.emplace_back(static_cast<int>(rng() % 2),
                          rng() % 2 == 0,
                          static_cast<GeoCoord>(rng() % 40));
    }
    const bool symmetric_via = it % 2 == 0;
    const std::vector<std::pair<int, Box>> via
        = symmetric_via ? std::vector<std::pair<int, Box>>{{0, {-3, -3, 3, 3}},
                                                           {1, {-2, -2, 2, 2}}}
                        : std::vector<std::pair<int, Box>>{{0, {-3, -5, 3, 5}},
                                                           {1, {-2, -1, 2, 1}}};
    auto signature = [&](const Isometry& place, const Box* perturb) {
      // Chip situation = place(local); signature with place^-1.
      CircuitSignature sig(place.inverse());
      for (const auto& [layer, b, label] : shapes) {
        sig.addShape(
            layer,
            place.apply(perturb != nullptr && label == 2 ? *perturb : b),
            label);
      }
      for (const auto& [layer, horizontal, c] : tracks) {
        const GridPoint p
            = place.apply(horizontal ? GridPoint{0, c} : GridPoint{c, 0});
        const bool chip_horizontal = horizontal != place.swapsAxes();
        sig.addTrack(layer, chip_horizontal, chip_horizontal ? p.y : p.x);
      }
      // The via model is the same in chip coordinates for every instance.
      sig.addVia(0, via);
      sig.addPin(0, 1);
      sig.addPin(1, 2);
      return sig.finish();
    };
    CircuitClassTable table;
    std::set<int> classes;
    for (int o = 0; o < 8; o++) {
      Isometry place = dihedral(o);
      place.dx = static_cast<GeoCoord>(rng() % 1000);
      place.dy = static_cast<GeoCoord>(rng() % 1000);
      classes.insert(table.classify(signature(place, nullptr)).first);
      // Mapping rep -> member: member(local) == M(rep(local)).
      Isometry rep = dihedral(static_cast<int>(rng() % 8));
      rep.dx = 17;
      rep.dy = -5;
      const Isometry m = place.compose(rep.inverse());
      const Box lb = std::get<1>(shapes.front());
      CHECK(m.apply(rep.apply(lb)) == place.apply(lb), "mapping");
    }
    if (symmetric_via) {
      CHECK(classes.size() == 1, "symmetric via: %zu classes", classes.size());
    } else {
      // The asymmetric via distinguishes orientations by its image.
      std::set<std::vector<std::pair<int, Box>>> images;
      for (int o = 0; o < 8; o++) {
        std::vector<std::pair<int, Box>> img;
        for (const auto& [layer, b] : via) {
          img.emplace_back(layer, dihedral(o).inverse().apply(b));
        }
        std::sort(img.begin(), img.end());
        images.insert(img);
      }
      CHECK(classes.size() == images.size(),
            "asymmetric via: %zu vs %zu",
            classes.size(),
            images.size());
    }
    // A perturbed shape gives a new class.
    bool has_label2 = false;
    for (const auto& s : shapes) {
      has_label2 = has_label2 || std::get<2>(s) == 2;
    }
    if (has_label2) {
      const Box shifted{-100, -100, -90, -90};
      const int before = table.numClasses();
      table.classify(signature(Isometry(), &shifted));
      CHECK(table.numClasses() == before + 1, "perturbation not detected");
    }
    cases++;
  }
  std::printf(
      "circuit signatures: %d situations x 8 isometries classified "
      "consistently\n",
      cases);
}

// ---------------------------------------------------------------------------
// Post-processing invariants.
void testCleanup(std::mt19937& rng)
{
  CleanupStats total;
  int cases = 0;
  int area_cases = 0;
  for (int it = 0; it < 3000; it++) {
    const GeoCoord half = 1 + static_cast<GeoCoord>(rng() % 2);
    const GeoCoord tau = 2 + static_cast<GeoCoord>(rng() % 6);
    CleanupRequest req;
    CleanupLayer layer;
    layer.tau = tau;
    layer.half_width = half;
    layer.rules.width = 2 * half;
    layer.rules.notch_spacing = rng() % 2 == 0 ? 2 * half : 0;
    layer.grid = 1 + static_cast<GeoCoord>(rng() % 3);
    layer.horizontal = rng() % 2 == 0;
    req.layers = {layer, layer, layer};
    const bool on_tracks = rng() % 4 == 0;
    req.layers[0].on_tracks = on_tracks;
    req.layers[1].rules.min_area = rng() % 2 == 0 ? 60 + rng() % 120 : 0;
    std::vector<Box> obstacles;
    req.allowed = [&obstacles](int, const Box& b) {
      for (const Box& o : obstacles) {
        if (b.xlo < o.xhi && o.xlo < b.xhi && b.ylo < o.yhi && o.ylo < b.yhi) {
          return false;
        }
      }
      return true;
    };
    req.via_pads = [half](
                       int, const GridPoint& p, bool, std::vector<Box>& out) {
      out.push_back(
          {p.x - half - 1, p.y - half - 1, p.x + half + 1, p.y + half + 1});
    };
    // Pins away from their access points: no new nodes there.
    std::vector<std::pair<int, Box>> forbidden;
    if (rng() % 2 == 0) {
      for (const Box& b : randomBoxes(rng, 3, 70, 12)) {
        forbidden.emplace_back(static_cast<int>(rng() % 3), b);
      }
    }
    req.node_allowed = [&forbidden](int layer, const GridPoint& p) {
      for (const auto& [l, b] : forbidden) {
        if (l == layer && b.containsClosed(p.x, p.y)) {
          return false;
        }
      }
      return true;
    };
    // Layer 0 run, via, layer 1 piece, via, layer 2 run.
    const bool via_chain = rng() % 3 == 0;
    std::vector<GridPoint> run0 = randomPath(rng, tau, 1, half, 3 + rng() % 6);
    for (GridPoint& p : run0) {
      p.x += 30;
      p.y += 30;
    }
    std::vector<RouteNode> nodes;
    for (const GridPoint& p : run0) {
      nodes.push_back({0, p});
    }
    if (via_chain) {
      GridPoint a = run0.back();
      nodes.push_back({1, a});
      const GeoCoord len = static_cast<GeoCoord>(rng() % 6);
      GridPoint b{a.x + len, a.y};
      if (len > 0) {
        nodes.push_back({1, b});
      }
      nodes.push_back({2, b});
      nodes.push_back({2, {b.x, b.y + 2 * tau}});
    }
    req.nodes = nodes;
    std::vector<Box> old_metal;
    for (size_t k = 0; k + 1 < nodes.size(); k++) {
      if (nodes[k].layer == nodes[k + 1].layer) {
        old_metal.push_back(wireBox(nodes[k].p, nodes[k + 1].p, half));
      }
    }
    // Obstacles around, but not on, the existing metal and via pads.
    for (const Box& o : randomBoxes(rng, 12, 100, 10)) {
      bool hits = false;
      for (const Box& m : old_metal) {
        hits = hits || m.bloated(half + 1).touches(o);
      }
      if (!hits) {
        obstacles.push_back(o);
      }
    }
    const CleanupResult r = cleanupConnection(req);
    total.add(r.stats);
    cases++;
    CHECK(r.nodes.front() == nodes.front() && r.nodes.back() == nodes.back(),
          "ends moved");
    CHECK(r.stats.violations_after <= r.stats.violations_before, "worse");
    // Via nodes preserved in order.
    std::vector<RouteNode> vias_in;
    std::vector<RouteNode> vias_out;
    for (size_t k = 0; k + 1 < nodes.size(); k++) {
      if (nodes[k].layer != nodes[k + 1].layer) {
        vias_in.push_back(nodes[k]);
      }
    }
    for (size_t k = 0; k + 1 < r.nodes.size(); k++) {
      if (r.nodes[k].layer != r.nodes[k + 1].layer) {
        vias_out.push_back(r.nodes[k]);
        CHECK(r.nodes[k].p == r.nodes[k + 1].p, "via moved");
      }
    }
    CHECK(vias_in == vias_out, "vias changed");
    for (size_t k = 0; k + 1 < r.nodes.size(); k++) {
      const RouteNode& a = r.nodes[k];
      const RouteNode& b = r.nodes[k + 1];
      if (a.layer != b.layer) {
        continue;
      }
      CHECK(a.p.x == b.p.x || a.p.y == b.p.y, "diagonal segment");
      CHECK(req.allowed(a.layer, wireBox(a.p, b.p, half)), "illegal new metal");
    }
    if (on_tracks) {
      // No jog shifts or collapses on an on-track-only layer.
      for (const RouteNode& n : r.nodes) {
        CHECK(n.layer != 0
                  || std::find(nodes.begin(), nodes.end(), n) != nodes.end(),
              "node moved on an on-track-only layer");
      }
    }
    for (const CleanupPatch& p : r.patches) {
      const CleanupLayer& pl = req.layers[p.layer];
      CHECK(pl.horizontal ? p.begin.y == p.end.y : p.begin.x == p.end.x,
            "patch not in the preferred direction");
      CHECK((std::abs(p.end.x - p.begin.x) + std::abs(p.end.y - p.begin.y))
                    % pl.grid
                == 0,
            "patch length off the manufacturing grid");
      CHECK(req.allowed(p.layer, wireBox(p.begin, p.end, half)),
            "illegal patch");
      CHECK(req.node_allowed(p.layer, p.end), "patch tip on a pin");
      area_cases++;
    }
    for (const RouteNode& n : r.nodes) {
      const bool existing
          = std::find(nodes.begin(), nodes.end(), n) != nodes.end();
      CHECK(existing || req.node_allowed(n.layer, n.p), "new node on a pin");
    }
  }
  std::printf(
      "post-processing: %d connections, %lld runs; violations %lld -> %lld "
      "(%lld jog shifts, %lld U collapses); min-area pieces %lld, extended "
      "%lld (%d patches checked)\n",
      cases,
      static_cast<long long>(total.runs),
      static_cast<long long>(total.violations_before),
      static_cast<long long>(total.violations_after),
      static_cast<long long>(total.moved_segments),
      static_cast<long long>(total.collapsed_notches),
      static_cast<long long>(total.area_violations),
      static_cast<long long>(total.area_fixed),
      area_cases);
}

// ---------------------------------------------------------------------------
// Pin contacts (DRT-1231 regression): FlexDR needs every pin contact of the
// initial routing at an access point. Paths from access-point sources on a
// pin layer without wires (points_only) touch the pin only there, also
// after post-processing with node_allowed; region sources (the paper's
// "anywhere on the pin") do not.
void testPinContacts(std::mt19937& rng)
{
  int64_t ap_paths = 0;
  int64_t ap_violations = 0;
  int64_t region_paths = 0;
  int64_t region_violating_paths = 0;
  int64_t ap_vertices = 0;
  int64_t region_vertices = 0;
  int64_t guard_cases = 0;
  OffTrackSearch search;
  for (int it = 0; it < 1500; it++) {
    const GeoCoord tau = 4 + static_cast<GeoCoord>(rng() % 6);
    const Box pin{20,
                  20,
                  20 + 4 + static_cast<GeoCoord>(rng() % 8),
                  20 + 10 + static_cast<GeoCoord>(rng() % 20)};
    std::vector<GridPoint> aps;
    const int num_aps = 1 + static_cast<int>(rng() % 3);
    for (int k = 0; k < num_aps; k++) {
      aps.push_back(
          {pin.xlo + static_cast<GeoCoord>(rng() % (pin.width() + 1)),
           pin.ylo + static_cast<GeoCoord>(rng() % (pin.height() + 1))});
    }
    // Other pins on layer 0 (via blockages) and obstacles on layer 1.
    std::vector<Box> neighbours = randomBoxes(rng, 6, 60, 10);
    std::erase_if(neighbours, [&](const Box& b) { return b.touches(pin); });
    const std::vector<Box> upper = randomBoxes(rng, 4, 60, 8);
    std::vector<OffTrackTarget> targets;
    for (GeoCoord x = 4; x <= 56; x += 9) {
      for (GeoCoord y = 4; y <= 56; y += 7) {
        targets.push_back({1, {x, y}, static_cast<int>(targets.size())});
      }
    }
    auto via_free = [&](const GridPoint& p) {
      for (const Box& b : neighbours) {
        if (b.bloated(2).containsOpen(p.x, p.y)) {
          return false;
        }
      }
      return true;
    };
    std::vector<ContactPin> pins(1);
    pins[0].shapes.emplace_back(0, pin);
    for (const GridPoint& a : aps) {
      pins[0].access_points.push_back({0, a});
    }
    auto contacts = [](const std::vector<RouteNode>& nodes) {
      std::vector<RouteNode> c{nodes.front(), nodes.back()};
      for (size_t k = 0; k + 1 < nodes.size(); k++) {
        if (nodes[k].layer != nodes[k + 1].layer) {
          c.push_back(nodes[k]);
          c.push_back(nodes[k + 1]);
        }
      }
      return c;
    };
    for (const bool from_aps : {true, false}) {
      BlockageGridInput gin;
      OffTrackRequest req;
      req.window = {-10, -10, 70, 70};
      OffTrackLayer l0;
      l0.wires = false;
      OffTrackLayer l1;
      l1.obstacles = upper;
      gin.layers.push_back({1, {}, from_aps});
      gin.layers.push_back({tau, upper});
      std::vector<Box> via_blocks;
      for (const Box& b : neighbours) {
        via_blocks.push_back(b.bloated(2));
      }
      gin.via_blockages.push_back(via_blocks);
      for (const OffTrackTarget& t : targets) {
        gin.points.push_back(t.p);
      }
      if (from_aps) {
        for (const GridPoint& a : aps) {
          gin.points.push_back(a);
          req.sources.push_back({0, {a.x, a.y, a.x, a.y}});
        }
      } else {
        gin.borders.push_back(pin);
        req.sources.push_back({0, pin});
      }
      gin.borders.push_back(req.window);
      const auto grid = buildBlockageGrid(gin);
      l0.grid = grid[0];
      l1.grid = grid[1];
      req.layers = {l0, l1};
      req.targets = targets;
      req.via_cost = 1;
      req.max_targets = 8;
      // As BonnRouteDR: from access points, vias into the pin layer only at
      // the access points (points_only).
      req.via_allowed = [&](int, const GridPoint& p) {
        return via_free(p)
               && (!from_aps
                   || std::find(aps.begin(), aps.end(), p) != aps.end());
      };
      const OffTrackResult res = search.run(req);
      (from_aps ? ap_vertices : region_vertices) += res.grid_vertices;
      for (const OffTrackPath& path : res.paths) {
        const int64_t v = pinContactViolations(contacts(path.nodes), pins);
        if (from_aps) {
          ap_paths++;
          ap_violations += v;
          CHECK(std::find(aps.begin(), aps.end(), path.nodes.front().p)
                    != aps.end(),
                "path does not start at an access point");
          // Post-processing the path with an on-track tail keeps it.
          CleanupRequest creq;
          creq.nodes = path.nodes;
          const RouteNode end = path.nodes.back();
          creq.nodes.push_back({end.layer, {end.p.x + 1, end.p.y}});
          creq.nodes.push_back({end.layer, {end.p.x + 1, end.p.y + 20}});
          CleanupLayer cl;
          cl.tau = tau;
          cl.half_width = 1;
          creq.layers = {cl, cl};
          creq.node_allowed = [&](int layer, const GridPoint& p) {
            return pinContactViolations({{layer, p}}, pins) == 0;
          };
          const CleanupResult cr = cleanupConnection(creq);
          ap_violations += pinContactViolations(contacts(cr.nodes), pins);
        } else {
          region_paths++;
          region_violating_paths += v > 0;
        }
      }
    }
    // The decision itself against a direct count.
    std::vector<RouteNode> probe;
    for (int k = 0; k < 6; k++) {
      probe.push_back(
          {static_cast<int>(rng() % 2),
           k < 2 ? aps[rng() % aps.size()]
                 : GridPoint{static_cast<GeoCoord>(18 + rng() % 16),
                             static_cast<GeoCoord>(18 + rng() % 34)}});
    }
    int64_t want = 0;
    for (const RouteNode& c : probe) {
      const bool on = c.layer == 0 && pin.containsClosed(c.p.x, c.p.y);
      const bool at_ap
          = c.layer == 0 && std::find(aps.begin(), aps.end(), c.p) != aps.end();
      want += on && !at_ap;
    }
    CHECK(pinContactViolations(probe, pins) == want, "guard count");
    guard_cases++;
  }
  CHECK(ap_violations == 0,
        "%lld contacts away from access points",
        static_cast<long long>(ap_violations));
  CHECK(region_violating_paths > 0, "control found no violating path");
  std::printf(
      "pin contacts: access-point sources %lld paths (post-processed), %lld "
      "contacts away from access points; region sources %lld/%lld paths "
      "touch the pin away from an access point; grid vertices %lld vs %lld "
      "(points_only vs region); %lld guard counts checked\n",
      static_cast<long long>(ap_paths),
      static_cast<long long>(ap_violations),
      static_cast<long long>(region_violating_paths),
      static_cast<long long>(region_paths),
      static_cast<long long>(ap_vertices),
      static_cast<long long>(region_vertices),
      static_cast<long long>(guard_cases));
}

// ---------------------------------------------------------------------------
// Layers whose wires must stay on tracks in the preferred direction (asap7
// M1-M3, unidirectional): with the arc filter BonnRouteDR uses, every
// returned segment lies on a track in the preferred direction, and the
// costs equal a brute force restricted the same way.
void testOnTrackLayers(std::mt19937& rng)
{
  int cases = 0;
  int found = 0;
  OffTrackSearch search;
  for (int it = 0; it < 1500; it++) {
    const bool horizontal = rng() % 2 == 0;
    const GeoCoord pitch = 3 + static_cast<GeoCoord>(rng() % 4);
    std::vector<GeoCoord> tracks;
    for (GeoCoord c = 1; c < 40; c += pitch) {
      tracks.push_back(c);
    }
    const std::vector<Box> obs
        = randomBoxes(rng, static_cast<int>(rng() % 6), 40, 8);
    auto on_track_point = [&]() {
      const GeoCoord c = tracks[rng() % tracks.size()];
      const GeoCoord a = static_cast<GeoCoord>(rng() % 40);
      return horizontal ? GridPoint{a, c} : GridPoint{c, a};
    };
    const GridPoint s = on_track_point();
    const GridPoint t
        = rng() % 2 == 0
              ? on_track_point()
              : (horizontal
                     ? GridPoint{static_cast<GeoCoord>(rng() % 40), s.y}
                     : GridPoint{s.x, static_cast<GeoCoord>(rng() % 40)});
    const GeoCoord tau = 1 + static_cast<GeoCoord>(rng() % 4);
    BlockageGridInput gin;
    gin.layers.push_back({tau, obs});
    gin.points = {s, t};
    OffTrackRequest req;
    OffTrackLayer l;
    l.grid = buildBlockageGrid(gin)[0];
    l.obstacles = obs;
    req.layers = {l};
    req.window = {-1000, -1000, 1000, 1000};
    req.sources.push_back({0, {s.x, s.y, s.x, s.y}});
    req.targets.push_back({0, t, 0});
    req.arc_filter =
        [&](int, const GridPoint& from, const GridPoint&, int, const int dir) {
          const bool along_x = dir == 0 || dir == 1;
          return along_x == horizontal
                 && std::binary_search(tracks.begin(),
                                       tracks.end(),
                                       horizontal ? from.y : from.x);
        };
    const OffTrackResult res = search.run(req);
    // Restricted brute force: a straight preferred segment along a track,
    // at least tau long unless s == t, avoiding the obstacles.
    int64_t brute = kInf;
    const bool same_line = horizontal ? s.y == t.y : s.x == t.x;
    if (s == t) {
      brute = 0;
    } else if (same_line
               && std::binary_search(
                   tracks.begin(), tracks.end(), horizontal ? s.y : s.x)) {
      const int64_t len
          = std::abs(int64_t{s.x} - t.x) + std::abs(int64_t{s.y} - t.y);
      if (len >= tau && !edgeBlocked(obs, s, t)) {
        brute = len;
      }
    }
    const int64_t got = res.paths.empty() ? kInf : res.paths.front().cost;
    CHECK(got == brute,
          "on-track layer it=%d G=%lld brute=%lld",
          it,
          static_cast<long long>(got == kInf ? -1 : got),
          static_cast<long long>(brute == kInf ? -1 : brute));
    for (const OffTrackPath& p : res.paths) {
      for (size_t k = 0; k + 1 < p.nodes.size(); k++) {
        const GridPoint& a = p.nodes[k].p;
        const GridPoint& b = p.nodes[k + 1].p;
        if (a == b) {
          continue;
        }
        const bool along_x = a.y == b.y;
        CHECK(along_x == horizontal, "wrong-way segment");
        CHECK(std::binary_search(
                  tracks.begin(), tracks.end(), horizontal ? a.y : a.x),
              "off-track segment");
      }
    }
    found += !res.paths.empty();
    cases++;
  }
  std::printf(
      "on-track-only layers: %d instances equal to the restricted brute "
      "force (%d with a path), no wrong-way or off-track segment\n",
      cases,
      found);
}

// ---------------------------------------------------------------------------
// Sources with initial costs (access points pin access did not choose start
// with a penalty): the cost equals the minimum over the sources of initial
// cost + brute-force distance; on an on-track-only layer a blocked chosen
// access point leaves paths from the others in the catalogue.
void testSourceCosts(std::mt19937& rng)
{
  int cases = 0;
  int fallback_paths = 0;
  OffTrackSearch search;
  for (int it = 0; it < 1000; it++) {
    Instance in;
    Instance::Layer l;
    l.tau = 1 + static_cast<GeoCoord>(rng() % 5);
    l.obstacles = randomBoxes(rng, static_cast<int>(rng() % 6), 30, 8);
    in.layers.push_back(l);
    in.t.p = freePoint(rng, l.obstacles, 30);
    const GridPoint s1 = freePoint(rng, l.obstacles, 30);
    const GridPoint s2 = freePoint(rng, l.obstacles, 30);
    const int64_t c2 = static_cast<int64_t>(rng() % 30);
    const GeoCoord margin = 4 * l.tau + 4;
    const Box dom{-margin, -margin, 42 + margin, 42 + margin};
    in.s.p = s1;
    const int64_t b1 = bruteForce(in, dom);
    in.s.p = s2;
    const int64_t b2 = bruteForce(in, dom);
    const int64_t want = std::min(b1, b2 >= kInf ? kInf : b2 + c2);
    BlockageGridInput gin;
    gin.layers.push_back({l.tau, l.obstacles});
    gin.points = {s1, s2, in.t.p};
    OffTrackRequest req;
    OffTrackLayer ol;
    ol.grid = buildBlockageGrid(gin)[0];
    ol.obstacles = l.obstacles;
    req.layers = {ol};
    req.window = {-1000000, -1000000, 1000000, 1000000};
    req.sources
        = {{0, {s1.x, s1.y, s1.x, s1.y}, 0}, {0, {s2.x, s2.y, s2.x, s2.y}, c2}};
    req.targets = {{0, in.t.p, 0}};
    const OffTrackResult res = search.run(req);
    const int64_t got = res.paths.empty() ? kInf : res.paths.front().cost;
    CHECK(got == want,
          "source costs it=%d G=%lld want=%lld",
          it,
          static_cast<long long>(got == kInf ? -1 : got),
          static_cast<long long>(want == kInf ? -1 : want));
    cases++;
  }
  // On-track-only layer: chosen access point on track 10, another on
  // track 20; an obstacle across track 10 next to the chosen one.
  for (int it = 0; it < 200; it++) {
    const std::vector<GeoCoord> tracks{0, 10, 20, 30};
    const GridPoint chosen{5, 10};
    const GridPoint other{5, 20};
    const bool blocked = it % 2 == 0;
    std::vector<Box> obs;
    if (blocked) {
      // Both ways along track 10 enter an obstacle interior at once.
      obs.push_back({5, 8, 9, 12});
      obs.push_back({1, 8, 5, 12});
    }
    BlockageGridInput gin;
    gin.layers.push_back({1, obs});
    gin.points = {chosen, other};
    OffTrackRequest req;
    std::vector<OffTrackTarget> targets;
    for (const GeoCoord y : tracks) {
      for (GeoCoord x = 0; x <= 30; x += 6) {
        targets.push_back({0, {x, y}, static_cast<int>(targets.size())});
        gin.points.push_back({x, y});
      }
    }
    OffTrackLayer ol;
    ol.grid = buildBlockageGrid(gin)[0];
    ol.obstacles = obs;
    req.layers = {ol};
    req.window = {-100, -100, 100, 100};
    req.sources = {{0, {chosen.x, chosen.y, chosen.x, chosen.y}, 0},
                   {0, {other.x, other.y, other.x, other.y}, 20}};
    req.targets = targets;
    req.arc_filter =
        [&](int, const GridPoint& from, const GridPoint&, int, const int dir) {
          return (dir == 0 || dir == 1)
                 && std::binary_search(tracks.begin(), tracks.end(), from.y);
        };
    const OffTrackResult res = search.run(req);
    int from_other = 0;
    for (const OffTrackPath& p : res.paths) {
      from_other += p.nodes.front().p == other;
    }
    if (blocked) {
      CHECK(!res.paths.empty()
                && from_other == static_cast<int>(res.paths.size()),
            "blocked chosen access point: %zu paths, %d from the other",
            res.paths.size(),
            from_other);
      fallback_paths += from_other;
    } else {
      CHECK(!res.paths.empty() && res.paths.front().nodes.front().p == chosen,
            "unblocked: the chosen access point must be first");
    }
  }
  std::printf(
      "source costs: %d instances equal to the brute-force minimum "
      "over penalised sources; %d catalogue paths from a non-chosen "
      "access point where the chosen one is blocked\n",
      cases,
      fallback_paths);
}

}  // namespace

int main()
{
  std::mt19937 rng(20261111);
  const auto t0 = std::chrono::steady_clock::now();
  testAlgorithm3(rng);
  testTheorem(rng, 1, 3000);
  testTheorem(rng, 2, 1000);
  reportVertexGrowth();
  testSameNetMapping(rng);
  testSelection(rng);
  testSignatures(rng);
  testCleanup(rng);
  testPinContacts(rng);
  testOnTrackLayers(rng);
  testSourceCosts(rng);
  const double secs
      = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
            .count();
  std::printf("%s: %d failures (%.1fs)\n",
              failures == 0 ? "PASS" : "FAIL",
              failures,
              secs);
  return failures == 0 ? 0 : 1;
}
