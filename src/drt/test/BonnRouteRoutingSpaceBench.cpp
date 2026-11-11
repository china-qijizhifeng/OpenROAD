// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Measurements of the BonnRoute routing-space data structures (Sec. 3.3,
// 3.4, 3.6) on random large instances with Nangate45-like rules: shape-grid
// interval counts, configuration-table size, fast-grid interval counts, the
// fraction of path-search legality queries answered by the fast grid
// (paper: 97.89%), and query / update times against an R-tree holding the
// same shapes with the same pairwise rule (the previous BonnRouteDR path).
//
// usage: bonnroute_routing_space_bench [seed] [extent_dbu] [searches]

#include <algorithm>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <utility>
#include <vector>

#include "BonnRouteRoutingSpaceTestUtil.h"
#include "br/BrCore.h"
#include "br/BrDistanceRules.h"
#include "br/BrFastGrid.h"
#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrTrackGraph.h"
#include "br/BrWireModel.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using namespace drt::br;
using namespace drt::br::test;

namespace {

using BPoint = bg::model::point<int, 2, bg::cs::cartesian>;
using BBox = bg::model::box<BPoint>;
using RTree = bgi::rtree<std::pair<BBox, int>, bgi::quadratic<16>>;

double seconds(std::chrono::steady_clock::time_point a)
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - a)
      .count();
}

BBox box(const Rect& r)
{
  return {{r.xlo, r.ylo}, {r.xhi, r.yhi}};
}

}  // namespace

int main(int argc, char** argv)
{
  const unsigned seed = argc > 1 ? std::atoi(argv[1]) : 1;
  const Coord extent = argc > 2 ? std::atoi(argv[2]) : 120000;
  const int num_searches = argc > 3 ? std::atoi(argv[3]) : 4000;
  std::mt19937 rng(seed);
  const int num_z = 4;
  const Coord pitch[num_z] = {280, 280, 280, 560};
  const Coord width[num_z] = {140, 140, 140, 280};
  // Rules: PRL / wide-metal table, line-end rule, cut spacing, an
  // inter-layer via rule between the upper via layers.
  std::vector<LayerRules> rules;
  for (int z = 0; z < num_z; z++) {
    LayerRules w;
    w.horizontal = z % 2 == 0;
    const Coord base = z == 3 ? 160 : 130;
    w.spacing = SpacingTable::sample(
        {600, 601, 2000, 2001}, {1000, 1001}, [&](Coord a, Coord b, Coord prl) {
          const Coord m = std::max(a, b);
          return base + (m > 600 && prl > 1000 ? 50 : 0) + (m > 2000 ? 150 : 0);
        });
    w.line_end = {width[z] + 40, base + 30, width[z] / 2};
    rules.push_back(w);
    if (z + 1 < num_z) {
      LayerRules v;
      v.via_layer = true;
      v.horizontal = w.horizontal;
      v.spacing = SpacingTable::constant(160);
      v.projection_spacing = z >= 1 ? 120 : 0;
      rules.push_back(v);
    }
  }
  ShapeClassTable classes;
  classes.setLayers(rules);
  // Track graph.
  TrackGraph graph;
  graph.layers.resize(num_z);
  for (int z = 0; z < num_z; z++) {
    graph.layers[z].horizontal = rules[2 * z].horizontal;
    for (Coord t = pitch[z] / 2; t < extent; t += pitch[z]) {
      graph.layers[z].tracks.push_back(t);
    }
  }
  for (int z = 0; z < num_z; z++) {
    std::vector<Coord> pos;
    for (const int nz : {z - 1, z + 1}) {
      if (nz >= 0 && nz < num_z) {
        pos.insert(pos.end(),
                   graph.layers[nz].tracks.begin(),
                   graph.layers[nz].tracks.end());
      }
    }
    std::ranges::sort(pos);
    pos.erase(std::ranges::unique(pos).begin(), pos.end());
    graph.layers[z].positions = pos;
    graph.layers[z].jog_factor = 4;
  }
  for (int z = 0; z + 1 < num_z; z++) {
    graph.via_cost.push_back(4 * pitch[z]);
  }
  const Coord tile = 15 * pitch[0];
  const ChunkAxis axis{0, tile, extent / tile + 1};
  graph.finalize([&](Coord x) { return axis.of(x); },
                 [&](Coord y) { return axis.of(y); });
  // Wire type 0.
  WireTypeGeometry geometry;
  for (int z = 0; z < num_z; z++) {
    geometry.width.push_back(width[z]);
  }
  for (int z = 0; z + 1 < num_z; z++) {
    const Coord c = 70;
    geometry.via.push_back({true,
                            {-c - 35, -c, c + 35, c},
                            {-c, -c, c, c},
                            {-c, -c - 35, c, c + 35}});
  }
  std::vector<WireType> types{makeWireType(geometry, classes)};
  const WireType& wt = types[0];
  // Shape grid geometry: cell from the spacing bound, dividing the pitch.
  std::vector<ShapeGridLayerGeometry> layers;
  for (int g = 0; g < classes.numLayers(); g++) {
    ShapeGridLayerGeometry geom;
    geom.horizontal = rules[g].horizontal;
    const Coord bound = maximalCellSize(classes, g);
    auto cell = [&](Coord p) {
      for (int k = 1; k <= 64; k++) {
        if (p % k == 0 && p / k <= bound) {
          return p / k >= bound / 2 ? p / k : bound;
        }
      }
      return bound;
    };
    const int z = g / 2;
    geom.cross_cell = cell(pitch[z]);
    geom.along_cell = cell(pitch[std::min(z + 1, num_z - 1)]);
    geom.cross_origin = -geom.cross_cell / 2 - 4 * geom.cross_cell;
    geom.along_origin = -geom.along_cell / 2 - 4 * geom.along_cell;
    geom.cross_cells = (extent - geom.cross_origin) / geom.cross_cell + 8;
    geom.along_cells = (extent - geom.along_origin) / geom.along_cell + 8;
    layers.push_back(geom);
  }
  ShapeGrid grid;
  grid.init(&classes, layers, axis, axis);
  ShapeGridWriteCache cache;
  std::vector<RTree> trees(classes.numLayers());
  struct Shape
  {
    int layer;
    Rect rect;
    ClassId cls;
    ShapeType type;
    int owner;
  };
  std::vector<Shape> shapes;
  auto add = [&](int layer, Rect r, ShapeRole role, ShapeType type, int owner) {
    const ClassId cls = classes.classify(layer, r, role);
    shapes.push_back({layer, r, cls, type, owner});
  };
  // Fixed shapes: M1 rails every 10 tracks, pins in the rows, OBS on M2.
  int net = 0;
  for (Coord y = 0; y < extent; y += 10 * pitch[0]) {
    add(0,
        {0, y - 85, extent, y + 85},
        ShapeRole::kFixed,
        ShapeType::kBlockage,
        kNoNet);
    for (Coord x = 0; x < extent;) {
      x += 280 + rng() % 1400;
      const Coord yy = y + 400 + rng() % 1800;
      add(0,
          {x, yy, x + 140, yy + 400 + (Coord) (rng() % 300)},
          ShapeRole::kFixed,
          ShapeType::kPin,
          net++);
    }
  }
  for (int k = 0; k < extent / 4000; k++) {
    const Coord x = rng() % (extent - 6000);
    const Coord y = rng() % (extent - 6000);
    add(2,
        {x,
         y,
         x + 1500 + (Coord) (rng() % 4000),
         y + 1500 + (Coord) (rng() % 4000)},
        ShapeRole::kFixed,
        ShapeType::kBlockage,
        kNoNet);
  }
  // Off-track pin-access shapes on M2 (Sec. 4.3): small rectangles at
  // arbitrary positions, the origin of edges whose legality the endpoints
  // do not imply (Fig. 4).
  for (int k = 0; k < extent / 60; k++) {
    const Coord x = 1000 + rng() % (extent - 2000);
    const Coord y = 1000 + rng() % (extent - 2000);
    add(2,
        {x, y, x + 60 + (Coord) (rng() % 60), y + 60 + (Coord) (rng() % 60)},
        ShapeRole::kFixed,
        ShapeType::kPin,
        (int) (rng() % 1000));
  }
  // Routed wiring: on-track segments with vias at their ends, ~35% of the
  // track length per layer, plus off-track access jogs on M2.
  const int pins = net;
  std::vector<std::pair<int, Rect>> routed_wires;
  for (int z = 1; z < num_z; z++) {
    const TrackLayer& l = graph.layers[z];
    for (int i = 0; i < l.numTracks(); i++) {
      for (int j = 0; j + 1 < l.numPositions();) {
        j += 1 + rng() % 40;
        const int len = 2 + rng() % 30;
        if (j + len >= l.numPositions()) {
          break;
        }
        const int owner = pins + static_cast<int>(rng() % 50000);
        const auto [x0, y0] = l.point(i, j);
        const auto [x1, y1] = l.point(i, j + len);
        const Rect& m = wt.preferred[z].rect;
        const Rect r{x0 + m.xlo, y0 + m.ylo, x1 + m.xhi, y1 + m.yhi};
        add(2 * z, r, ShapeRole::kPreferredWire, ShapeType::kWire, owner);
        routed_wires.emplace_back(2 * z, r);
        // Via down at the start.
        const ViaModel& v = wt.via[z - 1];
        const Point p{x0, y0};
        add(2 * z - 2,
            v.bottom.rect.translated(p),
            ShapeRole::kViaPad,
            ShapeType::kVia,
            owner);
        add(2 * z - 1,
            v.cut.rect.translated(p),
            ShapeRole::kCut,
            ShapeType::kVia,
            owner);
        add(2 * z,
            v.top.rect.translated(p),
            ShapeRole::kViaPad,
            ShapeType::kVia,
            owner);
        if (v.projection.valid) {
          shapes.push_back({2 * z + 1,
                            v.projection.rect.translated(p),
                            v.projection.cls,
                            ShapeType::kVia,
                            owner});
        }
        if (z == 1 && rng() % 4 == 0) {
          // Off-track access jog of half a pitch.
          const Coord h = 70;
          const Rect jr = l.horizontal
                              ? Rect{x0 - h, y0 - h - 140, x0 + h, y0 + h}
                              : Rect{x0 - h - 140, y0 - h, x0 + h, y0 + h};
          add(2, jr, ShapeRole::kNonPreferredWire, ShapeType::kWire, owner);
        }
        j += len;
      }
    }
  }
  std::printf(
      "instance: seed %u, extent %d DBU, %d layers, %zu shapes, %d pins\n",
      seed,
      extent,
      num_z,
      shapes.size(),
      pins);
  // 1. Shape grid and R-tree construction.
  auto t0 = std::chrono::steady_clock::now();
  for (const Shape& s : shapes) {
    grid.insert(s.layer, s.rect, s.cls, s.type, s.owner, cache);
  }
  grid.compact();
  const double t_grid = seconds(t0);
  t0 = std::chrono::steady_clock::now();
  {
    std::vector<std::vector<std::pair<BBox, int>>> per(classes.numLayers());
    for (size_t k = 0; k < shapes.size(); k++) {
      per[shapes[k].layer].emplace_back(box(shapes[k].rect),
                                        static_cast<int>(k));
    }
    for (int g = 0; g < classes.numLayers(); g++) {
      trees[g] = RTree(per[g].begin(), per[g].end());
    }
  }
  const double t_tree = seconds(t0);
  const ShapeGrid::Statistics sg = grid.statistics();
  int64_t cells_total = 0;
  for (const ShapeGridLayerGeometry& g : layers) {
    cells_total += int64_t{g.along_cells} * g.cross_cells;
  }
  std::printf("| shape grid | value |\n|---|---|\n");
  std::printf("| shapes inserted | %zu |\n", shapes.size());
  std::printf("| cells (all layers) | %lld |\n", (long long) cells_total);
  std::printf("| non-empty cells in intervals | %lld |\n",
              (long long) sg.interval_cells);
  std::printf("| intervals (AVL nodes) | %lld |\n", (long long) sg.intervals);
  std::printf("| cells per interval | %.2f |\n",
              sg.intervals ? double(sg.interval_cells) / sg.intervals : 0.0);
  std::printf("| shared cells (several owners) | %lld |\n",
              (long long) sg.shared_cells);
  std::printf("| cell configurations | %zu |\n", sg.configurations);
  std::printf("| shape classes | %zu |\n", classes.size());
  std::printf("| tile chunks with trees | %lld |\n", (long long) sg.chunks);
  std::printf("| memory (trees, cell maps, configurations) | %.1f MB |\n",
              sg.bytes / 1e6);
  std::printf("| bytes per interval / per shape | %.1f / %.1f |\n",
              double(sg.bytes) / std::max<int64_t>(1, sg.intervals),
              double(sg.bytes) / shapes.size());
  std::printf("| build time shape grid / R-tree (s) | %.3f / %.3f |\n\n",
              t_grid,
              t_tree);
  // 2. Fast grid.
  DistanceRuleChecker checker(&grid);
  FastGrid fast;
  fast.init(&graph, 1);
  FastGridBuilder builder;
  builder.init(&graph, &types, 1, &checker, &fast);
  t0 = std::chrono::steady_clock::now();
  for (int z = 0; z < num_z; z++) {
    builder.rebuildVertices(z, 0, graph.layers[z].numTracks() - 1);
  }
  for (int z = 0; z < num_z; z++) {
    builder.rebuildJogBits(z, 0, graph.layers[z].numTracks() - 1);
  }
  fast.compact();
  const double t_fast = seconds(t0);
  const FastGrid::Statistics fs = fast.statistics();
  std::printf("| fast grid | value |\n|---|---|\n");
  std::printf("| on-track vertices | %lld |\n", (long long) fs.vertices);
  std::printf("| wiring intervals | %lld |\n", (long long) fs.wiring_intervals);
  std::printf("| vertices per interval | %.2f |\n",
              double(fs.vertices) / fs.wiring_intervals);
  std::printf("| via-layer vertices / intervals | %lld / %lld |\n",
              (long long) fs.via_vertices,
              (long long) fs.via_intervals);
  std::printf("| intervals with an edge query bit | %lld |\n",
              (long long) fs.marked_edges);
  std::printf("| allocated tiles | %lld |\n", (long long) fs.tiles);
  std::printf("| memory | %.1f MB |\n", fs.bytes / 1e6);
  std::printf(
      "| bytes per on-track vertex | %.3f (dense FastGrid before: 26) |\n",
      double(fs.bytes) / fs.vertices);
  std::printf("| bytes per interval | %.2f |\n",
              double(fs.bytes) / (fs.wiring_intervals + fs.via_intervals));
  std::printf("| build time (s, 1 thread) | %.3f |\n\n", t_fast);
  // 3. Path searches: fraction answered by the fast grid.
  int64_t fast_q = 0;
  int64_t checker_q = 0;
  int64_t found = 0;
  int64_t settled = 0;
  Region region;
  const int tiles = axis.count;
  region.init(num_z, tiles, tiles);
  PathSearch search(graph, fast);
  IntervalSearch interval_search(graph, fast);
  double t_path = 0;
  double t_interval = 0;
  int64_t cost_sum = 0;
  struct Dumped
  {
    size_t source;
    size_t target;
    int box[4];
    int a0, a1, b0, b1;
  };
  std::vector<Dumped> dumped;
  int64_t settled_interval = 0;
  t0 = std::chrono::steady_clock::now();
  for (int k = 0; k < num_searches; k++) {
    const int z = rng() % num_z;
    const TrackLayer& l = graph.layers[z];
    const int i = rng() % l.numTracks();
    const int j = rng() % l.numPositions();
    const int z2 = rng() % num_z;
    const TrackLayer& l2 = graph.layers[z2];
    const auto [x, y] = l.point(i, j);
    const Coord span = 30 * pitch[0];
    const Coord tx = std::clamp<Coord>(
        x + (Coord) (rng() % (2 * span)) - span, 0, extent - 1);
    const Coord ty = std::clamp<Coord>(
        y + (Coord) (rng() % (2 * span)) - span, 0, extent - 1);
    const int i2 = std::min<int>(
        l2.numTracks() - 1,
        std::ranges::lower_bound(l2.tracks, l2.horizontal ? ty : tx)
            - l2.tracks.begin());
    const int j2 = std::min<int>(
        l2.numPositions() - 1,
        std::ranges::lower_bound(l2.positions, l2.horizontal ? tx : ty)
            - l2.positions.begin());
    region.reset();
    for (int zz = 0; zz < num_z; zz++) {
      for (int a = axis.of(std::min(x, tx)) - 1;
           a <= axis.of(std::max(x, tx)) + 1;
           a++) {
        for (int b = axis.of(std::min(y, ty)) - 1;
             b <= axis.of(std::max(y, ty)) + 1;
             b++) {
          region.allow(zz, a, b);
        }
      }
    }
    SearchRequest q;
    q.region = &region;
    q.sources = {{graph.vertex(z, i, j), 0}};
    const size_t tv = graph.vertex(z2, i2, j2);
    q.targets = {{tv, 0, 0}};
    const auto [px, py] = graph.point(tv);
    q.target_boxes = {{px, py, px, py, z2, z2}};
    dumped.push_back({graph.vertex(z, i, j),
                      tv,
                      {px, py, z2, 0},
                      axis.of(std::min(x, tx)) - 1,
                      axis.of(std::max(x, tx)) + 1,
                      axis.of(std::min(y, ty)) - 1,
                      axis.of(std::max(y, ty)) + 1});
    auto t1 = std::chrono::steady_clock::now();
    const SearchResult r = search.run(q);
    t_path += seconds(t1);
    t1 = std::chrono::steady_clock::now();
    const SearchResult ri = interval_search.run(q);
    t_interval += seconds(t1);
    settled_interval += ri.settled;
    if (ri.found != r.found || (r.found && ri.cost != r.cost)) {
      std::printf("interval / vertex search differ\n");
    }
    cost_sum += r.found ? r.cost : 0;
    fast_q += r.fast_grid_queries;
    checker_q += r.checker_queries;
    found += r.found;
    settled += r.settled;
  }
  const double t_search = seconds(t0);
  std::printf("| on-track path search | value |\n|---|---|\n");
  std::printf(
      "| searches (found) | %d (%lld) |\n", num_searches, (long long) found);
  std::printf("| settled vertices | %lld |\n", (long long) settled);
  std::printf("| cost checksum of found paths | %lld |\n",
              (long long) cost_sum);
  if (const char* path = std::getenv("BR_DUMP")) {
    // Graph, per-vertex levels of every kind and the requests, for the
    // original dense-grid search (ds-kit/dense_search.cpp).
    FILE* f = std::fopen(path, "wb");
    auto put = [&](const void* p, size_t n) { std::fwrite(p, 1, n, f); };
    auto put_int = [&](int64_t v) { put(&v, sizeof(v)); };
    put_int(num_z);
    for (const TrackLayer& l : graph.layers) {
      put_int(l.horizontal);
      put_int(l.jogs);
      put_int(l.jog_factor);
      put_int(l.numTracks());
      for (const Coord c : l.tracks) {
        put_int(c);
      }
      put_int(l.numPositions());
      for (const Coord c : l.positions) {
        put_int(c);
      }
    }
    for (const Cost c : graph.via_cost) {
      put_int(c);
    }
    put_int(tile);
    put_int(axis.count);
    for (int z = 0; z < num_z; z++) {
      const TrackLayer& l = graph.layers[z];
      std::vector<uint8_t> lv(
          static_cast<size_t>(l.numTracks()) * l.numPositions() * kNumKinds, 0);
      for (int i = 0; i < l.numTracks(); i++) {
        for (int j = 0; j < l.numPositions(); j++) {
          uint8_t* out = &lv[(static_cast<size_t>(i) * l.numPositions() + j)
                             * kNumKinds];
          for (int k = 0; k < kNumKinds; k++) {
            const Kind kind = static_cast<Kind>(k);
            if ((kind == kCutUp || kind == kViaUpPad) && z + 1 >= num_z) {
              continue;
            }
            if (kind == kViaDownPad && z == 0) {
              continue;
            }
            if (kind == kCutUp && l.up_track[j] < 0) {
              continue;
            }
            if (kind == kWireEdge && j + 1 >= l.numPositions()) {
              continue;
            }
            if (kind == kJogEdge && i + 1 >= l.numTracks()) {
              continue;
            }
            out[k] = fast.level(z, i, j, kind, 0);
          }
        }
      }
      put(lv.data(), lv.size());
    }
    put_int(static_cast<int64_t>(dumped.size()));
    for (const Dumped& d : dumped) {
      put(&d, sizeof(d));
    }
    std::fclose(f);
  }
  std::printf(
      "| vertex search time (s) / per settled vertex | %.3f / %.0f ns |\n",
      t_path,
      1e9 * t_path / std::max<int64_t>(1, settled));
  std::printf(
      "| interval search time (s) / per settled vertex | %.3f / %.0f ns |\n",
      t_interval,
      1e9 * t_interval / std::max<int64_t>(1, settled_interval));
  std::printf("| legality queries | %lld |\n",
              (long long) (fast_q + checker_q));
  std::printf("| answered by the fast grid | %.4f%% (paper 97.89%%) |\n",
              100.0 * fast_q / std::max<int64_t>(1, fast_q + checker_q));
  std::printf("| answered by the distance rule checker | %lld |\n",
              (long long) checker_q);
  std::printf("| search time (s) | %.3f |\n", t_search);
  std::printf("| worker scratch after the searches | %.2f MB |\n",
              search.scratchBytes() / 1e6);
  {
    // One search whose window is the whole die (a die-spanning net).
    region.reset();
    for (int zz = 0; zz < num_z; zz++) {
      for (int a = 0; a < tiles; a++) {
        for (int b = 0; b < tiles; b++) {
          region.allow(zz, a, b);
        }
      }
    }
    // From a legal vertex near the lower left corner to the upper right.
    size_t sv = 0;
    for (int i = 4; i < graph.layers[1].numTracks(); i++) {
      if (fast.level(graph.vertex(1, i, 4), kWire, 0) == kLegal) {
        sv = graph.vertex(1, i, 4);
        break;
      }
    }
    SearchRequest q;
    q.region = &region;
    q.sources = {{sv, 0}};
    const size_t tv = graph.vertex(
        1, graph.layers[1].numTracks() - 5, graph.layers[1].numPositions() - 5);
    q.targets = {{tv, 0, 0}};
    const auto [px, py] = graph.point(tv);
    q.target_boxes = {{px, py, px, py, 1, 1}};
    PathSearch big(graph, fast);
    const SearchResult r = big.run(q);
    std::printf(
        "| die-wide search: window vertices / settled / found | %zu / %lld / "
        "%d |\n",
        graph.numVertices(),
        (long long) r.settled,
        (int) r.found);
    std::printf(
        "| die-wide search scratch (paged) | %.2f MB (dense arrays: %.2f MB) "
        "|\n\n",
        big.scratchBytes() / 1e6,
        graph.numVertices() * 29.0 / 1e6);
  }
  // 4. Query times.
  const int queries = 200000;
  std::vector<std::pair<int, Rect>> probes;
  for (int k = 0; k < queries; k++) {
    const int z = rng() % num_z;
    const Coord len = 140 + rng() % 3000;
    const Coord x = rng() % (extent - len);
    const Coord y = rng() % (extent - len);
    const Rect& m = wt.preferred[z].rect;
    probes.emplace_back(
        2 * z,
        rules[2 * z].horizontal
            ? Rect{x + m.xlo, y + m.ylo, x + len + m.xhi, y + m.yhi}
            : Rect{x + m.xlo, y + m.ylo, x + m.xhi, y + len + m.yhi});
  }
  int64_t sum_grid = 0;
  t0 = std::chrono::steady_clock::now();
  for (const auto& [g, r] : probes) {
    sum_grid += checker.shapeLevel(g, r, wt.preferred[g / 2].cls, 3, kNoNet);
  }
  const double t_q_grid = seconds(t0);
  int64_t sum_tree = 0;
  int64_t mismatches = 0;
  t0 = std::chrono::steady_clock::now();
  for (size_t k = 0; k < probes.size(); k++) {
    const auto& [g, r] = probes[k];
    const Coord reach = classes.reach(g);
    RipupLevel level = 0;
    const PlacedShape probe{r, wt.preferred[g / 2].cls, 3};
    for (auto it
         = trees[g].qbegin(bgi::intersects(box(r.bloated(reach, reach))));
         it != trees[g].qend();
         ++it) {
      const Shape& s = shapes[it->second];
      if (shapesViolate(classes, probe, PlacedShape{s.rect, s.cls, 3})) {
        level = std::max(level, grid.level(s.type, s.owner));
      }
    }
    sum_tree += level;
  }
  const double t_q_tree = seconds(t0);
  for (size_t k = 0; k < 2000; k++) {
    const auto& [g, r] = probes[k];
    const Coord reach = classes.reach(g);
    RipupLevel level = 0;
    const PlacedShape probe{r, wt.preferred[g / 2].cls, 3};
    for (auto it
         = trees[g].qbegin(bgi::intersects(box(r.bloated(reach, reach))));
         it != trees[g].qend();
         ++it) {
      const Shape& s = shapes[it->second];
      if (shapesViolate(classes, probe, PlacedShape{s.rect, s.cls, 3})) {
        level = std::max(level, grid.level(s.type, s.owner));
      }
    }
    mismatches
        += level
           != checker.shapeLevel(g, r, wt.preferred[g / 2].cls, 3, kNoNet);
  }
  // Point queries with intervals (fast-grid update path).
  std::vector<LocatedModel> models;
  appendWireTypeModels(wt, 1, num_z, 0, models);
  CheckOptions opts;
  opts.groups = 4;
  CheckResult res;
  t0 = std::chrono::steady_clock::now();
  int64_t interval_len = 0;
  for (int k = 0; k < queries; k++) {
    checker.check({(Coord) (rng() % extent), (Coord) (rng() % extent)},
                  models,
                  opts,
                  res);
    // Layer 1 is vertical: its tracks run along y.
    interval_len += int64_t{res.y_hi} - res.y_lo;
  }
  const double t_check = seconds(t0);
  // Fast grid lookups vs a dense per-vertex byte array (previous FastGrid).
  std::vector<size_t> verts;
  for (int k = 0; k < queries; k++) {
    verts.push_back(rng() % graph.numVertices());
  }
  int64_t sum_fast = 0;
  t0 = std::chrono::steady_clock::now();
  for (int rep = 0; rep < 5; rep++) {
    for (const size_t v : verts) {
      sum_fast += fast.level(v, kWire, 0) + fast.level(v, kJog, 0);
    }
  }
  const double t_fast_q = seconds(t0) / 10;
  // Dense per-vertex levels (the previous FastGrid layout), filled with the
  // fast grid's values so its pages are real.
  std::vector<uint8_t> dense(graph.numVertices() * 2, 0);
  for (int z = 0; z < num_z; z++) {
    const TrackLayer& l = graph.layers[z];
    for (int i = 0; i < l.numTracks(); i++) {
      for (int j = 0; j < l.numPositions(); j++) {
        const size_t v = graph.vertex(z, i, j);
        dense[2 * v] = fast.level(z, i, j, kWire, 0);
        dense[2 * v + 1] = fast.level(z, i, j, kJog, 0);
      }
    }
  }
  int64_t sum_dense = 0;
  t0 = std::chrono::steady_clock::now();
  for (int rep = 0; rep < 5; rep++) {
    for (const size_t v : verts) {
      sum_dense += dense[2 * v] + dense[2 * v + 1];
    }
  }
  const double t_dense_q = seconds(t0) / 10;
  // Sequential access along the tracks of layer 1 and the search-like
  // pattern (the vertex and its two neighbouring tracks).
  const TrackLayer& l1 = graph.layers[1];
  const int64_t seq_n = int64_t{l1.numTracks()} * l1.numPositions() * 2;
  int64_t sum_seq_fast = 0;
  t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < l1.numTracks(); i++) {
    for (int j = 0; j < l1.numPositions(); j++) {
      sum_seq_fast
          += fast.level(1, i, j, kWire, 0) + fast.level(1, i, j, kJog, 0);
    }
  }
  const double t_seq_fast = seconds(t0);
  int64_t sum_seq_dense = 0;
  t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < l1.numTracks(); i++) {
    for (int j = 0; j < l1.numPositions(); j++) {
      const size_t v = graph.vertex(1, i, j);
      sum_seq_dense += dense[2 * v] + dense[2 * v + 1];
    }
  }
  const double t_seq_dense = seconds(t0);
  const int64_t nb_n = int64_t{l1.numTracks() - 2} * l1.numPositions() * 3;
  int64_t sum_nb_fast = 0;
  t0 = std::chrono::steady_clock::now();
  for (int i = 1; i + 1 < l1.numTracks(); i++) {
    for (int j = 0; j < l1.numPositions(); j++) {
      sum_nb_fast += fast.level(1, i - 1, j, kWire, 0)
                     + fast.level(1, i, j, kWire, 0)
                     + fast.level(1, i + 1, j, kWire, 0);
    }
  }
  const double t_nb_fast = seconds(t0);
  int64_t sum_nb_dense = 0;
  t0 = std::chrono::steady_clock::now();
  for (int i = 1; i + 1 < l1.numTracks(); i++) {
    for (int j = 0; j < l1.numPositions(); j++) {
      sum_nb_dense += dense[2 * graph.vertex(1, i - 1, j)]
                      + dense[2 * graph.vertex(1, i, j)]
                      + dense[2 * graph.vertex(1, i + 1, j)];
    }
  }
  const double t_nb_dense = seconds(t0);
#ifdef BR_FAST_GRID_READER
  FastGridReader reader;
  reader.reset(&fast, nullptr);
  int64_t sum_seq_reader = 0;
  t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < l1.numTracks(); i++) {
    for (int j = 0; j < l1.numPositions(); j++) {
      sum_seq_reader
          += reader.level(1, i, j, kWire, 0) + reader.level(1, i, j, kJog, 0);
    }
  }
  const double t_seq_reader = seconds(t0);
  int64_t sum_nb_reader = 0;
  t0 = std::chrono::steady_clock::now();
  for (int i = 1; i + 1 < l1.numTracks(); i++) {
    for (int j = 0; j < l1.numPositions(); j++) {
      sum_nb_reader += reader.level(1, i - 1, j, kWire, 0)
                       + reader.level(1, i, j, kWire, 0)
                       + reader.level(1, i + 1, j, kWire, 0);
    }
  }
  const double t_nb_reader = seconds(t0);
  if (sum_seq_reader != sum_seq_fast || sum_nb_reader != sum_nb_fast
      || sum_seq_dense != sum_seq_fast || sum_nb_dense != sum_nb_fast) {
    std::printf("lookup checksums differ\n");
  }
#endif
  std::printf("| query (%d random) | time per query |\n|---|---|\n", queries);
  std::printf("| placed wire, shape grid checker | %.0f ns |\n",
              1e9 * t_q_grid / queries);
  std::printf("| placed wire, R-tree + same rule | %.0f ns |\n",
              1e9 * t_q_tree / queries);
  std::printf("| answers differing (2000 compared) | %lld |\n",
              (long long) mismatches);
  std::printf(
      "| point query, wire type models + maximal intervals | %.0f ns (mean "
      "interval along the track %.0f DBU) |\n",
      1e9 * t_check / queries,
      double(interval_len) / queries);
  std::printf("| fast grid level, random vertex | %.1f ns |\n",
              1e9 * t_fast_q / queries);
  std::printf("| dense per-vertex array, random vertex | %.1f ns |\n\n",
              1e9 * t_dense_q / queries);
  std::printf(
      "| lookup per level (layer 1, %lld) | sequential along tracks | vertex + "
      "2 neighbour tracks |\n|---|---|---|\n",
      (long long) seq_n);
  std::printf(
      "| dense per-vertex array (previous FastGrid) | %.2f ns | %.2f ns |\n",
      1e9 * t_seq_dense / seq_n,
      1e9 * t_nb_dense / nb_n);
  std::printf("| FastGrid::level(z, i, j) | %.2f ns | %.2f ns |\n",
              1e9 * t_seq_fast / seq_n,
              1e9 * t_nb_fast / nb_n);
#ifdef BR_FAST_GRID_READER
  std::printf("| FastGridReader (run cursors) | %.2f ns | %.2f ns |\n",
              1e9 * t_seq_reader / seq_n,
              1e9 * t_nb_reader / nb_n);
#endif
  std::printf("\n");
  // 5. Updates: insert + remove a routed wire, with the fast-grid update.
  const int updates = 20000;
  t0 = std::chrono::steady_clock::now();
  for (int k = 0; k < updates; k++) {
    const auto& [g, r] = probes[k];
    const ClassId cls = wt.preferred[g / 2].cls;
    grid.insert(g, r, cls, ShapeType::kWire, 7, cache);
#ifdef BR_FAST_GRID_READER
    builder.shapeInserted(
        g, PlacedShape{r, cls, 3}, grid.level(ShapeType::kWire, 7));
#else
    builder.shapeChanged(g, r);
#endif
    grid.remove(g, r, cls, ShapeType::kWire, 7, cache);
#ifdef BR_FAST_GRID_READER
    builder.shapeRemoved(
        g, PlacedShape{r, cls, 3}, grid.level(ShapeType::kWire, 7));
#else
    builder.shapeChanged(g, r);
#endif
  }
  const double t_update = seconds(t0);
  // Insertions, then removals of the same shapes (fast grid part only).
  double t_ins = 0;
  double t_rem = 0;
  {
    const int m = 5000;
    for (int k = 0; k < m; k++) {
      const auto& [g, r] = probes[k];
      grid.insert(g, r, wt.preferred[g / 2].cls, ShapeType::kWire, 7, cache);
      auto t1 = std::chrono::steady_clock::now();
#ifdef BR_FAST_GRID_READER
      builder.shapeInserted(g,
                            PlacedShape{r, wt.preferred[g / 2].cls, 3},
                            grid.level(ShapeType::kWire, 7));
#else
      builder.shapeChanged(g, r);
#endif
      t_ins += seconds(t1);
    }
    for (int k = 0; k < m; k++) {
      const auto& [g, r] = probes[k];
      grid.remove(g, r, wt.preferred[g / 2].cls, ShapeType::kWire, 7, cache);
      auto t1 = std::chrono::steady_clock::now();
#ifdef BR_FAST_GRID_READER
      builder.shapeRemoved(g,
                           PlacedShape{r, wt.preferred[g / 2].cls, 3},
                           grid.level(ShapeType::kWire, 7));
#else
      builder.shapeChanged(g, r);
#endif
      t_rem += seconds(t1);
    }
    t_ins /= m;
    t_rem /= m;
  }
  t0 = std::chrono::steady_clock::now();
  for (int k = 0; k < updates; k++) {
    const auto& [g, r] = probes[k];
    grid.insert(g, r, wt.preferred[g / 2].cls, ShapeType::kWire, 7, cache);
    grid.remove(g, r, wt.preferred[g / 2].cls, ShapeType::kWire, 7, cache);
  }
  const double t_update_grid = seconds(t0);
  t0 = std::chrono::steady_clock::now();
  for (int k = 0; k < updates; k++) {
    const auto& [g, r] = probes[k];
    trees[g].insert({box(r), -1});
    trees[g].remove(std::make_pair(box(r), -1));
  }
  const double t_update_tree = seconds(t0);
  std::printf(
      "| update (%d insert + remove pairs) | time per pair |\n|---|---|\n",
      updates);
  std::printf("| shape grid + fast grid update | %.1f us |\n",
              1e6 * t_update / updates);
  std::printf("| shape grid only | %.1f us |\n", 1e6 * t_update_grid / updates);
  std::printf("| fast grid update per shape change | %.1f us |\n",
              1e6 * (t_update - t_update_grid) / (2 * updates));
  std::printf(
      "| fast grid update: insert / remove (5000 wires) | %.1f us / %.1f us "
      "|\n",
      1e6 * t_ins,
      1e6 * t_rem);
  std::printf("| R-tree only | %.1f us |\n", 1e6 * t_update_tree / updates);
  // 6. A net's own shapes around a connection (Sec. 4.4): pins as their
  // nets' own shapes, removed and reinserted physically, or (after) read
  // through the overlay.
  std::vector<size_t> pin_ids;
  for (size_t k = 0; k < shapes.size() && pin_ids.size() < 20000; k++) {
    if (shapes[k].type == ShapeType::kPin && shapes[k].layer == 0) {
      pin_ids.push_back(k);
    }
  }
  t0 = std::chrono::steady_clock::now();
  for (const size_t k : pin_ids) {
    const Shape& s = shapes[k];
    grid.remove(s.layer, s.rect, s.cls, s.type, s.owner, cache);
#ifdef BR_FAST_GRID_READER
    builder.shapeRemoved(s.layer, PlacedShape{s.rect, s.cls, 3}, kNotRemovable);
#else
    builder.shapeChanged(s.layer, s.rect);
#endif
    grid.insert(s.layer, s.rect, s.cls, s.type, s.owner, cache);
#ifdef BR_FAST_GRID_READER
    builder.shapeInserted(
        s.layer, PlacedShape{s.rect, s.cls, 3}, kNotRemovable);
#else
    builder.shapeChanged(s.layer, s.rect);
#endif
  }
  const double t_own_physical = seconds(t0);
  std::printf(
      "| own pin around a connection (%zu pins) | time per pin |\n|---|---|\n",
      pin_ids.size());
  std::printf("| remove + reinsert (shape grid and fast grid) | %.1f us |\n",
              1e6 * t_own_physical / std::max<size_t>(1, pin_ids.size()));
#ifdef BR_FAST_GRID_READER
  FastGridOverlay overlay;
  std::vector<FastGridBuilder::OwnShape> own;
  int64_t overlay_positions = 0;
  t0 = std::chrono::steady_clock::now();
  for (const size_t k : pin_ids) {
    const Shape& s = shapes[k];
    own.assign(1, {s.layer, PlacedShape{s.rect, s.cls, 3}, kNotRemovable});
    builder.buildOverlay(s.owner, own, overlay);
    for (const auto& [key, runs] : overlay.wiring) {
      for (const auto& run : runs) {
        overlay_positions += run.hi - run.lo + 1;
      }
    }
  }
  const double t_own_overlay = seconds(t0);
  std::printf(
      "| overlay with the pin removed (after) | %.1f us (%.1f positions) |\n",
      1e6 * t_own_overlay / std::max<size_t>(1, pin_ids.size()),
      double(overlay_positions) / std::max<size_t>(1, pin_ids.size()));
#endif
  std::printf("\n");
  std::printf("| total routing-space memory per on-track vertex | %.2f B |\n",
              double(sg.bytes + fs.bytes) / fs.vertices);
  std::printf("(checksums %lld %lld %lld %lld)\n",
              (long long) sum_grid,
              (long long) sum_tree,
              (long long) sum_fast,
              (long long) sum_dense);
  return mismatches != 0;
}
