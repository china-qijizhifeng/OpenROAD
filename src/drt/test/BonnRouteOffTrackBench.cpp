// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Measurement of the off-track pin access modules on generated standard
// cell rows (Nangate45-like dimensions: M1 horizontal pitch 140 / width 70 /
// spacing 65 with the pins, M2 vertical pitch 190 / width 70 / spacing 70,
// rows of height 1400 alternating R0 / MX). Mirrors BonnRouteDR's use:
// circuit classes by signature, a catalogue of shortest tau-feasible paths
// per pin on the blockage grid of M1 + M2 to on-track points within two
// pitches, a conflict-free primary path per class by destructive bounding,
// and the class solution mapped to every member and re-verified there.
//
// usage: bonnroute_offtrack_bench [rows] [cells_per_row] [seed] [patterns]
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "br/BrBlockageGrid.h"
#include "br/BrGeometry.h"
#include "br/BrOffTrackSearch.h"
#include "br/BrPinAccess.h"
#include "br/BrSameNetRules.h"

using namespace drt::br;

namespace {

constexpr GeoCoord kSite = 190;
constexpr GeoCoord kRow = 1400;
constexpr GeoCoord kM1Pitch = 140;
constexpr GeoCoord kM1Offset = 70;
constexpr GeoCoord kM2Pitch = 190;
constexpr GeoCoord kM2Offset = 95;
constexpr GeoCoord kWidth = 70;
constexpr GeoCoord kHalf = 35;
constexpr GeoCoord kSpaceM1 = 65;
constexpr GeoCoord kSpaceM2 = 70;
constexpr GeoCoord kRadius = 2;
// via1: bottom pad on M1 (asymmetric), cut, top pad on M2.
const Box kPadM1{-35, -70, 35, 70};
const Box kPadM2{-70, -35, 70, 35};

struct Master
{
  int sites = 1;
  std::vector<std::vector<Box>> pins;  // per pin, rects in the cell frame
};

struct Instance
{
  int master = 0;
  Isometry place;  // cell frame -> chip
  std::vector<int> nets;
};

struct Shape
{
  int layer;  // 0 M1, 1 M2
  Box box;
  int net;  // -1 blockage
};

std::vector<Master> makeMasters(std::mt19937& rng)
{
  std::vector<Master> masters;
  const int sites[] = {3, 4, 4, 5, 6, 8, 12, 19};
  const int pins[] = {2, 3, 3, 3, 4, 5, 5, 4};
  for (int m = 0; m < 8; m++) {
    Master ms;
    ms.sites = sites[m];
    const GeoCoord width = sites[m] * kSite;
    // Pins: vertical bars between the rails, x at multiples of 10.
    std::vector<GeoCoord> xs;
    for (int p = 0; p < pins[m]; p++) {
      const GeoCoord slot = (width - 200) / pins[m];
      xs.push_back(100 + p * slot + static_cast<GeoCoord>(rng() % 6) * 10);
    }
    for (const GeoCoord x : xs) {
      const GeoCoord ylo = 250 + static_cast<GeoCoord>(rng() % 30) * 10;
      const GeoCoord yhi = ylo + 300 + static_cast<GeoCoord>(rng() % 40) * 10;
      std::vector<Box> rects{{x, ylo, x + 90, std::min<GeoCoord>(yhi, 1150)}};
      if (rng() % 3 == 0) {
        // L-shaped pin: a horizontal foot.
        rects.push_back({x, ylo, x + 200, ylo + 90});
      }
      ms.pins.push_back(rects);
    }
    masters.push_back(ms);
  }
  return masters;
}

bool violates(const Box& a, const Box& b, GeoCoord s)
{
  const int64_t dx = std::max({0, a.xlo - b.xhi, b.xlo - a.xhi});
  const int64_t dy = std::max({0, a.ylo - b.yhi, b.ylo - a.yhi});
  return dx * dx + dy * dy < int64_t{s} * s;
}

struct Chip
{
  std::vector<Instance> instances;
  std::vector<Shape> shapes;
  GeoCoord width = 0;
  GeoCoord height = 0;
};

// `patterns` > 0: every row repeats one of that many master sequences
// (datapath-like regularity); 0: independent random masters.
Chip makeChip(std::mt19937& rng,
              const std::vector<Master>& masters,
              int rows,
              int cells,
              int patterns)
{
  Chip chip;
  int net = 0;
  std::vector<std::vector<std::pair<int, bool>>> sequences(patterns);
  for (auto& seq : sequences) {
    for (int c = 0; c < cells; c++) {
      seq.emplace_back(static_cast<int>(rng() % masters.size()),
                       rng() % 4 == 0);
    }
  }
  for (int r = 0; r < rows; r++) {
    GeoCoord x = 0;
    const int pattern = patterns > 0 ? static_cast<int>(rng() % patterns) : -1;
    for (int c = 0; c < cells; c++) {
      Instance inst;
      inst.master = pattern >= 0 ? sequences[pattern][c].first
                                 : static_cast<int>(rng() % masters.size());
      const bool gap
          = pattern >= 0 ? sequences[pattern][c].second : rng() % 4 == 0;
      const Master& m = masters[inst.master];
      if (r % 2 == 0) {
        inst.place = {1, 0, 0, 1, x, r * kRow};
      } else {
        inst.place = {1, 0, 0, -1, x, (r + 1) * kRow};  // MX
      }
      for (size_t p = 0; p < m.pins.size(); p++) {
        inst.nets.push_back(net);
        for (const Box& b : m.pins[p]) {
          chip.shapes.push_back({0, inst.place.apply(b), net});
        }
        net++;
      }
      chip.instances.push_back(inst);
      x += m.sites * kSite + (gap ? kSite : 0);
    }
    chip.width = std::max(chip.width, x);
  }
  chip.height = rows * kRow;
  for (int r = 0; r <= rows; r++) {
    chip.shapes.push_back(
        {0, {0, r * kRow - 85, chip.width, r * kRow + 85}, -1});
  }
  return chip;
}

std::vector<GeoCoord> tracks(GeoCoord offset,
                             GeoCoord pitch,
                             GeoCoord lo,
                             GeoCoord hi)
{
  std::vector<GeoCoord> t;
  GeoCoord c = offset + ((lo - offset) / pitch) * pitch;
  for (; c <= hi; c += pitch) {
    if (c >= lo) {
      t.push_back(c);
    }
  }
  return t;
}

struct Catalogue
{
  std::vector<OffTrackPath> paths;
  int64_t grid_vertices = 0;
  int64_t settled = 0;
  double micros = 0;         // whole catalogue
  double search_micros = 0;  // OffTrackSearch::run only
};

std::vector<Box> pathBoxes(const OffTrackPath& path, std::vector<int>& layers)
{
  std::vector<Box> boxes;
  layers.clear();
  for (size_t k = 0; k + 1 < path.nodes.size(); k++) {
    const RouteNode& a = path.nodes[k];
    const RouteNode& b = path.nodes[k + 1];
    if (a.layer != b.layer) {
      boxes.push_back(Isometry{1, 0, 0, 1, a.p.x, a.p.y}.apply(kPadM1));
      layers.push_back(0);
      boxes.push_back(Isometry{1, 0, 0, 1, a.p.x, a.p.y}.apply(kPadM2));
      layers.push_back(1);
    } else if (!(a.p == b.p)) {
      boxes.push_back(wireBox(a.p, b.p, kHalf));
      layers.push_back(a.layer);
    }
  }
  return boxes;
}

// Legal against the chip's shapes of other nets.
bool legalPath(const Chip& chip, const OffTrackPath& path, int net)
{
  std::vector<int> layers;
  const std::vector<Box> boxes = pathBoxes(path, layers);
  for (size_t k = 0; k < boxes.size(); k++) {
    for (const Shape& s : chip.shapes) {
      if (s.layer == layers[k] && s.net != net
          && violates(s.box, boxes[k], layers[k] == 0 ? kSpaceM1 : kSpaceM2)) {
        return false;
      }
    }
  }
  return true;
}

Catalogue catalogue(const Chip& chip,
                    const std::vector<Box>& pin,
                    int net,
                    GeoCoord tau,
                    int max_paths)
{
  Catalogue out;
  const auto t0 = std::chrono::steady_clock::now();
  Box pin_box = pin.front();
  for (const Box& b : pin) {
    pin_box = pin_box.merged(b);
  }
  const Box reach{pin_box.xlo - kRadius * kM2Pitch,
                  pin_box.ylo - kRadius * kM1Pitch,
                  pin_box.xhi + kRadius * kM2Pitch,
                  pin_box.yhi + kRadius * kM1Pitch};
  const std::vector<GeoCoord> m1y
      = tracks(kM1Offset, kM1Pitch, reach.ylo, reach.yhi);
  const std::vector<GeoCoord> m2x
      = tracks(kM2Offset, kM2Pitch, reach.xlo, reach.xhi);
  const Box window = reach.bloated(2 * tau + kM2Pitch);
  OffTrackRequest req;
  BlockageGridInput gin;
  std::vector<Box> via_blocks;
  for (int z = 0; z < 2; z++) {
    OffTrackLayer l;
    for (const Shape& s : chip.shapes) {
      if (s.layer != z || s.net == net || !s.box.bloated(400).touches(window)) {
        continue;
      }
      const GeoCoord sp = z == 0 ? kSpaceM1 : kSpaceM2;
      l.obstacles.push_back(s.box.bloated(sp + kHalf));
      const Box& pad = z == 0 ? kPadM1 : kPadM2;
      via_blocks.push_back({s.box.xlo - sp - pad.xhi,
                            s.box.ylo - sp - pad.yhi,
                            s.box.xhi + sp - pad.xlo,
                            s.box.yhi + sp - pad.ylo});
    }
    gin.layers.push_back({tau, l.obstacles});
    req.layers.push_back(l);
  }
  gin.via_blockages.push_back(via_blocks);
  for (const GeoCoord y : m1y) {
    for (const GeoCoord x : m2x) {
      gin.points.push_back({x, y});
      for (int z = 0; z < 2; z++) {
        req.targets.push_back(
            {z, {x, y}, static_cast<int>(req.targets.size())});
      }
    }
  }
  for (const Box& b : pin) {
    gin.borders.push_back(b);
    req.sources.push_back({0, b});
  }
  gin.borders.push_back(window);
  const auto grid = buildBlockageGrid(gin);
  for (int z = 0; z < 2; z++) {
    req.layers[z].grid = grid[z];
  }
  req.window = window;
  req.via_cost = kWidth;  // as BonnRouteDR: one wire width
  req.max_targets = 2 * max_paths;
  req.via_allowed = [&](int, const GridPoint& p) {
    for (const Box& b : via_blocks) {
      if (b.containsOpen(p.x, p.y)) {
        return false;
      }
    }
    return true;
  };
  static OffTrackSearch search;  // buffers reused across pins
  const auto s0 = std::chrono::steady_clock::now();
  const OffTrackResult res = search.run(req);
  out.search_micros = std::chrono::duration<double, std::micro>(
                          std::chrono::steady_clock::now() - s0)
                          .count();
  out.grid_vertices = res.grid_vertices;
  out.settled = res.settled;
  for (const OffTrackPath& p : res.paths) {
    if (legalPath(chip, p, net)) {
      out.paths.push_back(p);
      if (static_cast<int>(out.paths.size()) >= max_paths) {
        break;
      }
    }
  }
  out.micros = std::chrono::duration<double, std::micro>(
                   std::chrono::steady_clock::now() - t0)
                   .count();
  return out;
}

AccessCandidate candidateOf(const Chip& chip, const OffTrackPath& p, int net)
{
  AccessCandidate c;
  c.length = p.cost;
  c.end_layer = p.nodes.back().layer;
  c.end = p.nodes.back().p;
  std::vector<int> layers;
  const std::vector<Box> boxes = pathBoxes(p, layers);
  for (size_t k = 0; k < boxes.size(); k++) {
    const Box& b = boxes[k];
    if (layers[k] == 1) {
      for (const GeoCoord x : tracks(kM2Offset,
                                     kM2Pitch,
                                     b.xlo - kSpaceM2 - kHalf + 1,
                                     b.xhi + kSpaceM2 + kHalf - 1)) {
        c.blocked_tracks.push_back((int64_t{1} << 32) | (x / kM2Pitch));
      }
    } else {
      for (const GeoCoord y : tracks(kM1Offset,
                                     kM1Pitch,
                                     b.ylo - kSpaceM1 - kHalf + 1,
                                     b.yhi + kSpaceM1 + kHalf - 1)) {
        c.blocked_tracks.push_back(y / kM1Pitch);
      }
    }
  }
  std::sort(c.blocked_tracks.begin(), c.blocked_tracks.end());
  c.blocked_tracks.erase(
      std::unique(c.blocked_tracks.begin(), c.blocked_tracks.end()),
      c.blocked_tracks.end());
  // Continuations along the endpoint's track: one pitch both ways.
  int feasible = 1;  // via up (M3 is free)
  const GeoCoord step = c.end_layer == 0 ? kM2Pitch : kM1Pitch;
  for (const int sgn : {-1, 1}) {
    const GridPoint q = c.end_layer == 0
                            ? GridPoint{c.end.x + sgn * step, c.end.y}
                            : GridPoint{c.end.x, c.end.y + sgn * step};
    const Box wb = wireBox(c.end, q, kHalf);
    bool ok = true;
    for (const Shape& s : chip.shapes) {
      ok = ok
           && !(s.layer == c.end_layer && s.net != net
                && violates(s.box, wb, c.end_layer == 0 ? kSpaceM1 : kSpaceM2));
    }
    feasible += ok;
  }
  c.missing_directions = 3 - feasible;
  return c;
}

double percentile(std::vector<double> v, double q)
{
  if (v.empty()) {
    return 0;
  }
  std::sort(v.begin(), v.end());
  return v[static_cast<size_t>(q * (v.size() - 1))];
}

}  // namespace

int main(int argc, char** argv)
{
  const int rows = argc > 1 ? std::atoi(argv[1]) : 12;
  const int cells = argc > 2 ? std::atoi(argv[2]) : 40;
  const unsigned seed = argc > 3 ? std::atoi(argv[3]) : 1;
  const int patterns = argc > 4 ? std::atoi(argv[4]) : 0;
  std::mt19937 rng(seed);
  const std::vector<Master> masters = makeMasters(rng);
  const Chip chip = makeChip(rng, masters, rows, cells, patterns);
  SameNetLayerRules rules;
  rules.width = kWidth;
  const GeoCoord tau = minimumSegmentLength(rules, false, 0).tau;
  // Circuit classes.
  const auto t0 = std::chrono::steady_clock::now();
  CircuitClassTable table;
  std::vector<int> cls(chip.instances.size());
  std::vector<int> rep;
  const GeoCoord halo
      = kRadius * kM2Pitch + 2 * tau + 3 * kM2Pitch;  // as BonnRouteDR
  for (size_t k = 0; k < chip.instances.size(); k++) {
    const Instance& inst = chip.instances[k];
    const Master& m = masters[inst.master];
    const Box cell = inst.place.apply(Box{0, 0, m.sites * kSite, kRow});
    const Box window = cell.bloated(halo);
    CircuitSignature sig(inst.place.inverse());
    for (const Shape& s : chip.shapes) {
      if (!s.box.touches(window)) {
        continue;
      }
      int label = 0;
      for (size_t p = 0; p < inst.nets.size(); p++) {
        if (inst.nets[p] == s.net) {
          label = static_cast<int>(p) + 1;
        }
      }
      sig.addShape(s.layer, s.box.clipped(window), label);
    }
    for (const GeoCoord y :
         tracks(kM1Offset, kM1Pitch, window.ylo, window.yhi)) {
      sig.addTrack(0, true, y);
    }
    for (const GeoCoord x :
         tracks(kM2Offset, kM2Pitch, window.xlo, window.xhi)) {
      sig.addTrack(1, false, x);
    }
    sig.addVia(0, {{0, kPadM1}, {1, kPadM2}});
    for (size_t p = 0; p < inst.nets.size(); p++) {
      sig.addPin(static_cast<int>(p), static_cast<int>(p) + 1);
    }
    const auto [c, is_new] = table.classify(sig.finish());
    cls[k] = c;
    if (is_new) {
      rep.push_back(static_cast<int>(k));
    }
  }
  const double class_ms = std::chrono::duration<double, std::milli>(
                              std::chrono::steady_clock::now() - t0)
                              .count();
  // Representative catalogues and selection.
  std::vector<double> search_us;
  std::vector<double> run_us;
  std::vector<double> settled;
  std::vector<double> sizes;
  std::vector<double> vertices;
  std::vector<double> sel_nodes;
  std::vector<std::vector<Catalogue>> class_cat(rep.size());
  std::vector<std::vector<int>> class_choice(rep.size());
  int64_t empty = 0;
  int64_t pins = 0;
  int64_t optimal = 0;
  int64_t db_unassigned = 0;
  int64_t greedy_unassigned = 0;
  int64_t db_cost = 0;
  int64_t greedy_cost = 0;
  double sel_us = 0;
  for (size_t c = 0; c < rep.size(); c++) {
    const Instance& inst = chip.instances[rep[c]];
    const Master& m = masters[inst.master];
    SelectionInput in;
    std::vector<std::vector<std::vector<Box>>> boxes;
    std::vector<std::vector<std::vector<int>>> blayers;
    for (size_t p = 0; p < m.pins.size(); p++) {
      std::vector<Box> pin;
      for (const Box& b : m.pins[p]) {
        pin.push_back(inst.place.apply(b));
      }
      Catalogue cat = catalogue(chip, pin, inst.nets[p], tau, 8);
      search_us.push_back(cat.micros);
      run_us.push_back(cat.search_micros);
      settled.push_back(static_cast<double>(cat.settled));
      sizes.push_back(static_cast<double>(cat.paths.size()));
      vertices.push_back(static_cast<double>(cat.grid_vertices));
      empty += cat.paths.empty();
      pins++;
      in.pins.emplace_back();
      boxes.emplace_back();
      blayers.emplace_back();
      for (const OffTrackPath& path : cat.paths) {
        in.pins.back().push_back(candidateOf(chip, path, inst.nets[p]));
        blayers.back().emplace_back();
        boxes.back().push_back(pathBoxes(path, blayers.back().back()));
      }
      in.group.push_back(static_cast<int>(p));
      class_cat[c].push_back(std::move(cat));
    }
    in.conflict = [&](int pa, int ca, int pb, int cb) {
      const auto& A = boxes[pa][ca];
      const auto& B = boxes[pb][cb];
      for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < B.size(); j++) {
          const int la = blayers[pa][ca][i];
          if (la == blayers[pb][cb][j]
              && violates(A[i], B[j], la == 0 ? kSpaceM1 : kSpaceM2)) {
            return true;
          }
        }
      }
      return false;
    };
    in.weights.track = kM1Pitch;
    in.weights.direction = kM1Pitch;
    in.weights.spread = 1;
    in.weights.spread_distance = 2 * kM1Pitch;
    in.node_limit = 20000;
    const auto s0 = std::chrono::steady_clock::now();
    const SelectionResult r = selectConflictFree(in);
    sel_us += std::chrono::duration<double, std::micro>(
                  std::chrono::steady_clock::now() - s0)
                  .count();
    const SelectionResult g = selectGreedy(in);
    sel_nodes.push_back(static_cast<double>(r.nodes));
    optimal += r.optimal;
    for (size_t p = 0; p < r.choice.size(); p++) {
      db_unassigned += r.choice[p] < 0 && !in.pins[p].empty();
      greedy_unassigned += g.choice[p] < 0 && !in.pins[p].empty();
    }
    db_cost += r.cost % in.weights.unassigned;
    greedy_cost += g.cost % in.weights.unassigned;
    class_choice[c] = r.choice;
  }
  // Map class solutions to members and verify them in place.
  int64_t mapped = 0;
  int64_t mapped_illegal = 0;
  for (size_t k = 0; k < chip.instances.size(); k++) {
    const int c = cls[k];
    const Instance& inst = chip.instances[k];
    const Instance& r = chip.instances[rep[c]];
    const Isometry m = inst.place.compose(r.place.inverse());
    for (size_t p = 0; p < class_cat[c].size(); p++) {
      for (const OffTrackPath& path : class_cat[c][p].paths) {
        OffTrackPath mp = path;
        for (RouteNode& n : mp.nodes) {
          n.p = m.apply(n.p);
        }
        mapped++;
        mapped_illegal += !legalPath(chip, mp, inst.nets[p]);
      }
    }
  }
  double mean_size = 0;
  for (const double s : sizes) {
    mean_size += s;
  }
  mean_size /= std::max<size_t>(1, sizes.size());
  std::printf("bench rows=%d cells/row=%d seed=%u row patterns=%d tau=%d\n",
              rows,
              cells,
              seed,
              patterns,
              tau);
  std::printf(
      "  instances %zu, circuit classes %d (%.1f instances/class), "
      "signatures %.1f ms\n",
      chip.instances.size(),
      table.numClasses(),
      double(chip.instances.size()) / table.numClasses(),
      class_ms);
  std::printf(
      "  representative pins %lld, empty catalogues %lld, catalogue "
      "size mean %.2f\n",
      static_cast<long long>(pins),
      static_cast<long long>(empty),
      mean_size);
  std::printf(
      "  search per pin: mean %.0f us, p50 %.0f us, p95 %.0f us, max "
      "%.0f us; grid vertices p50 %.0f, max %.0f\n",
      [&] {
        double s = 0;
        for (const double v : search_us) {
          s += v;
        }
        return s / std::max<size_t>(1, search_us.size());
      }(),
      percentile(search_us, 0.5),
      percentile(search_us, 0.95),
      percentile(search_us, 1.0),
      percentile(vertices, 0.5),
      percentile(vertices, 1.0));
  std::printf(
      "  of which OffTrackSearch::run: p50 %.0f us, p95 %.0f us; settled "
      "states p50 %.0f, max %.0f\n",
      percentile(run_us, 0.5),
      percentile(run_us, 0.95),
      percentile(settled, 0.5),
      percentile(settled, 1.0));
  std::printf(
      "  selection: %lld/%zu optimal, nodes p50 %.0f max %.0f, total "
      "%.1f ms; unassigned pins destructive bounding %lld vs greedy "
      "%lld; cost (assigned part) %lld vs greedy %lld\n",
      static_cast<long long>(optimal),
      rep.size(),
      percentile(sel_nodes, 0.5),
      percentile(sel_nodes, 1.0),
      sel_us / 1000,
      static_cast<long long>(db_unassigned),
      static_cast<long long>(greedy_unassigned),
      static_cast<long long>(db_cost),
      static_cast<long long>(greedy_cost));
  std::printf(
      "  class mapping: %lld mapped member paths, %lld illegal at the "
      "member\n",
      static_cast<long long>(mapped),
      static_cast<long long>(mapped_illegal));
  return mapped_illegal == 0 ? 0 : 1;
}
