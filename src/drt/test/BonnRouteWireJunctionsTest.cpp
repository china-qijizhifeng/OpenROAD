// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

#include "br/BrWireJunctions.h"

namespace {
using namespace drt::br;

// Independent pairwise geometric oracle for the sweep implementation.
std::vector<std::vector<GridPoint>> reference(
    const std::vector<WireSegment>& wires)
{
  std::vector<std::vector<GridPoint>> result(wires.size());
  for (size_t a = 0; a < wires.size(); a++) {
    const WireSegment& h = wires[a];
    if (h.begin.y != h.end.y || h.begin.x == h.end.x) {
      continue;
    }
    for (size_t b = 0; b < wires.size(); b++) {
      const WireSegment& v = wires[b];
      if (h.layer != v.layer || v.begin.x != v.end.x || v.begin.y == v.end.y) {
        continue;
      }
      const GridPoint p{v.begin.x, h.begin.y};
      if (p.x > std::min(h.begin.x, h.end.x)
          && p.x < std::max(h.begin.x, h.end.x)
          && p.y > std::min(v.begin.y, v.end.y)
          && p.y < std::max(v.begin.y, v.end.y)) {
        result[a].push_back(p);
        result[b].push_back(p);
      }
    }
  }
  for (auto& points : result) {
    std::ranges::sort(points);
    points.erase(std::ranges::unique(points).begin(), points.end());
  }
  return result;
}

bool fixtures()
{
  // jr10 net100: the last vertical wire overlaps the previous short
  // vertical. FlexDR merges them and loses the endpoint at the trunk.
  std::vector<WireSegment> wires{{4, {29670, 558450}, {45310, 558450}},
                                 {4, {33350, 558110}, {33350, 558450}},
                                 {4, {33350, 559130}, {33350, 558110}}};
  const auto jr = interiorWireCrossings(wires);
  if (jr[0] != std::vector<GridPoint>{{33350, 558450}} || !jr[1].empty()
      || jr[2] != jr[0]) {
    return false;
  }
  // ORFS net4663: a horizontal backtrack overlaps the branch's short
  // horizontal, hiding its endpoint inside the long vertical trunk.
  wires = {{6, {62330, 507450}, {62330, 493170}},
           {6, {62790, 502010}, {62330, 502010}},
           {6, {60490, 502010}, {62790, 502010}}};
  const auto orfs = interiorWireCrossings(wires);
  if (orfs[0] != std::vector<GridPoint>{{62330, 502010}} || !orfs[1].empty()
      || orfs[2] != orfs[0]) {
    return false;
  }
  // T-junctions, abutting wires and crossings on different layers already
  // keep their connectivity; avoid changing their output representation.
  wires = {{4, {0, 10}, {20, 10}},
           {4, {10, 0}, {10, 10}},
           {4, {20, 10}, {30, 10}},
           {6, {10, 0}, {10, 20}},
           {4, {10, 10}, {10, 10}}};
  const auto untouched = interiorWireCrossings(wires);
  return std::ranges::all_of(untouched,
                             [](const auto& v) { return v.empty(); });
}

bool generated()
{
  std::mt19937 rng(9311);
  for (int run = 0; run < 1000; run++) {
    std::vector<WireSegment> wires;
    for (int k = 0; k < 64; k++) {
      auto coordinate = [&]() { return static_cast<int>(rng() % 41) - 20; };
      GridPoint a{coordinate(), coordinate()};
      GridPoint b = a;
      if (rng() % 2) {
        b.x = coordinate();
      } else {
        b.y = coordinate();
      }
      wires.push_back({static_cast<int>(rng() % 3), a, b});
    }
    if (interiorWireCrossings(wires) != reference(wires)) {
      std::printf("wire junction mismatch in case %d\n", run);
      return false;
    }
    std::ranges::reverse(wires);
    for (auto& w : wires) {
      std::swap(w.begin, w.end);
    }
    if (interiorWireCrossings(wires) != reference(wires)) {
      return false;
    }
  }
  return true;
}
}  // namespace

int main()
{
  return fixtures() && generated() ? 0 : 1;
}
