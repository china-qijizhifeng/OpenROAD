// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrBlockageGrid.h"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace drt::br {

namespace {

void sortUnique(std::vector<GeoCoord>& v)
{
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
}

// floor(a / b) and ceil(a / b) for b > 0.
int64_t floorDiv(int64_t a, int64_t b)
{
  return a >= 0 ? a / b : -((-a + b - 1) / b);
}
int64_t ceilDiv(int64_t a, int64_t b)
{
  return -floorDiv(-a, b);
}

}  // namespace

std::vector<GeoCoord> blockageGridCoordinates(
    const std::vector<GeoCoord>& hanan,
    const GeoCoord tau)
{
  // Lines 1-3: V' (the sentinels -inf, inf are implicit: the gap from the
  // first and to the last finite line is infinite, hence >= 4 tau).
  std::vector<GeoCoord> lines = hanan;
  sortUnique(lines);
  std::vector<GeoCoord> out;
  const int n = static_cast<int>(lines.size());
  const int64_t t = std::max<GeoCoord>(tau, 1);
  // Lines 5-6 for all i at once: runs of consecutive lines with gaps
  // < 4 tau. run_lo[i] = c_min, run_hi[i] = c_max.
  std::vector<int> run_lo(n);
  std::vector<int> run_hi(n);
  for (int i = 0; i < n; i++) {
    run_lo[i] = (i > 0 && int64_t{lines[i]} - lines[i - 1] < 4 * t)
                    ? run_lo[i - 1]
                    : i;
  }
  for (int i = n - 1; i >= 0; i--) {
    run_hi[i] = (i + 1 < n && int64_t{lines[i + 1]} - lines[i] < 4 * t)
                    ? run_hi[i + 1]
                    : i;
  }
  // Lines 4, 7-8.
  for (int i = 0; i < n; i++) {
    const int64_t xi = lines[i];
    const int64_t lo = int64_t{lines[run_lo[i]]} - 2 * t;
    const int64_t hi = int64_t{lines[run_hi[i]]} + 2 * t;
    const int64_t lambda_lo = ceilDiv(lo - xi, t);
    const int64_t lambda_hi = floorDiv(hi - xi, t);
    for (int64_t lambda = lambda_lo; lambda <= lambda_hi; lambda++) {
      out.push_back(static_cast<GeoCoord>(xi + lambda * t));
    }
  }
  // Line 11.
  sortUnique(out);
  return out;
}

std::vector<BlockageGridLayer> buildBlockageGrid(const BlockageGridInput& in)
{
  const int num = static_cast<int>(in.layers.size());
  std::vector<std::vector<GeoCoord>> own_x(num);
  std::vector<std::vector<GeoCoord>> own_y(num);
  auto add_box
      = [](std::vector<GeoCoord>& hx, std::vector<GeoCoord>& hy, const Box& b) {
          hx.push_back(b.xlo);
          hx.push_back(b.xhi);
          hy.push_back(b.ylo);
          hy.push_back(b.yhi);
        };
  auto pointsOnly = [&](const int k) { return in.layers[k].points_only; };
  for (int z = 0; z < num; z++) {
    std::vector<GeoCoord> hx;
    std::vector<GeoCoord> hy;
    if (!pointsOnly(z)) {
      for (int nz = z - 1; nz <= z + 1; nz++) {
        if (nz < 0 || nz >= num || pointsOnly(nz)) {
          continue;
        }
        for (const Box& b : in.layers[nz].blockages) {
          add_box(hx, hy, b);
        }
      }
      for (const int vz : {z - 1, z}) {
        if (vz >= 0 && vz < static_cast<int>(in.via_blockages.size())
            && !pointsOnly(vz) && !pointsOnly(vz + 1)) {
          for (const Box& b : in.via_blockages[vz]) {
            add_box(hx, hy, b);
          }
        }
      }
    }
    for (const GridPoint& p : in.points) {
      hx.push_back(p.x);
      hy.push_back(p.y);
    }
    for (const Box& b : in.borders) {
      add_box(hx, hy, b);
    }
    own_x[z] = blockageGridCoordinates(hx, in.layers[z].tau);
    own_y[z] = blockageGridCoordinates(hy, in.layers[z].tau);
  }
  std::vector<BlockageGridLayer> out(num);
  for (int z = 0; z < num; z++) {
    BlockageGridLayer& g = out[z];
    g.tau = in.layers[z].tau;
    for (int nz = z - 1; nz <= z + 1; nz++) {
      if (nz < 0 || nz >= num || (pointsOnly(z) && nz != z)) {
        continue;
      }
      g.xs.insert(g.xs.end(), own_x[nz].begin(), own_x[nz].end());
      g.ys.insert(g.ys.end(), own_y[nz].begin(), own_y[nz].end());
    }
    sortUnique(g.xs);
    sortUnique(g.ys);
  }
  return out;
}

}  // namespace drt::br
