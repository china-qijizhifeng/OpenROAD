// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Database-free tests of the Sec. 2.5 capacity helpers (BonnRouteCapacity.h).

#include <algorithm>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "BonnRouteCapacity.h"
#include "gtest/gtest.h"
#include "odb/geom.h"

namespace grt::bonnroute {
namespace {

// Brute force over every placement of k disjoint horizontal runs: the
// expected maximum column load and the number of placements.
std::pair<double, double> bruteForce(int rows, int columns, int length, int k)
{
  std::vector<std::pair<int, int>> slots;  // (row, first column)
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c + length <= columns; c++) {
      slots.emplace_back(r, c);
    }
  }
  double placements = 0.0;
  double sum_max = 0.0;
  std::vector<int> chosen;
  std::vector<std::vector<int>> used(rows, std::vector<int>(columns, 0));
  auto recurse = [&](auto&& self, size_t from) -> void {
    if (static_cast<int>(chosen.size()) == k) {
      int max_load = 0;
      for (int c = 0; c < columns; c++) {
        int load = 0;
        for (int r = 0; r < rows; r++) {
          load += used[r][c];
        }
        max_load = std::max(max_load, load);
      }
      placements += 1.0;
      sum_max += max_load;
      return;
    }
    for (size_t i = from; i < slots.size(); i++) {
      const auto [r, c] = slots[i];
      bool free = true;
      for (int j = c; j < c + length; j++) {
        free = free && used[r][j] == 0;
      }
      if (!free) {
        continue;
      }
      for (int j = c; j < c + length; j++) {
        used[r][j] = 1;
      }
      chosen.push_back(static_cast<int>(i));
      self(self, i + 1);
      chosen.pop_back();
      for (int j = c; j < c + length; j++) {
        used[r][j] = 0;
      }
    }
  };
  recurse(recurse, 0);
  return {placements > 0 ? sum_max / placements : -1.0, placements};
}

TEST(BonnRouteCapacity, ExpectedMaxColumnLoadMatchesBruteForce)
{
  for (int rows = 1; rows <= 3; rows++) {
    for (int columns = 1; columns <= 5; columns++) {
      for (int length = 1; length <= std::min(3, columns); length++) {
        const int max_k = 5;
        const std::vector<double> table
            = expectedMaxColumnLoad(rows, columns, length, max_k);
        ASSERT_EQ(table.size(), static_cast<size_t>(max_k + 1));
        EXPECT_DOUBLE_EQ(table[0], 0.0);
        for (int k = 1; k <= max_k; k++) {
          const auto [expected, placements]
              = bruteForce(rows, columns, length, k);
          if (placements == 0.0) {
            // Infeasible counts repeat the last feasible value.
            EXPECT_DOUBLE_EQ(table[k], table[k - 1]);
            continue;
          }
          EXPECT_NEAR(table[k], expected, 1e-9)
              << rows << "x" << columns << " p=" << length << " k=" << k;
        }
      }
    }
  }
}

TEST(BonnRouteCapacity, ExpectedMaxColumnLoadIsSublinear)
{
  // One pad always blocks a full column position; more pads add less than
  // one track each, and never less than the average column load.
  const int rows = 15;
  const int columns = 18;
  const int length = 2;
  const std::vector<double> table
      = expectedMaxColumnLoad(rows, columns, length, 40);
  EXPECT_NEAR(table[1], 1.0, 1e-12);
  for (int k = 2; k <= 40; k++) {
    EXPECT_GE(table[k], table[k - 1] - 1e-12);
    EXPECT_LE(table[k] - table[k - 1], 1.0 + 1e-12);
    EXPECT_GE(table[k], static_cast<double>(k) * length / columns - 1e-9);
    EXPECT_LE(table[k], std::min<double>(k, rows) + 1e-9);
  }
  EXPECT_LT(table[10], 10.0 * 0.5);
}

TEST(BonnRouteCapacity, Interpolate)
{
  const std::vector<double> table{0.0, 1.0, 1.5};
  EXPECT_DOUBLE_EQ(interpolate(table, -1.0), 0.0);
  EXPECT_DOUBLE_EQ(interpolate(table, 0.5), 0.5);
  EXPECT_DOUBLE_EQ(interpolate(table, 1.5), 1.25);
  EXPECT_DOUBLE_EQ(interpolate(table, 7.0), 1.5);
}

TEST(BonnRouteCapacity, RectilinearSpanningLength)
{
  EXPECT_EQ(rectilinearSpanningLength({}),
            std::make_pair(int64_t{0}, int64_t{0}));
  EXPECT_EQ(rectilinearSpanningLength({{0, 0}, {30, 40}}),
            std::make_pair(int64_t{30}, int64_t{40}));
  // Three collinear points and one off the line.
  const auto [h, v]
      = rectilinearSpanningLength({{0, 0}, {100, 0}, {200, 0}, {100, 50}});
  EXPECT_EQ(h, 200);
  EXPECT_EQ(v, 50);
}

TEST(BonnRouteCapacity, TrackLayerBitsMatchNaive)
{
  std::mt19937 rng(7);
  const int tracks = 7;
  const int vertices = 203;
  std::vector<int> track_coords;
  std::vector<int> vertex_coords;
  for (int i = 0; i < tracks; i++) {
    track_coords.push_back(i * 100);
  }
  for (int i = 0; i < vertices; i++) {
    vertex_coords.push_back(i * 37);
  }
  TrackLayer layer(0, track_coords, vertex_coords);
  std::vector<std::vector<char>> naive(tracks, std::vector<char>(vertices, 0));
  for (int op = 0; op < 400; op++) {
    const int t = static_cast<int>(rng() % tracks);
    int a = static_cast<int>(rng() % (vertices + 20)) - 10;
    int b = a + static_cast<int>(rng() % 90);
    layer.block(t, a, b);
    for (int v = std::max(a, 0); v <= std::min(b, vertices - 1); v++) {
      naive[t][v] = 1;
    }
    const int q = static_cast<int>(rng() % tracks);
    const int lo = static_cast<int>(rng() % vertices);
    const int hi = lo + static_cast<int>(rng() % (vertices - lo + 1));
    int expected = 0;
    for (int v = lo; v < hi; v++) {
      expected += naive[q][v] == 0;
    }
    ASSERT_EQ(layer.countUsable(q, lo, hi), expected);
    const int v = static_cast<int>(rng() % vertices);
    ASSERT_EQ(layer.usable(q, v), naive[q][v] == 0);
  }
  EXPECT_EQ(layer.trackRange(150, 420), std::make_pair(2, 4));
  EXPECT_EQ(layer.vertexRange(40, 73), std::make_pair(2, 1));  // empty
  EXPECT_EQ(layer.nearestVertex(55), 1);
  EXPECT_EQ(layer.nearestVertex(56), 2);
  EXPECT_EQ(layer.nearestTrack(-5), 0);
  EXPECT_EQ(layer.lowerVertex(74), 2);
}

// The estimator works on plain data: three layers (V, H, V) with tracks at
// multiples of 100, gcells x 0 / 1500 / 3000 / 4500 and y 0 / 1500 / 3000.
CapacityGrid makeGrid()
{
  CapacityGrid grid;
  grid.gridlines[0] = {0, 1500, 3000, 4500};
  grid.gridlines[1] = {0, 1500, 3000};
  grid.die = odb::Rect(0, 0, 4500, 3000);
  for (int l = 0; l < 3; l++) {
    CapacityLayer layer;
    layer.direction = l == 1 ? kHorizontal : kVertical;
    const int extent = l == 1 ? 3000 : 4500;
    for (int c = 0; c <= extent; c += 100) {
      layer.tracks.push_back(c);
    }
    layer.pitch = 100;
    layer.width = 40;
    layer.routable = l >= 1;
    grid.layers.push_back(layer);
  }
  return grid;
}

TEST(BonnRouteCapacity, EstimatorOnPlainGridData)
{
  CapacityGrid grid = makeGrid();
  // Keep-out 0 + 20 - 1 = 19: tracks y = 1700..1900 lose x = 2000..2300.
  grid.layers[1].obstacles.emplace_back(2000, 1700, 2300, 1900);
  grid.layers[1].adjustment = 0.5;
  CapacityEstimator estimator(grid, {});
  estimator.buildTrackGraph();
  const std::vector<double> capacity = estimator.capacities();
  // Row 1 (the last) holds tracks 1500..3000; the spans [750, 2250) and
  // [2250, 3750) hold 15 vertices per track. Adjustment 0.5 halves them.
  EXPECT_DOUBLE_EQ(capacity[grid.index(1, 0, 1)], 0.5 * (16.0 - 9.0 / 15));
  EXPECT_DOUBLE_EQ(capacity[grid.index(1, 1, 1)], 0.5 * (16.0 - 3.0 / 15));
  EXPECT_DOUBLE_EQ(capacity[grid.index(1, 0, 0)], 0.5 * 15.0);
  EXPECT_DOUBLE_EQ(capacity[grid.index(1, 2, 1)], 0.0);  // last column
  EXPECT_DOUBLE_EQ(capacity[grid.index(0, 0, 0)], 0.0);  // unroutable
  EXPECT_DOUBLE_EQ(capacity[grid.index(2, 0, 0)], 15.0);
  EXPECT_NEAR(
      estimator.totals().tracks[1], 15.0 * 2 + 16.0 * 2 - 12.0 / 15, 1e-12);
}

}  // namespace
}  // namespace grt::bonnroute
