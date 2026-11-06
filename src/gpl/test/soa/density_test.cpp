// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

#include "densitySoa.h"

int main()
{
  // Rectangular, rounded bins catch transposition and boundary mistakes that
  // a square, evenly divisible grid cannot detect.
  constexpr int nx = 8, ny = 16, left = -17, bottom = 9;
  constexpr double dx = 10.375, dy = 13.25;
  std::vector<int> lx, ly, ux, uy;
  std::vector<float> fx, fy;
  for (int x = 0; x < nx; ++x) {
    for (int y = 0; y < ny; ++y) {
      lx.push_back(left + std::lround(x * dx));
      ly.push_back(bottom + std::lround(y * dy));
      ux.push_back(left + std::lround((x + 1) * dx));
      uy.push_back(bottom + std::lround((y + 1) * dy));
      fx.push_back((x - 3 * y) * 0.125f);
      fy.push_back((y + 7 * x) * -0.0625f);
    }
  }
  std::vector<int> x_edges(nx + 1), y_edges(ny + 1);
  for (int x = 0; x <= nx; ++x) {
    x_edges[x] = left + std::lround(x * dx);
  }
  for (int y = 0; y <= ny; ++y) {
    y_edges[y] = bottom + std::lround(y * dy);
  }
  const gpl::DensityGridView grid{
      left, bottom, nx, ny, dx, dy, x_edges, y_edges, fx, fy};
  std::mt19937 random(42);
  std::vector<int64_t> actual(nx * ny, 0), expected(nx * ny, 0);
  for (int iteration = 0; iteration < 10000; ++iteration) {
    const int x = int(random() % 150) - 50;
    const int y = int(random() % 260) - 20;
    const gpl::DensityRect cell{x,
                                y,
                                x + int(random() % 100),
                                y + int(random() % 220),
                                float(random() % 100) / 113.0f};
    const bool macro = iteration % 3 == 0;
    float expected_x = 0, expected_y = 0;
    // Independent all-bin oracle, without the candidate's range calculation.
    for (int b = 0; b < nx * ny; ++b) {
      float area = 0;
      if (cell.lx < ux[b] && cell.ux > lx[b] && cell.ly < uy[b]
          && cell.uy > ly[b] && cell.lx < cell.ux && cell.ly < cell.uy) {
        const int w = std::min(ux[b], cell.ux) - std::max(lx[b], cell.lx);
        const int h = std::min(uy[b], cell.uy) - std::max(ly[b], cell.ly);
        area = float(w) * float(h) * cell.scale;
      }
      expected_x += area * fx[b];
      expected_y += area * fy[b];
      expected[b] += static_cast<int64_t>(macro ? area * 0.47f : area);
    }
    grid.forEachOverlap(cell, [&](size_t b, float area) {
      actual[b] += static_cast<int64_t>(macro ? area * 0.47f : area);
    });
    const auto [gx, gy] = grid.gather(cell);
    if (gx != expected_x || gy != expected_y || actual != expected) {
      throw std::runtime_error("resident density differs from all-bin oracle");
    }
  }
}
