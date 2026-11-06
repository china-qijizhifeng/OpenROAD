// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <utility>

namespace gpl {

// Call-local scalar registers, loaded once from the resident cell columns.
// There is no packed per-iteration copy of the cell population.
struct DensityRect
{
  int lx, ly, ux, uy;
  float scale;
};

// Non-owning views of resident columns. Recreate at a kernel boundary because
// topology changes can grow the columns. Grid storage is x * ny + y, matching
// the Poisson solver; bin boundaries retain OpenROAD's lround semantics.
struct DensityGridView
{
  int lx, ly, nx, ny;
  double bin_size_x, bin_size_y;
  std::span<const int> x_edges, y_edges;
  std::span<const float> field_x, field_y;

  template <typename Accumulate>
  void forEachOverlap(const DensityRect& cell, Accumulate&& add) const
  {
    const int x_begin = std::max(int((cell.lx - lx) / bin_size_x), 0);
    const int x_end = std::min(int(std::ceil((cell.ux - lx) / bin_size_x)), nx);
    const int y_begin = std::max(int((cell.ly - ly) / bin_size_y), 0);
    const int y_end = std::min(int(std::ceil((cell.uy - ly) / bin_size_y)), ny);
    // Most cells cover at most four bins per axis. Cache the y overlaps once;
    // large macros use the same separable calculation without bounded scratch.
    std::array<int, 8> heights;
    const bool small = y_end - y_begin <= static_cast<int>(heights.size());
    const auto height = [&](int y) {
      return std::min(y_edges[y + 1], cell.uy) - std::max(y_edges[y], cell.ly);
    };
    if (small) {
      for (int y = y_begin; y < y_end; ++y) {
        heights[y - y_begin] = height(y);
      }
    }
    for (int x = x_begin; x < x_end; ++x) {
      const int width
          = std::min(x_edges[x + 1], cell.ux) - std::max(x_edges[x], cell.lx);
      for (int y = y_begin; y < y_end; ++y) {
        const size_t b = static_cast<size_t>(x) * ny + y;
        const int h = small ? heights[y - y_begin] : height(y);
        const float area = width <= 0 || h <= 0 ? 0.0f
                                                : static_cast<float>(width)
                                                      * static_cast<float>(h);
        add(b, area * cell.scale);
      }
    }
  }

  std::pair<float, float> gather(const DensityRect& cell) const
  {
    float x = 0, y = 0;
    forEachOverlap(cell, [&](size_t b, float area) {
      // Keep the native x-outer/y-inner summation and untruncated weights.
      x += area * field_x[b];
      y += area * field_y[b];
    });
    return {x, y};
  }
};
}  // namespace gpl
