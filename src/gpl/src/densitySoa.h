// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once

#include <algorithm>
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
  std::span<const int> bin_lx, bin_ly, bin_ux, bin_uy;
  std::span<const float> field_x, field_y;

  template <typename Accumulate>
  void forEachOverlap(const DensityRect& cell, Accumulate&& add) const
  {
    const int x_begin = std::max(int((cell.lx - lx) / bin_size_x), 0);
    const int x_end = std::min(int(std::ceil((cell.ux - lx) / bin_size_x)), nx);
    const int y_begin = std::max(int((cell.ly - ly) / bin_size_y), 0);
    const int y_end = std::min(int(std::ceil((cell.uy - ly) / bin_size_y)), ny);
    for (int x = x_begin; x < x_end; ++x) {
      for (int y = y_begin; y < y_end; ++y) {
        const size_t b = static_cast<size_t>(x) * ny + y;
        const int width
            = std::min(bin_ux[b], cell.ux) - std::max(bin_lx[b], cell.lx);
        const int height
            = std::min(bin_uy[b], cell.uy) - std::max(bin_ly[b], cell.ly);
        const float area
            = width <= 0 || height <= 0
                  ? 0.0f
                  : static_cast<float>(width) * static_cast<float>(height);
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
