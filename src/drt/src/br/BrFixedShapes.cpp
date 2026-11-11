// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrFixedShapes.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace drt::br {

int64_t insertFixedShapes(
    ShapeGrid& grid,
    const std::vector<std::vector<FixedShapeInsert>>& inserts,
    const ChunkAxis x_chunks,
    const ChunkAxis y_chunks,
    const int bands_per_layer,
    const BandParallelFor& parallel_for)
{
  const ShapeClassTable& classes = grid.classes();
  const int num_layers
      = std::min(grid.numLayers(), static_cast<int>(inserts.size()));
  struct Band
  {
    int layer;
    Coord lo;  // cross range [lo, hi) of the band's rows
    Coord hi;
    std::vector<int> items;  // indices into inserts[layer]
  };
  std::vector<Band> bands;
  std::vector<int> first_band(num_layers + 1, 0);
  for (int g = 0; g < num_layers; g++) {
    first_band[g] = static_cast<int>(bands.size());
    const ShapeGridLayerGeometry& geom = grid.geometry(g);
    const ChunkAxis& cross = geom.horizontal ? y_chunks : x_chunks;
    // First row of each tile row (rows run in increasing tile rows).
    std::vector<int> chunk_row(cross.count + 1, geom.cross_cells);
    for (int r = geom.cross_cells - 1; r >= 0; r--) {
      const int chunk = cross.of(static_cast<Coord>(
          geom.cross_origin + static_cast<int64_t>(r) * geom.cross_cell));
      for (int k = chunk; k >= 0 && chunk_row[k] > r; k--) {
        chunk_row[k] = r;
      }
    }
    const int count = std::max(1, std::min(bands_per_layer, cross.count));
    auto coord = [&](const int row) {
      return static_cast<Coord>(geom.cross_origin
                                + static_cast<int64_t>(row) * geom.cross_cell);
    };
    for (int b = 0; b < count; b++) {
      const int row_lo = chunk_row[b * cross.count / count];
      const int row_hi = chunk_row[(b + 1) * cross.count / count];
      bands.push_back(
          {g,
           b == 0 ? std::numeric_limits<Coord>::min() : coord(row_lo),
           b + 1 == count ? std::numeric_limits<Coord>::max() : coord(row_hi),
           {}});
    }
  }
  first_band[num_layers] = static_cast<int>(bands.size());
  parallel_for(num_layers, [&](const int g) {
    const bool horizontal = grid.geometry(g).horizontal;
    for (int k = 0; k < static_cast<int>(inserts[g].size()); k++) {
      const FixedShapeInsert& in = inserts[g][k];
      if (classes.wide(in.cls)) {
        bands[first_band[g]].items.push_back(k);
        continue;
      }
      const Coord lo = horizontal ? in.rect.ylo : in.rect.xlo;
      const Coord hi = horizontal ? in.rect.yhi : in.rect.xhi;
      for (int b = first_band[g]; b < first_band[g + 1]; b++) {
        if (lo < bands[b].hi && hi > bands[b].lo) {
          bands[b].items.push_back(k);
        }
      }
    }
  });
  // Shapes whose insertion failed in some band.
  std::vector<std::vector<int>> failed(bands.size());
  parallel_for(static_cast<int>(bands.size()), [&](const int b) {
    const Band& band = bands[b];
    const bool horizontal = grid.geometry(band.layer).horizontal;
    ShapeGridWriteCache cache;
    for (const int k : band.items) {
      const FixedShapeInsert& in = inserts[band.layer][k];
      Rect r = in.rect;
      if (!classes.wide(in.cls)) {
        Coord& lo = horizontal ? r.ylo : r.xlo;
        Coord& hi = horizontal ? r.yhi : r.xhi;
        lo = std::max(lo, band.lo);
        hi = std::min(hi, band.hi);
      }
      if (!grid.insert(band.layer, r, in.cls, in.type, in.owner, cache)) {
        failed[b].push_back(k);
      }
    }
  });
  int64_t failures = 0;
  for (int g = 0; g < num_layers; g++) {
    std::vector<int> all;
    for (int b = first_band[g]; b < first_band[g + 1]; b++) {
      all.insert(all.end(), failed[b].begin(), failed[b].end());
    }
    std::ranges::sort(all);
    failures += std::ranges::unique(all).begin() - all.begin();
  }
  return failures;
}

}  // namespace drt::br
