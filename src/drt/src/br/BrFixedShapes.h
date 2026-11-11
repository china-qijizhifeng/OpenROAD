// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Parallel insertion of the fixed shapes (blockages, pins; Sec. 3.3) into
// the shape grid, with the result of inserting them one by one.
//
// The shape grid of a layer is chunked by global routing tiles (Sec. 5.1:
// writers of disjoint tile rows never touch the same tree or wide-shape
// list). The shapes of a layer are split into bands of tile rows; the
// thread of a band inserts, in the given order, every shape clipped to the
// band's rows. The pieces a shape leaves in a row depend only on its extent
// in that row, so every tile chunk sees the sequence of changes of one
// serial pass. Wide shapes (kept whole in the wide lists) are inserted
// unclipped, in order, by the layer's first band. Configuration ids may be
// numbered differently (the table is shared); configurations, intervals
// and wide lists are those of the serial pass.

#pragma once

#include <cstdint>
#include <functional>
#include <vector>

#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrWireModel.h"

namespace drt::br {

struct FixedShapeInsert
{
  Rect rect;
  ClassId cls = 0;
  ShapeType type = ShapeType::kBlockage;
  int32_t owner = kNoNet;
};

// Runs body(0) .. body(count - 1), possibly concurrently.
using BandParallelFor
    = std::function<void(int count, const std::function<void(int)>& body)>;

// inserts[g]: the shapes of grid layer g in insertion order. Chunk axes as
// passed to ShapeGrid::init; `bands_per_layer` tile-row bands per layer
// (at most one band per tile row). Returns the number of shapes whose
// insertion the shape grid reported as failed.
int64_t insertFixedShapes(
    ShapeGrid& grid,
    const std::vector<std::vector<FixedShapeInsert>>& inserts,
    ChunkAxis x_chunks,
    ChunkAxis y_chunks,
    int bands_per_layer,
    const BandParallelFor& parallel_for);

}  // namespace drt::br
