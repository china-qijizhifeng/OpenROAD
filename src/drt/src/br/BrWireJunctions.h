// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once

#include <span>
#include <vector>

#include "br/BrGeometry.h"

namespace drt::br {
struct WireSegment
{
  int layer;
  GridPoint begin;
  GridPoint end;
};

// Strict interior crossings between perpendicular wire centerlines on one
// layer. Split and truncate both wires here before FlexDR merges overlaps:
// its connectivity graph joins wire endpoints to other wires, not two
// interiors. Endpoint contacts already survive that graph unchanged.
// Each wire receives sorted unique cuts; input order and direction are kept.
std::vector<std::vector<GridPoint>> interiorWireCrossings(
    std::span<const WireSegment> wires);
}  // namespace drt::br
