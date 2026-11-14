// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once

#include <algorithm>

#include "frBaseTypes.h"

namespace drt {

// The four-corner spacing test separates into independent x and y minima.
// Use the endpoints, not the entire interval: an interval spanning the box
// can have both endpoints outside it.
inline frSquaredDistance minCornerAxisDistanceSquared(frCoord center,
                                                      frCoord half_width,
                                                      frCoord lower,
                                                      frCoord upper)
{
  const frCoord low_endpoint = center - half_width;
  const frCoord high_endpoint = center + half_width;
  const frCoord low_distance
      = std::max({lower - low_endpoint, low_endpoint - upper, 0});
  const frCoord high_distance
      = std::max({lower - high_endpoint, high_endpoint - upper, 0});
  const frCoord distance = std::min(low_distance, high_distance);
  return static_cast<frSquaredDistance>(distance) * distance;
}

}  // namespace drt
