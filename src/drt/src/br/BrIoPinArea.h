// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <optional>
#include <vector>

#include "br/BrGeometry.h"
#include "br/BrSameNetRules.h"

namespace drt::br {

// Routed metal on a top-level (IO) pin (BonnRouteDRParams::io_pin_area):
// sign-off keeps IO pin shapes on the pin purpose, so the routed metal of
// the net that touches the pin must meet the layer's minimum area by itself.

// Indices of the pieces of `metal` that touch `pin` and of every piece
// connected to them through touching pieces (one polygon without the pin).
inline std::vector<size_t> piecesOnPin(const std::vector<Box>& metal,
                                       const Box& pin)
{
  std::vector<char> in(metal.size(), 0);
  std::vector<size_t> group;
  for (size_t k = 0; k < metal.size(); k++) {
    if (metal[k].touches(pin)) {
      in[k] = 1;
      group.push_back(k);
    }
  }
  for (size_t head = 0; head < group.size(); head++) {
    const Box& g = metal[group[head]];
    for (size_t k = 0; k < metal.size(); k++) {
      if (!in[k] && metal[k].touches(g)) {
        in[k] = 1;
        group.push_back(k);
      }
    }
  }
  std::sort(group.begin(), group.end());
  return group;
}

struct IoPinPatchRequest
{
  std::vector<Box> group;  // the routed metal on the pin (piecesOnPin)
  GridPoint anchor;        // a via origin or wire end of the net on the pin
  bool horizontal = true;  // preferred direction of the layer
  GeoCoord half_width = 0;
  int64_t min_area = 0;
  GeoCoord grid = 1;  // manufacturing grid
  Box die;
  // True if new metal of the net covering the box is legal against other
  // nets.
  std::function<bool(const Box&)> legal;
};

// A patch from the anchor in the preferred direction, toward the die inside
// first (away from the nearer die edge), whose union with `group` meets
// min_area; lengths on the manufacturing grid. None when the group already
// meets it or no legal patch inside the die does.
inline std::optional<Box> ioPinPatch(const IoPinPatchRequest& req)
{
  int64_t area = unionArea(req.group);
  if (req.min_area <= 0 || area >= req.min_area) {
    return std::nullopt;
  }
  const GridPoint& p = req.anchor;
  const int64_t lo = req.horizontal ? int64_t{p.x} - req.die.xlo
                                    : int64_t{p.y} - req.die.ylo;
  const int64_t hi = req.horizontal ? int64_t{req.die.xhi} - p.x
                                    : int64_t{req.die.yhi} - p.y;
  const std::array<int, 2> directions
      = lo <= hi ? std::array<int, 2>{1, -1} : std::array<int, 2>{-1, 1};
  const int64_t width = std::max<int64_t>(1, 2 * int64_t{req.half_width});
  const int64_t grid = std::max<GeoCoord>(1, req.grid);
  for (const int dir : directions) {
    int64_t need = (req.min_area - area + width - 1) / width;
    for (int attempt = 0; attempt < 4 && need > 0; attempt++) {
      need = (need + grid - 1) / grid * grid;
      const GridPoint tip
          = req.horizontal
                ? GridPoint{static_cast<GeoCoord>(p.x + dir * need), p.y}
                : GridPoint{p.x, static_cast<GeoCoord>(p.y + dir * need)};
      const Box ext{std::min(p.x, tip.x) - req.half_width,
                    std::min(p.y, tip.y) - req.half_width,
                    std::max(p.x, tip.x) + req.half_width,
                    std::max(p.y, tip.y) + req.half_width};
      if (ext.xlo < req.die.xlo || ext.ylo < req.die.ylo
          || ext.xhi > req.die.xhi || ext.yhi > req.die.yhi
          || (req.legal && !req.legal(ext))) {
        break;
      }
      std::vector<Box> with = req.group;
      with.push_back(ext);
      const int64_t total = unionArea(with);
      if (total >= req.min_area) {
        return ext;
      }
      need += (req.min_area - total + width - 1) / width;
    }
  }
  return std::nullopt;
}

}  // namespace drt::br
