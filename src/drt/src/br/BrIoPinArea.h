// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <optional>
#include <utility>
#include <vector>

#include "br/BrGeometry.h"
#include "br/BrSameNetRules.h"

namespace drt::br {

// Routed metal on a top-level (IO) pin (BonnRouteDRParams::io_pin_area):
// sign-off keeps IO pin shapes on the pin purpose, so the routed metal of
// the net that touches the pin must meet the layer's minimum area by itself.

// Connected metal components that touch the pin. The pin purpose does not
// join disjoint routed components for minimum-area sign-off.
inline std::vector<std::vector<size_t>> metalOnPinComponents(
    const std::vector<Box>& metal,
    const Box& pin)
{
  std::vector<char> visited(metal.size(), 0);
  std::vector<std::vector<size_t>> groups;
  for (size_t seed = 0; seed < metal.size(); seed++) {
    if (visited[seed] || !metal[seed].touches(pin)) {
      continue;
    }
    std::vector<size_t> group{seed};
    visited[seed] = 1;
    for (size_t head = 0; head < group.size(); head++) {
      for (size_t k = 0; k < metal.size(); k++) {
        if (!visited[k] && metal[k].touches(metal[group[head]])) {
          visited[k] = 1;
          group.push_back(k);
        }
      }
    }
    std::sort(group.begin(), group.end());
    groups.push_back(std::move(group));
  }
  return groups;
}

// A wire on the pin layer that adds no metal beyond the pin and the net's
// via pads there: FlexDR's connectivity check removes it (a pin-layer
// segment next to a via whose pad already touches the pin), so it does not
// count toward the routed metal's area.
inline bool coveredByPinAndPads(const Box& wire,
                                const Box& pin,
                                const std::vector<Box>& pads)
{
  std::vector<Box> base = pads;
  base.push_back(pin);
  std::vector<Box> with = base;
  with.push_back(wire);
  return unionArea(with) == unionArea(base);
}

// Space between disjoint same-net metal still has to meet the layer's
// rule. Other-net legality deliberately ignores this net's own shapes.
inline bool patchLeavesNarrowGap(const Box& patch,
                                 const Box& other,
                                 const GeoCoord spacing)
{
  if (patch.touches(other)) {
    return false;
  }
  const int64_t dx = std::max<int64_t>(
      {other.xlo - int64_t{patch.xhi}, int64_t{patch.xlo} - other.xhi, 0});
  const int64_t dy = std::max<int64_t>(
      {other.ylo - int64_t{patch.yhi}, int64_t{patch.ylo} - other.yhi, 0});
  return dx * dx + dy * dy < int64_t{spacing} * spacing;
}

struct IoPinPatchRequest
{
  std::vector<Box> group;  // one routed metal component touching the pin
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
  // A narrow component near another net may have enough room only when
  // the extension is shared between both sides. Try the centred patch
  // after preserving the existing one-sided choices. Move it along the
  // preferred direction when needed to stay inside the die.
  int64_t need = (req.min_area - area + width - 1) / width;
  for (int attempt = 0; attempt < 4 && need > 0; attempt++) {
    need = (need + grid - 1) / grid * grid;
    const int64_t back = need / (2 * grid) * grid;
    const int64_t position = req.horizontal ? p.x : p.y;
    int64_t first = position - back - req.half_width;
    int64_t last = position + need - back + req.half_width;
    const int64_t die_lo = req.horizontal ? req.die.xlo : req.die.ylo;
    const int64_t die_hi = req.horizontal ? req.die.xhi : req.die.yhi;
    if (first < die_lo) {
      const int64_t shift = (die_lo - first + grid - 1) / grid * grid;
      first += shift;
      last += shift;
    }
    if (last > die_hi) {
      const int64_t shift = (last - die_hi + grid - 1) / grid * grid;
      first -= shift;
      last -= shift;
    }
    const Box ext = req.horizontal ? Box{static_cast<GeoCoord>(first),
                                         p.y - req.half_width,
                                         static_cast<GeoCoord>(last),
                                         p.y + req.half_width}
                                   : Box{p.x - req.half_width,
                                         static_cast<GeoCoord>(first),
                                         p.x + req.half_width,
                                         static_cast<GeoCoord>(last)};
    if (ext.xlo < req.die.xlo || ext.ylo < req.die.ylo || ext.xhi > req.die.xhi
        || ext.yhi > req.die.yhi || first > position || last < position) {
      break;
    }
    std::vector<Box> with = req.group;
    with.push_back(ext);
    const int64_t total = unionArea(with);
    if (total < req.min_area) {
      need += (req.min_area - total + width - 1) / width;
      continue;
    }
    if (!req.legal || req.legal(ext)) {
      return ext;
    }
    break;
  }
  return std::nullopt;
}

}  // namespace drt::br
