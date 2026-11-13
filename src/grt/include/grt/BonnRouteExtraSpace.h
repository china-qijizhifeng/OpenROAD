// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Extra space s(n, e) of BonnRoute global routing (Gester et al., TODAES
// 2013, Sec. 2.1) handed to detailed routing.
//
// After `global_route -bonnroute` with `set_bonnroute_param extra_space 1`,
// every net that received extra space on some global routing edge carries
// an odb::dbStringProperty named kExtraSpaceProperty whose value
// encodeExtraSpace() produced. Nets without the property have no extra
// space (global routing deletes stale properties of all nets it routes).
//
// An edge is identified by the routing level of its layer
// (odb::dbTechLayer::getRoutingLevel(), 1-based) and the indices (x, y) of
// its lower gcell in the block's dbGCellGrid: on a horizontal layer it
// joins gcells (x, y) and (x + 1, y), on a vertical layer (x, y) and
// (x, y + 1). `space` is the extra space in routing tracks of that layer,
// beyond the net's own wire width and minimum spacing, to be split between
// both sides of the wire. Global routing has already checked it against the
// edge capacity left after all nets' wires.
//
// Header-only and free of odb so that drt can include it alone.

#pragma once

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace grt::bonnroute {

inline constexpr const char* kExtraSpaceProperty = "bonnroute_extra_space";

struct ExtraSpaceEdge
{
  int level = 0;  // routing level, 1-based
  int x = 0;
  int y = 0;
  float space = 0.0f;  // tracks

  bool operator==(const ExtraSpaceEdge& other) const = default;
};

// "v1" followed by " level,x,y,space" per edge (space with 3 decimals).
inline std::string encodeExtraSpace(const std::vector<ExtraSpaceEdge>& edges)
{
  std::string out = "v1";
  char buffer[64];
  for (const ExtraSpaceEdge& e : edges) {
    std::snprintf(buffer,
                  sizeof(buffer),
                  " %d,%d,%d,%.3f",
                  e.level,
                  e.x,
                  e.y,
                  static_cast<double>(e.space));
    out += buffer;
  }
  return out;
}

// Returns false (and leaves `edges` empty) for an unknown format.
inline bool decodeExtraSpace(const std::string& text,
                             std::vector<ExtraSpaceEdge>& edges)
{
  edges.clear();
  if (text.rfind("v1", 0) != 0) {
    return false;
  }
  const char* p = text.c_str() + 2;
  while (*p == ' ') {
    ExtraSpaceEdge e;
    char* end = nullptr;
    e.level = static_cast<int>(std::strtol(p + 1, &end, 10));
    if (*end != ',') {
      edges.clear();
      return false;
    }
    e.x = static_cast<int>(std::strtol(end + 1, &end, 10));
    if (*end != ',') {
      edges.clear();
      return false;
    }
    e.y = static_cast<int>(std::strtol(end + 1, &end, 10));
    if (*end != ',') {
      edges.clear();
      return false;
    }
    e.space = std::strtof(end + 1, &end);
    edges.push_back(e);
    p = end;
  }
  if (*p != '\0') {
    edges.clear();
    return false;
  }
  return true;
}

}  // namespace grt::bonnroute
