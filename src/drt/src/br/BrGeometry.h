// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Plain geometry shared by the off-track modules (blockage grid, off-track
// path search, pin access, same-net post-processing). Independent of odb so
// the modules can be tested in isolation.

#pragma once

#include <algorithm>
#include <cstdint>
#include <tuple>

namespace drt::br {

using GeoCoord = int;

struct GridPoint
{
  GeoCoord x = 0;
  GeoCoord y = 0;
  bool operator==(const GridPoint& o) const { return x == o.x && y == o.y; }
  bool operator<(const GridPoint& o) const
  {
    return std::tie(x, y) < std::tie(o.x, o.y);
  }
};

// A node of a routed path: consecutive nodes on the same layer form a wire
// segment, consecutive nodes at the same point on adjacent layers a via.
// The layer numbering is the caller's.
struct RouteNode
{
  int layer = 0;
  GridPoint p;
  bool operator==(const RouteNode& o) const
  {
    return layer == o.layer && p == o.p;
  }
};

// Closed axis-parallel box [xlo, xhi] x [ylo, yhi]. Obstacles of the
// off-track search are the open interiors of boxes (Sec. 3.8: a path is
// tau-feasible if it "does not intersect the interior of any obstacle").
struct Box
{
  GeoCoord xlo = 0;
  GeoCoord ylo = 0;
  GeoCoord xhi = 0;
  GeoCoord yhi = 0;

  bool operator==(const Box& o) const
  {
    return xlo == o.xlo && ylo == o.ylo && xhi == o.xhi && yhi == o.yhi;
  }
  bool operator<(const Box& o) const
  {
    return std::tie(xlo, ylo, xhi, yhi) < std::tie(o.xlo, o.ylo, o.xhi, o.yhi);
  }
  int64_t width() const { return int64_t{xhi} - xlo; }
  int64_t height() const { return int64_t{yhi} - ylo; }
  int64_t area() const { return width() * height(); }
  bool containsClosed(GeoCoord x, GeoCoord y) const
  {
    return x >= xlo && x <= xhi && y >= ylo && y <= yhi;
  }
  bool containsOpen(GeoCoord x, GeoCoord y) const
  {
    return x > xlo && x < xhi && y > ylo && y < yhi;
  }
  // True if the closed boxes share a point.
  bool touches(const Box& o) const
  {
    return xlo <= o.xhi && o.xlo <= xhi && ylo <= o.yhi && o.ylo <= yhi;
  }
  Box bloated(GeoCoord d) const { return {xlo - d, ylo - d, xhi + d, yhi + d}; }
  Box merged(const Box& o) const
  {
    return {std::min(xlo, o.xlo),
            std::min(ylo, o.ylo),
            std::max(xhi, o.xhi),
            std::max(yhi, o.yhi)};
  }
  // Intersection; may be inverted when the boxes are disjoint.
  Box clipped(const Box& o) const
  {
    return {std::max(xlo, o.xlo),
            std::max(ylo, o.ylo),
            std::min(xhi, o.xhi),
            std::min(yhi, o.yhi)};
  }
  bool inverted() const { return xlo > xhi || ylo > yhi; }
};

// True if the closed segment (a, b) (axis-parallel or a point) meets the
// open interior of `box`.
inline bool segmentMeetsInterior(const GridPoint& a,
                                 const GridPoint& b,
                                 const Box& box)
{
  const GeoCoord xlo = std::min(a.x, b.x);
  const GeoCoord xhi = std::max(a.x, b.x);
  const GeoCoord ylo = std::min(a.y, b.y);
  const GeoCoord yhi = std::max(a.y, b.y);
  return xhi > box.xlo && xlo < box.xhi && yhi > box.ylo && ylo < box.yhi;
}

// Isometry of the integer lattice: p -> M p + offset with M one of the
// eight signed permutation matrices (rotations by multiples of 90 degrees
// and mirrors), as used for circuit classes (Sec. 4.3: "up to translation,
// mirroring, and rotation").
struct Isometry
{
  int m00 = 1;
  int m01 = 0;
  int m10 = 0;
  int m11 = 1;
  GeoCoord dx = 0;
  GeoCoord dy = 0;

  GridPoint apply(const GridPoint& p) const
  {
    return {m00 * p.x + m01 * p.y + dx, m10 * p.x + m11 * p.y + dy};
  }
  Box apply(const Box& b) const
  {
    const GridPoint p = apply(GridPoint{b.xlo, b.ylo});
    const GridPoint q = apply(GridPoint{b.xhi, b.yhi});
    return {std::min(p.x, q.x),
            std::min(p.y, q.y),
            std::max(p.x, q.x),
            std::max(p.y, q.y)};
  }
  // The linear part only (no translation).
  Isometry linear() const { return {m00, m01, m10, m11, 0, 0}; }
  // Signed permutation matrices are orthogonal: the inverse is the
  // transpose.
  Isometry inverse() const
  {
    Isometry r{m00, m10, m01, m11, 0, 0};
    const GridPoint o = r.apply(GridPoint{dx, dy});
    r.dx = -o.x;
    r.dy = -o.y;
    return r;
  }
  // (this o other)(p) = this(other(p)).
  Isometry compose(const Isometry& o) const
  {
    Isometry r;
    r.m00 = m00 * o.m00 + m01 * o.m10;
    r.m01 = m00 * o.m01 + m01 * o.m11;
    r.m10 = m10 * o.m00 + m11 * o.m10;
    r.m11 = m10 * o.m01 + m11 * o.m11;
    const GridPoint t = linear().apply(GridPoint{o.dx, o.dy});
    r.dx = t.x + dx;
    r.dy = t.y + dy;
    return r;
  }
  // True if the image of a horizontal line is vertical.
  bool swapsAxes() const { return m00 == 0; }
  bool operator==(const Isometry& o) const
  {
    return m00 == o.m00 && m01 == o.m01 && m10 == o.m10 && m11 == o.m11
           && dx == o.dx && dy == o.dy;
  }
};

// All eight linear isometries (dihedral group of the square).
inline Isometry dihedral(int k)
{
  static constexpr int kMats[8][4] = {{1, 0, 0, 1},
                                      {0, -1, 1, 0},
                                      {-1, 0, 0, -1},
                                      {0, 1, -1, 0},
                                      {-1, 0, 0, 1},
                                      {1, 0, 0, -1},
                                      {0, 1, 1, 0},
                                      {0, -1, -1, 0}};
  const int* m = kMats[k & 7];
  return {m[0], m[1], m[2], m[3], 0, 0};
}

}  // namespace drt::br
