// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Plain geometry shared by the off-track modules (blockage grid, off-track
// path search, pin access, same-net post-processing) and the routing-space
// modules (wire models, shape grid, distance rule checking, fast grid), plus
// an append-only table. Independent of odb so the modules can be tested in
// isolation.

#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <tuple>
#include <utility>

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

// --- Routing-space geometry (shape grid, distance rules, fast grid) ---

using Coord = int;

struct Point
{
  Coord x = 0;
  Coord y = 0;
  friend bool operator==(const Point&, const Point&) = default;
};

// Axis-parallel rectangle [xlo, xhi] x [ylo, yhi]. Two rectangles overlap if
// their interiors intersect; touching rectangles do not overlap but have
// distance 0.
struct Rect
{
  Coord xlo = 0;
  Coord ylo = 0;
  Coord xhi = 0;
  Coord yhi = 0;

  Coord dx() const { return xhi - xlo; }
  Coord dy() const { return yhi - ylo; }
  // Width of a rectangle: the edge length of its largest enclosed square
  // (Sec. 3.1).
  Coord width() const { return std::min(dx(), dy()); }
  Coord lo(const bool x_axis) const { return x_axis ? xlo : ylo; }
  Coord hi(const bool x_axis) const { return x_axis ? xhi : yhi; }
  Coord extent(const bool x_axis) const { return hi(x_axis) - lo(x_axis); }
  bool empty() const { return xlo >= xhi || ylo >= yhi; }
  Rect translated(const Point p) const
  {
    return {xlo + p.x, ylo + p.y, xhi + p.x, yhi + p.y};
  }
  Rect bloated(const Coord bx, const Coord by) const
  {
    return {xlo - bx, ylo - by, xhi + bx, yhi + by};
  }
  bool overlaps(const Rect& o) const
  {
    return xlo < o.xhi && o.xlo < xhi && ylo < o.yhi && o.ylo < yhi;
  }
  bool contains(const Point p) const
  {
    return p.x >= xlo && p.x <= xhi && p.y >= ylo && p.y <= yhi;
  }
  Rect hull(const Rect& o) const
  {
    return {std::min(xlo, o.xlo),
            std::min(ylo, o.ylo),
            std::max(xhi, o.xhi),
            std::max(yhi, o.yhi)};
  }
  // Transposed rectangle (x <-> y).
  Rect transposed() const { return {ylo, xlo, yhi, xhi}; }
  friend bool operator==(const Rect&, const Rect&) = default;
};

// Floor division for possibly negative numerators.
inline int64_t floorDivide(const int64_t a, const int64_t b)
{
  const int64_t q = a / b;
  return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}

// Largest g >= 0 with g * g < value, or -1 if value <= 0.
inline int64_t largestBelowSquare(const int64_t value)
{
  if (value <= 0) {
    return -1;
  }
  int64_t g = static_cast<int64_t>(std::sqrt(static_cast<double>(value)));
  while (g > 0 && g * g >= value) {
    g--;
  }
  while ((g + 1) * (g + 1) < value) {
    g++;
  }
  return g;
}

// Append-only table with stable element addresses. Appends are serialised by
// the caller's lock (see append()); readers may read any index that was
// published to them through that lock or a later synchronisation point
// without locking. Capacity is kBlocks * kBlockSize elements.
template <typename T>
class AppendOnlyTable
{
 public:
  static constexpr size_t kBlockSize = 4096;
  static constexpr size_t kBlocks = 4096;

  AppendOnlyTable() = default;
  AppendOnlyTable(const AppendOnlyTable&) = delete;
  AppendOnlyTable& operator=(const AppendOnlyTable&) = delete;

  // Requires the caller to hold the table's writer lock.
  size_t append(T&& value)
  {
    const size_t n = size_.load(std::memory_order_relaxed);
    const size_t block = n / kBlockSize;
    if (block >= kBlocks) {
      throw std::length_error("AppendOnlyTable capacity exceeded");
    }
    if (!blocks_[block]) {
      blocks_[block] = std::make_unique<T[]>(kBlockSize);
    }
    blocks_[block][n % kBlockSize] = std::move(value);
    size_.store(n + 1, std::memory_order_release);
    return n;
  }
  const T& operator[](const size_t index) const
  {
    return blocks_[index / kBlockSize][index % kBlockSize];
  }
  size_t size() const { return size_.load(std::memory_order_acquire); }
  void clear()
  {
    for (auto& b : blocks_) {
      b.reset();
    }
    size_.store(0, std::memory_order_release);
  }

 private:
  std::array<std::unique_ptr<T[]>, kBlocks> blocks_;
  std::atomic<size_t> size_{0};
};

}  // namespace drt::br
