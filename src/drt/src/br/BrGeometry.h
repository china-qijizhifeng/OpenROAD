// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Plain geometry and an append-only table shared by the BonnRoute routing
// space modules (wire models, shape grid, distance rule checking, fast grid).
// No drt dependencies.

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

namespace drt::br {

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
inline int64_t floorDiv(const int64_t a, const int64_t b)
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
