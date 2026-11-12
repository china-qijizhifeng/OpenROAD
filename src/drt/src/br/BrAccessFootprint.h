// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Access footprint of the routing space (Gester et al., TODAES 18(2), 2013,
// Sec. 5.1). No drt dependencies.
//
// "threads of disjoint regions never touch the same tile chunk": the shape
// grid keeps its intervals per global routing tile (ShapeGrid::Block), the
// fast grid its runs per tile (FastGrid::TileRuns). A thread that installs
// an AccessFootprint has every access to these per-tile containers recorded
// as a tile range, so a caller can verify that a unit of work (the search
// and commit of a net) stayed inside the tiles it may touch, and measure
// how far beyond its routing area it reached. Without an installed
// footprint a record is one thread-local load and a branch.

#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace drt::br {

// Global routing tiles [xlo, xhi] x [ylo, yhi]; empty while xlo > xhi.
struct TileBox
{
  int xlo = std::numeric_limits<int>::max();
  int ylo = std::numeric_limits<int>::max();
  int xhi = std::numeric_limits<int>::min();
  int yhi = std::numeric_limits<int>::min();

  bool empty() const { return xlo > xhi || ylo > yhi; }
  void add(const TileBox& other)
  {
    xlo = std::min(xlo, other.xlo);
    ylo = std::min(ylo, other.ylo);
    xhi = std::max(xhi, other.xhi);
    yhi = std::max(yhi, other.yhi);
  }
  bool contains(const TileBox& other) const
  {
    return other.empty()
           || (other.xlo >= xlo && other.ylo >= ylo && other.xhi <= xhi
               && other.yhi <= yhi);
  }
  // Tiles by which `other` extends beyond this box (0: contained).
  int reachOf(const TileBox& other) const
  {
    if (other.empty()) {
      return 0;
    }
    return std::max({0,
                     xlo - other.xlo,
                     ylo - other.ylo,
                     other.xhi - xhi,
                     other.yhi - yhi});
  }
  friend bool operator==(const TileBox&, const TileBox&) = default;
};

enum class AccessKind : uint8_t
{
  kShapeGridRead = 0,   // ShapeGrid blocks: collect(), test run of a removal
  kShapeGridWrite = 1,  // ShapeGrid blocks: insert(), remove()
  kWideShapeRead = 2,   // wide-shape lists of collect() (locked)
  kWideShapeWrite = 3,  // wide-shape lists of insert(), remove() (locked)
  kFastGridRead = 4,    // FastGrid tiles: words and runs
  kFastGridWrite = 5,   // FastGrid tiles: assign / write
};
constexpr int kNumAccessKinds = 6;

inline const char* accessKindName(const AccessKind kind)
{
  switch (kind) {
    case AccessKind::kShapeGridRead:
      return "shape grid read";
    case AccessKind::kShapeGridWrite:
      return "shape grid write";
    case AccessKind::kWideShapeRead:
      return "wide shape read";
    case AccessKind::kWideShapeWrite:
      return "wide shape write";
    case AccessKind::kFastGridRead:
      return "fast grid read";
    case AccessKind::kFastGridWrite:
      return "fast grid write";
  }
  return "";
}

inline bool isWriteAccess(const AccessKind kind)
{
  return kind == AccessKind::kShapeGridWrite
         || kind == AccessKind::kWideShapeWrite
         || kind == AccessKind::kFastGridWrite;
}

// Wide-shape lists are read under their lock by collect() windows that may
// reach into another region's tile; the pieces found there lie beyond the
// interaction distance of the query (BrShapeGrid.h), so these reads are
// recorded and left out of the box comparison.
inline bool isCheckedAccess(const AccessKind kind)
{
  return kind != AccessKind::kWideShapeRead;
}

class AccessFootprint
{
 public:
  // Steps of a unit of work; the caller names them.
  static constexpr int kMaxSteps = 8;
  // Violations kept with their details (all are counted).
  static constexpr size_t kMaxViolations = 16;

  struct Violation
  {
    AccessKind kind = AccessKind::kShapeGridRead;
    int layer = 0;  // shape-grid layer, or track-graph layer (fast grid)
    TileBox tiles;
    TileBox allowed;
    int32_t unit = -1;
    int step = 0;
  };

  // The footprint recording the calling thread's accesses, or nullptr.
  static AccessFootprint* current() { return current_; }
  // Installs `footprint` for the calling thread (nullptr: none).
  static void install(AccessFootprint* footprint) { current_ = footprint; }

  // Begins a unit of work: clears the accessed boxes; from now on a checked
  // access outside `allowed` is a violation.
  void begin(const int32_t unit, const TileBox& allowed)
  {
    unit_ = unit;
    allowed_ = allowed;
    step_ = 0;
    reads_.fill(TileBox{});
    writes_.fill(TileBox{});
    wide_reads_ = TileBox{};
    unit_violations_ = 0;
  }
  void setStep(const int step) { step_ = std::clamp(step, 0, kMaxSteps - 1); }
  int step() const { return step_; }

  void record(const AccessKind kind, const int layer, const TileBox& tiles)
  {
    records_++;
    if (!isCheckedAccess(kind)) {
      wide_reads_.add(tiles);
      return;
    }
    (isWriteAccess(kind) ? writes_ : reads_)[step_].add(tiles);
    if (!allowed_.contains(tiles)) {
      violation_count_++;
      unit_violations_++;
      if (violations_.size() < kMaxViolations) {
        violations_.push_back({kind, layer, tiles, allowed_, unit_, step_});
      }
    }
  }

  // Tiles read / written in `step` since begin() (checked kinds).
  const TileBox& reads(const int step) const { return reads_[step]; }
  const TileBox& writes(const int step) const { return writes_[step]; }
  const TileBox& wideReads() const { return wide_reads_; }
  // All checked accesses since begin().
  TileBox accessed() const
  {
    TileBox box;
    for (int s = 0; s < kMaxSteps; s++) {
      box.add(reads_[s]);
      box.add(writes_[s]);
    }
    return box;
  }
  int32_t unit() const { return unit_; }
  const TileBox& allowed() const { return allowed_; }
  // Violations since begin().
  int64_t unitViolations() const { return unit_violations_; }
  // Since construction.
  int64_t records() const { return records_; }
  int64_t violationCount() const { return violation_count_; }
  const std::vector<Violation>& violations() const { return violations_; }

 private:
  static inline thread_local AccessFootprint* current_ = nullptr;

  int32_t unit_ = -1;
  int step_ = 0;
  TileBox allowed_;
  std::array<TileBox, kMaxSteps> reads_;
  std::array<TileBox, kMaxSteps> writes_;
  TileBox wide_reads_;
  int64_t records_ = 0;
  int64_t violation_count_ = 0;
  int64_t unit_violations_ = 0;
  std::vector<Violation> violations_;
};

// Called by the routing space at every access to a per-tile container.
inline void recordAccess(const AccessKind kind,
                         const int layer,
                         const TileBox& tiles)
{
  if (AccessFootprint* footprint = AccessFootprint::current()) [[unlikely]] {
    footprint->record(kind, layer, tiles);
  }
}

// Installs a footprint for the calling thread within a scope.
class AccessFootprintScope
{
 public:
  explicit AccessFootprintScope(AccessFootprint* footprint)
      : previous_(AccessFootprint::current())
  {
    AccessFootprint::install(footprint);
  }
  ~AccessFootprintScope() { AccessFootprint::install(previous_); }
  AccessFootprintScope(const AccessFootprintScope&) = delete;
  AccessFootprintScope& operator=(const AccessFootprintScope&) = delete;

 private:
  AccessFootprint* previous_;
};

}  // namespace drt::br
