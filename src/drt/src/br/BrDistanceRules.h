// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Distance rule checking module of BonnRoute (Gester et al., TODAES 18(2),
// 2013, Sec. 3.4). No drt dependencies.
//
// "It gets a location l and a set W of wiretypes as input and returns the
// subset of wire- and via models of wiretypes in W that are allowed to be
// placed at l. [...] it computes this information by querying all
// intervals I from the shape grid that may contain shapes with distance
// rule violations to any of the models placed at location l. In addition,
// it returns a maximal horizontal and a maximal vertical interval of
// locations, each containing l, in which the same subset of wire and via
// models is allowed, and for which this can be determined without querying
// additional intervals I from the shape grid. Furthermore, the distance
// rule checking module can report nets whose (partial) removal increases
// the set of allowed wire or via models at the given location."
//
// The module answers only through ShapeGrid::collect(). Queries are const
// and thread-safe; scratch memory is per thread.

#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrWireModel.h"

namespace drt::br {

// A model to be placed at the query location: rect relative to the
// location on a shape-grid layer. Models with the same group form one
// answer (e.g. the shapes of a via model); group -1 = own index.
struct LocatedModel
{
  int layer = 0;
  Rect rect;
  ClassId cls = 0;
  uint8_t line_end_ends = 3;
  int group = -1;
};

struct CheckOptions
{
  // Answer = per group the lowest ripup level at which it becomes legal
  // (fast grid encoding, Sec. 3.6); false: answer = allowed or not.
  bool levels = true;
  bool report_nets = false;
  // Which maximal intervals to compute (horizontal / vertical line).
  bool intervals = true;
  bool x_interval = true;
  bool y_interval = true;
  // Shapes owned by this net count as removed (Sec. 4.4: shapes of the
  // connected components are removed temporarily).
  int32_t ignored_owner = kNoNet;
  // Additional area queried around the models to lengthen the intervals.
  Coord margin_x = 0;
  Coord margin_y = 0;
  // Minimum number of answer groups (groups without models are legal).
  int groups = 0;
};

struct BlockingNet
{
  int32_t net = kNoNet;
  RipupLevel level = 0;  // highest ripup level of its violating shapes
  uint64_t groups = 0;   // answer groups it blocks (bit per group < 64)
  friend bool operator==(const BlockingNet&, const BlockingNet&) = default;
};

struct CheckResult
{
  std::vector<RipupLevel> levels;  // per group; 0 = allowed
  // Maximal intervals of locations (inclusive) with the same answer.
  Coord x_lo = 0;
  Coord x_hi = 0;
  Coord y_lo = 0;
  Coord y_hi = 0;
  // The interval end is the end of the area known from the queried
  // intervals (otherwise the answer changes just beyond it).
  bool x_lo_known_bound = false;
  bool x_hi_known_bound = false;
  bool y_lo_known_bound = false;
  bool y_hi_known_bound = false;
  // Nets whose removal makes a blocked group legal, sorted by net; groups
  // blocked by non-removable shapes are not reported.
  std::vector<BlockingNet> nets;
  int64_t intervals = 0;  // shape-grid intervals queried
};

// Paper interface answer for one wire type at a wiring-layer location.
struct WireTypeAnswer
{
  RipupLevel preferred = 0;
  RipupLevel non_preferred = 0;
  RipupLevel via_up = 0;    // via model of via layer z (all its shapes)
  RipupLevel via_down = 0;  // via model of via layer z - 1
};

class DistanceRuleChecker
{
 public:
  DistanceRuleChecker() = default;
  explicit DistanceRuleChecker(const ShapeGrid* grid) : grid_(grid) {}
  void init(const ShapeGrid* grid) { grid_ = grid; }
  const ShapeGrid& grid() const { return *grid_; }

  void check(Point location,
             std::span<const LocatedModel> models,
             const CheckOptions& options,
             CheckResult& result) const;

  // query(l, W) of Sec. 3.4 at wiring layer z: preferred and non-preferred
  // wire model and the via models up and down of each wire type. The
  // result's levels are per (wire type, answer) in this order.
  void query(int z,
             Point location,
             std::span<const WireType* const> wire_types,
             const CheckOptions& options,
             std::vector<WireTypeAnswer>& answers,
             CheckResult& result) const;

  // Level of one placed shape (absolute rect), with the reference pairwise
  // rule on the stored pieces; optional blocking-net report.
  RipupLevel shapeLevel(int layer,
                        const Rect& rect,
                        ClassId cls,
                        uint8_t line_end_ends,
                        int32_t ignored_owner,
                        std::vector<BlockingNet>* nets = nullptr) const;

 private:
  const ShapeGrid* grid_ = nullptr;
};

// The models of `type` at wiring layer z, in WireTypeAnswer order as
// groups base + 0..3; invalid models are left out.
void appendWireTypeModels(const WireType& type,
                          int z,
                          int num_wiring_layers,
                          int base_group,
                          std::vector<LocatedModel>& models);

}  // namespace drt::br
