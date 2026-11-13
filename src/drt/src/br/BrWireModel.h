// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Diff-net distance rules, shape classes, wire models, via models and wire
// types of BonnRoute (Gester et al., TODAES 18(2), 2013, Sec. 3.1 and 3.2).
// No drt dependencies; BonnRouteDR fills the rules from the drt tech.
//
// Shape-grid layers: wiring layer z of the track graph is layer 2z, the via
// layer between wiring layers z and z + 1 is layer 2z + 1.

#pragma once

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "br/BrGeometry.h"

// Shape classes carry per-layer reaches from the interned classes and wide
// classes (ShapeClassTable::narrowReach, wide).
#define BR_SHAPE_CLASS_REACHES 1

namespace drt::br {

inline int wiringGridLayer(const int z)
{
  return 2 * z;
}
inline int viaGridLayer(const int z)
{
  return 2 * z + 1;
}

// Sec. 3.1: "The required distance between two shapes can be described as a
// nondecreasing function of their widths and common run-length."
// A table over width rows and run-length columns. Row k holds widths in
// [width_breaks[k], width_breaks[k + 1]), column c run-lengths in
// [prl_breaks[c], prl_breaks[c + 1]); column 0 is used for run-length <= 0.
class SpacingTable
{
 public:
  SpacingTable() = default;
  static SpacingTable constant(Coord spacing);
  // Samples f(width_a, width_b, prl) at the breakpoints (each break point
  // is a representative of its bucket; callers pass t and t + 1 for every
  // threshold t of a rule so both strict and non-strict thresholds are
  // bucket borders). The values are closed under maximum over smaller rows
  // and columns and symmetrised, which makes the table nondecreasing as
  // Sec. 3.1 states; a non-monotone rule is thereby modelled conservatively.
  template <typename F>
  static SpacingTable sample(std::vector<Coord> width_breaks,
                             std::vector<Coord> prl_breaks,
                             F&& f);
  int widthRow(Coord width) const;
  int prlColumn(Coord prl) const;
  int numRows() const { return static_cast<int>(width_breaks_.size()); }
  int numColumns() const { return static_cast<int>(prl_breaks_.size()); }
  Coord at(int row_a, int row_b, int column) const
  {
    return values_[(static_cast<size_t>(row_a) * numRows() + row_b)
                       * numColumns()
                   + column];
  }
  Coord maxValue() const;
  Coord minValue() const;
  const std::vector<Coord>& widthBreaks() const { return width_breaks_; }
  const std::vector<Coord>& prlBreaks() const { return prl_breaks_; }

 private:
  void close();

  std::vector<Coord> width_breaks_{0};
  std::vector<Coord> prl_breaks_{0};
  std::vector<Coord> values_{0};
};

template <typename F>
SpacingTable SpacingTable::sample(std::vector<Coord> width_breaks,
                                  std::vector<Coord> prl_breaks,
                                  F&& f)
{
  auto normalise = [](std::vector<Coord>& v) {
    v.push_back(0);
    std::erase_if(v, [](Coord c) { return c < 0; });
    std::ranges::sort(v);
    v.erase(std::ranges::unique(v).begin(), v.end());
  };
  normalise(width_breaks);
  normalise(prl_breaks);
  SpacingTable t;
  t.width_breaks_ = std::move(width_breaks);
  t.prl_breaks_ = std::move(prl_breaks);
  const int rows = t.numRows();
  const int cols = t.numColumns();
  t.values_.assign(static_cast<size_t>(rows) * rows * cols, 0);
  for (int a = 0; a < rows; a++) {
    for (int b = 0; b < rows; b++) {
      for (int c = 0; c < cols; c++) {
        t.values_[(static_cast<size_t>(a) * rows + b) * cols + c] = std::max(
            0,
            static_cast<Coord>(
                f(t.width_breaks_[a], t.width_breaks_[b], t.prl_breaks_[c])));
      }
    }
  }
  t.close();
  return t;
}

// Line-end rule (Sec. 3.1): an edge shorter than `max_edge` between two
// convex corners requires `space` in front of it, within `within` beyond its
// sides. Inactive if max_edge or space is 0.
struct LineEndRule
{
  Coord max_edge = 0;
  Coord space = 0;
  Coord within = 0;
  bool active() const { return max_edge > 0 && space > 0; }
};

// LEF58 rules that relate the sides and corners of two shapes
// (Lef58Rules.md). Each is evaluated from one shape's end-of-line sides:
// sides shorter than LayerRules::eol_widths[eol] (drt: an edge between two
// convex corners shorter than the rule's EOL width).
enum class SideRuleKind : uint8_t
{
  // EOLKEEPOUT: the other shape overlaps the keep-out region of the side
  // (`forward` beyond it, `backward` behind it, `side` beyond its ends).
  kRegionOverlap,
  // EOLKEEPOUT CORNERONLY: a corner of the other shape lies strictly
  // inside the keep-out region.
  kRegionCorner,
  // ENDOFLINE spacing: a side of the other shape facing the end-of-line
  // side is closer than `forward` (and not touching it) and overlaps it
  // laterally within `side`.
  kFacing,
  // ENDTOEND: as kFacing, but the facing side is an end-of-line side of the
  // other shape (same EOL width); `side` = 2 * the ENDTOEND extension.
  kFacingEol,
};

struct SideRule
{
  SideRuleKind kind = SideRuleKind::kRegionOverlap;
  uint8_t eol = 0;  // index into LayerRules::eol_widths
  Coord forward = 0;
  Coord backward = 0;
  Coord side = 0;
  friend bool operator==(const SideRule&, const SideRule&) = default;
};

// LEF58 cut class of a via layer (width = the smaller cut extent).
struct CutClassDims
{
  Coord width = 0;
  Coord length = 0;
  friend bool operator==(const CutClassDims&, const CutClassDims&) = default;
};

constexpr int kMaxEolWidths = 4;

// Rules of one shape-grid layer.
struct LayerRules
{
  bool via_layer = false;
  // Preferred routing direction; a via layer uses the direction of the
  // next lower wiring layer (Sec. 3.6).
  bool horizontal = true;
  // Wiring layer: metal-metal spacing; via layer: cut-cut spacing.
  SpacingTable spacing;
  LineEndRule line_end;
  // Via layer: inter-layer via rule (Sec. 3.1, 3.2) between this layer's
  // cuts and projections of cuts of the next lower via layer onto it; 0
  // means no rule (and no projections are stored).
  Coord projection_spacing = 0;
  // Wiring layer: minimum spacing of the nondefault rules (Sec. 3.2 wire
  // types); a metal class with ndr = k > 0 keeps at least ndr_spacing[k - 1]
  // to every other metal shape, as drt applies NDR spacing to the routed
  // shapes of NDR nets (max with the table value). 0: no NDR spacing.
  std::vector<Coord> ndr_spacing;
  Coord ndrSpacing(const int ndr) const
  {
    return ndr > 0 && ndr <= static_cast<int>(ndr_spacing.size())
               ? ndr_spacing[ndr - 1]
               : 0;
  }
  // Wiring layer: side rules (LEF58 EOL spacing, ENDTOEND, EOLKEEPOUT) with
  // their EOL widths (at most kMaxEolWidths), and the convex corner
  // spacing of LEF58_CORNERSPACING (largest x / y gap of two facing convex
  // corners, Chebyshev; 0: none). See Lef58Rules.md.
  std::vector<Coord> eol_widths;
  std::vector<SideRule> side_rules;
  Coord corner_spacing = 0;
  // Every routed preferred-direction wire ends inside a via pad of its wire
  // type (the pads cover the wire's end square), so its ends are no polygon
  // sides and its corners no polygon corners: the pads carry the side
  // rules. Otherwise preferred wires carry them at their ends.
  bool wire_ends_covered = false;
  // Via layer: LEF58 cut classes and the cut spacing table between them
  // ((n + 1)^2 entries, index 0 = no class, i.e. the table's DEFAULT),
  // edge to edge (a centre-to-centre entry is modelled edge to edge,
  // conservatively); the requirement is the max with `spacing`.
  std::vector<CutClassDims> cut_classes;
  std::vector<Coord> cut_class_spacing;
  bool sideRules() const { return !side_rules.empty() || corner_spacing > 0; }
  // 1 + index of the cut class of a dx x dy cut, 0 if none matches.
  int cutClassOf(Coord dx, Coord dy) const;
  Coord cutClassSpacing(int a, int b) const;
  // Largest distance at which side or corner rules relate two shapes.
  Coord sideReach() const;
  // Largest distance at which two shapes of this layer can violate a rule.
  Coord reach() const;
  // Smallest distance required between any two shapes of different nets on
  // this layer, ignoring projections (shape grid cell size, Sec. 3.3).
  Coord minimumSpacing() const;
};

enum class ShapeClassKind : uint8_t
{
  kMetal = 0,
  kCut = 1,
  kCutProjection = 2,
};

// Run-length column marking a preferred-direction wire: its run-length
// along the preferred direction is treated as unbounded (Sec. 3.1, rules
// whose run-lengths "are often not known during the construction of a
// route" are modelled conservatively).
constexpr uint8_t kUnboundedRun = 0xff;

// Sec. 3.2: a shape class is "used to determine the minimum distance
// requirements" of a shape. Two shapes of the same class have the same
// requirements to every other class.
struct ShapeClass
{
  int16_t layer = 0;
  ShapeClassKind kind = ShapeClassKind::kMetal;
  uint8_t width_row = 0;
  // Spacing-table column of the shape's own extent in x / y; the common
  // run-length of two shapes is bounded by the smaller of their extents.
  uint8_t run_column_x = 0;
  uint8_t run_column_y = 0;
  // Sec. 3.1: wires are assumed to have line-ends in preferred direction.
  bool line_end = false;
  // Nondefault rule of a routed metal shape (LayerRules::ndr_spacing); 0
  // for the default rule, pins, blockages and cuts.
  uint8_t ndr = 0;
  // Cut: 1 + LEF58 cut class (LayerRules::cut_classes), 0 = none.
  uint8_t cut_class = 0;
  // Side rules: its sides (bit per kSide*) that are polygon sides; its
  // rectangle corners between two such sides are polygon corners. A
  // preferred wire whose ends lie inside via pads keeps its flanks only
  // (LayerRules::wire_ends_covered).
  uint8_t sides_real = 0;
  // Side rules: per EOL width k, bits 4k..4k+3 = its sides (xlo, xhi, ylo,
  // yhi) that are end-of-line sides (shorter than the width).
  uint16_t eol_sides = 0;
  friend bool operator==(const ShapeClass&, const ShapeClass&) = default;
};

using ClassId = uint16_t;

struct ShapeClassHash
{
  size_t operator()(const ShapeClass& c) const;
};

// How a shape is classified (Sec. 3.1 modelling choices).
enum class ShapeRole : uint8_t
{
  kPreferredWire,     // line-end in preferred direction, unbounded run
  kNonPreferredWire,  // jog: no line-end ("optimistic", Sec. 3.1)
  kViaPad,            // line-end in preferred direction, actual run
  kCut,
  kCutProjection,
  kFixed,  // blockages and pins: actual geometry, no line-end
  // A preferred-direction wire whose ends are polygon sides even where
  // routed wires end in via pads (side rules: patches and off-track
  // pieces, Lef58Rules.md); otherwise as kPreferredWire.
  kRealWire,
};

// Minimum distance requirement between two shape classes (Sec. 3.1): the
// l2 distance must be at least `facing_x` for shapes separated in x whose
// y-projections overlap, `facing_y` for shapes separated in y whose
// x-projections overlap and `diagonal` otherwise; overlapping interiors
// violate if `overlap`.
struct Requirement
{
  Coord facing_x = 0;
  Coord facing_y = 0;
  Coord diagonal = 0;
  bool overlap = false;
  bool line_end = false;  // line-end rules apply between the two classes
  bool sides = false;     // side or corner rules apply (sideViolation)
};

// Interned shape classes (thread-safe interning, lock-free reads).
class ShapeClassTable
{
 public:
  ShapeClassTable() = default;
  ShapeClassTable(const ShapeClassTable&) = delete;
  ShapeClassTable& operator=(const ShapeClassTable&) = delete;

  void setLayers(std::vector<LayerRules> layers);
  int numLayers() const { return static_cast<int>(layers_.size()); }
  const LayerRules& layer(const int l) const { return layers_[l]; }
  ClassId intern(const ShapeClass& c);
  // The class of a shape (Sec. 3.1 modelling); classify() also interns it.
  // `ndr`: the nondefault rule of a routed shape of an NDR net (metal
  // only; 0 = default rule).
  ShapeClass describe(int layer,
                      const Rect& rect,
                      ShapeRole role,
                      int ndr = 0) const;
  ClassId classify(int layer, const Rect& rect, ShapeRole role, int ndr = 0)
  {
    return intern(describe(layer, rect, role, ndr));
  }
  const ShapeClass& get(const ClassId id) const { return classes_[id]; }
  size_t size() const { return classes_.size(); }
  Requirement requirement(ClassId a, ClassId b) const;
  // Largest distance at which shapes of the two classes violate a rule
  // (spacing, line-end extensions; 1 for an overlap rule), 0 for classes
  // of different layers.
  Coord distance(ClassId a, ClassId b) const;

  // Interaction distances of a layer (Sec. 3.3/3.4 windows). ruleReach()
  // (also reach()) bounds every pair the layer's rules allow, including
  // classes interned later (partition halos, grid borders, callers that
  // classify shapes after choosing a window). classReach() bounds every
  // pair of classes interned so far and narrowReach() every pair of
  // interned narrow classes; both only grow, when intern() adds a class. A
  // window around a model of class c must extend by reachOf(c); the shape
  // grid reports wide pieces from the wider area on its own
  // (ShapeGrid::collect). narrowRuleReach() bounds every pair of narrow
  // classes the layer's rules allow, including classes interned later:
  // narrowReach() <= narrowRuleReach() <= ruleReach() at all times.
  Coord ruleReach(const int layer) const { return rule_reach_[layer]; }
  Coord narrowRuleReach(const int layer) const
  {
    return narrow_rule_reach_[layer];
  }
  Coord reach(const int layer) const { return rule_reach_[layer]; }
  Coord classReach(const int layer) const
  {
    return class_reach_[layer].load(std::memory_order_relaxed);
  }
  Coord narrowReach(const int layer) const
  {
    return narrow_reach_[layer].load(std::memory_order_relaxed);
  }
  // Wide classes: shapes whose width row, or whose nondefault rule's
  // spacing, needs at least wideRatioPercent() of the narrowest row's
  // largest spacing (macro blockages, power straps, NDR wires with wide
  // spacing). The shape grid keeps their shapes in a separate index, so
  // that they do not widen every window of their layer.
  bool wide(const ClassId id) const
  {
    const ShapeClass& c = get(id);
    return c.kind != ShapeClassKind::kCutProjection
           && (c.width_row >= wide_row_[c.layer]
               || (c.ndr > 0 && c.ndr <= wide_ndr_[c.layer].size()
                   && wide_ndr_[c.layer][c.ndr - 1]));
  }
  Coord reachOf(const ClassId id) const
  {
    const int layer = get(id).layer;
    return wide(id) ? classReach(layer) : narrowReach(layer);
  }
  // Threshold of wide(); takes effect at the next setLayers(). 0: no class
  // is wide.
  void setWideRatioPercent(const int percent) { wide_percent_ = percent; }
  int wideRatioPercent() const { return wide_percent_; }
  // Unique per setLayers() call (keys of per-thread caches).
  uint64_t generation() const { return generation_; }

 private:
  std::vector<LayerRules> layers_;
  std::vector<Coord> rule_reach_;
  std::vector<Coord> narrow_rule_reach_;
  std::vector<int> wide_row_;
  std::vector<std::vector<bool>> wide_ndr_;  // per layer, per NDR
  std::unique_ptr<std::atomic<Coord>[]> class_reach_;
  std::unique_ptr<std::atomic<Coord>[]> narrow_reach_;
  // Interned classes per layer (guarded by mutex_).
  std::vector<std::vector<ClassId>> layer_classes_;
  int wide_percent_ = 200;
  uint64_t generation_ = 0;
  std::mutex mutex_;
  std::unordered_map<ShapeClass, ClassId, ShapeClassHash> index_;
  AppendOnlyTable<ShapeClass> classes_;
};

// A shape at an absolute position with its class. line_end_ends marks which
// ends along the preferred direction are real line-ends (bit 0: low end,
// bit 1: high end); only meaningful for line-end classes. A shape stored in
// pieces (shape grid cells) marks only the pieces at the shape's ends.
// Bits 2-5 (cutSideBit) mark sides of the rectangle that are no sides of
// the shape (a cell or clip border of a stored piece); they are kept on
// layers with side rules only, where they tell which sides and corners
// are real.
struct PlacedShape
{
  Rect rect;
  ClassId cls = 0;
  uint8_t line_end_ends = 3;
};

// Sides in the order of ShapeClass::eol_sides.
constexpr int kSideXlo = 0;
constexpr int kSideXhi = 1;
constexpr int kSideYlo = 2;
constexpr int kSideYhi = 3;
constexpr uint8_t kCutSides = 0x3c;
inline uint8_t cutSideBit(const int side)
{
  return static_cast<uint8_t>(4u << side);
}

// Side and corner rules of two shapes (Lef58Rules.md) as rectangle tests:
// a and b violate them iff some pair (from_a, from_b) overlaps strictly
// (Rect::overlaps; degenerate rectangles are points or segments). Each
// rectangle depends on its own shape only, so the pairs also give the
// placements of a moving shape that violate (BrDistanceRules). `rule` is
// the ViolationBits bit of the rule.
struct SideTest
{
  Rect from_a;
  Rect from_b;
  uint32_t rule = 0;
};
void appendSideTests(const ShapeClassTable& classes,
                     const PlacedShape& a,
                     const PlacedShape& b,
                     std::vector<SideTest>& out);
// Any side or corner rule violated.
bool sideViolation(const ShapeClassTable& classes,
                   const PlacedShape& a,
                   const PlacedShape& b);

// Line-end extension rectangles of a placed shape (Fig. 2, hatched): beyond
// each marked end in preferred direction, `space` deep and `within` wider
// than the shape on both sides. Returns the number written to out.
int lineEndExtensions(const ShapeClassTable& classes,
                      const PlacedShape& shape,
                      Rect out[2]);

// Reference pairwise check (Sec. 3.1 / 3.2): do two shapes of different nets
// on the same shape-grid layer violate a minimum distance rule?
bool shapesViolate(const ShapeClassTable& classes,
                   const PlacedShape& a,
                   const PlacedShape& b);

// Which rule families two shapes violate (tests and reports): spacing /
// overlap (requirement), line-end extensions, the side rule kinds and the
// corner spacing. shapesViolate() == (mask != 0).
enum ViolationBits : uint32_t
{
  kViolSpacing = 1u << 0,
  kViolLineEnd = 1u << 1,
  kViolRegionOverlap = 1u << 2,
  kViolRegionCorner = 1u << 3,
  kViolFacing = 1u << 4,
  kViolFacingEol = 1u << 5,
  kViolCorner = 1u << 6,
};
uint32_t violationMask(const ShapeClassTable& classes,
                       const PlacedShape& a,
                       const PlacedShape& b);

// Sec. 3.2: "Wire and via models map one-dimensional abstractions of
// wires, called stick figures, to the metal shapes induced by them [...]
// encoded by an axis-parallel rectangle and a shape class." The metal of a
// wire is the Minkowski sum of its stick figure and `rect`.
struct WireModel
{
  bool valid = false;
  Rect rect;
  ClassId cls = 0;
};

// "A via model consists of three axis-parallel rectangles and shape
// classes. If an inter-layer via rule applies, the projection of the via
// cut to the next higher via layer is added to the via model."
struct ViaModel
{
  bool valid = false;
  WireModel bottom;      // on wiring layer z
  WireModel cut;         // on via layer z
  WireModel top;         // on wiring layer z + 1
  WireModel projection;  // on via layer z + 1, if an inter-layer rule applies
};

// "A wire type is a function that maps wiring layers to pairs of wire
// models (for preferred and non-preferred routing direction), and via layers
// to via models."
struct WireType
{
  std::string name;
  std::vector<WireModel> preferred;      // per wiring layer z
  std::vector<WireModel> non_preferred;  // per wiring layer z
  std::vector<ViaModel> via;             // per via layer z (z -> z + 1)
};

// Geometry of a wire type before classification.
struct WireTypeGeometry
{
  std::string name;
  // Nondefault rule index of the wire type's metal models
  // (LayerRules::ndr_spacing); 0 for the default rule.
  int ndr = 0;
  std::vector<Coord> width;  // per wiring layer; square stick-figure model
  struct Via
  {
    bool valid = false;
    Rect bottom;
    Rect cut;
    Rect top;
  };
  std::vector<Via> via;  // per via layer
};

// Builds a wire type: wiring models are squares of the layer's width
// centred on the stick figure (a wire extends half its width beyond its
// end points); via projections are added where the next higher via layer
// has an inter-layer rule.
WireType makeWireType(const WireTypeGeometry& geometry,
                      ShapeClassTable& classes);

}  // namespace drt::br
