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
#include <cstdint>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "br/BrGeometry.h"

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
  ShapeClass describe(int layer, const Rect& rect, ShapeRole role) const;
  ClassId classify(int layer, const Rect& rect, ShapeRole role)
  {
    return intern(describe(layer, rect, role));
  }
  const ShapeClass& get(const ClassId id) const { return classes_[id]; }
  size_t size() const { return classes_.size(); }
  Requirement requirement(ClassId a, ClassId b) const;
  Coord reach(const int layer) const { return reach_[layer]; }
  // Unique per setLayers() call (keys of per-thread caches).
  uint64_t generation() const { return generation_; }

 private:
  std::vector<LayerRules> layers_;
  std::vector<Coord> reach_;
  uint64_t generation_ = 0;
  std::mutex mutex_;
  std::unordered_map<ShapeClass, ClassId, ShapeClassHash> index_;
  AppendOnlyTable<ShapeClass> classes_;
};

// A shape at an absolute position with its class. line_end_ends marks which
// ends along the preferred direction are real line-ends (bit 0: low end,
// bit 1: high end); only meaningful for line-end classes. A shape stored in
// pieces (shape grid cells) marks only the pieces at the shape's ends.
struct PlacedShape
{
  Rect rect;
  ClassId cls = 0;
  uint8_t line_end_ends = 3;
};

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
