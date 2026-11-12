// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Random rules, layouts and the brute-force distance-rule reference shared
// by the BonnRoute routing-space tests: every stored shape is kept in a
// vector and checked pairwise with the same shape classes.

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "br/BrDistanceRules.h"
#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrWireModel.h"

namespace drt::br::test {

inline int randomInt(std::mt19937& rng, int lo, int hi)
{
  return lo + static_cast<int>(rng() % static_cast<unsigned>(hi - lo + 1));
}

// Side rules (Lef58Rules.md): EOL widths, keep-out regions (whole and
// corner form), facing and end-to-end rules, convex corner spacing.
inline void randomSideRules(std::mt19937& rng, LayerRules& r, const Coord base)
{
  for (int k = randomInt(rng, 1, 3); k > 0; k--) {
    r.eol_widths.push_back(randomInt(rng, 6, 36));
  }
  const int n = static_cast<int>(r.eol_widths.size());
  for (int k = randomInt(rng, 0, 4); k > 0; k--) {
    SideRule rule;
    rule.kind = static_cast<SideRuleKind>(rng() % 4);
    rule.eol = static_cast<uint8_t>(randomInt(rng, 0, n - 1));
    rule.forward = randomInt(rng, 1, 2 * base);
    if (rule.kind == SideRuleKind::kRegionOverlap
        || rule.kind == SideRuleKind::kRegionCorner) {
      rule.backward = rng() % 3 == 0 ? 0 : randomInt(rng, 1, 2 * base);
    }
    rule.side = rng() % 4 == 0 ? 0 : randomInt(rng, 1, base);
    r.side_rules.push_back(rule);
  }
  if (rng() % 2 == 0) {
    r.corner_spacing = randomInt(rng, 2, 2 * base);
  }
  r.wire_ends_covered = rng() % 2 == 0;
}

// Spacing rules with run-length and width dependent steps, thresholds of
// both strict and non-strict kind, line-end rules, inter-layer via rules,
// side rules and cut class spacing tables (Lef58Rules.md).
inline std::vector<LayerRules> randomRules(std::mt19937& rng,
                                           const int num_wiring_layers,
                                           const bool first_horizontal)
{
  std::vector<LayerRules> rules;
  for (int l = 0; l < 2 * num_wiring_layers - 1; l++) {
    LayerRules r;
    const int z = l / 2;
    r.via_layer = (l % 2) == 1;
    r.horizontal = ((z % 2) == 0) == first_horizontal;
    if (r.via_layer) {
      r.spacing = SpacingTable::constant(randomInt(rng, 6, 16));
      if (l >= 3 && rng() % 2 == 0) {
        r.projection_spacing = randomInt(rng, 4, 14);
      }
      if (rng() % 2 == 0) {
        // Cut classes among the cut sizes of randomShape (4..12) and a
        // symmetric spacing table; index 0 = no class (the DEFAULT).
        for (int k = randomInt(rng, 1, 3); k > 0; k--) {
          const Coord a = randomInt(rng, 4, 12);
          const Coord b = randomInt(rng, 4, 12);
          r.cut_classes.push_back({std::min(a, b), std::max(a, b)});
        }
        const size_t n = r.cut_classes.size() + 1;
        r.cut_class_spacing.assign(n * n, 0);
        for (size_t a = 0; a < n; a++) {
          for (size_t b = a; b < n; b++) {
            const Coord v = randomInt(rng, 0, 24);
            r.cut_class_spacing[a * n + b] = v;
            r.cut_class_spacing[b * n + a] = v;
          }
        }
      }
    } else {
      const Coord base = randomInt(rng, 6, 14);
      const Coord wide = randomInt(rng, 14, 30);
      const Coord long_run = randomInt(rng, 20, 60);
      const Coord inc_w = randomInt(rng, 0, 8);
      const Coord inc_p = randomInt(rng, 0, 6);
      const bool strict = rng() % 2;
      r.spacing = SpacingTable::sample(
          {wide, wide + 1},
          {long_run, long_run + 1},
          [&](Coord wa, Coord wb, Coord prl) {
            const Coord w = std::max(wa, wb);
            const bool is_wide = strict ? w > wide : w >= wide;
            const bool is_long = strict ? prl > long_run : prl >= long_run;
            return base + (is_wide ? inc_w : 0)
                   + (is_wide && is_long ? inc_p : 0);
          });
      if (rng() % 3 != 0) {
        r.line_end.max_edge = randomInt(rng, 10, 24);
        r.line_end.space = randomInt(rng, base, base + 10);
        r.line_end.within = randomInt(rng, 1, base);
      }
      // Nondefault rules: spacing of NDR wires (0: none on this layer).
      for (int k = rng() % 3; k > 0; k--) {
        r.ndr_spacing.push_back(
            rng() % 4 == 0 ? 0 : randomInt(rng, base, 3 * base));
      }
      if (rng() % 2 == 0) {
        randomSideRules(rng, r, base);
      }
    }
    rules.push_back(r);
  }
  return rules;
}

struct StoredShape
{
  int layer;
  Rect rect;
  ClassId cls;
  ShapeType type;
  int32_t owner;
};

inline ShapeType randomType(std::mt19937& rng)
{
  return static_cast<ShapeType>(rng() % 5);
}

// Threshold of ShapeClassTable::wide() for a random case: none, every row
// above the first, or rows with clearly larger spacing.
inline int randomWideRatio(std::mt19937& rng)
{
  constexpr int kRatios[] = {0, 100, 130, 200};
  return kRatios[rng() % 4];
}

// A random shape on `layer` inside [0, extent)^2.
inline StoredShape randomShape(std::mt19937& rng,
                               ShapeClassTable& classes,
                               const int layer,
                               const Coord extent,
                               const int nets)
{
  const LayerRules& rules = classes.layer(layer);
  StoredShape s;
  s.layer = layer;
  s.type = randomType(rng);
  s.owner = isOwned(s.type) ? randomInt(rng, 0, nets - 1) : kNoNet;
  ShapeRole role = ShapeRole::kFixed;
  Coord dx;
  Coord dy;
  if (rules.via_layer) {
    dx = randomInt(rng, 4, 12);
    dy = rng() % 2 ? dx : randomInt(rng, 4, 12);
    role = rng() % 4 == 0 && layer >= 3 && rules.projection_spacing > 0
               ? ShapeRole::kCutProjection
               : ShapeRole::kCut;
  } else {
    const int kind = rng() % 4;
    const Coord w = randomInt(rng, 6, 34);
    const Coord len = randomInt(rng, 4, extent / 2);
    if (kind == 0) {
      role = ShapeRole::kPreferredWire;
      dx = rules.horizontal ? len : w;
      dy = rules.horizontal ? w : len;
    } else if (kind == 1) {
      role = ShapeRole::kNonPreferredWire;
      dx = rules.horizontal ? w : len / 3 + 4;
      dy = rules.horizontal ? len / 3 + 4 : w;
    } else if (kind == 2) {
      role = ShapeRole::kViaPad;
      dx = randomInt(rng, 8, 24);
      dy = randomInt(rng, 8, 24);
    } else {
      dx = randomInt(rng, 3, 80);
      dy = randomInt(rng, 3, 80);
    }
  }
  dx = std::min(dx, extent - 2);
  dy = std::min(dy, extent - 2);
  const Coord x = randomInt(rng, 1, extent - dx - 1);
  const Coord y = randomInt(rng, 1, extent - dy - 1);
  s.rect = {x, y, x + dx, y + dy};
  // Routed metal of an NDR net now and then.
  const int ndr
      = !rules.via_layer && !rules.ndr_spacing.empty()
            ? randomInt(rng, 0, static_cast<int>(rules.ndr_spacing.size()))
            : 0;
  s.cls = classes.classify(layer, s.rect, role, ndr);
  return s;
}

// Brute-force level of placed models in groups at a location.
inline std::vector<RipupLevel> bruteLevels(
    const ShapeClassTable& classes,
    const ShapeGrid& grid,
    const std::vector<StoredShape>& shapes,
    const std::vector<LocatedModel>& models,
    const int groups,
    const Point l,
    const int32_t ignored_owner,
    std::vector<BlockingNet>* nets = nullptr)
{
  std::vector<RipupLevel> levels(groups, 0);
  struct Hit
  {
    int group;
    int32_t net;
    RipupLevel level;
    bool removable;
  };
  std::vector<Hit> hits;
  for (size_t k = 0; k < models.size(); k++) {
    const LocatedModel& m = models[k];
    const int g = m.group < 0 ? static_cast<int>(k) : m.group;
    const PlacedShape placed{m.rect.translated(l), m.cls, m.line_end_ends};
    for (const StoredShape& s : shapes) {
      if (s.layer != m.layer
          || (s.owner == ignored_owner && s.owner != kNoNet)) {
        continue;
      }
      if (shapesViolate(classes, placed, PlacedShape{s.rect, s.cls, 3})) {
        const RipupLevel lv = grid.level(s.type, s.owner);
        levels[g] = std::max(levels[g], lv);
        hits.push_back({g, s.owner, lv, isRemovable(s.type)});
      }
    }
  }
  if (nets != nullptr) {
    nets->clear();
    for (const Hit& h : hits) {
      if (!h.removable || levels[h.group] >= kNotRemovable) {
        continue;
      }
      auto at = std::ranges::find(*nets, h.net, &BlockingNet::net);
      if (at == nets->end()) {
        nets->push_back({h.net, 0, 0});
        at = nets->end() - 1;
      }
      at->level = std::max(at->level, h.level);
      at->groups |= uint64_t{1} << h.group;
    }
    std::ranges::sort(*nets, {}, &BlockingNet::net);
  }
  return levels;
}

inline bool sameAnswer(const std::vector<RipupLevel>& a,
                       const std::vector<RipupLevel>& b,
                       const bool levels)
{
  if (a.size() != b.size()) {
    return false;
  }
  for (size_t k = 0; k < a.size(); k++) {
    if (levels ? a[k] != b[k] : (a[k] == 0) != (b[k] == 0)) {
      return false;
    }
  }
  return true;
}

// Shape-grid geometry over [0, extent)^2 with the paper's cell size bound.
inline std::vector<ShapeGridLayerGeometry> gridGeometry(
    std::mt19937& rng,
    const ShapeClassTable& classes,
    const Coord extent,
    const bool random_cells)
{
  std::vector<ShapeGridLayerGeometry> layers;
  for (int l = 0; l < classes.numLayers(); l++) {
    ShapeGridLayerGeometry g;
    g.horizontal = classes.layer(l).horizontal;
    const Coord bound = std::max<Coord>(1, maximalCellSize(classes, l));
    g.along_cell = random_cells ? randomInt(rng, 1, bound) : bound;
    g.cross_cell = random_cells ? randomInt(rng, 1, bound) : bound;
    g.along_origin = random_cells ? -randomInt(rng, 0, g.along_cell - 1) : 0;
    g.cross_origin = random_cells ? -randomInt(rng, 0, g.cross_cell - 1) : 0;
    g.along_cells = (extent - g.along_origin) / g.along_cell + 1;
    g.cross_cells = (extent - g.cross_origin) / g.cross_cell + 1;
    layers.push_back(g);
  }
  return layers;
}

}  // namespace drt::br::test
