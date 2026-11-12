// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrDrtRules.h"

#include <algorithm>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "br/BrWireModel.h"
#include "db/tech/frConstraint.h"
#include "db/tech/frLayer.h"
#include "frBaseTypes.h"
#include "odb/db.h"

namespace drt::br {

SpacingTable drtSpacingTable(frLayer* layer)
{
  std::vector<Coord> widths;
  std::vector<Coord> prls;
  auto add = [](std::vector<Coord>& v, const frCoord t) {
    if (t >= 0) {
      v.push_back(t);
      v.push_back(t + 1);
    }
  };
  frConstraint* con = layer->getMinSpacing();
  if (con != nullptr
      && con->typeId() == frConstraintTypeEnum::frcSpacingTablePrlConstraint) {
    auto table = static_cast<frSpacingTablePrlConstraint*>(con)->getLookupTbl();
    for (const frCoord w : table.getRows()) {
      add(widths, w);
    }
    for (const frCoord p : table.getCols()) {
      add(prls, p);
    }
  }
  odb::dbTechLayer* db_layer = layer->getDbLayer();
  if (db_layer != nullptr && db_layer->hasTwoWidthsSpacingRules()) {
    const uint32_t rows = db_layer->getTwoWidthsSpacingTableNumWidths();
    for (uint32_t row = 0; row < rows; row++) {
      add(widths,
          static_cast<frCoord>(db_layer->getTwoWidthsSpacingTableWidth(row)));
      if (db_layer->getTwoWidthsSpacingTableHasPRL(row)) {
        add(prls,
            static_cast<frCoord>(db_layer->getTwoWidthsSpacingTablePRL(row)));
      }
    }
  }
  if (layer->hasSpacingRangeConstraints()) {
    for (const auto* range : layer->getSpacingRangeConstraints()) {
      add(widths, range->getMinWidth());
      add(widths, range->getMaxWidth());
    }
  }
  return SpacingTable::sample(
      widths, prls, [layer](Coord a, Coord b, Coord prl) {
        return layer->getMinSpacingValue(a, b, prl, false);
      });
}

namespace {

// Index of an EOL width in rules.eol_widths. Beyond kMaxEolWidths widths
// the largest one is raised if needed: a larger width only makes more
// sides end-of-line sides (conservative).
uint8_t eolIndex(LayerRules& rules, const Coord width)
{
  std::vector<Coord>& w = rules.eol_widths;
  for (size_t k = 0; k < w.size(); k++) {
    if (w[k] == width) {
      return static_cast<uint8_t>(k);
    }
  }
  if (w.size() < kMaxEolWidths) {
    w.push_back(width);
    return static_cast<uint8_t>(w.size() - 1);
  }
  const auto it = std::ranges::max_element(w);
  *it = std::max(*it, width);
  return static_cast<uint8_t>(it - w.begin());
}

}  // namespace

LayerRules drtWiringRules(frLayer* layer, const bool lef58)
{
  LayerRules rules;
  rules.horizontal = layer->isHorizontal();
  rules.spacing = drtSpacingTable(layer);
  const drEolSpacingConstraint& eol = layer->getDrEolSpacingConstraint();
  rules.line_end = {eol.eolWidth, eol.eolSpace, eol.eolWithin};
  const auto& lef58_eols = layer->getLef58SpacingEndOfLineConstraints();
  const auto& keepouts = layer->getLef58EolKeepOutConstraints();
  const auto& corners = layer->getLef58CornerSpacingConstraints();
  if (!lef58 || (lef58_eols.empty() && keepouts.empty() && corners.empty())) {
    return rules;
  }
  // Lef58Rules.md: the end-of-line rules become side rules. The summary
  // line-end rule stays only for LEF58 EOL extension rules, which have no
  // side rule.
  if (layer->getLef58EolExtConstraints().empty()) {
    rules.line_end = {};
  }
  // ENDOFLINE spacing: a facing side closer than the space, within
  // `within` laterally; PARALLELEDGE, MINLENGTH, ENCLOSECUT and the
  // "nothing in between" exemption are ignored (conservative); the larger
  // of the regular, wrong-direction and end-PRL spacings.
  for (const frSpacingEndOfLineConstraint* con : layer->getEolSpacing()) {
    rules.side_rules.push_back({SideRuleKind::kFacing,
                                eolIndex(rules, con->getEolWidth()),
                                con->getMinSpacing(),
                                0,
                                con->getEolWithin()});
  }
  for (const frLef58SpacingEndOfLineConstraint* con : lef58_eols) {
    const uint8_t k = eolIndex(rules, con->getEolWidth());
    Coord space = con->getEolSpace();
    if (con->hasWrongDirSpacing()) {
      space = std::max(space, con->getWrongDirSpace());
    }
    Coord within = 0;
    if (con->hasWithinConstraint()) {
      const auto w = con->getWithinConstraint();
      within = std::max(w->getEolWithin(), w->getWrongDirWithin());
      if (w->hasEndPrlSpacing()) {
        space = std::max(space, w->getEndPrlSpacing());
      }
      if (w->hasEndToEndConstraint()) {
        // ENDTOEND: both sides end-of-line sides, overlapping laterally
        // after bloating each by the extension.
        const auto ete = w->getEndToEndConstraint();
        const Coord ext = ete->hasExtension() ? ete->getExtension() : 0;
        rules.side_rules.push_back({SideRuleKind::kFacingEol,
                                    k,
                                    ete->getEndToEndSpace(),
                                    0,
                                    2 * std::max<Coord>(ext, 0)});
      }
    }
    rules.side_rules.push_back({SideRuleKind::kFacing, k, space, 0, within});
  }
  // EOLKEEPOUT; EXCEPTWITHIN and CLASS are ignored (conservative).
  for (const frLef58EolKeepOutConstraint* con : keepouts) {
    rules.side_rules.push_back({con->isCornerOnly()
                                    ? SideRuleKind::kRegionCorner
                                    : SideRuleKind::kRegionOverlap,
                                eolIndex(rules, con->getEolWidth()),
                                con->getForwardExt(),
                                con->getBackwardExt(),
                                con->getSideExt()});
  }
  // CORNERSPACING CONVEXCORNER: the largest spacing of the width table, as
  // the larger of x and y (Chebyshev); EXCEPTEOL and CORNERTOCORNER
  // (Euclidean) are modelled conservatively. CONCAVECORNER is not
  // modelled.
  for (const frLef58CornerSpacingConstraint* con : corners) {
    if (con->getCornerType() == frCornerTypeEnum::CONVEX) {
      rules.corner_spacing = std::max(rules.corner_spacing, con->findMax());
    }
  }
  return rules;
}

LayerRules drtViaRules(frLayer* cut,
                       frLayer* below,
                       const frCoord cut_spacing,
                       const bool horizontal,
                       const bool lef58)
{
  LayerRules rules;
  rules.via_layer = true;
  rules.horizontal = horizontal;
  rules.spacing = SpacingTable::constant(cut_spacing);
  if (below != nullptr) {
    // Inter-layer via rule with the cut layer below (Sec. 3.1, 3.2).
    auto inter = [](frLayer* a, frLayer* b) -> frCoord {
      const auto& cons = a->getInterLayerCutSpacingConstraint(false);
      const frLayerNum n = b->getLayerNum();
      if (n < 0 || n >= static_cast<int>(cons.size()) || cons[n] == nullptr) {
        return 0;
      }
      return cons[n]->getCutSpacing();
    };
    rules.projection_spacing = std::max(inter(cut, below), inter(below, cut));
  }
  if (!lef58 || !cut->hasLef58DiffNetCutSpcTblConstraint()) {
    return rules;
  }
  // LEF58_SPACINGTABLE between cut classes (FlexGC_cut.cpp): the largest
  // entry over SIDE / END and both PRL entries, both orders; DEFAULT for
  // cuts of no class. Centre-to-centre entries are taken edge to edge
  // (conservative). It applies to cuts of one net as well unless the layer
  // has a SAMENET / SAMEMETAL table.
  odb::dbTechLayerCutSpacingTableDefRule* db
      = cut->getLef58DiffNetCutSpcTblConstraint()->getODBRule();
  std::vector<std::vector<std::pair<int, int>>> table;
  std::map<std::string, uint32_t> row_map;
  std::map<std::string, uint32_t> col_map;
  db->getSpacingTable(table, row_map, col_map);
  std::vector<std::string> names{""};
  for (const auto* map : {&row_map, &col_map}) {
    for (const auto& [key, index] : *map) {
      const std::string name = key.substr(0, key.find_last_of('/'));
      const int idx = cut->getCutClassIdx(name);
      if (idx < 0 || std::ranges::find(names, name) != names.end()) {
        continue;
      }
      const frLef58CutClass* cc = cut->getCutClass(idx);
      names.push_back(name);
      rules.cut_classes.push_back(
          {std::min(cc->getViaWidth(), cc->getViaLength()),
           std::max(cc->getViaWidth(), cc->getViaLength())});
    }
  }
  const size_t n = names.size();
  rules.cut_class_spacing.assign(n * n, 0);
  using Strategy = odb::dbTechLayerCutSpacingTableDefRule::LOOKUP_STRATEGY;
  for (size_t a = 0; a < n; a++) {
    for (size_t b = 0; b < n; b++) {
      rules.cut_class_spacing[a * n + b]
          = std::max(db->getMaxSpacing(names[a], names[b], Strategy::MAX),
                     db->getMaxSpacing(names[b], names[a], Strategy::MAX));
    }
  }
  return rules;
}

bool padsCoverWireEnd(const Coord width, const std::vector<Rect>& pads)
{
  const Coord half = width / 2;
  const Rect square{-half, -half, width - half, width - half};
  return std::ranges::all_of(pads, [&](const Rect& p) {
    return p.xlo <= square.xlo && p.ylo <= square.ylo && p.xhi >= square.xhi
           && p.yhi >= square.yhi;
  });
}

}  // namespace drt::br
