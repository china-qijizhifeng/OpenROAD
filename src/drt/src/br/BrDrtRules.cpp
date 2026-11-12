// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrDrtRules.h"

#include <algorithm>
#include <cstdint>
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

LayerRules drtWiringRules(frLayer* layer)
{
  LayerRules rules;
  rules.horizontal = layer->isHorizontal();
  rules.spacing = drtSpacingTable(layer);
  const drEolSpacingConstraint& eol = layer->getDrEolSpacingConstraint();
  rules.line_end = {eol.eolWidth, eol.eolSpace, eol.eolWithin};
  return rules;
}

LayerRules drtViaRules(frLayer* cut,
                       frLayer* below,
                       const frCoord cut_spacing,
                       const bool horizontal)
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
  return rules;
}

}  // namespace drt::br
