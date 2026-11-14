// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "db/tech/frLayer.h"

#include <algorithm>

#include "frBaseTypes.h"
#include "utl/Logger.h"

namespace drt {

void frLayer::printAllConstraints(utl::Logger* logger)
{
  for (auto& constraint : constraints_) {
    constraint->report(logger);
  }
  if (lef58CutSpacingTableSameNetMetalConstraint_) {
    lef58CutSpacingTableSameNetMetalConstraint_->report(logger);
  }
  if (lef58CutSpacingTableDiffNetConstraint_) {
    lef58CutSpacingTableDiffNetConstraint_->report(logger);
  }
  if (lef58SameNetInterCutSpacingTableConstraint_) {
    lef58SameNetInterCutSpacingTableConstraint_->report(logger);
  }
  if (lef58SameMetalInterCutSpacingTableConstraint_) {
    lef58SameMetalInterCutSpacingTableConstraint_->report(logger);
  }
  if (lef58DefaultInterCutSpacingTableConstraint_) {
    lef58DefaultInterCutSpacingTableConstraint_->report(logger);
  }
  for (auto& constraint : lef58SpacingEndOfLineConstraints_) {
    constraint->report(logger);
  }
  if (minSpc_) {
    minSpc_->report(logger);
  }
  if (spacingSamenet_) {
    spacingSamenet_->report(logger);
  }
  for (auto& constraint : eols_) {
    if (constraint) {
      constraint->report(logger);
    }
  }
  for (auto& constraint : cutConstraints_) {
    if (constraint) {
      constraint->report(logger);
    }
  }
  for (auto& constraint : cutSpacingSamenetConstraints_) {
    if (constraint) {
      constraint->report(logger);
    }
  }
  for (auto& constraint : interLayerCutSpacingConstraints_) {
    if (constraint) {
      constraint->report(logger);
    }
  }
  for (auto& constraint : interLayerCutSpacingSamenetConstraints_) {
    if (constraint) {
      constraint->report(logger);
    }
  }
  for (auto& constraint : lef58CutSpacingConstraints_) {
    if (constraint) {
      constraint->report(logger);
    }
  }
  for (auto& constraint : lef58CutSpacingSamenetConstraints_) {
    if (constraint) {
      constraint->report(logger);
    }
  }
  for (auto& constraint : minEnclosedAreaConstraints_) {
    if (constraint) {
      constraint->report(logger);
    }
  }
  if (areaConstraint_) {
    areaConstraint_->report(logger);
  }
  if (minStepConstraint_) {
    minStepConstraint_->report(logger);
  }
  for (auto& constraint : lef58MinStepConstraints_) {
    if (constraint) {
      constraint->report(logger);
    }
  }
  if (minWidthConstraint_) {
    minWidthConstraint_->report(logger);
  }
  for (auto& constraint : minimumcutConstraints_) {
    if (constraint) {
      constraint->report(logger);
    }
  }
  if (lef58RectOnlyConstraint_) {
    lef58RectOnlyConstraint_->report(logger);
  }
  if (lef58RightWayOnGridOnlyConstraint_) {
    lef58RightWayOnGridOnlyConstraint_->report(logger);
  }
  for (auto& constraint : lef58CornerSpacingConstraints_) {
    if (constraint) {
      constraint->report(logger);
    }
  }
}

void frLayer::setMinSpacing(frConstraint* in)
{
  minSpc_ = in;
  if (in == nullptr) {
    min_spacing_type_ = MinSpacingType::NONE;
    return;
  }
  switch (in->typeId()) {
    case frConstraintTypeEnum::frcSpacingConstraint:
      min_spacing_type_ = MinSpacingType::SPACING;
      break;
    case frConstraintTypeEnum::frcSpacingTablePrlConstraint:
      min_spacing_type_ = MinSpacingType::PRL_TABLE;
      break;
    case frConstraintTypeEnum::frcSpacingTableTwConstraint:
      min_spacing_type_ = MinSpacingType::TWO_WIDTHS_TABLE;
      break;
    default:
      min_spacing_type_ = MinSpacingType::OTHER;
      break;
  }
}

frCoord frLayer::getMinSpacingValue(frCoord width1,
                                    frCoord width2,
                                    frCoord prl,
                                    bool use_min_spacing)
{
  frCoord rangeSpc = -1;
  for (const auto* con : spacingRangeConstraints_) {
    if (use_min_spacing) {
      if (rangeSpc == -1) {
        rangeSpc = con->getMinSpacing();
      } else {
        rangeSpc = std::min(rangeSpc, con->getMinSpacing());
      }
    } else if (con->inRange(width1) || con->inRange(width2)) {
      rangeSpc = std::max(rangeSpc, con->getMinSpacing());
    }
  }

  frCoord minSpc = 0;
  switch (min_spacing_type_) {
    case MinSpacingType::NONE:
      return std::max(rangeSpc, 0);
    case MinSpacingType::SPACING:
      minSpc = static_cast<frSpacingConstraint*>(minSpc_)->getMinSpacing();
      break;
    case MinSpacingType::PRL_TABLE: {
      const auto* con = static_cast<frSpacingTablePrlConstraint*>(minSpc_);
      minSpc = use_min_spacing ? con->findMin()
                               : con->find(std::max(width1, width2), prl);
      break;
    }
    case MinSpacingType::TWO_WIDTHS_TABLE: {
      const auto* con = static_cast<frSpacingTableTwConstraint*>(minSpc_);
      minSpc
          = use_min_spacing ? con->findMin() : con->find(width1, width2, prl);
      break;
    }
    case MinSpacingType::OTHER:
      break;
  }
  if (use_min_spacing && rangeSpc != -1) {
    return std::min(rangeSpc, minSpc);
  }
  return std::max(rangeSpc, minSpc);
}

}  // namespace drt
