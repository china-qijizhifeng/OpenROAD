// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// The distance rules of the BonnRoute routing space (BrWireModel.h) from the
// drt tech: one place for BonnRouteDR and the differential test against
// drt's GC (BonnRouteLef58DiffTest).

#pragma once

#include "br/BrWireModel.h"
#include "db/tech/frLayer.h"
#include "frBaseTypes.h"

namespace drt::br {

// Sec. 3.1: distances are nondecreasing functions of the widths and the
// common run-length. drt's spacing lookup sampled at every threshold t and
// t + 1 of its PRL table, two-widths table and range rules, so buckets
// agree with the lookup for strict and non-strict thresholds.
SpacingTable drtSpacingTable(frLayer* layer);

// Rules of a routing layer: spacing table and the line-end rule from drt's
// summary of the EOL rules (frLayer::getDrEolSpacingConstraint, filled by
// FlexRP).
LayerRules drtWiringRules(frLayer* layer);

// Rules of a cut layer: `cut_spacing` between its cuts, the inter-layer
// via rule with the cut layer `below` (nullptr: none), `horizontal` = the
// direction of the wiring layer below.
LayerRules drtViaRules(frLayer* cut,
                       frLayer* below,
                       frCoord cut_spacing,
                       bool horizontal);

}  // namespace drt::br
