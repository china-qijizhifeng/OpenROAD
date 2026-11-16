// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>

#include "ant/AntennaChecker.hh"

namespace grt {

inline int antennaDiodeDemand(const ant::Violations& violations)
{
  int demand = 0;
  for (const auto& violation : violations) {
    // A violation without a diode solution must not look cheaper than a
    // repaired violation (also covers -jumper_only without a diode model).
    demand
        += std::max(1, violation.diode_count_per_gate) * violation.gates.size();
  }
  return demand;
}

inline bool antennaRatioDoesNotWorsen(const ant::RatioCheck& before,
                                      const ant::RatioCheck& after)
{
  if (!after.violated) {
    return true;
  }
  if (!before.violated || !std::isfinite(after.actual)
      || !std::isfinite(after.required)) {
    return false;
  }
  const double old_ratio = before.actual / std::max(before.required, 1e-30);
  const double new_ratio = after.actual / std::max(after.required, 1e-30);
  return new_ratio <= old_ratio + 1e-9 * std::max(1.0, old_ratio);
}

// Demand must fall and no gate/layer/ratio family may gain a violation or
// worsen an existing one. In particular, a PAR improvement cannot hide a new
// cumulative or cut-layer violation, or a gate lost by the trial wire builder.
inline bool localAntennaImproves(const ant::NetCheck& before,
                                 const ant::NetCheck& after,
                                 int demand_before,
                                 int demand_after)
{
  if (demand_after >= demand_before
      || before.pins.size() != after.pins.size()) {
    return false;
  }
  std::map<odb::dbITerm*, const ant::PinCheck*, odb::ODBPtrLess> old_pins;
  for (const auto& pin : before.pins) {
    old_pins.emplace(pin.gate, &pin);
  }
  for (const auto& pin : after.pins) {
    const auto old_pin = old_pins.find(pin.gate);
    if (old_pin == old_pins.end()) {
      return false;
    }
    std::map<odb::dbTechLayer*, const ant::LayerCheck*, odb::ODBPtrLess> layers;
    for (const auto& layer : old_pin->second->layers) {
      layers.emplace(layer.layer, &layer);
    }
    for (const auto& layer : pin.layers) {
      if (!layer.violated) {
        continue;
      }
      const auto old_layer = layers.find(layer.layer);
      if (old_layer == layers.end()) {
        return false;
      }
      const auto& old = *old_layer->second;
      if (!antennaRatioDoesNotWorsen(old.par, layer.par)
          || !antennaRatioDoesNotWorsen(old.car, layer.car)
          || (layer.psr
              && (!old.psr || !antennaRatioDoesNotWorsen(*old.psr, *layer.psr)))
          || (layer.csr
              && (!old.csr
                  || !antennaRatioDoesNotWorsen(*old.csr, *layer.csr)))) {
        return false;
      }
    }
    old_pins.erase(old_pin);
  }
  return old_pins.empty();
}

}  // namespace grt
