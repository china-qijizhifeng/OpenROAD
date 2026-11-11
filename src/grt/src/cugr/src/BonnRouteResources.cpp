// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "BonnRouteResources.h"

#include <algorithm>
#include <cmath>

namespace grt::bonnroute {

double SpaceLayerModel::powerPerLength(const double s) const
{
  return ground_capacitance + coupling_capacitance * gap / (gap + growth * s);
}

double SpaceLayerModel::yieldPerLength(const double s) const
{
  return defect_size
         * (defect_size / (gap + growth * s)
            + defect_size / (2.0 * wire_width));
}

SpaceLayerModel makeSpaceLayerModel(const int width,
                                    const int pitch,
                                    const int spacing,
                                    const double area_capacitance,
                                    const double edge_capacitance,
                                    const int dbu_per_micron,
                                    const double coupling_fraction,
                                    const double defect_size)
{
  SpaceLayerModel m;
  const double w = std::max(1, width);
  const double p = std::max(1, pitch);
  const double dbu = std::max(1, dbu_per_micron);
  m.gap = p - w > 0.0 ? p - w : std::max(1, spacing);
  m.growth = p / 2.0;
  m.wire_width = w;
  m.defect_size = defect_size;
  // pF per um of wire, then fF per DBU.
  const double per_um = std::max(0.0, area_capacitance) * w / dbu
                        + 2.0 * std::max(0.0, edge_capacitance);
  const double to_ff_per_dbu = 1000.0 / dbu;
  const double total = per_um > 0.0 ? per_um * to_ff_per_dbu : to_ff_per_dbu;
  const double k = std::clamp(coupling_fraction, 0.0, 1.0);
  m.coupling_capacitance = k * total;
  m.ground_capacitance = (1.0 - k) * total;
  return m;
}

int SpaceGrid::points() const
{
  if (step <= 0.0 || max <= 0.0) {
    return 1;
  }
  return 1 + static_cast<int>(std::floor(max / step + 1e-9));
}

namespace {
double spaceCost(const double price,
                 const double shape,
                 const double gap,
                 const double growth,
                 const double s)
{
  return price * s + shape / (gap + growth * s);
}
}  // namespace

SpaceChoice minimizeSpace(const double space_price,
                          const double shape,
                          const double gap,
                          const double growth,
                          const SpaceGrid& grid)
{
  const int n = grid.points();
  if (n == 1 || shape <= 0.0) {
    return {0.0, spaceCost(space_price, shape, gap, growth, 0.0)};
  }
  // d/ds = a_price - shape growth / (gap + growth s)^2 = 0.
  double index;
  if (space_price <= 0.0) {
    index = n - 1;
  } else {
    const double s = (std::sqrt(shape * growth / space_price) - gap) / growth;
    index = std::clamp(s / grid.step, 0.0, static_cast<double>(n - 1));
  }
  const int lo = static_cast<int>(std::floor(index));
  const int hi = std::min(n - 1, lo + 1);
  SpaceChoice best{grid.value(lo),
                   spaceCost(space_price, shape, gap, growth, grid.value(lo))};
  const double cost_hi
      = spaceCost(space_price, shape, gap, growth, grid.value(hi));
  if (cost_hi < best.cost) {
    best = {grid.value(hi), cost_hi};
  }
  return best;
}

SpaceChoice minimizeSpaceBruteForce(const double space_price,
                                    const double shape,
                                    const double gap,
                                    const double growth,
                                    const SpaceGrid& grid)
{
  SpaceChoice best{0.0, spaceCost(space_price, shape, gap, growth, 0.0)};
  for (int i = 1; i < grid.points(); i++) {
    const double cost
        = spaceCost(space_price, shape, gap, growth, grid.value(i));
    if (cost < best.cost) {
      best = {grid.value(i), cost};
    }
  }
  return best;
}

}  // namespace grt::bonnroute
