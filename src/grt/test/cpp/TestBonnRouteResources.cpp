// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Database-free tests of the Sec. 2.1 / Fig. 1 resource models
// (BonnRouteResources.h) and of the extra-space property format
// (grt/BonnRouteExtraSpace.h).

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "BonnRouteResources.h"
#include "grt/BonnRouteExtraSpace.h"
#include "gtest/gtest.h"

namespace grt::bonnroute {
namespace {

SpaceLayerModel nangateLikeModel()
{
  SpaceLayerModel m;
  m.gap = 70.0;  // pitch 140, width 70 (DBU 1000/um)
  m.growth = 70.0;
  m.ground_capacitance = 0.03;
  m.coupling_capacitance = 0.05;
  m.defect_size = 65.0;
  m.wire_width = 70.0;
  return m;
}

TEST(BonnRouteResources, GammaFunctionsDecreaseAndAreConvex)
{
  const SpaceLayerModel m = nangateLikeModel();
  const double h = 0.05;
  for (double s = 0.0; s <= 4.0; s += h) {
    for (auto gamma :
         {&SpaceLayerModel::powerPerLength, &SpaceLayerModel::yieldPerLength}) {
      const double f0 = (m.*gamma)(s);
      const double f1 = (m.*gamma)(s + h);
      const double f2 = (m.*gamma)(s + 2 * h);
      EXPECT_GT(f0, 0.0);
      EXPECT_LT(f1, f0);                    // decreasing in s
      EXPECT_GE(f0 - 2 * f1 + f2, -1e-15);  // convex
    }
  }
  // At s = 0 the coupling capacitance is fully present; far away only the
  // ground part (and for yield the opens part) is left.
  EXPECT_DOUBLE_EQ(m.powerPerLength(0.0),
                   m.ground_capacitance + m.coupling_capacitance);
  EXPECT_NEAR(m.powerPerLength(1e9), m.ground_capacitance, 1e-9);
  EXPECT_DOUBLE_EQ(m.yieldPerLength(0.0),
                   m.defect_size * m.defect_size / m.gap
                       + m.defect_size * m.defect_size / (2 * m.wire_width));
  // The shapes are the coefficients of 1 / (g + a s).
  for (const double s : {0.0, 0.5, 1.5}) {
    EXPECT_NEAR(m.powerPerLength(s),
                m.ground_capacitance + m.powerShape() / (m.gap + m.growth * s),
                1e-15);
    EXPECT_NEAR(m.yieldPerLength(s),
                m.defect_size * m.defect_size / (2 * m.wire_width)
                    + m.yieldShape() / (m.gap + m.growth * s),
                1e-12);
  }
}

TEST(BonnRouteResources, LayerModelUsesEitherCapacitanceSource)
{
  // Nangate45 metal2 at 2000 DBU/um: width 140, pitch 380. set_layer_rc
  // leaves an area capacitance (8.081e-5 pF/um / 0.07 um) and no edge
  // capacitance; the LEF has both. Either way the wire's total capacitance
  // per length is split into ground and coupling by coupling_fraction, so
  // the power model depends on s.
  const double rc_area = 8.081e-5 / 0.07;
  for (const auto& [area, edge] :
       {std::pair{rc_area, 0.0}, std::pair{4.0896e-5, 2.5157e-5}}) {
    const SpaceLayerModel m
        = makeSpaceLayerModel(140, 380, 140, area, edge, 2000, 0.5, 140.0);
    const double per_um = area * 0.07 + 2.0 * edge;
    EXPECT_NEAR((m.ground_capacitance + m.coupling_capacitance) * 2000.0,
                per_um * 1000.0,
                1e-12);
    EXPECT_GT(m.coupling_capacitance, 0.0);
    EXPECT_DOUBLE_EQ(m.coupling_capacitance, m.ground_capacitance);
    EXPECT_DOUBLE_EQ(m.gap, 240.0);
    EXPECT_DOUBLE_EQ(m.growth, 190.0);
    EXPECT_LT(m.powerPerLength(2.0), m.powerPerLength(0.0));
  }
  // No data: unit coupling capacitance; no gap from pitch: spacing.
  const SpaceLayerModel none
      = makeSpaceLayerModel(100, 100, 40, 0.0, 0.0, 1000, 0.5, 10.0);
  EXPECT_GT(none.coupling_capacitance, 0.0);
  EXPECT_DOUBLE_EQ(none.gap, 40.0);
}

TEST(BonnRouteResources, SpaceGridPoints)
{
  EXPECT_EQ((SpaceGrid{0.5, 2.0}).points(), 5);
  EXPECT_EQ((SpaceGrid{0.3, 1.0}).points(), 4);
  EXPECT_EQ((SpaceGrid{1.0, 0.0}).points(), 1);
  EXPECT_EQ((SpaceGrid{0.0, 2.0}).points(), 1);
}

TEST(BonnRouteResources, MinimizeSpaceMatchesBruteForce)
{
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  int nonzero = 0;
  int interior = 0;
  for (int i = 0; i < 20000; i++) {
    const SpaceGrid grid{0.25 * (1 + static_cast<int>(rng() % 4)),
                         0.5 * static_cast<int>(rng() % 9)};
    const double gap = 10.0 + 200.0 * unit(rng);
    const double growth = 10.0 + 200.0 * unit(rng);
    const double price = (rng() % 10 == 0) ? 0.0 : std::exp(8 * unit(rng) - 6);
    const double shape = (rng() % 10 == 0) ? 0.0 : std::exp(12 * unit(rng) - 4);
    const SpaceChoice fast = minimizeSpace(price, shape, gap, growth, grid);
    const SpaceChoice slow
        = minimizeSpaceBruteForce(price, shape, gap, growth, grid);
    ASSERT_NEAR(fast.cost, slow.cost, 1e-12 * std::max(1.0, slow.cost))
        << "price " << price << " shape " << shape << " gap " << gap
        << " growth " << growth << " step " << grid.step << " max " << grid.max;
    if (std::abs(fast.space - slow.space) > 1e-12) {
      // Only an exact tie may pick a different grid point.
      ASSERT_NEAR(fast.cost, slow.cost, 1e-14 * std::max(1.0, slow.cost));
    }
    nonzero += slow.space > 0.0;
    interior += slow.space > 0.0 && slow.space < grid.max;
  }
  // The sample covers s = 0, interior optima and the upper bound.
  EXPECT_GT(nonzero, 1000);
  EXPECT_GT(interior, 1000);
}

TEST(BonnRouteResources, MinimizeSpaceLimits)
{
  const SpaceGrid grid{0.5, 2.0};
  // No benefit from space: s = 0.
  EXPECT_EQ(minimizeSpace(1.0, 0.0, 50.0, 50.0, grid).space, 0.0);
  // Free space: the largest grid value.
  EXPECT_EQ(minimizeSpace(0.0, 1.0, 50.0, 50.0, grid).space, 2.0);
  // Expensive space: s = 0.
  EXPECT_EQ(minimizeSpace(1e6, 1.0, 50.0, 50.0, grid).space, 0.0);
}

TEST(BonnRouteExtraSpace, PropertyRoundTrip)
{
  const std::vector<ExtraSpaceEdge> edges{
      {2, 0, 1, 0.5f}, {3, 17, 4, 2.0f}, {10, 1234, 5678, 1.25f}};
  const std::string text = encodeExtraSpace(edges);
  EXPECT_EQ(text, "v1 2,0,1,0.500 3,17,4,2.000 10,1234,5678,1.250");
  std::vector<ExtraSpaceEdge> back;
  ASSERT_TRUE(decodeExtraSpace(text, back));
  EXPECT_EQ(back, edges);
  ASSERT_TRUE(decodeExtraSpace("v1", back));
  EXPECT_TRUE(back.empty());
  EXPECT_FALSE(decodeExtraSpace("v2 1,2,3,1.0", back));
  EXPECT_FALSE(decodeExtraSpace("v1 1,2", back));
  EXPECT_FALSE(decodeExtraSpace("v1 1,2,3,1.0x", back));
  EXPECT_TRUE(back.empty());
}

}  // namespace
}  // namespace grt::bonnroute
