// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include <algorithm>
#include <array>
#include <limits>
#include <random>

#include "dr/SpacingGeometry.h"
#include "frBaseTypes.h"
#include "gtest/gtest.h"
#include "odb/geom.h"

namespace drt {
namespace {

frSquaredDistance fourCornerDistance(const odb::Point& center,
                                     frCoord half_width,
                                     const odb::Rect& box)
{
  frSquaredDistance result = std::numeric_limits<frSquaredDistance>::max();
  for (const int sx : {-1, 1}) {
    for (const int sy : {-1, 1}) {
      const frCoord x = center.x() + sx * half_width;
      const frCoord y = center.y() + sy * half_width;
      const frCoord dx = std::max({box.xMin() - x, x - box.xMax(), 0});
      const frCoord dy = std::max({box.yMin() - y, y - box.yMax(), 0});
      result = std::min(result,
                        static_cast<frSquaredDistance>(dx) * dx
                            + static_cast<frSquaredDistance>(dy) * dy);
    }
  }
  return result;
}

TEST(SpacingGeometry, SpanningCornersRemainOutside)
{
  EXPECT_EQ(minCornerAxisDistanceSquared(0, 20, -5, 5), 225);
  EXPECT_EQ(minCornerAxisDistanceSquared(0, 5, -5, 5), 0);
  EXPECT_EQ(minCornerAxisDistanceSquared(-10, 0, -5, 5), 25);
}

TEST(SpacingGeometry, MatchesFourCornerOracle)
{
  std::mt19937 random(20261109);
  std::uniform_int_distribution<frCoord> coordinate(-1000000, 1000000);
  std::uniform_int_distribution<frCoord> extent(0, 10000);
  for (int i = 0; i < 100000; ++i) {
    const frCoord x = coordinate(random);
    const frCoord y = coordinate(random);
    const odb::Rect box(x, y, x + extent(random), y + extent(random));
    const odb::Point center(coordinate(random), coordinate(random));
    const frCoord half_width = extent(random);
    const auto distance = minCornerAxisDistanceSquared(
                              center.x(), half_width, box.xMin(), box.xMax())
                          + minCornerAxisDistanceSquared(
                              center.y(), half_width, box.yMin(), box.yMax());
    const auto expected = fourCornerDistance(center, half_width, box);
    ASSERT_EQ(distance, expected);
    for (const frCoord spacing : std::array<frCoord, 4>{0, 1, 100, 10000}) {
      const auto threshold = static_cast<frSquaredDistance>(spacing) * spacing;
      ASSERT_EQ(distance < threshold, expected < threshold);
    }
  }
}

TEST(SpacingGeometry, UsesWideSquaredDistances)
{
  EXPECT_EQ(minCornerAxisDistanceSquared(1000000000, 1, 0, 0),
            999999998000000001ULL);
}

}  // namespace
}  // namespace drt
