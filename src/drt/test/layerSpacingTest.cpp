// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include <algorithm>
#include <vector>

#include "db/tech/frConstraint.h"
#include "db/tech/frLayer.h"
#include "db/tech/frLookupTbl.h"
#include "gtest/gtest.h"
#include "odb/db.h"
#include "odb/dbTypes.h"

namespace drt {
namespace {

TEST(LayerSpacing, RuleReplacementAndMutableValues)
{
  frLayer layer;
  EXPECT_EQ(layer.getMinSpacingValue(100, 200, 300, false), 0);
  frSpacingConstraint spacing(45);
  layer.setMinSpacing(&spacing);
  EXPECT_EQ(layer.getMinSpacingValue(100, 200, 300, false), 45);
  spacing.setMinSpacing(70);
  EXPECT_EQ(layer.getMinSpacingValue(100, 200, 300, true), 70);

  frSpacingTablePrlConstraint table(fr2DLookupTbl<frCoord, frCoord, frCoord>(
      "WIDTH", {0}, "PRL", {0}, {{90}}));
  layer.setMinSpacing(&table);
  EXPECT_EQ(layer.getMinSpacingValue(100, 200, 300, false), 90);
  table.setLookupTbl(fr2DLookupTbl<frCoord, frCoord, frCoord>(
      "WIDTH", {0}, "PRL", {0}, {{110}}));
  EXPECT_EQ(layer.getMinSpacingValue(100, 200, 300, false), 110);
  EXPECT_EQ(layer.getMinSpacingValue(100, 200, 300, true), 110);

  // Unsupported constraint kinds retain the previous zero-base behavior.
  frSpacingSamenetConstraint same_net(200, false);
  layer.setMinSpacing(&same_net);
  EXPECT_EQ(layer.getMinSpacingValue(100, 200, 300, false), 0);
  layer.setMinSpacing(nullptr);
  EXPECT_FALSE(layer.hasMinSpacing());
  EXPECT_EQ(layer.getMinSpacingValue(100, 200, 300, true), 0);
  layer.setMinSpacing(&spacing);
  EXPECT_EQ(layer.getMinSpacingValue(100, 200, 300, false), 70);
}

TEST(LayerSpacing, OverlappingRangesAndMutation)
{
  frLayer layer;
  frSpacingRangeConstraint first;
  first.setMinWidth(100);
  first.setMaxWidth(200);
  first.setMinSpacing(70);
  frSpacingRangeConstraint second;
  second.setMinWidth(180);
  second.setMaxWidth(300);
  second.setMinSpacing(90);
  layer.addSpacingRangeConstraint(&first);
  layer.addSpacingRangeConstraint(&second);
  for (const auto width : {99, 100, 179, 180, 200, 201, 300, 301}) {
    const int expected = width < 100 || width > 300 ? 0 : width < 180 ? 70 : 90;
    EXPECT_EQ(layer.getMinSpacingValue(width, 0, 0, false), expected);
    EXPECT_EQ(layer.getMinSpacingValue(0, width, 0, false), expected);
    EXPECT_EQ(layer.getMinSpacingValue(width, 0, 0, true), 70);
  }
  frSpacingConstraint base(80);
  layer.setMinSpacing(&base);
  EXPECT_EQ(layer.getMinSpacingValue(150, 250, 0, false), 90);
  EXPECT_EQ(layer.getMinSpacingValue(0, 0, 0, false), 80);
  EXPECT_EQ(layer.getMinSpacingValue(0, 0, 0, true), 70);
  first.setMinSpacing(95);
  first.setMinWidth(200);
  first.setMaxWidth(220);
  EXPECT_EQ(layer.getMinSpacingValue(150, 0, 0, false), 80);
  EXPECT_EQ(layer.getMinSpacingValue(210, 0, 0, false), 95);
  EXPECT_EQ(layer.getMinSpacingValue(0, 0, 0, true), 80);
  layer.setMinSpacing(nullptr);
  EXPECT_EQ(layer.getMinSpacingValue(0, 0, 0, true), 90);
}

TEST(LayerSpacing, PrlThresholdsAndWidthSymmetry)
{
  const std::vector<frCoord> widths{0, 100, 200};
  const std::vector<frCoord> lengths{0, 50, 150};
  const std::vector<std::vector<frCoord>> values{
      {10, 20, 30}, {40, 50, 60}, {70, 80, 90}};
  frSpacingTablePrlConstraint rule(fr2DLookupTbl<frCoord, frCoord, frCoord>(
      "WIDTH", widths, "PRL", lengths, values));
  frLayer layer;
  layer.setMinSpacing(&rule);
  // LEF table thresholds select the previous row/column at equality.
  auto index = [](const auto& thresholds, int value) {
    int result = 0;
    for (int i = 0; i < static_cast<int>(thresholds.size()); ++i) {
      if (value > thresholds[i]) {
        result = i;
      }
    }
    return result;
  };
  for (const auto w1 : {-1, 0, 1, 99, 100, 101, 199, 200, 201, 1000}) {
    for (const auto w2 : {0, 100, 200, 201}) {
      for (const auto prl : {-1, 0, 1, 49, 50, 51, 149, 150, 151, 1000}) {
        const int expected
            = values[index(widths, std::max(w1, w2))][index(lengths, prl)];
        EXPECT_EQ(layer.getMinSpacingValue(w1, w2, prl, false), expected);
        EXPECT_EQ(layer.getMinSpacingValue(w2, w1, prl, false), expected);
        EXPECT_EQ(layer.getMinSpacingValue(w1, w2, prl, true), 10);
      }
    }
  }
}

TEST(LayerSpacing, TwoWidthsUsesLiveDatabaseTable)
{
  auto* db = odb::dbDatabase::create();
  auto* tech = odb::dbTech::create(db, "tech");
  auto* db_layer
      = odb::dbTechLayer::create(tech, "M1", odb::dbTechLayerType::ROUTING);
  db_layer->initTwoWidths(2);
  db_layer->addTwoWidthsIndexEntry(90);
  db_layer->addTwoWidthsIndexEntry(190, 100);
  db_layer->addTwoWidthsSpacingTableEntry(0, 0, 10);
  db_layer->addTwoWidthsSpacingTableEntry(0, 1, 50);
  db_layer->addTwoWidthsSpacingTableEntry(1, 0, 50);
  db_layer->addTwoWidthsSpacingTableEntry(1, 1, 150);
  frLayer layer;
  layer.setDbLayer(db_layer);
  frSpacingTableTwConstraint rule({{90, -1}, {190, 100}},
                                  {{10, 50}, {50, 150}});
  rule.setLayer(&layer);
  layer.setMinSpacing(&rule);
  for (const auto w1 : {0, 89, 90, 91, 189, 190, 191, 1000}) {
    for (const auto w2 : {0, 90, 190, 191, 1000}) {
      for (const auto prl : {-1, 0, 99, 100, 101, 1000}) {
        EXPECT_EQ(layer.getMinSpacingValue(w1, w2, prl, false),
                  db_layer->findTwSpacing(w1, w2, prl));
        EXPECT_EQ(layer.getMinSpacingValue(w1, w2, prl, true), 10);
      }
    }
  }
  // TW lookup reads OpenDB, while USEMINSPACING reads the rule's own minimum.
  db_layer->addTwoWidthsSpacingTableEntry(1, 1, 175);
  EXPECT_EQ(layer.getMinSpacingValue(1000, 1000, 1000, false), 175);
  rule.setSpacingTable({{90, -1}, {190, 100}}, {{20, 50}, {50, 175}});
  EXPECT_EQ(layer.getMinSpacingValue(1000, 1000, 1000, true), 20);
  layer.setMinSpacing(nullptr);
  EXPECT_EQ(layer.getMinSpacingValue(1000, 1000, 1000, false), 0);
  odb::dbDatabase::destroy(db);
}

}  // namespace
}  // namespace drt
