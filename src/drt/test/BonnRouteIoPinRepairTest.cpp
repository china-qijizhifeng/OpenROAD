// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include <memory>
#include <vector>

#include "br/BrIoPinArea.h"
#include "br/BrIoPinRepair.h"
#include "db/obj/frBPin.h"
#include "db/obj/frBTerm.h"
#include "db/obj/frBlock.h"
#include "db/obj/frBoundary.h"
#include "db/obj/frNet.h"
#include "db/obj/frShape.h"
#include "db/tech/frConstraint.h"
#include "db/tech/frLayer.h"
#include "fixture.h"
#include "gtest/gtest.h"

namespace drt {
namespace {

class IoPinRepairTest : public Fixture
{
 protected:
  void SetUp() override
  {
    frBoundary boundary;
    boundary.setPoints({{0, 0}, {5000, 0}, {5000, 5000}, {0, 5000}});
    design->getTopBlock()->setBoundaries({boundary});
    auto area = std::make_unique<frAreaConstraint>(20000);
    design->getTech()->getLayer(2)->setAreaConstraint(area.get());
    design->getTech()->addUConstraint(std::move(area));
    makeSimpleSpacingConstraint(2, 100);
  }

  void pin(frNet* net, const odb::Rect& box)
  {
    auto term = std::make_unique<frBTerm>("io");
    term->addToNet(net);
    auto port = std::make_unique<frBPin>();
    auto shape = std::make_unique<frRect>();
    shape->setLayerNum(2);
    shape->setBBox(box);
    port->addPinFig(std::move(shape));
    term->addPin(std::move(port));
    net->addBTerm(term.get());
    design->getTopBlock()->addTerm(std::move(term));
  }

  int64_t area(const frNet* net)
  {
    std::vector<br::Box> boxes;
    for (const auto& shape : net->getShapes()) {
      const odb::Rect b = shape->getBBox();
      boxes.push_back({b.xMin(), b.yMin(), b.xMax(), b.yMax()});
    }
    for (const auto& patch : net->getPatchWires()) {
      const odb::Rect b = patch->getBBox();
      boxes.push_back({b.xMin(), b.yMin(), b.xMax(), b.yMax()});
    }
    return br::unionArea(boxes);
  }
};

TEST_F(IoPinRepairTest, PreservesRoutedWiresAndIsIdempotent)
{
  frNet* io = makeNet("io");
  pin(io, odb::Rect(0, 950, 100, 1050));
  makePathsegExt(io, 2, {50, 1000}, {60, 1000});
  frNet* internal = makeNet("internal");
  makePathsegExt(internal, 2, {2000, 2000}, {3000, 2000});
  const odb::Rect io_before = io->getShapes().front()->getBBox();
  const odb::Rect internal_before = internal->getShapes().front()->getBBox();
  initRegionQuery();
  br::repairIoPinArea(design.get(), logger.get(), router_cfg.get());
  ASSERT_EQ(io->getPatchWires().size(), 1);
  EXPECT_GE(area(io), 20000);
  EXPECT_EQ(io->getShapes().front()->getBBox(), io_before);
  EXPECT_EQ(internal->getShapes().front()->getBBox(), internal_before);
  EXPECT_TRUE(internal->getPatchWires().empty());
  br::repairIoPinArea(design.get(), logger.get(), router_cfg.get());
  EXPECT_EQ(io->getPatchWires().size(), 1);
}

TEST_F(IoPinRepairTest, DoesNotJoinSeparateComponentsThroughPinPurpose)
{
  frNet* io = makeNet("io");
  pin(io, odb::Rect(0, 700, 100, 1300));
  makePathsegExt(io, 2, {50, 800}, {60, 800});
  makePathsegExt(io, 2, {50, 1200}, {60, 1200});
  initRegionQuery();
  br::repairIoPinArea(design.get(), logger.get(), router_cfg.get());
  ASSERT_EQ(io->getPatchWires().size(), 2);
  EXPECT_GE(area(io), 40000);
}

TEST_F(IoPinRepairTest, RejectsExtensionTowardDifferentNet)
{
  frNet* io = makeNet("io");
  pin(io, odb::Rect(0, 950, 100, 1050));
  makePathsegExt(io, 2, {50, 1000}, {60, 1000});
  frNet* blocker = makeNet("blocker");
  makePathsegExt(blocker, 2, {300, 1000}, {400, 1000});
  initRegionQuery();
  br::repairIoPinArea(design.get(), logger.get(), router_cfg.get());
  EXPECT_TRUE(io->getPatchWires().empty());
  EXPECT_TRUE(blocker->getPatchWires().empty());
  EXPECT_LT(area(io), 20000);
}

}  // namespace
}  // namespace drt
