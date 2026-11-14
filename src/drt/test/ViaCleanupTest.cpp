// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include <set>
#include <tuple>

#include "dr/ViaCleanup.h"
#include "fixture.h"
#include "gtest/gtest.h"

namespace drt {
namespace {

class ViaCleanupTest : public Fixture
{
 protected:
  ViaCleanupTest()
  {
    addLayer(design->getTech(), "v1", odb::dbTechLayerType::CUT);
    addLayer(design->getTech(),
             "m2",
             odb::dbTechLayerType::ROUTING,
             odb::dbTechLayerDir::VERTICAL);
    net_ = makeNet("signal");
    auto* def = makeViaDef("via", 3, {-80, -80}, {80, 80});
    makePathseg(net_,
                4,
                {400, 500},
                {1000, 500},
                120,
                frcExtendEndStyle,
                frcTruncateEndStyle);
    candidate_
        = {static_cast<frPathSeg*>(net_->getShapes().front().get()),
           {makeVia(def, net_, {400, 500}), makeVia(def, net_, {1000, 500})},
           2,
           600};
    net_->setModified(false);
  }

  using Geometry
      = std::set<std::tuple<int, int, int, int, int, frBlockObject*>>;

  Geometry geometry()
  {
    Geometry result;
    for (int layer = 2; layer <= 4; ++layer) {
      frRegionQuery::Objects<frBlockObject> objects;
      design->getRegionQuery()->queryDRObj({0, 0, 2000, 2000}, layer, objects);
      for (const auto& [box, object] : objects) {
        result.emplace(
            layer, box.xMin(), box.yMin(), box.xMax(), box.yMax(), object);
      }
    }
    return result;
  }

  int markerCount()
  {
    FlexGCWorker checker(design->getTech(), logger.get(), router_cfg.get());
    checker.setExtBox({0, 0, 2000, 2000});
    checker.setDrcBox({0, 0, 2000, 2000});
    checker.setEnableSurgicalFix(false);
    checker.init(design.get());
    checker.main();
    return checker.getMarkers().size();
  }

  frSegStyle targetStyle()
  {
    frSegStyle style;
    style.setWidth(100);
    style.setBeginStyle(frcTruncateEndStyle);
    style.setEndStyle(frcTruncateEndStyle);
    return style;
  }

  frNet* net_;
  via_cleanup::Candidate candidate_;
};

TEST_F(ViaCleanupTest, RejectedTrialRestoresGeometryAndOwnership)
{
  auto* obstruction = makeNet("other_signal");
  makePathseg(obstruction, 2, {700, 300}, {700, 700});
  initRegionQuery();
  ASSERT_EQ(markerCount(), 0);
  ASSERT_TRUE(via_cleanup::isolated(design->getRegionQuery(), candidate_));
  const auto before = geometry();
  const auto original_box = candidate_.segment->getBBox();
  const auto original_index = candidate_.segment->getIndexInOwner();
  {
    via_cleanup::Trial trial(
        design->getRegionQuery(), candidate_, targetStyle());
    EXPECT_GT(markerCount(), 0);
  }
  EXPECT_EQ(markerCount(), 0);
  EXPECT_EQ(geometry(), before);
  EXPECT_EQ(candidate_.segment->getLayerNum(), 4);
  EXPECT_EQ(candidate_.segment->getBBox(), original_box);
  EXPECT_EQ(candidate_.segment->getIndexInOwner(), original_index);
  EXPECT_EQ(net_->getShapes().front().get(), candidate_.segment);
  EXPECT_EQ(net_->getVias().front().get(), candidate_.vias[0]);
  EXPECT_EQ(net_->getVias().back().get(), candidate_.vias[1]);
  EXPECT_TRUE(net_->getPatchWires().empty());
  EXPECT_FALSE(net_->isModified());
}

TEST_F(ViaCleanupTest, RejectsInteriorSourceBranch)
{
  makePathseg(net_, 4, {700, 500}, {700, 900});
  initRegionQuery();
  EXPECT_FALSE(via_cleanup::isolated(design->getRegionQuery(), candidate_));
}

TEST_F(ViaCleanupTest, RejectsStackedViaOnOppositeCutLayer)
{
  addLayer(design->getTech(), "v2", odb::dbTechLayerType::CUT);
  addLayer(design->getTech(), "m3", odb::dbTechLayerType::ROUTING);
  auto* def = makeViaDef("upper_via", 5, {-80, -80}, {80, 80});
  makeVia(def, net_, {400, 500});
  initRegionQuery();
  EXPECT_FALSE(via_cleanup::isolated(design->getRegionQuery(), candidate_));
}

TEST_F(ViaCleanupTest, RejectsInteriorViaOnSameCutLayer)
{
  auto* def = makeViaDef("interior_via", 3, {-80, -80}, {80, 80});
  makeVia(def, net_, {700, 500});
  initRegionQuery();
  EXPECT_FALSE(via_cleanup::isolated(design->getRegionQuery(), candidate_));
}

TEST_F(ViaCleanupTest, RejectsOffCenterViaMetalContact)
{
  auto* def = makeViaDef("offset_via", 3, {-20, -20}, {20, 20});
  auto metal = std::make_unique<frRect>();
  metal->setLayerNum(4);
  metal->setBBox({-60, -160, 60, 20});
  def->addLayer2Fig(std::move(metal));
  // The origin and cut miss the source wire; only the enclosure touches it.
  makeVia(def, net_, {700, 700});
  initRegionQuery();
  EXPECT_FALSE(via_cleanup::isolated(design->getRegionQuery(), candidate_));
}

TEST_F(ViaCleanupTest, IgnoresUnrelatedRoutingLayer)
{
  makePathseg(makeNet("other_layer"), 2, {700, 300}, {700, 700});
  initRegionQuery();
  EXPECT_TRUE(via_cleanup::isolated(design->getRegionQuery(), candidate_));
}

TEST_F(ViaCleanupTest, RejectsFixedTargetInteriorContact)
{
  const auto [master, db_master] = makeMacro("pin", 0, 0, 2000, 2000);
  makeMacroPin(master, "interior_pin", 650, 450, 750, 550, 2);
  makeInst("pin_inst", master, db_master);
  initRegionQuery();
  EXPECT_FALSE(via_cleanup::isolated(design->getRegionQuery(), candidate_));
}

TEST_F(ViaCleanupTest, RejectsTargetInteriorJunction)
{
  makePathseg(net_, 2, {700, 300}, {700, 700});
  initRegionQuery();
  EXPECT_FALSE(via_cleanup::targetClear(design->getRegionQuery(), candidate_));
}

TEST_F(ViaCleanupTest, RejectsTargetCollinearOverlap)
{
  makePathseg(net_, 2, {700, 500}, {1200, 500});
  initRegionQuery();
  EXPECT_FALSE(via_cleanup::targetClear(design->getRegionQuery(), candidate_));
}

TEST_F(ViaCleanupTest, RejectsTargetInteriorVia)
{
  auto* def = makeViaDef("lower_via", 1, {-80, -80}, {80, 80});
  makeVia(def, net_, {700, 500});
  initRegionQuery();
  EXPECT_FALSE(via_cleanup::targetClear(design->getRegionQuery(), candidate_));
}

TEST_F(ViaCleanupTest, AllowsTargetEndpointJunctions)
{
  makePathseg(net_, 2, {400, 300}, {400, 700});
  makePathseg(net_, 2, {1000, 500}, {1200, 500});
  initRegionQuery();
  EXPECT_TRUE(via_cleanup::targetClear(design->getRegionQuery(), candidate_));
}

TEST_F(ViaCleanupTest, RejectsFixedSourceContact)
{
  const auto [master, db_master] = makeMacro("fixed", 0, 0, 2000, 2000);
  makeMacroObs(master, 650, 450, 750, 550, 4);
  makeInst("fixed_inst", master, db_master);
  initRegionQuery();
  EXPECT_FALSE(via_cleanup::isolated(design->getRegionQuery(), candidate_));
}

TEST_F(ViaCleanupTest, RejectsFixedTargetContact)
{
  const auto [master, db_master] = makeMacro("pin", 0, 0, 2000, 2000);
  makeMacroPin(master, "landing_pin", 350, 450, 450, 550, 2);
  makeInst("pin_inst", master, db_master);
  initRegionQuery();
  EXPECT_FALSE(via_cleanup::isolated(design->getRegionQuery(), candidate_));
}

TEST_F(ViaCleanupTest, CommitPreservesOffCenterLandingContacts)
{
  // These wires contact the old landing metal above the replacement wire.
  makePathseg(net_, 2, {400, 570}, {400, 900});
  makePathseg(net_, 2, {1000, 570}, {1000, 900});
  initRegionQuery();
  ASSERT_TRUE(via_cleanup::isolated(design->getRegionQuery(), candidate_));
  ASSERT_EQ(markerCount(), 0);
  const std::set<odb::Rect> landings{candidate_.vias[0]->getLayer1BBox(),
                                     candidate_.vias[1]->getLayer1BBox()};
  {
    via_cleanup::Trial trial(
        design->getRegionQuery(), candidate_, targetStyle());
    ASSERT_EQ(markerCount(), 0);
    trial.commit();
  }
  EXPECT_TRUE(net_->getVias().empty());
  EXPECT_EQ(candidate_.segment->getLayerNum(), 2);
  EXPECT_TRUE(net_->isModified());
  std::set<odb::Rect> retained;
  for (const auto& patch : net_->getPatchWires()) {
    EXPECT_EQ(patch->getNet(), net_);
    EXPECT_EQ(patch->getLayerNum(), 2);
    retained.insert(patch->getBBox());
  }
  EXPECT_EQ(retained, landings);
  frRegionQuery::Objects<frBlockObject> removed_layer;
  design->getRegionQuery()->queryDRObj({0, 0, 2000, 2000}, 4, removed_layer);
  EXPECT_TRUE(removed_layer.empty());
  EXPECT_EQ(markerCount(), 0);
}

}  // namespace
}  // namespace drt
