// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include <algorithm>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

#include "LocalAntennaRepair.h"
#include "RepairAntennas.h"
#include "ant/AntennaChecker.hh"
#include "db_sta/dbSta.hh"
#include "dpl/Opendp.h"
#include "grt/GRoute.h"
#include "grt/GlobalRouter.h"
#include "gtest/gtest.h"
#include "odb/db.h"
#include "odb/defin.h"
#include "stt/SteinerTreeBuilder.h"
#include "tst/db_fixture.h"
#include "tst/fixture.h"
#include "tst/physical_builder.h"
#include "utl/ServiceRegistry.h"

namespace grt {

class LocalAntennaRepairTestPeer
{
 public:
  static void refreshUsage(GlobalRouter& router)
  {
    router.updateDbCongestion();
  }

  static bool rejectedByAntenna(RepairAntennas& repair,
                                GRoute& route,
                                int segment,
                                int position,
                                odb::dbNet* net,
                                int tile_size)
  {
    repair.tile_size_ = tile_size;
    repair.jumper_size_ = 2 * tile_size;
    repair.smaller_seg_size_ = 5 * tile_size;
    auto before = repair.checkLocalAntenna(net, nullptr, 0);
    // Force rejection even if this is already a clean net: the candidate
    // must reach the real checker, which cannot improve on zero demand.
    before.violations.clear();
    return repair.tryLocalJumper(
               route, segment, position, net, nullptr, 0, before)
           == RepairAntennas::LocalJumperResult::Antenna;
  }
};

namespace {

class LocalAntennaCostTest : public tst::DbFixture
{
 protected:
  void SetUp() override
  {
    auto* tech = tst::makeTech(db_.get());
    auto layers = tst::makeLayerStack(tech, 3);
    auto* lib = tst::makeLib(db_.get(), tech);
    auto* block = tst::makeBlock(db_.get(), tech);
    auto* master = tst::makeMaster(lib,
                                   "gate",
                                   1000,
                                   1000,
                                   {{.name = "A",
                                     .layer = layers.metal(1),
                                     .rect = {0, 0, 100, 100},
                                     .gate_area = 1.0}});
    auto* gate = makeInst(block, master, "gate")->findITerm("A");
    auto* other_gate = makeInst(block, master, "other_gate")->findITerm("A");
    before_.pins = {{.gate = gate,
                     .layers = {{.layer = layers.metal(2),
                                 .par = {200.0, 100.0, true},
                                 .car = {},
                                 .psr = {},
                                 .csr = {},
                                 .violated = true}},
                     .violated = true},
                    {.gate = other_gate, .layers = {}, .violated = false}};
    before_.violated = true;
    upper_ = layers.metal(3);
  }

  ant::NetCheck before_;
  odb::dbTechLayer* upper_ = nullptr;
};

TEST_F(LocalAntennaCostTest, RequiresStrictDiodeDemandReduction)
{
  auto after = before_;
  after.pins[0].layers[0].par.actual = 150.0;
  EXPECT_FALSE(localAntennaImproves(before_, after, 2, 2));
  EXPECT_TRUE(localAntennaImproves(before_, after, 2, 1));
}

TEST_F(LocalAntennaCostTest, RejectsNewLayerViolation)
{
  auto after = before_;
  after.pins[0].layers.push_back({.layer = upper_,
                                  .par = {101.0, 100.0, true},
                                  .car = {},
                                  .psr = {},
                                  .csr = {},
                                  .violated = true});
  EXPECT_FALSE(localAntennaImproves(before_, after, 2, 1));
}

TEST_F(LocalAntennaCostTest, RejectsNewRatioFamilyAndWorseExistingRatio)
{
  auto after = before_;
  after.pins[0].layers[0].car = {101.0, 100.0, true};
  EXPECT_FALSE(localAntennaImproves(before_, after, 2, 1));
  after = before_;
  after.pins[0].layers[0].par.actual = 201.0;
  EXPECT_FALSE(localAntennaImproves(before_, after, 2, 1));
  after.pins[0].layers[0].par.actual = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(localAntennaImproves(before_, after, 2, 1));
}

TEST_F(LocalAntennaCostTest, RejectsMissingOrNewlyViolatingGate)
{
  auto after = before_;
  after.pins.pop_back();
  EXPECT_FALSE(localAntennaImproves(before_, after, 2, 0));
  after = before_;
  after.pins[1].layers = before_.pins[0].layers;
  EXPECT_FALSE(localAntennaImproves(before_, after, 2, 1));
}

TEST_F(LocalAntennaCostTest, AcceptsClearedViolationWithAllGatesPresent)
{
  auto after = before_;
  after.pins[0].layers.clear();
  after.pins[0].violated = false;
  after.violated = false;
  EXPECT_TRUE(localAntennaImproves(before_, after, 2, 0));
}

TEST_F(LocalAntennaCostTest, ZeroDiodeCountIsStillARealViolation)
{
  EXPECT_EQ(antennaDiodeDemand({{.routing_level = 2,
                                 .gates = {before_.pins[0].gate},
                                 .diode_count_per_gate = 0,
                                 .excess_ratio = 2.0}}),
            1);
  EXPECT_EQ(antennaDiodeDemand({}), 0);
}

using Guide = std::tuple<int, int, int, int, int, int, bool, bool>;

std::vector<Guide> guides(odb::dbNet* net)
{
  std::vector<Guide> result;
  for (auto* guide : net->getGuides()) {
    const auto box = guide->getBox();
    result.emplace_back(guide->getLayer()->getId(),
                        guide->getViaLayer()->getId(),
                        box.xMin(),
                        box.yMin(),
                        box.xMax(),
                        box.yMax(),
                        guide->isJumper(),
                        guide->isConnectedToTerm());
  }
  std::sort(result.begin(), result.end());
  return result;
}

class LocalAntennaRollbackTest : public tst::Fixture,
                                 public testing::WithParamInterface<bool>
{
 protected:
  void SetUp() override
  {
    readLiberty("_main/test/Nangate45/Nangate45_typ.lib");
    auto* lib = loadTechAndLib(
        "nangate45", "nangate45", "_main/test/Nangate45/Nangate45.lef");
    ASSERT_NE(lib, nullptr);
    auto* chip = odb::dbChip::create(db_.get(), lib->getTech());
    odb::defin reader(db_.get(), &logger_);
    std::vector<odb::dbLib*> libs{lib};
    ASSERT_TRUE(reader.readChip(
        libs,
        getFilePath("_main/src/grt/test/remove_buffers1.def").c_str(),
        chip));
    block_ = chip->getBlock();
    sta_->postReadDef(block_);
    router_.setUseCUGR(GetParam());
    router_.setMinRoutingLayer(2);
    router_.setMaxRoutingLayer(8);
    router_.globalRoute(/*save_guides=*/true);
  }

  std::vector<float> usage()
  {
    LocalAntennaRepairTestPeer::refreshUsage(router_);
    int nx, ny;
    router_.getGridSize(nx, ny);
    std::vector<float> result;
    for (int layer = 2; layer <= 8; layer++) {
      auto* metal = db_->getTech()->findRoutingLayer(layer);
      for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
          result.push_back(block_->getGCellGrid()->getUsage(metal, x, y));
        }
      }
    }
    return result;
  }

  odb::dbBlock* block_ = nullptr;
  utl::ServiceRegistry registry_{&logger_};
  stt::SteinerTreeBuilder stt_{&logger_};
  ant::AntennaChecker ant_{db_.get(), &logger_};
  dpl::Opendp dpl_{db_.get(), &logger_};
  GlobalRouter
      router_{&logger_, &registry_, &stt_, db_.get(), sta_.get(), &ant_, &dpl_};
};

TEST_P(LocalAntennaRollbackTest,
       AntennaRejectionPreservesRouteGuidesAndResources)
{
  RepairAntennas repair(&router_, &ant_, &dpl_, db_.get(), &logger_);
  const auto resources = usage();
  const int tile = router_.getTileSize();
  for (auto& [net, route] : router_.getRoutes()) {
    const GRoute original = route;
    const auto original_guides = guides(net);
    const bool had_jumpers = net->hasJumpers();
    for (int segment = 0; segment < original.size(); segment++) {
      const auto& seg = original[segment];
      if (seg.isVia() || seg.init_layer + 2 > 8 || seg.length() < 5 * tile) {
        continue;
      }
      const int lo = seg.init_x != seg.final_x ? seg.init_x : seg.init_y;
      const bool checked = LocalAntennaRepairTestPeer::rejectedByAntenna(
          repair, route, segment, lo + tile, net, tile);
      EXPECT_EQ(route, original);
      EXPECT_EQ(guides(net), original_guides);
      EXPECT_EQ(net->hasJumpers(), had_jumpers);
      EXPECT_EQ(net->getWire(), nullptr);
      EXPECT_EQ(usage(), resources);
      if (checked) {
        return;
      }
    }
  }
  FAIL() << "No legal trial reached the antenna checker";
}

INSTANTIATE_TEST_SUITE_P(RoutingEngines,
                         LocalAntennaRollbackTest,
                         testing::Bool(),
                         [](const testing::TestParamInfo<bool>& info) {
                           return info.param ? "CUGR" : "FastRoute";
                         });

}  // namespace
}  // namespace grt
