// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// BonnRoute global routing on a small synthetic design: the pin model of
// Gester et al., TODAES 2013, Sec. 2.1-2.2 (pins as vertex sets joined by
// zero-cost cliques) and the capacity estimate of Sec. 2.5.

#include <array>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "BonnRoute.h"
#include "BonnRouteCapacity.h"
#include "CUGR.h"
#include "Design.h"
#include "GRNet.h"
#include "GRTree.h"
#include "GridGraph.h"
#include "geo.h"
#include "gtest/gtest.h"
#include "odb/PtrSetMap.h"
#include "odb/db.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"
#include "utl/Logger.h"

namespace grt {
namespace {

// 9000 x 6000 die, M1 (V) / M2 (H) / M3 (V), all with a 100 pitch, tracks at
// multiples of 100. The gcell is 15 pitches; the last row and column absorb
// the remainder: x gridlines 0, 1500, 3000, 4500, 6000, 9000 and y gridlines
// 0, 1500, 3000, 6000 (5 x 3 cells, centers x 750 .. 5250, 7500 and y 750,
// 2250, 4500). Routing layers M2-M3.
class BonnRouteModels : public testing::Test
{
 protected:
  void SetUp() override
  {
    db_.reset(odb::dbDatabase::create());
    db_->setLogger(&logger_);
    tech_ = odb::dbTech::create(db_.get(), "test_tech");
    auto* chip = odb::dbChip::create(db_.get(), tech_);
    block_ = odb::dbBlock::create(chip, "test_block");
    block_->setDieArea(odb::Rect(0, 0, 9000, 6000));
    block_->setMinRoutingLayer(2);
    block_->setMaxRoutingLayer(3);
    for (int z = 0; z < 3; z++) {
      if (z > 0) {
        const std::string cut = "V" + std::to_string(z);
        odb::dbTechLayer::create(tech_, cut.c_str(), odb::dbTechLayerType::CUT);
      }
      const std::string name = "M" + std::to_string(z + 1);
      auto* metal = odb::dbTechLayer::create(
          tech_, name.c_str(), odb::dbTechLayerType::ROUTING);
      metal->setDirection(z == 1 ? odb::dbTechLayerDir::HORIZONTAL
                                 : odb::dbTechLayerDir::VERTICAL);
      metal->setWidth(40);
      metal->setSpacing(40);
      metal->setPitch(100);
      layers_.push_back(metal);
      auto* tracks = odb::dbTrackGrid::create(block_, metal);
      tracks->addGridPatternX(0, 91, 100);
      tracks->addGridPatternY(0, 61, 100);
    }
  }

  // Pins are (layer index, x0, y0, x1, y1) BTerm shapes; pin 0 drives.
  odb::dbNet* addNet(const std::string& name,
                     const std::vector<std::array<int, 5>>& pins)
  {
    auto* net = odb::dbNet::create(block_, name.c_str());
    for (size_t p = 0; p < pins.size(); p++) {
      const auto [layer, x0, y0, x1, y1] = pins[p];
      const std::string term_name = name + "_p" + std::to_string(p);
      auto* term = odb::dbBTerm::create(net, term_name.c_str());
      term->setIoType(p == 0 ? odb::dbIoType::INPUT : odb::dbIoType::OUTPUT);
      auto* pin = odb::dbBPin::create(term);
      pin->setPlacementStatus(odb::dbPlacementStatus::PLACED);
      odb::dbBox::create(pin, layers_[layer], x0, y0, x1, y1);
    }
    return net;
  }

  // Builds Design/GridGraph after the nets exist.
  void build()
  {
    design_ = std::make_unique<Design>(db_.get(),
                                       &logger_,
                                       constants_,
                                       2,
                                       3,
                                       odb::PtrSet<odb::dbNet>{},
                                       false);
    grid_ = std::make_unique<GridGraph>(design_.get(), constants_, &logger_);
    nets_.clear();
    for (const CUGRNet& base : design_->getAllNets()) {
      nets_.push_back(std::make_unique<GRNet>(base, grid_.get()));
    }
  }

  void route(const int capacity_model)
  {
    BonnRouteParams params;
    params.phases = 1;
    params.repair_passes = 0;
    params.verbose = false;
    params.full_graph = true;
    params.capacity_model = capacity_model;
    std::vector<GRNet*> nets;
    for (auto& net : nets_) {
      nets.push_back(net.get());
    }
    BonnRouteGlobal router(grid_.get(), nets, params, &logger_);
    router.run();
  }

  // Wire edges (layer, lower cell x, y) carrying demand.
  std::set<std::array<int, 3>> wireEdges(const GRNet& net) const
  {
    std::set<std::array<int, 3>> result;
    grid_->forEachTreeDemand(
        net.getRoutingTree(), {}, [&](int z, PointT p, double d) {
          if (d > 0) {
            result.insert({z, p.x(), p.y()});
          }
        });
    return result;
  }

  // GRNet pin index of the BTerm `name`.
  static int pinIndex(const GRNet& net, const std::string& name)
  {
    for (const auto& [index, bterm] : net.getBTermsByPinIndex()) {
      if (bterm->getName() == name) {
        return index;
      }
    }
    return -1;
  }

  int pinLinks(const GRNet& net, int* pin = nullptr) const
  {
    int links = 0;
    GRTreeNode::preorder(net.getRoutingTree(), [&](const auto& node) {
      for (const auto& child : node->getChildren()) {
        if (child->isPinLink()) {
          links++;
          if (pin != nullptr) {
            *pin = child->getLinkPin();
          }
        }
      }
    });
    return links;
  }

  utl::Logger logger_;
  std::unique_ptr<odb::dbDatabase, decltype(&odb::dbDatabase::destroy)> db_{
      nullptr,
      odb::dbDatabase::destroy};
  odb::dbTech* tech_ = nullptr;
  odb::dbBlock* block_ = nullptr;
  std::vector<odb::dbTechLayer*> layers_;
  Constants constants_;
  std::unique_ptr<Design> design_;
  std::unique_ptr<GridGraph> grid_;
  std::vector<std::unique_ptr<GRNet>> nets_;
};

class BonnRoutePinModel : public BonnRouteModels,
                          public testing::WithParamInterface<int>
{
};

TEST_P(BonnRoutePinModel, BranchesAttachToDifferentCellsOfOnePin)
{
  // The driver's M2 shape covers cells 1..3 of row 1. Each sink attaches to
  // the nearer end; the pin itself joins cells 1 and 3, so no wire runs
  // between them.
  addNet("n",
         {{1, 1600, 1600, 5900, 1640},
          {1, 100, 1600, 140, 1640},
          {1, 6200, 1600, 6240, 1640}});
  build();
  route(GetParam());
  const GRNet& net = *nets_.front();
  ASSERT_NE(net.getRoutingTree(), nullptr);
  const std::set<std::array<int, 3>> expected{{1, 0, 1}, {1, 3, 1}};
  EXPECT_EQ(wireEdges(net), expected);
  EXPECT_EQ(grid_->getTreeLength(net.getRoutingTree()), 2);
  int pin = -1;
  EXPECT_EQ(pinLinks(net, &pin), 1);
  EXPECT_EQ(pin, pinIndex(net, "n_p0"));
  EXPECT_FALSE(net.isLocal());

  // Demand round trip: the pin link commits nothing.
  const auto committed = grid_->snapshotDemand();
  grid_->removeTreeUsage(net);
  EXPECT_DOUBLE_EQ(grid_->getEdge(1, 0, 1).demand, 0.0);
  EXPECT_DOUBLE_EQ(grid_->getEdge(1, 1, 1).demand, 0.0);
  EXPECT_DOUBLE_EQ(grid_->getEdge(1, 2, 1).demand, 0.0);
  EXPECT_DOUBLE_EQ(grid_->getEdge(1, 3, 1).demand, 0.0);
  grid_->addTreeUsage(net);
  EXPECT_EQ(grid_->snapshotDemand(), committed);

  // Every pin link becomes a bridge whose corner points move along one
  // axis at a time from parent to child.
  GRTreeNode::preorder(net.getRoutingTree(), [&](const auto& node) {
    for (const auto& child : node->getChildren()) {
      if (!child->isPinLink()) {
        continue;
      }
      const std::vector<GRPoint> path = pinLinkPath(*node, *child, 1, 2);
      ASSERT_GE(path.size(), size_t{2});
      EXPECT_EQ(path.front().x(), node->x());
      EXPECT_EQ(path.back().x(), child->x());
      EXPECT_EQ(path.back().getLayerIdx(), child->getLayerIdx());
      for (size_t i = 0; i + 1 < path.size(); i++) {
        const int changes
            = (path[i].x() != path[i + 1].x())
              + (path[i].y() != path[i + 1].y())
              + (path[i].getLayerIdx() != path[i + 1].getLayerIdx());
        EXPECT_EQ(changes, 1);
      }
    }
  });
}

TEST_P(BonnRoutePinModel, UnconnectedPinBridgesTwoBranches)
{
  // The long pin is a sink; reaching it at cell 1 makes all of its cells
  // sources, so the last sink connects from cell 3.
  addNet("n",
         {{1, 100, 1600, 140, 1640},
          {1, 1600, 1600, 5900, 1640},
          {1, 6200, 1600, 6240, 1640}});
  build();
  route(GetParam());
  const GRNet& net = *nets_.front();
  const std::set<std::array<int, 3>> expected{{1, 0, 1}, {1, 3, 1}};
  EXPECT_EQ(wireEdges(net), expected);
  int pin = -1;
  EXPECT_EQ(pinLinks(net, &pin), 1);
  EXPECT_EQ(pin, pinIndex(net, "n_p1"));
}

TEST_P(BonnRoutePinModel, OverlappingPinsNeedNoEdge)
{
  addNet("n", {{1, 1600, 1600, 4400, 1640}, {1, 3100, 1600, 5900, 1640}});
  build();
  route(GetParam());
  const GRNet& net = *nets_.front();
  ASSERT_NE(net.getRoutingTree(), nullptr);
  EXPECT_TRUE(wireEdges(net).empty());
  EXPECT_EQ(pinLinks(net), 0);  // unused clique edges are pruned
  EXPECT_TRUE(net.getRoutingTree()->getChildren().empty());
}

TEST(BonnRoutePinBridge, PlanarRunsUseTheGivenRoutableLayers)
{
  // Two M1 cells of one pin, two cells apart in x and one in y: M1 takes
  // only via access (one gcell), the x run is on M2 and the y run on M3.
  const std::vector<GRPoint> path
      = pinLinkPath(GRPoint(0, 1, 1), GRPoint(0, 3, 2), 1, 2);
  const std::vector<std::array<int, 3>> expected{
      {0, 1, 1}, {1, 1, 1}, {1, 3, 1}, {2, 3, 1}, {2, 3, 2}, {0, 3, 2}};
  ASSERT_EQ(path.size(), expected.size());
  for (size_t i = 0; i < path.size(); i++) {
    EXPECT_EQ(path[i].getLayerIdx(), expected[i][0]);
    EXPECT_EQ(path[i].x(), expected[i][1]);
    EXPECT_EQ(path[i].y(), expected[i][2]);
  }
  for (size_t i = 0; i + 1 < path.size(); i++) {
    if (path[i].getLayerIdx() == path[i + 1].getLayerIdx()) {
      EXPECT_GE(path[i].getLayerIdx(), 1);
      // x runs on the horizontal layer, y runs on the vertical one.
      EXPECT_EQ(path[i].getLayerIdx(), path[i].x() != path[i + 1].x() ? 1 : 2);
    }
  }
  // Same cell, different layers: a via stack only.
  const std::vector<GRPoint> stack
      = pinLinkPath(GRPoint(0, 2, 2), GRPoint(1, 2, 2), 1, 2);
  ASSERT_EQ(stack.size(), size_t{2});
}

INSTANTIATE_TEST_SUITE_P(CapacityModels,
                         BonnRoutePinModel,
                         testing::Values(0, 1));

TEST_F(BonnRouteModels, UnblockedEstimateCountsTracks)
{
  build();
  bonnroute::CapacityEstimator estimator(capacityGridFromCugr(*grid_), {});
  estimator.buildTrackGraph();
  const std::vector<double> capacity = estimator.capacities();
  for (int l = 1; l < 3; l++) {
    for (int x = 0; x < grid_->getXSize(); x++) {
      for (int y = 0; y < grid_->getYSize(); y++) {
        EXPECT_DOUBLE_EQ(capacity[grid_->edgeFlatIndex(l, x, y)],
                         grid_->getEdge(l, x, y).capacity)
            << l << " " << x << " " << y;
      }
    }
  }
  EXPECT_DOUBLE_EQ(capacity[grid_->edgeFlatIndex(1, 0, 0)], 15.0);
  EXPECT_DOUBLE_EQ(capacity[grid_->edgeFlatIndex(0, 0, 0)], 0.0);
}

TEST_F(BonnRouteModels, BlockageRemovesTrackGraphVertices)
{
  // M2 blockage 2000..2300 x 1700..1900: keep-out 0 + 20 - 1 = 19 (no
  // spacing table) across and along the tracks, so tracks y = 1700..1900 lose
  // the vertices x = 2000..2300. Edge (0, 1) spans [750, 2250) with 15 vertices
  // per track (3 of them blocked on 3 tracks); edge (1, 1) loses x = 2300.
  odb::dbObstruction::create(block_, layers_[1], 2000, 1700, 2300, 1900);
  build();
  bonnroute::CapacityEstimator estimator(capacityGridFromCugr(*grid_), {});
  estimator.buildTrackGraph();
  const std::vector<double> capacity = estimator.capacities();
  EXPECT_DOUBLE_EQ(capacity[grid_->edgeFlatIndex(1, 0, 1)], 15.0 - 9.0 / 15);
  EXPECT_DOUBLE_EQ(capacity[grid_->edgeFlatIndex(1, 1, 1)], 15.0 - 3.0 / 15);
  EXPECT_DOUBLE_EQ(capacity[grid_->edgeFlatIndex(1, 2, 1)], 15.0);
  EXPECT_DOUBLE_EQ(capacity[grid_->edgeFlatIndex(1, 0, 0)], 15.0);
}

TEST_F(BonnRouteModels, SingleTileNetConsumesItsVertices)
{
  build();
  bonnroute::CapacityEstimator estimator(capacityGridFromCugr(*grid_), {});
  estimator.buildTrackGraph();
  bonnroute::EstimatorNet net;
  net.single_tile = true;
  net.tile_x = 1;
  net.tile_y = 1;
  net.min_layer = 1;
  net.max_layer = 2;
  net.pins.resize(2);
  net.pins[0].points = {{odb::Point(1600, 1600), 1}};
  net.pins[1].points = {{odb::Point(2400, 1600), 1}};
  estimator.addNets({net});
  const std::vector<double> capacity = estimator.capacities();
  EXPECT_EQ(estimator.totals().single_tile_routed, 1);
  EXPECT_EQ(estimator.totals().single_tile_vertices, 9);
  // Vertices 1600..2200 lie in edge (0, 1), 2300 and 2400 in edge (1, 1).
  EXPECT_DOUBLE_EQ(capacity[grid_->edgeFlatIndex(1, 0, 1)], 15.0 - 7.0 / 15);
  EXPECT_DOUBLE_EQ(capacity[grid_->edgeFlatIndex(1, 1, 1)], 15.0 - 2.0 / 15);
  EXPECT_NEAR(estimator.totals().single_tile[1], 9.0 / 15, 1e-12);
}

TEST_F(BonnRouteModels, IntraTileSteinerLengthReducesFlankEdges)
{
  build();
  bonnroute::CapacityEstimator estimator(capacityGridFromCugr(*grid_), {});
  estimator.buildTrackGraph();
  bonnroute::EstimatorNet net;
  net.min_layer = 1;
  net.max_layer = 2;
  net.pins.resize(2);
  net.pins[0].points = {{odb::Point(1600, 1600), 1}};
  net.pins[0].cell_x = 1;
  net.pins[0].cell_y = 1;
  net.pins[1].points = {{odb::Point(6800, 2200), 1}};
  net.pins[1].cell_x = 4;
  net.pins[1].cell_y = 1;
  estimator.addNets({net});
  const std::vector<double> capacity = estimator.capacities();
  // Tile (1, 1): pin to center (2250, 2250) is 650 + 650; half of each
  // goes to either flanking edge (M2 along x, M3 along y), divided by the
  // edge's length (1500, or 2250 next to the wide last row).
  EXPECT_NEAR(
      capacity[grid_->edgeFlatIndex(1, 0, 1)], 15.0 - 0.5 * 650 / 1500, 1e-12);
  EXPECT_NEAR(
      capacity[grid_->edgeFlatIndex(1, 1, 1)], 15.0 - 0.5 * 650 / 1500, 1e-12);
  EXPECT_NEAR(
      capacity[grid_->edgeFlatIndex(2, 1, 0)], 15.0 - 0.5 * 650 / 1500, 1e-12);
  EXPECT_NEAR(
      capacity[grid_->edgeFlatIndex(2, 1, 1)], 15.0 - 0.5 * 650 / 2250, 1e-12);
  // Tile (4, 1), center (7500, 2250): 700 + 50. The last column has one
  // flanking edge along x.
  EXPECT_NEAR(
      capacity[grid_->edgeFlatIndex(1, 3, 1)], 15.0 - 0.5 * 700 / 2250, 1e-12);
  EXPECT_NEAR(
      capacity[grid_->edgeFlatIndex(2, 4, 0)], 31.0 - 0.5 * 50 / 1500, 1e-12);
  EXPECT_NEAR(
      capacity[grid_->edgeFlatIndex(2, 4, 1)], 31.0 - 0.5 * 50 / 2250, 1e-12);
  EXPECT_NEAR(estimator.totals().intra_tile[1],
              2 * 0.5 * 650 / 1500 + 0.5 * 700 / 2250,
              1e-12);
}

TEST_F(BonnRouteModels, StackedViasReduceOnlyNearbyIntermediateEdges)
{
  build();
  bonnroute::CapacityEstimator estimator(capacityGridFromCugr(*grid_), {});
  estimator.buildTrackGraph();
  std::vector<double> stacks(static_cast<size_t>(grid_->getNumLayers())
                                 * grid_->getXSize() * grid_->getYSize(),
                             0.0);
  stacks[grid_->edgeFlatIndex(1, 2, 1)] = 6.0;  // M2 cell (2, 1)
  estimator.addStackedVias(stacks);
  const std::vector<double> capacity = estimator.capacities();
  const double left = capacity[grid_->edgeFlatIndex(1, 1, 1)];
  const double right = capacity[grid_->edgeFlatIndex(1, 2, 1)];
  EXPECT_LT(left, 15.0);
  EXPECT_DOUBLE_EQ(left, right);
  EXPECT_DOUBLE_EQ(capacity[grid_->edgeFlatIndex(1, 0, 1)], 15.0);
  EXPECT_DOUBLE_EQ(capacity[grid_->edgeFlatIndex(1, 3, 1)], 15.0);
  EXPECT_DOUBLE_EQ(capacity[grid_->edgeFlatIndex(1, 2, 0)], 15.0);
  EXPECT_GT(estimator.totals().stacked_vias[1], 0.0);
}

}  // namespace
}  // namespace grt
