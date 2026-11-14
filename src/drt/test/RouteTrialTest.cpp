// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include <set>
#include <tuple>

#include "dr/RouteTrial.h"
#include "dr/RouteTrialGraph.h"
#include "fixture.h"
#include "gtest/gtest.h"
#include "io/io.h"
#include "odb/dbWireCodec.h"

namespace drt {
namespace {

class RouteTrialTest : public Fixture
{
 protected:
  RouteTrialTest()
  {
    addLayer(design->getTech(), "v1", odb::dbTechLayerType::CUT);
    addLayer(design->getTech(),
             "m2",
             odb::dbTechLayerType::ROUTING,
             odb::dbTechLayerDir::HORIZONTAL);
    for (const auto layer : {2, 4}) {
      for (const bool x_axis : {true, false}) {
        design->getTopBlock()->addTrackPattern(
            std::make_unique<frTrackPattern>(x_axis, 0, 30, 100, layer));
      }
    }
    net_ = makeNet("signal");
    makePathseg(net_, 2, {400, 500}, {400, 900});
    makePathseg(net_, 4, {400, 900}, {1000, 900});
    makePathseg(net_, 2, {1000, 500}, {1000, 900});
    auto* def = makeViaDef("via", 3, {-50, -50}, {50, 50});
    makeVia(def, net_, {400, 900});
    makeVia(def, net_, {1000, 900});
    net_->setModified(false);
  }

  route_trial::Plan plan()
  {
    route_trial::Plan result{net_};
    result.anchors = {{2, 400, 500}, {2, 1000, 500}};
    for (const auto& shape : net_->getShapes()) {
      result.segments.push_back(static_cast<frPathSeg*>(shape.get()));
    }
    for (const auto& via : net_->getVias()) {
      result.vias.push_back(via.get());
    }
    auto replacement = std::make_unique<frPathSeg>();
    replacement->setLayerNum(2);
    replacement->setPoints({400, 500}, {1000, 500});
    replacement->setStyle(design->getTech()->getLayer(2)->getDefaultSegStyle());
    replacement->addToNet(net_);
    result.additions.push_back(std::move(replacement));
    result.old_length = 1400;
    result.new_length = 600;
    return result;
  }

  using Geometry
      = std::set<std::tuple<int, int, int, int, int, frBlockObject*>>;
  Geometry geometry()
  {
    Geometry result;
    for (int layer = 0; layer <= 4; ++layer) {
      frRegionQuery::Objects<frBlockObject> objects;
      design->getRegionQuery()->queryDRObj(
          odb::Rect(0, 0, 2000, 2000), layer, objects);
      for (const auto& [box, object] : objects) {
        result.emplace(
            layer, box.xMin(), box.yMin(), box.xMax(), box.yMax(), object);
      }
    }
    return result;
  }
  frNet* net_;
};

TEST_F(RouteTrialTest, AcceptsShorterMultiLayerChain)
{
  auto candidate = plan();
  initRegionQuery();
  ASSERT_TRUE(route_trial::apply(
      design.get(), logger.get(), router_cfg.get(), candidate));
  EXPECT_TRUE(net_->getVias().empty());
  ASSERT_EQ(net_->getShapes().size(), 1);
  auto* segment = static_cast<frPathSeg*>(net_->getShapes().front().get());
  EXPECT_EQ(segment->getLayerNum(), 2);
  EXPECT_EQ(segment->getBeginPoint(), odb::Point(400, 500));
  EXPECT_EQ(segment->getEndPoint(), odb::Point(1000, 500));
  EXPECT_TRUE(net_->isModified());
}

TEST_F(RouteTrialTest, DrcRejectionRestoresOriginalOwnersAndGeometry)
{
  auto candidate = plan();
  makeSpacingConstraint(2);
  auto* other = makeNet("other");
  makePathseg(other, 2, {600, 650}, {800, 650});
  initRegionQuery();
  ASSERT_TRUE(route_trial::isolated(design->getRegionQuery(), candidate));
  const auto before = geometry();
  EXPECT_FALSE(route_trial::apply(
      design.get(), logger.get(), router_cfg.get(), candidate));
  EXPECT_EQ(geometry(), before);
  EXPECT_EQ(net_->getShapes().size(), 3);
  EXPECT_EQ(net_->getVias().size(), 2);
  EXPECT_FALSE(net_->isModified());
}

TEST_F(RouteTrialTest, RejectsInteriorStackedVia)
{
  auto candidate = plan();
  auto* def = makeViaDef("lower", 1, {-50, -50}, {50, 50});
  makeVia(def, net_, {400, 700});
  initRegionQuery();
  EXPECT_FALSE(route_trial::isolated(design->getRegionQuery(), candidate));
}

TEST_F(RouteTrialTest, AllowsRetainedEndpointConnection)
{
  auto candidate = plan();
  makePathseg(net_, 2, {400, 100}, {400, 500});
  makePathseg(net_, 2, {1000, 100}, {1000, 500});
  initRegionQuery();
  EXPECT_TRUE(route_trial::isolated(design->getRegionQuery(), candidate));
}

TEST_F(RouteTrialTest, RejectsOverlapBeyondRetainedEndpoint)
{
  auto candidate = plan();
  makePathseg(net_, 2, {400, 100}, {400, 700});
  initRegionQuery();
  EXPECT_FALSE(route_trial::isolated(design->getRegionQuery(), candidate));
}

TEST_F(RouteTrialTest, RejectsNewInteriorJunction)
{
  auto candidate = plan();
  makePathseg(net_, 2, {700, 100}, {700, 500});
  initRegionQuery();
  EXPECT_FALSE(route_trial::isolated(design->getRegionQuery(), candidate));
}

TEST_F(RouteTrialTest, RejectsFixedContactBeforeMutation)
{
  auto candidate = plan();
  const auto [master, db_master] = makeMacro("fixed", 0, 0, 2000, 2000);
  makeMacroObs(master, 650, 450, 750, 550, 2);
  makeInst("fixed_instance", master, db_master);
  initRegionQuery();
  const auto before = geometry();
  EXPECT_FALSE(route_trial::apply(
      design.get(), logger.get(), router_cfg.get(), candidate));
  EXPECT_EQ(geometry(), before);
}

TEST_F(RouteTrialTest, ResolvesAndAppliesFourViasToTwo)
{
  auto* net = makeNet("four_vias");
  auto* def = design->getTech()->getVia("via");
  makePathseg(net, 2, {1400, 1500}, {1400, 1900});
  makePathseg(net, 4, {1400, 1900}, {1600, 1900});
  makePathseg(net, 2, {1600, 1900}, {1800, 1900});
  makePathseg(net, 4, {1800, 1900}, {2000, 1900});
  makePathseg(net, 2, {2000, 1500}, {2000, 1900});
  for (const auto x : {1400, 1600, 1800, 2000}) {
    makeVia(def, net, {x, 1900});
  }
  route_trial::Recipe recipe;
  recipe.net = "four_vias";
  recipe.anchors = {{"m1", 1400, 1500}, {"m1", 2000, 1500}};
  for (const auto& shape : net->getShapes()) {
    auto* segment = static_cast<frPathSeg*>(shape.get());
    recipe.remove_segments.push_back(
        {design->getTech()->getLayer(segment->getLayerNum())->getName(),
         segment->getBeginPoint(),
         segment->getEndPoint()});
  }
  for (const auto& via : net->getVias()) {
    recipe.remove_vias.push_back({def->getName(), via->getOrigin()});
  }
  recipe.add_segments = {{"m2", {1400, 1500}, {2000, 1500}}};
  recipe.add_vias
      = {{def->getName(), {1400, 1500}}, {def->getName(), {2000, 1500}}};
  auto candidate = route_trial::resolve(design.get(), recipe);
  EXPECT_EQ(candidate.old_length, 1400);
  EXPECT_EQ(candidate.new_length, 600);
  initRegionQuery();
  ASSERT_TRUE(route_trial::apply(
      design.get(), logger.get(), router_cfg.get(), candidate));
  EXPECT_EQ(net->getShapes().size(), 1);
  EXPECT_EQ(net->getVias().size(), 2);
  EXPECT_TRUE(net->isModified());
}

TEST_F(RouteTrialTest, ViaRequiresTrackOnBothMetalLayers)
{
  const auto* def = net_->getVias().front()->getViaDef();
  EXPECT_TRUE(route_trial::viaOnTracks(design.get(), def, {700, 500}));
  EXPECT_FALSE(route_trial::viaOnTracks(design.get(), def, {700, 511}));
}

TEST_F(RouteTrialTest, RemovesViasWithoutIncreasingLength)
{
  route_trial::Recipe recipe;
  recipe.net = "signal";
  recipe.anchors = {{"m1", 400, 900}, {"m1", 1000, 900}};
  recipe.remove_segments = {{"m2", {400, 900}, {1000, 900}}};
  recipe.remove_vias = {{"via", {400, 900}}, {"via", {1000, 900}}};
  recipe.add_segments = {{"m1", {400, 900}, {1000, 900}}};
  auto candidate = route_trial::resolve(design.get(), recipe);
  EXPECT_EQ(candidate.old_length, 600);
  EXPECT_EQ(candidate.new_length, 600);
  initRegionQuery();
  ASSERT_TRUE(route_trial::apply(
      design.get(), logger.get(), router_cfg.get(), candidate));
  EXPECT_TRUE(net_->getVias().empty());
  EXPECT_EQ(net_->getShapes().size(), 3);
}

TEST_F(RouteTrialTest, RejectsLongerRouteBeforeMutation)
{
  route_trial::Recipe recipe;
  recipe.net = "signal";
  recipe.anchors = {{"m1", 400, 900}, {"m1", 1000, 900}};
  recipe.remove_segments = {{"m2", {400, 900}, {1000, 900}}};
  recipe.remove_vias = {{"via", {400, 900}}, {"via", {1000, 900}}};
  recipe.add_segments = {{"m1", {400, 900}, {400, 1200}},
                         {"m1", {400, 1200}, {1000, 1200}},
                         {"m1", {1000, 900}, {1000, 1200}}};
  initRegionQuery();
  const auto before = geometry();
  EXPECT_THROW(route_trial::resolve(design.get(), recipe), std::runtime_error);
  EXPECT_EQ(geometry(), before);
  EXPECT_EQ(net_->getVias().size(), 2);
  EXPECT_FALSE(net_->isModified());
}

TEST_F(RouteTrialTest, ReconnectsThreePortsAndRetainsExternalBranch)
{
  makePathseg(net_, 4, {700, 900}, {700, 1300});
  route_trial::Recipe recipe;
  recipe.net = "signal";
  recipe.tree = true;
  recipe.anchors = {{"m1", 400, 500}, {"m1", 1000, 500}, {"m2", 700, 1300}};
  for (const auto& shape : net_->getShapes()) {
    const auto* segment = static_cast<frPathSeg*>(shape.get());
    recipe.remove_segments.push_back(
        {design->getTech()->getLayer(segment->getLayerNum())->getName(),
         segment->getBeginPoint(),
         segment->getEndPoint()});
  }
  recipe.remove_vias = {{"via", {400, 900}}, {"via", {1000, 900}}};
  recipe.add_segments
      = {{"m1", {400, 500}, {1000, 500}}, {"m1", {700, 500}, {700, 1300}}};
  recipe.add_vias = {{"via", {700, 1300}}};
  // Whole wire objects meet at an interior T junction in both source
  // and replacement. The third retained port must remain connected.
  makePathseg(net_, 4, {700, 1300}, {900, 1300});
  auto candidate = route_trial::resolve(design.get(), recipe);
  EXPECT_EQ(candidate.old_length, 1800);
  EXPECT_EQ(candidate.new_length, 1400);
  initRegionQuery();
  ASSERT_TRUE(route_trial::apply(
      design.get(), logger.get(), router_cfg.get(), candidate));
  EXPECT_EQ(net_->getVias().size(), 1);
  EXPECT_EQ(net_->getShapes().size(), 3);
  EXPECT_TRUE(net_->isModified());
}

TEST_F(RouteTrialTest, RejectsTreeThatLosesOnePortBeforeMutation)
{
  route_trial::Recipe recipe;
  recipe.net = "signal";
  recipe.tree = true;
  makePathseg(net_, 4, {700, 900}, {700, 1300});
  recipe.anchors = {{"m1", 400, 500}, {"m1", 1000, 500}, {"m2", 700, 1300}};
  for (const auto& shape : net_->getShapes()) {
    const auto* segment = static_cast<frPathSeg*>(shape.get());
    recipe.remove_segments.push_back(
        {design->getTech()->getLayer(segment->getLayerNum())->getName(),
         segment->getBeginPoint(),
         segment->getEndPoint()});
  }
  recipe.remove_vias = {{"via", {400, 900}}, {"via", {1000, 900}}};
  recipe.add_segments = {{"m1", {400, 500}, {1000, 500}}};
  initRegionQuery();
  const auto before = geometry();
  EXPECT_THROW(route_trial::resolve(design.get(), recipe), std::runtime_error);
  EXPECT_EQ(geometry(), before);
  EXPECT_FALSE(net_->isModified());
}

TEST_F(RouteTrialTest, TemporaryAddedViasRollbackExactly)
{
  auto candidate = plan();
  auto via = std::make_unique<frVia>(net_->getVias().front()->getViaDef(),
                                     odb::Point(700, 500));
  via->addToNet(net_);
  auto* added = via.get();
  candidate.added_vias.push_back(std::move(via));
  initRegionQuery();
  const auto before = geometry();
  {
    route_trial::Trial trial(design->getRegionQuery(), candidate);
    const auto during = geometry();
    EXPECT_NE(during, before);
    EXPECT_TRUE(
        std::any_of(during.begin(), during.end(), [added](const auto& item) {
          return std::get<5>(item) == added;
        }));
  }
  EXPECT_EQ(geometry(), before);
  EXPECT_EQ(net_->getVias().size(), 2);
  EXPECT_FALSE(net_->isModified());
}

TEST_F(RouteTrialTest, AddedViasCommitTransfersOwnership)
{
  auto candidate = plan();
  auto via = std::make_unique<frVia>(net_->getVias().front()->getViaDef(),
                                     odb::Point(700, 500));
  via->addToNet(net_);
  auto* added = via.get();
  candidate.added_vias.push_back(std::move(via));
  initRegionQuery();
  {
    route_trial::Trial trial(design->getRegionQuery(), candidate);
    trial.commit();
  }
  ASSERT_EQ(net_->getVias().size(), 1);
  EXPECT_EQ(net_->getVias().front().get(), added);
  const auto after = geometry();
  EXPECT_TRUE(
      std::any_of(after.begin(), after.end(), [added](const auto& item) {
        return std::get<5>(item) == added;
      }));
}

TEST_F(RouteTrialTest, GraphFiltersObstacleAndKeepsLegalDetourEdges)
{
  auto candidate = plan();
  auto* obstacle = makeNet("obstacle");
  makePathseg(obstacle, 2, {700, 400}, {700, 600});
  initRegionQuery();
  const auto before = geometry();
  EXPECT_FALSE(route_trial::graphSegment(
      design.get(), candidate, 2, {400, 500}, {1000, 500}));
  EXPECT_TRUE(route_trial::graphSegment(
      design.get(), candidate, 2, {400, 700}, {1000, 700}));
  EXPECT_FALSE(route_trial::graphSegment(
      design.get(), candidate, 2, {400, 711}, {1000, 711}));
  EXPECT_EQ(geometry(), before);
}

TEST_F(RouteTrialTest, GraphIncludesHaloTracksAndClipsToDie)
{
  frBoundary boundary;
  boundary.setPoints({{200, 200}, {1200, 200}, {1200, 1100}, {200, 1100}});
  design->getTopBlock()->setBoundaries({boundary});
  const auto candidate = plan();
  const auto [xs, ys]
      = route_trial::searchCoordinates(design.get(), candidate, {2, 4});
  // Source bounds are [400, 1000] x [500, 900]. Two 200-DBU
  // pitches add outside tracks, but nothing beyond the die is exported.
  ASSERT_FALSE(xs.empty());
  ASSERT_FALSE(ys.empty());
  EXPECT_EQ(xs.front(), 200);
  EXPECT_EQ(xs.back(), 1200);
  EXPECT_EQ(ys.front(), 200);
  EXPECT_EQ(ys.back(), 1100);
}

TEST_F(RouteTrialTest, GraphViaRejectsMetalObstacle)
{
  auto candidate = plan();
  const auto* def = net_->getVias().front()->getViaDef();
  auto* obstacle = makeNet("obstacle");
  makePathseg(obstacle, 4, {600, 500}, {800, 500});
  initRegionQuery();
  EXPECT_FALSE(route_trial::graphVia(design.get(), candidate, def, {700, 500}));
  EXPECT_TRUE(
      route_trial::graphVia(design.get(), candidate, def, {1100, 1200}));
}

TEST_F(RouteTrialTest, ImportedTaperSurvivesBendsViasAndRuleChanges)
{
  auto* db = db_tech->getDb();
  auto* chip = odb::dbChip::create(db, db_tech);
  auto* block = odb::dbBlock::create(chip, "import_block");
  auto* net = odb::dbNet::create(block, "imported");
  auto* ndr = odb::dbTechNonDefaultRule::create(db_tech, "import_ndr");
  net->setNonDefaultRule(ndr);
  auto* m1 = db_tech->findLayer("m1");
  auto* m2 = db_tech->findLayer("m2");
  auto* rule = odb::dbTechLayerRule::create(ndr, m1);
  rule->setWidth(100);
  rule->setSpacing(200);
  auto* upper_rule = odb::dbTechLayerRule::create(ndr, m2);
  upper_rule->setWidth(100);
  upper_rule->setSpacing(200);
  auto fr_rule = std::make_unique<frNonDefaultRule>();
  fr_rule->setName("import_ndr");
  design->getTech()->addNDR(std::move(fr_rule));
  auto* db_via = odb::dbTechVia::create(db_tech, "import_via");
  for (const auto* name : {"m1", "v1", "m2"}) {
    odb::dbBox::create(db_via, db_tech->findLayer(name), -50, -50, 50, 50);
  }
  makeViaDef("import_via", 3, {-50, -50}, {50, 50});
  odb::dbWireEncoder encoder;
  encoder.begin(odb::dbWire::create(net));
  encoder.newPath(m1, odb::dbWireType::ROUTED);
  encoder.addPoint(100, 100);
  encoder.addPoint(500, 100);
  encoder.addPoint(500, 500);
  encoder.addTechVia(db_via);
  encoder.newPath(m1, odb::dbWireType::ROUTED, rule);
  encoder.addPoint(100, 800);
  encoder.addPoint(500, 800);
  encoder.addTechVia(db_via);
  encoder.newPath(m1, odb::dbWireType::ROUTED);
  encoder.addPoint(100, 1200);
  encoder.addPoint(500, 1200);
  encoder.end();
  io::Parser parser(db, design.get(), logger.get(), router_cfg.get());
  const auto* imported = parser.addNet(net);
  std::vector<bool> segment_flags, via_flags;
  for (const auto& shape : imported->getShapes()) {
    segment_flags.push_back(static_cast<frPathSeg*>(shape.get())->isTapered());
  }
  for (const auto& via : imported->getVias()) {
    via_flags.push_back(via->isTapered());
  }
  EXPECT_EQ(segment_flags, (std::vector<bool>{true, true, false, true}));
  EXPECT_EQ(via_flags, (std::vector<bool>{true, false}));
}

TEST(RouteTrialTopology, RequiresOneConnectedChain)
{
  using route_trial::Node;
  const Node a{2, 0, 0}, b{2, 100, 0}, c{2, 200, 0};
  const Node x{4, 0, 100}, y{4, 100, 100}, z{4, 200, 100};
  EXPECT_TRUE(route_trial::isChain({{a, b}, {b, c}}, {a, c}));
  EXPECT_FALSE(route_trial::isChain({{a, b}, {b, c}, {b, x}}, {a, c}));
  EXPECT_FALSE(
      route_trial::isChain({{a, b}, {b, c}, {x, y}, {y, z}, {z, x}}, {a, c}));
  EXPECT_FALSE(route_trial::isChain({{a, b}, {b, c}}, {a, x}));
}

TEST(RouteTrialTopology, TreeSplitsCrossingsButRejectsCyclesAndOverlaps)
{
  using route_trial::Node;
  const Node left{2, 0, 100}, right{2, 200, 100};
  const Node low{2, 100, 0}, high{2, 100, 200}, middle{2, 100, 100};
  EXPECT_TRUE(route_trial::isTree({{left, right}, {low, high}},
                                  {left, right, low, high}));
  EXPECT_TRUE(route_trial::isTree({{left, right}, {middle, high}},
                                  {left, right, middle, high}));
  EXPECT_FALSE(
      route_trial::isTree({{left, right}, {low, high}}, {left, right, high}));
  EXPECT_FALSE(
      route_trial::isTree({{left, right}, {left, middle}}, {left, right}));
  const Node corner{2, 200, 200};
  EXPECT_FALSE(route_trial::isTree(
      {{left, right}, {middle, high}, {high, corner}, {corner, right}},
      {left, right}));
  const Node a{4, 0, 0}, b{4, 0, 100}, c{4, 100, 100}, d{4, 100, 0};
  EXPECT_FALSE(route_trial::isTree(
      {{left, right}, {a, b}, {b, c}, {c, d}, {d, a}}, {left, right}));
}

}  // namespace
}  // namespace drt
