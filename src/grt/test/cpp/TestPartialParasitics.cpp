// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Planar route keys: FastRoute re-estimates the parasitics of a net between
// slack updates exactly when its key changed, so equal keys must mean equal
// inputs of getPlanarRoute() and every change of those inputs must change the
// key.

#include <cstddef>
#include <vector>

#include "FastRoute.h"
#include "gtest/gtest.h"

namespace grt {
namespace {

TreeEdge makeEdge(const int len, const std::vector<GPoint3D>& grids)
{
  TreeEdge edge;
  edge.len = len;
  edge.route.type = RouteType::MazeRoute;
  edge.route.grids = grids;
  edge.route.routelen = static_cast<int>(grids.size()) - 1;
  return edge;
}

StTree makeTree()
{
  StTree tree;
  tree.edges.push_back(makeEdge(2, {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}));
  tree.edges.push_back(makeEdge(0, {{2, 0, 0}}));
  tree.edges.push_back(makeEdge(1, {{2, 0, 0}, {2, 1, 0}}));
  return tree;
}

std::vector<int> keyOf(const StTree& tree, const int min_layer)
{
  std::vector<int> values;
  appendPlanarRouteKey(tree, min_layer, values);
  return values;
}

TEST(PlanarRouteKey, EqualTreesHaveEqualKeys)
{
  EXPECT_EQ(keyOf(makeTree(), 1), keyOf(makeTree(), 1));
}

TEST(PlanarRouteKey, MinLayerIsPartOfTheKey)
{
  EXPECT_NE(keyOf(makeTree(), 1), keyOf(makeTree(), 2));
}

TEST(PlanarRouteKey, GridPointChangesTheKey)
{
  StTree tree = makeTree();
  tree.edges[0].route.grids[1] = {1, 1, 0};
  EXPECT_NE(keyOf(makeTree(), 1), keyOf(tree, 1));
}

TEST(PlanarRouteKey, LayerOfAGridPointIsNotPartOfTheKey)
{
  // getPlanarRoute() derives the layers from the net and the directions.
  StTree tree = makeTree();
  tree.edges[0].route.grids[1].layer = 3;
  EXPECT_EQ(keyOf(makeTree(), 1), keyOf(tree, 1));
}

TEST(PlanarRouteKey, EdgeWithoutLengthIsSkipped)
{
  StTree tree = makeTree();
  tree.edges[1].route.grids[0] = {5, 5, 0};
  EXPECT_EQ(keyOf(makeTree(), 1), keyOf(tree, 1));
}

TEST(PlanarRouteKey, EdgeGainingLengthChangesTheKey)
{
  StTree tree = makeTree();
  tree.edges[1] = makeEdge(1, {{2, 0, 0}, {3, 0, 0}});
  EXPECT_NE(keyOf(makeTree(), 1), keyOf(tree, 1));
}

TEST(PlanarRouteKey, SplittingAPathBetweenEdgesChangesTheKey)
{
  // The same grid points split differently between two edges give different
  // segments and bridges in getPlanarRoute().
  StTree a;
  a.edges.push_back(makeEdge(1, {{0, 0, 0}, {1, 0, 0}}));
  a.edges.push_back(makeEdge(2, {{1, 0, 0}, {2, 0, 0}, {3, 0, 0}}));
  StTree b;
  b.edges.push_back(makeEdge(2, {{0, 0, 0}, {1, 0, 0}, {1, 0, 0}}));
  b.edges.push_back(makeEdge(1, {{2, 0, 0}, {3, 0, 0}}));
  EXPECT_NE(keyOf(a, 1), keyOf(b, 1));
}

TEST(PlanarRouteKey, RouteLengthBelowTheGridsIsPartOfTheKey)
{
  // Only grids[0 .. routelen] are route points.
  StTree tree = makeTree();
  tree.edges[0].route.routelen = 1;
  EXPECT_NE(keyOf(makeTree(), 1), keyOf(tree, 1));
}

PlanarRouteKeys keysOf(const std::vector<StTree>& trees)
{
  PlanarRouteKeys keys;
  for (const StTree& tree : trees) {
    keys.begin.push_back(keys.values.size());
    appendPlanarRouteKey(tree, 1, keys.values);
  }
  keys.begin.push_back(keys.values.size());
  return keys;
}

TEST(PlanarRouteKeys, SameKeyComparesOneNetAtItsPosition)
{
  StTree changed = makeTree();
  changed.edges[2].route.grids[1] = {3, 0, 0};
  StTree longer = makeTree();
  longer.edges.push_back(makeEdge(1, {{0, 0, 0}, {0, 1, 0}}));

  const PlanarRouteKeys before = keysOf({makeTree(), makeTree(), makeTree()});
  const PlanarRouteKeys after = keysOf({longer, changed, makeTree()});
  ASSERT_EQ(before.netCount(), 3u);
  ASSERT_EQ(after.netCount(), 3u);
  EXPECT_FALSE(after.sameKey(before, 0));
  EXPECT_FALSE(after.sameKey(before, 1));
  // Net 2's key is equal although the keys before it are longer.
  EXPECT_TRUE(after.sameKey(before, 2));
  EXPECT_EQ(after.keySize(2), before.keySize(2));
}

TEST(PlanarRouteKeys, ClearEmptiesTheKeys)
{
  PlanarRouteKeys keys = keysOf({makeTree()});
  keys.clear();
  EXPECT_EQ(keys.netCount(), 0u);
  EXPECT_TRUE(keys.values.empty());
}

}  // namespace
}  // namespace grt
