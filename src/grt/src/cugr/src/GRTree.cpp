#include "GRTree.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

#include "utl/Logger.h"

namespace grt {

void GRTreeNode::preorder(
    const std::shared_ptr<GRTreeNode>& node,
    const std::function<void(const std::shared_ptr<GRTreeNode>&)>& visit)
{
  if (!node) {
    return;
  }
  visit(node);
  for (auto& child : node->children_) {
    preorder(child, visit);
  }
}

void GRTreeNode::print(const std::shared_ptr<GRTreeNode>& node,
                       utl::Logger* logger)
{
  preorder(node, [logger](const std::shared_ptr<GRTreeNode>& node) {
    std::stringstream ss;
    ss << *node << (!node->children_.empty() ? " -> " : "");
    for (auto& child : node->children_) {
      ss << *child << (child == node->children_.back() ? "" : ", ");
    }
    logger->report("{}", ss.str());
  });
}

std::vector<GRPoint> pinLinkPath(const GRPoint& from,
                                 const GRPoint& to,
                                 const int horizontal_layer,
                                 const int vertical_layer)
{
  std::vector<GRPoint> path{from};
  auto add = [&](const GRPoint& point) {
    const GRPoint& last = path.back();
    if (last.getLayerIdx() != point.getLayerIdx() || last.x() != point.x()
        || last.y() != point.y()) {
      path.push_back(point);
    }
  };
  int x = from.x();
  const int y = from.y();
  if (from.x() != to.x()) {
    add(GRPoint(horizontal_layer, x, y));
    x = to.x();
    add(GRPoint(horizontal_layer, x, y));
  }
  if (from.y() != to.y()) {
    add(GRPoint(vertical_layer, x, y));
    add(GRPoint(vertical_layer, x, to.y()));
  }
  add(to);
  return path;
}

}  // namespace grt
