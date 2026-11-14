#pragma once
#include <functional>
#include <memory>
#include <ostream>
#include <utility>
#include <vector>

#include "geo.h"
#include "utl/Logger.h"

namespace grt {

class GRPoint : public PointT
{
 public:
  GRPoint(int l, int x, int y) : PointT(x, y), layer_idx_(l) {}

  int getLayerIdx() const { return layer_idx_; }

  friend std::ostream& operator<<(std::ostream& os, const GRPoint& pt)
  {
    os << "(" << pt.layer_idx_ << ", " << pt.x() << ", " << pt.y() << ")";
    return os;
  }

 private:
  int layer_idx_;
};

class GRTreeNode : public GRPoint
{
 public:
  GRTreeNode(int l, int x, int y) : GRPoint(l, x, y) {}
  GRTreeNode(const GRPoint& point) : GRPoint(point) {}

  const std::vector<std::shared_ptr<GRTreeNode>>& getChildren() const
  {
    return children_;
  }

  // A nonnegative pin_index makes the edge to `child` a pin link: the two
  // cells are connected by the existing metal of that pin (the zero-cost
  // clique K(V_p) of Gester et al., TODAES 2013, Sec. 2.1/2.2). A pin link
  // is a topological edge only; it carries no wire or via demand.
  void addChild(std::shared_ptr<GRTreeNode> child, int pin_index = -1)
  {
    child->incoming_pin_ = pin_index;
    children_.push_back(std::move(child));
  }

  // True if the edge from this node's parent is a pin link.
  bool isPinLink() const { return incoming_pin_ >= 0; }
  // Pin index (GRNet numbering) of the incoming pin link, -1 otherwise.
  int getLinkPin() const { return incoming_pin_; }
  // True if the edge from this node's parent belongs to the wire of an
  // antenna jumper. On the jumper's layer the pads of the via stacks at its
  // two ends lie on the jumper's own track: they charge no edge the jumper
  // wire runs over (GridGraph::forEachTreeDemandImpl).
  bool isJumperWire() const { return incoming_jumper_; }
  void setJumperWire(const bool jumper) { incoming_jumper_ = jumper; }

  void removeChild(const std::shared_ptr<GRTreeNode>& child)
  {
    for (auto it = children_.begin(); it != children_.end(); ++it) {
      if (*it == child) {
        children_.erase(it);
        break;
      }
    }
  }

  static void preorder(
      const std::shared_ptr<GRTreeNode>& node,
      const std::function<void(const std::shared_ptr<GRTreeNode>&)>& visit);
  static void print(const std::shared_ptr<GRTreeNode>& node,
                    utl::Logger* logger);

 private:
  int incoming_pin_ = -1;
  bool incoming_jumper_ = false;
  std::vector<std::shared_ptr<GRTreeNode>> children_;
};

// Corner points of the bridge that represents a pin link between `from` and
// `to` in a segment-based route export: the x run on `horizontal_layer`, the
// y run on `vertical_layer` (routable layers of those preferred
// directions), joined to the two pin cells by via stacks. Consecutive points
// differ in exactly one coordinate. The bridge only keeps exported routes
// connected; the pin's own metal provides the connection, so it is never
// charged as routing demand.
std::vector<GRPoint> pinLinkPath(const GRPoint& from,
                                 const GRPoint& to,
                                 int horizontal_layer,
                                 int vertical_layer);

}  // namespace grt
