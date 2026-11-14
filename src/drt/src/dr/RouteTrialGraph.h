// SPDX-License-Identifier: BSD-3-Clause
// Bounded experimental track graph. Edges are geometry-filtered, NOT DRC
// proofs.
#pragma once

#include "dr/RouteTrial.h"

namespace drt::route_trial {

struct SearchGraph
{
  struct Edge
  {
    int from;
    int to;
    int64_t length;
    int64_t wrong_direction;
    std::string via;
  };
  std::vector<Node> nodes;
  std::map<Node, int> ids;
  std::vector<Edge> edges;
};

inline std::vector<frLayerNum> searchLayers(frDesign* design,
                                            const Plan& plan,
                                            RouterConfiguration* cfg)
{
  std::vector<frLayerNum> result;
  const auto first = std::get<0>(*plan.anchors.begin());
  const auto last = std::get<0>(*plan.anchors.rbegin());
  for (const auto& layer : design->getTech()->getLayers()) {
    const auto number = layer->getLayerNum();
    if (!layer->isRoutable() || layer->isBackside()
        || number < cfg->BOTTOM_ROUTING_LAYER
        || number > cfg->TOP_ROUTING_LAYER) {
      continue;
    }
    // Routing layers alternate with cut layers in the supported graph.
    const int minimum_vias
        = plan.tree ? (std::max(number, last) - std::min(number, first)) / 2
                    : (std::abs(number - first) + std::abs(number - last)) / 2;
    if (minimum_vias < static_cast<int>(plan.vias.size())) {
      result.push_back(number);
    }
  }
  return result;
}

inline std::pair<std::vector<frCoord>, std::vector<frCoord>> searchCoordinates(
    frDesign* design,
    const Plan& plan,
    const std::vector<frLayerNum>& layers)
{
  std::set<frCoord> xs, ys;
  odb::Rect box;
  box.mergeInit();
  for (const auto* segment : plan.segments) {
    for (const auto& point :
         {segment->getBeginPoint(), segment->getEndPoint()}) {
      xs.insert(point.x());
      ys.insert(point.y());
      box.merge(odb::Rect(point, point));
    }
  }
  for (const auto* via : plan.vias) {
    xs.insert(via->getOrigin().x());
    ys.insert(via->getOrigin().y());
    box.merge(odb::Rect(via->getOrigin(), via->getOrigin()));
  }
  for (const auto& [layer, x, y] : plan.anchors) {
    xs.insert(x);
    ys.insert(y);
  }
  frCoord halo = 0;
  for (const auto layer : layers) {
    halo = std::max(halo,
                    2
                        * static_cast<frCoord>(
                            design->getTech()->getLayer(layer)->getPitch()));
  }
  box = box.bloat(halo);
  const auto die = design->getTopBlock()->getDieBox();
  for (const auto layer : layers) {
    for (const auto& track : design->getTopBlock()->getTrackPatterns(layer)) {
      const bool x_axis = track->isHorizontal();
      const auto low = x_axis ? std::max(box.xMin(), die.xMin())
                              : std::max(box.yMin(), die.yMin());
      const auto high = x_axis ? std::min(box.xMax(), die.xMax())
                               : std::min(box.yMax(), die.yMax());
      const int64_t step = track->getTrackSpacing();
      if (step <= 0) {
        continue;
      }
      const int64_t start = track->getStartCoord();
      const int64_t first
          = std::max<int64_t>(0, (int64_t(low) - start + step - 1) / step);
      auto& coords = x_axis ? xs : ys;
      for (int64_t index = first; index < track->getNumTracks(); ++index) {
        const int64_t value = start + index * step;
        if (value > high) {
          break;
        }
        coords.insert(value);
      }
    }
  }
  return {{xs.begin(), xs.end()}, {ys.begin(), ys.end()}};
}

inline bool graphSegment(frDesign* design,
                         const Plan& plan,
                         frLayerNum layer,
                         const odb::Point& a,
                         const odb::Point& b)
{
  const auto* tech_layer = design->getTech()->getLayer(layer);
  frPathSeg segment;
  segment.setLayerNum(layer);
  segment.setPoints(a, b);
  segment.setStyle(tech_layer->getDefaultSegStyle());
  if (tech_layer->isUnidirectional()
      && tech_layer->isHorizontal() == segment.isVertical()) {
    return false;
  }
  return via_cleanup::onTrack(design->getTopBlock(), &segment, layer)
         && contactsAtAnchors(design->getRegionQuery(),
                              plan,
                              layer,
                              segment.getBBox(),
                              odb::Rect(a, b));
}

inline bool graphVia(frDesign* design,
                     const Plan& plan,
                     const frViaDef* def,
                     const odb::Point& origin)
{
  if (def->isMultiCut() || !viaOnTracks(design, def, origin)) {
    return false;
  }
  frVia via(def, origin);
  const odb::Rect point(origin, origin);
  for (const auto layer : {def->getLayer1Num(), def->getLayer2Num()}) {
    if (!contactsAtAnchors(design->getRegionQuery(),
                           plan,
                           layer,
                           via_cleanup::metalBox(&via, layer),
                           point)) {
      return false;
    }
  }
  frRegionQuery::Objects<frBlockObject> fixed;
  design->getRegionQuery()->query(
      via.getCutBBox(), def->getCutLayerNum(), fixed);
  if (!fixed.empty()) {
    return false;
  }
  std::vector<frBlockObject*> routed;
  design->getRegionQuery()->queryDRObj(via.getCutBBox(), routed);
  for (auto* object : routed) {
    if (object->typeId() != frcVia) {
      continue;
    }
    auto* other = static_cast<frVia*>(object);
    if (std::find(plan.vias.begin(), plan.vias.end(), other) == plan.vias.end()
        && other->getViaDef()->getCutLayerNum() == def->getCutLayerNum()
        && other->getCutBBox().intersects(via.getCutBBox())) {
      return false;
    }
  }
  return true;
}

inline void addGraphEdges(frDesign* design,
                          const Plan& plan,
                          const std::vector<frLayerNum>& layers,
                          const std::vector<frCoord>& xs,
                          const std::vector<frCoord>& ys,
                          SearchGraph& graph)
{
  for (const auto layer : layers) {
    const auto* tech_layer = design->getTech()->getLayer(layer);
    for (size_t x = 0; x < xs.size(); ++x) {
      for (size_t y = 0; y < ys.size(); ++y) {
        const int from = graph.ids.at({layer, xs[x], ys[y]});
        for (const bool horizontal : {true, false}) {
          if ((horizontal && x + 1 == xs.size())
              || (!horizontal && y + 1 == ys.size())) {
            continue;
          }
          const odb::Point a(xs[x], ys[y]);
          const odb::Point b(xs[x + horizontal], ys[y + !horizontal]);
          if (graphSegment(design, plan, layer, a, b)) {
            const int64_t length = odb::Point::manhattanDistance(a, b);
            graph.edges.push_back(
                {from,
                 graph.ids.at({layer, b.x(), b.y()}),
                 length,
                 tech_layer->isHorizontal() == horizontal ? 0 : length,
                 "-"});
          }
        }
        if (std::find(layers.begin(), layers.end(), layer + 2)
            == layers.end()) {
          continue;
        }
        const auto* def
            = design->getTech()->getLayer(layer + 1)->getDefaultViaDef();
        if (def && def->getLayer1Num() == layer
            && def->getLayer2Num() == layer + 2
            && graphVia(design, plan, def, {xs[x], ys[y]})) {
          graph.edges.push_back({from,
                                 graph.ids.at({layer + 2, xs[x], ys[y]}),
                                 0,
                                 0,
                                 def->getName()});
        }
      }
    }
  }
}

inline bool exportGraph(frDesign* design,
                        RouterConfiguration* cfg,
                        const Plan& plan,
                        const std::string& filename,
                        utl::Logger* logger)
{
  if (!sourceIsolated(design->getRegionQuery(), plan)) {
    logger->report("ROUTE_GRAPH rejected stage=source_contact");
    return false;
  }
  const auto layers = searchLayers(design, plan, cfg);
  const auto [xs, ys] = searchCoordinates(design, plan, layers);
  if (layers.empty() || xs.empty() || ys.empty()
      || layers.size() * xs.size() * ys.size() > 40000) {
    logger->report("ROUTE_GRAPH rejected stage=node_budget");
    return false;
  }
  SearchGraph graph;
  for (const auto layer : layers) {
    for (const auto x : xs) {
      for (const auto y : ys) {
        const Node node{layer, x, y};
        graph.ids.emplace(node, graph.nodes.size());
        graph.nodes.push_back(node);
      }
    }
  }
  for (const auto& anchor : plan.anchors) {
    if (!graph.ids.count(anchor)) {
      logger->report("ROUTE_GRAPH rejected stage=anchor_layer");
      return false;
    }
  }
  addGraphEdges(design, plan, layers, xs, ys, graph);
  std::ofstream file(filename);
  if (plan.tree) {
    file << "tree_graph_v1 " << plan.old_length << ' ' << plan.vias.size()
         << ' ' << plan.anchors.size();
    for (const auto& port : plan.anchors) {
      file << ' ' << graph.ids.at(port);
    }
    file << '\n';
  } else {
    file << "graph_v1 " << plan.old_length << ' ' << plan.vias.size() << ' '
         << graph.ids.at(*plan.anchors.begin()) << ' '
         << graph.ids.at(*plan.anchors.rbegin()) << '\n';
  }
  for (size_t id = 0; id < graph.nodes.size(); ++id) {
    const auto [layer, x, y] = graph.nodes[id];
    file << "node " << id << ' '
         << design->getTech()->getLayer(layer)->getName() << ' ' << layer << ' '
         << x << ' ' << y << '\n';
  }
  for (const auto& edge : graph.edges) {
    file << "edge " << edge.from << ' ' << edge.to << ' ' << edge.length << ' '
         << edge.wrong_direction << ' ' << edge.via << '\n';
  }
  file.close();
  if (!file) {
    throw std::runtime_error("Could not write route graph");
  }
  logger->report("ROUTE_GRAPH nodes={} edges={} geometry_only=1",
                 graph.nodes.size(),
                 graph.edges.size());
  return true;
}

}  // namespace drt::route_trial
