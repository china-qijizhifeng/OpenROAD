// SPDX-License-Identifier: BSD-3-Clause
// Explicit, checkpoint-based route experiments; never used by normal routing.
#pragma once

#include <algorithm>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "dr/ViaCleanup.h"

namespace drt::route_trial {

using Node = std::tuple<frLayerNum, frCoord, frCoord>;

struct Recipe
{
  std::string net;
  bool tree = false;
  struct Segment
  {
    std::string layer;
    odb::Point begin;
    odb::Point end;
  };
  struct Via
  {
    std::string name;
    odb::Point origin;
  };
  std::vector<Segment> remove_segments;
  std::vector<Via> remove_vias;
  std::vector<Via> add_vias;
  std::vector<Segment> add_segments;
  std::vector<std::tuple<std::string, frCoord, frCoord>> anchors;
};

inline Recipe read(const std::string& filename)
{
  std::ifstream file(filename);
  Recipe recipe;
  if (!std::getline(file, recipe.net) || recipe.net.empty()) {
    throw std::runtime_error("Route trial requires a net name");
  }
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) {
      continue;
    }
    std::istringstream stream(line);
    std::string operation, name, extra;
    frCoord x1, y1, x2, y2;
    if (!(stream >> operation >> name)) {
      throw std::runtime_error("Malformed route trial record");
    }
    if (operation == "topology") {
      if (name != "tree" || recipe.tree || (stream >> extra)) {
        throw std::runtime_error("Invalid route trial topology");
      }
      recipe.tree = true;
      continue;
    }
    if (!(stream >> x1 >> y1)) {
      throw std::runtime_error("Malformed route trial record");
    }
    if (operation == "anchor") {
      recipe.anchors.emplace_back(name, x1, y1);
    } else if (operation == "remove_via" || operation == "add_via") {
      auto& vias
          = operation == "remove_via" ? recipe.remove_vias : recipe.add_vias;
      vias.push_back({name, {x1, y1}});
    } else if (operation == "remove_segment" || operation == "add_segment") {
      if (!(stream >> x2 >> y2) || (x1 != x2 && y1 != y2)
          || (x1 == x2 && y1 == y2)) {
        throw std::runtime_error(
            "Route trial requires nonzero Manhattan segments");
      }
      auto& segments = operation == "remove_segment" ? recipe.remove_segments
                                                     : recipe.add_segments;
      odb::Point begin(x1, y1), end(x2, y2);
      if (end < begin) {
        std::swap(begin, end);
      }
      segments.push_back({name, begin, end});
    } else {
      throw std::runtime_error("Unsupported route trial operation");
    }
    if (stream >> extra) {
      throw std::runtime_error("Extra route trial fields");
    }
  }
  return recipe;
}

inline frLayerNum layerNumber(frTechObject* tech, const std::string& name)
{
  auto* layer = tech->getLayer(name);
  if (!layer || !layer->isRoutable() || layer->isBackside()) {
    throw std::runtime_error("Invalid route trial routing layer");
  }
  return layer->getLayerNum();
}

struct Plan
{
  explicit Plan(frNet* net_ptr) : net(net_ptr) {}

  frNet* net;
  bool tree = false;
  std::set<Node> anchors;
  std::vector<frPathSeg*> segments;
  std::vector<frVia*> vias;
  std::vector<std::unique_ptr<frShape>> additions;
  std::vector<std::unique_ptr<frVia>> added_vias;
  int64_t old_length = 0;
  int64_t new_length = 0;
};

inline bool isChain(const std::vector<std::pair<Node, Node>>& edges,
                    const std::set<Node>& anchors)
{
  if (anchors.size() != 2 || edges.empty()) {
    return false;
  }
  std::map<Node, std::vector<Node>> adjacency;
  for (const auto& [a, b] : edges) {
    adjacency[a].push_back(b);
    adjacency[b].push_back(a);
  }
  for (const auto& anchor : anchors) {
    if (adjacency[anchor].size() != 1) {
      return false;
    }
  }
  for (const auto& [node, neighbors] : adjacency) {
    if (neighbors.size() != (anchors.count(node) ? 1 : 2)) {
      return false;
    }
  }
  std::set<Node> visited;
  std::vector<Node> pending{*anchors.begin()};
  while (!pending.empty()) {
    const auto node = pending.back();
    pending.pop_back();
    if (visited.insert(node).second) {
      pending.insert(
          pending.end(), adjacency[node].begin(), adjacency[node].end());
    }
  }
  return visited.size() == adjacency.size();
}

// Multi-port trials preserve one connected tree. Split wires at geometric
// junctions, including retained ports inside whole source objects. Reject
// overlaps and cycles rather than relying on the writer to prune them.
inline bool isTree(const std::vector<std::pair<Node, Node>>& edges,
                   const std::set<Node>& ports)
{
  if (ports.size() < 2 || ports.size() > 8 || edges.empty()) {
    return false;
  }
  std::set<Node> points = ports;
  for (const auto& [a, b] : edges) {
    const auto [la, xa, ya] = a;
    const auto [lb, xb, yb] = b;
    if (a == b || (la == lb && xa != xb && ya != yb)
        || (la != lb && (xa != xb || ya != yb))) {
      return false;
    }
    points.insert(a);
    points.insert(b);
  }
  auto contains = [](const Node& a, const Node& b, const Node& p) {
    return std::get<0>(a) == std::get<0>(p) && std::get<0>(b) == std::get<0>(p)
           && std::get<1>(p) >= std::min(std::get<1>(a), std::get<1>(b))
           && std::get<1>(p) <= std::max(std::get<1>(a), std::get<1>(b))
           && std::get<2>(p) >= std::min(std::get<2>(a), std::get<2>(b))
           && std::get<2>(p) <= std::max(std::get<2>(a), std::get<2>(b));
  };
  for (const auto& [a, b] : edges) {
    if (std::get<0>(a) != std::get<0>(b)) {
      continue;
    }
    for (const auto& [c, d] : edges) {
      if (std::get<0>(a) != std::get<0>(c)
          || std::get<0>(c) != std::get<0>(d)) {
        continue;
      }
      for (const Node p :
           {Node{std::get<0>(a), std::get<1>(a), std::get<2>(c)},
            Node{std::get<0>(a), std::get<1>(c), std::get<2>(a)}}) {
        if (contains(a, b, p) && contains(c, d, p)) {
          points.insert(p);
        }
      }
    }
  }
  std::set<std::pair<Node, Node>> split_edges;
  auto add = [&](Node a, Node b) {
    if (b < a) {
      std::swap(a, b);
    }
    return split_edges.emplace(a, b).second;
  };
  for (const auto& [a, b] : edges) {
    if (std::get<0>(a) != std::get<0>(b)) {
      if (!add(a, b)) {
        return false;
      }
      continue;
    }
    std::vector<Node> on_wire;
    for (const auto& p : points) {
      if (contains(a, b, p)) {
        on_wire.push_back(p);
      }
    }
    for (size_t i = 1; i < on_wire.size(); ++i) {
      if (!add(on_wire[i - 1], on_wire[i])) {
        return false;
      }
    }
  }
  std::map<Node, std::vector<Node>> adjacency;
  for (const auto& [a, b] : split_edges) {
    adjacency[a].push_back(b);
    adjacency[b].push_back(a);
  }
  for (const auto& p : ports) {
    if (!adjacency.count(p)) {
      return false;
    }
  }
  for (const auto& [p, neighbors] : adjacency) {
    if (neighbors.size() == 1 && !ports.count(p)) {
      return false;
    }
  }
  if (split_edges.size() + 1 != adjacency.size()) {
    return false;
  }
  std::set<Node> visited;
  std::vector<Node> pending{*ports.begin()};
  while (!pending.empty()) {
    const auto node = pending.back();
    pending.pop_back();
    if (visited.insert(node).second) {
      pending.insert(
          pending.end(), adjacency[node].begin(), adjacency[node].end());
    }
  }
  return visited.size() == adjacency.size();
}

inline bool validTopology(const std::vector<std::pair<Node, Node>>& edges,
                          const Plan& plan)
{
  return plan.tree ? isTree(edges, plan.anchors) : isChain(edges, plan.anchors);
}

inline void resolveSegments(frTechObject* tech,
                            const Recipe& recipe,
                            Plan& plan)
{
  std::set<frPathSeg*> unique;
  for (const auto& request : recipe.remove_segments) {
    const auto layer = layerNumber(tech, request.layer);
    frPathSeg* match = nullptr;
    for (const auto& shape : plan.net->getShapes()) {
      if (shape->typeId() != frcPathSeg) {
        continue;
      }
      auto* segment = static_cast<frPathSeg*>(shape.get());
      if (segment->getLayerNum() == layer
          && segment->getBeginPoint() == request.begin
          && segment->getEndPoint() == request.end) {
        if (match) {
          throw std::runtime_error("Ambiguous route trial segment");
        }
        match = segment;
      }
    }
    if (!match || !unique.insert(match).second || match->isApPathSeg()
        || match->getStyle().getWidth() != tech->getLayer(layer)->getWidth()) {
      throw std::runtime_error("Unmatched or unsupported route trial segment");
    }
    plan.segments.push_back(match);
    plan.old_length
        += odb::Point::manhattanDistance(request.begin, request.end);
  }
}

inline void resolveVias(const Recipe& recipe, Plan& plan)
{
  std::set<frVia*> unique;
  for (const auto& request : recipe.remove_vias) {
    frVia* match = nullptr;
    for (const auto& via : plan.net->getVias()) {
      if (via->getOrigin() == request.origin
          && via->getViaDef()->getName() == request.name) {
        if (match) {
          throw std::runtime_error("Ambiguous route trial via");
        }
        match = via.get();
      }
    }
    if (!match || !unique.insert(match).second
        || match->getViaDef()->isMultiCut()) {
      throw std::runtime_error("Unmatched or unsupported route trial via");
    }
    plan.vias.push_back(match);
  }
}

inline Plan resolveSource(frDesign* design, const Recipe& recipe)
{
  auto* net = design->getTopBlock()->findNet(recipe.net);
  if (!net || net->isFixed() || net->isClock() || net->isSpecial()
      || net->hasNDR() || net->hasJumpers() || net->isFake()
      || net->isWatermark()) {
    throw std::runtime_error("Unsupported route trial net");
  }
  if ((!recipe.tree && recipe.anchors.size() != 2)
      || (recipe.tree
          && (recipe.anchors.size() < 2 || recipe.anchors.size() > 8
              || recipe.remove_segments.size() > 64
              || recipe.remove_vias.size() > 32))) {
    throw std::runtime_error("Route trial exceeds source topology bounds");
  }
  Plan plan{net};
  plan.tree = recipe.tree;
  auto* tech = design->getTech();
  for (const auto& [name, x, y] : recipe.anchors) {
    plan.anchors.emplace(layerNumber(tech, name), x, y);
  }
  resolveSegments(tech, recipe, plan);
  resolveVias(recipe, plan);
  std::vector<std::pair<Node, Node>> old_edges;
  for (auto* segment : plan.segments) {
    const auto layer = segment->getLayerNum();
    const auto [a, b] = segment->getPoints();
    old_edges.push_back({{layer, a.x(), a.y()}, {layer, b.x(), b.y()}});
  }
  for (auto* via : plan.vias) {
    const auto p = via->getOrigin();
    old_edges.push_back({{via->getViaDef()->getLayer1Num(), p.x(), p.y()},
                         {via->getViaDef()->getLayer2Num(), p.x(), p.y()}});
  }
  if (!validTopology(old_edges, plan) || plan.vias.size() < 2) {
    throw std::runtime_error(
        "Route trial source has invalid anchored topology");
  }
  return plan;
}

inline bool viaOnTracks(frDesign* design,
                        const frViaDef* def,
                        const odb::Point& origin)
{
  if (def->getLayer2Num() != def->getLayer1Num() + 2) {
    return false;
  }
  for (const auto number : {def->getLayer1Num(), def->getLayer2Num()}) {
    const auto* layer = design->getTech()->getLayer(number);
    if (!layer->isRoutable() || layer->isBackside()) {
      return false;
    }
    bool found = false;
    const auto coord = layer->isHorizontal() ? origin.y() : origin.x();
    for (const auto& track : design->getTopBlock()->getTrackPatterns(number)) {
      if (track->isHorizontal() == layer->isHorizontal()) {
        continue;
      }
      const int64_t delta = int64_t(coord) - track->getStartCoord();
      const int64_t step = track->getTrackSpacing();
      if (step > 0 && delta >= 0 && delta % step == 0
          && delta / step < track->getNumTracks()) {
        found = true;
        break;
      }
    }
    if (!found) {
      return false;
    }
  }
  return true;
}

inline Plan resolve(frDesign* design, const Recipe& recipe)
{
  auto plan = resolveSource(design, recipe);
  auto* tech = design->getTech();
  if (recipe.add_segments.empty() || recipe.add_segments.size() > 64
      || recipe.add_vias.size() >= plan.vias.size()) {
    throw std::runtime_error("Route trial exceeds bounded path or via budget");
  }
  std::vector<std::pair<Node, Node>> new_edges;
  for (const auto& request : recipe.add_segments) {
    const auto layer = layerNumber(tech, request.layer);
    auto segment = std::make_unique<frPathSeg>();
    segment->setLayerNum(layer);
    segment->setPoints(request.begin, request.end);
    segment->setStyle(tech->getLayer(layer)->getDefaultSegStyle());
    if (tech->getLayer(layer)->isUnidirectional()
        && tech->getLayer(layer)->isHorizontal() == segment->isVertical()) {
      throw std::runtime_error("Route trial violates layer direction");
    }
    if (!via_cleanup::onTrack(design->getTopBlock(), segment.get(), layer)) {
      throw std::runtime_error("Route trial is off track");
    }
    segment->addToNet(plan.net);
    const auto [a, b] = segment->getPoints();
    new_edges.push_back({{layer, a.x(), a.y()}, {layer, b.x(), b.y()}});
    plan.new_length += odb::Point::manhattanDistance(a, b);
    plan.additions.push_back(std::move(segment));
  }
  std::set<std::pair<std::string, odb::Point>> unique_vias;
  for (const auto& request : recipe.add_vias) {
    auto* def = tech->getVia(request.name);
    if (!def || def->isMultiCut() || !viaOnTracks(design, def, request.origin)
        || !unique_vias.emplace(request.name, request.origin).second) {
      throw std::runtime_error("Unsupported or duplicate added via");
    }
    layerNumber(tech, tech->getLayer(def->getLayer1Num())->getName());
    layerNumber(tech, tech->getLayer(def->getLayer2Num())->getName());
    auto via = std::make_unique<frVia>(def, request.origin);
    via->addToNet(plan.net);
    const auto p = request.origin;
    new_edges.push_back({{def->getLayer1Num(), p.x(), p.y()},
                         {def->getLayer2Num(), p.x(), p.y()}});
    plan.added_vias.push_back(std::move(via));
  }
  if (!validTopology(new_edges, plan) || plan.new_length > plan.old_length) {
    throw std::runtime_error(
        "Route trial must preserve topology, not increase length, and remove "
        "vias");
  }
  return plan;
}

inline bool atAnchor(const Plan& plan, frLayerNum layer, const odb::Point& p)
{
  return plan.anchors.count({layer, p.x(), p.y()});
}

// Both removed and new geometry may meet other routes only at retained ports.
// Query all layers because routed vias are indexed only on their cut layer.
inline bool contactsAtAnchors(frRegionQuery* rq,
                              const Plan& plan,
                              frLayerNum layer,
                              const odb::Rect& box,
                              const odb::Rect& centerline)
{
  frRegionQuery::Objects<frBlockObject> fixed;
  rq->query(box, layer, fixed);
  if (!fixed.empty()) {
    return false;
  }
  std::vector<frBlockObject*> objects;
  rq->queryDRObj(box, objects);
  for (auto* object : objects) {
    if (object->typeId() == frcVia) {
      auto* via = static_cast<frVia*>(object);
      if (std::find(plan.vias.begin(), plan.vias.end(), via) != plan.vias.end()
          || via_cleanup::otherLayer(via, layer) < 0
          || !via_cleanup::metalBox(via, layer).intersects(box)) {
        continue;
      }
      const auto p = via->getOrigin();
      if (via->getNet() != plan.net || !atAnchor(plan, layer, p)
          || !centerline.intersects(odb::Rect(p, p))) {
        return false;
      }
      continue;
    }
    auto* shape = static_cast<frShape*>(object);
    if (shape->getLayerNum() != layer) {
      continue;
    }
    if (shape->typeId() != frcPathSeg) {
      return false;
    }
    auto* segment = static_cast<frPathSeg*>(shape);
    if (std::find(plan.segments.begin(), plan.segments.end(), segment)
        != plan.segments.end()) {
      continue;
    }
    const odb::Rect other(segment->getBeginPoint(), segment->getEndPoint());
    if (segment->getNet() != plan.net || !other.intersects(centerline)) {
      return false;
    }
    const odb::Rect intersection(std::max(centerline.xMin(), other.xMin()),
                                 std::max(centerline.yMin(), other.yMin()),
                                 std::min(centerline.xMax(), other.xMax()),
                                 std::min(centerline.yMax(), other.yMax()));
    if (intersection.dx() != 0 || intersection.dy() != 0
        || !atAnchor(plan, layer, intersection.ll())) {
      return false;
    }
  }
  return true;
}

inline bool sourceIsolated(frRegionQuery* rq, const Plan& plan)
{
  for (auto* segment : plan.segments) {
    if (!contactsAtAnchors(
            rq,
            plan,
            segment->getLayerNum(),
            segment->getBBox(),
            odb::Rect(segment->getBeginPoint(), segment->getEndPoint()))) {
      return false;
    }
  }
  for (auto* via : plan.vias) {
    const odb::Rect point(via->getOrigin(), via->getOrigin());
    for (const auto layer :
         {via->getViaDef()->getLayer1Num(), via->getViaDef()->getLayer2Num()}) {
      if (!contactsAtAnchors(
              rq, plan, layer, via_cleanup::metalBox(via, layer), point)) {
        return false;
      }
    }
  }
  return true;
}

inline bool targetIsolated(frRegionQuery* rq, const Plan& plan)
{
  for (const auto& via : plan.added_vias) {
    const odb::Rect point(via->getOrigin(), via->getOrigin());
    for (const auto layer :
         {via->getViaDef()->getLayer1Num(), via->getViaDef()->getLayer2Num()}) {
      if (!contactsAtAnchors(rq,
                             plan,
                             layer,
                             via_cleanup::metalBox(via.get(), layer),
                             point)) {
        return false;
      }
    }
  }
  for (const auto& shape : plan.additions) {
    const auto* segment = static_cast<frPathSeg*>(shape.get());
    if (!contactsAtAnchors(
            rq,
            plan,
            segment->getLayerNum(),
            segment->getBBox(),
            odb::Rect(segment->getBeginPoint(), segment->getEndPoint()))) {
      return false;
    }
  }
  return true;
}

inline bool isolated(frRegionQuery* rq, const Plan& plan)
{
  return sourceIsolated(rq, plan) && targetIsolated(rq, plan);
}

class Trial
{
 public:
  Trial(frRegionQuery* rq, Plan& plan) : rq_(rq), plan_(plan)
  {
    for (auto* segment : plan_.segments) {
      rq_->removeDRObj(segment);
    }
    for (auto* via : plan_.vias) {
      rq_->removeDRObj(via);
    }
    for (auto& shape : plan_.additions) {
      rq_->addDRObj(shape.get());
    }
    for (auto& via : plan_.added_vias) {
      rq_->addDRObj(via.get());
    }
  }
  Trial(const Trial&) = delete;
  Trial& operator=(const Trial&) = delete;
  ~Trial()
  {
    if (!committed_) {
      for (auto& shape : plan_.additions) {
        rq_->removeDRObj(shape.get());
      }
      for (auto& via : plan_.added_vias) {
        rq_->removeDRObj(via.get());
      }
      for (auto* segment : plan_.segments) {
        rq_->addDRObj(segment);
      }
      for (auto* via : plan_.vias) {
        rq_->addDRObj(via);
      }
    }
  }
  void commit()
  {
    for (auto* segment : plan_.segments) {
      plan_.net->removeShape(segment);
    }
    for (auto* via : plan_.vias) {
      plan_.net->removeVia(via);
    }
    for (auto& shape : plan_.additions) {
      plan_.net->addShape(std::move(shape));
    }
    for (auto& via : plan_.added_vias) {
      plan_.net->addVia(std::move(via));
    }
    plan_.net->setModified(true);
    committed_ = true;
  }

 private:
  frRegionQuery* rq_;
  Plan& plan_;
  bool committed_ = false;
};

inline bool apply(frDesign* design,
                  utl::Logger* logger,
                  RouterConfiguration* cfg,
                  Plan& plan)
{
  auto* rq = design->getRegionQuery();
  if (!sourceIsolated(rq, plan)) {
    logger->report("ROUTE_TRIAL reject_stage=source_contact");
    return false;
  }
  if (!targetIsolated(rq, plan)) {
    logger->report("ROUTE_TRIAL reject_stage=target_contact");
    return false;
  }
  odb::Rect box;
  box.mergeInit();
  for (auto* segment : plan.segments) {
    box.merge(segment->getBBox());
  }
  for (auto* via : plan.vias) {
    box.merge(via->getBBox());
  }
  for (const auto& shape : plan.additions) {
    box.merge(shape->getBBox());
  }
  for (const auto& via : plan.added_vias) {
    box.merge(via->getBBox());
  }
  Trial trial(rq, plan);
  odb::Rect drc_box, ext_box;
  box.bloat(cfg->DRCSAFEDIST, drc_box);
  box.bloat(std::max(cfg->DRCSAFEDIST, cfg->MTSAFEDIST), ext_box);
  FlexGCWorker checker(design->getTech(), logger, cfg);
  checker.setDrcBox(drc_box);
  checker.setExtBox(ext_box);
  checker.setEnableSurgicalFix(false);
  checker.init(design);
  checker.main();
  if (!checker.getMarkers().empty()) {
    logger->report("ROUTE_TRIAL reject_stage=local_drc markers={}",
                   checker.getMarkers().size());
    return false;
  }
  trial.commit();
  return true;
}

}  // namespace drt::route_trial
