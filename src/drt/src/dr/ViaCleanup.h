// SPDX-License-Identifier: BSD-3-Clause
// Experimental, bounded post-route layer reassignment. Disabled by default.
#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <vector>

#include "db/obj/frNet.h"
#include "db/obj/frShape.h"
#include "db/obj/frVia.h"
#include "frDesign.h"
#include "frRegionQuery.h"
#include "gc/FlexGC.h"

namespace drt::via_cleanup {

struct Candidate
{
  frPathSeg* segment;
  std::array<frVia*, 2> vias;
  frLayerNum target;
  int64_t length;
};

inline frLayerNum otherLayer(const frVia* via, frLayerNum layer)
{
  const auto* def = via->getViaDef();
  if (def->getLayer1Num() == layer) {
    return def->getLayer2Num();
  }
  return def->getLayer2Num() == layer ? def->getLayer1Num() : -1;
}

inline odb::Rect metalBox(const frVia* via, frLayerNum layer)
{
  return via->getViaDef()->getLayer1Num() == layer ? via->getLayer1BBox()
                                                   : via->getLayer2BBox();
}

inline bool nearSeam(frBlock* block, const odb::Point& point)
{
  const auto& grids = block->getGCellPatterns();
  for (int axis = 0; axis < 2; ++axis) {
    const auto& grid = grids[axis];
    const int64_t pitch = grid.getSpacing();
    const int64_t coord = axis == 0 ? point.x() : point.y();
    const int64_t pos = coord - grid.getStartCoord();
    const int64_t period = 7 * pitch;
    const int64_t rem = (pos % period + period) % period;
    // Exclude artificial repetitions outside the die's internal window seams.
    const int64_t seam = pos - rem + (rem > period / 2 ? period : 0);
    if (seam > 0 && seam / pitch < grid.getCount()
        && 2 * std::min(rem, period - rem) < pitch) {
      return true;
    }
  }
  return false;
}

inline bool onTrack(frBlock* block, const frPathSeg* segment, frLayerNum layer)
{
  const auto point = segment->getBeginPoint();
  const int64_t coord = segment->isVertical() ? point.x() : point.y();
  for (const auto& track : block->getTrackPatterns(layer)) {
    if (track->isHorizontal() != segment->isVertical()) {
      continue;
    }
    const int64_t delta = coord - track->getStartCoord();
    const int64_t spacing = track->getTrackSpacing();
    if (spacing > 0 && delta >= 0 && delta % spacing == 0
        && delta / spacing < track->getNumTracks()) {
      return true;
    }
  }
  return false;
}

// The removed-layer metal must be a leaf chain with exactly the two vias.
// Conservatively reject ALL fixed pin/blockage contact and any other routed
// object, including same-net branches whose endpoints lie inside the segment.
inline bool isolated(frRegionQuery* rq, const Candidate& candidate)
{
  const auto layer = candidate.segment->getLayerNum();
  const std::array<odb::Rect, 3> boxes{candidate.segment->getBBox(),
                                       metalBox(candidate.vias[0], layer),
                                       metalBox(candidate.vias[1], layer)};
  for (const auto& box : boxes) {
    frRegionQuery::Objects<frBlockObject> objects;
    rq->query(box, layer, objects);
    if (!objects.empty()) {
      return false;
    }
    // Routed vias are indexed on their CUT layer, not either metal layer.
    // Query all layers, then test the metal actually removed by this trial.
    std::vector<frBlockObject*> routed;
    rq->queryDRObj(box, routed);
    for (auto* obj : routed) {
      if (obj == candidate.segment || obj == candidate.vias[0]
          || obj == candidate.vias[1]) {
        continue;
      }
      if (obj->typeId() == frcVia) {
        auto* via = static_cast<frVia*>(obj);
        if (otherLayer(via, layer) >= 0
            && metalBox(via, layer).intersects(box)) {
          return false;
        }
      } else if (static_cast<frShape*>(obj)->getLayerNum() == layer) {
        return false;
      }
    }
  }
  // The connectivity checker gives pin connections special endpoint styles.
  // Avoid fixed contacts anywhere on the replacement, including its interior.
  frPathSeg replacement(*candidate.segment);
  replacement.setLayerNum(candidate.target);
  replacement.setStyle(rq->getDesign()
                           ->getTech()
                           ->getLayer(candidate.target)
                           ->getDefaultSegStyle());
  frRegionQuery::Objects<frBlockObject> fixed;
  rq->query(replacement.getBBox(), candidate.target, fixed);
  if (!fixed.empty()) {
    return false;
  }
  for (const auto* via : candidate.vias) {
    frRegionQuery::Objects<frBlockObject> objects;
    rq->query(metalBox(via, candidate.target), candidate.target, objects);
    if (!objects.empty()) {
      return false;
    }
  }
  return true;
}

// Adding a junction inside the replacement creates a cycle that the global
// connectivity checker may prune, emitting new Recheck markers after local GC.
// Permit only target connections centered at the original via endpoints.
inline bool targetClear(frRegionQuery* rq, const Candidate& c)
{
  frPathSeg replacement(*c.segment);
  replacement.setLayerNum(c.target);
  replacement.setStyle(
      rq->getDesign()->getTech()->getLayer(c.target)->getDefaultSegStyle());
  const auto [begin, end] = replacement.getPoints();
  const odb::Rect line(begin, end);
  std::vector<frBlockObject*> objects;
  rq->queryDRObj(replacement.getBBox(), objects);
  for (auto* obj : objects) {
    if (obj == c.segment || obj == c.vias[0] || obj == c.vias[1]) {
      continue;
    }
    if (obj->typeId() == frcVia) {
      auto* via = static_cast<frVia*>(obj);
      if (otherLayer(via, c.target) >= 0
          && metalBox(via, c.target).intersects(replacement.getBBox())
          && via->getOrigin() != begin && via->getOrigin() != end) {
        return false;
      }
      continue;
    }
    auto* shape = static_cast<frShape*>(obj);
    if (shape->getLayerNum() != c.target) {
      continue;
    }
    if (shape->typeId() != frcPathSeg) {
      return false;
    }
    auto* segment = static_cast<frPathSeg*>(shape);
    const auto [bp, ep] = segment->getPoints();
    const odb::Rect other(bp, ep);
    // Reject side contacts, interior crossings, and collinear overlaps.
    if (!other.intersects(line)) {
      return false;
    }
    const odb::Rect intersection(std::max(begin.x(), bp.x()),
                                 std::max(begin.y(), bp.y()),
                                 std::min(end.x(), ep.x()),
                                 std::min(end.y(), ep.y()));
    if (intersection.dx() != 0 || intersection.dy() != 0
        || (intersection.ll() != begin && intersection.ll() != end)) {
      return false;
    }
  }
  return true;
}

inline std::vector<Candidate> candidates(frDesign* design,
                                         RouterConfiguration* cfg)
{
  auto* block = design->getTopBlock();
  auto* tech = design->getTech();
  const auto& grids = block->getGCellPatterns();
  std::vector<Candidate> result;
  if (grids.size() != 2 || grids[0].getSpacing() <= 0
      || grids[1].getSpacing() <= 0) {
    return result;
  }
  for (const auto& net : block->getNets()) {
    if (net->isFixed() || net->hasNDR() || net->isClock() || net->isSpecial()
        || net->isFake() || net->isWatermark() || net->hasJumpers()) {
      continue;
    }
    std::map<odb::Point, std::vector<frVia*>> by_point;
    for (const auto& via : net->getVias()) {
      if (!via->getViaDef()->isMultiCut()) {
        by_point[via->getOrigin()].push_back(via.get());
      }
    }
    for (const auto& shape : net->getShapes()) {
      if (shape->typeId() != frcPathSeg) {
        continue;
      }
      auto* seg = static_cast<frPathSeg*>(shape.get());
      const auto [begin, end] = seg->getPoints();
      const auto length = odb::Point::manhattanDistance(begin, end);
      if (seg->isApPathSeg() || length == 0
          || (begin.x() != end.x() && begin.y() != end.y())
          || length > grids[seg->isVertical() ? 1 : 0].getSpacing()
          || (!nearSeam(block, begin) && !nearSeam(block, end))) {
        continue;
      }
      const auto first = by_point.find(begin);
      const auto last = by_point.find(end);
      if (first == by_point.end() || last == by_point.end()) {
        continue;
      }
      for (auto* v1 : first->second) {
        const auto target = otherLayer(v1, seg->getLayerNum());
        if (target < cfg->BOTTOM_ROUTING_LAYER
            || target > cfg->TOP_ROUTING_LAYER) {
          continue;
        }
        auto* layer = tech->getLayer(target);
        if (!layer->isRoutable() || layer->isBackside()
            || (layer->isUnidirectional()
                && layer->isHorizontal() == seg->isVertical())
            || !onTrack(block, seg, target)) {
          continue;
        }
        for (auto* v2 : last->second) {
          if (otherLayer(v2, seg->getLayerNum()) == target) {
            result.push_back({seg, {v1, v2}, target, length});
          }
        }
      }
    }
  }
  // Stable tie order is original net/shape order; never pointer order.
  std::stable_sort(
      result.begin(), result.end(), [](const auto& a, const auto& b) {
        return a.length < b.length;
      });
  return result;
}

// Trial shapes live only in the region query. Original owners, IDs, list order,
// markers and modified flags remain untouched until acceptance. GC reads DR
// objects through the region query; it does not need trial objects in net
// lists.
class Trial
{
 public:
  Trial(frRegionQuery* rq, const Candidate& c, const frSegStyle& style)
      : rq_(rq),
        c_(c),
        old_layer_(c.segment->getLayerNum()),
        old_style_(c.segment->getStyle())
  {
    auto* net = c.segment->getNet();
    // Preserve target-layer landing metal exactly, including off-center
    // contacts.
    for (auto* via : c.vias) {
      const auto* def = via->getViaDef();
      const auto& figs = def->getLayer1Num() == c.target ? def->getLayer1Figs()
                                                         : def->getLayer2Figs();
      for (const auto& fig : figs) {
        auto patch = std::make_unique<frPatchWire>();
        patch->setOrigin(via->getOrigin());
        patch->setOffsetBox(fig->getBBox());
        patch->setLayerNum(c.target);
        patch->addToNet(net);
        patches_.push_back(std::move(patch));
      }
    }
    rq_->removeDRObj(c.segment);
    for (auto* via : c.vias) {
      rq_->removeDRObj(via);
    }
    c.segment->setLayerNum(c.target);
    c.segment->setStyle(style);
    rq_->addDRObj(c.segment);
    for (const auto& patch : patches_) {
      rq_->addDRObj(patch.get());
    }
  }
  Trial(const Trial&) = delete;
  Trial& operator=(const Trial&) = delete;
  ~Trial()
  {
    if (committed_) {
      return;
    }
    for (const auto& patch : patches_) {
      rq_->removeDRObj(patch.get());
    }
    rq_->removeDRObj(c_.segment);
    c_.segment->setLayerNum(old_layer_);
    c_.segment->setStyle(old_style_);
    rq_->addDRObj(c_.segment);
    for (auto* via : c_.vias) {
      rq_->addDRObj(via);
    }
  }
  void commit()
  {
    auto* net = c_.segment->getNet();
    for (auto& patch : patches_) {
      net->addPatchWire(std::move(patch));
    }
    for (auto* via : c_.vias) {
      net->removeVia(via);
    }
    net->setModified(true);
    committed_ = true;
  }

 private:
  frRegionQuery* rq_;
  Candidate c_;
  frLayerNum old_layer_;
  frSegStyle old_style_;
  std::vector<std::unique_ptr<frPatchWire>> patches_;
  bool committed_ = false;
};

inline int run(frDesign* design, utl::Logger* logger, RouterConfiguration* cfg)
{
  const auto start = std::chrono::steady_clock::now();
  auto* rq = design->getRegionQuery();
  const auto work = candidates(design, cfg);
  std::set<frVia*> removed;
  std::set<frPathSeg*> moved;
  int checked = 0;
  int accepted = 0;
  int rejected_topology = 0;
  constexpr int max_checks = 512;
  for (const auto& c : work) {
    if (checked == max_checks) {
      break;
    }
    if (moved.count(c.segment) || removed.count(c.vias[0])
        || removed.count(c.vias[1])) {
      continue;
    }
    // Bounding boxes are sufficient only for rectangular via metal.
    bool rectangles = true;
    for (const auto* via : c.vias) {
      for (const auto& fig : via->getViaDef()->getLayer1Figs()) {
        rectangles &= fig->typeId() == frcRect;
      }
      for (const auto& fig : via->getViaDef()->getLayer2Figs()) {
        rectangles &= fig->typeId() == frcRect;
      }
    }
    if (!rectangles || !isolated(rq, c) || !targetClear(rq, c)) {
      ++rejected_topology;
      continue;
    }
    odb::Rect box = c.segment->getBBox();
    for (const auto* via : c.vias) {
      box.merge(via->getBBox());
    }
    Trial trial(
        rq, c, design->getTech()->getLayer(c.target)->getDefaultSegStyle());
    box.merge(c.segment->getBBox());
    odb::Rect drc_box, ext_box;
    box.bloat(cfg->DRCSAFEDIST, drc_box);
    box.bloat(std::max(cfg->MTSAFEDIST, cfg->DRCSAFEDIST), ext_box);
    FlexGCWorker checker(design->getTech(), logger, cfg);
    checker.setDrcBox(drc_box);
    checker.setExtBox(ext_box);
    checker.setEnableSurgicalFix(false);
    checker.init(design);
    checker.main();
    ++checked;
    if (!checker.getMarkers().empty()) {
      continue;
    }
    removed.insert(c.vias[0]);
    removed.insert(c.vias[1]);
    moved.insert(c.segment);
    trial.commit();
    ++accepted;
  }
  const double seconds
      = std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
  logger->report(
      "VIA_CLEANUP candidates={} topology_rejected={} checked={} accepted={} "
      "removed={} seconds={:.6f}",
      work.size(),
      rejected_topology,
      checked,
      accepted,
      2 * accepted,
      seconds);
  logger->metric("route__via_cleanup__accepted", accepted);
  logger->metric("route__via_cleanup__seconds", seconds);
  return accepted;
}

}  // namespace drt::via_cleanup
