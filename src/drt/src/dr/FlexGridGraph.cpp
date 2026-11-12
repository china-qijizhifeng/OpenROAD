// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "dr/FlexGridGraph.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <vector>

#include "db/obj/frTrackPattern.h"
#include "db/tech/frViaDef.h"
#include "dr/FlexDR.h"
#include "frBaseTypes.h"
#include "frDesign.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"

using odb::dbTechLayerDir;
using odb::dbTechLayerType;

namespace drt {
void FlexGridGraph::addAccessPointLocation(frLayerNum layer_num,
                                           frCoord x_coord,
                                           frCoord y_coord)
{
  ap_locs_[layer_num].insert(odb::Point(x_coord, y_coord));
}

void initCoords(const frLayerCoordTrackPatternMap& map,
                std::vector<frCoord>& coords)
{
  coords.clear();
  for (auto& [l, m] : map) {
    coords.reserve(m.size());
    if (coords.empty()) {
      std::ranges::transform(m, std::back_inserter(coords), [](const auto& kv) {
        return kv.first;
      });
    } else {
      auto it = coords.begin();
      for (auto& [k, _] : m) {
        it = std::lower_bound(it, coords.end(), k);
        if (it == coords.end() || *it != k) {
          it = coords.insert(it, k);
        }
      }
    }
  }
}

bool FlexGridGraph::isAccessPointLocation(frLayerNum layer_num,
                                          frCoord x_coord,
                                          frCoord y_coord) const
{
  if (ap_locs_.size() <= layer_num) {
    return false;
  }
  const auto& layer_maze_locs = ap_locs_[layer_num];
  return layer_maze_locs.find(odb::Point(x_coord, y_coord))
         != layer_maze_locs.end();
}
void FlexGridGraph::initGrids(const frLayerCoordTrackPatternMap& xMap,
                              const frLayerCoordTrackPatternMap& yMap,
                              const frLayerDirMap& zMap,
                              bool followGuide)
{
  // initialize coord vectors
  initCoords(xMap, xCoords_);
  initCoords(yMap, yCoords_);
  zCoords_.clear();
  zHeights_.clear();
  layerRouteDirections_.clear();

  frCoord zHeight = 0;
  // std::vector<frCoord> via2viaMinLenTmp(4, 0);
  zCoords_.reserve(zMap.size());
  for (auto& [k, v] : zMap) {
    zCoords_.push_back(k);
    zHeight += getTech()->getLayer(k)->getPitch() * router_cfg_->VIACOST;
    zHeights_.push_back(zHeight);
    layerRouteDirections_.push_back(v);
  }
  // initialize all grids
  frMIdx xDim, yDim, zDim;
  getDim(xDim, yDim, zDim);
  const int capacity = xDim * yDim * zDim;

  nodes_.clear();
  nodes_.resize(capacity, Node());
  // new
  prevDirs_.clear();
  srcs_.clear();
  dsts_.clear();

  prevDirs_.resize(static_cast<std::size_t>(capacity) * 3, false);
  srcs_.resize(capacity, false);
  dsts_.resize(capacity, false);
  guides_.clear();
  if (followGuide) {
    guides_.resize(capacity, false);
  } else {
    guides_.resize(capacity, true);
  }
}

bool FlexGridGraph::outOfDieVia(frMIdx x,
                                frMIdx y,
                                frMIdx z,
                                const odb::Rect& dieBox)
{
  frLayerNum lNum = getLayerNum(z) + 1;
  if (lNum > getTech()->getTopLayerNum()) {
    return false;
  }
  const frViaDef* via = getTech()->getLayer(lNum)->getDefaultViaDef();
  if (!via) {
    return true;
  }
  odb::Rect viaBox(via->getLayer1ShapeBox());
  viaBox.merge(via->getLayer2ShapeBox());
  viaBox.moveDelta(xCoords_[x], yCoords_[y]);
  return !dieBox.contains(viaBox);
}
bool FlexGridGraph::hasOutOfDieViol(frMIdx x, frMIdx y, frMIdx z)
{
  const frLayerNum lNum = getLayerNum(z);
  if (!getTech()->getLayer(lNum)->isUnidirectional()) {
    return false;
  }
  odb::Rect testBoxUp;
  if (lNum + 1 <= getTech()->getTopLayerNum()) {
    const frViaDef* via = getTech()->getLayer(lNum + 1)->getDefaultViaDef();
    if (via) {
      testBoxUp = via->getLayer1ShapeBox();
      testBoxUp.merge(via->getLayer2ShapeBox());
      testBoxUp.moveDelta(xCoords_[x], yCoords_[y]);
    } else {
      // artificial value to indicate no via in test below
      dieBox_.bloat(1, testBoxUp);
    }
  }
  odb::Rect testBoxDown;
  if (lNum - 1 >= getTech()->getBottomLayerNum()) {
    const frViaDef* via = getTech()->getLayer(lNum - 1)->getDefaultViaDef();
    if (via) {
      testBoxDown = via->getLayer1ShapeBox();
      testBoxDown.merge(via->getLayer2ShapeBox());
      testBoxDown.moveDelta(xCoords_[x], yCoords_[y]);
    } else {
      // artificial value to indicate no via in test below
      dieBox_.bloat(1, testBoxDown);
    }
  }
  if (getTech()->getLayer(lNum)->isVertical()) {
    return (testBoxUp.xMax() > dieBox_.xMax()
            || testBoxUp.xMin() < dieBox_.xMin())
           && (testBoxDown.xMax() > dieBox_.xMax()
               || testBoxDown.xMin() < dieBox_.xMin());
  }
  // layer is horizontal
  return (testBoxUp.yMax() > dieBox_.yMax()
          || testBoxUp.yMin() < dieBox_.yMin())
         && (testBoxDown.yMax() > dieBox_.yMax()
             || testBoxDown.yMin() < dieBox_.yMin());
}

bool FlexGridGraph::isWorkerBorder(frMIdx v, bool isVert)
{
  if (isVert) {
    return xCoords_[v] == drWorker_->getRouteBox().xMin()
           || xCoords_[v] == drWorker_->getRouteBox().xMax();
  }
  return yCoords_[v] == drWorker_->getRouteBox().yMin()
         || yCoords_[v] == drWorker_->getRouteBox().yMax();
}
bool FlexGridGraph::hasAlignedUpDefTrack(
    frLayerNum layerNum,
    const std::map<frLayerNum, frTrackPattern*>& xSubMap,
    const std::map<frLayerNum, frTrackPattern*>& ySubMap) const
{
  for (frLayerNum lNum = layerNum + 2;
       lNum < (int) getTech()->getLayers().size();
       lNum += 2) {
    auto it = xSubMap.find(lNum);
    if (it != xSubMap.end()) {  // has x track in lNum
      if (it->second) {         // has track pattern, i.e., the track is default
        return true;
      }
    }
    it = ySubMap.find(lNum);
    if (it != ySubMap.end()) {  // has y track in lNum
      if (it->second) {         // has track pattern, i.e., the track is default
        return true;
      }
    }
  }
  return false;
}

void FlexGridGraph::initEdges(const frDesign* design,
                              frLayerCoordTrackPatternMap& xMap,
                              frLayerCoordTrackPatternMap& yMap,
                              const frLayerDirMap& zMap,
                              const odb::Rect& bbox,
                              bool initDR)
{
  frMIdx xDim, yDim, zDim;
  getDim(xDim, yDim, zDim);
  // initialize grid graph
  frMIdx zIdx = 0;
  dieBox_ = design->getTopBlock()->getDieBox();
  // Per-index parts of the tests made at every grid point. isEdgeInBox: both
  // end points of an edge lie in bbox, tested per x and per y index.
  // outOfDieVia: the die box contains the moved via box, also tested per x and
  // per y index. Track lookups on the layer two above depend on one index.
  std::vector<uint8_t> xInBox(xDim);
  std::vector<uint8_t> yInBox(yDim);
  for (frMIdx xIdx = 0; xIdx < xDim; xIdx++) {
    xInBox[xIdx]
        = bbox.xMin() <= xCoords_[xIdx] && xCoords_[xIdx] <= bbox.xMax();
  }
  for (frMIdx yIdx = 0; yIdx < yDim; yIdx++) {
    yInBox[yIdx]
        = bbox.yMin() <= yCoords_[yIdx] && yCoords_[yIdx] <= bbox.yMax();
  }
  // addEdge for the E, N and U edges at (x, y, z) without the per-call
  // direction dispatch: isValid of both end points and isEdgeInBox.
  auto addEdgeE = [&](frMIdx x, frMIdx y, frMIdx z) {
    if (x + 1 < xDim && xInBox[x] && xInBox[x + 1] && yInBox[y]) {
      nodes_[getIdx(x, y, z)].hasEastEdge = true;
    }
  };
  auto addEdgeN = [&](frMIdx x, frMIdx y, frMIdx z) {
    if (y + 1 < yDim && yInBox[y] && yInBox[y + 1] && xInBox[x]) {
      nodes_[getIdx(x, y, z)].hasNorthEdge = true;
    }
  };
  auto addEdgeU = [&](frMIdx x, frMIdx y, frMIdx z) {
    if (z + 1 < zDim && xInBox[x] && yInBox[y]) {
      nodes_[getIdx(x, y, z)].hasUpEdge = true;
    }
  };
  std::vector<uint8_t> xViaInDie(xDim);
  std::vector<uint8_t> yViaInDie(yDim);
  // Track at an index on the layer two above: 0 none, 1 track without
  // pattern, 2 track with pattern.
  std::vector<uint8_t> upTrack(std::max(xDim, yDim));
  for (const auto& [layerNum, dir] : zMap) {
    frLayerNum nonPrefLayerNum;
    const auto layer = getTech()->getLayer(layerNum);
    if (layerNum + 2 <= router_cfg_->TOP_ROUTING_LAYER) {
      nonPrefLayerNum = layerNum + 2;
    } else if (layerNum - 2 >= getTech()->getBottomLayerNum()) {
      nonPrefLayerNum = layerNum - 2;
    } else {
      nonPrefLayerNum = layerNum;
    }
    auto& yLayerMap = yMap[layerNum];
    auto& yLayer2Map = yMap[layerNum + 2];
    auto& yNonPrefLayerMap = yMap[nonPrefLayerNum];
    auto& xLayerMap = xMap[layerNum];
    auto& xLayer2Map = xMap[layerNum + 2];
    auto& xNonPrefLayerMap = xMap[nonPrefLayerNum];
    const bool inRoutingLayers = layerNum >= router_cfg_->BOTTOM_ROUTING_LAYER
                                 && layerNum <= router_cfg_->TOP_ROUTING_LAYER;
    const bool onGridOnly
        = layer->getLef58RightWayOnGridOnlyConstraint() != nullptr;
    // outOfDieVia(x, y, zIdx, dieBox_) for this layer.
    const frLayerNum viaLayerNum = getLayerNum(zIdx) + 1;
    const bool viaAboveTop = viaLayerNum > getTech()->getTopLayerNum();
    const frViaDef* dieVia
        = viaAboveTop ? nullptr
                      : getTech()->getLayer(viaLayerNum)->getDefaultViaDef();
    if (dieVia) {
      odb::Rect viaBox(dieVia->getLayer1ShapeBox());
      viaBox.merge(dieVia->getLayer2ShapeBox());
      for (frMIdx xIdx = 0; xIdx < xDim; xIdx++) {
        xViaInDie[xIdx] = dieBox_.xMin() <= viaBox.xMin() + xCoords_[xIdx]
                          && dieBox_.xMax() >= viaBox.xMax() + xCoords_[xIdx];
      }
      for (frMIdx yIdx = 0; yIdx < yDim; yIdx++) {
        yViaInDie[yIdx] = dieBox_.yMin() <= viaBox.yMin() + yCoords_[yIdx]
                          && dieBox_.yMax() >= viaBox.yMax() + yCoords_[yIdx];
      }
    }
    auto outOfDieViaAt = [&](frMIdx x, frMIdx y) {
      if (viaAboveTop) {
        return false;
      }
      if (dieVia == nullptr) {
        return true;
      }
      return !(xViaInDie[x] && yViaInDie[y]);
    };
    if (dir == dbTechLayerDir::HORIZONTAL) {
      for (frMIdx xIdx = 0; xIdx < xDim; xIdx++) {
        auto xIt2 = xLayer2Map.find(xCoords_[xIdx]);
        upTrack[xIdx] = xIt2 == xLayer2Map.end()  ? 0
                        : xIt2->second == nullptr ? 1
                                                  : 2;
      }
      for (frMIdx yIdx = 0; yIdx < yCoords_.size(); yIdx++) {
        auto yCoord = yCoords_[yIdx];
        auto yIt = yLayerMap.find(yCoord);
        if (yIt == yLayerMap.end()) {
          continue;
        }
        const bool rowHasPattern = yIt->second != nullptr;
        // add edge for preferred direction
        const bool prefEdges
            = inRoutingLayers && (!onGridOnly || rowHasPattern);
        const bool prefGridCost = !rowHasPattern || isWorkerBorder(yIdx, false);
        for (frMIdx xIdx = 0; xIdx < xCoords_.size(); xIdx++) {
          // add cost to out-of-die edge
          bool isOutOfDieVia = outOfDieViaAt(xIdx, yIdx);
          if (prefEdges
              && (!isOutOfDieVia || !hasOutOfDieViol(xIdx, yIdx, zIdx))) {
            addEdgeE(xIdx, yIdx, zIdx);
            if (prefGridCost) {
              setGridCostE(xIdx, yIdx, zIdx);
            }
          }
          if (isOutOfDieVia || upTrack[xIdx] == 0) {
            continue;
          }
          // via to upper layer
          const bool is_on_grid = rowHasPattern && upTrack[xIdx] == 2;
          if (is_on_grid || !onGridOnly
              || isAccessPointLocation(layerNum, xCoords_[xIdx], yCoord)) {
            addEdgeU(xIdx, yIdx, zIdx);
            if (!is_on_grid) {
              setGridCostU(xIdx, yIdx, zIdx);
            }
          }
        }
      }
      // get non pref track layer --> use upper layer pref dir track if
      // possible
      if (router_cfg_->USENONPREFTRACKS && !layer->isUnidirectional()) {
        for (frMIdx xIdx = 0; xIdx < xCoords_.size(); xIdx++) {
          auto xCoord = xCoords_[xIdx];
          auto xIt3 = xNonPrefLayerMap.find(xCoord);
          if (xIt3 == xNonPrefLayerMap.end()) {
            continue;
          }
          if (!inRoutingLayers) {
            continue;
          }
          // vertical non-pref track
          for (frMIdx yIdx = 0; yIdx < yCoords_.size(); yIdx++) {
            addEdgeN(xIdx, yIdx, zIdx);
            setGridCostN(xIdx, yIdx, zIdx);
          }
        }
      }
    } else if (dir == dbTechLayerDir::VERTICAL) {
      for (frMIdx yIdx = 0; yIdx < yDim; yIdx++) {
        auto yIt2 = yLayer2Map.find(yCoords_[yIdx]);
        upTrack[yIdx] = yIt2 == yLayer2Map.end()  ? 0
                        : yIt2->second == nullptr ? 1
                                                  : 2;
      }
      for (frMIdx xIdx = 0; xIdx < xCoords_.size(); xIdx++) {
        auto xCoord = xCoords_[xIdx];
        auto xIt = xLayerMap.find(xCoord);
        if (xIt == xLayerMap.end()) {
          continue;
        }
        const bool colHasPattern = xIt->second != nullptr;
        // add edge for preferred direction
        const bool prefEdges
            = inRoutingLayers && (!onGridOnly || colHasPattern);
        const bool prefGridCost = !colHasPattern || isWorkerBorder(xIdx, true);
        for (frMIdx yIdx = 0; yIdx < yCoords_.size(); yIdx++) {
          // add cost to out-of-die edge
          bool isOutOfDieVia = outOfDieViaAt(xIdx, yIdx);
          if (prefEdges
              && (!isOutOfDieVia || !hasOutOfDieViol(xIdx, yIdx, zIdx))) {
            addEdgeN(xIdx, yIdx, zIdx);
            if (prefGridCost) {
              setGridCostN(xIdx, yIdx, zIdx);
            }
          }
          if (isOutOfDieVia || upTrack[yIdx] == 0) {
            continue;
          }
          // via to upper layer
          const bool is_on_grid = colHasPattern && upTrack[yIdx] == 2;
          if (is_on_grid || !onGridOnly
              || isAccessPointLocation(layerNum, xCoord, yCoords_[yIdx])) {
            addEdgeU(xIdx, yIdx, zIdx);
            if (!is_on_grid) {
              setGridCostU(xIdx, yIdx, zIdx);
            }
          }
        }
      }
      // get non pref track layer --> use upper layer pref dir track if
      // possible
      if (router_cfg_->USENONPREFTRACKS && !layer->isUnidirectional()) {
        for (frMIdx yIdx = 0; yIdx < yCoords_.size(); yIdx++) {
          auto yCoord = yCoords_[yIdx];
          auto yIt3 = yNonPrefLayerMap.find(yCoord);
          if (yIt3 == yNonPrefLayerMap.end()) {
            continue;
          }
          if (!inRoutingLayers) {
            continue;
          }
          // horizontal non-pref track
          for (frMIdx xIdx = 0; xIdx < xCoords_.size(); xIdx++) {
            addEdgeE(xIdx, yIdx, zIdx);
            setGridCostE(xIdx, yIdx, zIdx);
          }
        }
      }
    }
    ++zIdx;
  }
  // this creates via edges over each ap until reaching a default track or a
  // layer with normal routing; in this case it creates jogs connections to the
  // neighboring tracks
  for (const odb::Point3D& apPt : drWorker_->getSpecialAccessAPs()) {
    for (int i = 0; i < 2; i++) {  // down and up
      bool up = (bool) i;
      int inc = up ? 1 : -1;
      frMIdx startZ = getMazeZIdx(apPt.z());
      frLayerNum nextLNum = getLayerNum(startZ) + 2 * inc;
      if (!up) {
        startZ--;
      }
      frMIdx xIdx = getMazeXIdx(apPt.x());
      frMIdx yIdx = getMazeYIdx(apPt.y());
      // create the edges
      for (int zIdx = startZ; zIdx >= 0 && zIdx < (int) zCoords_.size() - 1;
           zIdx += inc, nextLNum += inc * 2) {
        addEdge(xIdx, yIdx, zIdx, frDirEnum::U, bbox, initDR);
        frLayer* nextLayer = getTech()->getLayer(nextLNum);
        const bool restrictedRouting
            = nextLayer->isUnidirectional()
              || nextLNum < router_cfg_->BOTTOM_ROUTING_LAYER
              || nextLNum > router_cfg_->TOP_ROUTING_LAYER;
        if (!restrictedRouting || nextLayer->isVertical()) {
          auto& xSubMap = xMap[nextLNum];
          auto xTrack = xSubMap.find(apPt.x());
          if (xTrack != xSubMap.end() && xTrack->second != nullptr) {
            break;
          }
        }
        if (!restrictedRouting || nextLayer->isHorizontal()) {
          auto& ySubMap = yMap[nextLNum];
          auto yTrack = ySubMap.find(apPt.y());
          if (yTrack != ySubMap.end() && yTrack->second != nullptr) {
            break;
          }
        }
        // didn't find default track, then create tracks if possible
        if (!restrictedRouting
            && nextLNum >= router_cfg_->VIA_ACCESS_LAYERNUM) {
          dbTechLayerDir prefDir
              = design->getTech()->getLayer(nextLNum)->getDir();
          xMap[nextLNum][apPt.x()] = nullptr;  // to keep coherence
          yMap[nextLNum][apPt.y()] = nullptr;
          frMIdx nextZ = up ? zIdx + 1 : zIdx;
          // This is a value to make sure the edges we are adding will
          // reach a track on the layer of interest.  It is simpler to
          // be conservative than trying to figure out how many edges
          // to add to hit it precisely.  I intend to obviate the need
          // for this whole approach next.  Note that addEdge checks
          // for bounds so I don't.
          const int max_offset = 20;
          if (prefDir == dbTechLayerDir::HORIZONTAL) {
            for (int offset = 0; offset < max_offset; ++offset) {
              addEdge(xIdx, yIdx + offset, nextZ, frDirEnum::N, bbox, initDR);
              addEdge(xIdx, yIdx - offset, nextZ, frDirEnum::S, bbox, initDR);
            }
          } else {
            for (int offset = 0; offset < max_offset; ++offset) {
              addEdge(xIdx + offset, yIdx, nextZ, frDirEnum::E, bbox, initDR);
              addEdge(xIdx - offset, yIdx, nextZ, frDirEnum::W, bbox, initDR);
            }
          }
          break;
        }
      }
    }
  }
}

// initialization: update grid graph topology, does not assign edge cost
void FlexGridGraph::init(const frDesign* design,
                         const odb::Rect& routeBBox,
                         const odb::Rect& extBBox,
                         frLayerCoordTrackPatternMap& xMap,
                         frLayerCoordTrackPatternMap& yMap,
                         bool initDR,
                         bool followGuide)
{
  auto* via_data = getDRWorker()->getViaData();
  halfViaEncArea_ = &via_data->halfViaEncArea;

  // get tracks intersecting with the Maze bbox
  frLayerDirMap zMap;
  size_t layerCount = design->getTech()->getLayers().size();
  zMap.reserve(layerCount);

  initTracks(design, xMap, yMap, zMap, extBBox);
  initGrids(xMap, yMap, zMap, followGuide);  // buildGridGraph
  initEdges(
      design, xMap, yMap, zMap, routeBBox, initDR);  // add edges and edgeCost
  ap_locs_.clear();
}

// initialization helpers
// get all tracks intersecting with the Maze bbox, left/bottom are inclusive
void FlexGridGraph::initTracks(
    const frDesign* design,
    frLayerCoordTrackPatternMap& horLoc2TrackPatterns,
    frLayerCoordTrackPatternMap& vertLoc2TrackPatterns,
    frLayerDirMap& layerNum2PreRouteDir,
    const odb::Rect& bbox)
{
  for (auto& layer : getTech()->getLayers()) {
    if (layer->getType() != dbTechLayerType::ROUTING) {
      continue;
    }
    frLayerNum currLayerNum = layer->getLayerNum();
    // Layers above TOP_ROUTING_LAYER carry no routing edges or access vias
    // (initEdges and the special-access via loop already cap work at
    // TOP_ROUTING_LAYER), yet they were still added to the grid here, sizing
    // nodes_, the per-node arrays, and the maze search bounds across the full
    // layer stack. Skip them so the grid graph spans only the routing range,
    // cutting detailed-route memory and runtime when the design's top routing
    // layer is below the technology top. Default TOP_ROUTING_LAYER is INT_MAX,
    // so this is a no-op unless the top routing layer is explicitly lowered.
    if (currLayerNum > router_cfg_->TOP_ROUTING_LAYER) {
      continue;
    }
    dbTechLayerDir currPrefRouteDir = layer->getDir();
    for (auto& tp : design->getTopBlock()->getTrackPatterns(currLayerNum)) {
      // allow wrongway if global variable and design rule allow
      bool flag
          = (router_cfg_->USENONPREFTRACKS && !layer->isUnidirectional())
                ? (tp->isHorizontal()
                   && currPrefRouteDir == dbTechLayerDir::VERTICAL)
                      || (!tp->isHorizontal()
                          && currPrefRouteDir == dbTechLayerDir::HORIZONTAL)
                : true;
      if (flag) {
        int trackNum = ((tp->isHorizontal() ? bbox.xMin() : bbox.yMin())
                        - tp->getStartCoord())
                       / (int) tp->getTrackSpacing();
        trackNum = std::max(trackNum, 0);
        if (trackNum * (int) tp->getTrackSpacing() + tp->getStartCoord()
            < (tp->isHorizontal() ? bbox.xMin() : bbox.yMin())) {
          ++trackNum;
        }
        for (; trackNum < (int) tp->getNumTracks()
               && trackNum * (int) tp->getTrackSpacing() + tp->getStartCoord()
                      < (tp->isHorizontal() ? bbox.xMax() : bbox.yMax());
             ++trackNum) {
          frCoord trackLoc
              = trackNum * tp->getTrackSpacing() + tp->getStartCoord();
          if (tp->isHorizontal()) {
            horLoc2TrackPatterns[currLayerNum][trackLoc] = tp.get();
          } else {
            vertLoc2TrackPatterns[currLayerNum][trackLoc] = tp.get();
          }
        }
      }
    }
    layerNum2PreRouteDir[currLayerNum] = currPrefRouteDir;
  }
}

void FlexGridGraph::resetStatus()
{
  resetSrc();
  resetDst();
  resetPrevNodeDir();
}

void FlexGridGraph::resetSrc()
{
  srcs_.assign(srcs_.size(), false);
}

void FlexGridGraph::resetDst()
{
  dsts_.assign(dsts_.size(), false);
}

void FlexGridGraph::resetPrevNodeDir()
{
  prevDirs_.assign(prevDirs_.size(), false);
}

// print the grid graph with edge and vertex for debug purpose
void FlexGridGraph::print() const
{
  std::ofstream mazeLog(router_cfg_->OUT_MAZE_FILE.c_str());
  if (mazeLog.is_open()) {
    // print edges
    odb::Rect gridBBox;
    getBBox(gridBBox);
    mazeLog << "printing Maze grid (" << gridBBox.xMin() << ", "
            << gridBBox.yMin() << ") -- (" << gridBBox.xMax() << ", "
            << gridBBox.yMax() << ")\n";
    frMIdx xDim, yDim, zDim;
    getDim(xDim, yDim, zDim);

    if (xDim == 0 || yDim == 0 || zDim == 0) {
      std::cout << "Error: dimension == 0\n";
      return;
    }
    std::cout << "extBBox (xDim, yDim, zDim) = (" << xDim << ", " << yDim
              << ", " << zDim << ")\n";

    odb::Point p;
    for (frMIdx xIdx = 0; xIdx < xDim; ++xIdx) {
      for (frMIdx yIdx = 0; yIdx < yDim; ++yIdx) {
        for (frMIdx zIdx = 0; zIdx < zDim; ++zIdx) {
          if (hasEdge(xIdx, yIdx, zIdx, frDirEnum::N)) {
            if (yIdx + 1 >= yDim) {
              std::cout << "Error: no edge (" << xIdx << ", " << yIdx << ", "
                        << zIdx << ", N) " << yDim << '\n';
              continue;
            }
            mazeLog << "Edge: " << getPoint(p, xIdx, yIdx).x() << " "
                    << getPoint(p, xIdx, yIdx).y() << " " << zIdx << " "
                    << getPoint(p, xIdx, yIdx + 1).x() << " "
                    << getPoint(p, xIdx, yIdx + 1).y() << " " << zIdx << "\n";
          }
          if (hasEdge(xIdx, yIdx, zIdx, frDirEnum::E)) {
            if (xIdx + 1 >= xDim) {
              std::cout << "Error: no edge (" << xIdx << ", " << yIdx << ", "
                        << zIdx << ", E) " << xDim << '\n';
              continue;
            }
            mazeLog << "Edge: " << getPoint(p, xIdx, yIdx).x() << " "
                    << getPoint(p, xIdx, yIdx).y() << " " << zIdx << " "
                    << getPoint(p, xIdx + 1, yIdx).x() << " "
                    << getPoint(p, xIdx + 1, yIdx).y() << " " << zIdx << "\n";
          }
          if (hasEdge(xIdx, yIdx, zIdx, frDirEnum::U)) {
            if (zIdx + 1 >= zDim) {
              std::cout << "Error: no edge (" << xIdx << ", " << yIdx << ", "
                        << zIdx << ", U) " << zDim << '\n';
              continue;
            }
            mazeLog << "Edge: " << getPoint(p, xIdx, yIdx).x() << " "
                    << getPoint(p, xIdx, yIdx).y() << " " << zIdx << " "
                    << getPoint(p, xIdx, yIdx).x() << " "
                    << getPoint(p, xIdx, yIdx).y() << " " << zIdx + 1 << "\n";
          }
        }
      }
    }
  } else {
    std::cout << "Error: Fail to open maze log\n";
  }
}

}  // namespace drt
