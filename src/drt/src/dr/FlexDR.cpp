// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "dr/FlexDR.h"

#include <sys/stat.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"
#include "boost/io/ios_state.hpp"
#include "boost/polygon/polygon.hpp"
#include "db/infra/KDTree.hpp"
#include "db/infra/frTime.h"
#include "db/obj/frBlockObject.h"
#include "db/obj/frMarker.h"
#include "db/obj/frShape.h"
#include "db/obj/frVia.h"
#include "distributed/RoutingJobDescription.h"
#include "distributed/drUpdate.h"
#include "distributed/frArchive.h"
#include "dr/AbstractDRGraphics.h"
#include "dr/FlexDR_conn.h"
#include "drt/TritonRoute.h"
#include "dst/BalancerJobDescription.h"
#include "dst/Distributed.h"
#include "dst/JobMessage.h"
#include "frBaseTypes.h"
#include "frDesign.h"
#include "frProfileTask.h"
#include "gc/FlexGC.h"
#include "io/io.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"
#include "omp.h"
#include "serialization.h"
#include "utl/Logger.h"
#include "utl/Progress.h"
#include "utl/ScopedTemporaryFile.h"
#include "utl/exception.h"

using odb::dbTechLayerType;

BOOST_CLASS_EXPORT(drt::RoutingJobDescription)

namespace drt {

using utl::ThreadException;

enum class SerializationType
{
  READ,
  WRITE
};

void serializeWorker(FlexDRWorker* worker, std::string& workerStr)
{
  std::stringstream stream(std::ios_base::binary | std::ios_base::in
                           | std::ios_base::out);
  frOArchive ar(stream);
  registerTypes(ar);
  ar << *worker;
  workerStr = stream.str();
}

void deserializeWorker(FlexDRWorker* worker,
                       frDesign* design,
                       const std::string& workerStr)
{
  std::stringstream stream(
      workerStr,
      std::ios_base::binary | std::ios_base::in | std::ios_base::out);
  frIArchive ar(stream);
  ar.setDesign(design);
  registerTypes(ar);
  ar >> *worker;
}

void serializeViaData(const FlexDRViaData& viaData, std::string& serializedStr)
{
  std::stringstream stream(std::ios_base::binary | std::ios_base::in
                           | std::ios_base::out);
  frOArchive ar(stream);
  registerTypes(ar);
  ar << viaData;
  serializedStr = stream.str();
}

FlexDR::FlexDR(TritonRoute* router,
               frDesign* designIn,
               utl::Logger* loggerIn,
               odb::dbDatabase* dbIn,
               RouterConfiguration* router_cfg)
    : flow_state_machine_(std::make_unique<FlexDRFlow>()),
      router_(router),
      design_(designIn),
      logger_(loggerIn),
      db_(dbIn),
      router_cfg_(router_cfg),
      numWorkUnits_(0),
      dist_(nullptr),
      dist_on_(false),
      dist_port_(0),
      increaseClipsize_(false),
      clipSizeInc_(0),
      iter_(0)
{
}

FlexDR::~FlexDR() = default;

void FlexDR::setDebug(std::unique_ptr<AbstractDRGraphics> dr_graphics)
{
  graphics_ = std::move(dr_graphics);
}

std::string FlexDRWorker::reloadedMain()
{
  init(design_);
  debugPrint(logger_,
             utl::DRT,
             "autotuner",
             1,
             "Init number of markers {}",
             getInitNumMarkers());
  if (!skipRouting_) {
    route_queue();
  }
  setGCWorker(nullptr);
  cleanup();
  std::string workerStr;
  serializeWorker(this, workerStr);
  return workerStr;
}

void FlexDRWorker::writeUpdates(const std::string& file_name)
{
  std::vector<std::vector<drUpdate>> updates(1);
  for (const auto& net : getDesign()->getTopBlock()->getNets()) {
    for (const auto& guide : net->getGuides()) {
      for (const auto& connFig : guide->getRoutes()) {
        drUpdate update(drUpdate::ADD_GUIDE);
        frPathSeg* pathSeg = static_cast<frPathSeg*>(connFig.get());
        update.setPathSeg(*pathSeg);
        update.setIndexInOwner(guide->getIndexInOwner());
        update.setNet(net.get());
        updates.back().push_back(update);
      }
    }
    for (const auto& shape : net->getShapes()) {
      drUpdate update;
      auto pathSeg = static_cast<frPathSeg*>(shape.get());
      update.setPathSeg(*pathSeg);
      update.setNet(net.get());
      update.setUpdateType(drUpdate::ADD_SHAPE);
      updates.back().push_back(update);
    }
    for (const auto& shape : net->getVias()) {
      drUpdate update;
      auto via = static_cast<frVia*>(shape.get());
      update.setVia(*via);
      update.setNet(net.get());
      update.setUpdateType(drUpdate::ADD_SHAPE);
      updates.back().push_back(update);
    }
    for (const auto& shape : net->getPatchWires()) {
      drUpdate update;
      auto patch = static_cast<frPatchWire*>(shape.get());
      update.setPatchWire(*patch);
      update.setNet(net.get());
      update.setUpdateType(drUpdate::ADD_SHAPE);
      updates.back().push_back(update);
    }
  }
  for (const auto& marker : getDesign()->getTopBlock()->getMarkers()) {
    drUpdate update;
    update.setMarker(*marker);
    update.setUpdateType(drUpdate::ADD_SHAPE);
    updates.back().push_back(update);
  }
  std::ofstream file(file_name.c_str());
  frOArchive ar(file);
  registerTypes(ar);
  ar << updates;
  file.close();
}

int FlexDRWorker::main(frDesign* design)
{
  mainRead(design);
  return mainRoute();
}

void FlexDRWorker::mainRead(frDesign* design)
{
  ProfileTask profile("DRW:main");
  using std::chrono::high_resolution_clock;
  main_start_time_ = high_resolution_clock::now();
  auto micronPerDBU = 1.0 / getTech()->getDBUPerUU();
  if (router_cfg_->VERBOSE > 1) {
    logger_->report("start DR worker (BOX) ( {} {} ) ( {} {} )",
                    routeBox_.xMin() * micronPerDBU,
                    routeBox_.yMin() * micronPerDBU,
                    routeBox_.xMax() * micronPerDBU,
                    routeBox_.yMax() * micronPerDBU);
  }
  initMarkers(design);
  if (getDRIter() && getInitNumMarkers() == 0 && !needRecheck_) {
    skipRouting_ = true;
  }
  if (debugSettings_->debugDumpDR
      && (debugSettings_->box == odb::Rect(-1, -1, -1, -1)
          || routeBox_.intersects(debugSettings_->box))
      && !skipRouting_
      && (debugSettings_->iter == getDRIter()
          || debugSettings_->dumpLastWorker)) {
    std::string workerPath = fmt::format("{}/workerx{}_y{}",
                                         debugSettings_->dumpDir,
                                         routeBox_.xMin(),
                                         routeBox_.yMin());
    if (debugSettings_->dumpLastWorker) {
      workerPath = fmt::format("{}/workerx{}_y{}",
                               debugSettings_->dumpDir,
                               debugSettings_->box.xMin(),
                               debugSettings_->box.yMin());
    }
    if (mkdir(workerPath.c_str(), 0777) != 0) {
      logger_->error(
          DRT, 152, "Directory {} could not be created.", workerPath);
    }

    writeUpdates(fmt::format("{}/updates.bin", workerPath));
    {
      std::string viaDataStr;
      serializeViaData(*via_data_, viaDataStr);
      std::ofstream viaDataFile(
          fmt::format("{}/viadata.bin", workerPath).c_str());
      viaDataFile << viaDataStr;
      viaDataFile.close();
    }
    {
      std::string workerStr;
      serializeWorker(this, workerStr);
      std::ofstream workerFile(
          fmt::format("{}/worker.bin", workerPath).c_str());
      workerFile << workerStr;
      workerFile.close();
    }
    {
      std::ofstream router_cfgFile(
          fmt::format("{}/worker_router_cfg.bin", workerPath).c_str());
      frOArchive ar(router_cfgFile);
      registerTypes(ar);
      serializeGlobals(ar, router_cfg_);
      router_cfgFile.close();
    }
  }
  if (!skipRouting_) {
    readDesignNetObjects(design);
  }
}

int FlexDRWorker::mainRoute()
{
  ProfileTask profile("DRW:main");
  using std::chrono::high_resolution_clock;
  const high_resolution_clock::time_point t0 = main_start_time_;
  auto micronPerDBU = 1.0 / getTech()->getDBUPerUU();
  if (!skipRouting_) {
    init(getDesign());
  }
  const high_resolution_clock::time_point t1 = high_resolution_clock::now();
  if (!skipRouting_) {
    route_queue();
  }
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  const int num_markers = getNumMarkers();
  work_statistics_.init_seconds
      = std::chrono::duration<double>(t1 - t0).count();
  work_statistics_.total_seconds
      = std::chrono::duration<double>(t2 - t0).count();
  if (!skipRouting_
      && logger_->debugCheck(DRT, "bonnroute_cleanup_worker", 1)) {
    // Before cleanup(), which frees the pins of the worker's nets.
    reportWorkStatistics(num_markers);
  }
  // end() reads the worker's markers, nets and patches, not the GC state.
  // Releasing the GC worker here runs its destruction inside the parallel
  // loop of processWorkersBatch instead of in the serial endWorkersBatch,
  // and only the GC workers still running stay in memory together.
  setGCWorker(nullptr);
  cleanup();
  high_resolution_clock::time_point t3 = high_resolution_clock::now();
  work_statistics_.total_seconds
      = std::chrono::duration<double>(t3 - t0).count();

  using std::chrono::duration;
  using std::chrono::duration_cast;
  duration<double> time_span0 = duration_cast<duration<double>>(t1 - t0);
  duration<double> time_span1 = duration_cast<duration<double>>(t2 - t1);
  duration<double> time_span2 = duration_cast<duration<double>>(t3 - t2);

  if (router_cfg_->VERBOSE > 1) {
    std::stringstream ss;
    ss << "time (INIT/ROUTE/POST) " << time_span0.count() << " "
       << time_span1.count() << " " << time_span2.count() << " " << std::endl;
    std::cout << ss.str() << std::flush;
  }

  debugPrint(logger_,
             DRT,
             "autotuner",
             1,
             "worker ({:.3f} {:.3f}) ({:.3f} {:.3f}) time {} prev_#DRVs {} "
             "curr_#DRVs {}",
             routeBox_.xMin() * micronPerDBU,
             routeBox_.yMin() * micronPerDBU,
             routeBox_.xMax() * micronPerDBU,
             routeBox_.yMax() * micronPerDBU,
             duration_cast<duration<double>>(t3 - t0).count(),
             getInitNumMarkers(),
             num_markers);

  return 0;
}

void FlexDRWorker::reportWorkStatistics(const int num_markers) const
{
  // One line per worker that ran: what it was asked to repair, how much its
  // route queue did and where the seconds went. The nets with the most
  // reroutes follow as name:reroutes:pins in the worker:flags, with flags
  // c = clock net, n = non-default rule, l = left by the initial router
  // (no initial routing). `stopped 1` marks a worker whose route queue was
  // ended early (stubborn tiles after a clean worker, a discarded grid
  // worker started ahead); its numbers cover the part that ran.
  auto* block = design_->getTopBlock();
  const odb::Point gcell_min = block->getGCellIdx(routeBox_.ll());
  const odb::Point gcell_max = block->getGCellIdx(
      odb::Point(routeBox_.xMax() - 1, routeBox_.yMax() - 1));
  frMIdx x_dim = 0, y_dim = 0, z_dim = 0;
  if (!check_only_) {
    gridGraph_.getDim(x_dim, y_dim, z_dim);
  }
  std::vector<const drNet*> rerouted;
  int max_reroutes = 0;
  for (const auto& net : nets_) {
    if (net->getNumReroutes() > 0) {
      rerouted.push_back(net.get());
      max_reroutes = std::max(max_reroutes, net->getNumReroutes());
    }
  }
  std::ranges::stable_sort(rerouted, [](const drNet* lhs, const drNet* rhs) {
    return lhs->getNumReroutes() > rhs->getNumReroutes();
  });
  std::string top_nets;
  const size_t listed = std::min<size_t>(rerouted.size(), 6);
  for (size_t i = 0; i < listed; i++) {
    const drNet* net = rerouted[i];
    const frNet* fr_net = net->getFrNet();
    std::string flags;
    if (fr_net->isClock()) {
      flags += 'c';
    }
    if (fr_net->getNondefaultRule() != nullptr) {
      flags += 'n';
    }
    if (!fr_net->hasInitialRouting()) {
      flags += 'l';
    }
    top_nets += fmt::format("{}{}:{}:{}:{}",
                            i == 0 ? "" : " ",
                            fr_net->getName(),
                            net->getNumReroutes(),
                            net->getPins().size(),
                            flags.empty() ? "-" : flags);
  }
  logger_->debug(
      DRT,
      "bonnroute_cleanup_worker",
      "worker iter {} id {} box ( {} {} {} {} ) gcells ( {} {} {} {} ) mode {} "
      "maze_end_iter {} drc_cost {} marker_cost {} fixed_cost {} nets {} grid "
      "{}x{}x{} markers {} -> {} recheck {} routed {} checked {} dropped {} "
      "rerouted_nets {} max_reroutes {} init_s {:.3f} queue_s {:.3f} maze_s "
      "{:.3f} check_s {:.3f} route_end_s {:.3f} top_nets [{}] stopped {}",
      drIter_,
      worker_id_,
      routeBox_.xMin(),
      routeBox_.yMin(),
      routeBox_.xMax(),
      routeBox_.yMax(),
      gcell_min.x(),
      gcell_min.y(),
      gcell_max.x(),
      gcell_max.y(),
      static_cast<int>(ripupMode_),
      mazeEndIter_,
      workerDRCCost_,
      workerMarkerCost_,
      workerFixedShapeCost_,
      nets_.size(),
      x_dim,
      y_dim,
      z_dim,
      initNumMarkers_,
      num_markers,
      needRecheck_ ? 1 : 0,
      work_statistics_.routed_nets,
      work_statistics_.checked_objects,
      work_statistics_.dropped_entries,
      rerouted.size(),
      max_reroutes,
      work_statistics_.init_seconds,
      work_statistics_.queue_seconds,
      work_statistics_.maze_seconds,
      work_statistics_.check_seconds,
      work_statistics_.total_seconds,
      top_nets,
      stopped_ ? 1 : 0);
}

void FlexDRWorker::distributedMain(frDesign* design)
{
  ProfileTask profile("DR:main");
  if (router_cfg_->VERBOSE > 1) {
    logger_->report("start DR worker (BOX) ( {} {} ) ( {} {} )",
                    routeBox_.xMin() * 1.0 / getTech()->getDBUPerUU(),
                    routeBox_.yMin() * 1.0 / getTech()->getDBUPerUU(),
                    routeBox_.xMax() * 1.0 / getTech()->getDBUPerUU(),
                    routeBox_.yMax() * 1.0 / getTech()->getDBUPerUU());
  }
  initMarkers(design);
  if (getDRIter() && getInitNumMarkers() == 0 && !needRecheck_) {
    skipRouting_ = true;
    return;
  }
}

void FlexDR::initFromTA()
{
  // initialize lists
  for (auto& net : getDesign()->getTopBlock()->getNets()) {
    for (auto& guide : net->getGuides()) {
      for (auto& connFig : guide->getRoutes()) {
        if (connFig->typeId() == frcPathSeg) {
          std::unique_ptr<frShape> ps = std::make_unique<frPathSeg>(
              *(static_cast<frPathSeg*>(connFig.get())));
          auto [bp, ep] = static_cast<frPathSeg*>(ps.get())->getPoints();

          // skip TA dummy segment
          if (odb::Point::manhattanDistance(ep, bp) != 1) {
            net->addShape(std::move(ps));
          }
        } else {
          std::cout << "Error: initFromTA unsupported shape" << std::endl;
        }
      }
    }
  }
}

void FlexDR::initGCell2BoundaryPin()
{
  // initialize size
  auto topBlock = getDesign()->getTopBlock();
  auto gCellPatterns = topBlock->getGCellPatterns();
  auto& xgp = gCellPatterns.at(0);
  auto& ygp = gCellPatterns.at(1);
  auto tmpVec = std::vector<
      frOrderedIdMap<frNet*, std::set<std::pair<odb::Point, frLayerNum>>>>(
      (int) ygp.getCount());
  gcell2BoundaryPin_ = std::vector<std::vector<
      frOrderedIdMap<frNet*, std::set<std::pair<odb::Point, frLayerNum>>>>>(
      (int) xgp.getCount(), tmpVec);
  for (auto& net : topBlock->getNets()) {
    auto netPtr = net.get();
    for (auto& guide : net->getGuides()) {
      for (auto& connFig : guide->getRoutes()) {
        if (connFig->typeId() == frcPathSeg) {
          auto ps = static_cast<frPathSeg*>(connFig.get());
          auto [bp, ep] = ps->getPoints();
          frLayerNum layerNum = ps->getLayerNum();
          // skip TA dummy segment
          auto mdist = odb::Point::manhattanDistance(ep, bp);
          if (mdist == 1 || mdist == 0) {
            continue;
          }
          odb::Point idx1 = design_->getTopBlock()->getGCellIdx(bp);
          odb::Point idx2 = design_->getTopBlock()->getGCellIdx(ep);

          // update gcell2BoundaryPin
          // horizontal
          if (bp.y() == ep.y()) {
            int x1 = idx1.x();
            int x2 = idx2.x();
            int y = idx1.y();
            for (auto x = x1; x <= x2; ++x) {
              odb::Rect gcellBox = topBlock->getGCellBox(odb::Point(x, y));
              frCoord leftBound = gcellBox.xMin();
              frCoord rightBound = gcellBox.xMax();
              const bool hasLeftBound = bp.x() < leftBound;
              const bool hasRightBound = ep.x() >= rightBound;
              if (hasLeftBound) {
                odb::Point boundaryPt(leftBound, bp.y());
                gcell2BoundaryPin_[x][y][netPtr].emplace(boundaryPt, layerNum);
              }
              if (hasRightBound) {
                odb::Point boundaryPt(rightBound, ep.y());
                gcell2BoundaryPin_[x][y][netPtr].emplace(boundaryPt, layerNum);
              }
            }
          } else if (bp.x() == ep.x()) {
            int x = idx1.x();
            int y1 = idx1.y();
            int y2 = idx2.y();
            for (auto y = y1; y <= y2; ++y) {
              odb::Rect gcellBox = topBlock->getGCellBox(odb::Point(x, y));
              frCoord bottomBound = gcellBox.yMin();
              frCoord topBound = gcellBox.yMax();
              const bool hasBottomBound = bp.y() < bottomBound;
              const bool hasTopBound = ep.y() >= topBound;
              if (hasBottomBound) {
                odb::Point boundaryPt(bp.x(), bottomBound);
                gcell2BoundaryPin_[x][y][netPtr].emplace(boundaryPt, layerNum);
              }
              if (hasTopBound) {
                odb::Point boundaryPt(ep.x(), topBound);
                gcell2BoundaryPin_[x][y][netPtr].emplace(boundaryPt, layerNum);
              }
            }
          } else {
            std::cout
                << "Error: non-orthogonal pathseg in initGCell2BoundryPin\n";
          }
        }
      }
    }
  }
}

void FlexDR::init_halfViaEncArea()
{
  auto bottomLayerNum = getTech()->getBottomLayerNum();
  auto topLayerNum = getTech()->getTopLayerNum();
  auto& halfViaEncArea = via_data_.halfViaEncArea;
  for (int i = bottomLayerNum; i <= topLayerNum; i++) {
    if (getTech()->getLayer(i)->getType() != dbTechLayerType::ROUTING) {
      continue;
    }
    if (i + 1 <= topLayerNum
        && getTech()->getLayer(i + 1)->getType() == dbTechLayerType::CUT) {
      auto viaDef = getTech()->getLayer(i + 1)->getDefaultViaDef();
      if (viaDef) {
        frVia via(viaDef);
        odb::Rect layer1Box = via.getLayer1BBox();
        odb::Rect layer2Box = via.getLayer2BBox();
        auto layer1HalfArea = layer1Box.area() / 2;
        auto layer2HalfArea = layer2Box.area() / 2;
        halfViaEncArea.emplace_back(layer1HalfArea, layer2HalfArea);
      } else {
        halfViaEncArea.emplace_back(0, 0);
      }
    } else {
      halfViaEncArea.emplace_back(0, 0);
    }
  }
}

void FlexDR::init()
{
  ProfileTask profile("DR:init");
  frTime t;
  if (router_cfg_->VERBOSE > 0) {
    logger_->info(DRT, 187, "Start routing data preparation.");
  }
  initGCell2BoundaryPin();
  getRegionQuery()->initDRObj();  // first init in postProcess

  init_halfViaEncArea();

  if (router_cfg_->VERBOSE > 0) {
    t.print(logger_);
  }

  iter_ = 0;

  if (router_cfg_->VERBOSE > 0) {
    logger_->info(DRT, 194, "Start detail routing.");
  }
  for (const auto& net : getDesign()->getTopBlock()->getNets()) {
    if (!net->hasInitialRouting()) {
      continue;
    }
    for (const auto& via : net->getVias()) {
      const odb::Rect bottomBox = via->getLayer1BBox();
      const odb::Rect topBox = via->getLayer2BBox();
      // A via can land on an instance pin or on an IO pin.
      auto mark_via_connected = [&via, &bottomBox, &topBox](auto* term) {
        if (term->getBBox().intersects(bottomBox)) {
          std::vector<frRect> shapes;
          term->getShapes(shapes);
          for (const auto& shape : shapes) {
            if (shape.getLayerNum() != via->getViaDef()->getLayer1Num()) {
              continue;
            }
            if (!shape.intersects(bottomBox)) {
              continue;
            }
            via->setBottomConnected(true);
            break;
          }
        }
        if (via->isBottomConnected()) {
          return;
        }
        if (term->getBBox().intersects(topBox)) {
          std::vector<frRect> shapes;
          term->getShapes(shapes);
          for (const auto& shape : shapes) {
            if (shape.getLayerNum() != via->getViaDef()->getLayer2Num()) {
              continue;
            }
            if (!shape.intersects(topBox)) {
              continue;
            }
            via->setTopConnected(true);
            break;
          }
        }
      };
      for (auto term : net->getInstTerms()) {
        mark_via_connected(term);
      }
      for (auto term : net->getBTerms()) {
        mark_via_connected(term);
      }
    }
  }
}

void FlexDR::removeGCell2BoundaryPin()
{
  gcell2BoundaryPin_.clear();
  gcell2BoundaryPin_.shrink_to_fit();
}

frOrderedIdMap<frNet*, std::set<std::pair<odb::Point, frLayerNum>>>
FlexDR::initDR_mergeBoundaryPin(int startX,
                                int startY,
                                int size,
                                const odb::Rect& routeBox) const
{
  frOrderedIdMap<frNet*, std::set<std::pair<odb::Point, frLayerNum>>> bp;
  auto gCellPatterns = getDesign()->getTopBlock()->getGCellPatterns();
  auto& xgp = gCellPatterns.at(0);
  auto& ygp = gCellPatterns.at(1);
  for (int i = startX; i < (int) xgp.getCount() && i < startX + size; i++) {
    for (int j = startY; j < (int) ygp.getCount() && j < startY + size; j++) {
      auto& currBp = gcell2BoundaryPin_[i][j];
      for (auto& [net, s] : currBp) {
        for (auto& [pt, lNum] : s) {
          if (pt.x() == routeBox.xMin() || pt.x() == routeBox.xMax()
              || pt.y() == routeBox.yMin() || pt.y() == routeBox.yMax()) {
            bp[net].emplace(pt, lNum);
          }
        }
      }
    }
  }
  return bp;
}

void FlexDR::getBatchInfo(int& batchStepX, int& batchStepY)
{
  batchStepX = 2;
  batchStepY = 2;
}

std::unique_ptr<FlexDRWorker> FlexDR::createWorker(const int x_offset,
                                                   const int y_offset,
                                                   const SearchRepairArgs& args,
                                                   const odb::Rect& routeBoxIn)
{
  auto worker = std::make_unique<FlexDRWorker>(
      &via_data_, getDesign(), logger_, router_cfg_);
  odb::Rect route_box(routeBoxIn);
  if (route_box == odb::Rect(0, 0, 0, 0)) {
    auto gCellPatterns = getDesign()->getTopBlock()->getGCellPatterns();
    auto& xgp = gCellPatterns.at(0);
    auto& ygp = gCellPatterns.at(1);
    odb::Rect routeBox1 = getDesign()->getTopBlock()->getGCellBox(
        odb::Point(x_offset, y_offset));
    const int max_i
        = std::min((int) xgp.getCount() - 1, x_offset + args.size - 1);
    const int max_j = std::min((int) ygp.getCount(), y_offset + args.size - 1);
    odb::Rect routeBox2
        = getDesign()->getTopBlock()->getGCellBox(odb::Point(max_i, max_j));
    route_box.init(
        routeBox1.xMin(), routeBox1.yMin(), routeBox2.xMax(), routeBox2.yMax());
  }
  odb::Rect extBox;
  odb::Rect drcBox;
  route_box.bloat(router_cfg_->MTSAFEDIST, extBox);
  route_box.bloat(router_cfg_->DRCSAFEDIST, drcBox);
  worker->setRouteBox(route_box);
  worker->setExtBox(extBox);
  worker->setDrcBox(drcBox);
  worker->setMazeEndIter(args.mazeEndIter);
  worker->setDRIter(iter_);
  worker->setDebugSettings(router_->getDebugSettings());
  if (dist_on_) {
    worker->setDistributed(dist_, dist_ip_, dist_port_, dist_dir_);
  }
  if (!iter_) {
    auto bp = initDR_mergeBoundaryPin(x_offset, y_offset, args.size, route_box);
    worker->setDRIter(0, bp);
  }
  worker->setRipupMode(args.ripupMode);
  // Only iterations 0-2 (incremental because of initial routing): later
  // incremental iterations of designs with fixed nets also gauge worker
  // congestion on the grid graph.
  worker->setCheckOnlyWithoutReroutableNets(
      initial_routing_cleanup_.check_only_workers
      && args.ripupMode == RipUpMode::INCR && iter_ <= 2 && !dist_on_);
  worker->setFollowGuide(args.followGuide);
  // TODO: only pass to relevant workers
  worker->setGraphics(graphics_.get());
  worker->setCost(args.workerDRCCost,
                  args.workerMarkerCost,
                  args.workerFixedShapeCost,
                  args.workerMarkerDecay);
  return worker;
}

namespace {
void printIteration(utl::Logger* logger,
                    const int iter,
                    const std::string& flow_name)
{
  std::string suffix;
  if (iter == 1 || (iter > 20 && iter % 10 == 1)) {
    suffix = "st";
  } else if (iter == 2 || (iter > 20 && iter % 10 == 2)) {
    suffix = "nd";
  } else if (iter == 3 || (iter > 20 && iter % 10 == 3)) {
    suffix = "rd";
  } else {
    suffix = "th";
  }
  logger->info(DRT, 195, "Start {}{} {} iteration.", iter, suffix, flow_name);
}

void printIterationProgress(utl::Logger* logger,
                            FlexDR::IterationProgress& iter_prog,
                            const int num_markers,
                            const int max_perc = 90)
{
  int progress
      = (++iter_prog.cnt_done_workers * 100) / iter_prog.total_num_workers;
  progress = (progress / 10) * 10;  // clip it to multiples of 10%
  if (progress > iter_prog.last_reported_perc && progress <= max_perc) {
    iter_prog.last_reported_perc = progress;
    logger->report(
        "    Completing {}% with {} violations.", progress, num_markers);
    logger->report("    {}.", iter_prog.time);
  }
}
}  // namespace

void FlexDR::reportIterationViolations() const
{
  if (router_cfg_->VERBOSE > 0) {
    logger_->info(DRT,
                  199,
                  "  Number of violations = {}.",
                  getDesign()->getTopBlock()->getNumMarkers());
    if (getDesign()->getTopBlock()->getNumMarkers() > 0) {
      // report violations
      std::map<std::string, std::map<frLayerNum, uint32_t>> violations;
      std::set<frLayerNum> layers;
      const std::map<std::string, std::string> relabel
          = {{"Lef58SpacingEndOfLine", "EOL"},
             {"Lef58CutSpacingTable", "CutSpcTbl"},
             {"Lef58EolKeepOut", "eolKeepOut"}};
      for (const auto& marker : getDesign()->getTopBlock()->getMarkers()) {
        if (!marker->getConstraint()) {
          continue;
        }
        auto type = marker->getConstraint()->getViolName();
        if (relabel.find(type) != relabel.end()) {
          type = relabel.at(type);
        }
        violations[type][marker->getLayerNum()]++;
        layers.insert(marker->getLayerNum());
      }
      std::string line = fmt::format("{:<15}", "Viol/Layer");
      for (auto lNum : layers) {
        std::string lName = getTech()->getLayer(lNum)->getName();
        if (lName.size() >= 7) {
          lName = lName.substr(0, 2) + ".." + lName.substr(lName.size() - 2, 2);
        }
        line += fmt::format("{:>7}", lName);
      }
      logger_->report(line);
      for (auto& [type, typeViolations] : violations) {
        std::string typeName = type;
        if (typeName.size() >= 15) {
          typeName = typeName.substr(0, 12) + "..";
        }
        line = fmt::format("{:<15}", typeName);
        for (auto lNum : layers) {
          line += fmt::format("{:>7}", typeViolations[lNum]);
        }
        logger_->report(line);
      }
    }
  }
  if ((router_cfg_->DRC_RPT_ITER_STEP > 0 && iter_ > 0
       && iter_ % router_cfg_->DRC_RPT_ITER_STEP == 0)
      || logger_->debugCheck(DRT, "autotuner", 1)
      || logger_->debugCheck(DRT, "report", 1)) {
    router_->reportDRC(
        router_cfg_->DRC_RPT_FILE + '-' + std::to_string(iter_) + ".rpt",
        design_->getTopBlock()->getMarkers(),
        "DRC - iter " + std::to_string(iter_));
  }
}

void FlexDR::processWorkersBatch(
    std::vector<std::unique_ptr<FlexDRWorker>>& workers_batch,
    IterationProgress& iter_prog)
{
  const int num_markers = getDesign()->getTopBlock()->getNumMarkers();
  ThreadException exception;
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < (int) workers_batch.size(); i++) {  // NOLINT
    try {
      workers_batch[i]->main(getDesign());
#pragma omp critical
      {
        if (router_cfg_->VERBOSE > 0) {
          printIterationProgress(logger_, iter_prog, num_markers);
        }
      }
    } catch (...) {
      exception.capture();
    }
  }
  exception.rethrow();
}

std::vector<int> FlexDR::runStubbornWorkersUntilClean(
    std::vector<std::unique_ptr<FlexDRWorker>>& workers_batch,
    IterationProgress& iter_prog)
{
  const int num_markers = getDesign()->getTopBlock()->getNumMarkers();
  // Per cluster: the lowest creation index of a worker that ran to the end
  // with zero markers.  Workers at or below that index are never skipped or
  // stopped, so the first zero-marker worker always runs to the end.
  std::map<int, std::unique_ptr<std::atomic<int>>> first_clean;
  for (int i = 0; i < (int) workers_batch.size(); i++) {
    const int cluster = workers_batch[i]->getWorkerId();
    auto it = first_clean.find(cluster);
    if (it == first_clean.end()) {
      it = first_clean
               .emplace(cluster,
                        std::make_unique<std::atomic<int>>(
                            std::numeric_limits<int>::max()))
               .first;
    }
    workers_batch[i]->setFirstCleanCreationIndex(it->second.get(), i);
  }
  std::vector<char> finished(workers_batch.size(), 0);
  std::atomic<int> not_finished{0};
  ThreadException exception;
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < (int) workers_batch.size(); i++) {  // NOLINT
    try {
      auto& worker = workers_batch[i];
      std::atomic<int>& clean = *first_clean.at(worker->getWorkerId());
      if (clean.load() < i) {
        not_finished++;
        continue;
      }
      worker->main(getDesign());
      if (worker->isStopped()) {
        not_finished++;
        continue;
      }
      finished[i] = 1;
      if (worker->getBestNumMarkers() == 0) {
        int current = clean.load();
        while (i < current && !clean.compare_exchange_weak(current, i)) {
        }
      }
#pragma omp critical
      {
        if (router_cfg_->VERBOSE > 0) {
          printIterationProgress(logger_, iter_prog, num_markers);
        }
      }
    } catch (...) {
      exception.capture();
    }
  }
  // The shared indices live in this function; the workers outlive it.
  for (auto& worker : workers_batch) {
    worker->setFirstCleanCreationIndex(nullptr, 0);
  }
  exception.rethrow();
  iter_prog.total_num_workers -= not_finished.load();
  std::vector<int> ran;
  for (int i = 0; i < (int) workers_batch.size(); i++) {
    if (finished[i]) {
      ran.push_back(i);
    }
  }
  return ran;
}

void FlexDR::processWorkersBatchDistributed(
    std::vector<std::unique_ptr<FlexDRWorker>>& workers_batch,
    int& version,
    IterationProgress& iter_prog)
{
  if (router_->dist_pool_.has_value()) {
    router_->dist_pool_->join();
  }
  if (version++ == 0 && !design_->hasUpdates()) {
    std::string serializedViaData;
    serializeViaData(via_data_, serializedViaData);
    router_->sendGlobalsUpdates(router_cfg_path_, serializedViaData);
  } else {
    router_->sendDesignUpdates(router_cfg_path_, router_cfg_->MAX_THREADS);
  }

  ProfileTask task("DIST: PROCESS_BATCH");
  const int num_markers = getDesign()->getTopBlock()->getNumMarkers();
  // multi thread
  ThreadException exception;
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < (int) workers_batch.size(); i++) {  // NOLINT
    try {
      workers_batch[i]->distributedMain(getDesign());
    } catch (...) {
      exception.capture();
    }
  }
  exception.rethrow();
  int j = 0;
  std::vector<std::vector<std::pair<int, FlexDRWorker*>>> distWorkerBatches(
      router_->getCloudSize());
  for (int i = 0; i < workers_batch.size(); i++) {
    auto worker = workers_batch.at(i).get();
    if (!worker->isSkipRouting()) {
      distWorkerBatches[j].push_back({i, worker});
      j = (j + 1) % router_->getCloudSize();
    }
  }
  {
    ProfileTask task("DIST: SERIALIZE+SEND");
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < distWorkerBatches.size(); i++) {  // NOLINT
      sendWorkers(distWorkerBatches.at(i), workers_batch);
    }
  }
  std::vector<std::pair<int, std::string>> workers;
  router_->getWorkerResults(workers);
  {
    ProfileTask task("DIST: DESERIALIZING_BATCH");
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < workers.size(); i++) {  // NOLINT
      deserializeWorker(workers_batch.at(workers.at(i).first).get(),
                        design_,
                        workers.at(i).second);
#pragma omp critical
      {
        if (router_cfg_->VERBOSE > 0) {
          printIterationProgress(logger_, iter_prog, num_markers);
        }
      }
    }
  }
}

void FlexDR::normalizeInitialRouting()
{
  // The connectivity check processes the nets marked modified, which are the
  // nets the other router wrote; it clears the flag. Iteration 0 checks every
  // tile, so the Recheck markers of the removed and shrunk objects are
  // dropped.
  FlexDRConnectivityChecker checker(
      router_, logger_, router_cfg_, graphics_.get(), dist_on_);
  checker.check(0);
  std::vector<frMarker*> markers;
  for (const auto& marker : getDesign()->getTopBlock()->getMarkers()) {
    markers.push_back(marker.get());
  }
  for (frMarker* marker : markers) {
    getRegionQuery()->removeMarker(marker);
    getDesign()->getTopBlock()->removeMarker(marker);
  }
  if (router_cfg_->VERBOSE > 0) {
    logger_->info(DRT,
                  9150,
                  "Initial routing clean-up: connectivity normalisation "
                  "before iteration 0 removed or shrank {} objects.",
                  markers.size());
  }
}

void FlexDR::endWorker(FlexDRWorker& worker,
                       std::vector<odb::Rect>* changed_regions)
{
  cleanup_statistics_.workers++;
  cleanup_statistics_.workers_without_markers
      += worker.getInitNumMarkers() == 0 ? 1 : 0;
  cleanup_statistics_.workers_without_nets
      += worker.isSkippedWithoutNets() ? 1 : 0;
  cleanup_statistics_.workers_with_only_initial_routing
      += worker.hasOnlyNetsWithInitialRouting() ? 1 : 0;
  if (worker.end(getDesign())) {
    numWorkUnits_ += 1;
  }
  if (worker.isCheckOnly()) {
    num_check_only_workers_++;
  }
  if (changed_regions != nullptr && worker.hasChangedDesign()) {
    changed_regions->push_back(worker.getExtBox());
  }
  if (worker.isCongested()) {
    increaseClipsize_ = true;
  }
}

void FlexDR::endWorkersBatch(
    std::vector<std::unique_ptr<FlexDRWorker>>& workers_batch,
    std::vector<odb::Rect>* changed_regions)
{
  ProfileTask profile("DR:end_batch");
  // single thread
  for (auto& worker : workers_batch) {
    endWorker(*worker, changed_regions);
  }
  workers_batch.clear();
}

// The batch loop of optimizationFlow runs the batches one after the other:
// a batch starts when the slowest worker of the previous batch has ended.
// In the marker-driven iterations few workers run, so most threads idle
// behind one slow worker per batch.  This function keeps the results of that
// loop and fills the idle threads:
//  - A worker copies the changing design state in mainRead() and builds and
//    routes on its own data in mainRoute().  Every read and every write of
//    that design state happens under one mutex; mainRoute() runs outside it.
//  - Which workers of a batch run is decided when every earlier batch is
//    written back, on the design state the batch loop decides on: the
//    workers with a marker in their DRC box.  From iteration 2 on a worker
//    without one skips its routing and writes nothing (FlexDRWorker::main,
//    FlexDRWorker::end), so it needs no run.  The workers that run read the
//    design at that moment, before any result of their batch is written.
//  - Results are written back (end()) in the order of the batch loop: batch
//    by batch, creation order inside a batch; a result is written as soon as
//    every earlier result is written.
//  - A worker of a later batch may start ahead on an idle thread.  When its
//    batch starts, its design-read fingerprint is compared with the
//    fingerprint on the design state of that moment: equal means it read
//    exactly what it reads in the batch loop, so its result is kept; a
//    different one stops it and a fresh worker runs in the batch.
//  - When a result is written, the workers started ahead whose extended box
//    intersects the writer's are compared the same way at once; a worker
//    that read something else is stopped and may start ahead again.
// Start order (it does not enter the results): workers of the running batch,
// most markers first; then workers of later batches whose overlapping
// earlier workers are all written back; then the other workers of later
// batches.  The last group pays off when the overlapping earlier worker
// ends with more markers than it started with and writes nothing.
int FlexDR::runGridWorkersAhead(std::vector<GridWorkerSlot>& slots,
                                const int num_batches,
                                const SearchRepairArgs& args,
                                IterationProgress& iter_prog,
                                const bool skip_unaffected,
                                std::vector<odb::Rect>* changed_regions)
{
  enum class State
  {
    kWaiting,  // has not read the design
    kRead,     // read the design at the start of its batch
    kRunning,
    kDone,
    kEnded  // written back, or without markers at the start of its batch
  };
  struct Slot
  {
    State state = State::kWaiting;
    uint64_t fingerprint = 0;  // of a worker started ahead
    bool runs_in_batch = false;
    std::unique_ptr<std::atomic<bool>> stop_request;
    // Slots whose extended box intersects this one.
    std::vector<int> overlaps;
    // Markers in the DRC box when last counted.  The count orders the starts
    // and selects the workers that start ahead; it does not enter the
    // results.
    int num_markers = 0;
  };
  const int num_slots = slots.size();
  if (num_slots == 0) {
    return 0;
  }
  std::vector<Slot> states(num_slots);
  std::vector<int> batch_end(num_batches, 0);
  for (int i = 0; i < num_slots; i++) {
    batch_end[slots[i].batch] = i + 1;
  }
  for (int batch = 1; batch < num_batches; batch++) {
    batch_end[batch] = std::max(batch_end[batch], batch_end[batch - 1]);
  }
  {
    // Extended boxes intersect only between tiles next to each other.
    std::map<std::pair<int, int>, int> slot_at;
    for (int i = 0; i < num_slots; i++) {
      slot_at[{slots[i].x_offset, slots[i].y_offset}] = i;
    }
    for (int i = 0; i < num_slots; i++) {
      for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
          const auto it = slot_at.find({slots[i].x_offset + dx * args.size,
                                        slots[i].y_offset + dy * args.size});
          if (it == slot_at.end() || it->second == i
              || !slots[it->second].worker->getExtBox().intersects(
                  slots[i].worker->getExtBox())) {
            continue;
          }
          states[i].overlaps.push_back(it->second);
        }
      }
    }
  }

  std::mutex mutex;
  std::condition_variable wake;
  // Everything below is read and written with the mutex held.
  int front = 0;  // slots before it are ended
  int current_batch = -1;
  bool finished = false;
  bool failed = false;
  int num_markers_at_batch_start = 0;  // for the progress lines
  // Waiting slots of later batches with markers in their DRC box.
  std::set<int> marked_later_slots;
  bool no_worker_ahead = false;  // until the next result is written back
  int num_unaffected_workers = 0;
  int num_started_ahead = 0;
  int num_started_before_neighbours = 0;
  int num_kept = 0;
  int num_rerun = 0;
  int num_restarted_ahead = 0;
  int num_dropped = 0;
  std::priority_queue<std::pair<int, int>> batch_queue;  // markers, -slot
  std::vector<std::unique_ptr<FlexDRWorker>> discarded_workers;
  std::vector<std::unique_ptr<std::atomic<bool>>> discarded_stop_requests;
  ThreadException exception;

  // Needs a slot that still holds its worker.
  const auto count_markers = [&](const int i) {
    std::vector<frMarker*> markers;
    getDesign()->getRegionQuery()->queryMarker(slots[i].worker->getDrcBox(),
                                               markers);
    states[i].num_markers = markers.size();  // a Recheck marker counts
    if (slots[i].batch > current_batch && states[i].state == State::kWaiting
        && states[i].num_markers > 0) {
      marked_later_slots.insert(i);
    } else {
      marked_later_slots.erase(i);
    }
  };
  const auto discard = [&](const int i) {
    Slot& state = states[i];
    if (state.state == State::kRunning) {
      state.stop_request->store(true);
    }
    discarded_workers.push_back(std::move(slots[i].worker));
    discarded_stop_requests.push_back(std::move(state.stop_request));
  };
  const auto report_progress = [&]() {
    if (router_cfg_->VERBOSE > 0) {
      printIterationProgress(logger_, iter_prog, num_markers_at_batch_start);
    }
  };
  // Needs: every earlier batch is ended.
  const auto start_batch = [&](const int batch) {
    current_batch = batch;
    batch_queue = {};
    no_worker_ahead = false;
    num_markers_at_batch_start = getDesign()->getTopBlock()->getNumMarkers();
    for (int i = front; i < batch_end[batch]; i++) {
      Slot& state = states[i];
      const bool started = state.state != State::kWaiting;
      count_markers(i);
      if (state.num_markers == 0) {
        num_unaffected_workers++;
        if (skip_unaffected) {
          iter_prog.total_num_workers--;
        } else {
          cleanup_statistics_.workers++;
          cleanup_statistics_.workers_without_markers++;
          report_progress();
        }
        if (started) {
          num_dropped++;
          discard(i);
        }
        slots[i].worker.reset();
        state.state = State::kEnded;
        continue;
      }
      state.runs_in_batch = true;
      if (started) {
        if (state.fingerprint
            == slots[i].worker->getDesignReadFingerprint(getDesign())) {
          num_kept++;
          if (state.state == State::kDone) {
            report_progress();
          }
          continue;
        }
        num_rerun++;
        discard(i);
        slots[i].worker
            = createWorker(slots[i].x_offset, slots[i].y_offset, args);
      }
      slots[i].worker->mainRead(getDesign());
      state.state = State::kRead;
      batch_queue.emplace(state.num_markers, -i);
    }
  };
  // Writes back what the order allows.  Call after every change of a slot.
  const auto write_back = [&]() {
    while (!finished) {
      bool wrote = false;
      while (front < batch_end[current_batch]
             && (states[front].state == State::kDone
                 || states[front].state == State::kEnded)) {
        if (states[front].state == State::kDone) {
          endWorker(*slots[front].worker, changed_regions);
          slots[front].worker.reset();
          states[front].state = State::kEnded;
          for (const int other : states[front].overlaps) {
            if (slots[other].batch <= current_batch) {
              continue;
            }
            Slot& ahead = states[other];
            if ((ahead.state == State::kRunning
                 || ahead.state == State::kDone)
                && ahead.fingerprint
                       != slots[other].worker->getDesignReadFingerprint(
                           getDesign())) {
              num_restarted_ahead++;
              discard(other);
              slots[other].worker = createWorker(
                  slots[other].x_offset, slots[other].y_offset, args);
              ahead.state = State::kWaiting;
            }
            count_markers(other);
          }
          wrote = true;
        }
        front++;
      }
      if (wrote) {
        no_worker_ahead = false;
      }
      if (front < batch_end[current_batch]) {
        return;
      }
      if (current_batch + 1 == num_batches) {
        finished = true;
      } else {
        start_batch(current_batch + 1);
      }
    }
  };
  // A slot of a later batch to start ahead, or -1.
  const auto pick_ahead = [&]() {
    if (no_worker_ahead) {
      return -1;
    }
    int before_neighbours = -1;
    for (const int i : marked_later_slots) {
      // An overlapping worker of an earlier batch that is still to write
      // its result back changes what this worker reads, unless it ends
      // with more markers than it started with.
      const bool neighbours_pending
          = std::ranges::any_of(states[i].overlaps, [&](const int other) {
              if (slots[other].batch >= slots[i].batch
                  || states[other].state == State::kEnded) {
                return false;
              }
              return slots[other].batch == current_batch
                         ? states[other].runs_in_batch
                         : states[other].num_markers > 0;
            });
      if (!neighbours_pending) {
        return i;
      }
      if (before_neighbours < 0) {
        before_neighbours = i;
      }
    }
    if (before_neighbours >= 0) {
      num_started_before_neighbours++;
      return before_neighbours;
    }
    no_worker_ahead = true;
    return -1;
  };

  for (int i = 0; i < num_slots; i++) {
    count_markers(i);
  }
  start_batch(0);
  write_back();
#pragma omp parallel
  {
    std::unique_lock<std::mutex> lock(mutex);
    while (!finished && !failed) {
      try {
        int picked = -1;
        while (!batch_queue.empty() && picked < 0) {
          const int i = -batch_queue.top().second;
          batch_queue.pop();
          if (states[i].state == State::kRead) {
            picked = i;
          }
        }
        if (picked < 0) {
          picked = pick_ahead();
          if (picked >= 0) {
            marked_later_slots.erase(picked);
            num_started_ahead++;
            states[picked].fingerprint
                = slots[picked].worker->getDesignReadFingerprint(getDesign());
            slots[picked].worker->mainRead(getDesign());
          }
        }
        if (picked < 0) {
          wake.wait(lock);
          continue;
        }
        Slot& state = states[picked];
        FlexDRWorker* worker = slots[picked].worker.get();
        state.state = State::kRunning;
        state.stop_request = std::make_unique<std::atomic<bool>>(false);
        worker->setStopRequest(state.stop_request.get());
        lock.unlock();
        worker->mainRoute();
        lock.lock();
        // A worker started ahead may have been replaced or dropped.
        if (slots[picked].worker.get() == worker) {
          state.state = State::kDone;
          if (slots[picked].batch == current_batch) {
            report_progress();
          }
          write_back();
        }
        wake.notify_all();
      } catch (...) {
        if (!lock.owns_lock()) {
          lock.lock();
        }
        exception.capture();
        failed = true;
        wake.notify_all();
      }
    }
  }
  exception.rethrow();
  debugPrint(logger_,
             DRT,
             "bonnroute_cleanup",
             1,
             "grid workers ahead: iteration {} slots {} batches {} started "
             "ahead {} kept {} run again {} dropped without markers {} "
             "started before neighbours {} restarted ahead {}",
             iter_,
             num_slots,
             num_batches,
             num_started_ahead,
             num_kept,
             num_rerun,
             num_dropped,
             num_started_before_neighbours,
             num_restarted_ahead);
  return num_unaffected_workers;
}

void FlexDR::reportClusterResults(
    const std::vector<std::unique_ptr<FlexDRWorker>>& batch,
    const std::map<int, FlexDRWorker*>& chosen) const
{
  // One line per cluster of a guides-tiles or stubborn-tiles batch: the
  // workers that competed for it in creation order as
  // markers after:seconds:drc cost:marker cost:route box:state, and the one
  // whose result is written back. State: r = ran to the end, x = stopped
  // after an earlier worker of the cluster ended without markers, s = not
  // started for the same reason (stubborn tiles); only r workers compete.
  std::map<int, std::vector<const FlexDRWorker*>> clusters;
  for (const auto& worker : batch) {
    clusters[worker->getWorkerId()].push_back(worker.get());
  }
  for (const auto& [cluster_id, workers] : clusters) {
    const FlexDRWorker* best = chosen.at(cluster_id);
    int chosen_index = -1;
    std::string results;
    for (size_t i = 0; i < workers.size(); i++) {
      const FlexDRWorker* worker = workers[i];
      if (worker == best) {
        chosen_index = static_cast<int>(i);
      }
      const odb::Rect& box = worker->getRouteBox();
      const double seconds = worker->getWorkStatistics().total_seconds;
      results += fmt::format("{}{}:{:.3f}:{}:{}:{},{},{},{}:{}",
                             i == 0 ? "" : " ",
                             worker->getBestNumMarkers(),
                             seconds,
                             worker->getWorkerDRCCost(),
                             worker->getWorkerMarkerCost(),
                             box.xMin(),
                             box.yMin(),
                             box.xMax(),
                             box.yMax(),
                             worker->isStopped() ? 'x'
                                                 : (seconds > 0 ? 'r' : 's'));
    }
    logger_->debug(DRT,
                   "bonnroute_cleanup_worker",
                   "cluster iter {} flow {} id {} workers {} chosen {} "
                   "markers {} -> {} results [{}]",
                   iter_,
                   flow_state_machine_->getFlowName(),
                   cluster_id,
                   workers.size(),
                   chosen_index,
                   best->getInitNumMarkers(),
                   best->getBestNumMarkers(),
                   results);
  }
}

odb::Rect FlexDR::getDRVBBox(const odb::Rect& drv_rect) const
{
  odb::Rect route_box = drv_rect;
  auto min_idx = getDesign()->getTopBlock()->getGCellIdx(route_box.ll());
  auto max_idx = getDesign()->getTopBlock()->getGCellIdx(route_box.ur());
  return {min_idx, max_idx};
}

namespace stub_tiles {
/**
 * @brief Clusters DRV boxes.
 *
 * The function uses a greedy algorithm. It starts with a box and tries to merge
 * all other boxes with the following condition:
 * - The resulting merged box is smaller than or equal to 4x4 GCells.
 * It keeps merging untill there is no other possible merges.
 *
 * @param drv_boxes A list of gcell rectangles(with the gcell indices) that hold
 * the violations.
 * @returns a list of merged boxes.
 */
std::vector<odb::Rect> mergeBoxes(std::vector<odb::Rect>& drv_boxes)
{
  std::vector<odb::Rect> merge_boxes{drv_boxes};
  bool merged;
  do {
    merged = false;
    for (auto it1 = merge_boxes.begin(); it1 != merge_boxes.end(); ++it1) {
      for (auto it2 = it1 + 1; it2 != merge_boxes.end(); ++it2) {
        odb::Rect merge_box = (*it1);
        merge_box.merge((*it2));
        if (merge_box.dx() > 4 || merge_box.dy() > 4) {
          continue;
        }
        (*it1).merge((*it2));
        merge_boxes.erase(it2);
        merged = true;
        break;
      }
      if (merged) {
        break;
      }
    }
  } while (merged);
  return merge_boxes;
}
/**
 * @brief Checks if the passed rectangle intersects with any other rectangle in
 * the grid.
 *
 * @param grid A 2-d grid representing the current grid. (-1 value means
 * unoccupied).
 * @param rect_id the current rectangle id to disregard in the check.
 * @return True if the rectangle intersects with an occupied cell by another
 * rectangle.
 */
bool hasOtherRect(const std::vector<std::vector<int>>& grid,
                  const odb::Rect& rect,
                  const int rect_id)
{
  return std::any_of(grid.begin() + rect.xMin(),
                     grid.begin() + rect.xMax() + 1,
                     [&](const std::vector<int>& row) {
                       return std::any_of(row.begin() + rect.yMin(),
                                          row.begin() + rect.yMax() + 1,
                                          [rect_id](int id) {
                                            return id != -1 && id != rect_id;
                                          });
                     });
}
/**
 * @brief Expands the passed rectangle in the required direction.
 *
 * The function expands the passed rectangle in the required direction a gcell
 * at a time with the following conditions:
 * - The expansion is less than or equal to the max_expansion
 * - The expansion does not result a rectangle that intersects with another
 *   rectangle in the grid.
 * @param box The rectangle to be expanded.
 * @param id The rectangle's id.
 * @param grid A 2-d grid representing the current grid. (-1 value means
 * unoccupied).
 * @param dir The direction of the expansion (East, West, North, South).
 * @return The number of expanded gcells (less than or equal to max_expansion)
 */
int expandBox(odb::Rect& box,
              const int id,
              const std::vector<std::vector<int>>& grid,
              const frDirEnum dir,
              const int max_expansion)
{
  auto min_idx = box.ll();
  auto max_idx = box.ur();
  for (int i = 1; i <= max_expansion; i++) {
    switch (dir) {
      case frDirEnum::E:
        max_idx.addX(1);
        break;
      case frDirEnum::W:
        if (min_idx.x() == 0) {
          return i - 1;
        }
        min_idx.addX(-1);
        break;
      case frDirEnum::N:
        max_idx.addY(1);
        break;
      case frDirEnum::S:
        if (min_idx.y() == 0) {
          return i - 1;
        }
        min_idx.addY(-1);
        break;
      default:
        break;
    }
    odb::Rect expanded_box = {min_idx, max_idx};
    if (hasOtherRect(grid, expanded_box, id)) {
      return i - 1;
    }
    box = expanded_box;
  }
  return max_expansion;
}

void populateGrid(std::vector<std::vector<int>>& grid,
                  const odb::Rect& box,
                  const int id)
{
  for (auto x = box.xMin(); x <= box.xMax(); x++) {
    std::fill(
        grid[x].begin() + box.yMin(), grid[x].begin() + box.yMax() + 1, id);
  }
}

struct Wavefront
{
  int id;
  int expansions_done;
  int expansion_east;
  int expansion_west;
  int expansion_north;
  int expansion_south;
  bool operator<(const Wavefront& rhs) const
  {
    // Total order so that priority_queue pop order does not depend on the
    // standard library's heap implementation.
    if (expansions_done != rhs.expansions_done) {
      return expansions_done > rhs.expansions_done;
    }
    if (id != rhs.id) {
      return id > rhs.id;
    }
    return std::tie(
               expansion_east, expansion_west, expansion_north, expansion_south)
           > std::tie(rhs.expansion_east,
                      rhs.expansion_west,
                      rhs.expansion_north,
                      rhs.expansion_south);
  }
};

/**
 * @brief Expand the boxes to reach a max size of 7x7 gcells.
 *
 * The function tries to expand each of the passed merged_boxes into a 7x7 gcell
 * grids. There are 5 variantions of expansions that we consider in this
 * function:
 * - DRVs in Center-Center (Highest priority)
 * - DRVs in Center-East
 * - DRVs in Center-West
 * - DRVs in Center-North
 * - DRVs in Center-South
 * We first expand the first type and then work on the other 4 types of
 * expansions. We prioritize boxes with lower number of expansions done. That
 * way we try to balance the resulting boxes so that each DRV box optimally has
 * the same number of expanded boxes.
 */
std::vector<std::set<odb::Rect>> expandBoxes(
    std::vector<odb::Rect>& merged_boxes)
{
  frUInt4 min_x_idx = std::numeric_limits<frUInt4>::max();
  frUInt4 min_y_idx = std::numeric_limits<frUInt4>::max();
  frUInt4 max_x_idx = std::numeric_limits<frUInt4>::min();
  frUInt4 max_y_idx = std::numeric_limits<frUInt4>::min();
  for (const auto& box : merged_boxes) {
    min_x_idx = std::min(min_x_idx, (frUInt4) box.xMin());
    min_y_idx = std::min(min_y_idx, (frUInt4) box.yMin());
    max_x_idx = std::max(max_x_idx, (frUInt4) box.xMax());
    max_y_idx = std::max(max_y_idx, (frUInt4) box.yMax());
  }
  max_x_idx += 7;
  max_y_idx += 7;
  std::vector<std::vector<int>> grid(max_x_idx);
  std::fill(grid.begin(),
            grid.end(),
            std::vector<int>(max_y_idx, -1));  // -1 means unoccupied
  int id = 0;
  for (const auto& box : merged_boxes) {
    populateGrid(grid, box, id++);
  }

  std::vector<std::set<odb::Rect>> expanded_boxes;
  std::priority_queue<Wavefront> expansions;
  // first we create route boxes centering the violations
  id = 0;
  for (auto box : merged_boxes) {
    const int horizontal_expand = 6 - box.dx();
    const int half_horizontal_expand = horizontal_expand / 2;
    const int vertical_expand = 6 - box.dy();
    const int half_vertical_expand = vertical_expand / 2;
    const int expansion_east
        = expandBox(box, id, grid, frDirEnum::E, half_horizontal_expand);
    const int expansion_west = expandBox(
        box, id, grid, frDirEnum::W, horizontal_expand - expansion_east);
    const int expansion_north
        = expandBox(box, id, grid, frDirEnum::N, half_vertical_expand);
    const int expansion_south = expandBox(
        box, id, grid, frDirEnum::S, vertical_expand - expansion_north);

    expansions.push(
        {id, 1, 0, horizontal_expand, expansion_north, expansion_south});
    expansions.push(
        {id, 1, horizontal_expand, 0, expansion_north, expansion_south});
    expansions.push(
        {id, 1, expansion_east, expansion_west, 0, vertical_expand});
    expansions.push(
        {id, 1, expansion_east, expansion_west, vertical_expand, 0});
    expanded_boxes.push_back({box});
    populateGrid(grid, box, id++);
  }
  while (!expansions.empty()) {
    auto wavefront = expansions.top();
    expansions.pop();
    if (expanded_boxes[wavefront.id].size() != wavefront.expansions_done) {
      wavefront.expansions_done = expanded_boxes[wavefront.id].size();
      expansions.push(wavefront);
      continue;
    }
    auto box = merged_boxes[wavefront.id];
    expandBox(box, wavefront.id, grid, frDirEnum::E, wavefront.expansion_east);
    expandBox(box, wavefront.id, grid, frDirEnum::W, wavefront.expansion_west);
    expandBox(box, wavefront.id, grid, frDirEnum::N, wavefront.expansion_north);
    expandBox(box, wavefront.id, grid, frDirEnum::S, wavefront.expansion_south);
    expanded_boxes[wavefront.id].insert(box);
    populateGrid(grid, box, wavefront.id);
  }
  return expanded_boxes;
}

/**
 * @brief Distributes the workers into batches of non-intersecting borders.
 *
 * The function first calculates the external box to each of the passed route
 * boxes by bloating them by MTSAFEDIST. Then, it calculates the max box of
 * each worker by merging all possible variations of expanded boxes. Finally, it
 * distributes the workers into batches where each batch holds a set of workers
 * that do not intersect.
 */
std::vector<std::vector<int>> getWorkerBatchesBoxes(
    frDesign* design,
    std::vector<std::set<odb::Rect>>& expanded_boxes,
    const frCoord bloating_dist)
{
  if (expanded_boxes.empty()) {
    return {};
  }
  std::vector<odb::Rect> boxes_max;
  boxes_max.resize(expanded_boxes.size());
  int id = 0;
  for (auto& boxes : expanded_boxes) {
    boxes_max[id].mergeInit();
    for (auto rect : boxes) {
      rect.bloat(bloating_dist, rect);
      boxes_max[id].merge(rect);
    }
    ++id;
  }
  std::vector<std::vector<int>> batches;
  batches.push_back({0});
  for (int i = 1; i < boxes_max.size(); i++) {
    auto curr_box = boxes_max[i];
    bool added = false;
    for (auto& batch : batches) {
      bool intersects = false;
      for (const auto& id : batch) {
        if (curr_box.intersects(boxes_max[id])) {
          intersects = true;
          break;
        }
      }
      if (!intersects) {
        batch.emplace_back(i);
        added = true;
        break;
      }
    }
    if (!added) {
      batches.push_back({i});
    }
  }
  return batches;
}
}  // namespace stub_tiles

void FlexDR::stubbornTilesFlow(const SearchRepairArgs& args,
                               IterationProgress& iter_prog)
{
  if (graphics_) {
    graphics_->startIter(iter_, router_cfg_);
  }
  std::vector<odb::Rect> drv_boxes;
  for (const auto& marker : getDesign()->getTopBlock()->getMarkers()) {
    auto box = marker->getBBox();
    drv_boxes.push_back(getDRVBBox(box));
  }
  auto merged_boxes = stub_tiles::mergeBoxes(drv_boxes);
  auto expanded_boxes = stub_tiles::expandBoxes(merged_boxes);

  // Convert gcell indices to actual coordinates
  std::vector<std::set<odb::Rect>> expanded_boxes_coords;
  expanded_boxes_coords.reserve(expanded_boxes.size());
  std::ranges::transform(
      expanded_boxes,
      std::back_inserter(expanded_boxes_coords),
      [this](const auto& box_set) {
        std::set<odb::Rect> coord_set;
        std::ranges::transform(
            box_set,
            std::inserter(coord_set, coord_set.end()),
            [this](const auto& gcell_box) {
              return odb::Rect(
                  getDesign()->getTopBlock()->getGCellBox(gcell_box.ll()).ll(),
                  getDesign()->getTopBlock()->getGCellBox(gcell_box.ur()).ur());
            });
        return coord_set;
      });

  auto route_boxes_batches = stub_tiles::getWorkerBatchesBoxes(
      getDesign(), expanded_boxes_coords, router_cfg_->MTSAFEDIST);
  std::vector<frUInt4> drc_costs
      = {args.workerDRCCost, args.workerDRCCost / 2, args.workerDRCCost * 2};
  std::vector<frUInt4> marker_costs = {args.workerMarkerCost,
                                       args.workerMarkerCost / 2,
                                       args.workerMarkerCost * 2};
  std::vector<std::vector<std::unique_ptr<FlexDRWorker>>> workers_batches(
      route_boxes_batches.size());
  iter_prog.total_num_workers = 0;
  for (int batch_id = 0; batch_id < route_boxes_batches.size(); batch_id++) {
    auto& batch = route_boxes_batches[batch_id];
    for (const auto worker_id : batch) {
      for (const auto& route_box : expanded_boxes_coords[worker_id]) {
        for (auto drc_cost : drc_costs) {
          for (auto marker_cost : marker_costs) {
            auto worker_args = args;
            worker_args.workerDRCCost = drc_cost;
            worker_args.workerMarkerCost = marker_cost;
            workers_batches[batch_id].emplace_back(
                createWorker(0, 0, worker_args, route_box));
            workers_batches[batch_id].back()->setWorkerId(worker_id);
            iter_prog.total_num_workers++;
          }
        }
      }
    }
  }
  bool changed = false;
  omp_set_num_threads(router_cfg_->MAX_THREADS);
  for (auto& batch : workers_batches) {
    // The selection below keeps, per cluster (worker id), the first worker in
    // creation order with the fewest markers.  A worker that ends with zero
    // markers is therefore the selected one whatever the later workers of its
    // cluster find.  Workers only read the design until end(), so the later
    // workers of a cluster are skipped or stopped once an earlier one has
    // ended with zero markers; the selection runs over the workers that ran to
    // the end, a set that holds every worker up to the first zero-marker one,
    // so the selected worker is the same as when every worker runs, at any
    // thread count and timing.
    const std::vector<int> ran = runStubbornWorkersUntilClean(batch, iter_prog);
    std::map<int, FlexDRWorker*>
        worker_best_result;  // holds the best worker in results
    for (const int index : ran) {
      auto& worker = batch[index];
      const int worker_id = worker->getWorkerId();
      if (worker_best_result.find(worker_id) == worker_best_result.end()
          || worker->getBestNumMarkers()
                 < worker_best_result[worker_id]->getBestNumMarkers()) {
        worker_best_result[worker_id] = worker.get();
      }
    }
    if (logger_->debugCheck(DRT, "bonnroute_cleanup_worker", 1)) {
      reportClusterResults(batch, worker_best_result);
    }
    for (auto [_, worker] : worker_best_result) {
      changed
          |= (worker->end(getDesign())
              && worker->getBestNumMarkers() != worker->getInitNumMarkers());
    }
    batch.clear();
  }
  flow_state_machine_->setLastIterationEffective(changed);
}

namespace {

// Snap a segment to the gcells it crosses and grow it by one gcell at each end
// along its own routing axis (vertical -> up/down, horizontal -> left/right),
// giving the worker room to relocate the vias at the segment ends.
odb::Rect segmentBox(frBlock* block, frPathSeg* seg)
{
  const odb::Rect bbox = seg->getBBox();
  odb::Point ll = block->getGCellIdx(bbox.ll());
  odb::Point ur = block->getGCellIdx(bbox.ur());
  if (seg->isVertical()) {
    ll = {ll.x(), std::max(0, ll.y() - 1)};
    ur = {ur.x(), ur.y() + 1};
  } else {
    ll = {std::max(0, ll.x() - 1), ll.y()};
    ur = {ur.x() + 1, ur.y()};
  }
  return {block->getGCellBox(ll).ll(), block->getGCellBox(ur).ur()};
}

// Collect the net's routed objects (segments, patches, vias) overlapping a
// rect on the given layer. A via carries metal on more than one layer, so it
// is matched only when the enclosure (or cut) that actually sits on lnum
// intersects the rect.
void netObjsOverlapping(frNet* net,
                        frLayerNum lnum,
                        const odb::Rect& rect,
                        std::vector<frBlockObject*>& out)
{
  for (const auto& via : net->getVias()) {
    const frViaDef* via_def = via->getViaDef();
    const bool hit = (via_def->getLayer1Num() == lnum
                      && via->getLayer1BBox().intersects(rect))
                     || (via_def->getLayer2Num() == lnum
                         && via->getLayer2BBox().intersects(rect))
                     || (via_def->getCutLayerNum() == lnum
                         && via->getCutBBox().intersects(rect));
    if (hit) {
      out.push_back(via.get());
    }
  }
  for (const auto& shape : net->getShapes()) {
    if (shape->typeId() == frcPathSeg && shape->getLayerNum() == lnum
        && shape->getBBox().intersects(rect)) {
      out.push_back(shape.get());
    }
  }
  for (const auto& pwire : net->getPatchWires()) {
    if (pwire->getLayerNum() == lnum && pwire->getBBox().intersects(rect)) {
      out.push_back(pwire.get());
    }
  }
}

// Net path segments whose endpoint coincides with the via origin (on either
// the via's bottom or top routing layer).
void connectedSegments(frNet* net,
                       const odb::Point& origin,
                       std::vector<frPathSeg*>& out)
{
  for (const auto& shape : net->getShapes()) {
    if (shape->typeId() != frcPathSeg) {
      continue;
    }
    auto* seg = static_cast<frPathSeg*>(shape.get());
    if (seg->getBeginPoint() == origin || seg->getEndPoint() == origin) {
      out.push_back(seg);
    }
  }
}

// Reduce a src object to segment-anchored worker boxes:
//  - segment: its expanded gcell box
//  - via:     the expanded gcell box of every segment connected to it
//  - patch:   resolve the via/segment it patches and recurse once
void addObjBoxes(frBlock* block,
                 frBlockObject* obj,
                 std::vector<odb::Rect>& boxes,
                 int depth)
{
  switch (obj->typeId()) {
    case frcPathSeg:
      boxes.push_back(segmentBox(block, static_cast<frPathSeg*>(obj)));
      break;
    case frcVia: {
      auto* via = static_cast<frVia*>(obj);
      std::vector<frPathSeg*> segs;
      connectedSegments(via->getNet(), via->getOrigin(), segs);
      for (auto* seg : segs) {
        boxes.push_back(segmentBox(block, seg));
      }
      break;
    }
    case frcPatchWire: {
      if (depth > 0) {
        break;  // guard against patch-on-patch recursion
      }
      auto* pw = static_cast<frPatchWire*>(obj);
      std::vector<frBlockObject*> patched;
      netObjsOverlapping(
          pw->getNet(), pw->getLayerNum(), pw->getBBox(), patched);
      for (auto* p : patched) {
        if (p != obj) {
          addObjBoxes(block, p, boxes, depth + 1);
        }
      }
      break;
    }
    default:
      break;
  }
}

}  // namespace

// Builds worker boxes for a marker whose routing diverged off-guide (no orig
// guide covers it). The marker only records net owners and shape rects, so the
// actual src objects (segment/via/patch) are resolved by walking the shapes of
// each guided src net near the marker, then each is reduced to a segment gcell
// box that gives the worker room to reroute.
std::vector<odb::Rect> FlexDR::getOffGuideWorkerBoxes(frMarker* marker) const
{
  auto* block = getDesign()->getTopBlock();

  std::vector<frBlockObject*> src_objs;
  for (auto& obj : marker->getSrcs()) {
    if (obj->typeId() != frcNet) {
      continue;
    }
    auto* net = static_cast<frNet*>(obj);
    if (net->getOrigGuides().empty()) {
      continue;  // e.g. a power net whose rail can't be moved
    }
    netObjsOverlapping(net, marker->getLayerNum(), marker->getBBox(), src_objs);
  }

  std::vector<odb::Rect> boxes;
  for (auto* obj : src_objs) {
    addObjBoxes(block, obj, boxes, 0);
  }
  return boxes;
}

void FlexDR::guideTilesFlow(const SearchRepairArgs& args,
                            IterationProgress& iter_prog)
{
  if (graphics_) {
    graphics_->startIter(iter_, router_cfg_);
  }
  std::vector<odb::Rect> workers;
  for (const auto& marker : getDesign()->getTopBlock()->getMarkers()) {
    bool covered_by_guide = false;
    bool has_guided_net = false;
    for (auto src : marker->getSrcs()) {
      if (src->typeId() == frcNet) {
        auto net = static_cast<frNet*>(src);
        if (net->getOrigGuides().empty()) {
          continue;
        }
        has_guided_net = true;
        for (const auto& guide : net->getOrigGuides()) {
          if (!guide.getBBox().intersects(marker->getBBox())) {
            continue;
          }
          workers.push_back(guide.getBBox());
          covered_by_guide = true;
        }
      }
    }
    // The router can diverge from its guides, leaving a violation that no
    // guide covers. Guide-anchored workers miss it entirely, so build worker
    // boxes from the src objects behind the violation instead.
    if (has_guided_net && !covered_by_guide) {
      for (const auto& box : getOffGuideWorkerBoxes(marker.get())) {
        workers.push_back(box);
      }
    }
  }
  if (workers.empty()) {
    iter_prog.total_num_workers = 1;
    iter_prog.cnt_done_workers = 1;
    flow_state_machine_->setLastIterationEffective(false);
    return;
  }

  for (auto itr1 = workers.begin(); itr1 != workers.end(); ++itr1) {
    for (auto itr2 = itr1 + 1; itr2 != workers.end(); ++itr2) {
      if (itr1->getDir() == itr2->getDir() && itr1->intersects(*itr2)) {
        itr1->merge(*itr2);
        workers.erase(itr2);
        itr2 = itr1;
      }
    }
  }
  std::vector<std::set<odb::Rect>> boxes_set;
  boxes_set.reserve(workers.size());
  for (auto& worker : workers) {
    boxes_set.push_back({worker});
  }
  auto route_boxes_batches = stub_tiles::getWorkerBatchesBoxes(
      getDesign(), boxes_set, router_cfg_->MTSAFEDIST);
  std::vector<std::vector<std::unique_ptr<FlexDRWorker>>> workers_batches(
      route_boxes_batches.size());
  iter_prog.total_num_workers = 0;
  for (int batch_id = 0; batch_id < route_boxes_batches.size(); batch_id++) {
    auto& batch = route_boxes_batches[batch_id];
    for (const auto worker_id : batch) {
      for (auto route_box : boxes_set[worker_id]) {
        workers_batches[batch_id].emplace_back(
            createWorker(0, 0, args, route_box));
        workers_batches[batch_id].back()->setWorkerId(worker_id);
        iter_prog.total_num_workers++;
      }
    }
  }
  bool changed = false;
  omp_set_num_threads(router_cfg_->MAX_THREADS);
  for (auto& batch : workers_batches) {
    processWorkersBatch(batch, iter_prog);
    std::map<int, FlexDRWorker*>
        worker_best_result;  // holds the best worker in results
    for (auto& worker : batch) {
      const int worker_id = worker->getWorkerId();
      if (worker_best_result.find(worker_id) == worker_best_result.end()
          || worker->getBestNumMarkers()
                 < worker_best_result[worker_id]->getBestNumMarkers()) {
        worker_best_result[worker_id] = worker.get();
      }
    }
    if (logger_->debugCheck(DRT, "bonnroute_cleanup_worker", 1)) {
      reportClusterResults(batch, worker_best_result);
    }
    for (auto [_, worker] : worker_best_result) {
      changed
          |= (worker->end(getDesign())
              && worker->getBestNumMarkers() != worker->getInitNumMarkers());
    }
    batch.clear();
  }
  flow_state_machine_->setLastIterationEffective(changed);
}
void FlexDR::optimizationFlow(const SearchRepairArgs& args,
                              IterationProgress& iter_prog)
{
  if (graphics_) {
    graphics_->startIter(iter_, router_cfg_);
  }
  auto gCellPatterns = getDesign()->getTopBlock()->getGCellPatterns();
  auto& xgp = gCellPatterns.at(0);
  auto& ygp = gCellPatterns.at(1);
  const int size = args.size;
  const int offset = args.offset;
  iter_prog.total_num_workers
      = (((int) xgp.getCount() - 1 - offset) / size + 1)
        * (((int) ygp.getCount() - 1 - offset) / size + 1);

  std::vector<std::unique_ptr<FlexDRWorker>> uworkers;
  int batchStepX, batchStepY;

  getBatchInfo(batchStepX, batchStepY);

  std::vector<std::vector<std::vector<std::unique_ptr<FlexDRWorker>>>> workers(
      batchStepX * batchStepY);
  // The gcell offsets of the workers, in the shape of workers.
  std::vector<std::vector<std::vector<std::pair<int, int>>>> worker_offsets(
      batchStepX * batchStepY);
  const double create_start = omp_get_wtime();

  // In the incremental iterations a worker whose extended box sees no change
  // and no violation has no net to route and would re-check unchanged,
  // violation-free geometry, so it is skipped.  Workers are checked just
  // before their batch runs, after the earlier batches' changes.
  // After BonnRoute's incremental routing (arXiv:2111.06169 Sec. 5.2) the
  // repair iterations run only the workers with markers as well.
  const bool marker_cleanup
      = router_ != nullptr && router_->isBonnRouteCleanup();
  const bool skip_unaffected
      = (args.ripupMode == RipUpMode::INCR || marker_cleanup) && !dist_on_
        && incremental_changed_regions_.has_value();
  const bool markers_only = args.ripupMode != RipUpMode::INCR;
  std::vector<odb::Rect>* changed_regions
      = skip_unaffected ? &incremental_changed_regions_.value() : nullptr;
  const int num_workers = iter_prog.total_num_workers;
  int num_unaffected_workers = 0;

  int xIdx = 0, yIdx = 0;
  for (int i = offset; i < (int) xgp.getCount(); i += size) {
    for (int j = offset; j < (int) ygp.getCount(); j += size) {
      auto worker = createWorker(i, j, args);
      int batch_idx = (xIdx % batchStepX) * batchStepY + yIdx % batchStepY;
      const bool create_new_batch
          = workers[batch_idx].empty()
            || (!dist_on_
                && workers[batch_idx].back().size() >= router_cfg_->BATCHSIZE);
      if (create_new_batch) {
        workers[batch_idx].push_back(
            std::vector<std::unique_ptr<FlexDRWorker>>());
        worker_offsets[batch_idx].emplace_back();
      }
      workers[batch_idx].back().push_back(std::move(worker));
      worker_offsets[batch_idx].back().emplace_back(i, j);

      yIdx++;
    }
    yIdx = 0;
    xIdx++;
  }

  cleanup_statistics_.create_seconds += omp_get_wtime() - create_start;
  omp_set_num_threads(router_cfg_->MAX_THREADS);
  int version = 0;
  increaseClipsize_ = false;
  numWorkUnits_ = 0;
  num_check_only_workers_ = 0;
  // The marker-driven iterations of the clean-up after BonnRoute run few
  // workers per batch, so workers of later batches start ahead on idle
  // threads (runGridWorkersAhead; same results as the batch loop below).
  const bool start_ahead
      = initial_routing_cleanup_.start_grid_workers_ahead && markers_only
        && iter_ >= 2 && !dist_on_ && graphics_ == nullptr
        && (!skip_unaffected || marker_cleanup);
  if (start_ahead) {
    std::vector<GridWorkerSlot> slots;
    int num_batches = 0;
    for (int group = 0; group < (int) workers.size(); group++) {
      for (int batch = 0; batch < (int) workers[group].size(); batch++) {
        for (int k = 0; k < (int) workers[group][batch].size(); k++) {
          slots.push_back({.worker = std::move(workers[group][batch][k]),
                           .x_offset = worker_offsets[group][batch][k].first,
                           .y_offset = worker_offsets[group][batch][k].second,
                           .batch = num_batches});
        }
        num_batches++;
      }
    }
    workers.clear();
    const double batch_start = omp_get_wtime();
    num_unaffected_workers = runGridWorkersAhead(
        slots, num_batches, args, iter_prog, skip_unaffected, changed_regions);
    cleanup_statistics_.batch_seconds += omp_get_wtime() - batch_start;
  }
  // parallel execution
  for (auto& workerBatch : workers) {
    ProfileTask profile("DR:checkerboard");
    for (auto& workersInBatch : workerBatch) {
      if (skip_unaffected) {
        const int num_unaffected
            = std::erase_if(workersInBatch,
                            [this, markers_only, marker_cleanup](
                                const std::unique_ptr<FlexDRWorker>& worker) {
                              return !isAffectedByChanges(
                                  *worker, markers_only, marker_cleanup);
                            });
        num_unaffected_workers += num_unaffected;
        iter_prog.total_num_workers -= num_unaffected;
      }
      {
        const std::string batch_name = std::string("DR:batch<")
                                       + std::to_string(workersInBatch.size())
                                       + ">";
        ProfileTask profile(batch_name.c_str());
        const double batch_start = omp_get_wtime();
        if (dist_on_) {
          processWorkersBatchDistributed(workersInBatch, version, iter_prog);
        } else {
          processWorkersBatch(workersInBatch, iter_prog);
        }
        cleanup_statistics_.batch_seconds += omp_get_wtime() - batch_start;
      }
      const double end_batch_start = omp_get_wtime();
      endWorkersBatch(workersInBatch, changed_regions);
      cleanup_statistics_.end_batch_seconds
          += omp_get_wtime() - end_batch_start;
    }
  }

  if (skip_unaffected) {
    debugPrint(logger_,
               DRT,
               "incremental",
               1,
               "Skipped {} of {} workers.",
               num_unaffected_workers,
               num_workers);
    if (router_cfg_->VERBOSE > 0) {
      logger_->info(DRT,
                    629,
                    "Skipped {} of {} workers unaffected by changes since the "
                    "previous detailed routing.",
                    num_unaffected_workers,
                    num_workers);
    }
  }

  if (num_check_only_workers_ > 0 && router_cfg_->VERBOSE > 0) {
    logger_->info(DRT,
                  9121,
                  "{} of {} workers had only nets with initial routing and "
                  "ran the design-rule check only.",
                  num_check_only_workers_,
                  num_workers);
  }

  if (!iter_) {
    removeGCell2BoundaryPin();
  }
}

bool FlexDR::isAffectedByChanges(const FlexDRWorker& worker,
                                 const bool markers_only,
                                 const bool shape_regions) const
{
  odb::Rect affected_box;
  if (shape_regions) {
    // A changed shape can only cause markers within MTSAFEDIST of it, and
    // the worker reports the markers in its DRC box. Unrouted nets' regions
    // are whole gcell boxes, so the workers routing them intersect them.
    worker.getDrcBox().bloat(router_cfg_->MTSAFEDIST, affected_box);
  } else {
    // Guides are stored between gcell centers; keep one gcell of margin on
    // top of the region the worker reads.
    const auto& gcell_patterns = getDesign()->getTopBlock()->getGCellPatterns();
    const int margin
        = router_cfg_->MTSAFEDIST
          + static_cast<int>(std::max(gcell_patterns[0].getSpacing(),
                                      gcell_patterns[1].getSpacing()));
    worker.getExtBox().bloat(margin, affected_box);
  }
  if (!markers_only
      && std::ranges::any_of(*incremental_changed_regions_,
                             [&affected_box](const odb::Rect& region) {
                               return affected_box.intersects(region);
                             })) {
    return true;
  }
  std::vector<frMarker*> markers;
  getDesign()->getRegionQuery()->queryMarker(worker.getDrcBox(), markers);
  return !markers.empty();
}

void FlexDR::searchRepair(const SearchRepairArgs& args)
{
  // Calculate flow state
  const auto flow_state = flow_state_machine_->determineNextFlow(
      {.num_violations = getDesign()->getTopBlock()->getNumMarkers(),
       .args = args});

  const RipUpMode ripupMode = args.ripupMode;
  if ((ripupMode == RipUpMode::DRC || ripupMode == RipUpMode::NEARDRC)
      && getDesign()->getTopBlock()->getMarkers().empty()) {
    return;
  }
  ProfileTask profile(fmt::format("DR:searchRepair{}", iter_).c_str());
  const double iteration_start = omp_get_wtime();
  cleanup_statistics_ = {};

  if (dist_on_) {
    if ((iter_ % 10 == 0 && iter_ != 60) || iter_ == 3 || iter_ == 15) {
      router_cfg_path_ = fmt::format("{}globals.{}.ar", dist_dir_, iter_);
      router_->writeGlobals(router_cfg_path_);
    }
  }

  if (router_cfg_->VERBOSE > 0) {
    if (flow_state != FlexDRFlow::State::SKIP) {
      printIteration(logger_, iter_, flow_state_machine_->getFlowName());
    } else {
      logger_->info(DRT, 200, "Skipping iteration {}", iter_);
    }
  }
  // start timer for the current iteration
  IterationProgress iter_prog;

  switch (flow_state) {
    case FlexDRFlow::State::GUIDES:
      guideTilesFlow(args, iter_prog);
      break;

    case FlexDRFlow::State::STUBBORN:
      stubbornTilesFlow(args, iter_prog);
      break;

    case FlexDRFlow::State::OPTIMIZATION:
      optimizationFlow(args, iter_prog);
      break;

    case FlexDRFlow::State::SKIP:
      return;
  }

  if (router_cfg_->VERBOSE > 0) {
    iter_prog.cnt_done_workers--;  // decrement 1 and increment again in
                                   // printIterationProgress
    printIterationProgress(
        logger_, iter_prog, getDesign()->getTopBlock()->getNumMarkers(), 100);
  }
  FlexDRConnectivityChecker checker(
      router_, logger_, router_cfg_, graphics_.get(), dist_on_);
  const double check_start = omp_get_wtime();
  checker.check(iter_);
  const double check_seconds = omp_get_wtime() - check_start;
  flow_state_machine_->setFixingMaxSpacing(false);
  if (getDesign()->getTopBlock()->getNumMarkers() == 0
      && getTech()->hasMaxSpacingConstraints()) {
    fixMaxSpacing();
  }
  numViols_.push_back(getDesign()->getTopBlock()->getNumMarkers());
  if (logger_->debugCheck(DRT, "bonnroute_cleanup", 1)) {
    const auto& conn = checker.phaseSeconds();
    logger_->debug(
        DRT,
        "bonnroute_cleanup",
        "cleanup iter {} flow {} ripup_mode {} workers {} without_markers {} "
        "without_nets {} only_initial_routing {} check_only {} markers {} "
        "create_s {:.3f} batches_s {:.3f} end_batch_s {:.3f} check_s {:.3f} "
        "check_init_s {:.3f} check_merge_s {:.3f} check_astar_s {:.3f} "
        "check_finish_s {:.3f} iteration_s {:.3f}",
        iter_,
        flow_state_machine_->getFlowName(),
        static_cast<int>(args.ripupMode),
        cleanup_statistics_.workers,
        cleanup_statistics_.workers_without_markers,
        cleanup_statistics_.workers_without_nets,
        cleanup_statistics_.workers_with_only_initial_routing,
        num_check_only_workers_,
        getDesign()->getTopBlock()->getNumMarkers(),
        cleanup_statistics_.create_seconds,
        cleanup_statistics_.batch_seconds,
        cleanup_statistics_.end_batch_seconds,
        check_seconds,
        conn.init_parallel,
        conn.merge_serial,
        conn.astar_parallel,
        conn.finish_serial,
        omp_get_wtime() - iteration_start);
  }
  debugPrint(logger_,
             utl::DRT,
             "workers",
             1,
             "Number of work units = {}.",
             numWorkUnits_);
  reportIterationViolations();
  if (router_cfg_->VERBOSE > 0) {
    iter_prog.time.print(logger_);
    std::cout << std::flush;
  }
  end();
}

void FlexDR::end(bool done)
{
  if (done) {
    router_->reportDRC(
        router_cfg_->DRC_RPT_FILE, design_->getTopBlock()->getMarkers(), "DRC");
  }
  if (done && router_cfg_->VERBOSE > 0) {
    logger_->info(DRT, 198, "Complete detail routing.");
  }

  using ULL = uint64_t;
  const auto size = getTech()->getLayers().size();
  std::vector<ULL> wlen(size, 0);
  std::vector<ULL> sCut(size, 0);
  std::vector<ULL> mCut(size, 0);
  auto topBlock = getDesign()->getTopBlock();
  for (auto& net : topBlock->getNets()) {
    for (auto& shape : net->getShapes()) {
      if (shape->typeId() == frcPathSeg) {
        auto obj = static_cast<frPathSeg*>(shape.get());
        auto [bp, ep] = obj->getPoints();
        auto lNum = obj->getLayerNum();
        frCoord psLen = odb::Point::manhattanDistance(ep, bp);
        wlen[lNum] += psLen;
      }
    }
    for (auto& via : net->getVias()) {
      auto lNum = via->getViaDef()->getCutLayerNum();
      if (via->getViaDef()->isMultiCut()) {
        ++mCut[lNum];
      } else {
        ++sCut[lNum];
      }
    }
  }

  const ULL totWlen = std::accumulate(wlen.begin(), wlen.end(), ULL(0));
  const ULL totSCut = std::accumulate(sCut.begin(), sCut.end(), ULL(0));
  const ULL totMCut = std::accumulate(mCut.begin(), mCut.end(), ULL(0));

  if (done) {
    logger_->metric("route__drc_errors", topBlock->getNumMarkers());
    logger_->metric("route__wirelength", totWlen / topBlock->getDBUPerUU());
    logger_->metric("route__vias", totSCut + totMCut);
    logger_->metric("route__vias__singlecut", totSCut);
    logger_->metric("route__vias__multicut", totMCut);
  } else {
    logger_->metric(fmt::format("route__drc_errors__iter:{}", iter_),
                    topBlock->getNumMarkers());
    logger_->metric(fmt::format("route__wirelength__iter:{}", iter_),
                    totWlen / topBlock->getDBUPerUU());
  }

  if (router_cfg_->VERBOSE > 0) {
    logger_->report("Total wire length = {} um.",
                    totWlen / topBlock->getDBUPerUU());

    for (int i = getTech()->getBottomLayerNum();
         i <= getTech()->getTopLayerNum();
         i++) {
      if (getTech()->getLayer(i)->getType() == dbTechLayerType::ROUTING) {
        logger_->report("Total wire length on LAYER {} = {} um.",
                        getTech()->getLayer(i)->getName(),
                        wlen[i] / topBlock->getDBUPerUU());
      }
    }
    logger_->report("Total number of vias = {}.", totSCut + totMCut);
    if (totMCut > 0) {
      logger_->report("Total number of multi-cut vias = {} ({:5.1f}%).",
                      totMCut,
                      totMCut * 100.0 / (totSCut + totMCut));
      logger_->report("Total number of single-cut vias = {} ({:5.1f}%).",
                      totSCut,
                      totSCut * 100.0 / (totSCut + totMCut));
    }
    logger_->report("Up-via summary (total {}):", totSCut + totMCut);
    int nameLen = 0;
    for (int i = getTech()->getBottomLayerNum();
         i <= getTech()->getTopLayerNum();
         i++) {
      if (getTech()->getLayer(i)->getType() == dbTechLayerType::CUT) {
        nameLen = std::max(nameLen,
                           (int) getTech()->getLayer(i - 1)->getName().size());
      }
    }
    int maxL = 1 + nameLen + 4 + (int) std::to_string(totSCut).length();
    if (totMCut) {
      maxL += 9 + 4 + (int) std::to_string(totMCut).length() + 9 + 4
              + (int) std::to_string(totSCut + totMCut).length();
    }
    std::ostringstream msg;
    if (totMCut) {
      msg << " "
          << std::setw(nameLen + 4 + (int) std::to_string(totSCut).length() + 9)
          << "single-cut";
      msg << std::setw(4 + (int) std::to_string(totMCut).length() + 9)
          << "multi-cut"
          << std::setw(4 + (int) std::to_string(totSCut + totMCut).length())
          << "total";
    }
    msg << std::endl;
    for (int i = 0; i < maxL; i++) {
      msg << "-";
    }
    msg << std::endl;
    for (int i = getTech()->getBottomLayerNum();
         i <= getTech()->getTopLayerNum();
         i++) {
      if (getTech()->getLayer(i)->getType() == dbTechLayerType::CUT) {
        msg << " " << std::setw(nameLen)
            << getTech()->getLayer(i - 1)->getName() << "    "
            << std::setw((int) std::to_string(totSCut).length()) << sCut[i];
        if (totMCut) {
          msg << " (" << std::setw(5)
              << (double) ((sCut[i] + mCut[i])
                               ? sCut[i] * 100.0 / (sCut[i] + mCut[i])
                               : 0.0)
              << "%)";
          msg << "    " << std::setw((int) std::to_string(totMCut).length())
              << mCut[i] << " (" << std::setw(5)
              << (double) ((sCut[i] + mCut[i])
                               ? mCut[i] * 100.0 / (sCut[i] + mCut[i])
                               : 0.0)
              << "%)    "
              << std::setw((int) std::to_string(totSCut + totMCut).length())
              << sCut[i] + mCut[i];
        }
        msg << std::endl;
      }
    }
    for (int i = 0; i < maxL; i++) {
      msg << "-";
    }
    msg << std::endl;
    msg << " " << std::setw(nameLen) << "    "
        << std::setw((int) std::to_string(totSCut).length()) << totSCut;
    if (totMCut) {
      msg << " (" << std::setw(5)
          << ((totSCut + totMCut) ? totSCut * 100.0 / (totSCut + totMCut) : 0.0)
          << "%)";
      msg << "    " << std::setw((int) std::to_string(totMCut).length())
          << totMCut << " (" << std::setw(5)
          << ((totSCut + totMCut) ? totMCut * 100.0 / (totSCut + totMCut) : 0.0)
          << "%)    "
          << std::setw((int) std::to_string(totSCut + totMCut).length())
          << totSCut + totMCut;
    }
    msg << std::endl << std::endl;
    logger_->report("{}", msg.str());
  }
}

std::vector<FlexDR::SearchRepairArgs> strategy(const frUInt4 shapeCost,
                                               const frUInt4 markerCost)
{
  // clang-format off
  return {
    {7,  0,  3,      shapeCost,               0,       shapeCost, 0.950, RipUpMode::ALL    ,  true}, //  0
    {7, -2,  3,      shapeCost,       shapeCost,       shapeCost, 0.950, RipUpMode::ALL    ,  true}, //  1
    {7, -5,  3,      shapeCost,       shapeCost,       shapeCost, 0.950, RipUpMode::ALL    ,  true}, //  2
    {7,  0,  8,      shapeCost,      markerCost,   2 * shapeCost, 0.950, RipUpMode::DRC    , false}, //  3
    {7, -1,  8,      shapeCost,      markerCost,   2 * shapeCost, 0.950, RipUpMode::DRC    , false}, //  4
    {7, -2,  8,      shapeCost,      markerCost,   2 * shapeCost, 0.950, RipUpMode::DRC    , false}, //  5
    {7, -3,  8,      shapeCost,      markerCost,   2 * shapeCost, 0.950, RipUpMode::DRC    , false}, //  6
    {7, -4,  8,      shapeCost,      markerCost,   2 * shapeCost, 0.950, RipUpMode::DRC    , false}, //  7
    {7, -5,  8,      shapeCost,      markerCost,   2 * shapeCost, 0.950, RipUpMode::DRC    , false}, //  8
    {7, -6,  8,      shapeCost,      markerCost,   2 * shapeCost, 0.950, RipUpMode::DRC    , false}, //  9
    {7,  0,  8,  2 * shapeCost,      markerCost,   3 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 10
    {7, -1,  8,  2 * shapeCost,      markerCost,   3 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 11
    {7, -2,  8,  2 * shapeCost,      markerCost,   3 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 12
    {7, -3,  8,  2 * shapeCost,      markerCost,   3 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 13
    {7, -4,  8,  2 * shapeCost,      markerCost,   3 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 14
    {7, -5,  8,  2 * shapeCost,      markerCost,   4 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 15
    {7, -6,  8,  2 * shapeCost,      markerCost,   4 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 16
    {7, -3,  8,      shapeCost,      markerCost,   4 * shapeCost, 0.950, RipUpMode::ALL    , false}, // 17
    {7,  0,  8,  4 * shapeCost,      markerCost,   4 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 18
    {7, -1,  8,  4 * shapeCost,      markerCost,   4 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 19
    {7, -2,  8,  4 * shapeCost,      markerCost,  10 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 20
    {7, -3,  8,  4 * shapeCost,      markerCost,  10 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 21
    {7, -4,  8,  4 * shapeCost,      markerCost,  10 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 22
    {7, -5,  8,      shapeCost,      markerCost,  10 * shapeCost, 0.950, RipUpMode::NEARDRC, false}, // 23
    {7, -6,  8,  4 * shapeCost,      markerCost,  10 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 24
    {5, -2,  8,      shapeCost,      markerCost,  10 * shapeCost, 0.950, RipUpMode::ALL    , false}, // 25
    {7,  0,  8,  8 * shapeCost,  2 * markerCost,  10 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 26
    {7, -1,  8,  8 * shapeCost,  2 * markerCost,  10 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 27
    {7, -2,  8,  8 * shapeCost,  2 * markerCost,  10 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 28
    {7, -3,  8,  8 * shapeCost,  2 * markerCost,  10 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 29
    {7, -4,  8,      shapeCost,      markerCost,  50 * shapeCost, 0.950, RipUpMode::NEARDRC, false}, // 30
    {7, -5,  8,  8 * shapeCost,  2 * markerCost,  50 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 31
    {7, -6,  8,  8 * shapeCost,  2 * markerCost,  50 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 32
    {3, -1,  8,      shapeCost,      markerCost,  50 * shapeCost, 0.950, RipUpMode::ALL    , false}, // 33
    {7,  0,  8, 16 * shapeCost,  4 * markerCost,  50 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 34
    {7, -1,  8, 16 * shapeCost,  4 * markerCost,  50 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 35
    {7, -2,  8, 16 * shapeCost,  4 * markerCost,  50 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 36
    {7, -3,  8,      shapeCost,      markerCost,  50 * shapeCost, 0.950, RipUpMode::NEARDRC, false}, // 37
    {7, -4,  8, 16 * shapeCost,  4 * markerCost,  50 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 38
    {7, -5,  8, 16 * shapeCost,  4 * markerCost,  50 * shapeCost, 0.950, RipUpMode::DRC    , false}, // 39
    {7, -6,  8, 16 * shapeCost,  4 * markerCost, 100 * shapeCost, 0.990, RipUpMode::DRC    , false}, // 40
    {3, -2,  8,      shapeCost,      markerCost, 100 * shapeCost, 0.990, RipUpMode::ALL    , false}, // 41
    {7,  0, 16, 16 * shapeCost,  4 * markerCost, 100 * shapeCost, 0.990, RipUpMode::DRC    , false}, // 42
    {7, -1, 16, 16 * shapeCost,  4 * markerCost, 100 * shapeCost, 0.990, RipUpMode::DRC    , false}, // 43
    {7, -2, 16,      shapeCost,      markerCost, 100 * shapeCost, 0.990, RipUpMode::NEARDRC, false}, // 44
    {7, -3, 16, 16 * shapeCost,  4 * markerCost, 100 * shapeCost, 0.990, RipUpMode::DRC    , false}, // 45
    {7, -4, 16, 16 * shapeCost,  4 * markerCost, 100 * shapeCost, 0.990, RipUpMode::DRC    , false}, // 46
    {7, -5, 16, 16 * shapeCost,  4 * markerCost, 100 * shapeCost, 0.990, RipUpMode::DRC    , false}, // 47
    {7, -6, 16, 16 * shapeCost,  4 * markerCost, 100 * shapeCost, 0.990, RipUpMode::DRC    , false}, // 48
    {3, -0,  8,      shapeCost,      markerCost, 100 * shapeCost, 0.990, RipUpMode::ALL    , false}, // 49
    {7,  0, 32, 32 * shapeCost,  8 * markerCost, 100 * shapeCost, 0.999, RipUpMode::DRC    , false}, // 50
    {7, -1, 32,      shapeCost,      markerCost, 100 * shapeCost, 0.999, RipUpMode::NEARDRC, false}, // 51
    {7, -2, 32, 32 * shapeCost,  8 * markerCost, 100 * shapeCost, 0.999, RipUpMode::DRC    , false}, // 52
    {7, -3, 32, 32 * shapeCost,  8 * markerCost, 100 * shapeCost, 0.999, RipUpMode::DRC    , false}, // 53
    {7, -4, 32, 32 * shapeCost,  8 * markerCost, 100 * shapeCost, 0.999, RipUpMode::DRC    , false}, // 54
    {7, -5, 32, 32 * shapeCost,  8 * markerCost, 100 * shapeCost, 0.999, RipUpMode::DRC    , false}, // 55
    {7, -6, 32, 32 * shapeCost,  8 * markerCost, 100 * shapeCost, 0.999, RipUpMode::DRC    , false}, // 56
    {3, -1,  8,      shapeCost,      markerCost, 100 * shapeCost, 0.999, RipUpMode::ALL    , false}, // 57
    {7,  0, 64,      shapeCost,      markerCost, 100 * shapeCost, 0.999, RipUpMode::NEARDRC, false}, // 58
    {7, -1, 64, 64 * shapeCost, 16 * markerCost, 100 * shapeCost, 0.999, RipUpMode::DRC    , false}, // 59
    {7, -2, 64, 64 * shapeCost, 16 * markerCost, 100 * shapeCost, 0.999, RipUpMode::DRC    , false}, // 60
    {7, -3, 64, 64 * shapeCost, 16 * markerCost, 100 * shapeCost, 0.999, RipUpMode::DRC    , false}, // 61
    {7, -4, 64, 64 * shapeCost, 16 * markerCost, 100 * shapeCost, 0.999, RipUpMode::DRC    , false}, // 62
    {7, -5, 64, 64 * shapeCost, 16 * markerCost, 100 * shapeCost, 0.999, RipUpMode::DRC    , false}, // 63
    {7, -6, 64, 64 * shapeCost, 16 * markerCost, 100 * shapeCost, 0.999, RipUpMode::DRC    , false}  // 64
  };
  // clang-format on
}

void addRectToPolySet(gtl::polygon_90_set_data<frCoord>& polySet,
                      odb::Rect rect)
{
  gtl::polygon_90_data<frCoord> poly;
  std::vector<gtl::point_data<frCoord>> points;
  for (const auto& point : rect.getPoints()) {
    points.emplace_back(point.x(), point.y());
  }
  poly.set(points.begin(), points.end());
  using boost::polygon::operators::operator+=;
  polySet += poly;
}

void FlexDR::reportGuideCoverage()
{
  using boost::polygon::operators::operator&;

  const auto numLayers = getTech()->getLayers().size();
  std::vector<uint64_t> totalAreaByLayerNum(numLayers, 0);
  std::vector<uint64_t> totalCoveredAreaByLayerNum(numLayers, 0);
  std::map<frNet*, std::vector<float>> netsCoverage;
  const auto& nets = getDesign()->getTopBlock()->getNets();
  omp_set_num_threads(router_cfg_->MAX_THREADS);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < nets.size(); i++) {  // NOLINT
    const auto& net = nets.at(i);
    std::vector<gtl::polygon_90_set_data<frCoord>> routeSetByLayerNum(
        numLayers),
        guideSetByLayerNum(numLayers);
    for (const auto& shape : net->getShapes()) {
      odb::Rect rect = shape->getBBox();
      addRectToPolySet(routeSetByLayerNum[shape->getLayerNum()], rect);
    }
    for (const auto& pwire : net->getPatchWires()) {
      odb::Rect rect = pwire->getBBox();
      addRectToPolySet(routeSetByLayerNum[pwire->getLayerNum()], rect);
    }
    for (const auto& via : net->getVias()) {
      {
        odb::Rect rect = via->getLayer1BBox();
        addRectToPolySet(routeSetByLayerNum[via->getViaDef()->getLayer1Num()],
                         rect);
      }
      {
        odb::Rect rect = via->getLayer2BBox();
        addRectToPolySet(routeSetByLayerNum[via->getViaDef()->getLayer2Num()],
                         rect);
      }
    }

    for (const auto& shape : net->getOrigGuides()) {
      odb::Rect rect = shape.getBBox();
      addRectToPolySet(guideSetByLayerNum[shape.getLayerNum()], rect);
    }

    for (frLayerNum lNum = 0; lNum < numLayers; lNum++) {
      if (getTech()->getLayer(lNum)->getType() != dbTechLayerType::ROUTING
          || lNum > router_cfg_->TOP_ROUTING_LAYER) {
        continue;
      }
      float coveredPercentage = -1.0;
      uint64_t routingArea = 0;
      uint64_t coveredArea = 0;
      if (!routeSetByLayerNum[lNum].empty()) {
        routingArea = gtl::area(routeSetByLayerNum[lNum]);
        coveredArea
            = gtl::area(routeSetByLayerNum[lNum] & guideSetByLayerNum[lNum]);
        if (routingArea == 0) {
          coveredPercentage = -1.0;
        } else {
          coveredPercentage = (coveredArea / (double) routingArea) * 100;
        }
      }

#pragma omp critical
      {
        netsCoverage[net.get()].push_back(coveredPercentage);
        totalAreaByLayerNum[lNum] += routingArea;
        totalCoveredAreaByLayerNum[lNum] += coveredArea;
      }
    }
  }

  std::ofstream file(router_cfg_->GUIDE_REPORT_FILE);
  file << "Net,";
  for (const auto& layer : getTech()->getLayers()) {
    if (layer->getType() == dbTechLayerType::ROUTING
        && layer->getLayerNum() <= router_cfg_->TOP_ROUTING_LAYER) {
      file << layer->getName() << ",";
    }
  }
  file << std::endl;
  for (const auto& [net, coverage] : netsCoverage) {
    file << net->getName() << ",";
    for (auto coveredPercentage : coverage) {
      if (coveredPercentage < 0.0) {
        file << "NA,";
      } else {
        file << fmt::format("{:.2f}%,", coveredPercentage);
      }
    }
    file << std::endl;
  }
  file << "Total,";
  uint64_t totalArea = 0;
  uint64_t totalCoveredArea = 0;
  for (const auto& layer : getTech()->getLayers()) {
    if (layer->getType() == dbTechLayerType::ROUTING
        && layer->getLayerNum() <= router_cfg_->TOP_ROUTING_LAYER) {
      if (totalAreaByLayerNum[layer->getLayerNum()] == 0) {
        file << "NA,";
        continue;
      }
      totalArea += totalAreaByLayerNum[layer->getLayerNum()];
      totalCoveredArea += totalCoveredAreaByLayerNum[layer->getLayerNum()];
      auto coveredPercentage
          = (totalCoveredAreaByLayerNum[layer->getLayerNum()]
             / (double) totalAreaByLayerNum[layer->getLayerNum()])
            * 100;
      file << fmt::format("{:.2f}%,", coveredPercentage);
    }
  }
  if (totalArea == 0) {
    file << "NA";
  } else {
    auto totalCoveredPercentage = (totalCoveredArea / (double) totalArea) * 100;
    file << fmt::format("{:.2f}%,", totalCoveredPercentage);
  }
  file.close();
}
void FlexDR::fixMaxSpacing()
{
  logger_->info(DRT, 227, "Checking For LEF58_MAXSPACING violations");
  io::Parser parser(db_, getDesign(), logger_, router_cfg_);
  parser.initSecondaryVias();
  std::vector<frVia*> lonely_vias;
  for (const auto& layer : getTech()->getLayers()) {
    if (layer->getType() == odb::dbTechLayerType::CUT) {
      if (!layer->hasLef58MaxSpacingConstraints()) {
        continue;
      }
      for (const auto& rule : layer->getLef58MaxSpacingConstraints()) {
        auto result = getLonelyVias(
            layer.get(), rule->getMaxSpacing(), rule->getCutClassIdx());
        lonely_vias.insert(lonely_vias.end(), result.begin(), result.end());
        for (const auto via_def : layer->getViaDefs()) {
          if (via_def->getCutClassIdx() == rule->getCutClassIdx()) {
            router_->addAvoidViaDefPA(via_def);
          }
        }
      }
    }
  }
  std::vector<odb::Rect> lonely_vias_regions;
  lonely_vias_regions.reserve(lonely_vias.size());
  for (const auto& via : lonely_vias) {
    // Create LEF58_MAXSPACING Markers for the lonely vias
    auto marker = std::make_unique<frMarker>();
    marker->setBBox(via->getBBox());
    auto layer = getTech()->getLayer(via->getViaDef()->getCutLayerNum());
    marker->setLayerNum(layer->getLayerNum());
    marker->setConstraint(layer->getLef58MaxSpacingConstraints().at(0));
    marker->addSrc(via->getNet());
    marker->addVictim(
        via->getNet(),
        std::make_tuple(layer->getLayerNum(), via->getBBox(), false));
    marker->addAggressor(
        via->getNet(),
        std::make_tuple(layer->getLayerNum(), via->getBBox(), false));
    getRegionQuery()->addMarker(marker.get());
    getDesign()->getTopBlock()->addMarker(std::move(marker));
    // Define the lonely via region needed for drWorker creation
    // The region is of size 3x3 GCELLBOX with the via in the center GCELLBOX
    auto origin = via->getOrigin();
    auto block = getDesign()->getTopBlock();
    auto gcell_idx = block->getGCellIdx(origin);
    odb::Rect region, tmp_box;
    tmp_box = block->getGCellBox({gcell_idx.x() - 1, gcell_idx.y() - 1});
    region.set_xlo(tmp_box.xMin());
    region.set_ylo(tmp_box.yMin());
    tmp_box = block->getGCellBox({gcell_idx.x() + 1, gcell_idx.y() + 1});
    region.set_xhi(tmp_box.xMax());
    region.set_yhi(tmp_box.yMax());
    lonely_vias_regions.emplace_back(region);
    if (via->isBottomConnected() || via->isTopConnected()) {
      // get pins connected to the via
      frRegionQuery::Objects<frBlockObject> result;
      getRegionQuery()->query(
          via->isTopConnected() ? via->getLayer2BBox() : via->getLayer1BBox(),
          via->isTopConnected() ? via->getViaDef()->getLayer2Num()
                                : via->getViaDef()->getLayer1Num(),
          result);
      for (auto& [bx, obj] : result) {
        if (obj->typeId() == frcInstTerm) {
          auto inst_term = static_cast<frInstTerm*>(obj);
          if (inst_term->getNet() != via->getNet()) {
            continue;
          }
          inst_term->setStubborn(true);
          router_->addInstancePAData(inst_term->getInst());
        }
      }
    }
  }
  router_->updateDirtyPAData();
  // merge intersecting regions
  std::sort(lonely_vias_regions.begin(),
            lonely_vias_regions.end(),
            [](const odb::Rect& a, const odb::Rect& b) {
              return a.xMin() < b.xMin();
            });
  std::vector<odb::Rect> merged_regions;
  for (const auto& region : lonely_vias_regions) {
    bool found = false;
    for (auto& merged_region : merged_regions) {
      if (region.intersects(merged_region)) {
        merged_region.merge(region);
        found = true;
        break;
      }
    }
    if (!found) {
      merged_regions.emplace_back(region);
    }
  }
  // create drWorkers for the final regions
  omp_set_num_threads(router_cfg_->MAX_THREADS);
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < merged_regions.size(); i++) {
    auto route_box = merged_regions.at(i);
    auto worker = std::make_unique<FlexDRWorker>(
        &via_data_, design_, logger_, router_cfg_);
    odb::Rect ext_box;
    odb::Rect drc_box;
    route_box.bloat(router_cfg_->MTSAFEDIST, ext_box);
    route_box.bloat(router_cfg_->DRCSAFEDIST, drc_box);
    worker->setRouteBox(route_box);
    worker->setExtBox(ext_box);
    worker->setDrcBox(drc_box);
    worker->setDRIter(64);
    worker->setDebugSettings(router_->getDebugSettings());
    worker->setRipupMode(RipUpMode::VIASWAP);
    worker->setGraphics(graphics_.get());
    worker->main(getDesign());
#pragma omp critical
    {
      worker->end(getDesign());
    }
  }
  flow_state_machine_->setFixingMaxSpacing(true);
}

std::vector<frVia*> FlexDR::getLonelyVias(frLayer* layer,
                                          int max_spc,
                                          int cut_class)
{
  std::vector<frVia*> lonely_vias;
  if (layer->getSecondaryViaDefs().empty()) {
    return lonely_vias;
  }
  auto vias = getRegionQuery()->getVias(layer->getLayerNum());
  std::vector<odb::Point> via_positions;
  via_positions.reserve(vias.size());
  for (auto [obj, box] : vias) {
    via_positions.emplace_back(box.xCenter(), box.yCenter());
  }
  KDTree tree(via_positions);
  std::vector<std::atomic_bool> visited(via_positions.size());
  std::fill(visited.begin(), visited.end(), false);
  std::set<int> isolated_via_nodes;
  omp_set_num_threads(router_cfg_->MAX_THREADS);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < via_positions.size(); i++) {
    if (visited[i].load()) {
      continue;
    }
    visited[i].store(true);
    std::vector<int> neighbors = tree.radiusSearch(via_positions[i], max_spc);
    // Check if there are neighbors other than the point itself
    bool is_isolated = true;
    for (const auto& neighbor : neighbors) {
      if (neighbor != i) {
        visited[neighbor].store(true);
        is_isolated = false;
      }
    }
    if (is_isolated) {
#pragma omp critical
      {
        isolated_via_nodes.insert(i);
      }
    }
  }
  for (auto i : isolated_via_nodes) {
    auto obj = vias[i].first;
    if (obj->typeId() != frcVia) {
      continue;
    }
    auto via = static_cast<frVia*>(obj);
    auto net = via->getNet();
    if (net == nullptr || net->isFake() || net->isSpecial()) {
      continue;
    }
    if (via->getViaDef()->getCutClassIdx() != cut_class) {
      continue;
    }
    via->setIsLonely(true);
    lonely_vias.emplace_back(via);
  }
  return lonely_vias;
}

int FlexDR::main()
{
  ProfileTask profile("DR:main");
  auto reporter = logger_->progress()->startIterationReporting(
      "detailed routing", std::min(64, router_cfg_->END_ITERATION), {});

  init();
  incremental_changed_regions_ = router_->getIncrementalChangedRegions();
  frTime t;
  // If the loaded design already carries DRC markers (read from the database in
  // io::Parser::setMarkers), resume from the configured iteration, skipping the
  // initial full rip-up passes that are not useful when re-routing an existing
  // solution's violations.
  const int num_db_markers = getDesign()->getTopBlock()->getNumMarkers();
  if (num_db_markers > 0) {
    // Never resume past the configured end iteration, otherwise no repair
    // iteration would run at all.
    iter_ = std::max(0,
                     std::min(router_cfg_->REROUTE_VIOLATIONS_START_ITER,
                              router_cfg_->END_ITERATION));
    // gcell2BoundaryPin_ is only consumed by iteration 0 and is normally freed
    // at the end of that iteration. When resuming past iteration 0 that cleanup
    // never runs, so release it now to avoid holding the memory for the whole
    // routing.
    if (iter_ > 0) {
      removeGCell2BoundaryPin();
    }
    if (router_cfg_->VERBOSE > 0) {
      logger_->info(DRT,
                    627,
                    "Found {} DRC marker(s) in the database; resuming detailed "
                    "routing from iteration {}.",
                    num_db_markers,
                    iter_);
    }
  }
  bool incremental = false;
  bool hasFixed = false;
  for (const auto& net : getDesign()->getTopBlock()->getNets()) {
    incremental |= net->hasInitialRouting();
    hasFixed |= net->isFixed();
    if (incremental && hasFixed) {
      break;
    }
  }
  if (initial_routing_cleanup_.normalize_initial_routing && incremental
      && iter_ == 0 && num_db_markers == 0) {
    normalizeInitialRouting();
  }
  const auto search_repair_args
      = strategy(router_cfg_->ROUTESHAPECOST, router_cfg_->MARKERCOST);
  // iter_ may start > 0 when resuming from imported DRC markers; index the
  // strategy table by iter_ so the per-iteration parameters stay aligned.
  while (iter_ >= 0 && iter_ < static_cast<int>(search_repair_args.size())) {
    if (iter_ > router_cfg_->END_ITERATION) {
      break;
    }
    auto args = search_repair_args[iter_];
    // Clean-up of initial routing: a worker that has rerouted a net this
    // often without clearing its violations ends, and the next iteration
    // continues with shifted tiles.  The strategy's own value returns with
    // iteration 17, its first full rip-up iteration after the clean-up.
    if (incremental && initial_routing_cleanup_.reroutes_per_net > 0
        && iter_ >= 3 && iter_ <= 16) {
      args.mazeEndIter = std::min(args.mazeEndIter,
                                  initial_routing_cleanup_.reroutes_per_net);
    }
    int clipSize = args.size;
    if (args.ripupMode != RipUpMode::ALL) {
      if (increaseClipsize_) {
        clipSizeInc_ += 2;
      } else {
        clipSizeInc_ = std::max((float) 0, clipSizeInc_ - 0.2f);
      }
      clipSize += std::min(router_cfg_->MAX_CLIPSIZE_INCREASE,
                           (int) round(clipSizeInc_));
    }
    args.size = clipSize;
    if (args.ripupMode == RipUpMode::ALL) {
      if (hasFixed || (incremental && iter_ <= 2)) {
        args.ripupMode = RipUpMode::INCR;
      }
    }
    searchRepair(args);
    if (getDesign()->getTopBlock()->getNumMarkers() == 0) {
      break;
    }
    if (logger_->debugCheck(DRT, "snapshot", 1)) {
      io::Writer writer(getDesign(), logger_);
      writer.updateDb(db_, router_cfg_, false, true);

      std::string snapshotPath = fmt::format(
          "{}/drt_iter{}.odb", router_->getDebugSettings()->snapshotDir, iter_);
      db_->write(utl::OutStreamHandler(snapshotPath.c_str(), true).getStream());
    }
    if (reporter->incrementProgress()) {
      break;
    }
    // Iterations 1 and 2 keep every net with initial routing fixed
    // (incremental mode); after the configured number of them, the repair
    // that rips up such nets at violations (iteration 3) follows.
    if (incremental && iter_ < 2
        && iter_ == initial_routing_cleanup_.incremental_iterations) {
      if (router_cfg_->VERBOSE > 0) {
        logger_->info(DRT,
                      9120,
                      "Initial routing clean-up continues with iteration 3 "
                      "after iteration {}.",
                      iter_);
      }
      iter_ = 2;
    }
    ++iter_;
  }

  end(/* done */ true);
  reporter->end(true);

  if (!router_cfg_->GUIDE_REPORT_FILE.empty()) {
    reportGuideCoverage();
  }
  if (router_cfg_->VERBOSE > 0) {
    t.print(logger_);
    std::cout << std::endl;
  }
  return 0;
}

void FlexDR::sendWorkers(
    const std::vector<std::pair<int, FlexDRWorker*>>& remote_batch,
    std::vector<std::unique_ptr<FlexDRWorker>>& batch)
{
  if (remote_batch.empty()) {
    return;
  }
  std::vector<std::pair<int, std::string>> workers;
  {
    ProfileTask task("DIST: SERIALIZE_BATCH");
    for (auto& [idx, worker] : remote_batch) {
      std::string workerStr;
      serializeWorker(worker, workerStr);
      workers.emplace_back(idx, workerStr);
    }
  }
  std::string remote_ip = dist_ip_;
  uint16_t remote_port = dist_port_;
  if (router_->getCloudSize() > 1) {
    dst::JobMessage msg(dst::JobMessage::kBalancer),
        result(dst::JobMessage::kNone);
    bool ok = dist_->sendJob(msg, dist_ip_.c_str(), dist_port_, result);
    if (!ok) {
      logger_->error(utl::DRT, 7461, "Balancer failed");
    } else {
      dst::BalancerJobDescription* desc
          = static_cast<dst::BalancerJobDescription*>(
              result.getJobDescription());
      remote_ip = desc->getWorkerIP();
      remote_port = desc->getWorkerPort();
    }
  }
  {
    dst::JobMessage msg(dst::JobMessage::kRouting),
        result(dst::JobMessage::kNone);
    std::unique_ptr<dst::JobDescription> desc
        = std::make_unique<RoutingJobDescription>();
    RoutingJobDescription* rjd
        = static_cast<RoutingJobDescription*>(desc.get());
    rjd->setWorkers(workers);
    rjd->setSharedDir(dist_dir_);
    rjd->setSendEvery(20);
    msg.setJobDescription(std::move(desc));
    ProfileTask task("DIST: SENDJOB");
    bool ok = dist_->sendJobMultiResult(
        msg, remote_ip.c_str(), remote_port, result);
    if (!ok) {
      logger_->error(utl::DRT, 500, "Sending worker {} failed");
    }
    for (const auto& one_desc : result.getAllJobDescriptions()) {
      RoutingJobDescription* result_desc
          = static_cast<RoutingJobDescription*>(one_desc.get());
      router_->addWorkerResults(result_desc->getWorkers());
    }
  }
}

bool FlexDR::SearchRepairArgs::isEqualIgnoringSizeAndOffset(
    const SearchRepairArgs& other) const
{
  return (mazeEndIter == other.mazeEndIter
          && workerDRCCost == other.workerDRCCost
          && workerMarkerCost == other.workerMarkerCost
          && workerFixedShapeCost == other.workerFixedShapeCost
          && std::fabs(workerMarkerDecay - other.workerMarkerDecay) < 1e-6
          && ripupMode == other.ripupMode && followGuide == other.followGuide);
}

template <class Archive>
void FlexDRWorker::serialize(Archive& ar, const unsigned int version)
{
  // // We always serialize before calling main on the work unit so various
  // // fields are empty and don't need to be serialized.  I skip these to
  // // save having to write lots of serializers that will never be called.
  // if (!apSVia_.empty() || !nets_.empty() || !owner2nets_.empty()
  //     || !rq_.isEmpty() || gcWorker_) {
  //   logger_->error(DRT, 999, "Can't serialize used worker");
  // }

  // The logger_, graphics_ and debugSettings_ are handled by the caller to
  // use the current ones.
  if (is_loading(ar)) {
    frDesign* design = ar.getDesign();
    design_ = design;
  }
  (ar) & routeBox_;
  (ar) & extBox_;
  (ar) & drcBox_;
  (ar) & drIter_;
  (ar) & mazeEndIter_;
  (ar) & followGuide_;
  (ar) & needRecheck_;
  (ar) & skipRouting_;
  (ar) & ripupMode_;
  (ar) & workerDRCCost_;
  (ar) & workerMarkerCost_;
  (ar) & workerFixedShapeCost_;
  (ar) & workerMarkerDecay_;
  (ar) & initNumMarkers_;
  (ar) & nets_;
  (ar) & markers_;
  (ar) & bestMarkers_;
  (ar) & isCongested_;
  if (is_loading(ar)) {
    // boundaryPin_
    int sz = 0;
    (ar) & sz;
    while (sz--) {
      frBlockObject* obj;
      serializeBlockObject(ar, obj);
      std::set<std::pair<odb::Point, frLayerNum>> val;
      (ar) & val;
      boundaryPin_[(frNet*) obj] = std::move(val);
    }
    // owner2nets_
    for (auto& net : nets_) {
      owner2nets_[net->getFrNet()].push_back(net.get());
    }
    dist_on_ = true;
  } else {
    // boundaryPin_
    int sz = boundaryPin_.size();
    (ar) & sz;
    for (auto& [net, value] : boundaryPin_) {
      frBlockObject* obj = (frBlockObject*) net;
      serializeBlockObject(ar, obj);
      (ar) & value;
    }
  }
}

std::unique_ptr<FlexDRWorker> FlexDRWorker::load(
    const std::string& workerStr,
    FlexDRViaData* via_data,
    frDesign* design,
    utl::Logger* logger,
    RouterConfiguration* router_cfg)
{
  auto worker
      = std::make_unique<FlexDRWorker>(via_data, design, logger, router_cfg);
  deserializeWorker(worker.get(), design, workerStr);
  return worker;
}

// Explicit instantiations
template void FlexDRWorker::serialize<frIArchive>(
    frIArchive& ar,
    const unsigned int file_version);

template void FlexDRWorker::serialize<frOArchive>(
    frOArchive& ar,
    const unsigned int file_version);

}  // namespace drt
