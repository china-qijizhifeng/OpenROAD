// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#include "initialPlace.h"

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "AbstractGraphics.h"
#include "odb/dbTypes.h"
#include "placerBase.h"
#include "solver.h"
#include "utl/Logger.h"

namespace gpl {

using T = Eigen::Triplet<float>;

InitialPlaceVars::InitialPlaceVars(const PlaceOptions& options,
                                   const bool debug)
    : maxIter(options.initialPlaceMaxIter),
      minDiffLength(options.initialPlaceMinDiffLength),
      maxSolverIter(options.initialPlaceMaxSolverIter),
      maxFanout(options.initialPlaceMaxFanout),
      netWeightScale(options.initialPlaceNetWeightScale),
      debug(debug),
      forceCenter(options.forceCenterInitialPlace),
      placeIosMode(options.placeIosMode)
{
}

InitialPlace::InitialPlace(InitialPlaceVars ipVars,
                           std::shared_ptr<PlacerBaseCommon> pbc,
                           std::vector<std::shared_ptr<PlacerBase>>& pbVec,
                           std::unique_ptr<AbstractGraphics> graphics,
                           utl::Logger* log)
    : ipVars_(ipVars),
      pbc_(std::move(pbc)),
      pbVec_(pbVec),
      graphics_(std::move(graphics)),
      log_(log)
{
}

void InitialPlace::doBicgstabPlace(int threads)
{
  log_->info(utl::GPL, 5, "---- Execute Conjugate Gradient Initial Placement.");
  ResidualError error;

  if (graphics_) {
    graphics_->setDebugOn(ipVars_.debug);
  }
  const bool graphics_enabled = graphics_ && graphics_->enabled();
  if (graphics_enabled) {
    graphics_->debugForInitialPlace(pbc_, pbVec_);
  }

  placeInstsInitialPositions();

  // set ExtId for idx reference // easy recovery
  setPlaceInstExtId();

  if (graphics_enabled) {
    gif_key_ = graphics_->gifStart("initPlacement.gif");
  }

  for (size_t iter = 1; iter <= ipVars_.maxIter; iter++) {
    updatePinInfo();
    createSparseMatrix();
    error = cpuSparseSolve(ipVars_.maxSolverIter,
                           iter,
                           placeInstForceMatrixX_,
                           fixedInstForceVecX_,
                           instLocVecX_,
                           placeInstForceMatrixY_,
                           fixedInstForceVecY_,
                           instLocVecY_,
                           log_,
                           threads);

    if (graphics_enabled) {
      graphics_->cellPlot(true);

      odb::Rect region;
      odb::Rect bbox = pbc_->db()->getChip()->getBlock()->getBBox()->getBox();
      int max_dim = std::max(bbox.dx(), bbox.dy());
      double dbu_per_pixel = static_cast<double>(max_dim) / 1000.0;
      graphics_->gifAddFrame(gif_key_, region, 500, dbu_per_pixel, 20);
    }

    if (std::isnan(error.x) || std::isnan(error.y)) {
      log_->warn(utl::GPL,
                 325,
                 "Conjugate gradient initial placement solver failed at "
                 "iteration {}. ",
                 iter);
      break;
    }

    float error_max = std::max(error.x, error.y);
    log_->report(
        "[InitialPlace]  Iter: {} conjugate gradient residual: {:0.8f} HPWL: "
        "{}",
        iter,
        error_max,
        pbc_->getHpwl());
    updateCoordi();

    if (error_max <= 1e-5 && iter >= 5) {
      break;
    }
  }

  if (graphics_enabled) {
    graphics_->gifEnd(gif_key_);
    graphics_->setDebugOn(false);
    graphics_->cellPlot(false);
  }
}

// starting point of initial place is center.
void InitialPlace::placeInstsInitialPositions()
{
  const int core_center_x = pbc_->getDie().coreCx();
  const int core_center_y = pbc_->getDie().coreCy();

  int count_region_center = 0;
  int count_db_location = 0;
  int count_core_center = 0;

  for (auto inst : pbc_->placeInsts()) {
    if (inst->isLocked()) {
      continue;
    }

    const auto db_inst = inst->dbInst();
    const auto group = db_inst->getGroup();

    if (group && group->getRegion()) {
      auto region = group->getRegion();
      int region_x_min = std::numeric_limits<int>::max();
      int region_y_min = std::numeric_limits<int>::max();
      int region_x_max = std::numeric_limits<int>::min();
      int region_y_max = std::numeric_limits<int>::min();

      for (auto boundary : region->getBoundaries()) {
        region_x_min = std::min(region_x_min, boundary->xMin());
        region_y_min = std::min(region_y_min, boundary->yMin());
        region_x_max = std::max(region_x_max, boundary->xMax());
        region_y_max = std::max(region_y_max, boundary->yMax());
      }

      inst->setCenterLocation(region_x_max - (region_x_max - region_x_min) / 2,
                              region_y_max - (region_y_max - region_y_min) / 2);
      ++count_region_center;
    } else if (!ipVars_.forceCenter && db_inst->isPlaced()) {
      const auto bbox = db_inst->getBBox()->getBox();
      inst->setCenterLocation(bbox.xCenter(), bbox.yCenter());
      ++count_db_location;
    } else {
      inst->setCenterLocation(core_center_x, core_center_y);
      ++count_core_center;
    }
  }

  log_->info(utl::GPL,
             51,
             "Source of initial instance position counters:\n"
             "\tOdb location = {}"
             "\tCore center = {}"
             "\tRegion center = {}",
             count_db_location,
             count_core_center,
             count_region_center);
}

void InitialPlace::setPlaceInstExtId()
{
  // reset ExtId for all instances
  for (auto inst : pbc_->getInsts()) {
    inst->setExtId(INT_MAX);
  }
  // set index only with place-able instances
  for (size_t i = 0; i < pbc_->placeInsts().size(); ++i) {
    pbc_->placeInsts()[i]->setExtId(i);
  }
}

// Only movable IO port pins are excluded from the initial B2B solve;
// fixed ports (LOCKED/FIRM/COVER) keep anchoring cells here.
static bool isMovableIoPin(SoaRef<Pin> pin, bool place_ios_mode)
{
  return place_ios_mode && pin->isBTerm()
         && !pin->getDbBTerm()->getFirstPinPlacementStatus().isFixed();
}

void InitialPlace::updatePinInfo()
{
  auto& pins = *pbc_->pinColumns();
  auto& min_x = pins.column<PinColumns::minPinXField>();
  auto& min_y = pins.column<PinColumns::minPinYField>();
  auto& max_x = pins.column<PinColumns::maxPinXField>();
  auto& max_y = pins.column<PinColumns::maxPinYField>();
  std::fill(min_x.begin(), min_x.end(), 0);
  std::fill(min_y.begin(), min_y.end(), 0);
  std::fill(max_x.begin(), max_x.end(), 0);
  std::fill(max_y.begin(), max_y.end(), 0);
  const auto& x = pins.column<PinColumns::cx>();
  const auto& y = pins.column<PinColumns::cy>();
  const size_t none = x.size();
  triplet_capacity_x_ = triplet_capacity_y_ = pbc_->placeInsts().size();
  solve_pin_capacity_ = 0;
  for (const auto net : pbc_->getNets()) {
    size_t xmin = none, xmax = none, ymin = none, ymax = none, count = 0;
    int lx = INT_MAX, ly = INT_MAX, ux = INT_MIN, uy = INT_MIN;
    const auto& adjacency = net->getPins();
    for (size_t row : adjacency.rows()) {
      if (isMovableIoPin({adjacency.storage(), row}, ipVars_.placeIosMode)) {
        continue;
      }
      ++count;
      // Strict comparisons keep the first pin on equal coordinates.
      if (lx > x[row]) {
        lx = x[row];
        xmin = row;
      }
      if (ux < x[row]) {
        ux = x[row];
        xmax = row;
      }
      if (ly > y[row]) {
        ly = y[row];
        ymin = row;
      }
      if (uy < y[row]) {
        uy = y[row];
        ymax = row;
      }
    }
    if (xmin != none) {
      min_x[xmin] = 1;
    }
    if (xmax != none) {
      max_x[xmax] = 1;
    }
    if (ymin != none) {
      min_y[ymin] = 1;
    }
    if (ymax != none) {
      max_y[ymax] = 1;
    }
    if (count > 1 && adjacency.size() < ipVars_.maxFanout) {
      solve_pin_capacity_ = std::max(solve_pin_capacity_, count);
      // Each B2B pair contributes at most four triplets. Fixed and same-inst
      // pairs contribute fewer; this bound avoids growth without an O(n^2)
      // pass.
      auto bound = [count, none](size_t first, size_t last) {
        const size_t extrema
            = (first != none) + (last != none && last != first);
        return 4 * (extrema * (count - extrema) + extrema * (extrema - 1) / 2);
      };
      triplet_capacity_x_ += bound(xmin, xmax);
      triplet_capacity_y_ += bound(ymin, ymax);
    }
  }
}

// solve placeInstForceMatrixX_ * xcg_x_ = xcg_b_ and placeInstForceMatrixY_ *
// ycg_x_ = ycg_b_ eq.
void InitialPlace::createSparseMatrix()
{
  const int placeCnt = pbc_->placeInsts().size();
  instLocVecX_.resize(placeCnt);
  fixedInstForceVecX_.resize(placeCnt);
  instLocVecY_.resize(placeCnt);
  fixedInstForceVecY_.resize(placeCnt);

  placeInstForceMatrixX_.resize(placeCnt, placeCnt);
  placeInstForceMatrixY_.resize(placeCnt, placeCnt);

  //
  // listX and listY is a temporary vector that have tuples, (idx1, idx2, val)
  //
  // listX finally becomes placeInstForceMatrixX_
  // listY finally becomes placeInstForceMatrixY_
  //
  // The triplet vector is recommended usages
  // to fill in SparseMatrix from Eigen docs.
  //

  std::vector<T> listX, listY;
  listX.reserve(triplet_capacity_x_);
  listY.reserve(triplet_capacity_y_);

  // initialize vector
  for (auto inst : pbc_->placeInsts()) {
    int idx = inst->getExtId();

    instLocVecX_(idx) = inst->cx();
    instLocVecY_(idx) = inst->cy();

    fixedInstForceVecX_(idx) = fixedInstForceVecY_(idx) = 0;

    if (inst->isLocked()) {
      listX.emplace_back(idx, idx, 1.0f);
      listY.emplace_back(idx, idx, 1.0f);
      fixedInstForceVecX_(idx) = inst->cx();
      fixedInstForceVecY_(idx) = inst->cy();
    }
  }

  struct SolvePin
  {
    SoaRef<Instance> instance;
    int x, y, inst_x, inst_y, index;
    bool movable, locked, min_x, max_x, min_y, max_y;
  };
  std::vector<SolvePin> solve_pins;
  solve_pins.reserve(solve_pin_capacity_);
  const auto& pins = *pbc_->pinColumns();
  const auto& instances = *pbc_->instanceColumns();
  // for each net
  for (auto net : pbc_->getNets()) {
    // skip for small nets.
    if (net->getPins().size() <= 1) {
      continue;
    }

    // escape long time cals on huge fanout.
    //
    if (net->getPins().size() >= ipVars_.maxFanout) {
      continue;
    }

    float netWeight = ipVars_.netWeightScale / (net->getPins().size() - 1);
    solve_pins.clear();
    for (SoaRef<Pin> pin : net->getPins()) {
      if (!isMovableIoPin(pin, ipVars_.placeIosMode)) {
        const size_t row = pin.row();
        const auto instance = pins.at<PinColumns::inst>(row);
        const size_t ir = instance ? instance.row() : 0;
        const int inst_x = instance ? (instances.at<InstanceColumns::lx>(ir)
                                       + instances.at<InstanceColumns::ux>(ir))
                                          / 2
                                    : 0;
        const int inst_y = instance ? (instances.at<InstanceColumns::ly>(ir)
                                       + instances.at<InstanceColumns::uy>(ir))
                                          / 2
                                    : 0;
        solve_pins.push_back(
            {instance,
             pins.at<PinColumns::cx>(row),
             pins.at<PinColumns::cy>(row),
             inst_x,
             inst_y,
             instance ? instances.at<InstanceColumns::extId>(ir) : INT_MAX,
             pin->isPlaceInstConnected(),
             instance && instance->isLocked(),
             bool(pins.at<PinColumns::minPinXField>(row)),
             bool(pins.at<PinColumns::maxPinXField>(row)),
             bool(pins.at<PinColumns::minPinYField>(row)),
             bool(pins.at<PinColumns::maxPinYField>(row))});
      }
    }
    // Only pairs touching an x/y extreme can contribute B2B terms. Keep the
    // original lexicographic pin-pair order while skipping all other pairs.
    std::array<int, 4> extrema;
    size_t extreme_count = 0;
    for (size_t i = 0; i < solve_pins.size(); ++i) {
      const auto& pin = solve_pins[i];
      if (pin.min_x || pin.max_x || pin.min_y || pin.max_y) {
        extrema[extreme_count++] = i;
      }
    }
    for (int pinIdx1 = 1; pinIdx1 < solve_pins.size(); ++pinIdx1) {
      const auto& pin1 = solve_pins[pinIdx1];
      const bool all_predecessors
          = pin1.min_x || pin1.max_x || pin1.min_y || pin1.max_y;
      const size_t pair_count
          = all_predecessors
                ? pinIdx1
                : std::lower_bound(
                      extrema.begin(), extrema.begin() + extreme_count, pinIdx1)
                      - extrema.begin();
      for (size_t pair = 0; pair < pair_count; ++pair) {
        const int pinIdx2 = all_predecessors ? pair : extrema[pair];
        const auto& pin2 = solve_pins[pinIdx2];
        // no need to fill in when instance is same
        if (pin1.instance == pin2.instance) {
          continue;
        }

        // B2B modeling on min/maxX pins.
        if (pin1.min_x || pin1.max_x || pin2.min_x || pin2.max_x) {
          int diffX = abs(pin1.x - pin2.x);
          float weightX = 0;
          if (diffX > ipVars_.minDiffLength) {
            weightX = netWeight / diffX;
          } else {
            weightX = netWeight / ipVars_.minDiffLength;
          }

          // both pin cames from instance
          if (pin1.movable && pin2.movable) {
            const int inst1 = pin1.index;
            const int inst2 = pin2.index;
            const bool locked1 = pin1.locked;
            const bool locked2 = pin2.locked;

            if (!locked1 && !locked2) {
              // Both movable: standard B2B matrix terms
              listX.emplace_back(inst1, inst1, weightX);
              listX.emplace_back(inst2, inst2, weightX);
              listX.emplace_back(inst1, inst2, -weightX);
              listX.emplace_back(inst2, inst1, -weightX);

              fixedInstForceVecX_(inst1)
                  += -weightX
                     * ((pin1.x - pin1.inst_x) - (pin2.x - pin2.inst_x));
              fixedInstForceVecX_(inst2)
                  += -weightX
                     * ((pin2.x - pin2.inst_x) - (pin1.x - pin1.inst_x));
            } else if (locked1 && !locked2) {
              // inst1 is fixed: treat pin1 as a fixed terminal for inst2
              listX.emplace_back(inst2, inst2, weightX);
              fixedInstForceVecX_(inst2)
                  += weightX * (pin1.x - (pin2.x - pin2.inst_x));
            } else if (!locked1 && locked2) {
              // inst2 is fixed: treat pin2 as a fixed terminal for inst1
              listX.emplace_back(inst1, inst1, weightX);
              fixedInstForceVecX_(inst1)
                  += weightX * (pin2.x - (pin1.x - pin1.inst_x));
            }
            // both locked: no matrix or RHS contribution needed
          }
          // pin1 from IO port / pin2 from Instance
          else if (!pin1.movable && pin2.movable) {
            const int inst2 = pin2.index;
            listX.emplace_back(inst2, inst2, weightX);

            fixedInstForceVecX_(inst2)
                += weightX * (pin1.x - (pin2.x - pin2.inst_x));
          }
          // pin1 from Instance / pin2 from IO port
          else if (pin1.movable && !pin2.movable) {
            const int inst1 = pin1.index;
            listX.emplace_back(inst1, inst1, weightX);

            fixedInstForceVecX_(inst1)
                += weightX * (pin2.x - (pin1.x - pin1.inst_x));
          }
        }

        // B2B modeling on min/maxY pins.
        if (pin1.min_y || pin1.max_y || pin2.min_y || pin2.max_y) {
          int diffY = abs(pin1.y - pin2.y);
          float weightY = 0;
          if (diffY > ipVars_.minDiffLength) {
            weightY = netWeight / diffY;
          } else {
            weightY = netWeight / ipVars_.minDiffLength;
          }

          // both pin cames from instance
          if (pin1.movable && pin2.movable) {
            const int inst1 = pin1.index;
            const int inst2 = pin2.index;
            const bool locked1 = pin1.locked;
            const bool locked2 = pin2.locked;

            if (!locked1 && !locked2) {
              // Both movable: standard B2B matrix terms
              listY.emplace_back(inst1, inst1, weightY);
              listY.emplace_back(inst2, inst2, weightY);
              listY.emplace_back(inst1, inst2, -weightY);
              listY.emplace_back(inst2, inst1, -weightY);

              fixedInstForceVecY_(inst1)
                  += -weightY
                     * ((pin1.y - pin1.inst_y) - (pin2.y - pin2.inst_y));
              fixedInstForceVecY_(inst2)
                  += -weightY
                     * ((pin2.y - pin2.inst_y) - (pin1.y - pin1.inst_y));
            } else if (locked1 && !locked2) {
              // inst1 is fixed: treat pin1 as a fixed terminal for inst2
              listY.emplace_back(inst2, inst2, weightY);
              fixedInstForceVecY_(inst2)
                  += weightY * (pin1.y - (pin2.y - pin2.inst_y));
            } else if (!locked1 && locked2) {
              // inst2 is fixed: treat pin2 as a fixed terminal for inst1
              listY.emplace_back(inst1, inst1, weightY);
              fixedInstForceVecY_(inst1)
                  += weightY * (pin2.y - (pin1.y - pin1.inst_y));
            }
            // both locked: no matrix or RHS contribution needed
          }
          // pin1 from IO port / pin2 from Instance
          else if (!pin1.movable && pin2.movable) {
            const int inst2 = pin2.index;
            listY.emplace_back(inst2, inst2, weightY);

            fixedInstForceVecY_(inst2)
                += weightY * (pin1.y - (pin2.y - pin2.inst_y));
          }
          // pin1 from Instance / pin2 from IO port
          else if (pin1.movable && !pin2.movable) {
            const int inst1 = pin1.index;
            listY.emplace_back(inst1, inst1, weightY);

            fixedInstForceVecY_(inst1)
                += weightY * (pin2.y - (pin1.y - pin1.inst_y));
          }
        }
      }
    }
  }

  placeInstForceMatrixX_.setFromTriplets(listX.begin(), listX.end());
  placeInstForceMatrixY_.setFromTriplets(listY.begin(), listY.end());
}

void InitialPlace::updateCoordi()
{
  for (auto inst : pbc_->placeInsts()) {
    int idx = inst->getExtId();
    if (!inst->isLocked()) {
      int new_x = instLocVecX_(idx);
      int new_y = instLocVecY_(idx);

      // Constrain to core area
      const auto& die = pbc_->getDie();
      new_x = std::max(new_x, die.coreLx());
      new_x = std::min(new_x, die.coreUx());
      new_y = std::max(new_y, die.coreLy());
      new_y = std::min(new_y, die.coreUy());

      // If instance has a region constraint, use that instead
      const auto db_inst = inst->dbInst();
      const auto group = db_inst->getGroup();
      if (group && group->getRegion()) {
        auto region = group->getRegion();
        int region_x_min = std::numeric_limits<int>::max();
        int region_y_min = std::numeric_limits<int>::max();
        int region_x_max = std::numeric_limits<int>::min();
        int region_y_max = std::numeric_limits<int>::min();

        for (auto boundary : region->getBoundaries()) {
          region_x_min = std::min(region_x_min, boundary->xMin());
          region_y_min = std::min(region_y_min, boundary->yMin());
          region_x_max = std::max(region_x_max, boundary->xMax());
          region_y_max = std::max(region_y_max, boundary->yMax());
        }

        new_x = std::max(new_x, region_x_min);
        new_x = std::min(new_x, region_x_max);
        new_y = std::max(new_y, region_y_min);
        new_y = std::min(new_y, region_y_max);
      }

      inst->dbSetCenterLocation(new_x, new_y);
      inst->dbSetPlaced();
    }
  }
}

}  // namespace gpl
