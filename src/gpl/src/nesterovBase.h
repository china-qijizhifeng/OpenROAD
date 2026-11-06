// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#pragma once
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include "boost/random/mersenne_twister.hpp"
#include "boost/unordered/unordered_flat_map.hpp"
#include "densitySoa.h"
#include "gpl/Replace.h"
#include "hpwlBackend.h"
#include "odb/db.h"
#include "placementColumns.h"
#include "placerBase.h"
#include "point.h"
#include "routeBase.h"
#include "soaStorage.h"
#include "utl/Logger.h"
namespace odb {
class dbInst;
class dbITerm;
class dbBTerm;
class dbNet;
}  // namespace odb

namespace utl {
class Logger;
}

namespace gpl {

class Instance;
class Die;
class PlacerBaseCommon;
class PlacerBase;
class GCellHandle;

class Instance;
class Pin;
class Net;

class GPin;
class FFT;
class nesterovDbCbk;
class DeviceState;         // gpu/deviceState.h (GPU-only, forward decl here)
class RegionDensityField;  // gpu/regionDensityField.h (GPU-only)
class WirelengthGradientBackend;  // wirelengthGradientBackend.h
class DensityGradientBackend;     // densityGradientBackend.h
class NesterovDeviceContext;      // gpu/nesterovDeviceContext.h
class IoPinBase;                  // ioPinBase.h
enum class SlpSlot : int;         // gpu/nesterovDeviceContext.h
enum class SumGradSlot : int;     // gpu/nesterovDeviceContext.h

class GCell
{
 public:
  using NumericStore = GCellColumns::Store;
  using GCellChange = CellAreaChange;

  // instance cells
  GCell(SoaRef<Instance> inst,
        std::shared_ptr<NumericStore> storage
        = std::make_shared<NumericStore>());
  GCell(const SoaRows<Instance>& insts,
        std::shared_ptr<NumericStore> storage
        = std::make_shared<NumericStore>());

  // filler cells
  GCell(int cx,
        int cy,
        int dx,
        int dy,
        std::shared_ptr<NumericStore> storage
        = std::make_shared<NumericStore>());

  const SoaRows<Instance>& insts() const
  {
    return data_.get<GCellColumns::insts>();
  }
  const SoaRows<GPin>& gPins() const
  {
    return data_.get<GCellColumns::gPins>();
  }

  std::string getName() const;

  void addGPin(SoaRef<GPin> gPin);
  void clearGPins() { data_.get<GCellColumns::gPins>().clear(); }

  void updateLocations();

  bool isLocked() const;
  void lock();

  // normal coordinates
  int lx() const;
  int ly() const;
  int ux() const;
  int uy() const;
  int cx() const;
  int cy() const;
  int dx() const;
  int dy() const;

  // virtual density coordinates
  int dLx() const;
  int dLy() const;
  int dUx() const;
  int dUy() const;
  int dCx() const;
  int dCy() const;
  int dDx() const;
  int dDy() const;

  void setCenterLocation(int cx, int cy);
  // void setLocation(int x, int y);
  void setSize(int dx, int dy, GCellChange change = GCellChange::kNone);
  int64_t getInflationArea() const;
  void setAreaChangeType(GCellChange change)
  {
    data_.get<GCellColumns::change>() = change;
  }
  GCellChange changeType() const { return data_.get<GCellColumns::change>(); }
  void setAllLocations(int lx, int ly, int ux, int uy);

  void setDensityLocation(int dLx, int dLy);
  void setDensityCenterLocation(int dCx, int dCy);
  void setDensitySize(int dDx, int dDy);

  void setDensityScale(float densityScale);
  void setGradientX(float gradientX);
  void setGradientY(float gradientY);

  float getGradientX() const { return data_.get<GCellColumns::gradientX>(); }
  float getGradientY() const { return data_.get<GCellColumns::gradientY>(); }
  float getDensityScale() const
  {
    return data_.get<GCellColumns::densityScale>();
  }

  bool isInstance() const;
  bool isFiller() const;
  bool isMacroInstance() const;
  bool isStdInstance() const;
  bool contains(odb::dbInst* db_inst) const;

  void print(utl::Logger* logger, bool print_only_name) const;
  void writeAttributesToCSV(std::ostream& out) const;

  SoaRow<NumericStore>& numericRow() { return data_; }
  const SoaRow<NumericStore>& numericRow() const { return data_; }

  GCell(SoaViewTag tag, void* store, size_t row) : data_(tag, store, row) {}
  SoaRef<GCell> ref() const { return {data_.storageIdentity(), data_.slot()}; }
  using Store = NumericStore;

 private:
  SoaRow<NumericStore> data_;
};

inline int GCell::lx() const
{
  return data_.get<GCellColumns::lx>();
}
inline int GCell::ly() const
{
  return data_.get<GCellColumns::ly>();
}

inline int GCell::ux() const
{
  return data_.get<GCellColumns::ux>();
}

inline int GCell::uy() const
{
  return data_.get<GCellColumns::uy>();
}

inline int GCell::cx() const
{
  return (data_.get<GCellColumns::lx>() + data_.get<GCellColumns::ux>()) / 2;
}

inline int GCell::cy() const
{
  return (data_.get<GCellColumns::ly>() + data_.get<GCellColumns::uy>()) / 2;
}

inline int GCell::dx() const
{
  return data_.get<GCellColumns::ux>() - data_.get<GCellColumns::lx>();
}

inline int GCell::dy() const
{
  return data_.get<GCellColumns::uy>() - data_.get<GCellColumns::ly>();
}

inline int GCell::dLx() const
{
  return data_.get<GCellColumns::dLx>();
}

inline int GCell::dLy() const
{
  return data_.get<GCellColumns::dLy>();
}

inline int GCell::dUx() const
{
  return data_.get<GCellColumns::dUx>();
}

inline int GCell::dUy() const
{
  return data_.get<GCellColumns::dUy>();
}

inline int GCell::dCx() const
{
  return (data_.get<GCellColumns::dUx>() + data_.get<GCellColumns::dLx>()) / 2;
}

inline int GCell::dCy() const
{
  return (data_.get<GCellColumns::dUy>() + data_.get<GCellColumns::dLy>()) / 2;
}

inline int GCell::dDx() const
{
  return data_.get<GCellColumns::dUx>() - data_.get<GCellColumns::dLx>();
}

inline int GCell::dDy() const
{
  return data_.get<GCellColumns::dUy>() - data_.get<GCellColumns::dLy>();
}

class GNet
{
 public:
  using NumericStore = GNetColumns::Store;
  GNet(SoaRef<Net> net,
       std::shared_ptr<NumericStore> storage
       = std::make_shared<NumericStore>());
  GNet(const SoaRows<Net>& nets,
       std::shared_ptr<NumericStore> storage
       = std::make_shared<NumericStore>());

  SoaRef<Net> getPbNet() const;
  const SoaRows<Net>& getPbNets() const
  {
    return data_.get<GNetColumns::nets>();
  }
  const SoaRows<GPin>& getGPins() const
  {
    return data_.get<GNetColumns::gPins>();
  }

  int lx() const;
  int ly() const;
  int ux() const;
  int uy() const;

  void setTimingWeight(float timingWeight);
  void setCustomWeight(float customWeight);
  void setPowerWeight(float powerWeight);

  float getTotalWeight() const
  {
    return (data_.get<GNetColumns::timingWeight>()
            + data_.get<GNetColumns::powerWeight>())
           * data_.get<GNetColumns::customWeight>();
  }
  float getTimingWeight() const
  {
    return data_.get<GNetColumns::timingWeight>();
  }
  float getCustomWeight() const
  {
    return data_.get<GNetColumns::customWeight>();
  }

  void addGPin(SoaRef<GPin> gPin);
  void clearGPins() { data_.get<GNetColumns::gPins>().clear(); }
  void updateBox();
  // GPU path writes the device-computed bbox back so lx()/ly()/ux()/uy() stay
  // consistent with updateBox() without re-walking pins on the host. Caller
  // must pass updateBox()-equivalent values; no validation is done.
  void setBox(int lx, int ly, int ux, int uy);
  int64_t getHpwl() const;

  void setDontCare();
  bool isDontCare() const;

  // clear WA(Weighted Average) variables.
  void clearWaVars();

  void addWaExpMinSumX(float waExpMinSumX);
  void addWaXExpMinSumX(float waXExpMinSumX);

  void addWaExpMinSumY(float waExpMinSumY);
  void addWaYExpMinSumY(float waYExpMinSumY);

  void addWaExpMaxSumX(float waExpMaxSumX);
  void addWaXExpMaxSumX(float waXExpMaxSumX);

  void addWaExpMaxSumY(float waExpMaxSumY);
  void addWaYExpMaxSumY(float waYExpMaxSumY);

  float waExpMinSumX() const;
  float waXExpMinSumX() const;

  float waExpMinSumY() const;
  float waYExpMinSumY() const;

  float waExpMaxSumX() const;
  float waXExpMaxSumX() const;

  float waExpMaxSumY() const;
  float waYExpMaxSumY() const;

  void print(utl::Logger* log) const;

  SoaRow<NumericStore>& numericRow() { return data_; }
  const SoaRow<NumericStore>& numericRow() const { return data_; }

  GNet(SoaViewTag tag, void* store, size_t row) : data_(tag, store, row) {}
  SoaRef<GNet> ref() const { return {data_.storageIdentity(), data_.slot()}; }
  using Store = NumericStore;

 private:
  SoaRow<NumericStore> data_;

  //
  // weighted average WL model stor for better indexing
  // Please check the equation (4) in the ePlace-MS paper.
  //
  // WA: weighted Average
  // saving four variable will be helpful for
  // calculating the WA gradients/wirelengths.
  //
  // gamma: modeling accuracy.
  //
  // X forces.
  //
  // waExpMinSumX_: store sigma {exp(x_i/gamma)}
  // waXExpMinSumX_: store signa {x_i*exp(e_i/gamma)}
  // waExpMaxSumX_ : store sigma {exp(-x_i/gamma)}
  // waXExpMaxSumX_: store sigma {x_i*exp(-x_i/gamma)}
  //

  //
  // Y forces.
  //
  // waExpMinSumY_: store sigma {exp(y_i/gamma)}
  // waYExpMinSumY_: store signa {y_i*exp(e_i/gamma)}
  // waExpMaxSumY_ : store sigma {exp(-y_i/gamma)}
  // waYExpMaxSumY_: store sigma {y_i*exp(-y_i/gamma)}
  //
};

inline int GNet::lx() const
{
  return data_.get<GNetColumns::lx>();
}

inline int GNet::ly() const
{
  return data_.get<GNetColumns::ly>();
}

inline int GNet::ux() const
{
  return data_.get<GNetColumns::ux>();
}

inline int GNet::uy() const
{
  return data_.get<GNetColumns::uy>();
}

// eight add functions
inline void GNet::addWaExpMinSumX(float waExpMinSumX)
{
  data_.get<GNetColumns::waExpMinSumX>() += waExpMinSumX;
}

inline void GNet::addWaXExpMinSumX(float waXExpMinSumX)
{
  data_.get<GNetColumns::waXExpMinSumX>() += waXExpMinSumX;
}

inline void GNet::addWaExpMinSumY(float waExpMinSumY)
{
  data_.get<GNetColumns::waExpMinSumY>() += waExpMinSumY;
}

inline void GNet::addWaYExpMinSumY(float waYExpMinSumY)
{
  data_.get<GNetColumns::waYExpMinSumY>() += waYExpMinSumY;
}

inline void GNet::addWaExpMaxSumX(float waExpMaxSumX)
{
  data_.get<GNetColumns::waExpMaxSumX>() += waExpMaxSumX;
}

inline void GNet::addWaXExpMaxSumX(float waXExpMaxSumX)
{
  data_.get<GNetColumns::waXExpMaxSumX>() += waXExpMaxSumX;
}

inline void GNet::addWaExpMaxSumY(float waExpMaxSumY)
{
  data_.get<GNetColumns::waExpMaxSumY>() += waExpMaxSumY;
}

inline void GNet::addWaYExpMaxSumY(float waYExpMaxSumY)
{
  data_.get<GNetColumns::waYExpMaxSumY>() += waYExpMaxSumY;
}

inline float GNet::waExpMinSumX() const
{
  return data_.get<GNetColumns::waExpMinSumX>();
}

inline float GNet::waXExpMinSumX() const
{
  return data_.get<GNetColumns::waXExpMinSumX>();
}

inline float GNet::waExpMinSumY() const
{
  return data_.get<GNetColumns::waExpMinSumY>();
}

inline float GNet::waYExpMinSumY() const
{
  return data_.get<GNetColumns::waYExpMinSumY>();
}

inline float GNet::waExpMaxSumX() const
{
  return data_.get<GNetColumns::waExpMaxSumX>();
}

inline float GNet::waXExpMaxSumX() const
{
  return data_.get<GNetColumns::waXExpMaxSumX>();
}

inline float GNet::waExpMaxSumY() const
{
  return data_.get<GNetColumns::waExpMaxSumY>();
}

inline float GNet::waYExpMaxSumY() const
{
  return data_.get<GNetColumns::waYExpMaxSumY>();
}

class GPin
{
 public:
  using NumericStore = GPinColumns::Store;
  GPin(SoaRef<Pin> pin,
       std::shared_ptr<NumericStore> storage
       = std::make_shared<NumericStore>());
  GPin(const SoaRows<Pin>& pins,
       std::shared_ptr<NumericStore> storage
       = std::make_shared<NumericStore>());

  SoaRef<Pin> getPbPin() const;
  const SoaRows<Pin>& getPbPins() const
  {
    return data_.get<GPinColumns::pins>();
  }

  SoaRef<GCell> getGCell() const { return data_.get<GPinColumns::gCell>(); }
  SoaRef<GNet> getGNet() const { return data_.get<GPinColumns::gNet>(); }

  void setGCell(SoaRef<GCell> gCell);
  void setGNet(SoaRef<GNet> gNet);

  int cx() const { return data_.get<GPinColumns::cx>(); }
  int cy() const { return data_.get<GPinColumns::cy>(); }

  // Offset from the owning GCell's center. The absolute pin center
  // (cx_/cy_) is recomputed by updateLocation() as gCell->cx() + offsetCx_.
  // Exposed for GPU paths that maintain pin coordinates device-side from
  // inst centers + per-pin offsets (gpu/deviceState.cpp).
  int offsetCx() const { return data_.get<GPinColumns::offsetCx>(); }
  int offsetCy() const { return data_.get<GPinColumns::offsetCy>(); }

  // clear WA(Weighted Average) variables.
  void clearWaVars();

  void setMaxExpSumX(float maxExpSumX);
  void setMaxExpSumY(float maxExpSumY);
  void setMinExpSumX(float minExpSumX);
  void setMinExpSumY(float minExpSumY);

  float maxExpSumX() const
  {
    return data_.storage().at<GPinColumns::maxExpSumX>(
        data_.get<GPinColumns::waSlot>());
  }
  float maxExpSumY() const
  {
    return data_.storage().at<GPinColumns::maxExpSumY>(
        data_.get<GPinColumns::waSlot>());
  }
  float minExpSumX() const
  {
    return data_.storage().at<GPinColumns::minExpSumX>(
        data_.get<GPinColumns::waSlot>());
  }
  float minExpSumY() const
  {
    return data_.storage().at<GPinColumns::minExpSumY>(
        data_.get<GPinColumns::waSlot>());
  }

  bool hasMaxExpSumX() const
  {
    return data_.storage().at<GPinColumns::hasMaxExpSumX>(
        data_.get<GPinColumns::waSlot>());
  }
  bool hasMaxExpSumY() const
  {
    return data_.storage().at<GPinColumns::hasMaxExpSumY>(
        data_.get<GPinColumns::waSlot>());
  }
  bool hasMinExpSumX() const
  {
    return data_.storage().at<GPinColumns::hasMinExpSumX>(
        data_.get<GPinColumns::waSlot>());
  }
  bool hasMinExpSumY() const
  {
    return data_.storage().at<GPinColumns::hasMinExpSumY>(
        data_.get<GPinColumns::waSlot>());
  }

  void setCenterLocation(int cx, int cy);
  void updateLocation(SoaRef<GCell> gCell);
  void updateDensityLocation(SoaRef<GCell> gCell);

  void print(utl::Logger* log) const;
  void updateCoordi();

  SoaRow<NumericStore>& numericRow() { return data_; }
  const SoaRow<NumericStore>& numericRow() const { return data_; }

  GPin(SoaViewTag tag, void* store, size_t row) : data_(tag, store, row) {}
  SoaRef<GPin> ref() const { return {data_.storageIdentity(), data_.slot()}; }
  using Store = NumericStore;

 private:
  SoaRow<NumericStore> data_;

  // weighted average WL vals stor for better indexing
  // Please check the equation (4) in the ePlace-MS paper.
  //
  // maxExpSum_: holds exp(x_i/gamma)
  // minExpSum_: holds exp(-x_i/gamma)
  // the x_i is equal to cx_ variable.
  //

  // flag variables
  //
  // check whether
  // this pin is considered in a WA models.
};

class Bin
{
 public:
  using NumericStore = BinColumns::Store;
  Bin() : data_(std::make_shared<NumericStore>()) {}
  Bin(int x,
      int y,
      int lx,
      int ly,
      int ux,
      int uy,
      float targetDensity,
      std::shared_ptr<NumericStore> storage = std::make_shared<NumericStore>());

  int x() const { return data_.get<BinColumns::x>(); };
  int y() const { return data_.get<BinColumns::y>(); };

  int lx() const { return data_.get<BinColumns::lx>(); };
  int ly() const { return data_.get<BinColumns::ly>(); };
  int ux() const { return data_.get<BinColumns::ux>(); };
  int uy() const { return data_.get<BinColumns::uy>(); };
  int cx() const;
  int cy() const;
  int dx() const;
  int dy() const;

  float electroPhi() const;
  float electroFieldX() const;
  float electroFieldY() const;
  float getTargetDensity() const;
  float getDensity() const;

  void setDensity(float density);
  void setBinTargetDensity(float density);
  void setElectroField(float electroFieldX, float electroFieldY);
  void setElectroPhi(float phi);

  void setNonPlaceArea(int64_t area);
  void setInstPlacedArea(int64_t area);
  void setFillerArea(int64_t area);

  void setNonPlaceAreaUnscaled(int64_t area);
  void setInstPlacedAreaUnscaled(int64_t area);

  void addNonPlaceArea(int64_t area);
  void addInstPlacedArea(int64_t area);
  void addFillerArea(int64_t area);
  void atomicAddFillerArea(int64_t area);

  void addNonPlaceAreaUnscaled(int64_t area);
  void addInstPlacedAreaUnscaled(int64_t area);
  void atomicAddInstPlacedAreaUnscaled(int64_t area);

  int64_t getBinArea() const;
  int64_t getNonPlaceArea() const
  {
    return data_.get<BinColumns::nonPlaceArea>();
  }
  int64_t instPlacedArea() const
  {
    return data_.get<BinColumns::instPlacedArea>();
  }
  int64_t getNonPlaceAreaUnscaled() const
  {
    return data_.get<BinColumns::nonPlaceAreaUnscaled>();
  }
  int64_t getInstPlacedAreaUnscaled() const
  {
    return data_.get<BinColumns::instPlacedAreaUnscaled>();
  }

  int64_t getFillerArea() const { return data_.get<BinColumns::fillerArea>(); }

  Bin(SoaViewTag tag, void* store, size_t row) : data_(tag, store, row) {}
  SoaRef<Bin> ref() const { return {data_.storageIdentity(), data_.slot()}; }
  using Store = NumericStore;
  SoaRow<Store>& numericRow() { return data_; }
  const SoaRow<Store>& numericRow() const { return data_; }

 private:
  SoaRow<NumericStore> data_;
  // index

  // coordinate

  // will enable bin-wise density screening
};

inline int Bin::cx() const
{
  return (data_.get<BinColumns::ux>() + data_.get<BinColumns::lx>()) / 2;
}

inline int Bin::cy() const
{
  return (data_.get<BinColumns::uy>() + data_.get<BinColumns::ly>()) / 2;
}

inline int Bin::dx() const
{
  return (data_.get<BinColumns::ux>() - data_.get<BinColumns::lx>());
}

inline int Bin::dy() const
{
  return (data_.get<BinColumns::uy>() - data_.get<BinColumns::ly>());
}

inline void Bin::setNonPlaceArea(int64_t area)
{
  data_.get<BinColumns::nonPlaceArea>() = area;
}

inline void Bin::setNonPlaceAreaUnscaled(int64_t area)
{
  data_.get<BinColumns::nonPlaceAreaUnscaled>() = area;
}

inline void Bin::setInstPlacedArea(int64_t area)
{
  data_.get<BinColumns::instPlacedArea>() = area;
}

inline void Bin::setInstPlacedAreaUnscaled(int64_t area)
{
  data_.get<BinColumns::instPlacedAreaUnscaled>() = area;
}

inline void Bin::setFillerArea(int64_t area)
{
  data_.get<BinColumns::fillerArea>() = area;
}

inline void Bin::addNonPlaceArea(int64_t area)
{
  data_.get<BinColumns::nonPlaceArea>() += area;
}

inline void Bin::addInstPlacedArea(int64_t area)
{
  data_.get<BinColumns::instPlacedArea>() += area;
}

inline void Bin::addNonPlaceAreaUnscaled(int64_t area)
{
  data_.get<BinColumns::nonPlaceAreaUnscaled>() += area;
}

inline void Bin::addInstPlacedAreaUnscaled(int64_t area)
{
  data_.get<BinColumns::instPlacedAreaUnscaled>() += area;
}

inline void Bin::addFillerArea(int64_t area)
{
  data_.get<BinColumns::fillerArea>() += area;
}

// For the threaded density scatter. The accumulators are integers and each
// addend is truncated before it is added, so the sum is over a fixed multiset
// of integers -- associative and commutative, hence independent of the order
// threads reach it. Relaxed suffices: we need atomicity, not ordering, and the
// result is consumed only after the parallel region's implicit barrier.
// atomic_ref rather than an atomic member so the field stays a plain int64_t:
// Bin remains trivially copyable for its vector, and the serial path keeps its
// non-atomic add.
inline void Bin::atomicAddInstPlacedAreaUnscaled(int64_t area)
{
  std::atomic_ref<int64_t> ref(data_.get<BinColumns::instPlacedAreaUnscaled>());
  ref.fetch_add(area, std::memory_order_relaxed);
}

inline void Bin::atomicAddFillerArea(int64_t area)
{
  std::atomic_ref<int64_t> ref(data_.get<BinColumns::fillerArea>());
  ref.fetch_add(area, std::memory_order_relaxed);
}

//
// The bin can be non-uniform because of
// "integer" coordinates
//
class BinGrid
{
 public:
  BinGrid() = default;
  BinGrid(int lx, int ly, int ux, int uy);

  void setPlacerBase(std::shared_ptr<PlacerBase> pb);
  void setLogger(utl::Logger* log);
  void setRegionPoints(int lx, int ly, int ux, int uy);
  void setBinCnt(int binCntX, int binCntY);
  void setBinTargetDensity(float density);
  void updateBinsGCellDensityArea(const std::vector<GCellHandle>& cells,
                                  const std::vector<GCellRowList>& row_lists,
                                  int parallel_threads = 1);
  void scatterDensityAreaInPlace(const std::vector<GCellRowList>& row_lists,
                                 int parallel_threads);

  void initBins();

  // lx, ly, ux, uy will hold region area
  int lx() const;
  int ly() const;
  int ux() const;
  int uy() const;
  int cx() const;
  int cy() const;
  int dx() const;
  int dy() const;

  int getBinCntX() const;
  int getBinCntY() const;
  double getBinSizeX() const;
  double getBinSizeY() const;

  int64_t getOverflowArea() const;
  int64_t getOverflowAreaUnscaled() const;

  // GPU device-resident density path: the per-iteration scatter runs on
  // device, so the host bins are not updated; the device-reduced overflow
  // sums are stored here for the updateNextIter consumers.
  void setOverflowAreas(int64_t overflow_area, int64_t overflow_area_unscaled)
  {
    sumOverflowArea_ = overflow_area;
    sumOverflowAreaUnscaled_ = overflow_area_unscaled;
  }

  // return bins_ index with given gcell
  std::pair<int, int> getDensityMinMaxIdxX(SoaRef<GCell> gcell) const;
  std::pair<int, int> getDensityMinMaxIdxY(SoaRef<GCell> gcell) const;

  std::pair<int, int> getMinMaxIdxX(SoaRef<Instance> inst) const;
  std::pair<int, int> getMinMaxIdxY(SoaRef<Instance> inst) const;

  SoaPool<Bin>& getBins();
  const SoaPool<Bin>& getBinsConst() const { return bins_; };
  Bin::NumericStore& numericBins() { return *bin_columns_; }
  const Bin::NumericStore& numericBins() const { return *bin_columns_; }

  DensityGridView densityView() const;

  void updateBinsNonPlaceArea();

 private:
  std::shared_ptr<Bin::NumericStore> bin_columns_
      = std::make_shared<Bin::NumericStore>();
  SoaPool<Bin> bins_{bin_columns_};
  // Reused integer grids avoid atomics and preserve per-contribution
  // truncation.
  std::vector<int64_t> scatter_scratch_;
  std::vector<float> overflow_, overflow_unscaled_;
  std::vector<int> x_edges_, y_edges_;
  std::shared_ptr<PlacerBase> pb_;
  utl::Logger* log_ = nullptr;
  int lx_ = 0;
  int ly_ = 0;
  int ux_ = 0;
  int uy_ = 0;
  int binCntX_ = 0;
  int binCntY_ = 0;
  double binSizeX_ = 0;
  double binSizeY_ = 0;
  float targetDensity_ = 0;
  int64_t sumOverflowArea_ = 0;
  int64_t sumOverflowAreaUnscaled_ = 0;
  bool isSetBinCnt_ = false;
};

inline SoaPool<Bin>& BinGrid::getBins()
{
  return bins_;
}

struct NesterovBaseVars
{
  NesterovBaseVars(const PlaceOptions& options);

  const bool isSetBinCnt;
  const bool useUniformTargetDensity;
  bool placeIosMode = false;
  bool isMaxPhiCoefChanged = false;  // not user config
  const float targetDensity;
  const int binCntX;
  const int binCntY;

  const float minPhiCoef;
  float maxPhiCoef;  // may be updated after initialization
  const int initialPlacePerturbationSeed;
  const float initialPlacePerturbationDist;

  static constexpr float minWireLengthForceBar = -300;
};

struct NesterovPlaceVars
{
  // design_hpwl resolves referenceHpwl when the option is left at 0.
  NesterovPlaceVars(const PlaceOptions& options, int64_t design_hpwl = 0);

  int maxNesterovIter;
  static constexpr int maxBackTrack = 10;
  const float initDensityPenalty;                  // INIT_LAMBDA
  const float initWireLengthCoef;                  // base_wcof
  float targetOverflow;                            // overflow
  static constexpr float minPreconditioner = 1.0;  // MIN_PRE
  float initialPrevCoordiUpdateCoef = 100;         // z_ref_alpha
  const float referenceHpwl;                       // refDeltaHpwl

  // Floor for a derived referenceHpwl, the value this was hard-coded to
  // before it was derived at all. Designs at or below this scale keep exactly
  // the controller they had.
  static constexpr float kReferenceHpwlFloor = 446000000.0f;

  // Fraction of the design's HPWL a derived referenceHpwl spans. A wirelength
  // change of this size is what pulls the penalty ramp from its maximum to its
  // minimum, so it sets the controller's headroom before it saturates.
  static constexpr float kReferenceHpwlFraction = 0.01f;
  const float routability_end_overflow;
  const float routability_snapshot_overflow;
  const float keepResizeBelowOverflow;

  static constexpr int maxRecursionWlCoef = 10;
  static constexpr int maxRecursionInitSLPCoef = 10;

  bool timingDrivenMode;
  bool timingDrivenRepairTiming;
  float timingDrivenRepairTnsEndPercent;
  int timingDrivenIterCounter = 0;
  const bool routability_driven_mode;
  const bool disableRevertIfDiverge;

  bool debug = false;
  int debug_pause_iterations = 10;
  int debug_update_iterations = 10;
  bool debug_draw_bins = true;
  odb::dbInst* debug_inst = nullptr;
  odb::dbBTerm* debug_pin = nullptr;
  int debug_start_iter = 0;
  int debug_rudy_start = 5000;
  int debug_rudy_stride = 1;
  bool debug_generate_images = false;
  std::string debug_images_path = "REPORTS_DIR";
};

// Stores all pins, nets, and actual instances (static and movable)
// Used for calculating WL gradient
class NesterovBaseCommon
{
 public:
  NesterovBaseCommon(NesterovBaseVars nbVars,
                     std::shared_ptr<PlacerBaseCommon> pb,
                     utl::Logger* log,
                     int num_threads,
                     const Clusters& clusters);
  // Defined out-of-line (in nesterovBase.cpp) so the device_state_
  // std::unique_ptr<DeviceState> can default-destruct without exposing the
  // DeviceState definition (and its Kokkos types) in this header.
  ~NesterovBaseCommon();

  void reportInstanceExtensionByPinDensity() const;
  const SoaRows<GCell>& getGCells() const { return nbc_gcells_; }
  const SoaRows<GNet>& getGNets() const { return gNets_; }
  const SoaRows<GPin>& getGPins() const { return gPins_; }

  //
  // placerBase To NesterovBase functions
  //
  SoaRef<GCell> pbToNb(SoaRef<Instance> inst) const;
  SoaRef<GPin> pbToNb(SoaRef<Pin> pin) const;
  SoaRef<GNet> pbToNb(SoaRef<Net> net) const;

  //
  // OpenDB To NesterovBase functions
  //
  SoaRef<GCell> dbToNb(odb::dbInst* inst) const;
  SoaRef<GPin> dbToNb(odb::dbITerm* pin) const;
  SoaRef<GPin> dbToNb(odb::dbBTerm* pin) const;
  SoaRef<GNet> dbToNb(odb::dbNet* net) const;

  // WL force update based on WeightedAverage model
  // wlCoeffX : WireLengthCoefficient for X.
  //            equal to 1 / gamma_x
  // wlCoeffY : WireLengthCoefficient for Y.
  //            equal to 1 / gamma_y
  //
  // Gamma is described in the ePlaceMS paper.
  //
  // Public entry point — dispatches through wl_grad_backend_ (CPU or GPU).
  // Defined in wirelengthGradient.cpp.
  void updateWireLengthForceWA(float wlCoeffX, float wlCoeffY);

  // Native CPU body of updateWireLengthForceWA (the original OMP loop).
  // Called by CpuWirelengthGradientBackend; public so the backend in a
  // separate TU can dispatch into it. Defined in nesterovBase.cpp.
  void updateWireLengthForceWA_native(float wlCoeffX, float wlCoeffY);

  // Bulk per-cell wirelength gradient (hot path — replaces the
  // per-cell loop in NesterovBase::updateGradients). `out` is indexed
  // parallel to `gCells` (typically nb_gcells_, a per-NesterovBase view
  // into nbc gCellStor_). Defined in wirelengthGradient.cpp.
  void getAllWireLengthGradientsWA(const std::vector<GCellHandle>& gCells,
                                   PointColumns& out);

  // GPU path: run the per-inst gradient gather on device only (no host
  // copy) so the NB-level device scatter can consume it. No-op on the CPU
  // backend. Defined in wirelengthGradient.cpp.
  void prepareDeviceWlGradients();

  FloatPoint getWireLengthGradientPinWA(SoaRef<GPin> gPin,
                                        float wlCoeffX,
                                        float wlCoeffY) const;

  FloatPoint getWireLengthGradientWA(SoaRef<GCell> gCell,
                                     float wlCoeffX,
                                     float wlCoeffY) const;

  // for preconditioner
  FloatPoint getWireLengthPreconditioner(SoaRef<GCell> gCell) const;

  int64_t getHpwl();
  int64_t getHpwlColumns();

  // Refreshes host GNet boxes from the HPWL backend's last computeHpwl —
  // no-op on the CPU backend (always fresh). See hpwlBackend.h contract.
  void mirrorNetBoxesToHost();

  void updateDbGCells();

  // Device-resident state accessor (may be null when ENABLE_GPU is off).
  DeviceState* getDeviceState() { return device_state_.get(); }

  // Raw gCellStor_ accessor for DeviceState init (index correspondence).
  const SoaPool<GCell>& getGCellStor() const { return gCellStor_; }

  // Number of threads of execution
  size_t getNumThreads() { return num_threads_; }

  SoaRef<GCell> getGCellByIndex(size_t i);

  void setCbk(nesterovDbCbk* cbk) { db_cbk_ = cbk; }
  size_t createCbkGCell(odb::dbInst* db_inst);
  void createCbkGNet(odb::dbNet* net, bool skip_io_mode);
  void createCbkITerm(odb::dbITerm* iTerm);
  std::optional<std::pair<odb::dbInst*, size_t>> destroyCbkGCell(
      odb::dbInst* db_inst);
  void destroyCbkGNet(odb::dbNet*);
  void destroyCbkITerm(odb::dbITerm*);
  void resizeGCell(odb::dbInst* db_inst);
  void moveGCell(odb::dbInst* db_inst);
  void fixPointers();
  void updateDensityPinLocations(const GCellRowList& cells);
  void refreshCellLocks();
  void getWireLengthGradientsWA_cpu(const std::vector<GCellHandle>& cells,
                                    PointColumns& out,
                                    float coef_x,
                                    float coef_y);

  // Device-side sync for the timing-driven boundary; no-ops on CPU builds
  // and on the CPU runtime path. Non-virtual TD iterations create, destroy
  // (permuting storage indices), and resize instances — rebuildDeviceState()
  // reloads the DeviceState topology (sizes, CSRs, pin offsets, net weights)
  // from the repaired host storage. Virtual TD iterations only reweight
  // nets — refreshDeviceNetWeights() pushes the new weights.
  void rebuildDeviceState();
  void refreshDeviceNetWeights();

  void resetMinRcCellSize();
  void resizeMinRcCellSize();
  void updateMinRcCellSize();
  void revertGCellSizeToMinRc();

  SoaRef<GCell> getGCell(size_t index);
  size_t getGCellIndex(SoaRef<GCell> gCell) const;

  void printGCells();
  void printGPins();

  // TODO do this for each region? Also, manage this properly if other callbacks
  // are implemented.
  int64_t getDeltaArea() { return delta_area_; }
  void resetDeltaArea() { delta_area_ = 0; }
  int getNewGcellsCount() { return new_gcells_count_; }
  int getDeletedGcellsCount() { return deleted_gcells_count_; }
  void resetNewGcellsCount()
  {
    new_gcells_count_ = 0;
    deleted_gcells_count_ = 0;
  }

  NesterovBaseVars& getNbVars() { return nbVars_; }

 private:
  PointColumns pin_gradients_;
  bool pin_gradients_valid_ = false;
  uint64_t pin_gradient_weight_revision_ = 0;
  float pin_gradient_coef_x_ = 0, pin_gradient_coef_y_ = 0;
  void rebuildSoaTopology();
  void updateWireLengthForceWA_soa(float coef_x, float coef_y);
  FloatPoint getWireLengthGradientWA_soa(SoaRef<GCell> cell,
                                         float coef_x,
                                         float coef_y) const;

  NesterovBaseVars nbVars_;
  std::shared_ptr<PlacerBaseCommon> pbc_;
  utl::Logger* log_ = nullptr;

  std::shared_ptr<GCell::NumericStore> cell_columns_
      = std::make_shared<GCell::NumericStore>();
  std::shared_ptr<GPin::NumericStore> pin_columns_
      = std::make_shared<GPin::NumericStore>();
  std::shared_ptr<GNet::NumericStore> net_columns_
      = std::make_shared<GNet::NumericStore>();
  SoaPool<GCell> gCellStor_{cell_columns_};
  SoaPool<GNet> gNetStor_{net_columns_};
  SoaPool<GPin> gPinStor_{pin_columns_};

  SoaRows<GCell> nbc_gcells_;
  // For usage in routability mode, parallel to nbc_gcells_
  std::vector<int> min_rc_width_, min_rc_height_;
  SoaRows<GNet> gNets_;
  SoaRows<GPin> gPins_;

  std::vector<uint32_t> gCellMap_;
  std::vector<uint32_t> gPinMap_;
  std::vector<uint32_t> gNetMap_;

  boost::unordered::unordered_flat_map<odb::dbInst*, size_t>
      db_inst_to_nbc_index_map_;
  boost::unordered::unordered_flat_map<odb::dbNet*, size_t>
      db_net_to_index_map_;
  boost::unordered::unordered_flat_map<odb::dbITerm*, size_t>
      db_iterm_to_index_map_;
  boost::unordered::unordered_flat_map<odb::dbBTerm*, size_t>
      db_bterm_to_index_map_;

  // These three deques should not be required if placerBase allows for dynamic
  // modifications on its vectors.
  SoaPool<Instance> pb_insts_stor_;
  SoaPool<Net> pb_nets_stor_;
  SoaPool<Pin> pb_pins_stor_;

  int num_threads_;
  // Device-resident state for GPU backends (pin coords + per-net/per-pin
  // buffers; HPWL, WL grad, density gather all read from this).
  // Constructed in the ctor body after gCellStor_ / gPinStor_ / gNetStor_
  // are populated; null when ENABLE_GPU is off or gpl::gpuEnabled() returns
  // false. Must outlive hpwl_backend_ (backend borrows it), so it is
  // declared first and (since C++ destroys members in reverse declaration
  // order) destroyed last.
  std::unique_ptr<DeviceState> device_state_;
  std::unique_ptr<HpwlBackend> hpwl_backend_;
  // WA wirelength gradient dispatcher. CPU backend wraps the
  // updateWireLengthForceWA_native + per-cell helpers below; GPU backend
  // runs the 5-kernel Kokkos pipeline against device_state_'s pool.
  std::unique_ptr<WirelengthGradientBackend> wl_grad_backend_;
  int64_t delta_area_;
  int new_gcells_count_;
  int deleted_gcells_count_;
  nesterovDbCbk* db_cbk_{nullptr};
};

// Stores instances belonging to a specific power domain
// along with fillers and virtual blockages
// Also stores the bin grid for the power domain
// Used to calculate density gradient
class NesterovBase
{
 public:
  NesterovBase(NesterovBaseVars nbVars,
               std::shared_ptr<PlacerBase> pb,
               std::shared_ptr<NesterovBaseCommon> nbc,
               utl::Logger* log);
  ~NesterovBase();

  SoaRef<GCell> getFillerGCell(size_t index);

  NesterovBaseCommon* getNbc() { return nbc_.get(); }

  const std::vector<GCellHandle>& getGCells() const { return nb_gcells_; }

  // Live instances created since the last optimizer-state update. Their
  // callbacks have already built GCells; some may since have been resized.
  const std::unordered_set<odb::dbInst*>& getNewInstances() const
  {
    return new_instances_;
  }

  float getSumOverflow() const { return sum_overflow_; }
  float getSumOverflowUnscaled() const { return sum_overflow_unscaled_; }
  float getBaseWireLengthCoef() const { return baseWireLengthCoef_; }
  float getDensityPenalty() const { return densityPenalty_; }
  // Sets densityPenalty_ from the wirelength/density gradient ratio times
  // factor
  void updateDensityPenaltyFromRatio(float factor);

  float getWireLengthGradSum() const { return wireLengthGradSum_; }
  float getDensityGradSum() const { return densityGradSum_; }

  // update gCells with cx, cy
  void updateGCellCenterLocation(const PointColumns& coordis);

  void updateGCellDensityCenterLocation(const PointColumns& coordis);

  int getBinCntX() const;
  int getBinCntY() const;
  double getBinSizeX() const;
  double getBinSizeY() const;
  int64_t getOverflowArea() const;
  int64_t getOverflowAreaUnscaled() const;

  SoaPool<Bin>& getBins();
  const SoaPool<Bin>& getBinsConst() const { return bg_.getBinsConst(); };

  // filler cells / area control
  // will be used in Routability-driven loop
  int getFillerDx() const;
  int getFillerDy() const;
  int getFillerCnt() const;
  int64_t getFillerCellArea() const;
  int64_t getWhiteSpaceArea() const;
  int64_t getMovableArea() const;
  int64_t getTotalFillerArea() const;

  void setMovableArea(int64_t area) { movableArea_ = area; }

  // update
  // fillerArea, whiteSpaceArea, movableArea
  // and totalFillerArea after changing gCell's size
  void updateAreas();

  // update density sizes with changed dx and dy
  void updateDensitySize();

  // Account for retained timing repair and congestion inflation in this region.
  void updateTimingDrivenDensity();

  // should be separately defined.
  // This is mainly used for NesterovLoop
  int64_t getNesterovInstsArea() const;
  int64_t getInflationArea() const;
  int64_t getStdInstArea() const { return this->stdInstsArea_; }
  int64_t getMacroInstArea() const { return this->macroInstsArea_; }

  // sum phi and target density
  // used in NesterovPlace
  float getSumPhi() const;

  //
  // return uniform (lower bound) target density
  // LB of target density is required for massive runs.
  //
  float getUniformTargetDensity() const;

  // estimate target density based on the current placement
  float estimateTargetDensity(float overflow);

  // initTargetDensity is set by users
  // targetDensity is equal to initTargetDensity and
  // would be changed dynamically in RD loop
  //
  float initTargetDensity() const;
  float getTargetDensity() const;

  void setTargetDensity(float targetDensity);
  void checkConsistency();

  // RD can shrink the number of fillerCells.
  void cutFillerCells(int64_t targetFillerArea);

  void updateDensityCoordiLayoutInside(SoaRef<GCell> gcell);

  float getDensityCoordiLayoutInsideX(SoaRef<GCell> gCell, float cx) const;
  float getDensityCoordiLayoutInsideY(SoaRef<GCell> gCell, float cy) const;

  // FloatPoint getRegionGradient(SoaRef<GCell> gCell, FloatPoint nextLocation)
  // const;

  // for preconditioner
  FloatPoint getDensityPreconditioner(SoaRef<GCell> gCell) const;

  FloatPoint getDensityGradient(SoaRef<GCell> gCell) const;
  void getDensityGradientsColumns(PointColumns& out);

  // Fill out[i] with the density gradient for every filler cell in gCells, in
  // parallel. Used by the GPU density backend, whose per-inst gradients come
  // from the device but whose fillers (not in DeviceState) need the host
  // bin-overlap computation — a serial hotspot on the GPU path. Instance
  // entries (NesterovBaseCommon-backed) are left untouched.
  void fillFillerDensityGradients(const std::vector<GCellHandle>& gCells,
                                  PointColumns& out) const;

  // update electrostatic field within Bin
  void updateDensityFieldBin();

  BinGrid& getBinGrid() { return bg_; }

  // Nesterov Loop
  void initDensity1();
  float initDensity2(float wlCoeffX, float wlCoeffY);
  void setNpVars(const NesterovPlaceVars* npVars) { npVars_ = npVars; }
  void setIter(int iter) { iter_ = iter; }
  void setMaxPhiCoefChanged(bool maxPhiCoefChanged)
  {
    nbVars_.isMaxPhiCoefChanged = maxPhiCoefChanged;
  }

  void updateGradients(PointColumns& sumGrads,
                       PointColumns& wireLengthGrads,
                       PointColumns& densityGrads,
                       float wlCoeffX,
                       float wlCoeffY);

  void nbUpdatePrevGradient(float wlCoeffX, float wlCoeffY);
  void nbUpdateCurGradient(float wlCoeffX, float wlCoeffY);
  void nbUpdateNextGradient(float wlCoeffX, float wlCoeffY);

  void updateInitialPrevSLPCoordi();

  float getStepLength(const PointColumns& prevSLPCoordi_,
                      const PointColumns& prevSLPSumGrads_,
                      const PointColumns& curSLPCoordi_,
                      const PointColumns& curSLPSumGrads_);

  void updateNextIter(int iter);
  void setTrueReprintIterHeader() { reprint_iter_header_ = true; }
  float getPhiCoef(float scaledDiffHpwl) const;
  float getStoredPhiCoef() const { return phiCoef_; }
  float getStoredStepLength() const { return stepLength_; }
  float getStoredCoordiDistance() const { return coordiDistance_; }

  // True once the cells have stopped moving as fast as they were: the
  // per-iteration displacement has come down from its peak.
  //
  // Routability measures where wire will be from where the cells are, so it
  // needs a placement close to where it is going. That is a statement about
  // motion, not about density, and overflow does not capture it: the same
  // overflow can sit at very different displacements depending on how fast
  // the penalty schedule happens to be ramping.
  bool isSettled() const;
  float getSettleRatio() const;
  static constexpr float getSettleFraction() { return kSettleFraction; }

  float getStoredGradDistance() const { return gradDistance_; }

  bool checkConvergence(int gpl_iter_count,
                        int routability_gpl_iter_count,
                        RouteBase* rb);
  void resetConverged() { isConverged_ = false; }

  bool checkDivergence();
  void saveSnapshot();
  bool revertToSnapshot();

  // GPU device-resident density path: the hot loop no longer syncs coords
  // (or host GCell density centers) every iteration. Cold-path consumers
  // (snapshot save/revert, updateDb, routability filler cut/restore) call
  // this first to lazily refresh host state from the device cur slot.
  // No-op on the CPU path or when host state is already fresh. (TD and
  // routability runs never hold a device context, so they need no call.)
  void pullCoordsFromDevice();
  // pullCoordsFromDevice() plus the cur/prev SLP sum-grads, which the GPU
  // path keeps only on device. No-op on the CPU path.
  void pullSlpFromDevice();

  void updateDensityCenterCur();
  void updateDensityCenterCurSLP();
  void updateDensityCenterPrevSLP();
  void updateDensityCenterNextSLP();

  void nesterovUpdateCoordinates(float coeff);
  bool nesterovUpdateStepLength();
  void nesterovAdjustPhi();

  void resetMinSumOverflow();

  bool isDiverged() const { return isDiverged_; }
  // Resets isDiverged_ when no snapshot exists to revert to instead (the only
  // other place that clears it is revertToSnapshot()).
  void clearDivergence() { isDiverged_ = false; }

  void createCbkGCell(odb::dbInst* db_inst, size_t stor_index);
  std::optional<std::pair<odb::dbInst*, size_t>> destroyCbkGCell(
      odb::dbInst* db_inst);
  bool updateHandle(odb::dbInst* db_inst, size_t handle);

  // Must be called after fixPointers() to initialize internal values of gcells,
  // including parallel vectors.
  void updateGCellState();

  // Previous and current position and stored gradient of the instance's GCell
  // (the current gradient is the one the next Nesterov step uses).
  struct SlpState
  {
    FloatPoint prev_pos;
    FloatPoint prev_grad;
    FloatPoint cur_pos;
    FloatPoint cur_grad;
  };
  // std::nullopt when the instance is not in this region. Pulls from the
  // device first on the GPU path.
  std::optional<SlpState> getSlpState(odb::dbInst* db_inst);

  void destroyFillerGCell(size_t index_remove);
  void restoreRemovedFillers();
  void clearRemovedFillers() { removed_fillers_.clear(); }

  // Directly redistributes the existing filler into the free space.
  void redistributeFillerCells();

  void appendGCellCSVNote(const std::string& filename,
                          int iteration,
                          const std::string& message) const;
  // Helper to be used at nesterovPlace.cpp, inside core nesterov loop
  // Example:
  // for(auto& nb : nbVec_) {
  //   nb->writeGCellVectorsToCSV("gcells_vector.csv", nesterov_iter, 320, 1,
  //   1);
  // }
  void writeGCellVectorsToCSV(const std::string& filename,
                              int iteration,
                              int start_iteration = 0,
                              int iteration_stride = 50,
                              int gcell_index_stride = 10) const;

  std::shared_ptr<PlacerBase> getPb() const { return pb_; }

  odb::dbGroup* getGroup() const { return pb_->getGroup(); }

  std::pair<int, int> calculatePlacementPerturbationOffset(
      int dbu_per_micron) const;

  void updateDbIoPins();
  void refreshIoPins();
  // Null unless this region places the IO pins (-place_ios).
  IoPinBase* getIoPins() const { return io_.get(); }
  void setIoPins(std::unique_ptr<IoPinBase> io_pins);

 private:
  friend class TimingRepairIntegrationTest;

  template <bool GatherDensity>
  void combineGradientsColumns(PointColumns& sumGrads,
                               const PointColumns& wireLengthGrads,
                               PointColumns& densityGrads);

  NesterovBaseVars nbVars_;
  std::shared_ptr<PlacerBase> pb_;
  std::shared_ptr<NesterovBaseCommon> nbc_;
  utl::Logger* log_ = nullptr;
  mutable boost::random::mt19937 generator_;

  // Build (or rebuild) the GPU Nesterov device context against the current
  // nb_gcells_ size and sync host coords/grads into it. Called from
  // initDensity1 for the initial construction and from cutFillerCells /
  // restoreRemovedFillers after they resize nb_gcells_. No-op on CPU builds
  // and on GPU builds without a DeviceState (CPU runtime fallback).
  void rebuildNbDeviceCtx();

  // Scatter the named nb_device_ctx_ vector slot into DeviceState's per-inst
  // coord views, refresh device pin locations, and mark the DeviceState
  // coord flag fresh. Called after every GPU coord update (initDensity1,
  // updateInitialPrevSLPCoordi, nesterovUpdateCoordinates, revertToSnapshot,
  // rebuildNbDeviceCtx). No-op on CPU builds and when nb_device_ctx_ is null.
  void commitCoordsToDeviceState(SlpSlot source);

  BinGrid bg_;
  // Per-region FFT field Views (GPU). Declared before fft_ /
  // density_grad_backend_ / nb_device_ctx_ so it outlives the backends that
  // borrow it. Null on CPU builds and when the GPU path is off. Used by both
  // the device-resident density path (via nb_device_ctx_) and the host-staged
  // FFT path (fft_ / density_grad_backend_), so it lives on NesterovBase
  // rather than inside the device context (which is disabled in TD mode).
  std::unique_ptr<RegionDensityField> region_density_field_;
  std::unique_ptr<FFT> fft_;
  std::unique_ptr<DensityGradientBackend> density_grad_backend_;
  std::unique_ptr<NesterovDeviceContext> nb_device_ctx_;

  // True while the host coord vectors / GCell density centers mirror the
  // device state. The device-resident hot loop clears it each iteration;
  // pullCoordsFromDevice() restores it on demand. Always true on CPU.
  // Only ENABLE_GPU code reads it, but the field stays unconditional so
  // the class layout is identical in both builds (the ctor references it
  // to keep CPU-only clang builds quiet about the unused private field).
  bool host_coords_fresh_ = true;

  // Escape hatch for the device-resident density pipeline: set
  // GPL_GPU_HOST_DENSITY=1 to fall back to the per-iteration host
  // scatter/FFT staging path. Read in rebuildNbDeviceCtx() (reached from
  // initDensity1 and the filler cut/restore paths); always false in
  // TD/routability modes regardless of the env var. Same ENABLE_GPU-only
  // readers (and the same ctor reference) as host_coords_fresh_.
  bool use_device_density_ = false;

  int fillerDx_ = 0;
  int fillerDy_ = 0;
  int64_t whiteSpaceArea_ = 0;
  int64_t movableArea_ = 0;
  int64_t totalFillerArea_ = 0;
  int64_t initial_filler_area_ = 0;

  int64_t stdInstsArea_ = 0;
  int64_t macroInstsArea_ = 0;

  std::shared_ptr<GCell::NumericStore> filler_columns_
      = std::make_shared<GCell::NumericStore>();
  SoaPool<GCell> fillerStor_{filler_columns_};
  std::vector<GCellHandle> nb_gcells_;
  std::vector<GCellRowList> kernel_rows_;
  bool kernel_rows_dirty_ = true;
  PointColumns distance_scratch_;
  void appendCellHandle(GCellHandle handle);
  void eraseCellHandle(size_t index);
  void rebindCellHandle(size_t index, GCellHandle handle);
  void checkColumnState();
  const std::vector<GCellRowList>& kernelRows();
  void updateNextCoordinatesColumns(float coeff);

  std::unordered_map<odb::dbInst*, size_t> db_inst_to_nb_index_;
  std::unordered_map<size_t, size_t> filler_stor_index_to_nb_index_;

  // used to update gcell states after fixPointers() is called. Holds only
  // instances that are still alive: destroyCbkGCell() drops its entry, so a
  // deleted dbInst whose memory odb reuses cannot appear twice.
  std::unordered_set<odb::dbInst*> new_instances_;

  struct RemovedFillerColumns
  {
    explicit RemovedFillerColumns(std::shared_ptr<GCell::NumericStore> store)
        : gcells(std::move(store))
    {
    }
    SoaPool<GCell> gcells;
    PointColumns curSLPCoordi;
    PointColumns curSLPWireLengthGrads;
    PointColumns curSLPDensityGrads;
    PointColumns curSLPSumGrads;
    PointColumns nextSLPCoordi;
    PointColumns nextSLPWireLengthGrads;
    PointColumns nextSLPDensityGrads;
    PointColumns nextSLPSumGrads;
    PointColumns prevSLPCoordi;
    PointColumns prevSLPWireLengthGrads;
    PointColumns prevSLPDensityGrads;
    PointColumns prevSLPSumGrads;
    PointColumns curCoordi;
    PointColumns nextCoordi;
    PointColumns initCoordi;
    PointColumns snapshotCoordi;
    PointColumns snapshotSLPCoordi;
    PointColumns snapshotSLPSumGrads;
    PointColumns snapshotPrevSLPSumGrads;
    size_t size() const { return gcells.size(); }
    bool empty() const { return gcells.empty(); }
    void reserve(size_t count)
    {
      gcells.reserve(count);
      curSLPCoordi.reserve(count);
      curSLPWireLengthGrads.reserve(count);
      curSLPDensityGrads.reserve(count);
      curSLPSumGrads.reserve(count);
      nextSLPCoordi.reserve(count);
      nextSLPWireLengthGrads.reserve(count);
      nextSLPDensityGrads.reserve(count);
      nextSLPSumGrads.reserve(count);
      prevSLPCoordi.reserve(count);
      prevSLPWireLengthGrads.reserve(count);
      prevSLPDensityGrads.reserve(count);
      prevSLPSumGrads.reserve(count);
      curCoordi.reserve(count);
      nextCoordi.reserve(count);
      initCoordi.reserve(count);
      snapshotCoordi.reserve(count);
      snapshotSLPCoordi.reserve(count);
      snapshotSLPSumGrads.reserve(count);
      snapshotPrevSLPSumGrads.reserve(count);
    }
    void clear()
    {
      gcells.clear();
      curSLPCoordi.resize(0);
      curSLPWireLengthGrads.resize(0);
      curSLPDensityGrads.resize(0);
      curSLPSumGrads.resize(0);
      nextSLPCoordi.resize(0);
      nextSLPWireLengthGrads.resize(0);
      nextSLPDensityGrads.resize(0);
      nextSLPSumGrads.resize(0);
      prevSLPCoordi.resize(0);
      prevSLPWireLengthGrads.resize(0);
      prevSLPDensityGrads.resize(0);
      prevSLPSumGrads.resize(0);
      curCoordi.resize(0);
      nextCoordi.resize(0);
      initCoordi.resize(0);
      snapshotCoordi.resize(0);
      snapshotSLPCoordi.resize(0);
      snapshotSLPSumGrads.resize(0);
      snapshotPrevSLPSumGrads.resize(0);
    }
  };
  RemovedFillerColumns removed_fillers_{filler_columns_};

  float sumPhi_ = 0;
  float phiCoef_ = 0;
  float targetDensity_ = 0;
  float uniformTargetDensity_ = 0;

  // StepLength parameters (also included in the np debugPrint)
  // alpha
  float stepLength_ = 0;
  float coordiDistance_ = 0;
  float gradDistance_ = 0;

  // Largest per-iteration displacement seen so far, the yardstick isSettled()
  // measures the current one against.
  float peak_coordi_distance_ = 0;

  // Fraction of its own peak the per-iteration displacement must fall to for
  // the placement to count as settled. A ratio rather than a length, so it
  // carries no bin size, no DBU constant, and no dependence on the penalty
  // schedule's rate.
  static constexpr float kSettleFraction = 0.6f;

  // Nesterov loop data for each region, using parallel vectors
  // SLP is Step Length Prediction.
  //
  // y_st, y_dst, y_wdst, w_pdst
  PointColumns curSLPCoordi_;
  PointColumns curSLPWireLengthGrads_;
  PointColumns curSLPDensityGrads_;
  PointColumns curSLPSumGrads_;

  // y0_st, y0_dst, y0_wdst, y0_pdst
  PointColumns nextSLPCoordi_;
  PointColumns nextSLPWireLengthGrads_;
  PointColumns nextSLPDensityGrads_;
  PointColumns nextSLPSumGrads_;

  // z_st, z_dst, z_wdst, z_pdst
  PointColumns prevSLPCoordi_;
  PointColumns prevSLPWireLengthGrads_;
  PointColumns prevSLPDensityGrads_;
  PointColumns prevSLPSumGrads_;

  // x_st and x0_st
  PointColumns curCoordi_;
  PointColumns nextCoordi_;

  // save initial coordinates -- needed for RD
  PointColumns initCoordi_;

  // Snapshot data for routability, parallel vectors
  PointColumns snapshotCoordi_;
  PointColumns snapshotSLPCoordi_;
  PointColumns snapshotSLPSumGrads_;
  PointColumns snapshotPrevSLPSumGrads_;
  float snapshotDensityPenalty_ = 0;
  float snapshotStepLength_ = 0;

  // For destroying elements in parallel vectors
  void swapAndPop(PointColumns& vec, size_t remove_index, size_t last_index);
  void swapAndPopParallelVectors(size_t remove_index, size_t last_index);
  void appendParallelVectors();
  // Point whichever storage owns nb_gcells_[nb_index] back at nb_index.
  void rebindHandleIndex(size_t nb_index);

  float wireLengthGradSum_ = 0;
  float densityGradSum_ = 0;

  // opt_phi_cof
  float densityPenalty_ = 0;

  // base_wcof
  float baseWireLengthCoef_ = 0;

  // phi is described in ePlace paper.
  float sum_overflow_ = 0;
  float sum_overflow_unscaled_ = 0;
  float prev_reported_overflow_unscaled_ = 0;

  // half-parameter-wire-length
  int64_t prev_hpwl_ = 0;
  int64_t prev_reported_hpwl_ = 0;

  bool isDiverged_ = false;

  const NesterovPlaceVars* npVars_ = nullptr;

  float minSumOverflow_ = 1e30;
  float hpwlWithMinSumOverflow_ = 1e30;
  int iter_ = 0;
  bool isConverged_ = false;
  bool reprint_iter_header_ = false;

  void initFillerGCells();

  std::unique_ptr<IoPinBase> io_;
};

inline SoaPool<Bin>& NesterovBase::getBins()
{
  return bg_.getBins();
}

class biNormalParameters
{
 public:
  float meanX;
  float meanY;
  float sigmaX;
  float sigmaY;
  float lx;
  float ly;
  float ux;
  float uy;
};

class GCellHandle
{
 public:
  GCellHandle(NesterovBaseCommon* nbc, size_t idx)
      : storage_(nbc), storage_index_(idx)
  {
  }

  GCellHandle(NesterovBase* nb, size_t idx) : storage_(nb), storage_index_(idx)
  {
  }

  // Non-const versions
  SoaRef<GCell> operator->() { return getGCell(); }
  GCell operator*() { return *getGCell(); }
  operator SoaRef<GCell>() { return getGCell(); }

  // Const versions
  SoaRef<GCell> operator->() const { return getGCell(); }
  GCell operator*() const { return *getGCell(); }
  operator SoaRef<GCell>() const { return getGCell(); }

  bool isNesterovBaseCommon() const
  {
    return std::holds_alternative<NesterovBaseCommon*>(storage_);
  }

  void updateHandle(NesterovBaseCommon* nbc, size_t new_index)
  {
    storage_ = nbc;
    storage_index_ = new_index;
  }

  void updateHandle(NesterovBase* nb, size_t new_index)
  {
    storage_ = nb;
    storage_index_ = new_index;
  }

  size_t getStorageIndex() const { return storage_index_; }

 private:
  using StorageVariant = std::variant<NesterovBaseCommon*, NesterovBase*>;

  SoaRef<GCell> getGCell() const
  {
    if (std::holds_alternative<NesterovBaseCommon*>(storage_)) {
      return std::get<NesterovBaseCommon*>(storage_)->getGCell(storage_index_);
    }
    return std::get<NesterovBase*>(storage_)->getFillerGCell(storage_index_);
  }

  StorageVariant storage_;
  size_t storage_index_;
};

inline bool isValidSigType(const odb::dbSigType& db_type)
{
  return (db_type == odb::dbSigType::SIGNAL
          || db_type == odb::dbSigType::CLOCK);
}

inline constexpr const char* format_label_int = "{:27} {:10}";
inline constexpr const char* format_label_float = "{:27} {:10.4f}";
inline constexpr const char* format_label_um2 = "{:27} {:10.3f} um^2";
inline constexpr const char* format_label_percent = "{:27} {:10.2f} %";
inline constexpr const char* format_label_um2_with_delta
    = "{:27} {:10.3f} um^2 ({:+.2f}%)";

}  // namespace gpl
