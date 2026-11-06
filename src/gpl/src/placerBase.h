// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#pragma once

#include <climits>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "boost/unordered/unordered_flat_map.hpp"
#include "odb/geom.h"
#include "soaStorage.h"

namespace odb {
class dbDatabase;

class dbBTerm;
class dbGroup;
class dbITerm;
class dbInst;
class dbNet;
class dbObject;

class dbPlacementStatus;
class dbSigType;

class dbBox;

class Rect;
class Point;

}  // namespace odb

namespace utl {
class Logger;
}

namespace gpl {

class Pin;
class Net;
class GCell;
class PlacerBaseCommon;
struct PlaceOptions;

class Instance;
struct InstanceColumns
{
  enum Field : size_t
  {
    inst,
    pins,
    lx,
    ly,
    ux,
    uy,
    extId,
    is_macro,
    is_locked
  };
  using Store = SoaStorage<odb::dbInst*,
                           SoaRows<Pin>,
                           int,
                           int,
                           int,
                           int,
                           int,
                           uint8_t,
                           uint8_t>;
};
struct PinColumns
{
  enum Field : size_t
  {
    term,
    inst,
    net,
    cx,
    cy,
    offsetCx,
    offsetCy,
    iTermField,
    bTermField,
    minPinXField,
    minPinYField,
    maxPinXField,
    maxPinYField
  };
  using Store = SoaStorage<odb::dbObject*,
                           SoaRef<Instance>,
                           SoaRef<Net>,
                           int,
                           int,
                           int,
                           int,
                           uint8_t,
                           uint8_t,
                           uint8_t,
                           uint8_t,
                           uint8_t,
                           uint8_t>;
};
struct NetColumns
{
  enum Field : size_t
  {
    net,
    pins,
    lx,
    ly,
    ux,
    uy
  };
  using Store = SoaStorage<odb::dbNet*, SoaRows<Pin>, int, int, int, int>;
};

class Instance
{
 public:
  using Store = InstanceColumns::Store;
  explicit Instance(std::shared_ptr<Store> storage = std::make_shared<Store>());
  Instance(odb::dbInst* inst,
           PlacerBaseCommon* pbc,
           utl::Logger* logger,
           std::shared_ptr<Store> storage = std::make_shared<Store>());
  Instance(int lx,
           int ly,
           int ux,
           int uy,
           std::shared_ptr<Store> storage
           = std::make_shared<Store>());  // dummy instance
  ~Instance() = default;
  Instance(const Instance&) = default;
  Instance(Instance&&) noexcept = default;
  Instance& operator=(const Instance&) = default;
  Instance& operator=(Instance&&) noexcept = default;

  odb::dbInst* dbInst() const { return data_.get<InstanceColumns::inst>(); }

  // a cell that no need to be moved.
  bool isFixed() const;

  // a instance that need to be moved.
  bool isInstance() const;

  bool isPlaceInstance() const;

  bool isMacro() const;

  // A placeable instance may be fixed during part of incremental placement.
  // It remains in the set of placeable objects though so as to simplify
  // the unlocking process.
  bool isLocked() const;
  void lock();
  void unlock();

  // Dummy is virtual instance to fill in
  // unusable sites.  It will have inst_ as nullptr
  bool isDummy() const;

  void copyDbLocation(PlacerBaseCommon* pbc);
  void setLocation(int x, int y);
  void setCenterLocation(int x, int y);

  void dbSetPlaced();
  void dbSetPlacementStatus(const odb::dbPlacementStatus& ps);
  void dbSetLocation();
  void dbSetLocation(int x, int y);
  void dbSetCenterLocation(int x, int y);

  int lx() const;
  int ly() const;
  int ux() const;
  int uy() const;
  int cx() const;
  int cy() const;
  int dx() const;
  int dy() const;
  int64_t getArea() const;

  void setExtId(int extId);
  int getExtId() const { return data_.get<InstanceColumns::extId>(); }

  void addPin(SoaRef<Pin> pin);
  const SoaRows<Pin>& getPins() const
  {
    return data_.get<InstanceColumns::pins>();
  }
  void snapOutward(const odb::Point& origin, int step_x, int step_y);
  int64_t extendSizeByScale(double scale, utl::Logger* logger);

  Instance(SoaViewTag tag, void* store, size_t row) : data_(tag, store, row) {}
  SoaRef<Instance> ref() const
  {
    return {data_.storageIdentity(), data_.slot()};
  }
  SoaRow<Store>& numericRow() { return data_; }
  const SoaRow<Store>& numericRow() const { return data_; }

 private:
  SoaRow<Store> data_;
};

class Pin
{
 public:
  using Store = PinColumns::Store;
  explicit Pin(std::shared_ptr<Store> storage = std::make_shared<Store>());
  Pin(odb::dbITerm* iTerm,
      std::shared_ptr<Store> storage = std::make_shared<Store>());
  Pin(odb::dbBTerm* bTerm,
      utl::Logger* logger,
      bool placeIosMode = false,
      std::shared_ptr<Store> storage = std::make_shared<Store>());
  ~Pin() = default;
  Pin(const Pin&) = default;
  Pin(Pin&&) noexcept = default;
  Pin& operator=(const Pin&) = default;
  Pin& operator=(Pin&&) noexcept = default;

  odb::dbITerm* getDbITerm() const;
  odb::dbBTerm* getDbBTerm() const;

  bool isITerm() const;
  bool isBTerm() const;
  bool isMinPinX() const;
  bool isMaxPinX() const;
  bool isMinPinY() const;
  bool isMaxPinY() const;

  void setITerm();
  void setBTerm();
  void setMinPinX();
  void setMinPinY();
  void setMaxPinX();
  void setMaxPinY();
  void unsetMinPinX();
  void unsetMinPinY();
  void unsetMaxPinX();
  void unsetMaxPinY();

  int cx() const;
  int cy() const;

  int getOffsetCx() const;
  int getOffsetCy() const;

  void updateLocation(SoaRef<Instance> inst);

  void setInstance(SoaRef<Instance> inst);
  void setNet(SoaRef<Net> net);

  bool isPlaceInstConnected() const;

  SoaRef<Instance> getInstance() const { return data_.get<PinColumns::inst>(); }
  SoaRef<Net> getNet() const { return data_.get<PinColumns::net>(); }
  std::string getName() const;

  void updateCoordi(odb::dbITerm* iTerm);

  Pin(SoaViewTag tag, void* store, size_t row) : data_(tag, store, row) {}
  SoaRef<Pin> ref() const { return {data_.storageIdentity(), data_.slot()}; }
  SoaRow<Store>& numericRow() { return data_; }
  const SoaRow<Store>& numericRow() const { return data_; }

 private:
  SoaRow<Store> data_;

  // pin center coordinate is enough
  // Pins' placed location.

  // offset coordinates inside instance.
  // origin point is center point of instance.
  // (e.g. (DX/2,DY/2) )
  // This will increase efficiency for bloating

  void updateCoordi(odb::dbBTerm* bTerm,
                    utl::Logger* logger,
                    bool placeIosMode = false);
};

class Net
{
 public:
  using Store = NetColumns::Store;
  explicit Net(std::shared_ptr<Store> storage = std::make_shared<Store>());
  Net(odb::dbNet* net,
      bool skipIoMode,
      std::shared_ptr<Store> storage = std::make_shared<Store>());
  ~Net() = default;
  Net(const Net&) = default;
  Net(Net&&) noexcept = default;
  Net& operator=(const Net&) = default;
  Net& operator=(Net&&) noexcept = default;

  int lx() const;
  int ly() const;
  int ux() const;
  int uy() const;
  int cx() const;
  int cy() const;

  // HPWL: half-parameter-wire-length
  int64_t getHpwl() const;

  void updateBox(bool skipIoMode = false);

  const SoaRows<Pin>& getPins() const { return data_.get<NetColumns::pins>(); }

  odb::dbNet* getDbNet() const { return data_.get<NetColumns::net>(); }
  odb::dbSigType getSigType() const;

  void addPin(SoaRef<Pin> pin);

  Net(SoaViewTag tag, void* store, size_t row) : data_(tag, store, row) {}
  SoaRef<Net> ref() const { return {data_.storageIdentity(), data_.slot()}; }
  SoaRow<Store>& numericRow() { return data_; }
  const SoaRow<Store>& numericRow() const { return data_; }

 private:
  SoaRow<Store> data_;
};

class Die
{
 public:
  Die();
  Die(const odb::Rect& dieRect, const odb::Rect& coreRect);
  ~Die();

  void setDieBox(const odb::Rect& dieRect);
  void setCoreBox(const odb::Rect& coreRect);

  int dieLx() const { return dieLx_; }
  int dieLy() const { return dieLy_; }
  int dieUx() const { return dieUx_; }
  int dieUy() const { return dieUy_; }

  int coreLx() const { return coreLx_; }
  int coreLy() const { return coreLy_; }
  int coreUx() const { return coreUx_; }
  int coreUy() const { return coreUy_; }

  int dieCx() const;
  int dieCy() const;
  int dieDx() const;
  int dieDy() const;
  int coreCx() const;
  int coreCy() const;
  int coreDx() const;
  int coreDy() const;

  int64_t dieArea() const;
  int64_t coreArea() const;

 private:
  int dieLx_ = 0;
  int dieLy_ = 0;
  int dieUx_ = 0;
  int dieUy_ = 0;
  int coreLx_ = 0;
  int coreLy_ = 0;
  int coreUx_ = 0;
  int coreUy_ = 0;
};

struct PlacerBaseVars
{
  PlacerBaseVars(const PlaceOptions& options);

  const int padLeft;
  const int padRight;
  const bool skipIoMode;
  const bool placeIosMode;
  const bool disablePinDensityAdjust;
};

// Class includes everything from PlacerBase that is not region specific
class PlacerBaseCommon
{
 public:
  auto instanceColumns() const { return instance_columns_; }
  auto pinColumns() const { return pin_columns_; }
  auto netColumns() const { return net_columns_; }

  // temp padLeft/Right before OpenDB supporting...
  PlacerBaseCommon(odb::dbDatabase* db,
                   PlacerBaseVars pbVars,
                   utl::Logger* log);
  ~PlacerBaseCommon();

  const SoaRows<Instance>& placeInsts() const { return placeInsts_; }
  const SoaRows<Instance>& getInsts() const { return insts_; }
  const SoaRows<Pin>& getPins() const { return pins_; }
  const SoaRows<Net>& getNets() const { return nets_; }

  Die& getDie() { return die_; }

  // Pb : PlacerBase
  SoaRef<Instance> dbToPb(odb::dbInst* inst) const;
  SoaRef<Pin> dbToPb(odb::dbITerm* term) const;
  SoaRef<Pin> dbToPb(odb::dbBTerm* term) const;
  SoaRef<Net> dbToPb(odb::dbNet* net) const;

  int siteSizeX() const { return siteSizeX_; }
  int siteSizeY() const { return siteSizeY_; }

  int getPadLeft() const { return pbVars_.padLeft; }
  int getPadRight() const { return pbVars_.padRight; }
  bool isSkipIoMode() const { return pbVars_.skipIoMode; }

  int64_t getHpwl() const;
  void printInfo() const;

  int64_t getMacroInstsArea() const { return macroInstsArea_; }

  odb::dbDatabase* db() const { return db_; }

  void unlockAll();

 private:
  odb::dbDatabase* db_ = nullptr;
  utl::Logger* log_ = nullptr;

  PlacerBaseVars pbVars_;

  Die die_;

  std::shared_ptr<Instance::Store> instance_columns_
      = std::make_shared<Instance::Store>();
  SoaPool<Instance> instStor_{instance_columns_};
  std::shared_ptr<Pin::Store> pin_columns_ = std::make_shared<Pin::Store>();
  std::shared_ptr<Net::Store> net_columns_ = std::make_shared<Net::Store>();
  SoaPool<Pin> pinStor_{pin_columns_};
  SoaPool<Net> netStor_{net_columns_};

  SoaRows<Instance> insts_;
  SoaRows<Pin> pins_;
  SoaRows<Net> nets_;

  SoaRows<Instance> placeInsts_;

  boost::unordered::unordered_flat_map<odb::dbInst*, uint32_t> instMap_;
  // The key is a dbITerm or a dbBTerm
  boost::unordered::unordered_flat_map<odb::dbObject*, uint32_t> pinMap_;
  boost::unordered::unordered_flat_map<odb::dbNet*, uint32_t> netMap_;

  int siteSizeX_ = 0;
  int siteSizeY_ = 0;

  int64_t macroInstsArea_ = 0;

  void init();
  void reset();
};

class PlacerBase
{
 public:
  PlacerBase();
  // temp padLeft/Right before OpenDB supporting...
  PlacerBase(odb::dbDatabase* db,
             std::shared_ptr<PlacerBaseCommon> pbCommon,
             utl::Logger* log,
             bool check_density,
             odb::dbGroup* group = nullptr);
  ~PlacerBase();

  const SoaRows<Instance>& getInsts() const { return pb_insts_; }

  //
  // placeInsts : a real instance that need to be placed
  // fixedInsts : a real instance that is fixed (e.g. macros, tapcells)
  // dummyInsts : a fake instance that is for unusable site handling
  //
  // nonPlaceInsts : fixedInsts + dummyInsts to enable fast-iterate on Bin-init
  //
  const SoaRows<Instance>& placeInsts() const { return placeInsts_; }
  const SoaRows<Instance>& fixedInsts() const { return fixedInsts_; }
  const SoaRows<Instance>& dummyInsts() const { return dummyInsts_; }
  const SoaRows<Instance>& nonPlaceInsts() const { return nonPlaceInsts_; }

  Die& getDie() { return die_; }
  int64_t getRegionArea() const { return region_area_; }
  const odb::Rect& getRegionBBox() const { return region_bbox_; }

  int getSiteSizeX() const { return siteSizeX_; }
  int getSiteSizeY() const { return siteSizeY_; }

  int64_t getHpwl() const;
  void printInfo(bool check_density) const;

  int64_t placeInstsArea() const { return placeInstsArea_; }
  int64_t nonPlaceInstsArea() const { return nonPlaceInstsArea_; }
  int64_t macroInstsArea() const;
  int64_t stdInstsArea() const { return stdInstsArea_; }

  odb::dbDatabase* db() const { return db_; }
  odb::dbGroup* getGroup() const { return group_; }

  void unlockAll();

 private:
  odb::dbDatabase* db_ = nullptr;
  utl::Logger* log_ = nullptr;

  Die die_;
  int64_t region_area_;
  odb::Rect region_bbox_;

  std::shared_ptr<Instance::Store> instance_columns_
      = std::make_shared<Instance::Store>();
  SoaPool<Instance> instStor_;

  SoaRows<Instance> pb_insts_;

  SoaRows<Instance> placeInsts_;
  SoaRows<Instance> fixedInsts_;
  SoaRows<Instance> dummyInsts_;
  SoaRows<Instance> nonPlaceInsts_;

  int siteSizeX_ = 0;
  int siteSizeY_ = 9;

  int64_t placeInstsArea_ = 0;
  int64_t nonPlaceInstsArea_ = 0;

  // macroInstsArea_ + stdInstsArea_ = placeInstsArea_;
  // macroInstsArea_ should be separated
  // because of target_density tuning
  int64_t macroInstsArea_ = 0;
  int64_t stdInstsArea_ = 0;

  std::shared_ptr<PlacerBaseCommon> pbCommon_;
  odb::dbGroup* group_ = nullptr;

  void init(bool check_density);
  void initInstsForUnusableSites();

  void reset();
};

}  // namespace gpl
