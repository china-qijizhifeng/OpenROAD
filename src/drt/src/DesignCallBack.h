// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2022-2025, The OpenROAD Authors

#include "odb/db.h"
#include "odb/dbBlockCallBackObj.h"
namespace drt {
class TritonRoute;
}
namespace drt {
class DesignCallBack : public odb::dbBlockCallBackObj
{
 public:
  DesignCallBack(TritonRoute* router) : router_(router) {}
  void inDbInstCreate(odb::dbInst*) override;
  void inDbPreMoveInst(odb::dbInst* inst) override;
  void inDbPostMoveInst(odb::dbInst* inst) override;
  void inDbInstDestroy(odb::dbInst* inst) override;
  void inDbInstSwapMasterBefore(odb::dbInst*, odb::dbMaster*) override;
  void inDbInstSwapMasterAfter(odb::dbInst*) override;
  void inDbNetCreate(odb::dbNet* net) override;
  void inDbNetDestroy(odb::dbNet* net) override;
  void inDbITermPostDisconnect(odb::dbITerm* iterm, odb::dbNet* net) override;
  void inDbITermPostConnect(odb::dbITerm* iterm) override;
  void inDbInstPlacementStatusBefore(odb::dbInst* inst,
                                     const odb::dbPlacementStatus&) override;
  void inDbPostInstRename(odb::dbInst* inst, const char* old_name) override;
  void inDbNetPostMerge(odb::dbNet* preserved_net,
                        odb::dbNet* removed_net) override;

  // Routing changes made outside the router: the affected nets are re-checked
  // by the next incremental run.
  void inDbWireCreate(odb::dbWire* wire) override;
  void inDbWireDestroy(odb::dbWire* wire) override;
  void inDbWirePostModify(odb::dbWire* wire) override;
  void inDbWirePostAttach(odb::dbWire* wire) override;
  void inDbWirePostDetach(odb::dbWire* wire, odb::dbNet* net) override;
  void inDbWirePostAppend(odb::dbWire* src, odb::dbWire* dst) override;
  void inDbWirePostCopy(odb::dbWire* src, odb::dbWire* dst) override;

  // Changes the router does not track by region.
  void inDbBTermCreate(odb::dbBTerm*) override;
  void inDbBTermDestroy(odb::dbBTerm*) override;
  void inDbBTermPostConnect(odb::dbBTerm*) override;
  void inDbBTermPostDisConnect(odb::dbBTerm*, odb::dbNet*) override;
  void inDbBPinCreate(odb::dbBPin*) override;
  void inDbBPinAddBox(odb::dbBox*) override;
  void inDbBPinRemoveBox(odb::dbBox*) override;
  void inDbBPinDestroy(odb::dbBPin*) override;
  void inDbBPinPlacementStatusBefore(odb::dbBPin*,
                                     const odb::dbPlacementStatus&) override;
  void inDbBlockageCreate(odb::dbBlockage*) override;
  void inDbBlockageDestroy(odb::dbBlockage*) override;
  void inDbObstructionCreate(odb::dbObstruction*) override;
  void inDbObstructionDestroy(odb::dbObstruction*) override;
  void inDbRegionCreate(odb::dbRegion*) override;
  void inDbRegionAddBox(odb::dbRegion*, odb::dbBox*) override;
  void inDbRegionDestroy(odb::dbRegion*) override;
  void inDbRowCreate(odb::dbRow*) override;
  void inDbRowDestroy(odb::dbRow*) override;
  void inDbSWireCreate(odb::dbSWire*) override;
  void inDbSWireDestroy(odb::dbSWire*) override;
  void inDbSWireAddSBox(odb::dbSBox*) override;
  void inDbSWireRemoveSBox(odb::dbSBox*) override;
  void inDbSWirePreDestroySBoxes(odb::dbSWire*) override;
  void inDbFillCreate(odb::dbFill*) override;

 private:
  void recordInstRegion(odb::dbInst* inst);
  void recordWireNet(odb::dbWire* wire);
  void invalidate();

  TritonRoute* router_;
  // False while a master swap reuses the destroy and create handlers.
  bool recording_{true};
};
}  // namespace drt
