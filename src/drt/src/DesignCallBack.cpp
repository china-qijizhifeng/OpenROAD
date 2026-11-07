// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2022-2025, The OpenROAD Authors

#include "DesignCallBack.h"

#include "drt/TritonRoute.h"
#include "frDesign.h"
#include "io/io.h"
#include "odb/db.h"

namespace drt {

void DesignCallBack::inDbInstCreate(odb::dbInst* db_inst)
{
  if (recording_) {
    router_->recordInstCreated(db_inst);
  }
  auto design = router_->getDesign();
  if (design == nullptr || design->getTopBlock() == nullptr) {
    return;
  }
  io::Parser parser(router_->getDb(),
                    design,
                    router_->getLogger(),
                    router_->getRouterConfiguration());

  parser.setInst(db_inst);
}

void DesignCallBack::inDbPreMoveInst(odb::dbInst* db_inst)
{
  router_->recordInstTouched(db_inst);
  auto design = router_->getDesign();
  if (design == nullptr || design->getTopBlock() == nullptr) {
    return;
  }
  auto inst = design->getTopBlock()->getInst(db_inst->getName());
  if (inst == nullptr) {
    return;
  }
  if (design->getRegionQuery() != nullptr) {
    design->getRegionQuery()->removeBlockObj(inst);
  }
  router_->deleteInstancePAData(inst, false);
}

void DesignCallBack::inDbPostMoveInst(odb::dbInst* db_inst)
{
  auto design = router_->getDesign();
  if (design == nullptr || design->getTopBlock() == nullptr) {
    return;
  }
  auto inst = design->getTopBlock()->getInst(db_inst->getName());
  if (inst == nullptr) {
    return;
  }
  int x, y;
  db_inst->getLocation(x, y);
  inst->setOrigin({x, y});
  inst->setOrient(db_inst->getOrient());
  router_->addInstancePAData(inst);
  if (design->getRegionQuery() != nullptr) {
    design->getRegionQuery()->addBlockObj(inst);
  }
}

void DesignCallBack::inDbInstDestroy(odb::dbInst* db_inst)
{
  if (recording_) {
    router_->recordInstDestroyed(db_inst);
  }
  auto design = router_->getDesign();
  if (design == nullptr || design->getTopBlock() == nullptr) {
    return;
  }
  auto inst = design->getTopBlock()->getInst(db_inst->getName());
  if (inst == nullptr) {
    return;
  }
  if (design->getRegionQuery() != nullptr) {
    design->getRegionQuery()->removeBlockObj(inst);
  }
  design->getTopBlock()->removeInst(inst);
  router_->deleteInstancePAData(inst, true);
}

void DesignCallBack::inDbInstSwapMasterBefore(odb::dbInst* db_inst,
                                              odb::dbMaster* db_master)
{
  // The master change is compared with the snapshot, not recorded as a
  // destroyed and a created instance.
  router_->recordInstTouched(db_inst);
  recording_ = false;
  inDbInstDestroy(db_inst);
  recording_ = true;
}

void DesignCallBack::inDbInstSwapMasterAfter(odb::dbInst* db_inst)
{
  recording_ = false;
  inDbInstCreate(db_inst);
  recording_ = true;
  inDbPostMoveInst(db_inst);
  for (auto inst_term : db_inst->getITerms()) {
    if (inst_term->getNet()) {
      inDbITermPostConnect(inst_term);
    }
  }
}
void DesignCallBack::inDbNetCreate(odb::dbNet* db_net)
{
  router_->recordChangedNet(db_net->getName());
  auto design = router_->getDesign();
  if (design == nullptr || design->getTopBlock() == nullptr) {
    return;
  }
  io::Parser parser(router_->getDb(),
                    design,
                    router_->getLogger(),
                    router_->getRouterConfiguration());
  parser.addNet(db_net);
}

void DesignCallBack::inDbNetDestroy(odb::dbNet* db_net)
{
  auto design = router_->getDesign();
  if (design == nullptr || design->getTopBlock() == nullptr) {
    return;
  }
  auto net = design->getTopBlock()->findNet(db_net->getName());
  if (net == nullptr) {
    return;
  }
  router_->recordChangedNetRegions(net);
  design->getTopBlock()->removeNet(net);
}

void DesignCallBack::inDbITermPostDisconnect(odb::dbITerm* db_iterm,
                                             odb::dbNet* db_net)
{
  recordInstRegion(db_iterm->getInst());
  router_->recordChangedNet(db_net->getName());
  auto design = router_->getDesign();
  if (design == nullptr || design->getTopBlock() == nullptr) {
    return;
  }
  auto inst = design->getTopBlock()->findInst(db_iterm->getInst());
  auto net = design->getTopBlock()->findNet(db_net->getName());
  if (inst == nullptr || net == nullptr) {
    return;
  }
  auto inst_term = inst->getInstTerm(db_iterm->getMTerm()->getIndex());
  if (inst_term == nullptr) {
    return;
  }
  net->removeInstTerm(inst_term);
  inst_term->addToNet(nullptr);
}

void DesignCallBack::inDbITermPostConnect(odb::dbITerm* db_iterm)
{
  recordInstRegion(db_iterm->getInst());
  router_->recordChangedNet(db_iterm->getNet()->getName());
  auto design = router_->getDesign();
  if (design == nullptr || design->getTopBlock() == nullptr) {
    return;
  }
  auto inst = design->getTopBlock()->findInst(db_iterm->getInst());
  auto net = design->getTopBlock()->findNet(db_iterm->getNet()->getName());
  if (inst == nullptr || net == nullptr) {
    return;
  }
  auto inst_term = inst->getInstTerm(db_iterm->getMTerm()->getIndex());
  if (inst_term == nullptr) {
    return;
  }
  inst_term->addToNet(net);
  net->addInstTerm(inst_term);
  router_->addInstancePAData(inst);
}

void DesignCallBack::inDbInstPlacementStatusBefore(
    odb::dbInst* inst,
    const odb::dbPlacementStatus&)
{
  router_->recordInstTouched(inst);
}

void DesignCallBack::inDbPostInstRename(odb::dbInst*, const char*)
{
  invalidate();
}

void DesignCallBack::inDbNetPostMerge(odb::dbNet*, odb::dbNet*)
{
  invalidate();
}

void DesignCallBack::inDbWireCreate(odb::dbWire* wire)
{
  recordWireNet(wire);
}

void DesignCallBack::inDbWireDestroy(odb::dbWire* wire)
{
  recordWireNet(wire);
}

void DesignCallBack::inDbWirePostModify(odb::dbWire* wire)
{
  recordWireNet(wire);
}

void DesignCallBack::inDbWirePostAttach(odb::dbWire* wire)
{
  recordWireNet(wire);
}

void DesignCallBack::inDbWirePostDetach(odb::dbWire*, odb::dbNet* net)
{
  router_->recordChangedNet(net->getName());
}

void DesignCallBack::inDbWirePostAppend(odb::dbWire* src, odb::dbWire* dst)
{
  recordWireNet(src);
  recordWireNet(dst);
}

void DesignCallBack::inDbWirePostCopy(odb::dbWire* src, odb::dbWire* dst)
{
  recordWireNet(src);
  recordWireNet(dst);
}

void DesignCallBack::inDbBTermCreate(odb::dbBTerm*)
{
  invalidate();
}

void DesignCallBack::inDbBTermDestroy(odb::dbBTerm*)
{
  invalidate();
}

void DesignCallBack::inDbBTermPostConnect(odb::dbBTerm*)
{
  invalidate();
}

void DesignCallBack::inDbBTermPostDisConnect(odb::dbBTerm*, odb::dbNet*)
{
  invalidate();
}

void DesignCallBack::inDbBPinCreate(odb::dbBPin*)
{
  invalidate();
}

void DesignCallBack::inDbBPinAddBox(odb::dbBox*)
{
  invalidate();
}

void DesignCallBack::inDbBPinRemoveBox(odb::dbBox*)
{
  invalidate();
}

void DesignCallBack::inDbBPinDestroy(odb::dbBPin*)
{
  invalidate();
}

void DesignCallBack::inDbBPinPlacementStatusBefore(
    odb::dbBPin*,
    const odb::dbPlacementStatus&)
{
  invalidate();
}

void DesignCallBack::inDbBlockageCreate(odb::dbBlockage*)
{
  invalidate();
}

void DesignCallBack::inDbBlockageDestroy(odb::dbBlockage*)
{
  invalidate();
}

void DesignCallBack::inDbObstructionCreate(odb::dbObstruction*)
{
  invalidate();
}

void DesignCallBack::inDbObstructionDestroy(odb::dbObstruction*)
{
  invalidate();
}

void DesignCallBack::inDbRegionCreate(odb::dbRegion*)
{
  invalidate();
}

void DesignCallBack::inDbRegionAddBox(odb::dbRegion*, odb::dbBox*)
{
  invalidate();
}

void DesignCallBack::inDbRegionDestroy(odb::dbRegion*)
{
  invalidate();
}

void DesignCallBack::inDbRowCreate(odb::dbRow*)
{
  invalidate();
}

void DesignCallBack::inDbRowDestroy(odb::dbRow*)
{
  invalidate();
}

void DesignCallBack::inDbSWireCreate(odb::dbSWire*)
{
  invalidate();
}

void DesignCallBack::inDbSWireDestroy(odb::dbSWire*)
{
  invalidate();
}

void DesignCallBack::inDbSWireAddSBox(odb::dbSBox*)
{
  invalidate();
}

void DesignCallBack::inDbSWireRemoveSBox(odb::dbSBox*)
{
  invalidate();
}

void DesignCallBack::inDbSWirePreDestroySBoxes(odb::dbSWire*)
{
  invalidate();
}

void DesignCallBack::inDbFillCreate(odb::dbFill*)
{
  invalidate();
}

void DesignCallBack::recordInstRegion(odb::dbInst* inst)
{
  router_->recordChangedRegion(inst->getBBox()->getBox());
}

void DesignCallBack::recordWireNet(odb::dbWire* wire)
{
  if (odb::dbNet* net = wire->getNet()) {
    router_->recordChangedNet(net->getName());
  }
}

void DesignCallBack::invalidate()
{
  router_->invalidateChangeTracking();
}

}  // namespace drt
