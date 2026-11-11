// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Incremental BonnRoute (Ahrens, Henke, Rabenstein, Vygen,
// arXiv:2111.06169, Sec. 5.2): a detailed_route -bonnroute call on a design
// that already carries routing (the antenna-repair ECO after the first
// detailed routing) routes the nets that need routing against the existing
// routing, with rip-up and reroute, and writes them. FlexDR afterwards only
// checks the changed regions and repairs the markers there.
//
// The routing space (track graph, shape grid, fast grid, distance rules,
// pin access catalogues) stays alive in TritonRoute across detailed_route
// calls. Each call updates it with what changed since the previous call:
// instances (diodes inserted or moved), nets whose routing changed (wires
// destroyed by repair_antennas, routing changed by FlexDR) and nets whose
// pins changed. Existing routing is held as routed shapes of its net
// ("input routing"), which rip-up can remove like BonnRoute's own routes.

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "br/BonnRouteDR.h"
#include "br/BrIncrementalCore.h"
#include "db/obj/frBTerm.h"
#include "db/obj/frBlock.h"
#include "db/obj/frInst.h"
#include "db/obj/frInstBlockage.h"
#include "db/obj/frInstTerm.h"
#include "db/obj/frNet.h"
#include "db/obj/frShape.h"
#include "db/obj/frVia.h"
#include "db/tech/frLayer.h"
#include "db/tech/frViaDef.h"
#include "frDesign.h"
#include "frRegionQuery.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"
#include "utl/Logger.h"
#include "utl/timer.h"

namespace drt {

namespace {

uint64_t mix(uint64_t x)
{
  // splitmix64 finaliser.
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

uint64_t rectKey(const uint64_t kind, const int layer, const odb::Rect& r)
{
  uint64_t h = mix(kind * 1315423911ULL + static_cast<uint64_t>(layer));
  h = mix(h ^ static_cast<uint64_t>(static_cast<uint32_t>(r.xMin())));
  h = mix(h ^ static_cast<uint64_t>(static_cast<uint32_t>(r.yMin())));
  h = mix(h ^ static_cast<uint64_t>(static_cast<uint32_t>(r.xMax())));
  return mix(h ^ static_cast<uint64_t>(static_cast<uint32_t>(r.yMax())));
}

br::ShapeType routedShapeType(const br::ShapeRole role)
{
  return role == br::ShapeRole::kPreferredWire
                 || role == br::ShapeRole::kNonPreferredWire
             ? br::ShapeType::kWire
             : br::ShapeType::kVia;
}

const frInst* ownerInst(const frBlockObject* obj)
{
  if (obj->typeId() == frcInstTerm) {
    return static_cast<const frInstTerm*>(obj)->getInst();
  }
  if (obj->typeId() == frcInstBlockage) {
    return static_cast<const frInstBlockage*>(obj)->getInst();
  }
  return nullptr;
}

}  // namespace

bool BonnRouteDR::residentFor(const frDesign* design) const
{
  return design != nullptr && built_design_ == design
         && built_block_ == design->getTopBlock() && !layers_.empty();
}

uint64_t BonnRouteDR::routingFingerprint(const frNet* net)
{
  // Order-independent: FlexDR may reorder the routing objects.
  uint64_t h = 0;
  uint64_t n = 0;
  for (const auto& shape : net->getShapes()) {
    h += rectKey(1, shape->getLayerNum(), shape->getBBox());
    n++;
  }
  for (const auto& via : net->getVias()) {
    const odb::Point o = via->getOrigin();
    h += mix(mix(reinterpret_cast<uintptr_t>(via->getViaDef()))
             ^ rectKey(2, 0, odb::Rect(o, o)));
    n++;
  }
  for (const auto& patch : net->getPatchWires()) {
    h += rectKey(3, patch->getLayerNum(), patch->getBBox());
    n++;
  }
  return n == 0 ? 0 : mix(h ^ n);
}

uint64_t BonnRouteDR::termsFingerprint(const frNet* net)
{
  uint64_t h = 0;
  for (const frInstTerm* it : net->getInstTerms()) {
    h += mix(reinterpret_cast<uintptr_t>(it));
  }
  for (const frBTerm* bt : net->getBTerms()) {
    h += mix(reinterpret_cast<uintptr_t>(bt) ^ 0x5bd1e995ULL);
  }
  return mix(h ^ (net->getInstTerms().size() << 20)
             ^ net->getBTerms().size());
}

void BonnRouteDR::routingShapes(const frNet* net,
                                std::vector<RShape>& shapes) const
{
  shapes.clear();
  const int num_fr = static_cast<int>(z_of_fr_.size());
  auto wire_role = [&](const frLayerNum layer,
                       const odb::Rect& r,
                       const bool horizontal) {
    const int z = layer >= 0 && layer < num_fr ? z_of_fr_[layer] : -1;
    if (z < 0) {
      return br::ShapeRole::kFixed;
    }
    return horizontal == layers_[z].horizontal
               ? br::ShapeRole::kPreferredWire
               : br::ShapeRole::kNonPreferredWire;
  };
  for (const auto& shape : net->getShapes()) {
    const odb::Rect r = shape->getBBox();
    bool horizontal = r.dx() >= r.dy();
    if (shape->typeId() == frcPathSeg) {
      const auto [b, e] = static_cast<const frPathSeg*>(shape.get())->getPoints();
      if (b != e) {
        horizontal = b.y() == e.y();
      } else {
        const int z = shape->getLayerNum() >= 0
                              && shape->getLayerNum() < num_fr
                          ? z_of_fr_[shape->getLayerNum()]
                          : -1;
        horizontal = z >= 0 ? layers_[z].horizontal : horizontal;
      }
    }
    shapes.push_back({shape->getLayerNum(),
                      r,
                      wire_role(shape->getLayerNum(), r, horizontal)});
  }
  for (const auto& via : net->getVias()) {
    const frViaDef* def = via->getViaDef();
    shapes.push_back({def->getLayer1Num(),
                      via->getLayer1BBox(),
                      br::ShapeRole::kViaPad});
    shapes.push_back(
        {def->getCutLayerNum(), via->getCutBBox(), br::ShapeRole::kCut});
    shapes.push_back({def->getLayer2Num(),
                      via->getLayer2BBox(),
                      br::ShapeRole::kViaPad});
  }
  for (const auto& patch : net->getPatchWires()) {
    const odb::Rect r = patch->getBBox();
    shapes.push_back({patch->getLayerNum(),
                      r,
                      wire_role(patch->getLayerNum(), r, r.dx() >= r.dy())});
  }
}

void BonnRouteDR::setInputShapes(NetData& net,
                                 const bool present,
                                 br::FastGridBuilder::PendingRemovals* pending)
{
  if (net.input_present == present) {
    return;
  }
  if (present) {
    for (const RShape& s : net.input_shapes) {
      changeShape(s, routedShapeType(s.role), net.id, true, nullptr);
    }
    net.input_present = true;
    net.input_removed = false;
    return;
  }
  br::FastGridBuilder::PendingRemovals own;
  br::FastGridBuilder::PendingRemovals& batch
      = pending != nullptr ? *pending : own;
  for (const RShape& s : net.input_shapes) {
    changeShape(s, routedShapeType(s.role), net.id, false, nullptr, &batch);
  }
  if (pending == nullptr) {
    fast_builder_.flush(own);
  }
  net.input_present = false;
}

int BonnRouteDR::addNetData(frNet* net)
{
  NetData data;
  data.net = net;
  data.id = static_cast<int>(nets_.size());
  data.skip = true;
  odb::Rect box;
  box.mergeInit();
  for (const auto* it : net->getInstTerms()) {
    box.merge(it->getBBox());
  }
  for (const auto* bt : net->getBTerms()) {
    box.merge(bt->getBBox());
  }
  data.hpwl = box.isInverted() ? 0 : box.dx() + box.dy();
  net_index_[net] = data.id;
  nets_.push_back(std::move(data));
  return nets_.back().id;
}

bool BonnRouteDR::routableType(const NetData& data) const
{
  if (data.dead || data.net == nullptr) {
    return false;
  }
  const frNet* net = data.net;
  const auto type = net->getType();
  const size_t num_terms = net->getInstTerms().size() + net->getBTerms().size();
  return !(net->isSpecial() || net->isFake() || net->isFixed() || net->hasNDR()
           || type == odb::dbSigType::POWER || type == odb::dbSigType::GROUND
           || num_terms < 2 || !net->hasGuides());
}

void BonnRouteDR::recordFixedShape(const frBlockObject* obj,
                                   const RShape& shape,
                                   const br::ShapeType type,
                                   const int owner)
{
  const frInst* inst = ownerInst(obj);
  if (inst != nullptr) {
    inst_records_[inst].shapes.push_back({shape, type, owner});
  }
}

void BonnRouteDR::recordInstanceState(const frInst* inst,
                                      InstRecord& record) const
{
  record.name = inst->getName();
  record.bbox = inst->getBBox();
  record.orient = static_cast<int>(inst->getOrient().getValue());
  record.master = inst->getMaster();
  record.term_nets.clear();
  for (const auto& it : inst->getInstTerms()) {
    record.term_nets.push_back(it->getNet());
  }
}

void BonnRouteDR::readInstanceShapes(const frInst* inst, InstRecord& record)
{
  // The instance's pins and obstructions, as buildFixed reads them.
  frRegionQuery* rq = design_->getRegionQuery();
  odb::Rect window = inst->getBBox();
  window.bloat(layers_.front().pitch, window);
  const int num_fr = static_cast<int>(z_of_fr_.size());
  for (int fr = 0; fr < num_fr; fr++) {
    if (gridLayer(fr) < 0) {
      continue;
    }
    frRegionQuery::Objects<frBlockObject> result;
    rq->query(window, fr, result);
    for (const auto& [rect, obj] : result) {
      if (ownerInst(obj) != inst) {
        continue;
      }
      int owner = -1;
      if (obj->typeId() == frcInstTerm) {
        frNet* net = static_cast<frInstTerm*>(obj)->getNet();
        if (net != nullptr) {
          auto it = net_index_.find(net);
          owner = it != net_index_.end() ? it->second : addNetData(net);
        }
      }
      const RShape shape{static_cast<frLayerNum>(fr), rect};
      const br::ShapeType type
          = owner >= 0 ? br::ShapeType::kPin : br::ShapeType::kBlockage;
      changeShape(shape, type, owner, true, nullptr);
      if (owner >= 0) {
        nets_[owner].pins.push_back(shape);
      }
      record.shapes.push_back({shape, type, owner});
    }
  }
}

void BonnRouteDR::removeInstanceShapes(InstRecord& record, SyncStats& stats)
{
  br::FastGridBuilder::PendingRemovals pending;
  for (const FixedEntry& e : record.shapes) {
    changeShape(e.shape, e.type, e.owner, false, nullptr, &pending);
    stats.shapes_removed++;
    if (e.owner < 0 || e.owner >= static_cast<int>(nets_.size())) {
      continue;
    }
    auto& pins = nets_[e.owner].pins;
    auto it = std::ranges::find_if(pins, [&](const RShape& p) {
      return p.layer == e.shape.layer && p.rect == e.shape.rect;
    });
    if (it != pins.end()) {
      pins.erase(it);
    }
  }
  fast_builder_.flush(pending);
  record.shapes.clear();
}

void BonnRouteDR::syncInstances(SyncStats& stats)
{
  for (auto& [inst, record] : inst_records_) {
    record.seen = false;
  }
  InstRecord now;
  for (const auto& up : design_->getTopBlock()->getInsts()) {
    const frInst* inst = up.get();
    auto it = inst_records_.find(inst);
    if (it == inst_records_.end()) {
      InstRecord record;
      recordInstanceState(inst, record);
      const size_t before = record.shapes.size();
      readInstanceShapes(inst, record);
      stats.shapes_inserted += static_cast<int64_t>(record.shapes.size() - before);
      record.seen = true;
      stats.instances_changed++;
      stats.inst_regions.push_back(inst->getBBox());
      stats.changed_insts.insert(inst);
      inst_records_.emplace(inst, std::move(record));
      continue;
    }
    InstRecord& record = it->second;
    recordInstanceState(inst, now);
    if (record.name != now.name || record.bbox != now.bbox
        || record.orient != now.orient || record.master != now.master
        || record.term_nets != now.term_nets) {
      stats.inst_regions.push_back(record.bbox);
      stats.inst_regions.push_back(now.bbox);
      removeInstanceShapes(record, stats);
      record.name = std::move(now.name);
      record.bbox = now.bbox;
      record.orient = now.orient;
      record.master = now.master;
      record.term_nets = std::move(now.term_nets);
      readInstanceShapes(inst, record);
      stats.shapes_inserted += static_cast<int64_t>(record.shapes.size());
      stats.instances_changed++;
      stats.changed_insts.insert(inst);
    }
    record.seen = true;
  }
  for (auto it = inst_records_.begin(); it != inst_records_.end();) {
    if (it->second.seen) {
      ++it;
      continue;
    }
    stats.inst_regions.push_back(it->second.bbox);
    removeInstanceShapes(it->second, stats);
    stats.instances_changed++;
    it = inst_records_.erase(it);
  }
}

void BonnRouteDR::syncNets(SyncStats& stats)
{
  frBlock* block = design_->getTopBlock();
  // First pass: fingerprints; nets whose routing changed are re-read.
  std::vector<int> reread;
  std::vector<uint64_t> fp(nets_.size(), 0);
  std::vector<char> seen(nets_.size(), 0);
  int64_t changed_shapes = 0;
  for (const auto& up : block->getNets()) {
    frNet* net = up.get();
    auto it = net_index_.find(net);
    int id;
    if (it == net_index_.end()) {
      id = addNetData(net);
      stats.nets_added++;
      fp.resize(nets_.size(), 0);
      seen.resize(nets_.size(), 0);
    } else {
      id = it->second;
    }
    seen[id] = 1;
    NetData& data = nets_[id];
    fp[id] = routingFingerprint(net);
    if (fp[id] != data.routing_fp || data.routed || !data.shapes.empty()) {
      reread.push_back(id);
      changed_shapes += static_cast<int64_t>(data.input_shapes.size())
                        + static_cast<int64_t>(net->getShapes().size())
                        + 3 * static_cast<int64_t>(net->getVias().size());
    }
  }
  // Many changes (FlexDR's clean-up after the first BonnRoute run): update
  // the shape grid only and rebuild the fast grid once.
  const bool bulk = changed_shapes > 200000;
  bulk_insert_ = bulk;
  std::vector<RShape> shapes;
  for (const int id : reread) {
    NetData& data = nets_[id];
    br::FastGridBuilder::PendingRemovals pending;
    if (data.input_present) {
      stats.shapes_removed += static_cast<int64_t>(data.input_shapes.size());
      setInputShapes(data, false, &pending);
    }
    for (const RShape& s : data.shapes) {
      changeShape(s, routedShapeType(s.role), data.id, false, nullptr, &pending);
      stats.shapes_removed++;
    }
    fast_builder_.flush(pending);
    data.shapes.clear();
    data.wires.clear();
    data.vias.clear();
    data.routed = false;
    routingShapes(data.net, data.input_shapes);
    setInputShapes(data, true);
    stats.shapes_inserted += static_cast<int64_t>(data.input_shapes.size());
    data.routing_fp = fp[id];
    stats.nets_reread++;
  }
  bulk_insert_ = false;
  if (bulk) {
    rebuildFastGrid();
  }
  for (NetData& data : nets_) {
    if (data.dead || seen[data.id]) {
      continue;
    }
    // Deleted net: its routing leaves the routing space.
    setInputShapes(data, false);
    net_index_.erase(data.net);
    data.net = nullptr;
    data.dead = true;
    data.skip = true;
    stats.nets_removed++;
  }
}

void BonnRouteDR::insertDesignRouting()
{
  bulk_insert_ = true;
  for (NetData& data : nets_) {
    routingShapes(data.net, data.input_shapes);
    if (!data.input_shapes.empty()) {
      setInputShapes(data, true);
    }
    data.routing_fp = routingFingerprint(data.net);
    data.terms_fp = termsFingerprint(data.net);
  }
  bulk_insert_ = false;
  rebuildFastGrid();
  shape_grid_.compact();
  grid_.compact();
}

void BonnRouteDR::buildFromDesign()
{
  utl::Timer timer;
  if (modernSearch() && params_.use_intervals) {
    params_.use_intervals = false;
  }
  buildTrackGraph();
  if (layers_.empty()) {
    return;
  }
  buildRoutingSpace();
  if (params_.spreading) {
    buildSpreadTiles();
  }
  workers_.resize(std::max(1, router_cfg_->MAX_THREADS));
  for (Worker& w : workers_) {
    w.region.init(graph_.numLayers(), tiles_x_, tiles_y_);
    if (params_.use_intervals) {
      w.interval_search = std::make_unique<br::IntervalSearch>(graph_, grid_);
    } else {
      w.search = std::make_unique<br::PathSearch>(graph_, grid_);
    }
  }
  buildNets();
  record_fixed_ = true;
  buildFixed();
  buildSameNetRules();
  if (params_.spreading && params_.spread_from_extra_space) {
    loadExtraSpace();
  }
  for (const auto& up : design_->getTopBlock()->getInsts()) {
    recordInstanceState(up.get(), inst_records_[up.get()]);
  }
  insertDesignRouting();
  built_design_ = design_;
  built_block_ = design_->getTopBlock();
  if (params_.verbose) {
    logger_->report(
        "BonnRoute incremental: routing space built from the design in "
        "{:.2f}s ({} nets, {} instances).",
        timer.elapsed(),
        nets_.size(),
        inst_records_.size());
  }
}

void BonnRouteDR::keepRoutingSpace()
{
  for (const auto& up : design_->getTopBlock()->getInsts()) {
    recordInstanceState(up.get(), inst_records_[up.get()]);
  }
  for (NetData& data : nets_) {
    setReservations(data, false);
    if (data.routed) {
      // The written routes become input routing (their shapes stay).
      data.input_shapes = std::move(data.shapes);
      data.input_present = true;
    } else if (!data.shapes.empty()) {
      // Routed but dropped when written (contact guard): not in the design.
      br::FastGridBuilder::PendingRemovals pending;
      for (const RShape& s : data.shapes) {
        changeShape(s, routedShapeType(s.role), data.id, false, nullptr, &pending);
      }
      fast_builder_.flush(pending);
    }
    data.shapes.clear();
    data.wires.clear();
    data.vias.clear();
    data.pin_points.clear();
    data.routed = false;
    data.failed = false;
    data.routing_fp = routingFingerprint(data.net);
    data.terms_fp = termsFingerprint(data.net);
  }
  built_design_ = design_;
  built_block_ = design_->getTopBlock();
}

bool BonnRouteDR::prepareNet(NetData& net,
                             const std::unordered_set<const frInst*>& changed)
{
  if (net.prepared) {
    return !net.skip;
  }
  net.prepared = true;
  if (!routableType(net)) {
    net.skip = true;
    return false;
  }
  setReservations(net, false);
  // Catalogues of pins on unchanged instances are kept (Sec. 4.3); the
  // searches check their legality against the current routing space.
  // Only pins still on the net are looked at (a removed pin's instance
  // may be gone).
  std::unordered_set<const void*> current;
  for (const frInstTerm* it : net.net->getInstTerms()) {
    current.insert(it);
  }
  for (const frBTerm* bt : net.net->getBTerms()) {
    current.insert(bt);
  }
  std::unordered_map<const void*, Term> kept;
  for (Term& term : net.terms) {
    const void* key = term.inst_term != nullptr
                          ? static_cast<const void*>(term.inst_term)
                          : static_cast<const void*>(term.bterm);
    if (key == nullptr || term.options.empty() || !current.contains(key)
        || (term.inst_term != nullptr
            && changed.contains(term.inst_term->getInst()))) {
      continue;
    }
    kept.emplace(key, std::move(term));
  }
  net.terms.clear();
  net.input_reservation.clear();
  net.discount_runs.clear();
  net.skip = false;
  buildGuideTiles(net);
  buildAccess(net);
  size_t t = 0;
  for (frInstTerm* it : net.net->getInstTerms()) {
    net.terms[t++].inst_term = it;
  }
  for (frBTerm* bt : net.net->getBTerms()) {
    net.terms[t++].bterm = bt;
  }
  if (params_.offtrack_access) {
    AccessStats stats;
    for (Term& term : net.terms) {
      const void* key = term.inst_term != nullptr
                            ? static_cast<const void*>(term.inst_term)
                            : static_cast<const void*>(term.bterm);
      auto f = kept.find(key);
      if (f != kept.end()) {
        term = std::move(f->second);
        term.reservation_active = false;
        term.chosen = -1;
        if (term.primary >= 0 && term.primary < static_cast<int>(term.options.size())) {
          accessShapes(term.options[term.primary], term.reservation);
          setReservation(term, net.id, true);
        }
        continue;
      }
      std::vector<AccessOption> catalogue;
      offTrackCatalogue(net.id, term, nullptr, catalogue, stats);
      if (catalogue.empty()) {
        // The access stubs: the first conflict-free one is reserved.
        for (size_t o = 0; o < term.options.size(); o++) {
          std::vector<RShape> shapes;
          accessShapes(term.options[o], shapes);
          if (std::ranges::none_of(shapes, [&](const RShape& s) {
                return shapeConflicts(s, net.id);
              })) {
            term.primary = static_cast<int>(o);
            term.reservation = std::move(shapes);
            setReservation(term, net.id, true);
            break;
          }
        }
        continue;
      }
      std::vector<AccessOption> legal;
      std::vector<int> index;
      for (size_t o = 0; o < catalogue.size(); o++) {
        if (optionLegal(catalogue[o], net.id, nullptr)) {
          legal.push_back(catalogue[o]);
          index.push_back(static_cast<int>(o));
        }
      }
      const std::vector<int> choice = choosePrimaries({&legal}, {net.id}, stats);
      keepStubs(term, catalogue);
      if (!choice.empty() && choice[0] >= 0) {
        term.primary = index[choice[0]];
        accessShapes(term.options[term.primary], term.reservation);
        setReservation(term, net.id, true);
      }
    }
  } else {
    reserveAccess(net);
  }
  net.skip = std::ranges::any_of(
      net.terms, [](const Term& term) { return term.options.empty(); });
  return !net.skip;
}

bool BonnRouteDR::prepareVictim(const int id)
{
  if (!incremental_active_) {
    return true;
  }
  NetData& net = nets_[id];
  if (net.dead || victim_changed_insts_ == nullptr) {
    return false;
  }
  return prepareNet(net, *victim_changed_insts_);
}

void BonnRouteDR::routeIncremental(std::vector<odb::Rect>& regions)
{
  utl::Timer timer;
  regions.clear();
  incremental_active_ = true;
  SyncStats stats;
  syncInstances(stats);
  syncNets(stats);
  const double sync_s = timer.elapsed();
  victim_changed_insts_ = &stats.changed_insts;

  // Nets to route: without routing, or with routing whose pins changed (the
  // input routing is ripped up and the net routed again).
  std::vector<int> work;
  int64_t unrouted = 0;
  int64_t pins_changed = 0;
  for (NetData& data : nets_) {
    if (!routableType(data)) {
      continue;
    }
    const uint64_t tfp = termsFingerprint(data.net);
    const bool terms_changed = tfp != data.terms_fp;
    data.terms_fp = tfp;
    if (!data.input_present) {
      work.push_back(data.id);
      unrouted++;
    } else if (terms_changed) {
      setInputShapes(data, false);
      data.input_removed = true;
      work.push_back(data.id);
      pins_changed++;
    }
  }
  std::vector<int> order;
  std::vector<int> not_prepared;
  for (const int id : work) {
    if (prepareNet(nets_[id], stats.changed_insts)) {
      order.push_back(id);
    } else {
      not_prepared.push_back(id);
    }
  }
  const double access_s = timer.elapsed() - sync_s;
  std::ranges::stable_sort(
      order, [&](int a, int b) { return nets_[a].hpwl < nets_[b].hpwl; });
  const int64_t ripups_before = ripup_events_;
  std::vector<int> failed;
  routeParallel(order, failed);
  const int64_t first_failed = static_cast<int64_t>(failed.size());
  if (params_.ripup && !failed.empty()) {
    routeRipupQueue(failed);
  }
  const double route_s = timer.elapsed() - sync_s - access_s;

  // Write: nets routed now, and nets whose input routing was ripped up.
  std::vector<int> touched;
  for (NetData& data : nets_) {
    if (!data.dead && (data.routed || data.input_removed)) {
      touched.push_back(data.id);
    }
  }
  write_mask_.assign(nets_.size(), 0);
  for (const int id : touched) {
    NetData& data = nets_[id];
    if (data.input_removed) {
      data.net->clearRoutes();
      data.net->setHasInitialRouting(false);
      data.net->setModified(true);
    }
    write_mask_[id] = data.routed;
  }
  writeNets();
  write_mask_.clear();
  frBlock* block = design_->getTopBlock();
  auto add_guides = [&](const NetData& data) {
    for (const auto& guide : data.net->getGuides()) {
      odb::Rect box
          = block->getGCellBox(block->getGCellIdx(guide->getBeginPoint()));
      box.merge(block->getGCellBox(block->getGCellIdx(guide->getEndPoint())));
      regions.push_back(box);
    }
  };
  int64_t written = 0;
  int64_t left = 0;
  for (const int id : touched) {
    NetData& data = nets_[id];
    if (data.routed) {
      for (const RShape& s : data.shapes) {
        regions.push_back(s.rect);
      }
      data.input_shapes = std::move(data.shapes);
      data.input_present = true;
      written++;
    } else {
      // Not written (contact guard) or rip-up victim left unrouted:
      // FlexDR routes it from its guides.
      br::FastGridBuilder::PendingRemovals pending;
      for (const RShape& s : data.shapes) {
        changeShape(s, routedShapeType(s.role), data.id, false, nullptr, &pending);
      }
      fast_builder_.flush(pending);
      add_guides(data);
      left++;
    }
    data.shapes.clear();
    data.wires.clear();
    data.vias.clear();
    data.pin_points.clear();
    data.routed = false;
    data.input_removed = false;
    data.routing_fp = routingFingerprint(data.net);
  }
  for (const int id : order) {
    NetData& data = nets_[id];
    if (!data.input_present && !data.net->hasInitialRouting()) {
      // Failed in BonnRoute and never written.
      if (std::ranges::find(touched, id) == touched.end()) {
        add_guides(data);
        left++;
      }
    }
  }
  for (const int id : not_prepared) {
    add_guides(nets_[id]);
    left++;
  }
  regions.insert(
      regions.end(), stats.inst_regions.begin(), stats.inst_regions.end());
  for (NetData& data : nets_) {
    if (data.prepared) {
      setReservations(data, false);
      data.prepared = false;
      data.failed = false;
    }
  }
  victim_changed_insts_ = nullptr;
  incremental_active_ = false;
  const double write_s = timer.elapsed() - sync_s - access_s - route_s;
  if (params_.verbose) {
    logger_->report(
        "BonnRoute incremental: sync_s={:.3f} access_s={:.3f} route_s={:.3f} "
        "write_s={:.3f} total_s={:.3f} instances_changed={} nets_reread={} "
        "nets_added={} nets_removed={} shapes_removed={} shapes_inserted={} "
        "nets_to_route={} unrouted={} pins_changed={} not_prepared={} "
        "first_pass_failed={} ripups={} written={} left_to_flexdr={} "
        "regions={}",
        sync_s,
        access_s,
        route_s,
        write_s,
        timer.elapsed(),
        stats.instances_changed,
        stats.nets_reread,
        stats.nets_added,
        stats.nets_removed,
        stats.shapes_removed,
        stats.shapes_inserted,
        work.size(),
        unrouted,
        pins_changed,
        not_prepared.size(),
        first_failed,
        ripup_events_ - ripups_before,
        written,
        left,
        regions.size());
  }
}

}  // namespace drt
