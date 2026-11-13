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
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
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

void BonnRouteDR::gridOnlyChange(const RShape& shape,
                                 const br::ShapeType type,
                                 const int owner,
                                 const bool insert)
{
  const int g = gridLayer(shape.layer);
  if (g < 0) {
    return;
  }
  const br::Rect r{shape.rect.xMin(),
                   shape.rect.yMin(),
                   shape.rect.xMax(),
                   shape.rect.yMax()};
  auto apply = [&](const int layer, const br::ClassId cls) {
    const bool ok
        = insert
              ? shape_grid_.insert(layer, r, cls, type, owner, serial_cache_)
              : shape_grid_.remove(layer, r, cls, type, owner, serial_cache_);
    if (!ok) {
      inconsistent_updates_++;
    }
  };
  apply(g, shapeClass(shape, g));
  const int above = projectionLayer(g);
  if (above >= 0) {
    RShape projection = shape;
    projection.role = br::ShapeRole::kCutProjection;
    apply(above, shapeClass(projection, above));
  }
}

void BonnRouteDR::enableIncrementalDiscounts()
{
  if (incremental_discounts_ || !params_.incremental_reservations
      || params_.discount_denominator <= 1) {
    incremental_discounts_ = incremental_discounts_ || graph_.wire_unit > 1;
    return;
  }
  // The first call may have used unit wire costs (no reservations there):
  // scale every cost of the resident track graph and of the kept access
  // options to the discount denominator, so that discounted runs stay
  // integral (arXiv:2111.06169 Sec. 5.2). Path costs keep their order.
  if (graph_.wire_unit == 1) {
    const br::Cost f = params_.discount_denominator;
    graph_.wire_unit = f;
    for (br::TrackLayer& l : graph_.layers) {
      l.jog_factor *= f;
    }
    for (br::Cost& c : graph_.via_cost) {
      c *= f;
    }
    for (NetData& net : nets_) {
      for (Term& term : net.terms) {
        for (AccessOption& o : term.options) {
          o.cost *= f;
        }
      }
    }
  }
  // Discounted runs need the vertex search (the interval search assumes
  // unit wire costs).
  params_.use_intervals = false;
  for (Worker& w : workers_) {
    if (!w.search) {
      w.search = std::make_unique<br::PathSearch>(graph_, grid_);
    }
  }
  incremental_discounts_ = true;
}

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
  return mix(h ^ (net->getInstTerms().size() << 20) ^ net->getBTerms().size());
}

void BonnRouteDR::routingShapes(const frNet* net,
                                std::vector<RShape>& shapes) const
{
  shapes.clear();
  // Existing routing of an NDR net keeps the NDR spacing on its metal
  // (drt's checker applies it), except pieces tapered at pins.
  const int ndr = ndrType(net);
  const int num_fr = static_cast<int>(z_of_fr_.size());
  auto wire_role
      = [&](const frLayerNum layer, const odb::Rect& r, const bool horizontal) {
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
      const auto [b, e]
          = static_cast<const frPathSeg*>(shape.get())->getPoints();
      if (b != e) {
        horizontal = b.y() == e.y();
      } else {
        const int z = shape->getLayerNum() >= 0 && shape->getLayerNum() < num_fr
                          ? z_of_fr_[shape->getLayerNum()]
                          : -1;
        horizontal = z >= 0 ? layers_[z].horizontal : horizontal;
      }
    }
    const bool tapered
        = shape->typeId() == frcPathSeg
          && static_cast<const frPathSeg*>(shape.get())->isTapered();
    shapes.push_back({shape->getLayerNum(),
                      r,
                      wire_role(shape->getLayerNum(), r, horizontal),
                      tapered ? 0 : ndr});
  }
  for (const auto& via : net->getVias()) {
    const frViaDef* def = via->getViaDef();
    const int pad_ndr = via->isTapered() ? 0 : ndr;
    shapes.push_back({def->getLayer1Num(),
                      via->getLayer1BBox(),
                      br::ShapeRole::kViaPad,
                      pad_ndr});
    shapes.push_back(
        {def->getCutLayerNum(), via->getCutBBox(), br::ShapeRole::kCut});
    shapes.push_back({def->getLayer2Num(),
                      via->getLayer2BBox(),
                      br::ShapeRole::kViaPad,
                      pad_ndr});
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
      if (grid_only_) {
        gridOnlyChange(s, routedShapeType(s.role), net.id, true);
      } else {
        changeShape(s, routedShapeType(s.role), net.id, true, nullptr);
      }
    }
    // A net without routing has no input routing.
    net.input_present = !net.input_shapes.empty();
    net.input_removed = false;
    return;
  }
  br::FastGridBuilder::PendingRemovals own;
  br::FastGridBuilder::PendingRemovals& batch
      = pending != nullptr ? *pending : own;
  for (const RShape& s : net.input_shapes) {
    if (grid_only_) {
      gridOnlyChange(s, routedShapeType(s.role), net.id, false);
    } else {
      changeShape(s, routedShapeType(s.role), net.id, false, nullptr, &batch);
    }
  }
  if (pending == nullptr) {
    fast_builder_.flush(own);
  }
  net.input_present = false;
}

int BonnRouteDR::ndrType(const frNet* net) const
{
  if (!net->hasNDR()) {
    return 0;
  }
  const auto it = ndr_type_.find(net->getNondefaultRule());
  return it != ndr_type_.end() ? it->second : 0;
}

int BonnRouteDR::addNetData(frNet* net)
{
  NetData data;
  data.net = net;
  data.id = static_cast<int>(nets_.size());
  data.skip = true;
  if (params_.route_ndr_nets) {
    data.wire_type = ndrType(net);
  }
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
  // Nets with input routing need no guides: repair_antennas
  // -keep_detailed_wires leaves them without guides, and their routing
  // area comes from the kept routing (prepareNet).
  if (data.dead || data.net == nullptr) {
    return false;
  }
  const frNet* net = data.net;
  const auto type = net->getType();
  const size_t num_terms = net->getInstTerms().size() + net->getBTerms().size();
  // NDR nets with their wire type (route_ndr_nets), as in the first call.
  const bool ndr_routable = params_.route_ndr_nets && ndrType(net) > 0;
  return !(net->isSpecial() || net->isFake() || net->isFixed()
           || (net->hasNDR() && !ndr_routable) || type == odb::dbSigType::POWER
           || type == odb::dbSigType::GROUND || num_terms < 2
           || (!net->hasGuides() && !data.input_present));
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
      stats.shapes_inserted
          += static_cast<int64_t>(record.shapes.size() - before);
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
  // Nets whose routing fingerprint changed (or that BonnRoute routed in
  // the previous call) are compared shape by shape with the routing space;
  // only shapes that differ are removed or inserted.
  std::vector<int> reread;
  std::vector<uint64_t> fp(nets_.size(), 0);
  std::vector<char> seen(nets_.size(), 0);
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
    }
  }
  auto key = [](const RShape& r) {
    return std::make_tuple(r.layer,
                           r.rect.xMin(),
                           r.rect.yMin(),
                           r.rect.xMax(),
                           r.rect.yMax(),
                           static_cast<int>(r.role));
  };
  auto less = [&](const RShape& a, const RShape& b) { return key(a) < key(b); };
  // Per net: shapes in the space (old) and in the design (new), sorted.
  struct Delta
  {
    int id;
    std::vector<RShape> now;
    std::vector<RShape> remove;
    std::vector<RShape> insert;
  };
  std::vector<Delta> deltas;
  int64_t changed_shapes = 0;
  for (const int id : reread) {
    NetData& data = nets_[id];
    Delta d{id, {}, {}, {}};
    routingShapes(data.net, d.now);
    std::vector<RShape> old;
    if (data.input_present) {
      old = data.input_shapes;
    }
    old.insert(old.end(), data.shapes.begin(), data.shapes.end());
    std::ranges::sort(old, less);
    std::vector<RShape> now = d.now;
    std::ranges::sort(now, less);
    std::ranges::set_difference(old, now, std::back_inserter(d.remove), less);
    std::ranges::set_difference(now, old, std::back_inserter(d.insert), less);
    changed_shapes += static_cast<int64_t>(d.remove.size() + d.insert.size());
    deltas.push_back(std::move(d));
  }
  // Many changes: update the shape grid only and rebuild the fast grid once.
  const bool bulk = changed_shapes > 200000;
  for (Delta& d : deltas) {
    NetData& data = nets_[d.id];
    br::FastGridBuilder::PendingRemovals pending;
    for (const RShape& s : d.remove) {
      if (bulk) {
        gridOnlyChange(s, routedShapeType(s.role), data.id, false);
      } else {
        changeShape(
            s, routedShapeType(s.role), data.id, false, nullptr, &pending);
      }
    }
    if (!bulk) {
      fast_builder_.flush(pending);
    }
    for (const RShape& s : d.insert) {
      if (bulk) {
        gridOnlyChange(s, routedShapeType(s.role), data.id, true);
      } else {
        changeShape(s, routedShapeType(s.role), data.id, true, nullptr);
      }
    }
    stats.shapes_removed += static_cast<int64_t>(d.remove.size());
    stats.shapes_inserted += static_cast<int64_t>(d.insert.size());
    if (!d.remove.empty() || !d.insert.empty()) {
      stats.nets_reread++;
    }
    data.shapes.clear();
    data.wires.clear();
    data.vias.clear();
    data.routed = false;
    data.input_shapes = std::move(d.now);
    data.input_present = !data.input_shapes.empty();
    data.input_removed = false;
    data.routing_fp = fp[d.id];
  }
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
  grid_only_ = true;
  for (NetData& data : nets_) {
    routingShapes(data.net, data.input_shapes);
    if (!data.input_shapes.empty()) {
      setInputShapes(data, true);
    }
    data.routing_fp = routingFingerprint(data.net);
    data.terms_fp = termsFingerprint(data.net);
  }
  grid_only_ = false;
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
        changeShape(
            s, routedShapeType(s.role), data.id, false, nullptr, &pending);
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
  if (net.guide_tiles.empty()) {
    // No guides: the tiles of the net's routing and pins.
    frBlock* block = design_->getTopBlock();
    auto add = [&](const frLayerNum fr, const odb::Rect& r) {
      int z = fr >= 0 && fr < static_cast<int>(z_of_fr_.size()) ? z_of_fr_[fr]
                                                                : -1;
      if (z < 0) {
        z = fr < layers_.front().fr ? 0 : graph_.numLayers() - 1;
      }
      const odb::Point a = block->getGCellIdx(r.ll());
      const odb::Point b = block->getGCellIdx(r.ur());
      for (int x = a.x(); x <= b.x(); x++) {
        for (int y = a.y(); y <= b.y(); y++) {
          net.guide_tiles.push_back({z, x, y});
        }
      }
    };
    for (const RShape& sh : net.input_shapes) {
      add(sh.layer, sh.rect);
    }
    for (const RShape& sh : net.pins) {
      add(sh.layer, sh.rect);
    }
    std::ranges::sort(net.guide_tiles);
    net.guide_tiles.erase(std::ranges::unique(net.guide_tiles).begin(),
                          net.guide_tiles.end());
  }
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
        if (term.primary >= 0
            && term.primary < static_cast<int>(term.options.size())) {
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
      const std::vector<int> choice
          = choosePrimaries({&legal}, {net.id}, stats);
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

int BonnRouteDR::routingComponents(const frNet* net) const
{
  std::vector<RShape> shapes;
  routingShapes(net, shapes);
  // Objects: wires and patches one shape each, vias three; then pins.
  std::vector<int> object_of;
  int num_objects = 0;
  for (const auto& shape : net->getShapes()) {
    (void) shape;
    object_of.push_back(num_objects++);
  }
  for (const auto& via : net->getVias()) {
    (void) via;
    object_of.insert(object_of.end(), 3, num_objects++);
  }
  for (const auto& patch : net->getPatchWires()) {
    (void) patch;
    object_of.push_back(num_objects++);
  }
  struct Item
  {
    frLayerNum layer;
    odb::Rect rect;
    int object;
  };
  std::vector<Item> items;
  for (size_t k = 0; k < shapes.size(); k++) {
    items.push_back({shapes[k].layer, shapes[k].rect, object_of[k]});
  }
  const int first_pin = num_objects;
  std::vector<frRect> rects;
  for (const frInstTerm* it : net->getInstTerms()) {
    rects.clear();
    it->getShapes(rects);
    for (const frRect& r : rects) {
      items.push_back({r.getLayerNum(), r.getBBox(), num_objects});
    }
    num_objects++;
  }
  for (const frBTerm* bt : net->getBTerms()) {
    rects.clear();
    bt->getShapes(rects);
    for (const frRect& r : rects) {
      items.push_back({r.getLayerNum(), r.getBBox(), num_objects});
    }
    num_objects++;
  }
  std::vector<int> parent(num_objects);
  for (int k = 0; k < num_objects; k++) {
    parent[k] = k;
  }
  auto find = [&](int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };
  std::ranges::sort(items, [](const Item& a, const Item& b) {
    return std::make_pair(a.layer, a.rect.xMin())
           < std::make_pair(b.layer, b.rect.xMin());
  });
  for (size_t a = 0; a < items.size(); a++) {
    for (size_t b = a + 1; b < items.size() && items[b].layer == items[a].layer
                           && items[b].rect.xMin() <= items[a].rect.xMax();
         b++) {
      if (items[a].rect.intersects(items[b].rect)) {
        parent[find(items[a].object)] = find(items[b].object);
      }
    }
  }
  std::unordered_set<int> roots;
  for (int k = 0; k < num_objects; k++) {
    roots.insert(find(k));
  }
  (void) first_pin;
  return static_cast<int>(roots.size());
}

void BonnRouteDR::classifyInput(const NetData& net,
                                InputClassification& out) const
{
  out = InputClassification{};
  const frNet* fr = net.net;
  // Pieces and their shapes (as routingShapes builds them).
  std::vector<RShape> shapes;
  for (const auto& shape : fr->getShapes()) {
    InputPiece piece;
    piece.obj = shape.get();
    piece.kind = 0;
    out.pieces.push_back(std::move(piece));
  }
  for (const auto& via : fr->getVias()) {
    InputPiece piece;
    piece.obj = via.get();
    piece.kind = 1;
    out.pieces.push_back(std::move(piece));
  }
  for (const auto& patch : fr->getPatchWires()) {
    InputPiece piece;
    piece.obj = patch.get();
    piece.kind = 2;
    out.pieces.push_back(std::move(piece));
  }
  routingShapes(fr, shapes);
  {
    // routingShapes lists wires, then three shapes per via, then patches.
    size_t k = 0;
    for (InputPiece& piece : out.pieces) {
      const size_t n = piece.kind == 1 ? 3 : 1;
      piece.shapes.assign(shapes.begin() + k, shapes.begin() + k + n);
      k += n;
    }
  }
  const int num_pieces = static_cast<int>(out.pieces.size());
  // Pins in the order of net.terms (instance terms, then top-level terms).
  std::vector<std::vector<frRect>> pins;
  for (const frInstTerm* it : fr->getInstTerms()) {
    pins.emplace_back();
    it->getShapes(pins.back());
  }
  for (const frBTerm* bt : fr->getBTerms()) {
    pins.emplace_back();
    bt->getShapes(pins.back());
  }
  const int num_pins = static_cast<int>(pins.size());
  br::PieceGraph g;
  g.num_nodes = num_pieces + num_pins;
  g.is_pin.assign(g.num_nodes, 0);
  g.legal.assign(g.num_nodes, 1);
  g.is_short.assign(g.num_nodes, 0);
  for (int t = 0; t < num_pins; t++) {
    g.is_pin[num_pieces + t] = 1;
  }
  const int num_fr = static_cast<int>(z_of_fr_.size());
  for (int v = 0; v < num_pieces; v++) {
    const InputPiece& piece = out.pieces[v];
    for (const RShape& s : piece.shapes) {
      if (shapeConflicts(s, net.id)) {
        g.legal[v] = 0;
        out.any_illegal = true;
        break;
      }
    }
    if (piece.kind != 1) {
      const RShape& s = piece.shapes.front();
      const frLayer* layer = design_->getTech()->getLayer(s.layer);
      const int z = s.layer >= 0 && s.layer < num_fr ? z_of_fr_[s.layer] : -1;
      const int64_t pitch
          = z >= 0 ? layers_[z].pitch : static_cast<int64_t>(layer->getPitch());
      const int64_t length = std::max(s.rect.dx(), s.rect.dy());
      g.is_short[v] = length < params_.reserve_short_pitches * pitch;
    }
  }
  // Touching metal on a layer joins two nodes (sweep over x per layer).
  struct Item
  {
    odb::Rect rect;
    int node;
  };
  std::unordered_map<frLayerNum, std::vector<Item>> by_layer;
  for (int v = 0; v < num_pieces; v++) {
    for (const RShape& s : out.pieces[v].shapes) {
      by_layer[s.layer].push_back({s.rect, v});
    }
  }
  for (int t = 0; t < num_pins; t++) {
    for (const frRect& r : pins[t]) {
      by_layer[r.getLayerNum()].push_back({r.getBBox(), num_pieces + t});
    }
  }
  for (auto& [layer, items] : by_layer) {
    std::ranges::sort(items, [](const Item& a, const Item& b) {
      return a.rect.xMin() < b.rect.xMin();
    });
    for (size_t a = 0; a < items.size(); a++) {
      for (size_t b = a + 1;
           b < items.size() && items[b].rect.xMin() <= items[a].rect.xMax();
           b++) {
        if (items[a].node != items[b].node
            && items[a].rect.intersects(items[b].rect)) {
          g.edges.emplace_back(items[a].node, items[b].node);
        }
      }
    }
  }
  out.kept = br::keptPieces(g);
  out.reserved = br::reservedPieces(g, out.kept);
  const std::vector<int> label = br::keptComponents(g, out.kept);
  out.piece_label.assign(label.begin(), label.begin() + num_pieces);
  out.term_label.assign(label.begin() + num_pieces, label.end());
  out.complete
      = num_pins > 0 && std::ranges::all_of(out.term_label, [&](const int l) {
          return l == out.term_label.front();
        });
}

void BonnRouteDR::applyClassification(NetData& net,
                                      const InputClassification& c,
                                      ClassifyStats& stats)
{
  stats.nets++;
  // The net's input shapes are in the routing space; dropped and reserved
  // pieces leave it, reserved ones come back as reservations.
  br::FastGridBuilder::PendingRemovals pending;
  setInputShapes(net, false, &pending);
  fast_builder_.flush(pending);
  net.input_shapes.clear();
  net.input_reservation.clear();
  net.discount_runs.clear();
  const int num_fr = static_cast<int>(z_of_fr_.size());
  for (size_t v = 0; v < c.pieces.size(); v++) {
    const InputPiece& piece = c.pieces[v];
    stats.pieces++;
    if (c.kept[v]) {
      net.input_shapes.insert(
          net.input_shapes.end(), piece.shapes.begin(), piece.shapes.end());
      stats.kept++;
      continue;
    }
    if (!c.reserved[v]) {
      stats.dropped++;
      continue;
    }
    stats.reserved++;
    net.input_reservation.insert(
        net.input_reservation.end(), piece.shapes.begin(), piece.shapes.end());
    // Discounted runs need discount costs (enableIncrementalDiscounts);
    // otherwise the reservation only keeps the space for the net.
    if (piece.kind != 0 || graph_.wire_unit <= 1) {
      continue;
    }
    // Discount runs on the track the wire lies on (as reserveInputWires).
    const auto* seg = static_cast<const frPathSeg*>(piece.obj);
    if (seg->typeId() != frcPathSeg) {
      continue;
    }
    const frLayerNum fr = seg->getLayerNum();
    const int z = fr >= 0 && fr < num_fr ? z_of_fr_[fr] : -1;
    if (z < params_.reserve_min_layer || z < 0) {
      continue;
    }
    const auto [begin, end] = seg->getPoints();
    const br::TrackLayer& l = graph_.layers[z];
    const bool horizontal = begin.y() == end.y();
    if (horizontal != l.horizontal || begin == end) {
      continue;
    }
    const frCoord cross = horizontal ? begin.y() : begin.x();
    auto track = std::ranges::lower_bound(l.tracks, cross);
    if (track == l.tracks.end() || *track != cross) {
      continue;
    }
    const frCoord lo = horizontal ? std::min(begin.x(), end.x())
                                  : std::min(begin.y(), end.y());
    const frCoord hi = horizontal ? std::max(begin.x(), end.x())
                                  : std::max(begin.y(), end.y());
    const int jlo = static_cast<int>(std::ranges::lower_bound(l.positions, lo)
                                     - l.positions.begin());
    const int jhi = static_cast<int>(std::ranges::upper_bound(l.positions, hi)
                                     - l.positions.begin())
                    - 1;
    if (jhi > jlo) {
      net.discount_runs.push_back(
          {z, static_cast<int>(track - l.tracks.begin()), jlo, jhi});
      stats.discount_runs++;
    }
  }
  if (!net.input_shapes.empty()) {
    setInputShapes(net, true);
  }
  net.input_reserved = false;
  setInputReservation(net, true);
  net.keeps_input = true;
}

void BonnRouteDR::mergeKeptComponents(NetData& net,
                                      const InputClassification& c)
{
  // Kept pieces per component label.
  std::unordered_map<int, std::vector<size_t>> pieces_of;
  for (size_t v = 0; v < c.pieces.size(); v++) {
    if (c.kept[v]) {
      pieces_of[c.piece_label[v]].push_back(v);
    }
  }
  const int num_terms = static_cast<int>(net.terms.size());
  if (static_cast<int>(c.term_label.size()) != num_terms) {
    return;
  }
  std::unordered_map<int, int> merged_of;  // label -> index in out
  std::vector<Term> out;
  std::vector<Term> lone;
  for (int t = 0; t < num_terms; t++) {
    Term& term = net.terms[t];
    const int l = c.term_label[t];
    if (!pieces_of.contains(l)) {
      lone.push_back(std::move(term));
      continue;
    }
    // Reservations are re-made for the merged terminal.
    setReservation(term, net.id, false);
    auto it = merged_of.find(l);
    if (it == merged_of.end()) {
      merged_of.emplace(l, static_cast<int>(out.size()));
      out.push_back(std::move(term));
      continue;
    }
    Term& m = out[it->second];
    for (AccessOption& o : term.options) {
      m.options.push_back(std::move(o));
    }
    m.reservation.insert(
        m.reservation.end(), term.reservation.begin(), term.reservation.end());
  }
  // Track vertices on the kept pieces (cost 0): the kept routing is part of
  // the terminal (Fig. 9: new wires attach to the kept wires).
  const int num_fr = static_cast<int>(z_of_fr_.size());
  for (auto& [l, idx] : merged_of) {
    Term& m = out[idx];
    std::unordered_set<size_t> have;
    for (const AccessOption& o : m.options) {
      have.insert(o.vertex);
    }
    auto add_vertex = [&](const int z, const int i, const int j) {
      const size_t v = graph_.vertex(z, i, j);
      if (!have.insert(v).second) {
        return;
      }
      AccessOption opt;
      opt.vertex = v;
      opt.z = z;
      opt.cost = 0;
      m.options.push_back(std::move(opt));
    };
    for (const size_t p : pieces_of[l]) {
      const InputPiece& piece = c.pieces[p];
      if (piece.kind == 0 && piece.obj->typeId() == frcPathSeg) {
        const auto* seg = static_cast<const frPathSeg*>(piece.obj);
        const frLayerNum fr = seg->getLayerNum();
        const int z = fr >= 0 && fr < num_fr ? z_of_fr_[fr] : -1;
        if (z < 0) {
          continue;
        }
        const auto [b, e] = seg->getPoints();
        const br::TrackLayer& tl = graph_.layers[z];
        const bool horizontal = b == e ? tl.horizontal : b.y() == e.y();
        if (horizontal != tl.horizontal) {
          continue;
        }
        const frCoord cross = horizontal ? b.y() : b.x();
        auto track = std::ranges::lower_bound(tl.tracks, cross);
        if (track == tl.tracks.end() || *track != cross) {
          continue;
        }
        const int i = static_cast<int>(track - tl.tracks.begin());
        const frCoord lo
            = horizontal ? std::min(b.x(), e.x()) : std::min(b.y(), e.y());
        const frCoord hi
            = horizontal ? std::max(b.x(), e.x()) : std::max(b.y(), e.y());
        for (auto j = std::ranges::lower_bound(tl.positions, lo);
             j != tl.positions.end() && *j <= hi;
             ++j) {
          add_vertex(z, i, static_cast<int>(j - tl.positions.begin()));
        }
      } else if (piece.kind == 1) {
        const auto* via = static_cast<const frVia*>(piece.obj);
        const odb::Point o = via->getOrigin();
        for (const frLayerNum fr : {via->getViaDef()->getLayer1Num(),
                                    via->getViaDef()->getLayer2Num()}) {
          const int z = fr >= 0 && fr < num_fr ? z_of_fr_[fr] : -1;
          if (z < 0) {
            continue;
          }
          const br::TrackLayer& tl = graph_.layers[z];
          const frCoord cross = tl.horizontal ? o.y() : o.x();
          const frCoord along = tl.horizontal ? o.x() : o.y();
          auto track = std::ranges::lower_bound(tl.tracks, cross);
          auto pos = std::ranges::lower_bound(tl.positions, along);
          if (track != tl.tracks.end() && *track == cross
              && pos != tl.positions.end() && *pos == along) {
            add_vertex(z,
                       static_cast<int>(track - tl.tracks.begin()),
                       static_cast<int>(pos - tl.positions.begin()));
          }
        }
      }
    }
    m.chosen = -1;
    m.primary = -1;
    setReservation(m, net.id, true);
  }
  // The largest kept component is the source (term 0).
  std::ranges::stable_sort(out, [](const Term& a, const Term& b) {
    return a.options.size() > b.options.size();
  });
  for (Term& t : lone) {
    out.push_back(std::move(t));
  }
  net.terms = std::move(out);
  net.skip = std::ranges::any_of(
      net.terms, [](const Term& term) { return term.options.empty(); });
}

bool BonnRouteDR::prepareVictim(const int id)
{
  if (!incremental_active_) {
    return true;
  }
  NetData& net = nets_[id];
  if (net.dead || net.keeps_input || victim_changed_insts_ == nullptr) {
    return false;
  }
  return prepareNet(net, *victim_changed_insts_);
}

void BonnRouteDR::routeIncremental(std::vector<odb::Rect>& regions)
{
  utl::Timer timer;
  regions.clear();
  incremental_active_ = true;
  enableIncrementalDiscounts();
  // The paper's incremental runs use the general potential (Tables 4, 6).
  const int potential = params_.potential;
  params_.potential = params_.incremental_potential;
  SyncStats stats;
  syncInstances(stats);
  syncNets(stats);
  const double sync_s = timer.elapsed();
  victim_changed_insts_ = &stats.changed_insts;

  // Nets to route (Sec. 5.2): nets without routing (repair_antennas
  // destroyed it) are routed from their guides; nets with routing whose
  // pins changed, or with illegal pieces near changed instances, keep their
  // legal pin-to-pin paths and reserve their other legal pieces (Fig. 9).
  std::vector<int> work;
  std::vector<int> candidates;
  std::vector<char> terms_changed(nets_.size(), 0);
  int64_t unrouted = 0;
  int64_t pins_changed = 0;
  for (NetData& data : nets_) {
    if (!routableType(data)) {
      continue;
    }
    const uint64_t tfp = termsFingerprint(data.net);
    terms_changed[data.id] = tfp != data.terms_fp;
    data.terms_fp = tfp;
    if (!data.input_present) {
      work.push_back(data.id);
      unrouted++;
    } else if (terms_changed[data.id]) {
      candidates.push_back(data.id);
      pins_changed++;
    }
  }
  if (!stats.inst_regions.empty()) {
    // Routed nets whose shapes meet a changed instance: checked for pieces
    // that became illegal.
    for (NetData& data : nets_) {
      if (!data.input_present || terms_changed[data.id]
          || !routableType(data)) {
        continue;
      }
      const bool near
          = std::ranges::any_of(data.input_shapes, [&](const RShape& sh) {
              return std::ranges::any_of(
                  stats.inst_regions,
                  [&](const odb::Rect& r) { return r.intersects(sh.rect); });
            });
      if (near) {
        candidates.push_back(data.id);
      }
    }
  }
  if (logger_->debugCheck(utl::DRT, "incremental", 2)) {
    int reported = 0;
    for (const NetData& data : nets_) {
      if (data.dead || data.net == nullptr || reported >= 20) {
        continue;
      }
      if (terms_changed[data.id]
          || (!data.input_present && routableType(data))) {
        logger_->report(
            "  BonnRoute incremental net {}: routable {} input {} "
            "pins changed {} terms {}",
            data.net->getName(),
            routableType(data),
            data.input_present,
            static_cast<bool>(terms_changed[data.id]),
            data.net->getInstTerms().size() + data.net->getBTerms().size());
        reported++;
      }
    }
  }
  ClassifyStats cstats;
  int64_t complete = 0;
  std::vector<int> keep_work;
  for (const int id : candidates) {
    InputClassification c;
    classifyInput(nets_[id], c);
    if (c.complete && !(terms_changed[id] && c.any_illegal)) {
      // All pins joined by legal paths; illegal leftovers go to FlexDR.
      complete++;
      continue;
    }
    if (!c.complete && !terms_changed[id] && !c.any_illegal) {
      continue;
    }
    classified_.emplace(id, std::move(c));
    keep_work.push_back(id);
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
  for (const int id : keep_work) {
    NetData& data = nets_[id];
    if (!prepareNet(data, stats.changed_insts)) {
      // BonnRoute cannot route it: the input routing goes and FlexDR routes
      // the net from its guides.
      classified_.erase(id);
      setInputShapes(data, false);
      data.input_removed = true;
      not_prepared.push_back(id);
      continue;
    }
    const InputClassification& c = classified_.at(id);
    applyClassification(data, c, cstats);
    mergeKeptComponents(data, c);
    order.push_back(id);
  }
  const double access_s = timer.elapsed() - sync_s;
  std::ranges::stable_sort(
      order, [&](int a, int b) { return nets_[a].hpwl < nets_[b].hpwl; });
  const int64_t ripups_before = ripup_events_;
  std::vector<int> failed;
  routeParallel(order, failed);
  const int64_t first_failed = static_cast<int64_t>(failed.size());
  // Sec. 4.4 / arXiv:2111.06169 Sec. 5.2: the few nets of an incremental
  // call that fail get rip-up sequences in every profile; a net left to
  // FlexDR costs a track assignment and a FlexDR routing pass instead.
  if ((params_.ripup || params_.incremental_ripup) && !failed.empty()) {
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
  // Fig. 9 nets: the pieces that were not kept leave the frNet. A net that
  // could not be connected loses its input routing and goes to FlexDR when
  // it has guides; without guides (repair_antennas -keep_detailed_wires)
  // FlexDR cannot route it, so its input routing is restored unchanged.
  int64_t restored = 0;
  for (auto& [id, c] : classified_) {
    NetData& data = nets_[id];
    frNet* fr = data.net;
    if (data.routed) {
      for (size_t v = 0; v < c.pieces.size(); v++) {
        if (c.kept[v]) {
          continue;
        }
        auto* obj = const_cast<frBlockObject*>(c.pieces[v].obj);
        if (c.pieces[v].kind == 0) {
          fr->removeShape(static_cast<frShape*>(obj));
        } else if (c.pieces[v].kind == 1) {
          fr->removeVia(static_cast<frVia*>(obj));
        } else {
          fr->removePatchWire(static_cast<frShape*>(obj));
        }
      }
      fr->setModified(true);
    } else if (!fr->hasGuides()) {
      setReservations(data, false);
      setInputShapes(data, false);
      routingShapes(fr, data.input_shapes);
      setInputShapes(data, true);
      data.keeps_input = false;
      restored++;
    } else {
      setInputShapes(data, false);
      fr->clearRoutes();
      fr->setHasInitialRouting(false);
      fr->setModified(true);
      if (std::ranges::find(touched, id) == touched.end()) {
        touched.push_back(id);
      }
    }
  }
  writeNets();
  write_mask_.clear();
  // Debug check: nets written in this call whose routing and pins are not
  // one component.
  int64_t disconnected = 0;
  if (logger_->debugCheck(utl::DRT, "incremental", 1)) {
    for (const int id : touched) {
      const NetData& data = nets_[id];
      if (!data.routed) {
        continue;
      }
      const int components = routingComponents(data.net);
      if (components != 1) {
        disconnected++;
        if (disconnected <= 20) {
          logger_->report(
              "  BonnRoute incremental wrote net {} in {} components ({} "
              "wires, {} vias, {} patches, {} pins, input kept {}).",
              data.net->getName(),
              components,
              data.net->getShapes().size(),
              data.net->getVias().size(),
              data.net->getPatchWires().size(),
              data.net->getInstTerms().size() + data.net->getBTerms().size(),
              data.keeps_input);
        }
      }
    }
  }
  int64_t written = 0;
  int64_t left = 0;
  for (const int id : touched) {
    NetData& data = nets_[id];
    if (data.routed) {
      for (const RShape& s : data.shapes) {
        regions.push_back(s.rect);
      }
      if (data.keeps_input && data.input_present) {
        // Kept pieces stay; the new connections join them.
        data.input_shapes.insert(
            data.input_shapes.end(), data.shapes.begin(), data.shapes.end());
      } else {
        data.input_shapes = std::move(data.shapes);
      }
      data.input_present = true;
      written++;
    } else {
      // Not written (contact guard) or rip-up victim left unrouted:
      // FlexDR routes it from its guides.
      br::FastGridBuilder::PendingRemovals pending;
      for (const RShape& s : data.shapes) {
        changeShape(
            s, routedShapeType(s.role), data.id, false, nullptr, &pending);
      }
      fast_builder_.flush(pending);
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
        left++;
      }
    }
  }
  left += not_prepared.size();
  regions.insert(
      regions.end(), stats.inst_regions.begin(), stats.inst_regions.end());
  for (NetData& data : nets_) {
    if (data.prepared) {
      setReservations(data, false);
      data.prepared = false;
      data.failed = false;
      data.keeps_input = false;
      data.discount_runs.clear();
      data.input_reservation.clear();
    }
  }
  victim_changed_insts_ = nullptr;
  incremental_active_ = false;
  params_.potential = potential;
  classified_.clear();
  const double write_s = timer.elapsed() - sync_s - access_s - route_s;
  if (params_.verbose) {
    logger_->report(
        "BonnRoute incremental: sync_s={:.3f} access_s={:.3f} route_s={:.3f} "
        "write_s={:.3f} total_s={:.3f} instances_changed={} nets_reread={} "
        "nets_added={} nets_removed={} shapes_removed={} shapes_inserted={} "
        "nets_to_route={} unrouted={} pins_changed={} not_prepared={} "
        "first_pass_failed={} ripups={} written={} left_to_flexdr={} "
        "regions={} classified_nets={} complete_nets={} pieces={} kept={} "
        "reserved={} dropped={} discount_runs={} restored_unconnected={} "
        "disconnected_written={}",
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
        work.size() + keep_work.size(),
        unrouted,
        pins_changed,
        not_prepared.size(),
        first_failed,
        ripup_events_ - ripups_before,
        written,
        left,
        regions.size(),
        cstats.nets,
        complete,
        cstats.pieces,
        cstats.kept,
        cstats.reserved,
        cstats.dropped,
        cstats.discount_runs,
        restored,
        disconnected);
  }
}

}  // namespace drt
