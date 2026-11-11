// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BonnRouteDR.h"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "br/BrBlockageGrid.h"
#include "br/BrCore.h"
#include "br/BrGeometry.h"
#include "br/BrOffTrackSearch.h"
#include "br/BrPathCleanup.h"
#include "br/BrPinAccess.h"
#include "br/BrSameNetRules.h"
#include "db/obj/frAccess.h"
#include "db/obj/frBTerm.h"
#include "db/obj/frBlock.h"
#include "db/obj/frInst.h"
#include "db/obj/frInstTerm.h"
#include "db/obj/frNet.h"
#include "db/obj/frShape.h"
#include "db/obj/frVia.h"
#include "db/tech/frConstraint.h"
#include "db/tech/frLayer.h"
#include "db/tech/frTechObject.h"
#include "db/tech/frViaDef.h"
#include "drt-global.h"
#include "frDesign.h"
#include "frRegionQuery.h"
#include "odb/db.h"
#include "odb/dbTransform.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"
#include "utl/Logger.h"
#include "utl/timer.h"

namespace drt {

using br::Cost;
using br::Kind;

bool BonnRouteDRParams::set(const std::string& name, const double value)
{
  const int v = static_cast<int>(value);
  if (name == "via_cost") {
    via_cost = v;
  } else if (name == "jog_factor") {
    jog_factor = std::max(1, v);
  } else if (name == "ripup_cost") {
    ripup_cost = v;
  } else if (name == "max_ripups") {
    max_ripups = v;
  } else if (name == "max_access_jog") {
    max_access_jog = v;
  } else if (name == "max_access_options") {
    max_access_options = std::max(1, v);
  } else if (name == "pin_cost") {
    pin_cost = v;
  } else if (name == "layer_dilation") {
    layer_dilation = v;
  } else if (name == "retry_dilation") {
    retry_dilation = v;
  } else if (name == "use_intervals") {
    use_intervals = v != 0;
  } else if (name == "ripup") {
    ripup = v != 0;
  } else if (name == "interval_via_property") {
    interval_via_property = v != 0;
  } else if (name == "reservations") {
    reservations = v != 0;
  } else if (name == "discount_numerator") {
    discount_numerator = std::max(1, v);
  } else if (name == "discount_denominator") {
    discount_denominator = std::max(1, v);
  } else if (name == "reserve_min_layer") {
    reserve_min_layer = std::max(0, v);
  } else if (name == "potential") {
    potential = std::clamp(v, 0, 4);
  } else if (name == "potential_check") {
    potential_check = v != 0;
  } else if (name == "potential_detour_tiles") {
    potential_detour_tiles = std::max(0, v);
  } else if (name == "critical_nets") {
    critical_nets = v != 0;
  } else if (name == "partitions") {
    partitions = v != 0;
  } else if (name == "partition_regions") {
    partition_regions = std::max(0, v);
  } else if (name == "spreading") {
    spreading = v != 0;
  } else if (name == "spread_cost") {
    spread_cost = std::max(0, v);
  } else if (name == "spread_utilization") {
    spread_utilization = std::clamp(v, 0, 1000);
  } else if (name == "component_regions") {
    component_regions = v != 0;
  } else if (name == "ripup_depth") {
    ripup_depth = std::max(0, v);
  } else if (name == "ripup_passes") {
    ripup_passes = std::max(1, v);
  } else if (name == "offtrack_access") {
    offtrack_access = v != 0;
  } else if (name == "offtrack_radius") {
    offtrack_radius = std::max(0, v);
  } else if (name == "offtrack_regenerate") {
    offtrack_regenerate = std::max(0, v);
  } else if (name == "offtrack_tau") {
    offtrack_tau = std::max(0, v);
  } else if (name == "offtrack_notch_tau") {
    offtrack_notch_tau = v != 0;
  } else if (name == "access_classes") {
    access_classes = v != 0;
  } else if (name == "access_node_limit") {
    access_node_limit = std::max(1, v);
  } else if (name == "access_track_weight") {
    access_track_weight = std::max(0, v);
  } else if (name == "access_direction_weight") {
    access_direction_weight = std::max(0, v);
  } else if (name == "access_spread_weight") {
    access_spread_weight = std::max(0, v);
  } else if (name == "postprocess") {
    postprocess = v != 0;
  } else if (name == "verbose") {
    verbose = v != 0;
  } else {
    return false;
  }
  return true;
}

BonnRouteDR::BonnRouteDR(frDesign* design,
                         utl::Logger* logger,
                         RouterConfiguration* router_cfg,
                         const BonnRouteDRParams& params,
                         odb::dbBlock* block)
    : design_(design),
      logger_(logger),
      router_cfg_(router_cfg),
      params_(params),
      block_(block)
{
}

void BonnRouteDR::buildSpreadTiles()
{
  spread_tiles_.clear();
  spread_tile_count_ = 0;
  odb::dbGCellGrid* gcells
      = block_ != nullptr ? block_->getGCellGrid() : nullptr;
  if (gcells == nullptr) {
    logger_->warn(utl::DRT,
                  9104,
                  "BonnRoute: no global routing congestion in the database; "
                  "wire spreading is disabled.");
    return;
  }
  spread_tiles_.assign(
      graph_.numLayers(),
      std::vector<uint8_t>(static_cast<size_t>(tiles_x_) * tiles_y_, 0));
  for (int z = 0; z < graph_.numLayers(); z++) {
    odb::dbTechLayer* layer = layers_[z].layer->getDbLayer();
    for (int tx = 0; tx < tiles_x_; tx++) {
      for (int ty = 0; ty < tiles_y_; ty++) {
        const int cx = tile_x0_ + tx * tile_dx_ + tile_dx_ / 2;
        const int cy = tile_y0_ + ty * tile_dy_ + tile_dy_ / 2;
        const uint32_t gx = gcells->getXIdx(cx);
        const uint32_t gy = gcells->getYIdx(cy);
        const float capacity = gcells->getCapacity(layer, gx, gy);
        const float usage = gcells->getUsage(layer, gx, gy);
        if (capacity > 0
            && usage * 100 <= capacity * params_.spread_utilization) {
          spread_tiles_[z][static_cast<size_t>(tx) * tiles_y_ + ty] = 1;
          spread_tile_count_++;
        }
      }
    }
  }
}

BonnRouteDR::~BonnRouteDR() = default;

int BonnRouteDR::tileX(const frCoord x) const
{
  return std::clamp((x - tile_x0_) / tile_dx_, 0, tiles_x_ - 1);
}

int BonnRouteDR::tileY(const frCoord y) const
{
  return std::clamp((y - tile_y0_) / tile_dy_, 0, tiles_y_ - 1);
}

frCoord BonnRouteDR::spacing(const int z,
                             const frCoord width1,
                             const frCoord width2,
                             const frCoord prl) const
{
  return layers_[z].layer->getMinSpacingValue(width1, width2, prl, false);
}

void BonnRouteDR::buildTrackGraph()
{
  frTechObject* tech = design_->getTech();
  frBlock* block = design_->getTopBlock();
  const int num_fr = static_cast<int>(tech->getLayers().size());
  // Routing range: the layers the guides use, within drt's routing range.
  int lo = std::numeric_limits<int>::max();
  int hi = -1;
  for (const auto& net : block->getNets()) {
    for (const auto& guide : net->getGuides()) {
      lo = std::min(lo, static_cast<int>(guide->getLayerNum()));
      hi = std::max(hi, static_cast<int>(guide->getLayerNum()));
    }
  }
  lo = std::max(lo, router_cfg_->BOTTOM_ROUTING_LAYER);
  hi = std::min({hi, router_cfg_->TOP_ROUTING_LAYER, num_fr - 1});
  z_of_fr_.assign(num_fr, -1);
  cut_z_of_fr_.assign(num_fr, -1);
  layers_.clear();
  graph_.layers.clear();
  graph_.via_cost.clear();
  graph_.wire_unit = params_.reservations ? params_.discount_denominator : 1;
  for (int fr = lo; fr <= hi; fr++) {
    frLayer* layer = tech->getLayer(fr);
    if (layer->getType() != odb::dbTechLayerType::ROUTING) {
      continue;
    }
    LayerInfo info;
    info.fr = fr;
    info.layer = layer;
    info.width = static_cast<frCoord>(layer->getWidth());
    info.horizontal = layer->isHorizontal();
    const auto& eol = layer->getDrEolSpacingConstraint();
    info.eol_width = eol.eolWidth;
    info.eol_space = eol.eolSpace;
    br::TrackLayer tl;
    tl.horizontal = info.horizontal;
    for (const auto& tp : block->getTrackPatterns(fr)) {
      // X patterns (isHorizontal) hold x coordinates of vertical tracks.
      if (tp->isHorizontal() == info.horizontal) {
        continue;
      }
      for (frUInt4 k = 0; k < tp->getNumTracks(); k++) {
        tl.tracks.push_back(tp->getStartCoord()
                            + static_cast<frCoord>(k * tp->getTrackSpacing()));
      }
    }
    std::ranges::sort(tl.tracks);
    tl.tracks.erase(std::ranges::unique(tl.tracks).begin(), tl.tracks.end());
    if (tl.tracks.size() < 2) {
      continue;
    }
    info.pitch = layer->getPitch() > 0 ? static_cast<frCoord>(layer->getPitch())
                                       : tl.tracks[1] - tl.tracks[0];
    tl.jogs = !layer->isUnidirectional();
    tl.jog_factor = static_cast<Cost>(params_.jog_factor) * graph_.wire_unit;
    if (!layers_.empty() && layers_.back().horizontal == info.horizontal) {
      logger_->warn(utl::DRT,
                    9100,
                    "BonnRoute: layers {} and {} share a preferred direction; "
                    "detailed routing is skipped.",
                    layers_.back().layer->getName(),
                    layer->getName());
      layers_.clear();
      graph_.layers.clear();
      return;
    }
    z_of_fr_[fr] = static_cast<int>(layers_.size());
    layers_.push_back(info);
    graph_.layers.push_back(std::move(tl));
  }
  const int num_z = static_cast<int>(layers_.size());
  for (int z = 0; z < num_z; z++) {
    std::vector<br::Coord> pos;
    for (const int nz : {z - 1, z + 1}) {
      if (nz >= 0 && nz < num_z) {
        pos.insert(pos.end(),
                   graph_.layers[nz].tracks.begin(),
                   graph_.layers[nz].tracks.end());
      }
    }
    std::ranges::sort(pos);
    pos.erase(std::ranges::unique(pos).begin(), pos.end());
    graph_.layers[z].positions = std::move(pos);
  }
  for (int z = 0; z + 1 < num_z; z++) {
    LayerInfo& info = layers_[z];
    frLayer* cut = tech->getLayer(info.fr + 1);
    info.up_via = cut->getDefaultViaDef();
    if (info.up_via == nullptr) {
      logger_->warn(utl::DRT,
                    9101,
                    "BonnRoute: no default via on {}; detailed routing is "
                    "skipped.",
                    cut->getName());
      layers_.clear();
      graph_.layers.clear();
      return;
    }
    frVia via(info.up_via);
    info.up_pad = via.getLayer1BBox();
    info.up_cut = via.getCutBBox();
    layers_[z + 1].down_pad = via.getLayer2BBox();
    info.cut_spacing = cut->hasCutSpacing() ? cut->getCutSpacingValue() : 0;
    if (info.cut_spacing <= 0) {
      info.cut_spacing = info.up_cut.minDXDY();
    }
    cut_z_of_fr_[info.fr + 1] = z;
    graph_.via_cost.push_back(static_cast<Cost>(params_.via_cost)
                              * std::max(info.pitch, 1) * graph_.wire_unit);
  }
  die_ = block->getDieBox();
  const auto& gp = block->getGCellPatterns();
  tile_x0_ = gp[0].getStartCoord();
  tile_dx_ = std::max<frCoord>(1, gp[0].getSpacing());
  tiles_x_ = std::max(1, static_cast<int>(gp[0].getCount()));
  tile_y0_ = gp[1].getStartCoord();
  tile_dy_ = std::max<frCoord>(1, gp[1].getSpacing());
  tiles_y_ = std::max(1, static_cast<int>(gp[1].getCount()));
  graph_.tiles_x = tiles_x_;
  graph_.tiles_y = tiles_y_;
  graph_.finalize([this](br::Coord x) { return tileX(x); },
                  [this](br::Coord y) { return tileY(y); });
  fixed_tree_.clear();
  routed_tree_.clear();
  reservation_tree_.clear();
  fixed_tree_.resize(num_fr);
  routed_tree_.resize(num_fr);
  reservation_tree_.resize(num_fr);
  max_spacing_.assign(num_fr, 0);
  for (const LayerInfo& info : layers_) {
    max_spacing_[info.fr] = 4 * info.pitch;
    if (info.fr + 1 < num_fr) {
      max_spacing_[info.fr + 1] = 4 * info.pitch;
    }
  }
}

void BonnRouteDR::markShape(const RShape& shape,
                            const int mode,
                            const int owner,
                            const int delta)
{
  const odb::Rect& r = shape.rect;
  auto apply = [&](const int z,
                   const Kind kind,
                   const frCoord xlo,
                   const frCoord ylo,
                   const frCoord xhi,
                   const frCoord yhi) {
    br::forEachBlocked(
        graph_.layers[z], kind, xlo, ylo, xhi, yhi, [&](int i, int j) {
          const size_t v = graph_.vertex(z, i, j);
          if (mode == 0) {
            grid_.addFixed(v, kind);
          } else if (mode == 1) {
            grid_.addPin(v, kind, owner);
          } else if (mode == 2) {
            grid_.addRouted(v, kind, delta);
          } else {
            grid_.addReservation(v, kind, owner, delta);
          }
        });
    // The swept edges of wires and jogs (Sec. 3.6 edge exceptions).
    if (kind == br::kWire || kind == br::kJog) {
      const Kind edge = kind == br::kWire ? br::kWireEdge : br::kJogEdge;
      br::forEachBlocked(
          graph_.layers[z], edge, xlo, ylo, xhi, yhi, [&](int i, int j) {
            const size_t v = graph_.vertex(z, i, j);
            if (mode == 0) {
              grid_.addFixed(v, edge);
            } else if (mode == 1) {
              grid_.addPin(v, edge, owner);
            } else if (mode == 2) {
              grid_.addRouted(v, edge, delta);
            } else {
              grid_.addReservation(v, edge, owner, delta);
            }
          });
    }
  };
  // Origins p at which `model + p` is closer than (sx, sy) to r.
  auto block_model = [&](const int z,
                         const Kind kind,
                         const odb::Rect& model,
                         const frCoord sx,
                         const frCoord sy) {
    apply(z,
          kind,
          r.xMin() - sx - model.xMax() + 1,
          r.yMin() - sy - model.yMax() + 1,
          r.xMax() + sx - model.xMin() - 1,
          r.yMax() + sy - model.yMin() - 1);
  };
  const int num_fr = static_cast<int>(z_of_fr_.size());
  if (shape.layer < 0 || shape.layer >= num_fr) {
    return;
  }
  const int z = z_of_fr_[shape.layer];
  if (z >= 0) {
    const LayerInfo& info = layers_[z];
    const frCoord ws = r.minDXDY();
    const frCoord prl = r.maxDXDY();
    const frCoord w = info.width;
    const frCoord half = w / 2;
    const odb::Rect wire_model(-half, -half, w - half, w - half);
    const frCoord s_wire = spacing(z, ws, w, prl);
    // Wires are assumed to end in preferred direction (Sec. 3.1).
    const frCoord s_along = (info.eol_width > 0 && w < info.eol_width)
                                ? std::max(s_wire, info.eol_space)
                                : s_wire;
    block_model(z,
                br::kWire,
                wire_model,
                info.horizontal ? s_along : s_wire,
                info.horizontal ? s_wire : s_along);
    block_model(z, br::kJog, wire_model, s_wire, s_wire);
    if (z + 1 < static_cast<int>(layers_.size())) {
      const frCoord s = spacing(z, ws, info.up_pad.minDXDY(), prl);
      block_model(z, br::kViaUpPad, info.up_pad, s, s);
    }
    if (z > 0) {
      const frCoord s = spacing(z, ws, info.down_pad.minDXDY(), prl);
      block_model(z, br::kViaDownPad, info.down_pad, s, s);
    }
    return;
  }
  const int cz = cut_z_of_fr_[shape.layer];
  if (cz >= 0) {
    const LayerInfo& info = layers_[cz];
    block_model(
        cz, br::kCutUp, info.up_cut, info.cut_spacing, info.cut_spacing);
  }
}

void BonnRouteDR::buildNets()
{
  nets_.clear();
  net_index_.clear();
  for (const auto& uptr : design_->getTopBlock()->getNets()) {
    frNet* net = uptr.get();
    NetData data;
    data.net = net;
    data.id = static_cast<int>(nets_.size());
    const auto type = net->getType();
    const size_t num_terms
        = net->getInstTerms().size() + net->getBTerms().size();
    data.skip = net->isSpecial() || net->isFake() || net->isFixed()
                || net->hasInitialRouting() || net->hasNDR()
                || type == odb::dbSigType::POWER
                || type == odb::dbSigType::GROUND || num_terms < 2
                || !net->hasGuides();
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
  }
}

void BonnRouteDR::buildFixed()
{
  frRegionQuery* rq = design_->getRegionQuery();
  const odb::Rect die = design_->getTopBlock()->getDieBox();
  const int num_fr = static_cast<int>(z_of_fr_.size());
  for (int fr = 0; fr < num_fr; fr++) {
    if (z_of_fr_[fr] < 0 && cut_z_of_fr_[fr] < 0) {
      continue;
    }
    frRegionQuery::Objects<frBlockObject> result;
    rq->query(die, fr, result);
    std::vector<std::pair<odb::Rect, int>> entries;
    entries.reserve(result.size());
    for (const auto& [rect, obj] : result) {
      frNet* owner_net = nullptr;
      if (obj->typeId() == frcInstTerm) {
        owner_net = static_cast<frInstTerm*>(obj)->getNet();
      } else if (obj->typeId() == frcBTerm) {
        owner_net = static_cast<frBTerm*>(obj)->getNet();
      }
      int owner = -1;
      if (owner_net != nullptr) {
        auto it = net_index_.find(owner_net);
        if (it != net_index_.end()) {
          owner = it->second;
        }
      }
      markShape(
          {static_cast<frLayerNum>(fr), rect}, owner >= 0 ? 1 : 0, owner, 0);
      entries.emplace_back(rect, owner);
    }
    fixed_tree_[fr] = RTree<int>(entries.begin(), entries.end());
  }
}

bool BonnRouteDR::treeConflicts(const RTree<int>& tree,
                                const RShape& shape,
                                const int net,
                                const std::unordered_set<int>* skip) const
{
  const int num_fr = static_cast<int>(z_of_fr_.size());
  if (shape.layer < 0 || shape.layer >= num_fr) {
    return false;
  }
  const int z = z_of_fr_[shape.layer];
  const int cz = cut_z_of_fr_[shape.layer];
  if (z < 0 && cz < 0) {
    return false;
  }
  odb::Rect query = shape.rect;
  query.bloat(max_spacing_[shape.layer], query);
  for (auto it = tree.qbegin(bgi::intersects(query)); it != tree.qend(); ++it) {
    if ((it->second == net && net >= 0)
        || (skip != nullptr && skip->contains(it->second))) {
      continue;
    }
    const odb::Rect& other = it->first;
    frCoord s;
    if (z >= 0) {
      const frCoord prl = std::min(other.maxDXDY(), shape.rect.maxDXDY());
      s = spacing(z, other.minDXDY(), shape.rect.minDXDY(), prl);
    } else {
      s = layers_[cz].cut_spacing;
    }
    const int64_t dx = std::max({0,
                                 other.xMin() - shape.rect.xMax(),
                                 shape.rect.xMin() - other.xMax()});
    const int64_t dy = std::max({0,
                                 other.yMin() - shape.rect.yMax(),
                                 shape.rect.yMin() - other.yMax()});
    if (dx * dx + dy * dy < static_cast<int64_t>(s) * s) {
      return true;
    }
  }
  return false;
}

bool BonnRouteDR::shapeConflicts(const RShape& shape,
                                 const int net,
                                 const Worker* worker) const
{
  if (shape.layer < 0 || shape.layer >= static_cast<int>(z_of_fr_.size())) {
    return false;
  }
  const bool local = worker != nullptr && !worker->local_routed.empty();
  return treeConflicts(fixed_tree_[shape.layer], shape, net)
         || treeConflicts(routed_tree_[shape.layer], shape, net)
         || (local
             && treeConflicts(worker->local_routed[shape.layer], shape, net))
         || treeConflicts(reservation_tree_[shape.layer],
                          shape,
                          net,
                          worker != nullptr ? &worker->released : nullptr);
}

BonnRouteDR::RShape BonnRouteDR::wireShape(const OutWire& wire) const
{
  const int z = z_of_fr_[wire.layer];
  const frCoord half = layers_[z].width / 2;
  odb::Rect r(std::min(wire.begin.x(), wire.end.x()) - half,
              std::min(wire.begin.y(), wire.end.y()) - half,
              std::max(wire.begin.x(), wire.end.x()) + half,
              std::max(wire.begin.y(), wire.end.y()) + half);
  return {wire.layer, r};
}

void BonnRouteDR::viaShapes(const OutVia& via,
                            std::vector<RShape>& shapes) const
{
  frVia v(via.def, via.point);
  shapes.push_back({via.def->getLayer1Num(), v.getLayer1BBox()});
  shapes.push_back({via.def->getCutLayerNum(), v.getCutBBox()});
  shapes.push_back({via.def->getLayer2Num(), v.getLayer2BBox()});
}

void BonnRouteDR::shapesOf(const NetData& net,
                           std::vector<RShape>& shapes) const
{
  shapes.clear();
  for (const OutWire& w : net.wires) {
    shapes.push_back(wireShape(w));
  }
  for (const OutVia& v : net.vias) {
    viaShapes(v, shapes);
  }
}

void BonnRouteDR::buildAccess(NetData& net)
{
  struct Ap
  {
    const frAccessPoint* ap;
    odb::Point point;
    bool preferred;
  };
  auto add_term = [&](const std::vector<Ap>& aps) {
    Term term;
    for (const Ap& a : aps) {
      const frLayerNum la = a.ap->getLayerNum();
      const int num_fr = static_cast<int>(z_of_fr_.size());
      int z = -1;
      const frViaDef* via_def = nullptr;
      if (la >= 0 && la < num_fr && z_of_fr_[la] >= 0) {
        z = z_of_fr_[la];
      } else if (la + 2 < num_fr && z_of_fr_[la + 2] >= 0
                 && a.ap->hasAccess(frDirEnum::U) && a.ap->hasViaDef()) {
        z = z_of_fr_[la + 2];
        via_def = a.ap->getViaDefs(1).front();
      }
      if (z < 0) {
        continue;
      }
      const br::TrackLayer& l = graph_.layers[z];
      const LayerInfo& info = layers_[z];
      const frCoord c = l.horizontal ? a.point.y() : a.point.x();
      const frCoord along = l.horizontal ? a.point.x() : a.point.y();
      auto xy = [&](frCoord al, frCoord cr) {
        return l.horizontal ? odb::Point(al, cr) : odb::Point(cr, al);
      };
      const int it = static_cast<int>(std::ranges::lower_bound(l.tracks, c)
                                      - l.tracks.begin());
      const int jt = static_cast<int>(
          std::ranges::lower_bound(l.positions, along) - l.positions.begin());
      for (const int i : {it - 1, it}) {
        if (i < 0 || i >= l.numTracks()) {
          continue;
        }
        const frCoord jog = std::abs(l.tracks[i] - c);
        if (jog > params_.max_access_jog * info.pitch || (jog > 0 && !l.jogs)) {
          continue;
        }
        for (const int j : {jt - 1, jt}) {
          if (j < 0 || j >= l.numPositions()) {
            continue;
          }
          AccessOption opt;
          opt.z = z;
          opt.vertex = graph_.vertex(z, i, j);
          const frCoord run = std::abs(l.positions[j] - along);
          if (jog > 0) {
            opt.wires.push_back(
                {info.fr, xy(along, c), xy(along, l.tracks[i])});
          }
          if (run > 0) {
            opt.wires.push_back({info.fr,
                                 xy(along, l.tracks[i]),
                                 xy(l.positions[j], l.tracks[i])});
          }
          if (via_def != nullptr) {
            opt.vias.push_back({via_def, a.point});
          }
          if (via_def == nullptr) {
            opt.planar = true;
            opt.ap = a.point;
            opt.ap_layer = info.fr;
          }
          // The same stub as route nodes (pin -> vertex).
          if (via_def != nullptr) {
            opt.route.push_back({la, {a.point.x(), a.point.y()}});
          }
          opt.route.push_back({info.fr, {a.point.x(), a.point.y()}});
          for (const OutWire& w : opt.wires) {
            opt.route.push_back({info.fr, {w.end.x(), w.end.y()}});
          }
          opt.cost = (via_def != nullptr
                          ? (z > 0 ? graph_.via_cost[z - 1]
                                   : static_cast<Cost>(params_.via_cost)
                                         * info.pitch * unit())
                          : 0)
                     + (static_cast<Cost>(params_.jog_factor) * jog + run
                        + (a.preferred ? 0
                                       : static_cast<Cost>(params_.pin_cost)
                                             * info.pitch))
                           * unit();
          // Off-track stub against fixed shapes and other nets' pins.
          std::vector<RShape> stub;
          for (const OutWire& w : opt.wires) {
            stub.push_back(wireShape(w));
          }
          for (const OutVia& v : opt.vias) {
            std::vector<RShape> vs;
            viaShapes(v, vs);
            // The pin layer and cut were checked by pin access.
            stub.push_back(vs.back());
          }
          bool ok = true;
          for (const RShape& s : stub) {
            if (shapeConflicts(s, net.id)) {
              ok = false;
              break;
            }
          }
          if (ok) {
            term.options.push_back(std::move(opt));
          }
        }
      }
    }
    std::ranges::stable_sort(term.options, [](const auto& a, const auto& b) {
      return a.cost < b.cost;
    });
    // One option per vertex: the cheapest.
    std::vector<AccessOption> unique;
    std::unordered_set<size_t> seen;
    for (AccessOption& o : term.options) {
      if (seen.insert(o.vertex).second) {
        unique.push_back(std::move(o));
      }
      if (static_cast<int>(unique.size()) >= params_.max_access_options) {
        break;
      }
    }
    term.options = std::move(unique);
    if (term.options.empty()) {
      net.skip = true;
    }
    net.terms.push_back(std::move(term));
  };

  for (frInstTerm* it : net.net->getInstTerms()) {
    frInst* inst = it->getInst();
    const int pa_idx = inst->getPinAccessIdx();
    const odb::dbTransform xform = inst->getNoRotationTransform();
    const auto& prefs = it->getAccessPoints();
    std::vector<Ap> aps;
    int pin_idx = 0;
    for (const auto& pin : it->getTerm()->getPins()) {
      const frAccessPoint* pref
          = pin_idx < static_cast<int>(prefs.size()) ? prefs[pin_idx] : nullptr;
      pin_idx++;
      if (!pin->hasPinAccess() || pa_idx < 0
          || pa_idx >= pin->getNumPinAccess()) {
        continue;
      }
      for (const auto& ap : pin->getPinAccess(pa_idx)->getAccessPoints()) {
        odb::Point p = ap->getPoint();
        xform.apply(p);
        aps.push_back({ap.get(), p, ap.get() == pref});
      }
    }
    std::ranges::stable_sort(aps, [](const Ap& a, const Ap& b) {
      return a.preferred > b.preferred;
    });
    add_term(aps);
  }
  for (frBTerm* bt : net.net->getBTerms()) {
    std::vector<Ap> aps;
    for (const auto& pin : bt->getPins()) {
      if (!pin->hasPinAccess()) {
        continue;
      }
      for (const auto& ap : pin->getPinAccess(0)->getAccessPoints()) {
        aps.push_back({ap.get(), ap->getPoint(), true});
      }
    }
    add_term(aps);
  }
}

void BonnRouteDR::reserveAccess(NetData& net)
{
  // Select conflict-free primary paths from the access catalogues. Reserve
  // the complete geometry, not just the top pad of its access via.
  for (Term& term : net.terms) {
    for (const AccessOption& option : term.options) {
      std::vector<RShape> shapes;
      accessShapes(option, shapes);
      if (std::ranges::any_of(shapes, [&](const RShape& shape) {
            return shapeConflicts(shape, net.id);
          })) {
        continue;
      }
      term.reservation = std::move(shapes);
      setReservation(term, net.id, true);
      break;
    }
  }
}

// ---------------------------------------------------------------------------
// Off-track wiring and pin access (Sec. 3.7, 3.8, 4.3, 4.4); interpretation
// choices are documented in OffTrack.md.

namespace {

br::Box toBox(const odb::Rect& r)
{
  return {r.xMin(), r.yMin(), r.xMax(), r.yMax()};
}

odb::Rect toRect(const br::Box& b)
{
  return odb::Rect(b.xlo, b.ylo, b.xhi, b.yhi);
}

br::GridPoint toGrid(const odb::Point& p)
{
  return {p.x(), p.y()};
}

odb::Point toPoint(const br::GridPoint& p)
{
  return odb::Point(p.x, p.y);
}

// The lattice isometry of an odb transform (images of three points).
br::Isometry isometryOf(const odb::dbTransform& t)
{
  odb::Point o(0, 0);
  odb::Point ex(1, 0);
  odb::Point ey(0, 1);
  t.apply(o);
  t.apply(ex);
  t.apply(ey);
  br::Isometry m;
  m.m00 = ex.x() - o.x();
  m.m10 = ex.y() - o.y();
  m.m01 = ey.x() - o.x();
  m.m11 = ey.y() - o.y();
  m.dx = o.x();
  m.dy = o.y();
  return m;
}

}  // namespace

// The diff-net legality adapter (br::DiffNetLegality) over the current shape
// index: shapeConflicts for the track layers and their cuts; the same
// spacing semantics over low_tree_ for the pin layer and cut below the
// lowest track layer, which the shape index does not cover. The shape grid
// and distance rule checking module of Sec. 3.3-3.4 replace this class.
class BonnRouteDR::Legality final : public br::DiffNetLegality
{
 public:
  Legality(const BonnRouteDR& dr, const Worker* worker)
      : dr_(dr), worker_(worker)
  {
  }

  bool allowed(const int layer,
               const br::Box& box,
               const int net) const override
  {
    if (!valid(layer)) {
      return true;
    }
    const RShape shape{static_cast<frLayerNum>(layer), toRect(box)};
    if (covered(layer)) {
      return !dr_.shapeConflicts(shape, net, worker_);
    }
    bool hit = false;
    forEachForeign(layer, shape.rect, net, [&](const odb::Rect& other) {
      hit = violates(layer, other, shape.rect);
      return hit;
    });
    return !hit;
  }

  void blockages(const int layer,
                 const br::Box& window,
                 const int net,
                 const br::Box& model,
                 std::vector<br::Box>& out) const override
  {
    if (!valid(layer)) {
      return;
    }
    const odb::Rect m = toRect(model);
    odb::Rect query = toRect(window);
    query = odb::Rect(query.xMin() + m.xMin(),
                      query.yMin() + m.yMin(),
                      query.xMax() + m.xMax(),
                      query.yMax() + m.yMax());
    forEachForeign(layer, query, net, [&](const odb::Rect& o) {
      const frCoord s = spacingOf(layer, o.minDXDY(), m.minDXDY(), o.maxDXDY());
      out.push_back({o.xMin() - s - m.xMax(),
                     o.yMin() - s - m.yMax(),
                     o.xMax() + s - m.xMin(),
                     o.yMax() + s - m.yMin()});
      return false;
    });
  }

  bool conflict(const int layer,
                const br::Box& a,
                const br::Box& b) const override
  {
    return valid(layer) && violates(layer, toRect(a), toRect(b));
  }

 private:
  bool valid(const int layer) const
  {
    return layer >= 0 && layer < static_cast<int>(dr_.z_of_fr_.size());
  }
  bool covered(const int layer) const
  {
    return dr_.z_of_fr_[layer] >= 0 || dr_.cut_z_of_fr_[layer] >= 0;
  }
  // Spacing as treeConflicts applies it (width, width, parallel run).
  frCoord spacingOf(const int layer,
                    const frCoord w1,
                    const frCoord w2,
                    const frCoord prl) const
  {
    const int z = dr_.z_of_fr_[layer];
    if (z >= 0) {
      return dr_.spacing(z, w1, w2, prl);
    }
    const int cz = dr_.cut_z_of_fr_[layer];
    if (cz >= 0) {
      return dr_.layers_[cz].cut_spacing;
    }
    frLayer* l = dr_.design_->getTech()->getLayer(layer);
    if (l->getType() == odb::dbTechLayerType::ROUTING) {
      return l->getMinSpacingValue(w1, w2, prl, false);
    }
    return dr_.low_spacing_[layer];
  }
  bool violates(const int layer, const odb::Rect& a, const odb::Rect& b) const
  {
    const frCoord prl = std::min(a.maxDXDY(), b.maxDXDY());
    const frCoord s = spacingOf(layer, a.minDXDY(), b.minDXDY(), prl);
    const int64_t dx = std::max({0, a.xMin() - b.xMax(), b.xMin() - a.xMax()});
    const int64_t dy = std::max({0, a.yMin() - b.yMax(), b.yMin() - a.yMax()});
    return dx * dx + dy * dy < static_cast<int64_t>(s) * s;
  }
  frCoord halo(const int layer) const
  {
    return std::max(dr_.max_spacing_[layer], 4 * dr_.layers_[0].pitch);
  }
  // Calls fn(rect) for the shapes of other nets near `rect` (fixed shapes,
  // pins, routed shapes, reservations of other nets) with the skip rules of
  // shapeConflicts; stops when fn returns true.
  template <typename F>
  void forEachForeign(const int layer,
                      const odb::Rect& rect,
                      const int net,
                      F&& fn) const
  {
    odb::Rect query = rect;
    query.bloat(halo(layer), query);
    auto scan = [&](const RTree<int>& tree,
                    const std::unordered_set<int>* skip) {
      for (auto it = tree.qbegin(bgi::intersects(query)); it != tree.qend();
           ++it) {
        if ((it->second == net && net >= 0)
            || (skip != nullptr && skip->contains(it->second))) {
          continue;
        }
        if (fn(it->first)) {
          return true;
        }
      }
      return false;
    };
    if (scan(dr_.fixed_tree_[layer], nullptr)
        || (layer < static_cast<int>(dr_.low_tree_.size())
            && scan(dr_.low_tree_[layer], nullptr))
        || scan(dr_.routed_tree_[layer], nullptr)) {
      return;
    }
    if (worker_ != nullptr && !worker_->local_routed.empty()
        && scan(worker_->local_routed[layer], nullptr)) {
      return;
    }
    scan(dr_.reservation_tree_[layer],
         worker_ != nullptr ? &worker_->released : nullptr);
  }

  const BonnRouteDR& dr_;
  const Worker* worker_;
};

void BonnRouteDR::AccessStats::add(const AccessStats& o)
{
  pins += o.pins;
  searches += o.searches;
  paths += o.paths;
  empty += o.empty;
  fallback += o.fallback;
  drc_rejected += o.drc_rejected;
  same_net_rejected += o.same_net_rejected;
  grid_vertices += o.grid_vertices;
  settled += o.settled;
  search_seconds += o.search_seconds;
  instances += o.instances;
  classes += o.classes;
  class_reused += o.class_reused;
  reselected += o.reselected;
  mapping_failed += o.mapping_failed;
  selections += o.selections;
  selection_nodes += o.selection_nodes;
  selection_optimal += o.selection_optimal;
  primaries += o.primaries;
}

void BonnRouteDR::buildSameNetRules()
{
  const int num_z = static_cast<int>(layers_.size());
  same_net_rules_.assign(num_z, br::SameNetLayerRules());
  tau_.assign(num_z, 1);
  cleanup_layers_.assign(z_of_fr_.size(), br::CleanupLayer());
  for (int z = 0; z < num_z; z++) {
    const LayerInfo& info = layers_[z];
    frLayer* layer = info.layer;
    br::SameNetLayerRules rules;
    rules.width = info.width;
    frCoord step = 0;
    if (const frMinStepConstraint* c = layer->getMinStepConstraint()) {
      step = std::max(step, c->getMinStepLength());
    }
    for (const frLef58MinStepConstraint* c :
         layer->getLef58MinStepConstraints()) {
      step = std::max(step, c->getMinStepLength());
    }
    rules.min_step_length = step;
    rules.notch_spacing = spacing(z, info.width, info.width, 0);
    if (frAreaConstraint* area = layer->getAreaConstraint()) {
      rules.min_area = area->getMinArea();
    }
    const br::MinimumSegmentLength m = br::minimumSegmentLength(
        rules, params_.offtrack_notch_tau, params_.offtrack_tau);
    same_net_rules_[z] = rules;
    tau_[z] = m.tau;
    br::CleanupLayer& cl = cleanup_layers_[info.fr];
    cl.tau = m.tau;
    cl.half_width = info.width / 2;
    cl.rules = rules;
    cl.horizontal = info.horizontal;
    if (params_.verbose && (params_.offtrack_access || params_.postprocess)) {
      logger_->report(
          "BonnRoute same-net rules {}: width {}, min step {}, notch {}, "
          "min area {}, tau {} (short edge {}, notch {}).",
          layer->getName(),
          rules.width,
          rules.min_step_length,
          rules.notch_spacing,
          rules.min_area,
          m.tau,
          m.from_short_edge,
          m.from_notch);
    }
  }
}

void BonnRouteDR::buildLowLayerIndex()
{
  frTechObject* tech = design_->getTech();
  frRegionQuery* rq = design_->getRegionQuery();
  const odb::Rect die = design_->getTopBlock()->getDieBox();
  const int num_fr = static_cast<int>(z_of_fr_.size());
  low_tree_.assign(num_fr, RTree<int>());
  low_spacing_.assign(num_fr, 0);
  const frLayerNum bottom = layers_.front().fr;
  for (const frLayerNum fr : {bottom - 2, bottom - 1}) {
    if (fr < 0 || fr >= num_fr) {
      continue;
    }
    frLayer* layer = tech->getLayer(fr);
    const bool routing = layer->getType() == odb::dbTechLayerType::ROUTING;
    const bool cut = layer->getType() == odb::dbTechLayerType::CUT;
    if ((fr == bottom - 2 && !routing) || (fr == bottom - 1 && !cut)) {
      continue;
    }
    if (cut) {
      low_spacing_[fr]
          = layer->hasCutSpacing() ? layer->getCutSpacingValue() : 0;
      if (low_spacing_[fr] <= 0 && layer->getDefaultViaDef() != nullptr) {
        frVia via(layer->getDefaultViaDef());
        low_spacing_[fr] = via.getCutBBox().minDXDY();
      }
    }
    frRegionQuery::Objects<frBlockObject> result;
    rq->query(die, fr, result);
    std::vector<std::pair<odb::Rect, int>> entries;
    entries.reserve(result.size());
    for (const auto& [rect, obj] : result) {
      frNet* owner_net = nullptr;
      if (obj->typeId() == frcInstTerm) {
        owner_net = static_cast<frInstTerm*>(obj)->getNet();
      } else if (obj->typeId() == frcBTerm) {
        owner_net = static_cast<frBTerm*>(obj)->getNet();
      }
      int owner = -1;
      if (owner_net != nullptr) {
        auto it = net_index_.find(owner_net);
        if (it != net_index_.end()) {
          owner = it->second;
        }
      }
      entries.emplace_back(rect, owner);
    }
    low_tree_[fr] = RTree<int>(entries.begin(), entries.end());
  }
}

const frViaDef* BonnRouteDR::viaDefAbove(const frLayerNum lower) const
{
  const int num_fr = static_cast<int>(z_of_fr_.size());
  if (lower < 0 || lower >= num_fr) {
    return nullptr;
  }
  const int z = z_of_fr_[lower];
  if (z >= 0) {
    return z + 1 < static_cast<int>(layers_.size()) ? layers_[z].up_via
                                                    : nullptr;
  }
  if (lower + 1 < num_fr) {
    frLayer* cut = design_->getTech()->getLayer(lower + 1);
    if (cut->getType() == odb::dbTechLayerType::CUT) {
      return cut->getDefaultViaDef();
    }
  }
  return nullptr;
}

bool BonnRouteDR::pinShapes(const Term& term,
                            frLayerNum& layer,
                            std::vector<odb::Rect>& rects) const
{
  rects.clear();
  odb::Rect box;
  const frBlockObject* owner = nullptr;
  if (term.inst_term != nullptr) {
    box = term.inst_term->getBBox();
    owner = term.inst_term;
  } else if (term.bterm != nullptr) {
    box = term.bterm->getBBox();
    owner = term.bterm;
  } else {
    return false;
  }
  frRegionQuery* rq = design_->getRegionQuery();
  frTechObject* tech = design_->getTech();
  const int num_fr = static_cast<int>(z_of_fr_.size());
  const frLayerNum bottom = layers_.front().fr;
  for (frLayerNum fr = std::max(0, bottom - 2); fr < num_fr; fr++) {
    const bool usable = z_of_fr_[fr] >= 0
                        || (fr == bottom - 2
                            && tech->getLayer(fr)->getType()
                                   == odb::dbTechLayerType::ROUTING);
    if (!usable) {
      continue;
    }
    frRegionQuery::Objects<frBlockObject> result;
    rq->query(box, fr, result);
    for (const auto& [rect, obj] : result) {
      if (obj == owner) {
        rects.push_back(rect);
      }
    }
    if (!rects.empty()) {
      layer = fr;
      return true;
    }
  }
  return false;
}

bool BonnRouteDR::optionLegal(const AccessOption& option,
                              const int net,
                              const Worker* worker) const
{
  if (!grid_.usable(option.vertex, br::kWire, net)) {
    return false;
  }
  const Legality legality(*this, worker);
  std::vector<RShape> shapes;
  accessShapes(option, shapes);
  return std::ranges::all_of(shapes, [&](const RShape& s) {
    return legality.allowed(s.layer, toBox(s.rect), net);
  });
}

void BonnRouteDR::offTrackCatalogue(const int net,
                                    const Term& term,
                                    const Worker* worker,
                                    std::vector<AccessOption>& out,
                                    AccessStats& stats) const
{
  out.clear();
  stats.pins++;
  frLayerNum pin_layer = -1;
  std::vector<odb::Rect> pin;
  if (!pinShapes(term, pin_layer, pin)) {
    stats.empty++;
    return;
  }
  const int num_fr = static_cast<int>(z_of_fr_.size());
  const int num_z = static_cast<int>(layers_.size());
  const int pin_z = z_of_fr_[pin_layer];
  frLayerNum upper = -1;
  if (pin_z >= 0) {
    if (pin_z + 1 < num_z) {
      upper = layers_[pin_z + 1].fr;
    }
  } else if (pin_layer + 2 < num_fr && z_of_fr_[pin_layer + 2] >= 0) {
    upper = pin_layer + 2;
  }
  const frViaDef* via = upper >= 0 ? viaDefAbove(pin_layer) : nullptr;
  if (via != nullptr && via->getLayer2Num() != upper) {
    via = nullptr;
  }
  // Search layers: the pin layer and the next routing layer (Sec. 4.3:
  // access paths "connecting it to on-track points within a small radius").
  std::vector<frLayerNum> frs{pin_layer};
  if (via != nullptr) {
    frs.push_back(upper);
  }
  odb::Rect pin_box;
  pin_box.mergeInit();
  for (const odb::Rect& r : pin) {
    pin_box.merge(r);
  }
  struct Target
  {
    int s;  // search layer
    int z;
    size_t vertex;
    odb::Point p;
  };
  std::vector<Target> targets;
  frCoord max_pitch = 0;
  br::GeoCoord max_tau = 1;
  for (size_t s = 0; s < frs.size(); s++) {
    const int z = z_of_fr_[frs[s]];
    if (z < 0) {
      continue;
    }
    max_pitch = std::max(max_pitch, layers_[z].pitch);
    max_tau = std::max(max_tau, tau_[z]);
    const br::TrackLayer& l = graph_.layers[z];
    const frCoord r = params_.offtrack_radius * layers_[z].pitch;
    const frCoord clo = (l.horizontal ? pin_box.yMin() : pin_box.xMin()) - r;
    const frCoord chi = (l.horizontal ? pin_box.yMax() : pin_box.xMax()) + r;
    const frCoord alo = (l.horizontal ? pin_box.xMin() : pin_box.yMin()) - r;
    const frCoord ahi = (l.horizontal ? pin_box.xMax() : pin_box.yMax()) + r;
    const int i0 = static_cast<int>(std::ranges::lower_bound(l.tracks, clo)
                                    - l.tracks.begin());
    const int j0 = static_cast<int>(std::ranges::lower_bound(l.positions, alo)
                                    - l.positions.begin());
    for (int i = i0; i < l.numTracks() && l.tracks[i] <= chi; i++) {
      for (int j = j0; j < l.numPositions() && l.positions[j] <= ahi; j++) {
        const size_t v = graph_.vertex(z, i, j);
        if (!grid_.usable(v, br::kWire, net)) {
          continue;
        }
        const auto [x, y] = l.point(i, j);
        targets.push_back({static_cast<int>(s), z, v, odb::Point(x, y)});
      }
    }
  }
  if (targets.empty()) {
    stats.empty++;
    return;
  }
  odb::Rect area = pin_box;
  for (const Target& t : targets) {
    area.merge(odb::Rect(t.p, t.p));
  }
  const frCoord margin = 2 * max_tau + max_pitch;
  const br::Box window = toBox(area).bloated(margin);
  const Legality legality(*this, worker);
  br::BlockageGridInput grid_in;
  br::OffTrackRequest req;
  for (const frLayerNum fr : frs) {
    const int z = z_of_fr_[fr];
    br::OffTrackLayer layer;
    layer.wires = z >= 0;
    if (z >= 0) {
      const frCoord half = layers_[z].width / 2;
      legality.blockages(
          fr, window, net, {-half, -half, half, half}, layer.obstacles);
    }
    grid_in.layers.push_back({z >= 0 ? tau_[z] : 1, layer.obstacles});
    req.layers.push_back(std::move(layer));
  }
  std::vector<br::Box> via_blocks;
  std::array<RShape, 3> via_model{};
  if (via != nullptr) {
    const frVia v0(via, odb::Point(0, 0));
    via_model = {RShape{via->getLayer1Num(), v0.getLayer1BBox()},
                 RShape{via->getCutLayerNum(), v0.getCutBBox()},
                 RShape{via->getLayer2Num(), v0.getLayer2BBox()}};
    for (const RShape& m : via_model) {
      legality.blockages(m.layer, window, net, toBox(m.rect), via_blocks);
    }
    grid_in.via_blockages.push_back(via_blocks);
  }
  for (const Target& t : targets) {
    grid_in.points.push_back(toGrid(t.p));
  }
  for (const odb::Rect& r : pin) {
    grid_in.borders.push_back(toBox(r));
  }
  grid_in.borders.push_back(window);
  const std::vector<br::BlockageGridLayer> grid
      = br::buildBlockageGrid(grid_in);
  for (size_t s = 0; s < frs.size(); s++) {
    req.layers[s].grid = grid[s];
  }
  req.window = window;
  for (const odb::Rect& r : pin) {
    req.sources.push_back({0, toBox(r)});
  }
  for (size_t k = 0; k < targets.size(); k++) {
    req.targets.push_back(
        {targets[k].s, toGrid(targets[k].p), static_cast<int>(k)});
  }
  // The search ranks paths by length (Theorem 3.2); a via counts as one
  // wire width, a tie-breaker. The option's cost below uses the on-track
  // cost model.
  req.via_cost = layers_[z_of_fr_[frs.back()]].width;
  // The nearest on-track points; some may fail the exact checks below.
  req.max_targets = 2 * params_.max_access_options;
  req.via_allowed = [&](int, const br::GridPoint& p) {
    for (const br::Box& b : via_blocks) {
      if (b.containsOpen(p.x, p.y)) {
        return false;
      }
    }
    for (const RShape& m : via_model) {
      odb::Rect r = m.rect;
      r.moveDelta(p.x, p.y);
      if (!legality.allowed(m.layer, toBox(r), net)) {
        return false;
      }
    }
    return true;
  };
  utl::Timer timer;
  // One search per thread; its buffers are reused across pins.
  thread_local br::OffTrackSearch search;
  const br::OffTrackResult res = search.run(req);
  stats.search_seconds += timer.elapsed();
  stats.searches++;
  stats.grid_vertices += res.grid_vertices;
  stats.settled += res.settled;
  for (const br::OffTrackPath& path : res.paths) {
    const Target& target = targets[path.target_id];
    AccessOption opt;
    opt.offtrack = true;
    opt.vertex = target.vertex;
    opt.z = target.z;
    for (const br::RouteNode& n : path.nodes) {
      opt.route.push_back({frs[n.layer], n.p});
    }
    Cost cost = 0;
    for (size_t k = 0; k + 1 < opt.route.size(); k++) {
      const br::RouteNode& a = opt.route[k];
      const br::RouteNode& b = opt.route[k + 1];
      if (a.layer != b.layer) {
        const frLayerNum lower = std::min(a.layer, b.layer);
        opt.vias.push_back({viaDefAbove(lower), toPoint(a.p)});
        const int lz = z_of_fr_[lower];
        cost += lz >= 0 ? graph_.via_cost[lz]
                        : static_cast<Cost>(params_.via_cost) * layers_[0].pitch
                              * unit();
        continue;
      }
      if (a.p == b.p) {
        continue;
      }
      const int z = z_of_fr_[a.layer];
      opt.wires.push_back({a.layer, toPoint(a.p), toPoint(b.p)});
      const Cost len = std::abs(a.p.x - b.p.x) + std::abs(a.p.y - b.p.y);
      const bool preferred = (a.p.y == b.p.y) == layers_[z].horizontal;
      cost += len * (preferred ? 1 : params_.jog_factor) * unit();
    }
    opt.cost = cost;
    if (opt.route.size() == 1 || opt.route[1].layer == opt.route[0].layer) {
      opt.planar = true;
      opt.ap = toPoint(opt.route[0].p);
      opt.ap_layer = opt.route[0].layer;
    }
    // Catalogue paths are DRC clean (Sec. 4.3): every shape against the
    // routing space, and no notch within a layer.
    if (!optionLegal(opt, net, worker)) {
      stats.drc_rejected++;
      continue;
    }
    bool notch = false;
    for (size_t k = 0; k < opt.route.size() && !notch;) {
      size_t e = k;
      std::vector<br::GridPoint> run{opt.route[k].p};
      while (e + 1 < opt.route.size()
             && opt.route[e + 1].layer == opt.route[k].layer) {
        e++;
        run.push_back(opt.route[e].p);
      }
      const int z = z_of_fr_[opt.route[k].layer];
      if (z >= 0 && run.size() >= 4) {
        notch
            = br::pathViolations(run, tau_[z], same_net_rules_[z]).notches > 0;
      }
      k = e + 1;
    }
    if (notch) {
      stats.same_net_rejected++;
      continue;
    }
    out.push_back(std::move(opt));
    if (static_cast<int>(out.size()) >= params_.max_access_options) {
      break;
    }
  }
  stats.paths += static_cast<int64_t>(out.size());
  if (out.empty()) {
    stats.empty++;
  }
}

void BonnRouteDR::candidateOf(const AccessOption& option,
                              const int net,
                              br::AccessCandidate& candidate) const
{
  candidate = br::AccessCandidate();
  candidate.length = option.cost;
  candidate.end_layer = option.z;
  const auto [ex, ey] = graph_.point(option.vertex);
  candidate.end = {ex, ey};
  // Tracks the path blocks for wires of other nets.
  std::vector<RShape> shapes;
  accessShapes(option, shapes);
  for (const RShape& s : shapes) {
    const int z = s.layer >= 0 && s.layer < static_cast<int>(z_of_fr_.size())
                      ? z_of_fr_[s.layer]
                      : -1;
    if (z < 0) {
      continue;
    }
    const br::TrackLayer& l = graph_.layers[z];
    const frCoord w = layers_[z].width;
    const frCoord sp = spacing(z, s.rect.minDXDY(), w, s.rect.maxDXDY());
    const frCoord lo
        = (l.horizontal ? s.rect.yMin() : s.rect.xMin()) - sp - w / 2;
    const frCoord hi
        = (l.horizontal ? s.rect.yMax() : s.rect.xMax()) + sp + w / 2;
    auto it = std::ranges::upper_bound(l.tracks, lo);
    for (; it != l.tracks.end() && *it < hi; ++it) {
      candidate.blocked_tracks.push_back(
          (static_cast<int64_t>(z) << 32)
          | static_cast<int64_t>(it - l.tracks.begin()));
    }
  }
  std::ranges::sort(candidate.blocked_tracks);
  candidate.blocked_tracks.erase(
      std::ranges::unique(candidate.blocked_tracks).begin(),
      candidate.blocked_tracks.end());
  // On-track continuation directions at the endpoint: along the track both
  // ways and the via up.
  int z, i, j;
  graph_.decode(option.vertex, z, i, j);
  const br::TrackLayer& l = graph_.layers[z];
  int feasible = 0;
  const size_t v = option.vertex;
  if (j + 1 < l.numPositions() && grid_.usable(v, br::kWireEdge, net)
      && grid_.usable(graph_.vertex(z, i, j + 1), br::kWire, net)) {
    feasible++;
  }
  if (j > 0) {
    const size_t p = graph_.vertex(z, i, j - 1);
    if (grid_.usable(p, br::kWireEdge, net)
        && grid_.usable(p, br::kWire, net)) {
      feasible++;
    }
  }
  if (z + 1 < graph_.numLayers() && l.up_track[j] >= 0 && l.up_pos[i] >= 0) {
    const size_t u = graph_.vertex(z + 1, l.up_track[j], l.up_pos[i]);
    if (grid_.usable(v, br::kViaUpPad, net) && grid_.usable(v, br::kCutUp, net)
        && grid_.usable(u, br::kViaDownPad, net)) {
      feasible++;
    }
  }
  candidate.missing_directions = 3 - feasible;
}

std::vector<int> BonnRouteDR::choosePrimaries(
    const std::vector<const std::vector<AccessOption>*>& catalogues,
    const std::vector<int>& nets,
    AccessStats& stats) const
{
  const int n = static_cast<int>(catalogues.size());
  br::SelectionInput in;
  in.pins.resize(n);
  std::vector<std::vector<std::vector<RShape>>> shapes(n);
  for (int p = 0; p < n; p++) {
    int group = p;
    for (int q = 0; q < p; q++) {
      if (nets[q] == nets[p]) {
        group = in.group[q];
        break;
      }
    }
    in.group.push_back(group);
    for (const AccessOption& opt : *catalogues[p]) {
      br::AccessCandidate c;
      candidateOf(opt, nets[p], c);
      in.pins[p].push_back(std::move(c));
      shapes[p].emplace_back();
      accessShapes(opt, shapes[p].back());
    }
  }
  const Legality legality(*this, nullptr);
  in.conflict = [&](int pa, int ca, int pb, int cb) {
    for (const RShape& a : shapes[pa][ca]) {
      for (const RShape& b : shapes[pb][cb]) {
        if (a.layer == b.layer
            && legality.conflict(a.layer, toBox(a.rect), toBox(b.rect))) {
          return true;
        }
      }
    }
    return false;
  };
  const Cost pitch_unit = static_cast<Cost>(layers_[0].pitch) * unit();
  in.weights.track = params_.access_track_weight * pitch_unit;
  in.weights.direction = params_.access_direction_weight * pitch_unit;
  in.weights.spread
      = static_cast<int64_t>(params_.access_spread_weight) * unit();
  in.weights.spread_distance = 2 * layers_[0].pitch;
  in.node_limit = params_.access_node_limit;
  const br::SelectionResult r = br::selectConflictFree(in);
  stats.selections++;
  stats.selection_nodes += r.nodes;
  stats.selection_optimal += r.optimal;
  return r.choice;
}

void BonnRouteDR::buildOffTrackAccess(const std::vector<char>& pre_skip)
{
  utl::Timer timer;
  frBlock* block = design_->getTopBlock();
  frRegionQuery* rq = design_->getRegionQuery();
  AccessStats& stats = access_stats_;
  // Pins of the terms (same order as buildAccess).
  std::unordered_map<const frInstTerm*, std::pair<int, int>> term_of;
  std::vector<std::pair<int, int>> bterm_terms;
  for (NetData& net : nets_) {
    if (pre_skip[net.id]) {
      continue;
    }
    size_t t = 0;
    for (frInstTerm* it : net.net->getInstTerms()) {
      net.terms[t].inst_term = it;
      term_of[it] = {net.id, static_cast<int>(t)};
      t++;
    }
    for (frBTerm* bt : net.net->getBTerms()) {
      net.terms[t].bterm = bt;
      bterm_terms.emplace_back(net.id, static_cast<int>(t));
      t++;
    }
  }
  // Instances with access pins and their circuit signatures (Sec. 4.3:
  // "Expanding the bounding area of an instance of a circuit and taking
  // all shapes influencing this area and the track coordinates into
  // account, we can identify geometrically equal situations").
  struct InstanceAccess
  {
    frInst* inst = nullptr;
    std::vector<std::pair<int, int>> terms;  // (net, term) per access pin
    std::vector<int> pin_index;
    br::Isometry place;
    std::vector<int64_t> signature;
    int cls = -1;
  };
  std::vector<InstanceAccess> instances;
  for (const auto& inst : block->getInsts()) {
    InstanceAccess ia;
    ia.inst = inst.get();
    int k = 0;
    for (const auto& it : inst->getInstTerms()) {
      auto f = term_of.find(it.get());
      if (f != term_of.end()) {
        ia.terms.push_back(f->second);
        ia.pin_index.push_back(k);
      }
      k++;
    }
    if (!ia.terms.empty()) {
      ia.place = isometryOf(inst->getTransform());
      instances.push_back(std::move(ia));
    }
  }
  const int num_z = static_cast<int>(layers_.size());
  frCoord max_pitch = 0;
  br::GeoCoord max_tau = 1;
  for (int z = 0; z < std::min(num_z, 3); z++) {
    max_pitch = std::max(max_pitch, layers_[z].pitch);
    max_tau = std::max(max_tau, tau_[z]);
  }
  // Reach of a pin's search: targets within the radius, the window margin
  // (2 tau + pitch) and the spacing plus model extent of the obstacles that
  // can enter the window (at most two pitches).
  const frCoord halo
      = params_.offtrack_radius * max_pitch + 2 * max_tau + 3 * max_pitch;
  const frLayerNum bottom = layers_.front().fr;
  const frLayerNum sig_top = layers_[std::min(num_z, 3) - 1].fr;
  const int threads = static_cast<int>(workers_.size());
#pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (int k = 0; k < static_cast<int>(instances.size()); k++) {
    InstanceAccess& ia = instances[k];
    br::CircuitSignature sig(ia.place.inverse());
    odb::Rect window = ia.inst->getBBox();
    window.bloat(halo, window);
    std::vector<int> groups;
    for (size_t p = 0; p < ia.terms.size(); p++) {
      int g = static_cast<int>(p) + 1;
      for (size_t q = 0; q < p; q++) {
        if (ia.terms[q].first == ia.terms[p].first) {
          g = groups[q];
          break;
        }
      }
      groups.push_back(g);
      sig.addPin(ia.pin_index[p], g);
    }
    for (frLayerNum fr = std::max(0, bottom - 2); fr <= sig_top; fr++) {
      frRegionQuery::Objects<frBlockObject> result;
      rq->query(window, fr, result);
      for (const auto& [rect, obj] : result) {
        frNet* owner_net = nullptr;
        if (obj->typeId() == frcInstTerm) {
          owner_net = static_cast<frInstTerm*>(obj)->getNet();
        } else if (obj->typeId() == frcBTerm) {
          owner_net = static_cast<frBTerm*>(obj)->getNet();
        }
        int label = 0;
        if (owner_net != nullptr) {
          auto it = net_index_.find(owner_net);
          for (size_t p = 0; it != net_index_.end() && p < ia.terms.size();
               p++) {
            if (ia.terms[p].first == it->second) {
              label = groups[p];
              break;
            }
          }
        }
        const odb::Rect clipped = rect.intersect(window);
        sig.addShape(fr, toBox(clipped), label);
      }
      const frViaDef* via = viaDefAbove(fr);
      if (via != nullptr) {
        const frVia v0(via, odb::Point(0, 0));
        sig.addVia(fr,
                   {{via->getLayer1Num(), toBox(v0.getLayer1BBox())},
                    {via->getCutLayerNum(), toBox(v0.getCutBBox())},
                    {via->getLayer2Num(), toBox(v0.getLayer2BBox())}});
      }
    }
    for (int z = 0; z < std::min(num_z, 4); z++) {
      const br::TrackLayer& l = graph_.layers[z];
      const frCoord lo = l.horizontal ? window.yMin() : window.xMin();
      const frCoord hi = l.horizontal ? window.yMax() : window.xMax();
      for (auto it = std::ranges::lower_bound(l.tracks, lo);
           it != l.tracks.end() && *it <= hi;
           ++it) {
        sig.addTrack(z, l.horizontal, *it);
      }
    }
    sig.addValue(halo);
    sig.addValue(params_.offtrack_radius);
    if (!params_.access_classes) {
      sig.addValue(k);
    }
    ia.signature = sig.finish();
  }
  br::CircuitClassTable table;
  std::vector<int> representative;
  for (int k = 0; k < static_cast<int>(instances.size()); k++) {
    const auto [cls, is_new] = table.classify(instances[k].signature);
    instances[k].cls = cls;
    if (is_new) {
      representative.push_back(k);
    }
    instances[k].signature.clear();
    instances[k].signature.shrink_to_fit();
  }
  stats.instances = static_cast<int64_t>(instances.size());
  stats.classes = table.numClasses();
  // Catalogues and conflict-free solution of each class representative,
  // against fixed shapes and pins only (preprocessing).
  struct ClassData
  {
    std::vector<std::vector<AccessOption>> catalogues;
    std::vector<int> choice;
  };
  std::vector<ClassData> classes(representative.size());
  std::vector<AccessStats> class_stats(representative.size());
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
  for (int c = 0; c < static_cast<int>(representative.size()); c++) {
    const InstanceAccess& ia = instances[representative[c]];
    ClassData& cd = classes[c];
    cd.catalogues.resize(ia.terms.size());
    std::vector<const std::vector<AccessOption>*> cats;
    std::vector<int> nets;
    for (size_t p = 0; p < ia.terms.size(); p++) {
      const auto [net, t] = ia.terms[p];
      offTrackCatalogue(
          net, nets_[net].terms[t], nullptr, cd.catalogues[p], class_stats[c]);
      cats.push_back(&cd.catalogues[p]);
      nets.push_back(net);
    }
    cd.choice = choosePrimaries(cats, nets, class_stats[c]);
  }
  for (const AccessStats& s : class_stats) {
    stats.add(s);
  }
  // Per instance: map the class solution, check it against the routing
  // space (reservations of earlier instances), re-solve if needed, and
  // reserve the primary paths.
  auto map_option =
      [&](const AccessOption& in, const br::Isometry& m, AccessOption& out) {
        out = in;
        for (br::RouteNode& n : out.route) {
          n.p = m.apply(n.p);
        }
        for (OutWire& w : out.wires) {
          w.begin = toPoint(m.apply(toGrid(w.begin)));
          w.end = toPoint(m.apply(toGrid(w.end)));
        }
        for (OutVia& v : out.vias) {
          v.point = toPoint(m.apply(toGrid(v.point)));
        }
        out.ap = toPoint(m.apply(toGrid(in.ap)));
        const auto [x, y] = graph_.point(in.vertex);
        const br::GridPoint e = m.apply(br::GridPoint{x, y});
        const br::TrackLayer& l = graph_.layers[in.z];
        const frCoord cross = l.horizontal ? e.y : e.x;
        const frCoord along = l.horizontal ? e.x : e.y;
        auto ti = std::ranges::lower_bound(l.tracks, cross);
        auto pj = std::ranges::lower_bound(l.positions, along);
        if (ti == l.tracks.end() || *ti != cross || pj == l.positions.end()
            || *pj != along) {
          return false;
        }
        out.vertex = graph_.vertex(in.z,
                                   static_cast<int>(ti - l.tracks.begin()),
                                   static_cast<int>(pj - l.positions.begin()));
        return true;
      };
  auto reserve = [&](NetData& net, Term& term, const int primary) {
    term.primary = primary;
    accessShapes(term.options[primary], term.reservation);
    setReservation(term, net.id, true);
    stats.primaries++;
  };
  auto reserve_fallback = [&](NetData& net, Term& term) {
    // The previous access stubs: the first conflict-free one.
    stats.fallback++;
    for (size_t o = 0; o < term.options.size(); o++) {
      std::vector<RShape> shapes;
      accessShapes(term.options[o], shapes);
      if (std::ranges::any_of(shapes, [&](const RShape& shape) {
            return shapeConflicts(shape, net.id);
          })) {
        continue;
      }
      reserve(net, term, static_cast<int>(o));
      break;
    }
  };
  // Solves the pins against the current routing space: options that are
  // illegal now are kept in the catalogue (they may become feasible) but
  // cannot be primary.
  auto solve = [&](std::vector<std::vector<AccessOption>>& cats,
                   const std::vector<int>& nets) {
    const int n = static_cast<int>(cats.size());
    std::vector<std::vector<AccessOption>> legal(n);
    std::vector<std::vector<int>> index(n);
    std::vector<const std::vector<AccessOption>*> ptrs;
    for (int p = 0; p < n; p++) {
      for (size_t o = 0; o < cats[p].size(); o++) {
        if (optionLegal(cats[p][o], nets[p], nullptr)) {
          legal[p].push_back(cats[p][o]);
          index[p].push_back(static_cast<int>(o));
        }
      }
      ptrs.push_back(&legal[p]);
    }
    std::vector<int> choice = choosePrimaries(ptrs, nets, stats);
    for (int p = 0; p < n; p++) {
      choice[p] = choice[p] >= 0 ? index[p][choice[p]] : -1;
    }
    return choice;
  };
  for (int k = 0; k < static_cast<int>(instances.size()); k++) {
    const InstanceAccess& ia = instances[k];
    const int rep = representative[ia.cls];
    const ClassData& cd = classes[ia.cls];
    const int n = static_cast<int>(ia.terms.size());
    const br::Isometry m = ia.place.compose(instances[rep].place.inverse());
    std::vector<std::vector<AccessOption>> cats(n);
    std::vector<int> nets;
    bool mapped = true;
    for (int p = 0; p < n; p++) {
      const auto [net, t] = ia.terms[p];
      nets.push_back(net);
      if (k == rep) {
        cats[p] = cd.catalogues[p];
        continue;
      }
      for (const AccessOption& opt : cd.catalogues[p]) {
        AccessOption mo;
        if (!map_option(opt, m, mo)) {
          mapped = false;
          break;
        }
        cats[p].push_back(std::move(mo));
      }
      if (!mapped) {
        break;
      }
    }
    if (!mapped) {
      // Equal signatures with unequal track positions outside the window;
      // generate this instance's catalogues directly.
      stats.mapping_failed++;
      for (int p = 0; p < n; p++) {
        const auto [net, t] = ia.terms[p];
        offTrackCatalogue(net, nets_[net].terms[t], nullptr, cats[p], stats);
      }
    }
    std::vector<int> choice = mapped ? cd.choice : std::vector<int>(n, -1);
    bool reuse = mapped;
    for (int p = 0; p < n && reuse; p++) {
      reuse = choice[p] < 0 ? cats[p].empty()
                            : optionLegal(cats[p][choice[p]], nets[p], nullptr);
    }
    if (reuse) {
      stats.class_reused += k != rep;
    } else {
      stats.reselected++;
      choice = solve(cats, nets);
    }
    for (int p = 0; p < n; p++) {
      const auto [net_id, t] = ia.terms[p];
      NetData& net = nets_[net_id];
      Term& term = net.terms[t];
      if (cats[p].empty()) {
        reserve_fallback(net, term);
        continue;
      }
      term.options = std::move(cats[p]);
      if (choice[p] >= 0) {
        reserve(net, term, choice[p]);
      }
    }
  }
  // Top-level pins: catalogue and best feasible path each.
  for (const auto& [net_id, t] : bterm_terms) {
    NetData& net = nets_[net_id];
    Term& term = net.terms[t];
    std::vector<std::vector<AccessOption>> cats(1);
    offTrackCatalogue(net_id, term, nullptr, cats[0], stats);
    if (cats[0].empty()) {
      reserve_fallback(net, term);
      continue;
    }
    const std::vector<int> choice = solve(cats, {net_id});
    term.options = std::move(cats[0]);
    if (choice[0] >= 0) {
      reserve(net, term, choice[0]);
    }
  }
  for (NetData& net : nets_) {
    if (pre_skip[net.id]) {
      continue;
    }
    net.skip = std::ranges::any_of(
        net.terms, [](const Term& term) { return term.options.empty(); });
  }
  stats.seconds = timer.elapsed();
  if (params_.verbose) {
    logger_->report(
        "BonnRoute off-track access: {} instances in {} circuit classes, {} "
        "reused class solutions, {} re-solved, {} mapping failures; {} pin "
        "catalogues, {} paths ({:.2f} per pin), {} empty, {} fallback stubs; "
        "{} DRC-rejected, {} same-net-rejected; {} blockage-grid vertices, "
        "{} settled; search {:.2f}s; {} selections ({} optimal, {} nodes); "
        "{} primary paths; {:.2f}s.",
        stats.instances,
        stats.classes,
        stats.class_reused,
        stats.reselected,
        stats.mapping_failed,
        stats.searches,
        stats.paths,
        stats.searches > 0 ? double(stats.paths) / stats.searches : 0.0,
        stats.empty,
        stats.fallback,
        stats.drc_rejected,
        stats.same_net_rejected,
        stats.grid_vertices,
        stats.settled,
        stats.search_seconds,
        stats.selections,
        stats.selection_optimal,
        stats.selection_nodes,
        stats.primaries,
        stats.seconds);
  }
}

bool BonnRouteDR::regenerateAccess(NetData& net, const int t, Worker& worker)
{
  Term& term = net.terms[t];
  if (term.inst_term == nullptr && term.bterm == nullptr) {
    return false;
  }
  AccessStats stats;
  std::vector<AccessOption> fresh;
  offTrackCatalogue(net.id, term, &worker, fresh, stats);
  worker.regenerations++;
  int added = 0;
  for (AccessOption& opt : fresh) {
    if (term.options.size() >= 1000) {
      break;
    }
    const bool known
        = std::ranges::any_of(term.options, [&](const AccessOption& o) {
            return o.vertex == opt.vertex && o.route == opt.route;
          });
    if (!known) {
      term.options.push_back(std::move(opt));
      added++;
    }
  }
  worker.regenerated_paths += added;
  return added > 0;
}

void BonnRouteDR::routeToGeometry(const std::vector<br::RouteNode>& nodes,
                                  const std::vector<OutVia>& hints,
                                  NetData& net) const
{
  for (size_t k = 0; k + 1 < nodes.size(); k++) {
    const br::RouteNode& a = nodes[k];
    const br::RouteNode& b = nodes[k + 1];
    if (a.layer != b.layer) {
      const frLayerNum lower = std::min(a.layer, b.layer);
      const odb::Point p = toPoint(a.p);
      const frViaDef* def = nullptr;
      for (const OutVia& h : hints) {
        if (h.point == p && h.def != nullptr
            && h.def->getLayer1Num() == lower) {
          def = h.def;
          break;
        }
      }
      if (def == nullptr) {
        def = viaDefAbove(lower);
      }
      if (def != nullptr) {
        net.vias.push_back({def, p});
      }
      continue;
    }
    if (!(a.p == b.p)) {
      net.wires.push_back({a.layer, toPoint(a.p), toPoint(b.p)});
    }
  }
}

void BonnRouteDR::postprocessConnection(std::vector<br::RouteNode>& nodes,
                                        const std::vector<OutVia>& hints,
                                        const int net,
                                        Worker& worker,
                                        std::vector<OutWire>& patches) const
{
  const Legality legality(*this, &worker);
  br::CleanupRequest req;
  req.nodes = nodes;
  req.layers = cleanup_layers_;
  req.allowed = [&](const int layer, const br::Box& box) {
    return legality.allowed(layer, box, net);
  };
  req.via_pads = [&](const int layer,
                     const br::GridPoint& p,
                     const bool up,
                     std::vector<br::Box>& out) {
    const frLayerNum lower = up ? layer : layer - 2;
    const odb::Point point = toPoint(p);
    const frViaDef* def = nullptr;
    for (const OutVia& h : hints) {
      if (h.point == point && h.def != nullptr
          && h.def->getLayer1Num() == lower) {
        def = h.def;
        break;
      }
    }
    if (def == nullptr) {
      def = viaDefAbove(lower);
    }
    if (def != nullptr) {
      const frVia via(def, point);
      out.push_back(toBox(up ? via.getLayer1BBox() : via.getLayer2BBox()));
    }
  };
  const br::CleanupResult r = br::cleanupConnection(req);
  nodes = r.nodes;
  worker.cleanup.add(r.stats);
  for (const br::CleanupPatch& patch : r.patches) {
    patches.push_back({static_cast<frLayerNum>(patch.layer),
                       toPoint(patch.begin),
                       toPoint(patch.end)});
  }
}

void BonnRouteDR::accessShapes(const AccessOption& option,
                               std::vector<RShape>& shapes) const
{
  shapes.clear();
  for (const OutWire& wire : option.wires) {
    shapes.push_back(wireShape(wire));
  }
  for (const OutVia& via : option.vias) {
    viaShapes(via, shapes);
  }
}

void BonnRouteDR::setReservation(Term& term,
                                 const int net,
                                 const bool active,
                                 Worker* local)
{
  if (term.reservation_active == active) {
    return;
  }
  for (const RShape& shape : term.reservation) {
    if (shape.layer < 0
        || shape.layer >= static_cast<int>(reservation_tree_.size())) {
      continue;
    }
    markShape(shape, 3, net, active ? 1 : -1);
    auto& tree = reservation_tree_[shape.layer];
    if (local != nullptr) {
      local->released_shapes.emplace_back(shape, net);
    } else if (active) {
      tree.insert({shape.rect, net});
    } else {
      tree.remove(std::make_pair(shape.rect, net));
    }
  }
  term.reservation_active = active;
}

void BonnRouteDR::setReservations(NetData& net,
                                  const bool active,
                                  Worker* local)
{
  for (Term& term : net.terms) {
    setReservation(term, net.id, active, local);
  }
  if (!active) {
    setInputReservation(net, false, local);
  }
  if (local != nullptr) {
    local->released.insert(net.id);
  }
}

void BonnRouteDR::setInputReservation(NetData& net,
                                      const bool active,
                                      Worker* local)
{
  if (net.input_reserved == active) {
    return;
  }
  for (const RShape& shape : net.input_reservation) {
    markShape(shape, 3, net.id, active ? 1 : -1);
    auto& tree = reservation_tree_[shape.layer];
    if (local != nullptr) {
      local->released_shapes.emplace_back(shape, net.id);
    } else if (active) {
      tree.insert({shape.rect, net.id});
    } else {
      tree.remove(std::make_pair(shape.rect, net.id));
    }
  }
  net.input_reserved = active;
}

void BonnRouteDR::reserveInputWires()
{
  // Input wiring: the track assignment's path segments in the guides.
  const int num_fr = static_cast<int>(z_of_fr_.size());
  struct InputWire
  {
    int net;
    frLayerNum layer;
    odb::Point begin;
    odb::Point end;
  };
  std::vector<InputWire> wires;
  std::vector<std::vector<std::pair<odb::Rect, int>>> per_layer(num_fr);
  for (NetData& net : nets_) {
    if (net.skip) {
      continue;
    }
    for (const auto& guide : net.net->getGuides()) {
      for (const auto& route : guide->getRoutes()) {
        if (route->typeId() != frcPathSeg) {
          continue;
        }
        const auto* seg = static_cast<const frPathSeg*>(route.get());
        const frLayerNum fr = seg->getLayerNum();
        const auto [begin, end] = seg->getPoints();
        if (fr < 0 || fr >= num_fr || z_of_fr_[fr] < 0 || begin == end) {
          continue;
        }
        wires.push_back({net.id, fr, begin, end});
        per_layer[fr].emplace_back(wireShape({fr, begin, end}).rect, net.id);
      }
    }
  }
  std::vector<RTree<int>> input_tree(num_fr);
  for (int fr = 0; fr < num_fr; fr++) {
    input_tree[fr] = RTree<int>(per_layer[fr].begin(), per_layer[fr].end());
  }
  input_wires_ = static_cast<int64_t>(wires.size());
  // Reserve the parts (here: pieces between global routing tile
  // boundaries) of preferred-direction input wires on tracks that conflict
  // neither with fixed shapes, pin access reservations nor input wires of
  // other nets (Sec. 5.2).
  for (const InputWire& w : wires) {
    NetData& net = nets_[w.net];
    const int z = z_of_fr_[w.layer];
    if (z < params_.reserve_min_layer) {
      continue;
    }
    const br::TrackLayer& l = graph_.layers[z];
    const bool horizontal = w.begin.y() == w.end.y();
    if (horizontal != l.horizontal) {
      continue;
    }
    const frCoord cross = horizontal ? w.begin.y() : w.begin.x();
    auto track = std::ranges::lower_bound(l.tracks, cross);
    if (track == l.tracks.end() || *track != cross) {
      continue;
    }
    const int i = static_cast<int>(track - l.tracks.begin());
    const frCoord lo = horizontal ? std::min(w.begin.x(), w.end.x())
                                  : std::min(w.begin.y(), w.end.y());
    const frCoord hi = horizontal ? std::max(w.begin.x(), w.end.x())
                                  : std::max(w.begin.y(), w.end.y());
    const frCoord origin = horizontal ? tile_x0_ : tile_y0_;
    const frCoord pitch = horizontal ? tile_dx_ : tile_dy_;
    auto point = [&](frCoord along) {
      return horizontal ? odb::Point(along, cross) : odb::Point(cross, along);
    };
    auto reserve = [&](frCoord a, frCoord b) {
      if (b <= a) {
        return;
      }
      net.input_reservation.push_back(wireShape({w.layer, point(a), point(b)}));
      const int jlo = static_cast<int>(std::ranges::lower_bound(l.positions, a)
                                       - l.positions.begin());
      const int jhi = static_cast<int>(std::ranges::upper_bound(l.positions, b)
                                       - l.positions.begin())
                      - 1;
      if (jhi > jlo) {
        net.discount_runs.push_back({z, i, jlo, jhi});
      }
      input_reserved_++;
    };
    frCoord run_lo = lo;
    bool run_open = false;
    for (frCoord a = lo; a < hi;) {
      const frCoord next_boundary = origin + ((a - origin) / pitch + 1) * pitch;
      const frCoord b = std::min(hi, next_boundary);
      const RShape piece = wireShape({w.layer, point(a), point(b)});
      const bool legal = !shapeConflicts(piece, net.id)
                         && !treeConflicts(input_tree[w.layer], piece, net.id);
      if (legal && !run_open) {
        run_lo = a;
        run_open = true;
      } else if (!legal && run_open) {
        reserve(run_lo, a);
        run_open = false;
      }
      a = b;
    }
    if (run_open) {
      reserve(run_lo, hi);
    }
  }
  for (NetData& net : nets_) {
    if (!net.skip) {
      setInputReservation(net, true);
    }
  }
}

void BonnRouteDR::buildGuideTiles(NetData& net)
{
  net.guide_tiles.clear();
  const int num_z = graph_.numLayers();
  frBlock* block = design_->getTopBlock();
  for (const auto& guide : net.net->getGuides()) {
    const frLayerNum fr = guide->getLayerNum();
    int z = -1;
    if (fr >= 0 && fr < static_cast<int>(z_of_fr_.size())) {
      z = z_of_fr_[fr];
    }
    if (z < 0) {
      z = fr < layers_.front().fr ? 0 : num_z - 1;
    }
    const odb::Point a = block->getGCellIdx(guide->getBeginPoint());
    const odb::Point b = block->getGCellIdx(guide->getEndPoint());
    for (int x = std::min(a.x(), b.x()); x <= std::max(a.x(), b.x()); x++) {
      for (int y = std::min(a.y(), b.y()); y <= std::max(a.y(), b.y()); y++) {
        net.guide_tiles.push_back({z, x, y});
      }
    }
  }
  std::ranges::sort(net.guide_tiles);
  net.guide_tiles.erase(std::ranges::unique(net.guide_tiles).begin(),
                        net.guide_tiles.end());
}

void BonnRouteDR::buildRegion(const NetData& net,
                              const int dilation,
                              br::Region& region)
{
  region.reset();
  auto allow = [&](int z, int tx, int ty) {
    for (int dz = -params_.layer_dilation; dz <= params_.layer_dilation; dz++) {
      for (int dx = -dilation; dx <= dilation; dx++) {
        for (int dy = -dilation; dy <= dilation; dy++) {
          region.allow(z + dz, tx + dx, ty + dy);
        }
      }
    }
  };
  for (const auto& [z, x, y] : net.guide_tiles) {
    allow(z, x, y);
  }
  for (const Term& term : net.terms) {
    for (const AccessOption& o : term.options) {
      const auto [x, y] = graph_.point(o.vertex);
      allow(o.z, tileX(x), tileY(y));
    }
  }
}

bool BonnRouteDR::buildPathRegion(const NetData& net,
                                  const std::vector<char>& connected,
                                  const std::vector<size_t>& tree,
                                  const std::vector<std::vector<char>>& valid,
                                  Worker& worker)
{
  int xlo, ylo, xhi, yhi;
  netTileBox(net, xlo, ylo, xhi, yhi);
  xlo = std::max(xlo, 0);
  ylo = std::max(ylo, 0);
  xhi = std::min(xhi, tiles_x_ - 1);
  yhi = std::min(yhi, tiles_y_ - 1);
  if (xlo > xhi || ylo > yhi) {
    return false;
  }
  const int nx = xhi - xlo + 1;
  const int ny = yhi - ylo + 1;
  const int nz = graph_.numLayers();
  const size_t size = static_cast<size_t>(nx) * ny * nz;
  if (worker.tile_guide.size() < size) {
    worker.tile_guide.resize(size, 0);
    worker.tile_seen.resize(size, 0);
    worker.tile_target.resize(size, 0);
    worker.tile_parent.resize(size, -1);
  }
  if (++worker.tile_epoch == 0) {
    std::ranges::fill(worker.tile_guide, 0);
    std::ranges::fill(worker.tile_seen, 0);
    std::ranges::fill(worker.tile_target, 0);
    worker.tile_epoch = 1;
  }
  const uint32_t epoch = worker.tile_epoch;
  auto index = [&](int z, int x, int y) -> int64_t {
    if (z < 0 || z >= nz || x < xlo || x > xhi || y < ylo || y > yhi) {
      return -1;
    }
    return (static_cast<int64_t>(z) * nx + (x - xlo)) * ny + (y - ylo);
  };
  auto tile_index = [&](size_t v) {
    int z, i, j;
    graph_.decode(v, z, i, j);
    const auto [x, y] = graph_.layers[z].point(i, j);
    return index(z, tileX(x), tileY(y));
  };
  for (const auto& [z, x, y] : net.guide_tiles) {
    const int64_t t = index(z, x, y);
    if (t >= 0) {
      worker.tile_guide[t] = epoch;
    }
  }
  bool any_target = false;
  for (size_t t = 0; t < net.terms.size(); t++) {
    if (connected[t]) {
      continue;
    }
    for (size_t o = 0; o < net.terms[t].options.size(); o++) {
      const int64_t at
          = valid[t][o] ? tile_index(net.terms[t].options[o].vertex) : -1;
      if (at >= 0) {
        worker.tile_target[at] = epoch;
        worker.tile_guide[at] = epoch;
        any_target = true;
      }
    }
  }
  if (!any_target) {
    return false;
  }
  std::vector<int64_t> queue;
  auto push_source = [&](size_t v) {
    const int64_t at = tile_index(v);
    if (at >= 0 && worker.tile_seen[at] != epoch) {
      worker.tile_seen[at] = epoch;
      worker.tile_guide[at] = epoch;
      worker.tile_parent[at] = -1;
      queue.push_back(at);
    }
  };
  if (tree.empty()) {
    for (size_t o = 0; o < net.terms[0].options.size(); o++) {
      if (valid[0][o]) {
        push_source(net.terms[0].options[o].vertex);
      }
    }
  } else {
    for (const size_t v : tree) {
      push_source(v);
    }
  }
  // Breadth-first search over guide tiles: a shortest tile path along the
  // global routing (planar steps and layer changes).
  int64_t found = -1;
  for (size_t head = 0; head < queue.size() && found < 0; head++) {
    const int64_t at = queue[head];
    if (worker.tile_target[at] == epoch) {
      found = at;
      break;
    }
    const int z = static_cast<int>(at / (static_cast<int64_t>(nx) * ny));
    const int x = xlo + static_cast<int>((at / ny) % nx);
    const int y = ylo + static_cast<int>(at % ny);
    for (const auto& [dz, dx, dy] : {std::array<int, 3>{1, 0, 0},
                                     std::array<int, 3>{-1, 0, 0},
                                     std::array<int, 3>{0, 1, 0},
                                     std::array<int, 3>{0, -1, 0},
                                     std::array<int, 3>{0, 0, 1},
                                     std::array<int, 3>{0, 0, -1}}) {
      const int64_t next = index(z + dz, x + dx, y + dy);
      if (next < 0 || worker.tile_guide[next] != epoch
          || worker.tile_seen[next] == epoch) {
        continue;
      }
      worker.tile_seen[next] = epoch;
      worker.tile_parent[next] = at;
      queue.push_back(next);
    }
  }
  if (found < 0) {
    return false;
  }
  br::Region& region = worker.region;
  region.reset();
  int steps = 0;
  int first_x = 0;
  int first_y = 0;
  int last_x = 0;
  int last_y = 0;
  bool first = true;
  for (int64_t at = found; at >= 0; at = worker.tile_parent[at]) {
    const int z = static_cast<int>(at / (static_cast<int64_t>(nx) * ny));
    const int x = xlo + static_cast<int>((at / ny) % nx);
    const int y = ylo + static_cast<int>(at % ny);
    if (first) {
      first_x = x;
      first_y = y;
      first = false;
    } else if (x != last_x || y != last_y) {
      steps++;
    }
    last_x = x;
    last_y = y;
    for (int dz = -params_.layer_dilation; dz <= params_.layer_dilation; dz++) {
      region.allow(z + dz, x, y);
    }
  }
  worker.path_detour
      = steps - std::abs(last_x - first_x) - std::abs(last_y - first_y);
  return true;
}

void BonnRouteDR::pathToGeometry(const std::vector<size_t>& path,
                                 NetData& net) const
{
  // path runs target -> source; geometry is direction agnostic.
  size_t k = 0;
  while (k + 1 < path.size()) {
    int za, ia, ja, zb, ib, jb;
    graph_.decode(path[k], za, ia, ja);
    graph_.decode(path[k + 1], zb, ib, jb);
    if (za != zb) {
      const size_t lower = za < zb ? path[k] : path[k + 1];
      const auto [x, y] = graph_.point(lower);
      net.vias.push_back({layers_[std::min(za, zb)].up_via, odb::Point(x, y)});
      k++;
      continue;
    }
    // Merge a straight run of along (same track) or jog (same position)
    // moves.
    const bool along = ia == ib;
    size_t e = k + 1;
    while (e + 1 < path.size()) {
      int zc, ic, jc;
      graph_.decode(path[e + 1], zc, ic, jc);
      int ze, ie, je;
      graph_.decode(path[e], ze, ie, je);
      if (zc != za || (along ? ic != ia : jc != ja)) {
        break;
      }
      e++;
    }
    const auto [x0, y0] = graph_.point(path[k]);
    const auto [x1, y1] = graph_.point(path[e]);
    net.wires.push_back(
        {layers_[za].fr, odb::Point(x0, y0), odb::Point(x1, y1)});
    k = e;
  }
}

bool BonnRouteDR::connect(NetData& net,
                          const bool ripup,
                          Worker& worker,
                          const int extra_dilation)
{
  const int num_terms = static_cast<int>(net.terms.size());
  // Stub and attachment legality depend on the committed shapes, not on
  // the guide dilation. Check them once for both clean windows.
  std::vector<std::vector<char>> valid(num_terms);
  auto option_valid = [&](const AccessOption& opt) {
    bool ok = grid_.usable(opt.vertex, br::kWire, net.id)
              || (ripup && grid_.ripupable(opt.vertex, br::kWire, net.id));
    if (ok && !ripup && opt.offtrack) {
      // Blockage-grid paths: every shape against earlier routed nets
      // (Sec. 4.4 "these paths are first checked for diff-net rule
      // violations to earlier routed nets").
      return optionLegal(opt, net.id, &worker);
    }
    return ok;
  };
  for (int t = 0; t < num_terms; t++) {
    const Term& term = net.terms[t];
    valid[t].assign(term.options.size(), 0);
    int feasible = 0;
    for (size_t o = 0; o < term.options.size(); o++) {
      const AccessOption& opt = term.options[o];
      if (opt.offtrack) {
        valid[t][o] = option_valid(opt);
        feasible += valid[t][o];
        continue;
      }
      bool ok = grid_.usable(opt.vertex, br::kWire, net.id)
                || (ripup && grid_.ripupable(opt.vertex, br::kWire, net.id));
      if (ok && !ripup) {
        for (const OutWire& w : opt.wires) {
          if (shapeConflicts(wireShape(w), net.id, &worker)) {
            ok = false;
            break;
          }
        }
        // Landing pad against routed shapes and other pins' reservations.
        for (const OutVia& v : opt.vias) {
          if (!ok) {
            break;
          }
          const frVia via(v.def, v.point);
          ok = !shapeConflicts(
              {v.def->getLayer2Num(), via.getLayer2BBox()}, net.id, &worker);
        }
      }
      valid[t][o] = ok;
      feasible += ok;
    }
    // Sec. 4.4: "If no or only a small number of access paths remains
    // feasible, we dynamically generate new access paths."
    if (!ripup && params_.offtrack_access
        && feasible < params_.offtrack_regenerate) {
      const size_t before = net.terms[t].options.size();
      if (regenerateAccess(net, t, worker)) {
        const Term& grown = net.terms[t];
        valid[t].resize(grown.options.size(), 0);
        for (size_t o = before; o < grown.options.size(); o++) {
          valid[t][o] = option_valid(grown.options[o]);
        }
      }
    }
    if (std::ranges::find(valid[t], 1) == valid[t].end()) {
      return false;
    }
  }
  const int retry = params_.retry_dilation + extra_dilation;
  if (ripup) {
    return connectInRegion(net, retry, true, valid, worker);
  }
  return connectInRegion(net, extra_dilation, false, valid, worker)
         || (retry != extra_dilation
             && connectInRegion(net, retry, false, valid, worker));
}

bool BonnRouteDR::connectInRegion(NetData& net,
                                  const int dilation,
                                  const bool ripup,
                                  const std::vector<std::vector<char>>& valid,
                                  Worker& worker)
{
  net.wires.clear();
  net.vias.clear();
  net.pin_points.clear();
  for (Term& term : net.terms) {
    term.chosen = -1;
  }
  const int num_terms = static_cast<int>(net.terms.size());
  br::SearchRequest req;
  req.net = net.id;
  req.region = &worker.region;
  req.ripup = ripup;
  req.ripup_cost
      = static_cast<Cost>(params_.ripup_cost) * layers_[0].pitch * unit();
  req.discounts = net.discount_runs;
  if (params_.spreading && !spread_tiles_.empty()) {
    req.spread_tiles = &spread_tiles_;
    req.spread_cost
        = static_cast<Cost>(params_.spread_cost) * layers_[0].pitch * unit();
  }
  req.discount_rate = static_cast<Cost>(params_.discount_numerator) * unit()
                      / params_.discount_denominator;
  req.use_intervals = params_.use_intervals;
  req.interval_via_property = params_.interval_via_property;
  std::unordered_set<size_t> in_tree;
  std::vector<size_t> tree;
  std::vector<std::pair<size_t, Cost>> sources;
  std::vector<char> connected(num_terms, 0);
  connected[0] = 1;
  // Terms whose access path is part of a post-processed connection.
  std::vector<char> consumed(num_terms, 0);
  for (size_t o = 0; o < net.terms[0].options.size(); o++) {
    if (valid[0][o]) {
      sources.emplace_back(net.terms[0].options[o].vertex,
                           net.terms[0].options[o].cost);
    }
  }
  auto in_region = [&](size_t v) {
    int z, i, j;
    graph_.decode(v, z, i, j);
    return worker.region.contains(graph_.layers[z], z, i, j);
  };
  // Targets inside the current routing area and the future cost pi_H
  // (Sec. 4.1): l1 distance plus via lower bound to the nearest target;
  // per-term boxes keep it cheap for large nets.
  auto set_targets = [&](const bool restrict) {
    req.targets.clear();
    req.target_boxes.clear();
    std::vector<std::pair<int, size_t>> chosen;
    for (int t = 1; t < num_terms; t++) {
      if (connected[t]) {
        continue;
      }
      for (size_t o = 0; o < net.terms[t].options.size(); o++) {
        if (valid[t][o]
            && (!restrict || in_region(net.terms[t].options[o].vertex))) {
          chosen.emplace_back(t, o);
        }
      }
    }
    const bool per_point = chosen.size() <= 64;
    int box_term = -1;
    for (const auto& [t, o] : chosen) {
      const AccessOption& opt = net.terms[t].options[o];
      req.targets.push_back(
          {opt.vertex, opt.cost, t * 1024 + static_cast<int>(o)});
      const auto [x, y] = graph_.point(opt.vertex);
      if (per_point) {
        req.target_boxes.push_back({x, y, x, y, opt.z, opt.z});
        continue;
      }
      if (t != box_term) {
        req.target_boxes.push_back({x, y, x, y, opt.z, opt.z});
        box_term = t;
        continue;
      }
      br::TargetBox& box = req.target_boxes.back();
      box.xlo = std::min(box.xlo, x);
      box.ylo = std::min(box.ylo, y);
      box.xhi = std::max(box.xhi, x);
      box.yhi = std::max(box.yhi, y);
      box.zlo = std::min(box.zlo, opt.z);
      box.zhi = std::max(box.zhi, opt.z);
    }
    return !req.targets.empty();
  };
  // Future cost of one search (Sec. 4.1 / arXiv:2111.06169 Table 3).
  const br::GCellGrid gcells{tile_x0_,
                             tile_y0_,
                             tile_dx_,
                             tile_dy_,
                             die_.xMin(),
                             die_.yMin(),
                             die_.xMax(),
                             die_.yMax()};
  auto choose_potential = [&](const bool path_region) {
    req.future_cost = nullptr;
    using Model = br::PotentialOracle::Model;
    Model model;
    switch (params_.potential) {
      case 1:
        model = Model::kL1;
        break;
      case 2:
        model = Model::kSimple;
        break;
      case 3:
        model = Model::kGeneral;
        break;
      case 4:
        if (!path_region
            || worker.path_detour < params_.potential_detour_tiles) {
          return;
        }
        model = Model::kGeneral;
        break;
      default:
        return;
    }
    // The general model bounds distances in the routing area; a source
    // outside it could leave through edges the model does not contain.
    if (model == Model::kGeneral
        && !std::ranges::all_of(req.sources, [&](const auto& source) {
             return in_region(source.first);
           })) {
      worker.potential_fallbacks++;
      return;
    }
    worker.potential.prepare(model,
                             graph_,
                             req,
                             gcells,
                             params_.discount_numerator,
                             params_.discount_denominator);
    req.future_cost = &worker.potential;
  };
  auto run = [&](const bool path_region) {
    choose_potential(path_region);
    br::SearchResult res = params_.use_intervals
                               ? worker.interval_search->run(req)
                               : worker.search->run(req);
    if (params_.potential_check && req.future_cost != nullptr) {
      // A feasible future cost leaves the shortest path cost unchanged.
      const br::FutureCostOracle* oracle = req.future_cost;
      req.future_cost = nullptr;
      const br::SearchResult plain = worker.search->run(req);
      req.future_cost = oracle;
      worker.potential_checks++;
      if (plain.found != res.found || (res.found && plain.cost != res.cost)) {
        worker.potential_mismatches++;
      }
    }
    worker.searches++;
    worker.heap_pops += res.heap_pops;
    worker.settled += res.settled;
    return res;
  };
  bool full_region = false;
  int remaining = num_terms - 1;
  while (remaining > 0) {
    br::SearchResult res;
    // Sec. 4.4: route the connection inside the GR tiles of a guide path
    // to another component first; fall back to the net's whole area.
    if (params_.component_regions && dilation == 0 && !ripup
        && buildPathRegion(net, connected, tree, valid, worker)) {
      full_region = false;
      req.sources.clear();
      for (const auto& source : sources) {
        if (in_region(source.first)) {
          req.sources.push_back(source);
        }
      }
      if (!req.sources.empty() && set_targets(true)) {
        res = run(true);
        worker.path_region_searches++;
      }
    }
    if (!res.found) {
      if (!full_region) {
        buildRegion(net, dilation, worker.region);
        full_region = true;
      }
      req.sources = sources;
      if (!set_targets(false)) {
        return false;
      }
      res = run(false);
    }
    if (!res.found) {
      return false;
    }
    const int t = res.target_id / 1024;
    const int o = res.target_id % 1024;
    net.terms[t].chosen = o;
    connected[t] = 1;
    remaining--;
    const bool first_connection = net.terms[0].chosen < 0;
    if (net.terms[0].chosen < 0) {
      const size_t s = res.path.back();
      br::Cost best = br::kInfCost;
      for (size_t k = 0; k < net.terms[0].options.size(); k++) {
        const AccessOption& opt = net.terms[0].options[k];
        if (valid[0][k] && opt.vertex == s && opt.cost < best) {
          best = opt.cost;
          net.terms[0].chosen = static_cast<int>(k);
        }
      }
      sources.clear();
    }
    if (!params_.postprocess) {
      pathToGeometry(res.path, net);
      for (const size_t v : res.path) {
        if (in_tree.insert(v).second) {
          tree.push_back(v);
          sources.emplace_back(v, 0);
        }
      }
      continue;
    }
    // Sec. 4.4: combine the on-track path with the off-track paths at its
    // ends and post-process the result for same-net rules.
    std::vector<br::RouteNode> chain;
    std::vector<OutVia> hints;
    const AccessOption& target = net.terms[t].options[o];
    chain = target.route;
    hints = target.vias;
    for (size_t k = 0; k < res.path.size(); k++) {
      int z, i, j;
      graph_.decode(res.path[k], z, i, j);
      const auto [x, y] = graph_.layers[z].point(i, j);
      const br::RouteNode node{layers_[z].fr, {x, y}};
      if (chain.empty() || !(chain.back() == node)) {
        chain.push_back(node);
      }
    }
    consumed[t] = 1;
    if (first_connection && net.terms[0].chosen >= 0) {
      const AccessOption& source = net.terms[0].options[net.terms[0].chosen];
      for (auto it = source.route.rbegin(); it != source.route.rend(); ++it) {
        if (chain.empty() || !(chain.back() == *it)) {
          chain.push_back(*it);
        }
      }
      hints.insert(hints.end(), source.vias.begin(), source.vias.end());
      consumed[0] = 1;
    }
    std::vector<OutWire> patches;
    postprocessConnection(chain, hints, net.id, worker, patches);
    routeToGeometry(chain, hints, net);
    net.wires.insert(net.wires.end(), patches.begin(), patches.end());
    // Later connections attach to track vertices that still carry metal.
    auto on_chain = [&](const frLayerNum layer, const br::GridPoint& p) {
      for (size_t k = 0; k < chain.size(); k++) {
        const br::RouteNode& a = chain[k];
        if (a.layer == layer && a.p == p) {
          return true;
        }
        if (k + 1 < chain.size() && chain[k + 1].layer == layer
            && a.layer == layer) {
          const br::RouteNode& b = chain[k + 1];
          if (p.x >= std::min(a.p.x, b.p.x) && p.x <= std::max(a.p.x, b.p.x)
              && p.y >= std::min(a.p.y, b.p.y)
              && p.y <= std::max(a.p.y, b.p.y)) {
            return true;
          }
        }
      }
      return false;
    };
    for (const size_t v : res.path) {
      int z, i, j;
      graph_.decode(v, z, i, j);
      const auto [x, y] = graph_.layers[z].point(i, j);
      if (on_chain(layers_[z].fr, {x, y}) && in_tree.insert(v).second) {
        tree.push_back(v);
        sources.emplace_back(v, 0);
      }
    }
  }
  for (Term& term : net.terms) {
    if (term.chosen < 0) {
      return false;
    }
    const AccessOption& opt = term.options[term.chosen];
    const size_t t = &term - net.terms.data();
    if (!consumed[t]) {
      net.wires.insert(net.wires.end(), opt.wires.begin(), opt.wires.end());
      net.vias.insert(net.vias.end(), opt.vias.begin(), opt.vias.end());
    }
    if (opt.planar) {
      net.pin_points.emplace_back(opt.ap_layer, opt.ap);
    }
  }
  // History for rip-up (Sec. 4.2): vertices whose use forces a rip-up.
  if (ripup) {
    for (const size_t v : in_tree) {
      if (!grid_.usable(v, br::kWire, net.id)) {
        grid_.bumpHistory(v);
      }
    }
  }
  return true;
}

void BonnRouteDR::commitNet(NetData& net, Worker* local)
{
  shapesOf(net, net.shapes);
  for (const RShape& s : net.shapes) {
    markShape(s, 2, net.id, 1);
    if (s.layer >= 0 && s.layer < static_cast<int>(routed_tree_.size())) {
      if (local != nullptr) {
        local->local_routed[s.layer].insert({s.rect, net.id});
      } else {
        routed_tree_[s.layer].insert({s.rect, net.id});
      }
    }
  }
  if (local != nullptr) {
    local->committed.push_back(net.id);
  }
  net.routed = true;
}

void BonnRouteDR::mergeWorker(Worker& worker)
{
  for (const int id : worker.committed) {
    for (const RShape& s : nets_[id].shapes) {
      if (s.layer >= 0 && s.layer < static_cast<int>(routed_tree_.size())) {
        routed_tree_[s.layer].insert({s.rect, id});
      }
    }
  }
  for (const auto& [shape, net] : worker.released_shapes) {
    reservation_tree_[shape.layer].remove(std::make_pair(shape.rect, net));
  }
  for (RTree<int>& tree : worker.local_routed) {
    tree.clear();
  }
  worker.released.clear();
  worker.released_shapes.clear();
  worker.committed.clear();
}

void BonnRouteDR::ripupNet(NetData& net)
{
  for (const RShape& s : net.shapes) {
    markShape(s, 2, net.id, -1);
    if (s.layer >= 0 && s.layer < static_cast<int>(routed_tree_.size())) {
      routed_tree_[s.layer].remove(std::make_pair(s.rect, net.id));
    }
  }
  net.shapes.clear();
  net.wires.clear();
  net.vias.clear();
  net.routed = false;
}

void BonnRouteDR::collectVictims(const NetData& net,
                                 std::vector<int>& victims) const
{
  std::vector<RShape> shapes;
  shapesOf(net, shapes);
  std::unordered_set<int> found;
  for (const RShape& shape : shapes) {
    if (shape.layer < 0 || shape.layer >= static_cast<int>(z_of_fr_.size())) {
      continue;
    }
    const int z = z_of_fr_[shape.layer];
    const int cz = cut_z_of_fr_[shape.layer];
    if (z < 0 && cz < 0) {
      continue;
    }
    odb::Rect query = shape.rect;
    query.bloat(max_spacing_[shape.layer], query);
    const auto& tree = routed_tree_[shape.layer];
    for (auto it = tree.qbegin(bgi::intersects(query)); it != tree.qend();
         ++it) {
      if (it->second == net.id || found.contains(it->second)) {
        continue;
      }
      const odb::Rect& other = it->first;
      frCoord s;
      if (z >= 0) {
        s = spacing(z,
                    other.minDXDY(),
                    shape.rect.minDXDY(),
                    std::min(other.maxDXDY(), shape.rect.maxDXDY()));
      } else {
        s = layers_[cz].cut_spacing;
      }
      const int64_t dx = std::max({0,
                                   other.xMin() - shape.rect.xMax(),
                                   shape.rect.xMin() - other.xMax()});
      const int64_t dy = std::max({0,
                                   other.yMin() - shape.rect.yMax(),
                                   shape.rect.yMin() - other.yMax()});
      if (dx * dx + dy * dy < static_cast<int64_t>(s) * s) {
        found.insert(it->second);
      }
    }
  }
  victims.assign(found.begin(), found.end());
  std::ranges::sort(victims);
}

void BonnRouteDR::netTileBox(const NetData& net,
                             int& xlo,
                             int& ylo,
                             int& xhi,
                             int& yhi) const
{
  frBlock* block = design_->getTopBlock();
  xlo = ylo = std::numeric_limits<int>::max();
  xhi = yhi = -1;
  auto add = [&](int x, int y) {
    xlo = std::min(xlo, x);
    ylo = std::min(ylo, y);
    xhi = std::max(xhi, x);
    yhi = std::max(yhi, y);
  };
  for (const auto& guide : net.net->getGuides()) {
    const odb::Point a = block->getGCellIdx(guide->getBeginPoint());
    const odb::Point b = block->getGCellIdx(guide->getEndPoint());
    add(a.x(), a.y());
    add(b.x(), b.y());
  }
  for (const Term& term : net.terms) {
    for (const AccessOption& o : term.options) {
      const auto [x, y] = graph_.point(o.vertex);
      add(tileX(x), tileY(y));
    }
  }
  // Retry dilation plus one tile for the spacing halo of committed shapes.
  const int halo = params_.retry_dilation + 1;
  xlo -= halo;
  ylo -= halo;
  xhi += halo;
  yhi += halo;
}

void BonnRouteDR::routeParallel(const std::vector<int>& order,
                                std::vector<int>& failed)
{
  // Batches of nets whose tile boxes are pairwise disjoint are routed in
  // parallel against the committed state, then committed in net order.
  // Batch composition does not depend on the thread count.
  std::vector<uint32_t> stamp(static_cast<size_t>(tiles_x_) * tiles_y_, 0);
  uint32_t round = 0;
  std::vector<std::array<int, 4>> boxes(nets_.size());
  for (const int id : order) {
    auto& b = boxes[id];
    netTileBox(nets_[id], b[0], b[1], b[2], b[3]);
    b[0] = std::max(b[0], 0);
    b[1] = std::max(b[1], 0);
    b[2] = std::min(b[2], tiles_x_ - 1);
    b[3] = std::min(b[3], tiles_y_ - 1);
  }
  std::vector<int> pending = order;
  std::vector<int> batch;
  std::vector<int> next;
  std::vector<char> ok;
  const int threads = static_cast<int>(workers_.size());
  while (!pending.empty()) {
    round++;
    batch.clear();
    next.clear();
    for (const int id : pending) {
      const auto [xlo, ylo, xhi, yhi] = boxes[id];
      bool free = true;
      for (int x = xlo; x <= xhi && free; x++) {
        for (int y = ylo; y <= yhi; y++) {
          if (stamp[static_cast<size_t>(x) * tiles_y_ + y] == round) {
            free = false;
            break;
          }
        }
      }
      if (!free) {
        next.push_back(id);
        continue;
      }
      for (int x = xlo; x <= xhi; x++) {
        for (int y = ylo; y <= yhi; y++) {
          stamp[static_cast<size_t>(x) * tiles_y_ + y] = round;
        }
      }
      batch.push_back(id);
    }
    ok.assign(batch.size(), 0);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (int k = 0; k < static_cast<int>(batch.size()); k++) {
      Worker& w = workers_[omp_get_thread_num()];
      NetData& net = nets_[batch[k]];
      ok[k] = connect(net, false, w);
    }
    for (size_t k = 0; k < batch.size(); k++) {
      NetData& net = nets_[batch[k]];
      if (ok[k]) {
        commitNet(net);
        // Sec. 4.4: reservations of connected components are released.
        setReservations(net, false);
      } else {
        net.wires.clear();
        net.vias.clear();
        failed.push_back(net.id);
      }
    }
    pending.swap(next);
  }
}

void BonnRouteDR::partitionRegions(
    const std::vector<int>& nets,
    const std::vector<std::array<int, 4>>& boxes,
    const int count,
    const bool shifted,
    std::vector<std::array<int, 4>>& regions) const
{
  // Columns, then rows within each column, at weighted quantiles of the
  // pending nets' box centres, so the estimated workload (box half
  // perimeter) per region is balanced. Shifted partitions cut halfway
  // between the previous cuts, so nets crossing an earlier border fit.
  regions.clear();
  const double aspect = static_cast<double>(tiles_x_) / tiles_y_;
  const int kx = std::clamp(
      static_cast<int>(std::lround(std::sqrt(count * aspect))), 1, count);
  const int ky = std::max(1, (count + kx - 1) / kx);
  auto weight = [&](int id) {
    const auto& b = boxes[id];
    return static_cast<double>(b[2] - b[0] + b[3] - b[1] + 2);
  };
  auto cuts = [&](std::vector<std::pair<int, double>>& items,
                  int parts,
                  int lo,
                  int hi) {
    std::ranges::sort(items);
    double total = 0;
    for (const auto& item : items) {
      total += item.second;
    }
    std::vector<int> result{lo};
    double sum = 0;
    size_t k = 0;
    for (int part = 1; part <= parts; part++) {
      const double q = (shifted ? part - 0.5 : part) / parts;
      if (part == parts && !shifted) {
        break;
      }
      while (k < items.size() && sum + items[k].second <= q * total) {
        sum += items[k].second;
        k++;
      }
      const int c = k < items.size() ? items[k].first : hi + 1;
      if (c > result.back() && c <= hi) {
        result.push_back(c);
      }
    }
    result.push_back(hi + 1);
    return result;
  };
  std::vector<std::pair<int, double>> xs;
  for (const int id : nets) {
    xs.emplace_back((boxes[id][0] + boxes[id][2]) / 2, weight(id));
  }
  const std::vector<int> xcuts = cuts(xs, kx, 0, tiles_x_ - 1);
  for (size_t c = 0; c + 1 < xcuts.size(); c++) {
    const int x0 = xcuts[c];
    const int x1 = xcuts[c + 1] - 1;
    std::vector<std::pair<int, double>> ys;
    for (const int id : nets) {
      const int cx = (boxes[id][0] + boxes[id][2]) / 2;
      if (cx >= x0 && cx <= x1) {
        ys.emplace_back((boxes[id][1] + boxes[id][3]) / 2, weight(id));
      }
    }
    const std::vector<int> ycuts = cuts(ys, ky, 0, tiles_y_ - 1);
    for (size_t r = 0; r + 1 < ycuts.size(); r++) {
      regions.push_back({x0, ycuts[r], x1, ycuts[r + 1] - 1});
    }
  }
}

void BonnRouteDR::routePartitions(const std::vector<int>& order,
                                  std::vector<int>& failed,
                                  const int initial_regions,
                                  const char* stage)
{
  std::vector<std::array<int, 4>> boxes(nets_.size());
  for (const int id : order) {
    auto& b = boxes[id];
    netTileBox(nets_[id], b[0], b[1], b[2], b[3]);
    b[0] = std::max(b[0], 0);
    b[1] = std::max(b[1], 0);
    b[2] = std::min(b[2], tiles_x_ - 1);
    b[3] = std::min(b[3], tiles_y_ - 1);
  }
  const int threads = static_cast<int>(workers_.size());
  int count = initial_regions > 0             ? initial_regions
              : params_.partition_regions > 0 ? params_.partition_regions
                                              : 2 * threads;
  std::vector<int> pending = order;
  std::vector<std::array<int, 4>> regions;
  int partition = 0;
  while (!pending.empty()) {
    if (count <= 1) {
      regions.assign(1, {0, 0, tiles_x_ - 1, tiles_y_ - 1});
    } else {
      partitionRegions(pending, boxes, count, partition % 2 == 1, regions);
    }
    std::vector<std::vector<int>> members(regions.size());
    std::vector<std::vector<int>> region_failed(regions.size());
    std::vector<int> next;
    for (const int id : pending) {
      const auto& b = boxes[id];
      bool placed = false;
      for (size_t r = 0; r < regions.size() && !placed; r++) {
        const auto& g = regions[r];
        if (b[0] >= g[0] && b[1] >= g[1] && b[2] <= g[2] && b[3] <= g[3]) {
          members[r].push_back(id);
          placed = true;
        }
      }
      if (!placed) {
        next.push_back(id);
      }
    }
    if (count > 1 && next.size() == pending.size()) {
      // No pending net fits a region of this partition: continue with the
      // next partition of the sequence without routing this one.
      partition++;
      if (partition % 2 == 0) {
        count = std::max(1, count / 2);
      }
      continue;
    }
    // Each region is routed by one thread, which changes routing only
    // inside it; its nets are routed and committed in order.
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (int r = 0; r < static_cast<int>(regions.size()); r++) {
      Worker& w = workers_[omp_get_thread_num()];
      for (const int id : members[r]) {
        NetData& net = nets_[id];
        if (connect(net, false, w)) {
          commitNet(net, &w);
          setReservations(net, false, &w);
        } else {
          net.wires.clear();
          net.vias.clear();
          region_failed[r].push_back(id);
        }
      }
    }
    for (Worker& w : workers_) {
      mergeWorker(w);
    }
    for (const auto& f : region_failed) {
      failed.insert(failed.end(), f.begin(), f.end());
    }
    if (params_.verbose) {
      size_t routed = 0;
      for (const auto& m : members) {
        routed += m.size();
      }
      logger_->info(utl::DRT,
                    9103,
                    "BonnRoute {} {}: {} regions, {} nets inside, {} "
                    "deferred.",
                    stage,
                    partition,
                    regions.size(),
                    routed,
                    next.size());
    }
    pending.swap(next);
    partition++;
    // Fewer, larger regions in later partitions (Sec. 5.1).
    if (count > 1 && partition % 2 == 0) {
      count = std::max(1, count / 2);
    }
  }
}

BonnRouteDR::SavedRouting BonnRouteDR::saveRouting(const NetData& net) const
{
  SavedRouting saved;
  saved.net = net.id;
  saved.was_routed = net.routed;
  saved.wires = net.wires;
  saved.vias = net.vias;
  saved.pin_points = net.pin_points;
  saved.input_reserved = net.input_reserved;
  for (const Term& term : net.terms) {
    saved.chosen.push_back(term.chosen);
    saved.reserved.push_back(term.reservation_active);
  }
  return saved;
}

void BonnRouteDR::restoreRouting(const SavedRouting& saved)
{
  NetData& net = nets_[saved.net];
  if (net.routed) {
    ripupNet(net);
  }
  net.wires = saved.wires;
  net.vias = saved.vias;
  net.pin_points = saved.pin_points;
  for (size_t t = 0; t < net.terms.size(); t++) {
    net.terms[t].chosen = saved.chosen[t];
    setReservation(net.terms[t], net.id, saved.reserved[t] != 0);
  }
  setInputReservation(net, saved.input_reserved);
  if (saved.was_routed) {
    commitNet(net);
  } else {
    net.wires.clear();
    net.vias.clear();
  }
}

bool BonnRouteDR::ripupSequence(const int id,
                                const int depth,
                                const int dilation,
                                std::vector<SavedRouting>& journal)
{
  NetData& net = nets_[id];
  Worker& worker = workers_.front();
  journal.push_back(saveRouting(net));
  if (connect(net, false, worker, dilation)) {
    commitNet(net);
    setReservations(net, false);
    return true;
  }
  if (depth >= params_.ripup_depth || !connect(net, true, worker, dilation)) {
    return false;
  }
  std::vector<int> victims;
  collectVictims(net, victims);
  for (const int victim : victims) {
    if (nets_[victim].ripups >= params_.max_ripups) {
      return false;
    }
  }
  // Remove the victims, add the new route, then reconnect the victims,
  // again possibly with rip-up (Sec. 4.4).
  for (const int victim : victims) {
    journal.push_back(saveRouting(nets_[victim]));
    ripupNet(nets_[victim]);
    nets_[victim].ripups++;
    ripup_events_++;
  }
  commitNet(net);
  setReservations(net, false);
  for (const int victim : victims) {
    if (!nets_[victim].routed
        && !ripupSequence(victim, depth + 1, dilation, journal)) {
      return false;
    }
  }
  return true;
}

void BonnRouteDR::routeRipupQueue(const std::vector<int>& failed)
{
  std::vector<int> queue = failed;
  std::vector<SavedRouting> journal;
  for (int pass = 0; pass < params_.ripup_passes && !queue.empty(); pass++) {
    std::vector<int> unrouted;
    for (const int id : queue) {
      if (nets_[id].routed) {
        continue;
      }
      journal.clear();
      if (ripupSequence(id, 0, pass, journal)) {
        continue;
      }
      // Sec. 4.4: a failed sequence restores the original situation.
      for (auto it = journal.rbegin(); it != journal.rend(); ++it) {
        restoreRouting(*it);
      }
      rollbacks_++;
      unrouted.push_back(id);
    }
    queue.swap(unrouted);
  }
  for (const int id : queue) {
    NetData& net = nets_[id];
    net.wires.clear();
    net.vias.clear();
    net.failed = true;
  }
}

void BonnRouteDR::writeNets()
{
  frTechObject* tech = design_->getTech();
  for (NetData& data : nets_) {
    if (!data.routed) {
      continue;
    }
    frNet* net = data.net;
    // Split wires at planar pin points so the pin end is an endpoint.
    std::vector<OutWire> wires;
    for (const OutWire& w : data.wires) {
      std::vector<odb::Point> cuts{w.begin, w.end};
      for (const auto& [layer, p] : data.pin_points) {
        const bool on_x = w.begin.y() == w.end.y() && p.y() == w.begin.y()
                          && p.x() > std::min(w.begin.x(), w.end.x())
                          && p.x() < std::max(w.begin.x(), w.end.x());
        const bool on_y = w.begin.x() == w.end.x() && p.x() == w.begin.x()
                          && p.y() > std::min(w.begin.y(), w.end.y())
                          && p.y() < std::max(w.begin.y(), w.end.y());
        if (layer == w.layer && (on_x || on_y)) {
          cuts.push_back(p);
        }
      }
      std::ranges::sort(cuts, [](const odb::Point& a, const odb::Point& b) {
        return std::make_pair(a.x(), a.y()) < std::make_pair(b.x(), b.y());
      });
      for (size_t k = 0; k + 1 < cuts.size(); k++) {
        wires.push_back({w.layer, cuts[k], cuts[k + 1]});
      }
    }
    for (const OutWire& w : wires) {
      if (w.begin == w.end) {
        continue;
      }
      auto seg = std::make_unique<frPathSeg>();
      const odb::Point b(std::min(w.begin.x(), w.end.x()),
                         std::min(w.begin.y(), w.end.y()));
      const odb::Point e(std::max(w.begin.x(), w.end.x()),
                         std::max(w.begin.y(), w.end.y()));
      seg->setPoints(b, e);
      seg->setLayerNum(w.layer);
      frSegStyle style = tech->getLayer(w.layer)->getDefaultSegStyle();
      // A wire ending on a pin is truncated there, which is how FlexDR
      // recognises planar pin connections.
      for (const auto& [layer, p] : data.pin_points) {
        if (layer != w.layer) {
          continue;
        }
        if (p == b) {
          style.setBeginStyle(frcTruncateEndStyle, 0);
        }
        if (p == e) {
          style.setEndStyle(frcTruncateEndStyle, 0);
        }
      }
      seg->setStyle(style);
      net->addShape(std::move(seg));
    }
    for (const OutVia& v : data.vias) {
      net->addVia(std::make_unique<frVia>(v.def, v.point));
    }
    net->setHasInitialRouting(true);
    net->setModified(true);
  }
}

void BonnRouteDR::main()
{
  utl::Timer timer;
  if (modernSearch() && params_.use_intervals) {
    // The interval search assumes unit wire costs and the l1 future cost.
    logger_->warn(utl::DRT,
                  9102,
                  "BonnRoute: reservations and potentials use the vertex "
                  "search; use_intervals is ignored.");
    params_.use_intervals = false;
  }
  buildTrackGraph();
  if (layers_.empty()) {
    return;
  }
  grid_.init(graph_.numVertices());
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
  buildFixed();
  buildSameNetRules();
  buildLowLayerIndex();
  std::vector<char> pre_skip(nets_.size());
  for (const NetData& net : nets_) {
    pre_skip[net.id] = net.skip;
  }
  for (NetData& net : nets_) {
    if (!net.skip) {
      buildGuideTiles(net);
      buildAccess(net);
    }
  }
  if (params_.offtrack_access) {
    buildOffTrackAccess(pre_skip);
  } else {
    for (NetData& net : nets_) {
      if (!net.skip) {
        reserveAccess(net);
      }
    }
  }
  if (params_.reservations) {
    reserveInputWires();
  }
  const double setup_time = timer.elapsed();
  std::vector<int> order;
  for (const NetData& net : nets_) {
    if (!net.skip) {
      order.push_back(net.id);
    }
  }
  std::ranges::stable_sort(
      order, [&](int a, int b) { return nets_[a].hpwl < nets_[b].hpwl; });
  std::vector<int> failed_first;
  for (Worker& w : workers_) {
    w.local_routed.assign(z_of_fr_.size(), RTree<int>());
  }
  // Sec. 5.1: critical nets (here: clock nets) first, on the whole chip.
  std::vector<int> regular;
  std::vector<int> critical;
  for (const int id : order) {
    (params_.critical_nets && nets_[id].net->isClock() ? critical : regular)
        .push_back(id);
  }
  if (!critical.empty()) {
    routePartitions(critical, failed_first, 1, "critical-net partition");
  }
  if (params_.partitions) {
    routePartitions(regular, failed_first);
  } else {
    routeParallel(regular, failed_first);
  }
  const double parallel_time = timer.elapsed();
  // FlexDR handles unresolved nets unless serial rip-up is requested.
  if (params_.ripup) {
    routeRipupQueue(failed_first);
  }
  writeNets();
  int routed = 0;
  int failed = 0;
  int skipped = 0;
  for (const NetData& net : nets_) {
    routed += net.routed;
    failed += !net.routed && !net.skip;
    skipped += net.skip;
  }
  int64_t searches = 0;
  int64_t heap_pops = 0;
  int64_t settled = 0;
  int64_t intervals = 0;
  int64_t labels = 0;
  int64_t path_region_searches = 0;
  for (const Worker& w : workers_) {
    searches += w.searches;
    path_region_searches += w.path_region_searches;
    heap_pops += w.heap_pops;
    settled += w.settled;
    if (w.interval_search) {
      intervals += w.interval_search->intervals();
      labels += w.interval_search->labels();
    }
  }
  if (params_.verbose) {
    logger_->report(
        "BonnRoute DR: {} threads, parallel phase until {:.2f}s, {} nets "
        "to {}.",
        workers_.size(),
        parallel_time,
        failed_first.size(),
        params_.ripup ? "the rip-up phase" : "FlexDR");
    logger_->report(
        "BonnRoute DR: {} layers, {} track-graph vertices, setup {:.2f}s.",
        graph_.numLayers(),
        graph_.numVertices(),
        setup_time);
    if (params_.potential != 0) {
      int64_t prepares[3] = {0, 0, 0};
      int64_t fallbacks = 0;
      int64_t checks = 0;
      int64_t mismatches = 0;
      for (const Worker& w : workers_) {
        fallbacks += w.potential_fallbacks;
        checks += w.potential_checks;
        mismatches += w.potential_mismatches;
        for (int m = 0; m < 3; m++) {
          prepares[m] += w.potential.prepares(
              static_cast<br::PotentialOracle::Model>(m));
        }
      }
      logger_->report(
          "BonnRoute DR: potential {}: {} l1, {} simple, {} general "
          "preparations, {} general fallbacks to pi_H, {} of {} checked "
          "searches differ from pi_H.",
          params_.potential,
          prepares[0],
          prepares[1],
          prepares[2],
          fallbacks,
          mismatches,
          checks);
    }
    if (params_.spreading) {
      logger_->report(
          "BonnRoute DR: wire spreading in {} of {} layer tiles (GR "
          "utilization <= {}%).",
          spread_tile_count_,
          static_cast<int64_t>(graph_.numLayers()) * tiles_x_ * tiles_y_,
          params_.spread_utilization);
    }
    if (params_.reservations) {
      logger_->report(
          "BonnRoute DR: {} input wires, {} reserved runs, discount {}/{}.",
          input_wires_,
          input_reserved_,
          params_.discount_numerator,
          params_.discount_denominator);
    }
    logger_->report(
        "BonnRoute DR: routed {} nets, {} left to FlexDR, {} skipped; {} "
        "searches ({}), {} heap pops, {} settled vertices, {} intervals, "
        "{} interval labels, {} rip-ups, {} restored rip-up sequences, "
        "{} guide-path region searches, {:.2f}s.",
        routed,
        failed,
        skipped,
        searches,
        params_.use_intervals ? "interval labels" : "vertex labels",
        heap_pops,
        settled,
        intervals,
        labels,
        ripup_events_,
        rollbacks_,
        path_region_searches,
        timer.elapsed());
    br::CleanupStats cleanup;
    int64_t regenerations = 0;
    int64_t regenerated_paths = 0;
    for (const Worker& w : workers_) {
      cleanup.add(w.cleanup);
      regenerations += w.regenerations;
      regenerated_paths += w.regenerated_paths;
    }
    if (params_.offtrack_access) {
      logger_->report(
          "BonnRoute DR: {} dynamic access path generations added {} paths.",
          regenerations,
          regenerated_paths);
    }
    if (params_.postprocess) {
      logger_->report(
          "BonnRoute DR post-processing: {} layer runs, same-net violations "
          "{} -> {} ({} jog shifts, {} notch collapses), {} min-area pieces, "
          "{} extended.",
          cleanup.runs,
          cleanup.violations_before,
          cleanup.violations_after,
          cleanup.moved_segments,
          cleanup.collapsed_notches,
          cleanup.area_violations,
          cleanup.area_fixed);
    }
  }
}

}  // namespace drt
