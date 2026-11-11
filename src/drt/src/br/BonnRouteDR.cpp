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

#include "br/BrCore.h"
#include "db/obj/frAccess.h"
#include "db/obj/frBTerm.h"
#include "db/obj/frBlock.h"
#include "db/obj/frInst.h"
#include "db/obj/frInstTerm.h"
#include "db/obj/frNet.h"
#include "db/obj/frShape.h"
#include "db/obj/frVia.h"
#include "db/tech/frLayer.h"
#include "db/tech/frTechObject.h"
#include "db/tech/frViaDef.h"
#include "drt-global.h"
#include "frDesign.h"
#include "frRegionQuery.h"
#include "odb/db.h"
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
  for (int t = 0; t < num_terms; t++) {
    const Term& term = net.terms[t];
    valid[t].assign(term.options.size(), 0);
    for (size_t o = 0; o < term.options.size(); o++) {
      const AccessOption& opt = term.options[o];
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
    pathToGeometry(res.path, net);
    for (const size_t v : res.path) {
      if (in_tree.insert(v).second) {
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
    net.wires.insert(net.wires.end(), opt.wires.begin(), opt.wires.end());
    net.vias.insert(net.vias.end(), opt.vias.begin(), opt.vias.end());
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
                                  const int initial_regions)
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
                    "BonnRoute partition {}: {} regions, {} nets inside, {} "
                    "deferred.",
                    partition,
                    regions.size(),
                    routed,
                    next.size());
    }
    const bool placed_any = next.size() < pending.size();
    pending.swap(next);
    partition++;
    // Fewer, larger regions in later partitions (Sec. 5.1). A partition
    // that holds none of the pending nets is skipped in the sequence.
    if (count > 1 && (partition % 2 == 0 || !placed_any)) {
      count = std::max(1, count / 2);
      partition += partition % 2;
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
  for (NetData& net : nets_) {
    if (!net.skip) {
      buildGuideTiles(net);
      buildAccess(net);
    }
  }
  for (NetData& net : nets_) {
    if (!net.skip) {
      reserveAccess(net);
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
    routePartitions(critical, failed_first, 1);
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
      for (const Worker& w : workers_) {
        fallbacks += w.potential_fallbacks;
        for (int m = 0; m < 3; m++) {
          prepares[m] += w.potential.prepares(
              static_cast<br::PotentialOracle::Model>(m));
        }
      }
      logger_->report(
          "BonnRoute DR: potential {}: {} l1, {} simple, {} general "
          "preparations, {} general fallbacks to pi_H.",
          params_.potential,
          prepares[0],
          prepares[1],
          prepares[2],
          fallbacks);
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
  }
}

}  // namespace drt
