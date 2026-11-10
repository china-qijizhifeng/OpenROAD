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
                         const BonnRouteDRParams& params)
    : design_(design), logger_(logger), router_cfg_(router_cfg), params_(params)
{
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
    tl.jog_factor = params_.jog_factor;
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
                              * std::max(info.pitch, 1));
  }
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
  fixed_tree_.resize(num_fr);
  routed_tree_.resize(num_fr);
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
          } else {
            grid_.addRouted(v, kind, delta);
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
            } else {
              grid_.addRouted(v, edge, delta);
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

bool BonnRouteDR::shapeConflicts(const RShape& shape, const int net) const
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
  auto conflicts = [&](const odb::Rect& other) {
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
    return dx * dx + dy * dy < static_cast<int64_t>(s) * s;
  };
  for (const auto* tree :
       {&fixed_tree_[shape.layer], &routed_tree_[shape.layer]}) {
    for (auto it = tree->qbegin(bgi::intersects(query)); it != tree->qend();
         ++it) {
      if (it->second == net && net >= 0) {
        continue;
      }
      if (conflicts(it->first)) {
        return true;
      }
    }
  }
  return false;
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
                                         * info.pitch)
                          : 0)
                     + static_cast<Cost>(params_.jog_factor) * jog + run
                     + (a.preferred
                            ? 0
                            : static_cast<Cost>(params_.pin_cost) * info.pitch);
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
  // Sec. 4.3: the primary access via of each pin is reserved before
  // routing so later nets do not block it.
  for (const Term& term : net.terms) {
    if (term.options.empty()) {
      continue;
    }
    for (const OutVia& v : term.options.front().vias) {
      std::vector<RShape> vs;
      viaShapes(v, vs);
      const RShape& pad = vs.back();
      markShape(pad, 1, net.id, 0);
      fixed_tree_[pad.layer].insert({pad.rect, net.id});
    }
  }
}

void BonnRouteDR::buildRegion(const NetData& net,
                              const int dilation,
                              br::Region& region)
{
  region.reset();
  const int num_z = graph_.numLayers();
  frBlock* block = design_->getTopBlock();
  auto allow = [&](int z, int tx, int ty) {
    for (int dz = -params_.layer_dilation; dz <= params_.layer_dilation; dz++) {
      for (int dx = -dilation; dx <= dilation; dx++) {
        for (int dy = -dilation; dy <= dilation; dy++) {
          region.allow(z + dz, tx + dx, ty + dy);
        }
      }
    }
  };
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
        allow(z, x, y);
      }
    }
  }
  for (const Term& term : net.terms) {
    for (const AccessOption& o : term.options) {
      const auto [x, y] = graph_.point(o.vertex);
      allow(o.z, tileX(x), tileY(y));
    }
  }
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
                          const int dilation,
                          const bool ripup,
                          Worker& worker)
{
  net.wires.clear();
  net.vias.clear();
  net.pin_points.clear();
  for (Term& t : net.terms) {
    t.chosen = -1;
  }
  buildRegion(net, dilation, worker.region);
  const int num_terms = static_cast<int>(net.terms.size());
  // Options usable now: stub free of routed foreign shapes (unless rip-up)
  // and an enterable attachment vertex.
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
          if (shapeConflicts(wireShape(w), net.id)) {
            ok = false;
            break;
          }
        }
        // Landing pad against routed shapes and other pins' reservations.
        for (const OutVia& v : opt.vias) {
          if (!ok) {
            break;
          }
          std::vector<RShape> vs;
          viaShapes(v, vs);
          ok = !shapeConflicts(vs.back(), net.id);
        }
      }
      valid[t][o] = ok;
    }
    if (std::ranges::find(valid[t], 1) == valid[t].end()) {
      return false;
    }
  }
  br::SearchRequest req;
  req.net = net.id;
  req.region = &worker.region;
  req.ripup = ripup;
  req.ripup_cost = static_cast<Cost>(params_.ripup_cost) * layers_[0].pitch;
  req.use_intervals = params_.use_intervals;
  std::vector<size_t> tree;
  std::unordered_set<size_t> in_tree;
  std::vector<char> connected(num_terms, 0);
  connected[0] = 1;
  for (size_t o = 0; o < net.terms[0].options.size(); o++) {
    if (valid[0][o]) {
      req.sources.emplace_back(net.terms[0].options[o].vertex,
                               net.terms[0].options[o].cost);
    }
  }
  int remaining = num_terms - 1;
  while (remaining > 0) {
    req.targets.clear();
    req.target_boxes.clear();
    for (int t = 1; t < num_terms; t++) {
      if (connected[t]) {
        continue;
      }
      for (size_t o = 0; o < net.terms[t].options.size(); o++) {
        if (valid[t][o]) {
          const AccessOption& opt = net.terms[t].options[o];
          req.targets.push_back(
              {opt.vertex, opt.cost, t * 1024 + static_cast<int>(o)});
        }
      }
    }
    // Future cost pi_H (Sec. 4.1): l1 distance plus via lower bound to the
    // nearest target point; per-term boxes keep it cheap for large nets.
    const bool per_point = req.targets.size() <= 64;
    for (int t = 1; t < num_terms; t++) {
      if (connected[t]) {
        continue;
      }
      br::TargetBox box{std::numeric_limits<int>::max(),
                        std::numeric_limits<int>::max(),
                        std::numeric_limits<int>::min(),
                        std::numeric_limits<int>::min(),
                        std::numeric_limits<int>::max(),
                        -1};
      for (size_t o = 0; o < net.terms[t].options.size(); o++) {
        if (!valid[t][o]) {
          continue;
        }
        const AccessOption& opt = net.terms[t].options[o];
        const auto [x, y] = graph_.point(opt.vertex);
        if (per_point) {
          req.target_boxes.push_back({x, y, x, y, opt.z, opt.z});
          continue;
        }
        box.xlo = std::min(box.xlo, x);
        box.ylo = std::min(box.ylo, y);
        box.xhi = std::max(box.xhi, x);
        box.yhi = std::max(box.yhi, y);
        box.zlo = std::min(box.zlo, opt.z);
        box.zhi = std::max(box.zhi, opt.z);
      }
      if (!per_point) {
        req.target_boxes.push_back(box);
      }
    }
    br::SearchResult res = params_.use_intervals
                               ? worker.interval_search->run(req)
                               : worker.search->run(req);
    worker.searches++;
    worker.heap_pops += res.heap_pops;
    worker.settled += res.settled;
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
      req.sources.clear();
    }
    pathToGeometry(res.path, net);
    for (const size_t v : res.path) {
      if (in_tree.insert(v).second) {
        tree.push_back(v);
        req.sources.emplace_back(v, 0);
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
    for (const size_t v : tree) {
      if (!grid_.usable(v, br::kWire, net.id)) {
        grid_.bumpHistory(v);
      }
    }
  }
  return true;
}

void BonnRouteDR::commitNet(NetData& net)
{
  shapesOf(net, net.shapes);
  for (const RShape& s : net.shapes) {
    markShape(s, 2, net.id, 1);
    if (s.layer >= 0 && s.layer < static_cast<int>(routed_tree_.size())) {
      routed_tree_[s.layer].insert({s.rect, net.id});
    }
  }
  net.routed = true;
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

void BonnRouteDR::routeParallel(std::vector<int>& order,
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
      ok[k] = connect(net, 0, false, w)
              || connect(net, params_.retry_dilation, false, w);
    }
    for (size_t k = 0; k < batch.size(); k++) {
      NetData& net = nets_[batch[k]];
      if (ok[k]) {
        commitNet(net);
      } else {
        net.wires.clear();
        net.vias.clear();
        failed.push_back(net.id);
      }
    }
    pending.swap(next);
  }
}

bool BonnRouteDR::routeNet(NetData& net, const bool allow_ripup)
{
  if (net.skip || net.failed) {
    return false;
  }
  Worker& w = workers_[0];
  if (connect(net, 0, false, w)
      || connect(net, params_.retry_dilation, false, w)) {
    commitNet(net);
    return true;
  }
  if (!allow_ripup || !connect(net, params_.retry_dilation, true, w)) {
    net.wires.clear();
    net.vias.clear();
    net.failed = true;
    return false;
  }
  std::vector<int> victims;
  collectVictims(net, victims);
  for (const int v : victims) {
    ripupNet(nets_[v]);
    ripup_events_++;
  }
  commitNet(net);
  // Victims are rerouted by the caller.
  for (const int v : victims) {
    nets_[v].ripups++;
    if (nets_[v].ripups > params_.max_ripups) {
      nets_[v].failed = true;
    }
  }
  pending_victims_ = std::move(victims);
  return true;
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
  buildTrackGraph();
  if (layers_.empty()) {
    return;
  }
  grid_.init(graph_.numVertices());
  workers_.resize(std::max(1, router_cfg_->MAX_THREADS));
  for (Worker& w : workers_) {
    w.region.init(graph_.numLayers(), tiles_x_, tiles_y_);
    w.search = std::make_unique<br::PathSearch>(graph_, grid_);
    w.interval_search = std::make_unique<br::IntervalSearch>(graph_, grid_);
  }
  buildNets();
  buildFixed();
  for (NetData& net : nets_) {
    if (!net.skip) {
      buildAccess(net);
    }
  }
  for (NetData& net : nets_) {
    if (!net.skip) {
      reserveAccess(net);
    }
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
  routeParallel(order, failed_first);
  const double parallel_time = timer.elapsed();
  // Nets without a rip-up-free route: serial rip-up and reroute.
  std::deque<int> queue(failed_first.begin(), failed_first.end());
  while (!queue.empty()) {
    const int id = queue.front();
    queue.pop_front();
    NetData& net = nets_[id];
    if (net.routed || net.failed) {
      continue;
    }
    pending_victims_.clear();
    routeNet(net, params_.ripup);
    for (const int v : pending_victims_) {
      if (!nets_[v].failed) {
        queue.push_back(v);
      }
    }
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
  for (const Worker& w : workers_) {
    searches += w.searches;
    heap_pops += w.heap_pops;
    settled += w.settled;
    intervals += w.interval_search->intervals();
    labels += w.interval_search->labels();
  }
  if (params_.verbose) {
    logger_->report(
        "BonnRoute DR: {} threads, parallel phase until {:.2f}s, {} nets "
        "to the rip-up phase.",
        workers_.size(),
        parallel_time,
        failed_first.size());
    logger_->report(
        "BonnRoute DR: {} layers, {} track-graph vertices, setup {:.2f}s.",
        graph_.numLayers(),
        graph_.numVertices(),
        setup_time);
    logger_->report(
        "BonnRoute DR: routed {} nets, {} left to FlexDR, {} skipped; {} "
        "searches ({}), {} heap pops, {} settled vertices, {} intervals, "
        "{} interval labels, {} rip-ups, {:.2f}s.",
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
        timer.elapsed());
  }
}

}  // namespace drt
