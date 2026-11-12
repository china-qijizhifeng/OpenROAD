// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BonnRouteDR.h"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "br/BrBlockageGrid.h"
#include "br/BrCore.h"
#include "br/BrFixedShapes.h"
#include "br/BrGeometry.h"
#include "br/BrOffTrackSearch.h"
#include "br/BrPathCleanup.h"
#include "br/BrPinAccess.h"
#include "br/BrSameNetRules.h"
#include "br/BrTailPipeline.h"
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
#include "grt/BonnRouteExtraSpace.h"
#include "odb/db.h"
#include "odb/dbTransform.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"
#include "utl/Logger.h"
#include "utl/timer.h"

namespace drt {

using br::Cost;
using br::Kind;

void BonnRouteDRParams::applyProfile(const int which)
{
  profile = (which == 1 || which == 2) ? which : 0;
  if (profile >= 1) {
    // 2013-reimpl: Gester et al. 2013 as written. Interval search with
    // property (ii) for vias and the Algorithm 4 tie order (Sec. 4.1),
    // rip-up sequences (Sec. 4.4), the partition sequence to its end
    // (Sec. 5.1), pi_H with the exact via bound and the blockage-aware
    // pi_P on detours (Sec. 4.1), wire spreading (Sec. 4.2).
    use_intervals = true;
    interval_via_property = true;
    interval_nearest_first = false;
    ripup = true;
    tail_batch = 0;
    potential = 4;
    spreading = true;
  } else {
    // 2013-improved: measured optimizations (vertex search, speculative
    // tail batches, unresolved nets to FlexDR), and wire spreading (Sec.
    // 4.2): without it neighbours pack densely, and signal coupling
    // capacitance (+1.4..+4.0% against FlexDR) costs setup slack.
    use_intervals = false;
    interval_via_property = true;
    interval_nearest_first = false;
    ripup = false;
    tail_batch = 32;
    tail_order = true;
    tail_pipeline = false;
    tail_retry_delay = 1;
    tail_defer_settled = 300000;
    potential = 0;
    spreading = true;
  }
  partitions = true;
  critical_nets = true;
  // Sec. 3.2 / 5.1: NDR nets with their wire types in every profile. Left
  // to FlexDR, they are its largest clean-up (ibex: 17 clock nets) and, in
  // incremental calls, need track assignment.
  route_ndr_nets = true;
  component_regions = true;
  exact_future_cost = true;
  offtrack_access = true;
  postprocess = true;
  access_classes = true;
  reservations = false;
  if (profile == 2) {
    // 2021: the bulk routing experiment of Ahrens et al., arXiv:2111.06169
    // (Table 3) on top of 2013-reimpl: the general-model potential for every
    // search, in the interval search. Reservations (Sec. 5.2, Tables 4-6)
    // belong to incremental routing, whose input is a legal detailed
    // routing; `reservations 1` with discount 3/4 enables them separately.
    potential = 3;
    discount_numerator = 3;
    discount_denominator = 4;
  }
}

const char* BonnRouteDRParams::profileName() const
{
  return profile == 2 ? "2021"
                      : (profile == 1 ? "2013-reimpl" : "2013-improved");
}

bool BonnRouteDRParams::set(const std::string& name, const double value)
{
  const int v = static_cast<int>(value);
  if (name == "profile") {
    applyProfile(v);
  } else if (name == "via_cost") {
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
  } else if (name == "interval_nearest_first") {
    interval_nearest_first = v != 0;
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
  } else if (name == "exact_future_cost") {
    exact_future_cost = v != 0;
  } else if (name == "potential_check") {
    potential_check = v != 0;
  } else if (name == "potential_detour_tiles") {
    potential_detour_tiles = std::max(0, v);
  } else if (name == "critical_nets") {
    critical_nets = v != 0;
  } else if (name == "route_ndr_nets") {
    route_ndr_nets = v != 0;
  } else if (name == "fast_ndr_permille") {
    fast_ndr_permille = std::clamp(v, 0, 1000);
  } else if (name == "fast_wire_types") {
    fast_wire_types = std::clamp(v, 1, br::FastGrid::kMaxWireTypes);
  } else if (name == "tail_batch") {
    tail_batch = std::max(0, v);
  } else if (name == "tail_threads") {
    tail_threads = std::max(0, v);
  } else if (name == "tail_order") {
    tail_order = v != 0;
  } else if (name == "tail_pipeline") {
    tail_pipeline = v != 0;
  } else if (name == "tail_retry_delay") {
    tail_retry_delay = std::max(1, v);
  } else if (name == "tail_defer_settled") {
    tail_defer_settled = std::max(0, v);
  } else if (name == "partitions") {
    partitions = v != 0;
  } else if (name == "partition_regions") {
    partition_regions = std::max(0, v);
  } else if (name == "spread_from_extra_space") {
    spread_from_extra_space = v != 0;
  } else if (name == "spreading") {
    spreading = v != 0;
  } else if (name == "jog_factor_low") {
    jog_factor_low = std::max(0, v);
  } else if (name == "jog_low_utilization") {
    jog_low_utilization = std::clamp(v, 0, 100);
  } else if (name == "collapse_bridges") {
    collapse_bridges = v != 0;
  } else if (name == "collapse_lower_bridges") {
    collapse_lower_bridges = v != 0;
  } else if (name == "spread_unidirectional") {
    spread_unidirectional = v != 0;
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
  } else if (name == "access_class_margin") {
    access_class_margin = std::max(0, v);
  } else if (name == "offtrack_macro_pins") {
    offtrack_macro_pins = v != 0;
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
  } else if (name == "offtrack_from_access_points") {
    offtrack_from_access_points = v != 0;
  } else if (name == "offtrack_preferred_access_points") {
    offtrack_preferred_access_points = v != 0;
  } else if (name == "search_label_mb") {
    search_label_mb = std::max(0, v);
  } else if (name == "cleanup_incremental_iterations") {
    cleanup_incremental_iterations = std::clamp(v, 0, 2);
  } else if (name == "cleanup_check_only_workers") {
    cleanup_check_only_workers = v != 0;
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

void BonnRouteDR::loadExtraSpace()
{
  extra_space_nets_ = 0;
  extra_space_tiles_ = 0;
  odb::dbGCellGrid* gcells
      = block_ != nullptr ? block_->getGCellGrid() : nullptr;
  if (gcells == nullptr) {
    return;
  }
  std::vector<int> grid_x;
  std::vector<int> grid_y;
  gcells->getGridX(grid_x);
  gcells->getGridY(grid_y);
  odb::dbTech* tech = block_->getDb()->getTech();
  std::vector<grt::bonnroute::ExtraSpaceEdge> edges;
  for (NetData& net : nets_) {
    net.extra_space_tiles.clear();
    if (net.skip) {
      continue;
    }
    odb::dbNet* db_net = block_->findNet(net.net->getName().c_str());
    if (db_net == nullptr) {
      continue;
    }
    odb::dbStringProperty* property = odb::dbStringProperty::find(
        db_net, grt::bonnroute::kExtraSpaceProperty);
    if (property == nullptr
        || !grt::bonnroute::decodeExtraSpace(property->getValue(), edges)) {
      continue;
    }
    for (const grt::bonnroute::ExtraSpaceEdge& e : edges) {
      odb::dbTechLayer* layer = tech->findRoutingLayer(e.level);
      if (layer == nullptr || e.space <= 0) {
        continue;
      }
      frLayer* fr_layer = design_->getTech()->getLayer(layer->getName());
      if (fr_layer == nullptr) {
        continue;
      }
      const frLayerNum fr = fr_layer->getLayerNum();
      if (fr < 0 || fr >= static_cast<int>(z_of_fr_.size())
          || z_of_fr_[fr] < 0) {
        continue;
      }
      const int z = z_of_fr_[fr];
      if (!spreadLayer(z)) {
        continue;
      }
      // Both gcells of the edge, mapped through their centres.
      const bool horizontal
          = layer->getDirection() == odb::dbTechLayerDir::HORIZONTAL;
      for (int k = 0; k < 2; k++) {
        const int gx = e.x + (horizontal ? k : 0);
        const int gy = e.y + (horizontal ? 0 : k);
        if (gx < 0 || gy < 0 || gx >= static_cast<int>(grid_x.size())
            || gy >= static_cast<int>(grid_y.size())) {
          continue;
        }
        const int cx = gx + 1 < static_cast<int>(grid_x.size())
                           ? (grid_x[gx] + grid_x[gx + 1]) / 2
                           : grid_x[gx] + tile_dx_ / 2;
        const int cy = gy + 1 < static_cast<int>(grid_y.size())
                           ? (grid_y[gy] + grid_y[gy + 1]) / 2
                           : grid_y[gy] + tile_dy_ / 2;
        net.extra_space_tiles.emplace_back(z, tileX(cx) * tiles_y_ + tileY(cy));
      }
    }
    std::ranges::sort(net.extra_space_tiles);
    net.extra_space_tiles.erase(
        std::ranges::unique(net.extra_space_tiles).begin(),
        net.extra_space_tiles.end());
    if (!net.extra_space_tiles.empty()) {
      extra_space_nets_++;
      extra_space_tiles_ += static_cast<int64_t>(net.extra_space_tiles.size());
    }
  }
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
    if (!spreadLayer(z)) {
      continue;
    }
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

void BonnRouteDR::buildJogTiles()
{
  jog_low_tiles_.clear();
  jog_low_tile_count_ = 0;
  odb::dbGCellGrid* gcells
      = block_ != nullptr ? block_->getGCellGrid() : nullptr;
  if (gcells == nullptr) {
    return;
  }
  jog_low_tiles_.assign(
      graph_.numLayers(),
      std::vector<uint8_t>(static_cast<size_t>(tiles_x_) * tiles_y_, 0));
  for (int z = 0; z < graph_.numLayers(); z++) {
    if (!graph_.layers[z].jogs) {
      continue;
    }
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
            && usage * 100 <= capacity * params_.jog_low_utilization) {
          jog_low_tiles_[z][static_cast<size_t>(tx) * tiles_y_ + ty] = 1;
          jog_low_tile_count_++;
        }
      }
    }
  }
}

bool BonnRouteDR::spreadLayer(const int z) const
{
  // Layers whose wires stay on tracks in the preferred direction have no
  // jogs to leave a free track beside a wire; spreading there only moved
  // wires to other layers (aes_asap7: +7% vias, worse slack, more nets left
  // to FlexDR, no coupling gain).
  return params_.spread_unidirectional || !onTracksOnly(layers_[z].fr);
}

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

namespace {

br::Rect brRect(const odb::Rect& r)
{
  return {r.xMin(), r.yMin(), r.xMax(), r.yMax()};
}

}  // namespace

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
    tl.jog_factor_low = params_.jog_factor_low > 0
                            ? static_cast<Cost>(std::min(params_.jog_factor_low,
                                                         params_.jog_factor))
                                  * graph_.wire_unit
                            : tl.jog_factor;
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
}

br::SpacingTable BonnRouteDR::spacingTable(frLayer* layer) const
{
  // Sec. 3.1: distances are nondecreasing functions of the widths and the
  // common run-length. Sample drt's spacing lookup at every threshold t and
  // t + 1 of its PRL table, two-widths table and range rules, so buckets
  // agree with the lookup for strict and non-strict thresholds.
  std::vector<br::Coord> widths;
  std::vector<br::Coord> prls;
  auto add = [](std::vector<br::Coord>& v, const frCoord t) {
    if (t >= 0) {
      v.push_back(t);
      v.push_back(t + 1);
    }
  };
  frConstraint* con = layer->getMinSpacing();
  if (con != nullptr
      && con->typeId() == frConstraintTypeEnum::frcSpacingTablePrlConstraint) {
    auto table = static_cast<frSpacingTablePrlConstraint*>(con)->getLookupTbl();
    for (const frCoord w : table.getRows()) {
      add(widths, w);
    }
    for (const frCoord p : table.getCols()) {
      add(prls, p);
    }
  }
  odb::dbTechLayer* db_layer = layer->getDbLayer();
  if (db_layer != nullptr && db_layer->hasTwoWidthsSpacingRules()) {
    const uint32_t rows = db_layer->getTwoWidthsSpacingTableNumWidths();
    for (uint32_t row = 0; row < rows; row++) {
      add(widths,
          static_cast<frCoord>(db_layer->getTwoWidthsSpacingTableWidth(row)));
      if (db_layer->getTwoWidthsSpacingTableHasPRL(row)) {
        add(prls,
            static_cast<frCoord>(db_layer->getTwoWidthsSpacingTablePRL(row)));
      }
    }
  }
  if (layer->hasSpacingRangeConstraints()) {
    for (const auto* range : layer->getSpacingRangeConstraints()) {
      add(widths, range->getMinWidth());
      add(widths, range->getMaxWidth());
    }
  }
  return br::SpacingTable::sample(
      widths, prls, [layer](br::Coord a, br::Coord b, br::Coord prl) {
        return layer->getMinSpacingValue(a, b, prl, false);
      });
}

void BonnRouteDR::buildRoutingSpace()
{
  frTechObject* tech = design_->getTech();
  const int num_z = graph_.numLayers();
  // Shape-grid layers: wiring layer z is 2z, the via layer above it 2z + 1.
  std::vector<br::LayerRules> rules;
  for (int z = 0; z < num_z; z++) {
    const LayerInfo& info = layers_[z];
    br::LayerRules wiring;
    wiring.horizontal = info.horizontal;
    wiring.spacing = spacingTable(info.layer);
    // Line-end rule (Sec. 3.1), from drt's summary of the EOL rules.
    const drEolSpacingConstraint& eol = info.layer->getDrEolSpacingConstraint();
    wiring.line_end = {eol.eolWidth, eol.eolSpace, eol.eolWithin};
    rules.push_back(wiring);
    if (z + 1 == num_z) {
      continue;
    }
    br::LayerRules via;
    via.via_layer = true;
    via.horizontal = info.horizontal;
    via.spacing = br::SpacingTable::constant(info.cut_spacing);
    if (z > 0) {
      // Inter-layer via rule with the cut layer below (Sec. 3.1, 3.2).
      frLayer* cut = tech->getLayer(info.fr + 1);
      frLayer* below = tech->getLayer(layers_[z - 1].fr + 1);
      auto inter = [](frLayer* a, frLayer* b) -> frCoord {
        const auto& cons = a->getInterLayerCutSpacingConstraint(false);
        const frLayerNum n = b->getLayerNum();
        if (n < 0 || n >= static_cast<int>(cons.size()) || cons[n] == nullptr) {
          return 0;
        }
        return cons[n]->getCutSpacing();
      };
      via.projection_spacing = std::max(inter(cut, below), inter(below, cut));
    }
    rules.push_back(via);
  }
  // Wire types (Sec. 3.2): the default rule, then the NDRs by their number
  // of nets; the most frequent ones become fast wire types (Sec. 3.6).
  std::vector<frNonDefaultRule*> ndrs;
  std::unordered_map<const frNonDefaultRule*, int> ndr_nets;
  for (const auto& ndr : tech->getNondefaultRules()) {
    ndrs.push_back(ndr.get());
    ndr_nets[ndr.get()] = 0;
  }
  int64_t all_nets = 0;
  for (const auto& net : design_->getTopBlock()->getNets()) {
    all_nets++;
    if (net->getNondefaultRule() != nullptr) {
      ndr_nets[net->getNondefaultRule()]++;
    }
  }
  std::ranges::stable_sort(ndrs, [&](const auto* a, const auto* b) {
    return ndr_nets[a] > ndr_nets[b];
  });
  // Shape classes carry the NDR index in 8 bits; spacing of later NDRs is
  // not modelled (their nets are not routed, see buildNets).
  constexpr size_t kMaxNdrs = 255;
  for (int z = 0; z < num_z; z++) {
    std::vector<br::Coord>& spacing = rules[br::wiringGridLayer(z)].ndr_spacing;
    for (size_t k = 0; k < ndrs.size() && k < kMaxNdrs; k++) {
      spacing.push_back(ndrs[k]->getSpacing(layers_[z].fr / 2 - 1));
    }
  }
  track_grid_layers_ = static_cast<int>(rules.size());
  const int num_fr = static_cast<int>(tech->getLayers().size());
  grid_of_fr_.assign(num_fr, -1);
  for (int z = 0; z < num_z; z++) {
    grid_of_fr_[layers_[z].fr] = br::wiringGridLayer(z);
    if (z + 1 < num_z) {
      grid_of_fr_[layers_[z].fr + 1] = br::viaGridLayer(z);
    }
  }
  // The pin layer and the cut directly below the lowest track layer (Sec.
  // 4.3 off-track pin access and its vias): further shape-grid layers, not
  // part of the fast grid.
  const frLayerNum bottom = layers_.front().fr;
  for (const frLayerNum fr : {bottom - 2, bottom - 1}) {
    if (fr < 0 || fr >= num_fr) {
      continue;
    }
    frLayer* layer = tech->getLayer(fr);
    br::LayerRules low;
    if (fr == bottom - 2 && layer->getType() == odb::dbTechLayerType::ROUTING) {
      low.horizontal = layer->isHorizontal();
      low.spacing = spacingTable(layer);
      const drEolSpacingConstraint& eol = layer->getDrEolSpacingConstraint();
      low.line_end = {eol.eolWidth, eol.eolSpace, eol.eolWithin};
    } else if (fr == bottom - 1
               && layer->getType() == odb::dbTechLayerType::CUT) {
      low.via_layer = true;
      low.horizontal = layers_.front().horizontal;
      frCoord s = layer->hasCutSpacing() ? layer->getCutSpacingValue() : 0;
      if (s <= 0 && layer->getDefaultViaDef() != nullptr) {
        const frVia via(layer->getDefaultViaDef());
        s = via.getCutBBox().minDXDY();
      }
      low.spacing = br::SpacingTable::constant(s);
    } else {
      continue;
    }
    grid_of_fr_[fr] = static_cast<int>(rules.size());
    rules.push_back(low);
  }
  static std::atomic<uint64_t> generations{0};
  class_generation_ = ++generations;
  shape_classes_.setLayers(rules);
  // Wire types (Sec. 3.2): default widths and vias, then the NDRs.
  auto geometry
      = [&](const std::string& name, frNonDefaultRule* ndr, const int index) {
          br::WireTypeGeometry g;
          g.name = name;
          g.ndr = index;
          std::vector<frCoord> widths;
          std::vector<const frViaDef*> vias;
          for (int z = 0; z < num_z; z++) {
            frCoord w = layers_[z].width;
            if (ndr != nullptr) {
              w = std::max(w, ndr->getWidth(layers_[z].fr / 2 - 1));
            }
            g.width.push_back(w);
            widths.push_back(w);
          }
          for (int z = 0; z + 1 < num_z; z++) {
            const frViaDef* def = layers_[z].up_via;
            if (ndr != nullptr
                && ndr->getPrefVia(layers_[z].fr / 2 - 1) != nullptr) {
              def = ndr->getPrefVia(layers_[z].fr / 2 - 1);
            }
            vias.push_back(def);
            const frVia via(def);
            g.via.push_back({true,
                             brRect(via.getLayer1BBox()),
                             brRect(via.getCutBBox()),
                             brRect(via.getLayer2BBox())});
          }
          type_widths_.push_back(std::move(widths));
          type_vias_.push_back(std::move(vias));
          return g;
        };
  wire_types_.clear();
  type_widths_.clear();
  type_vias_.clear();
  type_nets_.clear();
  ndr_type_.clear();
  wire_types_.push_back(
      br::makeWireType(geometry("default", nullptr, 0), shape_classes_));
  type_nets_.push_back(0);
  for (size_t k = 0; k < ndrs.size() && k < kMaxNdrs; k++) {
    const int index = static_cast<int>(k) + 1;
    wire_types_.push_back(br::makeWireType(
        geometry(ndrs[k]->getName(), ndrs[k], index), shape_classes_));
    type_nets_.push_back(ndr_nets[ndrs[k]]);
    ndr_type_[ndrs[k]] = index;
  }
  // Fast wire types: the default rule and the NDRs with at least
  // fast_ndr_permille / 1000 of the nets, most frequent first, up to
  // fast_wire_types.
  fast_types_ = 1;
  while (fast_types_ < static_cast<int>(wire_types_.size())
         && fast_types_ < params_.fast_wire_types && type_nets_[fast_types_] > 0
         && int64_t{type_nets_[fast_types_]} * 1000
                >= int64_t{params_.fast_ndr_permille} * all_nets
         && params_.route_ndr_nets) {
    fast_types_++;
  }
  // Shape grid cells (Sec. 3.3): the largest edge whose diagonal stays
  // below the smallest diff-net spacing, reduced to a divisor of the track
  // pitch when one is at least half as large, and anchored so that tracks
  // lie at cell centres; equal geometry on different tracks then has equal
  // configurations.
  const odb::Rect die = design_->getTopBlock()->getDieBox();
  std::vector<br::ShapeGridLayerGeometry> geometries;
  frCoord max_pitch = 0;
  for (const LayerInfo& info : layers_) {
    max_pitch = std::max(max_pitch, info.pitch);
  }
  int64_t need = 0;
  for (int g = 0; g < shape_classes_.numLayers(); g++) {
    const bool low = g >= track_grid_layers_;
    const int z = low ? 0 : g / 2;
    const bool via_layer = shape_classes_.layer(g).via_layer;
    const br::TrackLayer& tl = graph_.layers[z];
    const frCoord bound = br::maximalCellSize(shape_classes_, g);
    auto cell = [&](frCoord pitch) -> frCoord {
      if (bound <= 0) {
        return std::max<frCoord>(1, pitch / 4);
      }
      for (int k = 1; k <= 64 && pitch > 0; k++) {
        if (pitch % k == 0 && pitch / k <= bound) {
          return pitch / k >= bound / 2 ? pitch / k : bound;
        }
      }
      return bound;
    };
    const frCoord cross_pitch = layers_[z].pitch;
    const frCoord along_pitch
        = z + 1 < num_z ? layers_[z + 1].pitch
                        : (z > 0 ? layers_[z - 1].pitch : cross_pitch);
    br::ShapeGridLayerGeometry geom;
    geom.horizontal = shape_classes_.layer(g).horizontal;
    geom.cross_cell = cell(cross_pitch);
    geom.along_cell = cell(along_pitch);
    const frCoord margin = shape_classes_.ruleReach(g) + 4 * max_pitch;
    const bool h = geom.horizontal;
    const frCoord cross_lo = (h ? die.yMin() : die.xMin()) - margin;
    const frCoord cross_hi = (h ? die.yMax() : die.xMax()) + margin;
    const frCoord along_lo = (h ? die.xMin() : die.yMin()) - margin;
    const frCoord along_hi = (h ? die.xMax() : die.yMax()) + margin;
    auto origin = [](frCoord anchor, frCoord lo, frCoord size) {
      const frCoord start = anchor - size / 2;
      return static_cast<frCoord>(
          start - br::floorDivide(start - lo + size - 1, size) * size);
    };
    // Anchors: the track layer's tracks across, its positions along (a low
    // layer of the other direction swaps them).
    const bool same = geom.horizontal == tl.horizontal;
    const std::vector<br::Coord>& across = same ? tl.tracks : tl.positions;
    const std::vector<br::Coord>& along = same ? tl.positions : tl.tracks;
    geom.cross_origin = origin(
        across.empty() ? cross_lo : across.front(), cross_lo, geom.cross_cell);
    geom.along_origin = origin(
        along.empty() ? along_lo : along.front(), along_lo, geom.along_cell);
    geom.cross_cells = static_cast<int>(
        (int64_t{cross_hi} - geom.cross_origin) / geom.cross_cell + 1);
    geom.along_cells = static_cast<int>(
        (int64_t{along_hi} - geom.along_origin) / geom.along_cell + 1);
    geometries.push_back(geom);
    // Largest distance from a location at which an update or query of this
    // layer can read or write (model, interaction distance, one cell).
    br::Coord model = 0;
    auto extent = [&](const br::WireModel& m) {
      if (m.valid) {
        model = std::max({model,
                          std::abs(m.rect.xlo),
                          std::abs(m.rect.ylo),
                          std::abs(m.rect.xhi),
                          std::abs(m.rect.yhi)});
      }
    };
    for (const br::WireType& w : wire_types_) {
      if (low) {
        model = std::max(model, layers_[0].pitch);
      } else if (!via_layer) {
        extent(w.preferred[z]);
        if (z + 1 < num_z) {
          extent(w.via[z].bottom);
        }
        if (z > 0) {
          extent(w.via[z - 1].top);
        }
      } else {
        extent(w.via[z].cut);
      }
    }
    need = std::max<int64_t>(need,
                             shape_classes_.ruleReach(g) + model
                                 + std::max(geom.cross_cell, geom.along_cell));
  }
  shape_grid_.init(&shape_classes_,
                   geometries,
                   br::ChunkAxis{tile_x0_, tile_dx_, tiles_x_},
                   br::ChunkAxis{tile_y0_, tile_dy_, tiles_y_});
  checker_.init(&shape_grid_);
  // Sec. 3.6: the fast grid holds the most frequent wire types.
  grid_.init(&graph_, fast_types_);
  fast_builder_.init(&graph_, &wire_types_, fast_types_, &checker_, &grid_);
  // Removal recomputations of a batch (flush) run their tile groups on the
  // router's threads; inside a parallel region (partition workers) serially.
  const int flush_threads = std::max(1, router_cfg_->MAX_THREADS);
  fast_builder_.setParallelFor(
      [flush_threads](const int count, const std::function<void(int)>& body) {
        if (omp_in_parallel()) {
          for (int k = 0; k < count; k++) {
            body(k);
          }
          return;
        }
#pragma omp parallel for num_threads(flush_threads) schedule(dynamic, 1)
        for (int k = 0; k < count; k++) {
          body(k);
        }
      });
  // A shape change writes up to need + one pitch beyond the shape and
  // reads up to 2 * need + one pitch (fast-grid recomputation); a query
  // reads up to need. netTileBox adds this many tiles beyond the largest
  // routing area, so threads of disjoint regions (Sec. 5.1) never touch
  // the same tile chunk.
  const int64_t tile = std::max<int64_t>(1, std::min(tile_dx_, tile_dy_));
  interaction_reach_ = need + max_pitch;
  halo_tiles_ = static_cast<int>(
      std::max<int64_t>(1, (2 * need + max_pitch + tile - 1) / tile));
  // Search label budget: a label count, search_label_mb MiB of 32-byte
  // labels (the vertex-search label size when the budget was introduced;
  // the count is kept so that searches stop at the same label).
  label_limit_ = static_cast<size_t>(std::max(0, params_.search_label_mb))
                 * (size_t{1} << 20) / 32;
}

void BonnRouteDR::rebuildFastGrid()
{
  // Whole tile rows per job: a tile's run lists are written by one thread.
  const std::vector<std::array<int, 3>> jobs = fast_builder_.rebuildJobs();
  const int threads = std::max(1, router_cfg_->MAX_THREADS);
  // Vertex data first; jog bits read both tracks of a pair.
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
  for (int k = 0; k < static_cast<int>(jobs.size()); k++) {
    fast_builder_.rebuildVertices(jobs[k][0], jobs[k][1], jobs[k][2]);
  }
  // Jog bits read the next tile row: compute all, then write.
  std::vector<std::vector<br::FastGridBuilder::JogBits>> bits(jobs.size());
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
  for (int k = 0; k < static_cast<int>(jobs.size()); k++) {
    fast_builder_.computeJogBits(jobs[k][0], jobs[k][1], jobs[k][2], bits[k]);
  }
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
  for (int k = 0; k < static_cast<int>(jobs.size()); k++) {
    fast_builder_.writeJogBits(bits[k]);
  }
}

int BonnRouteDR::gridLayer(const frLayerNum layer) const
{
  if (layer < 0 || layer >= static_cast<int>(grid_of_fr_.size())) {
    return -1;
  }
  return grid_of_fr_[layer];
}

int BonnRouteDR::projectionLayer(const int grid_layer) const
{
  const int above = grid_layer + 2;
  if (grid_layer < track_grid_layers_ && grid_layer % 2 == 1
      && above < track_grid_layers_
      && shape_classes_.layer(above).projection_spacing > 0) {
    return above;
  }
  return -1;
}

br::ClassId BonnRouteDR::shapeClass(const RShape& shape,
                                    const int grid_layer) const
{
  br::ShapeRole role = shape.role;
  if (shape_classes_.layer(grid_layer).via_layer
      && role != br::ShapeRole::kCutProjection) {
    role = br::ShapeRole::kCut;
  }
  const br::ShapeClass c = shape_classes_.describe(
      grid_layer, brRect(shape.rect), role, shape.ndr);
  // Per-thread cache; class ids belong to one routing-space build.
  struct Cache
  {
    uint64_t generation = 0;
    std::unordered_map<br::ShapeClass, br::ClassId, br::ShapeClassHash> ids;
  };
  thread_local Cache cache;
  if (cache.generation != class_generation_) {
    cache.ids.clear();
    cache.generation = class_generation_;
  }
  const auto it = cache.ids.find(c);
  if (it != cache.ids.end()) {
    return it->second;
  }
  const br::ClassId id = shape_classes_.intern(c);
  cache.ids.emplace(c, id);
  return id;
}

bool BonnRouteDR::usable(const size_t v,
                         const br::Kind kind,
                         const int net) const
{
  // The net's wire type (Sec. 3.2); its own shapes are ignored.
  const int type = net >= 0 ? nets_[net].wire_type : 0;
  if (type < fast_types_ && grid_.level(v, kind, type) == br::kLegal) {
    return true;
  }
  int z, i, j;
  graph_.decode(v, z, i, j);
  return fast_builder_.exactLevel(z, i, j, type, kind, net) == br::kLegal;
}

namespace {

br::ShapeType routedType(const br::ShapeRole role)
{
  return role == br::ShapeRole::kPreferredWire
                 || role == br::ShapeRole::kNonPreferredWire
             ? br::ShapeType::kWire
             : br::ShapeType::kVia;
}

}  // namespace

void BonnRouteDR::changeShape(const RShape& shape,
                              const br::ShapeType type,
                              const int owner,
                              const bool insert,
                              Worker* local,
                              br::FastGridBuilder::PendingRemovals* pending)
{
  const int g = gridLayer(shape.layer);
  if (g < 0) {
    return;
  }
  br::ShapeGridWriteCache& cache
      = local != nullptr ? local->cache : serial_cache_;
  const br::Rect r = brRect(shape.rect);
  auto apply = [&](const int layer, const br::ClassId cls) {
    const bool ok = insert
                        ? shape_grid_.insert(layer, r, cls, type, owner, cache)
                        : shape_grid_.remove(layer, r, cls, type, owner, cache);
    if (!ok) {
      inconsistent_updates_++;
    }
    // Sec. 3.6: the fast grid follows every change of the shape grid of
    // the track layers and their cuts, incrementally.
    if (layer < track_grid_layers_) {
      const br::PlacedShape placed{r, cls, 3};
      const br::RipupLevel level = shape_grid_.level(type, owner);
      if (insert) {
        fast_builder_.shapeInserted(layer, placed, level);
      } else if (pending != nullptr) {
        fast_builder_.shapeRemoved(layer, placed, level, *pending);
      } else {
        fast_builder_.shapeRemoved(layer, placed, level);
      }
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

void BonnRouteDR::ownShapes(
    const NetData& net,
    std::vector<br::FastGridBuilder::OwnShape>& out) const
{
  out.clear();
  auto add = [&](const RShape& shape, const br::ShapeType type) {
    const int g = gridLayer(shape.layer);
    if (g < 0 || g >= track_grid_layers_) {
      return;
    }
    const br::RipupLevel level = shape_grid_.level(type, net.id);
    out.push_back({g,
                   br::PlacedShape{brRect(shape.rect), shapeClass(shape, g), 3},
                   level});
    const int above = projectionLayer(g);
    if (above >= 0) {
      RShape projection = shape;
      projection.role = br::ShapeRole::kCutProjection;
      out.push_back({above,
                     br::PlacedShape{
                         brRect(shape.rect), shapeClass(projection, above), 3},
                     level});
    }
  };
  for (const RShape& pin : net.pins) {
    add(pin, br::ShapeType::kPin);
  }
  for (const Term& term : net.terms) {
    if (term.reservation_active) {
      for (const RShape& s : term.reservation) {
        add(s, br::ShapeType::kReservation);
      }
    }
  }
  if (net.input_reserved) {
    for (const RShape& s : net.input_reservation) {
      add(s, br::ShapeType::kReservation);
    }
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
    // NDR nets: with their wire type (Sec. 3.2) if route_ndr_nets and the
    // NDR has a wire type.
    const auto ndr = ndr_type_.find(net->getNondefaultRule());
    const bool ndr_routable = params_.route_ndr_nets && ndr != ndr_type_.end();
    if (net->hasNDR() && ndr_routable) {
      data.wire_type = ndr->second;
    }
    data.skip = net->isSpecial() || net->isFake() || net->isFixed()
                || net->hasInitialRouting() || (net->hasNDR() && !ndr_routable)
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
  // Blockages and pins (Sec. 3.3 shape types); pins are owned by their net
  // and removed temporarily while it is connected (Sec. 4.4).
  utl::Timer timer;
  frRegionQuery* rq = design_->getRegionQuery();
  const odb::Rect die = design_->getTopBlock()->getDieBox();
  const int num_fr = static_cast<int>(z_of_fr_.size());
  const int threads = std::max(1, router_cfg_->MAX_THREADS);
  // The fixed shapes of each drt layer in region-query order, with their
  // owner net (-1: blockage).
  struct FixedShape
  {
    odb::Rect rect;
    int owner;
    br::ClassId cls = 0;
    br::ClassId projection_cls = 0;
  };
  std::vector<std::vector<FixedShape>> fixed(num_fr);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
  for (int fr = 0; fr < num_fr; fr++) {
    if (gridLayer(fr) < 0) {
      continue;
    }
    frRegionQuery::Objects<frBlockObject> result;
    rq->query(die, fr, result);
    fixed[fr].reserve(result.size());
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
      fixed[fr].push_back({rect, owner});
    }
  }
  setup_phases_.fixed_query = timer.elapsed();
  // Shape classes in the order of the serial drt-layer loop (a shape, then
  // its cut projection), so that class ids are those of that loop.
  for (int fr = 0; fr < num_fr; fr++) {
    const int g = gridLayer(fr);
    if (g < 0) {
      continue;
    }
    const int above = projectionLayer(g);
    for (FixedShape& f : fixed[fr]) {
      RShape shape{static_cast<frLayerNum>(fr), f.rect};
      f.cls = shapeClass(shape, g);
      if (above >= 0) {
        shape.role = br::ShapeRole::kCutProjection;
        f.projection_cls = shapeClass(shape, above);
      }
    }
  }
  // Shape-grid layer g receives the shapes of its drt layer and the
  // projections of the cut layer below (changeShape), in the order of the
  // drt-layer loop: the cut below comes first. insertFixedShapes inserts
  // them in tile-row bands in parallel, with the result of that order.
  setup_phases_.fixed_classes = timer.elapsed() - setup_phases_.fixed_query;
  std::vector<std::vector<br::FixedShapeInsert>> inserts(
      shape_classes_.numLayers());
  for (int fr = 0; fr < num_fr; fr++) {
    const int g = gridLayer(fr);
    if (g < 0) {
      continue;
    }
    const int above = projectionLayer(g);
    for (const FixedShape& f : fixed[fr]) {
      const br::ShapeType type
          = f.owner >= 0 ? br::ShapeType::kPin : br::ShapeType::kBlockage;
      inserts[g].push_back({brRect(f.rect), f.cls, type, f.owner});
      if (above >= 0) {
        inserts[above].push_back(
            {brRect(f.rect), f.projection_cls, type, f.owner});
      }
    }
  }
  inconsistent_updates_ += br::insertFixedShapes(
      shape_grid_,
      inserts,
      br::ChunkAxis{tile_x0_, tile_dx_, tiles_x_},
      br::ChunkAxis{tile_y0_, tile_dy_, tiles_y_},
      2 * threads,
      [threads](const int count, const std::function<void(int)>& body) {
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
        for (int k = 0; k < count; k++) {
          body(k);
        }
      });
  for (int fr = 0; fr < num_fr; fr++) {
    for (const FixedShape& f : fixed[fr]) {
      if (f.owner >= 0) {
        nets_[f.owner].pins.push_back({static_cast<frLayerNum>(fr), f.rect});
      }
    }
  }
  setup_phases_.fixed_insert = timer.elapsed();
  rebuildFastGrid();
  setup_phases_.fast_grid_rebuild
      = timer.elapsed() - setup_phases_.fixed_insert;
  shape_grid_.compact();
  grid_.compact();
  setup_phases_.compaction = timer.elapsed() - setup_phases_.fixed_insert
                             - setup_phases_.fast_grid_rebuild;
}

bool BonnRouteDR::treeConflicts(const RTree<int>& tree,
                                const RShape& shape,
                                const int net) const
{
  const int g = gridLayer(shape.layer);
  if (g < 0) {
    return false;
  }
  // The pairwise rule of the shape grid on preferred-direction wires.
  const br::Rect r = brRect(shape.rect);
  const br::ClassId cls = shapeClass(shape, g);
  const br::Coord reach = shape_classes_.reach(g);
  odb::Rect query = shape.rect;
  query.bloat(reach, query);
  for (auto it = tree.qbegin(bgi::intersects(query)); it != tree.qend(); ++it) {
    if (it->second == net && net >= 0) {
      continue;
    }
    const br::Rect o = brRect(it->first);
    const br::ClassId other = shapeClass(
        RShape{shape.layer, it->first, br::ShapeRole::kPreferredWire}, g);
    if (br::shapesViolate(shape_classes_,
                          br::PlacedShape{r, cls, 3},
                          br::PlacedShape{o, other, 3})) {
      return true;
    }
  }
  return false;
}

bool BonnRouteDR::shapeConflicts(const RShape& shape,
                                 const int net,
                                 const Worker* worker) const
{
  const int g = gridLayer(shape.layer);
  if (g < 0) {
    return false;
  }
  // Sec. 3.4 / 4.4: everything but the net's own shapes.
  (void) worker;
  const br::Rect r = brRect(shape.rect);
  if (checker_.shapeLevel(g, r, shapeClass(shape, g), 3, net) != br::kLegal) {
    return true;
  }
  const int above = projectionLayer(g);
  if (above >= 0) {
    RShape projection = shape;
    projection.role = br::ShapeRole::kCutProjection;
    return checker_.shapeLevel(above, r, shapeClass(projection, above), 3, net)
           != br::kLegal;
  }
  return false;
}

BonnRouteDR::RShape BonnRouteDR::wireShape(const OutWire& wire) const
{
  // Sec. 3.2: the stick figure plus the wire model of wire type 0.
  const int z = z_of_fr_[wire.layer];
  const bool horizontal = wire.begin.y() == wire.end.y();
  if (z < 0) {
    // A layer below the track graph (pin access): the layer's width.
    const frCoord w = design_->getTech()->getLayer(wire.layer)->getWidth();
    const frCoord half = w / 2;
    return {wire.layer,
            odb::Rect(std::min(wire.begin.x(), wire.end.x()) - half,
                      std::min(wire.begin.y(), wire.end.y()) - half,
                      std::max(wire.begin.x(), wire.end.x()) + w - half,
                      std::max(wire.begin.y(), wire.end.y()) + w - half),
            br::ShapeRole::kFixed};
  }
  const bool preferred
      = wire.begin == wire.end || horizontal == layers_[z].horizontal;
  const br::WireType& type = wire_types_[wire.wire_type];
  const br::Rect& m
      = (preferred ? type.preferred[z] : type.non_preferred[z]).rect;
  odb::Rect r(std::min(wire.begin.x(), wire.end.x()) + m.xlo,
              std::min(wire.begin.y(), wire.end.y()) + m.ylo,
              std::max(wire.begin.x(), wire.end.x()) + m.xhi,
              std::max(wire.begin.y(), wire.end.y()) + m.yhi);
  return {wire.layer,
          r,
          preferred ? br::ShapeRole::kPreferredWire
                    : br::ShapeRole::kNonPreferredWire,
          wire.wire_type};
}

void BonnRouteDR::viaShapes(const OutVia& via,
                            std::vector<RShape>& shapes) const
{
  frVia v(via.def, via.point);
  shapes.push_back({via.def->getLayer1Num(),
                    v.getLayer1BBox(),
                    br::ShapeRole::kViaPad,
                    via.wire_type});
  shapes.push_back(
      {via.def->getCutLayerNum(), v.getCutBBox(), br::ShapeRole::kCut});
  shapes.push_back({via.def->getLayer2Num(),
                    v.getLayer2BBox(),
                    br::ShapeRole::kViaPad,
                    via.wire_type});
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

// The diff-net legality adapter (br::DiffNetLegality) of the off-track
// modules over the shape grid and the distance rule checker (Sec. 3.3, 3.4):
// shapes of the asking net count as removed (Sec. 4.4). Off-track shapes
// are classified by their geometry (no line-end of their own; stored
// shapes keep theirs).
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
    if (dr_.gridLayer(layer) < 0) {
      return true;
    }
    // A box on a track layer is a wire (preferred if it is not longer
    // across the preferred direction), on a cut layer a cut.
    RShape shape{static_cast<frLayerNum>(layer), toRect(box)};
    const int z = dr_.z_of_fr_[layer];
    if (z >= 0) {
      const bool along_x = shape.rect.dx() >= shape.rect.dy();
      shape.role = (shape.rect.dx() == shape.rect.dy()
                    || along_x == dr_.layers_[z].horizontal)
                       ? br::ShapeRole::kPreferredWire
                       : br::ShapeRole::kNonPreferredWire;
    } else if (dr_.design_->getTech()->getLayer(layer)->getType()
               == odb::dbTechLayerType::CUT) {
      shape.role = br::ShapeRole::kCut;
    }
    return allowedShape(shape, net);
  }

  // Exact check of a shape with its role (wire model, via pad, cut).
  bool allowedShape(const RShape& shape, const int net) const
  {
    return dr_.gridLayer(shape.layer) < 0
           || !dr_.shapeConflicts(shape, net, worker_);
  }

  void blockages(const int layer,
                 const br::Box& window,
                 const int net,
                 const br::Box& model,
                 std::vector<br::Box>& out) const override
  {
    const int g = dr_.gridLayer(layer);
    if (g < 0) {
      return;
    }
    // Reference points p with model + p closer than the requirement to a
    // stored piece: the piece bloated by the largest requirement (a box
    // around the l2 region) and its line-end extensions, minus the model.
    const br::ShapeClassTable& classes = dr_.shape_classes_;
    const br::Rect m{model.xlo, model.ylo, model.xhi, model.yhi};
    const br::ClassId mc
        = dr_.shapeClass(RShape{static_cast<frLayerNum>(layer),
                                odb::Rect(m.xlo, m.ylo, m.xhi, m.yhi)},
                         g);
    const br::Coord reach = classes.reach(g);
    thread_local br::GridWindow pieces;
    dr_.shape_grid_.collect(g,
                            br::Rect{window.xlo + m.xlo - reach,
                                     window.ylo + m.ylo - reach,
                                     window.xhi + m.xhi + reach,
                                     window.yhi + m.yhi + reach},
                            pieces);
    auto push = [&](const br::Rect& r, const br::Coord s) {
      out.push_back({r.xlo - s - m.xhi,
                     r.ylo - s - m.yhi,
                     r.xhi + s - m.xlo,
                     r.yhi + s - m.ylo});
    };
    for (const br::GridPiece& p : pieces.pieces) {
      if (p.owner == net && net >= 0 && br::isOwned(p.type)) {
        continue;
      }
      const br::Requirement req = classes.requirement(mc, p.cls);
      if (!req.overlap && !req.line_end
          && std::max({req.facing_x, req.facing_y, req.diagonal}) == 0) {
        continue;
      }
      push(p.rect, std::max({req.facing_x, req.facing_y, req.diagonal}));
      br::Rect ext[2];
      const int n = br::lineEndExtensions(
          classes, br::PlacedShape{p.rect, p.cls, p.line_end_ends}, ext);
      for (int k = 0; k < n; k++) {
        push(ext[k], 0);
      }
    }
  }

  bool conflict(const int layer,
                const br::Box& a,
                const br::Box& b) const override
  {
    const int g = dr_.gridLayer(layer);
    if (g < 0) {
      return false;
    }
    const RShape sa{static_cast<frLayerNum>(layer), toRect(a)};
    const RShape sb{static_cast<frLayerNum>(layer), toRect(b)};
    return br::shapesViolate(
        dr_.shape_classes_,
        br::PlacedShape{brRect(sa.rect), dr_.shapeClass(sa, g), 0},
        br::PlacedShape{brRect(sb.rect), dr_.shapeClass(sb, g), 0});
  }

 private:
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
  empty_no_pin += o.empty_no_pin;
  empty_no_targets += o.empty_no_targets;
  empty_no_path += o.empty_no_path;
  empty_rejected += o.empty_rejected;
  macro_pins += o.macro_pins;
  member_regenerated += o.member_regenerated;
  target_seconds += o.target_seconds;
  blockage_seconds += o.blockage_seconds;
  grid_seconds += o.grid_seconds;
  check_seconds += o.check_seconds;
  selection_seconds += o.selection_seconds;
  for (const std::string& e : o.empty_examples) {
    if (empty_examples.size() < 10) {
      empty_examples.push_back(e);
    }
  }
}

bool BonnRouteDR::onTracksOnly(const frLayerNum layer) const
{
  if (layer < 0 || layer >= static_cast<int>(z_of_fr_.size())) {
    return false;
  }
  frLayer* l = design_->getTech()->getLayer(layer);
  return l->getType() == odb::dbTechLayerType::ROUTING
         && (l->isUnidirectional()
             || l->getLef58RightWayOnGridOnlyConstraint() != nullptr);
}

bool BonnRouteDR::wireOnTracks(const OutWire& wire) const
{
  if (!onTracksOnly(wire.layer) || wire.begin == wire.end) {
    return true;
  }
  const int z = z_of_fr_[wire.layer];
  if (z < 0) {
    return false;
  }
  const br::TrackLayer& l = graph_.layers[z];
  const bool horizontal = wire.begin.y() == wire.end.y();
  if (horizontal != l.horizontal) {
    return false;
  }
  const frCoord cross = horizontal ? wire.begin.y() : wire.begin.x();
  return std::ranges::binary_search(l.tracks, cross);
}

void BonnRouteDR::accessSources(const Term& term,
                                std::vector<AccessSource>& out,
                                const bool preferred_only) const
{
  out.clear();
  auto add = [&](const frAccessPoint* ap,
                 const odb::Point& p,
                 const bool preferred = true) {
    const frViaDef* via = nullptr;
    if (ap->hasAccess(frDirEnum::U) && ap->hasViaDef()) {
      via = ap->getViaDefs(1).front();
    }
    out.push_back({ap->getLayerNum(), p, via, preferred});
  };
  if (term.inst_term != nullptr) {
    frInst* inst = term.inst_term->getInst();
    const int pa_idx = inst->getPinAccessIdx();
    const odb::dbTransform xform = inst->getNoRotationTransform();
    const auto& prefs = term.inst_term->getAccessPoints();
    int pin_idx = 0;
    for (const auto& pin : term.inst_term->getTerm()->getPins()) {
      const frAccessPoint* pref
          = pin_idx < static_cast<int>(prefs.size()) ? prefs[pin_idx] : nullptr;
      pin_idx++;
      if (!pin->hasPinAccess() || pa_idx < 0
          || pa_idx >= pin->getNumPinAccess()) {
        continue;
      }
      for (const auto& ap : pin->getPinAccess(pa_idx)->getAccessPoints()) {
        if (preferred_only && pref != nullptr && ap.get() != pref) {
          continue;
        }
        odb::Point p = ap->getPoint();
        xform.apply(p);
        add(ap.get(), p, pref == nullptr || ap.get() == pref);
      }
    }
  } else if (term.bterm != nullptr) {
    for (const auto& pin : term.bterm->getPins()) {
      if (!pin->hasPinAccess()) {
        continue;
      }
      for (const auto& ap : pin->getPinAccess(0)->getAccessPoints()) {
        add(ap.get(), ap->getPoint());
      }
    }
  }
}

void BonnRouteDR::accessPointsOf(const frBlockObject* pin,
                                 std::vector<br::RouteNode>& out) const
{
  out.clear();
  Term term;
  if (pin->typeId() == frcInstTerm) {
    term.inst_term
        = const_cast<frInstTerm*>(static_cast<const frInstTerm*>(pin));
  } else if (pin->typeId() == frcBTerm) {
    term.bterm = const_cast<frBTerm*>(static_cast<const frBTerm*>(pin));
  } else {
    return;
  }
  std::vector<AccessSource> sources;
  accessSources(term, sources);
  for (const AccessSource& s : sources) {
    out.push_back({s.layer, toGrid(s.point)});
  }
}

namespace {

// The net's pins (instance and top-level terms) whose shapes on `layer`
// contain p, with those shapes.
template <typename F>
void forEachNetPinAt(frRegionQuery* rq,
                     const frNet* net,
                     const frLayerNum layer,
                     const odb::Point& p,
                     F&& fn)
{
  frRegionQuery::Objects<frBlockObject> result;
  rq->query(odb::Rect(p, p), layer, result);
  for (const auto& [rect, obj] : result) {
    const frNet* owner = nullptr;
    if (obj->typeId() == frcInstTerm) {
      owner = static_cast<const frInstTerm*>(obj)->getNet();
    } else if (obj->typeId() == frcBTerm) {
      owner = static_cast<const frBTerm*>(obj)->getNet();
    }
    if (owner == net && rect.intersects(p)) {
      fn(obj, rect);
    }
  }
}

}  // namespace

int64_t BonnRouteDR::pinContactViolations(const NetData& net,
                                          const std::vector<OutWire>& wires,
                                          const std::vector<OutVia>& vias,
                                          const bool count_vias) const
{
  // FlexDR makes a pin a term of a worker (initNets_searchRepair_pin2epMap,
  // and initNets_segmentTerms for border segments) only at truncated
  // path-segment ends and at via pads marked bottom/top connected, and then
  // needs one of the pin's access points inside the worker's box
  // (initNet_term_helper, DRT-1231). The connected flags are set only by
  // FlexDR's maze on vias it creates (FlexDR_maze.cpp, at access points or
  // through checkViaConnectivityToAP, which requires an access point);
  // vias of the initial routing keep the default false (frVia). So only
  // the truncated ends, the planar pin points, are a risk; via pads on pins
  // are counted separately with `count_vias` for the log.
  std::vector<br::RouteNode> contacts;
  if (!count_vias) {
    for (const OutWire& w : wires) {
      for (const auto& [layer, p] : net.pin_points) {
        if (layer == w.layer && (p == w.begin || p == w.end)) {
          contacts.push_back({layer, toGrid(p)});
        }
      }
    }
  } else {
    for (const OutVia& v : vias) {
      contacts.push_back({v.def->getLayer1Num(), toGrid(v.point)});
      contacts.push_back({v.def->getLayer2Num(), toGrid(v.point)});
    }
  }
  std::ranges::sort(contacts,
                    [](const br::RouteNode& a, const br::RouteNode& b) {
                      return std::tie(a.layer, a.p.x, a.p.y)
                             < std::tie(b.layer, b.p.x, b.p.y);
                    });
  contacts.erase(std::ranges::unique(contacts).begin(), contacts.end());
  frRegionQuery* rq = design_->getRegionQuery();
  int64_t violations = 0;
  for (const br::RouteNode& c : contacts) {
    std::vector<br::ContactPin> pins;
    std::vector<const frBlockObject*> objs;
    forEachNetPinAt(rq,
                    net.net,
                    c.layer,
                    toPoint(c.p),
                    [&](const frBlockObject* obj, const odb::Rect& rect) {
                      auto it = std::ranges::find(objs, obj);
                      if (it == objs.end()) {
                        objs.push_back(obj);
                        pins.emplace_back();
                        accessPointsOf(obj, pins.back().access_points);
                        it = objs.end() - 1;
                      }
                      pins[it - objs.begin()].shapes.emplace_back(c.layer,
                                                                  toBox(rect));
                    });
    violations += br::pinContactViolations({c}, pins);
  }
  return violations;
}

bool BonnRouteDR::atAccessPoint(const NetData& net,
                                const frLayerNum layer,
                                const odb::Point& p) const
{
  bool found = false;
  forEachNetPinAt(design_->getRegionQuery(),
                  net.net,
                  layer,
                  p,
                  [&](const frBlockObject* obj, const odb::Rect&) {
                    std::vector<br::RouteNode> aps;
                    accessPointsOf(obj, aps);
                    found = found
                            || std::ranges::find(
                                   aps, br::RouteNode{layer, toGrid(p)})
                                   != aps.end();
                  });
  return found;
}

bool BonnRouteDR::nodeAllowed(const int net,
                              const frLayerNum layer,
                              const odb::Point& p) const
{
  bool ok = true;
  forEachNetPinAt(design_->getRegionQuery(),
                  nets_[net].net,
                  layer,
                  p,
                  [&](const frBlockObject* obj, const odb::Rect&) {
                    std::vector<br::RouteNode> aps;
                    accessPointsOf(obj, aps);
                    ok = ok
                         && std::ranges::find(aps,
                                              br::RouteNode{layer, toGrid(p)})
                                != aps.end();
                  });
  return ok;
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
    // tau on the manufacturing grid; without a mapped rule the grid itself.
    const frCoord grid
        = std::max<frCoord>(1, design_->getTech()->getManufacturingGrid());
    br::MinimumSegmentLength m = br::minimumSegmentLength(
        rules, params_.offtrack_notch_tau, params_.offtrack_tau, grid);
    m.tau = (m.tau + grid - 1) / grid * grid;
    same_net_rules_[z] = rules;
    tau_[z] = m.tau;
    br::CleanupLayer& cl = cleanup_layers_[info.fr];
    cl.tau = m.tau;
    cl.half_width = info.width / 2;
    cl.rules = rules;
    cl.horizontal = info.horizontal;
    cl.grid = grid;
    cl.on_tracks = onTracksOnly(info.fr);
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
  if (!usable(option.vertex, br::kWire, net)) {
    return false;
  }
  const Legality legality(*this, worker);
  std::vector<RShape> shapes;
  for (const OutWire& wire : option.wires) {
    shapes.push_back(wireShape(wire));
  }
  for (size_t k = 0; k < option.vias.size(); k++) {
    std::vector<RShape> via;
    viaShapes(option.vias[k], via);
    if (static_cast<int>(k) == option.pin_via) {
      // Pin access checked the pin-layer pad and the cut against the fixed
      // surroundings; they still meet other nets' routing and
      // reservations, and the top pad meets everything.
      if (routedConflicts(via[0], net) || routedConflicts(via[1], net)) {
        return false;
      }
      shapes.push_back(via.back());
    } else {
      shapes.insert(shapes.end(), via.begin(), via.end());
    }
  }
  return std::ranges::all_of(
      shapes, [&](const RShape& s) { return legality.allowedShape(s, net); });
}

bool BonnRouteDR::routedConflicts(const RShape& shape, const int net) const
{
  const int g = gridLayer(shape.layer);
  if (g < 0) {
    return false;
  }
  const br::ShapeClassTable& classes = shape_classes_;
  const br::Coord reach = classes.reach(g);
  const br::Rect r = brRect(shape.rect);
  thread_local br::GridWindow pieces;
  shape_grid_.collect(g, r.bloated(reach, reach), pieces);
  const br::PlacedShape placed{r, shapeClass(shape, g), 3};
  for (const br::GridPiece& p : pieces.pieces) {
    if (p.type == br::ShapeType::kBlockage || p.type == br::ShapeType::kPin
        || (p.owner == net && net >= 0)) {
      continue;
    }
    if (br::shapesViolate(
            classes, placed, br::PlacedShape{p.rect, p.cls, p.line_end_ends})) {
      return true;
    }
  }
  return false;
}

void BonnRouteDR::offTrackCatalogue(const int net,
                                    const Term& term,
                                    const Worker* worker,
                                    std::vector<AccessOption>& out,
                                    AccessStats& stats,
                                    const int extra_radius) const
{
  out.clear();
  stats.pins++;
  utl::Timer step_timer;
  auto empty = [&](int64_t& reason, const char* why) {
    reason++;
    stats.empty++;
    if (stats.empty_examples.size() < 10) {
      const std::string name = term.inst_term != nullptr
                                   ? term.inst_term->getName()
                               : term.bterm != nullptr ? term.bterm->getName()
                                                       : std::string("?");
      stats.empty_examples.push_back(name + " (" + why + ")");
    }
  };
  frLayerNum pin_layer = -1;
  std::vector<odb::Rect> pin;
  const bool from_aps = params_.offtrack_from_access_points;
  // Macro pins (e.g. fakeram45 metal3 pins) abut the macro's obstruction
  // on their layer and the next one: every path leaving the pin starts
  // inside the obstruction's spacing halo, so no tau-feasible path exists.
  // OpenROAD's pin access treats that abutment as intended; these pins keep
  // the access stubs unless offtrack_macro_pins is set.
  if (!params_.offtrack_macro_pins && term.inst_term != nullptr
      && term.inst_term->getInst()->getMaster()->getMasterType().isBlock()) {
    stats.pins--;
    stats.macro_pins++;
    return;
  }
  const int num_fr = static_cast<int>(z_of_fr_.size());
  const int num_z = static_cast<int>(layers_.size());
  // Sources (Sec. 4.3 paths "connecting it to on-track points"): the
  // pin's access points, or the whole pin with `from_aps` false.
  std::vector<AccessSource> aps;
  if (from_aps) {
    std::vector<AccessSource> all;
    accessSources(term, all);
    // Access points off the manufacturing grid would put the via and the
    // lines of the blockage grid off it.
    const frCoord mfg_grid
        = std::max<frCoord>(1, design_->getTech()->getManufacturingGrid());
    std::erase_if(all, [&](const AccessSource& a) {
      return a.point.x() % mfg_grid != 0 || a.point.y() % mfg_grid != 0;
    });
    const frLayerNum bottom = layers_.front().fr;
    for (const AccessSource& a : all) {
      const bool usable
          = (a.layer >= 0 && a.layer < num_fr && z_of_fr_[a.layer] >= 0)
            || (a.layer == bottom - 2 && a.via != nullptr
                && a.via->getLayer2Num() == bottom);
      if (usable && (pin_layer < 0 || a.layer < pin_layer)) {
        pin_layer = a.layer;
      }
    }
    for (const AccessSource& a : all) {
      if (a.layer == pin_layer
          && (z_of_fr_[a.layer] >= 0 || a.via != nullptr)) {
        aps.push_back(a);
        pin.emplace_back(a.point, a.point);
      }
    }
    if (aps.empty()) {
      empty(stats.empty_no_pin, "no usable access point");
      return;
    }
  } else if (!pinShapes(term, pin_layer, pin)) {
    empty(stats.empty_no_pin, "no pin shape");
    return;
  }
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
  // A pin layer without wires is used only at the access points, each with
  // its own via (points_only in the blockage grid).
  const bool points_only = from_aps && pin_z < 0;
  if (points_only) {
    via = nullptr;
  }
  const bool ap_vias
      = from_aps && upper >= 0 && std::ranges::any_of(aps, [&](const auto& a) {
          return a.via != nullptr && a.via->getLayer2Num() == upper;
        });
  // Search layers: the pin layer and the next routing layer (Sec. 4.3:
  // access paths "connecting it to on-track points within a small radius").
  std::vector<frLayerNum> frs{pin_layer};
  if (via != nullptr || ap_vias) {
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
    const frCoord r
        = (params_.offtrack_radius + extra_radius) * layers_[z].pitch;
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
        if (!usable(v, br::kWire, net)) {
          continue;
        }
        const auto [x, y] = l.point(i, j);
        targets.push_back({static_cast<int>(s), z, v, odb::Point(x, y)});
      }
    }
  }
  stats.target_seconds += step_timer.elapsed();
  if (targets.empty()) {
    empty(stats.empty_no_targets, "no usable on-track point");
    return;
  }
  step_timer.reset();
  odb::Rect area = pin_box;
  for (const Target& t : targets) {
    area.merge(odb::Rect(t.p, t.p));
  }
  const frCoord margin = 2 * max_tau + max_pitch;
  // Blockage-grid lines on the manufacturing grid: obstacles grow outward
  // to it (conservative), so every coordinate of Algorithm 3 (lines plus
  // multiples of tau, itself on the grid) is on it.
  const frCoord mfg
      = std::max<frCoord>(1, design_->getTech()->getManufacturingGrid());
  auto snap_out = [mfg](std::vector<br::Box>& boxes, const size_t from) {
    auto down = [mfg](const br::GeoCoord c) {
      return static_cast<br::GeoCoord>(c >= 0 ? c / mfg * mfg
                                              : -((-c + mfg - 1) / mfg) * mfg);
    };
    for (size_t k = from; k < boxes.size(); k++) {
      br::Box& b = boxes[k];
      b = {down(b.xlo), down(b.ylo), -down(-b.xhi), -down(-b.yhi)};
    }
  };
  std::vector<br::Box> window_box{toBox(area).bloated(margin)};
  snap_out(window_box, 0);
  const br::Box window = window_box.front();
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
      snap_out(layer.obstacles, 0);
      // The on-track cost model: non-preferred wire costs the jog factor.
      // (Theorem 3.2 is checked with such axis weights as well.)
      const int64_t jog = std::max(1, params_.jog_factor);
      layer.cost_x = layers_[z].horizontal ? 1 : jog;
      layer.cost_y = layers_[z].horizontal ? jog : 1;
    }
    grid_in.layers.push_back({z >= 0 ? tau_[z] : 1,
                              layer.obstacles,
                              points_only && fr == pin_layer});
    req.layers.push_back(std::move(layer));
  }
  std::vector<br::Box> via_blocks;
  std::array<RShape, 3> via_model{};
  if (via != nullptr) {
    const frVia v0(via, odb::Point(0, 0));
    via_model = {
        RShape{via->getLayer1Num(), v0.getLayer1BBox(), br::ShapeRole::kViaPad},
        RShape{via->getCutLayerNum(), v0.getCutBBox(), br::ShapeRole::kCut},
        RShape{
            via->getLayer2Num(), v0.getLayer2BBox(), br::ShapeRole::kViaPad}};
    for (const RShape& m : via_model) {
      legality.blockages(m.layer, window, net, toBox(m.rect), via_blocks);
    }
    snap_out(via_blocks, 0);
    grid_in.via_blockages.push_back(via_blocks);
  }
  for (const Target& t : targets) {
    grid_in.points.push_back(toGrid(t.p));
  }
  if (from_aps) {
    for (const AccessSource& a : aps) {
      grid_in.points.push_back(toGrid(a.point));
    }
  } else {
    for (const odb::Rect& r : pin) {
      grid_in.borders.push_back(toBox(r));
    }
  }
  grid_in.borders.push_back(window);
  stats.blockage_seconds += step_timer.elapsed();
  step_timer.reset();
  const std::vector<br::BlockageGridLayer> grid
      = br::buildBlockageGrid(grid_in);
  stats.grid_seconds += step_timer.elapsed();
  for (size_t s = 0; s < frs.size(); s++) {
    req.layers[s].grid = grid[s];
  }
  req.window = window;
  // Access points pin access did not choose start with the stubs' penalty
  // (pin_cost pitches), so their paths enter the catalogue only where they
  // reach targets the chosen one does not.
  const int64_t ap_penalty = static_cast<int64_t>(params_.pin_cost)
                             * std::max<frCoord>(max_pitch, 1);
  if (from_aps) {
    for (const AccessSource& a : aps) {
      const int64_t c = params_.offtrack_preferred_access_points && !a.preferred
                            ? ap_penalty
                            : 0;
      req.sources.push_back({0, toBox(odb::Rect(a.point, a.point)), c});
    }
  } else {
    for (const odb::Rect& r : pin) {
      req.sources.push_back({0, toBox(r)});
    }
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
  // On layers whose wires must stay on tracks in the preferred direction
  // (Sec. 3.8: "removing edges that induce a violation"), only planar arcs
  // along a track remain.
  for (size_t s = 0; s < frs.size(); s++) {
    if (onTracksOnly(frs[s]) && z_of_fr_[frs[s]] >= 0) {
      const br::TrackLayer& l = graph_.layers[z_of_fr_[frs[s]]];
      req.layers[s].tracks = &l.tracks;
      req.layers[s].horizontal = l.horizontal;
    }
  }
  req.via_allowed = [&](int lower, const br::GridPoint& p) {
    if (lower == 0) {
      // An access point's own via: pin access checked its pin-layer pad
      // and cut; its top pad must meet the routing space.
      for (const AccessSource& a : aps) {
        if (a.via != nullptr && a.via->getLayer2Num() == upper
            && toGrid(a.point) == p) {
          std::vector<RShape> shapes;
          viaShapes({a.via, a.point}, shapes);
          return !routedConflicts(shapes[0], net)
                 && !routedConflicts(shapes[1], net)
                 && legality.allowedShape(shapes[2], net);
        }
      }
    }
    if (via == nullptr) {
      return false;
    }
    for (const br::Box& b : via_blocks) {
      if (b.containsOpen(p.x, p.y)) {
        return false;
      }
    }
    for (RShape m : via_model) {
      m.rect.moveDelta(p.x, p.y);
      if (!legality.allowedShape(m, net)) {
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
  step_timer.reset();
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
        const frViaDef* def = viaDefAbove(lower);
        if (lower == pin_layer) {
          for (const AccessSource& ap : aps) {
            if (ap.via != nullptr && toGrid(ap.point) == a.p) {
              def = ap.via;
              opt.pin_via = static_cast<int>(opt.vias.size());
              break;
            }
          }
        }
        opt.vias.push_back({def, toPoint(a.p)});
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
    if (from_aps && params_.offtrack_preferred_access_points
        && !opt.route.empty()) {
      for (const AccessSource& a : aps) {
        if (toGrid(a.point) == opt.route.front().p && !a.preferred) {
          cost += static_cast<Cost>(params_.pin_cost) * layers_[0].pitch
                  * unit();
          break;
        }
      }
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
  stats.check_seconds += step_timer.elapsed();
  stats.paths += static_cast<int64_t>(out.size());
  if (out.empty()) {
    if (res.paths.empty()) {
      empty(stats.empty_no_path, "no tau-feasible path");
    } else {
      empty(stats.empty_rejected, "all paths rejected");
    }
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
  if (j + 1 < l.numPositions() && usable(v, br::kWireEdge, net)
      && usable(graph_.vertex(z, i, j + 1), br::kWire, net)) {
    feasible++;
  }
  if (j > 0) {
    const size_t p = graph_.vertex(z, i, j - 1);
    if (usable(p, br::kWireEdge, net) && usable(p, br::kWire, net)) {
      feasible++;
    }
  }
  if (z + 1 < graph_.numLayers() && l.up_track[j] >= 0 && l.up_pos[i] >= 0) {
    const size_t u = graph_.vertex(z + 1, l.up_track[j], l.up_pos[i]);
    if (usable(v, br::kViaUpPad, net) && usable(v, br::kCutUp, net)
        && usable(u, br::kViaDownPad, net)) {
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
  utl::Timer selection_timer;
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
  // Candidate shapes classified once, as Legality::conflict classifies
  // them (role kFixed, no line-end ends), with their grid layer and the
  // layer's rule reach: shapes farther apart than the reach in x or y
  // never violate a rule of the layer (ShapeClassTable::ruleReach).
  struct ConflictShape
  {
    frLayerNum layer;
    int grid_layer;
    br::Coord reach;
    br::PlacedShape placed;
  };
  std::vector<std::vector<std::vector<ConflictShape>>> classified(n);
  for (int p = 0; p < n; p++) {
    classified[p].resize(shapes[p].size());
    for (size_t c = 0; c < shapes[p].size(); c++) {
      for (const RShape& s : shapes[p][c]) {
        const int g = gridLayer(s.layer);
        if (g < 0) {
          continue;
        }
        const RShape fixed{s.layer, s.rect};
        classified[p][c].push_back(
            {s.layer,
             g,
             shape_classes_.reach(g),
             br::PlacedShape{brRect(s.rect), shapeClass(fixed, g), 0}});
      }
    }
  }
  in.conflict = [&](int pa, int ca, int pb, int cb) {
    for (const ConflictShape& a : classified[pa][ca]) {
      for (const ConflictShape& b : classified[pb][cb]) {
        if (a.layer != b.layer) {
          continue;
        }
        const br::Rect& ra = a.placed.rect;
        const br::Rect& rb = b.placed.rect;
        const int64_t gx = std::max<int64_t>(
            {0, int64_t{rb.xlo} - ra.xhi, int64_t{ra.xlo} - rb.xhi});
        const int64_t gy = std::max<int64_t>(
            {0, int64_t{rb.ylo} - ra.yhi, int64_t{ra.ylo} - rb.yhi});
        if (gx > a.reach || gy > a.reach) {
          continue;
        }
        if (br::shapesViolate(shape_classes_, a.placed, b.placed)) {
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
  stats.selection_seconds += selection_timer.elapsed();
  return r.choice;
}

void BonnRouteDR::buildOffTrackAccess(const std::vector<char>& pre_skip)
{
  utl::Timer timer;
  frBlock* block = design_->getTopBlock();
  frRegionQuery* rq = design_->getRegionQuery();
  AccessStats& stats = access_stats_;
  AccessPhases& phases = access_phases_;
  utl::Timer phase_timer;
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
  phases.pin_map = phase_timer.elapsed();
  phase_timer.reset();
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
  // Sec. 4.3: "Expanding the bounding area of an instance ... and taking
  // all shapes influencing this area and the track coordinates into
  // account". The area is the instance box plus access_class_margin
  // pitches; access paths reaching beyond it are re-checked at each member
  // (mapped paths illegal there are dropped and the pin's catalogue is
  // generated again if none remains).
  const frCoord halo = params_.access_class_margin * max_pitch;
  const frLayerNum bottom = layers_.front().fr;
  // Layers that influence access paths: up to the cut above the lowest
  // track layer (the search uses the pin layer and the next routing layer).
  const frLayerNum sig_top = layers_.front().fr + 1;
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
      if (params_.offtrack_from_access_points) {
        // The access points the catalogue paths start at (label -1 - p).
        std::vector<AccessSource> aps;
        accessSources(nets_[ia.terms[p].first].terms[ia.terms[p].second], aps);
        for (const AccessSource& a : aps) {
          sig.addShape(a.layer,
                       {a.point.x(), a.point.y(), a.point.x(), a.point.y()},
                       (a.preferred ? -1 : -1001) - static_cast<int>(p));
        }
      }
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
    // Tracks of the search layers and the layer above (its tracks are the
    // positions of the lowest track layer); higher layers' track phases
    // relative to the cell do not influence access and only split classes.
    for (int z = 0; z < std::min(num_z, 2); z++) {
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
  phases.signatures = phase_timer.elapsed();
  phase_timer.reset();
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
  phases.classification = phase_timer.elapsed();
  phase_timer.reset();
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
  AccessStats class_total;
  for (const AccessStats& s : class_stats) {
    stats.add(s);
    class_total.add(s);
  }
  phases.class_catalogues = phase_timer.elapsed();
  phase_timer.reset();
  utl::Timer step_timer;
  // Per instance: map the class solution, check it against the routing
  // space (reservations of earlier instances), re-solve if needed, and
  // reserve the primary paths.
  auto map_option = [&](const AccessOption& in,
                        const br::Isometry& m,
                        const std::vector<AccessSource>& member_aps,
                        AccessOption& out) {
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
    if (!params_.offtrack_from_access_points || out.route.empty()) {
      return true;
    }
    // The mapped start must be one of the member's access points (they are
    // in the signature; this also checks the via), with its via.
    const br::RouteNode& start = out.route.front();
    for (const AccessSource& a : member_aps) {
      if (a.layer == start.layer && toGrid(a.point) == start.p) {
        if (out.pin_via >= 0) {
          if (a.via == nullptr) {
            return false;
          }
          out.vias[out.pin_via].def = a.via;
        }
        return true;
      }
    }
    return false;
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
    // One net per pin, complete before mapping: solve() and the legality
    // checks below index nets by pin, also after a failed mapping.
    std::vector<int> nets;
    nets.reserve(n);
    for (const auto& [net, t] : ia.terms) {
      nets.push_back(net);
    }
    bool mapped = true;
    step_timer.reset();
    for (int p = 0; p < n; p++) {
      const auto [net, t] = ia.terms[p];
      if (k == rep) {
        cats[p] = cd.catalogues[p];
        continue;
      }
      std::vector<AccessSource> member_aps;
      accessSources(nets_[net].terms[t], member_aps);
      for (const AccessOption& opt : cd.catalogues[p]) {
        AccessOption mo;
        if (!map_option(opt, m, member_aps, mo)) {
          mapped = false;
          break;
        }
        cats[p].push_back(std::move(mo));
      }
      if (!mapped) {
        break;
      }
    }
    phases.member_mapping += step_timer.elapsed();
    if (!mapped) {
      // Equal signatures with unequal track positions outside the window or
      // a mapped start that is not an access point (other via); generate
      // this instance's catalogues directly.
      stats.mapping_failed++;
      step_timer.reset();
      for (int p = 0; p < n; p++) {
        const auto [net, t] = ia.terms[p];
        offTrackCatalogue(net, nets_[net].terms[t], nullptr, cats[p], stats);
      }
      phases.member_catalogues += step_timer.elapsed();
    }
    std::vector<int> choice = mapped ? cd.choice : std::vector<int>(n, -1);
    bool reuse = mapped;
    if (mapped && k != rep) {
      // Shapes outside the signature area can make mapped paths illegal at
      // this member; a pin left without a legal path gets its own
      // catalogue.
      for (int p = 0; p < n; p++) {
        step_timer.reset();
        const bool any_legal
            = cats[p].empty()
              || std::ranges::any_of(cats[p], [&](const AccessOption& o) {
                   return optionLegal(o, nets[p], nullptr);
                 });
        phases.member_legality += step_timer.elapsed();
        if (any_legal) {
          continue;
        }
        step_timer.reset();
        const auto [net, t] = ia.terms[p];
        offTrackCatalogue(net, nets_[net].terms[t], nullptr, cats[p], stats);
        stats.member_regenerated++;
        reuse = false;
        phases.member_catalogues += step_timer.elapsed();
      }
    }
    step_timer.reset();
    for (int p = 0; p < n && reuse; p++) {
      reuse = choice[p] < 0 ? cats[p].empty()
                            : optionLegal(cats[p][choice[p]], nets[p], nullptr);
    }
    phases.member_legality += step_timer.elapsed();
    if (reuse) {
      stats.class_reused += k != rep;
    } else {
      step_timer.reset();
      stats.reselected++;
      choice = solve(cats, nets);
      phases.member_solve += step_timer.elapsed();
    }
    step_timer.reset();
    for (int p = 0; p < n; p++) {
      const auto [net_id, t] = ia.terms[p];
      NetData& net = nets_[net_id];
      Term& term = net.terms[t];
      if (cats[p].empty()) {
        reserve_fallback(net, term);
        continue;
      }
      keepStubs(term, cats[p]);
      if (choice[p] >= 0) {
        reserve(net, term, choice[p]);
      }
    }
    phases.member_reservations += step_timer.elapsed();
  }
  phases.members = phase_timer.elapsed();
  phase_timer.reset();
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
    keepStubs(term, cats[0]);
    if (choice[0] >= 0) {
      reserve(net, term, choice[0]);
    }
  }
  phases.top_level_pins = phase_timer.elapsed();
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
        "{} primary paths; {:.2f}s. Empty catalogues: {} without usable "
        "access point or pin shape, {} without usable on-track point, {} "
        "without path, {} with all paths rejected, {} macro pins with "
        "access stubs; {} member pins regenerated{}{}.",
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
        stats.seconds,
        stats.empty_no_pin,
        stats.empty_no_targets,
        stats.empty_no_path,
        stats.empty_rejected,
        stats.macro_pins,
        stats.member_regenerated,
        stats.empty_examples.empty() ? "" : "; e.g. ",
        [&] {
          std::string list;
          for (const std::string& e : stats.empty_examples) {
            list += (list.empty() ? "" : ", ") + e;
          }
          return list;
        }());
    logger_->report(
        "BonnRoute off-track access phases: pin map {:.2f}s, signatures "
        "{:.2f}s, classification {:.2f}s, class catalogues {:.2f}s wall "
        "({:.2f}s thread: targets {:.2f}s, blockages {:.2f}s, blockage grid "
        "{:.2f}s, search {:.2f}s, path checks {:.2f}s, selection {:.2f}s), "
        "members {:.2f}s (mapping {:.2f}s, legality {:.2f}s, own catalogues "
        "{:.2f}s, re-solve {:.2f}s, reservations {:.2f}s), top-level pins "
        "{:.2f}s.",
        phases.pin_map,
        phases.signatures,
        phases.classification,
        phases.class_catalogues,
        class_total.target_seconds + class_total.blockage_seconds
            + class_total.grid_seconds + class_total.search_seconds
            + class_total.check_seconds + class_total.selection_seconds,
        class_total.target_seconds,
        class_total.blockage_seconds,
        class_total.grid_seconds,
        class_total.search_seconds,
        class_total.check_seconds,
        class_total.selection_seconds,
        phases.members,
        phases.member_mapping,
        phases.member_legality,
        phases.member_catalogues,
        phases.member_solve,
        phases.member_reservations,
        phases.top_level_pins);
  }
}

void BonnRouteDR::keepStubs(Term& term, std::vector<AccessOption>& catalogue)
{
  // The catalogue first (primary indices refer to it), then the access
  // stubs built from the same access points: a pin whose catalogue paths
  // are all blocked later still has the stubs' alternatives on other
  // tracks and positions (on unidirectional layers the catalogue of an
  // access point is confined to one track).
  std::vector<AccessOption> stubs = std::move(term.options);
  term.options = std::move(catalogue);
  for (AccessOption& stub : stubs) {
    if (term.options.size() >= 1000) {
      break;
    }
    term.options.push_back(std::move(stub));
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
  // Against the current routing space, two pitches farther out (Sec. 4.4
  // "we dynamically generate new access paths").
  offTrackCatalogue(net.id, term, &worker, fresh, stats, 2);
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
  req.node_allowed = [&](const int layer, const br::GridPoint& p) {
    return nodeAllowed(net, static_cast<frLayerNum>(layer), toPoint(p));
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
                       toPoint(patch.end),
                       true});
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
  // Reservations are shape-grid shapes owned by the net, non-removable for
  // the others (Sec. 4.3).
  for (const RShape& shape : term.reservation) {
    changeShape(shape, br::ShapeType::kReservation, net, active, local);
  }
  term.reservation_active = active;
}

void BonnRouteDR::setReservations(NetData& net,
                                  const bool active,
                                  Worker* local,
                                  br::FastGridBuilder::PendingRemovals* batch)
{
  if (active) {
    for (Term& term : net.terms) {
      setReservation(term, net.id, true, local);
    }
    return;
  }
  // Released together: one exact fast-grid recomputation.
  br::FastGridBuilder::PendingRemovals own;
  br::FastGridBuilder::PendingRemovals& pending
      = batch != nullptr ? *batch : own;
  for (Term& term : net.terms) {
    if (term.reservation_active) {
      for (const RShape& shape : term.reservation) {
        changeShape(
            shape, br::ShapeType::kReservation, net.id, false, local, &pending);
      }
      term.reservation_active = false;
    }
  }
  if (net.input_reserved) {
    for (const RShape& shape : net.input_reservation) {
      changeShape(
          shape, br::ShapeType::kReservation, net.id, false, local, &pending);
    }
    net.input_reserved = false;
  }
  if (batch == nullptr) {
    fast_builder_.flush(pending);
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
    changeShape(shape, br::ShapeType::kReservation, net.id, active, local);
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
      net.vias.push_back({type_vias_[net.wire_type][std::min(za, zb)],
                          odb::Point(x, y),
                          net.wire_type});
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
    net.wires.push_back({layers_[za].fr,
                         odb::Point(x0, y0),
                         odb::Point(x1, y1),
                         false,
                         net.wire_type});
    k = e;
  }
}

bool BonnRouteDR::connect(NetData& net,
                          const bool ripup,
                          Worker& worker,
                          const int extra_dilation)
{
  // Sec. 4.4: "The shapes of the components that contributed to S and T,
  // including reserved pin access paths, are removed temporarily from
  // routing space in order not to block access to these vertices." All
  // components of the net are sources or targets of its searches. The
  // searches read the fast grid through an exact overlay with the net's
  // shapes removed, and every checker query of the connection ignores them
  // (see DataStructures.md); the routing space itself is not changed.
  ownShapes(net, worker.own_shapes);
  fast_builder_.buildOverlay(net.id, worker.own_shapes, worker.overlay);
  worker.reader.reset(&grid_, &worker.overlay);
  worker.connects++;
  worker.own_shape_count += static_cast<int64_t>(worker.own_shapes.size());
  const bool ok = connectWithoutOwnShapes(net, ripup, worker, extra_dilation);
  worker.overlay.clear();
  worker.reader.reset(&grid_, nullptr);
  return ok;
}

bool BonnRouteDR::connectWithoutOwnShapes(NetData& net,
                                          const bool ripup,
                                          Worker& worker,
                                          const int extra_dilation)
{
  const int num_terms = static_cast<int>(net.terms.size());
  // Stub and attachment legality depend on the committed shapes, not on
  // the guide dilation. Check them once for both clean windows.
  std::vector<std::vector<char>> valid(num_terms);
  auto option_valid = [&](const AccessOption& opt) {
    const br::RipupLevel level
        = worker.reader.level(opt.vertex, br::kWire, net.wire_type);
    bool ok = level == br::kLegal || (ripup && level <= br::kMaxRemovableLevel);
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
      const br::RipupLevel level
          = worker.reader.level(opt.vertex, br::kWire, net.wire_type);
      bool ok
          = level == br::kLegal || (ripup && level <= br::kMaxRemovableLevel);
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
          ok = !shapeConflicts({v.def->getLayer2Num(),
                                via.getLayer2BBox(),
                                br::ShapeRole::kViaPad},
                               net.id,
                               &worker);
        }
      }
      valid[t][o] = ok;
      feasible += ok;
    }
    // Sec. 4.4: "If no or only a small number of access paths remains
    // feasible, we dynamically generate new access paths."
    if (!ripup && params_.offtrack_access && !worker.speculative
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
      worker.access_failures++;
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
  // Effort budget of a speculative tail attempt (tail_defer_settled):
  // checked before each search, so an attempt settles at most the budget
  // plus one search; the net is then deferred (see routeSpeculative).
  auto over_budget = [&worker]() {
    if (worker.effort_budget > 0
        && worker.settled - worker.effort_start >= worker.effort_budget) {
      worker.effort_exceeded = true;
    }
    return worker.effort_exceeded;
  };
  if (over_budget()) {
    return false;
  }
  const int num_terms = static_cast<int>(net.terms.size());
  br::SearchRequest req;
  req.net = net.id;
  req.region = &worker.region;
  req.ripup = ripup;
  req.ripup_cost
      = static_cast<Cost>(params_.ripup_cost) * layers_[0].pitch * unit();
  req.discounts = net.discount_runs;
  req.wire_type = net.wire_type;
  req.label_limit = label_limit_;
  req.overlay = &worker.overlay;
  // Clears the net's spreading mask when the connection ends.
  struct MaskGuard
  {
    Worker& worker;
    const NetData& net;
    ~MaskGuard()
    {
      for (const auto& [z, tile] : net.extra_space_tiles) {
        worker.spread_mask[z][tile] = 0;
      }
    }
  };
  std::optional<MaskGuard> mask_guard;
  if (params_.spreading && extra_space_nets_ > 0) {
    // Sec. 2.1 / 4.2: keep tracks free where GR gave this net extra space.
    if (worker.spread_mask.size() != static_cast<size_t>(graph_.numLayers())) {
      worker.spread_mask.assign(
          graph_.numLayers(),
          std::vector<uint8_t>(static_cast<size_t>(tiles_x_) * tiles_y_, 0));
    }
    for (const auto& [z, tile] : net.extra_space_tiles) {
      worker.spread_mask[z][tile] = 1;
    }
    mask_guard.emplace(MaskGuard{worker, net});
    req.spread_tiles = &worker.spread_mask;
    req.spread_cost
        = static_cast<Cost>(params_.spread_cost) * layers_[0].pitch * unit();
  } else if (params_.spreading && !spread_tiles_.empty()) {
    req.spread_tiles = &spread_tiles_;
    req.spread_cost
        = static_cast<Cost>(params_.spread_cost) * layers_[0].pitch * unit();
  }
  req.discount_rate = static_cast<Cost>(params_.discount_numerator) * unit()
                      / params_.discount_denominator;
  req.jog_low_tiles = jog_low_tiles_.empty() ? nullptr : &jog_low_tiles_;
  req.use_intervals = params_.use_intervals;
  req.interval_via_property = params_.interval_via_property;
  req.interval_nearest_first = params_.interval_nearest_first;
  req.exact_future_cost = params_.exact_future_cost;
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
        req.target_boxes.push_back({x, y, x, y, opt.z, opt.z, opt.cost});
        continue;
      }
      if (t != box_term) {
        req.target_boxes.push_back({x, y, x, y, opt.z, opt.z, opt.cost});
        box_term = t;
        continue;
      }
      br::TargetBox& box = req.target_boxes.back();
      box.cost = std::min(box.cost, opt.cost);
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
    const double prepare_start = omp_get_wtime();
    worker.potential.prepare(model,
                             graph_,
                             req,
                             gcells,
                             params_.discount_numerator,
                             params_.discount_denominator);
    worker.potential_seconds += omp_get_wtime() - prepare_start;
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
      const br::SearchResult plain = params_.use_intervals
                                         ? worker.interval_search->run(req)
                                         : worker.search->run(req);
      req.future_cost = oracle;
      worker.potential_checks++;
      if (plain.found != res.found || (res.found && plain.cost != res.cost)) {
        worker.potential_mismatches++;
      }
    }
    worker.searches++;
    worker.heap_pops += res.heap_pops;
    worker.settled += res.settled;
    worker.fast_grid_queries += res.fast_grid_queries;
    worker.checker_queries += res.checker_queries;
    worker.aborted_searches += res.aborted;
    return res;
  };
  bool full_region = false;
  int remaining = num_terms - 1;
  while (remaining > 0) {
    if (over_budget()) {
      return false;
    }
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
      if (over_budget()) {
        return false;
      }
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
    // NDR nets: on-track paths with the wire type, access paths of the
    // default rule (tapered at the pins, as drt does); the same-net
    // post-processing models the default rule only.
    if (!params_.postprocess || net.wire_type != 0) {
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
      if (worker.reader.level(v, br::kWire, net.wire_type) != br::kLegal) {
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
    changeShape(s, routedType(s.role), net.id, true, local);
  }
  net.routed = true;
}

int64_t BonnRouteDR::collapseBridges()
{
  // A bridge: vias a, b between the same layers L (below) and U (above) at
  // points A != B on one line, joined by exactly one wire of the net on U
  // from A to B, with nothing else of the net on U touching that wire (no
  // other wire, via, patch or pin point). The vias and the U wire are
  // replaced by a wire on L from A to B in L's non-preferred direction,
  // when that wire is legal against all other nets (the distance rule
  // checker, own shapes excluded) and L allows wrong-way wires.
  auto on_segment
      = [](const odb::Point& p, const odb::Point& a, const odb::Point& b) {
          return p.x() >= std::min(a.x(), b.x())
                 && p.x() <= std::max(a.x(), b.x())
                 && p.y() >= std::min(a.y(), b.y())
                 && p.y() <= std::max(a.y(), b.y())
                 && (a.x() == b.x() ? p.x() == a.x() : p.y() == a.y());
        };
  auto segments_touch
      = [&](const OutWire& w, const odb::Point& a, const odb::Point& b) {
          // Axis-parallel segments: their bounding boxes intersect.
          const odb::Rect s(w.begin, w.end);
          const odb::Rect t(a, b);
          return s.intersects(t);
        };
  int64_t collapsed = 0;
  for (NetData& net : nets_) {
    if (!net.routed || net.vias.size() < 2) {
      continue;
    }
    bool changed = true;
    while (changed) {
      changed = false;
      for (size_t ia = 0; ia < net.vias.size() && !changed; ia++) {
        for (size_t ib = ia + 1; ib < net.vias.size() && !changed; ib++) {
          const OutVia va = net.vias[ia];
          const OutVia vb = net.vias[ib];
          const frLayerNum lower = va.def->getLayer1Num();
          const frLayerNum upper = va.def->getLayer2Num();
          if (vb.def->getLayer1Num() != lower || vb.def->getLayer2Num() != upper
              || va.point == vb.point
              || (va.point.x() != vb.point.x()
                  && va.point.y() != vb.point.y())) {
            continue;
          }
          for (int mode = 0; mode < 2 && !changed; mode++) {
            // mode 0: the bridge wire is on the upper layer and the jog on the
            // lower one; mode 1 (collapse_lower_bridges): a wire on the lower
            // layer between two vias up becomes a wrong-way wire on the upper
            // layer, when the lower piece only joins two pieces of the upper
            // layer (sky130: most up-bridges cross met1 power rails and
            // cannot become met1 jogs).
            if (mode == 1 && !params_.collapse_lower_bridges) {
              continue;
            }
            const frLayerNum span = mode == 0 ? upper : lower;
            const frLayerNum jogl = mode == 0 ? lower : upper;
            BridgeStats& stats
                = mode == 0 ? bridge_stats_ : lower_bridge_stats_;
            const int zl = jogl >= 0 && jogl < static_cast<int>(z_of_fr_.size())
                               ? z_of_fr_[jogl]
                               : -1;
            if (zl < 0 || onTracksOnly(jogl)) {
              continue;
            }
            const odb::Point A = va.point;
            const odb::Point B = vb.point;
            int bridge = -1;
            int other_wires = 0;
            int patches = 0;
            for (size_t k = 0; k < net.wires.size(); k++) {
              const OutWire& w = net.wires[k];
              if (w.layer != span || !segments_touch(w, A, B)) {
                continue;
              }
              const bool same = !w.patch
                                && ((w.begin == A && w.end == B)
                                    || (w.begin == B && w.end == A));
              if (same && bridge < 0) {
                bridge = static_cast<int>(k);
              } else if (w.patch) {
                patches++;
              } else {
                other_wires++;
              }
            }
            if (bridge < 0) {
              continue;
            }
            stats.candidates++;
            if (other_wires > 0 || patches > 0) {
              (patches > 0 && other_wires == 0 ? stats.patch
                                               : stats.other_wire)++;
              continue;
            }
            bool clean = true;
            for (size_t k = 0; k < net.vias.size() && clean; k++) {
              if (k == ia || k == ib) {
                continue;
              }
              const OutVia& v = net.vias[k];
              if ((v.def->getLayer1Num() == span
                   || v.def->getLayer2Num() == span)
                  && on_segment(v.point, A, B)) {
                clean = false;
              }
            }
            for (const auto& [layer, p] : net.pin_points) {
              if (layer == span && on_segment(p, A, B)) {
                clean = false;
              }
            }
            if (!clean) {
              stats.other_via_or_pin++;
              continue;
            }
            const OutWire jog{jogl, A, B, false, net.wire_type};
            const RShape jog_shape = wireShape(jog);
            if (shapeConflicts(jog_shape, net.id, nullptr)) {
              stats.illegal++;
              // Blocked only by routed shapes (rip-up could free the jog) or
              // by fixed shapes (pins, blockages, power).
              const int g = gridLayer(jogl);
              const br::RipupLevel level
                  = g < 0 ? br::kNotRemovable
                          : checker_.shapeLevel(g,
                                                brRect(jog_shape.rect),
                                                shapeClass(jog_shape, g),
                                                3,
                                                net.id);
              (level == br::kNotRemovable ? stats.illegal_fixed
                                          : stats.illegal_routed)++;
              continue;
            }
            // Remove the bridge from the routing space, insert the jog.
            std::vector<RShape> old;
            old.push_back(wireShape(net.wires[bridge]));
            viaShapes(va, old);
            viaShapes(vb, old);
            br::FastGridBuilder::PendingRemovals pending;
            for (const RShape& s : old) {
              changeShape(
                  s, routedType(s.role), net.id, false, nullptr, &pending);
            }
            fast_builder_.flush(pending);
            changeShape(
                jog_shape, routedType(jog_shape.role), net.id, true, nullptr);
            net.wires.erase(net.wires.begin() + bridge);
            net.wires.push_back(jog);
            net.vias.erase(net.vias.begin() + ib);
            net.vias.erase(net.vias.begin() + ia);
            shapesOf(net, net.shapes);
            collapsed++;
            stats.collapsed++;
            changed = true;
          }
        }
      }
    }
  }
  return collapsed;
}

void BonnRouteDR::ripupNet(NetData& net)
{
  // One exact fast-grid recomputation for all shapes of the net.
  br::FastGridBuilder::PendingRemovals pending;
  for (const RShape& s : net.shapes) {
    changeShape(s, routedType(s.role), net.id, false, nullptr, &pending);
  }
  fast_builder_.flush(pending);
  net.shapes.clear();
  net.wires.clear();
  net.vias.clear();
  net.routed = false;
}

void BonnRouteDR::collectVictims(const NetData& net,
                                 std::vector<int>& victims) const
{
  // Sec. 3.4: "the distance rule checking module can report nets whose
  // (partial) removal increases the set of allowed wire or via models";
  // Sec. 4.2 / 4.4: those nets are ripped up.
  std::vector<RShape> shapes;
  shapesOf(net, shapes);
  std::unordered_set<int> found;
  std::vector<br::BlockingNet> nets;
  auto report = [&](const RShape& shape, const int g) {
    checker_.shapeLevel(
        g, brRect(shape.rect), shapeClass(shape, g), 3, net.id, &nets);
    for (const br::BlockingNet& b : nets) {
      if (b.net != net.id && b.net >= 0) {
        found.insert(b.net);
      }
    }
  };
  for (const RShape& shape : shapes) {
    const int g = gridLayer(shape.layer);
    if (g < 0) {
      continue;
    }
    report(shape, g);
    const int above = projectionLayer(g);
    if (above >= 0) {
      RShape projection = shape;
      projection.role = br::ShapeRole::kCutProjection;
      report(projection, above);
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
  // Pins are removed and reinserted while the net is connected (Sec. 4.4).
  for (const RShape& pin : net.pins) {
    add(tileX(pin.rect.xMin()), tileY(pin.rect.yMin()));
    add(tileX(pin.rect.xMax()), tileY(pin.rect.yMax()));
  }
  // Retry dilation plus the routing-space interaction halo (Sec. 5.1).
  const int halo = params_.retry_dilation + halo_tiles_;
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

void BonnRouteDR::routeSpeculative(const std::vector<int>& order,
                                   std::vector<int>& failed,
                                   const char* stage)
{
  const int threads
      = params_.tail_threads > 0
            ? std::min(params_.tail_threads, static_cast<int>(workers_.size()))
            : static_cast<int>(workers_.size());
  // Nets over the effort budget are routed after all others, in net order,
  // in a second pass without a budget.
  std::vector<int> deferred;
  std::vector<int> none;
  const int64_t budget = params_.tail_defer_settled;
  if (params_.tail_pipeline) {
    routeSpeculativePipeline(order, failed, deferred, stage, threads, budget);
  } else {
    routeSpeculativeBatches(order, failed, deferred, stage, threads, budget);
  }
  if (deferred.empty()) {
    return;
  }
  const std::string late = std::string(stage) + " deferred";
  if (params_.tail_pipeline) {
    routeSpeculativePipeline(deferred, failed, none, late, threads, 0);
  } else {
    routeSpeculativeBatches(deferred, failed, none, late, threads, 0);
  }
}

void BonnRouteDR::routeSpeculativeBatches(const std::vector<int>& order,
                                          std::vector<int>& failed,
                                          std::vector<int>& deferred,
                                          const std::string& stage,
                                          const int threads,
                                          const int64_t budget)
{
  const size_t batch_size = static_cast<size_t>(params_.tail_batch);
  std::deque<int> queue(order.begin(), order.end());
  std::vector<int> batch;
  std::vector<char> ok;
  std::vector<double> busy;
  std::vector<int64_t> settled;
  std::vector<char> over;
  int64_t deferred_count = 0;
  std::vector<int> slots;
  std::vector<double> estimate;
  std::vector<RShape> shapes;
  int64_t batches = 0;
  int64_t conflicts = 0;
  int64_t nets = 0;
  int64_t checked = 0;
  double commit_seconds = 0;
  double check_seconds = 0;
  double flush_seconds = 0;
  double wall_total = 0;
  double busy_total = 0;
  // Balance of the batches: the longest search of each batch, the thread
  // time idle at the batch barriers, and every search time of the stage.
  double longest_total = 0;
  double idle_total = 0;
  std::vector<double> search_times;
  while (!queue.empty()) {
    batch.clear();
    while (!queue.empty() && batch.size() < batch_size) {
      batch.push_back(queue.front());
      queue.pop_front();
    }
    ok.assign(batch.size(), 0);
    busy.assign(batch.size(), 0);
    settled.assign(batch.size(), 0);
    over.assign(batch.size(), 0);
    // Searches in dispatch order: by default the largest estimated effort
    // first, so a long search does not start last. Every search reads the
    // same routing space, so the order only changes the batch wall time.
    slots.resize(batch.size());
    std::iota(slots.begin(), slots.end(), 0);
    if (params_.tail_order) {
      estimate.resize(batch.size());
      for (size_t k = 0; k < batch.size(); k++) {
        estimate[k] = tailSearchEstimate(nets_[batch[k]]);
      }
      std::ranges::stable_sort(
          slots, [&](int a, int b) { return estimate[a] > estimate[b]; });
    }
    const double wall_start = omp_get_wtime();
    // Nothing writes the routing space while the batch is searched.
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (int i = 0; i < static_cast<int>(batch.size()); i++) {
      const int k = slots[i];
      Worker& w = workers_[omp_get_thread_num()];
      const double start = omp_get_wtime();
      const int64_t settled_before = w.settled;
      w.speculative = true;
      w.effort_budget = budget;
      w.effort_start = w.settled;
      w.effort_exceeded = false;
      ok[k] = connect(nets_[batch[k]], false, w);
      over[k] = !ok[k] && w.effort_exceeded;
      w.effort_budget = 0;
      w.effort_exceeded = false;
      w.speculative = false;
      settled[k] = w.settled - settled_before;
      busy[k] = omp_get_wtime() - start;
    }
    const double wall = omp_get_wtime() - wall_start;
    wall_total += wall;
    double batch_busy = 0;
    double longest = 0;
    for (const double b : busy) {
      batch_busy += b;
      longest = std::max(longest, b);
      search_times.push_back(b);
    }
    busy_total += batch_busy;
    longest_total += longest;
    idle_total += wall * threads - batch_busy;
    debugPrint(logger_,
               utl::DRT,
               "bonnroute_tail",
               1,
               "tail batch {} nets {} wall_ms {:.3f} longest_ms {:.3f} "
               "busy_ms {:.3f}",
               batches,
               batch.size(),
               1e3 * wall,
               1e3 * longest,
               1e3 * batch_busy);
    // Commit in net order. Each route was legal on the routing space at the
    // start of the batch, so only commits earlier in the batch can conflict
    // with it. The first committed net of a batch is never checked, which
    // guarantees progress.
    // Only shapes within the interaction distance of a shape committed
    // earlier in the batch can conflict: nets whose bounding boxes stay
    // farther apart skip the exact check.
    const double commit_start = omp_get_wtime();
    // The reservation releases of the whole batch are flushed together
    // after its commits. No search runs before the flush, inserts take the
    // maximum with stored levels, and the flush recomputes every removal
    // candidate exactly against the current shape grid, so the fast grid
    // is exact again before the next batch is searched.
    br::FastGridBuilder::PendingRemovals released;
    std::vector<int> retry;
    std::vector<odb::Rect> committed_boxes;
    for (size_t k = 0; k < batch.size(); k++) {
      NetData& net = nets_[batch[k]];
      // Per-net trace: 0 failed, 1 committed, 2 conflict, 3 deferred (over
      // the effort budget); exact check run.
      int outcome = 0;
      bool exact = false;
      if (!ok[k]) {
        net.wires.clear();
        net.vias.clear();
        if (over[k]) {
          deferred.push_back(net.id);
          deferred_count++;
          outcome = 3;
        } else {
          failed.push_back(net.id);
        }
      } else {
        bool clash = false;
        if (!committed_boxes.empty()) {
          shapesOf(net, shapes);
          odb::Rect box;
          box.mergeInit();
          for (const RShape& shape : shapes) {
            box.merge(shape.rect);
          }
          const bool near = !box.isInverted()
                            && std::ranges::any_of(
                                committed_boxes, [&](const odb::Rect& other) {
                                  return other.intersects(box);
                                });
          if (near) {
            const double check_start = omp_get_wtime();
            clash = std::ranges::any_of(shapes, [&](const RShape& shape) {
              return shapeConflicts(shape, net.id);
            });
            check_seconds += omp_get_wtime() - check_start;
            checked++;
            exact = true;
          }
        }
        if (clash) {
          net.wires.clear();
          net.vias.clear();
          net.pin_points.clear();
          retry.push_back(net.id);
          conflicts++;
          outcome = 2;
        } else {
          commitNet(net, nullptr);
          setReservations(net, false, nullptr, &released);
          odb::Rect box;
          box.mergeInit();
          for (const RShape& shape : net.shapes) {
            box.merge(shape.rect);
          }
          if (!box.isInverted()) {
            box.bloat(static_cast<int>(interaction_reach_), box);
            committed_boxes.push_back(box);
          }
          nets++;
          outcome = 1;
        }
      }
      if (logger_->debugCheck(utl::DRT, "bonnroute_tail", 1)) {
        int xlo, ylo, xhi, yhi;
        netTileBox(net, xlo, ylo, xhi, yhi);
        logger_->debug(utl::DRT,
                       "bonnroute_tail",
                       "tail net batch {} slot {} net {} search_ms {:.3f} "
                       "settled {} outcome {} exact {} box {} {} {} {} "
                       "guide_tiles {} terms {} hpwl {}",
                       batches,
                       k,
                       net.id,
                       1e3 * busy[k],
                       settled[k],
                       outcome,
                       exact ? 1 : 0,
                       std::max(xlo, 0),
                       std::max(ylo, 0),
                       std::min(xhi, tiles_x_ - 1),
                       std::min(yhi, tiles_y_ - 1),
                       net.guide_tiles.size(),
                       net.terms.size(),
                       net.hpwl);
      }
    }
    const double flush_start = omp_get_wtime();
    fast_builder_.flush(released);
    flush_seconds += omp_get_wtime() - flush_start;
    commit_seconds += omp_get_wtime() - commit_start;
    for (auto it = retry.rbegin(); it != retry.rend(); ++it) {
      queue.push_front(*it);
    }
    batches++;
  }
  partition_wall_ += wall_total;
  partition_busy_ += busy_total;
  partition_capacity_ += wall_total * threads;
  if (params_.verbose) {
    logger_->info(
        utl::DRT,
        9105,
        "BonnRoute {} speculative batches: {} batches of up to {} "
        "nets, {} nets committed, {} searched again after a "
        "conflict; parallel search {:.2f}s, thread utilization {:.0f}%; "
        "serial commit {:.2f}s ({} nets checked exactly, {:.2f}s; "
        "reservation release flush {:.2f}s).",
        stage,
        batches,
        batch_size,
        nets,
        conflicts,
        wall_total,
        wall_total > 0 ? 100.0 * busy_total / (wall_total * threads) : 0.0,
        commit_seconds,
        checked,
        check_seconds,
        flush_seconds);
    std::ranges::sort(search_times);
    auto quantile_ms = [&](const double q) {
      if (search_times.empty()) {
        return 0.0;
      }
      const size_t i = std::min(
          search_times.size() - 1,
          static_cast<size_t>(q * static_cast<double>(search_times.size())));
      return 1e3 * search_times[i];
    };
    logger_->info(
        utl::DRT,
        9125,
        "BonnRoute {} speculative batch balance on {} threads: longest "
        "search per batch {:.2f}s in total ({:.0f}% of the parallel search "
        "wall), idle thread time at the batch barriers {:.2f}s; {} searches, "
        "search time p50 {:.1f} ms, p90 {:.1f} ms, p99 {:.1f} ms, max {:.1f} "
        "ms; {} nets deferred over an effort of {} settled vertices; routing "
        "digest {:016x}.",
        stage,
        threads,
        longest_total,
        wall_total > 0 ? 100.0 * longest_total / wall_total : 0.0,
        idle_total,
        search_times.size(),
        quantile_ms(0.5),
        quantile_ms(0.9),
        quantile_ms(0.99),
        search_times.empty() ? 0.0 : 1e3 * search_times.back(),
        deferred_count,
        budget,
        tailRoutingDigest());
  }
}

double BonnRouteDR::tailSearchEstimate(const NetData& net) const
{
  // Search effort grows with the routing area and the number of
  // connections: guide tiles times terminals.
  return static_cast<double>(net.guide_tiles.size())
         * static_cast<double>(std::max<size_t>(1, net.terms.size()));
}

uint64_t BonnRouteDR::tailRoutingDigest() const
{
  uint64_t digest = 1469598103934665603ULL;
  auto add = [&digest](const int64_t v) {
    digest ^= static_cast<uint64_t>(v);
    digest *= 1099511628211ULL;
  };
  for (const NetData& net : nets_) {
    add(net.id);
    add(net.routed);
    for (const OutWire& w : net.wires) {
      add(w.layer);
      add(w.begin.x());
      add(w.begin.y());
      add(w.end.x());
      add(w.end.y());
      add(w.patch);
    }
    for (const OutVia& v : net.vias) {
      add(v.def->getCutLayerNum());
      add(static_cast<int64_t>(std::hash<std::string>()(v.def->getName())));
      add(v.point.x());
      add(v.point.y());
    }
    for (const auto& [layer, p] : net.pin_points) {
      add(layer);
      add(p.x());
      add(p.y());
    }
  }
  return digest;
}

void BonnRouteDR::routeSpeculativePipeline(const std::vector<int>& order,
                                           std::vector<int>& failed,
                                           std::vector<int>& deferred,
                                           const std::string& stage,
                                           const int threads,
                                           const int64_t budget)
{
  // Tile boxes hold every read and write of a net's search and commit
  // (netTileBox adds the interaction halo, Sec. 5.1), so nets with
  // disjoint boxes are searched and committed concurrently, as threads of
  // disjoint partition regions are.
  std::vector<br::TailPipeline::Box> boxes(nets_.size());
  for (const int id : order) {
    int xlo, ylo, xhi, yhi;
    netTileBox(nets_[id], xlo, ylo, xhi, yhi);
    boxes[id] = {std::max(xlo, 0),
                 std::max(ylo, 0),
                 std::min(xhi, tiles_x_ - 1),
                 std::min(yhi, tiles_y_ - 1)};
  }
  std::vector<double> search_seconds(nets_.size(), 0);
  std::vector<int64_t> settled(nets_.size(), 0);
  std::vector<char> over(nets_.size(), 0);
  // Shapes of committed nets, bloated by the interaction distance.
  std::vector<odb::Rect> committed_box(nets_.size());
  std::vector<std::vector<RShape>> shapes(threads);
  std::vector<int64_t> checked(threads, 0);
  std::vector<double> check_seconds(threads, 0);
  std::vector<std::vector<double>> times(threads);
  const bool trace = logger_->debugCheck(utl::DRT, "bonnroute_tail", 1);
  br::TailPipeline::Hooks hooks;
  hooks.parallel = [](const int count, const std::function<void(int)>& body) {
#pragma omp parallel num_threads(count)
    body(omp_get_thread_num());
  };
  hooks.search = [&](const int id, const int t) {
    Worker& w = workers_[t];
    const double start = omp_get_wtime();
    const int64_t before = w.settled;
    w.speculative = true;
    w.effort_budget = budget;
    w.effort_start = w.settled;
    w.effort_exceeded = false;
    const bool found = connect(nets_[id], false, w);
    over[id] = !found && w.effort_exceeded;
    w.effort_budget = 0;
    w.effort_exceeded = false;
    w.speculative = false;
    settled[id] = w.settled - before;
    search_seconds[id] = omp_get_wtime() - start;
    times[t].push_back(search_seconds[id]);
    return found;
  };
  hooks.decide = [&](const int id,
                     const bool found,
                     const std::vector<int>& earlier,
                     const int t) {
    using Decision = br::TailPipeline::Decision;
    Worker& w = workers_[t];
    NetData& net = nets_[id];
    Decision decision = Decision::kCommitted;
    bool exact = false;
    if (!found) {
      net.wires.clear();
      net.vias.clear();
      decision = Decision::kFailed;
    } else {
      // Only commits earlier in the batch can conflict (the route was
      // legal on the commits of the earlier batches), and only within the
      // interaction distance of their shapes.
      bool clash = false;
      if (!earlier.empty()) {
        shapesOf(net, shapes[t]);
        odb::Rect box;
        box.mergeInit();
        for (const RShape& shape : shapes[t]) {
          box.merge(shape.rect);
        }
        const bool near
            = !box.isInverted() && std::ranges::any_of(earlier, [&](int e) {
                return committed_box[e].intersects(box);
              });
        if (near) {
          const double check_start = omp_get_wtime();
          clash = std::ranges::any_of(shapes[t], [&](const RShape& shape) {
            return shapeConflicts(shape, net.id, &w);
          });
          check_seconds[t] += omp_get_wtime() - check_start;
          checked[t]++;
          exact = true;
        }
      }
      if (clash) {
        net.wires.clear();
        net.vias.clear();
        net.pin_points.clear();
        decision = Decision::kConflict;
      } else {
        commitNet(net, &w);
        // Released and recomputed exactly inside the net's box.
        setReservations(net, false, &w);
        odb::Rect box;
        box.mergeInit();
        for (const RShape& shape : net.shapes) {
          box.merge(shape.rect);
        }
        if (!box.isInverted()) {
          box.bloat(static_cast<int>(interaction_reach_), box);
          committed_box[id] = box;
        }
      }
    }
    if (trace) {
      const auto& b = boxes[id];
      logger_->debug(utl::DRT,
                     "bonnroute_tail",
                     "tail pipeline net {} search_ms {:.3f} settled {} "
                     "outcome {} exact {} box {} {} {} {} guide_tiles {} "
                     "terms {} hpwl {}",
                     id,
                     1e3 * search_seconds[id],
                     settled[id],
                     decision == Decision::kFailed
                         ? (over[id] ? 3 : 0)
                         : (decision == Decision::kCommitted ? 1 : 2),
                     exact ? 1 : 0,
                     b[0],
                     b[1],
                     b[2],
                     b[3],
                     net.guide_tiles.size(),
                     net.terms.size(),
                     net.hpwl);
    }
    return decision;
  };
  if (params_.tail_order) {
    hooks.cost = [&](const int id) { return tailSearchEstimate(nets_[id]); };
  }
  br::TailPipeline pipeline(params_.tail_batch, params_.tail_retry_delay);
  std::vector<int> unrouted;
  const br::TailPipeline::Stats stats
      = pipeline.run(order, boxes, threads, hooks, unrouted);
  int64_t deferred_count = 0;
  for (const int id : unrouted) {
    if (over[id]) {
      deferred.push_back(id);
      deferred_count++;
    } else {
      failed.push_back(id);
    }
  }
  partition_wall_ += stats.wall;
  partition_busy_ += stats.search_busy + stats.decide_busy;
  partition_capacity_ += stats.wall * threads;
  if (params_.verbose) {
    int64_t checked_total = 0;
    double check_total = 0;
    std::vector<double> all;
    for (int t = 0; t < threads; t++) {
      checked_total += checked[t];
      check_total += check_seconds[t];
      all.insert(all.end(), times[t].begin(), times[t].end());
    }
    std::ranges::sort(all);
    auto quantile_ms = [&](const double q) {
      if (all.empty()) {
        return 0.0;
      }
      const size_t i
          = std::min(all.size() - 1,
                     static_cast<size_t>(q * static_cast<double>(all.size())));
      return 1e3 * all[i];
    };
    const double capacity = stats.wall * threads;
    logger_->info(
        utl::DRT,
        9127,
        "BonnRoute {} speculative pipeline: {} batches of {} nets (retries "
        "{} batches later), {} nets committed, {} searched again after a "
        "conflict, {} without a route; {:.2f}s on {} threads, thread "
        "utilization {:.0f}% (searches {:.2f}s, decisions {:.2f}s incl. {} "
        "exact checks {:.2f}s, idle {:.2f}s); {} searches started before an "
        "earlier batch was decided; search time p50 {:.1f} ms, p90 {:.1f} "
        "ms, p99 {:.1f} ms, max {:.1f} ms; {} nets deferred over an effort "
        "of {} settled vertices; routing digest {:016x}.",
        stage,
        stats.batches,
        params_.tail_batch,
        params_.tail_retry_delay,
        stats.committed,
        stats.conflicts,
        stats.failed - deferred_count,
        stats.wall,
        threads,
        capacity > 0
            ? 100.0 * (stats.search_busy + stats.decide_busy) / capacity
            : 0.0,
        stats.search_busy,
        stats.decide_busy,
        checked_total,
        check_total,
        stats.idle,
        stats.overlapped_searches,
        quantile_ms(0.5),
        quantile_ms(0.9),
        quantile_ms(0.99),
        all.empty() ? 0.0 : 1e3 * all.back(),
        deferred_count,
        budget,
        tailRoutingDigest());
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
    if (params_.tail_batch > 0 && count < threads) {
      routeSpeculative(pending, failed, stage);
      return;
    }
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
    std::vector<double> busy(regions.size(), 0);
    const double wall_start = omp_get_wtime();
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (int r = 0; r < static_cast<int>(regions.size()); r++) {
      Worker& w = workers_[omp_get_thread_num()];
      const double region_start = omp_get_wtime();
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
      busy[r] = omp_get_wtime() - region_start;
    }
    const double wall = omp_get_wtime() - wall_start;
    double busy_sum = 0;
    double busy_max = 0;
    for (const double b : busy) {
      busy_sum += b;
      busy_max = std::max(busy_max, b);
    }
    partition_wall_ += wall;
    partition_busy_ += busy_sum;
    partition_capacity_
        += wall * std::min<int>(threads, static_cast<int>(regions.size()));
    for (const auto& f : region_failed) {
      failed.insert(failed.end(), f.begin(), f.end());
    }
    if (params_.verbose) {
      size_t routed = 0;
      for (const auto& m : members) {
        routed += m.size();
      }
      logger_->info(utl::DRT,
                    9107,
                    "BonnRoute {} {}: {} regions, {} nets inside, {} "
                    "deferred; {:.2f}s, longest region {:.2f}s, thread "
                    "utilization {:.0f}%.",
                    stage,
                    partition,
                    regions.size(),
                    routed,
                    next.size(),
                    wall,
                    busy_max,
                    wall > 0 ? 100.0 * busy_sum / (wall * threads) : 0.0);
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
    std::vector<OutWire> patches;
    for (const OutWire& w : data.wires) {
      if (w.patch) {
        patches.push_back(w);
        continue;
      }
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
        wires.push_back({w.layer, cuts[k], cuts[k + 1], w.patch, w.wire_type});
      }
    }
    // Contact guard: FlexDR's clean-up needs every pin contact at one of the
    // pin's access points; otherwise the net is left to FlexDR.
    guard_via_contacts_ += pinContactViolations(data, wires, data.vias, true);
    if (const int64_t bad = pinContactViolations(data, wires, data.vias)) {
      guard_nets_++;
      guard_contacts_ += bad;
      data.routed = false;
      data.wires.clear();
      data.vias.clear();
      continue;
    }
    // Wrong-way or off-track wires on layers that only allow on-track
    // preferred wires (asap7 M1-M3): FlexDR cannot reach a worker-boundary
    // point of such a wire (DRT-0255); leave the net to FlexDR.
    if (!std::ranges::all_of(wires,
                             [&](const OutWire& w) { return wireOnTracks(w); })
        || !std::ranges::all_of(
            patches, [&](const OutWire& w) { return wireOnTracks(w); })) {
      wrong_way_nets_++;
      data.routed = false;
      data.wires.clear();
      data.vias.clear();
      continue;
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
      // NDR wires (Sec. 3.2): the wire type's width, extended by half of
      // it, as FlexDR writes them.
      const int wz = z_of_fr_[w.layer];
      if (w.wire_type > 0 && wz >= 0
          && type_widths_[w.wire_type][wz]
                 > static_cast<int>(style.getWidth())) {
        const frCoord width = type_widths_[w.wire_type][wz];
        style.setWidth(width);
        style.setBeginExt(width / 2);
        style.setEndExt(width / 2);
      }
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
      // Access paths of an NDR net keep the default rule: tapered, as
      // FlexDR marks them, so that its checks apply no NDR spacing.
      seg->setTapered(data.wire_type != 0 && w.wire_type == 0);
      net->addShape(std::move(seg));
    }
    // Min-area extensions as patch wires anchored at their via (FlexDR's
    // own form; the connectivity check keeps patches at via origins).
    for (const OutWire& w : patches) {
      const RShape s = wireShape(w);
      auto pwire = std::make_unique<frPatchWire>();
      pwire->setLayerNum(w.layer);
      pwire->setOrigin(w.begin);
      odb::Rect offset = s.rect;
      offset.moveDelta(-w.begin.x(), -w.begin.y());
      pwire->setOffsetBox(offset);
      net->addPatchWire(std::move(pwire));
    }
    // A via pad at an access point of a pin of the net is that pin's
    // connection, as FlexDR marks its own vias (FlexDR_maze.cpp); unmarked,
    // FlexDR's connectivity check does not see the pin connected.
    for (const OutVia& v : data.vias) {
      auto via = std::make_unique<frVia>(v.def, v.point);
      if (atAccessPoint(data, v.def->getLayer1Num(), v.point)) {
        via->setBottomConnected(true);
        connected_vias_++;
      }
      if (atAccessPoint(data, v.def->getLayer2Num(), v.point)) {
        via->setTopConnected(true);
        connected_vias_++;
      }
      via->setTapered(data.wire_type != 0 && v.wire_type == 0);
      net->addVia(std::move(via));
    }
    net->setHasInitialRouting(true);
    net->setModified(true);
  }
}

namespace {

// FNV-1a style 64-bit digest of integer sequences (log fingerprints of the
// setup and routing results; equal runs print equal digests).
struct Digest
{
  uint64_t value = 1469598103934665603ULL;
  void add(const int64_t v)
  {
    value ^= static_cast<uint64_t>(v);
    value *= 1099511628211ULL;
  }
  void add(const odb::Point& p)
  {
    add(p.x());
    add(p.y());
  }
  void add(const odb::Rect& r)
  {
    add(r.xMin());
    add(r.yMin());
    add(r.xMax());
    add(r.yMax());
  }
};

}  // namespace

uint64_t BonnRouteDR::accessDigest() const
{
  Digest d;
  for (const NetData& net : nets_) {
    d.add(net.id);
    d.add(net.skip);
    for (const Term& term : net.terms) {
      d.add(static_cast<int64_t>(term.options.size()));
      d.add(term.primary);
      d.add(term.reservation_active);
      for (const AccessOption& o : term.options) {
        d.add(static_cast<int64_t>(o.vertex));
        d.add(o.cost);
        d.add(o.z);
        d.add(o.planar);
        d.add(o.ap);
        d.add(o.ap_layer);
        d.add(o.offtrack);
        d.add(o.pin_via);
        for (const OutWire& w : o.wires) {
          d.add(w.layer);
          d.add(w.begin);
          d.add(w.end);
          d.add(w.patch);
        }
        for (const OutVia& v : o.vias) {
          d.add(v.def->getCutLayerNum());
          d.add(
              static_cast<int64_t>(std::hash<std::string>()(v.def->getName())));
          d.add(v.point);
        }
        for (const br::RouteNode& n : o.route) {
          d.add(n.layer);
          d.add(n.p.x);
          d.add(n.p.y);
        }
      }
      for (const RShape& r : term.reservation) {
        d.add(r.layer);
        d.add(r.rect);
        d.add(static_cast<int64_t>(r.role));
      }
    }
  }
  return d.value;
}

uint64_t BonnRouteDR::routingDigest() const
{
  Digest d;
  for (const NetData& net : nets_) {
    d.add(net.id);
    d.add(net.routed);
    for (const OutWire& w : net.wires) {
      d.add(w.layer);
      d.add(w.begin);
      d.add(w.end);
      d.add(w.patch);
    }
    for (const OutVia& v : net.vias) {
      d.add(v.def->getCutLayerNum());
      d.add(static_cast<int64_t>(std::hash<std::string>()(v.def->getName())));
      d.add(v.point);
    }
    for (const auto& [layer, p] : net.pin_points) {
      d.add(layer);
      d.add(p);
    }
  }
  return d.value;
}

void BonnRouteDR::main()
{
  utl::Timer timer;
  if (params_.verbose) {
    logger_->report("BonnRoute DR profile {}.", params_.profileName());
  }
  // TODO(2021): discount rates in the interval search.
  if (modernSearch() && params_.use_intervals) {
    // The interval search assumes unit wire costs (no discounted runs).
    logger_->warn(utl::DRT,
                  9106,
                  "BonnRoute: reservations use the vertex search; "
                  "use_intervals is ignored.");
    params_.use_intervals = false;
  }
  SetupPhases& phases = setup_phases_;
  utl::Timer phase_timer;
  // Wall seconds since the previous lap.
  auto lap = [&phase_timer]() {
    const double t = phase_timer.elapsed();
    phase_timer.reset();
    return t;
  };
  buildTrackGraph();
  phases.track_graph = lap();
  if (layers_.empty()) {
    return;
  }
  buildRoutingSpace();
  phases.routing_space = lap();
  if (params_.spreading) {
    buildSpreadTiles();
  }
  if (params_.jog_factor_low > 0
      && params_.jog_factor_low < params_.jog_factor) {
    buildJogTiles();
    if (params_.verbose) {
      logger_->report(
          "BonnRoute DR jogs: factor {} in {} layer GCells at most {}% used, "
          "{} elsewhere.",
          params_.jog_factor_low,
          jog_low_tile_count_,
          params_.jog_low_utilization,
          params_.jog_factor);
    }
  }
  phases.spreading = lap();
  workers_.resize(std::max(1, router_cfg_->MAX_THREADS));
  for (Worker& w : workers_) {
    w.region.init(graph_.numLayers(), tiles_x_, tiles_y_);
    if (params_.use_intervals) {
      w.interval_search = std::make_unique<br::IntervalSearch>(graph_, grid_);
    } else {
      w.search = std::make_unique<br::PathSearch>(graph_, grid_);
    }
  }
  phases.workers = lap();
  buildNets();
  phases.nets = lap();
  buildFixed();
  lap();
  buildSameNetRules();
  phases.same_net_rules = lap();
  // Reads the extra space of the nets built above (grt/BonnRouteExtraSpace.h).
  if (params_.spreading && params_.spread_from_extra_space) {
    loadExtraSpace();
  }
  phases.spreading += lap();
  std::vector<char> pre_skip(nets_.size());
  for (const NetData& net : nets_) {
    pre_skip[net.id] = net.skip;
  }
  // Guide tiles and access stubs of each net read only the design, the
  // track graph and the routing space; buildAccess may mark its own net
  // skipped.
  std::vector<int> unskipped;
  for (const NetData& net : nets_) {
    if (!net.skip) {
      unskipped.push_back(net.id);
    }
  }
  const int setup_threads = static_cast<int>(workers_.size());
#pragma omp parallel for num_threads(setup_threads) schedule(dynamic, 64)
  for (int k = 0; k < static_cast<int>(unskipped.size()); k++) {
    NetData& net = nets_[unskipped[k]];
    buildGuideTiles(net);
    buildAccess(net);
  }
  phases.guide_tiles_and_stubs = lap();
  if (params_.offtrack_access) {
    buildOffTrackAccess(pre_skip);
    phases.offtrack_access = lap();
  } else {
    for (NetData& net : nets_) {
      if (!net.skip) {
        reserveAccess(net);
      }
    }
    phases.stub_reservations = lap();
  }
  if (params_.reservations) {
    reserveInputWires();
  }
  phases.input_wires = lap();
  const double setup_time = timer.elapsed();
  if (params_.verbose) {
    logger_->report(
        "BonnRoute DR setup phases: track graph {:.2f}s, routing space "
        "{:.2f}s, spreading tiles {:.2f}s, workers {:.2f}s, nets {:.2f}s, "
        "fixed shapes {:.2f}s (region query {:.2f}s, classes {:.2f}s), fast "
        "grid rebuild {:.2f}s, compaction "
        "{:.2f}s, same-net rules {:.2f}s, guide tiles and access stubs "
        "{:.2f}s, off-track access {:.2f}s, stub reservations {:.2f}s, input "
        "wire reservations {:.2f}s; setup {:.2f}s.",
        phases.track_graph,
        phases.routing_space,
        phases.spreading,
        phases.workers,
        phases.nets,
        phases.fixed_insert,
        phases.fixed_query,
        phases.fixed_classes,
        phases.fast_grid_rebuild,
        phases.compaction,
        phases.same_net_rules,
        phases.guide_tiles_and_stubs,
        phases.offtrack_access,
        phases.stub_reservations,
        phases.input_wires,
        setup_time);
    logger_->report("BonnRoute DR access digest {:016x}.", accessDigest());
  }
  std::vector<int> order;
  for (const NetData& net : nets_) {
    if (!net.skip) {
      order.push_back(net.id);
    }
  }
  std::ranges::stable_sort(
      order, [&](int a, int b) { return nets_[a].hpwl < nets_[b].hpwl; });
  std::vector<int> failed_first;
  // Sec. 5.1: critical nets (clock nets and NDR nets) first, on the whole
  // chip.
  std::vector<int> regular;
  std::vector<int> critical;
  for (const int id : order) {
    // Sec. 5.1: also nets "that use wide wires that would be hard to route
    // when the majority of nets has already been routed" (NDR nets).
    const bool is_critical
        = params_.critical_nets
          && (nets_[id].net->isClock() || nets_[id].wire_type != 0);
    (is_critical ? critical : regular).push_back(id);
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
  if (params_.collapse_bridges) {
    const double start = timer.elapsed();
    const int64_t collapsed = collapseBridges();
    if (params_.verbose) {
      logger_->report(
          "BonnRoute DR bridges: {} collapsed into wrong-way wires ({} vias "
          "fewer); kept: {} with other wires on the bridge, {} with "
          "patches, {} with other vias or pin points, {} illegal against "
          "other nets ({} only routed shapes, {} fixed shapes) (of {} "
          "candidate checks).",
          bridge_stats_.collapsed,
          2 * bridge_stats_.collapsed,
          bridge_stats_.other_wire,
          bridge_stats_.patch,
          bridge_stats_.other_via_or_pin,
          bridge_stats_.illegal,
          bridge_stats_.illegal_routed,
          bridge_stats_.illegal_fixed,
          bridge_stats_.candidates);
      if (params_.collapse_lower_bridges) {
        logger_->report(
            "BonnRoute DR lower bridges: {} collapsed into wrong-way wires ({} "
            "vias "
            "fewer); kept: {} with other wires on the bridge, {} with "
            "patches, {} with other vias or pin points, {} illegal against "
            "other nets ({} only routed shapes, {} fixed shapes) (of {} "
            "candidate checks).",
            lower_bridge_stats_.collapsed,
            2 * lower_bridge_stats_.collapsed,
            lower_bridge_stats_.other_wire,
            lower_bridge_stats_.patch,
            lower_bridge_stats_.other_via_or_pin,
            lower_bridge_stats_.illegal,
            lower_bridge_stats_.illegal_routed,
            lower_bridge_stats_.illegal_fixed,
            lower_bridge_stats_.candidates);
      }
      logger_->report("BonnRoute DR bridge collapse: {} bridges, {:.2f}s.",
                      collapsed,
                      timer.elapsed() - start);
    }
  }
  if (params_.verbose) {
    logger_->report("BonnRoute DR routing digest {:016x}.", routingDigest());
  }
  writeNets();
  if (params_.verbose) {
    // Fingerprint of the routes handed to FlexDR (FNV-1a over nets, wires
    // and vias in net order): equal routes give equal digests across builds.
    uint64_t digest = 14695981039346656037ULL;
    auto mix = [&digest](const int64_t value) {
      digest = (digest ^ static_cast<uint64_t>(value)) * 1099511628211ULL;
    };
    int64_t digest_wires = 0;
    int64_t digest_vias = 0;
    for (const NetData& net : nets_) {
      mix(net.id);
      mix(net.routed);
      for (const OutWire& w : net.wires) {
        mix(w.layer);
        mix(w.begin.x());
        mix(w.begin.y());
        mix(w.end.x());
        mix(w.end.y());
        mix(w.patch);
        digest_wires++;
      }
      for (const OutVia& v : net.vias) {
        for (const char c : v.def->getName()) {
          mix(c);
        }
        mix(v.point.x());
        mix(v.point.y());
        digest_vias++;
      }
    }
    logger_->report("BonnRoute DR output digest {:016x}: {} wires, {} vias.",
                    digest,
                    digest_wires,
                    digest_vias);
  }
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
  int64_t fast_grid_queries = 0;
  int64_t checker_queries = 0;
  int64_t aborted_searches = 0;
  int64_t connects = 0;
  int64_t own_shape_count = 0;
  size_t peak_scratch = 0;
  for (const Worker& w : workers_) {
    fast_grid_queries += w.fast_grid_queries;
    checker_queries += w.checker_queries;
    aborted_searches += w.aborted_searches;
    connects += w.connects;
    own_shape_count += w.own_shape_count;
    peak_scratch = std::max(
        peak_scratch,
        w.search ? w.search->peakScratchBytes()
                 : (w.interval_search ? w.interval_search->peakScratchBytes()
                                      : size_t{0}));
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
    const br::ShapeGrid::Statistics sg = shape_grid_.statistics();
    const br::FastGrid::Statistics fg = grid_.statistics();
    logger_->report(
        "BonnRoute DR routing space: {} shape-grid intervals ({} cells, {} "
        "shared cells), {} cell configurations, {} shape classes, {} wire "
        "types, {:.1f} MB; fast grid {} intervals for {} vertices, {} via "
        "intervals, {} marked edges, {:.1f} MB; partition halo {} tiles; {} "
        "legality queries, {:.2f}% answered by the fast grid; {} "
        "inconsistent updates.",
        sg.intervals,
        sg.interval_cells,
        sg.shared_cells,
        sg.configurations,
        shape_classes_.size(),
        wire_types_.size(),
        sg.bytes / 1e6,
        fg.wiring_intervals,
        fg.vertices,
        fg.via_intervals,
        fg.marked_edges,
        fg.bytes / 1e6,
        halo_tiles_,
        fast_grid_queries + checker_queries,
        fast_grid_queries + checker_queries > 0
            ? 100.0 * fast_grid_queries / (fast_grid_queries + checker_queries)
            : 100.0,
        inconsistent_updates_.load());
    logger_->report(
        "BonnRoute DR search scratch: peak {:.1f} MB per worker, budget {} "
        "MB, {} searches stopped at the budget.",
        peak_scratch / 1e6,
        params_.search_label_mb,
        aborted_searches);
    logger_->report(
        "BonnRoute DR own shapes: {} connects, {:.1f} own shapes (pins, "
        "reservations, cut projections) per connect, read through the "
        "fast-grid overlay.",
        connects,
        connects > 0 ? static_cast<double>(own_shape_count) / connects : 0.0);
    if (wire_types_.size() > 1) {
      // Sec. 3.2 / 3.6: wire types, nets routed with each, fast ones.
      std::vector<int> routed_by_type(wire_types_.size(), 0);
      std::vector<int> nets_by_type(wire_types_.size(), 0);
      int ndr_skipped = 0;
      for (const NetData& net : nets_) {
        nets_by_type[net.wire_type] += !net.skip;
        routed_by_type[net.wire_type] += net.routed;
        ndr_skipped += net.skip && net.net->hasNDR();
      }
      std::string types;
      for (size_t t = 0; t < wire_types_.size(); t++) {
        types += (t == 0 ? "" : ", ") + wire_types_[t].name + " "
                 + std::to_string(routed_by_type[t]) + "/"
                 + std::to_string(nets_by_type[t]) + " routed"
                 + (static_cast<int>(t) < fast_types_ ? " (fast)" : "");
      }
      logger_->report(
          "BonnRoute DR wire types: {}; {} fast; {} NDR nets skipped.",
          types,
          fast_types_,
          ndr_skipped);
    }
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
    if (params_.use_intervals) {
      std::array<int64_t, br::IntervalSearch::kNumIntervalEnds> ends{};
      int64_t vertices = 0;
      int64_t count = 0;
      for (const Worker& w : workers_) {
        if (!w.interval_search) {
          continue;
        }
        for (int e = 0; e < br::IntervalSearch::kNumIntervalEnds; e++) {
          ends[e] += w.interval_search->intervalEnds(
              static_cast<br::IntervalSearch::IntervalEnd>(e));
        }
        vertices += w.interval_search->intervalVertices();
        count += w.interval_search->intervals();
      }
      logger_->report(
          "BonnRoute DR intervals: {:.2f} vertices per interval; interval "
          "ends: track {}, start vertex {}, neighbour {}, wire edge {}, jog "
          "property {}, via property {}, future cost slope {}.",
          count > 0 ? static_cast<double>(vertices) / count : 0.0,
          ends[0],
          ends[1],
          ends[2],
          ends[3],
          ends[4],
          ends[5],
          ends[6]);
    }
    logger_->report(
        "BonnRoute DR partitions: {:.2f}s wall, {:.1f}% of {} threads busy "
        "({:.1f}% of the usable regions).",
        partition_wall_,
        partition_wall_ > 0
            ? 100.0 * partition_busy_ / (partition_wall_ * workers_.size())
            : 0.0,
        workers_.size(),
        partition_capacity_ > 0 ? 100.0 * partition_busy_ / partition_capacity_
                                : 0.0);
    if (params_.potential != 0) {
      double potential_seconds = 0;
      for (const Worker& w : workers_) {
        potential_seconds += w.potential_seconds;
      }
      logger_->report("BonnRoute DR: {:.2f}s thread time preparing potentials.",
                      potential_seconds);
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
      if (prepares[2] > 0) {
        double build_seconds = 0;
        double propagation_seconds = 0;
        int64_t columns = 0;
        int64_t rows = 0;
        int64_t elements = 0;
        int64_t pops = 0;
        for (const Worker& w : workers_) {
          build_seconds += w.potential.buildSeconds();
          propagation_seconds += w.potential.propagationSeconds();
          columns += w.potential.tileColumns();
          rows += w.potential.tileRows();
          elements += w.potential.generalCounters().boundary_elements;
          pops += w.potential.generalCounters().heap_pops;
        }
        const double n = static_cast<double>(prepares[2]);
        logger_->report(
            "BonnRoute DR: general potential: tile models {:.2f}s, "
            "propagation {:.2f}s thread time; per preparation {:.1f} x {:.1f} "
            "tiles, {:.0f} boundary elements, {:.0f} heap pops.",
            build_seconds,
            propagation_seconds,
            columns / n,
            rows / n,
            elements / n,
            pops / n);
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
    if (params_.spreading && extra_space_nets_ > 0) {
      logger_->report(
          "BonnRoute DR: wire spreading from global routing extra space: "
          "{} nets, {} net layer tiles.",
          extra_space_nets_,
          extra_space_tiles_);
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
    int64_t access_failures = 0;
    for (const Worker& w : workers_) {
      access_failures += w.access_failures;
      cleanup.add(w.cleanup);
      regenerations += w.regenerations;
      regenerated_paths += w.regenerated_paths;
    }
    if (params_.offtrack_access) {
      logger_->report(
          "BonnRoute DR: {} dynamic access path generations added {} paths; "
          "{} connections failed for a pin without feasible access option.",
          regenerations,
          regenerated_paths,
          access_failures);
    }
    logger_->report(
        "BonnRoute DR contact guard: {} nets with {} truncated pin contacts "
        "away from access points left to FlexDR; {} via pads on pins away "
        "from access points (not FlexDR terms, kept); {} via pads at access "
        "points marked connected; {} nets with wrong-way or off-track wires "
        "on on-track-only layers left to FlexDR.",
        guard_nets_,
        guard_contacts_,
        guard_via_contacts_,
        connected_vias_,
        wrong_way_nets_);
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
