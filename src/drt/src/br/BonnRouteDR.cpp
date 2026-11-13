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
#include <iterator>
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
#include "br/BrDrtRules.h"
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
#include "db/obj/frInstBlockage.h"
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
#include "gc/FlexGC.h"
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
    tail_lookahead = true;
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
  // Sec. 5.2 incremental routing of the later detailed_route -bonnroute
  // calls (the post-antenna ECO routes), in every profile.
  incremental = true;
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
  } else if (name == "critical_ripup_level") {
    critical_ripup_level = std::clamp<int>(v, 1, br::kMaxRemovableLevel);
  } else if (name == "first_pass_ripup_level") {
    first_pass_ripup_level = std::clamp<int>(v, 1, br::kMaxRemovableLevel);
  } else if (name == "region_ripup_effort") {
    region_ripup_effort = std::max(0, v);
  } else if (name == "keep_partial") {
    keep_partial = v != 0;
  } else if (name == "ripup_partitions") {
    ripup_partitions = std::max(0, v);
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
  } else if (name == "guide_dilation") {
    guide_dilation = std::max(0, v);
  } else if (name == "corridor_cost") {
    corridor_cost = std::max(0, v);
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
  } else if (name == "lef58_rules") {
    lef58_rules = v != 0;
  } else if (name == "same_net_vias") {
    same_net_vias = v != 0;
  } else if (name == "notch_fill") {
    notch_fill = v != 0;
  } else if (name == "final_ripup") {
    final_ripup = v != 0;
  } else if (name == "final_ripup_max_share") {
    final_ripup_max_share = std::clamp(v, 0, 100);
  } else if (name == "final_ripup_effort") {
    final_ripup_effort = std::max(0, v);
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
  } else if (name == "tail_batch_min") {
    tail_batch_min = std::max(0, v);
  } else if (name == "tail_batch_high") {
    tail_batch_high = std::clamp(v, 0, 100);
  } else if (name == "tail_batch_low") {
    tail_batch_low = std::clamp(v, 0, 100);
  } else if (name == "tail_threads") {
    tail_threads = std::max(0, v);
  } else if (name == "tail_order") {
    tail_order = v != 0;
  } else if (name == "tail_pipeline") {
    tail_pipeline = v != 0;
  } else if (name == "tail_lookahead") {
    tail_lookahead = v != 0;
  } else if (name == "tail_retry_delay") {
    tail_retry_delay = std::max(1, v);
  } else if (name == "tail_deferred_pass") {
    tail_deferred_pass = v != 0;
  } else if (name == "tail_deferred_settled") {
    tail_deferred_settled = std::max(0, v);
  } else if (name == "tail_defer_settled") {
    tail_defer_settled = std::max(0, v);
  } else if (name == "partitions") {
    partitions = v != 0;
  } else if (name == "partition_regions") {
    partition_regions = std::max(0, v);
  } else if (name == "tail_regions") {
    tail_regions = std::max(0, v);
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
  } else if (name == "offtrack_hanan_grid") {
    offtrack_hanan_grid = v != 0;
  } else if (name == "offtrack_on_track_ties") {
    offtrack_on_track_ties = v != 0;
  } else if (name == "offtrack_pin_classes") {
    offtrack_pin_classes = v != 0;
  } else if (name == "offtrack_access_stats") {
    offtrack_access_stats = v != 0;
  } else if (name == "offtrack_catalogue_sample") {
    offtrack_catalogue_sample = std::max(0, v);
  } else if (name == "offtrack_access_report") {
    offtrack_access_report = std::max(0, v);
  } else if (name == "offtrack_algorithm3_fallback") {
    offtrack_algorithm3_fallback = v != 0;
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
  } else if (name == "fast_grid_dense") {
    fast_grid_dense = v != 0;
  } else if (name == "fast_grid_memory_mb") {
    fast_grid_memory_mb = std::max(0, v);
  } else if (name == "cleanup_incremental_iterations") {
    cleanup_incremental_iterations = std::clamp(v, 0, 2);
  } else if (name == "cleanup_check_only_workers") {
    cleanup_check_only_workers = v != 0;
  } else if (name == "partition_check") {
    partition_check = v != 0;
  } else if (name == "start_check") {
    start_check = std::clamp(v, 0, 1);
  } else if (name == "planar_first") {
    planar_first = v != 0;
  } else if (name == "cleanup_normalize_initial_routing") {
    cleanup_normalize_initial_routing = v != 0;
  } else if (name == "left_net_track_assignment") {
    left_net_track_assignment = v != 0;
  } else if (name == "incremental") {
    incremental = v != 0;
  } else if (name == "incremental_resident") {
    incremental_resident = v != 0;
  } else if (name == "reserve_short_pitches") {
    reserve_short_pitches = std::max(0, v);
  } else if (name == "incremental_reservations") {
    incremental_reservations = v != 0;
  } else if (name == "incremental_potential") {
    incremental_potential = std::clamp(v, 0, 4);
  } else if (name == "incremental_ripup") {
    incremental_ripup = v != 0;
  } else if (name == "incremental_check") {
    incremental_check = std::clamp(v, 0, 2);
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
  // A routing space kept for incremental calls records instance shapes.
  record_fixed_ = params_.incremental && params_.incremental_resident;
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

void BonnRouteDR::buildRoutingSpace()
{
  frTechObject* tech = design_->getTech();
  const int num_z = graph_.numLayers();
  // Shape-grid layers: wiring layer z is 2z, the via layer above it 2z + 1.
  std::vector<br::LayerRules> rules;
  for (int z = 0; z < num_z; z++) {
    const LayerInfo& info = layers_[z];
    // Spacing and line-end rules (Sec. 3.1) from the drt tech (BrDrtRules).
    rules.push_back(br::drtWiringRules(info.layer, params_.lef58_rules));
    if (z + 1 == num_z) {
      continue;
    }
    // Cut spacing and the inter-layer via rule with the cut layer below
    // (Sec. 3.1, 3.2).
    rules.push_back(
        br::drtViaRules(tech->getLayer(info.fr + 1),
                        z > 0 ? tech->getLayer(layers_[z - 1].fr + 1) : nullptr,
                        info.cut_spacing,
                        info.horizontal,
                        params_.lef58_rules));
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
      low = br::drtWiringRules(layer, params_.lef58_rules);
    } else if (fr == bottom - 1
               && layer->getType() == odb::dbTechLayerType::CUT) {
      frCoord s = layer->hasCutSpacing() ? layer->getCutSpacingValue() : 0;
      if (s <= 0 && layer->getDefaultViaDef() != nullptr) {
        const frVia via(layer->getDefaultViaDef());
        s = via.getCutBBox().minDXDY();
      }
      low = br::drtViaRules(
          layer, nullptr, s, layers_.front().horizontal, params_.lef58_rules);
    } else {
      continue;
    }
    grid_of_fr_[fr] = static_cast<int>(rules.size());
    rules.push_back(low);
  }
  // Side rules (Lef58Rules.md): preferred wires carry no line-end sides and
  // corners when every wire type's via pads on the layer cover their ends.
  for (int z = 0; z < num_z; z++) {
    br::LayerRules& r = rules[br::wiringGridLayer(z)];
    if (!r.sideRules()) {
      continue;
    }
    bool covered = true;
    for (size_t k = 0; k <= std::min(ndrs.size(), kMaxNdrs); k++) {
      frNonDefaultRule* ndr = k == 0 ? nullptr : ndrs[k - 1];
      frCoord w = layers_[z].width;
      if (ndr != nullptr) {
        w = std::max(w, ndr->getWidth(layers_[z].fr / 2 - 1));
      }
      auto via_of = [&](const int vz) {
        const frViaDef* def = layers_[vz].up_via;
        if (ndr != nullptr
            && ndr->getPrefVia(layers_[vz].fr / 2 - 1) != nullptr) {
          def = ndr->getPrefVia(layers_[vz].fr / 2 - 1);
        }
        return def;
      };
      std::vector<br::Rect> pads;
      if (z + 1 < num_z) {
        pads.push_back(brRect(frVia(via_of(z)).getLayer1BBox()));
      }
      if (z > 0) {
        pads.push_back(brRect(frVia(via_of(z - 1)).getLayer2BBox()));
      }
      covered = covered && br::padsCoverWireEnd(w, pads);
    }
    r.wire_ends_covered = covered;
  }
  if (params_.lef58_rules) {
    for (int g = 0; g < static_cast<int>(rules.size()); g++) {
      const br::LayerRules& r = rules[g];
      if (!r.sideRules() && r.cut_classes.empty()) {
        continue;
      }
      std::string text;
      for (const br::SideRule& rule : r.side_rules) {
        static const char* kKinds[]
            = {"keepout", "keepout-corner", "eol", "end-to-end"};
        text += fmt::format(" {}(eol<{} fwd {} back {} side {})",
                            kKinds[static_cast<int>(rule.kind)],
                            r.eol_widths[rule.eol],
                            rule.forward,
                            rule.backward,
                            rule.side);
      }
      if (r.corner_spacing > 0) {
        text += fmt::format(" corner {}", r.corner_spacing);
      }
      if (r.sideRules()) {
        text += r.wire_ends_covered ? " wire-ends-covered" : " wire-ends-real";
      }
      if (!r.cut_classes.empty()) {
        text += fmt::format(" cut classes {} spacing", r.cut_classes.size());
        for (const br::Coord c : r.cut_class_spacing) {
          text += fmt::format(" {}", c);
        }
      }
      frLayerNum fr = -1;
      for (int l = 0; l < static_cast<int>(grid_of_fr_.size()); l++) {
        if (grid_of_fr_[l] == g) {
          fr = l;
        }
      }
      logger_->report("BonnRoute DR LEF58 rules {}:{}",
                      fr >= 0 ? tech->getLayer(fr)->getName() : "?",
                      text);
    }
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
  need_by_grid_layer_.assign(shape_classes_.numLayers(), 0);
  narrow_need_by_grid_layer_.assign(shape_classes_.numLayers(), 0);
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
    need_by_grid_layer_[g] = int64_t{shape_classes_.ruleReach(g)} + model
                             + std::max(geom.cross_cell, geom.along_cell);
    narrow_need_by_grid_layer_[g]
        = int64_t{shape_classes_.narrowRuleReach(g)} + model
          + std::max(geom.cross_cell, geom.along_cell);
    need = std::max(need, need_by_grid_layer_[g]);
  }
  shape_grid_.init(&shape_classes_,
                   geometries,
                   br::ChunkAxis{tile_x0_, tile_dx_, tiles_x_},
                   br::ChunkAxis{tile_y0_, tile_dy_, tiles_y_});
  checker_.init(&shape_grid_);
  // Sec. 3.6: the fast grid holds the most frequent wire types.
  grid_.init(&graph_, fast_types_);
  if (params_.fast_grid_dense) {
    // Dense layers bottom-up while their arrays fit the budget.
    const size_t budget
        = static_cast<size_t>(params_.fast_grid_memory_mb) * (size_t{1} << 20);
    size_t used = 0;
    std::string dense_layers;
    std::string run_layers;
    for (int z = 0; z < graph_.numLayers(); z++) {
      const size_t b = grid_.denseBytes(z);
      std::string& list = used + b <= budget ? dense_layers : run_layers;
      if (used + b <= budget) {
        grid_.setDense(z, true);
        used += b;
      }
      list += (list.empty() ? "" : " ") + std::to_string(z);
    }
    logger_->report(
        "BonnRoute DR fast grid dense: track layers [{}] dense ({:.1f} MB of "
        "{} MB), layers [{}] in runs.",
        dense_layers,
        used / 1e6,
        params_.fast_grid_memory_mb,
        run_layers);
  }
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
  // Per net halos follow in buildPartitionHalos(), once tau is known.
  area_halo_tiles_.assign(num_z, halo_tiles_);
  access_halo_tiles_.assign(num_z, halo_tiles_);
  narrow_area_halo_tiles_.assign(num_z, halo_tiles_);
  narrow_access_halo_tiles_.assign(num_z, halo_tiles_);
  // Search label budget: a label count, search_label_mb MiB of 32-byte
  // labels (the vertex-search label size when the budget was introduced;
  // the count is kept so that searches stop at the same label).
  label_limit_ = static_cast<size_t>(std::max(0, params_.search_label_mb))
                 * (size_t{1} << 20) / 32;
}

void BonnRouteDR::buildPartitionHalos()
{
  const int num_z = graph_.numLayers();
  const int64_t tile = std::max<int64_t>(1, std::min(tile_dx_, tile_dy_));
  auto tiles = [tile](const int64_t distance) {
    return static_cast<int>(
        std::max<int64_t>(1, (distance + tile - 1) / tile));
  };
  // Access paths start within this many pitches of the pin (catalogue
  // radius, stub jogs); dynamic generation looks two pitches farther and
  // collects obstacles one pitch plus 2 tau beyond its targets.
  const int64_t radius
      = std::max(params_.offtrack_radius, params_.max_access_jog);
  int64_t reach = 0;
  int64_t narrow_reach = 0;
  frCoord pitch = 0;
  int64_t tau = 0;
  for (int z = 0; z < num_z; z++) {
    // Grid layers of a net up to track layer z: wires up to z, vias up to
    // the one between z - 1 and z, whose cut projects into the cut layer
    // above z; the low layers hold pins and access vias.
    for (int g = 0; g < shape_classes_.numLayers(); g++) {
      if (g >= track_grid_layers_ || g <= 2 * z + 1) {
        reach = std::max(reach, need_by_grid_layer_[g]);
        // The narrow distance holds where the fast grid's models are narrow
        // (FastGridBuilder::modelReach); the low layers have no models.
        const bool narrow_models
            = g >= track_grid_layers_ || !fast_builder_.wideModels(g);
        narrow_reach = std::max(narrow_reach,
                                narrow_models ? narrow_need_by_grid_layer_[g]
                                              : need_by_grid_layer_[g]);
      }
    }
    // Pitches of the jog pairs and edges the updates recompute: layer z and
    // its positions, the tracks of layer z + 1.
    pitch = std::max(pitch, layers_[std::min(z + 1, num_z - 1)].pitch);
    pitch = std::max(pitch, layers_[z].pitch);
    tau = std::max<int64_t>(tau, z < static_cast<int>(tau_.size()) ? tau_[z] : 1);
    // A shape change writes up to reach + one pitch beyond the shape and
    // reads up to 2 * reach + one pitch; a query reads up to reach.
    area_halo_tiles_[z] = tiles(2 * reach + pitch);
    narrow_area_halo_tiles_[z] = tiles(2 * narrow_reach + pitch);
    // Around a pin. Its access paths, also the dynamically generated ones,
    // lie within (radius + 3) pitches + 2 tau of it (targets up to radius +
    // 2 pitches, paths inside the window one pitch + 2 tau beyond them) and
    // are committed and released like any shape: 2 * distance + pitch
    // beyond them. Dynamic generation collects obstacles up to the layer's
    // whole rule distance beyond that window (Legality::blockages).
    const int64_t paths = (radius + 3) * pitch + 2 * tau;
    access_halo_tiles_[z] = tiles(paths + 2 * reach + pitch);
    narrow_access_halo_tiles_[z]
        = tiles(std::max(paths + 2 * narrow_reach + pitch, paths + reach));
  }
}

void BonnRouteDR::markWideNets()
{
  // Nets whose own shapes can be wide: nondefault wire types; pins and
  // reservations (the overlay of the net's own shapes reads up to
  // classReach() around a wide one); access path shapes, which a partition
  // thread commits and releases.
  const int threads = static_cast<int>(workers_.size());
#pragma omp parallel for num_threads(threads) schedule(dynamic, 256)
  for (int k = 0; k < static_cast<int>(nets_.size()); k++) {
    nets_[k].wide_shapes = netHasWideShapes(nets_[k]);
  }
}

bool BonnRouteDR::netHasWideShapes(const NetData& net) const
{
  {
    bool wide = net.wire_type != 0;
    auto check = [&](const RShape& shape) {
      const int g = gridLayer(shape.layer);
      if (g < 0) {
        return;
      }
      wide = wide || shape_classes_.wide(shapeClass(shape, g));
      const int above = projectionLayer(g);
      if (above >= 0) {
        RShape projection = shape;
        projection.role = br::ShapeRole::kCutProjection;
        wide = wide || shape_classes_.wide(shapeClass(projection, above));
      }
    };
    for (const RShape& pin : net.pins) {
      check(pin);
    }
    for (const RShape& shape : net.input_reservation) {
      check(shape);
    }
    std::vector<RShape> shapes;
    for (const Term& term : net.terms) {
      for (const RShape& shape : term.reservation) {
        check(shape);
      }
      for (size_t o = 0; o < term.options.size() && !wide; o++) {
        accessShapes(term.options[o], shapes);
        for (const RShape& shape : shapes) {
          check(shape);
        }
      }
    }
    return wide;
  }
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
                 || role == br::ShapeRole::kRealWire
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
    // Sec. 5.1: the halo of a net without wide shapes covers narrow classes.
    if (local != nullptr && owner >= 0 && !nets_[owner].wide_shapes
        && shape_classes_.wide(cls)) {
      unexpected_wide_updates_++;
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
  // D1 S2: the shapes of the connections the net closed so far.
  if (net.partial) {
    for (const RShape& s : net.shapes) {
      add(s, routedType(s.role));
    }
  }
  // Incremental mode: kept input routing is a component of the net.
  if (net.input_present) {
    for (const RShape& s : net.input_shapes) {
      add(s, routedType(s.role));
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
    const frBlockObject* obj = nullptr;
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
      fixed[fr].push_back({rect, owner, obj});
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
      // Incremental mode: the fixed shapes of each instance.
      if (record_fixed_) {
        recordFixedShape(
            f.obj,
            {static_cast<frLayerNum>(fr), f.rect},
            f.owner >= 0 ? br::ShapeType::kPin : br::ShapeType::kBlockage,
            f.owner);
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
  br::ShapeRole role = preferred ? br::ShapeRole::kPreferredWire
                                 : br::ShapeRole::kNonPreferredWire;
  if (wire.patch && preferred
      && shape_classes_.layer(br::wiringGridLayer(z)).sideRules()) {
    // Side rules: a patch ends where no via pad covers its end; its ends
    // and corners are real (Lef58Rules.md).
    role = br::ShapeRole::kRealWire;
  }
  return {wire.layer, r, role, wire.wire_type};
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
      if (shape.role == br::ShapeRole::kPreferredWire
          && dr_.shape_classes_.layer(br::wiringGridLayer(z)).sideRules()) {
        // Side rules: an off-track piece or patch ends where no via pad
        // covers its end (Lef58Rules.md).
        shape.role = br::ShapeRole::kRealWire;
      }
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
      if (!req.overlap && !req.line_end && !req.sides
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
      if (req.sides) {
        // Side and corner rules: reference points where a rectangle of
        // the model meets one of the piece (closed, Lef58Rules.md).
        thread_local std::vector<br::SideTest> tests;
        tests.clear();
        br::appendSideTests(classes,
                            br::PlacedShape{m, mc, 3},
                            br::PlacedShape{p.rect, p.cls, p.line_end_ends},
                            tests);
        for (const br::SideTest& t : tests) {
          out.push_back({t.from_b.xlo - t.from_a.xhi,
                         t.from_b.ylo - t.from_a.yhi,
                         t.from_b.xhi - t.from_a.xlo,
                         t.from_b.yhi - t.from_a.ylo});
        }
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
  offgrid_rejected += o.offgrid_rejected;
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
  fallback_instances += o.fallback_instances;
  fallback_improved += o.fallback_improved;
  pin_classes += o.pin_classes;
  mapped_options += o.mapped_options;
  mapped_dropped += o.mapped_dropped;
  sample_pins += o.sample_pins;
  sample_equal += o.sample_equal;
  sample_mapped_only += o.sample_mapped_only;
  sample_own_only += o.sample_own_only;
  sample_same_environment += o.sample_same_environment;
  sample_same_environment_translated += o.sample_same_environment_translated;
  sample_same_environment_translated_unequal
      += o.sample_same_environment_translated_unequal;
  grid_lines += o.grid_lines;
  grid_lines_other += o.grid_lines_other;
  grid_vertices_other += o.grid_vertices_other;
  search_states += o.search_states;
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

const frAccessPoint* BonnRouteDR::accessPointAt(const frBlockObject* term,
                                                const frLayerNum layer,
                                                const odb::Point& p) const
{
  if (term->typeId() == frcInstTerm) {
    const auto* inst_term = static_cast<const frInstTerm*>(term);
    frInst* inst = inst_term->getInst();
    const int pa_idx = inst->getPinAccessIdx();
    const odb::dbTransform xform = inst->getNoRotationTransform();
    for (const auto& pin : inst_term->getTerm()->getPins()) {
      if (!pin->hasPinAccess() || pa_idx < 0
          || pa_idx >= pin->getNumPinAccess()) {
        continue;
      }
      for (const auto& ap : pin->getPinAccess(pa_idx)->getAccessPoints()) {
        odb::Point q = ap->getPoint();
        xform.apply(q);
        if (q == p && ap->getLayerNum() == layer) {
          return ap.get();
        }
      }
    }
  } else if (term->typeId() == frcBTerm) {
    for (const auto& pin : static_cast<const frBTerm*>(term)->getPins()) {
      if (!pin->hasPinAccess()) {
        continue;
      }
      for (const auto& ap : pin->getPinAccess(0)->getAccessPoints()) {
        if (ap->getPoint() == p && ap->getLayerNum() == layer) {
          return ap.get();
        }
      }
    }
  }
  return nullptr;
}

bool BonnRouteDR::viaCheckClean(frBlockObject* term,
                                frNet* net,
                                const frViaDef* def,
                                const odb::Point& p,
                                const bool own_net_only) const
{
  FlexGCWorker gc(design_->getTech(), logger_, router_cfg_);
  gc.setIgnoreMinArea();
  gc.setIgnoreLongSideEOL();
  gc.setIgnoreCornerSpacing();
  const frCoord pitch
      = design_->getTech()->getLayer(def->getLayer1Num())->getPitch();
  odb::Rect ext_box;
  odb::Rect(p, p).bloat(5 * pitch, ext_box);
  gc.setExtBox(ext_box);
  gc.setDrcBox(ext_box);
  if (!net->hasNDR()
      || net->autoTaperEnabled(router_cfg_->AUTO_TAPER_NDR_NETS)) {
    if (term->typeId() == frcInstTerm) {
      gc.addTargetObj(static_cast<frInstTerm*>(term)->getInst());
    } else {
      gc.addTargetObj(term);
    }
  }
  gc.initPA0(design_);
  frVia via(def, p);
  gc.addPAObj(&via, net);
  gc.initPA1();
  gc.main();
  for (const auto& marker : gc.getMarkers()) {
    if (!own_net_only) {
      return false;
    }
    const auto& srcs = marker->getSrcs();
    if (!srcs.empty()
        && std::ranges::all_of(
            srcs, [net](const frBlockObject* s) { return s == net; })) {
      return false;
    }
  }
  return true;
}

void BonnRouteDR::buildPinViaLegality()
{
  // Every lower vertex of a search via whose bottom pad overlaps a pin shape
  // of the net. Approved as is: an access point whose FlexPA via definitions
  // include the via. Otherwise the via is checked once per instance class
  // (master, orientation, offset in master coordinates, definition) on the
  // first instance of the class in net order; same-net markers only, so the
  // verdict depends on the class alone (different-net spacing is the routing
  // space's).
  const double start = omp_get_wtime();
  struct ClassKey
  {
    const void* owner;
    int orient;
    frCoord x, y;
    const frViaDef* def;
    bool operator==(const ClassKey&) const = default;
  };
  struct ClassKeyHash
  {
    size_t operator()(const ClassKey& k) const
    {
      size_t h = std::hash<const void*>()(k.owner);
      for (const size_t v : {static_cast<size_t>(k.orient),
                             static_cast<size_t>(k.x),
                             static_cast<size_t>(k.y),
                             std::hash<const void*>()(k.def)}) {
        h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
      }
      return h;
    }
  };
  struct ClassCheck
  {
    frBlockObject* term;
    frNet* net;
    const frViaDef* def;
    odb::Point point;
    bool clean = true;
  };
  struct Position
  {
    int net;
    size_t vertex;
    int check;
  };
  std::unordered_map<ClassKey, int, ClassKeyHash> class_index;
  std::vector<ClassCheck> checks;
  std::vector<Position> positions;
  int64_t approved = 0;
  const int num_fr = static_cast<int>(z_of_fr_.size());
  std::vector<frRect> shapes;
  for (int n = 0; n < static_cast<int>(nets_.size()); n++) {
    NetData& data = nets_[n];
    if (data.skip) {
      continue;
    }
    const br::WireType& type = wire_types_[data.wire_type];
    auto visit = [&](frBlockObject* term) {
      shapes.clear();
      if (term->typeId() == frcInstTerm) {
        static_cast<frInstTerm*>(term)->getShapes(shapes);
      } else {
        static_cast<frBTerm*>(term)->getShapes(shapes);
      }
      for (const frRect& shape : shapes) {
        const frLayerNum fr = shape.getLayerNum();
        const int z = fr >= 0 && fr < num_fr ? z_of_fr_[fr] : -1;
        if (z < 0 || z + 1 >= graph_.numLayers()
            || z >= static_cast<int>(type.via.size()) || !type.via[z].valid
            || z >= static_cast<int>(type_vias_[data.wire_type].size())
            || type_vias_[data.wire_type][z] == nullptr) {
          continue;
        }
        const frViaDef* def = type_vias_[data.wire_type][z];
        const br::Rect pad = type.via[z].bottom.rect;
        const odb::Rect r = shape.getBBox();
        const br::TrackLayer& layer = graph_.layers[z];
        // Vertices v with (v + pad) overlapping r with positive area.
        const br::Coord c_lo = layer.horizontal ? r.yMin() - pad.yhi
                                                : r.xMin() - pad.xhi;
        const br::Coord c_hi = layer.horizontal ? r.yMax() - pad.ylo
                                                : r.xMax() - pad.xlo;
        const br::Coord a_lo = layer.horizontal ? r.xMin() - pad.xhi
                                                : r.yMin() - pad.yhi;
        const br::Coord a_hi = layer.horizontal ? r.xMax() - pad.xlo
                                                : r.yMax() - pad.ylo;
        const int i0 = static_cast<int>(
            std::ranges::upper_bound(layer.tracks, c_lo) - layer.tracks.begin());
        const int i1 = static_cast<int>(
            std::ranges::lower_bound(layer.tracks, c_hi) - layer.tracks.begin());
        const int j0 = static_cast<int>(
            std::ranges::upper_bound(layer.positions, a_lo)
            - layer.positions.begin());
        const int j1 = static_cast<int>(
            std::ranges::lower_bound(layer.positions, a_hi)
            - layer.positions.begin());
        for (int i = i0; i < i1; i++) {
          for (int j = j0; j < j1; j++) {
            const auto [x, y] = layer.point(i, j);
            const odb::Point p(x, y);
            const size_t v = graph_.vertex(z, i, j);
            if (const frAccessPoint* ap = accessPointAt(term, fr, p)) {
              bool listed = false;
              for (int cuts = 1; ap->hasViaDef(cuts) && !listed; cuts++) {
                listed = std::ranges::find(ap->getViaDefs(cuts), def)
                         != ap->getViaDefs(cuts).end();
              }
              if (listed) {
                approved++;
                continue;
              }
            }
            ClassKey key{term, 0, p.x(), p.y(), def};
            if (term->typeId() == frcInstTerm) {
              frInst* inst = static_cast<frInstTerm*>(term)->getInst();
              odb::dbTransform xform = inst->getTransform();
              xform.invert();
              odb::Point q = p;
              xform.apply(q);
              key = {inst->getMaster(),
                     static_cast<int>(inst->getOrient().getValue()),
                     q.x(),
                     q.y(),
                     def};
            }
            auto [it, inserted]
                = class_index.try_emplace(key, static_cast<int>(checks.size()));
            if (inserted) {
              checks.push_back({term, data.net, def, p});
            }
            positions.push_back({n, v, it->second});
          }
        }
      }
    };
    for (frInstTerm* term : data.net->getInstTerms()) {
      visit(term);
    }
    for (frBTerm* term : data.net->getBTerms()) {
      visit(term);
    }
  }
#pragma omp parallel for num_threads(router_cfg_->MAX_THREADS) \
    schedule(dynamic, 64)
  for (int k = 0; k < static_cast<int>(checks.size()); k++) {
    ClassCheck& c = checks[k];
    c.clean = viaCheckClean(c.term, c.net, c.def, c.point, true);
  }
  int64_t forbidden = 0;
  int64_t forbidden_nets = 0;
  int64_t unclean_classes = 0;
  for (const ClassCheck& c : checks) {
    unclean_classes += !c.clean;
  }
  for (const Position& pos : positions) {
    if (!checks[pos.check].clean) {
      nets_[pos.net].forbidden_pin_vias.push_back(pos.vertex);
    }
  }
  for (NetData& data : nets_) {
    std::vector<size_t>& f = data.forbidden_pin_vias;
    std::ranges::sort(f);
    f.erase(std::unique(f.begin(), f.end()), f.end());
    forbidden += static_cast<int64_t>(f.size());
    forbidden_nets += !f.empty();
  }
  if (params_.verbose) {
    logger_->report(
        "BonnRoute DR planar-first rule: {} via positions on own pins ({} at "
        "approved access points), {} instance-class checks ({} unclean), {} "
        "forbidden positions on {} nets; {:.2f}s.",
        positions.size() + approved,
        approved,
        checks.size(),
        unclean_classes,
        forbidden,
        forbidden_nets,
        omp_get_wtime() - start);
  }
}

void BonnRouteDR::checkPathStarts()
{
  // Every via whose bottom pad lies on a pin shape of its own net, checked
  // against the fixed design (pins, obstructions, special nets) as FlexPA
  // checks an access via (FlexPA::isViaViolationFree): FlexGC with min area,
  // long-side EOL and corner spacing ignored, ext and DRC box 5 pitches
  // around the via, the instance (or top-level term) as target.
  struct StartCheck
  {
    const NetData* net;
    const frViaDef* def;
    odb::Point point;
    frBlockObject* term;
    const frAccessPoint* ap;
    bool clean = true;
    int alternatives = 0;
    const frViaDef* clean_alternative = nullptr;
    // Other one-cut definitions of the cut layer (technology).
    int technology_alternatives = 0;
    const frViaDef* clean_technology_alternative = nullptr;
  };
  const double start = omp_get_wtime();
  std::vector<StartCheck> checks;
  frRegionQuery* rq = design_->getRegionQuery();
  for (const NetData& data : nets_) {
    if (!data.routed) {
      continue;
    }
    for (const OutVia& v : data.vias) {
      const frLayerNum layer = v.def->getLayer1Num();
      frBlockObject* term = nullptr;
      forEachNetPinAt(
          rq, data.net, layer, v.point, [&](const frBlockObject* obj, const odb::Rect&) {
            if (term == nullptr) {
              term = const_cast<frBlockObject*>(obj);
            }
          });
      if (term == nullptr) {
        continue;
      }
      checks.push_back(
          {&data, v.def, v.point, term, accessPointAt(term, layer, v.point)});
    }
  }
  auto via_clean = [this](const StartCheck& c, const frViaDef* def) {
    return viaCheckClean(c.term, c.net->net, def, c.point, false);
  };
#pragma omp parallel for num_threads(router_cfg_->MAX_THREADS) \
    schedule(dynamic, 64)
  for (int i = 0; i < static_cast<int>(checks.size()); i++) {
    StartCheck& c = checks[i];
    c.clean = via_clean(c, c.def);
    if (c.clean) {
      continue;
    }
    for (const frViaDef* def :
         design_->getTech()->getLayer(c.def->getCutLayerNum())->getViaDefs()) {
      if (def == c.def || def->getNumCut() != 1) {
        continue;
      }
      c.technology_alternatives++;
      if (c.clean_technology_alternative == nullptr && via_clean(c, def)) {
        c.clean_technology_alternative = def;
      }
    }
    if (c.ap == nullptr) {
      continue;
    }
    for (int cuts = 1; c.ap->hasViaDef(cuts); cuts++) {
      for (const frViaDef* def : c.ap->getViaDefs(cuts)) {
        if (def == c.def) {
          continue;
        }
        c.alternatives++;
        if (c.clean_alternative == nullptr && via_clean(c, def)) {
          c.clean_alternative = def;
        }
      }
    }
  }
  int64_t at_access_point = 0;
  int64_t unclean = 0;
  int64_t unclean_at_access_point = 0;
  int64_t unclean_first_definition = 0;
  int64_t with_clean_alternative = 0;
  int64_t without_alternative = 0;
  int64_t with_clean_technology_alternative = 0;
  for (const StartCheck& c : checks) {
    at_access_point += c.ap != nullptr;
    if (c.clean) {
      continue;
    }
    unclean++;
    if (c.ap != nullptr) {
      unclean_at_access_point++;
      unclean_first_definition
          += c.ap->hasViaDef() && c.ap->getViaDefs(1).front() == c.def;
    }
    with_clean_alternative += c.clean_alternative != nullptr;
    without_alternative += c.alternatives == 0;
    with_clean_technology_alternative
        += c.clean_technology_alternative != nullptr;
    logger_->report(
        "BonnRoute DR start check unclean: net {} via {} at ( {} {} ) access "
        "point {} alternatives {} clean alternative {} technology "
        "alternatives {} clean technology alternative {}",
        c.net->net->getName(),
        c.def->getName(),
        c.point.x(),
        c.point.y(),
        c.ap != nullptr ? 1 : 0,
        c.alternatives,
        c.clean_alternative != nullptr ? c.clean_alternative->getName()
                                       : std::string("-"),
        c.technology_alternatives,
        c.clean_technology_alternative != nullptr
            ? c.clean_technology_alternative->getName()
            : std::string("-"));
  }
  logger_->report(
      "BonnRoute DR start check: {} vias on own pins checked ({} at access "
      "points); {} unclean ({} at access points, {} with the access point's "
      "first one-cut definition); {} with a clean alternative definition, {} "
      "without an alternative definition; {} with a clean one-cut definition "
      "of the cut layer; {:.2f}s.",
      checks.size(),
      at_access_point,
      unclean,
      unclean_at_access_point,
      unclean_first_definition,
      with_clean_alternative,
      without_alternative,
      with_clean_technology_alternative,
      omp_get_wtime() - start);
}

int BonnRouteDR::pinsWithoutContact(const NetData& net,
                                    const std::vector<OutWire>& wires) const
{
  // FlexDR's pin contacts (dr/FlexDR_conn.cpp buildPin2epMap, and the
  // worker's initNets_searchRepair_pin2epMap): every pin of the net with a
  // shape at a truncated path-segment end, or at a via marked bottom / top
  // connected, on that layer. writeNets truncates wire ends at the net's
  // pin points and marks vias at access points connected.
  std::vector<const frBlockObject*> contacted;
  frRegionQuery* rq = design_->getRegionQuery();
  auto add_pins = [&](const frLayerNum layer, const odb::Point& p) {
    forEachNetPinAt(
        rq, net.net, layer, p, [&](const frBlockObject* obj, const odb::Rect&) {
          contacted.push_back(obj);
        });
  };
  for (const OutWire& w : wires) {
    if (w.begin == w.end) {
      continue;
    }
    for (const auto& [layer, p] : net.pin_points) {
      if (layer == w.layer && (p == w.begin || p == w.end)) {
        add_pins(layer, p);
      }
    }
  }
  for (const OutVia& v : net.vias) {
    for (const frLayerNum layer :
         {v.def->getLayer1Num(), v.def->getLayer2Num()}) {
      if (atAccessPoint(net, layer, v.point)) {
        add_pins(layer, v.point);
      }
    }
  }
  if (net.keeps_input) {
    // Incremental calls (Fig. 9): the kept pin-to-pin paths stay in the
    // frNet and carry the contacts of the pins they join, by the same rule.
    for (const auto& shape : net.net->getShapes()) {
      if (shape->typeId() != frcPathSeg) {
        continue;
      }
      const auto* seg = static_cast<const frPathSeg*>(shape.get());
      const auto [b, e] = seg->getPoints();
      if (seg->getBeginStyle() == frcTruncateEndStyle) {
        add_pins(seg->getLayerNum(), b);
      }
      if (seg->getEndStyle() == frcTruncateEndStyle) {
        add_pins(seg->getLayerNum(), e);
      }
    }
    for (const auto& via : net.net->getVias()) {
      const frViaDef* def = via->getViaDef();
      const odb::Point o = via->getOrigin();
      if (via->isBottomConnected()
          || atAccessPoint(net, def->getLayer1Num(), o)) {
        add_pins(def->getLayer1Num(), o);
      }
      if (via->isTopConnected() || atAccessPoint(net, def->getLayer2Num(), o)) {
        add_pins(def->getLayer2Num(), o);
      }
    }
  }
  int missing = 0;
  auto check = [&](const frBlockObject* pin) {
    if (std::ranges::find(contacted, pin) == contacted.end()) {
      missing++;
    }
  };
  for (const frInstTerm* it : net.net->getInstTerms()) {
    check(it);
  }
  for (const frBTerm* bt : net.net->getBTerms()) {
    check(bt);
  }
  return missing;
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
                                    const int extra_radius,
                                    const bool algorithm3) const
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
    // Layers without a track layer (the pin layer below the routing
    // layers) take the manufacturing grid: their lines at multiples of tau
    // are also lines of the adjacent layers (Algorithm 3), and tau = 1 put
    // wires of those layers off the grid.
    grid_in.layers.push_back({z >= 0 ? tau_[z] : mfg,
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
  // Every line is on the manufacturing grid (see snap_out); with tau at most
  // the grid unit the Hanan lines suffice (D9). A layer whose wires stay on
  // tracks needs its tracks inside the window as lines: a path may run
  // along a track that holds no target.
  // `algorithm3`: the grid, states and ties of offtrack_hanan_grid 0.
  const bool hanan = params_.offtrack_hanan_grid && !algorithm3;
  grid_in.grid_unit = hanan ? mfg : 0;
  if (hanan) {
    for (const frLayerNum fr : frs) {
      const int z = z_of_fr_[fr];
      if (z < 0 || !onTracksOnly(fr)) {
        continue;
      }
      const br::TrackLayer& l = graph_.layers[z];
      const br::GeoCoord lo = l.horizontal ? window.ylo : window.xlo;
      const br::GeoCoord hi = l.horizontal ? window.yhi : window.xhi;
      for (auto it = std::ranges::lower_bound(l.tracks, lo);
           it != l.tracks.end() && *it <= hi;
           ++it) {
        grid_in.borders.push_back(
            l.horizontal ? br::Box{window.xlo, *it, window.xhi, *it}
                         : br::Box{*it, window.ylo, *it, window.yhi});
      }
    }
  }
  stats.blockage_seconds += step_timer.elapsed();
  step_timer.reset();
  const std::vector<br::BlockageGridLayer> grid
      = br::buildBlockageGrid(grid_in);
  stats.grid_seconds += step_timer.elapsed();
  if (params_.offtrack_access_stats) {
    // M3 (c): lines and vertices inside the window, of this grid and of the
    // other one (Hanan lines / Algorithm 3).
    auto count = [&window](const std::vector<br::BlockageGridLayer>& layers,
                           int64_t& lines,
                           int64_t& vertices) {
      for (const br::BlockageGridLayer& l : layers) {
        const int64_t nx = std::ranges::count_if(l.xs, [&](br::GeoCoord x) {
          return x >= window.xlo && x <= window.xhi;
        });
        const int64_t ny = std::ranges::count_if(l.ys, [&](br::GeoCoord y) {
          return y >= window.ylo && y <= window.yhi;
        });
        lines += nx + ny;
        vertices += nx * ny;
      }
    };
    int64_t unused = 0;
    count(grid, stats.grid_lines, unused);
    br::BlockageGridInput other = grid_in;
    other.grid_unit = grid_in.grid_unit > 0 ? 0 : mfg;
    count(br::buildBlockageGrid(other),
          stats.grid_lines_other,
          stats.grid_vertices_other);
  }
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
  req.merge_direction_states = hanan;
  req.on_track_ties = params_.offtrack_on_track_ties && !algorithm3;
  for (size_t s = 0; s < frs.size(); s++) {
    const int z = z_of_fr_[frs[s]];
    if (z >= 0) {
      req.layers[s].tie_tracks = &graph_.layers[z].tracks;
      req.layers[s].horizontal = graph_.layers[z].horizontal;
    }
  }
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
  stats.search_states += res.states;
  stats.settled += res.settled;
  step_timer.reset();
  for (const br::OffTrackPath& path : res.paths) {
    const Target& target = targets[path.target_id];
    // Every node is on the manufacturing grid (all blockage-grid lines are);
    // a path that is not would be an Off Grid violation for drt.
    if (std::ranges::any_of(path.nodes, [mfg](const br::RouteNode& n) {
          return n.p.x % mfg != 0 || n.p.y % mfg != 0;
        })) {
      stats.offgrid_rejected++;
      continue;
    }
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
  // Signature components left out by the M3 ablation (bits): pin net
  // groups, preferred access point marks, shapes not of the instance,
  // tracks.
  enum SignatureAblation
  {
    kAblateNetGroups = 1,
    kAblatePreferred = 2,
    kAblateForeignShapes = 4,
    kAblateTracks = 8
  };
  auto owned_by = [](const frBlockObject* obj, const frInst* inst) {
    if (obj->typeId() == frcInstTerm) {
      return static_cast<const frInstTerm*>(obj)->getInst() == inst;
    }
    if (obj->typeId() == frcInstBlockage) {
      return static_cast<const frInstBlockage*>(obj)->getInst() == inst;
    }
    return false;
  };
  auto instance_signature = [&](const InstanceAccess& ia,
                                const int k,
                                const int ablation) {
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
      if (ablation & kAblateNetGroups) {
        g = 1;
      }
      groups.push_back(g);
      sig.addPin(ia.pin_index[p], g);
      if (params_.offtrack_from_access_points) {
        // The access points the catalogue paths start at (label -1 - p).
        std::vector<AccessSource> aps;
        accessSources(nets_[ia.terms[p].first].terms[ia.terms[p].second], aps);
        for (const AccessSource& a : aps) {
          const bool preferred = a.preferred || (ablation & kAblatePreferred);
          sig.addShape(a.layer,
                       {a.point.x(), a.point.y(), a.point.x(), a.point.y()},
                       (preferred ? -1 : -1001) - static_cast<int>(p));
        }
      }
    }
    for (frLayerNum fr = std::max(0, bottom - 2); fr <= sig_top; fr++) {
      frRegionQuery::Objects<frBlockObject> result;
      rq->query(window, fr, result);
      for (const auto& [rect, obj] : result) {
        if ((ablation & kAblateForeignShapes) && !owned_by(obj, ia.inst)) {
          continue;
        }
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
    for (int z = 0; z < std::min(num_z, 2) && !(ablation & kAblateTracks);
         z++) {
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
    return sig.finish();
  };
  const int threads = static_cast<int>(workers_.size());
  // D10: catalogues per pin class (key of the pin's own search inputs)
  // instead of per circuit class.
  const bool pin_classes
      = params_.offtrack_pin_classes && params_.access_classes;
  // M3 (offtrack_access_stats): circuit classes with one signature
  // component left out, and the pin-key classes, without changing results.
  const bool shadow_stats = params_.offtrack_access_stats;
  std::array<int64_t, 4> ablated_classes{};
  if (!pin_classes || shadow_stats) {
#pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
    for (int k = 0; k < static_cast<int>(instances.size()); k++) {
      instances[k].signature = instance_signature(instances[k], k, 0);
    }
  }
  if (shadow_stats) {
    for (int a = 0; a < 4; a++) {
      std::vector<std::vector<int64_t>> ablated(instances.size());
#pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
      for (int k = 0; k < static_cast<int>(instances.size()); k++) {
        ablated[k] = instance_signature(instances[k], k, 1 << a);
      }
      br::CircuitClassTable table;
      for (const std::vector<int64_t>& sig : ablated) {
        table.classify(sig);
      }
      ablated_classes[a] = table.numClasses();
    }
  }
  phases.signatures = phase_timer.elapsed();
  phase_timer.reset();
  br::CircuitClassTable table;
  std::vector<int> representative;
  for (int k = 0; k < static_cast<int>(instances.size()) && !pin_classes;
       k++) {
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
  if (shadow_stats && pin_classes) {
    br::CircuitClassTable full;
    for (const InstanceAccess& ia : instances) {
      full.classify(ia.signature);
    }
    stats.classes = full.numClasses();
  }
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
  // offtrack_algorithm3_fallback: an instance whose selection leaves a pin
  // with a catalogue unassigned tries Algorithm 3 catalogues (equal costs,
  // other equal-cost paths) and keeps them if they assign more pins.
  auto algorithm3_fallback = [&](const InstanceAccess& ia,
                                 std::vector<std::vector<AccessOption>>& cats,
                                 std::vector<int>& choice,
                                 const std::vector<int>& nets) {
    if (!params_.offtrack_algorithm3_fallback
        || !params_.offtrack_hanan_grid) {
      return;
    }
    const int n = static_cast<int>(cats.size());
    int assigned = 0;
    bool missing = false;
    for (int p = 0; p < n; p++) {
      assigned += choice[p] >= 0;
      missing = missing || (choice[p] < 0 && !cats[p].empty());
    }
    if (!missing) {
      return;
    }
    stats.fallback_instances++;
    std::vector<std::vector<AccessOption>> alt(n);
    for (int p = 0; p < n; p++) {
      const auto [net, t] = ia.terms[p];
      offTrackCatalogue(
          net, nets_[net].terms[t], nullptr, alt[p], stats, 0, true);
    }
    std::vector<int> alt_choice = solve(alt, nets);
    const int alt_assigned = static_cast<int>(
        std::ranges::count_if(alt_choice, [](const int c) { return c >= 0; }));
    if (alt_assigned > assigned) {
      cats = std::move(alt);
      choice = std::move(alt_choice);
      stats.fallback_improved++;
    }
  };
  // D10 pin keys (also for the M3 statistics): per access pin, the inputs
  // of its own catalogue search that are equal between instances of a
  // master, in the instance's local frame. `environment` adds every other
  // shape in the pin's reach (the sampled comparison's test whether a
  // member sees exactly the representative's surroundings).
  std::vector<std::vector<int>> pin_class(instances.size());
  std::vector<std::pair<int, int>> pin_representative;
  std::vector<int> master_of(instances.size());
  // Search reach of a pin: targets within the radius, the window margin
  // (2 tau + pitch) and two pitches for obstacles entering the window.
  const frCoord pin_reach
      = (params_.offtrack_radius + 3) * max_pitch + 2 * max_tau;
  auto pin_key = [&](const InstanceAccess& ia,
                     const size_t p,
                     const int master,
                     const bool environment) {
    const auto [net, t] = ia.terms[p];
    const Term& term = nets_[net].terms[t];
    br::CircuitSignature key(ia.place.inverse());
    key.addValue(master);
    key.addValue(ia.pin_index[p]);
    key.addValue(environment ? 1 : 0);
    key.addValue(params_.offtrack_radius);
    key.addValue(params_.max_access_options);
    key.addValue(params_.jog_factor);
    key.addValue(params_.pin_cost);
    key.addPin(ia.pin_index[p], 1);
    odb::Rect box = ia.inst->getBBox();
    if (params_.offtrack_from_access_points) {
      std::vector<AccessSource> aps;
      accessSources(term, aps);
      bool first = true;
      for (const AccessSource& a : aps) {
        const bool marked
            = params_.offtrack_preferred_access_points && !a.preferred;
        key.addShape(a.layer,
                     {a.point.x(), a.point.y(), a.point.x(), a.point.y()},
                     marked ? -2 : -1);
        if (a.via != nullptr) {
          const frVia v(a.via, a.point);
          key.addShape(a.via->getLayer1Num(), toBox(v.getLayer1BBox()), -3);
          key.addShape(a.via->getCutLayerNum(), toBox(v.getCutBBox()), -3);
          key.addShape(a.via->getLayer2Num(), toBox(v.getLayer2BBox()), -3);
        }
        const odb::Rect point(a.point, a.point);
        if (first) {
          box = point;
          first = false;
        } else {
          box.merge(point);
        }
      }
    }
    odb::Rect window = box;
    window.bloat(pin_reach, window);
    for (frLayerNum fr = std::max(0, bottom - 2); fr <= sig_top; fr++) {
      frRegionQuery::Objects<frBlockObject> result;
      rq->query(window, fr, result);
      for (const auto& [rect, obj] : result) {
        const bool own = owned_by(obj, ia.inst);
        if (!own && !environment) {
          continue;
        }
        const frNet* owner_net = nullptr;
        if (obj->typeId() == frcInstTerm) {
          owner_net = static_cast<const frInstTerm*>(obj)->getNet();
        } else if (obj->typeId() == frcBTerm) {
          owner_net = static_cast<const frBTerm*>(obj)->getNet();
        }
        const bool same_net
            = owner_net != nullptr && owner_net == nets_[net].net;
        key.addShape(fr,
                     toBox(rect.intersect(window)),
                     (own ? 1 : 3) + (same_net ? 1 : 0));
      }
      const frViaDef* via = viaDefAbove(fr);
      if (via != nullptr) {
        const frVia v0(via, odb::Point(0, 0));
        key.addVia(fr,
                   {{via->getLayer1Num(), toBox(v0.getLayer1BBox())},
                    {via->getCutLayerNum(), toBox(v0.getCutBBox())},
                    {via->getLayer2Num(), toBox(v0.getLayer2BBox())}});
      }
    }
    for (int z = 0; z < std::min(num_z, 2); z++) {
      const br::TrackLayer& l = graph_.layers[z];
      const frCoord lo = l.horizontal ? window.yMin() : window.xMin();
      const frCoord hi = l.horizontal ? window.yMax() : window.xMax();
      for (auto it = std::ranges::lower_bound(l.tracks, lo);
           it != l.tracks.end() && *it <= hi;
           ++it) {
        key.addTrack(z, l.horizontal, *it);
      }
    }
    return key.finish();
  };
  if (pin_classes || shadow_stats) {
    phase_timer.reset();
    std::unordered_map<const frMaster*, int> master_ids;
    for (size_t k = 0; k < instances.size(); k++) {
      const auto [it, inserted] = master_ids.emplace(
          instances[k].inst->getMaster(), static_cast<int>(master_ids.size()));
      master_of[k] = it->second;
    }
    std::vector<std::vector<std::vector<int64_t>>> keys(instances.size());
#pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
    for (int k = 0; k < static_cast<int>(instances.size()); k++) {
      const InstanceAccess& ia = instances[k];
      for (size_t p = 0; p < ia.terms.size(); p++) {
        keys[k].push_back(pin_key(ia, p, master_of[k], false));
      }
    }
    phases.signatures += phase_timer.elapsed();
    phase_timer.reset();
    br::CircuitClassTable pin_table;
    for (size_t k = 0; k < instances.size(); k++) {
      for (size_t p = 0; p < keys[k].size(); p++) {
        const auto [cls, is_new] = pin_table.classify(keys[k][p]);
        pin_class[k].push_back(cls);
        if (is_new) {
          pin_representative.emplace_back(static_cast<int>(k),
                                          static_cast<int>(p));
        }
      }
      keys[k].clear();
      keys[k].shrink_to_fit();
    }
    stats.pin_classes = pin_table.numClasses();
    phases.classification += phase_timer.elapsed();
    phase_timer.reset();
  }
  if (pin_classes) {
    // Catalogue of each pin class at its first pin, in that pin's own
    // surroundings (fixed shapes and pins), in parallel.
    phase_timer.reset();
    std::vector<std::vector<AccessOption>> pin_catalogue(
        pin_representative.size());
    std::vector<AccessStats> pin_stats(pin_representative.size());
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (int c = 0; c < static_cast<int>(pin_representative.size()); c++) {
      const auto [k, p] = pin_representative[c];
      const auto [net, t] = instances[k].terms[p];
      offTrackCatalogue(
          net, nets_[net].terms[t], nullptr, pin_catalogue[c], pin_stats[c]);
    }
    for (const AccessStats& s : pin_stats) {
      stats.add(s);
      class_total.add(s);
    }
    phases.class_catalogues = phase_timer.elapsed();
    phase_timer.reset();
    // Per instance, in parallel and before any reservation: each pin maps
    // its class catalogue and keeps the paths legal here (or generates its
    // own catalogue), then the instance's conflict-free selection.
    struct MemberData
    {
      std::vector<std::vector<AccessOption>> catalogues;
      std::vector<int> choice;
    };
    std::vector<MemberData> members(instances.size());
    std::vector<AccessStats> member_stats(instances.size());
    const int sample = params_.offtrack_catalogue_sample;
#pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
    for (int k = 0; k < static_cast<int>(instances.size()); k++) {
      const InstanceAccess& ia = instances[k];
      MemberData& md = members[k];
      AccessStats& ms = member_stats[k];
      const int n = static_cast<int>(ia.terms.size());
      md.catalogues.resize(n);
      std::vector<int> nets;
      for (const auto& [net, t] : ia.terms) {
        nets.push_back(net);
      }
      for (int p = 0; p < n; p++) {
        const auto [net, t] = ia.terms[p];
        const Term& term = nets_[net].terms[t];
        const int c = pin_class[k][p];
        const auto [rk, rp] = pin_representative[c];
        const std::vector<AccessOption>& source = pin_catalogue[c];
        std::vector<AccessOption>& cat = md.catalogues[p];
        if (rk == k && rp == p) {
          cat = source;
          continue;
        }
        const br::Isometry m
            = ia.place.compose(instances[rk].place.inverse());
        std::vector<AccessSource> member_aps;
        accessSources(term, member_aps);
        bool mapped = true;
        for (const AccessOption& opt : source) {
          AccessOption mo;
          if (!map_option(opt, m, member_aps, mo)) {
            mapped = false;
            break;
          }
          ms.mapped_options++;
          if (optionLegal(mo, net, nullptr)) {
            cat.push_back(std::move(mo));
          } else {
            ms.mapped_dropped++;
          }
        }
        const size_t needed = std::min<size_t>(
            static_cast<size_t>(params_.offtrack_regenerate), source.size());
        if (!mapped) {
          ms.mapping_failed++;
          offTrackCatalogue(net, term, nullptr, cat, ms);
          continue;
        }
        if (cat.size() < needed) {
          ms.member_regenerated++;
          offTrackCatalogue(net, term, nullptr, cat, ms);
          continue;
        }
        if (sample > 0 && (static_cast<int64_t>(k) * 31 + p) % sample == 0) {
          // Sampled comparison with the catalogue generated here.
          AccessStats own_stats;
          std::vector<AccessOption> own;
          offTrackCatalogue(net, term, nullptr, own, own_stats);
          auto pairs = [](const std::vector<AccessOption>& options) {
            std::vector<std::pair<size_t, Cost>> out;
            for (const AccessOption& o : options) {
              out.emplace_back(o.vertex, o.cost);
            }
            std::ranges::sort(out);
            return out;
          };
          const auto a = pairs(cat);
          const auto b = pairs(own);
          std::vector<std::pair<size_t, Cost>> only;
          std::ranges::set_difference(a, b, std::back_inserter(only));
          const bool mapped_only = !only.empty();
          only.clear();
          std::ranges::set_difference(b, a, std::back_inserter(only));
          const bool own_only = !only.empty();
          ms.sample_pins++;
          ms.sample_equal += a == b;
          ms.sample_mapped_only += mapped_only;
          ms.sample_own_only += own_only;
          const InstanceAccess& rep = instances[rk];
          if (pin_key(ia, p, master_of[k], true)
              == pin_key(rep, rp, master_of[rk], true)) {
            ms.sample_same_environment++;
            const br::Isometry linear = m.linear();
            if (linear.m00 == 1 && linear.m11 == 1 && linear.m01 == 0
                && linear.m10 == 0) {
              ms.sample_same_environment_translated++;
              ms.sample_same_environment_translated_unequal += a != b;
            }
          }
        }
      }
      std::vector<const std::vector<AccessOption>*> ptrs;
      for (const std::vector<AccessOption>& cat : md.catalogues) {
        ptrs.push_back(&cat);
      }
      md.choice = choosePrimaries(ptrs, nets, ms);
    }
    for (const AccessStats& s : member_stats) {
      stats.add(s);
    }
    phases.member_mapping = phase_timer.elapsed();
    phase_timer.reset();
    // In instance order: keep the selection where its primary paths are
    // legal against the reservations of earlier instances, else re-solve;
    // reserve.
    for (int k = 0; k < static_cast<int>(instances.size()); k++) {
      const InstanceAccess& ia = instances[k];
      MemberData& md = members[k];
      const int n = static_cast<int>(ia.terms.size());
      std::vector<int> nets;
      for (const auto& [net, t] : ia.terms) {
        nets.push_back(net);
      }
      std::vector<std::vector<AccessOption>>& cats = md.catalogues;
      std::vector<int> choice = md.choice;
      step_timer.reset();
      bool reuse = true;
      for (int p = 0; p < n && reuse; p++) {
        reuse = choice[p] < 0
                    ? cats[p].empty()
                    : optionLegal(cats[p][choice[p]], nets[p], nullptr);
      }
      phases.member_legality += step_timer.elapsed();
      if (reuse) {
        stats.class_reused++;
      } else {
        step_timer.reset();
        stats.reselected++;
        choice = solve(cats, nets);
        phases.member_solve += step_timer.elapsed();
      }
      step_timer.reset();
      algorithm3_fallback(ia, cats, choice, nets);
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
      md = MemberData();
    }
  }
  for (int k = 0; k < static_cast<int>(instances.size()) && !pin_classes;
       k++) {
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
    algorithm3_fallback(ia, cats, choice, nets);
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
    std::vector<int> choice = solve(cats, {net_id});
    if (choice[0] < 0 && params_.offtrack_algorithm3_fallback
        && params_.offtrack_hanan_grid) {
      // A top-level pin without a legal path: Algorithm 3 catalogue.
      stats.fallback_instances++;
      std::vector<std::vector<AccessOption>> alt(1);
      offTrackCatalogue(net_id, term, nullptr, alt[0], stats, 0, true);
      const std::vector<int> alt_choice = solve(alt, {net_id});
      if (alt_choice[0] >= 0) {
        cats = std::move(alt);
        choice = alt_choice;
        stats.fallback_improved++;
      }
    }
    keepStubs(term, cats[0]);
    if (choice[0] >= 0) {
      reserve(net, term, choice[0]);
    }
  }
  phases.top_level_pins = phase_timer.elapsed();
  if (params_.verbose && params_.offtrack_algorithm3_fallback) {
    logger_->info(utl::DRT,
                  9165,
                  "BonnRoute off-track Algorithm 3 fallback: {} instances with "
                  "a pin without a primary path on the Hanan grid, {} assign "
                  "more pins with Algorithm 3 catalogues.",
                  stats.fallback_instances,
                  stats.fallback_improved);
  }
  for (NetData& net : nets_) {
    if (pre_skip[net.id]) {
      continue;
    }
    net.skip = std::ranges::any_of(
        net.terms, [](const Term& term) { return term.options.empty(); });
  }
  stats.seconds = timer.elapsed();
  if (params_.verbose && params_.offtrack_access_report > 0) {
    reportOffTrackAccess();
  }
  if (params_.verbose) {
    logger_->report(
        "BonnRoute off-track access: {} instances in {} circuit classes, {} "
        "reused class solutions, {} re-solved, {} mapping failures; {} pin "
        "catalogues, {} paths ({:.2f} per pin), {} empty, {} fallback stubs; "
        "{} DRC-rejected, {} same-net-rejected, {} off-grid-rejected; {} "
        "blockage-grid vertices, "
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
        stats.offgrid_rejected,
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
    if (pin_classes) {
      logger_->report(
          "BonnRoute off-track pin classes: {} pin classes; {} mapped paths, "
          "{} illegal at the member, {} pins with a failed mapping, {} pins "
          "regenerated; {} instances kept their selection, {} re-solved; "
          "sample: {} pins, {} equal, {} with mapped-only (target, cost), {} "
          "with own-only, {} with the representative's environment ({} "
          "translated, {} of those unequal).",
          stats.pin_classes,
          stats.mapped_options,
          stats.mapped_dropped,
          stats.mapping_failed,
          stats.member_regenerated,
          stats.class_reused,
          stats.reselected,
          stats.sample_pins,
          stats.sample_equal,
          stats.sample_mapped_only,
          stats.sample_own_only,
          stats.sample_same_environment,
          stats.sample_same_environment_translated,
          stats.sample_same_environment_translated_unequal);
    }
    if (shadow_stats) {
      logger_->report(
          "BonnRoute off-track access stats: {} instances, {} circuit "
          "classes; without net groups {}, without preferred marks {}, "
          "without foreign shapes {}, without tracks {}; {} pin-key classes "
          "of {} access pins; blockage grid lines {} (other grid {}), "
          "vertices {} (other grid {}), search states {}, settled {}.",
          stats.instances,
          stats.classes,
          ablated_classes[0],
          ablated_classes[1],
          ablated_classes[2],
          ablated_classes[3],
          stats.pin_classes,
          [&] {
            int64_t pins = 0;
            for (const std::vector<int>& c : pin_class) {
              pins += static_cast<int64_t>(c.size());
            }
            return pins;
          }(),
          stats.grid_lines,
          stats.grid_lines_other,
          stats.grid_vertices,
          stats.grid_vertices_other,
          stats.search_states,
          stats.settled);
    }
  }
}

std::string BonnRouteDR::termName(const Term& term) const
{
  if (term.inst_term != nullptr) {
    return term.inst_term->getName() + " ("
           + term.inst_term->getInst()->getMaster()->getName() + ")";
  }
  if (term.bterm != nullptr) {
    return term.bterm->getName() + " (top)";
  }
  return "?";
}

void BonnRouteDR::reportOffTrackAccess() const
{
  int64_t pins = 0;
  int64_t options = 0;
  int64_t pairs = 0;
  int64_t sharing_pairs = 0;
  int64_t common = 0;
  int64_t wire = 0;
  int64_t distinct_wire = 0;
  int64_t catalogue_preferred = 0;
  int64_t catalogue_off_track = 0;
  int64_t primaries = 0;
  int64_t primary_preferred = 0;
  int64_t primary_off_track = 0;
  int64_t primary_tracks = 0;
  int64_t primary_missing = 0;
  int64_t primary_common = 0;
  std::vector<std::string> without;
  for (const NetData& net : nets_) {
    for (const Term& term : net.terms) {
      std::vector<std::vector<br::CatalogueWire>> paths;
      std::vector<int64_t> off_track_of;
      for (const AccessOption& opt : term.options) {
        if (!opt.offtrack) {
          continue;
        }
        std::vector<br::CatalogueWire>& path = paths.emplace_back();
        int64_t off = 0;
        for (const OutWire& w : opt.wires) {
          const bool horizontal = w.begin.y() == w.end.y();
          const br::GeoCoord line = horizontal ? w.begin.y() : w.begin.x();
          const br::GeoCoord a = horizontal ? w.begin.x() : w.begin.y();
          const br::GeoCoord b = horizontal ? w.end.x() : w.end.y();
          path.push_back(
              {w.layer, horizontal, line, std::min(a, b), std::max(a, b)});
          const int z = w.layer >= 0
                                && w.layer < static_cast<int>(z_of_fr_.size())
                            ? z_of_fr_[w.layer]
                            : -1;
          if (z < 0 || horizontal != layers_[z].horizontal) {
            continue;
          }
          const int64_t len = std::abs(int64_t{b} - a);
          catalogue_preferred += len;
          if (!std::ranges::binary_search(graph_.layers[z].tracks, line)) {
            catalogue_off_track += len;
            off += len;
          }
        }
        off_track_of.push_back(off);
      }
      if (paths.empty()) {
        continue;
      }
      pins++;
      options += static_cast<int64_t>(paths.size());
      const br::CatalogueSharing s = br::catalogueSharing(paths);
      pairs += s.pairs;
      sharing_pairs += s.sharing_pairs;
      common += s.common;
      wire += s.wire;
      distinct_wire += s.distinct_wire;
      if (term.primary >= 0 && term.options[term.primary].offtrack) {
        const AccessOption& opt = term.options[term.primary];
        primaries++;
        primary_common += s.common;
        for (const OutWire& w : opt.wires) {
          const int z = w.layer >= 0
                                && w.layer < static_cast<int>(z_of_fr_.size())
                            ? z_of_fr_[w.layer]
                            : -1;
          const bool horizontal = w.begin.y() == w.end.y();
          if (z >= 0 && horizontal == layers_[z].horizontal) {
            primary_preferred += std::abs(int64_t{w.end.x()} - w.begin.x())
                                 + std::abs(int64_t{w.end.y()} - w.begin.y());
          }
        }
        primary_off_track += off_track_of[term.primary];
        br::AccessCandidate c;
        candidateOf(opt, net.id, c);
        primary_tracks += static_cast<int64_t>(c.blocked_tracks.size());
        primary_missing += c.missing_directions;
      } else if (term.primary < 0) {
        int64_t legal = 0;
        for (const AccessOption& opt : term.options) {
          legal += opt.offtrack && optionLegal(opt, net.id, nullptr);
        }
        without.push_back(fmt::format("{}: net {}, {} catalogue paths, {} "
                                      "legal after setup, common wire {}",
                                      termName(term),
                                      net.net->getName(),
                                      paths.size(),
                                      legal,
                                      s.common ? 1 : 0));
      }
    }
  }
  std::ranges::sort(without);
  logger_->info(
      utl::DRT,
      9160,
      "BonnRoute off-track access report: {} pins with catalogues, {} paths; "
      "path pairs sharing wire {} of {} ({:.1f}%); pins whose paths all "
      "share a wire piece {} ({:.1f}%); distinct / total catalogue wire "
      "{:.3f}; preferred-direction catalogue wire off track {:.1f}%; {} "
      "primary paths: preferred-direction wire off track {:.1f}%, {:.2f} "
      "blocked tracks and {:.2f} missing directions per path, {} of them on "
      "pins whose paths all share a wire piece; {} pins with catalogues and "
      "no primary path.",
      pins,
      options,
      sharing_pairs,
      pairs,
      pairs > 0 ? 100.0 * sharing_pairs / pairs : 0.0,
      common,
      pins > 0 ? 100.0 * common / pins : 0.0,
      wire > 0 ? double(distinct_wire) / wire : 0.0,
      catalogue_preferred > 0 ? 100.0 * catalogue_off_track / catalogue_preferred
                              : 0.0,
      primaries,
      primary_preferred > 0 ? 100.0 * primary_off_track / primary_preferred
                            : 0.0,
      primaries > 0 ? double(primary_tracks) / primaries : 0.0,
      primaries > 0 ? double(primary_missing) / primaries : 0.0,
      primary_common,
      without.size());
  for (size_t k = 0;
       k < without.size()
       && k < static_cast<size_t>(params_.offtrack_access_report);
       k++) {
    logger_->info(utl::DRT,
                  9161,
                  "BonnRoute pin without a primary path: {}.",
                  without[k]);
  }
}

void BonnRouteDR::reportAccessFailures() const
{
  std::vector<std::pair<int, int>> failed;
  for (const Worker& w : workers_) {
    failed.insert(failed.end(),
                  w.access_failure_terms.begin(),
                  w.access_failure_terms.end());
  }
  std::ranges::sort(failed);
  std::vector<std::pair<std::pair<int, int>, int>> counted;
  for (const auto& f : failed) {
    if (!counted.empty() && counted.back().first == f) {
      counted.back().second++;
    } else {
      counted.push_back({f, 1});
    }
  }
  std::vector<std::string> pins;
  for (const auto& [f, count] : counted) {
    const NetData& net = nets_[f.first];
    const Term& term = net.terms[f.second];
    pins.push_back(fmt::format("{}: net {}, {} failed connects, primary {}, "
                               "routed {}",
                               termName(term),
                               net.net->getName(),
                               count,
                               term.primary >= 0 ? 1 : 0,
                               net.routed ? 1 : 0));
  }
  std::ranges::sort(pins);
  std::vector<std::string> left;
  int64_t left_without_primary = 0;
  int64_t left_with_failure = 0;
  for (const NetData& net : nets_) {
    if (net.routed || net.skip) {
      continue;
    }
    int without_primary = 0;
    int with_failure = 0;
    for (size_t t = 0; t < net.terms.size(); t++) {
      const Term& term = net.terms[t];
      without_primary += term.primary < 0
                         && std::ranges::any_of(term.options,
                                                [](const AccessOption& o) {
                                                  return o.offtrack;
                                                });
      with_failure += std::ranges::binary_search(
          failed, std::pair<int, int>(net.id, static_cast<int>(t)));
    }
    left_without_primary += without_primary > 0;
    left_with_failure += with_failure > 0;
    left.push_back(fmt::format("{}: {} pins, {} with catalogues and no "
                               "primary path, {} with failed connects",
                               net.net->getName(),
                               net.terms.size(),
                               without_primary,
                               with_failure));
  }
  std::ranges::sort(left);
  logger_->info(utl::DRT,
                9162,
                "BonnRoute access failures: {} pins in {} failed connects; "
                "{} nets left to FlexDR, {} with a pin without a primary "
                "path, {} with a pin of a failed connect.",
                pins.size(),
                failed.size(),
                left.size(),
                left_without_primary,
                left_with_failure);
  const size_t limit = static_cast<size_t>(params_.offtrack_access_report);
  for (size_t k = 0; k < pins.size() && k < limit; k++) {
    logger_->info(
        utl::DRT, 9163, "BonnRoute pin of failed connects: {}.", pins[k]);
  }
  for (size_t k = 0; k < left.size() && k < limit; k++) {
    logger_->info(utl::DRT, 9164, "BonnRoute net left to FlexDR: {}.", left[k]);
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
  const int around = dilation + params_.guide_dilation;
  // corridor_cost: the tiles around the guide and option tiles are the
  // widened part of the area.
  const bool priced = params_.corridor_cost > 0;
  auto allow = [&](int z, int tx, int ty) {
    for (int dz = -params_.layer_dilation; dz <= params_.layer_dilation; dz++) {
      for (int dx = -around; dx <= around; dx++) {
        for (int dy = -around; dy <= around; dy++) {
          if (priced && (dx != 0 || dy != 0)) {
            region.allowWidened(z + dz, tx + dx, ty + dy);
          } else {
            region.allow(z + dz, tx + dx, ty + dy);
          }
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
  // corridor_cost: a cell beside the path is plain when its planar tile is
  // a guide, source or target tile of the net on a layer within
  // layer_dilation of the cell's layer (the rule of buildRegion).
  auto plain = [&](const int zc, const int xc, const int yc) {
    for (int dz = -params_.layer_dilation; dz <= params_.layer_dilation; dz++) {
      const int64_t at = index(zc + dz, xc, yc);
      if (at >= 0 && worker.tile_guide[at] == epoch) {
        return true;
      }
    }
    return false;
  };
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
    const int around = params_.guide_dilation;
    for (int dz = -params_.layer_dilation; dz <= params_.layer_dilation; dz++) {
      for (int dx = -around; dx <= around; dx++) {
        for (int dy = -around; dy <= around; dy++) {
          // corridor_cost: tiles beside the path are the widened part,
          // unless they are tiles of the net themselves.
          if (params_.corridor_cost > 0 && (dx != 0 || dy != 0)
              && !plain(z + dz, x + dx, y + dy)) {
            region.allowWidened(z + dz, x + dx, y + dy);
          } else {
            region.allow(z + dz, x + dx, y + dy);
          }
        }
      }
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

namespace {

// Same-net via rules of one connection (Lef58Rules.md): the vias of the
// net's wire type against the net's cuts outside the search and against the
// via that entered the current layer (a bridge: down, along one layer, up
// again on the same via layer, or the reverse).
class NetViaRules final : public br::SameNetVias
{
 public:
  NetViaRules(const br::ShapeClassTable& classes,
              const br::TrackGraph& graph,
              const br::WireType& type)
      : classes_(classes), graph_(graph), type_(type)
  {
    const int num_vias = std::max(graph.numLayers() - 1, 0);
    cuts_.resize(num_vias);
    sorted_.assign(num_vias, 1);
    max_dx_.assign(num_vias, 0);
    layers.assign(num_vias, 0);
    cut_rule_.assign(num_vias, 0);
    for (int zv = 0; zv < num_vias && zv < static_cast<int>(type.via.size());
         zv++) {
      const br::ViaModel& m = type.via[zv];
      if (!m.valid) {
        continue;
      }
      // Rules in scope (Lef58Rules.md): the LEF58 cut spacing table of the
      // cut layer, the side rules of the pad layers.
      cut_rule_[zv] = !classes.layer(br::viaGridLayer(zv)).cut_classes.empty();
      const bool pads
          = classes.layer(br::wiringGridLayer(zv)).sideRules()
            || classes.layer(br::wiringGridLayer(zv + 1)).sideRules();
      layers[zv] = cut_rule_[zv] || pads;
      br::Coord extent = 0;
      for (const br::WireModel* w : {&m.bottom, &m.cut, &m.top}) {
        extent = std::max({extent, w->rect.dx(), w->rect.dy()});
      }
      reach
          = std::max(reach,
                     std::max({classes.ruleReach(br::viaGridLayer(zv)),
                               classes.ruleReach(br::wiringGridLayer(zv)),
                               classes.ruleReach(br::wiringGridLayer(zv + 1))})
                         + 2 * extent);
    }
  }

  bool any() const { return std::ranges::find(layers, 1) != layers.end(); }

  void addCut(const int zv, const br::Rect& rect, const br::ClassId cls)
  {
    if (zv >= 0 && zv < static_cast<int>(cuts_.size()) && cut_rule_[zv]) {
      cuts_[zv].emplace_back(rect, cls);
      sorted_[zv] = 0;
      max_dx_[zv] = std::max(max_dx_[zv], rect.dx());
    }
  }

  // Before a search: the cuts of each layer sorted by xlo.
  void sort()
  {
    for (size_t zv = 0; zv < cuts_.size(); zv++) {
      if (!sorted_[zv]) {
        std::ranges::sort(cuts_[zv], {}, [](const auto& c) {
          return std::pair(c.first.xlo, c.first.ylo);
        });
        sorted_[zv] = 1;
      }
    }
  }

  // planar_first: vias on the net's own pins that the searches may not place
  // (sorted lower vertices); their via layers become active.
  void setForbiddenPinVias(const std::vector<size_t>& vertices)
  {
    forbidden_ = &vertices;
    for (const size_t v : vertices) {
      int z, i, j;
      graph_.decode(v, z, i, j);
      if (z >= 0 && z < static_cast<int>(layers.size())) {
        layers[z] = 1;
      }
    }
  }

  bool allowed(const int zv, const int i, const int j) const override
  {
    if (forbidden_ != nullptr
        && std::ranges::binary_search(*forbidden_, graph_.vertex(zv, i, j))) {
      return false;
    }
    const br::ViaModel* m = model(zv);
    if (m == nullptr || !cut_rule_[zv] || cuts_[zv].empty()) {
      return true;
    }
    const br::Rect cut = m->cut.rect.translated(at(zv, i, j));
    const br::Coord d = classes_.ruleReach(br::viaGridLayer(zv));
    const std::vector<std::pair<br::Rect, br::ClassId>>& cuts = cuts_[zv];
    // Cuts whose xlo is within reach (cut extents bound the rest).
    auto it = std::ranges::lower_bound(
        cuts, cut.xlo - d - max_dx_[zv], {}, [](const auto& c) {
          return c.first.xlo;
        });
    for (; it != cuts.end() && it->first.xlo <= cut.xhi + d; ++it) {
      const auto& [rect, cls] = *it;
      if (rect == cut) {
        continue;  // the same via
      }
      if (br::shapesViolate(classes_,
                            br::PlacedShape{cut, m->cut.cls, 3},
                            br::PlacedShape{rect, cls, 3})) {
        return false;
      }
    }
    return true;
  }

  bool pairAllowed(const int zv,
                   const int i,
                   const int j,
                   const int i0,
                   const int j0,
                   const int turn_layer) const override
  {
    const br::ViaModel* m = model(zv);
    if (m == nullptr) {
      return true;
    }
    const br::Point a = at(zv, i, j);
    const br::Point b = at(zv, i0, j0);
    if (a == b) {
      return true;
    }
    if (cut_rule_[zv]
        && br::shapesViolate(
            classes_,
            br::PlacedShape{m->cut.rect.translated(a), m->cut.cls, 3},
            br::PlacedShape{m->cut.rect.translated(b), m->cut.cls, 3})) {
      return false;
    }
    // The pads on the layer the path does not turn on belong to different
    // polygons unless they touch.
    const int far_layer = turn_layer == zv ? zv + 1 : zv;
    if (!classes_.layer(br::wiringGridLayer(far_layer)).sideRules()) {
      return true;
    }
    const br::WireModel& far = turn_layer == zv ? m->top : m->bottom;
    const br::Rect pa = far.rect.translated(a);
    const br::Rect pb = far.rect.translated(b);
    if (pa.xlo <= pb.xhi && pb.xlo <= pa.xhi && pa.ylo <= pb.yhi
        && pb.ylo <= pa.yhi) {
      return true;
    }
    return !br::shapesViolate(classes_,
                              br::PlacedShape{pa, far.cls, 3},
                              br::PlacedShape{pb, far.cls, 3});
  }

 private:
  const br::ViaModel* model(const int zv) const
  {
    if (zv < 0 || zv >= static_cast<int>(type_.via.size())
        || !type_.via[zv].valid) {
      return nullptr;
    }
    return &type_.via[zv];
  }
  br::Point at(const int zv, const int i, const int j) const
  {
    const auto [x, y] = graph_.layers[zv].point(i, j);
    return {x, y};
  }

  const br::ShapeClassTable& classes_;
  const br::TrackGraph& graph_;
  const br::WireType& type_;
  std::vector<std::vector<std::pair<br::Rect, br::ClassId>>> cuts_;
  std::vector<uint8_t> sorted_;
  std::vector<uint8_t> cut_rule_;
  // Largest x extent of the cuts per layer: the window of allowed() starts
  // this much before the reach.
  std::vector<br::Coord> max_dx_;
  const std::vector<size_t>* forbidden_ = nullptr;
};

}  // namespace

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
  net.left_step = kLeftNone;
  net.left_term = -1;
  // keep_partial: no progress of an earlier net survives into this call.
  worker.progress.valid = false;
  worker.progress.net = -1;
  worker.footprint.setStep(kStepAccess);
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
  bool options_grown = false;
  auto option_valid = [&](const AccessOption& opt) {
    const br::RipupLevel level
        = worker.reader.level(opt.vertex, br::kWire, net.wire_type);
    bool ok = level == br::kLegal || (ripup && level <= worker.max_ripup_level);
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
    // D1 S2: a joined term keeps its chosen access path.
    if (net.partial && net.partial_state.connected[t]) {
      valid[t][net.partial_state.chosen[t]] = 1;
      continue;
    }
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
          = level == br::kLegal || (ripup && level <= worker.max_ripup_level);
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
      worker.footprint.setStep(kStepRegenerate);
      const bool regenerated = regenerateAccess(net, t, worker);
      worker.footprint.setStep(kStepAccess);
      if (regenerated) {
        options_grown = true;
        const Term& grown = net.terms[t];
        valid[t].resize(grown.options.size(), 0);
        for (size_t o = before; o < grown.options.size(); o++) {
          valid[t][o] = option_valid(grown.options[o]);
        }
      }
    }
    if (std::ranges::find(valid[t], 1) == valid[t].end()) {
      worker.access_failures++;
      net.left_step = kLeftAccess;
      net.left_term = t;
      if (params_.offtrack_access_report > 0) {
        worker.access_failure_terms.emplace_back(net.id, t);
      }
      return false;
    }
  }
  // Sec. 5.1: generated access paths extend the net's area. The thread
  // continues only where the area's tile box still lies inside its region.
  if (options_grown && worker.limit != nullptr) {
    const br::TileBox box = netAreaTileBox(net, extra_dilation);
    if (!worker.limit->contains(box)) {
      worker.area_deferred = true;
      worker.area_deferrals++;
      return false;
    }
    worker.footprint.widen(box);
  }
  const int retry = params_.retry_dilation + extra_dilation;
  worker.footprint.setStep(kStepSearch);
  if (ripup) {
    // A region's rip-up sequence searches only inside the region.
    if (worker.limit != nullptr) {
      const br::TileBox box = netAreaTileBox(net, retry);
      if (!worker.limit->contains(box)) {
        worker.area_deferred = true;
        return false;
      }
      worker.footprint.widen(box);
    }
    return connectInRegion(net, retry, true, valid, worker);
  }
  if (!net.dilated_area) {
    if (connectInRegion(net, extra_dilation, false, valid, worker)) {
      return true;
    }
    // A search skipped for the effort budget says nothing about the area.
    if (retry == extra_dilation || worker.effort_exceeded) {
      return false;
    }
    net.dilated_area = true;
  }
  // Sec. 5.1: the dilated area is searched where its tile box lies inside
  // the thread's region. Otherwise the net continues there in a later
  // partition, whose regions are larger, or in the tail.
  if (worker.limit != nullptr) {
    const br::TileBox box = netAreaTileBox(net, retry);
    if (!worker.limit->contains(box)) {
      worker.area_deferred = true;
      worker.area_deferrals++;
      return false;
    }
    worker.footprint.widen(box);
  }
  worker.footprint.setStep(kStepRetrySearch);
  return connectInRegion(net, retry, false, valid, worker);
}

bool BonnRouteDR::connectInRegion(NetData& net,
                                  const int dilation,
                                  const bool ripup,
                                  const std::vector<std::vector<char>>& valid,
                                  Worker& worker)
{
  // D1 S2: a partial net resumes from the connections it closed.
  const bool resume = net.partial;
  worker.progress.valid = false;
  if (resume) {
    net.wires = net.partial_state.wires;
    net.vias = net.partial_state.vias;
    worker.partial_resumes++;
  } else {
    net.wires.clear();
    net.vias.clear();
  }
  net.pin_points.clear();
  for (size_t t = 0; t < net.terms.size(); t++) {
    net.terms[t].chosen = resume ? net.partial_state.chosen[t] : -1;
  }
  // Effort budget of a speculative tail attempt (tail_defer_settled):
  // checked before each search, so an attempt settles at most the budget
  // plus one search; the net is then deferred (see routeSpeculative).
  auto over_budget = [&worker, &net]() {
    if (worker.effort_budget > 0
        && worker.settled - worker.effort_start >= worker.effort_budget) {
      worker.effort_exceeded = true;
    }
    if (worker.effort_exceeded) {
      net.left_step = kLeftSearchBudget;
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
  req.max_ripup_level = worker.max_ripup_level;
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
  req.corridor_cost_percent = params_.corridor_cost;
  req.use_intervals = params_.use_intervals;
  req.interval_via_property = params_.interval_via_property;
  req.interval_nearest_first = params_.interval_nearest_first;
  req.exact_future_cost = params_.exact_future_cost;
  // Same-net via rules (Lef58Rules.md): the net's cuts so far (earlier
  // connections and the access vias of connected terms) are synced before
  // every search.
  std::optional<NetViaRules> net_vias;
  size_t synced_vias = 0;
  std::vector<char> synced_terms(net.terms.size(), 0);
  if (params_.same_net_vias && !params_.use_intervals) {
    net_vias.emplace(shape_classes_, graph_, wire_types_[net.wire_type]);
    if (!net.forbidden_pin_vias.empty()) {
      net_vias->setForbiddenPinVias(net.forbidden_pin_vias);
    }
    if (net_vias->any()) {
      req.same_net_vias = &*net_vias;
    } else {
      net_vias.reset();
    }
  }
  auto add_cut = [&](const OutVia& v) {
    const frVia via(v.def, v.point);
    const frLayerNum layer = v.def->getCutLayerNum();
    const int g = gridLayer(layer);
    if (g < 0 || g >= track_grid_layers_ || g % 2 == 0) {
      return;
    }
    net_vias->addCut(
        (g - 1) / 2,
        brRect(via.getCutBBox()),
        shapeClass(RShape{layer, via.getCutBBox(), br::ShapeRole::kCut}, g));
  };
  auto sync_cuts = [&]() {
    if (!net_vias) {
      return;
    }
    for (; synced_vias < net.vias.size(); synced_vias++) {
      add_cut(net.vias[synced_vias]);
    }
    for (size_t t = 0; t < net.terms.size(); t++) {
      const Term& term = net.terms[t];
      if (synced_terms[t] || term.chosen < 0) {
        continue;
      }
      synced_terms[t] = 1;
      for (const OutVia& v : term.options[term.chosen].vias) {
        add_cut(v);
      }
    }
    net_vias->sort();
  };
  std::unordered_set<size_t> in_tree;
  std::vector<size_t> tree;
  std::vector<std::pair<size_t, Cost>> sources;
  std::vector<char> connected(num_terms, 0);
  connected[0] = 1;
  // Terms whose access path is part of a post-processed connection.
  std::vector<char> consumed(num_terms, 0);
  if (resume) {
    connected = net.partial_state.connected;
    consumed = net.partial_state.consumed;
    for (const size_t v : net.partial_state.tree) {
      if (in_tree.insert(v).second) {
        tree.push_back(v);
        sources.emplace_back(v, 0);
      }
    }
  } else {
    for (size_t o = 0; o < net.terms[0].options.size(); o++) {
      if (valid[0][o]) {
        sources.emplace_back(net.terms[0].options[o].vertex,
                             net.terms[0].options[o].cost);
      }
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
    sync_cuts();
    choose_potential(path_region);
    const double search_start = trace_connections_ ? omp_get_wtime() : 0.0;
    br::SearchResult res = params_.use_intervals
                               ? worker.interval_search->run(req)
                               : worker.search->run(req);
    if (trace_connections_) {
      // The area of the search: bounding box of its allowed tiles, their
      // number over all layers and the number of distinct planar tiles.
      std::vector<size_t> planar;
      planar.reserve(worker.region.touched.size());
      for (const auto& [z, index] : worker.region.touched) {
        planar.push_back(index);
      }
      std::ranges::sort(planar);
      const size_t planar_tiles
          = std::ranges::unique(planar).begin() - planar.begin();
      logger_->debug(
          utl::DRT,
          "bonnroute_connect",
          "connect net {} phase {} dilation {} path {} area {} {} {} {} tiles "
          "{} planar {} settled {} found {} ms {:.3f}",
          net.id,
          worker.speculative ? "tail"
                             : (worker.limit != nullptr ? "region" : "serial"),
          dilation,
          path_region ? 1 : 0,
          worker.region.bx_lo,
          worker.region.by_lo,
          worker.region.bx_hi,
          worker.region.by_hi,
          worker.region.touched.size(),
          planar_tiles,
          res.settled,
          res.found ? 1 : 0,
          1e3 * (omp_get_wtime() - search_start));
    }
    if (params_.potential_check && req.future_cost != nullptr) {
      // A feasible future cost leaves the shortest path cost unchanged.
      const br::FutureCostOracle* oracle = req.future_cost;
      req.future_cost = nullptr;
      const br::SearchResult plain = params_.use_intervals
                                         ? worker.interval_search->run(req)
                                         : worker.search->run(req);
      req.future_cost = oracle;
      worker.potential_checks++;
      worker.potential_check_pops += res.heap_pops;
      worker.potential_check_settled += res.settled;
      worker.potential_check_plain_pops += plain.heap_pops;
      worker.potential_check_plain_settled += plain.settled;
      if (plain.found != res.found || (res.found && plain.cost != res.cost)) {
        worker.potential_mismatches++;
      }
    }
    worker.searches++;
    worker.heap_pops += res.heap_pops;
    worker.held_pops += res.held_pops;
    worker.settled += res.settled;
    worker.fast_grid_queries += res.fast_grid_queries;
    worker.checker_queries += res.checker_queries;
    worker.aborted_searches += res.aborted;
    logSearchStatistics(net, worker, req, res, path_region, dilation, ripup);
    return res;
  };
  bool full_region = false;
  int remaining = static_cast<int>(
      std::ranges::count(connected, static_cast<char>(0)));
  // D1 S2: on failure the progress of this call (connections closed so
  // far) is kept for the caller.
  bool succeeded = false;
  struct OnExit
  {
    std::function<void()> f;
    ~OnExit() { f(); }
  };
  const OnExit keep_progress{[&]() {
    if (succeeded || !params_.keep_partial) {
      return;
    }
    PartialState& p = worker.progress;
    p.valid = true;
    p.net = net.id;
    p.joined = static_cast<int>(std::ranges::count(connected, char{1}));
    p.connected = connected;
    p.consumed = consumed;
    p.tree = tree;
    p.chosen.clear();
    for (const Term& term : net.terms) {
      p.chosen.push_back(term.chosen);
    }
    p.wires = net.wires;
    p.vias = net.vias;
  }};
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
        net.left_step = kLeftNoPath;
        return false;
      }
      res = run(false);
    }
    if (!res.found) {
      net.left_step = res.aborted ? kLeftSearchBudget : kLeftNoPath;
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
      net.left_step = kLeftNoPath;
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
        if (worker.limit != nullptr) {
          worker.history_bumps.push_back(v);
        } else {
          grid_.bumpHistory(v);
        }
      }
    }
  }
  succeeded = true;
  return true;
}

void BonnRouteDR::logSearchStatistics(const NetData& net,
                                      const Worker& worker,
                                      const br::SearchRequest& request,
                                      const br::SearchResult& result,
                                      const bool path_region,
                                      const int dilation,
                                      const bool ripup) const
{
  if (!logger_->debugCheck(utl::DRT, "bonnroute_search", 1)) {
    return;
  }
  // Mode: tail (speculative with an effort budget), tail-unbudgeted
  // (speculative without one: the deferred pass), serial-ripup (a search
  // that may rip up), region (any other: partitions and the rip-up phase's
  // plain connects); the stage follows from the next stage summary line.
  const char* mode = worker.speculative
                         ? (worker.effort_budget > 0 ? "tail"
                                                     : "tail-unbudgeted")
                         : (ripup ? "serial-ripup" : "region");
  const br::Region& region = *request.region;
  logger_->debug(
      utl::DRT,
      "bonnroute_search",
      "search net {} worker {} mode {} level {} dilation {} found {} "
      "aborted {} settled {} pops {} pushes {} fast_grid {} checker {} "
      "run_pops {} run_continuations {} sources {} targets {} target_boxes "
      "{} potential {} region_layer_tiles {} region_box {} {} {} {} terms "
      "{}",
      net.id,
      &worker - workers_.data(),
      mode,
      path_region ? "path" : "full",
      dilation,
      result.found ? 1 : 0,
      result.aborted ? 1 : 0,
      result.settled,
      result.heap_pops,
      result.heap_pushes,
      result.fast_grid_queries,
      result.checker_queries,
      result.run_label_pops,
      result.run_continuations,
      request.sources.size(),
      request.targets.size(),
      request.target_boxes.size(),
      request.future_cost != nullptr ? 1 : 0,
      region.touched.size(),
      region.bx_lo,
      region.by_lo,
      region.bx_hi,
      region.by_hi,
      net.terms.size());
}

void BonnRouteDR::uncommitShapes(NetData& net, Worker* local)
{
  // One exact fast-grid recomputation for all shapes of the net.
  br::FastGridBuilder::PendingRemovals pending;
  for (const RShape& s : net.shapes) {
    changeShape(s, routedType(s.role), net.id, false, local, &pending);
  }
  fast_builder_.flush(pending);
  net.shapes.clear();
}

void BonnRouteDR::commitPartial(NetData& net,
                                const PartialState& state,
                                Worker* local)
{
  if (state.net != net.id || state.connected.size() != net.terms.size()
      || state.chosen.size() != net.terms.size()) {
    throw std::logic_error("BonnRouteDR: partial routing of another net");
  }
  if (net.partial) {
    uncommitShapes(net, local);
  }
  net.partial_state = state;
  net.partial = true;
  net.routed = false;
  net.wires = state.wires;
  net.vias = state.vias;
  shapesOf(net, net.shapes);
  for (const RShape& s : net.shapes) {
    changeShape(s, routedType(s.role), net.id, true, local);
  }
  net.wires.clear();
  net.vias.clear();
}

void BonnRouteDR::commitNet(NetData& net, Worker* local)
{
  // D1 S2: the whole net replaces the shapes of its partial routing.
  if (net.partial) {
    uncommitShapes(net, local);
    net.partial = false;
    net.partial_state = PartialState{};
  }
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
            // Bridges collapse only away from the net's pins. A bridge via
            // whose point lies on a pin of the net is a pin contact (FlexDR
            // finds pins at via points and at truncated wire ends,
            // dr/FlexDR_conn.cpp); the jog would end inside that pin with
            // an extended end, which FlexDR does not see as connected
            // (ibex: a met3 port at the die edge, DRT-0206).
            bool on_pin = false;
            for (const odb::Point& p : {A, B}) {
              for (const frLayerNum layer : {lower, upper}) {
                forEachNetPinAt(design_->getRegionQuery(),
                                net.net,
                                layer,
                                p,
                                [&](const frBlockObject*, const odb::Rect&) {
                                  on_pin = true;
                                });
              }
            }
            if (on_pin) {
              stats.pin_contact++;
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

void BonnRouteDR::ripupNet(NetData& net, Worker* local)
{
  // One exact fast-grid recomputation for all shapes of the net.
  br::FastGridBuilder::PendingRemovals pending;
  for (const RShape& s : net.shapes) {
    changeShape(s, routedType(s.role), net.id, false, local, &pending);
  }
  // Incremental mode: input routing is ripped up with the net.
  if (net.input_present) {
    setInputShapes(net, false, &pending);
    net.input_removed = true;
  }
  fast_builder_.flush(pending);
  net.shapes.clear();
  net.wires.clear();
  net.vias.clear();
  net.routed = false;
  net.partial = false;
  net.partial_state = PartialState{};
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

br::TileBox BonnRouteDR::netCoreTileBox(const NetData& net) const
{
  frBlock* block = design_->getTopBlock();
  br::TileBox box;
  box.xhi = box.yhi = -1;
  auto add = [&box](int x, int y) {
    box.xlo = std::min(box.xlo, x);
    box.ylo = std::min(box.ylo, y);
    box.xhi = std::max(box.xhi, x);
    box.yhi = std::max(box.yhi, y);
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
  return box;
}

int BonnRouteDR::netTopLayer(const NetData& net) const
{
  const int num_z = graph_.numLayers();
  // Routing area (buildRegion, buildPathRegion): guide and option layers
  // plus layer_dilation.
  int top = -1;
  for (const auto& tile : net.guide_tiles) {
    top = std::max(top, tile[0]);
  }
  for (const Term& term : net.terms) {
    for (const AccessOption& option : term.options) {
      top = std::max(top, option.z);
    }
  }
  if (top < 0) {
    return num_z - 1;
  }
  return std::min(top + params_.layer_dilation, num_z - 1);
}

int BonnRouteDR::netAccessTopLayer(const NetData& net) const
{
  const int num_z = graph_.numLayers();
  int top = 0;
  for (const RShape& pin : net.pins) {
    if (pin.layer >= 0 && pin.layer < static_cast<int>(z_of_fr_.size())
        && z_of_fr_[pin.layer] >= 0) {
      top = std::max(top, z_of_fr_[pin.layer] + 1);
    }
  }
  return std::min(top, num_z - 1);
}

int BonnRouteDR::netHaloTiles(const NetData& net) const
{
  const int top = netTopLayer(net);
  const int access_top = netAccessTopLayer(net);
  return net.wide_shapes ? std::max(area_halo_tiles_[top],
                                    access_halo_tiles_[access_top])
                         : std::max(narrow_area_halo_tiles_[top],
                                    narrow_access_halo_tiles_[access_top]);
}

void BonnRouteDR::netTileBox(const NetData& net,
                             int& xlo,
                             int& ylo,
                             int& xhi,
                             int& yhi) const
{
  const br::TileBox core = netCoreTileBox(net);
  // Corridor width, retry dilation and the routing-space interaction halo
  // (Sec. 5.1).
  const int halo
      = params_.guide_dilation + params_.retry_dilation + netHaloTiles(net);
  xlo = core.xlo - halo;
  ylo = core.ylo - halo;
  xhi = core.xhi + halo;
  yhi = core.yhi + halo;
}

br::TileBox BonnRouteDR::netAreaTileBox(const NetData& net,
                                        const int dilation) const
{
  br::TileBox box = netCoreTileBox(net);
  const int margin = params_.guide_dilation + dilation + netHaloTiles(net);
  box.xlo = std::max(box.xlo - margin, 0);
  box.ylo = std::max(box.ylo - margin, 0);
  box.xhi = std::min(box.xhi + margin, tiles_x_ - 1);
  box.yhi = std::min(box.yhi + margin, tiles_y_ - 1);
  return box;
}

br::TileBox BonnRouteDR::netClippedTileBox(const NetData& net) const
{
  br::TileBox box;
  netTileBox(net, box.xlo, box.ylo, box.xhi, box.yhi);
  box.xlo = std::max(box.xlo, 0);
  box.ylo = std::max(box.ylo, 0);
  box.xhi = std::min(box.xhi, tiles_x_ - 1);
  box.yhi = std::min(box.yhi, tiles_y_ - 1);
  return box;
}

void BonnRouteDR::beginFootprint(Worker& worker,
                                 const NetData& net,
                                 const FootprintPhase phase,
                                 const br::TileBox& box)
{
  worker.footprint_phase = phase;
  worker.footprint.begin(net.id, box);
}

void BonnRouteDR::endFootprint(Worker& worker,
                               const NetData& net,
                               const br::TileBox* region)
{
  const br::AccessFootprint& footprint = worker.footprint;
  FootprintStats& stats = worker.footprint_stats;
  const int phase = worker.footprint_phase;
  stats.units[phase]++;
  if (footprint.unitViolations() > 0) {
    stats.outside_box_units[phase]++;
  }
  const br::TileBox accessed = footprint.accessed();
  if (region != nullptr && !region->contains(accessed)) {
    stats.outside_region_units[phase]++;
    if (stats.region_examples.size() < 8) {
      stats.region_examples.push_back(fmt::format(
          "net {} ({}): tiles [{}, {}] - [{}, {}], region [{}, {}] - [{}, {}]",
          net.id,
          net.net->getName(),
          accessed.xlo,
          accessed.ylo,
          accessed.xhi,
          accessed.yhi,
          region->xlo,
          region->ylo,
          region->xhi,
          region->yhi));
    }
  }
  // Reach beyond the core box, per step: reads and writes.
  const br::TileBox core = netCoreTileBox(net);
  for (int step = 0; step < kNumFootprintSteps; step++) {
    for (int write = 0; write < 2; write++) {
      const br::TileBox& box
          = write != 0 ? footprint.writes(step) : footprint.reads(step);
      if (!box.empty()) {
        const int reach
            = std::min(core.reachOf(box), FootprintStats::kReachBins - 1);
        stats.reach[phase][step][write][reach]++;
      }
    }
  }
}

void BonnRouteDR::reportFootprints()
{
  static constexpr std::array<const char*, kNumFootprintPhases> kPhaseNames{
      "region", "tail search", "tail commit"};
  static constexpr std::array<const char*, kNumFootprintSteps> kStepNames{
      "access", "regenerate", "search", "retry search", "commit", "release"};
  FootprintStats total;
  int64_t records = 0;
  int64_t violations = 0;
  std::vector<br::AccessFootprint::Violation> examples;
  for (const Worker& w : workers_) {
    const FootprintStats& s = w.footprint_stats;
    for (int phase = 0; phase < kNumFootprintPhases; phase++) {
      total.units[phase] += s.units[phase];
      total.outside_box_units[phase] += s.outside_box_units[phase];
      total.outside_region_units[phase] += s.outside_region_units[phase];
      for (int step = 0; step < kNumFootprintSteps; step++) {
        for (int write = 0; write < 2; write++) {
          for (int bin = 0; bin < FootprintStats::kReachBins; bin++) {
            total.reach[phase][step][write][bin]
                += s.reach[phase][step][write][bin];
          }
        }
      }
    }
    total.region_examples.insert(total.region_examples.end(),
                                 s.region_examples.begin(),
                                 s.region_examples.end());
    records += w.footprint.records();
    violations += w.footprint.violationCount();
    examples.insert(examples.end(),
                    w.footprint.violations().begin(),
                    w.footprint.violations().end());
  }
  int64_t outside_box = 0;
  int64_t outside_region = 0;
  for (int phase = 0; phase < kNumFootprintPhases; phase++) {
    outside_box += total.outside_box_units[phase];
    outside_region += total.outside_region_units[phase];
  }
  logger_->report(
      "BonnRoute DR access footprint check: {} region units, {} tail search "
      "units, {} tail commit units, {} recorded accesses; box margin {} "
      "tiles (retry dilation {} + halo {}) at most; {} accesses outside the "
      "net box in {} units, {} units outside their region.",
      total.units[kPhaseRegion],
      total.units[kPhaseTailSearch],
      total.units[kPhaseTailCommit],
      records,
      params_.guide_dilation + params_.retry_dilation + halo_tiles_,
      params_.guide_dilation + params_.retry_dilation,
      halo_tiles_,
      violations,
      outside_box,
      outside_region);
  for (int phase = 0; phase < kNumFootprintPhases; phase++) {
    for (int step = 0; step < kNumFootprintSteps; step++) {
      for (int write = 0; write < 2; write++) {
        const FootprintStats::Histogram& h = total.reach[phase][step][write];
        int64_t units = 0;
        std::string bins;
        for (int bin = 0; bin < FootprintStats::kReachBins; bin++) {
          units += h[bin];
          bins += fmt::format(" {}:{}", bin, h[bin]);
        }
        if (units == 0) {
          continue;
        }
        logger_->report(
            "BonnRoute DR access footprint reach: phase {} step {} {} units "
            "{} tiles beyond the core box{}.",
            kPhaseNames[phase],
            kStepNames[step],
            write != 0 ? "writes" : "reads",
            units,
            bins);
      }
    }
  }
  if (violations == 0 && outside_region == 0) {
    return;
  }
  std::string detail;
  for (size_t k = 0; k < examples.size() && k < 8; k++) {
    const br::AccessFootprint::Violation& v = examples[k];
    detail += fmt::format(
        "\n  net {} ({}) step {} {} layer {}: tiles [{}, {}] - [{}, {}], "
        "box [{}, {}] - [{}, {}]",
        v.unit,
        v.unit >= 0 && v.unit < static_cast<int>(nets_.size())
            ? nets_[v.unit].net->getName()
            : std::string("-"),
        kStepNames[std::clamp(v.step, 0, kNumFootprintSteps - 1)],
        br::accessKindName(v.kind),
        v.layer,
        v.tiles.xlo,
        v.tiles.ylo,
        v.tiles.xhi,
        v.tiles.yhi,
        v.allowed.xlo,
        v.allowed.ylo,
        v.allowed.xhi,
        v.allowed.yhi);
  }
  for (size_t k = 0; k < total.region_examples.size() && k < 8; k++) {
    detail += "\n  " + total.region_examples[k];
  }
  logger_->warn(utl::DRT,
                9130,
                "BonnRoute access footprint check failed: {} accesses outside "
                "the net box in {} units, {} units outside their region.{}",
                violations,
                outside_box,
                outside_region,
                detail);
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
  const bool scheduled = params_.tail_lookahead || params_.tail_pipeline;
  if (scheduled) {
    routeSpeculativePipeline(order, failed, deferred, stage, threads, budget);
  } else {
    routeSpeculativeBatches(order, failed, deferred, stage, threads, budget);
  }
  if (deferred.empty()) {
    return;
  }
  if (!params_.tail_deferred_pass) {
    // To the rip-up phase, in net order.
    failed.insert(failed.end(), deferred.begin(), deferred.end());
    return;
  }
  const std::string late = std::string(stage) + " deferred";
  const int64_t late_budget = params_.tail_deferred_settled;
  if (scheduled) {
    routeSpeculativePipeline(
        deferred, failed, none, late, threads, late_budget);
  } else {
    routeSpeculativeBatches(
        deferred, failed, none, late, threads, late_budget);
  }
  // Over the deferred pass's budget: to the rip-up phase.
  failed.insert(failed.end(), none.begin(), none.end());
}

br::TailBatchPolicy BonnRouteDR::tailBatchPolicy() const
{
  return {params_.tail_batch,
          params_.tail_batch_min,
          params_.tail_batch_high,
          params_.tail_batch_low};
}

void BonnRouteDR::routeSpeculativeBatches(const std::vector<int>& order,
                                          std::vector<int>& failed,
                                          std::vector<int>& deferred,
                                          const std::string& stage,
                                          const int threads,
                                          const int64_t budget)
{
  const size_t batch_size = static_cast<size_t>(params_.tail_batch);
  br::TailBatchSizes sizes(tailBatchPolicy());
  // Retries of the previous batch at the queue's front: all join the batch.
  size_t pending_retries = 0;
  int64_t previous_attempts = 0;
  int64_t previous_conflicts = 0;
  int64_t smallest_batch = 0;
  int64_t largest_batch = 0;
  std::deque<int> queue(order.begin(), order.end());
  std::vector<int> batch;
  std::vector<char> ok;
  std::vector<double> busy;
  std::vector<int64_t> settled;
  std::vector<char> over;
  // keep_partial: progress of the failed attempts of the batch.
  std::vector<PartialState> batch_progress;
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
    if (batches > 0) {
      sizes.record(
          static_cast<int>(batches) - 1, previous_attempts, previous_conflicts);
    }
    const size_t size = std::max<size_t>(
        pending_retries, static_cast<size_t>(sizes.size(batches)));
    while (!queue.empty() && batch.size() < size) {
      batch.push_back(queue.front());
      queue.pop_front();
    }
    smallest_batch = smallest_batch == 0
                         ? static_cast<int64_t>(batch.size())
                         : std::min<int64_t>(smallest_batch, batch.size());
    largest_batch = std::max<int64_t>(largest_batch, batch.size());
    ok.assign(batch.size(), 0);
    busy.assign(batch.size(), 0);
    settled.assign(batch.size(), 0);
    over.assign(batch.size(), 0);
    batch_progress.assign(batch.size(), PartialState{});
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
      const br::AccessFootprintScope footprint_scope(
          params_.partition_check ? &w.footprint : nullptr);
      if (params_.partition_check) {
        beginFootprint(w,
                       nets_[batch[k]],
                       kPhaseTailSearch,
                       netClippedTileBox(nets_[batch[k]]));
      }
      w.speculative = true;
      w.effort_budget = budget;
      w.effort_start = w.settled;
      w.effort_exceeded = false;
      ok[k] = connect(nets_[batch[k]], false, w);
      over[k] = !ok[k] && w.effort_exceeded;
      if (!ok[k] && params_.keep_partial
          && progressAdds(nets_[batch[k]], w.progress)) {
        batch_progress[k] = std::move(w.progress);
      }
      w.effort_budget = 0;
      w.effort_exceeded = false;
      w.speculative = false;
      settled[k] = w.settled - settled_before;
      if (params_.partition_check) {
        endFootprint(w, nets_[batch[k]], nullptr);
      }
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
      // partition_check: the decision and commit of the net, recorded on
      // the first worker (no search runs during the commits).
      Worker& commit_worker = workers_.front();
      // keep_partial: a failed attempt that closed connections commits
      // them like a route (checked against the batch's commits).
      const bool partial = !ok[k] && batch_progress[k].valid;
      const bool writes = ok[k] || partial;
      const br::AccessFootprintScope footprint_scope(
          params_.partition_check && writes ? &commit_worker.footprint
                                            : nullptr);
      if (params_.partition_check && writes) {
        beginFootprint(
            commit_worker, net, kPhaseTailCommit, netClippedTileBox(net));
        commit_worker.footprint.setStep(kStepCommit);
      }
      auto unrouted = [&]() {
        if (over[k]) {
          deferred.push_back(net.id);
          deferred_count++;
          outcome = 3;
        } else {
          failed.push_back(net.id);
        }
      };
      if (!writes) {
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
        if (clash && partial) {
          net.wires.clear();
          net.vias.clear();
          net.pin_points.clear();
          unrouted();
        } else if (clash) {
          net.wires.clear();
          net.vias.clear();
          net.pin_points.clear();
          retry.push_back(net.id);
          conflicts++;
          outcome = 2;
        } else if (partial) {
          commitPartial(net, batch_progress[k], nullptr);
          commit_worker.tail_partial_commits++;
          odb::Rect box;
          box.mergeInit();
          for (const RShape& shape : net.shapes) {
            box.merge(shape.rect);
          }
          if (!box.isInverted()) {
            box.bloat(static_cast<int>(interaction_reach_), box);
            committed_boxes.push_back(box);
          }
          unrouted();
          if (outcome == 0) {
            outcome = 4;
          }
        } else {
          commitNet(net, nullptr);
          commit_worker.footprint.setStep(kStepRelease);
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
        if (params_.partition_check) {
          endFootprint(commit_worker, net, nullptr);
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
    pending_retries = retry.size();
    previous_attempts = static_cast<int64_t>(batch.size());
    previous_conflicts = static_cast<int64_t>(retry.size());
    batches++;
  }
  if (params_.verbose && tailBatchPolicy().adaptive()) {
    logger_->info(utl::DRT,
                  9132,
                  "BonnRoute {} adaptive batches: {} batches of {} to {} nets "
                  "(tail_batch_min {}, tail_batch {}).",
                  stage,
                  batches,
                  smallest_batch,
                  largest_batch,
                  params_.tail_batch_min,
                  params_.tail_batch);
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
  // keep_partial (lookahead only): a failed attempt that closed connections
  // is reported as found and commits them at its decision.
  const bool tail_partial = params_.keep_partial && params_.tail_lookahead;
  std::vector<char> route_found(nets_.size(), 0);
  std::vector<char> partial_attempt(nets_.size(), 0);
  std::vector<PartialState> tail_progress(tail_partial ? nets_.size() : 0);
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
    const br::AccessFootprintScope footprint_scope(
        params_.partition_check ? &w.footprint : nullptr);
    if (params_.partition_check) {
      const auto& b = boxes[id];
      beginFootprint(
          w, nets_[id], kPhaseTailSearch, {b[0], b[1], b[2], b[3]});
    }
    w.speculative = true;
    w.effort_budget = budget;
    w.effort_start = w.settled;
    w.effort_exceeded = false;
    const bool found = connect(nets_[id], false, w);
    over[id] = !found && w.effort_exceeded;
    route_found[id] = found;
    bool progress = false;
    if (!found && tail_partial && progressAdds(nets_[id], w.progress)) {
      tail_progress[id] = std::move(w.progress);
      progress = true;
    }
    w.effort_budget = 0;
    w.effort_exceeded = false;
    w.speculative = false;
    settled[id] = w.settled - before;
    if (params_.partition_check) {
      endFootprint(w, nets_[id], nullptr);
    }
    search_seconds[id] = omp_get_wtime() - start;
    times[t].push_back(search_seconds[id]);
    return found || progress;
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
    const br::AccessFootprintScope footprint_scope(
        params_.partition_check && found ? &w.footprint : nullptr);
    if (params_.partition_check && found) {
      const auto& b = boxes[id];
      beginFootprint(w, net, kPhaseTailCommit, {b[0], b[1], b[2], b[3]});
      w.footprint.setStep(kStepCommit);
    }
    // keep_partial: found without a route = connections to commit.
    const bool partial = found && !route_found[id];
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
      if (clash && partial) {
        net.wires.clear();
        net.vias.clear();
        net.pin_points.clear();
        decision = Decision::kFailed;
      } else if (clash) {
        net.wires.clear();
        net.vias.clear();
        net.pin_points.clear();
        decision = Decision::kConflict;
      } else if (partial) {
        // Committed for the batch (later decisions see it); the net stays
        // unrouted (see the failed list below).
        commitPartial(net, tail_progress[id], &w);
        partial_attempt[id] = 1;
        w.tail_partial_commits++;
        odb::Rect box;
        box.mergeInit();
        for (const RShape& shape : net.shapes) {
          box.merge(shape.rect);
        }
        if (!box.isInverted()) {
          box.bloat(static_cast<int>(interaction_reach_), box);
          committed_box[id] = box;
        }
      } else {
        commitNet(net, &w);
        // Released and recomputed exactly inside the net's box.
        w.footprint.setStep(kStepRelease);
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
      if (params_.partition_check) {
        endFootprint(w, net, nullptr);
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
  std::vector<int> unrouted;
  br::TailPipeline::Stats stats;
  if (params_.tail_lookahead) {
    br::TailLookahead lookahead(tailBatchPolicy());
    std::vector<br::TailLookahead::Attempt> attempts;
    stats = lookahead.run(order,
                          boxes,
                          threads,
                          hooks,
                          unrouted,
                          tail_partial ? &attempts : nullptr);
    if (tail_partial) {
      // Unrouted: failed attempts and partial commits, in (batch, position)
      // order, as the barrier batches list them.
      unrouted.clear();
      for (const br::TailLookahead::Attempt& a : attempts) {
        if (a.decision == br::TailPipeline::Decision::kFailed
            || (a.decision == br::TailPipeline::Decision::kCommitted
                && partial_attempt[a.net])) {
          unrouted.push_back(a.net);
        }
      }
    }
    if (params_.verbose && tailBatchPolicy().adaptive()) {
      logger_->info(utl::DRT,
                    9132,
                    "BonnRoute {} adaptive batches: {} batches of {} to {} "
                    "nets (tail_batch_min {}, tail_batch {}).",
                    stage,
                    stats.batches,
                    stats.smallest_batch,
                    stats.largest_batch,
                    params_.tail_batch_min,
                    params_.tail_batch);
    }
  } else {
    br::TailPipeline pipeline(params_.tail_batch, params_.tail_retry_delay);
    stats = pipeline.run(order, boxes, threads, hooks, unrouted);
  }
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
    if (params_.tail_lookahead) {
      logger_->info(
          utl::DRT,
          9140,
          "BonnRoute {} speculative lookahead: {} batches of up to {} nets, "
          "{} nets committed, {} searched again after a conflict, {} "
          "without a route; {:.2f}s on {} threads, thread utilization "
          "{:.0f}% (searches {:.2f}s, decisions {:.2f}s incl. {} exact "
          "checks {:.2f}s, idle {:.2f}s); {} searches started before their "
          "batch; search time p50 {:.1f} ms, p90 {:.1f} ms, p99 {:.1f} ms, "
          "max {:.1f} ms; {} nets deferred over an effort of {} settled "
          "vertices; routing digest {:016x}.",
          stage,
          stats.batches,
          params_.tail_batch,
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
      return;
    }
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
  // Sec. 5.1: a net enters a region that contains the tile box of its
  // guide area (netAreaTileBox). A net whose connections failed there
  // enters with the box of the dilated area; connect() tests that box
  // against the region before it searches the dilated area.
  std::vector<std::array<int, 4>> boxes(nets_.size());
  auto set_box = [&](const int id) {
    const NetData& net = nets_[id];
    const br::TileBox box
        = netAreaTileBox(net, net.dilated_area ? params_.retry_dilation : 0);
    boxes[id] = {box.xlo, box.ylo, box.xhi, box.yhi};
  };
  for (const int id : order) {
    set_box(id);
  }
  // Marks of the nets that stay pending after a partition.
  std::vector<char> stays(nets_.size(), 0);
  const int threads = static_cast<int>(workers_.size());
  int count = initial_regions > 0             ? initial_regions
              : params_.partition_regions > 0 ? params_.partition_regions
                                              : 2 * threads;
  // The tail takes over below this many regions (a partition with fewer
  // regions than threads leaves threads idle).
  const int tail_below
      = params_.tail_regions > 0 ? params_.tail_regions : threads;
  std::vector<int> pending = order;
  std::vector<std::array<int, 4>> regions;
  int partition = 0;
  while (!pending.empty()) {
    if (params_.tail_batch > 0 && count < tail_below) {
      routeSpeculative(pending, failed, stage);
      return;
    }
    if (count <= 1) {
      regions.assign(1, {0, 0, tiles_x_ - 1, tiles_y_ - 1});
    } else {
      partitionRegions(pending, boxes, count, partition % 2 == 1, regions);
    }
    if (trace_connections_) {
      for (const auto& g : regions) {
        logger_->debug(utl::DRT,
                       "bonnroute_connect",
                       "partition {} {} region {} {} {} {}",
                       stage,
                       partition,
                       g[0],
                       g[1],
                       g[2],
                       g[3]);
      }
    }
    std::vector<std::vector<int>> members(regions.size());
    std::vector<std::vector<int>> region_failed(regions.size());
    // Nets that need the dilated area, whose box leaves the region.
    std::vector<std::vector<int>> region_again(regions.size());
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
      const int64_t region_settled_start = w.settled;
      w.region_ripup_settled = 0;
      const br::AccessFootprintScope footprint_scope(
          params_.partition_check ? &w.footprint : nullptr);
      const br::TileBox region{
          regions[r][0], regions[r][1], regions[r][2], regions[r][3]};
      for (const int id : members[r]) {
        NetData& net = nets_[id];
        if (params_.partition_check) {
          const auto& b = boxes[id];
          beginFootprint(w, net, kPhaseRegion, {b[0], b[1], b[2], b[3]});
        }
        w.limit = &region;
        w.area_deferred = false;
        const bool connected = connect(net, false, w);
        // D1 S2: a failed net keeps the connections it closed, when it
        // joined more terms than before.
        if (!connected && params_.keep_partial
            && progressAdds(net, w.progress)) {
          w.footprint.setStep(kStepCommit);
          commitPartial(net, w.progress, &w);
          w.partial_commits++;
        }
        // D3: a failure inside the region gets a rip-up sequence there,
        // which commits the net when it succeeds.
        const bool by_ripup = !connected && !w.area_deferred
                              && params_.region_ripup_effort > 0
                              && regionRipup(net, w, region, region_settled_start);
        w.limit = nullptr;
        if (connected) {
          w.footprint.setStep(kStepCommit);
          commitNet(net, &w);
          w.footprint.setStep(kStepRelease);
          setReservations(net, false, &w);
        } else if (!by_ripup) {
          net.wires.clear();
          net.vias.clear();
          (w.area_deferred ? region_again : region_failed)[r].push_back(id);
        }
        if (params_.partition_check) {
          endFootprint(w, net, &region);
        }
      }
      busy[r] = omp_get_wtime() - region_start;
    }
    const double wall = omp_get_wtime() - wall_start;
    // Rip-up history of the regions' sequences (Sec. 4.2).
    for (Worker& w : workers_) {
      for (const size_t v : w.history_bumps) {
        grid_.bumpHistory(v);
      }
      w.history_bumps.clear();
    }
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
    // The nets that continue in the dilated area stay pending with that
    // area's box, at their positions in the net order.
    size_t again = 0;
    for (const auto& a : region_again) {
      again += a.size();
      for (const int id : a) {
        set_box(id);
        stays[id] = 1;
      }
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
                    "utilization {:.0f}%; {} nets inside continue in the "
                    "dilated area later.",
                    stage,
                    partition,
                    regions.size(),
                    routed,
                    next.size(),
                    wall,
                    busy_max,
                    wall > 0 ? 100.0 * busy_sum / (wall * threads) : 0.0,
                    again);
    }
    if (again > 0) {
      for (const int id : next) {
        stays[id] = 1;
      }
      next.clear();
      for (const int id : pending) {
        if (stays[id] != 0) {
          next.push_back(id);
          stays[id] = 0;
        }
      }
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
  saved.input_present = net.input_present;
  for (const Term& term : net.terms) {
    saved.chosen.push_back(term.chosen);
    saved.reserved.push_back(term.reservation_active);
  }
  saved.was_partial = net.partial;
  saved.partial_state = net.partial_state;
  return saved;
}

void BonnRouteDR::restoreRouting(const SavedRouting& saved, Worker* local)
{
  NetData& net = nets_[saved.net];
  if (net.routed || net.partial) {
    ripupNet(net, local);
  }
  net.wires = saved.wires;
  net.vias = saved.vias;
  net.pin_points = saved.pin_points;
  for (size_t t = 0; t < net.terms.size(); t++) {
    net.terms[t].chosen = saved.chosen[t];
    setReservation(net.terms[t], net.id, saved.reserved[t] != 0, local);
  }
  setInputReservation(net, saved.input_reserved, local);
  if (saved.input_present && !net.input_present) {
    setInputShapes(net, true);
  }
  if (saved.was_routed) {
    commitNet(net, local);
  } else {
    net.wires.clear();
    net.vias.clear();
    if (saved.was_partial) {
      commitPartial(net, saved.partial_state, local);
    }
  }
}

bool BonnRouteDR::ripupSequence(const int id,
                                const int depth,
                                const int dilation,
                                std::vector<SavedRouting>& journal,
                                const bool try_plain,
                                Worker& worker,
                                RipupState& state,
                                const br::TileBox* resident)
{
  NetData& net = nets_[id];
  Worker* local = resident != nullptr ? &worker : nullptr;
  journal.push_back(saveRouting(net));
  if (try_plain) {
    const double start = omp_get_wtime();
    const bool ok = connect(net, false, worker, dilation);
    state.plain_seconds += omp_get_wtime() - start;
    state.plain_connects++;
    if (ok) {
      commitNet(net, local);
      setReservations(net, false, local);
      return true;
    }
  }
  if (depth >= params_.ripup_depth) {
    state.outcome = kRipupDepth;
    state.outcome_depth = depth;
    return false;
  }
  const double start = omp_get_wtime();
  const bool found = connect(net, true, worker, dilation);
  state.mode_seconds += omp_get_wtime() - start;
  state.mode_connects++;
  if (!found) {
    state.outcome = kRipupNoPath;
    state.outcome_depth = depth;
    return false;
  }
  std::vector<int> victims;
  collectVictims(net, victims);
  for (const int victim : victims) {
    if (nets_[victim].ripups >= params_.max_ripups) {
      state.outcome = kRipupVictimLimit;
      state.outcome_depth = depth;
      return false;
    }
    // Sec. 3.3: shapes above the pass's level stay (an access path or a
    // post-processing shape of the new route can meet them).
    if (shape_grid_.netRipupLevel(victim) > worker.max_ripup_level) {
      state.outcome = kRipupVictimLevel;
      state.outcome_depth = depth;
      state.victim_level_failures++;
      return false;
    }
    // Sec. 5.1: a region's sequence changes routing only inside the region.
    if (resident != nullptr
        && !resident->contains(
            netAreaTileBox(nets_[victim], params_.retry_dilation))) {
      state.outcome = kRipupVictimOutside;
      state.outcome_depth = depth;
      state.victim_outside_failures++;
      return false;
    }
    if (!prepareVictim(victim)) {
      state.outcome = kRipupVictimFailed;
      state.outcome_depth = depth;
      return false;
    }
  }
  // Remove the victims, add the new route, then reconnect the victims,
  // again possibly with rip-up (Sec. 4.4).
  for (const int victim : victims) {
    journal.push_back(saveRouting(nets_[victim]));
    ripupNet(nets_[victim], local);
    nets_[victim].ripups++;
    state.events++;
    state.victims_by_level[shape_grid_.netRipupLevel(victim)]++;
  }
  commitNet(net, local);
  setReservations(net, false, local);
  for (const int victim : victims) {
    if (!nets_[victim].routed
        && !ripupSequence(victim,
                          depth + 1,
                          dilation,
                          journal,
                          true,
                          worker,
                          state,
                          resident)) {
      return false;
    }
  }
  return true;
}

bool BonnRouteDR::regionRipup(NetData& net,
                              Worker& worker,
                              const br::TileBox& region,
                              const int64_t region_settled_start)
{
  RipupState& state = worker.ripup;
  // Budget from this region's own work only (deterministic).
  const int64_t own
      = worker.settled - region_settled_start - worker.region_ripup_settled;
  if (worker.region_ripup_settled * 100 >= own * params_.region_ripup_effort) {
    state.budget_skips++;
    return false;
  }
  const int64_t settled_before = worker.settled;
  const br::RipupLevel level = worker.max_ripup_level;
  worker.max_ripup_level
      = static_cast<br::RipupLevel>(params_.first_pass_ripup_level);
  // The sequence may change routing anywhere in the region.
  worker.footprint.widen(region);
  std::vector<SavedRouting> journal;
  state.outcome = kRipupNotStarted;
  state.outcome_depth = 0;
  state.sequences++;
  const bool ok
      = ripupSequence(net.id, 0, 0, journal, false, worker, state, &region);
  if (ok) {
    state.routed++;
  } else {
    // Sec. 4.4: a failed sequence restores the original situation.
    for (auto it = journal.rbegin(); it != journal.rend(); ++it) {
      restoreRouting(*it, &worker);
    }
    state.rollbacks++;
  }
  worker.max_ripup_level = level;
  worker.area_deferred = false;
  worker.region_ripup_settled += worker.settled - settled_before;
  state.settled += worker.settled - settled_before;
  return ok;
}


void BonnRouteDR::routeRipupQueue(const std::vector<int>& failed,
                                  const int64_t settled_budget)
{
  std::vector<int> queue = failed;
  std::vector<SavedRouting> journal;
  Worker& worker = workers_.front();
  // Vertices settled by the serial loop's own sequences.
  int64_t serial_settled = 0;
  for (int pass = 0; pass < params_.ripup_passes && !queue.empty(); pass++) {
    // Region-parallel rounds first; the serial loop takes the rest.
    // Incremental calls prepare victims and move input routing serially.
    for (int round = 0;
         !incremental_active_ && round < params_.ripup_partitions;
         round++) {
      std::vector<int> pending;
      for (const int id : queue) {
        if (!nets_[id].routed) {
          pending.push_back(id);
        }
      }
      if (pending.empty()) {
        break;
      }
      parallelRipupRound(pending, pass, round, settled_budget);
    }
    // Sec. 4.4: later passes with higher ripup effort.
    worker.max_ripup_level
        = pass == 0 ? static_cast<br::RipupLevel>(params_.first_pass_ripup_level)
                    : br::kMaxRemovableLevel;
    std::vector<int> unrouted;
    for (const int id : queue) {
      if (nets_[id].routed) {
        continue;
      }
      if (settled_budget > 0 && serial_settled >= settled_budget) {
        unrouted.push_back(id);
        ripup_budget_skipped_++;
        if (nets_[id].ripup_outcome == kRipupNotStarted) {
          nets_[id].ripup_outcome = kRipupBudget;
        }
        continue;
      }
      journal.clear();
      serial_ripup_.outcome = kRipupNotStarted;
      serial_ripup_.outcome_depth = 0;
      const int64_t sequence_start = worker.settled;
      const bool routed = ripupSequence(
          id, 0, pass, journal, false, worker, serial_ripup_, nullptr);
      serial_settled += worker.settled - sequence_start;
      if (routed) {
        continue;
      }
      // The step the net's own rip-up mode search failed at, before the
      // rollback restores the routings.
      const uint8_t step = nets_[id].left_step;
      const int term = nets_[id].left_term;
      // Sec. 4.4: a failed sequence restores the original situation.
      for (auto it = journal.rbegin(); it != journal.rend(); ++it) {
        restoreRouting(*it);
      }
      serial_ripup_.rollbacks++;
      unrouted.push_back(id);
      // A failure below the net's own call is a victim's.
      nets_[id].ripup_outcome = serial_ripup_.outcome_depth > 0
                                    ? uint8_t{kRipupVictimFailed}
                                    : serial_ripup_.outcome;
      nets_[id].left_step = step;
      nets_[id].left_term = term;
    }
    queue.swap(unrouted);
  }
  worker.max_ripup_level = br::kMaxRemovableLevel;
  for (const int id : queue) {
    NetData& net = nets_[id];
    net.wires.clear();
    net.vias.clear();
    net.failed = true;
  }
}

int64_t BonnRouteDR::parallelRipupRound(const std::vector<int>& queue,
                                        const int pass,
                                        const int round,
                                        const int64_t settled_budget)
{
  const double wall_start = omp_get_wtime();
  std::vector<std::array<int, 4>> boxes(nets_.size());
  const int dilation = params_.retry_dilation + pass;
  for (const int id : queue) {
    const br::TileBox b = netAreaTileBox(nets_[id], dilation);
    boxes[id] = {b.xlo, b.ylo, b.xhi, b.yhi};
  }
  const int threads = static_cast<int>(workers_.size());
  std::vector<std::array<int, 4>> regions;
  partitionRegions(queue, boxes, 2 * threads, round % 2 == 1, regions);
  std::vector<std::vector<int>> members(regions.size());
  for (const int id : queue) {
    const auto& b = boxes[id];
    for (size_t r = 0; r < regions.size(); r++) {
      const auto& g = regions[r];
      if (b[0] >= g[0] && b[1] >= g[1] && b[2] <= g[2] && b[3] <= g[3]) {
        members[r].push_back(id);
        break;
      }
    }
  }
  std::vector<int64_t> settled(regions.size(), 0);
  std::vector<int64_t> sequences(regions.size(), 0);
  std::vector<int64_t> routed(regions.size(), 0);
  std::vector<int64_t> rollbacks(regions.size(), 0);
  std::vector<int64_t> stops(regions.size(), 0);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
  for (int r = 0; r < static_cast<int>(regions.size()); r++) {
    if (members[r].empty()) {
      continue;
    }
    Worker& w = workers_[omp_get_thread_num()];
    const br::TileBox region{
        regions[r][0], regions[r][1], regions[r][2], regions[r][3]};
    const br::AccessFootprintScope footprint_scope(
        params_.partition_check ? &w.footprint : nullptr);
    // The region's share of the serial budget (deterministic).
    const int64_t budget
        = settled_budget > 0
              ? settled_budget * static_cast<int64_t>(members[r].size())
                    / static_cast<int64_t>(queue.size())
              : 0;
    const int64_t start = w.settled;
    const br::RipupLevel level = w.max_ripup_level;
    w.max_ripup_level
        = pass == 0 ? static_cast<br::RipupLevel>(params_.first_pass_ripup_level)
                    : br::kMaxRemovableLevel;
    std::vector<SavedRouting> journal;
    for (const int id : members[r]) {
      NetData& net = nets_[id];
      if (net.routed) {
        continue;
      }
      if (budget > 0 && w.settled - start >= budget) {
        stops[r]++;
        break;
      }
      if (params_.partition_check) {
        const auto& b = boxes[id];
        beginFootprint(w, net, kPhaseRegion, {b[0], b[1], b[2], b[3]});
        w.footprint.widen(region);
      }
      w.limit = &region;
      w.area_deferred = false;
      journal.clear();
      w.ripup.outcome = kRipupNotStarted;
      w.ripup.outcome_depth = 0;
      sequences[r]++;
      if (ripupSequence(id, 0, pass, journal, false, w, w.ripup, &region)) {
        routed[r]++;
      } else {
        const uint8_t step = net.left_step;
        const int term = net.left_term;
        for (auto it = journal.rbegin(); it != journal.rend(); ++it) {
          restoreRouting(*it, &w);
        }
        rollbacks[r]++;
        net.ripup_outcome = w.ripup.outcome_depth > 0
                                ? uint8_t{kRipupVictimFailed}
                                : w.ripup.outcome;
        net.left_step = step;
        net.left_term = term;
      }
      w.limit = nullptr;
      w.area_deferred = false;
      if (params_.partition_check) {
        endFootprint(w, net, &region);
      }
    }
    w.max_ripup_level = level;
    settled[r] = w.settled - start;
  }
  // Rip-up history of the regions' sequences (Sec. 4.2).
  for (Worker& w : workers_) {
    for (const size_t v : w.history_bumps) {
      grid_.bumpHistory(v);
    }
    w.history_bumps.clear();
  }
  int64_t total = 0;
  for (size_t r = 0; r < regions.size(); r++) {
    total += settled[r];
    parallel_ripup_sequences_ += sequences[r];
    parallel_ripup_routed_ += routed[r];
    parallel_ripup_rollbacks_ += rollbacks[r];
    parallel_ripup_budget_stops_ += stops[r];
  }
  parallel_ripup_rounds_++;
  parallel_ripup_seconds_ += omp_get_wtime() - wall_start;
  return total;
}

int BonnRouteDR::blockerKinds(const RShape& shape, const int net) const
{
  int kinds = 0;
  auto scan = [&](const int g, const RShape& s) {
    const br::Coord reach = shape_classes_.reach(g);
    const br::Rect r = brRect(s.rect);
    br::GridWindow pieces;
    shape_grid_.collect(g, r.bloated(reach, reach), pieces);
    const br::PlacedShape placed{r, shapeClass(s, g), 3};
    for (const br::GridPiece& p : pieces.pieces) {
      if (p.owner == net && net >= 0 && br::isOwned(p.type)) {
        continue;
      }
      if (!br::shapesViolate(shape_classes_,
                             placed,
                             br::PlacedShape{p.rect, p.cls, p.line_end_ends})) {
        continue;
      }
      switch (p.type) {
        case br::ShapeType::kBlockage:
          kinds |= 1;
          break;
        case br::ShapeType::kPin:
          kinds |= 2;
          break;
        case br::ShapeType::kReservation:
          kinds |= 4;
          break;
        default:
          kinds |= 8;
          break;
      }
    }
  };
  const int g = gridLayer(shape.layer);
  if (g < 0) {
    return 0;
  }
  scan(g, shape);
  const int above = projectionLayer(g);
  if (above >= 0) {
    RShape projection = shape;
    projection.role = br::ShapeRole::kCutProjection;
    scan(above, projection);
  }
  return kinds;
}

BonnRouteDR::AccessBlock BonnRouteDR::accessBlock(const NetData& net,
                                                  const Term& term) const
{
  AccessBlock best = kAccessNoOption;
  bool any = false;
  std::vector<RShape> via;
  for (const AccessOption& option : term.options) {
    int kinds = 0;
    // The wire at the option's track vertex.
    const auto [x, y] = graph_.point(option.vertex);
    const frCoord half = layers_[option.z].width / 2;
    kinds |= blockerKinds(RShape{layers_[option.z].fr,
                                 odb::Rect(x - half, y - half, x + half, y + half),
                                 br::ShapeRole::kPreferredWire},
                          net.id);
    for (const OutWire& wire : option.wires) {
      kinds |= blockerKinds(wireShape(wire), net.id);
    }
    for (size_t k = 0; k < option.vias.size(); k++) {
      via.clear();
      viaShapes(option.vias[k], via);
      for (size_t s = 0; s < via.size(); s++) {
        int v = blockerKinds(via[s], net.id);
        // Pin access checked the pin-layer pad and the cut of the access
        // point's own via against the fixed surroundings (optionLegal).
        if (static_cast<int>(k) == option.pin_via && s + 1 < via.size()) {
          v &= 4 | 8;
        }
        kinds |= v;
      }
    }
    const AccessBlock block = kinds == 0             ? kAccessFeasible
                              : (kinds & (1 | 2)) != 0 ? kAccessFixed
                              : (kinds & 4) != 0       ? kAccessReservation
                                                       : kAccessRouted;
    best = any ? std::min(best, block) : block;
    any = true;
  }
  return best;
}

void BonnRouteDR::reportLeftNets() const
{
  static constexpr std::array<const char*, kNumLeftSteps> kSteps{
      "none", "access", "no path", "search budget"};
  static constexpr std::array<const char*, kNumRipupOutcomes> kOutcomes{
      "not started",
      "effort budget",
      "depth",
      "no path with rip-up",
      "victim at its rip-up limit",
      "victim not routed again",
      "victim above the rip-up level",
      "victim outside the region"};
  static constexpr std::array<const char*, kNumAccessBlocks> kBlocks{
      "no option",
      "feasible at the end",
      "routed shapes",
      "reservations",
      "fixed shapes or pins"};
  int64_t left = 0;
  int64_t guard_contact = 0;
  int64_t guard_on_track = 0;
  int64_t guard_missing_contact = 0;
  std::array<int64_t, kNumLeftSteps> by_step{};
  std::array<int64_t, kNumRipupOutcomes> by_outcome{};
  std::array<int64_t, kNumAccessBlocks> by_block{};
  const bool trace = logger_->debugCheck(utl::DRT, "bonnroute_left", 1);
  for (const NetData& net : nets_) {
    if (net.routed || net.skip) {
      continue;
    }
    left++;
    AccessBlock block = kAccessNoOption;
    const bool access = net.left_guard == kGuardNone
                        && net.left_step == kLeftAccess && net.left_term >= 0
                        && net.left_term < static_cast<int>(net.terms.size());
    if (net.left_guard == kGuardPinContact) {
      guard_contact++;
    } else if (net.left_guard == kGuardOnTrack) {
      guard_on_track++;
    } else if (net.left_guard == kGuardMissingContact) {
      guard_missing_contact++;
    } else {
      by_step[net.left_step]++;
      by_outcome[net.ripup_outcome]++;
      if (access) {
        block = accessBlock(net, net.terms[net.left_term]);
        by_block[block]++;
      }
    }
    if (trace) {
      logger_->debug(
          utl::DRT,
          "bonnroute_left",
          "left net {} ({}) terms {} hpwl {} guide_tiles {} guard {} step {} "
          "term {} access {} ripup {} ripups {} dilated {}",
          net.id,
          net.net->getName(),
          net.terms.size(),
          net.hpwl,
          net.guide_tiles.size(),
          static_cast<int>(net.left_guard),
          kSteps[net.left_step],
          net.left_term,
          access ? kBlocks[block] : "-",
          kOutcomes[net.ripup_outcome],
          net.ripups,
          net.dilated_area ? 1 : 0);
    }
  }
  logger_->report(
      "BonnRoute DR left nets: {}; output guards: {} pin contact, {} "
      "on-track-only layer, {} without a FlexDR pin contact; last failed "
      "step of the others: {} access ({} pin without an option, {} blocked "
      "by fixed shapes or pins, {} by "
      "reservations, {} by routed shapes, {} feasible at the end), {} no "
      "path, {} search budget, {} none; their rip-up sequence: {} not "
      "started, {} over the effort budget, {} no path with rip-up, {} victim "
      "at its rip-up limit, {} victim not routed again, {} depth, {} victim "
      "above the rip-up level, {} victim outside the region.",
      left,
      guard_contact,
      guard_on_track,
      guard_missing_contact,
      by_step[kLeftAccess],
      by_block[kAccessNoOption],
      by_block[kAccessFixed],
      by_block[kAccessReservation],
      by_block[kAccessRouted],
      by_block[kAccessFeasible],
      by_step[kLeftNoPath],
      by_step[kLeftSearchBudget],
      by_step[kLeftNone],
      by_outcome[kRipupNotStarted],
      by_outcome[kRipupBudget],
      by_outcome[kRipupNoPath],
      by_outcome[kRipupVictimLimit],
      by_outcome[kRipupVictimFailed],
      by_outcome[kRipupDepth],
      by_outcome[kRipupVictimLevel],
      by_outcome[kRipupVictimOutside]);
}

void BonnRouteDR::writeNets()
{
  frTechObject* tech = design_->getTech();
  std::vector<std::pair<frLayerNum, odb::Rect>> notch_fill_rects;
  // (net, index in notch_fill_rects) of every fill written.
  std::vector<std::pair<int, size_t>> written_fills;
  for (NetData& data : nets_) {
    if (!data.routed || (!write_mask_.empty() && !write_mask_[data.id])) {
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
      data.left_guard = kGuardPinContact;
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
      data.left_guard = kGuardOnTrack;
      data.routed = false;
      data.wires.clear();
      data.vias.clear();
      continue;
    }
    // Every pin of the net needs a contact FlexDR recognises; a pin without
    // one stops FlexDR's connectivity check (DRT-0206). Leave the net to
    // FlexDR instead.
    if (const int missing = pinsWithoutContact(data, wires)) {
      uncontacted_nets_++;
      uncontacted_pins_ += missing;
      data.left_guard = kGuardMissingContact;
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
    if (params_.notch_fill) {
      const size_t first_fill = notch_fill_rects.size();
      fillNotches(net, data.id, notch_fill_rects);
      for (size_t k = first_fill; k < notch_fill_rects.size(); k++) {
        written_fills.emplace_back(data.id, k);
      }
    }
    net->setHasInitialRouting(true);
    net->setModified(true);
  }
  // A routing space kept for incremental calls holds the notch fills as
  // routed shapes of their net (patch wires, as routingShapes reads them);
  // inserted after all nets are written, so the written routing is the
  // same with and without a kept space.
  if (record_fixed_ || incremental_active_) {
    for (const auto& [id, k] : written_fills) {
      const RShape s
          = patchShape(notch_fill_rects[k].first, notch_fill_rects[k].second);
      nets_[id].shapes.push_back(s);
      changeShape(s, routedType(s.role), id, true, nullptr);
    }
  }
}

void BonnRouteDR::fillNotches(
    frNet* net,
    const int id,
    std::vector<std::pair<frLayerNum, odb::Rect>>& fills)
{
  // drt checks the spacing between separated pieces of one net's metal on
  // a layer (same-net Metal Spacing): a via pad next to a jog of its own
  // wire, two via pads one pitch apart, a pad next to a pin shape. Filling
  // the gap joins the pieces into one polygon.
  struct Piece
  {
    frLayerNum layer;
    odb::Rect rect;
    bool routed;
  };
  std::vector<Piece> pieces;
  std::vector<std::pair<frLayerNum, odb::Point>> anchors;
  for (const auto& shape : net->getShapes()) {
    if (shape->typeId() != frcPathSeg) {
      continue;
    }
    const auto* seg = static_cast<const frPathSeg*>(shape.get());
    pieces.push_back({seg->getLayerNum(), seg->getBBox(), true});
    const auto [b, e] = seg->getPoints();
    anchors.emplace_back(seg->getLayerNum(), b);
    anchors.emplace_back(seg->getLayerNum(), e);
  }
  for (const auto& via : net->getVias()) {
    const frViaDef* def = via->getViaDef();
    pieces.push_back({def->getLayer1Num(), via->getLayer1BBox(), true});
    pieces.push_back({def->getLayer2Num(), via->getLayer2BBox(), true});
    anchors.emplace_back(def->getLayer1Num(), via->getOrigin());
    anchors.emplace_back(def->getLayer2Num(), via->getOrigin());
  }
  for (const auto& patch : net->getPatchWires()) {
    pieces.push_back({patch->getLayerNum(), patch->getBBox(), true});
  }
  std::vector<frRect> pin_shapes;
  for (frInstTerm* it : net->getInstTerms()) {
    it->getShapes(pin_shapes);
  }
  for (frBTerm* bt : net->getBTerms()) {
    bt->getShapes(pin_shapes);
  }
  for (const frRect& r : pin_shapes) {
    pieces.push_back({r.getLayerNum(), r.getBBox(), false});
  }
  const size_t original = pieces.size();
  for (size_t a = 0; a < original; a++) {
    if (!pieces[a].routed) {
      continue;
    }
    const frLayerNum layer = pieces[a].layer;
    const int z = layer >= 0 && layer < static_cast<int>(z_of_fr_.size())
                      ? z_of_fr_[layer]
                      : -1;
    if (z < 0) {
      continue;
    }
    for (size_t b = 0; b < original; b++) {
      if (b == a || pieces[b].layer != layer || (pieces[b].routed && b < a)) {
        continue;
      }
      const odb::Rect& ra = pieces[a].rect;
      const odb::Rect& rb = pieces[b].rect;
      const int gx = std::max(rb.xMin() - ra.xMax(), ra.xMin() - rb.xMax());
      const int gy = std::max(rb.yMin() - ra.yMax(), ra.yMin() - rb.yMax());
      odb::Rect fill;
      int gap = 0;
      int prl = 0;
      if (gx > 0 && gy < 0) {
        gap = gx;
        prl = -gy;
        fill.init(std::min(ra.xMax(), rb.xMax()),
                  std::max(ra.yMin(), rb.yMin()),
                  std::max(ra.xMin(), rb.xMin()),
                  std::min(ra.yMax(), rb.yMax()));
      } else if (gy > 0 && gx < 0) {
        gap = gy;
        prl = -gx;
        fill.init(std::max(ra.xMin(), rb.xMin()),
                  std::min(ra.yMax(), rb.yMax()),
                  std::min(ra.xMax(), rb.xMax()),
                  std::max(ra.yMin(), rb.yMin()));
      } else {
        continue;
      }
      if (gap >= spacing(z, ra.minDXDY(), rb.minDXDY(), prl)) {
        continue;
      }
      // The joined metal is as wide as the overlap: below the minimum width
      // the fill would be a Min Width violation (jpeg: 143 new markers).
      if (prl < static_cast<int>(
              design_->getTech()->getLayer(layer)->getMinWidth())) {
        continue;
      }
      // Already joined by another piece.
      if (std::ranges::any_of(pieces, [&](const Piece& c) {
            return c.layer == layer && c.rect.contains(fill);
          })) {
        continue;
      }
      // Legal against other nets: the routing space (their routing, pins
      // and blockages) and the fills written for earlier nets.
      const bool along_x = fill.dx() >= fill.dy();
      RShape shape{layer, fill};
      shape.role = (fill.dx() == fill.dy() || along_x == layers_[z].horizontal)
                       ? br::ShapeRole::kPreferredWire
                       : br::ShapeRole::kNonPreferredWire;
      if (shapeConflicts(shape, id, nullptr)) {
        continue;
      }
      const frCoord need = spacing(z, fill.minDXDY(), fill.minDXDY(), 0);
      if (std::ranges::any_of(fills, [&](const auto& other) {
            if (other.first != layer) {
              return false;
            }
            const odb::Rect& o = other.second;
            const int dx
                = std::max(o.xMin() - fill.xMax(), fill.xMin() - o.xMax());
            const int dy
                = std::max(o.yMin() - fill.yMax(), fill.yMin() - o.yMax());
            return std::max(dx, dy) < need;
          })) {
        continue;
      }
      // Anchored at a via origin or a wire end of the net on this layer
      // (FlexDR's connectivity check keeps patches only there).
      const odb::Point center((fill.xMin() + fill.xMax()) / 2,
                              (fill.yMin() + fill.yMax()) / 2);
      const std::pair<frLayerNum, odb::Point>* anchor = nullptr;
      int64_t best = std::numeric_limits<int64_t>::max();
      for (const auto& candidate : anchors) {
        if (candidate.first != layer) {
          continue;
        }
        const int64_t d
            = odb::Point::manhattanDistance(candidate.second, center);
        if (d < best) {
          best = d;
          anchor = &candidate;
        }
      }
      if (anchor == nullptr) {
        continue;
      }
      auto pwire = std::make_unique<frPatchWire>();
      pwire->setLayerNum(layer);
      pwire->setOrigin(anchor->second);
      odb::Rect offset = fill;
      offset.moveDelta(-anchor->second.x(), -anchor->second.y());
      pwire->setOffsetBox(offset);
      net->addPatchWire(std::move(pwire));
      pieces.push_back({layer, fill, true});
      fills.emplace_back(layer, fill);
      notch_fills_++;
    }
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
  if (params_.use_intervals && params_.corridor_cost > 0) {
    // The interval search has no corridor price (BrCore.h).
    logger_->warn(utl::DRT,
                  9131,
                  "BonnRoute: corridor_cost needs the vertex search "
                  "(use_intervals 0); corridor_cost is ignored.");
    params_.corridor_cost = 0;
  }
  trace_connections_ = logger_->debugCheck(utl::DRT, "bonnroute_connect", 1);
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
  if (params_.planar_first) {
    buildPinViaLegality();
  }
  buildPartitionHalos();
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
  markWideNets();
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
  // Sec. 3.3: rip-up level of the critical nets, set before any of their
  // shapes is committed; the parallel phases only read the level table.
  if (params_.critical_ripup_level > 1 && !critical.empty()) {
    shape_grid_.setNetRipupLevel(static_cast<int32_t>(nets_.size()) - 1, 1);
    for (const int id : critical) {
      shape_grid_.setNetRipupLevel(
          id, static_cast<br::RipupLevel>(params_.critical_ripup_level));
    }
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
  if (params_.region_ripup_effort > 0 && params_.verbose) {
    RipupState total;
    for (const Worker& w : workers_) {
      const RipupState& s = w.ripup;
      total.sequences += s.sequences;
      total.routed += s.routed;
      total.budget_skips += s.budget_skips;
      total.events += s.events;
      total.rollbacks += s.rollbacks;
      total.victim_level_failures += s.victim_level_failures;
      total.victim_outside_failures += s.victim_outside_failures;
      total.mode_connects += s.mode_connects;
      total.plain_connects += s.plain_connects;
      total.mode_seconds += s.mode_seconds;
      total.plain_seconds += s.plain_seconds;
      total.settled += s.settled;
    }
    logger_->report(
        "BonnRoute DR region rip-up: {} sequences, {} nets routed, {} rolled "
        "back, {} over the region budget; {} rip-ups; failed at a victim "
        "above the level {}, outside the region {}; {} connects with rip-up "
        "{:.2f} thread-s, {} without {:.2f} thread-s; {} settled vertices.",
        total.sequences,
        total.routed,
        total.rollbacks,
        total.budget_skips,
        total.events,
        total.victim_level_failures,
        total.victim_outside_failures,
        total.mode_connects,
        total.mode_seconds,
        total.plain_connects,
        total.plain_seconds,
        total.settled);
  }
  if (params_.partition_check) {
    reportFootprints();
  }
  // Unresolved nets get the serial rip-up sequences of Sec. 4.4 (ripup,
  // or final_ripup for the nets the parallel phases left). A net left to
  // FlexDR is routed by its iteration 0 over the fixed BonnRoute routing,
  // which is where its shorts come from.
  // Many failed nets point at rules the routing space does not model
  // (aes_asap7: 1,698 of ~17k nets; cut spacing tables, EOL keep-out):
  // rip-up then routes nets FlexDR must clean up anyway (DR stage +3.8%
  // to +162% for budgets 5% to unlimited), so the final phase is skipped
  // above final_ripup_max_share percent of the routed nets.
  const bool final_phase = params_.final_ripup
                           && static_cast<int64_t>(failed_first.size()) * 100
                                  <= static_cast<int64_t>(order.size())
                                         * params_.final_ripup_max_share;
  if (params_.final_ripup && !final_phase && params_.verbose) {
    logger_->report(
        "BonnRoute DR rip-up phase skipped: {} of {} nets failed (more than "
        "{}%).",
        failed_first.size(),
        order.size(),
        params_.final_ripup_max_share);
  }
  if (params_.ripup || final_phase) {
    const double ripup_start = timer.elapsed();
    const size_t before = failed_first.size();
    // final_ripup_effort: the serial phase may settle this percentage of
    // the vertices the parallel phases settled (deterministic). Smaller
    // nets first, so the budget goes to the most nets.
    int64_t budget = 0;
    if (!params_.ripup && params_.final_ripup_effort > 0) {
      int64_t parallel_settled = 0;
      for (const Worker& w : workers_) {
        parallel_settled += w.settled;
      }
      budget = parallel_settled * params_.final_ripup_effort / 100;
      std::ranges::stable_sort(failed_first, [&](int a, int b) {
        return nets_[a].hpwl < nets_[b].hpwl;
      });
    }
    routeRipupQueue(failed_first, budget);
    const auto still = std::ranges::count_if(
        failed_first, [&](const int id) { return !nets_[id].routed; });
    if (params_.verbose) {
      logger_->report(
          "BonnRoute DR rip-up phase: {} nets, {} routed, {} left ({} over "
          "the effort budget); {} rip-ups, {} rollbacks; {:.2f}s ({} connects "
          "without rip-up {:.2f}s, {} with rip-up {:.2f}s).",
          before,
          static_cast<int64_t>(before) - still,
          still,
          ripup_budget_skipped_,
          serial_ripup_.events,
          serial_ripup_.rollbacks,
          timer.elapsed() - ripup_start,
          serial_ripup_.plain_connects,
          serial_ripup_.plain_seconds,
          serial_ripup_.mode_connects,
          serial_ripup_.mode_seconds);
      if (params_.ripup_partitions > 0) {
        logger_->report(
            "BonnRoute DR parallel rip-up: {} rounds, {} sequences, {} nets "
            "routed, {} rolled back, {} regions stopped at their budget; "
            "{:.2f}s.",
            parallel_ripup_rounds_,
            parallel_ripup_sequences_,
            parallel_ripup_routed_,
            parallel_ripup_rollbacks_,
            parallel_ripup_budget_stops_,
            parallel_ripup_seconds_);
      }
      logger_->report(
          "BonnRoute DR rip-up victims by level: 1: {}, 2: {}, 3: {}, 4: {}, "
          "5: {}, 6: {}; critical_ripup_level {}, first_pass_ripup_level {}; "
          "{} sequences failed at a victim above the pass's level.",
          serial_ripup_.victims_by_level[1],
          serial_ripup_.victims_by_level[2],
          serial_ripup_.victims_by_level[3],
          serial_ripup_.victims_by_level[4],
          serial_ripup_.victims_by_level[5],
          serial_ripup_.victims_by_level[6],
          params_.critical_ripup_level,
          params_.first_pass_ripup_level,
          serial_ripup_.victim_level_failures);
    }
  }
  if (params_.keep_partial) {
    int64_t commits = 0;
    int64_t resumes = 0;
    int64_t tail_commits = 0;
    for (const Worker& w : workers_) {
      commits += w.partial_commits;
      resumes += w.partial_resumes;
      tail_commits += w.tail_partial_commits;
    }
    int64_t removed = 0;
    int64_t completed = 0;
    for (NetData& net : nets_) {
      if (net.partial && !net.routed) {
        ripupNet(net);
        removed++;
      }
    }
    for (const NetData& net : nets_) {
      completed += net.routed && !net.skip;
    }
    if (params_.verbose) {
      logger_->report(
          "BonnRoute DR partial routing: {} partial commits in the regions, "
          "{} in the tail, {} connects resumed a partial net; {} nets still "
          "partial at the end, removed and left to FlexDR whole; {} nets "
          "routed.",
          commits,
          tail_commits,
          resumes,
          removed,
          completed);
    }
  }
  if (params_.collapse_bridges) {
    const double start = timer.elapsed();
    const int64_t collapsed = collapseBridges();
    if (params_.verbose) {
      logger_->report(
          "BonnRoute DR bridges: {} collapsed into wrong-way wires ({} vias "
          "fewer); kept: {} with other wires on the bridge, {} with "
          "patches, {} with other vias or pin points, {} with a via on a "
          "pin of the net, {} illegal against other nets ({} only routed "
          "shapes, {} fixed shapes) (of {} candidate checks).",
          bridge_stats_.collapsed,
          2 * bridge_stats_.collapsed,
          bridge_stats_.other_wire,
          bridge_stats_.patch,
          bridge_stats_.other_via_or_pin,
          bridge_stats_.pin_contact,
          bridge_stats_.illegal,
          bridge_stats_.illegal_routed,
          bridge_stats_.illegal_fixed,
          bridge_stats_.candidates);
      if (params_.collapse_lower_bridges) {
        logger_->report(
            "BonnRoute DR lower bridges: {} collapsed into wrong-way wires ({} "
            "vias "
            "fewer); kept: {} with other wires on the bridge, {} with "
            "patches, {} with other vias or pin points, {} with a via on a "
            "pin of the net, {} illegal against other nets ({} only routed "
            "shapes, {} fixed shapes) (of {} candidate checks).",
            lower_bridge_stats_.collapsed,
            2 * lower_bridge_stats_.collapsed,
            lower_bridge_stats_.other_wire,
            lower_bridge_stats_.patch,
            lower_bridge_stats_.other_via_or_pin,
            lower_bridge_stats_.pin_contact,
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
  if (params_.start_check > 0) {
    checkPathStarts();
  }
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
  int64_t held_pops = 0;
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
    held_pops += w.held_pops;
    settled += w.settled;
    if (w.interval_search) {
      intervals += w.interval_search->intervals();
      labels += w.interval_search->labels();
    }
  }
  if (params_.verbose) {
    // Sec. 5.1: halos of the nets' boxes.
    std::vector<int64_t> nets_by_top(graph_.numLayers(), 0);
    std::vector<int64_t> nets_by_access(graph_.numLayers(), 0);
    std::vector<int64_t> nets_by_halo;
    for (const NetData& net : nets_) {
      if (net.skip) {
        continue;
      }
      nets_by_top[netTopLayer(net)]++;
      nets_by_access[netAccessTopLayer(net)]++;
      const size_t halo = static_cast<size_t>(netHaloTiles(net));
      if (nets_by_halo.size() <= halo) {
        nets_by_halo.resize(halo + 1, 0);
      }
      nets_by_halo[halo]++;
    }
    std::string by_layer;
    int64_t wide_nets = 0;
    for (const NetData& net : nets_) {
      wide_nets += !net.skip && net.wide_shapes;
    }
    for (int z = 0; z < graph_.numLayers(); z++) {
      by_layer += fmt::format(
          "{}{} area {} (narrow {}) tiles ({} nets) access {} (narrow {}) "
          "tiles ({} nets)",
          z == 0 ? "" : ", ",
          layers_[z].layer->getName(),
          area_halo_tiles_[z],
          narrow_area_halo_tiles_[z],
          nets_by_top[z],
          access_halo_tiles_[z],
          narrow_access_halo_tiles_[z],
          nets_by_access[z]);
    }
    std::string by_halo;
    for (size_t h = 0; h < nets_by_halo.size(); h++) {
      if (nets_by_halo[h] > 0) {
        by_halo += fmt::format(" {}:{}", h, nets_by_halo[h]);
      }
    }
    // narrowReach() grows with the classes interned during routing; the
    // halos assume it stays within narrowRuleReach().
    int reach_layers = 0;
    for (int g = 0; g < shape_classes_.numLayers(); g++) {
      reach_layers += shape_classes_.narrowReach(g)
                      > shape_classes_.narrowRuleReach(g);
    }
    logger_->report(
        "BonnRoute DR partition halo by top layer: {}; nets by halo tiles{}; "
        "{} nets with wide shapes use the rule distance; {} wide shape "
        "updates for other nets by partition threads, {} layers whose "
        "narrow class distance exceeds its rule bound.",
        by_layer,
        by_halo,
        wide_nets,
        unexpected_wide_updates_.load(),
        reach_layers);
    // Sec. 4.4 / 5.1: nets that left their guide area.
    int64_t dilated = 0;
    int64_t dilated_routed = 0;
    for (const NetData& net : nets_) {
      if (net.dilated_area) {
        dilated++;
        dilated_routed += net.routed;
      }
    }
    int64_t area_deferrals = 0;
    for (const Worker& w : workers_) {
      area_deferrals += w.area_deferrals;
    }
    logger_->report(
        "BonnRoute DR routing areas: {} nets continued in the area dilated "
        "by {} tiles ({} of them routed); {} times a net's next area did "
        "not fit its partition region and the net continued in a later "
        "partition or the tail.",
        dilated,
        params_.retry_dilation,
        dilated_routed,
        area_deferrals);
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
      int64_t check_pops = 0;
      int64_t check_settled = 0;
      int64_t check_plain_pops = 0;
      int64_t check_plain_settled = 0;
      for (const Worker& w : workers_) {
        fallbacks += w.potential_fallbacks;
        checks += w.potential_checks;
        mismatches += w.potential_mismatches;
        check_pops += w.potential_check_pops;
        check_settled += w.potential_check_settled;
        check_plain_pops += w.potential_check_plain_pops;
        check_plain_settled += w.potential_check_plain_settled;
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
      if (checks > 0) {
        logger_->report(
            "BonnRoute DR: potential check: {} searches repeated with pi_H in "
            "the same state; heap pops {} with the potential, {} with pi_H; "
            "settled vertices {} with the potential, {} with pi_H.",
            checks,
            check_pops,
            check_plain_pops,
            check_settled,
            check_plain_settled);
      }
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
    reportLeftNets();
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
        serial_ripup_.events,
        serial_ripup_.rollbacks,
        path_region_searches,
        timer.elapsed());
    if (!params_.use_intervals) {
      logger_->report(
          "BonnRoute DR vertex search queue: {} of {} heap pops answered by "
          "the last pushed item ({:.1f}%).",
          held_pops,
          heap_pops,
          heap_pops > 0 ? 100.0 * held_pops / heap_pops : 0.0);
    }
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
      if (params_.offtrack_access_report > 0) {
        reportAccessFailures();
      }
    }
    logger_->report(
        "BonnRoute DR contact guard: {} nets with {} truncated pin contacts "
        "away from access points left to FlexDR; {} via pads on pins away "
        "from access points (not FlexDR terms, kept); {} via pads at access "
        "points marked connected; {} nets with wrong-way or off-track wires "
        "on on-track-only layers left to FlexDR; {} same-net notch fills.",
        guard_nets_,
        guard_contacts_,
        guard_via_contacts_,
        connected_vias_,
        wrong_way_nets_,
        notch_fills_);
    logger_->report(
        "BonnRoute DR pin contacts: {} nets with {} pins without a FlexDR "
        "pin contact (truncated wire end or connected via) left to FlexDR.",
        uncontacted_nets_,
        uncontacted_pins_);
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
