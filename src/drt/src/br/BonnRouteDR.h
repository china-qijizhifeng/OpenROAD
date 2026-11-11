// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// BonnRoute detailed routing on the drt database (Gester et al., TODAES
// 18(2), 2013, Sec. 3-4). Runs after pin access, guide processing and track
// assignment; writes complete net routings as initial routing, which
// FlexDR then cleans up (the paper's BR + external DRC clean-up flow).

#pragma once

#include <array>
#include <atomic>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "br/BrCore.h"
#include "br/BrDistanceRules.h"
#include "br/BrFastGrid.h"
#include "br/BrGeometry.h"
#include "br/BrPathCleanup.h"
#include "br/BrPinAccess.h"
#include "br/BrPotentialOracle.h"
#include "br/BrSameNetRules.h"
#include "br/BrShapeGrid.h"
#include "br/BrWireModel.h"
#include "frBaseTypes.h"
#include "frRTree.h"
#include "odb/geom.h"

namespace odb {
class dbBlock;
}

namespace utl {
class Logger;
}

namespace drt {

class frDesign;
class frNet;
class frInst;
class frInstTerm;
class frBTerm;
class frLayer;
class frViaDef;
struct RouterConfiguration;

struct BonnRouteDRParams
{
  // Via cost in pitches of the lower layer (gamma of Sec. 4.1).
  int via_cost = 4;
  // Jog cost factor (beta of Sec. 4.1).
  int jog_factor = 4;
  // Rip-up penalty per vertex, in pitches of the bottom layer.
  int ripup_cost = 40;
  int max_ripups = 3;
  // Off-track pin access: largest jog to a track, in pitches.
  int max_access_jog = 2;
  int max_access_options = 8;
  // Extra cost of a non-preferred access point, in pitches.
  int pin_cost = 2;
  // Routing area: guide tiles plus this many adjacent layers / tiles.
  int layer_dilation = 1;
  int retry_dilation = 2;
  bool use_intervals = false;
  bool ripup = false;
  // Sec. 4.4: search area of a connection = GR tiles of a shortest path
  // along the guides from the source component to another component.
  bool component_regions = true;
  // Sec. 5.1 partition sequence; 0 selects batches of nets with disjoint
  // bounding boxes. partition_regions = 0 uses twice the thread count.
  bool partitions = true;
  int partition_regions = 0;
  // Sec. 5.1: a critical net routing step before the partitions; clock
  // nets are routed first, serially, on the whole chip.
  bool critical_nets = true;
  // Interval property (ii) also for vias (Sec. 4.1); 0 relaxes vias per
  // vertex and lets intervals ignore via availability.
  bool interval_via_property = true;
  // Sec. 4.4: bound on the depth of a rip-up sequence; a failed sequence
  // restores the original routing.
  int ripup_depth = 3;
  // Later passes over unrouted nets use more rip-up effort and a larger
  // routing area (Sec. 4.4, last sentence).
  int ripup_passes = 2;
  // Ahrens et al., arXiv:2111.06169, Sec. 5.2: legal track-assignment wires
  // become reservations of their net, usable by it at cost factor
  // discount_numerator / discount_denominator and blocked for other nets.
  bool reservations = false;
  int discount_numerator = 3;
  int discount_denominator = 4;
  // Lowest track-graph layer (0 = bottom routing layer) with reservations.
  int reserve_min_layer = 0;
  // Future cost: 0 = pi_H of the 2013 paper (Sec. 4.1); arXiv:2111.06169
  // potentials 1 = l1 (Table 3), 2 = simple model (Theorem 3), 3 = general
  // model (Section 4, guide area and discounts), 4 = Gester et al. 2013,
  // Sec. 4.1: the general (blockage-aware) model only for connections
  // whose guide path detours by at least potential_detour_tiles, pi_H
  // otherwise. Reservations or a potential select the vertex search.
  int potential = 0;
  int potential_detour_tiles = 2;
  // Debug: repeat every search that used a potential with pi_H and count
  // cost differences (a feasible potential has none).
  bool potential_check = false;
  // pi_H with the exact via-stack bound lb_via of Sec. 4.1 and the
  // smallest terminal cost per target box; 0 keeps dz * (cheapest via).
  bool exact_future_cost = true;
  // Sec. 4.2 wire spreading: extra cost (pitches of the layer) of a wire
  // vertex beside an occupied track, in GCells whose global routing
  // utilization (dbGCellGrid usage / capacity) is at most
  // spread_utilization percent. Selects the vertex search.
  bool spreading = false;
  int spread_cost = 1;
  int spread_utilization = 70;
  // Off-track pin access (Sec. 3.7, 3.8, 4.3, 4.4): per pin a catalogue of
  // shortest tau-feasible paths on the blockage grid to on-track points
  // within `offtrack_radius` pitches, primary paths per circuit class by
  // destructive bounding. 0 selects the previous access stubs from
  // OpenROAD access points and greedy reservations.
  bool offtrack_access = true;
  int offtrack_radius = 2;
  // Regenerate a pin's catalogue when fewer of its paths remain feasible
  // (Sec. 4.4 "If no or only a small number of access paths remains
  // feasible, we dynamically generate new access paths").
  int offtrack_regenerate = 2;
  // Minimum segment length tau in DBU; 0 derives it from the same-net
  // rules (Sec. 3.7). `offtrack_notch_tau` adds notch spacing + width.
  int offtrack_tau = 0;
  bool offtrack_notch_tau = false;
  // Reuse a circuit class's catalogues and conflict-free solution for all
  // its instances; 0 solves every instance separately.
  bool access_classes = true;
  int access_node_limit = 20000;
  // Weights of the conflict-free solution (Sec. 4.3), in pitches of the
  // bottom track layer: per blocked track, per missing on-track
  // continuation direction; spreading per DBU two endpoints are closer than
  // two pitches.
  int access_track_weight = 1;
  int access_direction_weight = 1;
  int access_spread_weight = 1;
  // Same-net post-processing of each connection (Sec. 4.4).
  bool postprocess = true;
  // Sec. 4.1 searches: per-worker label budget in MB (0: unlimited); a
  // search that would exceed it stops without a path (counted).
  int search_label_mb = 512;
  bool verbose = true;

  bool set(const std::string& name, double value);
};

class BonnRouteDR
{
 public:
  // block: source of the global routing congestion (dbGCellGrid) used by
  // wire spreading; may be nullptr.
  BonnRouteDR(frDesign* design,
              utl::Logger* logger,
              RouterConfiguration* router_cfg,
              const BonnRouteDRParams& params,
              odb::dbBlock* block = nullptr);
  ~BonnRouteDR();

  void main();

 private:
  struct LayerInfo
  {
    frLayerNum fr = 0;
    frLayer* layer = nullptr;
    frCoord width = 0;
    frCoord pitch = 0;
    bool horizontal = true;
    // Via to the next track layer and its shapes relative to the origin.
    const frViaDef* up_via = nullptr;
    odb::Rect up_pad;    // bottom pad on this layer
    odb::Rect up_cut;    // cut on fr + 1
    odb::Rect down_pad;  // top pad of the via coming from below
    frCoord cut_spacing = 0;
    frCoord eol_width = 0;
    frCoord eol_space = 0;
  };
  // Routed geometry of one net, as written to frNet.
  struct OutWire
  {
    frLayerNum layer;
    odb::Point begin;
    odb::Point end;
  };
  struct OutVia
  {
    const frViaDef* def;
    odb::Point point;
  };
  // A shape on a drt layer; role selects its shape class (Sec. 3.1).
  struct RShape
  {
    frLayerNum layer;
    odb::Rect rect;
    br::ShapeRole role = br::ShapeRole::kFixed;
  };
  struct AccessOption
  {
    size_t vertex = 0;
    br::Cost cost = 0;
    int z = 0;
    std::vector<OutWire> wires;
    std::vector<OutVia> vias;
    // Planar access: the stub ends on the pin at `ap` (truncated end).
    bool planar = false;
    odb::Point ap;
    frLayerNum ap_layer = 0;
    // The path from the pin to the vertex as route nodes (drt layer
    // numbers); the post-processing combines it with on-track paths.
    std::vector<br::RouteNode> route;
    // Generated on the blockage grid (legality is checked for all shapes).
    bool offtrack = false;
  };
  struct Term
  {
    std::vector<AccessOption> options;
    std::vector<RShape> reservation;
    bool reservation_active = false;
    int chosen = -1;
    // The pin, for dynamic generation of access paths.
    frInstTerm* inst_term = nullptr;
    frBTerm* bterm = nullptr;
    int primary = -1;
  };
  struct NetData
  {
    frNet* net = nullptr;
    int id = 0;
    std::vector<Term> terms;
    std::vector<OutWire> wires;
    std::vector<OutVia> vias;
    std::vector<RShape> shapes;  // committed shapes
    std::vector<RShape> pins;    // pin shapes in the shape grid
    std::vector<std::pair<frLayerNum, odb::Point>> pin_points;
    // Legal input wires reserved for this net (arXiv:2111.06169 Sec. 5.2)
    // and the track-graph edges they cover.
    std::vector<RShape> input_reservation;
    bool input_reserved = false;
    std::vector<br::DiscountRun> discount_runs;
    // Guide tiles (z, tile x, tile y), deduplicated.
    std::vector<std::array<int, 3>> guide_tiles;
    int ripups = 0;
    bool routed = false;
    bool failed = false;
    bool skip = false;
    int64_t hpwl = 0;
  };

  // Per-thread search state (Sec. 5.1: nets in disjoint areas are routed
  // in parallel).
  struct Worker
  {
    br::Region region;
    std::unique_ptr<br::PathSearch> search;
    std::unique_ptr<br::IntervalSearch> interval_search;
    br::PotentialOracle potential;
    int64_t potential_fallbacks = 0;
    double potential_seconds = 0;
    int64_t potential_checks = 0;
    int64_t potential_mismatches = 0;
    // Guide path of the last buildPathRegion: planar tile steps minus the
    // tile distance of its ends.
    int path_detour = 0;
    // Scratch of the guide-path search for connection regions.
    std::vector<uint32_t> tile_guide;
    std::vector<uint32_t> tile_seen;
    std::vector<uint32_t> tile_target;
    std::vector<int64_t> tile_parent;
    uint32_t tile_epoch = 0;
    // Sec. 5.1: a thread changes the routing space only inside its region;
    // shape grid and fast grid are chunked by global routing tile, so the
    // worker updates them in place. Its configuration transitions are
    // memoised here.
    br::ShapeGridWriteCache cache;
    int64_t fast_grid_queries = 0;
    int64_t checker_queries = 0;
    int64_t aborted_searches = 0;
    br::CleanupStats cleanup;
    int64_t regenerations = 0;
    int64_t regenerated_paths = 0;
    int64_t searches = 0;
    int64_t path_region_searches = 0;
    int64_t heap_pops = 0;
    int64_t settled = 0;
  };
  // Saved routing of a net, for restoring a failed rip-up sequence.
  struct SavedRouting
  {
    int net;
    bool was_routed;
    std::vector<OutWire> wires;
    std::vector<OutVia> vias;
    std::vector<std::pair<frLayerNum, odb::Point>> pin_points;
    std::vector<int> chosen;
    std::vector<char> reserved;
    bool input_reserved;
  };

  void buildTrackGraph();
  void buildFixed();
  void buildNets();
  void buildAccess(NetData& net);
  void reserveAccess(NetData& net);

  // Off-track pin access (Sec. 3.7-3.8, 4.3-4.4); see OffTrack.md.
  class Legality;
  struct AccessStats
  {
    int64_t pins = 0;
    int64_t searches = 0;
    int64_t paths = 0;
    int64_t empty = 0;
    int64_t fallback = 0;
    int64_t drc_rejected = 0;
    int64_t same_net_rejected = 0;
    int64_t grid_vertices = 0;
    int64_t settled = 0;
    double search_seconds = 0;
    int64_t instances = 0;
    int64_t classes = 0;
    int64_t class_reused = 0;
    int64_t reselected = 0;
    int64_t mapping_failed = 0;
    int64_t selections = 0;
    int64_t selection_nodes = 0;
    int64_t selection_optimal = 0;
    int64_t primaries = 0;
    double seconds = 0;
    void add(const AccessStats& o);
  };
  void buildSameNetRules();
  // pre_skip: net skip flags before buildAccess.
  void buildOffTrackAccess(const std::vector<char>& pre_skip);
  // Rectangles of the pin on its access layer (the lowest routing layer
  // that is a track layer or directly below the lowest one).
  bool pinShapes(const Term& term,
                 frLayerNum& layer,
                 std::vector<odb::Rect>& rects) const;
  // Catalogue of the pin's off-track access paths against the current
  // routing space (Sec. 4.3).
  void offTrackCatalogue(int net,
                         const Term& term,
                         const Worker* worker,
                         std::vector<AccessOption>& out,
                         AccessStats& stats) const;
  void candidateOf(const AccessOption& option,
                   int net,
                   br::AccessCandidate& candidate) const;
  // Chooses primary paths for `pins` (net ids per pin) by destructive
  // bounding; -1 for pins without a feasible path.
  std::vector<int> choosePrimaries(
      const std::vector<const std::vector<AccessOption>*>& catalogues,
      const std::vector<int>& nets,
      AccessStats& stats) const;
  bool regenerateAccess(NetData& net, int t, Worker& worker);
  bool optionLegal(const AccessOption& option,
                   int net,
                   const Worker* worker) const;
  const frViaDef* viaDefAbove(frLayerNum lower) const;
  // Route nodes -> wires and vias of the net; via defs of `hints` at their
  // points take precedence over the defaults.
  void routeToGeometry(const std::vector<br::RouteNode>& nodes,
                       const std::vector<OutVia>& hints,
                       NetData& net) const;
  // Sec. 4.4 post-processing of one connection.
  void postprocessConnection(std::vector<br::RouteNode>& nodes,
                             const std::vector<OutVia>& hints,
                             int net,
                             Worker& worker,
                             std::vector<OutWire>& patches) const;
  // `local`: the worker whose region contains the net (Sec. 5.1).
  void setReservations(NetData& net, bool active, Worker* local = nullptr);
  void reserveInputWires();
  void setInputReservation(NetData& net, bool active, Worker* local = nullptr);
  // Diff-net conflict with the shapes of an R-tree (input wires only).
  bool treeConflicts(const RTree<int>& tree,
                     const RShape& shape,
                     int net) const;
  // Cost unit of the search: discount_denominator when discounts are used.
  br::Cost unit() const { return graph_.wire_unit; }
  bool modernSearch() const
  {
    return params_.reservations || params_.potential != 0 || params_.spreading;
  }
  // Tiles where global routing left room for spreading (Sec. 4.2).
  void buildSpreadTiles();
  void setReservation(Term& term,
                      int net,
                      bool active,
                      Worker* local = nullptr);
  void accessShapes(const AccessOption& option,
                    std::vector<RShape>& shapes) const;
  // extra_dilation enlarges the routing area of later rip-up passes.
  bool connect(NetData& net,
               bool ripup,
               Worker& worker,
               int extra_dilation = 0);
  bool connectWithoutOwnShapes(NetData& net,
                               bool ripup,
                               Worker& worker,
                               int extra_dilation);
  bool connectInRegion(NetData& net,
                       int dilation,
                       bool ripup,
                       const std::vector<std::vector<char>>& valid,
                       Worker& worker);
  void buildRegion(const NetData& net, int dilation, br::Region& region);
  void buildGuideTiles(NetData& net);
  // Region of the guide path from the routed component to the nearest
  // unconnected terminal; false if the guides do not connect them.
  bool buildPathRegion(const NetData& net,
                       const std::vector<char>& connected,
                       const std::vector<size_t>& tree,
                       const std::vector<std::vector<char>>& valid,
                       Worker& worker);
  // Tile box of everything the net's searches may touch, with halo.
  void netTileBox(const NetData& net,
                  int& xlo,
                  int& ylo,
                  int& xhi,
                  int& yhi) const;
  void routeParallel(const std::vector<int>& order, std::vector<int>& failed);
  void routeRipupQueue(const std::vector<int>& failed);
  bool ripupSequence(int net,
                     int depth,
                     int dilation,
                     std::vector<SavedRouting>& journal);
  SavedRouting saveRouting(const NetData& net) const;
  int64_t input_wires_ = 0;
  int64_t input_reserved_ = 0;
  void restoreRouting(const SavedRouting& saved);
  void pathToGeometry(const std::vector<size_t>& path, NetData& net) const;
  void shapesOf(const NetData& net, std::vector<RShape>& shapes) const;
  void viaShapes(const OutVia& via, std::vector<RShape>& shapes) const;
  RShape wireShape(const OutWire& wire) const;
  // Distance rules, shape classes and wire types from the drt tech (Sec.
  // 3.1, 3.2), the shape grid (Sec. 3.3), checker (3.4) and fast grid (3.6).
  void buildRoutingSpace();
  br::SpacingTable spacingTable(frLayer* layer) const;
  void rebuildFastGrid();
  // Shape-grid layer of a drt layer, or -1: track layers, their cuts, and
  // the pin layer and cut directly below the lowest track layer.
  int gridLayer(frLayerNum layer) const;
  // Thread-safe (per-thread class cache).
  br::ClassId shapeClass(const RShape& shape, int grid_layer) const;
  // Grid layer above a cut layer that holds its inter-layer projection, or
  // -1 (Sec. 3.2).
  int projectionLayer(int grid_layer) const;
  // Inserts / removes a shape in the shape grid (with the inter-layer
  // projection of a cut, Sec. 3.2) and updates the fast grid (Sec. 3.6).
  void changeShape(const RShape& shape,
                   br::ShapeType type,
                   int owner,
                   bool insert,
                   Worker* local);
  // Sec. 4.4: the shapes of the components being connected (pins,
  // reservations) are removed temporarily from routing space.
  void setOwnShapes(NetData& net, bool present, Worker* local);
  // Track-graph vertex kind usable by `net`: the fast grid, and where it
  // reports a blockage, the distance rule checker with the net's own
  // shapes treated as removed (Sec. 4.4).
  bool usable(size_t v, br::Kind kind, int net) const;
  void commitNet(NetData& net, Worker* local = nullptr);
  // Sec. 5.1: a sequence of balanced partitions into thread regions.
  // initial_regions: regions of the first partition, 0 = parameter;
  // stage names the step in the log.
  void routePartitions(const std::vector<int>& order,
                       std::vector<int>& failed,
                       int initial_regions = 0,
                       const char* stage = "partition");
  void partitionRegions(const std::vector<int>& nets,
                        const std::vector<std::array<int, 4>>& boxes,
                        int count,
                        bool shifted,
                        std::vector<std::array<int, 4>>& regions) const;
  void ripupNet(NetData& net);
  // Rip-up victims: nets reported by the distance rule checker as blocking
  // the net's new shapes (Sec. 3.4, 4.2).
  void collectVictims(const NetData& net, std::vector<int>& victims) const;
  // Distance rule check of a shape against everything but `net`'s own
  // shapes (the adapter for off-track paths and pin access).
  bool shapeConflicts(const RShape& shape,
                      int net,
                      const Worker* worker = nullptr) const;
  frCoord spacing(int z, frCoord width1, frCoord width2, frCoord prl) const;
  void writeNets();

  int tileX(frCoord x) const;
  int tileY(frCoord y) const;

  frDesign* design_;
  utl::Logger* logger_;
  RouterConfiguration* router_cfg_;
  BonnRouteDRParams params_;

  br::TrackGraph graph_;
  // Interning is internally locked (classes appear during routing).
  mutable br::ShapeClassTable shape_classes_;
  uint64_t class_generation_ = 0;
  // Wire type 0: default widths and vias (the fast-grid wire type); then
  // one wire type per non-default rule (Sec. 3.2).
  std::vector<br::WireType> wire_types_;
  br::ShapeGrid shape_grid_;
  br::DistanceRuleChecker checker_;
  br::FastGrid grid_;
  br::FastGridBuilder fast_builder_;
  br::ShapeGridWriteCache serial_cache_;
  // Shape-grid layers per drt layer (-1: not covered); number of track
  // grid layers (2 z + 1 for z track layers); the low layers follow.
  std::vector<int> grid_of_fr_;
  int track_grid_layers_ = 0;
  // During buildFixed the fast grid is rebuilt once at the end.
  bool bulk_insert_ = false;
  // Tiles between a net's routing area and its partition box border, so
  // that its routing-space updates and queries stay inside (Sec. 5.1).
  int halo_tiles_ = 1;
  size_t label_limit_ = 0;
  std::atomic<int64_t> inconsistent_updates_{0};
  std::vector<Worker> workers_;
  std::vector<LayerInfo> layers_;  // per z
  std::vector<int> z_of_fr_;       // routing layer -> z, -1
  std::vector<int> cut_z_of_fr_;   // cut layer -> lower z, -1
  frCoord tile_x0_ = 0;
  frCoord tile_y0_ = 0;
  frCoord tile_dx_ = 1;
  frCoord tile_dy_ = 1;
  int tiles_x_ = 1;
  odb::Rect die_;
  odb::dbBlock* block_ = nullptr;
  std::vector<std::vector<uint8_t>> spread_tiles_;
  int64_t spread_tile_count_ = 0;
  int tiles_y_ = 1;

  std::vector<NetData> nets_;
  std::unordered_map<const frNet*, int> net_index_;
  // Same-net rules and tau per track layer z; cleanup layers per drt layer.
  std::vector<br::SameNetLayerRules> same_net_rules_;
  std::vector<br::GeoCoord> tau_;
  std::vector<br::CleanupLayer> cleanup_layers_;
  AccessStats access_stats_;
  // Statistics.
  int64_t ripup_events_ = 0;
  int64_t rollbacks_ = 0;
  double partition_wall_ = 0;
  double partition_busy_ = 0;
  double partition_capacity_ = 0;
};

}  // namespace drt
