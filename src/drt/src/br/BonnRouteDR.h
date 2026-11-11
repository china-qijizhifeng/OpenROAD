// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// BonnRoute detailed routing on the drt database (Gester et al., TODAES
// 18(2), 2013, Sec. 3-4). Runs after pin access, guide processing and track
// assignment; writes complete net routings as initial routing, which
// FlexDR then cleans up (the paper's BR + external DRC clean-up flow).

#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "br/BrCore.h"
#include "br/BrPotentialOracle.h"
#include "br/BrGeometry.h"
#include "br/BrPathCleanup.h"
#include "br/BrPinAccess.h"
#include "br/BrSameNetRules.h"
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
  struct RShape
  {
    frLayerNum layer;
    odb::Rect rect;
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
    // Sec. 5.1 partitions: a thread commits only inside its region. The
    // shared shape index is updated after each partition; until then the
    // worker sees its own commits and releases through these.
    std::vector<RTree<int>> local_routed;
    std::unordered_set<int> released;
    std::vector<std::pair<RShape, int>> released_shapes;
    std::vector<int> committed;
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
  void buildLowLayerIndex();
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
  // With `local`, shape-index updates are deferred (see Worker).
  void setReservations(NetData& net, bool active, Worker* local = nullptr);
  void reserveInputWires();
  void setInputReservation(NetData& net, bool active, Worker* local = nullptr);
  bool treeConflicts(const RTree<int>& tree,
                     const RShape& shape,
                     int net,
                     const std::unordered_set<int>* skip = nullptr) const;
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
  // Fast grid update: 0 fixed, 1 pin, 2 routed, 3 owned reservation.
  void markShape(const RShape& shape, int mode, int owner, int delta);
  void commitNet(NetData& net, Worker* local = nullptr);
  // Applies a worker's deferred shape-index updates (serial).
  void mergeWorker(Worker& worker);
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
  void collectVictims(const NetData& net, std::vector<int>& victims) const;
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
  br::FastGrid grid_;
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

  // Shape index (shape grid of Sec. 3.3, reduced to what the off-track
  // checks and rip-up need): fixed shapes and pins, and routed shapes.
  std::vector<RTree<int>> fixed_tree_;        // owner net or -1
  std::vector<RTree<int>> routed_tree_;       // owner net
  std::vector<RTree<int>> reservation_tree_;  // primary access paths
  std::vector<frCoord> max_spacing_;

  std::vector<NetData> nets_;
  std::unordered_map<const frNet*, int> net_index_;
  // Same-net rules and tau per track layer z; cleanup layers per drt layer.
  std::vector<br::SameNetLayerRules> same_net_rules_;
  std::vector<br::GeoCoord> tau_;
  std::vector<br::CleanupLayer> cleanup_layers_;
  // Fixed shapes of the pin layer and cut directly below the lowest track
  // layer (not covered by the shape index), per drt layer.
  std::vector<RTree<int>> low_tree_;
  std::vector<frCoord> low_spacing_;
  AccessStats access_stats_;
  // Statistics.
  int64_t ripup_events_ = 0;
  int64_t rollbacks_ = 0;
};

}  // namespace drt
