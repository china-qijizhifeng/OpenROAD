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
#include "frBaseTypes.h"
#include "frRTree.h"
#include "odb/geom.h"

namespace utl {
class Logger;
}

namespace drt {

class frDesign;
class frNet;
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
  bool verbose = true;

  bool set(const std::string& name, double value);
};

class BonnRouteDR
{
 public:
  BonnRouteDR(frDesign* design,
              utl::Logger* logger,
              RouterConfiguration* router_cfg,
              const BonnRouteDRParams& params);
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
  };
  struct Term
  {
    std::vector<AccessOption> options;
    std::vector<RShape> reservation;
    bool reservation_active = false;
    int chosen = -1;
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
    return params_.reservations || params_.potential != 0;
  }
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
  // initial_regions: regions of the first partition, 0 = parameter.
  void routePartitions(const std::vector<int>& order,
                       std::vector<int>& failed,
                       int initial_regions = 0);
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
  int tiles_y_ = 1;

  // Shape index (shape grid of Sec. 3.3, reduced to what the off-track
  // checks and rip-up need): fixed shapes and pins, and routed shapes.
  std::vector<RTree<int>> fixed_tree_;        // owner net or -1
  std::vector<RTree<int>> routed_tree_;       // owner net
  std::vector<RTree<int>> reservation_tree_;  // primary access paths
  std::vector<frCoord> max_spacing_;

  std::vector<NetData> nets_;
  std::unordered_map<const frNet*, int> net_index_;
  // Statistics.
  int64_t ripup_events_ = 0;
  int64_t rollbacks_ = 0;
};

}  // namespace drt
