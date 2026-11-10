// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// BonnRoute detailed routing on the drt database (Gester et al., TODAES
// 18(2), 2013, Sec. 3-4). Runs after pin access, guide processing and track
// assignment; writes complete net routings as initial routing, which
// FlexDR then cleans up (the paper's BR + external DRC clean-up flow).

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "br/BrCore.h"
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
  bool use_intervals = true;
  bool ripup = true;
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
  };
  struct Term
  {
    std::vector<AccessOption> options;
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
    int ripups = 0;
    bool routed = false;
    bool failed = false;
    bool skip = false;
    int64_t hpwl = 0;
  };

  void buildTrackGraph();
  void buildFixed();
  void buildNets();
  void buildAccess(NetData& net);
  void reserveAccess(NetData& net);
  bool routeNet(NetData& net, bool allow_ripup);
  bool connect(NetData& net, int dilation, bool ripup);
  void buildRegion(const NetData& net, int dilation);
  void pathToGeometry(const std::vector<size_t>& path, NetData& net) const;
  void shapesOf(const NetData& net, std::vector<RShape>& shapes) const;
  void viaShapes(const OutVia& via, std::vector<RShape>& shapes) const;
  RShape wireShape(const OutWire& wire) const;
  // Fast grid update for a shape; mode 0 fixed, 1 pin (owner), 2 routed.
  void markShape(const RShape& shape, int mode, int owner, int delta);
  void commitNet(NetData& net);
  void ripupNet(NetData& net);
  void collectVictims(const NetData& net, std::vector<int>& victims) const;
  bool shapeConflicts(const RShape& shape, int net) const;
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
  br::Region region_;
  std::unique_ptr<br::PathSearch> search_;
  std::vector<LayerInfo> layers_;  // per z
  std::vector<int> z_of_fr_;       // routing layer -> z, -1
  std::vector<int> cut_z_of_fr_;   // cut layer -> lower z, -1
  frCoord tile_x0_ = 0;
  frCoord tile_y0_ = 0;
  frCoord tile_dx_ = 1;
  frCoord tile_dy_ = 1;
  int tiles_x_ = 1;
  int tiles_y_ = 1;

  // Shape index (shape grid of Sec. 3.3, reduced to what the off-track
  // checks and rip-up need): fixed shapes and pins, and routed shapes.
  std::vector<RTree<int>> fixed_tree_;   // owner net or -1
  std::vector<RTree<int>> routed_tree_;  // owner net
  std::vector<frCoord> max_spacing_;

  std::vector<NetData> nets_;
  std::unordered_map<const frNet*, int> net_index_;
  // Nets ripped up by the last routeNet().
  std::vector<int> pending_victims_;

  // Statistics.
  int64_t searches_ = 0;
  int64_t heap_pops_ = 0;
  int64_t settled_ = 0;
  int64_t ripup_events_ = 0;
};

}  // namespace drt
