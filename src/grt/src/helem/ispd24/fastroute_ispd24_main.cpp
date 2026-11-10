// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// fastroute_ispd24: runs OpenROAD's FastRoute core on ISPD 2024 contest
// inputs, so that it routes the same GCell grid, capacities and pin access
// points as helem_route. Capacities are floored to whole tracks (FastRoute
// stores integer capacities); pin access uses grt::helem::choosePinAccessPoints.
//
//   fastroute_ispd24 -cap F -net F -output F [-threads N]
//                    [-congestion-iterations N] [-pin-access box|center]

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "FastRoute.h"
#include "core/GCellGraph.h"
#include "db_sta/dbSta.hh"
#include "core/HelemRouter.h"
#include "grt/GRoute.h"
#include "ispd24/Ispd24Io.h"
#include "odb/db.h"
#include "odb/dbTypes.h"
#include "sta/Sta.hh"
#include "stt/SteinerTreeBuilder.h"
#include "tcl.h"
#include "utl/Logger.h"
#include "utl/ServiceRegistry.h"

namespace {

// FastRoute reports positions as tile_size * (gcell + 0.5) + origin; with a
// tile size of 2 and origin 0 the gcell index is (position - 1) / 2.
constexpr int kTileSize = 2;

int toGCell(int position)
{
  return (position - 1) / kTileSize;
}

void usage()
{
  std::fprintf(stderr,
               "usage: fastroute_ispd24 -cap F -net F -output F [-threads N] "
               "[-congestion-iterations N] [-pin-access box|center]\n");
}

}  // namespace

int main(int argc, char** argv)
{
  namespace helem = grt::helem;
  std::string cap_path, net_path, output_path;
  int threads = 8;
  int congestion_iterations = 50;
  helem::PinAccessMode pin_access = helem::PinAccessMode::kOtherPinsBox;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto value = [&]() -> std::string {
      if (i + 1 >= argc) {
        usage();
        std::exit(2);
      }
      return argv[++i];
    };
    if (arg == "-cap") {
      cap_path = value();
    } else if (arg == "-net") {
      net_path = value();
    } else if (arg == "-output" || arg == "-out") {
      output_path = value();
    } else if (arg == "-threads") {
      threads = std::stoi(value());
    } else if (arg == "-congestion-iterations") {
      congestion_iterations = std::stoi(value());
    } else if (arg == "-pin-access") {
      pin_access = value() == "center" ? helem::PinAccessMode::kNetCenter
                                       : helem::PinAccessMode::kOtherPinsBox;
    } else {
      usage();
      return 2;
    }
  }
  if (cap_path.empty() || net_path.empty() || output_path.empty()) {
    usage();
    return 2;
  }

  const auto start = std::chrono::steady_clock::now();
  auto elapsed = [&]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now()
                                         - start)
        .count();
  };
  std::string error;
  helem::GCellGraph graph;
  std::vector<helem::RoutingNet> nets;
  if (!helem::readIspd24Resource(cap_path, graph, error)
      || !helem::readIspd24Nets(net_path, nets, error)) {
    std::fprintf(stderr, "error: %s\n", error.c_str());
    return 1;
  }
  const double read_seconds = elapsed();

  utl::Logger logger(nullptr);
  odb::dbDatabase* db = odb::dbDatabase::create();
  db->setLogger(&logger);
  odb::dbTech* tech = odb::dbTech::create(db, "ispd24");
  // FastRoute reads routing layers (width, resistance, upper via layer) in
  // run(); mirror the contest layers as routing layers with cut layers.
  for (int l = 0; l < graph.num_layers; ++l) {
    if (l > 0) {
      odb::dbTechLayer::create(
          tech, ("via" + std::to_string(l)).c_str(), odb::dbTechLayerType::CUT);
    }
    odb::dbTechLayer* layer = odb::dbTechLayer::create(
        tech, graph.layers[l].name.c_str(), odb::dbTechLayerType::ROUTING);
    layer->setWidth(100);
    layer->setDirection(graph.layers[l].horizontal
                            ? odb::dbTechLayerDir::HORIZONTAL
                            : odb::dbTechLayerDir::VERTICAL);
  }
  odb::dbChip* chip = odb::dbChip::create(db, tech);
  odb::dbBlock* block = odb::dbBlock::create(chip, "top");
  stt::SteinerTreeBuilder stt_builder(&logger);
  // FastRoute queries the timer for slacks (updateSlacks); without a liberty
  // library it skips timing, as in GlobalRouter runs without timing. The
  // setup mirrors the OpenROAD test fixture (tst::Fixture).
  sta::initSta();
  Tcl_Interp* interp = Tcl_CreateInterp();
  auto sta = std::make_unique<sta::dbSta>(interp, db, &logger);
  utl::ServiceRegistry service_registry(&logger);

  auto fastroute = std::make_unique<grt::FastRouteCore>(
      db, &logger, &service_registry, &stt_builder, sta.get());
  const int layers = graph.num_layers;
  fastroute->clear();
  fastroute->setVerbose(true);
  fastroute->setOverflowIterations(congestion_iterations);
  fastroute->setCriticalNetsPercentage(0);
  fastroute->setResistanceAware(false);
  fastroute->setNumThreads(threads);
  fastroute->setRegularX(true);
  fastroute->setRegularY(true);
  fastroute->setLowerLeft(0, 0);
  fastroute->setTileSize(kTileSize);
  fastroute->setGridsAndLayers(graph.x_size, graph.y_size, layers);
  fastroute->setGridMax(graph.x_size * kTileSize, graph.y_size * kTileSize);
  for (int l = 0; l < layers; ++l) {
    fastroute->addLayerDirection(l,
                                 graph.layers[l].horizontal
                                     ? odb::dbTechLayerDir::HORIZONTAL
                                     : odb::dbTechLayerDir::VERTICAL);
  }

  // Capacities, following GlobalRouter's contract: every edge of a layer
  // first gets the layer's nominal capacity (here the largest contest
  // capacity on that layer, floored to whole tracks), which also becomes the
  // layer's base capacity; edges with less capacity are then reduced with
  // addAdjustment, as GlobalRouter does for blockages. FastRoute layer
  // numbers are 1-based.
  auto edge_cap = [&](int l, int x, int y) {
    return graph.isRoutingLayer(l)
               ? static_cast<int>(std::floor(
                     graph.capacity[l][graph.layerEdgeIndex(l, x, y)]))
               : 0;
  };
  std::vector<int> nominal(layers, 0);
  for (int l = 0; l < layers; ++l) {
    for (int y = 0; y < graph.y_size; ++y) {
      for (int x = 0; x < graph.x_size; ++x) {
        const bool edge = graph.layers[l].horizontal ? x + 1 < graph.x_size
                                                     : y + 1 < graph.y_size;
        if (edge) {
          nominal[l] = std::max(nominal[l], edge_cap(l, x, y));
        }
      }
    }
  }
  fastroute->initEdges();
  for (int l = 0; l < layers; ++l) {
    if (graph.layers[l].horizontal) {
      for (int y = 0; y < graph.y_size; ++y) {
        for (int x = 0; x + 1 < graph.x_size; ++x) {
          fastroute->setEdgeCapacity(x, y, x + 1, y, l + 1, nominal[l]);
        }
      }
      fastroute->addHCapacity(nominal[l], l + 1);
    } else {
      for (int x = 0; x < graph.x_size; ++x) {
        for (int y = 0; y + 1 < graph.y_size; ++y) {
          fastroute->setEdgeCapacity(x, y, x, y + 1, l + 1, nominal[l]);
        }
      }
      fastroute->addVCapacity(nominal[l], l + 1);
    }
  }
  fastroute->initLowerBoundCapacities();
  std::vector<int> track_space(layers, 1);
  fastroute->initBlockedIntervals(track_space);
  fastroute->saveResourcesBeforeAdjustments();
  int64_t adjusted = 0;
  for (int l = 0; l < layers; ++l) {
    const bool horizontal = graph.layers[l].horizontal;
    for (int y = 0; y < graph.y_size; ++y) {
      for (int x = 0; x < graph.x_size; ++x) {
        if (horizontal ? x + 1 >= graph.x_size : y + 1 >= graph.y_size) {
          continue;
        }
        const int cap = edge_cap(l, x, y);
        if (cap < nominal[l]) {
          fastroute->addAdjustment(x,
                                   y,
                                   horizontal ? x + 1 : x,
                                   horizontal ? y : y + 1,
                                   l + 1,
                                   static_cast<uint16_t>(cap),
                                   true);
          ++adjusted;
        }
      }
    }
  }
  fastroute->initEdgesCapacityPerLayer();
  std::printf("[fastroute] nominal capacity per layer:");
  for (int l = 0; l < layers; ++l) {
    std::printf(" %d", nominal[l]);
  }
  std::printf("; adjusted edges=%lld\n", static_cast<long long>(adjusted));

  // Nets.
  std::vector<odb::dbNet*> db_nets(nets.size(), nullptr);
  // Pin layer interval per net for nets whose pins share one GCell.
  std::vector<std::tuple<int, int, int, int>> local_vias;
  int max_degree = 0;
  int routed_nets = 0;
  for (size_t n = 0; n < nets.size(); ++n) {
    const helem::RoutingNet& net = nets[n];
    odb::dbNet* db_net = odb::dbNet::create(block, net.name.c_str());
    db_nets[n] = db_net;
    const std::vector<const helem::AccessPoint*> choice
        = helem::choosePinAccessPoints(net, pin_access, false);
    std::vector<std::tuple<int, int, int>> pins;
    for (const helem::AccessPoint* ap : choice) {
      if (ap != nullptr) {
        pins.emplace_back(ap->x, ap->y, ap->layer);
      }
    }
    std::sort(pins.begin(), pins.end());
    pins.erase(std::unique(pins.begin(), pins.end()), pins.end());
    if (pins.size() < 2) {
      continue;
    }
    bool local = true;
    for (const auto& pin : pins) {
      local = local && std::get<0>(pin) == std::get<0>(pins.front())
              && std::get<1>(pin) == std::get<1>(pins.front());
    }
    if (local) {
      int low = layers, high = -1;
      for (const auto& pin : pins) {
        low = std::min(low, std::get<2>(pin));
        high = std::max(high, std::get<2>(pin));
      }
      if (low < high) {
        local_vias.emplace_back(static_cast<int>(n), low, high, 0);
      }
      continue;
    }
    grt::FrNet* fr_net = fastroute->addNet(db_net,
                                           /*is_clock=*/false,
                                           /*is_local=*/false,
                                           /*driver_idx=*/0,
                                           /*cost=*/1,
                                           /*min_layer=*/graph.min_routing_layer,
                                           /*max_layer=*/graph.topRoutingLayer(),
                                           /*slack=*/0.0f,
                                           /*edge_cost_per_layer=*/nullptr);
    for (const auto& [x, y, layer] : pins) {
      fr_net->addPin(x, y, layer);
    }
    max_degree = std::max(max_degree, static_cast<int>(pins.size()));
    ++routed_nets;
  }
  fastroute->setMaxNetDegree(max_degree);
  fastroute->initAuxVar();
  fastroute->setAllowCongestion(true);
  const double setup_seconds = elapsed() - read_seconds;
  std::printf("[fastroute] nets=%zu routed=%d local_vias=%zu setup_s=%.3f\n",
              nets.size(),
              routed_nets,
              local_vias.size(),
              setup_seconds);
  std::fflush(stdout);

  const double route_start = elapsed();
  grt::NetRouteMap routes = fastroute->run();
  const double route_seconds = elapsed() - route_start;

  // Convert to contest segments.
  std::vector<helem::NetRoute3D> out(nets.size());
  // odb specializes std::less for its pointer types; hash by address.
  std::unordered_map<odb::dbNet*, int> net_index;
  for (size_t n = 0; n < nets.size(); ++n) {
    net_index[db_nets[n]] = static_cast<int>(n);
  }
  for (const auto& [db_net, route] : routes) {
    helem::NetRoute3D& net_route = out[net_index.at(db_net)];
    for (const grt::GSegment& seg : route) {
      const int x0 = toGCell(seg.init_x), y0 = toGCell(seg.init_y);
      const int x1 = toGCell(seg.final_x), y1 = toGCell(seg.final_y);
      // FastRoute layers in GSegment are 1-based.
      const int l0 = seg.init_layer - 1, l1 = seg.final_layer - 1;
      if (l0 != l1) {
        net_route.segments.push_back(
            {x0, y0, std::min(l0, l1), x0, y0, std::max(l0, l1)});
      } else if (x0 != x1 || y0 != y1) {
        net_route.segments.push_back({std::min(x0, x1),
                                      std::min(y0, y1),
                                      l0,
                                      std::max(x0, x1),
                                      std::max(y0, y1),
                                      l0});
      }
    }
  }
  for (const auto& [n, low, high, unused] : local_vias) {
    const helem::AccessPoint* ap = nullptr;
    for (const helem::AccessPoint* choice :
         helem::choosePinAccessPoints(nets[n], pin_access, false)) {
      ap = ap != nullptr ? ap : choice;
    }
    out[n].segments.push_back({ap->x, ap->y, low, ap->x, ap->y, high});
  }
  if (!helem::writeIspd24Solution(output_path, nets, out, error)) {
    std::fprintf(stderr, "error: %s\n", error.c_str());
    return 1;
  }
  std::printf(
      "{\"read_s\":%.3f,\"setup_s\":%.3f,\"route_s\":%.3f,\"total_s\":%.3f,"
      "\"fastroute_overflow\":%d}\n",
      read_seconds,
      setup_seconds,
      route_seconds,
      elapsed(),
      fastroute->totalOverflow());
  return 0;
}
