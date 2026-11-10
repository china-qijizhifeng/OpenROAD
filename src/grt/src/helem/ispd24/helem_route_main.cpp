// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// helem_route: runs the HeLEM global router on ISPD 2024 contest inputs.
//
//   helem_route -cap design.cap -net design.net -output design.out
//               [-threads N] [-lr N] [-lem N] [-batch N] [-expand R]
//               [-wire-weight normalized|score] [-la-order short|long]
//               [-la-batch N] [-cap2d-scale S] [-tau T] [-rho R]
//               [-exclude-own] [-verify] [-quiet]
//
// The last line on stdout is a JSON summary for experiment scripts.

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "core/GCellGraph.h"
#include "core/HelemRouter.h"
#include "ispd24/Ispd24Io.h"
#include "stt/flute.h"

namespace {


void usage()
{
  std::fprintf(stderr,
               "usage: helem_route -cap F -net F -output F [options]\n");
}

}  // namespace

int main(int argc, char** argv)
{
  using grt::helem::HelemOptions;
  std::string cap_path, net_path, output_path;
  HelemOptions options;
  int flute_accuracy = 10;
  bool flute_in_dbu = true;
  options.threads = 8;
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
      options.threads = std::stoi(value());
    } else if (arg == "-lr") {
      options.lr_iterations = std::stoi(value());
    } else if (arg == "-lem") {
      options.lem_iterations = std::stoi(value());
    } else if (arg == "-bcd") {
      options.bcd_passes = std::stoi(value());
    } else if (arg == "-batch") {
      options.batch_size = std::stoi(value());
    } else if (arg == "-expand") {
      options.bbox_expand_ratio = std::stod(value());
    } else if (arg == "-expand-max") {
      options.bbox_expand_max = std::stoi(value());
    } else if (arg == "-tau") {
      options.tau = std::stod(value());
    } else if (arg == "-rho") {
      options.rho_initial = std::stod(value());
    } else if (arg == "-sigma") {
      options.rho_growth = std::stod(value());
    } else if (arg == "-y0") {
      options.multiplier_initial = std::stod(value());
    } else if (arg == "-wire-weight") {
      const std::string mode = value();
      options.wire_weight_mode
          = mode == "score" ? grt::helem::WireWeightMode::kScoreLength
                            : grt::helem::WireWeightMode::kNormalizedLength;
    } else if (arg == "-la-order") {
      options.layer_assignment_order
          = value() == "long" ? grt::helem::LayerAssignmentOrder::kLongFirst
                              : grt::helem::LayerAssignmentOrder::kShortFirst;
    } else if (arg == "-la-batch") {
      options.layer_assignment_batch_size = std::stoi(value());
    } else if (arg == "-cap2d-scale") {
      options.capacity_2d_scale = std::stod(value());
    } else if (arg == "-flute-acc") {
      flute_accuracy = std::stoi(value());
    } else if (arg == "-flute-grid") {
      flute_in_dbu = false;
    } else if (arg == "-bend") {
      options.bend_cost_scale = std::stod(value());
    } else if (arg == "-pin-access") {
      options.pin_access_mode
          = value() != "center" ? grt::helem::PinAccessMode::kOtherPinsBox
                             : grt::helem::PinAccessMode::kNetCenter;
    } else if (arg == "-pin-high") {
      options.pin_access_prefer_high_layer = true;
    } else if (arg == "-exclude-own") {
      options.price_includes_own_route = false;
    } else if (arg == "-verify") {
      options.verify_ledger = true;
    } else if (arg == "-quiet") {
      options.verbose = false;
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
  grt::helem::GCellGraph graph;
  if (!grt::helem::readIspd24Resource(cap_path, graph, error)) {
    std::fprintf(stderr, "error: %s\n", error.c_str());
    return 1;
  }
  std::vector<grt::helem::RoutingNet> nets;
  if (!grt::helem::readIspd24Nets(net_path, nets, error)) {
    std::fprintf(stderr, "error: %s\n", error.c_str());
    return 1;
  }
  const double read_seconds = elapsed();
  std::printf("[helem] read layers=%d grid=%dx%d nets=%zu seconds=%.3f\n",
              graph.num_layers,
              graph.x_size,
              graph.y_size,
              nets.size(),
              read_seconds);

  stt::flt::Flute flute;
  flute.completeLUT();
  // GCell centres in DBU. The contest grid is anisotropic (e.g. 3800 x 2800),
  // so the Steiner tree is built on centre coordinates, not on indices.
  std::vector<int> center_x(graph.x_size, 0), center_y(graph.y_size, 0);
  for (int x = 0; x + 1 < graph.x_size; ++x) {
    center_x[x + 1] = center_x[x] + graph.horizontal_edge_length[x];
  }
  for (int y = 0; y + 1 < graph.y_size; ++y) {
    center_y[y + 1] = center_y[y] + graph.vertical_edge_length[y];
  }
  auto to_index = [](const std::vector<int>& centers, int coordinate) {
    return static_cast<int>(
        std::lower_bound(centers.begin(), centers.end(), coordinate)
        - centers.begin());
  };
  auto steiner = [&](const std::vector<grt::helem::GCellPoint>& points) {
    std::vector<int> xs, ys;
    xs.reserve(points.size());
    ys.reserve(points.size());
    for (const auto& p : points) {
      xs.push_back(flute_in_dbu ? center_x[p.x] : p.x);
      ys.push_back(flute_in_dbu ? center_y[p.y] : p.y);
    }
    stt::Tree tree = flute.flute(xs, ys, flute_accuracy);
    if (flute_in_dbu) {
      // Steiner points lie on the Hanan grid of the inputs, i.e. on centres.
      for (stt::Branch& branch : tree.branch) {
        branch.x = to_index(center_x, branch.x);
        branch.y = to_index(center_y, branch.y);
      }
    }
    std::vector<grt::helem::SteinerEdge> edges;
    const int branches = tree.branchCount();
    edges.reserve(branches);
    for (int i = 0; i < branches; ++i) {
      const stt::Branch& branch = tree.branch[i];
      const stt::Branch& neighbor = tree.branch[branch.n];
      edges.push_back({{branch.x, branch.y}, {neighbor.x, neighbor.y}});
    }
    return edges;
  };
  auto logger = [](const std::string& message) {
    std::printf("%s\n", message.c_str());
    std::fflush(stdout);
  };
  grt::helem::HelemRouter router(graph, nets, options, steiner, logger);
  const std::vector<grt::helem::NetRoute3D> routes = router.run();
  const double route_seconds = elapsed() - read_seconds;
  if (!grt::helem::writeIspd24Solution(output_path, nets, routes, error)) {
    std::fprintf(stderr, "error: %s\n", error.c_str());
    return 1;
  }
  const double total_seconds = elapsed();
  const grt::helem::HelemScore score = router.score();
  std::printf(
      "{\"read_s\":%.3f,\"route_s\":%.3f,\"total_s\":%.3f,"
      "\"wire\":%.4f,\"via\":%.4f,\"overflow\":%.4f,\"score\":%.4f,"
      "\"vias\":%lld,\"stages\":[",
      read_seconds,
      route_seconds,
      total_seconds,
      score.wire_cost,
      score.via_cost,
      score.overflow_cost,
      score.total(),
      static_cast<long long>(score.via_count));
  bool first = true;
  for (const auto& report : router.stageReports()) {
    std::printf("%s{\"stage\":\"%s\",\"s\":%.3f,\"of_total\":%lld,"
                "\"of_max\":%d,\"of_edges\":%lld}",
                first ? "" : ",",
                report.stage.c_str(),
                report.seconds,
                static_cast<long long>(report.total_overflow_2d),
                report.max_overflow_2d,
                static_cast<long long>(report.overflowed_edges_2d));
    first = false;
  }
  std::printf("]}\n");
  return 0;
}
