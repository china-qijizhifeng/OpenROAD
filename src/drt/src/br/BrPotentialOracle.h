// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Future cost of the vertex search from the potentials of Ahrens et al.,
// arXiv:2111.06169 (BrPotential.h). The track graph is a subgraph of the
// grid graph G of Section 1.1 with the same per-DBU rates: preferred wires
// cost wire_unit (discount_rate on reserved runs), jogs jog_factor, vias
// via_cost. Distances in G are therefore feasible potentials for the
// track graph search. With the general model the routing area of the search
// (Region) gives infinite costs outside, as in the blockage-aware future
// cost of Peyer et al. that Gester et al. 2013, Sec. 4.1, use on detours.

#pragma once

#include <cstdint>
#include <vector>

#include "br/BrCore.h"
#include "br/BrPotential.h"

namespace drt::br {

class PotentialOracle : public FutureCostOracle
{
 public:
  enum class Model : uint8_t
  {
    kL1,       // Table 3 baseline
    kSimple,   // Theorem 3
    kGeneral,  // Section 4 with the search's Region and discounts
  };

  void prepare(Model model,
               const TrackGraph& graph,
               const SearchRequest& request,
               const GCellGrid& grid,
               int64_t discount_numerator,
               int64_t discount_denominator)
  {
    model_ = model;
    const int num_layers = graph.numLayers();
    rates_.horizontal.assign(num_layers, kInfCost);
    rates_.vertical.assign(num_layers, kInfCost);
    rates_.via_up.assign(std::max(0, num_layers - 1), kInfCost);
    // Discounts apply to preferred-direction runs only. The simple models
    // have one rate per layer and direction, so they take the discounted
    // rate where the request has discounted runs (arXiv:2111.06169,
    // Sec. 5.2: the l1 bound has to be multiplied by delta).
    std::vector<char> discounted(num_layers, 0);
    for (const DiscountRun& run : request.discounts) {
      discounted[run.z] = 1;
    }
    for (int z = 0; z < num_layers; z++) {
      const TrackLayer& layer = graph.layers[z];
      Cost preferred = graph.wire_unit;
      if (model != Model::kGeneral && discounted[z]) {
        preferred = std::min(preferred, request.discount_rate);
      }
      const Cost jog = layer.jogs ? layer.jog_factor : kInfCost;
      rates_.horizontal[z] = layer.horizontal ? preferred : jog;
      rates_.vertical[z] = layer.horizontal ? jog : preferred;
      if (z + 1 < num_layers) {
        rates_.via_up[z] = graph.via_cost[z];
      }
    }
    targets_.clear();
    for (const SearchRequest::Target& target : request.targets) {
      int z, i, j;
      graph.decode(target.vertex, z, i, j);
      const auto [x, y] = graph.layers[z].point(i, j);
      targets_.push_back({x, y, x, y, z, target.cost});
    }
    switch (model) {
      case Model::kL1:
        l1_.prepare(rates_, targets_);
        break;
      case Model::kSimple:
        simple_.prepare(rates_, targets_);
        break;
      case Model::kGeneral:
        wires_.clear();
        for (const DiscountRun& run : request.discounts) {
          const TrackLayer& layer = graph.layers[run.z];
          const auto [x1, y1] = layer.point(run.i, run.jlo);
          const auto [x2, y2] = layer.point(run.i, run.jhi);
          wires_.push_back({run.z,
                            x1,
                            y1,
                            x2,
                            y2,
                            discount_numerator,
                            discount_denominator});
        }
        builder_.build(*request.region, grid, rates_, wires_, tiles_);
        general_.prepare(tiles_, targets_);
        break;
    }
    prepares_[static_cast<int>(model)]++;
  }

  Cost query(const int z, const Coord x, const Coord y) const override
  {
    switch (model_) {
      case Model::kL1:
        return l1_.query(z, x, y);
      case Model::kSimple:
        return simple_.query(z, x, y);
      case Model::kGeneral:
        return general_.query(z, x, y);
    }
    return 0;
  }

  int64_t prepares(Model model) const
  {
    return prepares_[static_cast<int>(model)];
  }

 private:
  Model model_ = Model::kL1;
  LayerCostRates rates_;
  std::vector<TargetRectangle> targets_;
  std::vector<DiscountedWire> wires_;
  L1Potential l1_;
  SimplePotential simple_;
  TileCostModelBuilder builder_;
  TileCostModel tiles_;
  GeneralPotential general_;
  int64_t prepares_[3] = {0, 0, 0};
};

}  // namespace drt::br
