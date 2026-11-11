// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Blockage grid for off-track wiring (Gester et al., TODAES 18(2), 2013,
// Sec. 3.8, Algorithm 3; Maßberg and Nieberg, Discrete Appl. Math. 161(12),
// 2013). Coordinates are the Hanan lines of the blockage borders and of
// s, t, plus lines at multiples of tau around groups of Hanan lines whose
// consecutive distances are below 4 tau.

#pragma once

#include <cstdint>
#include <vector>

#include "br/BrGeometry.h"

namespace drt::br {

// Algorithm 3 ("Blockage Grid Vertical"): from the blockage border
// coordinates x_j, x'_j (j = 1..n) and x_s, x_t (all passed in `hanan`, in
// any order, duplicates allowed) and tau > 0, the list V of coordinates:
//
//   1: V' := {x_s, x_t, -inf, inf}; V := {}
//   2: add x_j, x'_j to V'
//   3: sort V' and remove duplicates
//   4: for i = 1 .. |V'|:
//   5:   c_min := i - min{k in N | x_{i-k} - x_{i-k-1} >= 4 tau}
//   6:   c_max := i + min{k in N | x_{i+k+1} - x_{i+k} >= 4 tau}
//   7:   for lambda in Z with x_i + lambda tau in [x_cmin - 2 tau,
//                                                 x_cmax + 2 tau]:
//   8:     add x_i + lambda tau to V
//  11: sort V and remove duplicates
//
// (The paper's line 7 reads "x_i + lambda * u"; u is tau.) N contains 0, so
// [c_min, c_max] is the maximal run of Hanan lines around i whose
// consecutive gaps are below 4 tau. The sentinels -inf and inf only bound
// these runs; they generate no coordinates themselves.
std::vector<GeoCoord> blockageGridCoordinates(
    const std::vector<GeoCoord>& hanan,
    GeoCoord tau);

// One layer of a multi-layer blockage grid.
struct BlockageGridLayer
{
  GeoCoord tau = 1;
  std::vector<GeoCoord> xs;  // vertical lines, sorted and unique
  std::vector<GeoCoord> ys;  // horizontal lines, sorted and unique
  size_t numVertices() const { return xs.size() * ys.size(); }
};

// Input of a multi-layer blockage grid: per wiring layer its blockages and
// tau; per via layer (between wiring layers z and z + 1) its blockages.
struct BlockageGridInput
{
  struct Layer
  {
    GeoCoord tau = 1;
    std::vector<Box> blockages;
    // The caller uses this layer only at `points` (see buildBlockageGrid).
    bool points_only = false;
  };
  std::vector<Layer> layers;
  std::vector<std::vector<Box>> via_blockages;  // size layers.size() - 1
  // Coordinates every layer must contain: s and t, pin borders, on-track
  // target points, the search window.
  std::vector<GridPoint> points;
  std::vector<Box> borders;
};

// Builds the coordinates of each wiring layer (Sec. 3.8: "we have a
// blockage grid for each wiring and via layer. Coordinates are also induced
// by blockages of neighboring layers."):
//   H_z = borders of the blockages of wiring layers z - 1, z, z + 1 and of
//         the via layers z - 1 | z and z | z + 1, plus `points`/`borders`;
//   X_z = Algorithm 3 (H_z, tau_z), separately for x and y;
//   via layer z | z + 1 grid = X_z u X_{z+1};
//   final wiring layer grid = X_{z-1} u X_z u X_{z+1}
// so every via-layer vertex is a vertex of both wiring layers it connects.
// A `points_only` layer is one the caller connects only at `points`: a pin
// layer without wires whose vias are allowed only at the access points
// (sources). It keeps H_z = `points`/`borders` and its own coordinates, and
// neither its blockages nor the via blockages next to it induce coordinates
// elsewhere. (If vias were allowed elsewhere, a down-and-up via pair would
// reset the incoming direction and could shorten paths, so the reduction
// needs that restriction.)
// Adding lines keeps Theorem 3.2's guarantee (a superset of the grid still
// contains the endpoints of some shortest tau-feasible path).
std::vector<BlockageGridLayer> buildBlockageGrid(const BlockageGridInput& in);

}  // namespace drt::br
