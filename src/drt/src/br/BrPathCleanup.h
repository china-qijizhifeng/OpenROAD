// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Same-net post-processing of a routed connection (Gester et al., TODAES
// 18(2), 2013, Sec. 4.4): "we combine it with the corresponding off-track
// path if present, and then perform a postprocessing step on the resulting
// path to clean up same-net rule violations which occur particularly where
// on-track and off-track paths meet, but also within on-track paths".

#pragma once

#include <cstdint>
#include <functional>
#include <vector>

#include "br/BrGeometry.h"
#include "br/BrSameNetRules.h"

namespace drt::br {

struct CleanupLayer
{
  GeoCoord tau = 1;
  GeoCoord half_width = 0;
  SameNetLayerRules rules;
  bool horizontal = true;  // preferred direction (min-area extensions)
  // Manufacturing grid: extension lengths are multiples of it.
  GeoCoord grid = 1;
  // Wires only on tracks in the preferred direction (unidirectional
  // layers): no jog shifts or U collapses (they move wires off track).
  bool on_tracks = false;
};

struct CleanupRequest
{
  // The connection: consecutive nodes on one layer are wire segments,
  // consecutive nodes at one point on adjacent layers are vias. The first
  // and last nodes (pin / tree attachments) and all via nodes stay fixed.
  std::vector<RouteNode> nodes;
  std::vector<CleanupLayer> layers;  // indexed by RouteNode::layer
  // True if new metal of the net covering `box` on `layer` is legal
  // (diff-net rules; DiffNetLegality::allowed).
  std::function<bool(int layer, const Box& box)> allowed;
  // Pads on `layer` of the via at `p` going up (true) or coming from below
  // (false), for minimum area.
  std::function<void(int layer, const GridPoint& p, bool up, std::vector<Box>&)>
      via_pads;
  // False for points where the post-processing must not create a node
  // (wire end, bend or extension tip), e.g. a pin of the net away from its
  // access points: FlexDR treats pin-touching ends as pin connections.
  // Existing nodes stay where they are.
  std::function<bool(int layer, const GridPoint& p)> node_allowed;
};

struct CleanupPatch
{
  int layer = 0;
  GridPoint begin;
  GridPoint end;
};

struct CleanupStats
{
  int64_t runs = 0;
  int64_t violations_before = 0;  // short segments and notches
  int64_t violations_after = 0;
  int64_t moved_segments = 0;     // jog shifts (fix A / B)
  int64_t collapsed_notches = 0;  // U-leg collapses (fix C)
  int64_t area_violations = 0;    // isolated via-to-via pieces below min area
  int64_t area_fixed = 0;
  void add(const CleanupStats& o)
  {
    runs += o.runs;
    violations_before += o.violations_before;
    violations_after += o.violations_after;
    moved_segments += o.moved_segments;
    collapsed_notches += o.collapsed_notches;
    area_violations += o.area_violations;
    area_fixed += o.area_fixed;
  }
};

struct CleanupResult
{
  std::vector<RouteNode> nodes;
  // Min-area extensions (stubs beyond a via), emitted as extra wires.
  std::vector<CleanupPatch> patches;
  CleanupStats stats;
};

// Per layer run between vias / ends:
//  1. drop repeated points and interior collinear points (spikes included);
//  2. remove short segments s = (P_k, P_k+1) by shifting the perpendicular
//     neighbour a = (P_k-1, P_k) by s (fix A: s merges into the segment
//     before a) or b = (P_k+1, P_k+2) by -s (fix B), and collapse the
//     shorter leg of a U whose base is short (fix C); a change is kept only
//     if its new metal is `allowed` and it reduces the run's inner
//     violations without adding end violations;
//  3. extend via-to-via pieces below the minimum area beyond a via, in the
//     preferred direction and by a multiple of the manufacturing grid.
// Chain ends and vias never move, and no change creates a node where
// `node_allowed` is false.
CleanupResult cleanupConnection(const CleanupRequest& request);

}  // namespace drt::br
