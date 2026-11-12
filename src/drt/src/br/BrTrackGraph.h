// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Track graph of BonnRoute (Gester et al., TODAES 18(2), 2013, Sec. 3.5)
// and the shape kinds whose legality the fast grid answers per vertex
// (Sec. 3.6). No drt dependencies.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "br/BrGeometry.h"

namespace drt::br {

using Cost = int64_t;
constexpr Cost kInfCost = std::numeric_limits<Cost>::max() / 4;

// One wiring layer of the track graph (Sec. 3.5). Vertices are the
// intersections of the layer's tracks with the tracks of the adjacent
// layers, which have the perpendicular preferred direction.
struct TrackLayer
{
  bool horizontal = true;
  std::vector<Coord> tracks;     // cross coordinates, sorted and unique
  std::vector<Coord> positions;  // along coordinates, sorted and unique
  // Per position: track index of layer z +/- 1 at this coordinate, or -1.
  std::vector<int> up_track;
  std::vector<int> down_track;
  // Per track: position index of this track's coordinate in layer z +/- 1.
  std::vector<int> up_pos;
  std::vector<int> down_pos;
  // Per track / position: global routing tile row / column (Region lookup).
  std::vector<int> track_tile;
  std::vector<int> pos_tile;
  bool jogs = true;
  Cost jog_factor = 4;  // beta_z of Sec. 4.1
  // Jog cost in tiles of SearchRequest::jog_low_tiles (<= jog_factor).
  Cost jog_factor_low = 4;
  size_t base = 0;

  int numTracks() const { return static_cast<int>(tracks.size()); }
  int numPositions() const { return static_cast<int>(positions.size()); }
  // (x, y) of vertex (i, j).
  std::pair<Coord, Coord> point(int i, int j) const
  {
    return horizontal ? std::make_pair(positions[j], tracks[i])
                      : std::make_pair(tracks[i], positions[j]);
  }
};

class TrackGraph
{
 public:
  std::vector<TrackLayer> layers;
  std::vector<Cost> via_cost;  // via between z and z + 1
  // Cost of one DBU of preferred-direction wire. Values above 1 let a
  // reservation discount (arXiv:2111.06169, Sec. 5.2) stay integral; jog
  // factors, via costs and penalties are expressed in the same unit.
  Cost wire_unit = 1;
  int tiles_x = 1;
  int tiles_y = 1;

  // Fills the up/down maps and vertex bases; tile_x / tile_y map a
  // coordinate to its tile index.
  template <typename TileX, typename TileY>
  void finalize(TileX tile_x, TileY tile_y);

  int numLayers() const { return static_cast<int>(layers.size()); }
  size_t numVertices() const { return num_vertices_; }
  size_t vertex(int z, int i, int j) const
  {
    const TrackLayer& l = layers[z];
    return l.base + static_cast<size_t>(i) * l.positions.size() + j;
  }
  // Inverse of vertex().
  void decode(size_t v, int& z, int& i, int& j) const
  {
    z = static_cast<int>(
            std::upper_bound(layer_bases_.begin(), layer_bases_.end(), v)
            - layer_bases_.begin())
        - 1;
    const TrackLayer& l = layers[z];
    const size_t rem = v - l.base;
    i = static_cast<int>(rem / l.positions.size());
    j = static_cast<int>(rem % l.positions.size());
  }
  std::pair<Coord, Coord> point(size_t v) const
  {
    int z, i, j;
    decode(v, z, i, j);
    return layers[z].point(i, j);
  }
  Cost minViaCost() const
  {
    Cost c = kInfCost;
    for (const Cost v : via_cost) {
      c = std::min(c, v);
    }
    return c == kInfCost ? 0 : c;
  }

 private:
  size_t num_vertices_ = 0;
  std::vector<size_t> layer_bases_;
};

template <typename TileX, typename TileY>
void TrackGraph::finalize(TileX tile_x, TileY tile_y)
{
  auto index_of = [](const std::vector<Coord>& v, Coord c) {
    auto it = std::lower_bound(v.begin(), v.end(), c);
    return (it != v.end() && *it == c) ? static_cast<int>(it - v.begin()) : -1;
  };
  num_vertices_ = 0;
  layer_bases_.clear();
  for (int z = 0; z < numLayers(); z++) {
    TrackLayer& l = layers[z];
    l.base = num_vertices_;
    layer_bases_.push_back(num_vertices_);
    num_vertices_ += l.tracks.size() * l.positions.size();
    l.up_track.assign(l.positions.size(), -1);
    l.down_track.assign(l.positions.size(), -1);
    l.up_pos.assign(l.tracks.size(), -1);
    l.down_pos.assign(l.tracks.size(), -1);
    for (int dz : {-1, 1}) {
      const int nz = z + dz;
      if (nz < 0 || nz >= numLayers()) {
        continue;
      }
      const TrackLayer& n = layers[nz];
      for (int j = 0; j < l.numPositions(); j++) {
        (dz > 0 ? l.up_track : l.down_track)[j]
            = index_of(n.tracks, l.positions[j]);
      }
      for (int i = 0; i < l.numTracks(); i++) {
        (dz > 0 ? l.up_pos : l.down_pos)[i]
            = index_of(n.positions, l.tracks[i]);
      }
    }
    l.track_tile.resize(l.tracks.size());
    l.pos_tile.resize(l.positions.size());
    for (int i = 0; i < l.numTracks(); i++) {
      l.track_tile[i]
          = l.horizontal ? tile_y(l.tracks[i]) : tile_x(l.tracks[i]);
    }
    for (int j = 0; j < l.numPositions(); j++) {
      l.pos_tile[j]
          = l.horizontal ? tile_x(l.positions[j]) : tile_y(l.positions[j]);
    }
  }
}

// Shape kinds whose legality the on-track path search asks for (Sec. 3.6,
// 4.1). Vertex kinds are the zero-dimensional stick figures of the fast
// grid: preferred-direction wire, jog (non-preferred wire), bottom pad of
// the via going up, the via layer above (cut and inter-layer projection)
// and top pad of the via coming from below.
enum Kind : uint8_t
{
  kWire = 0,
  kJog = 1,
  kViaUpPad = 2,
  kCutUp = 3,
  kViaDownPad = 4,
  // Edge kinds, addressed by the lower endpoint: the preferred-direction
  // edge to j + 1 and the jog to track i + 1. An edge's level is the
  // maximum of its endpoints' wire (jog) levels and of the level of the
  // swept wire on the edge; the latter can exceed both endpoints (a shape
  // between them, Fig. 4 zigzag edge).
  kWireEdge = 5,
  kJogEdge = 6,
  kNumKinds = 7
};

// Calls fn(i, j) for every vertex of layer `l` inside the box [xlo, xhi] x
// [ylo, yhi] (vertex kinds), or for every edge whose segment meets it:
// kWireEdge (j, j + 1) on track i, kJogEdge (i, i + 1) at position j.
template <typename F>
void forEachBlocked(const TrackLayer& l,
                    const Kind kind,
                    const Coord xlo,
                    const Coord ylo,
                    const Coord xhi,
                    const Coord yhi,
                    F&& fn)
{
  if (xlo > xhi || ylo > yhi) {
    return;
  }
  const Coord clo = l.horizontal ? ylo : xlo;
  const Coord chi = l.horizontal ? yhi : xhi;
  const Coord alo = l.horizontal ? xlo : ylo;
  const Coord ahi = l.horizontal ? xhi : yhi;
  auto first_ge = [](const std::vector<Coord>& v, Coord c) {
    return static_cast<int>(std::lower_bound(v.begin(), v.end(), c)
                            - v.begin());
  };
  auto last_le = [](const std::vector<Coord>& v, Coord c) {
    return static_cast<int>(std::upper_bound(v.begin(), v.end(), c) - v.begin())
           - 1;
  };
  int i0 = first_ge(l.tracks, clo);
  int i1 = last_le(l.tracks, chi);
  int j0 = first_ge(l.positions, alo);
  int j1 = last_le(l.positions, ahi);
  if (kind == kWireEdge) {
    // Edges (j, j + 1) with positions[j] <= ahi and positions[j + 1] >= alo.
    j0 = std::max(0, j0 - 1);
    j1 = std::min(j1, l.numPositions() - 2);
  } else if (kind == kJogEdge) {
    i0 = std::max(0, i0 - 1);
    i1 = std::min(i1, l.numTracks() - 2);
  }
  for (int i = i0; i <= i1; i++) {
    for (int j = j0; j <= j1; j++) {
      fn(i, j);
    }
  }
}

}  // namespace drt::br
