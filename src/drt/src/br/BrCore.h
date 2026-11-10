// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// BonnRoute detailed routing core (Gester et al., TODAES 18(2), 2013,
// Sec. 3.5, 3.6, 4.1, 4.2): track graph, fast grid and the on-track
// interval path search with future cost. Independent of the drt database so
// it can be tested in isolation; BonnRouteDR adapts it to frDesign.

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace drt::br {

using Coord = int;
using Cost = int64_t;
constexpr Cost kInfCost = std::numeric_limits<Cost>::max() / 4;

// One wiring layer of the track graph (Sec. 3.5). Vertices are the
// intersections of the layer's tracks with the tracks of the adjacent
// layers, which have the perpendicular preferred direction.
struct TrackLayer
{
  bool horizontal = true;
  std::vector<Coord> tracks;     // cross coordinates, sorted
  std::vector<Coord> positions;  // along coordinates, sorted
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

// Shape kinds whose legality the fast grid stores per vertex (Sec. 3.6):
// preferred-direction wire, jog, bottom pad and cut of the via going up,
// top pad of the via coming from below.
enum Kind : uint8_t
{
  kWire = 0,
  kJog = 1,
  kViaUpPad = 2,
  kCutUp = 3,
  kViaDownPad = 4,
  kNumKinds = 5
};

class FastGrid
{
 public:
  static constexpr int32_t kNoOwner = -1;
  static constexpr int32_t kManyOwners = -2;

  void init(size_t num_vertices)
  {
    fixed_.assign(num_vertices, 0);
    pin_mask_.assign(num_vertices, 0);
    pin_owner_.assign(num_vertices, kNoOwner);
    count_.assign(num_vertices, {});
    history_.assign(num_vertices, 0);
  }
  void addFixed(size_t v, Kind k) { fixed_[v] |= (1u << k); }
  // A pin shape blocks every net except its owner.
  void addPin(size_t v, Kind k, int32_t owner)
  {
    if (owner < 0) {
      addFixed(v, k);
      return;
    }
    pin_mask_[v] |= (1u << k);
    if (pin_owner_[v] == kNoOwner) {
      pin_owner_[v] = owner;
    } else if (pin_owner_[v] != owner) {
      pin_owner_[v] = kManyOwners;
    }
  }
  // Routed (removable) shapes; saturating counter.
  void addRouted(size_t v, Kind k, int delta)
  {
    int c = count_[v][k] + delta;
    count_[v][k] = static_cast<uint8_t>(std::clamp(c, 0, 255));
  }
  bool hardBlocked(size_t v, Kind k, int32_t net) const
  {
    if ((fixed_[v] >> k) & 1u) {
      return true;
    }
    return ((pin_mask_[v] >> k) & 1u) && pin_owner_[v] != net;
  }
  bool usable(size_t v, Kind k, int32_t net) const
  {
    return !hardBlocked(v, k, net) && count_[v][k] == 0;
  }
  // Blocked only by routed shapes, which a rip-up can remove (Sec. 4.2).
  bool ripupable(size_t v, Kind k, int32_t net) const
  {
    return !hardBlocked(v, k, net) && count_[v][k] > 0;
  }
  uint8_t routedCount(size_t v, Kind k) const { return count_[v][k]; }
  uint16_t history(size_t v) const { return history_[v]; }
  void bumpHistory(size_t v)
  {
    if (history_[v] < 0xffff) {
      history_[v]++;
    }
  }

 private:
  std::vector<uint8_t> fixed_;
  std::vector<uint8_t> pin_mask_;
  std::vector<int32_t> pin_owner_;
  std::vector<std::array<uint8_t, kNumKinds>> count_;
  std::vector<uint16_t> history_;
};

// Routing area R of one search: allowed global routing tiles per layer.
struct Region
{
  int tiles_x = 1;
  int tiles_y = 1;
  std::vector<std::vector<uint8_t>> allowed;  // [z][tx * tiles_y + ty]

  void init(int num_layers, int tx, int ty)
  {
    tiles_x = tx;
    tiles_y = ty;
    allowed.assign(num_layers, std::vector<uint8_t>(size_t(tx) * ty, 0));
  }
  std::vector<std::pair<int, size_t>> touched;
  // Clears the allowed tiles in O(#allowed).
  void reset()
  {
    for (const auto& [z, index] : touched) {
      allowed[z][index] = 0;
    }
    touched.clear();
    bx_lo = by_lo = std::numeric_limits<int>::max();
    bx_hi = by_hi = -1;
  }
  bool allowedTile(int z, int tx, int ty) const
  {
    return z >= 0 && z < static_cast<int>(allowed.size()) && tx >= 0
           && tx < tiles_x && ty >= 0 && ty < tiles_y
           && allowed[z][size_t(tx) * tiles_y + ty] != 0;
  }
  // Bounding box of all allowed tiles.
  int bx_lo = std::numeric_limits<int>::max();
  int by_lo = std::numeric_limits<int>::max();
  int bx_hi = -1;
  int by_hi = -1;

  void allow(int z, int tx, int ty)
  {
    if (z >= 0 && z < static_cast<int>(allowed.size()) && tx >= 0
        && tx < tiles_x && ty >= 0 && ty < tiles_y) {
      uint8_t& cell = allowed[z][size_t(tx) * tiles_y + ty];
      if (!cell) {
        cell = 1;
        touched.emplace_back(z, size_t(tx) * tiles_y + ty);
      }
      bx_lo = std::min(bx_lo, tx);
      by_lo = std::min(by_lo, ty);
      bx_hi = std::max(bx_hi, tx);
      by_hi = std::max(by_hi, ty);
    }
  }
  bool contains(const TrackLayer& l, int z, int i, int j) const
  {
    const int t = l.track_tile[i];
    const int p = l.pos_tile[j];
    const int tx = l.horizontal ? p : t;
    const int ty = l.horizontal ? t : p;
    return allowed[z][size_t(tx) * tiles_y + ty] != 0;
  }
};

// Target cluster for the future cost: a rectangle on a layer range.
struct TargetBox
{
  Coord xlo, ylo, xhi, yhi;
  int zlo, zhi;
};

struct SearchRequest
{
  int32_t net = -1;
  const Region* region = nullptr;
  // Source vertices with their initial cost (0 on the routed tree, stub
  // cost at pin attachments).
  std::vector<std::pair<size_t, Cost>> sources;
  // Target vertices: (vertex, terminal cost, target id).
  struct Target
  {
    size_t vertex;
    Cost cost;
    int id;
  };
  std::vector<Target> targets;
  std::vector<TargetBox> target_boxes;
  // Rip-up mode: vertices blocked only by routed shapes are usable at
  // ripup_cost * (1 + history).
  bool ripup = false;
  Cost ripup_cost = 0;
  // Vertices never usable in this search (e.g. own pin landing conflicts).
  bool use_intervals = true;
};

struct SearchResult
{
  bool found = false;
  Cost cost = kInfCost;
  int target_id = -1;
  size_t target_vertex = 0;
  std::vector<size_t> path;  // from target back to a source
  // Statistics.
  int64_t heap_pops = 0;
  int64_t heap_pushes = 0;
  int64_t settled = 0;
};

// Dijkstra with future cost on the track graph (Sec. 4.1). With
// use_intervals, a preferred-direction run whose reduced cost is zero is
// labelled in one sweep instead of vertex by vertex (Algorithm 4); the
// result cost equals the vertex-by-vertex search.
class PathSearch
{
 public:
  PathSearch(const TrackGraph& graph, const FastGrid& grid)
      : graph_(graph), grid_(grid)
  {
  }
  SearchResult run(const SearchRequest& request);

 private:
  enum Parent : uint8_t
  {
    kSource = 0,
    kFromLower = 1,       // along: came from j - 1
    kFromHigher = 2,      // along: came from j + 1
    kFromTrackBelow = 3,  // jog: came from i - 1
    kFromTrackAbove = 4,  // jog: came from i + 1
    kFromLayerBelow = 5,  // via from z - 1
    kFromLayerAbove = 6,  // via from z + 1
  };
  struct Item
  {
    Cost key;
    Cost g;
    size_t v;
    int8_t dir;  // 0 vertex label, +1 / -1 ray, 2 sink
    int32_t aux;
    bool operator>(const Item& o) const
    {
      if (key != o.key) {
        return key > o.key;
      }
      if (g != o.g) {
        return g < o.g;  // deeper first on ties
      }
      return v > o.v;
    }
  };

  Cost futureCost(int z, Coord x, Coord y) const;
  Cost vertexPenalty(size_t v, Kind k) const;
  bool enter(size_t v, Kind k) const;
  void push(const Item& item);
  bool relax(size_t from_v, size_t to_v, Cost g, Parent parent, int8_t dir);
  void settle(int z, int i, int j, size_t v, Cost g);

  const TrackGraph& graph_;
  const FastGrid& grid_;
  const SearchRequest* req_ = nullptr;
  SearchResult* res_ = nullptr;
  Cost min_via_cost_ = 0;

  // Search window: per layer track range [i0, i0 + ni) and position range
  // [j0, j0 + nj) covering the region; local index = wbase + di * nj + dj.
  struct Window
  {
    int i0 = 0, ni = 0, j0 = 0, nj = 0;
    size_t wbase = 0;
  };
  bool local(int z, int i, int j, size_t& index) const
  {
    const Window& w = window_[z];
    const int di = i - w.i0;
    const int dj = j - w.j0;
    if (di < 0 || di >= w.ni || dj < 0 || dj >= w.nj) {
      return false;
    }
    index = w.wbase + static_cast<size_t>(di) * w.nj + dj;
    return true;
  }
  void setupWindow(const SearchRequest& request);
  std::vector<Window> window_;

  std::vector<Cost> dist_;
  std::vector<uint8_t> parent_;
  std::vector<uint32_t> stamp_;
  std::vector<uint32_t> source_stamp_;
  std::vector<uint32_t> target_stamp_;
  std::vector<int32_t> target_index_;
  uint32_t epoch_ = 0;
  std::vector<Item> heap_;
};

}  // namespace drt::br
