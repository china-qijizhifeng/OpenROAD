// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "BonnRouteCapacity.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

#include "odb/geom.h"

namespace grt::bonnroute {

std::vector<double> expectedMaxColumnLoad(int rows,
                                          int columns,
                                          int length,
                                          int max_count)
{
  rows = std::max(rows, 1);
  columns = std::max(columns, 1);
  length = std::clamp(length, 1, columns);
  max_count = std::max(max_count, 0);
  // A state counts the rows by how many more columns their current run
  // covers, a[0] = free rows and a[j] (j >= 1) = rows covered in the
  // current column and j - 1 columns after it. Rows are symmetric, so the
  // number of placements follows from these counts and binomials.
  std::vector<std::vector<int>> states;
  std::map<std::vector<int>, int> index;
  std::vector<int> current(length, 0);
  std::function<void(int, int)> enumerate = [&](int part, int left) {
    if (part == length - 1) {
      current[part] = left;
      index[current] = static_cast<int>(states.size());
      states.push_back(current);
      return;
    }
    for (int v = 0; v <= left; v++) {
      current[part] = v;
      enumerate(part + 1, left - v);
    }
  };
  enumerate(0, rows);
  std::vector<std::vector<double>> binom(rows + 1,
                                         std::vector<double>(rows + 1, 0.0));
  for (int n = 0; n <= rows; n++) {
    binom[n][0] = 1.0;
    for (int k = 1; k <= n; k++) {
      binom[n][k] = binom[n - 1][k - 1] + (k <= n - 1 ? binom[n - 1][k] : 0.0);
    }
  }
  std::vector<int> free_state(length, 0);
  free_state[0] = rows;
  const int all_free = index.at(free_state);
  const size_t num_states = states.size();
  // count[m][k]: placements of k runs with at most m chosen vertices in
  // every column; count[rows] counts all placements.
  std::vector<std::vector<double>> count(
      rows + 1, std::vector<double>(max_count + 1, 0.0));
  std::vector<double> cur(num_states * (max_count + 1));
  std::vector<double> next(cur.size());
  std::vector<int> successor(length);
  for (int bound = 0; bound <= rows; bound++) {
    std::ranges::fill(cur, 0.0);
    cur[static_cast<size_t>(all_free) * (max_count + 1)] = 1.0;
    for (int col = 0; col < columns; col++) {
      std::ranges::fill(next, 0.0);
      const bool can_start = col + length <= columns;
      for (size_t s = 0; s < num_states; s++) {
        const std::vector<int>& a = states[s];
        const int covered = rows - a[0];
        if (covered > bound) {
          continue;
        }
        const int max_start = can_start ? std::min(a[0], bound - covered) : 0;
        for (int start = 0; start <= max_start; start++) {
          if (length == 1) {
            successor[0] = a[0];
          } else {
            successor[0] = a[0] - start + a[1];
            for (int j = 1; j + 1 < length; j++) {
              successor[j] = a[j + 1];
            }
            successor[length - 1] = start;
          }
          const size_t t = index.at(successor);
          const double ways = binom[a[0]][start];
          for (int k = 0; k + start <= max_count; k++) {
            const double v = cur[s * (max_count + 1) + k];
            if (v != 0.0) {
              next[t * (max_count + 1) + k + start] += v * ways;
            }
          }
        }
      }
      std::swap(cur, next);
    }
    for (int k = 0; k <= max_count; k++) {
      count[bound][k]
          = cur[static_cast<size_t>(all_free) * (max_count + 1) + k];
    }
  }
  std::vector<double> result(max_count + 1, 0.0);
  for (int k = 1; k <= max_count; k++) {
    const double total = count[rows][k];
    if (total <= 0.0) {
      result[k] = result[k - 1];
      continue;
    }
    // E[max] = sum over m >= 0 of P(max > m).
    double expected = 0.0;
    for (int bound = 0; bound < rows; bound++) {
      expected += 1.0 - count[bound][k] / total;
    }
    result[k] = expected;
  }
  return result;
}

double interpolate(const std::vector<double>& table, const double count)
{
  if (table.empty() || count <= 0.0) {
    return table.empty() ? 0.0 : table.front();
  }
  const auto last = static_cast<double>(table.size() - 1);
  if (count >= last) {
    return table.back();
  }
  const auto i = static_cast<size_t>(count);
  const double f = count - static_cast<double>(i);
  return table[i] * (1.0 - f) + table[i + 1] * f;
}

std::pair<int64_t, int64_t> rectilinearSpanningLength(
    const std::vector<odb::Point>& points)
{
  const size_t n = points.size();
  if (n < 2) {
    return {0, 0};
  }
  std::vector<int64_t> best(n, std::numeric_limits<int64_t>::max());
  std::vector<size_t> from(n, 0);
  std::vector<char> done(n, 0);
  best[0] = 0;
  int64_t horizontal = 0;
  int64_t vertical = 0;
  for (size_t step = 0; step < n; step++) {
    size_t u = n;
    for (size_t i = 0; i < n; i++) {
      if (!done[i] && (u == n || best[i] < best[u])) {
        u = i;
      }
    }
    done[u] = 1;
    if (step > 0) {
      horizontal += std::abs(static_cast<int64_t>(points[u].x())
                             - points[from[u]].x());
      vertical += std::abs(static_cast<int64_t>(points[u].y())
                           - points[from[u]].y());
    }
    for (size_t i = 0; i < n; i++) {
      if (done[i]) {
        continue;
      }
      const int64_t d
          = std::abs(static_cast<int64_t>(points[u].x()) - points[i].x())
            + std::abs(static_cast<int64_t>(points[u].y()) - points[i].y());
      if (d < best[i]) {
        best[i] = d;
        from[i] = u;
      }
    }
  }
  return {horizontal, vertical};
}

TrackLayer::TrackLayer(const int direction,
                       std::vector<int> tracks,
                       std::vector<int> vertices)
    : direction_(direction),
      tracks_(std::move(tracks)),
      vertices_(std::move(vertices))
{
  const size_t bits = tracks_.size() * vertices_.size();
  blocked_.assign((bits + 63) / 64, 0);
}

bool TrackLayer::usable(const int track, const int vertex) const
{
  const size_t bit = static_cast<size_t>(track) * vertices_.size() + vertex;
  return ((blocked_[bit / 64] >> (bit % 64)) & 1) == 0;
}

void TrackLayer::block(const int track, int first, int last)
{
  first = std::max(first, 0);
  last = std::min(last, numVertices() - 1);
  if (first > last) {
    return;
  }
  const size_t row = static_cast<size_t>(track) * vertices_.size();
  size_t begin = row + first;
  const size_t end = row + last + 1;
  while (begin < end) {
    const size_t word = begin / 64;
    const size_t offset = begin % 64;
    const size_t count = std::min<size_t>(64 - offset, end - begin);
    const uint64_t mask
        = (count == 64 ? ~uint64_t{0} : ((uint64_t{1} << count) - 1)) << offset;
    blocked_[word] |= mask;
    begin += count;
  }
}

int TrackLayer::countUsable(const int track, int first, int last) const
{
  first = std::max(first, 0);
  last = std::min(last, numVertices());
  if (first >= last) {
    return 0;
  }
  const size_t row = static_cast<size_t>(track) * vertices_.size();
  size_t begin = row + first;
  const size_t end = row + last;
  int blocked = 0;
  while (begin < end) {
    const size_t word = begin / 64;
    const size_t offset = begin % 64;
    const size_t count = std::min<size_t>(64 - offset, end - begin);
    const uint64_t mask
        = (count == 64 ? ~uint64_t{0} : ((uint64_t{1} << count) - 1)) << offset;
    blocked += std::popcount(blocked_[word] & mask);
    begin += count;
  }
  return (last - first) - blocked;
}

namespace {
std::pair<int, int> coordRange(const std::vector<int>& coords,
                               const int lo,
                               const int hi)
{
  const auto first = std::ranges::lower_bound(coords, lo) - coords.begin();
  const auto last = std::ranges::upper_bound(coords, hi) - coords.begin() - 1;
  return {static_cast<int>(first), static_cast<int>(last)};
}

int nearest(const std::vector<int>& coords, const int coord)
{
  if (coords.empty()) {
    return -1;
  }
  const auto it = std::ranges::lower_bound(coords, coord);
  if (it == coords.end()) {
    return static_cast<int>(coords.size()) - 1;
  }
  const int i = static_cast<int>(it - coords.begin());
  if (i > 0 && coord - coords[i - 1] <= *it - coord) {
    return i - 1;
  }
  return i;
}

int exactIndex(const std::vector<int>& coords, const int coord)
{
  const auto it = std::ranges::lower_bound(coords, coord);
  return (it != coords.end() && *it == coord)
             ? static_cast<int>(it - coords.begin())
             : -1;
}
}  // namespace

std::pair<int, int> TrackLayer::trackRange(const int lo, const int hi) const
{
  return coordRange(tracks_, lo, hi);
}

std::pair<int, int> TrackLayer::vertexRange(const int lo, const int hi) const
{
  return coordRange(vertices_, lo, hi);
}

int TrackLayer::lowerVertex(const int coord) const
{
  return static_cast<int>(std::ranges::lower_bound(vertices_, coord)
                          - vertices_.begin());
}

int TrackLayer::nearestTrack(const int coord) const
{
  return nearest(tracks_, coord);
}

int TrackLayer::nearestVertex(const int coord) const
{
  return nearest(vertices_, coord);
}

int CapacityLayer::parallelSpacing(const int wire_width, const int length) const
{
  int w = static_cast<int>(spacing_widths.size()) - 1;
  while (w > 0 && spacing_widths[w] >= wire_width) {
    w--;
  }
  if (length == 0) {
    return spacing[w][0];
  }
  int l = static_cast<int>(spacing_lengths.size()) - 1;
  while (l > 0 && spacing_lengths[l] >= length) {
    l--;
  }
  return spacing[w][l];
}

CapacityEstimator::CapacityEstimator(CapacityGrid grid,
                                     const CapacityEstimatorOptions& options)
    : grid_(std::move(grid)), options_(options)
{
  num_layers_ = grid_.numLayers();
  x_size_ = grid_.size(0);
  y_size_ = grid_.size(1);
  min_layer_ = num_layers_;
  for (int l = 0; l < num_layers_; l++) {
    if (grid_.layers[l].routable) {
      min_layer_ = l;
      break;
    }
  }
  const size_t n = static_cast<size_t>(num_layers_) * x_size_ * y_size_;
  stacked_.assign(n, 0.0);
  intra_.assign(n, 0.0);
}

// Sec. 2.5, paragraph 3: "we count the usable vertices of the track graph
// (cf. Section 3.5) in the tile areas of v and w between the coordinates of
// the tile centers c_v and c_w in preferred direction after extending each
// blockage by a small constant value in preferred direction."
void CapacityEstimator::buildTrackGraph()
{
  layers_.clear();
  layers_.reserve(num_layers_);
  for (int l = 0; l < num_layers_; l++) {
    const CapacityLayer& layer = grid_.layers[l];
    const int direction = layer.direction;
    // "The intersection points of routing tracks with tracks projected from
    // neighboring wiring layers define the vertices of our track graph."
    std::vector<int> vertices;
    for (const int n : {l - 1, l + 1}) {
      if (n >= 0 && n < num_layers_ && grid_.layers[n].direction != direction) {
        const std::vector<int>& c = grid_.layers[n].tracks;
        vertices.insert(vertices.end(), c.begin(), c.end());
      }
    }
    if (vertices.empty()) {
      const int lo
          = direction == kHorizontal ? grid_.die.xMin() : grid_.die.yMin();
      const int hi
          = direction == kHorizontal ? grid_.die.xMax() : grid_.die.yMax();
      const int pitch = std::max(1, layer.pitch);
      for (int c = lo + pitch / 2; c <= hi; c += pitch) {
        vertices.push_back(c);
      }
    }
    std::ranges::sort(vertices);
    vertices.erase(std::ranges::unique(vertices).begin(), vertices.end());
    layers_.emplace_back(direction, layer.tracks, std::move(vertices));
  }

  for (int l = min_layer_; l < num_layers_; l++) {
    const CapacityLayer& layer = grid_.layers[l];
    const int direction = layer.direction;
    TrackLayer& tracks = layers_[l];
    int min_edge_length = std::numeric_limits<int>::max();
    for (int e = 0; e + 1 < grid_.size(direction); e++) {
      min_edge_length
          = std::min(min_edge_length, grid_.edgeLength(direction, e));
    }
    for (const odb::Rect& obs : layer.obstacles) {
      const int width = std::min(obs.dx(), obs.dy());
      const int along_lo = direction == kHorizontal ? obs.xMin() : obs.yMin();
      const int along_hi = direction == kHorizontal ? obs.xMax() : obs.yMax();
      const int across_lo = direction == kHorizontal ? obs.yMin() : obs.xMin();
      const int across_hi = direction == kHorizontal ? obs.yMax() : obs.xMax();
      // Keep-out of a standard wire's stick figure, as CUGR uses it across
      // the tracks; extended by the same amount along them.
      const int keep_out
          = layer.parallelSpacing(
                width, std::min(min_edge_length, along_hi - along_lo))
            + layer.width / 2 - 1;
      const int extension = static_cast<int>(
          std::lround(options_.blockage_extension * std::max(0, keep_out)));
      const auto [t0, t1]
          = tracks.trackRange(across_lo - keep_out, across_hi + keep_out);
      const auto [v0, v1]
          = tracks.vertexRange(along_lo - extension, along_hi + extension);
      if (v0 > v1) {
        continue;
      }
      for (int t = t0; t <= t1; t++) {
        tracks.block(t, v0, v1);
      }
    }
  }
  totals_ = CapacityTotals();
  before_single_tile_ = vertexCapacities();
}

std::vector<double> CapacityEstimator::vertexCapacities() const
{
  std::vector<double> capacity(
      static_cast<size_t>(num_layers_) * x_size_ * y_size_, 0.0);
  for (int l = min_layer_; l < num_layers_; l++) {
    const TrackLayer& tracks = layers_[l];
    const int direction = tracks.direction();
    const int rows = grid_.size(1 - direction);
    const int edges = grid_.size(direction) - 1;
#pragma omp parallel for schedule(dynamic, 8)
    for (int r = 0; r < rows; r++) {
      const int lo = grid_.gridline(1 - direction, r);
      const int hi = grid_.gridline(1 - direction, r + 1);
      // Tracks on the upper gridline belong to the next row.
      const auto [t0, t1] = tracks.trackRange(lo, r + 1 == rows ? hi : hi - 1);
      for (int e = 0; e < edges; e++) {
        const int first = tracks.lowerVertex(grid_.center(direction, e));
        const int last = tracks.lowerVertex(grid_.center(direction, e + 1));
        const int n = last - first;
        double cap = 0.0;
        if (n <= 0) {
          cap = std::max(0, t1 - t0 + 1);
        } else {
          // "Then we divide the result by the number of vertices on a track
          // between c_v and c_w."
          int usable = 0;
          for (int t = t0; t <= t1; t++) {
            usable += tracks.countUsable(t, first, last);
          }
          cap = static_cast<double>(usable) / n;
        }
        const int x = direction == kHorizontal ? e : r;
        const int y = direction == kHorizontal ? r : e;
        capacity[flat(l, x, y)] = cap;
      }
    }
  }
  return capacity;
}

int CapacityEstimator::lowestLayer(const int direction,
                                   const int min_layer,
                                   const int max_layer) const
{
  for (int l = std::max(min_layer, min_layer_);
       l <= std::min(max_layer, num_layers_ - 1);
       l++) {
    if (layers_[l].direction() == direction) {
      return l;
    }
  }
  return -1;
}

void CapacityEstimator::addNets(const std::vector<EstimatorNet>& nets)
{
  for (const EstimatorNet& net : nets) {
    if (net.single_tile) {
      if (!options_.single_tile_nets) {
        continue;
      }
      if (routeSingleTile(net)) {
        totals_.single_tile_routed++;
      } else {
        totals_.single_tile_failed++;
        // Fall back to the Steiner-length estimate of (c).
        addIntraTileNet(net, false);
      }
    } else if (options_.intra_tile) {
      addIntraTileNet(net, true);
    }
  }
}

// Sec. 2.5, paragraph 4: "prior to global routing capacity estimation we
// run our detailed router on nets whose pins are contained in a single
// tile, allowing it to route within a slightly larger area than that of the
// tile." Here a shortest-path (path composition, Algorithm 1) router on the
// track graph of the lowest routable layers stands in for the detailed
// router; the vertices of its route become unusable.
bool CapacityEstimator::routeSingleTile(const EstimatorNet& net)
{
  const int low = std::max(min_layer_, net.min_layer);
  const int high = std::min({low + std::max(1, options_.single_tile_layers) - 1,
                             net.max_layer,
                             num_layers_ - 1});
  if (high < low || net.tile_x < 0 || net.tile_y < 0) {
    return false;
  }
  const int x0 = grid_.gridline(0, net.tile_x);
  const int x1 = grid_.gridline(0, net.tile_x + 1);
  const int y0 = grid_.gridline(1, net.tile_y);
  const int y1 = grid_.gridline(1, net.tile_y + 1);
  const int margin = static_cast<int>(
      std::lround(options_.single_tile_margin * std::max(x1 - x0, y1 - y0)));
  const int region[2][2]
      = {{x0 - margin, x1 + margin}, {y0 - margin, y1 + margin}};

  const int num = high - low + 1;
  std::vector<int> t_first(num), t_count(num), v_first(num), v_count(num);
  std::vector<int64_t> offset(num + 1, 0);
  for (int j = 0; j < num; j++) {
    const TrackLayer& tracks = layers_[low + j];
    const int d = tracks.direction();
    const auto [t0, t1] = tracks.trackRange(region[1 - d][0], region[1 - d][1]);
    const auto [v0, v1] = tracks.vertexRange(region[d][0], region[d][1]);
    if (t0 > t1 || v0 > v1) {
      return false;
    }
    t_first[j] = t0;
    t_count[j] = t1 - t0 + 1;
    v_first[j] = v0;
    v_count[j] = v1 - v0 + 1;
    offset[j + 1] = offset[j] + static_cast<int64_t>(t_count[j]) * v_count[j];
  }
  const int64_t num_nodes = offset[num];
  struct Node
  {
    int j;
    int t;
    int v;
  };
  auto decode = [&](int64_t id) {
    int j = 0;
    while (offset[j + 1] <= id) {
      j++;
    }
    const int64_t rem = id - offset[j];
    return Node{j,
                t_first[j] + static_cast<int>(rem / v_count[j]),
                v_first[j] + static_cast<int>(rem % v_count[j])};
  };
  auto encode = [&](int j, int t, int v) -> int64_t {
    if (t < t_first[j] || t >= t_first[j] + t_count[j] || v < v_first[j]
        || v >= v_first[j] + v_count[j]) {
      return -1;
    }
    return offset[j] + static_cast<int64_t>(t - t_first[j]) * v_count[j]
           + (v - v_first[j]);
  };

  // Terminals: each pin is the set of its access points' nearest vertices
  // (the zero-cost clique of Sec. 2.1 applies inside the tile as well).
  std::unordered_map<int64_t, std::vector<int>> pins_at;
  std::vector<std::vector<int64_t>> terminals(net.pins.size());
  for (size_t p = 0; p < net.pins.size(); p++) {
    for (const auto& [point, ap_layer] : net.pins[p].points) {
      if (ap_layer > high) {
        continue;
      }
      const int j = std::max(ap_layer, low) - low;
      const TrackLayer& tracks = layers_[low + j];
      const int d = tracks.direction();
      const int across = d == kHorizontal ? point.y() : point.x();
      const int along = d == kHorizontal ? point.x() : point.y();
      const int64_t id
          = encode(j, tracks.nearestTrack(across), tracks.nearestVertex(along));
      if (id < 0) {
        continue;
      }
      if (std::ranges::find(terminals[p], id) == terminals[p].end()) {
        terminals[p].push_back(id);
        pins_at[id].push_back(static_cast<int>(p));
      }
    }
    if (terminals[p].empty()) {
      return false;
    }
  }

  const int64_t via_cost = 4LL * std::max(1, grid_.layers[low].pitch);
  constexpr int64_t kInf = std::numeric_limits<int64_t>::max();
  std::vector<int64_t> dist(num_nodes, kInf);
  std::vector<int64_t> parent(num_nodes, -1);
  std::vector<char> in_tree(num_nodes, 0);
  std::vector<char> connected(net.pins.size(), 0);
  size_t remaining = net.pins.size();
  using Item = std::pair<int64_t, int64_t>;
  std::priority_queue<Item, std::vector<Item>, std::greater<>> heap;
  auto add_source = [&](int64_t id) {
    in_tree[id] = 1;
    dist[id] = 0;
    parent[id] = -1;
    heap.emplace(0, id);
  };
  std::function<void(int64_t)> connect_at = [&](int64_t id) {
    auto it = pins_at.find(id);
    if (it == pins_at.end()) {
      return;
    }
    for (const int p : it->second) {
      if (connected[p]) {
        continue;
      }
      connected[p] = 1;
      remaining--;
      for (const int64_t w : terminals[p]) {
        if (!in_tree[w]) {
          add_source(w);
          connect_at(w);
        }
      }
    }
  };
  add_source(terminals[0].front());
  connect_at(terminals[0].front());
  auto passable = [&](int64_t id) {
    const Node n = decode(id);
    return layers_[low + n.j].usable(n.t, n.v) || pins_at.contains(id);
  };
  std::vector<int64_t> used;
  while (!heap.empty() && remaining > 0) {
    const auto [d, u] = heap.top();
    heap.pop();
    if (d > dist[u]) {
      continue;
    }
    if (!in_tree[u] && pins_at.contains(u)
        && std::ranges::any_of(pins_at[u],
                               [&](int p) { return !connected[p]; })) {
      // Algorithm 1, step 5: add the path to the tree.
      std::vector<int64_t> path;
      int64_t cur = u;
      while (cur >= 0 && !in_tree[cur]) {
        path.push_back(cur);
        cur = parent[cur];
      }
      if (cur >= 0) {
        used.push_back(cur);  // the attachment point carries wire too
      }
      for (const int64_t id : path) {
        used.push_back(id);
        add_source(id);
      }
      for (const int64_t id : path) {
        connect_at(id);
      }
      continue;
    }
    const Node n = decode(u);
    const TrackLayer& tracks = layers_[low + n.j];
    auto relax = [&](int64_t to, int64_t cost) {
      if (to < 0 || !passable(to)) {
        return;
      }
      if (d + cost < dist[to]) {
        dist[to] = d + cost;
        parent[to] = u;
        heap.emplace(dist[to], to);
      }
    };
    for (const int dv : {-1, 1}) {
      const int v = n.v + dv;
      if (v >= 0 && v < tracks.numVertices()) {
        relax(encode(n.j, n.t, v),
              std::abs(tracks.vertices()[v] - tracks.vertices()[n.v]));
      }
    }
    for (const int dj : {-1, 1}) {
      const int j = n.j + dj;
      if (j < 0 || j >= num) {
        continue;
      }
      const TrackLayer& other = layers_[low + j];
      if (other.direction() == tracks.direction()) {
        continue;
      }
      // The vertex lies on a track of the other layer, and this layer's
      // track is one of the other layer's vertex coordinates.
      const int t = exactIndex(other.tracks(), tracks.vertices()[n.v]);
      const int v = exactIndex(other.vertices(), tracks.tracks()[n.t]);
      if (t >= 0 && v >= 0) {
        relax(encode(j, t, v), via_cost);
      }
    }
  }
  if (remaining > 0) {
    return false;
  }
  std::ranges::sort(used);
  used.erase(std::ranges::unique(used).begin(), used.end());
  for (const int64_t id : used) {
    const Node n = decode(id);
    layers_[low + n.j].block(n.t, n.v, n.v);
  }
  totals_.single_tile_vertices += static_cast<int64_t>(used.size());
  return true;
}

void CapacityEstimator::addIntraTile(const int layer,
                                     const int tile_x,
                                     const int tile_y,
                                     const double length)
{
  if (layer < 0 || length <= 0.0) {
    return;
  }
  const int direction = layers_[layer].direction();
  const int i = direction == kHorizontal ? tile_x : tile_y;
  const int size = grid_.size(direction);
  // Half of the tile lies in the span of each flanking edge.
  for (const int e : {i - 1, i}) {
    if (e < 0 || e + 1 >= size) {
      continue;
    }
    const int edge_length = grid_.edgeLength(direction, e);
    if (edge_length <= 0) {
      continue;
    }
    const int x = direction == kHorizontal ? e : tile_x;
    const int y = direction == kHorizontal ? tile_y : e;
    intra_[flat(layer, x, y)] += 0.5 * length / edge_length;
  }
}

// Sec. 2.5, paragraph 4: "intra-tile connection lengths of longer nets are
// estimated by their Steiner lengths and edge capacities are reduced
// accordingly, in a similar way as in Wei et al. [2012]." The pins of a net
// inside one tile are joined with the tile center, where the global route
// attaches; the horizontal and vertical lengths go to the lowest routable
// layer of that direction.
void CapacityEstimator::addIntraTileNet(const EstimatorNet& net,
                                        const bool include_center)
{
  std::map<std::pair<int, int>, std::vector<odb::Point>> tiles;
  for (const EstimatorPin& pin : net.pins) {
    int cx = pin.cell_x;
    int cy = pin.cell_y;
    if (net.single_tile) {
      cx = net.tile_x;
      cy = net.tile_y;
    }
    if (cx < 0 || cy < 0 || pin.points.empty()) {
      continue;
    }
    const odb::Point center(grid_.center(0, cx), grid_.center(1, cy));
    const odb::Point* best = nullptr;
    int64_t best_distance = std::numeric_limits<int64_t>::max();
    for (const auto& [point, layer] : pin.points) {
      const int64_t d
          = std::abs(static_cast<int64_t>(point.x()) - center.x())
            + std::abs(static_cast<int64_t>(point.y()) - center.y());
      if (d < best_distance) {
        best_distance = d;
        best = &point;
      }
    }
    tiles[{cx, cy}].push_back(*best);
  }
  for (auto& [tile, points] : tiles) {
    if (include_center) {
      points.emplace_back(grid_.center(0, tile.first),
                          grid_.center(1, tile.second));
    }
    const auto [horizontal, vertical] = rectilinearSpanningLength(points);
    addIntraTile(lowestLayer(kHorizontal, net.min_layer, net.max_layer),
                 tile.first,
                 tile.second,
                 static_cast<double>(horizontal));
    addIntraTile(lowestLayer(kVertical, net.min_layer, net.max_layer),
                 tile.first,
                 tile.second,
                 static_cast<double>(vertical));
  }
}

// Sec. 2.5, paragraphs 5-6: a stacked via from l - 1 to l + 1 consumes space
// on l. The expected maximum column load of k random disjoint pads in the
// lattice of one edge span approximates the capacity it costs ("the
// expected capacity reduction is sublinear in the number of stacked
// vias"); pads blocking several tracks scale it. CUGR's demand model
// already charges every pad linearly (GridGraph's via-pad flank demand), so
// only the excess of the statistical estimate over that linear charge
// reduces the capacity here.
void CapacityEstimator::addStackedVias(const std::vector<double>& stack_counts)
{
  if (!options_.stacked_vias) {
    return;
  }
  const bool linear_charged = grid_.linear_via_demand;
  for (int l = std::max(1, min_layer_); l + 1 < num_layers_; l++) {
    const TrackLayer& tracks = layers_[l];
    const int direction = tracks.direction();
    const int rows = grid_.size(1 - direction);
    const int cells = grid_.size(direction);
    if (cells < 2 || tracks.numVertices() < 2) {
      continue;
    }
    double max_count = 0.0;
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cells; c++) {
        const int x = direction == kHorizontal ? c : r;
        const int y = direction == kHorizontal ? r : c;
        max_count = std::max(max_count, stack_counts[flat(l, x, y)]);
      }
    }
    if (max_count <= 0.0) {
      continue;
    }
    const double pad_length = grid_.layers[l].through_pad_length;
    const double along = grid_.layers[l].through_pad_along;
    const double scale = std::max(1.0, grid_.layers[l].through_pad_tracks);
    // Normalized region: typical tracks per row and vertices per span.
    std::vector<int> row_tracks;
    for (int r = 0; r < rows; r++) {
      const int lo = grid_.gridline(1 - direction, r);
      const int hi = grid_.gridline(1 - direction, r + 1);
      const auto [t0, t1] = tracks.trackRange(lo, r + 1 == rows ? hi : hi - 1);
      row_tracks.push_back(std::max(0, t1 - t0 + 1));
    }
    std::vector<int> span_vertices;
    for (int e = 0; e + 1 < cells; e++) {
      span_vertices.push_back(tracks.lowerVertex(grid_.center(direction, e + 1))
                              - tracks.lowerVertex(grid_.center(direction, e)));
    }
    std::ranges::nth_element(row_tracks,
                             row_tracks.begin() + row_tracks.size() / 2);
    std::ranges::nth_element(span_vertices,
                             span_vertices.begin() + span_vertices.size() / 2);
    const int typical_tracks = row_tracks[row_tracks.size() / 2];
    const int columns
        = std::clamp(span_vertices[span_vertices.size() / 2], 1, 40);
    const double spacing = static_cast<double>(tracks.vertices().back()
                                               - tracks.vertices().front())
                           / (tracks.numVertices() - 1);
    const int length = std::clamp(
        static_cast<int>(std::lround(along / std::max(1.0, spacing))), 1, 3);
    const int lattice_rows = std::clamp(
        static_cast<int>(std::lround(typical_tracks / scale)), 1, 20);
    const int table_max = std::min(48, static_cast<int>(std::ceil(max_count)));
    const std::vector<double> table
        = expectedMaxColumnLoad(lattice_rows, columns, length, table_max);
    auto expected_loss = [&](double k) {
      double load = interpolate(table, k);
      if (k > table_max) {
        load += (k - table_max) * length / columns;
      }
      return scale * std::min<double>(load, lattice_rows);
    };
    auto flank_sum = [&](int c) {
      int sum = 0;
      if (c > 0) {
        sum += grid_.edgeLength(direction, c - 1);
      }
      if (c + 1 < cells) {
        sum += grid_.edgeLength(direction, c);
      }
      return sum;
    };
    for (int r = 0; r < rows; r++) {
      for (int e = 0; e + 1 < cells; e++) {
        const int x = direction == kHorizontal ? e : r;
        const int y = direction == kHorizontal ? r : e;
        const int xn = direction == kHorizontal ? e + 1 : r;
        const int yn = direction == kHorizontal ? r : e + 1;
        const double k0 = stack_counts[flat(l, x, y)];
        const double k1 = stack_counts[flat(l, xn, yn)];
        if (k0 + k1 <= 0.0) {
          continue;
        }
        // Half of each tile lies in this edge's span.
        const double statistical = expected_loss(0.5 * (k0 + k1));
        double linear = 0.0;
        if (linear_charged) {
          const int s0 = flank_sum(e);
          const int s1 = flank_sum(e + 1);
          linear += s0 > 0 ? k0 * pad_length / s0 : 0.0;
          linear += s1 > 0 ? k1 * pad_length / s1 : 0.0;
        }
        stacked_[flat(l, x, y)] = std::max(0.0, statistical - linear);
      }
    }
  }
}

std::vector<double> CapacityEstimator::capacities()
{
  std::vector<double> capacity = vertexCapacities();
  totals_.tracks.assign(num_layers_, 0.0);
  totals_.single_tile.assign(num_layers_, 0.0);
  totals_.intra_tile.assign(num_layers_, 0.0);
  totals_.stacked_vias.assign(num_layers_, 0.0);
  totals_.estimate.assign(num_layers_, 0.0);
  for (int l = min_layer_; l < num_layers_; l++) {
    const double adjustment = grid_.layers[l].adjustment;
    for (int x = 0; x < x_size_; x++) {
      for (int y = 0; y < y_size_; y++) {
        const size_t f = flat(l, x, y);
        totals_.tracks[l] += before_single_tile_[f];
        totals_.single_tile[l] += before_single_tile_[f] - capacity[f];
        double cap = capacity[f];
        const double intra = std::min(cap, intra_[f]);
        cap -= intra;
        const double stacked = std::min(cap, stacked_[f]);
        cap -= stacked;
        totals_.intra_tile[l] += intra;
        totals_.stacked_vias[l] += stacked;
        // The user layer adjustment, applied as GridGraph applies it to
        // CUGR's capacities (including its floor of one track).
        if (adjustment != 0.0) {
          const double before = cap;
          cap *= (1.0 - adjustment);
          if (before > 0.0 && adjustment != 1.0) {
            cap = std::max(cap, 1.0);
          }
        }
        capacity[f] = cap;
        totals_.estimate[l] += cap;
      }
    }
  }
  return capacity;
}

}  // namespace grt::bonnroute
