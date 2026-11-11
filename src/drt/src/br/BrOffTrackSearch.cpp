// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrOffTrackSearch.h"

#include <algorithm>
#include <cstdint>
#include <functional>
#include <queue>
#include <utility>
#include <vector>

namespace drt::br {

namespace {

constexpr int kFresh = 4;
constexpr int kStates = 5;
constexpr int64_t kInf = std::numeric_limits<int64_t>::max() / 4;

int exactIndex(const std::vector<GeoCoord>& v, const GeoCoord c)
{
  auto it = std::lower_bound(v.begin(), v.end(), c);
  return (it != v.end() && *it == c) ? static_cast<int>(it - v.begin()) : -1;
}

// Open intervals of the obstacles whose open interior meets the line
// `along` = c (horizontal line y = c if `rows`), sorted and merged where
// they overlap (touching open intervals stay apart: the common end point is
// free).
std::vector<std::pair<GeoCoord, GeoCoord>> lineIntervals(
    const std::vector<Box>& obstacles,
    const GeoCoord c,
    const bool rows)
{
  std::vector<std::pair<GeoCoord, GeoCoord>> iv;
  for (const Box& o : obstacles) {
    if (rows ? (c > o.ylo && c < o.yhi) : (c > o.xlo && c < o.xhi)) {
      iv.emplace_back(rows ? o.xlo : o.ylo, rows ? o.xhi : o.yhi);
    }
  }
  std::sort(iv.begin(), iv.end());
  std::vector<std::pair<GeoCoord, GeoCoord>> merged;
  for (const auto& [lo, hi] : iv) {
    if (lo >= hi) {
      continue;
    }
    if (!merged.empty() && lo < merged.back().second) {
      merged.back().second = std::max(merged.back().second, hi);
    } else {
      merged.emplace_back(lo, hi);
    }
  }
  return merged;
}

// True if [a, b] meets one of the open intervals.
bool meets(const std::vector<std::pair<GeoCoord, GeoCoord>>& iv,
           const GeoCoord a,
           const GeoCoord b)
{
  // First interval with hi > a.
  auto it = std::upper_bound(
      iv.begin(), iv.end(), a, [](GeoCoord v, const auto& p) {
        return v < p.second;
      });
  return it != iv.end() && it->first < b;
}

}  // namespace

bool OffTrackSearch::freeH(LayerIndex& l,
                           const int iy,
                           const int ix0,
                           const int ix1)
{
  if (!l.row_built[iy]) {
    l.rows[iy] = lineIntervals(*l.obstacles, l.ys[iy], true);
    l.row_built[iy] = 1;
  }
  return !meets(l.rows[iy], l.xs[ix0], l.xs[ix1]);
}

bool OffTrackSearch::freeV(LayerIndex& l,
                           const int ix,
                           const int iy0,
                           const int iy1)
{
  if (!l.col_built[ix]) {
    l.cols[ix] = lineIntervals(*l.obstacles, l.xs[ix], false);
    l.col_built[ix] = 1;
  }
  return !meets(l.cols[ix], l.ys[iy0], l.ys[iy1]);
}

OffTrackResult OffTrackSearch::run(const OffTrackRequest& req)
{
  OffTrackResult res;
  const int num_layers = static_cast<int>(req.layers.size());
  index_.assign(num_layers, LayerIndex());
  size_t total = 0;
  const Box& win = req.window;
  for (int z = 0; z < num_layers; z++) {
    const OffTrackLayer& in = req.layers[z];
    LayerIndex& l = index_[z];
    for (const GeoCoord x : in.grid.xs) {
      if (x >= win.xlo && x <= win.xhi) {
        l.xs.push_back(x);
      }
    }
    for (const GeoCoord y : in.grid.ys) {
      if (y >= win.ylo && y <= win.yhi) {
        l.ys.push_back(y);
      }
    }
    l.base = total;
    total += static_cast<size_t>(l.nx()) * l.ny();
    l.obstacles = &in.obstacles;
    l.rows.assign(l.ny(), {});
    l.cols.assign(l.nx(), {});
    l.row_built.assign(l.ny(), 0);
    l.col_built.assign(l.nx(), 0);
  }
  res.grid_vertices = static_cast<int64_t>(total);
  res.states = static_cast<int64_t>(total) * kStates;
  if (dist_.size() < total * kStates) {
    dist_.resize(total * kStates);
    parent_.resize(total * kStates);
    stamp_.resize(total * kStates, 0);
  }
  if (via_memo_.size() < total) {
    via_memo_.resize(total);
    via_stamp_.resize(total, 0);
    future_memo_.resize(total);
    future_stamp_.resize(total, 0);
  }
  if (++epoch_ == 0) {
    std::fill(stamp_.begin(), stamp_.end(), 0);
    std::fill(via_stamp_.begin(), via_stamp_.end(), 0);
    std::fill(future_stamp_.begin(), future_stamp_.end(), 0);
    epoch_ = 1;
  }
  for (int z = 0; z < num_layers; z++) {
    LayerIndex& l = index_[z];
    auto map_lines = [](const std::vector<GeoCoord>& from,
                        const std::vector<GeoCoord>& to,
                        std::vector<int>& out) {
      out.assign(from.size(), -1);
      size_t k = 0;
      for (size_t i = 0; i < from.size(); i++) {
        while (k < to.size() && to[k] < from[i]) {
          k++;
        }
        if (k < to.size() && to[k] == from[i]) {
          out[i] = static_cast<int>(k);
        }
      }
    };
    if (z + 1 < num_layers) {
      map_lines(l.xs, index_[z + 1].xs, l.up_x);
      map_lines(l.ys, index_[z + 1].ys, l.up_y);
    }
    if (z > 0) {
      map_lines(l.xs, index_[z - 1].xs, l.down_x);
      map_lines(l.ys, index_[z - 1].ys, l.down_y);
    }
  }
  auto dist = [&](const int64_t st) {
    return stamp_[st] == epoch_ ? dist_[st] : kInf;
  };
  auto set = [&](const int64_t st, const int64_t d, const int64_t parent) {
    stamp_[st] = epoch_;
    dist_[st] = d;
    parent_[st] = parent;
  };

  auto layer_of = [&](const size_t v) {
    int z = num_layers - 1;
    while (z > 0 && index_[z].base > v) {
      z--;
    }
    return z;
  };
  // Future cost (A*, as the on-track search of Sec. 4.1): l1 distance to
  // the lattice spanned by the target coordinates of a layer plus the via
  // costs to reach it. The lattice contains every target, so the bound is
  // admissible; it changes by at most an arc's cost along any arc, so it is
  // consistent, every state is settled with its exact distance, and targets
  // (future cost 0) are settled in increasing distance.
  std::vector<std::vector<GeoCoord>> target_xs(num_layers);
  std::vector<std::vector<GeoCoord>> target_ys(num_layers);
  int64_t unit_cost = std::numeric_limits<int64_t>::max();
  for (const OffTrackLayer& l : req.layers) {
    unit_cost = std::min({unit_cost, l.cost_x, l.cost_y});
  }
  unit_cost = std::max<int64_t>(0, unit_cost);
  for (const OffTrackTarget& t : req.targets) {
    if (t.layer >= 0 && t.layer < num_layers) {
      target_xs[t.layer].push_back(t.p.x);
      target_ys[t.layer].push_back(t.p.y);
    }
  }
  for (int z = 0; z < num_layers; z++) {
    std::sort(target_xs[z].begin(), target_xs[z].end());
    target_xs[z].erase(std::unique(target_xs[z].begin(), target_xs[z].end()),
                       target_xs[z].end());
    std::sort(target_ys[z].begin(), target_ys[z].end());
    target_ys[z].erase(std::unique(target_ys[z].begin(), target_ys[z].end()),
                       target_ys[z].end());
  }
  auto nearest = [](const std::vector<GeoCoord>& v, const GeoCoord c) {
    auto it = std::lower_bound(v.begin(), v.end(), c);
    int64_t d = std::numeric_limits<int64_t>::max() / 8;
    if (it != v.end()) {
      d = int64_t{*it} - c;
    }
    if (it != v.begin()) {
      d = std::min(d, int64_t{c} - *(it - 1));
    }
    return d;
  };
  auto future = [&](const int z, const GridPoint& p) -> int64_t {
    if (!req.future_cost) {
      return 0;
    }
    int64_t h = kInf;
    for (int t = 0; t < num_layers; t++) {
      if (target_xs[t].empty()) {
        continue;
      }
      h = std::min(
          h,
          (nearest(target_xs[t], p.x) + nearest(target_ys[t], p.y)) * unit_cost
              + req.via_cost * std::abs(t - z));
    }
    return h == kInf ? 0 : h;
  };
  struct Item
  {
    int64_t f;
    int64_t g;
    int64_t state;
    bool operator>(const Item& o) const { return f > o.f; }
  };
  std::priority_queue<Item, std::vector<Item>, std::greater<>> heap;
  for (const OffTrackSource& s : req.sources) {
    if (s.layer < 0 || s.layer >= num_layers) {
      continue;
    }
    const LayerIndex& l = index_[s.layer];
    const int x0 = static_cast<int>(
        std::lower_bound(l.xs.begin(), l.xs.end(), s.region.xlo)
        - l.xs.begin());
    const int y0 = static_cast<int>(
        std::lower_bound(l.ys.begin(), l.ys.end(), s.region.ylo)
        - l.ys.begin());
    for (int ix = x0; ix < l.nx() && l.xs[ix] <= s.region.xhi; ix++) {
      for (int iy = y0; iy < l.ny() && l.ys[iy] <= s.region.yhi; iy++) {
        const size_t v = l.base + static_cast<size_t>(ix) * l.ny() + iy;
        const int64_t st = static_cast<int64_t>(v) * kStates + kFresh;
        if (dist(st) > 0) {
          set(st, 0, -1);
          heap.push({future(s.layer, {l.xs[ix], l.ys[iy]}), 0, st});
        }
      }
    }
  }
  // Target vertices -> target ids.
  std::vector<std::pair<size_t, int>> targets;
  for (const OffTrackTarget& t : req.targets) {
    if (t.layer < 0 || t.layer >= num_layers) {
      continue;
    }
    const LayerIndex& l = index_[t.layer];
    const int ix = exactIndex(l.xs, t.p.x);
    const int iy = exactIndex(l.ys, t.p.y);
    if (ix < 0 || iy < 0) {
      continue;
    }
    targets.emplace_back(l.base + static_cast<size_t>(ix) * l.ny() + iy, t.id);
  }
  std::sort(targets.begin(), targets.end());
  std::vector<int64_t> reached(targets.size(), -1);
  size_t distinct = 0;
  for (size_t k = 0; k < targets.size(); k++) {
    distinct += (k == 0 || targets[k].first != targets[k - 1].first);
  }
  size_t remaining = distinct;

  auto relax = [&](const int64_t from,
                   const int64_t to,
                   const int64_t cost,
                   const int to_layer,
                   const GridPoint& to_point) {
    res.arcs++;
    if (cost < dist(to)) {
      set(to, cost, from);
      const size_t v = static_cast<size_t>(to / kStates);
      if (future_stamp_[v] != epoch_) {
        future_stamp_[v] = epoch_;
        future_memo_[v] = future(to_layer, to_point);
      }
      heap.push({cost + future_memo_[v], cost, to});
    }
  };
  auto via_ok
      = [&](const int lower_z, const size_t lower_v, const GridPoint p) {
          if (via_stamp_[lower_v] != epoch_) {
            via_stamp_[lower_v] = epoch_;
            via_memo_[lower_v]
                = (req.via_allowed && req.via_allowed(lower_z, p)) ? 1 : 0;
          }
          return via_memo_[lower_v] == 1;
        };

  while (!heap.empty() && remaining > 0) {
    const Item top = heap.top();
    heap.pop();
    const int64_t cost = top.g;
    const int64_t st = top.state;
    if (cost > dist(st)) {
      continue;
    }
    if (cost > req.max_cost) {
      break;
    }
    res.settled++;
    const size_t v = static_cast<size_t>(st / kStates);
    const int d = static_cast<int>(st % kStates);
    const int z = layer_of(v);
    LayerIndex& l = index_[z];
    const OffTrackLayer& in = req.layers[z];
    const int ix = static_cast<int>((v - l.base) / l.ny());
    const int iy = static_cast<int>((v - l.base) % l.ny());
    const GridPoint p{l.xs[ix], l.ys[iy]};
    {
      auto it = std::lower_bound(
          targets.begin(), targets.end(), std::make_pair(v, -1));
      if (it != targets.end() && it->first == v) {
        const size_t k = it - targets.begin();
        if (reached[k] < 0) {
          for (size_t q = k; q < targets.size() && targets[q].first == v; q++) {
            reached[q] = st;
          }
          remaining--;
          if (req.max_targets > 0
              && distinct - remaining >= static_cast<size_t>(req.max_targets)) {
            remaining = 0;
          }
        }
      }
    }
    auto state = [&](const int jx, const int jy, const int dir) {
      return static_cast<int64_t>(l.base + static_cast<size_t>(jx) * l.ny()
                                  + jy)
                 * kStates
             + dir;
    };
    auto planar = [&](const int jx, const int jy, const int dir) {
      const GridPoint q{l.xs[jx], l.ys[jy]};
      if (req.arc_filter && !req.arc_filter(z, p, q, d, dir)) {
        return;
      }
      const int64_t c = std::abs(int64_t{q.x} - p.x) * in.cost_x
                        + std::abs(int64_t{q.y} - p.y) * in.cost_y;
      relax(st, state(jx, jy, dir), cost + c, z, q);
    };
    if (in.wires) {
      // Straight arc to the neighbour in the incoming direction.
      if (d == 0 && ix + 1 < l.nx() && freeH(l, iy, ix, ix + 1)) {
        planar(ix + 1, iy, 0);
      } else if (d == 1 && ix > 0 && freeH(l, iy, ix - 1, ix)) {
        planar(ix - 1, iy, 1);
      } else if (d == 2 && iy + 1 < l.ny() && freeV(l, ix, iy, iy + 1)) {
        planar(ix, iy + 1, 2);
      } else if (d == 3 && iy > 0 && freeV(l, ix, iy - 1, iy)) {
        planar(ix, iy - 1, 3);
      }
      // Turn arcs to the nearest vertices at distance >= tau.
      const GeoCoord tau = in.grid.tau;
      const bool turn_x = d == kFresh || d == 2 || d == 3;
      const bool turn_y = d == kFresh || d == 0 || d == 1;
      if (turn_x) {
        const int kp = static_cast<int>(
            std::lower_bound(l.xs.begin(), l.xs.end(), p.x + tau)
            - l.xs.begin());
        if (kp < l.nx() && freeH(l, iy, ix, kp)) {
          planar(kp, iy, 0);
        }
        const int km = static_cast<int>(
                           std::upper_bound(l.xs.begin(), l.xs.end(), p.x - tau)
                           - l.xs.begin())
                       - 1;
        if (km >= 0 && freeH(l, iy, km, ix)) {
          planar(km, iy, 1);
        }
      }
      if (turn_y) {
        const int kp = static_cast<int>(
            std::lower_bound(l.ys.begin(), l.ys.end(), p.y + tau)
            - l.ys.begin());
        if (kp < l.ny() && freeV(l, ix, iy, kp)) {
          planar(ix, kp, 2);
        }
        const int km = static_cast<int>(
                           std::upper_bound(l.ys.begin(), l.ys.end(), p.y - tau)
                           - l.ys.begin())
                       - 1;
        if (km >= 0 && freeV(l, ix, km, iy)) {
          planar(ix, km, 3);
        }
      }
    }
    // Via arcs to the "no incoming direction" state of the adjacent layer.
    for (const int nz : {z + 1, z - 1}) {
      if (nz < 0 || nz >= num_layers) {
        continue;
      }
      const LayerIndex& n = index_[nz];
      const int jx = nz > z ? l.up_x[ix] : l.down_x[ix];
      const int jy = nz > z ? l.up_y[iy] : l.down_y[iy];
      if (jx < 0 || jy < 0) {
        continue;
      }
      const size_t w = n.base + static_cast<size_t>(jx) * n.ny() + jy;
      const bool ok = nz > z ? via_ok(z, v, p) : via_ok(nz, w, p);
      if (ok) {
        relax(st,
              static_cast<int64_t>(w) * kStates + kFresh,
              cost + req.via_cost,
              nz,
              p);
      }
    }
  }

  // Paths to the reached targets.
  for (size_t k = 0; k < targets.size(); k++) {
    if (reached[k] < 0) {
      continue;
    }
    OffTrackPath path;
    path.target_id = targets[k].second;
    path.cost = dist_[reached[k]];
    std::vector<RouteNode> full;
    for (int64_t s = reached[k]; s >= 0; s = parent_[s]) {
      const size_t v = static_cast<size_t>(s / kStates);
      const int z = layer_of(v);
      const LayerIndex& l = index_[z];
      const int ix = static_cast<int>((v - l.base) / l.ny());
      const int iy = static_cast<int>((v - l.base) % l.ny());
      full.push_back({z, {l.xs[ix], l.ys[iy]}});
    }
    std::reverse(full.begin(), full.end());
    const size_t n = full.size();
    for (size_t q = 0; q < n; q++) {
      bool keep = q == 0 || q + 1 == n;
      if (!keep) {
        const RouteNode& a = full[q - 1];
        const RouteNode& b = full[q];
        const RouteNode& c = full[q + 1];
        keep = a.layer != b.layer || c.layer != b.layer
               || !((a.p.x == b.p.x && b.p.x == c.p.x)
                    || (a.p.y == b.p.y && b.p.y == c.p.y));
      }
      if (keep) {
        path.nodes.push_back(full[q]);
      }
    }
    res.paths.push_back(std::move(path));
  }
  std::stable_sort(res.paths.begin(),
                   res.paths.end(),
                   [](const OffTrackPath& a, const OffTrackPath& b) {
                     return a.cost < b.cost;
                   });
  return res;
}

}  // namespace drt::br
