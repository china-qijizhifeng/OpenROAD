// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrPathCleanup.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include "br/BrSameNetRules.h"

namespace drt::br {

namespace {

struct RunPoint
{
  GridPoint p;
  bool fixed = false;
};

bool collinear(const GridPoint& a, const GridPoint& b, const GridPoint& c)
{
  return (a.x == b.x && b.x == c.x) || (a.y == b.y && b.y == c.y);
}

void normalize(std::vector<RunPoint>& r)
{
  bool changed = true;
  while (changed) {
    changed = false;
    std::vector<RunPoint> out;
    for (const RunPoint& q : r) {
      if (!out.empty() && out.back().p == q.p) {
        out.back().fixed = out.back().fixed || q.fixed;
        changed = true;
        continue;
      }
      out.push_back(q);
    }
    for (size_t k = 1; k + 1 < out.size(); k++) {
      if (!out[k].fixed && collinear(out[k - 1].p, out[k].p, out[k + 1].p)) {
        out.erase(out.begin() + static_cast<long>(k));
        changed = true;
        break;
      }
    }
    r.swap(out);
  }
}

std::vector<GridPoint> points(const std::vector<RunPoint>& r)
{
  std::vector<GridPoint> p;
  p.reserve(r.size());
  for (const RunPoint& q : r) {
    p.push_back(q.p);
  }
  return p;
}

struct Score
{
  int inner = 0;  // short inner segments + notches
  int end = 0;
  bool betterThan(const Score& o) const
  {
    return inner < o.inner && end <= o.end;
  }
};

Score score(const std::vector<RunPoint>& r, const CleanupLayer& layer)
{
  const SameNetViolations v = pathViolations(points(r), layer.tau, layer.rules);
  return {v.short_inner + v.notches, v.short_end};
}

int64_t segLength(const GridPoint& a, const GridPoint& b)
{
  return std::abs(int64_t{a.x} - b.x) + std::abs(int64_t{a.y} - b.y);
}

GridPoint plus(const GridPoint& a, const GridPoint& v)
{
  return {a.x + v.x, a.y + v.y};
}
GridPoint minus(const GridPoint& a, const GridPoint& v)
{
  return {a.x - v.x, a.y - v.y};
}

class RunCleaner
{
 public:
  RunCleaner(const CleanupRequest& req, const int layer, CleanupStats& stats)
      : req_(req), layer_(layer), rules_(req.layers[layer]), stats_(stats)
  {
  }

  void clean(std::vector<RunPoint>& r)
  {
    normalize(r);
    stats_.violations_before += score(r, rules_).inner + score(r, rules_).end;
    const int max_changes = 2 * static_cast<int>(r.size()) + 8;
    int changes = 0;
    bool progress = true;
    while (progress && changes < max_changes) {
      progress = false;
      const int m = static_cast<int>(r.size()) - 1;
      for (int k = 0; k < m && !progress; k++) {
        if (segLength(r[k].p, r[k + 1].p) >= rules_.tau) {
          continue;
        }
        progress = tryFixes(r, k);
      }
      changes += progress;
    }
    const Score s = score(r, rules_);
    stats_.violations_after += s.inner + s.end;
  }

 private:
  bool legal(const GridPoint& a, const GridPoint& b) const
  {
    return !req_.allowed
           || req_.allowed(layer_, wireBox(a, b, rules_.half_width));
  }

  bool accept(std::vector<RunPoint>& r, std::vector<RunPoint>& cand)
  {
    normalize(cand);
    if (score(cand, rules_).betterThan(score(r, rules_))) {
      r.swap(cand);
      return true;
    }
    return false;
  }

  bool tryFixes(std::vector<RunPoint>& r, const int k)
  {
    const int m = static_cast<int>(r.size()) - 1;
    const GridPoint sv{r[k + 1].p.x - r[k].p.x, r[k + 1].p.y - r[k].p.y};
    // Fix A: move a = (P_k-1, P_k) by sv; P_k-1 must be interior.
    if (k >= 2 && !r[k].fixed && !r[k - 1].fixed) {
      const GridPoint moved = plus(r[k - 1].p, sv);
      if (legal(r[k - 1].p, moved) && legal(moved, r[k + 1].p)) {
        std::vector<RunPoint> cand = r;
        cand[k - 1].p = moved;
        cand.erase(cand.begin() + k);
        if (accept(r, cand)) {
          stats_.moved_segments++;
          return true;
        }
      }
    }
    // Fix B: move b = (P_k+1, P_k+2) by -sv; P_k+2 must be interior.
    if (k + 2 <= m - 1 && !r[k + 1].fixed && !r[k + 2].fixed) {
      const GridPoint moved = minus(r[k + 2].p, sv);
      if (legal(r[k].p, moved) && legal(moved, r[k + 2].p)) {
        std::vector<RunPoint> cand = r;
        cand[k + 2].p = moved;
        cand.insert(cand.begin() + k + 3, RunPoint{r[k + 2].p, false});
        cand.erase(cand.begin() + k + 1);
        if (accept(r, cand)) {
          stats_.moved_segments++;
          return true;
        }
      }
    }
    // Fix C: a U with short base s; collapse the shorter leg.
    if (k >= 1 && k + 2 <= m && !r[k].fixed && !r[k + 1].fixed) {
      const GridPoint av{r[k].p.x - r[k - 1].p.x, r[k].p.y - r[k - 1].p.y};
      const GridPoint bv{r[k + 2].p.x - r[k + 1].p.x,
                         r[k + 2].p.y - r[k + 1].p.y};
      const bool opposite
          = (av.x != 0 && bv.x != 0 && (av.x > 0) != (bv.x > 0))
            || (av.y != 0 && bv.y != 0 && (av.y > 0) != (bv.y > 0));
      if (opposite) {
        const int64_t la = segLength(r[k - 1].p, r[k].p);
        const int64_t lb = segLength(r[k + 1].p, r[k + 2].p);
        std::vector<RunPoint> cand = r;
        bool ok = false;
        if (la <= lb) {
          const GridPoint q = plus(r[k - 1].p, sv);
          if (legal(r[k - 1].p, q)) {
            cand[k].p = q;
            cand.erase(cand.begin() + k + 1);
            ok = true;
          }
        } else {
          const GridPoint q = plus(r[k].p, bv);
          if (legal(q, r[k + 2].p)) {
            cand[k + 1].p = q;
            cand.erase(cand.begin() + k);
            ok = true;
          }
        }
        if (ok && accept(r, cand)) {
          stats_.collapsed_notches++;
          return true;
        }
      }
    }
    return false;
  }

  const CleanupRequest& req_;
  const int layer_;
  const CleanupLayer& rules_;
  CleanupStats& stats_;
};

void fixMinArea(const CleanupRequest& req,
                const int layer,
                const std::vector<RunPoint>& r,
                const bool start_up,
                const bool end_up,
                std::vector<CleanupPatch>& patches,
                CleanupStats& stats)
{
  const CleanupLayer& l = req.layers[layer];
  if (l.rules.min_area <= 0 || r.empty()) {
    return;
  }
  std::vector<Box> boxes;
  for (size_t k = 0; k + 1 < r.size(); k++) {
    boxes.push_back(wireBox(r[k].p, r[k + 1].p, l.half_width));
  }
  if (req.via_pads) {
    req.via_pads(layer, r.front().p, start_up, boxes);
    req.via_pads(layer, r.back().p, end_up, boxes);
  }
  if (unionArea(boxes) >= l.rules.min_area) {
    return;
  }
  stats.area_violations++;
  const GeoCoord width = std::max<GeoCoord>(1, 2 * l.half_width);
  // Candidate extensions: beyond the last node along the last segment,
  // beyond the first node along the first segment, then both preferred
  // directions.
  std::vector<std::pair<GridPoint, GridPoint>> dirs;  // (anchor, unit)
  auto unit = [](const GridPoint& from, const GridPoint& to) {
    return GridPoint{(to.x > from.x) - (to.x < from.x),
                     (to.y > from.y) - (to.y < from.y)};
  };
  if (r.size() >= 2) {
    dirs.emplace_back(r.back().p, unit(r[r.size() - 2].p, r.back().p));
    dirs.emplace_back(r.front().p, unit(r[1].p, r.front().p));
  } else {
    const GridPoint pref = l.horizontal ? GridPoint{1, 0} : GridPoint{0, 1};
    dirs.emplace_back(r.front().p, pref);
    dirs.emplace_back(r.front().p, GridPoint{-pref.x, -pref.y});
  }
  for (const auto& [anchor, u] : dirs) {
    int64_t need = (l.rules.min_area - unionArea(boxes) + width - 1) / width;
    for (int attempt = 0; attempt < 4 && need > 0; attempt++) {
      const GridPoint tip{anchor.x + static_cast<GeoCoord>(u.x * need),
                          anchor.y + static_cast<GeoCoord>(u.y * need)};
      const Box ext = wireBox(anchor, tip, l.half_width);
      if (req.allowed && !req.allowed(layer, ext)) {
        break;
      }
      std::vector<Box> with = boxes;
      with.push_back(ext);
      const int64_t area = unionArea(with);
      if (area >= l.rules.min_area) {
        patches.push_back({layer, anchor, tip});
        stats.area_fixed++;
        return;
      }
      need += (l.rules.min_area - area + width - 1) / width;
    }
  }
}

}  // namespace

CleanupResult cleanupConnection(const CleanupRequest& req)
{
  CleanupResult out;
  const std::vector<RouteNode>& nodes = req.nodes;
  const size_t n = nodes.size();
  size_t s = 0;
  while (s < n) {
    size_t e = s;
    while (e + 1 < n && nodes[e + 1].layer == nodes[s].layer) {
      e++;
    }
    const int layer = nodes[s].layer;
    std::vector<RunPoint> run;
    for (size_t k = s; k <= e; k++) {
      // Run ends are vias or chain ends; both stay.
      run.push_back({nodes[k].p, k == s || k == e});
    }
    out.stats.runs++;
    if (layer >= 0 && layer < static_cast<int>(req.layers.size())) {
      RunCleaner cleaner(req, layer, out.stats);
      cleaner.clean(run);
      // Via-to-via piece (neither end is a chain end).
      if (s > 0 && e + 1 < n) {
        fixMinArea(req,
                   layer,
                   run,
                   nodes[s - 1].layer > layer,
                   nodes[e + 1].layer > layer,
                   out.patches,
                   out.stats);
      }
    }
    for (const RunPoint& q : run) {
      out.nodes.push_back({layer, q.p});
    }
    s = e + 1;
  }
  return out;
}

}  // namespace drt::br
