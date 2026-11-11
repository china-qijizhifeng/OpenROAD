// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Same-net rules (Gester et al., TODAES 18(2), 2013, Sec. 3.7) and their
// mapping to minimum segment lengths tau (Sec. 3.8: "Most same-net rules can
// be mapped to requirements on minimum segment lengths in the constructed
// path [Nieberg 2011]"). Nieberg 2011 was not available; the mapping below
// is derived from the boundary polygon of a centre-line path of uniform
// width and is documented in OffTrack.md.

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include "br/BrGeometry.h"

namespace drt::br {

// The same-net rules of one wiring layer that the off-track search and the
// post-processing take into account.
struct SameNetLayerRules
{
  // Width of the wire model (both directions; Sec. 3.2 default wire type).
  GeoCoord width = 0;
  // Short edge rule (LEF MINSTEP): of two adjacent boundary edges at least
  // one must be this long; 0 = no rule.
  GeoCoord min_step_length = 0;
  // Notch rule: distance between non-adjacent parts of the same metal; 0 =
  // no rule. Taken from the layer's minimum spacing.
  GeoCoord notch_spacing = 0;
  // Minimum area rule (LEF AREA) in DBU^2; 0 = no rule.
  int64_t min_area = 0;
};

struct MinimumSegmentLength
{
  GeoCoord tau = 1;
  // The individual lower bounds (0 where the rule is absent or not mapped).
  GeoCoord from_short_edge = 0;
  GeoCoord from_notch = 0;
};

// Mapping of the rules to one minimum segment length tau per layer.
//
// A centre-line segment of length L between two bends of a wire of width w
// carries boundary edges of length L - w, L or L + w (L - w on the inside of
// a U, L on both sides of a Z). Requiring every segment to have L - w >=
// min_step_length makes every boundary edge except the line-end edges (of
// length w) at least min_step_length long, and line-end edges are adjacent
// only to such edges, so the short edge rule holds: tau_step =
// min_step_length + w. The two legs of a U whose base has length L are L - w
// apart; tau_notch = notch_spacing + w. Notch avoidance through tau is
// optional (`include_notch`) because it also forbids short Z-jogs, which
// have no notch; without it, U-notches are filtered by `pathViolations`
// and repaired by the post-processing.
//
// Minimum area is not mapped to tau: it constrains a whole connected
// polygon, not each segment, and the off-track paths are attached to a pin
// or continued by on-track wiring. The post-processing handles it.
//
// Without any mapped rule tau is the wire width (the smallest segment that
// does not create a boundary edge shorter than the line end).
inline MinimumSegmentLength minimumSegmentLength(const SameNetLayerRules& r,
                                                 const bool include_notch,
                                                 const GeoCoord override_tau)
{
  MinimumSegmentLength m;
  if (r.min_step_length > 0) {
    m.from_short_edge = r.min_step_length + r.width;
  }
  if (r.notch_spacing > 0) {
    m.from_notch = r.notch_spacing + r.width;
  }
  m.tau = std::max<GeoCoord>(m.from_short_edge, r.width);
  if (include_notch) {
    m.tau = std::max(m.tau, m.from_notch);
  }
  if (override_tau > 0) {
    m.tau = override_tau;
  }
  m.tau = std::max<GeoCoord>(m.tau, 1);
  return m;
}

// Exact area of a union of boxes (coordinate compression).
inline int64_t unionArea(const std::vector<Box>& boxes)
{
  std::vector<GeoCoord> xs;
  std::vector<GeoCoord> ys;
  for (const Box& b : boxes) {
    if (b.inverted()) {
      continue;
    }
    xs.push_back(b.xlo);
    xs.push_back(b.xhi);
    ys.push_back(b.ylo);
    ys.push_back(b.yhi);
  }
  std::sort(xs.begin(), xs.end());
  xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
  std::sort(ys.begin(), ys.end());
  ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
  int64_t area = 0;
  for (size_t i = 0; i + 1 < xs.size(); i++) {
    for (size_t j = 0; j + 1 < ys.size(); j++) {
      const GeoCoord cx = xs[i];
      const GeoCoord cy = ys[j];
      for (const Box& b : boxes) {
        if (!b.inverted() && cx >= b.xlo && xs[i + 1] <= b.xhi && cy >= b.ylo
            && ys[j + 1] <= b.yhi) {
          area += (int64_t{xs[i + 1]} - cx) * (int64_t{ys[j + 1]} - cy);
          break;
        }
      }
    }
  }
  return area;
}

// Box of a centre-line wire piece with square ends (half-width extension).
inline Box wireBox(const GridPoint& a, const GridPoint& b, GeoCoord half_width)
{
  return {std::min(a.x, b.x) - half_width,
          std::min(a.y, b.y) - half_width,
          std::max(a.x, b.x) + half_width,
          std::max(a.y, b.y) + half_width};
}

// Same-net violations of a planar centre-line path (one layer, consecutive
// points form axis-parallel segments).
struct SameNetViolations
{
  // Segments shorter than tau between two bends.
  int short_inner = 0;
  // First or last segments shorter than tau (their outer end is a path
  // end: pin, via or attachment).
  int short_end = 0;
  // U-shapes whose legs are closer than the notch spacing.
  int notches = 0;
  int total() const { return short_inner + short_end + notches; }
};

// Points with zero-length and collinear interior points removed (reversals
// included: the dropped point lies on the remaining metal or only adds a
// spike).
inline std::vector<GridPoint> simplifyPolyline(const std::vector<GridPoint>& in)
{
  std::vector<GridPoint> out = in;
  bool changed = true;
  while (changed) {
    changed = false;
    std::vector<GridPoint> next;
    for (const GridPoint& p : out) {
      if (!next.empty() && next.back() == p) {
        changed = true;
        continue;
      }
      if (next.size() >= 2) {
        const GridPoint& a = next[next.size() - 2];
        const GridPoint& b = next.back();
        if ((a.x == b.x && b.x == p.x) || (a.y == b.y && b.y == p.y)) {
          next.back() = p;
          changed = true;
          continue;
        }
      }
      next.push_back(p);
    }
    out.swap(next);
  }
  return out;
}

inline SameNetViolations pathViolations(const std::vector<GridPoint>& raw,
                                        const GeoCoord tau,
                                        const SameNetLayerRules& rules)
{
  SameNetViolations v;
  const std::vector<GridPoint> p = simplifyPolyline(raw);
  const int n = static_cast<int>(p.size()) - 1;  // segments
  auto len = [&](int k) {
    return std::abs(p[k + 1].x - p[k].x) + std::abs(p[k + 1].y - p[k].y);
  };
  for (int k = 0; k < n; k++) {
    if (len(k) >= tau) {
      continue;
    }
    if (k == 0 || k == n - 1) {
      v.short_end++;
    } else {
      v.short_inner++;
    }
  }
  if (rules.notch_spacing > 0) {
    for (int k = 1; k + 1 < n; k++) {
      // Legs k - 1 and k + 1 opposite: a U with base k.
      const int ax = p[k].x - p[k - 1].x;
      const int ay = p[k].y - p[k - 1].y;
      const int bx = p[k + 2].x - p[k + 1].x;
      const int by = p[k + 2].y - p[k + 1].y;
      const bool opposite = (ax != 0 && bx != 0 && (ax > 0) != (bx > 0))
                            || (ay != 0 && by != 0 && (ay > 0) != (by > 0));
      if (opposite && len(k) - rules.width < rules.notch_spacing) {
        v.notches++;
      }
    }
  }
  return v;
}

}  // namespace drt::br
