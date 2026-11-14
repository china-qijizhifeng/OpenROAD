// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrWireJunctions.h"

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <tuple>
#include <vector>

namespace drt::br {
std::vector<std::vector<GridPoint>> interiorWireCrossings(
    const std::span<const WireSegment> wires)
{
  enum Kind
  {
    kStart,
    kVertical,
    kEnd
  };
  struct Event
  {
    int layer;
    GeoCoord x;
    Kind kind;
    size_t wire;
  };
  std::vector<Event> events;
  for (size_t k = 0; k < wires.size(); k++) {
    const WireSegment& w = wires[k];
    if (w.begin == w.end) {
      continue;
    }
    if (w.begin.y == w.end.y) {
      events.push_back({w.layer, std::min(w.begin.x, w.end.x), kStart, k});
      events.push_back({w.layer, std::max(w.begin.x, w.end.x), kEnd, k});
    } else if (w.begin.x == w.end.x) {
      events.push_back({w.layer, w.begin.x, kVertical, k});
    }
  }
  std::ranges::sort(events, [](const Event& a, const Event& b) {
    return std::tie(a.layer, a.x, a.kind, a.wire)
           < std::tie(b.layer, b.x, b.kind, b.wire);
  });
  std::vector<std::vector<GridPoint>> cuts(wires.size());
  // Sweep x; only horizontal wires spanning this x can cross a vertical.
  std::map<GeoCoord, std::set<size_t>> active;
  int layer = events.empty() ? 0 : events.front().layer;
  for (const Event& event : events) {
    if (event.layer != layer) {
      active.clear();
      layer = event.layer;
    }
    const WireSegment& w = wires[event.wire];
    if (event.kind == kStart) {
      active[w.begin.y].insert(event.wire);
    } else if (event.kind == kEnd) {
      auto row = active.find(w.begin.y);
      row->second.erase(event.wire);
      if (row->second.empty()) {
        active.erase(row);
      }
    } else {
      const GeoCoord low = std::min(w.begin.y, w.end.y);
      const GeoCoord high = std::max(w.begin.y, w.end.y);
      for (auto row = active.upper_bound(low);
           row != active.end() && row->first < high;
           row++) {
        for (const size_t k : row->second) {
          const WireSegment& h = wires[k];
          if (event.x <= std::min(h.begin.x, h.end.x)
              || event.x >= std::max(h.begin.x, h.end.x)) {
            continue;
          }
          const GridPoint cross{event.x, row->first};
          cuts[k].push_back(cross);
          cuts[event.wire].push_back(cross);
        }
      }
    }
  }
  for (auto& wire_cuts : cuts) {
    std::ranges::sort(wire_cuts);
    wire_cuts.erase(std::ranges::unique(wire_cuts).begin(), wire_cuts.end());
  }
  return cuts;
}
}  // namespace drt::br
