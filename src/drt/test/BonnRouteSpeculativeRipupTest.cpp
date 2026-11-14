// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include <cstdint>
#include <cstdio>
#include <memory>
#include <vector>

#include "br/BrCore.h"
#include "br/BrSpeculativeRipup.h"

namespace {
using namespace drt::br;

struct SearchSnapshot
{
  SearchResult result;
  TileBox reads;
};

TileBox written(const AccessFootprint& footprint)
{
  TileBox box;
  for (int step = 0; step < AccessFootprint::kMaxSteps; step++) {
    box.add(footprint.writes(step));
  }
  return box;
}

struct Fixture
{
  TrackGraph graph;
  FastGrid grid;
  Region region;
  SearchRequest request;
  std::unique_ptr<PathSearch> vertex;
  std::unique_ptr<IntervalSearch> interval;
  bool intervals;

  Fixture(const bool dense, const bool use_intervals) : intervals(use_intervals)
  {
    graph.layers.resize(1);
    graph.layers[0].tracks = {0, 10, 20};
    for (int j = 0; j <= 32; j++) {
      graph.layers[0].positions.push_back(j * 10);
    }
    graph.finalize([](Coord x) { return x / 40; },
                   [](Coord y) { return y / 40; });
    grid.init(&graph, 1);
    grid.setDense(0, dense);
    region.init(1, graph.tiles_x, graph.tiles_y);
    // Leave distant tiles outside the search region for the reuse case.
    for (int x = 0; x <= 3; x++) {
      region.allow(0, x, 0);
    }
    request.region = &region;
    request.sources = {{graph.vertex(0, 1, 0), 0}};
    request.targets = {{graph.vertex(0, 1, 12), 0, 0}};
    request.target_boxes = {{120, 10, 120, 10, 0, 0}};
    vertex = std::make_unique<PathSearch>(graph, grid);
    interval = std::make_unique<IntervalSearch>(graph, grid);
  }

  SearchSnapshot search()
  {
    AccessFootprint footprint;
    footprint.begin(0, {0, 0, graph.tiles_x - 1, graph.tiles_y - 1});
    const AccessFootprintScope scope(&footprint);
    SearchSnapshot snapshot;
    snapshot.result = intervals ? interval->run(request) : vertex->run(request);
    snapshot.reads = footprint.accessed();
    snapshot.reads.add(footprint.wideReads());
    return snapshot;
  }

  TileBox edge(const int track, const int position, const RipupLevel level)
  {
    AccessFootprint footprint;
    footprint.begin(0, {0, 0, graph.tiles_x - 1, graph.tiles_y - 1});
    const AccessFootprintScope scope(&footprint);
    grid.setLevel(graph.vertex(0, track, position), kWire, 0, level);
    return written(footprint);
  }
};

bool cachedReads(const bool dense, const bool overlaid)
{
  TrackGraph graph;
  graph.layers.resize(2);
  graph.layers[0].tracks = {0, 10, 20};
  graph.layers[1].horizontal = false;
  graph.layers[1].positions = {0, 10, 20};
  for (int j = 0; j <= 32; j++) {
    graph.layers[0].positions.push_back(j * 10);
    graph.layers[1].tracks.push_back(j * 10);
  }
  graph.finalize([](Coord x) { return x / 40; },
                 [](Coord y) { return y / 40; });
  FastGrid grid;
  grid.init(&graph, 1);
  grid.setDense(0, dense);
  FastGridOverlay overlay;
  if (overlaid) {
    const auto key = FastGridOverlay::key(0, 1);
    overlay.wiring[key] = {{0, 32, 1}};
    overlay.via[key] = {{0, 32, 1}};
  }
  FastGridReader reader;
  reader.reset(&grid, overlaid ? &overlay : nullptr);
  const uint64_t expected = overlaid ? 1 : 0;
  // This first load warms a run or a whole dense row without recording.
  if (reader.wiringWord(0, 1, 1) != expected
      || reader.viaWord(0, 1, grid.viaIndex(0, 1)) != expected) {
    return false;
  }
  for (const int position : {2, 31, 30, 1}) {
    for (const bool via : {false, true}) {
      AccessFootprint footprint;
      footprint.begin(0, {0, 0, graph.tiles_x - 1, graph.tiles_y - 1});
      const AccessFootprintScope scope(&footprint);
      const uint64_t word
          = via ? reader.viaWord(0, 1, grid.viaIndex(0, position))
                : reader.wiringWord(0, 1, position);
      const TileBox tile{position / 4, 0, position / 4, 0};
      if (word != expected || footprint.reads(0) != tile) {
        return false;
      }
    }
  }
  return true;
}

bool removal(const bool dense, const bool intervals)
{
  Fixture f(dense, intervals);
  f.edge(1, 5, kNotRemovable);
  const auto before = f.search();
  const TileBox writes = f.edge(1, 5, kLegal);
  const auto after = f.search();
  // Removing the obstacle cannot invalidate the old detour's legality,
  // yet its cost is no longer optimal. Path legality alone cannot reuse it.
  return before.result.found && after.result.found
         && before.result.cost > after.result.cost
         && !ripupSearchReusable(before.reads, writes, true);
}

bool failedSearch(const bool dense, const bool intervals)
{
  Fixture f(dense, intervals);
  for (int i = 0; i < 3; i++) {
    f.edge(i, 5, kNotRemovable);
  }
  const auto before = f.search();
  const TileBox writes = f.edge(1, 5, kLegal);
  const auto after = f.search();
  return !before.result.found && after.result.found
         && !ripupSearchReusable(before.reads, writes, true);
}

bool distantAndRollback(const bool dense, const bool intervals)
{
  Fixture f(dense, intervals);
  const auto before = f.search();
  const TileBox distant = f.edge(1, 31, kNotRemovable);
  const auto after = f.search();
  if (!before.result.found || before.result.path != after.result.path
      || before.result.cost != after.result.cost
      || !ripupSearchReusable(before.reads, distant, true)) {
    return false;
  }
  TileBox rollback = f.edge(1, 5, kNotRemovable);
  rollback.add(f.edge(1, 5, kLegal));
  const auto restored = f.search();
  // A rollback still records its writes; a changed root revision also
  // invalidates presearch when no shared geometry changed.
  return before.result.path == restored.result.path
         && !ripupSearchReusable(before.reads, rollback, true)
         && !ripupSearchReusable(before.reads, {}, false);
}

bool history(const bool dense, const bool intervals)
{
  Fixture f(dense, intervals);
  f.request.ripup = true;
  f.request.ripup_cost = 10;
  const size_t conflict = f.graph.vertex(0, 1, 6);
  f.grid.setLevel(conflict, kWire, 0, 1);
  const auto before = f.search();
  AccessFootprint footprint;
  footprint.begin(0, {0, 0, f.graph.tiles_x - 1, f.graph.tiles_y - 1});
  {
    const AccessFootprintScope scope(&footprint);
    for (int i = 0; i < 20; i++) {
      f.grid.bumpHistory(conflict);
    }
  }
  const auto after = f.search();
  return before.result.found && after.result.found
         && before.result.cost < after.result.cost
         && before.result.path != after.result.path
         && !ripupSearchReusable(before.reads, written(footprint), true);
}
}  // namespace

int main()
{
  bool passed = true;
  for (const bool dense : {false, true}) {
    for (const bool overlaid : {false, true}) {
      if (!cachedReads(dense, overlaid)) {
        std::printf(
            "cached footprint failed: dense=%d overlay=%d\n", dense, overlaid);
        passed = false;
      }
    }
    for (const bool intervals : {false, true}) {
      const bool removed = removal(dense, intervals);
      const bool failed = failedSearch(dense, intervals);
      const bool rollback = distantAndRollback(dense, intervals);
      const bool bumped = history(dense, intervals);
      if (!removed || !failed || !rollback || !bumped) {
        std::printf(
            "speculative invalidation failed: dense=%d intervals=%d "
            "removal=%d failure=%d rollback=%d history=%d\n",
            dense,
            intervals,
            removed,
            failed,
            rollback,
            bumped);
        passed = false;
      }
    }
  }
  return passed ? 0 : 1;
}
