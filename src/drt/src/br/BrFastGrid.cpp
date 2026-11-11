// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrFastGrid.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "br/BrDistanceRules.h"
#include "br/BrGeometry.h"
#include "br/BrShapeGrid.h"
#include "br/BrTrackGraph.h"
#include "br/BrWireModel.h"

namespace drt::br {

namespace {

int wiringKindIndex(const Kind kind)
{
  switch (kind) {
    case kWire:
      return 0;
    case kJog:
      return 1;
    case kViaUpPad:
      return 2;
    case kViaDownPad:
      return 3;
    default:
      return -1;
  }
}

}  // namespace

void FastGrid::init(const TrackGraph* graph, const int num_wire_types)
{
  graph_ = graph;
  num_wire_types_ = std::clamp(num_wire_types, 0, kMaxWireTypes);
  history_.clear();
  layers_.clear();
  layers_.resize(graph->numLayers());
  // Consecutive equal tile indices form one chunk.
  auto chunks = [](const std::vector<int>& tile,
                   const int n,
                   std::vector<int>& chunk_of,
                   std::vector<int>& first) {
    chunk_of.assign(n, 0);
    first.assign(1, 0);
    int chunk = 0;
    for (int k = 1; k < n; k++) {
      if (tile.size() == static_cast<size_t>(n) && tile[k] != tile[k - 1]) {
        chunk++;
        first.push_back(k);
      }
      chunk_of[k] = chunk;
    }
    first.push_back(n);
    return chunk + 1;
  };
  for (int z = 0; z < graph->numLayers(); z++) {
    const TrackLayer& l = graph->layers[z];
    Layer& fl = layers_[z];
    const int np = l.numPositions();
    fl.pos_chunks = chunks(l.pos_tile, np, fl.pos_chunk, fl.chunk_first_pos);
    fl.track_chunks = chunks(
        l.track_tile, l.numTracks(), fl.track_chunk, fl.chunk_first_track);
    const size_t tiles = static_cast<size_t>(fl.pos_chunks) * fl.track_chunks;
    fl.wiring.clear();
    fl.wiring.resize(tiles);
    fl.via_index.assign(np, -1);
    fl.via_position.clear();
    fl.via_chunk.clear();
    fl.via.clear();
    if (z + 1 < graph->numLayers()) {
      for (int j = 0; j < np; j++) {
        if (l.up_track[j] >= 0) {
          fl.via_index[j] = static_cast<int>(fl.via_position.size());
          fl.via_position.push_back(j);
          fl.via_chunk.push_back(fl.pos_chunk[j]);
        }
      }
      fl.via.resize(tiles);
    }
    fl.via_chunk_first.assign(fl.pos_chunks + 1,
                              static_cast<int>(fl.via_position.size()));
    for (int c = static_cast<int>(fl.via_chunk.size()) - 1; c >= 0; c--) {
      for (int k = fl.via_chunk[c]; k >= 0 && fl.via_chunk_first[k] > c; k--) {
        fl.via_chunk_first[k] = c;
      }
    }
  }
}

uint64_t FastGrid::lookup(const TileRuns* tile,
                          const int track,
                          const int index)
{
  if (tile == nullptr) {
    return 0;
  }
  const uint32_t a = tile->track_first[track];
  const uint32_t b = tile->track_first[track + 1];
  if (a == b) {
    return 0;
  }
  const auto begin = tile->starts.begin() + a;
  const auto it = std::upper_bound(begin, tile->starts.begin() + b, index);
  return tile->words[(it - tile->starts.begin()) - 1];
}

void FastGrid::assign(std::unique_ptr<TileRuns>& tile,
                      const int tracks,
                      const int track,
                      const int count,
                      const int lo,
                      const int hi,
                      const uint64_t word,
                      const uint64_t mask)
{
  if (!tile) {
    if ((word & mask) == 0) {
      return;
    }
    tile = std::make_unique<TileRuns>();
    tile->track_first.assign(tracks + 1, 0);
  }
  TileRuns& t = *tile;
  const uint32_t a = t.track_first[track];
  const uint32_t b = t.track_first[track + 1];
  thread_local std::vector<uint64_t> words;
  words.assign(count, 0);
  for (uint32_t r = a; r < b; r++) {
    const int end = r + 1 < b ? t.starts[r + 1] : count;
    std::fill(words.begin() + t.starts[r], words.begin() + end, t.words[r]);
  }
  for (int k = lo; k <= hi; k++) {
    words[k] = (words[k] & ~mask) | (word & mask);
  }
  thread_local std::vector<uint16_t> new_starts;
  thread_local std::vector<uint64_t> new_words;
  new_starts.clear();
  new_words.clear();
  for (int k = 0; k < count; k++) {
    if (k == 0 || words[k] != words[k - 1]) {
      new_starts.push_back(static_cast<uint16_t>(k));
      new_words.push_back(words[k]);
    }
  }
  if (new_words.size() == 1 && new_words.front() == 0) {
    new_starts.clear();
    new_words.clear();
  }
  const int64_t delta = static_cast<int64_t>(new_words.size()) - (b - a);
  t.starts.erase(t.starts.begin() + a, t.starts.begin() + b);
  t.starts.insert(t.starts.begin() + a, new_starts.begin(), new_starts.end());
  t.words.erase(t.words.begin() + a, t.words.begin() + b);
  t.words.insert(t.words.begin() + a, new_words.begin(), new_words.end());
  for (int k = track + 1; k <= tracks; k++) {
    t.track_first[k] = static_cast<uint32_t>(t.track_first[k] + delta);
  }
}

uint64_t FastGrid::wiringWord(const int z, const int i, const int j) const
{
  const Layer& fl = layers_[z];
  const int k = fl.pos_chunk[j];
  const int tc = fl.track_chunk[i];
  return lookup(fl.wiring[static_cast<size_t>(tc) * fl.pos_chunks + k].get(),
                i - fl.chunk_first_track[tc],
                j - fl.chunk_first_pos[k]);
}

uint64_t FastGrid::viaWord(const int z, const int i, const int c) const
{
  const Layer& fl = layers_[z];
  const int k = fl.via_chunk[c];
  const int tc = fl.track_chunk[i];
  return lookup(fl.via[static_cast<size_t>(tc) * fl.pos_chunks + k].get(),
                i - fl.chunk_first_track[tc],
                c - fl.via_chunk_first[k]);
}

void FastGrid::assignWiring(const int z,
                            const int i,
                            const int j0,
                            const int j1,
                            const uint64_t word,
                            const uint64_t mask)
{
  Layer& fl = layers_[z];
  const int tc = fl.track_chunk[i];
  const int tracks = fl.chunk_first_track[tc + 1] - fl.chunk_first_track[tc];
  for (int j = j0; j <= j1;) {
    const int k = fl.pos_chunk[j];
    const int first = fl.chunk_first_pos[k];
    const int last = fl.chunk_first_pos[k + 1] - 1;
    const int end = std::min(j1, last);
    assign(fl.wiring[static_cast<size_t>(tc) * fl.pos_chunks + k],
           tracks,
           i - fl.chunk_first_track[tc],
           last - first + 1,
           j - first,
           end - first,
           word,
           mask);
    j = end + 1;
  }
}

void FastGrid::assignVia(const int z,
                         const int i,
                         const int c0,
                         const int c1,
                         const uint64_t word,
                         const uint64_t mask)
{
  Layer& fl = layers_[z];
  const int tc = fl.track_chunk[i];
  const int tracks = fl.chunk_first_track[tc + 1] - fl.chunk_first_track[tc];
  for (int c = c0; c <= c1;) {
    const int k = fl.via_chunk[c];
    const int first = fl.via_chunk_first[k];
    const int last = fl.via_chunk_first[k + 1] - 1;
    const int end = std::min(c1, last);
    assign(fl.via[static_cast<size_t>(tc) * fl.pos_chunks + k],
           tracks,
           i - fl.chunk_first_track[tc],
           last - first + 1,
           c - first,
           end - first,
           word,
           mask);
    c = end + 1;
  }
}

bool FastGrid::edgeQuery(const int z,
                         const int i,
                         const int j,
                         const Kind edge_kind) const
{
  const int bit = edge_kind == kWireEdge ? kWireEdgeBit : kJogEdgeBit;
  return ((wiringWord(z, i, j) >> bit) & 1) != 0;
}

RipupLevel FastGrid::level(const int z,
                           const int i,
                           const int j,
                           const Kind kind,
                           const int wire_type,
                           LegalityCounters* counters) const
{
  const bool edge = kind == kWireEdge || kind == kJogEdge;
  if (wire_type >= num_wire_types_) {
    if (counters != nullptr) {
      counters->checker++;
    }
    if (oracle_ == nullptr) {
      return kNotRemovable;
    }
    return edge ? oracle_->edgeLevel(z, i, j, wire_type, kind)
                : oracle_->vertexLevel(z, i, j, wire_type, kind);
  }
  if (edge) {
    const uint64_t w = wiringWord(z, i, j);
    const bool along = kind == kWireEdge;
    if ((w >> (along ? kWireEdgeBit : kJogEdgeBit)) & 1) {
      if (counters != nullptr) {
        counters->checker++;
      }
      return oracle_ != nullptr ? oracle_->edgeLevel(z, i, j, wire_type, kind)
                                : kNotRemovable;
    }
    if (counters != nullptr) {
      counters->fast_grid++;
    }
    // Deduced from the vertex data at both endpoints.
    const int shift = wiringShift(wire_type, along ? 0 : 1);
    const uint64_t other
        = along ? wiringWord(z, i, j + 1) : wiringWord(z, i + 1, j);
    return std::max(field(w, shift), field(other, shift));
  }
  if (counters != nullptr) {
    counters->fast_grid++;
  }
  if (kind == kCutUp) {
    const int c = layers_[z].via_index[j];
    if (c < 0) {
      return kNotRemovable;
    }
    const uint64_t w = viaWord(z, i, c);
    return std::max(field(w, viaShift(wire_type, 0)),
                    field(w, viaShift(wire_type, 1)));
  }
  return field(wiringWord(z, i, j),
               wiringShift(wire_type, wiringKindIndex(kind)));
}

void FastGrid::setLevel(const size_t v,
                        const Kind kind,
                        const int wire_type,
                        const RipupLevel level)
{
  int z, i, j;
  graph_->decode(v, z, i, j);
  if (kind == kCutUp) {
    const int c = layers_[z].via_index[j];
    if (c >= 0) {
      const int shift = viaShift(wire_type, 0);
      assignVia(z, i, c, c, uint64_t{level} << shift, uint64_t{7} << shift);
    }
    return;
  }
  const int k = wiringKindIndex(kind);
  if (k < 0) {
    return;
  }
  const int shift = wiringShift(wire_type, k);
  assignWiring(z, i, j, j, uint64_t{level} << shift, uint64_t{7} << shift);
}

void FastGrid::setEdgeQuery(const size_t v,
                            const Kind edge_kind,
                            const bool query)
{
  int z, i, j;
  graph_->decode(v, z, i, j);
  const uint64_t bit = uint64_t{1}
                       << (edge_kind == kWireEdge ? kWireEdgeBit : kJogEdgeBit);
  assignWiring(z, i, j, j, query ? bit : 0, bit);
}

void FastGrid::bumpHistory(const size_t v)
{
  uint16_t& h = history_[v];
  if (h < 0xffff) {
    h++;
  }
}

size_t FastGrid::bytes(const std::vector<std::unique_ptr<TileRuns>>& tiles)
{
  size_t n = tiles.capacity() * sizeof(tiles[0]);
  for (const auto& t : tiles) {
    if (t) {
      n += sizeof(TileRuns) + t->track_first.capacity() * sizeof(uint32_t)
           + t->starts.capacity() * sizeof(uint16_t)
           + t->words.capacity() * sizeof(uint64_t);
    }
  }
  return n;
}

void FastGrid::compact()
{
  for (Layer& fl : layers_) {
    for (auto* tiles : {&fl.wiring, &fl.via}) {
      for (auto& t : *tiles) {
        if (t) {
          t->starts.shrink_to_fit();
          t->words.shrink_to_fit();
        }
      }
    }
  }
}

FastGrid::Statistics FastGrid::statistics() const
{
  Statistics s;
  for (size_t z = 0; z < layers_.size(); z++) {
    const Layer& fl = layers_[z];
    const TrackLayer& l = graph_->layers[z];
    s.vertices += static_cast<int64_t>(l.numTracks()) * l.numPositions();
    s.via_vertices
        += static_cast<int64_t>(l.numTracks()) * fl.via_position.size();
    s.bytes += bytes(fl.wiring) + bytes(fl.via)
               + (fl.pos_chunk.capacity() + fl.track_chunk.capacity()
                  + fl.via_index.capacity() + fl.via_position.capacity()
                  + fl.via_chunk.capacity())
                     * sizeof(int);
    // A track of a tile without runs is one all-zero interval.
    auto count = [&](const std::vector<std::unique_ptr<TileRuns>>& tiles,
                     int64_t& intervals,
                     bool marks) {
      for (int tc = 0; tc < fl.track_chunks; tc++) {
        const int tracks
            = fl.chunk_first_track[tc + 1] - fl.chunk_first_track[tc];
        for (int k = 0; k < fl.pos_chunks; k++) {
          const TileRuns* t
              = tiles.empty()
                    ? nullptr
                    : tiles[static_cast<size_t>(tc) * fl.pos_chunks + k].get();
          if (t == nullptr) {
            intervals += tracks;
            continue;
          }
          s.tiles++;
          for (int r = 0; r < tracks; r++) {
            intervals += std::max<uint32_t>(
                1, t->track_first[r + 1] - t->track_first[r]);
          }
          if (marks) {
            for (const uint64_t w : t->words) {
              s.marked_edges += (w & kEdgeBits) != 0;
            }
          }
        }
      }
    };
    count(fl.wiring, s.wiring_intervals, true);
    if (!fl.via.empty()) {
      count(fl.via, s.via_intervals, false);
    }
  }
  s.bytes += history_.size() * (sizeof(size_t) + sizeof(uint16_t) + 16);
  return s;
}

// --- FastGridBuilder ---------------------------------------------------

void FastGridBuilder::init(const TrackGraph* graph,
                           const std::vector<WireType>* types,
                           const int fast,
                           const DistanceRuleChecker* checker,
                           FastGrid* grid)
{
  graph_ = graph;
  types_ = types;
  fast_ = std::min(
      {fast, static_cast<int>(types->size()), FastGrid::kMaxWireTypes});
  checker_ = checker;
  grid_ = grid;
  grid_->setOracle(this);
  const int num_z = graph->numLayers();
  wiring_model_box_.assign(num_z, Rect{});
  cut_model_box_.assign(num_z, Rect{});
  projection_model_box_.assign(num_z, Rect{});
  has_projection_.assign(num_z, false);
  std::vector<LocatedModel> models;
  for (int z = 0; z < num_z; z++) {
    models.clear();
    pointModels(z, models);
    for (size_t k = 0; k < models.size(); k++) {
      wiring_model_box_[z]
          = k == 0 ? models[k].rect : wiring_model_box_[z].hull(models[k].rect);
    }
    models.clear();
    viaModels(z, models);
    bool first_cut = true;
    bool first_projection = true;
    for (const LocatedModel& m : models) {
      if (m.layer == viaGridLayer(z)) {
        cut_model_box_[z] = first_cut ? m.rect : cut_model_box_[z].hull(m.rect);
        first_cut = false;
      } else {
        projection_model_box_[z]
            = first_projection ? m.rect : projection_model_box_[z].hull(m.rect);
        first_projection = false;
        has_projection_[z] = true;
      }
    }
  }
}

Point FastGridBuilder::location(const int z, const int i, const int j) const
{
  const auto [x, y] = graph_->layers[z].point(i, j);
  return {x, y};
}

void FastGridBuilder::pointModels(const int z,
                                  std::vector<LocatedModel>& out) const
{
  const int num_z = graph_->numLayers();
  const int layer = wiringGridLayer(z);
  for (int t = 0; t < fast_; t++) {
    const WireType& w = (*types_)[t];
    auto add = [&](const WireModel& m, int kind) {
      if (m.valid) {
        out.push_back({layer, m.rect, m.cls, 3, 4 * t + kind});
      }
    };
    add(w.preferred[z], 0);
    add(w.non_preferred[z], 1);
    if (z + 1 < num_z && w.via[z].valid) {
      add(w.via[z].bottom, 2);
    }
    if (z > 0 && w.via[z - 1].valid) {
      add(w.via[z - 1].top, 3);
    }
  }
}

void FastGridBuilder::viaModels(const int z,
                                std::vector<LocatedModel>& out) const
{
  if (z + 1 >= graph_->numLayers()) {
    return;
  }
  for (int t = 0; t < fast_; t++) {
    const ViaModel& v = (*types_)[t].via[z];
    if (!v.valid) {
      continue;
    }
    out.push_back({viaGridLayer(z), v.cut.rect, v.cut.cls, 3, 2 * t});
    if (v.projection.valid) {
      out.push_back({viaGridLayer(z + 1),
                     v.projection.rect,
                     v.projection.cls,
                     3,
                     2 * t + 1});
    }
  }
}

namespace {

// Writes runs of equal flags over [lo, lo + flags.size()) of track i.
template <typename Assign>
void assignFlags(const std::vector<uint8_t>& flags,
                 const int lo,
                 const uint64_t bit,
                 Assign&& assign)
{
  for (size_t k = 0; k < flags.size();) {
    size_t e = k;
    while (e + 1 < flags.size() && flags[e + 1] == flags[k]) {
      e++;
    }
    assign(lo + static_cast<int>(k),
           lo + static_cast<int>(e),
           flags[k] ? bit : uint64_t{0},
           bit);
    k = e + 1;
  }
}

}  // namespace

void FastGridBuilder::recomputeWiring(const int z,
                                      const Range& r,
                                      const bool vertices,
                                      const bool jogs) const
{
  const TrackLayer& layer = graph_->layers[z];
  const int np = layer.numPositions();
  thread_local std::vector<LocatedModel> models;
  thread_local std::vector<int> run_of;
  thread_local std::vector<std::pair<Coord, Coord>> run_span;
  thread_local std::vector<uint8_t> flags;
  CheckResult res;
  CheckOptions opts;
  opts.levels = true;
  // Only the intervals along the tracks are needed.
  opts.x_interval = layer.horizontal;
  opts.y_interval = !layer.horizontal;
  const uint64_t level_mask
      = fast_ == 0 ? 0 : (uint64_t{1} << (12 * fast_)) - 1;
  if (vertices) {
    models.clear();
    pointModels(z, models);
    opts.groups = 4 * fast_;
    for (int i = r.i0; i <= r.i1; i++) {
      // Vertex data from maximal intervals of equal answers.
      run_of.assign(std::max(0, r.j1 - r.j0 + 1), -1);
      run_span.clear();
      int run = 0;
      for (int j = r.j0; j <= r.j1; run++) {
        checker_->check(location(z, i, j), models, opts, res);
        uint64_t word = 0;
        for (int g = 0; g < opts.groups; g++) {
          word |= uint64_t{res.levels[g]}
                  << FastGrid::wiringShift(g / 4, g % 4);
        }
        const Coord hi = layer.horizontal ? res.x_hi : res.y_hi;
        run_span.emplace_back(layer.horizontal ? res.x_lo : res.y_lo, hi);
        int end = static_cast<int>(std::upper_bound(layer.positions.begin(),
                                                    layer.positions.end(),
                                                    hi)
                                   - layer.positions.begin())
                  - 1;
        end = std::clamp(end, j, r.j1);
        grid_->assignWiring(z, i, j, end, word, level_mask);
        for (int k = j; k <= end; k++) {
          run_of[k - r.j0] = run;
        }
        j = end + 1;
      }
      // Preferred-direction edges: inside one interval, or across two
      // intervals that together cover the edge, the point levels on the
      // edge are those of its endpoints, and every violation of the swept
      // wire is one of a point model on it, so the deduction holds;
      // otherwise the swept wire is checked.
      const int e1 = std::min(r.e1, np - 2);
      if (r.e0 > e1) {
        continue;
      }
      flags.assign(e1 - r.e0 + 1, 0);
      for (int j = r.e0; j <= e1; j++) {
        if (j >= r.j0 && j + 1 <= r.j1) {
          const int a = run_of[j - r.j0];
          const int b = run_of[j + 1 - r.j0];
          if (a == b
              || (b == a + 1
                  && int64_t{run_span[a].second} + 1 >= run_span[b].first)) {
            continue;
          }
        }
        const uint64_t wa = grid_->wiringWord(z, i, j);
        const uint64_t wb = grid_->wiringWord(z, i, j + 1);
        for (int t = 0; t < fast_; t++) {
          const int shift = FastGrid::wiringShift(t, 0);
          const RipupLevel ends = std::max(FastGrid::field(wa, shift),
                                           FastGrid::field(wb, shift));
          if (sweptLevel(z, i, j, t, kWireEdge, kNoNet) > ends) {
            flags[j - r.e0] = 1;
            break;
          }
        }
      }
      assignFlags(flags,
                  r.e0,
                  uint64_t{1} << FastGrid::kWireEdgeBit,
                  [&](int lo, int hi, uint64_t w, uint64_t m) {
                    grid_->assignWiring(z, i, lo, hi, w, m);
                  });
    }
  }
  if (jogs && layer.jogs && r.j0 <= r.j1) {
    for (int p = std::max(r.p0, 0); p <= std::min(r.p1, layer.numTracks() - 2);
         p++) {
      jogFlags(z, p, r.j0, r.j1, flags);
      writeFlags(z, p, r.j0, flags, FastGrid::kJogEdgeBit);
    }
  }
}

void FastGridBuilder::jogFlags(const int z,
                               const int p,
                               const int j0,
                               const int j1,
                               std::vector<uint8_t>& flags) const
{
  // Jog edges (p, p + 1): the swept non-preferred wire as one model.
  const TrackLayer& layer = graph_->layers[z];
  thread_local std::vector<LocatedModel> models;
  models.clear();
  const Coord d = layer.tracks[p + 1] - layer.tracks[p];
  for (int t = 0; t < fast_; t++) {
    const WireModel& m = (*types_)[t].non_preferred[z];
    if (!m.valid) {
      continue;
    }
    const Rect shifted
        = m.rect.translated(layer.horizontal ? Point{0, d} : Point{d, 0});
    models.push_back({wiringGridLayer(z), m.rect.hull(shifted), m.cls, 3, t});
  }
  CheckOptions opts;
  opts.levels = true;
  opts.groups = fast_;
  opts.x_interval = layer.horizontal;
  opts.y_interval = !layer.horizontal;
  CheckResult res;
  flags.assign(j1 - j0 + 1, 0);
  for (int j = j0; j <= j1;) {
    checker_->check(location(z, p, j), models, opts, res);
    const Coord hi = layer.horizontal ? res.x_hi : res.y_hi;
    int end = static_cast<int>(std::upper_bound(layer.positions.begin(),
                                                layer.positions.end(),
                                                hi)
                               - layer.positions.begin())
              - 1;
    end = std::clamp(end, j, j1);
    for (int k = j; k <= end; k++) {
      const uint64_t wa = grid_->wiringWord(z, p, k);
      const uint64_t wb = grid_->wiringWord(z, p + 1, k);
      for (int t = 0; t < fast_; t++) {
        const int shift = FastGrid::wiringShift(t, 1);
        if (res.levels[t] > std::max(FastGrid::field(wa, shift),
                                     FastGrid::field(wb, shift))) {
          flags[k - j0] = 1;
          break;
        }
      }
    }
    j = end + 1;
  }
}

void FastGridBuilder::writeFlags(const int z,
                                 const int p,
                                 const int j0,
                                 const std::vector<uint8_t>& flags,
                                 const int bit) const
{
  assignFlags(flags,
              j0,
              uint64_t{1} << bit,
              [&](int lo, int hi, uint64_t w, uint64_t m) {
                grid_->assignWiring(z, p, lo, hi, w, m);
              });
}

void FastGridBuilder::recomputeVia(const int z,
                                   const int i0,
                                   const int i1,
                                   const int c0,
                                   const int c1) const
{
  if (z + 1 >= graph_->numLayers() || c0 > c1) {
    return;
  }
  thread_local std::vector<LocatedModel> models;
  models.clear();
  viaModels(z, models);
  const TrackLayer& layer = graph_->layers[z];
  CheckOptions opts;
  opts.levels = true;
  opts.groups = 2 * fast_;
  opts.x_interval = layer.horizontal;
  opts.y_interval = !layer.horizontal;
  const uint64_t mask = fast_ == 0 ? 0 : (uint64_t{1} << (6 * fast_)) - 1;
  CheckResult res;
  for (int i = i0; i <= i1; i++) {
    for (int c = c0; c <= c1;) {
      const int j = grid_->viaPosition(z, c);
      checker_->check(location(z, i, j), models, opts, res);
      uint64_t word = 0;
      for (int g = 0; g < opts.groups; g++) {
        word |= uint64_t{res.levels[g]} << FastGrid::viaShift(g / 2, g % 2);
      }
      const Coord hi = layer.horizontal ? res.x_hi : res.y_hi;
      int end = c;
      while (end + 1 <= c1
             && layer.positions[grid_->viaPosition(z, end + 1)] <= hi) {
        end++;
      }
      grid_->assignVia(z, i, c, end, word, mask);
      c = end + 1;
    }
  }
}

void FastGridBuilder::rebuildVertices(const int z, const int i0, const int i1)
{
  const TrackLayer& layer = graph_->layers[z];
  Range r{i0,
          i1,
          0,
          layer.numPositions() - 1,
          i0,
          std::min(i1, layer.numTracks() - 2),
          0,
          layer.numPositions() - 2};
  recomputeWiring(z, r, true, false);
  if (z + 1 < graph_->numLayers()) {
    recomputeVia(z, i0, i1, 0, grid_->numViaPositions(z) - 1);
  }
}

std::vector<std::array<int, 3>> FastGridBuilder::rebuildJobs() const
{
  std::vector<std::array<int, 3>> jobs;
  for (int z = 0; z < graph_->numLayers(); z++) {
    const TrackLayer& l = graph_->layers[z];
    for (int i = 0; i < l.numTracks();) {
      int e = i;
      while (e + 1 < l.numTracks() && l.track_tile.size() == l.tracks.size()
             && l.track_tile[e + 1] == l.track_tile[i]) {
        e++;
      }
      jobs.push_back({z, i, e});
      i = e + 1;
    }
  }
  return jobs;
}

void FastGridBuilder::computeJogBits(const int z,
                                     const int i0,
                                     const int i1,
                                     std::vector<JogBits>& out) const
{
  const TrackLayer& layer = graph_->layers[z];
  if (!layer.jogs || layer.numPositions() == 0) {
    return;
  }
  for (int p = i0; p <= std::min(i1, layer.numTracks() - 2); p++) {
    JogBits bits{z, p, {}};
    jogFlags(z, p, 0, layer.numPositions() - 1, bits.flags);
    out.push_back(std::move(bits));
  }
}

void FastGridBuilder::writeJogBits(const std::vector<JogBits>& bits) const
{
  for (const JogBits& b : bits) {
    writeFlags(b.z, b.p, 0, b.flags, FastGrid::kJogEdgeBit);
  }
}

void FastGridBuilder::rebuildJogBits(const int z, const int i0, const int i1)
{
  std::vector<JogBits> bits;
  computeJogBits(z, i0, i1, bits);
  writeJogBits(bits);
}

namespace {

// Index range of sorted coordinates inside [lo, hi].
std::pair<int, int> indexRange(const std::vector<Coord>& v,
                               const int64_t lo,
                               const int64_t hi)
{
  const int a = static_cast<int>(
      std::lower_bound(
          v.begin(), v.end(), lo, [](Coord c, int64_t x) { return c < x; })
      - v.begin());
  const int b = static_cast<int>(
                    std::upper_bound(v.begin(),
                                     v.end(),
                                     hi,
                                     [](int64_t x, Coord c) { return x < c; })
                    - v.begin())
                - 1;
  return {a, b};
}

}  // namespace

FastGridBuilder::Range FastGridBuilder::wiringRange(const int z,
                                                    const Rect& rect,
                                                    const int grid_layer) const
{
  const TrackLayer& layer = graph_->layers[z];
  const Rect& box = wiring_model_box_[z];
  const int64_t reach = checker_->grid().classes().reach(grid_layer) + 1;
  const int64_t x_lo = int64_t{rect.xlo} - box.xhi - reach;
  const int64_t x_hi = int64_t{rect.xhi} - box.xlo + reach;
  const int64_t y_lo = int64_t{rect.ylo} - box.yhi - reach;
  const int64_t y_hi = int64_t{rect.yhi} - box.ylo + reach;
  const auto [i0, i1] = layer.horizontal ? indexRange(layer.tracks, y_lo, y_hi)
                                         : indexRange(layer.tracks, x_lo, x_hi);
  const auto [j0, j1] = layer.horizontal
                            ? indexRange(layer.positions, x_lo, x_hi)
                            : indexRange(layer.positions, y_lo, y_hi);
  Range r;
  r.i0 = i0;
  r.i1 = i1;
  r.j0 = j0;
  r.j1 = j1;
  // Jog pairs and preferred edges straddling the box.
  r.p0 = std::max(0, i0 - 1);
  r.p1 = std::min(i1, layer.numTracks() - 2);
  r.e0 = std::max(0, j0 - 1);
  r.e1 = std::min(j1, layer.numPositions() - 2);
  return r;
}

void FastGridBuilder::shapeChanged(const int grid_layer, const Rect& rect)
{
  if (grid_layer % 2 == 0) {
    const int z = grid_layer / 2;
    if (z < graph_->numLayers()) {
      const Range r = wiringRange(z, rect, grid_layer);
      recomputeWiring(z, r, true, true);
    }
    return;
  }
  const int64_t reach = checker_->grid().classes().reach(grid_layer) + 1;
  auto via_update = [&](int z, const Rect& box) {
    const TrackLayer& layer = graph_->layers[z];
    const int64_t x_lo = int64_t{rect.xlo} - box.xhi - reach;
    const int64_t x_hi = int64_t{rect.xhi} - box.xlo + reach;
    const int64_t y_lo = int64_t{rect.ylo} - box.yhi - reach;
    const int64_t y_hi = int64_t{rect.yhi} - box.ylo + reach;
    const auto [i0, i1] = layer.horizontal
                              ? indexRange(layer.tracks, y_lo, y_hi)
                              : indexRange(layer.tracks, x_lo, x_hi);
    const auto [j0, j1] = layer.horizontal
                              ? indexRange(layer.positions, x_lo, x_hi)
                              : indexRange(layer.positions, y_lo, y_hi);
    if (j0 > j1) {
      return;
    }
    int c0 = 0;
    while (c0 < grid_->numViaPositions(z) && grid_->viaPosition(z, c0) < j0) {
      c0++;
    }
    int c1 = c0 - 1;
    while (c1 + 1 < grid_->numViaPositions(z)
           && grid_->viaPosition(z, c1 + 1) <= j1) {
      c1++;
    }
    recomputeVia(z, i0, i1, c0, c1);
  };
  const int z = (grid_layer - 1) / 2;
  if (z + 1 < graph_->numLayers()) {
    via_update(z, cut_model_box_[z]);
  }
  if (z >= 1 && has_projection_[z - 1]) {
    via_update(z - 1, projection_model_box_[z - 1]);
  }
}

RipupLevel FastGridBuilder::vertexLevel(const int z,
                                        const int i,
                                        const int j,
                                        const int wire_type,
                                        const Kind kind) const
{
  return exactLevel(z, i, j, wire_type, kind, kNoNet);
}

RipupLevel FastGridBuilder::edgeLevel(const int z,
                                      const int i,
                                      const int j,
                                      const int wire_type,
                                      const Kind kind) const
{
  return exactLevel(z, i, j, wire_type, kind, kNoNet);
}

RipupLevel FastGridBuilder::exactLevel(const int z,
                                       const int i,
                                       const int j,
                                       const int wire_type,
                                       const Kind kind,
                                       const int32_t ignored_owner) const
{
  const WireType& w = (*types_)[wire_type];
  const int num_z = graph_->numLayers();
  thread_local std::vector<LocatedModel> models;
  models.clear();
  auto add = [&](int layer, const WireModel& m) {
    if (m.valid) {
      models.push_back({layer, m.rect, m.cls, 3, 0});
    }
  };
  switch (kind) {
    case kWire:
      add(wiringGridLayer(z), w.preferred[z]);
      break;
    case kJog:
      add(wiringGridLayer(z), w.non_preferred[z]);
      break;
    case kViaUpPad:
      if (z + 1 < num_z && w.via[z].valid) {
        add(wiringGridLayer(z), w.via[z].bottom);
      }
      break;
    case kViaDownPad:
      if (z > 0 && w.via[z - 1].valid) {
        add(wiringGridLayer(z), w.via[z - 1].top);
      }
      break;
    case kCutUp:
      if (z + 1 < num_z && w.via[z].valid) {
        add(viaGridLayer(z), w.via[z].cut);
        add(viaGridLayer(z + 1), w.via[z].projection);
      }
      break;
    default: {
      // Edge kinds: both endpoints and the swept wire.
      const bool along = kind == kWireEdge;
      const Kind vertex_kind = along ? kWire : kJog;
      return std::max(
          {exactLevel(z, i, j, wire_type, vertex_kind, ignored_owner),
           along
               ? exactLevel(z, i, j + 1, wire_type, vertex_kind, ignored_owner)
               : exactLevel(z, i + 1, j, wire_type, vertex_kind, ignored_owner),
           sweptLevel(z, i, j, wire_type, kind, ignored_owner)});
    }
  }
  if (models.empty()) {
    return kNotRemovable;
  }
  CheckOptions opts;
  opts.intervals = false;
  opts.ignored_owner = ignored_owner;
  CheckResult res;
  checker_->check(location(z, i, j), models, opts, res);
  return res.levels[0];
}

RipupLevel FastGridBuilder::sweptLevel(const int z,
                                       const int i,
                                       const int j,
                                       const int wire_type,
                                       const Kind kind,
                                       const int32_t ignored_owner) const
{
  const TrackLayer& layer = graph_->layers[z];
  const bool along = kind == kWireEdge;
  const WireModel& m = along ? (*types_)[wire_type].preferred[z]
                             : (*types_)[wire_type].non_preferred[z];
  if (!m.valid) {
    return kNotRemovable;
  }
  const Coord d = along ? layer.positions[j + 1] - layer.positions[j]
                        : layer.tracks[i + 1] - layer.tracks[i];
  // Along the preferred direction for kWireEdge, across it for kJogEdge.
  const bool shift_x = along == layer.horizontal;
  const Rect swept
      = m.rect.hull(m.rect.translated(shift_x ? Point{d, 0} : Point{0, d}));
  const LocatedModel model{wiringGridLayer(z), swept, m.cls, 3, 0};
  CheckOptions opts;
  opts.intervals = false;
  opts.ignored_owner = ignored_owner;
  CheckResult res;
  checker_->check(location(z, i, j), std::span(&model, 1), opts, res);
  return res.levels[0];
}

}  // namespace drt::br
