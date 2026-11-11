// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrFastGrid.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <tuple>
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

void FastGrid::assignWords(std::unique_ptr<TileRuns>& tile,
                           const int tracks,
                           const int track,
                           const int count,
                           const int lo,
                           const uint64_t* values,
                           const int n,
                           const uint64_t mask)
{
  if (!tile) {
    bool zero = true;
    for (int k = 0; k < n && zero; k++) {
      zero = (values[k] & mask) == 0;
    }
    if (zero) {
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
  bool changed = false;
  for (int k = 0; k < n; k++) {
    const uint64_t w = (words[lo + k] & ~mask) | (values[k] & mask);
    changed |= w != words[lo + k];
    words[lo + k] = w;
  }
  if (!changed) {
    return;
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

void FastGrid::writeWiring(const int z,
                           const int i,
                           const int j0,
                           const uint64_t* values,
                           const int n,
                           const uint64_t mask)
{
  Layer& fl = layers_[z];
  const int tc = fl.track_chunk[i];
  const int tracks = fl.chunk_first_track[tc + 1] - fl.chunk_first_track[tc];
  for (int j = j0; j < j0 + n;) {
    const int k = fl.pos_chunk[j];
    const int first = fl.chunk_first_pos[k];
    const int last = fl.chunk_first_pos[k + 1] - 1;
    const int end = std::min(j0 + n - 1, last);
    assignWords(fl.wiring[static_cast<size_t>(tc) * fl.pos_chunks + k],
                tracks,
                i - fl.chunk_first_track[tc],
                last - first + 1,
                j - first,
                values + (j - j0),
                end - j + 1,
                mask);
    j = end + 1;
  }
}

void FastGrid::writeVia(const int z,
                        const int i,
                        const int c0,
                        const uint64_t* values,
                        const int n,
                        const uint64_t mask)
{
  Layer& fl = layers_[z];
  const int tc = fl.track_chunk[i];
  const int tracks = fl.chunk_first_track[tc + 1] - fl.chunk_first_track[tc];
  for (int c = c0; c < c0 + n;) {
    const int k = fl.via_chunk[c];
    const int first = fl.via_chunk_first[k];
    const int last = fl.via_chunk_first[k + 1] - 1;
    const int end = std::min(c0 + n - 1, last);
    assignWords(fl.via[static_cast<size_t>(tc) * fl.pos_chunks + k],
                tracks,
                i - fl.chunk_first_track[tc],
                last - first + 1,
                c - first,
                values + (c - c0),
                end - c + 1,
                mask);
    c = end + 1;
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
  return detail::levelOf(
      *this,
      z,
      i,
      j,
      kind,
      wire_type,
      kNoNet,
      counters,
      [this](int zz, int ii, int jj) { return wiringWord(zz, ii, jj); },
      [this](int zz, int ii, int cc) { return viaWord(zz, ii, cc); });
}

uint64_t FastGrid::lookupRun(const TileRuns* tile,
                             const int track,
                             const int index,
                             const int count,
                             int& lo,
                             int& hi)
{
  if (tile == nullptr) {
    lo = 0;
    hi = count - 1;
    return 0;
  }
  const uint32_t a = tile->track_first[track];
  const uint32_t b = tile->track_first[track + 1];
  if (a == b) {
    lo = 0;
    hi = count - 1;
    return 0;
  }
  const auto begin = tile->starts.begin() + a;
  const auto it = std::upper_bound(begin, tile->starts.begin() + b, index);
  const size_t r = (it - tile->starts.begin()) - 1;
  lo = tile->starts[r];
  hi = (r + 1 < b ? tile->starts[r + 1] : count) - 1;
  return tile->words[r];
}

uint64_t FastGrid::wiringRun(const int z,
                             const int i,
                             const int j,
                             int& lo,
                             int& hi) const
{
  const Layer& fl = layers_[z];
  const int k = fl.pos_chunk[j];
  const int tc = fl.track_chunk[i];
  const int first = fl.chunk_first_pos[k];
  const uint64_t w
      = lookupRun(fl.wiring[static_cast<size_t>(tc) * fl.pos_chunks + k].get(),
                  i - fl.chunk_first_track[tc],
                  j - first,
                  fl.chunk_first_pos[k + 1] - first,
                  lo,
                  hi);
  lo += first;
  hi += first;
  return w;
}

uint64_t FastGrid::viaRun(const int z,
                          const int i,
                          const int c,
                          int& lo,
                          int& hi) const
{
  const Layer& fl = layers_[z];
  const int k = fl.via_chunk[c];
  const int tc = fl.track_chunk[i];
  const int first = fl.via_chunk_first[k];
  const uint64_t w
      = lookupRun(fl.via[static_cast<size_t>(tc) * fl.pos_chunks + k].get(),
                  i - fl.chunk_first_track[tc],
                  c - first,
                  fl.via_chunk_first[k + 1] - first,
                  lo,
                  hi);
  lo += first;
  hi += first;
  return w;
}

// --- FastGridReader ----------------------------------------------------

void FastGridReader::reset(const FastGrid* grid, const FastGridOverlay* overlay)
{
  grid_ = grid;
  overlay_ = overlay != nullptr && overlay->owner != kNoNet ? overlay : nullptr;
  for (Cursor& c : wiring_) {
    c.z = -1;
  }
  for (Cursor& c : via_) {
    c.z = -1;
  }
}

RipupLevel FastGridReader::level(const int z,
                                 const int i,
                                 const int j,
                                 const Kind kind,
                                 const int wire_type,
                                 LegalityCounters* counters)
{
  return detail::levelOf(
      *grid_,
      z,
      i,
      j,
      kind,
      wire_type,
      ignoredOwner(),
      counters,
      [this](int zz, int ii, int jj) { return wiringWord(zz, ii, jj); },
      [this](int zz, int ii, int cc) { return viaWord(zz, ii, cc); });
}

uint64_t FastGridReader::applyOverlay(
    const std::unordered_map<uint64_t, std::vector<FastGridOverlay::Run>>& runs,
    const int z,
    const int i,
    const int j,
    const uint64_t word,
    int& lo,
    int& hi)
{
  const auto it = runs.find(FastGridOverlay::key(z, i));
  if (it == runs.end()) {
    return word;
  }
  const std::vector<FastGridOverlay::Run>& r = it->second;
  // First overlay run ending at or after j.
  const auto at = std::lower_bound(
      r.begin(), r.end(), j, [](const FastGridOverlay::Run& run, int x) {
        return run.hi < x;
      });
  if (at != r.end() && at->lo <= j) {
    lo = std::max(lo, at->lo);
    hi = std::min(hi, at->hi);
    return at->word;
  }
  if (at != r.end()) {
    hi = std::min(hi, at->lo - 1);
  }
  if (at != r.begin()) {
    lo = std::max(lo, (at - 1)->hi + 1);
  }
  return word;
}

FastGrid::RunsView FastGrid::wiringRuns(const int z,
                                        const int i,
                                        const int j) const
{
  const Layer& fl = layers_[z];
  const int k = fl.pos_chunk[j];
  const int tc = fl.track_chunk[i];
  RunsView v;
  v.first = fl.chunk_first_pos[k];
  v.count = fl.chunk_first_pos[k + 1] - v.first;
  const TileRuns* t
      = fl.wiring[static_cast<size_t>(tc) * fl.pos_chunks + k].get();
  if (t != nullptr) {
    const int track = i - fl.chunk_first_track[tc];
    const uint32_t a = t->track_first[track];
    v.n = static_cast<int>(t->track_first[track + 1] - a);
    v.starts = t->starts.data() + a;
    v.words = t->words.data() + a;
  }
  return v;
}

uint64_t FastGridReader::fillWiring(Cursor& c,
                                    const int z,
                                    const int i,
                                    const int j)
{
  // Same track and tile: the next run, or a search in the cached runs.
  if (c.z == z && c.i == i && c.run >= 0 && j >= c.runs.first
      && j < c.runs.first + c.runs.count) {
    const FastGrid::RunsView& v = c.runs;
    const int x = j - v.first;
    int r = c.run;
    if (r + 1 < v.n && x >= v.starts[r + 1]
        && (r + 2 >= v.n || x < v.starts[r + 2])) {
      r++;
    } else {
      r = static_cast<int>(std::upper_bound(v.starts, v.starts + v.n, x)
                           - v.starts)
          - 1;
    }
    c.run = r;
    c.lo = v.first + v.starts[r];
    c.hi = v.first + (r + 1 < v.n ? v.starts[r + 1] : v.count) - 1;
    c.word = v.words[r];
    return c.word;
  }
  const bool overlaid
      = overlay_ != nullptr
        && overlay_->wiring.contains(FastGridOverlay::key(z, i));
  if (!overlaid) {
    c.z = z;
    c.i = i;
    c.runs = grid_->wiringRuns(z, i, j);
    const FastGrid::RunsView& v = c.runs;
    if (v.n == 0) {
      c.run = -1;
      c.lo = v.first;
      c.hi = v.first + v.count - 1;
      c.word = 0;
      return 0;
    }
    const int x = j - v.first;
    const int r = static_cast<int>(std::upper_bound(v.starts, v.starts + v.n, x)
                                   - v.starts)
                  - 1;
    c.run = r;
    c.lo = v.first + v.starts[r];
    c.hi = v.first + (r + 1 < v.n ? v.starts[r + 1] : v.count) - 1;
    c.word = v.words[r];
    return c.word;
  }
  int lo;
  int hi;
  uint64_t w = grid_->wiringRun(z, i, j, lo, hi);
  w = applyOverlay(overlay_->wiring, z, i, j, w, lo, hi);
  c.z = z;
  c.i = i;
  c.lo = lo;
  c.hi = hi;
  c.word = w;
  c.run = -1;
  return w;
}

uint64_t FastGridReader::fillVia(Cursor& c,
                                 const int z,
                                 const int i,
                                 const int v)
{
  int lo;
  int hi;
  uint64_t w = grid_->viaRun(z, i, v, lo, hi);
  if (overlay_ != nullptr) {
    w = applyOverlay(overlay_->via, z, i, v, w, lo, hi);
  }
  c.z = z;
  c.i = i;
  c.lo = lo;
  c.hi = hi;
  c.word = w;
  c.run = -1;
  return w;
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

inline void raise(uint64_t& word, const int shift, const RipupLevel level)
{
  if (FastGrid::field(word, shift) < level) {
    word = (word & ~(uint64_t{7} << shift)) | (uint64_t{level} << shift);
  }
}

Rect windowOf(const Point a, const Point b, const Rect& box, const Coord reach)
{
  return {std::min(a.x, b.x) + box.xlo - reach,
          std::min(a.y, b.y) + box.ylo - reach,
          std::max(a.x, b.x) + box.xhi + reach,
          std::max(a.y, b.y) + box.yhi + reach};
}

// Largest distance at which shapes of the two classes interact (spacing
// and line-end extensions), cached per thread.
Coord interaction(const ShapeClassTable& classes,
                  const ClassId model,
                  const ClassId piece)
{
  struct Entry
  {
    uint64_t generation = 0;
    uint32_t key = 0xffffffff;
    Coord d = 0;
  };
  thread_local std::array<Entry, 4096> cache;
  const uint32_t key = (static_cast<uint32_t>(model) << 16) | piece;
  Entry& e = cache[(model * 131u + piece) & 4095u];
  if (e.generation != classes.generation() || e.key != key) {
    const Requirement r = classes.requirement(model, piece);
    Coord d = std::max({r.facing_x, r.facing_y, r.diagonal});
    if (r.line_end) {
      const LineEndRule& le = classes.layer(classes.get(model).layer).line_end;
      d = std::max({d, le.space, le.within});
    }
    if (r.overlap) {
      d = std::max<Coord>(d, 1);
    }
    e = {classes.generation(), key, d};
  }
  return e.d;
}

// Exact prefilter: the model on the line (cross range [c_lo, c_hi] and along
// range [a_lo, a_hi] of all its placements) cannot violate the piece.
inline bool farApart(const ShapeClassTable& classes,
                     const LocatedModel& m,
                     const GridPiece& p,
                     const bool horizontal,
                     const Coord c_lo,
                     const Coord c_hi,
                     const Coord a_lo,
                     const Coord a_hi)
{
  const Coord d = interaction(classes, m.cls, p.cls);
  const Coord pc_lo = horizontal ? p.rect.ylo : p.rect.xlo;
  const Coord pc_hi = horizontal ? p.rect.yhi : p.rect.xhi;
  const Coord pa_lo = horizontal ? p.rect.xlo : p.rect.ylo;
  const Coord pa_hi = horizontal ? p.rect.xhi : p.rect.yhi;
  return pc_lo >= c_hi + d || pc_hi <= c_lo - d || pa_lo >= a_hi + d
         || pa_hi <= a_lo - d;
}

}  // namespace

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

LocatedModel FastGridBuilder::jogHull(const int z,
                                      const int p,
                                      const int t) const
{
  // The swept non-preferred wire from track p to p + 1 as one model.
  const TrackLayer& layer = graph_->layers[z];
  const WireModel& m = (*types_)[t].non_preferred[z];
  const Coord d = layer.tracks[p + 1] - layer.tracks[p];
  const Rect shifted
      = m.rect.translated(layer.horizontal ? Point{0, d} : Point{d, 0});
  return {wiringGridLayer(z), m.rect.hull(shifted), m.cls, 3, t};
}

Rect FastGridBuilder::prefHull(const int z,
                               const int i,
                               const int j,
                               const int t) const
{
  const TrackLayer& layer = graph_->layers[z];
  const Rect& m = (*types_)[t].preferred[z].rect;
  const Coord d = layer.positions[j + 1] - layer.positions[j];
  return m.hull(m.translated(layer.horizontal ? Point{d, 0} : Point{0, d}))
      .translated(location(z, i, j));
}

void FastGridBuilder::trackLevels(const int z,
                                  const int i,
                                  const int j0,
                                  const int j1,
                                  const int32_t ignored,
                                  TrackLevels& out) const
{
  const TrackLayer& layer = graph_->layers[z];
  const std::vector<Coord>& pos = layer.positions;
  const int n = j1 - j0 + 1;
  const int edges = std::max(0, n - 1);
  out.words.assign(n, 0);
  out.interior.assign(static_cast<size_t>(fast_) * edges, 0);
  out.pieces.clear();
  thread_local std::vector<LocatedModel> models;
  models.clear();
  pointModels(z, models);
  if (models.empty() || n <= 0) {
    return;
  }
  const ShapeGrid& sg = checker_->grid();
  const ShapeClassTable& classes = sg.classes();
  const int g = wiringGridLayer(z);
  const Point a = location(z, i, j0);
  const Point b = location(z, i, j1);
  sg.collect(
      g, windowOf(a, b, wiring_model_box_[z], classes.reach(g)), out.pieces);
  thread_local std::vector<LineSpan> spans;
  for (const GridPiece& p : out.pieces.pieces) {
    if (p.owner == ignored && p.owner != kNoNet && isOwned(p.type)) {
      continue;
    }
    const RipupLevel lv = sg.level(p.type, p.owner);
    const PlacedShape ps{p.rect, p.cls, p.line_end_ends};
    for (const LocatedModel& m : models) {
      spans.clear();
      appendViolationSpans(
          classes, m.rect, m.cls, 3, a, layer.horizontal, ps, spans);
      const int shift = FastGrid::wiringShift(m.group / 4, m.group % 4);
      for (const LineSpan& sp : spans) {
        const auto [ja, jb] = indexRange(pos, sp.lo, sp.hi);
        for (int k = std::max(ja, j0); k <= std::min(jb, j1); k++) {
          raise(out.words[k - j0], shift, lv);
        }
        if (m.group % 4 == 0 && edges > 0) {
          // Edges (e, e + 1) whose segment meets the span.
          const int e_lo = std::max(j0, ja - 1);
          const int e_hi = std::min(j1 - 1, jb);
          RipupLevel* in
              = &out.interior[static_cast<size_t>(m.group / 4) * edges];
          for (int e = e_lo; e <= e_hi; e++) {
            in[e - j0] = std::max(in[e - j0], lv);
          }
        }
      }
    }
  }
}

void FastGridBuilder::jogLevels(const int z,
                                const int p,
                                const int j0,
                                const int j1,
                                const int32_t ignored,
                                std::vector<RipupLevel>& out) const
{
  const TrackLayer& layer = graph_->layers[z];
  const int n = j1 - j0 + 1;
  out.assign(static_cast<size_t>(fast_) * std::max(0, n), 0);
  if (n <= 0) {
    return;
  }
  thread_local std::vector<LocatedModel> models;
  models.clear();
  Rect box;
  for (int t = 0; t < fast_; t++) {
    if ((*types_)[t].non_preferred[z].valid) {
      models.push_back(jogHull(z, p, t));
      box = models.size() == 1 ? models.back().rect
                               : box.hull(models.back().rect);
    }
  }
  if (models.empty()) {
    return;
  }
  const ShapeGrid& sg = checker_->grid();
  const ShapeClassTable& classes = sg.classes();
  const int g = wiringGridLayer(z);
  const Point a = location(z, p, j0);
  const Point b = location(z, p, j1);
  thread_local GridWindow window;
  sg.collect(g, windowOf(a, b, box, classes.reach(g)), window);
  thread_local std::vector<LineSpan> spans;
  for (const GridPiece& piece : window.pieces) {
    if (piece.owner == ignored && piece.owner != kNoNet
        && isOwned(piece.type)) {
      continue;
    }
    const RipupLevel lv = sg.level(piece.type, piece.owner);
    const PlacedShape ps{piece.rect, piece.cls, piece.line_end_ends};
    for (const LocatedModel& m : models) {
      spans.clear();
      appendViolationSpans(
          classes, m.rect, m.cls, 3, a, layer.horizontal, ps, spans);
      for (const LineSpan& sp : spans) {
        const auto [ja, jb] = indexRange(layer.positions, sp.lo, sp.hi);
        for (int k = std::max(ja, j0); k <= std::min(jb, j1); k++) {
          RipupLevel& s = out[static_cast<size_t>(m.group) * n + (k - j0)];
          s = std::max(s, lv);
        }
      }
    }
  }
}

void FastGridBuilder::computeWiring(const int z,
                                    const int i0,
                                    const int i1,
                                    const int j0,
                                    const int j1,
                                    const int32_t ignored,
                                    std::vector<uint64_t>& words) const
{
  const TrackLayer& layer = graph_->layers[z];
  const int np = layer.numPositions();
  const int nt = layer.numTracks();
  const int n = j1 - j0 + 1;
  words.assign(static_cast<size_t>(std::max(0, i1 - i0 + 1)) * std::max(0, n),
               0);
  if (n <= 0 || i1 < i0) {
    return;
  }
  const int je = std::min(j1 + 1, np - 1);
  const int edges = je - j0;
  const bool jogs = layer.jogs;
  const int last = jogs && i1 + 1 < nt ? i1 + 1 : i1;
  thread_local std::vector<TrackLevels> levels;
  if (levels.size() < static_cast<size_t>(last - i0 + 1)) {
    levels.resize(last - i0 + 1);
  }
  for (int i = i0; i <= last; i++) {
    trackLevels(z, i, j0, je, ignored, levels[i - i0]);
  }
  const ShapeGrid& sg = checker_->grid();
  const ShapeClassTable& classes = sg.classes();
  thread_local std::vector<RipupLevel> jog;
  for (int i = i0; i <= i1; i++) {
    const TrackLevels& tl = levels[i - i0];
    uint64_t* out = &words[static_cast<size_t>(i - i0) * n];
    for (int j = j0; j <= j1; j++) {
      out[j - j0] = tl.words[j - j0];
    }
    // Preferred-direction edges: the swept wire is checked only where a
    // point level inside the edge exceeds both endpoints; every violation
    // of the swept wire is one of a point model on it.
    for (int j = j0; j <= j1 && j + 1 < np; j++) {
      for (int t = 0; t < fast_; t++) {
        const int shift = FastGrid::wiringShift(t, 0);
        const RipupLevel ends
            = std::max(FastGrid::field(tl.words[j - j0], shift),
                       FastGrid::field(tl.words[j + 1 - j0], shift));
        if (tl.interior[static_cast<size_t>(t) * edges + (j - j0)] <= ends) {
          continue;
        }
        const PlacedShape hull{
            prefHull(z, i, j, t), (*types_)[t].preferred[z].cls, 3};
        RipupLevel swept = 0;
        for (const GridPiece& p : tl.pieces.pieces) {
          if ((p.owner == ignored && p.owner != kNoNet && isOwned(p.type))
              || sg.level(p.type, p.owner) <= std::max(swept, ends)) {
            continue;
          }
          if (shapesViolate(
                  classes, hull, PlacedShape{p.rect, p.cls, p.line_end_ends})) {
            swept = sg.level(p.type, p.owner);
          }
        }
        if (swept > ends) {
          out[j - j0] |= uint64_t{1} << FastGrid::kWireEdgeBit;
          break;
        }
      }
    }
    // Jog edges (i, i + 1).
    if (jogs && i + 1 < nt) {
      const TrackLevels& above = levels[i + 1 - i0];
      jogLevels(z, i, j0, j1, ignored, jog);
      for (int j = j0; j <= j1; j++) {
        for (int t = 0; t < fast_; t++) {
          const int shift = FastGrid::wiringShift(t, 1);
          if (jog[static_cast<size_t>(t) * n + (j - j0)]
              > std::max(FastGrid::field(tl.words[j - j0], shift),
                         FastGrid::field(above.words[j - j0], shift))) {
            out[j - j0] |= uint64_t{1} << FastGrid::kJogEdgeBit;
            break;
          }
        }
      }
    }
  }
}

void FastGridBuilder::computeVia(const int z,
                                 const int i0,
                                 const int i1,
                                 const int c0,
                                 const int c1,
                                 const int32_t ignored,
                                 std::vector<uint64_t>& words) const
{
  const int n = c1 - c0 + 1;
  words.assign(static_cast<size_t>(std::max(0, i1 - i0 + 1)) * std::max(0, n),
               0);
  if (n <= 0 || i1 < i0 || z + 1 >= graph_->numLayers()) {
    return;
  }
  const TrackLayer& layer = graph_->layers[z];
  thread_local std::vector<LocatedModel> models;
  models.clear();
  viaModels(z, models);
  if (models.empty()) {
    return;
  }
  thread_local std::vector<Coord> coords;
  coords.resize(n);
  for (int c = c0; c <= c1; c++) {
    coords[c - c0] = layer.positions[grid_->viaPosition(z, c)];
  }
  const ShapeGrid& sg = checker_->grid();
  const ShapeClassTable& classes = sg.classes();
  thread_local GridWindow window;
  thread_local std::vector<LineSpan> spans;
  const int jv0 = grid_->viaPosition(z, c0);
  const int jv1 = grid_->viaPosition(z, c1);
  for (int i = i0; i <= i1; i++) {
    uint64_t* out = &words[static_cast<size_t>(i - i0) * n];
    const Point a = location(z, i, jv0);
    const Point b = location(z, i, jv1);
    for (const int layer_index : {viaGridLayer(z), viaGridLayer(z + 1)}) {
      if (layer_index == viaGridLayer(z + 1) && !has_projection_[z]) {
        continue;
      }
      const Rect& box = layer_index == viaGridLayer(z)
                            ? cut_model_box_[z]
                            : projection_model_box_[z];
      sg.collect(
          layer_index, windowOf(a, b, box, classes.reach(layer_index)), window);
      for (const GridPiece& p : window.pieces) {
        if (p.owner == ignored && p.owner != kNoNet && isOwned(p.type)) {
          continue;
        }
        const RipupLevel lv = sg.level(p.type, p.owner);
        const PlacedShape ps{p.rect, p.cls, p.line_end_ends};
        for (const LocatedModel& m : models) {
          if (m.layer != layer_index) {
            continue;
          }
          spans.clear();
          appendViolationSpans(
              classes, m.rect, m.cls, 3, a, layer.horizontal, ps, spans);
          const int shift = FastGrid::viaShift(m.group / 2, m.group % 2);
          for (const LineSpan& sp : spans) {
            const auto [ca, cb] = indexRange(coords, sp.lo, sp.hi);
            for (int k = std::max(ca, 0); k <= std::min(cb, n - 1); k++) {
              raise(out[k], shift, lv);
            }
          }
        }
      }
    }
  }
}

void FastGridBuilder::writeWiring(const int z,
                                  const int i0,
                                  const int i1,
                                  const int j0,
                                  const int j1,
                                  const std::vector<uint64_t>& words,
                                  const uint64_t mask) const
{
  const int n = j1 - j0 + 1;
  for (int i = i0; i <= i1 && n > 0; i++) {
    grid_->writeWiring(
        z, i, j0, &words[static_cast<size_t>(i - i0) * n], n, mask);
  }
}

void FastGridBuilder::writeVia(const int z,
                               const int i0,
                               const int i1,
                               const int c0,
                               const int c1,
                               const std::vector<uint64_t>& words) const
{
  const int n = c1 - c0 + 1;
  for (int i = i0; i <= i1 && n > 0; i++) {
    grid_->writeVia(
        z, i, c0, &words[static_cast<size_t>(i - i0) * n], n, ~uint64_t{0});
  }
}

void FastGridBuilder::rebuildVertices(const int z, const int i0, const int i1)
{
  const TrackLayer& layer = graph_->layers[z];
  const int np = layer.numPositions();
  // Blocks of positions bound the windows of the per-track computation.
  constexpr int kBlock = 256;
  thread_local std::vector<uint64_t> words;
  const uint64_t mask = ~(uint64_t{1} << FastGrid::kJogEdgeBit);
  for (int i = i0; i <= i1; i++) {
    for (int j = 0; j < np; j += kBlock) {
      const int j1 = std::min(np - 1, j + kBlock - 1);
      // Levels and preferred-edge bits; jog bits follow in the second phase.
      thread_local TrackLevels tl;
      const int je = std::min(j1 + 1, np - 1);
      trackLevels(z, i, j, je, kNoNet, tl);
      words.assign(j1 - j + 1, 0);
      const int edges = je - j;
      const ShapeGrid& sg = checker_->grid();
      for (int k = j; k <= j1; k++) {
        words[k - j] = tl.words[k - j];
        if (k + 1 >= np) {
          continue;
        }
        for (int t = 0; t < fast_; t++) {
          const int shift = FastGrid::wiringShift(t, 0);
          const RipupLevel ends
              = std::max(FastGrid::field(tl.words[k - j], shift),
                         FastGrid::field(tl.words[k + 1 - j], shift));
          if (tl.interior[static_cast<size_t>(t) * edges + (k - j)] <= ends) {
            continue;
          }
          const PlacedShape hull{
              prefHull(z, i, k, t), (*types_)[t].preferred[z].cls, 3};
          bool marked = false;
          for (const GridPiece& p : tl.pieces.pieces) {
            if (sg.level(p.type, p.owner) > ends
                && shapesViolate(sg.classes(),
                                 hull,
                                 PlacedShape{p.rect, p.cls, p.line_end_ends})) {
              marked = true;
              break;
            }
          }
          if (marked) {
            words[k - j] |= uint64_t{1} << FastGrid::kWireEdgeBit;
            break;
          }
        }
      }
      grid_->writeWiring(z, i, j, words.data(), j1 - j + 1, mask);
    }
  }
  if (z + 1 < graph_->numLayers() && grid_->numViaPositions(z) > 0) {
    const int nv = grid_->numViaPositions(z);
    for (int c = 0; c < nv; c += kBlock) {
      const int c1 = std::min(nv - 1, c + kBlock - 1);
      computeVia(z, i0, i1, c, c1, kNoNet, words);
      writeVia(z, i0, i1, c, c1, words);
    }
  }
}

void FastGridBuilder::computeJogBits(const int z,
                                     const int i0,
                                     const int i1,
                                     std::vector<JogBits>& out) const
{
  const TrackLayer& layer = graph_->layers[z];
  const int np = layer.numPositions();
  if (!layer.jogs || np == 0) {
    return;
  }
  thread_local std::vector<RipupLevel> jog;
  for (int p = i0; p <= std::min(i1, layer.numTracks() - 2); p++) {
    JogBits bits{z, p, std::vector<uint8_t>(np, 0)};
    constexpr int kBlock = 256;
    for (int j0 = 0; j0 < np; j0 += kBlock) {
      const int j1 = std::min(np - 1, j0 + kBlock - 1);
      const int n = j1 - j0 + 1;
      jogLevels(z, p, j0, j1, kNoNet, jog);
      for (int j = j0; j <= j1; j++) {
        const uint64_t wa = grid_->wiringWord(z, p, j);
        const uint64_t wb = grid_->wiringWord(z, p + 1, j);
        for (int t = 0; t < fast_; t++) {
          const int shift = FastGrid::wiringShift(t, 1);
          if (jog[static_cast<size_t>(t) * n + (j - j0)] > std::max(
                  FastGrid::field(wa, shift), FastGrid::field(wb, shift))) {
            bits.flags[j] = 1;
            break;
          }
        }
      }
    }
    out.push_back(std::move(bits));
  }
}

void FastGridBuilder::writeJogBits(const std::vector<JogBits>& bits) const
{
  thread_local std::vector<uint64_t> words;
  for (const JogBits& b : bits) {
    words.resize(b.flags.size());
    for (size_t k = 0; k < b.flags.size(); k++) {
      words[k] = b.flags[k] ? uint64_t{1} << FastGrid::kJogEdgeBit : 0;
    }
    grid_->writeWiring(b.z,
                       b.p,
                       0,
                       words.data(),
                       static_cast<int>(words.size()),
                       uint64_t{1} << FastGrid::kJogEdgeBit);
  }
}

void FastGridBuilder::rebuildJogBits(const int z, const int i0, const int i1)
{
  std::vector<JogBits> bits;
  computeJogBits(z, i0, i1, bits);
  writeJogBits(bits);
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

void FastGridBuilder::viaRange(const int z,
                               const int grid_layer,
                               const Rect& rect,
                               const Rect& box,
                               int& i0,
                               int& i1,
                               int& c0,
                               int& c1) const
{
  const TrackLayer& layer = graph_->layers[z];
  const int64_t reach = checker_->grid().classes().reach(grid_layer) + 1;
  const int64_t x_lo = int64_t{rect.xlo} - box.xhi - reach;
  const int64_t x_hi = int64_t{rect.xhi} - box.xlo + reach;
  const int64_t y_lo = int64_t{rect.ylo} - box.yhi - reach;
  const int64_t y_hi = int64_t{rect.yhi} - box.ylo + reach;
  std::tie(i0, i1) = layer.horizontal ? indexRange(layer.tracks, y_lo, y_hi)
                                      : indexRange(layer.tracks, x_lo, x_hi);
  const auto [j0, j1] = layer.horizontal
                            ? indexRange(layer.positions, x_lo, x_hi)
                            : indexRange(layer.positions, y_lo, y_hi);
  const int nv = grid_->numViaPositions(z);
  c0 = 0;
  c1 = -1;
  if (j0 > j1 || nv == 0) {
    return;
  }
  // Via positions are increasing positions.
  int lo = 0;
  int hi = nv;
  while (lo < hi) {
    const int mid = (lo + hi) / 2;
    if (grid_->viaPosition(z, mid) < j0) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  c0 = lo;
  c1 = c0 - 1;
  while (c1 + 1 < nv && grid_->viaPosition(z, c1 + 1) <= j1) {
    c1++;
  }
}

void FastGridBuilder::insertWiring(const int z,
                                   const int g,
                                   const PlacedShape& s,
                                   const RipupLevel lv)
{
  const TrackLayer& layer = graph_->layers[z];
  const int np = layer.numPositions();
  const int nt = layer.numTracks();
  const Range r = wiringRange(z, s.rect, g);
  const int tb0 = std::max(0, std::min(r.i0, r.p0));
  const int tb1 = std::min(nt - 1, std::max(r.i1, r.p1 + 1));
  const int jb0 = std::max(0, std::min(r.j0, r.e0));
  const int jb1 = std::min(np - 1, std::max(r.j1, r.e1 + 1));
  if (tb0 > tb1 || jb0 > jb1) {
    return;
  }
  const int n = jb1 - jb0 + 1;
  thread_local std::vector<uint64_t> before;
  thread_local std::vector<uint64_t> after;
  before.resize(static_cast<size_t>(tb1 - tb0 + 1) * n);
  for (int i = tb0; i <= tb1; i++) {
    for (int j = jb0; j <= jb1;) {
      int lo;
      int hi;
      const uint64_t w = grid_->wiringRun(z, i, j, lo, hi);
      for (; j <= std::min(hi, jb1); j++) {
        before[static_cast<size_t>(i - tb0) * n + (j - jb0)] = w;
      }
    }
  }
  after = before;
  auto at = [&](std::vector<uint64_t>& v, int i, int j) -> uint64_t& {
    return v[static_cast<size_t>(i - tb0) * n + (j - jb0)];
  };
  const ShapeClassTable& classes = checker_->grid().classes();
  thread_local std::vector<LocatedModel> models;
  models.clear();
  pointModels(z, models);
  thread_local std::vector<LineSpan> spans;
  // Point levels: max(old, the new shape's level where it violates).
  for (int i = r.i0; i <= r.i1; i++) {
    const Point a = location(z, i, std::max(r.j0, 0));
    for (const LocatedModel& m : models) {
      spans.clear();
      appendViolationSpans(
          classes, m.rect, m.cls, 3, a, layer.horizontal, s, spans);
      const int shift = FastGrid::wiringShift(m.group / 4, m.group % 4);
      for (const LineSpan& sp : spans) {
        const auto [ja, jb] = indexRange(layer.positions, sp.lo, sp.hi);
        for (int k = std::max(ja, r.j0); k <= std::min(jb, r.j1); k++) {
          raise(at(after, i, k), shift, lv);
        }
      }
    }
  }
  // Edge bits: with S' = max(S, lv if hit), M' >= M and S <= M where the bit
  // was clear, the bit becomes set exactly if hit and lv > M'; a set bit
  // stays set while M' == M, and is rechecked otherwise.
  auto decide = [&](bool old_bit,
                    const uint64_t wa0,
                    const uint64_t wb0,
                    const uint64_t wa1,
                    const uint64_t wb1,
                    const int kind,
                    auto&& hit,
                    auto&& exact) {
    bool changed = false;
    for (int t = 0; t < fast_; t++) {
      const int shift = FastGrid::wiringShift(t, kind);
      const RipupLevel m0
          = std::max(FastGrid::field(wa0, shift), FastGrid::field(wb0, shift));
      const RipupLevel m1
          = std::max(FastGrid::field(wa1, shift), FastGrid::field(wb1, shift));
      if (lv > m1 && hit(t)) {
        return true;
      }
      changed |= m0 != m1;
    }
    if (!old_bit) {
      return false;
    }
    return changed ? exact() : true;
  };
  for (int i = r.i0; i <= r.i1; i++) {
    for (int e = r.e0; e <= r.e1; e++) {
      const uint64_t bit = uint64_t{1} << FastGrid::kWireEdgeBit;
      const bool set = decide(
          (at(before, i, e) & bit) != 0,
          at(before, i, e),
          at(before, i, e + 1),
          at(after, i, e),
          at(after, i, e + 1),
          0,
          [&](int t) {
            return shapesViolate(
                classes,
                PlacedShape{
                    prefHull(z, i, e, t), (*types_)[t].preferred[z].cls, 3},
                s);
          },
          [&]() {
            const uint64_t wa = at(after, i, e);
            const uint64_t wb = at(after, i, e + 1);
            for (int t = 0; t < fast_; t++) {
              const int shift = FastGrid::wiringShift(t, 0);
              if (sweptLevel(z, i, e, t, kWireEdge, kNoNet) > std::max(
                      FastGrid::field(wa, shift), FastGrid::field(wb, shift))) {
                return true;
              }
            }
            return false;
          });
      at(after, i, e) = set ? at(after, i, e) | bit : at(after, i, e) & ~bit;
    }
  }
  if (layer.jogs) {
    for (int p = r.p0; p <= r.p1; p++) {
      const Point a = location(z, p, std::max(r.j0, 0));
      thread_local std::vector<std::vector<LineSpan>> hits;
      hits.assign(fast_, {});
      for (int t = 0; t < fast_; t++) {
        if ((*types_)[t].non_preferred[z].valid) {
          const LocatedModel m = jogHull(z, p, t);
          appendViolationSpans(
              classes, m.rect, m.cls, 3, a, layer.horizontal, s, hits[t]);
        }
      }
      for (int j = r.j0; j <= r.j1; j++) {
        const uint64_t bit = uint64_t{1} << FastGrid::kJogEdgeBit;
        const Coord c = layer.positions[j];
        const bool set
            = decide((at(before, p, j) & bit) != 0,
                     at(before, p, j),
                     at(before, p + 1, j),
                     at(after, p, j),
                     at(after, p + 1, j),
                     1,
                     [&](int t) {
                       for (const LineSpan& sp : hits[t]) {
                         if (sp.lo <= c && c <= sp.hi) {
                           return true;
                         }
                       }
                       return false;
                     },
                     [&]() {
                       const uint64_t wa = at(after, p, j);
                       const uint64_t wb = at(after, p + 1, j);
                       for (int t = 0; t < fast_; t++) {
                         const int shift = FastGrid::wiringShift(t, 1);
                         if (sweptLevel(z, p, j, t, kJogEdge, kNoNet)
                             > std::max(FastGrid::field(wa, shift),
                                        FastGrid::field(wb, shift))) {
                           return true;
                         }
                       }
                       return false;
                     });
        at(after, p, j) = set ? at(after, p, j) | bit : at(after, p, j) & ~bit;
      }
    }
  }
  for (int i = tb0; i <= tb1; i++) {
    const size_t off = static_cast<size_t>(i - tb0) * n;
    if (!std::equal(before.begin() + off,
                    before.begin() + off + n,
                    after.begin() + off)) {
      grid_->writeWiring(z, i, jb0, &after[off], n, ~uint64_t{0});
    }
  }
}

void FastGridBuilder::insertVia(const int z,
                                const int g,
                                const PlacedShape& s,
                                const RipupLevel lv,
                                const bool projection)
{
  // Via layer z data: cut models live on layer g = 2z + 1, projection
  // models of via layer z on layer 2z + 3.
  int i0, i1, c0, c1;
  viaRange(z,
           g,
           s.rect,
           projection ? projection_model_box_[z] : cut_model_box_[z],
           i0,
           i1,
           c0,
           c1);
  if (i0 > i1 || c0 > c1) {
    return;
  }
  const TrackLayer& layer = graph_->layers[z];
  const ShapeClassTable& classes = checker_->grid().classes();
  thread_local std::vector<LocatedModel> models;
  models.clear();
  viaModels(z, models);
  const int n = c1 - c0 + 1;
  thread_local std::vector<uint64_t> words;
  thread_local std::vector<LineSpan> spans;
  thread_local std::vector<Coord> coords;
  coords.resize(n);
  for (int c = c0; c <= c1; c++) {
    coords[c - c0] = layer.positions[grid_->viaPosition(z, c)];
  }
  for (int i = i0; i <= i1; i++) {
    words.resize(n);
    for (int c = c0; c <= c1; c++) {
      words[c - c0] = grid_->viaWord(z, i, c);
    }
    bool changed = false;
    const Point a = location(z, i, grid_->viaPosition(z, c0));
    for (const LocatedModel& m : models) {
      if (m.layer != g) {
        continue;
      }
      spans.clear();
      appendViolationSpans(
          classes, m.rect, m.cls, 3, a, layer.horizontal, s, spans);
      const int shift = FastGrid::viaShift(m.group / 2, m.group % 2);
      for (const LineSpan& sp : spans) {
        const auto [ca, cb] = indexRange(coords, sp.lo, sp.hi);
        for (int k = std::max(ca, 0); k <= std::min(cb, n - 1); k++) {
          const uint64_t old = words[k];
          raise(words[k], shift, lv);
          changed |= old != words[k];
        }
      }
    }
    if (changed) {
      grid_->writeVia(z, i, c0, words.data(), n, ~uint64_t{0});
    }
  }
}

void FastGridBuilder::shapeInserted(const int grid_layer,
                                    const PlacedShape& shape,
                                    const RipupLevel level)
{
  if (grid_layer % 2 == 0) {
    const int z = grid_layer / 2;
    if (z < graph_->numLayers()) {
      insertWiring(z, grid_layer, shape, level);
    }
    return;
  }
  const int z = (grid_layer - 1) / 2;
  if (z + 1 < graph_->numLayers()) {
    insertVia(z, grid_layer, shape, level, false);
  }
  if (z >= 1 && has_projection_[z - 1]) {
    insertVia(z - 1, grid_layer, shape, level, true);
  }
}

void FastGridBuilder::wiringCandidates(const int z,
                                       const int g,
                                       const PlacedShape& s,
                                       const RipupLevel lv,
                                       Candidates& out) const
{
  const TrackLayer& layer = graph_->layers[z];
  const Range r = wiringRange(z, s.rect, g);
  const ShapeClassTable& classes = checker_->grid().classes();
  thread_local std::vector<LocatedModel> models;
  models.clear();
  pointModels(z, models);
  thread_local std::vector<LineSpan> spans;
  // Point levels equal to the shape's level where it violates.
  for (int i = r.i0; i <= r.i1; i++) {
    const Point a = location(z, i, std::max(r.j0, 0));
    for (const LocatedModel& m : models) {
      spans.clear();
      appendViolationSpans(
          classes, m.rect, m.cls, 3, a, layer.horizontal, s, spans);
      const int shift = FastGrid::wiringShift(m.group / 4, m.group % 4);
      for (const LineSpan& sp : spans) {
        const auto [ja, jb] = indexRange(layer.positions, sp.lo, sp.hi);
        for (int k = std::max(ja, r.j0); k <= std::min(jb, r.j1); k++) {
          if (FastGrid::field(grid_->wiringWord(z, i, k), shift) == lv) {
            out.add(i, k);
          }
        }
      }
    }
    // Marked preferred edges whose swept wire the shape violates.
    for (int e = r.e0; e <= r.e1; e++) {
      if (!grid_->edgeQuery(z, i, e, kWireEdge)) {
        continue;
      }
      for (int t = 0; t < fast_; t++) {
        if (shapesViolate(
                classes,
                PlacedShape{
                    prefHull(z, i, e, t), (*types_)[t].preferred[z].cls, 3},
                s)) {
          out.add(i, e);
          break;
        }
      }
    }
  }
  if (!layer.jogs) {
    return;
  }
  for (int p = r.p0; p <= r.p1; p++) {
    const Point a = location(z, p, std::max(r.j0, 0));
    for (int t = 0; t < fast_; t++) {
      if (!(*types_)[t].non_preferred[z].valid) {
        continue;
      }
      const LocatedModel m = jogHull(z, p, t);
      spans.clear();
      appendViolationSpans(
          classes, m.rect, m.cls, 3, a, layer.horizontal, s, spans);
      for (const LineSpan& sp : spans) {
        const auto [ja, jb] = indexRange(layer.positions, sp.lo, sp.hi);
        for (int k = std::max(ja, r.j0); k <= std::min(jb, r.j1); k++) {
          if (grid_->edgeQuery(z, p, k, kJogEdge)) {
            out.add(p, k);
          }
        }
      }
    }
  }
}

void FastGridBuilder::viaCandidates(const int z,
                                    const int g,
                                    const PlacedShape& s,
                                    const RipupLevel lv,
                                    const bool projection,
                                    Candidates& out) const
{
  int i0, i1, c0, c1;
  viaRange(z,
           g,
           s.rect,
           projection ? projection_model_box_[z] : cut_model_box_[z],
           i0,
           i1,
           c0,
           c1);
  if (i0 > i1 || c0 > c1) {
    return;
  }
  const TrackLayer& layer = graph_->layers[z];
  const ShapeClassTable& classes = checker_->grid().classes();
  thread_local std::vector<LocatedModel> models;
  models.clear();
  viaModels(z, models);
  thread_local std::vector<LineSpan> spans;
  thread_local std::vector<Coord> coords;
  const int n = c1 - c0 + 1;
  coords.resize(n);
  for (int c = c0; c <= c1; c++) {
    coords[c - c0] = layer.positions[grid_->viaPosition(z, c)];
  }
  for (int i = i0; i <= i1; i++) {
    const Point a = location(z, i, grid_->viaPosition(z, c0));
    for (const LocatedModel& m : models) {
      if (m.layer != g) {
        continue;
      }
      spans.clear();
      appendViolationSpans(
          classes, m.rect, m.cls, 3, a, layer.horizontal, s, spans);
      const int shift = FastGrid::viaShift(m.group / 2, m.group % 2);
      for (const LineSpan& sp : spans) {
        const auto [ca, cb] = indexRange(coords, sp.lo, sp.hi);
        for (int k = std::max(ca, 0); k <= std::min(cb, n - 1); k++) {
          if (FastGrid::field(grid_->viaWord(z, i, c0 + k), shift) == lv) {
            out.add(i, c0 + k);
          }
        }
      }
    }
  }
}

void FastGridBuilder::recomputeCandidates(
    const int z,
    Candidates& cand,
    const int32_t ignored,
    std::vector<std::tuple<int, int, uint64_t>>& out) const
{
  const TrackLayer& layer = graph_->layers[z];
  const int np = layer.numPositions();
  const int nt = layer.numTracks();
  std::ranges::sort(cand.at);
  cand.at.erase(std::ranges::unique(cand.at).begin(), cand.at.end());
  struct Track
  {
    int i;
    int lo;
    int hi;
    std::vector<uint64_t> words;       // level fields, positions lo .. hi
    std::vector<RipupLevel> interior;  // [t * (n - 1) + (j - lo)]
    std::vector<uint32_t> near;        // pieces that can meet the track
  };
  thread_local std::vector<Track> tracks;
  tracks.clear();
  for (size_t k = 0; k < cand.at.size();) {
    size_t e = k;
    while (e + 1 < cand.at.size() && cand.at[e + 1].first == cand.at[k].first) {
      e++;
    }
    tracks.push_back({cand.at[k].first,
                      std::max(0, cand.at[k].second - 1),
                      std::min(np - 1, cand.at[e].second + 1),
                      {},
                      {},
                      {}});
    k = e + 1;
  }
  const ShapeGrid& sg = checker_->grid();
  const ShapeClassTable& classes = sg.classes();
  const int g = wiringGridLayer(z);
  const Coord reach = classes.reach(g);
  // One window for all tracks, including the neighbouring tracks of the
  // jog pairs.
  Rect window;
  bool first = true;
  for (const Track& t : tracks) {
    for (const int i : {t.i - 1, t.i, t.i + 1}) {
      if (i < 0 || i >= nt) {
        continue;
      }
      const Rect r = windowOf(location(z, i, t.lo),
                              location(z, i, t.hi),
                              wiring_model_box_[z],
                              reach);
      window = first ? r : window.hull(r);
      first = false;
    }
  }
  thread_local GridWindow pieces;
  sg.collect(g, window, pieces);
  thread_local std::vector<LocatedModel> models;
  models.clear();
  pointModels(z, models);
  const bool h = layer.horizontal;
  const Rect& box = wiring_model_box_[z];
  auto cross_meets = [&](const GridPiece& p, Coord lo, Coord hi) {
    const Coord p_lo = h ? p.rect.ylo : p.rect.xlo;
    const Coord p_hi = h ? p.rect.yhi : p.rect.xhi;
    return p_hi > lo - reach && p_lo < hi + reach;
  };
  auto skip = [&](const GridPiece& p) {
    return p.owner == ignored && p.owner != kNoNet && isOwned(p.type);
  };
  thread_local std::vector<LineSpan> spans;
  thread_local std::vector<uint32_t> near;
  for (Track& t : tracks) {
    const int n = t.hi - t.lo + 1;
    t.words.assign(n, 0);
    t.interior.assign(static_cast<size_t>(fast_) * std::max(0, n - 1), 0);
    const Coord c = layer.tracks[t.i];
    const Coord c_lo = c + (h ? box.ylo : box.xlo);
    const Coord c_hi = c + (h ? box.yhi : box.xhi);
    const Point a = location(z, t.i, t.lo);
    const Coord along_lo = layer.positions[t.lo];
    const Coord along_hi = layer.positions[t.hi];
    near.clear();
    for (const GridPiece& p : pieces.pieces) {
      const Coord pa_lo = h ? p.rect.xlo : p.rect.ylo;
      const Coord pa_hi = h ? p.rect.xhi : p.rect.yhi;
      if (skip(p) || !cross_meets(p, c_lo, c_hi)
          || pa_lo >= along_hi + (h ? box.xhi : box.yhi) + reach
          || pa_hi <= along_lo + (h ? box.xlo : box.ylo) - reach) {
        continue;
      }
      near.push_back(static_cast<uint32_t>(&p - pieces.pieces.data()));
    }
    t.near = near;
    for (const uint32_t pi : near) {
      const GridPiece& p = pieces.pieces[pi];
      const RipupLevel lv = sg.level(p.type, p.owner);
      const PlacedShape ps{p.rect, p.cls, p.line_end_ends};
      for (const LocatedModel& m : models) {
        if (farApart(classes,
                     m,
                     p,
                     h,
                     c + (h ? m.rect.ylo : m.rect.xlo),
                     c + (h ? m.rect.yhi : m.rect.xhi),
                     along_lo + (h ? m.rect.xlo : m.rect.ylo),
                     along_hi + (h ? m.rect.xhi : m.rect.yhi))) {
          continue;
        }
        spans.clear();
        appendViolationSpans(classes, m.rect, m.cls, 3, a, h, ps, spans);
        const int shift = FastGrid::wiringShift(m.group / 4, m.group % 4);
        for (const LineSpan& sp : spans) {
          const auto [ja, jb] = indexRange(layer.positions, sp.lo, sp.hi);
          for (int k = std::max(ja, t.lo); k <= std::min(jb, t.hi); k++) {
            raise(t.words[k - t.lo], shift, lv);
          }
          if (m.group % 4 == 0 && n > 1) {
            RipupLevel* in
                = &t.interior[static_cast<size_t>(m.group / 4) * (n - 1)];
            for (int e = std::max(t.lo, ja - 1); e <= std::min(t.hi - 1, jb);
                 e++) {
              in[e - t.lo] = std::max(in[e - t.lo], lv);
            }
          }
        }
      }
    }
  }
  auto find_track = [&](int i) -> const Track* {
    const auto it = std::ranges::lower_bound(tracks, i, {}, &Track::i);
    return it != tracks.end() && it->i == i ? &*it : nullptr;
  };
  constexpr uint64_t kWireBit = uint64_t{1} << FastGrid::kWireEdgeBit;
  constexpr uint64_t kJogBit = uint64_t{1} << FastGrid::kJogEdgeBit;
  auto levels_at = [&](int i, int j) -> uint64_t {
    const Track* t = find_track(i);
    if (t != nullptr && j >= t->lo && j <= t->hi) {
      return t->words[j - t->lo];
    }
    return grid_->wiringWord(z, i, j) & ~FastGrid::kEdgeBits;
  };
  // Jog bits to recompute: pairs below and above every candidate.
  thread_local std::vector<std::pair<int, int>> jog;
  jog.clear();
  if (layer.jogs) {
    for (const auto& [i, j] : cand.at) {
      if (i + 1 < nt) {
        jog.emplace_back(i, j);
      }
      if (i >= 1) {
        jog.emplace_back(i - 1, j);
      }
    }
    std::ranges::sort(jog);
    jog.erase(std::ranges::unique(jog).begin(), jog.end());
  }
  thread_local std::vector<std::pair<std::pair<int, int>, bool>> jog_bits;
  jog_bits.clear();
  thread_local std::vector<RipupLevel> swept;
  for (size_t k = 0; k < jog.size();) {
    const int p = jog[k].first;
    size_t e = k;
    while (e + 1 < jog.size() && jog[e + 1].first == p) {
      e++;
    }
    const int m = static_cast<int>(e - k + 1);
    swept.assign(static_cast<size_t>(fast_) * m, 0);
    const Coord c = layer.tracks[p];
    const Coord c2 = layer.tracks[p + 1];
    const Point a = location(z, p, jog[k].second);
    thread_local std::vector<std::pair<int, LocatedModel>> hulls;
    hulls.clear();
    for (int t = 0; t < fast_; t++) {
      if ((*types_)[t].non_preferred[z].valid) {
        hulls.emplace_back(t, jogHull(z, p, t));
      }
    }
    for (const GridPiece& piece : pieces.pieces) {
      if (skip(piece)
          || !cross_meets(piece,
                          c + (h ? box.ylo : box.xlo),
                          c2 + (h ? box.yhi : box.xhi))) {
        continue;
      }
      const RipupLevel lv = sg.level(piece.type, piece.owner);
      const PlacedShape ps{piece.rect, piece.cls, piece.line_end_ends};
      for (const auto& [t, hull] : hulls) {
        if (farApart(classes,
                     hull,
                     piece,
                     h,
                     c + (h ? hull.rect.ylo : hull.rect.xlo),
                     c + (h ? hull.rect.yhi : hull.rect.xhi),
                     layer.positions[jog[k].second]
                         + (h ? hull.rect.xlo : hull.rect.ylo),
                     layer.positions[jog[e].second]
                         + (h ? hull.rect.xhi : hull.rect.yhi))) {
          continue;
        }
        spans.clear();
        appendViolationSpans(classes, hull.rect, hull.cls, 3, a, h, ps, spans);
        for (const LineSpan& sp : spans) {
          for (int q = 0; q < m; q++) {
            const Coord x = layer.positions[jog[k + q].second];
            if (sp.lo <= x && x <= sp.hi) {
              RipupLevel& s = swept[static_cast<size_t>(t) * m + q];
              s = std::max(s, lv);
            }
          }
        }
      }
    }
    for (int q = 0; q < m; q++) {
      const int j = jog[k + q].second;
      const uint64_t wa = levels_at(p, j);
      const uint64_t wb = levels_at(p + 1, j);
      bool bit = false;
      for (int t = 0; t < fast_ && !bit; t++) {
        const int shift = FastGrid::wiringShift(t, 1);
        bit = swept[static_cast<size_t>(t) * m + q] > std::max(
                  FastGrid::field(wa, shift), FastGrid::field(wb, shift));
      }
      jog_bits.push_back({{p, j}, bit});
    }
    k = e + 1;
  }
  auto jog_bit = [&](int p, int j, bool stored) {
    const auto it = std::ranges::lower_bound(
        jog_bits, std::make_pair(p, j), {}, [](const auto& x) {
          return x.first;
        });
    return it != jog_bits.end() && it->first == std::make_pair(p, j)
               ? it->second
               : stored;
  };
  // Words of the candidate tracks: levels, preferred bits of the edges
  // with a changed endpoint, jog bits.
  for (const Track& t : tracks) {
    const int n = t.hi - t.lo + 1;
    for (int j = t.lo; j <= t.hi; j++) {
      const uint64_t stored = grid_->wiringWord(z, t.i, j);
      uint64_t w = t.words[j - t.lo];
      bool wire_bit = (stored & kWireBit) != 0;
      if (j < t.hi) {
        wire_bit = false;
        for (int ty = 0; ty < fast_ && !wire_bit; ty++) {
          const int shift = FastGrid::wiringShift(ty, 0);
          const RipupLevel ends
              = std::max(FastGrid::field(t.words[j - t.lo], shift),
                         FastGrid::field(t.words[j + 1 - t.lo], shift));
          if (t.interior[static_cast<size_t>(ty) * (n - 1) + (j - t.lo)]
              <= ends) {
            continue;
          }
          const PlacedShape hull{
              prefHull(z, t.i, j, ty), (*types_)[ty].preferred[z].cls, 3};
          for (const uint32_t pi : t.near) {
            const GridPiece& p = pieces.pieces[pi];
            if (sg.level(p.type, p.owner) > ends
                && shapesViolate(classes,
                                 hull,
                                 PlacedShape{p.rect, p.cls, p.line_end_ends})) {
              wire_bit = true;
              break;
            }
          }
        }
      }
      if (wire_bit) {
        w |= kWireBit;
      }
      if (jog_bit(t.i, j, (stored & kJogBit) != 0)) {
        w |= kJogBit;
      }
      out.emplace_back(t.i, j, w);
    }
  }
  // Jog bits of tracks without own candidates.
  for (const auto& [pj, bit] : jog_bits) {
    const auto [p, j] = pj;
    const Track* t = find_track(p);
    if (t != nullptr && j >= t->lo && j <= t->hi) {
      continue;
    }
    const uint64_t stored = grid_->wiringWord(z, p, j);
    out.emplace_back(p, j, bit ? stored | kJogBit : stored & ~kJogBit);
  }
}

void FastGridBuilder::shapeRemoved(const int grid_layer,
                                   const PlacedShape& shape,
                                   const RipupLevel level)
{
  thread_local std::vector<uint64_t> words;
  if (grid_layer % 2 == 0) {
    const int z = grid_layer / 2;
    if (z >= graph_->numLayers()) {
      return;
    }
    Candidates c;
    wiringCandidates(z, grid_layer, shape, level, c);
    if (c.empty()) {
      return;
    }
    thread_local std::vector<std::tuple<int, int, uint64_t>> changed;
    changed.clear();
    recomputeCandidates(z, c, kNoNet, changed);
    writeChanged(z, changed);
    return;
  }
  const int z = (grid_layer - 1) / 2;
  for (const bool projection : {false, true}) {
    const int vz = projection ? z - 1 : z;
    if (vz < 0 || vz + 1 >= graph_->numLayers()
        || (projection && !has_projection_[vz])) {
      continue;
    }
    Candidates c;
    viaCandidates(vz, grid_layer, shape, level, projection, c);
    if (!c.empty()) {
      computeVia(vz, c.i0, c.i1, c.j0, c.j1, kNoNet, words);
      writeVia(vz, c.i0, c.i1, c.j0, c.j1, words);
    }
  }
}

void FastGridBuilder::writeChanged(
    const int z,
    std::vector<std::tuple<int, int, uint64_t>>& changed) const
{
  // Runs of consecutive positions per track.
  thread_local std::vector<uint64_t> words;
  std::ranges::sort(changed);
  for (size_t k = 0; k < changed.size();) {
    size_t e = k;
    words.assign(1, std::get<2>(changed[k]));
    while (e + 1 < changed.size()
           && std::get<0>(changed[e + 1]) == std::get<0>(changed[k])
           && std::get<1>(changed[e + 1]) == std::get<1>(changed[e]) + 1) {
      e++;
      words.push_back(std::get<2>(changed[e]));
    }
    grid_->writeWiring(z,
                       std::get<0>(changed[k]),
                       std::get<1>(changed[k]),
                       words.data(),
                       static_cast<int>(words.size()),
                       ~uint64_t{0});
    k = e + 1;
  }
}

void FastGridBuilder::shapeRemoved(const int grid_layer,
                                   const PlacedShape& shape,
                                   const RipupLevel level,
                                   PendingRemovals& pending)
{
  if (grid_layer % 2 == 1) {
    shapeRemoved(grid_layer, shape, level);
    return;
  }
  const int z = grid_layer / 2;
  if (z >= graph_->numLayers()) {
    return;
  }
  if (pending.wiring.size() < static_cast<size_t>(graph_->numLayers())) {
    pending.wiring.resize(graph_->numLayers());
  }
  // Stored levels are unchanged until flush(), so the candidates of every
  // removed shape are found against the levels before the group.
  wiringCandidates(z, grid_layer, shape, level, pending.wiring[z]);
}

void FastGridBuilder::flush(PendingRemovals& pending)
{
  thread_local std::vector<std::tuple<int, int, uint64_t>> changed;
  for (size_t z = 0; z < pending.wiring.size(); z++) {
    if (pending.wiring[z].empty()) {
      continue;
    }
    changed.clear();
    recomputeCandidates(
        static_cast<int>(z), pending.wiring[z], kNoNet, changed);
    writeChanged(static_cast<int>(z), changed);
    pending.wiring[z] = Candidates{};
  }
}

void FastGridBuilder::shapeChanged(const int grid_layer, const Rect& rect)
{
  thread_local std::vector<uint64_t> words;
  if (grid_layer % 2 == 0) {
    const int z = grid_layer / 2;
    if (z >= graph_->numLayers()) {
      return;
    }
    const Range r = wiringRange(z, rect, grid_layer);
    const int nt = graph_->layers[z].numTracks();
    const int np = graph_->layers[z].numPositions();
    const int i0 = std::max(0, std::min(r.i0, r.p0));
    const int i1 = std::min(nt - 1, std::max(r.i1, r.p1));
    const int j0 = std::max(0, std::min(r.j0, r.e0));
    const int j1 = std::min(np - 1, std::max(r.j1, r.e1));
    if (i0 <= i1 && j0 <= j1) {
      computeWiring(z, i0, i1, j0, j1, kNoNet, words);
      writeWiring(z, i0, i1, j0, j1, words, ~uint64_t{0});
    }
    return;
  }
  const int z = (grid_layer - 1) / 2;
  for (const bool projection : {false, true}) {
    const int vz = projection ? z - 1 : z;
    if (vz < 0 || vz + 1 >= graph_->numLayers()
        || (projection && !has_projection_[vz])) {
      continue;
    }
    int i0, i1, c0, c1;
    viaRange(vz,
             grid_layer,
             rect,
             projection ? projection_model_box_[vz] : cut_model_box_[vz],
             i0,
             i1,
             c0,
             c1);
    if (i0 <= i1 && c0 <= c1) {
      computeVia(vz, i0, i1, c0, c1, kNoNet, words);
      writeVia(vz, i0, i1, c0, c1, words);
    }
  }
}

void FastGridBuilder::buildOverlay(const int32_t owner,
                                   const std::vector<OwnShape>& shapes,
                                   FastGridOverlay& overlay) const
{
  overlay.clear();
  overlay.owner = owner;
  // Positions that can change when the owner's shapes are removed, as for
  // a removal, recomputed exactly with the owner's shapes ignored.
  thread_local std::vector<uint64_t> words;
  thread_local std::map<uint64_t, std::map<int, uint64_t>> wiring;
  thread_local std::map<uint64_t, std::map<int, uint64_t>> via;
  wiring.clear();
  via.clear();
  thread_local std::vector<Candidates> per_layer;
  per_layer.assign(graph_->numLayers(), Candidates{});
  for (const OwnShape& own : shapes) {
    const int g = own.grid_layer;
    if (g % 2 == 0) {
      const int z = g / 2;
      if (z >= graph_->numLayers()) {
        continue;
      }
      wiringCandidates(z, g, own.shape, own.level, per_layer[z]);
      continue;
    }
    const int z = (g - 1) / 2;
    for (const bool projection : {false, true}) {
      const int vz = projection ? z - 1 : z;
      if (vz < 0 || vz + 1 >= graph_->numLayers()
          || (projection && !has_projection_[vz])) {
        continue;
      }
      Candidates c;
      viaCandidates(vz, g, own.shape, own.level, projection, c);
      if (c.empty()) {
        continue;
      }
      computeVia(vz, c.i0, c.i1, c.j0, c.j1, owner, words);
      const int n = c.j1 - c.j0 + 1;
      for (int i = c.i0; i <= c.i1; i++) {
        auto& track = via[FastGridOverlay::key(vz, i)];
        for (int k = c.j0; k <= c.j1; k++) {
          track[k] = words[static_cast<size_t>(i - c.i0) * n + (k - c.j0)];
        }
      }
    }
  }
  // All candidates of a layer in one recomputation.
  thread_local std::vector<std::tuple<int, int, uint64_t>> changed;
  for (int z = 0; z < graph_->numLayers(); z++) {
    if (per_layer[z].empty()) {
      continue;
    }
    changed.clear();
    recomputeCandidates(z, per_layer[z], owner, changed);
    for (const auto& [i, j, w] : changed) {
      wiring[FastGridOverlay::key(z, i)][j] = w;
    }
  }
  auto runs =
      [](const std::map<uint64_t, std::map<int, uint64_t>>& in,
         std::unordered_map<uint64_t, std::vector<FastGridOverlay::Run>>& out) {
        for (const auto& [key, track] : in) {
          std::vector<FastGridOverlay::Run>& r = out[key];
          for (const auto& [j, w] : track) {
            if (!r.empty() && r.back().hi + 1 == j && r.back().word == w) {
              r.back().hi = j;
            } else {
              r.push_back({j, j, w});
            }
          }
        }
      };
  runs(wiring, overlay.wiring);
  runs(via, overlay.via);
}

RipupLevel FastGridBuilder::vertexLevel(const int z,
                                        const int i,
                                        const int j,
                                        const int wire_type,
                                        const Kind kind,
                                        const int32_t ignored_owner) const
{
  return exactLevel(z, i, j, wire_type, kind, ignored_owner);
}

RipupLevel FastGridBuilder::edgeLevel(const int z,
                                      const int i,
                                      const int j,
                                      const int wire_type,
                                      const Kind kind,
                                      const int32_t ignored_owner) const
{
  return exactLevel(z, i, j, wire_type, kind, ignored_owner);
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
