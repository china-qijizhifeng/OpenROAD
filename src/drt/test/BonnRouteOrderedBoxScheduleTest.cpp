// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Ordered box schedule (BrOrderedBoxSchedule.h): random instances on a toy
// routing space of tiles. A task reads the tiles of its box and writes a
// value derived from them and its index to every tile of the box. On 1-8
// threads with random task jitter the final tiles and the value each task
// read must equal the serial execution in index order. Every tile also
// records the task working on it: a second task entering a tile in use is a
// violation even when the results happen to agree. Also run it with
// -fsanitize=thread.
//
// Also: disjointLayerOrder() against a definition-based first-fit layering
// (same order, chain at most the number of layers); retries in the same run
// (task n + i repeats task i with a larger box when task i asks for it, as
// the rip-up boxes do with ripup_box_chunk_retry) against the serial order;
// the busy-thread time line adds up to the wall time and never shows more
// busy threads than threads; the critical path is at least the longest task.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <memory>
#include <random>
#include <thread>
#include <vector>

#include "br/BrOrderedBoxSchedule.h"

using drt::br::OrderedBoxSchedule;
using Box = OrderedBoxSchedule::Box;

namespace {

uint64_t mix(uint64_t x)
{
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ULL;
  x ^= x >> 33;
  return x;
}

struct Space
{
  int tiles_x;
  int tiles_y;
  std::vector<uint64_t> tiles;
  std::unique_ptr<std::atomic<int>[]> user;
  std::atomic<int64_t> violations{0};
  Space(const int x, const int y)
      : tiles_x(x), tiles_y(y), tiles(x * y, 0), user(new std::atomic<int>[x * y])
  {
    for (int i = 0; i < x * y; i++) {
      user[i] = -1;
    }
  }
  // The work of task `index`: returns the value it read.
  uint64_t work(const Box& box, const int index, const bool check, const int jitter_us)
  {
    if (check) {
      for (int y = box[1]; y <= box[3]; y++) {
        for (int x = box[0]; x <= box[2]; x++) {
          int expected = -1;
          if (!user[y * tiles_x + x].compare_exchange_strong(expected, index)) {
            violations++;
          }
        }
      }
    }
    uint64_t read = 0;
    for (int y = box[1]; y <= box[3]; y++) {
      for (int x = box[0]; x <= box[2]; x++) {
        read = mix(read ^ tiles[y * tiles_x + x]);
      }
    }
    if (jitter_us > 0) {
      std::this_thread::sleep_for(std::chrono::microseconds(jitter_us));
    }
    const uint64_t value = mix(read + static_cast<uint64_t>(index) + 1);
    for (int y = box[1]; y <= box[3]; y++) {
      for (int x = box[0]; x <= box[2]; x++) {
        tiles[y * tiles_x + x] = value ^ static_cast<uint64_t>(y * tiles_x + x);
      }
    }
    if (check) {
      for (int y = box[1]; y <= box[3]; y++) {
        for (int x = box[0]; x <= box[2]; x++) {
          int expected = index;
          if (!user[y * tiles_x + x].compare_exchange_strong(expected, -1)) {
            violations++;
          }
        }
      }
    }
    return read;
  }
};

// By definition: some tile lies in both boxes.
bool shareTile(const Box& a, const Box& b)
{
  for (int y = a[1]; y <= a[3]; y++) {
    for (int x = a[0]; x <= a[2]; x++) {
      if (x >= b[0] && x <= b[2] && y >= b[1] && y <= b[3]) {
        return true;
      }
    }
  }
  return false;
}

void parallel(const int threads, const std::function<void(int)>& body)
{
  std::vector<std::thread> pool;
  for (int t = 1; t < threads; t++) {
    pool.emplace_back(body, t);
  }
  body(0);
  for (std::thread& thread : pool) {
    thread.join();
  }
}

// First-fit layering by definition: task j joins the first layer in which
// no box shares a tile with box j.
std::vector<int> referenceLayerOrder(const std::vector<Box>& boxes, int& layer_count)
{
  std::vector<std::vector<int>> layers;
  for (int j = 0; j < static_cast<int>(boxes.size()); j++) {
    size_t l = 0;
    for (; l < layers.size(); l++) {
      bool free = true;
      for (const int i : layers[l]) {
        free = free && !shareTile(boxes[i], boxes[j]);
      }
      if (free) {
        break;
      }
    }
    if (l == layers.size()) {
      layers.emplace_back();
    }
    layers[l].push_back(j);
  }
  layer_count = static_cast<int>(layers.size());
  std::vector<int> order;
  for (const auto& layer : layers) {
    order.insert(order.end(), layer.begin(), layer.end());
  }
  return order;
}

int64_t longestChain(const std::vector<Box>& boxes)
{
  std::vector<int> chain(boxes.size(), 1);
  int64_t longest = 0;
  for (size_t j = 0; j < boxes.size(); j++) {
    for (size_t i = 0; i < j; i++) {
      if (shareTile(boxes[i], boxes[j])) {
        chain[j] = std::max(chain[j], chain[i] + 1);
      }
    }
    longest = std::max<int64_t>(longest, chain[j]);
  }
  return longest;
}

// Retries in the same run: tasks 0 .. n-1 with the boxes, tasks n .. 2n-1
// with the boxes widened by `margin` tiles; task n + i works only when task
// i read a value that asks for a retry. Returns the number of mismatches
// against the serial order over 1, 2, 3 and 8 threads.
int64_t retryCase(std::mt19937_64& rng,
                  const int tiles_x,
                  const int tiles_y,
                  const std::vector<Box>& boxes,
                  int64_t& retried,
                  int64_t& violations)
{
  const int count = static_cast<int>(boxes.size());
  std::vector<Box> all = boxes;
  const int margin = 1 + static_cast<int>(rng() % 3);
  for (const Box& b : boxes) {
    all.push_back({std::max(b[0] - margin, 0),
                   std::max(b[1] - margin, 0),
                   std::min(b[2] + margin, tiles_x - 1),
                   std::min(b[3] + margin, tiles_y - 1)});
  }
  auto asks = [](const uint64_t read) { return read % 3 == 0; };
  Space serial(tiles_x, tiles_y);
  std::vector<uint64_t> serial_read(2 * count, 0);
  std::vector<char> serial_ran(2 * count, 0);
  for (int i = 0; i < 2 * count; i++) {
    if (i >= count && !asks(serial_read[i - count])) {
      continue;
    }
    serial_read[i] = serial.work(all[i], i, false, 0);
    serial_ran[i] = 1;
    retried += i >= count;
  }
  int64_t mismatches = 0;
  for (const int threads : {1, 2, 3, 8}) {
    Space space(tiles_x, tiles_y);
    std::vector<uint64_t> read(2 * count, 0);
    std::vector<char> ran(2 * count, 0);
    std::vector<int> jitter(2 * count, 0);
    for (int& j : jitter) {
      j = static_cast<int>(rng() % 2) * 50;
    }
    OrderedBoxSchedule::run(all, threads, parallel, [&](const int index, const int) {
      // Task i returned before task n + i starts (their boxes meet).
      if (index >= count && !asks(read[index - count])) {
        return;
      }
      read[index] = space.work(all[index], index, true, jitter[index]);
      ran[index] = 1;
    });
    if (space.tiles != serial.tiles || read != serial_read || ran != serial_ran) {
      mismatches++;
    }
    violations += space.violations;
  }
  return mismatches;
}

}  // namespace

int main()
{
  std::mt19937_64 rng(20261113);
  int64_t instances = 0;
  int64_t mismatches = 0;
  int64_t violations = 0;
  int64_t chain_errors = 0;
  int64_t concurrent_instances = 0;
  int64_t layer_errors = 0;
  int64_t layered_shorter = 0;
  int64_t line_errors = 0;
  int64_t retry_mismatches = 0;
  int64_t retried = 0;
  for (int round = 0; round < 400; round++) {
    const int tiles_x = 4 + static_cast<int>(rng() % 28);
    const int tiles_y = 4 + static_cast<int>(rng() % 28);
    const int count = static_cast<int>(rng() % 97);
    // Mostly small boxes; some instances hold boxes over the whole space.
    const int max_span = 1 + static_cast<int>(rng() % 6);
    const bool whole = rng() % 4 == 0;
    std::vector<Box> boxes;
    for (int i = 0; i < count; i++) {
      if (whole && rng() % 16 == 0) {
        boxes.push_back({0, 0, tiles_x - 1, tiles_y - 1});
        continue;
      }
      const int x = static_cast<int>(rng() % tiles_x);
      const int y = static_cast<int>(rng() % tiles_y);
      const int w = static_cast<int>(rng() % max_span);
      const int h = static_cast<int>(rng() % max_span);
      boxes.push_back({x, y, std::min(x + w, tiles_x - 1), std::min(y + h, tiles_y - 1)});
    }
    Space serial(tiles_x, tiles_y);
    std::vector<uint64_t> serial_read(count, 0);
    for (int i = 0; i < count; i++) {
      serial_read[i] = serial.work(boxes[i], i, false, 0);
    }
    // Longest chain by definition (a schedule with more waits than the
    // boxes require differs here).
    std::vector<int> chain(count, 1);
    int64_t longest = 0;
    for (int j = 0; j < count; j++) {
      for (int i = 0; i < j; i++) {
        if (shareTile(boxes[i], boxes[j])) {
          chain[j] = std::max(chain[j], chain[i] + 1);
        }
      }
      longest = std::max<int64_t>(longest, chain[j]);
    }
    {
      int layer_count = 0;
      const std::vector<int> expected = referenceLayerOrder(boxes, layer_count);
      const std::vector<int> order = OrderedBoxSchedule::disjointLayerOrder(boxes);
      std::vector<Box> permuted;
      for (const int i : order) {
        permuted.push_back(boxes[i]);
      }
      const int64_t layered = longestChain(permuted);
      if (order != expected || layered > layer_count || layered > longest) {
        layer_errors++;
      }
      layered_shorter += layered < longest;
    }
    if (round % 4 == 0) {
      retry_mismatches += retryCase(rng, tiles_x, tiles_y, boxes, retried, violations);
    }
    for (const int threads : {1, 2, 3, 8}) {
      Space space(tiles_x, tiles_y);
      std::vector<uint64_t> read(count, 0);
      std::vector<int> jitter(count, 0);
      for (int i = 0; i < count; i++) {
        jitter[i] = rng() % 3 == 0 ? static_cast<int>(rng() % 300) : 0;
      }
      std::atomic<int> running{0};
      std::atomic<int> peak{0};
      const OrderedBoxSchedule::Stats stats = OrderedBoxSchedule::run(
          boxes, threads, parallel, [&](const int index, const int thread) {
            const int now = ++running;
            int seen = peak.load();
            while (now > seen && !peak.compare_exchange_weak(seen, now)) {
            }
            if (thread < 0 || thread >= threads) {
              space.violations++;
            }
            read[index] = space.work(boxes[index], index, true, jitter[index]);
            --running;
          });
      instances++;
      if (space.tiles != serial.tiles || read != serial_read) {
        mismatches++;
      }
      violations += space.violations;
      if (stats.tasks != count || stats.longest_chain != longest || peak > threads) {
        chain_errors++;
      }
      double line_total = 0;
      for (int c = 0; c < drt::br::BusyThreads::kClasses; c++) {
        line_total += stats.busy_threads[c];
        if (stats.busy_threads[c] > 0 && c > drt::br::BusyThreads::classOf(threads)) {
          line_errors++;
        }
      }
      if (std::abs(line_total - stats.wall) > 1e-6 + 1e-9 * stats.wall
          || stats.critical_path > stats.wall + 1e-6
          || (count > 0 && stats.critical_path <= 0 && stats.busy > 0)) {
        line_errors++;
      }
      if (peak > 1) {
        concurrent_instances++;
      }
    }
  }
  std::printf("ordered box schedule: %lld runs, %lld mismatches against the serial order, %lld tile violations, "
              "%lld statistics errors, %lld runs with concurrent tasks\n",
              static_cast<long long>(instances),
              static_cast<long long>(mismatches),
              static_cast<long long>(violations),
              static_cast<long long>(chain_errors),
              static_cast<long long>(concurrent_instances));
  std::printf("disjoint layer order: %lld errors, %lld instances with a shorter chain; busy-thread line / critical "
              "path: %lld errors; retries in the same run: %lld retries, %lld mismatches against the serial order\n",
              static_cast<long long>(layer_errors),
              static_cast<long long>(layered_shorter),
              static_cast<long long>(line_errors),
              static_cast<long long>(retried),
              static_cast<long long>(retry_mismatches));
  const bool ok = mismatches == 0 && violations == 0 && chain_errors == 0 && concurrent_instances > 100
                  && layer_errors == 0 && layered_shorter > 10 && line_errors == 0 && retry_mismatches == 0
                  && retried > 100;
  std::printf(ok ? "PASS\n" : "FAIL\n");
  return ok ? 0 : 1;
}
