// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Speculative tail pipeline (BrTailPipeline.h): random instances on a toy
// routing space of tiles. A search reads the tiles of its net's box, a
// decision detects conflicts with earlier commits of its batch and writes
// the tiles of the box. The pipeline on 1-8 threads (random task jitter,
// random cost estimates) must equal a serial execution of the semantics:
// final tiles, per-net decision sequences, failed nets and counts. Also run
// it with -fsanitize=thread.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <map>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

#include "br/BrTailPipeline.h"

using drt::br::TailPipeline;
using Box = TailPipeline::Box;
using Decision = TailPipeline::Decision;

namespace {

struct Instance
{
  int tiles_x;
  int tiles_y;
  std::vector<Box> boxes;
  std::vector<int> order;
  int batch;
  int retry_delay;
  uint64_t salt;
};

struct Outcome
{
  std::vector<uint64_t> tiles;
  std::map<int, std::vector<int>> decisions;  // per net, in attempt order
  std::vector<int> failed;
  int64_t committed = 0;
  int64_t conflicts = 0;
};

uint64_t mix(uint64_t h, const uint64_t v)
{
  h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  return h * 0xff51afd7ed558ccdULL;
}

bool meets(const Box& a, const Box& b)
{
  return a[0] <= b[2] && b[0] <= a[2] && a[1] <= b[3] && b[1] <= a[3];
}

// The route of a net: a function of the tiles of its box only.
uint64_t searchRoute(const Instance& in,
                     const std::vector<uint64_t>& tiles,
                     const int net,
                     const int attempt)
{
  const Box& b = in.boxes[net];
  uint64_t h = mix(in.salt, net * 7919 + attempt);
  for (int x = b[0]; x <= b[2]; x++) {
    for (int y = b[1]; y <= b[3]; y++) {
      h = mix(h, tiles[static_cast<size_t>(x) * in.tiles_y + y]);
    }
  }
  return h;
}

bool routeFound(const uint64_t route)
{
  return route % 11 != 0;
}

bool routesConflict(const uint64_t a, const uint64_t b)
{
  return mix(a, b) % 4 == 0;
}

void commitRoute(const Instance& in,
                 std::vector<uint64_t>& tiles,
                 const int net,
                 const uint64_t route)
{
  const Box& b = in.boxes[net];
  for (int x = b[0]; x <= b[2]; x++) {
    for (int y = b[1]; y <= b[3]; y++) {
      uint64_t& t = tiles[static_cast<size_t>(x) * in.tiles_y + y];
      t = mix(t, route + net);
    }
  }
}

// Serial execution of the semantics.
Outcome reference(const Instance& in)
{
  Outcome out;
  out.tiles.assign(static_cast<size_t>(in.tiles_x) * in.tiles_y, 1);
  struct Attempt
  {
    int net;
    int attempt;
  };
  std::vector<std::vector<Attempt>> batches;
  for (size_t i = 0; i < in.order.size(); i++) {
    const size_t b = i / in.batch;
    if (batches.size() <= b) {
      batches.resize(b + 1);
    }
    batches[b].push_back({in.order[i], 0});
  }
  std::vector<uint64_t> route(in.boxes.size(), 0);
  std::vector<char> found(in.boxes.size(), 0);
  for (size_t b = 0; b < batches.size(); b++) {
    const std::vector<uint64_t> snapshot = out.tiles;
    for (const Attempt& a : batches[b]) {
      route[a.net] = searchRoute(in, snapshot, a.net, a.attempt);
      found[a.net] = routeFound(route[a.net]);
    }
    std::vector<int> committed;
    for (size_t k = 0; k < batches[b].size(); k++) {
      const Attempt a = batches[b][k];
      Decision d = Decision::kCommitted;
      if (!found[a.net]) {
        d = Decision::kFailed;
        out.failed.push_back(a.net);
      } else {
        for (const int e : committed) {
          if (meets(in.boxes[e], in.boxes[a.net])
              && routesConflict(route[a.net], route[e])) {
            d = Decision::kConflict;
          }
        }
        if (d == Decision::kConflict) {
          const size_t t = b + in.retry_delay;
          if (batches.size() <= t) {
            batches.resize(t + 1);
          }
          batches[t].push_back({a.net, a.attempt + 1});
          out.conflicts++;
        } else {
          commitRoute(in, out.tiles, a.net, route[a.net]);
          committed.push_back(a.net);
          out.committed++;
        }
      }
      out.decisions[a.net].push_back(static_cast<int>(d));
    }
  }
  return out;
}

Outcome pipelined(const Instance& in, const int threads, const uint64_t seed)
{
  Outcome out;
  out.tiles.assign(static_cast<size_t>(in.tiles_x) * in.tiles_y, 1);
  std::vector<uint64_t> route(in.boxes.size(), 0);
  std::vector<int> attempt(in.boxes.size(), 0);
  std::mutex log_mutex;
  std::vector<std::mt19937_64> jitter;
  for (int t = 0; t < threads; t++) {
    jitter.emplace_back(seed * 131 + t);
  }
  auto pause = [&](const int thread) {
    const int us = static_cast<int>(jitter[thread]() % 200);
    if (us < 60) {
      std::this_thread::sleep_for(std::chrono::microseconds(us));
    }
  };
  std::mt19937_64 cost_rng(seed);
  std::vector<double> cost(in.boxes.size());
  for (double& c : cost) {
    c = static_cast<double>(cost_rng() % 1000);
  }
  TailPipeline::Hooks hooks;
  hooks.parallel = [](const int count, const std::function<void(int)>& body) {
    std::vector<std::thread> pool;
    for (int t = 0; t < count; t++) {
      pool.emplace_back(body, t);
    }
    for (std::thread& t : pool) {
      t.join();
    }
  };
  hooks.search = [&](const int net, const int thread) {
    pause(thread);
    route[net] = searchRoute(in, out.tiles, net, attempt[net]);
    pause(thread);
    return routeFound(route[net]);
  };
  hooks.decide = [&](const int net,
                     const bool found,
                     const std::vector<int>& earlier,
                     const int thread) {
    pause(thread);
    Decision d = Decision::kCommitted;
    if (!found) {
      d = Decision::kFailed;
    } else {
      for (const int e : earlier) {
        if (routesConflict(route[net], route[e])) {
          d = Decision::kConflict;
        }
      }
      if (d == Decision::kCommitted) {
        commitRoute(in, out.tiles, net, route[net]);
      } else {
        attempt[net]++;
      }
    }
    {
      std::lock_guard<std::mutex> lock(log_mutex);
      out.decisions[net].push_back(static_cast<int>(d));
    }
    return d;
  };
  hooks.cost = [&](const int net) { return cost[net]; };
  TailPipeline pipeline(in.batch, in.retry_delay);
  const TailPipeline::Stats stats
      = pipeline.run(in.order, in.boxes, threads, hooks, out.failed);
  out.committed = stats.committed;
  out.conflicts = stats.conflicts;
  return out;
}

Instance randomInstance(std::mt19937_64& rng)
{
  auto uniform = [&](const int lo, const int hi) {
    return lo + static_cast<int>(rng() % static_cast<uint64_t>(hi - lo + 1));
  };
  Instance in;
  in.tiles_x = uniform(3, 24);
  in.tiles_y = uniform(3, 24);
  const int nets = uniform(1, 220);
  const int max_side = uniform(1, std::max(in.tiles_x, in.tiles_y));
  for (int n = 0; n < nets; n++) {
    const int w = uniform(1, std::min(max_side, in.tiles_x));
    const int h = uniform(1, std::min(max_side, in.tiles_y));
    const int x = uniform(0, in.tiles_x - w);
    const int y = uniform(0, in.tiles_y - h);
    in.boxes.push_back({x, y, x + w - 1, y + h - 1});
  }
  // Every net once, in random order; some nets left out.
  for (int n = 0; n < nets; n++) {
    if (rng() % 8 != 0) {
      in.order.push_back(n);
    }
  }
  std::shuffle(in.order.begin(), in.order.end(), rng);
  in.batch = uniform(1, 40);
  in.retry_delay = uniform(1, 3);
  in.salt = rng();
  return in;
}

}  // namespace

int main(int argc, char** argv)
{
  const int cases = argc > 1 ? std::atoi(argv[1]) : 300;
  std::mt19937_64 rng(20261111);
  int failures = 0;
  int64_t attempts = 0;
  int64_t conflicts = 0;
  for (int c = 0; c < cases; c++) {
    const Instance in = randomInstance(rng);
    const Outcome ref = reference(in);
    for (const int threads : {1, 2, 3, 8}) {
      const Outcome got = pipelined(in, threads, rng());
      const bool equal
          = got.tiles == ref.tiles && got.decisions == ref.decisions
            && got.failed == ref.failed && got.committed == ref.committed
            && got.conflicts == ref.conflicts;
      if (!equal) {
        failures++;
        std::printf(
            "case %d threads %d: mismatch (nets %zu, batch %d, delay %d; "
            "committed %lld vs %lld, conflicts %lld vs %lld)\n",
            c,
            threads,
            in.order.size(),
            in.batch,
            in.retry_delay,
            static_cast<long long>(got.committed),
            static_cast<long long>(ref.committed),
            static_cast<long long>(got.conflicts),
            static_cast<long long>(ref.conflicts));
      }
    }
    attempts += ref.committed + ref.conflicts
                + static_cast<int64_t>(ref.failed.size());
    conflicts += ref.conflicts;
  }
  std::printf(
      "tail pipeline: %d cases x 4 thread counts, %lld attempts, %lld "
      "conflicts, %d failures\n",
      cases,
      static_cast<long long>(attempts),
      static_cast<long long>(conflicts),
      failures);
  return failures == 0 ? 0 : 1;
}
