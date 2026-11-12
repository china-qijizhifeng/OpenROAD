// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Speculative tail pipeline (BrTailPipeline.h): random instances on a toy
// routing space of tiles. A search reads the tiles of its net's box, a
// decision detects conflicts with earlier commits of its batch and writes
// the tiles of the box. The pipeline on 1-8 threads (random task jitter,
// random cost estimates) must equal a serial execution of the semantics:
// final tiles, per-net decision sequences, failed nets and counts. Also run
// it with -fsanitize=thread.
//
// Tail lookahead (TailLookahead, design D4): random instances against a
// serial execution of the barrier semantics S, the main pass with an
// effort budget followed by the deferred pass, as routeSpeculative runs
// them. Compared item by item: every attempt (batch, position, net,
// decision), the commits each search sees (tiles of its box and the set of
// committed nets whose boxes meet it), the failed and deferred lists in
// order, and the final tiles. Each tile also counts the searches reading it
// and the commits writing it: a commit during a search of a meeting box is
// a violation even when the results happen to agree.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include "br/BrTailPipeline.h"

using drt::br::TailLookahead;
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


// ---------------------------------------------------------------------------
// Tail lookahead against the barrier semantics S.

struct LookaheadInstance
{
  int tiles_x;
  int tiles_y;
  std::vector<Box> boxes;
  std::vector<int> order;
  int batch;
  uint64_t salt;
  int found_modulus;     // a route is missing when route % found_modulus == 0
  int conflict_modulus;  // routes conflict when mix % conflict_modulus == 0
  int over_permille;     // share of missing routes over the effort budget
  int long_permille;     // share of searches that sleep long
};

// (batch, position, net, decision) of every attempt of one pass.
using AttemptList = std::vector<std::tuple<int, int, int, int>>;

struct LookaheadOutcome
{
  std::vector<uint64_t> tiles;
  AttemptList main_pass;
  AttemptList deferred_pass;
  std::vector<int> failed;
  std::vector<int> deferred;
  // Per (net, attempt): hash of the committed nets that meet the net's box,
  // as seen by the search.
  std::map<std::pair<int, int>, uint64_t> visible;
  int64_t violations = 0;
  int64_t early = 0;
};

bool lookaheadFound(const LookaheadInstance& in, const uint64_t route)
{
  return route % static_cast<uint64_t>(in.found_modulus) != 0;
}

bool lookaheadConflict(const LookaheadInstance& in,
                       const uint64_t a,
                       const uint64_t b)
{
  return mix(a, b) % static_cast<uint64_t>(in.conflict_modulus) == 0;
}

bool overBudget(const LookaheadInstance& in, const int net, const int attempt)
{
  return mix(in.salt ^ 0x5bd1e995ULL, net * 131 + attempt) % 1000
         < static_cast<uint64_t>(in.over_permille);
}

uint64_t lookaheadRoute(const LookaheadInstance& in,
                        const std::vector<uint64_t>& tiles,
                        const int net,
                        const int attempt,
                        const bool budget)
{
  const Box& b = in.boxes[net];
  uint64_t h = mix(in.salt, net * 7919 + attempt * 2 + (budget ? 1 : 0));
  for (int x = b[0]; x <= b[2]; x++) {
    for (int y = b[1]; y <= b[3]; y++) {
      h = mix(h, tiles[static_cast<size_t>(x) * in.tiles_y + y]);
    }
  }
  return h;
}

void lookaheadCommit(const LookaheadInstance& in,
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

template <class Committed>
uint64_t visibleCommits(const LookaheadInstance& in,
                        const int net,
                        const Committed& committed)
{
  uint64_t h = 7;
  for (int m = 0; m < static_cast<int>(in.boxes.size()); m++) {
    if (m != net && meets(in.boxes[m], in.boxes[net]) && committed(m)) {
      h = mix(h, m);
    }
  }
  return h;
}

// Serial S, as routeSpeculativeBatches: the main pass with the budget, then
// the deferred nets in a pass without one.
LookaheadOutcome lookaheadReference(const LookaheadInstance& in)
{
  LookaheadOutcome out;
  out.tiles.assign(static_cast<size_t>(in.tiles_x) * in.tiles_y, 1);
  std::vector<int> attempt(in.boxes.size(), 0);
  std::vector<char> committed(in.boxes.size(), 0);
  std::vector<uint64_t> route(in.boxes.size(), 0);
  auto pass = [&](const std::vector<int>& order,
                  const bool budget,
                  AttemptList& list,
                  std::vector<int>& deferred) {
    std::deque<int> queue(order.begin(), order.end());
    int b = 0;
    while (!queue.empty()) {
      std::vector<int> batch;
      while (!queue.empty() && static_cast<int>(batch.size()) < in.batch) {
        batch.push_back(queue.front());
        queue.pop_front();
      }
      std::vector<char> found(batch.size(), 0);
      for (size_t k = 0; k < batch.size(); k++) {
        const int n = batch[k];
        route[n] = lookaheadRoute(in, out.tiles, n, attempt[n], budget);
        found[k] = lookaheadFound(in, route[n]);
        out.visible[{n, attempt[n]}] = visibleCommits(
            in, n, [&](const int m) { return committed[m] != 0; });
      }
      std::vector<int> retry;
      std::vector<int> earlier;
      for (size_t k = 0; k < batch.size(); k++) {
        const int n = batch[k];
        Decision d = Decision::kCommitted;
        if (!found[k]) {
          d = Decision::kFailed;
          if (budget && overBudget(in, n, attempt[n])) {
            deferred.push_back(n);
          } else {
            out.failed.push_back(n);
          }
        } else {
          for (const int e : earlier) {
            if (meets(in.boxes[e], in.boxes[n])
                && lookaheadConflict(in, route[n], route[e])) {
              d = Decision::kConflict;
            }
          }
          if (d == Decision::kConflict) {
            retry.push_back(n);
          } else {
            lookaheadCommit(in, out.tiles, n, route[n]);
            committed[n] = 1;
            earlier.push_back(n);
          }
        }
        attempt[n]++;
        list.emplace_back(b, static_cast<int>(k), n, static_cast<int>(d));
      }
      for (auto it = retry.rbegin(); it != retry.rend(); ++it) {
        queue.push_front(*it);
      }
      b++;
    }
  };
  pass(in.order, true, out.main_pass, out.deferred);
  std::vector<int> none;
  pass(out.deferred, false, out.deferred_pass, none);
  return out;
}

LookaheadOutcome lookaheadScheduled(const LookaheadInstance& in,
                                    const int threads,
                                    const uint64_t seed)
{
  LookaheadOutcome out;
  const size_t tile_count = static_cast<size_t>(in.tiles_x) * in.tiles_y;
  out.tiles.assign(tile_count, 1);
  std::vector<int> attempt(in.boxes.size(), 0);
  std::vector<uint64_t> route(in.boxes.size(), 0);
  std::unique_ptr<std::atomic<int>[]> committed(
      new std::atomic<int>[in.boxes.size()]);
  std::unique_ptr<std::atomic<int>[]> readers(new std::atomic<int>[tile_count]);
  std::unique_ptr<std::atomic<int>[]> writers(new std::atomic<int>[tile_count]);
  for (size_t n = 0; n < in.boxes.size(); n++) {
    committed[n] = 0;
  }
  for (size_t t = 0; t < tile_count; t++) {
    readers[t] = 0;
    writers[t] = 0;
  }
  std::atomic<int64_t> violations{0};
  std::mutex log_mutex;
  std::vector<std::mt19937_64> jitter;
  for (int t = 0; t < threads; t++) {
    jitter.emplace_back(seed * 131 + t);
  }
  auto pause = [&](const int thread, const bool long_task) {
    const uint64_t r = jitter[thread]();
    if (long_task) {
      std::this_thread::sleep_for(std::chrono::microseconds(300 + r % 900));
    } else if (r % 8 == 0) {
      std::this_thread::sleep_for(std::chrono::microseconds(r % 40));
    } else if (r % 8 == 1) {
      std::this_thread::yield();
    }
  };
  // Marks the tiles of a box: delta to `mine`, and counts a violation when
  // the other counter is set.
  auto mark = [&](const int net,
                  std::atomic<int>* mine,
                  const std::atomic<int>* other,
                  const int delta) {
    const Box& b = in.boxes[net];
    for (int x = b[0]; x <= b[2]; x++) {
      for (int y = b[1]; y <= b[3]; y++) {
        const size_t t = static_cast<size_t>(x) * in.tiles_y + y;
        mine[t] += delta;
        if (other[t].load() != 0) {
          violations++;
        }
      }
    }
  };
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
  bool budget = true;
  hooks.search = [&](const int net, const int thread) {
    mark(net, readers.get(), writers.get(), 1);
    const bool long_task
        = mix(in.salt ^ 77, net * 31 + attempt[net]) % 1000
          < static_cast<uint64_t>(in.long_permille);
    pause(thread, long_task);
    route[net] = lookaheadRoute(in, out.tiles, net, attempt[net], budget);
    const uint64_t seen = visibleCommits(
        in, net, [&](const int m) { return committed[m].load() != 0; });
    {
      std::lock_guard<std::mutex> lock(log_mutex);
      out.visible[{net, attempt[net]}] = seen;
    }
    pause(thread, false);
    mark(net, readers.get(), writers.get(), -1);
    return lookaheadFound(in, route[net]);
  };
  hooks.decide = [&](const int net,
                     const bool found,
                     const std::vector<int>& earlier,
                     const int thread) {
    pause(thread, false);
    Decision d = Decision::kCommitted;
    if (!found) {
      d = Decision::kFailed;
    } else {
      for (const int e : earlier) {
        if (lookaheadConflict(in, route[net], route[e])) {
          d = Decision::kConflict;
        }
      }
      if (d == Decision::kCommitted) {
        mark(net, writers.get(), readers.get(), 1);
        lookaheadCommit(in, out.tiles, net, route[net]);
        committed[net] = 1;
        pause(thread, false);
        mark(net, writers.get(), readers.get(), -1);
      }
    }
    attempt[net]++;
    return d;
  };
  std::mt19937_64 cost_rng(seed);
  std::vector<double> cost(in.boxes.size());
  for (double& c : cost) {
    c = static_cast<double>(cost_rng() % 1000);
  }
  hooks.cost = [&](const int net) { return cost[net]; };
  auto pass = [&](const std::vector<int>& order,
                  AttemptList& list,
                  std::vector<int>& deferred) {
    TailLookahead lookahead(in.batch);
    std::vector<int> unrouted;
    std::vector<TailLookahead::Attempt> attempts;
    const TailLookahead::Stats stats
        = lookahead.run(order, in.boxes, threads, hooks, unrouted, &attempts);
    out.early += stats.overlapped_searches;
    for (const TailLookahead::Attempt& a : attempts) {
      list.emplace_back(
          a.batch, a.position, a.net, static_cast<int>(a.decision));
    }
    // As routeSpeculativePipeline: unrouted nets over the budget are
    // deferred, the others failed; the attempt counter already advanced.
    for (const int n : unrouted) {
      if (budget && overBudget(in, n, attempt[n] - 1)) {
        deferred.push_back(n);
      } else {
        out.failed.push_back(n);
      }
    }
  };
  pass(in.order, out.main_pass, out.deferred);
  budget = false;
  std::vector<int> none;
  pass(out.deferred, out.deferred_pass, none);
  out.violations = violations.load();
  return out;
}

LookaheadInstance randomLookaheadInstance(std::mt19937_64& rng)
{
  auto uniform = [&](const int lo, const int hi) {
    return lo + static_cast<int>(rng() % static_cast<uint64_t>(hi - lo + 1));
  };
  LookaheadInstance in;
  in.tiles_x = uniform(2, 28);
  in.tiles_y = uniform(2, 28);
  const int nets = uniform(1, 200);
  const int max_side = uniform(1, std::max(in.tiles_x, in.tiles_y));
  for (int n = 0; n < nets; n++) {
    const int w = uniform(1, std::min(max_side, in.tiles_x));
    const int h = uniform(1, std::min(max_side, in.tiles_y));
    const int x = uniform(0, in.tiles_x - w);
    const int y = uniform(0, in.tiles_y - h);
    in.boxes.push_back({x, y, x + w - 1, y + h - 1});
  }
  for (int n = 0; n < nets; n++) {
    if (rng() % 8 != 0) {
      in.order.push_back(n);
    }
  }
  std::shuffle(in.order.begin(), in.order.end(), rng);
  in.batch = uniform(1, 40);
  in.salt = rng();
  in.found_modulus = uniform(2, 16);
  in.conflict_modulus = uniform(1, 8);
  in.over_permille = uniform(0, 3) == 0 ? 0 : uniform(0, 1000);
  in.long_permille = uniform(0, 40);
  return in;
}

// Returns the number of failing runs.
int runLookaheadCases(const int cases, int64_t& attempts, int64_t& early)
{
  std::mt19937_64 rng(20261113);
  int failures = 0;
  for (int c = 0; c < cases; c++) {
    const LookaheadInstance in = randomLookaheadInstance(rng);
    const LookaheadOutcome ref = lookaheadReference(in);
    // Thread counts 1-8, each in a quarter of the cases.
    for (const int threads : {1 + c % 8, 1 + (c + 4) % 8}) {
      const LookaheadOutcome got = lookaheadScheduled(in, threads, rng());
      const bool equal
          = got.tiles == ref.tiles && got.main_pass == ref.main_pass
            && got.deferred_pass == ref.deferred_pass
            && got.failed == ref.failed && got.deferred == ref.deferred
            && got.visible == ref.visible && got.violations == 0;
      early += got.early;
      if (!equal) {
        failures++;
        std::printf(
            "lookahead case %d threads %d: mismatch (nets %zu, batch %d; "
            "attempts %zu+%zu vs %zu+%zu, failed %zu vs %zu, deferred %zu "
            "vs %zu, visible %s, tiles %s, violations %lld)\n",
            c,
            threads,
            in.order.size(),
            in.batch,
            got.main_pass.size(),
            got.deferred_pass.size(),
            ref.main_pass.size(),
            ref.deferred_pass.size(),
            got.failed.size(),
            ref.failed.size(),
            got.deferred.size(),
            ref.deferred.size(),
            got.visible == ref.visible ? "same" : "DIFF",
            got.tiles == ref.tiles ? "same" : "DIFF",
            static_cast<long long>(got.violations));
      }
    }
    attempts += static_cast<int64_t>(ref.main_pass.size()
                                     + ref.deferred_pass.size());
  }
  return failures;
}

}  // namespace

int main(int argc, char** argv)
{
  const int cases = argc > 1 ? std::atoi(argv[1]) : 300;
  const int lookahead_cases = argc > 2 ? std::atoi(argv[2]) : 10000;
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
  int64_t lookahead_attempts = 0;
  int64_t early = 0;
  const int lookahead_failures
      = runLookaheadCases(lookahead_cases, lookahead_attempts, early);
  std::printf(
      "tail lookahead: %d cases x 2 thread counts (1-8), %lld attempts, "
      "%lld searches started early, %d failures\n",
      lookahead_cases,
      static_cast<long long>(lookahead_attempts),
      static_cast<long long>(early),
      lookahead_failures);
  // The lookahead must actually start searches early on these instances.
  const bool exercised = lookahead_cases == 0 || early > 0;
  return failures == 0 && lookahead_failures == 0 && exercised ? 0 : 1;
}
