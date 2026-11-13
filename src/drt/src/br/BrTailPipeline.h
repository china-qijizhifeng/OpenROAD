// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

// Deterministic pipeline for the speculative tail of BonnRoute DR
// (Sec. 5.1 extension, profile 2013-improved): nets are searched in
// parallel on a committed routing state that is defined by the net order
// alone, while earlier nets are still being committed elsewhere.
//
// Semantics (independent of threads and timing):
//  * The nets of `order` form fixed batches of `batch` nets ("fresh"
//    attempts). A net whose route conflicts at its decision is appended as
//    a new attempt to the end of batch b + retry_delay.
//  * An attempt of batch b is searched on the routing state holding
//    exactly the commits of the attempts of batches < b.
//  * The attempts of a batch are decided in list order. A decision sees
//    the commits of the attempts earlier in its batch, so it can detect a
//    conflict with them (`earlier`), and then commits or rejects the route.
//
// Each net has a box (inclusive global routing tiles) that holds all its
// reads and writes of the routing space; attempts with disjoint boxes do not
// influence each other. The pipeline uses this to run them concurrently:
//  * a search of an attempt in batch b waits for the decisions of the
//    attempts of earlier batches whose boxes meet its box;
//  * a decision waits for the searches of the other attempts of its batch
//    whose boxes meet its box and for the decisions of those earlier in the
//    batch.
// Batch b's member list is final once batch b - 1 is final and batch
// b - retry_delay is fully decided; searches of fresh attempts of batch b
// start once batch b - 1 is final, decisions once batch b is final.
// Every search and decision thereby reads exactly the state the semantics
// define, so results equal a serial execution of the semantics.

#include <array>
#include <condition_variable>
#include <algorithm>
#include <cstdint>
#include <functional>
#include <mutex>
#include <queue>
#include <tuple>
#include <utility>
#include <optional>
#include <stdexcept>
#include <vector>

#include "br/BrBusyThreads.h"

namespace drt::br {

class TailPipeline
{
 public:
  // xlo, ylo, xhi, yhi in global routing tiles, inclusive.
  using Box = std::array<int, 4>;
  enum class Decision
  {
    kCommitted,
    kConflict,
    kFailed
  };
  struct Hooks
  {
    // Runs body(thread) for thread = 0 .. threads - 1 concurrently and
    // returns when all of them returned.
    std::function<void(int threads, const std::function<void(int)>& body)>
        parallel;
    // Searches a route for the net on thread `thread`; true when found.
    std::function<bool(int net, int thread)> search;
    // Decides the attempt: `earlier` lists the nets committed earlier in
    // the same batch whose boxes meet the net's box, in list order.
    std::function<Decision(int net,
                           bool found,
                           const std::vector<int>& earlier,
                           int thread)>
        decide;
    // Optional estimated search cost: within a batch, searches with a
    // larger estimate start first (any order gives the same results).
    std::function<double(int net)> cost;
  };
  struct Stats
  {
    int64_t batches = 0;
    int64_t attempts = 0;
    int64_t committed = 0;
    int64_t conflicts = 0;
    int64_t failed = 0;
    double wall = 0;
    double search_busy = 0;
    double decide_busy = 0;
    // Thread seconds spent waiting for a ready search or decision.
    double idle = 0;
    // Searches that started while an earlier batch still had undecided
    // attempts (overlap gained over a barrier per batch).
    int64_t overlapped_searches = 0;
    // TailLookahead: smallest and largest batch (attempts, retries
    // included).
    int64_t smallest_batch = 0;
    int64_t largest_batch = 0;
    // TailLookahead: wall seconds by the number of threads inside a search
    // or decision.
    BusyThreads::Seconds busy_threads{};
  };

  TailPipeline(int batch, int retry_delay)
      : batch_(batch < 1 ? 1 : batch),
        retry_delay_(retry_delay < 1 ? 1 : retry_delay)
  {
  }

  // boxes: indexed by net id. failed receives the nets without a route in
  // (batch, position) order.
  Stats run(const std::vector<int>& order,
            const std::vector<Box>& boxes,
            int threads,
            const Hooks& hooks,
            std::vector<int>& failed);

 private:
  struct Item
  {
    int net;
    int batch;
    int position;
    bool fresh;
    bool searched = false;
    bool decided = false;
    bool search_started = false;
    bool decide_started = false;
    bool found = false;
    Decision decision = Decision::kFailed;
    int search_deps = -1;  // -1: not computed yet
    int decide_deps = -1;
    double cost = 0;
    std::vector<int> on_searched;  // items whose decide_deps count it
    std::vector<int> on_decided;   // items whose search/decide deps count it
    std::vector<char> on_decided_kind;  // 0 search, 1 decide
  };
  struct Batch
  {
    std::vector<int> items;
    int decided = 0;
    // Conflicting attempts of batch - retry_delay: (source position, net).
    std::vector<std::pair<int, int>> retries;
  };

  static bool meets(const Box& a, const Box& b)
  {
    return a[0] <= b[2] && b[0] <= a[2] && a[1] <= b[3] && b[1] <= a[3];
  }
  bool meets(int a, int b) const
  {
    return meets((*boxes_)[items_[a].net], (*boxes_)[items_[b].net]);
  }
  void addItem(int net, int batch, bool fresh);
  void computeSearchDeps(int item);
  void computeDecideDeps(int item);
  void finalizeBatches();
  void readySearch(int item);
  void readyDecide(int item);
  void searched(int item);
  void decided(int item);
  bool finished() const;

  const int batch_;
  const int retry_delay_;
  const std::vector<Box>* boxes_ = nullptr;
  const Hooks* hooks_ = nullptr;
  std::vector<Item> items_;
  std::vector<Batch> batches_;
  int final_ = 0;  // batches < final_ have their final member list
  int low_ = 0;    // batches < low_ are fully decided
  int open_items_ = 0;
  int running_ = 0;
  bool deadlock_ = false;
  // Ready queues: decisions by (batch, position); searches by (batch,
  // -cost, position).
  std::priority_queue<std::pair<int64_t, int>,
                      std::vector<std::pair<int64_t, int>>,
                      std::greater<>>
      ready_decide_;
  std::priority_queue<std::tuple<int, double, int, int>,
                      std::vector<std::tuple<int, double, int, int>>,
                      std::greater<>>
      ready_search_;
  std::mutex mutex_;
  std::condition_variable wake_;
};

// The speculative tail with the batch semantics of the barrier version
// (BonnRouteDR::routeSpeculativeBatches), where idle threads search nets of
// the next batch early (design D4, tail lookahead). Results equal a serial
// execution of these semantics for every thread count and timing.
//
// Semantics (S):
//  * The queue Q starts as `order`. Batch k takes the first
//    min(batch, |Q|) nets of Q.
//  * Every attempt of batch k is searched on the routing state holding the
//    commits of batches < k.
//  * The attempts of batch k are decided in list order. A decision sees the
//    commits earlier in the batch (`earlier`), and commits, rejects
//    (conflict) or fails the attempt.
//  * The conflicting nets of batch k go back to the front of Q in list
//    order, so batch k + 1 holds them first, followed by fresh nets of Q.
//
// Schedule (each net's box holds all its reads and writes of the routing
// space; boxes "meet" when they share a tile):
//  * R1: batch k + 1 is formed once every attempt of batch k is decided.
//  * R2: the decision of y waits for y's search, for the searches of the
//    other attempts of its batch whose boxes meet y's box, and for the
//    decisions of the attempts earlier in the batch whose boxes meet it.
//  * R3: while batch k is decided, let r be the attempts of batch k decided
//    as conflicts and U the attempts of batch k that are undecided and not
//    known to be without a route (searching, or found and undecided). The
//    p-th net after batch k in Q (p from 0) starts its search when
//    (a) p < batch - r - |U| and (b) no box of U meets its box. Such a
//    search has the lowest priority, candidates in increasing p.
//  * R4: that search is the net's search in batch k + 1.
// Batch k produces at most r + |U| conflicts, so by (a) the net belongs to
// batch k + 1. By (b), no attempt of batch k writes its box after the search
// starts, so the search reads the state of the commits of batches <= k.
// D6: size of the speculative batches. Batch k holds every retry of batch
// k-1 first, then fresh nets of the queue up to size(k) in total. size(0) =
// size(1) = the start size (start, clamped to [min, max]; 0: max);
// size(k+1) follows the conflict share of batch k-1 (the last batch decided
// when batch k is formed, so the size of the next batch is known while
// batch k is decided, as TailLookahead's R3 needs): above high_percent
// conflicts per attempt it halves, below low_percent it doubles, within
// [min, max]. min 0 or >= max: every batch has max nets. A start below max
// lets the batches grow while conflicts are rare (tail_batch_max).
struct TailBatchPolicy
{
  int max = 32;
  int min = 0;
  int high_percent = 50;
  int low_percent = 20;
  int start = 0;
  bool adaptive() const { return min > 0 && min < max; }
  int startSize() const
  {
    if (!adaptive() || start <= 0) {
      return max;
    }
    return std::clamp(start, min, max);
  }
  int next(int size, int64_t attempts, int64_t conflicts) const
  {
    if (!adaptive() || attempts <= 0) {
      return adaptive() ? size : max;
    }
    if (conflicts * 100 > static_cast<int64_t>(high_percent) * attempts) {
      return std::max(min, size / 2);
    }
    if (conflicts * 100 < static_cast<int64_t>(low_percent) * attempts) {
      return std::min(max, size * 2);
    }
    return size;
  }
};

// Sizes of the batches of one pass under a TailBatchPolicy.
class TailBatchSizes
{
 public:
  explicit TailBatchSizes(const TailBatchPolicy& policy)
      : policy_(policy), max_(std::max(1, policy.max))
  {
  }
  // size(k) needs the outcomes of the batches before k - 1.
  int size(const int k)
  {
    while (static_cast<int>(sizes_.size()) <= k) {
      const int j = static_cast<int>(sizes_.size());
      if (!policy_.adaptive()) {
        sizes_.push_back(max_);
        continue;
      }
      if (j < 2) {
        sizes_.push_back(policy_.startSize());
        continue;
      }
      if (static_cast<int>(outcomes_.size()) < j - 1) {
        throw std::logic_error("TailBatchSizes: outcome of batch missing");
      }
      const auto [attempts, conflicts] = outcomes_[j - 2];
      sizes_.push_back(policy_.next(sizes_[j - 1], attempts, conflicts));
    }
    return sizes_[k];
  }
  // Batch k is decided: its attempts and conflicts.
  void record(const int k, const int64_t attempts, const int64_t conflicts)
  {
    if (static_cast<int>(outcomes_.size()) != k) {
      throw std::logic_error("TailBatchSizes: outcomes out of order");
    }
    outcomes_.emplace_back(attempts, conflicts);
  }

 private:
  TailBatchPolicy policy_;
  int max_;
  std::vector<int> sizes_;
  std::vector<std::pair<int64_t, int64_t>> outcomes_;
};

// Outcomes of one speculative tail pass (log statistics only). Every
// attempt is decided once: a complete route committed, a partial commit
// (keep_partial: its connections stay, the net stays unrouted), a conflict
// (searched again in a later batch) or failed. The unrouted nets, partial
// commits and failed attempts, go to the deferred pass (over the effort
// budget) or to the rip-up phase.
struct TailOutcomes
{
  int64_t attempts = 0;
  int64_t complete = 0;
  int64_t partial = 0;
  int64_t conflicts = 0;
  int64_t failed = 0;
  int64_t deferred = 0;
  int64_t deferred_partial = 0;
  int64_t to_ripup = 0;
  int64_t to_ripup_partial = 0;

  void decided(const TailPipeline::Decision decision, const bool partial_commit)
  {
    attempts++;
    switch (decision) {
      case TailPipeline::Decision::kCommitted:
        (partial_commit ? partial : complete)++;
        break;
      case TailPipeline::Decision::kConflict:
        conflicts++;
        break;
      case TailPipeline::Decision::kFailed:
        failed++;
        break;
    }
  }
  void unrouted(const bool over_budget, const bool partial_commit)
  {
    if (over_budget) {
      deferred++;
      deferred_partial += partial_commit;
    } else {
      to_ripup++;
      to_ripup_partial += partial_commit;
    }
  }
  // Failed attempts (no route) that go to the rip-up phase.
  int64_t failedToRipup() const { return to_ripup - to_ripup_partial; }
  // Every attempt has one outcome; every unrouted net (partial commit or
  // failed attempt) goes to exactly one of the deferred pass and the
  // rip-up phase.
  bool conserved() const
  {
    const int64_t counts[] = {attempts,
                              complete,
                              partial,
                              conflicts,
                              failed,
                              deferred,
                              deferred_partial,
                              to_ripup,
                              to_ripup_partial};
    for (const int64_t c : counts) {
      if (c < 0) {
        return false;
      }
    }
    return attempts == complete + partial + conflicts + failed
           && partial + failed == deferred + to_ripup
           && deferred_partial + to_ripup_partial == partial
           && deferred_partial <= deferred && to_ripup_partial <= to_ripup;
  }
};

class TailLookahead
{
 public:
  using Box = TailPipeline::Box;
  using Decision = TailPipeline::Decision;
  using Hooks = TailPipeline::Hooks;
  // overlapped_searches counts the searches started under R3.
  using Stats = TailPipeline::Stats;
  // One attempt of the semantics, in (batch, position) order.
  struct Attempt
  {
    int net;
    int batch;
    int position;
    Decision decision;
    bool early;  // searched under R3
  };

  explicit TailLookahead(int batch) : TailLookahead(TailBatchPolicy{batch}) {}
  explicit TailLookahead(const TailBatchPolicy& policy)
      : batch_(policy.max < 1 ? 1 : policy.max), policy_(policy)
  {
  }

  // boxes: indexed by net id. failed receives the nets without a route in
  // (batch, position) order. attempts (optional) receives every attempt in
  // (batch, position) order.
  Stats run(const std::vector<int>& order,
            const std::vector<Box>& boxes,
            int threads,
            const Hooks& hooks,
            std::vector<int>& failed,
            std::vector<Attempt>* attempts = nullptr);

 private:
  struct Item
  {
    int net;
    int batch = -1;  // -1: searched under R3 before its batch is formed
    int position = -1;
    bool early = false;
    bool search_started = false;
    bool searched = false;
    bool decided = false;
    bool found = false;
    Decision decision = Decision::kFailed;
    int decide_deps = -1;  // -1: not computed yet
    double cost = 0;
    std::vector<int> on_searched;  // items whose decide_deps count it
    std::vector<int> on_decided;   // items whose decide_deps count it
  };

  static bool meets(const Box& a, const Box& b)
  {
    return a[0] <= b[2] && b[0] <= a[2] && a[1] <= b[3] && b[1] <= a[3];
  }
  bool meets(int a, int b) const
  {
    return meets((*boxes_)[items_[a].net], (*boxes_)[items_[b].net]);
  }
  int addItem(int net);
  void formBatch();
  void computeDecideDeps(int item);
  void readySearch(int item);
  void readyDecide(int item);
  int lookaheadCandidate();
  void searched(int item);
  void decided(int item);
  bool finished() const { return current_.empty(); }

  const int batch_;  // largest batch size
  const TailBatchPolicy policy_;
  std::optional<TailBatchSizes> sizes_;
  const std::vector<int>* order_ = nullptr;
  const std::vector<Box>* boxes_ = nullptr;
  const Hooks* hooks_ = nullptr;
  std::vector<Item> items_;
  // Items of the formed batches, in (batch, position) order.
  std::vector<std::vector<int>> batches_;
  // Items of the batch being decided (batches_.back()).
  std::vector<int> current_;
  int current_decided_ = 0;
  int current_conflicts_ = 0;
  // Conflicts of the current batch: (position, net).
  std::vector<std::pair<int, int>> retries_;
  // First net of order_ that is in no formed batch.
  size_t queue_head_ = 0;
  // Items searched under R3, by p = index in order_ - queue_head_ (-1: none).
  std::vector<int> early_;
  int running_ = 0;
  bool deadlock_ = false;
  // Ready queues of the current batch: decisions by position; searches by
  // (-cost, position).
  std::priority_queue<std::pair<int, int>,
                      std::vector<std::pair<int, int>>,
                      std::greater<>>
      ready_decide_;
  std::priority_queue<std::tuple<double, int, int>,
                      std::vector<std::tuple<double, int, int>>,
                      std::greater<>>
      ready_search_;
  std::mutex mutex_;
  std::condition_variable wake_;
};

}  // namespace drt::br
