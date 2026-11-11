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
#include <cstdint>
#include <functional>
#include <mutex>
#include <queue>
#include <tuple>
#include <utility>
#include <vector>

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

}  // namespace drt::br
