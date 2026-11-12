// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrTailPipeline.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <mutex>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace drt::br {

namespace {

double now()
{
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

}  // namespace

void TailPipeline::addItem(const int net, const int batch, const bool fresh)
{
  while (static_cast<int>(batches_.size()) <= batch) {
    batches_.emplace_back();
  }
  Item item;
  item.net = net;
  item.batch = batch;
  item.position = static_cast<int>(batches_[batch].items.size());
  item.fresh = fresh;
  item.cost = hooks_->cost ? hooks_->cost(net) : 0.0;
  batches_[batch].items.push_back(static_cast<int>(items_.size()));
  items_.push_back(std::move(item));
  open_items_++;
}

void TailPipeline::computeSearchDeps(const int x)
{
  // The commits of interacting attempts of earlier batches.
  int deps = 0;
  const int b = items_[x].batch;
  for (int bb = low_; bb < b; bb++) {
    for (const int y : batches_[bb].items) {
      if (!items_[y].decided && meets(x, y)) {
        deps++;
        items_[y].on_decided.push_back(x);
        items_[y].on_decided_kind.push_back(0);
      }
    }
  }
  items_[x].search_deps = deps;
  if (deps == 0) {
    readySearch(x);
  }
}

void TailPipeline::computeDecideDeps(const int x)
{
  // Its own search, the searches of interacting attempts of its batch (they
  // must not see its commit) and the decisions of interacting attempts
  // earlier in the batch (its decision sees their commits).
  Item& item = items_[x];
  int deps = item.searched ? 0 : 1;
  for (const int y : batches_[item.batch].items) {
    if (y == x || !meets(x, y)) {
      continue;
    }
    if (!items_[y].searched) {
      deps++;
      items_[y].on_searched.push_back(x);
    }
    if (items_[y].position < item.position && !items_[y].decided) {
      deps++;
      items_[y].on_decided.push_back(x);
      items_[y].on_decided_kind.push_back(1);
    }
  }
  items_[x].decide_deps = deps;
  if (deps == 0) {
    readyDecide(x);
  }
}

void TailPipeline::finalizeBatches()
{
  while (final_ < static_cast<int>(batches_.size())) {
    const int b = final_;
    if (b >= retry_delay_) {
      const Batch& source = batches_[b - retry_delay_];
      if (source.decided < static_cast<int>(source.items.size())) {
        return;
      }
    }
    final_++;
    std::vector<std::pair<int, int>> retries;
    retries.swap(batches_[b].retries);
    std::ranges::sort(retries);
    for (const auto& [position, net] : retries) {
      addItem(net, b, false);
    }
    for (const int x : batches_[b].items) {
      if (!items_[x].fresh) {
        computeSearchDeps(x);
      }
    }
    for (const int x : batches_[b].items) {
      computeDecideDeps(x);
    }
    if (b + 1 < static_cast<int>(batches_.size())) {
      for (const int x : batches_[b + 1].items) {
        if (items_[x].fresh) {
          computeSearchDeps(x);
        }
      }
    }
  }
}

void TailPipeline::readySearch(const int x)
{
  const Item& item = items_[x];
  ready_search_.emplace(item.batch, -item.cost, item.position, x);
}

void TailPipeline::readyDecide(const int x)
{
  const Item& item = items_[x];
  ready_decide_.emplace(
      (static_cast<int64_t>(item.batch) << 32) | item.position, x);
}

void TailPipeline::searched(const int x)
{
  Item& item = items_[x];
  item.searched = true;
  if (item.decide_deps > 0 && --item.decide_deps == 0) {
    readyDecide(x);
  }
  for (const int y : item.on_searched) {
    if (--items_[y].decide_deps == 0) {
      readyDecide(y);
    }
  }
}

void TailPipeline::decided(const int x)
{
  items_[x].decided = true;
  // Copies: addItem below may reallocate items_.
  const std::vector<int> waiters = items_[x].on_decided;
  const std::vector<char> kinds = items_[x].on_decided_kind;
  for (size_t k = 0; k < waiters.size(); k++) {
    Item& y = items_[waiters[k]];
    if (kinds[k] == 0) {
      if (--y.search_deps == 0) {
        readySearch(waiters[k]);
      }
    } else if (--y.decide_deps == 0) {
      readyDecide(waiters[k]);
    }
  }
  const int b = items_[x].batch;
  if (items_[x].decision == Decision::kConflict) {
    // Appended when the target batch becomes final, in the position order
    // of their source batch (decisions complete in any order).
    const int target = b + retry_delay_;
    while (static_cast<int>(batches_.size()) <= target) {
      batches_.emplace_back();
    }
    batches_[target].retries.emplace_back(items_[x].position, items_[x].net);
  }
  batches_[b].decided++;
  open_items_--;
  while (low_ < final_
         && batches_[low_].decided
                == static_cast<int>(batches_[low_].items.size())) {
    low_++;
  }
  finalizeBatches();
  while (low_ < final_
         && batches_[low_].decided
                == static_cast<int>(batches_[low_].items.size())) {
    low_++;
  }
}

bool TailPipeline::finished() const
{
  return open_items_ == 0 && final_ == static_cast<int>(batches_.size());
}

TailPipeline::Stats TailPipeline::run(const std::vector<int>& order,
                                      const std::vector<Box>& boxes,
                                      const int threads,
                                      const Hooks& hooks,
                                      std::vector<int>& failed)
{
  boxes_ = &boxes;
  hooks_ = &hooks;
  items_.clear();
  batches_.clear();
  final_ = 0;
  low_ = 0;
  open_items_ = 0;
  running_ = 0;
  deadlock_ = false;
  ready_decide_ = {};
  ready_search_ = {};
  Stats stats;
  items_.reserve(order.size() * 2);
  for (size_t i = 0; i < order.size(); i++) {
    addItem(order[i], static_cast<int>(i / batch_), true);
  }
  if (!batches_.empty()) {
    for (const int x : batches_[0].items) {
      computeSearchDeps(x);
    }
  }
  finalizeBatches();
  const int count = threads < 1 ? 1 : threads;
  std::vector<Stats> local(count);
  const double start = now();
  hooks.parallel(count, [&](const int thread) {
    Stats& mine = local[thread];
    std::vector<int> earlier;
    std::unique_lock<std::mutex> lock(mutex_);
    while (!finished() && !deadlock_) {
      int x = -1;
      bool decide = false;
      if (!ready_decide_.empty()) {
        x = ready_decide_.top().second;
        ready_decide_.pop();
        decide = true;
      } else if (!ready_search_.empty()) {
        x = std::get<3>(ready_search_.top());
        ready_search_.pop();
      } else {
        if (running_ == 0) {
          // Nothing runs and nothing is ready: a dependency cycle (bug).
          deadlock_ = true;
          wake_.notify_all();
          break;
        }
        const double wait_start = now();
        wake_.wait(lock);
        mine.idle += now() - wait_start;
        continue;
      }
      const int net = items_[x].net;
      const bool found = items_[x].found;
      earlier.clear();
      if (decide) {
        items_[x].decide_started = true;
        // Every interacting attempt earlier in the batch is decided.
        for (const int y : batches_[items_[x].batch].items) {
          if (items_[y].position < items_[x].position && meets(x, y)
              && items_[y].decision == Decision::kCommitted) {
            earlier.push_back(items_[y].net);
          }
        }
      } else {
        items_[x].search_started = true;
        if (items_[x].batch > low_) {
          mine.overlapped_searches++;
        }
      }
      running_++;
      lock.unlock();
      const double task_start = now();
      bool result_found = false;
      Decision decision = Decision::kFailed;
      if (decide) {
        decision = hooks.decide(net, found, earlier, thread);
      } else {
        result_found = hooks.search(net, thread);
      }
      const double task_seconds = now() - task_start;
      lock.lock();
      running_--;
      if (decide) {
        mine.decide_busy += task_seconds;
        items_[x].decision = decision;
        decided(x);
      } else {
        mine.search_busy += task_seconds;
        items_[x].found = result_found;
        searched(x);
      }
      wake_.notify_all();
    }
    wake_.notify_all();
  });
  stats.wall = now() - start;
  if (deadlock_) {
    throw std::logic_error("TailPipeline: no ready search or decision");
  }
  for (const Stats& s : local) {
    stats.search_busy += s.search_busy;
    stats.decide_busy += s.decide_busy;
    stats.idle += s.idle;
    stats.overlapped_searches += s.overlapped_searches;
  }
  stats.batches = static_cast<int64_t>(batches_.size());
  stats.attempts = static_cast<int64_t>(items_.size());
  for (const Batch& batch : batches_) {
    for (const int x : batch.items) {
      switch (items_[x].decision) {
        case Decision::kCommitted:
          stats.committed++;
          break;
        case Decision::kConflict:
          stats.conflicts++;
          break;
        case Decision::kFailed:
          stats.failed++;
          failed.push_back(items_[x].net);
          break;
      }
    }
  }
  return stats;
}


int TailLookahead::addItem(const int net)
{
  Item item;
  item.net = net;
  item.cost = hooks_->cost ? hooks_->cost(net) : 0.0;
  items_.push_back(std::move(item));
  return static_cast<int>(items_.size()) - 1;
}

void TailLookahead::formBatch()
{
  // R1: the conflicts of the previous batch in list order, then fresh nets
  // of the queue; R4: a fresh net searched under R3 keeps that search.
  std::ranges::sort(retries_);
  std::vector<int> members;
  members.reserve(batch_);
  for (const auto& [position, net] : retries_) {
    members.push_back(addItem(net));
  }
  retries_.clear();
  const size_t room = static_cast<size_t>(batch_) - members.size();
  const size_t fresh = std::min(room, order_->size() - queue_head_);
  for (size_t p = 0; p < early_.size(); p++) {
    if (p < fresh) {
      members.push_back(early_[p] >= 0
                            ? early_[p]
                            : addItem((*order_)[queue_head_ + p]));
    } else if (early_[p] >= 0) {
      throw std::logic_error(
          "TailLookahead: a search under R3 outside the next batch");
    }
  }
  queue_head_ += fresh;
  early_.assign(batch_, -1);
  current_ = std::move(members);
  current_decided_ = 0;
  current_conflicts_ = 0;
  if (current_.empty()) {
    return;
  }
  const int b = static_cast<int>(batches_.size());
  for (size_t k = 0; k < current_.size(); k++) {
    items_[current_[k]].batch = b;
    items_[current_[k]].position = static_cast<int>(k);
  }
  batches_.push_back(current_);
  for (const int x : current_) {
    if (!items_[x].search_started) {
      readySearch(x);
    }
  }
  for (const int x : current_) {
    computeDecideDeps(x);
  }
}

void TailLookahead::computeDecideDeps(const int x)
{
  // R2: its own search, the searches of interacting attempts of its batch
  // (they must not see its commit) and the decisions of interacting
  // attempts earlier in the batch (its decision sees their commits).
  int deps = items_[x].searched ? 0 : 1;
  for (const int y : current_) {
    if (y == x || !meets(x, y)) {
      continue;
    }
    if (!items_[y].searched) {
      deps++;
      items_[y].on_searched.push_back(x);
    }
    if (items_[y].position < items_[x].position && !items_[y].decided) {
      deps++;
      items_[y].on_decided.push_back(x);
    }
  }
  items_[x].decide_deps = deps;
  if (deps == 0) {
    readyDecide(x);
  }
}

void TailLookahead::readySearch(const int x)
{
  ready_search_.emplace(-items_[x].cost, items_[x].position, x);
}

void TailLookahead::readyDecide(const int x)
{
  ready_decide_.emplace(items_[x].position, x);
}

int TailLookahead::lookaheadCandidate()
{
  if (queue_head_ >= order_->size()) {
    return -1;
  }
  // U: attempts of the current batch that may still write the routing
  // space (search running or pending, or a route found and undecided).
  std::vector<int> writers;
  for (const int x : current_) {
    const Item& item = items_[x];
    if (!item.decided && (!item.searched || item.found)) {
      writers.push_back(x);
    }
  }
  // R3 (a): p < batch - r - |U|.
  const int64_t limit = std::min<int64_t>(
      static_cast<int64_t>(batch_) - current_conflicts_
          - static_cast<int64_t>(writers.size()),
      static_cast<int64_t>(order_->size() - queue_head_));
  for (int64_t p = 0; p < limit; p++) {
    if (early_[p] >= 0) {
      continue;
    }
    const int net = (*order_)[queue_head_ + p];
    // R3 (b): no box of U meets its box.
    const Box& box = (*boxes_)[net];
    const bool blocked = std::ranges::any_of(writers, [&](const int y) {
      return meets((*boxes_)[items_[y].net], box);
    });
    if (blocked) {
      continue;
    }
    const int x = addItem(net);
    items_[x].early = true;
    items_[x].search_started = true;
    early_[p] = x;
    return x;
  }
  return -1;
}

void TailLookahead::searched(const int x)
{
  items_[x].searched = true;
  if (items_[x].decide_deps > 0 && --items_[x].decide_deps == 0) {
    readyDecide(x);
  }
  for (const int y : items_[x].on_searched) {
    if (--items_[y].decide_deps == 0) {
      readyDecide(y);
    }
  }
}

void TailLookahead::decided(const int x)
{
  items_[x].decided = true;
  for (const int y : items_[x].on_decided) {
    if (--items_[y].decide_deps == 0) {
      readyDecide(y);
    }
  }
  if (items_[x].decision == Decision::kConflict) {
    retries_.emplace_back(items_[x].position, items_[x].net);
    current_conflicts_++;
  }
  current_decided_++;
  if (current_decided_ == static_cast<int>(current_.size())) {
    formBatch();
  }
}

TailLookahead::Stats TailLookahead::run(const std::vector<int>& order,
                                        const std::vector<Box>& boxes,
                                        const int threads,
                                        const Hooks& hooks,
                                        std::vector<int>& failed,
                                        std::vector<Attempt>* attempts)
{
  order_ = &order;
  boxes_ = &boxes;
  hooks_ = &hooks;
  items_.clear();
  batches_.clear();
  current_.clear();
  retries_.clear();
  queue_head_ = 0;
  early_.assign(batch_, -1);
  running_ = 0;
  deadlock_ = false;
  ready_decide_ = {};
  ready_search_ = {};
  items_.reserve(order.size() + order.size() / 4);
  formBatch();
  Stats stats;
  const int count = threads < 1 ? 1 : threads;
  std::vector<Stats> local(count);
  const double start = now();
  hooks.parallel(count, [&](const int thread) {
    Stats& mine = local[thread];
    std::vector<int> earlier;
    std::unique_lock<std::mutex> lock(mutex_);
    while (!finished() && !deadlock_) {
      int x = -1;
      bool decide = false;
      if (!ready_decide_.empty()) {
        x = ready_decide_.top().second;
        ready_decide_.pop();
        decide = true;
      } else if (!ready_search_.empty()) {
        x = std::get<2>(ready_search_.top());
        ready_search_.pop();
        items_[x].search_started = true;
      } else if ((x = lookaheadCandidate()) >= 0) {
        mine.overlapped_searches++;
      } else {
        if (running_ == 0) {
          // Nothing runs and nothing is ready: a dependency cycle (bug).
          deadlock_ = true;
          wake_.notify_all();
          break;
        }
        const double wait_start = now();
        wake_.wait(lock);
        mine.idle += now() - wait_start;
        continue;
      }
      const int net = items_[x].net;
      const bool found = items_[x].found;
      earlier.clear();
      if (decide) {
        // Every interacting attempt earlier in the batch is decided (R2).
        for (const int y : current_) {
          if (items_[y].position < items_[x].position && meets(x, y)
              && items_[y].decision == Decision::kCommitted) {
            earlier.push_back(items_[y].net);
          }
        }
      }
      running_++;
      lock.unlock();
      const double task_start = now();
      bool result_found = false;
      Decision decision = Decision::kFailed;
      if (decide) {
        decision = hooks.decide(net, found, earlier, thread);
      } else {
        result_found = hooks.search(net, thread);
      }
      const double task_seconds = now() - task_start;
      lock.lock();
      running_--;
      if (decide) {
        mine.decide_busy += task_seconds;
        items_[x].decision = decision;
        decided(x);
      } else {
        mine.search_busy += task_seconds;
        items_[x].found = result_found;
        searched(x);
      }
      wake_.notify_all();
    }
    wake_.notify_all();
  });
  stats.wall = now() - start;
  if (deadlock_) {
    throw std::logic_error("TailLookahead: no ready search or decision");
  }
  for (const Stats& s : local) {
    stats.search_busy += s.search_busy;
    stats.decide_busy += s.decide_busy;
    stats.idle += s.idle;
    stats.overlapped_searches += s.overlapped_searches;
  }
  stats.batches = static_cast<int64_t>(batches_.size());
  for (const std::vector<int>& batch : batches_) {
    for (const int x : batch) {
      const Item& item = items_[x];
      stats.attempts++;
      switch (item.decision) {
        case Decision::kCommitted:
          stats.committed++;
          break;
        case Decision::kConflict:
          stats.conflicts++;
          break;
        case Decision::kFailed:
          stats.failed++;
          failed.push_back(item.net);
          break;
      }
      if (attempts != nullptr) {
        attempts->push_back(
            {item.net, item.batch, item.position, item.decision, item.early});
      }
    }
  }
  return stats;
}

}  // namespace drt::br
