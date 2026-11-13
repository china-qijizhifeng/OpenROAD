// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>

namespace drt {

// Outcome of one maze search and of routing one net.  kCancelled means the
// worker observed its cancel request and the search ended without a result;
// it is a separate outcome from kNoPath.
enum class MazeSearchStatus
{
  kFound,
  kNoPath,
  kCancelled
};

// Cancel request of one run of a FlexDR worker (FlexDR::runGridWorkersAhead).
// The scheduler requests a cancel only for a run whose worker it has detached
// from its slot: such a result is never written back, so stopping the run
// changes no output.  The scheduler owns the request until the thread of the
// run has returned.
class WorkerCancelRequest
{
 public:
  using Clock = std::chrono::steady_clock;

  // Requester side.  Called with the scheduler's mutex held; the request time
  // and the pop count are read back under the same mutex after the run's
  // thread returned.
  void request()
  {
    requested_at_ = Clock::now();
    pops_at_request_ = published_pops_.load(std::memory_order_relaxed);
    requested_.store(true, std::memory_order_release);
  }
  bool isRequested() const
  {
    return requested_.load(std::memory_order_acquire);
  }
  Clock::time_point getRequestTime() const { return requested_at_; }
  // A* pops the worker had published when the cancel was requested; the
  // worker publishes at every check, so this lags the true count by less
  // than WorkerCancelPoll::kPopsPerCheck pops.
  uint64_t getPopsAtRequest() const { return pops_at_request_; }

  // Worker side.
  void publishPops(const uint64_t pops)
  {
    published_pops_.store(pops, std::memory_order_relaxed);
  }

 private:
  std::atomic<bool> requested_{false};
  std::atomic<uint64_t> published_pops_{0};
  Clock::time_point requested_at_;
  uint64_t pops_at_request_ = 0;
};

// Worker side of the cancel request: counts the A* pops of the worker and
// polls the request at the coarse check points and after every
// kPopsPerCheck pops.  Without a request object every poll returns false.
class WorkerCancelPoll
{
 public:
  static constexpr uint64_t kPopsPerCheck = 1024;  // a power of two
  static_assert((kPopsPerCheck & (kPopsPerCheck - 1)) == 0);

  void attach(WorkerCancelRequest* request) { request_ = request; }
  // Coarse check point.  True once the request is observed; it stays true.
  bool poll()
  {
    if (cancelled_) {
      return true;
    }
    if (request_ == nullptr) {
      return false;
    }
    request_->publishPops(pops_);
    cancelled_ = request_->isRequested();
    return cancelled_;
  }
  // After every A* pop.  True once the request is observed.
  bool afterPop()
  {
    ++pops_;
    if ((pops_ & (kPopsPerCheck - 1)) != 0) {
      return cancelled_;
    }
    return poll();
  }
  bool isCancelled() const { return cancelled_; }
  uint64_t getPops() const { return pops_; }

 private:
  WorkerCancelRequest* request_ = nullptr;
  uint64_t pops_ = 0;
  bool cancelled_ = false;
};

// Log-only counters of the cancelled runs of one clean-up iteration.
struct WorkerCancelStatistics
{
  using Clock = WorkerCancelRequest::Clock;

  int requests = 0;   // runs detached from their slot while running
  int cancelled = 0;  // of those: ended at a cancel check
  int completed = 0;  // of those: ended their routing before a check
  double exit_latency_sum_ms = 0;
  double exit_latency_max_ms = 0;
  uint64_t pops_after_request_sum = 0;
  uint64_t pops_after_request_max = 0;
  // Thread seconds of the detached runs, from their start to their exit.
  double detached_run_seconds = 0;
  double detached_run_seconds_after_request = 0;

  // A detached run whose thread returned at exit_time.
  void addExit(const WorkerCancelRequest& request,
               const bool observed_cancel,
               const uint64_t pops_at_exit,
               const Clock::time_point run_start,
               const Clock::time_point exit_time)
  {
    if (observed_cancel) {
      cancelled++;
    } else {
      completed++;
    }
    const double latency_ms = std::max(
        0.0,
        std::chrono::duration<double, std::milli>(exit_time
                                                  - request.getRequestTime())
            .count());
    exit_latency_sum_ms += latency_ms;
    exit_latency_max_ms = std::max(exit_latency_max_ms, latency_ms);
    const uint64_t pops_after
        = pops_at_exit > request.getPopsAtRequest()
              ? pops_at_exit - request.getPopsAtRequest()
              : 0;
    pops_after_request_sum += pops_after;
    pops_after_request_max = std::max(pops_after_request_max, pops_after);
    detached_run_seconds
        += std::chrono::duration<double>(exit_time - run_start).count();
    detached_run_seconds_after_request += latency_ms / 1000.0;
  }
  int exits() const { return cancelled + completed; }
};

}  // namespace drt
