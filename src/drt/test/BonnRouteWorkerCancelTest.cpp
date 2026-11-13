// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Cancel of detached FlexDR clean-up workers (dr/FlexDRWorkerCancel.h).
// The worker side polls the request at coarse check points and after every
// WorkerCancelPoll::kPopsPerCheck A* pops.  Checked here:
//  - without a request object, and with a request that is never set, no poll
//    reports a cancel, and the pops are counted;
//  - a coarse poll observes a set request at once; the A* poll observes it
//    at the next check, so at most kPopsPerCheck pops follow the request;
//    the observed state is sticky;
//  - a search loop written like FlexGridGraph::search reports kCancelled for
//    a cancelled search, never kNoPath, and kNoPath / kFound without a
//    cancel exactly as without the poll (same pops, same result);
//  - a request from another thread stops a busy worker thread; the counters
//    of WorkerCancelStatistics add up.
// Also run it with -fsanitize=thread.

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <queue>
#include <thread>
#include <utility>
#include <vector>

#include "dr/FlexDRWorkerCancel.h"

using drt::MazeSearchStatus;
using drt::WorkerCancelPoll;
using drt::WorkerCancelRequest;
using drt::WorkerCancelStatistics;

namespace {

int failures = 0;

void check(const bool condition, const char* what, const int line)
{
  if (!condition) {
    std::printf("FAIL line %d: %s\n", line, what);
    failures++;
  }
}
#define CHECK(condition) check((condition), #condition, __LINE__)

constexpr uint64_t kPopsPerCheck = WorkerCancelPoll::kPopsPerCheck;

void testNoRequest()
{
  WorkerCancelPoll detached_never;
  for (uint64_t i = 0; i < 10 * kPopsPerCheck + 7; i++) {
    CHECK(!detached_never.afterPop());
  }
  CHECK(!detached_never.poll());
  CHECK(!detached_never.isCancelled());
  CHECK(detached_never.getPops() == 10 * kPopsPerCheck + 7);

  WorkerCancelRequest request;
  WorkerCancelPoll poll;
  poll.attach(&request);
  for (uint64_t i = 0; i < 5 * kPopsPerCheck; i++) {
    CHECK(!poll.afterPop());
  }
  CHECK(!poll.poll());
  CHECK(!poll.isCancelled());
  CHECK(!request.isRequested());
}

void testCoarsePoll()
{
  WorkerCancelRequest request;
  WorkerCancelPoll poll;
  poll.attach(&request);
  CHECK(!poll.poll());
  request.request();
  CHECK(request.isRequested());
  CHECK(poll.poll());
  CHECK(poll.isCancelled());
  CHECK(poll.poll());  // sticky
  CHECK(poll.afterPop());
}

void testAstarPollBound()
{
  for (const uint64_t pops_before :
       {uint64_t{0}, uint64_t{1}, kPopsPerCheck - 1, kPopsPerCheck,
        kPopsPerCheck + 1, 7 * kPopsPerCheck + 513}) {
    WorkerCancelRequest request;
    WorkerCancelPoll poll;
    poll.attach(&request);
    for (uint64_t i = 0; i < pops_before; i++) {
      CHECK(!poll.afterPop());
    }
    request.request();
    // The published count lags the true count by less than one check.
    CHECK(request.getPopsAtRequest() <= pops_before);
    CHECK(pops_before - request.getPopsAtRequest() < kPopsPerCheck);
    uint64_t pops_after = 0;
    bool cancelled = false;
    while (!cancelled) {
      pops_after++;
      cancelled = poll.afterPop();
      CHECK(pops_after <= kPopsPerCheck);
      if (pops_after > kPopsPerCheck) {
        break;
      }
    }
    CHECK(cancelled);
    CHECK(poll.isCancelled());
    CHECK(poll.getPops() == pops_before + pops_after);
    CHECK(poll.getPops() % kPopsPerCheck == 0);
  }
}

// A* on a W x H 4-neighbour grid with blocked cells, shaped like
// FlexGridGraph::search: the poll runs right after every pop.
struct Grid
{
  int width;
  int height;
  std::vector<char> blocked;
};

MazeSearchStatus search(const Grid& grid,
                        const int source,
                        const int target,
                        WorkerCancelPoll& poll,
                        const bool with_poll,
                        uint64_t& pops)
{
  std::vector<int> cost(grid.width * grid.height, -1);
  std::priority_queue<std::pair<int, int>,
                      std::vector<std::pair<int, int>>,
                      std::greater<>>
      wavefront;
  wavefront.emplace(0, source);
  std::vector<char> done(grid.width * grid.height, 0);
  cost[source] = 0;
  while (!wavefront.empty()) {
    const auto [current_cost, node] = wavefront.top();
    wavefront.pop();
    pops++;
    if (with_poll && poll.afterPop()) {
      return MazeSearchStatus::kCancelled;
    }
    if (done[node]) {
      continue;
    }
    done[node] = 1;
    if (node == target) {
      return MazeSearchStatus::kFound;
    }
    const int x = node % grid.width;
    const int y = node / grid.width;
    const int neighbours[4][2] = {{x + 1, y}, {x - 1, y}, {x, y + 1}, {x, y - 1}};
    for (const auto& [nx, ny] : neighbours) {
      if (nx < 0 || ny < 0 || nx >= grid.width || ny >= grid.height) {
        continue;
      }
      const int next = ny * grid.width + nx;
      if (grid.blocked[next]
          || (cost[next] >= 0 && cost[next] <= current_cost + 1)) {
        continue;
      }
      cost[next] = current_cost + 1;
      wavefront.emplace(current_cost + 1, next);
    }
  }
  return MazeSearchStatus::kNoPath;
}

void testSearchStatus()
{
  // 100 x 100 grid; the target in the corner is walled in.
  Grid grid{100, 100, std::vector<char>(100 * 100, 0)};
  const int target = 99 * 100 + 99;
  grid.blocked[98 * 100 + 99] = 1;
  grid.blocked[99 * 100 + 98] = 1;
  const int reachable = 50 * 100 + 60;
  for (const int goal : {target, reachable}) {
    const MazeSearchStatus expected = goal == target
                                          ? MazeSearchStatus::kNoPath
                                          : MazeSearchStatus::kFound;
    // Without a request the poll changes neither the result nor the pops.
    uint64_t pops_plain = 0;
    WorkerCancelPoll unused;
    CHECK(search(grid, 0, goal, unused, false, pops_plain) == expected);
    WorkerCancelRequest request;
    WorkerCancelPoll poll;
    poll.attach(&request);
    uint64_t pops_polled = 0;
    CHECK(search(grid, 0, goal, poll, true, pops_polled) == expected);
    CHECK(pops_polled == pops_plain);
    CHECK(!poll.isCancelled());
  }
  // A cancel requested before or during an unreachable search ends it as
  // kCancelled, never as kNoPath.
  for (const uint64_t pops_before : {uint64_t{0}, uint64_t{3000}}) {
    WorkerCancelRequest request;
    WorkerCancelPoll poll;
    poll.attach(&request);
    for (uint64_t i = 0; i < pops_before; i++) {
      poll.afterPop();
    }
    request.request();
    uint64_t pops = 0;
    CHECK(search(grid, 0, target, poll, true, pops)
          == MazeSearchStatus::kCancelled);
    CHECK(pops <= kPopsPerCheck);
  }
}

void testCrossThread()
{
  WorkerCancelStatistics statistics;
  for (int round = 0; round < 20; round++) {
    WorkerCancelRequest request;
    WorkerCancelPoll poll;
    poll.attach(&request);
    std::atomic<bool> started{false};
    std::atomic<bool> exited_by_cancel{false};
    uint64_t pops_at_exit = 0;
    const auto run_start = WorkerCancelRequest::Clock::now();
    std::thread worker([&]() {
      started.store(true);
      // Busy until the cancel is observed; bounded for a broken poll.
      for (uint64_t i = 0; i < (uint64_t{1} << 36); i++) {
        if (poll.afterPop()) {
          exited_by_cancel.store(true);
          break;
        }
      }
      pops_at_exit = poll.getPops();
    });
    while (!started.load()) {
      std::this_thread::yield();
    }
    std::this_thread::sleep_for(std::chrono::microseconds(100 * (round % 5)));
    request.request();
    worker.join();
    const auto exit_time = WorkerCancelRequest::Clock::now();
    CHECK(exited_by_cancel.load());
    CHECK(pops_at_exit >= request.getPopsAtRequest());
    CHECK(pops_at_exit - request.getPopsAtRequest() <= 2 * kPopsPerCheck);
    statistics.requests++;
    statistics.addExit(request, poll.isCancelled(), pops_at_exit, run_start,
                       exit_time);
  }
  CHECK(statistics.requests == 20);
  CHECK(statistics.cancelled == 20);
  CHECK(statistics.completed == 0);
  CHECK(statistics.exits() == 20);
  CHECK(statistics.pops_after_request_max <= 2 * kPopsPerCheck);
  CHECK(statistics.exit_latency_max_ms >= 0);
  CHECK(statistics.exit_latency_sum_ms >= statistics.exit_latency_max_ms);
  CHECK(statistics.detached_run_seconds
        >= statistics.detached_run_seconds_after_request);
}

void testStatisticsCompleted()
{
  // A run that ended its routing before any check after the request.
  WorkerCancelStatistics statistics;
  WorkerCancelRequest request;
  WorkerCancelPoll poll;
  poll.attach(&request);
  for (int i = 0; i < 100; i++) {
    poll.afterPop();
  }
  const auto run_start = WorkerCancelRequest::Clock::now();
  request.request();
  statistics.requests++;
  statistics.addExit(request, poll.isCancelled(), poll.getPops(), run_start,
                     WorkerCancelRequest::Clock::now());
  CHECK(statistics.cancelled == 0);
  CHECK(statistics.completed == 1);
  CHECK(statistics.pops_after_request_sum == 100);  // published at 0
}

}  // namespace

int main()
{
  testNoRequest();
  testCoarsePoll();
  testAstarPollBound();
  testSearchStatus();
  testCrossThread();
  testStatisticsCompleted();
  if (failures != 0) {
    std::printf("bonnroute_worker_cancel_test: %d failures\n", failures);
    return EXIT_FAILURE;
  }
  std::printf("bonnroute_worker_cancel_test: all checks passed\n");
  return EXIT_SUCCESS;
}
