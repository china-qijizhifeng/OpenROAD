// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

// Ordered schedule of tasks with boxes (BonnRoute DR, Sec. 5.1 extension).
//
// Task i has a box (inclusive global routing tiles) that holds every read
// and write of its work on the routing space. Semantics: the tasks run one
// after the other in index order. Schedule: task j starts once every task
// i < j whose box shares a tile with box j has returned; tasks with
// disjoint boxes run concurrently. Each task thereby reads exactly the
// state the serial order defines, so the results equal the serial
// execution for every thread count and timing. Among the ready tasks the
// smallest index starts first.

#include <algorithm>
#include <array>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <mutex>
#include <queue>
#include <vector>

namespace drt::br {

class OrderedBoxSchedule
{
 public:
  // xlo, ylo, xhi, yhi in global routing tiles, inclusive.
  using Box = std::array<int, 4>;
  // Runs body(thread) for thread = 0 .. threads - 1 concurrently and
  // returns when all of them returned.
  using Parallel
      = std::function<void(int threads, const std::function<void(int)>& body)>;
  using Task = std::function<void(int index, int thread)>;
  struct Stats
  {
    int64_t tasks = 0;
    // Tasks on the longest chain of tasks that wait for each other.
    int64_t longest_chain = 0;
    int threads = 0;
    double wall = 0;
    // Thread seconds inside tasks, and waiting for a ready task.
    double busy = 0;
    double idle = 0;
  };

  static bool meets(const Box& a, const Box& b)
  {
    return a[0] <= b[2] && b[0] <= a[2] && a[1] <= b[3] && b[1] <= a[3];
  }

  // Runs task(i, thread) for every box i; at most `threads` at a time.
  static Stats run(const std::vector<Box>& boxes,
                   const int threads,
                   const Parallel& parallel,
                   const Task& task)
  {
    Stats stats;
    const int count = static_cast<int>(boxes.size());
    stats.tasks = count;
    if (count == 0) {
      return stats;
    }
    // Earlier tasks not returned yet whose boxes meet the task's box, and
    // the later tasks that wait for a task.
    std::vector<int> waiting(count, 0);
    std::vector<std::vector<int>> later(count);
    std::vector<int> chain(count, 1);
    for (int j = 0; j < count; j++) {
      for (int i = 0; i < j; i++) {
        if (meets(boxes[i], boxes[j])) {
          waiting[j]++;
          later[i].push_back(j);
          chain[j] = std::max(chain[j], chain[i] + 1);
        }
      }
      stats.longest_chain = std::max<int64_t>(stats.longest_chain, chain[j]);
    }
    std::priority_queue<int, std::vector<int>, std::greater<>> ready;
    for (int j = 0; j < count; j++) {
      if (waiting[j] == 0) {
        ready.push(j);
      }
    }
    int open = count;
    std::mutex mutex;
    std::condition_variable wake;
    const int used = std::clamp(threads, 1, count);
    stats.threads = used;
    std::vector<double> busy(used, 0);
    std::vector<double> idle(used, 0);
    const double start = now();
    parallel(used, [&](const int thread) {
      std::unique_lock<std::mutex> lock(mutex);
      while (open > 0) {
        if (ready.empty()) {
          const double wait_start = now();
          wake.wait(lock);
          idle[thread] += now() - wait_start;
          continue;
        }
        const int x = ready.top();
        ready.pop();
        lock.unlock();
        const double task_start = now();
        task(x, thread);
        const double task_seconds = now() - task_start;
        lock.lock();
        busy[thread] += task_seconds;
        open--;
        bool released = false;
        for (const int y : later[x]) {
          if (--waiting[y] == 0) {
            ready.push(y);
            released = true;
          }
        }
        if (released || open == 0) {
          wake.notify_all();
        }
      }
    });
    stats.wall = now() - start;
    for (int t = 0; t < used; t++) {
      stats.busy += busy[t];
      stats.idle += idle[t];
    }
    return stats;
  }

 private:
  static double now()
  {
    return std::chrono::duration<double>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
  }
};

}  // namespace drt::br
