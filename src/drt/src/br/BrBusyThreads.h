// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

// Busy-thread time line of a scheduler (BonnRoute DR statistics, log only):
// the wall seconds during which 0, 1, 2, 3-4, 5-8, 9-16, 17-32 and 33 or
// more threads were inside a task. The owner calls change() with its mutex
// held, with the number of busy threads before the change.

#include <array>
#include <cstddef>
#include <cstdio>
#include <string>

namespace drt::br {

class BusyThreads
{
 public:
  static constexpr int kClasses = 8;
  using Seconds = std::array<double, kClasses>;

  static int classOf(const int busy)
  {
    if (busy <= 2) {
      return busy < 0 ? 0 : busy;
    }
    int c = 3;
    for (int high = 4; busy > high && c < kClasses - 1; high *= 2) {
      c++;
    }
    return c;
  }

  void start(const double now) { last_ = now; }
  // `busy` threads were inside a task from the previous call until `now`.
  void change(const int busy, const double now)
  {
    seconds_[classOf(busy)] += now - last_;
    last_ = now;
  }
  const Seconds& seconds() const { return seconds_; }

  static void add(Seconds& total, const Seconds& part)
  {
    for (size_t c = 0; c < total.size(); c++) {
      total[c] += part[c];
    }
  }

  // "0: 0.10s, 1: 2.00s, 2: 0.50s, 3-4: 0.20s, ..., 33+: 0.00s"
  static std::string text(const Seconds& seconds)
  {
    static constexpr std::array<const char*, kClasses> kNames
        = {"0", "1", "2", "3-4", "5-8", "9-16", "17-32", "33+"};
    std::string out;
    for (int c = 0; c < kClasses; c++) {
      char buffer[48];
      std::snprintf(buffer,
                    sizeof(buffer),
                    "%s%s: %.2fs",
                    c == 0 ? "" : ", ",
                    kNames[c],
                    seconds[c]);
      out += buffer;
    }
    return out;
  }

 private:
  Seconds seconds_{};
  double last_ = 0;
};

}  // namespace drt::br
