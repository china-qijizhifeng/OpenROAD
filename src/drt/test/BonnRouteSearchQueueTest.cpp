// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// SearchQueue (aligned 4-ary heap plus the held last push) against a
// reference priority queue: random push / pop sequences shaped like the
// vertex search (monotone keys with many ties, bursts of pushes, the last
// push often the next minimum) must pop the same items in the same order.
// VertexOrder: pack / unpack round trips and the packed order equals the
// lexicographic order of (z, i, j, dir).

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <queue>
#include <random>
#include <tuple>
#include <vector>

#include "br/BrSearchQueue.h"

using namespace drt::br;

namespace {

struct TestItem
{
  int64_t key;
  int64_t g;
  uint64_t order;
  int32_t aux;
  bool operator>(const TestItem& o) const
  {
    if (key != o.key) {
      return key > o.key;
    }
    if (g != o.g) {
      return g < o.g;
    }
    return order > o.order;
  }
  bool operator==(const TestItem& o) const
  {
    return key == o.key && g == o.g && order == o.order && aux == o.aux;
  }
};

struct Greater
{
  bool operator()(const TestItem& a, const TestItem& b) const { return a > b; }
};

int fails = 0;
void fail(const char* what, const int c)
{
  if (fails++ < 12) {
    std::printf("case %d: %s\n", c, what);
  }
}

}  // namespace

int main()
{
  std::mt19937_64 rng(5);
  int64_t pops = 0;
  int64_t held = 0;
  int64_t max_size = 0;
  SearchQueue<TestItem> queue;  // reused across cases, as in PathSearch
  for (int c = 0; c < 3000; c++) {
    queue.clear();
    std::priority_queue<TestItem, std::vector<TestItem>, Greater> ref;
    uint64_t next_order = 0;
    int64_t floor = 0;  // keys at least the last popped key (monotone)
    const int key_range = 1 + static_cast<int>(rng() % (c % 3 == 0 ? 3 : 40));
    const int steps = 1 + static_cast<int>(rng() % (c % 10 == 0 ? 20000 : 600));
    const int64_t held_before = queue.heldPops();
    for (int s = 0; s < steps; s++) {
      const int burst = static_cast<int>(rng() % 5);
      for (int b = 0; b < burst; b++) {
        TestItem item;
        item.key = floor + static_cast<int64_t>(rng() % key_range);
        item.g = static_cast<int64_t>(rng() % 4);
        // Unique orders keep the order strict (as vertices are queued at
        // most once per g in the search).
        next_order += 1 + rng() % 3;
        item.order = rng() % 2 ? next_order : ~next_order;
        item.aux = static_cast<int32_t>(rng() % 1000);
        queue.push(item);
        ref.push(item);
      }
      max_size = std::max<int64_t>(max_size, static_cast<int64_t>(ref.size()));
      if (queue.size() != ref.size() || queue.empty() != ref.empty()) {
        fail("size differs", c);
      }
      const int npop = static_cast<int>(rng() % 3);
      for (int p = 0; p < npop && !ref.empty(); p++) {
        const TestItem a = queue.pop();
        const TestItem b = ref.top();
        ref.pop();
        pops++;
        if (!(a == b)) {
          fail("pop differs", c);
        }
        floor = b.key;
      }
    }
    while (!ref.empty()) {
      const TestItem a = queue.pop();
      const TestItem b = ref.top();
      ref.pop();
      pops++;
      if (!(a == b)) {
        fail("drain pop differs", c);
      }
    }
    if (!queue.empty()) {
      fail("queue not empty after the drain", c);
    }
    held += queue.heldPops() - held_before;
  }
  // VertexOrder.
  int64_t order_checks = 0;
  for (int c = 0; c < 2000; c++) {
    VertexOrder order;
    const int layers = 1 + static_cast<int>(rng() % 12);
    const int tracks = 1 + static_cast<int>(rng() % (c % 2 ? 300000 : 70));
    const int positions = 1 + static_cast<int>(rng() % (c % 3 ? 300000 : 90));
    if (!order.setup(layers, tracks, positions)) {
      fail("order setup rejected a small graph", c);
      continue;
    }
    auto random_vertex = [&]() {
      return std::make_tuple(static_cast<int>(rng() % layers),
                             static_cast<int>(rng() % tracks),
                             static_cast<int>(rng() % positions),
                             static_cast<int>(rng() % 4) - 1);
    };
    for (int k = 0; k < 200; k++) {
      const auto a = random_vertex();
      auto b = random_vertex();
      if (rng() % 3 == 0) {
        std::get<0>(b) = std::get<0>(a);
        if (rng() % 2 == 0) {
          std::get<1>(b) = std::get<1>(a);
        }
      }
      const uint64_t pa = order.pack(
          std::get<0>(a), std::get<1>(a), std::get<2>(a), std::get<3>(a));
      const uint64_t pb = order.pack(
          std::get<0>(b), std::get<1>(b), std::get<2>(b), std::get<3>(b));
      order_checks++;
      if (order.z(pa) != std::get<0>(a) || order.i(pa) != std::get<1>(a)
          || order.j(pa) != std::get<2>(a)
          || VertexOrder::dir(pa) != std::get<3>(a)) {
        fail("pack / unpack round trip", c);
      }
      if ((pa < pb) != (a < b) || (pa == pb) != (a == b)) {
        fail("packed order differs from (z, i, j, dir)", c);
      }
    }
  }
  VertexOrder big;
  if (big.setup(64, 1 << 30, 1 << 30)) {
    fail("order setup accepted 70 bits", -1);
  }
  std::printf(
      "search queue: %lld pops (%lld from the held item), largest queue "
      "%lld, %lld order checks, %d failures\n",
      static_cast<long long>(pops),
      static_cast<long long>(held),
      static_cast<long long>(max_size),
      static_cast<long long>(order_checks),
      fails);
  return fails != 0;
}
