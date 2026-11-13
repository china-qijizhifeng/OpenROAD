// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Incremental routing (arXiv:2111.06169, Sec. 5.2, Fig. 9): the
// classification of input pieces into kept pin-to-pin paths, reservations
// and dropped pieces, against brute force on random small piece graphs.
// Brute force enumerates every simple path between two distinct pins whose
// inner nodes are legal pieces or pins and marks its pieces kept.

#include <cstdint>
#include <cstdio>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include "br/BrIncrementalCore.h"

using namespace drt::br;

namespace {

struct BruteForce
{
  const PieceGraph& g;
  std::vector<std::vector<int>> adj;
  std::vector<char> on_path;
  std::vector<int> path;
  std::vector<char> kept;

  explicit BruteForce(const PieceGraph& graph)
      : g(graph),
        adj(graph.num_nodes),
        on_path(graph.num_nodes, 0),
        kept(graph.num_nodes, 0)
  {
    for (const auto& [a, b] : g.edges) {
      if (a != b) {
        adj[a].push_back(b);
        adj[b].push_back(a);
      }
    }
  }

  bool allowed(const int u) const { return g.is_pin[u] || g.legal[u]; }

  void extend(const int start, const int u)
  {
    if (u != start && g.is_pin[u]) {
      for (const int x : path) {
        if (!g.is_pin[x]) {
          kept[x] = 1;
        }
      }
    }
    for (const int w : adj[u]) {
      if (on_path[w] || !allowed(w)) {
        continue;
      }
      on_path[w] = 1;
      path.push_back(w);
      extend(start, w);
      path.pop_back();
      on_path[w] = 0;
    }
  }

  std::vector<char> run()
  {
    for (int p = 0; p < g.num_nodes; p++) {
      if (!g.is_pin[p]) {
        continue;
      }
      on_path[p] = 1;
      path.assign(1, p);
      extend(p, p);
      on_path[p] = 0;
    }
    return kept;
  }
};

// Labels by brute force: nodes connected through kept pieces and pins.
std::vector<int> bruteComponents(const PieceGraph& g,
                                 const std::vector<char>& kept)
{
  std::vector<int> label(g.num_nodes, -1);
  auto member = [&](int u) { return g.is_pin[u] || kept[u]; };
  int next = 0;
  for (int s = 0; s < g.num_nodes; s++) {
    if (!member(s) || label[s] >= 0) {
      continue;
    }
    std::vector<int> stack{s};
    label[s] = next;
    while (!stack.empty()) {
      const int u = stack.back();
      stack.pop_back();
      for (const auto& [a, b] : g.edges) {
        const int w = a == u ? b : (b == u ? a : -1);
        if (w >= 0 && member(w) && label[w] < 0) {
          label[w] = next;
          stack.push_back(w);
        }
      }
    }
    next++;
  }
  return label;
}

bool samePartition(const std::vector<int>& a, const std::vector<int>& b)
{
  if (a.size() != b.size()) {
    return false;
  }
  std::set<std::pair<int, int>> pairs;
  std::set<int> la;
  std::set<int> lb;
  for (size_t i = 0; i < a.size(); i++) {
    if ((a[i] < 0) != (b[i] < 0)) {
      return false;
    }
    if (a[i] >= 0) {
      pairs.insert({a[i], b[i]});
      la.insert(a[i]);
      lb.insert(b[i]);
    }
  }
  return pairs.size() == la.size() && pairs.size() == lb.size();
}

}  // namespace

int main()
{
  std::mt19937 rng(20261112);
  int fails = 0;
  int cases = 0;
  int64_t kept_total = 0;
  int64_t reserved_total = 0;
  for (int c = 0; c < 20000; c++) {
    PieceGraph g;
    g.num_nodes = 2 + static_cast<int>(rng() % 9);
    g.is_pin.assign(g.num_nodes, 0);
    g.legal.assign(g.num_nodes, 0);
    g.is_short.assign(g.num_nodes, 0);
    for (int u = 0; u < g.num_nodes; u++) {
      g.is_pin[u] = rng() % 4 == 0;
      g.legal[u] = rng() % 5 != 0;
      g.is_short[u] = rng() % 3 == 0;
    }
    const int num_edges = static_cast<int>(rng() % (2 * g.num_nodes + 1));
    for (int e = 0; e < num_edges; e++) {
      g.edges.emplace_back(static_cast<int>(rng() % g.num_nodes),
                           static_cast<int>(rng() % g.num_nodes));
    }
    cases++;
    const std::vector<char> kept = keptPieces(g);
    const std::vector<char> expect = BruteForce(g).run();
    if (kept != expect) {
      std::printf("case %d: kept pieces differ from brute force\n", c);
      fails++;
    }
    if (!samePartition(keptComponents(g, kept), bruteComponents(g, kept))) {
      std::printf("case %d: kept components differ\n", c);
      fails++;
    }
    // Reservation choice from the definition, neighbour by neighbour.
    const std::vector<char> reserve = reservedPieces(g, kept);
    for (int v = 0; v < g.num_nodes; v++) {
      bool want = !g.is_pin[v] && g.legal[v] && !kept[v];
      if (want && g.is_short[v]) {
        int neighbours = 0;
        bool other = false;
        for (const auto& [a, b] : g.edges) {
          if (a == b || (a != v && b != v)) {
            continue;
          }
          const int w = a == v ? b : a;
          neighbours++;
          other |= g.is_pin[w] || g.legal[w];
        }
        want = neighbours == 0 || other;
      }
      if (static_cast<bool>(reserve[v]) != want) {
        std::printf("case %d node %d: reservation choice differs\n", c, v);
        fails++;
      }
      kept_total += kept[v];
      reserved_total += reserve[v];
    }
  }
  std::printf(
      "incremental classification: %d cases, %lld kept, %lld reserved, %d "
      "failures\n",
      cases,
      static_cast<long long>(kept_total),
      static_cast<long long>(reserved_total),
      fails);
  return fails != 0;
}
