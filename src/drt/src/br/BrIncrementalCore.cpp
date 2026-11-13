// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrIncrementalCore.h"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <utility>
#include <vector>

namespace drt::br {

namespace {

// Unit-capacity flow network with vertex splitting: node u becomes the arc
// in(u) = 2u -> out(u) = 2u + 1; sink = 2 * num_nodes.
class SplitNetwork
{
 public:
  explicit SplitNetwork(const PieceGraph& g) : n_(2 * g.num_nodes + 1)
  {
    head_.assign(n_, -1);
    for (int u = 0; u < g.num_nodes; u++) {
      split_.push_back(addArc(2 * u, 2 * u + 1, allowed(g, u) ? 1 : 0));
    }
    for (const auto& [a, b] : g.edges) {
      if (a == b || !allowed(g, a) || !allowed(g, b)) {
        continue;
      }
      addArc(2 * a + 1, 2 * b, 1);
      addArc(2 * b + 1, 2 * a, 1);
    }
    for (int u = 0; u < g.num_nodes; u++) {
      if (g.is_pin[u]) {
        addArc(2 * u + 1, sink(), 1);
      }
    }
    base_cap_ = cap_;
  }

  int sink() const { return n_ - 1; }

  // Number of vertex-disjoint paths from node v to distinct pins, up to
  // `limit` (v counts with capacity `limit`).
  int paths(const int v, const int limit)
  {
    cap_ = base_cap_;
    cap_[split_[v]] = limit;
    int flow = 0;
    std::vector<int> parent_arc(n_);
    std::vector<int> queue;
    while (flow < limit) {
      std::fill(parent_arc.begin(), parent_arc.end(), -2);
      queue.assign(1, 2 * v);
      parent_arc[2 * v] = -1;
      for (size_t k = 0; k < queue.size() && parent_arc[sink()] == -2; k++) {
        const int x = queue[k];
        for (int a = head_[x]; a >= 0; a = next_[a]) {
          const int y = to_[a];
          if (cap_[a] > 0 && parent_arc[y] == -2) {
            parent_arc[y] = a;
            queue.push_back(y);
          }
        }
      }
      if (parent_arc[sink()] == -2) {
        break;
      }
      for (int x = sink(); parent_arc[x] >= 0; x = to_[parent_arc[x] ^ 1]) {
        cap_[parent_arc[x]]--;
        cap_[parent_arc[x] ^ 1]++;
      }
      flow++;
    }
    return flow;
  }

 private:
  static bool allowed(const PieceGraph& g, const int u)
  {
    return g.is_pin[u] || g.legal[u];
  }

  int addArc(const int from, const int to, const int cap)
  {
    const int a = static_cast<int>(to_.size());
    to_.push_back(to);
    cap_.push_back(cap);
    next_.push_back(head_[from]);
    head_[from] = a;
    to_.push_back(from);
    cap_.push_back(0);
    next_.push_back(head_[to]);
    head_[to] = a + 1;
    return a;
  }

  int n_;
  std::vector<int> head_;
  std::vector<int> to_;
  std::vector<int> cap_;
  std::vector<int> base_cap_;
  std::vector<int> next_;
  std::vector<int> split_;
};

int find(std::vector<int>& parent, int x)
{
  while (parent[x] != x) {
    parent[x] = parent[parent[x]];
    x = parent[x];
  }
  return x;
}

}  // namespace

std::vector<char> keptPieces(const PieceGraph& g)
{
  std::vector<char> kept(g.num_nodes, 0);
  if (g.num_nodes == 0) {
    return kept;
  }
  SplitNetwork network(g);
  for (int v = 0; v < g.num_nodes; v++) {
    if (!g.is_pin[v] && g.legal[v]) {
      kept[v] = network.paths(v, 2) >= 2;
    }
  }
  return kept;
}

std::vector<int> keptComponents(const PieceGraph& g,
                                const std::vector<char>& kept)
{
  std::vector<int> parent(g.num_nodes);
  std::iota(parent.begin(), parent.end(), 0);
  auto member = [&](const int u) { return g.is_pin[u] || kept[u]; };
  for (const auto& [a, b] : g.edges) {
    if (member(a) && member(b)) {
      parent[find(parent, a)] = find(parent, b);
    }
  }
  std::vector<int> label(g.num_nodes, -1);
  std::vector<int> root_label(g.num_nodes, -1);
  int next = 0;
  for (int u = 0; u < g.num_nodes; u++) {
    if (!member(u)) {
      continue;
    }
    const int r = find(parent, u);
    if (root_label[r] < 0) {
      root_label[r] = next++;
    }
    label[u] = root_label[r];
  }
  return label;
}

std::vector<char> reservedPieces(const PieceGraph& g,
                                 const std::vector<char>& kept)
{
  std::vector<std::vector<int>> adj(g.num_nodes);
  for (const auto& [a, b] : g.edges) {
    if (a != b) {
      adj[a].push_back(b);
      adj[b].push_back(a);
    }
  }
  std::vector<char> reserve(g.num_nodes, 0);
  for (int v = 0; v < g.num_nodes; v++) {
    if (g.is_pin[v] || !g.legal[v] || kept[v]) {
      continue;
    }
    bool only_illegal = !adj[v].empty();
    for (const int w : adj[v]) {
      if (g.is_pin[w] || g.legal[w]) {
        only_illegal = false;
        break;
      }
    }
    reserve[v] = !(g.is_short[v] && only_illegal);
  }
  return reserve;
}

}  // namespace drt::br
