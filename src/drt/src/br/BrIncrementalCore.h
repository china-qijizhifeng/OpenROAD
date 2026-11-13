// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Classification of a net's input routing for incremental routing
// (Ahrens, Henke, Rabenstein, Vygen, arXiv:2111.06169, Sec. 5.2, Fig. 9).
// Independent of the drt database: the adapter (BrIncremental.cpp) builds
// the piece graph of a net from its detailed wires, vias and pins.

#pragma once

#include <utility>
#include <vector>

namespace drt::br {

// Pieces of one net's input routing (wires, vias, patches) and its pins as
// nodes; an edge joins two nodes whose metal touches.
struct PieceGraph
{
  int num_nodes = 0;
  std::vector<std::pair<int, int>> edges;
  // Node is a pin of the net (pins are never kept, reserved or dropped).
  std::vector<char> is_pin;
  // Piece without a conflict with other nets' shapes (ignored for pins).
  std::vector<char> legal;
  // Piece shorter than the short-wire threshold (ignored for pins).
  std::vector<char> is_short;
};

// Sec. 5.2: "we keep pin-to-pin paths fixed if they consist only of
// detailed wires" (and only of legal ones). kept[v] holds for a legal piece
// v that lies on a simple path between two distinct pins whose inner nodes
// are legal pieces or pins. Exact: two vertex-disjoint paths from v to two
// distinct pins (Menger), found by unit-capacity augmenting paths.
std::vector<char> keptPieces(const PieceGraph& g);

// Components of the kept pieces together with the pins they touch: label
// per node (pins included), -1 for nodes outside every kept component. A
// pin without kept pieces gets its own label.
std::vector<int> keptComponents(const PieceGraph& g,
                                const std::vector<char>& kept);

// Sec. 5.2: legal input wires become reservations of their net, "except for
// short wires that connect only to an illegal input wire". Kept pieces stay
// fixed and are not reservations. reserve[v] holds for a legal piece v that
// is not kept and is not a short piece whose neighbours are all illegal
// pieces (at least one neighbour, no pin).
std::vector<char> reservedPieces(const PieceGraph& g,
                                 const std::vector<char>& kept);

}  // namespace drt::br
