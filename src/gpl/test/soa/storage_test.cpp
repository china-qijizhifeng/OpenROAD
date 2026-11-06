// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include <memory>
#include <stdexcept>
#include <vector>

#include "soaStorage.h"

namespace {
void require(bool condition)
{
  if (!condition) {
    throw std::runtime_error("resident SoA lifecycle invariant failed");
  }
}
}  // namespace

int main()
{
  using Store = gpl::SoaStorage<int, float, int64_t>;
  using Row = gpl::SoaRow<Store>;
  auto storage = std::make_shared<Store>();
  std::vector<Row> cells;
  for (int i = 0; i < 10000; ++i) {
    cells.emplace_back(storage);
    cells.back().get<0>() = i;
    cells.back().get<1>() = i * 0.25f;
    cells.back().get<2>() = int64_t(i) * i;
  }
  for (int i = 0; i < 10000; ++i) {
    require(cells[i].get<0>() == i);
    require(cells[i].get<1>() == i * 0.25f);
    require(std::as_const(cells[i]).get<2>() == int64_t(i) * i);
  }
  // A filler-removal snapshot retains values, not an alias of the live slot.
  Row snapshot(cells[3]);
  cells[3].get<0>() = -1;
  require(snapshot.get<0>() == 3);
  // Resize/insert/delete callbacks use swap-remove and may reallocate columns.
  cells[3] = std::move(cells.back());
  cells.pop_back();
  require(cells[3].get<0>() == 9999);
  for (int i = 0; i < 5000; ++i) {
    cells.emplace_back(storage);
    require(cells.back().get<0>() == 0);
    require(cells.back().get<1>() == 0);
    require(cells.back().get<2>() == 0);
  }
  require(snapshot.get<0>() == 3);
  cells[3] = snapshot;
  snapshot.get<0>() = 22;
  require(cells[3].get<0>() == 3);
  // Cross-region value assignment must retain the source snapshot's values.
  auto other_store = std::make_shared<Store>();
  Row other(other_store);
  other = snapshot;
  require(other.get<0>() == 22);
  snapshot.get<0>() = 23;
  require(other.get<0>() == 22);
  // Independent owner lifetime and self-assignment.
  cells.clear();
  storage.reset();
  const Row& same_row = other;
  other = same_row;
  require(other.get<0>() == 22);
  Row moved(std::move(other));
  require(moved.get<0>() == 22);
  require(snapshot.get<0>() == 23);
}
