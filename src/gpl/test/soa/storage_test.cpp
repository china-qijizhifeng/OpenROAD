// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include <limits>
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
struct Cell
{
  struct Store : gpl::SoaStorage<int>
  {
    uint64_t revision = 0;
  };
  explicit Cell(std::shared_ptr<Store> store) : row(std::move(store)) {}
  Cell(gpl::SoaViewTag tag, void* store, size_t slot) : row(tag, store, slot) {}
  gpl::SoaRow<Store>& numericRow() { return row; }
  const gpl::SoaRow<Store>& numericRow() const { return row; }
  int value() const { return row.get<0>(); }
  void set(int value) { row.get<0>() = value; }
  gpl::SoaRow<Store> row;
};

void checkCompactLinks()
{
  using Links = gpl::SoaStorage<gpl::SoaRef<Cell>>;
  auto cells = std::make_shared<Cell::Store>();
  gpl::SoaPool<Cell> pool(cells);
  pool.emplace_back(cells);
  pool.back().set(91);
  gpl::SoaRef<Cell> reference(cells.get(), pool.front().numericRow().slot());
  auto links = std::make_shared<Links>();
  gpl::SoaRow<Links> row(links);
  row.get<0>() = reference;
  require(row.get<0>()->value() == 91);
  gpl::SoaRow<Links> copy(row);
  row.get<0>() = nullptr;
  require(!row.get<0>());
  require(std::as_const(copy).get<0>() == reference);
  pool.reserve(10000);
  require(copy.get<0>()->value() == 91);
}

void checkRowReferences()
{
  auto store = std::make_shared<Cell::Store>();
  gpl::SoaPool<Cell> pool(store);
  pool.reserve(100);
  require(store->capacity() == 100);
  require(store->allocatedBytes() >= 100 * sizeof(int));
  for (int i = 0; i < 100; ++i) {
    pool.emplace_back(store);
    pool.back().set(i);
  }
  const auto reference = pool.ref(3);
  gpl::SoaRows<Cell> adjacency{reference};
  pool.reserve(10000);
  for (int i = 100; i < 10000; ++i) {
    pool.emplace_back(store);
  }
  require(adjacency[0]->value() == 3);
  pool.swapRows(3, pool.size() - 1);
  require(pool.ref(pool.size() - 1) == reference);
  require(pool.indexOf(reference) == pool.size() - 1);
  // Snapshots clone state while ordinary references keep the row identity.
  gpl::SoaPool<Cell> snapshot(store);
  snapshot.push_back(*reference);
  reference->set(77);
  require(snapshot.front().value() == 3);
  pool.pop_back();
  pool.emplace_back(store);
  require(pool.back().value() == 0);
  pool.push_back(snapshot.front());
  snapshot.clear();
  require(pool.back().value() == 3);
  // Iteration creates views: writing through one must update the column.
  for (auto cell : pool) {
    cell.set(42);
  }
  require(pool.front().value() == 42);
  require(pool.back().value() == 42);
  // OpenMP may dereference a temporary iterator obtained by begin + index.
  auto view = *(pool.begin() + 3);
  view.set(17);
  require(pool[3].value() == 17);
}
}  // namespace

int main()
{
  checkRowReferences();
  checkCompactLinks();
  using Store = gpl::SoaStorage<int, float, int64_t>;
  using Row = gpl::SoaRow<Store>;
  auto storage = std::make_shared<Store>();
  // Reject unsupported batch sizes before trying to allocate any columns.
  bool rejected = false;
  try {
    storage->reserveRows(
        static_cast<size_t>(std::numeric_limits<uint32_t>::max()) + 1);
  } catch (const std::length_error&) {
    rejected = true;
  }
  require(rejected);
  require(storage->size() == 0 && storage->capacity() == 0);
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
