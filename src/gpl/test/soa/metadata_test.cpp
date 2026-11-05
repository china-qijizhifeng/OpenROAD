// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include <memory>
#include <stdexcept>
#include <vector>

#include "placementColumns.h"

int main()
{
  using Columns = gpl::GCellColumns;
  using Row = gpl::SoaRow<Columns::Store>;
  auto pool = std::make_shared<Columns::Store>();
  Row live(pool);
  live.get<Columns::insts>().push_back(nullptr);
  live.get<Columns::gPins>().resize(3, nullptr);
  live.get<Columns::change>() = gpl::CellAreaChange::kRoutability;
  Row snapshot(live);
  live.get<Columns::insts>().clear();
  live.get<Columns::gPins>().pop_back();
  live.get<Columns::change>() = gpl::CellAreaChange::kNone;
  std::vector<Row> others;
  for (int i = 0; i < 10000; ++i) {
    others.emplace_back(pool);
  }
  if (snapshot.get<Columns::insts>().size() != 1
      || snapshot.get<Columns::gPins>().size() != 3
      || snapshot.get<Columns::change>() != gpl::CellAreaChange::kRoutability) {
    throw std::runtime_error("metadata snapshot changed after topology growth");
  }
  const size_t freed = others.back().slot();
  others.back().get<Columns::insts>().resize(5);
  others.pop_back();
  Row reuse(pool);
  if (reuse.slot() != freed || !reuse.get<Columns::insts>().empty()
      || reuse.get<Columns::change>() != gpl::CellAreaChange::kNone) {
    throw std::runtime_error("reused metadata row retains deleted topology");
  }
}
