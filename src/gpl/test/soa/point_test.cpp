// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include "point.h"

#include <stdexcept>
#include <utility>

int main()
{
  gpl::PointColumns cur;
  cur.resize(1000);
  for (size_t i = 0; i < cur.size(); ++i) {
    cur[i] = gpl::FloatPoint(i * 0.25f, -float(i));
  }
  auto snapshot = cur;
  cur[3] = cur[7];
  cur[7].x += 1;
  if (cur[3].x != 1.75f || snapshot[3].x != 0.75f) {
    throw std::runtime_error("row or rollback snapshot aliases a source row");
  }
  cur.swapRows(3, 999);
  cur.pop_back();
  cur.emplace_back();
  if (cur[3].x != 249.75f || cur[999].x != 0 || cur[999].y != 0) {
    throw std::runtime_error("optimizer swap-remove/insert lost row values");
  }
  auto next = cur;
  next.x()[12] = 123;
  std::swap(cur, next);
  const auto& read = cur;
  if (read[12].x != 123 || next[12].x != 3) {
    throw std::runtime_error("backtracking state swap lost column ownership");
  }
  // IO placement uses a conditional value/reference expression.
  cur[0] = true ? cur[1] : gpl::FloatPoint{};
  cur[1] = false ? cur[0] : gpl::FloatPoint(7, 8);
  if (cur[0].x != 0.25f || cur[1].y != 8) {
    throw std::runtime_error("conditional row assignment failed");
  }
}
