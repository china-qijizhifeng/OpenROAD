// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include "wirelengthSoa.h"

#include <algorithm>
#include <numeric>
#include <span>
#include <vector>

#include "nesterovBase.h"

namespace gpl {
namespace {
// CSR is keyed by stable numeric row, not the movable object-vector index.
// Fill order is the original pin order to retain floating-point summation.
template <typename Objects, typename Pins>
void buildPinIndex(const Objects& objects,
                   size_t row_count,
                   Pins get_pins,
                   std::vector<size_t>& offsets,
                   std::vector<size_t>& entries)
{
  offsets.assign(row_count + 1, 0);
  for (const auto& object : objects) {
    offsets[object.numericRow().slot() + 1] = get_pins(object).size();
  }
  std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
  entries.resize(offsets.back());
  for (const auto& object : objects) {
    size_t index = offsets[object.numericRow().slot()];
    for (const GPin* pin : get_pins(object)) {
      entries[index++] = pin->numericRow().slot();
    }
  }
}
}  // namespace

void NesterovBaseCommon::rebuildSoaTopology()
{
  buildPinIndex(
      gNetStor_,
      net_columns_->column<GNetColumns::lx>().size(),
      [](const GNet& net) -> const auto& { return net.getGPins(); },
      net_pin_offsets_,
      net_pin_rows_);
  buildPinIndex(
      gCellStor_,
      cell_columns_->column<GCellColumns::lx>().size(),
      [](const GCell& cell) -> const auto& { return cell.gPins(); },
      cell_pin_offsets_,
      cell_pin_rows_);
  net_rows_.clear();
  net_rows_.reserve(gNetStor_.size());
  for (const auto& net : gNetStor_) {
    net_rows_.push_back(net.numericRow().slot());
  }
  pin_net_rows_.resize(pin_columns_->column<GPinColumns::cx>().size());
  for (const auto& pin : gPinStor_) {
    if (pin.getGNet()) {
      pin_net_rows_[pin.numericRow().slot()]
          = pin.getGNet()->numericRow().slot();
    }
  }
}

void NesterovBaseCommon::updateWireLengthForceWA_soa(float coef_x, float coef_y)
{
  // Disconnected pins retain the native cleared-WA state. Clearing columns
  // also covers holes left by swap-remove; no per-pin object pass is needed.
  auto clear = []<typename T>(std::vector<T>& column) {
    std::fill(column.begin(), column.end(), T{});
  };
  clear(pin_columns_->column<GPinColumns::minExpSumX>());
  clear(pin_columns_->column<GPinColumns::maxExpSumX>());
  clear(pin_columns_->column<GPinColumns::minExpSumY>());
  clear(pin_columns_->column<GPinColumns::maxExpSumY>());
  clear(pin_columns_->column<GPinColumns::hasMinExpSumX>());
  clear(pin_columns_->column<GPinColumns::hasMaxExpSumX>());
  clear(pin_columns_->column<GPinColumns::hasMinExpSumY>());
  clear(pin_columns_->column<GPinColumns::hasMaxExpSumY>());
#pragma omp parallel for num_threads(num_threads_)
  for (size_t i = 0; i < net_rows_.size(); ++i) {
    const size_t n = net_rows_[i];
    const auto rows
        = std::span<const size_t>(net_pin_rows_)
              .subspan(net_pin_offsets_[n],
                       net_pin_offsets_[n + 1] - net_pin_offsets_[n]);
    updateNetForceColumns(*net_columns_,
                          *pin_columns_,
                          n,
                          rows,
                          coef_x,
                          coef_y,
                          nbVars_.minWireLengthForceBar);
  }
}

FloatPoint NesterovBaseCommon::getWireLengthGradientWA_soa(const GCell* cell,
                                                           float coef_x,
                                                           float coef_y) const
{
  if (cell->isFiller()) {
    return {};
  }
  const size_t c = cell->numericRow().slot();
  const auto& nets = *net_columns_;
  FloatPoint sum;
  for (size_t i = cell_pin_offsets_[c]; i < cell_pin_offsets_[c + 1]; ++i) {
    const size_t p = cell_pin_rows_[i], n = pin_net_rows_[p];
    const auto gradient
        = pinGradientColumns(nets, *pin_columns_, n, p, coef_x, coef_y);
    const float weight = (nets.at<GNetColumns::timingWeight>(n)
                          + nets.at<GNetColumns::powerWeight>(n))
                         * nets.at<GNetColumns::customWeight>(n);
    sum.x += gradient.x * weight;
    sum.y += gradient.y * weight;
  }
  return sum;
}
int64_t NesterovBaseCommon::getHpwlColumns()
{
  int64_t hpwl = 0;
#pragma omp parallel for num_threads(num_threads_) reduction(+ : hpwl)
  for (size_t i = 0; i < net_rows_.size(); ++i) {
    const size_t n = net_rows_[i];
    const auto rows
        = std::span<const size_t>(net_pin_rows_)
              .subspan(net_pin_offsets_[n],
                       net_pin_offsets_[n + 1] - net_pin_offsets_[n]);
    auto& nets = *net_columns_;
    updateNetBoxColumns(nets, *pin_columns_, n, rows);
    if (nets.at<GNetColumns::ux>(n) >= nets.at<GNetColumns::lx>(n)) {
      const int64_t lx = nets.at<GNetColumns::lx>(n);
      const int64_t ly = nets.at<GNetColumns::ly>(n);
      const int64_t ux = nets.at<GNetColumns::ux>(n);
      const int64_t uy = nets.at<GNetColumns::uy>(n);
      hpwl += ux - lx + uy - ly;
    }
  }
  return hpwl;
}
}  // namespace gpl
