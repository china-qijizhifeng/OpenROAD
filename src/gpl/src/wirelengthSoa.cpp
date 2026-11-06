// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include "wirelengthSoa.h"

#include <algorithm>
#include <bit>
#include <cassert>
#include <limits>
#include <numeric>
#include <span>
#include <vector>

#include "nesterovBase.h"
#include "wirelengthSimd.h"

namespace gpl {
void NesterovBaseCommon::rebuildSoaTopology()
{
  // Adjacency row lists are the sole topology. Only WA data is laid out in
  // net traversal order; the temporary permutation is released at this
  // boundary.
  std::vector<size_t> wa_order;
  wa_order.reserve(gPinStor_.size());
  for (auto net : gNetStor_) {
    const auto& rows = net.getGPins().rows();
    wa_order.insert(wa_order.end(), rows.begin(), rows.end());
  }
  resetPinWaColumns(*pin_columns_, wa_order);
  auto& gradient_slots = pin_columns_->column<GPinColumns::gradientSlot>();
  std::fill(gradient_slots.begin(),
            gradient_slots.end(),
            std::numeric_limits<uint32_t>::max());
  size_t count = 0;
  for (auto cell : gCellStor_) {
    cell_columns_->at<GCellColumns::gradientBegin>(cell.numericRow().slot())
        = count;
    for (size_t p : cell.gPins().rows()) {
      if (gradient_slots[p] != std::numeric_limits<uint32_t>::max()) {
        throw std::logic_error("GPL pin belongs to multiple gradient outputs");
      }
      gradient_slots[p] = count++;
    }
  }
  pin_gradients_.resize(count);
  pin_gradients_valid_ = false;
  if (log_->debugCheck(utl::GPL, "soaMemory", 1)) {
    auto report = [&](const char* name, const auto& store, size_t index_bytes) {
      debugPrint(
          log_,
          utl::GPL,
          "soaMemory",
          1,
          "{} slots={} live_rows={} capacity={} column_bytes={} index_bytes={}",
          name,
          store.size(),
          store.liveSize(),
          store.capacity(),
          store.allocatedBytes(),
          index_bytes);
    };
    report("Instance", *pbc_->instanceColumns(), 0);
    report("Pin", *pbc_->pinColumns(), 0);
    report("Net", *pbc_->netColumns(), 0);
    report("GCell", *cell_columns_, gCellStor_.allocatedBytes());
    report("GPin", *pin_columns_, gPinStor_.allocatedBytes());
    report("GNet", *net_columns_, gNetStor_.allocatedBytes());
    debugPrint(
        log_,
        utl::GPL,
        "soaMemory",
        1,
        "consumer_pins={} gradient_bytes={} physical_row_map_bytes={}",
        count,
        (pin_gradients_.x().capacity() + pin_gradients_.y().capacity())
            * sizeof(float),
        (gCellMap_.capacity() + gPinMap_.capacity() + gNetMap_.capacity())
            * sizeof(uint32_t));
  }
}

void NesterovBaseCommon::updateDensityPinLocations(const GCellRowList& cells)
{
  pin_gradients_valid_ = false;
  const auto& c = *cells.columns;
  const auto& lx = c.column<GCellColumns::dLx>();
  const auto& ly = c.column<GCellColumns::dLy>();
  const auto& ux = c.column<GCellColumns::dUx>();
  const auto& uy = c.column<GCellColumns::dUy>();
  auto& px = pin_columns_->column<GPinColumns::cx>();
  auto& py = pin_columns_->column<GPinColumns::cy>();
  const auto& ox = pin_columns_->column<GPinColumns::offsetCx>();
  const auto& oy = pin_columns_->column<GPinColumns::offsetCy>();
#pragma omp parallel for num_threads(num_threads_)
  for (size_t r : cells.rows) {
    const int x = (ux[r] + lx[r]) / 2, y = (uy[r] + ly[r]) / 2;
    for (size_t p : c.at<GCellColumns::gPins>(r).rows()) {
      px[p] = x + ox[p];
      py[p] = y + oy[p];
    }
  }
}

void NesterovBaseCommon::updateWireLengthForceWA_soa(float coef_x, float coef_y)
{
  const bool use_avx2 = useAvx2Wirelength();
#pragma omp parallel for num_threads(num_threads_)
  for (size_t i = 0; i < gNetStor_.size(); ++i) {
    const size_t n = gNetStor_.rows()[i];
    const auto& rows = net_columns_->at<GNetColumns::gPins>(n).rows();
    updateNetForceColumns(*net_columns_,
                          *pin_columns_,
                          n,
                          rows,
                          coef_x,
                          coef_y,
                          nbVars_.minWireLengthForceBar);
    // Reuse the just-written net/pin data while it is local. Cell gather then
    // reads two gradient columns instead of repeating cross-net arithmetic.
    const NetGradientTerms terms(*net_columns_, n);
    const float weight = (net_columns_->at<GNetColumns::timingWeight>(n)
                          + net_columns_->at<GNetColumns::powerWeight>(n))
                         * net_columns_->at<GNetColumns::customWeight>(n);
    scatterNetPinGradients(terms,
                           *pin_columns_,
                           rows,
                           pin_gradients_,
                           coef_x,
                           coef_y,
                           weight,
                           use_avx2);
  }
  pin_gradient_weight_revision_ = net_columns_->weight_revision;
  pin_gradient_coef_x_ = coef_x;
  pin_gradient_coef_y_ = coef_y;
  pin_gradients_valid_ = true;
}

void NesterovBaseCommon::getWireLengthGradientsWA_cpu(
    const std::vector<GCellHandle>& cells,
    PointColumns& out,
    float coef_x,
    float coef_y)
{
  const bool debug = log_->debugCheck(utl::GPL, "getGradientWA", 1)
                     || log_->debugCheck(utl::GPL, "getGradientWAPin", 1);
  const bool use_cached
      = !debug && pin_gradients_valid_
        && pin_gradient_weight_revision_ == net_columns_->weight_revision
        && std::bit_cast<uint32_t>(coef_x)
               == std::bit_cast<uint32_t>(pin_gradient_coef_x_)
        && std::bit_cast<uint32_t>(coef_y)
               == std::bit_cast<uint32_t>(pin_gradient_coef_y_);
#pragma omp parallel for num_threads(num_threads_)
  for (size_t i = 0; i < cells.size(); ++i) {
    if (!use_cached) {
      out[i] = getWireLengthGradientWA(cells[i], coef_x, coef_y);
      continue;
    }
    const SoaRef<GCell> cell = cells[i];
    FloatPoint sum;
    if (cell.storage() == cell_columns_.get()) {
      const size_t begin
          = cell_columns_->at<GCellColumns::gradientBegin>(cell.row());
      const size_t end
          = begin + cell_columns_->at<GCellColumns::gPins>(cell.row()).size();
      for (size_t p = begin; p < end; ++p) {
        const FloatPoint gradient = pin_gradients_[p];
        sum.x += gradient.x;
        sum.y += gradient.y;
      }
    }
    out[i] = sum;
  }
}

FloatPoint NesterovBaseCommon::getWireLengthGradientWA_soa(SoaRef<GCell> cell,
                                                           float coef_x,
                                                           float coef_y) const
{
  if (cell->isFiller()) {
    return {};
  }
  const size_t c = cell->numericRow().slot();
  const auto& nets = *net_columns_;
  FloatPoint sum;
  for (size_t p : cell_columns_->at<GCellColumns::gPins>(c).rows()) {
    const size_t n = pin_columns_->at<GPinColumns::gNet>(p).row();
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
  for (size_t i = 0; i < gNetStor_.size(); ++i) {
    const size_t n = gNetStor_.rows()[i];
    const auto& rows = net_columns_->at<GNetColumns::gPins>(n).rows();
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
