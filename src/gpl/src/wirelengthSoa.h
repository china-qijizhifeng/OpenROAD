// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once

#include <algorithm>
#include <cassert>
#include <climits>
#include <span>

#include "placementColumns.h"
#include "point.h"

namespace gpl {

// Reindex only the WA intermediates, leaving stable geometry rows untouched.
// Called at topology boundaries; every connected pin is overwritten per force.
inline void resetPinWaColumns(GPinColumns::Store& pins,
                              std::span<const size_t> net_pin_rows)
{
  auto& wa_rows = pins.column<GPinColumns::waSlot>();
  std::fill(wa_rows.begin(), wa_rows.end(), wa_rows.size());
  for (size_t i = 0; i < net_pin_rows.size(); ++i) {
    assert(net_pin_rows[i] < wa_rows.size());
    assert(wa_rows[net_pin_rows[i]] == wa_rows.size());
    wa_rows[net_pin_rows[i]] = i;
  }
  size_t disconnected = net_pin_rows.size();
  for (size_t& row : wa_rows) {
    if (row == wa_rows.size()) {
      row = disconnected++;
    }
  }
  assert(disconnected == wa_rows.size());
  // Topology changes remap all WA rows, including disconnected pins and holes.
  // Connected pins overwrite every WA value/flag on each force update.
  auto clear = []<typename T>(std::vector<T>& column) {
    std::fill(column.begin(), column.end(), T{});
  };
  clear(pins.column<GPinColumns::minExpSumX>());
  clear(pins.column<GPinColumns::maxExpSumX>());
  clear(pins.column<GPinColumns::minExpSumY>());
  clear(pins.column<GPinColumns::maxExpSumY>());
  clear(pins.column<GPinColumns::hasMinExpSumX>());
  clear(pins.column<GPinColumns::hasMaxExpSumX>());
  clear(pins.column<GPinColumns::hasMinExpSumY>());
  clear(pins.column<GPinColumns::hasMaxExpSumY>());
}

template <typename Rows>
inline void updateNetBoxColumns(GNetColumns::Store& nets,
                                const GPinColumns::Store& pins,
                                size_t n,
                                const Rows& pin_rows)
{
  int lx = INT_MAX, ly = INT_MAX, ux = INT_MIN, uy = INT_MIN;
  const auto& x = pins.column<GPinColumns::cx>();
  const auto& y = pins.column<GPinColumns::cy>();
  for (size_t p : pin_rows) {
    lx = std::min(lx, x[p]);
    ly = std::min(ly, y[p]);
    ux = std::max(ux, x[p]);
    uy = std::max(uy, y[p]);
  }
  nets.at<GNetColumns::lx>(n) = lx;
  nets.at<GNetColumns::ly>(n) = ly;
  nets.at<GNetColumns::ux>(n) = ux;
  nets.at<GNetColumns::uy>(n) = uy;
}

// Use exactly OpenROAD's exponential approximation and operation ordering.
// DREAMPlace's exp()/coefficient rearrangement is intentionally not
// substituted.
inline float wirelengthFastExp(float value)
{
  value = 1.0f + value / 1024.0f;
  for (int i = 0; i < 10; ++i) {
    value *= value;
  }
  return value;
}

template <typename Rows>
inline void updateNetForceColumns(GNetColumns::Store& nets,
                                  GPinColumns::Store& pins,
                                  size_t n,
                                  const Rows& rows,
                                  float coef_x,
                                  float coef_y,
                                  float min_force)
{
  updateNetBoxColumns(nets, pins, n, rows);
  // Keep each net's ordered reductions in registers, then publish once.
  // Pin columns and net columns may otherwise alias through Store accessors.
  const int lx = nets.at<GNetColumns::lx>(n);
  const int ly = nets.at<GNetColumns::ly>(n);
  const int ux = nets.at<GNetColumns::ux>(n);
  const int uy = nets.at<GNetColumns::uy>(n);
  float min_sum_x = 0;
  float weighted_min_x = 0;
  float max_sum_x = 0;
  float weighted_max_x = 0;
  float min_sum_y = 0;
  float weighted_min_y = 0;
  float max_sum_y = 0;
  float weighted_max_y = 0;
  for (size_t p : rows) {
    const size_t wa = pins.at<GPinColumns::waSlot>(p);
    const int x = pins.at<GPinColumns::cx>(p);
    const int y = pins.at<GPinColumns::cy>(p);
    {
      const float exponent = (lx - x) * coef_x;
      const bool active = exponent > min_force;
      pins.at<GPinColumns::hasMinExpSumX>(wa) = active;
      const float value = active ? wirelengthFastExp(exponent) : 0.0f;
      pins.at<GPinColumns::minExpSumX>(wa) = value;
      if (active) {
        min_sum_x += value;
        weighted_min_x += x * value;
      }
    }
    {
      const float exponent = (x - ux) * coef_x;
      const bool active = exponent > min_force;
      pins.at<GPinColumns::hasMaxExpSumX>(wa) = active;
      const float value = active ? wirelengthFastExp(exponent) : 0.0f;
      pins.at<GPinColumns::maxExpSumX>(wa) = value;
      if (active) {
        max_sum_x += value;
        weighted_max_x += x * value;
      }
    }
    {
      const float exponent = (ly - y) * coef_y;
      const bool active = exponent > min_force;
      pins.at<GPinColumns::hasMinExpSumY>(wa) = active;
      const float value = active ? wirelengthFastExp(exponent) : 0.0f;
      pins.at<GPinColumns::minExpSumY>(wa) = value;
      if (active) {
        min_sum_y += value;
        weighted_min_y += y * value;
      }
    }
    {
      const float exponent = (y - uy) * coef_y;
      const bool active = exponent > min_force;
      pins.at<GPinColumns::hasMaxExpSumY>(wa) = active;
      const float value = active ? wirelengthFastExp(exponent) : 0.0f;
      pins.at<GPinColumns::maxExpSumY>(wa) = value;
      if (active) {
        max_sum_y += value;
        weighted_max_y += y * value;
      }
    }
  }
  nets.at<GNetColumns::waExpMinSumX>(n) = min_sum_x;
  nets.at<GNetColumns::waXExpMinSumX>(n) = weighted_min_x;
  nets.at<GNetColumns::waExpMaxSumX>(n) = max_sum_x;
  nets.at<GNetColumns::waXExpMaxSumX>(n) = weighted_max_x;
  nets.at<GNetColumns::waExpMinSumY>(n) = min_sum_y;
  nets.at<GNetColumns::waYExpMinSumY>(n) = weighted_min_y;
  nets.at<GNetColumns::waExpMaxSumY>(n) = max_sum_y;
  nets.at<GNetColumns::waYExpMaxSumY>(n) = weighted_max_y;
}

// Net-invariant loads and float denominator products are shared by all pins.
// The pin numerator retains the native float/double evaluation order.
struct NetGradientTerms
{
  float min_x, weighted_min_x, max_x, weighted_max_x;
  float min_y, weighted_min_y, max_y, weighted_max_y;
  float min_x_squared, max_x_squared, min_y_squared, max_y_squared;
  NetGradientTerms(const GNetColumns::Store& nets, size_t n)
      : min_x(nets.at<GNetColumns::waExpMinSumX>(n)),
        weighted_min_x(nets.at<GNetColumns::waXExpMinSumX>(n)),
        max_x(nets.at<GNetColumns::waExpMaxSumX>(n)),
        weighted_max_x(nets.at<GNetColumns::waXExpMaxSumX>(n)),
        min_y(nets.at<GNetColumns::waExpMinSumY>(n)),
        weighted_min_y(nets.at<GNetColumns::waYExpMinSumY>(n)),
        max_y(nets.at<GNetColumns::waExpMaxSumY>(n)),
        weighted_max_y(nets.at<GNetColumns::waYExpMaxSumY>(n)),
        min_x_squared(min_x * min_x),
        max_x_squared(max_x * max_x),
        min_y_squared(min_y * min_y),
        max_y_squared(max_y * max_y)
  {
  }
};

inline FloatPoint pinGradientColumns(const NetGradientTerms& terms,
                                     const GPinColumns::Store& pins,
                                     size_t p,
                                     float wlCoeffX,
                                     float wlCoeffY)
{
  const size_t wa = pins.at<GPinColumns::waSlot>(p);
  float gradientMinX = 0, gradientMinY = 0;
  float gradientMaxX = 0, gradientMaxY = 0;

  // min x
  if (pins.at<GPinColumns::hasMinExpSumX>(wa)) {
    // from Net.
    float waExpMinSumX = terms.min_x;
    float waXExpMinSumX = terms.weighted_min_x;

    gradientMinX
        = (waExpMinSumX
               * (pins.at<GPinColumns::minExpSumX>(wa)
                  * (1.0 - wlCoeffX * pins.at<GPinColumns::cx>(p)))
           + wlCoeffX * pins.at<GPinColumns::minExpSumX>(wa) * waXExpMinSumX)
          / terms.min_x_squared;
  }

  // max x
  if (pins.at<GPinColumns::hasMaxExpSumX>(wa)) {
    float waExpMaxSumX = terms.max_x;
    float waXExpMaxSumX = terms.weighted_max_x;

    gradientMaxX
        = (waExpMaxSumX
               * (pins.at<GPinColumns::maxExpSumX>(wa)
                  * (1.0 + wlCoeffX * pins.at<GPinColumns::cx>(p)))
           - wlCoeffX * pins.at<GPinColumns::maxExpSumX>(wa) * waXExpMaxSumX)
          / terms.max_x_squared;
  }

  // min y
  if (pins.at<GPinColumns::hasMinExpSumY>(wa)) {
    float waExpMinSumY = terms.min_y;
    float waYExpMinSumY = terms.weighted_min_y;

    gradientMinY
        = (waExpMinSumY
               * (pins.at<GPinColumns::minExpSumY>(wa)
                  * (1.0 - wlCoeffY * pins.at<GPinColumns::cy>(p)))
           + wlCoeffY * pins.at<GPinColumns::minExpSumY>(wa) * waYExpMinSumY)
          / terms.min_y_squared;
  }

  // max y
  if (pins.at<GPinColumns::hasMaxExpSumY>(wa)) {
    float waExpMaxSumY = terms.max_y;
    float waYExpMaxSumY = terms.weighted_max_y;

    gradientMaxY
        = (waExpMaxSumY
               * (pins.at<GPinColumns::maxExpSumY>(wa)
                  * (1.0 + wlCoeffY * pins.at<GPinColumns::cy>(p)))
           - wlCoeffY * pins.at<GPinColumns::maxExpSumY>(wa) * waYExpMaxSumY)
          / terms.max_y_squared;
  }

  return {gradientMinX - gradientMaxX, gradientMinY - gradientMaxY};
}
inline FloatPoint pinGradientColumns(const GNetColumns::Store& nets,
                                     const GPinColumns::Store& pins,
                                     size_t n,
                                     size_t p,
                                     float coef_x,
                                     float coef_y)
{
  return pinGradientColumns(NetGradientTerms(nets, n), pins, p, coef_x, coef_y);
}
}  // namespace gpl
