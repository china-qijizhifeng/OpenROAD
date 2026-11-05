// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once

#include <algorithm>
#include <climits>
#include <span>

#include "placementColumns.h"
#include "point.h"

namespace gpl {

inline void updateNetBoxColumns(GNetColumns::Store& nets,
                                const GPinColumns::Store& pins,
                                size_t n,
                                std::span<const size_t> pin_rows)
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

inline void updateNetForceColumns(GNetColumns::Store& nets,
                                  GPinColumns::Store& pins,
                                  size_t n,
                                  std::span<const size_t> rows,
                                  float coef_x,
                                  float coef_y,
                                  float min_force)
{
  updateNetBoxColumns(nets, pins, n, rows);
  nets.at<GNetColumns::waExpMinSumX>(n) = 0;
  nets.at<GNetColumns::waXExpMinSumX>(n) = 0;
  nets.at<GNetColumns::waExpMaxSumX>(n) = 0;
  nets.at<GNetColumns::waXExpMaxSumX>(n) = 0;
  nets.at<GNetColumns::waExpMinSumY>(n) = 0;
  nets.at<GNetColumns::waYExpMinSumY>(n) = 0;
  nets.at<GNetColumns::waExpMaxSumY>(n) = 0;
  nets.at<GNetColumns::waYExpMaxSumY>(n) = 0;
  for (size_t p : rows) {
    {
      const float exponent
          = (nets.at<GNetColumns::lx>(n) - pins.at<GPinColumns::cx>(p))
            * coef_x;
      const bool active = exponent > min_force;
      pins.at<GPinColumns::hasMinExpSumX>(p) = active;
      const float value = active ? wirelengthFastExp(exponent) : 0.0f;
      pins.at<GPinColumns::minExpSumX>(p) = value;
      if (active) {
        nets.at<GNetColumns::waExpMinSumX>(n) += value;
        nets.at<GNetColumns::waXExpMinSumX>(n)
            += pins.at<GPinColumns::cx>(p) * value;
      }
    }
    {
      const float exponent
          = (pins.at<GPinColumns::cx>(p) - nets.at<GNetColumns::ux>(n))
            * coef_x;
      const bool active = exponent > min_force;
      pins.at<GPinColumns::hasMaxExpSumX>(p) = active;
      const float value = active ? wirelengthFastExp(exponent) : 0.0f;
      pins.at<GPinColumns::maxExpSumX>(p) = value;
      if (active) {
        nets.at<GNetColumns::waExpMaxSumX>(n) += value;
        nets.at<GNetColumns::waXExpMaxSumX>(n)
            += pins.at<GPinColumns::cx>(p) * value;
      }
    }
    {
      const float exponent
          = (nets.at<GNetColumns::ly>(n) - pins.at<GPinColumns::cy>(p))
            * coef_y;
      const bool active = exponent > min_force;
      pins.at<GPinColumns::hasMinExpSumY>(p) = active;
      const float value = active ? wirelengthFastExp(exponent) : 0.0f;
      pins.at<GPinColumns::minExpSumY>(p) = value;
      if (active) {
        nets.at<GNetColumns::waExpMinSumY>(n) += value;
        nets.at<GNetColumns::waYExpMinSumY>(n)
            += pins.at<GPinColumns::cy>(p) * value;
      }
    }
    {
      const float exponent
          = (pins.at<GPinColumns::cy>(p) - nets.at<GNetColumns::uy>(n))
            * coef_y;
      const bool active = exponent > min_force;
      pins.at<GPinColumns::hasMaxExpSumY>(p) = active;
      const float value = active ? wirelengthFastExp(exponent) : 0.0f;
      pins.at<GPinColumns::maxExpSumY>(p) = value;
      if (active) {
        nets.at<GNetColumns::waExpMaxSumY>(n) += value;
        nets.at<GNetColumns::waYExpMaxSumY>(n)
            += pins.at<GPinColumns::cy>(p) * value;
      }
    }
  }
}

inline FloatPoint pinGradientColumns(const GNetColumns::Store& nets,
                                     const GPinColumns::Store& pins,
                                     size_t n,
                                     size_t p,
                                     float wlCoeffX,
                                     float wlCoeffY)
{
  float gradientMinX = 0, gradientMinY = 0;
  float gradientMaxX = 0, gradientMaxY = 0;

  // min x
  if (pins.at<GPinColumns::hasMinExpSumX>(p)) {
    // from Net.
    float waExpMinSumX = nets.at<GNetColumns::waExpMinSumX>(n);
    float waXExpMinSumX = nets.at<GNetColumns::waXExpMinSumX>(n);

    gradientMinX
        = (waExpMinSumX
               * (pins.at<GPinColumns::minExpSumX>(p)
                  * (1.0 - wlCoeffX * pins.at<GPinColumns::cx>(p)))
           + wlCoeffX * pins.at<GPinColumns::minExpSumX>(p) * waXExpMinSumX)
          / (waExpMinSumX * waExpMinSumX);
  }

  // max x
  if (pins.at<GPinColumns::hasMaxExpSumX>(p)) {
    float waExpMaxSumX = nets.at<GNetColumns::waExpMaxSumX>(n);
    float waXExpMaxSumX = nets.at<GNetColumns::waXExpMaxSumX>(n);

    gradientMaxX
        = (waExpMaxSumX
               * (pins.at<GPinColumns::maxExpSumX>(p)
                  * (1.0 + wlCoeffX * pins.at<GPinColumns::cx>(p)))
           - wlCoeffX * pins.at<GPinColumns::maxExpSumX>(p) * waXExpMaxSumX)
          / (waExpMaxSumX * waExpMaxSumX);
  }

  // min y
  if (pins.at<GPinColumns::hasMinExpSumY>(p)) {
    float waExpMinSumY = nets.at<GNetColumns::waExpMinSumY>(n);
    float waYExpMinSumY = nets.at<GNetColumns::waYExpMinSumY>(n);

    gradientMinY
        = (waExpMinSumY
               * (pins.at<GPinColumns::minExpSumY>(p)
                  * (1.0 - wlCoeffY * pins.at<GPinColumns::cy>(p)))
           + wlCoeffY * pins.at<GPinColumns::minExpSumY>(p) * waYExpMinSumY)
          / (waExpMinSumY * waExpMinSumY);
  }

  // max y
  if (pins.at<GPinColumns::hasMaxExpSumY>(p)) {
    float waExpMaxSumY = nets.at<GNetColumns::waExpMaxSumY>(n);
    float waYExpMaxSumY = nets.at<GNetColumns::waYExpMaxSumY>(n);

    gradientMaxY
        = (waExpMaxSumY
               * (pins.at<GPinColumns::maxExpSumY>(p)
                  * (1.0 + wlCoeffY * pins.at<GPinColumns::cy>(p)))
           - wlCoeffY * pins.at<GPinColumns::maxExpSumY>(p) * waYExpMaxSumY)
          / (waExpMaxSumY * waExpMaxSumY);
  }

  return {gradientMinX - gradientMaxX, gradientMinY - gradientMaxY};
}
}  // namespace gpl
