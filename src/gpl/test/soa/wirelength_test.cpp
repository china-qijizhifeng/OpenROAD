// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include <algorithm>
#include <array>
#include <bit>
#include <climits>
#include <cstdint>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include "wirelengthSoa.h"

namespace {
using namespace gpl;

void require(bool value)
{
  if (!value) {
    throw std::runtime_error("WA topology or numerical mismatch");
  }
}

bool same(float a, float b)
{
  return std::bit_cast<uint32_t>(a) == std::bit_cast<uint32_t>(b);
}

struct PinForce
{
  std::array<float, 4> value{};
  std::array<bool, 4> active{};
};

// Scalar native equations use geometry-row pins and individual net
// accumulators. They intentionally know nothing about WA storage slots or CSR
// remapping.
struct NativeNet
{
  int lx = INT_MAX, ly = INT_MAX, ux = INT_MIN, uy = INT_MIN;
  std::array<float, 4> sum{}, weighted{};
};

NativeNet nativeForce(const std::vector<int>& x,
                      const std::vector<int>& y,
                      std::span<const size_t> rows,
                      float cx,
                      float cy,
                      float minimum,
                      std::vector<PinForce>& forces)
{
  NativeNet net;
  for (size_t p : rows) {
    net.lx = std::min(net.lx, x[p]);
    net.ly = std::min(net.ly, y[p]);
    net.ux = std::max(net.ux, x[p]);
    net.uy = std::max(net.uy, y[p]);
  }
  for (size_t p : rows) {
    const std::array exponent{(net.lx - x[p]) * cx,
                              (x[p] - net.ux) * cx,
                              (net.ly - y[p]) * cy,
                              (y[p] - net.uy) * cy};
    for (size_t axis = 0; axis < 4; ++axis) {
      auto& pin = forces[p];
      pin.active[axis] = exponent[axis] > minimum;
      pin.value[axis] = 0;
      if (pin.active[axis]) {
        float value = 1.0f + exponent[axis] / 1024.0f;
        for (int square = 0; square < 10; ++square) {
          value *= value;
        }
        pin.value[axis] = value;
        net.sum[axis] += value;
        net.weighted[axis] += (axis < 2 ? x[p] : y[p]) * value;
      }
    }
  }
  return net;
}

FloatPoint nativeGradient(const NativeNet& net,
                          const PinForce& pin,
                          int x,
                          int y,
                          float cx,
                          float cy)
{
  std::array<float, 4> gradient{};
  for (size_t axis = 0; axis < 4; ++axis) {
    if (!pin.active[axis]) {
      continue;
    }
    const float coefficient = axis < 2 ? cx : cy;
    const int coordinate = axis < 2 ? x : y;
    const float sum = net.sum[axis], weighted = net.weighted[axis];
    const float value = pin.value[axis];
    if (axis % 2 == 0) {
      gradient[axis] = (sum * (value * (1.0 - coefficient * coordinate))
                        + coefficient * value * weighted)
                       / (sum * sum);
    } else {
      gradient[axis] = (sum * (value * (1.0 + coefficient * coordinate))
                        - coefficient * value * weighted)
                       / (sum * sum);
    }
  }
  return {gradient[0] - gradient[1], gradient[2] - gradient[3]};
}

void comparePin(const GPinColumns::Store& pins,
                size_t p,
                const PinForce& expected)
{
  const size_t wa = pins.at<GPinColumns::waSlot>(p);
  const std::array values{pins.at<GPinColumns::minExpSumX>(wa),
                          pins.at<GPinColumns::maxExpSumX>(wa),
                          pins.at<GPinColumns::minExpSumY>(wa),
                          pins.at<GPinColumns::maxExpSumY>(wa)};
  const std::array flags{pins.at<GPinColumns::hasMinExpSumX>(wa),
                         pins.at<GPinColumns::hasMaxExpSumX>(wa),
                         pins.at<GPinColumns::hasMinExpSumY>(wa),
                         pins.at<GPinColumns::hasMaxExpSumY>(wa)};
  for (size_t axis = 0; axis < 4; ++axis) {
    require(same(values[axis], expected.value[axis]));
    require(bool(flags[axis]) == expected.active[axis]);
  }
}

void compareNet(const GNetColumns::Store& nets,
                size_t n,
                const NativeNet& expected)
{
  require(nets.at<GNetColumns::lx>(n) == expected.lx);
  require(nets.at<GNetColumns::ly>(n) == expected.ly);
  require(nets.at<GNetColumns::ux>(n) == expected.ux);
  require(nets.at<GNetColumns::uy>(n) == expected.uy);
  const std::array sums{nets.at<GNetColumns::waExpMinSumX>(n),
                        nets.at<GNetColumns::waExpMaxSumX>(n),
                        nets.at<GNetColumns::waExpMinSumY>(n),
                        nets.at<GNetColumns::waExpMaxSumY>(n)};
  const std::array weighted{nets.at<GNetColumns::waXExpMinSumX>(n),
                            nets.at<GNetColumns::waXExpMaxSumX>(n),
                            nets.at<GNetColumns::waYExpMinSumY>(n),
                            nets.at<GNetColumns::waYExpMaxSumY>(n)};
  for (size_t axis = 0; axis < 4; ++axis) {
    require(same(sums[axis], expected.sum[axis]));
    require(same(weighted[axis], expected.weighted[axis]));
  }
}
}  // namespace

int main()
{
  using namespace gpl;
  constexpr size_t pin_count = 4096, net_count = 137;
  GPinColumns::Store pins;
  GNetColumns::Store nets;
  std::mt19937 rng(42);
  for (size_t p = 0; p < pin_count; ++p) {
    pins.allocate();
    pins.at<GPinColumns::cx>(p) = int(rng() % 2000001) - 1000000;
    pins.at<GPinColumns::cy>(p) = int(rng() % 2000001) - 1000000;
  }
  const auto x = pins.column<GPinColumns::cx>();
  const auto y = pins.column<GPinColumns::cy>();
  for (size_t n = 0; n < net_count; ++n) {
    nets.allocate();
  }
  // Released slots remain allocated columns, as in real topology changes.
  for (size_t p = 0; p < pin_count; p += 13) {
    pins.release(p);
  }
  std::vector<size_t> geometry_rows;
  for (size_t p = 0; p < pin_count; ++p) {
    if (p % 13 != 0) {
      geometry_rows.push_back(p);
    }
  }
  for (int topology = 0; topology < 5; ++topology) {
    std::shuffle(geometry_rows.begin(), geometry_rows.end(), rng);
    // Alternately disconnect/reconnect previously active pins; include an empty
    // topology to prove reset clears all rows without an ensuing force update.
    const size_t connected
        = topology == 2 ? 0 : geometry_rows.size() - topology * 127;
    const std::vector<size_t> entries(geometry_rows.begin(),
                                      geometry_rows.begin() + connected);
    resetPinWaColumns(pins, entries);
    std::vector<bool> seen(pin_count);
    for (size_t p = 0; p < pin_count; ++p) {
      const size_t wa = pins.at<GPinColumns::waSlot>(p);
      require(wa < pin_count && !seen[wa]);
      seen[wa] = true;
      comparePin(pins, p, {});
    }
    for (size_t i = 0; i < entries.size(); ++i) {
      require(pins.at<GPinColumns::waSlot>(entries[i]) == i);
    }
    std::vector<size_t> offsets{0};
    for (size_t n = 0; n < net_count; ++n) {
      // Empty, singleton, large and irregular nets all use permuted geometry.
      const size_t count = n + 1 == net_count
                               ? connected - offsets.back()
                               : std::min(n % 41, connected - offsets.back());
      offsets.push_back(offsets.back() + count);
    }
    for (int iteration = 0; iteration < 8; ++iteration) {
      const float cx = iteration % 3 == 0 ? 0.0f : 0.00001f;
      const float cy = cx * 2;
      const float minimum = iteration % 2 == 0 ? -10.0f : 0.0f;
      std::vector<PinForce> forces(pin_count);
      for (size_t n = 0; n < net_count; ++n) {
        auto rows = std::span<const size_t>(entries).subspan(
            offsets[n], offsets[n + 1] - offsets[n]);
        const auto expected = nativeForce(x, y, rows, cx, cy, minimum, forces);
        updateNetForceColumns(nets, pins, n, rows, cx, cy, minimum);
        compareNet(nets, n, expected);
        for (size_t p : rows) {
          comparePin(pins, p, forces[p]);
          const auto actual = pinGradientColumns(nets, pins, n, p, cx, cy);
          const auto gradient
              = nativeGradient(expected, forces[p], x[p], y[p], cx, cy);
          require(same(actual.x, gradient.x) && same(actual.y, gradient.y));
        }
      }
      for (size_t p = 0; p < pin_count; ++p) {
        comparePin(pins, p, forces[p]);
      }
      require(pins.column<GPinColumns::cx>() == x);
      require(pins.column<GPinColumns::cy>() == y);
    }
  }
}
