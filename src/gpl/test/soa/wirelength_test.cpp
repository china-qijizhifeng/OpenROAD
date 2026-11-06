// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include <algorithm>
#include <array>
#include <bit>
#include <cfenv>
#include <climits>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "wirelengthSimd.h"
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

size_t simd_calls = 0;

// The oracle uses geometry rows; production kernels read contiguous WA slots
// and scatter to a separately permuted output. Untouched output slots retain
// distinctive values, including negative zero, to detect stray stores.
void compareScatter(const GPinColumns::Store& pins,
                    const GNetColumns::Store& nets,
                    size_t n,
                    std::span<const size_t> rows,
                    const NativeNet& net,
                    const std::vector<PinForce>& forces,
                    float cx,
                    float cy,
                    float weight,
                    bool use_simd)
{
  std::vector<uint32_t> scatter_rows(rows.begin(), rows.end());
  PointColumns expected;
  expected.resize(pins.size() + 17, {-0.0f, 123.25f});
  for (size_t p : rows) {
    const uint32_t slot = pins.at<GPinColumns::gradientSlot>(p);
    if (slot == UINT32_MAX) {
      continue;
    }
    const auto gradient = nativeGradient(net,
                                         forces[p],
                                         pins.at<GPinColumns::cx>(p),
                                         pins.at<GPinColumns::cy>(p),
                                         cx,
                                         cy);
    expected[slot] = {gradient.x * weight, gradient.y * weight};
  }
  const NetGradientTerms terms(nets, n);
  auto check = [&](auto kernel, const char* name) {
    PointColumns actual;
    actual.resize(expected.size(), {-0.0f, 123.25f});
    std::fenv_t saved;
    require(std::fegetenv(&saved) == 0);
    require(std::feclearexcept(FE_ALL_EXCEPT) == 0);
    kernel(terms, pins, scatter_rows, actual, cx, cy, weight);
    const int exceptions = std::fetestexcept(FE_INVALID | FE_DIVBYZERO);
    require(std::fesetenv(&saved) == 0);
    if (exceptions != 0) {
      throw std::runtime_error(std::string(name)
                               + " raised invalid/divide-by-zero, pins="
                               + std::to_string(rows.size()));
    }
    for (size_t slot = 0; slot < actual.size(); ++slot) {
      if (!same(actual[slot].x, expected[slot].x)
          || !same(actual[slot].y, expected[slot].y)) {
        throw std::runtime_error(std::string(name) + " scatter mismatch, pins="
                                 + std::to_string(rows.size())
                                 + " slot=" + std::to_string(slot));
      }
    }
  };
  check(scatterNetPinGradientsScalar, "scalar");
#if GPL_WIRELENGTH_AVX2
  if (use_simd) {
    check(scatterNetPinGradientsAvx2, "AVX2");
    ++simd_calls;
  }
#else
  require(!use_simd);
#endif
}

void scatterBoundaryCases(bool use_simd)
{
  std::mt19937 rng(739);
  for (size_t count : {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                       13, 14, 15, 16, 17, 31, 32, 33, 63, 64, 65, 129}) {
    GPinColumns::Store pins;
    GNetColumns::Store nets;
    nets.allocate();
    const size_t geometry_count = count + 11;
    std::vector<size_t> permutation(geometry_count);
    std::iota(permutation.begin(), permutation.end(), 0);
    std::shuffle(permutation.begin(), permutation.end(), rng);
    for (size_t p = 0; p < geometry_count; ++p) {
      pins.allocate();
      // Odd coordinates beyond float's exact integer range expose premature
      // conversion to double or a different coefficient multiplication order.
      pins.at<GPinColumns::cx>(p)
          = p % 3 == 0 ? 20000001 : int(rng() % 2000001) - 1000000;
      pins.at<GPinColumns::cy>(p)
          = p % 3 == 1 ? -20000001 : int(rng() % 2000001) - 1000000;
      pins.at<GPinColumns::gradientSlot>(p) = permutation[p] + 8;
    }
    std::shuffle(permutation.begin(), permutation.end(), rng);
    resetPinWaColumns(pins, permutation);
    // An unaligned WA start catches assumptions hidden by aligned test data.
    const auto rows = std::span<const size_t>(permutation).subspan(3, count);
    const auto x = pins.column<GPinColumns::cx>();
    const auto y = pins.column<GPinColumns::cy>();
    for (int consumers = 0; consumers < 3; ++consumers) {
      for (size_t p = 0; p < geometry_count; ++p) {
        pins.at<GPinColumns::gradientSlot>(p)
            = consumers == 2 || (consumers == 1 && p % 3 == 0) ? UINT32_MAX
                                                               : p + 8;
      }
      for (const auto coefficient : {FloatPoint{0, 0},
                                     FloatPoint{0.00001f, 0.00002f},
                                     FloatPoint{0.00000001f, 0.00000003f},
                                     FloatPoint{0, 0.00003f},
                                     FloatPoint{0.00003f, 0}}) {
        for (float minimum : {-10.0f, -1.0f, 0.0f}) {
          std::vector<PinForce> forces(geometry_count);
          const auto net = nativeForce(
              x, y, rows, coefficient.x, coefficient.y, minimum, forces);
          updateNetForceColumns(
              nets, pins, 0, rows, coefficient.x, coefficient.y, minimum);
          for (float weight : {0.0f, -0.0f, 1.0f, -1.0f, 0.375f, 7.25f}) {
            compareScatter(pins,
                           nets,
                           0,
                           rows,
                           net,
                           forces,
                           coefficient.x,
                           coefficient.y,
                           weight,
                           use_simd);
          }
        }
      }
    }
  }
}
}  // namespace

int main()
{
  using namespace gpl;
  const bool use_simd = hasAvx2Wirelength();
  scatterBoundaryCases(use_simd);
  constexpr size_t pin_count = 4096, net_count = 137;
  GPinColumns::Store pins;
  GNetColumns::Store nets;
  std::mt19937 rng(42);
  for (size_t p = 0; p < pin_count; ++p) {
    pins.allocate();
    pins.at<GPinColumns::cx>(p) = int(rng() % 2000001) - 1000000;
    pins.at<GPinColumns::cy>(p) = int(rng() % 2000001) - 1000000;
  }
  for (size_t n = 0; n < net_count; ++n) {
    nets.allocate();
  }
  std::vector<uint32_t> gradient_slots(pin_count);
  std::iota(gradient_slots.begin(), gradient_slots.end(), 8);
  std::mt19937 slot_rng(918);
  std::shuffle(gradient_slots.begin(), gradient_slots.end(), slot_rng);
  for (size_t p = 0; p < pin_count; ++p) {
    pins.at<GPinColumns::gradientSlot>(p)
        = p % 7 == 0 ? UINT32_MAX : gradient_slots[p];
  }
  // Deleted rows release their topology; the column allocation remains
  // reusable.
  for (size_t p = 0; p < pin_count; p += 13) {
    pins.release(p);
  }
  const auto x = pins.column<GPinColumns::cx>();
  const auto y = pins.column<GPinColumns::cy>();
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
        compareScatter(pins,
                       nets,
                       n,
                       rows,
                       expected,
                       forces,
                       cx,
                       cy,
                       iteration % 2 == 0 ? 0.375f : 7.25f,
                       use_simd);
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
  if (use_simd) {
    require(simd_calls > 0);
    std::cout << "PASS AVX2 native-oracle bitwise scatter: " << simd_calls
              << " calls\n";
  } else {
    std::cout << "SKIP AVX2 scatter: runtime support unavailable; "
                 "scalar native-oracle bitwise scatter passed\n";
  }
}
