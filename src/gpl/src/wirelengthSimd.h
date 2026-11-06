// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <span>

#include "wirelengthSoa.h"

#if (defined(__x86_64__) || defined(__i386__)) && defined(__GNUC__)
#include <immintrin.h>
#define GPL_WIRELENGTH_AVX2 1
#else
#define GPL_WIRELENGTH_AVX2 0
#endif

namespace gpl {

inline bool hasAvx2Wirelength()
{
#if GPL_WIRELENGTH_AVX2
  static const bool available = __builtin_cpu_supports("avx2");
  return available;
#else
  return false;
#endif
}

// Explicit selection until full-placement measurements justify a default.
// Unsupported CPUs always retain the portable scalar path.
inline bool useAvx2Wirelength()
{
  static const bool enabled = [] {
    const char* mode = std::getenv("GPL_WIRELENGTH_SIMD");
    return mode != nullptr && std::strcmp(mode, "avx2") == 0
           && hasAvx2Wirelength();
  }();
  return enabled;
}

inline void scatterNetPinGradientsScalar(const NetGradientTerms& terms,
                                         const GPinColumns::Store& pins,
                                         std::span<const uint32_t> rows,
                                         PointColumns& out,
                                         float coef_x,
                                         float coef_y,
                                         float weight)
{
  for (size_t p : rows) {
    const uint32_t output = pins.at<GPinColumns::gradientSlot>(p);
    if (output != std::numeric_limits<uint32_t>::max()) {
      const auto gradient = pinGradientColumns(terms, pins, p, coef_x, coef_y);
      out[output] = FloatPoint{gradient.x * weight, gradient.y * weight};
    }
  }
}

#if GPL_WIRELENGTH_AVX2
namespace detail {

// Four independent pins, with precisely the scalar mixed-precision sequence:
// coef*coordinate and coef*exp*weighted are float, the 1.0 expression and
// division are double, then each directional derivative rounds to float.
// Inactive lanes use zero inputs and denominator one, never an invalid divide.
template <bool Maximum>
__attribute__((target("avx2"), always_inline)) inline __m128 gradientAxis4(
    __m128 coordinate,
    __m128 exponential,
    __m128 active,
    float coefficient,
    float sum,
    float weighted,
    float squared)
{
  if (_mm_movemask_ps(active) == 0) {
    return _mm_setzero_ps();
  }
  const __m128 coef = _mm_and_ps(active, _mm_set1_ps(coefficient));
  const __m128 value = _mm_and_ps(active, exponential);
  const __m128 coord = _mm_and_ps(active, coordinate);
  const __m128 net_sum = _mm_and_ps(active, _mm_set1_ps(sum));
  const __m128 net_weighted = _mm_and_ps(active, _mm_set1_ps(weighted));
  const __m128 denominator
      = _mm_or_ps(_mm_and_ps(active, _mm_set1_ps(squared)),
                  _mm_andnot_ps(active, _mm_set1_ps(1.0f)));
  const __m256d scaled_coord = _mm256_cvtps_pd(_mm_mul_ps(coef, coord));
  const __m256d inner = Maximum
                            ? _mm256_add_pd(_mm256_set1_pd(1.0), scaled_coord)
                            : _mm256_sub_pd(_mm256_set1_pd(1.0), scaled_coord);
  const __m256d left = _mm256_mul_pd(
      _mm256_cvtps_pd(net_sum), _mm256_mul_pd(_mm256_cvtps_pd(value), inner));
  const __m256d right
      = _mm256_cvtps_pd(_mm_mul_ps(_mm_mul_ps(coef, value), net_weighted));
  const __m256d numerator
      = Maximum ? _mm256_sub_pd(left, right) : _mm256_add_pd(left, right);
  return _mm_and_ps(
      active,
      _mm256_cvtpd_ps(_mm256_div_pd(numerator, _mm256_cvtps_pd(denominator))));
}

__attribute__((target("avx2"), always_inline)) inline __m128 activePins4(
    const uint8_t* flags,
    __m128 consumers)
{
  uint32_t bytes;
  std::memcpy(&bytes, flags, sizeof(bytes));
  const __m128i values = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(bytes));
  return _mm_and_ps(
      consumers,
      _mm_castsi128_ps(_mm_cmpgt_epi32(values, _mm_setzero_si128())));
}

}  // namespace detail

// rows is a net's canonical pin order. WA slots are contiguous in this order;
// geometry rows and output slots may be scattered. No net or cell reduction
// is vectorized. Only this function requires AVX2, not its callers or Eigen.
__attribute__((target("avx2"))) inline void scatterNetPinGradientsAvx2(
    const NetGradientTerms& terms,
    const GPinColumns::Store& pins,
    std::span<const uint32_t> rows,
    PointColumns& out,
    float coef_x,
    float coef_y,
    float weight)
{
  const auto& x = pins.column<GPinColumns::cx>();
  const auto& y = pins.column<GPinColumns::cy>();
  const auto& outputs = pins.column<GPinColumns::gradientSlot>();
  const auto& wa_slots = pins.column<GPinColumns::waSlot>();
  size_t i = 0;
  for (; i + 4 <= rows.size(); i += 4) {
    const uint32_t p0 = rows[i], p1 = rows[i + 1], p2 = rows[i + 2],
                   p3 = rows[i + 3];
    const size_t wa = wa_slots[p0];
    assert(wa_slots[p1] == wa + 1 && wa_slots[p2] == wa + 2
           && wa_slots[p3] == wa + 3);
    const __m128 consumers
        = _mm_castsi128_ps(_mm_setr_epi32(outputs[p0] != UINT32_MAX ? -1 : 0,
                                          outputs[p1] != UINT32_MAX ? -1 : 0,
                                          outputs[p2] != UINT32_MAX ? -1 : 0,
                                          outputs[p3] != UINT32_MAX ? -1 : 0));
    if (_mm_movemask_ps(consumers) == 0) {
      continue;
    }
    const __m128 cx
        = _mm_cvtepi32_ps(_mm_setr_epi32(x[p0], x[p1], x[p2], x[p3]));
    const __m128 cy
        = _mm_cvtepi32_ps(_mm_setr_epi32(y[p0], y[p1], y[p2], y[p3]));
    const __m128 min_x = detail::gradientAxis4<false>(
        cx,
        _mm_loadu_ps(pins.column<GPinColumns::minExpSumX>().data() + wa),
        detail::activePins4(
            pins.column<GPinColumns::hasMinExpSumX>().data() + wa, consumers),
        coef_x,
        terms.min_x,
        terms.weighted_min_x,
        terms.min_x_squared);
    const __m128 max_x = detail::gradientAxis4<true>(
        cx,
        _mm_loadu_ps(pins.column<GPinColumns::maxExpSumX>().data() + wa),
        detail::activePins4(
            pins.column<GPinColumns::hasMaxExpSumX>().data() + wa, consumers),
        coef_x,
        terms.max_x,
        terms.weighted_max_x,
        terms.max_x_squared);
    const __m128 min_y = detail::gradientAxis4<false>(
        cy,
        _mm_loadu_ps(pins.column<GPinColumns::minExpSumY>().data() + wa),
        detail::activePins4(
            pins.column<GPinColumns::hasMinExpSumY>().data() + wa, consumers),
        coef_y,
        terms.min_y,
        terms.weighted_min_y,
        terms.min_y_squared);
    const __m128 max_y = detail::gradientAxis4<true>(
        cy,
        _mm_loadu_ps(pins.column<GPinColumns::maxExpSumY>().data() + wa),
        detail::activePins4(
            pins.column<GPinColumns::hasMaxExpSumY>().data() + wa, consumers),
        coef_y,
        terms.max_y,
        terms.weighted_max_y,
        terms.max_y_squared);
    float gx[4], gy[4];
    _mm_storeu_ps(gx,
                  _mm_mul_ps(_mm_sub_ps(min_x, max_x), _mm_set1_ps(weight)));
    _mm_storeu_ps(gy,
                  _mm_mul_ps(_mm_sub_ps(min_y, max_y), _mm_set1_ps(weight)));
    for (size_t lane = 0; lane < 4; ++lane) {
      const uint32_t output = outputs[rows[i + lane]];
      if (output != UINT32_MAX) {
        out[output] = FloatPoint{gx[lane], gy[lane]};
      }
    }
  }
  scatterNetPinGradientsScalar(
      terms, pins, rows.subspan(i), out, coef_x, coef_y, weight);
}
#endif

inline void scatterNetPinGradients(const NetGradientTerms& terms,
                                   const GPinColumns::Store& pins,
                                   std::span<const uint32_t> rows,
                                   PointColumns& out,
                                   float coef_x,
                                   float coef_y,
                                   float weight,
                                   [[maybe_unused]] bool use_avx2)
{
#if GPL_WIRELENGTH_AVX2
  if (use_avx2 && rows.size() >= 4) {
    scatterNetPinGradientsAvx2(terms, pins, rows, out, coef_x, coef_y, weight);
    return;
  }
#endif
  scatterNetPinGradientsScalar(terms, pins, rows, out, coef_x, coef_y, weight);
}

}  // namespace gpl
