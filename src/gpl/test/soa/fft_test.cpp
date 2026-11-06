// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include "fft.h"

#include <array>
#include <cstring>
#include <random>
#include <stdexcept>
#include <vector>

int main()
{
  std::mt19937 random(42);
  // Cover both axes, tiny fallback, and parallel narrow rectangular grids.
  for (const auto [nx, ny] : {std::pair{2, 2},
                              {8, 4},
                              {64, 64},
                              {128, 32},
                              {32, 128},
                              {256, 256},
                              {2048, 2},
                              {2, 2048}}) {
    gpl::FFT reference(nx, ny, 10.375f, 13.25f);
    gpl::FFT candidate(nx, ny, 10.375f, 13.25f);
    const size_t count = nx * ny;
    std::vector<float> density(count);
    std::array<std::vector<float>, 3> expected, actual;
    for (auto& c : expected) {
      c.resize(count);
    }
    for (auto& c : actual) {
      c.resize(count);
    }
    // Change thread counts and inputs between calls to exercise table/scratch
    // reuse, including initialization on the parallel path.
    for (int threads : {8, 1, 2, 8}) {
      for (float& v : density) {
        v = static_cast<int>(random() % 20001) / 97.0f - 100.0f;
      }
      candidate.setNumThreads(threads);
      reference.solveColumns(density, expected[0], expected[1], expected[2]);
      candidate.solveColumns(density, actual[0], actual[1], actual[2]);
      for (int i = 0; i < 3; ++i) {
        if (std::memcmp(
                expected[i].data(), actual[i].data(), count * sizeof(float))
            != 0) {
          throw std::runtime_error(
              "parallel Ooura field differs from serial library");
        }
      }
    }
  }
}
