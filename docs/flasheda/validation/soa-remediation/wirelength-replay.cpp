#include <bit>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include "wirelengthSoa.h"
#include "previousWirelength.h"

using Clock = std::chrono::steady_clock;
template <size_t... F, typename Store>
void compare(const Store& a, const Store& b, std::index_sequence<F...>)
{
  const auto same = []<typename T>(const std::vector<T>& x,
                                 const std::vector<T>& y) {
    if (x.size() != y.size()) return false;
    for (size_t i = 0; i < x.size(); ++i) {
      if constexpr (std::is_same_v<T, float>) {
        if (std::bit_cast<uint32_t>(x[i]) != std::bit_cast<uint32_t>(y[i]))
          return false;
      } else if (x[i] != y[i]) return false;
    }
    return true;
  };
  if (!(same(a.template column<F>(), b.template column<F>()) && ...))
    throw std::runtime_error("wirelength columns differ");
}

int main()
{
  using namespace gpl;
  GNetColumns::Store nets;
  GPinColumns::Store pins;
  std::vector<size_t> offsets{0}, entries;
  std::mt19937 rng(42);
  constexpr size_t count = 10000;
  for (size_t n = 0; n < count; ++n) {
    nets.allocate();
    const size_t degree = n % 251;
    for (size_t k = 0; k < degree; ++k) {
      const size_t p = pins.allocate();
      pins.at<GPinColumns::cx>(p) = int(rng() % 2000001) - 1000000;
      pins.at<GPinColumns::cy>(p) = int(rng() % 2000001) - 1000000;
      entries.push_back(p);
    }
    offsets.push_back(entries.size());
  }
  auto old_nets = nets;
  auto old_pins = pins;
  for (int sample = 0; sample < 11; ++sample) {
    for (int order = 0; order < 2; ++order) {
      const bool previous = (order + sample) % 2 == 0;
      const float coefficient = sample % 3 == 0 ? 0.0f : 0.00001f;
      const float minimum = sample % 5 == 0 ? 0.0f : -10.0f;
      const auto start = Clock::now();
      for (size_t n = 0; n < count; ++n) {
        const auto rows = std::span<const size_t>(entries).subspan(
            offsets[n], offsets[n + 1] - offsets[n]);
        if (previous)
          updateNetForceColumnsPrevious(old_nets, old_pins, n, rows,
                                        coefficient, coefficient * 2, minimum);
        else
          updateNetForceColumns(nets, pins, n, rows,
                                coefficient, coefficient * 2, minimum);
      }
      const double ms = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
      std::cout << "{\"variant\":\"" << (previous ? "previous" : "registers")
                << "\",\"sample\":" << sample << ",\"ms\":" << ms << "}\n";
    }
    compare(nets, old_nets, std::make_index_sequence<15>{});
    compare(pins, old_pins, std::make_index_sequence<12>{});
  }
}
