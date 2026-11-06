#pragma once
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace flasheda::density {
struct Cell {
  int lx, ly, ux, uy, kind; // 0 ignored, 1 standard, 2 macro, 3 filler
  float scale, gx, gy;
};
struct Bin {
  int lx, ly, ux, uy;
  int64_t fixed, fixed_unscaled, inst, filler;
  float target, density, phi, ex, ey;
};
struct Snapshot {
  int nx = 0, ny = 0;
  double sx = 0, sy = 0;
  std::vector<Cell> cells;
  std::vector<Bin> bins;
};
inline void write(const Snapshot &s, const std::string &path) {
  std::ofstream f(path);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  f << std::setprecision(17) << "FLASHEDA_DENSITY_V1 " << s.nx << ' ' << s.ny
    << ' ' << s.sx << ' ' << s.sy << ' ' << s.cells.size() << '\n';
  for (const auto &c : s.cells)
    f << c.lx << ' ' << c.ly << ' ' << c.ux << ' ' << c.uy << ' ' << c.kind
      << ' ' << c.scale << ' ' << c.gx << ' ' << c.gy << '\n';
  for (const auto &b : s.bins)
    f << b.lx << ' ' << b.ly << ' ' << b.ux << ' ' << b.uy << ' ' << b.fixed
      << ' ' << b.fixed_unscaled << ' ' << b.inst << ' ' << b.filler << ' '
      << b.target << ' ' << b.density << ' ' << b.phi << ' ' << b.ex << ' '
      << b.ey << '\n';
}
inline Snapshot read(const std::string &path) {
  std::ifstream f(path);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  Snapshot s;
  std::string magic;
  size_t count;
  f >> magic >> s.nx >> s.ny >> s.sx >> s.sy >> count;
  if (magic != "FLASHEDA_DENSITY_V1" || s.nx < 2 || s.ny < 2 || s.nx > 4096 ||
      s.ny > 4096 || count > 10000000 || s.sx <= 0 || s.sy <= 0)
    throw std::runtime_error("Invalid density snapshot header");
  s.cells.resize(count);
  s.bins.resize(size_t(s.nx) * s.ny);
  for (auto &c : s.cells)
    f >> c.lx >> c.ly >> c.ux >> c.uy >> c.kind >> c.scale >> c.gx >> c.gy;
  for (auto &b : s.bins)
    f >> b.lx >> b.ly >> b.ux >> b.uy >> b.fixed >> b.fixed_unscaled >>
        b.inst >> b.filler >> b.target >> b.density >> b.phi >> b.ex >> b.ey;
  return s;
}
} // namespace flasheda::density
