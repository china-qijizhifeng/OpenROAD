// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrPotential.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <tuple>
#include <utility>
#include <vector>

#include "br/BrCore.h"

namespace drt::br {

namespace {

Cost addCost(const Cost a, const Cost b)
{
  if (a >= kInfCost || b >= kInfCost) {
    return kInfCost;
  }
  const Cost sum = a + b;
  return sum >= kInfCost ? kInfCost : sum;
}

// rate * length for length >= 0; a forbidden direction costs nothing over a
// zero length.
Cost scaleCost(const Cost rate, const Cost length)
{
  if (length == 0) {
    return 0;
  }
  if (rate >= kInfCost) {
    return kInfCost;
  }
  const Cost product = rate * length;
  return product >= kInfCost ? kInfCost : product;
}

Cost clampCost(const Cost c)
{
  return c >= kInfCost ? kInfCost : c;
}

// floor(a / b) for b > 0.
Cost floorDiv(const Cost a, const Cost b)
{
  const Cost q = a / b;
  return (a % b != 0 && a < 0) ? q - 1 : q;
}

Cost axisGap(const Coord c, const Coord lo, const Coord hi)
{
  if (c < lo) {
    return static_cast<Cost>(lo) - c;
  }
  if (c > hi) {
    return static_cast<Cost>(c) - hi;
  }
  return 0;
}

// Prefix sums of the finite via costs and prefix counts of forbidden vias:
// the via stack between z0 and z1 is forbidden iff a forbidden via lies in
// between.
void buildViaPrefix(const std::vector<Cost>& via_up,
                    const int num_layers,
                    std::vector<Cost>& prefix,
                    std::vector<int>& forbidden)
{
  prefix.assign(std::max(num_layers, 1), 0);
  forbidden.assign(std::max(num_layers, 1), 0);
  for (int z = 1; z < num_layers; z++) {
    const Cost v
        = z - 1 < static_cast<int>(via_up.size()) ? via_up[z - 1] : kInfCost;
    const bool blocked = v >= kInfCost;
    prefix[z] = prefix[z - 1] + (blocked ? 0 : v);
    forbidden[z] = forbidden[z - 1] + (blocked ? 1 : 0);
  }
}

Cost viaStack(const std::vector<Cost>& prefix,
              const std::vector<int>& forbidden,
              int z0,
              int z1)
{
  if (z0 > z1) {
    std::swap(z0, z1);
  }
  if (forbidden[z1] != forbidden[z0]) {
    return kInfCost;
  }
  return prefix[z1] - prefix[z0];
}

}  // namespace

// ---------------------------------------------------------------------------
// L1Potential

void L1Potential::prepare(const LayerCostRates& rates,
                          const std::vector<TargetRectangle>& targets)
{
  counters_.prepares++;
  const int num_layers = rates.numLayers();
  min_horizontal_ = kInfCost;
  min_vertical_ = kInfCost;
  for (int z = 0; z < num_layers; z++) {
    min_horizontal_ = std::min(min_horizontal_, rates.horizontal[z]);
    min_vertical_ = std::min(min_vertical_, rates.vertical[z]);
  }
  buildViaPrefix(rates.via_up, num_layers, via_prefix_, via_forbidden_prefix_);
  targets_.assign(targets.begin(), targets.end());
}

Cost L1Potential::viaLowerBound(const int z0, const int z1) const
{
  return viaStack(via_prefix_, via_forbidden_prefix_, z0, z1);
}

Cost L1Potential::query(const int z, const Coord x, const Coord y) const
{
  counters_.queries++;
  Cost best = kInfCost;
  for (const TargetRectangle& t : targets_) {
    Cost d = addCost(scaleCost(min_horizontal_, axisGap(x, t.xlo, t.xhi)),
                     scaleCost(min_vertical_, axisGap(y, t.ylo, t.yhi)));
    d = addCost(d, viaLowerBound(z, t.z));
    best = std::min(best, addCost(d, t.terminal_cost));
  }
  return best;
}

// ---------------------------------------------------------------------------
// SimplePotential

void SimplePotential::prepare(const LayerCostRates& rates,
                              const std::vector<TargetRectangle>& targets)
{
  counters_.prepares++;
  // The table depends on the rates only.
  const bool same_rates = table_valid_ && num_layers_ == rates.numLayers()
                          && rates_.horizontal == rates.horizontal
                          && rates_.vertical == rates.vertical
                          && std::equal(rates.via_up.begin(),
                                        rates.via_up.end(),
                                        rates_.via_up.begin(),
                                        rates_.via_up.end());
  num_layers_ = rates.numLayers();
  assert(num_layers_ <= kMaxPotentialLayers);
  rates_.horizontal.assign(rates.horizontal.begin(), rates.horizontal.end());
  rates_.vertical.assign(rates.vertical.begin(), rates.vertical.end());
  rates_.via_up.assign(rates.via_up.begin(), rates.via_up.end());
  rates_.via_up.resize(std::max(num_layers_ - 1, 0), kInfCost);
  buildViaPrefix(
      rates_.via_up, num_layers_, via_prefix_, via_forbidden_prefix_);
  min_horizontal_ = kInfCost;
  min_vertical_ = kInfCost;
  for (int z = 0; z < num_layers_; z++) {
    min_horizontal_ = std::min(min_horizontal_, rates_.horizontal[z]);
    min_vertical_ = std::min(min_vertical_, rates_.vertical[z]);
  }
  targets_.assign(targets.begin(), targets.end());
  if (!same_rates) {
    buildTable();
  }
}

Cost SimplePotential::viaCost(const int z0, const int z1) const
{
  return viaStack(via_prefix_, via_forbidden_prefix_, z0, z1);
}

void SimplePotential::buildTable()
{
  const int n = num_layers_;
  functions_.clear();
  first_.assign(static_cast<size_t>(n) * n + 1, 0);
  std::vector<LinearCost> all;
  for (int tz = 0; tz < n; tz++) {
    for (int z = 0; z < n; z++) {
      all.clear();
      for (int order = 0; order < 2; order++) {
        const std::vector<Cost>& first_rate
            = order == 0 ? rates_.horizontal : rates_.vertical;
        const std::vector<Cost>& second_rate
            = order == 0 ? rates_.vertical : rates_.horizontal;
        for (int a = 0; a < n; a++) {
          for (int b = 0; b < n; b++) {
            const Cost c = addCost(addCost(viaCost(tz, a), viaCost(a, b)),
                                   viaCost(b, z));
            if (c >= kInfCost) {
              continue;
            }
            const Cost h = order == 0 ? first_rate[a] : second_rate[b];
            const Cost v = order == 0 ? second_rate[b] : first_rate[a];
            all.push_back({c, h, v});
          }
        }
      }
      // Keep the functions not bounded from below by another one (ties:
      // the first in a fixed order).
      std::ranges::sort(all, [](const LinearCost& p, const LinearCost& q) {
        return std::tie(p.constant, p.horizontal, p.vertical)
               < std::tie(q.constant, q.horizontal, q.vertical);
      });
      const size_t begin = functions_.size();
      for (const LinearCost& f : all) {
        bool dominated = false;
        for (size_t k = begin; k < functions_.size() && !dominated; k++) {
          const LinearCost& g = functions_[k];
          dominated = g.constant <= f.constant && g.horizontal <= f.horizontal
                      && g.vertical <= f.vertical;
        }
        if (!dominated) {
          functions_.push_back(f);
        }
      }
      first_[static_cast<size_t>(tz) * n + z + 1]
          = static_cast<int>(functions_.size());
    }
  }
  table_valid_ = true;
}

Cost SimplePotential::pointDistance(const int target_z,
                                    const int z,
                                    const Cost dx,
                                    const Cost dy) const
{
  const size_t k = static_cast<size_t>(target_z) * num_layers_ + z;
  Cost best = kInfCost;
  for (int f = first_[k]; f < first_[k + 1]; f++) {
    const LinearCost& c = functions_[f];
    best = std::min(best,
                    addCost(addCost(c.constant, scaleCost(c.horizontal, dx)),
                            scaleCost(c.vertical, dy)));
  }
  return best;
}

// Theorem 3. Paths from r = (target point, target_z) to s = (query, z) of
// the form vias, horizontal run on layer a, vias, vertical run on layer b,
// vias (and the same with the two runs swapped) contain a shortest path
// (Lemma 1). With first[z] = via(r_z, z) + first_rate[z] * first_length,
// d[z] = min_a first[a] + via(a, z) is the distance from r to the corner
// vertex on layer z; two sweeps over the via chain compute it in O(l), and
// the answer is min_b d[b] + second_rate[b] * second_length + via(b, s_z).
Cost SimplePotential::pointDistanceSweep(const int target_z,
                                         const int z,
                                         const Cost dx,
                                         const Cost dy) const
{
  if (dx == 0 && dy == 0) {
    return viaCost(target_z, z);
  }
  const int num_layers = num_layers_;
  const std::vector<Cost>& via = rates_.via_up;
  std::array<Cost, kMaxPotentialLayers> corner;
  Cost best = kInfCost;
  for (int order = 0; order < 2; order++) {
    const std::vector<Cost>& first_rate
        = order == 0 ? rates_.horizontal : rates_.vertical;
    const std::vector<Cost>& second_rate
        = order == 0 ? rates_.vertical : rates_.horizontal;
    const Cost first_length = order == 0 ? dx : dy;
    const Cost second_length = order == 0 ? dy : dx;
    // Downward sweep: corner[a] = min over a' >= a of first[a'] + via(a', a).
    Cost carry = kInfCost;
    for (int a = num_layers - 1; a >= 0; a--) {
      if (a + 1 < num_layers) {
        carry = addCost(carry, via[a]);
      }
      const Cost first = addCost(viaCost(target_z, a),
                                 scaleCost(first_rate[a], first_length));
      carry = std::min(carry, first);
      corner[a] = carry;
    }
    // Upward sweep completes the minimum over all a and evaluates the
    // second run on each layer b.
    carry = kInfCost;
    for (int b = 0; b < num_layers; b++) {
      if (b > 0) {
        carry = addCost(carry, via[b - 1]);
      }
      carry = std::min(carry, corner[b]);
      if (carry >= best) {
        continue;
      }
      const Cost total
          = addCost(addCost(carry, scaleCost(second_rate[b], second_length)),
                    viaCost(b, z));
      best = std::min(best, total);
    }
  }
  return best;
}

Cost SimplePotential::query(const int z, const Coord x, const Coord y) const
{
  counters_.queries++;
  // The L1 bound (minimum rates times the gaps plus the via stack, which
  // every path between the two layers uses) orders the targets; the exact
  // Theorem 3 evaluation stops once the bound reaches the best distance.
  constexpr int kOrderedTargets = 16;
  const int num_targets = static_cast<int>(targets_.size());
  std::array<std::pair<Cost, int>, kOrderedTargets> order;
  const int ordered = std::min(num_targets, kOrderedTargets);
  auto bound = [&](const TargetRectangle& t) {
    const Cost gaps
        = addCost(scaleCost(min_horizontal_, axisGap(x, t.xlo, t.xhi)),
                  scaleCost(min_vertical_, axisGap(y, t.ylo, t.yhi)));
    return addCost(addCost(gaps, viaCost(t.z, z)), t.terminal_cost);
  };
  for (int k = 0; k < ordered; k++) {
    order[k] = {bound(targets_[k]), k};
    for (int m = k; m > 0 && order[m].first < order[m - 1].first; m--) {
      std::swap(order[m], order[m - 1]);
    }
  }
  Cost best = kInfCost;
  auto evaluate = [&](const TargetRectangle& t) {
    const Cost d = pointDistance(
        t.z, z, axisGap(x, t.xlo, t.xhi), axisGap(y, t.ylo, t.yhi));
    best = std::min(best, addCost(d, t.terminal_cost));
  };
  for (int k = 0; k < ordered && order[k].first < best; k++) {
    evaluate(targets_[order[k].second]);
  }
  for (int k = ordered; k < num_targets; k++) {
    if (bound(targets_[k]) < best) {
      evaluate(targets_[k]);
    }
  }
  return best;
}

// ---------------------------------------------------------------------------
// TileCostModel

void TileCostModel::resize(const int layers, const Cost rate)
{
  num_layers = layers;
  const size_t n = static_cast<size_t>(layers) * numColumns() * numRows();
  horizontal.assign(n, rate);
  vertical.assign(n, rate);
  via_up.assign(n, rate);
}

void TileCostModel::assignSimple(const LayerCostRates& rates,
                                 const Coord xlo,
                                 const Coord ylo,
                                 const Coord xhi,
                                 const Coord yhi)
{
  x_cuts.clear();
  y_cuts.clear();
  domain_xlo = xlo;
  domain_ylo = ylo;
  domain_xhi = xhi;
  domain_yhi = yhi;
  resize(rates.numLayers(), kInfCost);
  for (int z = 0; z < num_layers; z++) {
    horizontal[index(z, 0, 0)] = rates.horizontal[z];
    vertical[index(z, 0, 0)] = rates.vertical[z];
    if (z + 1 < num_layers) {
      via_up[index(z, 0, 0)] = rates.via_up[z];
    }
  }
}

LayerCostRates minimumLayerRates(const TileCostModel& model)
{
  LayerCostRates rates;
  const int num_layers = model.num_layers;
  rates.horizontal.assign(num_layers, kInfCost);
  rates.vertical.assign(num_layers, kInfCost);
  rates.via_up.assign(std::max(num_layers - 1, 0), kInfCost);
  for (int z = 0; z < num_layers; z++) {
    for (int i = 0; i < model.numColumns(); i++) {
      for (int j = 0; j < model.numRows(); j++) {
        const size_t k = model.index(z, i, j);
        rates.horizontal[z]
            = std::min(rates.horizontal[z], model.horizontal[k]);
        rates.vertical[z] = std::min(rates.vertical[z], model.vertical[k]);
        if (z + 1 < num_layers) {
          rates.via_up[z] = std::min(rates.via_up[z], model.via_up[k]);
        }
      }
    }
  }
  return rates;
}

// ---------------------------------------------------------------------------
// TileCostModelBuilder

void TileCostModelBuilder::collectLines(std::vector<Line>& lines,
                                        std::vector<Coord>& cuts)
{
  std::sort(lines.begin(), lines.end(), [](const Line& a, const Line& b) {
    return a.position < b.position;
  });
  cuts.clear();
  for (size_t k = 0; k < lines.size();) {
    bool doubled = false;
    size_t e = k;
    for (; e < lines.size() && lines[e].position == lines[k].position; e++) {
      doubled = doubled || lines[e].doubled;
    }
    cuts.push_back(lines[k].position);
    if (doubled) {
      cuts.push_back(lines[k].position);
    }
    k = e;
  }
}

void TileCostModelBuilder::build(const Region& region,
                                 const GCellGrid& grid,
                                 const LayerCostRates& base,
                                 const std::vector<DiscountedWire>& wires,
                                 TileCostModel& model)
{
  const int num_layers = base.numLayers();
  model.x_cuts.clear();
  model.y_cuts.clear();
  if (region.bx_hi < 0 || region.by_hi < 0) {
    model.domain_xlo = grid.die_xlo;
    model.domain_ylo = grid.die_ylo;
    model.domain_xhi = grid.die_xhi;
    model.domain_yhi = grid.die_yhi;
    model.resize(num_layers, kInfCost);
    return;
  }
  const int tiles_x = region.tiles_x;
  const int tiles_y = region.tiles_y;
  const Coord dom_xlo = region.bx_lo == 0
                            ? std::min(grid.die_xlo, grid.origin_x)
                            : grid.origin_x + region.bx_lo * grid.pitch_x;
  const Coord dom_xhi
      = region.bx_hi == tiles_x - 1
            ? std::max(grid.die_xhi, grid.origin_x + tiles_x * grid.pitch_x)
            : grid.origin_x + (region.bx_hi + 1) * grid.pitch_x;
  const Coord dom_ylo = region.by_lo == 0
                            ? std::min(grid.die_ylo, grid.origin_y)
                            : grid.origin_y + region.by_lo * grid.pitch_y;
  const Coord dom_yhi
      = region.by_hi == tiles_y - 1
            ? std::max(grid.die_yhi, grid.origin_y + tiles_y * grid.pitch_y)
            : grid.origin_y + (region.by_hi + 1) * grid.pitch_y;
  model.domain_xlo = dom_xlo;
  model.domain_ylo = dom_ylo;
  model.domain_xhi = dom_xhi;
  model.domain_yhi = dom_yhi;

  x_lines_.clear();
  y_lines_.clear();
  for (int tx = region.bx_lo; tx <= region.bx_hi + 1; tx++) {
    if (tx > 0 && tx < tiles_x) {
      x_lines_.push_back({grid.origin_x + tx * grid.pitch_x, false});
    }
  }
  for (int ty = region.by_lo; ty <= region.by_hi + 1; ty++) {
    if (ty > 0 && ty < tiles_y) {
      y_lines_.push_back({grid.origin_y + ty * grid.pitch_y, false});
    }
  }
  horizontal_wires_.clear();
  vertical_wires_.clear();
  for (const DiscountedWire& w : wires) {
    if (w.z < 0 || w.z >= num_layers || w.denominator <= 0 || w.numerator < 0) {
      continue;
    }
    if (w.y1 == w.y2 && w.x1 != w.x2) {
      DiscountedWire c = w;
      c.x1 = std::max(std::min(w.x1, w.x2), dom_xlo);
      c.x2 = std::min(std::max(w.x1, w.x2), dom_xhi);
      if (c.x1 >= c.x2 || w.y1 < dom_ylo || w.y1 > dom_yhi) {
        continue;
      }
      horizontal_wires_.push_back(c);
      x_lines_.push_back({c.x1, false});
      x_lines_.push_back({c.x2, false});
      y_lines_.push_back({c.y1, true});
    } else if (w.x1 == w.x2 && w.y1 != w.y2) {
      DiscountedWire c = w;
      c.y1 = std::max(std::min(w.y1, w.y2), dom_ylo);
      c.y2 = std::min(std::max(w.y1, w.y2), dom_yhi);
      if (c.y1 >= c.y2 || w.x1 < dom_xlo || w.x1 > dom_xhi) {
        continue;
      }
      vertical_wires_.push_back(c);
      y_lines_.push_back({c.y1, false});
      y_lines_.push_back({c.y2, false});
      x_lines_.push_back({c.x1, true});
    }
  }
  collectLines(x_lines_, model.x_cuts);
  collectLines(y_lines_, model.y_cuts);
  model.resize(num_layers, kInfCost);

  const int p = static_cast<int>(model.x_cuts.size());
  const int q = static_cast<int>(model.y_cuts.size());
  const auto& xc = model.x_cuts;
  const auto& yc = model.y_cuts;
  auto column_degenerate
      = [&](int i) { return i >= 1 && i < p && xc[i - 1] == xc[i]; };
  auto row_degenerate
      = [&](int j) { return j >= 1 && j < q && yc[j - 1] == yc[j]; };
  auto gcell_x = [&](Coord x) {
    return std::clamp(
        static_cast<int>(floorDiv(Cost{x} - grid.origin_x, grid.pitch_x)),
        0,
        tiles_x - 1);
  };
  auto gcell_y = [&](Coord y) {
    return std::clamp(
        static_cast<int>(floorDiv(Cost{y} - grid.origin_y, grid.pitch_y)),
        0,
        tiles_y - 1);
  };
  // Regular tiles: each lies in one GCell (all GCell boundaries are cuts).
  for (int i = 0; i <= p; i++) {
    if (column_degenerate(i)) {
      continue;
    }
    const int tx = i == 0 ? gcell_x(p > 0 ? xc[0] - 1 : grid.origin_x)
                          : gcell_x(xc[i - 1]);
    for (int j = 0; j <= q; j++) {
      if (row_degenerate(j)) {
        continue;
      }
      const int ty = j == 0 ? gcell_y(q > 0 ? yc[0] - 1 : grid.origin_y)
                            : gcell_y(yc[j - 1]);
      for (int z = 0; z < num_layers; z++) {
        if (!region.allowedTile(z, tx, ty)) {
          continue;
        }
        const size_t k = model.index(z, i, j);
        model.horizontal[k] = base.horizontal[z];
        model.vertical[k] = base.vertical[z];
        if (z + 1 < num_layers && region.allowedTile(z + 1, tx, ty)) {
          model.via_up[k] = base.via_up[z];
        }
      }
    }
  }
  auto discounted = [](Cost rate, const DiscountedWire& w) {
    if (rate >= kInfCost) {
      return kInfCost;
    }
    return rate * w.numerator / w.denominator;
  };
  // Degenerate rows: horizontal reservations on the line.
  for (int j = 1; j < q; j++) {
    if (!row_degenerate(j)) {
      continue;
    }
    const Coord line = yc[j - 1];
    for (const DiscountedWire& w : horizontal_wires_) {
      if (w.y1 != line) {
        continue;
      }
      for (int i = 1; i < p; i++) {
        if (column_degenerate(i) || xc[i - 1] < w.x1 || xc[i] > w.x2) {
          continue;
        }
        const Cost neighbor
            = std::min(model.horizontal[model.index(w.z, i, j - 1)],
                       model.horizontal[model.index(w.z, i, j + 1)]);
        Cost& rate = model.horizontal[model.index(w.z, i, j)];
        rate = std::min(rate, discounted(neighbor, w));
      }
    }
  }
  // Degenerate columns: vertical reservations on the line.
  for (int i = 1; i < p; i++) {
    if (!column_degenerate(i)) {
      continue;
    }
    const Coord line = xc[i - 1];
    for (const DiscountedWire& w : vertical_wires_) {
      if (w.x1 != line) {
        continue;
      }
      for (int j = 1; j < q; j++) {
        if (row_degenerate(j) || yc[j - 1] < w.y1 || yc[j] > w.y2) {
          continue;
        }
        const Cost neighbor
            = std::min(model.vertical[model.index(w.z, i - 1, j)],
                       model.vertical[model.index(w.z, i + 1, j)]);
        Cost& rate = model.vertical[model.index(w.z, i, j)];
        rate = std::min(rate, discounted(neighbor, w));
      }
    }
  }
  mergeEqualColumns(model);
  mergeEqualRows(model);
}

// Removes a single (not doubled) cut between two columns with identical
// rates: every edge on the cut line had the minimum of two equal rates.
void TileCostModelBuilder::mergeEqualColumns(TileCostModel& model)
{
  const int p = static_cast<int>(model.x_cuts.size());
  const int rows = model.numRows();
  const int num_layers = model.num_layers;
  auto same = [&](int a, int b) {
    for (int z = 0; z < num_layers; z++) {
      for (int j = 0; j < rows; j++) {
        const size_t ka = model.index(z, a, j);
        const size_t kb = model.index(z, b, j);
        if (model.horizontal[ka] != model.horizontal[kb]
            || model.vertical[ka] != model.vertical[kb]
            || model.via_up[ka] != model.via_up[kb]) {
          return false;
        }
      }
    }
    return true;
  };
  keep_.assign(p, 1);
  bool any = false;
  for (int c = 0; c < p; c++) {
    const bool doubled
        = (c > 0 && model.x_cuts[c - 1] == model.x_cuts[c])
          || (c + 1 < p && model.x_cuts[c + 1] == model.x_cuts[c]);
    if (!doubled && same(c, c + 1)) {
      keep_[c] = 0;
      any = true;
    }
  }
  if (!any) {
    return;
  }
  scratch_cuts_.clear();
  std::vector<int>& columns = kept_indices_;  // kept old column indices
  columns.assign(1, 0);
  for (int c = 0; c < p; c++) {
    if (keep_[c]) {
      scratch_cuts_.push_back(model.x_cuts[c]);
      columns.push_back(c + 1);
    }
  }
  const int new_columns = static_cast<int>(columns.size());
  for (std::vector<Cost>* rates :
       {&model.horizontal, &model.vertical, &model.via_up}) {
    scratch_.resize(static_cast<size_t>(num_layers) * new_columns * rows);
    for (int z = 0; z < num_layers; z++) {
      for (int i = 0; i < new_columns; i++) {
        for (int j = 0; j < rows; j++) {
          scratch_[(static_cast<size_t>(z) * new_columns + i) * rows + j]
              = (*rates)[model.index(z, columns[i], j)];
        }
      }
    }
    rates->swap(scratch_);
  }
  model.x_cuts.swap(scratch_cuts_);
}

void TileCostModelBuilder::mergeEqualRows(TileCostModel& model)
{
  const int q = static_cast<int>(model.y_cuts.size());
  const int columns = model.numColumns();
  const int num_layers = model.num_layers;
  auto same = [&](int a, int b) {
    for (int z = 0; z < num_layers; z++) {
      for (int i = 0; i < columns; i++) {
        const size_t ka = model.index(z, i, a);
        const size_t kb = model.index(z, i, b);
        if (model.horizontal[ka] != model.horizontal[kb]
            || model.vertical[ka] != model.vertical[kb]
            || model.via_up[ka] != model.via_up[kb]) {
          return false;
        }
      }
    }
    return true;
  };
  keep_.assign(q, 1);
  bool any = false;
  for (int c = 0; c < q; c++) {
    const bool doubled
        = (c > 0 && model.y_cuts[c - 1] == model.y_cuts[c])
          || (c + 1 < q && model.y_cuts[c + 1] == model.y_cuts[c]);
    if (!doubled && same(c, c + 1)) {
      keep_[c] = 0;
      any = true;
    }
  }
  if (!any) {
    return;
  }
  scratch_cuts_.clear();
  std::vector<int>& rows = kept_indices_;
  rows.assign(1, 0);
  for (int c = 0; c < q; c++) {
    if (keep_[c]) {
      scratch_cuts_.push_back(model.y_cuts[c]);
      rows.push_back(c + 1);
    }
  }
  const int new_rows = static_cast<int>(rows.size());
  for (std::vector<Cost>* rates :
       {&model.horizontal, &model.vertical, &model.via_up}) {
    scratch_.resize(static_cast<size_t>(num_layers) * columns * new_rows);
    for (int z = 0; z < num_layers; z++) {
      for (int i = 0; i < columns; i++) {
        for (int j = 0; j < new_rows; j++) {
          scratch_[(static_cast<size_t>(z) * columns + i) * new_rows + j]
              = (*rates)[model.index(z, i, rows[j])];
        }
      }
    }
    rates->swap(scratch_);
  }
  model.y_cuts.swap(scratch_cuts_);
}

// ---------------------------------------------------------------------------
// GeneralPotential

void GeneralPotential::Axis::build(const std::vector<Coord>& cuts,
                                   const std::vector<Coord>& extra,
                                   const Coord lo,
                                   const Coord hi)
{
  positions.assign(cuts.begin(), cuts.end());
  positions.insert(positions.end(), extra.begin(), extra.end());
  std::sort(positions.begin(), positions.end());
  positions.erase(std::unique(positions.begin(), positions.end()),
                  positions.end());
  domain_lo = lo;
  domain_hi = hi;
  if (!positions.empty()) {
    domain_lo = std::min(domain_lo, positions.front());
    domain_hi = std::max(domain_hi, positions.back());
  }
  const int num_positions = numPositions();
  first_tile.resize(numSlots());
  last_tile.resize(numSlots());
  for (int k = 0; k <= num_positions; k++) {
    // Open interval below positions[k]: the tile column after the last cut
    // <= positions[k - 1].
    const int tile
        = k == 0
              ? 0
              : static_cast<int>(
                    std::upper_bound(cuts.begin(), cuts.end(), positions[k - 1])
                    - cuts.begin());
    first_tile[2 * k] = tile;
    last_tile[2 * k] = tile;
    if (k < num_positions) {
      const Coord c = positions[k];
      first_tile[2 * k + 1] = static_cast<int>(
          std::lower_bound(cuts.begin(), cuts.end(), c) - cuts.begin());
      last_tile[2 * k + 1] = static_cast<int>(
          std::upper_bound(cuts.begin(), cuts.end(), c) - cuts.begin());
    }
  }
}

int GeneralPotential::Axis::locate(const Coord c) const
{
  const int k
      = static_cast<int>(std::lower_bound(positions.begin(), positions.end(), c)
                         - positions.begin());
  return (k < numPositions() && positions[k] == c) ? 2 * k + 1 : 2 * k;
}

Coord GeneralPotential::Axis::slotLo(const int slot) const
{
  const int k = slot / 2;
  if (slot & 1) {
    return positions[k];
  }
  return k == 0 ? domain_lo : positions[k - 1] + 1;
}

Coord GeneralPotential::Axis::slotHi(const int slot) const
{
  const int k = slot / 2;
  if (slot & 1) {
    return positions[k];
  }
  return k == numPositions() ? domain_hi : positions[k] - 1;
}

Cost GeneralPotential::rate(const RateKind kind,
                            const int x_slot,
                            const int y_slot,
                            const int z) const
{
  const std::vector<Cost>& rates = rates_[kind];
  Cost best = kInfCost;
  const int i0 = x_axis_.first_tile[x_slot];
  const int i1 = x_axis_.last_tile[x_slot];
  const int j0 = y_axis_.first_tile[y_slot];
  const int j1 = y_axis_.last_tile[y_slot];
  for (int i = i0; i <= i1; i++) {
    const size_t row_base
        = (static_cast<size_t>(z) * model_columns_ + i) * model_rows_;
    for (int j = j0; j <= j1; j++) {
      best = std::min(best, rates[row_base + j]);
    }
  }
  return best;
}

Cost GeneralPotential::runRate(const LineRun& run,
                               const RateKind kind,
                               const int along_slot,
                               const int line_slot,
                               const int z) const
{
  if (run.horizontal_lines) {
    return rate(kind, along_slot, line_slot, z);
  }
  const RateKind swapped = kind == kHorizontalRate ? kVerticalRate
                           : kind == kVerticalRate ? kHorizontalRate
                                                   : kViaRate;
  return rate(swapped, line_slot, along_slot, z);
}

void GeneralPotential::prepare(const TileCostModel& model,
                               const std::vector<TargetRectangle>& targets)
{
  counters_.prepares++;
  num_layers_ = model.num_layers;
  assert(num_layers_ <= kMaxPotentialLayers);
  model_columns_ = model.numColumns();
  model_rows_ = model.numRows();
  rates_[kHorizontalRate].assign(model.horizontal.begin(),
                                 model.horizontal.end());
  rates_[kVerticalRate].assign(model.vertical.begin(), model.vertical.end());
  rates_[kViaRate].assign(model.via_up.begin(), model.via_up.end());
  targets_.clear();
  for (const TargetRectangle& t : targets) {
    if (t.z >= 0 && t.z < num_layers_ && t.xlo <= t.xhi && t.ylo <= t.yhi
        && t.terminal_cost < kInfCost) {
      targets_.push_back(t);
    }
  }
  has_targets_ = !targets_.empty();
  functions_.clear();
  heap_.clear();
  tile_functions_.clear();
  if (++tile_epoch_ == 0) {
    for (TileListEntry& e : tile_lists_) {
      e.stamp = 0;
    }
    tile_epoch_ = 1;
  }
  if (!has_targets_) {
    return;
  }

  // Refine the grid with respect to the targets (Corollary 15).
  Coord xlo = model.domain_xlo;
  Coord xhi = model.domain_xhi;
  Coord ylo = model.domain_ylo;
  Coord yhi = model.domain_yhi;
  scratch_positions_.clear();
  for (const TargetRectangle& t : targets_) {
    scratch_positions_.push_back(t.xlo);
    scratch_positions_.push_back(t.xhi);
    xlo = std::min(xlo, t.xlo);
    xhi = std::max(xhi, t.xhi);
  }
  x_axis_.build(model.x_cuts, scratch_positions_, xlo, xhi);
  scratch_positions_.clear();
  for (const TargetRectangle& t : targets_) {
    scratch_positions_.push_back(t.ylo);
    scratch_positions_.push_back(t.yhi);
    ylo = std::min(ylo, t.ylo);
    yhi = std::max(yhi, t.yhi);
  }
  y_axis_.build(model.y_cuts, scratch_positions_, ylo, yhi);

  tile_lists_.resize(static_cast<size_t>(num_layers_)
                     * (x_axis_.numPositions() + 1)
                     * (y_axis_.numPositions() + 1));

  for (LineRun* run : {&horizontal_run_, &vertical_run_}) {
    run->horizontal_lines = run == &horizontal_run_;
    const Axis& along = alongAxis(*run);
    const Axis& lines = lineAxis(*run);
    run->num_lines = lines.numPositions();
    run->num_slots = along.numSlots();
    run->head.assign(
        static_cast<size_t>(num_layers_) * run->num_lines * run->num_slots, -1);
    counters_.boundary_elements += static_cast<int64_t>(run->head.size());
    // Initialization: every element inside a target gets the constant
    // terminal cost; elements touching a target get their functions by
    // propagation from the target's grid points.
    for (const TargetRectangle& t : targets_) {
      const Coord a_lo = run->horizontal_lines ? t.xlo : t.ylo;
      const Coord a_hi = run->horizontal_lines ? t.xhi : t.yhi;
      const Coord l_lo = run->horizontal_lines ? t.ylo : t.xlo;
      const Coord l_hi = run->horizontal_lines ? t.yhi : t.xhi;
      const int s0 = along.locate(a_lo);
      const int s1 = along.locate(a_hi);
      const int m0 = lines.locate(l_lo) / 2;
      const int m1 = lines.locate(l_hi) / 2;
      for (int m = m0; m <= m1; m++) {
        for (int s = s0; s <= s1; s++) {
          const Coord lo = along.slotLo(s);
          const Coord hi = along.slotHi(s);
          if (lo <= hi) {
            offer(*run, run->element(t.z, m, s), lo, hi, t.terminal_cost, 0);
          }
        }
      }
    }
    runPropagation(*run);
  }
}

// Lemma 9 over integer points: the new function g takes over the integer
// interval [u, w] where it is strictly below the envelope (an interval by
// concavity); functions owning only points of [u, w] are removed, the two
// neighbours are trimmed (their keys can only increase).
void GeneralPotential::offer(LineRun& run,
                             const size_t element,
                             const Coord lo,
                             const Coord hi,
                             const Cost offset,
                             const Cost slope)
{
  counters_.functions_offered++;
  int32_t head = run.head[element];
  Cost u = Cost{hi} + 1;
  Cost w = Cost{lo} - 1;
  if (head < 0) {
    u = lo;
    w = hi;
  } else {
    for (int32_t p = head; p >= 0; p = functions_[p].next) {
      const LineFunction& h = functions_[p];
      const Cost d_slope = slope - h.slope;
      const Cost d_offset = offset - h.offset;  // g - h = d_offset + d_slope*x
      Cost a = h.owned_lo;
      Cost b = h.owned_hi;
      if (d_slope == 0) {
        if (d_offset >= 0) {
          continue;
        }
      } else if (d_slope > 0) {
        b = std::min(b, floorDiv(-d_offset - 1, d_slope));
      } else {
        a = std::max(a, floorDiv(d_offset, -d_slope) + 1);
      }
      if (a <= b) {
        u = std::min(u, a);
        w = std::max(w, b);
      }
    }
    if (u > w) {
      return;  // dominated at every integer point
    }
  }
  counters_.functions_inserted++;
  const int32_t id = static_cast<int32_t>(functions_.size());
  // Trim or remove the functions overlapping [u, w]; find the predecessor.
  int32_t prev = -1;
  int32_t before = -1;  // last surviving function left of u
  for (int32_t p = head; p >= 0;) {
    LineFunction& h = functions_[p];
    const int32_t next = h.next;
    if (h.owned_lo >= u && h.owned_hi <= w) {
      h.alive = false;
      counters_.functions_removed++;
      if (prev < 0) {
        head = next;
      } else {
        functions_[prev].next = next;
      }
      p = next;
      continue;
    }
    bool trimmed = false;
    if (h.owned_lo < u && h.owned_hi >= u) {
      h.owned_hi = static_cast<Coord>(u - 1);
      trimmed = true;
    }
    if (h.owned_lo <= w && h.owned_hi > w) {
      h.owned_lo = static_cast<Coord>(w + 1);
      trimmed = true;
    }
    if (trimmed) {
      h.key = std::min(h.at(h.owned_lo), h.at(h.owned_hi));
    }
    if (h.owned_hi < u) {
      before = p;
    }
    prev = p;
    p = next;
  }
  LineFunction g;
  g.offset = offset;
  g.slope = slope;
  g.owned_lo = static_cast<Coord>(u);
  g.owned_hi = static_cast<Coord>(w);
  g.key = std::min(g.at(g.owned_lo), g.at(g.owned_hi));
  g.element = static_cast<int32_t>(element);
  g.alive = true;
  g.processed = false;
  if (before < 0) {
    g.next = head;
    head = id;
  } else {
    g.next = functions_[before].next;
    functions_[before].next = id;
  }
  functions_.push_back(g);
  run.head[element] = head;
  heap_.push_back({g.key, id});
  std::push_heap(heap_.begin(), heap_.end(), std::greater<HeapEntry>());
}

void GeneralPotential::runPropagation(LineRun& run)
{
  while (!heap_.empty()) {
    std::pop_heap(heap_.begin(), heap_.end(), std::greater<HeapEntry>());
    const HeapEntry entry = heap_.back();
    heap_.pop_back();
    const LineFunction& f = functions_[entry.function];
    if (!f.alive || f.processed) {
      counters_.stale_heap_entries++;
      continue;
    }
    if (f.key != entry.key) {
      // The key increased after trimming (Lemma 9): reinsert lazily.
      counters_.stale_heap_entries++;
      heap_.push_back({f.key, entry.function});
      std::push_heap(heap_.begin(), heap_.end(), std::greater<HeapEntry>());
      continue;
    }
    process(run, entry.function);
  }
}

// Propagation of Theorem 13 from a function of a line element to the
// adjacent elements: vias (down / up), the straight run to the neighbouring
// grid line (south / north), and the horizontal run to the next grid point
// (west / east). Grid points and open intervals are separate elements, so
// every move uses the exact edge rate of the closed-tile minimum rule.
void GeneralPotential::process(LineRun& run, const int32_t function)
{
  functions_[function].processed = true;
  counters_.heap_pops++;
  const LineFunction f = functions_[function];
  const Axis& along = alongAxis(run);
  const Axis& lines = lineAxis(run);
  const size_t element = static_cast<size_t>(f.element);
  const int s = static_cast<int>(element % run.num_slots);
  const size_t rest = element / run.num_slots;
  const int m = static_cast<int>(rest % run.num_lines);
  const int z = static_cast<int>(rest / run.num_lines);
  const int line_slot = 2 * m + 1;
  const Coord lo = along.slotLo(s);
  const Coord hi = along.slotHi(s);

  for (const int dz : {-1, 1}) {
    const int z2 = z + dz;
    if (z2 < 0 || z2 >= num_layers_) {
      continue;
    }
    const Cost via = runRate(run, kViaRate, s, line_slot, std::min(z, z2));
    if (via < kInfCost) {
      offer(run, run.element(z2, m, s), lo, hi, f.offset + via, f.slope);
    }
  }
  for (const int dm : {-1, 1}) {
    const int m2 = m + dm;
    if (m2 < 0 || m2 >= run.num_lines) {
      continue;
    }
    const int between = 2 * std::max(m, m2);
    const Cost r = runRate(run, kVerticalRate, s, between, z);
    if (r < kInfCost) {
      const Cost length = static_cast<Cost>(lines.positions[std::max(m, m2)])
                          - lines.positions[std::min(m, m2)];
      offer(run, run.element(z, m2, s), lo, hi, f.offset + r * length, f.slope);
    }
  }
  const int k = s / 2;
  if ((s & 1) == 0) {
    // Open interval: continue to its end points, from the owner of the end.
    const Cost c = runRate(run, kHorizontalRate, s, line_slot, z);
    if (c >= kInfCost) {
      return;
    }
    if (k > 0 && f.owned_lo == lo) {
      const Coord u = along.positions[k - 1];
      offer(run, run.element(z, m, s - 1), u, u, f.at(lo) + c, 0);
    }
    if (k < along.numPositions() && f.owned_hi == hi) {
      const Coord u = along.positions[k];
      offer(run, run.element(z, m, s + 1), u, u, f.at(hi) + c, 0);
    }
    return;
  }
  // Grid point at positions[k].
  const Coord u = along.positions[k];
  const Cost value = f.at(u);
  for (const int side : {-1, 1}) {
    const int interval = s + side;
    const Cost c = runRate(run, kHorizontalRate, interval, line_slot, z);
    if (c >= kInfCost) {
      continue;
    }
    const Coord ilo = along.slotLo(interval);
    const Coord ihi = along.slotHi(interval);
    if (ilo <= ihi) {
      // value + c * |x - u| = value + side * c * (x - u) on this side.
      offer(run,
            run.element(z, m, interval),
            ilo,
            ihi,
            value - side * c * u,
            side * c);
    } else {
      // Adjacent grid points (distance 1).
      const int point = s + 2 * side;
      if (point >= 0 && point < along.numSlots()) {
        const Coord v = along.positions[point / 2];
        offer(run, run.element(z, m, point), v, v, value + c, 0);
      }
    }
  }
}

Cost GeneralPotential::evaluate(const LineRun& run,
                                const int z,
                                const int line,
                                const int slot,
                                const Coord u) const
{
  for (int32_t p = run.head[run.element(z, line, slot)]; p >= 0;
       p = functions_[p].next) {
    const LineFunction& f = functions_[p];
    if (u <= f.owned_hi) {
      return clampCost(f.at(u));
    }
  }
  return kInfCost;
}

Cost GeneralPotential::query(const int z, Coord x, Coord y) const
{
  counters_.queries++;
  if (!has_targets_ || z < 0 || z >= num_layers_) {
    return kInfCost;
  }
  x = std::clamp(x, x_axis_.domain_lo, x_axis_.domain_hi);
  y = std::clamp(y, y_axis_.domain_lo, y_axis_.domain_hi);
  const int x_slot = x_axis_.locate(x);
  const int y_slot = y_axis_.locate(y);
  if (y_slot & 1) {
    return evaluate(horizontal_run_, z, y_slot / 2, x_slot, x);
  }
  if (x_slot & 1) {
    return evaluate(vertical_run_, z, x_slot / 2, y_slot, y);
  }
  if (mode_ == QueryMode::kBoundary) {
    return tileInteriorBoundary(z, x_slot, y_slot, x, y);
  }
  return tileInteriorFunctions(z, x_slot, y_slot, x, y);
}

namespace {

// Via stack costs from layer z to every layer inside one tile interior.
template <typename ViaRate>
void viaChain(const int z,
              const int num_layers,
              ViaRate via_rate,
              std::array<Cost, kMaxPotentialLayers>& chain)
{
  chain[z] = 0;
  for (int a = z + 1; a < num_layers; a++) {
    chain[a] = addCost(chain[a - 1], via_rate(a - 1));
  }
  for (int a = z - 1; a >= 0; a--) {
    chain[a] = addCost(chain[a + 1], via_rate(a));
  }
}

}  // namespace

// Post-Theorem-13 query for a point in the open interior of a tile: a
// shortest path starts with vias inside the tile and then runs straight to
// one of the four sides (Lemma 8), unless it ends in a target by vias only.
Cost GeneralPotential::tileInteriorBoundary(const int z,
                                            const int x_slot,
                                            const int y_slot,
                                            const Coord x,
                                            const Coord y) const
{
  const int kx = x_slot / 2;
  const int ky = y_slot / 2;
  std::array<Cost, kMaxPotentialLayers> chain;
  viaChain(
      z,
      num_layers_,
      [&](int a) { return rate(kViaRate, x_slot, y_slot, a); },
      chain);
  Cost best = kInfCost;
  for (const TargetRectangle& t : targets_) {
    if (t.xlo <= x && x <= t.xhi && t.ylo <= y && y <= t.yhi) {
      best = std::min(best, addCost(chain[t.z], t.terminal_cost));
    }
  }
  const auto& xs = x_axis_.positions;
  const auto& ys = y_axis_.positions;
  for (int a = 0; a < num_layers_; a++) {
    const Cost via = chain[a];
    if (via >= best) {
      continue;
    }
    const Cost rv = rate(kVerticalRate, x_slot, y_slot, a);
    if (rv < kInfCost) {
      if (ky > 0) {
        const Cost e = evaluate(horizontal_run_, a, ky - 1, x_slot, x);
        best = std::min(best,
                        addCost(via, addCost(e, rv * (Cost{y} - ys[ky - 1]))));
      }
      if (ky < y_axis_.numPositions()) {
        const Cost e = evaluate(horizontal_run_, a, ky, x_slot, x);
        best
            = std::min(best, addCost(via, addCost(e, rv * (Cost{ys[ky]} - y))));
      }
    }
    const Cost rh = rate(kHorizontalRate, x_slot, y_slot, a);
    if (rh < kInfCost) {
      if (kx > 0) {
        const Cost e = evaluate(vertical_run_, a, kx - 1, y_slot, y);
        best = std::min(best,
                        addCost(via, addCost(e, rh * (Cost{x} - xs[kx - 1]))));
      }
      if (kx < x_axis_.numPositions()) {
        const Cost e = evaluate(vertical_run_, a, kx, y_slot, y);
        best
            = std::min(best, addCost(via, addCost(e, rh * (Cost{xs[kx]} - x))));
      }
    }
  }
  return best;
}

void GeneralPotential::addTileCandidate(const TileFunction& f,
                                        const Coord xa,
                                        const Coord ya,
                                        const Coord xb,
                                        const Coord yb) const
{
  counters_.tile_candidates++;
  TileCandidate c;
  c.function = f;
  c.corner[0] = f.offset + f.slope_x * xa + f.slope_y * ya;
  c.corner[1] = f.offset + f.slope_x * xa + f.slope_y * yb;
  c.corner[2] = f.offset + f.slope_x * xb + f.slope_y * ya;
  c.corner[3] = f.offset + f.slope_x * xb + f.slope_y * yb;
  auto dominates = [](const TileCandidate& g, const TileCandidate& h) {
    return g.corner[0] <= h.corner[0] && g.corner[1] <= h.corner[1]
           && g.corner[2] <= h.corner[2] && g.corner[3] <= h.corner[3];
  };
  for (const TileCandidate& kept : tile_scratch_) {
    if (dominates(kept, c)) {
      return;
    }
  }
  for (size_t k = 0; k < tile_scratch_.size();) {
    if (dominates(c, tile_scratch_[k])) {
      tile_scratch_[k] = tile_scratch_.back();
      tile_scratch_.pop_back();
    } else {
      k++;
    }
  }
  tile_scratch_.push_back(c);
}

// Sec. 5.1 variant: the candidates of tileInteriorBoundary() as 2-D affine
// functions over the closed tile interior box; an affine function that is
// >= another one at the four corners is >= it on the whole box.
void GeneralPotential::buildTileList(const int z,
                                     const int x_slot,
                                     const int y_slot,
                                     TileListEntry& entry) const
{
  counters_.tile_lists_built++;
  tile_scratch_.clear();
  const int kx = x_slot / 2;
  const int ky = y_slot / 2;
  const Coord xa = x_axis_.slotLo(x_slot);
  const Coord xb = x_axis_.slotHi(x_slot);
  const Coord ya = y_axis_.slotLo(y_slot);
  const Coord yb = y_axis_.slotHi(y_slot);
  std::array<Cost, kMaxPotentialLayers> chain;
  viaChain(
      z,
      num_layers_,
      [&](int a) { return rate(kViaRate, x_slot, y_slot, a); },
      chain);
  for (const TargetRectangle& t : targets_) {
    if (t.xlo <= xa && xa <= t.xhi && t.ylo <= ya && ya <= t.yhi
        && chain[t.z] < kInfCost) {
      addTileCandidate({chain[t.z] + t.terminal_cost, 0, 0}, xa, ya, xb, yb);
    }
  }
  const auto& xs = x_axis_.positions;
  const auto& ys = y_axis_.positions;
  for (int a = 0; a < num_layers_; a++) {
    const Cost via = chain[a];
    if (via >= kInfCost) {
      continue;
    }
    const Cost rv = rate(kVerticalRate, x_slot, y_slot, a);
    if (rv < kInfCost) {
      for (const int side : {0, 1}) {
        const int line = ky - 1 + side;
        if (line < 0 || line >= y_axis_.numPositions()) {
          continue;
        }
        // South (side 0): rv * (y - Y); north (side 1): rv * (Y - y).
        const Cost sign = side == 0 ? 1 : -1;
        for (int32_t p
             = horizontal_run_.head[horizontal_run_.element(a, line, x_slot)];
             p >= 0;
             p = functions_[p].next) {
          const LineFunction& g = functions_[p];
          addTileCandidate(
              {via + g.offset - sign * rv * ys[line], g.slope, sign * rv},
              xa,
              ya,
              xb,
              yb);
        }
      }
    }
    const Cost rh = rate(kHorizontalRate, x_slot, y_slot, a);
    if (rh < kInfCost) {
      for (const int side : {0, 1}) {
        const int line = kx - 1 + side;
        if (line < 0 || line >= x_axis_.numPositions()) {
          continue;
        }
        const Cost sign = side == 0 ? 1 : -1;
        for (int32_t p
             = vertical_run_.head[vertical_run_.element(a, line, y_slot)];
             p >= 0;
             p = functions_[p].next) {
          const LineFunction& g = functions_[p];
          addTileCandidate(
              {via + g.offset - sign * rh * xs[line], sign * rh, g.slope},
              xa,
              ya,
              xb,
              yb);
        }
      }
    }
  }
  entry.stamp = tile_epoch_;
  entry.begin = static_cast<int32_t>(tile_functions_.size());
  entry.count = static_cast<int32_t>(tile_scratch_.size());
  for (const TileCandidate& c : tile_scratch_) {
    tile_functions_.push_back(c.function);
  }
  counters_.tile_functions += entry.count;
}

Cost GeneralPotential::tileInteriorFunctions(const int z,
                                             const int x_slot,
                                             const int y_slot,
                                             const Coord x,
                                             const Coord y) const
{
  const size_t index
      = (static_cast<size_t>(z) * (x_axis_.numPositions() + 1) + x_slot / 2)
            * (y_axis_.numPositions() + 1)
        + y_slot / 2;
  TileListEntry& entry = tile_lists_[index];
  if (entry.stamp != tile_epoch_) {
    buildTileList(z, x_slot, y_slot, entry);
  }
  Cost best = kInfCost;
  const TileFunction* f = tile_functions_.data() + entry.begin;
  for (int32_t k = 0; k < entry.count; k++) {
    best = std::min(best, f[k].offset + f[k].slope_x * x + f[k].slope_y * y);
  }
  return clampCost(best);
}

}  // namespace drt::br
