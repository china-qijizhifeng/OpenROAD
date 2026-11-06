// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once
#include <cstdint>

#include "soaStorage.h"

namespace gpl {
class Instance;
class Pin;
class Net;
class GCell;
class GPin;
class GNet;
enum class CellAreaChange : uint8_t
{
  kNone,
  kRoutability,
  kTimingDriven,
  kNewInstance,
  kDownsize,
  kUpsize,
  kResizeNoChange
};

enum class CellKind : uint8_t
{
  kFiller,
  kStandard,
  kMacro
};

struct GCellColumns
{
  enum Field : size_t
  {
    lx,
    ly,
    ux,
    uy,
    dLx,
    dLy,
    dUx,
    dUy,
    densityScale,
    gradientX,
    gradientY,
    insts,
    gPins,
    change,
    kind,
    locked,
  };
  using Store = SoaStorage<int,
                           int,
                           int,
                           int,
                           int,
                           int,
                           int,
                           int,
                           float,
                           float,
                           float,
                           std::vector<Instance*>,
                           std::vector<GPin*>,
                           CellAreaChange,
                           CellKind,
                           uint8_t>;
};
// Stable row/optimizer-index maps, rebuilt only when the region topology
// changes. They own no numerical copies or pointers into reallocatable columns.
struct GCellRowList
{
  GCellColumns::Store* columns = nullptr;
  std::vector<size_t> rows;
  std::vector<size_t> indices;
};

struct GNetColumns
{
  enum Field : size_t
  {
    lx,
    ly,
    ux,
    uy,
    timingWeight,
    customWeight,
    powerWeight,
    waExpMinSumX,
    waXExpMinSumX,
    waExpMaxSumX,
    waXExpMaxSumX,
    waExpMinSumY,
    waYExpMinSumY,
    waExpMaxSumY,
    waYExpMaxSumY,
    gPins,
    nets,
    isDontCare,
  };
  using Store = SoaStorage<int,
                           int,
                           int,
                           int,
                           float,
                           float,
                           float,
                           float,
                           float,
                           float,
                           float,
                           float,
                           float,
                           float,
                           float,
                           std::vector<GPin*>,
                           std::vector<Net*>,
                           uint8_t>;
};
struct GPinColumns
{
  enum Field : size_t
  {
    offsetCx,
    offsetCy,
    cx,
    cy,
    maxExpSumX,
    maxExpSumY,
    minExpSumX,
    minExpSumY,
    hasMaxExpSumX,
    hasMaxExpSumY,
    hasMinExpSumX,
    hasMinExpSumY,
    gCell,
    gNet,
    pins,
    // WA intermediates are indexed in net CSR order; geometry keeps its row.
    waSlot,
  };
  using Store = SoaStorage<int,
                           int,
                           int,
                           int,
                           float,
                           float,
                           float,
                           float,
                           uint8_t,
                           uint8_t,
                           uint8_t,
                           uint8_t,
                           GCell*,
                           GNet*,
                           std::vector<Pin*>,
                           size_t>;
};
struct BinColumns
{
  enum Field : size_t
  {
    x,
    y,
    lx,
    ly,
    ux,
    uy,
    nonPlaceArea,
    instPlacedArea,
    instPlacedAreaUnscaled,
    nonPlaceAreaUnscaled,
    fillerArea,
    density,
    targetDensity,
    electroPhi,
    electroFieldX,
    electroFieldY,
  };
  using Store = SoaStorage<int,
                           int,
                           int,
                           int,
                           int,
                           int,
                           int64_t,
                           int64_t,
                           int64_t,
                           int64_t,
                           int64_t,
                           float,
                           float,
                           float,
                           float,
                           float>;
};
}  // namespace gpl
