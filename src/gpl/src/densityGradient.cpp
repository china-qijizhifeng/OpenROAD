// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Density gradient backends + dispatch. Mirrors wirelengthGradient.cpp.

#include <cassert>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

#include "backendContext.h"
#include "densityGradientBackend.h"
#include "nesterovBase.h"
#include "point.h"

#ifdef ENABLE_GPU
#include "gpu/deviceState.h"
#include "gpu/gpuDensityGradientBackend.h"
#include "gpu/gpuRuntime.h"
#include "gpu/regionDensityField.h"
#endif

namespace gpl {

namespace {

class CpuDensityGradientBackend : public DensityGradientBackend
{
 public:
  explicit CpuDensityGradientBackend(NesterovBase* nb) : nb_(nb) {}

  void getCellGradients(const std::vector<GCellHandle>& gCells,
                        PointColumns& out) override
  {
    if (out.size() != gCells.size()) {
      throw std::invalid_argument(
          "density gradient output size differs from cells");
    }
    if (&gCells == &nb_->getGCells()) {
      nb_->getDensityGradientsColumns(out);
      return;
    }
    // The backend contract accepts subsets as well as the owning region.
#pragma omp parallel for num_threads(nb_->getNbc()->getNumThreads())
    for (size_t i = 0; i < gCells.size(); ++i) {
      out[i] = nb_->getDensityGradient(gCells[i]);
    }
  }

  FloatPoint getCellGradient(const GCell* gCell) override
  {
    return nb_->getDensityGradient(gCell);
  }

  const char* name() const override { return "CPU"; }

 private:
  NesterovBase* nb_;
};

}  // namespace

std::unique_ptr<DensityGradientBackend> makeDensityGradientBackend(
    const BackendContext& ctx)
{
#ifdef ENABLE_GPU
  if (gpuEnabled() && ctx.device_state && ctx.region_field
      && ctx.region_field->numBins() > 0) {
    return std::make_unique<GpuDensityGradientBackend>(
        ctx.nb, ctx.device_state, ctx.region_field);
  }
#endif
  return std::make_unique<CpuDensityGradientBackend>(ctx.nb);
}

}  // namespace gpl
