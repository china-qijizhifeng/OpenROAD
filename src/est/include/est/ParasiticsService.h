// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026-2026, The OpenROAD Authors

#pragma once

#include <utility>
#include <vector>

#include "grt/GRoute.h"

namespace est {

// Service interface published by the parasitics estimator. Consumers
// look it up through utl::ServiceRegistry and do not depend on the
// concrete est::EstimateParasitics type.
class ParasiticsService
{
 public:
  virtual ~ParasiticsService() = default;

  // Clear all existing parasitics and re-estimate them from the
  // current global routing state (partial routes). With threads > 1 the
  // nets are estimated in parallel; the result is identical to one
  // thread. Only a caller that knows nothing else reads or writes
  // parasitics or delays meanwhile may pass threads > 1.
  virtual void estimateAllGlobalRouteParasitics(int threads = 1) = 0;

  // Re-estimate one net's parasitics from the given route; no-op if empty.
  virtual void updateGlobalRouteParasitics(odb::dbNet* net, grt::GRoute& route)
      = 0;

  // Replace the parasitics of the given nets with estimates from the given
  // partial routes, keep the parasitics of every other net and invalidate
  // all delays, as estimateAllGlobalRouteParasitics() does. With every net
  // whose partial route changed since the previous estimate in `routes`, the
  // parasitics and the timing equal those of
  // estimateAllGlobalRouteParasitics(). Works only while the partial routes
  // are FastRoute's routes of every net
  // (grt::GlobalRouter::partialRoutesAreFastRouteRoutes()); returns whether
  // it did.
  virtual bool updatePartialRouteParasitics(
      std::vector<std::pair<odb::dbNet*, grt::GRoute>>& routes,
      int threads)
      = 0;
};

}  // namespace est
