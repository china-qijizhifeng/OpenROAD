// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once

#include "br/BrAccessFootprint.h"

namespace drt::br {
// A search may change its optimum after either insertion or removal anywhere
// it read. Checking only the selected path's legality is insufficient.
inline bool ripupSearchReusable(const TileBox& reads,
                                const TileBox& writes,
                                const bool same_net_state)
{
  return same_net_state
         && (reads.empty() || writes.empty() || reads.xhi < writes.xlo
             || writes.xhi < reads.xlo || reads.yhi < writes.ylo
             || writes.yhi < reads.ylo);
}
}  // namespace drt::br
