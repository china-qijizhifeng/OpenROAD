// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

namespace utl {
class Logger;
}

namespace drt {
class frDesign;
struct RouterConfiguration;

namespace br {

// Add legal IO metal area patches after FlexDR has fixed the routing.
// The drawn metal must meet area rules without the pin-purpose geometry.
void repairIoPinArea(frDesign* design,
                     utl::Logger* logger,
                     RouterConfiguration* router_cfg);

}  // namespace br
}  // namespace drt
