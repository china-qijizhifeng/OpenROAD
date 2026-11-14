// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

#include <string>
#include <vector>

#include "core/GCellGraph.h"

namespace grt::helem {

// Readers and writer for the ISPD 2024 global routing contest formats
// (resource file .cap, net file .net, solution file).
bool readIspd24Resource(const std::string& path,
                        GCellGraph& graph,
                        std::string& error);
bool readIspd24Nets(const std::string& path,
                    std::vector<RoutingNet>& nets,
                    std::string& error);
bool writeIspd24Solution(const std::string& path,
                         const std::vector<RoutingNet>& nets,
                         const std::vector<NetRoute3D>& routes,
                         std::string& error);

}  // namespace grt::helem
