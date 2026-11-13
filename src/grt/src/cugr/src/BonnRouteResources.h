// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Resources of BonnRoute global routing beyond edge capacity and the
// objective, after Gester et al., "BonnRoute: Algorithms and Data Structures
// for Fast and Good VLSI Routing", TODAES 18(2), 2013, Sec. 2.1-2.2 and
// Fig. 1: a net n using edge e may get extra space s(n, e) >= 0 (in tracks)
// next to its wire. Each resource r then has a convex use gamma_{n,e}^r(s):
//
//  - space on e:      gamma(s) = w(n, e) + s                 (linear),
//  - power:           gamma(s) = L (C_ground + C_coupling * g / (g + a s)),
//  - yield loss:      gamma(s) = L x0 (x0 / (g + a s) + x0 / (2 w_e)),
//
// with L the edge length, g the gap to a neighbor wire on the next track
// (pitch - width), a = pitch / 2 the gap growth per track of extra space
// (the space is split evenly between both sides), x0 the minimum defect
// size. Power is the switching capacitance of the wire (activity 1): ground
// capacitance plus coupling to both neighbors, inversely proportional to
// the gap (parallel-plate model). Yield loss is the critical area of the
// wire for a defect density proportional to x0^2 / x^3 (x >= x0): shorts to
// both neighbors give L x0^2 / d (d >= x0, the classical result for this
// density), opens of a wire of width w give L x0^2 / (2 w), which does not
// depend on s. Both decrease and are convex in s, as in Fig. 1.
//
// Eq. (1) prices an edge by min over s of sum_r y_r gamma_r(s) / u_r. For
// the forms above this is a s + K / (g + a s) + C with a, K >= 0, minimized
// in closed form and then over the grid s in {0, step, ..., s_max}; for a
// convex function the best grid point is a neighbor of the continuous
// minimum.
//
// Everything here is plain data and arithmetic, unit tested without a
// database (src/grt/test/cpp/TestBonnRouteResources.cpp).

#pragma once

#include <vector>

namespace grt::bonnroute {

// Per wiring layer. Lengths in DBU, capacitance in fF per DBU of wire.
struct SpaceLayerModel
{
  double gap = 1.0;     // g: gap to the neighbor track at s = 0
  double growth = 1.0;  // a: gap growth per track of extra space
  double ground_capacitance = 0.0;
  double coupling_capacitance = 0.0;  // both neighbors together, at gap g
  double defect_size = 1.0;           // x0
  double wire_width = 1.0;            // w_e for the opens term

  // gamma per unit length of wire (multiply by L).
  double powerPerLength(double s) const;
  double yieldPerLength(double s) const;
  // The coefficient k with powerPerLength(s) = const + k / (g + a s).
  double powerShape() const { return coupling_capacitance * gap; }
  double yieldShape() const { return defect_size * defect_size; }
};

// Model of one layer from its geometry (DBU) and capacitance data as odb
// holds it: `area_capacitance` in pF per um^2 (LEF CPERSQDIST, or the
// per-area value that set_layer_rc writes) and `edge_capacitance` in pF per
// um and edge (LEF EDGECAPACITANCE; 0 after set_layer_rc). The wire's total
// capacitance per length is area * width + 2 * edge; `coupling_fraction` of
// it couples to the neighbor wires at the minimum gap, the rest goes to
// ground. Without any capacitance data the layer gets a unit coupling
// capacitance. `defect_size` is x0.
SpaceLayerModel makeSpaceLayerModel(int width,
                                    int pitch,
                                    int spacing,
                                    double area_capacitance,
                                    double edge_capacitance,
                                    int dbu_per_micron,
                                    double coupling_fraction,
                                    double defect_size);

struct SpaceGrid
{
  double step = 0.5;   // tracks
  double max = 2.0;    // tracks
  int points() const;  // number of grid values, >= 1
  double value(int i) const { return i * step; }
};

struct SpaceChoice
{
  double space = 0.0;  // s
  double cost = 0.0;   // a s + K / (g + a s), without the constant parts
};

// argmin over the grid of space_price * s + shape / (gap + growth * s)
// (space_price, shape >= 0); ties keep the smaller s.
SpaceChoice minimizeSpace(double space_price,
                          double shape,
                          double gap,
                          double growth,
                          const SpaceGrid& grid);

// The same by enumerating every grid point (reference for tests).
SpaceChoice minimizeSpaceBruteForce(double space_price,
                                    double shape,
                                    double gap,
                                    double growth,
                                    const SpaceGrid& grid);

}  // namespace grt::bonnroute
