// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// io_pin_area (ORFS finding F2, asap7 KLayout M4.W.5): routed metal on a
// top-level pin must meet the layer's minimum area without the pin shape.
#include <cstdint>
#include <cstdio>
#include <optional>
#include <vector>

#include "br/BrIoPinArea.h"

using namespace drt::br;

namespace {

int fails = 0;

void check(const bool ok, const char* what)
{
  if (!ok) {
    std::printf("FAIL: %s\n", what);
    fails++;
  }
}

}  // namespace

int main()
{
  // asap7 aes key[61]: M4 IO pin (0, 36096)-(84, 36120) at the left die
  // edge, VIA34 top pad 40 x 24 at x = 45; M4 min area 2000, width 24.
  const Box die{0, 0, 48709, 48709};
  const Box pin{0, 36096, 84, 36120};
  const Box pad{25, 36096, 65, 36120};
  const Box far_wire{45, 36240, 2937, 36264};  // the net's M4 wire 6 tracks up
  {
    const std::vector<Box> metal{pad, far_wire};
    const std::vector<size_t> group = piecesOnPin(metal, pin);
    check(group.size() == 1 && group[0] == 0, "only the pad touches the pin");
    IoPinPatchRequest req;
    req.group = {pad};
    req.anchor = {45, 36108};
    req.horizontal = true;
    req.half_width = 12;
    req.min_area = 2000;
    req.grid = 1;
    req.die = die;
    const std::optional<Box> patch = ioPinPatch(req);
    check(patch.has_value(), "a patch for the lone pad");
    if (patch) {
      check(patch->xlo >= die.xlo, "patch inside the die");
      check(patch->xhi > 65, "patch extends toward the die inside");
      check(patch->ylo == 36096 && patch->yhi == 36120, "patch on the track");
      const int64_t area = unionArea({pad, *patch});
      check(area >= 2000, "pad and patch meet the minimum area");
      // KLayout M4.W.5: horizontal width of the M4 polygon at least 44.
      check(patch->merged(pad).width() >= 44, "polygon at least 44 wide");
    }
    // Blocked inside: the other direction leaves the die, so no patch.
    req.legal = [](const Box& b) { return b.xhi <= 65; };
    check(!ioPinPatch(req).has_value(), "no patch when blocked and at edge");
  }
  {
    // A pin at the right die edge extends to the left.
    IoPinPatchRequest req;
    req.group = {{48644, 100, 48684, 124}};
    req.anchor = {48664, 112};
    req.horizontal = true;
    req.half_width = 12;
    req.min_area = 2000;
    req.die = die;
    const std::optional<Box> patch = ioPinPatch(req);
    check(patch.has_value() && patch->xlo < 48644, "right edge extends left");
  }
  {
    // Metal already large enough: no patch.
    IoPinPatchRequest req;
    req.group = {{60, 36096, 2997, 36120}};
    req.anchor = {60, 36108};
    req.half_width = 12;
    req.min_area = 2000;
    req.die = die;
    check(!ioPinPatch(req).has_value(), "long wire needs no patch");
  }
  {
    // Vertical layer, manufacturing grid 5: patch length on the grid.
    IoPinPatchRequest req;
    req.group = {{990, 0, 1010, 30}};
    req.anchor = {1000, 15};
    req.horizontal = false;
    req.half_width = 10;
    req.min_area = 2000;
    req.grid = 5;
    req.die = die;
    const std::optional<Box> patch = ioPinPatch(req);
    check(patch.has_value(), "vertical patch");
    if (patch) {
      check((patch->yhi - 10 - 15) % 5 == 0, "length on the grid");
      check(unionArea({req.group[0], *patch}) >= 2000, "vertical area");
    }
  }
  {
    // Connected pieces: a wire touching the pad counts with it.
    const std::vector<Box> metal{
        pad, {55, 36096, 120, 36120}, {200, 36096, 260, 36120}};
    const std::vector<size_t> group = piecesOnPin(metal, pin);
    check(group.size() == 2, "touching wire joins the group");
  }
  {
    // asap7 aes text_in[25]: the access wire from the pin centre to the via
    // at x = 81 lies inside pin and pad; FlexDR removes it, so it adds no
    // area and the lone pad still needs its patch.
    const Box pin25{0, 37536, 84, 37560};
    const Box pad25{61, 37536, 101, 37560};
    const Box wire25{30, 37536, 93, 37560};
    check(coveredByPinAndPads(wire25, pin25, {pad25}), "access wire covered");
    check(!coveredByPinAndPads({30, 37536, 140, 37560}, pin25, {pad25}),
          "longer wire counts");
  }
  std::printf("io pin area: %d failures\n", fails);
  return fails != 0;
}
