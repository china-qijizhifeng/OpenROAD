#!/bin/bash
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, The OpenROAD Authors
#
# Builds the database-independent off-track modules with the host compiler,
# runs the randomized test and the measurement binary on generated standard
# cell rows, and writes everything to the output file.
# usage: bonnroute_offtrack_measure.sh <output file>
set -euo pipefail
OUT=${1:-offtrack-measurements.txt}
DRT=$(cd "$(dirname "$0")/.." && pwd)
BIN=$(mktemp -d)
SRCS="$DRT/src/br/BrBlockageGrid.cpp $DRT/src/br/BrOffTrackSearch.cpp
      $DRT/src/br/BrPinAccess.cpp $DRT/src/br/BrPathCleanup.cpp"
CXX=${CXX:-g++}
$CXX -std=c++20 -O2 -I"$DRT/src" -o "$BIN/test" \
  "$DRT/test/BonnRouteOffTrackTest.cpp" $SRCS
$CXX -std=c++20 -O2 -I"$DRT/src" -o "$BIN/bench" \
  "$DRT/test/BonnRouteOffTrackBench.cpp" $SRCS
{
  echo "# $(date -u +%FT%TZ) $(git -C "$DRT" rev-parse --short HEAD 2>/dev/null) $($CXX --version | head -1)"
  echo "## randomized test"
  "$BIN/test"
  echo "## bench: rows cells/row seed row-patterns"
  for args in "12 40 1 0" "12 40 2 0" "24 60 3 0" "12 40 1 4" "24 60 3 8"; do
    "$BIN/bench" $args
  done
  echo "## gcd-like pin access: region sources vs access points"
  "$BIN/bench" gcd 12 30 1
  "$BIN/bench" gcd 24 40 2
} | tee "$OUT"
