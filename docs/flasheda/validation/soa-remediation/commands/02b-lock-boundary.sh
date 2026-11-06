#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
LOCAL=/tmp/flasheda-gpl-timing
test "$(cat "$TASK/commands/02-build.sh.exit")" = 0
export PATH="/usr/local/bin:$PATH"
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
export OMP_PROC_BIND=FALSE OMP_DYNAMIC=FALSE
unset OMP_PLACES
tar -xzf "$TASK/source-v2.tar.gz" -C "$LOCAL/OpenROAD"
for unit in storage density point metadata; do
 c++ -std=c++20 -O2 -I "$LOCAL/OpenROAD/src/gpl/src" "$LOCAL/OpenROAD/src/gpl/test/soa/${unit}_test.cpp" -o "$TASK/$unit-test"
 "$TASK/$unit-test"
done
c++ -std=c++20 -O3 -ffp-contract=off -fopenmp -I "$LOCAL/OpenROAD/src/gpl/src" "$LOCAL/OpenROAD/src/gpl/test/soa/fft_test.cpp" "$LOCAL/OpenROAD/src/gpl/src/fft.cpp" "$LOCAL/OpenROAD/src/gpl/src/fftsg.cpp" "$LOCAL/OpenROAD/src/gpl/src/fftsg2d.cpp" -o "$TASK/fft-test"
"$TASK/fft-test"
echo 'All standalone column and Ooura differential tests passed.'
cmake --build "$LOCAL/build" --target openroad -j20
cp "$LOCAL/build/bin/openroad" "$TASK/soa-final-openroad"
sha256sum "$TASK/soa-final-openroad" > "$TASK/binary.sha256"
