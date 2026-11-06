#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
test "$(cat "$TASK/commands/08-rebuild-inputs.sh.exit")" = 0
c++ -std=c++20 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -ffp-contract=off -fopenmp -I /tmp/flasheda-gpl-timing/OpenROAD/src/gpl/src /tmp/flasheda-gpl-timing/OpenROAD/src/gpl/test/soa/fft_test.cpp /tmp/flasheda-gpl-timing/OpenROAD/src/gpl/src/fft.cpp /tmp/flasheda-gpl-timing/OpenROAD/src/gpl/src/fftsg.cpp /tmp/flasheda-gpl-timing/OpenROAD/src/gpl/src/fftsg2d.cpp -o "$TASK/fft-sanitized-test"
"$TASK/fft-sanitized-test"
echo 'Parallel Ooura ASan/UBSan passed.'
c++ -std=c++20 -O3 -ffp-contract=off -I /tmp/flasheda-gpl-timing/OpenROAD/src/gpl/src -I "$TASK" "$TASK/density-replay.cpp" -o "$TASK/density-replay"
for call in 2 400 1100; do
 "$TASK/density-replay" "/gpfs/users/linlizhi/flasheda-density-cpu-20261106/results/snapshots/aes-$call.txt" > "$TASK/replay-$call.jsonl"
done
