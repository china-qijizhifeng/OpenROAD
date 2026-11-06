#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
LOCAL=/tmp/flasheda-gpl-timing
# Command 15 built successfully; its post-build checker mistook fft_test for GPL.
test "$(cat "$TASK/commands/18-isa-audit.sh.exit")" = 0
python "$TASK/check-kernel-isa.py"
cp "$LOCAL/build/bin/openroad" "$TASK/tuned-openroad"
cp "$LOCAL/build/src/gpl/CMakeFiles/gpl_lib.dir/flags.make" "$TASK/tuned-flags.make"
cp "$LOCAL/build/CMakeCache.txt" "$TASK/tuned-CMakeCache.txt"
cp "$LOCAL/build/compile_commands.json" "$TASK/tuned-compile-commands.json"
sha256sum "$TASK/tuned-openroad" > "$TASK/tuned-binary.sha256"
