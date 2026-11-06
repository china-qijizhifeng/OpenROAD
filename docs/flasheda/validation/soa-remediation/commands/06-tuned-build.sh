#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
LOCAL=/tmp/flasheda-gpl-timing
test "$(cat "$TASK/commands/05-benchmark.sh.exit")" = 0
export PATH="/usr/local/bin:$PATH"
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
cmake -S "$LOCAL/OpenROAD" -B "$LOCAL/build" -DGPL_ENABLE_IPO=ON -DGPL_CPU_ARCH=icelake-server -DFETCHCONTENT_FULLY_DISCONNECTED=ON
cmake --build "$LOCAL/build" --target openroad -j20
cp "$LOCAL/build/bin/openroad" "$TASK/tuned-openroad"
cp "$LOCAL/build/src/gpl/CMakeFiles/gpl_lib.dir/flags.make" "$TASK/tuned-flags.make"
cp "$LOCAL/build/CMakeCache.txt" "$TASK/tuned-CMakeCache.txt"
sha256sum "$TASK/tuned-openroad" > "$TASK/tuned-binary.sha256"
python - "$TASK/tuned-flags.make" <<'PYFLAGS'
import sys
from pathlib import Path
flags=Path(sys.argv[1]).read_text()
for required in ['-flto', '-march=icelake-server', '-ffp-contract=off']:
 if required not in flags: raise RuntimeError(f'Missing tuned build flag: {required}')
print('Verified GPL-only IPO, architecture and float-contract flags.')
PYFLAGS
