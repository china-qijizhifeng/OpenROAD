#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
LOCAL=/tmp/flasheda-gpl-timing
test "$(cat "$TASK/commands/14-eigen-audit.sh.exit")" = 0
export PATH="/usr/local/bin:$PATH"
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
python - <<'PY'
from pathlib import Path
names=['nesterovBase.cpp','fft.cpp','fftsg.cpp','fftsg2d.cpp','hpwl.cpp','wirelengthGradient.cpp','wirelengthSoa.cpp','densityGradient.cpp']
root=Path('/tmp/flasheda-gpl-timing/build/src/gpl/CMakeFiles/gpl_lib.dir/src')
for name in names:
 assert 'eigen3/' not in (root/(name+'.o.d')).read_text().lower(), name
print('Static dependency audit: all specialized kernels are Eigen-free.')
PY
cp "$TASK/tuned-openroad" "$TASK/unsafe-wide-isa-openroad"
cp "$TASK/tuned-flags.make" "$TASK/unsafe-wide-isa-flags.make"
cp "$TASK/tuned-CMakeCache.txt" "$TASK/unsafe-wide-isa-CMakeCache.txt"
cp "$TASK/CMakeLists-kernel-isa.txt" "$LOCAL/OpenROAD/src/gpl/CMakeLists.txt"
python - <<'PY'
import hashlib,json
from pathlib import Path
task=Path('/gpfs/users/linlizhi/flasheda-soa-remediation-20261106')
root=Path('/tmp/flasheda-gpl-timing/OpenROAD')
for name, expected in json.loads((task/'source-sha256-kernel-isa.json').read_text()).items():
 assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected,name
print('All final checkout source hashes verified.')
PY
cmake -S "$LOCAL/OpenROAD" -B "$LOCAL/build" -DGPL_ENABLE_IPO=ON -DGPL_CPU_ARCH=icelake-server -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DFETCHCONTENT_FULLY_DISCONNECTED=ON
cmake --build "$LOCAL/build" --target openroad -j20
python "$TASK/check-kernel-isa.py"
cp "$LOCAL/build/bin/openroad" "$TASK/tuned-openroad"
cp "$LOCAL/build/src/gpl/CMakeFiles/gpl_lib.dir/flags.make" "$TASK/tuned-flags.make"
cp "$LOCAL/build/CMakeCache.txt" "$TASK/tuned-CMakeCache.txt"
cp "$LOCAL/build/compile_commands.json" "$TASK/tuned-compile-commands.json"
sha256sum "$TASK/tuned-openroad" > "$TASK/tuned-binary.sha256"
