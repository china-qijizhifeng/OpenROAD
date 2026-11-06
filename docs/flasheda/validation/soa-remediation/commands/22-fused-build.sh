#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
LOCAL=/tmp/flasheda-gpl-timing
test -f "$TASK/validation-v1-complete"
export PATH="/usr/local/bin:$PATH"
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
tar -xzf "$TASK/source-fused.tar.gz" -C "$LOCAL/OpenROAD"
python - <<'PY'
import hashlib,json,os
from pathlib import Path
task=Path('/gpfs/users/linlizhi/flasheda-soa-remediation-20261106')
root=Path('/tmp/flasheda-gpl-timing/OpenROAD')
for name,expected in json.loads((task/'source-sha256-fused.json').read_text()).items():
 assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected,name
 os.utime(root/name,None)
print('Final fused source input hashes verified and mtimes refreshed.')
PY
cmake -S "$LOCAL/OpenROAD" -B "$LOCAL/build" -DGPL_ENABLE_IPO=OFF -DGPL_CPU_ARCH= -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DFETCHCONTENT_FULLY_DISCONNECTED=ON
cmake --build "$LOCAL/build" --target openroad -j20
cp "$LOCAL/build/bin/openroad" "$TASK/soa-fused-openroad"
cp "$LOCAL/build/src/gpl/CMakeFiles/gpl_lib.dir/flags.make" "$TASK/fused-flags.make"
cp "$LOCAL/build/CMakeCache.txt" "$TASK/fused-CMakeCache.txt"
sha256sum "$TASK/soa-fused-openroad" > "$TASK/fused-binary.sha256"
python - "$TASK/fused-flags.make" <<'PY'
import sys
from pathlib import Path
s=Path(sys.argv[1]).read_text()
assert '-ffp-contract=off' in s and '-march=' not in s and '-flto' not in s
print('Verified portable default GPL build flags; no IPO or architecture override.')
PY
c++ -std=c++20 -O3 -ffp-contract=off -I "$LOCAL/OpenROAD/src/gpl/src" -I "$TASK" "$TASK/wirelength-replay.cpp" -o "$TASK/wirelength-replay"
"$TASK/wirelength-replay" > "$TASK/wirelength-replay.jsonl"
