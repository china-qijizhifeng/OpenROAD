#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
LOCAL=/tmp/flasheda-gpl-timing
test "$(cat "$TASK/commands/26b-profile-fused.sh.exit")" = 0
export PATH="/usr/local/bin:$PATH"
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
tar -xzf "$TASK/source-v3.tar.gz" -C "$LOCAL/OpenROAD"
python - <<'PY'
import hashlib,json,os
from pathlib import Path
task=Path('/gpfs/users/linlizhi/flasheda-soa-remediation-20261106')
root=Path('/tmp/flasheda-gpl-timing/OpenROAD')
for name,expected in json.loads((task/'source-sha256-v3.json').read_text()).items():
 assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected,name
 os.utime(root/name,None)
print('All v3 source inputs verified; mtimes refreshed.')
PY
for unit in storage point density metadata wirelength; do
 c++ -std=c++20 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -ffp-contract=off -I "$LOCAL/OpenROAD/src/gpl/src" "$LOCAL/OpenROAD/src/gpl/test/soa/${unit}_test.cpp" -o "$TASK/v3-$unit-sanitized-test"
 "$TASK/v3-$unit-sanitized-test"
 c++ -std=c++20 -O3 -DNDEBUG -ffp-contract=off -I "$LOCAL/OpenROAD/src/gpl/src" "$LOCAL/OpenROAD/src/gpl/test/soa/${unit}_test.cpp" -o "$TASK/v3-$unit-release-test"
 "$TASK/v3-$unit-release-test"
 echo "v3 $unit ASan/UBSan and NDEBUG passed."
done
c++ -std=c++20 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -ffp-contract=off -fopenmp -I "$LOCAL/OpenROAD/src/gpl/src" "$LOCAL/OpenROAD/src/gpl/test/soa/fft_test.cpp" "$LOCAL/OpenROAD/src/gpl/src/fft.cpp" "$LOCAL/OpenROAD/src/gpl/src/fftsg.cpp" "$LOCAL/OpenROAD/src/gpl/src/fftsg2d.cpp" -o "$TASK/v3-fft-sanitized-test"
"$TASK/v3-fft-sanitized-test"
echo 'v3 FFT ASan/UBSan byte comparisons passed.'
cmake -S "$LOCAL/OpenROAD" -B "$LOCAL/build" -U GPL_ENABLE_IPO -U GPL_CPU_ARCH -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DFETCHCONTENT_FULLY_DISCONNECTED=ON
cmake --build "$LOCAL/build" --target openroad -j20
cp "$LOCAL/build/bin/openroad" "$TASK/soa-v3-openroad"
cp "$LOCAL/build/src/gpl/CMakeFiles/gpl_lib.dir/flags.make" "$TASK/v3-flags.make"
cp "$LOCAL/build/CMakeCache.txt" "$TASK/v3-CMakeCache.txt"
cp "$LOCAL/build/compile_commands.json" "$TASK/v3-compile-commands.json"
sha256sum "$TASK/soa-v3-openroad" > "$TASK/v3-binary.sha256"
python - "$TASK/v3-flags.make" <<'PY'
import sys
from pathlib import Path
s=Path(sys.argv[1]).read_text()
assert '-ffp-contract=off' in s and '-march=' not in s and '-flto' not in s
print('Verified default v3 flags: no IPO/ISA override.')
PY
