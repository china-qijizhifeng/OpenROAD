#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
printf 'addr2line='; command -v addr2line || true
printf 'g++='; command -v c++
find /tmp/flasheda-gpl-timing/OpenROAD/src/gpl/test -maxdepth 2 -iname '*timing*' -o -iname '*integration*'
python - <<'PYCOPY'
from pathlib import Path
import shutil
root=Path('/tmp/flasheda-gpl-timing/OpenROAD/src/gpl/test')
out=Path('/gpfs/users/linlizhi/flasheda-soa-remediation-20261106/fixture-audit');out.mkdir(exist_ok=True)
for f in root.rglob('*.cpp'):
 if 'TimingRepairIntegrationTest' in f.read_text(errors='replace'):
  print('TimingRepairIntegrationTest fixture:',f)
  shutil.copy2(f,out/f.name)
for name in ['CMakeLists.txt','BUILD']:
 f=root/name
 if f.exists():shutil.copy2(f,out/name)
PYCOPY
