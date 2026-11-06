#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
python - <<'PY'
from pathlib import Path
import shutil
root=Path('/tmp/flasheda-gpl-timing/OpenROAD/src/tst')
out=Path('/gpfs/users/linlizhi/flasheda-soa-remediation-20261106/fixture-audit/tst');out.mkdir(parents=True,exist_ok=True)
for f in root.rglob('*'):
 if f.is_file() and (f.name in ['CMakeLists.txt','BUILD'] or 'IntegratedFixture' in f.name):
  print(f)
  d=out/f.relative_to(root);d.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,d)
for f in Path('/tmp/flasheda-gpl-timing/build/src/tst').rglob('*.a'):
 print('built tst library:',f)
PY
