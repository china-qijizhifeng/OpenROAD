#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
python - <<'PY'
from pathlib import Path
import subprocess
p=Path('/gpfs/users/linlizhi/flasheda-soa-remediation-20261106')
symbols=subprocess.check_output(['nm','-C',str(p/'density-subset-adapter.o')],text=True)
assert 'DensityBackendFullSubsetAndReleaseSizes_Test::TestBody()' in symbols
root=p/'density-subset-runfiles';assert (root/'_main').resolve()==Path('/tmp/flasheda-gpl-timing/OpenROAD')
r=subprocess.run([str(p/'density-subset-fixture'),'--gtest_filter=TimingRepairIntegrationTest.DensityBackendFullSubsetAndReleaseSizes'],cwd=root,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
print(r.stdout)
assert r.returncode==0 and '[  PASSED  ] 1 test.' in r.stdout and 'Running 1 test from 1 test suite' in r.stdout
PY
