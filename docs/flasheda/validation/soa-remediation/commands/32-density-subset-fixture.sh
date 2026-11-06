#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
LOCAL=/tmp/flasheda-gpl-timing
test "$(cat "$TASK/commands/31-v3-summary.sh.exit")" = 0
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
python - <<'PY'
from pathlib import Path
import json,shlex,subprocess
task=Path('/gpfs/users/linlizhi/flasheda-soa-remediation-20261106')
build=Path('/tmp/flasheda-gpl-timing/build')
root=Path('/tmp/flasheda-gpl-timing/OpenROAD')
commands=json.loads((build/'compile_commands.json').read_text())
target='src/gpl/CMakeFiles/gpl_lib.dir/src/nesterovBase.cpp.o'
row=next(r for r in commands if r.get('output')==target)
args=shlex.split(row['command']);args=args[:args.index('-o')]
args += ['-I'+str(root),'-I'+str(root/'src/tst/include')]
args += ['-isystem','/usr/local/include','-c',str(task/'density-subset-fixture-test.cpp'),'-o',str(task/'density-subset-fixture.o')]
print('Compile existing TimingRepairIntegrationTest adapter with Release/NDEBUG flags.',flush=True)
subprocess.run(args,cwd=row['directory'],check=True)
link=shlex.split((build/'src/CMakeFiles/openroad.dir/link.txt').read_text())
libraries=link[link.index('-o')+2:]
for source in (root/'src/tst/src').glob('*.cpp'):
 if source.name not in ['db_fixture.cpp','logger_checks.cpp','fixture.cpp','loaded_design.cpp','IntegratedFixture.cpp']:
  continue
 command=shlex.split(row['command']);command=command[:command.index('-o')]
 obj=task/('density-subset-'+source.stem+'.o')
 command += ['-I'+str(root/'src/tst/include'),'-c',str(source),'-o',str(obj)]
 subprocess.run(command,cwd=row['directory'],check=True)
objects=[str(task/'density-subset-fixture.o')]+[str(p) for p in task.glob('density-subset-*.o') if p.name!='density-subset-fixture.o']
gtest=['/usr/local/lib/libgtest_main.a','/usr/local/lib/libgtest.a']
assert all(Path(p).exists() for p in gtest)
subprocess.run(['/usr/bin/c++','-o',str(task/'density-subset-fixture'),*objects,'-Wl,--start-group',*libraries,*gtest,'-Wl,--end-group'],cwd=build/'src',check=True)
runfiles=task/'density-subset-runfiles';runfiles.mkdir(exist_ok=True)
(runfiles/'_main').symlink_to(root, target_is_directory=True)
subprocess.run([str(task/'density-subset-fixture'),'--gtest_filter=TimingRepairIntegrationTest.DensityBackendFullSubsetAndReleaseSizes'],cwd=runfiles,check=True)
PY
