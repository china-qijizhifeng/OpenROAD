#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
test "$(cat "$TASK/commands/16b-kernel-isa-acceptance.sh.exit")" = 0
python "$TASK/summarize.py" --task "$TASK"
python - "$TASK" <<'PY'
import json,os,platform,subprocess,sys
from pathlib import Path
p=Path(sys.argv[1])
d={'python':sys.version,'platform':platform.platform(),'hostname':platform.node(),'affinity':sorted(os.sched_getaffinity(0)),'compiler':subprocess.check_output(['c++','--version'],text=True),'cmake':subprocess.check_output(['cmake','--version'],text=True)}
try:
 import google.protobuf
 from google.protobuf.internal import api_implementation
 d.update(protobuf=google.protobuf.__version__,protobuf_implementation=api_implementation.Type())
except ImportError: d['protobuf']='unavailable (not used by C++ validation)'
for f in ['/sys/fs/cgroup/cpu.max','/sys/fs/cgroup/memory.max']:
 if Path(f).exists():d[f]=Path(f).read_text().strip()
(p/'runtime.json').write_text(json.dumps(d,indent=2)+'\n')
PY
touch "$TASK/validation-v1-complete"
