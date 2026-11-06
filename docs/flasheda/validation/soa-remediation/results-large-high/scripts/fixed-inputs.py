import hashlib
import json
import os
import platform
from pathlib import Path
import subprocess
import sys

task=Path('/gpfs/users/linlizhi/flasheda-soa-large-high-remediation-20261106')
source=Path('/tmp/flasheda-gpl-timing/OpenROAD')
paths=[task/'baseline',task/'candidate',task/'state.tcl',task/'gpl-only.tcl']
paths += [source/'src/gpl/test'/name for name in ['nangate45.lef','large01.def','large02.def']]
hashes={str(p):hashlib.file_digest(p.open('rb'),'sha256').hexdigest() for p in paths}
assert hashes[str(task/'candidate')]=='4d627682245a2896c22ce73248c5e122ad5b6f01b3f73b387f8bb5f1c07a4252'
assert hashes[str(task/'baseline')]=='d3e694dc50bde8d8c848fe153ff7d64e55ce6a96d9d9b4548cc6484ab07f1c0f'
(task/'input-binary-sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
d={'python':sys.version,'platform':platform.platform(),'hostname':platform.node(),'affinity':sorted(os.sched_getaffinity(0)),'compiler':subprocess.check_output(['c++','--version'],text=True)}
try:
 import google.protobuf
 from google.protobuf.internal import api_implementation
 d.update(protobuf=google.protobuf.__version__,protobuf_implementation=api_implementation.Type())
except ImportError:d['protobuf']='unavailable (unused by CPU C++ verification)'
for f in ['/sys/fs/cgroup/cpu.max','/sys/fs/cgroup/memory.max']:
 if Path(f).exists():d[f]=Path(f).read_text().strip()
(task/'runtime.json').write_text(json.dumps(d,indent=2)+'\n')
print('Binary hashes and all design input hashes fixed.')
