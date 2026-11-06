import argparse
import os
import hashlib
import json
from pathlib import Path
parser = argparse.ArgumentParser()
parser.add_argument('--touch', action='store_true')
args = parser.parse_args()
root = Path('/tmp/flasheda-gpl-timing/OpenROAD')
task = Path('/gpfs/users/linlizhi/flasheda-soa-remediation-20261106')
for name, expected in json.loads((task/'source-sha256-final.json').read_text()).items():
    actual = hashlib.sha256((root/name).read_bytes()).hexdigest()
    if actual != expected:
        raise RuntimeError(f'Source mismatch: {name}')
    if args.touch:
        os.utime(root/name, None)
print('All production and test source hashes match the local checkout.')
for name in ['src/gpl/CMakeFiles/gpl_lib.dir/flags.make', 'src/CMakeFiles/openroad.dir/link.txt', 'CMakeCache.txt']:
    f = Path('/tmp/flasheda-gpl-timing/build')/name
    if f.exists(): (task/f.name).write_bytes(f.read_bytes())
