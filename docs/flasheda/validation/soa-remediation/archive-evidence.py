"""Archive only this task's receipts; large placement outputs/binaries stay on GPFS."""
import hashlib
import json
from pathlib import Path
import tarfile

def digest_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


task = Path('/gpfs/users/linlizhi/flasheda-soa-remediation-20261106')
assert (task / 'worker.exit').read_text().strip() == '0'
for prefix in ['08-', '09-', '10-', '11-', '14-', '15b-', '16b-', '17b-', '22-', '23-', '24b-', '25b-', '26b-', '27-', '28-', '29-', '30-', '31-', '32d-']:
    matches = list((task / 'commands').glob(f'{prefix}*.sh.exit'))
    assert len(matches) == 1 and matches[0].read_text().strip() == '0', matches

hashes = {}
for name in ['baseline', 'previous-soa-openroad', 'soa-final-openroad',
             'tuned-openroad', 'soa-fused-openroad', 'soa-v3-openroad', 'stale-overlay-openroad', 'unsafe-wide-isa-openroad']:
    f = task / name
    hashes[str(f)] = digest_file(f)
inputs = [Path('/gpfs/users/linlizhi/flasheda-gpl-timing-20261105/source/gpl-timing/validate.py'),
          Path('/gpfs/users/linlizhi/flasheda-openroad-rudy-20261105/source/rudy/place.tcl')]
inputs.extend(Path('/gpfs/users/linlizhi/flasheda-density-cpu-20261106/results/snapshots').glob('aes-*.txt'))
for f in inputs:
    hashes[str(f)] = digest_file(f)
(task / 'input-binary-sha256.json').write_text(json.dumps(hashes, indent=2) + '\n')

selected = set()
for f in task.iterdir():
    if f.is_file() and (f.suffix in {'.py', '.sh', '.cpp', '.h', '.yaml', '.json',
                                    '.jsonl', '.log', '.exit', '.sha256', '.make', '.txt', '.md'}
                        or f.name.startswith('source') and f.name.endswith('.tar.gz')
                        or f.name == 'STOP'):
        selected.add(f)
for directory in ['commands', 'fork-results', 'profile-tools', 'profile-fused', 'profile-fused-host']:
    for f in (task / directory).rglob('*'):
        if f.is_file() and f.suffix in {'.py', '.sh', '.log', '.exit', '.json', '.tcl', '.txt', '.c', '.md'}:
            selected.add(f)
with tarfile.open(task / 'evidence.tar.gz', 'w:gz', compresslevel=6) as archive:
    for f in sorted(selected):
        archive.add(f, arcname=str(f.relative_to(task)), recursive=False)
    for f in inputs[:2]:
        archive.add(f, arcname='dependencies/' + f.name, recursive=False)
f = task / 'evidence.tar.gz'
digest = digest_file(f)
(task / 'evidence.sha256').write_text(f'{digest}  evidence.tar.gz\n')
print(f'{len(selected)} receipts; {f.stat().st_size} bytes; SHA-256 {digest}')
