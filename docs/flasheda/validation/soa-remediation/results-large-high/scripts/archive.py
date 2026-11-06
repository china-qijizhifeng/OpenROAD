"""Archive receipts and exact input bytes, preserving binaries/raw artifacts on GPFS."""
import hashlib
import json
from pathlib import Path
import tarfile

def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()

task=Path('/gpfs/users/linlizhi/flasheda-soa-large-high-remediation-20261106')
assert (task/'worker.exit').read_text().strip()=='0'
for prefix in ['01-','02-','03-']:
 files=list((task/'commands').glob(prefix+'*.sh.exit'))
 assert len(files)==1 and files[0].read_text().strip()=='0',files
data=json.loads((task/'results/summary.json').read_text())
assert data['complete'] and len(data['runs'])==36
for r in data['runs']:
 f=Path(r['log']);assert digest(f)==r['log_sha256']
for name,expected_digest in data['binary_sha256'].items():
 assert digest(task/name)==expected_digest,name
files=[]
for f in task.rglob('*'):
 if not f.is_file():continue
 if f.suffix in {'.py','.sh','.log','.exit','.json','.tcl','.txt','.make','.yaml','.md','.sha256'} or f.name in ['STOP','frozen-inputs.tar.gz','source-v2.tar.gz','source-v3.tar.gz']:
  files.append(f)
with tarfile.open(task/'evidence.tar.gz','w:gz',compresslevel=6) as archive:
 for f in sorted(files):archive.add(f,arcname=str(f.relative_to(task)),recursive=False)
archive_digest=digest(task/'evidence.tar.gz')
(task/'evidence.sha256').write_text(f'{archive_digest}  evidence.tar.gz\n')
print(len(files),'receipts;',archive_digest)
