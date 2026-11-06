"""Persist the exact input bytes before the disposable worker is released."""
import hashlib
import json
from pathlib import Path
import tarfile

task=Path('/gpfs/users/linlizhi/flasheda-soa-large-high-remediation-20261106')
source=Path('/tmp/flasheda-gpl-timing/OpenROAD')
manifest=json.loads((task/'input-binary-sha256.json').read_text())
with tarfile.open(task/'frozen-inputs.tar.gz','w:gz',compresslevel=6) as archive:
 for name,expected in manifest.items():
  path=Path(name)
  assert hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==expected,name
  if path in [task/'baseline',task/'candidate']:continue
  relative='source/'+str(path.relative_to(source)) if path.is_relative_to(source) else 'task/'+path.name
  archive.add(path,arcname=relative,recursive=False)
digest=hashlib.file_digest((task/'frozen-inputs.tar.gz').open('rb'),'sha256').hexdigest()
(task/'frozen-inputs.sha256').write_text(f'{digest}  frozen-inputs.tar.gz\n')
print('Frozen exact public design inputs and Tcl bytes:',digest)
