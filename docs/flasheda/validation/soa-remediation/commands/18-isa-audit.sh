#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
cp /tmp/flasheda-gpl-timing/build/compile_commands.json "$TASK/compile-commands-audit.json"
python - <<'PY'
import json
from pathlib import Path
for row in json.loads(Path('/tmp/flasheda-gpl-timing/build/compile_commands.json').read_text()):
 if '/src/gpl/' in row['file'] and Path(row['file']).name in {'fft.cpp','nesterovBase.cpp','initialPlace.cpp','solver.cpp','replace.cpp'}:
  print(json.dumps(row))
PY
