#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
test "$(cat "$TASK/commands/30-v3-benchmark.sh.exit")" = 0
python - "$TASK" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);d={}
b=json.loads((p/'fork-results/benchmark-v3/summary.json').read_text());assert b['complete'] and len(b['runs'])==12
d['benchmark']=b
for label,size in [('v3',4),('cases-v3',10),('v3-consistency',4),('cases-v3-consistency',3)]:
 results=json.loads((p/'fork-results'/label/'summary.json').read_text())
 assert len(results)==size
 if label.startswith('cases-'):assert all(len(v)==2 and all(r['exact'] for r in v.values()) for v in results.values())
 d[label]=results
(p/'summary-v3.json').write_text(json.dumps(d,indent=2)+'\n')
print(json.dumps(b['summary'],indent=2))
PY
touch "$TASK/validation-v3-complete"
