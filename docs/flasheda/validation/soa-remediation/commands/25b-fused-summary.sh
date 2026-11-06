#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
test "$(cat "$TASK/commands/24b-fused-benchmark.sh.exit")" = 0
python - "$TASK" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);d={}
b=json.loads((p/'fork-results/benchmark-fused/summary.json').read_text());assert b['complete'] and len(b['runs'])==12
d['benchmark']=b
for label,size in [('fused',4),('cases-fused',10)]:
 results=json.loads((p/'fork-results'/label/'summary.json').read_text())
 assert len(results)==size
 if label=='cases-fused':assert all(len(v)==2 and all(r['exact'] for r in v.values()) for v in results.values())
 d[label]=results
(p/'summary-fused.json').write_text(json.dumps(d,indent=2)+'\n')
print(json.dumps(b['summary'],indent=2))
PY
touch "$TASK/validation-fused-complete"
