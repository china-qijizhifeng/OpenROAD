#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-large-high-remediation-20261106
test "$(cat "$TASK/commands/02-acceptance.sh.exit")" = 0
export PATH="/usr/local/bin:$PATH"
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
export OMP_PROC_BIND=FALSE OMP_DYNAMIC=FALSE
unset OMP_PLACES LD_PRELOAD
python3 "$TASK/run-large.py" --task "$TASK" --phase benchmark
python3 "$TASK/freeze-inputs.py"
python3 "$TASK/report.py" --task "$TASK"
touch "$TASK/validation-complete"
touch "$TASK/STOP"
