#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
test "$(cat "$TASK/commands/23-fused-acceptance.sh.exit")" = 0
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
export OMP_PROC_BIND=FALSE OMP_DYNAMIC=FALSE
unset OMP_PLACES
python "$TASK/benchmark-fused.py" --task "$TASK"
