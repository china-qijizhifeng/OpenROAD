#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
test "$(cat "$TASK/commands/02b-lock-boundary.sh.exit")" = 0
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
export OMP_PROC_BIND=FALSE OMP_DYNAMIC=FALSE
unset OMP_PLACES
python "$TASK/verify-source.py"
python "$TASK/validate-resident.py" --binary "$TASK/soa-final-openroad" --label first
python "$TASK/validate-cases.py" --baseline "$TASK/baseline" --candidate "$TASK/soa-final-openroad" --source /tmp/flasheda-gpl-timing/OpenROAD --output "$TASK/fork-results/cases"
