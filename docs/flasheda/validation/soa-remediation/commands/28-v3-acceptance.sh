#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
test "$(cat "$TASK/commands/27-v3-build.sh.exit")" = 0
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
export OMP_PROC_BIND=FALSE OMP_DYNAMIC=FALSE
unset OMP_PLACES
python "$TASK/validate-resident-v3.py" --binary "$TASK/soa-v3-openroad" --label v3
python "$TASK/validate-cases-v3.py" --baseline "$TASK/baseline" --candidate "$TASK/soa-v3-openroad" --source /tmp/flasheda-gpl-timing/OpenROAD --output "$TASK/fork-results/cases-v3" --baseline-results "$TASK/fork-results/cases-final"
