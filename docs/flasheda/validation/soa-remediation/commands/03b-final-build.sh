#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
LOCAL=/tmp/flasheda-gpl-timing
test "$(cat "$TASK/commands/03-acceptance.sh.exit")" = 0
export PATH="/usr/local/bin:$PATH"
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
export OMP_PROC_BIND=FALSE OMP_DYNAMIC=FALSE
unset OMP_PLACES
tar -xzf "$TASK/source-final.tar.gz" -C "$LOCAL/OpenROAD"
cmake --build "$LOCAL/build" --target openroad -j20
cp "$LOCAL/build/bin/openroad" "$TASK/soa-final-openroad"
sha256sum "$TASK/soa-final-openroad" > "$TASK/binary.sha256"
python "$TASK/verify-source-final.py"
python "$TASK/validate-resident.py" --binary "$TASK/soa-final-openroad" --label final
python "$TASK/validate-cases.py" --baseline "$TASK/baseline" --candidate "$TASK/soa-final-openroad" --source "$LOCAL/OpenROAD" --output "$TASK/fork-results/cases-final" --baseline-results "$TASK/fork-results/cases"
