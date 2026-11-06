#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
LOCAL=/tmp/flasheda-gpl-timing
# Previous acceptance stopped at incremental02: overlay mtimes reused stale objects.
# Keep every failed log and rebuild all overlaid inputs before resuming.
export PATH="/usr/local/bin:$PATH"
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
export OMP_PROC_BIND=FALSE OMP_DYNAMIC=FALSE
unset OMP_PLACES
tar -xzf "$TASK/source-final.tar.gz" -C "$LOCAL/OpenROAD"
python "$TASK/verify-source-final.py" --touch
cmake --build "$LOCAL/build" --target openroad -j20
cp "$LOCAL/build/bin/openroad" "$TASK/soa-final-openroad"
sha256sum "$TASK/soa-final-openroad" > "$TASK/binary.sha256"
python "$TASK/verify-source-final.py"
python - "$TASK/soa-final-openroad" <<'PYSYMBOL'
import subprocess,sys
symbols=subprocess.check_output(['nm','-C',sys.argv[1]], text=True)
for required in ['gpl::NesterovBaseCommon::refreshCellLocks()', 'gpl::NesterovBase::updateNextCoordinatesColumns(float)']:
 if required not in symbols: raise RuntimeError(f'Missing final implementation: {required}')
print('Final default binary contains lock-refresh and column-coordinate implementations.')
PYSYMBOL
python "$TASK/validate-resident.py" --binary "$TASK/soa-final-openroad" --label final
python "$TASK/validate-cases.py" --baseline "$TASK/baseline" --candidate "$TASK/soa-final-openroad" --source "$LOCAL/OpenROAD" --output "$TASK/fork-results/cases-final" --baseline-results "$TASK/fork-results/cases"
