#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
test -f "$TASK/validation-fused-complete"
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
export OMP_PROC_BIND=FALSE OMP_DYNAMIC=FALSE
unset OMP_PLACES
cd "$TASK/profile-tools"
gcc -O2 -Wall -fPIC -shared -o sampler.so sampler.c -lrt
gcc -O2 -Wall -o pcount pcount.c
python "$TASK/profile-fused.py" --old "$TASK" --task "$TASK" --parser /gpfs/users/linlizhi/flasheda-gpl-timing-20261105/source/gpl-timing/validate.py --tcl /gpfs/users/linlizhi/flasheda-openroad-rudy-20261105/source/rudy/place.tcl
python analyze.py "fused-1=$TASK/fork-results/profile-fused/soa-fused-openroad-t1/folded.txt" "fused-8=$TASK/fork-results/profile-fused/soa-fused-openroad-t8/folded.txt" --json "$TASK/phases-fused.json" > "$TASK/phases-fused.md"
