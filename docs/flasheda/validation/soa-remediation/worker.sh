#!/bin/bash
set -uo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
exec > "$TASK/worker.log" 2>&1
trap 'echo "$?" > "$TASK/worker.exit"' EXIT
date -u
while test ! -f "$TASK/STOP"; do
  for command in "$TASK"/commands/*.sh; do
    test -f "$command" || continue
    test ! -f "$command.exit" || continue
    bash "$command" > "$command.log" 2>&1
    echo "$?" > "$command.exit"
  done
  sleep 3
done
