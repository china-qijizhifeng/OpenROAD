#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-large-high-remediation-20261106
DEPS=/gpfs/users/linlizhi/flasheda-sta-online-20261102
for archive in cmake-deps sysdeps flex-header; do tar -xzf "$DEPS/$archive.tar.gz" -C /; done
ldconfig
tar -xzf /gpfs/users/linlizhi/flasheda-density-cpu-20261106/fork-build-cache.tar.gz -C /tmp flasheda-gpl-timing/OpenROAD
lscpu > "$TASK/lscpu.txt"
python3 "$TASK/fixed-inputs.py"
