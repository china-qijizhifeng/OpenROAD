#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
DEPS=/gpfs/users/linlizhi/flasheda-sta-online-20261102
for archive in cmake-deps sysdeps flex-header; do tar -xzf "$DEPS/$archive.tar.gz" -C /; done
ldconfig
tar -xzf /gpfs/users/linlizhi/flasheda-density-cpu-20261106/fork-build-cache.tar.gz -C /tmp
lscpu > "$TASK/lscpu.txt"
cp /gpfs/users/linlizhi/flasheda-density-cpu-20261106/baseline "$TASK/baseline"
cp /gpfs/users/linlizhi/flasheda-density-cpu-20261106/soa-final-openroad "$TASK/previous-soa-openroad"
