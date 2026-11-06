#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
nm -C "$TASK/density-subset-fixture.o" | grep -E 'DensityBackend|RegisterTest|GLOBAL__sub|test_info' | head -25 || true
nm -C "$TASK/density-subset-fixture" | grep -E 'DensityBackend|RegisterTest|TimingRepairIntegrationTest|GLOBAL__sub_I' | head -30 || true
cat /tmp/flasheda-gpl-timing/build/src/CMakeFiles/openroad.dir/link.txt
