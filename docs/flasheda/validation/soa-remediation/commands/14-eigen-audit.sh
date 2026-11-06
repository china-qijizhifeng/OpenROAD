#!/bin/bash
set -euo pipefail
TASK=/gpfs/users/linlizhi/flasheda-soa-remediation-20261106
LOCAL=/tmp/flasheda-gpl-timing
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:/opt/or-tools/lib:${LD_LIBRARY_PATH:-}"
echo 'Static Eigen alignment contract and source users'
for arch in default icelake-server; do
  flags=()
  if test "$arch" != default; then flags=(-march="$arch"); fi
  echo "$arch"
  printf '#include <Eigen/Core>\n' | c++ -x c++ -std=c++20 -O3 "${flags[@]}" -I/usr/local/include/eigen3 -dM -E - | grep -E '^#define EIGEN_(MAX_ALIGN_BYTES|DEFAULT_ALIGN_BYTES|IDEAL_MAX_ALIGN_BYTES|MALLOC_ALREADY_ALIGNED|GLIBC_MALLOC_ALREADY_ALIGNED|VECTORIZE_AVX|VECTORIZE_SSE)'
done
sed -n '1,225p' /usr/local/include/eigen3/Eigen/src/Core/util/Memory.h
grep -R -n -E 'Eigen::|Eigen/|EIGEN_MAX_ALIGN_BYTES|EIGEN_MALLOC_ALREADY_ALIGNED' "$LOCAL/OpenROAD/src/gpl" "$LOCAL/OpenROAD/src/rsz" "$LOCAL/OpenROAD/src/dpl" "$LOCAL/OpenROAD/src/grt" | head -120 || true
objdump -d -C --start-address=0x14f0880 --stop-address=0x14f0910 "$TASK/tuned-openroad"
