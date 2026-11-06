# Compact 3D maze scratch validation

The 3D distance, path-length, and destination-mark arrays now use the actual
`layers × y_grid × x_grid` dimensions. Previously they used the legacy 2D
scratch square, whose side is `max(1000, x_grid, y_grid)`. The 3D heap-offset
decoder uses the corresponding actual row and layer strides. The routing costs,
pointer heap, tie order, net order, and congestion policies are unchanged.

The three integer arrays and destination bitmap previously reserved roughly
`12.125 × layers × padded_side²` bytes, excluding allocator overhead. Their new
size is approximately `12.125 × layers × x_grid × y_grid` bytes. The remaining
2D scratch arrays are unchanged. The reduction depends on the design's grid
dimensions; a square grid with a side of at least 1000 does not shrink.

## Reproduce

Build the baseline commit `2847efec57c70e532ca40cf7969e044541a25a94` and the PR
head separately with the same Bazel `--config=opt` toolchain. The benchmark is
Linux-only because it reads CPU ticks from `/proc/self/stat`.

From the repository root, prepare the fixed placement input:

```sh
bench_dir="$PWD/docs/flasheda/validation/grt-maze3d-scratch"
bench_run=/absolute/path/to/new-benchmark-directory
mkdir -p "$bench_run"
gzip -dc "$bench_dir/aes-placed.def.gz" > "$bench_run/aes-placed.def"
```

Run each executable with a separate output directory:

```sh
mkdir -p "$bench_run/baseline-1"
GRT_THREADS=1 \
GRT_LEF="$PWD/test/Nangate45/Nangate45.lef" \
GRT_DEF="$bench_run/aes-placed.def" \
GRT_OUT="$bench_run/baseline-1" \
OMP_DYNAMIC=FALSE OMP_PROC_BIND=FALSE OMP_TOOL=disabled \
/usr/bin/time -v /absolute/path/to/baseline-openroad \
  -no_init -exit "$bench_dir/route.tcl" \
  > "$bench_run/baseline-1/run.log" 2> "$bench_run/baseline-1/resource.txt"
```

Repeat with the candidate executable. Warm up each once, then run baseline,
candidate, candidate, baseline twice serially; repeat with `GRT_THREADS=8`.
Use `GRT_MEASURE wall_s` for `global_route` time and `time -v` for process peak
RSS. File loading and output serialization are outside the GRT time window.
Keep profiler runs separate from this timing sequence.

The frozen AES placement contains 20,297 components and 18,632 nets. It comes
from the public OpenROAD AES placement fixture used in the earlier GPL
experiments. The LEF is the repository's unchanged Nangate45 test LEF. This
experiment has no Liberty/SDC or CTS input and is not a full ORFS timing-driven
flow. It routes metal2–metal6 with a 50% capacity reduction and a maximum of 50
congestion iterations. The expected successful routing has 18,269 routed nets,
zero overflow, wirelength 655,594 µm, and 152,307 vias.

Input SHA-256:

```text
aes-placed.def 596b5a0b5a5f56604b411c525dd61e18034fcccfc331c38267dde3ce24840b86
Nangate45.lef  caa4e9b45e081f93720ce2c4d581497321ddc6fbcb1ffe1b384e324bbd5d43a7
```

Compare the actual `route.guide` and `routed.odb` files for every run, not only
the final wirelength. Existing native GRT regressions additionally cover
rectangular search grids, recovery after heap exhaustion, colocated pins,
overlapping tree edges, resistance-aware routing, and incremental rerouting.
