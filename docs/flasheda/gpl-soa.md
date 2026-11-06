# Resident GPL SoA

This experimental fork starts from OpenROAD
`18216f3506f090aa7e12f83620b47f21ddff4f40`.
The migration changes storage and access, not the placement objective or the
numerical rules. CPU acceptance is required before publishing performance claims.

## Ownership

All fields of the seven repeated placement entity types (`Instance`, `Pin`,
`Net`, `GCell`, `GPin`, `GNet`, `Bin`) live in typed column stores. This includes
geometry, flags, weights, OpenDB associations, ownership pointers, and adjacency
lists. Entity objects carry only a row handle. Physical wrappers use a
PlacerBaseCommon pool, virtual cells/pins/nets use NesterovBaseCommon pools,
fillers use a region pool, and bins use a grid pool. Standalone temporary objects
can own a private store; population constructors explicitly use their owner pool.

`SoaRow` has value semantics: a copy allocates an independent row and deep-copies
list fields; a move transfers the row. This preserves filler snapshot/restore and
swap-remove semantics. Columns may reallocate at topology boundaries, so a handle
holds a row index, never an element pointer. Allocation and release are serial;
parallel kernels operate on distinct rows, private integer grids, or existing
synchronized outputs. Field vectors containing pointers preserve OpenDB boundary
associations without copying the database into the numerical core.

Net-to-pin and cell-to-pin numerical access uses flat pin-row arrays plus CSR
start offsets, rebuilt after construction and `fixPointers()` following committed
timing repair. Row IDs are stable across object vector relocation; callback
swap-removes therefore do not silently change CSR meaning. Pin order inside each
net/cell is retained. Debug tracing retains the original object-facing algorithm,
but its fields are backed by the same authoritative columns.

Nesterov coordinates, gradients, backtracking states, routability snapshots,
removed-filler snapshots, and IO optimizer coordinates use separate x/y columns.
Minimum routability cell sizes use width/height columns. Initial placement's
Eigen sparse matrices and separate x/y vectors remain the solver representation.
Singleton configuration, local scalar return values and OpenDB-owned structures
are not converted into artificial one-row arrays. IO constraint/slot geometry
and auxiliary routing data are not covered by the seven-entity migration.

## Density pipeline

Resident bins use `x * ny + y`, matching the existing Poisson implementation.
Legacy Bin iteration views remain y-major to preserve scalar reduction order.
Scatter and gather load cell rectangle fields once, then visit bin columns
without Bin objects or a per-iteration population pack. One thread uses ordinary
integer addition. Multiple threads use reusable private integer grids, with a
64 MiB scratch budget, dynamic cell scheduling and deterministic integer merge.

Each scatter contribution is computed with the original float operation order
and cast to int64 **before** addition. Macro target-density scaling remains before
the cast. Gather uses untruncated weights and the original x-outer/y-inner order.
Rounded bin edges and clipping follow OpenROAD, not a reconstructed uniform
floating-point grid.

Poisson writes phi/field outputs directly into the resident bin columns. Density
is still copied into scratch because the Ooura DCT is destructive. This is one
necessary scratch copy, not a layout conversion. FFT stages remain separate
because their global transforms form synchronization boundaries.

## DREAMPlace reference

Inspected revision: `6627f3327e6cc17db7782c0b90073a498531ca3c`.

- [BasicPlace.py](https://github.com/limbo018/DREAMPlace/blob/6627f3327e6cc17db7782c0b90073a498531ca3c/dreamplace/BasicPlace.py): persistent coordinate/size/offset tensors and flat node/net-to-pin maps.
- [electric_density_map.cpp](https://github.com/limbo018/DREAMPlace/blob/6627f3327e6cc17db7782c0b90073a498531ca3c/dreamplace/ops/electric_potential/src/electric_density_map.cpp): direct-array density operators and x-major bins.
- [electric_force.cpp](https://github.com/limbo018/DREAMPlace/blob/6627f3327e6cc17db7782c0b90073a498531ca3c/dreamplace/ops/electric_potential/src/electric_force.cpp): field gather over continuous arrays.
- [weighted_average_wirelength_merged.cpp](https://github.com/limbo018/DREAMPlace/blob/6627f3327e6cc17db7782c0b90073a498531ca3c/dreamplace/ops/weighted_average_wirelength/src/weighted_average_wirelength_merged.cpp): net-parallel CSR traversal and operator-level fusion.

These are architectural references; DREAMPlace operator code is not imported.
Its scaled fixed-point deterministic accumulation differs from OpenROAD's
unscaled integer truncation. Its exp/coefficient rearrangements and optional
skipping of electrostatic energy are also not substituted. CPU/GPU scheduling
must be measured on the actual design distribution; width sorting is not assumed
to solve large macro imbalance.

## Validation

Standalone Bazel tests are under `src/gpl/test/soa`. They cover row lifetime,
metadata growth/reuse, optimizer copy/swap/remove, and an independent all-bin
oracle for density on a rectangular rounded grid including macro scaling.
Local Clang ASan/UBSan and worker GCC checks supplement full OpenROAD runs.

Acceptance must compare real placement output against the pinned native binary,
including coordinates/orientations/masters, HPWL, WNS/TNS, area and legality.
GCD/AES alone do not establish macro, IO, region or timing-repair correctness.
GPU interfaces are migrated for the new coordinate type but are not GPU-validated;
CPU results must not be advertised as GPU acceleration.

## Existing controls

`set_thread_count N` selects GPL CPU parallelism; the density scatter uses that
count subject to its private-grid scratch budget, and gather/wirelength remain
parallel over cells/nets. Keep `sta::set_thread_count 1` when isolating GPL changes.
The 64 MiB scratch cap is an implementation constant, not a new Tcl option.
The fork enables resident SoA directly; no per-call AoS-to-SoA switch is required.

Final parallel acceptance uses `OMP_PROC_BIND=FALSE`, no `OMP_PLACES`, and
`OMP_DYNAMIC=FALSE`. The earlier close/cores setting binds the main thread before
OpenROAD samples sched_getaffinity(), causing physicalCoreCount() to report one
core and silently cap requested GPL parallelism to one. Those earlier requested
8-thread case results are serial evidence only. Final runners reject any reported
placement thread cap. This is an existing runtime/configuration interaction,
not a SoA kernel limitation. Setting more threads can lose on small designs because startup,
barriers and private-grid reduction remain real work.

Measured results and exact comparison coverage: [验收报告](validation/RESULTS.md).

## Diagnosis remediation (2026-11-06)

The iteration kernels now use topology-boundary pool/row maps for coordinates,
center and pin synchronization, density and preconditioning. Bin overflow is
computed in contiguous FFT order with a separate native-order scalar reduction.
Density overlap uses separable rounded edge arrays. CPU Ooura transforms can
parallelize independent rows and columns with the configured GPL worker count.
The cached cell lock flag is refreshed after incremental phase 2 unlocks physical
instances. Optimizer membership changes go through append, swap-remove and
rebind helpers; the consistency debug channel checks the row maps and cached
flags against live instances. WA intermediate columns use net-CSR slots while
pin geometry retains stable rows, and are cleared only at topology rebuilds.
CPU density gather and preconditioning share a traversal. The measured IPO/ISA
experiment regressed and its options were removed from the final source.
See [remediation and acceptance](validation/SOA_REMEDIATION.md) for scope,
build experiments and differential results. AES measurements are specific to
that configuration; larger designs are validated separately.
