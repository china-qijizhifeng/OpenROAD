# Nesterov coordinate SIMD investigation (read-only production)

Source: current working `src/gpl/src/nesterovBase.cpp` at 2026-11-06 17:52 UTC; parent is freezing snapshot10a independently. No production edits or additional full-placement sample requested by this investigation.

## Contract established from source

- `updateNextCoordinatesColumns` (4170-4223): per mapped cell `r=cells.rows[j]`, `k=cells.indices[j]`; optimizer arrays are x/y contiguous `PointColumns`, geometry and locked byte belong to storage row r.
- Locked cell writes `(nx,ny)=(cx,cy)` and `(nsx,nsy)=(sx,sy)` exactly and skips arithmetic/bounds.
- Unlocked operation order: `x=sx+step*gx`; `slp=x+coeff*(x-cx)`. y follows the same expression. No reassociation/FMA is allowed under the existing contract.
- `hx=(ux[r]-lx[r])/2.0f`; geometry subtraction is integer first, then float division. The inside function tests original value-minus-half against lower, then original value-plus-half against upper; both tests use the original value. If both hold (oversized cell), upper assignment wins. A nested min/max clamp or mutually exclusive conditions can change this result. Ordered comparisons also preserve NaN/unordered behavior.
- `kernelRows` (2924-2949) groups by storage while appending optimizer indices in increasing original order. Row slots can be recycled/permuted after topology changes; do not assume r==k or global physical-row contiguity. `checkColumnState` verifies a permutation, not contiguity.
- `getStepLength` (3694-3738) calls `getDistance` twice; the scratch square pass (5805-5824) has independent contiguous writes. Its final `sum += x[i]; sum += y[i]` is an ordered scalar reduction. Do not change this order to make step length SIMD; iteration/backtrack trajectory depends on the resulting bits.
- `updateNextIter` locked-restore loop (4020+) separately copies 5 pairs before vector swaps, and still resolves handles/isInstance/isLocked. It is separate from coordinate arithmetic and should not be conflated with the function under study.

## Existing measured evidence (not a new profile)

AES 1T frozen original targeted sampling: 99Hz cycles:u DWARF8192; perf lost=0. `targeted-reports/{current,new}/flat.txt` coordinate outlined self is 130 samples each: 4,728,920,686 period (3.35%) current and 4,775,368,956 (3.30%) snapshot8a. `getDistance` outlined self is 16 / 15 samples (0.42% / 0.38%). These samples support a small isolated target, not a claim that coordinates dominate whole-placement time. No exact prediction of future speedup follows from the percentages.

## Minimal candidate and limits

1. First inspect final ELF outlined assembly. A packed x/y instruction alone is not cross-cell SIMD. Look for four/eight different cell lanes and load/store address progression, actual conditional masks and tail handling. Watch_build has been asked to export existing real artifacts; no independent compile is authorized here.
2. For explicit AVX2, keep the original mapped domain and process 4/8 independent cells: gather optimizer and geometry columns, ordered comparison plus blend for lower, then ordered comparison plus blend for upper, then preserve original current coordinates for locked lanes. A final output blend alone does not preserve the skip-arithmetic contract: packed operations may raise invalid/overflow on a locked lane first. Minimal exact policy: any block containing a locked lane takes the unchanged scalar per-cell fallback, so no locked lane enters SIMD arithmetic; only all-unlocked blocks vectorize. A more elaborate masked-input policy requires its own proof for nonfinite step/coefficients and extreme inputs, and fenv regression tests, before selection. Scatter outputs with scalar lane stores because AVX2 has no packed scatter. Preserve float integer-conversion/subtraction order, round points and no FMA; do not mix in horizontal reductions. Retain scalar tails and non-AVX2 fallback.
3. A stronger layout fast path can reuse contiguous optimizer-index runs (mapping metadata, no numeric shadow columns), retaining generic mapped fallback. This is optional only after real run-length evidence: topology changes can fragment row/order, and gather/scatter cost may exceed arithmetic saved.
4. Splitting into a dense arithmetic pass and indexed boundary/lock pass would enable auto-vectorization but adds memory traffic and speculatively computes locked lanes; treat as a separately measured alternative, not an assumed free simplification.

## Acceptance before production selection

- Differential oracle against unchanged scalar for all 4 output columns bitwise; locked and mixed masks; sizes 0/1/3/4/7/8/9 and tails; permuted/recycled rows and nonidentity optimizer indices; odd integer sizes, exact lower/upper equality, both-conditions-true oversized cells, signed zero, NaN/infinity where supported by scalar contract; compare FE_INVALID/FE_OVERFLOW/FE_DIVBYZERO and retain skip-arithmetic on locked lanes (mixed-lock blocks can use scalar fallback). No change to step-length serial reduction.
- Existing full placement state/QoR/iteration/backtrack/gradient trace must remain exact for AES and both large designs at actual 1/4 threads.
- Real final ELF must contain cross-cell packed arithmetic in the reachable runtime-dispatched helper, generic fallback baseline ISA, no FMA or fastmath.
- Benchmark unchanged scalar and candidate in the same E2B sandbox serially with three alternating repetitions. Require isolated coordinate-kernel benefit and no meaningful full-placement regression; current ~3.3% self samples set modest expectations. Do not add full-placement profiling just to answer this static question.

## Verified actual snapshot9 final ELF assembly

Watch_build export: `build/soa/row-storage-20261106/snapshot9-loop-disassembly-v2/receipt.json`; formal uninstrumented ELF SHA `d7ddd6f2280b9842076ce6f9c38edf40e6351945d1515c8d5621fc8968392e7f`. Each exported .s SHA matches its receipt.

- `2.s`, coordinate outlined 0x203d3b0: loop advances `inc %rax` at 0x203d4d1, then loads one row (`mov (%rdx,%rax,8),%rdi`) and one optimizer index (`mov 0x0(%r13,%rax,8),%rcx`) at 0x203d4dd. Locked branch is `cmpb`/`jne` at 0x203d4e6/4eb. Each output uses `movss` to that one `k`.
- Actual x/y step and momentum at 0x203d529-0x203d5b3 are scalar `mulss/addss/subss`. No cross-cell arithmetic SIMD.
- Existing partial SLP packs four values *within the same cell*: x/y and their momentum values are assembled by `unpcklps`/`movlhps` at 0x203d5c9-5d5, geometry half sizes and bounds use `cvtdq2ps`/`mulps`; lower and upper comparisons use `subps/cmpltps/movmskps` (0x203d5d9-5e6 and 0x203d61c-625). This is useful SSE intra-cell clipping, not four cells at once. Individual mask bits then drive a branch ladder 0x203d5ee-6ee.
- The coordinate outlined has only two calls, `__kmpc_for_static_init_8u` before the numerical loop and `__kmpc_for_static_fini` after it. There is no per-cell accessor/helper call in this kernel; calls are not its remaining obstacle. Double indirect row/index access, scalar locked choice, noncontiguous output stores and clipping branch ladder are the concrete remaining obstacles.
- `6.s`, getDistance scratch outlined, has actual cross-cell SSE4-lane arithmetic at 0x203a1f0-0x203a21e: `movups` four x, four y; `subps/mulps`; `movups` four independent scratch outputs; `add $0x4,%rsi`. There are runtime alias checks and scalar/tail fallback. This pass is already auto-vectorized; merely adding an omp simd pragma does not establish a new gain.
- `5.s`, getDistance reduction, deliberately remains ordered `addss` x0/y0/x1/y1/x2/y2/x3/y3 at 0x2039f50-0x2039f78; loop unrolls four but has one scalar accumulator. Keep this order.

The evidence supports a local AVX2 cross-cell coordinate helper or removing the existing intra-cell mask-to-branch ladder as separate experiments. It does not support a new rewrite of getDistance reduction. Snapshot10a source has no intended coordinate arithmetic change; final 10a assembly may be checked when exported, without delaying this snapshot9 conclusion.


## Actual10a cross-check (distinct binary)

Final uninstrumented 10a ELF SHA `89f9d9777b01434c9d7c1a744148a16446bd3ea4490fa61fb1010ca21e6941f1`, outlined export `snapshot10a-loop-disassembly/2.s` SHA `c50c3748170ce7c12cf09581f52662520631002b8b699b8c3f946bc14eaa07f7`. Its coordinate outlined starts 0x000000000203d820 and has the same mnemonic counts as snapshot9: 225 instructions, 4 mulss, 4 addss, 2 subss, 2 addps, 2 cmpltps, 2 movmskps, 16 movss, 2 OMP calls. This confirms the same intra-cell packing/scalar per-cell update structure in actual10a; snapshot9 addresses above describe snapshot9 only.

Diagnostic prelink remarks `snapshot10a/nesterovBase-vectorization.log` explicitly report scratch loop at5811 vector width4/interleave1, and ordered reduction5819 not vectorized. The loop4202 pragma is reported “call instruction cannot be vectorized” and the outer storage-list loop4195 not vectorized. Because final outlined machine code contains no inner numerical calls, this pragma/outer-loop remark is not evidence of a remaining per-cell helper call. Actual final code is the stronger evidence for the inner loop.
