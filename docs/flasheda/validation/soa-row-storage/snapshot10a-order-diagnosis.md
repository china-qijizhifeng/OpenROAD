# Snapshot10a whole-placement regression checks

No extra sampling or workload was run for this analysis. Existing full36 formal results are read in fixed dispatch order; large02/1T/r2/10a is retained as an excluded polluted sample, with one authorized later replacement.

## What the existing numbers establish

- large01 1T candidate is slower in all three matched pairs: +3.631082s/+4.0276% (9→10a), +0.586456s/+0.6537% (10a→9), +2.197202s/+2.4828% (9→10a). Therefore simply blaming later-in-time execution is unsupported. Candidate first repeat is unusually slow, but removing it alone would not establish a benefit; do not drop samples post hoc.
- Whole-process user CPU rises +4.21/+0.90/+2.90s and sys changes -0.01/-0.02/-0.01s. All six1T runs have (user+sys)/wall0.9996–0.9997. The evidence shows CPU execution time growth, not additional measured I/O/sys wait or an idle scheduler gap. CPU frequency, memory stalls and codegen can all change user CPU time; no periodic frequency/steal/cache trace was recorded, so none is identified as root cause.
- Whole-process nonplacement time also increases +0.568918/+0.293544/+0.692798s. The placement-only slowdown is not the only elapsed difference. Since this region is not phase/profiled in the current formal round, assigning it to WA/cache is unsupported.
- large01 4T matched placement changes are -0.176%/-5.723%/+4.588%; ranges overlap and the median effect is small. AES4T all matched pairs are slightly faster while AES1T changes are mixed. large02 original4T changes are -3.271%/-4.804%/-0.155%; full1T final statistics await/require the single replacement.

## Verified static equivalence and concrete candidates

44 current live source hashes equal the10a final source receipt and formal ELF receipt. Of42 common recorded files from9,40 are byte-identical, including nesterovBase.cpp/.h, coordinate, density, initialPlace, SoA storage, and optimizer columns. The only changed common files are wirelengthSoa.cpp and the exact SIMD unit test; newly recorded BUILD/Simd.h supply the helper and build dependency. Fixed compiler SHA, O3/ThinLTO/C++20/ffp-contract-off, no global AVX2/FMA/fastmath, same Tcl/input SHA and same state/algorithm trace are checked. Formal minimal environments exclude observation/replay controls. There is no evidence of a changed iteration/backtrack algorithm.

The9 WA force outline inlines scalar derivative arithmetic.10a actual final outline has one numerical call to the target-attributed `scatterNetPinGradientsAvx2` for eligible per-net spans>=4; it loads/repackages irregular geometry, masks four flags/direction, stores lane values and handles scalar tails. These costs can dominate small eligible nets and cannot be amortized from the mere presence of AVX2 opcodes. In the existing real large02 replay only41.45%pins are in4pin groups and58.55%tail. No captured large01 per-net degree/flag/consumer distribution exists in this experiment, so its exact SIMD coverage or call-cost share is not known.

D/E replay compares identical actual production SoA stores/helper with a warmed, repeated single frozen numeric snapshot and excludes decode/other placement. Its local improvement is compatible with whole-placement regression when dispatch/cache/working-set/adjacent kernels differ. This is an explanation of why both observations can coexist, not a proof of cache cause. Preserve both observations.

## Falsifiable next decision

Do not change the4pin threshold or cross-net batching based on the local D/E results alone. The current experiment has not proved10a a whole-placement winner on all designs; large01 1T is a measured regression. A later controlled experiment can compare exactly the same10a ELF path with forced scalar versus AVX2 only, and real per-net4/8/16+ classes, then decide a threshold/default. That would separate target dispatch/tail/mask cost from cross-ELF codegen or host drift. It requires a frozen instrumentation/selection design and separate formal timing; no such extra experiment is run or claimed here. Keep9 as the known comparison and expose SIMD as optional only if parent chooses that product tradeoff.

Machine-readable paired wall/user/sys/RSS/order evidence: `summary/snapshot10a-pair-order-resources.json`. All raw results remain immutable; polluted sample exclusion and exact replacement are independently recorded.
