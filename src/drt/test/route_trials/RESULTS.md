# Recorded routing experiments

This branch is a research draft. It provides mechanisms to propose and validate
edits; it does not establish an optimizer that preserves every original quality
metric. No trial layout was adopted or combined. The initial unoptimized
checkpoint remained the baseline and fallback.

The older adaptive strategy comparison imports
[upstream PR #11868](https://github.com/The-OpenROAD-Project/OpenROAD/pull/11868)
at `7c73eec6d5f2e6294da8f5bd2c0d2f8c8d9f1ea3`, with explicit experiment switches.
Its policy-specific tests were not imported. The validation below must not be
presented as comprehensive adaptive-strategy coverage.

## Native C++ and baseline fidelity

The v9 binary, built from commit `38ab64d3740f273ab6fae09ec9a70859878e9b66`
plus the recorded working-tree source hashes, passed 21 RouteTrial and 14
ViaCleanup unit tests, and the `ndr_vias1` / `ndr_no_auto_taper` native golden
regressions. Its SHA256 is
`f29211a1b7fd54475e0749d15d9b4d8fc249cea5cf06cccb166ea50b4e918067`.
The ViaCleanup result was a valid Bazel cached result; the two NDR regressions
executed in that validation run.

Two-design original/noop controls compared DEF, date-normalized SPEF, endpoint,
quality and identity reports exactly. The NDR import fix preserves taper state
across bends/vias and resets it on path/RULE boundaries; it applies regardless
of whether the experimental optimizer is enabled.

Final PR cleanup uses the repository CI formatter. This is separately recorded
from the historical v9 build: include sorting, adjacent string-literal splitting,
comment reflow and whitespace normalization are not a new post-format C++ build.
Exact tested-source hashes and binary identities are in [validation.json](validation.json).

The portable package passed three Bazel targets containing 26 synthetic tests
in an isolated Bazel 9.2.0 workspace using rules_python 2.3.2, Python 3.13 and
the repository Python configuration. All ten historical selector snapshots
replayed identically, and the historical proposal produced the same complete
adapter transaction. This was not a whole-repository build.

A separate smoke run of the portable Tcl helper on the historical v9 binary
passed four controls: sourcing without database changes, exact three-net
roundtrip fidelity, a known proposal modifying three nets while preserving
untouched-net tokens, and a deliberate error in the second net. The failure
control confirmed the first net had already changed, exited without exporting
a database, and preserved original checkpoint hashes. This demonstrates
private-process isolation, not live rollback. The three successful controls
had zero fresh DRC and antenna violations and reported routed status. No new
search or electrical-quality acceptance was performed in this smoke run.

## Structure-guided single-net pilot

Full-design features covered tinyRocket/Nangate45, Ibex/Sky130HD and
JPEG/Sky130HD; native candidates covered tinyRocket and Ibex. Native input
stages matched the original post-antenna-repair checkpoint. JPEG's historical
original already had one antenna violation; it was profiled only.

Twenty-nine selected sources produced 22 native graphs. Five source-contact
rejections and two node-budget rejections were not timing failures or proofs of
physical impossibility. Three candidates survived native physical checks; all
three failed the strict endpoint gate. No candidate was adopted.

Direct-fanout and transitive-cone selection found smaller apparent impact, but
also less geometric opportunity. Exact audits found 14 unproductive chains and
four three-terminal trees infeasible under the via/length budgets **in those
fixed graphs**. This did not establish that annealing was trapped, or that a
larger window, other layer allocation or joint-net move was impossible.

The measured lesson was to screen geometric opportunity before spending search
budget, then use circuit and electrical structure to estimate impact. A small
logical cone is neither a guarantee of physical feasibility nor a complete
coupling/timing influence set.

## Multi-net import and feedback comparison

The first independent multi-net interface run evaluated 12 three-net proposals
on the v7 binary (`c06955ea92a17500cc0b7d1d85f7ad3c20aae387934d9a14d097937587a80c1f`).
All passed physical checks and removed two vias, but each regressed 46 setup
and two or three hold endpoints. Private-process roundtrip and deliberate
second-net-error rollback controls passed. This is checkpoint-level isolation,
not a shared C++ live-database multi-net transaction.

The subsequent feedback pilot reused that v7 binary and compared a fixed
geometry/diversity control with mixed timing/physical feedback. It used a new
26-candidate pool in the original region and a 30-candidate pool in an independent
two-net region of the same design. All twelve historical layouts were excluded;
their timing results did not train either arm.

Before any new native timing result, eight independently ranked regions were
screened with equal geometric budgets. Candidate counts were
`0 / 10 / 30 / 102 / 31 / 0 / 0 / 54`; the first region with at least twelve
candidates was chosen. This pre-native protocol amendment is recorded as
exploratory work, not a pristine held-out confirmation trial.

Each arm had fourteen logical evaluations, including four common seeds.
Twenty-four requests reduced to nineteen actual layouts through five cross-arm
cache hits. All choices and training snapshots were saved before native
execution, then independently replayed. No control-only or future labels leaked
into feedback selection.

| Region | Arm | Logical evaluations | Physical pass | All endpoint slacks nonworse | Full PASS |
| --- | --- | ---: | ---: | ---: | ---: |
| Original region | Feedback | 7 | 6 | 0 | 0 |
| Original region | Geometry/diversity | 7 | 5 | 0 | 0 |
| Independent region | Feedback | 7 | 2 | 2 | 0 |
| Independent region | Geometry/diversity | 7 | 1 | 1 | 0 |

Counts include shared seeds and cannot be summed as independent layouts.
Among nineteen actual layouts, ten passed physical checks and nine were rejected
for spacing. The frozen comparator reported eight FAIL, two INCOMPLETE and zero
PASS. Default power activity did not converge for original or candidate; its
numeric deltas were not interpreted as reliable power gains or regressions.

Both endpoint-nonworse layouts removed two vias and 760 DBU (0.38 micrometers).
Each of the 4061 finite setup and 4061 finite hold endpoints was exactly equal
to the original, as was the complete exported endpoint file. The first such
layout appeared at the same regional evaluation budget in both arms, during
exploration before usable timing labels. The extra feedback layout was selected
with only one usable timing sample: every predicted timing loss was zero, so
physical-failure estimates and geometry ranking still distinguished candidates.
This does not isolate or prove the benefit of STA labels alone.

In the original region, feedback found a lower hold-loss candidate at the cost
of slightly greater setup loss. The compared actual layouts were:

| Representative | Setup maximum loss (ps) | Hold maximum loss (ps) | Setup/hold regressed endpoints |
| --- | ---: | ---: | ---: |
| Control's shared seed | 0.045537949 | 0.011429191 | 1 / 5 |
| Feedback choice | 0.048995018 | 0.001937151 | 1 / 4 |

Both failed the original zero-regression gate. Weighted or normalized search
scores never compensated for a failed acceptance metric. The geometric KNN also
made substantial errors: one round predicted one setup regression and produced
34; another returned to 46. No model or pool was retuned after observing results.

## Why equal endpoint slack was insufficient

Additional through-net reports exposed setup regressions that the frozen
endpoint-worst-slack metric hid. The same startpoint, full pin/net sequence,
clock paths and required time were preserved in every paired comparison:

| Layout | Setup slack delta through first edited net (ps) | Through second edited net (ps) |
| --- | ---: | ---: |
| First endpoint-nonworse layout | -0.747957252 | -0.450528503 |
| Second endpoint-nonworse layout | -0.739186490 | -0.433209024 |

These branches still had about 246–353 ps more setup margin than the endpoint's
current worst branch, so its reported worst slack did not change. The affected
endpoint was constrained and had finite min/max slack; this was not an
unconstrained-net explanation. Minimum-delay through-net paths became slower,
improving their hold margin. At least fourteen nets changed SPEF aggregates.

The frozen statistics remain INCOMPLETE under their original metric definition;
the posthoc path regression is reported separately, not silently inserted into
the experiment's acceptance statistics. Nevertheless, endpoint equality cannot
be advertised as all-path nonregression, and supplying only missing power
results would not establish that stronger claim.

In the original region, remaining regressed paths did not traverse the edited
nets. Changed loads on unedited neighboring nets affected launch/propagation
behavior. Circuit structure should therefore guide evaluation of paths through
both edited nets and electrically affected neighbors, rather than only the
edited net's logical cone or endpoint worst slack.

## Evidence and limits

The native candidate phases totaled 432.595 seconds, serially with four threads
on one 8-CPU/16-GiB guest. This excludes preparation, generation, controls,
transport and posthoc diagnosis. No speedup against full rerouting is claimed.

The recovered archive is 676,294,242 bytes with SHA256
`68641fba5d3b894f67b3d6a3ccfec1715be2bad169918ebd291f38df2d0c8ce7`.
All 814 manifest files and 71 distinct observation evidence files passed an
independent byte-level audit. Large DB/DEF/SPEF archives, machine-specific
controllers and credentials are intentionally outside Git. The compact
[validation manifest](validation.json) records the identities and outcomes.
All experiment sandboxes were released and original checkpoints retained.

The repository package contains portable preparation, import and selection
helpers plus synthetic tests. It does not include a portable end-to-end
production flow, a new multi-design/corner acceptance suite, complete coupling
model, or a calibrated timing surrogate. Broader default-flow/performance,
adaptive-policy and distributed-mode coverage remains required before this
research branch should be considered ready for production.
