# Experimental route trials and joint checkpoint edits

These tools support bounded routing experiments, not a production via optimizer.
No layout from the recorded experiments was adopted. Fewer vias, shorter wires,
clean DRC, and unchanged endpoint worst slack do not establish complete timing
or power nonregression. See [RESULTS.md](RESULTS.md) for the measured limits.

## Interfaces and isolation

| Interface | Scope | What success means |
| --- | --- | --- |
| `DRT_VIA_CLEANUP=1` | Conservative seam-adjacent via cleanup after a clean fresh local route | Local geometry/GC and connectivity checks; no RCX/STA acceptance |
| `drt::experimental_route_trial_cmd recipe` | One net; a two-port chain or a 2–8-port tree | Requested replacement passed local checks and connectivity/writeback |
| `drt::experimental_route_graph_cmd source.recipe graph.txt` | Bounded single-net search graph | Geometry-filtered graph exported; edges are not full DRC proofs |
| `joint_adapter.py` + `route_trials::import_joint` | Up to 16 nets, 2–64 boundary ports per edited net | Prepared proposal and private-process import; no complete QoR acceptance |
| `feedback.py` + `select_feedback.py` | Frozen-pool sequential research selection | A reproducible next choice, never permission to adopt a layout |

Use a fresh private OpenROAD process loaded from an immutable original checkpoint
for every trial. A failed process is discarded. The C++ local edit can roll back
its local geometry/GC trial, but subsequent global connectivity processing and
writeback do not provide a whole-design live-database transaction. The multi-net
Tcl helper likewise has no live-database rollback and must not be caught and
continued on a shared working database.

The optimization switches are off unless explicitly requested. This does not
mean the entire branch has no default effects: the NDR taper-state import fix is
unconditional, and the imported incremental probe still increments maze-search
counters when reporting is disabled. No default-performance claim is made.

## Existing routing experiment switches

| Switch | Behavior |
| --- | --- |
| `detailed_route -adaptive_strategy` | Enable the imported adaptive search-and-repair strategy |
| `DRT_ADAPTIVE` | Presence enables the same strategy; unset it to disable the environment override |
| `DRT_ADAPTIVE_FORCE_UNIDIR` | Presence selects the unidirectional-policy branch during experiments |
| `DRT_OR_SEED`, `DRT_OR_K` | Override net-order seed and perturbation strength; K defaults to 0.1 when the seed override is set |
| `DRT_ITER_MODE` | Override selected ALL iterations, e.g. `1:NEARDRC,2:DRC` |
| `DRT_NEARDRC_BLOAT_MULT` | Marker-neighborhood multiplier; default 2 |
| `DRT_NETCHG` | Append incremental-routing probe records to the supplied filename |
| `DRT_VIA_CLEANUP` | Cleanup runs only when the value is exactly `1` |

The presence-based switches remain enabled if their value is `0`. Several probe
settings are cached in process-local statics; start a new process when changing
experiment settings. The adaptive/probe commits are included for research
continuity and are not a claim of improved QoR.

## Single-net recipes and graph export

The first line is the exact database net name, including any DEF escaping.
Remaining records use database units and actual technology layer/via names:

```text
<net-name>
topology tree
anchor <routing-layer> <x> <y>
remove_segment <routing-layer> <x1> <y1> <x2> <y2>
remove_via <via-definition> <x> <y>
add_segment <routing-layer> <x1> <y1> <x2> <y2>
add_via <via-definition> <x> <y>
```

`topology tree` is optional; omit it for a two-port chain. Supply every boundary
port and exact whole source objects. Tree sources are capped at 64 segments and
32 vias; replacement segments are capped at 64. The replacement must preserve
connected, acyclic topology, strictly reduce vias, and not increase centerline
length. Internal ports and crossings are normalized before topology checks.
Clock, fixed, special, NDR and unsupported routing configurations are rejected.
A saved regular GCELLGRID is required; MAXSPACING technology is unsupported.

A source-only recipe can be exported with:

```tcl
# Run in a private process after loading the original checkpoint and technology.
set ok [drt::experimental_route_graph_cmd source.recipe graph.txt]
```

The graph covers source geometry plus a two-pitch halo and has a 40,000-node
budget. Its first record is `graph_v1 old_length old_vias start end`, or
`tree_graph_v1 old_length old_vias terminal_count terminal_ids...`. Subsequent
records describe nodes and candidate edges. A graph rejection or absence of a
path in this fixed graph does not prove global routing impossibility.

## Joint proposals

`joint_adapter.py` uses the Python standard library and accepts explicit paths:

```sh
python3 src/drt/test/route_trials/joint_adapter.py   original.def tech.lef prepared --proposal proposal.json
```

The proposal schema is `joint_route_proposal_v1`. It binds the original DEF SHA256
and a canonical JSON content hash. Each edit has `net`, `ports`,
`remove_segments`, `remove_vias`, `add_atomic_segments`, and `add_vias`.
Segment coordinates are `[layer_number, [x1,y1], [x2,y2]]`; a via is
`[via_name, [x,y]]`; ports are `[layer_number,x,y]`. Hash the complete proposal
object except `proposal_sha256` and `acceptance` with `joint_adapter.digest`.

The adapter recomputes boundary ports from retained geometry, verifies whole
object removal, tracks, topology and the original identity, then writes
`patch.def` and `transaction.json` into a previously nonexistent directory.
Every edit is checked before output is written. Vias must decrease across the
whole proposal, and total wire length must not increase; an individual net may
become longer. Output status is `PREPARED_NOT_ACCEPTED`.

The parser intentionally supports only `metalN` layers, supported explicit LEF
via definitions, a rectangular die, and a limited regular-net DEF grammar. It
rejects clock/NDR/fixed/cover/taper/unsupported syntax for edited nets. It is not
a general DEF reader and does not model complete spacing, pins or obstructions.
Native physical checks remain necessary.

In a private process with libraries, the original DB and SDC already loaded:

```tcl
source src/drt/test/route_trials/import_joint.tcl
route_trials::import_joint prepared/patch.def
# The caller checks fidelity, fresh DRC, antenna and routed status before
# extracting parasitics and evaluating the complete joint candidate.
```

`read_def -incremental` skips NETS in the tested OpenDB version. The helper uses
`-floorplan_initialize` for an adapter-produced partial DEF containing only the
original header and existing NETS. It verifies net and terminal membership,
restores ITerm ordering using `disconnectDbNet`/`connect` (preserving modnets),
and invokes the experimental native connectivity pass. It does not save or
publish a DB. Never pass an arbitrary floorplan DEF to this helper.

Use `joint_adapter.audit(original_def, candidate_def, tech_lef, transaction_json)`
to check the exported candidate against the exact proposal and unchanged nets.
Fidelity is separate from physical legality and electrical quality. A roundtrip
control can be prepared with `--roundtrip NET_A NET_B`; run it and an unmodified
control through the same downstream flow before testing changes.

## Feedback selection

`feedback.endpoint_feedback` compares identical endpoint sets, rejects changed
applicability/nonfinite values, and returns separate setup/hold positive losses
(maximum, total and count). Improvements never cancel regressions. Missing
measurements remain missing, not zero. This endpoint feedback does not include
every path through a net; the measured branch-masking example in RESULTS.md is
an explicit limitation.

```sh
python3 src/drt/test/route_trials/select_feedback.py freeze config.json plan.json
python3 src/drt/test/route_trials/select_feedback.py next   plan.json observations.json REGION ROUND choice.json
```

A config lists one or two regions as `{"regions":[{"region":"example","pool":"pool.json"}]}`.
Pool paths resolve relative to the config; proposal paths resolve relative to the
pool. A pool contains canonical geometry IDs, route-edge sets, geometric feature
vectors, via/wire metrics and proposal hashes. The selector is not a candidate
generator or a physical router.

Freeze the pool, model and complete control sequence before native outcomes.
The plan gives four shared observations and ten later choices per arm (split
equally for two regions): fourteen logical evaluations per arm, at most 24
new distinct layouts. Native callers must enforce that budget and cache exact
geometry across arms. The helper alone does not execute or authenticate STA.

The feedback arm reads only shared and its own earlier observations in the
same region. It uses up to three neighbors, combining route-edge Jaccard and
standardized geometric distance. Every third round explores geometrically;
other rounds use independent setup/hold predictions, a Pareto filter and a
separately normalized minimax ranking. A physical-failure estimate also affects
selection: this is mixed feedback, not a clean ablation of timing labels alone.
It is not a calibrated probability of safety. Save the full prediction receipt
before evaluating the selected layout, then independently replay it for audit.

## Validation and acceptance

```sh
bazelisk test //src/drt/test:route_trial_unittest   //src/drt/test:via_cleanup_unittest   //src/drt/test:ndr_vias1-tcl_test   //src/drt/test:ndr_no_auto_taper-tcl_test
bazelisk test //src/drt/test/route_trials:all
```

The portable Python tests use small synthetic DEF/LEF fixtures; they do not
replace native import/rollback, physical, RCX or timing validation. Historical
native versions and evidence hashes are recorded in [validation.json](validation.json).
The repository helpers are extracted from those experiments; their portability
checks are recorded separately from the original native runs. The extracted
Tcl helper also passed source-only, three-net roundtrip, known-proposal and
deliberate later-net-failure smoke checks on the historical v9 binary. This
checks the portable helper; it is not a rebuild of the final formatted C++.

Keep the initial unoptimized checkpoint as the comparison baseline and fallback.
Check the complete joint layout, not a sum of separately accepted edits. Reject
regressions and distinguish missing evidence from a pass. Unconverged default
power activity is not evidence of power nonregression. No helper in this folder
provides automatic adoption, and no current result establishes all-path or
cross-corner nonregression.
