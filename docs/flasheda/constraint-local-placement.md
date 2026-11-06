# Local placement with CP-SAT and spatial constraint solver research

Tracking issue: [#5](https://github.com/china-qijizhifeng/OpenROAD/issues/5).

## Implemented slice

`openroad.cp_sat_placement.optimize_local_placement` operates on a loaded OpenROAD
`Design`. It jointly chooses horizontal positions and legal horizontal mirrors
for explicitly selected single-height standard cells. Every selected cell stays
on its original contiguous horizontal row. All other instances, ports, masters,
net connections, placement statuses and rows remain unchanged.

This is a pre-route local placement operation, intended after legalization and
before pin access. It is not a standalone JSON demonstration or a replacement
for the global placer. The optional OR-Tools Python dependency is imported only when solving;
the default OpenROAD flow and wheel do not require OR-Tools Python.

The finite model contains:

- actual row sites and bounds, cell widths/heights and existing placement padding;
- immutable neighbors, including movable cells outside the explicit selection;
- conservative exclusion of placement blockage rectangles;
- pairwise non-overlap of the selected cells' padded rectangles;
- position-dependent pin centroids and orientation-dependent pin locations;
- all affected non-supply nets, including fixed endpoints outside the selection
  and placed block ports (using OpenDB's first physical pin location);
- a displacement limit, cell-count limit, candidate-count limit and solve timeout.

The solver uses a finite integer candidate index per cell, `Element` constraints
for coordinates and pin positions, and the native `NoOverlap2D` constraint for
selected padded rectangles. HPWL extrema use exact `MinEquality`/`MaxEquality`.
Lexicographic optimization runs three solves under a shared time budget, fixing
an objective only after CP-SAT returns OPTIMAL. The solver uses one worker and a
fixed seed. FEASIBLE is reported as UNKNOWN for the adapter's acceptance contract;
this first version does not apply a partially optimized incumbent.

CP-SAT was selected for its finite-domain and native rectangle constraints. The
SMT/OMT alternatives remain part of the research and future backend work. The
shared finite model and independent checker do not depend on OR-Tools objects.

Each connected ITerm uses the average of its transformed shape centers, with
the same truncation toward zero as `dbITerm::getAvgXY`. Averaging is performed
before rounding, including for mirrored pins. Thus this metric is explicitly
**affected-net pin-centroid HPWL**, not cell-center distance, routed wirelength
or a timing metric. Fixed endpoint extrema can be compressed to two corners
without changing HPWL. Unaffected nets contribute a constant and are omitted.

The lexicographic objective minimizes HPWL, total displacement, then mirror
count. A candidate is applied only when the finite optimization has completed
with OPTIMAL status for each objective, its decoded model passes an independent Python
overlap/objective check, and HPWL strictly improves. OpenDB computes HPWL again
after the update, and real `Opendp.checkPlacement` checks the resulting design.
If placement validation fails, original origins and orientations are restored.
Unexpected exceptions also restore placement and propagate. No model is applied
for UNSAT, UNKNOWN or an incomplete optimization result.

The baseline must already pass `checkPlacement`. Fixed/dont-touch/unplaced,
multi-height, rotated or region/group-constrained selected instances are
rejected. Affected signal nets with wires, guides or extracted cap nodes, and
selected pins with preferred access points, are rejected. Supply nets are not
optimized. The caller must use this before routing/pin access and refresh
placement parasitics before subsequent timing analysis; ODB guards cannot detect
all external tool caches or analysis data loaded only into STA.

Padding and blockages are modeled conservatively. Advanced placement rules such
as one-site gaps, edge spacing and pin/obstruction interactions are checked by DPL at acceptance,
not encoded completely in CP-SAT. A rejected candidate is not automatically refined
and re-solved in this version. UNSAT means only that this bounded finite model
has no solution, not that the chip cannot be placed. A timeout means UNKNOWN.

Sizing, buffer insertion, cross-row movement, clock synthesis, routing choices,
timing/power objectives and automatic window selection are follow-ups in #5.
An HPWL improvement is not a claim of routability, setup/hold, power or signoff
improvement. Existing connectivity is preserved by construction; this operation
does not rewrite logic.

## Use in the actual placement flow

Build the fork's Python wheel and install the optional dependency in its venv:

```sh
bazelisk build //python/openroad:openroad_wheel
python3.13 -m venv .venv-cp-sat
. .venv-cp-sat/bin/activate
python -m pip install bazel-bin/python/openroad/*.whl
python -m pip install --require-hashes -r bazel/requirements_cp_sat_lock.txt
```

Call the function in the same Python flow that owns placement/padding state:

```python
from dataclasses import asdict
from openroad import Design, Tech
from openroad.cp_sat_placement import optimize_local_placement

tech = Tech()  # Keep Tech alive while its Design is in use.
design = Design(tech)
design.readDb("post_legalization.odb")
# Reapply flow-specific DPL padding/settings when restoring an ODB checkpoint.
# For example: design.getOpendp().setPaddingGlobal(1, 1)
result = optimize_local_placement(
    design,
    ["u_path/g1", "u_path/g2"],  # Replace with actual instance names.
    max_displacement_sites=4,
    timeout_ms=10_000,
    problem_path="local-placement.json",
)
print(asdict(result))
if result.status in ("applied", "unchanged"):
    design.writeDb("post_cp_sat.odb")
# In a flow with loaded timing libraries/SDC and configured wire RC:
# design.evalTclString("estimate_parasitics -placement")
```

For a checkpoint with no additional in-memory DPL settings, a CLI is also
available. It never overwrites its input checkpoint, and writes no output ODB on
UNKNOWN, UNSAT or rejected placement:

```sh
python -m openroad.cp_sat_placement post_legalization.odb post_cp_sat.odb \
  --instances u_path/g1 u_path/g2 --max-displacement-sites 4 \
  --timeout-ms 10000 --problem local-placement.json
```

ODB does not preserve every flow setting. Use the Python API and restore the
original DPL settings when the flow uses per-instance/master padding or other
configuration. The output must continue through the normal parasitic, timing,
routing and physical validation flow.

## Which solvers understand spatial constraints?

For this problem, "spatial SMT" means an encoding choice rather than a standard
native VLSI geometry theory. A rectangle can be represented by integer DBU
coordinates and dimensions. Non-overlap is a disjunction:

```text
x_i + w_i <= x_j OR x_j + w_j <= x_i
OR y_i + h_i <= y_j OR y_j + h_j <= y_i
```

With fixed dimensions, individual separation inequalities are difference
constraints. Site-index coordinates, selectable widths, HPWL sums and conditional
decisions generally move the whole problem beyond pure difference logic.
Unquantified linear integer arithmetic plus Boolean choices is a natural first
encoding. Avoid nonlinear transistor/electrical models in this first spatial
layer. Integer DBU avoids floating-point tolerances at legal boundaries.

| Solver / approach | Useful capability | Boundary and recommendation |
| --- | --- | --- |
| Z3 / OMT | Boolean choices, integer arithmetic, conditional expressions, lexicographic objectives | Alternative research backend, not implemented in this PR. Solve bounded local problems; record UNKNOWN and bounds honestly. |
| cvc5 | Integer/bit-vector theories, incremental satisfiability and proof facilities | Useful independent satisfiability backend. Use explicit bounded objective queries if needed; do not assume compatibility with Z3's Optimize API. Not benchmarked in this PR. |
| OptiMathSAT | Incremental multi-objective OMT over linear arithmetic and other theories | Relevant comparison for arithmetic-heavy objectives. Check its current MathSAT license terms before distributing an integration. Not integrated here. |
| OR-Tools CP-SAT | Finite-domain integer choices, interval and `NoOverlap2D` constraints | Implemented backend. Joint candidate indices use `Element`; padded cell rectangles use `NoOverlap2D`. Native global constraints are a good structural fit, but speed advantages require measurement. |
| MILP | Binary choices, linear costs, resource budgets and disjunctive formulations | Suitable alternative for the master/resource problem. Big-M bounds must be valid and tight; numerical feasibility tolerances require physical validation. |
| dReal | Nonlinear real constraints with delta-complete reasoning | Relevant to approximate continuous electrical models, not the default exact DBU legality engine. Delta-SAT is not exact satisfiability of the original rules. |
| SAT / FlashEDA | Finite Boolean selection, conflict clauses, bounded pseudo-Boolean objectives | Planned backend for this finite candidate model. Requires an encoder and an optimization controller, not just replacing `CpSolver.solve()`. |

Spatial *indexing* (R-trees, interval trees) is complementary: it finds which
objects can interact, reducing candidate-generation and conflict construction
cost. It does not solve the joint choices. This first bounded implementation
uses straightforward geometric filtering; it is not a scalable full-chip spatial
index or a claim of a new geometry decision procedure.

## Research precedents and scope

- [SP&R, ASP-DAC 2020](https://doi.org/10.1109/ASP-DAC47756.2020.9045729)
  jointly solves transistor placement and routing inside a standard cell using
  OMT and dynamic pin allocation. Its granularity is not a full OpenROAD block.
- [SMTCell](https://github.com/ckchengucsd/SMTCellUCSD) provides a concrete Z3
  encoding and layout conversion flow, plus pre-partitioning to control scale.
  Its documented solver version is not the version required by this new API.
- [ROAD, ISPD 2019](https://doi.org/10.1145/3299902.3309752) uses SAT to diagnose
  local routing conflicts. It motivates explanation-based feedback, not treating
  a heuristic router's failure as an UNSAT proof.
- [SO3-Cell, ICCAD 2025](https://doi.org/10.1109/ICCAD66269.2025.11240677)
  jointly optimizes topology, transistor placement and internal routing with
  **MILP**, not SMT. It is a useful formulation comparison at cell granularity.
- [Simultaneous gate sizing and placement, TCAD 2000](https://doi.org/10.1109/43.828549)
  focuses joint optimization on critical paths and neighboring cells. This
  motivates bounded windows rather than a monolithic full-chip formula.

Primary solver documentation consulted on 2026-11-07:
[Z3 arithmetic](https://microsoft.github.io/z3guide/docs/theories/Arithmetic/),
[Z3 optimization](https://microsoft.github.io/z3guide/docs/optimization/intro/),
[cvc5 bit-vectors](https://cvc5.github.io/docs/latest/examples/bitvectors.html),
[OptiMathSAT](https://optimathsat.disi.unitn.it/),
[CP-SAT model API source](https://github.com/google/or-tools/blob/stable/ortools/sat/python/cp_model.py),
[dReal](https://dreal.github.io/).
These references establish capabilities and research scope; this PR does not
claim a measured Z3-versus-cvc5/CP-SAT/FlashEDA performance ranking.

## Hierarchy and future refinement

The proposed hierarchy is:

1. Existing OpenROAD algorithms choose global placement and promising regions.
2. A bounded local problem coordinates strongly coupled placement decisions.
3. Detailed physical and timing tools validate the resulting concrete state.

The local boundary must retain fixed neighbors, ports, row/site rules, padding
and geometric resource ownership. Adding sizing requires compatible cell
functions and pin mappings, upstream load/slew effects and downstream required
times. Adding buffers also requires connectivity/equivalence checks and updated
parasitics. A timing budget shared by several windows must not be spent repeatedly.

A later logic-based Benders loop can return proven conflict cuts to a master
problem. An abstraction-refinement loop can add missing geometric detail where
the coarse model admitted a spurious candidate. These mechanisms are design
directions, not implemented solver claims. Cuts must include their boundary
assumptions, and cannot outlive a change to those assumptions without revalidation.

## FlashEDA SAT boundary

The exported `openroad.local-placement.v1` JSON contains integer geometry,
ordered cells, ordered candidates, per-candidate pin coordinates, per-cell net
identities and fixed endpoint extrema. The export is a reproducibility artifact;
there is no external-solver execution hook or imported-model apply API in this PR.

The future SAT backend should use the same finite choices:

1. Give every `(cell_index, candidate_index)` a stable Boolean ID and impose
   exactly-one selection per cell.
2. Encode incompatible padded rectangles as conflict clauses. Static obstacles
   have already filtered candidates, so retain the boundary fingerprint too.
3. Encode bounded integer HPWL, displacement and mirror objectives with proven
   ranges. Integer coordinates may be negative; offset them or use explicitly
   signed encodings. Do not permit overflow or silently change comparison semantics.
4. Export DIMACS plus a mapping/manifest. Prefer shared Boolean/AIG structure
   when the selected FlashEDA runner can consume it; plain CNF remains the
   baseline interchange format.
5. Run the pinned FlashEDA frontend/portfolio and record actual winning backend,
   result, limits, model and proof artifacts. Do not attribute a portfolio result
   to a particular internal solver without evidence.
6. Check SAT assignments against every CNF clause, decode choices, then call the
   independent finite-model evaluator and the same OpenDB acceptance path.
7. Implement lexicographic optimization as bounded feasibility queries. Tighten
   bounds only after trusted/verified SAT or UNSAT results; UNKNOWN neither proves
   infeasibility nor optimality. Fix an objective only after its optimum is proved.

This separation allows future FlashEDA experiments without teaching the SAT
engine OpenDB, Liberty or SMT arithmetic. No FlashEDA speedup is claimed here.
Encoding, every objective query, validation and process overhead must be included
in end-to-end comparisons.

## Validation and acceptance

Registered tests:

```sh
bazelisk test --test_output=errors \
  //python/openroad:placement_model_test \
  //python/openroad:cp_sat_placement_test
```

The finite model is checked against independent exhaustive enumeration on tiny
cases. Integration tests use real OpenROAD/OpenDB with asymmetric physical pins,
fixed neighbors and block ports, exercise joint swaps/mirroring and rollback,
and validate persisted placement. A solver timeout status is simulated only to test the
adapter's no-write contract; physical DB objects and placement validation are real.

Execution versions, source fingerprints and results are recorded separately in
the validation report accompanying this change. Tiny HPWL examples establish
functionality, not chip-scale runtime or PPA improvements. A later comparison
should fix designs, constraints, tool versions and compute budgets, then compare
ordinary DPL, local CP-SAT and an SMT/SAT alternative using final routing and STA.
