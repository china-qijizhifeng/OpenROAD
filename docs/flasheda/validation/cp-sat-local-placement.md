# CP-SAT local placement validation

The implementation is described in [the integration and solver research note](../constraint-local-placement.md).
Tracking issue: [#5](https://github.com/china-qijizhifeng/OpenROAD/issues/5).

## Results (2026-11-07)

Validated in an isolated Linux E2B environment, using fork base
`2847efec57c70e532ca40cf7969e044541a25a94` plus this change, and OR-Tools
`9.14.6206`. The OpenROAD Python wheel was built from the fork with Bazel; its
installed Python modules were checked byte-for-byte against the tested sources.
The [machine-readable receipt](cp-sat-local-placement.json) records source,
dependency, wheel and output fingerprints. Bazel used hermetic Python 3.13;
the isolated installed-wheel check used Python 3.12.3. The task's sandbox was
released after collecting the evidence.

| Check | Result |
| --- | --- |
| `//python/openroad:placement_model_test` | Passed; 45 finite problems compared with exhaustive enumeration |
| `//python/openroad:cp_sat_placement_test` | Passed; 7 real OpenDB integration methods with subcases |
| Build `//python/openroad:openroad_wheel` | Passed |
| Import installed OpenROAD and the new module without OR-Tools | Passed; solver is optional |
| API suite from isolated installed wheel | 7 methods passed |
| Installed CLI, input ODB to new output ODB | `applied`, HPWL 10,000 → 6,000 DBU, displacement 4,000 DBU, 6 total candidates |
| Fresh Design reads CLI output and runs DPL `checkPlacement` | Passed; joint swap, fixed instances and connections verified |
| Real one-site-gap rejection | DPL rejected the candidate and the original placement was restored |

The joint-swap case first tries one movable cell alone and verifies that it stays
unchanged. Selecting both cells then permits the swap and the strict geometric
improvement. No cell master or net connection is changed.

During acceptance, two fixture calls used DEF's `FIXED` spelling instead of
OpenDB's `LOCKED`, and the pure Bazel test lacked its source import path. These
test/registration errors were corrected before the passing runs. A one-site-gap
candidate was then correctly rejected by DPL; it is now an explicit rejection
and rollback case, alongside a two-site-padding case that succeeds. The installed
wheel suite also uses an independent directory and imports OpenROAD before `odb`
to exercise the actual module-loader contract without source-file shadowing.

## Scope

The acceptance target is a real OpenROAD/OpenDB placement operation, not just a
standalone constraint model. The fixture has two movable cells in a completely
occupied row between fixed neighbors, plus asymmetric-pin cells on both row
orientations. One-at-a-time movement cannot swap the packed cells. Joint placement
can reduce their affected-net pin-centroid HPWL from 10,000 to 6,000 DBU while
preserving connectivity and passing DPL checks. These values are fixture
acceptance expectations, not full-chip PPA measurements.

Validation also covers fixed external instance pins and ports, placement padding,
mirroring, exact origin/orientation rollback, unsupported input rejection, and
ODB write/read persistence. The pure finite model has an independent enumeration
oracle covering 45 tiny problems, including actual overlap/padding infeasibility,
empty domains, empty nets, negative coordinates, multiple pins per net and fixed
boundary points.

## Reproduction

```sh
bazelisk test --test_output=errors \
  //python/openroad:placement_model_test \
  //python/openroad:cp_sat_placement_test
bazelisk build //python/openroad:openroad_wheel
```

The installed-wheel check additionally runs without the optional OR-Tools Python
package to verify that importing OpenROAD still succeeds, then installs the
pinned dependencies and executes the ODB-to-ODB CLI. The output is read by a fresh
Design and checked with `checkPlacement`; installed module bytes are compared
with the tested source so a source-directory import cannot mask packaging errors.

## Result boundaries

The model tests compare final lexicographic objectives with exhaustive enumeration.
Transaction tests inject solver UNKNOWN/UNSAT and downstream rejection to exercise
no-write/rollback behavior; these injections are not measured solver timeouts or
new physical infeasibility proofs. Geometry and placement checks use actual
OpenDB/DPL objects.

The flow has no Liberty/SDC timing acceptance or detailed routing benchmark in this
first slice. It does not establish timing, power, routability, cell-sizing or
buffer-insertion improvement. No cross-solver performance ranking is reported.
Raw build time includes dependency fetching/cache behavior and is not a placement
performance result.
