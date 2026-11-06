# Validation runners

These scripts preserve the executed worker commands and comparison criteria.
They use the existing Northjob experiment's absolute paths, parser and Tcl fixture;
they are **receipts**, not a portable public test harness. The fixture data are
upstream OpenROAD GPL tests. Standalone portable C++ tests are registered under
`//src/gpl/test/soa/...`.

- `validate-resident.py`: GCD and AES against the accepted native reference.
- `validate-cases.py`: original macro/region/cluster/timing/IO/incremental Tcl
  tests, followed by exact instance/port/net-connectivity comparison against the
  native binary at the same 1 or 8 threads. Upstream golden DEF comparison is
  replaced by this same-environment differential comparison; other Tcl assertions
  are retained. A mismatch stops the batch.
- `benchmark-resident.py`: alternating-order AES baseline/candidate measurements,
  three samples each at 1 and 8 threads, with exact coordinate/QoR checks.

- `measure-memory.py` / `measure-resources.py`: separate 8-thread runs with child-process peak RSS and CPU time.

The final source receipt includes SHA-256 hashes of every modified production
file and binary. Earlier timing/RUDY experimental overlays are removed from the
final candidate before acceptance. No source changes are made to OpenSTA.
