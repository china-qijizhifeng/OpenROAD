# Final large-design validation index

84 final measurement records, seven configurations, 12 records per configuration. Every configuration has one worker and one exact final state.

| Configuration | Worker identity | Records | Exact state |
|---|---|---:|---|
| ibex_core-td-rd | job-388e9a7b6ce7-20261106113510 | 12 | yes |
| medium04 | job-388e9a7b6ce7-20261106113510 | 12 | yes |
| macro03 | job-388e9a7b6ce7-20261106113510 | 12 | yes |
| macro03-rd | job-388e9a7b6ce7-20261106113510 | 12 | yes |
| large01 | job-bb6ab222440c-20261106125349 | 12 | yes |
| large02 | job-bb6ab222440c-20261106125349 | 12 | yes |
| large01-rd | job-bb6ab222440c-20261106125349 | 12 | yes |

Native sample 0 for the original four configurations is an earlier-window sample from the same worker, verified against raw log SHA and parsed state. The high worker executes all 36 samples fresh. Timing samples are never pooled across workers.

| Configuration | GPL threads | Native median seconds [min, max] | v3 median seconds [min, max] | Reduction |
|---|---:|---:|---:|---:|
| ibex_core-td-rd | 1 | 175.590 [172.270, 176.125] | 130.857 [130.706, 131.060] | 25.5% |
| ibex_core-td-rd | 8 | 92.610 [91.218, 97.427] | 61.482 [61.344, 62.515] | 33.6% |
| medium04 | 1 | 221.268 [220.708, 223.214] | 149.112 [149.105, 150.757] | 32.6% |
| medium04 | 8 | 87.042 [84.673, 90.042] | 46.403 [45.169, 47.421] | 46.7% |
| macro03 | 1 | 189.959 [189.556, 192.995] | 135.379 [135.088, 135.940] | 28.7% |
| macro03 | 8 | 72.266 [70.912, 74.659] | 40.494 [39.952, 41.651] | 44.0% |
| macro03-rd | 1 | 456.460 [454.330, 460.040] | 322.789 [322.226, 322.888] | 29.3% |
| macro03-rd | 8 | 185.937 [180.459, 192.001] | 91.989 [91.353, 92.918] | 50.5% |
| large01 | 1 | 215.435 [212.947, 218.313] | 180.010 [178.097, 180.205] | 16.4% |
| large01 | 8 | 93.309 [91.529, 97.161] | 51.322 [50.725, 54.237] | 45.0% |
| large02 | 1 | 243.114 [241.486, 243.539] | 214.071 [210.978, 216.140] | 11.9% |
| large02 | 8 | 106.412 [100.685, 110.652] | 65.464 [64.084, 68.374] | 38.5% |
| large01-rd | 1 | 468.669 [468.176, 468.887] | 374.864 [373.615, 382.126] | 20.0% |
| large01-rd | 8 | 194.173 [177.380, 202.310] | 105.706 [104.200, 105.746] | 45.6% |

ibex_core includes detailed placement, legality, WNS/TNS, area and HPWL. Other configurations use LEF/DEF only; no timing or legalization claim.

Both jobs reached Completed, compact evidence was locally SHA-256 verified, exact aliases were deleted, and subsequent queries returned empty lists. Binary, input, source, archive and job receipts are indexed in cross-shard-summary.json.

Local evidence archives are stored outside the documentation directory:

- [flasheda-soa-large-remediation-20261106 evidence](/Users/shom/Downloads/openroad-gpl-soa/build/soa/remediation-large/evidence.tar.gz)
- [flasheda-soa-large-high-remediation-20261106 evidence](/Users/shom/Downloads/openroad-gpl-soa/build/soa/remediation-large-high/evidence.tar.gz)