# K290/SPGLIB reduction audit

This package tests what the current qualified periodic g-XTB implementation
actually reduces, rather than inferring memory or timing from the number of
irreducible k points.  A periodic CH4 cell on a regular `6 x 6 x 6` mesh was
run three times through five routes with one exactly pinned CPU and one thread
for OpenMP and every recorded BLAS runtime.

All 15 calculations use CP2K executable SHA-256
`b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f`,
converge in 11 SCC steps, terminate normally, and agree with the explicit
216-point result within `7e-15` hartree.  The verifier independently checks
the executable and input hashes, affinity proof, thread environment, selected
symmetry backend, k-point count, normal termination, and fused-storage marker.

Median results over the three interleaved rounds are:

| Route | SCF k points | Elapsed (s) | Change vs full | Peak RSS (KiB) | Change vs full |
|---|---:|---:|---:|---:|---:|
| Explicit full / complete array | 216 | 64.0341 | 0.00% | 1,062,768 | 0.00% |
| K290 / complete array | 10 | 57.5782 | -10.08% | 1,054,904 | -0.74% |
| SPGLIB / complete array | 10 | 57.1455 | -10.76% | 1,055,780 | -0.66% |
| K290 / symmetry fused | 10 | 58.1330 | -9.22% | 1,053,336 | -0.89% |
| SPGLIB / symmetry fused | 10 | 58.1174 | -9.24% | 1,051,892 | -1.02% |

K290 and SPGLIB therefore reduce the SCF/diagonalization set from 216 to 10
points (95.37%) and save about 10% elapsed time in this small serial case.
They do **not** reduce whole-process memory by a comparable factor with the
default `COMPLETE_ARRAY` exchange backend, because the coupled 216-point
exchange mesh is still reconstructed and materialized.

The optional `SYMMETRY_FUSED` route prints and passes
`persistentFullExchangeStarMatrices=0`, with an eight-member batch and only
ten cached representative overlaps.  This is structural evidence that the
persistent full exchange-star arrays have been removed.  The total-process
RSS falls by only about 1% here because other CP2K and model-library storage
dominates, and the fused path is slightly slower than the default reduced
route.  No general memory or speedup claim should be extrapolated from this
single small system.

`online_head_delta.json` records the online-source reconciliation performed
before the audit.  The two later CP2K commits affect restart alias projection
and regular-mesh inference for explicit `GENERAL` lists; these cold-start
`MONKHORST-PACK` inputs enter neither path.  The online save_tblite production
source tree is identical to the qualified provider tree.

Reproduce the machine-readable decision with:

```bash
python3 verify_kpoint_reduction.py --output verification.reproduced.json
```

The archived first launcher attempt is retained remotely as a setup failure:
`/usr/bin/time` was unavailable and no CP2K calculation was started.  It is
not included in the 15 scientific runs or the reported medians.
