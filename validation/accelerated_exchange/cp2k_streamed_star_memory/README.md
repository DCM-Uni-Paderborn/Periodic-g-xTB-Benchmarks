# CP2K streamed symmetry-star memory qualification

Date: 2026-07-17

This archive qualifies an isolated CP2K memory optimization for the native
g-xTB Fock mixer.  It contains implementation-oracle evidence for Part II and
no LC10, DMC13, or other application benchmark.

Only the release qualification labeled V4 is authoritative.  V1 lacked the
MPI launcher in its noninteractive path, V2 was rejected after PRRTE overrode
the requested affinity, and V3 passed numerically but exercised the superseded
four-work-block implementation.  Those attempts are excluded from every table
and archive here.  `selector_edges_v2` and `debug_subset_v4` are component-run
labels for the final three-block source and are part of the authoritative
evidence bundle.

## Optimization and selectors

After native Fock mixing on a symmetry-reduced mesh, the reference path
materializes `mixed_fock_full(N,N,S,K)` only to validate symmetry covariance.
The new streamed path expands one member of one symmetry star, validates every
operation-coset image, and immediately applies the weighted real adjoint.  It
therefore preserves the calibrated CP2K `U_g`, `w_full/w_irred` fold factor,
and antiunitary time reversal without storing the complete mesh.

`CP2K_GXTB_SYMMETRY_STAR_CONTRACTION` is independent of the existing forward
exchange and reverse-derivative selectors:

- `DENSE` is the default and permanent complete-mesh oracle;
- `STREAMED` uses three reusable `N x N` complex work blocks;
- `QUALIFY` runs both paths in the same mixer call and aborts above the
  existing covariance or weighted-round-trip gates.

Unknown and overlong values abort.  With `N=nao`, `S=nspin`, and `K=nfull`, the
targeted dense temporary contains exactly `N^2 S K` complex elements, whereas
the streamed temporary contains exactly `3 N^2`.  The reduction factor is
therefore `S K / 3`.  This is an exact allocation claim for this post-mixer
gate, not a whole-process RSS or speedup claim.

The allocation/lifetime audit in `source_snapshot/STREAMED_STAR_MEMORY_AUDIT.md`
also records that the existing forward exchange, Fock foldback, and bounded
reverse exchange already stream full-star members; they were not changed by
this patch.

## Release qualification

The authoritative Linux matrix contains 48 clean calculations: 16 case/rank
configurations, each run with `DENSE`, `STREAMED`, and `QUALIFY`.  Coverage is:

- RKS and UKS;
- K290 and SPGLIB symmetry reduction;
- unitary and antiunitary time-reversal stars;
- true 1D, 2D, and 3D periodicity;
- MPI `P=1,2,4`, including `P > nred`;
- energy, analytical force, and analytical stress;
- the existing dense forward/reverse oracles through
  `KGROUP_PARTIAL_ROOT`, derivative `QUALIFY`, and the iteration-one
  full-mesh oracle.

All 16 triples have zero difference at printed precision for energy, every
force component, and every stress component.  The largest internally printed
covariance or weighted expand/fold residual is `6.661338e-16`, versus the
`1.0e-10` gate.

Representative exact memory counters are:

| Case | `S` | `K` | Dense complex elements | Streamed complex elements | Factor |
|---|---:|---:|---:|---:|---:|
| RKS K290/SPGLIB | 1 | 8 | 512 | 192 | 2.667 |
| RKS time reversal | 1 | 3 | 12 | 12 | 1.000 |
| UKS time reversal | 2 | 3 | 384 | 192 | 2.000 |
| RKS 2D time reversal | 1 | 9 | 11664 | 3888 | 3.000 |

Every run used `mpiexec --bind-to none` inside one of eight disjoint four-CPU
sets (`192-195` through `220-223`).  The runner sampled every live CP2K rank
from `/proc`, required its exact `Cpus_allowed_list`, and required its current
processor (`/proc/PID/stat` field 39) to be inside that set.  Missing or escaped
rank affinity is a hard failure.  The complete per-rank proof is stored in each
`run.json` inside the raw archive.

The same single-shot small-system matrix records wall time but was not designed
as a performance benchmark.  `STREAMED/DENSE` ranges from `0.944` to `1.080`,
with median `0.971`; this spread is launch and SCF noise, so no speedup is
claimed.  The exact allocation counters, rather than process RSS or these wall
times, support the memory statement.  The individual values are frozen in
`timing_summary.tsv`.

## Debug and fail-closed gates

Four Debug/check calculations passed: K290 RKS at `P=2`, UKS time reversal at
`P=4`, true 1D time reversal at `P=4`, and true 2D time reversal at `P=4`.
They retained bounds, FPE, undefined-behavior, and Fortran runtime checks.
Leak reporting was disabled only for known Open MPI/PMIx process-finalization
allocations; the failed preliminary run retaining those third-party reports is
not part of the authoritative archive.

The default-without-environment-variable run is exactly identical to explicit
`DENSE` for energy, force, and stress.  `BOGUS` and a 33-character value both
terminate nonzero with their specific fail-closed diagnostic.  These runs also
carry live rank-affinity proofs.

## Frozen provenance

| Artifact | SHA-256 |
|---|---|
| Tested `src/tblite_interface.F` | `996bb0cd04b4f8e8d560886a2963d90efd075af39b2f4bb9bd9e0ad930c2703a` |
| Release `libcp2k.so.2026.2` | `f22762b82b1c6be6adf4f4fde8aa476a1d350d1283112ad9507e360d0f8ff861` |
| Debug `libcp2k.so.2026.2` | `2d3acc1d47427a0df7ef4b689b885baa427b044424a14538101cbf2884c6fb10` |
| Release `cp2k.psmp` launcher | `16c4b8905943cb30b11425b56b10b98af62e538a93776e817d6ae230eeceeef2` |
| Debug `cp2k.pdbg` launcher | `9d983212849417543180ca6b560607c1abe9f2b8d8fc021ac01160208a31d92b` |
| Frozen `save_tblite` `libtblite.a` | `20c74bf3272a229e125893956880757bca365c0dc1b54fa8e892bf99a67e7760` |
| SPGLIB shared library | `8fdf95db83704cfc8da870062a1292380038d9d864a90b4a8b495d5b239572ce` |
| Complete raw evidence archive | `5f8aaf36a25f96048850868d97ca3719cfa5a55005390ed9f251c3f4acba7745` |

The source starts from CP2K g-xTB-pbc base
`68f677114f5829a32292171251150dd8e00ce458`; repository commit identifiers are
kept here for reproducibility and are not intended for the manuscript or SI.

Run `scripts/verify_test_matrix.py` from this archive directory to recheck
hashes, termination, live affinity, observables, internal residuals, and exact
allocation formulas.  The verifier rejects archive links and paths escaping its
temporary extraction root, reads the frozen release qualification directly
from `raw_archive/`, and verifies the recorded CP2K executable and shared-library
hashes against the curated provenance manifests; no original Terok path is
required.  `SHA256SUMS` covers every curated artifact in this directory except
itself.
