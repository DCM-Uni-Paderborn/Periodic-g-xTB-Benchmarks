# Bounded reverse exchange qualification

This directory preserves the provider-level qualification of the bounded
reverse/derivative transaction used by periodic g-xTB.  The dense complete-mesh
implementation remains the numerical oracle.

The frozen source was tested independently on macOS and on `terok`
(Debian GNU/Linux, x86-64, GNU Fortran 15.2.0).  On terok the release exchange
suite passes 31/31 tests and the complete g-xTB suite reports 40 ordinary
passes plus four expected-failure diagnostics.

## Dense-oracle residuals on terok

| Case | Max overlap-adjoint error | Max force-response error | Max stress-response error | Peak complex scalars |
|---|---:|---:|---:|---:|
| RKS 3x1x1, batch 2 | 1.3878e-17 | 0 | 0 | 113 |
| UKS shifted/permuted 9x9x1, batch 7 | 5.5511e-17 | 0 | 6.7763e-21 | 585 |
| RKS shifted/permuted 2x2x2, batch 3 | 2.7756e-17 | 0 | 0 | 144 |

The transaction-owned complex high-water count is

`6*B*nAO^2*nspin + 5*nAO*nspin + 2*Nk + 3*max(nmesh) + 10*nAO^2`.

Thus the AO-matrix part is bounded by the image-batch size rather than by the
full k-point count.  Shared immutable real model data and compiler-created
temporaries are excluded from this exact provider counter; process RSS is
recorded separately in the JSON files.

## Controlled-thread Linux runs

| OMP threads | Exchange elapsed / s | Max RSS / KiB | Result |
|---:|---:|---:|---|
| 1 | 1.069 | 12,528 | pass |
| 2 | 1.495 | 12,852 | pass |
| 4 | 2.815 | 56,756 | pass |
| 8 | 12.814 | 188,596 | pass |

All four logs are byte-identical and therefore also establish deterministic
thread-count-independent numerical output for this test.  The small unit
suite is not a scaling benchmark: OpenMP overhead and per-thread runtime
storage dominate above one thread, so these timings must not be presented as
a speedup result.  Production MPI/CP2K timing is qualified separately.

Raw transcripts and `resource.getrusage` records are under
`terok_20260716/evidence-linux/`.  `streamed_reverse_qualification_20260716.md`
describes the API, state-machine negatives, and memory accounting.

