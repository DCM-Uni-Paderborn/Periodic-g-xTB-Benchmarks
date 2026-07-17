# Accelerated periodic g-xTB exchange validation

This directory contains implementation-oracle evidence for Part II of the
periodic g-xTB work.  Application benchmarks are intentionally kept separate.
The explicit expanded-full-mesh implementation remains the numerical oracle.

| Archive | Validated scope | Performance claim |
|---|---|---|
| `provider_streamed_reverse/` | Bounded save_tblite reverse transaction; overlap adjoint, force, stress, state-machine and exact workspace tests | Provider AO-matrix storage is bounded by the declared image-batch size; the unit suite is not a speedup benchmark |
| `cp2k_streamed_reverse_consumer/` | CP2K bounded derivative consumer; 21 Linux RKS/UKS, full/TR/K290/SPGLIB, shifted, 1D/2D/3D cases at MPI P=1/2/4 plus force/stress finite differences | Correctness and exact provider high-water are qualified; the small-system process-RSS comparison does not establish a material whole-process reduction or speedup |
| `cp2k_kgroup_owner/` | CP2K owner/communicator precursor; full/TR/K290/SPGLIB, shifted, 1D/2D/3D, RKS/UKS, and Linux MPI P=1/2/4 oracle comparisons | No coupled-kernel speedup: save_tblite still has one global stream state |
| `cp2k_distributed_images_20260717/` | Replicated reduced source plus disjoint nonlinear BvK-image kernels and all-k additive response folding; 30 dense/distributed pairs and six fail-closed cases spanning RKS/UKS, Gamma/shifted, K290/symmetry, 1D/2D/3D, MPI P=1/2/4, forces, and stress | Numerically equivalent to the full-mesh oracle and capable of coupled-kernel parallelism; the archived single-shot timings are not yet a repeated scaling benchmark |
| `crossmesh_restart/` | Opt-in, metadata-validated regular-mesh density restart; official regression, malformed-file fallback, RKS/UKS, force/stress, density and Linux MPI-2 tests | Reduces SCF iterations in the qualified cases; the interpolated initial guess is not guaranteed N-representable, so default activation awaits a spectral projection or accept/fallback gate |

The distributed-image implementation preserves the required global coupling:
it first merges additive partial k-to-R accumulators, then evaluates disjoint
subsets of the nonlinear image kernel, and finally sums every image-local
contribution to the full Bloch response.  Independently finalized k-group
calculations must still never be summed because that would omit cross-group
terms in Brillouin-zone-coupled nonlocal exchange.  The provider boundary and
its validation rules are recorded in
`cp2k_kgroup_owner/provider_partial_accumulator_abi.md`.

Each subdirectory carries raw output, verification scripts, provenance, and
SHA-256 manifests.  Run those scripts from within the corresponding archive
after copying or publishing the data.
