# Accelerated periodic g-xTB exchange validation

This directory contains implementation-oracle evidence for Part II of the
periodic g-xTB work.  Application benchmarks are intentionally kept separate.
The explicit expanded-full-mesh implementation remains the numerical oracle.

| Archive | Validated scope | Performance claim |
|---|---|---|
| `provider_streamed_reverse/` | Bounded save_tblite reverse transaction; overlap adjoint, force, stress, state-machine and exact workspace tests | Provider AO-matrix storage is bounded by the declared image-batch size; the unit suite is not a speedup benchmark |
| `cp2k_streamed_reverse_consumer/` | CP2K bounded derivative consumer; 21 Linux RKS/UKS, full/TR/K290/SPGLIB, shifted, 1D/2D/3D cases at MPI P=1/2/4 plus force/stress finite differences | Correctness and exact provider high-water are qualified; the small-system process-RSS comparison does not establish a material whole-process reduction or speedup |
| `cp2k_streamed_star_memory/` | Native-mixer symmetry-star gate; 48 Linux DENSE/STREAMED/QUALIFY runs spanning RKS/UKS, K290/SPGLIB/TR, 1D/2D/3D, MPI P=1/2/4, force, and stress | Removes the remaining `N^2 S K` post-mixer full-mesh temporary in favor of exactly `3 N^2` work elements; this is not a whole-process RSS or speedup claim |
| `cp2k_kgroup_owner/` | CP2K owner/communicator precursor; full/TR/K290/SPGLIB, shifted, 1D/2D/3D, RKS/UKS, and Linux MPI P=1/2/4 oracle comparisons | No coupled-kernel speedup: save_tblite still has one global stream state |
| `crossmesh_restart/` | Opt-in, metadata-validated regular-mesh density restart; official regression, malformed-file fallback, RKS/UKS, force/stress, density and Linux MPI-2 tests | Reduces SCF iterations in the qualified cases; the interpolated initial guess is not guaranteed N-representable, so default activation awaits a spectral projection or accept/fallback gate |

Genuine k-group scaling requires additive partial k-to-R accumulators to be
merged before the nonlinear provider kernel is applied.  The exact proposed
boundary and its validation rules are frozen in
`cp2k_kgroup_owner/provider_partial_accumulator_abi.md`.  Independently
finalized group calculations must never be summed because that would omit
cross-group terms in Brillouin-zone-coupled nonlocal exchange.

Each subdirectory carries raw output, verification scripts, provenance, and
SHA-256 manifests.  Run those scripts from within the corresponding archive
after copying or publishing the data.
