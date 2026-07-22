# Accelerated periodic g-xTB exchange validation

This directory contains implementation-oracle evidence for Part II of the
periodic g-xTB work.  Application benchmarks are intentionally kept separate.
The explicit expanded-full-mesh implementation remains the numerical oracle.

| Archive | Validated scope | Performance claim |
|---|---|---|
| `acp_mesh_contraction_20260722/` | Bounded ACP Bloch batches and sparse projector-image reverse contraction; dense/streamed/qualify energy, force, stress, response-oracle, open-shell-neutrality, and fail-closed selector tests on macOS and Linux | Removes the complete ACP Bloch tensor and quadratic projector-difference set from production; scoped allocation claim only, with no total-RSS or speedup claim |
| `provider_streamed_reverse/` | Bounded save_tblite reverse transaction; overlap adjoint, force, stress, state-machine and exact workspace tests | Provider AO-matrix storage is bounded by the declared image-batch size; the unit suite is not a speedup benchmark |
| `cp2k_streamed_reverse_consumer/` | CP2K bounded derivative consumer; 21 Linux RKS/UKS, full/TR/K290/SPGLIB, shifted, 1D/2D/3D cases at MPI P=1/2/4 plus force/stress finite differences | Correctness and exact provider high-water are qualified; historical P=2/P=4 shared-mask timings are legacy/non-scaling and only numerical equivalence remains usable |
| `cp2k_streamed_star_memory/` | Native-mixer symmetry-star gate; 48 Linux DENSE/STREAMED/QUALIFY runs spanning RKS/UKS, K290/SPGLIB/TR, 1D/2D/3D, MPI P=1/2/4, force, and stress | Removes the remaining `N^2 S K` post-mixer full-mesh temporary in favor of exactly `3 N^2` work elements; this is not a whole-process RSS or speedup claim |
| `cp2k_kgroup_owner/` | CP2K owner/communicator precursor; full/TR/K290/SPGLIB, shifted, 1D/2D/3D, RKS/UKS, and Linux MPI P=1/2/4 oracle comparisons | No coupled-kernel speedup: save_tblite still has one global stream state |
| `mixed_radix_fft_20260717/` | Dense-oracle equivalence of separable and mixed-radix FFT regular-mesh exchange transforms for RKS/UKS, shifted/unshifted K290/TR/SPGLIB, 1D/2D/3D, energy, force and stress | Correctness qualification only; the short serial matrix is not a speedup benchmark |
| `symmetry_phase_cache_20260717/` | Cached-versus-uncached symmetry-star atom phases across RKS/UKS, K290/TR/SPGLIB, shifted meshes and derivative diagnostics | Correctness qualification only; the short serial matrix is not a speedup benchmark |
| `cp2k_distributed_images_20260717/` | Replicated reduced source plus disjoint nonlinear BvK-image kernels and all-k additive response folding; 30 dense/distributed pairs and six fail-closed cases spanning RKS/UKS, Gamma/shifted, K290/symmetry, 1D/2D/3D, MPI P=1/2/4, forces, and stress | Numerically equivalent to the full-mesh oracle and capable of coupled-kernel parallelism; the archived single-shot timings are not yet a repeated scaling benchmark |
| `crossmesh_restart/` | Opt-in, metadata-validated regular-mesh density restart with electron-conserving overlap-metric spectral projection; eight accepted and six fail-closed RKS/UKS, K290/SPGLIB, shifted, 1D/2D/3D, force/stress, density, and MPI P=1/2/4 controls | Saves one to five SCC iterations in the qualified cases; this is not a repeated wall-time or process-memory claim |

The distributed-image implementation preserves the required global coupling:
it first merges additive partial k-to-R accumulators, then evaluates disjoint
subsets of the nonlinear image kernel, and finally sums every image-local
contribution to the full Bloch response.  Independently finalized k-group
calculations must still never be summed because that would omit cross-group
terms in Brillouin-zone-coupled nonlocal exchange.  The provider boundary and
its validation rules are recorded in
`cp2k_kgroup_owner/provider_partial_accumulator_abi.md`.

The frozen integrated P=2/P=4 timing matrix and other multi-rank runs launched
with taskset plus `--bind-to none` are preserved byte-for-byte, but are
classified `legacy_timing_non_scaling`. Their energy, Fock, force, and stress
comparisons remain valid; no speedup or scalability claim may use their times.

Each subdirectory carries raw output, verification scripts, provenance, and
SHA-256 manifests.  Run those scripts from within the corresponding archive
after copying or publishing the data.

`fft_phase_cache_crossmesh_20260717/README.md` is the compact index for the
three corresponding acceleration archives and their deliberately separate
qualification boundaries.
