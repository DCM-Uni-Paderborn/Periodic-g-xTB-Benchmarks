# Periodic g-xTB in CP2K: Part II validation data

This branch contains the reproducibility record for *Periodic g-xTB in CP2K.
II. Exact, Memory-Bounded, and Distributed Brillouin-Zone-Coupled Nonlocal
Exchange*.
Application-accuracy benchmarks and other Part-I material live on `main`; they
are intentionally not duplicated here.

Any benchmark retained on this branch measures a numerical acceleration
module---for example transform cost, bounded workspace, or MPI scaling---and
not the predictive accuracy of g-xTB.  DMC13 and LC10 are therefore outside
the Part-II data model.

Every optimized path is tested against the explicit expanded-full-mesh
implementation. That reference remains selectable as the numerical oracle.
An SCC endpoint alone is insufficient: energy, Fock response, overlap adjoint,
forces, and stress must agree before timing or memory results are accepted.

## Contents

- `campaigns/exchange-acceleration-20260716/`: machine-readable module and
  case matrix, raw outputs, deterministic comparisons, and exact source/build
  provenance for the reference and accelerated paths.
- `campaigns/cp2k-batched-exchange-runtime-20260716T160444Z/`: bounded
  image-batch correctness and memory/work trade-off tests.
- `campaigns/gxtb-exchange-distributed-separable-dft-20260716T170019Z/`:
  distributed image-range and separable direct-transform qualification.
- `validation/accelerated_exchange/`: immutable reference-equivalence,
  force/stress, 0D--3D PBC, K290/SPGLIB/time-reversal, MPI, cache, transform,
  restart, timing, and memory archives with SHA-256 manifests.
- `validation/accelerated_exchange/auto_policy_regression_20260722/`: final
  public-source/exact-binary regression for the default production policy;
  48 calculations and all 196 matcher results pass on the final CP2K and
  save_tblite revisions.
- `validation/accelerated_exchange/figures/`: publication-ready Part-II
  figures whose generating scripts and checksums are retained in this branch.
- `scripts/benchmark_execution.py`: fail-closed CPU reservation and affinity
  helper used for reproducible timing runs.
- `tests/`: repository-level archive, derivation, and CPU-affinity checks.

The retained acceleration components cover streamed symmetry-star
contractions, bounded provider image batches, bounded-forward/sparse-reverse
ACP mesh contractions, invariant phase/symmetry caches,
regular-mesh dense-oracle and separable/mixed-radix transforms, distributed
image kernels and MPI ownership, and metadata-validated cross-mesh restarts.
The automatic policy selects the exact qualified combination only after the
complete mesh and MPI layout are known; one-point ACP remains dense, while
periodic multi-point ACP uses bounded forward batches and sparse reverse.
Each component has an independent qualification boundary; a component-level
result is never promoted to an end-to-end scaling claim.

## Acceptance order

1. fixed-matrix expansion, transform, foldback, and adjoint identities;
2. identical-input energy, Fock, overlap-adjoint, force, and stress agreement;
3. converged-state identity and analytical derivatives versus finite
   differences;
4. shifted/unshifted meshes, RKS/UKS, 0D--3D PBC, time reversal, K290 and
   SPGLIB reduction, and multiple MPI layouts;
5. timing and peak/aggregate memory only after all applicable correctness
   gates pass.

The physical Brillouin-zone sum is unchanged. Streaming, caching,
transform factorization, and distribution alter storage, ordering, and
communication only; they do not replace transformed star members by scalar
weights or omit cross-k-point exchange terms.

## Repository checks

Run the self-contained standard-library checks with:

```bash
python3 -m unittest discover -s tests
```

Each immutable archive also documents its own checksum and regeneration
commands. Start with `validation/accelerated_exchange/README.md` and the three
campaign READMEs above. Historical source paths in manifests are provenance
strings only; all inputs used for the retained Part-II tests are archived
inside this branch.

## Scope boundary

This branch contains no DMC-ICE13, lattice-constant/cohesive-energy, or
molecular-crystal application benchmark. Those datasets remain on `main` and
in Git history. Complete periodic GFN1/GFN2 benchmark data remain canonical in
[`DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks`](https://github.com/DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks).
