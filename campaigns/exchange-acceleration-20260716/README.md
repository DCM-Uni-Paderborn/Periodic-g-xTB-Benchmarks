# Brillouin-zone-coupled exchange acceleration campaign

Status: active implementation and qualification campaign.

This campaign supplies the machine-readable evidence for
`main_jcp_II_v1.tex`, *Periodic g-xTB in CP2K. II. Accelerated and Scalable
Brillouin-Zone-Coupled Nonlocal Exchange*.  The explicit expanded-full-mesh
implementation described in Part I is the permanent numerical oracle.  An
optimized path is not accepted from a successful SCC endpoint alone: its
energy, Fock response, overlap adjoint, forces, and stress must all reproduce
the oracle before timing or memory data may enter the paper.

## Independently selectable modules

1. **Canonical regular-grid and Born--von Karman cache.**  Prove the complete
   shifted product mesh in linear time, cache exact order maps, phases and
   static exchange kernels, and invalidate on every geometry, lattice, mesh,
   basis, topology, or model change.
2. **Streamed symmetry-star contraction.**  Generate a full-mesh density or
   overlap block transiently from an irreducible representative, push it into
   the coupled exchange accumulator, and fold each returned response block
   immediately with the weighted real adjoint.  Time reversal must use the
   existing anti-linear CP2K transformation; scalar star weights are not a
   replacement for transformed matrices.
3. **Regular-grid FFT.**  Replace dense k-to-R and R-to-k transformations only
   after the grid planner has established a complete product mesh and its
   common twist.  Keep the dense transform selectable as the oracle and as the
   fallback for irregular inputs.
4. **MPI k-point groups.**  Distribute star, real-space, or AO-pair tiles with
   explicit ownership and reductions for energy, irreducible Fock blocks,
   overlap adjoints, forces, and stress.
5. **Validated cross-mesh restart.**  Transfer the density through a common
   real-space representation, restore symmetry, Hermiticity and electron
   number, reject incompatible metadata, and fall back to a cold start when
   the initial residual/energy gate fails.

The modules are tested separately and in all meaningful combinations.  The
full physical Brillouin-zone sum is never reduced by assumption; only its
storage, ordering, transformation and distribution are changed.

## Qualification order

1. fixed-matrix forward/inverse and weighted-adjoint identities;
2. single-step exchange energy, Fock and overlap response;
3. converged SCC endpoint and state identity;
4. analytical forces and stress against both the oracle and central finite
   differences;
5. 0D, 1D, 2D and 3D boundary conditions, shifted and unshifted meshes,
   time reversal, K290 and SPGLIB reduction, and a nonsymmorphic operation;
6. one versus several MPI k-point groups and repeated rank layouts;
7. cold versus transferred-density SCC runs;
8. wall time and peak/aggregate memory only after all applicable correctness
   gates have passed.

## Evidence layout

- `validation_matrix.json` is the authoritative module/case status index.
- `provenance/` contains immutable source, build, binary, input and launcher
  hashes.
- `raw/` contains unedited stdout, stderr, timing and maximum-RSS records.
- `derived/` contains deterministic comparisons generated from `raw/`.
- `paper/` contains only table/figure data whose source rows have passed the
  gates in `validation_matrix.json`.

No raw result is overwritten.  A rerun receives a new timestamped directory
and an independent checksum manifest.

## Qualified CP2K block helpers

The CP2K one-block expansion and weighted real-adjoint foldback helpers are now
qualified as a separately bounded component.  Their once-per-map runtime check
uses deterministic complex-Hermitian probes and enforces the following gates:

- physical overlap expansion versus the unchanged full-array path at
  `max(1e-6, 100*eps_geo)` relative residual;
- the weighted variational adjoint identity for every star member at `1e-10`,
  including antiunitary negative-`rotp` operations;
- the accumulated blockwise fold versus the unchanged full-array oracle at
  `1e-12` relative residual.

Five current/baseline run pairs cover K290 and SPGLIB `2x2x2` reductions and a
`3x1x1` time-reversal reduction.  All ten launchers returned zero and reached
`PROGRAM ENDED`.  The maximum current/baseline differences are
`7.105e-15 Ha` over all printed energy evaluations,
`1.024e-16 Ha/bohr` for analytical force components,
`2.633e-10 bar` for printed analytical stress, and `5.370e-10 Ha` for the
finite-difference numerical virial.  The largest current finite-difference
sum is `1.893e-9 Ha`, well inside the campaign derivative gate.

The deterministic derivation is
`scripts/compare_cp2k_block_helpers.py`; its JSON and compact text products are
under `derived/`.  Unedited outputs, input/output manifests, the exact source
patch, build logs, launchers, audit, and reference runtime harness are preserved
under `raw/cp2k_block_helpers/` and `provenance/cp2k_block_helpers/`.

This result does **not** qualify the production streamed exchange module.  It
establishes only the CP2K block operators that the streamed module will call;
the latter remains `implementation_in_progress` until its full provider and
CP2K forward/reverse paths pass the complete oracle matrix.

## Qualified save_tblite provider cache and matrix-lean forward stream

The provider cache/planner and matrix-lean forward transaction are now
qualified as two narrowly scoped components.  An exact-source GNU Fortran
Debug build (`-O0 -fcheck=all -fbacktrace`) passes the focused
`bvk_exchange_supercell` test with return code zero.  That single test is a
deliberately dense functional qualification rather than a smoke test: it
compares the BvK provider path with the explicit supercell and unchanged dense
full-mesh oracles, exercises cache hits and exact invalidation, rejects and
recovers from incomplete and duplicate grids, and checks a regular `9x9x1`
mesh beyond the small dense phase-orthogonality oracle.

For the reduced forward stream the same passing test verifies energy, shell
potential and every Fock block, arbitrary block arrival order, a common mesh
twist, physical k-point permutation, duplicate/missing-block recovery and
charge-dependent onsite-state invalidation/recovery.  The public storage query
is asserted six times in reduced mode (before and after application, across
ordered, permuted and `9x9x1` transactions) and confirms only that the stream
does not retain the complete k-space density and overlap input arrays.  A
separate assertion confirms that the deliberately retained oracle mode does
own those arrays.

This storage query is deliberately narrow: its implementation tests only
whether `stream%density` or `stream%overlap` is allocated.  It does not inspect
the remaining stream or cache allocations and therefore cannot establish
bounded memory.  The current forward path retains three complete complex
`nao x nao x Nk x nspin` BvK-image tensors (`amat_r`, `cmat_r`, and `vmat_r`),
while the provider cache retains two dense `Nk x Nk` phase tables.  Consequently
the present result is classified as **matrix-lean**, with no claim about total
peak-memory reduction or asymptotically bounded storage.

The deterministic derivation is
`scripts/summarize_save_tblite_provider_forward.py`; its JSON and compact text
products are under `derived/`.  Unedited stdout/stderr, return code and command
are under `raw/save_tblite_provider_forward/`; the cleanly applicable exact
patch, complete modified-source snapshot, CMake cache, build commands and
checksums are under `provenance/save_tblite_provider_forward/`.  An earlier
terok 30/30 record is preserved separately as historical raw data but is not a
qualification basis because its executable path was replaced during that
long-running test.

This passed status is intentionally limited to the provider cache/planner and
matrix-lean **forward** stream, specifically the absence of retained full
k-space density/overlap input arrays.  True bounded-memory R/image batching,
including batched or on-demand phase generation and measured peak-memory
scaling, remains `implementation_in_progress`.  The tested reduced stream also
rejects the reverse call by contract; reduced-memory overlap adjoints, forces
and stress therefore remain `implementation_in_progress`.  Likewise,
separately passing provider and CP2K block-helper tests do not yet qualify the
production CP2K stream consumer, which also remains
`implementation_in_progress`.
