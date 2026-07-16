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

