# CP2K consumer qualification: separable direct DFT backend

This archive qualifies the selectable `SEPARABLE_DFT` backend of the periodic
g-xTB Brillouin-zone-coupled nonlocal exchange path against the permanent
`DENSE` transform oracle in CP2K/Quickstep.

The accelerated backend is a separable **direct DFT**, not an FFT.  On a
machine-regular `n1 x n2 x n3` mesh it applies one-dimensional character
transforms along the active directions and costs
`O(N_AO^2 N_k (n1+n2+n3))`.  The dense path remains independently selectable.

## Scope

- full and symmetry-reduced meshes (SPGLIB and K290);
- unshifted and shifted regular meshes;
- 1D, 2D, and 3D periodic boundary conditions;
- RKS and UKS/time-reversal reduction;
- total energy, folded Fock response, analytical forces, and analytical stress;
- a 23-energy finite-difference force/stress DEBUG trajectory;
- post-install MPI-distributed smoke tests with two ranks, including the
  combined distributed-forward/separable-derivative selection.

All paired production runs contain `PROGRAM ENDED`.  The maximum paired
deviations are tabulated in `qualification_summary.tsv`.  The complete-mesh
DEBUG trajectory differs by at most `7.105427357601002e-15` Ha over all 23
energies.  Its force difference is `7.58372569e-17` Ha/bohr and its analytical
stress difference is `3.9999940781854093e-07` bar.  The finite-difference
stress residual remains `8.22e-10` a.u. or smaller for the separable backend.

The completed symmetry-reduced dense-oracle diagnostics report at most
`5.551115e-17` Ha in both the full and folded Fock responses.  The provider's
independent transform-only timing and exact storage results are archived in
the companion repository campaign.  A separate release CP2K `6x6x6` CH4
probe with 11 KS builds is stored under `release_bench/`.  Across three
sequential runs the median exclusive `build_tblite_ks_matrix` timer decreases
from 0.185 to 0.154 s (1.20x), while its inclusive timer changes from 1.131 to
1.103 s (1.03x) and total wall time is indistinguishable at about 7.18 s.
Thus this small-system probe demonstrates reduced transform work but does not
support an end-to-end SCF speedup claim.

The MPI-distributed forward contraction remains a separate algorithm from the
transform backend.  Its post-install two-rank smoke test is included here, but
the analytical derivative contraction still follows the complete-mesh path;
no MPI derivative scaling claim is made.

## Reproduction environment

- CP2K binary: `/private/tmp/cp2k_gxtb_batch_build/bin/cp2k.pdbg`
- CP2K data: `/private/tmp/cp2k_gxtb_symmetry_v2/data`
- provider installation: `/private/tmp/save_tblite_batched_install`
- `OMP_NUM_THREADS=1`
- primary mode: `CP2K_GXTB_EXCHANGE_STREAM_MODE=legacy`
- oracle: `CP2K_GXTB_EXCHANGE_TRANSFORM_MODE=dense`
- candidate: `CP2K_GXTB_EXCHANGE_TRANSFORM_MODE=separable_dft`

Inputs are under `inputs/`; raw transcripts are under `runs/`.
Release timing transcripts and their machine-readable summary are under
`release_bench/`.
