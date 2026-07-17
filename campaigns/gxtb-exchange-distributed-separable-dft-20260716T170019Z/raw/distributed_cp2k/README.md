# CP2K distributed g-xTB exchange qualification

Date: 2026-07-16

This directory records the first end-to-end qualification of
`CP2K_GXTB_EXCHANGE_STREAM_MODE=distributed` against the pre-existing
`legacy` path and the explicit full-mesh exchange oracle.  All production
comparisons used the debug CP2K executable, one OpenMP thread, and either one
or two MPI ranks.  The distributed runs use a contiguous, exactly-once
partition of the full Born--von Karman image set, followed by global sums of
the partial exchange energy, shell potential, and folded Fock response.

## Covered gates

- K=1 with P=2, including the inactive-rank (`P>K`) provider-cache path.
- K=8 full mesh and reduced SPGLIB/K290 paths in three dimensions.
- Uneven K=3/P=2 image ownership, time-reversal reduction, and UKS.
- True one-dimensional (`PERIODIC X`) and two-dimensional (`PERIODIC XZ`)
  cells, with printed force and analytical-stress components.
- A shifted 2x2x2 mesh.
- Componentwise force and stress comparisons and full DEBUG finite-
  difference force/stress checks for the 3D full, SPGLIB, and K290 cases.
- Per-iteration oracle checks of exchange energy, shell potential, folded
  Fock response, Hermiticity, symmetry covariance, and adjoint duality.

All completed production controls ended normally.  Final energies agree
exactly in the K=1, uneven K=3/P=2, one-dimensional, two-dimensional,
shifted-mesh, and UKS pairs.  Across all 23 displaced geometries of each 3D
DEBUG pair, the largest legacy/distributed energy difference is
7.11e-15 hartree.  The largest printed componentwise force difference is
3.33e-16 hartree/bohr, and the largest analytical-stress difference is
1.40e-6 bar.  The full-mesh oracle bounds the exchange-only differences by
2.66e-15 hartree for energy and 3.05e-16 hartree for the folded Fock response.

Symmetry-covariance residuals are normally at roundoff.  During finite-
difference cell/coordinate displacements in the reduced SPGLIB and K290
DEBUG runs they reach 8.85e-12, while energy, Fock, duality, force, and stress
comparisons remain within the bounds above.  This is retained explicitly in
the tabular record rather than hidden by a looser aggregate tolerance.

The full lower-dimensional `RUN_TYPE DEBUG` numerical-stress driver traps
after its initial analytical evaluation in both the pre-existing legacy and
new distributed modes.  The corresponding `ENERGY_FORCE` calculations,
which exercise and print analytical forces and stress, complete normally and
match.  Therefore the trap is recorded as a pre-existing DEBUG deformation
diagnostic limitation, not as a distributed-exchange regression.

`qualification_summary.tsv` is the machine-readable comparison table.
`SHA256SUMS` covers all preserved inputs, outputs, and summaries except the
manifest itself.

