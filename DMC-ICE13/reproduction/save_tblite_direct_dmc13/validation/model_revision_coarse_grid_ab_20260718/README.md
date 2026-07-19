# Coarse-grid g-xTB model/revision A/B diagnostic

This diagnostic separates CP2K integration effects from differences between
the concurrently developed g-xTB and periodic `save_tblite` source histories.
All calculations use the same absolute Cartesian DMC-ICE13 structures and the
same direct-CLI protocol:

```text
tblite run --method gxtb --acc 0.01 --iterations 300 --no-restart
```

The `1 x 1 x 1` and `2 x 2 x 2` inputs are explicit Gamma-centred Born--von Karman
supercells.  Energies are divided by the replication factor before relative
energies per water molecule are formed with respect to ice Ih.

The compared executables are:

- `current`: the CP2K-integration provider used by the native-k-point work;
- `legacy_mstore_inorganic`: Leopold Seidler's obsolete `mstore-inorganic`
  history after the periodic exchange-cutoff and electrostatics changes;
- `gxtb_v201`: the separate post-March molecular g-xTB model history;
- `dcm_main`: the then-current DCM `main` merge on the newer tblite base.

The raw result is unambiguous: the coarse-grid DMC13 MAE depends strongly on
which source history is called “the g-xTB CLI.”  The current integration and
DCM `main` are close, whereas the two historical development lines differ by
tens to hundreds of kJ mol^-1 per water molecule at these deliberately
unconverged grids.  Therefore a smaller (1^3) or (2^3) MAE from one of
those lines is not evidence for a CP2K-native k-point error.

The complete Linux `3 x 3 x 3` repetition adds the final author `pbc` branch
to this controlled comparison.  The current integration and final `pbc` MAEs
are `34.04849` and `34.07059` kJ mol^-1 per water, respectively.  The obsolete
`mstore-inorganic` history instead gives `17.83062` kJ mol^-1 per water.  Its
apparently better value is therefore source-history dependent; the largest
legacy/final relative-energy shift is `76.52863` kJ mol^-1 for ice VII.

These values are diagnostic only.  They must not replace the converged
native-k-point benchmark.  The subsequent
`../wigner_seitz_branch_diagnosis_20260718` gate shows that `authors_exchange`
is not the final periodic author reference: it retains a Wigner--Seitz
self-image indexing artifact fixed by Leopold Seidler in the newer `pbc`
branch.  The final `pbc` executable is compared directly in
`../provider_revision_bvk_ab_20260718` and agrees closely with the current
provider.

The selected `2 x 2 x 2` legacy results were repeated on Linux.  Ice Ih, VII,
and XVII agree with the macOS build to within
9.1e-13 hartree in the supercell energy, excluding a platform
or BLAS explanation.

Run `python3 generate_comparison_tables.py` to regenerate the tables and
`python3 verify_comparison.py` to recompute every reported MAE independently,
verify complete `3 x 3 x 3` coverage, normal termination, singleton-CPU
affinity, source and executable identities, all input hashes, and the
cross-platform sentinel.
