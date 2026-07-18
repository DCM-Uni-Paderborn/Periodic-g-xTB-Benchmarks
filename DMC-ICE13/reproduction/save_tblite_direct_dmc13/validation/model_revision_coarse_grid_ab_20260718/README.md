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
- `authors_exchange`: Leopold Seidler's `mstore-inorganic` history after
  the periodic exchange-cutoff and electrostatics changes;
- `gxtb_v201`: the separate post-March molecular g-xTB model history;
- `dcm_main`: the then-current DCM `main` merge on the newer tblite base.

The raw result is unambiguous: the coarse-grid DMC13 MAE depends strongly on
which source history is called “the g-xTB CLI.”  The current integration and
DCM `main` are close, whereas the two historical development lines differ by
tens to hundreds of kJ mol^-1 per water molecule at these deliberately
unconverged grids.  Therefore a smaller (1^3) or (2^3) MAE from one of
those lines is not evidence for a CP2K-native k-point error.

These values are diagnostic only.  They must not replace the converged
native-k-point benchmark.  The decisive follow-up is the full `3 x 3 x 3` and
denser comparison using the exact authors' executable, followed by a
component-level port only if the difference survives convergence and passes
energy, force, and stress validation.

The selected `2 x 2 x 2` `authors_exchange` results were repeated on Linux.  Ice
Ih, VII, and XVII agree with the macOS build to within
9.1e-13 hartree in the supercell energy, excluding a platform
or BLAS explanation.

Run `python3 verify_comparison.py` to recompute every reported MAE from the
raw JSON files and to verify the cross-platform sentinel.
