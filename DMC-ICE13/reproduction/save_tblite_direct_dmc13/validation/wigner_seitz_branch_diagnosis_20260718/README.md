# Wigner--Seitz branch diagnosis

This diagnostic resolves why a direct CLI calculation initially attributed to
the save_tblite authors produced a much smaller coarse-grid DMC-ICE13 error than
the CP2K-integrated provider.

## Branch identity

The compared author histories are not successive names for the same branch.
They diverge after a common ancestor:

- `mstore-inorganic` ends at a 24 June 2026 revision and retains compact
  Wigner--Seitz image positions after filtering the zero translation.
- `pbc` ends at a 16 July 2026 revision and contains Leopold Seidler's 10 July
  fix that maps each compact distance entry back to its original lattice
  translation.

`source_identity.txt` records the exact revisions. The final author reference
for periodic calculations is therefore `pbc`, not `mstore-inorganic`.

## Controlled source A/B test

Only two lines in the current integration source were changed for the decisive
diagnostic: `orig(pos)` was replaced by `pos` in the Wigner--Seitz image list.
This reproduces the old compact-position semantics while leaving the model,
input structures, compiler, SCF settings, and all other source code unchanged.
A separate threshold-only build changed the image-distance tolerance without
changing the index mapping.

For the explicit Gamma-centred `2 x 2 x 2` supercells, the compact-index build
moves ice VII from about `-300` to about `-146` kJ mol-1 per water relative to
ice Ih, close to the `mstore-inorganic` value of about `-152`. The threshold-only
build remains at the current value. With all Coulomb terms removed, the
compact-index build and `mstore-inorganic` agree to 0.08 kJ mol-1 for ice VII
and below 0.001 kJ mol-1 for ice XVII. This isolates the large discrepancy to
the legacy self-image mapping; the remaining full-model difference is the
separate hard-versus-smooth Wigner--Seitz electrostatics update.

The legacy mapping is not a valid accuracy improvement. Running the current
Wigner--Seitz unit tests against it fails three invariants: a cubic self-image
points to the origin, and orthorhombic and skew cells omit a nearest image. The
fixed mapping passes the same suite. Raw test logs and controlled CLI outputs
are archived under `raw/`.

## Consequence for DMC-ICE13

The existing `provider_revision_bvk_ab_20260718` gate compares the actual final
`pbc` executable with the current provider. Across all 13 phases at `2 x 2 x 2`,
their MAEs are 88.70687 and 88.68138 kJ mol-1 per water, respectively. The
author branch changes the MAE by only +0.02549 kJ mol-1 and is therefore
slightly worse at this deliberately unconverged mesh. Selected `3 x 3 x 3`
relative energies differ by at most 0.0772 kJ mol-1 for ice VII and 0.00003
kJ mol-1 for ice XVII.

Thus the previously reported lower coarse-grid error belongs to the obsolete
`mstore-inorganic` self-image bug, not to the final author `pbc` implementation
and not to a CP2K-native k-point normalization error. Production calculations
must retain the corrected original-translation mapping.

`relative_energy_comparison.csv`, `no_coulomb_relative_energy_comparison.csv`,
and `wsc_self_image_fix.patch` provide the numerical and source-level evidence.
