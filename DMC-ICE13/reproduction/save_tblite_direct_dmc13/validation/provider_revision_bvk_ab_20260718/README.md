# Direct BvK provider-revision comparison

This gate compares two standalone `save_tblite` executables on identical
explicit Born--von-Karman supercells.  It therefore does not involve CP2K,
Bloch transforms, symmetry reduction, or a native-k normalization choice.

- `current_save_tblite_cli` contains the current provider results.  The same
  energies agree with the CP2K-native primitive-cell energies to about
  `1e-7` Hartree per primitive cell over the complete 13-phase `2^3` and
  `3^3` sets.
- `seidler_pbc_cli` contains independent recalculations with the unchanged
  Seidler `pbc` reference executable for ice Ih, VII, and XVII at `2^3` and
  `3^3`.

The Seidler executable changes the difficult ice-VII relative energy by only
`-0.09212` kJ mol^-1 per water at `2^3` and `-0.07717` kJ mol^-1 per water at
`3^3` relative to the current provider.  Ice XVII changes by at most
`0.00041` kJ mol^-1 per water.  These differences are far too small to explain
the earlier large DMC-ICE13 discrepancy or a substantially denser native-k
convergence tail.

`absolute_energy_comparison.csv` reports primitive-cell energies and
`relative_energy_comparison.csv` reports Ih-referenced values.  Positive
absolute-energy differences mean Seidler `pbc` minus current `save_tblite`.
The archived process outputs, JSON files, executable hashes, structure hashes,
and integrity manifest provide the complete provenance.

The Linux author-branch repetition in `seidler_pbc_cli_linux` extends both the
`2 x 2 x 2` and `3 x 3 x 3` comparisons to all 13 phases.  Its source tree is
byte-for-byte identical to revision `c932120`, and the three calculations
repeated on both macOS and Linux agree within `1.9e-12` Hartree for the complete
supercell.  At `2 x 2 x 2`, the current and author-branch MAEs are `88.68138`
and `88.70687` kJ mol^-1 per water, respectively, a change of only `+0.02549`
kJ mol^-1.  At `3 x 3 x 3`, they are `34.04849` and `34.07059` kJ mol^-1,
respectively, a change of `+0.02210` kJ mol^-1; the largest individual
author-minus-current relative-energy shift is `0.10048` kJ mol^-1 for ice VI.
Thus the final author branch is marginally worse on both complete meshes and
does not account for an accuracy or k-point-convergence improvement.
`full_k222_relative_comparison.csv`, `full_k333_relative_comparison.csv`, and
`full_k333_absolute_energy_comparison.csv` contain the phase-resolved evidence.

`compare_complete_mesh.py` applies the same absolute-energy, Ih-referenced
relative-energy, and MAE analysis to any complete 13-phase direct-BvK mesh.
It refuses missing, malformed, or non-finite JSON results and writes a
SHA-256 input manifest alongside the generated tables, so a denser author
comparison cannot silently mix executables or incomplete phases.
