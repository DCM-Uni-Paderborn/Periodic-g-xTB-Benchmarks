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

