# Direct save_tblite reproduction of DMC-ICE13

This package permits an implementation-independent comparison between native
CP2K k-point calculations and direct `save_tblite` calculations on the exactly
commensurate Born--von-Karman supercells.

## Structures

`structures/kNNN/<phase>/POSCAR` contains Cartesian coordinates in Angstrom.
For mesh size `N`, the lattice and atom list are explicit `N x N x N`
replications of the primitive DMC-ICE13 cell.  Hence a Gamma-only direct CLI
calculation on that POSCAR is the real-space counterpart of the native,
Gamma-centred `N x N x N` Brillouin-zone mesh.

The 13 phases are `Ih`, `II`, `III`, `IV`, `VI`, `VII`, `VIII`, `IX`, `XI`,
`XIII`, `XIV`, `XV`, and `XVII`.  The structures were generated from the exact
inputs used for the current CP2K-native comparison.  An independent geometric
audit found them equivalent to the supplied `gamma_only_dmc_ice13` archive to
better than `1.2e-9` Angstrom in the periodic pair-distance spectra.

`provenance/structure_hashes.csv` records the SHA-256 digest of every POSCAR
and its source CP2K input.

The package is intentionally self-contained: the structures contain absolute
Cartesian positions rather than scaled coordinates, and no CP2K parser is
needed for a direct recalculation.

## Recalculation

Set the direct CLI executable and run, for example,

```bash
export TBLITE_EXE=/path/to/save_tblite/build/app/tblite
MESHES="1 2 3" ACCURACY=0.1 ./run_save_tblite.sh
```

`PHASES="Ih VII XVII"` restricts a run to selected polymorphs, while
`RESULT_ROOT=/path/to/results` keeps results from different executables fully
separate.  The driver records the executable hash, reported version, numerical
settings, structure hashes, JSON results, and process output.  Set
`SKIP_EXISTING=1` to resume an interrupted series without repeating completed
JSON results.

The archived Gamma CP2K inputs used `ACCURACY 0.01`; a direct A/B test showed
that changing `0.1` to `0.01` has no material effect on the tested relative
energies.  Set `ACCURACY=0.01` to reproduce that input choice literally.

For a single phase and mesh, the equivalent command is

```bash
tblite run --method gxtb --acc 0.1 --iterations 300 --no-restart \
  --json result.json structures/k333/VII/POSCAR
```

Add `--grad gradient.txt` to obtain the molecular gradient and virial.

## Absolute-energy comparison

`tables/absolute_energies_vs_mesh.csv` contains absolute total energies in
Hartree.  CP2K-native values refer to one primitive cell.  Direct CLI values
refer to the explicit supercell and are also reported after division by
`N^3`.  The strict comparison is therefore

```text
E(CP2K native, N x N x N) = E(save_tblite CLI, explicit N^3 cell) / N^3.
```

Relative DMC-ICE13 energies and MAEs are deliberately kept in separate tables:
this file is an absolute-energy oracle intended to expose unit, normalization,
cell, occupation, and k-point inconsistencies directly.

## Supplied reference coverage

The current package contains direct `save_tblite` absolute energies for all 13
phases at Gamma, `2 x 2 x 2`, and `3 x 3 x 3`.  CP2K-native absolute energies
are supplied for all phases at `2 x 2 x 2`, `3 x 3 x 3`, and `4 x 4 x 4`.
Independent direct `4 x 4 x 4` supercell gates are supplied for Ih, VII, and
XVII.

Empty CSV fields mean "not yet calculated" and are never interpolated or
replaced by a value from another implementation.

The direct/current and native/current columns agree within
`1.06e-7` Hartree per primitive cell over the complete `2 x 2 x 2` set.  The
complete 13-phase `3 x 3 x 3` set agrees within `1.03e-7` Hartree per primitive
cell (RMS `2.96e-8` Hartree).  The three `4 x 4 x 4` gates agree within
`1.12e-7` Hartree per primitive cell.  This makes the table suitable as a
strict absolute-energy comparison for another `save_tblite` build.

## Derivative validation

`validation/k222_XVII_derivatives` contains a native-versus-explicit-supercell
comparison for energy, forces, and stress, together with independent central
finite differences.  The direct differences are `4.30e-9` Hartree per
primitive cell for the energy, `2.85e-7` Hartree/bohr for the largest force
component, and `0.113` bar for the largest stress component.  Raw inputs,
outputs, gradients, displaced structures, and energies are retained beside the
machine-readable summaries.

`validation/native_derivative_hardening_20260718` contains the final exact
hybrid ACP reverse-path gates: a real-space image contraction for self-inverse
meshes and a direct Bloch-space contraction for genuinely complex meshes.  It
also archives all force/stress ablations and failed representation trials with
SHA-256 provenance.

`validation/provider_revision_bvk_ab_20260718` compares the unchanged Seidler
`pbc` executable and the current `save_tblite` executable directly on the same
explicit `2 x 2 x 2` and `3 x 3 x 3` BvK supercells.  For the difficult ice-VII
gate, the provider revision changes the Ih-referenced energy by less than
`0.1` kJ mol^-1 per water at either mesh; it therefore cannot explain the
earlier large implementation discrepancy or a materially shorter convergence
tail.  An independent Linux repetition covers all 13 phases at `2 x 2 x 2`:
the author-branch and current-provider MAEs differ by only
`+0.02549` kJ mol^-1 per water, with the author branch being marginally worse
at this unconverged mesh.
