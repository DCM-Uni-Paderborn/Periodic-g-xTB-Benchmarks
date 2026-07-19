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

The archived Gamma CP2K inputs used `ACCURACY 0.01`.
`validation/accuracy_sensitivity_20260718` contains the complete direct and
CP2K-native A/B gate: changing `0.1` to `0.01` alters the tested relative
energies by less than `6e-8` kJ mol^-1 per water.  Set `ACCURACY=0.01` to
reproduce that input choice literally.

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
are supplied for all phases at the same three meshes and at `4 x 4 x 4`.
Independent direct `4 x 4 x 4` supercell gates are supplied for Ih, VII, and
XVII.

Empty CSV fields mean "not yet calculated" and are never interpolated or
replaced by a value from another implementation.

The direct/current and native/current columns agree within `2.97e-8` Hartree
per primitive cell over the complete Gamma set (RMS `1.57e-8` Hartree) and
within `1.06e-7` Hartree over the complete `2 x 2 x 2` set.  The complete
13-phase `3 x 3 x 3` set agrees within `1.03e-7` Hartree per primitive cell
(RMS `2.96e-8` Hartree).  The three `4 x 4 x 4` gates agree within `1.12e-7`
Hartree per primitive cell.  This makes the table suitable as a strict
absolute-energy comparison for another `save_tblite` build.  Run
`tools/verify_absolute_energy_parity.py` to check required coverage, hashes,
normal termination, and the numerical tolerance.

`validation/kpoint_grid_bvk_oracle_20260718` independently checks the point
sets themselves.  The complete printed `2 x 2 x 2`, `3 x 3 x 3`, and
`4 x 4 x 4` native grids are equal-weight realizations of the corresponding
Gamma-centred BvK reciprocal grids, including the even-mesh MacDonald shifts.
Thus the direct-supercell comparison does not hide a shifted-mesh mismatch.

For the explicit CP2K Gamma-supercell oracle, use
`tools/build_cp2k_gamma_supercell.py` to generate the input and
`tools/verify_gamma_supercell_input.py` to compare its cell, atom ordering,
and every Cartesian coordinate against the archived POSCAR.  The verifier
also rejects a non-XYZ cell or any explicit KPOINTS section, so the energy
gate cannot silently compare different structures or boundary conditions.
`tools/compare_native_symmetry_cli.py` complements this input gate by requiring
normal termination and a selected CP2K executable hash before it compares an
explicit full native-k mesh, its symmetry-reduced counterpart, and the direct
BvK CLI energy on the common primitive-cell normalization.

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
at this unconverged mesh.  The complete all-phase `3 x 3 x 3` repetition gives
the same conclusion: the corresponding MAE change is `+0.02210` kJ mol^-1 per
water, and the largest phase-resolved relative-energy shift is `0.10048`
kJ mol^-1 per water for ice VI.

`validation/model_revision_coarse_grid_ab_20260718` extends that audit to the
separate `mstore-inorganic`, post-March molecular g-xTB, and DCM `main` source
histories.  On identical Cartesian structures, these histories produce very
different Gamma and `2 x 2 x 2` DMC-ICE13 MAEs.  The current integration and
DCM `main` remain close, while the two historical development lines yield much
smaller coarse-grid MAEs.  This establishes a source-revision effect, not a
CP2K-native k-point discrepancy.  The complete `3 x 3 x 3` extension sharpens
that conclusion: current integration and final author `pbc` give MAEs of
`34.04849` and `34.07059` kJ mol^-1 per water, whereas the obsolete
`mstore-inorganic` history gives `17.83062` kJ mol^-1 per water.

`validation/wigner_seitz_branch_diagnosis_20260718` resolves the dominant
source-revision effect.  The older `mstore-inorganic` history retained a
compact-position Wigner--Seitz self-image mapping that fails the cubic,
orthorhombic, and skew-cell invariants.  Leopold Seidler fixed that mapping in
the newer, final `pbc` branch.  A controlled two-line A/B build reproduces the
large coarse-grid shift, whereas changing only the image-distance threshold
does not.  The apparently improved `mstore-inorganic` MAE is therefore a
legacy indexing artifact and must not be used as the periodic author
reference.

`validation/cp2k_response_fix_ab_20260719` closes the independent CP2K-side
audit.  Its complete 13-structure `2 x 2 x 2` gate compares native,
SPGLIB-reduced primitive-cell energies with explicit standalone-CLI BvK
supercells and finds a maximum Ih-referenced difference of only
`2.1414e-5` kJ mol^-1 per water.  Controlled before/after and no-ACP tests
assign the earlier CP2K energy shift predominantly to the corrected periodic
BvK Coulomb response rather than to smearing, symmetry acceleration, or the
provider frontend.  The package also contains focused GENERAL-grid and
density-mixer restart gates plus the complete final CP2K regression result:
78 correct, zero wrong, and zero failed checks.  Run its
`verify_response_fix.py` script to recompute the tables and verify the full
SHA-256 manifest.

`validation/qualified_energy_sentinels_20260719` then compares that
response-corrected reference build with the subsequently hardened source on
ice Ih and ice VII at `2 x 2 x 2` and on ice VII at `3 x 3 x 3`. All three
qualified-minus-reference energies are exactly zero at the precision printed
by CP2K. This proves that the regular-mesh/restart hardening is inactive for
the frozen DMC energy inputs; previously completed response-corrected points
therefore do not require a blanket rerun merely because of that hardening.

For the final adaptive table, `tools/select_adaptive_endpoints.py` enforces the
phase-local one-step rule directly on normally terminated calculations from a
required CP2K executable hash.  It rejects a missing earlier adjacent pair,
selects the denser member of the first pair with an absolute relative-energy
change no larger than the chosen threshold, and computes aggregate error
statistics only after all twelve non-reference phases pass.  No aggregate MAE
or RMS condition participates in endpoint selection.

`tools/dmc_phase_convergence.py` and `tools/dmc_mixed_mae.py` apply the same
qualification to incremental decisions and progress reports.  When an
execution hash is required, they also verify normal exit and the hash of the
actual input.  They additionally require the directory mesh, the dimensions
inside the input, and the canonical even/odd MacDonald shift to agree.  A
denser result produced by another executable, from a modified input, or with a
noncanonical shift is skipped instead of masking a valid lower-mesh pair.
`tools/verify_adaptive_dmc13.py` is an independent final oracle: it reparses
the raw CP2K energies and MacDonald meshes, verifies binary and input
provenance, proves that every reported endpoint is the first passing adjacent
pair, and recomputes all aggregate statistics.  The positive and negative
end-to-end checks in `tools/tests/test_adaptive_reporting.py` can be run with
`python3 -m unittest -v tools/tests/test_adaptive_reporting.py`.
All portable integrity manifests, including their relative paths and complete
file coverage, can be checked after a fresh clone with
`python3 tools/verify_sha256_manifests.py`.  The separately named
`SOURCE_SHA256SUMS` files preserve historical full-worktree inventories and
are intentionally not interpreted as manifests for files shipped in the
compact archive.

If an unresolved phase needs a denser regular mesh,
`tools/build_native_mesh_input.py` rewrites exactly the single cubic
`SCHEME MACDONALD` line of a frozen input.  Along with the mesh dimensions it
recomputes the canonical Gamma-centred BvK MacDonald shift: zero for odd meshes
and `(N-1)/(2N)` for even meshes.  It refuses anisotropic or ambiguous sources
and records both input hashes, both shifts, and the sole changed line,
preventing a dynamic extension from silently inheriting the wrong even/odd
shift or changing any other numerical setting.
