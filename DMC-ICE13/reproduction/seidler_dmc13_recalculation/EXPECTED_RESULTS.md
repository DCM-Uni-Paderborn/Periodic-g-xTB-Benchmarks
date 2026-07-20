# Expected checkpoints and interpretation

This file separates three questions that should not be conflated.

## 1. Does CP2K-native reproduce `pbc`?

Use the current pbc-derived integration source state `15915c943564...`
recorded in `sources.json`, not the older `c932120d258...` author snapshot.
The direct
`save_tblite` executable evaluates explicit Gamma-centred Born--von Karman
supercells, whereas CP2K evaluates native Bloch k points with symmetry
reduction.  All 52 required points from `1 x 1 x 1` through `4 x 4 x 4`
pass.  Across the complete matrix, the largest absolute primitive-cell
difference is `1.12e-7 Eh`; after referencing every phase to same-mesh ice
Ih, the largest difference is `2.14e-5 kJ mol-1 H2O-1`.  These are numerical
differences, not a benchmark-accuracy difference.

The primary parity matrix uses `--acc 0.1`, exactly matching the qualified
CP2K inputs.  The effective value is not inferred from a directory name: the
assembler checks the printed `1e-7 Eh` energy and `2e-6 e` density thresholds
for every direct run.  A separately labelled `--acc 0.01` calculation is a
useful tighter convergence test, but it is not substituted for the
same-setting comparison.

The primary files are:

- `tables/three_route_absolute_energies_k333.csv`;
- `tables/three_route_relative_energies_k333.csv`;
- `tables/pbc_cli_vs_cp2k_native_absolute_parity.csv`.

The repository-level gate
`validation/binary_provider_identity_20260720/` additionally verifies that
the direct CLI link rule and the qualified CP2K build use the same static
provider archive, not merely source trees bearing the same revision label.
The independent
`validation/relative_energy_postprocessing_20260720/` gate reconstructs the
BvK normalization, per-water same-mesh Ih referencing, unit conversion, DMC
errors, and aggregate statistics directly from the raw absolute outputs.
The source-level and exact-arithmetic gate
`validation/macdonald_bvk_mesh_equivalence_20260720/` proves that the native
MacDonald grids in all archived inputs are precisely the reciprocal-space
folding grids of the corresponding Gamma-only BvK supercells before symmetry
reduction.  Thus the direct/native comparison does not contrast different
twists or different Brillouin-zone meshes.

An independent clean build of the exact `15915c...` source should reproduce
the direct CLI column before CP2K is involved.  A deviation larger than approximately `1e-7 Eh`
per primitive cell should be investigated from the absolute energy and exact
input, not from a relative MAE.

## 2. Are `mstore-inorganic` and `pbc` the same model state?

No.  The supplied complete `2^3` and `3^3` diagnostic calculations show
phase-dependent differences far above numerical noise.  At `3^3`, the MAEs
are `17.8306` for `mstore-inorganic` and `34.0485 kJ mol-1 H2O-1` for the
current pbc-derived source, and the largest individual relative-energy shift
between those branches is `76.4515 kJ mol-1 H2O-1`.  The most direct view is
`tables/mstore_vs_pbc_relative_differences.csv`; the underlying absolute
`mstore-inorganic` energies are in
`tables/mstore_inorganic_absolute_energies.csv`.

The complete `2^3` and `3^3` comparisons are now strictly same-setting at CLI
accuracy `0.1`.  A separate independently rebuilt `mstore-inorganic` `3^3`
sensitivity matrix at `0.01` differs from the `0.1` matrix by at most
`2.41e-12 Eh` per explicit supercell.  Controlled pbc-derived runs likewise
show a sub-`2e-10 Eh` sensitivity over that interval.  The tens-of-kJ/mol
branch shift is therefore not an SCC-accuracy artifact; an author rerun is
still requested as an independent provenance check rather than as closure of
a missing common-setting comparison.

The separate `c932120...` `pbc` snapshot is retained in the author-pbc tables.
It differs slightly from `15915c...`; that same-source distinction must not be
misclassified as a CP2K interface error.  Across the twelve relative energies,
the largest `c932120...` versus `15915c...` shift is
`0.1355 kJ mol-1 H2O-1` at `2^3` and `0.1005 kJ mol-1 H2O-1` at `3^3`.
Those values are far larger than the current-CLI/CP2K numerical residual, but
far smaller than the historical `mstore-inorganic`/`pbc` model shift.

The `mstore-inorganic` numbers were produced from the historical source
revision in `sources.json`.  Its obsolete dependency locator was repaired to
permit compilation, but no `save_tblite` source file was changed.  We ask the
authors to repeat this comparison with their own clean builds so that branch
history, dependency state, and executable provenance are independently
controlled.

## 3. Which source state produced the previously quoted lower DMC error?

This is still open.  A lower error from `mstore-inorganic` at a sparse mesh
does not establish its converged accuracy, and neither sparse-mesh value is a
substitute for the ongoing dense-grid CP2K-native series.  The requested
author rerun should identify the exact source revision, executable hash,
structures, `--acc` value, mesh convention, and post-processing used for the
previously quoted result.

The decisive diagnostic sequence is therefore:

1. reproduce the supplied absolute `pbc` CLI energies;
2. verify that the independently built `pbc` CLI agrees with CP2K-native;
3. rerun the identical structures with `mstore-inorganic`;
4. identify which branch and settings generated the earlier benchmark data.

All relative energies in this package use ice Ih from the same model and the
same mesh.  No value is compared against an Ih energy from another mesh.
