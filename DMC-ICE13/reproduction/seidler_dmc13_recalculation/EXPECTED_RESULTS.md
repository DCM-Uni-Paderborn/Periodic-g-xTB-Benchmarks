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
`6.5e-11 Eh` per explicit supercell and `2.26e-10 kJ mol-1 H2O-1` after
same-mesh Ih referencing.  Controlled pbc-derived runs likewise
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

The partial historical `4^3` extension further shows that the mean absolute
branch gap on the same eleven completed phases shrinks from `49.2110` at
`2^3` to `18.0207` at `3^3` and `5.2456 kJ mol-1 H2O-1` at `4^3`.  Its
historical same-subset MAE is `7.0043`, but the mean absolute
`3^3`-to-`4^3` change is still `12.6064 kJ mol-1 H2O-1`.  Because XIII was
terminated before the first SCC result, `evidence/mstore_inorganic_k444_partial/`
sets `full_matrix_complete` and `usable_for_full_benchmark_statistics` to
false.  These data diagnose finite-size convergence; they must not be quoted
as a complete `4^3` DMC-ICE13 MAE.

The large historical branch shift is now component-classified rather than
merely observed.  For ice VII relative to same-mesh Ih at `2^3`, the full
author-`pbc` minus `mstore-inorganic` gap is
`-148.1194 kJ mol-1 H2O-1`.  Disabling exchange reduces its magnitude by
`98.57%`, while disabling ACP alone reduces it by only `4.28%`.  The
corresponding exact-source, self-consistent runs are archived in
`evidence/mstore_pbc_component_ablation/`.  This identifies exchange-path
changes between the branches as the origin of their different sparse-mesh
behavior; it does not implicate the CP2K interface.

Reciprocal one-patch builds identify the dominant exchange-path change more
precisely.  The `30b04691e0af` Wigner--Seitz fix preserves the original
translation index after the zero-distance origin is removed from a compacted
distance list.  Switching only that routine in the `pbc` build shifts the
ice-VII-minus-Ih `2^3` energy by `-153.5880 kJ mol-1 H2O-1`; switching only the
same routine in the historical `mstore-inorganic` build shifts it by
`-153.9405 kJ mol-1 H2O-1`.  The reciprocal shifts agree within
`0.3524 kJ mol-1 H2O-1` and explain more than 95% of the full branch gap.
The corrected same-build `pbc` executable reproduces the archived author-`pbc`
Ih and VII absolute energies within `2e-12 Eh`.  Full raw evidence and the
verifier are under `evidence/wigner_seitz_self_image_attribution/`.

The smaller post-Wigner--Seitz residual is now independently source-attributed
as well.  Starting from the author `pbc` snapshot and reverting only the later
minimum-image form of the periodic second-order Coulomb term changes the same
ice-VII-minus-Ih `2^3` relative energy from `-300.0673203` to
`-305.8884691 kJ mol-1 H2O-1`.  The independently Wigner--Seitz-corrected
`mstore-inorganic` result is `-305.8884581 kJ mol-1 H2O-1`; their difference is
only `1.10e-5 kJ mol-1 H2O-1`.  This reciprocal one-patch experiment explains
`99.99981%` of the entire `5.82114 kJ mol-1 H2O-1` post-Wigner--Seitz residual.
The raw inputs, converged outputs, inverse source patch, matched build options,
and verifier are under `evidence/second_order_mic_attribution/`.

The same three source states were also evaluated for the smaller ice-XVII
cell.  There the post-Wigner--Seitz residual is
`0.16420623 kJ mol-1 H2O-1`; removing only the minimum-image second-order
variant leaves `7.90e-7 kJ mol-1 H2O-1` and independently explains
`99.99952%`.  This rules out an ice-VII-specific cancellation.

Together, the two reciprocal source-patch experiments classify the complete
historical branch separation at this diagnostic point: the dominant term is
the corrected periodic exchange self-image mapping, and the remaining term is
the later minimum-image second-order Coulomb form.  Neither is generated by
the CP2K-native interface.

The `mstore-inorganic` numbers were produced from the source revision in
`sources.json`.  A direct remote-head check on 2026-07-20 confirmed that this
revision is still the tip of `lmseidler/save_tblite:mstore-inorganic`; the word
"historical" elsewhere in this package describes its older code lineage, not
a stale local branch or an invented reconstruction.  The same check identifies
`c932120...` as the tip of `lmseidler/save_tblite:pbc`.  The obsolete dependency
locator of the former was repaired to permit compilation, but no `save_tblite`
source file was changed.  We ask the authors to repeat this comparison with
their own clean builds so that dependency state and executable provenance are
independently controlled.

## 3. Which source state produced the previously quoted lower DMC error?

The lower sparse-mesh trend is now consistent with the historical
`mstore-inorganic` Wigner--Seitz and second-order Coulomb behavior, but the exact provenance of the
previously quoted benchmark number is still open.  A lower error from that
source state does not establish its converged accuracy, and neither
sparse-mesh value substitutes for the ongoing dense-grid CP2K-native series.
The requested author rerun should identify the exact source revision,
executable hash, structures, `--acc` value, mesh convention, and
post-processing used for the quoted result.

The decisive diagnostic sequence is therefore:

1. reproduce the supplied absolute `pbc` CLI energies;
2. verify that the independently built `pbc` CLI agrees with CP2K-native;
3. rerun the identical structures with `mstore-inorganic`;
4. identify which branch and settings generated the earlier benchmark data.

All relative energies in this package use ice Ih from the same model and the
same mesh.  No value is compared against an Ih energy from another mesh.
