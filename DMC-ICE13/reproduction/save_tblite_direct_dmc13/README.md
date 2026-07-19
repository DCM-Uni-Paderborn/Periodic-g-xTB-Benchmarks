# Direct save_tblite reproduction of DMC-ICE13

This package permits an implementation-independent comparison between native
CP2K k-point calculations and direct `save_tblite` calculations on the exactly
commensurate Born--von-Karman supercells.

## Structures

`structures/kNNN/<phase>/POSCAR` and `structure.xyz` contain Cartesian
coordinates in Angstrom.  The extended-XYZ comment line carries the complete
lattice and periodicity, so the geometry can be inspected without interpreting
POSCAR scaling conventions.
For mesh size `N`, the lattice and atom list are explicit `N x N x N`
replications of the primitive DMC-ICE13 cell.  Hence a Gamma-only direct CLI
calculation on that POSCAR is the real-space counterpart of the native,
Gamma-centred `N x N x N` Brillouin-zone mesh.

The 13 phases are `Ih`, `II`, `III`, `IV`, `VI`, `VII`, `VIII`, `IX`, `XI`,
`XIII`, `XIV`, `XV`, and `XVII`.  The structures were generated from the exact
inputs used for the current CP2K-native comparison.  An independent geometric
audit found them equivalent to the supplied `gamma_only_dmc_ice13` archive to
better than `1.2e-9` Angstrom in the periodic pair-distance spectra.

The supplied archive is not an independent direct-CLI result from the model
authors.  Its embedded README identifies it as a Gamma-only exchange package
created on 6 June 2026 from the local CP2K working tree, and its energy table
points back to the corresponding early CP2K output directories.  It is used
here only as an independently packaged geometry source and historical CP2K
snapshot.  `provenance/gamma_exchange_archive_20260606.md` records the archive
identity and the evidence for this classification.

`provenance/structure_hashes.csv` records the SHA-256 digest of every POSCAR,
extended XYZ, and source CP2K input.

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
settings, structure hashes, JSON results, process output, and exit status.
Each phase directory also carries the exact executable and input hashes.  A
nonzero CLI exit or a successful process that produces no JSON is retained as
an explicit failed result.  Set `SKIP_EXISTING=1` to resume an interrupted
series; an existing result is reused only when the executable and input hashes,
exit status, JSON, and normal output markers all qualify.

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

`tables/dmc_reference_relative_energies.csv` gives the DMC reference energies
relative to ice Ih.  `tables/relative_energies_vs_mesh.csv` combines those
references with our absolute CP2K-native and direct-CLI energies at every
available mesh.  It reports the Ih-referenced energy and signed DMC error per
water molecule while retaining the absolute energy and source-file hash in
the same row.  Missing provider/mesh combinations are omitted, never
interpolated.  Regenerate the derived table with
`tools/collect_relative_energies.py`.
`tools/verify_recalculation_package.py` independently checks all Cartesian
structure hashes, atom counts and cell metadata, the complete phase set, both
energy providers, and every relative-energy/error conversion.

Keeping the absolute and relative tables separate makes the former an oracle
for unit, normalization, cell, occupation, and k-point inconsistencies while
the latter is the compact phase-by-phase comparison requested for independent
recalculation.

## Supplied reference coverage

The current package contains direct `save_tblite` absolute energies for all 13
phases at Gamma, `2 x 2 x 2`, and `3 x 3 x 3`.  CP2K-native absolute energies
are supplied for all phases at the same three meshes and at `4 x 4 x 4`.
Independent direct `4 x 4 x 4` supercell gates are supplied for Ih, II, VII,
XI, XIV, and XVII. The final qualified five-phase subset Ih/II/VII/XI/XIV is
also retained with exact executable/input hashes and affinity proofs under
`validation/cli_native_k444_sentinels_20260719`.

Empty CSV fields mean "not yet calculated" and are never interpolated or
replaced by a value from another implementation.

The direct/current and native/current columns agree within `2.97e-8` Hartree
per primitive cell over the complete Gamma set (RMS `1.57e-8` Hartree) and
within `1.06e-7` Hartree over the complete `2 x 2 x 2` set.  The complete
13-phase `3 x 3 x 3` set agrees within `1.03e-7` Hartree per primitive cell
(RMS `2.96e-8` Hartree).  The six `4 x 4 x 4` gates agree within `1.12e-7`
Hartree per primitive cell.  This makes the table suitable as a strict
absolute-energy comparison for another `save_tblite` build.  Run
`tools/verify_absolute_energy_parity.py` to check required coverage, hashes,
normal termination, and the numerical tolerance.  On the Ih-referenced,
per-water scale the maximum native/CLI differences are
`1.89e-6`, `2.14e-5`, and `1.71e-5` kJ mol^-1 for the complete Gamma,
`2 x 2 x 2`, and `3 x 3 x 3` sets, respectively; the independent
`4 x 4 x 4` gates give `1.87e-5` kJ mol^-1.  These differences are many
orders of magnitude below the DMC benchmark errors and cannot explain a
different apparent k-point convergence curve.

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
normal CP2K and direct-CLI termination before it compares an explicit full
native-k mesh, its symmetry-reduced counterpart, and the direct BvK CLI energy
on the common primitive-cell normalization.  Selected CP2K/CLI executable and
input hashes can additionally be made mandatory.  The Gamma-supercell
comparator applies the same termination and optional provenance gates.  The
positive and negative checks in `tools/tests/test_oracle_qualification.py`
ensure that successful parity is accepted while a failed CLI run or changed
CP2K input is rejected before any energy tolerance is considered.

`tools/verify_gamma_cli_requalification.py` independently qualifies a fresh
all-phase Gamma repetition.  It requires a singleton-CPU affinity proof,
controller and phase exit status zero, the requested source revision and
executable hash, exact structure hashes, numerical settings, JSON and
process-completion markers.
Only then does it compare every new absolute energy both with the archived
direct result and CP2K-native, followed by the stricter Ih-referenced,
per-water comparison.  Its positive and deliberately corrupted input, exit,
and affinity cases are covered by
`tools/tests/test_gamma_cli_requalification.py`.

`tools/verify_k222_cli_native_requalification.py` applies the stricter
same-host provenance gate to a fresh all-phase `2 x 2 x 2` repetition.  Every
direct and native phase must carry the selected Linux executable hash, the
exact input hash, normal termination, and a singleton-CPU affinity proof on
both sides.  The native input is independently required to contain the canonical
shifted MacDonald mesh with symmetry reduction enabled.  Only after all
thirteen phases satisfy these conditions are absolute and Ih-referenced
CLI/native differences evaluated.  Positive and deliberately altered binary
and mesh cases are covered by
`tools/tests/test_k222_cli_native_requalification.py`.

`tools/verify_cli_native_mesh_requalification.py` generalizes the same gate
to any cubic MacDonald mesh.  It derives the odd/even mesh shift, verifies the
explicit Born--von--Karman replication and Cartesian coordinates, accepts only
declared singleton CPUs, and compares both total and non-self-consistent
dispersion energies.  It can qualify either all thirteen phases or a declared
sentinel subset that includes ice Ih.  Dedicated positive and negative tests
cover a complete `3 x 3 x 3` fixture, the even-mesh shift rule, altered shifts,
and undeclared CPU assignments.

`tools/verify_no_acp_cli_native.py` supplies a separate component-ablation
gate for ice Ih and XVII on the `2 x 2 x 2` mesh.  It compares the final
CP2K-native energies with direct explicit-BvK CLI energies generated from the
identical ACP-disabled parameter file.  Before any numerical comparison, it
requires normal termination on both sides, exact executable, structure, and
parameter hashes, the selected provider revision, and a singleton-CPU
affinity proof for each CLI calculation.  It also parses both the internal
full-model export and the No-ACP parameter file.  The latter is accepted only
if the global `[acp]` activation table has been removed and every other parsed
parameter is unchanged; the retained elemental projector records are then
inactive by construction.  The gate checks both absolute
primitive-cell energies and the Ih-referenced energy per water.  Positive and
deliberately corrupted provenance, affinity, and energy cases are covered by
`tools/tests/test_no_acp_cli_native.py`.

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

The remaining small final-`pbc`/current-provider difference is resolved by the
top-level `validation/provider_component_attribution_20260719` and
`validation/pbc_h0_anisotropy_attribution_20260719` archives.  A self-consistent
exchange ablation first localized the amplification to the exchange-containing
SCC path.  A controlled source A/B build then proved the actual origin:
restoring the historical central-cell H0 anisotropy reproduces the final-`pbc`
phase-VII `2 x 2 x 2` supercell energy within `1.91e-11` hartree and accounts
for more than 99.999999% of the provider gap.  The current image-complete H0
definition is retained because an equivalent lattice-image representation is
invariant to `2.73e-12` hartree, whereas the historical treatment changes by
`5.60e-8` hartree.

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

`tools/verify_k222_full_reduced_set.py` provides the complementary fresh-build
symmetry gate for all 13 DMC-ICE13 structures.  It accepts only normally
terminated runs from the requested CP2K executable hash, proves singleton CPU
affinity and exact input provenance, and requires the full-grid and reduced
inputs to differ only in their project and symmetry-control lines.  Total and
non-self-consistent dispersion energies are then compared component by
component.  The negative tests reject altered inputs, binaries, affinity,
mesh flags, total energies, and dispersion components.

For a single reproducible implementation audit,
`tools/verify_part_i_implementation.py` executes the complete archived gate
set: absolute CLI/native energy parity, numerical-accuracy sensitivity,
periodic-response A/B, energy/force/stress derivatives, k-point/BvK grid
identity, provider and model revisions, derivative hardening, unchanged-build
sentinels, Wigner--Seitz diagnosis, final-build low-k/partial-PBC derivatives,
exchange/ACP component ablations, periodic-H0 source/invariance attribution,
direct periodic source tests for H0, Wigner--Seitz, exchange, Fock response,
forces, stress, and transform oracles, the complete author-`pbc`/current-CLI/
CP2K-native `3 x 3 x 3` energy closure, the sibling
`../seidler_dmc13_recalculation` author-facing rerun package, and all portable
SHA-256 manifests.  The archived `tools/controllers` suite additionally
enforces the VII--Ih--Gamma-oracle priority, exact executable hashes,
memory-safe serialization, and phase-local `0.10` adaptive pruning in a
synthetic end-to-end dry run.  Its
JSON report records every return code and hashes the script and captured
output of each gate; any failed subordinate verifier fails the aggregate
audit.

For the final adaptive table, `tools/select_adaptive_endpoints.py` enforces the
phase-local one-step rule directly on normally terminated calculations from a
required CP2K executable hash.  It rejects a missing earlier adjacent pair,
selects the denser member of the first pair with an absolute relative-energy
change no larger than the user-selected default threshold of `0.10`
kJ mol^-1 per water molecule, and computes aggregate error
statistics only after all twelve non-reference phases pass.  No aggregate MAE
or RMS condition participates in endpoint selection.

`tools/dmc_phase_convergence.py` and `tools/dmc_mixed_mae.py` apply the same
qualification to incremental decisions and progress reports.  When an
execution hash is required, they also verify normal exit and the hash of the
actual input.  They additionally require the directory mesh, the dimensions
inside the input, and the canonical even/odd MacDonald shift to agree.  A
denser result produced by another executable, from a modified input, or with a
noncanonical shift is skipped instead of masking a valid lower-mesh pair.
`tools/monitor_qualified_mixed_mae.sh` wraps the latter evaluator for a live
calculation campaign.  It writes a result only when all twelve phases have a
normally terminated same-mesh phase/Ih pair from the required executable and
canonical input; otherwise it records an explicit `NOT_READY` reason.  Its
history includes the selected mesh vector and the same-mesh paper comparator,
and `ONCE=1` provides a nonpolling reproducibility gate.
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
