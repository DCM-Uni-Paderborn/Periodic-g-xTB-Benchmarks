# DMC-ICE13 recalculation package for the save_tblite authors

This compact package is intended for an independent rerun with any selected
`save_tblite` executable.  It contains all 13 DMC-ICE13 cells in absolute
Cartesian Angstrom coordinates, the published Ih-referenced DMC energies, and
the absolute energies obtained with the author `pbc` snapshot, the current
CP2K-integration provider, CP2K-native Bloch k points, and the historical
`mstore-inorganic` branch.  `REQUEST_TO_AUTHORS.md` is a ready-to-send request
for an independent two-branch rerun; `EXPECTED_RESULTS.md` gives the acceptance
checks and explains which conclusions follow from each comparison.
`MODEL_SOURCE_DIFFERENCES.md` inventories the model-relevant source history
that separates `mstore-inorganic` from `pbc`.  The reciprocal one-patch
evidence under `evidence/wigner_seitz_self_image_attribution/` further
identifies the Wigner--Seitz self-image-index correction as the dominant
causal source of their different sparse-mesh exchange energies.  The
reciprocal test under `evidence/second_order_mic_attribution/` attributes the
entire remaining source-state residual to the later minimum-image form of the
periodic second-order Coulomb term, within SCC numerical resolution.
The directory `evidence/mstore_inorganic_k444_partial/` adds an independently
qualified historical-source convergence check at `4 x 4 x 4`: Ih and eleven
benchmark phases completed, while phase XIII is retained with its exact
memory-termination record and is excluded from all reported partial metrics.

The author `pbc` snapshot (`c932120...`) and the later pbc-derived integration
provider (`15915c...`) are deliberately separate entries.  CP2K-native parity
must be judged against the latter exact source state; the former is retained
to measure the smaller intervening provider-revision shift.

The primitive structures live under `structures/primitive`.  Explicit
Gamma-centred Born--von Karman cells are generated without changing atom
order by `scripts/build_bvk_from_poscar.py`.  Every actually evaluated current
CLI input from `1 x 1 x 1` through `4 x 4 x 4` is retained under `raw`, checked
byte-for-byte against its recorded input SHA-256, and compared with the
deterministic generator using exact species/order and a maximum cell/coordinate
tolerance of `5e-12` Angstrom.  The historical `2 x 2 x 2` and `3 x 3 x 3`
branch-comparison POSCARs are checked independently against the same generator.
This retains both the exact evaluated inputs and a compact deterministic route
to recreate them.

Run one phase with a chosen executable, for example:

```text
python3 scripts/run_save_tblite.py /path/to/tblite VII 3 results/my-build \
  --require-binary-sha256 <sha256>
```

The command is exactly equivalent to:

```text
tblite run --method gxtb --acc 0.1 --iterations 300 --no-restart \
  --json result.json POSCAR
```

`0.1` is the primary parity setting because it is the value used by the
qualified CP2K-native DMC-ICE13 inputs.  A second run with `--accuracy 0.01`
is an optional tighter SCC-convergence sensitivity check; it must not be mixed
silently into the same-setting parity matrix.  The assembler reparses the
printed energy and density thresholds from every direct-CLI text output and
rejects a parity mesh unless all thirteen runs correspond to `0.1`.

Run the requested complete matrix for one executable sequentially with:

```text
python3 scripts/run_branch_matrix.py /path/to/tblite results/pbc
python3 scripts/run_branch_matrix.py /path/to/tblite results/mstore-inorganic
```

After both branches finish, create the direct return table with:

```text
python3 scripts/summarize_author_results.py \
  results/pbc results/mstore-inorganic results/independent_branch_comparison.csv
```

## Comparison tables

- `tables/cp2k_native_absolute_energies_by_mesh.csv` contains every presently
  completed, hash-qualified CP2K-native absolute energy through `8^3`, the
  exact executed input, and its independently verified settings state.
- `tables/cp2k_native_relative_energies_by_mesh.csv` applies the same-mesh ice
  Ih reference only when both required absolute outputs are qualified.
- `tables/pbc_cli_vs_cp2k_native_absolute_parity.csv` is the direct absolute
  current-`pbc` CLI versus CP2K-native comparison for the complete 52-point
  matrix from `1^3` through `4^3`.
- `tables/current_cli_convergence_provenance.csv` records the convergence
  thresholds printed by every direct-CLI run and infers the effective
  `--acc` value independently of filenames or launcher descriptions.
- `tables/author_pbc_absolute_energies.csv` and its relative companion contain
  the complete author-`pbc` snapshot series at `2^3` and `3^3`.
- `tables/mstore_inorganic_relative_energies_by_mesh.csv` contains the
  independently rebuilt historical `mstore-inorganic` results at Gamma,
  `2^3`, and `3^3`.  A mesh enters
  the statistics table only after all twelve non-reference phases and Ih are
  present with one binary hash.
- `tables/mstore_inorganic_absolute_energies.csv` retains the corresponding
  absolute supercell and primitive-cell energies together with every input,
  output, and executable hash.
- `evidence/mstore_inorganic_k444_partial/` contains the raw historical
  `4 x 4 x 4` outputs for the eleven completed non-reference phases plus Ih,
  absolute and relative-energy tables, the failed XIII record, and a
  self-contained verifier.  Its same-eleven-phase MAE is a convergence
  diagnostic only and is explicitly not a complete DMC-ICE13 statistic.
- `tables/mstore_vs_pbc_relative_differences.csv` is the shortest direct view
  of the branch effect for every phase and common mesh; it also states both
  effective CLI accuracies and whether they are identical for that row.
- `tables/all_branch_relative_energy_comparison.csv` and
  `tables/branch_comparison_statistics.csv` collect the branch-resolved values
  without hiding incomplete meshes.
- `tables/three_route_absolute_energies_k333.csv` and its relative companion
  give the direct current-CLI, author-`pbc`, and CP2K-native `3^3` closure.

The `raw` directory retains the source outputs, exact evaluated CLI POSCARs,
exact executed CP2K inputs, input/binary hashes, and exit states used to
assemble the new tables.  Running `scripts/assemble_comparison_tables.py`
rejects incomplete CP2K outputs, any run without an archived integer exit
status of zero, a non-matching input or executable hash, and inputs that do not
use the qualified `ACCURACY`, SCC mixer, SCF threshold, MacDonald shift, and
SPGLIB-reduced native-k settings.  The generated absolute-energy table exposes
the exit status and normal-termination qualification for every admitted row.
For concatenated output, the last program segment itself must contain the
ordered start marker, final energy, and end marker; a new incomplete segment
cannot inherit an older successful end marker.
`prepare_package.py` independently converts scaled CP2K coordinates where
needed and verifies cell vectors, species order, and Cartesian positions
against the canonical primitive POSCAR for every admitted native endpoint.

The compact `evidence` directory contains the independently reproducible
three-route `3^3` closure, the CP2K native-k versus explicit Gamma-BvK oracle,
the exchange/ACP component ablation, and reciprocal Wigner--Seitz one-patch
and second-order minimum-image one-patch tests, together with the partial
historical `4^3` convergence matrix, including raw inputs, outputs,
verification summaries, and hash manifests.  Historical table rows whose original raw output was not retained
remain explicitly identified by the hashes stored in the tables.

`sources.json` records the exact source states and executable hashes used in
this package.  It also records that the tested `mstore-inorganic` and `pbc`
revisions were still the respective upstream branch tips when checked on
2026-07-20; the literal remote query is archived in
`author_branch_heads_20260720.txt`.  Here, "historical" denotes the older code
lineage rather than a stale local branch.  The `mstore-inorganic` build required only a
dependency-fetch repair for the obsolete `mctc-lib` wrap; the `save_tblite`
source revision itself was not modified.  This diagnostic build is therefore
suitable for locating the model-revision difference, while the requested clean
author builds remain the decisive independent check.

`comparison_workbook.xlsx` provides the same data in a convenient review
workbook; CSV files remain the machine-readable source of truth.
`dmc_ice13_small_mesh_energy_matrix.xlsx` is the compact two-panel view: one
worksheet lists the published DMC cohesive energies beside qualified
CP2K-native, current-`pbc` CLI, and historical `mstore-inorganic` absolute
energies; the second lists their Ih-referenced relative energies.  The matrix
covers Gamma through `4^3` for CP2K-native/current `pbc`, includes the
historical `mstore-inorganic` Gamma, `2^3`, and `3^3` results, and leaves only
the uncomputed historical `mstore-inorganic` `4^3` column explicitly blank.
`comparison_summary.json` gives the compact machine-readable conclusion,
source states, numerical parity checks, and complete-mesh statistics.

All generated files are covered by `SHA256SUMS`.  Rebuild and verify the
package from its authoritative parent archive with:

```text
python3 scripts/assemble_comparison_tables.py
python3 prepare_package.py --refresh-manifest-only
shasum -a 256 -c SHA256SUMS
```
