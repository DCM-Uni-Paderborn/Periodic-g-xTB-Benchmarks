# X23b Periodic GFN Benchmark

This directory contains CP2K/tblite calculations for the X23b molecular-crystal
benchmark of Dolgonos, Hoja, and Boese. The reference lattice energies are the
recommended experimental back-corrected values from Table 5, and the reference
cell volumes are the electronic reference volumes from Table 2 of that work.
The primary relaxed-cell benchmark in the current manuscript revision uses
native Bloch 2x2x2 CP2K `&KPOINTS` cell optimizations with full SPGLIB
symmetry reduction. Reported lattice energies are reevaluated on those final
geometries with a Gamma-centered 3x3x3 mesh. This is neither a Gamma-only cell
optimization nor a Born-von-Karman supercell calculation.

The crystal structures are taken from the open X23 `refdata` set. Hexamine is
the only special case: the open experimental CIF contains only heavy atoms, so
the complete X23 Quantum ESPRESSO crystal input is used for that system.

## Contents

- `structures/`: P1 CIF crystal structures and gas-phase molecular starting
  geometries.
- `inputs/`: CP2K input files for crystal single points, gas-phase molecular
  optimizations, and retained Gamma-point crystal cell optimizations.
- `runs/`: generated CP2K working directories, ignored by Git.
- `data/`: metadata, reference values, extracted energies, volume errors, and
  aggregate statistics, including the DMC-X23 comparison values used for the
  system-resolved lattice-energy figure. The final-geometry mesh convergence
  is retained in `x23b_final_geometry_kpoint_{rows,summary}.csv`. The
  exact source, patch, executable, and protocol record is
  `data/build_provenance.json`; the
  `x23b_reference_cli_{gfn1,gfn2}_{rows,summary}.csv` files record direct
  CP2K-native versus tblite CLI checks on all initial and final Gamma
  geometries.
- `figures/`: PDF, SVG, and PNG versions of the three X23b plots used in the
  revised manuscript and Supporting Information.
- `scripts/`: input generation, analysis, plotting, and run scripts.

## Run Defaults

The run script expects the CP2K executable through the `CP2K` environment
variable, or otherwise falls back to `cp2k.psmp`. The default execution mode is
many independent single-core jobs:

- `OMP_NUM_THREADS=1`
- `OPENBLAS_NUM_THREADS=1`
- `MKL_NUM_THREADS=1`
- `VECLIB_MAXIMUM_THREADS=1`
- `OMP_WAIT_POLICY=PASSIVE`
- `CP2K_PARALLEL_JOBS=20`

This was faster for the small DMC-ICE13 and X23b-style xTB jobs than hybrid
MPI/OpenMP execution.

## Current primary result

The final 23/23 converged X23b relaxed-cell data per method are stored in
`data/x23b_lattice_energies.csv`, `data/x23b_cell_volumes.csv`, and
`data/x23b_summary.csv`. The volume rows are `cell_opt,k222`; the manuscript
lattice-energy rows are `cell_opt_single_point,k333`. The raw k222 energies
and the k444 convergence checks remain in the same files.

| Quantity | Method | ME | MAE | RMSE | MaxAE |
|---|---|---:|---:|---:|---:|
| Lattice energy / kJ mol-1, k333 on k222 geometry | GFN1-xTB | 0.258871 | 11.345702 | 14.019344 | 30.935058 |
| Lattice energy / kJ mol-1, k333 on k222 geometry | GFN2-xTB | -12.018989 | 14.092104 | 21.341752 | 77.785392 |
| Cell volume / percent, k222 optimization | GFN1-xTB | -5.960071 | 7.514116 | 9.019708 | 19.236681 |
| Cell volume / percent, k222 optimization | GFN2-xTB | -1.657324 | 5.842296 | 7.530373 | 19.952589 |

The k333-to-k444 mean absolute energy changes on the final geometries are
0.079329 kJ mol-1 for GFN1-xTB and 0.084265 kJ mol-1 for GFN2-xTB. The
fixed-reference-geometry single-point rows remain as a separate diagnostic.

Recreate the curated tables and figures after collecting the cell
optimizations and final-geometry single points:

```bash
python3 X23b/scripts/x23b_pipeline.py analyse \
  --cellopt-csv /path/to/x23b_k222_cellopt_results.csv \
  --final-kpoint-csv X23b/data/x23b_final_geometry_kpoint_rows.csv
```

## Adaptive all-method k-point convergence

`scripts/x23b_adaptive_kpoint_convergence.py` defines a separate, exact
23-crystal convergence campaign for `GFN1`, `GFN2`, and `GXTB`.  All three
methods are executed by the same frozen CP2K executable and build manifest.
For every method/crystal pair it evaluates two independent native-Bloch
series, always with SPGLIB symmetry reduction and `FULL_GRID F`:

- fixed experimental-cell single points, starting with k111 and k222; and
- full `CELL_OPT` calculations independently restarted from the experimental
  crystal at each mesh, starting with k111, k222, and k333.

One consecutive passing step is sufficient and the denser result is retained.
The single-point threshold is an absolute cohesive-energy change of at most
0.05 kJ/mol.  A cell step passes only when the maximum relative change of the
three lengths is at most 0.05%, the relative volume change is at most 0.10%,
and the maximum absolute angle change is at most 0.05 degrees.  Failed tracks
advance one mesh at a time through k888.  The workflow stores every raw delta
and every individual pass flag.

Production is deliberately held until a new manifest is both
`production_ready` and identifies the final post-5582 build.  The currently
archived pre-post-5582 manifest is rejected even though it is marked
`production_ready`.  Once the final build exists, use one empty output root:

```bash
SCRIPT=X23b/scripts/x23b_adaptive_kpoint_convergence.py
KCONV_ROOT=/path/to/x23b-adaptive-post5582
CP2K=/path/to/final-post5582/bin/cp2k.psmp
CP2K_SOURCE=/path/to/final-post5582/cp2k
SAVE_TBLITE=/path/to/final-post5582/bin/tblite
SAVE_TBLITE_SOURCE=/path/to/final-post5582/save_tblite
CAMPAIGN_MANIFEST=campaigns/gxtb-pbc-v1-post5582-YYYYMMDD/build_manifest.json

python3 "$SCRIPT" prepare --output-root "$KCONV_ROOT" \
  --cp2k "$CP2K" --cp2k-source "$CP2K_SOURCE" \
  --save-tblite "$SAVE_TBLITE" --save-tblite-source "$SAVE_TBLITE_SOURCE" \
  --campaign-manifest "$CAMPAIGN_MANIFEST"

# Copy the second line printed by prepare; do not recompute it after editing files.
WORKFLOW_SHA256=<exact-sha256-printed-by-prepare>
python3 "$SCRIPT" run \
  --workflow-manifest "$KCONV_ROOT/workflow_manifest.json" \
  --workflow-manifest-sha256 "$WORKFLOW_SHA256" \
  --cp2k "$CP2K" --cp2k-source "$CP2K_SOURCE" \
  --save-tblite "$SAVE_TBLITE" --save-tblite-source "$SAVE_TBLITE_SOURCE" \
  --campaign-manifest "$CAMPAIGN_MANIFEST" \
  --jobs 23 --threads-per-job 1

python3 "$SCRIPT" finalize \
  --workflow-manifest "$KCONV_ROOT/workflow_manifest.json" \
  --workflow-manifest-sha256 "$WORKFLOW_SHA256" \
  --output-dir X23b/data/kpoint_convergence \
  --cp2k "$CP2K" --cp2k-source "$CP2K_SOURCE" \
  --save-tblite "$SAVE_TBLITE" --save-tblite-source "$SAVE_TBLITE_SOURCE" \
  --campaign-manifest "$CAMPAIGN_MANIFEST"
```

Preparation freezes all 1104 possible inputs and their exact order, but `run`
submits only the current adaptive frontier.  Every completed result is bound
to the input, output, direct CP2K command, thread count, build/campaign,
experimental CIF, generator, and workflow-manifest hashes.  Partial, stale,
duplicated, reordered, post-convergence, or foreign results abort collection.
The final CSV, JSON, and TeX artifacts are replaced atomically only after all
138 method/crystal/series tracks have converged.

## Additive g-xTB production workflow

`GXTB` is an additive third method.  Its gas molecule, Gamma preoptimization,
k222 cell, and k333/k444 energies are all method-owned; no GFN1/GFN2 molecular
energy or geometry is accepted as a source.  The generated CP2K inputs contain
`METHOD GXTB` and exactly `SCC_MIXER TBLITE`, selecting save_tblite's native
complete-Fock potential DIIS.  Production inputs contain no alternative mixer
override. All native-Bloch production meshes use the same full SPGLIB symmetry
reduction as GFN1/GFN2 (`SYMMETRY T`, `FULL_GRID F`, and the SPGLIB backend
and reduction method). CP2K expands the irreducible matrices to the complete
mesh for the coupled save_tblite evaluation and folds the response back.

The V1 campaign paths used below are:

```bash
CP2K=/tmp/cp2k_gxtb_18d37c_bench_build/bin/cp2k.psmp
CP2K_SOURCE=/tmp/cp2k_gxtb
SAVE_TBLITE=/tmp/save_tblite_1449feb_bench_install/bin/tblite
SAVE_TBLITE_SOURCE=/tmp/save_tblite_cp2k
CAMPAIGN_MANIFEST=campaigns/gxtb-pbc-v1-20260714/build_manifest.json
```

### Required shifted-k222 force/stress gate

Before changing the campaign state to `production_ready`, run the versioned
Ammonia comparison in `validation/gxtb_k222_force_stress`. Both inputs use the
frozen X23b reference structure and identical `GXTB`, native TBLITE mixer, SCF,
shifted `2x2x2`, force-print, and analytical-stress settings. The diagnostic
reference uses all eight points (`SYMMETRY F`, `FULL_GRID T`); the candidate
uses the production SPGLIB reduction. This full-grid exception is isolated
from and cannot be collected by the production workflow.

```bash
GATE=X23b/scripts/x23b_k222_force_stress_gate.py
python3 "$GATE" prepare \
  --cp2k "$CP2K" --cp2k-source "$CP2K_SOURCE" \
  --save-tblite "$SAVE_TBLITE" --save-tblite-source "$SAVE_TBLITE_SOURCE" \
  --campaign-manifest "$CAMPAIGN_MANIFEST"
python3 "$GATE" run --variant full \
  --cp2k "$CP2K" --cp2k-source "$CP2K_SOURCE" \
  --save-tblite "$SAVE_TBLITE" --save-tblite-source "$SAVE_TBLITE_SOURCE" \
  --campaign-manifest "$CAMPAIGN_MANIFEST"
python3 "$GATE" run --variant spglib \
  --cp2k "$CP2K" --cp2k-source "$CP2K_SOURCE" \
  --save-tblite "$SAVE_TBLITE" --save-tblite-source "$SAVE_TBLITE_SOURCE" \
  --campaign-manifest "$CAMPAIGN_MANIFEST"
python3 "$GATE" check --campaign-manifest "$CAMPAIGN_MANIFEST"
```

Each `run` command revalidates the exact launcher, loaded `libcp2k`, CP2K and
save_tblite source revisions, save_tblite CLI, and static `libtblite.a`, then
stamps the input and output hashes. `check` derives the full count of eight
from the hashed `MACDONALD 2 2 2` plus `SYMMETRY F`/`FULL_GRID T` contract.
The SPGLIB output must report fewer than eight special points and all eight
rows of its explicit `2x2x2` mesh mapping. It passes only for
`|Delta E| <= 1e-9 Hartree`, maximum force-component difference
`<= 1e-6 Hartree/Bohr`, and maximum stress-component difference
`<= 1e-5 GPa`; RMS differences are also recorded. The machine-readable result
is `X23b/runs/validation/gxtb_k222_force_stress_v1_final/gate_result.json`. The script
accepts `validation_in_progress` only for this pre-production gate and never
promotes the central manifest itself.

Once that result passes and the central manifest records the validation and is
explicitly promoted to `production_ready`, prepare and run the 23 gas jobs.
The previous GXTB Gamma CELL_OPT outputs are quarantined diagnostics and are
not a source for the experimental-reference V1 path:

```bash
python3 X23b/scripts/x23b_pipeline.py prepare --method GXTB --production-only
python3 X23b/scripts/x23b_pipeline.py run --method GXTB \
  --phase molecule_geoopt --cp2k "$CP2K" --jobs 6 --threads-per-job 1 \
  --cp2k-source "$CP2K_SOURCE" --save-tblite "$SAVE_TBLITE" \
  --save-tblite-source "$SAVE_TBLITE_SOURCE" \
  --campaign-manifest "$CAMPAIGN_MANIFEST" \
  --prune-transients
```

The default `--campaign-manifest` is the V1 freeze above and is the single
source of truth. Production is blocked until its state is `production_ready`.
The gate requires clean CP2K and save_tblite source trees and verifies the exact
CP2K launcher, loaded `libcp2k`, embedded/full CP2K revisions, save_tblite CLI,
static `libtblite.a`, CMake fingerprints, and fetched-dependency lock. A
mismatch is rejected rather than stamped; no runtime-link inference is made
for the static archive.

Before any cell optimization, measure energy, forces, and analytical stress on
all 23 frozen X23 reference structures with the shifted k222 SPGLIB contract.
This additive phase does not start a CELL_OPT and does not touch any Gamma run:

```bash
PREFLIGHT_ROOT=X23b/runs/gxtb_native/experimental_k222_preflight
PREFLIGHT=X23b/scripts/x23b_experimental_k222_preflight.py
python3 "$PREFLIGHT" prepare --output-root "$PREFLIGHT_ROOT"
python3 "$PREFLIGHT" run --output-root "$PREFLIGHT_ROOT" \
  --cp2k "$CP2K" --cp2k-source "$CP2K_SOURCE" \
  --save-tblite "$SAVE_TBLITE" --save-tblite-source "$SAVE_TBLITE_SOURCE" \
  --campaign-manifest "$CAMPAIGN_MANIFEST" \
  --jobs 6 --threads-per-job 1
python3 "$PREFLIGHT" collect --output-root "$PREFLIGHT_ROOT" \
  --csv X23b/data/gxtb_staging/x23b_experimental_k222_preflight.csv
```

The preflight manifest contains all 23 source-input, structure, generated-input,
campaign hashes, and output paths. Completed job stamps and the collector add
the output hashes. The collector records energy, maximum and RMS
force, all nine stress components, maximum absolute stress, pressure, and start
volume. A fully parsed run is marked `measured_not_approved`: large finite
derivatives remain reportable scientific findings and are never silently
converted into an approval.

Before using those derivatives in a cell optimization, run the additive finite-
difference pilot on the exact repository IDs `ammonia`,
`14-cyclohexanedione`, `acetic_acid`, and `ethylcarbamate`. The last system has
a fully oblique triclinic reference cell and is the deliberately low-symmetry
case. The default pilot has only 36 shifted-k222 calculations: one analytical
force/stress baseline, two central collective-coordinate pairs, one central
isotropic-strain pair, and one central symmetric-xy-shear pair per crystal.
The collective Cartesian directions are deterministic, translation-free, and
normalized over all `3N` components; this is not a full `6N` coordinate scan.

```bash
FD_ROOT=X23b/runs/gxtb_native/k222_fd_gate_v1
FD_GATE=X23b/scripts/x23b_k222_fd_gate.py
python3 "$FD_GATE" prepare --output-root "$FD_ROOT" \
  --coordinate-step-bohr 1.0e-3 --strain-step 5.0e-4 \
  --coordinate-directions 2
python3 "$FD_GATE" run --output-root "$FD_ROOT" \
  --cp2k "$CP2K" --cp2k-source "$CP2K_SOURCE" \
  --save-tblite "$SAVE_TBLITE" --save-tblite-source "$SAVE_TBLITE_SOURCE" \
  --campaign-manifest "$CAMPAIGN_MANIFEST" \
  --jobs 2 --threads-per-job 1
python3 "$FD_GATE" collect --output-root "$FD_ROOT" \
  --csv X23b/data/gxtb_staging/x23b_k222_fd_measured.csv \
  --json X23b/data/gxtb_staging/x23b_k222_fd_measured.json
```

For a normalized Cartesian direction `d`, the report compares the central
finite difference with `-sum_i F_i dot d_i`. For a symmetric strain generator
`G`, it compares `dE/dh` with
`-V (sigma:G) * 2.2937122783963248e-4 Hartree/(GPa Angstrom^3)`, including the
CP2K stress sign, the reference volume, the full double contraction, and the
unit factor as explicit report fields. Preparation freezes the steps,
directions, campaign/build identity, source structures, generated inputs, and
all hashes in an immutable manifest. A matching successful stamp is the only
resumable state; foreign partial output or a changed input is reported as
`STALE_OUTPUT`/`STALE_STAMP` and is never overwritten.

Collection writes only `measured_not_approved` CSV/JSON artifacts. Scientific
approval is a separate immutable JSON and requires a named reviewer plus two
explicitly chosen tolerances; the measurement files are not modified:

```bash
: "${COORD_TOL_HA_PER_BOHR:?set reviewed coordinate tolerance}"
: "${STRESS_TOL_GPA:?set reviewed stress-conjugation tolerance}"
python3 "$FD_GATE" approve \
  --report-json X23b/data/gxtb_staging/x23b_k222_fd_measured.json \
  --approval-json X23b/data/gxtb_staging/x23b_k222_fd_approval.json \
  --reviewer "<name>" \
  --coordinate-abs-tolerance-hartree-per-bohr "$COORD_TOL_HA_PER_BOHR" \
  --stress-abs-tolerance-gpa "$STRESS_TOL_GPA"
```

After reviewing the complete preflight table, prepare the 23 native-Bloch k222
cell optimizations from the same frozen experimental-reference inputs and the
matching stamped preflight outputs, then run and collect them:

```bash
KROOT=X23b/runs/gxtb_native
python3 X23b/scripts/x23b_kpoint_cellopt.py prepare \
  --source-policy experimental_reference --preflight-root "$PREFLIGHT_ROOT" \
  --output-root "$KROOT/k222" --method GXTB
python3 X23b/scripts/x23b_kpoint_cellopt.py run \
  --output-root "$KROOT/k222" --method GXTB --cp2k "$CP2K" \
  --cp2k-source "$CP2K_SOURCE" --save-tblite "$SAVE_TBLITE" \
  --save-tblite-source "$SAVE_TBLITE_SOURCE" \
  --campaign-manifest "$CAMPAIGN_MANIFEST" \
  --jobs 6 --threads-per-job 1 --prune-transients
# Only if RUN reports MAX_ITER; this keeps the same model and native mixer:
python3 X23b/scripts/x23b_kpoint_cellopt.py continue-maxiter \
  --output-root "$KROOT/k222" --method GXTB --cp2k "$CP2K" \
  --cp2k-source "$CP2K_SOURCE" --save-tblite "$SAVE_TBLITE" \
  --save-tblite-source "$SAVE_TBLITE_SOURCE" \
  --campaign-manifest "$CAMPAIGN_MANIFEST" \
  --jobs 2 --threads-per-job 1 --prune-transients
python3 X23b/scripts/x23b_kpoint_cellopt.py collect \
  --output-root "$KROOT/k222" --method GXTB --molecule-run-root X23b \
  --csv X23b/data/gxtb_staging/x23b_k222_cellopt_results.csv
```

`experimental_reference` is GXTB-only, uses the separate
`k222_cellopt_keep_angles_from_experimental_reference` variant, and forbids
`--override`. Policy, frozen structure/input hashes, and preflight input/output
hashes are immutable parts of every k222 job stamp. Thus a Gamma result, a
policy-less old manifest, or a modified source artifact is rejected as stale.
The legacy matching-Gamma route remains available without changing GFN1/GFN2:
use an independent output root and explicitly select
`--source-policy gamma_cellopt_restart --gamma-root X23b`.

Run the 23 k333 and 23 k444 single points on each system's own final GXTB-k222
restart and build the strict 23-system mesh table:

```bash
for MESH in 3 4; do
  python3 X23b/scripts/x23b_final_kpoint_sp.py prepare \
    --cellopt-root "$KROOT/k222" --output-root "$KROOT/k${MESH}${MESH}${MESH}" \
    --mesh "$MESH" --method GXTB
  python3 X23b/scripts/x23b_final_kpoint_sp.py run \
    --output-root "$KROOT/k${MESH}${MESH}${MESH}" --mesh "$MESH" \
    --method GXTB --cp2k "$CP2K" --jobs 6 --threads-per-job 1 \
    --cp2k-source "$CP2K_SOURCE" --save-tblite "$SAVE_TBLITE" \
    --save-tblite-source "$SAVE_TBLITE_SOURCE" \
    --campaign-manifest "$CAMPAIGN_MANIFEST" \
    --prune-transients
  python3 X23b/scripts/x23b_final_kpoint_sp.py collect \
    --output-root "$KROOT/k${MESH}${MESH}${MESH}" --mesh "$MESH" --method GXTB \
    --molecule-run-root X23b \
    --csv "X23b/data/gxtb_staging/x23b_k${MESH}${MESH}${MESH}_results.csv"
done
python3 X23b/scripts/summarize_x23b_final_kpoints.py \
  --method GXTB \
  --k333-csv X23b/data/gxtb_staging/x23b_k333_results.csv \
  --k444-csv X23b/data/gxtb_staging/x23b_k444_results.csv \
  --rows-csv X23b/data/gxtb_staging/x23b_final_geometry_kpoint_rows.csv \
  --summary-csv X23b/data/gxtb_staging/x23b_final_geometry_kpoint_summary.csv
python3 X23b/scripts/x23b_pipeline.py analyse --method GXTB --skip-plots \
  --cellopt-csv X23b/data/gxtb_staging/x23b_k222_cellopt_results.csv \
  --final-kpoint-csv X23b/data/gxtb_staging/x23b_final_geometry_kpoint_rows.csv \
  --output-dir X23b/data/gxtb_staging
```

This V1 production protocol is 115 CP2K jobs: 23 gas optimizations, 23 frozen-
structure k222 derivative preflights, 23 k222 optimizations, and 23 single
points at each of k333 and k444. The focused finite-difference validation adds
36 separate pilot single points; it is not hidden in that production count.
All GXTB resumptions require a per-job
`job_provenance.json` whose
input/output hashes and complete campaign fingerprint match the frozen
manifest. A missing, stale, or mixed-build stamp produces the non-success
`STALE_OUTPUT`; use the GXTB-only `--force` only after reviewing the mismatch.
A concurrent lock produces non-success `BUSY`.
The k222 run, continuation, and polish commands select only the primary inputs
owned by `x23b_k222_cellopt_manifest.csv`, never generated continuation inputs.
`--force`, `--clean`, and `--prune-transients` are guarded so
they cannot delete the published GFN1/GFN2 trees.  Pruning occurs only after a
phase-specific success check; it retains the final continuation restart,
current WFN, final structures, inputs, outputs, and a JSON deletion record.

The GXTB provenance file is updated atomically by every production batch and
collector. It records a derived snapshot of the central campaign identity,
artifact/source hashes, workflow roots, and campaign-matched completion counts.
Collectors reject stamps from any other build. The fixed filename is
`data/build_provenance_gxtb.json`; `data/build_provenance.json` remains frozen.
Full-grid GXTB outputs made before the SPGLIB production path was available are
retained only as diagnostics. They have no matching production stamp and are
not eligible for collection, continuation, or publication tables.

### Fail-closed paper artifact

The staging analysis above is not a publication boundary. After the final
direct-ACP k-point build has been validated against the frozen predecessor,
freeze the GFN1/GFN2/GXTB comparison with:

```bash
python3 X23b/scripts/finalize_x23b_paper_summary.py
```

The command atomically creates
`X23b/data/x23b_gfn_gxtb_paper_summary.csv` and the complete hash/lineage
manifest `X23b/data/x23b_gfn_gxtb_paper_summary.json`. It requires exact common
23-system coverage, all stamped gas/preflight/k222/k333/k444 raw jobs, the
approved FD pilot, the exact preflight-to-CELL_OPT and final-restart lineage,
the production-ready build manifest, and two additional reviewed artifacts:

- `X23b/data/gxtb_staging/x23b_direct_acp_cross_build_approval.json`, binding
  the final direct-ACP build and frozen predecessor manifests to 23 paired raw
  input/output/stamp comparisons;
- `X23b/data/gxtb_staging/x23b_k333_k444_convergence_approval.json`, recording
  a named reviewer's explicit mean-absolute and maximum-absolute convergence
  thresholds and passing checks for GFN1, GFN2, and GXTB.

Missing, incomplete, rejected, or hash-inconsistent gates remove both paper
outputs rather than leaving a partial or stale publication bundle. The frozen
GFN1/GFN2 source tables and `data/build_provenance.json` are read and hashed,
never rewritten. Non-default gate or output paths can be supplied with
`--cross-build-approval`, `--kpoint-approval`, `--output-csv`, and
`--output-json`.

### Optional fixed-reference SI diagnostic

The relaxed-cell benchmark above remains the default. To reproduce the
additional GFN1/GFN2-style SI diagnostic on the frozen reference crystals,
explicitly opt in to 23 GXTB single points at each of Gamma, k111, k222, and
k333 (92 additional jobs). Its explicit k-point inputs also use SPGLIB
reduction; a full-grid comparison, if desired, is a separate diagnostic and
must not be mixed into these tables:

```bash
python3 X23b/scripts/x23b_pipeline.py prepare --method GXTB --production-only \
  --include-fixed-reference
python3 X23b/scripts/x23b_pipeline.py run --method GXTB --phase crystal_sp \
  --cp2k "$CP2K" --cp2k-source "$CP2K_SOURCE" \
  --save-tblite "$SAVE_TBLITE" --save-tblite-source "$SAVE_TBLITE_SOURCE" \
  --campaign-manifest "$CAMPAIGN_MANIFEST" \
  --jobs 4 --threads-per-job 1 --prune-transients
python3 X23b/scripts/x23b_pipeline.py analyse --method GXTB --skip-plots \
  --include-fixed-reference \
  --cellopt-csv X23b/data/gxtb_staging/x23b_k222_cellopt_results.csv \
  --final-kpoint-csv X23b/data/gxtb_staging/x23b_final_geometry_kpoint_rows.csv \
  --output-dir X23b/data/gxtb_staging
```

Supplying `--phase crystal_sp` is the opt-in boundary and defaults to all four
diagnostic meshes. Individual meshes can be selected by repeating
`--fixed-reference-mesh gamma|k111|k222|k333`. GXTB analysis with
`--include-fixed-reference` refuses partial 92-case coverage.

The optional Gamma CP2K-native versus save_tblite-CLI energy/gradient/virial
diagnostic is deliberately separate from production:

```bash
python3 scripts/run_x23b_reference_cli_checks.py \
  --benchmark-root . --method GXTB --cp2k "$CP2K" \
  --tblite "$SAVE_TBLITE" \
  --out X23b/runs/reference_cli_gxtb --jobs 6 --resume
```

The candidate build and current stamp-validated completion state are recorded in
`data/build_provenance_gxtb.json`. Curated GFN1/GFN2 tables are not touched by
the GXTB commands: analysis defaults to/stays in `data/gxtb_staging` until all
five 23/23 coverage checks have passed and the provenance counts are updated.
