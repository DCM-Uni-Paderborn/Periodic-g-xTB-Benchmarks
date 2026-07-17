# DMC-ICE13 Periodic GFN Benchmark

This directory contains CP2K/tblite single-point calculations for the
DMC-ICE13 ice polymorph benchmark. The calculations compare periodic
GFN1-xTB and GFN2-xTB relative energies against the diffusion Monte Carlo
reference values of Della Pia, Zen, Alfe, and Michaelides,
J. Chem. Phys. 157, 134701 (2022), DOI: 10.1063/5.0102645.

## Data included

- `poscars/`: POSCAR geometries for the 13 DMC-ICE13 polymorphs.
- `inputs/`: Gamma-only CP2K input files for GFN1-xTB and GFN2-xTB.
- `kpoint_inputs/`: explicit native Bloch 1x1x1, 2x2x2, 3x3x3, 4x4x4, and
  5x5x5 MacDonald k-point CP2K input files.
- `runs/`: generated Gamma-only CP2K working directories, ignored by Git.
- `runs_kpoints/`: generated k-point CP2K working directories, ignored by Git.
- `data/results.json`: raw CP2K total energies, per-water energies, relative
  energies with respect to ice Ih, and error statistics for the Gamma-only
  calculations.
- `data/kpoint_results.json`: raw and relative energies for the k-point
  dependent calculations.
- `data/dmc_ice13_relative_energies.csv`: 3x3x3 relative energies and GFN
  errors used as the primary manuscript values.
- `data/dmc_ice13_kpoint_stats.csv`: aggregate DMC-ICE13 error statistics as a
  function of k-point mesh.
- `data/dmc_ice13_kpoint_relative_energies.csv`: phase-resolved relative
  energies and errors as a function of k-point mesh.
- `data/dmc_ice13_gxtb_adaptive_frontier.{csv,json}`: the diagnostic adaptive
  frontier from explicit unshifted `1x1x1` (presented as Gamma after the
  separate route-equivalence check) through `11x11x11`. At each frontier,
  accepted phases retain their denser endpoint and unresolved phases use the
  current mesh; this diagnostic is not an additional stopping gate.
- `data/dmc_ice13_gxtb_phase_vii_kpoint_convergence.csv`: the complete ice-VII
  relative-energy and adjacent-mesh sequence used in the new Supporting
  Information.
- `data/previous_vs_full_pr350_mae.csv` and the companion Markdown file:
  explicit comparison with the earlier partial-PR350 manuscript stack.
- `data/dmc_ice13_relative_mae_comparison.csv`: comparison with the published
  DFT data from the DMC-ICE13 paper.
- `data/dmc_ice13_published_dft_absolute_energies.csv`: published DMC and DFT
  absolute lattice energies from the DMC-ICE13 paper, used to compute the
  relative-energy MAE ranking.
- `data/build_provenance.json`: source revisions, executable and shared-library
  hashes, patch hashes, build flags, and the completed-calculation count.
- `data/dmc_ice13_reference_cli_rows.csv` and
  `data/dmc_ice13_reference_cli_summary.csv`: direct CP2K-native versus tblite
  CLI energy, gradient, and virial checks for all 26 Gamma calculations.
- `figures/`: PDF, SVG, and PNG versions of the three DMC-ICE13 plots used in
  the revised manuscript and Supporting Information.
- `scripts/`: input generation, extraction, analysis, plotting, and run scripts.

The original PDF and Supporting Information are not redistributed here. The
geometries and DMC reference values are documented through the paper DOI above.

## CP2K setup used

The calculations were run from CP2K development trunk, not from a numbered
release. The executable reports `2026.1 (Development Version)` and is
interfaced to tblite:

- CP2K source revision: `faf9aae91266170dfee8a9f7171a5135bc5eb368`
- CP2K flags reported by the executable: `omp no_statm_access spglib libdftd4
  dftd4_v4_2 s_dftd3 mctc-lib tblite`
- tblite: `8a9d09474b93d25c044d6f46ce920750c7fe4cf7` (`tblite` 0.6.0),
  merging current `main` with the complete tblite PR 350 series through
  `8c5e562`; the earlier PR 343 changes are also included
- CP2K working-tree additions: overlap-covariant native-Bloch full symmetry,
  analytical force/stress checks, separate Bloch-wavefunction and multipolar
  SCC restart handling, and Broyden-mixer history fixes
- `TBLITE/ACCURACY`: `0.1`
- `EPS_SCF`: `1.0E-9`
- run-script defaults: `OMP_NUM_THREADS=1`,
  `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and
  `CP2K_PARALLEL_JOBS=20`, i.e. independent single-core CP2K jobs are launched
  concurrently.

The primary comparison uses the Gamma-centered 3x3x3 k-point mesh, matching
the non-hybrid DFT single-point setup in the DMC-ICE13 reference. The explicit
1x1x1 mesh verifies equivalence to the Gamma-only calculation, the 2x2x2 mesh
documents the approach to convergence, and the 4x4x4 and 5x5x5 checks confirm
that the 3x3x3 aggregate statistics are converged. All energies in the CSV
summaries are relative to ice Ih and reported in kJ mol-1 per water molecule.
The independent reference-CLI checks use analytical CP2K stress tensors; all
26 calculations complete, with maximum force-component differences of
`3.57e-8` (GFN1) and `1.51e-7` atomic units and maximum virial-component
differences of `1.26e-6` and `2.29e-6` atomic units, respectively.

Current aggregate MAEs:

| Mesh | GFN1-xTB | GFN2-xTB |
|---|---:|---:|
| Gamma | 6.694624 | 5.578897 |
| 2x2x2 | 7.956838 | 3.510100 |
| 3x3x3 | 8.005255 | 3.462919 |
| 4x4x4 | 8.006494 | 3.461424 |
| 5x5x5 | 8.006485 | 3.461353 |

## Additive g-xTB production workflow

The follow-on g-xTB benchmark is kept separate from the frozen GFN1-xTB and
GFN2-xTB paper data above.  Running the g-xTB workflow does not replace the
versioned GFN result JSON/CSV files, their figures, or their raw run
directories.

The production contract is:

- the same 13 fixed DMC-ICE13 reference geometries and six sampling members as
  for GFN1/GFN2, for 78 g-xTB single points in total;
- implicit Gamma without a `&KPOINTS` section, followed by an independent
  explicit `MACDONALD 1 1 1 0.0 0.0 0.0` calculation;
- for every larger mesh, `SYMMETRY T`, `FULL_GRID F`,
  `SYMMETRY_BACKEND SPGLIB`, and `SYMMETRY_REDUCTION_METHOD SPGLIB`;
- CP2K expands the irreducible density and overlap to the complete coupled
  mesh internally for the single save_tblite g-xTB evaluation, and folds the
  response back to the irreducible representation;
- `TBLITE/ACCURACY 0.1`, `EPS_SCF 1.0E-9`, and the native save_tblite Fock-DIIS
  path selected by `SCC_MIXER TBLITE` with an explicit 300-iteration limit.

The isolated production locations are:

- `gxtb_spglib_inputs/`: generated inputs;
- `runs_gxtb_spglib/`: raw outputs and per-job hash stamps (ignored by Git);
- `data/dmc_ice13_gxtb_spglib_*`: additive JSON/CSV analysis products;
- `data/build_provenance_gxtb_spglib.json`: executable, source, input, output,
  invocation, and protocol provenance;
- `figures/dmc_ice13_gxtb_spglib_*`: additive plots once all meshes are
  complete.

Earlier g-xTB files below `kpoint_inputs/`, `runs/GXTB/`, and
`runs_kpoints/*/GXTB/` used the pre-symmetry full-grid development route.
They remain available only as diagnostics.  The production runner never reads
them as benchmark results and records their hashes separately under
`legacy_diagnostics` in the new provenance file.

After building the intended clean CP2K and save_tblite revisions, the complete
matrix can be run with, for example:

```bash
python3 scripts/run_dmc13_kpoint_jobs.py \
  --root DMC-ICE13 \
  --cp2k /path/to/cp2k.psmp \
  --tblite /path/to/tblite \
  --tblite-static-library /path/to/lib/libtblite.a \
  --cp2k-source /path/to/cp2k-source \
  --tblite-source /path/to/save_tblite-source \
  --method GXTB \
  --jobs 4
```

For MPI production on Terok, each worker instead receives a literal ordered
PE list with exactly one unique logical CPU per MPI rank.  The runner injects
Open MPI's `--map-by pe-list=...:ordered --bind-to core --report-bindings`,
sets one OpenMP and BLAS thread per rank, removes inherited OMPI/PRTE binding
overrides, and verifies rank-numbered singleton masks through `/proc` before it
can finalize a schema-v2 execution record.  It also holds one host-local
`flock` per requested logical CPU, so a second independently started
production driver fails before launch instead of sharing that CPU. Earlier
misbound samples remain sticky even if a later sample looks correct.  After
locking and again immediately before each launch, a Linux `/proc` preflight
also rejects selected CPUs already present in any live non-zombie CP2K or MPI
rank mask, including jobs started by older non-locking launchers.  Sequential
same-rank/same-mask PID generations created during sanitizer teardown are
folded into one rank proof; concurrently live duplicate ranks, rank migration,
or any successor mask change remain fatal. Extra MPI
launcher arguments are not accepted because Open MPI aliases, appfiles, and
MCA parameter files could otherwise bypass the immutable mapping contract. For
example:

```bash
python3 scripts/run_dmc13_kpoint_jobs.py \
  --root DMC-ICE13 \
  --cp2k /path/to/cp2k.psmp \
  --tblite /path/to/tblite \
  --tblite-static-library /path/to/lib/libtblite.a \
  --cp2k-source /path/to/cp2k-source \
  --tblite-source /path/to/save_tblite-source \
  --method GXTB --jobs 2 --threads-per-job 1 \
  --mpi-ranks-per-job 4 --mpi-launcher /path/to/mpirun \
  --pe-list 96,97,98,99 --pe-list 100,101,102,103
```

If a physical core exposes multiple hardware threads and `--bind-to core`
therefore produces a non-singleton `/proc` mask, the run is rejected. No SMT
exception is inferred from the requested PE list.

Resume is bound to the exact execution mode, MPI-rank count, thread settings,
and execution-contract hash. A direct invocation therefore cannot silently
reuse an MPI-produced stamp, or vice versa. Aggregate timing provenance is
scaling-eligible only if every included job has return code zero and a fully
revalidated schema-v2 execution sidecar; mixed, legacy, and failed populations
are classified `timing_non_scaling`.

The older shell launchers and any schema-v1 taskset/`--bind-to none` records
are retained unchanged as historical raw provenance.  Their energies, forces,
and stresses remain usable after the normal hash and scientific validation,
but all multi-rank timings from that shared-mask policy are classified
`legacy_timing_non_scaling` and must not be used in speedup or scaling plots.

This production runner is intentionally g-xTB-only: `--method GXTB` must be
given exactly once, so it cannot regenerate or overwrite the frozen GFN1/GFN2
inputs and results.  It additionally requires the exact additive analysis
prefix `gxtb_spglib` and keeps both generated-input and run roots below this
DMC-ICE13 directory.  Unprefixed analysis remains restricted to GFN1/GFN2.

A result is admitted to analysis only if the output completed with a converged
SCF, the input satisfies the explicit SPGLIB/implicit-Gamma contract, and the
input, output, CP2K launcher, the `libcp2k` selected through `otool`/RPATH on
macOS or `ldd` on Linux, and
the statically linked `libtblite.a` hashes match its stamp.  The default
`campaigns/gxtb-pbc-v1-20260714/build_manifest.json` additionally freezes all
four build artifacts; the CP2K-reported revision must resolve exactly to the
CP2K source HEAD, and the save_tblite source HEAD must equal the manifest
revision.  `--force` archives prior production files instead of deleting a run
directory.  Unstamped output is never retroactively adopted as production,
because its producing executable cannot be established from the output alone.

### Optional dense convergence extensions

The frozen production core remains exactly the six meshes from implicit Gamma
through `k555` (78 jobs).  `k666` through `k131313` are opt-in
convergence extensions; they are never added by the runner default and do not
change the meaning of a completed core campaign or any existing stamp.  Their
MacDonald definitions are `6 6 6 5/12 5/12 5/12`, `7 7 7 0 0 0`,
`8 8 8 7/16 7/16 7/16`, `9 9 9 0 0 0`, `10 10 10 9/20 9/20 9/20`,
`11 11 11 0 0 0`, `12 12 12 11/24 11/24 11/24`, and
`13 13 13 0 0 0`, respectively, with the same SPGLIB-reduced/full-coupled-mesh
contract as the core.

An eight-phase `k666` triage pilot can be selected by repeating `--phase`:

```bash
python3 scripts/run_dmc13_kpoint_jobs.py \
  --root DMC-ICE13 \
  --cp2k /path/to/cp2k.psmp \
  --cp2k-library /path/to/libcp2k.dylib \
  --tblite /path/to/tblite \
  --tblite-static-library /path/to/lib/libtblite.a \
  --cp2k-source /path/to/cp2k-source \
  --tblite-source /path/to/save_tblite-source \
  --method GXTB --mesh k666 \
  --phase Ih --phase VII --phase XIV --phase XI --phase VIII \
  --phase VI --phase XV --phase XVII \
  --jobs 1
```

Every invocation holds a nonblocking campaign lock, verifies both the generated
input and the frozen copy that CP2K actually executed, and writes a deterministic
content-addressed validation snapshot from currently valid stamps.  The named
validation index is only an atomic current-copy; analysis receives the immutable
snapshot path.  The run then emits the validation index, the DMC-independent
fixed-mesh convergence report, and the DMC-referenced phase-wise result:

- `data/dmc_ice13_gxtb_spglib_validation_index.json`;
- `data/dmc_ice13_gxtb_spglib_kpoint_convergence.json` and the companion
  phase-level CSV;
- `data/dmc_ice13_gxtb_spglib_phasewise_kpoint_convergence.json` and its CSV.

Once all twelve non-reference phases are marked `phasewise_kpoint_converged`,
freeze the three-method paper table and its complete raw-energy/provenance
lineage with:

```bash
python3 DMC-ICE13/scripts/finalize_dmc13_phasewise_summary.py \
  --root DMC-ICE13
```

This refuses an incomplete GFN1/GFN2/g-XTB result and removes stale summaries.
On success it writes exactly one row per method to
`data/dmc_ice13_gfn_gxtb_phasewise_summary.csv` and the phase-resolved energies,
same-mesh Ih references, mesh choices, source hashes, build identities, and DMC
reference sensitivity to the companion JSON file.  The same three rows also
carry a strictly identical `k333` fixed-mesh comparison for GFN1-xTB,
GFN2-xTB, and g-xTB.  This fixed-mesh block is a cross-method diagnostic, not
the primary converged result; the g-xTB `k333` value is explicitly labelled
`numerically_unconverged_same_mesh_comparator` and must not be substituted for
the phase-wise k-point-converged value.  The raw g-xTB matrix was executed
before CP2K PR #5582 and retains that build identity in every raw record.  It
is not relabelled as post-#5582 data.  The energy-only sentinel qualification
below now authorizes its use in the DMC-ICE13 paper table while preserving the
pre-#5582 origin explicitly.

### Same-source builds on another architecture

Schema-v1 validation snapshots remain immutable and readable.  New snapshots
use schema v2: every `(mesh, phase)` record names its own content-addressed
execution-build identity, so records made by the frozen macOS/ARM64 build and a
qualified Linux/x86_64 build can coexist without weakening any per-file or
per-stamp check.  All build identities in one index must use the exact frozen
CP2K and save_tblite source revisions.

A build whose executable/library hashes differ from `build_manifest.json` is
accepted only with an explicit `--execution-build-manifest`.  This additive
JSON file must contain schema version 1, the campaign and protocol IDs, the
SHA256 of the unchanged campaign manifest, the deterministic `build_id`, a
`build_identity` block with the five stamp-bound fields plus
`tblite_cli_sha256`.  `qualification` must use evidence schema 3, declare
finite total- and relative-energy tolerances no looser than `1e-10` Eh and
`0.001` kJ mol-1 per H2O, report matching observed maxima within those
tolerances, and contain at least one counted same-mesh dense relative-energy
sentinel.  Each sentinel binds the non-Ih phase and Ih inputs, four completed
outputs (remote phase, remote Ih, reference phase, and reference Ih), and all
four corresponding run stamps by safe relative paths and 64-hex SHA256 hashes,
and binds the remote and reference sides to the execution-manifest and
frozen-campaign `build_id`,
respectively.  The inputs must have the exact
`ice_<phase>_GXTB_<mesh>.inp` names, projects, and complete frozen g-xTB/SPGLIB
contract.  Their declared positive water counts are checked against the
explicit O atoms in the hashed `&COORD` blocks.  Phase and Ih inputs must be
distinct; all four output paths and all four stamp paths must be distinct.
Phase/Ih hashes within one build and remote/reference hashes for one system
must differ.  The readers hash and parse the same bytes.  Each output must
report exactly one matching input, project, CP2K revision, tblite revision,
and energy header.  Alternate outputs require the full 40-hex save_tblite
revision; frozen legacy outputs may report `unknown` because their exact bytes
and stamps are bound by the pinned reference record.  Fresh remote stamps must
be schema v2, bind the five-field build identity and source/frozen-input/output
hashes, and state `adopted_existing_output: false`.  Frozen reference outputs
and stamps must exactly match the SHA256-pinned base index.
The sentinel also records the fixed conversion
`hartree_to_kjmol = 2625.499638`.  The readers parse all four hashed outputs,
recompute and verify the separate phase and Ih cross-build total-energy
deltas, recompute both remote and reference Ih-referenced relative energies
from the raw totals and water counts, and then verify the declared relative
delta.  `observed_max_abs_total_energy_delta_hartree` must equal the maximum
over both the phase and Ih total-energy deltas of every sentinel.  The manifest
and every evidence path must be relative to this benchmark root; absolute,
`..`, and escaping symlink paths are rejected.  The schema-v2 index itself
content-addresses the qualification manifest.  Both source worktrees must be
clean, and a qualified alternate library must embed the full save_tblite
revision rather than `unknown`.

The qualification declares the expected remote `PROGRAM STARTED ON`,
`Program compiled on`, and `Program compiled for` values.  Both remote outputs
must match them; this ARM64/macOS-to-x86_64/Terok campaign additionally expects
the reference and remote host/platform fingerprints to differ.  This is an
accidental-copy defense, not cryptographic execution attestation.  The record
establishes hash-bound internal consistency and still relies on the manually
trusted SSH transfer.  Cryptographic proof would require a signature or
trusted scheduler receipt whose trust root is outside this artifact tree.

The compact qualification block is:

```json
{
  "status": "passed",
  "evidence_schema_version": 3,
  "remote_execution_environment": {
    "program_started_on": "terok",
    "program_compiled_on": "terok",
    "program_compiled_for": "x86_64"
  },
  "total_energy_tolerance_hartree": 1e-10,
  "relative_energy_tolerance_kjmol_per_h2o": 0.001,
  "observed_max_abs_total_energy_delta_hartree": 0.0,
  "observed_max_abs_relative_energy_delta_kjmol_per_h2o": 0.0,
  "same_mesh_dense_relative_sentinel_count": 1,
  "same_mesh_dense_relative_sentinels": [{
    "kind": "same_mesh_dense_relative_energy",
    "mesh": "k666",
    "phase": "VII",
    "remote_build_id": "sha256:<64 lowercase hex>",
    "reference_build_id": "sha256:<64 lowercase hex>",
    "phase_input": "data/execution_builds/<build>/requalification/k666/VII/ice_VII_GXTB_k666.inp",
    "phase_input_sha256": "<64 lowercase hex>",
    "ih_input": "data/execution_builds/<build>/requalification/k666/Ih/ice_Ih_GXTB_k666.inp",
    "ih_input_sha256": "<64 lowercase hex>",
    "remote_phase_output": "data/execution_builds/<build>/requalification/k666/VII/ice_VII_GXTB_k666.out",
    "remote_phase_output_sha256": "<64 lowercase hex>",
    "remote_ih_output": "data/execution_builds/<build>/requalification/k666/Ih/ice_Ih_GXTB_k666.out",
    "remote_ih_output_sha256": "<64 lowercase hex>",
    "reference_phase_output": "runs_gxtb_spglib/k666/VII/ice_VII_GXTB_k666.out",
    "reference_phase_output_sha256": "<64 lowercase hex>",
    "reference_ih_output": "runs_gxtb_spglib/k666/Ih/ice_Ih_GXTB_k666.out",
    "reference_ih_output_sha256": "<64 lowercase hex>",
    "remote_phase_stamp": "data/execution_builds/<build>/requalification/k666/VII/ice_VII_GXTB_k666.run.json",
    "remote_phase_stamp_sha256": "<64 lowercase hex>",
    "remote_ih_stamp": "data/execution_builds/<build>/requalification/k666/Ih/ice_Ih_GXTB_k666.run.json",
    "remote_ih_stamp_sha256": "<64 lowercase hex>",
    "reference_phase_stamp": "runs_gxtb_spglib/k666/VII/ice_VII_GXTB_k666.run.json",
    "reference_phase_stamp_sha256": "<64 lowercase hex>",
    "reference_ih_stamp": "runs_gxtb_spglib/k666/Ih/ice_Ih_GXTB_k666.run.json",
    "reference_ih_stamp_sha256": "<64 lowercase hex>",
    "phase_water_count": 12,
    "ih_water_count": 12,
    "hartree_to_kjmol": 2625.499638,
    "phase_total_energy_delta_hartree": 0.0,
    "ih_total_energy_delta_hartree": 0.0,
    "remote_relative_energy_kjmol_per_h2o": 0.0,
    "reference_relative_energy_kjmol_per_h2o": 0.0,
    "relative_energy_delta_kjmol_per_h2o": 0.0
  }]
}
```

Pass the immutable local snapshot with `--base-validation-index` and its exact
`--base-validation-index-sha256` on the remote run.  It is fully reverified,
its valid logical jobs are not rerun or archived, and only missing jobs are
appended.  A different record for an already indexed
`(mesh, phase)` is a hard collision, including under `--force`.  The resulting
schema-v2 snapshot records the base snapshot as its parent.  Transfer only the
new generated source inputs, frozen executed inputs, outputs, stamps, the
content-addressed qualification manifest, and all evidence artifacts into a
local staging directory.  Verify their hashes before placing them at their
relative paths, and never use an unqualified `rsync --delete` or overwrite an
existing production record.  The base snapshot is checked before input
generation and checked again immediately before the merged snapshot is
written; only missing selected phases are prepared.  Analysis similarly takes
the immutable snapshot with `--validation-index-sha256`; a named current-copy
or whitespace-modified copy cannot silently replace it.

For each pair of meshes the report evaluates changes of the Ih-referenced
relative energies, not changes of the DMC error.  The frozen stopping rule is a
maximum phase change of at most 0.10 and an RMS change of at most 0.05
kJ mol-1 per H2O over two consecutive, fully covered mesh refinements.  A
partial pilot is reported as `coverage: pilot`.  It rejects a candidate only if
an observed maximum exceeds 0.10 or the full-set RMS lower bound
`sqrt(sum(observed_delta^2)/12)` exceeds 0.05.  Otherwise it is explicitly
inconclusive, even when the observed-subset RMS exceeds 0.05.  A pilot is always
`eligible_for_stopping: false`; only 13/13 hash-valid meshes can establish
formal convergence, and only the trailing contiguous sequence of full adjacent
mesh pairs determines the stopping state.

The separate phase-wise result uses one direct rule: a phase is k-point
converged at `N x N x N` when its same-mesh-Ih-referenced relative energy
changes by at most 0.05 kJ mol-1 per H2O from `(N-1) x (N-1) x (N-1)`.  The
smallest such mesh not contradicted by already available denser evidence is
reported for every phase.  Once all 12 non-reference phases are selected, the
report emits the **phase-wise k-point-converged MAE** together with each
relative energy, DMC reference, error, last delta, mesh label, `N`, and total
k-point count.  RMS, mean absolute value, and maximum of the 12 last deltas are
diagnostics only; they do not add a second acceptance condition.

The JCP draft convergence figures and their exact plotted tables are regenerated
from the hash-pinned raw and phase-wise JSON files with

```bash
python3 DMC-ICE13/scripts/make_dmc_adaptive_convergence_plots.py
```

The main-text figure shows the adaptive-frontier MAE through `11x11x11`; the
SI figure resolves ice VII and demonstrates that `9x9x9 -> 10x10x10` does not
pass, whereas `10x10x10 -> 11x11x11` changes the relative energy by only
0.0311 kJ mol-1 per H2O and therefore satisfies the 0.05 criterion.

The independent Gamma CLI comparison is likewise additive.  Run
`scripts/dmc_gxtb_gamma_cli_check.py` with both `--cp2k` and `--tblite` only
after the Gamma production jobs exist.  By default it reads CP2K references
from `runs_gxtb_spglib/`, requires their V1 protocol stamps to match the
current launcher, loaded `libcp2k`, static `libtblite.a`, source revisions, and
central campaign manifest, and writes stamped CLI jobs to
`runs_cli_gxtb_spglib/`.  Its CSV/JSON products use the
`dmc_ice13_gxtb_spglib_*` prefix.  Earlier `runs_cli/GXTB/` files and the
unprefixed CLI tables are not adopted or modified.

### Post-upstream cross-build requalification

After a CP2K integration commit is merged upstream, the frozen benchmark
energies are reusable only after a new common production build reproduces a
hash-pinned cross-build matrix.  The dedicated
`scripts/requalify_dmc13_cross_build.py` runner consumes both the immutable
reference validation snapshot and the common candidate-build manifest by exact
SHA256.  It verifies clean source trees, requires the candidate CP2K revision
to contain upstream `c92cc08b45378b85150447011b5a4bb552f5b797`, executes the
exact frozen input bytes, and compares total energies within `1e-10` Eh and
same-mesh-Ih relative energies within `0.001` kJ mol-1 per H2O.

There are two deliberately separate scopes.  A `--scope sentinel` run takes
explicit same-mesh Ih/phase selections and can establish only
`sentinel_passed`; it always writes `old_results_reusable: false`.  This is the
fast gate for implicit Gamma, representative even/odd meshes, or the final
adjacent Ih/VII pair.

The publication gate uses `--scope full-publication-matrix` and does not accept
free-form `--selection` arguments.  It derives its exact matrix from a final
summary supplied with `--final-summary` and `--final-summary-sha256`.  The
matrix contains the selected and immediately previous mesh for all twelve
non-reference phases, each corresponding same-mesh Ih energy, and all thirteen
raw `k333` energies used by the fixed-mesh comparator.  Only a hash-bound
`full_publication_matrix_passed` report from this automatically derived matrix
sets `old_results_reusable: true` and `paper_freeze_authorized: true` in this
generic runner.  The report path is removed before every attempt and remains
absent after any failure, preventing a stale passing report from surviving a
failed rerun.

#### Accepted DMC-ICE13 energy-only sentinel qualification

For the frozen DMC-ICE13 energy benchmark, an explicitly reviewed narrower
policy was accepted after the generic runner was written.  The existing 62/62
successful pre-#5582 production matrix is qualified by rerunning the exact
frozen ice-VII `k666` input with the production-ready post-#5582 build.  The
old and post-#5582 energies are respectively
`-917.472187146001147` and `-917.472187146001261` Ha, an absolute difference
of `1.14e-13` Ha; both calculations required 12 SCC iterations.  This is well
inside the accepted `1e-10` Ha energy tolerance.

The immutable candidate manifest is archived at
`campaigns/gxtb-pbc-v1-post5582-20260714/build_manifest.json` with SHA256
`b0feea6a411f02dedb1eb57190092e35d38f4c5705a985893d9f97070ddb1d51`.
`data/dmc_ice13_gxtb_post5582_energy_sentinel_qualification.json` binds that
manifest, the exact sentinel input, both output hashes, the 113-record frozen
validation index, the phase-wise source tables, and the original build
provenance.  The finalizer revalidates every binding before it can emit
`qualified_by_post5582_energy_sentinel`.  This exception qualifies only
DMC-ICE13 total and relative energies; it does not qualify forces, stress,
X23b, or LC10 data, and it never changes the pre-#5582 source revisions stored
with the raw results.
