LC10 native-Bloch CP2K/tblite paper benchmark
===============================================

This directory retains the 12 cubic solids studied by Goldzak, Wang, Ye, and
Berkelbach, J. Chem. Phys. 157, 174112 (2022).  The paper benchmark is the
fixed identical ten-system subset `C, Si, SiC, BN, BP, AlN, AlP, MgS, LiF,
LiCl` for GFN1-xTB, GFN2-xTB, and g-xTB.  LiH and MgO remain diagnostic data
only and never enter a reported statistic for any of the three methods.

Current production run
----------------------

The 2026-07-12 rerun used:

- CP2K trunk revision `faf9aae91266170dfee8a9f7171a5135bc5eb368` with the
  local CP2K/tblite interface patch recorded by hash in
  `data/build_provenance.json`;
- tblite revision `8a9d09474b93d25c044d6f46ce920750c7fe4cf7`, which combines current
  `main` with PR #350 and includes the previously merged PR #343;
- conventional cubic eight-atom cells and CP2K native Bloch sampling through
  `&KPOINTS` with `SCHEME MACDONALD`, `SYMMETRY T`, and `FULL_GRID F`, using
  full SPGLIB symmetry reduction; no Born-von-Karman supercells;
- a `k444` cubic equation of state, `k333/k444/k555` final single points, and
  `k555` as the reported cohesive-energy mesh;
- CP2K energies extrapolated to electronic temperature `T -> 0`;
- matching tblite CLI isolated-atom references with explicit atomic spins and
  `ACCURACY 0.05`.

Only completed SCF points enter an EOS fit. A quadratic fit is rejected when
its local RMSE exceeds 0.02 hartree or when its fitted minimum lies more than
0.02 hartree above the sampled local minimum.  All ten paper systems must pass
for every method.  Existing LiH/MgO branch and multistart material is retained
for diagnostics but is not a production gate and needs only a short methods
note in the manuscript.

Current versus previous results
-------------------------------

The values in this section are the frozen GFN1/GFN2 production results. Adding
GXTB is method-selective and does not replace these rows or their raw outputs.

| method | fixed coverage | lattice MAE (A) | cohesive-energy MAE (eV/atom) |
|---|---:|---:|---:|
| GFN1-xTB | 10/10 | 0.145118 | 1.543851 |
| GFN2-xTB | 10/10 | 0.062410 | 1.299325 |

These values are recomputed from the existing raw result table on precisely
the same ten systems. Historical 12-system and old/new diagnostic tables are
preserved under `data/`, but are not direct paper comparisons.

Reproduction
------------

Run from the repository root:

```bash
python3 Goldzak12/scripts/run_goldzak12_eos_benchmark.py \
  --cp2k /path/to/cp2k.ssmp \
  --tblite /path/to/tblite \
  --cp2k-source /path/to/cp2k \
  --tblite-source /path/to/tblite \
  --jobs 10 --threads 1 --eos-mesh k444 \
  --energy-mesh k333 --energy-mesh k444 --energy-mesh k555 \
  --result-mesh k555

python3 Goldzak12/scripts/compare_goldzak12_results.py --mesh k555
python3 Goldzak12/scripts/plot_literature_comparison.py
python3 Goldzak12/scripts/validate_goldzak12_results.py \
  --eos-mesh k444 --energy-mesh k333 --energy-mesh k444 \
  --energy-mesh k555 --result-mesh k555

python3 Goldzak12/scripts/continue_goldzak12_eos.py \
  --solid LiH --method GFN2 --mesh k444 --start-scale 0.94 \
  --scale 0.93 --scale 0.92 --variant lih_scc_continuation \
  --mixer tblite --memory 250 --damping 0.4 --promote
```

On Terok, use the shared fail-closed MPI/affinity path.  The manifest hash must
come from the completed build qualification; do not compute and accept a new
value inside the launch command.  Repeat `--solid` and `--cpu-set` as needed;
the number of CPU sets must equal `--jobs`, the sets must be disjoint, and each
must contain at least `MPI ranks times OpenMP threads` CPUs.  This compact
two-worker example demonstrates the exact CLI spelling:

```bash
MANIFEST=campaigns/gxtb-pbc-v1-post5582-20260714/build_manifest.json
MANIFEST_SHA256=REPLACE_WITH_EXTERNALLY_REVIEWED_SHA256
MPI=/path/to/qualified/environment/bin/mpirun

python3 Goldzak12/scripts/run_goldzak12_eos_benchmark.py \
  --method GXTB \
  --campaign-manifest "$MANIFEST" \
  --campaign-manifest-sha256 "$MANIFEST_SHA256" \
  --cp2k-source /path/to/clean/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/clean/save_tblite-pbc \
  --jobs 2 --mpi-ranks-per-job 8 --threads 1 \
  --mpi-launcher "$MPI" \
  --mpi-launcher-arg=--bind-to --mpi-launcher-arg=none \
  --taskset /usr/bin/taskset --cpu-set 96-103 --cpu-set 104-111 \
  --eos-mesh k444 --energy-mesh k333 --energy-mesh k444 \
  --energy-mesh k555 --result-mesh k555 --stop-after-eos
```

GXTB production extension
-------------------------

GXTB is run selectively with the save_tblite-enabled CP2K executable and the
matching save_tblite CLI. The generated CP2K inputs pin `METHOD GXTB`,
`SCC_MIXER TBLITE` (the native g-XTB FDIIS potential/Fock mixer), and CP2K
`DIRECT_P_MIXING`; no CP2K-Fock, CP2K-density, or modified-Broyden retry is
used as an alternative GXTB production mixer. The 11 LC10 elements use
`save_tblite run --method gxtb --spin 2S`, where `2S = multiplicity - 1` is
recorded alongside every energy.

All three methods now use the identical native-Bloch k-point contract:
`SCHEME MACDONALD`, `SYMMETRY T`, `FULL_GRID F`, and the SPGLIB backend and
reduction method. For GXTB, CP2K expands the irreducible density/overlap data
to the complete mesh before the coupled save_tblite evaluation and folds the
response back afterwards. Older GXTB inputs or outputs with `SYMMETRY F` and
`FULL_GRID T` are diagnostics only and are never accepted as LC10 production
data.

LC10 obtains its EOS and cohesive energies entirely from `RUN_TYPE ENERGY`.
GXTB energy and isolated-atom inputs therefore make no analytical-stress
request; this keeps the V1 energy benchmark independent of the later
force/stress implementation. The frozen GFN1/GFN2 inputs retain their existing
`STRESS_TENSOR ANALYTICAL` setting.

The diagnostic-only wavefunction-continuation results are summarized in
`data/gxtb_wfn_hysteresis.csv` and documented in
`data/gxtb_wfn_hysteresis.md`. They identify multiple self-consistent LiH
roots with strong path hysteresis and confirm that MgO at scale 0.85 enters a
symmetry-breaking/collapsed branch. These diagnostics are not included in the
production lattice-constant or cohesive-energy MAE.

```bash
python3 Goldzak12/scripts/run_goldzak12_eos_benchmark.py \
  --method GXTB \
  --campaign-manifest campaigns/gxtb-pbc-v1-post5582-20260714/build_manifest.json \
  --cp2k-source /path/to/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/save_tblite-pbc \
  --jobs 3 --threads 1 --eos-mesh k444 \
  --energy-mesh k333 --energy-mesh k444 --energy-mesh k555 \
  --result-mesh k555 --stop-after-eos

# After reviewing eos_fits.csv, gxtb_eos_branch_diagnostics.csv, and any
# classification/adaptive follow-up, approve this exact fit fingerprint:
python3 Goldzak12/scripts/run_goldzak12_eos_benchmark.py \
  --method GXTB \
  --campaign-manifest campaigns/gxtb-pbc-v1-post5582-20260714/build_manifest.json \
  --cp2k-source /path/to/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/save_tblite-pbc \
  --jobs 3 --threads 1 --eos-mesh k444 \
  --energy-mesh k333 --energy-mesh k444 --energy-mesh k555 \
  --result-mesh k555 --approve-fits --prune-transients

python3 Goldzak12/scripts/validate_goldzak12_results.py \
  --method GXTB --eos-mesh k444 \
  --energy-mesh k333 --energy-mesh k444 --energy-mesh k555 \
  --result-mesh k555
```

After all three methods have been collected, validate the combined archive and
freeze the publication table and its complete raw-output lineage:

```bash
python3 Goldzak12/scripts/validate_goldzak12_results.py \
  --method GFN1 --method GFN2 --method GXTB --eos-mesh k444 \
  --energy-mesh k333 --energy-mesh k444 --energy-mesh k555 \
  --result-mesh k555

python3 Goldzak12/scripts/finalize_goldzak12_paper_summary.py
```

The finalizer atomically writes `data/lc10_gfn_gxtb_paper_summary.csv`,
`data/lc10_gfn_gxtb_paper_summary.json`, and
`data/lc10_gfn_gxtb_paper_summary.tex`. The CSV has exactly one ten-system row
per method with ME, MAE, RMSE, and MaxAE for lattice constants and cohesive
energies. The JSON adds raw EOS/final-energy lineage, hashes, direct g-xTB/GFN
comparisons, and build provenance; the TeX file exports paper macros. Missing
or substituted systems, reduced coverage, stale stamps, or tampered raw data
remove all three outputs and make finalization fail.

`k555` is the runner and validator default result mesh, matching the frozen
GFN1/GFN2 report. The runner fixes `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`,
and `VECLIB_MAXIMUM_THREADS` to 1 and sets `OMP_WAIT_POLICY=PASSIVE`; CP2K's
outer `OMP_NUM_THREADS` remains controlled by `--threads`.

`scripts/benchmark_execution.py` leaves the scientific job-stamp schema and
matcher unchanged.  It writes an additive atomic `*.execution.json` record
that binds the exact taskset mask, launcher hash and command, `--bind-to none`,
observed CP2K child/rank PIDs and kernel CPU masks, input hash, output hash, and
scientific-stamp hash.  Kernel-normalized equivalent CPU-set spellings are
compared as sets.  A completed scientific output with missing or invalid
execution evidence is never deleted or rerun implicitly; the driver stops and
requires explicit review before `--force` can replace it.

GXTB final single points are never launched implicitly. `--stop-after-eos`
runs the atom/EOS stage and stops after writing the fits; `--fit-only`
recollects already stamped EOS outputs without launching an executable; and
only `--approve-fits` records the exact current fit-table fingerprint and
permits the `k333/k444/k555` stage. For the 38-GB development Mac,
`--jobs 3 --threads 1` is the conservative unmeasured default because the final queue
contains `k555` jobs; increase concurrency only after measuring their resident
set size. The 11 isolated-atom checks are small and can separately use more
workers.

The versioned `campaigns/gxtb-pbc-v1-post5582-20260714/build_manifest.json` is the
single source of truth for the CP2K launcher, its actually loaded `libcp2k`,
the save_tblite CLI, `libtblite.a`, their hashes, and the frozen source
revisions. The runner refuses every production invocation until this manifest
has `campaign_state: production_ready`. Its path above is also the runner
default. Explicit executable or library overrides are accepted only when they
resolve to the exact manifest artifacts. The runner additionally verifies the
dynamic loader resolution,
the CP2K embedded revision against both the manifest and clean source HEAD,
the clean save_tblite source HEAD, and then embeds the complete campaign
identity in every GXTB job stamp and collector. This path-independent identity
also freezes both CMake-cache hashes and the fetched-dependency lock; the
manifest path and complete file hash remain separate provenance records.
MPI/affinity launches additionally require the externally reviewed manifest
file hash through `--campaign-manifest-sha256`, and every GXTB EOS path requires
the CP2K source to descend from merged PR #5582 commit
`c92cc08b45378b85150447011b5a4bb552f5b797`.

An unsuccessful save_tblite atom job or CP2K job makes the command return
nonzero only after the concurrent batch has finished, so successful jobs are
kept. GXTB resume skips are accepted only when the prior output has a per-job
stamp whose input and complete manifest-derived campaign identity still matches
(launcher, loaded `libcp2k`, CLI, `libtblite.a`, both source revisions, both
CMake-cache hashes, and the dependency lock).
An atom JSON by itself is deliberately insufficient. GFN1/GFN2 keep their
existing retry protocol, but exhausted retries are also fatal instead of
silently producing partial tables.

Before accepting any GXTB atomic cohesive-energy reference, run the independent
CP2K-versus-save_tblite check (the 11 elements used by LC10, no solid jobs):

```bash
python3 Goldzak12/scripts/run_goldzak12_benchmark.py atom-check \
  --method GXTB \
  --campaign-manifest campaigns/gxtb-pbc-v1-post5582-20260714/build_manifest.json \
  --cp2k-source /path/to/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/save_tblite-pbc \
  --jobs 10 --threads 1 --tolerance-hartree 1e-6
```

This writes `data/atom_reference_cp2k_vs_save_tblite_gxtb.csv` and fails if an
atom is missing or exceeds the selected tolerance.

The default GXTB production plan is 11 save_tblite atom jobs, 110 CP2K k444
EOS jobs (ten solids times 11 standard scales), and 30 final single points.
Thus exact LC10 production has 151 jobs before the independent CP2K atom gate.
Every `--adaptive-scale SOLID=SCALE` adds one CP2K EOS job.

```bash
python3 Goldzak12/scripts/run_goldzak12_eos_benchmark.py \
  --method GXTB --adaptive-scale MgO=0.92000 --adaptive-scale MgO=0.90000 \
  --campaign-manifest campaigns/gxtb-pbc-v1-post5582-20260714/build_manifest.json \
  --cp2k-source /path/to/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/save_tblite-pbc --stop-after-eos
```

The separate continuation script may still be used for diagnostics, but GXTB
`--promote` is rejected: an independently copied continuation output cannot
satisfy the canonical input plus complete campaign-stamp contract.

LiH/MgO multi-start branch qualification
-----------------------------------------

LiH and MgO are outside the fixed LC10 paper benchmark. Their versioned
multi-start map remains an optional diagnostic and is not required before
LC10 production, validation, or finalization. It never contributes a reported
GFN1/GFN2/g-xTB statistic.

The runner rejects a CP2K source that is not descended from upstream commit
`c92cc08b45378b85150447011b5a4bb552f5b797` (merged PR #5582).  In particular,
the earlier `18d37c` build remains frozen diagnostic provenance and must not be
used for this map.  Pre-production execution also needs an explicit state
argument; the production default remains `production_ready`:

```bash
python3 Goldzak12/scripts/run_gxtb_multistart_branches.py \
  --campaign-manifest campaigns/gxtb-pbc-v1-post5582-20260714/build_manifest.json \
  --campaign-manifest-sha256 REPLACE_WITH_EXTERNALLY_REVIEWED_SHA256 \
  --cp2k-source /path/to/clean/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/clean/save_tblite-pbc \
  --cold-workers 2 --mpi-ranks-per-job 8 --threads 1 \
  --mpi-launcher /path/to/qualified/environment/bin/mpirun \
  --mpi-launcher-arg=--bind-to --mpi-launcher-arg=none \
  --taskset /usr/bin/taskset --cpu-set 112-119 --cpu-set 120-127

python3 Goldzak12/scripts/classify_gxtb_multistart_branches.py \
  --campaign-root Goldzak12/runs/gxtb_multistart_branches/FINGERPRINT
```

Classification re-hashes every input, output, WFN restart, job stamp, parent
restart, and parent candidate manifest.  It rejects missing scales, failed SCC
or symmetry gates, negative Mulliken populations, inverted LiH/MgO polarity,
inequivalent atoms, charge/electron-count failures, and discontinuous adjacent
charge/Fermi descriptors.  Among full physical paths it selects the lowest
relative-energy continuous branch and then applies the existing single-well
and quadratic-fit gates.  The selection remains diagnostic: it never copies an
output into `runs/eos`, never approves an LC10 fit, and never starts `k555`.

The independent CP2K-versus-CLI atom check adds 11 CP2K atom jobs. Counting
that mandatory acceptance gate, exact LC10 publication therefore contains 162
calculations before adaptive additions.

Final-single-point inputs are runnable only after a valid quadratic EOS fit.
Each generated input has an adjacent `*.inp.eos.json` lineage record containing
the fitted lattice constant, EOS mesh, energy mesh, input hash, and SPGLIB
contract. A pre-generated input at the experimental lattice constant, an input
without this valid lineage, or a lineage attached to an invalid fit is treated
as stale. Once the fit is accepted, the runner regenerates the input at the
actual EOS minimum before CP2K can run it.

The runner never forces a fit. A missing or discontinuous EOS minimum for any
of the ten paper systems stops before final single points and is recorded in
`data/gxtb_adaptive_followup.csv`/`.md`; the paper artifact still requires
10/10. LiH/MgO diagnostic behavior has no effect on this gate.

The exact per-system grid is persisted in
`data/gxtb_eos_scale_manifest.json`, restored on resume, and used by the
validator. A large local energy discontinuity is reported as a numerical SCC
branch *candidate*, never immediately as a physical failure. Each candidate
must be reviewed in `data/gxtb_eos_classifications.json` with a per-point
`action` (`exclude` or `retain`), `classification`, and nonempty `rationale`;
the generated `data/gxtb_eos_classification_candidates.json` is the template.

GXTB additionally has a topology gate before its quadratic fit. Across the
complete reviewed point set, the energy must decrease monotonically towards
the sampled global minimum and increase monotonically away from it, allowing
`1e-4` hartree per eight-atom LC12 cell as numerical tolerance. A lower
compressed branch plus a higher local well is recorded as
`nonmonotonic_branch`, with no fitted lattice constant. Adding points near the
local well cannot clear this status because the original reversal remains in
the persisted scale manifest; only an explicit, scientifically reviewed point
exclusion can change the set admitted to the gate. Such exclusions still need
the existing classification and nonempty rationale.

Reduced paper coverage is not supported: every one of the ten fixed systems
must have a valid EOS and all three final meshes for every method. GXTB rows
are merged additively into the dynamic working tables, while the publication
finalizer independently enforces the exact identical set. GXTB build and
protocol metadata are
kept separately in `data/build_provenance_gxtb.json`; the existing
`build_provenance.json`, `atom_energies_tblite_cli.csv`, and GFN1/GFN2 run
trees remain untouched by a selective GXTB invocation. The new atomic values
are written to `data/atom_energies_save_tblite_cli_gxtb.csv`.

GXTB inputs set `BACKUP_COPIES 0` and disable SCF restart printing. The
explicit `--prune-transients` safeguard removes only large `RESTART.kp`/WFN
transients located below a GXTB run tree and only when a successful output is
present in the same directory. It never traverses GFN1/GFN2 run trees.

Raw calculations and generated inputs are kept below `Goldzak12/runs` and
`Goldzak12/inputs` and are ignored by Git. Curated CSV/Markdown tables and
provenance are versioned in `Goldzak12/data`; `Goldzak12/figures` contains the
PNG/PDF figure used in the revised manuscript.

Literature comparison
---------------------

`scripts/plot_literature_comparison.py` augments the current EOS results with
published DFT and post-HF values. All errors are recomputed against the same
zero-point-corrected experimental values from Goldzak et al.; experimental
columns from other sources are retained only for provenance.

The comparison includes SCAN, SCAN-L, r2SCAN, and r2SCAN-L for 11 common
solids from Mejia-Rodriguez and Trickey (2020), and LSDA, PBE, PBEsol, TPSS,
revTPSS, TM, HSE06, and optB86b-vdW data from Mo et al. (2017). Coverage is
reported explicitly because the literature subsets differ.
