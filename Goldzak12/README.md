LC10 native-Bloch CP2K/tblite paper benchmark
===============================================

This directory retains the 12 cubic solids studied by Goldzak, Wang, Ye, and
Berkelbach, J. Chem. Phys. 157, 174112 (2022).  The paper benchmark is the
fixed identical ten-system subset `C, Si, SiC, BN, BP, AlN, AlP, MgS, LiF,
LiCl` for GFN1-xTB, GFN2-xTB, and g-xTB.  LiH and MgO remain diagnostic data
only and never enter a reported statistic for any of the three methods.

Complete periodic GFN1/GFN2 inputs and paper datasets are maintained in
[`Periodic-GFN2-Benchmarks`](https://github.com/DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks).
This directory retains g-xTB data, mixed comparison tables, shared references,
and a historical LC12 snapshot whose complete contents are not present in the
cited canonical source revision; see `../GFN_BASELINE_SOURCE.md`.

Publication protocol
--------------------

GFN1-xTB, GFN2-xTB, and g-xTB are rerun with the same post-#5582 CP2K binary
and the same `save_tblite` provider.  The externally reviewed
`production_ready` manifest and its independently supplied SHA256 pin are
mandatory.  Conventional cubic eight-atom cells use native Bloch sampling
with `SCHEME MACDONALD`, `SYMMETRY T`, `FULL_GRID F`, and full SPGLIB symmetry
reduction; no Born-von-Karman supercells are used.

K-point convergence is assessed independently for every method/system pair.
Separate EOS fits and equilibrium single points are first evaluated at
`k333`, `k444`, and `k555`.  The earliest single adjacent step
`n^3 -> (n+1)^3` is accepted when both
`|Delta a0| <= 0.001 A` and
`|Delta Ecoh| <= 0.05 kJ mol-1 atom-1`
(`0.000518213... eV atom-1`).  Exactly one passing interval suffices: there is
no RMS gate and no two-step requirement.  The denser `(n+1)^3` values are
reported.  An unresolved pair alone proceeds to `k666`, `k777`, `k888`,
`k999`, `k101010`, and onward until the same two criteria pass.
There is no scientific maximum mesh.  The optional integer
`--maximum-mesh N` is only a technical resource guard; reaching it without a
passing adjacent step stops with an error and never selects a value.

The historical `k333/k444/k555` energy series on a common `k444` EOS geometry
is still generated, but is explicitly a fixed-geometry diagnostic.  It never
selects `a0` or `Ecoh` and cannot establish lattice-constant convergence.
CP2K energies are extrapolated to electronic temperature `T -> 0`; isolated
atoms use the same save_tblite CLI with explicit spins and `ACCURACY 0.05`.

Only completed SCF points enter an EOS fit. A quadratic fit is rejected when
its local RMSE exceeds 0.02 hartree or when its fitted minimum lies more than
0.02 hartree above the sampled local minimum.  All ten paper systems must pass
for every method.  Existing LiH/MgO branch and multistart material is retained
for diagnostics but is not a production gate and needs only a short methods
note in the manuscript.

Current versus previous results
-------------------------------

The values below are compact copies of the earlier fixed-`k444` GFN1/GFN2
baseline used only for comparison. The complete baseline is canonical in the
GFN repository. They are not the new adaptive, same-binary LC10 paper
comparison; that table is emitted only after the complete campaign passes.

| method | fixed coverage | lattice MAE (A) | cohesive-energy MAE (eV/atom) |
|---|---:|---:|---:|
| GFN1-xTB | 10/10 | 0.145118 | 1.543851 |
| GFN2-xTB | 10/10 | 0.062410 | 1.299325 |

These values are recomputed from the existing raw result table on precisely
the same ten systems. Historical 12-system and old/new diagnostic tables are
preserved under `data/`, but are not direct paper comparisons.

Reproduction
------------

Run the publication campaign from the repository root.  Do not compute and
accept the manifest hash inside the launch command; it must be supplied from
the completed external build review.

```bash
python3 Goldzak12/scripts/run_goldzak12_k_convergence.py \
  --method GXTB \
  --campaign-manifest /path/to/build_manifest.json \
  --campaign-manifest-sha256 REVIEWED_64_HEX_SHA256 \
  --cp2k-source /path/to/clean/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/clean/save_tblite-pbc \
  --jobs 10 --threads 1 --stop-after-convergence

# Review eos_fits.csv, lc10_adaptive_k_steps.csv, every branch classification,
# and the printed fit fingerprint, then run the same command with:
#   --approve-fits

python3 Goldzak12/scripts/finalize_goldzak12_paper_summary.py
```

`--method` is repeatable and defaults to all three methods for backwards
compatibility.  The production GXTB continuation uses `--method GXTB` without
`--maximum-mesh`; it leaves the frozen external GFN1/GFN2 provenance, atom
references, and run trees untouched.  A selective run records only its
execution scope; the publication finalizer remains fail-closed and accepts
only the complete three-method by LC10 bundle.

On Terok, use the shared fail-closed MPI/affinity path.  The manifest hash must
come from the completed build qualification; do not compute and accept a new
value inside the launch command.  Repeat `--solid` and `--pe-list` as needed;
the number of ordered PE lists must equal `--jobs`, the lists must be disjoint,
and each must contain exactly one explicit logical CPU per MPI rank.  Production
MPI requires one OpenMP thread per rank.  This compact two-worker example
demonstrates the exact CLI spelling:

```bash
MANIFEST=campaigns/gxtb-pbc-v1-post5582-20260714/build_manifest.json
MANIFEST_SHA256=REPLACE_WITH_EXTERNALLY_REVIEWED_SHA256
MPI=/path/to/qualified/environment/bin/mpirun

python3 Goldzak12/scripts/run_goldzak12_k_convergence.py \
  --method GXTB \
  --campaign-manifest "$MANIFEST" \
  --campaign-manifest-sha256 "$MANIFEST_SHA256" \
  --cp2k-source /path/to/clean/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/clean/save_tblite-pbc \
  --jobs 2 --mpi-ranks-per-job 8 --threads 1 \
  --mpi-launcher "$MPI" \
  --pe-list 96,97,98,99,100,101,102,103 \
  --pe-list 104,105,106,107,108,109,110,111 \
  --eos-mesh k444 --energy-mesh k333 --energy-mesh k444 \
  --energy-mesh k555 --result-mesh k555 --stop-after-eos
```

Each logical CPU is additionally protected by a host-local `flock` for the
life of the driver. Independent production invocations therefore cannot reuse
the same CPU.  A Linux `/proc` preflight additionally rejects overlap with live
non-zombie CP2K/MPI ranks from non-locking launchers, both when the pool is
created and immediately before launch. Runtime verification is sticky across
all `/proc` samples; sequential same-rank/same-mask sanitizer PID generations
are aggregated, while concurrent duplicate ranks, rank migration, and changed
successor masks fail closed.  Resume stamps bind the direct/MPI mode, rank
count, thread settings, and exact execution-contract hash. User-supplied MPI
launcher arguments are rejected; the driver alone injects mapping, binding,
rank count, and binding reports.

Single-mesh diagnostics and g-XTB input contract
------------------------------------------------

`run_goldzak12_eos_benchmark.py` remains available for single-mesh diagnostics;
it is not the publication runner.  The g-XTB CP2K inputs pin `METHOD GXTB`,
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

After the adaptive runner has completed and its exact fit fingerprint has been
approved, freeze the publication table and its complete raw-output lineage:

```bash
python3 Goldzak12/scripts/finalize_goldzak12_paper_summary.py
```

The finalizer atomically writes `data/lc10_gfn_gxtb_paper_summary.csv`,
`data/lc10_gfn_gxtb_paper_summary.json`, and
`data/lc10_gfn_gxtb_paper_summary.tex`. The CSV has exactly one ten-system row
per method with ME, MAE, RMSE, and MaxAE for lattice constants and cohesive
energies. Each solid may have a different selected dense mesh from `k444`
upward without a fixed cap; the JSON records the raw adjacent criteria and
selected mesh. It
also adds complete EOS/equilibrium-output lineage, hashes, direct g-xTB/GFN
comparisons, and build provenance; the TeX file exports paper macros. Missing
or substituted systems, reduced coverage, stale stamps, a second/RMS gate, or
tampered raw data remove all three outputs and make finalization fail.

The runner fixes `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`,
and `VECLIB_MAXIMUM_THREADS` to 1 and sets `OMP_WAIT_POLICY=PASSIVE`; CP2K's
outer `OMP_NUM_THREADS` remains controlled by `--threads`.

`scripts/benchmark_execution.py` leaves the scientific job-stamp schema and
matcher unchanged.  It writes an additive atomic schema-v2
`*.execution.json` record that binds the launcher hash and exact command,
Open MPI's injected `--map-by pe-list=...:ordered --bind-to core
--report-bindings` policy, the hash of the preserved launcher/binding log, and
the input, output, and scientific-stamp hashes.  Rank identity comes from
`OMPI_COMM_WORLD_RANK`, not PID order; every `/proc` mask must be the singleton
CPU assigned to that rank. Sequential process generations retain their full
PID history in the record and may be combined only when rank and singleton mask
are unchanged and no two generation PIDs are concurrently live.  Binding/mapping CLI overrides, `--bind-to none`,
outer taskset, duplicate CPUs, wrong list lengths, unavailable CPUs, and
multi-threaded MPI ranks fail closed.  Schema-v1 shared-taskset records remain
readable and their numerical energies/forces/stresses retain provenance, but
their timings are explicitly `legacy_timing_non_scaling` and cannot enter a
scaling comparison. On an SMT host where `--bind-to core` exposes more than
one logical CPU in a rank mask, the singleton gate also fails closed; such a
machine needs a separately qualified binding policy rather than an exception.

The separate fixed-geometry diagnostic single points are never launched
implicitly. `--stop-after-convergence` completes the adaptive independent-EOS
stage and stops; `--fit-only` recollects already stamped outputs without
launching an executable; and only `--approve-fits` records the exact current
multi-mesh fit-table fingerprint and permits the diagnostic series. For the
38-GB development Mac,
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

The initial same-binary three-method publication plan contains 33 save_tblite
atom jobs, 990 EOS jobs (3 methods x 10 solids x 3 meshes x 11 scales), and 90
own-minimum equilibrium single points.  After approval, the separate
fixed-geometry diagnostic adds 90 single points.  The independent g-XTB CP2K
atom gate adds 11 jobs, for 1214 calculations before any adaptive `k666+`
extension.  Each unresolved method/solid/mesh adds 11 EOS points plus one
own-minimum single point.  A targeted extra EOS scale adds one job at each mesh
where it is requested.

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
  --pe-list 112,113,114,115,116,117,118,119 \
  --pe-list 120,121,122,123,124,125,126,127

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

Equilibrium and fixed-diagnostic inputs are runnable only after a valid
quadratic EOS fit. Each generated input has an adjacent `*.inp.eos.json`
lineage record containing the fitted lattice constant, EOS mesh, energy mesh,
input hash, and SPGLIB contract. A pre-generated input at the experimental
lattice constant, an input without this valid lineage, or a lineage attached
to an invalid fit is stale. Every adaptive Ecoh value is evaluated at its own
mesh's EOS minimum; the fixed-geometry series is stored in a different run
tree and marked ineligible for paper selection.

The runner never forces a fit. A missing or discontinuous EOS minimum for any
of the ten paper systems stops before final single points and is recorded in
`data/gxtb_adaptive_followup.csv`/`.md`; the paper artifact still requires
10/10. LiH/MgO diagnostic behavior has no effect on this gate.

The exact method/system/mesh grid is persisted in
`data/lc10_k_convergence_scale_manifest.json` and used by the finalizer. A
large local energy discontinuity is reported as a numerical SCC
branch *candidate*, never immediately as a physical failure. Starting with the
second available mesh, the runner also compares the pointwise
`E_N(scale)-E_(N-1)(scale)` correction with its robust scale trend. An isolated
cross-mesh residual of at least `0.01` hartree (and eight robust MADs) is
reported fail-closed; this catches a converged alternate SCC root that can look
locally smooth within one EOS curve. The cross-mesh gate never excludes a
point automatically. Each candidate
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
