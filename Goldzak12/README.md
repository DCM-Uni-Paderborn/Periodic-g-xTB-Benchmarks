LC12 (Goldzak12) native-Bloch CP2K/tblite benchmark
===================================================

This directory contains the 12 cubic covalent and ionic solids studied by
Goldzak, Wang, Ye, and Berkelbach, J. Chem. Phys. 157, 174112 (2022). It
compares CP2K/tblite GFN1-xTB and GFN2-xTB and CP2K/save_tblite g-XTB with the
reported HF, MP2, SCS-MP2, SOS-MP2, and zero-point-corrected experimental
lattice constants and cohesive energies.

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
0.02 hartree above the sampled local minimum. The current run gives 12/12
valid GFN1 fits and 10/12 valid GFN2 fits. GFN2/MgO has no bracketed stable
minimum on the compressed branch, while GFN2/LiH has a discontinuous EOS and
fails the general fit-quality criterion.

Volume continuation with separate CP2K Bloch-wavefunction and native tblite
SCC restarts removes the earlier independent-start failures. LiH converges for
all 32 sampled points down to scale 0.71, but its energy continues to decrease
until the electronic branch collapses below that range. MgO can be followed in
fine steps to scale 0.926, where its energy is still decreasing; the additional
0.90 and 0.88 points enter the same charge-collapse branch even after damped
2400-step SCC retries. The 10/12 coverage therefore reflects missing physical
EOS minima, not unfinished production jobs.

Current versus previous results
-------------------------------

The values in this section are the frozen GFN1/GFN2 production results. Adding
GXTB is method-selective and does not replace these rows or their raw outputs.

| method | coverage | lattice MAE (A) | cohesive-energy MAE (eV/atom) |
|---|---:|---:|---:|
| GFN1 current | 12/12 | 0.136650 | 1.457694 |
| GFN1 previous | 12/12 | 0.164341 | 1.457325 |
| GFN2 current | 10/12 | 0.062410 | 1.299325 |
| GFN2 previous | 11/12 | 0.147638 | 1.731839 |

On the identical ten-system GFN2 subset, the lattice-constant MAE decreases
from 0.133264 to 0.062410 A and the cohesive-energy MAE decreases from 1.534521
to 1.299325 eV/atom. The frozen previous tables are in
`data/baseline_20260710`; `data/old_vs_new.md` and the associated CSV files
contain the complete per-system comparison.

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

GXTB production extension
-------------------------

GXTB is run selectively with the save_tblite-enabled CP2K executable and the
matching save_tblite CLI. The generated CP2K inputs pin `METHOD GXTB`,
`SCC_MIXER TBLITE` (the native g-XTB FDIIS potential/Fock mixer), and CP2K
`DIRECT_P_MIXING`; no CP2K-Fock, CP2K-density, or modified-Broyden retry is
used as an alternative GXTB production mixer. The 13 isolated atoms use
`save_tblite run --method gxtb --spin 2S`, where `2S = multiplicity - 1` is
recorded alongside every energy.

All three methods now use the identical native-Bloch k-point contract:
`SCHEME MACDONALD`, `SYMMETRY T`, `FULL_GRID F`, and the SPGLIB backend and
reduction method. For GXTB, CP2K expands the irreducible density/overlap data
to the complete mesh before the coupled save_tblite evaluation and folds the
response back afterwards. Older GXTB inputs or outputs with `SYMMETRY F` and
`FULL_GRID T` are diagnostics only and are never accepted as LC12 production
data.

LC12 obtains its EOS and cohesive energies entirely from `RUN_TYPE ENERGY`.
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
  --campaign-manifest campaigns/gxtb-pbc-v1-20260714/build_manifest.json \
  --cp2k-source /path/to/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/save_tblite-pbc \
  --jobs 3 --threads 1 --eos-mesh k444 \
  --energy-mesh k333 --energy-mesh k444 --energy-mesh k555 \
  --result-mesh k555 --stop-after-eos

# After reviewing eos_fits.csv, gxtb_eos_branch_diagnostics.csv, and any
# classification/adaptive follow-up, approve this exact fit fingerprint:
python3 Goldzak12/scripts/run_goldzak12_eos_benchmark.py \
  --method GXTB \
  --campaign-manifest campaigns/gxtb-pbc-v1-20260714/build_manifest.json \
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

The finalizer writes `data/lc12_gfn_gxtb_paper_summary.csv` and
`data/lc12_gfn_gxtb_paper_summary.json` atomically.  The CSV has one reported-
coverage and one identical-three-method-common-subset row per method, including
ME, MAE, RMSE, and MaxAE for both lattice constants and cohesive energies.  The
JSON additionally records every accepted EOS point, every `k333/k444/k555`
single point, isolated-atom references, input/output/stamp hashes, the exact
fit approval, and both build-provenance records.  It removes stale publication
files and fails without replacement when a raw energy, input lineage, campaign
stamp, fit fingerprint, or required mesh is incomplete.

Full 12/12 GXTB coverage is the default.  A reduced GXTB set can be frozen only
when the already approved GXTB provenance explicitly records
`allow_reduced_coverage`, the approved minimum is met, and every omitted system
has a completed adaptive investigation.  Such a system remains in the JSON
with its fit status, interpretation, EOS raw hashes, and invalidated final-input
lineage, but it receives no fabricated lattice or cohesive-energy result.  Run
the validator with the same `--allow-reduced-coverage` and
`--minimum-valid-fits` values before finalizing such a campaign.

`k555` is the runner and validator default result mesh, matching the frozen
GFN1/GFN2 report. The runner fixes `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`,
and `VECLIB_MAXIMUM_THREADS` to 1 and sets `OMP_WAIT_POLICY=PASSIVE`; CP2K's
outer `OMP_NUM_THREADS` remains controlled by `--threads`.

GXTB final single points are never launched implicitly. `--stop-after-eos`
runs the atom/EOS stage and stops after writing the fits; `--fit-only`
recollects already stamped EOS outputs without launching an executable; and
only `--approve-fits` records the exact current fit-table fingerprint and
permits the `k333/k444/k555` stage. For the 38-GB development Mac,
`--jobs 3 --threads 1` is the conservative unmeasured default because the final queue
contains `k555` jobs; increase concurrency only after measuring their resident
set size. The 13 isolated-atom checks are small and can separately use more
workers.

The versioned `campaigns/gxtb-pbc-v1-20260714/build_manifest.json` is the
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
CP2K-versus-save_tblite check (13 atoms, no solid production jobs):

```bash
python3 Goldzak12/scripts/run_goldzak12_benchmark.py atom-check \
  --method GXTB \
  --campaign-manifest campaigns/gxtb-pbc-v1-20260714/build_manifest.json \
  --cp2k-source /path/to/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/save_tblite-pbc \
  --jobs 10 --threads 1 --tolerance-hartree 1e-6
```

This writes `data/atom_reference_cp2k_vs_save_tblite_gxtb.csv` and fails if an
atom is missing or exceeds the selected tolerance.

The default GXTB production plan is 13 save_tblite atom jobs, 132 CP2K k444 EOS jobs
(12 solids times 11 standard scales), and three final single points per valid
EOS minimum. Thus a 12/12 run has 181 jobs in total (168 CP2K plus 13 CLI);
with `n` valid minima it has `145 + 3n` jobs. Every
`--adaptive-scale SOLID=SCALE` adds one CP2K EOS job. A targeted point can also
be added through the stamped main runner:

```bash
python3 Goldzak12/scripts/run_goldzak12_eos_benchmark.py \
  --method GXTB --adaptive-scale MgO=0.92000 --adaptive-scale MgO=0.90000 \
  --campaign-manifest campaigns/gxtb-pbc-v1-20260714/build_manifest.json \
  --cp2k-source /path/to/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/save_tblite-pbc --stop-after-eos
```

The separate continuation script may still be used for diagnostics, but GXTB
`--promote` is rejected: an independently copied continuation output cannot
satisfy the canonical input plus complete campaign-stamp contract.

LiH/MgO multi-start branch qualification
-----------------------------------------

LiH and MgO are not eligible for an EOS fit until the versioned multi-start
map in `data/gxtb_multistart_plan.json` has been completed and classified.  It
contains 18 LiH and 20 MgO scales.  At every scale the protocol requests an
independent cold start plus ascending and descending WFN-continuation chains
(52 and 58 calculations, respectively).  A failed chain is archived and
stopped; it is never retried implicitly.  Every input pins the reduced shifted
`k444` SPGLIB contract, native save_tblite FDIIS, and explicit restart logging.

The runner rejects a CP2K source that is not descended from upstream commit
`c92cc08b45378b85150447011b5a4bb552f5b797` (merged PR #5582).  In particular,
the earlier `18d37c` build remains frozen diagnostic provenance and must not be
used for this map.  Pre-production execution also needs an explicit state
argument; the production default remains `production_ready`:

```bash
python3 Goldzak12/scripts/run_gxtb_multistart_branches.py \
  --campaign-manifest /path/to/post-c92cc08/build_manifest.json \
  --cp2k-source /path/to/clean/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/clean/save_tblite-pbc \
  --campaign-state qualification_pending \
  --cold-workers 8 --threads 1

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
output into `runs/eos`, never approves a fit, and never starts `k555`.  Manual
promotion requires an identity-identical `production_ready` manifest and a
second hash review.

The independent CP2K-versus-CLI atom check adds 13 CP2K atom jobs. Counting
that mandatory acceptance gate, a full 12/12 publication campaign therefore
contains 194 calculations before adaptive additions (`158 + 3n` for `n`
valid fits), not 181.

Final-single-point inputs are runnable only after a valid quadratic EOS fit.
Each generated input has an adjacent `*.inp.eos.json` lineage record containing
the fitted lattice constant, EOS mesh, energy mesh, input hash, and SPGLIB
contract. A pre-generated input at the experimental lattice constant, an input
without this valid lineage, or a lineage attached to an invalid fit is treated
as stale. Once the fit is accepted, the runner regenerates the input at the
actual EOS minimum before CP2K can run it.

Missing or discontinuous EOS minima remain explicit reduced coverage; the
runner does not force a fit. It first writes
`data/gxtb_adaptive_followup.csv`/`.md` with suggested scales and stops before
final single points. Add the suggested points, then rerun, for example:

```bash
python3 Goldzak12/scripts/run_goldzak12_eos_benchmark.py \
  --method GXTB \
  --campaign-manifest campaigns/gxtb-pbc-v1-20260714/build_manifest.json \
  --cp2k-source /path/to/cp2k-g-xTB-pbc \
  --save-tblite-source /path/to/save_tblite-pbc \
  --adaptive-scale MgO=0.91000 --adaptive-scale MgO=0.93000 \
  --stop-after-eos
```

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

Only after adaptive investigation may an incomplete, but still meaningful,
subset be accepted explicitly (default minimum 8 of 12 quadratic fits):

```bash
python3 Goldzak12/scripts/run_goldzak12_eos_benchmark.py \
  --method GXTB --allow-reduced-coverage --minimum-valid-fits 8 [same run options]
python3 Goldzak12/scripts/validate_goldzak12_results.py \
  --method GXTB --allow-reduced-coverage --minimum-valid-fits 8
```

Without this flag both runner and validator require 12/12 valid GXTB fits;
zero or merely incidental fit coverage can no longer pass. Every requested
EOS point must be completed or explicitly classified with a rationale; a
classified failed point still makes its runner invocation return nonzero and
can only enter a separately validated reduced-coverage result. GXTB rows are
merged additively into the dynamic
three-method tables, and a separate common-valid-subset table is written to
`data/eos_common_subset_summary.csv`. GXTB build and protocol metadata are
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
