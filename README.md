# Periodic g-xTB Benchmarks

This public repository collects the paper-relevant inputs, validation gates,
curated output data, analysis scripts, and figures for the periodic g-xTB
implementation in CP2K with save_tblite.

Complete GFN1-xTB/GFN2-xTB inputs and result datasets are not duplicated here.
Their canonical public home is
[`DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks`](https://github.com/DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks).
Only compact values needed for explicit g-xTB comparisons are retained; the
verified split and source revision are recorded in
[`GFN_BASELINE_SOURCE.md`](GFN_BASELINE_SOURCE.md). All g-xTB workflows,
inputs, results, and Part-I/Part-II paper artifacts belong exclusively here.

## Contents

- `DMC-ICE13/`: phase-wise adaptive native-Bloch k-point convergence for
  periodic g-xTB relative ice-polymorph energies, with compact GFN1/GFN2 and
  diffusion Monte Carlo comparison values.
- `X23b/`: correctness-gated g-xTB molecular-crystal workflow covering
  molecular optimization, shifted-k preflights, force/stress finite
  differences, native-Bloch cell optimization, and final k-point energies.
- `Goldzak12/`: LC12 source calculations and the fixed LC10 paper subset
  (C, Si, SiC, BN, BP, AlN, AlP, MgS, LiF, and LiCl), including equations of
  state, cohesive energies, literature comparisons, and SCC-root diagnostics.
- `campaigns/`: immutable g-xTB build, protocol, and qualification manifests.
- `validation/`: molecular, primitive-cell/supercell, symmetry, force, and
  stress regression inputs for the CP2K/save_tblite bridge. The canonical
  artifact digests removed from the manuscript and Supplementary Material are
  stored in `validation/paper_artifact_sha256.json`. The complete post-#5582
  analytical-force/stress qualification, including raw outputs and its
  discarded CPU-affinity diagnostic, is frozen in
  `validation/gxtb_derivative_regtests_post5582_20260716/`. The corresponding
  full-grid/K290/SPGLIB, force, and periodic-virial gates for representative
  one- and two-dimensional PBC are in
  `validation/gxtb_partial_pbc_post5582_20260716/`. The same 23 low-k,
  force, stress, K290/SPGLIB, and partial-PBC cases were repeated with the
  final clean DMC build in
  `validation/gxtb_final_lowk_derivatives_20260719/`. This repetition closes
  the formerly open Born--von Karman supercell residual: the one-dimensional
  difference is zero at printed precision and the two-dimensional difference
  is `4.070e-9` hartree per primitive cell.
- The final-provider component and source attribution is archived in
  `validation/provider_component_attribution_20260719/` and
  `validation/pbc_h0_anisotropy_attribution_20260719/`.  The controlled source
  A/B reproduces the final `pbc` phase-VII energy within `1.91e-11` hartree and
  assigns the residual provider difference to the historical central-cell H0
  anisotropy.  An equivalent lattice-image test supports retaining the current
  image-complete periodic H0 definition.
- Exchange/ACP energy and derivative ablations for the final CP2K build are in
  `validation/dmc13_k222_viii_component_ablation_20260719/` and
  `validation/dmc13_k222_xvii_derivative_component_ablation_20260719/`.
- Direct source-level qualification of the final `save_tblite` periodic H0,
  Wigner--Seitz, exchange, Fock-response, force, stress, and transform paths is
  retained in `validation/save_tblite_periodic_source_tests_20260719/`.  The
  periodic subsets pass completely; the sole 74/75 Hamiltonian-suite outlier
  is a nonperiodic CeCl3 finite-difference threshold case reproduced with
  identical component errors by the historical final `pbc` source.
- `patches/`: local CP2K and tblite patches used for the final benchmark
  revision.
- `scripts/`: helper scripts used for the final k-point, cell-optimization,
  and CP2K-native-vs-tblite-CLI checks.
- `scripts/finalize_paper_benchmark_bundle.py`: fail-closed aggregation of the
  three completed benchmark-specific publication bundles into one CSV, JSON
  lineage record, and set of TeX number macros under `paper/`. The JSON and
  TeX exports additionally contain like-for-like g-XTB-minus-GFN1 and
  g-XTB-minus-GFN2 deltas, error ratios, and percentage changes. The LC result
  is emitted only for the exact LC10 set shared by GFN1-xTB, GFN2-xTB, and
  g-XTB; method-specific LC12 coverage is never exported or compared.
- `CODE_PATCHES.md`: compact implementation provenance retained where it is
  also relevant to the g-xTB integration.

Routine generated CP2K working directories, raw standard-output files, and
optional diagnostic plots are not tracked. Explicitly named immutable
qualification archives are the exception: they retain the raw evidence needed
to reproduce a paper table or diagnose an invalidated launch. All other runs
can be recreated from the versioned inputs and scripts.

## Current g-xTB campaign

The active campaign is deliberately fail-closed: provisional numbers are not
promoted to final paper summaries until the stored build identity, raw
artifacts, k-point convergence, symmetry, force/stress, and cross-build gates
pass.  The current same-build implementation qualification uses CP2K source
revision `8520b2e592cd04d35081ab4ad46d92c606071e23` and `save_tblite` source
revision `15915c9435644eb257178ca8f8bf7220c38b1a84`; the corresponding Linux
executables have SHA-256 digests
`b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f`
and `f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a`.
The complete implementation audit and exact provenance are retained under
`DMC-ICE13/reproduction/save_tblite_direct_dmc13/validation/`.

The older phase-wise DMC-ICE13 trace remains frozen only as a pre-response
baseline.  Its response-corrected replacement is promoted only after every
phase has a qualifying adjacent-mesh endpoint from the same immutable build.
LC10 remains explicitly provisional until all ten solids pass the per-solid
lattice-constant and cohesive-energy gates.  Other retained files outside the
current manuscript benchmark scope are technical regression material and not
accuracy claims.  The publication finalizers refuse incomplete or unhashed
data.

The older three-benchmark bundle command remains available for archival
campaigns:

```bash
python3 scripts/finalize_paper_benchmark_bundle.py
```

The command deletes stale aggregate outputs and fails without emitting a
replacement if any child bundle is incomplete, has changed hashes, lacks the
three-method comparison, or has inconsistent coverage.

## External GFN1/GFN2 comparison source

The comparison calculations use DCM-Uni-Paderborn CP2K development trunk revision
`faf9aae91266170dfee8a9f7171a5135bc5eb368` with tblite support. The tblite
build combines `main` revision `eb50bbfbe1c0869e2e18c9b7cc13144e5130b6df`
with PR 350 head `8c5e56255dc0f7001615489f24162ed770888d8b` in local merge
`8a9d09474b93d25c044d6f46ce920750c7fe4cf7`; PR 343 is already in the base.
The frozen CP2K and tblite executable SHA-256 hashes are
`f2b8e6e516b60d49af722997dd0bf06c10b54b2a2a221f786e5eaea38cccd8a5`
and `d50145af569a6ce4ea4e73e68d1cb004c3ca240105deb941c0244b7d431ed47f`.

The raw inputs and complete result tables for these values are in the canonical
GFN repository linked above. Primary aggregate comparison values retained here
are:

| Benchmark | Setup | Method | MAE |
|---|---|---|---:|
| DMC-ICE13 relative energies | native Bloch 3x3x3 | GFN1-xTB | 8.005255 kJ mol-1 |
| DMC-ICE13 relative energies | native Bloch 3x3x3 | GFN2-xTB | 3.462919 kJ mol-1 |
| X23b lattice energies | k333 SP on native Bloch k222 cell opt | GFN1-xTB | 11.345702 kJ mol-1 |
| X23b lattice energies | k333 SP on native Bloch k222 cell opt | GFN2-xTB | 14.092104 kJ mol-1 |
| X23b cell volumes | native Bloch k222 cell opt | GFN1-xTB | 7.514116 percent |
| X23b cell volumes | native Bloch k222 cell opt | GFN2-xTB | 5.842296 percent |
| LC10 lattice constants | exact common set, k444 EOS | GFN1-xTB | 0.145118 A |
| LC10 lattice constants | exact common set, k444 EOS | GFN2-xTB | 0.062410 A |
| LC10 cohesive energies | exact common set, k555 on k444 EOS minima | GFN1-xTB | 1.543851 eV atom-1 |
| LC10 cohesive energies | exact common set, k555 on k444 EOS minima | GFN2-xTB | 1.299325 eV atom-1 |

All production k-point calculations use native Bloch sampling with full
SPGLIB symmetry reduction. The completed production counts are 156/156 for
DMC-ICE13, 46/46 for X23b k222 cell optimization, and 46/46 each for the
X23b k333 and k444 final-geometry single points.

The new LC10 paper table is not populated from those historical rows.  All
three methods share one post-#5582 CP2K/save_tblite build and independently
converge both a0 and Ecoh per solid from k333 upward without a fixed maximum
mesh. One passing adjacent step is sufficient, both thresholds must pass, and
the denser value is retained; no RMS or two-step gate is applied.
