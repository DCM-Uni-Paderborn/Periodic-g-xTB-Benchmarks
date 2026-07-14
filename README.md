# Periodic g-xTB Benchmarks

This private repository collects the paper-relevant inputs, validation gates,
curated output data, analysis scripts, and figures for the periodic g-xTB
implementation in CP2K with save_tblite.

Frozen GFN1-xTB/GFN2-xTB values are retained only as comparison baselines.
Their canonical public home remains
[`DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks`](https://github.com/DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks).
New g-xTB workflows, inputs, results, and paper artifacts belong exclusively
in this repository.

## Contents

- `DMC-ICE13/`: phase-wise adaptive native-Bloch k-point convergence for
  periodic g-xTB relative ice-polymorph energies, with frozen GFN1/GFN2 and
  diffusion Monte Carlo comparison data.
- `X23b/`: correctness-gated g-xTB molecular-crystal workflow covering
  molecular optimization, shifted-k preflights, force/stress finite
  differences, native-Bloch cell optimization, and final k-point energies.
- `Goldzak12/`: LC12 equations of state, cohesive energies, literature
  comparisons, SCC-root/hysteresis diagnostics, and lattice-constant fits for
  cubic covalent and ionic solids.
- `campaigns/`: immutable g-xTB build, protocol, and qualification manifests.
- `validation/`: molecular, primitive-cell/supercell, symmetry, force, and
  stress regression inputs for the CP2K/save_tblite bridge.
- `patches/`: local CP2K and tblite patches used for the final benchmark
  revision.
- `scripts/`: helper scripts used for the final k-point, cell-optimization,
  and CP2K-native-vs-tblite-CLI checks.
- `scripts/finalize_paper_benchmark_bundle.py`: fail-closed aggregation of the
  three completed benchmark-specific publication bundles into one CSV, JSON
  lineage record, and set of TeX number macros under `paper/`. The JSON and
  TeX exports additionally contain like-for-like g-XTB-minus-GFN2 deltas,
  error ratios, and percentage changes. LC12 comparisons are emitted only for
  the identical three-method common subset, never across differing coverage.
- `FINAL_RESULTS.md`, `CODE_PATCHES.md`, and `paper_revision_numbers.csv`:
  compact provenance for the current paper revision.

Generated CP2K working directories, raw standard-output files, and optional
diagnostic plots are not tracked. They can be recreated from the versioned
inputs and scripts; the curated CSV, JSON, plotting data, and manuscript figure
files are the benchmark data used in the paper.

## Current g-xTB campaign

The active campaign is deliberately fail-closed: provisional numbers are not
promoted to paper summaries until the stored build identity, raw artifacts,
k-point convergence, symmetry, force/stress, and cross-build gates pass.
Current production calculations use CP2K source revision
`18d37c946413...` with save_tblite revision `1449febde312...`; the direct
ACP k-space successor build is still in qualification. See
`campaigns/gxtb-pbc-v1-20260714/` for the complete machine-readable identity.

The current g-xTB DMC-ICE13 and LC12 values are therefore explicitly
provisional, and X23b has no g-xTB accuracy result yet. The publication
finalizers refuse incomplete or unhashed data.

After the DMC-ICE13, X23b, and LC12 finalizers have all succeeded, create the
single manuscript-facing bundle with

```bash
python3 scripts/finalize_paper_benchmark_bundle.py
```

The command deletes stale aggregate outputs and fails without emitting a
replacement if any child bundle is incomplete, has changed hashes, lacks the
three-method comparison, or has inconsistent coverage.

## Frozen GFN1/GFN2 comparison snapshot

The imported comparison calculations use DCM-Uni-Paderborn CP2K development trunk revision
`faf9aae91266170dfee8a9f7171a5135bc5eb368` with tblite support. The tblite
build combines `main` revision `eb50bbfbe1c0869e2e18c9b7cc13144e5130b6df`
with PR 350 head `8c5e56255dc0f7001615489f24162ed770888d8b` in local merge
`8a9d09474b93d25c044d6f46ce920750c7fe4cf7`; PR 343 is already in the base.
The frozen CP2K and tblite executable SHA-256 hashes are
`f2b8e6e516b60d49af722997dd0bf06c10b54b2a2a221f786e5eaea38cccd8a5`
and `d50145af569a6ce4ea4e73e68d1cb004c3ca240105deb941c0244b7d431ed47f`.

Primary aggregate results:

| Benchmark | Setup | Method | MAE |
|---|---|---|---:|
| DMC-ICE13 relative energies | native Bloch 3x3x3 | GFN1-xTB | 8.005255 kJ mol-1 |
| DMC-ICE13 relative energies | native Bloch 3x3x3 | GFN2-xTB | 3.462919 kJ mol-1 |
| X23b lattice energies | k333 SP on native Bloch k222 cell opt | GFN1-xTB | 11.345702 kJ mol-1 |
| X23b lattice energies | k333 SP on native Bloch k222 cell opt | GFN2-xTB | 14.092104 kJ mol-1 |
| X23b cell volumes | native Bloch k222 cell opt | GFN1-xTB | 7.514116 percent |
| X23b cell volumes | native Bloch k222 cell opt | GFN2-xTB | 5.842296 percent |
| LC12 lattice constants | k444 EOS, 12/12 | GFN1-xTB | 0.136650 A |
| LC12 lattice constants | k444 EOS, 10/12 | GFN2-xTB | 0.062410 A |
| LC12 cohesive energies | k555 on k444 EOS minima, 12/12 | GFN1-xTB | 1.457694 eV atom-1 |
| LC12 cohesive energies | k555 on k444 EOS minima, 10/12 | GFN2-xTB | 1.299325 eV atom-1 |

All production k-point calculations use native Bloch sampling with full
SPGLIB symmetry reduction. The completed production counts are 156/156 for
DMC-ICE13, 46/46 for X23b k222 cell optimization, and 46/46 each for the
X23b k333 and k444 final-geometry single points.
