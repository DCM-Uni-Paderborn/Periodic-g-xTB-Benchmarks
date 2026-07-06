# Periodic GFN2 Benchmarks

This repository collects the paper-relevant benchmark inputs, curated output
data, analysis scripts, and figures for periodic GFN calculations in CP2K.

## Contents

- `DMC-ICE13/`: CP2K/tblite single-point benchmark for the DMC-ICE13 ice
  polymorph data set, comparing periodic GFN1-xTB and GFN2-xTB relative
  energies with the diffusion Monte Carlo reference energies. The current
  manuscript data use native Bloch k-point calculations through CP2K
  `&KPOINTS`.
- `X23b/`: CP2K/tblite molecular-crystal benchmark with gas-phase molecular
  optimizations, crystal single-point k-point tests, native Bloch 2x2x2
  crystal cell optimizations, extracted X23b lattice energies, volume errors,
  and summary plots.
- `patches/`: local CP2K and tblite patches used for the final benchmark
  revision.
- `scripts/`: helper scripts used for the final k-point, cell-optimization,
  and CP2K-native-vs-tblite-CLI checks.
- `FINAL_RESULTS.md`, `CODE_PATCHES.md`, and `paper_revision_numbers.csv`:
  compact provenance for the current paper revision.

Generated CP2K working directories and raw standard-output files are not
tracked. They can be recreated from the versioned inputs and scripts; the
curated CSV, JSON, and plotting data files are the benchmark data used in the
manuscript.

## Current revision snapshot

The final calculations use CP2K trunk revision
`518a50992f009b083c127372f294e6485306c05b` with tblite support and tblite
revision `5b14b8430bb2ffb3c96808466ad670821f81f745` (`tblite` 0.6.0),
including the changes corresponding to tblite PRs 343 and 350.

Primary aggregate results:

| Benchmark | Setup | Method | MAE |
|---|---|---|---:|
| DMC-ICE13 relative energies | native Bloch 3x3x3 | GFN1-xTB | 8.008187 kJ mol-1 |
| DMC-ICE13 relative energies | native Bloch 3x3x3 | GFN2-xTB | 3.185301 kJ mol-1 |
| X23b lattice energies | native Bloch 2x2x2 cell opt | GFN1-xTB | 11.129018 kJ mol-1 |
| X23b lattice energies | native Bloch 2x2x2 cell opt | GFN2-xTB | 14.459836 kJ mol-1 |
| X23b cell volumes | native Bloch 2x2x2 cell opt | GFN1-xTB | 7.914787 percent |
| X23b cell volumes | native Bloch 2x2x2 cell opt | GFN2-xTB | 5.616637 percent |
