# DMC-ICE13 Periodic GFN Benchmark

This directory contains CP2K/tblite single-point calculations for the
DMC-ICE13 ice polymorph benchmark. The calculations compare periodic
GFN1-xTB and GFN2-xTB relative energies against the diffusion Monte Carlo
reference values of Della Pia, Zen, Alfe, and Michaelides,
J. Chem. Phys. 157, 134701 (2022), DOI: 10.1063/5.0102645.

## Data included

- `poscars/`: POSCAR geometries for the 13 DMC-ICE13 polymorphs.
- `inputs/`: Gamma-only CP2K input files for GFN1-xTB and GFN2-xTB.
- `kpoint_inputs/`: explicit 1x1x1, 2x2x2, 3x3x3, 4x4x4, and 5x5x5
  MacDonald k-point CP2K input files.
- `runs/`: Gamma-only CP2K inputs and outputs for each method and polymorph.
- `runs_kpoints/`: CP2K inputs and outputs for the k-point benchmark.
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
- `data/dmc_ice13_relative_mae_comparison.csv`: comparison with the published
  DFT data from the DMC-ICE13 paper.
- `data/dmc_ice13_published_dft_absolute_energies.csv`: published DMC and DFT
  absolute lattice energies from the DMC-ICE13 paper, used to compute the
  relative-energy MAE ranking.
- `figures/`: SVG and PNG plots generated from the benchmark data.
- `scripts/`: extraction, analysis, plotting, and remote run scripts.

The original PDF and Supporting Information are not redistributed here. The
geometries and DMC reference values are documented through the paper DOI above.

## CP2K setup used

The calculations were run with a CP2K 2026.1 development build interfaced to
tblite. The current Spark default is available as
`/home/kuehne88/bin/cp2k-current-tblite.psmp`:

- CP2K source revision: `0622d442e4`
- CP2K flags: `omp fftw3 libxc parallel scalapack mpi_f08 xsmm spglib
  libdftd4 dftd4_v3 mctc-lib tblite`
- tblite: `0.5.0`
- `TBLITE/ACCURACY`: `0.1`
- `EPS_SCF`: `1.0E-9`
- `OMP_NUM_THREADS`: `4`

The primary comparison uses the Gamma-centered 3x3x3 k-point mesh, matching
the non-hybrid DFT single-point setup in the DMC-ICE13 reference. The explicit
1x1x1 mesh verifies equivalence to the Gamma-only calculation, the 2x2x2 mesh
documents the approach to convergence, and the 4x4x4 and 5x5x5 checks confirm
that the 3x3x3 aggregate statistics are converged. All energies in the CSV
summaries are relative to ice Ih and reported in kJ mol-1 per water molecule.
