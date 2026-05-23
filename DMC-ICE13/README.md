# DMC-ICE13 Periodic GFN Benchmark

This directory contains CP2K/tblite single-point calculations for the
DMC-ICE13 ice polymorph benchmark. The calculations compare periodic
GFN1-xTB and GFN2-xTB relative energies against the diffusion Monte Carlo
reference values of Della Pia, Zen, Alfe, and Michaelides,
J. Chem. Phys. 157, 134701 (2022), DOI: 10.1063/5.0102645.

## Data included

- `poscars/`: POSCAR geometries for the 13 DMC-ICE13 polymorphs.
- `inputs/`: CP2K input files for GFN1-xTB and GFN2-xTB.
- `runs/`: CP2K inputs and outputs for each method and polymorph.
- `data/results.json`: raw CP2K total energies, per-water energies, relative
  energies with respect to ice Ih, and error statistics.
- `data/dmc_ice13_relative_energies.csv`: relative energies and GFN errors.
- `data/dmc_ice13_relative_mae_comparison.csv`: comparison with selected
  published DFT data from the DMC-ICE13 paper.
- `figures/`: SVG and PNG plots generated from the benchmark data.
- `scripts/`: extraction, analysis, plotting, and remote run scripts.

The original PDF and Supporting Information are not redistributed here. The
geometries and DMC reference values are documented through the paper DOI above.

## CP2K setup used

The calculations were run with a CP2K 2026.1 development build interfaced to
tblite:

- CP2K source revision: `95e0cafc31`
- CP2K flags: `omp libint fftw3 libxc parallel scalapack mpi_f08 xsmm spglib
  libdftd4 dftd4_v3 mctc-lib tblite`
- tblite: `0.5.0`
- `TBLITE/ACCURACY`: `0.1`
- `EPS_SCF`: `1.0E-9`
- `OMP_NUM_THREADS`: `4`

All energies in the CSV summaries are relative to ice Ih and reported in
kJ mol-1 per water molecule.

