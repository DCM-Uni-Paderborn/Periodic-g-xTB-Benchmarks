# Code Patch Inventory

## tblite

Patch: `patches/tblite_wsc_multipole_ewald_local.patch`

Base/HEAD: `5b14b8430bb2ffb3c96808466ad670821f81f745`.

Files changed:

- `src/tblite/coulomb/multipole.f90`
- `src/tblite/coulomb/ewald.f90`
- `src/tblite/wignerseitz.f90`
- `src/tblite/cutoff.f90`

Purpose:

- Correct Wigner-Seitz image indexing and image weighting for multipolar electrostatics.
- Use multipole-aware Ewald real/reciprocal cutoff estimates.
- Use WSC images consistently in multipole matrix and gradient/virial paths.
- Respect directional periodicity masks for cutoff and central-cell wrapping.

## CP2K

Patch: `patches/cp2k_tblite_interface_local.patch`

Base/HEAD: `518a50992f009b083c127372f294e6485306c05b`.

File changed:

- `src/tblite_interface.F`

Purpose:

- Forward mixer settings into `tb%calc%mixer_input` for the current tblite API.

## Benchmark Scripts

Full helper scripts used for the final revision are in `scripts/`.
Some runner defaults preserve the local production paths used for the paper
revision; override `--benchmark-root`, `--out`, or `CP2K` when replaying them
in a different checkout. `scripts/update_x23b_k222_figures.py` is repo-relative
and regenerates the X23b figures from the versioned `X23b/data` files.

- `run_dmc13_kpoint_jobs.py`: DMC13 native Bloch k-point benchmark runner.
- `run_x23b_cellopt_variant_matrix.py`: X23b Gamma cellopt variant runner, including `cg_2pnt_keep_angles`.
- `run_x23b_cellopt_final_kpoint_sp.py`: X23b final-cellopt native Bloch k222/k333 single-point runner.
- `run_x23b_reference_cli_checks.py`: CP2K-native vs tblite CLI reference comparison.
