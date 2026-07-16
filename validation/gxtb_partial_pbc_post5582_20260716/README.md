# Post-#5582 one- and two-dimensional PBC qualification

This archive freezes the 1D/2D periodic g-xTB qualification performed on
`terok.casus.science` on 2026-07-16. It uses the same post-CP2K-#5582 build as
`../gxtb_derivative_regtests_post5582_20260716/`.

The twelve successful cases cover a one-dimensional Ar chain and a
two-dimensional Ar layer. For each dimensionality the archive contains native
Bloch sampling, a Gamma-centered native mesh, the commensurate Gamma
supercell, and force/virial finite differences with a full mesh, explicit K290
reduction, and explicit SPGLIB reduction.

All production runs have return code zero and `PROGRAM ENDED AT`. K290 and
SPGLIB reduce the 1D mesh from 2 to 1 point and the 2D mesh from 4 to 2 points.
Both reproduce the corresponding full-grid energy and analytical derivative
components at the printed precision. `partial_pbc_summary.csv` contains the
finite-difference residuals.

The Gamma-supercell energies per primitive cell differ from the native
Gamma-centered Bloch calculations by `+1.53572793806234e-5 Eh` (1D) and
`-4.71936336907675e-5 Eh` (2D), or `+0.0403205` and `-0.123907 kJ mol-1`.
These residuals are small but not machine-precision identities. They remain an
open diagnostic of lower-dimensional image sums, electrostatics, and the
environment-dependent q-vSZP state; this archive therefore does not claim
general 1D/2D production readiness.

The first detached launch used an unavailable bare `mpirun` command and ended
immediately with return code 127. Its metadata and empty run shells are kept in
`runs_failed_missing_mpirun/`. The production relaunch uses the absolute MPI
path and `--bind-to none`, so it is independent of login-shell state and CPU
affinity defaults.

`SHA256SUMS` covers the complete archive. Each production run also contains
its own input, output, return code, and initial/final digest manifests.
