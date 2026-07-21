# Qualified traceless-QQ molecular-limit series

This directory archives the qualified CP2K-native molecular-limit rerun after
the zero-dimensional quadrupole--quadrupole contraction was made consistent
with the traceless atomic-quadrupole convention used by the periodic Ewald
route.

The requested H2O geometry is given in bohr in every archived input. The
series contains the molecular zero-dimensional reference and cubic periodic
cells of 8, 10, 12, 15, 20, 30, 40, 50, 60, 80, and 100 A.

- `cp2k_native/raw/` contains the CP2K-native energy and force inputs, outputs,
  hashes, exit status, qualification status, and pre-exec affinity evidence.
- `stress/raw/` contains the corresponding CP2K stress and numerical-virial
  evidence; direct-CLI diagnostic files are deliberately excluded.
- `build_provenance/` contains build, link, library, and regression-test logs.
- `source/` contains the compressed signed source patch at save_tblite commit
  `fad7fe4b188f99794d7c047d5b710667c3a2ce84`.
- `scripts/analyze_native_forces.py` regenerates the qualified table from the
  raw CP2K outputs.

All accepted calculations terminated normally on exact singleton CPU 42 with
OMP/BLAS thread counts fixed to one. They carry CP2K executable SHA-256
`a606cb0ff838dc1a5f967238154d5c0892da5f5b63488d2c10959ec12d6e4d7c`.
The force comparison is reported in the native CP2K unit `Eh/a0`; stress is
reported in GPa. The direct-CLI debugging comparison is not part of the Part-I
manuscript or Supporting Information.

The prelaunch resource records retain `MemAvailable` and the PID/UID/RSS/CPU/
command-name fields for every live process. Command arguments were omitted
from the public archive because they are not needed for the memory-admission
proof.
