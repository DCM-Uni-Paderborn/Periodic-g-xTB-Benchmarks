# Explicit CP2K Gamma-supercell oracle

This archive qualifies the independent CP2K side of the native-Bloch versus
explicit Born--von Karman comparison for DMC-ICE13 ice XVII on the `2 x 2 x 2`
mesh.  The generated 144-atom input is an explicit Gamma-only supercell: it
has `PERIODIC XYZ`, contains no `KPOINTS` section, and preserves the POSCAR
cell, atom order, and every Cartesian coordinate exactly.

The input gate is already complete in `inputs/XVII/input-verification.json`.
The energy run must use the qualified CP2K executable with SHA-256
`b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f`,
singleton CPU affinity, and one OpenMP/BLAS thread.  Its result is accepted
only after `compare_gamma_supercell_oracle.py` independently reparses normal
termination and compares, per primitive cell, all three routes:

1. symmetry-reduced native `2 x 2 x 2` Bloch k points;
2. this explicit CP2K Gamma supercell;
3. the direct `save_tblite` explicit-BvK CLI calculation.

No energy result is recorded before all execution and input hashes are
available.
