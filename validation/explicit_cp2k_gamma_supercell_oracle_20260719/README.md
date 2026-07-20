# Explicit CP2K Gamma-supercell oracle

This archive qualifies the independent CP2K side of the native-Bloch versus
explicit Born--von Karman comparison for DMC-ICE13 ice XVII on the `2 x 2 x 2`
mesh.  The generated 144-atom input is an explicit Gamma-only supercell: it
has `PERIODIC XYZ`, contains no `KPOINTS` section, and preserves the POSCAR
cell, atom order, and every Cartesian coordinate exactly.

The input gate is recorded in `inputs/XVII/input-verification.json`.  The
energy calculations used the qualified CP2K executable with SHA-256
`b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f`
and the direct `save_tblite` executable with SHA-256
`f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a`.
Every run has a singleton-CPU affinity proof and used one OpenMP/BLAS thread.
The independent verifier in `scripts/compare_gamma_supercell_oracle.py`
reparses normal termination and compares, per primitive cell, all three
routes:

1. symmetry-reduced native `2 x 2 x 2` Bloch k points;
2. this explicit CP2K Gamma supercell;
3. the direct `save_tblite` explicit-BvK CLI calculation.

The verification passes.  Native Bloch sampling and the explicit CP2K Gamma
supercell differ by only `1.1254996934439987e-11` hartree per primitive cell.
The maximum pairwise residual among all three routes is
`4.298783551348606e-9` hartree per primitive cell, far below the conservative
`2e-7`-hartree gate.  Thus this independent Born--von Karman oracle excludes
the native k-space transform and k-point symmetry reduction as the source of
the much larger DMC-ICE13 deviations.

The complete machine-readable comparison is in `verification.json`; raw
inputs, outputs, executable and input hashes, and affinity evidence are kept
under `inputs/` and `results/`.

The reusable `scripts/make_gamma_supercell_input.py` generator reproduces the
archived ice-XVII input and its verification record byte for byte.  It also
provides the independently gated 288-atom ice-VII input under `inputs/VII/`.
Ice VII is the phase with the largest residual in the complete direct-CLI
versus native-CP2K matrix.  `scripts/run_vii_gamma_oracle.sh` stages the
corresponding explicit-Gamma calculation with the same qualified binary,
input-hash gate, singleton CPU affinity, and conservative memory gate.  The
calculation terminated normally and the independent three-route check in
`verification-vii.json` passes: native Bloch sampling and the explicit CP2K
Gamma supercell differ by `1.192184981846367e-7` hartree per primitive cell,
while the maximum pairwise difference among native CP2K, Gamma-supercell CP2K,
and direct `save_tblite` is the same value and remains below the conservative
`2e-7`-hartree gate.  Ice VII therefore independently confirms the conclusion
from ice XVII that neither the native k-space transform nor k-point symmetry
reduction explains the much larger benchmark deviation.
