# Expanded derivative validation

This archive contains independent central-finite-difference checks for two
SPGLIB-reduced g-xTB calculations in CP2K.

- The nonsymmorphic ammonia crystal uses a MacDonald `2 x 2 x 2` mesh. SPGLIB
  retains four of eight k points.
- Cubic methane uses a Monkhorst-Pack `6 x 6 x 6` mesh. SPGLIB retains ten of
  216 k points.

Both calculations use `DX 0.0001`, `EPS_SCF 1.0E-10`, one checked Cartesian
force component, and the complete nine-component numerical cell derivative.
The CP2K DEBUG driver performs 23 converged SCF calculations per case.

| Case | Checked force sum / Eh bohr^-1 | Virial sum / Eh | Exit status |
|---|---:|---:|---:|
| NH3, SPGLIB 8 -> 4 | 0 | 4.850e-9 | 0 |
| CH4, SPGLIB 216 -> 10 | 0 | 1.130e-9 | 0 |

`run_identity.sha256`, the per-result binary and input hashes, affinity proofs,
raw CP2K outputs, and exit-status files record the calculation provenance.
