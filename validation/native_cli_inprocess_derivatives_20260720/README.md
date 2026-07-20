# Native/CLI in-process derivative parity

This gate exercises CP2K's built-in `REFERENCE_CLI` path with the exact
qualified CP2K binary and direct `save_tblite` CLI used for the Part-I
DMC-ICE13 comparison.  Unlike the complete 52-point multi-k energy matrix,
this diagnostic compares energy, Cartesian gradients, and, where available,
the cell virial inside one CP2K execution.

Two independent cases are retained:

- fully periodic cubic H2O at Gamma, including energy, forces, and virial;
- nonperiodic spin-polarized triplet O2, including energy and forces.

Both calculations use `STOP_ON_ERROR T` and an a priori component limit of
`1e-5`.  The verifier additionally enforces a tighter `1e-9 Eh` energy limit,
checks the executable and input identities, and independently confirms the O2
CLI energy from a separately pinned direct run.

The first H2O preflight used a program path longer than CP2K's 80-character
input-field limit and therefore stopped before calculation.  It is retained
under `raw/preflight_program_path_failure` as provenance but is excluded from
the scientific gate.  The successful inputs use a short symbolic path whose
target has the qualified CLI SHA-256.

Run from the repository root with

```bash
python3 validation/native_cli_inprocess_derivatives_20260720/verify_inprocess_derivative_parity.py
```
