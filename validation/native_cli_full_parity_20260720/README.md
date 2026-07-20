# Complete direct-CLI/native parity gate

This gate verifies the complete same-provider energy matrix used to decide
whether the CP2K-native periodic g-xTB implementation reproduces direct
`save_tblite` calculations.

It requires all thirteen structures at every mesh from `1 x 1 x 1` through
`4 x 4 x 4` (52 unique points), effective direct-CLI accuracy `0.1`, passing
input qualification, and three SHA-256 provenance fields per point.  It then
recomputes every primitive-cell absolute difference and, independently from
the precomputed relative-energy tables, references both routes to same-mesh
ice Ih per water molecule.

The exact acceptance limits are

- `2e-7 Eh` per primitive cell for every absolute energy; and
- `5e-5 kJ mol-1 H2O-1` for every same-mesh Ih-referenced relative energy.

Run from the repository root with

```bash
python3 validation/native_cli_full_parity_20260720/verify_native_cli_full_parity.py
```

The generated `verification.json` and `verification.stdout` contain all
individual checks and per-mesh maxima.  This is an implementation-parity
gate, not a claim that a sparse DMC-ICE13 mesh is physically converged.
