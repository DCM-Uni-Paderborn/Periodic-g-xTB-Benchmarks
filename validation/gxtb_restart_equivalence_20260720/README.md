# Native g-xTB restart equivalence

This gate tests two restart modes of the CP2K-native periodic g-xTB path with
the exact qualified production executable:

1. a strict `3 x 3 x 3` same-mesh restart from a converged WFN payload; and
2. a validated BvK transfer from `3 x 3 x 3` to `4 x 4 x 4`, compared with an
   independent cold `4 x 4 x 4` calculation.

The four calculations use a periodic methane cell, the native Monkhorst--Pack
path, the default g-xTB Fock mixer, one BLAS/OpenMP thread, and an independently
recorded singleton CPU affinity. The same-mesh restart reproduces the source
energy exactly and requires one rather than twelve SCF steps. The cross-mesh
transfer differs from the cold `4 x 4 x 4` result by `7e-15` hartree and
requires seven rather than eleven SCF steps.

The restart is deliberately a density/Fock warm start. Mixer history is reset
when the WFN payload is read, whereas static ACP data may remain cached. This
avoids treating a mesh-dependent DIIS/Broyden history as transferable state.

The corresponding permanent CP2K regression inputs are
`CH4_gxtb_kp_restart_3_same.inp`, `CH4_gxtb_kp_restart_4.inp`, and
`CH4_gxtb_kp_cold_4.inp` on the `g-xTB-pbc` branch. Raw inputs, outputs,
affinity proofs, exit states, the restart payload, and their hashes are under
`raw/`.

Run from the repository root with

```bash
python3 validation/gxtb_restart_equivalence_20260720/verify_restart_equivalence.py
```

This test proves restart equivalence for the tested same-mesh and validated
regular-mesh transfer paths. It does not imply that a run configured with
restart writing disabled can be recovered after interruption.
