# Periodic g-xTB Part-I implementation audit

`summary.json` is the machine-readable result of running
`tools/verify_part_i_implementation.py` from the root of the direct
`save_tblite` reproduction package on 19 July 2026.

The audit reruns twenty independent archived gates covering absolute and
Ih-referenced CLI/native energies, numerical accuracy, the periodic-response
correction, energy/force/stress derivatives, native-k/BvK grid identity,
provider and model revisions, final-source retention sentinels, the
Wigner--Seitz branch diagnosis, final-build low-k and partial-PBC derivatives,
exchange/ACP component ablations, the periodic-H0 source attribution and
equivalent-image invariant, direct periodic source tests for H0,
Wigner--Seitz, exchange, force, and stress paths, the complete all-phase
author-`pbc`/current-CLI/CP2K-native `3 x 3 x 3` closure, the portable
author-facing DMC-ICE13 recalculation package, the fail-closed adaptive
production-controller dry run, and every portable SHA-256 manifest.  All
twenty gates pass.  The same-host, all-phase CLI/native repetitions at Gamma,
`2 x 2 x 2`, and `3 x 3 x 3` are retained by their dedicated provenance gates;
the independent `4 x 4 x 4` sentinel comparison is retained separately.

Reproduce the report with:

```bash
python3 tools/verify_part_i_implementation.py \
  --output-json validation/part_i_implementation_audit_20260719/summary.json
```
