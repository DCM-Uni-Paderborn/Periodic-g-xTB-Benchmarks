# Periodic g-xTB Part-I implementation audit

`summary.json` is the machine-readable result of running
`tools/verify_part_i_implementation.py` from the root of the direct
`save_tblite` reproduction package on 19 July 2026.

The audit reruns eleven independent archived gates covering absolute and
Ih-referenced CLI/native energies, numerical accuracy, the periodic-response
correction, energy/force/stress derivatives, native-k/BvK grid identity,
provider and model revisions, final-source retention sentinels, the
Wigner--Seitz branch diagnosis, and every portable SHA-256 manifest.  All
eleven gates pass.  This report does not replace the fresh same-host all-phase
`2 x 2 x 2` repetition; that stricter provenance gate is archived separately
once its queued production calculations complete.

Reproduce the report with:

```bash
python3 tools/verify_part_i_implementation.py \
  --output-json validation/part_i_implementation_audit_20260719/summary.json
```
