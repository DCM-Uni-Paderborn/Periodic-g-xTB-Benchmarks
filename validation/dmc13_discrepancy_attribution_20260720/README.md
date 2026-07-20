# DMC-ICE13 discrepancy attribution

This gate independently reconstructs the numerical hierarchy needed to
separate a CP2K-interface error from a provider-source or benchmark difference.
It compares four quantities without treating sparse-mesh DMC statistics as
converged accuracy estimates:

- current pbc-derived direct CLI versus CP2K-native Bloch calculations;
- the upstream author `pbc` snapshot versus the later pbc-derived integration
  provider;
- `mstore-inorganic` versus the current pbc-derived provider at identical
  structures, meshes, and CLI accuracy;
- the still-provisional dense-mesh adaptive DMC-ICE13 statistic.

Run from the repository root with

```bash
python3 validation/dmc13_discrepancy_attribution_20260720/verify_discrepancy_attribution.py
```

The resulting `verification.json` deliberately leaves the provenance of the
previously quoted lower author result open.  Its purpose is narrower and
testable: to show which numerical scales can and cannot originate in the
CP2K-native interface.
