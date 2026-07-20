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
- reciprocal one-patch Wigner--Seitz builds that causally isolate the dominant
  historical periodic self-image exchange difference;
- a linked second reciprocal source-patch test that attributes the complete
  post-Wigner--Seitz residual to the later minimum-image form of the
  second-order Coulomb term, independently for ice VII and XVII;
- the still-provisional dense-mesh adaptive DMC-ICE13 statistic.

Run from the repository root with

```bash
python3 validation/dmc13_discrepancy_attribution_20260720/verify_discrepancy_attribution.py
```

The resulting `verification.json` closes the tested source-state chain: more
than 95% of the original sparse-mesh `mstore-inorganic`/`pbc` gap is caused by
the Wigner--Seitz self-image-index correction, and the remaining residual is
removed to below 5e-5 kJ mol\(^{-1}\) H2O\(^{-1}\) by reverting only the later
second-order minimum-image change.  The latter attribution is reproduced for
ice VII and the independent ice XVII cross-check.  Exact source and executable
provenance is still required before assigning a separately quoted author
benchmark to either source state.  The gate therefore shows both which
numerical scales cannot originate in the CP2K-native interface and which
source changes produce the historical CLI trend.
