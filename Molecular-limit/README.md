# H2O molecular/periodic-limit check

This directory contains the CP2K data behind the Part-I molecular-limit test.
The same water geometry, charge, spin, model parameters, and SCC threshold are
used in the zero-dimensional route and in progressively enlarged
three-dimensional cubic cells.

## Contents

- `results.csv`: the complete archived 8--200 A CP2K energy/force sequence;
- `results_energy_force_stress_8_50.csv`: the current qualified 8--50 A
  three-panel data used by Fig. 2, including analytical/numerical virial
  agreement and the exact executable hash;
- `figures/h2o_molecular_periodic_limit.tex`: current energy/force/stress
  three-panel figure source, with all panels restricted to 50 A;
- `current_build_20260721/raw_energy_force/`: accepted CP2K inputs, outputs,
  executable hashes, termination status, and pre-exec affinity records for the
  0D and 8--50 A energy/force sequence;
- `current_build_20260721/raw_stress/`: corresponding accepted stress and
  numerical-virial evidence;
- `scripts/analyze_three_panel_series.py`: strict same-build/termination
  verifier and table regenerator;
- `raw/*/input.inp`: the corresponding CP2K inputs;
- `H2O_gxtb_molecular_limit.inp`: the canonical production template;
- `H2O_gxtb_molecular_limit_tight.inp` and
  `H2O_gxtb_molecular_limit_shifted.inp`: the tighter-threshold and translated
  controls cited in the manuscript.

The reported energy difference is

```text
Delta E(L) = E_CP2K,3D(L) - E_CP2K,0D .
```

The force quantity is the largest Cartesian component difference relative to
the zero-dimensional CP2K result. The energy crosses the molecular value near
40 A; over the displayed range the force difference reaches about
`4.9e-4 eV/A`, while the largest absolute stress component decreases from
`530.84 MPa` at 8 A to `3.80 MPa` at 50 A.  The maximum analytical-minus-
numerical virial difference is `1.184e-8 Eh`.

All accepted displayed cases terminate normally and carry CP2K executable
SHA-256 `b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f`.
The direct-CLI component-ablation diagnostic is intentionally not part of the
manuscript or Supporting Information. Source and executable provenance is
recorded in [`../PART_I_PROVENANCE.md`](../PART_I_PROVENANCE.md), and curated
artifact hashes are listed in
[`../validation/paper_artifact_sha256.json`](../validation/paper_artifact_sha256.json).
