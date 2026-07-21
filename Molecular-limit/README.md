# H2O molecular/periodic-limit check

This directory contains the CP2K data behind the Part-I molecular-limit test.
The same water geometry, charge, spin, model parameters, and SCC threshold are
used in the zero-dimensional route and in progressively enlarged
three-dimensional cubic cells.

## Contents

- `results.csv`: the complete archived 8--200 A CP2K energy/force sequence;
- `results_energy_force_stress_8_100.csv`: the current qualified 8--100 A
  three-panel data used by Fig. 2, including analytical/numerical virial
  agreement and the exact executable hash;
- `figures/h2o_molecular_periodic_limit.tex`: current energy/force/stress
  three-panel figure source, with all panels extended to 100 A;
- `traceless_qq_fix_20260721/`: the accepted pure CP2K-native 0D and 8--100 A
  energy/force outputs, the corresponding stress and numerical-virial
  evidence, executable and library hashes, CPU-affinity records, build/test
  logs, and signed source patch;
- `current_build_20260721/`: superseded 8--50 A evidence retained unchanged for
  historical traceability;
- `traceless_qq_fix_20260721/scripts/analyze_native_forces.py`: strict
  same-build/termination verifier and current table regenerator;
- `scripts/analyze_three_panel_series.py`: superseded 8--50 A table generator
  retained for historical traceability;
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
the zero-dimensional CP2K result and is reported in the native CP2K unit
`Eh/a0`; stress is reported in GPa. After correcting the zero-dimensional
quadrupole--quadrupole contraction, energy, force, and stress converge
monotonically through 100 A and the former force plateau is absent. The
zero-dimensional and 100 A CP2K energies are `-76.437385109217445 Eh` and
`-76.437385305373553 Eh`, respectively. Their difference is
`-1.96156108e-7 Eh` (`-5.15007791e-4 kJ/mol`); at 100 A the largest force
component difference is `2.6020e-7 Eh/a0` and the largest stress component is
`6.8918e-5 GPa`. The maximum analytical-minus-numerical virial difference is
`1.0858e-8 Eh`.

All accepted displayed cases terminate normally and carry CP2K executable
SHA-256 `a606cb0ff838dc1a5f967238154d5c0892da5f5b63488d2c10959ec12d6e4d7c`.
The direct-CLI component-ablation diagnostic is intentionally not part of the
manuscript or Supporting Information. Source and executable provenance is
recorded in [`../PART_I_PROVENANCE.md`](../PART_I_PROVENANCE.md), and curated
artifact hashes are listed in
[`../validation/paper_artifact_sha256.json`](../validation/paper_artifact_sha256.json).
