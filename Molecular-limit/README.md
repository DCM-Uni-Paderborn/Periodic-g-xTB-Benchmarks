# H2O molecular/periodic-limit check

This directory contains the CP2K data behind the Part-I molecular-limit test.
The same water geometry, charge, spin, model parameters, and SCC threshold are
used in the zero-dimensional route and in progressively enlarged
three-dimensional cubic cells.

## Contents

- `results.csv`: the complete 8--200 A CP2K sequence used by the manuscript
  table and convergence figure;
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
40 A and then approaches a small positive offset; the force difference reaches
a plateau near `4.86e-4 eV/A`. Tightening the numerical thresholds and
translating the molecule do not remove this behavior.

Only the CP2K quantities reported in Part I are retained here. Source and
executable provenance is recorded in [`../PART_I_PROVENANCE.md`](../PART_I_PROVENANCE.md),
and the curated result-file hash is listed in
[`../validation/paper_artifact_sha256.json`](../validation/paper_artifact_sha256.json).
