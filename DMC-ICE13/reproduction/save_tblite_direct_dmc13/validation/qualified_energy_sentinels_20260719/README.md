# Qualified-build energy sentinels

This archive determines whether the regular-mesh and restart hardening applied
after the periodic-response correction changes the DMC-ICE13 energy path. It
does not infer equivalence from source inspection: the reference and qualified
executables are run independently on the same frozen inputs.

The gates comprise ice Ih and ice VII on the native `2 x 2 x 2` mesh and ice
VII on the native, SPGLIB-reduced `3 x 3 x 3` mesh. The first two reference
outputs are copied from the final response-corrected campaign; the `3 x 3 x 3`
reference is repeated beside its qualified counterpart. All calculations use
one thread and an archived singleton-CPU affinity proof.

The three qualified-minus-reference energy differences are exactly zero at
the precision printed by CP2K. The maximum absolute difference is therefore
`0.0` hartree, well below the `1e-12` hartree gate. Consequently, the later
hardening is inactive for these DMC energy inputs and does not require a
blanket repetition of already completed response-corrected points.

## Contents

- `inputs/`: the three frozen CP2K inputs;
- `raw/reference-k222-*`: the two response-corrected `2 x 2 x 2` references;
- `raw/baseline-k333-VII`: the repeated response-corrected `3 x 3 x 3` run;
- `raw/qualified-*`: the three runs from the qualified source tree;
- `raw/comparison.json`: machine-readable energies and differences;
- `provenance/`: source and executable identities for both builds;
- `run_qualified_energy_sentinels.sh`: persistent Terok driver;
- `SHA256SUMS`: complete archive integrity manifest.

Recompute the energy, hash, termination, and affinity gates with:

```bash
python3 verify_qualified_energy_sentinels.py
```
