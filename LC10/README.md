# LC10 periodic g-xTB benchmark

This directory contains the fixed ten-solid Part-I subset:

```text
C, Si, SiC, BN, BP, AlN, AlP, MgS, LiF, LiCl
```

MgO and LiH are not members of the reported benchmark and no longer appear in
the current data tree.

## Current paper snapshot

The current Part-I paper reports one qualified, common native-Bloch `7x7x7`
mesh for all ten solids.  It is deliberately described as a common-mesh
snapshot rather than as a k-point-converged adaptive result: every `6x6x6` to
`7x7x7` lattice-constant increment is at most `0.005 A`, whereas none of the
cohesive-energy increments is at most `0.05 kJ mol-1 atom-1`.

## Paper data

- `data/lc10_gxtb_uniform_k777_snapshot.csv`: current per-solid `7x7x7`
  lattice constants, cohesive energies, signed errors, and adjacent-mesh
  changes at manuscript precision;
- `data/lc10_method_comparison.csv`: signed mean errors and mean absolute
  errors for current g-xTB and the literature GFN/post-HF context over the
  exact common ten-solid set;
- `figures/lc10_method_comparison.*`: current method-context figure generated
  by `scripts/generate_lc10_method_comparison.py`;
- `data/reference_goldzak2022.csv`: zero-point-corrected experimental and
  high-level reference values.

For the current uniform `7x7x7` g-xTB snapshot, ME/MAE are
`-0.132181178/0.134252148 A` for lattice constants and
`+0.144612412/0.190758163 eV atom-1` for cohesive energies.  The literature
rows use different numerical protocols (`4x4x4` for the CP2K/tblite GFN rows
and the protocols of Goldzak et al. for HF/MP2 variants); they provide method
context and are not same-mesh implementation comparisons.

The plotting script regenerates the repository figure from the versioned
table. Complete GFN1-xTB/GFN2-xTB raw data remain in the separate
`Periodic-GFN2-Benchmarks` repository.
