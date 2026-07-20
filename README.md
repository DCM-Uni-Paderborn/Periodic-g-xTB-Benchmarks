# Periodic g-xTB in CP2K — Part I data

This `main` branch is the public data and validation companion to
**Periodic g-xTB in CP2K. I. Analytic Derivatives, Brillouin-Zone-Coupled
Nonlocal Exchange, and Symmetry Reduction**.

It contains only material used directly in Part I or required to reproduce a
Part-I statement. Numerical-method performance material for Part II is kept on
the separate `part-II` branch. Historical files removed from the present tree
remain available through the Git history.

## Repository map

- [`DMC-ICE13/`](DMC-ICE13/): primitive structures, native-Bloch inputs,
  DMC references, the frozen-build uniform and phase-wise g-xTB results, plots,
  and the compact independent recalculation package.
- [`LC10/`](LC10/): the ten-solid equation-of-state benchmark used in the
  manuscript, including final per-solid values, first-passing k-point
  intervals, aggregate comparisons, and plot data.
- [`Molecular-limit/`](Molecular-limit/): the molecular/periodic large-cell
  inputs and the complete CP2K energy/force sequence used in Part I.
- [`validation/`](validation/): energy, force, stress, symmetry, K290/SPGLIB,
  lower-dimensional PBC, primitive-cell/supercell, and source-level regression
  evidence for the Part-I reference implementation.
- [`PART_I_PROVENANCE.md`](PART_I_PROVENANCE.md): source and executable
  identities for the frozen qualification build.
- [`GFN_BASELINE_SOURCE.md`](GFN_BASELINE_SOURCE.md): provenance of the compact
  GFN1-xTB/GFN2-xTB comparison values imported from the separate periodic-GFN
  project.

## Current headline values

### DMC-ICE13

The uniform same-build MAE sequence is:

| Mesh | ME | MAE | RMSE | MaxAE |
|---|---:|---:|---:|---:|
| Gamma | -155.6376 | 163.8345 | 218.1678 | 496.8631 |
| 2x2x2 | -86.0399 | 88.6814 | 125.5465 | 304.9652 |
| 3x3x3 | -34.0485 | 34.0485 | 57.0019 | 158.8881 |
| 4x4x4 | -11.3631 | 11.3655 | 21.9579 | 68.9453 |
| 5x5x5 | -4.2818 | 4.3464 | 8.2379 | 26.8510 |

All entries are in kJ mol-1 per H2O over the twelve non-reference phases and
use ice Ih on the same mesh. The current phase-wise set through at most 8x8x8
has ME -1.6775, MAE 1.7478, RMSE 2.1950, and MaxAE 5.4315 kJ mol-1 per H2O.
It is not the final adaptive result: eight of twelve phases currently satisfy
the declared one-step condition
`|R(N)-R(N-1)| <= 0.10 kJ mol-1 per H2O`.

Machine-readable values are in
[`DMC-ICE13/data/`](DMC-ICE13/data/). The recalculation package for an
independent `save_tblite` evaluation is in
[`DMC-ICE13/reproduction/seidler_dmc13_recalculation/`](DMC-ICE13/reproduction/seidler_dmc13_recalculation/).
The aggregate implementation gate is
[`validation/implementation_audit_20260720/`](validation/implementation_audit_20260720/).

### LC10

All ten solids pass both per-solid criteria,
`|Delta a0| <= 0.025 A` and
`|Delta Ecoh| <= 0.25 kJ mol-1 atom-1`, by at most 9x9x9. The final g-xTB
MAEs are 0.1434 A for lattice constants and 0.2947 eV atom-1 for cohesive
energies. The exact per-solid endpoints and signed changes are in
[`LC10/data/lc10_gxtb_final.csv`](LC10/data/lc10_gxtb_final.csv).

## Data policy

Only compact GFN1-xTB/GFN2-xTB values required for explicit manuscript
comparisons are stored here. Their complete inputs and raw data remain in
[`DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks`](https://github.com/DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks).

Generated working directories are not tracked. Named validation archives are
tracked when they substantiate a manuscript table, figure, or implementation
claim. Integrity values are kept in this repository rather than in the
manuscript or Supporting Information.
