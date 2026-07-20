# Second-order minimum-image attribution

This reciprocal one-patch test resolves the small source-state residual that
remains after correcting the historical Wigner--Seitz self-image index.
It compares ice VII and, independently, ice XVII with same-source ice Ih in
explicit `2 x 2 x 2` BvK supercells.  The VII and Ih cells contain 96 water
molecules; the smaller XVII cell contains 48.

The three independently converged source states are:

- author `pbc` with the corrected Wigner--Seitz indexing;
- exact historical `mstore-inorganic` plus only that Wigner--Seitz correction;
- author `pbc` with only the later minimum-image second-order change reverted.

All runs use g-xTB, `--acc 0.1`, 300 allowed SCC iterations, no restart, the
same structures, and matching Release compiler/dependency settings.  The
staged inverse patch is archived as `source/revert-083f220.patch`; its changed
file set and resulting Git tree are recorded separately.

| Source state | VII minus Ih / kJ mol-1 H2O-1 |
|---|---:|
| corrected `pbc` | -300.067320330933 |
| WSC-corrected `mstore-inorganic` | -305.888458118205 |
| `pbc` without the minimum-image second-order variant | -305.888469079012 |

Thus the corrected-`pbc`/corrected-`mstore-inorganic` residual is
`5.821137787272 kJ mol-1 H2O-1`.  Reverting only the minimum-image
second-order change leaves `-0.000010960807 kJ mol-1 H2O-1`, explaining
`99.9998117%` of that residual.  The remaining value is below the SCC-level
numerical resolution of this diagnostic.

The independent ice-XVII cross-check gives:

| Source state | XVII minus Ih / kJ mol-1 H2O-1 |
|---|---:|
| corrected `pbc` | -51.559121449977 |
| WSC-corrected `mstore-inorganic` | -51.723327684319 |
| `pbc` without the minimum-image second-order variant | -51.723328473957 |

Here the original residual is `0.164206234342 kJ mol-1 H2O-1`, whereas the
no-MIC `pbc` result differs from WSC-corrected `mstore-inorganic` by only
`0.000000789638 kJ mol-1 H2O-1`, independently explaining `99.9995191%`.

Together with the independent Wigner--Seitz reciprocal-patch test, this
classifies the complete historical sparse-mesh branch separation: the large
part comes from misassigned periodic exchange self-images, and the small
post-WSC part comes from the subsequently changed second-order periodic
Coulomb kernel.  Neither effect is a CP2K-native/CLI interface discrepancy.

Run `python3 evaluate_second_order_mic_attribution.py` to reproduce all
numerical, input, source, build-option, and raw-file checks.
