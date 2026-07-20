# Direct-CLI accuracy and CP2K parity control

This validation resolves the convergence-setting ambiguity in the DMC-ICE13
direct-CLI archive.  The qualified CP2K inputs use `TBLITE/ACCURACY 0.1`.
Fresh, exactly pinned direct `save_tblite` calculations were therefore run for
ice Ih and the particularly sensitive ice VII cell at `1 x 1 x 1` with both
`--acc 0.1` and the ten-times tighter `--acc 0.01`.

The text output independently proves the effective setting:

| `--acc` | printed energy threshold | printed density threshold | iterations |
|---:|---:|---:|---:|
| 0.1 | `1e-7 Eh` | `2e-6 e` | 8 |
| 0.01 | `1e-8 Eh` | `2e-7 e` | 10 |

Tightening the direct calculation changes the Ih energy by only
`-1.9065e-10 Eh` and the ice VII energy by `-1.8986e-11 Eh`.  At the identical
`0.1` setting, CP2K-native differs from the direct CLI by `1.2696e-8 Eh` for
Ih and `2.0150e-8 Eh` for ice VII.  Both are far below the conservative
`2e-7 Eh` parity gate.

The common basis-cutoff formula saturates at its maximum integral threshold
for both `0.1` and `0.01`; only the SCC convergence criteria become tighter in
this range.  The current upstream implementation also contains a reversed
clamping-bound order, but correcting that order would not distinguish these
two settings because both still lie in the saturated high-accuracy range.
This shared cutoff issue therefore does not explain any CP2K/direct-CLI or
`mstore-inorganic`/`pbc` energy difference.

The previous direct-CLI archive was described as `0.01`, but its complete text
outputs show that all `1^3`, `2^3`, and `3^3` calculations and ten of twelve
completed `4^3` calculations actually used `0.1`.  The final parity matrix is
now gated by the printed thresholds.  The two tighter `4^3` outputs are
preserved here before their same-setting replacements are installed.

Reproduce the controlled comparison with:

```text
python3 verify_accuracy_equivalence.py --output verification.reproduced.json
```
