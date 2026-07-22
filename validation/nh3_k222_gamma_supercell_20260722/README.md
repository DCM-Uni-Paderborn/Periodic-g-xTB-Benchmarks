# NH3 k-point reduction and Gamma-supercell comparison

This archive supports the ammonia-crystal entries in Table I and Table SII of
Part I. The nonsymmorphic P2_13 primitive cell uses a Gamma-compatible
MacDonald 2 x 2 x 2 mesh. K290 and SPGLIB both retain four of the eight k
points identified from twelve space-group operations.

| Route | Total energy / Eh | Energy per primitive cell / Eh |
|---|---:|---:|
| Native full mesh | -226.435296880878070 | -226.435296880878070 |
| Native K290 | -226.435296880878042 | -226.435296880878042 |
| Native SPGLIB | -226.435296880878042 | -226.435296880878042 |
| Eightfold Gamma supercell | -1811.482375046768539 | -226.435296880846067 |

The Gamma-supercell-minus-native-full difference is
`+3.2003e-11 Eh` per primitive cell. K290 and SPGLIB differ from the explicit
full mesh by `+2.8e-14 Eh`.

The input and output pairs retain the complete CP2K records used for these
comparisons. `SHA256SUMS` provides file-level integrity checks.
