# Direct CLI/native identity at 4 x 4 x 4

This archive is the final same-host energy-identity check between the direct
periodic `save_tblite` executable and the CP2K-native periodic g-xTB path at
`4 x 4 x 4`. It covers ice Ih plus phases II, VII, XI, and XIV. The set was
chosen to include the reference phase, the largest DMC-ICE13 outlier, and
different primitive-cell sizes and symmetry reductions.

Every calculation is tied to the exact executable and input SHA-256 recorded
beside its output. The pre-exec records prove disjoint singleton CPU affinity,
both executables terminated normally, and the direct supercell energy is
divided by exactly `4^3` before comparison with the native primitive-cell
energy. The Ih-referenced comparison also accounts for the number of water
molecules in each primitive cell.

The largest absolute native-minus-direct difference is
`1.1152530987601494e-7` hartree per primitive cell. The largest Ih-referenced
difference is `1.8722655024642987e-5` kJ mol-1 per water molecule. Both pass
the frozen gates of `2e-7` hartree and `5e-5` kJ mol-1 per water molecule.
Together with the complete Gamma, `2 x 2 x 2`, and `3 x 3 x 3` comparisons,
this excludes the CP2K/CLI interface, units, primitive-cell normalization, and
native-k/BvK mapping as explanations for the much larger DMC benchmark error.

This is an implementation-validation record. It is intentionally retained in
the repository rather than presented as a Part-I benchmark result in the
manuscript.

## Contents

- `inputs/native/`: frozen CP2K inputs; the matching direct POSCAR files are
  the canonical `structures/k444` files in the reproduction package;
- `raw/cli/`: direct executable outputs, results, exit states, hashes, and
  affinity proofs;
- `raw/native/`: CP2K outputs, exit states, hashes, and affinity proofs;
- `report.json`: the machine-readable five-phase comparison;
- `verify_archive.py`: fail-closed recomputation of all energy and provenance
  gates;
- `SHA256SUMS`: complete integrity manifest for this directory.

Reproduce the archived verdict with:

```bash
python3 verify_archive.py
```
