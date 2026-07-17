# Symmetry-star atom-phase cache validation

Verdict: **PASS for numerical equivalence to the uncached reference in the
archived scope**.  The candidate reuses atom-phase factors that are invariant
within the repeated symmetry-star transformations.  This archive is a
correctness qualification and does not claim a measured speedup.

## Qualified scope

Seven reference/candidate pairs cover RKS and UKS, K290, time-reversal and
SPGLIB handling, shifted and unshifted meshes, one- and three-dimensional
periodicity, and both ordinary energy and derivative-validation paths.  Every
output must contain `PROGRAM ENDED`.

For every pair, all filtered scientific lines (`ENERGY|` and `DEBUG|`) are
byte-identical.  This includes 50- and 81-line finite-difference force/stress
diagnostics in the corresponding cases.  The only separately parsed printed
force/stress comparison is the shifted Si/SPGLIB case; its maximum residuals
are `1.12727326999999988e-16` hartree/bohr and
`2.66995176391900044e-10` bar.  All final-energy differences are exactly zero
at printed precision.  Complete case-wise results are in `summary.tsv`.

The verifier requires `1e-10` hartree energy agreement, `1e-8`
hartree/bohr force agreement, `1e-3` bar stress agreement, and exact filtered
scientific lines.

## Source and build identity

- uncached reference revision:
  `eb9b157c728f040ed4f0297e994ee9078c41cfc5`
- cache implementation commit:
  `828918efd961d24a6491e483a8d6a68595138e8f`
- reproducible implementation patch:
  `source/0001-Cache-atom-phases-for-g-xTB-symmetry-stars.patch.gz`
- tested candidate executable SHA-256:
  `22a5339f49316e286b3bd67a782ad1be9d2a31ea5e4040764867b18bd5053abe`

The candidate executable was built from the exact cache worktree before its
change was committed, so its printed revision is its base
`09ac4390801888c7c16b780ca68da24622665728`; the archived patch freezes the
tested cache change.

## Reproduction and integrity

The qualification used one MPI rank and one OpenMP thread on
`MacBook-Pro-von-Thomas-3.local`.  Run from this directory:

```text
python3 verify.py
shasum -a 256 -c SHA256SUMS
```

The verifier must reproduce `summary.tsv` and exit zero.  `SHA256SUMS`
authenticates every archived file except the manifest itself.
