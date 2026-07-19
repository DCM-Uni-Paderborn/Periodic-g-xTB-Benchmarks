# Periodic `save_tblite` source qualification

This archive qualifies the source-level periodic paths used by the CP2K-native
g-xTB integration.  It complements the end-to-end CLI/native comparisons: the
latter prove that both front ends return the same energy, while these unit tests
exercise the periodic H0, Wigner--Seitz, and Brillouin-zone-coupled exchange
kernels directly.

## Qualified revisions

- current integration source: `15915c9435644eb257178ca8f8bf7220c38b1a84`
- historical `pbc` baseline: `c932120d2580811901de6a1fe3f89b943c251766`
- compiler: GNU Fortran 16.1.0, Release build, one OpenMP/BLAS thread

The historical baseline did not build its unit-test executable as committed.
Two test-registration repairs were applied in its detached diagnostic worktree:
the missing comma in `test/unit/test_gxtb.f90` and the omitted
`wignerseitz` source entry in `test/unit/CMakeLists.txt`.  No library source was
changed.  `results.json` records these two and only these two modified files.

## Results

| Test group | Result | Numerically relevant bound |
|---|---:|---:|
| periodic H0 anisotropy and gradient | 3/3 pass | test tolerances satisfied |
| Wigner--Seitz images, including partial periodicity | 6/6 pass | test tolerances satisfied |
| exchange, Fock response, gradients, stress, and transforms | 36/36 pass | explicit-BvK residual: `4.6134320322299693e-13` Ha |
| ACP | 22/22 pass | image, k-space, gradient, and stress paths included |
| charge and multipole Coulomb | 109/109 pass | periodic gradients and stress included |
| dispersion and repulsion | 34/34 pass | periodic energy/gradient/stress cases included |
| complete g-xTB model suite | 40/40 pass | periodic ACP/H0/k-mesh cases included |
| full Hamiltonian suite | 74/75 pass | sole failure described below |

The mixed-radix exchange transform differs from its dense oracle by at most
`9.7144514654701197e-17` Ha.  Thus neither the regular-mesh acceleration nor the
whole-mesh exchange contraction introduces a material energy difference.

The only full-suite failure is the nonperiodic CeCl3 numerical-gradient test.
Its maximum analytic/finite-difference difference is
`2.87804443471762e-10` Ha/bohr, narrowly above the test threshold
`2.220446049250313e-10` Ha/bohr.  All twelve difference components are
bit-for-bit identical in the historical `pbc` baseline and the current source
under the same compiler and dependencies.  It is therefore an inherited,
compiler-sensitive finite-difference threshold case, not a regression of the
periodic implementation.

## Reproduction and verification

With both build worktrees present, rerun the archive with:

```sh
python3 run_source_tests.py \
  --current-root /path/to/current-save_tblite \
  --pbc-root /path/to/pbc-baseline
python3 verify_source_tests.py
shasum -a 256 -c SHA256SUMS
```

`results.json` contains executable hashes, commands, return codes, pass/fail
counts, residuals, and the complete CeCl3 difference vectors.  The raw command
output is retained in the corresponding `.log` files.
