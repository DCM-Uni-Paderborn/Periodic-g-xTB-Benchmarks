# KGROUP_PARTIAL_ROOT independent oracle harness

This temporary, uncommitted harness compares the selectable CP2K
`KGROUP_PARTIAL_ROOT` bridge with the unchanged `DENSE` bridge at identical MPI
rank counts.  It freezes the already-qualified KGROUP_OWNER inputs and covers:

- CH4 full `2x2x2`, P=1/2/4;
- CH4 K290 and SPGLIB reductions, P=1/2;
- H2 and UKS O2 time-reversal `3x1x1`, P=1/2/4, including P>nred;
- Ar2 1D and Ar4 2D periodicity, P=1/2;
- shifted-grid Si, P=1/2/4.

Every partial-root output must contain the in-process complete-mesh forward
oracle (`dE`, `dVsh`, folded Fock response) and reverse oracle (overlap
adjoint, direct force, direct stress) with each residual at most `1e-10`.
Final energy, atomic forces, and analytical stress are parsed independently
from DENSE and partial-root outputs.  Return code zero and exactly one
`PROGRAM ENDED` marker are mandatory.

The internal forward/reverse oracle ceiling is `1e-10`.  The independently
printed DENSE-vs-partial-root observables use the SI acceptance limits:
`1e-9` Ha for the total energy, `1e-7` Ha/bohr for atomic forces, and
`1e-5` GPa = `0.1` bar for analytical stress.  Actual residuals are retained
at full parsed precision even when they are far below these limits.  A true
one-point Gamma mesh is intentionally not claimed as a Mode-6 reverse test:
CP2K dispatches that mathematically trivial derivative through the dense
one-point fallback.  The present matrix begins with nontrivial meshes and
therefore exercises the partial-root reverse bridge itself.

Threading is forced to one for OpenMP, OpenBLAS, MKL, BLIS, GotoBLAS, and
Accelerate.  On terok the intended invocation pins every sequential MPI case
to the reserved CPU set 128--143; no cases are run concurrently.

Example:

```sh
./run_matrix.py --cp2k /path/to/cp2k.psmp \
  --rank-prefix 'taskset -c 128-143'
./verify_matrix.py
./freeze_manifest.py
sha256sum -c SHA256SUMS
```

The runner refuses to overwrite a nonempty run directory and the verifier
fails on missing pairs, missing markers, metadata/hash drift, non-finite
values, or a numerical gate violation.
