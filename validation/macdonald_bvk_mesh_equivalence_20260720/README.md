# MacDonald mesh / Gamma-supercell BvK equivalence

This gate proves that the native CP2K meshes used for the Part-I DMC-ICE13
production data are exactly the reciprocal grids obtained by folding a
Gamma-only `N x N x N` Born--von Karman supercell into the primitive-cell
Brillouin zone.  It is not an inference from close energies.

For the production path, `GAMMA_CENTERED` is off and CP2K's
`full_grid_gen` constructs the one-dimensional reduced coordinate

```text
k_i = (2 i - N - 1)/(2 N) + s,  i = 1,...,N.
```

The frozen inputs use `s = 0` for odd `N` and
`s = (N-1)/(2N)` for even `N`.  Hence

- odd `N`: `k_i = (i-(N+1)/2)/N`, whose values modulo a reciprocal lattice
  vector are `{0,1/N,...,(N-1)/N}`;
- even `N`: `k_i = (i-1)/N` directly.

The Cartesian product in three dimensions is therefore precisely the
Gamma-supercell BvK grid for every mesh used here.  Space-group and
time-reversal reduction only select representatives and weights from this
same full grid; the g-xTB coupling reconstructs the full mesh internally.

`verify_macdonald_bvk_mesh.py` performs the mathematical proof with exact
rational arithmetic for `N=1,...,9` and parses every archived production
input to verify its mesh, shift, symmetry setting, and generated point set.
For non-terminating decimal shifts such as `5/12`, it additionally verifies
that the written decimal and the exact fraction round to the same IEEE-754
binary64 value used by CP2K; the rational text-level residual is retained in
the machine-readable output.  It also checks the source excerpt against the
qualified CP2K source identity recorded in this repository.

Run with:

```bash
python3 validation/macdonald_bvk_mesh_equivalence_20260720/verify_macdonald_bvk_mesh.py
```
