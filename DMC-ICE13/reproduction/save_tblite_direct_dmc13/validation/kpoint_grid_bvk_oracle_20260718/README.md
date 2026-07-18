# Native k-point/BvK grid oracle

This gate verifies the regular meshes used for the native CP2K calculations
against the reciprocal grid represented by an explicit Born--von-Karman
supercell.  It is independent of the energy comparison: CP2K prints every
full-grid coordinate and weight, and `verify_grid.py` compares the canonical
point set with

```text
{(i/N, j/N, k/N) mod 1 | i,j,k = 0,...,N-1}
```

for `N = 2, 3, 4`.

The even MacDonald grids use the explicit shift `(N-1)/(2N)`.  Applied to
the underlying Monkhorst--Pack coordinate
`(2*i-N-1)/(2*N)`, this gives `(i-1)/N` exactly.  Odd grids use zero shift;
their negative representatives are equivalent to the same BvK points modulo
a reciprocal lattice vector.  Thus the alternating input notation does not
alternate between different physical meshes.

All three small Si g-xTB jobs ended normally.  They use the same k-point
generator and current CP2K build as the DMC-ICE13 calculations.  The separate
DMC absolute-energy oracle then tests the complete Hamiltonian and
normalization on the actual ice structures.

Run the machine-readable gate with

```bash
python3 verify_grid.py
```

