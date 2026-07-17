# Separable direct-DFT qualification

The new backend is a factorized **separable DFT**, not an FFT.  No FFT
dependency is linked by this save_tblite build.  It applies one direct DFT per
regular-grid pencil and changes the character-product cost from
`O(Nk^2)` to `O(Nk * sum(nmesh))`; the dense phase-table/BLAS backend remains
selectable as the numerical oracle.

## Correctness

- Debug build (`-O0 -fcheck=all -fbacktrace`): exchange suite 31/31 passed.
- Debug build: g-xTB suite 44/44 passed (40 ordinary passes and four expected
  failure diagnostics).
- The focused explicit transform oracle (7x1x1, 3x4x1, and 2x3x4): arbitrary
  k ordering, arbitrary representative ordering, a common shifted-grid twist,
  and exact normalized inverse/forward adjoint identity passed at an absolute
  tolerance of `1e-11`.
- End-to-end dense-versus-separable comparisons cover RKS and UKS energy,
  Fock, shell potential, overlap adjoint, atomic force response, and stress at
  `1e-11`.
- A near-regular mesh that lies within the dense diagnostic tolerance but is
  not machine regular is explicitly rejected by the separable backend.  The
  same applies to non-machine-uniform weights.  The dense route retains its
  prior behavior.

## Timing

The accompanying TSV contains release (`-O3`) transform-only timings on the
local Apple-silicon host, averaged over 20 paired k-to-R/R-to-k calls with 64
complex rows.  At fixed `Nk=216`, the separable backend is expectedly slower
for a purely 1D `216x1x1` mesh (0.72x), but is 4.00x faster for `18x12x1` and
5.75x faster for `6x6x6`.  At `Nk=512`, it reaches 5.74x for `32x16x1` and
10.08x for `8x8x8`.  This is the expected dimensionality dependence of a
separable direct DFT; no `O(Nk log Nk)` claim is made.  The deliberately
ill-shaped 512x1x1 stress case accumulates an absolute dense/separable
difference of `1.3911e-11` (against an unnormalised forward result whose scale
grows with Nk); it is listed transparently in the TSV and is not part of the
focused `1e-11` gate.

## Memory

The dense plan stores two `Nk x Nk` complex phase tables, i.e. exactly
`2*Nk^2` complex numbers in addition to the common kernel and ordering data.
The compact plan stores `Nk + 3*max(nmesh)` complex numbers in the current
rectangular roots/twist representation, plus two `Nk` integer permutations.
Each separable transform allocates two `Nrow x Nk` complex work arrays; these
are `O(Nrow*Nk)` and independent of `Nk^2`.  For the timed 8x8x8 case this is
524,288 complex phase-table elements for dense versus 536 compact-plan
complex elements and 65,536 complex scratch elements for 64 transformed rows.
