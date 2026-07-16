# Final verification report

## Verdict

**PASS** for the reviewed correctness-first Root-A bridge at the exact source
and build identities recorded in this bundle.

The verdict means that the new selectable `KGROUP_PARTIAL_ROOT` implementation
is statically consistent with the provider ABI, fails closed at its reviewed
protocol boundaries, and is numerically equivalent to the unchanged explicit
full-mesh dense reference for the frozen validation matrix.  It does not mean
that the root-serial nonlinear kernel and R-to-k work have been accelerated.

## Reperformed bundle checks

- CP2K source SHA-256: exact expected match
- CP2K patience-patch SHA-256: exact expected match
- Provider source-tar SHA-256: exact expected match
- Provider CTest classification manifest: all entries OK
- Final-smoke remapped manifests: 62 checked, 0 missing, 0 mismatched
- Oracle manifest: 181 entries OK
- Oracle fail-closed verifier: 20/20 DENSE/partial-root pairs PASS
- Oracle tamper self-test: forged output digest rejected
- Provider focused exchange suite: 31/31 Release, 31/31 Debug
- CP2K positive/fault/Gamma status tables: 22/22 rows PASS
- CP2K Release and Debug builds: both complete through `[4159/4159]`
- CP2K `git diff --check`: clean

## Oracle numerical envelope

Across the 20 frozen comparison pairs:

- maximum externally parsed total-energy difference:
  `1.4210854715202004e-14` Ha
- maximum externally parsed force-component difference:
  `4.7081310160000002e-16` Ha/bohr
- maximum externally parsed analytical-stress component difference:
  `3.9999995351536199e-06` bar
- maximum of all recorded in-process dense-oracle residual fields:
  `2.220446e-15`

These are below the frozen external gates (`1e-9` Ha, `1e-7` Ha/bohr,
`0.1` bar) and the internal `1e-10` forward/reverse gate.

## Coverage

The independent matrix exercises full-grid 3D, K290 and SPGLIB symmetry
reduction, shifted 3D, time reversal, P>nred, RKS, UKS, 1D, and 2D at selected
MPI sizes P=1,2,4.  Each pair compares final energy, atomic forces, and
analytical stress and also checks the internal energy, shell-potential, folded
Fock, overlap-adjoint, direct-force, and direct-stress dense oracles.

Explicit one-point Gamma is deliberately classified separately: its Mode-6
forward result agrees with explicit dense, while CP2K routes the derivative
through the established dense one-point fallback.  The small explicit-vs-
implicit Gamma baseline remains unchanged under tighter settings and is not a
partial-root error.
