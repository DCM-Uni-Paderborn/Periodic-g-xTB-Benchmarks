# CP2K periodic-response and direct-CLI validation

This package closes two independent questions for the periodic g-xTB
implementation used in Part I:

1. Does native CP2K k-point sampling reproduce a direct `save_tblite` CLI
   Born--von Karman calculation on the same structures?
2. Which source change caused the previously observed reduction of the
   coarse-grid DMC-ICE13 error?

## CLI/native gate

All thirteen absolute Cartesian DMC-ICE13 structures are evaluated on a
Gamma-centred `2 x 2 x 2` grid.  CP2K uses native k points with SPGLIB
reduction.  The CLI uses the mathematically equivalent explicit `2 x 2 x 2`
Born--von Karman supercell; its energy is divided by eight before comparison.

Across the complete set, the largest absolute-energy difference is
`1.0520830e-7` hartree (ice VII), the RMS difference is `3.0359973e-8`
hartree, and the largest relative-energy difference is only
`2.1413856e-5` kJ mol^-1 per water molecule.  Thus the native and CLI energy
paths are numerically equivalent.  The small residual is consistent with the
different matrix representations and summation orders and is treated as a
numerical parity bound rather than bitwise identity.  The independent input
accuracy gate is stored in `../accuracy_sensitivity_20260718`.

The largest residual was challenged further on both sides.  For ice Ih and
VII, tightening the direct CLI energy and density criteria to `1e-10`
hartree and `2e-9` e changes a primitive-cell energy by at most
`9.49e-11` hartree.  Independently, tightening native ice VII from
`ACCURACY 0.1` and `EPS_SCF 1e-9` to `ACCURACY 1e-4` and `EPS_SCF 1e-12`
changes its energy by only `1.14e-13` hartree.  The tight-on-both-sides
native/CLI difference remains `1.0520e-7` hartree.  The residual is therefore
not caused by premature SCC termination in either front end.  Exact outputs,
inputs, affinity proofs, and exit records are retained in
`raw/cli_tight_scc` and `raw/cp2k_native_tight_scc`.

## Response-fix A/B gate

The old qualified CP2K build and the final clean build link the identical
`save_tblite` static library.  On the deliberately unconverged `2 x 2 x 2`
grid, the DMC-ICE13 MAE changes from `90.8922187` to `88.6813751`
kJ mol^-1, an improvement of `2.2108436` kJ mol^-1 or `2.4324%`.

An exact build of the CP2K parent immediately before the response fix
reproduces the older ice-Ih and ice-XVII energies to the printed digit.  This
excludes the exact symmetry-star/transform accelerations developed earlier as
the origin of the shift.  Explicit 300 K smearing gives exactly the same final
Ih and XVII energies as the zero-temperature default, so the changed default
occupation setting is not responsible either.

Removing ACP from both Ih and XVII retains `5.2328450` of the full
`5.8196408` kJ mol^-1 relative-energy shift (`89.917%`).  The dominant cause is
therefore the corrected periodic Born--von Karman Coulomb response.  This
ablation does not imply that ACP is generally negligible; it only assigns the
observed Ih--XVII A/B shift.

The coarse-grid MAEs are diagnostic and must not replace the adaptively
converged benchmark.

## Hardening and complete regression gate

The final response-corrected source was additionally hardened in two places.
An explicitly listed regular `GENERAL` mesh is now subjected to the same
strict Cartesian-product inference as the coupled exchange map, so it creates
the identical Born--von Karman Coulomb cell as the equivalent shifted
MacDonald mesh.  A transiently inconsistent set of redundant real-space
density images is projected by group averaging onto the canonical finite-mesh
subspace when a CP2K density-mixer restart is written, rather than aborting a
healthy SCF cycle.

Focused tests give exactly zero printed energy difference both between the
regular `GENERAL` and shifted-grid paths and between CP2K density mixing and
Fock-DIIS.  Reading the resulting same-mesh density-mixer restart succeeds.
The complete `xTB/regtest-tblite-gxtb` suite then passes with 78 correct,
zero wrong, and zero failed checks.  The pre-refresh run is retained to show
that its 17 wrong checks were exclusively stale numerical references: it had
zero runtime failures, and all 17 become correct after refreshing the expected
response-corrected values.

## Contents and verification

- `cli_native_k222.csv`: all thirteen absolute and relative CLI/native values;
- `response_fix_k222.csv`: per-phase old/final relative energies and errors;
- `ablation_summary.csv`: exact-parent, smearing, full-model, and no-ACP gates;
- `summary.json`: machine-readable aggregate metrics;
- `source_identity.txt` and `provenance/`: source and binary identities;
- `raw/`: CP2K outputs, CLI JSON results, inputs, affinity proofs, and exit
  statuses;
- `hardening_validation/`: focused GENERAL-grid, mixer/restart, unit-test, and
  complete before/after regression evidence.  Each regression directory has
  a portable `SHA256SUMS` for the archived proof files; its
  `SOURCE_SHA256SUMS` additionally preserves the original full temporary
  regression-worktree inventory, including files that are not duplicated in
  this compact archive;
- `SHA256SUMS`: complete integrity manifest for every archived file.

Recompute all tables and enforce the numerical gates with:

```bash
python3 verify_response_fix.py
```
