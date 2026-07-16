# save_tblite CTest failure classification

## Verdict

No regression from `257ba442684c39454175e5192c8a2342b4c6380f` to
`35e7942b60edd89bb407ab3da5768d3410af83f5` is exposed by the complete current
failure set. The byte-matched base reproduces every selected test name, CTest
status, and diagnostic signature.

| Configuration | 35e7942 | 257ba442 | Classification |
|---|---|---|---|
| Release | C-API failed; tblite/xtbml failed | identical | missing ddX (`WITH_DDX=OFF`) |
| Debug | C-API failed; coulomb-multipole failed; gfn1-xtb NUMERICAL; xtb-param NUMERICAL; xtbml failed | identical | missing ddX plus pre-existing Debug checks/traps |

The Release diagnostic union additionally proves that `coulomb-multipole`,
`gfn1-xtb`, and `xtb-param` pass in both revisions.

## Failure signatures

- `C-API`: four ddCOSMO/ddCPCM/ddPCM calls report that ddX support is unavailable.
- `tblite/xtbml`: only `xtbml-energy-sum-up-gfn1` fails for unavailable ddX.
- Debug `tblite/coulomb-multipole`: the `energy-gfn2-zno-sc` test indexes atomic
  number 30 into test parameter array `p_dkernel(20)` at
  `test_coulomb_multipole.f90:258`; this is exposed by `-fcheck=all`.
- Debug `tblite/gfn1-xtb`: SIGFPE in the unchanged OpenBLAS/LAPACK `syevr` call at
  `src/tblite/lapack/sygvr.f90:280`, exposed by
  `-ffpe-trap=invalid,zero,overflow`.
- Debug `tblite/xtb-param`: SIGFPE at the unchanged generic harmonic mean in
  `src/tblite/xtb/calculator.f90:1619`; a zero pi-channel scale is divided into.

The failure test sources and the two relevant implementation sources are
byte-identical between the revisions. Although `src/tblite/xtb/gxtb.f90` changed,
the `xtb-param` backtrace is in the unchanged generic calculator path.

## Reproducibility

- Base checkout is clean at exact commit `257ba442...`.
- Compiler, OpenBLAS, installed dependencies, CMake flags/options, and environment
  are shared with the current builds.
- Normalized selected CMake-cache diffs are empty for Release and Debug.
- Every test used `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`,
  `MKL_NUM_THREADS=1`, `BLIS_NUM_THREADS=1`, and
  `VECLIB_MAXIMUM_THREADS=1`.
- Raw commands are in `provenance/commands.txt`; raw configure/build/CTest logs
  are in `logs/`; hashes are in `provenance/artifacts.sha256`.

Recommended follow-ups are independent of the partial-k-to-R bridge: either
enable ddX or feature-gate the ddX-only tests; repair the Zn test parameter
indexing; make the generic diatomic scale mean zero-safe; and avoid trapping
floating exceptions raised internally by the selected LAPACK backend unless
that backend behavior is itself under test.
