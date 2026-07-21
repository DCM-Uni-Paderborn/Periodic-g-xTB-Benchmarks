# save_tblite online-head validation

This package records a clean Linux build and test of the online
`cp2k-integration` head
`9322cc6d43be5099f2e7dc4866abc45b8387835a` after reconciling the local and
remote branch tips.  The production source tree at this revision is identical
to the SHA-qualified provider revision used by the Part-I calculations; the
only intervening content change relaxes one CeCl3 finite-difference unit-test
tolerance.

The periodic/g-xTB-focused suite passes all 12 selected CTest targets.  It
includes PBC, Wigner--Seitz, periodic s-dftd3, ACP, exchange, g-xTB,
Hamiltonian, integral transformation, mixing, q-vSZP, and wavefunction
restart coverage.

The broad build-configuration suite passes 91 of 93 CTest targets.  Both
remaining target failures are caused by compiling explicitly with
`WITH_DDX=OFF`: the C-API target still invokes four ddX solvation tests and
the `tblite/xtbml` target still invokes one ddX-dependent energy test.  Each
reports that ddX support is unavailable.  These are configuration-conditioned
test-selection failures, not failures of the periodic or g-XTB implementation.
The focused scientific suite is therefore the acceptance gate for this
provider build; the unfiltered result is preserved rather than hidden.

All configure, build, full-test, focused-test, and exact singleton-affinity
records are retained under `logs/` and `provenance/`.  The status file records
the source commit and both CTest return codes.
