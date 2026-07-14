# g-xTB PBC benchmark campaign V1

This directory freezes the first publication benchmark campaign for the
CP2K/save_tblite periodic g-xTB integration. V1 is the explicitly expanded
full-mesh implementation behind CP2K's ordinary SPGLIB-reduced user interface.
The later memory-oriented V2 implementation must use a separate campaign and
must reproduce the frozen V1 results before it can replace them.
The deferred symmetry, FFT, cache, MPI, restart, validation, and manuscript
contract is recorded separately in `../gxtb-pbc-v2/README.md`; it does not
alter this frozen V1 campaign.

The existing Overleaf project *Advancing GFN2-xTB for Periodic Systems via
multipolar Ewald Summation in CP2K* is a read-only protocol reference. Its
checkout at `/tmp/gfn2_overleaf_reference_20260714` (commit `eb9f475`) is not a
result destination and was not modified.

## Frozen software

- CP2K `g-xTB-pbc`: `18d37c946413dba1b848f57563c46d16b866ce20`
- save_tblite `cp2k-integration`:
  `1449febde312874cd0fac4227919f5ba4e4b69b8`
- CP2K launcher SHA256:
  `c3c19763eef2fcd51b0266db811be20d0a3d13a586c1bbd4eb0e83c59cca2c8a`
- loaded `libcp2k.2026.1.dylib` SHA256:
  `a1ab82829a43d872b32d2a4f9b929c394f96e5dd0d293db163716c946ad3b8b3`
- save_tblite CLI SHA256:
  `8b0a1d5acd36df23efeb9078dceced61adb95a3cc8523f25ff76114da602cfa5`
- installed static `libtblite.a` SHA256:
  `8ce501b698e1b7ea5db03bc4db979302dc08253e843ee44cb5ef2aecb5a4668c`

Both source worktrees were clean at build time. The remote DCM branches point
to these exact commits. `DCM-Uni-Paderborn/save_tblite:pbc` and
`lmseidler/save_tblite:pbc` both point to `083f22030f0be7abff3e3d27b35c141c04c2aa6d`.

## Production contract

All explicit g-xTB k meshes use `SYMMETRY T`, `FULL_GRID F`,
`SYMMETRY_BACKEND SPGLIB`, and `SYMMETRY_REDUCTION_METHOD SPGLIB`. CP2K
expands the density and overlap to the coupled full mesh internally, invokes
save_tblite once for that mesh, then folds the response back to the irreducible
representation. Full-grid calculations are validation oracles only.

Production SCC uses save_tblite's native Fock DIIS (`SCC_MIXER TBLITE`). The
CP2K density and Fock mixers remain regression-tested alternatives, not the
benchmark default. DDX is disabled, so continuum-solvation calculations are
outside this campaign.

The benchmark matrices are:

- DMC-ICE13: 13 phases at implicit Gamma and k111/k222/k333/k444/k555;
  k333 is the primary comparison and denser meshes establish convergence.
- X23b: 23 gas GEO_OPT, 23 Gamma CELL_OPT, 23 k222 SPGLIB CELL_OPT, then
  k333 and k444 single points on the k222 geometry.
- LC12: 13 save_tblite atomic references plus 13 CP2K atom checks, 12 solids
  on an 11-point k444 EOS, then k333/k444/k555 energies at each accepted EOS
  minimum.

The DMC runner additionally supports explicit `k666` through `k131313`
convergence extensions.  They are outside the frozen six-mesh/78-job
core, do not change the campaign manifest or any prior stamp, and are admitted
only through per-record build identities and per-output hash gates.  A
same-source build on another architecture requires an explicit qualified
execution-build manifest; the default remains the exact frozen artifact gate.
The manifest is itself hash-bound and must carry hashed, finite, numerically
consistent same-mesh dense sentinel evidence schema v3 at tolerances no looser
than `1e-10` Eh and `0.001` kJ mol-1 per H2O.  Every sentinel binds the phase
and Ih inputs, four completed outputs (remote/reference phase and Ih), and all
four corresponding run stamps with safe relative paths and SHA256 hashes.
`remote_build_id` must equal the
qualified execution-manifest identity and `reference_build_id` the exact
frozen campaign identity.  The readers hash and parse the same bytes, require
exact phase/mesh/project g-xTB/SPGLIB input semantics, derive both water counts
from the explicit O atoms, and bind every output to its exact CP2K input and
project header, the `tblite_gxtb` capability, a matching 7--40-hex CP2K source
revision prefix, and the exact alternate-build save_tblite revision.  Remote
stamps are schema v2 and bind the full build, input, frozen-input, and output
identity; frozen reference artifacts must come from the explicitly
SHA256-pinned base index.  Phase and Ih evidence must be distinct, as must each
remote/reference output hash.  The manifest-declared Terok start host, compile
host, and x86_64 platform are checked from the same cached output bytes.  These
environment fields are accidental-copy evidence, not cryptographic execution
attestation; the SSH transfer remains manually trusted.  The sentinel records
the fixed `hartree_to_kjmol = 2625.499638`; the readers parse the four outputs,
verify separate phase and Ih cross-build total-energy deltas, and recompute the
remote and reference Ih-referenced relative energies and their delta from raw
totals.  The declared observed total-energy maximum covers both the phase and
Ih deltas across all sentinels.  Partial dense pilots remain
diagnostic for the separate fixed-mesh stopping report.  The
phase-wise result selects, for each non-reference phase, the smallest adjacent
`N^3` mesh whose same-mesh-Ih-referenced relative energy changes by at most
0.05 kJ mol-1 per H2O from `(N-1)^3`; it is complete when all 12 phases are
selected.  RMS, mean absolute value, and maximum of those 12 last deltas are
reported as diagnostics only.  The runner serializes invocations with a
nonblocking lock and feeds analysis an immutable content-addressed validation
snapshot that binds both the generated and actually executed input copies.

Every accepted output must carry one homogeneous fingerprint for the CP2K
launcher, loaded libcp2k, embedded/source revision, and static libtblite
archive.  A schema-v2 validation index may combine records from qualified
ARM64 and x86_64 builds only because each record independently names and
verifies that fingerprint while all records retain the exact frozen source
revisions.  Schema-v1 snapshots remain byte-for-byte immutable and readable.
Stale, unqualified mixed-build, old full-grid, or unstamped outputs are
diagnostics and never enter the paper tables.

## Result locations

Raw outputs remain in the ignored run trees so they are preserved locally
without replacing the published GFN1/GFN2 tables. Curated CSV/JSON results and
their provenance are additive files in the corresponding `data` directories:

- `DMC-ICE13/runs_gxtb_spglib` and `DMC-ICE13/data/gxtb_spglib_*`
- `X23b/runs` and `X23b/data/gxtb_staging`
- `Goldzak12/runs` and `Goldzak12/data`

`build_manifest.json` is the machine-readable software and protocol freeze.
The per-benchmark provenance files remain authoritative for individual input
and output hashes.
