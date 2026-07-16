# Final CP2K Root-A partial-k-to-R bridge evidence

This is the immutable evidence bundle for the final `KGROUP_PARTIAL_ROOT`
bridge between CP2K and `save_tblite`.  It was assembled without changing any
pre-existing raw archive, verified recursively, published first on its named
evidence branch, and then fast-forwarded into the private benchmark
repository's default branch.

## Qualified source identities

- CP2K base: `0a1f7e3329a3e6c2a6accff28617af53fb9943b4`
- CP2K review branch: `codex/gxtb-partial-root-bridge-clean`
- Only changed CP2K file: `src/tblite_interface.F`
- CP2K source SHA-256:
  `47c9b039b2e0d081f1ac3688f29f5c75ffed9a60acbb490f7bee1ae99593dd5d`
- Patience diff SHA-256:
  `8763981ef9c6ba7e9db26a11f4245e11fbcc52ac393a6b6421b09e8e7628156f`
- `save_tblite` provider commit:
  `35e7942b60edd89bb407ab3da5768d3410af83f5`
- Provider source-tar SHA-256:
  `2535519767302bc851c30512852e4ca031fadc0f97d0ece860e983effedbfd28`

Exact build-product identities and original paths are recorded in
`provenance/SOURCE_AND_BUILD_IDENTITY.md`.  The CP2K and provider CMake caches,
Ninja link graphs, and complete build logs retained here prove the selected
configuration and final successful links.

## Validation verdict

1. Static review found no remaining correctness or fail-closed defect in the
   bridge diff.  The review and its explicit scope are in
   `provenance/static-review/STATIC_REVIEW.md`.
2. The focused provider exchange suite passed all 31 tests in both Release and
   Debug builds.
3. The provider-wide CTest failures were reproduced byte-for-byte at the exact
   base revision `257ba442684c39454175e5192c8a2342b4c6380f`.  They are therefore
   classified as pre-existing or caused by `WITH_DDX=OFF`, not by the partial
   bridge.  Raw logs and a complete manifest are under
   `provider/full-ctest-base-classification/`.
4. Release and Debug CP2K positive smokes passed for a nontrivial `2x2x2` mesh
   at MPI P=1,2,4.  Explicit-Gamma Mode 6 agrees with explicit-Gamma dense;
   its derivative deliberately uses CP2K's one-point dense fallback.  The tiny
   explicit-vs-implicit Gamma baseline is unchanged by tighter settings and is
   not a Mode-6 discrepancy.
5. All five Release fault-hook cases aborted without timeout and with the
   required diagnostic: nonleader failure, anti-Hermitian reverse response,
   non-finite forward result, unknown selector, and a truncated selector.
6. The frozen independent oracle matrix passed 20/20 DENSE-vs-partial-root
   pairs: 40/40 CP2K outputs returned zero and contain exactly one
   `PROGRAM ENDED`.  Coverage includes RKS and UKS, full-grid, K290 and SPGLIB
   symmetry reduction, time reversal including P>nred, 1D, 2D, shifted 3D
   meshes, forces, and analytical stress.  Every in-process forward/reverse
   dense-oracle residual is at most `1e-10`; independently parsed observables
   meet the stated energy, force, and stress gates.  See
   `validation/oracle-harness/summary.tsv` and its raw runs.

## Directory map

- `cp2k/source/`: reviewed source and replayable patch
- `cp2k/build/{release,debug}/`: configure/build logs, CMake caches, link graphs
- `provider/source/`: qualified provider source tar and changed-file hashes
- `provider/build/{release,debug}/`: provider build identities and 31-test logs
- `provider/full-ctest-base-classification/`: base/current classification, raw
  logs, normalized caches, commands, and internal SHA-256 manifest
- `validation/final-smokes/`: Release/Debug positive runs, Gamma diagnostics,
  fault hooks, inputs, scripts, raw output, and per-campaign manifests
- `validation/oracle-harness/`: frozen 40-run oracle campaign, verifier,
  tamper self-test, raw outputs, summary, provenance, and SHA-256 manifest
- `provenance/`: static review, build/link proof, source/build identities
- `SHA256SUMS`: top-level manifest for the complete bundle

## Integrity verification

From this directory:

```sh
shasum -a 256 -c SHA256SUMS
(cd provider/full-ctest-base-classification && \
  shasum -a 256 -c provenance/artifacts.sha256)
(cd validation/oracle-harness && shasum -a 256 -c SHA256SUMS)
```

The final-smoke manifests preserve their original absolute Terok paths.  Their
62 entries rooted at
`/home/kuehne88/work/gxtb-partial-root-final-35e7942-20260716` were remapped to
`validation/final-smokes/` during bundle verification; all 62 matched, with no
missing entry.

## Proposed archive destination

If this staging bundle is approved for archival, use a new branch in the
private `DCM-Uni-Paderborn/Periodic-g-xTB-Benchmarks` repository:

`evidence/root-a-partial-bridge-final-35e7942-20260717`

The branch has intentionally not been created, committed, or pushed here.
