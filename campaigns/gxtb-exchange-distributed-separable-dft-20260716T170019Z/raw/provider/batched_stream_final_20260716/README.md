# Bounded image-batch exchange evidence (2026-07-16)

This directory records the final local qualification of the independent
`save_tblite` worktree at
`/private/tmp/save_tblite_batched_stream_independent_20260716`.

## Build configuration

- Branch: `codex/periodic-exchange-cache`
- Base revision: `257ba442684c39454175e5192c8a2342b4c6380f`
- Compiler: `/opt/homebrew/bin/gfortran`
- CMake build type: `Debug`
- Debug flags: `-g -O0 -fcheck=all -fbacktrace`
- `WITH_API=OFF`, `WITH_OpenMP=OFF`, `WITH_DDX=OFF`
- All FetchContent dependencies were supplied from the existing disconnected
  source cache under `/private/tmp/save_tblite_kmesh_memory/_build-memory/_deps`.

The complete CMake build, including dependency test executables, completed.

## Targeted result

`build/test/unit/tblite-tester exchange` passed 30/30 tests after the final
provider hardening.  The 9x9x1 qualification case has 81 k points, a nonzero
common twist, a coprime permutation of external k-point order, two spin
channels, k-dependent Hermitian density and overlap blocks, explicitly
nonzero R != 0 Fourier images, batch size 7, reverse block arrival, and a
four-image final batch.
The final targeted CTest transcript is preserved verbatim in
`exchange_ctest.log` and included in the SHA-256 manifest.

Measured against the unchanged complete-mesh oracle before removing the
temporary diagnostic print:

- maximum absolute Fock difference: `1.1102230246251565e-16`
- maximum absolute shell-potential difference: `0`
- absolute energy difference: `1.3877787807814457e-16`

## Transaction semantics exercised

Each public batched stream is strictly sequential:

1. `cp2k_exchange_stream_begin(..., mode=cp2k_exchange_stream_batched,
   batch_size=B)`
2. `cp2k_exchange_stream_batch_bounds`
3. exactly one `cp2k_exchange_stream_push` for every k block
4. `cp2k_exchange_stream_batch_apply`
5. exactly one `cp2k_exchange_stream_batch_pull` for every k block
6. `cp2k_exchange_stream_batch_advance`
7. repeat steps 2--6, then call `cp2k_exchange_stream_batch_result` exactly
   once and `cp2k_exchange_stream_end`

An optional `image_first,image_count` pair restricts one stream to a
contiguous additive image range.  Independent streams or MPI ranks can cover
disjoint ranges concurrently.  Each rank then submits the K Bloch blocks only
once when its local batch holds its complete range; the caller sums energy,
shell potential, and Fock contributions over an exact nonoverlapping cover of
images 1..K.  The unit oracle splits the image mesh into two disjoint ranges,
evaluates each with one K-block submission, and verifies their external sum
against the unchanged full-mesh result.  No thread-safety claim is implied for
concurrent use of the same stream object.

The tests cover missing and duplicate push/pull/get operations, pull before
apply, advance before apply, duplicate apply, operations after the final
advance, result before completion, duplicate result, mandatory-result close
failure and recovery, wrong pull dimensions, changed density, changed push
overlap, changed pull overlap, stale `onecxints`, stale `kq`, stale onsite
state, changed image kernel, changed phase table, and changed input-to-grid
ordering.  A nonzero caller-owned Fock sentinel verifies additive pull
semantics.  Plan generation, exact representative order, twist, kernel
fingerprints, phase fingerprints, model signature, and onsite snapshots gate
every apply.  Pulls use only the already validated stream response and
therefore do not repeat an O(Nk) model/plan scan.

Batch mode rejects weights or grid characters that differ from an exact
uniform finite Fourier group by more than `64*epsilon(1.0_wp)`.  This is
stricter than the general whole-mesh parser and prevents the inverse-image
energy contraction from dropping accepted nonzero off-pair Gram coefficients.
The retained inverse-pair coefficient is evaluated from the compact regular
grid character plan.  The grid part cancels exactly for an inverse-image pair,
leaving only the supplied weight sum and its cached twist character.

The compact phase plan stores one roots-of-unity vector per mesh direction and
one twist character per image.  Repeated batch push/pull loops therefore use
only integer modular indexing and complex multiplication: they contain no
transcendental phase evaluations and no dense K-by-K phase table.  The plan is
cached with the geometry/model/kernel fingerprint and invalidated with it.

## Complex-storage high-water mark

For `nao=N`, spin count `S`, k-point count `K`, and allocated batch size `B`,
the exact explicitly allocated complex-scalar count reported by
`cp2k_exchange_stream_peak_complex_elements` is

`4*N^2*B*S + 3*N*S + 2*K + 3*max(nmesh) + max(2*N^2, 5*N^2 + 6*N,
N^2*S + 4*N^2 + B)`.

The test asserts this value for 3x1x1/RKS transactions with B < K, B = K, and
B > K, and for the 9x9x1/UKS transaction at every batch, including the short
final batch.  Dense KxK phase tables are absent in batched mode.  The
conservative full-mesh query is tested false for B < K and true for both B = K
and B > K, where the allocated batch size is clamped to K.

The count deliberately excludes shared immutable real model/image-kernel
caches and compiler-created expression temporaries.  It is consequently a
transaction-owned complex-storage claim, not a total-process RSS claim.

## Remaining boundaries

- A single stream is sequential and no measured MPI scaling claim is made.
  The disjoint image-range API is only the provider primitive needed for
  caller-owned MPI distribution and an additive reduction.
- Every image batch in a serial full-range stream still resubmits all k blocks.
  With O(B*N^2) AO storage and no persistent full Bloch/image mesh, avoiding
  those resubmissions is impossible under the current caller-driven protocol:
  after a block has been released, later image sums contain information that
  cannot be reconstructed from the bounded state.  The remaining matrix-product
  work is therefore proportional to `K*ceil(K/B)`.  A disjoint one-range-per-rank
  decomposition changes the wall-time critical path to one K submission per
  rank without changing that memory bound; it does not constitute a serial
  speedup.  Compact phase and inverse-coefficient setup is O(K+max(nmesh)) in
  time and storage.
- The bounded path covers forward energy, shell potential, and Fock response.
  Force/stress reverse response still uses the complete-mesh oracle and is not
  claimed to be bounded by B.
- Repeated-block consistency uses deterministic 64-bit storage-bit
  fingerprints.  They reliably catch the exercised mutations but are not a
  cryptographic identity proof.
- The high-water API returns `-1` outside an active batched transaction rather
  than under-reporting reduced/oracle modes.

## Broader CTest inventory

After rebuilding every target, `ctest --test-dir build --output-on-failure -j 8`
completed the current 43-test inventory: 39 passed and 4 unrelated tests
failed.  The exchange test passed.  The four failures are documented in
`full_ctest_summary.txt`; none of their diagnostics exercises the new public
image-batch transaction.

## SHA-256

The final source, test binary, focused log, README, and broader CTest summary
hashes are in `source_sha256.txt`.
