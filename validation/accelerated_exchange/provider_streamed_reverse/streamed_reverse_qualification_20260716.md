# Streamed reverse/gradient qualification (2026-07-16)

## Scope

This qualification covers the memory-bounded reverse path for the
Brillouin-zone-coupled g-xTB exchange provider.  The unchanged complete-mesh
`cp2k_exchange_kmesh_gradient` routine is the dense oracle.

The caller-driven transaction is:

1. `cp2k_exchange_stream_begin(..., reverse_only=.true.)`, or a normal
   reduced forward `push/apply` followed by
   `cp2k_exchange_stream_reverse_begin` without forward Fock pulls;
2. for each image batch, resubmit disjoint indexed k blocks with
   `cp2k_exchange_stream_reverse_push`;
3. call `cp2k_exchange_stream_reverse_batch_apply` once;
4. accumulate one overlap-adjoint block at a time with
   `cp2k_exchange_stream_reverse_batch_pull`;
5. call `cp2k_exchange_stream_reverse_batch_advance` and repeat;
6. obtain forces and the positive homogeneous-strain response exactly once
   with `cp2k_exchange_stream_reverse_result`.

No complete Bloch density, overlap, Fock, or overlap-adjoint mesh is owned by
a batched reverse transaction.  Image-local kernel adjoints are contracted
directly to forces and stress before bounded storage is reused.  Density and
overlap bit signatures guard every resubmission and pull.  Duplicate push,
apply, pull, and result transitions are negative-tested.

## Correctness gates

All commands were run from
`/private/tmp/save_tblite_batched_stream_independent_20260716` with a
controlled thread count:

```text
cmake --build build -j 8
OMP_NUM_THREADS=1 build/test/unit/tblite-tester exchange
cmake --build build-release -j 8
OMP_NUM_THREADS=1 build-release/test/unit/tblite-tester exchange
OMP_NUM_THREADS=1 build-release/test/unit/tblite-tester gxtb
cmake --install build-release
```

Both exchange suites passed 31/31 tests.  The release g-xTB suite passed 40
ordinary tests plus four expected-failure diagnostics.

Maximum absolute streamed-minus-dense residuals are:

| Build | Spin/grid | overlap adjoint | force | stress/strain | peak complex scalars |
|---|---:|---:|---:|---:|---:|
| Release | RKS, 3x1x1, B=2 | 1.3877787807814457e-17 | 0 | 0 | 113 |
| Release | UKS, shifted/permuted 9x9x1, B=7 | 5.5511151231257827e-17 | 0 | 1.3552527156068805e-20 | 585 |
| Release | RKS, shifted/permuted 2x2x2, B=3 | 2.7755575615628914e-17 | 0 | 8.2718061255302767e-25 | 144 |
| Debug/check | RKS, 3x1x1, B=2 | 1.3877787807814457e-17 | 1.3234889800848443e-23 | 2.6469779601696886e-23 | 113 |
| Debug/check | UKS, shifted/permuted 9x9x1, B=7 | 8.3266726846886741e-17 | 5.2939559203393771e-23 | 6.7762635780344027e-21 | 585 |
| Debug/check | RKS, shifted/permuted 2x2x2, B=3 | 2.7755575615628914e-17 | 0 | 1.0339757656912846e-24 | 144 |

The unit-test hard gate is 1e-11 independently of the older supercell test's
1e-9 tolerance.

## Storage bound

For batch size `B`, AO dimension `n`, spin count `s`, mesh size `K`, and
`m=max(nmesh)`, the exact queried peak in explicitly allocated complex
scalars is

```text
Ncomplex_peak = 6 B n^2 s + 5 n s + 2 K + 3 m + 10 n^2.
```

The first two terms are the six bounded image matrices and five onsite
diagonal vectors.  `2 K + 3 m` is the compact regular-grid phase plan.  The
last term is the largest blockwise pull workspace.  The image-kernel reverse
apply additionally owns `2 B n^2` real scalars temporarily for the two kernel
parameter responses.  Thus all AO-matrix storage is O(B n^2 s), while the
only K-dependent transaction metadata is O(K); no O(K n^2) AO mesh exists
when B is fixed.

For the 9x9x1 UKS test (`n=2`, `s=2`, `B=7`, `K=81`, `m=9`), the exact query
is 585 complex scalars and the additional real response workspace is 56 real
scalars.  Increasing K at fixed B changes only the compact O(K) maps.

## Frozen hashes

```text
d3791e769816e28c70b8091b501a187fb3a59083c80588aa2603ff7b86e25a8f  src/tblite/cp2k_compat.f90
17ec6f792419bf060334f04628bfcc151f89520bc73acade6b8682d81ab8a121  src/tblite/exchange/fock.f90
840a32e2fb48084f425b9006c537908a3d89a96456ce1dc530d91effa881ccdf  test/unit/test_exchange.f90
ee665d8de68f5efa21e4b6a7b192f8d0c5b4cc3934d5e65d8c618301c122dd0b  build-release/test/unit/tblite-tester
201f04f51e739d635ca6463ce25cd906135c66e4034261a86752b7be452aa130  build/test/unit/tblite-tester
10d5a11663b5664715f45abfbfa8c1998c10464d31b1e06c25e40dea2a0ddf33  /private/tmp/save_tblite_batched_install/lib/libtblite.a
236722ff47dfa458895213516deb368eb52c35fea5d7f4b1b130a588f30ae2a1  /private/tmp/save_tblite_batched_install/include/tblite/GNU-16.1.0/tblite_cp2k_compat.mod
```

The installed provider prefix is `/private/tmp/save_tblite_batched_install`.
