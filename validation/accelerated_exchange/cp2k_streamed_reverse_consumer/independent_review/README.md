# Independent review: CP2K streamed g-xTB reverse consumer

Date: 2026-07-16

## Verdict

**PASS, with a documentation-only memory-scope correction.** No functional source defect was found
and no executable statement was changed by this review. The selectable streamed reverse path agrees
with the unchanged dense derivative oracle for energy response, overlap adjoint, forces, and stress
in all frozen Linux cases. The original wording could be read as claiming `K`-independent total
CP2K memory; that claim is false and has been corrected in the source comments and evidence README.

## Code and ABI audit

- `DENSE` remains the default and unchanged numerical oracle.
- `STREAMED` keeps CP2K's irreducible density and response matrices, reconstructs one full-star
  provider batch, pulls its overlap adjoints, and folds them immediately. It does not retain an
  additional expanded full-star density/overlap/adjoint mesh.
- `QUALIFY` applies the streamed response exactly once and separately constructs the dense meshes
  only for comparison; it does not double-apply a force, stress, or Fock response.
- Direct atomic and homogeneous-strain derivatives are retrieved once after the final reverse batch.
- The same CP2K symmetry adjoint and k-to-R/R-to-k machinery is used by the streamed and dense paths.
- Spin-summed exchange response is added to CP2K channel 1 in both RKS and UKS, matching the dense
  implementation; the UKS full-mesh case passes.
- Every participating rank follows the same collective sequence. The P=4 true-1D case has more MPI
  ranks than full k points and completes without deadlock.
- The frozen public `save_tblite` reverse-stream state machine rejects duplicate push, apply, pull,
  result, missing-pull, and stale-state misuse; Debug and Release provider suites both pass 31/31.

## Exact numerical maxima

The independent verifier checks every per-run SHA-256 manifest, return code, `PROGRAM ENDED` record,
input hash, qualification record, high-water equality, finite-difference result, and available
P=1/P=2/P=4 printed observable.

| Check | Maximum |
| --- | ---: |
| Streamed vs dense overlap adjoint | `1.85962356624714e-15` |
| Streamed vs dense direct force | `1.08420217248550e-19` |
| Streamed vs dense direct stress | `3.46944695195361e-18` |
| Finite-difference force sum | `2.0e-8` hartree/bohr |
| Finite-difference stress sum | `1.30618e-7` a.u. |
| P=1/P=2/P=4 printed energy delta | `7.105427357601002e-15` hartree |
| P=1/P=2/P=4 printed force delta | `1.59159765e-16` hartree/bohr |
| P=1/P=2/P=4 printed stress delta | `2.07317494193e-10` bar |

The canonical matrix contains 21 passing runs and 39 passing dense-oracle evaluations: RKS and UKS,
K290, shifted and time-reversal meshes, explicit full meshes, SPGLIB, true 1D and 2D periodicity,
forces, stress, and finite differences, each at P=1,2,4. The old orphan directory `k290_fd_p1`
has no canonical manifest and is intentionally excluded; it is not required by or counted in the
21-run status matrix.

## Provider memory bound

For batch size `B`, AO dimension `n`, spin dimension `s`, complete mesh size `K`, and regular mesh
dimensions `nmesh`, the audited complex-scalar high-water invariant is

```
6 B n^2 s + 5 n s + 2 K + 3 max(nmesh) + 10 n^2.
```

This exactly covers the explicitly allocated complex arrays owned by the provider transaction plus
its largest explicit complex pull workspace. It excludes existing CP2K matrices, provider real and
integer metadata, immutable shared caches, compiler-created temporaries, and a transient real apply
workspace. Therefore the supported statement is:

- CP2K irreducible input/output storage remains `O(nred n^2 s)`; with `FULL_GRID`, `nred = K`.
- The consumer avoids an *additional* expanded full-star mesh.
- The provider AO-image workspace is fixed-`B`, `O(B n^2 s)`.
- Total CP2K process memory is not proven independent of `K`.

The largest single-process RSS in the correctness matrix was 603144 KiB. In the like-for-like
shifted-SPGLIB P=1 run, dense and streamed used 601676 and 601548 KiB process RSS, respectively
(625152 and 624752 KiB process-tree RSS). The streamed wall time was 37.751557909 s versus
37.276223755 s for dense. This small case supports neither a material whole-process memory-saving nor
a speedup claim; it only confirms numerical equivalence and the exact provider-workspace invariant.

## Frozen artifacts and hashes

| Artifact | SHA-256 |
| --- | --- |
| Terok-tested `src/tblite_interface.F` | `3ac492c0070137547bbd897fe82eeab312953385de0a1701799abd386b46ac93` |
| Reviewed source after comment-only correction | `bb7c755478333be0113b280ee198bf8d02ce90b4bbf1c5113387195d6dfefeb1` |
| `src/tblite_types.F` | `e848fd796e69021a1ba21473d0fd7e2da3dbe485fe0e799e495ecd916889398e` |
| Frozen pre-format source copy | `890135629af5ed9dad8f702b8a42d8b5add30ab2dfb78a12221f0d9d72f5a4fb` |
| Terok `cp2k.pdbg` | `8cdecb1e925ef1a2e6391da7298bd720d7bf4f335ab8d405127097db9e53033d` |
| Terok `libtblite.a` | `0095c215ed46318cf3e368f06b425829a82fdd40e4a8ceed6661749a496adb2b` |
| Provider `cp2k_compat.f90` | `d3791e769816e28c70b8091b501a187fb3a59083c80588aa2603ff7b86e25a8f` |
| Provider Debug tester | `201f04f51e739d635ca6463ce25cd906135c66e4034261a86752b7be452aa130` |
| Provider Release tester | `ee665d8de68f5efa21e4b6a7b192f8d0c5b4cc3934d5e65d8c618301c122dd0b` |
| Independent verifier | `9959696258459be38e7328078a3c63b99b49d78a8399bc1e223b756db31c37db` |

The portable verifier reads
`../raw_archive/cp2k_gxtb_streamed_reverse_consumer_evidence_20260716.tar.gz`,
rejects archive links or paths that escape its temporary extraction root, and
then verifies the canonical `linux_matrix_terok/` matrix and
`linux_mode_rss_terok/` dense/streamed pair directly.  No original
`/private/tmp` directory is required.

The final local Fortran formatter reported 0 files changed, both source property checks passed, and
`git diff --check` is clean. The remote precommit service was unreachable during the last re-run;
the complete precommit pass had already succeeded before the final comment-only correction.
