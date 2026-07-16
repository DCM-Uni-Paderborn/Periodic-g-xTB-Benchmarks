# Provider-bounded g-xTB reverse consumer qualification

Date: 2026-07-16

This directory qualifies the CP2K consumer of the frozen `save_tblite` bounded reverse API. It is
method-development evidence only: no LC10, DMC13, or other application benchmark is included.

## Implementation under test

`CP2K_GXTB_EXCHANGE_GRADIENT_MODE` selects three deliberately separate paths:

- `DENSE` retains the Part-I complete-mesh derivative unchanged.
- `STREAMED` stores the irreducible CP2K density and overlap response, regenerates one bounded batch
  of full-star density/overlap blocks, and folds every returned overlap adjoint immediately.
- `QUALIFY` applies the streamed response once, then constructs the complete meshes solely for a
  dense numerical oracle and aborts above a relative/absolute `1.0e-10` gate.

For `K > 1`, the batch is forcibly smaller than the complete mesh. The provider must report that it
owns no complete AO mesh. Its exact complex-element high-water mark is checked at runtime against

```
6 B n^2 s + 5 n s + 2 K + 3 max(nmesh) + 10 n^2,
```

where `B` is the active image batch, `n` the AO dimension, `s` the spin dimension, and `K` the full
k mesh. Atomic and homogeneous-strain derivatives are retrieved exactly once after the last image
batch. The established one-point Gamma derivative remains on the dense path.

This is an exact bound for the explicitly allocated complex workspace owned by the provider; it is
not a claim that the complete CP2K process is independent of `K`. CP2K retains its irreducible input
and response matrices, `O(nred n^2 s)`, and `FULL_GRID` has `nred = K`. The streamed path avoids an
*additional* expanded full-star mesh in the consumer and bounds the provider AO-image workspace by
`O(B n^2 s)`. Existing CP2K k-point/real-space output matrices, real and integer provider metadata,
immutable caches, and compiler temporaries are outside the complex-element counter above.

The helper scripts are:

- `run_linux_matrix.sh`: disjoint CPU sets, MPI `P=1,2,4`, and the full correctness matrix;
- `run_with_rss.py`: 20-ms `/proc` sampling of the complete process-tree resident set;
- `run_linux_mode_rss.sh`: like-for-like dense/streamed process-level RSS comparison.

## Frozen sources and build

- Local consumer worktree: `/private/tmp/cp2k_gxtb_streamed_reverse_consumer`
- Local build: `/private/tmp/cp2k_gxtb_streamed_reverse_consumer_build`
- Terok consumer root: `/home/kuehne88/work/cp2k-gxtb-streamed-reverse-consumer-20260716`
- Frozen provider source: `/private/tmp/save_tblite_batched_stream_independent_20260716`
- Terok provider root: `/home/kuehne88/work/codex-gxtb-streamed-reverse-frozen-20260716`
- Provider qualification:
  `/private/tmp/save_tblite_batched_stream_independent_20260716/evidence/streamed_reverse_qualification_20260716.md`

Terok used GCC/GFortran 15.2.0, Open MPI, ScaLAPACK, SPGLIB, one OpenMP thread per MPI rank, and a
Debug/check build with leak reporting disabled for the known Open-MPI finalization allocations. The
build advertises `tblite_gxtb` and `spglib`.

Frozen SHA-256 values:

| Artifact                 | SHA-256                                                            |
| ------------------------ | ------------------------------------------------------------------ |
| Terok `cp2k.pdbg`        | `8cdecb1e925ef1a2e6391da7298bd720d7bf4f335ab8d405127097db9e53033d` |
| Terok-tested `src/tblite_interface.F` | `3ac492c0070137547bbd897fe82eeab312953385de0a1701799abd386b46ac93` |
| Reviewed source after documentation-only correction | `bb7c755478333be0113b280ee198bf8d02ce90b4bbf1c5113387195d6dfefeb1` |
| `src/tblite_types.F`     | `e848fd796e69021a1ba21473d0fd7e2da3dbe485fe0e799e495ecd916889398e` |
| Terok `libtblite.a`      | `0095c215ed46318cf3e368f06b425829a82fdd40e4a8ceed6661749a496adb2b` |

The final precommit pass changed only continuation indentation in two oracle-validation calls. An
incremental Terok rebuild retained the byte-identical `cp2k.pdbg` hash shown above. The subsequent
review correction changed comments and this README only; no executable statement changed, so the
numerical evidence still applies to the reviewed source.

## Linux MPI result matrix

All 21 runs returned zero, printed `PROGRAM ENDED`, and printed the bounded qualification record.
The 39 derivative evaluations all satisfied the exact provider high-water formula.

| Reference case        | Coverage                                  |    P=1 wall / tree RSS |    P=2 wall / tree RSS |    P=4 wall / tree RSS |
| --------------------- | ----------------------------------------- | ---------------------: | ---------------------: | ---------------------: |
| K290 RKS 3D FD        | K290, force, stress                       |  62.899 s / 220032 KiB |  55.628 s / 410976 KiB |  53.915 s / 791080 KiB |
| Shifted SPGLIB RKS 3D | shifted mesh, SPGLIB, force, stress       |  40.782 s / 626464 KiB | 32.989 s / 1096852 KiB | 24.400 s / 2002620 KiB |
| TR RKS 3D FD          | time reversal, force, stress              | 124.978 s / 247520 KiB | 102.869 s / 470844 KiB |  84.779 s / 914284 KiB |
| Full-mesh RKS 3D FD   | explicit full mesh, force, stress         | 124.091 s / 248556 KiB | 104.530 s / 469396 KiB |  85.157 s / 911196 KiB |
| Full-mesh UKS 3D FD   | UKS, explicit full mesh, force, stress    |  27.469 s / 207880 KiB |  27.789 s / 385752 KiB |  27.721 s / 739820 KiB |
| RKS 1D FD             | true `PERIODIC X`, force, diagonal stress |   5.362 s / 190440 KiB |   6.102 s / 355432 KiB |   6.636 s / 686408 KiB |
| SPGLIB RKS 2D FD      | true `PERIODIC YZ`, SPGLIB, force, stress |  39.894 s / 224908 KiB |  33.961 s / 413640 KiB |  35.430 s / 785036 KiB |

Across every rank count and case, the largest streamed-versus-dense residuals were:

- overlap adjoint: `1.85962356624714e-15`;
- direct force derivative: `1.08420217248550e-19`;
- direct stress derivative: `3.46944695195361e-18`.

The largest finite-difference sum was `1.30618e-7` for stress and `2.0e-8` for force; the printed
force relative error was `0.00%`. The RSS values above are sums over each complete process tree, so
their near-linear increase with `P` reflects replicated CP2K rank state rather than a larger
per-rank streamed transaction. The maximum single-process RSS over the matrix was 603144 KiB.

Raw Linux evidence, including input hashes, return codes, selected records, complete output, per-run
SHA-256 manifests, and RSS JSON, is frozen locally at
`/private/tmp/cp2k_gxtb_streamed_reverse_consumer_evidence/linux_matrix_terok` and remotely below
`.../evidence/linux_matrix` in the Terok consumer root.

The complete frozen evidence archive is
`/private/tmp/cp2k_gxtb_streamed_reverse_consumer_evidence_20260716.tar.gz`, SHA-256
`c934c892f9ba105deb270b69339c9d18433e4ff84cc391d4beea718b36ae3d0f`. The identical archive is
stored in the Terok consumer root.

## Like-for-like dense versus streamed RSS

The shifted-SPGLIB P=1 `ENERGY_FORCE` case was repeated without the qualification oracle:

| Mode          |     Energy / hartree |     Wall |   Tree RSS | Single-process RSS |
| ------------- | -------------------: | -------: | ---------: | -----------------: |
| Dense         | -579.050838553817130 | 37.276 s | 625152 KiB |         601676 KiB |
| Streamed, B=2 | -579.050838553817130 | 37.752 s | 624752 KiB |         601548 KiB |

The printed energy is identical. Maximum printed differences are `6.11e-17` hartree/bohr in force
components and `4.37e-10` bar in off-diagonal stresses whose reference values are themselves at
roundoff; diagonal stress is identical. The process-level RSS reduction is only 400 KiB in this
small test because the CP2K SCF/full-matrix state dominates. Provider-workspace boundedness is
therefore established by the exact high-water invariant. The complete CP2K process still retains
`O(nred n^2 s)` irreducible matrices; a larger method benchmark is required before claiming a
material whole-process memory reduction.

Raw data is frozen at
`/private/tmp/cp2k_gxtb_streamed_reverse_consumer_evidence/linux_mode_rss_terok`.

## Excluded negative evidence

Two failures are retained but are not mixed into the positive statistics:

1. The first Terok configure selected `CP2K_TBLITE_PROVIDER=UPSTREAM`. Consequently `__TBLITE_GXTB`
   was absent and the incompatible legacy interface failed to compile. Reconfiguring the unchanged
   source with `CP2K_TBLITE_PROVIDER=SAVE` produced `g-xTB support: 1` and the clean 2710/2710 link
   above. The discarded and accepted build logs are both archived in
   `/private/tmp/cp2k_gxtb_streamed_reverse_consumer_evidence/terok_build`.
1. A macOS H2 time-reversal debug run converged its SCF, then received SIGILL in the runtime before
   the derivative qualifier. It has no `PROGRAM ENDED` record and is excluded. Its raw output is
   retained at `/private/tmp/cp2k_gxtb_streamed_reverse_consumer_evidence/tr_311_fd/run.out`,
   SHA-256 `42dc232aeaa5d3dae73c8e01f3310fca0e93743cebea2c688a857b53d36a5c40`. The same TR input
   passed on Linux for P=1,2,4 and is the authoritative result.

## Scope statement

This matrix proves correctness and fixed-`B` provider ownership for the selectable CP2K consumer,
and proves that the streamed consumer does not materialize an additional expanded full-star mesh.
It does not prove `K`-independent total CP2K memory: the irreducible CP2K matrices remain resident,
and `FULL_GRID` makes that set the complete mesh. It also does not by itself establish distributed
nonlocal-exchange scaling across k groups: in this path the provider transaction remains replicated
per participating MPI rank. True k-group decomposition and partial-k-to-R accumulator merging
require a separate provider/CP2K qualification and must not be inferred from the P=1/2/4 correctness
runs above.
