# Independent adversarial review evidence

## Scope and verdict

This review covered the staged cross-mesh k-point restart implementation in the CP2K worktree `/private/tmp/cp2k_gxtb_crossmesh_restart`, branch `codex/gxtb-crossmesh-restart`.

Verdict: **PASS for the opt-in `VALIDATED_BVK_TRANSFER` initial-guess path**, subject to one explicitly documented limitation: the interpolated metric density is not mathematically guaranteed to be N-representable. The implementation must therefore not be described as positivity-preserving or unconditionally safe, and broad/default activation should wait for a spectral projection or a conservative accept/fallback gate.

No manuscript or SI files were changed during this review.

## Reviewed identity

- base source commit: `28df9380abb327d56bbf216d2469a1fd8c953fc0`
- original staged-diff SHA-256 supplied to the reviewer: `7a33a193cfe2202b1c4cccc2b595b8067bc305294714c053955cb13253422432`
- final staged-diff SHA-256 after review fixes, documentation clarification, and formatting: `d933b478a69ca0db5566d87a4d5b132455ec263081da5cf315de70e469c1add4`
- final local debug binary: `/private/tmp/cp2k_gxtb_crossmesh_restart_build/bin/cp2k.pdbg`
- final local debug binary SHA-256: `c55fd95cc3a28e33cca3a01b5e27aa1baa7fa72973b49dc50a2adcbc788e86f8`
- final patch size: 15 files, 2480 insertions, 270 deletions

The final staged tree is clean apart from the intended staged changes, and `git diff --cached --check` reports no whitespace errors.

## Review fixes applied

The adversarial review found and fixed the following file-validation and compatibility issues:

1. Reject NaN/Inf in version-2 shift, cell, positions, tolerances, smearing, and checksum metadata before any comparisons.
2. Reject non-finite reconstructed target density matrices.
3. Prevent checksum overflow for finite but extremely large payload values on both write and read.
4. Reject trailing records/data after the declared version-2 payload.
5. Preserve legacy behavior through the default `COMPATIBLE` mode:
   - legacy version-1, same-mesh restarts remain accepted by default;
   - version-2, same-mesh restarts remain accepted by default;
   - explicit `EXACT` accepts only self-describing version-2 same-mesh restarts;
   - explicit `VALIDATED_BVK_TRANSFER` is required for cross-mesh interpolation.
6. Clarify log messages and input documentation for the three modes.
7. State explicitly that the accepted-transfer checks cover spin/electron counts,
   smearing metadata, Fourier representation, Hermiticity, and electron trace,
   but do not constitute an N-representability guarantee.

The validation remains transactional: phase state is restored on rejection, and accepted transfers clear CP2K SCF history plus the tblite/g-xTB mixer history rather than transferring incompatible Hamiltonian or DIIS state.

## Build and static checks

Final exact-source rebuild:

```text
cmake --build /private/tmp/cp2k_gxtb_crossmesh_restart_build -j 8 --target cp2k.pdbg
[80/80] Linking Fortran executable ../bin/cp2k.pdbg
```

Checks:

- `git diff --cached --check`: PASS
- CP2K precommit, first uncached pass: 15 files found, 8 skipped, 7 checked, 0 failed, status OK; formatter output was restaged.
- CP2K precommit, final cached pass: 15 files found, 15 skipped from cache, 0 failed, status OK.
- `tools/regtesting/check_inputs.py` for `tests/QS/regtest-kp-restart`: PASS, exit 0.

## Official focused regression

Archived under `official_regtest/`; original run directory:

`/private/tmp/cp2k-crossmesh-independent-regtest-final/TEST-2026-07-16_22-13-32`

Result: 3/3 matchers correct, 0 failed tests, 0 wrong results, status OK.

| Case | Mesh | Transfer | SCF steps | Energy / Ha | Result |
|---|---:|---|---:|---:|---|
| source | `5x5x5` | none | 10 | -2.951453440859938 | PROGRAM ENDED |
| target | `6x6x6` | accepted | 2 | -2.951453440027636 | PROGRAM ENDED |
| marker | -- | `Validated BvK mesh transfer accepted` | -- | -- | matched |

## Five malformed-file adversarial tests

Mutation generator: `mutate_restart.py`. Each final exact-source run is archived as `post_fix/<mode>/test_final_exact.out`.

| Mutation | Required behavior | Observed behavior | Cold SCF steps | Energy / Ha |
|---|---|---|---:|---:|
| `nan_shift` | reject | `version-2 metadata contains non-finite values` | 10 | -2.951453440874063 |
| `nan_hmat` | reject | `version-2 metadata contains non-finite values` | 10 | -2.951453440874063 |
| `nan_checksum` | reject | `version-2 metadata contains non-finite values` | 10 | -2.951453440874063 |
| `huge_payload` (`1e200`) | reject | `payload magnitude would overflow its checksum` | 10 | -2.951453440874063 |
| `append_record` | reject | `trailing data after version-2 payload` | 10 | -2.951453440874063 |

All five completed with `PROGRAM ENDED`; none accepted the corrupted restart. Pre-fix evidence is retained under `pre_fix/`: trailing data was accepted, while NaN metadata reached IEEE-trapping debug failures instead of a controlled fallback.

## Default and explicit-mode compatibility

The archived legacy fixture has SHA-256 `deaca45d972f70e2fb2e4da70ef80d40655c6ebdac5edbd1870d398372c4c687`.

| Fixture/mode | Expected | Observed | SCF steps | Energy / Ha | Wall / s | Max RSS / B |
|---|---|---|---:|---:|---:|---:|
| legacy v1 / default | accept | `Legacy-compatible restart accepted` | 1 | -4.349166473124775 | 0.78 | 119177216 |
| legacy v1 / `EXACT` | reject/fallback | rejected: only `COMPATIBLE` accepts v1 | 12 | -4.349166468899968 | 3.78 | 146751488 |
| v2 / default | accept | `Strict same-mesh restart accepted` | 1 | -4.349166473124775 | 0.77 | 119259136 |
| v2 / `EXACT` | accept | `Strict same-mesh restart accepted` | 1 | -4.349166473124775 | 0.77 | 121470976 |

This verifies that merely adding the new feature does not silently disable legacy same-mesh restart files.

## RKS energy, force, stress, and final-density comparison

System: two asymmetric He atoms in a `4.0 x 4.3 x 4.7` Angstrom orthorhombic cell, PBE/K290, source `3x3x3`, target `4x4x4`.

Primary outputs:

- transferred: `force_stress/target_restart_force_mokp.out`
- cold: `force_stress/target_cold_force.out`
- exact parsed comparison: `force_stress/comparison.txt`

Results:

- transferred run accepted the restart and converged in 7 SCF steps;
- cold run converged in 15 SCF steps;
- energy delta: `-8.881784197001252e-16 Ha`;
- maximum absolute force delta: `9.999999717180685e-10 Ha/bohr`;
- maximum absolute stress delta: `1.499999780207872e-05 bar`;
- final restart-file density comparison: maximum absolute delta `1.738517871330281e-10`, RMS `4.161419153312078e-12`, relative Frobenius delta `1.177179419484646e-10`.

Transferred force vector, in Ha/bohr:

```text
atom 1  -0.412085241  -0.831141456   0.312564075
atom 2  -0.407316372   0.000115389907  -0.191595378
```

Transferred stress tensor, in bar:

```text
 -2210.51997427    1873.01224982     130.496085162
  1873.01224982  221427.910070      -910.207082495
   130.496085162   -910.207082495  110329.619479
```

## UKS and final-density comparison

System: one He atom, UKS singlet (`1` alpha and `1` beta electron), source `3x3x3`, target `4x4x4`.

| Case | SCF steps | Energy / Ha |
|---|---:|---:|
| source | 11 | -2.951453440459917 |
| target transferred | 6 | -2.951453435945473 |
| target cold | 11 | -2.951453436494849 |

The transferred run accepted the restart and completed normally. Its final density versus the cold target has maximum absolute delta `8.456397040945696e-12`, RMS `4.911276352010702e-13`, and relative Frobenius delta `1.949211790267738e-11`.

## Two-rank Linux MPI validation on terok

An isolated exact-base-plus-patch tree was built on terok:

- source: `/home/kuehne88/work/cp2k-crossmesh-independent-review`
- binary: `/home/kuehne88/work/cp2k-crossmesh-independent-review-build/bin/cp2k.psmp`
- remote evidence root: `/home/kuehne88/work/cp2k-crossmesh-independent-review-mpi2`
- local copied logs: `mpi2/`

| Case | MPI ranks | Mesh | Transfer | SCF steps | Energy / Ha | Wall / s |
|---|---:|---:|---|---:|---:|---:|
| source | 2 | `5x5x5` | none | 10 | -2.951453440859939 | 3.274605726 |
| target | 2 | `6x6x6` | accepted | 2 | -2.951453440027636 | 2.974859933 |

Both logs contain `PROGRAM ENDED`; the target contains `Validated BvK mesh transfer accepted`.

## Metric-density occupation audit and remaining limitation

`diagnose_occupations.py` uses CP2K's high-precision `MO_KP ... OVERLAP_MATRIX` output, reconstructs the source BvK density at every target irreducible k point, applies the implementation's scalar electron-number normalization, and diagonalizes `S^(1/2) P_sigma S^(1/2)`.

RKS audit (`force_stress/occupation_diagnostics.txt`):

- Hermiticity relative residual: `6.960137652469458e-18`
- minimum overlap eigenvalue: `0.2982695380098867`
- trace before normalization: `3.999999613804153`
- normalization factor: `1.000000096548971`
- normalized metric occupations: `[-1.109384066673194e-5, 2.000959228752370]`
- bound violations: lower `1.109384066673194e-5`, upper `9.592287523703114e-4`

UKS audit (`uks_density/occupation_diagnostics.txt`):

- Hermiticity relative residual: `5.421010862427522e-20`
- minimum overlap eigenvalue: `0.7071537542418023`
- trace before normalization: `1.000000007788061` per spin
- normalization factor: `0.9999999922119395`
- normalized metric occupations: `[-4.885708932167015e-7, 1.000050668119747]`
- bound violations: lower `4.885708932167015e-7`, upper `5.066811974718810e-5`

Therefore scalar trace normalization alone is insufficient to guarantee positivity and Pauli bounds. This did not affect the tested final observables, because the SCF healed the interpolated initial guess, but it is a real algorithmic qualification. Recommended next step before broad activation: project occupations in the overlap metric while conserving the electron count, or reject/fall back when an inexpensive spectral bound exceeds a documented tolerance.

## Code-path review notes

The following invariants were inspected in addition to dynamic testing:

- fingerprinting covers basis topology, primitive exponents and contractions, species, geometry, cell, symmetry, mesh/shift, spin, smearing, and payload checksums;
- regular Monkhorst-Pack/MacDonald meshes, canonical cells/aliases, Hermiticity, electron trace, and Nyquist representability are guarded;
- restart rejection is a cold-start fallback, not a fatal error;
- only the density is transferred; no Fock, DIIS, or mixer history is reused;
- CP2K and tblite/g-xTB histories are cleared after acceptance;
- pseudopotential and XC identity are not explicitly fingerprinted, which is acceptable for a non-authoritative initial guess but should be mentioned if the semantics are ever strengthened beyond that.

## Reproduction helpers and archive

- malformed restart generator: `mutate_restart.py`
- restart parser/density comparator: `parse_restart.py`
- force/stress comparator: `compare_force_stress.py`
- metric-occupation audit: `diagnose_occupations.py`
- archive integrity: `SHA256SUMS`

The absolute paths above describe the review environment. The hash manifest provides content identity if the evidence is moved into the benchmark repository.
