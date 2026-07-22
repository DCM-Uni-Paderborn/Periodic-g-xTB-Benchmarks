# Exact BvK cache-identity qualification

This archive closes the persistent Born--von Karman (BvK) exchange-plan
identity gate for the Part-II acceleration tree.  The production cache already
stores and compares the complete model signature and the ordered active BvK
plan.  Provider commit `5c430c5aa32ce84ca940b26d66f0482c517b3805`
adds a field-by-field negative regression for that contract.

## Validated identity

The test first requires a freshly built plan to match.  It then changes exactly
one member at a time and requires cache rejection for all of the following:

- integer topology: `nao`, `nsh`, `maxsh`, `nsh_id`, `nao_sh`, `ish_at`, and
  `iao_sh`;
- scalar model parameters: `frscale`, `omega`, `lrscale`, `ondiag_scale`,
  `hubbard_exp`, `hubbard_exp_r0`, `gexp`, and `corr_exp`;
- array model parameters: `hubbard`, `onecxints`, `offdiag_scale`, `rad`,
  `kq`, `corr_scale`, and `corr_rad`.

Every real field is restored from the byte-identical cached signature before
the next probe.  The test then copies the live BvK kernel, swaps representative
columns 1 and 2 without changing the representative set, and requires
`bvk_plan_matches` to reject the altered order.  This distinguishes exact plan
identity from a weaker unordered-set comparison.

## Results

| Platform / build | `tblite/exchange` | `tblite/gxtb` | Result |
|---|---:|---:|---|
| macOS Debug | 3.62 s | 14.47 s | 2/2 passed; 0 failed |
| macOS Release | 1.76 s | 6.21 s | 2/2 passed; 0 failed |
| Terok Linux Release | 2.61 s | 12.40 s | 2/2 passed; 0 failed |

The exchange suite contains 36 subtests and the g-xTB suite 44 subtests.  The
largest reported BvK-supercell residual was
`4.6134320322299693e-13`, within the established regression tolerance.  The
Linux test executable SHA-256 was
`4c4cc95530e6db966118e5abea5dff8870e84ce213cb6fbd1dd5cf5f43a539b1`.

The Linux build and tests ran under the CPU-141 reservation lock with
`Cpus_allowed_list: 141` and all OMP/BLAS thread limits set to one.  Before the
launch, `MemAvailable` was 518074280 KiB.  One healthy CP2K calculation was
live; its remaining-growth allowance was 103957604 KiB.  After the 16-GiB
candidate allowance, the computed margin was 397339460 KiB, above the required
134217728 KiB.  The unrelated healthy CP2K process on CPU 76 was not touched.

## Contents

- `local/`: raw macOS Debug and Release CTest logs and executable identities.
- `linux/`: raw Terok build/CTest logs, complete live-RSS inventories, memory
  arithmetic, exact affinity proof, source hashes, source state/history, and
  platform provenance.
- `source/`: the signed provider patch, source history, and exact launch and
  evidence-collection scripts.
- `verify_archive.sh`: integrity and acceptance checks.
- `SHA256SUMS`: portable SHA-256 manifest over the complete archive.

This is a correctness qualification, not a timing or whole-process memory
benchmark.  It closes exact model-signature and active-plan-order coverage;
broader repeated performance matrices remain a separate gate.
