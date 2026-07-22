# Complex-density ACP reverse-orientation qualification

Date: 2026-07-22

Verdict: **PASS**.  This bundle closes the independent H2 sparse-ACP
derivative finding retained by the compact partial-transform campaign.  The
accepted result is a correction to the AO/image orientation of the CP2K
real-space density in the compact ACP reverse contraction; it does not change
the ACP energy or introduce an approximation.

## Symptom and independent diagnosis

The original H2 time-reversal `3 x 1 x 1` qualification aborted as intended:
the provider-level response oracle was zero, but the final sparse reverse
residual was `3.9903e-5`, above its `1e-7` gate.  The channel maxima were
`2.1848e-5` for `dEdCN`, `4.2713e-6` for `dEdq`, `4.5013e-6` for the
Cartesian gradient, and `3.9903e-5` for strain.  The exact energy was unchanged
at `-1.154127291280940` hartree, whereas the first-atom x force was
`5.05908009e-2` hartree/bohr in the dense route and `5.05964691e-2`
hartree/bohr in the pre-fix streamed route.

Central energy differences independently select the dense force:

| displacement / angstrom | central-difference force / hartree bohr^-1 |
|---:|---:|
| 1e-3 | 5.0590838571e-2 |
| 1e-4 | 5.0590801516e-2 |

The `1e-4` result differs from the dense analytical force by about `6.2e-10`
hartree/bohr, while the old streamed force differs by `5.6682e-6`
hartree/bohr.  This finite-difference result is independent of the dense
derivative oracle.

## Corrected orientation

Let `C_V` be a compact projector image, `L` the diagonal ACP channel strength,
and `P_R` the real-space density gathered with CP2K's k-to-image phase
convention.  The response of projector image `U` is

```
D_U = sum_V [ L C_V P_(U-V) + L C_V P_(V-U)^T ].
```

The pre-fix implementation attached the transpose to the opposite translation.
That is invisible for Gamma-only or real symmetric image blocks but is wrong
for a complex time-reversal mesh.  The implementation now uses the expression
above in both the production contraction and its internal complete-image
oracle.  The provider unit test retains its general image-density convention
and separately constructs the CP2K gathered density, so the two conventions
cannot be conflated again.

## Same-build qualification

All four corrected qualification calculations terminated normally.  The
provider compact-response residual is exactly zero in every row.  `Final`,
`gradient`, and `sigma` are the absolute channel maxima printed by the
independently allocated dense/sparse reverse comparison; the gate is `1e-7`.

| Host | Mesh | Energy / hartree | Final | gradient | sigma |
|---|---|---:|---:|---:|---:|
| macOS arm64 | 1D TR `3 x 1 x 1` | -1.154127291280940 | 5.5511e-17 | 2.7756e-17 | 0 |
| macOS arm64 | 2D TR `3 x 1 x 3` | -1.137599991534640 | 6.6613e-16 | 5.5511e-17 | 5.5511e-17 |
| Terok Linux | 1D TR `3 x 1 x 1` | -1.154127291280940 | 1.1102e-16 | 2.7756e-17 | 1.1102e-16 |
| Terok Linux | 2D TR `3 x 1 x 3` | -1.137599991534640 | 5.5511e-16 | 5.5511e-17 | 8.3267e-17 |

The macOS qualification binary is SHA-256
`13691b6015ff4a4e1573df1a05cfb6fac6ebf94efba2b486e84f2210227450dd`;
its linked provider archive is
`ff940532a44c5392e4b6bc52aced52aefb8355b88e16020bf6db05f6a240b248`.
The accepted Terok binary is
`df6466552a495e94e710174dbf468ec1765f1a4690ef955f9129ef983f35790b`;
its provider archive is
`fe210c64a4c4fa6897668a8657dd234046143c55bda2f7c9279d24108f2f152a`.
Every Linux case records these identities before launch.

The Linux Release provider selection passed `tblite/exchange` and
`tblite/gxtb` (2/2).  The CP2K 1D and 2D regression inputs now select
`ACP_MESH_CONTRACTION QUALIFY`, and their matchers require both the provider
response and final sparse-reverse qualifications.

## Source and execution provenance

- save_tblite commits:
  `331757ad2d4ac814521911a620bc724ad753bc57` and
  `5c25e149ad10f85bb12f31b0f50860cd3589a1ce`
- CP2K commit: `e01b895346411dc17e1df865ac0384dd48a0efdb`
- all three commits carry
  `Signed-off-by: Thomas D. Kühne <tkuehne@cp2k.org>`

Before transfer, every Terok source file was checked byte-for-byte against the
expected parent commit and backed up under `linux/provenance/`.  The replacement
files were then checked against the published branch versions.  Provider and
CP2K builds ran on reserved singleton CPU 90 with OMP/BLAS thread counts equal
to one.  Each CP2K run held disjoint reservation locks for CPUs 90 and 91 and
recorded pre-exec `/proc` masks for both MPI ranks.  `MemAvailable` was about
494 GiB before each launch.  The unrelated healthy CP2K calculation on CPU 76
was not touched.

## Contents

- `local/pre_gate/`: the preserved gate failure.
- `local/diagnostic/`: dense, old streamed, corrected, two-dimensional, and
  four finite-difference raw calculations.
- `linux/cases/`: accepted 1D/2D inputs, outputs, hashes, complete prelaunch
  RSS inventories, and rank-affinity proofs.
- `linux/build/`: provider tests and provider/CP2K build logs.
- `linux/provenance/`: pre-update sources, source hashes, build-affinity proofs,
  and qualified binary identities.
- `source/`: exact signed source patches, snapshots, regression inputs, and
  commit histories.
- `summary.tsv`: machine-readable numerical summary.
- `verify_archive.sh` and `SHA256SUMS`: semantic and byte-integrity checks.

