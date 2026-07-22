# Bounded ACP mesh contraction qualification

This archive qualifies the independently selectable periodic atom-centered-potential (ACP) mesh contraction added to the Part-II acceleration tree. The physical ACP operator and derivatives are unchanged.

## Implemented boundary

- `DENSE` retains the complete `ACP(AO,AO,Nk_full)` Bloch operator and the established dense derivative oracle.
- `STREAMED` evaluates at most the negotiated full-mesh k-point batch, folds each block immediately, and never allocates the complete ACP Bloch tensor.
- The streamed reverse path contracts compact projector images directly with CP2K's existing sparse real-space density images. It forms neither the quadratic projector-translation difference set nor a zero-padded density tensor.
- `QUALIFY` evaluates both representations, gates the compact provider response against the complete image formula at relative tolerance `1e-12`, gates the final coordination/charge/Cartesian/strain response at `1e-7`, and applies the streamed result once.
- Exact integer primitive-cell labels are used for projector and density translations. Missing density translations are exact zero blocks.

The scoped forward workspace changes from `NAO^2 * Nk_full` to `NAO^2 * B` complex values, apart from the irreducible response already required by CP2K. The reverse path retains compact projector images, their response, and existing sparse density images. This is not a claim that total process RSS is independent of the k-point mesh.

## Acceptance gates

Every accepted run must terminate normally and use the archived executable/provider hashes. Dense, streamed, and qualify modes are compared for total energy, force components, and stress components when printed. The open-shell O2 control additionally compares the alpha/beta populations. An unknown environment selector must fail before an SCC response is exposed.

The provider `tblite/acp` and `tblite/gxtb` tests pass in both local Debug and Release builds. They cover exactly-once streamed transactions, incomplete/duplicate/nonfinite pushes, geometry invalidation, compact forward batches, sparse density support, dense-response equality, and optional-output safety.

The archived source patches are byte-identical to provider commit
`294099dab74fc0a49e58f42e789ed10b961907ca` and CP2K commit
`448a5e4d0316775a272f941a7f54434a5c9d6ce5`; both commits carry the required
developer sign-off.  Their parent histories are included under `source/`.

## Same-build matrix

All dense, streamed, and qualify calculations terminated normally.  The macOS
matrix used CP2K executable SHA-256
`9a35911c4e7d2a95e120310962c4184659d95388fdd02ab4588ca5fa68dfb967`
and provider archive SHA-256
`7c82e18b5bd14409849f3766c3d6e0d686478674c825d8ae7a292613aa693ca4`.
The independent Linux Release matrix used CP2K executable SHA-256
`18f29832f1044e918dff7c02f07b3e2862db322e2735332362c7c6ddf82e3335`,
linked CP2K library SHA-256
`384149b9a7011bbb1a47d9fa7b13602ea9910fb28d046bc7f418214575c1b201`,
and provider archive SHA-256
`d627125671914f59628a24b8b725e62f43f1175eba9b264a9d4182c4ff64197f`.

| System | Platform / mesh / route | Energy (Ha) | max abs streamed-dense force (Ha/bohr) | max abs streamed-dense stress (bar) | Sparse counts / gate |
|---|---|---:|---:|---:|---|
| CH4 | macOS K290 2x2x2, complex | -40.473748967057020 | 3.6000e-9 | 1.5109e-3 | 23 projector, 443 density; provider 0; final 2.5685e-8 <= 1e-7 |
| CH4 | Linux K290 2x2x2, complex | -40.473748967057020 | 2.0000e-9 | 1.5108e-3 | 23 projector, 443 density; provider 0; final 2.5685e-8 <= 1e-7 |
| CH4 | macOS K290 3x3x3, exact symmetry-fused combination | -40.468866070692428 | 1.8000e-9 | 7.5070e-4 | 27 full points in batches of 2; provider 0; final 1.2767e-8 <= 1e-7 |
| CH4 | Linux K290 3x3x3, exact symmetry-fused combination | -40.468866070692428 | 1.0000e-9 | 7.5070e-4 | 27 full points in batches of 2; provider 0; final 1.2767e-8 <= 1e-7 |
| Si | macOS shifted full 2x2x2 | -579.050928767344203 | 0 | 0 | 231 projector, 3025 density; provider 0; final 0 |
| Si | Linux shifted full 2x2x2 | -579.050928767344203 | 0 | 0 | 231 projector, 3025 density; provider 0; final 0 |
| O2 UKS | macOS and Linux 3x1x1 | -150.541892482180998 | not printed | not printed | no ACP channel; alpha/beta populations 7/5 in every mode |

For CH4, both platforms used batch 2 with `fullStorage=0` on the 8- and
27-point physical meshes; for Si, batch 8 with `fullStorage=0`.  The CH4 final
residual is roundoff from separately allocated small-matrix derivative work.
The provider response itself is exactly zero against its stricter `1e-12`
gate.  The invalid selector controls exited nonzero with
`Unknown CP2K_GXTB_ACP_MESH_CONTRACTION value`.

The exact symmetry-fused K290 combination is also permanent regression input,
not only an environment-override campaign case: its production and
qualification files select `STREAMED` and `QUALIFY` directly, respectively,
and dedicated matchers require both the bounded-forward/sparse-reverse markers
and the two response-oracle markers.

## Linux release qualification

The Terok release build and exact-singleton qualification are archived under
`linux/`.  The build and every matrix case ran under the CPU-141 reservation
lock with `Cpus_allowed_list: 141` and all OMP/BLAS thread limits set to one.
Before the final launch and before every case, the inventory records every live
RSS value, available memory, remaining-growth allowance, candidate peak, and
the 128-GiB safety margin.  The minimum computed margin was 306127200 KiB.
The unrelated healthy CP2K calculation on CPU 76 was left untouched.

## Contents

- `inputs/`: official CP2K regression inputs.
- `local/`: raw macOS outputs, exit code, and binary hashes.
- `linux/`: raw Terok build/test/run evidence and exact affinity/memory provenance.
- `source/`: provider and CP2K source patches, build/final launchers, and
  parent history.
- `summary.tsv`: machine-readable paired numerical results.
- `verify_archive.sh`: checksum, normal-end, storage-marker, response-gate,
  and invalid-selector verification.
- `SHA256SUMS`: portable integrity manifest over the complete archive.

The O2 calculation contains no ACP parameter channel and is therefore a selector-neutrality/open-shell regression, not evidence for the sparse ACP reverse contraction. The CH4 and Si rows are the ACP-specific energy/force/stress qualification.
