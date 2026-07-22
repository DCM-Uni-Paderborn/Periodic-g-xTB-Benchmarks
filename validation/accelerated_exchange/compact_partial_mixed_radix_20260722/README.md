# Range-local partial exchange with compact mixed-radix transforms

Date: 2026-07-22

Verdict: **PASS for the archived correctness and targeted-allocation scope**.
This bundle qualifies the direct composition of the range-local partial-image
exchange route with compact mixed-radix AO-pair transforms, streamed
post-mixer symmetry-star contraction, and streamed analytical reverse
contraction. It is method-development evidence for Part II; it is not an
application benchmark and makes no wall-time or whole-process RSS claim.

## Frozen implementation

- `save_tblite` commit:
  `62e28eb7b89af51dd497fb780bc24cd19679ce14`
- CP2K commit:
  `bd5b7aed038c3964d6dd58f2b742d8acd8d9d784`
- Terok CP2K executable SHA-256:
  `df6466552a495e94e710174dbf468ec1765f1a4690ef955f9129ef983f35790b`
- Terok linked `libtblite.a` SHA-256:
  `d55a9cdfd52c80e955b6dcb8854e993b8f26b0c629a3ed4e0028547d9ed44976`
- macOS CP2K executable SHA-256:
  `180fd8eb055d46547cf2936f7df197b80b01681e8275d3d8186b466d6fad64d2`
- macOS linked `libtblite.a` SHA-256:
  `85389e07c58fa1e390a975a9dfd687e0d04fc24cc1aa4df42525993a1a418d88`

The exact committed patches are under `source/provider/` and `source/cp2k/`.
The precommit patches are retained separately because they bind the tested
Terok source hashes before the commit objects were created. Committing did not
change the source bytes.

## Main end-to-end gate

The CH4 K290-reduced `2 x 2 x 2` calculation uses two CP2K k-point groups,
image batch two, `K_GROUP_PARTIAL_DISTRIBUTED_IMAGES`, `MIXED_RADIX_FFT`,
streamed post-mixer contraction, compact ACP response, and streamed analytical
exchange reverse in the same SCC sequence. The production calculation and the
qualification calculation both terminate normally at
`-40.473748967057020` hartree on macOS and Terok Linux. Every forward and
reverse transform marker reports `fallback=0`.

| Host | Max forward dE | Max forward dVsh | Max forward dFfold | Reverse dOverlap | Reverse dForce | Reverse dStress |
|---|---:|---:|---:|---:|---:|---:|
| macOS | 5.5511e-16 | 3.4694e-18 | 7.2165e-16 | 7.4940e-16 | 3.2526e-19 | 1.7347e-17 |
| Terok Linux | 8.8818e-16 | 3.4694e-18 | 8.3267e-16 | 6.3838e-16 | 3.7947e-19 | 1.5613e-17 |

The streamed post-mixer temporary contains 192 complex elements versus 512
for its dense full-star oracle. The independent ACP sparse reverse residual is
`2.5685e-8`, below its architecture-portable `1e-7` gate.

## Breadth matrix

The compact partial transform and analytical reverse were also exercised for
true one- and two-dimensional periodicity, RKS and UKS, time reversal, a
shifted explicit full grid, and shifted SPGLIB reduction. Shared macOS/Linux
cases have identical final energies at the printed precision.

| Case | Platform coverage | Final energy / hartree | Max forward residual | Max reverse residual |
|---|---|---:|---:|---:|
| H2, 1D TR `3 x 1 x 1` | macOS + Linux | -1.154127291280940 | 5.5511e-17 | 2.7759e-17 |
| H2, 2D TR `3 x 1 x 3` | macOS + Linux | -1.137599991534640 | 8.3267e-17 | 5.5512e-17 |
| O2, UKS TR `3 x 1 x 1` | macOS + Linux | -150.541892482180998 | 8.8818e-16 | 1.7365e-16 |
| Si, shifted full `2 x 2 x 2` | macOS + Linux | -579.050928767344203 | 9.9920e-16 | 8.3267e-17 |
| CH4, shifted SPGLIB `4 x 4 x 4` | Linux | -40.468551982577388 | 3.2196e-15 | 1.5543e-15 |

All nine positive breadth runs terminate normally and report the exact
mixed-radix forward/reverse markers with zero fallback. `matrix_summary.tsv`
is generated directly from the raw outputs and retains each observable
separately.

## Provider and affinity qualification

The Terok provider Release suites pass 36/36 exchange tests and 44/44 g-XTB
tests. The compact partial-transform provider residual is at most
`1.3878e-16` in the shifted `9 x 9 x 1` UKS case. Its targeted compact/full
high-water counters are 3700/12144 complex elements for that case and 236/368
for the shifted `2 x 2 x 2` case.

Every CP2K run uses two disjoint singleton CPU reservations, ranks 0 and 1 on
logical CPUs 141 and 142. The launcher sets OMP, OpenBLAS, MKL, BLIS, and
VECLIB thread counts to one, holds one lock per CPU, and records the pre-exec
`/proc` mask. `linux/results/*/provenance/` stores the exact binary, provider,
input, and source hashes as well as the complete prelaunch memory/RSS
inventory.

## Preserved fail-closed and independent findings

Negative records are retained but excluded from positive statistics:

1. `local/pre_guard_fix` preserves the original combination guard before the
   compact partial route was enabled.
2. `matrix/local_missing_spglib` is the expected fail-closed result from the
   local binary that was built without SPGLIB; Linux provides the qualified
   SPGLIB result.
3. `matrix/local_pre_fullgrid_guard` shows that streamed symmetry-star
   contraction correctly rejects an unreduced full grid. The accepted full
   grid case uses the dense star selector while retaining compact partial
   forward/reverse transforms.
4. `matrix/local_pre_acp_gate` records an independent H2 ACP sparse-gradient
   mismatch of about `3.9903e-5`, above the ACP-only `1e-7` gate. The compact
   exchange forward/reverse checks in that run had already passed. This ACP
   finding remains open and is not mixed into the compact-transform
   acceptance claim.

## Verification

From this directory run:

```sh
./verify_archive.sh
```

`verify.py` reconstructs both TSV summaries from raw outputs, checks all
normal-termination and exact-selector markers, applies the `1e-12` numerical
gates, checks the main ACP and mixer counters, verifies cross-platform
energies, and requires every preserved fail-closed record. `SHA256SUMS`
authenticates every curated file except the manifest itself.
