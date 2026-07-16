# Independently reviewed CP2K cross-mesh k-point restart

Verdict: **PASS for the opt-in `VALIDATED_BVK_TRANSFER` initial-guess
path**.  This is not the default restart mode and it does not transfer a
Hamiltonian, Fock matrix, DIIS state, or tblite/g-xTB mixer history.

## Frozen source identity

- CP2K branch: `codex/gxtb-crossmesh-restart`
- reviewed base: `28df9380abb327d56bbf216d2469a1fd8c953fc0`
- published implementation: `8f10e4bed387ae3a430fc9fab0be158970b0f20b`
- reviewed staged-diff SHA-256:
  `d933b478a69ca0db5566d87a4d5b132455ec263081da5cf315de70e469c1add4`
- reviewed debug-binary SHA-256:
  `c55fd95cc3a28e33cca3a01b5e27aa1baa7fa72973b49dc50a2adcbc788e86f8`
- reproducible source patch:
  `source/0001-Add-validated-cross-mesh-k-point-restarts.patch.gz`

## Qualified behavior

The independently reviewed build, the focused official regression, five
malformed-file fallbacks, legacy/default-mode compatibility, RKS and UKS
transfers, force/stress evaluation, final-density comparison, and a two-rank
Linux MPI run all passed.  Rejected files fall back transactionally to a cold
initial guess.  In particular, NaN metadata, non-finite or checksum-overflowing
payloads, and trailing data are rejected without aborting CP2K.

The principal numerical residuals against independently converged cold target
runs are:

| Test | SCF steps, transfer / cold | Final energy difference / Ha | Maximum final residual |
|---|---:|---:|---:|
| RKS, `3x3x3` to `4x4x4` | 7 / 15 | `8.881784197001252e-16` | force `9.999999717180685e-10` Ha/bohr; stress `1.499999780207872e-05` bar |
| RKS final density | -- | -- | max `1.738517871330281e-10`; relative Frobenius `1.177179419484646e-10` |
| UKS, `3x3x3` to `4x4x4` | 6 / 11 | `5.49376e-10` | final-density max `8.456397040945696e-12`; relative Frobenius `1.949211790267738e-11` |
| MPI-2, `5x5x5` to `6x6x6` | 2 / 10 (source) | target energy `-2.951453440027636` | accepted transfer; both runs ended normally |

The UKS energy entry compares the separately converged transferred and cold
runs reported in `EVIDENCE.md`; their density residual is the more direct
same-target diagnostic.  The complete precision, inputs, raw outputs, scripts,
and provenance are retained in this directory.

## Explicit limitation

Fourier interpolation followed by scalar electron-number normalization does
**not** guarantee an N-representable or positive-semidefinite metric density.
The measured occupation ranges are
`[-1.1093841e-5, 2.0009592]` for the RKS test and
`[-4.8857089e-7, 1.0000507]` for the UKS test.  The tested SCFs healed those
small initial-guess excursions and converged to the cold references, but the
path must not be described as positivity preserving or unconditionally safe.
Before broad or default activation, add either an electron-conserving spectral
projection in the overlap metric or a conservative spectral accept/fallback
gate.

`EVIDENCE.md` is the complete review report.  `SHA256SUMS` authenticates every
archived file except the manifest itself; verify it from this directory with
`shasum -a 256 -c SHA256SUMS`.
