# Mixed-radix FFT exchange-transform validation

Verdict: **PASS for numerical equivalence to the dense transform oracle in
the archived scope**.  This archive qualifies the selectable separable and
mixed-radix FFT implementations of the regular-mesh exchange transform.  It
does not by itself establish a speedup claim.

## Qualified scope

Every candidate is compared with the dense transform using the complete
printed energy sequence and the final energy.  Printed forces and stress are
also compared where the input requests them.  The matrix covers:

- RKS and UKS;
- K290, time-reversal, and SPGLIB k-point handling;
- shifted and unshifted regular meshes;
- one-, two-, and three-dimensional periodicity; and
- energy, force, and stress paths.

The verifier requires `PROGRAM ENDED` in every output.  Its acceptance limits
are `1e-10` hartree for energy, `1e-8` hartree/bohr for force, and `1e-3` bar
for stress.  The largest archived differences are:

| Observable | Maximum absolute difference |
|---|---:|
| printed energy sequence | `1.42108547152020037e-14` hartree |
| final energy | `7.10542735760100186e-15` hartree |
| force | `1.74607420800000003e-16` hartree/bohr |
| stress | `5.00000169267877936e-7` bar |

`summary.tsv` contains all case-wise residuals.  Both the separable and FFT
paths pass against the same dense oracle.

## Source and build identity

- CP2K implementation commit: `cc782e71dce6d7f4808a2e0db209c41adbecf0ed`
- CP2K implementation patch:
  `source/0001-Select-mixed-radix-FFT-for-g-xTB-exchange.patch.gz`
- save_tblite implementation commit:
  `bef49671d1a238bc9f4098565f4d22b77050b657`
- save_tblite implementation patch:
  `source/0001-Accelerate-regular-exchange-transforms-with-mixed-ra.patch.gz`
- tested CP2K executable SHA-256:
  `759b85c7aacf4d49617cbbfb42b63b8d0b10a8b49c5e1b479689ab67b1aead7b`
- linked static `libtblite.a` SHA-256:
  `5b623106986393b91f965405da651933960f747968516925a548a298e2983bca`

The executable was built from the exact implementation worktree before the
selector change was committed, so its printed CP2K revision is the parent
`09ac4390801888c7c16b780ca68da24622665728`; the two archived patches freeze
the tested source changes.  `CP2K_GXTB_EXCHANGE_TRANSFORM_MODE` selected the
`dense`, `separable`, or mixed-radix FFT path for each run.

## Reproduction and integrity

The short qualification runs used one MPI rank and one OpenMP thread on
`MacBook-Pro-von-Thomas-3.local`.  Absolute build paths and original output
locations are retained in `provenance.txt`; they are provenance, not required
locations for verification.

From this directory run:

```text
python3 verify.py
shasum -a 256 -c SHA256SUMS
```

The first command must reproduce `summary.tsv` and exit zero.  `SHA256SUMS`
authenticates every archived file except the manifest itself.
