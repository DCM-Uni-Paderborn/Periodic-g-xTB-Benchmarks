# save_tblite pre-kernel partial k-to-R accumulator ABI

Date: 2026-07-16

Source worktree: `/private/tmp/save_tblite_partial_k_to_r_20260716`

Branch: `codex/partial-k-to-r-accumulator`

Base revision: `257ba442684c39454175e5192c8a2342b4c6380f`

Qualified commit: `e57095e15e88a3cf5cecb70cf1379f68ff05867c`

## Contract

The public API is MPI-free.  Every host partition pushes only its uniquely
owned Bloch blocks into a bounded stream, exports a contiguous *linear
pre-kernel* source payload and an integer ownership vector, and may then be
discarded.  The host sums those payloads and ownership vectors.  An importer
proves that every k point has global owner count one and only then invokes the
unchanged nonlinear exchange kernel.  The strict designated-root mode uses one
importer.  The optional replicated-kernel mode imports the same exact sources
on several ranks and supplies one common global `owner_rank(k)` map: each
Bloch Fock/overlap-adjoint block is pulled on exactly one importer, including
the valid empty-importer case when `P > Nk`.  Exactly one declared result owner
may retrieve replicated scalar energy/shell-potential or direct force/stress
results.  Finalized results are never summed across partitions.

Every logical forward or reverse transaction requires a fresh nonzero
host-generated 128-bit identity, broadcast unchanged to all contributors.
Descriptors expose a canonical padding-free 30-word `int64` encoding for host
allgather.  The host must unpack and compare every field exactly before payload
reduction/import; the internal two-word buffer digest is not an acceptance
criterion.  Descriptor and raw payload are one atomic transport record:
rebinding saved raw values to a fresh descriptor is explicitly forbidden.

Forward and reverse transactions have separate descriptors, generations,
payload types, and state-machine flags.  A descriptor binds the ABI and scalar
width, phase, generation, `Nk`, `nAO`, `nspin`, shell count, regular mesh,
current image batch and capacity, image range, transaction identity, k coordinates and ordering,
weights, representatives and input mapping, compact Fourier plan and twist,
exchange-kernel/model state, charge-dependent onsite state, geometry, lattice,
species, and periodicity.  Imports rebuild this identity from the live stream
and also rerun the calculator/cache/geometry state check.

## Exact payload storage

For current image batch width `B`, the exported forward payload contains

```
N_forward = 4 * nAO^2 * B * nspin
```

complex scalars: `A_R`, partner `A_R`, `C_R`, and `V_R`.

The exported reverse payload contains

```
N_reverse = 6 * nAO^2 * B * nspin + 5 * nAO * nspin
```

complex scalars: the six batched AO-image sources (`A_R`, `C_R`, `V_R`,
partner `A_R`, reverse `A_R`, reverse `V_R`) and the five onsite vectors
(`gdiagP`, `gdiagSP`, `gdiagSPS`, reverse `bdiagSP`, reverse `bdiagSPS`).
With binary64 complex scalars these payload byte counts are exactly
`16*N_forward` and `16*N_reverse`.  Each export separately carries `Nk`
default-integer ownership counts and one fixed-size opaque descriptor.  Unit
tests query the public size functions and assert both complex-scalar formulas.
All extent arithmetic and slicing use checked `int64` values.  The provider
fails before allocation if a payload would exceed `huge(default integer)`;
the CP2K host can independently chunk MPI reductions below that indexing cap.

## Qualification matrix

- `P = 1, 2, 4` on an RKS `3x1x1` mesh; `P=4 > Nk=3` deliberately contains
  one empty contributor.
- A separate `P=4`, `B=1` run traverses all three image batches, imports on all
  four simulated ranks, distributes forward and reverse pulls with one exact
  global owner map, and deliberately leaves rank four empty (`P=4 > Nk=3`).
  Only rank one retrieves scalar/direct results.  An additional empty exported
  contributor exercises the public contributor-advance operation in every
  generation.  The test checks exact bounded payload formulas and rejects stale
  forward and reverse descriptors from generation one in generation two.
- Permuted k-input order with a common mesh twist on `3x1x1` (`P=4`).
- Shifted/noncanonical lower-dimensional `9x9x1` UKS mesh (`P=4`).
- Shifted three-dimensional `2x2x2` mesh (`P=2`).
- Forward energy, Fock response, and onsite response compared with the dense
  complete-mesh oracle.
- Reverse overlap adjoint, analytic force, and analytic stress compared with
  the dense complete-mesh oracle.
- Negative cases cover duplicate local push, out-of-range k index,
  apply-before-import, missing/duplicate/negative ownership, descriptor/image
  range mismatch, cross-transaction payload replay, a divergent descriptor on
  an empty rank, descriptor pack/unpack corruption, oversized payload arithmetic,
  non-finite payload, second import, push-after-import, discard-before-export,
  attempts to discard an imported accumulator, invalid/duplicate pull-owner
  setup, missing owned pulls, unowned pulls, and non-owner scalar/direct-result
  retrieval.  Unconfigured designated-root cases in the Debug build also guard
  against non-short-circuit access to an absent distributed-pull mask.

## Local results (macOS arm64)

Compilers/tools: GNU Fortran 16.1.0, CMake 4.3.4, Ninja 1.13.2.

- Debug build with `-O0 -g -fcheck=all -fbacktrace` and floating-point traps:
  `tblite-tester exchange`: 31/31 passed.
- Release build: `tblite-tester exchange`: 31/31 passed.
- Debug and Release `ctest -R exchange`: 1/1 passed each.
- Maximum partial-ABI residual over all forward/reverse cases:
  `9.7144514654701197E-017`.
- Representative reverse residuals (response, force, stress):
  - RKS `3x1x1`: `1.3877787807814457E-017`, `0`, `0` (Release).
  - UKS batched `9x9x1`: `5.5511151231257827E-017`, `0`,
    `1.3552527156068805E-020` (Release).
  - shifted 3D `2x2x2`: `2.7755575615628914E-017`, `0`,
    `2.0679515313825692E-025` (Release).

The broader local Release CTest run passes the exchange target and all tests
relevant to this ABI.  A few unrelated finite-difference threshold/platform
tests and the DDX-disabled xtbml case remain baseline/environment failures and
are recorded separately; they do not exercise the partial accumulator.

## Linux reproduction on terok

Isolated source:
`/home/kuehne88/work/codex-save-tblite-partial-k-to-r-20260716`

Release build:
`/home/kuehne88/work/codex-save-tblite-partial-k-to-r-release-20260716`

Debug build:
`/home/kuehne88/work/codex-save-tblite-partial-k-to-r-debug-20260716`

GNU Fortran 15.2.0 and the existing conda/OpenBLAS dependency environment are
used.  OpenMP linkage is enabled because the preinstalled dependency archives
were built with OpenMP; the partial accumulator API itself contains neither
MPI nor OpenMP collectives and leaves all reductions to the host.

All controlled qualification commands set `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `GOTO_NUM_THREADS=1` to
avoid host-default thread oversubscription.  Earlier unrestricted builds were
left running and were not used as qualification evidence.

- Linux Release `tblite-tester exchange`: 31/31 passed; maximum partial-ABI
  residual `1.2490009027033011E-016`.
- Linux Debug with `-O0 -g -fcheck=all -fbacktrace` and floating-point traps:
  `tblite-tester exchange`: 31/31 passed; maximum partial-ABI residual
  `9.7144514654701197E-017`.
- MPI-freedom audit of the Release archive/executable/source reports zero
  undefined MPI symbols, zero linked MPI libraries, and zero MPI source
  imports.  This is intentional: the CP2K host owns the collective SUM.
- A two-rank `mpirun -np 2` Release smoke launches two independent copies of
  the complete exchange qualification: both reach 31/31 passed, neither emits
  a failure, and both report `1.2490009027033011E-016`.  The simulated
  `P=1/2/4` host reductions inside each qualification are the numerical test
  of the MPI collective contract; no provider-internal collective is used.
