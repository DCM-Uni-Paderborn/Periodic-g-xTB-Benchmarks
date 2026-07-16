# save_tblite pre-kernel partial k-to-R accumulator ABI: revised qualification

Date: 2026-07-16

Branch: `codex/partial-k-to-r-accumulator`

Base revision: `257ba442684c39454175e5192c8a2342b4c6380f`

Qualified commit: `35e7942b60edd89bb407ab3da5768d3410af83f5`

The qualified commit is local and intentionally unpushed.  This bundle is a
new revision; the evidence bundle for `e57095e15e88a3cf5cecb70cf1379f68ff05867c`
is retained unchanged as historical evidence.

## Contract under review

The provider API is MPI-free.  Each host partition pushes only its uniquely
owned Bloch blocks and exports an additive, pre-kernel source payload plus a
local ownership vector.  CP2K sums payloads and ownership vectors, proves that
every k point has exactly one global owner, and imports the result before one
nonlinear exchange apply.  Replicated importers may partition Fock or overlap-
adjoint pulls with one common owner map; exactly one declared rank retrieves
the replicated shell-potential/energy or force/stress result.  Empty ranks are
valid.

Each forward or reverse transaction uses a fresh, nonzero host-generated
128-bit identity.  Its canonical 30-word int64 descriptor binds ABI/scalar
width, phase, generation, mesh and dimensions, current batch and image range,
k ordering and weights, geometry and periodicity, compact Fourier data,
exchange model/kernel state, and charge-dependent onsite state.  Imports
rebuild the descriptor from live state and compare every field exactly.

## Adversarial hardening added during review

The review found one reproducible default-integer overflow alias in the public
mesh gate: on 32-bit default integers, `PRODUCT([65536,65537,1])` can wrap to
65536.  The three public whole-mesh/stream/gradient gates now use checked int64
arithmetic and reject this exact case before plan construction.  Since the
accepted product is then at most `Nk <= huge(default integer)`, the remaining
internal lexicographic mesh index products are bounded.

Two modular phase products were independently widened to int64.  Image-range
validation no longer relies on Fortran logical short-circuit behavior, and
current-batch end/advance arithmetic avoids the `huge(default integer)+1`
intermediate at the final image.  Tests cover the overflow-alias mesh and an
extreme invalid image start.

## Exact source and Linux reproduction

`source-worktree.tar.gz` is `git archive HEAD` for the qualified commit.  Its
SHA-256 is recorded in `SHA256SUMS`.  The uploaded archive had the same hash on
Terok before extraction.  Hashes of all nine changed implementation/test files
matched locally and after extraction.

Terok paths:

- source: `/home/kuehne88/work/codex-save-tblite-partial-k-to-r-35e7942`
- Debug: `/home/kuehne88/work/codex-save-tblite-partial-k-to-r-debug-35e7942`
- Release: `/home/kuehne88/work/codex-save-tblite-partial-k-to-r-release-35e7942`

Linux used GNU Fortran 15.2.0 and the existing conda/OpenBLAS dependency
environment.  All qualification runs set `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `GOTO_NUM_THREADS=1`.

## Results

- macOS Debug exchange suite: 31/31 passed.
- macOS Release exchange suite: 31/31 passed.
- macOS Debug and Release CTest exchange target: 1/1 passed each.
- Terok Linux Debug exchange suite: 31/31 passed; maximum partial residual
  `9.7144514654701197E-017`.
- Terok Linux Release exchange suite: 31/31 passed; maximum partial residual
  `1.2490009027033011E-016`.
- Terok two-rank Release smoke: 62 passed test cases, zero failures, and two
  complete residual reports of `1.2490009027033011E-016`.
- MPI-freedom audit: zero undefined MPI symbols, zero linked MPI libraries,
  and zero MPI source imports.

The qualification covers energy and Fock response, shell/onsite response,
overlap adjoints, analytic force and stress, one/two/four simulated host
partitions, an empty contributor (`P > Nk`), multi-batch generation changes,
permuted and shifted meshes, RKS and UKS, and 1D/2D/3D regular grids.  Negative
tests cover stale/cross-transaction descriptors, descriptor corruption,
missing/duplicate/negative owners, invalid pull ownership, missing/unowned
pulls, illegal state transitions, non-finite payloads, and oversized extents.

The remaining requirements are explicit host responsibilities, not hidden
provider guarantees: globally fresh transaction identities, exact descriptor
consensus before reduction, one common pull-owner map/result owner, and atomic
descriptor/payload transport.  The storage-bit fingerprints are deterministic
state-change guards, not cryptographic authentication.

