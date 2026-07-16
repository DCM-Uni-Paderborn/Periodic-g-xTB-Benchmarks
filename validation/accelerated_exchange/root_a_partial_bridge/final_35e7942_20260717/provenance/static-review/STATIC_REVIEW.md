# Static review: CP2K `KGROUP_PARTIAL_ROOT`

## Scope and result

The reviewed CP2K change is the sole diff in `src/tblite_interface.F` relative
to base `0a1f7e3329a3e6c2a6accff28617af53fb9943b4`.  Reviewed source SHA-256 is
`47c9b039b2e0d081f1ac3688f29f5c75ffed9a60acbb490f7bee1ae99593dd5d`;
the patience patch SHA-256 is
`8763981ef9c6ba7e9db26a11f4245e11fbcc52ac393a6b6421b09e8e7628156f`.

Result: **PASS, no remaining finding in the reviewed bridge diff.**  This
conclusion is for correctness, collective safety, ABI conformance, and
fail-closed validation of the new Root-A path.  Runtime evidence is separate
and retained in this bundle; the static result does not substitute for it.

## Data-flow and collective review

- Mode 6 distributes only the additive pre-kernel k-to-R construction among
  CP2K's official k groups.  One service leader performs the nonlinear
  provider apply and pulls the R-to-k responses.  The code and diagnostics do
  not overclaim distribution of the root-serial kernel or reverse transform.
- Forward and reverse paths verify the official k-group range, group index,
  leader communicator size, and unique service-root relationship before
  entering the transaction.
- Every stage that may fail locally synchronizes failure first within each
  k group and then across the leader communicator before any later collective.
  This prevents a nonleader or leader error from stranding peers in a mismatched
  broadcast, reduction, or provider call.
- Global owner counts are reduced alongside the additive payload.  The
  provider accepts only exact-one ownership for each full-mesh point, including
  the zero contributions from leaders that do not own a given source.
- Only the service leader applies the coupled kernel and consumes the unique
  scalar/direct result.  Pulled response blocks are broadcast among leaders and
  folded by the unique irreducible source owner, so force/stress direct terms
  are not multiply counted.

## Provider ABI and transaction review

- CP2K requires partial ABI version 1, exactly 30 canonical descriptor words,
  and descriptor size 30 before using Mode 6.
- Each leader packs the full canonical descriptor, receives the service-root
  descriptor, unpacks it, and compares every field before payload reduction or
  import.  Forward and reverse payload sizes are cross-checked in int64 and
  bounded by the largest representable default integer before allocation.
- Each logical transaction receives a fresh nonzero two-word host identity,
  shared across leaders.  Provider-side validation rechecks transaction,
  phase/generation, geometry, model, kernel, mesh, compact phase, and onsite
  state.  A payload cannot silently be rebound to a different live descriptor.
- Integer arithmetic used by the new bridge checks regular-mesh products,
  batch/image ranges, and modular phase products before conversion to default
  integer indexing.

## Numerical and fail-closed review

- Local density, shell charge, packed payload, scalar/vector result, folded
  Fock response, overlap adjoint, direct gradient, and stress are all rejected
  if non-finite.
- Overlap and reverse-response blocks are checked for Hermiticity before any
  symmetric projection.  The folded Fock response is likewise checked before
  applying it to CP2K work matrices.
- Qualification mode retains hard dense oracles for energy, shell potential,
  folded Fock response, overlap adjoint, direct force, and direct stress.  Each
  residual has a fail-closed ceiling of `1e-10`.
- The one-point mesh preserves CP2K's established Gamma symmetry/weight path.
  Forward Mode 6 can be compared against explicit dense Gamma, while the
  derivative intentionally dispatches to the dense one-point fallback instead
  of claiming a partial-root reverse test.
- Qualification-only fault selection is strict.  Unknown and overlength
  selectors abort, as do injected nonleader failure, non-finite forward output,
  and anti-Hermitian reverse output.  Runtime fault-hook evidence covers all of
  these paths.

## Provider-side independent review

The retained provider adversarial review covers descriptor completeness and
corruption rejection, structural and state revalidation, mandatory transaction
identity, distributed ownership, finite payloads, int64 bounds, and dense
oracles.  It reports no unresolved defect at provider commit
`35e7942b60edd89bb407ab3da5768d3410af83f5`.

## Explicit limitation and unrelated follow-up

This Root-A bridge is a correctness-first first consumer: the nonlinear apply
and R-to-k reverse remain service-root serial.  Performance/scalability claims
must therefore be limited to the distributed additive k-to-R construction
until further work is implemented and measured.

One pre-existing, generic follow-up lies outside the bridge diff:
`src/cryssym.F:337` forms `nk(1)*nk(2)*nk(3)` in default integer arithmetic.
It should eventually receive the same checked-product treatment for malformed
extreme input.  It is not introduced by this change and does not alter the
review verdict for the bridge's own guarded mesh arithmetic.
