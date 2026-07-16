# Minimal save_tblite ABI for partial k-to-R accumulation

## Mathematical boundary

The Brillouin-zone transform is linear in each owned Bloch density/overlap
block, but the coupled exchange kernel is not additive in independently
completed k-point subsets.  Partial data must therefore be merged **before**
`cp2k_exchange_stream_batch_apply`.  Summing energy, shell potential, Fock,
force, or stress from independently finalized group streams is invalid because
it omits cross-group terms.

For one current image batch of capacity `B`, the additive forward state is:

- complex `A_R(nAO,nAO,B,nspin)`;
- complex `C_R(nAO,nAO,B,nspin)`;
- complex `V_R(nAO,nAO,B,nspin)`;
- real/complex onsite accumulators `gdiagP`, `gdiagSP`, and `gdiagSPS` with the
  provider's existing dimensions and normalization;
- integer `owner_count(Nk)`.

The reverse/derivative state needs the corresponding additive `B_R`,
`A_seed,R`, and `V_seed,R` buffers before the adjoint R-to-k transform.  These
names describe the mathematical roles; the public object may keep them opaque.

## Required public operations

A minimal Fortran-facing contract is:

```fortran
call cp2k_exchange_partial_begin(partial, calc, mol, ecache, wfn, &
   nmesh, xkp, weight, batch_first, batch_last, error)

call cp2k_exchange_partial_push(partial, ik, density, overlap, error)

call cp2k_exchange_partial_export_forward(partial, descriptor, &
   owner_count, forward_buffer, error)

! CP2K performs SUM over its inter-k communicator here.

call cp2k_exchange_partial_import_forward(stream, descriptor, &
   owner_count, forward_buffer, error)

call cp2k_exchange_stream_batch_apply(stream, calc, mol, ecache, error)

call cp2k_exchange_partial_export_reverse(stream, descriptor, &
   reverse_buffer, error)

! CP2K performs the required SUM/reduce-scatter here.

call cp2k_exchange_partial_import_reverse(stream, descriptor, &
   reverse_buffer, error)
```

`forward_buffer` and `reverse_buffer` should be opaque, contiguous pack objects
with public size/pack/unpack queries; CP2K must never access private stream
components.  save_tblite must not call MPI.  CP2K owns `MPI_Allreduce(SUM)`,
`MPI_Reduce_scatter`, or point-to-point scheduling as appropriate.

## Validation and state machine

The descriptor must be immutable and carry (or cryptographically identify):

- mesh dimensions, canonical k ordering, coordinates, weights, common twist,
  and input-to-grid permutation;
- `Nk`, `nAO`, `nspin`, shell dimensions, transform backend and normalization;
- current image-batch bounds and capacity;
- model/calculator identity, geometry/cell generation, charge-state generation,
  and kernel/ACP signature;
- buffer-layout ABI version and scalar precision.

`import_forward` must reject a descriptor mismatch, non-finite data, a stale
generation, or any globally summed `owner_count(ik) /= 1`.  It is legal to
export after an arbitrary nonempty set of uniquely owned k blocks; it is not
legal to apply until a complete validated import has occurred.  A second
import, push-after-import, duplicate owner, mismatched batch, or apply-before-
import must be a hard error.  Reverse export is legal only after a successful
apply of that same descriptor; reverse import/finalization must obey the same
single-use generation rules.

The forward collective volume per current batch is
`O(B*nAO^2*nspin)` complex values plus `O(Nk)` owner integers.  That is the
minimal public boundary needed to distribute the k-to-R construction exactly.
Scaling the coupled provider kernel itself additionally requires disjoint
image/AO-pair work ownership after this merge (or a reduce-scatter schedule),
with additive outputs reduced only after the kernel stage.

## Mandatory provider tests

1. Two partial streams whose disjoint k sets form the full mesh must reproduce
   the existing dense stream for forward buffers, energy, shell potential, and
   every returned Fock block.
2. Duplicate, missing, or out-of-range k ownership must fail.
3. Permuted input ordering must pass only with a matching descriptor and must
   reproduce the canonical result.
4. Mismatched mesh shift, weights, batch bounds, geometry generation, model,
   spin, precision, or transform normalization must fail before apply.
5. Reverse buffers, analytic forces, and stress must match the dense stream and
   central finite differences for 1D, 2D, and 3D periodicity, time reversal,
   UKS, SPGLIB/K290, and shifted meshes.
6. Results must be invariant to one, two, and four CP2K k groups and to a rank
   count greater than `Nk`.
