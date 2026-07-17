# g-xTB symmetry-star allocation/lifetime audit

Base: `g-xTB-pbc` at `68f677114f5829a32292171251150dd8e00ce458`.
Scope: `src/tblite_interface.F`.  Let `N=nao`, `S=nspin`, `K=nfull`,
`R=nred`, `L=nlocal`, and `B=image batch size`.  Counts below refer to dense
complex AO elements and exclude CP2K DBCSR/full-matrix storage that exists on
both sides of the selector.

## Existing exchange paths

| Path | Long-lived CP2K-side mesh storage during transaction | Full mesh? | Lifetime |
|---|---:|---|---|
| Dense forward oracle | `N^2 K (2S+1)` from density, overlap, full Fock | yes | one exchange build |
| Streamed forward (provider modes 1--4) | `N^2 S R` reduced density plus `N^2 (2S+2)` star work blocks | no, except explicit qualification oracle | one exchange build |
| K-group owner / partial-root (modes 5--6) | `N^2 S L` owner-local density plus star work blocks | no, except explicit qualification oracle | one exchange build |
| Dense reverse oracle | `N^2 K (S+2)` from density, overlap, full overlap adjoint | yes | one derivative build |
| Bounded streamed reverse | `N^2 R (S+1)` reduced density/adjoint plus star work blocks and provider `O(B N^2 S)` | no, except `QUALIFY` | one derivative build |

Thus forward exchange, Fock foldback, and reverse exchange already expand and
fold one complete-mesh star member at a time.  Those selectors are unchanged.

## Remaining unconditional symmetry full-mesh temporary

After native Fock mixing on a symmetry-reduced mesh,
`tb_prepare_gxtb_kpoint_mixer` unconditionally allocated
`mixed_fock_full(N,N,S,K)`, expanded every irreducible representative, checked
covariance, and immediately deallocated it.  The array is not consumed by the
provider or SCF; it is solely a post-mixer symmetry gate.  Its targeted peak is
exactly `N^2 S K` complex elements on every participating rank, once per native
mixer call.

## Isolated optimization

`CP2K_GXTB_SYMMETRY_STAR_CONTRACTION` is fail-closed and selects:

- `DENSE` (default, permanent oracle): unchanged full-mesh gate.
- `STREAMED`: for each irreducible source and each member of its star, expand
  one AO block, compare every operation-coset image, and immediately apply the
  weighted real adjoint via `tb_accumulate_gxtb_kpoint_response`.  This includes
  `fold_factor = w_full/w_irred`, calibrated `U_g`, and antiunitary time
  reversal.  The folded star is Hermitized and compared with the projected
  irreducible input at `1e-10`.
- `QUALIFY`: execute both gates in the same mixer call and print their maximum
  covariance residuals plus the streamed weighted round-trip residual.

The streamed gate has exactly three additional `N x N` complex work arrays,
independent of `K`, `R`, and `S`; it eliminates the targeted `N^2 S K` array.
The existing covariance tolerance remains
`max(1e-10, 100*eps_geo)`.  The weighted expand/fold round trip uses the same
`1e-10` gate as the existing symmetry/duality qualification paths.

## Execution-affinity qualification

Historical schema-1 matrix outputs used an outer shared `taskset` mask and
`mpiexec --bind-to none`. Their energy, force, stress, covariance, round-trip,
and allocation-counter comparisons remain numerical evidence after hash
verification, but their wall times are `legacy_timing_non_scaling` and are not
speedup evidence.

The clean 48-run rerun described in `EXACT_BINDING_RERUN.md` uses an explicit
ordered PE list with one logical CPU per MPI rank, Open MPI core binding and a
complete binding report, sticky lifetime sampling keyed by
`OMPI_COMM_WORLD_RANK`, and host-wide per-CPU reservation locks. The dense
implementation remains the numerical oracle; only fully revalidated schema-2
runs can be classified `production_scaling_eligible`.
