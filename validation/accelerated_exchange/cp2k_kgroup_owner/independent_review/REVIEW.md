# Independent review: CP2K `KGROUP_OWNER` precursor

Date: 2026-07-16

## Verdict

**PASS for the explicitly scoped correctness/ownership precursor.**  The
review found no source-level defect requiring a patch.  This verdict is not a
performance or distributed-provider-kernel qualification.

## Frozen-source identity and automated evidence

- The reviewed worktree files are byte-identical to the evidence snapshot:
  - `tblite_interface.F`: `be2dc79a785758fb712df108a75f6e750682dacfb7664919409e81b554ff58b5`
  - `tblite_types.F`: `e848fd796e69021a1ba21473d0fd7e2da3dbe485fe0e799e495ecd916889398e`
- `sha256sum -c SHA256SUMS` passes for the local archive and for
  `terok_linux`.
- `verify.py` passes all 11 frozen local cases.
- `verify_linux.py` passes frozen Linux P=1,2,4 cases.
- `git diff --check` is clean.
- A forced fresh CP2K precommit run checked both source files and returned
  `failed 0`.

Maxima recomputed from the frozen local matrix are
`dE=3.219647e-15 Ha`, `dVsh=1.387779e-17 Ha`, and
`dFfold=1.665335e-16 Ha`.  Final dense-reference differences are at most
`7.105427357601002e-15 Ha`, `7.16760127e-16 Ha/bohr`, and
`1.20000186143443e-6 bar`.

The frozen Terok P=1,2,4 matrix gives maximum in-process residuals
`dE=2.553513e-15 Ha`, `dVsh=0`, and `dFfold=6.941217e-18 Ha`; final
rank-count differences are at most `7.105427357601002e-15 Ha`,
`1.910821622e-16 Ha/bohr`, and `8.000006346e-7 bar`.

## Additional independent Linux cases

Two additional runs used the frozen Release binary whose SHA-256 is
`0d8c3257280cc27b09c08d5f2b3da7080ef68587996813a1871f728cffd27f90`.

1. CH4 Gamma, P=4, one k group of size four (`P > nred=1`): `PROGRAM
   ENDED`; maximum oracle residuals `dE=1.110223e-16 Ha`, `dVsh=0`,
   `dFfold=3.469447e-18 Ha`.  Versus the frozen dense Gamma reference:
   `dEtot=0`, `dFmax=4.52210711e-16 Ha/bohr`,
   `dStressMax=1.0000003385e-6 bar`.
2. H2 3x1x1 with time reversal, P=6, two k groups of size three
   (`P > nred=2`) and unequal star weights 2/3 and 1/3: `PROGRAM ENDED`;
   maximum oracle residuals `dE=5.551115e-17 Ha`, `dVsh=0`,
   `dFfold=1.387779e-17 Ha`.  Versus the accepted P=2 result:
   `dEtot=0`, `dFmax=0`, `dStressMax=5.0000380725e-7 bar`.

Independent-output hashes:

- Gamma P=4 output: `617d03e20cb66f76475e5b0ad530a1880107041d96cee01941d95dc088d8ffc2`
- H2 P=6 output: `e612353fd40f6de39907d3228a42615051e96edd96cdbc80ee2ea1c08dc7030d`

## Static MPI/weighting review

- The path obtains `kp_range`, `kp_dist`, `nkp_groups`, `para_env_kp`, and
  `para_env_inter_kp` only through `get_kpoint_info` and checks their mutual
  consistency (`tblite_interface.F:3701-3719`).  No direct MPI handle or
  newly constructed communicator appears in either reviewed file.
- CP2K itself constructs equal-sized, nonempty k groups and asserts
  `MOD(nkp,nkp_grp)==0` (`kpoint_methods.F:713-768`).  Consequently inactive
  or unequal-sized communicator groups are not legal states of this official
  metadata.  P>Nred is represented by multiple ranks inside each nonempty
  k group, as covered by the P=4 and P=6 Linux runs above.
- Every rank traverses every full-k block in identical order.  Exactly one
  owner-group leader sets `owner_count=1`; the leader inter-k communicator
  sums and checks that count before provider push
  (`tblite_interface.F:3771-3823`).  The unchanged overlap construction is
  then entered collectively by every global rank in the same order.
- One global source alone owns and applies the provider stream
  (`tblite_interface.F:3759-3768,3826-3830`).  Energy and shell potential have
  exactly one nonzero source before inter-group and intra-group propagation
  (`3831-3834`).  Each returned full-k Fock block has one source, is folded
  only by the irreducible owner, then distributed by the same two-level
  collective order (`3836-3866`).
- The provider receives full-mesh weights once.  The returned unweighted
  Fock response is folded with `w_full/w_irred` once
  (`tblite_interface.F:2519-2583`); no owner transfer applies another weight.
  The dense gradient path supplies full weights to the provider once and
  does not reweight the provider's total direct force/stress response
  (`5233-5239`).  The exchange force and stress are added once to the total
  response (`8021-8035`).

## Memory claim

The exact claim that can be accepted is limited to the CP2K owner-local
density allocation.  `density_local(nao,nao,nspin,nlocal)` is allocated with
`nlocal=nred/G` under CP2K's equal-group invariant
(`tblite_interface.F:3720-3733`; `kpoint_methods.F:735-768`).  Thus that one
allocation contains `nAO^2*nspin*nred/G` complex values per rank.

This is not a total-memory scaling result: the folded CP2K response remains
`nAO^2*nspin*nred` per rank (`tblite_interface.F:3748-3749`), the current
gradient path still materializes full-mesh arrays (`5121-5124`), and the
single provider stream retains full image-space state.  The frozen README
correctly scopes its formula to the CP2K density allocation.

## Nonblocking limitations that must remain explicit

1. The provider kernel is single-state on one global source.  There is no
   k-group speedup.  Frozen Linux wall time is 4.43--4.60 s for P=1--4 while
   aggregate child CPU time grows from 3.47 to 11.49 s.
2. The in-process oracle uses an independent dense provider evaluation and a
   separate full-array fold, but shares the CP2K star-block expansion,
   overlap construction, and orbit metadata with the candidate.  Frozen
   external dense runs independently cover final energy/force/stress; they do
   not constitute a separately stored full-Fock oracle.
3. "Unequal" coverage means unequal star sizes/weights, not unequal-sized or
   inactive communicator groups; the latter are excluded by CP2K's official
   k-group constructor.
4. The frozen verifier scripts reproduce their tables and enforce the
   in-process `1e-10` gate, but do not independently impose a numerical gate
   on final external observables or parse the launched MPI rank count.  The
   present review recomputed and inspected those values directly.
5. Genuine distributed provider work still requires the documented
   pre-kernel partial-accumulator merge.  Independently finalizing group
   streams and summing energies/Fock/forces/stress would be mathematically
   invalid.
