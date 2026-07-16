# Periodic g-xTB symmetry/k-point V2 audit

Audit baseline:

- CP2K: `28df9380abb327d56bbf216d2469a1fd8c953fc0`
- save_tblite: `257ba442684c39454175e5192c8a2342b4c6380f`
- V2 worktree: `/private/tmp/cp2k_gxtb_symmetry_v2`
- V2 branch: `gxtb-symmetry-v2-audit`
- The qualified V1 source, binary, and benchmark data are out of scope and were not modified.

## 1. Main conclusion

The complete Brillouin-zone sum is physically coupled for g-xTB exchange. It is not merely a
CP2K storage convention that can be removed by treating irreducible k-points as independent
GFN1/GFN2 contributions.

CP2K already has the necessary orbit map and calibrated real-linear AO transforms:

- `tb_build_gxtb_kmesh_map` builds representative, star-member, coset-operation, and fold-weight
  maps.
- `tb_transform_gxtb_kpoint_matrix` implements the forward `T_g` and the real adjoint `T_g*`,
  including atom permutations, AO rotations, Bloch phases, and antiunitary time reversal.
- `tb_expand_gxtb_kpoint_density` and `tb_fold_gxtb_kpoint_matrix` are the V1 full-mesh reference
  projector and adjoint.

However, the provider API currently requires complete arrays. In save_tblite,
`cp2k_exchange_kmesh` rejects anything except `PRODUCT(nmesh) == nk`, uniform weights, and a
complete Fourier mesh. `exchange_fock%get_KFock_kmesh` then couples all k-points through the dense
`k -> R -> k` transforms in `apply_bvk_Kkernel_inplace`. The reverse path in
`get_KGrad_kmesh` has the same global coupling.

Therefore, the correct V2 target is a streamed orbit-aware provider interface: CP2K generates one
star member from its irreducible representative with `T_g`, save_tblite accumulates the coupled
real-space intermediates, and CP2K folds each returned response immediately with `T_g*`. The
physical sum over every star member remains; the persistent/full dense k-mesh arrays do not.

## 2. Current materialization and implementation points

| Role | Current implementation | V2 insertion point |
|---|---|---|
| Orbit topology and weights | `src/tblite_interface.F::tb_build_gxtb_kmesh_map` | Keep as canonical CP2K orbit descriptor; extend only with immutable transform-cache metadata. |
| Phase calibration | `src/kpoint_methods.F::calibrate_symmetry_phases` | Keep as the only gauge authority; cache derived descriptors only after this routine succeeds. |
| AO star transform / adjoint | `src/tblite_interface.F::tb_transform_gxtb_kpoint_matrix` | Factor blockwise `T_g`/`T_g*` application for streamed provider calls. |
| Full density expansion | `src/tblite_interface.F::tb_expand_gxtb_kpoint_density` | Retain as qualification oracle; production V2 calls a one-star-member generator. |
| Full response folding | `src/tblite_interface.F::tb_fold_gxtb_kpoint_matrix` | Retain as qualification oracle; add `tb_accumulate_gxtb_kpoint_response` for immediate adjoint accumulation. |
| SCF exchange | `src/tblite_interface.F::tb_build_gxtb_kpoint_exchange` | Replace simultaneous `density`, `overlap`, and `exchange_fock_full` arrays with begin/push/finalize/pull provider phases. |
| Exchange gradient | `src/tblite_interface.F::tb_gxtb_kpoint_exchange_gradient` | Stream full-star density/overlap into the reverse plan and fold each overlap response immediately. |
| ACP Fock | `src/tblite_interface.F::tb_build_gxtb_kpoint_acp` | ACP is separable by k; expose provider block evaluation and fold each block immediately. |
| ACP gradient | `src/tblite_interface.F::tb_gxtb_kpoint_gradient` | Replace `acp_density_k(:,:,:,nfull)` with a provider accumulator accepting one weighted k block at a time. |
| Provider ABI | `save_tblite/src/tblite/cp2k_compat.f90` | Add stable plan/stream routines while retaining full-mesh calls as reference. |
| Coupled exchange | `save_tblite/src/tblite/exchange/fock.f90` | Split `build_KFock_kmesh` and `get_KGrad_kmesh` into accumulation, real-space kernel, and output-block stages. |
| Provider invariant cache | `save_tblite/src/tblite/exchange/cache.f90` | Cache BvK kernel, mesh permutation/twist, phase/FFT plan, and validation fingerprint. |
| ACP implementation | `save_tblite/src/tblite/acp/type.f90` | Split `get_acp_kmesh[_gradient]` into reusable plan plus per-k accumulation. |
| MPI k groups | `src/kpoint_methods.F::kpoint_env_initialize`, `kpoint_type%para_env_kp`, `%para_env_inter_kp` | Reuse existing communicators after the serial streamed reference is correct. |
| Mesh-changing restart | `src/kpoint_io.F`, call in `src/qs_initial_guess.F` | Introduce a versioned mesh-aware restart and a validated source-R/target-k interpolation path. |

The largest V1 CP2K allocations occur in `tb_build_gxtb_kpoint_exchange`: `density` is
`nao^2*nspin*nfull`, `overlap` is `nao^2*nfull`, and `exchange_fock_full` is
`nao^2*nspin*nfull`, all complex. The gradient path materializes analogous input and adjoint
arrays. save_tblite additionally allocates several `nao^2*nk` forward/reverse intermediates.

## 3. Prioritized implementation plan

### P0: Freeze reference oracles and add timers

1. Keep all V1 full expansion/folding routines callable as a qualification backend.
2. Add separate timers and peak-workspace counters for CP2K star expansion, dense overlap build,
   provider exchange, response folding, ACP, and provider reverse mode.
3. Freeze full-grid versus reduced-grid numerical artifacts before enabling any optimized path.

Go/no-go: no optimized default until all energy, duality, force, and stress oracles below are
green on the same binary.

### P1: Cache immutable metadata only

1. CP2K: cache the resolved AO-rotation slot for each `(full star member, coset operation)`. The
   current V2 micro-patch does exactly this and changes no floating-point arithmetic.
2. CP2K: next cache compact transform descriptors after phase calibration: atom permutation,
   signed operation, selected gauge/mode, and atom phase arguments. Do not cache dense `U_g`.
3. save_tblite: cache `exchange_bvk_kernel`, Fourier mesh permutation, common twist, and dense
   phase tables (or FFT plan) in `exchange_cache`, keyed by geometry generation, `nmesh`, k-list,
   and weights. The current code rebuilds the BvK kernel and phases on every SCC evaluation and
   again for the gradient.
4. ACP cache keys must include the basis/charge generation as well as geometry. q-vSZP makes a
   geometry-only ACP cache unsafe unless basis invariance is demonstrated.

Go/no-go: cache invalidation tests must cover coordinate move, homogeneous cell strain, mesh or
shift change, symmetry backend change, and every charge-dependent basis update.

### P2: Stream ACP and simple star folding

1. Add a save_tblite ACP plan that computes projector cells/phases once and returns one `H_ACP(k)`
   block at a time.
2. Add `tb_accumulate_gxtb_kpoint_response` in CP2K. It applies every real adjoint in the operation
   coset and adds `fold_factor/op_count` directly to the irreducible destination.
3. Keep `tb_fold_gxtb_kpoint_matrix` as an independent full-array oracle.
4. Split ACP gradient into begin/accumulate/finalize so CP2K never allocates
   `acp_density_k(:,:,:,nfull)`.

ACP is the safe first streaming target because its k blocks are separable before the final weighted
gradient contraction.

### P3: Stream the coupled SCF exchange

Add a provider plan with a stable reference implementation along these lines:

1. `exchange_plan_init`: validate the complete regular mesh once; build/cache BvK kernel and
   Fourier descriptor.
2. `exchange_plan_push_k`: CP2K generates a star block `(P_k,S_k)` from one irreducible
   representative. save_tblite forms the local nonlinear intermediates (`S P`, `S P S`, and `P`)
   and accumulates their inverse Fourier transforms and Brillouin-zone diagonal terms.
3. `exchange_plan_apply_r`: apply the image-resolved exchange kernels after all physical star
   members have been accumulated.
4. `exchange_plan_get_k`: reconstruct one unweighted `F_k`; CP2K immediately applies `T_g*` and
   accumulates it into the irreducible Fock matrix. Energy and shell potential are accumulated
   exactly once with the complete-mesh weights.
5. Retain the current `cp2k_exchange_kmesh` call as the independent reference backend until the
   streamed path passes permutation, duality, and finite-difference gates.

This removes simultaneous CP2K full-mesh inputs/outputs. A first implementation may still keep a
small number of provider `nao^2*nR` work arrays; AO tiling can reduce that peak later.

### P4: Stream reverse mode for forces and stress

1. Split `get_KGrad_kmesh` at the same mathematical boundaries as P3.
2. Regenerate star blocks from irreducible inputs if retaining all forward intermediates would cost
   more memory than a second pass.
3. Emit each unweighted overlap-gradient block and fold it immediately with the real adjoint.
4. All atomic and strain derivatives returned by the provider remain complete primitive-cell BZ
   contractions and must not receive a second k weight.
5. Keep `tb_qualification_gxtb_duality`, full-response covariance, folded-response covariance, and
   the permuted-full-mesh oracle active in qualification builds.

### P5: Regular-mesh FFT backend

The dense transforms are in save_tblite, not CP2K:

- `build_KFock_kmesh` constructs `phase_forward`/`phase_inverse`.
- `apply_bvk_Kkernel_inplace` executes both directions as dense complex matrix products.
- `reverse_bvk_Kkernel` and `get_KGrad_kmesh` require the exact adjoint normalization.

Implementation:

1. Map arbitrary input ordering to integer `(ix,iy,iz)` mesh indices and map BvK representatives
   to the same canonical order.
2. Factor a common shifted-mesh twist into diagonal pre/post phases.
3. Add a batched 1D/2D/3D complex FFT backend over the mesh dimension for AO-pair tiles.
4. Keep dense GEMM as the reference and as the small-mesh/irregular-GENERAL fallback. Select the
   crossover by measured timings, not mesh size intuition.
5. Verify the reverse transform as the real adjoint, including normalization and antiunitary star
   folding.

This changes the transform cost per AO pair from `O(Nk^2)` to `O(Nk log Nk)` for regular meshes;
it does not remove the physical star sum.

### P6: MPI k groups

The current `tb_collect_kpoint_density` gathers every irreducible density through the global BLACS
environment, after which every rank executes the replicated dense provider calculation. Existing
CP2K k groups therefore help diagonalization but do not yet parallelize the coupled g-xTB provider
work.

The selectable `KGROUP_OWNER` precursor now validates the CP2K side of that decomposition without
claiming provider-kernel speedup:

1. It obtains `kp_range`, `kp_dist`, `nkp_groups`, `para_env_kp`, and
   `para_env_inter_kp` exclusively through `get_kpoint_info` and checks their consistency.
2. A rank stores only its group's `nred/nkp_groups` irreducible density blocks. The owner-group
   leader expands each full-mesh star member exactly once; the inter-group leader communicator
   transfers one Bloch block at a time and verifies an exact owner count of one.
3. Because the public provider stream cannot merge partial k-to-R accumulators, exactly one global
   source owns the save_tblite stream and coupled kernel. It receives all blocks sequentially, so
   neither a complete CP2K density mesh nor replicated provider states exist. This is an exact
   ownership and bounded-communication baseline, not scalable k-group execution.
4. The source returns one full-mesh Fock block at a time. Only the corresponding irreducible owner
   folds it; a leader-inter-group reduction followed by an intra-group reduction distributes the
   complete irreducible response. Energy and shell potential follow the same communicator order.
5. A qualification-only switch evaluates the dense full-mesh provider in the same SCF iteration
   and compares energy, shell potential, and the folded Fock response at `1e-10`.

The remaining public save_tblite boundary is deliberately small. For each current image batch the
provider must expose an opaque, additive partial-accumulator object containing the forward k-to-R
buffers `A_R`, `C_R`, and `V_R`, plus `gdiagP`, `gdiagSP`, and `gdiagSPS`. It must export and import
those buffers before `cp2k_exchange_stream_batch_apply`; CP2K then performs `MPI_Allreduce(SUM)`
over `O(B*nAO^2*nspin)` complex values and `O(Nk)` integer owner counts. The provider must reject a
merge unless every k point has exactly one owner and the mesh/order/weights/twist, model and
geometry generation, kernel signature, current image-batch bounds, AO/spin dimensions, and
transform convention all match. The reverse stage analogously exports/imports additive `B_R`,
`A_seed,R`, and `V_seed,R` buffers before their adjoint transform.

Only after this pre-kernel merge may groups execute disjoint image or AO-pair tiles and sum their
additive scalar, shell, Fock, force, and stress contributions. Independently finalized group
results are mathematically invalid because the coupled nonlocal kernel contains cross-group terms.
save_tblite should remain MPI-free: CP2K owns all collectives, while the provider owns the state
machine, validation, packing, and pure-compute stages.

### P7: Validated `n^3 -> (n+1)^3` density restart

The current `.kp` format (`src/kpoint_io.F`, version 1) stores real-space density images and cell
labels, but no source mesh, shift, irreducible/full-k map, or alias convention. The reader copies
only cells that exist in the target neighbor-list map. That is not a sufficient contract for a
validated 5^3-to-6^3 transfer.

1. Add a version-2 restart header containing source `nmesh`, shift, gamma-centering, full/reduced
   k coordinates and weights, spin, and canonical BvK cell representatives.
2. Expand the source irreducible density with the source symmetry map, inverse-transform to a
   canonical source-R representation, embed common centered R coefficients into the target BvK
   lattice, zero-pad unresolved coefficients, and forward-transform to the target mesh.
3. Apply target little-group projection and Hermitization.
4. Validate electron count `sum_k w_k Tr[P(k)S(k)]`, finiteness, Hermiticity, star covariance, and
   absence of BvK alias collisions. Fall back to the normal atomic/RS guess if any gate fails.
5. Do not require identical final SCF iteration histories. Require identical converged observables
   and demonstrate a statistically useful iteration/time reduction on the benchmark set.

## 4. Generic CP2K boundary and HF/post-HF reuse

The V2 work should not become a second, g-xTB-private k-point subsystem. CP2K already contains
the beginnings of a reusable implementation in `src/kpoint_methods.F`:

- `kpoint_env_initialize` owns `kp_range`, `para_env_kp`, and `para_env_inter_kp`;
- `kpoint_density_transform` expands irreducible AO densities and immediately accumulates the
  weighted real-space density;
- `symtrans_phase` applies the atom permutation, AO rotation, Bloch phase, and antiunitary time
  reversal;
- `rskp_transform` transforms real-space DBCSR blocks to one k point.

The current g-xTB code needed a dense-array form of the same real-linear orbit operation because
the save_tblite provider accepts dense full-mesh arrays. The long-term refactor should move the
immutable topology and transform descriptors to a method-independent layer, while leaving storage
and physical kernels behind explicit adapters.

### 4.1 Proposed layers

| Layer | Reusable contract | Must not leak into the contract |
|---|---|---|
| `kmesh_orbit_plan` | Periodic dimensions, integer regular-grid coordinates, shift/twist, canonical k and BvK-R order, full-to-irrep stars, operation cosets, time reversal, weights, and a strict invalidation fingerprint. | g-xTB exchange parameters, RI metrics, Coulomb-head models, orbital occupations. |
| `krep_transform_plan` | Cached representation-specific permutation/rotation/phase descriptors and exact forward/real-adjoint calls. Start with AO one-body matrices and dense/DBCSR backends. | An assumption that AO, RI auxiliary, and MO-band representations share the same `U_g`. They do not. |
| `kfourier_plan` | Dense reference and batched FFT implementations of regular-mesh `k <-> R`, exact adjoint normalization, arbitrary input ordering, compatible shifted meshes, and 1D/2D/3D partial PBC. | Coulomb singularity treatment, minimum-image tie rules, screening, or nonlinear model kernels. |
| `kspace_work_plan` | Ownership of k blocks, R blocks, and AO/RI-pair tiles using the existing k-group communicators, with deterministic reductions and a serial reference schedule. | Method-specific tensor contractions or a direct MPI dependency in save_tblite. |
| Method adapter | Streamed begin/push/apply/get interface and method-specific force/stress reverse path. | Reimplementation of mesh canonicalization, phase gauges, or communicator construction. |

`U_g` is reusable only at a fixed representation level. AO density and Fock matrices can use the
current atom/AO transform. RI-HFX needs a separate auxiliary-basis representation for two- and
three-center tensors. MO-based MP2/RPA/GW additionally needs band sewing matrices and a policy for
degenerate subspaces; a raw AO rotation is not a valid replacement. Momentum-transfer and
two-particle objects live on pair or quadruple k orbits, not on the single-k stars used by a
one-body density.

### 4.2 Concrete reuse points found in CP2K

| Method/path | Existing code | Reusable part | Method-specific remainder |
|---|---|---|---|
| Ordinary SCF/DFT | `src/kpoint_methods.F::kpoint_density_transform`, `symtrans_phase`, `rskp_transform` | Orbit descriptors, AO `U_g/U_g*`, phase cache, Fourier plan, existing k-group communicators. | Density construction, smearing/occupations, local XC/grid integration. Local/semi-local DFT may sum irreducible densities directly and does not require the coupled exchange protocol. |
| RI-HFX with k points | `src/hfx_ri_kp.F::get_pmat_images`, `hfx_ri_update_ks_kp`, `hfx_ri_update_forces_kp` | Source AO-density expansion, k/R transform, R/image ownership, deterministic k/R/pair scheduling. `KP_NGROUPS` is an existing method-local prototype for subgroup work distribution. | Local RI domains, 2c/3c integrals, truncated/short-range potential, bump functions, sparsity screening, delta-P update, and all integral/bump derivatives. |
| RI metric and Coulomb matrix | `src/mp2_ri_2c.F::RI_2c_integral_mat`, `compute_V_by_lattice_sum`; `src/kpoint_coulomb_2c.F::lattice_sum` | Batched regular R-to-k transform and R/pair distribution. The present code explicitly evaluates cos/sin for each residual cell and k point. | Coulomb/truncated-Coulomb operator, lattice-sum convergence, even-mesh restriction, extrapolation, and RI factorization. |
| Cubic RPA/GW | `src/rpa_gw_kpoints_util.F::real_space_to_kpoint_transform_rpa`, `invert_eps_compute_W_and_Erpa_kp`, `compute_Wc_real_space_tau_GW` | Fourier plan, canonical cell ordering, k/R ownership, and subgroup redistribution. `real_space_to_kpoint_transform_rpa` is a second explicit dense phase loop that should eventually use the same backend. | Polarizability, dielectric inversion, frequency/time transforms, screened interaction, and atom-pair minimum-image tie averaging. |
| GW/RPA long-wavelength terms | `src/rpa_gw_kpoints_util.F::compute_wkp_W`; `src/rpa_gw.F::compute_eps_head_Berry`, `kpoint_sum_for_eps_inv_head_Berry` | Mesh coordinates and orbit-invariant bookkeeping only. | Tailored `1/k` or `1/k^2` integration weights, Berry head/wing limits, dimensional corrections, and extrapolation. These must remain explicit kernel-head policies. |
| Periodic MP2 or other two-particle methods | Current k-point infrastructure plus future method code; the present `mp2_gpw.F` k-point branch primarily serves cubic RPA/GW rather than a general periodic canonical-MP2 implementation. | Mesh canonicalization, representation transforms, FFT, and MPI pair-group scheduling. | Momentum-conserving pair/quadruple enumeration, occupied/virtual band gauges, denominators, exchange antisymmetry, q=0 corrections, amplitudes, and their derivatives. |

For RI-HFX, the useful decomposition is therefore not "make HFX call save_tblite". Both methods
can consume a common orbit/Fourier/work plan, but their kernel stages remain independent. The
existing RI-HFX implementation already works largely in real-space image tensors and distributes
`(iatom,jatom,R)` work; it may benefit more from a shared transform and communicator layer than
from materializing full k stars.

### 4.3 Kernel heads, staggered meshes, and correctness boundary

An FFT or symmetry orbit engine is algebraic. It must never silently decide how a singular or
nearly singular physical kernel is integrated. In particular:

- HF/RI-HFX chooses Coulomb versus truncated/short-range operators and lattice cutoffs;
- RPA/GW uses tailored k weights, dielectric head/wing limits, Berry matrix elements, and
  dimensional extrapolations;
- a staggered occupied/virtual or k/q construction needs two mesh descriptors and an explicit
  proof that differences and momentum-conserving sums close on the target grid;
- shifted meshes require a common twist for FFT use; incompatible or irregular pairs fall back to
  the dense reference;
- force and stress paths need derivatives of the method kernel, basis, metric, and any head
  correction. The adjoint Fourier/orbit transform supplies only the algebraic pullback.

Accordingly, the reusable API should carry opaque caller metadata or callbacks for kernel-head
terms; it should not contain a universal `q=0` value. A method adapter must declare whether its
weights are physical BZ weights, tailored quadrature weights, or an unweighted response before a
transform is allowed.

### 4.4 Safe implementation sequence for the shared layer

1. Extract a read-only orbit descriptor from the already validated `kpoint_type/kpoint_sym_type`
   data. Make `kpoint_density_transform` and the g-xTB qualification path consume it without
   changing arithmetic.
2. Extract a phase iterator/dense-reference Fourier plan and route `rskp_transform` plus one
   RPA transform through it. Preserve the existing routines as bitwise/reference backends.
3. Add AO dense and DBCSR representation adapters; add RI auxiliary and MO sewing adapters only
   with independent covariance tests.
4. Add FFT behind a per-call capability check. Then connect RI Coulomb/RPA transforms one at a
   time against their existing explicit cos/sin loops.
5. Unify scheduler metadata with `para_env_kp/para_env_inter_kp`; retain RI-HFX `KP_NGROUPS` and
   current RPA subgroup routes as independent MPI oracles until their ownership/reduction tests
   pass.
6. Only after single-k one-body paths are qualified, add pair/transfer-momentum orbit enumeration
   for HF/post-HF. Do not infer two-particle star weights from one-body weights.

This is an engineering reuse opportunity, not by itself a scientific novelty claim. Symmetry
stars, FFTs, and MPI decomposition are standard components. A paper must tie any novelty claim to
the model-specific coupled formulation and its validated periodic energies/forces/stress (or to a
demonstrably new algorithm), not to the mere combination of those standard techniques.

### 4.5 Additional HF/post-HF qualification gates

- `kpoint_density_transform`: old versus plan-backed full/reduced AO density for 1D/2D/3D,
  time reversal, shifted meshes, UKS, and rank-count changes;
- RI-HFX: energy and Fock with delta-P on/off, full versus reduced mesh, `KP_NGROUPS=1/2/4`, then
  analytic force/stress versus finite differences;
- RI Coulomb: explicit cos/sin versus FFT for Coulomb and truncated operators, even/odd supported
  meshes, partial PBC, and lattice-sum-size convergence;
- RPA/GW: correlation energy and selected `W(k,iw)` blocks with uniform and tailored weights,
  head correction on/off, full versus reduced symmetry, and dense versus FFT;
- future two-particle orbit code: brute-force versus orbit-reduced momentum quadruplets on 2x2x1
  and 2x2x2 meshes, invariance to MO phase and rotations within degenerate bands, and staggered
  mesh closure/fallback tests.

## 5. Numerical oracle matrix

Every optimized result is compared to a full-grid, dense-transform, serial V1/reference path from
the same source geometry and model parameters.

| Periodicity / case | Mesh variants | Symmetry variants | Required oracles |
|---|---|---|---|
| 0D H2O in a box | implicit Gamma; explicit 1x1x1 | full grid; Gamma symmetry | Energy, charges, SCF convergence, analytic/FD force and stress-path consistency. |
| 1D H2 chain, `PERIODIC X` | 3x1x1 and 4x1x1; Gamma-centred and one compatible shift | full; inversion/time reversal; atomic symmetry where supported | Primitive versus 3x1x1 Gamma supercell energy/forces/axial stress; full/reduced Fock duality; force FD. |
| 2D BN or molecular sheet, `PERIODIC XY` | 2x2x1 and 3x3x1; shifted 2x2x1 | full; SPGLIB; inversion/time reversal | Primitive versus 2x2x1 Gamma supercell energy/forces/in-plane stress; full/reduced equality; strain FD for xx, yy, xy. |
| 3D Si | 2x2x2, 3x3x3; shifted 2x2x2 | full; K290; SPGLIB; inversion/time reversal | Energy, SCF, density/Fock covariance and duality, forces, full stress, dense-versus-FFT. |
| 3D low-symmetry molecular crystal | 2x2x2 and 3x3x3 | full; SPGLIB; permuted explicit GENERAL full mesh | Ordering invariance, energy, forces, stress, ACP image coverage. |
| UKS O2 periodic cell | 3x1x1 | full versus time reversal | Alpha/beta populations, energy, spin-resolved duality, Fock Hermiticity, forces/stress. |
| Restart | 2->3 small tests; production 5->6 and optionally 6->7 | full and symmetry-reduced; Gamma-centred and shifted | Seed charge/Hermiticity/covariance; cold-versus-restart converged energy, forces, stress; iterations and wall time. |
| MPI | representative 1D/2D/3D cases | 1, 2, and 4 k groups; multiple ranks per group | Same converged observables and duality as serial; no rank-count-dependent cache state. |

Minimum numerical gates:

- full versus reduced total energy: `1e-10 Ha` for fixed arithmetic; `1e-9 Ha` for FFT/MPI paths;
- full/reduced variational duality and covariance residuals: `<= 1e-10`;
- dense versus FFT transform relative residual: `<= 1e-11` before SCF;
- analytic force versus central FD: `<= 1e-6 Ha/bohr`;
- analytic strain derivative versus central FD: `<= 1e-6 Ha` in the tested component;
- primitive/supercell energy per primitive cell: `<= 1e-8 Ha`, with correspondingly mapped forces
  and stress;
- restart and cold-run converged energy: `<= 1e-9 Ha`; force/stress differences within their FD
  tolerances;
- every run must terminate with the same electron/spin populations and a genuine SCF convergence
  criterion, not merely a similar energy.

Negative tests are required for irregular/incomplete GENERAL meshes, replication in a nonperiodic
direction, inconsistent mesh weights, stale geometry/basis caches, and a shifted mesh that cannot
be mapped to a regular Fourier grid.

## 6. Current isolated V2 micro-patch

The worktree currently adds `tb_gxtb_kmesh_map_type%rotation_slot`. The slot is resolved once when
the immutable star map is built and is passed to every forward/adjoint AO transform. The old code
linearly searched `kpoints%ibrot` on every transform call. No phase, matrix multiplication, weight,
or accumulation order was changed.

Current local checks:

- `./make_pretty.sh src/tblite_interface.F src/tblite_types.F`: pass, 2/2 files;
- `python3 tools/precommit/precommit.py src/tblite_interface.F src/tblite_types.F`: pass;
- `git diff --check`: pass;
- isolated GNU Fortran compilation against the pinned Terok V1 modules and save_tblite provider:
  pass for the exact formatted `tblite_types.F` and `tblite_interface.F` sources. Because
  `tb_gxtb_kmesh_map_type` is public, the isolated module set correctly rebuilt its two direct
  consumers, `qs_environment_types.F` and `cp_dbcsr_output.F`, without modifying the V1 build.
  The isolated artifacts are under
  `/home/kuehne88/work/codex-gxtb-symmetry-v2-syntax-20260714`.

It must remain non-default/unmerged until a separately linked V2 binary passes the full/reduced
energy, duality, force, stress, time-reversal, and shifted-SPGLIB runtime tests.
