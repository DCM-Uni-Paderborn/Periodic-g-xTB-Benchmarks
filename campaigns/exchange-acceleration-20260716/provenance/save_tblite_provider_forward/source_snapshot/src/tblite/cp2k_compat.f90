! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later

!> Compatibility helpers for CP2K's low-level tblite integration.
module tblite_cp2k_compat
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use multicharge, only : get_charges
   use tblite_acp_cache, only : acp_cache
   use tblite_acp_type, only : acp_image_type, get_acp, get_acp_images, &
      & get_acp_kmesh, get_acp_gradient, get_acp_image_gradient, &
      & get_acp_kmesh_gradient
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_cache, only : basis_cache
   use tblite_basis_type, only : basis_type, get_cutoff, get_pair_cutoff
   use tblite_blas, only : gemv
   use tblite_container_cache, only : container_cache
   use tblite_cutoff, only : get_lattice_points
   use tblite_exchange_cache, only : exchange_cache
   use tblite_exchange_fock, only : exchange_fock
   use tblite_integral_multipole, only : multipole_cgto, multipole_grad_cgto
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_scf_potential, only : potential_type, add_pot_to_h1
   use tblite_scf_info, only : scf_info, atom_resolved, shell_resolved
   use tblite_scf_mixer_cache, only : mixer_cache_container
   use tblite_scf_mixer_diis, only : diis_input, diis_mixer, new_diis
   use tblite_scf_mixer_input, only : mixer_mode
   use tblite_scf_mixer_simple, only : new_simple, simple_input, simple_mixer
   use tblite_wavefunction_spin, only : magnet_to_updown, updown_to_magnet
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_anisotropy, get_anisotropy_gradient, &
      & get_hamiltonian, get_hamiltonian_gradient, hamiltonian_image_type
   implicit none
   private

   public :: cp2k_get_mixer_dimension, cp2k_get_scc_mixer_dimension
   public :: cp2k_get_acp_cutoff
   public :: cp2k_potential_mixer_type, new_cp2k_potential_mixer
   public :: cp2k_set_scf_options
   public :: cp2k_update_basis, cp2k_get_cgto
   public :: cp2k_multipole_cgto, cp2k_multipole_grad_cgto
   public :: cp2k_build_gamma_h0, cp2k_build_image_h0, cp2k_build_image_acp, &
      & cp2k_build_acp_projectors
   public :: cp2k_prepare_exchange, cp2k_exchange_kpoint, cp2k_exchange_kmesh, &
      & cp2k_exchange_kmesh_gradient, cp2k_exchange_kpoint_gradient, &
      & cp2k_gamma_gradient
   public :: cp2k_exchange_stream_type, cp2k_exchange_stream_begin, &
      & cp2k_exchange_stream_push, cp2k_exchange_stream_apply, &
      & cp2k_exchange_stream_get, cp2k_exchange_stream_reverse_apply, &
      & cp2k_exchange_stream_reverse_get, cp2k_exchange_stream_end, &
      & cp2k_exchange_stream_has_full_mesh_storage
   public :: cp2k_image_h0_gradient
   public :: cp2k_acp_image_gradient, cp2k_acp_kmesh, &
      & cp2k_acp_kmesh_gradient, cp2k_apply_basis_gradient
   public :: acp_image_type, hamiltonian_image_type

   integer, parameter, public :: cp2k_exchange_stream_reduced = 1
   integer, parameter, public :: cp2k_exchange_stream_oracle = 2

   !> CP2K-facing wrapper for g-xTB's native potential mixer.
   !>
   !> CP2K owns the packing of all real Fock degrees of freedom into one flat
   !> vector.  This wrapper deliberately owns only the mixer state and mirrors
   !> g-xTB's iterator schedule: simple mixing with damping 0.2 for the first
   !> three Fock builds, DIIS precollection on the third build, and DIIS from
   !> the fourth build onward.
   type, public :: cp2k_potential_mixer_type
      private
      type(simple_mixer) :: simple
      type(diis_mixer) :: diis
      type(mixer_cache_container) :: simple_cache
      type(mixer_cache_container) :: diis_cache
      integer :: ndim = 0
      integer :: iteration = 0
      real(wp) :: residual = 0.0_wp
      logical :: initialized = .false.
   contains
      procedure, public :: reset => cp2k_potential_mixer_reset
      procedure, public :: advance => cp2k_potential_mixer_advance
      procedure, public :: get_error => cp2k_potential_mixer_get_error
   end type cp2k_potential_mixer_type

   !> Stateful collector for Brillouin-zone-coupled exchange.
   !>
   !> The default reduced mode retains only three weighted k-to-R
   !> intermediates per spin and reconstructs requested Fock blocks on demand.
   !> The separately selectable oracle mode retains all Bloch blocks and
   !> dispatches the unchanged complete-mesh evaluator, including reverse mode.
   type, public :: cp2k_exchange_stream_type
      private
      logical :: active = .false.
      logical :: applied = .false.
      logical :: reverse_applied = .false.
      integer :: nk = 0
      integer :: nao = 0
      integer :: nspin = 0
      integer :: nsh = 0
      integer :: mode = 0
      integer :: nmesh(3) = 0
      integer, allocatable :: reps(:, :)
      logical, allocatable :: seen(:), reverse_pulled(:)
      real(wp), allocatable :: kfrac(:, :), weights(:), vsh(:)
      real(wp), allocatable :: g_onsfx(:, :, :), g_onsri(:, :), &
         & dgdq_onsfx(:, :, :), dgdq_onsri(:, :)
      complex(wp), allocatable :: density(:, :, :, :), overlap(:, :, :), &
         & fock(:, :, :, :), overlap_adjoint(:, :, :)
      complex(wp), allocatable :: amat_r(:, :, :, :), cmat_r(:, :, :, :), &
         & vmat_r(:, :, :, :), gdiagP(:, :), gdiagSP(:, :), gdiagSPS(:, :)
      real(wp) :: energy = 0.0_wp
   end type cp2k_exchange_stream_type

contains

!> Return a real-space cutoff that covers both the valence basis and the
!> auxiliary ACP projectors used by g-xTB.
pure function cp2k_get_acp_cutoff(calc, accuracy) result(cutoff)
   type(xtb_calculator), intent(in) :: calc
   real(wp), intent(in) :: accuracy
   real(wp) :: cutoff

   cutoff = get_cutoff(calc%bas, accuracy)
   if (allocated(calc%acp)) then
      cutoff = max(cutoff, get_pair_cutoff(calc%bas, calc%acp%auxbas, accuracy))
   end if
end function cp2k_get_acp_cutoff

!> Construct a CP2K-facing g-xTB potential mixer.
subroutine new_cp2k_potential_mixer(self, ndim, error)
   type(cp2k_potential_mixer_type), intent(out) :: self
   integer, intent(in) :: ndim
   type(error_type), allocatable, intent(out) :: error

   call self%reset(ndim, error)
end subroutine new_cp2k_potential_mixer

!> Reset a CP2K-facing g-xTB potential mixer for a flat Fock vector.
subroutine cp2k_potential_mixer_reset(self, ndim, error)
   class(cp2k_potential_mixer_type), intent(inout) :: self
   integer, intent(in) :: ndim
   type(error_type), allocatable, intent(out) :: error

   if (allocated(self%simple_cache%raw)) then
      call self%simple%cleanup(self%simple_cache)
      deallocate(self%simple_cache%raw)
   end if
   if (allocated(self%diis_cache%raw)) then
      call self%diis%cleanup(self%diis_cache)
      deallocate(self%diis_cache%raw)
   end if

   self%ndim = 0
   self%iteration = 0
   self%residual = 0.0_wp
   self%initialized = .false.

   if (ndim <= 0) then
      call fatal_error(error, "CP2K potential mixer dimension must be positive")
      return
   end if

   call new_simple(self%simple, simple_input(mode=mixer_mode%potential, &
      & start=1, damp=0.2_wp))
   call new_diis(self%diis, diis_input(mode=mixer_mode%potential, start=4, &
      & precollect=1, memory=7, damp=1.0_wp, output_fraction=1.0_wp))

   self%ndim = ndim
   self%initialized = .true.
end subroutine cp2k_potential_mixer_reset

!> Advance g-xTB's native potential mixer from a raw to a mixed Fock vector.
subroutine cp2k_potential_mixer_advance(self, raw_fock, mixed_fock, error)
   class(cp2k_potential_mixer_type), intent(inout) :: self
   real(wp), intent(in) :: raw_fock(:)
   real(wp), intent(out) :: mixed_fock(:)
   type(error_type), allocatable, intent(out) :: error

   integer :: next_iteration

   mixed_fock = 0.0_wp
   if (.not. self%initialized) then
      call fatal_error(error, "CP2K potential mixer must be reset before use")
      return
   end if
   if (size(raw_fock) /= self%ndim .or. size(mixed_fock) /= self%ndim) then
      call fatal_error(error, "CP2K potential mixer Fock vector dimension mismatch")
      return
   end if
   if (.not. all(ieee_is_finite(raw_fock))) then
      call fatal_error(error, "CP2K potential mixer received a non-finite Fock vector")
      return
   end if

   next_iteration = self%iteration + 1
   select case(next_iteration)
   case(1)
      call self%simple%update(self%simple_cache, self%ndim)
      self%simple_cache%raw%iset = 0
      call self%simple%set(self%simple_cache, raw_fock)
      self%simple_cache%raw%initialized = .true.
      mixed_fock = raw_fock
      self%residual = 0.0_wp

   case(2)
      ! The DIIS cache is prepared one step before its first precollection.
      call self%diis%update(self%diis_cache, self%ndim)

      self%simple_cache%raw%idif = 0
      call self%simple%diff(self%simple_cache, raw_fock)
      self%residual = self%simple%get_error(self%simple_cache)
      call self%simple%collect(self%simple_cache)
      call self%simple%next(self%simple_cache, error)
      if (allocated(error)) return
      self%simple_cache%raw%iget = 0
      call self%simple%get(self%simple_cache, mixed_fock)

      call set_simple_input(self, mixed_fock)
      call set_diis_input(self, mixed_fock)

   case(3)
      self%simple_cache%raw%idif = 0
      call self%simple%diff(self%simple_cache, raw_fock)
      self%diis_cache%raw%idif = 0
      call self%diis%diff(self%diis_cache, raw_fock)
      self%residual = self%simple%get_error(self%simple_cache)

      call self%simple%collect(self%simple_cache)
      call self%diis%collect(self%diis_cache)
      call self%simple%next(self%simple_cache, error)
      if (allocated(error)) return
      self%simple_cache%raw%iget = 0
      call self%simple%get(self%simple_cache, mixed_fock)

      call set_simple_input(self, mixed_fock)
      call set_diis_input(self, mixed_fock)

   case default
      if (next_iteration == 4) then
         ! This is the iterator's simple-to-DIIS transition point.
         call self%simple%cleanup(self%simple_cache)
         deallocate(self%simple_cache%raw)
         call self%diis%update(self%diis_cache, self%ndim)
      end if

      self%diis_cache%raw%idif = 0
      call self%diis%diff(self%diis_cache, raw_fock)
      self%residual = self%diis%get_error(self%diis_cache)
      call self%diis%collect(self%diis_cache)
      call self%diis%next(self%diis_cache, error)
      if (allocated(error)) return
      self%diis_cache%raw%iget = 0
      call self%diis%get(self%diis_cache, mixed_fock)
      call set_diis_input(self, mixed_fock)
   end select

   if (.not. ieee_is_finite(self%residual)) then
      call fatal_error(error, "CP2K potential mixer produced a non-finite residual")
      return
   end if
   if (.not. all(ieee_is_finite(mixed_fock))) then
      call fatal_error(error, "CP2K potential mixer produced a non-finite Fock vector")
      return
   end if
   self%iteration = next_iteration
end subroutine cp2k_potential_mixer_advance

!> Return the RMS residual from the most recent raw Fock vector.
pure function cp2k_potential_mixer_get_error(self) result(residual)
   class(cp2k_potential_mixer_type), intent(in) :: self
   real(wp) :: residual

   residual = self%residual
end function cp2k_potential_mixer_get_error

!> Store the current extrapolated vector in the simple mixer cache.
subroutine set_simple_input(self, fock)
   class(cp2k_potential_mixer_type), intent(inout) :: self
   real(wp), intent(in) :: fock(:)

   self%simple_cache%raw%iset = 0
   call self%simple%set(self%simple_cache, fock)
   self%simple_cache%raw%initialized = .true.
end subroutine set_simple_input

!> Store the current extrapolated vector in the DIIS mixer cache.
subroutine set_diis_input(self, fock)
   class(cp2k_potential_mixer_type), intent(inout) :: self
   real(wp), intent(in) :: fock(:)

   self%diis_cache%raw%iset = 0
   call self%diis%set(self%diis_cache, fock)
   self%diis_cache%raw%initialized = .true.
end subroutine set_diis_input

!> Return the CP2K-owned SCC mixer dimension without the spin multiplicity.
!>
!> The orbital density is deliberately not included here, even when a
!> contribution advertises an orbital-resolved density potential.  CP2K
!> diagonalizes and stores that density itself; only charges and multipoles
!> are exchanged with its SCC mixer through the compatibility interface.
pure function cp2k_get_scc_mixer_dimension(mol, bas, info) result(ndim)
   type(structure_type), intent(in) :: mol
   class(basis_type), intent(in) :: bas
   type(scf_info), intent(in) :: info
   integer :: ndim

   ndim = 0
   select case(info%charge)
   case(atom_resolved)
      ndim = ndim + mol%nat
   case(shell_resolved)
      ndim = ndim + bas%nsh
   end select
   if (info%dipole == atom_resolved) ndim = ndim + 3*mol%nat
   if (info%quadrupole == atom_resolved) ndim = ndim + 6*mol%nat
end function cp2k_get_scc_mixer_dimension

!> Backward-compatible name for the CP2K SCC mixer dimension.
pure function cp2k_get_mixer_dimension(mol, bas, info) result(ndim)
   type(structure_type), intent(in) :: mol
   class(basis_type), intent(in) :: bas
   type(scf_info), intent(in) :: info
   integer :: ndim

   ndim = cp2k_get_scc_mixer_dimension(mol, bas, info)
end function cp2k_get_mixer_dimension

!> Apply SCC settings using the iterator API used by save_tblite.
subroutine cp2k_set_scf_options(calc, max_iter, damping)
   type(xtb_calculator), intent(inout) :: calc
   integer, intent(in) :: max_iter
   real(wp), intent(in) :: damping

   calc%iterator%max_iter = max_iter
   call calc%iterator%set_mixer_damping(simple_damping=damping, &
      & broyden_damping=damping)
end subroutine cp2k_set_scf_options

!> Construct the immutable auxiliary charges and cache the effective basis.
subroutine cp2k_update_basis(calc, mol, cache, error, grad, wfn_aux)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(basis_cache), intent(inout) :: cache
   type(error_type), allocatable, intent(out) :: error
   logical, intent(in), optional :: grad
   type(wavefunction_type), intent(inout), optional :: wfn_aux
   type(wavefunction_type) :: local_wfn_aux
   logical :: do_grad

   do_grad = .false.
   if (present(grad)) do_grad = grad

   if (allocated(calc%charge_model)) then
      if (present(wfn_aux)) then
         call new_wavefunction(wfn_aux, mol%nat, calc%bas%nsh, 0, 1, 0.0_wp, do_grad)
         if (do_grad) then
            call get_charges(calc%charge_model, mol, error, wfn_aux%qat(:, 1), &
               & dqdr=wfn_aux%dqatdr(:, :, :, 1), dqdL=wfn_aux%dqatdL(:, :, :, 1))
         else
            call get_charges(calc%charge_model, mol, error, wfn_aux%qat(:, 1))
         end if
         if (allocated(error)) return
         call calc%bas%update(mol, cache, do_grad, wfn_aux)
      else
         call new_wavefunction(local_wfn_aux, mol%nat, calc%bas%nsh, 0, 1, 0.0_wp, do_grad)
         if (do_grad) then
            call get_charges(calc%charge_model, mol, error, local_wfn_aux%qat(:, 1), &
               & dqdr=local_wfn_aux%dqatdr(:, :, :, 1), &
               & dqdL=local_wfn_aux%dqatdL(:, :, :, 1))
         else
            call get_charges(calc%charge_model, mol, error, local_wfn_aux%qat(:, 1))
         end if
         if (allocated(error)) return
         call calc%bas%update(mol, cache, do_grad, local_wfn_aux)
      end if
   else
      call calc%bas%update(mol, cache, do_grad)
   end if
end subroutine cp2k_update_basis

!> Build the complete Gamma-point overlap and one-electron g-xTB Hamiltonian.
!>
!> The dense matrices are deliberately returned in tblite's atom/AO ordering.
!> CP2K can scatter them into its block-sparse matrices while retaining the
!> element-wise basis metadata used for block dimensions and cutoffs.
subroutine cp2k_build_gamma_h0(calc, mol, bcache, acache, accuracy, selfenergy, &
   & overlap, dipole, quadrupole, hamiltonian)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(basis_cache), intent(in) :: bcache
   type(acp_cache), intent(inout) :: acache
   real(wp), intent(in) :: accuracy
   real(wp), intent(in) :: selfenergy(:)
   real(wp), intent(out) :: overlap(:, :), dipole(:, :, :), quadrupole(:, :, :)
   real(wp), intent(out) :: hamiltonian(:, :)

   real(wp) :: cutoff
   real(wp), allocatable :: aniso_dip(:, :), lattr(:, :)
   type(adjacency_list) :: list

   allocate(aniso_dip(3, mol%nat))
   cutoff = cp2k_get_acp_cutoff(calc, accuracy)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)
   if (any(mol%periodic)) then
      call get_anisotropy(calc%h0, mol, lattr, list, aniso_dip)
   else
      call get_anisotropy(calc%h0, mol, aniso_dip)
   end if
   call get_hamiltonian(mol, lattr, list, calc%bas, bcache, calc%h0, &
      & selfenergy, aniso_dip, overlap, dipole, quadrupole, hamiltonian)

   if (allocated(calc%acp)) then
      call calc%acp%update(mol, acache)
      call get_acp(mol, lattr, list, calc%bas, bcache, calc%acp, acache, hamiltonian)
   end if
end subroutine cp2k_build_gamma_h0

!> Build image-resolved q-vSZP overlap, multipoles, and H0 matrices.
!>
!> This routine deliberately covers the one-electron H0 only.  The atomic
!> correction potential is a projector convolution and therefore requires a
!> separate image-resolved implementation before a complete k-point
!> Hamiltonian can be assembled.
subroutine cp2k_build_image_h0(calc, mol, bcache, accuracy, selfenergy, &
   & translations, images, overlap_gamma, dipole_gamma, quadrupole_gamma, &
   & hamiltonian_gamma)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(basis_cache), intent(in) :: bcache
   real(wp), intent(in) :: accuracy
   real(wp), intent(in) :: selfenergy(:)
   real(wp), intent(in) :: translations(:, :)
   type(hamiltonian_image_type), intent(out) :: images
   real(wp), intent(out), optional :: overlap_gamma(:, :)
   real(wp), intent(out), optional :: dipole_gamma(:, :, :)
   real(wp), intent(out), optional :: quadrupole_gamma(:, :, :)
   real(wp), intent(out), optional :: hamiltonian_gamma(:, :)

   real(wp) :: cutoff
   real(wp), allocatable :: aniso_dip(:, :), dipole(:, :, :), &
      & hamiltonian(:, :), overlap(:, :), quadrupole(:, :, :)
   type(adjacency_list) :: list

   allocate(aniso_dip(3, mol%nat), overlap(calc%bas%nao, calc%bas%nao), &
      & dipole(3, calc%bas%nao, calc%bas%nao), &
      & quadrupole(6, calc%bas%nao, calc%bas%nao), &
      & hamiltonian(calc%bas%nao, calc%bas%nao))
   cutoff = get_cutoff(calc%bas, accuracy)
   call new_adjacency_list(list, mol, translations, cutoff)
   if (any(mol%periodic)) then
      call get_anisotropy(calc%h0, mol, translations, list, aniso_dip)
   else
      call get_anisotropy(calc%h0, mol, aniso_dip)
   end if
   call get_hamiltonian(mol, translations, list, calc%bas, bcache, calc%h0, &
      & selfenergy, aniso_dip, overlap, dipole, quadrupole, hamiltonian, images)

   if (present(overlap_gamma)) overlap_gamma = overlap
   if (present(dipole_gamma)) dipole_gamma = dipole
   if (present(quadrupole_gamma)) quadrupole_gamma = quadrupole
   if (present(hamiltonian_gamma)) hamiltonian_gamma = hamiltonian
end subroutine cp2k_build_image_h0

!> Build the image-resolved separable ACP contribution.
!>
!> The returned Hamiltonian translations are generally a larger difference
!> set than the translations supplied for the projector integrals.  External
!> drivers must therefore scatter this set independently of the pair-local H0
!> image set.
subroutine cp2k_build_image_acp(calc, mol, bcache, acache, translations, &
   & images, hamiltonian_gamma)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(basis_cache), intent(in) :: bcache
   type(acp_cache), intent(inout) :: acache
   real(wp), intent(in) :: translations(:, :)
   type(acp_image_type), intent(out) :: images
   real(wp), intent(out), optional :: hamiltonian_gamma(:, :)

   call get_acp_images(mol, translations, calc%bas, bcache, calc%acp, acache, &
      & images)
   if (present(hamiltonian_gamma)) then
      hamiltonian_gamma = sum(images%hamiltonian, dim=3)
   end if
end subroutine cp2k_build_image_acp

!> Build only the compact projector-overlap images required by the direct
!> Bloch-space ACP forward and reverse sweeps.  Unlike
!> `cp2k_build_image_acp`, this deliberately does not materialize the
!> quadratic translation-difference set or H(R).
subroutine cp2k_build_acp_projectors(calc, mol, bcache, acache, translations, &
   & images)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(basis_cache), intent(in) :: bcache
   type(acp_cache), intent(inout) :: acache
   real(wp), intent(in) :: translations(:, :)
   type(acp_image_type), intent(out) :: images

   call get_acp_images(mol, translations, calc%bas, bcache, calc%acp, acache, &
      & images, build_hamiltonian=.false.)
end subroutine cp2k_build_acp_projectors

!> Image-resolved ACP response for CP2K k-point density matrices.
!>
!> `density(:, :, R, spin)` must be in save_tblite's AO image orientation and
!> dual to the Hamiltonian returned by `cp2k_build_image_acp`, so that its
!> energy contribution is the element-wise contraction over equal images.
!> CP2K's irreducible-star weights must already be included in these real-space
!> densities.
subroutine cp2k_acp_image_gradient(calc, mol, bcache, acache, images, &
   & density_translation, density, dEdcnbas, dEdqbas, gradient, sigma, error)
   type(xtb_calculator), intent(in) :: calc
   type(structure_type), intent(in) :: mol
   type(basis_cache), intent(in) :: bcache
   type(acp_cache), intent(in) :: acache
   type(acp_image_type), intent(in) :: images
   real(wp), intent(in) :: density_translation(:, :)
   real(wp), intent(in) :: density(:, :, :, :)
   real(wp), intent(inout) :: dEdcnbas(:), dEdqbas(:)
   real(wp), intent(inout) :: gradient(:, :), sigma(:, :)
   type(error_type), allocatable, intent(out) :: error

   if (.not. allocated(calc%acp)) then
      call fatal_error(error, "g-xTB ACP image gradient requires an ACP model")
      return
   end if
   call get_acp_image_gradient(mol, calc%bas, bcache, calc%acp, acache, &
      & images, density_translation, density, dEdcnbas, dEdqbas, gradient, &
      & sigma, error)
end subroutine cp2k_acp_image_gradient

!> Evaluate the separable ACP directly on CP2K Bloch points.
subroutine cp2k_acp_kmesh(calc, mol, images, kfrac, hamiltonian, error)
   type(xtb_calculator), intent(in) :: calc
   type(structure_type), intent(in) :: mol
   type(acp_image_type), intent(in) :: images
   real(wp), intent(in) :: kfrac(:, :)
   complex(wp), intent(out) :: hamiltonian(:, :, :)
   type(error_type), allocatable, intent(out) :: error

   if (.not.allocated(calc%acp)) then
      call fatal_error(error, "g-xTB ACP container is not allocated")
      return
   end if
   call get_acp_kmesh(mol, images, kfrac, hamiltonian, error)
end subroutine cp2k_acp_kmesh

!> Differentiate the complete weighted Bloch contraction of the ACP.
subroutine cp2k_acp_kmesh_gradient(calc, mol, bcache, acache, images, &
   & kfrac, weights, density, dEdcnbas, dEdqbas, gradient, sigma, error)
   type(xtb_calculator), intent(in) :: calc
   type(structure_type), intent(in) :: mol
   type(basis_cache), intent(in) :: bcache
   type(acp_cache), intent(in) :: acache
   type(acp_image_type), intent(in) :: images
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :)
   real(wp), intent(inout) :: dEdcnbas(:), dEdqbas(:)
   real(wp), intent(inout) :: gradient(:, :), sigma(:, :)
   type(error_type), allocatable, intent(out) :: error

   if (.not.allocated(calc%acp)) then
      call fatal_error(error, "g-xTB ACP container is not allocated")
      return
   end if
   call get_acp_kmesh_gradient(mol, calc%bas, bcache, calc%acp, acache, &
      & images, kfrac, weights, density, dEdcnbas, dEdqbas, gradient, &
      & sigma, error)
end subroutine cp2k_acp_kmesh_gradient

!> Apply a component's q-vSZP basis response and external-charge chain rule.
!>
!> Image-resolved one-electron components return their own additive
!> `dE/dCN_basis` and `dE/dq_basis`.  This helper maps exactly that component
!> onto Cartesian and strain derivatives.  Calling it separately for H0 and
!> ACP therefore remains additive and avoids applying either response twice.
subroutine cp2k_apply_basis_gradient(calc, mol, wfn_aux, dEdcnbas, dEdqbas, &
   & gradient, sigma)
   type(xtb_calculator), intent(in) :: calc
   type(structure_type), intent(in) :: mol
   type(wavefunction_type), intent(in) :: wfn_aux
   real(wp), intent(in) :: dEdcnbas(:), dEdqbas(:)
   real(wp), intent(inout) :: gradient(:, :), sigma(:, :)

   real(wp), allocatable :: dEdq(:)

   allocate(dEdq(mol%nat), source=0.0_wp)
   if (calc%bas%charge_dependent) then
      call calc%bas%get_basis_gradient(mol, dEdcnbas, dEdqbas, dEdq, &
         & gradient, sigma)
   end if
   if (allocated(wfn_aux%dqatdr)) then
      call gemv(wfn_aux%dqatdr(:, :, :, 1), dEdq, gradient, beta=1.0_wp)
   end if
   if (allocated(wfn_aux%dqatdL)) then
      call gemv(wfn_aux%dqatdL(:, :, :, 1), dEdq, sigma, beta=1.0_wp)
   end if
end subroutine cp2k_apply_basis_gradient

!> Prepare geometry- and charge-dependent exchange kernels.
!>
!> This routine must be called once, serially, before evaluating the k-point
!> blocks of every SCC iteration.  Rebuilding here also makes changes of the
!> geometry or lattice visible to the Wigner-Seitz and bond-order kernels.
subroutine cp2k_prepare_exchange(calc, mol, ecache, wfn, error)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: ecache
   type(wavefunction_type), intent(in) :: wfn
   type(error_type), allocatable, intent(out) :: error

   if (.not. allocated(calc%exchange)) then
      call fatal_error(error, "g-xTB exchange container is not allocated")
      return
   end if

   call calc%exchange%update(mol, ecache)
   select type (cache => ecache%raw)
   type is (exchange_cache)
      call calc%exchange%get_onsite_Kmatrix(mol, wfn, cache)
   class default
      call fatal_error(error, "unexpected g-xTB exchange cache type")
   end select
end subroutine cp2k_prepare_exchange

!> Evaluate the g-xTB exchange operator for one complex k-point block.
!>
!> The input P(k) must not contain a k-point or star weight.  The caller owns
!> all Brillouin-zone weights and symmetry expansion.  This routine returns
!> an unweighted Hermitian F(k), an atom/shell-resolved potential, and an
!> unweighted exchange energy.  For a reduced k mesh the shell potential must
!> be permuted with every star operation before the weighted star sum is made.
!> Call cp2k_prepare_exchange serially before entering the k-point loop; this
!> evaluator then only reads the shared cache and can be called concurrently.
subroutine cp2k_exchange_kpoint(calc, mol, ecache, density, overlap, &
   & fock, vsh, energy, error)
   type(xtb_calculator), intent(in) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(in) :: ecache
   complex(wp), intent(in) :: density(:, :, :), overlap(:, :)
   complex(wp), intent(out) :: fock(:, :, :)
   real(wp), intent(out) :: vsh(:), energy
   type(error_type), allocatable, intent(out) :: error

   if (.not. allocated(calc%exchange)) then
      call fatal_error(error, "g-xTB exchange container is not allocated")
      return
   end if
   if (.not. allocated(ecache%raw)) then
      call fatal_error(error, &
         & "k-point exchange cache is not prepared; call cp2k_prepare_exchange")
      return
   end if

   select type (cache => ecache%raw)
   type is (exchange_cache)
      select type (exchange => calc%exchange)
      type is (exchange_fock)
         call exchange%get_KFock_kpoint(mol, cache, density, overlap, fock, &
            & vsh, energy)
      class default
         call fatal_error(error, &
            & "k-point exchange requires the exchange_fock implementation")
      end select
   class default
      call fatal_error(error, "unexpected g-xTB exchange cache type")
   end select
end subroutine cp2k_exchange_kpoint

!> Build or reuse the geometry- and mesh-static BvK kernel/Fourier plan.
subroutine prepare_bvk_plan(exchange, mol, cache, nmesh, kfrac, weights, error)
   type(exchange_fock), intent(in) :: exchange
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(inout) :: cache
   integer, intent(in) :: nmesh(3)
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   type(error_type), allocatable, intent(out) :: error

   integer :: grid(3), icell, idim, igrid, ik, jcell, nearest, nk
   logical, allocatable :: seen(:)
   real(wp) :: angle, scaled, target
   real(wp), allocatable :: kred(:, :)
   complex(wp) :: orthogonality, phase

   nk = size(weights)
   if (cache%bvk_matches(mol, nmesh, kfrac, weights) &
      & .and. exchange%bvk_model_matches(cache)) then
      if (any(shape(cache%bvk_kernel%g_mulliken_r) /= &
         & [exchange%nsh, exchange%nsh, nk])) then
         cache%bvk_plan_valid = .false.
      else
         return
      end if
   end if
   cache%bvk_plan_valid = .false.

   if (allocated(cache%bvk_input_to_grid)) &
      & deallocate(cache%bvk_input_to_grid)
   if (allocated(cache%bvk_grid_to_input)) &
      & deallocate(cache%bvk_grid_to_input)
   allocate(cache%bvk_input_to_grid(nk), cache%bvk_grid_to_input(nk), source=0)
   allocate(seen(nk), source=.false.)
   allocate(kred(3, nk))
   kred = modulo(kfrac, 1.0_wp)
   do idim = 1, 3
      cache%bvk_twist(idim) = kred(idim, 1) &
         & -real(floor(kred(idim, 1)*real(nmesh(idim), wp)), wp) &
         & /real(nmesh(idim), wp)
   end do
   do ik = 1, nk
      do idim = 1, 3
         scaled = (kred(idim, ik)-cache%bvk_twist(idim)) &
            & *real(nmesh(idim), wp)
         nearest = nint(scaled)
         if (abs(scaled-real(nearest, wp)) > 1.0e-10_wp) then
            call fatal_error(error, &
               & "g-xTB k points do not share a regular-grid twist")
            return
         end if
         grid(idim) = modulo(nearest, nmesh(idim))
      end do
      igrid = 1 + grid(1) + nmesh(1)*(grid(2)+nmesh(2)*grid(3))
      if (seen(igrid)) then
         call fatal_error(error, "g-xTB k-point mesh contains a duplicate grid point")
         return
      end if
      seen(igrid) = .true.
      cache%bvk_input_to_grid(ik) = igrid
      cache%bvk_grid_to_input(igrid) = ik
   end do
   if (.not.all(seen)) then
      call fatal_error(error, "g-xTB k-point mesh is incomplete")
      return
   end if

   ! Construct the expensive image kernel only after the O(Nk) mesh proof.
   call exchange%get_bvk_Kmatrix(mol, nmesh, cache%bvk_kernel)
   if (any(shape(cache%bvk_kernel%reps) /= [3, nk]) &
      & .or. any(shape(cache%bvk_kernel%g_mulliken_r) /= &
      & [exchange%nsh, exchange%nsh, nk]) &
      & .or. any(shape(cache%bvk_kernel%g_bocorr_r) /= &
      & [mol%nat, mol%nat, nk])) then
      call fatal_error(error, "g-xTB BvK kernel has inconsistent dimensions")
      return
   end if

   if (allocated(cache%bvk_phase_forward)) &
      & deallocate(cache%bvk_phase_forward)
   if (allocated(cache%bvk_phase_inverse)) &
      & deallocate(cache%bvk_phase_inverse)
   allocate(cache%bvk_phase_forward(size(cache%bvk_kernel%reps, 2), nk), &
      & cache%bvk_phase_inverse(nk, size(cache%bvk_kernel%reps, 2)))
   do ik = 1, nk
      do icell = 1, size(cache%bvk_kernel%reps, 2)
         angle = 2.0_wp*pi*dot_product(kred(:, ik), &
            & real(cache%bvk_kernel%reps(:, icell), wp))
         phase = exp(cmplx(0.0_wp, angle, wp))
         cache%bvk_phase_forward(icell, ik) = phase
         cache%bvk_phase_inverse(ik, icell) = weights(ik)*conjg(phase)
      end do
   end do

   ! Keep a dense phase-orthogonality oracle for small meshes.  The unique
   ! integer-grid map above is the O(Nk) production proof for larger meshes
   ! and is also the ordering descriptor required by a future FFT backend.
   if (nk <= 64) then
      do icell = 1, nk
         do jcell = 1, nk
            orthogonality = sum(cache%bvk_phase_forward(icell, :) &
               & *cache%bvk_phase_inverse(:, jcell))
            target = merge(1.0_wp, 0.0_wp, icell == jcell)
            if (abs(orthogonality-target) > 1.0e-10_wp) then
               call fatal_error(error, &
                  & "g-xTB k points are not a complete Fourier mesh")
               return
            end if
         end do
      end do
   end if
   call exchange%set_bvk_model_signature(cache)
   call cache%set_bvk_signature(mol, nmesh, kfrac, weights)
end subroutine prepare_bvk_plan

!> Evaluate BvK-supercell-equivalent exchange on a complete regular k mesh.
!>
!> Density, overlap, and Fock contain unweighted Bloch blocks.  Energy and
!> shell potential are returned as complete Brillouin-zone contractions per
!> primitive cell and must not be weighted a second time by CP2K.
subroutine cp2k_exchange_kmesh(calc, mol, ecache, wfn, nmesh, kfrac, weights, &
   & density, overlap, fock, vsh, energy, error)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: ecache
   type(wavefunction_type), intent(in) :: wfn
   integer, intent(in) :: nmesh(3)
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   complex(wp), intent(out) :: fock(:, :, :, :)
   real(wp), intent(out) :: vsh(:), energy
   type(error_type), allocatable, intent(out) :: error

   integer :: idim, nk

   nk = size(weights)
   if (any(nmesh < 1) .or. product(nmesh) /= nk) then
      call fatal_error(error, "g-xTB exchange requires a complete regular k mesh")
      return
   end if
   do idim = 1, 3
      if (.not.mol%periodic(idim) .and. nmesh(idim) /= 1) then
         call fatal_error(error, &
            & "nonperiodic directions cannot be replicated by the g-xTB k mesh")
         return
      end if
   end do
   if (size(kfrac, 1) /= 3 .or. size(kfrac, 2) /= nk &
      & .or. size(overlap, 3) /= nk .or. size(density, 4) /= nk &
      & .or. size(fock, 4) /= nk) then
      call fatal_error(error, "inconsistent g-xTB whole-mesh array dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(kfrac))) then
      call fatal_error(error, "g-xTB whole-mesh k points must be finite")
      return
   end if
   if (.not.all(ieee_is_finite(weights)) .or. &
      & abs(sum(weights)-1.0_wp) > 1.0e-12_wp &
      & .or. maxval(abs(weights-1.0_wp/real(nk, wp))) > 1.0e-12_wp) then
      call fatal_error(error, "g-xTB whole-mesh exchange requires uniform weights")
      return
   end if

   call cp2k_prepare_exchange(calc, mol, ecache, wfn, error)
   if (allocated(error)) return
   select type (cache => ecache%raw)
   type is (exchange_cache)
      select type (exchange => calc%exchange)
      type is (exchange_fock)
         if (size(overlap, 1) /= exchange%nao &
            & .or. size(overlap, 2) /= exchange%nao &
            & .or. size(density, 1) /= exchange%nao &
            & .or. size(density, 2) /= exchange%nao &
            & .or. size(fock, 1) /= exchange%nao &
            & .or. size(fock, 2) /= exchange%nao &
            & .or. (size(density, 3) /= 1 .and. size(density, 3) /= 2) &
            & .or. size(density, 3) /= wfn%nspin &
            & .or. size(fock, 3) /= size(density, 3) &
            & .or. size(vsh) /= exchange%nsh) then
            call fatal_error(error, &
               & "inconsistent g-xTB whole-mesh AO, spin, or shell dimensions")
            return
         end if
         call prepare_bvk_plan(exchange, mol, cache, nmesh, kfrac, weights, error)
         if (allocated(error)) return
         call exchange%get_KFock_kmesh(mol, cache, cache%bvk_kernel, &
            & kfrac, weights, density, overlap, fock, vsh, energy)
      class default
         call fatal_error(error, &
            & "whole-mesh exchange requires the exchange_fock implementation")
      end select
   class default
      call fatal_error(error, "unexpected g-xTB exchange cache type")
   end select

end subroutine cp2k_exchange_kmesh


!> Release all storage owned by a Brillouin-zone exchange stream.
subroutine clear_exchange_stream(stream)
   type(cp2k_exchange_stream_type), intent(inout) :: stream

   if (allocated(stream%seen)) deallocate(stream%seen)
   if (allocated(stream%reverse_pulled)) deallocate(stream%reverse_pulled)
   if (allocated(stream%reps)) deallocate(stream%reps)
   if (allocated(stream%kfrac)) deallocate(stream%kfrac)
   if (allocated(stream%weights)) deallocate(stream%weights)
   if (allocated(stream%vsh)) deallocate(stream%vsh)
   if (allocated(stream%g_onsfx)) deallocate(stream%g_onsfx)
   if (allocated(stream%g_onsri)) deallocate(stream%g_onsri)
   if (allocated(stream%dgdq_onsfx)) deallocate(stream%dgdq_onsfx)
   if (allocated(stream%dgdq_onsri)) deallocate(stream%dgdq_onsri)
   if (allocated(stream%density)) deallocate(stream%density)
   if (allocated(stream%overlap)) deallocate(stream%overlap)
   if (allocated(stream%fock)) deallocate(stream%fock)
   if (allocated(stream%overlap_adjoint)) deallocate(stream%overlap_adjoint)
   if (allocated(stream%amat_r)) deallocate(stream%amat_r)
   if (allocated(stream%cmat_r)) deallocate(stream%cmat_r)
   if (allocated(stream%vmat_r)) deallocate(stream%vmat_r)
   if (allocated(stream%gdiagP)) deallocate(stream%gdiagP)
   if (allocated(stream%gdiagSP)) deallocate(stream%gdiagSP)
   if (allocated(stream%gdiagSPS)) deallocate(stream%gdiagSPS)
   stream%active = .false.
   stream%applied = .false.
   stream%reverse_applied = .false.
   stream%nk = 0
   stream%nao = 0
   stream%nspin = 0
   stream%nsh = 0
   stream%mode = 0
   stream%nmesh = 0
   stream%energy = 0.0_wp
end subroutine clear_exchange_stream


!> Report whether a stream owns complete density or overlap k-mesh storage.
logical function cp2k_exchange_stream_has_full_mesh_storage(stream) &
   & result(has_storage)
   type(cp2k_exchange_stream_type), intent(in) :: stream

   has_storage = allocated(stream%density) .or. allocated(stream%overlap)
end function cp2k_exchange_stream_has_full_mesh_storage


!> Check the exact charge-dependent onsite state captured by `begin`.
logical function exchange_stream_onsite_matches(stream, cache) result(matches)
   type(cp2k_exchange_stream_type), intent(in) :: stream
   type(exchange_cache), intent(in) :: cache

   matches = .false.
   if (.not.allocated(stream%g_onsfx) &
      & .or. .not.allocated(stream%g_onsri) &
      & .or. .not.allocated(stream%dgdq_onsfx) &
      & .or. .not.allocated(stream%dgdq_onsri) &
      & .or. .not.allocated(cache%g_onsfx) &
      & .or. .not.allocated(cache%g_onsri) &
      & .or. .not.allocated(cache%dgdq_onsfx) &
      & .or. .not.allocated(cache%dgdq_onsri)) return
   if (any(shape(stream%g_onsfx) /= shape(cache%g_onsfx)) &
      & .or. any(shape(stream%g_onsri) /= shape(cache%g_onsri)) &
      & .or. any(shape(stream%dgdq_onsfx) /= shape(cache%dgdq_onsfx)) &
      & .or. any(shape(stream%dgdq_onsri) /= shape(cache%dgdq_onsri))) return
   if (any(stream%g_onsfx /= cache%g_onsfx) &
      & .or. any(stream%g_onsri /= cache%g_onsri) &
      & .or. any(stream%dgdq_onsfx /= cache%dgdq_onsfx) &
      & .or. any(stream%dgdq_onsri /= cache%dgdq_onsri)) return
   matches = .true.
end function exchange_stream_onsite_matches


!> Begin collecting the unweighted Bloch blocks of one complete regular mesh.
!>
!> Only one transaction may be active in a stream object.  The expensive BvK
!> plan is prepared here, before blocks are accepted, so that `apply` can also
!> detect an intervening geometry, mesh, or model change.  The optional mode
!> selects the memory-reduced production path or the full-mesh oracle; reduced
!> mode is the default.
subroutine cp2k_exchange_stream_begin(stream, calc, mol, ecache, wfn, nmesh, &
   & kfrac, weights, error, mode)
   type(cp2k_exchange_stream_type), intent(inout) :: stream
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: ecache
   type(wavefunction_type), intent(in) :: wfn
   integer, intent(in) :: nmesh(3)
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   type(error_type), allocatable, intent(out) :: error
   integer, intent(in), optional :: mode

   integer :: idim, nk, selected_mode

   if (stream%active) then
      call fatal_error(error, "g-xTB exchange stream is already active")
      return
   end if

   selected_mode = cp2k_exchange_stream_reduced
   if (present(mode)) selected_mode = mode
   if (selected_mode /= cp2k_exchange_stream_reduced &
      & .and. selected_mode /= cp2k_exchange_stream_oracle) then
      call fatal_error(error, "unknown g-xTB exchange stream mode")
      return
   end if

   nk = size(weights)
   if (any(nmesh < 1) .or. product(nmesh) /= nk) then
      call fatal_error(error, &
         & "g-xTB exchange stream requires a complete regular k mesh")
      return
   end if
   do idim = 1, 3
      if (.not.mol%periodic(idim) .and. nmesh(idim) /= 1) then
         call fatal_error(error, &
            & "nonperiodic directions cannot be replicated by the g-xTB k mesh")
         return
      end if
   end do
   if (size(kfrac, 1) /= 3 .or. size(kfrac, 2) /= nk) then
      call fatal_error(error, &
         & "inconsistent g-xTB exchange stream mesh dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(kfrac))) then
      call fatal_error(error, "g-xTB exchange stream k points must be finite")
      return
   end if
   if (.not.all(ieee_is_finite(weights)) .or. &
      & abs(sum(weights)-1.0_wp) > 1.0e-12_wp &
      & .or. maxval(abs(weights-1.0_wp/real(nk, wp))) > 1.0e-12_wp) then
      call fatal_error(error, &
         & "g-xTB exchange stream requires uniform k-point weights")
      return
   end if

   call cp2k_prepare_exchange(calc, mol, ecache, wfn, error)
   if (allocated(error)) return
   select type (cache => ecache%raw)
   type is (exchange_cache)
      select type (exchange => calc%exchange)
      type is (exchange_fock)
         if (wfn%nspin /= 1 .and. wfn%nspin /= 2) then
            call fatal_error(error, &
               & "g-xTB exchange stream requires one or two spin channels")
            return
         end if
         call prepare_bvk_plan(exchange, mol, cache, nmesh, kfrac, weights, error)
         if (allocated(error)) return
         if (.not.allocated(cache%g_onsfx) &
            & .or. .not.allocated(cache%g_onsri) &
            & .or. .not.allocated(cache%dgdq_onsfx) &
            & .or. .not.allocated(cache%dgdq_onsri)) then
            call fatal_error(error, &
               & "g-xTB exchange stream onsite state is not prepared")
            return
         end if

         call clear_exchange_stream(stream)
         stream%nk = nk
         stream%nao = exchange%nao
         stream%nspin = wfn%nspin
         stream%nsh = exchange%nsh
         stream%mode = selected_mode
         stream%nmesh = nmesh
         allocate(stream%seen(nk), source=.false.)
         allocate(stream%reverse_pulled(nk), source=.false.)
         allocate(stream%reps, source=cache%bvk_kernel%reps)
         allocate(stream%kfrac(3, nk), source=kfrac)
         allocate(stream%weights(nk), source=weights)
         allocate(stream%g_onsfx, source=cache%g_onsfx)
         allocate(stream%g_onsri, source=cache%g_onsri)
         allocate(stream%dgdq_onsfx, source=cache%dgdq_onsfx)
         allocate(stream%dgdq_onsri, source=cache%dgdq_onsri)
         if (selected_mode == cp2k_exchange_stream_oracle) then
            allocate(stream%density(exchange%nao, exchange%nao, wfn%nspin, nk), &
               & stream%overlap(exchange%nao, exchange%nao, nk), &
               & stream%fock(exchange%nao, exchange%nao, wfn%nspin, nk), &
               & source=(0.0_wp, 0.0_wp))
            allocate(stream%overlap_adjoint(exchange%nao, exchange%nao, nk), &
               & source=(0.0_wp, 0.0_wp))
         else
            allocate(stream%amat_r(exchange%nao, exchange%nao, nk, wfn%nspin), &
               & stream%cmat_r(exchange%nao, exchange%nao, nk, wfn%nspin), &
               & stream%vmat_r(exchange%nao, exchange%nao, nk, wfn%nspin), &
               & source=(0.0_wp, 0.0_wp))
            allocate(stream%gdiagP(exchange%nao, wfn%nspin), &
               & stream%gdiagSP(exchange%nao, wfn%nspin), &
               & stream%gdiagSPS(exchange%nao, wfn%nspin), &
               & source=(0.0_wp, 0.0_wp))
         end if
         allocate(stream%vsh(exchange%nsh), source=0.0_wp)
         stream%energy = 0.0_wp
         stream%active = .true.
      class default
         call fatal_error(error, &
            & "exchange stream requires the exchange_fock implementation")
      end select
   class default
      call fatal_error(error, "unexpected g-xTB exchange cache type")
   end select
end subroutine cp2k_exchange_stream_begin


!> Push one unweighted density/overlap Bloch block into an active stream.
subroutine cp2k_exchange_stream_push(stream, ik, density, overlap, error)
   type(cp2k_exchange_stream_type), intent(inout) :: stream
   integer, intent(in) :: ik
   complex(wp), intent(in) :: density(:, :, :), overlap(:, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: icell, spin
   real(wp) :: angle
   complex(wp) :: phase
   complex(wp), allocatable :: amat(:, :), cmat(:, :)

   if (.not.stream%active) then
      call fatal_error(error, "g-xTB exchange stream is not active")
      return
   end if
   if (stream%applied) then
      call fatal_error(error, "g-xTB exchange stream was already applied")
      return
   end if
   if (ik < 1 .or. ik > stream%nk) then
      call fatal_error(error, "g-xTB exchange stream block index is out of range")
      return
   end if
   if (stream%seen(ik)) then
      call fatal_error(error, "duplicate g-xTB exchange stream block")
      return
   end if
   if (size(density, 1) /= stream%nao &
      & .or. size(density, 2) /= stream%nao &
      & .or. size(density, 3) /= stream%nspin &
      & .or. size(overlap, 1) /= stream%nao &
      & .or. size(overlap, 2) /= stream%nao) then
      call fatal_error(error, &
         & "inconsistent g-xTB exchange stream block dimensions")
      return
   end if

   if (stream%mode == cp2k_exchange_stream_oracle) then
      stream%density(:, :, :, ik) = density
      stream%overlap(:, :, ik) = overlap
   else
      allocate(amat(stream%nao, stream%nao), cmat(stream%nao, stream%nao))
      do spin = 1, stream%nspin
         amat = matmul(overlap, density(:, :, spin))
         cmat = 0.5_wp*matmul(amat, overlap)
         do icell = 1, stream%nk
            angle = 2.0_wp*pi*dot_product(stream%kfrac(:, ik), &
               & real(stream%reps(:, icell), wp))
            phase = stream%weights(ik)*exp(cmplx(0.0_wp, -angle, wp))
            stream%amat_r(:, :, icell, spin) = &
               & stream%amat_r(:, :, icell, spin) + phase*amat
            stream%cmat_r(:, :, icell, spin) = &
               & stream%cmat_r(:, :, icell, spin) + phase*cmat
            stream%vmat_r(:, :, icell, spin) = &
               & stream%vmat_r(:, :, icell, spin) + phase*density(:, :, spin)
         end do
      end do
   end if
   stream%seen(ik) = .true.
end subroutine cp2k_exchange_stream_push


!> Apply exchange after all blocks have been collected.
!>
!> Reduced mode applies the kernels directly to the accumulated BvK images;
!> oracle mode dispatches the unchanged complete-mesh evaluator.
subroutine cp2k_exchange_stream_apply(stream, calc, mol, ecache, vsh, energy, &
   & error)
   type(cp2k_exchange_stream_type), intent(inout) :: stream
   type(xtb_calculator), intent(in) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: ecache
   real(wp), intent(out) :: vsh(:), energy
   type(error_type), allocatable, intent(out) :: error

   if (.not.stream%active) then
      call fatal_error(error, "g-xTB exchange stream is not active")
      return
   end if
   if (stream%applied) then
      call fatal_error(error, "g-xTB exchange stream was already applied")
      return
   end if
   if (.not.all(stream%seen)) then
      call fatal_error(error, "g-xTB exchange stream has missing blocks")
      return
   end if
   if (size(vsh) /= stream%nsh) then
      call fatal_error(error, &
         & "inconsistent g-xTB exchange stream shell-potential dimensions")
      return
   end if
   if (.not.allocated(calc%exchange) .or. .not.allocated(ecache%raw)) then
      call fatal_error(error, &
         & "g-xTB exchange stream state was invalidated before apply")
      return
   end if

   select type (cache => ecache%raw)
   type is (exchange_cache)
      select type (exchange => calc%exchange)
      type is (exchange_fock)
         if (exchange%nao /= stream%nao .or. exchange%nsh /= stream%nsh &
            & .or. .not.cache%bvk_matches(mol, stream%nmesh, &
            & stream%kfrac, stream%weights) &
            & .or. .not.exchange%bvk_model_matches(cache) &
            & .or. .not.exchange_stream_onsite_matches(stream, cache)) then
            call fatal_error(error, &
               & "g-xTB exchange stream state was invalidated before apply")
            return
         end if
         if (stream%mode == cp2k_exchange_stream_oracle) then
            call exchange%get_KFock_kmesh(mol, cache, cache%bvk_kernel, &
               & stream%kfrac, stream%weights, stream%density, stream%overlap, &
               & stream%fock, stream%vsh, stream%energy)
         else
            call exchange%get_KFock_stream_apply(mol, cache, &
               & cache%bvk_kernel, cache%bvk_phase_forward, stream%weights, &
               & stream%amat_r, stream%cmat_r, stream%vmat_r, &
               & stream%gdiagP, stream%gdiagSP, stream%gdiagSPS, &
               & stream%vsh, stream%energy)
         end if
         vsh = stream%vsh
         energy = stream%energy
         stream%applied = .true.
      class default
         call fatal_error(error, &
            & "exchange stream requires the exchange_fock implementation")
      end select
   class default
      call fatal_error(error, "unexpected g-xTB exchange cache type")
   end select
end subroutine cp2k_exchange_stream_apply


!> Retrieve one unweighted Fock block after a successful stream application.
!>
!> Reduced mode requires the corresponding overlap block again and
!> reconstructs only this requested k point.  Oracle mode reads the retained
!> complete-mesh result and does not require the optional overlap argument.
subroutine cp2k_exchange_stream_get(stream, ik, fock, error, overlap)
   type(cp2k_exchange_stream_type), intent(in) :: stream
   integer, intent(in) :: ik
   complex(wp), intent(out) :: fock(:, :, :)
   type(error_type), allocatable, intent(out) :: error
   complex(wp), intent(in), optional :: overlap(:, :)

   integer :: iao, icell, ii, jj, spin
   real(wp) :: angle, spin_factor
   complex(wp) :: phase, tmp
   complex(wp), allocatable :: amat(:, :), cmat(:, :), vmat(:, :), work(:, :)

   if (.not.stream%active) then
      call fatal_error(error, "g-xTB exchange stream is not active")
      return
   end if
   if (.not.stream%applied) then
      call fatal_error(error, "g-xTB exchange stream was not applied")
      return
   end if
   if (ik < 1 .or. ik > stream%nk) then
      call fatal_error(error, "g-xTB exchange stream block index is out of range")
      return
   end if
   if (size(fock, 1) /= stream%nao .or. size(fock, 2) /= stream%nao &
      & .or. size(fock, 3) /= stream%nspin) then
      call fatal_error(error, &
         & "inconsistent g-xTB exchange stream Fock dimensions")
      return
   end if

   if (stream%mode == cp2k_exchange_stream_oracle) then
      fock = stream%fock(:, :, :, ik)
      return
   end if
   if (.not.present(overlap)) then
      call fatal_error(error, &
         & "memory-reduced g-xTB exchange stream get requires overlap block")
      return
   end if
   if (size(overlap, 1) /= stream%nao .or. size(overlap, 2) /= stream%nao) then
      call fatal_error(error, &
         & "inconsistent g-xTB exchange stream get overlap dimensions")
      return
   end if

   allocate(amat(stream%nao, stream%nao), cmat(stream%nao, stream%nao), &
      & vmat(stream%nao, stream%nao), work(stream%nao, stream%nao))
   spin_factor = 0.5_wp
   if (stream%nspin > 1) spin_factor = 1.0_wp
   do spin = 1, stream%nspin
      amat = (0.0_wp, 0.0_wp)
      cmat = (0.0_wp, 0.0_wp)
      vmat = (0.0_wp, 0.0_wp)
      do icell = 1, stream%nk
         angle = 2.0_wp*pi*dot_product(stream%kfrac(:, ik), &
            & real(stream%reps(:, icell), wp))
         phase = exp(cmplx(0.0_wp, angle, wp))
         amat = amat + phase*stream%amat_r(:, :, icell, spin)
         cmat = cmat + phase*stream%cmat_r(:, :, icell, spin)
         vmat = vmat + phase*stream%vmat_r(:, :, icell, spin)
      end do
      work = amat
      do iao = 1, stream%nao
         work(:, iao) = work(:, iao) &
            & +0.25_wp*stream%gdiagP(iao, spin)*overlap(:, iao)
      end do
      work = work + 0.5_wp*matmul(overlap, vmat)
      fock(:, :, spin) = cmat + matmul(work, overlap)
      do iao = 1, stream%nao
         fock(:, iao, spin) = fock(:, iao, spin) &
            & +0.5_wp*stream%gdiagSP(iao, spin)*overlap(:, iao)
         fock(iao, iao, spin) = fock(iao, iao, spin) &
            & +0.25_wp*stream%gdiagSPS(iao, spin)
      end do
      do ii = 1, stream%nao
         fock(ii, ii, spin) = cmplx(-0.5_wp*spin_factor &
            & *real(fock(ii, ii, spin), wp), 0.0_wp, wp)
         do jj = 1, ii-1
            tmp = -0.25_wp*spin_factor*(fock(jj, ii, spin) &
               & +conjg(fock(ii, jj, spin)))
            fock(jj, ii, spin) = tmp
            fock(ii, jj, spin) = conjg(tmp)
         end do
      end do
   end do
end subroutine cp2k_exchange_stream_get


!> Close a successfully applied exchange stream and release all storage.
subroutine cp2k_exchange_stream_end(stream, error)
   type(cp2k_exchange_stream_type), intent(inout) :: stream
   type(error_type), allocatable, intent(out) :: error

   if (.not.stream%active) then
      call fatal_error(error, "g-xTB exchange stream is not active")
      return
   end if
   if (.not.stream%applied) then
      call clear_exchange_stream(stream)
      call fatal_error(error, &
         & "g-xTB exchange stream ended before successful apply")
      return
   end if
   if (stream%reverse_applied .and. .not.all(stream%reverse_pulled)) then
      call clear_exchange_stream(stream)
      call fatal_error(error, &
         & "g-xTB exchange stream ended with missing reverse pulls")
      return
   end if
   call clear_exchange_stream(stream)
end subroutine cp2k_exchange_stream_end


!> Differentiate BvK-supercell-equivalent exchange on a complete k mesh.
!>
!> The overlap response is spin-summed, Hermitian, and unweighted:
!> ``dE = sum_k w_k Re sum_ij conjg(overlap_grad_ij(k))*dS_ij(k)``.
!> Atomic and positive homogeneous-strain derivatives already contain the
!> complete primitive-cell Brillouin-zone contraction and must not be weighted
!> again by CP2K.  The strain response has no inverse-volume factor.
subroutine cp2k_exchange_kmesh_gradient(calc, mol, ecache, nmesh, kfrac, &
   & weights, density, overlap, overlap_grad, gradient, sigma, error)
   type(xtb_calculator), intent(in) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: ecache
   integer, intent(in) :: nmesh(3)
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   complex(wp), intent(out) :: overlap_grad(:, :, :)
   real(wp), intent(out) :: gradient(:, :), sigma(:, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: idim, nk
   real(wp), allocatable :: bocorr_grad_r(:, :, :), &
      & mulliken_grad_r(:, :, :)

   gradient = 0.0_wp
   sigma = 0.0_wp
   overlap_grad = (0.0_wp, 0.0_wp)
   nk = size(weights)

   if (.not.allocated(calc%exchange)) then
      call fatal_error(error, "g-xTB exchange container is not allocated")
      return
   end if
   if (.not.allocated(ecache%raw)) then
      call fatal_error(error, &
         & "k-point exchange cache is not prepared; call cp2k_prepare_exchange")
      return
   end if
   if (any(nmesh < 1) .or. product(nmesh) /= nk) then
      call fatal_error(error, "g-xTB exchange requires a complete regular k mesh")
      return
   end if
   do idim = 1, 3
      if (.not.mol%periodic(idim) .and. nmesh(idim) /= 1) then
         call fatal_error(error, &
            & "nonperiodic directions cannot be replicated by the g-xTB k mesh")
         return
      end if
   end do
   if (size(kfrac, 1) /= 3 .or. size(kfrac, 2) /= nk &
      & .or. size(overlap, 3) /= nk .or. size(overlap_grad, 3) /= nk &
      & .or. size(density, 4) /= nk) then
      call fatal_error(error, &
         & "inconsistent g-xTB whole-mesh gradient array dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(kfrac))) then
      call fatal_error(error, "g-xTB whole-mesh k points must be finite")
      return
   end if
   if (.not.all(ieee_is_finite(weights)) .or. &
      & abs(sum(weights)-1.0_wp) > 1.0e-12_wp &
      & .or. maxval(abs(weights-1.0_wp/real(nk, wp))) > 1.0e-12_wp) then
      call fatal_error(error, "g-xTB whole-mesh exchange requires uniform weights")
      return
   end if

   select type (cache => ecache%raw)
   type is (exchange_cache)
      select type (exchange => calc%exchange)
      type is (exchange_fock)
         if (size(overlap, 1) /= exchange%nao &
            & .or. size(overlap, 2) /= exchange%nao &
            & .or. size(overlap_grad, 1) /= exchange%nao &
            & .or. size(overlap_grad, 2) /= exchange%nao &
            & .or. size(density, 1) /= exchange%nao &
            & .or. size(density, 2) /= exchange%nao &
            & .or. (size(density, 3) /= 1 .and. size(density, 3) /= 2) &
            & .or. size(gradient, 1) /= 3 &
            & .or. size(gradient, 2) /= mol%nat &
            & .or. size(sigma, 1) /= 3 .or. size(sigma, 2) /= 3) then
            call fatal_error(error, &
               & "inconsistent g-xTB whole-mesh gradient AO or atom dimensions")
            return
         end if

         call prepare_bvk_plan(exchange, mol, cache, nmesh, kfrac, weights, error)
         if (allocated(error)) return

         allocate(mulliken_grad_r(exchange%nao, exchange%nao, nk), &
            & bocorr_grad_r(exchange%nao, exchange%nao, nk))
         call exchange%get_KGrad_kmesh(mol, cache, cache%bvk_kernel, &
            & kfrac, weights, density, overlap, overlap_grad, &
            & mulliken_grad_r, bocorr_grad_r)
         call exchange%get_bvk_Kmatrix_derivs(mol, cache%bvk_kernel, &
            & mulliken_grad_r, &
            & bocorr_grad_r, gradient, sigma)
      class default
         call fatal_error(error, &
            & "whole-mesh exchange gradient requires exchange_fock")
      end select
   class default
      call fatal_error(error, "unexpected g-xTB exchange cache type")
   end select

end subroutine cp2k_exchange_kmesh_gradient


!> Apply the exact complete-mesh reverse exchange evaluator once.
!>
!> This transition is legal only after the forward stream was applied.  The
!> complete primitive-cell atomic and positive-strain derivatives are returned
!> directly exactly once.  The unweighted overlap adjoints remain in the
!> stream and must subsequently be pulled blockwise with `reverse_get`.
subroutine cp2k_exchange_stream_reverse_apply(stream, calc, mol, ecache, &
   & gradient, sigma, error)
   type(cp2k_exchange_stream_type), intent(inout) :: stream
   type(xtb_calculator), intent(in) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: ecache
   real(wp), intent(out) :: gradient(:, :), sigma(:, :)
   type(error_type), allocatable, intent(out) :: error

   if (.not.stream%active) then
      call fatal_error(error, "g-xTB exchange stream is not active")
      return
   end if
   if (.not.stream%applied) then
      call fatal_error(error, &
         & "g-xTB exchange stream reverse requires forward apply")
      return
   end if
   if (stream%mode /= cp2k_exchange_stream_oracle) then
      call fatal_error(error, &
         & "g-xTB exchange stream reverse requires oracle mode")
      return
   end if
   if (stream%reverse_applied) then
      call fatal_error(error, "g-xTB exchange stream reverse was already applied")
      return
   end if
   if (size(gradient, 1) /= 3 .or. size(gradient, 2) /= mol%nat &
      & .or. size(sigma, 1) /= 3 .or. size(sigma, 2) /= 3) then
      call fatal_error(error, &
         & "inconsistent g-xTB exchange stream reverse dimensions")
      return
   end if
   if (.not.allocated(calc%exchange) .or. .not.allocated(ecache%raw)) then
      call fatal_error(error, &
         & "g-xTB exchange stream state was invalidated before reverse apply")
      return
   end if

   select type (cache => ecache%raw)
   type is (exchange_cache)
      select type (exchange => calc%exchange)
      type is (exchange_fock)
         if (exchange%nao /= stream%nao .or. exchange%nsh /= stream%nsh &
            & .or. .not.cache%bvk_matches(mol, stream%nmesh, &
            & stream%kfrac, stream%weights) &
            & .or. .not.exchange%bvk_model_matches(cache) &
            & .or. .not.exchange_stream_onsite_matches(stream, cache)) then
            call fatal_error(error, &
               & "g-xTB exchange stream state was invalidated before reverse apply")
            return
         end if
      class default
         call fatal_error(error, &
            & "exchange stream reverse requires the exchange_fock implementation")
         return
      end select
   class default
      call fatal_error(error, "unexpected g-xTB exchange cache type")
      return
   end select

   call cp2k_exchange_kmesh_gradient(calc, mol, ecache, stream%nmesh, &
      & stream%kfrac, stream%weights, stream%density, stream%overlap, &
      & stream%overlap_adjoint, gradient, sigma, error)
   if (allocated(error)) return
   stream%reverse_pulled = .false.
   stream%reverse_applied = .true.
end subroutine cp2k_exchange_stream_reverse_apply


!> Pull one unweighted overlap-adjoint block from the reverse stream.
subroutine cp2k_exchange_stream_reverse_get(stream, ik, overlap_adjoint, error)
   type(cp2k_exchange_stream_type), intent(inout) :: stream
   integer, intent(in) :: ik
   complex(wp), intent(out) :: overlap_adjoint(:, :)
   type(error_type), allocatable, intent(out) :: error

   if (.not.stream%active) then
      call fatal_error(error, "g-xTB exchange stream is not active")
      return
   end if
   if (.not.stream%reverse_applied) then
      call fatal_error(error, "g-xTB exchange stream reverse was not applied")
      return
   end if
   if (ik < 1 .or. ik > stream%nk) then
      call fatal_error(error, &
         & "g-xTB exchange stream reverse block index is out of range")
      return
   end if
   if (stream%reverse_pulled(ik)) then
      call fatal_error(error, "duplicate g-xTB exchange stream reverse pull")
      return
   end if
   if (size(overlap_adjoint, 1) /= stream%nao &
      & .or. size(overlap_adjoint, 2) /= stream%nao) then
      call fatal_error(error, &
         & "inconsistent g-xTB exchange stream overlap-adjoint dimensions")
      return
   end if

   overlap_adjoint = stream%overlap_adjoint(:, :, ik)
   stream%reverse_pulled(ik) = .true.
end subroutine cp2k_exchange_stream_reverse_get


!> Differentiate the g-xTB exchange functional for one complex k-point block.
!>
!> All returned quantities are unweighted.  For a reduced mesh the external
!> driver must expand `gradient` and `sigma` over the complete k star.  The
!> complex Hermitian `overlap_grad` transforms like a density matrix and is
!> intended to be folded to real-space overlap images by the same weighted
!> star transform as an ordinary k-point density.
subroutine cp2k_exchange_kpoint_gradient(calc, mol, ecache, density, overlap, &
   & overlap_grad, gradient, sigma, error)
   type(xtb_calculator), intent(in) :: calc
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: ecache
   complex(wp), intent(in) :: density(:, :, :), overlap(:, :)
   complex(wp), intent(out) :: overlap_grad(:, :)
   real(wp), intent(out) :: gradient(:, :), sigma(:, :)
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: mulliken_grad(:, :), bocorr_grad(:, :)

   gradient = 0.0_wp
   sigma = 0.0_wp
   overlap_grad = (0.0_wp, 0.0_wp)
   if (.not. allocated(calc%exchange)) then
      call fatal_error(error, "g-xTB exchange container is not allocated")
      return
   end if
   if (.not. allocated(ecache%raw)) then
      call fatal_error(error, &
         & "k-point exchange cache is not prepared; call cp2k_prepare_exchange")
      return
   end if

   select type (cache => ecache%raw)
   type is (exchange_cache)
      select type (exchange => calc%exchange)
      type is (exchange_fock)
         allocate(mulliken_grad(exchange%nao, exchange%nao), &
            & bocorr_grad(exchange%nao, exchange%nao))
         call exchange%get_KGrad_kpoint(mol, cache, density, overlap, &
            & mulliken_grad, bocorr_grad, overlap_grad)
         call exchange%get_mulliken_derivs_direct(mol, cache, mulliken_grad, &
            & gradient, sigma)
         call exchange%get_bocorr_derivs_direct(mol, cache, bocorr_grad, &
            & gradient, sigma)
      class default
         call fatal_error(error, &
            & "k-point exchange gradient requires the exchange_fock implementation")
      end select
   class default
      call fatal_error(error, "unexpected g-xTB exchange cache type")
   end select
end subroutine cp2k_exchange_kpoint_gradient

!> Complete Gamma-point H0/ACP/exchange and q-vSZP response for CP2K.
subroutine cp2k_gamma_gradient(calc, mol, bcache, acache, ecache, wfn_aux, &
   & accuracy, selfenergy, dsedcn, overlap, pot, wfn, wdensity, gradient, sigma)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(basis_cache), intent(in) :: bcache
   type(acp_cache), intent(inout) :: acache
   type(container_cache), intent(inout) :: ecache
   type(wavefunction_type), intent(in) :: wfn_aux
   real(wp), intent(in) :: accuracy
   real(wp), intent(in) :: selfenergy(:), dsedcn(:), overlap(:, :)
   type(potential_type), intent(in) :: pot
   type(wavefunction_type), intent(inout) :: wfn
   real(wp), intent(inout) :: wdensity(:, :, :)
   real(wp), intent(out) :: gradient(:, :), sigma(:, :)

   real(wp) :: cutoff
   real(wp), allocatable :: aniso_dip(:, :), cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: dEdcn(:), dEdad(:, :), dEdcnbas(:), dEdqbas(:)
   real(wp), allocatable :: h1(:, :, :), lattr(:, :)
   type(adjacency_list) :: list
   type(integral_type) :: ints
   type(potential_type) :: pot_work

   gradient = 0.0_wp
   sigma = 0.0_wp
   allocate(aniso_dip(3, mol%nat), dEdcn(mol%nat), dEdad(3, mol%nat), &
      & dEdcnbas(mol%nat), dEdqbas(mol%nat), source=0.0_wp)

   cutoff = cp2k_get_acp_cutoff(calc, accuracy)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)
   if (any(mol%periodic)) then
      call get_anisotropy(calc%h0, mol, lattr, list, aniso_dip)
   else
      call get_anisotropy(calc%h0, mol, aniso_dip)
   end if

   ! CP2K constructs the potential contribution to its KS matrix directly.
   ! The native gradient kernel additionally expects the atom/shell shifts
   ! spread to pot%vao.  Prepare that representation on a private copy so the
   ! persistent CP2K potential is neither changed nor accumulated twice.
   pot_work = pot
   call new_integral(ints, size(overlap, 1))
   ints%overlap = overlap
   allocate(h1(size(overlap, 1), size(overlap, 2), wfn%nspin))
   call add_pot_to_h1(calc%bas, ints, pot_work, h1)

   if (allocated(calc%ncoord)) then
      allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
      call calc%ncoord%get_cn(mol, cn, dcndr, dcndL)
   end if

   if (allocated(calc%exchange)) then
      call calc%exchange%get_gradient_w_overlap(mol, ecache, wfn, overlap, &
         & wdensity(:, :, 1), gradient, sigma)
   end if

   if (allocated(calc%acp)) then
      call get_acp_gradient(mol, lattr, list, calc%bas, bcache, calc%acp, acache, &
         & wfn, dEdcnbas, dEdqbas, gradient, sigma)
   end if

   call updown_to_magnet(wfn%density)
   call updown_to_magnet(wdensity)
   call get_hamiltonian_gradient(mol, lattr, list, calc%bas, bcache, calc%h0, &
      & selfenergy, dsedcn, aniso_dip, pot_work, wfn%density, wdensity, dEdcn, dEdad, &
      & dEdcnbas, dEdqbas, gradient, sigma)
   call magnet_to_updown(wfn%density)

   if (allocated(dcndr)) call gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
   if (allocated(dcndL)) call gemv(dcndL, dEdcn, sigma, beta=1.0_wp)

   call cp2k_apply_basis_gradient(calc, mol, wfn_aux, dEdcnbas, dEdqbas, &
      & gradient, sigma)

   if (any(mol%periodic)) then
      call get_anisotropy_gradient(calc%h0, mol, lattr, list, dEdad, gradient, sigma)
   else
      call get_anisotropy_gradient(calc%h0, mol, dEdad, gradient, sigma)
   end if
end subroutine cp2k_gamma_gradient


!> Image-resolved H0/Pulay and q-vSZP response for CP2K.
!>
!> `translations` must be the exact translation list used to construct the
!> supplied real-space matrices.  Their layout is
!> `(nao,nao,ntranslation,nspin)`.  CP2K supplies the spin channels in
!> up/down representation; this wrapper converts private copies to tblite's
!> charge--magnetization convention before the contraction.  The matrices
!> must satisfy P(R)^T=P(-R) and W(R)^T=W(-R).
!>
!> This routine deliberately contains only the one-electron H0, overlap
!> (Pulay), CN, anisotropy, and charge-dependent q-vSZP basis response.  ACP
!> and exchange image responses are independent contributions and must be
!> added by their respective CP2K bridges.
subroutine cp2k_image_h0_gradient(calc, mol, bcache, wfn_aux, accuracy, &
   & translations, selfenergy, dsedcn, overlap_gamma, pot, density_images, &
   & wdensity_images, gradient, sigma)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   type(basis_cache), intent(in) :: bcache
   type(wavefunction_type), intent(in) :: wfn_aux
   real(wp), intent(in) :: accuracy
   real(wp), intent(in) :: translations(:, :)
   real(wp), intent(in) :: selfenergy(:), dsedcn(:), overlap_gamma(:, :)
   type(potential_type), intent(in) :: pot
   real(wp), intent(in) :: density_images(:, :, :, :)
   real(wp), intent(in) :: wdensity_images(:, :, :, :)
   real(wp), intent(out) :: gradient(:, :), sigma(:, :)

   real(wp) :: cutoff
   real(wp), allocatable :: aniso_dip(:, :), cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: dEdcn(:), dEdad(:, :), dEdcnbas(:), dEdqbas(:), dEdq(:)
   real(wp), allocatable :: h1(:, :, :), density_work(:, :, :, :), &
      & wdensity_work(:, :, :, :)
   type(adjacency_list) :: list
   type(integral_type) :: ints
   type(potential_type) :: pot_work

   gradient = 0.0_wp
   sigma = 0.0_wp
   allocate(aniso_dip(3, mol%nat), dEdcn(mol%nat), dEdad(3, mol%nat), &
      & dEdcnbas(mol%nat), dEdqbas(mol%nat), dEdq(mol%nat), source=0.0_wp)

   cutoff = get_cutoff(calc%bas, accuracy)
   call new_adjacency_list(list, mol, translations, cutoff)
   if (any(mol%periodic)) then
      call get_anisotropy(calc%h0, mol, translations, list, aniso_dip)
   else
      call get_anisotropy(calc%h0, mol, aniso_dip)
   end if

   ! Match the Gamma bridge: spread the atom/shell shifts to AO potentials on
   ! a private copy, without changing CP2K's persistent potential object.
   pot_work = pot
   call new_integral(ints, size(overlap_gamma, 1))
   ints%overlap = overlap_gamma
   allocate(h1(size(overlap_gamma, 1), size(overlap_gamma, 2), &
      & size(density_images, 4)))
   call add_pot_to_h1(calc%bas, ints, pot_work, h1)

   if (allocated(calc%ncoord)) then
      allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
      call calc%ncoord%get_cn(mol, cn, dcndr, dcndL)
   end if

   density_work = density_images
   wdensity_work = wdensity_images
   call updown_to_magnet(density_work)
   call updown_to_magnet(wdensity_work)
   call get_hamiltonian_gradient(mol, translations, list, calc%bas, bcache, &
      & calc%h0, selfenergy, dsedcn, aniso_dip, pot_work, density_work, &
      & wdensity_work, dEdcn, dEdad, dEdcnbas, dEdqbas, gradient, sigma)

   if (allocated(dcndr)) call gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
   if (allocated(dcndL)) call gemv(dcndL, dEdcn, sigma, beta=1.0_wp)

   if (calc%bas%charge_dependent) then
      call calc%bas%get_basis_gradient(mol, dEdcnbas, dEdqbas, dEdq, gradient, sigma)
   end if

   if (allocated(wfn_aux%dqatdr)) then
      call gemv(wfn_aux%dqatdr(:, :, :, 1), dEdq, gradient, beta=1.0_wp)
   end if
   if (allocated(wfn_aux%dqatdL)) then
      call gemv(wfn_aux%dqatdL(:, :, :, 1), dEdq, sigma, beta=1.0_wp)
   end if

   if (any(mol%periodic)) then
      call get_anisotropy_gradient(calc%h0, mol, translations, list, dEdad, &
         & gradient, sigma)
   else
      call get_anisotropy_gradient(calc%h0, mol, dEdad, gradient, sigma)
   end if
end subroutine cp2k_image_h0_gradient

!> Extract one effective CGTO in the representation expected by CP2K.
subroutine cp2k_get_cgto(bas, cache, ish, izp, iat, ang, nprim, alpha, coeff)
   class(basis_type), intent(in) :: bas
   type(basis_cache), intent(in) :: cache
   integer, intent(in) :: ish, izp, iat
   integer, intent(out) :: ang, nprim
   real(wp), intent(out) :: alpha(:), coeff(:)

   associate(cgto => bas%cgto(ish, izp)%raw)
      ang = cgto%ang
      nprim = cgto%nprim
      alpha = 0.0_wp
      coeff = 0.0_wp
      alpha(:nprim) = cgto%alpha(:nprim)
      call cgto%get_coeffs(cache%cgto(ish, iat), coeff(:nprim))
   end associate
end subroutine cp2k_get_cgto

!> Evaluate multipole integrals with save_tblite's explicit basis cache.
pure subroutine cp2k_multipole_cgto(bas, cache, jsh, jzp, jat, ish, izp, iat, &
   & r2, vec, overlap, dpint, qpint)
   class(basis_type), intent(in) :: bas
   type(basis_cache), intent(in) :: cache
   integer, intent(in) :: jsh, jzp, jat, ish, izp, iat
   real(wp), intent(in) :: r2, vec(3)
   real(wp), intent(out) :: overlap(*), dpint(3, *), qpint(6, *)

   call multipole_cgto(bas%cgto(jsh, jzp)%raw, bas%cgto(ish, izp)%raw, &
      & cache%cgto(jsh, jat), cache%cgto(ish, iat), r2, vec, bas%intcut, &
      & overlap, dpint, qpint)
end subroutine cp2k_multipole_cgto

!> Evaluate multipole derivatives with save_tblite's explicit basis cache.
pure subroutine cp2k_multipole_grad_cgto(bas, cache, jsh, jzp, jat, ish, izp, iat, &
   & r2, vec, overlap, dpint, qpint, doverlap, ddpintj, dqpintj, ddpinti, dqpinti)
   class(basis_type), intent(in) :: bas
   type(basis_cache), intent(in) :: cache
   integer, intent(in) :: jsh, jzp, jat, ish, izp, iat
   real(wp), intent(in) :: r2, vec(3)
   real(wp), intent(out) :: overlap(*), dpint(3, *), qpint(6, *)
   real(wp), intent(out) :: doverlap(3, *)
   real(wp), intent(out) :: ddpintj(3, 3, *), dqpintj(3, 6, *)
   real(wp), intent(out) :: ddpinti(3, 3, *), dqpinti(3, 6, *)

   call multipole_grad_cgto(bas%cgto(jsh, jzp)%raw, bas%cgto(ish, izp)%raw, &
      & cache%cgto(jsh, jat), cache%cgto(ish, iat), r2, vec, bas%intcut, &
      & overlap, dpint, qpint, doverlap, ddpintj, dqpintj, ddpinti, dqpinti)
end subroutine cp2k_multipole_grad_cgto

end module tblite_cp2k_compat
