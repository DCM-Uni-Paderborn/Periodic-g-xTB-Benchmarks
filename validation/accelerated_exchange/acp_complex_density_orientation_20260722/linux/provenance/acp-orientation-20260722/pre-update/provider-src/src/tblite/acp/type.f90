! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/acp/type.f90
!> Provides the atomic correction potential for xTB.

!> Implementation of the atomic correction potential
module tblite_acp_type
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matinv_3x3
   use tblite_acp_cache, only : acp_cache
   use tblite_adjlist, only : adjacency_list
   use tblite_basis_cache, only : basis_cache
   use tblite_basis_type, only : basis_type, new_basis, cgto_container
   use tblite_blas, only : gemm
   use tblite_integral_overlap, only : overlap_cgto, overlap_grad_cgto, &
      & maxl, msao, smap
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_acp, get_acp, get_acp_images, get_acp_kmesh, &
      & get_acp_gradient, get_acp_image_gradient, &
      & get_acp_projector_image_gradient, get_acp_kmesh_gradient, &
      & acp_kmesh_gradient_stream_begin, &
      & acp_kmesh_gradient_stream_push, &
      & acp_kmesh_gradient_stream_end, &
      & acp_kmesh_gradient_stream_discard, &
      & acp_kmesh_gradient_stream_real_elements, &
      & acp_kmesh_gradient_stream_complex_elements

   type, public :: acp_type
      !> Auxiliary basis set for ACP projector
      class(basis_type), allocatable :: auxbas
      !> Energy levels of the auxiliary basis functions
      real(wp), allocatable :: levels(:, :)
   contains
      !> Update the ACP cache
      procedure :: update
      !> Set the integral cutoffs for a given accuracy
      procedure :: set_cutoff
   end type acp_type

   !> Image representation of the separable ACP operator.
   !>
   !> `projector_overlap(:, :, T)` stores
   !> \f$C_T=\langle p_T|\mu_0\rangle\f$.  Consequently its Bloch
   !> representation is \f$C(k)=\sum_T C_T\exp(-i kT)\f$.  The Hamiltonian
   !> images use the convention expected by CP2K,
   !> \f$H(k)=\sum_R H_R\exp(+i kR)\f$, and are obtained from the exact
   !> finite-image convolution
   !> \f$H_R=\sum_{T-V=R}C_T^T L C_V\f$.
   type, public :: acp_image_type
      real(wp), allocatable :: projector_translation(:, :)
      real(wp), allocatable :: projector_overlap(:, :, :)
      real(wp), allocatable :: scaled_projector_overlap(:, :, :)
      real(wp), allocatable :: level(:)
      real(wp), allocatable :: translation(:, :)
      integer, allocatable :: inverse(:)
      integer :: origin = 0
      logical :: inversion_closed = .false.
      real(wp), allocatable :: hamiltonian(:, :, :)
   end type acp_image_type

   !> Bounded accumulator for the weighted Bloch-space ACP derivative.
   !>
   !> The complete `(AO,AO,spin,k)` density is never retained.  The stream
   !> stores only the response of every compact projector image together with
   !> two auxiliary-by-valence work matrices and enforces a sequential,
   !> exactly-once k-point transaction.
   type, public :: acp_kmesh_gradient_stream_type
      private
      logical :: active = .false.
      integer :: nao = 0
      integer :: naux = 0
      integer :: nk_expected = 0
      integer :: nk_pushed = 0
      integer :: nspin = 0
      integer :: ntr = 0
      real(wp) :: weight_sum = 0.0_wp
      real(wp), allocatable :: projector_cell(:, :)
      real(wp), allocatable :: response(:, :, :, :)
      complex(wp), allocatable :: csk(:, :)
      complex(wp), allocatable :: work(:, :)
   end type acp_kmesh_gradient_stream_type

contains


!> Factory for a new atomic correction potential object
subroutine new_acp(self, mol, nproj, cgtp, levels, accuracy)
   !> ACP object
   type(acp_type), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of projectors per species
   integer, intent(in) :: nproj(:)
   !> Contracte gaussian type projector functions for each projector and species
   type(cgto_container), allocatable, intent(in) :: cgtp(:, :)
   !> Energy levels of the ACP auxiliary basis functions
   real(wp), intent(in) :: levels(:, :)
   !> Optional accuracy specification
   real(wp), intent(in), optional :: accuracy

   allocate(self%auxbas)
   call new_basis(self%auxbas, mol, nproj, cgtp, accuracy=accuracy)

   self%levels = levels

end subroutine new_acp


! Update ACP cache
subroutine update(self, mol, cache)
   !> Instance of the basis type
   class(acp_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(acp_cache), intent(inout) :: cache
   call self%auxbas%update(mol, cache%auxbas, .false.)

end subroutine update


!> Set the integral cutoffs for the auxiliary basis for a given accuracy
subroutine set_cutoff(self, accuracy)
   !> Instance of the basis set data
   class(acp_type), intent(inout) :: self
   !> Accuracy factor for electronic convergence thresholds
   real(wp), intent(in) :: accuracy

   call self%auxbas%set_cutoff(accuracy)

end subroutine set_cutoff


subroutine get_acp(mol, trans, list, bas, bcache, acp, acache, hamiltonian)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(inout) :: acache
   !> Effective Hamiltonian
   real(wp), intent(inout) :: hamiltonian(:, :)

   real(wp), allocatable :: pv_overlap(:, :)

   if (.not. allocated(acache%scaled_pv_overlap)) then
      allocate(acache%scaled_pv_overlap(acp%auxbas%nao, bas%nao), source=0.0_wp)
   else
      acache%scaled_pv_overlap(:, :) = 0.0_wp
   end if
   allocate(pv_overlap(acp%auxbas%nao, bas%nao), source=0.0_wp)

   ! Obtain the (scaled) projector-valence overlap matrix
   if (any(mol%periodic)) then
      call get_pv_overlap_3d(mol, trans, list, bas, bcache, acp, acache, &
         & pv_overlap)
   else
      call get_pv_overlap_0d(mol, bas, bcache, acp, acache, pv_overlap)
   end if

   ! Contract over the auxiliary projectors and add the ACP to the Hamiltonian
   call gemm(amat=pv_overlap, bmat=acache%scaled_pv_overlap, cmat=hamiltonian, &
      & transa='T', beta=1.0_wp)

end subroutine get_acp


!> Build image-resolved projector overlaps and their ACP convolution.
!>
!> Unlike a pair-local one-electron term, the ACP range is the difference set
!> of the projector-overlap translations.  Returning that complete difference
!> set is essential: assigning the Gamma-point ACP matrix to the origin image
!> would give an incorrect, k-independent operator.
subroutine get_acp_images(mol, trans, bas, bcache, acp, acache, images, &
   & build_hamiltonian, use_projector_cache)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points used for the projector-valence overlap
   real(wp), intent(in) :: trans(:, :)
   !> Valence basis set information
   class(basis_type), intent(in) :: bas
   !> Valence basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(inout) :: acache
   !> Image-resolved separable ACP operator
   type(acp_image_type), intent(out) :: images
   !> Form the quadratic projector-translation difference set and H(R).
   !> Direct Bloch-space callers only need the compact projector images.
   logical, intent(in), optional :: build_hamiltonian
   !> Reuse the geometry-static projector payload.  Disabling this provides
   !> an explicit direct-build oracle and invalidates any previous payload.
   logical, intent(in), optional :: use_projector_cache

   integer :: iat, jat, isp, jsp, itr, jtr, ktr, ninput, ntr, nhtr
   integer :: ish, jproj, is, js, ii, jj, iao, jao, nao, ij
   logical :: cache_enabled, cache_hit, with_hamiltonian
   real(wp) :: diff(3), level, r2, scale, tol2, vec(3)
   integer, allocatable :: aux_ang(:), aux_nprim(:), bas_ang(:), bas_nprim(:)
   real(wp), allocatable :: aux_alpha(:), aux_coeff(:), bas_alpha(:), &
      & bas_coeff(:)
   real(wp), allocatable :: htrans_work(:, :), projector_overlap(:, :, :), &
      & scaled_projector_overlap(:, :, :), stmp(:)

   call acp%update(mol, acache)

   call get_acp_projector_basis_signature(mol, bas, bcache, acp, acache, &
      & bas_ang, bas_nprim, bas_alpha, bas_coeff, aux_ang, aux_nprim, &
      & aux_alpha, aux_coeff)
   cache_enabled = .true.
   if (present(use_projector_cache)) cache_enabled = use_projector_cache
   cache_hit = cache_enabled .and. &
      & acp_projector_cache_matches(acache, mol, trans, bas, acp, &
      & bas_ang, bas_nprim, bas_alpha, bas_coeff, aux_ang, aux_nprim, &
      & aux_alpha, aux_coeff)

   if (cache_hit) then
      acache%projector_hits = acache%projector_hits + 1
      allocate(images%projector_translation, &
         & source=acache%projector_translation)
      allocate(images%projector_overlap, source=acache%projector_overlap)
      allocate(images%level, source=acache%projector_level)
      ntr = size(images%projector_translation, 2)
      allocate(images%scaled_projector_overlap(acp%auxbas%nao, bas%nao, ntr))
      do itr = 1, ntr
         do iao = 1, bas%nao
            images%scaled_projector_overlap(:, iao, itr) = &
               & images%level*images%projector_overlap(:, iao, itr)
         end do
      end do
      if (allocated(acache%scaled_pv_overlap)) then
         if (any(shape(acache%scaled_pv_overlap) /= &
            & [acp%auxbas%nao, bas%nao])) then
            deallocate(acache%scaled_pv_overlap)
         end if
      end if
      if (.not.allocated(acache%scaled_pv_overlap)) then
         allocate(acache%scaled_pv_overlap(acp%auxbas%nao, bas%nao))
      end if
      acache%scaled_pv_overlap = sum(images%scaled_projector_overlap, dim=3)
   else
      if (cache_enabled) then
         acache%projector_misses = acache%projector_misses + 1
      else
         acache%projector_bypasses = acache%projector_bypasses + 1
      end if
      ! Invalidate before any allocation or integral evaluation.  A failed or
      ! interrupted rebuild can therefore never expose the previous payload.
      acache%projector_valid = .false.

      ninput = size(trans, 2)
      allocate(projector_overlap(acp%auxbas%nao, bas%nao, ninput), &
         & scaled_projector_overlap(acp%auxbas%nao, bas%nao, ninput), &
         & images%level(acp%auxbas%nao), source=0.0_wp)

      if (allocated(acache%scaled_pv_overlap)) then
         if (any(shape(acache%scaled_pv_overlap) /= &
            & [acp%auxbas%nao, bas%nao])) then
            deallocate(acache%scaled_pv_overlap)
         end if
      end if
      if (.not.allocated(acache%scaled_pv_overlap)) then
         allocate(acache%scaled_pv_overlap(acp%auxbas%nao, bas%nao), &
            & source=0.0_wp)
      else
         acache%scaled_pv_overlap(:, :) = 0.0_wp
      end if

      ! Expand the shell levels once to the auxiliary-AO representation.
      do jat = 1, mol%nat
         jsp = mol%id(jat)
         js = acp%auxbas%ish_at(jat)
         do jproj = 1, acp%auxbas%nsh_id(jsp)
            jj = acp%auxbas%iao_sh(js+jproj)
            level = acp%levels(jproj, jsp)
            nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
            images%level(jj+1:jj+nao) = level
         end do
      end do

      allocate(stmp(msao(bas%maxl)*msao(acp%auxbas%maxl)))

      ! C_T = <p_T|mu_0>.  This is the same integral and finite image set used
      ! by get_pv_overlap_3d, only without summing over T prematurely.
      do iat = 1, mol%nat
         isp = mol%id(iat)
         is = bas%ish_at(iat)
         do jat = 1, mol%nat
            jsp = mol%id(jat)
            js = acp%auxbas%ish_at(jat)
            do itr = 1, ninput
               vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
               r2 = dot_product(vec, vec)
               do ish = 1, bas%nsh_id(isp)
                  ii = bas%iao_sh(is+ish)
                  do jproj = 1, acp%auxbas%nsh_id(jsp)
                     jj = acp%auxbas%iao_sh(js+jproj)
                     call overlap_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                        & bas%cgto(ish, isp)%raw, &
                        & acache%auxbas%cgto(jproj, jat), &
                        & bcache%cgto(ish, iat), r2, vec, &
                        & acp%auxbas%intcut, stmp)

                     level = acp%levels(jproj, jsp)
                     nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
                     do iao = 1, msao(bas%cgto(ish, isp)%raw%ang)
                        do jao = 1, nao
                           ij = jao + nao*(iao-1)
                           projector_overlap(jj+jao, ii+iao, itr) = &
                              & stmp(ij)
                           scaled_projector_overlap(jj+jao, ii+iao, itr) = &
                              & level*stmp(ij)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! The supplied CP2K image range can be wider than the ACP integral range.
      ! overlap_cgto applies the integral cutoff by returning exact zeros.  Drop
      ! precisely those empty projector translations before forming the O(ntr^2)
      ! difference set; no additional numerical/physical threshold is introduced.
      ntr = 0
      do itr = 1, ninput
         if (any(projector_overlap(:, :, itr) /= 0.0_wp)) ntr = ntr + 1
      end do
      allocate(images%projector_translation(3, ntr), &
         & images%projector_overlap(acp%auxbas%nao, bas%nao, ntr), &
         & images%scaled_projector_overlap(acp%auxbas%nao, bas%nao, ntr))
      ktr = 0
      do itr = 1, ninput
         if (.not.any(projector_overlap(:, :, itr) /= 0.0_wp)) cycle
         ktr = ktr + 1
         images%projector_translation(:, ktr) = trans(:, itr)
         images%projector_overlap(:, :, ktr) = &
            & projector_overlap(:, :, itr)
         images%scaled_projector_overlap(:, :, ktr) = &
            & scaled_projector_overlap(:, :, itr)
      end do

      acache%scaled_pv_overlap = &
         & sum(images%scaled_projector_overlap, dim=3)
      if (cache_enabled) then
         call store_acp_projector_cache(acache, mol, trans, bas, acp, &
            & bas_ang, bas_nprim, bas_alpha, bas_coeff, aux_ang, aux_nprim, &
            & aux_alpha, aux_coeff, images)
      end if
   end if

   with_hamiltonian = .true.
   if (present(build_hamiltonian)) with_hamiltonian = build_hamiltonian
   if (.not.with_hamiltonian) return

   ! First construct the unique translation-difference set.  The tolerance
   ! only absorbs floating-point noise from equivalent lattice differences.
   allocate(htrans_work(3, max(1, ntr*ntr)), source=0.0_wp)
   scale = max(1.0_wp, maxval(abs(trans)))
   tol2 = (1024.0_wp*epsilon(1.0_wp)*scale)**2
   nhtr = 0
   if (ntr == 0) then
      nhtr = 1
      htrans_work(:, nhtr) = 0.0_wp
   else
      do itr = 1, ntr
         do jtr = 1, ntr
            diff = images%projector_translation(:, itr) &
               & - images%projector_translation(:, jtr)
            ktr = find_translation(htrans_work(:, :nhtr), diff, tol2)
            if (ktr == 0) then
               nhtr = nhtr + 1
               htrans_work(:, nhtr) = diff
            end if
         end do
      end do
   end if

   allocate(images%translation(3, nhtr), &
      & images%hamiltonian(bas%nao, bas%nao, nhtr), source=0.0_wp)
   allocate(images%inverse(nhtr), source=0)
   images%translation = htrans_work(:, :nhtr)
   images%origin = find_translation(images%translation, &
      & [0.0_wp, 0.0_wp, 0.0_wp], tol2)
   images%inversion_closed = .true.
   do itr = 1, nhtr
      images%inverse(itr) = find_translation(images%translation, &
         & -images%translation(:, itr), tol2)
      images%inversion_closed = images%inversion_closed .and. &
         & images%inverse(itr) > 0
   end do

   ! H_R = sum_(T-V=R) C_T^T L C_V.  Each ordered translation pair is
   ! included, so H_R^T=H_-R and the Bloch sum is Hermitian at arbitrary k.
   do itr = 1, ntr
      do jtr = 1, ntr
         diff = images%projector_translation(:, itr) &
            & - images%projector_translation(:, jtr)
         ktr = find_translation(images%translation, diff, tol2)
         call gemm(amat=images%projector_overlap(:, :, itr), &
            & bmat=images%scaled_projector_overlap(:, :, jtr), &
            & cmat=images%hamiltonian(:, :, ktr), transa='T', beta=1.0_wp)
      end do
   end do

end subroutine get_acp_images


!> Build an exact signature of every primitive and effective contraction
!> coefficient entering the ACP projector integrals.
subroutine get_acp_projector_basis_signature(mol, bas, bcache, acp, acache, &
   & bas_ang, bas_nprim, bas_alpha, bas_coeff, aux_ang, aux_nprim, &
   & aux_alpha, aux_coeff)
   type(structure_type), intent(in) :: mol
   class(basis_type), intent(in) :: bas
   type(basis_cache), intent(in) :: bcache
   type(acp_type), intent(in) :: acp
   type(acp_cache), intent(in) :: acache
   integer, allocatable, intent(out) :: bas_ang(:), bas_nprim(:)
   integer, allocatable, intent(out) :: aux_ang(:), aux_nprim(:)
   real(wp), allocatable, intent(out) :: bas_alpha(:), bas_coeff(:)
   real(wp), allocatable, intent(out) :: aux_alpha(:), aux_coeff(:)

   integer :: iat, isp, ish, offset, shell, total

   total = 0
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ish = 1, bas%nsh_at(iat)
         total = total + bas%cgto(ish, isp)%raw%nprim
      end do
   end do
   allocate(bas_ang(bas%nsh), bas_nprim(bas%nsh), &
      & bas_alpha(total), bas_coeff(total))
   shell = 0
   offset = 0
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ish = 1, bas%nsh_at(iat)
         associate(cgto => bas%cgto(ish, isp)%raw)
            shell = shell + 1
            bas_ang(shell) = cgto%ang
            bas_nprim(shell) = cgto%nprim
            bas_alpha(offset+1:offset+cgto%nprim) = &
               & cgto%alpha(1:cgto%nprim)
            call cgto%get_coeffs(bcache%cgto(ish, iat), &
               & bas_coeff(offset+1:offset+cgto%nprim))
            offset = offset + cgto%nprim
         end associate
      end do
   end do

   total = 0
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ish = 1, acp%auxbas%nsh_at(iat)
         total = total + acp%auxbas%cgto(ish, isp)%raw%nprim
      end do
   end do
   allocate(aux_ang(acp%auxbas%nsh), aux_nprim(acp%auxbas%nsh), &
      & aux_alpha(total), aux_coeff(total))
   shell = 0
   offset = 0
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ish = 1, acp%auxbas%nsh_at(iat)
         associate(cgto => acp%auxbas%cgto(ish, isp)%raw)
            shell = shell + 1
            aux_ang(shell) = cgto%ang
            aux_nprim(shell) = cgto%nprim
            aux_alpha(offset+1:offset+cgto%nprim) = &
               & cgto%alpha(1:cgto%nprim)
            call cgto%get_coeffs(acache%auxbas%cgto(ish, iat), &
               & aux_coeff(offset+1:offset+cgto%nprim))
            offset = offset + cgto%nprim
         end associate
      end do
   end do
end subroutine get_acp_projector_basis_signature


!> Check whether a compact ACP projector payload was formed from exactly the
!> current geometry, cell, translation range, basis and model parameters.
logical function acp_projector_cache_matches(cache, mol, trans, bas, acp, &
   & bas_ang, bas_nprim, bas_alpha, bas_coeff, aux_ang, aux_nprim, &
   & aux_alpha, aux_coeff) result(matches)
   type(acp_cache), intent(in) :: cache
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: trans(:, :)
   class(basis_type), intent(in) :: bas
   type(acp_type), intent(in) :: acp
   integer, intent(in) :: bas_ang(:), bas_nprim(:), aux_ang(:), aux_nprim(:)
   real(wp), intent(in) :: bas_alpha(:), bas_coeff(:), aux_alpha(:), &
      & aux_coeff(:)

   matches = .false.
   if (.not.cache%projector_valid) return
   if (.not.allocated(cache%projector_id) .or. &
      & .not.allocated(cache%projector_xyz) .or. &
      & .not.allocated(cache%projector_input_translation) .or. &
      & .not.allocated(cache%projector_bas_ang) .or. &
      & .not.allocated(cache%projector_bas_nprim) .or. &
      & .not.allocated(cache%projector_bas_alpha) .or. &
      & .not.allocated(cache%projector_bas_coeff) .or. &
      & .not.allocated(cache%projector_aux_ang) .or. &
      & .not.allocated(cache%projector_aux_nprim) .or. &
      & .not.allocated(cache%projector_aux_alpha) .or. &
      & .not.allocated(cache%projector_aux_coeff) .or. &
      & .not.allocated(cache%projector_levels) .or. &
      & .not.allocated(cache%projector_translation) .or. &
      & .not.allocated(cache%projector_overlap) .or. &
      & .not.allocated(cache%projector_level)) return

   if (size(cache%projector_id) /= size(mol%id)) return
   if (any(shape(cache%projector_xyz) /= shape(mol%xyz))) return
   if (any(shape(cache%projector_input_translation) /= shape(trans))) return
   if (size(cache%projector_bas_ang) /= size(bas_ang)) return
   if (size(cache%projector_bas_nprim) /= size(bas_nprim)) return
   if (size(cache%projector_bas_alpha) /= size(bas_alpha)) return
   if (size(cache%projector_bas_coeff) /= size(bas_coeff)) return
   if (size(cache%projector_aux_ang) /= size(aux_ang)) return
   if (size(cache%projector_aux_nprim) /= size(aux_nprim)) return
   if (size(cache%projector_aux_alpha) /= size(aux_alpha)) return
   if (size(cache%projector_aux_coeff) /= size(aux_coeff)) return
   if (any(shape(cache%projector_levels) /= shape(acp%levels))) return
   if (any(shape(cache%projector_overlap) /= &
      & [acp%auxbas%nao, bas%nao, &
      & size(cache%projector_translation, 2)])) return
   if (size(cache%projector_level) /= acp%auxbas%nao) return

   if (any(cache%projector_periodic .neqv. mol%periodic)) return
   if (any(cache%projector_id /= mol%id)) return
   if (any(cache%projector_lattice /= mol%lattice)) return
   if (any(cache%projector_xyz /= mol%xyz)) return
   if (any(cache%projector_input_translation /= trans)) return
   if (any(cache%projector_intcut /= [bas%intcut, acp%auxbas%intcut])) return
   if (any(cache%projector_bas_ang /= bas_ang)) return
   if (any(cache%projector_bas_nprim /= bas_nprim)) return
   if (any(cache%projector_bas_alpha /= bas_alpha)) return
   if (any(cache%projector_bas_coeff /= bas_coeff)) return
   if (any(cache%projector_aux_ang /= aux_ang)) return
   if (any(cache%projector_aux_nprim /= aux_nprim)) return
   if (any(cache%projector_aux_alpha /= aux_alpha)) return
   if (any(cache%projector_aux_coeff /= aux_coeff)) return
   if (any(cache%projector_levels /= acp%levels)) return
   matches = .true.
end function acp_projector_cache_matches


!> Replace the compact geometry-static ACP projector cache after invalidation.
subroutine store_acp_projector_cache(cache, mol, trans, bas, acp, bas_ang, &
   & bas_nprim, bas_alpha, bas_coeff, aux_ang, aux_nprim, aux_alpha, &
   & aux_coeff, images)
   type(acp_cache), intent(inout) :: cache
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: trans(:, :)
   class(basis_type), intent(in) :: bas
   type(acp_type), intent(in) :: acp
   integer, intent(in) :: bas_ang(:), bas_nprim(:), aux_ang(:), aux_nprim(:)
   real(wp), intent(in) :: bas_alpha(:), bas_coeff(:), aux_alpha(:), &
      & aux_coeff(:)
   type(acp_image_type), intent(in) :: images

   cache%projector_valid = .false.
   cache%projector_periodic = mol%periodic
   cache%projector_lattice = mol%lattice
   cache%projector_intcut = [bas%intcut, acp%auxbas%intcut]
   cache%projector_id = mol%id
   cache%projector_xyz = mol%xyz
   cache%projector_input_translation = trans
   cache%projector_bas_ang = bas_ang
   cache%projector_bas_nprim = bas_nprim
   cache%projector_bas_alpha = bas_alpha
   cache%projector_bas_coeff = bas_coeff
   cache%projector_aux_ang = aux_ang
   cache%projector_aux_nprim = aux_nprim
   cache%projector_aux_alpha = aux_alpha
   cache%projector_aux_coeff = aux_coeff
   cache%projector_levels = acp%levels
   cache%projector_translation = images%projector_translation
   cache%projector_overlap = images%projector_overlap
   cache%projector_level = images%level
   cache%projector_valid = .true.
end subroutine store_acp_projector_cache


!> Evaluate the separable ACP directly on arbitrary Bloch points.
!>
!> With `C(k)=sum_T C_T exp(-2*pi*i*k.T)`, the returned unweighted
!> Hamiltonian is `H(k)=C(k)^H L C(k)`.  Direct evaluation avoids both the
!> quadratic translation-difference set and any real-space image truncation.
subroutine get_acp_kmesh(mol, images, kfrac, hamiltonian, error)
   type(structure_type), intent(in) :: mol
   type(acp_image_type), intent(in) :: images
   real(wp), intent(in) :: kfrac(:, :)
   complex(wp), intent(out) :: hamiltonian(:, :, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: ik, itr, nao, nk, ntr
   real(wp) :: angle
   real(wp), allocatable :: projector_cell(:, :)
   complex(wp) :: phase
   complex(wp), allocatable :: ck(:, :), csk(:, :)

   if (.not.allocated(images%projector_translation) .or. &
      & .not.allocated(images%projector_overlap) .or. &
      & .not.allocated(images%scaled_projector_overlap)) then
      call fatal_error(error, "ACP k-mesh evaluation requires projector images")
      return
   end if
   nao = size(images%projector_overlap, 2)
   ntr = size(images%projector_overlap, 3)
   nk = size(kfrac, 2)
   if (size(kfrac, 1) /= 3 .or. nk < 1 .or. &
      & any(shape(images%scaled_projector_overlap) /= &
      & shape(images%projector_overlap)) .or. &
      & any(shape(images%projector_translation) /= [3, ntr]) .or. &
      & any(shape(hamiltonian) /= [nao, nao, nk])) then
      call fatal_error(error, "inconsistent ACP k-mesh dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(kfrac))) then
      call fatal_error(error, "ACP k points must be finite")
      return
   end if

   call get_acp_projector_cells(mol, images, projector_cell, error)
   if (allocated(error)) return
   allocate(ck(size(images%projector_overlap, 1), nao), &
      & csk(size(images%projector_overlap, 1), nao))
   hamiltonian = (0.0_wp, 0.0_wp)
   do ik = 1, nk
      ck = (0.0_wp, 0.0_wp)
      csk = (0.0_wp, 0.0_wp)
      do itr = 1, ntr
         angle = -2.0_wp*pi*dot_product(kfrac(:, ik), projector_cell(:, itr))
         phase = exp(cmplx(0.0_wp, angle, wp))
         ck = ck + phase*images%projector_overlap(:, :, itr)
         csk = csk + phase*images%scaled_projector_overlap(:, :, itr)
      end do
      hamiltonian(:, :, ik) = matmul(conjg(transpose(ck)), csk)
   end do
end subroutine get_acp_kmesh


!> Differentiate a complete weighted Bloch contraction of the separable ACP.
!>
!> The density blocks are unweighted and Hermitian.  The real response of a
!> compact projector image is
!>
!> `D_T = 2 Re sum_k w_k exp(+2*pi*i*k.T) L C(k) P(k)`.
!>
!> Contracting `D_T` with `dC_T` is exactly dual to the direct k-space ACP
!> energy, including complex shifted meshes, without requiring density images
!> for the global `T-V` difference set.
subroutine get_acp_kmesh_gradient(mol, bas, bcache, acp, acache, images, &
   & kfrac, weights, density, dEdcnbas, dEdqbas, gradient, sigma, error)
   type(structure_type), intent(in) :: mol
   class(basis_type), intent(in) :: bas
   type(basis_cache), intent(in) :: bcache
   type(acp_type), intent(in) :: acp
   type(acp_cache), intent(in) :: acache
   type(acp_image_type), intent(in) :: images
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :)
   real(wp), intent(inout) :: dEdcnbas(:), dEdqbas(:)
   real(wp), intent(inout) :: gradient(:, :), sigma(:, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: ik, itr, nao, nk, nspin, ntr, spin
   real(wp) :: angle
   real(wp), allocatable :: projector_cell(:, :), response(:, :, :, :)
   complex(wp) :: phase
   complex(wp), allocatable :: csk(:, :), work(:, :)

   if (.not.allocated(images%projector_translation) .or. &
      & .not.allocated(images%projector_overlap) .or. &
      & .not.allocated(images%scaled_projector_overlap)) then
      call fatal_error(error, "ACP k-mesh gradient requires projector images")
      return
   end if
   nao = size(images%projector_overlap, 2)
   ntr = size(images%projector_overlap, 3)
   nk = size(weights)
   nspin = size(density, 3)
   if (size(kfrac, 1) /= 3 .or. size(kfrac, 2) /= nk .or. nk < 1 .or. &
      & size(density, 1) /= nao .or. size(density, 2) /= nao .or. &
      & size(density, 4) /= nk .or. nspin < 1 .or. &
      & any(shape(images%scaled_projector_overlap) /= &
      & shape(images%projector_overlap)) .or. &
      & any(shape(images%projector_translation) /= [3, ntr])) then
      call fatal_error(error, "inconsistent ACP k-mesh gradient dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(kfrac)) .or. &
      & .not.all(ieee_is_finite(weights)) .or. &
      & .not.all(ieee_is_finite(real(density, wp))) .or. &
      & .not.all(ieee_is_finite(aimag(density))) .or. &
      & any(weights < 0.0_wp) .or. &
      & abs(sum(weights)-1.0_wp) > 1.0e-12_wp) then
      call fatal_error(error, "invalid ACP k-mesh density or weights")
      return
   end if

   call get_acp_projector_cells(mol, images, projector_cell, error)
   if (allocated(error)) return
   allocate(response(size(images%projector_overlap, 1), nao, ntr, nspin), &
      & source=0.0_wp)
   allocate(csk(size(images%projector_overlap, 1), nao), &
      & work(size(images%projector_overlap, 1), nao))
   do ik = 1, nk
      csk = (0.0_wp, 0.0_wp)
      do itr = 1, ntr
         angle = -2.0_wp*pi*dot_product(kfrac(:, ik), projector_cell(:, itr))
         phase = exp(cmplx(0.0_wp, angle, wp))
         csk = csk + phase*images%scaled_projector_overlap(:, :, itr)
      end do
      do spin = 1, nspin
         work = matmul(csk, density(:, :, spin, ik))
         do itr = 1, ntr
            angle = -2.0_wp*pi*dot_product(kfrac(:, ik), projector_cell(:, itr))
            phase = exp(cmplx(0.0_wp, angle, wp))
            response(:, :, itr, spin) = response(:, :, itr, spin) + &
               & 2.0_wp*weights(ik)*real(conjg(phase)*work, wp)
         end do
      end do
   end do

   call get_pv_overlap_deriv_images(mol, images%projector_translation, bas, &
      & bcache, acp, acache, response, dEdcnbas, dEdqbas, gradient, sigma)
end subroutine get_acp_kmesh_gradient


!> Start a bounded exactly-once ACP k-mesh derivative transaction.
subroutine acp_kmesh_gradient_stream_begin(stream, mol, images, nspin, &
   & nk_expected, error)
   type(acp_kmesh_gradient_stream_type), intent(inout) :: stream
   type(structure_type), intent(in) :: mol
   type(acp_image_type), intent(in) :: images
   integer, intent(in) :: nspin, nk_expected
   type(error_type), allocatable, intent(out) :: error

   integer :: nao, naux, ntr

   if (stream%active) then
      call fatal_error(error, "ACP k-mesh gradient stream is already active")
      return
   end if
   call clear_acp_kmesh_gradient_stream(stream)
   if (.not.allocated(images%projector_translation) .or. &
      & .not.allocated(images%projector_overlap) .or. &
      & .not.allocated(images%scaled_projector_overlap)) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream requires projector images")
      return
   end if
   naux = size(images%projector_overlap, 1)
   nao = size(images%projector_overlap, 2)
   ntr = size(images%projector_overlap, 3)
   if (naux < 1 .or. nao < 1 .or. ntr < 1 .or. nspin < 1 .or. &
      & nk_expected < 1 .or. &
      & any(shape(images%scaled_projector_overlap) /= &
      & shape(images%projector_overlap)) .or. &
      & any(shape(images%projector_translation) /= [3, ntr])) then
      call fatal_error(error, &
         & "inconsistent ACP k-mesh gradient stream dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(images%projector_translation)) .or. &
      & .not.all(ieee_is_finite(images%projector_overlap)) .or. &
      & .not.all(ieee_is_finite(images%scaled_projector_overlap))) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream projector data must be finite")
      return
   end if

   call get_acp_projector_cells(mol, images, stream%projector_cell, error)
   if (allocated(error)) return
   allocate(stream%response(naux, nao, ntr, nspin), source=0.0_wp)
   allocate(stream%csk(naux, nao), stream%work(naux, nao))
   stream%nao = nao
   stream%naux = naux
   stream%ntr = ntr
   stream%nspin = nspin
   stream%nk_expected = nk_expected
   stream%nk_pushed = 0
   stream%weight_sum = 0.0_wp
   stream%active = .true.
end subroutine acp_kmesh_gradient_stream_begin


!> Add one ordered, unweighted Hermitian density block to an ACP derivative.
subroutine acp_kmesh_gradient_stream_push(stream, images, ik, kfrac, weight, &
   & density, error)
   type(acp_kmesh_gradient_stream_type), intent(inout) :: stream
   type(acp_image_type), intent(in) :: images
   integer, intent(in) :: ik
   real(wp), intent(in) :: kfrac(:), weight
   complex(wp), intent(in) :: density(:, :, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: itr, spin
   real(wp) :: angle
   complex(wp) :: phase

   if (.not.stream%active) then
      call fatal_error(error, "ACP k-mesh gradient stream is not active")
      return
   end if
   if (ik /= stream%nk_pushed+1 .or. ik > stream%nk_expected) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream k-point order changed")
      return
   end if
   if (.not.allocated(images%projector_translation) .or. &
      & .not.allocated(images%scaled_projector_overlap)) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream lost its projector images")
      return
   end if
   if (size(kfrac) /= 3 .or. &
      & any(shape(density) /= [stream%nao, stream%nao, stream%nspin]) .or. &
      & any(shape(images%scaled_projector_overlap) /= &
      & [stream%naux, stream%nao, stream%ntr]) .or. &
      & any(shape(images%projector_translation) /= [3, stream%ntr])) then
      call fatal_error(error, &
         & "inconsistent ACP k-mesh gradient stream block dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(kfrac)) .or. &
      & .not.ieee_is_finite(weight) .or. weight < 0.0_wp .or. &
      & .not.all(ieee_is_finite(real(density, wp))) .or. &
      & .not.all(ieee_is_finite(aimag(density)))) then
      call fatal_error(error, "invalid ACP k-mesh gradient stream block")
      return
   end if

   stream%csk = (0.0_wp, 0.0_wp)
   do itr = 1, stream%ntr
      angle = -2.0_wp*pi*dot_product(kfrac, &
         & stream%projector_cell(:, itr))
      phase = exp(cmplx(0.0_wp, angle, wp))
      stream%csk = stream%csk + &
         & phase*images%scaled_projector_overlap(:, :, itr)
   end do
   do spin = 1, stream%nspin
      stream%work = matmul(stream%csk, density(:, :, spin))
      do itr = 1, stream%ntr
         angle = -2.0_wp*pi*dot_product(kfrac, &
            & stream%projector_cell(:, itr))
         phase = exp(cmplx(0.0_wp, angle, wp))
         stream%response(:, :, itr, spin) = &
            & stream%response(:, :, itr, spin) + &
            & 2.0_wp*weight*real(conjg(phase)*stream%work, wp)
      end do
   end do
   stream%nk_pushed = ik
   stream%weight_sum = stream%weight_sum+weight
end subroutine acp_kmesh_gradient_stream_push


!> Finish one ACP derivative stream and contract projector derivatives once.
subroutine acp_kmesh_gradient_stream_end(stream, mol, bas, bcache, acp, &
   & acache, images, dEdcnbas, dEdqbas, gradient, sigma, error)
   type(acp_kmesh_gradient_stream_type), intent(inout) :: stream
   type(structure_type), intent(in) :: mol
   class(basis_type), intent(in) :: bas
   type(basis_cache), intent(in) :: bcache
   type(acp_type), intent(in) :: acp
   type(acp_cache), intent(in) :: acache
   type(acp_image_type), intent(in) :: images
   real(wp), intent(inout) :: dEdcnbas(:), dEdqbas(:)
   real(wp), intent(inout) :: gradient(:, :), sigma(:, :)
   type(error_type), allocatable, intent(out) :: error

   if (.not.stream%active) then
      call fatal_error(error, "ACP k-mesh gradient stream is not active")
      return
   end if
   if (stream%nk_pushed /= stream%nk_expected) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream ended with missing k points")
      return
   end if
   if (abs(stream%weight_sum-1.0_wp) > 1.0e-12_wp) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream weights are not normalized")
      return
   end if
   if (.not.allocated(images%projector_translation)) then
      call fatal_error(error, &
         & "ACP k-mesh gradient stream lost its projector translations")
      return
   end if
   if (bas%nao /= stream%nao .or. acp%auxbas%nao /= stream%naux .or. &
      & any(shape(images%projector_translation) /= [3, stream%ntr]) .or. &
      & size(dEdcnbas) /= mol%nat .or. size(dEdqbas) /= mol%nat .or. &
      & any(shape(gradient) /= [3, mol%nat]) .or. &
      & any(shape(sigma) /= [3, 3])) then
      call fatal_error(error, &
         & "inconsistent ACP k-mesh gradient stream output dimensions")
      return
   end if

   call get_pv_overlap_deriv_images(mol, images%projector_translation, bas, &
      & bcache, acp, acache, stream%response, dEdcnbas, dEdqbas, gradient, &
      & sigma)
   call clear_acp_kmesh_gradient_stream(stream)
end subroutine acp_kmesh_gradient_stream_end


!> Discard an ACP derivative stream without producing a result.
subroutine acp_kmesh_gradient_stream_discard(stream)
   type(acp_kmesh_gradient_stream_type), intent(inout) :: stream

   call clear_acp_kmesh_gradient_stream(stream)
end subroutine acp_kmesh_gradient_stream_discard


!> Count real scalar elements owned by an ACP derivative stream.
integer function acp_kmesh_gradient_stream_real_elements(stream) &
   & result(nelements)
   type(acp_kmesh_gradient_stream_type), intent(in) :: stream

   nelements = 0
   if (allocated(stream%projector_cell)) &
      & nelements = nelements+size(stream%projector_cell)
   if (allocated(stream%response)) nelements = nelements+size(stream%response)
end function acp_kmesh_gradient_stream_real_elements


!> Count complex scalar elements owned by an ACP derivative stream.
integer function acp_kmesh_gradient_stream_complex_elements(stream) &
   & result(nelements)
   type(acp_kmesh_gradient_stream_type), intent(in) :: stream

   nelements = 0
   if (allocated(stream%csk)) nelements = nelements+size(stream%csk)
   if (allocated(stream%work)) nelements = nelements+size(stream%work)
end function acp_kmesh_gradient_stream_complex_elements


!> Reset all storage and state of an ACP derivative stream.
subroutine clear_acp_kmesh_gradient_stream(stream)
   type(acp_kmesh_gradient_stream_type), intent(inout) :: stream

   if (allocated(stream%projector_cell)) deallocate(stream%projector_cell)
   if (allocated(stream%response)) deallocate(stream%response)
   if (allocated(stream%csk)) deallocate(stream%csk)
   if (allocated(stream%work)) deallocate(stream%work)
   stream%active = .false.
   stream%nao = 0
   stream%naux = 0
   stream%nk_expected = 0
   stream%nk_pushed = 0
   stream%nspin = 0
   stream%ntr = 0
   stream%weight_sum = 0.0_wp
end subroutine clear_acp_kmesh_gradient_stream


!> Differentiate compact ACP projector images against sparse density images.
!>
!> Missing translation differences are treated as exact zero-density blocks.
!> Unlike `get_acp_image_gradient`, neither the quadratic Hamiltonian
!> difference set nor a zero-padded density array is required.
subroutine get_acp_projector_image_gradient(mol, bas, bcache, acp, acache, &
   & images, density_translation, density, dEdcnbas, dEdqbas, gradient, &
   & sigma, error)
   type(structure_type), intent(in) :: mol
   class(basis_type), intent(in) :: bas
   type(basis_cache), intent(in) :: bcache
   type(acp_type), intent(in) :: acp
   type(acp_cache), intent(in) :: acache
   type(acp_image_type), intent(in) :: images
   real(wp), intent(in) :: density_translation(:, :)
   real(wp), intent(in) :: density(:, :, :, :)
   real(wp), intent(inout) :: dEdcnbas(:), dEdqbas(:)
   real(wp), intent(inout) :: gradient(:, :), sigma(:, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: dtr, itr, jtr, ktr, nspin, ntr, spin
   integer, allocatable :: density_map(:)
   real(wp) :: diff(3), response_residual, scale, tol2
   real(wp), allocatable :: density_cell(:, :), projector_cell(:, :), &
      & response(:, :, :, :), response_oracle(:, :, :, :)
   character(len=256) :: message

   if (.not.allocated(images%projector_translation) .or. &
      & .not.allocated(images%projector_overlap) .or. &
      & .not.allocated(images%scaled_projector_overlap)) then
      call fatal_error(error, &
         & "ACP sparse image gradient requires compact projector images")
      return
   end if
   ntr = size(images%projector_translation, 2)
   nspin = size(density, 4)
   if (size(images%projector_overlap, 1) /= acp%auxbas%nao .or. &
      & size(images%projector_overlap, 2) /= bas%nao .or. &
      & size(images%projector_overlap, 3) /= ntr .or. &
      & any(shape(images%scaled_projector_overlap) /= &
      & shape(images%projector_overlap)) .or. &
      & any(shape(images%projector_translation) /= [3, ntr]) .or. &
      & size(density_translation, 1) /= 3 .or. &
      & size(density, 1) /= bas%nao .or. size(density, 2) /= bas%nao .or. &
      & size(density, 3) /= size(density_translation, 2) .or. &
      & size(density, 3) < 1 .or. nspin < 1 .or. &
      & size(dEdcnbas) /= mol%nat .or. size(dEdqbas) /= mol%nat .or. &
      & any(shape(gradient) /= [3, mol%nat]) .or. &
      & any(shape(sigma) /= [3, 3])) then
      call fatal_error(error, &
         & "inconsistent ACP sparse image gradient dimensions")
      return
   end if
   if (.not.all(ieee_is_finite(density_translation)) .or. &
      & .not.all(ieee_is_finite(density))) then
      call fatal_error(error, &
         & "ACP sparse image gradient inputs must be finite")
      return
   end if

   call get_acp_translation_cells(mol, images%projector_translation, &
      & projector_cell, error)
   if (allocated(error)) return
   call get_acp_translation_cells(mol, density_translation, density_cell, error)
   if (allocated(error)) return
   allocate(response(acp%auxbas%nao, bas%nao, ntr, nspin), source=0.0_wp)
   do itr = 1, ntr
      do jtr = 1, ntr
         diff = projector_cell(:, itr)-projector_cell(:, jtr)
         dtr = find_translation(density_cell, diff, 0.0_wp)
         if (dtr > 0) then
            do spin = 1, nspin
               call gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                  & bmat=density(:, :, dtr, spin), &
                  & cmat=response(:, :, itr, spin), transb='T', beta=1.0_wp)
            end do
         end if

         diff = projector_cell(:, jtr)-projector_cell(:, itr)
         dtr = find_translation(density_cell, diff, 0.0_wp)
         if (dtr > 0) then
            do spin = 1, nspin
               call gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                  & bmat=density(:, :, dtr, spin), &
                  & cmat=response(:, :, itr, spin), beta=1.0_wp)
            end do
         end if
      end do
   end do

   if (allocated(images%translation) .and. allocated(images%hamiltonian)) then
      scale = max(1.0_wp, maxval(abs(images%translation)), &
         & maxval(abs(density_translation)))
      tol2 = (2048.0_wp*epsilon(1.0_wp)*scale)**2
      allocate(density_map(size(images%translation, 2)), source=0)
      do ktr = 1, size(images%translation, 2)
         density_map(ktr) = find_translation(density_translation, &
            & images%translation(:, ktr), tol2)
         if (density_map(ktr) == 0) then
            call fatal_error(error, &
               & "ACP sparse qualification lacks a dense-oracle translation")
            return
         end if
      end do
      allocate(response_oracle, mold=response)
      response_oracle = 0.0_wp
      do itr = 1, ntr
         do jtr = 1, ntr
            diff = images%projector_translation(:, itr) &
               & - images%projector_translation(:, jtr)
            ktr = find_translation(images%translation, diff, tol2)
            dtr = density_map(ktr)
            do spin = 1, nspin
               call gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                  & bmat=density(:, :, dtr, spin), &
                  & cmat=response_oracle(:, :, itr, spin), transb='T', &
                  & beta=1.0_wp)
            end do
            diff = images%projector_translation(:, jtr) &
               & - images%projector_translation(:, itr)
            ktr = find_translation(images%translation, diff, tol2)
            dtr = density_map(ktr)
            do spin = 1, nspin
               call gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
                  & bmat=density(:, :, dtr, spin), &
                  & cmat=response_oracle(:, :, itr, spin), beta=1.0_wp)
            end do
         end do
      end do
      if (.not.all(ieee_is_finite(response)) .or. &
         & .not.all(ieee_is_finite(response_oracle))) then
         call fatal_error(error, &
            & "ACP sparse qualification produced a non-finite response")
         return
      end if
      scale = max(1.0_wp, maxval(abs(response_oracle)))
      response_residual = maxval(abs(response-response_oracle))/scale
      write(*, '(a,es12.4)') &
         & " ACP_SPARSE_RESPONSE_QUALIFY residual=", response_residual
      if (response_residual > 1.0e-12_wp) then
         write(message, '(a,es24.16,a,es24.16)') &
            & "ACP sparse response exceeded its dense-oracle gate: relative=", &
            & response_residual, ", scale=", scale
         call fatal_error(error, trim(message))
         return
      end if
   end if

   call get_pv_overlap_deriv_images(mol, images%projector_translation, bas, &
      & bcache, acp, acache, response, dEdcnbas, dEdqbas, gradient, sigma)
end subroutine get_acp_projector_image_gradient


!> Convert Cartesian projector translations to integer primitive-cell labels.
subroutine get_acp_projector_cells(mol, images, cells, error)
   type(structure_type), intent(in) :: mol
   type(acp_image_type), intent(in) :: images
   real(wp), allocatable, intent(out) :: cells(:, :)
   type(error_type), allocatable, intent(out) :: error

   call get_acp_translation_cells(mol, images%projector_translation, cells, &
      & error)
end subroutine get_acp_projector_cells


!> Convert Cartesian translations to exact integer primitive-cell labels.
subroutine get_acp_translation_cells(mol, translation, cells, error)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: translation(:, :)
   real(wp), allocatable, intent(out) :: cells(:, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: itr
   real(wp) :: scale
   real(wp) :: invlat(3, 3)

   allocate(cells(3, size(translation, 2)), source=0.0_wp)
   if (.not.any(mol%periodic)) return
   invlat = matinv_3x3(mol%lattice)
   do itr = 1, size(cells, 2)
      cells(:, itr) = matmul(invlat, translation(:, itr))
   end do
   scale = max(1.0_wp, maxval(abs(cells)))
   if (maxval(abs(cells-real(nint(cells), wp))) > &
      & 4096.0_wp*epsilon(1.0_wp)*scale) then
      call fatal_error(error, "ACP translation is not a lattice vector")
      return
   end if
   cells = real(nint(cells), wp)
end subroutine get_acp_translation_cells


!> Locate a Cartesian lattice translation within a numerical tolerance.
pure function find_translation(trans, target, tol2) result(index)
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: target(:)
   real(wp), intent(in) :: tol2
   integer :: index
   integer :: itr

   index = 0
   do itr = 1, size(trans, 2)
      if (sum((trans(:, itr) - target)**2) <= tol2) then
         index = itr
         exit
      end if
   end do
end function find_translation


subroutine get_pv_overlap_0d(mol, bas, bcache, acp, acache, pv_overlap)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(inout) :: acache
   !> Projector-valence overlap matrix
   real(wp), intent(out) :: pv_overlap(:, :)

   integer :: iat, jat, isp, jsp, k, img, inl
   integer :: ish, jproj, is, js, ii, jj, iao, jao, nao, ij, iaosh, jaoproj
   real(wp) :: r2, vec(3), level
   real(wp), allocatable :: stmp(:)

   allocate(stmp(msao(bas%maxl)*msao(acp%auxbas%maxl)))
   
   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, bcache, acp, acache, pv_overlap) &
   !$omp private(iat, jat, isp, jsp, inl, img, is, js, ish, jproj, ii, jj) &
   !$omp private(iao, jao, iaosh, jaoproj, nao, ij, r2, vec, stmp, level)
   do iat = 1, mol%nat
      isp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jsp = mol%id(jat)
         js = acp%auxbas%ish_at(jat)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         ! Loop over valence shells
         do ish = 1, bas%nsh_id(isp)
            ii = bas%iao_sh(is+ish)
            ! Loop over projectors
            do jproj = 1, acp%auxbas%nsh_id(jsp)
               jj = acp%auxbas%iao_sh(js+jproj)

               call overlap_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                  & bas%cgto(ish, isp)%raw, acache%auxbas%cgto(jproj, jat), &
                  & bcache%cgto(ish, iat), r2, vec, acp%auxbas%intcut, stmp)

               level = acp%levels(jproj, jsp)
               
               nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
               do iao = 1, msao(bas%cgto(ish, isp)%raw%ang)
                  do jao = 1, nao
                     ij = jao + nao*(iao-1)

                     ! Store scaled overlap for later gradient calculation
                     acache%scaled_pv_overlap(jj+jao, ii+iao) = level * stmp(ij)

                     pv_overlap(jj+jao, ii+iao) = stmp(ij)
                  end do
               end do

            end do
         end do

      end do
   end do

end subroutine get_pv_overlap_0d


subroutine get_pv_overlap_3d(mol, trans, list, bas, bcache, acp, acache, &
   & pv_overlap)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(inout) :: acache
   !> Projector-valence overlap matrix
   real(wp), intent(out) :: pv_overlap(:, :)

   integer :: iat, jat, isp, jsp, itr, k, img, inl
   integer :: ish, jproj, is, js, ii, jj, iao, jao, nao, ij, iaosh, jaoproj
   real(wp) :: r2, vec(3), level
   real(wp), allocatable :: stmp(:)

   allocate(stmp(msao(bas%maxl)*msao(acp%auxbas%maxl)))
   
   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, bcache, acp, acache, pv_overlap, trans) &
   !$omp private(iat, jat, isp, jsp, itr, inl, img, is, js, ish, jproj, ii, jj) &
   !$omp private(iao, jao, iaosh, jaoproj, nao, ij, r2, vec, stmp, level)
   do iat = 1, mol%nat
      isp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jsp = mol%id(jat)
         js = acp%auxbas%ish_at(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            ! Loop over valence shells
            do ish = 1, bas%nsh_id(isp)
               ii = bas%iao_sh(is+ish)
               ! Loop over projectors
               do jproj = 1, acp%auxbas%nsh_id(jsp)
                  jj = acp%auxbas%iao_sh(js+jproj)

                  call overlap_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                     & bas%cgto(ish, isp)%raw, acache%auxbas%cgto(jproj, jat), &
                     & bcache%cgto(ish, iat), r2, vec, acp%auxbas%intcut, stmp)

                  level = acp%levels(jproj, jsp)
                  
                  nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
                  do iao = 1, msao(bas%cgto(ish, isp)%raw%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)

                        ! Store scaled overlap for later gradient calculation
                        !$omp atomic
                        acache%scaled_pv_overlap(jj+jao, ii+iao) = acache%scaled_pv_overlap(jj+jao, ii+iao) + level * stmp(ij)
                        !$omp atomic
                        pv_overlap(jj+jao, ii+iao) = pv_overlap(jj+jao, ii+iao) + stmp(ij)
                     end do
                  end do

               end do
            end do
         end do
      end do
   end do


end subroutine get_pv_overlap_3d


subroutine get_acp_gradient(mol, trans, list, bas, bcache, acp, acache, wfn, &
   & dEdcnbas, dEdqbas, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(in) :: acache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Derivative of the electronic energy w.r.t. basis set coordination numbers
   real(wp), intent(inout) :: dEdcnbas(:)
   !> Derivative of the electronic energy w.r.t. basis set atomic charges
   real(wp), intent(inout) :: dEdqbas(:)
   !> Derivative of the electronic energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the electronic energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)

   if (any(mol%periodic)) then
      call get_pv_overlap_deriv_3d(mol, trans, list, bas, bcache, acp, acache, &
         & wfn%density, dEdcnbas, dEdqbas, gradient, sigma)
   else
      call get_pv_overlap_deriv_0d(mol, bas, bcache, acp, acache, wfn%density, &
         & dEdcnbas, dEdqbas, gradient, sigma)
   end if

end subroutine get_acp_gradient


!> Contract the ACP response with image-resolved real-space densities.
!>
!> The density images use the dual convention of `acp_image_type`: the ACP
!> energy is the element-wise sum of `P(R)*H(R)` over images, spins and AOs.
!>
!> With `H_R = sum_(T-V=R) C_T^T L C_V`, the coefficient of an overlap
!> derivative `dC_U` is
!>
!> `D_U = sum_V L C_V P_(U-V)^T + sum_T L C_T P_(T-U)`.
!>
!> Both terms are evaluated explicitly.  For inversion-compatible density
!> images, `P_(-R)=P_R^T`, they are identical.  If every density image is the
!> same symmetric Gamma-point density, this expression reduces exactly to the
!> `get_acp_gradient` intermediate `2 L C P`.
subroutine get_acp_image_gradient(mol, bas, bcache, acp, acache, images, &
   & density_translation, density, dEdcnbas, dEdqbas, gradient, sigma, error)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Valence basis set information
   class(basis_type), intent(in) :: bas
   !> Valence basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(in) :: acache
   !> Projector and Hamiltonian image data at the current geometry
   type(acp_image_type), intent(in) :: images
   !> Cartesian translations of the density images
   real(wp), intent(in) :: density_translation(:, :)
   !> Density images `(AO, AO, image, spin)` dual to `images%hamiltonian`
   real(wp), intent(in) :: density(:, :, :, :)
   !> Derivative w.r.t. basis-set coordination numbers
   real(wp), intent(inout) :: dEdcnbas(:)
   !> Derivative w.r.t. basis-set atomic charges
   real(wp), intent(inout) :: dEdqbas(:)
   !> Cartesian ACP gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> ACP strain derivative
   real(wp), intent(inout) :: sigma(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: itr, jtr, ktr, dtr, ntr, nspin, nearest
   integer, allocatable :: density_map(:)
   real(wp) :: diff(3), nearest2, scale, tol2
   real(wp), allocatable :: response(:, :, :, :)
   character(len=512) :: message

   if (.not. allocated(images%projector_translation) .or. &
      & .not. allocated(images%projector_overlap) .or. &
      & .not. allocated(images%scaled_projector_overlap) .or. &
      & .not. allocated(images%translation) .or. &
      & .not. allocated(images%hamiltonian)) then
      call fatal_error(error, "ACP image gradient requires initialized image data")
      return
   end if

   ntr = size(images%projector_translation, 2)
   nspin = size(density, 4)
   if (size(images%projector_overlap, 1) /= acp%auxbas%nao .or. &
      & size(images%projector_overlap, 2) /= bas%nao .or. &
      & size(images%projector_overlap, 3) /= ntr .or. &
      & any(shape(images%scaled_projector_overlap) /= &
      & shape(images%projector_overlap))) then
      call fatal_error(error, "inconsistent ACP projector image dimensions")
      return
   end if
   if (size(density, 1) /= bas%nao .or. size(density, 2) /= bas%nao .or. &
      & size(density, 3) /= size(density_translation, 2) .or. &
      & size(density, 3) < 1 .or. nspin < 1) then
      call fatal_error(error, "inconsistent ACP density image dimensions")
      return
   end if
   if (size(density_translation, 1) /= 3 .or. &
      & size(images%translation, 1) /= 3 .or. &
      & size(images%projector_translation, 1) /= 3) then
      call fatal_error(error, "ACP image translations must be three-dimensional")
      return
   end if
   if (size(images%hamiltonian, 1) /= bas%nao .or. &
      & size(images%hamiltonian, 2) /= bas%nao .or. &
      & size(images%hamiltonian, 3) /= size(images%translation, 2)) then
      call fatal_error(error, "inconsistent ACP Hamiltonian image dimensions")
      return
   end if
   if (size(dEdcnbas) /= mol%nat .or. size(dEdqbas) /= mol%nat .or. &
      & any(shape(gradient) /= [3, mol%nat]) .or. &
      & any(shape(sigma) /= [3, 3])) then
      call fatal_error(error, "inconsistent ACP image gradient output dimensions")
      return
   end if

   scale = max(1.0_wp, maxval(abs(images%translation)), &
      & maxval(abs(density_translation)))
   tol2 = (2048.0_wp*epsilon(1.0_wp)*scale)**2
   allocate(density_map(size(images%translation, 2)), source=0)
   do ktr = 1, size(images%translation, 2)
      density_map(ktr) = find_translation(density_translation, &
         & images%translation(:, ktr), tol2)
      if (density_map(ktr) == 0) then
         nearest = 0
         nearest2 = huge(1.0_wp)
         do dtr = 1, size(density_translation, 2)
            diff = density_translation(:, dtr) - images%translation(:, ktr)
            if (sum(diff**2) < nearest2) then
               nearest = dtr
               nearest2 = sum(diff**2)
            end if
         end do
         write(message, '(a,i0,a,3(es16.8,1x),a,i0,a,i0,a,es16.8)') &
            & "density image set does not cover ACP translation ", ktr, &
            & " at [", images%translation(:, ktr), "]; density images=", &
            & size(density_translation, 2), ", nearest=", nearest, &
            & ", distance=", sqrt(nearest2)
         call fatal_error(error, trim(message))
         return
      end if
   end do

   allocate(response(acp%auxbas%nao, bas%nao, ntr, nspin), source=0.0_wp)
   do itr = 1, ntr
      do jtr = 1, ntr
         ! dC_U from the left projector factor:
         !       L C_V P_(U-V)^T
         diff = images%projector_translation(:, itr) &
            & - images%projector_translation(:, jtr)
         ktr = find_translation(images%translation, diff, tol2)
         if (ktr == 0) then
            call fatal_error(error, "incomplete ACP Hamiltonian difference set")
            return
         end if
         dtr = density_map(ktr)
         do ktr = 1, nspin
            call gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
               & bmat=density(:, :, dtr, ktr), &
               & cmat=response(:, :, itr, ktr), transb='T', beta=1.0_wp)
         end do

         ! dC_U from the right projector factor:
         !       L C_T P_(T-U)
         diff = images%projector_translation(:, jtr) &
            & - images%projector_translation(:, itr)
         ktr = find_translation(images%translation, diff, tol2)
         if (ktr == 0) then
            call fatal_error(error, "incomplete ACP Hamiltonian difference set")
            return
         end if
         dtr = density_map(ktr)
         do ktr = 1, nspin
            call gemm(amat=images%scaled_projector_overlap(:, :, jtr), &
               & bmat=density(:, :, dtr, ktr), &
               & cmat=response(:, :, itr, ktr), beta=1.0_wp)
         end do
      end do
   end do

   call get_pv_overlap_deriv_images(mol, images%projector_translation, bas, &
      & bcache, acp, acache, response, dEdcnbas, dEdqbas, gradient, sigma)
end subroutine get_acp_image_gradient


!> Contract derivatives of every compact projector image with its response.
subroutine get_pv_overlap_deriv_images(mol, trans, bas, bcache, acp, acache, &
   & response, dEdcnbas, dEdqbas, gradient, sigma)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: trans(:, :)
   class(basis_type), intent(in) :: bas
   type(basis_cache), intent(in) :: bcache
   type(acp_type), intent(in) :: acp
   type(acp_cache), intent(in) :: acache
   !> Coefficient of every projector-overlap derivative `(aux, AO, T, spin)`
   real(wp), intent(in) :: response(:, :, :, :)
   real(wp), intent(inout) :: dEdcnbas(:)
   real(wp), intent(inout) :: dEdqbas(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, isp, jsp, spin, nspin, itr
   integer :: ish, jproj, is, js, ii, jj, iao, jao, nao, ij
   real(wp) :: r2, vec(3), dG(3), dcnbasi, dqbasi, tmp
   real(wp), allocatable :: stmp(:), dstmp(:, :)
   real(wp), allocatable :: dstmpdqeffi(:), dstmpdqeffj(:)
   logical :: compute_qeff_grad
   real(wp), allocatable :: dEdcnbas_local(:), dEdqbas_local(:), &
      & gradient_local(:, :), sigma_local(:, :)

   compute_qeff_grad = bas%charge_dependent .or. acp%auxbas%charge_dependent
   nspin = size(response, 4)
   allocate(stmp(msao(bas%maxl)*msao(acp%auxbas%maxl)), &
      & dstmp(3, msao(bas%maxl)*msao(acp%auxbas%maxl)))
   if (compute_qeff_grad) then
      allocate(dstmpdqeffi(msao(bas%maxl)*msao(acp%auxbas%maxl)), &
         & dstmpdqeffj(msao(bas%maxl)*msao(acp%auxbas%maxl)))
   end if

   do spin = 1, nspin
      !$omp parallel default(none) shared(dEdcnbas, dEdqbas, gradient, sigma, trans) &
      !$omp shared(mol, bas, bcache, acp, acache, response, compute_qeff_grad, spin) &
      !$omp private(iat, jat, itr, isp, jsp, is, js, ish, jproj, ii, jj, iao, jao, nao, ij) &
      !$omp private(r2, vec, stmp, dstmp, dstmpdqeffi, dstmpdqeffj, dcnbasi, dqbasi) &
      !$omp private(tmp, dG, dEdcnbas_local, dEdqbas_local, gradient_local, sigma_local)
      if (compute_qeff_grad) then
         allocate(dEdcnbas_local(size(dEdcnbas)), source=0.0_wp)
         allocate(dEdqbas_local(size(dEdqbas)), source=0.0_wp)
      end if
      allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source=0.0_wp)
      allocate(sigma_local(size(sigma, 1), size(sigma, 2)), source=0.0_wp)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         isp = mol%id(iat)
         is = bas%ish_at(iat)
         dcnbasi = 0.0_wp
         dqbasi = 0.0_wp
         do jat = 1, mol%nat
            jsp = mol%id(jat)
            js = acp%auxbas%ish_at(jat)
            do itr = 1, size(trans, 2)
               vec = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
               r2 = dot_product(vec, vec)
               dG = 0.0_wp
               do ish = 1, bas%nsh_id(isp)
                  ii = bas%iao_sh(is+ish)
                  do jproj = 1, acp%auxbas%nsh_id(jsp)
                     jj = acp%auxbas%iao_sh(js+jproj)
                     if (compute_qeff_grad) then
                        call overlap_grad_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                           & bas%cgto(ish, isp)%raw, &
                           & acache%auxbas%cgto(jproj, jat), &
                           & bcache%cgto(ish, iat), r2, vec, &
                           & acp%auxbas%intcut, stmp, dstmp, &
                           & dstmpdqeffj, dstmpdqeffi)
                     else
                        call overlap_grad_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                           & bas%cgto(ish, isp)%raw, &
                           & acache%auxbas%cgto(jproj, jat), &
                           & bcache%cgto(ish, iat), r2, vec, &
                           & acp%auxbas%intcut, stmp, dstmp)
                     end if

                     nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
                     do iao = 1, msao(bas%cgto(ish, isp)%raw%ang)
                        do jao = 1, nao
                           ij = jao + nao*(iao-1)
                           dG = dG + response(jj+jao, ii+iao, itr, spin) &
                              & * dstmp(:, ij)
                           if (compute_qeff_grad) then
                              tmp = response(jj+jao, ii+iao, itr, spin) &
                                 & * dstmpdqeffi(ij)
                              dcnbasi = dcnbasi &
                                 & + bcache%cgto(ish, iat)%dqeffdcn*tmp
                              dqbasi = dqbasi &
                                 & + bcache%cgto(ish, iat)%dqeffdq*tmp
                           end if
                        end do
                     end do
                  end do
               end do

               gradient_local(:, iat) = gradient_local(:, iat) + dG
               gradient_local(:, jat) = gradient_local(:, jat) - dG
               sigma_local = sigma_local + 0.5_wp*(spread(vec, 1, 3) &
                  & * spread(dG, 2, 3) + spread(dG, 1, 3)*spread(vec, 2, 3))
            end do
         end do
         if (compute_qeff_grad) then
            dEdcnbas_local(iat) = dEdcnbas_local(iat) + dcnbasi
            dEdqbas_local(iat) = dEdqbas_local(iat) + dqbasi
         end if
      end do
      !$omp end do

      !$omp critical
      if (compute_qeff_grad) then
         dEdcnbas = dEdcnbas + dEdcnbas_local
         dEdqbas = dEdqbas + dEdqbas_local
      end if
      gradient = gradient + gradient_local
      sigma = sigma + sigma_local
      !$omp end critical
      if (compute_qeff_grad) then
         deallocate(dEdcnbas_local, dEdqbas_local)
      end if
      deallocate(gradient_local, sigma_local)
      !$omp end parallel
   end do
end subroutine get_pv_overlap_deriv_images


subroutine get_pv_overlap_deriv_0d(mol, bas, bcache, acp, acache, pmat, &
   & dEdcnbas, dEdqbas, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(in) :: acache
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :, :)
   !> Derivative of the electronic energy w.r.t. basis set coordination numbers
   real(wp), intent(inout) :: dEdcnbas(:)
   !> Derivative of the electronic energy w.r.t. basis set atomic charges
   real(wp), intent(inout) :: dEdqbas(:)
   !> Derivative of the electronic energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the electronic energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, isp, jsp, spin, nspin
   integer :: ish, jproj, is, js, ii, jj, iao, jao, nao, ij
   real(wp) :: r2, vec(3), dG(3), dcnbasi, dqbasi, tmp
   real(wp), allocatable :: stmp(:), dstmp(:, :), spmat(:, :)
   real(wp), allocatable :: dstmpdqeffi(:), dstmpdqeffj(:)
   logical :: compute_qeff_grad

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: dEdcnbas_local(:), dEdqbas_local(:), &
      & gradient_local(:, :), sigma_local(:, :)

   ! Determine before the loop if effective charge derivatives are needed
   compute_qeff_grad = bas%charge_dependent .or. acp%auxbas%charge_dependent

   nspin = size(pmat, 3)

   allocate(spmat(acp%auxbas%nao, bas%nao), stmp(msao(bas%maxl)*msao(acp%auxbas%maxl)), &
      & dstmp(3, msao(bas%maxl)*msao(acp%auxbas%maxl)))
   
   if (compute_qeff_grad) then
      allocate(dstmpdqeffi(msao(bas%maxl)*msao(acp%auxbas%maxl)), &
         & dstmpdqeffj(msao(bas%maxl)*msao(acp%auxbas%maxl)))
   end if

   do spin = 1, nspin
      ! Precalculate scaled projector-density matrix product
      call gemm(amat=acache%scaled_pv_overlap, bmat=pmat(:, :, spin), &
         & cmat=spmat(:, :), beta=0.0_wp)

      !$omp parallel default(none) shared(dEdcnbas, dEdqbas, gradient, sigma) &
      !$omp shared(mol, bas, bcache, acp, acache, spmat, compute_qeff_grad) &
      !$omp private(iat, jat, isp, jsp, is, js, ish, jproj, ii, jj, iao, jao, nao, ij) &
      !$omp private(r2, vec, stmp, dstmp, dstmpdqeffi, dstmpdqeffj, dcnbasi, dqbasi) &
      !$omp private(tmp, dG, dEdcnbas_local, dEdqbas_local, gradient_local, sigma_local)
      if (compute_qeff_grad) then
         allocate(dEdcnbas_local(size(dEdcnbas)), source=0.0_wp)
         allocate(dEdqbas_local(size(dEdqbas)), source=0.0_wp)
      end if
      allocate(gradient_local(size(gradient,1), size(gradient,2)), source=0.0_wp)
      allocate(sigma_local(size(sigma,1), size(sigma,2)), source=0.0_wp)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         isp = mol%id(iat)
         is = bas%ish_at(iat)
         dcnbasi = 0.0_wp
         dqbasi = 0.0_wp
         do jat = 1, mol%nat
            !if (iat == jat) cycle
            jsp = mol%id(jat)
            js = acp%auxbas%ish_at(jat)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            ! Loop over valence shells
            dG(:) = 0.0_wp
            do ish = 1, bas%nsh_id(isp)
               ii = bas%iao_sh(is+ish)
               ! Loop over projectors
               do jproj = 1, acp%auxbas%nsh_id(jsp)
                  jj = acp%auxbas%iao_sh(js+jproj)

                  ! Calculate overlap integral derivatives
                  if (compute_qeff_grad) then
                     call overlap_grad_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                        & bas%cgto(ish, isp)%raw, &
                        & acache%auxbas%cgto(jproj, jat), &
                        & bcache%cgto(ish, iat), r2, vec, &
                        & acp%auxbas%intcut, stmp, dstmp, &
                        & dstmpdqeffj, dstmpdqeffi)
                  else
                     call overlap_grad_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                        & bas%cgto(ish, isp)%raw, &
                        & acache%auxbas%cgto(jproj, jat), &
                        & bcache%cgto(ish, iat), r2, vec, &
                        & acp%auxbas%intcut, stmp, dstmp)
                  end if

                  ! Contract derivatives with the SxP intermediate
                  nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
                  do iao = 1, msao(bas%cgto(ish, isp)%raw%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        dG(:) = dG + 2.0_wp * spmat(jj+jao, ii+iao) * dstmp(:, ij)
                        if (compute_qeff_grad) then
                           tmp = 2.0_wp * spmat(jj+jao, ii+iao) * dstmpdqeffi(ij)
                           dcnbasi = dcnbasi + bcache%cgto(ish, iat)%dqeffdcn * tmp
                           dqbasi = dqbasi + bcache%cgto(ish, iat)%dqeffdq * tmp
                        end if
                     end do
                  end do

               end do
            end do
            ! Collect gradient and sigma once per atom pair
            gradient_local(:, iat) = gradient_local(:, iat) + dG
            gradient_local(:, jat) = gradient_local(:, jat) - dG
            sigma_local(:, :) = sigma_local + 0.5_wp * (spread(vec, 1, 3) &
               & * spread(dG, 2, 3) + spread(dG, 1, 3) * spread(vec, 2, 3))
         end do
         ! Collect effective charge derivatives once per atom
         if (compute_qeff_grad) then
            dEdcnbas_local(iat) = dEdcnbas_local(iat) + dcnbasi
            dEdqbas_local(iat) = dEdqbas_local(iat) + dqbasi
         end if
      end do
      !$omp end do
      !$omp critical
      if (compute_qeff_grad) then
         dEdcnbas(:) = dEdcnbas + dEdcnbas_local
         dEdqbas(:) = dEdqbas + dEdqbas_local
      end if
      gradient(:, :) = gradient + gradient_local
      sigma(:, :) = sigma + sigma_local
      !$omp end critical
      if (compute_qeff_grad) then
         deallocate(dEdcnbas_local, dEdqbas_local)
      end if
      deallocate(gradient_local, sigma_local)
      !$omp end parallel
   end do

end subroutine get_pv_overlap_deriv_0d


subroutine get_pv_overlap_deriv_3d(mol, trans, list, bas, bcache, acp, acache, &
   & pmat, dEdcnbas, dEdqbas, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(in) :: acache
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :, :)
   !> Derivative of the electronic energy w.r.t. basis set coordination numbers
   real(wp), intent(inout) :: dEdcnbas(:)
   !> Derivative of the electronic energy w.r.t. basis set atomic charges
   real(wp), intent(inout) :: dEdqbas(:)
   !> Derivative of the ACP energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the ACP energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, isp, jsp, spin, nspin, itr
   integer :: ish, jproj, is, js, ii, jj, iao, jao, nao, ij
   real(wp) :: r2, vec(3), dG(3), dcnbasi, dqbasi, tmp
   real(wp), allocatable :: stmp(:), dstmp(:, :), spmat(:, :)
   real(wp), allocatable :: dstmpdqeffi(:), dstmpdqeffj(:)
   logical :: compute_qeff_grad

   ! Thread-private array for reduction
   real(wp), allocatable :: dEdcnbas_local(:), dEdqbas_local(:), &
      & gradient_local(:, :), sigma_local(:, :)

   ! Determine before the loop if effective charge derivatives are needed
   compute_qeff_grad = bas%charge_dependent .or. acp%auxbas%charge_dependent

   nspin = size(pmat, 3)

   allocate(spmat(acp%auxbas%nao, bas%nao), stmp(msao(bas%maxl)*msao(acp%auxbas%maxl)), &
      & dstmp(3, msao(bas%maxl)*msao(acp%auxbas%maxl)))
   
   if (compute_qeff_grad) then
      allocate(dstmpdqeffi(msao(bas%maxl)*msao(acp%auxbas%maxl)), &
         & dstmpdqeffj(msao(bas%maxl)*msao(acp%auxbas%maxl)))
   end if

   do spin = 1, nspin
      ! Precalculate scaled projector-density matrix product
      call gemm(amat=acache%scaled_pv_overlap, bmat=pmat(:, :, spin), &
         & cmat=spmat(:, :), beta=0.0_wp)

      !$omp parallel default(none) shared(dEdcnbas, dEdqbas, gradient, sigma, trans) &
      !$omp shared(mol, bas, bcache, acp, acache, spmat, compute_qeff_grad) &
      !$omp private(iat, jat, itr, isp, jsp, is, js, ish, jproj, ii, jj, iao, jao, nao, ij) &
      !$omp private(r2, vec, stmp, dstmp, dstmpdqeffi, dstmpdqeffj, dcnbasi, dqbasi) &
      !$omp private(tmp, dG, dEdcnbas_local, dEdqbas_local, gradient_local, sigma_local)
      if (compute_qeff_grad) then
         allocate(dEdcnbas_local(size(dEdcnbas)), source=0.0_wp)
         allocate(dEdqbas_local(size(dEdqbas)), source=0.0_wp)
      end if
      allocate(gradient_local(size(gradient,1), size(gradient,2)), source=0.0_wp)
      allocate(sigma_local(size(sigma,1), size(sigma,2)), source=0.0_wp)
      
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         isp = mol%id(iat)
         is = bas%ish_at(iat)
         dcnbasi = 0.0_wp
         dqbasi = 0.0_wp
         
         do jat = 1, mol%nat
            jsp = mol%id(jat)
            js = acp%auxbas%ish_at(jat)
            
            ! Loop over translations
            do itr = 1, size(trans, 2)
               vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
               r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
               
               dG(:) = 0.0_wp
               
               ! Loop over valence shells
               do ish = 1, bas%nsh_id(isp)
                  ii = bas%iao_sh(is+ish)
                  ! Loop over projectors
                  do jproj = 1, acp%auxbas%nsh_id(jsp)
                     jj = acp%auxbas%iao_sh(js+jproj)

                     ! Calculate overlap integral derivatives
                     if (compute_qeff_grad) then
                        call overlap_grad_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                           & bas%cgto(ish, isp)%raw, &
                           & acache%auxbas%cgto(jproj, jat), &
                           & bcache%cgto(ish, iat), r2, vec, &
                           & acp%auxbas%intcut, stmp, dstmp, &
                           & dstmpdqeffj, dstmpdqeffi)
                     else
                        call overlap_grad_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                           & bas%cgto(ish, isp)%raw, &
                           & acache%auxbas%cgto(jproj, jat), &
                           & bcache%cgto(ish, iat), r2, vec, &
                           & acp%auxbas%intcut, stmp, dstmp)
                     end if

                     ! Contract derivatives with the SxP intermediate
                     nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
                     do iao = 1, msao(bas%cgto(ish, isp)%raw%ang)
                        do jao = 1, nao
                           ij = jao + nao*(iao-1)
                           dG(:) = dG + 2.0_wp * spmat(jj+jao, ii+iao) * dstmp(:, ij)
                           if (compute_qeff_grad) then
                              tmp = 2.0_wp * spmat(jj+jao, ii+iao) * dstmpdqeffi(ij)
                              dcnbasi = dcnbasi + bcache%cgto(ish, iat)%dqeffdcn * tmp
                              dqbasi = dqbasi + bcache%cgto(ish, iat)%dqeffdq * tmp
                           end if
                        end do
                     end do

                  end do
               end do
               
               ! Collect gradient and sigma once per atom-pair translation
               ! this implicitly also handles self-images correctly
               gradient_local(:, iat) = gradient_local(:, iat) + dG
               gradient_local(:, jat) = gradient_local(:, jat) - dG
               sigma_local(:, :) = sigma_local + 0.5_wp * (spread(vec, 1, 3) &
                  & * spread(dG, 2, 3) + spread(dG, 1, 3) * spread(vec, 2, 3))
               
            end do ! itr
         end do ! jat
         
         ! Collect effective charge derivatives once per atom
         if (compute_qeff_grad) then
            dEdcnbas_local(iat) = dEdcnbas_local(iat) + dcnbasi
            dEdqbas_local(iat) = dEdqbas_local(iat) + dqbasi
         end if
      end do
      !$omp end do
      
      !$omp critical
      if (compute_qeff_grad) then
         dEdcnbas(:) = dEdcnbas + dEdcnbas_local
         dEdqbas(:) = dEdqbas + dEdqbas_local
      end if
      gradient(:, :) = gradient + gradient_local
      sigma(:, :) = sigma + sigma_local
      !$omp end critical
      
      if (compute_qeff_grad) then
         deallocate(dEdcnbas_local, dEdqbas_local)
      end if
      deallocate(gradient_local, sigma_local)
      !$omp end parallel
   end do

end subroutine get_pv_overlap_deriv_3d


end module tblite_acp_type
