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

!> @file tblite/exchange/fock.f90
!> Provides matrix multiplication algorithm for approximated Fock exchange

!> Approximated Fock exchange based on matrix multiplication algorithm
module tblite_exchange_fock
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_blas, only: gemm, symm, axpy
   use tblite_basis_type, only : basis_type
   use tblite_container_cache, only : container_cache
   use tblite_cutoff, only : get_lattice_points
   use tblite_exchange_cache, only : exchange_bvk_kernel, exchange_cache
   use tblite_exchange_type, only : exchange_type
   use tblite_utils_average, only : average_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wignerseitz, only : wignerseitz_cell, &
      & get_wignerseitz_pairs, wignerseitz_threshold
   implicit none
   private

   public :: new_exchange_fock

   public :: exchange_bvk_kernel

   type, public, extends(exchange_type) :: exchange_fock
      !> Averaged Hubbard parameter for each shell and species
      real(wp), allocatable :: hubbard(:, :, :, :)
      !> One center exchange integrals
      real(wp), allocatable :: onecxints(:, :, :)
      !> Diagonal scaling of the Fock exchange
      real(wp) :: ondiag_scale
      !> Off-diagonal scaling of the Fock exchange
      real(wp), allocatable :: offdiag_scale(:, :, :, :)
      !> Exponent of radius dependent hubbard scaling
      real(wp) :: hubbard_exp
      !> Radius prefactor of radius dependent hubbard scaling
      real(wp) :: hubbard_exp_r0
      !> Smoothening exponent (1 = Mataga-Nishimoto, 2 = Klopman-Ohno)
      real(wp) :: gexp
      !> Pairwise radii for approximate exchange integrals
      real(wp), allocatable :: rad(:, :)
      !> Charge-dependence of the onsite Fock exchange
      real(wp), allocatable :: kq(:, :)
      !> Bond-order correlation scaling factor for each atom pair
      real(wp), allocatable :: corr_scale(:, :)
      !> Bond-order correlation damping exponent
      real(wp) :: corr_exp
      !> Bond-order correlation radius for each atom pair
      real(wp), allocatable :: corr_rad(:, :)
   contains
      !> Evaluate Mulliken Fock exchange gamma matrix
      procedure :: get_mulliken_Kmatrix
      !> Evaluate onsite Fock exchange gamma matrix
      procedure :: get_onsite_Kmatrix
      !> Evaluate bond-order correlation correction gamma matrix
      procedure :: get_bocorr_Kmatrix
      !> Evaluate image-resolved exchange kernels on a regular BvK mesh
      procedure :: get_bvk_Kmatrix
      !> Check/set the exact model signature of a persistent BvK plan
      procedure :: bvk_model_matches
      procedure :: set_bvk_model_signature
      !> Check a complete cached plan, including kernel representative order
      procedure :: bvk_plan_matches
      !> Contract image-resolved kernel responses with geometry derivatives
      procedure :: get_bvk_Kmatrix_derivs
      !> Evaluate the gradient of the Mulliken exchange energy
      procedure :: get_mulliken_derivs
      !> Evaluate Mulliken derivatives from direct independent-kernel responses
      procedure :: get_mulliken_derivs_direct
      !> Evaluate the gradient of the bond-order correlation correction energy
      procedure :: get_bocorr_derivs
      !> Evaluate bond-order derivatives from direct independent-kernel responses
      procedure :: get_bocorr_derivs_direct
      !> Calculate exchange contribution to the Fock matrix
      procedure :: get_KFock
      !> Calculate exchange contribution for a complex Hermitian k-point block
      procedure :: get_KFock_kpoint
      !> Calculate exchange contribution for a complete regular k-point mesh
      procedure :: get_KFock_kmesh
      !> Apply/reconstruct the memory-reduced BvK streaming forward path
      procedure :: get_KFock_stream_apply
      procedure :: get_KFock_stream_block
      !> Calculate exchange response for a complex Hermitian k-point block
      procedure :: get_KGrad_kpoint
      !> Calculate exchange response for a complete regular k-point mesh
      procedure :: get_KGrad_kmesh
      !> Calculate exchange contribution to the gradient
      procedure :: get_KGrad
   end type exchange_fock

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   character(len=*), parameter :: label = "Mulliken and onsite Fock exchange + bond-order correlation"

contains


!> Matrix product for contiguous complex arrays.
!>
!> Keeping the MATMUL behind an assumed-shape contiguous boundary prevents
!> compilers from packing the rank-remapped k-mesh views below.
pure subroutine matmul_complex_contiguous(amat, bmat, cmat)
   complex(wp), contiguous, intent(in) :: amat(:, :), bmat(:, :)
   complex(wp), contiguous, intent(out) :: cmat(:, :)

   cmat = matmul(amat, bmat)
end subroutine matmul_complex_contiguous


!> Create a new approximate exchange container
subroutine new_exchange_fock(self, mol, bas, hubbard, hubbard_average, &
   & avg_exponents, ondiag_scale, offdiag_scale, offdiag_average, hubbard_exp, &
   & hubbard_exp_r0, rad, gexp, onecxints, kq, corr_scale, corr_scale_average, &
   & corr_exp, corr_rad, corr_rad_average, frscale, omega, lrscale)
   !> Instance of the Fock exchange container
   type(exchange_fock), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Description of the basis set
   type(basis_type), intent(in) :: bas
   !> Hubbard parameter for all shells and species
   real(wp), intent(in) :: hubbard(:, :)
   !> Averaging function for Hubbard parameter of a shell-pair
   type(average_type), intent(in) :: hubbard_average
   !> Averaging exponents for all shells and species
   real(wp), intent(in) :: avg_exponents(:, :)
   !> Diagonal scaling of the Fock exchange
   real(wp), intent(in) :: ondiag_scale
   !> Off-diagonal scaling of the Fock exchange
   real(wp), intent(in) :: offdiag_scale(:, :)
   !> Averaging function for the off-diagonal scaling of a shell-pair
   type(average_type), intent(in) :: offdiag_average
   !> Exponent of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Radius prefactor of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp_r0
   !> Radius for hubbard scaling
   real(wp), intent(in) :: rad(:, :)
   !> Smoothening exponent; 1 = Mataga-Nishimoto, 2 = Klopman-Ohno
   real(wp), intent(in) :: gexp
   !> One center exchange integrals
   real(wp), intent(in) :: onecxints(:, :, :)
   !> Charge dependence of the onsite Fock exchange 
   real(wp), intent(in) :: kq(:, :)
   !> Bond-order correlation scaling factor for each atoa
   real(wp), intent(in) :: corr_scale(:)
   !> Averaging function for correlation scaling of an atom-pair
   type(average_type), intent(in) :: corr_scale_average
   !> Bond-order correlation damping exponent
   real(wp), intent(in) :: corr_exp
   !> Bond-order correlation radius for each atom
   real(wp), intent(in) :: corr_rad(:)
   !> Averaging function for correlation radius of an atom-pair
   type(average_type), intent(in) :: corr_rad_average
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Optional range separation parameter
   real(wp), intent(in), optional :: omega
   !> Optional long range scaling for range seperated exchange
   real(wp), intent(in), optional :: lrscale

   integer :: isp, jsp, ish, jsh
   real(wp) :: inter_avg

   self%label = label

   allocate(self%nsh_id(mol%nid), self%nao_sh(bas%nsh), self%ish_at(mol%nat), &
      & self%iao_sh(bas%nsh))
   self%nao = bas%nao
   self%nsh = bas%nsh
   self%nsh_id = bas%nsh_id
   self%nao_sh = bas%nao_sh
   self%ish_at = bas%ish_at
   self%iao_sh = bas%iao_sh
   self%maxsh = maxval(bas%nsh_id)

   ! Global and optional range-separated exchange scaling
   self%frscale = frscale
   if (present(omega).and.present(lrscale)) then
      self%omega = omega
      self%lrscale = lrscale
   end if

   ! Diagonal and off-diagonal scaling in the Mulliken exchange matrix
   allocate(self%offdiag_scale(self%maxsh, self%maxsh, mol%nid, mol%nid))
   self%offdiag_scale = 0.0_wp
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         do ish = 1, bas%nsh_id(isp)
            do jsh = 1, bas%nsh_id(jsp)
               self%offdiag_scale(jsh, ish, jsp, isp) = offdiag_average%value( &
                  & offdiag_scale(ish, isp), offdiag_scale(jsh, jsp))
            end do
         end do
      end do
   end do
   self%ondiag_scale  = ondiag_scale

   ! Error-function damped Klopman-Ohno-Mataga kernel for Mulliken exchange
   self%gexp = gexp
   self%hubbard_exp = hubbard_exp
   self%hubbard_exp_r0 = hubbard_exp_r0
   self%rad = rad

   ! Pairwise averaged Hubbard parameters for Mulliken exchange
   allocate(self%hubbard(self%maxsh, self%maxsh, mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%hubbard(:, :, jsp, isp) = 0.0_wp
         do ish = 1, bas%nsh_id(isp)
            do jsh = 1, bas%nsh_id(jsp)
               inter_avg = max(avg_exponents(ish, isp), avg_exponents(jsh, jsp))
               self%hubbard(jsh, ish, jsp, isp) = hubbard_average%value(&
                  & hubbard(ish, isp), hubbard(jsh, jsp), inter_avg)
            end do
         end do
      end do
   end do

   ! Charge-dependence and onsite exchange integrals for onsite Fock exchange
   self%kq = kq
   self%onecxints = onecxints

   ! Bond-order based correlation scaling and radii
   allocate(self%corr_scale(mol%nid, mol%nid), self%corr_rad(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%corr_scale(isp, jsp) = &
            & corr_scale_average%value(corr_scale(isp), corr_scale(jsp))
         self%corr_rad(isp, jsp) = &
            & corr_rad_average%value(corr_rad(isp), corr_rad(jsp))
         end do 
   end do 
   self%corr_exp = corr_exp

end subroutine new_exchange_fock


!> Check whether a persistent BvK plan belongs to this exact exchange model.
logical function bvk_model_matches(self, cache) result(matches)
   class(exchange_fock), intent(in) :: self
   type(exchange_cache), intent(in) :: cache

   matches = .false.
   if (.not.cache%bvk_model_valid) return
   if (.not.allocated(cache%bvk_model_nsh_id) &
      & .or. .not.allocated(cache%bvk_model_nao_sh) &
      & .or. .not.allocated(cache%bvk_model_ish_at) &
      & .or. .not.allocated(cache%bvk_model_iao_sh) &
      & .or. .not.allocated(cache%bvk_model_hubbard) &
      & .or. .not.allocated(cache%bvk_model_offdiag_scale) &
      & .or. .not.allocated(cache%bvk_model_rad) &
      & .or. .not.allocated(cache%bvk_model_corr_scale) &
      & .or. .not.allocated(cache%bvk_model_corr_rad)) return
   if (cache%bvk_model_nao /= self%nao &
      & .or. cache%bvk_model_nsh /= self%nsh &
      & .or. cache%bvk_model_maxsh /= self%maxsh) return
   if (size(cache%bvk_model_nsh_id) /= size(self%nsh_id) &
      & .or. size(cache%bvk_model_nao_sh) /= size(self%nao_sh) &
      & .or. size(cache%bvk_model_ish_at) /= size(self%ish_at) &
      & .or. size(cache%bvk_model_iao_sh) /= size(self%iao_sh)) return
   if (any(shape(cache%bvk_model_hubbard) /= shape(self%hubbard)) &
      & .or. any(shape(cache%bvk_model_offdiag_scale) /= &
      & shape(self%offdiag_scale)) &
      & .or. any(shape(cache%bvk_model_rad) /= shape(self%rad)) &
      & .or. any(shape(cache%bvk_model_corr_scale) /= &
      & shape(self%corr_scale)) &
      & .or. any(shape(cache%bvk_model_corr_rad) /= &
      & shape(self%corr_rad))) return
   if (any(cache%bvk_model_nsh_id /= self%nsh_id) &
      & .or. any(cache%bvk_model_nao_sh /= self%nao_sh) &
      & .or. any(cache%bvk_model_ish_at /= self%ish_at) &
      & .or. any(cache%bvk_model_iao_sh /= self%iao_sh)) return
   if (cache%bvk_model_frscale /= self%frscale &
      & .or. cache%bvk_model_omega /= self%omega &
      & .or. cache%bvk_model_lrscale /= self%lrscale &
      & .or. cache%bvk_model_ondiag_scale /= self%ondiag_scale &
      & .or. cache%bvk_model_hubbard_exp /= self%hubbard_exp &
      & .or. cache%bvk_model_hubbard_exp_r0 /= self%hubbard_exp_r0 &
      & .or. cache%bvk_model_gexp /= self%gexp &
      & .or. cache%bvk_model_corr_exp /= self%corr_exp) return
   if (any(cache%bvk_model_hubbard /= self%hubbard) &
      & .or. any(cache%bvk_model_offdiag_scale /= self%offdiag_scale) &
      & .or. any(cache%bvk_model_rad /= self%rad) &
      & .or. any(cache%bvk_model_corr_scale /= self%corr_scale) &
      & .or. any(cache%bvk_model_corr_rad /= self%corr_rad)) return
   matches = .true.
end function bvk_model_matches


!> Store the exact static exchange-model inputs used by a BvK kernel.
subroutine set_bvk_model_signature(self, cache)
   class(exchange_fock), intent(in) :: self
   type(exchange_cache), intent(inout) :: cache

   cache%bvk_model_valid = .false.
   cache%bvk_model_nao = self%nao
   cache%bvk_model_nsh = self%nsh
   cache%bvk_model_maxsh = self%maxsh
   cache%bvk_model_frscale = self%frscale
   cache%bvk_model_omega = self%omega
   cache%bvk_model_lrscale = self%lrscale
   cache%bvk_model_ondiag_scale = self%ondiag_scale
   cache%bvk_model_hubbard_exp = self%hubbard_exp
   cache%bvk_model_hubbard_exp_r0 = self%hubbard_exp_r0
   cache%bvk_model_gexp = self%gexp
   cache%bvk_model_corr_exp = self%corr_exp
   if (allocated(cache%bvk_model_nsh_id)) deallocate(cache%bvk_model_nsh_id)
   if (allocated(cache%bvk_model_nao_sh)) deallocate(cache%bvk_model_nao_sh)
   if (allocated(cache%bvk_model_ish_at)) deallocate(cache%bvk_model_ish_at)
   if (allocated(cache%bvk_model_iao_sh)) deallocate(cache%bvk_model_iao_sh)
   if (allocated(cache%bvk_model_hubbard)) deallocate(cache%bvk_model_hubbard)
   if (allocated(cache%bvk_model_offdiag_scale)) &
      & deallocate(cache%bvk_model_offdiag_scale)
   if (allocated(cache%bvk_model_rad)) deallocate(cache%bvk_model_rad)
   if (allocated(cache%bvk_model_corr_scale)) &
      & deallocate(cache%bvk_model_corr_scale)
   if (allocated(cache%bvk_model_corr_rad)) deallocate(cache%bvk_model_corr_rad)
   allocate(cache%bvk_model_nsh_id, source=self%nsh_id)
   allocate(cache%bvk_model_nao_sh, source=self%nao_sh)
   allocate(cache%bvk_model_ish_at, source=self%ish_at)
   allocate(cache%bvk_model_iao_sh, source=self%iao_sh)
   allocate(cache%bvk_model_hubbard, source=self%hubbard)
   allocate(cache%bvk_model_offdiag_scale, source=self%offdiag_scale)
   allocate(cache%bvk_model_rad, source=self%rad)
   allocate(cache%bvk_model_corr_scale, source=self%corr_scale)
   allocate(cache%bvk_model_corr_rad, source=self%corr_rad)
   cache%bvk_model_valid = .true.
end subroutine set_bvk_model_signature


!> Check the model, mesh signature, and representative ordering of a plan.
logical function bvk_plan_matches(self, mol, cache, kernel, kfrac, weights) &
   & result(matches)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(in) :: cache
   type(exchange_bvk_kernel), intent(in) :: kernel
   real(wp), intent(in) :: kfrac(:, :), weights(:)

   matches = .false.
   if (.not.cache%bvk_matches(mol, kernel%nmesh, kfrac, weights)) return
   if (.not.self%bvk_model_matches(cache)) return
   if (.not.allocated(kernel%reps)) return
   if (any(shape(kernel%reps) /= shape(cache%bvk_kernel%reps))) return
   if (any(kernel%reps /= cache%bvk_kernel%reps)) return
   matches = .true.
end function bvk_plan_matches


!> Build image-resolved exchange kernels for a regular Born-von Karman mesh
!>
!> This routine evaluates the Wigner-Seitz kernels of the corresponding BvK
!> supercell without constructing that supercell.  Images are indexed by the
!> representatives ``0 <= reps(:, icell) < nmesh``.  The returned block has
!> its row in the home cell and its column in the represented cell, i.e.
!> ``G(k) = sum_R exp(+i k.R) G_R``.
subroutine get_bvk_Kmatrix(self, mol, nmesh, kernel)
   !> Instance of the Fock exchange container
   class(exchange_fock), intent(in) :: self
   !> Primitive-cell structure
   type(structure_type), intent(in) :: mol
   !> Full regular k-point mesh; non-periodic directions must have extent one
   integer, intent(in) :: nmesh(3)
   !> Image-resolved exchange kernels
   type(exchange_bvk_kernel), intent(out) :: kernel

   integer :: icell, ix, iy, iz, iat, jat, izp, jzp, is, js
   integer :: ish, jsh, img, nimg, ncell, idim
   integer, allocatable :: tridx(:)
   real(wp) :: bvk_lattice(3, 3), cellvec(3), vec(3), r1
   real(wp) :: rsh, gam, scale, wsw, arg, damp, corr
   real(wp), allocatable :: trans(:, :)

   kernel%nmesh = nmesh
   ncell = product(nmesh)
   allocate(kernel%reps(3, ncell))
   allocate(kernel%g_mulliken_r(self%nsh, self%nsh, ncell), &
      & kernel%g_bocorr_r(mol%nat, mol%nat, ncell), source=0.0_wp)

   icell = 0
   do iz = 0, nmesh(3)-1
      do iy = 0, nmesh(2)-1
         do ix = 0, nmesh(1)-1
            icell = icell + 1
            kernel%reps(:, icell) = [ix, iy, iz]
         end do
      end do
   end do

   bvk_lattice = mol%lattice
   do idim = 1, 3
      bvk_lattice(:, idim) = bvk_lattice(:, idim)*real(nmesh(idim), wp)
   end do
   call get_lattice_points(mol%periodic, bvk_lattice, &
      & wignerseitz_threshold, trans)
   allocate(tridx(size(trans, 2)))

   do icell = 1, ncell
      cellvec = matmul(mol%lattice, real(kernel%reps(:, icell), wp))
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = self%ish_at(iat)
         do jat = 1, mol%nat
            jzp = mol%id(jat)
            js = self%ish_at(jat)

            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - cellvec
            call get_wignerseitz_pairs(nimg, trans, vec, tridx)
            if (nimg > 0) then
               wsw = 1.0_wp/real(nimg, wp)
               do img = 1, nimg
                  r1 = norm2(vec - trans(:, tridx(img)))
                  if (r1 < eps) cycle

                  do ish = 1, self%nsh_id(izp)
                     do jsh = 1, self%nsh_id(jzp)
                        rsh = wsw*get_gmulliken_pair(r1, &
                           & self%hubbard(jsh, ish, jzp, izp), &
                           & self%offdiag_scale(jsh, ish, jzp, izp), &
                           & self%hubbard_exp, self%hubbard_exp_r0, &
                           & self%rad(izp, jzp), self%gexp, self%frscale, &
                           & self%lrscale, self%omega)
                        kernel%g_mulliken_r(is+ish, js+jsh, icell) = &
                           & kernel%g_mulliken_r(is+ish, js+jsh, icell) + rsh
                     end do
                  end do

                  arg = self%corr_exp*(r1-self%corr_rad(izp, jzp)) &
                     & /self%rad(izp, jzp)
                  damp = 0.5_wp*(1.0_wp + erf(-arg))
                  corr = wsw*self%corr_scale(izp, jzp)*damp
                  kernel%g_bocorr_r(iat, jat, icell) = &
                     & kernel%g_bocorr_r(iat, jat, icell) + corr
               end do
            end if

            ! The analytic one-center Mulliken term belongs to R=0 only.
            ! Periodic self-image terms were included by the pair loop above.
            if (icell == 1 .and. iat == jat) then
               do ish = 1, self%nsh_id(izp)
                  do jsh = 1, ish-1
                     scale = self%offdiag_scale(jsh, ish, izp, izp)
                     gam = self%hubbard(jsh, ish, izp, izp)*scale &
                        & *self%frscale
                     kernel%g_mulliken_r(is+jsh, is+ish, icell) = &
                        & kernel%g_mulliken_r(is+jsh, is+ish, icell) + gam
                     kernel%g_mulliken_r(is+ish, is+jsh, icell) = &
                        & kernel%g_mulliken_r(is+ish, is+jsh, icell) + gam
                  end do
                  gam = self%hubbard(ish, ish, izp, izp) &
                     & *self%ondiag_scale*self%frscale
                  kernel%g_mulliken_r(is+ish, is+ish, icell) = &
                     & kernel%g_mulliken_r(is+ish, is+ish, icell) + gam
               end do
            end if
         end do
      end do
   end do

end subroutine get_bvk_Kmatrix


!> Contract image-resolved BvK kernel adjoints with geometry derivatives.
!>
!> The response matrices use the same oriented home-row/translated-column
!> convention as get_bvk_Kmatrix and must not be symmetrized between images.
!> The returned stress-like tensor is the positive homogeneous-strain
!> derivative, ``sigma(a,b) = dE/d epsilon(a,b)``; no volume factor is applied.
subroutine get_bvk_Kmatrix_derivs(self, mol, kernel, mulliken_grad_r, &
   & bocorr_grad_r, gradient, sigma)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_bvk_kernel), intent(in) :: kernel
   real(wp), intent(in) :: mulliken_grad_r(:, :, :), bocorr_grad_r(:, :, :)
   real(wp), intent(inout) :: gradient(:, :), sigma(:, :)

   integer :: icell, iat, jat, izp, jzp, is, js, ish, jsh
   integer :: iao, jao, ii, jj, ni, nj, img, nimg, idim
   integer, allocatable :: tridx(:)
   real(wp) :: atom_grad, bvk_lattice(3, 3), cellvec(3), coeff
   real(wp) :: dcorr, dkernel, dpair(3), pairvec(3), r1, shell_grad, wsw
   real(wp) :: arg
   real(wp), allocatable :: trans(:, :)

   bvk_lattice = mol%lattice
   do idim = 1, 3
      bvk_lattice(:, idim) = bvk_lattice(:, idim) &
         & *real(kernel%nmesh(idim), wp)
   end do
   call get_lattice_points(mol%periodic, bvk_lattice, &
      & wignerseitz_threshold, trans)
   allocate(tridx(size(trans, 2)))

   do icell = 1, size(kernel%reps, 2)
      cellvec = matmul(mol%lattice, real(kernel%reps(:, icell), wp))
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = self%ish_at(iat)
         do jat = 1, mol%nat
            jzp = mol%id(jat)
            js = self%ish_at(jat)

            pairvec = mol%xyz(:, iat) - mol%xyz(:, jat) - cellvec
            call get_wignerseitz_pairs(nimg, trans, pairvec, tridx)
            if (nimg <= 0) cycle
            wsw = 1.0_wp/real(nimg, wp)

            do img = 1, nimg
               pairvec = mol%xyz(:, iat) - mol%xyz(:, jat) - cellvec &
                  & - trans(:, tridx(img))
               r1 = norm2(pairvec)
               if (r1 < eps) cycle

               coeff = 0.0_wp
               do ish = 1, self%nsh_id(izp)
                  ii = self%iao_sh(is+ish)
                  ni = self%nao_sh(is+ish)
                  do jsh = 1, self%nsh_id(jzp)
                     jj = self%iao_sh(js+jsh)
                     nj = self%nao_sh(js+jsh)
                     shell_grad = 0.0_wp
                     do iao = 1, ni
                        do jao = 1, nj
                           shell_grad = shell_grad &
                              & + mulliken_grad_r(ii+iao, jj+jao, icell)
                        end do
                     end do
                     call get_gmulliken_pair_deriv(r1, &
                        & self%hubbard(jsh, ish, jzp, izp), &
                        & self%offdiag_scale(jsh, ish, jzp, izp), &
                        & self%hubbard_exp, self%hubbard_exp_r0, &
                        & self%rad(izp, jzp), self%gexp, self%frscale, &
                        & self%lrscale, self%omega, dkernel)
                     coeff = coeff + wsw*shell_grad*dkernel
                  end do
               end do

               atom_grad = 0.0_wp
               do ish = 1, self%nsh_id(izp)
                  ii = self%iao_sh(is+ish)
                  ni = self%nao_sh(is+ish)
                  do jsh = 1, self%nsh_id(jzp)
                     jj = self%iao_sh(js+jsh)
                     nj = self%nao_sh(js+jsh)
                     do iao = 1, ni
                        do jao = 1, nj
                           atom_grad = atom_grad &
                              & + bocorr_grad_r(ii+iao, jj+jao, icell)
                        end do
                     end do
                  end do
               end do
               arg = self%corr_exp*(r1-self%corr_rad(izp, jzp)) &
                  & /self%rad(izp, jzp)
               dcorr = -self%corr_scale(izp, jzp)*self%corr_exp &
                  & *exp(-arg**2)/(sqrtpi*self%rad(izp, jzp))
               coeff = coeff + wsw*atom_grad*dcorr

               dpair = coeff*pairvec/r1
               gradient(:, iat) = gradient(:, iat) + dpair
               gradient(:, jat) = gradient(:, jat) - dpair
               sigma = sigma + spread(dpair, 1, 3)*spread(pairvec, 2, 3)
            end do
         end do
      end do
   end do

end subroutine get_bvk_Kmatrix_derivs


!> Calculate exchange contribution to the Fock matrix and atomic potential
subroutine get_KFock(self, mol, cache, density, overlap)
   !> Instance of the Fock exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container with intermediates and the final Fock matrix/potential
   type(exchange_cache), intent(inout) :: cache
   !> Density matrix
   real(wp), intent(in) :: density(:, :, :)
   !> Overlap matrix
   real(wp), intent(in) :: overlap(:, :)

   integer :: spin, iat, izp, ish, is, ii, iao, jsh, jj, jao
   real(wp) :: gfx, gri, tmp
   real(wp), allocatable :: tmpA(:,:), tmpB(:,:)
   real(wp), allocatable :: diagP(:), gdiagP(:)
   real(wp), allocatable :: tmpSP(:), gdiagSP(:)
   real(wp), allocatable :: tmpSPS(:), gdiagSPS(:)

   real(wp) :: spin_factor

   allocate(tmpA(self%nao, self%nao), tmpB(self%nao, self%nao), diagP(self%nao), &
      & gdiagP(self%nao), tmpSP(self%nao), gdiagSP(self%nao), tmpSPS(self%nao), &
      & gdiagSPS(self%nao), source = 0.0_wp)

   ! Select spin factor to cancel the quadratic dependence of the exchange energy
   ! on the occupation numbers (0.5 for restricted, and 1.0 for unrestricted)
   spin_factor = 0.5_wp
   if(size(density, 3) .gt. 1) then
      spin_factor = 1.0_wp
   end if

   cache%prev_F(:, :, :) = 0.0_wp
   cache%prev_vsh(:, :) = 0.0_wp

   ! Evaluate the Fock matrix contribution for Mulliken
   ! and onsite approximated Fock exchange for a symmetric density matrix
   do spin = 1, size(density, 3)

      ! Intermediate A = S x P
      call symm(amat=overlap, bmat=density(:, :, spin), cmat=tmpA)

      ! Collect P diagonal onsite correction
      do iao = 1, self%nao
         diagP(iao) = density(iao, iao, spin)
      end do
      call onsite_fx_symv(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_onsfx, diagP, gdiagP)

      ! Collect S x P diagonal onsite correction
      do iao = 1, self%nao
         tmpSP(iao) = tmpA(iao, iao)
      end do
      call onsite_fx_symv(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_onsfx, tmpSP, gdiagSP)

      ! Collect S * P summation
      tmpSP = 0.0_wp
      tmpSP(:) = sum(density(:, :, spin) * overlap, dim=2)

      ! Collect S x P x S diagonal onsite correction
      tmpB(:, :) = tmpA * overlap
      tmpSPS(:) = sum(tmpB, dim=2)
      call onsite_fx_symv(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_onsfx, tmpSPS, gdiagSPS)

      ! Intermediate F = 1/2 * (S x P) x S
      call gemm(amat=tmpA, bmat=overlap, alpha = 0.5_wp, &
         & cmat=cache%prev_F(:, :, spin))

      ! Shell potential from self-consistent onsite term
      !$omp parallel do default(none) schedule(runtime) shared(mol, self, cache) &
      !$omp shared(density, tmpA, diagP, tmpSP, tmpSPS, spin, spin_factor) &
      !$omp private(iat, izp, is, ish, jsh, ii, jj, iao, jao, gfx, gri)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is  = self%ish_at(iat)
         do ish = 1, self%nsh_id(izp)
            ii = self%iao_sh(is+ish)
            do jsh = 1, self%nsh_id(izp)
               jj = self%iao_sh(is+jsh)

               gfx = 0.0_wp
               gri = 0.0_wp
               do iao = 1, self%nao_sh(is+ish)
                  do jao = 1, self%nao_sh(is+jsh)
                     gfx = gfx - spin_factor * 0.25_wp * ( &
                        & + (density(jj+jao, ii+iao, spin) * cache%prev_F(jj+jao, ii+iao, spin)) &
                        & + 0.5_wp * (tmpA(ii+iao, jj+jao) * tmpA(ii+iao, jj+jao)) &
                        & + 0.25_wp * (tmpSPS(jj+jao) * diagP(ii+iao)) &
                        & + 0.5_wp * (tmpSP(jj+jao) * tmpSP(ii+iao)) &
                        & + 0.25_wp * (diagP(jj+jao) * tmpSPS(ii+iao)) )

                     if (ish == jsh) then
                        gri = gri + spin_factor * ( &
                           & + (density(jj+jao, ii+iao, spin) * cache%prev_F(jj+jao, ii+iao, spin)) &
                           & + 0.5_wp * (tmpA(ii+iao, jj+jao) * tmpA(jj+jao, ii+iao)) )
                     end if
                  end do
               end do

               if (ish == jsh) then
                  cache%prev_vsh(is+ish, 1) = cache%prev_vsh(is+ish, 1) &
                     & + cache%dgdq_onsri(ish, is+ish) * gri & 
                     & + cache%dgdq_onsfx(jsh, ish, is+ish) * gfx
               else
                  cache%prev_vsh(is+ish, 1) = cache%prev_vsh(is+ish, 1) + &
                     & cache%dgdq_onsfx(jsh, ish, is+ish) * gfx
                  cache%prev_vsh(is+jsh, 1) = cache%prev_vsh(is+jsh, 1) + &
                     & cache%dgdq_onsfx(jsh, ish, is+jsh) * gfx
               end if
            end do
         end do
      end do

      ! Apply Mulliken, onsite, and bond-order correlation matrices as g * (S x P)
      tmpB = 0.0_wp
      call shell_hadamard_add(self%nsh, self%nao_sh, self%iao_sh, cache%g_mulliken, &
         & tmpA, 1.0_wp, tmpB)
      call atom_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_bocorr, tmpA, -4.0_wp, tmpB)
      call onsite_fx_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsfx, tmpA, 0.5_wp, tmpB, trans_src=.true.)
      call onsite_ri_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsri, tmpA, -2.0_wp, tmpB)
      tmpA = tmpB

      ! Apply Mulliken, onsite, and bond-order correlation matrices as g * (S x P x S)
      tmpB = 0.0_wp
      call onsite_fx_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsfx, cache%prev_F(:, :, spin), &
         & 0.5_wp, tmpB)
      call onsite_ri_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsri, cache%prev_F(:, :, spin), &
         & -2.0_wp, tmpB)
      call shell_hadamard_add(self%nsh, self%nao_sh, self%iao_sh, cache%g_mulliken, &
         & cache%prev_F(:, :, spin), 1.0_wp, tmpB)
      cache%prev_F(:, :, spin) = tmpB

      ! Apply Mulliken, onsite, and bond-order correlation matrices as g * P
      tmpB = 0.0_wp
      call shell_hadamard_add(self%nsh, self%nao_sh, self%iao_sh, cache%g_mulliken, &
         & density(:, :, spin), 1.0_wp, tmpB)
      call onsite_fx_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsfx, density(:, :, spin), &
         & 0.5_wp, tmpB)
      call onsite_ri_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsri, density(:, :, spin), &
         & -2.0_wp, tmpB)

      ! Add P diagonal onsite correction
      do iao = 1, self%nao
         call axpy(xvec=overlap(:, iao), yvec=tmpA(:, iao), &
            & alpha=0.25_wp*gdiagP(iao))
      end do 

      ! Intermediate A += 1/2 * S x (g * P)
      call symm(amat=overlap, bmat=tmpB, alpha=0.5_wp, &
         & cmat=tmpA, beta=1.0_wp)

      ! Add intermediate F += (S x X) x S
      call gemm(amat=tmpA, bmat=overlap, cmat=cache%prev_F(:, :, spin), &
         & beta = 1.0_wp)

      ! Add S x P and S x P x S diagonal onsite corrections
      do iao = 1, self%nao
         call axpy(xvec=overlap(:, iao), yvec=cache%prev_F(:, iao, spin), &
            & alpha=0.5_wp*gdiagSP(iao))
         cache%prev_F(iao, iao, spin) = cache%prev_F(iao, iao, spin) &
            & + 0.25_wp * gdiagSPS(iao)
      end do

      ! Save symmetrized Fock matrix for energy evaluation
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(self, cache, spin_factor, spin) &
      !$omp private(ii, jj, tmp)
      do ii = 1, self%nao
         cache%prev_F(ii, ii, spin) = -0.5_wp * spin_factor &
            & * cache%prev_F(ii, ii, spin)
         do jj = 1, ii-1
            tmp = -0.25_wp * spin_factor * &
               & (cache%prev_F(jj, ii, spin) + cache%prev_F(ii, jj, spin))
            cache%prev_F(jj, ii, spin) = tmp
            cache%prev_F(ii, jj, spin) = tmp
         end do
      end do
   end do

end subroutine get_KFock


!> Calculate exchange for one complex Hermitian k-point block.
!>
!> The real Gamma implementation historically assumes symmetric S and P.
!> Here every transpose implied by that algebra is promoted to a Hermitian
!> adjoint.  The returned matrix is the Hermitian projection of the same
!> linear Fock map.  The charge-dependent shell potential is evaluated as
!> 1/2 Tr[P dF/dq], which both preserves the Gamma functional and avoids
!> phase-dependent products of individual complex matrix elements.
subroutine get_KFock_kpoint(self, mol, cache, density, overlap, fock, vsh, energy)
   !> Instance of the Fock exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable exchange kernels and charge derivatives
   type(exchange_cache), intent(in) :: cache
   !> Complex Hermitian density matrix for all spin channels
   complex(wp), intent(in) :: density(:, :, :)
   !> Complex Hermitian overlap matrix
   complex(wp), intent(in) :: overlap(:, :)
   !> Complex Hermitian exchange Fock matrix
   complex(wp), intent(out) :: fock(:, :, :)
   !> Charge-dependent shell potential
   real(wp), intent(out) :: vsh(:)
   !> Optional exchange energy for this k-point block
   real(wp), intent(out), optional :: energy

   integer :: iat, izp, is, ish, iq, spin
   real(wp), allocatable :: dg_onsfx(:, :, :), dg_onsri(:, :)
   real(wp), allocatable :: zero_at(:, :), zero_sh(:, :)
   complex(wp), allocatable :: dfock(:, :, :)
   real(wp) :: tmp_energy

   call build_KFock_complex(self, mol, density, overlap, cache%g_mulliken, &
      & cache%g_bocorr, cache%g_onsfx, cache%g_onsri, fock)

   tmp_energy = 0.0_wp
   do spin = 1, size(density, 3)
      tmp_energy = tmp_energy + 0.5_wp*real(sum(conjg(density(:, :, spin)) &
         & * fock(:, :, spin)), wp)
   end do
   if (present(energy)) energy = tmp_energy

   allocate(zero_sh(self%nsh, self%nsh), zero_at(mol%nat, mol%nat), &
      & dg_onsfx(self%maxsh, self%maxsh, mol%nat), &
      & dg_onsri(self%maxsh, mol%nat), source=0.0_wp)
   allocate(dfock(self%nao, self%nao, size(density, 3)), &
      & source=(0.0_wp, 0.0_wp))
   vsh = 0.0_wp

   ! The onsite kernels are the only exchange kernels that depend explicitly
   ! on shell charges.  Apply their sparse derivative through the same linear
   ! Fock map, then contract with P.  This is intentionally an exact reference
   ! implementation; the contractions can be fused after validation in CP2K.
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = self%ish_at(iat)
      do ish = 1, self%nsh_id(izp)
         iq = is + ish
         dg_onsfx = 0.0_wp
         dg_onsri = 0.0_wp
         dg_onsfx(:, :, iat) = cache%dgdq_onsfx(:, :, iq)
         dg_onsri(:, iat) = cache%dgdq_onsri(:, iq)
         call build_KFock_complex(self, mol, density, overlap, zero_sh, &
            & zero_at, dg_onsfx, dg_onsri, dfock)
         do spin = 1, size(density, 3)
            vsh(iq) = vsh(iq) + 0.5_wp*real(sum( &
               & conjg(density(:, :, spin))*dfock(:, :, spin)), wp)
         end do
      end do
   end do

end subroutine get_KFock_kpoint


!> Calculate exchange for a complete regular Born-von Karman k-point mesh.
!>
!> Matrix products are local in k space, while every distance-dependent
!> Hadamard kernel is applied to image matrices.  This is the Fourier-dual
!> form of the Gamma-point exchange functional of the corresponding BvK
!> supercell and therefore couples all k points on the mesh.
subroutine get_KFock_kmesh(self, mol, cache, kernel, kfrac, weights, &
   & density, overlap, fock, vsh, energy)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(in) :: cache
   type(exchange_bvk_kernel), intent(in) :: kernel
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   complex(wp), intent(out) :: fock(:, :, :, :)
   real(wp), intent(out) :: vsh(:), energy

   integer :: spin, ik

   if (self%bvk_plan_matches(mol, cache, kernel, kfrac, weights)) then
      call build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
         & kernel, cache%g_onsfx, cache%g_onsri, fock, .true., &
         & cache%bvk_phase_forward, cache%bvk_phase_inverse)
   else
      call build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
         & kernel, cache%g_onsfx, cache%g_onsri, fock, .true.)
   end if

   energy = 0.0_wp
   do ik = 1, size(weights)
      do spin = 1, size(density, 3)
         energy = energy + 0.5_wp*weights(ik)*real(sum( &
            & conjg(density(:, :, spin, ik))*fock(:, :, spin, ik)), wp)
      end do
   end do

   call get_KFock_kmesh_vsh(self, mol, cache, weights, density, overlap, vsh)

end subroutine get_KFock_kmesh


!> Contract the complete-mesh energy with charge-dependent onsite kernels.
!>
!> The onsite kernels live only in the home image.  Reverse the whole-mesh
!> Fock assembly once to obtain dE/dG, then contract that compact gradient
!> with dG/dq.  This is algebraically identical to rebuilding the complete
!> Fock mesh for every shell charge, but avoids all repeated S*P products and
!> Fourier transforms.
subroutine get_KFock_kmesh_vsh(self, mol, cache, weights, density, overlap, vsh)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(in) :: cache
   real(wp), intent(in) :: weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   real(wp), intent(out) :: vsh(:)

   integer :: iao, iat, ik, iq, is, ish, izp, spin
   real(wp) :: spin_factor
   real(wp), allocatable :: grad_onsfx(:, :, :), grad_onsri(:, :)
   complex(wp), allocatable :: bq(:, :), bw(:, :), bt0(:, :), bu0(:, :), &
      & bv0(:, :), d0(:, :), dmat(:, :), p0(:, :), x0(:, :), xmat(:, :)
   complex(wp), allocatable :: bgp(:), bgsp(:), bgsps(:), diag2d(:), &
      & diagp(:), diagx(:)

   allocate(grad_onsfx(self%maxsh, self%maxsh, mol%nat), &
      & grad_onsri(self%maxsh, mol%nat), source=0.0_wp)
   allocate(bq(self%nao, self%nao), bw(self%nao, self%nao), &
      & bt0(self%nao, self%nao), bu0(self%nao, self%nao), &
      & bv0(self%nao, self%nao), d0(self%nao, self%nao), &
      & dmat(self%nao, self%nao), p0(self%nao, self%nao), &
      & x0(self%nao, self%nao), xmat(self%nao, self%nao), &
      & bgp(self%nao), bgsp(self%nao), bgsps(self%nao), &
      & diag2d(self%nao), diagp(self%nao), diagx(self%nao), &
      & source=(0.0_wp, 0.0_wp))

   spin_factor = 0.5_wp
   if (size(density, 3) > 1) spin_factor = 1.0_wp

   do spin = 1, size(density, 3)
      p0 = (0.0_wp, 0.0_wp)
      x0 = (0.0_wp, 0.0_wp)
      d0 = (0.0_wp, 0.0_wp)
      do ik = 1, size(weights)
         xmat = matmul(overlap(:, :, ik), density(:, :, spin, ik))
         dmat = 0.5_wp*matmul(xmat, overlap(:, :, ik))
         p0 = p0 + weights(ik)*density(:, :, spin, ik)
         x0 = x0 + weights(ik)*xmat
         d0 = d0 + weights(ik)*dmat
      end do

      bt0 = (0.0_wp, 0.0_wp)
      bu0 = (0.0_wp, 0.0_wp)
      bv0 = (0.0_wp, 0.0_wp)
      bgp = (0.0_wp, 0.0_wp)
      bgsp = (0.0_wp, 0.0_wp)
      bgsps = (0.0_wp, 0.0_wp)
      do ik = 1, size(weights)
         bq = -0.25_wp*spin_factor*weights(ik)*density(:, :, spin, ik)
         bu0 = bu0 + bq
         bw = matmul(bq, conjg(transpose(overlap(:, :, ik))))
         bt0 = bt0 + bw
         bv0 = bv0 + 0.5_wp*matmul( &
            & conjg(transpose(overlap(:, :, ik))), bw)
         do iao = 1, self%nao
            bgp(iao) = bgp(iao) + 0.25_wp*sum( &
               & conjg(overlap(:, iao, ik))*bw(:, iao))
            bgsp(iao) = bgsp(iao) + 0.5_wp*sum( &
               & conjg(overlap(:, iao, ik))*bq(:, iao))
            bgsps(iao) = bgsps(iao) + 0.25_wp*bq(iao, iao)
         end do
      end do

      call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, x0, &
         & 0.5_wp, grad_onsfx, adjoint_src=.true.)
      call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, d0, &
         & 0.5_wp, grad_onsfx)
      call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, p0, &
         & 0.5_wp, grad_onsfx)

      call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, x0, &
         & -2.0_wp, grad_onsri)
      call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, d0, &
         & -2.0_wp, grad_onsri)
      call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, p0, &
         & -2.0_wp, grad_onsri)

      do iao = 1, self%nao
         diagp(iao) = p0(iao, iao)
         diagx(iao) = x0(iao, iao)
         diag2d(iao) = 2.0_wp*d0(iao, iao)
      end do
      call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgp, &
         & diagp, grad_onsfx)
      call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsp, &
         & diagx, grad_onsfx)
      call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsps, &
         & diag2d, grad_onsfx)
   end do

   vsh = 0.0_wp
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = self%ish_at(iat)
      do ish = 1, self%nsh_id(izp)
         iq = is + ish
         vsh(iq) = sum(grad_onsfx(:, :, iat) &
            & *cache%dgdq_onsfx(:, :, iq)) + sum(grad_onsri(:, iat) &
            & *cache%dgdq_onsri(:, iq))
      end do
   end do

end subroutine get_KFock_kmesh_vsh


!> Build the whole-mesh exchange Fock matrix for fixed exchange kernels.
subroutine build_KFock_kmesh(self, mol, density, overlap, kfrac, weights, &
   & kernel, g_onsfx, g_onsri, fock, use_bvk_kernel, &
   & phase_forward_cached, phase_inverse_cached)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   type(exchange_bvk_kernel), intent(in) :: kernel
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   complex(wp), intent(out) :: fock(:, :, :, :)
   logical, intent(in) :: use_bvk_kernel
   complex(wp), contiguous, optional, target, intent(in) :: phase_forward_cached(:, :)
   complex(wp), contiguous, optional, target, intent(in) :: phase_inverse_cached(:, :)

   integer :: iao, icell, ii, jj, ik, spin
   logical :: use_cached_phases
   real(wp) :: angle, spin_factor
   complex(wp) :: phase, tmp
   complex(wp), allocatable, target :: amat(:, :, :), cmat(:, :, :)
   complex(wp), allocatable :: diagP(:), diagSP(:), diagSPS(:)
   complex(wp), allocatable :: gdiagP(:), gdiagSP(:), gdiagSPS(:)
   complex(wp), allocatable, target :: vmat(:, :, :)
   complex(wp), allocatable, target :: source_r(:, :, :), result_r(:, :, :)
   complex(wp), allocatable :: work(:, :)
   complex(wp), allocatable, target :: phase_forward_local(:, :), phase_inverse_local(:, :)
   complex(wp), contiguous, pointer :: phase_forward(:, :), phase_inverse(:, :)

   allocate(amat(self%nao, self%nao, size(weights)), &
      & cmat(self%nao, self%nao, size(weights)), &
      & vmat(self%nao, self%nao, size(weights)), &
      & source_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & result_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & work(self%nao, self%nao), diagP(self%nao), diagSP(self%nao), &
      & diagSPS(self%nao), gdiagP(self%nao), gdiagSP(self%nao), &
      & gdiagSPS(self%nao))
   use_cached_phases = .false.
   if (present(phase_forward_cached) .and. present(phase_inverse_cached)) then
      use_cached_phases = all(shape(phase_forward_cached) == &
         & [size(kernel%reps, 2), size(weights)]) .and. &
         & all(shape(phase_inverse_cached) == &
         & [size(weights), size(kernel%reps, 2)])
   end if
   if (use_cached_phases) then
      phase_forward => phase_forward_cached
      phase_inverse => phase_inverse_cached
   else
      allocate(phase_forward_local(size(kernel%reps, 2), size(weights)), &
         & phase_inverse_local(size(weights), size(kernel%reps, 2)))
      phase_forward => phase_forward_local
      phase_inverse => phase_inverse_local

      ! No ordering of the k points is assumed, and a common
      ! Monkhorst-Pack twist is retained explicitly in every phase.
      do ik = 1, size(weights)
         do icell = 1, size(kernel%reps, 2)
            angle = 2.0_wp*pi*dot_product(kfrac(:, ik), &
               & real(kernel%reps(:, icell), wp))
            phase = exp(cmplx(0.0_wp, angle, wp))
            phase_forward(icell, ik) = phase
            phase_inverse(ik, icell) = weights(ik)*conjg(phase)
         end do
      end do
   end if

   spin_factor = 0.5_wp
   if (size(density, 3) > 1) spin_factor = 1.0_wp
   fock = (0.0_wp, 0.0_wp)

   do spin = 1, size(density, 3)
      diagP = (0.0_wp, 0.0_wp)
      diagSP = (0.0_wp, 0.0_wp)
      diagSPS = (0.0_wp, 0.0_wp)
      do ik = 1, size(weights)
         amat(:, :, ik) = matmul(overlap(:, :, ik), density(:, :, spin, ik))
         cmat(:, :, ik) = 0.5_wp*matmul(amat(:, :, ik), overlap(:, :, ik))
         do iao = 1, self%nao
            diagP(iao) = diagP(iao) + weights(ik)*density(iao, iao, spin, ik)
            diagSP(iao) = diagSP(iao) + weights(ik)*amat(iao, iao, ik)
            diagSPS(iao) = diagSPS(iao) + 2.0_wp*weights(ik)*cmat(iao, iao, ik)
         end do
      end do
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagP, gdiagP)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagSP, gdiagSP)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagSPS, gdiagSPS)

      call apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, phase_forward, &
         & amat, source_r, result_r, &
         & g_onsfx, g_onsri, 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, &
         & use_bvk_kernel, .true.)
      call apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, phase_forward, &
         & cmat, source_r, result_r, &
         & g_onsfx, g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, &
         & use_bvk_kernel, .false.)

      ! A fixed-spin rank-four section is non-contiguous for UKS.  Copy it
      ! explicitly into the third reusable mesh buffer instead of relying on
      ! a compiler-created pack temporary at the contiguous kernel boundary.
      do ik = 1, size(weights)
         vmat(:, :, ik) = density(:, :, spin, ik)
      end do
      call apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, phase_forward, &
         & vmat, source_r, result_r, g_onsfx, g_onsri, 1.0_wp, 0.0_wp, &
         & 0.5_wp, -2.0_wp, use_bvk_kernel, .false.)

      do ik = 1, size(weights)
         work = amat(:, :, ik)
         do iao = 1, self%nao
            work(:, iao) = work(:, iao) &
               & + 0.25_wp*gdiagP(iao)*overlap(:, iao, ik)
         end do
         work = work + 0.5_wp*matmul(overlap(:, :, ik), vmat(:, :, ik))
         fock(:, :, spin, ik) = cmat(:, :, ik) &
            & + matmul(work, overlap(:, :, ik))
         do iao = 1, self%nao
            fock(:, iao, spin, ik) = fock(:, iao, spin, ik) &
               & + 0.5_wp*gdiagSP(iao)*overlap(:, iao, ik)
            fock(iao, iao, spin, ik) = fock(iao, iao, spin, ik) &
               & + 0.25_wp*gdiagSPS(iao)
         end do

         do ii = 1, self%nao
            fock(ii, ii, spin, ik) = cmplx(-0.5_wp*spin_factor &
               & *real(fock(ii, ii, spin, ik), wp), 0.0_wp, wp)
            do jj = 1, ii-1
               tmp = -0.25_wp*spin_factor*(fock(jj, ii, spin, ik) &
                  & + conjg(fock(ii, jj, spin, ik)))
               fock(jj, ii, spin, ik) = tmp
               fock(ii, jj, spin, ik) = conjg(tmp)
            end do
         end do
      end do
   end do

end subroutine build_KFock_kmesh


!> Apply one exchange kernel directly to already transformed BvK images.
subroutine apply_bvk_Kkernel_realspace(self, mol, kernel, source_r, result_r, &
   & g_onsfx, g_onsri, mulliken_scale, bocorr_scale, onsite_fx_scale, &
   & onsite_ri_scale, use_bvk_kernel, adjoint_src)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_bvk_kernel), intent(in) :: kernel
   complex(wp), intent(in) :: source_r(:, :, :)
   complex(wp), intent(out) :: result_r(:, :, :)
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   real(wp), intent(in) :: mulliken_scale, bocorr_scale
   real(wp), intent(in) :: onsite_fx_scale, onsite_ri_scale
   logical, intent(in) :: use_bvk_kernel, adjoint_src

   integer :: icell, origin

   result_r = (0.0_wp, 0.0_wp)
   origin = 0
   do icell = 1, size(kernel%reps, 2)
      if (all(kernel%reps(:, icell) == 0)) origin = icell
      if (use_bvk_kernel) then
         if (mulliken_scale /= 0.0_wp) &
            call shell_hadamard_add_complex(self%nsh, self%nao_sh, &
               & self%iao_sh, kernel%g_mulliken_r(:, :, icell), &
               & source_r(:, :, icell), mulliken_scale, &
               & result_r(:, :, icell))
         if (bocorr_scale /= 0.0_wp) &
            call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, &
               & kernel%g_bocorr_r(:, :, icell), source_r(:, :, icell), &
               & bocorr_scale, result_r(:, :, icell))
      end if
   end do

   if (origin > 0) then
      if (onsite_fx_scale /= 0.0_wp) &
         call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
            & source_r(:, :, origin), onsite_fx_scale, &
            & result_r(:, :, origin), adjoint_src=adjoint_src)
      if (onsite_ri_scale /= 0.0_wp) &
         call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
            & source_r(:, :, origin), onsite_ri_scale, &
            & result_r(:, :, origin))
   end if
end subroutine apply_bvk_Kkernel_realspace


!> Apply the forward exchange map to memory-reduced BvK image accumulators.
!>
!> On entry `amat_r`, `cmat_r`, and `vmat_r` contain the weighted inverse
!> transforms of S*P, 1/2*S*P*S, and P.  They are overwritten by their
!> kernel-applied image matrices.  No complete density, overlap, or Fock mesh
!> is materialized.  Energy and the charge-dependent shell potential are
!> contracted before the source images are overwritten.
subroutine get_KFock_stream_apply(self, mol, cache, kernel, phase_forward, &
   & weights, amat_r, cmat_r, vmat_r, gdiagP, gdiagSP, gdiagSPS, vsh, energy)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(in) :: cache
   type(exchange_bvk_kernel), intent(in) :: kernel
   complex(wp), intent(in) :: phase_forward(:, :)
   real(wp), intent(in) :: weights(:)
   complex(wp), intent(inout) :: amat_r(:, :, :, :), cmat_r(:, :, :, :), &
      & vmat_r(:, :, :, :)
   complex(wp), intent(out) :: gdiagP(:, :), gdiagSP(:, :), gdiagSPS(:, :)
   real(wp), intent(out) :: vsh(:), energy

   integer :: iao, iat, icell, iq, is, ish, izp, jcell, origin, spin
   real(wp) :: spin_factor
   complex(wp) :: coeff, energy_raw
   real(wp), allocatable :: grad_onsfx(:, :, :), grad_onsri(:, :)
   complex(wp), allocatable :: bgp(:), bgsp(:), bgsps(:), bmat(:, :), &
      & bt0(:, :), bu0(:, :), bv0(:, :), diag2d(:), diagp(:), diagx(:), &
      & result_r(:, :, :), result2_r(:, :, :)

   origin = 0
   do icell = 1, size(kernel%reps, 2)
      if (all(kernel%reps(:, icell) == 0)) origin = icell
   end do

   allocate(grad_onsfx(self%maxsh, self%maxsh, mol%nat), &
      & grad_onsri(self%maxsh, mol%nat), source=0.0_wp)
   allocate(bgp(self%nao), bgsp(self%nao), bgsps(self%nao), &
      & bmat(self%nao, self%nao), bt0(self%nao, self%nao), &
      & bu0(self%nao, self%nao), bv0(self%nao, self%nao), &
      & diag2d(self%nao), diagp(self%nao), diagx(self%nao), &
      & result_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & result2_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & source=(0.0_wp, 0.0_wp))

   energy = 0.0_wp
   vsh = 0.0_wp
   spin_factor = 0.5_wp
   if (size(amat_r, 4) > 1) spin_factor = 1.0_wp

   do spin = 1, size(amat_r, 4)
      do iao = 1, self%nao
         diagp(iao) = vmat_r(iao, iao, origin, spin)
         diagx(iao) = amat_r(iao, iao, origin, spin)
         diag2d(iao) = 2.0_wp*cmat_r(iao, iao, origin, spin)
      end do
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & diagp, gdiagP(:, spin))
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & diagx, gdiagSP(:, spin))
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & diag2d, gdiagSPS(:, spin))

      bt0 = -0.25_wp*spin_factor*conjg( &
         & transpose(amat_r(:, :, origin, spin)))
      bu0 = -0.25_wp*spin_factor*vmat_r(:, :, origin, spin)
      bv0 = -0.25_wp*spin_factor*cmat_r(:, :, origin, spin)
      do iao = 1, self%nao
         bgp(iao) = -0.125_wp*spin_factor &
            & *cmat_r(iao, iao, origin, spin)
         bgsp(iao) = -0.125_wp*spin_factor &
            & *amat_r(iao, iao, origin, spin)
         bgsps(iao) = -0.0625_wp*spin_factor &
            & *vmat_r(iao, iao, origin, spin)
      end do

      call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, &
         & amat_r(:, :, origin, spin), 0.5_wp, grad_onsfx, &
         & adjoint_src=.true.)
      call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, &
         & cmat_r(:, :, origin, spin), 0.5_wp, grad_onsfx)
      call onsite_fx_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, &
         & vmat_r(:, :, origin, spin), 0.5_wp, grad_onsfx)
      call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bt0, &
         & amat_r(:, :, origin, spin), -2.0_wp, grad_onsri)
      call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bu0, &
         & cmat_r(:, :, origin, spin), -2.0_wp, grad_onsri)
      call onsite_ri_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bv0, &
         & vmat_r(:, :, origin, spin), -2.0_wp, grad_onsri)
      call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgp, &
         & diagp, grad_onsfx)
      call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsp, &
         & diagx, grad_onsfx)
      call onsite_fx_symv_parameter_gradient_complex(mol%nat, mol%id, &
         & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, bgsps, &
         & diag2d, grad_onsfx)

      energy_raw = (0.0_wp, 0.0_wp)
      do iao = 1, self%nao
         energy_raw = energy_raw &
            & +0.5_wp*cmat_r(iao, iao, origin, spin)*gdiagP(iao, spin) &
            & +0.5_wp*conjg(amat_r(iao, iao, origin, spin)) &
            & *gdiagSP(iao, spin) &
            & +0.25_wp*vmat_r(iao, iao, origin, spin)*gdiagSPS(iao, spin)
      end do

      call apply_bvk_Kkernel_realspace(self, mol, kernel, &
         & amat_r(:, :, :, spin), result_r, cache%g_onsfx, cache%g_onsri, &
         & 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, .true., .true.)
      do icell = 1, size(kernel%reps, 2)
         bmat = (0.0_wp, 0.0_wp)
         do jcell = 1, size(kernel%reps, 2)
            coeff = sum(weights*conjg(phase_forward(icell, :) &
               & *phase_forward(jcell, :)))
            bmat = bmat + coeff*conjg( &
               & transpose(amat_r(:, :, jcell, spin)))
         end do
         energy_raw = energy_raw + sum(conjg(bmat)*result_r(:, :, icell))
      end do
      amat_r(:, :, :, spin) = result_r

      call apply_bvk_Kkernel_realspace(self, mol, kernel, &
         & cmat_r(:, :, :, spin), result_r, cache%g_onsfx, cache%g_onsri, &
         & 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .true., .false.)
      energy_raw = energy_raw + sum(conjg(vmat_r(:, :, :, spin))*result_r)

      call apply_bvk_Kkernel_realspace(self, mol, kernel, &
         & vmat_r(:, :, :, spin), result2_r, cache%g_onsfx, &
         & cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .true., .false.)
      energy_raw = energy_raw + sum(conjg(cmat_r(:, :, :, spin)) &
         & *result2_r)
      vmat_r(:, :, :, spin) = result2_r
      cmat_r(:, :, :, spin) = result_r
      energy = energy - 0.25_wp*spin_factor*real(energy_raw, wp)
   end do

   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = self%ish_at(iat)
      do ish = 1, self%nsh_id(izp)
         iq = is + ish
         vsh(iq) = sum(grad_onsfx(:, :, iat) &
            & *cache%dgdq_onsfx(:, :, iq)) + sum(grad_onsri(:, iat) &
            & *cache%dgdq_onsri(:, iq))
      end do
   end do
end subroutine get_KFock_stream_apply


!> Reconstruct one unweighted Fock block from kernel-applied BvK images.
subroutine get_KFock_stream_block(self, overlap, phase_forward, amat_r, cmat_r, &
   & vmat_r, gdiagP, gdiagSP, gdiagSPS, fock)
   class(exchange_fock), intent(in) :: self
   complex(wp), intent(in) :: overlap(:, :), phase_forward(:)
   complex(wp), intent(in) :: amat_r(:, :, :, :), cmat_r(:, :, :, :), &
      & vmat_r(:, :, :, :), gdiagP(:, :), gdiagSP(:, :), gdiagSPS(:, :)
   complex(wp), intent(out) :: fock(:, :, :)

   integer :: iao, icell, ii, jj, spin
   real(wp) :: spin_factor
   complex(wp) :: tmp
   complex(wp), allocatable :: amat(:, :), cmat(:, :), vmat(:, :), work(:, :)

   allocate(amat(self%nao, self%nao), cmat(self%nao, self%nao), &
      & vmat(self%nao, self%nao), work(self%nao, self%nao))
   spin_factor = 0.5_wp
   if (size(fock, 3) > 1) spin_factor = 1.0_wp
   do spin = 1, size(fock, 3)
      amat = (0.0_wp, 0.0_wp)
      cmat = (0.0_wp, 0.0_wp)
      vmat = (0.0_wp, 0.0_wp)
      do icell = 1, size(phase_forward)
         amat = amat + phase_forward(icell)*amat_r(:, :, icell, spin)
         cmat = cmat + phase_forward(icell)*cmat_r(:, :, icell, spin)
         vmat = vmat + phase_forward(icell)*vmat_r(:, :, icell, spin)
      end do
      work = amat
      do iao = 1, self%nao
         work(:, iao) = work(:, iao) &
            & + 0.25_wp*gdiagP(iao, spin)*overlap(:, iao)
      end do
      work = work + 0.5_wp*matmul(overlap, vmat)
      fock(:, :, spin) = cmat + matmul(work, overlap)
      do iao = 1, self%nao
         fock(:, iao, spin) = fock(:, iao, spin) &
            & + 0.5_wp*gdiagSP(iao, spin)*overlap(:, iao)
         fock(iao, iao, spin) = fock(iao, iao, spin) &
            & + 0.25_wp*gdiagSPS(iao, spin)
      end do
      do ii = 1, self%nao
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
end subroutine get_KFock_stream_block


!> Apply image-resolved and onsite Hadamard exchange kernels to a k-mesh.
subroutine apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, phase_forward, &
   & mesh_k, source_r, result_r, g_onsfx, g_onsri, mulliken_scale, &
   & bocorr_scale, onsite_fx_scale, onsite_ri_scale, use_bvk_kernel, &
   & adjoint_src)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_bvk_kernel), intent(in) :: kernel
   complex(wp), contiguous, intent(in) :: phase_inverse(:, :), phase_forward(:, :)
   complex(wp), contiguous, target, intent(inout) :: mesh_k(:, :, :)
   complex(wp), contiguous, target, intent(inout) :: source_r(:, :, :)
   complex(wp), contiguous, target, intent(inout) :: result_r(:, :, :)
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   real(wp), intent(in) :: mulliken_scale, bocorr_scale
   real(wp), intent(in) :: onsite_fx_scale, onsite_ri_scale
   logical, intent(in) :: use_bvk_kernel, adjoint_src

   integer :: icell, npair, origin
   complex(wp), pointer, contiguous :: mesh_k_2d(:, :), source_r_2d(:, :), &
      & result_r_2d(:, :)

   npair = self%nao*self%nao
   mesh_k_2d(1:npair, 1:size(mesh_k, 3)) => mesh_k
   source_r_2d(1:npair, 1:size(source_r, 3)) => source_r
   result_r_2d(1:npair, 1:size(result_r, 3)) => result_r

   ! Treat all AO pairs as independent rows of one dense transform.  This
   ! replaces the scalar cell-by-k loops and lets the compiler/BLAS backend
   ! execute the inverse transform as a complex matrix multiplication.
   call matmul_complex_contiguous(mesh_k_2d, phase_inverse, source_r_2d)
   result_r = (0.0_wp, 0.0_wp)

   origin = 0
   do icell = 1, size(kernel%reps, 2)
      if (all(kernel%reps(:, icell) == 0)) origin = icell
   end do

   if (use_bvk_kernel) then
      do icell = 1, size(kernel%reps, 2)
         if (mulliken_scale /= 0.0_wp) &
            call shell_hadamard_add_complex(self%nsh, self%nao_sh, &
               & self%iao_sh, kernel%g_mulliken_r(:, :, icell), &
               & source_r(:, :, icell), mulliken_scale, result_r(:, :, icell))
         if (bocorr_scale /= 0.0_wp) &
            call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, &
               & kernel%g_bocorr_r(:, :, icell), source_r(:, :, icell), &
               & bocorr_scale, result_r(:, :, icell))
      end do
   end if

   if (origin > 0) then
      if (onsite_fx_scale /= 0.0_wp) &
         call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
            & source_r(:, :, origin), onsite_fx_scale, result_r(:, :, origin), &
            & adjoint_src=adjoint_src)
      if (onsite_ri_scale /= 0.0_wp) &
         call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
            & source_r(:, :, origin), onsite_ri_scale, result_r(:, :, origin))
   end if

   ! The forward transform uses the same phase table and likewise batches all
   ! AO pairs in one complex matrix multiplication.
   call matmul_complex_contiguous(result_r_2d, phase_forward, mesh_k_2d)

end subroutine apply_bvk_Kkernel_inplace


!> Apply a BvK kernel out of place while retaining the source for reverse mode.
!>
!> The SCF path calls apply_bvk_Kkernel_inplace directly and reuses its two
!> real-space work arrays.  Geometry response needs both the original source
!> and the transformed result, so this wrapper copies into the caller-provided
!> result before invoking the same in-place implementation.
subroutine apply_bvk_Kkernel(self, mol, kernel, phase_inverse, phase_forward, &
   & source_k, g_onsfx, g_onsri, mulliken_scale, bocorr_scale, &
   & onsite_fx_scale, onsite_ri_scale, use_bvk_kernel, adjoint_src, result_k)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_bvk_kernel), intent(in) :: kernel
   complex(wp), intent(in) :: phase_inverse(:, :), phase_forward(:, :)
   complex(wp), contiguous, intent(in) :: source_k(:, :, :)
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   real(wp), intent(in) :: mulliken_scale, bocorr_scale
   real(wp), intent(in) :: onsite_fx_scale, onsite_ri_scale
   logical, intent(in) :: use_bvk_kernel, adjoint_src
   complex(wp), contiguous, target, intent(out) :: result_k(:, :, :)

   complex(wp), allocatable, target :: result_r(:, :, :), source_r(:, :, :)

   allocate(source_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & result_r(self%nao, self%nao, size(kernel%reps, 2)))
   result_k = source_k
   call apply_bvk_Kkernel_inplace(self, mol, kernel, phase_inverse, &
      & phase_forward, result_k, source_r, result_r, g_onsfx, g_onsri, &
      & mulliken_scale, bocorr_scale, onsite_fx_scale, onsite_ri_scale, &
      & use_bvk_kernel, adjoint_src)

end subroutine apply_bvk_Kkernel


!> Reverse an image-resolved and onsite Hadamard exchange map.
!>
!> Forward Fourier convention:
!> ``X_R = sum_k w_k exp(-i k.R) X_k`` and
!> ``Y_k = sum_R exp(+i k.R) Y_R``.  Consequently the result adjoint is
!> transformed without k weights, while the returned source adjoint carries
!> one weight for every k block.  Image-kernel responses remain oriented and
!> are deliberately not symmetrized between R and -R.
subroutine reverse_bvk_Kkernel(self, mol, kernel, weights, phase_inverse, &
   & phase_forward, source_k, result_adjoint_k, g_onsfx, g_onsri, &
   & mulliken_scale, bocorr_scale, onsite_fx_scale, onsite_ri_scale, &
   & use_bvk_kernel, adjoint_src, source_adjoint_k, mulliken_grad_r, &
   & bocorr_grad_r)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_bvk_kernel), intent(in) :: kernel
   real(wp), intent(in) :: weights(:)
   complex(wp), intent(in) :: phase_inverse(:, :), phase_forward(:, :)
   complex(wp), contiguous, target, intent(in) :: source_k(:, :, :)
   complex(wp), contiguous, target, intent(in) :: result_adjoint_k(:, :, :)
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   real(wp), intent(in) :: mulliken_scale, bocorr_scale
   real(wp), intent(in) :: onsite_fx_scale, onsite_ri_scale
   logical, intent(in) :: use_bvk_kernel, adjoint_src
   complex(wp), contiguous, target, intent(out) :: source_adjoint_k(:, :, :)
   real(wp), intent(inout) :: mulliken_grad_r(:, :, :), bocorr_grad_r(:, :, :)

   integer :: icell, ik, npair, origin
   complex(wp), allocatable, target :: result_adjoint_r(:, :, :), &
      & source_adjoint_r(:, :, :), source_r(:, :, :)
   complex(wp), pointer :: result_adjoint_k_2d(:, :), &
      & result_adjoint_r_2d(:, :), source_adjoint_k_2d(:, :), &
      & source_adjoint_r_2d(:, :), source_k_2d(:, :), source_r_2d(:, :)

   allocate(source_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & result_adjoint_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & source_adjoint_r(self%nao, self%nao, size(kernel%reps, 2)), &
      & source=(0.0_wp, 0.0_wp))
   npair = self%nao*self%nao
   source_k_2d(1:npair, 1:size(source_k, 3)) => source_k
   source_r_2d(1:npair, 1:size(source_r, 3)) => source_r
   result_adjoint_k_2d(1:npair, 1:size(result_adjoint_k, 3)) &
      & => result_adjoint_k
   result_adjoint_r_2d(1:npair, 1:size(result_adjoint_r, 3)) &
      & => result_adjoint_r
   source_adjoint_r_2d(1:npair, 1:size(source_adjoint_r, 3)) &
      & => source_adjoint_r
   source_adjoint_k_2d(1:npair, 1:size(source_adjoint_k, 3)) &
      & => source_adjoint_k

   source_r_2d = matmul(source_k_2d, phase_inverse)
   result_adjoint_r_2d = matmul(result_adjoint_k_2d, &
      & conjg(transpose(phase_forward)))

   origin = 0
   do icell = 1, size(kernel%reps, 2)
      if (all(kernel%reps(:, icell) == 0)) origin = icell
   end do

   if (use_bvk_kernel) then
      do icell = 1, size(kernel%reps, 2)
         if (mulliken_scale /= 0.0_wp) then
            call shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
               & self%iao_sh, result_adjoint_r(:, :, icell), &
               & source_r(:, :, icell), mulliken_scale, &
               & mulliken_grad_r(:, :, icell))
            call shell_hadamard_add_complex(self%nsh, self%nao_sh, &
               & self%iao_sh, kernel%g_mulliken_r(:, :, icell), &
               & result_adjoint_r(:, :, icell), mulliken_scale, &
               & source_adjoint_r(:, :, icell))
         end if
         if (bocorr_scale /= 0.0_wp) then
            call atom_parameter_gradient_complex(mol%nat, mol%id, &
               & self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, &
               & result_adjoint_r(:, :, icell), source_r(:, :, icell), &
               & bocorr_scale, bocorr_grad_r(:, :, icell))
            call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
               & self%nao_sh, self%ish_at, self%iao_sh, &
               & kernel%g_bocorr_r(:, :, icell), &
               & result_adjoint_r(:, :, icell), bocorr_scale, &
               & source_adjoint_r(:, :, icell))
         end if
      end do
   end if

   if (origin > 0) then
      if (onsite_fx_scale /= 0.0_wp) &
         call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
            & result_adjoint_r(:, :, origin), onsite_fx_scale, &
            & source_adjoint_r(:, :, origin), adjoint_src=adjoint_src)
      if (onsite_ri_scale /= 0.0_wp) &
         call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
            & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
            & result_adjoint_r(:, :, origin), onsite_ri_scale, &
            & source_adjoint_r(:, :, origin))
   end if

   source_adjoint_k_2d = matmul(source_adjoint_r_2d, phase_forward)
   do ik = 1, size(weights)
      source_adjoint_k(:, :, ik) = weights(ik)*source_adjoint_k(:, :, ik)
   end do

end subroutine reverse_bvk_Kkernel


!> Differentiate exchange on a complete regular BvK k-point mesh.
!>
!> ``overlap_grad`` is returned in the unweighted public convention
!> ``dE = sum_k w_k Re sum_ij conjg(overlap_grad_ij(k))*dS_ij(k)``.
!> The image-kernel adjoints already contain the complete primitive-cell BZ
!> contraction and must not be weighted again by the caller.
subroutine get_KGrad_kmesh(self, mol, cache, kernel, kfrac, weights, density, &
   & overlap, overlap_grad, mulliken_grad_r, bocorr_grad_r)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), target, intent(in) :: cache
   type(exchange_bvk_kernel), intent(in) :: kernel
   real(wp), intent(in) :: kfrac(:, :), weights(:)
   complex(wp), intent(in) :: density(:, :, :, :), overlap(:, :, :)
   complex(wp), intent(out) :: overlap_grad(:, :, :)
   real(wp), intent(out) :: mulliken_grad_r(:, :, :), bocorr_grad_r(:, :, :)

   integer :: iao, icell, ik, spin
   logical :: use_cached_phases
   real(wp) :: angle, spin_factor
   complex(wp) :: phase
   complex(wp), allocatable :: amat(:, :, :), &
      & bdiag(:), bdiag_p(:), bdiag_sp(:), bdiag_sps(:), &
      & btmp(:, :, :), dmat(:, :, :), ddensity(:, :, :), &
      & gdiag_p(:), gdiag_sp(:), gdiag_sps(:), tmat(:, :, :), &
      & umat(:, :, :), vmat(:, :, :), work(:, :, :), xmat(:, :, :)
   complex(wp), allocatable, target :: phase_forward_local(:, :), phase_inverse_local(:, :)
   complex(wp), contiguous, pointer :: phase_forward(:, :), phase_inverse(:, :)

   allocate(amat(self%nao, self%nao, size(weights)), &
      & btmp(self%nao, self%nao, size(weights)), &
      & dmat(self%nao, self%nao, size(weights)), &
      & ddensity(self%nao, self%nao, size(weights)), &
      & tmat(self%nao, self%nao, size(weights)), &
      & umat(self%nao, self%nao, size(weights)), &
      & vmat(self%nao, self%nao, size(weights)), &
      & work(self%nao, self%nao, size(weights)), &
      & xmat(self%nao, self%nao, size(weights)), &
      & bdiag(self%nao), bdiag_p(self%nao), bdiag_sp(self%nao), &
      & bdiag_sps(self%nao), gdiag_p(self%nao), gdiag_sp(self%nao), &
      & gdiag_sps(self%nao), source=(0.0_wp, 0.0_wp))
   use_cached_phases = self%bvk_plan_matches(mol, cache, kernel, kfrac, weights)
   if (use_cached_phases) then
      phase_forward => cache%bvk_phase_forward
      phase_inverse => cache%bvk_phase_inverse
   else
      allocate(phase_forward_local(size(kernel%reps, 2), size(weights)), &
         & phase_inverse_local(size(weights), size(kernel%reps, 2)))
      phase_forward => phase_forward_local
      phase_inverse => phase_inverse_local
      do ik = 1, size(weights)
         do icell = 1, size(kernel%reps, 2)
            angle = 2.0_wp*pi*dot_product(kfrac(:, ik), &
               & real(kernel%reps(:, icell), wp))
            phase = exp(cmplx(0.0_wp, angle, wp))
            phase_forward(icell, ik) = phase
            phase_inverse(ik, icell) = weights(ik)*conjg(phase)
         end do
      end do
   end if

   spin_factor = 0.5_wp
   if (size(density, 3) > 1) spin_factor = 1.0_wp
   overlap_grad = (0.0_wp, 0.0_wp)
   mulliken_grad_r = 0.0_wp
   bocorr_grad_r = 0.0_wp

   do spin = 1, size(density, 3)
      bdiag_p = (0.0_wp, 0.0_wp)
      bdiag_sp = (0.0_wp, 0.0_wp)
      bdiag_sps = (0.0_wp, 0.0_wp)
      gdiag_p = (0.0_wp, 0.0_wp)
      gdiag_sp = (0.0_wp, 0.0_wp)
      gdiag_sps = (0.0_wp, 0.0_wp)

      do ik = 1, size(weights)
         xmat(:, :, ik) = matmul(overlap(:, :, ik), &
            & density(:, :, spin, ik))
         dmat(:, :, ik) = 0.5_wp*matmul(xmat(:, :, ik), overlap(:, :, ik))
         do iao = 1, self%nao
            gdiag_p(iao) = gdiag_p(iao) &
               & + weights(ik)*density(iao, iao, spin, ik)
            gdiag_sp(iao) = gdiag_sp(iao) &
               & + weights(ik)*xmat(iao, iao, ik)
            gdiag_sps(iao) = gdiag_sps(iao) &
               & + 2.0_wp*weights(ik)*dmat(iao, iao, ik)
         end do
      end do
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & gdiag_p, bdiag_p)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & gdiag_sp, bdiag_sp)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & gdiag_sps, bdiag_sps)
      gdiag_p = bdiag_p
      gdiag_sp = bdiag_sp
      gdiag_sps = bdiag_sps
      bdiag_sp = (0.0_wp, 0.0_wp)
      bdiag_sps = (0.0_wp, 0.0_wp)

      call apply_bvk_Kkernel(self, mol, kernel, phase_inverse, phase_forward, &
         & xmat, cache%g_onsfx, cache%g_onsri, 1.0_wp, -4.0_wp, &
         & 0.5_wp, -2.0_wp, .true., .true., tmat)
      call apply_bvk_Kkernel(self, mol, kernel, phase_inverse, phase_forward, &
         & dmat, cache%g_onsfx, cache%g_onsri, 1.0_wp, 0.0_wp, &
         & 0.5_wp, -2.0_wp, .true., .false., umat)
      call apply_bvk_Kkernel(self, mol, kernel, phase_inverse, phase_forward, &
         & density(:, :, spin, :), cache%g_onsfx, cache%g_onsri, &
         & 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .true., .false., vmat)

      do ik = 1, size(weights)
         work(:, :, ik) = tmat(:, :, ik)
         do iao = 1, self%nao
            work(:, iao, ik) = work(:, iao, ik) &
               & + 0.25_wp*gdiag_p(iao)*overlap(:, iao, ik)
         end do
         work(:, :, ik) = work(:, :, ik) &
            & + 0.5_wp*matmul(overlap(:, :, ik), vmat(:, :, ik))

         ! Seed E = -(spin_factor/4) sum_k w_k Re<P_k,A_k>.
         amat(:, :, ik) = -0.25_wp*spin_factor*weights(ik) &
            & *density(:, :, spin, ik)
         umat(:, :, ik) = amat(:, :, ik)
         btmp(:, :, ik) = matmul(amat(:, :, ik), &
            & conjg(transpose(overlap(:, :, ik))))
         overlap_grad(:, :, ik) = overlap_grad(:, :, ik) &
            & + matmul(conjg(transpose(work(:, :, ik))), amat(:, :, ik))
         do iao = 1, self%nao
            overlap_grad(:, iao, ik) = overlap_grad(:, iao, ik) &
               & + 0.5_wp*conjg(gdiag_sp(iao))*amat(:, iao, ik)
            bdiag_sp(iao) = bdiag_sp(iao) &
               & + 0.5_wp*dot_product(overlap(:, iao, ik), &
               & amat(:, iao, ik))
            bdiag_sps(iao) = bdiag_sps(iao) &
               & + 0.25_wp*amat(iao, iao, ik)
         end do

         do iao = 1, self%nao
            overlap_grad(:, iao, ik) = overlap_grad(:, iao, ik) &
               & + 0.25_wp*conjg(gdiag_p(iao))*btmp(:, iao, ik)
         end do
         overlap_grad(:, :, ik) = overlap_grad(:, :, ik) &
            & + 0.5_wp*matmul(btmp(:, :, ik), &
            & conjg(transpose(vmat(:, :, ik))))
         vmat(:, :, ik) = 0.5_wp*matmul( &
            & conjg(transpose(overlap(:, :, ik))), btmp(:, :, ik))
      end do

      call reverse_bvk_Kkernel(self, mol, kernel, weights, phase_inverse, &
         & phase_forward, dmat, umat, cache%g_onsfx, cache%g_onsri, &
         & 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .true., .false., tmat, &
         & mulliken_grad_r, bocorr_grad_r)
      call reverse_bvk_Kkernel(self, mol, kernel, weights, phase_inverse, &
         & phase_forward, xmat, btmp, cache%g_onsfx, cache%g_onsri, &
         & 1.0_wp, -4.0_wp, 0.5_wp, -2.0_wp, .true., .true., umat, &
         & mulliken_grad_r, bocorr_grad_r)
      call reverse_bvk_Kkernel(self, mol, kernel, weights, phase_inverse, &
         & phase_forward, density(:, :, spin, :), vmat, cache%g_onsfx, &
         & cache%g_onsri, 1.0_wp, 0.0_wp, 0.5_wp, -2.0_wp, .true., &
         & .false., ddensity, mulliken_grad_r, bocorr_grad_r)

      ! Reverse the three BZ-averaged onsite diagonal maps.  diag(P) has no
      ! overlap response; its kernel response belongs to the shell potential.
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & bdiag_sp, bdiag)
      do ik = 1, size(weights)
         do iao = 1, self%nao
            umat(iao, iao, ik) = umat(iao, iao, ik) &
               & + weights(ik)*bdiag(iao)
         end do
      end do
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & bdiag_sps, bdiag)
      do ik = 1, size(weights)
         do iao = 1, self%nao
            umat(iao, :, ik) = umat(iao, :, ik) &
               & + weights(ik)*bdiag(iao)*overlap(iao, :, ik)
            overlap_grad(:, iao, ik) = overlap_grad(:, iao, ik) &
               & + weights(ik)*bdiag(iao)*conjg(xmat(iao, :, ik))
         end do

         ! D_k = X_k S_k/2 and X_k = S_k P_k.
         umat(:, :, ik) = umat(:, :, ik) &
            & + 0.5_wp*matmul(tmat(:, :, ik), &
            & conjg(transpose(overlap(:, :, ik))))
         overlap_grad(:, :, ik) = overlap_grad(:, :, ik) &
            & + 0.5_wp*matmul(conjg(transpose(xmat(:, :, ik))), &
            & tmat(:, :, ik))
         overlap_grad(:, :, ik) = overlap_grad(:, :, ik) &
            & + matmul(umat(:, :, ik), &
            & conjg(transpose(density(:, :, spin, ik))))
      end do
   end do

   do ik = 1, size(weights)
      overlap_grad(:, :, ik) = 0.5_wp*(overlap_grad(:, :, ik) &
         & + conjg(transpose(overlap_grad(:, :, ik))))/weights(ik)
   end do

end subroutine get_KGrad_kmesh


!> Differentiate the complex k-point exchange functional.
!>
!> The returned overlap response obeys
!> `dE = real(sum(conjg(overlap_grad)*dS))`.  The two real operator responses
!> are direct derivatives with respect to the independent elements of the
!> symmetric kernels.  They can be accumulated over a full k star and must
!> then be contracted with get_mulliken_derivs_direct and
!> get_bocorr_derivs_direct.  All inputs and outputs are unweighted.
subroutine get_KGrad_kpoint(self, mol, cache, density, overlap, &
   & mulliken_grad, bocorr_grad, overlap_grad)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(in) :: cache
   complex(wp), intent(in) :: density(:, :, :), overlap(:, :)
   real(wp), intent(out) :: mulliken_grad(:, :), bocorr_grad(:, :)
   complex(wp), intent(out) :: overlap_grad(:, :)

   integer :: iao, spin
   real(wp) :: spin_factor
   real(wp), allocatable :: mulliken_raw(:, :), bocorr_raw(:, :)
   complex(wp), allocatable :: amat(:, :), bamat(:, :), bdiag(:), &
      & bdiag_p(:), bdiag_sp(:), bdiag_sps(:), bmat(:, :), btmp(:, :), &
      & dmat(:, :), gdiag_p(:), gdiag_sp(:), gdiag_sps(:), tmat(:, :), &
      & t2mat(:, :), umat(:, :), vmat(:, :), xmat(:, :)

   allocate(amat(self%nao, self%nao), bamat(self%nao, self%nao), &
      & bmat(self%nao, self%nao), btmp(self%nao, self%nao), &
      & dmat(self%nao, self%nao), tmat(self%nao, self%nao), &
      & t2mat(self%nao, self%nao), umat(self%nao, self%nao), &
      & vmat(self%nao, self%nao), xmat(self%nao, self%nao), &
      & bdiag(self%nao), bdiag_p(self%nao), bdiag_sp(self%nao), &
      & bdiag_sps(self%nao), gdiag_p(self%nao), gdiag_sp(self%nao), &
      & gdiag_sps(self%nao), source=(0.0_wp, 0.0_wp))
   allocate(mulliken_raw(self%nao, self%nao), &
      & bocorr_raw(self%nao, self%nao), source=0.0_wp)

   spin_factor = 0.5_wp
   if (size(density, 3) > 1) spin_factor = 1.0_wp
   overlap_grad = (0.0_wp, 0.0_wp)

   do spin = 1, size(density, 3)
      ! Forward intermediates of build_KFock_complex, kept explicitly for the
      ! reverse sweep below.
      xmat = matmul(overlap, density(:, :, spin))
      dmat = 0.5_wp*matmul(xmat, overlap)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & [(density(iao, iao, spin), iao=1, self%nao)], gdiag_p)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & [(xmat(iao, iao), iao=1, self%nao)], gdiag_sp)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & [(2.0_wp*dmat(iao, iao), iao=1, self%nao)], gdiag_sps)

      tmat = (0.0_wp, 0.0_wp)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & cache%g_mulliken, xmat, 1.0_wp, tmat)
      call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_bocorr, xmat, &
         & -4.0_wp, tmat)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, xmat, &
         & 0.5_wp, tmat, adjoint_src=.true.)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, xmat, &
         & -2.0_wp, tmat)

      umat = (0.0_wp, 0.0_wp)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & cache%g_mulliken, dmat, 1.0_wp, umat)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, dmat, &
         & 0.5_wp, umat)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, dmat, &
         & -2.0_wp, umat)

      vmat = (0.0_wp, 0.0_wp)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & cache%g_mulliken, density(:, :, spin), 1.0_wp, vmat)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & density(:, :, spin), 0.5_wp, vmat)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, &
         & density(:, :, spin), -2.0_wp, vmat)

      t2mat = tmat + 0.5_wp*matmul(overlap, vmat)
      do iao = 1, self%nao
         t2mat(:, iao) = t2mat(:, iao) &
            & + 0.25_wp*gdiag_p(iao)*overlap(:, iao)
      end do
      amat = umat + matmul(t2mat, overlap)
      do iao = 1, self%nao
         amat(:, iao) = amat(:, iao) &
            & + 0.5_wp*gdiag_sp(iao)*overlap(:, iao)
         amat(iao, iao) = amat(iao, iao) + 0.25_wp*gdiag_sps(iao)
      end do

      ! Reverse sweep, seeded by E=-(spin_factor/4)*Re<P,A>.
      bamat = -0.25_wp*spin_factor*density(:, :, spin)
      bmat = bamat  ! adjoint of U
      btmp = matmul(bamat, conjg(transpose(overlap)))  ! adjoint of T2
      overlap_grad = overlap_grad + matmul(conjg(transpose(t2mat)), bamat)
      bdiag_sp = (0.0_wp, 0.0_wp)
      bdiag_sps = (0.0_wp, 0.0_wp)
      do iao = 1, self%nao
         overlap_grad(:, iao) = overlap_grad(:, iao) &
            & + 0.5_wp*conjg(gdiag_sp(iao))*bamat(:, iao)
         bdiag_sp(iao) = 0.5_wp*dot_product(overlap(:, iao), bamat(:, iao))
         bdiag_sps(iao) = 0.25_wp*bamat(iao, iao)
      end do

      ! T2 = T + 1/4 S diag(gdiagP) + 1/2 S V.
      tmat = btmp  ! adjoint of T
      bdiag_p = (0.0_wp, 0.0_wp)
      do iao = 1, self%nao
         overlap_grad(:, iao) = overlap_grad(:, iao) &
            & + 0.25_wp*conjg(gdiag_p(iao))*btmp(:, iao)
         bdiag_p(iao) = 0.25_wp*dot_product(overlap(:, iao), btmp(:, iao))
      end do
      overlap_grad = overlap_grad &
         & + 0.5_wp*matmul(btmp, conjg(transpose(vmat)))
      vmat = 0.5_wp*matmul(conjg(transpose(overlap)), btmp) ! adjoint of V

      ! U response: U = Gm(D) + 1/2 Gfx(D) - 2 Gri(D).
      dmat = (0.0_wp, 0.0_wp)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & cache%g_mulliken, bmat, 1.0_wp, dmat)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, bmat, &
         & 0.5_wp, dmat)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, bmat, &
         & -2.0_wp, dmat)
      call shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
         & self%iao_sh, bmat, 0.5_wp*matmul(xmat, overlap), 1.0_wp, &
         & mulliken_raw)

      ! T response, including the adjoint of Gfx(X^H).
      bamat = (0.0_wp, 0.0_wp) ! adjoint of X
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & cache%g_mulliken, tmat, 1.0_wp, bamat)
      call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_bocorr, tmat, &
         & -4.0_wp, bamat)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, tmat, &
         & 0.5_wp, bamat, adjoint_src=.true.)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsri, tmat, &
         & -2.0_wp, bamat)
      call shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
         & self%iao_sh, tmat, xmat, 1.0_wp, mulliken_raw)
      call atom_parameter_gradient_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, tmat, xmat, -4.0_wp, &
         & bocorr_raw)

      ! V response contributes only to the Mulliken geometry kernel.
      call shell_parameter_gradient_complex(self%nsh, self%nao_sh, &
         & self%iao_sh, vmat, density(:, :, spin), 1.0_wp, mulliken_raw)

      ! Reverse the three onsite diagonal maps.
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & bdiag_sp, bdiag)
      do iao = 1, self%nao
         bamat(iao, iao) = bamat(iao, iao) + bdiag(iao)
      end do
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, cache%g_onsfx, &
         & bdiag_sps, bdiag)
      do iao = 1, self%nao
         bamat(iao, :) = bamat(iao, :) + bdiag(iao)*overlap(iao, :)
         overlap_grad(:, iao) = overlap_grad(:, iao) &
            & + bdiag(iao)*conjg(xmat(iao, :))
      end do

      ! D = 1/2 X S and X = S P.
      bamat = bamat + 0.5_wp*matmul(dmat, conjg(transpose(overlap)))
      overlap_grad = overlap_grad &
         & + 0.5_wp*matmul(conjg(transpose(xmat)), dmat)
      overlap_grad = overlap_grad &
         & + matmul(bamat, conjg(transpose(density(:, :, spin))))
   end do

   ! Only Hermitian overlap variations are physical at a full k point.
   overlap_grad = 0.5_wp*(overlap_grad+conjg(transpose(overlap_grad)))
   call symmetrize_shell_parameter_gradient(self%nsh, self%nao_sh, &
      & self%iao_sh, mulliken_raw, mulliken_grad)
   call symmetrize_atom_parameter_gradient(mol%nat, mol%id, self%nsh_id, &
      & self%nao_sh, self%ish_at, self%iao_sh, bocorr_raw, bocorr_grad)

end subroutine get_KGrad_kpoint


!> Apply the g-xTB exchange Fock map to complex Hermitian matrices.
subroutine build_KFock_complex(self, mol, density, overlap, g_mulliken, &
   & g_bocorr, g_onsfx, g_onsri, fock)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   complex(wp), intent(in) :: density(:, :, :), overlap(:, :)
   real(wp), intent(in) :: g_mulliken(:, :), g_bocorr(:, :)
   real(wp), intent(in) :: g_onsfx(:, :, :), g_onsri(:, :)
   complex(wp), intent(out) :: fock(:, :, :)

   integer :: iao, ii, jj, spin
   real(wp) :: spin_factor
   complex(wp) :: tmp
   complex(wp), allocatable :: diagP(:), diagSP(:), gdiagP(:), &
      & gdiagSP(:), gdiagSPS(:), tmpA(:, :), tmpB(:, :), tmpSPS(:)

   allocate(tmpA(self%nao, self%nao), tmpB(self%nao, self%nao), &
      & diagP(self%nao), diagSP(self%nao), gdiagP(self%nao), &
      & gdiagSP(self%nao), tmpSPS(self%nao), gdiagSPS(self%nao), &
      & source=(0.0_wp, 0.0_wp))

   spin_factor = 0.5_wp
   if (size(density, 3) > 1) spin_factor = 1.0_wp
   fock = (0.0_wp, 0.0_wp)

   do spin = 1, size(density, 3)
      ! A = S P
      tmpA = matmul(overlap, density(:, :, spin))

      do iao = 1, self%nao
         diagP(iao) = density(iao, iao, spin)
         diagSP(iao) = tmpA(iao, iao)
         tmpSPS(iao) = sum(tmpA(iao, :)*conjg(overlap(iao, :)))
      end do
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagP, gdiagP)

      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, diagSP, gdiagSP)
      call onsite_fx_symv_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, tmpSPS, gdiagSPS)

      ! Initial 1/2 S P S intermediate.
      fock(:, :, spin) = 0.5_wp*matmul(tmpA, overlap)

      ! X = g*(S P), including the adjoint P S onsite contribution.
      tmpB = (0.0_wp, 0.0_wp)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & g_mulliken, tmpA, 1.0_wp, tmpB)
      call atom_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_bocorr, tmpA, &
         & -4.0_wp, tmpB)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, tmpA, &
         & 0.5_wp, tmpB, adjoint_src=.true.)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, tmpA, &
         & -2.0_wp, tmpB)
      tmpA = tmpB

      ! g*(1/2 S P S)
      tmpB = (0.0_wp, 0.0_wp)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
         & fock(:, :, spin), 0.5_wp, tmpB)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
         & fock(:, :, spin), -2.0_wp, tmpB)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & g_mulliken, fock(:, :, spin), 1.0_wp, tmpB)
      fock(:, :, spin) = tmpB

      ! g*P
      tmpB = (0.0_wp, 0.0_wp)
      call shell_hadamard_add_complex(self%nsh, self%nao_sh, self%iao_sh, &
         & g_mulliken, density(:, :, spin), 1.0_wp, tmpB)
      call onsite_fx_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsfx, &
         & density(:, :, spin), 0.5_wp, tmpB)
      call onsite_ri_hadamard_add_complex(mol%nat, mol%id, self%nsh_id, &
         & self%nao_sh, self%ish_at, self%iao_sh, g_onsri, &
         & density(:, :, spin), -2.0_wp, tmpB)

      do iao = 1, self%nao
         tmpA(:, iao) = tmpA(:, iao) &
            & + 0.25_wp*gdiagP(iao)*overlap(:, iao)
      end do
      tmpA = tmpA + 0.5_wp*matmul(overlap, tmpB)
      fock(:, :, spin) = fock(:, :, spin) + matmul(tmpA, overlap)

      do iao = 1, self%nao
         fock(:, iao, spin) = fock(:, iao, spin) &
            & + 0.5_wp*gdiagSP(iao)*overlap(:, iao)
         fock(iao, iao, spin) = fock(iao, iao, spin) &
            & + 0.25_wp*gdiagSPS(iao)
      end do

      ! Hermitian projection and exchange sign/spin scaling.
      do ii = 1, self%nao
         fock(ii, ii, spin) = cmplx(-0.5_wp*spin_factor &
            & * real(fock(ii, ii, spin), wp), 0.0_wp, wp)
         do jj = 1, ii-1
            tmp = -0.25_wp*spin_factor*(fock(jj, ii, spin) &
               & + conjg(fock(ii, jj, spin)))
            fock(jj, ii, spin) = tmp
            fock(ii, jj, spin) = conjg(tmp)
         end do
      end do
   end do

end subroutine build_KFock_complex


!> Calculate exchange contribution to the gradient
subroutine get_KGrad(self, mol, cache, density, overlap, mulliken_grad, &
   & bocorr_grad, ao_grad)
   !> Instance of the exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container with intermediates
   type(exchange_cache), intent(inout) :: cache
   !> Density matrix
   real(wp), intent(in) :: density(:, :, :)
   !> Overlap matrix
   real(wp), intent(in) :: overlap(:, :)
   !> Operator gradient w.r.t. the Mulliken gamma matrix
   real(wp), contiguous, intent(out) :: mulliken_grad(:, :)
   !> Operator gradient w.r.t. the bond-order correlation matrix
   real(wp), contiguous, intent(out) :: bocorr_grad(:, :)
   !> Orbital gradient contribution to the energy-weighted density matrix
   real(wp), contiguous, intent(inout) :: ao_grad(:, :)

   integer :: spin, iao
   real(wp), allocatable :: tmpA(:, :), tmpB(:, :), tmpC(:, :), tmpD(:, :)
   real(wp), allocatable :: tmpVec(:), diagP(:), diagSP(:), tmpSPSvec(:)
   real(wp) :: spin_factor

   mulliken_grad = 0.0_wp
   bocorr_grad = 0.0_wp

   allocate(tmpA(self%nao, self%nao), tmpB(self%nao, self%nao), &
      & tmpC(self%nao, self%nao), tmpD(self%nao, self%nao), &
      & tmpVec(self%nao), tmpSPSvec(self%nao), diagP(self%nao), &
      & diagSP(self%nao), source = 0.0_wp)

   ! Select spin factor to cancel the quadratic dependence of the exchange energy
   ! on the occupation numbers (0.5 for restricted, and 1.0 for unrestricted)
   spin_factor = 0.5_wp
   if(size(density, 3) .gt. 1) then
      spin_factor = 1.0_wp
   end if

   ! Evaluate the operator and overlap derivative matrices for Mulliken
   ! and onsite approximated Fock exchange for a symmetric density matrix
   do spin = 1, size(density, 3)
      
      ! Intermediate A = S x P
      call symm(amat=overlap, bmat=density(:, :, spin), cmat=tmpA)

      ! Collect S x P diagonal onsite correction
      tmpVec(:) = 0.0_wp
      do iao = 1, self%nao
         tmpVec(iao) = tmpA(iao, iao)
      end do
      call onsite_fx_symv(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_onsfx, tmpVec, diagSP)

      ! Collect P diagonal onsite correction
      tmpVec(:) = 0.0_wp
      do iao = 1, self%nao
         tmpVec(iao) = density(iao, iao, spin)
      end do
      call onsite_fx_symv(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_onsfx, tmpVec, diagP)

      ! Out-of-place transpose of A for onsite correction
      tmpC = transpose(tmpA)

      ! Operator derivatives (P x S) * (S x P)
      mulliken_grad = mulliken_grad + 0.5_wp * spin_factor * tmpA * tmpC
      bocorr_grad = bocorr_grad + 2.0_wp * spin_factor * tmpA * tmpC

      ! Apply Mulliken, onsite, bond-order correlation matrices
      ! as g * (S x P) and g * (P x S)
      tmpB = 0.0_wp
      call shell_hadamard_add(self%nsh, self%nao_sh, self%iao_sh, cache%g_mulliken, &
         & tmpA, 1.0_wp, tmpB)
      call atom_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_bocorr, tmpA, -4.0_wp, tmpB)
      call onsite_fx_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsfx, tmpC, 0.5_wp, tmpB)
      call onsite_ri_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsri, tmpA, -2.0_wp, tmpB)

      ! Add diagP * S x P and (P x S) * diagP onsite corrections
      tmpD = 0.0_wp
      do iao = 1, self%nao
         call axpy(xvec=density(:, iao, spin), yvec=tmpD(:, iao), &
            & alpha=0.5_wp * diagSP(iao))
         call axpy(xvec=tmpC(:, iao), yvec=tmpD(:, iao), &
            & alpha=0.5_wp * diagP(iao))
      end do

      ! Apply Mulliken and onsite matrices as g * P
      tmpC = 0.0_wp
      call shell_hadamard_add(self%nsh, self%nao_sh, self%iao_sh, cache%g_mulliken, &
         & density(:, :, spin), 1.0_wp, tmpC)
      call onsite_fx_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsfx, density(:, :, spin), 0.5_wp, tmpC)
      call onsite_ri_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsri, density(:, :, spin), -2.0_wp, tmpC)

      ! Add intermediate B += S x (g * P)
      call symm(amat=overlap, bmat=tmpC, cmat=tmpB, beta=1.0_wp)

      ! Add intermediate D += P x (S x X)
      call symm(amat=density(:, :, spin), bmat=tmpB, cmat=tmpD, beta=1.0_wp)

      ! Out-of-place transpose of D for gradient contribution
      tmpB = transpose(tmpD)

      ! Symmetrized overlap derivative contribution
      ! 0.25 arises since we need to treat the upper and lower triangular overlap
      ! as indpendent variables to get the correct energy weighed density matrix
      ao_grad = ao_grad + 0.25_wp * spin_factor * (tmpD + tmpB)

      ! Intermediate B = (S x P) x S
      call gemm(amat=tmpA, bmat=overlap, cmat=tmpB)

      ! Operator derivatives P * (S x P) x S
      mulliken_grad = mulliken_grad + 0.5_wp * spin_factor * density(:, :, spin) * tmpB

      ! Collect S x P x S diagonal vector
      tmpSPSvec(:) = 0.0_wp
      do iao = 1, self%nao
        tmpSPSvec(iao) = tmpB(iao, iao)
      end do
   end do

end subroutine get_KGrad


!> Evaluate Mulliken exchange gamma matrix
subroutine get_mulliken_Kmatrix(self, mol, cache)
   !> Instance of the exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container with intermediates and the final exchange matrix
   type(exchange_cache), intent(inout) :: cache

   if (any(mol%periodic)) then
      call get_gmulliken_3d(mol, self%nsh_id, self%ish_at, self%hubbard, &
         & self%ondiag_scale, self%offdiag_scale, self%hubbard_exp, &
         & self%hubbard_exp_r0, self%rad, self%gexp, self%frscale, self%omega, &
         & self%lrscale, cache%wsc, cache%g_mulliken)
   else
      call get_gmulliken_0d(mol, self%nsh_id, self%ish_at, self%hubbard, &
         & self%ondiag_scale, self%offdiag_scale, self%hubbard_exp, &
         & self%hubbard_exp_r0, self%rad, self%gexp, self%frscale, self%omega, &
         & self%lrscale, cache%g_mulliken)
   end if

end subroutine get_mulliken_Kmatrix


!> Evaluate range separated exchange matrix for finite systems
subroutine get_gmulliken_0d(mol, nsh_id, ish_at, hubbard, ondiag_scale, &
   & offdiag_scale, hubbard_exp, hubbard_exp_r0, rad, gexp, frscale, &
   & omega, lrscale, g_mulliken)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Hubbard parameter parameter for each shell
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Diagonal scaling of the Fock exchange
   real(wp), intent(in) :: ondiag_scale
   !> Off-diagonal scaling of the Fock exchange
   real(wp), intent(in) :: offdiag_scale(:, :, :, :)
   !> Exponent of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Radius prefactor of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp_r0
   !> Radius for hubbard scaling
   real(wp), intent(in) :: rad(:, :)
   !> Exponent of exchange kernel
   real(wp), intent(in) :: gexp
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Range separation parameter
   real(wp), intent(in) :: omega
   !> Long-range scaling factor
   real(wp), intent(in) :: lrscale
   !> Mulliken exchange matrix
   real(wp), intent(out) :: g_mulliken(:, :)

   integer :: iat, jat, izp, jzp, is, js, ish, jsh
   real(wp) :: vec(3), r1, r1g, gam, denom, rsh, scale, damp

   g_mulliken(:, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) shared(g_mulliken) &
   !$omp shared(mol, nsh_id, ish_at, hubbard, gexp, hubbard_exp, hubbard_exp_r0) &
   !$omp shared(rad, offdiag_scale, ondiag_scale, frscale, omega, lrscale) &
   !$omp private(iat, izp, is, ish, jat, jzp, js, jsh) &
   !$omp private(gam, vec, r1, r1g, denom, rsh, scale, damp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = ish_at(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         js = ish_at(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         r1g = r1**gexp
         damp = exp (-(hubbard_exp + hubbard_exp_r0 * rad(izp, jzp)) * r1)
         do ish = 1, nsh_id(izp)
            do jsh = 1, nsh_id(jzp)
               scale = offdiag_scale(jsh, ish, jzp, izp) / damp
               gam = hubbard(jsh, ish, jzp, izp) * scale
               denom = (r1g + gam**(-gexp))**(1.0_wp/gexp)
               rsh = (frscale+lrscale*erf(omega*r1)) / denom

               g_mulliken(js+jsh, is+ish) = rsh
               g_mulliken(is+ish, js+jsh) = rsh
            end do
         end do
      end do
      ! Onsite terms
      do ish = 1, nsh_id(izp)
         do jsh = 1, ish-1
            scale = offdiag_scale(jsh, ish, izp, izp)
            gam = hubbard(jsh, ish, izp, izp) * scale * frscale

            g_mulliken(is+jsh, is+ish) = gam
            g_mulliken(is+ish, is+jsh) = gam
         end do
         ! Diagonal elements
         gam = hubbard(ish, ish, izp, izp) * ondiag_scale * frscale

         g_mulliken(is+ish, is+ish) = gam
      end do
   end do 

end subroutine get_gmulliken_0d


!> Evaluate range separated exchange matrix for periodic systems
subroutine get_gmulliken_3d(mol, nsh_id, ish_at, hubbard, &
   & ondiag_scale, offdiag_scale, hubbard_exp, hubbard_exp_r0, rad, &
   & gexp, frscale, omega, lrscale, wsc, g_mulliken)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Hubbard parameter parameter for each shell
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Diagonal scaling of the Fock exchange
   real(wp), intent(in) :: ondiag_scale
   !> Off-diagonal scaling of the Fock exchange
   real(wp), intent(in) :: offdiag_scale(:, :, :, :)
   !> Exponent of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Radius prefactor of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp_r0
   !> Radius for hubbard scaling
   real(wp), intent(in) :: rad(:, :)
   !> Exponent of exchange kernel
   real(wp), intent(in) :: gexp
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Range separation parameter
   real(wp), intent(in) :: omega
   !> Long-range scaling factor
   real(wp), intent(in) :: lrscale
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Mulliken exchange matrix
   real(wp), intent(out) :: g_mulliken(:, :)

   integer :: iat, jat, izp, jzp, is, js, ish, jsh, img
   real(wp) :: vec(3), r1, rsh, gam, scale, wsw

   g_mulliken(:, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(g_mulliken, mol, nsh_id, ish_at, hubbard) &
   !$omp shared(gexp, hubbard_exp, hubbard_exp_r0, rad, offdiag_scale) &
   !$omp shared(ondiag_scale, frscale, omega, lrscale, wsc) &
   !$omp private(iat, izp, is, ish, jat, jzp, js, jsh) &
   !$omp private(vec, r1, rsh, gam, scale, wsw, img)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = ish_at(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         js = ish_at(jat)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            r1 = norm2(vec)
            if (r1 < eps) cycle
            do ish = 1, nsh_id(izp)
               do jsh = 1, nsh_id(jzp)
                  rsh = wsw * get_gmulliken_pair(r1, hubbard(jsh, ish, jzp, izp), &
                     & offdiag_scale(jsh, ish, jzp, izp), hubbard_exp, &
                     & hubbard_exp_r0, rad(izp, jzp), gexp, frscale, lrscale, omega)

                  !$omp atomic
                  g_mulliken(js+jsh, is+ish) = g_mulliken(js+jsh, is+ish) + rsh
                  !$omp atomic
                  g_mulliken(is+ish, js+jsh) = g_mulliken(is+ish, js+jsh) + rsh
               end do
            end do
         end do
      end do
      ! Onsite terms
      do ish = 1, nsh_id(izp)
         do jsh = 1, ish-1
            scale = offdiag_scale(jsh, ish, izp, izp)
            gam = hubbard(jsh, ish, izp, izp) * scale * frscale

            !$omp atomic
            g_mulliken(is+jsh, is+ish) = g_mulliken(is+jsh, is+ish) + gam
            !$omp atomic
            g_mulliken(is+ish, is+jsh) = g_mulliken(is+ish, is+jsh) + gam
         end do
         ! Diagonal elements
         gam = hubbard(ish, ish, izp, izp) * ondiag_scale * frscale

         !$omp atomic
         g_mulliken(is+ish, is+ish) = g_mulliken(is+ish, is+ish) + gam
      end do

      ! Self-interaction with periodic images
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         r1 = norm2(vec)
         if (r1 < eps) cycle

         do ish = 1, nsh_id(izp)
            ! Off-diagonal shell pairs
            do jsh = 1, ish-1
               rsh = wsw * get_gmulliken_pair(r1, hubbard(jsh, ish, izp, izp), &
                  & offdiag_scale(jsh, ish, izp, izp), hubbard_exp, &
                  & hubbard_exp_r0, rad(izp, izp), gexp, frscale, lrscale, omega)

               !$omp atomic
               g_mulliken(is+jsh, is+ish) = g_mulliken(is+jsh, is+ish) + rsh
               !$omp atomic
               g_mulliken(is+ish, is+jsh) = g_mulliken(is+ish, is+jsh) + rsh
            end do

            ! Same-shell block
            rsh = wsw * get_gmulliken_pair(r1, hubbard(ish, ish, izp, izp), &
               & offdiag_scale(ish, ish, izp, izp), hubbard_exp, &
               & hubbard_exp_r0, rad(izp, izp), gexp, frscale, lrscale, omega)

            !$omp atomic
            g_mulliken(is+ish, is+ish) = g_mulliken(is+ish, is+ish) + rsh
         end do
      end do
   end do

end subroutine get_gmulliken_3d


!> Evaluate the gradient of the Mulliken exchange energy
subroutine get_mulliken_derivs(self, mol, cache, mulliken_grad, gradient, sigma)
   !> Instance of the exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(exchange_cache), intent(inout) :: cache
   !> Operator gradient w.r.t. the gamma Mulliken matrix
   real(wp), contiguous, intent(in) :: mulliken_grad(:, :)
   !> Molecular gradient of the exchange energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   if (any(mol%periodic)) then
      call get_gmulliken_derivs_3d(mol, self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, &
         & self%hubbard, self%offdiag_scale, self%hubbard_exp, &
         & self%hubbard_exp_r0, self%rad, self%gexp, self%frscale, self%omega, &
         & self%lrscale, cache%wsc, mulliken_grad, gradient, sigma)
   else
      call get_gmulliken_derivs_0d(mol, self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, &
         & self%hubbard, self%ondiag_scale, self%offdiag_scale, self%hubbard_exp, &
         & self%hubbard_exp_r0, self%rad, self%gexp, self%frscale, self%omega, &
         & self%lrscale, mulliken_grad, gradient, sigma)
   end if

end subroutine get_mulliken_derivs


!> Evaluate Mulliken derivatives from direct symmetric-kernel responses.
!>
!> get_KGrad_kpoint returns one derivative for each independent element of the
!> symmetric shell kernel.  The established real-space derivative routines
!> predate that convention: their off-diagonal shell blocks carry the opposite
!> sign and their diagonal shell blocks additionally contain both occurrences
!> of the same kernel element.  Convert that representation here so callers do
!> not have to reproduce this implementation detail.
subroutine get_mulliken_derivs_direct(self, mol, cache, mulliken_grad, &
   & gradient, sigma)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(inout) :: cache
   real(wp), contiguous, intent(in) :: mulliken_grad(:, :)
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: ii, ish, ni
   real(wp), allocatable :: legacy_grad(:, :)

   allocate(legacy_grad(self%nao, self%nao), source=-mulliken_grad)
   do ish = 1, self%nsh
      ii = self%iao_sh(ish)
      ni = self%nao_sh(ish)
      legacy_grad(ii+1:ii+ni, ii+1:ii+ni) = &
         & -2.0_wp*mulliken_grad(ii+1:ii+ni, ii+1:ii+ni)
   end do
   call self%get_mulliken_derivs(mol, cache, legacy_grad, gradient, sigma)
end subroutine get_mulliken_derivs_direct


!> Evaluate derivatives of Mulliken exchange matrix for finite systems (0D)
subroutine get_gmulliken_derivs_0d(mol, nsh_id, nao_sh, ish_at, iao_sh, hubbard, &
   & ondiag_scale, offdiag_scale, hubbard_exp, hubbard_exp_r0, rad, &
   & gexp, frscale, omega, lrscale, mulliken_grad, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Hubbard parameter parameter for each shell
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Diagonal scaling of the Fock exchange (unused here)
   real(wp), intent(in) :: ondiag_scale
   !> Off-diagonal scaling of the Fock exchange
   real(wp), intent(in) :: offdiag_scale(:, :, :, :)
   !> Exponent of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Radius prefactor of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp_r0
   !> Radius for hubbard scaling
   real(wp), intent(in) :: rad(:, :)
   !> Exponent of exchange kernel
   real(wp), intent(in) :: gexp
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Range separation parameter
   real(wp), intent(in) :: omega
   !> Long-range scaling factor
   real(wp), intent(in) :: lrscale
   !> Operator gradient w.r.t. the Mulliken gamma matrix
   real(wp), intent(in) :: mulliken_grad(:, :)
   !> Molecular gradient of the exchange energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp, is, js, ii, jj, ish, jsh, iao, jao
   real(wp) :: vec(3), r1, r1g, scale, damp, exparg, rsh, drsh, gam
   real(wp) :: denom, denom_pow, denom_deriv, tmp, shell_grad, dG(3)

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)

   !$omp parallel default(none) &
   !$omp shared(mol, nsh_id, nao_sh, ish_at, iao_sh, hubbard) &
   !$omp shared(gexp, hubbard_exp, hubbard_exp_r0, rad, offdiag_scale) &
   !$omp shared(frscale, omega, lrscale, mulliken_grad, gradient, sigma) &
   !$omp private(iat, izp, is, ii, ish, iao, jat, jzp, js, jj, jsh, jao, vec) &
   !$omp private(r1, r1g, scale, damp, exparg, rsh, drsh, gam, denom, denom_pow) &
   !$omp private(denom_deriv, tmp, shell_grad, dG, gradient_local, sigma_local)
   allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source = 0.0_wp)
   allocate(sigma_local(size(sigma, 1), size(sigma, 2)), source = 0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = ish_at(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         js = ish_at(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         r1g = r1**gexp

         ! Radius dependent hardness damping factor
         exparg = (hubbard_exp + hubbard_exp_r0 * rad(izp, jzp))
         damp = exp(-exparg * r1)
         
         ! Range-separation factor and derivative
         rsh = (frscale+lrscale*erf(omega*r1)) 
         drsh = lrscale * 2.0_wp*omega / sqrtpi * exp(-(omega*r1)**2)

         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js+jsh)

               ! Diagonal/off-diagonal scaling factor for hubbard parameter
               scale = offdiag_scale(jsh, ish, jzp, izp) / damp

               ! Scaled hubbard parameter
               gam = hubbard(jsh, ish, jzp, izp) * scale

               ! Damped coulomb interaction denominator
               denom = r1g + gam**(-gexp)
               denom_pow = denom**(1.0_wp/gexp)
               denom_deriv = denom_pow * denom
               
               ! Derivative of range-separation factor
               tmp = drsh / denom_pow

               ! Derivative of damped coulomb interaction denominator
               tmp = tmp + rsh * ( -r1g/r1 + exparg * gam**(-gexp)) / denom_deriv

               ! Collect all operator gradient contributions per shell pair
               shell_grad = 0.0_wp
               do iao = 1, nao_sh(is + ish)
                  do jao = 1, nao_sh(js + jsh)
                     shell_grad = shell_grad + mulliken_grad(ii+iao, jj+jao)
                  end do
               end do
               
               ! Add operator contribution to the gradient
               dG = shell_grad * tmp * vec(:)/r1
               gradient_local(:, iat) = gradient_local(:, iat) + dG
               gradient_local(:, jat) = gradient_local(:, jat) - dG
               sigma_local(:, :) = sigma_local - 0.5_wp * (spread(vec, 1, 3) &
                  & * spread(dG, 2, 3) + spread(dG, 1, 3) * spread(vec, 2, 3))
            end do
         end do
      end do
   end do
   !$omp critical (get_gmulliken_derivs_0d_)
   gradient(:, :) = gradient + gradient_local
   sigma(:, :) = sigma + sigma_local
   !$omp end critical (get_gmulliken_derivs_0d_)
   deallocate(gradient_local, sigma_local)
   !$omp end parallel

end subroutine get_gmulliken_derivs_0d


!> Evaluate derivatives of Mulliken exchange matrix for periodic systems (3D)
subroutine get_gmulliken_derivs_3d(mol, nsh_id, nao_sh, ish_at, iao_sh, hubbard, &
   & offdiag_scale, hubbard_exp, hubbard_exp_r0, rad, &
   & gexp, frscale, omega, lrscale, wsc, mulliken_grad, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Hubbard parameter parameter for each shell
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Off-diagonal scaling of the Fock exchange
   real(wp), intent(in) :: offdiag_scale(:, :, :, :)
   !> Exponent of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Radius prefactor of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp_r0
   !> Radius for hubbard scaling
   real(wp), intent(in) :: rad(:, :)
   !> Exponent of exchange kernel
   real(wp), intent(in) :: gexp
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Range separation parameter
   real(wp), intent(in) :: omega
   !> Long-range scaling factor
   real(wp), intent(in) :: lrscale
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Operator gradient w.r.t. the Mulliken gamma matrix
   real(wp), intent(in) :: mulliken_grad(:, :)
   !> Molecular gradient of the exchange energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp, is, js, ii, jj, ish, jsh, iao, jao, img
   real(wp) :: vec(3), r1, drsh, shell_grad, dG(3), wsw

   ! Thread-private arrays for reduction
   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)

   !$omp parallel default(none) &
   !$omp shared(mol, nsh_id, nao_sh, ish_at, iao_sh, hubbard) &
   !$omp shared(gexp, hubbard_exp, hubbard_exp_r0, rad, offdiag_scale) &
   !$omp shared(frscale, omega, lrscale, wsc, mulliken_grad) &
   !$omp shared(gradient, sigma) &
   !$omp private(iat, izp, is, ii, ish, jat, jzp, js, jj, jsh, iao, jao, img) &
   !$omp private(vec, r1, drsh, shell_grad, dG, wsw) &
   !$omp private(gradient_local, sigma_local)
   allocate(gradient_local(3, mol%nat), source = 0.0_wp)
   allocate(sigma_local(3, 3), source = 0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = ish_at(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         js = ish_at(jat)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            r1 = norm2(vec)
            if (r1 < eps) cycle

            do ish = 1, nsh_id(izp)
               ii = iao_sh(is+ish)
               do jsh = 1, nsh_id(jzp)
                  jj = iao_sh(js+jsh)
                  call get_gmulliken_pair_deriv(r1, hubbard(jsh, ish, jzp, izp), &
                     & offdiag_scale(jsh, ish, jzp, izp), hubbard_exp, &
                     & hubbard_exp_r0, rad(izp, jzp), gexp, frscale, lrscale, omega, &
                     & drsh)

                  shell_grad = 0.0_wp
                  do iao = 1, nao_sh(is + ish)
                     do jao = 1, nao_sh(js + jsh)
                        shell_grad = shell_grad + mulliken_grad(ii+iao, jj+jao)
                     end do
                  end do
                  dG(:) = shell_grad * wsw * drsh * vec(:) / r1

                  gradient_local(:, iat) = gradient_local(:, iat) - dG
                  gradient_local(:, jat) = gradient_local(:, jat) + dG

                  sigma_local(:, :) = sigma_local - spread(dG, 1, 3) &
                     & * spread(vec, 2, 3)
               end do
            end do
         end do
      end do

      ! Self-interaction periodic images
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         r1 = norm2(vec)
         if (r1 < eps) cycle

         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            do jsh = 1, ish-1
               jj = iao_sh(is+jsh)
               call get_gmulliken_pair_deriv(r1, hubbard(jsh, ish, izp, izp), &
                  & offdiag_scale(jsh, ish, izp, izp), hubbard_exp, &
                  & hubbard_exp_r0, rad(izp, izp), gexp, frscale, lrscale, omega, &
                  & drsh)

               shell_grad = 0.0_wp
               do iao = 1, nao_sh(is + ish)
                  do jao = 1, nao_sh(is + jsh)
                     shell_grad = shell_grad + mulliken_grad(ii+iao, jj+jao)
                  end do
               end do

               dG(:) = shell_grad * wsw * drsh * vec(:) / r1

               ! Both symmetric shell blocks depend on the self-image kernel.
               sigma_local(:, :) = sigma_local - spread(dG, 1, 3) &
                  & * spread(vec, 2, 3)
            end do

            ! Same-shell block
            call get_gmulliken_pair_deriv(r1, hubbard(ish, ish, izp, izp), &
               & offdiag_scale(ish, ish, izp, izp), hubbard_exp, &
               & hubbard_exp_r0, rad(izp, izp), gexp, frscale, lrscale, omega, &
               & drsh)

            shell_grad = 0.0_wp
            do iao = 1, nao_sh(is + ish)
               do jao = 1, nao_sh(is + ish)
                  shell_grad = shell_grad + mulliken_grad(ii+iao, ii+jao)
               end do
            end do

            dG(:) = shell_grad * wsw * drsh * vec(:) / r1
            sigma_local(:, :) = sigma_local - 0.5_wp * spread(dG, 1, 3) &
               & * spread(vec, 2, 3)
         end do
      end do
   end do
   !$omp critical (get_gmulliken_derivs_3d_)
   gradient(:, :) = gradient + gradient_local
   sigma(:, :) = sigma + sigma_local
   !$omp end critical (get_gmulliken_derivs_3d_)
   deallocate(gradient_local, sigma_local)
   !$omp end parallel
end subroutine get_gmulliken_derivs_3d

!> Evaluate onsite exchange gamma matrix
subroutine get_onsite_Kmatrix(self, mol, wfn, cache)
   !> Instance of the exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Reusable data container with intermediates and the final exchange matrices
   type(exchange_cache), intent(inout) :: cache

   call get_gons(mol, self%nsh_id, self%ish_at, self%onecxints, self%frscale, &
      & self%kq, wfn%qsh(:, 1), cache%g_onsfx, cache%g_onsri, cache%dgdq_onsfx, &
      & cache%dgdq_onsri)

end subroutine get_onsite_Kmatrix


!> Evaluate onsite exchange and rotational invariance correction matrices
subroutine get_gons(mol, nsh_id, ish_at, onecxints, frscale, kq, qsh, &
   & g_onsfx, g_onsri, dgdq_onsfx, dgdq_onsri)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> One center exchange integrals
   real(wp), intent(in) :: onecxints(:, :, :)
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Charge-dependence of the effective Fock exchange 
   real(wp), intent(in) :: kq(:, :)
   !> Shell-resolved charges
   real(wp), intent(in) :: qsh(:)
   !> Onsite exchange matrix
   real(wp), intent(out) :: g_onsfx(:, :, :)
   !> Onsite rotational invariance correction matrix
   real(wp), intent(out) :: g_onsri(:, :)
   !> Charge-derivative of the onsite exchange matrix
   real(wp), intent(out) :: dgdq_onsfx(:, :, :)
   !> Charge-derivative of the onsite rotational invariance correction matrix
   real(wp), intent(out) :: dgdq_onsri(:, :)

   integer :: iat, izp, is, ish, jsh
   real(wp) :: gam, dgami, dgamj, denom

   g_onsfx(:, :, :) = 0.0_wp
   g_onsri(:, :) = 0.0_wp
   dgdq_onsfx(:, :, :) = 0.0_wp
   dgdq_onsri(:, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(g_onsfx, g_onsri, dgdq_onsfx, dgdq_onsri) &
   !$omp shared(onecxints, mol, frscale, kq, qsh, nsh_id, ish_at) &
   !$omp private(iat, izp, is, ish, jsh, gam, dgami, dgamj, denom)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         do jsh = 1, nsh_id(izp)

            ! Scaled onsite exchange integral
            gam = frscale * onecxints(jsh, ish, izp) * (1.0_wp - 0.5_wp * &
               & (kq(ish, izp) * qsh(is+ish) + kq(jsh, izp) * qsh(is+jsh)))
            dgami = -0.5_wp * frscale * kq(ish, izp) * onecxints(jsh, ish, izp)
            dgamj = -0.5_wp * frscale * kq(jsh, izp) * onecxints(jsh, ish, izp)

            ! Onsite correction exchange matrix
            g_onsfx(jsh, ish, iat) = gam
            if (ish == jsh) then
               dgdq_onsfx(jsh, ish, is+ish) = dgami + dgamj
            else
               dgdq_onsfx(jsh, ish, is+ish) = dgami
               dgdq_onsfx(jsh, ish, is+jsh) = dgamj
            end if

            ! Rotational invariance correction matrix
            if (ish == jsh .and. ish > 1) then
               denom = 1.0_wp / real(4*ish - 2, wp)
               g_onsri(ish, iat) = gam * denom
               dgdq_onsri(ish, is+ish) = (dgami + dgamj) * denom
            end if
         end do
      end do
   end do 

end subroutine get_gons


!> Evaluate bond-order correlation correction matrix
subroutine get_bocorr_Kmatrix(self, mol, cache)
   !> Instance of the exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(exchange_cache), intent(inout) :: cache

   if (any(mol%periodic)) then
      call get_gbocorr_3d(mol, self%corr_scale, self%corr_exp, self%corr_rad, &
         & self%rad, cache%wsc, cache%g_bocorr)
   else
      call get_gbocorr_0d(mol, self%corr_scale, self%corr_exp, self%corr_rad, &
         & self%rad, cache%g_bocorr)
   end if

end subroutine get_bocorr_Kmatrix


!> Evaluate bond-order correlation correction matrix for finite systems
subroutine get_gbocorr_0d(mol, corr_scale, corr_exp, corr_rad, &
   & rad, g_bocorr)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Bond-order correlation scaling factor for each atom pair
   real(wp), intent(in) :: corr_scale(:, :)
   !> Bond-order correlation damping exponent
   real(wp), intent(in) :: corr_exp
   !> Bond-order correlation radius for each atom pair
   real(wp), intent(in) :: corr_rad(:, :)
   !> Reference van-der-Waals radius
   real(wp), intent(in) :: rad(:, :)
   !> bond_order correlation correction matrix
   real(wp), intent(out) :: g_bocorr(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r1, arg, damp, corr

   g_bocorr(:, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) shared(g_bocorr) &
   !$omp shared(mol, corr_rad, corr_exp, corr_scale, rad) &
   !$omp private(iat, izp, jat, jzp, vec, r1, arg, damp, corr)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         arg = corr_exp * (r1 - corr_rad(izp, jzp)) / rad(izp, jzp)
         damp = 0.5_wp * (1.0_wp + erf(-arg))
         corr = corr_scale(izp, jzp) * damp

         g_bocorr(jat, iat) = corr
         g_bocorr(iat, jat) = corr
      end do
   end do

end subroutine get_gbocorr_0d


!> Evaluate bond-order correlation correction matrix for periodic systems (3D)
subroutine get_gbocorr_3d(mol, corr_scale, corr_exp, corr_rad, &
   & rad, wsc, g_bocorr)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Bond-order correlation scaling factor for each atom pair
   real(wp), intent(in) :: corr_scale(:, :)
   !> Bond-order correlation damping exponent
   real(wp), intent(in) :: corr_exp
   !> Bond-order correlation radius for each atom pair
   real(wp), intent(in) :: corr_rad(:, :)
   !> Reference van-der-Waals radius
   real(wp), intent(in) :: rad(:, :)
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Bond-order correlation correction matrix
   real(wp), intent(out) :: g_bocorr(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vec(3), r1, arg, damp, wsw, corr

   g_bocorr(:, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(g_bocorr, mol) &
   !$omp shared(corr_rad, corr_exp, corr_scale, rad, wsc) &
   !$omp private(iat, izp, jat, jzp) &
   !$omp private(vec, r1, arg, damp, wsw, corr, img)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) &
               & - wsc%trans(:, wsc%tridx(img, jat, iat))
            r1 = norm2(vec)
            if (r1 < eps) cycle
            arg = corr_exp * (r1 - corr_rad(izp, jzp)) / rad(izp, jzp)
            damp = 0.5_wp * (1.0_wp + erf(-arg))
            corr = corr_scale(izp, jzp) * damp * wsw

            !$omp atomic
            g_bocorr(jat, iat) = g_bocorr(jat, iat) + corr
            !$omp atomic
            g_bocorr(iat, jat) = g_bocorr(iat, jat) + corr
         end do
      end do

      ! Self-image periodic interactions
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         r1 = norm2(vec)
         if (r1 < eps) cycle
         arg = corr_exp * (r1 - corr_rad(izp, izp)) / rad(izp, izp)
         damp = 0.5_wp * (1.0_wp + erf(-arg))
         corr = corr_scale(izp, izp) * damp * wsw

         !$omp atomic
         g_bocorr(iat, iat) = g_bocorr(iat, iat) + corr
      end do
   end do

end subroutine get_gbocorr_3d


!> Evaluate the gradient of the bond-order correlation energy
subroutine get_bocorr_derivs(self, mol, cache, bocorr_grad, gradient, sigma)
   !> Instance of the exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(exchange_cache), intent(inout) :: cache
   !> Operator gradient w.r.t. the bond-order correlation matrix
   real(wp), contiguous, intent(in) :: bocorr_grad(:, :)
   !> Molecular gradient of the exchange energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   if (any(mol%periodic)) then
      call get_gbocorr_derivs_3d(mol, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, self%corr_scale, self%corr_exp, self%corr_rad, &
         & self%rad, cache%wsc, bocorr_grad, gradient, sigma)
   else
      call get_gbocorr_derivs_0d(mol, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, self%corr_scale, self%corr_exp, self%corr_rad, &
         & self%rad, bocorr_grad, gradient, sigma)
   end if

end subroutine get_bocorr_derivs


!> Evaluate bond-order derivatives from direct symmetric-kernel responses.
!>
!> Off-diagonal atom blocks already use the direct dE/dg convention.  A
!> diagonal atom block represents one independent kernel element but occurs
!> twice in the historical derivative contraction, and therefore needs the
!> factor of two applied below.
subroutine get_bocorr_derivs_direct(self, mol, cache, bocorr_grad, gradient, &
   & sigma)
   class(exchange_fock), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(exchange_cache), intent(inout) :: cache
   real(wp), contiguous, intent(in) :: bocorr_grad(:, :)
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: iat, ii, is, izp, ni
   real(wp), allocatable :: legacy_grad(:, :)

   allocate(legacy_grad(self%nao, self%nao), source=bocorr_grad)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = self%ish_at(iat)
      ii = self%iao_sh(is+1)
      ni = sum(self%nao_sh(is+1:is+self%nsh_id(izp)))
      legacy_grad(ii+1:ii+ni, ii+1:ii+ni) = &
         & 2.0_wp*bocorr_grad(ii+1:ii+ni, ii+1:ii+ni)
   end do
   call self%get_bocorr_derivs(mol, cache, legacy_grad, gradient, sigma)
end subroutine get_bocorr_derivs_direct


!> Evaluate derivatives of Mulliken exchange matrix for finite systems (0D)
subroutine get_gbocorr_derivs_0d(mol, nsh_id, nao_sh, ish_at, iao_sh, corr_scale, &
   & corr_exp, corr_rad, rad, bocorr_grad, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Bond-order correlation scaling factor for each atom pair
   real(wp), intent(in) :: corr_scale(:, :)
   !> Bond-order correlation damping exponent
   real(wp), intent(in) :: corr_exp
   !> Bond-order correlation radius for each atom pair
   real(wp), intent(in) :: corr_rad(:, :)
   !> Reference van-der-Waals radius
   real(wp), intent(in) :: rad(:, :)
   !> Operator gradient w.r.t. the bond-order correlation matrix
   real(wp), intent(in) :: bocorr_grad(:, :)
   !> Molecular gradient of the exchange energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), intent(inout) :: sigma(:, :)
   
   integer :: iat, jat, izp, jzp, is, js, ii, jj, ish, jsh, iao, jao
   real(wp) :: vec(3), r1, arg, ddamp, dcorr, atom_grad, dG(3)

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)

   !$omp parallel default(none) &
   !$omp shared(mol, nsh_id, nao_sh, ish_at, iao_sh, corr_rad, corr_exp) &
   !$omp shared(corr_scale, rad, bocorr_grad, gradient, sigma) &
   !$omp private(iat, izp, is, ii, ish, iao, jat, jzp, js, jj, jsh, jao, vec, r1) &
   !$omp private(arg, ddamp, dcorr, atom_grad, dG, gradient_local, sigma_local)
   allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source = 0.0_wp)
   allocate(sigma_local(size(sigma, 1), size(sigma, 2)), source = 0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = ish_at(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         js = ish_at(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         arg = corr_exp * (r1 - corr_rad(izp, jzp)) / rad(izp, jzp)
         ddamp = corr_exp * exp(-arg**2) / (sqrtpi * rad(izp, jzp))
         dcorr = corr_scale(izp, jzp) * ddamp

         ! Collect all operator gradient contributions per atom pair
         atom_grad = 0.0_wp
         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js+jsh)
               do iao = 1, nao_sh(is + ish)
                  do jao = 1, nao_sh(js + jsh)
                     atom_grad = atom_grad + bocorr_grad(ii+iao, jj+jao)
                  end do
               end do
            end do
         end do
         ! Add operator contribution to the gradient
         dG(:) = atom_grad * dcorr * vec/r1
         gradient_local(:, iat) = gradient_local(:, iat) + dG
         gradient_local(:, jat) = gradient_local(:, jat) - dG
         sigma_local(:, :) = sigma_local - 0.5_wp * (spread(vec, 1, 3) &
            & * spread(dG, 2, 3) + spread(dG, 1, 3) * spread(vec, 2, 3))
      end do
   end do
   !$omp critical (get_gbocorr_derivs_0d_)
   gradient(:, :) = gradient + gradient_local
   sigma(:, :) = sigma + sigma_local
   !$omp end critical (get_gbocorr_derivs_0d_)
   deallocate(gradient_local, sigma_local)
   !$omp end parallel

end subroutine get_gbocorr_derivs_0d


!> Evaluate derivatives of bond-order correlation matrix for periodic systems (3D)
subroutine get_gbocorr_derivs_3d(mol, nsh_id, nao_sh, ish_at, iao_sh, corr_scale, &
   & corr_exp, corr_rad, rad, wsc, bocorr_grad, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Bond-order correlation scaling factor for each atom pair
   real(wp), intent(in) :: corr_scale(:, :)
   !> Bond-order correlation damping exponent
   real(wp), intent(in) :: corr_exp
   !> Bond-order correlation radius for each atom pair
   real(wp), intent(in) :: corr_rad(:, :)
   !> Reference van-der-Waals radius
   real(wp), intent(in) :: rad(:, :)
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Operator gradient w.r.t. the bond-order correlation matrix
   real(wp), intent(in) :: bocorr_grad(:, :)
   !> Molecular gradient of the exchange energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp, is, js, ii, jj, ish, jsh, iao, jao, img
   real(wp) :: vec(3), r1, arg, ddamp, dcorr, atom_grad, dG(3), wsw

   ! Thread-private arrays for reduction
   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)

   !$omp parallel default(none) &
   !$omp shared(mol, nsh_id, nao_sh, ish_at, iao_sh) &
   !$omp shared(corr_rad, corr_exp) &
   !$omp shared(corr_scale, rad, wsc, bocorr_grad, gradient, sigma) &
   !$omp private(iat, izp, is, ii, jat, jzp, js, jj, img, ish, jsh, iao, jao) &
   !$omp private(vec, r1, arg, ddamp, dcorr, atom_grad, dG, wsw) &
   !$omp private(gradient_local, sigma_local)
   allocate(gradient_local(3, mol%nat), source = 0.0_wp)
   allocate(sigma_local(3, 3), source = 0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = ish_at(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         js = ish_at(jat)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            r1 = norm2(vec)
            if (r1 < eps) cycle

            arg = corr_exp * (r1 - corr_rad(izp, jzp)) / rad(izp, jzp)
            ddamp = -corr_exp * exp(-arg**2) / (sqrtpi * rad(izp, jzp))
            dcorr = corr_scale(izp, jzp) * ddamp

            atom_grad = 0.0_wp
            do ish = 1, nsh_id(izp)
               ii = iao_sh(is+ish)
               do jsh = 1, nsh_id(jzp)
                  jj = iao_sh(js+jsh)
                  do iao = 1, nao_sh(is + ish)
                     do jao = 1, nao_sh(js + jsh)
                        atom_grad = atom_grad + bocorr_grad(ii+iao, jj+jao)
                     end do
                  end do
               end do
            end do

            dG(:) = atom_grad * wsw * dcorr * vec(:)/r1

            gradient_local(:, iat) = gradient_local(:, iat) + dG
            gradient_local(:, jat) = gradient_local(:, jat) - dG

            ! Here vec = r_i - r_j - T and gradient_i = dG.
            sigma_local(:, :) = sigma_local + spread(dG, 1, 3) &
               & * spread(vec, 2, 3)
         end do
      end do

      ! Self-interaction periodic images
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         r1 = norm2(vec)
         if (r1 < eps) cycle

         arg = corr_exp * (r1 - corr_rad(izp, izp)) / rad(izp, izp)
         ddamp = -corr_exp * exp(-arg**2) / (sqrtpi * rad(izp, izp))
         dcorr = corr_scale(izp, izp) * ddamp

         atom_grad = 0.0_wp
         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            do jsh = 1, nsh_id(izp)
               jj = iao_sh(is+jsh)
               do iao = 1, nao_sh(is + ish)
                  do jao = 1, nao_sh(is + jsh)
                     atom_grad = atom_grad + bocorr_grad(ii+iao, jj+jao)
                  end do
               end do
            end do
         end do

         dG(:) = atom_grad * wsw * dcorr * vec(:)/r1

         ! The diagonal self-image contribution is counted only once.
         sigma_local(:, :) = sigma_local + 0.5_wp * spread(dG, 1, 3) &
            & * spread(vec, 2, 3)
      end do
   end do
   !$omp critical (get_gbocorr_derivs_3d_)
   gradient(:, :) = gradient + gradient_local
   sigma(:, :) = sigma + sigma_local
   !$omp end critical (get_gbocorr_derivs_3d_)
   deallocate(gradient_local, sigma_local)
   !$omp end parallel
end subroutine get_gbocorr_derivs_3d


subroutine shell_hadamard_add(nsh, nao_sh, iao_sh, g_sh, src, alpha, dst, trans_src)
   !> Number of shells in the molecule
   integer, intent(in) :: nsh
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Shell-resolved exchange matrix: [nsh, nsh]
   real(wp), intent(in) :: g_sh(:, :)
   !> Source AO matrix
   real(wp), intent(in) :: src(:, :)
   !> Prefactor
   real(wp), intent(in) :: alpha
   !> Destination AO matrix
   real(wp), intent(inout) :: dst(:, :)
   !> Whether src should be read as transposed
   logical, intent(in), optional :: trans_src

   integer :: ish, jsh, ii, jj, ni, nj, iao, jao
   real(wp) :: scale
   logical :: trans

   trans = .false.
   if (present(trans_src)) trans = trans_src

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(nsh, nao_sh, iao_sh, g_sh, src, alpha, dst, trans) &
   !$omp private(ish, jsh, ii, jj, ni, nj, iao, jao, scale)
   do ish = 1, nsh
      ii = iao_sh(ish)
      ni = nao_sh(ish)

      do jsh = 1, nsh
         jj = iao_sh(jsh)
         nj = nao_sh(jsh)

         scale = alpha * g_sh(jsh, ish)
         if (abs(scale) < epsilon(1.0_wp)) cycle

         if (.not.trans) then
            dst(jj+1:jj+nj, ii+1:ii+ni) = dst(jj+1:jj+nj, ii+1:ii+ni) + &
               scale * src(jj+1:jj+nj, ii+1:ii+ni)
         else
            do iao = 1, ni
               do jao = 1, nj
                  dst(jj+jao, ii+iao) = dst(jj+jao, ii+iao) + &
                     scale * src(ii+iao, jj+jao)
               end do
            end do
         end if
      end do
   end do

end subroutine shell_hadamard_add


subroutine atom_hadamard_add(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
   & g_at, src, alpha, dst, trans_src)
   !> Number of atoms in the molecule
   integer, intent(in) :: nat
   !> Species identifier for each atom
   integer, intent(in) :: id(:)
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Atom-resolved matrix: [nat, nat]
   real(wp), intent(in) :: g_at(:, :)
   !> Source AO matrix
   real(wp), intent(in) :: src(:, :)
   !> Prefactor
   real(wp), intent(in) :: alpha
   !> Destination AO matrix
   real(wp), intent(inout) :: dst(:, :)
   !> Whether src should be read as transposed
   logical, intent(in), optional :: trans_src

   integer :: iat, jat, izp, jzp, is, js, ish, jsh, ii, jj, ni, nj, iao, jao
   real(wp) :: scale
   logical :: trans

   trans = .false.
   if (present(trans_src)) trans = trans_src

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(nat, id, nsh_id, nao_sh, ish_at, iao_sh, g_at, src, alpha, dst, trans) &
   !$omp private(iat, jat, izp, jzp, is, js, ish, jsh, ii, jj, ni, nj, iao, jao, scale)
   do iat = 1, nat
      izp = id(iat)
      is  = ish_at(iat)

      do jat = 1, nat
         jzp = id(jat)
         js  = ish_at(jat)

         scale = alpha * g_at(jat, iat)
         if (abs(scale) < epsilon(1.0_wp)) cycle

         do ish = 1, nsh_id(izp)
            ii = iao_sh(is + ish)
            ni = nao_sh(is + ish)

            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js + jsh)
               nj = nao_sh(js + jsh)

               if (.not.trans) then
                  dst(jj+1:jj+nj, ii+1:ii+ni) = dst(jj+1:jj+nj, ii+1:ii+ni) + &
                     & scale * src(jj+1:jj+nj, ii+1:ii+ni)
               else
                  do iao = 1, ni
                     do jao = 1, nj
                        dst(jj+jao, ii+iao) = dst(jj+jao, ii+iao) + &
                           & scale * src(ii+iao, jj+jao)
                     end do
                  end do
               end if
            end do
         end do
      end do
   end do

end subroutine atom_hadamard_add


subroutine onsite_fx_hadamard_add(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
   & g_onsfx, src, alpha, dst, trans_src)
   !> Number of atoms in the molecule
   integer, intent(in) :: nat
   !> Species identifier for each atom
   integer, intent(in) :: id(:)
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Compact shell-block onsite exchange matrix: [maxsh, maxsh, nat]
   real(wp), intent(in) :: g_onsfx(:, :, :)
   !> Source AO matrix
   real(wp), intent(in) :: src(:, :)
   !> Prefactor
   real(wp), intent(in) :: alpha
   !> Destination AO matrix
   real(wp), intent(inout) :: dst(:, :)
   !> Whether src should be read as transposed
   logical, intent(in), optional :: trans_src

   integer :: iat, izp, is, ish, jsh, ii, jj, ni, nj, iao, jao
   real(wp) :: scale
   logical :: trans

   trans = .false.
   if (present(trans_src)) trans = trans_src

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(nat, id, ish_at, nsh_id, iao_sh, nao_sh, alpha, g_onsfx, src, dst, trans) &
   !$omp private(iat, izp, is, ish, jsh, ii, jj, ni, nj, iao, jao, scale)
   do iat = 1, nat
      izp = id(iat)
      is  = ish_at(iat)

      do ish = 1, nsh_id(izp)
         ii = iao_sh(is + ish)
         ni = nao_sh(is + ish)

         do jsh = 1, nsh_id(izp)
            jj = iao_sh(is + jsh)
            nj = nao_sh(is + jsh)

            scale = alpha * g_onsfx(jsh, ish, iat)
            if (abs(scale) < epsilon(1.0_wp)) cycle

            if (.not.trans) then
               dst(jj+1:jj+nj, ii+1:ii+ni) = dst(jj+1:jj+nj, ii+1:ii+ni) + &
                  & scale * src(jj+1:jj+nj, ii+1:ii+ni)
            else
               do iao = 1, ni
                  do jao = 1, nj
                     dst(jj+jao, ii+iao) = dst(jj+jao, ii+iao) + &
                        & scale * src(ii+iao, jj+jao)
                  end do
               end do
            end if
         end do
      end do
   end do

end subroutine onsite_fx_hadamard_add


subroutine onsite_fx_symv(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
   & g_onsfx, xvec, yvec)
   !> Number of atoms in the molecule
   integer, intent(in) :: nat
   !> Species identifier for each atom
   integer, intent(in) :: id(:)
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Compact diagonal storage of the onsite exchange matrix: [maxsh, maxsh, nat]
   real(wp), intent(in) :: g_onsfx(:, :, :)
   !> Input vector to be contracted: [nao]
   real(wp), intent(in) :: xvec(:)
   !> Output vector: [nao]
   real(wp), intent(out) :: yvec(:)

   integer :: iat, izp, is, ish, jsh, ii, jj, ni, nj
   real(wp) :: shsum

   yvec = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(nat, id, ish_at, nsh_id, iao_sh, nao_sh, g_onsfx, xvec, yvec) &
   !$omp private(iat, izp, is, ish, jsh, ii, jj, ni, nj, shsum)
   do iat = 1, nat
      izp = id(iat)
      is  = ish_at(iat)

      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)

         shsum = sum(xvec(ii+1:ii+ni))
         if (abs(shsum) < epsilon(1.0_wp)) cycle

         do jsh = 1, nsh_id(izp)
            jj = iao_sh(is+jsh)
            nj = nao_sh(is+jsh)

            yvec(jj+1:jj+nj) = yvec(jj+1:jj+nj) + g_onsfx(jsh, ish, iat) * shsum
         end do
      end do
   end do

end subroutine onsite_fx_symv


subroutine onsite_ri_hadamard_add(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
   & g_onsri, src, alpha, dst)
   !> Number of atoms in the molecule
   integer, intent(in) :: nat
   !> Species identifier for each atom
   integer, intent(in) :: id(:)
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Compact shell-diagonal onsite rotational invariance matrix: [maxsh, nat]
   real(wp), intent(in) :: g_onsri(:, :)
   !> Source AO matrix
   real(wp), intent(in) :: src(:, :)
   !> Prefactor
   real(wp), intent(in) :: alpha
   !> Destination AO matrix
   real(wp), intent(inout) :: dst(:, :)

   integer :: iat, izp, is, ish, ii, ni
   real(wp) :: scale

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(nat, id, ish_at, nsh_id, iao_sh, nao_sh, alpha, g_onsri, src, dst) &
   !$omp private(iat, izp, is, ish, ii, ni, scale)
   do iat = 1, nat
      izp = id(iat)
      is  = ish_at(iat)

      do ish = 1, nsh_id(izp)
         ii = iao_sh(is + ish)
         ni = nao_sh(is + ish)

         scale = alpha * g_onsri(ish, iat)
         if (abs(scale) < epsilon(1.0_wp)) cycle

         dst(ii+1:ii+ni, ii+1:ii+ni) = dst(ii+1:ii+ni, ii+1:ii+ni) + &
            & scale * src(ii+1:ii+ni, ii+1:ii+ni)
      end do
   end do

end subroutine onsite_ri_hadamard_add


!> Complex counterpart of shell_hadamard_add for Hermitian k-point blocks.
subroutine shell_hadamard_add_complex(nsh, nao_sh, iao_sh, g_sh, src, &
   & alpha, dst, adjoint_src)
   integer, intent(in) :: nsh, nao_sh(:), iao_sh(:)
   real(wp), intent(in) :: g_sh(:, :), alpha
   complex(wp), intent(in) :: src(:, :)
   complex(wp), intent(inout) :: dst(:, :)
   logical, intent(in), optional :: adjoint_src

   integer :: iao, ii, ish, jao, jj, jsh, ni, nj
   real(wp) :: scale
   logical :: adjoint

   adjoint = .false.
   if (present(adjoint_src)) adjoint = adjoint_src
   do ish = 1, nsh
      ii = iao_sh(ish)
      ni = nao_sh(ish)
      do jsh = 1, nsh
         jj = iao_sh(jsh)
         nj = nao_sh(jsh)
         scale = alpha*g_sh(jsh, ish)
         if (abs(scale) < epsilon(1.0_wp)) cycle
         if (.not. adjoint) then
            dst(jj+1:jj+nj, ii+1:ii+ni) = dst(jj+1:jj+nj, ii+1:ii+ni) &
               & + scale*src(jj+1:jj+nj, ii+1:ii+ni)
         else
            do iao = 1, ni
               do jao = 1, nj
                  dst(jj+jao, ii+iao) = dst(jj+jao, ii+iao) &
                     & + scale*conjg(src(ii+iao, jj+jao))
               end do
            end do
         end if
      end do
   end do
end subroutine shell_hadamard_add_complex


!> Accumulate the unconstrained real derivative of a shell-block kernel.
subroutine shell_parameter_gradient_complex(nsh, nao_sh, iao_sh, adjoint, &
   & source, alpha, gradient)
   integer, intent(in) :: nsh, nao_sh(:), iao_sh(:)
   complex(wp), intent(in) :: adjoint(:, :), source(:, :)
   real(wp), intent(in) :: alpha
   real(wp), intent(inout) :: gradient(:, :)

   integer :: iao, ii, ish, jao, jj, jsh, ni, nj

   do ish = 1, nsh
      ii = iao_sh(ish)
      ni = nao_sh(ish)
      do jsh = 1, nsh
         jj = iao_sh(jsh)
         nj = nao_sh(jsh)
         do iao = 1, ni
            do jao = 1, nj
               gradient(jj+jao, ii+iao) = gradient(jj+jao, ii+iao) &
                  & + alpha*real(conjg(adjoint(jj+jao, ii+iao)) &
                  & * source(jj+jao, ii+iao), wp)
            end do
         end do
      end do
   end do
end subroutine shell_parameter_gradient_complex


!> Convert independent shell-kernel derivatives to the symmetric convention
!> consumed by get_mulliken_derivs.
subroutine symmetrize_shell_parameter_gradient(nsh, nao_sh, iao_sh, raw, &
   & gradient)
   integer, intent(in) :: nsh, nao_sh(:), iao_sh(:)
   real(wp), intent(in) :: raw(:, :)
   real(wp), intent(out) :: gradient(:, :)

   integer :: iao, ii, ish, jao, jj, jsh, ni, nj
   real(wp) :: tmp

   gradient = raw
   do ish = 1, nsh
      ii = iao_sh(ish)
      ni = nao_sh(ish)
      do jsh = 1, ish-1
         jj = iao_sh(jsh)
         nj = nao_sh(jsh)
         do iao = 1, ni
            do jao = 1, nj
               tmp = raw(jj+jao, ii+iao) + raw(ii+iao, jj+jao)
               gradient(jj+jao, ii+iao) = tmp
               gradient(ii+iao, jj+jao) = tmp
            end do
         end do
      end do
   end do
end subroutine symmetrize_shell_parameter_gradient


!> Complex counterpart of atom_hadamard_add for Hermitian k-point blocks.
subroutine atom_hadamard_add_complex(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
   & g_at, src, alpha, dst, adjoint_src)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   real(wp), intent(in) :: g_at(:, :), alpha
   complex(wp), intent(in) :: src(:, :)
   complex(wp), intent(inout) :: dst(:, :)
   logical, intent(in), optional :: adjoint_src

   integer :: iao, iat, ii, is, ish, izp, jao, jat, jj, js, jsh, jzp, ni, nj
   real(wp) :: scale
   logical :: adjoint

   adjoint = .false.
   if (present(adjoint_src)) adjoint = adjoint_src
   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do jat = 1, nat
         jzp = id(jat)
         js = ish_at(jat)
         scale = alpha*g_at(jat, iat)
         if (abs(scale) < epsilon(1.0_wp)) cycle
         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            ni = nao_sh(is+ish)
            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js+jsh)
               nj = nao_sh(js+jsh)
               if (.not. adjoint) then
                  dst(jj+1:jj+nj, ii+1:ii+ni) = &
                     & dst(jj+1:jj+nj, ii+1:ii+ni) &
                     & + scale*src(jj+1:jj+nj, ii+1:ii+ni)
               else
                  do iao = 1, ni
                     do jao = 1, nj
                        dst(jj+jao, ii+iao) = dst(jj+jao, ii+iao) &
                           & + scale*conjg(src(ii+iao, jj+jao))
                     end do
                  end do
               end if
            end do
         end do
      end do
   end do
end subroutine atom_hadamard_add_complex


!> Accumulate the unconstrained real derivative of an atom-block kernel.
subroutine atom_parameter_gradient_complex(nat, id, nsh_id, nao_sh, ish_at, &
   & iao_sh, adjoint, source, alpha, gradient)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   complex(wp), intent(in) :: adjoint(:, :), source(:, :)
   real(wp), intent(in) :: alpha
   real(wp), intent(inout) :: gradient(:, :)

   integer :: iao, iat, ii, is, ish, izp, jao, jat, jj, js, jsh, jzp, ni, nj

   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do jat = 1, nat
         jzp = id(jat)
         js = ish_at(jat)
         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            ni = nao_sh(is+ish)
            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js+jsh)
               nj = nao_sh(js+jsh)
               do iao = 1, ni
                  do jao = 1, nj
                     gradient(jj+jao, ii+iao) = gradient(jj+jao, ii+iao) &
                        & + alpha*real(conjg(adjoint(jj+jao, ii+iao)) &
                        & * source(jj+jao, ii+iao), wp)
                  end do
               end do
            end do
         end do
      end do
   end do
end subroutine atom_parameter_gradient_complex


!> Convert independent atom-kernel derivatives to the symmetric convention
!> consumed by get_bocorr_derivs.
subroutine symmetrize_atom_parameter_gradient(nat, id, nsh_id, nao_sh, ish_at, &
   & iao_sh, raw, gradient)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   real(wp), intent(in) :: raw(:, :)
   real(wp), intent(out) :: gradient(:, :)

   integer :: iao, iat, ii, is, ish, izp, jao, jat, jj, js, jsh, jzp, ni, nj
   real(wp) :: tmp

   gradient = raw
   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do jat = 1, iat-1
         jzp = id(jat)
         js = ish_at(jat)
         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            ni = nao_sh(is+ish)
            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js+jsh)
               nj = nao_sh(js+jsh)
               do iao = 1, ni
                  do jao = 1, nj
                     tmp = raw(jj+jao, ii+iao) + raw(ii+iao, jj+jao)
                     gradient(jj+jao, ii+iao) = tmp
                     gradient(ii+iao, jj+jao) = tmp
                  end do
               end do
            end do
         end do
      end do
   end do
end subroutine symmetrize_atom_parameter_gradient


!> Complex counterpart of onsite_fx_hadamard_add.
subroutine onsite_fx_hadamard_add_complex(nat, id, nsh_id, nao_sh, ish_at, &
   & iao_sh, g_onsfx, src, alpha, dst, adjoint_src)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   real(wp), intent(in) :: g_onsfx(:, :, :), alpha
   complex(wp), intent(in) :: src(:, :)
   complex(wp), intent(inout) :: dst(:, :)
   logical, intent(in), optional :: adjoint_src

   integer :: iao, iat, ii, is, ish, izp, jao, jj, jsh, ni, nj
   real(wp) :: scale
   logical :: adjoint

   adjoint = .false.
   if (present(adjoint_src)) adjoint = adjoint_src
   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)
         do jsh = 1, nsh_id(izp)
            jj = iao_sh(is+jsh)
            nj = nao_sh(is+jsh)
            scale = alpha*g_onsfx(jsh, ish, iat)
            if (abs(scale) < epsilon(1.0_wp)) cycle
            if (.not. adjoint) then
               dst(jj+1:jj+nj, ii+1:ii+ni) = &
                  & dst(jj+1:jj+nj, ii+1:ii+ni) &
                  & + scale*src(jj+1:jj+nj, ii+1:ii+ni)
            else
               do iao = 1, ni
                  do jao = 1, nj
                     dst(jj+jao, ii+iao) = dst(jj+jao, ii+iao) &
                        & + scale*conjg(src(ii+iao, jj+jao))
                  end do
               end do
            end if
         end do
      end do
   end do
end subroutine onsite_fx_hadamard_add_complex


!> Complex counterpart of onsite_fx_symv.
subroutine onsite_fx_symv_complex(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
   & g_onsfx, xvec, yvec)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   real(wp), intent(in) :: g_onsfx(:, :, :)
   complex(wp), intent(in) :: xvec(:)
   complex(wp), intent(out) :: yvec(:)

   integer :: iat, ii, is, ish, izp, jj, jsh, ni, nj
   complex(wp) :: shsum

   yvec = (0.0_wp, 0.0_wp)
   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)
         shsum = sum(xvec(ii+1:ii+ni))
         if (abs(shsum) < epsilon(1.0_wp)) cycle
         do jsh = 1, nsh_id(izp)
            jj = iao_sh(is+jsh)
            nj = nao_sh(is+jsh)
            yvec(jj+1:jj+nj) = yvec(jj+1:jj+nj) &
               & + g_onsfx(jsh, ish, iat)*shsum
         end do
      end do
   end do
end subroutine onsite_fx_symv_complex


!> Complex counterpart of onsite_ri_hadamard_add.
subroutine onsite_ri_hadamard_add_complex(nat, id, nsh_id, nao_sh, ish_at, &
   & iao_sh, g_onsri, src, alpha, dst)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   real(wp), intent(in) :: g_onsri(:, :), alpha
   complex(wp), intent(in) :: src(:, :)
   complex(wp), intent(inout) :: dst(:, :)

   integer :: iat, ii, is, ish, izp, ni
   real(wp) :: scale

   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)
         scale = alpha*g_onsri(ish, iat)
         if (abs(scale) < epsilon(1.0_wp)) cycle
         dst(ii+1:ii+ni, ii+1:ii+ni) = &
            & dst(ii+1:ii+ni, ii+1:ii+ni) &
            & + scale*src(ii+1:ii+ni, ii+1:ii+ni)
      end do
   end do
end subroutine onsite_ri_hadamard_add_complex


!> Reverse an onsite exchange Hadamard map with respect to its real kernel.
subroutine onsite_fx_parameter_gradient_complex(nat, id, nsh_id, nao_sh, &
   & ish_at, iao_sh, adjoint, src, alpha, gradient, adjoint_src)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   complex(wp), intent(in) :: adjoint(:, :), src(:, :)
   real(wp), intent(in) :: alpha
   real(wp), intent(inout) :: gradient(:, :, :)
   logical, intent(in), optional :: adjoint_src

   integer :: iao, iat, ii, is, ish, izp, jao, jj, jsh, ni, nj
   real(wp) :: value
   logical :: take_adjoint

   take_adjoint = .false.
   if (present(adjoint_src)) take_adjoint = adjoint_src
   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)
         do jsh = 1, nsh_id(izp)
            jj = iao_sh(is+jsh)
            nj = nao_sh(is+jsh)
            value = 0.0_wp
            if (.not.take_adjoint) then
               value = real(sum(conjg(adjoint(jj+1:jj+nj, ii+1:ii+ni)) &
                  & *src(jj+1:jj+nj, ii+1:ii+ni)), wp)
            else
               do iao = 1, ni
                  do jao = 1, nj
                     value = value + real(conjg(adjoint(jj+jao, ii+iao)) &
                        & *conjg(src(ii+iao, jj+jao)), wp)
                  end do
               end do
            end if
            gradient(jsh, ish, iat) = gradient(jsh, ish, iat) + alpha*value
         end do
      end do
   end do
end subroutine onsite_fx_parameter_gradient_complex


!> Reverse an onsite rotational-invariance Hadamard map.
subroutine onsite_ri_parameter_gradient_complex(nat, id, nsh_id, nao_sh, &
   & ish_at, iao_sh, adjoint, src, alpha, gradient)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   complex(wp), intent(in) :: adjoint(:, :), src(:, :)
   real(wp), intent(in) :: alpha
   real(wp), intent(inout) :: gradient(:, :)

   integer :: iat, ii, is, ish, izp, ni

   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)
         gradient(ish, iat) = gradient(ish, iat) + alpha*real(sum( &
            & conjg(adjoint(ii+1:ii+ni, ii+1:ii+ni)) &
            & *src(ii+1:ii+ni, ii+1:ii+ni)), wp)
      end do
   end do
end subroutine onsite_ri_parameter_gradient_complex


!> Reverse the onsite shell-sum map with respect to its real kernel.
subroutine onsite_fx_symv_parameter_gradient_complex(nat, id, nsh_id, &
   & nao_sh, ish_at, iao_sh, adjoint, src, gradient)
   integer, intent(in) :: nat, id(:), nsh_id(:), nao_sh(:), ish_at(:), iao_sh(:)
   complex(wp), intent(in) :: adjoint(:), src(:)
   real(wp), intent(inout) :: gradient(:, :, :)

   integer :: iat, ii, is, ish, izp, jj, jsh, ni, nj
   complex(wp) :: adjoint_sum, source_sum

   do iat = 1, nat
      izp = id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)
         source_sum = sum(src(ii+1:ii+ni))
         do jsh = 1, nsh_id(izp)
            jj = iao_sh(is+jsh)
            nj = nao_sh(is+jsh)
            adjoint_sum = sum(adjoint(jj+1:jj+nj))
            gradient(jsh, ish, iat) = gradient(jsh, ish, iat) &
               & + real(conjg(adjoint_sum)*source_sum, wp)
         end do
      end do
   end do
end subroutine onsite_fx_symv_parameter_gradient_complex


!> Compute the kernel pair value for two shells.
pure function get_gmulliken_pair(r1, gam_ij, offdiag_scale_ij, &
   & hubbard_exp, hubbard_exp_r0, rad_ij, gexp, frscale, lrscale, omega) result(rsh)
   real(wp), intent(in) :: r1
   real(wp), intent(in) :: gam_ij
   real(wp), intent(in) :: offdiag_scale_ij
   real(wp), intent(in) :: hubbard_exp
   real(wp), intent(in) :: hubbard_exp_r0
   real(wp), intent(in) :: rad_ij
   real(wp), intent(in) :: gexp
   real(wp), intent(in) :: frscale
   real(wp), intent(in) :: lrscale
   real(wp), intent(in) :: omega
   real(wp) :: rsh

   real(wp) :: damp, gam, denom, r1g

   r1g = r1**gexp
   damp = exp(-(hubbard_exp + hubbard_exp_r0 * rad_ij) * r1)
   gam = gam_ij * offdiag_scale_ij / damp
   denom = (r1g + gam**(-gexp))**(1.0_wp/gexp)
   rsh = (frscale + lrscale * erf(omega * r1)) / denom

end function get_gmulliken_pair


!> Compute the derivative of the kernel pair value w.r.t. distance.
pure subroutine get_gmulliken_pair_deriv(r1, gam_ij, offdiag_scale_ij, &
   & hubbard_exp, hubbard_exp_r0, rad_ij, gexp, frscale, lrscale, omega, &
   & drsh)
   real(wp), intent(in) :: r1
   real(wp), intent(in) :: gam_ij
   real(wp), intent(in) :: offdiag_scale_ij
   real(wp), intent(in) :: hubbard_exp
   real(wp), intent(in) :: hubbard_exp_r0
   real(wp), intent(in) :: rad_ij
   real(wp), intent(in) :: gexp
   real(wp), intent(in) :: frscale
   real(wp), intent(in) :: lrscale
   real(wp), intent(in) :: omega
   real(wp), intent(out) :: drsh

   real(wp) :: damp, gam, exparg, r1g
   real(wp) :: rsh_factor, drsh_factor, denom, denom_pow, denom_deriv

   r1g = r1**gexp
   exparg = hubbard_exp + hubbard_exp_r0 * rad_ij
   damp = exp(-exparg * r1)
   gam = gam_ij * offdiag_scale_ij / damp

   ! Range-separation factor and its derivative
   rsh_factor = frscale + lrscale * erf(omega * r1)
   drsh_factor = lrscale * 2.0_wp * omega / sqrtpi * exp(-(omega * r1)**2)

   ! Damped coulomb interaction denominator
   denom = r1g + gam**(-gexp)
   denom_pow = denom**(1.0_wp / gexp)
   denom_deriv = denom_pow * denom

   drsh = drsh_factor / denom_pow &
      & + rsh_factor * (-r1g / r1 + exparg * gam**(-gexp)) / denom_deriv

end subroutine get_gmulliken_pair_deriv


end module tblite_exchange_fock
