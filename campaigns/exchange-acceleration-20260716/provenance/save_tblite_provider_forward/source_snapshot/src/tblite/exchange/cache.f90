
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

!> @file tblite/exchange/cache.f90
!> Provides a cache specific for all exchange interactions

!> Data container for mutable data in exchange calculations
module tblite_exchange_cache
   use mctc_env, only : wp, dp
   use mctc_io, only : structure_type
   use tblite_wignerseitz, only : wignerseitz_cell, new_wignerseitz_cell
   implicit none
   private

   public :: exchange_bvk_kernel, exchange_cache

   !> Image-resolved exchange kernels for a regular Born-von Karman mesh.
   !>
   !> This neutral cache representation is shared by the exchange evaluator
   !> and the CP2K compatibility plan without introducing module-global state.
   type :: exchange_bvk_kernel
      integer :: nmesh(3) = 1
      integer, allocatable :: reps(:, :)
      real(wp), allocatable :: g_mulliken_r(:, :, :)
      real(wp), allocatable :: g_bocorr_r(:, :, :)
   end type exchange_bvk_kernel

   type :: exchange_cache
      !> Wigner-Seitz cell for periodic calculations
      type(wignerseitz_cell) :: wsc

      !> Mulliken approximate exchange matrix: [nsh, nsh]
      real(wp), allocatable :: g_mulliken(:, :)
      !> Integrals for the onsite exchange correction: [maxsh, maxsh, nat]
      real(wp), allocatable :: g_onsfx(:, :, :)
      !> Charge-derivative of the integrals for the onsite exchange correction: [maxsh, maxsh, nsh]
      real(wp), allocatable :: dgdq_onsfx(:, :, :)
      !> Integrals for the onsite rotational invariance correction: [maxsh, nat]
      real(wp), allocatable :: g_onsri(:, :)
      !> Charge-derivative of the onsite rotational invariance correction matrix: [maxsh, nsh]
      real(wp), allocatable :: dgdq_onsri(:, :)
      !> Integrals for the bond-order correlation correction: [nat, nat]
      real(wp), allocatable :: g_bocorr(:, :)

      !> Previously calculated Fock matrix contribution 
      real(wp), allocatable :: prev_F(:, :, :)
      !> Previously calculated shell potential contribution
      real(wp), allocatable :: prev_vsh(:, :)

      !> Geometry- and mesh-static BvK kernel and Fourier plan
      logical :: bvk_plan_valid = .false.
      !> Number of successfully constructed BvK plans (diagnostic/test oracle)
      integer :: bvk_plan_builds = 0
      !> Exact exchange-model signature associated with the cached kernel
      logical :: bvk_model_valid = .false.
      integer :: bvk_model_nao = 0, bvk_model_nsh = 0, bvk_model_maxsh = 0
      integer, allocatable :: bvk_model_nsh_id(:), bvk_model_nao_sh(:)
      integer, allocatable :: bvk_model_ish_at(:), bvk_model_iao_sh(:)
      real(wp) :: bvk_model_frscale = 0.0_wp
      real(wp) :: bvk_model_omega = 0.0_wp
      real(wp) :: bvk_model_lrscale = 0.0_wp
      real(wp) :: bvk_model_ondiag_scale = 0.0_wp
      real(wp) :: bvk_model_hubbard_exp = 0.0_wp
      real(wp) :: bvk_model_hubbard_exp_r0 = 0.0_wp
      real(wp) :: bvk_model_gexp = 0.0_wp
      real(wp) :: bvk_model_corr_exp = 0.0_wp
      real(wp), allocatable :: bvk_model_hubbard(:, :, :, :)
      real(wp), allocatable :: bvk_model_offdiag_scale(:, :, :, :)
      real(wp), allocatable :: bvk_model_rad(:, :)
      real(wp), allocatable :: bvk_model_corr_scale(:, :)
      real(wp), allocatable :: bvk_model_corr_rad(:, :)
      logical :: bvk_periodic(3) = .false.
      integer :: bvk_nmesh(3) = 0
      integer, allocatable :: bvk_id(:)
      real(wp), allocatable :: bvk_xyz(:, :), bvk_lattice(:, :)
      real(wp), allocatable :: bvk_kfrac(:, :), bvk_weights(:)
      real(wp) :: bvk_twist(3) = 0.0_wp
      integer, allocatable :: bvk_input_to_grid(:), bvk_grid_to_input(:)
      complex(wp), allocatable :: bvk_phase_forward(:, :)
      complex(wp), allocatable :: bvk_phase_inverse(:, :)
      type(exchange_bvk_kernel) :: bvk_kernel
   contains
      procedure :: update
      procedure :: bvk_matches
      procedure :: set_bvk_signature
   end type exchange_cache


contains


subroutine update(self, mol)
   !> Instance of the electrostatic container
   class(exchange_cache), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   if (any(mol%periodic)) then
      call new_wignerseitz_cell(self%wsc, mol)
   end if

end subroutine update


!> Check whether the persistent BvK plan belongs to this exact geometry and mesh.
logical function bvk_matches(self, mol, nmesh, kfrac, weights) result(matches)
   class(exchange_cache), intent(in) :: self
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: nmesh(3)
   real(wp), intent(in) :: kfrac(:, :), weights(:)

   integer :: ik, nk

   matches = .false.
   if (.not.self%bvk_plan_valid) return
   if (.not.self%bvk_model_valid) return
   nk = size(weights)
   if (any(self%bvk_nmesh /= nmesh)) return
   if (any(self%bvk_periodic .neqv. mol%periodic)) return
   if (.not.allocated(self%bvk_id) .or. .not.allocated(self%bvk_xyz) &
      & .or. .not.allocated(self%bvk_lattice) &
      & .or. .not.allocated(self%bvk_kfrac) &
      & .or. .not.allocated(self%bvk_weights) &
      & .or. .not.allocated(self%bvk_input_to_grid) &
      & .or. .not.allocated(self%bvk_grid_to_input) &
      & .or. .not.allocated(self%bvk_phase_forward) &
      & .or. .not.allocated(self%bvk_phase_inverse) &
      & .or. .not.allocated(self%bvk_kernel%reps) &
      & .or. .not.allocated(self%bvk_kernel%g_mulliken_r) &
      & .or. .not.allocated(self%bvk_kernel%g_bocorr_r)) return
   if (size(self%bvk_id) /= size(mol%id)) return
   if (any(shape(self%bvk_xyz) /= shape(mol%xyz))) return
   if (any(shape(self%bvk_lattice) /= shape(mol%lattice))) return
   if (any(shape(self%bvk_kfrac) /= shape(kfrac))) return
   if (size(self%bvk_weights) /= size(weights)) return
   if (size(self%bvk_input_to_grid) /= nk) return
   if (size(self%bvk_grid_to_input) /= nk) return
   if (any(self%bvk_kernel%nmesh /= nmesh)) return
   if (any(shape(self%bvk_kernel%reps) /= [3, nk])) return
   if (size(self%bvk_kernel%g_mulliken_r, 1) /= &
      & size(self%bvk_kernel%g_mulliken_r, 2)) return
   if (size(self%bvk_kernel%g_mulliken_r, 3) /= nk) return
   if (any(shape(self%bvk_kernel%g_bocorr_r) /= [mol%nat, mol%nat, nk])) return
   if (any(shape(self%bvk_phase_forward) /= &
      & [nk, nk])) return
   if (any(shape(self%bvk_phase_inverse) /= &
      & [nk, nk])) return
   if (any(self%bvk_input_to_grid < 1) &
      & .or. any(self%bvk_input_to_grid > nk)) return
   if (any(self%bvk_grid_to_input < 1) &
      & .or. any(self%bvk_grid_to_input > nk)) return
   do ik = 1, nk
      if (self%bvk_grid_to_input(self%bvk_input_to_grid(ik)) /= ik) return
   end do
   if (any(self%bvk_id /= mol%id)) return
   if (any(self%bvk_xyz /= mol%xyz)) return
   if (any(self%bvk_lattice /= mol%lattice)) return
   if (any(self%bvk_kfrac /= kfrac)) return
   if (any(self%bvk_weights /= weights)) return
   matches = .true.
end function bvk_matches


!> Record the exact fingerprint after a BvK kernel and Fourier plan were validated.
subroutine set_bvk_signature(self, mol, nmesh, kfrac, weights)
   class(exchange_cache), intent(inout) :: self
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: nmesh(3)
   real(wp), intent(in) :: kfrac(:, :), weights(:)

   self%bvk_plan_valid = .false.
   self%bvk_nmesh = nmesh
   self%bvk_periodic = mol%periodic
   if (allocated(self%bvk_id)) deallocate(self%bvk_id)
   if (allocated(self%bvk_xyz)) deallocate(self%bvk_xyz)
   if (allocated(self%bvk_lattice)) deallocate(self%bvk_lattice)
   if (allocated(self%bvk_kfrac)) deallocate(self%bvk_kfrac)
   if (allocated(self%bvk_weights)) deallocate(self%bvk_weights)
   allocate(self%bvk_id, source=mol%id)
   allocate(self%bvk_xyz, source=mol%xyz)
   allocate(self%bvk_lattice, source=mol%lattice)
   allocate(self%bvk_kfrac, source=kfrac)
   allocate(self%bvk_weights, source=weights)
   self%bvk_plan_builds = self%bvk_plan_builds + 1
   self%bvk_plan_valid = .true.
end subroutine set_bvk_signature

end module tblite_exchange_cache
