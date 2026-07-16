!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief types for tblite
!> \author JVP
!> \history creation 09.2024
! **************************************************************************************************

MODULE tblite_types

#if defined(__TBLITE)
   USE mctc_io, ONLY: structure_type, new
#if defined(__TBLITE_GXTB)
   USE tblite_acp_cache, ONLY: acp_cache
   USE tblite_basis_cache, ONLY: basis_cache
   USE tblite_cp2k_compat, ONLY: cp2k_potential_mixer_type
#endif
   USE tblite_xtb_calculator, ONLY: xtb_calculator
   USE tblite_container, ONLY: container_cache
   USE tblite_scc_mixer, ONLY: cp2k_tblite_broyden_mixer_type
   USE tblite_wavefunction_type, ONLY: wavefunction_type
   USE tblite_scf_potential, ONLY: potential_type
   USE tblite_param, ONLY: param_record
#endif

   USE cp_dbcsr_api, ONLY: dbcsr_p_type
   USE cp_dbcsr_operations, ONLY: dbcsr_deallocate_matrix_set
   USE input_constants, ONLY: tblite_mixer_damping_default, &
                              tblite_mixer_iterations_default, &
                              tblite_mixer_max_weight_default, &
                              tblite_mixer_min_weight_default, &
                              tblite_mixer_omega0_default, &
                              tblite_mixer_weight_factor_default, &
                              tblite_solver_gvd
   USE kinds, ONLY: dp

#include "./base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_types'

   PUBLIC :: tblite_type, deallocate_tblite_type, allocate_tblite_type
#if defined(__TBLITE_GXTB)
   PUBLIC :: tb_gxtb_kmesh_map_type

   TYPE tb_gxtb_kmesh_map_type
      CHARACTER(LEN=32)                                :: scheme = ""
      INTEGER                                          :: nfull = 0, nred = 0
      LOGICAL                                          :: explicit_general = .FALSE.
      LOGICAL                                          :: full_grid = .FALSE., gamma_centered = .FALSE.
      LOGICAL                                          :: identity = .FALSE., inversion_only = .FALSE.
      INTEGER                                          :: symmetry_backend = 0, symmetry_method = 0
      INTEGER, DIMENSION(3)                            :: nmesh = 1
      REAL(KIND=dp), DIMENSION(3)                      :: kp_shift = 0.0_dp
      INTEGER, ALLOCATABLE, DIMENSION(:)               :: op_count, representative, source, symop
      INTEGER, ALLOCATABLE, DIMENSION(:, :)            :: rotation_slot, symops
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)      :: reduced_xkp, xkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)         :: reduced_weight, weight, fold_factor
   END TYPE tb_gxtb_kmesh_map_type
#endif

   TYPE tblite_type

      LOGICAL                                            :: use_virial = .FALSE.
      INTEGER, ALLOCATABLE                               :: el_num(:)

      REAL(KIND=dp)                                      :: accuracy = 1.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)                     :: sigma = -1.0_dp
      REAL(KIND=dp), ALLOCATABLE                         :: e_hal(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_rep(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_disp(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_scd(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_es(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_int(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_inc(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_exchange(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_spin(:)
      REAL(KIND=dp), ALLOCATABLE                         :: selfenergy(:)
      REAL(KIND=dp), ALLOCATABLE                         :: cn(:)

      REAL(KIND=dp), ALLOCATABLE                         :: grad(:, :)
      !derivatives w. r. t. coordination number
      REAL(KIND=dp), ALLOCATABLE                         :: dsedcn(:)
      REAL(KIND=dp), ALLOCATABLE                         :: dcndr(:, :, :)
      REAL(KIND=dp), ALLOCATABLE                         :: dcndL(:, :, :)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipbra => Null(), &
                                                            dipket => Null()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: quadbra => Null(), &
                                                            quadket => Null()
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp_ref => Null()
      INTEGER                                            :: mixer_memory = tblite_mixer_iterations_default
      INTEGER                                            :: mixer_solver = tblite_solver_gvd
      REAL(KIND=dp)                                      :: mixer_damping = tblite_mixer_damping_default
      REAL(KIND=dp)                                      :: mixer_max_weight = tblite_mixer_max_weight_default
      REAL(KIND=dp)                                      :: mixer_min_weight = tblite_mixer_min_weight_default
      REAL(KIND=dp)                                      :: mixer_omega0 = tblite_mixer_omega0_default
      REAL(KIND=dp)                                      :: mixer_weight_factor = tblite_mixer_weight_factor_default

#if defined(__TBLITE_GXTB)
      TYPE(basis_cache)                                  :: bcache = basis_cache()
      TYPE(acp_cache)                                    :: acache = acp_cache()
      REAL(KIND=dp), ALLOCATABLE                         :: overlap_dense(:, :)
      TYPE(tb_gxtb_kmesh_map_type), ALLOCATABLE          :: gxtb_kmesh_map
      COMPLEX(KIND=dp), ALLOCATABLE                      :: acp_fock_kp(:, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: qualification_acp_fock_full(:, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: exchange_fock_kp(:, :, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: mixed_fock_kp(:, :, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: previous_density_kp(:, :, :, :)
      LOGICAL                                            :: gxtb_density_map_validated = .FALSE.
      LOGICAL                                            :: gxtb_transform_validated = .FALSE.
      LOGICAL                                            :: acp_fock_kp_valid = .FALSE.
      LOGICAL                                            :: exchange_fock_kp_valid = .FALSE.
      LOGICAL                                            :: mixed_fock_kp_valid = .FALSE.
      LOGICAL                                            :: previous_density_kp_valid = .FALSE.
      REAL(KIND=dp)                                      :: e_acp_kp = 0.0_dp
      INTEGER                                            :: potential_mixer_dimension = 0
      TYPE(cp2k_potential_mixer_type), ALLOCATABLE       :: potential_mixer
      TYPE(wavefunction_type)                            :: wfn_aux = wavefunction_type()
      TYPE(container_cache)                              :: ecache = container_cache()
      TYPE(container_cache)                              :: rcache = container_cache()
      TYPE(container_cache)                              :: scache = container_cache()
#endif
#if defined(__TBLITE)
      TYPE(structure_type)                               :: mol = structure_type()
      TYPE(xtb_calculator)                               :: calc = xtb_calculator()
      TYPE(container_cache)                              :: cache = container_cache()
      TYPE(container_cache)                              :: dcache = container_cache()
      TYPE(container_cache)                              :: icache = container_cache()
      TYPE(wavefunction_type)                            :: wfn = wavefunction_type()
      TYPE(potential_type)                               :: pot = potential_type()
      TYPE(cp2k_tblite_broyden_mixer_type), ALLOCATABLE  :: mixer
      TYPE(param_record), ALLOCATABLE                    :: param
#endif

   END TYPE tblite_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param tb_tblite ...
! **************************************************************************************************
   SUBROUTINE allocate_tblite_type(tb_tblite)

      TYPE(tblite_type), POINTER                         :: tb_tblite

      CALL deallocate_tblite_type(tb_tblite)

      ALLOCATE (tb_tblite)

   END SUBROUTINE allocate_tblite_type

! **************************************************************************************************
!> \brief ...
!> \param tb_tblite ...
! **************************************************************************************************
   SUBROUTINE deallocate_tblite_type(tb_tblite)

      TYPE(tblite_type), POINTER                         :: tb_tblite

      IF (ASSOCIATED(tb_tblite)) THEN
         IF (ALLOCATED(tb_tblite%el_num)) DEALLOCATE (tb_tblite%el_num)

         IF (ALLOCATED(tb_tblite%e_hal)) DEALLOCATE (tb_tblite%e_hal)
         IF (ALLOCATED(tb_tblite%e_rep)) DEALLOCATE (tb_tblite%e_rep)
         IF (ALLOCATED(tb_tblite%e_disp)) DEALLOCATE (tb_tblite%e_disp)
         IF (ALLOCATED(tb_tblite%e_scd)) DEALLOCATE (tb_tblite%e_scd)
         IF (ALLOCATED(tb_tblite%e_es)) DEALLOCATE (tb_tblite%e_es)
         IF (ALLOCATED(tb_tblite%e_int)) DEALLOCATE (tb_tblite%e_int)
         IF (ALLOCATED(tb_tblite%e_inc)) DEALLOCATE (tb_tblite%e_inc)
         IF (ALLOCATED(tb_tblite%e_exchange)) DEALLOCATE (tb_tblite%e_exchange)
         IF (ALLOCATED(tb_tblite%e_spin)) DEALLOCATE (tb_tblite%e_spin)
         IF (ALLOCATED(tb_tblite%selfenergy)) DEALLOCATE (tb_tblite%selfenergy)
         IF (ALLOCATED(tb_tblite%cn)) DEALLOCATE (tb_tblite%cn)

         IF (ALLOCATED(tb_tblite%grad)) DEALLOCATE (tb_tblite%grad)
         IF (ALLOCATED(tb_tblite%dsedcn)) DEALLOCATE (tb_tblite%dsedcn)
         IF (ALLOCATED(tb_tblite%dcndr)) DEALLOCATE (tb_tblite%dcndr)
         IF (ALLOCATED(tb_tblite%dcndL)) DEALLOCATE (tb_tblite%dcndL)
#if defined(__TBLITE_GXTB)
         IF (ALLOCATED(tb_tblite%overlap_dense)) DEALLOCATE (tb_tblite%overlap_dense)
         IF (ALLOCATED(tb_tblite%gxtb_kmesh_map)) DEALLOCATE (tb_tblite%gxtb_kmesh_map)
         IF (ALLOCATED(tb_tblite%acp_fock_kp)) DEALLOCATE (tb_tblite%acp_fock_kp)
         IF (ALLOCATED(tb_tblite%qualification_acp_fock_full)) &
            DEALLOCATE (tb_tblite%qualification_acp_fock_full)
         IF (ALLOCATED(tb_tblite%exchange_fock_kp)) DEALLOCATE (tb_tblite%exchange_fock_kp)
         IF (ALLOCATED(tb_tblite%mixed_fock_kp)) DEALLOCATE (tb_tblite%mixed_fock_kp)
         IF (ALLOCATED(tb_tblite%previous_density_kp)) DEALLOCATE (tb_tblite%previous_density_kp)
         tb_tblite%acp_fock_kp_valid = .FALSE.
         tb_tblite%exchange_fock_kp_valid = .FALSE.
         tb_tblite%mixed_fock_kp_valid = .FALSE.
         tb_tblite%previous_density_kp_valid = .FALSE.
         IF (ALLOCATED(tb_tblite%potential_mixer)) DEALLOCATE (tb_tblite%potential_mixer)
         tb_tblite%potential_mixer_dimension = 0
#endif

         IF (ASSOCIATED(tb_tblite%dipbra)) &
            CALL dbcsr_deallocate_matrix_set(tb_tblite%dipbra)
         IF (ASSOCIATED(tb_tblite%dipket)) &
            CALL dbcsr_deallocate_matrix_set(tb_tblite%dipket)
         IF (ASSOCIATED(tb_tblite%quadbra)) &
            CALL dbcsr_deallocate_matrix_set(tb_tblite%quadbra)
         IF (ASSOCIATED(tb_tblite%quadket)) &
            CALL dbcsr_deallocate_matrix_set(tb_tblite%quadket)
         IF (ASSOCIATED(tb_tblite%rho_ao_kp_ref)) &
            CALL dbcsr_deallocate_matrix_set(tb_tblite%rho_ao_kp_ref)
#if defined(__TBLITE)
         IF (ALLOCATED(tb_tblite%mixer)) DEALLOCATE (tb_tblite%mixer)
         IF (ALLOCATED(tb_tblite%param)) DEALLOCATE (tb_tblite%param)
#endif

         DEALLOCATE (tb_tblite)
      END IF

   END SUBROUTINE deallocate_tblite_type

END MODULE tblite_types
