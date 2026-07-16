!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief interface to tblite
!> \author JVP
!> \history creation 09.2024
! **************************************************************************************************

MODULE tblite_interface

   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite

#if defined(__TBLITE)
   USE mctc_env, ONLY: error_type
   USE mctc_io, ONLY: structure_type, new
   USE mctc_io_symbols, ONLY: symbol_to_number
   USE tblite_adjlist, ONLY: adjacency_list, new_adjacency_list
   USE tblite_basis_type, ONLY: get_cutoff
   USE tblite_container, ONLY: container_cache
   USE tblite_container_type, ONLY: container_type
   USE tblite_cutoff, ONLY: get_lattice_points
   USE tblite_data_spin, ONLY: get_spin_constant
#if defined(__TBLITE_GXTB)
   USE tblite_cp2k_compat, ONLY: get_mixer_dimension => cp2k_get_scc_mixer_dimension, &
                                 acp_image_type, hamiltonian_image_type, &
                                 cp2k_acp_kmesh, cp2k_acp_kmesh_gradient, &
                                 cp2k_apply_basis_gradient, cp2k_build_acp_projectors, &
                                 cp2k_build_gamma_h0, &
                                 cp2k_build_image_h0, cp2k_exchange_kmesh, &
                                 cp2k_exchange_kmesh_gradient, &
                                 cp2k_exchange_kpoint_gradient, &
                                 cp2k_exchange_stream_batch_advance, &
                                 cp2k_exchange_stream_batch_apply, &
                                 cp2k_exchange_stream_batch_bounds, &
                                 cp2k_exchange_stream_batch_pull, &
                                 cp2k_exchange_stream_batch_result, &
                                 cp2k_exchange_stream_apply, cp2k_exchange_stream_begin, &
                                 cp2k_exchange_stream_batched, &
                                 cp2k_exchange_stream_end, cp2k_exchange_stream_get, &
                                 cp2k_exchange_stream_has_full_mesh_storage, &
                                 cp2k_exchange_stream_peak_complex_elements, &
                                 cp2k_exchange_stream_oracle, cp2k_exchange_stream_push, &
                                 cp2k_exchange_stream_reduced, cp2k_exchange_stream_type, &
                                 cp2k_exchange_stream_reverse_batch_advance, &
                                 cp2k_exchange_stream_reverse_batch_apply, &
                                 cp2k_exchange_stream_reverse_batch_pull, &
                                 cp2k_exchange_stream_reverse_push, &
                                 cp2k_exchange_stream_reverse_result, &
                                 cp2k_exchange_transform_dense, &
                                 cp2k_exchange_transform_separable_dft, &
                                 cp2k_gamma_gradient, cp2k_image_h0_gradient, &
                                 cp2k_get_acp_cutoff, cp2k_get_cgto, cp2k_multipole_cgto, &
                                 cp2k_multipole_grad_cgto, cp2k_set_scf_options, &
                                 cp2k_prepare_exchange, cp2k_update_basis, &
                                 new_cp2k_potential_mixer
#else
   USE tblite_integral_multipole, ONLY: multipole_cgto, multipole_grad_cgto
#endif
   USE tblite_integral_multipole, ONLY: maxl, msao
   USE tblite_integral_type, ONLY: integral_type, new_integral
#if !defined(__TBLITE_GXTB)
   USE tblite_scf, ONLY: get_mixer_dimension
#endif
   USE tblite_scf_info, ONLY: scf_info, atom_resolved, shell_resolved, not_used
   USE tblite_scf_potential, ONLY: potential_type, new_potential, add_pot_to_h1
   USE tblite_spin, ONLY: spin_polarization, new_spin_polarization
   USE tblite_wavefunction_type, ONLY: wavefunction_type, new_wavefunction
   USE tblite_xtb_calculator, ONLY: xtb_calculator, new_xtb_calculator
   USE tblite_xtb_gfn1, ONLY: new_gfn1_calculator
   USE tblite_xtb_gfn2, ONLY: new_gfn2_calculator
#if defined(__TBLITE_GXTB)
   USE tblite_xtb_gxtb, ONLY: new_gxtb_calculator
#endif
   USE tblite_xtb_h0, ONLY: get_selfenergy, get_hamiltonian, get_occupation, &
                            get_hamiltonian_gradient, tb_hamiltonian
   USE tblite_xtb_ipea1, ONLY: new_ipea1_calculator
#endif
   USE ai_contraction, ONLY: block_add, &
                             contraction
   USE ai_overlap, ONLY: overlap_ab
   USE atomic_kind_types, ONLY: atomic_kind_type, get_atomic_kind, get_atomic_kind_set
   USE atprop_types, ONLY: atprop_type
   USE basis_set_types, ONLY: gto_basis_set_type, gto_basis_set_p_type, &
                              & allocate_gto_basis_set, write_gto_basis_set, process_gto_basis
   USE block_p_types, ONLY: block_p_type
   USE cell_types, ONLY: cell_type, get_cell
   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_control_types, ONLY: dft_control_type, xtb_reference_cli_type
   USE cp_dbcsr_api, ONLY: dbcsr_type, dbcsr_p_type, dbcsr_create, dbcsr_add, dbcsr_copy, &
                           dbcsr_get_block_p, dbcsr_finalize, &
                           dbcsr_deallocate_matrix, dbcsr_set, &
                           dbcsr_type_antisymmetric, dbcsr_type_symmetric, &
                           dbcsr_iterator_type, dbcsr_iterator_blocks_left, &
                           dbcsr_iterator_start, dbcsr_iterator_stop, &
                           dbcsr_iterator_next_block
   USE cp_dbcsr_contrib, ONLY: dbcsr_dot, dbcsr_print
   USE cp_dbcsr_cp2k_link, ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations, ONLY: dbcsr_allocate_matrix_set, dbcsr_deallocate_matrix_set
   USE cp_fm_struct, ONLY: cp_fm_struct_type
   USE cp_fm_types, ONLY: copy_info_type, cp_fm_cleanup_copy_general, &
                          cp_fm_create, cp_fm_finish_copy_general, cp_fm_get_info, &
                          cp_fm_release, cp_fm_start_copy_general, cp_fm_type
   USE cp_fm_pool_types, ONLY: cp_fm_pool_p_type, &
                               fm_pool_create_fm, fm_pool_give_back_fm
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_type, cp_logger_get_default_io_unit
   USE cp_output_handling, ONLY: cp_print_key_should_output, &
                                 cp_print_key_unit_nr, cp_print_key_finished_output, &
                                 debug_print_level, high_print_level, silent_print_level
   USE cp_units, ONLY: cp_unit_from_cp2k
#if defined(__TBLITE_GXTB)
   USE cryssym, ONLY: full_grid_gen
#endif
   USE input_constants, ONLY: gfn1xtb, gfn2xtb, gxtb, ipea1xtb, smear_energy_window, &
                              smear_fermi_dirac, smear_gaussian, smear_list, smear_mp, smear_mv, &
                              tblite_cli_born_kernel_auto, tblite_cli_born_kernel_p16, &
                              tblite_cli_born_kernel_still, tblite_cli_solution_state_bar1mol, &
                              tblite_cli_solution_state_gsolv, tblite_cli_solution_state_reference, &
                              tblite_cli_solvation_alpb, tblite_cli_solvation_cpcm, &
                              tblite_cli_solvation_gb, tblite_cli_solvation_gbe, &
                              tblite_cli_solvation_gbsa, &
                              tblite_guess_auto, tblite_guess_ceh, tblite_guess_eeq, tblite_guess_sad, &
                              tblite_mixer_damping_default, &
                              tblite_mixer_max_weight_default, tblite_mixer_min_weight_default, &
                              tblite_mixer_omega0_default, tblite_mixer_weight_factor_default, &
                              tblite_scc_mixer_auto, tblite_scc_mixer_cp2k, &
                              tblite_scc_mixer_none, tblite_scc_mixer_tblite, &
                              tblite_solver_gvd, tblite_solver_gvr
   USE input_section_types, ONLY: section_vals_val_get
   USE kinds, ONLY: default_path_length, default_string_length, dp, int_8
   USE kpoint_methods, ONLY: calibrate_symmetry_phases, kpoint_density_transform, rskp_transform
   USE kpoint_types, ONLY: get_kpoint_info, kind_rotmat_type, kpoint_env_type, &
                           kpoint_phase_gauge_input, kpoint_phase_gauge_pbc, &
                           kpoint_sym_type, kpoint_type
#if defined(__TBLITE_GXTB)
   USE mathconstants, ONLY: twopi
#endif
   USE memory_utilities, ONLY: reallocate
   USE message_passing, ONLY: mp_para_env_type
   USE mulliken, ONLY: ao_charges
   USE orbital_pointers, ONLY: ncoset
   USE particle_types, ONLY: particle_type
   USE qs_charge_mixing, ONLY: charge_mixing
   USE qs_condnum, ONLY: overlap_condnum
   USE qs_density_mixing_types, ONLY: direct_mixing_nr
   USE qs_energy_types, ONLY: qs_energy_type
   USE qs_environment_types, ONLY: get_qs_env, qs_environment_type
   USE qs_force_types, ONLY: qs_force_type, total_qs_force
   USE qs_integral_utils, ONLY: basis_set_list_setup, get_memory_usage
   USE qs_kind_types, ONLY: get_qs_kind, qs_kind_type, get_qs_kind_set
   USE qs_ks_types, ONLY: get_ks_env, qs_ks_env_type, set_ks_env
   USE qs_matrix_pools, ONLY: mpools_get
   USE qs_neighbor_list_types, ONLY: neighbor_list_iterator_create, neighbor_list_iterate, &
                                     get_iterator_info, neighbor_list_set_p_type, &
                                     neighbor_list_iterator_p_type, neighbor_list_iterator_release
   USE qs_overlap, ONLY: build_overlap_matrix, create_sab_matrix
   USE cp_dbcsr_output, ONLY: cp_dbcsr_write_sparse_matrix
   USE qs_rho_types, ONLY: qs_rho_get, qs_rho_type
   USE qs_scf_types, ONLY: qs_scf_env_type
   USE scf_control_types, ONLY: scf_control_type
   USE input_section_types, ONLY: section_vals_get_subs_vals, section_vals_type
   USE string_utilities, ONLY: integer_to_string
#if defined(__TBLITE)
   USE tblite_scc_mixer, ONLY: new_cp2k_tblite_mixer
#endif
   USE tblite_types, ONLY: tblite_type, deallocate_tblite_type, allocate_tblite_type
#if defined(__TBLITE_GXTB)
   USE tblite_types, ONLY: tb_gxtb_kmesh_map_type
#endif
   USE virial_types, ONLY: virial_type
   USE xtb_types, ONLY: get_xtb_atom_param, xtb_atom_type
   USE xtb_types, ONLY: xtb_atom_type

!$ USE OMP_LIB, ONLY: omp_destroy_lock, omp_init_lock, omp_set_lock, &
!$                    omp_unset_lock, omp_lock_kind

#include "./base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_interface'

   INTEGER, PARAMETER                                 :: dip_n = 3
   INTEGER, PARAMETER                                 :: quad_n = 6
#if defined(__TBLITE_GXTB)
   ! CP2K-only qualification mode.  The provider sees an ordinary bounded
   ! transaction, while CP2K assigns one disjoint additive BvK image range to
   ! every participating MPI rank and reduces the folded response afterwards.
   INTEGER, PARAMETER                                 :: tb_gxtb_exchange_stream_distributed = 4
   ! Exact CP2K k-group ownership precursor.  Only the owner of an irreducible
   ! k point expands its star blocks; official inter-k-point communicators feed
   ! one nonreplicated provider stream, and owner groups perform the weighted
   ! foldback.  The coupled provider kernel itself is intentionally not claimed
   ! to scale over k groups until save_tblite exposes a public partial-k-to-R
   ! accumulator merge/finalize interface.
   INTEGER, PARAMETER                                 :: tb_gxtb_exchange_stream_kgroup_owner = 5
   INTEGER, PARAMETER                                 :: tb_gxtb_exchange_gradient_dense = 0
   INTEGER, PARAMETER                                 :: tb_gxtb_exchange_gradient_streamed = 1
   INTEGER, PARAMETER                                 :: tb_gxtb_exchange_gradient_qualify = 2
#endif
   REAL(KIND=dp), PARAMETER                           :: same_atom = 0.00001_dp
   REAL(KIND=dp), PARAMETER                           :: tblite_scc_pconv = 2.0E-5_dp
   PUBLIC :: tb_set_calculator, tb_init_geometry, tb_init_wf
   PUBLIC :: tb_get_basis, build_tblite_matrices
   PUBLIC :: tb_get_energy, tb_update_charges, tb_ham_add_coulomb
   PUBLIC :: tb_native_scc_mixer_active
   PUBLIC :: tb_scf_mixer_error
   PUBLIC :: tb_get_multipole
   PUBLIC :: tb_derive_dH_off
   PUBLIC :: tb_gxtb_gamma_gradient, tb_gxtb_kpoint_gradient
#if defined(__TBLITE_GXTB)
   PUBLIC :: tb_apply_gxtb_kpoint_fm, tb_mix_gxtb_gamma_fm, tb_prepare_gxtb_kpoint_mixer, &
             tb_sync_gxtb_kpoint_density, tb_gxtb_kpoint_exchange_gradient
#endif
   PUBLIC :: tb_reference_cli_compare

CONTAINS

#if defined(__TBLITE)
! **************************************************************************************************
!> \brief Project a tblite charge/magnetization quantity to a CP2K spin channel.
!> \param values ...
!> \param ispin ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION tb_spin_project(values, ispin) RESULT(value)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: values
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp)                                      :: value

      value = values(1)
      IF (SIZE(values) > 1) THEN
         SELECT CASE (ispin)
         CASE (1)
            value = values(1) + values(2)
         CASE (2)
            value = values(1) - values(2)
         CASE DEFAULT
            value = values(1)
         END SELECT
      END IF

   END FUNCTION tb_spin_project

! **************************************************************************************************
!> \brief Store a private copy of the converged density for tblite force derivatives.
!> \param tb ...
!> \param matrix_p ...
! **************************************************************************************************
   SUBROUTINE tb_store_density_ref(tb, matrix_p)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p

      INTEGER                                            :: img, ispin, nimg, nspin

      nspin = SIZE(matrix_p, 1)
      nimg = SIZE(matrix_p, 2)
      IF (ASSOCIATED(tb%rho_ao_kp_ref)) THEN
         IF (SIZE(tb%rho_ao_kp_ref, 1) /= nspin .OR. SIZE(tb%rho_ao_kp_ref, 2) /= nimg) &
            CALL dbcsr_deallocate_matrix_set(tb%rho_ao_kp_ref)
      END IF
      IF (.NOT. ASSOCIATED(tb%rho_ao_kp_ref)) THEN
         CALL dbcsr_allocate_matrix_set(tb%rho_ao_kp_ref, nspin, nimg)
         DO img = 1, nimg
            DO ispin = 1, nspin
               ALLOCATE (tb%rho_ao_kp_ref(ispin, img)%matrix)
            END DO
         END DO
      END IF
      DO img = 1, nimg
         DO ispin = 1, nspin
            CALL dbcsr_copy(tb%rho_ao_kp_ref(ispin, img)%matrix, matrix_p(ispin, img)%matrix)
         END DO
      END DO

   END SUBROUTINE tb_store_density_ref
#endif

! **************************************************************************************************
!> \brief intialize geometry objects ...
!> \param qs_env ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_init_geometry(qs_env, tb)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb

#if defined(__TBLITE)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tblite_init_geometry'

      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER                                            :: iatom, natom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xyz
      INTEGER                                            :: handle, ikind
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL, DIMENSION(3)                              :: lperiod

      CALL timeset(routineN, handle)

      !get info from environment vaiarable
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, qs_kind_set=qs_kind_set)

      !get information about particles
      natom = SIZE(particle_set)
      ALLOCATE (xyz(3, natom))
      CALL allocate_tblite_type(tb)
      ALLOCATE (tb%el_num(natom))
      tb%el_num = -9
      DO iatom = 1, natom
         xyz(:, iatom) = particle_set(iatom)%r(:)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), zatom=tb%el_num(iatom))
         IF (tb%el_num(iatom) < 1 .OR. tb%el_num(iatom) > 85) THEN
            CPABORT("only elements 1-85 are supported by tblite")
         END IF
      END DO

      !get information about cell / lattice
      CALL get_cell(cell=cell, periodic=periodic)
      lperiod(1) = periodic(1) == 1
      lperiod(2) = periodic(2) == 1
      lperiod(3) = periodic(3) == 1

      !prepare for the call to the dispersion function
      CALL new(tb%mol, tb%el_num, xyz, lattice=cell%hmat, periodic=lperiod)

      DEALLOCATE (xyz)

      CALL timestop(handle)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_init_geometry

! **************************************************************************************************
!> \brief updating coordinates...
!> \param qs_env ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_update_geometry(qs_env, tb)

      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tblite_type)                                  :: tb

#if defined(__TBLITE)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tblite_update_geometry'

      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER                                            :: iatom, natom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xyz
      INTEGER                                            :: handle
#if defined(__TBLITE_GXTB)
      LOGICAL                                            :: geometry_changed
#endif

      CALL timeset(routineN, handle)

      !get info from environment vaiarable
      NULLIFY (cell)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell)

      !get information about particles
      natom = SIZE(particle_set)
      ALLOCATE (xyz(3, natom))
      DO iatom = 1, natom
         xyz(:, iatom) = particle_set(iatom)%r(:)
      END DO
#if defined(__TBLITE_GXTB)
      geometry_changed = ANY(tb%mol%xyz /= xyz) .OR. &
                         ANY(tb%mol%lattice /= cell%hmat)
      IF (geometry_changed) CALL tb_reset_gxtb_kpoint_cache(tb)
#endif
      tb%mol%xyz(:, :) = xyz
      tb%mol%lattice(:, :) = cell%hmat

      DEALLOCATE (xyz)

      CALL timestop(handle)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_update_geometry

! **************************************************************************************************
!> \brief initialize wavefunction ...
!> \param tb ...
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE tb_init_wf(tb, dft_control)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(dft_control_type), POINTER                    :: dft_control

#if defined(__TBLITE)

      INTEGER                                            :: nSpin

      TYPE(scf_info)                                     :: info

      nSpin = dft_control%nspins
      IF (nSpin /= 1 .AND. nSpin /= 2) CPABORT("tblite supports only one or two spin channels")

      tb%mol%charge = dft_control%charge
      tb%mol%uhf = MAX(0, dft_control%multiplicity - 1)
      IF (nSpin == 2) THEN
#if defined(__TBLITE_GXTB)
         ! g-xTB carries its fitted spin-polarization container directly on
         ! the calculator. Other methods still use CP2K's legacy wrapper.
         IF (.NOT. ALLOCATED(tb%calc%spin_polarization)) CALL tb_add_spin_polarization(tb)
#else
         CALL tb_add_spin_polarization(tb)
#endif
      END IF

      info = tb%calc%variable_info()
      IF (info%charge > shell_resolved) CPABORT("tblite: no support for orbital resolved charge")
      IF (info%dipole > atom_resolved) CPABORT("tblite: no support for shell resolved dipole moment")
      IF (info%quadrupole > atom_resolved) &
         CPABORT("tblite: no support shell resolved quadrupole moment")

      CALL new_wavefunction(tb%wfn, tb%mol%nat, tb%calc%bas%nsh, tb%calc%bas%nao, nSpin, 0.0_dp)
      CALL get_occupation(tb%mol, tb%calc%bas, tb%calc%h0, tb%wfn%nocc, tb%wfn%n0at, tb%wfn%n0sh)
#if defined(__TBLITE_GXTB)
      IF (dft_control%qs_control%xtb_control%tblite_method == gxtb) THEN
         IF (ALLOCATED(tb%mixer)) DEALLOCATE (tb%mixer)
      ELSE
         CALL tb_reset_mixer(tb)
      END IF
#else
      CALL tb_reset_mixer(tb)
#endif

      CALL new_potential(tb%pot, tb%mol, tb%calc%bas, tb%wfn%nspin)

      !allocate quantities later required
      ALLOCATE (tb%e_hal(tb%mol%nat), tb%e_rep(tb%mol%nat), tb%e_disp(tb%mol%nat))
      ALLOCATE (tb%e_scd(tb%mol%nat), tb%e_es(tb%mol%nat), tb%e_int(tb%mol%nat))
      ALLOCATE (tb%e_inc(tb%mol%nat), tb%e_exchange(tb%mol%nat), tb%e_spin(tb%mol%nat))
      ALLOCATE (tb%selfenergy(tb%calc%bas%nsh))
      IF (ALLOCATED(tb%calc%ncoord)) ALLOCATE (tb%cn(tb%mol%nat))

#else
      MARK_USED(tb)
      MARK_USED(dft_control)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_init_wf

#if defined(__TBLITE)
! **************************************************************************************************
!> \brief Add tblite's on-site spin-polarization interaction.
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_add_spin_polarization(tb)

      TYPE(tblite_type), POINTER                         :: tb

      CLASS(container_type), ALLOCATABLE                 :: cont
      TYPE(spin_polarization), ALLOCATABLE               :: spin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: wll

      ALLOCATE (spin)
      CALL tb_get_spin_constants(tb, wll)
      CALL new_spin_polarization(spin, tb%mol, wll, tb%calc%bas%nsh_id)
      CALL MOVE_ALLOC(spin, cont)
      CALL tb%calc%push_back(cont)

   END SUBROUTINE tb_add_spin_polarization

! **************************************************************************************************
!> \brief Build tblite spin constants for the current basis.
!> \param tb ...
!> \param wll ...
! **************************************************************************************************
   SUBROUTINE tb_get_spin_constants(tb, wll)

      TYPE(tblite_type), POINTER                         :: tb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: wll

      INTEGER                                            :: il, ish, izp, jl, jsh

      ALLOCATE (wll(tb%calc%bas%nsh, tb%calc%bas%nsh, tb%mol%nid))
      wll = 0.0_dp
      DO izp = 1, tb%mol%nid
         DO ish = 1, tb%calc%bas%nsh_id(izp)
#if defined(__TBLITE_GXTB)
            il = tb%calc%bas%cgto(ish, izp)%raw%ang
#else
            il = tb%calc%bas%cgto(ish, izp)%ang
#endif
            DO jsh = 1, tb%calc%bas%nsh_id(izp)
#if defined(__TBLITE_GXTB)
               jl = tb%calc%bas%cgto(jsh, izp)%raw%ang
#else
               jl = tb%calc%bas%cgto(jsh, izp)%ang
#endif
               wll(jsh, ish, izp) = get_spin_constant(jl, il, tb%mol%num(izp))
            END DO
         END DO
      END DO

   END SUBROUTINE tb_get_spin_constants
#endif

! **************************************************************************************************
!> \brief Reset tblite's internal SCC mixer for a new CP2K SCF cycle.
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_reset_mixer(tb)

      TYPE(tblite_type), POINTER                         :: tb

#if defined(__TBLITE)

      TYPE(scf_info)                                     :: info

      info = tb%calc%variable_info()
      IF (ALLOCATED(tb%mixer)) DEALLOCATE (tb%mixer)
      CALL new_cp2k_tblite_mixer(tb%mixer, tb%mixer_memory, &
                                 tb%wfn%nspin*get_mixer_dimension(tb%mol, tb%calc%bas, info), &
                                 tb%mixer_damping, tb%mixer_omega0, tb%mixer_min_weight, &
                                 tb%mixer_max_weight, tb%mixer_weight_factor)

#else
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_reset_mixer

! **************************************************************************************************
!> \brief Configure tblite's internal SCC mixer from CP2K input.
!> \param tb ...
!> \param iterations ...
!> \param memory ...
!> \param damping ...
!> \param omega0 ...
!> \param min_weight ...
!> \param max_weight ...
!> \param weight_factor ...
!> \param solver ...
! **************************************************************************************************
   SUBROUTINE tb_configure_mixer(tb, iterations, memory, damping, omega0, min_weight, max_weight, &
                                 weight_factor, solver)

      TYPE(tblite_type), POINTER                         :: tb
      INTEGER, INTENT(IN)                                :: iterations, memory, solver
      REAL(KIND=dp), INTENT(IN)                          :: damping, max_weight, min_weight, omega0, &
                                                            weight_factor

#if defined(__TBLITE)

      IF (iterations < 1) CPABORT("tblite SCC mixer ITERATIONS must be positive")
      IF (memory < 1) CPABORT("tblite SCC mixer MEMORY must be positive")
      IF (damping <= 0.0_dp) CPABORT("tblite SCC mixer damping must be positive")
      IF (omega0 <= 0.0_dp) CPABORT("tblite SCC mixer OMEGA0 must be positive")
      IF (min_weight <= 0.0_dp) CPABORT("tblite SCC mixer MIN_WEIGHT must be positive")
      IF (max_weight <= 0.0_dp) CPABORT("tblite SCC mixer MAX_WEIGHT must be positive")
      IF (max_weight < min_weight) &
         CPABORT("tblite SCC mixer MAX_WEIGHT must not be smaller than MIN_WEIGHT")
      IF (weight_factor <= 0.0_dp) CPABORT("tblite SCC mixer WEIGHT_FACTOR must be positive")
      SELECT CASE (solver)
      CASE (tblite_solver_gvd, tblite_solver_gvr)
      CASE DEFAULT
         CPABORT("Unknown tblite SCC mixer SOLVER")
      END SELECT

      tb%mixer_memory = memory
      tb%mixer_solver = solver
      tb%mixer_damping = damping
#if defined(__TBLITE_GXTB)
      CALL cp2k_set_scf_options(tb%calc, iterations, damping)
#else
      tb%calc%max_iter = iterations
      tb%calc%mixer_damping = damping
#endif
      tb%mixer_omega0 = omega0
      tb%mixer_min_weight = min_weight
      tb%mixer_max_weight = max_weight
      tb%mixer_weight_factor = weight_factor

#else
      MARK_USED(tb)
      MARK_USED(iterations)
      MARK_USED(memory)
      MARK_USED(damping)
      MARK_USED(omega0)
      MARK_USED(min_weight)
      MARK_USED(max_weight)
      MARK_USED(weight_factor)
      MARK_USED(solver)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_configure_mixer

! **************************************************************************************************
!> \brief Return whether the tblite native SCC mixer is active for this run.
!> \param dft_control ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_native_scc_mixer_active(dft_control) RESULT(use_native_mixer)

      TYPE(dft_control_type), POINTER                    :: dft_control
      LOGICAL                                            :: use_native_mixer

      use_native_mixer = .FALSE.
      IF (.NOT. ASSOCIATED(dft_control)) RETURN
      IF (dft_control%qs_control%do_ls_scf) RETURN

      SELECT CASE (dft_control%qs_control%xtb_control%tblite_scc_mixer)
      CASE (tblite_scc_mixer_auto)
         use_native_mixer = .TRUE.
      CASE (tblite_scc_mixer_tblite)
         use_native_mixer = .TRUE.
      CASE (tblite_scc_mixer_cp2k, tblite_scc_mixer_none)
         use_native_mixer = .FALSE.
      CASE DEFAULT
         CPABORT("Unknown tblite SCC mixer")
      END SELECT

   END FUNCTION tb_native_scc_mixer_active

! **************************************************************************************************
!> \brief Return the native tblite SCC mixer residual on the CP2K iter_delta scale.
!> \param dft_control ...
!> \param tb ...
!> \param eps_scf CP2K reporting scale for the native residual.
!> \return ...
! **************************************************************************************************
   FUNCTION tb_scf_mixer_error(dft_control, tb, eps_scf) RESULT(mixer_error)

      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(tblite_type), POINTER                         :: tb
      REAL(KIND=dp), INTENT(IN)                          :: eps_scf
      REAL(KIND=dp)                                      :: mixer_error

#if defined(__TBLITE)
      REAL(KIND=dp)                                      :: raw_error
#endif

      mixer_error = 0.0_dp

#if defined(__TBLITE)
      IF (.NOT. ASSOCIATED(tb)) RETURN
      IF (.NOT. tb_native_scc_mixer_active(dft_control)) RETURN
      raw_error = 0.0_dp
#if defined(__TBLITE_GXTB)
      IF (dft_control%qs_control%xtb_control%tblite_method == gxtb) THEN
         IF (ALLOCATED(tb%potential_mixer)) &
            raw_error = REAL(tb%potential_mixer%get_error(), KIND=dp)
      ELSEIF (ALLOCATED(tb%mixer)) THEN
         raw_error = REAL(tb%mixer%get_error(), KIND=dp)
      END IF
#else
      IF (ALLOCATED(tb%mixer)) raw_error = REAL(tb%mixer%get_error(), KIND=dp)
#endif
      IF (raw_error > 0.0_dp) THEN
         IF (eps_scf > 0.0_dp) THEN
            mixer_error = eps_scf*raw_error/ &
                          (tblite_scc_pconv*dft_control%qs_control%xtb_control%tblite_accuracy)
         ELSE
            mixer_error = raw_error
         END IF
      END IF
#else
      MARK_USED(dft_control)
      MARK_USED(tb)
      MARK_USED(eps_scf)
#endif

   END FUNCTION tb_scf_mixer_error

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param typ ...
!> \param accuracy ...
!> \param param_file ...
! **************************************************************************************************
   SUBROUTINE tb_set_calculator(tb, typ, accuracy, param_file)

      TYPE(tblite_type), POINTER                         :: tb
      INTEGER                                            :: typ
      REAL(KIND=dp), INTENT(IN)                          :: accuracy
      CHARACTER(LEN=*), INTENT(IN)                       :: param_file

#if defined(__TBLITE)

      TYPE(error_type), ALLOCATABLE                      :: error

      IF (ALLOCATED(tb%param)) DEALLOCATE (tb%param)
      IF (LEN_TRIM(param_file) > 0) THEN
         ALLOCATE (tb%param)
         CALL tb%param%load(TRIM(param_file), error)
         IF (ALLOCATED(error)) CPABORT("Could not load tblite PARAM file: "//TRIM(param_file))
         CALL new_xtb_calculator(tb%calc, tb%mol, tb%param, error)
      ELSE
         SELECT CASE (typ)
         CASE default
            CPABORT("Unknown xtb type")
         CASE (gfn1xtb)
            CALL new_gfn1_calculator(tb%calc, tb%mol, error)
         CASE (gfn2xtb)
            CALL new_gfn2_calculator(tb%calc, tb%mol, error)
         CASE (ipea1xtb)
            CALL new_ipea1_calculator(tb%calc, tb%mol, error)
         CASE (gxtb)
#if defined(__TBLITE_GXTB)
            CALL new_gxtb_calculator(tb%calc, tb%mol, error)
#else
            CPABORT("METHOD GXTB requires CP2K_TBLITE_PROVIDER=SAVE")
#endif
         END SELECT
      END IF
      IF (ALLOCATED(error)) CPABORT("tblite calculator setup failed")

      tb%accuracy = accuracy
#if defined(__TBLITE_GXTB)
      CALL cp2k_update_basis(tb%calc, tb%mol, tb%bcache, error)
      IF (ALLOCATED(error)) CPABORT("save_tblite basis setup failed")
#endif

#else
      MARK_USED(tb)
      MARK_USED(typ)
      MARK_USED(accuracy)
      MARK_USED(param_file)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_set_calculator

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_init_ham(qs_env, tb, para_env)

      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tblite_type)                                  :: tb
      TYPE(mp_para_env_type)                             :: para_env

#if defined(__TBLITE)

      TYPE(container_cache)                              :: hcache, incache, rcache

      tb%e_hal = 0.0_dp
      tb%e_rep = 0.0_dp
      tb%e_disp = 0.0_dp
      tb%e_int = 0.0_dp
      tb%e_inc = 0.0_dp
      tb%e_exchange = 0.0_dp
      tb%e_spin = 0.0_dp
      IF (ALLOCATED(tb%grad)) THEN
         tb%grad = 0.0_dp
         CALL tb_zero_force(qs_env)
      END IF
      tb%sigma = 0.0_dp

#if defined(__TBLITE_GXTB)
      IF (ALLOCATED(tb%calc%increment)) THEN
         CALL tb%calc%increment%get_engrad(tb%mol, incache, tb%e_inc)
      END IF
#endif

      IF (ALLOCATED(tb%calc%halogen)) THEN
         CALL tb%calc%halogen%update(tb%mol, hcache)
         IF (ALLOCATED(tb%grad)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%halogen%get_engrad(tb%mol, hcache, tb%e_hal, &
            & tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_halogen", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 0)
         ELSE
            CALL tb%calc%halogen%get_engrad(tb%mol, hcache, tb%e_hal)
         END IF
      END IF

      IF (ALLOCATED(tb%calc%repulsion)) THEN
         CALL tb%calc%repulsion%update(tb%mol, rcache)
         IF (ALLOCATED(tb%grad)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%repulsion%get_engrad(tb%mol, rcache, tb%e_rep, &
            & tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_repulsion", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 1)
         ELSE
            CALL tb%calc%repulsion%get_engrad(tb%mol, rcache, tb%e_rep)
         END IF
      END IF

      IF (ALLOCATED(tb%calc%dispersion)) THEN
         CALL tb%calc%dispersion%update(tb%mol, tb%dcache)
         IF (ALLOCATED(tb%grad)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%dispersion%get_engrad(tb%mol, tb%dcache, tb%e_disp, &
            & tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_dispersion_static", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 2)
         ELSE
            CALL tb%calc%dispersion%get_engrad(tb%mol, tb%dcache, tb%e_disp)
         END IF
      END IF

      IF (ALLOCATED(tb%calc%interactions)) THEN
         CALL tb%calc%interactions%update(tb%mol, tb%icache)
      END IF

      CALL new_potential(tb%pot, tb%mol, tb%calc%bas, tb%wfn%nspin)
      IF (ALLOCATED(tb%calc%coulomb)) THEN
         CALL tb%calc%coulomb%update(tb%mol, tb%cache)
      END IF

      IF (ALLOCATED(tb%grad)) THEN
         IF (ALLOCATED(tb%calc%ncoord)) THEN
            CALL tb%calc%ncoord%get_cn(tb%mol, tb%cn, tb%dcndr, tb%dcndL)
         END IF
         CALL get_selfenergy(tb%calc%h0, tb%mol%id, tb%calc%bas%ish_at, &
         &  tb%calc%bas%nsh_id, cn=tb%cn, selfenergy=tb%selfenergy, dsedcn=tb%dsedcn)
      ELSE
         IF (ALLOCATED(tb%calc%ncoord)) THEN
            CALL tb%calc%ncoord%get_cn(tb%mol, tb%cn)
         END IF
         CALL get_selfenergy(tb%calc%h0, tb%mol%id, tb%calc%bas%ish_at, &
         &  tb%calc%bas%nsh_id, cn=tb%cn, selfenergy=tb%selfenergy, dsedcn=tb%dsedcn)
      END IF

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(para_env)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_init_ham

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE tb_get_energy(qs_env, tb, energy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      TYPE(qs_energy_type), POINTER                      :: energy

#if defined(__TBLITE)

      INTEGER                                            :: iounit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: scf_section
      REAL(KIND=dp)                                      :: xtb_inter
      NULLIFY (scf_section, logger)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      scf_section => section_vals_get_subs_vals(qs_env%input, "DFT%SCF")

      energy%repulsive = SUM(tb%e_rep)
      energy%el_stat = SUM(tb%e_es)
      energy%dispersion = SUM(tb%e_disp)
      energy%dispersion_sc = SUM(tb%e_scd)
      energy%xtb_xb_inter = SUM(tb%e_hal)
      xtb_inter = SUM(tb%e_int) + SUM(tb%e_inc) + SUM(tb%e_exchange) + SUM(tb%e_spin)
#if defined(__TBLITE_GXTB)
      xtb_inter = xtb_inter + tb%e_acp_kp
#endif

      energy%total = energy%core + energy%repulsive + energy%el_stat + energy%dispersion &
                     + energy%dispersion_sc + energy%xtb_xb_inter + xtb_inter &
                     + energy%kTS + energy%efield + energy%qmmm_el

      iounit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DETAILED_ENERGY", &
                                    extension=".scfLog")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,(T9,A,T60,F20.10))") &
            "Repulsive pair potential energy:               ", energy%repulsive, &
            "Zeroth order Hamiltonian energy:               ", energy%core, &
            "Electrostatic energy:                          ", energy%el_stat, &
            "Self-consistent dispersion energy:             ", energy%dispersion_sc, &
            "Non-self consistent dispersion energy:         ", energy%dispersion
         IF (ABS(energy%xtb_xb_inter) > 1.e-9_dp) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "Correction for halogen bonding:                ", energy%xtb_xb_inter
         END IF
         IF (ABS(xtb_inter) > 1.e-9_dp) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "Additional interaction (e.g. spin):            ", xtb_inter
         END IF
         IF (ABS(energy%efield) > 1.e-9_dp) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "Electric field interaction energy:             ", energy%efield
         END IF
         IF (qs_env%qmmm) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "QM/MM Electrostatic energy:                    ", energy%qmmm_el
         END IF
      END IF
      CALL cp_print_key_finished_output(iounit, logger, scf_section, &
                                        "PRINT%DETAILED_ENERGY")

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(energy)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_get_energy

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param gto_basis_set ...
!> \param element_symbol ...
!> \param param ...
!> \param occ ...
! **************************************************************************************************
   SUBROUTINE tb_get_basis(tb, gto_basis_set, element_symbol, param, occ)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(gto_basis_set_type), POINTER                  :: gto_basis_set
      CHARACTER(len=2), INTENT(IN)                       :: element_symbol
      TYPE(xtb_atom_type), POINTER                       :: param
      INTEGER, DIMENSION(5), INTENT(out)                 :: occ

#if defined(__TBLITE)

      REAL(KIND=dp)                                      :: docc
      CHARACTER(LEN=default_string_length)               :: sng
      INTEGER                                            :: ang, i_type, id_atom, ind_ao, ipgf, iSH, &
                                                            ishell, ityp, maxl, mprim, natorb, &
                                                            nset, nshell
      LOGICAL                                            :: do_ortho
#if defined(__TBLITE_GXTB)
      REAL(KIND=dp), ALLOCATABLE                         :: tb_alpha(:), tb_coeff(:)
#endif

      CALL allocate_gto_basis_set(gto_basis_set)

      !identifying element in the bas data
      CALL symbol_to_number(i_type, element_symbol)
      DO id_atom = 1, tb%mol%nat
         IF (i_type == tb%el_num(id_atom)) EXIT
      END DO
      param%z = i_type
      param%symbol = element_symbol
      param%defined = .TRUE.
      ityp = tb%mol%id(id_atom)

      !getting size information
      nset = tb%calc%bas%nsh_id(ityp)
      nshell = 1
      mprim = 0
      DO ishell = 1, nset
#if defined(__TBLITE_GXTB)
         mprim = MAX(mprim, tb%calc%bas%cgto(ishell, ityp)%raw%nprim)
#else
         mprim = MAX(mprim, tb%calc%bas%cgto(ishell, ityp)%nprim)
#endif
      END DO
      param%nshell = nset
      natorb = 0

      !write basis set information
      CALL integer_to_string(mprim, sng)
      gto_basis_set%name = element_symbol//"_STO-"//TRIM(sng)//"G"
      gto_basis_set%nset = nset
      CALL reallocate(gto_basis_set%lmax, 1, nset)
      CALL reallocate(gto_basis_set%lmin, 1, nset)
      CALL reallocate(gto_basis_set%npgf, 1, nset)
      CALL reallocate(gto_basis_set%nshell, 1, nset)
      CALL reallocate(gto_basis_set%n, 1, 1, 1, nset)
      CALL reallocate(gto_basis_set%l, 1, 1, 1, nset)
      CALL reallocate(gto_basis_set%zet, 1, mprim, 1, nset)
      CALL reallocate(gto_basis_set%gcc, 1, mprim, 1, 1, 1, nset)

      ind_ao = 0
      maxl = 0
#if defined(__TBLITE_GXTB)
      ALLOCATE (tb_alpha(mprim), tb_coeff(mprim))
#endif
      DO ishell = 1, nset
#if defined(__TBLITE_GXTB)
         CALL cp2k_get_cgto(tb%calc%bas, tb%bcache, ishell, ityp, id_atom, &
                            ang, gto_basis_set%npgf(ishell), tb_alpha, tb_coeff)
#else
         ang = tb%calc%bas%cgto(ishell, ityp)%ang
#endif
         natorb = natorb + (2*ang + 1)
         param%lval(ishell) = ang
         maxl = MAX(ang, maxl)
         gto_basis_set%lmax(ishell) = ang
         gto_basis_set%lmin(ishell) = ang
#if !defined(__TBLITE_GXTB)
         gto_basis_set%npgf(ishell) = tb%calc%bas%cgto(ishell, ityp)%nprim
#endif
         gto_basis_set%nshell(ishell) = nshell
         gto_basis_set%n(1, ishell) = ang + 1
         gto_basis_set%l(1, ishell) = ang
         DO ipgf = 1, gto_basis_set%npgf(ishell)
#if defined(__TBLITE_GXTB)
            gto_basis_set%gcc(ipgf, 1, ishell) = tb_coeff(ipgf)
            gto_basis_set%zet(ipgf, ishell) = tb_alpha(ipgf)
#else
            gto_basis_set%gcc(ipgf, 1, ishell) = tb%calc%bas%cgto(ishell, ityp)%coeff(ipgf)
            gto_basis_set%zet(ipgf, ishell) = tb%calc%bas%cgto(ishell, ityp)%alpha(ipgf)
#endif
         END DO
         DO ipgf = 1, (2*ang + 1)
            ind_ao = ind_ao + 1
            param%lao(ind_ao) = ang
            param%nao(ind_ao) = ishell
         END DO
      END DO

      do_ortho = .FALSE.
      CALL process_gto_basis(gto_basis_set, do_ortho, nset, maxl)

      !setting additional values in parameter
#if defined(__TBLITE_GXTB)
      param%rcut = cp2k_get_acp_cutoff(tb%calc, tb%accuracy)
#else
      param%rcut = get_cutoff(tb%calc%bas, tb%accuracy)
#endif
      param%natorb = natorb
      param%lmax = maxl !max angular momentum

      !getting occupation
      occ = 0
      docc = 0.0_dp
      IF (tb%calc%bas%nsh_at(id_atom) > 5) CPABORT("too many shells in tblite")
      DO iSh = 1, tb%calc%bas%nsh_at(id_atom)
         occ(iSh) = NINT(tb%calc%h0%refocc(iSh, ityp) + docc)
         docc = docc + tb%calc%h0%refocc(iSh, ityp) - REAL(occ(iSh))
         param%occupation(iSh) = occ(iSh)
      END DO
      IF (ABS(docc) > 0.1_dp) CPABORT("Getting occupation numbers from tblite fails")
      param%zeff = SUM(occ) !effective core charge

      !set normalization process
      gto_basis_set%norm_type = 3

#else
      occ = 0
      MARK_USED(tb)
      MARK_USED(gto_basis_set)
      MARK_USED(element_symbol)
      MARK_USED(param)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_get_basis

#if defined(__TBLITE_GXTB)
! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param iatom ...
!> \param first_ao ...
!> \param nao ...
!> Return the dense tblite AO range belonging to one CP2K atom block.
! **************************************************************************************************
   SUBROUTINE tb_atom_ao_range(tb, iatom, first_ao, nao)

      TYPE(tblite_type), POINTER                         :: tb
      INTEGER, INTENT(IN)                                :: iatom
      INTEGER, INTENT(OUT)                               :: first_ao, nao

      INTEGER                                            :: first_shell, last_shell

      first_shell = tb%calc%bas%ish_at(iatom) + 1
      last_shell = first_shell + tb%calc%bas%nsh_at(iatom) - 1
      first_ao = tb%calc%bas%iao_sh(first_shell) + 1
      nao = SUM(tb%calc%bas%nao_sh(first_shell:last_shell))

   END SUBROUTINE tb_atom_ao_range

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param dense ...
!> \param matrix ...
!> Scatter a dense tblite Gamma matrix into locally owned CP2K atom blocks.
! **************************************************************************************************
   SUBROUTINE tb_dense_to_dbcsr(tb, dense, matrix)

      TYPE(tblite_type), POINTER                         :: tb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: dense
      TYPE(dbcsr_type), POINTER                          :: matrix

      INTEGER                                            :: col, col_first, col_nao, row, row_first, &
                                                            row_nao
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         NULLIFY (block)
         CALL dbcsr_iterator_next_block(iter, row, col, block)
         CALL tb_atom_ao_range(tb, row, row_first, row_nao)
         CALL tb_atom_ao_range(tb, col, col_first, col_nao)
         CPASSERT(SIZE(block, 1) == row_nao)
         CPASSERT(SIZE(block, 2) == col_nao)
         block(:, :) = dense(row_first:row_first + row_nao - 1, &
                             col_first:col_first + col_nao - 1)
      END DO
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE tb_dense_to_dbcsr

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param matrix ...
!> \param dense ...
!> \param para_env ...
!> Gather a symmetric CP2K Gamma matrix into tblite atom/AO ordering.
! **************************************************************************************************
   SUBROUTINE tb_dbcsr_to_dense(tb, matrix, dense, para_env)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(dbcsr_type), POINTER                          :: matrix
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: dense
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: col, col_first, col_nao, row, row_first, &
                                                            row_nao
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      dense = 0.0_dp
      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         NULLIFY (block)
         CALL dbcsr_iterator_next_block(iter, row, col, block)
         CALL tb_atom_ao_range(tb, row, row_first, row_nao)
         CALL tb_atom_ao_range(tb, col, col_first, col_nao)
         CPASSERT(SIZE(block, 1) == row_nao)
         CPASSERT(SIZE(block, 2) == col_nao)
         dense(row_first:row_first + row_nao - 1, &
               col_first:col_first + col_nao - 1) = block
         IF (row /= col) THEN
            dense(col_first:col_first + col_nao - 1, &
                  row_first:row_first + row_nao - 1) = TRANSPOSE(block)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL para_env%sum(dense)

   END SUBROUTINE tb_dbcsr_to_dense

! **************************************************************************************************
!> \brief ...
!> \param index_to_cell ...
!> \param cell_to_index ...
!> \param closed_index_to_cell ...
!> \param closed_cell_to_index ...
!> \param nimg_closed ...
!> Close CP2K's compact k-point image set under lattice inversion for save_tblite.
!>
!> Original cells retain their CP2K image indices.  Missing inverse cells are appended in original
!> order, yielding a deterministic response/build view on every MPI rank without changing CP2K's
!> DBCSR storage.  The returned dense lookup preserves the possibly negative bounds of CP2K's map.
! **************************************************************************************************
   SUBROUTINE tb_close_gxtb_image_cells(index_to_cell, cell_to_index, &
                                        closed_index_to_cell, closed_cell_to_index, nimg_closed)

      INTEGER, DIMENSION(:, :), INTENT(IN), POINTER      :: index_to_cell
      INTEGER, DIMENSION(:, :, :), INTENT(IN), POINTER   :: cell_to_index
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: closed_index_to_cell
      INTEGER, DIMENSION(:, :, :), INTENT(OUT), POINTER  :: closed_cell_to_index
      INTEGER, INTENT(OUT)                               :: nimg_closed

      INTEGER                                            :: ic, img, inv, nimg
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: image_work
      INTEGER, DIMENSION(3)                              :: image_cell
      LOGICAL                                            :: image_in_bounds

      nimg = SIZE(index_to_cell, 2)
      CPASSERT(nimg > 0)
      ALLOCATE (image_work(3, 2*nimg))
      image_work(:, 1:nimg) = index_to_cell
      nimg_closed = nimg
      DO img = 1, nimg
         image_cell = index_to_cell(:, img)
         image_in_bounds = &
            cell_to_index(image_cell(1), image_cell(2), image_cell(3)) == img
         CPASSERT(image_in_bounds)
         image_cell = -image_cell
         IF (cell_to_index(image_cell(1), image_cell(2), image_cell(3)) == 0) THEN
            nimg_closed = nimg_closed + 1
            CPASSERT(nimg_closed <= SIZE(image_work, 2))
            image_work(:, nimg_closed) = image_cell
         END IF
      END DO
      ALLOCATE (closed_index_to_cell(3, nimg_closed))
      closed_index_to_cell(:, :) = image_work(:, 1:nimg_closed)
      DEALLOCATE (image_work)

      ALLOCATE (closed_cell_to_index( &
                LBOUND(cell_to_index, 1):UBOUND(cell_to_index, 1), &
                LBOUND(cell_to_index, 2):UBOUND(cell_to_index, 2), &
                LBOUND(cell_to_index, 3):UBOUND(cell_to_index, 3)))
      closed_cell_to_index = 0
      DO ic = 1, nimg_closed
         image_cell = closed_index_to_cell(:, ic)
         image_in_bounds = &
            image_cell(1) >= LBOUND(closed_cell_to_index, 1) .AND. &
            image_cell(1) <= UBOUND(closed_cell_to_index, 1) .AND. &
            image_cell(2) >= LBOUND(closed_cell_to_index, 2) .AND. &
            image_cell(2) <= UBOUND(closed_cell_to_index, 2) .AND. &
            image_cell(3) >= LBOUND(closed_cell_to_index, 3) .AND. &
            image_cell(3) <= UBOUND(closed_cell_to_index, 3)
         CPASSERT(image_in_bounds)
         closed_cell_to_index(image_cell(1), image_cell(2), image_cell(3)) = ic
      END DO
      CPASSERT(MAXVAL(closed_cell_to_index) == nimg_closed)
      DO ic = 1, nimg_closed
         image_cell = closed_index_to_cell(:, ic)
         image_in_bounds = &
            closed_cell_to_index(image_cell(1), image_cell(2), image_cell(3)) == ic
         CPASSERT(image_in_bounds)
         inv = closed_cell_to_index(-image_cell(1), -image_cell(2), -image_cell(3))
         image_in_bounds = inv > 0 .AND. inv <= nimg_closed
         CPASSERT(image_in_bounds)
         image_in_bounds = ALL(closed_index_to_cell(:, inv) == -image_cell)
         CPASSERT(image_in_bounds)
      END DO

   END SUBROUTINE tb_close_gxtb_image_cells

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param matrices ...
!> \param sab_kp ...
!> \param cell_to_index ...
!> \param dense_cell_to_index ...
!> \param dense ...
!> \param para_env ...
!> Gather CP2K real-space image matrices in save_tblite's AO/image orientation.
!>
!> This restores the complete dual representation of `tb_scatter_gxtb_images`.  An individual
!> CP2K image matrix is not a symmetric dense matrix: its upper-triangular atom blocks encode the atom/cell
!> canonicalization of the k-point neighbour list.  Expanding it with `tb_dbcsr_to_dense` would
!> therefore put the transpose into the wrong image.  The returned arrays instead obey the dual
!> convention used by save_tblite, E=sum_R,ij P(i,j,R)*H(i,j,R).  For different atoms CP2K's
!> symmetric neighbour list retains only one member of the (i,j,R)/(j,i,-R) pair, so the missing
!> inverse-cell transpose is reconstructed explicitly.  Same-atom +/-R entries are already
!> present explicitly and need no reconstruction.
! **************************************************************************************************
   SUBROUTINE tb_gather_gxtb_images(tb, matrices, sab_kp, cell_to_index, &
                                    dense_cell_to_index, dense, para_env)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrices
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_kp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index, dense_cell_to_index
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(OUT)  :: dense
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: col_first, col_nao, iatom, ic, ic_dense, &
                                                            ic_inv, icol, irow, ispin, jatom, &
                                                            row_first, row_nao
      INTEGER, DIMENSION(3)                              :: cell
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      CPASSERT(SIZE(dense, 1) == tb%calc%bas%nao)
      CPASSERT(SIZE(dense, 2) == tb%calc%bas%nao)
      found = MAXVAL(dense_cell_to_index) == SIZE(dense, 3)
      CPASSERT(found)
      CPASSERT(SIZE(dense, 4) == SIZE(matrices, 1))

      dense = 0.0_dp
      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, sab_kp)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, cell=cell)
         ic = cell_to_index(cell(1), cell(2), cell(3))
         CPASSERT(ic > 0 .AND. ic <= SIZE(matrices, 2))
         ic_dense = dense_cell_to_index(cell(1), cell(2), cell(3))
         CPASSERT(ic_dense > 0 .AND. ic_dense <= SIZE(dense, 3))
         IF (iatom /= jatom) THEN
            ic_inv = dense_cell_to_index(-cell(1), -cell(2), -cell(3))
            CPASSERT(ic_inv > 0 .AND. ic_inv <= SIZE(dense, 3))
         ELSE
            ic_inv = cell_to_index(-cell(1), -cell(2), -cell(3))
            CPASSERT(ic_inv > 0 .AND. ic_inv <= SIZE(matrices, 2))
         END IF
         irow = MIN(iatom, jatom)
         icol = MAX(iatom, jatom)
         CALL tb_atom_ao_range(tb, iatom, row_first, row_nao)
         CALL tb_atom_ao_range(tb, jatom, col_first, col_nao)

         DO ispin = 1, SIZE(matrices, 1)
            NULLIFY (block)
            CALL dbcsr_get_block_p(matrix=matrices(ispin, ic)%matrix, row=irow, col=icol, &
                                   BLOCK=block, found=found)
            CPASSERT(found)
            IF (iatom <= jatom) THEN
               CPASSERT(SIZE(block, 1) == row_nao .AND. SIZE(block, 2) == col_nao)
               dense(col_first:col_first + col_nao - 1, &
                     row_first:row_first + row_nao - 1, ic_dense, ispin) = TRANSPOSE(block)
               IF (iatom /= jatom) &
                  dense(row_first:row_first + row_nao - 1, &
                        col_first:col_first + col_nao - 1, ic_inv, ispin) = block
            ELSE
               CPASSERT(SIZE(block, 1) == col_nao .AND. SIZE(block, 2) == row_nao)
               dense(col_first:col_first + col_nao - 1, &
                     row_first:row_first + row_nao - 1, ic_dense, ispin) = block
               IF (iatom /= jatom) &
                  dense(row_first:row_first + row_nao - 1, &
                        col_first:col_first + col_nao - 1, ic_inv, ispin) = TRANSPOSE(block)
            END IF
         END DO
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)
      CALL para_env%sum(dense)

   END SUBROUTINE tb_gather_gxtb_images

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param dense ...
!> \param para_env ...
!> Gather one distributed full matrix into a dense array replicated on the global communicator.
! **************************************************************************************************
   SUBROUTINE tb_fm_to_dense(matrix, dense, para_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: matrix
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: dense
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: iloc, jloc, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data

      NULLIFY (col_indices, row_indices, local_data)
      CALL cp_fm_get_info(matrix, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, local_data=local_data)
      dense = 0.0_dp
      DO jloc = 1, ncol_local
         DO iloc = 1, nrow_local
            dense(row_indices(iloc), col_indices(jloc)) = local_data(iloc, jloc)
         END DO
      END DO
      CALL para_env%sum(dense)

   END SUBROUTINE tb_fm_to_dense

#if defined(__TBLITE_GXTB)
! **************************************************************************************************
!> \brief ...
!> \param dense ...
!> \param matrix ...
!> Scatter a replicated dense matrix into one distributed full matrix.
! **************************************************************************************************
   SUBROUTINE tb_dense_to_fm(dense, matrix)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: dense
      TYPE(cp_fm_type), INTENT(INOUT)                    :: matrix

      INTEGER                                            :: iloc, jloc, ncol_global, ncol_local, &
                                                            nrow_global, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data

      NULLIFY (col_indices, row_indices, local_data)
      CALL cp_fm_get_info(matrix, nrow_global=nrow_global, ncol_global=ncol_global, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, local_data=local_data)
      CPASSERT(nrow_global == SIZE(dense, 1))
      CPASSERT(ncol_global == SIZE(dense, 2))
      DO jloc = 1, ncol_local
         DO iloc = 1, nrow_local
            local_data(iloc, jloc) = dense(row_indices(iloc), col_indices(jloc))
         END DO
      END DO

   END SUBROUTINE tb_dense_to_fm

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param raw_fock ...
!> \param mixed_fock ...
!> \param iter_count ...
!> Advance save_tblite's native g-xTB potential mixer on one replicated real vector.
! **************************************************************************************************
   SUBROUTINE tb_advance_gxtb_potential_mixer(tb, raw_fock, mixed_fock, iter_count)

      TYPE(tblite_type), POINTER                         :: tb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: raw_fock
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: mixed_fock
      INTEGER, INTENT(IN)                                :: iter_count

      TYPE(error_type), ALLOCATABLE                      :: error

      IF (iter_count < 1) CPABORT("Invalid SCF iteration for the g-xTB potential mixer")
      IF (SIZE(raw_fock) < 1 .OR. SIZE(mixed_fock) /= SIZE(raw_fock)) &
         CPABORT("Invalid g-xTB potential mixer vector dimensions")
      IF (.NOT. ALL(ieee_is_finite(raw_fock))) &
         CPABORT("CP2K supplied a non-finite Fock matrix to the g-xTB potential mixer")

      IF (.NOT. ALLOCATED(tb%potential_mixer)) THEN
         ALLOCATE (tb%potential_mixer)
         CALL new_cp2k_potential_mixer(tb%potential_mixer, SIZE(raw_fock), error)
         IF (ALLOCATED(error)) CPABORT("Could not initialize the g-xTB potential mixer")
         tb%potential_mixer_dimension = SIZE(raw_fock)
      ELSEIF (iter_count == 1 .OR. tb%potential_mixer_dimension /= SIZE(raw_fock)) THEN
         CALL tb%potential_mixer%reset(SIZE(raw_fock), error)
         IF (ALLOCATED(error)) CPABORT("Could not reset the g-xTB potential mixer")
         tb%potential_mixer_dimension = SIZE(raw_fock)
      END IF

      CALL tb%potential_mixer%advance(raw_fock, mixed_fock, error)
      IF (ALLOCATED(error)) CPABORT("The g-xTB potential mixer failed")
      IF (.NOT. ALL(ieee_is_finite(mixed_fock))) &
         CPABORT("save_tblite returned a non-finite mixed g-xTB Fock matrix")

   END SUBROUTINE tb_advance_gxtb_potential_mixer
#endif

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param matrix ...
!> \param dense ...
!> \param para_env ...
!> \param antisymmetric ...
!> Gather a symmetric or antisymmetric DBCSR k-point component into a replicated dense array.
! **************************************************************************************************
   SUBROUTINE tb_dbcsr_component_to_dense(tb, matrix, dense, para_env, antisymmetric)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(dbcsr_type), POINTER                          :: matrix
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: dense
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: antisymmetric

      INTEGER                                            :: col, col_first, col_nao, row, row_first, &
                                                            row_nao
      REAL(KIND=dp)                                      :: transpose_sign
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      transpose_sign = 1.0_dp
      IF (antisymmetric) transpose_sign = -1.0_dp
      dense = 0.0_dp
      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         NULLIFY (block)
         CALL dbcsr_iterator_next_block(iter, row, col, block)
         CALL tb_atom_ao_range(tb, row, row_first, row_nao)
         CALL tb_atom_ao_range(tb, col, col_first, col_nao)
         dense(row_first:row_first + row_nao - 1, &
               col_first:col_first + col_nao - 1) = block
         IF (row /= col) THEN
            dense(col_first:col_first + col_nao - 1, &
                  row_first:row_first + row_nao - 1) = transpose_sign*TRANSPOSE(block)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL para_env%sum(dense)

   END SUBROUTINE tb_dbcsr_component_to_dense

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> \param kpoints ...
!> \param density ...
!> \param target_index optional output slot for each irreducible k-point
!> Collect the unweighted complex density P(k) of every irreducible k-point.
!>
!> CP2K stores these matrices in the owning k-point group.  This follows the communication pattern
!> used by kpoint_density_transform, but gathers the matrices without applying w(k).  In particular,
!> it never reconstructs P(k) from the real-space density images, whose mesh-equivalent image cells
!> could otherwise introduce a multiplicity alias.
! **************************************************************************************************
   SUBROUTINE tb_collect_kpoint_density(qs_env, tb, kpoints, density, target_index)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      TYPE(kpoint_type), POINTER                         :: kpoints
      COMPLEX(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(OUT)                                     :: density
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: target_index

      INTEGER                                            :: ic, ik, ikk, ispin, itarget, nc, nkp, &
                                                            nspin
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: my_kpgrp, use_real_wfn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dense_part
      TYPE(copy_info_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: info
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_ao_fm_pools
      TYPE(cp_fm_type)                                   :: fmdummy
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fmwork
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (ao_ao_fm_pools, kp, para_env)
      CALL mpools_get(qs_env%mpools, ao_ao_fm_pools=ao_ao_fm_pools)
      CPASSERT(ASSOCIATED(ao_ao_fm_pools))
      CPASSERT(SIZE(ao_ao_fm_pools) >= 1)
      CPASSERT(ASSOCIATED(ao_ao_fm_pools(1)%pool))
      para_env => kpoints%blacs_env_all%para_env
      CALL get_kpoint_info(kpoints, nkp=nkp, kp_range=kp_range, use_real_wfn=use_real_wfn)
      nspin = tb%wfn%nspin
      nc = 1
      IF (.NOT. use_real_wfn) nc = 2
      CPASSERT(SIZE(density, 3) == nspin)
      IF (PRESENT(target_index)) THEN
         CPASSERT(SIZE(target_index) == nkp)
         CPASSERT(ALL(target_index >= 1 .AND. target_index <= SIZE(density, 4)))
      ELSE
         CPASSERT(SIZE(density, 4) == nkp)
      END IF
      ALLOCATE (fmwork(nc))
      DO ic = 1, nc
         CALL fm_pool_create_fm(ao_ao_fm_pools(1)%pool, fmwork(ic))
      END DO
      ALLOCATE (dense_part(SIZE(density, 1), SIZE(density, 2)))
      ALLOCATE (info(nc, nkp, nspin))

      DO ispin = 1, nspin
         DO ik = 1, nkp
            my_kpgrp = ik >= kp_range(1) .AND. ik <= kp_range(2)
            IF (my_kpgrp) THEN
               ikk = ik - kp_range(1) + 1
               kp => kpoints%kp_env(ikk)%kpoint_env
            ELSE
               NULLIFY (kp)
            END IF
            DO ic = 1, nc
               IF (my_kpgrp) THEN
                  CALL cp_fm_start_copy_general(kp%pmat(ic, ispin), fmwork(ic), para_env, info(ic, ik, ispin))
               ELSE
                  CALL cp_fm_start_copy_general(fmdummy, fmwork(ic), para_env, info(ic, ik, ispin))
               END IF
            END DO
         END DO
      END DO

      density = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO ispin = 1, nspin
         DO ik = 1, nkp
            itarget = ik
            IF (PRESENT(target_index)) itarget = target_index(ik)
            DO ic = 1, nc
               CALL cp_fm_finish_copy_general(fmwork(ic), info(ic, ik, ispin))
               CALL tb_fm_to_dense(fmwork(ic), dense_part, para_env)
               IF (ic == 1) THEN
                  density(:, :, ispin, itarget) = CMPLX(dense_part, 0.0_dp, KIND=dp)
               ELSE
                  density(:, :, ispin, itarget) = density(:, :, ispin, itarget) + &
                                                  CMPLX(0.0_dp, dense_part, KIND=dp)
               END IF
            END DO
         END DO
      END DO

      DO ispin = 1, nspin
         DO ik = 1, nkp
            my_kpgrp = ik >= kp_range(1) .AND. ik <= kp_range(2)
            IF (.NOT. my_kpgrp) CYCLE
            DO ic = 1, nc
               CALL cp_fm_cleanup_copy_general(info(ic, ik, ispin))
            END DO
         END DO
      END DO
      DO ic = 1, nc
         CALL fm_pool_give_back_fm(ao_ao_fm_pools(1)%pool, fmwork(ic))
      END DO
      DEALLOCATE (fmwork)
      DEALLOCATE (dense_part, info)

   END SUBROUTINE tb_collect_kpoint_density

! **************************************************************************************************
!> \brief Collect only the irreducible densities owned by this CP2K k-point group.
!> \param tb ...
!> \param kpoints ...
!> \param density_local local irreducible P(k), in increasing global-k order
!>
!> Unlike tb_collect_kpoint_density, this routine never enters the global
!> communicator and never materializes another k group's density.  Each pmat
!> is already distributed over para_env_kp; the official group communicator
!> therefore suffices to replicate one local dense block inside its owner
!> group.  Inter-group transfer remains blockwise in the caller.
! **************************************************************************************************
   SUBROUTINE tb_collect_local_kpoint_density(tb, kpoints, density_local)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(kpoint_type), POINTER                         :: kpoints
      COMPLEX(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(OUT)                                     :: density_local

      INTEGER                                            :: ic, ik_local, ispin, nc, nlocal, nspin
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: use_real_wfn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dense_part
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mp_para_env_type), POINTER                    :: para_env_kp

      NULLIFY (kp, para_env_kp)
      CALL get_kpoint_info(kpoints, kp_range=kp_range, use_real_wfn=use_real_wfn, &
                           para_env_kp=para_env_kp)
      CPASSERT(ASSOCIATED(para_env_kp))
      nlocal = kp_range(2) - kp_range(1) + 1
      nspin = tb%wfn%nspin
      nc = 1
      IF (.NOT. use_real_wfn) nc = 2
      CPASSERT(SIZE(density_local, 3) == nspin)
      CPASSERT(SIZE(density_local, 4) == nlocal)
      ALLOCATE (dense_part(SIZE(density_local, 1), SIZE(density_local, 2)))
      density_local = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      DO ik_local = 1, nlocal
         kp => kpoints%kp_env(ik_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(kp%nkpoint == kp_range(1) + ik_local - 1)
         DO ispin = 1, nspin
            DO ic = 1, nc
               CALL tb_fm_to_dense(kp%pmat(ic, ispin), dense_part, para_env_kp)
               IF (ic == 1) THEN
                  density_local(:, :, ispin, ik_local) = &
                     CMPLX(dense_part, 0.0_dp, KIND=dp)
               ELSE
                  density_local(:, :, ispin, ik_local) = &
                     density_local(:, :, ispin, ik_local) + &
                     CMPLX(0.0_dp, dense_part, KIND=dp)
               END IF
            END DO
         END DO
      END DO
      DEALLOCATE (dense_part)

   END SUBROUTINE tb_collect_local_kpoint_density

! **************************************************************************************************
!> \brief Drop geometry-dependent g-xTB k-point maps and dense overlap data.
!> \param tb ...
!>
!> The cache is retained across SCC iterations, but never across an actual coordinate or cell
!> change. Dynamic density and Fock data are deliberately not part of this static cache.
! **************************************************************************************************
   SUBROUTINE tb_reset_gxtb_kpoint_cache(tb)

      TYPE(tblite_type), INTENT(INOUT)                   :: tb

      IF (ALLOCATED(tb%gxtb_kmesh_map)) DEALLOCATE (tb%gxtb_kmesh_map)
      IF (ALLOCATED(tb%acp_fock_kp)) DEALLOCATE (tb%acp_fock_kp)
      IF (ALLOCATED(tb%qualification_acp_fock_full)) &
         DEALLOCATE (tb%qualification_acp_fock_full)
      tb%gxtb_density_map_validated = .FALSE.
      tb%gxtb_transform_validated = .FALSE.
      tb%acp_fock_kp_valid = .FALSE.
      tb%e_acp_kp = 0.0_dp

   END SUBROUTINE tb_reset_gxtb_kpoint_cache

! **************************************************************************************************
!> \brief Compare two fractional k-points modulo reciprocal lattice vectors.
!> \param xkp_a ...
!> \param xkp_b ...
!> \return true when both points are equivalent
! **************************************************************************************************
   LOGICAL FUNCTION tb_same_periodic_kpoint(xkp_a, xkp_b) RESULT(same)

      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp_a, xkp_b

      REAL(KIND=dp), DIMENSION(3)                        :: diff

      diff = xkp_a - xkp_b
      diff = diff - REAL(NINT(diff), KIND=dp)
      same = SUM(ABS(diff)) < 1.0E-8_dp

   END FUNCTION tb_same_periodic_kpoint

! **************************************************************************************************
!> \brief Find a fractional k-point in a complete periodic mesh.
!> \param xkp_mesh ...
!> \param xkp_search ...
!> \return matching mesh index, or zero
! **************************************************************************************************
   INTEGER FUNCTION tb_find_periodic_kpoint(xkp_mesh, xkp_search) RESULT(ik_match)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp_mesh
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp_search

      INTEGER                                            :: ik

      ik_match = 0
      DO ik = 1, SIZE(xkp_mesh, 2)
         IF (tb_same_periodic_kpoint(xkp_mesh(1:3, ik), xkp_search)) THEN
            ik_match = ik
            RETURN
         END IF
      END DO

   END FUNCTION tb_find_periodic_kpoint

! **************************************************************************************************
!> \brief Recover a complete uniform Cartesian mesh from explicit GENERAL k-points.
!> \param kpoints ...
!> \param nmesh number of distinct equally spaced points along each reciprocal axis
!>
!> save_tblite's coupled exchange needs the Born-von Karman mesh dimensions.  CP2K deliberately
!> leaves nkp_grid unset for GENERAL input, so accept that input only when its original, unreduced
!> list is a complete Cartesian product with uniform weights.  The caller retains the original
!> input ordering; this routine supplies only the validated mesh descriptor.
! **************************************************************************************************
   SUBROUTINE tb_infer_gxtb_general_mesh(kpoints, nmesh)

      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: nmesh

      INTEGER                                            :: iaxis, ik, inext, inode, nfull, slot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: occupied
      INTEGER, DIMENSION(3)                              :: node_index
      REAL(KIND=dp)                                      :: distance, expected_spacing, gap, q, &
                                                            swap, tol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: nodes

      IF (TRIM(kpoints%kp_scheme) /= "GENERAL" .OR. &
          .NOT. ASSOCIATED(kpoints%xkp_input) .OR. &
          .NOT. ASSOCIATED(kpoints%wkp_input)) &
         CPABORT("g-XTB GENERAL meshes require the original explicit k-point list")
      nfull = SIZE(kpoints%wkp_input)
      IF (nfull < 1 .OR. SIZE(kpoints%xkp_input, 1) < 3 .OR. &
          SIZE(kpoints%xkp_input, 2) /= nfull) &
         CPABORT("Incomplete explicit GENERAL k-point list for g-XTB exchange")
      IF (.NOT. ALL(ieee_is_finite(kpoints%xkp_input(1:3, 1:nfull))) .OR. &
          .NOT. ALL(ieee_is_finite(kpoints%wkp_input(1:nfull)))) &
         CPABORT("Non-finite explicit GENERAL k-point list for g-XTB exchange")

      ! This descriptor is passed directly to save_tblite's Fourier mesh, whose
      ! periodic indexing is defined at 1e-10.  Reject a coarser approximate
      ! GENERAL grid here instead of accepting it under CP2K's geometry tolerance.
      tol = 1.0E-10_dp
      IF (MAXVAL(ABS(kpoints%wkp_input(1:nfull) - &
                     1.0_dp/REAL(nfull, KIND=dp))) > tol) &
         CPABORT("g-XTB GENERAL meshes require uniform full-mesh weights")

      ALLOCATE (nodes(nfull, 3), source=0.0_dp)
      nmesh = 0
      DO iaxis = 1, 3
         DO ik = 1, nfull
            q = MODULO(kpoints%xkp_input(iaxis, ik), 1.0_dp)
            IF (ABS(q - 1.0_dp) < tol) q = 0.0_dp
            DO inode = 1, nmesh(iaxis)
               distance = ABS(q - nodes(inode, iaxis))
               distance = MIN(distance, ABS(distance - 1.0_dp))
               IF (distance < tol) EXIT
            END DO
            IF (inode > nmesh(iaxis)) THEN
               nmesh(iaxis) = nmesh(iaxis) + 1
               nodes(nmesh(iaxis), iaxis) = q
            END IF
         END DO

         ! Sort the distinct coordinates before checking all circular spacings.
         DO inode = 2, nmesh(iaxis)
            inext = inode
            DO WHILE (inext > 1 .AND. nodes(inext, iaxis) < nodes(inext - 1, iaxis))
               swap = nodes(inext, iaxis)
               nodes(inext, iaxis) = nodes(inext - 1, iaxis)
               nodes(inext - 1, iaxis) = swap
               inext = inext - 1
            END DO
         END DO
         expected_spacing = 1.0_dp/REAL(nmesh(iaxis), KIND=dp)
         DO inode = 1, nmesh(iaxis)
            inext = MOD(inode, nmesh(iaxis)) + 1
            gap = nodes(inext, iaxis) - nodes(inode, iaxis)
            IF (inext == 1) gap = gap + 1.0_dp
            IF (ABS(gap - expected_spacing) > tol) &
               CPABORT("g-XTB GENERAL k-points are not an equally spaced Cartesian mesh")
         END DO
      END DO
      IF (PRODUCT(nmesh) /= nfull) &
         CPABORT("g-XTB GENERAL k-points do not form a complete Cartesian mesh")

      ALLOCATE (occupied(nfull), source=0)
      DO ik = 1, nfull
         DO iaxis = 1, 3
            q = MODULO(kpoints%xkp_input(iaxis, ik), 1.0_dp)
            IF (ABS(q - 1.0_dp) < tol) q = 0.0_dp
            node_index(iaxis) = 0
            DO inode = 1, nmesh(iaxis)
               distance = ABS(q - nodes(inode, iaxis))
               distance = MIN(distance, ABS(distance - 1.0_dp))
               IF (distance < tol) THEN
                  node_index(iaxis) = inode
                  EXIT
               END IF
            END DO
            IF (node_index(iaxis) == 0) &
               CPABORT("Failed to index an explicit g-XTB GENERAL k-point")
         END DO
         slot = node_index(1) + nmesh(1)*((node_index(2) - 1) + &
                                          nmesh(2)*(node_index(3) - 1))
         occupied(slot) = occupied(slot) + 1
      END DO
      IF (ANY(occupied /= 1)) &
         CPABORT("g-XTB GENERAL k-points contain duplicate or missing Cartesian mesh points")
      DEALLOCATE (nodes, occupied)

   END SUBROUTINE tb_infer_gxtb_general_mesh

! **************************************************************************************************
!> \brief Build the complete regular mesh represented by CP2K's irreducible k-points.
!> \param kpoints ...
!> \param kmap ...
!>
!> symop=0 is a direct representative, symop=-1 is pure time reversal, and a positive value is
!> the corresponding entry of source(k)%kp_sym.  Direct and pure-time-reversal matches take
!> precedence over atom/AO operations, following the full-mesh Wannier90 reconstruction.
! **************************************************************************************************
   SUBROUTINE tb_build_gxtb_kmesh_map(kpoints, kmap)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), INTENT(OUT)          :: kmap

      INTEGER                                            :: expected_ops, ifull, ikred, imatch, iop, &
                                                            ir, ira, isym, jr, max_ops, nred
      LOGICAL                                            :: explicit_general_mesh
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: reduced_weight
      TYPE(kpoint_sym_type), POINTER                     :: kpsym

      NULLIFY (kpsym)
      nred = kpoints%nkp
      kmap%scheme = TRIM(kpoints%kp_scheme)
      kmap%nred = nred
      kmap%full_grid = kpoints%full_grid
      kmap%inversion_only = kpoints%inversion_symmetry_only
      kmap%symmetry_backend = kpoints%symmetry_backend
      kmap%symmetry_method = kpoints%symmetry_reduction_method
      kmap%kp_shift = kpoints%kp_shift
      kmap%gamma_centered = kpoints%gamma_centered
      IF (nred < 1 .OR. .NOT. ASSOCIATED(kpoints%xkp) .OR. &
          .NOT. ASSOCIATED(kpoints%wkp)) &
         CPABORT("Incomplete CP2K k-point set for g-xTB exchange")
      IF (ANY(kpoints%wkp(1:nred) <= 0.0_dp) .OR. &
          ABS(SUM(kpoints%wkp(1:nred)) - 1.0_dp) > 1.0E-10_dp) &
         CPABORT("Invalid irreducible CP2K k-point weights for g-xTB exchange")

      explicit_general_mesh = TRIM(kpoints%kp_scheme) == "GENERAL" .AND. &
                              ASSOCIATED(kpoints%xkp_input) .AND. &
                              ASSOCIATED(kpoints%wkp_input)
      kmap%explicit_general = explicit_general_mesh
      IF (explicit_general_mesh) THEN
         CALL tb_infer_gxtb_general_mesh(kpoints, kmap%nmesh)
      ELSE
         kmap%nmesh = kpoints%nkp_grid
      END IF
      IF (ANY(kmap%nmesh < 1)) THEN
         IF (nred /= 1 .OR. &
             .NOT. tb_same_periodic_kpoint(kpoints%xkp(1:3, 1), [0.0_dp, 0.0_dp, 0.0_dp]) .OR. &
             ABS(kpoints%wkp(1) - 1.0_dp) > 1.0E-10_dp) &
            CPABORT("Descriptorless g-XTB k-point input requires a unit-weight Gamma point")
         kmap%nmesh = 1
      END IF
      kmap%nfull = PRODUCT(kmap%nmesh)
      IF (kmap%nfull < 1) CPABORT("Invalid complete g-xTB k-point mesh")
      ALLOCATE (kmap%xkp(3, kmap%nfull), kmap%weight(kmap%nfull), &
                kmap%source(kmap%nfull), kmap%symop(kmap%nfull), &
                kmap%fold_factor(kmap%nfull), kmap%representative(nred), &
                kmap%reduced_xkp(3, nred), kmap%reduced_weight(nred), reduced_weight(nred))
      kmap%reduced_xkp(:, :) = kpoints%xkp(1:3, 1:nred)
      kmap%reduced_weight(:) = kpoints%wkp(1:nred)

      IF (explicit_general_mesh) THEN
         IF (SIZE(kpoints%xkp_input, 2) /= kmap%nfull) &
            CPABORT("Explicit GENERAL mesh size changed during g-XTB reconstruction")
         kmap%xkp(1:3, 1:kmap%nfull) = kpoints%xkp_input(1:3, 1:kmap%nfull)
      ELSEIF (kmap%nfull == 1) THEN
         kmap%xkp(1:3, 1) = kpoints%xkp(1:3, 1)
      ELSE
         CALL full_grid_gen(kmap%nmesh, kmap%xkp, kmap%weight, kpoints%kp_shift, &
                            kpoints%gamma_centered)
      END IF
      kmap%weight = 1.0_dp/REAL(kmap%nfull, KIND=dp)
      kmap%source = 0
      kmap%symop = 0
      kmap%representative = 0

      ! Representatives themselves.
      DO ikred = 1, nred
         imatch = tb_find_periodic_kpoint(kmap%xkp, kpoints%xkp(1:3, ikred))
         kmap%representative(ikred) = imatch
         IF (imatch > 0 .AND. kmap%source(imatch) == 0) THEN
            kmap%source(imatch) = ikred
            kmap%symop(imatch) = 0
         END IF
      END DO
      IF (ANY(kmap%representative == 0)) &
         CPABORT("An irreducible g-xTB k-point is absent from its complete mesh")

      ! CP2K applies this reduction even when atom symmetry is disabled.
      DO ikred = 1, nred
         imatch = tb_find_periodic_kpoint(kmap%xkp, -kpoints%xkp(1:3, ikred))
         IF (imatch > 0 .AND. kmap%source(imatch) == 0) THEN
            kmap%source(imatch) = ikred
            kmap%symop(imatch) = -1
         END IF
      END DO

      ! Atomic stars (including antiunitary operations) complete the remaining mesh.
      IF (ASSOCIATED(kpoints%kp_sym)) THEN
         DO ikred = 1, nred
            kpsym => kpoints%kp_sym(ikred)%kpoint_sym
            IF (.NOT. ASSOCIATED(kpsym)) CYCLE
            IF (kpsym%nwred < 1) CYCLE
            DO isym = 1, kpsym%nwred
               imatch = tb_find_periodic_kpoint(kmap%xkp, kpsym%xkp(1:3, isym))
               IF (imatch > 0 .AND. kmap%source(imatch) == 0) THEN
                  kmap%source(imatch) = ikred
                  kmap%symop(imatch) = isym
               END IF
            END DO
         END DO
      END IF

      reduced_weight = 0.0_dp
      DO ifull = 1, kmap%nfull
         IF (kmap%source(ifull) == 0) &
            CPABORT("CP2K symmetry orbits do not span the complete g-xTB k-point mesh")
         ikred = kmap%source(ifull)
         reduced_weight(ikred) = reduced_weight(ikred) + kmap%weight(ifull)
         kmap%fold_factor(ifull) = kmap%weight(ifull)/kpoints%wkp(ikred)
      END DO
      IF (MAXVAL(ABS(reduced_weight - kpoints%wkp(1:nred))) > 1.0E-8_dp) &
         CPABORT("CP2K irreducible k-point weights disagree with their reconstructed stars")

      ! Retain every group operation that maps a representative onto a full point.  The first
      ! nwght entries contain one operation per distinct star member; nwred also contains the
      ! alternative little-group operations.  Averaging all of them projects an SCF density that
      ! is not itself little-group invariant onto the symmetry-reduced complete mesh.
      ALLOCATE (kmap%op_count(kmap%nfull), source=0)
      max_ops = 1
      DO ifull = 1, kmap%nfull
         ikred = kmap%source(ifull)
         IF (ASSOCIATED(kpoints%kp_sym)) THEN
            kpsym => kpoints%kp_sym(ikred)%kpoint_sym
         ELSE
            NULLIFY (kpsym)
         END IF
         IF (ASSOCIATED(kpsym)) THEN
            IF (kpsym%nwred > 0) THEN
               IF (kpsym%nwght < 1 .OR. MOD(kpsym%nwred, kpsym%nwght) /= 0) &
                  CPABORT("Inconsistent CP2K g-XTB k-point orbit dimensions")
               expected_ops = kpsym%nwred/kpsym%nwght
               DO isym = 1, kpsym%nwred
                  IF (tb_same_periodic_kpoint(kpsym%xkp(1:3, isym), &
                                              kmap%xkp(1:3, ifull))) &
                     kmap%op_count(ifull) = kmap%op_count(ifull) + 1
               END DO
               IF (kmap%op_count(ifull) /= expected_ops) &
                  CPABORT("Incomplete CP2K operation coset for a g-XTB k-point star member")
            END IF
         END IF
         IF (kmap%op_count(ifull) == 0) kmap%op_count(ifull) = 1
         max_ops = MAX(max_ops, kmap%op_count(ifull))
      END DO
      ALLOCATE (kmap%symops(max_ops, kmap%nfull), source=0)
      ALLOCATE (kmap%rotation_slot(max_ops, kmap%nfull), source=0)
      DO ifull = 1, kmap%nfull
         ikred = kmap%source(ifull)
         iop = 0
         IF (ASSOCIATED(kpoints%kp_sym)) THEN
            kpsym => kpoints%kp_sym(ikred)%kpoint_sym
         ELSE
            NULLIFY (kpsym)
         END IF
         IF (ASSOCIATED(kpsym)) THEN
            IF (kpsym%nwred > 0) THEN
               DO isym = 1, kpsym%nwred
                  IF (tb_same_periodic_kpoint(kpsym%xkp(1:3, isym), &
                                              kmap%xkp(1:3, ifull))) THEN
                     iop = iop + 1
                     kmap%symops(iop, ifull) = isym
                  END IF
               END DO
            END IF
         END IF
         IF (iop == 0) THEN
            iop = 1
            IF (tb_same_periodic_kpoint(kmap%xkp(1:3, ifull), &
                                        kpoints%xkp(1:3, ikred))) THEN
               kmap%symops(iop, ifull) = 0
            ELSEIF (tb_same_periodic_kpoint(kmap%xkp(1:3, ifull), &
                                            -kpoints%xkp(1:3, ikred))) THEN
               kmap%symops(iop, ifull) = -1
            ELSE
               CPABORT("No CP2K operation maps an irreducible g-xTB k-point onto its full mesh")
            END IF
         END IF
         CPASSERT(iop == kmap%op_count(ifull))
         kmap%symop(ifull) = kmap%symops(1, ifull)
         DO iop = 1, kmap%op_count(ifull)
            isym = kmap%symops(iop, ifull)
            IF (isym <= 0) CYCLE
            CPASSERT(ASSOCIATED(kpsym))
            CPASSERT(isym <= kpsym%nwred)
            CPASSERT(ASSOCIATED(kpoints%ibrot))
            ir = ABS(kpsym%rotp(isym))
            ira = 0
            DO jr = 1, SIZE(kpoints%ibrot)
               IF (ir == kpoints%ibrot(jr)) ira = jr
            END DO
            IF (ira == 0) &
               CPABORT("Missing CP2K AO rotation for a cached g-xTB symmetry operation")
            kmap%rotation_slot(iop, ifull) = ira
         END DO
      END DO

      kmap%identity = kmap%nfull == nred
      IF (kmap%identity) THEN
         DO ifull = 1, kmap%nfull
            IF (kmap%source(ifull) /= ifull .OR. kmap%op_count(ifull) /= 1 .OR. &
                kmap%symops(1, ifull) /= 0) THEN
               kmap%identity = .FALSE.
               EXIT
            END IF
         END DO
      END IF
      DEALLOCATE (reduced_weight)

   END SUBROUTINE tb_build_gxtb_kmesh_map

! **************************************************************************************************
!> \brief Check whether a cached mesh descriptor still belongs to this k-point object.
!> \param kpoints ...
!> \param kmap ...
!> \return true when all cached topology inputs and referenced symmetry operations still match
! **************************************************************************************************
   LOGICAL FUNCTION tb_gxtb_kmesh_map_matches(kpoints, kmap) RESULT(matches)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), INTENT(IN)           :: kmap

      INTEGER                                            :: current_ops, ifull, ikred, iop, isym
      LOGICAL                                            :: explicit_general
      TYPE(kpoint_sym_type), POINTER                     :: kpsym

      NULLIFY (kpsym)
      matches = .FALSE.
      explicit_general = TRIM(kpoints%kp_scheme) == "GENERAL" .AND. &
                         ASSOCIATED(kpoints%xkp_input) .AND. ASSOCIATED(kpoints%wkp_input)
      IF (kmap%nred /= kpoints%nkp .OR. kmap%nfull /= PRODUCT(kmap%nmesh)) RETURN
      IF (kmap%explicit_general .NEQV. explicit_general) RETURN
      IF (TRIM(kmap%scheme) /= TRIM(kpoints%kp_scheme)) RETURN
      IF (kmap%full_grid .NEQV. kpoints%full_grid) RETURN
      IF (kmap%gamma_centered .NEQV. kpoints%gamma_centered) RETURN
      IF (kmap%inversion_only .NEQV. kpoints%inversion_symmetry_only) RETURN
      IF (kmap%symmetry_backend /= kpoints%symmetry_backend .OR. &
          kmap%symmetry_method /= kpoints%symmetry_reduction_method) RETURN
      IF (ANY(kmap%kp_shift /= kpoints%kp_shift)) RETURN
      IF (.NOT. ALLOCATED(kmap%reduced_xkp) .OR. .NOT. ALLOCATED(kmap%reduced_weight) .OR. &
          .NOT. ALLOCATED(kmap%representative) .OR. .NOT. ALLOCATED(kmap%source) .OR. &
          .NOT. ALLOCATED(kmap%op_count) .OR. .NOT. ALLOCATED(kmap%symops) .OR. &
          .NOT. ALLOCATED(kmap%rotation_slot) .OR. &
          .NOT. ALLOCATED(kmap%xkp) .OR. .NOT. ALLOCATED(kmap%weight)) RETURN
      IF (ANY(SHAPE(kmap%reduced_xkp) /= [3, kmap%nred]) .OR. &
          SIZE(kmap%reduced_weight) /= kmap%nred .OR. &
          SIZE(kmap%representative) /= kmap%nred .OR. &
          ANY(SHAPE(kmap%xkp) /= [3, kmap%nfull]) .OR. &
          SIZE(kmap%weight) /= kmap%nfull .OR. SIZE(kmap%source) /= kmap%nfull .OR. &
          SIZE(kmap%op_count) /= kmap%nfull .OR. SIZE(kmap%symops, 2) /= kmap%nfull .OR. &
          ANY(SHAPE(kmap%rotation_slot) /= SHAPE(kmap%symops))) RETURN
      IF (.NOT. ASSOCIATED(kpoints%xkp) .OR. .NOT. ASSOCIATED(kpoints%wkp)) RETURN
      IF (ANY(kmap%reduced_xkp /= kpoints%xkp(1:3, 1:kmap%nred)) .OR. &
          ANY(kmap%reduced_weight /= kpoints%wkp(1:kmap%nred))) RETURN

      IF (kmap%explicit_general) THEN
         IF (.NOT. ASSOCIATED(kpoints%xkp_input) .OR. .NOT. ASSOCIATED(kpoints%wkp_input)) RETURN
         IF (SIZE(kpoints%xkp_input, 1) < 3 .OR. &
             SIZE(kpoints%xkp_input, 2) /= kmap%nfull .OR. &
             SIZE(kpoints%wkp_input) /= kmap%nfull) RETURN
         IF (ANY(kmap%xkp /= kpoints%xkp_input(1:3, 1:kmap%nfull)) .OR. &
             MAXVAL(ABS(kmap%weight - kpoints%wkp_input(1:kmap%nfull))) > 1.0E-10_dp) RETURN
      ELSE
         IF (ANY(kpoints%nkp_grid > 0 .AND. kpoints%nkp_grid /= kmap%nmesh)) RETURN
      END IF

      DO ikred = 1, kmap%nred
         ifull = kmap%representative(ikred)
         IF (ifull < 1 .OR. ifull > kmap%nfull) RETURN
         IF (.NOT. tb_same_periodic_kpoint(kmap%xkp(1:3, ifull), &
                                           kpoints%xkp(1:3, ikred))) RETURN
      END DO
      DO ifull = 1, kmap%nfull
         ikred = kmap%source(ifull)
         IF (ikred < 1 .OR. ikred > kmap%nred) RETURN
         IF (kmap%op_count(ifull) < 1 .OR. &
             kmap%op_count(ifull) > SIZE(kmap%symops, 1)) RETURN
         IF (ASSOCIATED(kpoints%kp_sym)) THEN
            kpsym => kpoints%kp_sym(ikred)%kpoint_sym
         ELSE
            NULLIFY (kpsym)
         END IF
         current_ops = 0
         IF (ASSOCIATED(kpsym)) THEN
            IF (kpsym%nwred > 0 .AND. ASSOCIATED(kpsym%xkp)) THEN
               DO isym = 1, kpsym%nwred
                  IF (tb_same_periodic_kpoint(kpsym%xkp(1:3, isym), &
                                              kmap%xkp(1:3, ifull))) current_ops = current_ops + 1
               END DO
            END IF
         END IF
         IF (current_ops == 0) current_ops = 1
         IF (current_ops /= kmap%op_count(ifull)) RETURN
         DO iop = 1, kmap%op_count(ifull)
            isym = kmap%symops(iop, ifull)
            IF (isym == 0) THEN
               IF (kmap%rotation_slot(iop, ifull) /= 0) RETURN
               IF (.NOT. tb_same_periodic_kpoint(kpoints%xkp(1:3, ikred), &
                                                 kmap%xkp(1:3, ifull))) RETURN
               CYCLE
            ELSEIF (isym == -1) THEN
               IF (kmap%rotation_slot(iop, ifull) /= 0) RETURN
               IF (.NOT. tb_same_periodic_kpoint(-kpoints%xkp(1:3, ikred), &
                                                 kmap%xkp(1:3, ifull))) RETURN
               CYCLE
            ELSEIF (isym < -1) THEN
               RETURN
            END IF
            IF (.NOT. ASSOCIATED(kpsym)) RETURN
            IF (isym > kpsym%nwred) RETURN
            IF (.NOT. ASSOCIATED(kpsym%xkp) .OR. .NOT. ASSOCIATED(kpsym%rotp)) RETURN
            IF (.NOT. tb_same_periodic_kpoint(kpsym%xkp(1:3, isym), &
                                              kmap%xkp(1:3, ifull))) RETURN
            IF (.NOT. ASSOCIATED(kpoints%ibrot)) RETURN
            IF (kmap%rotation_slot(iop, ifull) < 1 .OR. &
                kmap%rotation_slot(iop, ifull) > SIZE(kpoints%ibrot)) RETURN
            IF (ABS(kpsym%rotp(isym)) /= &
                kpoints%ibrot(kmap%rotation_slot(iop, ifull))) RETURN
         END DO
      END DO
      matches = .TRUE.

   END FUNCTION tb_gxtb_kmesh_map_matches

! **************************************************************************************************
!> \brief Build or return the geometry-static g-xTB complete-mesh descriptor.
!> \param tb ...
!> \param kpoints ...
!> \param matrix_s ...
!> \param sab_nl ...
!> \param cell_to_index ...
!>
!> The GENERAL validator and star map are invariant over an SCC cycle. They are invalidated when
!> tb_update_geometry observes changed coordinates or cell, and a cheap fingerprint also guards
!> against a replaced/reinitialized k-point object. Phase calibration has its own inexpensive
!> readiness test and is therefore called on the concrete current symmetry arrays every time.
! **************************************************************************************************
   SUBROUTINE tb_prepare_gxtb_kmesh_map(tb, kpoints, matrix_s, sab_nl, cell_to_index)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index

      IF (ALLOCATED(tb%gxtb_kmesh_map)) THEN
         IF (.NOT. tb_gxtb_kmesh_map_matches(kpoints, tb%gxtb_kmesh_map)) &
            CALL tb_reset_gxtb_kpoint_cache(tb)
      END IF
      IF (.NOT. ALLOCATED(tb%gxtb_kmesh_map)) THEN
         ALLOCATE (tb%gxtb_kmesh_map)
         CALL tb_build_gxtb_kmesh_map(kpoints, tb%gxtb_kmesh_map)
         tb%gxtb_density_map_validated = .FALSE.
         tb%gxtb_transform_validated = .FALSE.
      END IF

      IF (ANY(tb%gxtb_kmesh_map%symop > 0)) &
         CALL calibrate_symmetry_phases(kpoints, matrix_s, matrix_s(1, 1)%matrix, &
                                        sab_nl, cell_to_index, include_inactive=.TRUE.)

   END SUBROUTINE tb_prepare_gxtb_kmesh_map

! **************************************************************************************************
!> \brief Return the dense tblite AO layout for all primitive-cell atoms.
!> \param tb ...
!> \param ao_first ...
!> \param ao_size ...
! **************************************************************************************************
   SUBROUTINE tb_gxtb_ao_layout(tb, ao_first, ao_size)

      TYPE(tblite_type), POINTER                         :: tb
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: ao_first, ao_size

      INTEGER                                            :: iatom, nat

      nat = tb%mol%nat
      ALLOCATE (ao_first(nat), ao_size(nat))
      DO iatom = 1, nat
         CALL tb_atom_ao_range(tb, iatom, ao_first(iatom), ao_size(iatom))
         IF (iatom > 1) THEN
            CPASSERT(ao_first(iatom) == ao_first(iatom - 1) + ao_size(iatom - 1))
         END IF
      END DO
      CPASSERT(ao_first(nat) + ao_size(nat) - 1 == tb%calc%bas%nao)

   END SUBROUTINE tb_gxtb_ao_layout

! **************************************************************************************************
!> \brief Apply a CP2K atom/AO k-point symmetry transform or its real adjoint.
!> \param kpoints ...
!> \param ikred ...
!> \param isym 0 direct, -1 pure time reversal, positive kp_sym entry
!> \param rotation_slot cached index into kpoints%kind_rotmat; zero for identity/time reversal
!> \param ao_first ...
!> \param ao_size ...
!> \param source ...
!> \param TARGET ...
!> \param adjoint apply the adjoint needed to fold a full-mesh response
!>
!> This is the dense Hermitian-matrix analogue of kpoint_methods::symtrans_phase.  In particular,
!> phase_mode and phase_gauge are selected by CP2K's overlap-covariance calibration; no tblite
!> phase convention is introduced here.  The kgphase used by CP2K's Wannier-MO transform belongs
!> to a different gauge and is deliberately not part of this AO-matrix transform.
! **************************************************************************************************
   SUBROUTINE tb_transform_gxtb_kpoint_matrix(kpoints, ikred, isym, rotation_slot, ao_first, ao_size, &
                                              source, TARGET, adjoint)

      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: ikred, isym, rotation_slot
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_first, ao_size
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: source
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: target
      LOGICAL, INTENT(IN)                                :: adjoint

      COMPLEX(KIND=dp)                                   :: phase
      INTEGER                                            :: i_first, i_last, iatom, ikind, ip, &
                                                            ip_first, ip_last, j_first, j_last, &
                                                            jatom, jkind, jp, jp_first, jp_last
      INTEGER, DIMENSION(3)                              :: shift
      INTEGER, DIMENSION(:, :), POINTER                  :: phase_fcell
      LOGICAL                                            :: reverse_phase, time_reversal
      REAL(KIND=dp)                                      :: arg
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: kroti, krotj
      TYPE(kind_rotmat_type), DIMENSION(:), POINTER      :: kind_rot
      TYPE(kpoint_sym_type), POINTER                     :: kpsym

      CPASSERT(SIZE(source, 1) == SIZE(source, 2))
      CPASSERT(ALL(SHAPE(source) == SHAPE(TARGET)))
      IF (isym == 0) THEN
         CPASSERT(rotation_slot == 0)
         TARGET = source
         RETURN
      ELSEIF (isym == -1) THEN
         CPASSERT(rotation_slot == 0)
         TARGET = CONJG(source)
         RETURN
      END IF

      CPASSERT(ikred >= 1 .AND. ikred <= kpoints%nkp)
      CPASSERT(ASSOCIATED(kpoints%kp_sym))
      kpsym => kpoints%kp_sym(ikred)%kpoint_sym
      CPASSERT(ASSOCIATED(kpsym))
      CPASSERT(kpsym%nwred > 0)
      CPASSERT(isym >= 1 .AND. isym <= kpsym%nwred)
      CPASSERT(ASSOCIATED(kpsym%phase_mode) .AND. ASSOCIATED(kpsym%phase_gauge))
      CPASSERT(kpsym%phase_mode(isym) > 0 .AND. kpsym%phase_gauge(isym) > 0)
      CPASSERT(ASSOCIATED(kpoints%atype) .AND. ASSOCIATED(kpoints%ibrot))
      CPASSERT(ASSOCIATED(kpoints%kind_rotmat))
      CPASSERT(SIZE(ao_first) == SIZE(kpoints%atype))
      CPASSERT(SIZE(ao_size) == SIZE(kpoints%atype))

      reverse_phase = kpsym%phase_mode(isym) == 2
      SELECT CASE (kpsym%phase_gauge(isym))
      CASE (kpoint_phase_gauge_pbc)
         phase_fcell => kpsym%fcell_gauge(:, :, isym)
      CASE (kpoint_phase_gauge_input)
         phase_fcell => kpsym%fcell(:, :, isym)
      CASE DEFAULT
         CPABORT("Unknown calibrated CP2K k-point phase gauge")
      END SELECT

      CPASSERT(rotation_slot >= 1 .AND. rotation_slot <= SIZE(kpoints%ibrot))
      CPASSERT(ABS(kpsym%rotp(isym)) == kpoints%ibrot(rotation_slot))
      kind_rot => kpoints%kind_rotmat(rotation_slot, :)
      time_reversal = kpsym%rotp(isym) < 0
      TARGET = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      DO iatom = 1, SIZE(ao_first)
         ip = kpsym%f0(iatom, isym)
         CPASSERT(ip >= 1 .AND. ip <= SIZE(ao_first))
         ikind = kpoints%atype(iatom)
         kroti => kind_rot(ikind)%rmat
         CPASSERT(ASSOCIATED(kroti))
         i_first = ao_first(iatom)
         i_last = i_first + ao_size(iatom) - 1
         ip_first = ao_first(ip)
         ip_last = ip_first + ao_size(ip) - 1
         CPASSERT(SIZE(kroti, 1) == ao_size(ip) .AND. SIZE(kroti, 2) == ao_size(iatom))
         DO jatom = 1, SIZE(ao_first)
            jp = kpsym%f0(jatom, isym)
            CPASSERT(jp >= 1 .AND. jp <= SIZE(ao_first))
            jkind = kpoints%atype(jatom)
            krotj => kind_rot(jkind)%rmat
            CPASSERT(ASSOCIATED(krotj))
            j_first = ao_first(jatom)
            j_last = j_first + ao_size(jatom) - 1
            jp_first = ao_first(jp)
            jp_last = jp_first + ao_size(jp) - 1
            CPASSERT(SIZE(krotj, 1) == ao_size(jp) .AND. SIZE(krotj, 2) == ao_size(jatom))

            IF (reverse_phase) THEN
               shift = phase_fcell(1:3, iatom) - phase_fcell(1:3, jatom)
               arg = DOT_PRODUCT(REAL(shift, KIND=dp), kpsym%xkp(1:3, isym))
            ELSE
               shift = phase_fcell(1:3, jatom) - phase_fcell(1:3, iatom)
               arg = DOT_PRODUCT(REAL(shift, KIND=dp), kpsym%xkp(1:3, isym))
            END IF
            phase = CMPLX(COS(twopi*arg), -SIN(twopi*arg), KIND=dp)

            IF (.NOT. adjoint) THEN
               IF (time_reversal) THEN
                  TARGET(ip_first:ip_last, jp_first:jp_last) = phase* &
                                                               MATMUL(kroti, MATMUL(CONJG(source(i_first:i_last, j_first:j_last)), &
                                                                                    TRANSPOSE(krotj)))
               ELSE
                  TARGET(ip_first:ip_last, jp_first:jp_last) = phase* &
                                                               MATMUL(kroti, MATMUL(source(i_first:i_last, j_first:j_last), &
                                                                                    TRANSPOSE(krotj)))
               END IF
            ELSE
               IF (time_reversal) THEN
                  TARGET(i_first:i_last, j_first:j_last) = phase*CONJG( &
                                                           MATMUL(TRANSPOSE(kroti), &
                                                                  MATMUL(source(ip_first:ip_last, jp_first:jp_last), krotj)))
               ELSE
                  TARGET(i_first:i_last, j_first:j_last) = CONJG(phase)* &
                                                           MATMUL(TRANSPOSE(kroti), &
                                                                  MATMUL(source(ip_first:ip_last, jp_first:jp_last), krotj))
               END IF
            END IF
         END DO
      END DO

   END SUBROUTINE tb_transform_gxtb_kpoint_matrix

! **************************************************************************************************
!> \brief Expand one irreducible matrix into one complete-mesh star member.
!> \param kpoints ...
!> \param kmap ...
!> \param ao_first ...
!> \param ao_size ...
!> \param ifull complete-mesh destination index
!> \param representative matrix at the irreducible source k point
!> \param full_block projected matrix at ifull
!> \param work caller-owned transform workspace
!>
!> The operation-coset average is the same projector used by the full-array density expansion.
!> Keeping the workspace outside this helper permits a future streamed caller to reuse it without
!> allocating for every star member.  The cached CP2K rotation and calibrated phase gauge are used
!> unchanged, including negative rotp entries which request the antiunitary transform.
! **************************************************************************************************
   SUBROUTINE tb_expand_gxtb_kpoint_block(kpoints, kmap, ao_first, ao_size, ifull, &
                                          representative, full_block, work)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), INTENT(IN)           :: kmap
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_first, ao_size
      INTEGER, INTENT(IN)                                :: ifull
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: representative
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: full_block, work

      INTEGER                                            :: ikred, iop

      CPASSERT(kmap%nred == kpoints%nkp)
      CPASSERT(ifull >= 1)
      CPASSERT(ifull <= kmap%nfull)
      CPASSERT(ALLOCATED(kmap%source))
      CPASSERT(ALLOCATED(kmap%op_count))
      CPASSERT(ALLOCATED(kmap%symops))
      CPASSERT(ALLOCATED(kmap%rotation_slot))
      CPASSERT(SIZE(kmap%source) == kmap%nfull)
      CPASSERT(SIZE(kmap%op_count) == kmap%nfull)
      CPASSERT(SIZE(kmap%symops, 2) == kmap%nfull)
      CPASSERT(SIZE(kmap%rotation_slot, 1) == SIZE(kmap%symops, 1))
      CPASSERT(SIZE(kmap%rotation_slot, 2) == SIZE(kmap%symops, 2))
      CPASSERT(SIZE(representative, 1) == SIZE(representative, 2))
      CPASSERT(SIZE(full_block, 1) == SIZE(representative, 1))
      CPASSERT(SIZE(full_block, 2) == SIZE(representative, 2))
      CPASSERT(SIZE(work, 1) == SIZE(representative, 1))
      CPASSERT(SIZE(work, 2) == SIZE(representative, 2))
      ikred = kmap%source(ifull)
      CPASSERT(ikred >= 1)
      CPASSERT(ikred <= kmap%nred)
      CPASSERT(kmap%op_count(ifull) >= 1)
      CPASSERT(kmap%op_count(ifull) <= SIZE(kmap%symops, 1))

      full_block = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO iop = 1, kmap%op_count(ifull)
         CALL tb_transform_gxtb_kpoint_matrix(kpoints, ikred, kmap%symops(iop, ifull), &
                                              kmap%rotation_slot(iop, ifull), ao_first, ao_size, &
                                              representative, work, .FALSE.)
         full_block(:, :) = full_block + work/REAL(kmap%op_count(ifull), KIND=dp)
      END DO
      full_block(:, :) = 0.5_dp*(full_block + CONJG(TRANSPOSE(full_block)))

   END SUBROUTINE tb_expand_gxtb_kpoint_block

! **************************************************************************************************
!> \brief Expand irreducible spin densities to the complete regular k mesh.
!> \param kpoints ...
!> \param kmap ...
!> \param ao_first ...
!> \param ao_size ...
!> \param density_full full-sized buffer with irreducible inputs in kmap%representative slots
!> \param validate ...
! **************************************************************************************************
   SUBROUTINE tb_expand_gxtb_kpoint_density(kpoints, kmap, ao_first, ao_size, &
                                            density_full, validate)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), INTENT(IN)           :: kmap
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_first, ao_size
      COMPLEX(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT)                                   :: density_full
      LOGICAL, INTENT(IN)                                :: validate

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: transformed
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: representative_density
      INTEGER                                            :: ifull, ikred, iop, irepresentative, ispin
      REAL(KIND=dp)                                      :: residual, scale, symmetry_tolerance

      CPASSERT(SIZE(density_full, 4) == kmap%nfull)
      ALLOCATE (transformed(SIZE(density_full, 1), SIZE(density_full, 2)), &
                representative_density(SIZE(density_full, 1), SIZE(density_full, 2), &
                                       SIZE(density_full, 3)))
      symmetry_tolerance = MAX(1.0E-10_dp, 100.0_dp*kpoints%eps_geo)
      DO ikred = 1, kpoints%nkp
         irepresentative = kmap%representative(ikred)
         representative_density(:, :, :) = density_full(:, :, :, irepresentative)
         DO ifull = 1, kmap%nfull
            IF (kmap%source(ifull) /= ikred) CYCLE
            DO ispin = 1, SIZE(density_full, 3)
               density_full(:, :, ispin, ifull) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
               DO iop = 1, kmap%op_count(ifull)
                  CALL tb_transform_gxtb_kpoint_matrix(kpoints, ikred, &
                                                       kmap%symops(iop, ifull), &
                                                       kmap%rotation_slot(iop, ifull), ao_first, ao_size, &
                                                       representative_density(:, :, ispin), &
                                                       transformed, .FALSE.)
                  density_full(:, :, ispin, ifull) = density_full(:, :, ispin, ifull) + &
                                                     transformed/REAL(kmap%op_count(ifull), KIND=dp)
               END DO
               density_full(:, :, ispin, ifull) = &
                  0.5_dp*(density_full(:, :, ispin, ifull) + &
                          CONJG(TRANSPOSE(density_full(:, :, ispin, ifull))))
            END DO
         END DO
      END DO

      IF (validate) THEN
         ! The projector and static operation map do not change between SCC iterations. Validate
         ! one physical density after each geometry change, rather than repeating this O(Nk*Na^2)
         ! diagnostic for every subsequent SCC step.
         DO ifull = 1, kmap%nfull
            irepresentative = kmap%representative(kmap%source(ifull))
            CPASSERT(irepresentative > 0)
            DO iop = 1, kmap%op_count(ifull)
               DO ispin = 1, SIZE(density_full, 3)
                  CALL tb_transform_gxtb_kpoint_matrix(kpoints, kmap%source(ifull), &
                                                       kmap%symops(iop, ifull), &
                                                       kmap%rotation_slot(iop, ifull), ao_first, ao_size, &
                                                       density_full(:, :, ispin, irepresentative), &
                                                       transformed, .FALSE.)
                  scale = MAX(1.0_dp, MAXVAL(ABS(density_full(:, :, ispin, ifull))))
                  residual = MAXVAL(ABS(transformed - density_full(:, :, ispin, ifull)))/scale
                  IF (residual > symmetry_tolerance) &
                     CPABORT("Projected g-xTB density is not invariant under its CP2K little group")
               END DO
            END DO
         END DO
      END IF
      DEALLOCATE (representative_density, transformed)

   END SUBROUTINE tb_expand_gxtb_kpoint_density

! **************************************************************************************************
!> \brief Accumulate the weighted real adjoint of one complete-mesh response block.
!> \param kpoints ...
!> \param kmap ...
!> \param ao_first ...
!> \param ao_size ...
!> \param ifull complete-mesh source index
!> \param full_block unweighted response at ifull
!> \param reduced_block accumulator for kmap%source(ifull)
!> \param work caller-owned transform workspace
!>
!> fold_factor is w_full/w_irred.  Averaging over the complete operation coset and applying the
!> real adjoint therefore preserves the weighted variational contraction for both unitary and
!> antiunitary star members.  The caller hermitizes reduced_block after all members belonging to
!> the irreducible destination have been accumulated, as done by the full-array oracle below.
! **************************************************************************************************
   SUBROUTINE tb_accumulate_gxtb_kpoint_response(kpoints, kmap, ao_first, ao_size, ifull, &
                                                 full_block, reduced_block, work)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), INTENT(IN)           :: kmap
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_first, ao_size
      INTEGER, INTENT(IN)                                :: ifull
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: full_block
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT)   :: reduced_block
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: work

      INTEGER                                            :: ikred, iop

      CPASSERT(kmap%nred == kpoints%nkp)
      CPASSERT(ifull >= 1)
      CPASSERT(ifull <= kmap%nfull)
      CPASSERT(ALLOCATED(kmap%source))
      CPASSERT(ALLOCATED(kmap%op_count))
      CPASSERT(ALLOCATED(kmap%symops))
      CPASSERT(ALLOCATED(kmap%rotation_slot))
      CPASSERT(ALLOCATED(kmap%fold_factor))
      CPASSERT(SIZE(kmap%source) == kmap%nfull)
      CPASSERT(SIZE(kmap%op_count) == kmap%nfull)
      CPASSERT(SIZE(kmap%fold_factor) == kmap%nfull)
      CPASSERT(SIZE(kmap%symops, 2) == kmap%nfull)
      CPASSERT(SIZE(kmap%rotation_slot, 1) == SIZE(kmap%symops, 1))
      CPASSERT(SIZE(kmap%rotation_slot, 2) == SIZE(kmap%symops, 2))
      CPASSERT(SIZE(full_block, 1) == SIZE(full_block, 2))
      CPASSERT(SIZE(reduced_block, 1) == SIZE(full_block, 1))
      CPASSERT(SIZE(reduced_block, 2) == SIZE(full_block, 2))
      CPASSERT(SIZE(work, 1) == SIZE(full_block, 1))
      CPASSERT(SIZE(work, 2) == SIZE(full_block, 2))
      CPASSERT(ieee_is_finite(kmap%fold_factor(ifull)))
      CPASSERT(kmap%fold_factor(ifull) >= 0.0_dp)
      ikred = kmap%source(ifull)
      CPASSERT(ikred >= 1)
      CPASSERT(ikred <= kmap%nred)
      CPASSERT(kmap%op_count(ifull) >= 1)
      CPASSERT(kmap%op_count(ifull) <= SIZE(kmap%symops, 1))

      DO iop = 1, kmap%op_count(ifull)
         CALL tb_transform_gxtb_kpoint_matrix(kpoints, ikred, kmap%symops(iop, ifull), &
                                              kmap%rotation_slot(iop, ifull), ao_first, ao_size, &
                                              full_block, work, .TRUE.)
         reduced_block(:, :) = reduced_block + &
                               kmap%fold_factor(ifull)*work/ &
                               REAL(kmap%op_count(ifull), KIND=dp)
      END DO

   END SUBROUTINE tb_accumulate_gxtb_kpoint_response

! **************************************************************************************************
!> \brief Fold an unweighted complete-mesh Hermitian response to irreducible k-points.
!> \param kpoints ...
!> \param kmap ...
!> \param ao_first ...
!> \param ao_size ...
!> \param full ...
!> \param reduced ...
! **************************************************************************************************
   SUBROUTINE tb_fold_gxtb_kpoint_matrix(kpoints, kmap, ao_first, ao_size, full, reduced)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), INTENT(IN)           :: kmap
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_first, ao_size
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: full
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: reduced

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: work
      INTEGER                                            :: ifull, ikred, iop

      CPASSERT(SIZE(full, 3) == kmap%nfull)
      CPASSERT(SIZE(reduced, 3) == kpoints%nkp)
      ALLOCATE (work(SIZE(full, 1), SIZE(full, 2)))
      reduced = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO ifull = 1, kmap%nfull
         ikred = kmap%source(ifull)
         DO iop = 1, kmap%op_count(ifull)
            CALL tb_transform_gxtb_kpoint_matrix(kpoints, ikred, kmap%symops(iop, ifull), &
                                                 kmap%rotation_slot(iop, ifull), ao_first, ao_size, &
                                                 full(:, :, ifull), work, .TRUE.)
            reduced(:, :, ikred) = reduced(:, :, ikred) + &
                                   kmap%fold_factor(ifull)*work/ &
                                   REAL(kmap%op_count(ifull), KIND=dp)
         END DO
      END DO
      DO ikred = 1, kpoints%nkp
         reduced(:, :, ikred) = 0.5_dp*(reduced(:, :, ikred) + &
                                        CONJG(TRANSPOSE(reduced(:, :, ikred))))
      END DO
      DEALLOCATE (work)

   END SUBROUTINE tb_fold_gxtb_kpoint_matrix

! **************************************************************************************************
!> \brief Independently test the weighted variational duality of full and folded responses.
!>
!> This scalar check deliberately does not call the folding helper.  It contracts the current
!> irreducible density and its separately expanded complete-mesh image with the two response
!> representations, so that star weights, adjoints, and antiunitary operations are all tested.
!> \param kpoints ...
!> \param kmap ...
!> \param density_full ...
!> \param density_reduced ...
!> \param response_full ...
!> \param response_reduced ...
!> \param full_contraction ...
!> \param reduced_contraction ...
!> \param relative_residual ...
! **************************************************************************************************
   SUBROUTINE tb_qualification_gxtb_duality(kpoints, kmap, density_full, density_reduced, &
                                            response_full, response_reduced, full_contraction, &
                                            reduced_contraction, relative_residual)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), INTENT(IN)           :: kmap
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: density_full, density_reduced, &
                                                            response_full, response_reduced
      REAL(KIND=dp), INTENT(OUT)                         :: full_contraction, reduced_contraction, &
                                                            relative_residual

      INTEGER                                            :: ifull, ikred
      REAL(KIND=dp)                                      :: scale

      CPASSERT(SIZE(density_full, 3) == kmap%nfull)
      CPASSERT(SIZE(response_full, 3) == kmap%nfull)
      CPASSERT(SIZE(density_reduced, 3) == kmap%nred)
      CPASSERT(SIZE(response_reduced, 3) == kmap%nred)
      CPASSERT(ALL(SHAPE(density_full) == SHAPE(response_full)))
      CPASSERT(ALL(SHAPE(density_reduced) == SHAPE(response_reduced)))

      full_contraction = 0.0_dp
      DO ifull = 1, kmap%nfull
         full_contraction = full_contraction + kmap%weight(ifull)* &
                            REAL(SUM(CONJG(density_full(:, :, ifull))* &
                                     response_full(:, :, ifull)), KIND=dp)
      END DO
      reduced_contraction = 0.0_dp
      DO ikred = 1, kmap%nred
         reduced_contraction = reduced_contraction + kpoints%wkp(ikred)* &
                               REAL(SUM(CONJG(density_reduced(:, :, ikred))* &
                                        response_reduced(:, :, ikred)), KIND=dp)
      END DO
      scale = MAX(1.0E-30_dp, ABS(full_contraction), ABS(reduced_contraction))
      relative_residual = ABS(full_contraction - reduced_contraction)/scale

   END SUBROUTINE tb_qualification_gxtb_duality

! **************************************************************************************************
!> \brief Qualify the one-block expansion and folding helpers against the full-array oracle.
!> \param kpoints ...
!> \param kmap ...
!> \param ao_first ...
!> \param ao_size ...
!> \param overlap complete-mesh Hermitian test matrices
!> \param require_negative_rotp require a signed antiunitary space-group operation
!>
!> This check is called only from the once-per-map symmetry qualification.  It deliberately leaves
!> the production full-array expansion and folding routines unchanged.  Besides comparing the
!> accumulated blockwise fold to tb_fold_gxtb_kpoint_matrix, it checks the weighted real-adjoint
!> identity for each star member separately.  A dedicated qualification case can additionally
!> require a positive symmetry index whose signed rotp entry is negative.
! **************************************************************************************************
   SUBROUTINE tb_validate_gxtb_kpoint_block_helpers(kpoints, kmap, ao_first, ao_size, overlap, &
                                                    require_negative_rotp)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), INTENT(IN)           :: kmap
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_first, ao_size
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: overlap
      LOGICAL, INTENT(IN)                                :: require_negative_rotp

      CHARACTER(LEN=256)                                 :: error_message
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: expanded, folded, probe_representative, &
                                                            work
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: blockwise, oracle, probe_full
      INTEGER                                            :: i, ifull, iop, irepresentative, isym, j, &
                                                            negative_rotp_seen, source
      REAL(KIND=dp)                                      :: lhs, residual, rhs, scale, &
                                                            symmetry_tolerance
      TYPE(kpoint_sym_type), POINTER                     :: kpsym
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (kpsym, para_env)
      para_env => kpoints%blacs_env_all%para_env
      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(SIZE(overlap, 3) == kmap%nfull)
      ALLOCATE (expanded(SIZE(overlap, 1), SIZE(overlap, 2)), &
                folded(SIZE(overlap, 1), SIZE(overlap, 2)), &
                probe_representative(SIZE(overlap, 1), SIZE(overlap, 2)), &
                work(SIZE(overlap, 1), SIZE(overlap, 2)), &
                blockwise(SIZE(overlap, 1), SIZE(overlap, 2), kmap%nred), &
                oracle(SIZE(overlap, 1), SIZE(overlap, 2), kmap%nred), &
                probe_full(SIZE(overlap, 1), SIZE(overlap, 2), kmap%nfull))
      blockwise = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      symmetry_tolerance = MAX(1.0E-6_dp, 100.0_dp*kpoints%eps_geo)
      negative_rotp_seen = 0
      DO ifull = 1, kmap%nfull
         source = kmap%source(ifull)
         IF (.NOT. ASSOCIATED(kpoints%kp_sym)) CYCLE
         kpsym => kpoints%kp_sym(source)%kpoint_sym
         IF (.NOT. ASSOCIATED(kpsym) .OR. .NOT. ASSOCIATED(kpsym%rotp)) CYCLE
         DO iop = 1, kmap%op_count(ifull)
            isym = kmap%symops(iop, ifull)
            IF (isym > 0) THEN
               CPASSERT(isym <= SIZE(kpsym%rotp))
               IF (kpsym%rotp(isym) < 0) negative_rotp_seen = negative_rotp_seen + 1
            END IF
         END DO
      END DO
      IF (require_negative_rotp .AND. negative_rotp_seen == 0) &
         CPABORT("Block-helper qualification did not exercise a negative signed rotp entry")
      IF (para_env%is_source()) &
         WRITE (*, '(A,I0)') " GXTB-QUALIFICATION_ONLY BLOCK-HELPERS negative_rotp=", &
         negative_rotp_seen
      DO ifull = 1, kmap%nfull
         DO j = 1, SIZE(overlap, 2)
            DO i = 1, j
               probe_full(i, j, ifull) = &
                  CMPLX(0.01_dp*REAL(i + 2*j + 3*ifull, KIND=dp), &
                        0.003_dp*REAL((i - j)*(ifull + 1), KIND=dp), KIND=dp)
               probe_full(j, i, ifull) = CONJG(probe_full(i, j, ifull))
            END DO
         END DO
      END DO

      DO ifull = 1, kmap%nfull
         source = kmap%source(ifull)
         CPASSERT(source >= 1)
         CPASSERT(source <= kmap%nred)
         irepresentative = kmap%representative(source)
         CPASSERT(irepresentative >= 1)
         CPASSERT(irepresentative <= kmap%nfull)
         CALL tb_expand_gxtb_kpoint_block(kpoints, kmap, ao_first, ao_size, ifull, &
                                          overlap(:, :, irepresentative), expanded, work)
         scale = MAX(1.0_dp, MAXVAL(ABS(overlap(:, :, ifull))))
         residual = MAXVAL(ABS(expanded - overlap(:, :, ifull)))/scale
         IF (residual > symmetry_tolerance) THEN
            WRITE (error_message, '(A,I0,A,I0,A,ES12.4)') &
               "Blockwise CP2K symmetry expansion failed for full k-point ", ifull, &
               ", irreducible source ", source, "; residual=", residual
            CALL cp_abort(__LOCATION__, TRIM(error_message))
         END IF

         DO j = 1, SIZE(overlap, 2)
            DO i = 1, j
               probe_representative(i, j) = &
                  CMPLX(0.02_dp*REAL(2*i + j + source, KIND=dp), &
                        0.007_dp*REAL((i - j)*(source + 2), KIND=dp), KIND=dp)
               probe_representative(j, i) = CONJG(probe_representative(i, j))
            END DO
         END DO
         CALL tb_expand_gxtb_kpoint_block(kpoints, kmap, ao_first, ao_size, ifull, &
                                          probe_representative, expanded, work)
         folded = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         CALL tb_accumulate_gxtb_kpoint_response(kpoints, kmap, ao_first, ao_size, ifull, &
                                                 probe_full(:, :, ifull), folded, work)
         blockwise(:, :, source) = blockwise(:, :, source) + folded

         lhs = kmap%weight(ifull)* &
               REAL(SUM(CONJG(probe_full(:, :, ifull))*expanded), KIND=dp)
         rhs = kpoints%wkp(source)* &
               REAL(SUM(CONJG(probe_representative)*folded), KIND=dp)
         residual = ABS(lhs - rhs)/MAX(1.0_dp, ABS(lhs), ABS(rhs))
         IF (residual > 1.0E-10_dp) THEN
            WRITE (error_message, '(A,I0,A,I0,A,ES12.4)') &
               "Weighted blockwise CP2K symmetry adjoint failed for full k-point ", ifull, &
               ", irreducible source ", source, "; residual=", residual
            CALL cp_abort(__LOCATION__, TRIM(error_message))
         END IF
      END DO

      DO source = 1, kmap%nred
         blockwise(:, :, source) = 0.5_dp*(blockwise(:, :, source) + &
                                           CONJG(TRANSPOSE(blockwise(:, :, source))))
      END DO
      CALL tb_fold_gxtb_kpoint_matrix(kpoints, kmap, ao_first, ao_size, probe_full, oracle)
      scale = MAX(1.0_dp, MAXVAL(ABS(oracle)))
      residual = MAXVAL(ABS(blockwise - oracle))/scale
      IF (residual > 1.0E-12_dp) THEN
         WRITE (error_message, '(A,ES12.4)') &
            "Blockwise CP2K symmetry fold disagrees with the full-array oracle; residual=", residual
         CALL cp_abort(__LOCATION__, TRIM(error_message))
      END IF
      DEALLOCATE (blockwise, expanded, folded, oracle, probe_full, probe_representative, work)

   END SUBROUTINE tb_validate_gxtb_kpoint_block_helpers

! **************************************************************************************************
!> \brief Validate the dense symmetry transform against direct CP2K overlap matrices.
!>
!> Besides checking S(k_full)=T[S(k_irred)], this verifies the real adjoint used for antiunitary
!> fold-back through a Frobenius dot-product identity before save_tblite sees a reduced mesh.
!> \param kpoints ...
!> \param kmap ...
!> \param ao_first ...
!> \param ao_size ...
!> \param overlap ...
! **************************************************************************************************
   SUBROUTINE tb_validate_gxtb_kmesh_transform(kpoints, kmap, ao_first, ao_size, overlap)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), INTENT(IN)           :: kmap
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_first, ao_size
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: overlap

      CHARACTER(LEN=256)                                 :: error_message
      CHARACTER(LEN=32)                                  :: qualification_value
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: adjoint, transformed
      INTEGER                                            :: ifull, iop, irepresentative, &
                                                            qualification_status, source
      LOGICAL                                            :: require_negative_rotp, &
                                                            run_block_qualification
      REAL(KIND=dp)                                      :: adjoint_residual, covariance_residual, &
                                                            lhs, rhs, scale

      CPASSERT(SIZE(overlap, 3) == kmap%nfull)
      qualification_status = 1
      qualification_value = ""
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_GXTB_QUALIFICATION_BLOCK_HELPERS", &
                                    qualification_value, STATUS=qualification_status)
      run_block_qualification = .FALSE.
      require_negative_rotp = .FALSE.
      IF (qualification_status == 0) THEN
         SELECT CASE (TRIM(ADJUSTL(qualification_value)))
         CASE ("", "0", "OFF", "off", "FALSE", "false")
            run_block_qualification = .FALSE.
         CASE ("1", "ON", "on", "TRUE", "true")
            run_block_qualification = .TRUE.
         CASE ("REQUIRE_NEGATIVE_ROTP", "require_negative_rotp")
            run_block_qualification = .TRUE.
            require_negative_rotp = .TRUE.
         CASE DEFAULT
            CPABORT("Unknown CP2K_GXTB_QUALIFICATION_BLOCK_HELPERS value")
         END SELECT
      END IF
      ALLOCATE (adjoint(SIZE(overlap, 1), SIZE(overlap, 2)), &
                transformed(SIZE(overlap, 1), SIZE(overlap, 2)))
      DO ifull = 1, kmap%nfull
         source = kmap%source(ifull)
         irepresentative = kmap%representative(source)
         IF (irepresentative == 0) &
            CPABORT("Could not locate an irreducible representative in the complete g-xTB mesh")
         DO iop = 1, kmap%op_count(ifull)
            CALL tb_transform_gxtb_kpoint_matrix(kpoints, source, kmap%symops(iop, ifull), &
                                                 kmap%rotation_slot(iop, ifull), ao_first, ao_size, &
                                                 overlap(:, :, irepresentative), &
                                                 transformed, .FALSE.)
            scale = MAX(1.0_dp, MAXVAL(ABS(overlap(:, :, ifull))))
            covariance_residual = MAXVAL(ABS(transformed - overlap(:, :, ifull)))/scale
            IF (covariance_residual > MAX(1.0E-6_dp, 100.0_dp*kpoints%eps_geo)) THEN
               WRITE (error_message, '(A,I0,A,I0,A,I0,A,ES12.4)') &
                  "Dense CP2K overlap symmetry failed for full k-point ", ifull, &
                  ", irreducible source ", source, ", symmetry operation ", &
                  kmap%symops(iop, ifull), "; residual=", covariance_residual
               CALL cp_abort(__LOCATION__, TRIM(error_message))
            END IF

            CALL tb_transform_gxtb_kpoint_matrix(kpoints, source, kmap%symops(iop, ifull), &
                                                 kmap%rotation_slot(iop, ifull), ao_first, ao_size, &
                                                 overlap(:, :, ifull), &
                                                 adjoint, .TRUE.)
            lhs = REAL(SUM(CONJG(overlap(:, :, ifull))*transformed), KIND=dp)
            rhs = REAL(SUM(CONJG(adjoint)*overlap(:, :, irepresentative)), KIND=dp)
            adjoint_residual = ABS(lhs - rhs)/MAX(1.0_dp, ABS(lhs), ABS(rhs))
            IF (adjoint_residual > 1.0E-10_dp) THEN
               WRITE (error_message, '(A,I0,A,I0,A,I0,A,ES12.4)') &
                  "Dense CP2K symmetry adjoint failed for full k-point ", ifull, &
                  ", irreducible source ", source, ", symmetry operation ", &
                  kmap%symops(iop, ifull), "; residual=", adjoint_residual
               CALL cp_abort(__LOCATION__, TRIM(error_message))
            END IF
         END DO
      END DO
      IF (run_block_qualification) &
         CALL tb_validate_gxtb_kpoint_block_helpers(kpoints, kmap, ao_first, ao_size, overlap, &
                                                    require_negative_rotp)
      DEALLOCATE (adjoint, transformed)

   END SUBROUTINE tb_validate_gxtb_kmesh_transform

! **************************************************************************************************
!> \brief Check that a complete-mesh response is covariant before symmetry folding.
!> \param kpoints ...
!> \param kmap ...
!> \param ao_first ...
!> \param ao_size ...
!> \param full ...
!> \param label ...
!> \param max_residual ...
! **************************************************************************************************
   SUBROUTINE tb_validate_gxtb_full_response(kpoints, kmap, ao_first, ao_size, full, label, max_residual)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), INTENT(IN)           :: kmap
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_first, ao_size
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: full
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: max_residual

      CHARACTER(LEN=256)                                 :: error_message
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: transformed
      INTEGER                                            :: ifull, iop, irepresentative, source
      REAL(KIND=dp)                                      :: residual, scale, symmetry_tolerance

      CPASSERT(SIZE(full, 3) == kmap%nfull)
      ALLOCATE (transformed(SIZE(full, 1), SIZE(full, 2)))
      symmetry_tolerance = MAX(1.0E-10_dp, 100.0_dp*kpoints%eps_geo)
      IF (PRESENT(max_residual)) max_residual = 0.0_dp
      DO ifull = 1, kmap%nfull
         source = kmap%source(ifull)
         irepresentative = kmap%representative(source)
         IF (irepresentative == 0) &
            CPABORT("Could not locate an irreducible representative in the complete g-xTB mesh")
         DO iop = 1, kmap%op_count(ifull)
            CALL tb_transform_gxtb_kpoint_matrix(kpoints, source, kmap%symops(iop, ifull), &
                                                 kmap%rotation_slot(iop, ifull), ao_first, ao_size, &
                                                 full(:, :, irepresentative), &
                                                 transformed, .FALSE.)
            scale = MAX(1.0_dp, MAXVAL(ABS(full(:, :, ifull))))
            residual = MAXVAL(ABS(transformed - full(:, :, ifull)))/scale
            IF (PRESENT(max_residual)) max_residual = MAX(max_residual, residual)
            IF (residual > symmetry_tolerance) THEN
               WRITE (error_message, '(A,A,A,I0,A,I0,A,I0,A,I0,A,ES12.4)') &
                  "Complete-mesh g-xTB ", TRIM(label), " violates CP2K symmetry at full k-point ", &
                  ifull, ", representative ", irepresentative, ", irreducible source ", source, &
                  ", symmetry operation ", kmap%symops(iop, ifull), "; residual=", residual
               CALL cp_abort(__LOCATION__, TRIM(error_message))
            END IF
         END DO
      END DO
      DEALLOCATE (transformed)

   END SUBROUTINE tb_validate_gxtb_full_response

! **************************************************************************************************
!> \brief Check that a folded response reproduces every complete-mesh symmetry image.
!> \param kpoints ...
!> \param kmap ...
!> \param ao_first ...
!> \param ao_size ...
!> \param full ...
!> \param reduced ...
!> \param label ...
!> \param max_residual ...
! **************************************************************************************************
   SUBROUTINE tb_validate_gxtb_folded_response(kpoints, kmap, ao_first, ao_size, full, reduced, label, &
                                               max_residual)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), INTENT(IN)           :: kmap
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_first, ao_size
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: full, reduced
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: max_residual

      CHARACTER(LEN=256)                                 :: error_message
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: transformed
      INTEGER                                            :: ifull, iop, source
      REAL(KIND=dp)                                      :: residual, scale, symmetry_tolerance

      CPASSERT(SIZE(full, 3) == kmap%nfull)
      CPASSERT(SIZE(reduced, 3) == kpoints%nkp)
      ALLOCATE (transformed(SIZE(full, 1), SIZE(full, 2)))
      symmetry_tolerance = MAX(1.0E-10_dp, 100.0_dp*kpoints%eps_geo)
      IF (PRESENT(max_residual)) max_residual = 0.0_dp
      DO ifull = 1, kmap%nfull
         source = kmap%source(ifull)
         DO iop = 1, kmap%op_count(ifull)
            CALL tb_transform_gxtb_kpoint_matrix(kpoints, source, kmap%symops(iop, ifull), &
                                                 kmap%rotation_slot(iop, ifull), ao_first, ao_size, &
                                                 reduced(:, :, source), &
                                                 transformed, .FALSE.)
            scale = MAX(1.0_dp, MAXVAL(ABS(full(:, :, ifull))))
            residual = MAXVAL(ABS(transformed - full(:, :, ifull)))/scale
            IF (PRESENT(max_residual)) max_residual = MAX(max_residual, residual)
            IF (residual > symmetry_tolerance) THEN
               WRITE (error_message, '(A,A,A,I0,A,I0,A,I0,A,ES12.4)') &
                  "Folded g-xTB ", TRIM(label), " violates CP2K symmetry at full k-point ", &
                  ifull, ", irreducible source ", source, ", symmetry operation ", &
                  kmap%symops(iop, ifull), "; residual=", residual
               CALL cp_abort(__LOCATION__, TRIM(error_message))
            END IF
         END DO
      END DO
      DEALLOCATE (transformed)

   END SUBROUTINE tb_validate_gxtb_folded_response

! **************************************************************************************************
!> \brief Project irreducible Hermitian responses into their CP2K little-group invariant spaces.
!>
!> A symmetry-reduced SCF represents every star by one matrix.  Numerical diagonalization can
!> introduce a small component that is not invariant under operations which leave that
!> representative k point fixed.  Such a component is outside the physical reduced variational
!> space and must not enter the model-native Fock-DIIS history.  Averaging the already calibrated
!> CP2K AO transforms is the real-orthogonal projector for unitary and antiunitary operations.
!> \param kpoints ...
!> \param kmap ...
!> \param ao_first ...
!> \param ao_size ...
!> \param reduced ...
! **************************************************************************************************
   SUBROUTINE tb_project_gxtb_reduced_response(kpoints, kmap, ao_first, ao_size, reduced)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), INTENT(IN)           :: kmap
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_first, ao_size
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: reduced

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: projected, source, transformed
      INTEGER                                            :: ifull, ikred, iop

      CPASSERT(SIZE(reduced, 3) == kpoints%nkp)
      ALLOCATE (projected(SIZE(reduced, 1), SIZE(reduced, 2)), &
                source(SIZE(reduced, 1), SIZE(reduced, 2)), &
                transformed(SIZE(reduced, 1), SIZE(reduced, 2)))
      DO ikred = 1, kpoints%nkp
         ifull = kmap%representative(ikred)
         CPASSERT(ifull >= 1 .AND. ifull <= kmap%nfull)
         source(:, :) = reduced(:, :, ikred)
         projected = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO iop = 1, kmap%op_count(ifull)
            CALL tb_transform_gxtb_kpoint_matrix(kpoints, ikred, &
                                                 kmap%symops(iop, ifull), &
                                                 kmap%rotation_slot(iop, ifull), ao_first, ao_size, &
                                                 source, transformed, .FALSE.)
            projected(:, :) = projected + transformed/REAL(kmap%op_count(ifull), KIND=dp)
         END DO
         reduced(:, :, ikred) = 0.5_dp*(projected + CONJG(TRANSPOSE(projected)))
      END DO
      DEALLOCATE (projected, source, transformed)

   END SUBROUTINE tb_project_gxtb_reduced_response

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_ks ...
!> \param fmwork ...
!> Apply save_tblite's native potential mixer to the complete Gamma-point Fock matrices.
! **************************************************************************************************
   SUBROUTINE tb_mix_gxtb_gamma_fm(qs_env, matrix_ks, fmwork)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: fmwork

      INTEGER                                            :: ispin, nao, ndim, nspin
      REAL(KIND=dp)                                      :: hermitian_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mixed_vector, raw_vector
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: mixed_fock, raw_fock
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(tblite_type), POINTER                         :: tb

      NULLIFY (dft_control, para_env, scf_env, tb)
      IF (.NOT. ASSOCIATED(qs_env)) RETURN
      CALL get_qs_env(qs_env, dft_control=dft_control, para_env=para_env, scf_env=scf_env)
      tb => qs_env%tb_tblite
      IF (.NOT. ASSOCIATED(tb) .OR. .NOT. ASSOCIATED(dft_control) .OR. &
          .NOT. ASSOCIATED(para_env) .OR. .NOT. ASSOCIATED(scf_env)) RETURN
      IF (dft_control%qs_control%xtb_control%tblite_method /= gxtb) RETURN
      IF (.NOT. tb_native_scc_mixer_active(dft_control)) RETURN

      nao = tb%calc%bas%nao
      nspin = SIZE(matrix_ks)
      CPASSERT(SIZE(fmwork) >= nspin)
      CPASSERT(nspin == tb%wfn%nspin)
      ndim = nao*nao*nspin
      ALLOCATE (raw_fock(nao, nao, nspin), mixed_fock(nao, nao, nspin), &
                raw_vector(ndim), mixed_vector(ndim))
      DO ispin = 1, nspin
         CALL tb_dbcsr_to_dense(tb, matrix_ks(ispin)%matrix, raw_fock(:, :, ispin), para_env)
         hermitian_error = MAXVAL(ABS(raw_fock(:, :, ispin) - &
                                      TRANSPOSE(raw_fock(:, :, ispin))))
         IF (hermitian_error > 1.0E-10_dp) &
            CPABORT("CP2K supplied a non-symmetric Gamma Fock matrix to the g-xTB mixer")
      END DO

      raw_vector(:) = RESHAPE(raw_fock, [ndim])
      CALL tb_advance_gxtb_potential_mixer(tb, raw_vector, mixed_vector, scf_env%iter_count)
      mixed_fock(:, :, :) = RESHAPE(mixed_vector, SHAPE(mixed_fock))
      DO ispin = 1, nspin
         hermitian_error = MAXVAL(ABS(mixed_fock(:, :, ispin) - &
                                      TRANSPOSE(mixed_fock(:, :, ispin))))
         IF (hermitian_error > 1.0E-10_dp) &
            CPABORT("save_tblite returned a non-symmetric mixed Gamma Fock matrix")
         CALL tb_dense_to_fm(mixed_fock(:, :, ispin), fmwork(ispin))
      END DO
      DEALLOCATE (raw_fock, mixed_fock, raw_vector, mixed_vector)

   END SUBROUTINE tb_mix_gxtb_gamma_fm

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_ks ...
!> \param kpoints ...
!> \param scf_env ...
!> Build and mix all complete irreducible-k g-xTB Fock matrices as one weighted vector.
! **************************************************************************************************
   SUBROUTINE tb_prepare_gxtb_kpoint_mixer(qs_env, matrix_ks, kpoints, scf_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: mixed_fock_full, raw_fock
      INTEGER                                            :: block_size, ik, ispin, nao, nc, ndim, &
                                                            nkp, nspin, offset
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_first, ao_size
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: use_real_wfn
      REAL(KIND=dp)                                      :: hermitian_error, scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mixed_vector, raw_vector
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: fock_imag, fock_real
      TYPE(dbcsr_type), POINTER                          :: cmatrix, rmatrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(tb_gxtb_kmesh_map_type), POINTER              :: kmap
      TYPE(tblite_type), POINTER                         :: tb

      NULLIFY (cell_to_index, cmatrix, dft_control, kmap, para_env, rmatrix, sab_nl, tb)
      IF (.NOT. ASSOCIATED(qs_env) .OR. .NOT. ASSOCIATED(kpoints) .OR. &
          .NOT. ASSOCIATED(scf_env)) RETURN
      CALL get_qs_env(qs_env, dft_control=dft_control)
      tb => qs_env%tb_tblite
      IF (.NOT. ASSOCIATED(tb) .OR. .NOT. ASSOCIATED(dft_control)) RETURN
      tb%mixed_fock_kp_valid = .FALSE.
      tb%previous_density_kp_valid = .FALSE.
      IF (dft_control%qs_control%xtb_control%tblite_method /= gxtb) RETURN

      CALL get_kpoint_info(kpoints, nkp=nkp, use_real_wfn=use_real_wfn, &
                           cell_to_index=cell_to_index, sab_nl=sab_nl)
      CPASSERT(ASSOCIATED(sab_nl))
      IF (ABS(SUM(kpoints%wkp(1:nkp)) - 1.0_dp) > 1.0E-10_dp) &
         CPABORT("g-xTB k-point weights do not sum to one")
      IF (ANY(kpoints%wkp(1:nkp) <= 0.0_dp)) &
         CPABORT("g-xTB potential mixing requires positive k-point weights")
      IF (.NOT. tb_native_scc_mixer_active(dft_control)) THEN
         nao = tb%calc%bas%nao
         nspin = SIZE(matrix_ks, 1)
         IF (ALLOCATED(tb%previous_density_kp)) DEALLOCATE (tb%previous_density_kp)
         ALLOCATE (tb%previous_density_kp(nao, nao, nspin, nkp))
         CALL tb_collect_kpoint_density(qs_env, tb, kpoints, tb%previous_density_kp)
         tb%previous_density_kp_valid = .TRUE.
         RETURN
      END IF
      IF (ALLOCATED(tb%calc%exchange) .AND. .NOT. tb%exchange_fock_kp_valid) &
         CPABORT("Missing g-xTB exchange Fock matrices before potential mixing")
      IF (ALLOCATED(tb%calc%acp) .AND. .NOT. tb%acp_fock_kp_valid) &
         CPABORT("Missing direct k-space ACP matrices before potential mixing")
      nao = tb%calc%bas%nao
      nspin = SIZE(matrix_ks, 1)
      CPASSERT(nspin == tb%wfn%nspin)
      nc = 2
      IF (use_real_wfn) nc = 1
      block_size = nao*nao
      ndim = block_size*nspin*nkp*nc
      ALLOCATE (raw_fock(nao, nao, nspin, nkp), fock_real(nao, nao), &
                fock_imag(nao, nao), raw_vector(ndim), mixed_vector(ndim))

      para_env => kpoints%blacs_env_all%para_env
      ALLOCATE (rmatrix, cmatrix)
      CALL dbcsr_create(rmatrix, template=matrix_ks(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix, template=matrix_ks(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_nl)

      DO ik = 1, nkp
         DO ispin = 1, nspin
            CALL dbcsr_set(rmatrix, 0.0_dp)
            CALL dbcsr_set(cmatrix, 0.0_dp)
            CALL rskp_transform(rmatrix, cmatrix, matrix_ks, ispin, kpoints%xkp(1:3, ik), &
                                cell_to_index, sab_nl)
            CALL tb_dbcsr_component_to_dense(tb, rmatrix, fock_real, para_env, .FALSE.)
            IF (use_real_wfn) THEN
               fock_imag = 0.0_dp
            ELSE
               CALL tb_dbcsr_component_to_dense(tb, cmatrix, fock_imag, para_env, .TRUE.)
            END IF
            raw_fock(:, :, ispin, ik) = CMPLX(fock_real, fock_imag, KIND=dp)
            IF (tb%acp_fock_kp_valid) &
               raw_fock(:, :, ispin, ik) = raw_fock(:, :, ispin, ik) + &
                                           tb%acp_fock_kp(:, :, ik)
            IF (tb%exchange_fock_kp_valid) &
               raw_fock(:, :, ispin, ik) = raw_fock(:, :, ispin, ik) + &
                                           tb%exchange_fock_kp(:, :, ispin, ik)
            hermitian_error = MAXVAL(ABS(raw_fock(:, :, ispin, ik) - &
                                         CONJG(TRANSPOSE(raw_fock(:, :, ispin, ik)))))
            IF (hermitian_error > 1.0E-10_dp) &
               CPABORT("CP2K supplied a non-Hermitian k-point Fock matrix to the g-xTB mixer")
         END DO
      END DO
      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)

      IF (.NOT. ALLOCATED(tb%gxtb_kmesh_map)) &
         CPABORT("Missing complete g-xTB k-mesh map before potential mixing")
      kmap => tb%gxtb_kmesh_map
      IF (.NOT. kmap%identity) THEN
         CALL tb_gxtb_ao_layout(tb, ao_first, ao_size)
         DO ispin = 1, nspin
            CALL tb_project_gxtb_reduced_response(kpoints, kmap, ao_first, ao_size, &
                                                  raw_fock(:, :, ispin, :))
         END DO
      END IF

      offset = 0
      DO ik = 1, nkp
         scale = SQRT(REAL(nkp*nc, KIND=dp)*kpoints%wkp(ik))
         DO ispin = 1, nspin
            raw_vector(offset + 1:offset + block_size) = &
               scale*RESHAPE(REAL(raw_fock(:, :, ispin, ik), KIND=dp), [block_size])
            offset = offset + block_size
            IF (.NOT. use_real_wfn) THEN
               raw_vector(offset + 1:offset + block_size) = &
                  scale*RESHAPE(AIMAG(raw_fock(:, :, ispin, ik)), [block_size])
               offset = offset + block_size
            END IF
         END DO
      END DO
      CPASSERT(offset == ndim)
      CALL tb_advance_gxtb_potential_mixer(tb, raw_vector, mixed_vector, scf_env%iter_count)

      IF (ALLOCATED(tb%mixed_fock_kp)) DEALLOCATE (tb%mixed_fock_kp)
      ALLOCATE (tb%mixed_fock_kp(nao, nao, nspin, nkp), &
                source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      offset = 0
      DO ik = 1, nkp
         scale = SQRT(REAL(nkp*nc, KIND=dp)*kpoints%wkp(ik))
         DO ispin = 1, nspin
            fock_real(:, :) = RESHAPE(mixed_vector(offset + 1:offset + block_size)/scale, &
                                      SHAPE(fock_real))
            offset = offset + block_size
            IF (use_real_wfn) THEN
               fock_imag = 0.0_dp
            ELSE
               fock_imag(:, :) = RESHAPE(mixed_vector(offset + 1:offset + block_size)/scale, &
                                         SHAPE(fock_imag))
               offset = offset + block_size
            END IF
            tb%mixed_fock_kp(:, :, ispin, ik) = CMPLX(fock_real, fock_imag, KIND=dp)
            hermitian_error = MAXVAL(ABS(tb%mixed_fock_kp(:, :, ispin, ik) - &
                                         CONJG(TRANSPOSE(tb%mixed_fock_kp(:, :, ispin, ik)))))
            IF (hermitian_error > 1.0E-10_dp) &
               CPABORT("save_tblite returned a non-Hermitian mixed k-point Fock matrix")
         END DO
      END DO
      CPASSERT(offset == ndim)
      IF (.NOT. kmap%identity) THEN
         DO ispin = 1, nspin
            CALL tb_project_gxtb_reduced_response(kpoints, kmap, ao_first, ao_size, &
                                                  tb%mixed_fock_kp(:, :, ispin, :))
         END DO
         ALLOCATE (mixed_fock_full(nao, nao, nspin, kmap%nfull), &
                   source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))
         DO ik = 1, nkp
            mixed_fock_full(:, :, :, kmap%representative(ik)) = &
               tb%mixed_fock_kp(:, :, :, ik)
         END DO
         CALL tb_expand_gxtb_kpoint_density(kpoints, kmap, ao_first, ao_size, &
                                            mixed_fock_full, .FALSE.)
         DO ispin = 1, nspin
            CALL tb_validate_gxtb_folded_response(kpoints, kmap, ao_first, ao_size, &
                                                  mixed_fock_full(:, :, ispin, :), &
                                                  tb%mixed_fock_kp(:, :, ispin, :), &
                                                  "mixed total Fock")
         END DO
         DEALLOCATE (mixed_fock_full, ao_first, ao_size)
      END IF
      tb%mixed_fock_kp_valid = .TRUE.
      DEALLOCATE (raw_fock, fock_real, fock_imag, raw_vector, mixed_vector)

   END SUBROUTINE tb_prepare_gxtb_kpoint_mixer

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param diis_step ...
!> Keep P(k) consistent with CP2K's real-space DIRECT_P mixing for a non-native g-xTB backend.
!>
!> CP2K's standard density mixer updates only the real-space image matrices, while g-xTB exchange
!> consumes the irreducible-k kp%pmat matrices directly.  Mix the latter with the same scalar step
!> after each diagonalization so the next F(k) build sees exactly the selected CP2K density state.
! **************************************************************************************************
   SUBROUTINE tb_sync_gxtb_kpoint_density(qs_env, diis_step)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: diis_step

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: current_density, mixed_density
      INTEGER                                            :: ik, ikk, ispin, nao, nkp, nspin
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: use_real_wfn
      REAL(KIND=dp)                                      :: alpha, hermitian_error
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(tblite_type), POINTER                         :: tb

      NULLIFY (dft_control, kp, kpoints, scf_env, tb)
      IF (.NOT. ASSOCIATED(qs_env)) RETURN
      CALL get_qs_env(qs_env, dft_control=dft_control, kpoints=kpoints, scf_env=scf_env)
      tb => qs_env%tb_tblite
      IF (.NOT. ASSOCIATED(tb) .OR. .NOT. ASSOCIATED(dft_control) .OR. &
          .NOT. ASSOCIATED(kpoints) .OR. .NOT. ASSOCIATED(scf_env)) RETURN
      IF (dft_control%qs_control%xtb_control%tblite_method /= gxtb) RETURN
      IF (tb_native_scc_mixer_active(dft_control)) RETURN
      IF (scf_env%mixing_method /= direct_mixing_nr) &
         CPABORT("g-xTB with the CP2K/NONE backend currently requires DIRECT_P_MIXING")
      IF (.NOT. tb%previous_density_kp_valid .OR. &
          .NOT. ALLOCATED(tb%previous_density_kp)) &
         CPABORT("Missing input P(k) for CP2K g-xTB density mixing")

      CALL get_kpoint_info(kpoints, nkp=nkp, kp_range=kp_range, use_real_wfn=use_real_wfn)
      nao = tb%calc%bas%nao
      nspin = tb%wfn%nspin
      CPASSERT(ALL(SHAPE(tb%previous_density_kp) == [nao, nao, nspin, nkp]))
      ALLOCATE (current_density(nao, nao, nspin, nkp), mixed_density(nao, nao, nspin, nkp))
      CALL tb_collect_kpoint_density(qs_env, tb, kpoints, current_density)

      alpha = scf_env%mixing_store%alpha
      IF (diis_step .OR. scf_env%iter_count < scf_env%mixing_store%nskip_mixing) alpha = 1.0_dp
      IF (alpha < 0.0_dp .OR. alpha > 1.0_dp) &
         CPABORT("g-xTB k-point DIRECT_P mixing requires 0 <= ALPHA <= 1")
      mixed_density(:, :, :, :) = &
         alpha*current_density + (1.0_dp - alpha)*tb%previous_density_kp
      IF (.NOT. ALL(ieee_is_finite(REAL(mixed_density, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(mixed_density)))) &
         CPABORT("CP2K g-xTB k-point density mixing produced a non-finite P(k)")

      DO ik = kp_range(1), kp_range(2)
         ikk = ik - kp_range(1) + 1
         kp => kpoints%kp_env(ikk)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         DO ispin = 1, nspin
            hermitian_error = MAXVAL(ABS(mixed_density(:, :, ispin, ik) - &
                                         CONJG(TRANSPOSE(mixed_density(:, :, ispin, ik)))))
            IF (hermitian_error > 1.0E-10_dp) &
               CPABORT("CP2K g-xTB density mixing produced a non-Hermitian P(k)")
            CALL tb_dense_to_fm(REAL(mixed_density(:, :, ispin, ik), KIND=dp), &
                                kp%pmat(1, ispin))
            IF (use_real_wfn) THEN
               IF (MAXVAL(ABS(AIMAG(mixed_density(:, :, ispin, ik)))) > 1.0E-10_dp) &
                  CPABORT("Complex mixed g-xTB P(k) cannot use real k-point wavefunctions")
            ELSE
               CALL tb_dense_to_fm(AIMAG(mixed_density(:, :, ispin, ik)), kp%pmat(2, ispin))
            END IF
         END DO
      END DO

      tb%previous_density_kp(:, :, :, :) = mixed_density(:, :, :, :)
      tb%previous_density_kp_valid = .FALSE.
      DEALLOCATE (current_density, mixed_density)

   END SUBROUTINE tb_sync_gxtb_kpoint_density

! **************************************************************************************************
!> \brief Build the geometry-static separable ACP directly on the complete k mesh.
!> \param qs_env ...
!> \param tb ...
!> \param matrix_s ...
!> \param sab_nl ...
!> \param cell_to_index ...
!> \param projectors compact save_tblite projector images (no H(R) difference set)
!>
!> The full-mesh Bloch operator is folded with the real adjoint of CP2K's density expansion.
!> It is therefore an unweighted irreducible-k response and can be added directly to every spin
!> Fock matrix.  No Born-von-Karman modulo folding or real-space image truncation is involved.
! **************************************************************************************************
   SUBROUTINE tb_build_gxtb_kpoint_acp(qs_env, tb, matrix_s, sab_nl, &
                                       cell_to_index, projectors)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(acp_image_type), INTENT(IN)                   :: projectors

      CHARACTER(LEN=32)                                  :: qualification_value
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: acp_full
      INTEGER                                            :: ik, nao, nred, qualification_status
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_first, ao_size
      LOGICAL                                            :: use_real_wfn
      REAL(KIND=dp)                                      :: hermitian_error
      TYPE(error_type), ALLOCATABLE                      :: error
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(tb_gxtb_kmesh_map_type), POINTER              :: kmap

      NULLIFY (kmap, kpoints)
      tb%acp_fock_kp_valid = .FALSE.
      tb%e_acp_kp = 0.0_dp
      IF (ALLOCATED(tb%acp_fock_kp)) DEALLOCATE (tb%acp_fock_kp)
      CALL get_qs_env(qs_env, kpoints=kpoints)
      CPASSERT(ASSOCIATED(kpoints))
      CALL get_kpoint_info(kpoints, nkp=nred, use_real_wfn=use_real_wfn)
      CALL tb_prepare_gxtb_kmesh_map(tb, kpoints, matrix_s, sab_nl, cell_to_index)
      kmap => tb%gxtb_kmesh_map
      nao = tb%calc%bas%nao
      CALL tb_gxtb_ao_layout(tb, ao_first, ao_size)
      ALLOCATE (acp_full(nao, nao, kmap%nfull), &
                tb%acp_fock_kp(nao, nao, nred), &
                source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))

      CALL cp2k_acp_kmesh(tb%calc, tb%mol, projectors, kmap%xkp, acp_full, error)
      IF (ALLOCATED(error)) &
         CPABORT("save_tblite direct k-space ACP evaluation failed: "//TRIM(error%message))
      IF (.NOT. ALL(ieee_is_finite(REAL(acp_full, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(acp_full)))) &
         CPABORT("save_tblite returned a non-finite k-space ACP matrix")
      DO ik = 1, kmap%nfull
         hermitian_error = MAXVAL(ABS(acp_full(:, :, ik) - &
                                      CONJG(TRANSPOSE(acp_full(:, :, ik)))))
         IF (hermitian_error > 1.0E-10_dp) &
            CPABORT("save_tblite returned a non-Hermitian k-space ACP matrix")
      END DO

      IF (.NOT. kmap%identity) &
         CALL tb_validate_gxtb_full_response(kpoints, kmap, ao_first, ao_size, &
                                             acp_full, "ACP")
      CALL tb_fold_gxtb_kpoint_matrix(kpoints, kmap, ao_first, ao_size, &
                                      acp_full, tb%acp_fock_kp)
      IF (.NOT. kmap%identity) &
         CALL tb_validate_gxtb_folded_response(kpoints, kmap, ao_first, ao_size, &
                                               acp_full, tb%acp_fock_kp, "ACP")
      qualification_value = ""
      qualification_status = 1
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION", &
                                    qualification_value, STATUS=qualification_status)
      IF (ALLOCATED(tb%qualification_acp_fock_full)) &
         DEALLOCATE (tb%qualification_acp_fock_full)
      IF (qualification_status == 0) THEN
         ALLOCATE (tb%qualification_acp_fock_full(nao, nao, kmap%nfull))
         tb%qualification_acp_fock_full(:, :, :) = acp_full(:, :, :)
      END IF
      DO ik = 1, nred
         hermitian_error = MAXVAL(ABS(tb%acp_fock_kp(:, :, ik) - &
                                      CONJG(TRANSPOSE(tb%acp_fock_kp(:, :, ik)))))
         IF (hermitian_error > 1.0E-10_dp) &
            CPABORT("Folded g-xTB ACP matrix is not Hermitian")
         IF (use_real_wfn .AND. MAXVAL(ABS(AIMAG(tb%acp_fock_kp(:, :, ik)))) > 1.0E-10_dp) &
            CPABORT("Complex g-xTB ACP requires complex k-point wavefunctions")
      END DO
      tb%acp_fock_kp_valid = .TRUE.
      DEALLOCATE (acp_full, ao_first, ao_size)

   END SUBROUTINE tb_build_gxtb_kpoint_acp

! **************************************************************************************************
!> \brief Contract the current irreducible density with the folded k-space ACP.
!> \param qs_env ...
!> \param tb ...
!> \param kpoints ...
!>
!> The usual real-space core-energy contraction no longer contains ACP once H(R) scattering is
!> removed.  This exact Bloch contraction restores that one-electron energy without changing the
!> Fock operator or applying a mesh weight twice.
! **************************************************************************************************
   SUBROUTINE tb_update_gxtb_kpoint_acp_energy(qs_env, tb, kpoints)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      TYPE(kpoint_type), POINTER                         :: kpoints

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: density
      INTEGER                                            :: ik, ispin, nao, nkp, nspin

      IF (.NOT. tb%acp_fock_kp_valid .OR. .NOT. ALLOCATED(tb%acp_fock_kp)) &
         CPABORT("Missing direct k-space ACP before its energy contraction")
      CALL get_kpoint_info(kpoints, nkp=nkp)
      nao = tb%calc%bas%nao
      nspin = tb%wfn%nspin
      CPASSERT(ALL(SHAPE(tb%acp_fock_kp) == [nao, nao, nkp]))
      ALLOCATE (density(nao, nao, nspin, nkp))
      CALL tb_collect_kpoint_density(qs_env, tb, kpoints, density)
      IF (.NOT. ALL(ieee_is_finite(REAL(density, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(density)))) &
         CPABORT("CP2K supplied a non-finite density to the k-space ACP energy")
      tb%e_acp_kp = 0.0_dp
      DO ik = 1, nkp
         DO ispin = 1, nspin
            tb%e_acp_kp = tb%e_acp_kp + kpoints%wkp(ik)* &
                          REAL(SUM(CONJG(density(:, :, ispin, ik))*tb%acp_fock_kp(:, :, ik)), KIND=dp)
         END DO
      END DO
      IF (.NOT. ieee_is_finite(tb%e_acp_kp)) &
         CPABORT("Non-finite direct k-space ACP energy")
      DEALLOCATE (density)

   END SUBROUTINE tb_update_gxtb_kpoint_acp_energy

! **************************************************************************************************
!> \brief Select the qualification-only CP2K bridge for streamed g-xTB exchange.
!>
!> Zero retains the Part-I complete-array implementation.  The nonzero values expose the
!> provider's matrix-lean, full-mesh, bounded image-batch, and MPI-distributed image-range
!> transaction modes for end-to-end qualification without prematurely changing the public CP2K
!> input contract.  Once the complete oracle matrix has passed, a qualified reduced path can
!> become input-selectable.
!> \return selected provider stream mode, or zero for the complete-array bridge
! **************************************************************************************************
   INTEGER FUNCTION tb_gxtb_exchange_stream_mode() RESULT(stream_mode)

      CHARACTER(LEN=32)                                  :: value
      INTEGER                                            :: status

      stream_mode = 0
      status = 1
      value = ""
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_GXTB_EXCHANGE_STREAM_MODE", value, STATUS=status)
      IF (status /= 0 .OR. LEN_TRIM(value) == 0) RETURN

      SELECT CASE (TRIM(ADJUSTL(value)))
      CASE ("LEGACY", "legacy", "0")
         stream_mode = 0
      CASE ("REDUCED", "reduced", "1")
         stream_mode = cp2k_exchange_stream_reduced
      CASE ("ORACLE", "oracle", "2")
         stream_mode = cp2k_exchange_stream_oracle
      CASE ("BATCHED", "batched", "3")
         stream_mode = cp2k_exchange_stream_batched
      CASE ("DISTRIBUTED", "distributed", "MPI", "mpi", "4")
         stream_mode = tb_gxtb_exchange_stream_distributed
      CASE ("KGROUP_OWNER", "kgroup_owner", "K_GROUP_OWNER", "k_group_owner", "5")
         stream_mode = tb_gxtb_exchange_stream_kgroup_owner
      CASE DEFAULT
         CPABORT("Unknown CP2K_GXTB_EXCHANGE_STREAM_MODE value")
      END SELECT

   END FUNCTION tb_gxtb_exchange_stream_mode

! **************************************************************************************************
!> \brief Return the qualification-only bounded image batch size.
!>
!> The default keeps eight BvK image representatives at a time.  A positive environment override
!> permits memory/scaling qualification without extending the public input contract prematurely.
!> \return requested positive image batch size
! **************************************************************************************************
   INTEGER FUNCTION tb_gxtb_exchange_stream_batch_size() RESULT(batch_size)

      CHARACTER(LEN=32)                                  :: value
      INTEGER                                            :: status

      batch_size = 8
      status = 1
      value = ""
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE", value, STATUS=status)
      IF (status /= 0 .OR. LEN_TRIM(value) == 0) RETURN
      READ (value, *, IOSTAT=status) batch_size
      IF (status /= 0 .OR. batch_size < 1) &
         CPABORT("CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE must be a positive integer")

   END FUNCTION tb_gxtb_exchange_stream_batch_size

! **************************************************************************************************
!> \brief Select the CP2K consumer used for the coupled-exchange derivative.
!>
!> DENSE preserves the Part-I complete-mesh implementation.  STREAMED uses the
!> provider's bounded reverse-only transaction and folds each returned full-star
!> overlap adjoint immediately.  QUALIFY executes STREAMED, then constructs the
!> unchanged complete-mesh derivative only as a numerical oracle before applying
!> the streamed response exactly once.
!> \return selected derivative bridge
! **************************************************************************************************
   INTEGER FUNCTION tb_gxtb_exchange_gradient_mode() RESULT(gradient_mode)

      CHARACTER(LEN=32)                                  :: value
      INTEGER                                            :: status

      gradient_mode = tb_gxtb_exchange_gradient_dense
      status = 1
      value = ""
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_GXTB_EXCHANGE_GRADIENT_MODE", value, STATUS=status)
      IF (status /= 0 .OR. LEN_TRIM(value) == 0) RETURN

      SELECT CASE (TRIM(ADJUSTL(value)))
      CASE ("DENSE", "dense", "LEGACY", "legacy", "0")
         gradient_mode = tb_gxtb_exchange_gradient_dense
      CASE ("STREAMED", "streamed", "BATCHED", "batched", "1")
         gradient_mode = tb_gxtb_exchange_gradient_streamed
      CASE ("QUALIFY", "qualify", "ORACLE", "oracle", "2")
         gradient_mode = tb_gxtb_exchange_gradient_qualify
      CASE DEFAULT
         CPABORT("Unknown CP2K_GXTB_EXCHANGE_GRADIENT_MODE value")
      END SELECT

   END FUNCTION tb_gxtb_exchange_gradient_mode

! **************************************************************************************************
!> \brief Select the regular-mesh Fourier backend used by complete-mesh g-xTB exchange.
!>
!> Dense phase-table multiplication remains the default and permanent numerical oracle.  The
!> separable backend is a factorized direct DFT, not an FFT: it applies one-dimensional transforms
!> along each regular mesh direction and therefore scales as Nk*sum(nmesh) character products.
!> The provider rejects meshes or weights that are not machine-regular instead of silently falling
!> back, so qualification runs always exercise the backend that was requested.
!> \return provider transform-backend identifier
! **************************************************************************************************
   INTEGER FUNCTION tb_gxtb_exchange_transform_backend() RESULT(transform_backend)

      CHARACTER(LEN=32)                                  :: value
      INTEGER                                            :: status

      transform_backend = cp2k_exchange_transform_dense
      status = 1
      value = ""
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_GXTB_EXCHANGE_TRANSFORM_MODE", value, STATUS=status)
      IF (status /= 0 .OR. LEN_TRIM(value) == 0) RETURN

      SELECT CASE (TRIM(ADJUSTL(value)))
      CASE ("DENSE", "dense", "1")
         transform_backend = cp2k_exchange_transform_dense
      CASE ("SEPARABLE", "separable", "SEPARABLE_DFT", "separable_dft", "2")
         transform_backend = cp2k_exchange_transform_separable_dft
      CASE DEFAULT
         CPABORT("Unknown CP2K_GXTB_EXCHANGE_TRANSFORM_MODE value")
      END SELECT

   END FUNCTION tb_gxtb_exchange_transform_backend

! **************************************************************************************************
!> \brief Exact owner-generated g-xTB exchange bridge using CP2K k-point groups.
!> \param qs_env ...
!> \param tb ...
!> \param matrix_s ...
!> \param run_qualification_oracle ...
!> \param iter_count ...
!>
!> This qualification path removes the globally replicated irreducible-density
!> array from the streamed bridge.  Every CP2K k-point group gathers only its
!> own P(k).  Its group leader expands exactly the full-mesh star members whose
!> irreducible source lies in kp_range; para_env_inter_kp transfers one block at
!> a time to the unique provider rank.  The provider state and coupled kernel
!> are intentionally nonreplicated.  Returned full-k Fock blocks travel only
!> among group leaders, are folded by the irreducible owner, and are finally
!> reduced over para_env_inter_kp and broadcast inside para_env_kp.
!>
!> The current public save_tblite stream API requires every Bloch block in one
!> stream before apply and does not expose mergeable partial k-to-R
!> accumulators.  Consequently this exact precursor demonstrates CP2K
!> ownership and bounded block communication, but does not claim parallel
!> provider-kernel scaling.  A future distributed variant must use a public
!> provider accumulator create/push/merge/finalize API, never private stream
!> components.
! **************************************************************************************************
   SUBROUTINE tb_build_gxtb_kpoint_exchange_kgroup_owner(qs_env, tb, matrix_s, &
                                                         run_qualification_oracle, iter_count)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      LOGICAL, INTENT(IN)                                :: run_qualification_oracle
      INTEGER, INTENT(IN)                                :: iter_count

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: overlap_block, transform_work
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: density_block, fock_block
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: density_local
      INTEGER                                            :: group_index, ifull, iop, ispin, &
                                                            local_source, nao, nkp_groups, &
                                                            nlocal, nred, nspin, owner_count, &
                                                            source
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_first, ao_size
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: owns_source, producer, service_rank, &
                                                            use_real_wfn, validate_density
      REAL(KIND=dp)                                      :: density_norm, exchange_energy, &
                                                            hermitian_error, residual, scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: exchange_vsh
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: overlap_part
      TYPE(cp2k_exchange_stream_type)                    :: stream
      TYPE(dbcsr_type), POINTER                          :: cmatrix, rmatrix
      TYPE(error_type), ALLOCATABLE                      :: error
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env_global, para_env_inter_kp, &
                                                            para_env_kp
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(tb_gxtb_kmesh_map_type), POINTER              :: kmap

      NULLIFY (cell_to_index, cmatrix, kmap, kp_dist, kpoints, para_env_global, &
               para_env_inter_kp, para_env_kp, rmatrix, sab_nl)
      CALL get_qs_env(qs_env, kpoints=kpoints)
      CPASSERT(ASSOCIATED(kpoints))
      CALL get_kpoint_info(kpoints, nkp=nred, kp_range=kp_range, nkp_groups=nkp_groups, &
                           kp_dist=kp_dist, para_env_kp=para_env_kp, &
                           para_env_inter_kp=para_env_inter_kp, &
                           para_env=para_env_global, &
                           use_real_wfn=use_real_wfn, cell_to_index=cell_to_index, &
                           sab_nl=sab_nl)
      IF (.NOT. ASSOCIATED(para_env_global) .OR. .NOT. ASSOCIATED(para_env_kp) .OR. &
          .NOT. ASSOCIATED(para_env_inter_kp)) &
         CPABORT("Missing official CP2K k-point communicators for KGROUP_OWNER")
      IF (.NOT. ASSOCIATED(kp_dist) .OR. .NOT. ASSOCIATED(cell_to_index) .OR. &
          .NOT. ASSOCIATED(sab_nl)) &
         CPABORT("Incomplete CP2K k-point ownership or image map for KGROUP_OWNER")
      IF (nkp_groups /= para_env_inter_kp%num_pe) &
         CPABORT("CP2K inter-k communicator is inconsistent with nkp_groups")
      group_index = para_env_inter_kp%mepos + 1
      IF (group_index < 1 .OR. group_index > nkp_groups) &
         CPABORT("Invalid CP2K k-point group index")
      IF (ANY(kp_range /= kp_dist(:, group_index))) &
         CPABORT("CP2K k-point range differs from its official group distribution")
      nlocal = kp_range(2) - kp_range(1) + 1
      IF (nlocal < 1) CPABORT("KGROUP_OWNER requires at least one irreducible k point per group")

      CPASSERT(ALLOCATED(tb%gxtb_kmesh_map))
      kmap => tb%gxtb_kmesh_map
      CPASSERT(kmap%nred == nred)
      nao = tb%calc%bas%nao
      nspin = tb%wfn%nspin
      CALL tb_gxtb_ao_layout(tb, ao_first, ao_size)
      ALLOCATE (density_local(nao, nao, nspin, nlocal), &
                density_block(nao, nao, nspin), fock_block(nao, nao, nspin), &
                overlap_block(nao, nao), transform_work(nao, nao), &
                overlap_part(nao, nao), exchange_vsh(tb%calc%bas%nsh))
      CALL tb_collect_local_kpoint_density(tb, kpoints, density_local)
      IF (.NOT. ALL(ieee_is_finite(REAL(density_local, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(density_local)))) &
         CPABORT("CP2K supplied a non-finite owner-local density to KGROUP_OWNER")
      IF (.NOT. ALL(ieee_is_finite(tb%wfn%qsh))) &
         CPABORT("CP2K supplied non-finite shell charges to KGROUP_OWNER")

      density_norm = MAXVAL(ABS(density_local))
      CALL para_env_global%max(density_norm)
      validate_density = .NOT. kmap%identity .AND. &
                         .NOT. tb%gxtb_density_map_validated .AND. &
                         density_norm > 1.0E-14_dp
      service_rank = para_env_global%is_source()
      exchange_energy = 0.0_dp
      exchange_vsh = 0.0_dp
      ALLOCATE (tb%exchange_fock_kp(nao, nao, nspin, nred), &
                source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))

      ALLOCATE (rmatrix, cmatrix)
      CALL dbcsr_create(rmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_nl)

      IF (service_rank) THEN
         CALL cp2k_exchange_stream_begin(stream, tb%calc, tb%mol, tb%ecache, tb%wfn, &
                                         kmap%nmesh, kmap%xkp, kmap%weight, error, &
                                         mode=cp2k_exchange_stream_reduced)
      ELSE
         ! Every replicated CP2K model copy needs the current geometry- and
         ! charge-dependent derivative cache although only one stream owns the
         ! forward coupled kernel.
         CALL cp2k_prepare_exchange(tb%calc, tb%mol, tb%ecache, tb%wfn, error)
      END IF
      IF (ALLOCATED(error)) CPABORT("save_tblite KGROUP_OWNER initialization failed")

      DO ifull = 1, kmap%nfull
         source = kmap%source(ifull)
         owns_source = source >= kp_range(1) .AND. source <= kp_range(2)
         producer = owns_source .AND. para_env_kp%is_source()
         density_block = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         owner_count = 0
         IF (producer) THEN
            local_source = source - kp_range(1) + 1
            CPASSERT(local_source >= 1 .AND. local_source <= nlocal)
            DO ispin = 1, nspin
               CALL tb_expand_gxtb_kpoint_block( &
                  kpoints, kmap, ao_first, ao_size, ifull, &
                  density_local(:, :, ispin, local_source), &
                  density_block(:, :, ispin), transform_work)
               IF (validate_density) THEN
                  DO iop = 1, kmap%op_count(ifull)
                     CALL tb_transform_gxtb_kpoint_matrix( &
                        kpoints, source, kmap%symops(iop, ifull), &
                        kmap%rotation_slot(iop, ifull), ao_first, ao_size, &
                        density_local(:, :, ispin, local_source), transform_work, .FALSE.)
                     scale = MAX(1.0_dp, MAXVAL(ABS(density_block(:, :, ispin))))
                     residual = MAXVAL(ABS(transform_work - density_block(:, :, ispin)))/scale
                     IF (residual > MAX(1.0E-10_dp, 100.0_dp*kpoints%eps_geo)) &
                        CPABORT("KGROUP_OWNER density is not invariant under its CP2K little group")
                  END DO
               END IF
            END DO
            owner_count = 1
         END IF
         ! Only the communicator containing the k-group leaders carries the
         ! block.  Other inter-k communicators execute the same collective on
         ! zeros, avoiding assumptions about underlying MPI handles.
         CALL para_env_inter_kp%sum(owner_count)
         CALL para_env_inter_kp%sum(density_block)
         IF (para_env_kp%is_source() .AND. owner_count /= 1) &
            CPABORT("KGROUP_OWNER star block does not have exactly one owner")
         IF (service_rank) THEN
            hermitian_error = 0.0_dp
            DO ispin = 1, nspin
               hermitian_error = MAX(hermitian_error, &
                                     MAXVAL(ABS(density_block(:, :, ispin) - &
                                                CONJG(TRANSPOSE(density_block(:, :, ispin))))))
            END DO
            IF (hermitian_error > 1.0E-10_dp) &
               CPABORT("KGROUP_OWNER transferred a non-Hermitian density block")
         END IF

         CALL prepare_overlap_block(ifull)
         IF (service_rank) THEN
            CALL cp2k_exchange_stream_push(stream, ifull, density_block, overlap_block, error)
            IF (ALLOCATED(error)) CPABORT("save_tblite rejected a KGROUP_OWNER Bloch block")
         END IF
      END DO
      IF (validate_density) tb%gxtb_density_map_validated = .TRUE.

      IF (service_rank) THEN
         CALL cp2k_exchange_stream_apply(stream, tb%calc, tb%mol, tb%ecache, &
                                         exchange_vsh, exchange_energy, error)
         IF (ALLOCATED(error)) CPABORT("save_tblite KGROUP_OWNER apply failed")
      END IF
      CALL para_env_inter_kp%sum(exchange_energy)
      CALL para_env_inter_kp%sum(exchange_vsh)
      CALL para_env_kp%sum(exchange_energy)
      CALL para_env_kp%sum(exchange_vsh)

      DO ifull = 1, kmap%nfull
         source = kmap%source(ifull)
         owns_source = source >= kp_range(1) .AND. source <= kp_range(2)
         CALL prepare_overlap_block(ifull)
         fock_block = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         IF (service_rank) THEN
            CALL cp2k_exchange_stream_get(stream, ifull, fock_block, error, overlap_block)
            IF (ALLOCATED(error)) CPABORT("save_tblite KGROUP_OWNER Fock retrieval failed")
            IF (.NOT. ALL(ieee_is_finite(REAL(fock_block, KIND=dp))) .OR. &
                .NOT. ALL(ieee_is_finite(AIMAG(fock_block)))) &
               CPABORT("save_tblite returned a non-finite KGROUP_OWNER Fock block")
         END IF
         CALL para_env_inter_kp%sum(fock_block)
         IF (owns_source .AND. para_env_kp%is_source()) THEN
            DO ispin = 1, nspin
               CALL tb_accumulate_gxtb_kpoint_response( &
                  kpoints, kmap, ao_first, ao_size, ifull, fock_block(:, :, ispin), &
                  tb%exchange_fock_kp(:, :, ispin, source), transform_work)
            END DO
         END IF
      END DO
      IF (service_rank) THEN
         CALL cp2k_exchange_stream_end(stream, error)
         IF (ALLOCATED(error)) CPABORT("save_tblite KGROUP_OWNER finalization failed")
      END IF

      ! The leader-only owner pieces become a complete irreducible response on
      ! the leader inter-k communicator, then para_env_kp broadcasts that
      ! result within each CP2K k-point group.
      CALL para_env_inter_kp%sum(tb%exchange_fock_kp)
      CALL para_env_kp%sum(tb%exchange_fock_kp)
      DO source = 1, nred
         DO ispin = 1, nspin
            hermitian_error = &
               MAXVAL(ABS(tb%exchange_fock_kp(:, :, ispin, source) - &
                          CONJG(TRANSPOSE(tb%exchange_fock_kp(:, :, ispin, source)))))
            IF (hermitian_error > 1.0E-10_dp) &
               CPABORT("Folded KGROUP_OWNER exchange matrix is not Hermitian")
            tb%exchange_fock_kp(:, :, ispin, source) = &
               0.5_dp*(tb%exchange_fock_kp(:, :, ispin, source) + &
                       CONJG(TRANSPOSE(tb%exchange_fock_kp(:, :, ispin, source))))
         END DO
         IF (use_real_wfn .AND. &
             MAXVAL(ABS(AIMAG(tb%exchange_fock_kp(:, :, :, source)))) > 1.0E-10_dp) &
            CPABORT("Complex KGROUP_OWNER exchange cannot use real wavefunctions")
      END DO

      tb%e_exchange = 0.0_dp
      tb%e_exchange(1) = exchange_energy
      tb%pot%vsh(:, 1) = tb%pot%vsh(:, 1) + exchange_vsh
      tb%exchange_fock_kp_valid = .TRUE.
      IF (run_qualification_oracle) CALL qualify_kgroup_oracle()
      IF (service_rank .AND. iter_count == 1) THEN
         WRITE (*, '(A,I0,A,I0,A,I0,A)') &
            " GXTB-KGROUP-OWNER groups=", nkp_groups, ", nred=", nred, &
            ", nfull=", kmap%nfull, &
            "; provider kernel remains single-state until public accumulator merge ABI"
      END IF

      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)
      DEALLOCATE (ao_first, ao_size, density_block, density_local, exchange_vsh, &
                  fock_block, overlap_block, overlap_part, transform_work)

   CONTAINS

! **************************************************************************************************
!> \brief Build one replicated S(k) block through the unchanged CP2K image transform.
!> \param ifull_local ...
! **************************************************************************************************
      SUBROUTINE prepare_overlap_block(ifull_local)

      INTEGER, INTENT(IN)                                :: ifull_local

      REAL(KIND=dp)                                      :: local_hermitian_error

         CALL dbcsr_set(rmatrix, 0.0_dp)
         CALL dbcsr_set(cmatrix, 0.0_dp)
         CALL rskp_transform(rmatrix, cmatrix, matrix_s, 1, kmap%xkp(1:3, ifull_local), &
                             cell_to_index, sab_nl)
         CALL tb_dbcsr_component_to_dense(tb, rmatrix, overlap_part, para_env_global, .FALSE.)
         overlap_block(:, :) = CMPLX(overlap_part, 0.0_dp, KIND=dp)
         CALL tb_dbcsr_component_to_dense(tb, cmatrix, overlap_part, para_env_global, .TRUE.)
         overlap_block(:, :) = CMPLX(REAL(overlap_block, KIND=dp), overlap_part, KIND=dp)
         IF (.NOT. ALL(ieee_is_finite(REAL(overlap_block, KIND=dp))) .OR. &
             .NOT. ALL(ieee_is_finite(AIMAG(overlap_block)))) &
            CPABORT("CP2K supplied a non-finite overlap to KGROUP_OWNER")
         local_hermitian_error = &
            MAXVAL(ABS(overlap_block - CONJG(TRANSPOSE(overlap_block))))
         IF (local_hermitian_error > 1.0E-10_dp) &
            CPABORT("CP2K supplied a non-Hermitian overlap to KGROUP_OWNER")

      END SUBROUTINE prepare_overlap_block

! **************************************************************************************************
!> \brief Qualification-only complete-array dense oracle for energy and folded Fock.
! **************************************************************************************************
      SUBROUTINE qualify_kgroup_oracle()

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: overlap_full
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: density_full, fock_full, &
                                                            fock_reduced_oracle
      INTEGER                                            :: ifull_local, ispin_local, &
                                                            local_source_local, source_local
      REAL(KIND=dp)                                      :: energy_delta, energy_oracle, fock_delta, &
                                                            vsh_delta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vsh_oracle

         IF (service_rank) THEN
            ALLOCATE (density_full(nao, nao, nspin, kmap%nfull), &
                      fock_full(nao, nao, nspin, kmap%nfull), &
                      fock_reduced_oracle(nao, nao, nspin, nred), &
                      overlap_full(nao, nao, kmap%nfull), &
                      source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))
            ALLOCATE (vsh_oracle(SIZE(exchange_vsh)), source=0.0_dp)
         END IF
         DO ifull_local = 1, kmap%nfull
            source_local = kmap%source(ifull_local)
            owns_source = source_local >= kp_range(1) .AND. source_local <= kp_range(2)
            producer = owns_source .AND. para_env_kp%is_source()
            density_block = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            IF (producer) THEN
               local_source_local = source_local - kp_range(1) + 1
               DO ispin_local = 1, nspin
                  CALL tb_expand_gxtb_kpoint_block( &
                     kpoints, kmap, ao_first, ao_size, ifull_local, &
                     density_local(:, :, ispin_local, local_source_local), &
                     density_block(:, :, ispin_local), transform_work)
               END DO
            END IF
            CALL para_env_inter_kp%sum(density_block)
            CALL prepare_overlap_block(ifull_local)
            IF (service_rank) THEN
               density_full(:, :, :, ifull_local) = density_block
               overlap_full(:, :, ifull_local) = overlap_block
            END IF
         END DO
         IF (service_rank) THEN
            CALL cp2k_exchange_kmesh(tb%calc, tb%mol, tb%ecache, tb%wfn, kmap%nmesh, &
                                     kmap%xkp, kmap%weight, density_full, overlap_full, &
                                     fock_full, vsh_oracle, energy_oracle, error, &
                                     transform_backend=cp2k_exchange_transform_dense)
            IF (ALLOCATED(error)) CPABORT("KGROUP_OWNER dense qualification oracle failed")
            fock_reduced_oracle = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            DO ispin_local = 1, nspin
               CALL tb_fold_gxtb_kpoint_matrix( &
                  kpoints, kmap, ao_first, ao_size, fock_full(:, :, ispin_local, :), &
                  fock_reduced_oracle(:, :, ispin_local, :))
            END DO
            energy_delta = ABS(exchange_energy - energy_oracle)
            vsh_delta = MAXVAL(ABS(exchange_vsh - vsh_oracle))
            fock_delta = MAXVAL(ABS(tb%exchange_fock_kp - fock_reduced_oracle))
            WRITE (*, '(A,I0,3(A,ES14.6))') &
               " GXTB-QUALIFICATION_ONLY KGROUP-ORACLE iter=", iter_count, &
               " dE=", energy_delta, " dVsh=", vsh_delta, " dFfold=", fock_delta
            IF (energy_delta > 1.0E-10_dp .OR. vsh_delta > 1.0E-10_dp .OR. &
                fock_delta > 1.0E-10_dp) &
               CPABORT("KGROUP_OWNER dense qualification oracle exceeded 1.0E-10")
            DEALLOCATE (density_full, fock_full, fock_reduced_oracle, overlap_full, vsh_oracle)
         END IF

      END SUBROUTINE qualify_kgroup_oracle

   END SUBROUTINE tb_build_gxtb_kpoint_exchange_kgroup_owner

! **************************************************************************************************
!> \brief Evaluate complete-mesh g-xTB exchange without extra expanded full-star host arrays.
!> \param qs_env ...
!> \param tb ...
!> \param matrix_s ...
!> \param para_env ...
!> \param stream_mode provider reduced, oracle, bounded image-batch, or MPI image-range mode
!>
!> Irreducible CP2K densities are expanded one star member at a time and pushed sequentially into
!> save_tblite; one stream object is never entered concurrently by k-point groups or threads.
!> After the coupled contraction, one unweighted Fock block at a time is reconstructed and folded
!> immediately with the weighted real adjoint.  In bounded image-batch mode, CP2K regenerates the
!> Bloch inputs for every batch and folds each additive Fock contribution before the batch storage
!> is reused.  MPI image-range mode gives every active rank exactly one disjoint contiguous image
!> tile, then sums energy, shell potential, and the already folded irreducible Fock response.  All
!> ranks still enter the CP2K dense-matrix collectives in identical order.  The provider oracle mode
!> retains its complete arrays internally as a reference.
!> \param run_qualification_oracle ...
!> \param iter_count ...
! **************************************************************************************************
   SUBROUTINE tb_build_gxtb_kpoint_exchange_streamed(qs_env, tb, matrix_s, para_env, stream_mode, &
                                                     run_qualification_oracle, iter_count)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: stream_mode
      LOGICAL, INTENT(IN)                                :: run_qualification_oracle
      INTEGER, INTENT(IN)                                :: iter_count

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: overlap_block, transform_work
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: density_block, fock_block
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: density_reduced
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: overlap_full_qualification
      INTEGER                                            :: batch_size, first_image, ifull, ispin, &
                                                            last_image, local_image_count, &
                                                            local_image_first, nao, nowners, nred, &
                                                            nspin, source
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_first, ao_size, image_owner_count
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: batches_done, distributed_mode, &
                                                            expected_full_mesh_storage, &
                                                            has_full_kspace_inputs, &
                                                            owns_image_range, &
                                                            run_block_qualification, use_real_wfn, &
                                                            validate_density
      REAL(KIND=dp)                                      :: exchange_energy, hermitian_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: exchange_vsh
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: overlap_part
      TYPE(cp2k_exchange_stream_type)                    :: stream
      TYPE(dbcsr_type), POINTER                          :: cmatrix, rmatrix
      TYPE(error_type), ALLOCATABLE                      :: error
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(tb_gxtb_kmesh_map_type), POINTER              :: kmap

      NULLIFY (cell_to_index, cmatrix, kmap, kpoints, rmatrix, sab_nl)
      IF (stream_mode /= cp2k_exchange_stream_reduced .AND. &
          stream_mode /= cp2k_exchange_stream_oracle .AND. &
          stream_mode /= cp2k_exchange_stream_batched .AND. &
          stream_mode /= tb_gxtb_exchange_stream_distributed) &
         CPABORT("Invalid save_tblite exchange stream mode")
      CALL get_qs_env(qs_env, kpoints=kpoints)
      CPASSERT(ASSOCIATED(kpoints))
      CALL get_kpoint_info(kpoints, nkp=nred, use_real_wfn=use_real_wfn, &
                           cell_to_index=cell_to_index, sab_nl=sab_nl)
      IF (.NOT. ASSOCIATED(cell_to_index) .OR. .NOT. ASSOCIATED(sab_nl)) &
         CPABORT("Incomplete CP2K image map for streamed g-xTB exchange")
      CPASSERT(ALLOCATED(tb%gxtb_kmesh_map))
      kmap => tb%gxtb_kmesh_map
      CPASSERT(kmap%nred == nred)

      nao = tb%calc%bas%nao
      nspin = tb%wfn%nspin
      CALL tb_gxtb_ao_layout(tb, ao_first, ao_size)
      ALLOCATE (density_reduced(nao, nao, nspin, nred), &
                density_block(nao, nao, nspin), fock_block(nao, nao, nspin), &
                overlap_block(nao, nao), transform_work(nao, nao), &
                overlap_part(nao, nao))
      ALLOCATE (exchange_vsh(tb%calc%bas%nsh), source=0.0_dp)
      exchange_energy = 0.0_dp
      CALL tb_collect_kpoint_density(qs_env, tb, kpoints, density_reduced)
      IF (.NOT. ALL(ieee_is_finite(REAL(density_reduced, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(density_reduced)))) &
         CPABORT("CP2K supplied a non-finite irreducible density to streamed save_tblite")
      IF (.NOT. ALL(ieee_is_finite(tb%wfn%qsh))) &
         CPABORT("CP2K supplied non-finite shell charges to streamed save_tblite exchange")

      ALLOCATE (rmatrix, cmatrix)
      CALL dbcsr_create(rmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_nl)

      validate_density = .NOT. kmap%identity .AND. &
                         .NOT. tb%gxtb_density_map_validated .AND. &
                         MAXVAL(ABS(density_reduced)) > 1.0E-14_dp
      IF (validate_density) THEN
         CALL validate_stream_density()
         tb%gxtb_density_map_validated = .TRUE.
      END IF
      IF (.NOT. kmap%identity .AND. .NOT. tb%gxtb_transform_validated) THEN
         run_block_qualification = block_helper_qualification_requested()
         IF (run_block_qualification) THEN
            ! This explicit qualification-only gate intentionally materializes S(k) so that the
            ! block helpers are compared with the unchanged full-array oracle.  Normal streamed
            ! operation uses the reduced-storage validation below and never enters this branch.
            ALLOCATE (overlap_full_qualification(nao, nao, kmap%nfull))
            DO ifull = 1, kmap%nfull
               CALL prepare_stream_block(ifull, .FALSE.)
               overlap_full_qualification(:, :, ifull) = overlap_block
            END DO
            CALL tb_validate_gxtb_kmesh_transform(kpoints, kmap, ao_first, ao_size, &
                                                  overlap_full_qualification)
            DEALLOCATE (overlap_full_qualification)
         ELSE
            CALL validate_stream_overlap()
         END IF
         tb%gxtb_transform_validated = .TRUE.
      END IF

      distributed_mode = stream_mode == tb_gxtb_exchange_stream_distributed
      owns_image_range = .TRUE.
      local_image_first = 1
      local_image_count = kmap%nfull
      IF (distributed_mode) THEN
         nowners = MIN(para_env%num_pe, kmap%nfull)
         owns_image_range = para_env%mepos < nowners
         IF (owns_image_range) THEN
            local_image_count = kmap%nfull/nowners
            IF (para_env%mepos < MOD(kmap%nfull, nowners)) &
               local_image_count = local_image_count + 1
            local_image_first = 1 + para_env%mepos*(kmap%nfull/nowners) + &
                                MIN(para_env%mepos, MOD(kmap%nfull, nowners))
            CPASSERT(local_image_count >= 1)
            CPASSERT(local_image_first >= 1)
            CPASSERT(local_image_first + local_image_count - 1 <= kmap%nfull)
         ELSE
            local_image_first = 0
            local_image_count = 0
         END IF
         ALLOCATE (image_owner_count(kmap%nfull), source=0)
         IF (owns_image_range) &
            image_owner_count(local_image_first:local_image_first + local_image_count - 1) = 1
         CALL para_env%sum(image_owner_count)
         IF (ANY(image_owner_count /= 1)) &
            CPABORT("Distributed g-xTB image ranges do not form an exact nonoverlapping cover")
         DEALLOCATE (image_owner_count)
      END IF

      batch_size = 0
      IF (distributed_mode) THEN
         IF (owns_image_range) THEN
            batch_size = local_image_count
            CALL cp2k_exchange_stream_begin(stream, tb%calc, tb%mol, tb%ecache, tb%wfn, &
                                            kmap%nmesh, kmap%xkp, kmap%weight, error, &
                                            mode=cp2k_exchange_stream_batched, &
                                            batch_size=batch_size, &
                                            image_first=local_image_first, &
                                            image_count=local_image_count)
         ELSE
            ! Ranks beyond the number of BvK images still need the current
            ! geometry/charge-dependent exchange cache for the replicated
            ! analytical derivative path used after the SCF calculation.
            CALL cp2k_prepare_exchange(tb%calc, tb%mol, tb%ecache, tb%wfn, error)
         END IF
      ELSE IF (stream_mode == cp2k_exchange_stream_batched) THEN
         batch_size = tb_gxtb_exchange_stream_batch_size()
         CALL cp2k_exchange_stream_begin(stream, tb%calc, tb%mol, tb%ecache, tb%wfn, &
                                         kmap%nmesh, kmap%xkp, kmap%weight, error, &
                                         mode=stream_mode, batch_size=batch_size)
      ELSE
         CALL cp2k_exchange_stream_begin(stream, tb%calc, tb%mol, tb%ecache, tb%wfn, &
                                         kmap%nmesh, kmap%xkp, kmap%weight, error, &
                                         mode=stream_mode)
      END IF
      IF (ALLOCATED(error)) CPABORT("save_tblite exchange stream initialization failed")
      has_full_kspace_inputs = .FALSE.
      IF (owns_image_range) THEN
         has_full_kspace_inputs = cp2k_exchange_stream_has_full_mesh_storage(stream)
         IF (distributed_mode) THEN
            expected_full_mesh_storage = local_image_count >= kmap%nfull
         ELSE
            expected_full_mesh_storage = stream_mode /= cp2k_exchange_stream_batched .OR. &
                                         MIN(batch_size, kmap%nfull) >= kmap%nfull
         END IF
         IF (has_full_kspace_inputs .NEQV. expected_full_mesh_storage) &
            CPABORT("Unexpected save_tblite exchange-stream storage class")
      END IF
      ALLOCATE (tb%exchange_fock_kp(nao, nao, nspin, nred), &
                source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))

      IF (distributed_mode) THEN
         ! Every rank enters prepare_stream_block in the same order because it
         ! contains collectives over CP2K's distributed AO matrices.  Only the
         ! owner ranks call their independent provider transaction.
         IF (owns_image_range) THEN
            CALL cp2k_exchange_stream_batch_bounds(stream, first_image, last_image, error)
            IF (ALLOCATED(error)) CPABORT("save_tblite distributed image bounds failed")
            IF (first_image /= local_image_first .OR. &
                last_image /= local_image_first + local_image_count - 1) &
               CPABORT("save_tblite returned a different distributed image range")
         END IF
         DO ifull = 1, kmap%nfull
            CALL prepare_stream_block(ifull, .TRUE.)
            IF (owns_image_range) THEN
               CALL cp2k_exchange_stream_push(stream, ifull, density_block, overlap_block, error)
               IF (ALLOCATED(error)) &
                  CPABORT("save_tblite distributed image range rejected a Bloch block")
            END IF
         END DO
         IF (owns_image_range) THEN
            CALL cp2k_exchange_stream_batch_apply(stream, tb%calc, tb%mol, tb%ecache, error)
            IF (ALLOCATED(error)) CPABORT("save_tblite distributed image evaluation failed")
         END IF
         DO ifull = 1, kmap%nfull
            CALL prepare_stream_block(ifull, .FALSE.)
            IF (owns_image_range) THEN
               fock_block = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
               CALL cp2k_exchange_stream_batch_pull(stream, ifull, overlap_block, &
                                                    fock_block, error)
               IF (ALLOCATED(error)) &
                  CPABORT("save_tblite distributed image Fock retrieval failed")
               CALL fold_stream_fock_block(ifull, .FALSE.)
            END IF
         END DO
         IF (owns_image_range) THEN
            CALL cp2k_exchange_stream_batch_advance(stream, batches_done, error)
            IF (ALLOCATED(error)) CPABORT("save_tblite distributed image advance failed")
            IF (.NOT. batches_done) CPABORT("save_tblite distributed image range was not complete")
            CALL cp2k_exchange_stream_batch_result(stream, tb%calc, tb%mol, tb%ecache, &
                                                   exchange_vsh, exchange_energy, error)
            IF (ALLOCATED(error)) CPABORT("save_tblite distributed image result failed")
         END IF
         CALL para_env%sum(exchange_energy)
         CALL para_env%sum(exchange_vsh)
         CALL para_env%sum(tb%exchange_fock_kp)
      ELSE IF (stream_mode == cp2k_exchange_stream_batched) THEN
         batches_done = .FALSE.
         DO WHILE (.NOT. batches_done)
            CALL cp2k_exchange_stream_batch_bounds(stream, first_image, last_image, error)
            IF (ALLOCATED(error)) CPABORT("save_tblite exchange image-batch bounds failed")
            IF (first_image < 1 .OR. last_image < first_image .OR. &
                last_image > kmap%nfull .OR. last_image - first_image + 1 > batch_size) &
               CPABORT("save_tblite returned invalid exchange image-batch bounds")

            DO ifull = 1, kmap%nfull
               CALL prepare_stream_block(ifull, .TRUE.)
               CALL cp2k_exchange_stream_push(stream, ifull, density_block, overlap_block, error)
               IF (ALLOCATED(error)) &
                  CPABORT("save_tblite exchange image batch rejected a Bloch block")
            END DO
            CALL cp2k_exchange_stream_batch_apply(stream, tb%calc, tb%mol, tb%ecache, error)
            IF (ALLOCATED(error)) CPABORT("save_tblite exchange image-batch evaluation failed")

            DO ifull = 1, kmap%nfull
               CALL prepare_stream_block(ifull, .FALSE.)
               fock_block = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
               CALL cp2k_exchange_stream_batch_pull(stream, ifull, overlap_block, &
                                                    fock_block, error)
               IF (ALLOCATED(error)) &
                  CPABORT("save_tblite exchange image-batch Fock retrieval failed")
               CALL fold_stream_fock_block(ifull, .FALSE.)
            END DO
            CALL cp2k_exchange_stream_batch_advance(stream, batches_done, error)
            IF (ALLOCATED(error)) CPABORT("save_tblite exchange image-batch advance failed")
         END DO
         CALL cp2k_exchange_stream_batch_result(stream, tb%calc, tb%mol, tb%ecache, &
                                                exchange_vsh, exchange_energy, error)
         IF (ALLOCATED(error)) CPABORT("save_tblite exchange image-batch result failed")
      ELSE
         DO ifull = 1, kmap%nfull
            CALL prepare_stream_block(ifull, .TRUE.)
            CALL cp2k_exchange_stream_push(stream, ifull, density_block, overlap_block, error)
            IF (ALLOCATED(error)) CPABORT("save_tblite exchange stream rejected a Bloch block")
         END DO
         CALL cp2k_exchange_stream_apply(stream, tb%calc, tb%mol, tb%ecache, exchange_vsh, &
                                         exchange_energy, error)
         IF (ALLOCATED(error)) CPABORT("save_tblite exchange stream evaluation failed")
         DO ifull = 1, kmap%nfull
            CALL prepare_stream_block(ifull, .FALSE.)
            CALL cp2k_exchange_stream_get(stream, ifull, fock_block, error, overlap_block)
            IF (ALLOCATED(error)) CPABORT("save_tblite exchange stream Fock retrieval failed")
            CALL fold_stream_fock_block(ifull, .TRUE.)
         END DO
      END IF

      IF (.NOT. ALL(ieee_is_finite(exchange_vsh)) .OR. &
          .NOT. ieee_is_finite(exchange_energy)) &
         CPABORT("save_tblite returned a non-finite streamed exchange response")

      DO source = 1, nred
         DO ispin = 1, nspin
            hermitian_error = MAXVAL(ABS(tb%exchange_fock_kp(:, :, ispin, source) - &
                                         CONJG(TRANSPOSE( &
                                               tb%exchange_fock_kp(:, :, ispin, source)))))
            IF (hermitian_error > 1.0E-10_dp) &
               CPABORT("Folded streamed g-xTB exchange matrix is not Hermitian")
            tb%exchange_fock_kp(:, :, ispin, source) = &
               0.5_dp*(tb%exchange_fock_kp(:, :, ispin, source) + &
                       CONJG(TRANSPOSE(tb%exchange_fock_kp(:, :, ispin, source))))
         END DO
         IF (use_real_wfn .AND. &
             MAXVAL(ABS(AIMAG(tb%exchange_fock_kp(:, :, :, source)))) > 1.0E-10_dp) &
            CPABORT("Complex streamed exchange cannot be represented by real wavefunctions")
      END DO

      IF (owns_image_range) THEN
         CALL cp2k_exchange_stream_end(stream, error)
         IF (ALLOCATED(error)) CPABORT("save_tblite exchange stream finalization failed")
      END IF
      IF (run_qualification_oracle) CALL qualify_stream_oracle()
      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)

      tb%e_exchange = 0.0_dp
      tb%e_exchange(1) = exchange_energy
      tb%pot%vsh(:, 1) = tb%pot%vsh(:, 1) + exchange_vsh
      tb%exchange_fock_kp_valid = .TRUE.
      DEALLOCATE (ao_first, ao_size, density_block, density_reduced, exchange_vsh, &
                  fock_block, overlap_block, overlap_part, transform_work)

   CONTAINS

! **************************************************************************************************
!> \brief Return whether the explicit block-helper/full-array qualification was requested.
!> \return ...
! **************************************************************************************************
      LOGICAL FUNCTION block_helper_qualification_requested() RESULT(requested)

      CHARACTER(LEN=32)                                  :: value
      INTEGER                                            :: status

         requested = .FALSE.
         value = ""
         status = 1
         CALL GET_ENVIRONMENT_VARIABLE("CP2K_GXTB_QUALIFICATION_BLOCK_HELPERS", &
                                       value, STATUS=status)
         IF (status /= 0) RETURN
         SELECT CASE (TRIM(ADJUSTL(value)))
         CASE ("", "0", "OFF", "off", "FALSE", "false")
            requested = .FALSE.
         CASE DEFAULT
            ! The full validator interprets valid nonzero modes and diagnoses unknown values.
            requested = .TRUE.
         END SELECT

      END FUNCTION block_helper_qualification_requested

! **************************************************************************************************
!> \brief Validate streamed density expansion without materializing a complete density mesh.
! **************************************************************************************************
      SUBROUTINE validate_stream_density()

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: expanded, transformed
      INTEGER                                            :: ifull_local, iop_local, ispin_local, &
                                                            source_local
      REAL(KIND=dp)                                      :: residual, scale, symmetry_tolerance

         ALLOCATE (expanded(nao, nao), transformed(nao, nao))
         symmetry_tolerance = MAX(1.0E-10_dp, 100.0_dp*kpoints%eps_geo)
         DO ifull_local = 1, kmap%nfull
            source_local = kmap%source(ifull_local)
            DO ispin_local = 1, nspin
               CALL tb_expand_gxtb_kpoint_block(kpoints, kmap, ao_first, ao_size, ifull_local, &
                                                density_reduced(:, :, ispin_local, source_local), &
                                                expanded, transform_work)
               DO iop_local = 1, kmap%op_count(ifull_local)
                  CALL tb_transform_gxtb_kpoint_matrix( &
                     kpoints, source_local, kmap%symops(iop_local, ifull_local), &
                     kmap%rotation_slot(iop_local, ifull_local), ao_first, ao_size, &
                     density_reduced(:, :, ispin_local, source_local), transformed, .FALSE.)
                  scale = MAX(1.0_dp, MAXVAL(ABS(expanded)))
                  residual = MAXVAL(ABS(transformed - expanded))/scale
                  IF (residual > symmetry_tolerance) &
                     CPABORT("Streamed g-xTB density is not invariant under its CP2K little group")
               END DO
            END DO
         END DO
         DEALLOCATE (expanded, transformed)

      END SUBROUTINE validate_stream_density

! **************************************************************************************************
!> \brief Validate direct-overlap covariance and its real adjoint with reduced storage.
! **************************************************************************************************
      SUBROUTINE validate_stream_overlap()

      CHARACTER(LEN=256)                                 :: error_message
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: adjoint, transformed
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: overlap_reduced
      INTEGER                                            :: ifull_local, iop_local, source_local
      REAL(KIND=dp)                                      :: adjoint_residual, covariance_residual, &
                                                            lhs, rhs, scale

         ALLOCATE (adjoint(nao, nao), transformed(nao, nao), &
                   overlap_reduced(nao, nao, kmap%nred))
         DO source_local = 1, kmap%nred
            ifull_local = kmap%representative(source_local)
            CPASSERT(ifull_local >= 1 .AND. ifull_local <= kmap%nfull)
            CALL prepare_stream_block(ifull_local, .FALSE.)
            overlap_reduced(:, :, source_local) = overlap_block
         END DO

         DO ifull_local = 1, kmap%nfull
            source_local = kmap%source(ifull_local)
            CALL prepare_stream_block(ifull_local, .FALSE.)
            DO iop_local = 1, kmap%op_count(ifull_local)
               CALL tb_transform_gxtb_kpoint_matrix( &
                  kpoints, source_local, kmap%symops(iop_local, ifull_local), &
                  kmap%rotation_slot(iop_local, ifull_local), ao_first, ao_size, &
                  overlap_reduced(:, :, source_local), transformed, .FALSE.)
               scale = MAX(1.0_dp, MAXVAL(ABS(overlap_block)))
               covariance_residual = MAXVAL(ABS(transformed - overlap_block))/scale
               IF (covariance_residual > MAX(1.0E-6_dp, 100.0_dp*kpoints%eps_geo)) THEN
                  WRITE (error_message, '(A,I0,A,I0,A,I0,A,ES12.4)') &
                     "Streamed CP2K overlap symmetry failed for full k-point ", ifull_local, &
                     ", irreducible source ", source_local, ", symmetry operation ", &
                     kmap%symops(iop_local, ifull_local), "; residual=", covariance_residual
                  CALL cp_abort(__LOCATION__, TRIM(error_message))
               END IF

               CALL tb_transform_gxtb_kpoint_matrix( &
                  kpoints, source_local, kmap%symops(iop_local, ifull_local), &
                  kmap%rotation_slot(iop_local, ifull_local), ao_first, ao_size, &
                  overlap_block, adjoint, .TRUE.)
               lhs = REAL(SUM(CONJG(overlap_block)*transformed), KIND=dp)
               rhs = REAL(SUM(CONJG(adjoint)*overlap_reduced(:, :, source_local)), KIND=dp)
               adjoint_residual = ABS(lhs - rhs)/MAX(1.0_dp, ABS(lhs), ABS(rhs))
               IF (adjoint_residual > 1.0E-10_dp) THEN
                  WRITE (error_message, '(A,I0,A,I0,A,I0,A,ES12.4)') &
                     "Streamed CP2K symmetry adjoint failed for full k-point ", ifull_local, &
                     ", irreducible source ", source_local, ", symmetry operation ", &
                     kmap%symops(iop_local, ifull_local), "; residual=", adjoint_residual
                  CALL cp_abort(__LOCATION__, TRIM(error_message))
               END IF
            END DO
         END DO
         DEALLOCATE (adjoint, overlap_reduced, transformed)

      END SUBROUTINE validate_stream_overlap

! **************************************************************************************************
!> \brief Compare the streamed folded response with the unchanged full-mesh provider oracle.
! **************************************************************************************************
      SUBROUTINE qualify_stream_oracle()

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: overlap_full
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: density_full, fock_full, fock_reduced
      INTEGER                                            :: ifull_local, ispin_local, source_local
      REAL(KIND=dp) :: duality_full, duality_reduced, duality_relative, energy_delta, &
         energy_oracle, fock_delta, hermitian_residual, oracle_duality_residual, &
         oracle_fold_covariance_residual, oracle_full_covariance_residual, &
         stream_covariance_residual, stream_duality_residual, vsh_delta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vsh_oracle

         ALLOCATE (density_full(nao, nao, nspin, kmap%nfull), &
                   fock_full(nao, nao, nspin, kmap%nfull), &
                   fock_reduced(nao, nao, nspin, nred), &
                   overlap_full(nao, nao, kmap%nfull), &
                   source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))
         ALLOCATE (vsh_oracle(SIZE(exchange_vsh)), source=0.0_dp)
         DO source_local = 1, nred
            ifull_local = kmap%representative(source_local)
            CPASSERT(ifull_local >= 1 .AND. ifull_local <= kmap%nfull)
            density_full(:, :, :, ifull_local) = density_reduced(:, :, :, source_local)
         END DO
         CALL tb_expand_gxtb_kpoint_density(kpoints, kmap, ao_first, ao_size, density_full, .TRUE.)
         DO ifull_local = 1, kmap%nfull
            CALL prepare_stream_block(ifull_local, .FALSE.)
            overlap_full(:, :, ifull_local) = overlap_block
         END DO
         CALL cp2k_exchange_kmesh(tb%calc, tb%mol, tb%ecache, tb%wfn, kmap%nmesh, &
                                  kmap%xkp, kmap%weight, density_full, overlap_full, fock_full, &
                                  vsh_oracle, energy_oracle, error, &
                                  transform_backend=cp2k_exchange_transform_dense)
         IF (ALLOCATED(error)) &
            CPABORT("QUALIFICATION_ONLY streamed full-mesh exchange oracle failed")
         IF (.NOT. ieee_is_finite(energy_oracle)) &
            CPABORT("QUALIFICATION_ONLY streamed oracle returned a non-finite energy")
         IF (.NOT. ALL(ieee_is_finite(vsh_oracle))) &
            CPABORT("QUALIFICATION_ONLY streamed oracle returned a non-finite shell response")
         IF (.NOT. ALL(ieee_is_finite(REAL(fock_full, KIND=dp))) .OR. &
             .NOT. ALL(ieee_is_finite(AIMAG(fock_full)))) &
            CPABORT("QUALIFICATION_ONLY streamed oracle returned a non-finite Fock response")
         hermitian_residual = 0.0_dp
         DO ifull_local = 1, kmap%nfull
            DO ispin_local = 1, nspin
               hermitian_residual = MAX(hermitian_residual, &
                                        MAXVAL(ABS(fock_full(:, :, ispin_local, ifull_local) - &
                                                   CONJG(TRANSPOSE( &
                                                         fock_full(:, :, ispin_local, ifull_local))))))
            END DO
         END DO
         IF (hermitian_residual > 1.0E-10_dp) &
            CPABORT("QUALIFICATION_ONLY streamed oracle returned a non-Hermitian full-k Fock response")
         fock_reduced = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         oracle_full_covariance_residual = 0.0_dp
         oracle_fold_covariance_residual = 0.0_dp
         stream_covariance_residual = 0.0_dp
         oracle_duality_residual = 0.0_dp
         stream_duality_residual = 0.0_dp
         DO ispin_local = 1, nspin
            CALL tb_validate_gxtb_full_response(kpoints, kmap, ao_first, ao_size, &
                                                fock_full(:, :, ispin_local, :), &
                                                "QUALIFICATION_ONLY streamed-oracle Fock", &
                                                duality_relative)
            oracle_full_covariance_residual = &
               MAX(oracle_full_covariance_residual, duality_relative)
            CALL tb_fold_gxtb_kpoint_matrix(kpoints, kmap, ao_first, ao_size, &
                                            fock_full(:, :, ispin_local, :), &
                                            fock_reduced(:, :, ispin_local, :))
            CALL tb_validate_gxtb_folded_response(kpoints, kmap, ao_first, ao_size, &
                                                  fock_full(:, :, ispin_local, :), &
                                                  fock_reduced(:, :, ispin_local, :), &
                                                  "QUALIFICATION_ONLY streamed-oracle Fock", &
                                                  duality_relative)
            oracle_fold_covariance_residual = &
               MAX(oracle_fold_covariance_residual, duality_relative)
            CALL tb_qualification_gxtb_duality( &
               kpoints, kmap, density_full(:, :, ispin_local, :), &
               density_reduced(:, :, ispin_local, :), fock_full(:, :, ispin_local, :), &
               fock_reduced(:, :, ispin_local, :), duality_full, duality_reduced, &
               duality_relative)
            oracle_duality_residual = MAX(oracle_duality_residual, duality_relative)
            IF (duality_relative > 1.0E-10_dp) &
               CPABORT("QUALIFICATION_ONLY streamed oracle-fold duality exceeded 1.0E-10")
            CALL tb_validate_gxtb_folded_response( &
               kpoints, kmap, ao_first, ao_size, fock_full(:, :, ispin_local, :), &
               tb%exchange_fock_kp(:, :, ispin_local, :), &
               "QUALIFICATION_ONLY direct streamed Fock", duality_relative)
            stream_covariance_residual = MAX(stream_covariance_residual, duality_relative)
            CALL tb_qualification_gxtb_duality( &
               kpoints, kmap, density_full(:, :, ispin_local, :), &
               density_reduced(:, :, ispin_local, :), fock_full(:, :, ispin_local, :), &
               tb%exchange_fock_kp(:, :, ispin_local, :), duality_full, duality_reduced, &
               duality_relative)
            stream_duality_residual = MAX(stream_duality_residual, duality_relative)
            IF (duality_relative > 1.0E-10_dp) &
               CPABORT("QUALIFICATION_ONLY direct streamed duality exceeded 1.0E-10")
         END DO
         energy_delta = ABS(exchange_energy - energy_oracle)
         vsh_delta = MAXVAL(ABS(exchange_vsh - vsh_oracle))
         fock_delta = MAXVAL(ABS(tb%exchange_fock_kp - fock_reduced))
         IF (para_env%is_source()) THEN
            WRITE (*, '(A,I0,9(A,ES14.6))') &
               " GXTB-QUALIFICATION_ONLY STREAM-ORACLE iter=", iter_count, &
               " dE=", energy_delta, " dVsh=", vsh_delta, &
               " dFfold=", fock_delta, " hermFull=", hermitian_residual, &
               " covOracleFull=", oracle_full_covariance_residual, &
               " covOracleFold=", oracle_fold_covariance_residual, &
               " dualOracleFold=", oracle_duality_residual, &
               " covStream=", stream_covariance_residual, &
               " dualStream=", stream_duality_residual
         END IF
         IF (energy_delta > 1.0E-10_dp .OR. vsh_delta > 1.0E-10_dp .OR. &
             fock_delta > 1.0E-10_dp .OR. hermitian_residual > 1.0E-10_dp .OR. &
             oracle_full_covariance_residual > 1.0E-10_dp .OR. &
             oracle_fold_covariance_residual > 1.0E-10_dp .OR. &
             oracle_duality_residual > 1.0E-10_dp .OR. &
             stream_covariance_residual > 1.0E-10_dp .OR. &
             stream_duality_residual > 1.0E-10_dp) &
            CPABORT("QUALIFICATION_ONLY streamed full-mesh oracle exceeded its 1.0E-10 gate")
         DEALLOCATE (density_full, fock_full, fock_reduced, overlap_full, vsh_oracle)

      END SUBROUTINE qualify_stream_oracle

! **************************************************************************************************
!> \brief ...
!> \param ifull_local ...
!> \param include_density ...
! **************************************************************************************************
      SUBROUTINE prepare_stream_block(ifull_local, include_density)

      INTEGER, INTENT(IN)                                :: ifull_local
      LOGICAL, INTENT(IN)                                :: include_density

      INTEGER                                            :: ispin_local, source_local
      REAL(KIND=dp)                                      :: local_hermitian_error

         source_local = kmap%source(ifull_local)
         IF (include_density) THEN
            DO ispin_local = 1, nspin
               CALL tb_expand_gxtb_kpoint_block(kpoints, kmap, ao_first, ao_size, ifull_local, &
                                                density_reduced(:, :, ispin_local, source_local), &
                                                density_block(:, :, ispin_local), transform_work)
               local_hermitian_error = &
                  MAXVAL(ABS(density_block(:, :, ispin_local) - &
                             CONJG(TRANSPOSE(density_block(:, :, ispin_local)))))
               IF (local_hermitian_error > 1.0E-10_dp) &
                  CPABORT("Streamed g-xTB density block is not Hermitian")
            END DO
         END IF

         CALL dbcsr_set(rmatrix, 0.0_dp)
         CALL dbcsr_set(cmatrix, 0.0_dp)
         CALL rskp_transform(rmatrix, cmatrix, matrix_s, 1, kmap%xkp(1:3, ifull_local), &
                             cell_to_index, sab_nl)
         CALL tb_dbcsr_component_to_dense(tb, rmatrix, overlap_part, para_env, .FALSE.)
         overlap_block(:, :) = CMPLX(overlap_part, 0.0_dp, KIND=dp)
         CALL tb_dbcsr_component_to_dense(tb, cmatrix, overlap_part, para_env, .TRUE.)
         overlap_block(:, :) = CMPLX(REAL(overlap_block, KIND=dp), overlap_part, KIND=dp)
         IF (.NOT. ALL(ieee_is_finite(REAL(overlap_block, KIND=dp))) .OR. &
             .NOT. ALL(ieee_is_finite(AIMAG(overlap_block)))) &
            CPABORT("CP2K supplied a non-finite overlap to streamed save_tblite")
         local_hermitian_error = &
            MAXVAL(ABS(overlap_block - CONJG(TRANSPOSE(overlap_block))))
         IF (local_hermitian_error > 1.0E-10_dp) &
            CPABORT("CP2K supplied a non-Hermitian overlap to streamed save_tblite")

      END SUBROUTINE prepare_stream_block

! **************************************************************************************************
!> \brief ...
!> \param ifull_local ...
!> \param require_hermitian ...
! **************************************************************************************************
      SUBROUTINE fold_stream_fock_block(ifull_local, require_hermitian)

      INTEGER, INTENT(IN)                                :: ifull_local
      LOGICAL, INTENT(IN)                                :: require_hermitian

      INTEGER                                            :: ispin_local, source_local
      REAL(KIND=dp)                                      :: local_hermitian_error

         IF (.NOT. ALL(ieee_is_finite(REAL(fock_block, KIND=dp))) .OR. &
             .NOT. ALL(ieee_is_finite(AIMAG(fock_block)))) &
            CPABORT("save_tblite returned a non-finite streamed Fock block")
         source_local = kmap%source(ifull_local)
         DO ispin_local = 1, nspin
            local_hermitian_error = &
               MAXVAL(ABS(fock_block(:, :, ispin_local) - &
                          CONJG(TRANSPOSE(fock_block(:, :, ispin_local)))))
            IF (require_hermitian .AND. local_hermitian_error > 1.0E-10_dp) &
               CPABORT("save_tblite returned a non-Hermitian streamed Fock block")
            CALL tb_accumulate_gxtb_kpoint_response(kpoints, kmap, ao_first, ao_size, &
                                                    ifull_local, fock_block(:, :, ispin_local), &
                                                    tb%exchange_fock_kp(:, :, ispin_local, &
                                                                        source_local), &
                                                    transform_work)
         END DO

      END SUBROUTINE fold_stream_fock_block

   END SUBROUTINE tb_build_gxtb_kpoint_exchange_streamed

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> \param matrix_s ...
!> \param para_env ...
!> Build g-xTB exchange from the complete regular mesh represented by CP2K's k-points.
!> \param iter_count ...
! **************************************************************************************************
   SUBROUTINE tb_build_gxtb_kpoint_exchange(qs_env, tb, matrix_s, para_env, iter_count)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: iter_count

      CHARACTER(LEN=32)                                  :: qualification_value
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: overlap, overlap_oracle
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :) :: density, density_oracle, &
                                              density_reduced, exchange_fock_full, exchange_fock_oracle, exchange_fock_unpermuted, &
                                                              folded_fock_oracle
      INTEGER                                            :: ik, ispin, nao, nred, nspin, &
                                                            qualification_iteration, &
                                                            qualification_status, source_index, &
                                                            stream_mode, transform_backend
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_first, ao_size
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: run_qualification_duality, &
                                                            run_qualification_oracle, &
                                                            use_real_wfn, validate_density
      REAL(KIND=dp) :: acp_duality_full, acp_duality_reduced, acp_duality_relative, &
                       acp_duality_total_full, acp_duality_total_reduced, duality_tolerance, &
                       exchange_duality_full, exchange_duality_reduced, exchange_duality_relative, &
                       exchange_duality_total_full, exchange_duality_total_reduced, exchange_energy, &
                       exchange_energy_oracle, hermitian_error, oracle_energy_delta, oracle_fold_abs, &
                       oracle_fold_rel, oracle_full_abs, oracle_full_rel, oracle_scale, oracle_sym_fold, &
                       oracle_sym_full, oracle_tolerance, oracle_vsh_delta, primary_sym_fold, primary_sym_full, &
                       symmetry_residual
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: exchange_vsh, exchange_vsh_oracle, &
                                                            weight_oracle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: overlap_part, xkp_oracle
      TYPE(dbcsr_type), POINTER                          :: cmatrix, rmatrix
      TYPE(error_type), ALLOCATABLE                      :: error
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(tb_gxtb_kmesh_map_type), POINTER              :: kmap

      NULLIFY (cell_to_index, cmatrix, kmap, kpoints, rmatrix, sab_nl)
      CALL get_qs_env(qs_env, kpoints=kpoints)
      CPASSERT(ASSOCIATED(kpoints))
      CALL get_kpoint_info(kpoints, nkp=nred, &
                           use_real_wfn=use_real_wfn, cell_to_index=cell_to_index, sab_nl=sab_nl)
      IF (.NOT. ASSOCIATED(cell_to_index) .OR. .NOT. ASSOCIATED(sab_nl)) &
         CPABORT("Incomplete CP2K image map for g-xTB exchange")
      CALL tb_prepare_gxtb_kmesh_map(tb, kpoints, matrix_s, sab_nl, cell_to_index)
      kmap => tb%gxtb_kmesh_map
      qualification_iteration = -1
      qualification_status = 1
      qualification_value = ""
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_GXTB_QUALIFICATION_FULLMESH_ORACLE_ITERATION", &
                                    qualification_value, STATUS=qualification_status)
      IF (qualification_status == 0) &
         READ (qualification_value, *, IOSTAT=qualification_status) qualification_iteration
      run_qualification_oracle = qualification_status == 0 .AND. &
                                 qualification_iteration == iter_count
      stream_mode = tb_gxtb_exchange_stream_mode()
      transform_backend = tb_gxtb_exchange_transform_backend()
      IF (run_qualification_oracle .AND. kmap%identity .AND. &
          stream_mode /= cp2k_exchange_stream_batched .AND. &
          stream_mode /= tb_gxtb_exchange_stream_distributed .AND. &
          stream_mode /= tb_gxtb_exchange_stream_kgroup_owner .AND. &
          .NOT. (stream_mode == 0 .AND. &
                 transform_backend == cp2k_exchange_transform_separable_dft)) &
         CPABORT("QUALIFICATION_ONLY identity oracle requires batched, distributed, or separable mode")
      tb%exchange_fock_kp_valid = .FALSE.
      tb%mixed_fock_kp_valid = .FALSE.
      IF (ALLOCATED(tb%exchange_fock_kp)) DEALLOCATE (tb%exchange_fock_kp)
      IF (ALLOCATED(tb%mixed_fock_kp)) DEALLOCATE (tb%mixed_fock_kp)

      IF (stream_mode == tb_gxtb_exchange_stream_kgroup_owner) THEN
         CALL tb_build_gxtb_kpoint_exchange_kgroup_owner( &
            qs_env, tb, matrix_s, run_qualification_oracle, iter_count)
         RETURN
      ELSE IF (stream_mode /= 0) THEN
         CALL tb_build_gxtb_kpoint_exchange_streamed(qs_env, tb, matrix_s, para_env, stream_mode, &
                                                     run_qualification_oracle, iter_count)
         RETURN
      END IF

      nao = tb%calc%bas%nao
      nspin = tb%wfn%nspin
      CALL tb_gxtb_ao_layout(tb, ao_first, ao_size)
      ALLOCATE (density(nao, nao, nspin, kmap%nfull), overlap(nao, nao, kmap%nfull), &
                exchange_fock_full(nao, nao, nspin, kmap%nfull), &
                overlap_part(nao, nao), &
                exchange_vsh(tb%calc%bas%nsh))
      IF (run_qualification_oracle) ALLOCATE (density_reduced(nao, nao, nspin, nred))
      ! Before the first diagonalization of a freshly initialized k-point SCF (including an
      ! RS-only restart), CP2K has no P(k) yet.  The zero pmat is a deliberate initial Hamiltonian;
      ! every subsequent SCC step uses the exact P(k) produced by the preceding diagonalization,
      ! without an alias-prone inverse image transform.
      IF (kmap%identity) THEN
         CALL tb_collect_kpoint_density(qs_env, tb, kpoints, density)
      ELSE
         CALL tb_collect_kpoint_density(qs_env, tb, kpoints, density, kmap%representative)
         IF (run_qualification_oracle) THEN
            DO ik = 1, nred
               density_reduced(:, :, :, ik) = density(:, :, :, kmap%representative(ik))
            END DO
         END IF
         validate_density = .NOT. tb%gxtb_density_map_validated .AND. &
                            MAXVAL(ABS(density)) > 1.0E-14_dp
         CALL tb_expand_gxtb_kpoint_density(kpoints, kmap, ao_first, ao_size, density, &
                                            validate_density)
         IF (validate_density) tb%gxtb_density_map_validated = .TRUE.
      END IF
      IF (.NOT. ALL(ieee_is_finite(REAL(density, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(density)))) &
         CPABORT("CP2K supplied a non-finite k-point density matrix to save_tblite")

      IF (.NOT. ALL(ieee_is_finite(tb%wfn%qsh))) &
         CPABORT("CP2K supplied non-finite shell charges to save_tblite exchange")

      ALLOCATE (rmatrix, cmatrix)
      CALL dbcsr_create(rmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_nl)

      DO ik = 1, kmap%nfull
         CALL dbcsr_set(rmatrix, 0.0_dp)
         CALL dbcsr_set(cmatrix, 0.0_dp)
         CALL rskp_transform(rmatrix, cmatrix, matrix_s, 1, kmap%xkp(1:3, ik), &
                             cell_to_index, sab_nl)
         CALL tb_dbcsr_component_to_dense(tb, rmatrix, overlap_part, para_env, .FALSE.)
         IF (.NOT. ALL(ieee_is_finite(overlap_part))) &
            CPABORT("CP2K supplied a non-finite k-point overlap matrix to save_tblite")
         overlap(:, :, ik) = CMPLX(overlap_part, 0.0_dp, KIND=dp)
         CALL tb_dbcsr_component_to_dense(tb, cmatrix, overlap_part, para_env, .TRUE.)
         IF (.NOT. ALL(ieee_is_finite(overlap_part))) &
            CPABORT("CP2K supplied a non-finite k-point overlap matrix to save_tblite")
         overlap(:, :, ik) = CMPLX(REAL(overlap(:, :, ik), KIND=dp), overlap_part, KIND=dp)

         hermitian_error = MAXVAL(ABS(overlap(:, :, ik) - &
                                      CONJG(TRANSPOSE(overlap(:, :, ik)))))
         IF (hermitian_error > 1.0E-10_dp) &
            CPABORT("CP2K supplied a non-Hermitian k-point overlap matrix to save_tblite")
         DO ispin = 1, nspin
            hermitian_error = MAXVAL(ABS(density(:, :, ispin, ik) - &
                                         CONJG(TRANSPOSE(density(:, :, ispin, ik)))))
            IF (hermitian_error > 1.0E-10_dp) &
               CPABORT("CP2K supplied a non-Hermitian k-point density matrix to save_tblite")
         END DO

      END DO

      IF (.NOT. kmap%identity .AND. .NOT. tb%gxtb_transform_validated) THEN
         CALL tb_validate_gxtb_kmesh_transform(kpoints, kmap, ao_first, ao_size, overlap)
         tb%gxtb_transform_validated = .TRUE.
      END IF
      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)
      DEALLOCATE (overlap_part)
      CALL cp2k_exchange_kmesh(tb%calc, tb%mol, tb%ecache, tb%wfn, kmap%nmesh, &
                               kmap%xkp, kmap%weight, density, overlap, exchange_fock_full, exchange_vsh, &
                               exchange_energy, error, transform_backend=transform_backend)
      IF (ALLOCATED(error)) CPABORT("save_tblite whole-mesh exchange evaluation failed")
      IF (.NOT. ALL(ieee_is_finite(REAL(exchange_fock_full, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(exchange_fock_full)))) &
         CPABORT("save_tblite returned a non-finite whole-mesh exchange matrix")
      IF (.NOT. ALL(ieee_is_finite(exchange_vsh)) .OR. &
          .NOT. ieee_is_finite(exchange_energy)) &
         CPABORT("save_tblite returned a non-finite whole-mesh exchange response")
      DO ik = 1, kmap%nfull
         DO ispin = 1, nspin
            hermitian_error = MAXVAL(ABS(exchange_fock_full(:, :, ispin, ik) - &
                                         CONJG(TRANSPOSE(exchange_fock_full(:, :, ispin, ik)))))
            IF (hermitian_error > 1.0E-10_dp) &
               CPABORT("save_tblite returned a non-Hermitian k-point exchange matrix")
         END DO
      END DO
      ALLOCATE (tb%exchange_fock_kp(nao, nao, nspin, nred), &
                source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      primary_sym_full = 0.0_dp
      primary_sym_fold = 0.0_dp
      DO ispin = 1, nspin
         IF (.NOT. kmap%identity) THEN
            CALL tb_validate_gxtb_full_response(kpoints, kmap, ao_first, ao_size, &
                                                exchange_fock_full(:, :, ispin, :), "Fock", &
                                                symmetry_residual)
            primary_sym_full = MAX(primary_sym_full, symmetry_residual)
         END IF
         CALL tb_fold_gxtb_kpoint_matrix(kpoints, kmap, ao_first, ao_size, &
                                         exchange_fock_full(:, :, ispin, :), &
                                         tb%exchange_fock_kp(:, :, ispin, :))
         IF (.NOT. kmap%identity) THEN
            CALL tb_validate_gxtb_folded_response(kpoints, kmap, ao_first, ao_size, &
                                                  exchange_fock_full(:, :, ispin, :), &
                                                  tb%exchange_fock_kp(:, :, ispin, :), "Fock", &
                                                  symmetry_residual)
            primary_sym_fold = MAX(primary_sym_fold, symmetry_residual)
         END IF
      END DO
      DO ik = 1, nred
         DO ispin = 1, nspin
            hermitian_error = MAXVAL(ABS(tb%exchange_fock_kp(:, :, ispin, ik) - &
                                         CONJG(TRANSPOSE(tb%exchange_fock_kp(:, :, ispin, ik)))))
            IF (hermitian_error > 1.0E-10_dp) &
               CPABORT("Folded g-xTB exchange matrix is not Hermitian")
         END DO
         IF (use_real_wfn .AND. &
             MAXVAL(ABS(AIMAG(tb%exchange_fock_kp(:, :, :, ik)))) > 1.0E-10_dp) &
            CPABORT("Complex g-xTB exchange cannot be represented by real k-point wavefunctions")
      END DO

      run_qualification_duality = .FALSE.
      IF (run_qualification_oracle) &
         run_qualification_duality = MAXVAL(ABS(density_reduced)) > 1.0E-14_dp
      IF (run_qualification_duality) THEN
         IF (.NOT. ALLOCATED(tb%qualification_acp_fock_full) .OR. &
             .NOT. ALLOCATED(tb%acp_fock_kp)) &
            CPABORT("QUALIFICATION_ONLY duality requires full and folded direct k-space ACP")
         CPASSERT(ALL(SHAPE(tb%qualification_acp_fock_full) == [nao, nao, kmap%nfull]))
         CPASSERT(ALL(SHAPE(tb%acp_fock_kp) == [nao, nao, nred]))
         duality_tolerance = 1.0E-10_dp
         exchange_duality_total_full = 0.0_dp
         exchange_duality_total_reduced = 0.0_dp
         acp_duality_total_full = 0.0_dp
         acp_duality_total_reduced = 0.0_dp
         DO ispin = 1, nspin
            CALL tb_qualification_gxtb_duality(kpoints, kmap, density(:, :, ispin, :), &
                                               density_reduced(:, :, ispin, :), &
                                               exchange_fock_full(:, :, ispin, :), &
                                               tb%exchange_fock_kp(:, :, ispin, :), &
                                               exchange_duality_full, exchange_duality_reduced, &
                                               exchange_duality_relative)
            CALL tb_qualification_gxtb_duality(kpoints, kmap, density(:, :, ispin, :), &
                                               density_reduced(:, :, ispin, :), &
                                               tb%qualification_acp_fock_full, tb%acp_fock_kp, &
                                               acp_duality_full, acp_duality_reduced, &
                                               acp_duality_relative)
            exchange_duality_total_full = exchange_duality_total_full + exchange_duality_full
            exchange_duality_total_reduced = exchange_duality_total_reduced + &
                                             exchange_duality_reduced
            acp_duality_total_full = acp_duality_total_full + acp_duality_full
            acp_duality_total_reduced = acp_duality_total_reduced + acp_duality_reduced
            IF (para_env%is_source()) THEN
               WRITE (*, '(A,I0,A,I0,6(A,ES14.6))') &
                  " GXTB-QUALIFICATION_ONLY DUALITY iter=", iter_count, " spin=", ispin, &
                  " exFull=", exchange_duality_full, " exRed=", exchange_duality_reduced, &
                  " exRel=", exchange_duality_relative, " acpFull=", acp_duality_full, &
                  " acpRed=", acp_duality_reduced, " acpRel=", acp_duality_relative
            END IF
            IF (exchange_duality_relative > duality_tolerance .OR. &
                acp_duality_relative > duality_tolerance) &
               CPABORT("QUALIFICATION_ONLY spin-resolved full/reduced duality exceeded 1.0E-10")
         END DO
         exchange_duality_relative = ABS(exchange_duality_total_full - &
                                         exchange_duality_total_reduced)/ &
                                     MAX(1.0E-30_dp, ABS(exchange_duality_total_full), &
                                         ABS(exchange_duality_total_reduced))
         acp_duality_relative = ABS(acp_duality_total_full - acp_duality_total_reduced)/ &
                                MAX(1.0E-30_dp, ABS(acp_duality_total_full), &
                                    ABS(acp_duality_total_reduced))
         IF (para_env%is_source()) THEN
            WRITE (*, '(A,I0,A,6(A,ES14.6))') &
               " GXTB-QUALIFICATION_ONLY DUALITY iter=", iter_count, " total", &
               " exFull=", exchange_duality_total_full, &
               " exRed=", exchange_duality_total_reduced, &
               " exRel=", exchange_duality_relative, " acpFull=", acp_duality_total_full, &
               " acpRed=", acp_duality_total_reduced, " acpRel=", acp_duality_relative
         END IF
         IF (exchange_duality_relative > duality_tolerance .OR. &
             acp_duality_relative > duality_tolerance) &
            CPABORT("QUALIFICATION_ONLY total full/reduced duality exceeded 1.0E-10")
      END IF

      ! QUALIFICATION_ONLY: this is deliberately dormant unless an exact SCC iteration is named
      ! through the environment.  It re-evaluates the already star-expanded density on a reversed
      ! complete mesh, restores the original ordering, folds independently, and applies strict
      ! numerical gates.  Production calculations must never set this variable; the hook is
      ! removed again after qualification artifacts have been frozen.
      IF (run_qualification_oracle) THEN
         ALLOCATE (density_oracle(nao, nao, nspin, kmap%nfull), &
                   overlap_oracle(nao, nao, kmap%nfull), &
                   exchange_fock_oracle(nao, nao, nspin, kmap%nfull), &
                   exchange_fock_unpermuted(nao, nao, nspin, kmap%nfull), &
                   folded_fock_oracle(nao, nao, nspin, nred), &
                   exchange_vsh_oracle(SIZE(exchange_vsh)), &
                   weight_oracle(kmap%nfull), xkp_oracle(3, kmap%nfull))
         DO ik = 1, kmap%nfull
            source_index = kmap%nfull + 1 - ik
            density_oracle(:, :, :, ik) = density(:, :, :, source_index)
            overlap_oracle(:, :, ik) = overlap(:, :, source_index)
            weight_oracle(ik) = kmap%weight(source_index)
            xkp_oracle(:, ik) = kmap%xkp(:, source_index)
         END DO
         CALL cp2k_exchange_kmesh(tb%calc, tb%mol, tb%ecache, tb%wfn, kmap%nmesh, &
                                  xkp_oracle, weight_oracle, density_oracle, overlap_oracle, &
                                  exchange_fock_oracle, exchange_vsh_oracle, &
                                  exchange_energy_oracle, error, &
                                  transform_backend=cp2k_exchange_transform_dense)
         IF (ALLOCATED(error)) &
            CPABORT("QUALIFICATION_ONLY permuted full-mesh exchange evaluation failed")
         DO ik = 1, kmap%nfull
            source_index = kmap%nfull + 1 - ik
            exchange_fock_unpermuted(:, :, :, source_index) = exchange_fock_oracle(:, :, :, ik)
         END DO
         folded_fock_oracle = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         oracle_sym_full = 0.0_dp
         oracle_sym_fold = 0.0_dp
         DO ispin = 1, nspin
            CALL tb_validate_gxtb_full_response(kpoints, kmap, ao_first, ao_size, &
                                                exchange_fock_unpermuted(:, :, ispin, :), &
                                                "QUALIFICATION_ONLY oracle Fock", symmetry_residual)
            oracle_sym_full = MAX(oracle_sym_full, symmetry_residual)
            CALL tb_fold_gxtb_kpoint_matrix(kpoints, kmap, ao_first, ao_size, &
                                            exchange_fock_unpermuted(:, :, ispin, :), &
                                            folded_fock_oracle(:, :, ispin, :))
            CALL tb_validate_gxtb_folded_response(kpoints, kmap, ao_first, ao_size, &
                                                  exchange_fock_unpermuted(:, :, ispin, :), &
                                                  folded_fock_oracle(:, :, ispin, :), &
                                                  "QUALIFICATION_ONLY oracle Fock", symmetry_residual)
            oracle_sym_fold = MAX(oracle_sym_fold, symmetry_residual)
         END DO
         oracle_energy_delta = ABS(exchange_energy_oracle - exchange_energy)
         oracle_vsh_delta = MAXVAL(ABS(exchange_vsh_oracle - exchange_vsh))
         oracle_full_abs = MAXVAL(ABS(exchange_fock_unpermuted - exchange_fock_full))
         oracle_scale = MAX(1.0_dp, MAXVAL(ABS(exchange_fock_full)))
         oracle_full_rel = oracle_full_abs/oracle_scale
         oracle_fold_abs = MAXVAL(ABS(folded_fock_oracle - tb%exchange_fock_kp))
         oracle_scale = MAX(1.0_dp, MAXVAL(ABS(tb%exchange_fock_kp)))
         oracle_fold_rel = oracle_fold_abs/oracle_scale
         oracle_tolerance = 1.0E-10_dp
         IF (para_env%is_source()) THEN
            WRITE (*, '(A,I0,10(A,ES14.6))') &
               " GXTB-QUALIFICATION_ONLY FULLMESH-ORACLE iter=", iter_count, &
               " dE=", oracle_energy_delta, " dVsh=", oracle_vsh_delta, &
               " dFfull=", oracle_full_abs, " rFfull=", oracle_full_rel, &
               " dFfold=", oracle_fold_abs, " rFfold=", oracle_fold_rel, &
               " symPfull=", primary_sym_full, " symPfold=", primary_sym_fold, &
               " symOfull=", oracle_sym_full, " symOfold=", oracle_sym_fold
         END IF
         IF (oracle_energy_delta > oracle_tolerance .OR. &
             oracle_vsh_delta > oracle_tolerance .OR. &
             oracle_full_rel > oracle_tolerance .OR. oracle_fold_rel > oracle_tolerance .OR. &
             primary_sym_full > oracle_tolerance .OR. primary_sym_fold > oracle_tolerance .OR. &
             oracle_sym_full > oracle_tolerance .OR. oracle_sym_fold > oracle_tolerance) &
            CPABORT("QUALIFICATION_ONLY permuted full-mesh oracle exceeded its 1.0E-10 gate")
         DEALLOCATE (density_oracle, overlap_oracle, exchange_fock_oracle, &
                     exchange_fock_unpermuted, folded_fock_oracle, exchange_vsh_oracle, &
                     weight_oracle, xkp_oracle)
      END IF

      ! F(k) is unweighted.  The energy and shell response returned by the
      ! whole-mesh API already contain the complete BZ contraction.
      tb%e_exchange = 0.0_dp
      tb%e_exchange(1) = exchange_energy
      tb%pot%vsh(:, 1) = tb%pot%vsh(:, 1) + exchange_vsh

      DEALLOCATE (ao_first, ao_size, density, exchange_fock_full, overlap, exchange_vsh)
      IF (ALLOCATED(density_reduced)) DEALLOCATE (density_reduced)
      tb%exchange_fock_kp_valid = .TRUE.

   END SUBROUTINE tb_build_gxtb_kpoint_exchange

! **************************************************************************************************
!> \brief Dispatch the coupled g-xTB exchange derivative to the dense oracle or bounded consumer.
!> \param qs_env ...
!> \param tb ...
!> \param gradient ...
!> \param sigma ...
! **************************************************************************************************
   SUBROUTINE tb_gxtb_kpoint_exchange_gradient(qs_env, tb, gradient, sigma)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: gradient, sigma

      INTEGER                                            :: gradient_mode

      gradient_mode = tb_gxtb_exchange_gradient_mode()
      SELECT CASE (gradient_mode)
      CASE (tb_gxtb_exchange_gradient_dense)
         CALL tb_gxtb_kpoint_exchange_gradient_dense(qs_env, tb, gradient, sigma)
      CASE (tb_gxtb_exchange_gradient_streamed)
         CALL tb_gxtb_kpoint_exchange_gradient_streamed(qs_env, tb, gradient, sigma, .FALSE.)
      CASE (tb_gxtb_exchange_gradient_qualify)
         CALL tb_gxtb_kpoint_exchange_gradient_streamed(qs_env, tb, gradient, sigma, .TRUE.)
      CASE DEFAULT
         CPABORT("Invalid g-xTB exchange-gradient bridge")
      END SELECT

   END SUBROUTINE tb_gxtb_kpoint_exchange_gradient

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> \param gradient ...
!> \param sigma ...
!> Differentiate the explicit g-xTB exchange contribution on a CP2K k-point mesh.
!>
!> A complete multi-k mesh is differentiated by one coupled save_tblite reverse sweep.  Its direct
!> atomic gradient and strain derivative are already complete Brillouin-zone contractions per
!> primitive cell.  The returned unweighted overlap response transforms like a density matrix and
!> is folded back to CP2K's real-space image convention by kpoint_density_transform, which applies
!> the k-point weights, Bloch phases, and image orientation exactly once.
!>
!> This isolated bridge puts the spin-summed exchange overlap response into matrix_w(1,:).  The
!> ordinary alpha/beta channels must still be separate at entry; a later up/down-to-charge/magnet
!> conversion adds matrix_w(2,:) to channel one.  H0/ACP response and final force/virial dispatch
!> deliberately remain outside.  The established independent-k derivative remains the one-point
!> fallback so the Gamma path retains its existing response and symmetry handling.
! **************************************************************************************************
   SUBROUTINE tb_gxtb_kpoint_exchange_gradient_dense(qs_env, tb, gradient, sigma)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: gradient, sigma

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: overlap, overlap_grad, overlap_grad_full
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: density
      INTEGER                                            :: ic, ik, ikk, img, ispin, isym, kplocal, &
                                                            nao, nat, nc, nimg, nred, nspin, &
                                                            source_atom, target_atom, &
                                                            transform_backend
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_first, ao_size
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: full_grid, use_real_wfn, validate_density
      REAL(KIND=dp)                                      :: hermitian_error, star_weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: k_gradient, k_sigma, overlap_part
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rotated_sigma
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_ao_fm_pools
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fmwork
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: pmat_ext
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: exchange_w_rs, matrix_s, matrix_w
      TYPE(dbcsr_type), POINTER                          :: cmatrix, rmatrix
      TYPE(error_type), ALLOCATABLE                      :: error
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_sym_type), POINTER                     :: kpsym
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(tb_gxtb_kmesh_map_type), POINTER              :: kmap

      NULLIFY (ao_ao_fm_pools, cell_to_index, cmatrix, exchange_w_rs, kmap, kp, kpoints, &
               kpsym, matrix_s, matrix_struct, matrix_w, para_env, rmatrix, sab_nl)
      gradient = 0.0_dp
      sigma = 0.0_dp

      CALL get_qs_env(qs_env, kpoints=kpoints, matrix_s_kp=matrix_s, matrix_w_kp=matrix_w, &
                      para_env=para_env)
      IF (.NOT. ASSOCIATED(kpoints)) CPABORT("Missing k-point environment for g-xTB exchange gradient")
      IF (.NOT. ASSOCIATED(matrix_s)) CPABORT("Missing overlap images for g-xTB exchange gradient")
      IF (.NOT. ASSOCIATED(matrix_w)) &
         CPABORT("Missing energy-weighted density images for g-xTB exchange gradient")
      IF (.NOT. ASSOCIATED(para_env)) &
         CPABORT("Missing parallel environment for g-xTB exchange gradient")
      CALL mpools_get(qs_env%mpools, ao_ao_fm_pools=ao_ao_fm_pools)
      IF (.NOT. ASSOCIATED(ao_ao_fm_pools)) &
         CPABORT("Missing AO full-matrix pool for g-xTB exchange gradient")
      CPASSERT(SIZE(ao_ao_fm_pools) >= 1)
      IF (.NOT. ASSOCIATED(ao_ao_fm_pools(1)%pool)) &
         CPABORT("Uninitialized AO full-matrix pool for g-xTB exchange gradient")

      CALL get_kpoint_info(kpoints, nkp=nred, kp_range=kp_range, &
                           full_grid=full_grid, use_real_wfn=use_real_wfn, &
                           cell_to_index=cell_to_index, sab_nl=sab_nl)
      IF (.NOT. ASSOCIATED(cell_to_index) .OR. .NOT. ASSOCIATED(sab_nl)) &
         CPABORT("Incomplete k-point image map for g-xTB exchange gradient")
      CALL tb_prepare_gxtb_kmesh_map(tb, kpoints, matrix_s, sab_nl, cell_to_index)
      kmap => tb%gxtb_kmesh_map
      transform_backend = tb_gxtb_exchange_transform_backend()

      nao = tb%calc%bas%nao
      nat = tb%mol%nat
      nspin = tb%wfn%nspin
      nimg = SIZE(matrix_s, 2)
      nc = MERGE(1, 2, use_real_wfn)
      kplocal = kp_range(2) - kp_range(1) + 1
      CALL tb_gxtb_ao_layout(tb, ao_first, ao_size)
      CPASSERT(SIZE(gradient, 1) == 3 .AND. SIZE(gradient, 2) == nat)
      CPASSERT(SIZE(sigma, 1) == 3 .AND. SIZE(sigma, 2) == 3)
      CPASSERT(SIZE(matrix_s, 1) >= 1)
      CPASSERT(SIZE(matrix_w, 1) >= nspin .AND. SIZE(matrix_w, 2) == nimg)
      CPASSERT(kplocal > 0)

      ALLOCATE (density(nao, nao, nspin, kmap%nfull), overlap(nao, nao, kmap%nfull), &
                overlap_grad_full(nao, nao, kmap%nfull), &
                overlap_part(nao, nao), &
                k_gradient(3, nat), k_sigma(3, 3))
      overlap_grad_full = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      IF (kmap%identity) THEN
         CALL tb_collect_kpoint_density(qs_env, tb, kpoints, density)
      ELSE
         CALL tb_collect_kpoint_density(qs_env, tb, kpoints, density, kmap%representative)
         validate_density = .NOT. tb%gxtb_density_map_validated .AND. &
                            MAXVAL(ABS(density)) > 1.0E-14_dp
         CALL tb_expand_gxtb_kpoint_density(kpoints, kmap, ao_first, ao_size, density, &
                                            validate_density)
         IF (validate_density) tb%gxtb_density_map_validated = .TRUE.
      END IF
      IF (.NOT. ALL(ieee_is_finite(REAL(density, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(density)))) &
         CPABORT("CP2K supplied a non-finite k-point density to the exchange gradient")

      ALLOCATE (rmatrix, cmatrix)
      CALL dbcsr_create(rmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_nl)

      DO ik = 1, kmap%nfull
         CALL dbcsr_set(rmatrix, 0.0_dp)
         CALL dbcsr_set(cmatrix, 0.0_dp)
         CALL rskp_transform(rmatrix, cmatrix, matrix_s, 1, kmap%xkp(1:3, ik), &
                             cell_to_index, sab_nl)
         CALL tb_dbcsr_component_to_dense(tb, rmatrix, overlap_part, para_env, .FALSE.)
         IF (.NOT. ALL(ieee_is_finite(overlap_part))) &
            CPABORT("CP2K supplied a non-finite overlap to the exchange gradient")
         overlap(:, :, ik) = CMPLX(overlap_part, 0.0_dp, KIND=dp)
         CALL tb_dbcsr_component_to_dense(tb, cmatrix, overlap_part, para_env, .TRUE.)
         IF (.NOT. ALL(ieee_is_finite(overlap_part))) &
            CPABORT("CP2K supplied a non-finite overlap to the exchange gradient")
         overlap(:, :, ik) = CMPLX(REAL(overlap(:, :, ik), KIND=dp), overlap_part, KIND=dp)
         hermitian_error = MAXVAL(ABS(overlap(:, :, ik) - &
                                      CONJG(TRANSPOSE(overlap(:, :, ik)))))
         IF (hermitian_error > 1.0E-10_dp) &
            CPABORT("CP2K supplied a non-Hermitian overlap to the exchange gradient")
         DO ispin = 1, nspin
            hermitian_error = MAXVAL(ABS(density(:, :, ispin, ik) - &
                                         CONJG(TRANSPOSE(density(:, :, ispin, ik)))))
            IF (hermitian_error > 1.0E-10_dp) &
               CPABORT("CP2K supplied a non-Hermitian density to the exchange gradient")
         END DO
      END DO

      IF (.NOT. kmap%identity .AND. .NOT. tb%gxtb_transform_validated) THEN
         CALL tb_validate_gxtb_kmesh_transform(kpoints, kmap, ao_first, ao_size, overlap)
         tb%gxtb_transform_validated = .TRUE.
      END IF
      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)
      DEALLOCATE (overlap_part)

      IF (kmap%nfull == 1) THEN
         ! Preserve the established Gamma/one-point response, including its existing symmetry
         ! expansion.  Both the response and the direct derivatives are unweighted here.
         CALL cp2k_exchange_kpoint_gradient(tb%calc, tb%mol, tb%ecache, &
                                            density(:, :, :, 1), overlap(:, :, 1), &
                                            overlap_grad_full(:, :, 1), &
                                            k_gradient, k_sigma, error)
         IF (ALLOCATED(error)) CPABORT("save_tblite k-point exchange gradient failed")
         ik = 1
         IF (full_grid) THEN
            gradient = gradient + kpoints%wkp(ik)*k_gradient
            sigma = sigma + kpoints%wkp(ik)*k_sigma
         ELSE
            kpsym => kpoints%kp_sym(ik)%kpoint_sym
            CPASSERT(ASSOCIATED(kpsym))
            IF (kpsym%nwred > 0) THEN
               CPASSERT(ASSOCIATED(kpsym%f0))
               CPASSERT(ASSOCIATED(kpsym%rot))
               CPASSERT(ASSOCIATED(kpsym%rotp))
               CPASSERT(SIZE(kpsym%f0, 1) == nat)
               CPASSERT(SIZE(kpsym%f0, 2) >= kpsym%nwred)
               CPASSERT(SIZE(kpsym%rot, 3) >= kpsym%nwred)
               CPASSERT(SIZE(kpsym%rotp) >= kpsym%nwred)
               star_weight = kpoints%wkp(ik)/REAL(kpsym%nwred, KIND=dp)
               DO isym = 1, kpsym%nwred
                  ! A negative rotp marks time reversal.  These real nuclear and strain
                  ! derivatives transform with the spatial rotation only: there is no extra sign.
                  DO source_atom = 1, nat
                     target_atom = kpsym%f0(source_atom, isym)
                     CPASSERT(target_atom >= 1 .AND. target_atom <= nat)
                     gradient(:, target_atom) = gradient(:, target_atom) + star_weight* &
                                                MATMUL(kpsym%rot(1:3, 1:3, isym), &
                                                       k_gradient(:, source_atom))
                  END DO
                  rotated_sigma = MATMUL(kpsym%rot(1:3, 1:3, isym), &
                                         MATMUL(k_sigma, &
                                                TRANSPOSE(kpsym%rot(1:3, 1:3, isym))))
                  sigma = sigma + star_weight*rotated_sigma
               END DO
            ELSE
               ! Inversion/time-reversal-only reduction has equal real derivatives at k and -k;
               ! its representative already carries the combined CP2K mesh weight.
               gradient = gradient + kpoints%wkp(ik)*k_gradient
               sigma = sigma + kpoints%wkp(ik)*k_sigma
            END IF
         END IF
      ELSE
         ! Expected save_tblite ABI:
         ! cp2k_exchange_kmesh_gradient(calc, mol, ecache, nmesh, kfrac, weights, density,
         !                              overlap, overlap_grad, gradient, sigma, error)
         ! overlap_grad is spin-summed and unweighted; k_gradient/k_sigma are already total
         ! primitive-cell responses and therefore must not receive another k-point weight here.
         CALL cp2k_exchange_kmesh_gradient(tb%calc, tb%mol, tb%ecache, kmap%nmesh, &
                                           kmap%xkp, kmap%weight, density, overlap, overlap_grad_full, &
                                           k_gradient, k_sigma, error, &
                                           transform_backend=transform_backend)
         IF (ALLOCATED(error)) CPABORT("save_tblite whole-mesh exchange gradient failed")
         gradient = k_gradient
         sigma = k_sigma
      END IF
      DEALLOCATE (density, overlap)

      IF (.NOT. ALL(ieee_is_finite(REAL(overlap_grad_full, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(overlap_grad_full))) .OR. &
          .NOT. ALL(ieee_is_finite(k_gradient)) .OR. &
          .NOT. ALL(ieee_is_finite(k_sigma))) &
         CPABORT("save_tblite returned a non-finite exchange gradient")
      DO ik = 1, kmap%nfull
         hermitian_error = MAXVAL(ABS(overlap_grad_full(:, :, ik) - &
                                      CONJG(TRANSPOSE(overlap_grad_full(:, :, ik)))))
         IF (hermitian_error > 1.0E-10_dp) &
            CPABORT("save_tblite returned a non-Hermitian overlap gradient")
      END DO
      IF (.NOT. kmap%identity) &
         CALL tb_validate_gxtb_full_response(kpoints, kmap, ao_first, ao_size, &
                                             overlap_grad_full, "overlap gradient")
      ALLOCATE (overlap_grad(nao, nao, nred), source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      CALL tb_fold_gxtb_kpoint_matrix(kpoints, kmap, ao_first, ao_size, &
                                      overlap_grad_full, overlap_grad)
      IF (.NOT. kmap%identity) &
         CALL tb_validate_gxtb_folded_response(kpoints, kmap, ao_first, ao_size, &
                                               overlap_grad_full, overlap_grad, "overlap gradient")

      ! pmat_ext is local to each k-point group, just like kp%pmat. Allocate it only after the
      ! provider reverse sweep, so it does not add to the unavoidable full-mesh API peak.
      ALLOCATE (pmat_ext(kplocal, nc, nspin))
      DO ikk = 1, kplocal
         kp => kpoints%kp_env(ikk)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(SIZE(kp%pmat, 1) >= nc .AND. SIZE(kp%pmat, 2) >= nspin)
         DO ispin = 1, nspin
            DO ic = 1, nc
               CALL cp_fm_get_info(kp%pmat(ic, ispin), matrix_struct=matrix_struct)
               CPASSERT(ASSOCIATED(matrix_struct))
               CALL cp_fm_create(pmat_ext(ikk, ic, ispin), matrix_struct, set_zero=.TRUE.)
            END DO
         END DO
      END DO
      DO ik = 1, nred
         hermitian_error = MAXVAL(ABS(overlap_grad(:, :, ik) - &
                                      CONJG(TRANSPOSE(overlap_grad(:, :, ik)))))
         IF (hermitian_error > 1.0E-10_dp) &
            CPABORT("Folded g-xTB overlap gradient is not Hermitian")
         IF (use_real_wfn .AND. MAXVAL(ABS(AIMAG(overlap_grad(:, :, ik)))) > 1.0E-10_dp) &
            CPABORT("Complex exchange overlap gradient cannot use real k-point wavefunctions")

         IF (ik >= kp_range(1) .AND. ik <= kp_range(2)) THEN
            ikk = ik - kp_range(1) + 1
            CALL tb_dense_to_fm(REAL(overlap_grad(:, :, ik), KIND=dp), pmat_ext(ikk, 1, 1))
            IF (.NOT. use_real_wfn) &
               CALL tb_dense_to_fm(AIMAG(overlap_grad(:, :, ik)), pmat_ext(ikk, 2, 1))
         END IF
      END DO

      CALL dbcsr_allocate_matrix_set(exchange_w_rs, nspin, nimg)
      DO ispin = 1, nspin
         DO img = 1, nimg
            ALLOCATE (exchange_w_rs(ispin, img)%matrix)
            CALL dbcsr_create(exchange_w_rs(ispin, img)%matrix, &
                              template=matrix_s(1, 1)%matrix, &
                              name="g-XTB EXCHANGE OVERLAP RESPONSE", &
                              matrix_type=dbcsr_type_symmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(exchange_w_rs(ispin, img)%matrix, sab_nl)
            CALL dbcsr_set(exchange_w_rs(ispin, img)%matrix, 0.0_dp)
         END DO
      END DO

      ! SCF work matrices are released before force evaluation.  Borrow a force-safe AO full-
      ! matrix workspace from the persistent pool for the inverse k-point transform instead.
      ALLOCATE (fmwork(nc))
      DO ic = 1, nc
         CALL fm_pool_create_fm(ao_ao_fm_pools(1)%pool, fmwork(ic))
      END DO
      CALL kpoint_density_transform(kpoints, exchange_w_rs, .TRUE., &
                                    matrix_s(1, 1)%matrix, sab_nl, fmwork, &
                                    pmat_ext=pmat_ext, overlap_rs=matrix_s)
      DO ic = 1, nc
         CALL fm_pool_give_back_fm(ao_ao_fm_pools(1)%pool, fmwork(ic))
      END DO
      DEALLOCATE (fmwork)
      DO img = 1, nimg
         ! Both APIs give the unweighted G_k=dE_exchange/dS.  At Gamma this is
         ! -ao_grad_gamma, while the established path adds ao_grad_gamma to W; hence the common
         ! image convention is W_exchange=-G_k.  kpoint_density_transform has already supplied the
         ! one and only mesh weight.  Avoiding a dense round trip preserves image orientation.
         CALL dbcsr_add(matrix_w(1, img)%matrix, exchange_w_rs(1, img)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      END DO

      CALL dbcsr_deallocate_matrix_set(exchange_w_rs)
      DO ispin = 1, nspin
         DO ic = 1, nc
            DO ikk = 1, kplocal
               CALL cp_fm_release(pmat_ext(ikk, ic, ispin))
            END DO
         END DO
      END DO
      DEALLOCATE (ao_first, ao_size, pmat_ext, overlap_grad, &
                  overlap_grad_full, k_gradient, k_sigma)

   END SUBROUTINE tb_gxtb_kpoint_exchange_gradient_dense

! **************************************************************************************************
!> \brief Memory-bounded reverse consumer for Brillouin-zone-coupled g-xTB exchange.
!> \param qs_env ...
!> \param tb ...
!> \param gradient ...
!> \param sigma ...
!> \param qualify compare with the unchanged complete-mesh derivative before applying the response
!>
!> CP2K retains dense irreducible density and response arrays, requiring O(nred*nao**2) host
!> storage; for FULL_GRID, nred equals the complete-mesh size.  It does not additionally retain an
!> expanded full-star density, overlap, or response mesh.  Instead, every full-star density and
!> overlap block is reconstructed for the current provider image batch, and every returned overlap
!> adjoint is folded immediately with the calibrated CP2K symmetry adjoint.  The provider's AO-image
!> workspace is O(batch_size*nao**2*nspin).  The direct atomic and homogeneous-strain derivatives are
!> retrieved once after all image batches.  QUALIFY alone allocates additional complete meshes for
!> the dense numerical oracle.
! **************************************************************************************************
   SUBROUTINE tb_gxtb_kpoint_exchange_gradient_streamed(qs_env, tb, gradient, sigma, qualify)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: gradient, sigma
      LOGICAL, INTENT(IN)                                :: qualify

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: density_reduced
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: density_block, overlap_grad
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: overlap_adjoint, overlap_block, &
                                                            transform_work
      INTEGER                                            :: actual_batch_size, batch_size, ic, ifull, &
                                                            ik, ikk, img, ispin, kplocal, nao, nat, &
                                                            nc, nimg, nred, nspin
      INTEGER(KIND=int_8)                                :: expected_peak, peak_complex
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_first, ao_size
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: batch_done, use_real_wfn
      REAL(KIND=dp)                                      :: hermitian_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: k_gradient, k_sigma, overlap_part
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_ao_fm_pools
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fmwork
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: pmat_ext
      TYPE(cp2k_exchange_stream_type)                    :: stream
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: exchange_w_rs, matrix_s, matrix_w
      TYPE(dbcsr_type), POINTER                          :: cmatrix, rmatrix
      TYPE(error_type), ALLOCATABLE                      :: error
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(tb_gxtb_kmesh_map_type), POINTER              :: kmap

      NULLIFY (ao_ao_fm_pools, cell_to_index, cmatrix, exchange_w_rs, kmap, kp, kpoints, &
               matrix_s, matrix_struct, matrix_w, para_env, rmatrix, sab_nl)
      gradient = 0.0_dp
      sigma = 0.0_dp

      CALL get_qs_env(qs_env, kpoints=kpoints, matrix_s_kp=matrix_s, matrix_w_kp=matrix_w, &
                      para_env=para_env)
      IF (.NOT. ASSOCIATED(kpoints)) &
         CPABORT("Missing k-point environment for streamed g-xTB exchange gradient")
      IF (.NOT. ASSOCIATED(matrix_s)) &
         CPABORT("Missing overlap images for streamed g-xTB exchange gradient")
      IF (.NOT. ASSOCIATED(matrix_w)) &
         CPABORT("Missing energy-weighted density images for streamed g-xTB exchange gradient")
      IF (.NOT. ASSOCIATED(para_env)) &
         CPABORT("Missing parallel environment for streamed g-xTB exchange gradient")
      CALL mpools_get(qs_env%mpools, ao_ao_fm_pools=ao_ao_fm_pools)
      IF (.NOT. ASSOCIATED(ao_ao_fm_pools)) &
         CPABORT("Missing AO full-matrix pool for streamed g-xTB exchange gradient")
      CPASSERT(SIZE(ao_ao_fm_pools) >= 1)
      IF (.NOT. ASSOCIATED(ao_ao_fm_pools(1)%pool)) &
         CPABORT("Uninitialized AO full-matrix pool for streamed g-xTB exchange gradient")

      CALL get_kpoint_info(kpoints, nkp=nred, kp_range=kp_range, use_real_wfn=use_real_wfn, &
                           cell_to_index=cell_to_index, sab_nl=sab_nl)
      IF (.NOT. ASSOCIATED(cell_to_index) .OR. .NOT. ASSOCIATED(sab_nl)) &
         CPABORT("Incomplete k-point image map for streamed g-xTB exchange gradient")
      CALL tb_prepare_gxtb_kmesh_map(tb, kpoints, matrix_s, sab_nl, cell_to_index)
      kmap => tb%gxtb_kmesh_map

      ! Keep the established one-point derivative, including its special CP2K symmetry handling.
      IF (kmap%nfull == 1) THEN
         CALL tb_gxtb_kpoint_exchange_gradient_dense(qs_env, tb, gradient, sigma)
         RETURN
      END IF

      nao = tb%calc%bas%nao
      nat = tb%mol%nat
      nspin = tb%wfn%nspin
      nimg = SIZE(matrix_s, 2)
      nc = MERGE(1, 2, use_real_wfn)
      kplocal = kp_range(2) - kp_range(1) + 1
      CALL tb_gxtb_ao_layout(tb, ao_first, ao_size)
      CPASSERT(SIZE(gradient, 1) == 3 .AND. SIZE(gradient, 2) == nat)
      CPASSERT(SIZE(sigma, 1) == 3 .AND. SIZE(sigma, 2) == 3)
      CPASSERT(SIZE(matrix_s, 1) >= 1)
      CPASSERT(SIZE(matrix_w, 1) >= nspin .AND. SIZE(matrix_w, 2) == nimg)
      CPASSERT(kplocal > 0)

      ALLOCATE (density_reduced(nao, nao, nspin, nred), density_block(nao, nao, nspin), &
                overlap_adjoint(nao, nao), overlap_block(nao, nao), &
                overlap_grad(nao, nao, nred), transform_work(nao, nao), &
                source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      ALLOCATE (overlap_part(nao, nao), k_gradient(3, nat), k_sigma(3, 3), source=0.0_dp)
      CALL tb_collect_kpoint_density(qs_env, tb, kpoints, density_reduced)
      IF (.NOT. ALL(ieee_is_finite(REAL(density_reduced, KIND=dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(density_reduced)))) &
         CPABORT("CP2K supplied a non-finite irreducible density to streamed save_tblite")

      ALLOCATE (rmatrix, cmatrix)
      CALL dbcsr_create(rmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_nl)

      batch_size = tb_gxtb_exchange_stream_batch_size()
      ! A bounded qualification must not silently become one complete AO image mesh.
      batch_size = MIN(batch_size, kmap%nfull - 1)
      actual_batch_size = MIN(batch_size, kmap%nfull)
      CALL cp2k_exchange_stream_begin(stream, tb%calc, tb%mol, tb%ecache, tb%wfn, &
                                      kmap%nmesh, kmap%xkp, kmap%weight, error, &
                                      mode=cp2k_exchange_stream_batched, batch_size=batch_size, &
                                      reverse_only=.TRUE.)
      IF (ALLOCATED(error)) &
         CPABORT("save_tblite bounded reverse begin failed: "//TRIM(error%message))
      IF (cp2k_exchange_stream_has_full_mesh_storage(stream)) &
         CPABORT("Bounded g-xTB reverse stream unexpectedly owns a complete AO mesh")
      peak_complex = cp2k_exchange_stream_peak_complex_elements(stream)
      expected_peak = 6_int_8*INT(actual_batch_size, int_8)*INT(nao, int_8)**2* &
                      INT(nspin, int_8) + 5_int_8*INT(nao, int_8)*INT(nspin, int_8) + &
                      2_int_8*INT(kmap%nfull, int_8) + &
                      3_int_8*INT(MAXVAL(kmap%nmesh), int_8) + &
                      10_int_8*INT(nao, int_8)**2
      IF (peak_complex /= expected_peak) &
         CPABORT("save_tblite bounded reverse high-water counter violated its exact formula")

      batch_done = .FALSE.
      DO WHILE (.NOT. batch_done)
         DO ifull = 1, kmap%nfull
            CALL prepare_reverse_block(ifull)
            CALL cp2k_exchange_stream_reverse_push(stream, ifull, density_block, &
                                                   overlap_block, error)
            IF (ALLOCATED(error)) &
               CPABORT("save_tblite bounded reverse push failed: "//TRIM(error%message))
         END DO
         CALL cp2k_exchange_stream_reverse_batch_apply(stream, tb%calc, tb%mol, tb%ecache, error)
         IF (ALLOCATED(error)) &
            CPABORT("save_tblite bounded reverse apply failed: "//TRIM(error%message))

         DO ifull = 1, kmap%nfull
            CALL prepare_reverse_block(ifull)
            overlap_adjoint = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            CALL cp2k_exchange_stream_reverse_batch_pull(stream, ifull, density_block, &
                                                         overlap_block, overlap_adjoint, error)
            IF (ALLOCATED(error)) &
               CPABORT("save_tblite bounded reverse pull failed: "//TRIM(error%message))
            IF (.NOT. ALL(ieee_is_finite(REAL(overlap_adjoint, KIND=dp))) .OR. &
                .NOT. ALL(ieee_is_finite(AIMAG(overlap_adjoint)))) &
               CPABORT("save_tblite returned a non-finite bounded overlap adjoint")
            hermitian_error = MAXVAL(ABS(overlap_adjoint - CONJG(TRANSPOSE(overlap_adjoint))))
            IF (hermitian_error > 1.0E-10_dp) &
               CPABORT("save_tblite returned a non-Hermitian bounded overlap adjoint")
            ik = kmap%source(ifull)
            CALL tb_accumulate_gxtb_kpoint_response(kpoints, kmap, ao_first, ao_size, ifull, &
                                                    overlap_adjoint, overlap_grad(:, :, ik), &
                                                    transform_work)
         END DO
         CALL cp2k_exchange_stream_reverse_batch_advance(stream, batch_done, error)
         IF (ALLOCATED(error)) &
            CPABORT("save_tblite bounded reverse advance failed: "//TRIM(error%message))
      END DO
      CALL cp2k_exchange_stream_reverse_result(stream, tb%calc, tb%mol, tb%ecache, &
                                               k_gradient, k_sigma, error)
      IF (ALLOCATED(error)) &
         CPABORT("save_tblite bounded reverse result failed: "//TRIM(error%message))
      CALL cp2k_exchange_stream_end(stream, error)
      IF (ALLOCATED(error)) &
         CPABORT("save_tblite bounded reverse end failed: "//TRIM(error%message))
      gradient = k_gradient
      sigma = k_sigma

      DO ik = 1, nred
         overlap_grad(:, :, ik) = 0.5_dp*(overlap_grad(:, :, ik) + &
                                          CONJG(TRANSPOSE(overlap_grad(:, :, ik))))
         hermitian_error = MAXVAL(ABS(overlap_grad(:, :, ik) - &
                                      CONJG(TRANSPOSE(overlap_grad(:, :, ik)))))
         IF (hermitian_error > 1.0E-10_dp) &
            CPABORT("Folded bounded g-xTB overlap gradient is not Hermitian")
         IF (use_real_wfn .AND. MAXVAL(ABS(AIMAG(overlap_grad(:, :, ik)))) > 1.0E-10_dp) &
            CPABORT("Complex bounded exchange overlap gradient cannot use real wavefunctions")
      END DO
      IF (.NOT. ALL(ieee_is_finite(k_gradient)) .OR. .NOT. ALL(ieee_is_finite(k_sigma))) &
         CPABORT("save_tblite returned a non-finite bounded direct exchange derivative")

      IF (qualify) CALL qualify_bounded_reverse()
      IF (para_env%is_source()) THEN
         WRITE (*, '(A,I0,A,I0,A,I0)') " GXTB STREAMED-REVERSE nfull=", kmap%nfull, &
            " batch=", actual_batch_size, " peak_complex=", peak_complex
      END IF

      ! Convert only the irreducible response to CP2K's force-safe full matrices.
      ALLOCATE (pmat_ext(kplocal, nc, nspin))
      DO ikk = 1, kplocal
         kp => kpoints%kp_env(ikk)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         CPASSERT(SIZE(kp%pmat, 1) >= nc .AND. SIZE(kp%pmat, 2) >= nspin)
         DO ispin = 1, nspin
            DO ic = 1, nc
               CALL cp_fm_get_info(kp%pmat(ic, ispin), matrix_struct=matrix_struct)
               CPASSERT(ASSOCIATED(matrix_struct))
               CALL cp_fm_create(pmat_ext(ikk, ic, ispin), matrix_struct, set_zero=.TRUE.)
            END DO
         END DO
      END DO
      DO ik = 1, nred
         IF (ik >= kp_range(1) .AND. ik <= kp_range(2)) THEN
            ikk = ik - kp_range(1) + 1
            CALL tb_dense_to_fm(REAL(overlap_grad(:, :, ik), KIND=dp), pmat_ext(ikk, 1, 1))
            IF (.NOT. use_real_wfn) &
               CALL tb_dense_to_fm(AIMAG(overlap_grad(:, :, ik)), pmat_ext(ikk, 2, 1))
         END IF
      END DO

      CALL dbcsr_allocate_matrix_set(exchange_w_rs, nspin, nimg)
      DO ispin = 1, nspin
         DO img = 1, nimg
            ALLOCATE (exchange_w_rs(ispin, img)%matrix)
            CALL dbcsr_create(exchange_w_rs(ispin, img)%matrix, &
                              template=matrix_s(1, 1)%matrix, &
                              name="BOUNDED g-XTB EXCHANGE OVERLAP RESPONSE", &
                              matrix_type=dbcsr_type_symmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(exchange_w_rs(ispin, img)%matrix, sab_nl)
            CALL dbcsr_set(exchange_w_rs(ispin, img)%matrix, 0.0_dp)
         END DO
      END DO
      ALLOCATE (fmwork(nc))
      DO ic = 1, nc
         CALL fm_pool_create_fm(ao_ao_fm_pools(1)%pool, fmwork(ic))
      END DO
      CALL kpoint_density_transform(kpoints, exchange_w_rs, .TRUE., &
                                    matrix_s(1, 1)%matrix, sab_nl, fmwork, &
                                    pmat_ext=pmat_ext, overlap_rs=matrix_s)
      DO ic = 1, nc
         CALL fm_pool_give_back_fm(ao_ao_fm_pools(1)%pool, fmwork(ic))
      END DO
      DEALLOCATE (fmwork)
      DO img = 1, nimg
         CALL dbcsr_add(matrix_w(1, img)%matrix, exchange_w_rs(1, img)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      END DO

      CALL dbcsr_deallocate_matrix_set(exchange_w_rs)
      DO ispin = 1, nspin
         DO ic = 1, nc
            DO ikk = 1, kplocal
               CALL cp_fm_release(pmat_ext(ikk, ic, ispin))
            END DO
         END DO
      END DO
      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)
      DEALLOCATE (ao_first, ao_size, density_block, density_reduced, k_gradient, k_sigma, &
                  overlap_adjoint, overlap_block, overlap_grad, overlap_part, pmat_ext, &
                  transform_work)

   CONTAINS

! **************************************************************************************************
!> \brief Reconstruct one exact full-star density and overlap block.
!> \param ifull_local ...
! **************************************************************************************************
      SUBROUTINE prepare_reverse_block(ifull_local)

      INTEGER, INTENT(IN)                                :: ifull_local

      INTEGER                                            :: ispin_local, source_local
      REAL(KIND=dp)                                      :: local_hermitian_error

         source_local = kmap%source(ifull_local)
         DO ispin_local = 1, nspin
            CALL tb_expand_gxtb_kpoint_block(kpoints, kmap, ao_first, ao_size, ifull_local, &
                                             density_reduced(:, :, ispin_local, source_local), &
                                             density_block(:, :, ispin_local), transform_work)
            local_hermitian_error = &
               MAXVAL(ABS(density_block(:, :, ispin_local) - &
                          CONJG(TRANSPOSE(density_block(:, :, ispin_local)))))
            IF (local_hermitian_error > 1.0E-10_dp) &
               CPABORT("Streamed reverse g-xTB density block is not Hermitian")
         END DO

         CALL dbcsr_set(rmatrix, 0.0_dp)
         CALL dbcsr_set(cmatrix, 0.0_dp)
         CALL rskp_transform(rmatrix, cmatrix, matrix_s, 1, kmap%xkp(1:3, ifull_local), &
                             cell_to_index, sab_nl)
         CALL tb_dbcsr_component_to_dense(tb, rmatrix, overlap_part, para_env, .FALSE.)
         overlap_block(:, :) = CMPLX(overlap_part, 0.0_dp, KIND=dp)
         CALL tb_dbcsr_component_to_dense(tb, cmatrix, overlap_part, para_env, .TRUE.)
         overlap_block(:, :) = CMPLX(REAL(overlap_block, KIND=dp), overlap_part, KIND=dp)
         IF (.NOT. ALL(ieee_is_finite(REAL(overlap_block, KIND=dp))) .OR. &
             .NOT. ALL(ieee_is_finite(AIMAG(overlap_block)))) &
            CPABORT("CP2K supplied a non-finite overlap to streamed reverse save_tblite")
         local_hermitian_error = &
            MAXVAL(ABS(overlap_block - CONJG(TRANSPOSE(overlap_block))))
         IF (local_hermitian_error > 1.0E-10_dp) &
            CPABORT("CP2K supplied a non-Hermitian overlap to streamed reverse save_tblite")

      END SUBROUTINE prepare_reverse_block

! **************************************************************************************************
!> \brief Compare the bounded result with the unchanged complete-mesh provider derivative.
! **************************************************************************************************
      SUBROUTINE qualify_bounded_reverse()

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: overlap_full, overlap_grad_full, &
                                                            overlap_grad_oracle
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: density_full
      INTEGER                                            :: ifull_local
      REAL(KIND=dp)                                      :: force_residual, response_residual, &
                                                            response_scale, stress_residual
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gradient_oracle, sigma_oracle

         ALLOCATE (density_full(nao, nao, nspin, kmap%nfull), &
                   overlap_full(nao, nao, kmap%nfull), &
                   overlap_grad_full(nao, nao, kmap%nfull), &
                   overlap_grad_oracle(nao, nao, nred), &
                   source=CMPLX(0.0_dp, 0.0_dp, KIND=dp))
         ALLOCATE (gradient_oracle(3, nat), sigma_oracle(3, 3), source=0.0_dp)
         DO ifull_local = 1, kmap%nfull
            CALL prepare_reverse_block(ifull_local)
            density_full(:, :, :, ifull_local) = density_block
            overlap_full(:, :, ifull_local) = overlap_block
         END DO
         CALL cp2k_exchange_kmesh_gradient(tb%calc, tb%mol, tb%ecache, kmap%nmesh, &
                                           kmap%xkp, kmap%weight, density_full, overlap_full, &
                                           overlap_grad_full, gradient_oracle, sigma_oracle, error, &
                                           transform_backend=tb_gxtb_exchange_transform_backend())
         IF (ALLOCATED(error)) &
            CPABORT("save_tblite dense reverse oracle failed: "//TRIM(error%message))
         IF (.NOT. kmap%identity) &
            CALL tb_validate_gxtb_full_response(kpoints, kmap, ao_first, ao_size, &
                                                overlap_grad_full, &
                                                "bounded-reverse oracle overlap gradient")
         CALL tb_fold_gxtb_kpoint_matrix(kpoints, kmap, ao_first, ao_size, &
                                         overlap_grad_full, overlap_grad_oracle)
         IF (.NOT. kmap%identity) &
            CALL tb_validate_gxtb_folded_response(kpoints, kmap, ao_first, ao_size, &
                                                  overlap_grad_full, overlap_grad_oracle, &
                                                  "bounded-reverse oracle overlap gradient")
         response_scale = MAX(1.0_dp, MAXVAL(ABS(overlap_grad_oracle)))
         response_residual = MAXVAL(ABS(overlap_grad - overlap_grad_oracle))/response_scale
         force_residual = MAXVAL(ABS(gradient - gradient_oracle))
         stress_residual = MAXVAL(ABS(sigma - sigma_oracle))
         IF (para_env%is_source()) THEN
            WRITE (*, '(A,3(A,ES22.14),2(A,I0))') &
               " GXTB-QUALIFICATION_ONLY STREAMED-REVERSE", &
               " dOverlap=", response_residual, " dForce=", force_residual, &
               " dStress=", stress_residual, " peak=", peak_complex, &
               " expected=", expected_peak
         END IF
         IF (response_residual > 1.0E-10_dp .OR. force_residual > 1.0E-10_dp .OR. &
             stress_residual > 1.0E-10_dp) &
            CPABORT("QUALIFICATION_ONLY bounded reverse exceeded its 1.0E-10 dense-oracle gate")
         DEALLOCATE (density_full, gradient_oracle, overlap_full, overlap_grad_full, &
                     overlap_grad_oracle, sigma_oracle)

      END SUBROUTINE qualify_bounded_reverse

   END SUBROUTINE tb_gxtb_kpoint_exchange_gradient_streamed

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ik ...
!> \param ispin ...
!> \param fm_real ...
!> \param fm_imag ...
!> Apply the complete mixed g-xTB F(k), or add explicit exchange for a non-native backend.
!>
!> This hook is called after the ordinary real-space KS matrix has been Bloch summed, but before it
!> is redistributed to a k-point group.  AUTO/TBLITE replaces the raw matrix by the complete result
!> of save_tblite's potential mixer.  CP2K/NONE retain CP2K's matrix and only add exchange.
! **************************************************************************************************
   SUBROUTINE tb_apply_gxtb_kpoint_fm(qs_env, ik, ispin, fm_real, fm_imag)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ik, ispin
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_real
      TYPE(cp_fm_type), INTENT(INOUT), OPTIONAL          :: fm_imag

      INTEGER                                            :: iloc, jloc, ncol_global, ncol_local, &
                                                            nrow_global, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: use_mixed
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(tblite_type), POINTER                         :: tb

      NULLIFY (tb)
      IF (.NOT. ASSOCIATED(qs_env)) RETURN
      tb => qs_env%tb_tblite
      IF (.NOT. ASSOCIATED(tb)) RETURN
      use_mixed = tb%mixed_fock_kp_valid
      IF (.NOT. use_mixed .AND. .NOT. tb%exchange_fock_kp_valid .AND. &
          .NOT. tb%acp_fock_kp_valid) RETURN
      IF (use_mixed) THEN
         IF (.NOT. ALLOCATED(tb%mixed_fock_kp)) &
            CPABORT("Missing mixed g-xTB k-point Fock matrix cache")
         CPASSERT(ik >= 1 .AND. ik <= SIZE(tb%mixed_fock_kp, 4))
         CPASSERT(ispin >= 1 .AND. ispin <= SIZE(tb%mixed_fock_kp, 3))
      ELSE
         IF (tb%exchange_fock_kp_valid) THEN
            IF (.NOT. ALLOCATED(tb%exchange_fock_kp)) &
               CPABORT("Missing g-xTB k-point exchange matrix cache")
            CPASSERT(ik >= 1 .AND. ik <= SIZE(tb%exchange_fock_kp, 4))
            CPASSERT(ispin >= 1 .AND. ispin <= SIZE(tb%exchange_fock_kp, 3))
         END IF
         IF (tb%acp_fock_kp_valid) THEN
            IF (.NOT. ALLOCATED(tb%acp_fock_kp)) &
               CPABORT("Missing direct g-xTB k-point ACP matrix cache")
            CPASSERT(ik >= 1 .AND. ik <= SIZE(tb%acp_fock_kp, 3))
         END IF
      END IF

      NULLIFY (col_indices, row_indices, local_data)
      CALL cp_fm_get_info(fm_real, nrow_global=nrow_global, ncol_global=ncol_global, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, local_data=local_data)
      IF (use_mixed) THEN
         CPASSERT(nrow_global == SIZE(tb%mixed_fock_kp, 1))
         CPASSERT(ncol_global == SIZE(tb%mixed_fock_kp, 2))
      ELSE
         CPASSERT(nrow_global == tb%calc%bas%nao)
         CPASSERT(ncol_global == tb%calc%bas%nao)
      END IF
      IF (ANY(row_indices(1:nrow_local) < 1) .OR. &
          ANY(row_indices(1:nrow_local) > nrow_global)) &
         CPABORT("g-xTB exchange row index is outside the AO basis")
      IF (ANY(col_indices(1:ncol_local) < 1) .OR. &
          ANY(col_indices(1:ncol_local) > ncol_global)) &
         CPABORT("g-xTB exchange column index is outside the AO basis")
      DO jloc = 1, ncol_local
         DO iloc = 1, nrow_local
            IF (use_mixed) THEN
               local_data(iloc, jloc) = &
                  REAL(tb%mixed_fock_kp(row_indices(iloc), col_indices(jloc), ispin, ik), KIND=dp)
            ELSE
               IF (tb%acp_fock_kp_valid) &
                  local_data(iloc, jloc) = local_data(iloc, jloc) + &
                                           REAL(tb%acp_fock_kp(row_indices(iloc), col_indices(jloc), ik), KIND=dp)
               IF (tb%exchange_fock_kp_valid) &
                  local_data(iloc, jloc) = local_data(iloc, jloc) + &
                                           REAL(tb%exchange_fock_kp(row_indices(iloc), col_indices(jloc), ispin, ik), KIND=dp)
            END IF
         END DO
      END DO

      IF (PRESENT(fm_imag)) THEN
         NULLIFY (col_indices, row_indices, local_data)
         CALL cp_fm_get_info(fm_imag, nrow_global=nrow_global, ncol_global=ncol_global, &
                             nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices, local_data=local_data)
         IF (ANY(row_indices(1:nrow_local) < 1) .OR. &
             ANY(row_indices(1:nrow_local) > nrow_global)) &
            CPABORT("g-xTB exchange row index is outside the AO basis")
         IF (ANY(col_indices(1:ncol_local) < 1) .OR. &
             ANY(col_indices(1:ncol_local) > ncol_global)) &
            CPABORT("g-xTB exchange column index is outside the AO basis")
         DO jloc = 1, ncol_local
            DO iloc = 1, nrow_local
               IF (use_mixed) THEN
                  local_data(iloc, jloc) = &
                     AIMAG(tb%mixed_fock_kp(row_indices(iloc), col_indices(jloc), ispin, ik))
               ELSE
                  IF (tb%acp_fock_kp_valid) &
                     local_data(iloc, jloc) = local_data(iloc, jloc) + &
                                              AIMAG(tb%acp_fock_kp(row_indices(iloc), col_indices(jloc), ik))
                  IF (tb%exchange_fock_kp_valid) &
                     local_data(iloc, jloc) = local_data(iloc, jloc) + &
                                              AIMAG(tb%exchange_fock_kp(row_indices(iloc), col_indices(jloc), ispin, ik))
               END IF
            END DO
         END DO
      ELSEIF (use_mixed) THEN
         IF (MAXVAL(ABS(AIMAG(tb%mixed_fock_kp(:, :, ispin, ik)))) > 1.0E-10_dp) &
            CPABORT("Complex mixed g-xTB Fock matrix requires complex k-point work matrices")
      ELSE
         IF (tb%acp_fock_kp_valid) THEN
            IF (MAXVAL(ABS(AIMAG(tb%acp_fock_kp(:, :, ik)))) > 1.0E-10_dp) &
               CPABORT("Complex g-xTB ACP requires complex k-point work matrices")
         END IF
         IF (tb%exchange_fock_kp_valid) THEN
            IF (MAXVAL(ABS(AIMAG(tb%exchange_fock_kp(:, :, ispin, ik)))) > 1.0E-10_dp) &
               CPABORT("Complex g-xTB exchange requires complex k-point work matrices")
         END IF
      END IF

   END SUBROUTINE tb_apply_gxtb_kpoint_fm

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param matrix_s ...
!> \param sab_orb ...
!> \param nimg ...
!> Allocate the image-resolved multipole matrices used by the SCC potential.
! **************************************************************************************************
   SUBROUTINE tb_allocate_image_multipoles(tb, matrix_s, sab_orb, nimg)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      INTEGER, INTENT(IN)                                :: nimg

      INTEGER                                            :: i, idx, img

      IF (ASSOCIATED(tb%dipbra)) CALL dbcsr_deallocate_matrix_set(tb%dipbra)
      IF (ASSOCIATED(tb%dipket)) CALL dbcsr_deallocate_matrix_set(tb%dipket)
      IF (ASSOCIATED(tb%quadbra)) CALL dbcsr_deallocate_matrix_set(tb%quadbra)
      IF (ASSOCIATED(tb%quadket)) CALL dbcsr_deallocate_matrix_set(tb%quadket)

      CALL dbcsr_allocate_matrix_set(tb%dipbra, dip_n*nimg)
      CALL dbcsr_allocate_matrix_set(tb%dipket, dip_n*nimg)
      CALL dbcsr_allocate_matrix_set(tb%quadbra, quad_n*nimg)
      CALL dbcsr_allocate_matrix_set(tb%quadket, quad_n*nimg)
      DO img = 1, nimg
         DO i = 1, dip_n
            idx = i + dip_n*(img - 1)
            ALLOCATE (tb%dipbra(idx)%matrix, tb%dipket(idx)%matrix)
            CALL dbcsr_create(tb%dipbra(idx)%matrix, template=matrix_s(1, img)%matrix, &
                              name="DIPOLE BRA MATRIX")
            CALL dbcsr_create(tb%dipket(idx)%matrix, template=matrix_s(1, img)%matrix, &
                              name="DIPOLE KET MATRIX")
            CALL cp_dbcsr_alloc_block_from_nbl(tb%dipbra(idx)%matrix, sab_orb)
            CALL cp_dbcsr_alloc_block_from_nbl(tb%dipket(idx)%matrix, sab_orb)
         END DO
         DO i = 1, quad_n
            idx = i + quad_n*(img - 1)
            ALLOCATE (tb%quadbra(idx)%matrix, tb%quadket(idx)%matrix)
            CALL dbcsr_create(tb%quadbra(idx)%matrix, template=matrix_s(1, img)%matrix, &
                              name="QUADRUPOLE BRA MATRIX")
            CALL dbcsr_create(tb%quadket(idx)%matrix, template=matrix_s(1, img)%matrix, &
                              name="QUADRUPOLE KET MATRIX")
            CALL cp_dbcsr_alloc_block_from_nbl(tb%quadbra(idx)%matrix, sab_orb)
            CALL cp_dbcsr_alloc_block_from_nbl(tb%quadket(idx)%matrix, sab_orb)
         END DO
      END DO

   END SUBROUTINE tb_allocate_image_multipoles

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param images ...
!> \param translations ...
!> \param matrix_s ...
!> \param matrix_h ...
!> \param sab_orb ...
!> \param cell_to_index ...
!> \param nimg ...
!> Scatter save_tblite's full AO image matrices into CP2K's symmetric DBCSR representation.
!>
!> save_tblite returns A(R) with the translated AO on the row and the home-cell AO on the
!> column. CP2K stores S_mu^0,nu^R in an upper-triangular atom block. Hence entries whose
!> neighbour-list atom order is already upper triangular are transposed. The bra/ket gauges
!> of the multipoles have to be exchanged under the same transpose.
!>
!> rskp_transform subsequently uses H(k)=sum_R H(R)*exp(+2*pi*i*k.R).  There is no k-point
!> weight in H0, S, or the multipole operators; irreducible-star weights belong only to
!> density/energy contractions.  ACP is deliberately absent here and is evaluated directly on
!> the complete Bloch mesh.  Every MPI rank evaluates the same small dense save_tblite images
!> and writes only its locally owned DBCSR blocks.
! **************************************************************************************************
   SUBROUTINE tb_scatter_gxtb_images(tb, images, translations, &
                                     matrix_s, matrix_h, sab_orb, cell_to_index, nimg)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(hamiltonian_image_type), INTENT(IN)           :: images
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: translations
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_h
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, INTENT(IN)                                :: nimg

      INTEGER                                            :: i, i_first, i_nao, iatom, ic, icol, id, &
                                                            iq, irow, j_first, j_nao, jatom
      INTEGER, DIMENSION(3)                              :: cell
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block, hblock, sblock
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      CPASSERT(SIZE(images%overlap, 3) == SIZE(translations, 2))
      CPASSERT(SIZE(images%translation, 2) == SIZE(translations, 2))
      CPASSERT(SIZE(translations, 2) >= nimg)

      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, cell=cell)
         ic = cell_to_index(cell(1), cell(2), cell(3))
         CPASSERT(ic > 0 .AND. ic <= nimg)
         irow = MIN(iatom, jatom)
         icol = MAX(iatom, jatom)
         CALL tb_atom_ao_range(tb, iatom, i_first, i_nao)
         CALL tb_atom_ao_range(tb, jatom, j_first, j_nao)

         NULLIFY (sblock, hblock)
         CALL dbcsr_get_block_p(matrix_s(1, ic)%matrix, irow, icol, sblock, found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix_h(1, ic)%matrix, irow, icol, hblock, found)
         CPASSERT(found)

         IF (iatom <= jatom) THEN
            CPASSERT(SIZE(sblock, 1) == i_nao .AND. SIZE(sblock, 2) == j_nao)
            sblock = TRANSPOSE(images%overlap(j_first:j_first + j_nao - 1, &
                                              i_first:i_first + i_nao - 1, ic))
            hblock = TRANSPOSE(images%hamiltonian(j_first:j_first + j_nao - 1, &
                                                  i_first:i_first + i_nao - 1, ic))

            DO i = 1, dip_n
               id = i + dip_n*(ic - 1)
               CALL dbcsr_get_block_p(tb%dipket(id)%matrix, irow, icol, block, found)
               CPASSERT(found)
               block = TRANSPOSE(images%dipole_bra(i, j_first:j_first + j_nao - 1, &
                                                   i_first:i_first + i_nao - 1, ic))
               CALL dbcsr_get_block_p(tb%dipbra(id)%matrix, irow, icol, block, found)
               CPASSERT(found)
               block = TRANSPOSE(images%dipole_ket(i, j_first:j_first + j_nao - 1, &
                                                   i_first:i_first + i_nao - 1, ic))
            END DO
            DO i = 1, quad_n
               iq = i + quad_n*(ic - 1)
               CALL dbcsr_get_block_p(tb%quadket(iq)%matrix, irow, icol, block, found)
               CPASSERT(found)
               block = TRANSPOSE(images%quadrupole_bra(i, j_first:j_first + j_nao - 1, &
                                                       i_first:i_first + i_nao - 1, ic))
               CALL dbcsr_get_block_p(tb%quadbra(iq)%matrix, irow, icol, block, found)
               CPASSERT(found)
               block = TRANSPOSE(images%quadrupole_ket(i, j_first:j_first + j_nao - 1, &
                                                       i_first:i_first + i_nao - 1, ic))
            END DO
         ELSE
            CPASSERT(SIZE(sblock, 1) == j_nao .AND. SIZE(sblock, 2) == i_nao)
            sblock = images%overlap(j_first:j_first + j_nao - 1, &
                                    i_first:i_first + i_nao - 1, ic)
            hblock = images%hamiltonian(j_first:j_first + j_nao - 1, &
                                        i_first:i_first + i_nao - 1, ic)

            DO i = 1, dip_n
               id = i + dip_n*(ic - 1)
               CALL dbcsr_get_block_p(tb%dipbra(id)%matrix, irow, icol, block, found)
               CPASSERT(found)
               block = images%dipole_bra(i, j_first:j_first + j_nao - 1, &
                                         i_first:i_first + i_nao - 1, ic)
               CALL dbcsr_get_block_p(tb%dipket(id)%matrix, irow, icol, block, found)
               CPASSERT(found)
               block = images%dipole_ket(i, j_first:j_first + j_nao - 1, &
                                         i_first:i_first + i_nao - 1, ic)
            END DO
            DO i = 1, quad_n
               iq = i + quad_n*(ic - 1)
               CALL dbcsr_get_block_p(tb%quadbra(iq)%matrix, irow, icol, block, found)
               CPASSERT(found)
               block = images%quadrupole_bra(i, j_first:j_first + j_nao - 1, &
                                             i_first:i_first + i_nao - 1, ic)
               CALL dbcsr_get_block_p(tb%quadket(iq)%matrix, irow, icol, block, found)
               CPASSERT(found)
               block = images%quadrupole_ket(i, j_first:j_first + j_nao - 1, &
                                             i_first:i_first + i_nao - 1, ic)
            END DO
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO i = 1, SIZE(tb%dipbra)
         CALL dbcsr_finalize(tb%dipbra(i)%matrix)
         CALL dbcsr_finalize(tb%dipket(i)%matrix)
      END DO
      DO i = 1, SIZE(tb%quadbra)
         CALL dbcsr_finalize(tb%quadbra(i)%matrix)
         CALL dbcsr_finalize(tb%quadket(i)%matrix)
      END DO
   END SUBROUTINE tb_scatter_gxtb_images
#endif

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE build_tblite_matrices(qs_env, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

#if defined(__TBLITE)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_tblite_matrices'

      INTEGER                                            :: handle, maxder, nderivatives, nimg, img, nkind, i, &
                                                            ic, iw, iatom, jatom, ikind, jkind, iset, jset, n1, n2, icol, &
                                                            irow, ia, ib, sgfa, sgfb, ldsab, nseta, nsetb, &
                                                            natorb_a, natorb_b, raw_iatom, raw_jatom, &
                                                            slot, idim, jdim
      LOGICAL                                            :: do_kpoints, found, norml1, norml2, use_arnoldi, &
                                                            use_gxtb_gamma, use_gxtb_images
      REAL(KIND=dp)                                      :: dr, dshpoly, ff, hij_base, r2, rr
      REAL(KIND=dp)                                      :: native_dot_tmp, native_cn_icol, native_cn_irow, &
                                                            native_dot_weight
      INTEGER, DIMENSION(3)                              :: cell
      REAL(KIND=dp)                                      :: hij, shpoly
      REAL(KIND=dp), DIMENSION(2)                        :: condnum
      REAL(KIND=dp), DIMENSION(3)                        :: native_h0_overlap_force, native_radial_force, raw_rij, rij
      REAL(KIND=dp), DIMENSION(3, 3)                     :: native_radial_dot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(:, :, :), POINTER                :: cell_to_index
      INTEGER, DIMENSION(:, :), POINTER                   :: index_to_cell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: owork
#if defined(__TBLITE_GXTB)
      INTEGER                                            :: nimg_gxtb
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: gxtb_index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: gxtb_cell_to_index
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dense_h, projector_translations, translations
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dense_dipole, dense_quadrupole
      TYPE(error_type), ALLOCATABLE                      :: error
      TYPE(acp_image_type)                               :: acp_images
      TYPE(hamiltonian_image_type)                       :: h0_images
#endif
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: native_cn_deriv_thread
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: native_radial_force_thread
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: oint, sint, hint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: radial_hint
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min
      INTEGER, DIMENSION(:), POINTER                     :: npgfa, npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, zeta, zetb, scon_a, scon_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fblock, pblock, sblock
      TYPE(block_p_type), DIMENSION(3, 3)                 :: radial_blocks
!$    INTEGER                                            :: hash, hash1, lock_num
!$    INTEGER(KIND=int_8)                                :: iatom8
!$    INTEGER(kind=omp_lock_kind), ALLOCATABLE, DIMENSION(:) :: locks
      INTEGER, PARAMETER                                :: nlock = 501

      TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
      TYPE(atprop_type), POINTER                            :: atprop
      TYPE(cell_type), POINTER                              :: cell_obj
      TYPE(cp_blacs_env_type), POINTER                      :: blacs_env
      TYPE(cp_logger_type), POINTER                         :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER          :: matrix_dh_native, matrix_h, matrix_p, &
                                                               matrix_q_native, matrix_s, matrix_s_native, matrix_w
      TYPE(dft_control_type), POINTER                       :: dft_control
      TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
      TYPE(gto_basis_set_type), POINTER                     :: basis_set_a, basis_set_b
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER     :: basis_set_list
      TYPE(kpoint_type), POINTER                            :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_kp
      TYPE(mp_para_env_type), POINTER                       :: para_env
      TYPE(qs_energy_type), POINTER                         :: energy
      TYPE(qs_ks_env_type), POINTER                         :: ks_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER             :: qs_kind_set
      TYPE(qs_rho_type), POINTER                            :: rho
      TYPE(tblite_type), POINTER                            :: tb
      TYPE(tb_hamiltonian), POINTER                         :: h0
      TYPE(virial_type), POINTER                            :: virial

      CALL timeset(routineN, handle)

      NULLIFY (ks_env, energy, atomic_kind_set, qs_kind_set)
      NULLIFY (matrix_dh_native, matrix_h, matrix_q_native, matrix_s, matrix_s_native, atprop, dft_control)
      NULLIFY (sab_orb, sab_kp, rho, tb, kpoints, cell_to_index, index_to_cell, cell_obj)
#if defined(__TBLITE_GXTB)
      NULLIFY (gxtb_cell_to_index)
#endif

      CALL get_qs_env(qs_env=qs_env, &
                      ks_env=ks_env, para_env=para_env, &
                      energy=energy, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_h_kp=matrix_h, &
                      matrix_s_kp=matrix_s, &
                      atprop=atprop, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      sab_orb=sab_orb, &
                      sab_kp=sab_kp, &
                      rho=rho, tb_tblite=tb, cell=cell_obj)
      h0 => tb%calc%h0

      !update geometry (required for debug / geometry optimization)
      CALL tb_update_geometry(qs_env, tb)
      nkind = SIZE(atomic_kind_set)
      nderivatives = 0
      IF (calculate_forces) THEN
         nderivatives = 1
         IF (ALLOCATED(tb%grad)) DEALLOCATE (tb%grad)
         ALLOCATE (tb%grad(3, tb%mol%nat))
         IF (ALLOCATED(tb%dsedcn)) DEALLOCATE (tb%dsedcn)
         ALLOCATE (tb%dsedcn(tb%calc%bas%nsh))
         IF (ALLOCATED(tb%calc%ncoord)) THEN
            IF (ALLOCATED(tb%dcndr)) DEALLOCATE (tb%dcndr)
            ALLOCATE (tb%dcndr(3, tb%mol%nat, tb%mol%nat))
            IF (ALLOCATED(tb%dcndL)) DEALLOCATE (tb%dcndL)
            ALLOCATE (tb%dcndL(3, 3, tb%mol%nat))
         END IF
      ELSE
         IF (ALLOCATED(tb%grad)) DEALLOCATE (tb%grad)
         IF (ALLOCATED(tb%dcndr)) DEALLOCATE (tb%dcndr)
         IF (ALLOCATED(tb%dcndL)) DEALLOCATE (tb%dcndL)
      END IF
      maxder = ncoset(nderivatives)
      nimg = dft_control%nimages
      use_gxtb_gamma = .FALSE.
      use_gxtb_images = .FALSE.
#if defined(__TBLITE_GXTB)
      use_gxtb_gamma = dft_control%qs_control%xtb_control%tblite_method == gxtb .AND. .NOT. do_kpoints
      use_gxtb_images = dft_control%qs_control%xtb_control%tblite_method == gxtb .AND. do_kpoints
      IF (dft_control%qs_control%xtb_control%tblite_method == gxtb) THEN
         CALL cp2k_update_basis(tb%calc, tb%mol, tb%bcache, error, calculate_forces, tb%wfn_aux)
         IF (ALLOCATED(error)) CPABORT("save_tblite q-vSZP basis update failed")
      END IF
#endif
      IF (do_kpoints) THEN
         IF (.NOT. ASSOCIATED(sab_kp)) CPABORT("Missing k-point neighbor list for tblite")
         sab_orb => sab_kp
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)
      END IF

      !intialise hamiltonian
      CALL tb_init_ham(qs_env, tb, para_env)

      ! get density matrtix
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)

      ! set up matrices for force calculations
      IF (calculate_forces) THEN
         NULLIFY (force, matrix_w, virial)
         CALL get_qs_env(qs_env=qs_env, &
                         matrix_w_kp=matrix_w, &
                         virial=virial, force=force)

         IF (SIZE(matrix_p, 1) == 2 .AND. .NOT. (use_gxtb_gamma .OR. use_gxtb_images)) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
               CALL dbcsr_add(matrix_w(1, img)%matrix, matrix_w(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            END DO
         END IF
         tb%use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      END IF

      IF (calculate_forces .AND. .NOT. (use_gxtb_gamma .OR. use_gxtb_images)) THEN
         CALL dbcsr_allocate_matrix_set(matrix_q_native, 1, nimg)
         DO img = 1, nimg
            ALLOCATE (matrix_q_native(1, img)%matrix)
            CALL dbcsr_copy(matrix_q_native(1, img)%matrix, matrix_w(1, img)%matrix, &
                            name="TBLITE NATIVE OVERLAP FORCE MATRIX")
         END DO
      END IF

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
      IF (calculate_forces .AND. .NOT. (use_gxtb_gamma .OR. use_gxtb_images)) THEN
         ALLOCATE (native_radial_force_thread(3, SIZE(atom_of_kind)))
         ALLOCATE (native_cn_deriv_thread(SIZE(atom_of_kind)))
         native_cn_deriv_thread = 0.0_dp
         native_radial_force_thread = 0.0_dp
      END IF

      ! set up basis set lists
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ! allocate overlap matrix
      CALL dbcsr_allocate_matrix_set(matrix_s, 1, nimg)
      CALL create_sab_matrix(ks_env, matrix_s, "OVERLAP MATRIX", basis_set_list, basis_set_list, &
                             sab_orb, .TRUE.)
      CALL set_ks_env(ks_env, matrix_s_kp=matrix_s)

      ! initialize H matrix
      CALL dbcsr_allocate_matrix_set(matrix_h, 1, nimg)
      DO img = 1, nimg
         ALLOCATE (matrix_h(1, img)%matrix)
         CALL dbcsr_create(matrix_h(1, img)%matrix, template=matrix_s(1, img)%matrix, &
                           name="HAMILTONIAN MATRIX")
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_h(1, img)%matrix, sab_orb)
      END DO
      CALL set_ks_env(ks_env, matrix_h_kp=matrix_h)
#if defined(__TBLITE_GXTB)
      IF (use_gxtb_gamma) THEN
         n1 = tb%calc%bas%nao
         IF (ALLOCATED(tb%overlap_dense)) DEALLOCATE (tb%overlap_dense)
         ALLOCATE (tb%overlap_dense(n1, n1), dense_h(n1, n1), &
                   dense_dipole(3, n1, n1), dense_quadrupole(6, n1, n1))
         CALL cp2k_build_gamma_h0(tb%calc, tb%mol, tb%bcache, tb%acache, tb%accuracy, &
                                  tb%selfenergy, tb%overlap_dense, dense_dipole, &
                                  dense_quadrupole, dense_h)
         CALL tb_dense_to_dbcsr(tb, tb%overlap_dense, matrix_s(1, 1)%matrix)
         CALL tb_dense_to_dbcsr(tb, dense_h, matrix_h(1, 1)%matrix)
         DEALLOCATE (dense_h, dense_dipole, dense_quadrupole)
         GOTO 900
      ELSEIF (use_gxtb_images) THEN
         ! Preserve CP2K's compact image order as a prefix and append any missing inverse cells.
         ! Converting these exact integer cells with hmat preserves the +2*pi*k.R convention for
         ! arbitrary GENERAL k points as well as full or symmetry-reduced regular meshes; no
         ! modulo folding or grid-dependent aliasing is used here.
         CPASSERT(ASSOCIATED(index_to_cell))
         CPASSERT(ASSOCIATED(cell_obj))
         CPASSERT(SIZE(index_to_cell, 2) == nimg)
         CALL tb_close_gxtb_image_cells(index_to_cell, cell_to_index, &
                                        gxtb_index_to_cell, gxtb_cell_to_index, nimg_gxtb)
         ALLOCATE (translations(3, nimg_gxtb))
         DO img = 1, nimg_gxtb
            translations(:, img) = MATMUL(cell_obj%hmat, REAL(gxtb_index_to_cell(:, img), dp))
         END DO
         CALL cp2k_build_image_h0(tb%calc, tb%mol, tb%bcache, tb%accuracy, tb%selfenergy, &
                                  translations, h0_images)
         CALL tb_allocate_image_multipoles(tb, matrix_s, sab_orb, nimg)
         IF (ALLOCATED(tb%calc%acp)) THEN
            ! ACP is a convolution of projector/valence overlaps.  Build C(T)
            ! only on save_tblite's one-overlap cutoff, while SAB_KP above
            ! represents its larger T-V difference set.  Passing the expanded
            ! CP2K image set here would allocate an unnecessary difference of
            ! an already doubled range.
            CALL get_lattice_points(tb%mol%periodic, tb%mol%lattice, &
                                    cp2k_get_acp_cutoff(tb%calc, tb%accuracy), &
                                    projector_translations)
            CALL cp2k_build_acp_projectors(tb%calc, tb%mol, tb%bcache, tb%acache, &
                                           projector_translations, acp_images)
            CALL tb_scatter_gxtb_images(tb, h0_images, translations=translations, &
                                        matrix_s=matrix_s, matrix_h=matrix_h, sab_orb=sab_orb, &
                                        cell_to_index=cell_to_index, nimg=nimg)
            CALL tb_build_gxtb_kpoint_acp(qs_env, tb, matrix_s, sab_orb, &
                                          cell_to_index, acp_images)
            DEALLOCATE (projector_translations)
         ELSE
            CALL tb_scatter_gxtb_images(tb, h0_images, translations=translations, &
                                        matrix_s=matrix_s, matrix_h=matrix_h, sab_orb=sab_orb, &
                                        cell_to_index=cell_to_index, nimg=nimg)
         END IF
         DEALLOCATE (gxtb_cell_to_index, gxtb_index_to_cell, translations)
         GOTO 900
      END IF
#endif
      IF (calculate_forces .AND. nimg > 1) THEN
         CALL dbcsr_allocate_matrix_set(matrix_dh_native, 9, nimg)
         DO img = 1, nimg
            DO idim = 1, 9
               ALLOCATE (matrix_dh_native(idim, img)%matrix)
               CALL dbcsr_create(matrix_dh_native(idim, img)%matrix, template=matrix_s(1, img)%matrix, &
                                 name="TBLITE H0 STRAIN DERIVATIVE")
               CALL cp_dbcsr_alloc_block_from_nbl(matrix_dh_native(idim, img)%matrix, sab_orb)
            END DO
         END DO
      END IF
      ldsab = get_memory_usage(qs_kind_set, "ORB", "ORB")

      native_radial_dot = 0.0_dp

      ! loop over all atom pairs with a non-zero overlap (sab_orb)
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (calculate_forces, basis_set_list, sab_orb, matrix_dh_native, matrix_s, matrix_h, matrix_p, &
!$OMP        cell_to_index, nimg, atom_of_kind, qs_kind_set, tb, h0, ldsab, maxder, locks, &
!$OMP        native_cn_deriv_thread, native_radial_force_thread, ncoset) &
!$OMP PRIVATE (slot, hash, hash1, iatom8, lock_num, ikind, jkind, iatom, jatom, cell, rij, ic, irow, &
!$OMP        icol, dr, dshpoly, ff, hij_base, r2, raw_iatom, raw_jatom, &
!$OMP        raw_rij, native_cn_icol, native_cn_irow, native_dot_weight, native_h0_overlap_force, &
!$OMP        native_radial_force, &
!$OMP        n1, n2, ia, ib, i, idim, jdim, iset, jset, sgfa, sgfb, &
!$OMP        nseta, nsetb, natorb_a, &
!$OMP        natorb_b, found, zeta, first_sgfa, la_max, la_min, npgfa, nsgfa, rpgfa, set_radius_a, &
!$OMP        scon_a, first_sgfb, lb_max, lb_min, npgfb, nsgfb, rpgfb, set_radius_b, scon_b, zetb, rr, hij, shpoly, &
!$OMP        owork, oint, sint, hint, radial_hint, basis_set_a, basis_set_b, sblock, fblock, &
!$OMP        pblock, radial_blocks) &
!$OMP REDUCTION (+ : native_radial_dot)

!$OMP SINGLE
!$    ALLOCATE (locks(nlock))
!$OMP END SINGLE
!$OMP DO
!$    DO lock_num = 1, nlock
!$       CALL omp_init_lock(locks(lock_num))
!$    END DO
!$OMP END DO

      ALLOCATE (oint(ldsab, ldsab, maxder), owork(ldsab, ldsab))

!$OMP DO SCHEDULE(GUIDED)
      DO slot = 1, sab_orb(1)%nl_size

         ikind = sab_orb(1)%nlist_task(slot)%ikind
         jkind = sab_orb(1)%nlist_task(slot)%jkind
         iatom = sab_orb(1)%nlist_task(slot)%iatom
         jatom = sab_orb(1)%nlist_task(slot)%jatom
         cell(:) = sab_orb(1)%nlist_task(slot)%cell(:)
         rij(1:3) = sab_orb(1)%nlist_task(slot)%r(1:3)
         raw_iatom = iatom
         raw_jatom = jatom
         raw_rij = rij
         native_cn_icol = 0.0_dp
         native_cn_irow = 0.0_dp
         native_h0_overlap_force = 0.0_dp
         native_radial_force = 0.0_dp

         ! canonicalize pair ordering as in current serial flow
         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)
         IF (iatom < jatom) THEN
            rij = -rij
            i = ikind
            ikind = jkind
            jkind = i
         END IF

         dr = NORM2(rij(:))
         r2 = dr*dr

         IF (nimg == 1) THEN
            ic = 1
         ELSE
            ic = cell_to_index(cell(1), cell(2), cell(3))
            CPASSERT(ic > 0)
         END IF

         NULLIFY (sblock)
         CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
         CPASSERT(found)
         NULLIFY (fblock)
         CALL dbcsr_get_block_p(matrix=matrix_h(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=fblock, found=found)
         CPASSERT(found)
         IF (calculate_forces) THEN
            NULLIFY (pblock)
            CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                   row=irow, col=icol, BLOCK=pblock, found=found)
            CPASSERT(found)
            IF (nimg > 1) THEN
               DO jdim = 1, 3
                  DO idim = 1, 3
                     NULLIFY (radial_blocks(idim, jdim)%block)
                     CALL dbcsr_get_block_p(matrix=matrix_dh_native(idim + 3*(jdim - 1), ic)%matrix, &
                                            row=irow, col=icol, BLOCK=radial_blocks(idim, jdim)%block, &
                                            found=found)
                     CPASSERT(found)
                  END DO
               END DO
            END IF
         END IF

         ! --------- Overlap
         !get basis information
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis a
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         scon_a => basis_set_a%scon
         zeta => basis_set_a%zet
         ! basis b
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         scon_b => basis_set_b%scon
         zetb => basis_set_b%zet

         natorb_a = 0
         DO iset = 1, nseta
            natorb_a = natorb_a + (2*basis_set_a%l(1, iset) + 1)
         END DO
         natorb_b = 0
         DO iset = 1, nsetb
            natorb_b = natorb_b + (2*basis_set_b%l(1, iset) + 1)
         END DO
         ALLOCATE (sint(natorb_a, natorb_b, maxder))
         sint = 0.0_dp
         ALLOCATE (hint(natorb_a, natorb_b, maxder))
         hint = 0.0_dp
         IF (calculate_forces .AND. nimg > 1) THEN
            ALLOCATE (radial_hint(natorb_a, natorb_b, 3, 3))
            radial_hint = 0.0_dp
         END IF

         !----------------- overlap integrals
         DO iset = 1, nseta
            n1 = npgfa(iset)*(ncoset(la_max(iset)) - ncoset(la_min(iset) - 1))
            sgfa = first_sgfa(1, iset)
            DO jset = 1, nsetb
               IF (set_radius_a(iset) + set_radius_b(jset) < dr) CYCLE
               n2 = npgfb(jset)*(ncoset(lb_max(jset)) - ncoset(lb_min(jset) - 1))
               sgfb = first_sgfb(1, jset)
               IF (calculate_forces) THEN
                  CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                  lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                                  rij, sab=oint(:, :, 1), dab=oint(:, :, 2:4))
               ELSE
                  CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                  lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                                  rij, sab=oint(:, :, 1))
               END IF
               ! Contraction
               CALL contraction(oint(:, :, 1), owork, ca=scon_a(:, sgfa:), na=n1, ma=nsgfa(iset), &
                                cb=scon_b(:, sgfb:), nb=n2, mb=nsgfb(jset), fscale=1.0_dp, trans=.FALSE.)
               CALL block_add("IN", owork, nsgfa(iset), nsgfb(jset), sint(:, :, 1), sgfa, sgfb, trans=.FALSE.)
               IF (calculate_forces) THEN
                  DO i = 2, 4
                     CALL contraction(oint(:, :, i), owork, ca=scon_a(:, sgfa:), na=n1, ma=nsgfa(iset), &
                                      cb=scon_b(:, sgfb:), nb=n2, mb=nsgfb(jset), fscale=1.0_dp, trans=.FALSE.)
                     CALL block_add("IN", owork, nsgfa(iset), nsgfb(jset), sint(:, :, i), sgfa, sgfb, trans=.FALSE.)
                  END DO
               END IF
            END DO
         END DO

!$       iatom8 = INT(iatom - 1, int_8)*INT(SIZE(atom_of_kind), int_8) + INT(jatom, int_8)
!$       hash1 = INT(MOD(iatom8, INT(nlock, int_8)) + 1)
         ! update S matrix
!$       hash = hash1
!$       CALL omp_set_lock(locks(hash))
         IF (icol <= irow) THEN
            sblock(:, :) = sblock(:, :) + sint(:, :, 1)
         ELSE
            sblock(:, :) = sblock(:, :) + TRANSPOSE(sint(:, :, 1))
         END IF
!$       CALL omp_unset_lock(locks(hash))

         ! --------- Hamiltonian
         IF (icol == irow .AND. dr < same_atom) THEN
            !get diagonal F matrix from selfenergy
            n1 = tb%calc%bas%ish_at(icol)
            DO iset = 1, nseta
               sgfa = first_sgfa(1, iset)
               hij = tb%selfenergy(n1 + iset)
               DO ia = sgfa, sgfa + nsgfa(iset) - 1
                  hint(ia, ia, 1) = hij
                  IF (calculate_forces) THEN
                     native_cn_icol = native_cn_icol + tb%dsedcn(n1 + iset)*pblock(ia, ia)
                     DO jdim = 1, 3
                        DO idim = 1, 3
                           IF (nimg == 1) THEN
                              native_radial_dot(idim, jdim) = native_radial_dot(idim, jdim) + &
                                                              tb%dsedcn(n1 + iset)*tb%dcndL(idim, jdim, icol)*pblock(ia, ia)
                           ELSE
                              radial_hint(ia, ia, idim, jdim) = radial_hint(ia, ia, idim, jdim) + &
                                                                tb%dsedcn(n1 + iset)*tb%dcndL(idim, jdim, icol)
                           END IF
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         ELSE
            !get off-diagonal F matrix
            rr = SQRT(dr/(h0%rad(jkind) + h0%rad(ikind)))
            n1 = tb%calc%bas%ish_at(icol)
            DO iset = 1, nseta
               sgfa = first_sgfa(1, iset)
               n2 = tb%calc%bas%ish_at(irow)
               DO jset = 1, nsetb
                  sgfb = first_sgfb(1, jset)
                  shpoly = (1.0_dp + h0%shpoly(iset, ikind)*rr) &
                           *(1.0_dp + h0%shpoly(jset, jkind)*rr)
                  dshpoly = ((1.0_dp + h0%shpoly(iset, ikind)*rr)*h0%shpoly(jset, jkind)*rr &
                             + (1.0_dp + h0%shpoly(jset, jkind)*rr)*h0%shpoly(iset, ikind)*rr) &
                            *0.5_dp/r2
                  hij_base = 0.5_dp*(tb%selfenergy(n1 + iset) + tb%selfenergy(n2 + jset)) &
                             *h0%hscale(iset, jset, ikind, jkind)
                  hij = hij_base*shpoly
                  DO ia = sgfa, sgfa + nsgfa(iset) - 1
                     DO ib = sgfb, sgfb + nsgfb(jset) - 1
                        hint(ia, ib, 1) = hij*sint(ia, ib, 1)
                        IF (calculate_forces) THEN
                           native_dot_weight = 2.0_dp
                           IF (icol == irow) native_dot_weight = 1.0_dp
                           native_cn_icol = native_cn_icol + native_dot_weight*0.5_dp* &
                                            h0%hscale(iset, jset, ikind, jkind)*shpoly* &
                                            tb%dsedcn(n1 + iset)*pblock(ib, ia)*sint(ia, ib, 1)
                           native_cn_irow = native_cn_irow + native_dot_weight*0.5_dp* &
                                            h0%hscale(iset, jset, ikind, jkind)*shpoly* &
                                            tb%dsedcn(n2 + jset)*pblock(ib, ia)*sint(ia, ib, 1)
                           native_radial_force = native_radial_force + &
                                                 hij_base*dshpoly*pblock(ib, ia)*sint(ia, ib, 1)*raw_rij
                           native_h0_overlap_force = native_h0_overlap_force + &
                                                     hij*pblock(ib, ia)*sint(ia, ib, 2:4)
                           DO jdim = 1, 3
                              DO idim = 1, 3
                                 IF (nimg == 1) THEN
                                    native_radial_dot(idim, jdim) = native_radial_dot(idim, jdim) + &
                                                                    native_dot_weight*pblock(ib, ia)*( &
                                                                    (-hij)*sint(ia, ib, idim + 1)*rij(jdim) + &
                                                                    hij_base*dshpoly*sint(ia, ib, 1)*rij(idim)*rij(jdim) + &
                                                                0.5_dp*h0%hscale(iset, jset, ikind, jkind)*shpoly*sint(ia, ib, 1)* &
                                                                    (tb%dsedcn(n1 + iset)*tb%dcndL(idim, jdim, icol) + &
                                                                     tb%dsedcn(n2 + jset)*tb%dcndL(idim, jdim, irow)))
                                 ELSE
                                    radial_hint(ia, ib, idim, jdim) = radial_hint(ia, ib, idim, jdim) + &
                                                                      (-hij)*sint(ia, ib, idim + 1)*rij(jdim) + &
                                                                      hij_base*dshpoly*sint(ia, ib, 1)*rij(idim)*rij(jdim) + &
                                                                0.5_dp*h0%hscale(iset, jset, ikind, jkind)*shpoly*sint(ia, ib, 1)* &
                                                                      (tb%dsedcn(n1 + iset)*tb%dcndL(idim, jdim, icol) + &
                                                                       tb%dsedcn(n2 + jset)*tb%dcndL(idim, jdim, irow))
                                 END IF
                              END DO
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END IF

         ! update F matrix
!$       CALL omp_set_lock(locks(hash))
         IF (icol <= irow) THEN
            fblock(:, :) = fblock(:, :) + hint(:, :, 1)
            IF (calculate_forces .AND. nimg > 1) THEN
               DO jdim = 1, 3
                  DO idim = 1, 3
                     radial_blocks(idim, jdim)%block(:, :) = radial_blocks(idim, jdim)%block(:, :) + &
                                                             radial_hint(:, :, idim, jdim)
                  END DO
               END DO
            END IF
         ELSE
            fblock(:, :) = fblock(:, :) + TRANSPOSE(hint(:, :, 1))
            IF (calculate_forces .AND. nimg > 1) THEN
               DO jdim = 1, 3
                  DO idim = 1, 3
                     radial_blocks(idim, jdim)%block(:, :) = radial_blocks(idim, jdim)%block(:, :) + &
                                                             TRANSPOSE(radial_hint(:, :, idim, jdim))
                  END DO
               END DO
            END IF
         END IF
!$       CALL omp_unset_lock(locks(hash))

         IF (calculate_forces) THEN
            ff = 2.0_dp
!$OMP CRITICAL(tblite_native_radial_force)
            native_radial_force_thread(:, raw_iatom) = &
               native_radial_force_thread(:, raw_iatom) - ff*native_radial_force
            native_radial_force_thread(:, raw_jatom) = &
               native_radial_force_thread(:, raw_jatom) + ff*native_radial_force
            native_radial_force_thread(:, icol) = &
               native_radial_force_thread(:, icol) + ff*native_h0_overlap_force
            native_radial_force_thread(:, irow) = &
               native_radial_force_thread(:, irow) - ff*native_h0_overlap_force
            native_cn_deriv_thread(icol) = native_cn_deriv_thread(icol) + native_cn_icol
            native_cn_deriv_thread(irow) = native_cn_deriv_thread(irow) + native_cn_irow
!$OMP END CRITICAL(tblite_native_radial_force)
         END IF

         DEALLOCATE (sint, hint)
         IF (ALLOCATED(radial_hint)) DEALLOCATE (radial_hint)

      END DO
!$OMP END DO
      DEALLOCATE (oint, owork)
!$OMP DO
!$    DO lock_num = 1, nlock
!$       CALL omp_destroy_lock(locks(lock_num))
!$    END DO
!$OMP END DO
!$OMP SINGLE
!$    DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL

      IF (calculate_forces) THEN
         IF (nimg > 1) THEN
            native_radial_dot = 0.0_dp
            DO img = 1, nimg
               DO jdim = 1, 3
                  DO idim = 1, 3
                     CALL dbcsr_finalize(matrix_dh_native(idim + 3*(jdim - 1), img)%matrix)
                     CALL dbcsr_dot(matrix_dh_native(idim + 3*(jdim - 1), img)%matrix, &
                                    matrix_p(1, img)%matrix, native_dot_tmp)
                     native_radial_dot(idim, jdim) = native_radial_dot(idim, jdim) + native_dot_tmp
                  END DO
               END DO
            END DO
            CALL dbcsr_deallocate_matrix_set(matrix_dh_native)
         ELSE
            CALL para_env%sum(native_radial_dot)
         END IF
         CALL para_env%sum(native_cn_deriv_thread)
         tb%grad = 0.0_dp
         CALL tb_add_grad(tb%grad, tb%dcndr, native_cn_deriv_thread, tb%mol%nat)
         CALL tb_grad2force(qs_env, tb, para_env, 4)
         DO iatom = 1, SIZE(atom_of_kind)
            ikind = kind_of(iatom)
            force(ikind)%overlap(:, atom_of_kind(iatom)) = &
               force(ikind)%overlap(:, atom_of_kind(iatom)) + native_radial_force_thread(:, iatom)
         END DO
         IF (tb%use_virial) THEN
            virial%pv_overlap = virial%pv_overlap - native_radial_dot/para_env%num_pe
            virial%pv_virial = virial%pv_virial - native_radial_dot/para_env%num_pe
         END IF
      END IF

900   CONTINUE
      DO img = 1, nimg
         DO i = 1, SIZE(matrix_s, 1)
            CALL dbcsr_finalize(matrix_s(i, img)%matrix)
         END DO
         DO i = 1, SIZE(matrix_h, 1)
            CALL dbcsr_finalize(matrix_h(i, img)%matrix)
         END DO
      END DO

      IF (calculate_forces .AND. .NOT. (use_gxtb_gamma .OR. use_gxtb_images)) THEN
         CALL build_overlap_matrix(ks_env, matrixkp_s=matrix_s_native, &
                                   matrix_name="OVERLAP MATRIX", nderivative=nderivatives, &
                                   basis_type_a="ORB", basis_type_b="ORB", sab_nl=sab_orb, &
                                   calculate_forces=.TRUE., matrixkp_p=matrix_q_native)
         CALL dbcsr_deallocate_matrix_set(matrix_s_native)
         CALL dbcsr_deallocate_matrix_set(matrix_q_native)
      END IF

      !compute multipole moments for gfn2
      IF (ANY(dft_control%qs_control%xtb_control%tblite_method == [gfn2xtb, gxtb])) THEN
         IF (.NOT. use_gxtb_images) CALL tb_get_multipole(qs_env, tb)
      END IF

      ! output overlap information
      NULLIFY (logger)
      logger => cp_get_default_logger()
      IF (.NOT. calculate_forces) THEN
         IF (cp_print_key_should_output(logger%iter_info, qs_env%input, &
                                        "DFT%PRINT%OVERLAP_CONDITION") /= 0) THEN
            iw = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%OVERLAP_CONDITION", &
                                      extension=".Log")
            CALL section_vals_val_get(qs_env%input, "DFT%PRINT%OVERLAP_CONDITION%1-NORM", l_val=norml1)
            CALL section_vals_val_get(qs_env%input, "DFT%PRINT%OVERLAP_CONDITION%DIAGONALIZATION", l_val=norml2)
            CALL section_vals_val_get(qs_env%input, "DFT%PRINT%OVERLAP_CONDITION%ARNOLDI", l_val=use_arnoldi)
            CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env)
            CALL overlap_condnum(matrix_s, condnum, iw, norml1, norml2, use_arnoldi, blacs_env)
         END IF
      END IF

      DEALLOCATE (basis_set_list)
      IF (ALLOCATED(native_cn_deriv_thread)) DEALLOCATE (native_cn_deriv_thread)
      IF (ALLOCATED(native_radial_force_thread)) DEALLOCATE (native_radial_force_thread)

      CALL timestop(handle)

#else
      MARK_USED(qs_env)
      MARK_USED(calculate_forces)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE build_tblite_matrices

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dft_control ...
!> \param tb ...
!> \param calculate_forces ...
!> \param use_rho ...
! **************************************************************************************************
   SUBROUTINE tb_update_charges(qs_env, dft_control, tb, calculate_forces, use_rho)

      TYPE(qs_environment_type), POINTER                                 :: qs_env
      TYPE(dft_control_type), POINTER                                    :: dft_control
      TYPE(tblite_type), POINTER                                         :: tb
      LOGICAL, INTENT(IN)                                                :: calculate_forces
      LOGICAL, INTENT(IN)                                                :: use_rho

#if defined(__TBLITE)

      INTEGER                                            :: iatom, ikind, is, ns, atom_a, ii, im
      INTEGER                                            :: ispin, nspin
      INTEGER                                            :: nimg, nkind, nsgf, natorb, na, n_mix_cols, mix_offset
      INTEGER                                            :: n_atom, max_orb, max_shell
      INTEGER                                            :: raw_state_status, raw_state_unit
      LOGICAL                                            :: advance_native_mixer, discard_mixed_output, do_combined_mixing, &
                                                            do_dipole, do_kpoints, do_quadrupole, native_sign_mixing, &
                                                            skip_charge_mixing, reuse_native_input, skip_scf_dispersion, &
                                                            skip_scf_dispersion_energy, seed_native_from_rho, &
                                                            skip_scf_dispersion_gradient, skip_scf_dispersion_potential, &
                                                            use_gxtb_potential_mixer, use_native_mixer, use_no_mixer
      REAL(KIND=dp)                                      :: native_seed_charge, norm, new_charge, pao
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      INTEGER                                            :: debug_status
      CHARACTER(LEN=32)                                  :: debug_value
#endif
      CHARACTER(LEN=default_path_length)                 :: raw_state_file
      INTEGER, DIMENSION(5)                              :: occ
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(25)                             :: nao
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ch_atom, ch_shell, ch_ref, ch_orb, mix_vars
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aocg, ao_dip, ao_quad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: aocg_spin, ao_dip_spin, ao_quad_spin, &
                                                            ch_orb_spin, ch_shell_spin

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_p
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_matrix
      TYPE(dbcsr_type), POINTER                          :: s_matrix
      TYPE(error_type), ALLOCATABLE                      :: error
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      ! compute mulliken charges required for charge update
      NULLIFY (kpoints, particle_set, qs_kind_set, atomic_kind_set, scf_control)
      CALL get_qs_env(qs_env=qs_env, scf_env=scf_env, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, matrix_s_kp=matrix_s, rho=rho, para_env=para_env, &
                      scf_control=scf_control, do_kpoints=do_kpoints, kpoints=kpoints)

      ! also compute multipoles needed by GFN2
      do_dipole = .FALSE.
      do_quadrupole = .FALSE.
      skip_scf_dispersion = .FALSE.
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DEBUG_SKIP_SCF_DISPERSION", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) skip_scf_dispersion
#endif
      skip_scf_dispersion_energy = skip_scf_dispersion
      skip_scf_dispersion_gradient = skip_scf_dispersion
      skip_scf_dispersion_potential = skip_scf_dispersion
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DEBUG_SKIP_SCF_DISPERSION_ENERGY", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) skip_scf_dispersion_energy
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DEBUG_SKIP_SCF_DISPERSION_GRADIENT", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) skip_scf_dispersion_gradient
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DEBUG_SKIP_SCF_DISPERSION_POTENTIAL", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) skip_scf_dispersion_potential
#endif
      IF (dft_control%qs_control%do_ls_scf .OR. scf_control%use_ot) THEN
         use_native_mixer = .FALSE.
         use_no_mixer = .TRUE.
      ELSE
         SELECT CASE (dft_control%qs_control%xtb_control%tblite_scc_mixer)
         CASE (tblite_scc_mixer_auto)
            use_native_mixer = tb_native_scc_mixer_active(dft_control)
            use_no_mixer = .FALSE.
         CASE (tblite_scc_mixer_tblite)
            use_native_mixer = .TRUE.
            use_no_mixer = .FALSE.
         CASE (tblite_scc_mixer_cp2k)
            use_native_mixer = .FALSE.
            use_no_mixer = .FALSE.
         CASE (tblite_scc_mixer_none)
            use_native_mixer = .FALSE.
            use_no_mixer = .TRUE.
         CASE DEFAULT
            CPABORT("Unknown tblite SCC mixer")
         END SELECT
      END IF
      use_gxtb_potential_mixer = .FALSE.
#if defined(__TBLITE_GXTB)
      use_gxtb_potential_mixer = use_native_mixer .AND. &
                                 dft_control%qs_control%xtb_control%tblite_method == gxtb
#endif
      IF (use_native_mixer) THEN
         IF (.NOT. ASSOCIATED(scf_env)) CPABORT("tblite SCC mixer requires a QS SCF environment")
         IF (use_rho .AND. (.NOT. calculate_forces)) THEN
            IF (scf_env%iter_count > dft_control%qs_control%xtb_control%tblite_mixer_iterations) THEN
               CPABORT("tblite SCC mixer exceeded TBLITE_MIXER/ITERATIONS")
            END IF
            IF (scf_env%iter_count == 1 .AND. .NOT. use_gxtb_potential_mixer) THEN
               CALL tb_configure_mixer(tb, dft_control%qs_control%xtb_control%tblite_mixer_iterations, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_memory, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_damping, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_omega0, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_min_weight, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_max_weight, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_weight_factor, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_solver)
               CALL tb_reset_mixer(tb)
            END IF
         END IF
      END IF
#if defined(__TBLITE_GXTB)
      IF (dft_control%qs_control%xtb_control%tblite_method == gxtb) THEN
         IF (.NOT. use_gxtb_potential_mixer .AND. &
             .NOT. (dft_control%qs_control%do_ls_scf .OR. scf_control%use_ot) .AND. &
             scf_env%mixing_method /= direct_mixing_nr) &
            CPABORT("g-xTB with the CP2K/NONE backend currently requires DIRECT_P_MIXING")
         ! g-xTB's model-native iteration mixes the complete potential.  For AUTO/TBLITE that
         ! happens immediately before diagonalization; CP2K/NONE leave the raw SCC response here
         ! and let the ordinary CP2K SCF density/Fock controls define the external backend.
         use_native_mixer = .FALSE.
         use_no_mixer = .TRUE.
      END IF
#endif
      nspin = dft_control%nspins
      IF (nspin /= tb%wfn%nspin) CPABORT("CP2K/tblite spin channel mismatch")

      NULLIFY (matrix_p)
      IF (use_rho) THEN
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ELSEIF (calculate_forces .AND. nspin > 1) THEN
         IF (.NOT. ASSOCIATED(tb%rho_ao_kp_ref)) &
            CPABORT("Missing converged tblite density for UKS/LSD forces")
         matrix_p => tb%rho_ao_kp_ref
      ELSE
         matrix_p => scf_env%p_mix_new
      END IF
      nimg = dft_control%nimages
#if defined(__TBLITE_GXTB)
      IF (dft_control%qs_control%xtb_control%tblite_method == gxtb) THEN
         IF (.NOT. do_kpoints) THEN
            DO ispin = 1, nspin
               CALL tb_dbcsr_to_dense(tb, matrix_p(ispin, 1)%matrix, &
                                      tb%wfn%density(:, :, ispin), para_env)
            END DO
         END IF
      END IF
#endif
      IF (nspin > 1 .AND. (.NOT. calculate_forces)) CALL tb_store_density_ref(tb, matrix_p)
      IF (ASSOCIATED(tb%dipbra)) do_dipole = .TRUE.
      IF (ASSOCIATED(tb%quadbra)) do_quadrupole = .TRUE.
      reuse_native_input = .FALSE.
      IF (use_native_mixer) THEN
         IF (scf_env%iter_count == 1) THEN
            reuse_native_input = ANY(ABS(tb%wfn%qsh) > 1.0E-14_dp)
            IF (do_dipole) reuse_native_input = reuse_native_input .OR. &
                                                ANY(ABS(tb%wfn%dpat) > 1.0E-14_dp)
            IF (do_quadrupole) reuse_native_input = reuse_native_input .OR. &
                                                    ANY(ABS(tb%wfn%qpat) > 1.0E-14_dp)
         END IF
      END IF
      n_atom = SIZE(particle_set)
      nkind = SIZE(atomic_kind_set)
      CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)
      ALLOCATE (aocg(nsgf, n_atom))
      ALLOCATE (aocg_spin(nsgf, n_atom, nspin))
      aocg = 0.0_dp
      aocg_spin = 0.0_dp
      IF (do_dipole) THEN
         ALLOCATE (ao_dip(n_atom, dip_n))
         ALLOCATE (ao_dip_spin(n_atom, dip_n, nspin))
         ao_dip_spin = 0.0_dp
      END IF
      IF (do_quadrupole) THEN
         ALLOCATE (ao_quad(n_atom, quad_n))
         ALLOCATE (ao_quad_spin(n_atom, quad_n, nspin))
         ao_quad_spin = 0.0_dp
      END IF
      max_orb = 0
      max_shell = 0
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
         CALL get_xtb_atom_param(xtb_kind, natorb=natorb)
         max_orb = MAX(max_orb, natorb)
      END DO
      DO is = 1, n_atom
         max_shell = MAX(max_shell, tb%calc%bas%nsh_at(is))
      END DO
      ALLOCATE (ch_atom(n_atom, nspin), ch_shell(n_atom, max_shell))
      ALLOCATE (ch_orb(max_orb, n_atom), ch_ref(max_orb, n_atom))
      ALLOCATE (ch_orb_spin(max_orb, n_atom, nspin), ch_shell_spin(n_atom, max_shell, nspin))
      ch_atom = 0.0_dp
      ch_shell = 0.0_dp
      ch_orb = 0.0_dp
      ch_orb_spin = 0.0_dp
      ch_shell_spin = 0.0_dp
      ch_ref = 0.0_dp
      IF (do_kpoints) THEN
         DO ispin = 1, nspin
            CALL tb_ao_charges_kp_spin(matrix_p, matrix_s, aocg_spin(:, :, ispin), ispin, para_env)
            IF (do_dipole) THEN
               DO im = 1, dip_n
                  CALL tb_contract_dens_kp_spin(matrix_p, tb%dipbra, tb%dipket, im, dip_n, &
                                                ao_dip_spin(:, im, ispin), ispin, para_env)
               END DO
            END IF
            IF (do_quadrupole) THEN
               DO im = 1, quad_n
                  CALL tb_contract_dens_kp_spin(matrix_p, tb%quadbra, tb%quadket, im, quad_n, &
                                                ao_quad_spin(:, im, ispin), ispin, para_env)
               END DO
            END IF
         END DO
      ELSE
         NULLIFY (p_matrix, s_matrix)
         p_matrix => matrix_p(:, 1)
         s_matrix => matrix_s(1, 1)%matrix
         DO ispin = 1, nspin
            CALL tb_ao_charges_matrix(matrix_p(ispin, 1)%matrix, s_matrix, aocg_spin(:, :, ispin), para_env)
            IF (do_dipole) THEN
               DO im = 1, dip_n
                  CALL tb_contract_dens_matrix(matrix_p(ispin, 1)%matrix, tb%dipbra(im)%matrix, &
                                               tb%dipket(im)%matrix, ao_dip_spin(:, im, ispin), para_env)
               END DO
            END IF
            IF (do_quadrupole) THEN
               DO im = 1, quad_n
                  CALL tb_contract_dens_matrix(matrix_p(ispin, 1)%matrix, tb%quadbra(im)%matrix, &
                                               tb%quadket(im)%matrix, ao_quad_spin(:, im, ispin), para_env)
               END DO
            END IF
         END DO
      END IF
      IF (nspin == 1) THEN
         aocg(:, :) = aocg_spin(:, :, 1)
         IF (do_dipole) ao_dip(:, :) = ao_dip_spin(:, :, 1)
         IF (do_quadrupole) ao_quad(:, :) = ao_quad_spin(:, :, 1)
      ELSE
         aocg(:, :) = aocg_spin(:, :, 1) + aocg_spin(:, :, 2)
         IF (do_dipole) THEN
            DO im = 1, dip_n
               DO iatom = 1, n_atom
                  pao = ao_dip_spin(iatom, im, 1)
                  ao_dip_spin(iatom, im, 1) = pao + ao_dip_spin(iatom, im, 2)
                  ao_dip_spin(iatom, im, 2) = pao - ao_dip_spin(iatom, im, 2)
               END DO
            END DO
            ao_dip(:, :) = ao_dip_spin(:, :, 1)
         END IF
         IF (do_quadrupole) THEN
            DO im = 1, quad_n
               DO iatom = 1, n_atom
                  pao = ao_quad_spin(iatom, im, 1)
                  ao_quad_spin(iatom, im, 1) = pao + ao_quad_spin(iatom, im, 2)
                  ao_quad_spin(iatom, im, 2) = pao - ao_quad_spin(iatom, im, 2)
               END DO
            END DO
            ao_quad(:, :) = ao_quad_spin(:, :, 1)
         END IF
      END IF
      NULLIFY (xtb_kind)
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=na)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
         CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, nao=nao, occupation=occ)
         DO iatom = 1, na
            atom_a = atomic_kind_set(ikind)%atom_list(iatom)
            DO is = 1, natorb
               ns = lao(is) + 1
               norm = 2*lao(is) + 1
               ch_ref(is, atom_a) = tb%calc%h0%refocc(nao(is), ikind)/norm
               ch_orb(is, atom_a) = aocg(is, atom_a) - ch_ref(is, atom_a)
               ch_orb_spin(is, atom_a, 1) = ch_orb(is, atom_a)
               IF (nspin == 2) ch_orb_spin(is, atom_a, 2) = &
                  aocg_spin(is, atom_a, 1) - aocg_spin(is, atom_a, 2)
               ch_shell(atom_a, ns) = ch_orb(is, atom_a) + ch_shell(atom_a, ns)
               DO ispin = 1, nspin
                  ch_shell_spin(atom_a, ns, ispin) = ch_orb_spin(is, atom_a, ispin) + &
                                                     ch_shell_spin(atom_a, ns, ispin)
               END DO
            END DO
            DO ispin = 1, nspin
               ch_atom(atom_a, ispin) = SUM(ch_orb_spin(:, atom_a, ispin))
            END DO
         END DO
      END DO
      native_seed_charge = -SUM(ch_atom(:, 1))
      seed_native_from_rho = SUM(ABS(aocg_spin)) > 1.0E-10_dp .AND. &
                             ABS(native_seed_charge - REAL(dft_control%charge, dp)) < 1.0E-5_dp
      DEALLOCATE (aocg, aocg_spin)

      raw_state_status = 1
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_RAW_STATE_DUMP", raw_state_file, STATUS=raw_state_status)
#endif
      IF (raw_state_status == 0) THEN
         OPEN (NEWUNIT=raw_state_unit, FILE=TRIM(raw_state_file), STATUS="REPLACE", ACTION="WRITE")
         WRITE (raw_state_unit, *) "qat"
         DO iatom = 1, n_atom
            WRITE (raw_state_unit, "(I0,1X,ES24.16)") iatom, -ch_atom(iatom, 1)
         END DO
         WRITE (raw_state_unit, *) "qsh"
         DO iatom = 1, n_atom
            DO is = 1, tb%calc%bas%nsh_at(iatom)
               WRITE (raw_state_unit, "(I0,1X,ES24.16)") tb%calc%bas%ish_at(iatom) + is, -ch_shell(iatom, is)
            END DO
         END DO
         IF (do_dipole) THEN
            WRITE (raw_state_unit, *) "dpat"
            DO iatom = 1, n_atom
               WRITE (raw_state_unit, "(I0,3(1X,ES24.16))") iatom, -ao_dip(iatom, :)
            END DO
         END IF
         IF (do_quadrupole) THEN
            WRITE (raw_state_unit, *) "qpat"
            DO iatom = 1, n_atom
               WRITE (raw_state_unit, "(I0,6(1X,ES24.16))") iatom, -ao_quad(iatom, :)
            END DO
         END IF
         CLOSE (raw_state_unit)
      END IF

      IF (use_native_mixer) THEN
         IF (.NOT. ALLOCATED(tb%mixer)) CPABORT("tblite mixer not initialized")
         advance_native_mixer = .FALSE.
         IF (use_rho .AND. (.NOT. calculate_forces)) THEN
            advance_native_mixer = scf_env%iter_count > 1
         END IF
         IF (advance_native_mixer) THEN
            CALL tb%mixer%next(error)
            IF (ALLOCATED(error)) CPABORT("tblite native mixer failed")
            CALL tb%mixer%get(tb%wfn%qsh)
            tb%wfn%qat(:, :) = 0.0_dp
            DO iatom = 1, n_atom
               ii = tb%calc%bas%ish_at(iatom)
               DO ispin = 1, nspin
                  tb%wfn%qat(iatom, ispin) = &
                     SUM(tb%wfn%qsh(ii + 1:ii + tb%calc%bas%nsh_at(iatom), ispin))
               END DO
            END DO
            IF (do_dipole) THEN
               CALL tb%mixer%get(tb%wfn%dpat)
               DEALLOCATE (ao_dip)
            END IF
            IF (do_quadrupole) THEN
               CALL tb%mixer%get(tb%wfn%qpat)
               DEALLOCATE (ao_quad)
            END IF
         ELSE
            IF (use_rho .AND. (.NOT. calculate_forces)) THEN
               IF (.NOT. reuse_native_input) THEN
                  IF (seed_native_from_rho) THEN
                     ! Seed the native tblite mixer from CP2K's current density so SCF_GUESS/RESTART
                     ! controls the initial SCC state.
                     DO iatom = 1, n_atom
                        ii = tb%calc%bas%ish_at(iatom)
                        DO ispin = 1, nspin
                           DO is = 1, tb%calc%bas%nsh_at(iatom)
                              tb%wfn%qsh(ii + is, ispin) = -ch_shell_spin(iatom, is, ispin)
                           END DO
                           tb%wfn%qat(iatom, ispin) = &
                              SUM(tb%wfn%qsh(ii + 1:ii + tb%calc%bas%nsh_at(iatom), ispin))
                        END DO
                     END DO
                     IF (do_dipole) THEN
                        DO iatom = 1, n_atom
                           DO ispin = 1, nspin
                              tb%wfn%dpat(:, iatom, ispin) = -ao_dip_spin(iatom, :, ispin)
                           END DO
                        END DO
                     END IF
                     IF (do_quadrupole) THEN
                        DO iatom = 1, n_atom
                           DO ispin = 1, nspin
                              tb%wfn%qpat(:, iatom, ispin) = -ao_quad_spin(iatom, :, ispin)
                           END DO
                        END DO
                     END IF
                  ELSE
                     tb%wfn%qsh(:, :) = 0.0_dp
                     tb%wfn%qat(:, :) = 0.0_dp
                     IF (do_dipole) tb%wfn%dpat(:, :, :) = 0.0_dp
                     IF (do_quadrupole) tb%wfn%qpat(:, :, :) = 0.0_dp
                  END IF
               END IF
               IF (do_dipole) DEALLOCATE (ao_dip)
               IF (do_quadrupole) DEALLOCATE (ao_quad)
            ELSE
               IF ((.NOT. use_rho) .AND. (.NOT. calculate_forces)) THEN
                  CALL tb%mixer%set(tb%wfn%qsh)
                  IF (do_dipole) CALL tb%mixer%set(tb%wfn%dpat)
                  IF (do_quadrupole) CALL tb%mixer%set(tb%wfn%qpat)
               END IF
               DO iatom = 1, n_atom
                  ii = tb%calc%bas%ish_at(iatom)
                  DO ispin = 1, nspin
                     DO is = 1, tb%calc%bas%nsh_at(iatom)
                        tb%wfn%qsh(ii + is, ispin) = -ch_shell_spin(iatom, is, ispin)
                     END DO
                     tb%wfn%qat(iatom, ispin) = &
                        SUM(tb%wfn%qsh(ii + 1:ii + tb%calc%bas%nsh_at(iatom), ispin))
                  END DO
               END DO
               IF (do_dipole) THEN
                  DO iatom = 1, n_atom
                     DO ispin = 1, nspin
                        tb%wfn%dpat(:, iatom, ispin) = -ao_dip_spin(iatom, :, ispin)
                     END DO
                  END DO
                  DEALLOCATE (ao_dip)
               END IF
               IF (do_quadrupole) THEN
                  DO iatom = 1, n_atom
                     DO ispin = 1, nspin
                        tb%wfn%qpat(:, iatom, ispin) = -ao_quad_spin(iatom, :, ispin)
                     END DO
                  END DO
                  DEALLOCATE (ao_quad)
               END IF
               IF ((.NOT. use_rho) .AND. (.NOT. calculate_forces)) THEN
                  CALL tb%mixer%diff(tb%wfn%qsh)
                  IF (do_dipole) CALL tb%mixer%diff(tb%wfn%dpat)
                  IF (do_quadrupole) CALL tb%mixer%diff(tb%wfn%qpat)
               END IF
            END IF
         END IF
      ELSE
         ! charge mixing
         native_sign_mixing = .FALSE.
         IF (dft_control%qs_control%do_ls_scf .OR. scf_control%use_ot) THEN
            ! LS_SCF and OT optimize the electronic variables directly. Use the
            ! current shell/multipole response from that density without an
            ! extra SCC variable mixing step.
         ELSEIF (nspin > 1) THEN
            n_mix_cols = nspin*max_shell
            IF (do_dipole) n_mix_cols = n_mix_cols + nspin*dip_n
            IF (do_quadrupole) n_mix_cols = n_mix_cols + nspin*quad_n
            ALLOCATE (mix_vars(n_atom, n_mix_cols))
            mix_vars = 0.0_dp

            mix_offset = 0
            DO ispin = 1, nspin
               mix_vars(:, mix_offset + 1:mix_offset + max_shell) = &
                  -ch_shell_spin(:, 1:max_shell, ispin)
               mix_offset = mix_offset + max_shell
            END DO
            IF (do_dipole) THEN
               DO ispin = 1, nspin
                  mix_vars(:, mix_offset + 1:mix_offset + dip_n) = &
                     -ao_dip_spin(:, 1:dip_n, ispin)
                  mix_offset = mix_offset + dip_n
               END DO
            END IF
            IF (do_quadrupole) THEN
               DO ispin = 1, nspin
                  mix_vars(:, mix_offset + 1:mix_offset + quad_n) = &
                     -ao_quad_spin(:, 1:quad_n, ispin)
                  mix_offset = mix_offset + quad_n
               END DO
            END IF

            IF (.NOT. use_no_mixer) THEN
               CALL charge_mixing(scf_env%mixing_method, scf_env%mixing_store, &
                                  mix_vars, para_env, scf_env%iter_count)
            END IF

            mix_offset = 0
            DO ispin = 1, nspin
               ch_shell_spin(:, 1:max_shell, ispin) = &
                  -mix_vars(:, mix_offset + 1:mix_offset + max_shell)
               mix_offset = mix_offset + max_shell
            END DO
            ch_shell(:, 1:max_shell) = ch_shell_spin(:, 1:max_shell, 1)
            IF (do_dipole) THEN
               DO ispin = 1, nspin
                  ao_dip_spin(:, 1:dip_n, ispin) = &
                     -mix_vars(:, mix_offset + 1:mix_offset + dip_n)
                  mix_offset = mix_offset + dip_n
               END DO
               ao_dip(:, 1:dip_n) = ao_dip_spin(:, 1:dip_n, 1)
            END IF
            IF (do_quadrupole) THEN
               DO ispin = 1, nspin
                  ao_quad_spin(:, 1:quad_n, ispin) = &
                     -mix_vars(:, mix_offset + 1:mix_offset + quad_n)
                  mix_offset = mix_offset + quad_n
               END DO
               ao_quad(:, 1:quad_n) = ao_quad_spin(:, 1:quad_n, 1)
            END IF
            DEALLOCATE (mix_vars)
         ELSE
            do_combined_mixing = do_dipole .OR. do_quadrupole
            native_sign_mixing = do_dipole .OR. do_quadrupole
            discard_mixed_output = .FALSE.
            skip_charge_mixing = use_no_mixer
            IF (skip_charge_mixing) THEN
               !
            ELSEIF (do_combined_mixing) THEN
               n_mix_cols = max_shell
               IF (do_dipole) n_mix_cols = n_mix_cols + dip_n
               IF (do_quadrupole) n_mix_cols = n_mix_cols + quad_n
               ALLOCATE (mix_vars(n_atom, n_mix_cols))
               IF (native_sign_mixing) THEN
                  mix_vars(:, 1:max_shell) = -ch_shell(:, 1:max_shell)
               ELSE
                  mix_vars(:, 1:max_shell) = ch_shell(:, 1:max_shell)
               END IF
               mix_offset = max_shell
               IF (do_dipole) THEN
                  IF (native_sign_mixing) THEN
                     mix_vars(:, mix_offset + 1:mix_offset + dip_n) = -ao_dip(:, 1:dip_n)
                  ELSE
                     mix_vars(:, mix_offset + 1:mix_offset + dip_n) = ao_dip(:, 1:dip_n)
                  END IF
                  mix_offset = mix_offset + dip_n
               END IF
               IF (do_quadrupole) THEN
                  IF (native_sign_mixing) THEN
                     mix_vars(:, mix_offset + 1:mix_offset + quad_n) = -ao_quad(:, 1:quad_n)
                  ELSE
                     mix_vars(:, mix_offset + 1:mix_offset + quad_n) = ao_quad(:, 1:quad_n)
                  END IF
               END IF
               CALL charge_mixing(scf_env%mixing_method, scf_env%mixing_store, &
                                  mix_vars, para_env, scf_env%iter_count)
               IF (.NOT. discard_mixed_output) THEN
                  IF (native_sign_mixing) THEN
                     ch_shell(:, 1:max_shell) = -mix_vars(:, 1:max_shell)
                  ELSE
                     ch_shell(:, 1:max_shell) = mix_vars(:, 1:max_shell)
                  END IF
                  mix_offset = max_shell
                  IF (do_dipole) THEN
                     IF (native_sign_mixing) THEN
                        ao_dip(:, 1:dip_n) = -mix_vars(:, mix_offset + 1:mix_offset + dip_n)
                     ELSE
                        ao_dip(:, 1:dip_n) = mix_vars(:, mix_offset + 1:mix_offset + dip_n)
                     END IF
                     mix_offset = mix_offset + dip_n
                  END IF
                  IF (do_quadrupole) THEN
                     IF (native_sign_mixing) THEN
                        ao_quad(:, 1:quad_n) = -mix_vars(:, mix_offset + 1:mix_offset + quad_n)
                     ELSE
                        ao_quad(:, 1:quad_n) = mix_vars(:, mix_offset + 1:mix_offset + quad_n)
                     END IF
                  END IF
               END IF
               DEALLOCATE (mix_vars)
            ELSE
               CALL charge_mixing(scf_env%mixing_method, scf_env%mixing_store, &
                                  ch_shell, para_env, scf_env%iter_count)
            END IF
         END IF

         !setting new wave function
         CALL tb%pot%reset
         tb%e_es = 0.0_dp
         tb%e_scd = 0.0_dp
         IF (nspin > 1) THEN
            DO iatom = 1, n_atom
               ii = tb%calc%bas%ish_at(iatom)
               DO ispin = 1, nspin
                  DO is = 1, tb%calc%bas%nsh_at(iatom)
                     tb%wfn%qsh(ii + is, ispin) = -ch_shell_spin(iatom, is, ispin)
                  END DO
                  tb%wfn%qat(iatom, ispin) = &
                     SUM(tb%wfn%qsh(ii + 1:ii + tb%calc%bas%nsh_at(iatom), ispin))
               END DO
            END DO
            IF (do_dipole) THEN
               DO iatom = 1, n_atom
                  DO ispin = 1, nspin
                     tb%wfn%dpat(:, iatom, ispin) = -ao_dip_spin(iatom, :, ispin)
                  END DO
               END DO
               DEALLOCATE (ao_dip)
            END IF
            IF (do_quadrupole) THEN
               DO iatom = 1, n_atom
                  DO ispin = 1, nspin
                     tb%wfn%qpat(:, iatom, ispin) = -ao_quad_spin(iatom, :, ispin)
                  END DO
               END DO
               DEALLOCATE (ao_quad)
            END IF
         ELSE
            DO iatom = 1, n_atom
               ii = tb%calc%bas%ish_at(iatom)
               DO is = 1, tb%calc%bas%nsh_at(iatom)
                  new_charge = -ch_shell(iatom, is)
                  tb%wfn%qsh(ii + is, 1) = new_charge
               END DO
               IF (native_sign_mixing) THEN
                  tb%wfn%qat(iatom, 1) = SUM(tb%wfn%qsh(ii + 1:ii + tb%calc%bas%nsh_at(iatom), 1))
               ELSE
                  tb%wfn%qat(iatom, 1) = -ch_atom(iatom, 1)
               END IF
            END DO

            IF (do_dipole) THEN
               DO iatom = 1, n_atom
                  DO im = 1, dip_n
                     tb%wfn%dpat(im, iatom, 1) = -ao_dip(iatom, im)
                  END DO
               END DO
               DEALLOCATE (ao_dip)
            END IF
            IF (do_quadrupole) THEN
               DO iatom = 1, n_atom
                  DO im = 1, quad_n
                     tb%wfn%qpat(im, iatom, 1) = -ao_quad(iatom, im)
                  END DO
               END DO
               DEALLOCATE (ao_quad)
            END IF
         END IF
      END IF

      CALL tb%pot%reset
      tb%e_es = 0.0_dp
      tb%e_scd = 0.0_dp
      tb%e_int = 0.0_dp
      IF (ALLOCATED(tb%calc%coulomb)) THEN
         CALL tb%calc%coulomb%get_potential(tb%mol, tb%cache, tb%wfn, tb%pot)
         CALL tb%calc%coulomb%get_energy(tb%mol, tb%cache, tb%wfn, tb%e_es)
      END IF
      IF (ALLOCATED(tb%calc%dispersion)) THEN
         IF (.NOT. skip_scf_dispersion_potential) &
            CALL tb%calc%dispersion%get_potential(tb%mol, tb%dcache, tb%wfn, tb%pot)
         IF (.NOT. skip_scf_dispersion_energy) &
            CALL tb%calc%dispersion%get_energy(tb%mol, tb%dcache, tb%wfn, tb%e_scd)
      END IF
      IF (ALLOCATED(tb%calc%interactions)) THEN
         CALL tb%calc%interactions%get_potential(tb%mol, tb%icache, tb%wfn, tb%pot)
         CALL tb%calc%interactions%get_energy(tb%mol, tb%icache, tb%wfn, tb%e_int)
      END IF

#if defined(__TBLITE_GXTB)
      IF (dft_control%qs_control%xtb_control%tblite_method == gxtb) THEN
         tb%e_exchange = 0.0_dp
         tb%e_spin = 0.0_dp
         tb%e_acp_kp = 0.0_dp
         tb%exchange_fock_kp_valid = .FALSE.
         tb%mixed_fock_kp_valid = .FALSE.
         IF (ALLOCATED(tb%exchange_fock_kp)) DEALLOCATE (tb%exchange_fock_kp)

         IF (ALLOCATED(tb%calc%exchange)) THEN
            IF (.NOT. do_kpoints) THEN
               CALL tb%calc%exchange%update(tb%mol, tb%ecache)
               CALL tb%calc%exchange%get_potential_w_overlap(tb%mol, tb%ecache, tb%wfn, &
                                                             tb%overlap_dense, tb%pot)
               CALL tb%calc%exchange%get_energy_w_overlap(tb%mol, tb%ecache, tb%wfn, &
                                                          tb%overlap_dense, tb%e_exchange)
            ELSE
               CALL tb_build_gxtb_kpoint_exchange(qs_env, tb, matrix_s, para_env, scf_env%iter_count)
            END IF
         END IF
         IF (do_kpoints .AND. ALLOCATED(tb%calc%acp)) THEN
            CPASSERT(ASSOCIATED(kpoints))
            CALL tb_update_gxtb_kpoint_acp_energy(qs_env, tb, kpoints)
         END IF

         IF (ALLOCATED(tb%calc%repulsion)) THEN
            tb%e_rep = 0.0_dp
            CALL tb%calc%repulsion%update(tb%mol, tb%rcache)
            CALL tb%calc%repulsion%get_potential(tb%mol, tb%rcache, tb%wfn, tb%pot)
            CALL tb%calc%repulsion%get_energy(tb%mol, tb%rcache, tb%wfn, tb%e_rep)
         END IF

         IF (ALLOCATED(tb%calc%spin_polarization)) THEN
            CALL tb%calc%spin_polarization%get_potential(tb%mol, tb%scache, tb%wfn, tb%pot)
            CALL tb%calc%spin_polarization%get_energy(tb%mol, tb%scache, tb%wfn, tb%e_spin)
         END IF
      END IF
#endif

      IF (calculate_forces) THEN
         IF (ALLOCATED(tb%calc%coulomb)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%coulomb%get_gradient(tb%mol, tb%cache, tb%wfn, tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_coulomb", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 3)
         END IF

         IF (ALLOCATED(tb%calc%dispersion) .AND. .NOT. skip_scf_dispersion_gradient) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%dispersion%get_gradient(tb%mol, tb%dcache, tb%wfn, tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_dispersion_scf", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 2)
         END IF

         IF (ALLOCATED(tb%calc%interactions)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%interactions%get_gradient(tb%mol, tb%icache, tb%wfn, tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_interactions_scf", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 3)
         END IF
#if defined(__TBLITE_GXTB)
         IF (dft_control%qs_control%xtb_control%tblite_method == gxtb .AND. &
             ALLOCATED(tb%calc%repulsion)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%repulsion%get_gradient(tb%mol, tb%rcache, tb%wfn, tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_gxtb_repulsion", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 1)
         END IF
#endif
      END IF

      IF (ALLOCATED(ao_dip_spin)) DEALLOCATE (ao_dip_spin)
      IF (ALLOCATED(ao_quad_spin)) DEALLOCATE (ao_quad_spin)
      DEALLOCATE (ch_atom, ch_shell, ch_orb, ch_ref, ch_orb_spin, ch_shell_spin)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(dft_control)
      MARK_USED(calculate_forces)
      MARK_USED(use_rho)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_update_charges

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE tb_ham_add_coulomb(qs_env, tb, dft_control)

      TYPE(qs_environment_type), POINTER                       :: qs_env
      TYPE(tblite_type), POINTER                               :: tb
      TYPE(dft_control_type), POINTER                          :: dft_control

#if defined(__TBLITE)

      INTEGER                                            :: ikind, jkind, iatom, jatom, icol, irow
      INTEGER                                            :: ic, id1, id2, id3, iq1, iq2, iq3, iq4, iq5, iq6, &
                                                            is, nimg, ni, nj, i, j, nspin, first_i, first_j
      INTEGER                                            :: la, lb, za, zb
      LOGICAL                                            :: do_kpoints, found
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(25)                             :: naoa, naob
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp)                                      :: mpfac
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      INTEGER                                            :: debug_status
      CHARACTER(LEN=32)                                  :: debug_value
#endif
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of, sum_shell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ashift, bshift
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ksblock, sblock
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dip_ket1, dip_ket2, dip_ket3, &
                                                            dip_bra1, dip_bra2, dip_bra3
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: quad_ket1, quad_ket2, quad_ket3, &
                                                            quad_ket4, quad_ket5, quad_ket6, &
                                                            quad_bra1, quad_bra2, quad_bra3, &
                                                            quad_bra4, quad_bra5, quad_bra6

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: kp_list
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b

      nimg = dft_control%nimages
      mpfac = -0.5_dp

      NULLIFY (matrix_s, ks_matrix, n_list, kp_list, qs_kind_set)
      CALL get_qs_env(qs_env=qs_env, sab_orb=n_list, sab_kp=kp_list, &
                      matrix_s_kp=matrix_s, matrix_ks_kp=ks_matrix, qs_kind_set=qs_kind_set, &
                      do_kpoints=do_kpoints)
      IF (do_kpoints) THEN
         IF (.NOT. ASSOCIATED(kp_list)) CPABORT("Missing k-point neighbor list for tblite Hamiltonian")
         n_list => kp_list
      END IF
      nspin = SIZE(ks_matrix, 1)

      !creating sum of shell lists
      ALLOCATE (sum_shell(tb%mol%nat))
      i = 0
      DO j = 1, tb%mol%nat
         sum_shell(j) = i
         i = i + tb%calc%bas%nsh_at(j)
      END DO

      IF (.NOT. do_kpoints) THEN
         ! no k-points; all matrices have been transformed to periodic bsf
         CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  kind_of=kind_of)
         CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)

            ikind = kind_of(irow)
            jkind = kind_of(icol)

            ! atomic parameters
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa)
            CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob)

            ni = SIZE(sblock, 1)
            ALLOCATE (ashift(ni, ni))

            nj = SIZE(sblock, 2)
            ALLOCATE (bshift(nj, nj))

            DO is = 1, nspin
               ashift = 0.0_dp
               DO i = 1, ni
                  la = naoa(i) + sum_shell(irow)
                  ashift(i, i) = tb_spin_project(tb%pot%vsh(la, :), is)
               END DO
               bshift = 0.0_dp
               DO j = 1, nj
                  lb = naob(j) + sum_shell(icol)
                  bshift(j, j) = tb_spin_project(tb%pot%vsh(lb, :), is)
               END DO
               NULLIFY (ksblock)
               CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                      row=irow, col=icol, block=ksblock, found=found)
               CPASSERT(found)
               ksblock = ksblock - 0.5_dp*(MATMUL(ashift, sblock) &
                                           + MATMUL(sblock, bshift))
               ksblock = ksblock - 0.5_dp*(tb_spin_project(tb%pot%vat(irow, :), is) &
                                           + tb_spin_project(tb%pot%vat(icol, :), is))*sblock
            END DO
            DEALLOCATE (ashift, bshift)
         END DO
         CALL dbcsr_iterator_stop(iter)

         IF (ASSOCIATED(tb%dipbra)) THEN
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)

               NULLIFY (dip_bra1, dip_bra2, dip_bra3)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(1)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(2)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(3)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra3, found=found)
               CPASSERT(found)
               NULLIFY (dip_ket1, dip_ket2, dip_ket3)
               CALL dbcsr_get_block_p(matrix=tb%dipket(1)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipket(2)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipket(3)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket3, found=found)
               CPASSERT(found)

               DO is = 1, nspin
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock + mpfac*(dip_ket1*tb_spin_project(tb%pot%vdp(1, irow, :), is) &
                                             + dip_ket2*tb_spin_project(tb%pot%vdp(2, irow, :), is) &
                                             + dip_ket3*tb_spin_project(tb%pot%vdp(3, irow, :), is) &
                                             + dip_bra1*tb_spin_project(tb%pot%vdp(1, icol, :), is) &
                                             + dip_bra2*tb_spin_project(tb%pot%vdp(2, icol, :), is) &
                                             + dip_bra3*tb_spin_project(tb%pot%vdp(3, icol, :), is))
               END DO
            END DO
            CALL dbcsr_iterator_stop(iter)
         END IF

         IF (ASSOCIATED(tb%quadbra)) THEN
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)

               NULLIFY (quad_bra1, quad_bra2, quad_bra3, quad_bra4, quad_bra5, quad_bra6)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(1)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(2)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(3)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra3, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(4)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra4, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(5)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra5, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(6)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra6, found=found)
               CPASSERT(found)

               NULLIFY (quad_ket1, quad_ket2, quad_ket3, quad_ket4, quad_ket5, quad_ket6)
               CALL dbcsr_get_block_p(matrix=tb%quadket(1)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(2)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(3)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket3, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(4)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket4, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(5)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket5, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(6)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket6, found=found)
               CPASSERT(found)

               DO is = 1, nspin
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)

                  ksblock = ksblock + mpfac*(quad_ket1*tb_spin_project(tb%pot%vqp(1, irow, :), is) &
                                             + quad_ket2*tb_spin_project(tb%pot%vqp(2, irow, :), is) &
                                             + quad_ket3*tb_spin_project(tb%pot%vqp(3, irow, :), is) &
                                             + quad_ket4*tb_spin_project(tb%pot%vqp(4, irow, :), is) &
                                             + quad_ket5*tb_spin_project(tb%pot%vqp(5, irow, :), is) &
                                             + quad_ket6*tb_spin_project(tb%pot%vqp(6, irow, :), is) &
                                             + quad_bra1*tb_spin_project(tb%pot%vqp(1, icol, :), is) &
                                             + quad_bra2*tb_spin_project(tb%pot%vqp(2, icol, :), is) &
                                             + quad_bra3*tb_spin_project(tb%pot%vqp(3, icol, :), is) &
                                             + quad_bra4*tb_spin_project(tb%pot%vqp(4, icol, :), is) &
                                             + quad_bra5*tb_spin_project(tb%pot%vqp(5, icol, :), is) &
                                             + quad_bra6*tb_spin_project(tb%pot%vqp(6, icol, :), is))
               END DO
            END DO
            CALL dbcsr_iterator_stop(iter)
         END IF

#if defined(__TBLITE_GXTB)
         IF (dft_control%qs_control%xtb_control%tblite_method == gxtb .AND. &
             ALLOCATED(tb%calc%exchange)) THEN
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)
               CALL tb_atom_ao_range(tb, irow, first_i, ni)
               CALL tb_atom_ao_range(tb, icol, first_j, nj)
               DO is = 1, nspin
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock + tb%pot%kao(first_i:first_i + ni - 1, &
                                                 first_j:first_j + nj - 1, is)
               END DO
            END DO
            CALL dbcsr_iterator_stop(iter)
         END IF
#endif

      ELSE
         NULLIFY (kpoints)
         CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
         NULLIFY (cell_to_index)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)

         NULLIFY (nl_iterator)
         CALL neighbor_list_iterator_create(nl_iterator, n_list)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                   iatom=iatom, jatom=jatom, r=rij, cell=cellind)

            icol = MAX(iatom, jatom)
            irow = MIN(iatom, jatom)

            IF (iatom > jatom) THEN
               i = ikind
               ikind = jkind
               jkind = i
            END IF

            ic = cell_to_index(cellind(1), cellind(2), cellind(3))
            CPASSERT(ic > 0)

            NULLIFY (sblock)
            CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                   row=irow, col=icol, block=sblock, found=found)
            CPASSERT(found)

            ! atomic parameters
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa)
            CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob)

            ni = SIZE(sblock, 1)
            ALLOCATE (ashift(ni, ni))

            nj = SIZE(sblock, 2)
            ALLOCATE (bshift(nj, nj))

            DO is = 1, nspin
               ashift = 0.0_dp
               DO i = 1, ni
                  la = naoa(i) + sum_shell(irow)
                  ashift(i, i) = tb_spin_project(tb%pot%vsh(la, :), is)
               END DO
               bshift = 0.0_dp
               DO j = 1, nj
                  lb = naob(j) + sum_shell(icol)
                  bshift(j, j) = tb_spin_project(tb%pot%vsh(lb, :), is)
               END DO
               NULLIFY (ksblock)
               CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                      row=irow, col=icol, block=ksblock, found=found)
               CPASSERT(found)
               ksblock = ksblock - 0.5_dp*(MATMUL(ashift, sblock) &
                                           + MATMUL(sblock, bshift))
               ksblock = ksblock - 0.5_dp*(tb_spin_project(tb%pot%vat(irow, :), is) &
                                           + tb_spin_project(tb%pot%vat(icol, :), is))*sblock
            END DO
            DEALLOCATE (ashift, bshift)
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

         IF (ASSOCIATED(tb%dipbra)) THEN
            NULLIFY (nl_iterator)
            CALL neighbor_list_iterator_create(nl_iterator, n_list)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
               CALL get_iterator_info(nl_iterator, &
                                      iatom=iatom, jatom=jatom, cell=cellind)
               icol = MAX(iatom, jatom)
               irow = MIN(iatom, jatom)
               ic = cell_to_index(cellind(1), cellind(2), cellind(3))
               CPASSERT(ic > 0)
               id1 = 1 + dip_n*(ic - 1)
               id2 = 2 + dip_n*(ic - 1)
               id3 = 3 + dip_n*(ic - 1)

               NULLIFY (dip_bra1, dip_bra2, dip_bra3)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(id1)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(id2)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(id3)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra3, found=found)
               CPASSERT(found)
               NULLIFY (dip_ket1, dip_ket2, dip_ket3)
               CALL dbcsr_get_block_p(matrix=tb%dipket(id1)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipket(id2)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipket(id3)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket3, found=found)
               CPASSERT(found)

               DO is = 1, nspin
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock + mpfac*(dip_ket1*tb_spin_project(tb%pot%vdp(1, irow, :), is) &
                                             + dip_ket2*tb_spin_project(tb%pot%vdp(2, irow, :), is) &
                                             + dip_ket3*tb_spin_project(tb%pot%vdp(3, irow, :), is) &
                                             + dip_bra1*tb_spin_project(tb%pot%vdp(1, icol, :), is) &
                                             + dip_bra2*tb_spin_project(tb%pot%vdp(2, icol, :), is) &
                                             + dip_bra3*tb_spin_project(tb%pot%vdp(3, icol, :), is))
               END DO
            END DO
            CALL neighbor_list_iterator_release(nl_iterator)
         END IF

         IF (ASSOCIATED(tb%quadbra)) THEN
            NULLIFY (nl_iterator)
            CALL neighbor_list_iterator_create(nl_iterator, n_list)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
               CALL get_iterator_info(nl_iterator, &
                                      iatom=iatom, jatom=jatom, cell=cellind)
               icol = MAX(iatom, jatom)
               irow = MIN(iatom, jatom)
               ic = cell_to_index(cellind(1), cellind(2), cellind(3))
               CPASSERT(ic > 0)
               iq1 = 1 + quad_n*(ic - 1)
               iq2 = 2 + quad_n*(ic - 1)
               iq3 = 3 + quad_n*(ic - 1)
               iq4 = 4 + quad_n*(ic - 1)
               iq5 = 5 + quad_n*(ic - 1)
               iq6 = 6 + quad_n*(ic - 1)

               NULLIFY (quad_bra1, quad_bra2, quad_bra3, quad_bra4, quad_bra5, quad_bra6)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(iq1)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(iq2)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(iq3)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra3, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(iq4)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra4, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(iq5)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra5, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(iq6)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra6, found=found)
               CPASSERT(found)

               NULLIFY (quad_ket1, quad_ket2, quad_ket3, quad_ket4, quad_ket5, quad_ket6)
               CALL dbcsr_get_block_p(matrix=tb%quadket(iq1)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(iq2)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(iq3)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket3, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(iq4)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket4, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(iq5)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket5, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(iq6)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket6, found=found)
               CPASSERT(found)

               DO is = 1, nspin
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)

                  ksblock = ksblock + mpfac*(quad_ket1*tb_spin_project(tb%pot%vqp(1, irow, :), is) &
                                             + quad_ket2*tb_spin_project(tb%pot%vqp(2, irow, :), is) &
                                             + quad_ket3*tb_spin_project(tb%pot%vqp(3, irow, :), is) &
                                             + quad_ket4*tb_spin_project(tb%pot%vqp(4, irow, :), is) &
                                             + quad_ket5*tb_spin_project(tb%pot%vqp(5, irow, :), is) &
                                             + quad_ket6*tb_spin_project(tb%pot%vqp(6, irow, :), is) &
                                             + quad_bra1*tb_spin_project(tb%pot%vqp(1, icol, :), is) &
                                             + quad_bra2*tb_spin_project(tb%pot%vqp(2, icol, :), is) &
                                             + quad_bra3*tb_spin_project(tb%pot%vqp(3, icol, :), is) &
                                             + quad_bra4*tb_spin_project(tb%pot%vqp(4, icol, :), is) &
                                             + quad_bra5*tb_spin_project(tb%pot%vqp(5, icol, :), is) &
                                             + quad_bra6*tb_spin_project(tb%pot%vqp(6, icol, :), is))
               END DO
            END DO
            CALL neighbor_list_iterator_release(nl_iterator)
         END IF

      END IF

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(dft_control)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_ham_add_coulomb

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_get_multipole(qs_env, tb)

      TYPE(qs_environment_type), POINTER                       :: qs_env
      TYPE(tblite_type), POINTER                               :: tb

#if defined(__TBLITE)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tb_get_multipole'

      INTEGER                                     :: ikind, jkind, iatom, jatom, icol, irow, iset, jset, ityp, jtyp
      INTEGER                                     :: ic, idx, id1, id2, id3, img, iq1, iq2, iq3, iq4, iq5, iq6
      INTEGER                                     :: nkind, natom, handle, nimg, i, inda, indb
      INTEGER                                     :: atom_a, atom_b, nseta, nsetb, ia, ib, ij
      LOGICAL                                     :: do_kpoints, found
      REAL(KIND=dp)                               :: r2
      INTEGER, DIMENSION(3)                       :: cell
      INTEGER, DIMENSION(:, :, :), POINTER        :: cell_to_index
      REAL(KIND=dp), DIMENSION(3)                 :: rij
      INTEGER, DIMENSION(:), POINTER              :: la_max, lb_max
      INTEGER, DIMENSION(:), POINTER              :: nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER           :: first_sgfa, first_sgfb
      INTEGER, ALLOCATABLE                        :: atom_of_kind(:)
      REAL(KIND=dp), ALLOCATABLE                  :: stmp(:)
      REAL(KIND=dp), ALLOCATABLE                  :: dtmp(:, :), qtmp(:, :), dtmpj(:, :), qtmpj(:, :)
      REAL(KIND=dp), DIMENSION(:, :), POINTER     :: dip_ket1, dip_ket2, dip_ket3, &
                                                     dip_bra1, dip_bra2, dip_bra3
      REAL(KIND=dp), DIMENSION(:, :), POINTER     :: quad_ket1, quad_ket2, quad_ket3, &
                                                     quad_ket4, quad_ket5, quad_ket6, &
                                                     quad_bra1, quad_bra2, quad_bra3, &
                                                     quad_bra4, quad_bra5, quad_bra6

      TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                       :: dft_control
      TYPE(gto_basis_set_type), POINTER                     :: basis_set_a, basis_set_b
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER     :: basis_set_list
      TYPE(kpoint_type), POINTER                            :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_kp
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                              :: nl_iterator
      TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER             :: qs_kind_set

      CALL timeset(routineN, handle)

      !get info from environment vaiarable
      NULLIFY (atomic_kind_set, qs_kind_set, sab_orb, sab_kp, particle_set)
      NULLIFY (dft_control, matrix_s, kpoints, cell_to_index)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      sab_orb=sab_orb, &
                      sab_kp=sab_kp, &
                      particle_set=particle_set, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      kpoints=kpoints, &
                      matrix_s_kp=matrix_s)
      natom = SIZE(particle_set)
      nkind = SIZE(atomic_kind_set)
      nimg = dft_control%nimages
      IF (do_kpoints) THEN
         IF (.NOT. ASSOCIATED(sab_kp)) CPABORT("Missing k-point neighbor list for tblite multipoles")
         sab_orb => sab_kp
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind)

      !set up basis set lists
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ALLOCATE (stmp(msao(tb%calc%bas%maxl)**2))
      ALLOCATE (dtmp(dip_n, msao(tb%calc%bas%maxl)**2))
      ALLOCATE (qtmp(quad_n, msao(tb%calc%bas%maxl)**2))
      ALLOCATE (dtmpj(dip_n, msao(tb%calc%bas%maxl)**2))
      ALLOCATE (qtmpj(quad_n, msao(tb%calc%bas%maxl)**2))

      ! allocate dipole/quadrupole moment matrix elemnts
      CALL dbcsr_allocate_matrix_set(tb%dipbra, dip_n*nimg)
      CALL dbcsr_allocate_matrix_set(tb%dipket, dip_n*nimg)
      CALL dbcsr_allocate_matrix_set(tb%quadbra, quad_n*nimg)
      CALL dbcsr_allocate_matrix_set(tb%quadket, quad_n*nimg)
      DO img = 1, nimg
         DO i = 1, dip_n
            idx = i + dip_n*(img - 1)
            ALLOCATE (tb%dipbra(idx)%matrix)
            ALLOCATE (tb%dipket(idx)%matrix)
            CALL dbcsr_create(tb%dipbra(idx)%matrix, template=matrix_s(1, img)%matrix, &
                              name="DIPOLE BRAMATRIX")
            CALL dbcsr_create(tb%dipket(idx)%matrix, template=matrix_s(1, img)%matrix, &
                              name="DIPOLE KETMATRIX")
            CALL cp_dbcsr_alloc_block_from_nbl(tb%dipbra(idx)%matrix, sab_orb)
            CALL cp_dbcsr_alloc_block_from_nbl(tb%dipket(idx)%matrix, sab_orb)
         END DO
         DO i = 1, quad_n
            idx = i + quad_n*(img - 1)
            ALLOCATE (tb%quadbra(idx)%matrix)
            ALLOCATE (tb%quadket(idx)%matrix)
            CALL dbcsr_create(tb%quadbra(idx)%matrix, template=matrix_s(1, img)%matrix, &
                              name="QUADRUPOLE BRAMATRIX")
            CALL dbcsr_create(tb%quadket(idx)%matrix, template=matrix_s(1, img)%matrix, &
                              name="QUADRUPOLE KETMATRIX")
            CALL cp_dbcsr_alloc_block_from_nbl(tb%quadbra(idx)%matrix, sab_orb)
            CALL cp_dbcsr_alloc_block_from_nbl(tb%quadket(idx)%matrix, sab_orb)
         END DO
      END DO

      !loop over all atom pairs with a non-zero overlap (sab_orb)
      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)

         r2 = NORM2(rij(:))**2

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)

         IF (iatom < jatom) THEN
            rij = -rij
            i = ikind
            ikind = jkind
            jkind = i
         END IF

         IF (nimg == 1) THEN
            ic = 1
         ELSE
            ic = cell_to_index(cell(1), cell(2), cell(3))
            CPASSERT(ic > 0)
         END IF
         id1 = 1 + dip_n*(ic - 1)
         id2 = 2 + dip_n*(ic - 1)
         id3 = 3 + dip_n*(ic - 1)
         iq1 = 1 + quad_n*(ic - 1)
         iq2 = 2 + quad_n*(ic - 1)
         iq3 = 3 + quad_n*(ic - 1)
         iq4 = 4 + quad_n*(ic - 1)
         iq5 = 5 + quad_n*(ic - 1)
         iq6 = 6 + quad_n*(ic - 1)

         ityp = tb%mol%id(icol)
         jtyp = tb%mol%id(irow)

         NULLIFY (dip_bra1, dip_bra2, dip_bra3)
         CALL dbcsr_get_block_p(matrix=tb%dipbra(id1)%matrix, &
                                row=irow, col=icol, BLOCK=dip_bra1, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%dipbra(id2)%matrix, &
                                row=irow, col=icol, BLOCK=dip_bra2, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%dipbra(id3)%matrix, &
                                row=irow, col=icol, BLOCK=dip_bra3, found=found)
         CPASSERT(found)

         NULLIFY (dip_ket1, dip_ket2, dip_ket3)
         CALL dbcsr_get_block_p(matrix=tb%dipket(id1)%matrix, &
                                row=irow, col=icol, BLOCK=dip_ket1, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%dipket(id2)%matrix, &
                                row=irow, col=icol, BLOCK=dip_ket2, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%dipket(id3)%matrix, &
                                row=irow, col=icol, BLOCK=dip_ket3, found=found)
         CPASSERT(found)

         NULLIFY (quad_bra1, quad_bra2, quad_bra3, quad_bra4, quad_bra5, quad_bra6)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(iq1)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra1, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(iq2)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra2, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(iq3)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra3, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(iq4)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra4, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(iq5)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra5, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(iq6)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra6, found=found)
         CPASSERT(found)

         NULLIFY (quad_ket1, quad_ket2, quad_ket3, quad_ket4, quad_ket5, quad_ket6)
         CALL dbcsr_get_block_p(matrix=tb%quadket(iq1)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket1, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(iq2)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket2, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(iq3)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket3, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(iq4)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket4, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(iq5)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket5, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(iq6)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket6, found=found)
         CPASSERT(found)

         !get basis information
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         atom_a = atom_of_kind(icol)
         atom_b = atom_of_kind(irow)
         ! basis a
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         ! basis b
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set

         ! --------- Hamiltonian
         ! Periodic self-images are off-diagonal lattice contributions, not the on-site block.
         IF (icol == irow .AND. r2 < same_atom**2) THEN
            DO iset = 1, nseta
               DO jset = 1, nsetb
#if defined(__TBLITE_GXTB)
                  CALL cp2k_multipole_cgto(tb%calc%bas, tb%bcache, jset, ityp, irow, &
                       & iset, ityp, icol, r2, -rij, stmp, dtmp, qtmp)
#else
                  CALL multipole_cgto(tb%calc%bas%cgto(jset, ityp), tb%calc%bas%cgto(iset, ityp), &
                       & r2, -rij, tb%calc%bas%intcut, stmp, dtmp, qtmp)
#endif

                  DO inda = 1, nsgfa(iset)
                     ia = first_sgfa(1, iset) - first_sgfa(1, 1) + inda
                     DO indb = 1, nsgfb(jset)
                        ib = first_sgfb(1, jset) - first_sgfb(1, 1) + indb
                        ij = indb + nsgfb(jset)*(inda - 1)

                        dip_ket1(ib, ia) = dip_ket1(ib, ia) + dtmp(1, ij)
                        dip_ket2(ib, ia) = dip_ket2(ib, ia) + dtmp(2, ij)
                        dip_ket3(ib, ia) = dip_ket3(ib, ia) + dtmp(3, ij)

                        quad_ket1(ib, ia) = quad_ket1(ib, ia) + qtmp(1, ij)
                        quad_ket2(ib, ia) = quad_ket2(ib, ia) + qtmp(2, ij)
                        quad_ket3(ib, ia) = quad_ket3(ib, ia) + qtmp(3, ij)
                        quad_ket4(ib, ia) = quad_ket4(ib, ia) + qtmp(4, ij)
                        quad_ket5(ib, ia) = quad_ket5(ib, ia) + qtmp(5, ij)
                        quad_ket6(ib, ia) = quad_ket6(ib, ia) + qtmp(6, ij)

                        dip_bra1(ib, ia) = dip_bra1(ib, ia) + dtmp(1, ij)
                        dip_bra2(ib, ia) = dip_bra2(ib, ia) + dtmp(2, ij)
                        dip_bra3(ib, ia) = dip_bra3(ib, ia) + dtmp(3, ij)

                        quad_bra1(ib, ia) = quad_bra1(ib, ia) + qtmp(1, ij)
                        quad_bra2(ib, ia) = quad_bra2(ib, ia) + qtmp(2, ij)
                        quad_bra3(ib, ia) = quad_bra3(ib, ia) + qtmp(3, ij)
                        quad_bra4(ib, ia) = quad_bra4(ib, ia) + qtmp(4, ij)
                        quad_bra5(ib, ia) = quad_bra5(ib, ia) + qtmp(5, ij)
                        quad_bra6(ib, ia) = quad_bra6(ib, ia) + qtmp(6, ij)
                     END DO
                  END DO
               END DO
            END DO
         ELSE
            DO iset = 1, nseta
               DO jset = 1, nsetb
#if defined(__TBLITE_GXTB)
                  CALL cp2k_multipole_cgto(tb%calc%bas, tb%bcache, jset, jtyp, irow, &
                      & iset, ityp, icol, r2, -rij, stmp, dtmp, qtmp)
#else
                  CALL multipole_cgto(tb%calc%bas%cgto(jset, jtyp), tb%calc%bas%cgto(iset, ityp), &
                      & r2, -rij, tb%calc%bas%intcut, stmp, dtmp, qtmp)
#endif

                  DO inda = 1, nsgfa(iset)
                     ia = first_sgfa(1, iset) - first_sgfa(1, 1) + inda
                     DO indb = 1, nsgfb(jset)
                        ib = first_sgfb(1, jset) - first_sgfb(1, 1) + indb

                        ij = indb + nsgfb(jset)*(inda - 1)
                        CALL tb_shift_multipole(-rij, stmp(ij), dtmp(:, ij), qtmp(:, ij), &
                                                dtmpj(:, ij), qtmpj(:, ij))

                        dip_bra1(ib, ia) = dip_bra1(ib, ia) + dtmp(1, ij)
                        dip_bra2(ib, ia) = dip_bra2(ib, ia) + dtmp(2, ij)
                        dip_bra3(ib, ia) = dip_bra3(ib, ia) + dtmp(3, ij)

                        quad_bra1(ib, ia) = quad_bra1(ib, ia) + qtmp(1, ij)
                        quad_bra2(ib, ia) = quad_bra2(ib, ia) + qtmp(2, ij)
                        quad_bra3(ib, ia) = quad_bra3(ib, ia) + qtmp(3, ij)
                        quad_bra4(ib, ia) = quad_bra4(ib, ia) + qtmp(4, ij)
                        quad_bra5(ib, ia) = quad_bra5(ib, ia) + qtmp(5, ij)
                        quad_bra6(ib, ia) = quad_bra6(ib, ia) + qtmp(6, ij)

                        dip_ket1(ib, ia) = dip_ket1(ib, ia) + dtmpj(1, ij)
                        dip_ket2(ib, ia) = dip_ket2(ib, ia) + dtmpj(2, ij)
                        dip_ket3(ib, ia) = dip_ket3(ib, ia) + dtmpj(3, ij)

                        quad_ket1(ib, ia) = quad_ket1(ib, ia) + qtmpj(1, ij)
                        quad_ket2(ib, ia) = quad_ket2(ib, ia) + qtmpj(2, ij)
                        quad_ket3(ib, ia) = quad_ket3(ib, ia) + qtmpj(3, ij)
                        quad_ket4(ib, ia) = quad_ket4(ib, ia) + qtmpj(4, ij)
                        quad_ket5(ib, ia) = quad_ket5(ib, ia) + qtmpj(5, ij)
                        quad_ket6(ib, ia) = quad_ket6(ib, ia) + qtmpj(6, ij)
                     END DO
                  END DO
               END DO
            END DO
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO i = 1, SIZE(tb%dipbra)
         CALL dbcsr_finalize(tb%dipbra(i)%matrix)
         CALL dbcsr_finalize(tb%dipket(i)%matrix)
      END DO
      DO i = 1, SIZE(tb%quadbra)
         CALL dbcsr_finalize(tb%quadbra(i)%matrix)
         CALL dbcsr_finalize(tb%quadket(i)%matrix)
      END DO

      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_get_multipole

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> Complete save_tblite Gamma-point response for g-xTB.
! **************************************************************************************************
   SUBROUTINE tb_gxtb_gamma_gradient(qs_env, tb)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb

#if defined(__TBLITE_GXTB)
      INTEGER                                            :: ispin, nao, nspin
      REAL(KIND=dp), ALLOCATABLE                         :: dense_w(:, :, :), gamma_sigma(:, :)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_w
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (dft_control, matrix_w, para_env)
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, matrix_w_kp=matrix_w, &
                      para_env=para_env)
      IF (dft_control%nimages /= 1) &
         CPABORT("g-xTB Gamma gradient received k-point image matrices")

      nao = tb%calc%bas%nao
      nspin = tb%wfn%nspin
      ALLOCATE (dense_w(nao, nao, nspin), gamma_sigma(3, 3))
      DO ispin = 1, nspin
         CALL tb_dbcsr_to_dense(tb, matrix_w(ispin, 1)%matrix, dense_w(:, :, ispin), para_env)
      END DO
      tb%grad = 0.0_dp
      CALL cp2k_gamma_gradient(tb%calc, tb%mol, tb%bcache, tb%acache, tb%ecache, &
                               tb%wfn_aux, tb%accuracy, tb%selfenergy, tb%dsedcn, &
                               tb%overlap_dense, tb%pot, tb%wfn, dense_w, tb%grad, gamma_sigma)
      CALL tb_grad2force(qs_env, tb, para_env, 4)
      tb%sigma = tb%sigma + gamma_sigma
      CALL tb_dump_sigma_component("after_gxtb_gamma_response", tb%sigma, para_env)
      IF (tb%use_virial) CALL tb_add_stress(qs_env, tb, para_env)

      DEALLOCATE (dense_w, gamma_sigma)
#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      CPABORT("Built without save_tblite g-xTB support")
#endif

   END SUBROUTINE tb_gxtb_gamma_gradient

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> Complete the image-resolved, non-exchange g-xTB response for CP2K k points.
!>
!> CP2K's real-space P(R) and W(R) already contain the irreducible-star weights and symmetry
!> expansion.  They are gathered without symmetrizing individual images, then passed in their
!> original up/down spin representation to save_tblite.  The explicit nonlinear exchange response
!> first augments matrix_w with its overlap derivative; its direct atomic and stress derivatives
!> are then accumulated with the image-resolved H0/Pulay and ACP response.
! **************************************************************************************************
   SUBROUTINE tb_gxtb_kpoint_gradient(qs_env, tb)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb

#if defined(__TBLITE_GXTB)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: acp_density_k(:, :, :, :)
      INTEGER                                            :: img, nao, nimg, nimg_response, nspin
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index, response_cell_to_index
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_first, ao_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: response_index_to_cell
      LOGICAL                                            :: validate_density
      REAL(KIND=dp), ALLOCATABLE                         :: acp_gradient(:, :), acp_sigma(:, :), &
                                                            dEdcnbas(:), dEdqbas(:), &
                                                            density_images(:, :, :, :), &
                                                            exchange_gradient(:, :), &
                                                            exchange_sigma(:, :), &
                                                            gradient(:, :), overlap_gamma(:, :), &
                                                            overlap_images(:, :, :, :), &
                                                            projector_translations(:, :), &
                                                            sigma(:, :), translations(:, :), &
                                                            wdensity_images(:, :, :, :)
      TYPE(acp_image_type)                               :: acp_images
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s, matrix_w
      TYPE(error_type), ALLOCATABLE                      :: error
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_kp
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(tb_gxtb_kmesh_map_type), POINTER              :: kmap

      NULLIFY (cell, cell_to_index, index_to_cell, kmap, kpoints, matrix_p, matrix_s, matrix_w, &
               para_env, response_cell_to_index, rho, sab_kp)
      CALL get_qs_env(qs_env=qs_env, cell=cell, kpoints=kpoints, matrix_s_kp=matrix_s, &
                      matrix_w_kp=matrix_w, para_env=para_env, rho=rho, sab_kp=sab_kp)
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index, index_to_cell=index_to_cell)

      IF (.NOT. ASSOCIATED(cell) .OR. .NOT. ASSOCIATED(sab_kp) .OR. &
          .NOT. ASSOCIATED(cell_to_index) .OR. .NOT. ASSOCIATED(index_to_cell)) &
         CPABORT("Incomplete CP2K k-point image data for the g-xTB gradient")

      nao = tb%calc%bas%nao
      nspin = SIZE(matrix_p, 1)
      nimg = SIZE(matrix_p, 2)
      CPASSERT(nspin == tb%wfn%nspin)
      CPASSERT(SIZE(matrix_s, 1) == 1 .AND. SIZE(matrix_s, 2) == nimg)
      CPASSERT(SIZE(matrix_w, 1) == nspin .AND. SIZE(matrix_w, 2) == nimg)
      CPASSERT(SIZE(index_to_cell, 2) == nimg)

      CALL tb_close_gxtb_image_cells(index_to_cell, cell_to_index, &
                                     response_index_to_cell, response_cell_to_index, nimg_response)

      ALLOCATE (density_images(nao, nao, nimg_response, nspin), &
                wdensity_images(nao, nao, nimg_response, nspin), &
                overlap_images(nao, nao, nimg_response, 1), overlap_gamma(nao, nao), &
                translations(3, nimg_response), gradient(3, tb%mol%nat), sigma(3, 3))
      DO img = 1, nimg_response
         translations(:, img) = MATMUL(cell%hmat, REAL(response_index_to_cell(:, img), dp))
      END DO

      CALL tb_gather_gxtb_images(tb, matrix_p, sab_kp, cell_to_index, &
                                 response_cell_to_index, density_images, para_env)
      CALL tb_gather_gxtb_images(tb, matrix_s, sab_kp, cell_to_index, &
                                 response_cell_to_index, overlap_images, para_env)
      overlap_gamma(:, :) = SUM(overlap_images(:, :, :, 1), DIM=3)

      IF (ALLOCATED(tb%calc%exchange)) THEN
         ALLOCATE (exchange_gradient(3, tb%mol%nat), exchange_sigma(3, 3))
         CALL tb_gxtb_kpoint_exchange_gradient(qs_env, tb, exchange_gradient, exchange_sigma)
      END IF
      CALL tb_gather_gxtb_images(tb, matrix_w, sab_kp, cell_to_index, &
                                 response_cell_to_index, wdensity_images, para_env)

      CALL cp2k_image_h0_gradient(tb%calc, tb%mol, tb%bcache, tb%wfn_aux, tb%accuracy, &
                                  translations, tb%selfenergy, tb%dsedcn, overlap_gamma, &
                                  tb%pot, density_images, wdensity_images, gradient, sigma)

      IF (ALLOCATED(exchange_gradient)) THEN
         gradient(:, :) = gradient(:, :) + exchange_gradient(:, :)
         sigma(:, :) = sigma(:, :) + exchange_sigma(:, :)
         DEALLOCATE (exchange_gradient, exchange_sigma)
      END IF

      IF (ALLOCATED(tb%calc%acp)) THEN
         CALL get_lattice_points(tb%mol%periodic, tb%mol%lattice, &
                                 cp2k_get_acp_cutoff(tb%calc, tb%accuracy), &
                                 projector_translations)
         CALL cp2k_build_acp_projectors(tb%calc, tb%mol, tb%bcache, tb%acache, &
                                        projector_translations, acp_images)
         CALL tb_prepare_gxtb_kmesh_map(tb, kpoints, matrix_s, sab_kp, cell_to_index)
         kmap => tb%gxtb_kmesh_map
         CALL tb_gxtb_ao_layout(tb, ao_first, ao_size)
         ALLOCATE (acp_density_k(nao, nao, nspin, kmap%nfull), &
                   dEdcnbas(tb%mol%nat), dEdqbas(tb%mol%nat), &
                   acp_gradient(3, tb%mol%nat), acp_sigma(3, 3))
         IF (kmap%identity) THEN
            CALL tb_collect_kpoint_density(qs_env, tb, kpoints, acp_density_k)
         ELSE
            CALL tb_collect_kpoint_density(qs_env, tb, kpoints, acp_density_k, &
                                           kmap%representative)
            validate_density = .NOT. tb%gxtb_density_map_validated .AND. &
                               MAXVAL(ABS(acp_density_k)) > 1.0E-14_dp
            CALL tb_expand_gxtb_kpoint_density(kpoints, kmap, ao_first, ao_size, &
                                               acp_density_k, validate_density)
            IF (validate_density) tb%gxtb_density_map_validated = .TRUE.
         END IF
         IF (.NOT. ALL(ieee_is_finite(REAL(acp_density_k, KIND=dp))) .OR. &
             .NOT. ALL(ieee_is_finite(AIMAG(acp_density_k)))) &
            CPABORT("CP2K supplied a non-finite k-point density to the ACP gradient")
         dEdcnbas = 0.0_dp
         dEdqbas = 0.0_dp
         acp_gradient = 0.0_dp
         acp_sigma = 0.0_dp
         CALL cp2k_acp_kmesh_gradient(tb%calc, tb%mol, tb%bcache, tb%acache, &
                                      acp_images, kmap%xkp, kmap%weight, acp_density_k, &
                                      dEdcnbas, dEdqbas, acp_gradient, acp_sigma, error)
         IF (ALLOCATED(error)) &
            CPABORT("save_tblite g-xTB ACP k-mesh gradient failed: "//TRIM(error%message))
         CALL cp2k_apply_basis_gradient(tb%calc, tb%mol, tb%wfn_aux, dEdcnbas, dEdqbas, &
                                        acp_gradient, acp_sigma)
         gradient(:, :) = gradient(:, :) + acp_gradient(:, :)
         sigma(:, :) = sigma(:, :) + acp_sigma(:, :)
         DEALLOCATE (acp_density_k, acp_gradient, acp_sigma, ao_first, ao_size, &
                     dEdcnbas, dEdqbas, projector_translations)
      END IF

      tb%grad(:, :) = gradient(:, :)
      CALL tb_grad2force(qs_env, tb, para_env, 4)
      tb%sigma = tb%sigma + sigma
      CALL tb_dump_sigma_component("after_gxtb_kpoint_response", tb%sigma, para_env)
      IF (tb%use_virial) CALL tb_add_stress(qs_env, tb, para_env)

      DEALLOCATE (density_images, gradient, overlap_gamma, overlap_images, sigma, &
                  translations, wdensity_images, response_cell_to_index, response_index_to_cell)
#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      CPABORT("Built without save_tblite g-xTB support")
#endif

   END SUBROUTINE tb_gxtb_kpoint_gradient

! **************************************************************************************************
!> \brief Shift a multipole operator from one center to the other.
!> \param vec displacement vector between the two centers
!> \param s overlap integral
!> \param di dipole integral on the original center
!> \param qi quadrupole integral on the original center
!> \param dj dipole integral on the shifted center
!> \param qj quadrupole integral on the shifted center
! **************************************************************************************************
   PURE SUBROUTINE tb_shift_multipole(vec, s, di, qi, dj, qj)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vec
      REAL(KIND=dp), INTENT(IN)                          :: s
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: di, qi
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: dj, qj

      REAL(KIND=dp)                                      :: tr

      dj(1) = di(1) + vec(1)*s
      dj(2) = di(2) + vec(2)*s
      dj(3) = di(3) + vec(3)*s

      qj(1) = 2*vec(1)*di(1) + vec(1)**2*s
      qj(3) = 2*vec(2)*di(2) + vec(2)**2*s
      qj(6) = 2*vec(3)*di(3) + vec(3)**2*s
      qj(2) = vec(1)*di(2) + vec(2)*di(1) + vec(1)*vec(2)*s
      qj(4) = vec(1)*di(3) + vec(3)*di(1) + vec(1)*vec(3)*s
      qj(5) = vec(2)*di(3) + vec(3)*di(2) + vec(2)*vec(3)*s
      tr = 0.5_dp*(qj(1) + qj(3) + qj(6))

      qj(1) = qi(1) + 1.5_dp*qj(1) - tr
      qj(2) = qi(2) + 1.5_dp*qj(2)
      qj(3) = qi(3) + 1.5_dp*qj(3) - tr
      qj(4) = qi(4) + 1.5_dp*qj(4)
      qj(5) = qi(5) + 1.5_dp*qj(5)
      qj(6) = qi(6) + 1.5_dp*qj(6) - tr

   END SUBROUTINE tb_shift_multipole

! **************************************************************************************************
!> \brief compute the mulliken properties (AO resolved)
!> \param p_mat ...
!> \param s_matrix ...
!> \param charges ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_ao_charges_matrix(p_mat, s_matrix, charges, para_env)
      TYPE(dbcsr_type), POINTER                          :: p_mat, s_matrix
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: charges
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: i, iblock_col, iblock_row, j
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_block, s_block
      TYPE(dbcsr_iterator_type)                          :: iter

      charges = 0.0_dp
      CALL dbcsr_iterator_start(iter, s_matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         NULLIFY (s_block, p_block)
         CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, s_block)
         CALL dbcsr_get_block_p(matrix=p_mat, row=iblock_row, col=iblock_col, BLOCK=p_block, found=found)
         IF (.NOT. found) CYCLE
         IF (.NOT. (ASSOCIATED(s_block) .AND. ASSOCIATED(p_block))) CYCLE

         DO j = 1, SIZE(p_block, 2)
            DO i = 1, SIZE(p_block, 1)
               charges(i, iblock_row) = charges(i, iblock_row) + p_block(i, j)*s_block(i, j)
            END DO
         END DO
         IF (iblock_col /= iblock_row) THEN
            DO j = 1, SIZE(p_block, 2)
               DO i = 1, SIZE(p_block, 1)
                  charges(j, iblock_col) = charges(j, iblock_col) + p_block(i, j)*s_block(i, j)
               END DO
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL para_env%sum(charges)

   END SUBROUTINE tb_ao_charges_matrix

! **************************************************************************************************
!> \brief compute the AO-resolved Mulliken charges for one k-point spin channel.
!> \param p_matrix_kp ...
!> \param s_matrix_kp ...
!> \param charges ...
!> \param ispin ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_ao_charges_kp_spin(p_matrix_kp, s_matrix_kp, charges, ispin, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: p_matrix_kp, s_matrix_kp
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: charges
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: ic
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: image_charges
      TYPE(dbcsr_type), POINTER                          :: p_mat, s_mat

      charges = 0.0_dp
      ALLOCATE (image_charges(SIZE(charges, 1), SIZE(charges, 2)))
      DO ic = 1, SIZE(s_matrix_kp, 2)
         NULLIFY (p_mat, s_mat)
         p_mat => p_matrix_kp(ispin, ic)%matrix
         s_mat => s_matrix_kp(1, ic)%matrix
         IF (ASSOCIATED(p_mat) .AND. ASSOCIATED(s_mat)) THEN
            image_charges = 0.0_dp
            CALL tb_ao_charges_matrix(p_mat, s_mat, image_charges, para_env)
            charges(:, :) = charges(:, :) + image_charges(:, :)
         END IF
      END DO
      DEALLOCATE (image_charges)

   END SUBROUTINE tb_ao_charges_kp_spin

! **************************************************************************************************
!> \brief compute the mulliken properties (AO resolved)
!> \param p_mat ...
!> \param bra_mat ...
!> \param ket_mat ...
!> \param output ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_contract_dens_matrix(p_mat, bra_mat, ket_mat, output, para_env)
      TYPE(dbcsr_type), POINTER                          :: p_mat, bra_mat, ket_mat
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: output
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: i, iblock_col, iblock_row, j
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: bra, ket, p_block
      TYPE(dbcsr_iterator_type)                          :: iter

      output = 0.0_dp
      CALL dbcsr_iterator_start(iter, bra_mat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         NULLIFY (p_block, bra, ket)
         CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, bra)
         CALL dbcsr_get_block_p(matrix=p_mat, row=iblock_row, col=iblock_col, BLOCK=p_block, found=found)
         IF (.NOT. found) CYCLE
         CALL dbcsr_get_block_p(matrix=ket_mat, row=iblock_row, col=iblock_col, BLOCK=ket, found=found)
         IF (.NOT. found) CPABORT("missing block")

         IF (.NOT. (ASSOCIATED(bra) .AND. ASSOCIATED(p_block))) CYCLE
         IF (iblock_col == iblock_row) THEN
            DO j = 1, SIZE(p_block, 1)
               DO i = 1, SIZE(p_block, 2)
                  output(iblock_row) = output(iblock_row) + p_block(j, i)*bra(j, i)
               END DO
            END DO
         ELSE
            DO j = 1, SIZE(p_block, 1)
               DO i = 1, SIZE(p_block, 2)
                  output(iblock_row) = output(iblock_row) + p_block(j, i)*ket(j, i)
               END DO
            END DO
            DO j = 1, SIZE(p_block, 1)
               DO i = 1, SIZE(p_block, 2)
                  output(iblock_col) = output(iblock_col) + p_block(j, i)*bra(j, i)
               END DO
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL para_env%sum(output)

   END SUBROUTINE tb_contract_dens_matrix

! **************************************************************************************************
!> \brief compute the AO-resolved density contraction for one k-point spin channel.
!> \param p_matrix ...
!> \param bra_mat ...
!> \param ket_mat ...
!> \param iop ...
!> \param nops ...
!> \param output ...
!> \param ispin ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_contract_dens_kp_spin(p_matrix, bra_mat, ket_mat, iop, nops, output, ispin, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: p_matrix
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: bra_mat, ket_mat
      INTEGER, INTENT(IN)                                :: iop, nops
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: output
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: ic, idx, nimg
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: image_output
      TYPE(dbcsr_type), POINTER                          :: p_mat

      nimg = SIZE(p_matrix, 2)
      output = 0.0_dp
      ALLOCATE (image_output(SIZE(output)))
      DO ic = 1, nimg
         idx = iop + nops*(ic - 1)
         CPASSERT(idx <= SIZE(bra_mat))
         CPASSERT(idx <= SIZE(ket_mat))
         NULLIFY (p_mat)
         p_mat => p_matrix(ispin, ic)%matrix
         image_output = 0.0_dp
         CALL tb_contract_dens_matrix(p_mat, bra_mat(idx)%matrix, ket_mat(idx)%matrix, image_output, para_env)
         output = output + image_output
      END DO
      DEALLOCATE (image_output)

   END SUBROUTINE tb_contract_dens_kp_spin

! **************************************************************************************************
!> \brief compute the mulliken properties (AO resolved)
!> \param p_matrix ...
!> \param bra_mat ...
!> \param ket_mat ...
!> \param output ...
!> \param para_env ...
!> \par History
!>      adapted from ao_charges_2
!> \note
! **************************************************************************************************
   SUBROUTINE contract_dens(p_matrix, bra_mat, ket_mat, output, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_matrix
      TYPE(dbcsr_type), POINTER                          :: bra_mat, ket_mat
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: output
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'contract_dens'

      INTEGER                                            :: handle, i, iblock_col, iblock_row, &
                                                            ispin, j, nspin
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: bra, ket, p_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      nspin = SIZE(p_matrix)
      output = 0.0_dp
      DO ispin = 1, nspin
         CALL dbcsr_iterator_start(iter, bra_mat)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            NULLIFY (p_block, bra, ket)

            CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, bra)
            CALL dbcsr_get_block_p(matrix=p_matrix(ispin)%matrix, &
                                   row=iblock_row, col=iblock_col, BLOCK=p_block, found=found)
            IF (.NOT. found) CYCLE
            CALL dbcsr_get_block_p(matrix=ket_mat, &
                                   row=iblock_row, col=iblock_col, BLOCK=ket, found=found)
            IF (.NOT. found) CPABORT("missing block")

            IF (.NOT. (ASSOCIATED(bra) .AND. ASSOCIATED(p_block))) CYCLE
            IF (iblock_col == iblock_row) THEN
               DO j = 1, SIZE(p_block, 1)
                  DO i = 1, SIZE(p_block, 2)
                     output(iblock_row) = output(iblock_row) + p_block(j, i)*bra(j, i)
                  END DO
               END DO
            ELSE
               DO j = 1, SIZE(p_block, 1)
                  DO i = 1, SIZE(p_block, 2)
                     output(iblock_row) = output(iblock_row) + p_block(j, i)*ket(j, i)
                  END DO
               END DO
               DO j = 1, SIZE(p_block, 1)
                  DO i = 1, SIZE(p_block, 2)
                     output(iblock_col) = output(iblock_col) + p_block(j, i)*bra(j, i)
                  END DO
               END DO
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
      END DO

      CALL para_env%sum(output)
      CALL timestop(handle)

   END SUBROUTINE contract_dens

! **************************************************************************************************
!> \brief compute the AO-resolved density contraction for real-space k-point image matrices
!> \param p_matrix ...
!> \param bra_mat ...
!> \param ket_mat ...
!> \param iop ...
!> \param nops ...
!> \param output ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE contract_dens_kp(p_matrix, bra_mat, ket_mat, iop, nops, output, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: p_matrix
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: bra_mat, ket_mat
      INTEGER, INTENT(IN)                                :: iop, nops
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: output
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: ic, idx, nimg
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: image_output
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_image

      nimg = SIZE(p_matrix, 2)
      output = 0.0_dp
      ALLOCATE (image_output(SIZE(output)))
      DO ic = 1, nimg
         idx = iop + nops*(ic - 1)
         CPASSERT(idx <= SIZE(bra_mat))
         CPASSERT(idx <= SIZE(ket_mat))
         NULLIFY (p_image)
         p_image => p_matrix(:, ic)
         image_output = 0.0_dp
         CALL contract_dens(p_image, bra_mat(idx)%matrix, ket_mat(idx)%matrix, image_output, para_env)
         output = output + image_output
      END DO
      DEALLOCATE (image_output)

   END SUBROUTINE contract_dens_kp

! **************************************************************************************************
!> \brief save gradient to force
!> \param qs_env ...
!> \param tb ...
!> \param para_env ...
!> \param ityp ...
!> \note
! **************************************************************************************************
   SUBROUTINE tb_grad2force(qs_env, tb, para_env, ityp)

      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tblite_type)                                  :: tb
      TYPE(mp_para_env_type)                             :: para_env
      INTEGER                                            :: ityp

      CHARACTER(len=*), PARAMETER                        :: routineN = 'tb_grad2force'

      CHARACTER(LEN=default_path_length)                 :: dump_file
      INTEGER                                            :: atoma, dump_status, dump_unit, handle, &
                                                            iatom, ikind, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force

      CALL timeset(routineN, handle)

      NULLIFY (force, atomic_kind_set)
      CALL get_qs_env(qs_env=qs_env, force=force, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of)

      natom = SIZE(particle_set)

      dump_status = 1
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_FORCE_DUMP", dump_file, STATUS=dump_status)
#endif
      IF (dump_status == 0) THEN
         OPEN (NEWUNIT=dump_unit, FILE=TRIM(dump_file), STATUS="UNKNOWN", &
               POSITION="APPEND", ACTION="WRITE")
         WRITE (dump_unit, "(A,1X,I0)") "component", ityp
         DO iatom = 1, natom
            WRITE (dump_unit, "(I0,3(1X,ES24.16))") iatom, tb%grad(:, iatom)/para_env%num_pe
         END DO
         CLOSE (dump_unit)
      END IF

      SELECT CASE (ityp)
      CASE DEFAULT
         CPABORT("unknown force type")
      CASE (0)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%all_potential(:, atoma) = &
               force(ikind)%all_potential(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (1)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%repulsive(:, atoma) = &
               force(ikind)%repulsive(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (2)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%dispersion(:, atoma) = &
               force(ikind)%dispersion(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (3)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%rho_elec(:, atoma) = &
               force(ikind)%rho_elec(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (4)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%overlap(:, atoma) = &
               force(ikind)%overlap(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (5)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%efield(:, atoma) = &
               force(ikind)%efield(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE tb_grad2force

! **************************************************************************************************
!> \brief set gradient to zero
!> \param qs_env ...
!> \note
! **************************************************************************************************
   SUBROUTINE tb_zero_force(qs_env)

      TYPE(qs_environment_type)                          :: qs_env

      INTEGER                                            :: iatom, ikind, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force

      NULLIFY (force, atomic_kind_set)
      CALL get_qs_env(qs_env=qs_env, force=force, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               kind_of=kind_of)

      natom = SIZE(particle_set)

      DO iatom = 1, natom
         ikind = kind_of(iatom)
         force(ikind)%all_potential = 0.0_dp
         force(ikind)%repulsive = 0.0_dp
         force(ikind)%dispersion = 0.0_dp
         force(ikind)%rho_elec = 0.0_dp
         force(ikind)%overlap = 0.0_dp
         force(ikind)%efield = 0.0_dp
      END DO

   END SUBROUTINE tb_zero_force

! **************************************************************************************************
!> \brief Add SCC-overlap and direct multipole Hamiltonian derivatives.
!> \param qs_env ...
!> \param use_rho ...
!> \param nimg ...
! **************************************************************************************************
   SUBROUTINE tb_derive_dH_off(qs_env, use_rho, nimg)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: use_rho
      INTEGER, INTENT(IN)                                :: nimg

#if defined(__TBLITE)
      INTEGER                                            :: i, idim, ij, iatom, ic, icol, ikind, img, ispin, &
                                                            jdim, ni, nj, nkind, nel, &
                                                            ityp, jatom, jkind, jrow, jtyp, iset, jset, nseti, nsetj, &
                                                            ia, ib, inda, indb, sampled_axes, sampled_even_axes, &
                                                            sampled_gamma_axes, nspin, ikp_axis
      INTEGER, DIMENSION(3)                              :: cellind, nkp_cellind, nkp_grid
      INTEGER, DIMENSION(:), POINTER                     :: nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found, gamma_centered, gamma_sampled_image_pair, &
                                                            has_multipole_response, sampled_image_pair, &
                                                            use_matrix_scc_stress
      LOGICAL, DIMENSION(3)                              :: mesh_has_gamma
      REAL(KIND=dp)                                      :: r2, dr, i_a_shift, j_a_shift, i_a_shift_mag, j_a_shift_mag, &
                                                            ishift, jshift, ishift_mag, jshift_mag, pij_charge, &
                                                            pij_magnet, mp_pair_scale, kpoint_coordinate, native_dot_tmp
      REAL(KIND=dp), DIMENSION(3)                        :: kp_shift
      REAL(KIND=dp), DIMENSION(3)                        :: rij, dgrad, dhgrad_charge, dhgrad_magnet, &
                                                            mpgrad_charge, mpgrad_magnet
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hsigma
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: t_ov, idip, jdip, idip_mag, jdip_mag, &
                                                            iquad, jquad, iquad_mag, jquad_mag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: t_dip, t_quad, t_d_ov
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: t_i_dip, t_i_quad, t_j_dip, t_j_quad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :, :) :: scc_strain_hint
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pblock, pblock_beta
      TYPE(block_p_type), DIMENSION(3, 3, 2)             :: scc_strain_blocks

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_dh_scc, matrix_p
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER            :: basis_set_a, basis_set_b
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_kp
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(tblite_type), POINTER                         :: tb
      TYPE(virial_type), POINTER                         :: virial

      ! compute mulliken charges required for charge update
      NULLIFY (scf_env, rho, tb, sab_orb, sab_kp, para_env, kpoints, matrix_dh_scc, virial)
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      scf_env=scf_env, &
                      rho=rho, &
                      tb_tblite=tb, &
                      sab_orb=sab_orb, &
                      sab_kp=sab_kp, &
                      para_env=para_env, &
                      qs_kind_set=qs_kind_set, &
                      kpoints=kpoints)

      NULLIFY (cell_to_index)
      IF (nimg > 1) THEN
         IF (.NOT. ASSOCIATED(sab_kp)) CPABORT("Missing tblite k-point neighbor list")
         sab_orb => sab_kp
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index, nkp_grid=nkp_grid, &
                              kp_shift=kp_shift, gamma_centered=gamma_centered)
      ELSE
         nkp_grid = 1
         kp_shift = 0.0_dp
         gamma_centered = .TRUE.
      END IF
      mesh_has_gamma = .FALSE.
      DO i = 1, 3
         DO ikp_axis = 1, nkp_grid(i)
            IF (gamma_centered .AND. MODULO(nkp_grid(i), 2) == 0) THEN
               kpoint_coordinate = REAL(2*ikp_axis - nkp_grid(i), KIND=dp)/ &
                                   REAL(2*nkp_grid(i), KIND=dp) + kp_shift(i)
            ELSE
               kpoint_coordinate = REAL(2*ikp_axis - nkp_grid(i) - 1, KIND=dp)/ &
                                   REAL(2*nkp_grid(i), KIND=dp) + kp_shift(i)
            END IF
            IF (ABS(kpoint_coordinate - ANINT(kpoint_coordinate)) < 1.0e-12_dp) &
               mesh_has_gamma(i) = .TRUE.
         END DO
      END DO
      has_multipole_response = ASSOCIATED(tb%dipbra) .OR. ASSOCIATED(tb%quadbra)

      NULLIFY (matrix_p)
      IF (use_rho) THEN
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ELSEIF (ASSOCIATED(tb%rho_ao_kp_ref)) THEN
         matrix_p => tb%rho_ao_kp_ref
      ELSE
         matrix_p => scf_env%p_mix_new
      END IF
      nspin = SIZE(matrix_p, 1)
      use_matrix_scc_stress = tb%use_virial .AND. nimg > 1 .AND. ALL(tb%mol%periodic)

      IF (use_matrix_scc_stress) THEN
         CALL dbcsr_allocate_matrix_set(matrix_dh_scc, 9*nspin, nimg)
         DO img = 1, nimg
            DO ispin = 1, nspin
               DO jdim = 1, 3
                  DO idim = 1, 3
                     i = idim + 3*(jdim - 1) + 9*(ispin - 1)
                     ALLOCATE (matrix_dh_scc(i, img)%matrix)
                     CALL dbcsr_create(matrix_dh_scc(i, img)%matrix, template=matrix_p(ispin, img)%matrix, &
                                       name="TBLITE SCC STRAIN DERIVATIVE")
                     CALL cp_dbcsr_alloc_block_from_nbl(matrix_dh_scc(i, img)%matrix, sab_orb)
                  END DO
               END DO
            END DO
         END DO
      END IF

      ! set up basis set lists
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      nel = msao(tb%calc%bas%maxl)**2
      ALLOCATE (t_ov(nel))
      ALLOCATE (t_d_ov(3, nel))
      ALLOCATE (t_dip(dip_n, nel))
      ALLOCATE (t_i_dip(3, dip_n, nel), t_j_dip(3, dip_n, nel))
      ALLOCATE (t_quad(quad_n, nel))
      ALLOCATE (t_i_quad(3, quad_n, nel), t_j_quad(3, quad_n, nel))

      ALLOCATE (idip(dip_n), jdip(dip_n), idip_mag(dip_n), jdip_mag(dip_n))
      ALLOCATE (iquad(quad_n), jquad(quad_n), iquad_mag(quad_n), jquad_mag(quad_n))

      tb%grad = 0.0_dp
      hsigma = 0.0_dp
      ! loop over all atom pairs with a non-zero overlap (sab_orb)
      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cellind)

         icol = MAX(iatom, jatom)
         jrow = MIN(iatom, jatom)

         IF (iatom < jatom) THEN
            rij = -rij
            i = ikind
            ikind = jkind
            jkind = i
         END IF

         ityp = tb%mol%id(icol)
         jtyp = tb%mol%id(jrow)

         r2 = DOT_PRODUCT(rij, rij)
         dr = SQRT(r2)
         IF (icol == jrow .AND. dr < same_atom) CYCLE

         !get basis information
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         first_sgfa => basis_set_a%first_sgf
         nsgfa => basis_set_a%nsgf_set
         nseti = basis_set_a%nset
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         first_sgfb => basis_set_b%first_sgf
         nsgfb => basis_set_b%nsgf_set
         nsetj = basis_set_b%nset

         IF (nimg == 1) THEN
            ic = 1
         ELSE
            ic = cell_to_index(cellind(1), cellind(2), cellind(3))
            CPASSERT(ic > 0)
         END IF
         nkp_cellind = 0
         DO i = 1, 3
            IF (nkp_grid(i) > 1) THEN
               nkp_cellind(i) = MODULO(cellind(i), nkp_grid(i))
               IF (2*nkp_cellind(i) > nkp_grid(i)) nkp_cellind(i) = nkp_cellind(i) - nkp_grid(i)
            END IF
         END DO
         sampled_axes = 0
         IF (nkp_cellind(1) /= 0) sampled_axes = sampled_axes + 1
         IF (nkp_cellind(2) /= 0) sampled_axes = sampled_axes + 1
         IF (nkp_cellind(3) /= 0) sampled_axes = sampled_axes + 1
         sampled_even_axes = 0
         sampled_gamma_axes = 0
         DO i = 1, 3
            IF (nkp_grid(i) > 1 .AND. MODULO(nkp_grid(i), 2) == 0 .AND. &
                ABS(2*nkp_cellind(i)) == nkp_grid(i)) THEN
               sampled_even_axes = sampled_even_axes + 1
               IF (mesh_has_gamma(i)) sampled_gamma_axes = sampled_gamma_axes + 1
            END IF
         END DO
         sampled_image_pair = sampled_axes > 0 .AND. sampled_even_axes > 0
         gamma_sampled_image_pair = sampled_image_pair .AND. sampled_gamma_axes == sampled_even_axes
         mp_pair_scale = 1.0_dp
         IF (icol == jrow .AND. sampled_image_pair .AND. has_multipole_response .AND. &
             .NOT. gamma_sampled_image_pair) mp_pair_scale = -1.0_dp

         NULLIFY (pblock, pblock_beta)
         CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                row=jrow, col=icol, BLOCK=pblock, found=found)
         IF (.NOT. found) CPABORT("pblock not found")
         IF (nspin > 1) THEN
            CALL dbcsr_get_block_p(matrix=matrix_p(2, ic)%matrix, &
                                   row=jrow, col=icol, BLOCK=pblock_beta, found=found)
            IF (.NOT. found) CPABORT("pblock beta not found")
         END IF
         IF (use_matrix_scc_stress) THEN
            ALLOCATE (scc_strain_hint(SIZE(pblock, 2), SIZE(pblock, 1), 3, 3, nspin))
            scc_strain_hint = 0.0_dp
            DO ispin = 1, nspin
               DO jdim = 1, 3
                  DO idim = 1, 3
                     NULLIFY (scc_strain_blocks(idim, jdim, ispin)%block)
                     i = idim + 3*(jdim - 1) + 9*(ispin - 1)
                     CALL dbcsr_get_block_p(matrix=matrix_dh_scc(i, ic)%matrix, &
                                            row=jrow, col=icol, &
                                            BLOCK=scc_strain_blocks(idim, jdim, ispin)%block, found=found)
                     IF (.NOT. found) CPABORT("SCC strain derivative block not found")
                  END DO
               END DO
            END DO
         END IF
         i_a_shift = tb%pot%vat(icol, 1)
         j_a_shift = tb%pot%vat(jrow, 1)
         i_a_shift_mag = 0.0_dp
         j_a_shift_mag = 0.0_dp
         IF (SIZE(tb%pot%vat, 2) > 1) THEN
            i_a_shift_mag = tb%pot%vat(icol, 2)
            j_a_shift_mag = tb%pot%vat(jrow, 2)
         END IF
         idip(:) = tb%pot%vdp(:, icol, 1)
         jdip(:) = tb%pot%vdp(:, jrow, 1)
         idip_mag(:) = 0.0_dp
         jdip_mag(:) = 0.0_dp
         IF (SIZE(tb%pot%vdp, 3) > 1) THEN
            idip_mag(:) = tb%pot%vdp(:, icol, 2)
            jdip_mag(:) = tb%pot%vdp(:, jrow, 2)
         END IF
         iquad(:) = tb%pot%vqp(:, icol, 1)
         jquad(:) = tb%pot%vqp(:, jrow, 1)
         iquad_mag(:) = 0.0_dp
         jquad_mag(:) = 0.0_dp
         IF (SIZE(tb%pot%vqp, 3) > 1) THEN
            iquad_mag(:) = tb%pot%vqp(:, icol, 2)
            jquad_mag(:) = tb%pot%vqp(:, jrow, 2)
         END IF
         ni = tb%calc%bas%ish_at(icol)
         DO iset = 1, nseti
            ishift = i_a_shift + tb%pot%vsh(ni + iset, 1)
            ishift_mag = 0.0_dp
            IF (SIZE(tb%pot%vsh, 2) > 1) ishift_mag = i_a_shift_mag + tb%pot%vsh(ni + iset, 2)
            nj = tb%calc%bas%ish_at(jrow)
            DO jset = 1, nsetj
               jshift = j_a_shift + tb%pot%vsh(nj + jset, 1)
               jshift_mag = 0.0_dp
               IF (SIZE(tb%pot%vsh, 2) > 1) jshift_mag = j_a_shift_mag + tb%pot%vsh(nj + jset, 2)

               !get integrals and derivatives
#if defined(__TBLITE_GXTB)
               CALL cp2k_multipole_grad_cgto(tb%calc%bas, tb%bcache, iset, ityp, icol, &
                  & jset, jtyp, jrow, r2, rij, t_ov, t_dip, t_quad, t_d_ov, &
                  & t_i_dip, t_i_quad, t_j_dip, t_j_quad)
#else
               CALL multipole_grad_cgto(tb%calc%bas%cgto(iset, ityp), tb%calc%bas%cgto(jset, jtyp), &
                  & r2, rij, tb%calc%bas%intcut, t_ov, t_dip, t_quad, t_d_ov, t_i_dip, t_i_quad, &
                  & t_j_dip, t_j_quad)
#endif

               dgrad = 0.0_dp
               DO inda = 1, nsgfa(iset)
                  ia = first_sgfa(1, iset) - first_sgfa(1, 1) + inda
                  DO indb = 1, nsgfb(jset)
                     ib = first_sgfb(1, jset) - first_sgfb(1, 1) + indb

                     ij = inda + nsgfa(iset)*(indb - 1)

                     pij_charge = pblock(ib, ia)
                     pij_magnet = 0.0_dp
                     IF (nspin > 1) THEN
                        pij_charge = pij_charge + pblock_beta(ib, ia)
                        pij_magnet = pblock(ib, ia) - pblock_beta(ib, ia)
                     END IF
                     mpgrad_charge = MATMUL(t_i_dip(:, :, ij), idip) &
                                     + MATMUL(t_j_dip(:, :, ij), jdip) &
                                     + MATMUL(t_i_quad(:, :, ij), iquad) &
                                     + MATMUL(t_j_quad(:, :, ij), jquad)
                     mpgrad_magnet = MATMUL(t_i_dip(:, :, ij), idip_mag) &
                                     + MATMUL(t_j_dip(:, :, ij), jdip_mag) &
                                     + MATMUL(t_i_quad(:, :, ij), iquad_mag) &
                                     + MATMUL(t_j_quad(:, :, ij), jquad_mag)
                     dhgrad_charge = -(ishift + jshift)*t_d_ov(:, ij) - mp_pair_scale*mpgrad_charge
                     dhgrad_magnet = -(ishift_mag + jshift_mag)*t_d_ov(:, ij) - mp_pair_scale*mpgrad_magnet
                     dgrad(:) = dgrad(:) - &
                                ((ishift + jshift)*pij_charge + &
                                 (ishift_mag + jshift_mag)*pij_magnet)*t_d_ov(:, ij) - &
                                mp_pair_scale*(pij_charge*mpgrad_charge + pij_magnet*mpgrad_magnet)

                     IF (ALLOCATED(scc_strain_hint)) THEN
                        DO jdim = 1, 3
                           DO idim = 1, 3
                              scc_strain_hint(ia, ib, idim, jdim, 1) = &
                                 scc_strain_hint(ia, ib, idim, jdim, 1) &
                                 + (dhgrad_charge(idim) + MERGE(dhgrad_magnet(idim), 0.0_dp, nspin > 1))*rij(jdim)
                              IF (nspin > 1) THEN
                                 scc_strain_hint(ia, ib, idim, jdim, 2) = &
                                    scc_strain_hint(ia, ib, idim, jdim, 2) &
                                    + (dhgrad_charge(idim) - dhgrad_magnet(idim))*rij(jdim)
                              END IF
                           END DO
                        END DO
                     END IF

                  END DO
               END DO
               tb%grad(:, icol) = tb%grad(:, icol) - dgrad
               tb%grad(:, jrow) = tb%grad(:, jrow) + dgrad
               IF (tb%use_virial .AND. .NOT. use_matrix_scc_stress) THEN
                  IF (icol == jrow) THEN
                     DO ia = 1, 3
                        DO ib = 1, 3
                           IF (sampled_image_pair .AND. .NOT. gamma_sampled_image_pair) THEN
                              hsigma(ia, ib) = hsigma(ia, ib) - 0.25_dp* &
                                               (rij(ia)*dgrad(ib) + rij(ib)*dgrad(ia))
                           ELSE
                              hsigma(ia, ib) = hsigma(ia, ib) + 0.25_dp* &
                                               (rij(ia)*dgrad(ib) + rij(ib)*dgrad(ia))
                           END IF
                        END DO
                     END DO
                  ELSE
                     DO ia = 1, 3
                        DO ib = 1, 3
                           hsigma(ia, ib) = hsigma(ia, ib) + 0.50_dp*(rij(ia)*dgrad(ib) + rij(ib)*dgrad(ia))
                        END DO
                     END DO
                  END IF
               END IF
            END DO
         END DO
         IF (ALLOCATED(scc_strain_hint)) THEN
            DO ispin = 1, nspin
               DO jdim = 1, 3
                  DO idim = 1, 3
                     IF (icol <= jrow) THEN
                        scc_strain_blocks(idim, jdim, ispin)%block(:, :) = &
                           scc_strain_blocks(idim, jdim, ispin)%block(:, :) &
                           + scc_strain_hint(:, :, idim, jdim, ispin)
                     ELSE
                        scc_strain_blocks(idim, jdim, ispin)%block(:, :) = &
                           scc_strain_blocks(idim, jdim, ispin)%block(:, :) &
                           + TRANSPOSE(scc_strain_hint(:, :, idim, jdim, ispin))
                     END IF
                  END DO
               END DO
            END DO
            DEALLOCATE (scc_strain_hint)
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      IF (use_matrix_scc_stress) THEN
         hsigma = 0.0_dp
         DO img = 1, nimg
            DO ispin = 1, nspin
               DO jdim = 1, 3
                  DO idim = 1, 3
                     i = idim + 3*(jdim - 1) + 9*(ispin - 1)
                     CALL dbcsr_finalize(matrix_dh_scc(i, img)%matrix)
                     CALL dbcsr_dot(matrix_dh_scc(i, img)%matrix, matrix_p(ispin, img)%matrix, native_dot_tmp)
                     ! Symmetric DBCSR storage represents both orientations of each assembled pair.
                     hsigma(idim, jdim) = hsigma(idim, jdim) + 0.5_dp*native_dot_tmp
                  END DO
               END DO
            END DO
         END DO
         CALL dbcsr_deallocate_matrix_set(matrix_dh_scc)
      ELSE
         CALL para_env%sum(hsigma)
      END IF
      CALL para_env%sum(tb%grad)
      CALL tb_grad2force(qs_env, tb, para_env, 4)

      IF (.NOT. use_matrix_scc_stress) tb%sigma = tb%sigma + hsigma

      DEALLOCATE (basis_set_list)
      DEALLOCATE (t_ov, t_d_ov)
      DEALLOCATE (t_dip, t_i_dip, t_j_dip)
      DEALLOCATE (t_quad, t_i_quad, t_j_quad)
      DEALLOCATE (idip, jdip, idip_mag, jdip_mag, iquad, jquad, iquad_mag, jquad_mag)

      IF (tb%use_virial) THEN
         CALL tb_add_stress(qs_env, tb, para_env)
         IF (use_matrix_scc_stress) THEN
            CALL get_qs_env(qs_env=qs_env, virial=virial)
            virial%pv_virial = virial%pv_virial - hsigma/para_env%num_pe
         END IF
      END IF

#else
      MARK_USED(qs_env)
      MARK_USED(use_rho)
      MARK_USED(nimg)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_derive_dH_off

! **************************************************************************************************
!> \brief Run native tblite CLI and compare against CP2K/tblite.
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE tb_reference_cli_compare(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tb_reference_cli_compare'

      CHARACTER(LEN=16)                                  :: solvation_model_name
      CHARACTER(LEN=32) :: acc_str, charge_str, efield_x_str, efield_y_str, efield_z_str, &
         etemp_guess_val_str, etemp_str, iter_str, spin_str, spinpol_str
      CHARACTER(LEN=4*default_path_length+16) :: efield_str, etemp_guess_str, guess_str, &
         param_str, post_processing_output_str, post_processing_str, restart_str, solvation_str, &
         verbosity_str
      CHARACTER(LEN=8)                                   :: guess, method, solver
      CHARACTER(LEN=8*default_path_length)               :: command
      CHARACTER(LEN=default_path_length)                 :: file_base, gen_file, grad_file, &
                                                            json_file, log_file, &
                                                            post_processing_output_file
      INTEGER                                            :: cmdstat, exitstat, handle, iounit, &
                                                            n_periodic, natom, nkp, &
                                                            reference_iterations, spin
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: do_kpoints, have_energy, have_gradient, &
                                                            have_virial, too_large, &
                                                            unsupported_kpoints
      REAL(KIND=dp)                                      :: cli_energy, cp_energy, ediff, etemp, &
                                                            etemp_guess, fmax, fsum, vmax, vsum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cli_gradient, cli_virial, cp_gradient
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_reference_cli_type)                       :: ref

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, cell, dft_control, energy, force, kpoints, logger, para_env, particle_set, &
               scf_control, virial, xkp)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                      dft_control=dft_control, energy=energy, force=force, &
                      para_env=para_env, particle_set=particle_set, scf_control=scf_control, &
                      virial=virial, do_kpoints=do_kpoints, kpoints=kpoints)

      ref = dft_control%qs_control%xtb_control%reference_cli
      IF (.NOT. ref%enabled) THEN
         CALL timestop(handle)
         RETURN
      END IF
      IF (.NOT. para_env%is_source()) THEN
         CALL timestop(handle)
         RETURN
      END IF

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      verbosity_str = ""
      IF (logger%iter_info%print_level == silent_print_level) THEN
         verbosity_str = " --silent"
      ELSE IF (logger%iter_info%print_level == high_print_level .OR. &
               logger%iter_info%print_level == debug_print_level) THEN
         verbosity_str = " --verbose"
      END IF
      periodic = 0
      IF (ASSOCIATED(cell)) CALL get_cell(cell=cell, periodic=periodic)
      n_periodic = COUNT(periodic == 1)
      IF (n_periodic == 1 .OR. n_periodic == 2) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "tblite reference CLI check skipped: partial periodicity cannot be represented faithfully."
         WRITE (UNIT=iounit, FMT="(T2,A)") &
            "DFTB+ GEN distinguishes only molecular and fully three-dimensional periodic geometries."
         IF (ref%stop_on_error) &
            CPABORT("REFERENCE_CLI cannot represent one- or two-dimensional periodicity")
         CALL timestop(handle)
         RETURN
      END IF
      IF (ref%solvation_active) THEN
         IF (n_periodic == 3) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") &
               "tblite reference CLI implicit solvation is not supported for PERIODIC XYZ."
            WRITE (UNIT=iounit, FMT="(T2,A)") &
               "Use PERIODIC NONE for molecular solvation diagnostics, or remove IMPLICIT_SOLVATION."
            CPABORT("REFERENCE_CLI implicit solvation is incompatible with PERIODIC XYZ")
         ELSE IF (n_periodic > 0) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") &
               "WARNING: tblite reference CLI implicit solvation with finite periodicity is diagnostic only."
            WRITE (UNIT=iounit, FMT="(T2,A,I0,A)") &
               "The generated native tblite reference geometry has ", n_periodic, &
               " periodic direction(s); continuum-solvation conventions are primarily molecular."
         END IF
      END IF
      unsupported_kpoints = .FALSE.
      IF (do_kpoints .AND. ASSOCIATED(kpoints)) THEN
         nkp = 0
         CALL get_kpoint_info(kpoint=kpoints, nkp=nkp, xkp=xkp)
         unsupported_kpoints = nkp > 1
         IF (nkp == 1 .AND. ASSOCIATED(xkp)) unsupported_kpoints = ANY(ABS(xkp(:, 1)) > 1.0E-12_dp)
      END IF
      IF (unsupported_kpoints) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "tblite reference CLI check skipped: CP2K KPOINTS are active."
         WRITE (UNIT=iounit, FMT="(T2,A)") &
            "The native tblite CLI reference path does not reproduce CP2K multi-k-point sampling."
         IF (ref%stop_on_error) CPABORT("tblite reference CLI cannot check CP2K k-point calculations")
         CALL timestop(handle)
         RETURN
      END IF
      IF (dft_control%qs_control%xtb_control%tblite_scc_mixer == tblite_scc_mixer_cp2k) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "WARNING: tblite reference CLI cannot reproduce XTB/SCC_MIXER CP2K."
         IF (dft_control%qs_control%xtb_control%tblite_method == gxtb) THEN
            WRITE (UNIT=iounit, FMT="(T2,A)") &
               "The external run uses g-xTB potential DIIS, not CP2K P/Fock mixing; compare converged results only."
         ELSE
            WRITE (UNIT=iounit, FMT="(T2,A)") &
               "The external native tblite run uses tblite's own SCC mixer; only the converged result is compared."
         END IF
         IF (ref%stop_on_error) CPABORT("tblite reference CLI cannot reproduce SCC_MIXER CP2K")
      END IF
      natom = SIZE(particle_set)
      method = tb_reference_method_name(dft_control%qs_control%xtb_control%tblite_method)
      guess = tb_reference_guess_name(ref%guess)
      guess_str = ""
      IF (ref%guess /= tblite_guess_auto) guess_str = " --guess "//TRIM(guess)
      solver = tb_reference_solver_name(dft_control%qs_control%xtb_control%tblite_mixer_solver)
      file_base = tb_join_path(ref%work_directory, ref%prefix)
      gen_file = TRIM(file_base)//".gen"
      grad_file = TRIM(file_base)//".grad"
      json_file = TRIM(file_base)//".json"
      log_file = TRIM(file_base)//".log"
      post_processing_output_file = ""
      IF (LEN_TRIM(ref%grad_file) > 0) grad_file = ref%grad_file
      IF (LEN_TRIM(ref%json_file) > 0) json_file = ref%json_file
      IF (LEN_TRIM(ref%post_processing_output_file) > 0) &
         post_processing_output_file = ref%post_processing_output_file

      WRITE (charge_str, "(I0)") dft_control%charge
      spin = MAX(0, dft_control%multiplicity - 1)
      WRITE (spin_str, "(I0)") spin
      WRITE (acc_str, "(ES16.8)") dft_control%qs_control%xtb_control%tblite_accuracy
      reference_iterations = dft_control%qs_control%xtb_control%tblite_mixer_iterations
      WRITE (iter_str, "(I0)") reference_iterations
      efield_str = ""
      IF (ref%efield_active) THEN
         WRITE (efield_x_str, "(ES16.8)") ref%efield(1)
         WRITE (efield_y_str, "(ES16.8)") ref%efield(2)
         WRITE (efield_z_str, "(ES16.8)") ref%efield(3)
         efield_str = " --efield "//TRIM(ADJUSTL(efield_x_str))//","// &
                      TRIM(ADJUSTL(efield_y_str))//","//TRIM(ADJUSTL(efield_z_str))
      END IF
      solvation_str = ""
      solvation_model_name = ""
      IF (ref%solvation_active) THEN
         SELECT CASE (ref%solvation_model)
         CASE (tblite_cli_solvation_alpb)
            solvation_model_name = "ALPB"
            solvation_str = " --alpb "
         CASE (tblite_cli_solvation_gbsa)
            solvation_model_name = "GBSA"
            solvation_str = " --gbsa "
         CASE (tblite_cli_solvation_gbe)
            solvation_model_name = "GBE"
            solvation_str = " --gbe "
         CASE (tblite_cli_solvation_gb)
            solvation_model_name = "GB"
            solvation_str = " --gb "
         CASE (tblite_cli_solvation_cpcm)
            solvation_model_name = "CPCM"
            solvation_str = " --cpcm "
         CASE DEFAULT
            CPABORT("Unknown tblite reference CLI implicit-solvation model")
         END SELECT
         solvation_str = TRIM(solvation_str)//" "//TRIM(tb_shell_quote(ref%solvation_solvent))
         SELECT CASE (ref%solvation_born_kernel)
         CASE (tblite_cli_born_kernel_auto)
         CASE (tblite_cli_born_kernel_p16)
            solvation_str = TRIM(solvation_str)//" --born-kernel p16"
         CASE (tblite_cli_born_kernel_still)
            solvation_str = TRIM(solvation_str)//" --born-kernel still"
         CASE DEFAULT
            CPABORT("Unknown tblite reference CLI Born kernel")
         END SELECT
         SELECT CASE (ref%solvation_state)
         CASE (tblite_cli_solution_state_gsolv)
         CASE (tblite_cli_solution_state_bar1mol)
            solvation_str = TRIM(solvation_str)//" --solv-state bar1mol"
         CASE (tblite_cli_solution_state_reference)
            solvation_str = TRIM(solvation_str)//" --solv-state reference"
         CASE DEFAULT
            CPABORT("Unknown tblite reference CLI solution state")
         END SELECT
      END IF
      IF (dft_control%qs_control%xtb_control%tblite_method /= gxtb) THEN
         IF (dft_control%qs_control%xtb_control%tblite_mixer_memory /= reference_iterations) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") &
               "WARNING: tblite reference CLI cannot reproduce XTB/TBLITE_MIXER/MEMORY."
            WRITE (UNIT=iounit, FMT="(T2,A,I0,A,I0,A)") &
               "The native reference run uses tblite's internal mixer memory tied to --iterations (", &
               reference_iterations, "), while CP2K uses MEMORY ", &
               dft_control%qs_control%xtb_control%tblite_mixer_memory, "."
            IF (ref%stop_on_error) &
               CPABORT("tblite reference CLI cannot reproduce TBLITE_MIXER/MEMORY")
         END IF
         IF (ABS(dft_control%qs_control%xtb_control%tblite_mixer_damping - &
                 tblite_mixer_damping_default) > 10.0_dp*EPSILON(1.0_dp)) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") &
               "WARNING: tblite reference CLI cannot reproduce XTB/TBLITE_MIXER/DAMPING."
            WRITE (UNIT=iounit, FMT="(T2,A,F8.4,A,F8.4,A)") &
               "The native reference run uses tblite's library default ", tblite_mixer_damping_default, &
               ", while CP2K uses DAMPING ", dft_control%qs_control%xtb_control%tblite_mixer_damping, "."
            IF (ref%stop_on_error) &
               CPABORT("tblite reference CLI cannot reproduce TBLITE_MIXER/DAMPING")
         END IF
         IF (ABS(dft_control%qs_control%xtb_control%tblite_mixer_omega0 - &
                 tblite_mixer_omega0_default) > 10.0_dp*EPSILON(1.0_dp)) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") &
               "WARNING: tblite reference CLI cannot reproduce XTB/TBLITE_MIXER/OMEGA0."
            WRITE (UNIT=iounit, FMT="(T2,A,ES12.4,A,ES12.4,A)") &
               "The native reference run uses tblite's library default ", tblite_mixer_omega0_default, &
               ", while CP2K uses OMEGA0 ", dft_control%qs_control%xtb_control%tblite_mixer_omega0, "."
            IF (ref%stop_on_error) &
               CPABORT("tblite reference CLI cannot reproduce TBLITE_MIXER/OMEGA0")
         END IF
         IF (ABS(dft_control%qs_control%xtb_control%tblite_mixer_min_weight - &
                 tblite_mixer_min_weight_default) > 10.0_dp*EPSILON(1.0_dp)) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") &
               "WARNING: tblite reference CLI cannot reproduce XTB/TBLITE_MIXER/MIN_WEIGHT."
            WRITE (UNIT=iounit, FMT="(T2,A,ES12.4,A,ES12.4,A)") &
               "The native reference run uses tblite's library default ", tblite_mixer_min_weight_default, &
               ", while CP2K uses MIN_WEIGHT ", dft_control%qs_control%xtb_control%tblite_mixer_min_weight, "."
            IF (ref%stop_on_error) &
               CPABORT("tblite reference CLI cannot reproduce TBLITE_MIXER/MIN_WEIGHT")
         END IF
         IF (ABS(dft_control%qs_control%xtb_control%tblite_mixer_max_weight - &
                 tblite_mixer_max_weight_default) > 10.0_dp*EPSILON(1.0_dp)) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") &
               "WARNING: tblite reference CLI cannot reproduce XTB/TBLITE_MIXER/MAX_WEIGHT."
            WRITE (UNIT=iounit, FMT="(T2,A,ES12.4,A,ES12.4,A)") &
               "The native reference run uses tblite's library default ", tblite_mixer_max_weight_default, &
               ", while CP2K uses MAX_WEIGHT ", dft_control%qs_control%xtb_control%tblite_mixer_max_weight, "."
            IF (ref%stop_on_error) &
               CPABORT("tblite reference CLI cannot reproduce TBLITE_MIXER/MAX_WEIGHT")
         END IF
         IF (ABS(dft_control%qs_control%xtb_control%tblite_mixer_weight_factor - &
                 tblite_mixer_weight_factor_default) > 10.0_dp*EPSILON(1.0_dp)) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") &
               "WARNING: tblite reference CLI cannot reproduce XTB/TBLITE_MIXER/WEIGHT_FACTOR."
            WRITE (UNIT=iounit, FMT="(T2,A,ES12.4,A,ES12.4,A)") &
               "The native reference run uses tblite's library default ", tblite_mixer_weight_factor_default, &
               ", while CP2K uses WEIGHT_FACTOR ", &
               dft_control%qs_control%xtb_control%tblite_mixer_weight_factor, "."
            IF (ref%stop_on_error) &
               CPABORT("tblite reference CLI cannot reproduce TBLITE_MIXER/WEIGHT_FACTOR")
         END IF
      END IF
      etemp = 300.0_dp
      IF (ASSOCIATED(scf_control) .AND. ASSOCIATED(scf_control%smear)) THEN
         IF (scf_control%smear%do_smear) THEN
            etemp = cp_unit_from_cp2k(scf_control%smear%electronic_temperature, "K")
            IF (scf_control%smear%method /= smear_fermi_dirac) THEN
               WRITE (UNIT=iounit, FMT="(/,T2,A,A,A)") &
                  "WARNING: tblite reference CLI cannot reproduce CP2K smearing method ", &
                  TRIM(tb_reference_smear_method_name(scf_control%smear%method)), "."
               WRITE (UNIT=iounit, FMT="(T2,A,F12.3,A)") &
                  "The native reference run uses Fermi-Dirac electronic temperature ", etemp, " K instead."
            END IF
         END IF
      END IF
      WRITE (etemp_str, "(ES16.8)") etemp
      etemp_guess = 0.0_dp
      etemp_guess_str = ""
      IF (ref%electronic_temperature_guess > 0.0_dp) THEN
         etemp_guess = cp_unit_from_cp2k(ref%electronic_temperature_guess, "K")
         WRITE (etemp_guess_val_str, "(ES16.8)") etemp_guess
         etemp_guess_str = " --etemp-guess "//TRIM(ADJUSTL(etemp_guess_val_str))
      END IF
      param_str = ""
      IF (LEN_TRIM(dft_control%qs_control%xtb_control%tblite_param_file) > 0) &
         param_str = " --param "//TRIM(tb_shell_quote(dft_control%qs_control%xtb_control%tblite_param_file))
      spinpol_str = ""
      IF (dft_control%lsd) spinpol_str = " --spin-polarized"
      post_processing_str = ""
      IF (LEN_TRIM(ref%post_processing) > 0) &
         post_processing_str = " --post-processing "//TRIM(tb_shell_quote(ref%post_processing))
      post_processing_output_str = ""
      IF (LEN_TRIM(post_processing_output_file) > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "WARNING: tblite reference CLI POST_PROCESSING_OUTPUT was requested explicitly."
         WRITE (UNIT=iounit, FMT="(T2,A)") &
            "Some tblite 0.5.0 command-line builds document --post-processing-output but do not parse it."
         post_processing_output_str = " --post-processing-output "// &
                                      TRIM(tb_shell_quote(post_processing_output_file))
      END IF
      restart_str = " --no-restart"
      IF (LEN_TRIM(ref%restart_file) > 0) &
         restart_str = " --restart "//TRIM(tb_shell_quote(ref%restart_file))

      CALL tb_write_reference_gen(qs_env, TRIM(gen_file))

      command = TRIM(tb_shell_quote(ref%program_name))//" run --method "//TRIM(method)// &
                TRIM(param_str)// &
                TRIM(spinpol_str)// &
                " --charge "//TRIM(ADJUSTL(charge_str))// &
                " --spin "//TRIM(ADJUSTL(spin_str))// &
                " --acc "//TRIM(ADJUSTL(acc_str))// &
                TRIM(guess_str)// &
                " --solver "//TRIM(solver)// &
                " --iterations "//TRIM(ADJUSTL(iter_str))// &
                " --etemp "//TRIM(ADJUSTL(etemp_str))// &
                TRIM(etemp_guess_str)// &
                TRIM(efield_str)// &
                TRIM(solvation_str)// &
                TRIM(post_processing_str)// &
                TRIM(post_processing_output_str)// &
                TRIM(restart_str)// &
                TRIM(verbosity_str)//" --input "//TRIM(tb_shell_quote(ref%input_format))// &
                " --grad "//TRIM(tb_shell_quote(grad_file))// &
                " --json "//TRIM(tb_shell_quote(json_file))//" "//TRIM(tb_shell_quote(gen_file))// &
                " > "//TRIM(tb_shell_quote(log_file))//" 2>&1"

      cmdstat = 0
      exitstat = 0
      CALL execute_command_line(TRIM(command), exitstat=exitstat, cmdstat=cmdstat)
      IF (cmdstat /= 0 .OR. exitstat /= 0) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") "tblite reference CLI check failed to run."
         WRITE (UNIT=iounit, FMT="(T2,A,A)") "Command: ", TRIM(command)
         WRITE (UNIT=iounit, FMT="(T2,A,I0,T32,A,I0)") "cmdstat:", cmdstat, "exitstat:", exitstat
         IF (ref%stop_on_error) CPABORT("tblite reference CLI command failed")
         CALL tb_reference_cleanup(ref, gen_file, grad_file, json_file, log_file, post_processing_output_file)
         CALL timestop(handle)
         RETURN
      END IF

      ALLOCATE (cli_gradient(3, natom), cli_virial(3, 3))
      CALL tb_read_reference_grad(TRIM(grad_file), natom, cli_energy, cli_gradient, cli_virial, &
                                  have_energy, have_gradient, have_virial)

      WRITE (UNIT=iounit, FMT="(/,T2,A)") "tblite reference CLI check"
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Executable: ", TRIM(ref%program_name)
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Method:     ", TRIM(method)
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Guess:      ", TRIM(guess)
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Solver:     ", TRIM(solver)
      WRITE (UNIT=iounit, FMT="(T2,A,L1)") "Spin-pol.:  ", dft_control%lsd
      IF (ref%efield_active) &
         WRITE (UNIT=iounit, FMT="(T2,A,3ES16.8,A)") "Efield:     ", ref%efield, " V/Angstrom"
      IF (ref%solvation_active) &
         WRITE (UNIT=iounit, FMT="(T2,A,A,1X,A)") "Solvation:  ", TRIM(solvation_model_name), &
         TRIM(ref%solvation_solvent)
      IF (LEN_TRIM(ref%post_processing) > 0) &
         WRITE (UNIT=iounit, FMT="(T2,A,A)") "Post proc.: ", TRIM(ref%post_processing)
      IF (LEN_TRIM(post_processing_output_file) > 0) &
         WRITE (UNIT=iounit, FMT="(T2,A,A)") "PP output:  ", TRIM(post_processing_output_file)
      IF (ref%electronic_temperature_guess > 0.0_dp) &
         WRITE (UNIT=iounit, FMT="(T2,A,F12.3,A)") "Guess etemp:", etemp_guess, " K"
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Grad file:  ", TRIM(grad_file)
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "JSON file:  ", TRIM(json_file)
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Log file:   ", TRIM(log_file)

      too_large = .FALSE.
      IF (ref%check_energy) THEN
         IF (have_energy) THEN
            cp_energy = energy%total
            ediff = ABS(cp_energy - cli_energy)
            WRITE (UNIT=iounit, FMT="(T2,A,3ES22.12)") &
               "Energy CP2K/CLI/absdiff:", cp_energy, cli_energy, ediff
            too_large = too_large .OR. ediff > ref%error_limit
         ELSE
            WRITE (UNIT=iounit, FMT="(T2,A)") "Energy check skipped: no CLI energy found."
            too_large = .TRUE.
         END IF
      END IF

      IF (ref%check_forces) THEN
         IF (have_gradient .AND. ASSOCIATED(force)) THEN
            ALLOCATE (cp_gradient(3, natom))
            CALL total_qs_force(cp_gradient, force, atomic_kind_set)
            fsum = SUM(ABS(cp_gradient - cli_gradient))
            fmax = MAXVAL(ABS(cp_gradient - cli_gradient))
            WRITE (UNIT=iounit, FMT="(T2,A,2ES22.12)") "Gradient diff sum/max:", fsum, fmax
            too_large = too_large .OR. fmax > ref%error_limit
            DEALLOCATE (cp_gradient)
         ELSE IF (ASSOCIATED(force)) THEN
            WRITE (UNIT=iounit, FMT="(T2,A)") "Gradient check failed: no CLI gradient found."
            too_large = .TRUE.
         ELSE
            WRITE (UNIT=iounit, FMT="(T2,A)") &
               "Gradient check skipped: this CP2K calculation did not request forces."
         END IF
      END IF

      IF (ref%check_virial) THEN
         IF (have_virial .AND. ASSOCIATED(virial)) THEN
            ! Native tblite prints the positive cell derivative; CP2K stores the PV virial
            ! with the opposite sign.
            vsum = SUM(ABS(-virial%pv_virial - cli_virial))
            vmax = MAXVAL(ABS(-virial%pv_virial - cli_virial))
            WRITE (UNIT=iounit, FMT="(T2,A,2ES22.12)") "Virial diff sum/max:", vsum, vmax
            too_large = too_large .OR. vmax > ref%error_limit
         ELSE IF (ASSOCIATED(virial)) THEN
            WRITE (UNIT=iounit, FMT="(T2,A)") "Virial check failed: no CLI virial found."
            too_large = .TRUE.
         ELSE
            WRITE (UNIT=iounit, FMT="(T2,A)") &
               "Virial check skipped: this CP2K calculation did not request a virial."
         END IF
      END IF

      IF (too_large) THEN
         WRITE (UNIT=iounit, FMT="(T2,A,ES12.4)") &
            "tblite reference CLI deviation exceeded ERROR_LIMIT = ", ref%error_limit
         IF (ref%stop_on_error) CPABORT("tblite reference CLI deviation exceeded ERROR_LIMIT")
      END IF

      CALL tb_reference_cli_aux_commands(ref, dft_control, gen_file, file_base, verbosity_str, iounit)

      CALL tb_reference_cleanup(ref, gen_file, grad_file, json_file, log_file, post_processing_output_file)
      DEALLOCATE (cli_gradient, cli_virial)

      CALL timestop(handle)

   END SUBROUTINE tb_reference_cli_compare

! **************************************************************************************************
!> \brief Map CP2K tblite method id to native tblite CLI method name.
!> \param method_id ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_reference_method_name(method_id) RESULT(method)
      INTEGER, INTENT(IN)                                :: method_id
      CHARACTER(LEN=8)                                   :: method

      SELECT CASE (method_id)
      CASE (gfn1xtb)
         method = "gfn1"
      CASE (gfn2xtb)
         method = "gfn2"
      CASE (ipea1xtb)
         method = "ipea1"
      CASE (gxtb)
         method = "gxtb"
      CASE DEFAULT
         CPABORT("Unknown tblite reference CLI method")
      END SELECT

   END FUNCTION tb_reference_method_name

! **************************************************************************************************
!> \brief Map CP2K tblite reference CLI guess id to native tblite CLI guess name.
!> \param guess_id ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_reference_guess_name(guess_id) RESULT(guess)
      INTEGER, INTENT(IN)                                :: guess_id
      CHARACTER(LEN=8)                                   :: guess

      SELECT CASE (guess_id)
      CASE (tblite_guess_auto)
         guess = "auto"
      CASE (tblite_guess_sad)
         guess = "sad"
      CASE (tblite_guess_eeq)
         guess = "eeq"
      CASE (tblite_guess_ceh)
         guess = "ceh"
      CASE DEFAULT
         CPABORT("Unknown tblite reference CLI guess")
      END SELECT

   END FUNCTION tb_reference_guess_name

! **************************************************************************************************
!> \brief Map CP2K tblite solver id to native tblite CLI solver name.
!> \param solver_id ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_reference_solver_name(solver_id) RESULT(solver)
      INTEGER, INTENT(IN)                                :: solver_id
      CHARACTER(LEN=8)                                   :: solver

      SELECT CASE (solver_id)
      CASE (tblite_solver_gvd)
         solver = "gvd"
      CASE (tblite_solver_gvr)
         solver = "gvr"
      CASE DEFAULT
         CPABORT("Unknown tblite reference CLI solver")
      END SELECT

   END FUNCTION tb_reference_solver_name

! **************************************************************************************************
!> \brief Map CP2K smearing method id to a diagnostic label.
!> \param method_id ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_reference_smear_method_name(method_id) RESULT(method)
      INTEGER, INTENT(IN)                                :: method_id
      CHARACTER(LEN=24)                                  :: method

      SELECT CASE (method_id)
      CASE (smear_fermi_dirac)
         method = "FERMI_DIRAC"
      CASE (smear_energy_window)
         method = "ENERGY_WINDOW"
      CASE (smear_list)
         method = "LIST"
      CASE (smear_gaussian)
         method = "GAUSSIAN"
      CASE (smear_mp)
         method = "METHFESSEL_PAXTON"
      CASE (smear_mv)
         method = "MARZARI_VANDERBILT"
      CASE DEFAULT
         method = "UNKNOWN"
      END SELECT

   END FUNCTION tb_reference_smear_method_name

! **************************************************************************************************
!> \brief Run optional native tblite auxiliary subcommands.
!> \param ref ...
!> \param dft_control ...
!> \param gen_file ...
!> \param file_base ...
!> \param verbosity_str ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE tb_reference_cli_aux_commands(ref, dft_control, gen_file, file_base, verbosity_str, iounit)

      TYPE(xtb_reference_cli_type), INTENT(IN)           :: ref
      TYPE(dft_control_type), INTENT(IN)                 :: dft_control
      CHARACTER(LEN=*), INTENT(IN)                       :: gen_file, file_base, verbosity_str
      INTEGER, INTENT(IN)                                :: iounit

      CHARACTER(LEN=32)                                  :: charge_str, efield_x_str, efield_y_str, &
                                                            efield_z_str, etemp_guess_val_str, &
                                                            spin_str
      CHARACTER(LEN=4*default_path_length+16)            :: copy_str, dry_run_str, efield_str, &
                                                            etemp_guess_str, grad_str, json_str, &
                                                            method_str, output_str
      CHARACTER(LEN=8*default_path_length)               :: command
      CHARACTER(LEN=default_path_length)                 :: guess_input, log_file
      INTEGER                                            :: spin
      REAL(KIND=dp)                                      :: etemp_guess

      WRITE (charge_str, "(I0)") dft_control%charge
      spin = MAX(0, dft_control%multiplicity - 1)
      WRITE (spin_str, "(I0)") spin

      IF (ref%guess_cli%enabled) THEN
         guess_input = ref%guess_cli%input_file
         IF (LEN_TRIM(guess_input) == 0) guess_input = gen_file
         etemp_guess_str = ""
         IF (ref%guess_cli%electronic_temperature_guess > 0.0_dp) THEN
            etemp_guess = cp_unit_from_cp2k(ref%guess_cli%electronic_temperature_guess, "K")
            WRITE (etemp_guess_val_str, "(ES16.8)") etemp_guess
            etemp_guess_str = " --etemp-guess "//TRIM(ADJUSTL(etemp_guess_val_str))
         END IF
         efield_str = ""
         IF (ref%guess_cli%efield_active) THEN
            WRITE (efield_x_str, "(ES16.8)") ref%guess_cli%efield(1)
            WRITE (efield_y_str, "(ES16.8)") ref%guess_cli%efield(2)
            WRITE (efield_z_str, "(ES16.8)") ref%guess_cli%efield(3)
            efield_str = " --efield "//TRIM(ADJUSTL(efield_x_str))//","// &
                         TRIM(ADJUSTL(efield_y_str))//","//TRIM(ADJUSTL(efield_z_str))
         END IF
         grad_str = ""
         IF (ref%guess_cli%grad) grad_str = " --grad"
         json_str = ""
         IF (LEN_TRIM(ref%guess_cli%json_file) > 0) &
            json_str = " --json "//TRIM(tb_shell_quote(ref%guess_cli%json_file))
         log_file = TRIM(file_base)//".guess.log"
         command = TRIM(tb_shell_quote(ref%program_name))// &
                   " guess --charge "//TRIM(ADJUSTL(charge_str))// &
                   " --spin "//TRIM(ADJUSTL(spin_str))// &
                   " --method "//TRIM(tb_reference_guess_name(ref%guess_cli%method))// &
                   " --solver "//TRIM(tb_reference_solver_name(ref%guess_cli%solver))// &
                   TRIM(etemp_guess_str)// &
                   TRIM(efield_str)// &
                   TRIM(grad_str)// &
                   TRIM(json_str)// &
                   TRIM(verbosity_str)//" --input "//TRIM(tb_shell_quote(ref%guess_cli%input_format))// &
                   " "//TRIM(tb_shell_quote(guess_input))// &
                   " > "//TRIM(tb_shell_quote(log_file))//" 2>&1"
         CALL tb_reference_cli_execute(ref, "guess", command, log_file, iounit)
         IF (.NOT. ref%keep_files .AND. LEN_TRIM(ref%guess_cli%json_file) > 0) &
            CALL tb_delete_file(ref%guess_cli%json_file)
      END IF

      IF (ref%param_cli%enabled) THEN
         method_str = ""
         IF (ref%param_cli%method_explicit .OR. LEN_TRIM(ref%param_cli%input_file) == 0) &
            method_str = " --method "// &
                         TRIM(tb_reference_method_name(MERGE(ref%param_cli%method, &
                                                             dft_control%qs_control%xtb_control%tblite_method, &
                                                             ref%param_cli%method_explicit)))
         output_str = ""
         IF (LEN_TRIM(ref%param_cli%output_file) > 0) &
            output_str = " --output "//TRIM(tb_shell_quote(ref%param_cli%output_file))
         log_file = TRIM(file_base)//".param.log"
         command = TRIM(tb_shell_quote(ref%program_name))//" param"// &
                   TRIM(method_str)// &
                   TRIM(output_str)
         IF (LEN_TRIM(ref%param_cli%input_file) > 0) &
            command = TRIM(command)//" "//TRIM(tb_shell_quote(ref%param_cli%input_file))
         command = TRIM(command)//" > "//TRIM(tb_shell_quote(log_file))//" 2>&1"
         CALL tb_reference_cli_execute(ref, "param", command, log_file, iounit)
         IF (.NOT. ref%keep_files .AND. LEN_TRIM(ref%param_cli%output_file) > 0) &
            CALL tb_delete_file(ref%param_cli%output_file)
      END IF

      IF (ref%fit_cli%enabled) THEN
         dry_run_str = ""
         IF (ref%fit_cli%dry_run) dry_run_str = " --dry-run"
         copy_str = ""
         IF (LEN_TRIM(ref%fit_cli%copy_file) > 0) &
            copy_str = " --copy "//TRIM(tb_shell_quote(ref%fit_cli%copy_file))
         log_file = TRIM(file_base)//".fit.log"
         command = TRIM(tb_shell_quote(ref%program_name))//" fit"// &
                   TRIM(dry_run_str)// &
                   TRIM(copy_str)// &
                   TRIM(verbosity_str)//" "//TRIM(tb_shell_quote(ref%fit_cli%param_file))// &
                   " "//TRIM(tb_shell_quote(ref%fit_cli%input_file))// &
                   " > "//TRIM(tb_shell_quote(log_file))//" 2>&1"
         CALL tb_reference_cli_execute(ref, "fit", command, log_file, iounit)
         IF (.NOT. ref%keep_files .AND. LEN_TRIM(ref%fit_cli%copy_file) > 0) &
            CALL tb_delete_file(ref%fit_cli%copy_file)
      END IF

      IF (ref%tagdiff_cli%enabled) THEN
         method_str = ""
         IF (ref%tagdiff_cli%fit) method_str = " --fit"
         log_file = TRIM(file_base)//".tagdiff.log"
         command = TRIM(tb_shell_quote(ref%program_name))//" tagdiff"// &
                   TRIM(method_str)//" "//TRIM(tb_shell_quote(ref%tagdiff_cli%actual_file))// &
                   " "//TRIM(tb_shell_quote(ref%tagdiff_cli%reference_file))// &
                   " > "//TRIM(tb_shell_quote(log_file))//" 2>&1"
         CALL tb_reference_cli_execute(ref, "tagdiff", command, log_file, iounit)
      END IF

   END SUBROUTINE tb_reference_cli_aux_commands

! **************************************************************************************************
!> \brief Execute one native tblite REFERENCE_CLI auxiliary command.
!> \param ref ...
!> \param label ...
!> \param command ...
!> \param log_file ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE tb_reference_cli_execute(ref, label, command, log_file, iounit)

      TYPE(xtb_reference_cli_type), INTENT(IN)           :: ref
      CHARACTER(LEN=*), INTENT(IN)                       :: label, command, log_file
      INTEGER, INTENT(IN)                                :: iounit

      INTEGER                                            :: cmdstat, exitstat

      cmdstat = 0
      exitstat = 0
      CALL execute_command_line(TRIM(command), exitstat=exitstat, cmdstat=cmdstat)
      IF (cmdstat /= 0 .OR. exitstat /= 0) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A,A)") "tblite reference CLI auxiliary command failed: ", &
            TRIM(label)
         WRITE (UNIT=iounit, FMT="(T2,A,A)") "Command: ", TRIM(command)
         WRITE (UNIT=iounit, FMT="(T2,A,I0,T32,A,I0)") "cmdstat:", cmdstat, "exitstat:", exitstat
         IF (ref%stop_on_error) CPABORT("tblite reference CLI auxiliary command failed")
      ELSE
         WRITE (UNIT=iounit, FMT="(/,T2,A,A)") "tblite reference CLI auxiliary command completed: ", &
            TRIM(label)
         WRITE (UNIT=iounit, FMT="(T2,A,A)") "Log file: ", TRIM(log_file)
      END IF
      IF (.NOT. ref%keep_files) CALL tb_delete_file(log_file)

   END SUBROUTINE tb_reference_cli_execute

! **************************************************************************************************
!> \brief Join directory and filename.
!> \param directory ...
!> \param filename ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_join_path(directory, filename) RESULT(path)
      CHARACTER(LEN=*), INTENT(IN)                       :: directory, filename
      CHARACTER(LEN=default_path_length)                 :: path

      IF (LEN_TRIM(directory) == 0 .OR. TRIM(directory) == ".") THEN
         path = TRIM(filename)
      ELSE IF (directory(LEN_TRIM(directory):LEN_TRIM(directory)) == "/") THEN
         path = TRIM(directory)//TRIM(filename)
      ELSE
         path = TRIM(directory)//"/"//TRIM(filename)
      END IF

   END FUNCTION tb_join_path

! **************************************************************************************************
!> \brief Shell-quote a filename or executable path.
!> \param text ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_shell_quote(text) RESULT(quoted)
      CHARACTER(LEN=*), INTENT(IN)                       :: text
      CHARACTER(LEN=4*default_path_length)               :: quoted

      INTEGER                                            :: i

      quoted = "'"
      DO i = 1, LEN_TRIM(text)
         IF (text(i:i) == "'") THEN
            quoted = TRIM(quoted)//"'\\''"
         ELSE
            quoted = TRIM(quoted)//text(i:i)
         END IF
      END DO
      quoted = TRIM(quoted)//"'"

   END FUNCTION tb_shell_quote

! **************************************************************************************************
!> \brief Write current CP2K geometry as DFTB+ gen format for native tblite.
!> \param qs_env ...
!> \param filename ...
! **************************************************************************************************
   SUBROUTINE tb_write_reference_gen(qs_env, filename)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(LEN=*), INTENT(IN)                       :: filename

      CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:)        :: symbols, unique_symbols
      INTEGER                                            :: iatom, ikind, ios, natom, nuniq, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: species
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: to_angstrom
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env=qs_env, cell=cell, particle_set=particle_set)

      natom = SIZE(particle_set)
      to_angstrom = cp_unit_from_cp2k(1.0_dp, "angstrom")
      ALLOCATE (symbols(natom), unique_symbols(natom), species(natom))
      nuniq = 0
      DO iatom = 1, natom
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, element_symbol=symbols(iatom))
         found = .FALSE.
         DO ikind = 1, nuniq
            IF (TRIM(unique_symbols(ikind)) == TRIM(symbols(iatom))) THEN
               found = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT. found) THEN
            nuniq = nuniq + 1
            unique_symbols(nuniq) = symbols(iatom)
            ikind = nuniq
         END IF
         species(iatom) = ikind
      END DO

      OPEN (NEWUNIT=unit_nr, FILE=TRIM(filename), STATUS="REPLACE", ACTION="WRITE", &
            FORM="FORMATTED", IOSTAT=ios)
      IF (ios /= 0) CPABORT("Could not open tblite reference CLI geometry file")

      CALL get_cell(cell=cell, periodic=periodic)
      IF (ANY(periodic == 1)) THEN
         WRITE (UNIT=unit_nr, FMT="(I0,1X,A)") natom, "S"
      ELSE
         WRITE (UNIT=unit_nr, FMT="(I0,1X,A)") natom, "C"
      END IF
      WRITE (UNIT=unit_nr, FMT="(*(A,1X))") (TRIM(unique_symbols(ikind)), ikind=1, nuniq)
      DO iatom = 1, natom
         WRITE (UNIT=unit_nr, FMT="(I0,1X,I0,3(1X,ES24.16))") &
            iatom, species(iatom), particle_set(iatom)%r(:)*to_angstrom
      END DO
      IF (ANY(periodic == 1)) THEN
         WRITE (UNIT=unit_nr, FMT="(3(1X,ES24.16))") 0.0_dp, 0.0_dp, 0.0_dp
         DO ikind = 1, 3
            WRITE (UNIT=unit_nr, FMT="(3(1X,ES24.16))") cell%hmat(:, ikind)*to_angstrom
         END DO
      END IF
      CLOSE (unit_nr)

      DEALLOCATE (symbols, unique_symbols, species)

   END SUBROUTINE tb_write_reference_gen

! **************************************************************************************************
!> \brief Read native tblite gradient file.
!> \param filename ...
!> \param natom ...
!> \param energy ...
!> \param gradient ...
!> \param virial ...
!> \param have_energy ...
!> \param have_gradient ...
!> \param have_virial ...
! **************************************************************************************************
   SUBROUTINE tb_read_reference_grad(filename, natom, energy, gradient, virial, &
                                     have_energy, have_gradient, have_virial)

      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      INTEGER, INTENT(IN)                                :: natom
      REAL(KIND=dp), INTENT(OUT)                         :: energy
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: gradient, virial
      LOGICAL, INTENT(OUT)                               :: have_energy, have_gradient, have_virial

      CHARACTER(LEN=1024)                                :: line
      INTEGER                                            :: ios, nread, unit_nr
      LOGICAL                                            :: exists
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: values

      have_energy = .FALSE.
      have_gradient = .FALSE.
      have_virial = .FALSE.
      energy = 0.0_dp
      gradient = 0.0_dp
      virial = 0.0_dp

      INQUIRE (FILE=TRIM(filename), EXIST=exists)
      IF (.NOT. exists) RETURN

      OPEN (NEWUNIT=unit_nr, FILE=TRIM(filename), STATUS="OLD", ACTION="READ", &
            FORM="FORMATTED", IOSTAT=ios)
      IF (ios /= 0) RETURN

      DO
         READ (UNIT=unit_nr, FMT="(A)", IOSTAT=ios) line
         IF (ios /= 0) EXIT
         IF (INDEX(line, "energy             :real:0:") > 0) THEN
            READ (UNIT=unit_nr, FMT="(A)", IOSTAT=ios) line
            IF (ios == 0) THEN
               READ (line, *, IOSTAT=ios) energy
               have_energy = ios == 0
            END IF
         ELSE IF (INDEX(line, "gradient           :real:2:3,") > 0) THEN
            ALLOCATE (values(3*natom))
            CALL tb_read_real_values(unit_nr, values, nread)
            IF (nread == 3*natom) THEN
               CALL tb_values_to_matrix(values, gradient)
               have_gradient = .TRUE.
            END IF
            DEALLOCATE (values)
         ELSE IF (INDEX(line, "virial             :real:2:3,3") > 0) THEN
            ALLOCATE (values(9))
            CALL tb_read_real_values(unit_nr, values, nread)
            IF (nread == 9) THEN
               CALL tb_values_to_matrix(values, virial)
               have_virial = .TRUE.
            END IF
            DEALLOCATE (values)
         END IF
      END DO
      CLOSE (unit_nr)

   END SUBROUTINE tb_read_reference_grad

! **************************************************************************************************
!> \brief Read a fixed number of real values from following lines.
!> \param unit_nr ...
!> \param values ...
!> \param nread ...
! **************************************************************************************************
   SUBROUTINE tb_read_real_values(unit_nr, values, nread)

      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: values
      INTEGER, INTENT(OUT)                               :: nread

      CHARACTER(LEN=1024)                                :: line
      INTEGER                                            :: ios

      nread = 0
      DO WHILE (nread < SIZE(values))
         READ (UNIT=unit_nr, FMT="(A)", IOSTAT=ios) line
         IF (ios /= 0) EXIT
         CALL tb_parse_real_line(line, values, nread)
      END DO

   END SUBROUTINE tb_read_real_values

! **************************************************************************************************
!> \brief Parse real values from one text line.
!> \param line ...
!> \param values ...
!> \param nread ...
! **************************************************************************************************
   SUBROUTINE tb_parse_real_line(line, values, nread)

      CHARACTER(LEN=*), INTENT(IN)                       :: line
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: values
      INTEGER, INTENT(INOUT)                             :: nread

      CHARACTER(LEN=128)                                 :: token
      INTEGER                                            :: first, ios, last, pos

      pos = 1
      DO WHILE (pos <= LEN_TRIM(line) .AND. nread < SIZE(values))
         DO WHILE (pos <= LEN_TRIM(line) .AND. INDEX(" ,[]", line(pos:pos)) > 0)
            pos = pos + 1
         END DO
         IF (pos > LEN_TRIM(line)) EXIT
         first = pos
         DO WHILE (pos <= LEN_TRIM(line) .AND. INDEX(" ,[]", line(pos:pos)) == 0)
            pos = pos + 1
         END DO
         last = pos - 1
         token = line(first:last)
         READ (token, *, IOSTAT=ios) values(nread + 1)
         IF (ios == 0) nread = nread + 1
      END DO

   END SUBROUTINE tb_parse_real_line

! **************************************************************************************************
!> \brief Convert flat native tblite values to CP2K atom-major matrix layout.
!> \param values ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE tb_values_to_matrix(values, matrix)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: values
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: matrix

      INTEGER                                            :: i, j, n

      n = 0
      DO j = 1, SIZE(matrix, 2)
         DO i = 1, SIZE(matrix, 1)
            n = n + 1
            matrix(i, j) = values(n)
         END DO
      END DO

   END SUBROUTINE tb_values_to_matrix

! **************************************************************************************************
!> \brief Delete temporary files unless requested otherwise.
!> \param ref ...
!> \param gen_file ...
!> \param grad_file ...
!> \param json_file ...
!> \param log_file ...
!> \param post_processing_output_file ...
! **************************************************************************************************
   SUBROUTINE tb_reference_cleanup(ref, gen_file, grad_file, json_file, log_file, post_processing_output_file)

      TYPE(xtb_reference_cli_type), INTENT(IN)           :: ref
      CHARACTER(LEN=*), INTENT(IN)                       :: gen_file, grad_file, json_file, &
                                                            log_file, post_processing_output_file

      IF (ref%keep_files) RETURN
      CALL tb_delete_file(gen_file)
      CALL tb_delete_file(grad_file)
      CALL tb_delete_file(json_file)
      CALL tb_delete_file(log_file)
      IF (LEN_TRIM(post_processing_output_file) > 0) &
         CALL tb_delete_file(post_processing_output_file)

   END SUBROUTINE tb_reference_cleanup

! **************************************************************************************************
!> \brief Delete a file if it exists.
!> \param filename ...
! **************************************************************************************************
   SUBROUTINE tb_delete_file(filename)

      CHARACTER(LEN=*), INTENT(IN)                       :: filename

      INTEGER                                            :: ios, unit_nr
      LOGICAL                                            :: exists

      INQUIRE (FILE=TRIM(filename), EXIST=exists)
      IF (.NOT. exists) RETURN
      OPEN (NEWUNIT=unit_nr, FILE=TRIM(filename), STATUS="OLD", IOSTAT=ios)
      IF (ios == 0) CLOSE (unit_nr, STATUS="DELETE")

   END SUBROUTINE tb_delete_file

! **************************************************************************************************
!> \brief Dump cumulative tblite virial pieces for local debugging.
!> \param label ...
!> \param sigma ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_dump_sigma_component(label, sigma, para_env)

      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sigma
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=default_path_length)                 :: dump_file
      INTEGER                                            :: dump_status, dump_unit, i

      dump_status = 1
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_SIGMA_COMPONENT_DUMP", dump_file, STATUS=dump_status)
#endif
      IF (dump_status /= 0) RETURN

      OPEN (NEWUNIT=dump_unit, FILE=TRIM(dump_file), STATUS="UNKNOWN", &
            POSITION="APPEND", ACTION="WRITE")
      WRITE (dump_unit, "(A)") TRIM(label)
      DO i = 1, 3
         WRITE (dump_unit, "(3(1X,ES24.16))") sigma(i, :)/para_env%num_pe
      END DO
      CLOSE (dump_unit)

   END SUBROUTINE tb_dump_sigma_component

! **************************************************************************************************
!> \brief save stress tensor
!> \param qs_env ...
!> \param tb ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_add_stress(qs_env, tb, para_env)

      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tblite_type)                                  :: tb
      TYPE(mp_para_env_type)                             :: para_env

      CHARACTER(LEN=default_path_length)                 :: dump_file
      INTEGER                                            :: dump_status, dump_unit, i
      INTEGER, DIMENSION(3)                              :: periodic
      TYPE(cell_type), POINTER                           :: cell
      TYPE(virial_type), POINTER                         :: virial

      NULLIFY (virial, cell)
      CALL get_qs_env(qs_env=qs_env, virial=virial, cell=cell)
      CALL get_cell(cell=cell, periodic=periodic)

      IF (ALL(periodic == 0)) THEN
         CALL cp_warn(__LOCATION__, &
                      "tblite stress tensor requested for an isolated system. "// &
                      "The reported virial is useful for finite-difference checks, "// &
                      "but it is not a physically meaningful bulk stress for an isolated molecule.")
      END IF

      dump_status = 1
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_VIRIAL_DUMP", dump_file, STATUS=dump_status)
#endif
      IF (dump_status == 0) THEN
         OPEN (NEWUNIT=dump_unit, FILE=TRIM(dump_file), STATUS="UNKNOWN", &
               POSITION="APPEND", ACTION="WRITE")
         WRITE (dump_unit, "(A)") "sigma"
         DO i = 1, 3
            WRITE (dump_unit, "(3(1X,ES24.16))") tb%sigma(i, :)/para_env%num_pe
         END DO
         CLOSE (dump_unit)
      END IF

      virial%pv_virial = virial%pv_virial - tb%sigma/para_env%num_pe

   END SUBROUTINE tb_add_stress

! **************************************************************************************************
!> \brief add contrib. to gradient
!> \param grad ...
!> \param deriv ...
!> \param dE ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE tb_add_grad(grad, deriv, dE, natom)

      REAL(KIND=dp), DIMENSION(:, :)                     :: grad
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: deriv
      REAL(KIND=dp), DIMENSION(:)                        :: dE
      INTEGER                                            :: natom

      INTEGER                                            :: i, j

      DO i = 1, natom
         DO j = 1, natom
            grad(:, i) = grad(:, i) + deriv(:, i, j)*dE(j)
         END DO
      END DO

   END SUBROUTINE tb_add_grad

! **************************************************************************************************
!> \brief add contrib. to sigma
!> \param sig ...
!> \param deriv ...
!> \param dE ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE tb_add_sig(sig, deriv, dE, natom)

      REAL(KIND=dp), DIMENSION(:, :)                     :: sig
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: deriv
      REAL(KIND=dp), DIMENSION(:)                        :: dE
      INTEGER                                            :: natom

      INTEGER                                            :: i, j

      DO i = 1, 3
         DO j = 1, natom
            sig(:, i) = sig(:, i) + deriv(:, i, j)*dE(j)
         END DO
      END DO

   END SUBROUTINE tb_add_sig

END MODULE tblite_interface
