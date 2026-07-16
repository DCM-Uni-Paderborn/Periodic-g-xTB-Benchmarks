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

module test_exchange
   use mctc_data_vdwrad, only : get_vdw_rad
   use mctc_env, only : wp
   use mctc_io_convert, only : autoaa
   use mctc_io_constants, only : pi
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use mstore, only : get_structure
   use multicharge, only : get_eeqbc_charges
   use tblite_basis_cache, only : basis_cache
   use tblite_basis_qvszp, only : qvszp_basis_type, qvszp_cgto_type, &
      & new_qvszp_cgto, new_qvszp_basis
   use tblite_basis_type, only : basis_type, cgto_container, cgto_type
   use tblite_blas, only : gemv
   use tblite_container_cache, only : container_cache
   use tblite_cp2k_compat, only : cp2k_exchange_kmesh, &
      & cp2k_exchange_kmesh_gradient, cp2k_exchange_kpoint_gradient, &
      & cp2k_prepare_exchange, cp2k_exchange_stream_type, &
      & cp2k_exchange_stream_begin, cp2k_exchange_stream_push, &
      & cp2k_exchange_stream_apply, cp2k_exchange_stream_get, &
      & cp2k_exchange_stream_reverse_apply, &
      & cp2k_exchange_stream_reverse_get, &
      & cp2k_exchange_stream_end, cp2k_exchange_stream_has_full_mesh_storage, &
      & cp2k_exchange_stream_reduced, cp2k_exchange_stream_oracle
   use tblite_cutoff, only : get_lattice_points
   use tblite_data_onecxints, only : get_onecxints
   use tblite_exchange, only : exchange_type, new_exchange_fock, exchange_fock, &
      & exchange_bvk_kernel, exchange_cache
   use tblite_integral_overlap, only : get_overlap
   use tblite_scf, only: new_potential, potential_type
   use tblite_utils_average, only : average_type, new_average, average_id
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: collect_exchange

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine exchange_maker(exchange, mol, bas)
         import :: exchange_type, structure_type, basis_type
         class(exchange_type), allocatable, intent(out) :: exchange
         type(structure_type), intent(in) :: mol
         !> Basis set information
         class(basis_type), intent(in) :: bas
      end subroutine exchange_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_exchange(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gamma_fock_h2", test_gamma_fock_h2), &
      new_unittest("gamma_fock_lih", test_gamma_fock_lih), &
      new_unittest("gamma_fock_no", test_gamma_fock_no), &
      new_unittest("gamma_fock_s2", test_gamma_fock_s2), &
      new_unittest("gamma_fock_cecl3", test_gamma_fock_cecl3), &
      new_unittest("energy_fock_h2", test_e_fock_h2), &
      new_unittest("energy_fock_lih", test_e_fock_lih), &
      new_unittest("energy_fock_no", test_e_fock_no), &
      new_unittest("energy_fock_n2", test_e_fock_n2), &
      new_unittest("energy_fock_h2o", test_e_fock_h2o), &
      new_unittest("energy_fock_s2", test_e_fock_s2), &
      new_unittest("energy_fock_sih4", test_e_fock_sih4), &
      new_unittest("energy_fock_cecl3", test_e_fock_cecl3), &
      new_unittest("energy_fock_ce2", test_e_fock_ce2), &
      new_unittest("potential_fock_h2", test_p_fock_h2), &
      new_unittest("potential_fock_lih", test_p_fock_lih), &
      new_unittest("potential_fock_no", test_p_fock_no), &
      new_unittest("potential_fock_s2", test_p_fock_s2), &
      new_unittest("potential_fock_cecl3", test_p_fock_cecl3), &
      new_unittest("potential_fock_ce2", test_p_fock_ce2), &
      new_unittest("gradient_op_fock_h2", test_op_g_fock_h2), &
      new_unittest("gradient_op_fock_lih", test_op_g_fock_lih), &
      new_unittest("gradient_op_fock_no", test_op_g_fock_no), &
      new_unittest("gradient_op_fock_s2", test_op_g_fock_s2), &
      new_unittest("gradient_op_fock_cecl3", test_op_g_fock_cecl3), &
      new_unittest("gradient_op_fock_pyrazine", test_op_g_fock_pyrazine), &
      new_unittest("gradient_op_fock_oxacb", test_op_g_fock_oxacb), &
      new_unittest("gradient_op_fock_co2", test_op_g_fock_co2), &
      new_unittest("bvk_exchange_supercell", test_bvk_exchange_kernel), &
      new_unittest("kpoint_fock_complex", test_kpoint_fock_complex) &
   ]

end subroutine collect_exchange


subroutine make_qvszp_basis(bas, mol, error)
   type(qvszp_basis_type), intent(out) :: bas
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error

   !> Parameter: Number of shells selected from the q-vSZP basis set
   integer, parameter :: pa_nshell(60) = [ &
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, & !1-20
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & !21-40
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 4, 4, 4] !41-60

   integer :: isp, izp, ish
   integer, allocatable :: nshell(:)
   type(cgto_container), allocatable :: cgto(:, :)
   type(cgto_container), allocatable :: cgto_h0(:, :)
   type(qvszp_cgto_type), allocatable :: cgto_qvszp
   real(wp) :: alpha(12), k0, k2, k3

   nshell = pa_nshell(mol%num)
   allocate(cgto(maxval(nshell), mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nshell(isp)
         allocate(cgto_qvszp)
         call new_qvszp_cgto(cgto_qvszp, izp, ish, .true., error)
         if (allocated(error)) return
         call move_alloc(cgto_qvszp, cgto(ish, isp)%raw)
      end do
   end do

   call new_qvszp_basis(bas, mol, nshell, cgto, error, accuracy=1.0_wp, &
      & cgto_h0=cgto_h0)
   if (allocated(error)) return

end subroutine make_qvszp_basis


!> Factory to create exchange object
subroutine make_exchange_gxtb(exchange, mol, bas)

   !> New exchange object
   class(exchange_type), allocatable, intent(out) :: exchange

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Basis set information
   class(basis_type), intent(in) :: bas

   !> Atomic Hubbard parameters or chemical hardness
   real(wp), parameter :: hubbard_parameter(60) = [&
      &  0.4725928800_wp,  0.9220339100_wp,  0.1745288800_wp,  0.2570073300_wp, & !1-4
      &  0.3394908600_wp,  0.4219541200_wp,  0.5043819300_wp,  0.5869186300_wp, & !5-8
      &  0.6693135100_wp,  0.7519160700_wp,  0.1796410500_wp,  0.2215727600_wp, & !9-12
      &  0.2634857800_wp,  0.3053964500_wp,  0.3473401400_wp,  0.3892472500_wp, & !13-16
      &  0.4311567000_wp,  0.4730826900_wp,  0.1710546900_wp,  0.2027624400_wp, & !17-20
      &  0.2100732200_wp,  0.2173964700_wp,  0.2247103900_wp,  0.2320150100_wp, & !21-24
      &  0.2393396900_wp,  0.2466563800_wp,  0.2539825500_wp,  0.2612886300_wp, & !25-28
      &  0.2685947600_wp,  0.2759256500_wp,  0.3076299900_wp,  0.3393158000_wp, & !29-32
      &  0.3723598500_wp,  0.4027354900_wp,  0.4344577600_wp,  0.4661170800_wp, & !33-36
      &  0.1558507900_wp,  0.1864932400_wp,  0.1935621000_wp,  0.2006331100_wp, & !37-40
      &  0.2077052200_wp,  0.2147725400_wp,  0.2218461400_wp,  0.2289187200_wp, & !41-44
      &  0.2359862100_wp,  0.2430561200_wp,  0.2501301800_wp,  0.2571993700_wp, & !45-48
      &  0.2878478000_wp,  0.3184867300_wp,  0.3491243100_wp,  0.3797659300_wp, & !49-52
      &  0.4104080800_wp,  0.4410577700_wp,  0.0501933200_wp,  0.0676257000_wp, & !53-56
      &  0.0850444500_wp,  0.1024773600_wp,  0.1199110500_wp,  0.1373277200_wp] !57-60

   real(wp), parameter :: fock_shell_hubbard(4, 60) = reshape([&
      &  3.3736927441_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !1
      &  5.7821690333_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !2
      &  1.3351178406_wp,  6.8872623116_wp,  0.0000000000_wp,  0.0000000000_wp, & !3
      &  2.4034872257_wp,  4.2333715067_wp,  0.0000000000_wp,  0.0000000000_wp, & !4
      &  1.4293672476_wp,  3.4588078363_wp,  0.0000000000_wp,  0.0000000000_wp, & !5
      &  2.5669958364_wp,  2.6678235131_wp,  0.0000000000_wp,  0.0000000000_wp, & !6
      &  4.8091413037_wp,  2.4087151818_wp,  0.0000000000_wp,  0.0000000000_wp, & !7
      &  5.1575887922_wp,  1.7604336086_wp,  0.0000000000_wp,  0.0000000000_wp, & !8
      &  5.6943326551_wp,  1.8352088241_wp,  0.0000000000_wp,  0.0000000000_wp, & !9
      &  4.5498436963_wp,  2.2916791150_wp,  0.0000000000_wp,  0.0000000000_wp, & !10
      &  3.7580992801_wp,  7.9514333677_wp,  0.0000000000_wp,  0.0000000000_wp, & !11
      &  2.0745094615_wp,  0.6401688716_wp,  0.0717621465_wp,  0.0000000000_wp, & !12
      &  2.6629478250_wp,  2.7525805658_wp,  0.0405648634_wp,  0.0000000000_wp, & !13
      &  4.3696365329_wp,  3.1617167429_wp,  0.0503738814_wp,  0.0000000000_wp, & !14
      &  5.6145647692_wp,  2.8978111700_wp,  0.0537795667_wp,  0.0000000000_wp, & !15
      &  3.8873165963_wp,  2.4254422115_wp,  0.0442026749_wp,  0.0000000000_wp, & !16
      &  3.9645174695_wp,  2.6867457351_wp,  1.8118397987_wp,  0.0000000000_wp, & !17
      &  4.0741332202_wp,  3.1424848876_wp,  0.1430915244_wp,  0.0000000000_wp, & !18
      &  1.3072532381_wp,  6.3378448661_wp,  0.0000000000_wp,  0.0000000000_wp, & !19
      &  1.2750980853_wp,  0.0635898847_wp,  0.2566121433_wp,  0.0000000000_wp, & !20
      &  3.2177234627_wp,  0.1959633878_wp,  5.7183716039_wp,  0.0000000000_wp, & !21
      &  4.1533395668_wp,  0.4152874467_wp,  6.2572212867_wp,  0.0000000000_wp, & !22
      &  4.1591909613_wp,  0.1390620807_wp,  5.9663550535_wp,  0.0000000000_wp, & !23
      &  2.7046153001_wp,  0.3415990565_wp,  5.7762402164_wp,  0.0000000000_wp, & !24
      &  3.6572739230_wp,  0.3129330436_wp,  5.8444982170_wp,  0.0000000000_wp, & !25
      &  2.2651998337_wp,  1.3653240301_wp,  6.0069435652_wp,  0.0000000000_wp, & !26
      &  2.5801471333_wp,  1.8964342295_wp,  5.4116718965_wp,  0.0000000000_wp, & !27
      &  2.5360729882_wp,  0.7118663442_wp,  5.7481530472_wp,  0.0000000000_wp, & !28
      &  2.9461034490_wp,  0.1230794473_wp,  5.8710029264_wp,  0.0000000000_wp, & !29
      &  2.6265625768_wp,  0.1448974804_wp,  0.0000000000_wp,  0.0000000000_wp, & !30
      &  4.0208890030_wp,  3.8443882747_wp,  0.0550326561_wp,  0.0000000000_wp, & !31
      &  3.9416820971_wp,  2.8623156888_wp,  0.0604375873_wp,  0.0000000000_wp, & !32
      &  3.8763117016_wp,  2.3413349488_wp,  0.1378292565_wp,  0.0000000000_wp, & !33
      &  3.8302236000_wp,  2.1240501008_wp,  0.0522070041_wp,  0.0000000000_wp, & !34
      &  2.6511323619_wp,  2.3146852635_wp,  0.0477540983_wp,  0.0000000000_wp, & !35
      &  4.0236949882_wp,  1.9479648586_wp,  0.0528299056_wp,  0.0000000000_wp, & !36
      &  3.8028426024_wp,  6.1510800760_wp,  0.0000000000_wp,  0.0000000000_wp, & !37
      &  1.5398833615_wp,  0.0469962124_wp,  0.4618960428_wp,  0.0000000000_wp, & !38
      &  2.6192646364_wp,  0.1774815842_wp,  5.7992744018_wp,  0.0000000000_wp, & !39
      &  2.9771698621_wp,  0.3273649537_wp,  6.1815149222_wp,  0.0000000000_wp, & !40
      &  4.0363248368_wp,  0.1815082362_wp,  5.8443575218_wp,  0.0000000000_wp, & !41
      &  3.0398497575_wp,  0.2478952306_wp,  4.5246140125_wp,  0.0000000000_wp, & !42
      &  3.1729738779_wp,  0.2323279988_wp,  4.4717438114_wp,  0.0000000000_wp, & !43
      &  3.3330011595_wp,  1.2564839962_wp,  5.6480626609_wp,  0.0000000000_wp, & !44
      &  3.6797088328_wp,  1.3554372375_wp,  4.6463903741_wp,  0.0000000000_wp, & !45
      &  2.3183681380_wp,  0.7409404169_wp,  3.8240765933_wp,  0.0000000000_wp, & !46
      &  4.1192819848_wp,  0.0967114383_wp,  4.4370977082_wp,  0.0000000000_wp, & !47
      &  3.0835232605_wp,  0.0438366605_wp,  0.0000000000_wp,  0.0000000000_wp, & !48
      &  4.7502780901_wp,  3.8043909165_wp,  0.0549850116_wp,  0.0000000000_wp, & !49
      &  4.1139475134_wp,  2.7119042380_wp,  0.0646933085_wp,  0.0000000000_wp, & !50
      &  4.9780075612_wp,  2.6973509435_wp,  0.0635386575_wp,  0.0000000000_wp, & !51
      &  3.1439645685_wp,  1.9915060640_wp,  0.0430816035_wp,  0.0000000000_wp, & !52
      &  3.0968244165_wp,  2.2532539703_wp,  0.0437247742_wp,  0.0000000000_wp, & !53
      &  3.4990871432_wp,  1.5345773312_wp,  0.1367881276_wp,  0.0000000000_wp, & !54
      &  3.0594065734_wp,  7.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !55
      &  2.7240224409_wp,  0.2030444630_wp,  0.4354451565_wp,  0.0000000000_wp, & !56
      &  3.7394472644_wp,  0.2280241628_wp, 13.4231022158_wp,  0.0000000000_wp, & !57
      &  3.9877249366_wp,  0.4186506926_wp,  3.1603362554_wp,  9.5757084786_wp, & !58
      &  4.1588690283_wp,  0.3985958915_wp,  2.1669384435_wp,  8.9609255762_wp, & !59
      &  4.2724917144_wp,  1.5603613289_wp,  2.6049529753_wp, 10.5049543896_wp],& !60
      & shape(fock_shell_hubbard))

   real(wp), parameter :: fock_avg_exp(4, 60) = reshape([&
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !1
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !2
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !3
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !4
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !5
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !6
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !7
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !8
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !9
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !10
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !11
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !12
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !13
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !14
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !15
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !16
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !17
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !18
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !19
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !20
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !21
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !22
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !23
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !24
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !25
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !26
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !27
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !28
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !29
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !30
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !31
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !32
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !33
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !34
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !35
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !36
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !37
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !38
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !39
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !40
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !41
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !42
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !43
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !44
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !45
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !46
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !47
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !48
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !49
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !50
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !51
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !52
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !53
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !54
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !55
      &  1.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp, & !56
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !57
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !58
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !59
      &  1.0000000000_wp,  2.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp],& !60
      & shape(fock_avg_exp))

   !> Parameter: Correlation radius for bond-order correlation correction to the Fock exchange
   real(wp), parameter :: fock_corr_rad(60) = [&
      & 14.3839259001_wp, 69.8731739578_wp,  6.5320635735_wp,  7.3163209686_wp, & !1-4
      & 17.1096864921_wp, 14.0786602534_wp,  4.8617350446_wp,  2.7996336557_wp, & !5-8
      &  0.7023326119_wp,  2.5626334368_wp, 14.0657663817_wp,  3.7434249280_wp, & !9-12
      &  6.0091669156_wp, 11.7789515108_wp, 18.2551554677_wp,  6.7406985032_wp, & !13-16
      &  5.3208032300_wp,  8.5000000000_wp,  9.0000000000_wp,  2.6999106381_wp, & !17-20
      &  4.2067112352_wp,  6.2106126930_wp,  9.4638003846_wp,  6.2838931204_wp, & !21-24
      &  4.0356555486_wp,  7.7605187678_wp,  2.1902457130_wp,  5.8666053908_wp, & !25-28
      &  5.8303790215_wp,  6.8955319781_wp, 12.2944956646_wp,  7.2446956658_wp, & !29-32
      & 11.5934354745_wp,  7.1345083204_wp,  4.2685057643_wp,  6.3145368912_wp, & !33-36
      & 11.0760022868_wp,  4.4062242213_wp,  3.2423160508_wp,  7.8875122192_wp, & !37-40
      &  7.5872400810_wp,  9.1209497744_wp,  4.3428939689_wp, 10.1832745738_wp, & !41-44
      &  5.7030273921_wp, 10.9046610938_wp, 11.1582685110_wp,  5.7040786241_wp, & !45-48
      & 11.7537065100_wp,  8.3351454926_wp, 12.9946420528_wp, 11.5449196095_wp, & !49-52
      &  4.0944483072_wp,  5.2170336023_wp, 15.0000000000_wp,  8.0046362470_wp, & !53-56
      &  5.8601951026_wp,  2.6802487350_wp,  4.2753129680_wp,  5.0680352245_wp] !57-60
      
   !> Parameter: Scaling factor for bond-order correlation correction to the Fock exchange
   real(wp), parameter :: fock_corr_scale(60) = [&
      & -0.0098632492_wp, -0.0017675886_wp, -0.0007595460_wp, -0.0001650446_wp, & !1-4
      & -0.0100209651_wp, -0.0126422276_wp, -0.0342892622_wp, -0.0178181598_wp, & !5-8
      & -0.0662165049_wp, -0.0689017240_wp, -0.0518985878_wp, -0.0217939514_wp, & !9-12
      & -0.0074794756_wp, -0.0130603938_wp, -0.0114161245_wp, -0.0239163432_wp, & !13-16
      & -0.0073912488_wp, -0.0392779244_wp, -0.0121215583_wp, -0.0124350881_wp, & !17-20
      & -0.0015511412_wp, -0.0033761025_wp, -0.0163814942_wp, -0.0200187053_wp, & !21-24
      & -0.0150071461_wp, -0.0010523703_wp, -0.0040408510_wp, -0.0006406616_wp, & !25-28
      & -0.0176195804_wp, -0.0073902932_wp, -0.0161133134_wp, -0.0012004437_wp, & !29-32
      & -0.0057631315_wp, -0.0086071526_wp, -0.0062089443_wp, -0.0184032543_wp, & !33-36
      & -0.0349279664_wp, -0.0020233113_wp, -0.0075322731_wp, -0.0004720921_wp, & !37-40
      & -0.0015572573_wp, -0.0074100261_wp, -0.0066068332_wp, -0.0006283098_wp, & !41-44
      & -0.0084164850_wp, -0.0030546848_wp, -0.0039673973_wp, -0.0104996025_wp, & !45-48
      & -0.0031252618_wp, -0.0028410123_wp, -0.0114280369_wp, -0.0044901228_wp, & !49-52
      & -0.0080558853_wp, -0.0243832808_wp, -0.0416556059_wp, -0.0071439456_wp, & !53-56
      & -0.0290053165_wp, -0.0100176448_wp, -0.0029583603_wp, -0.0105390843_wp] !57-60

   !> Parameter: van-der-Waals radii scaling
   real(wp), parameter :: rvdw_scale(60) = [&
      &  1.0000000000_wp,  1.0000000000_wp,  1.0600000000_wp,  1.0000000000_wp, & !1-4
      &  1.0200000000_wp,  1.0000000000_wp,  1.0500000000_wp,  1.0800000000_wp, & !5-8
      &  1.1500000000_wp,  0.8000000000_wp,  1.2000000000_wp,  1.1500000000_wp, & !9-12
      &  1.0700000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !13-16
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !17-20
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !21-24
      &  1.0000000000_wp,  1.1000000000_wp,  1.1000000000_wp,  1.1000000000_wp, & !25-28
      &  1.0000000000_wp,  1.1500000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !29-32
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !33-36
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !37-40
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !41-44
      &  1.1000000000_wp,  1.1000000000_wp,  1.1000000000_wp,  1.1000000000_wp, & !45-48
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !49-52
      &  1.0000000000_wp,  1.0000000000_wp,  0.6500000000_wp,  1.0000000000_wp, & !53-56
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp] !57-60

   real(wp), parameter :: p_kq(0:3) = &
      & [1.1000000000_wp, 0.5500000000_wp, 0.2750000000_wp, 0.1375000000_wp]
   real(wp), parameter :: p_ondiag = 1.3826597204_wp
   real(wp), parameter :: p_offdiag_l(0:3) = &
      & [0.0500000000_wp, 2.5000000000_wp, 2.5000000000_wp, 2.5000000000_wp]
   real(wp), parameter :: p_hubbard_exp = -0.2960502355_wp
   real(wp), parameter :: p_hubbard_exp_r0 = 0.0349817377_wp * autoaa
   real(wp), parameter :: p_omega = 0.2000000000_wp
   real(wp), parameter :: p_a0 = 0.1500000000_wp
   real(wp), parameter :: p_corr_exp = 1.0000000000_wp / autoaa

   integer :: izp, jzp, isp, jsp, ish, jsh, il, jl
   real(wp), allocatable :: hardness_fx(:, :), onecxints(:, :, :), offdiag(:, :), kq(:, :)
   real(wp), allocatable :: avg_exponents(:, :), rad(:, :)
   real(wp), allocatable :: corr_scale(:), corr_rad(:)
   type(exchange_fock), allocatable :: fock
   type(average_type), allocatable :: hubbard_average
   type(average_type), allocatable :: offdiag_average
   type(average_type), allocatable :: corr_scale_average
   type(average_type), allocatable :: corr_rad_average
   type(average_type), allocatable :: rvdw_average
   
   ! Obtain radii for radius-dependent hubbard scaling
   allocate(rvdw_average)
   call new_average(rvdw_average, average_id%arithmetic)
   allocate(rad(mol%nid, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, isp
         jzp = mol%num(jsp)
         rad(jsp, isp) = get_vdw_rad(jzp, izp) &
            & * rvdw_average%value(rvdw_scale(izp), rvdw_scale(jzp))
         rad(isp, jsp) = rad(jsp, isp)
      end do
   end do

   ! Obtain the shell Hubbard parameters
   allocate(hardness_fx(maxval(bas%nsh_id), mol%nid))
   hardness_fx(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         hardness_fx(ish, isp) = hubbard_parameter(izp) &
            & * fock_shell_hubbard(ish, izp)
      end do
   end do

   ! Set averager for hubbard and correlation parameters
   allocate(hubbard_average, offdiag_average, corr_scale_average, corr_rad_average)
   call new_average(hubbard_average, average_id%general)
   call new_average(offdiag_average, average_id%geometric)
   call new_average(corr_scale_average, average_id%arithmetic)
   call new_average(corr_rad_average, average_id%geometric)

   ! Obtain the averaging exponents for the shell Hubbard parameters
   allocate(avg_exponents(maxval(bas%nsh_id), mol%nid))
   avg_exponents(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         avg_exponents(ish, isp) = fock_avg_exp(ish, izp)
      end do
   end do

   ! Obtain the one-center exchange integrals
   allocate(onecxints(maxval(bas%nsh_id), maxval(bas%nsh_id), mol%nid))
   onecxints(:, :, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         il = bas%cgto(ish, isp)%raw%ang
         do jsh = 1, bas%nsh_id(isp)
            jl = bas%cgto(jsh, isp)%raw%ang
            onecxints(jsh, ish, isp) = get_onecxints(jl, il, izp)
         end do
      end do
   end do

   ! Obtain the off-diagaonal scaling exchange factors
   allocate(offdiag(maxval(bas%nsh_id), mol%nid))
   offdiag(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         il = bas%cgto(ish, isp)%raw%ang
         offdiag(ish, isp) = p_offdiag_l(il)
      end do
   end do

   ! Obtain the charge scaling of the exchange fraction
   allocate(kq(maxval(bas%nsh_id), mol%nid))
   kq(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         il = bas%cgto(ish, isp)%raw%ang
         kq(ish, isp) = p_kq(il)
      end do
   end do

   ! Obtain the correlation scaling factor and radii
   corr_scale = fock_corr_scale(mol%num)
   corr_rad = fock_corr_rad(mol%num)

   allocate(fock)
   call new_exchange_fock(fock, mol, bas, hardness_fx, hubbard_average, &
      & avg_exponents, p_ondiag, offdiag, offdiag_average, p_hubbard_exp, &
      & p_hubbard_exp_r0, rad, 1.0_wp, onecxints, kq, corr_scale, &
      & corr_scale_average, p_corr_exp, corr_rad, corr_rad_average, p_a0, &
      & p_omega, 1.0_wp-p_a0)
   call move_alloc(fock, exchange)

end subroutine make_exchange_gxtb


subroutine test_gamma_generic(error, mol, make_exchange, qsh, ref_g_mulliken, &
   & ref_g_onsfx, ref_dgdq_onsfx, ref_g_onsri, ref_dgdq_onsri, ref_g_bocorr, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Factory to create new exchange objects
   procedure(exchange_maker) :: make_exchange

   !> Shell-resolved charges
   real(wp), intent(in) :: qsh(:)

   !> Reference matrix for Mulliken to check against
   real(wp), intent(in) :: ref_g_mulliken(:, :)
   
   !> Compressed reference matrix for onsite to check against
   real(wp), intent(in) :: ref_g_onsfx(:, :, :)

   !> Compressed reference matrix for charge derivative of onsite to check against
   real(wp), intent(in) :: ref_dgdq_onsfx(:, :, :)

   !> Compressed reference matrix for rotational invariance to check against
   real(wp), intent(in) :: ref_g_onsri(:, :)

   !> Compressed reference matrix for charge derivative of rotational invariance to check against
   real(wp), intent(in) :: ref_dgdq_onsri(:, :)

   !> Reference matrix for bond-order correlation to check against
   real(wp), intent(in) :: ref_g_bocorr(:, :)

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   class(basis_type), allocatable :: bas
   type(basis_cache) :: bcache
   class(exchange_type), allocatable :: exchange
   type(container_cache) :: cache
   type(exchange_cache), pointer :: ptr
   type(wavefunction_type) :: wfn, wfn_aux
   class(qvszp_basis_type), allocatable :: qvszp_bas
   real(wp) :: thr_
   real(wp), allocatable :: lattr(:, :)
   real(wp), allocatable :: overlap(:,:)
   integer :: ii, jj

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup basis
   allocate(qvszp_bas)
   call make_qvszp_basis(qvszp_bas, mol, error)
   call move_alloc(qvszp_bas, bas)
   if (allocated(error)) return

   ! Obtain EEQBC charges adapt the basis set
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   ! Setup the wavefunction with the provided charges
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 1, &
      & 0.0_wp, .false.)
   wfn%qsh(:, 1) = qsh

   ! Get lattice points
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   ! Setup and calculate exchange
   call make_exchange(exchange, mol, bas)

   call exchange%update(mol, cache)

   call view(cache, ptr)
   ! Calculate the onsite exchange matrix with charge dependence
   call exchange%get_onsite_Kmatrix(mol, wfn, ptr)

   ! Check Mulliken matrix
   if (any(abs(ptr%g_mulliken - ref_g_mulliken) > thr_)) then
      call test_failed(error, "Mulliken matrix does not match")
      write(*,*) "Reference:"
      print'(3es21.14)', ref_g_mulliken
      write(*,*) "Mulliken matrix:"
      print'(3es21.14)', ptr%g_mulliken
      write(*,*) "difference:"
      print'(3es21.14)', ptr%g_mulliken - ref_g_mulliken
   end if

   ! Check onsite matrix
   if (any(abs(ptr%g_onsfx - ref_g_onsfx) > thr_)) then
      call test_failed(error, "Onsite matrix does not match")
      write(*,*) "Reference:"
      print'(3es21.14)', ref_g_onsfx
      write(*,*) "Onsite matrix:"
      print'(3es21.14)', ptr%g_onsfx
      write(*,*) "difference:"
      print'(3es21.14)', ptr%g_onsfx - ref_g_onsfx
   end if

   ! Check charge-derivative of onsite matrix
   if (any(abs(ptr%dgdq_onsfx - ref_dgdq_onsfx) > thr_)) then
      call test_failed(error, "Charge-derivative of onsite matrix does not match")
      write(*,*) "Reference:"
      print'(3es21.14)', ref_dgdq_onsfx
      write(*,*) "Charge-derivative of onsite matrix:"
      print'(3es21.14)', ptr%dgdq_onsfx
      write(*,*) "difference:"
      print'(3es21.14)', ptr%dgdq_onsfx - ref_dgdq_onsfx
   end if

   ! Check rotational invariance matrix
   if (any(abs(ptr%g_onsri - ref_g_onsri) > thr_)) then
      call test_failed(error, "Rotational invariance matrix does not match")
      write(*,*) "Reference:"
      print'(3es21.14)', ref_g_onsri
      write(*,*) "Rotational invariance matrix:"
      print'(3es21.14)', ptr%g_onsri
      write(*,*) "difference:"
      print'(3es21.14)', ptr%g_onsri - ref_g_onsri
   end if

   ! Check charge-derivative of rotational invariance matrix
   if (any(abs(ptr%dgdq_onsri - ref_dgdq_onsri) > thr_)) then
      call test_failed(error, "Charge-derivative of rotational invariance matrix does not match")
      write(*,*) "Reference:"
      print'(3es21.14)', ref_dgdq_onsri
      write(*,*) "Charge-derivative rotational invariance matrix:"
      print'(3es21.14)', ptr%dgdq_onsri
      write(*,*) "difference:"
      print'(3es21.14)', ptr%dgdq_onsri - ref_dgdq_onsri
   end if
   
   ! Check bond-order correlation matrix
   if (any(abs(ptr%g_bocorr - ref_g_bocorr) > thr_)) then
      call test_failed(error, "Bond-order correlation matrix does not match")
      write(*,*) "Reference:"
      print'(3es21.14)', ref_g_bocorr
      write(*,*) "Bond-order correlation matrix:"
      print'(3es21.14)', ptr%g_bocorr
      write(*,*) "difference:"
      print'(3es21.14)', ptr%g_bocorr - ref_g_bocorr
   end if

end subroutine test_gamma_generic

!> Return reference to container cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(exchange_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(exchange_cache)
      ptr => target
   end select
end subroutine view


subroutine test_energy_generic(error, mol, density, qsh, make_exchange, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Density matrix for this structure
   real(wp), intent(in) :: density(:, :, :)

   !> Shell partial charges for this structure
   real(wp), intent(in) :: qsh(:)

   !> Factory to create new exchange objects
   procedure(exchange_maker) :: make_exchange

   !> Reference value to check against
   real(wp), intent(in) :: ref

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   class(basis_type), allocatable :: bas
   class(qvszp_basis_type), allocatable :: qvszp_bas
   type(basis_cache) :: bcache
   class(exchange_type), allocatable :: exchange
   type(potential_type) :: pot
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn, wfn_aux
   integer :: nspin
   real(wp) :: energy(mol%nat)
   real(wp), allocatable :: lattr(:, :)
   real(wp), allocatable :: overlap(:,:)
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in
   nspin = size(density, dim=3)

   ! Setup basis set
   allocate(qvszp_bas)
   call make_qvszp_basis(qvszp_bas, mol, error)
   call move_alloc(qvszp_bas, bas)
   if (allocated(error)) return

   call check(error, bas%nao, size(density, 1))
   if (allocated(error)) return

   ! Obtain EEQBC charges adapt the basis set
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   ! Get lattice points
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   ! Calculate the overlap integrals
   allocate(overlap(bas%nao, bas%nao))
   call get_overlap(mol, lattr, cutoff, bas, bcache, overlap)

   ! Setup potential and wavefunction
   call new_potential(pot, mol, bas, nspin, .false.)
   call pot%reset()
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, nspin, &
      & 0.0_wp, .false.)
   wfn%density(:, :, :) = density
   wfn%qsh(:, 1) = qsh

   ! Setup and calculate exchange
   call make_exchange(exchange, mol, bas)

   energy = 0.0_wp
   call exchange%update(mol, cache)
   call exchange%get_energy_w_overlap(mol, cache, wfn, overlap, energy)

   call check(error, sum(energy), ref, thr=thr_)
   if (allocated(error)) then
      print*,ref, sum(energy), sum(energy) - ref
   end if

end subroutine test_energy_generic


subroutine test_numpot(error, mol, density, qsh, make_exchange, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Density matrix for this structure
   real(wp), intent(in) :: density(:, :, :)

   !> Shell-resolved charges for this structure
   real(wp), intent(in) :: qsh(:)

   !> Factory to create new exchange objects
   procedure(exchange_maker) :: make_exchange

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   class(basis_type), allocatable :: bas
   class(qvszp_basis_type), allocatable :: qvszp_bas
   type(basis_cache) :: bcache
   class(exchange_type), allocatable :: exchange
   type(potential_type) :: pot
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn, wfn_aux
   integer :: ish, nspin
   real(wp) :: energy(mol%nat)
   real(wp), allocatable :: lattr(:, :)
   real(wp), allocatable :: overlap(:,:)
   real(wp), allocatable :: numvsh(:)
   real(wp) :: er(mol%nat), el(mol%nat)
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in
   nspin = size(density, dim=3)

   ! Setup basis set
   allocate(qvszp_bas)
   call make_qvszp_basis(qvszp_bas, mol, error)
   call move_alloc(qvszp_bas, bas)
   if (allocated(error)) return

   call check(error, bas%nao, size(density, 1))
   if (allocated(error)) return

   ! Obtain EEQBC charges adapt the basis set
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   ! Get lattice points
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   ! Calculate the overlap integrals
   allocate(overlap(bas%nao, bas%nao))
   call get_overlap(mol, lattr, cutoff, bas, bcache, overlap)

   ! Setup potential and wavefunction
   call new_potential(pot, mol, bas, nspin, .false.)
   call pot%reset()
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, nspin, &
      & 0.0_wp, .false.)
   wfn%density(:, :, :) = density
   wfn%qsh(:, 1) = qsh

   ! Setup and calculate exchange
   call make_exchange(exchange, mol, bas)
   call exchange%update(mol, cache)
      
   ! Numerical atomic potential
   allocate(numvsh(bas%nsh), source=0.0_wp)
   do ish = 1, bas%nsh
      er = 0.0_wp
      el = 0.0_wp
      ! Right hand side
      wfn%qsh(ish, 1) = wfn%qsh(ish, 1) + step
      call exchange%get_energy_w_overlap(mol, cache, wfn, overlap, er)
      ! Left hand side
      wfn%qsh(ish, 1) = wfn%qsh(ish, 1) - 2*step
      call exchange%get_energy_w_overlap(mol, cache, wfn, overlap, el)

      wfn%qsh(ish, 1) = wfn%qsh(ish, 1) + step
      numvsh(ish) = 0.5_wp*(sum(er) - sum(el))/step
   end do

   ! Analytic potential after updating the cached fock matrix
   call pot%reset()
   call exchange%get_energy_w_overlap(mol, cache, wfn, overlap, el) 
   call exchange%get_potential_w_overlap(mol, cache, wfn, overlap, pot)

   if (any(abs(pot%vsh(:, 1) - numvsh) > thr_)) then
      call test_failed(error, "Shell-resolved potential does not match")
      write(*,*) "numerical potential:"
      print'(3es21.14)', numvsh
      write(*,*) "analytical potential:"
      print'(3es21.14)', pot%vsh(: ,1)
      write(*,*) "difference:"
      print'(3es21.14)', pot%vsh(: ,1) - numvsh
   end if

end subroutine test_numpot


subroutine test_num_op_grad(error, mol, density, make_exchange, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Density matrix for this structure
   real(wp), intent(in) :: density(:, :, :)

   !> Factory to create new exchange objects
   procedure(exchange_maker) :: make_exchange

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic, nspin
   class(basis_type), allocatable :: bas
   class(qvszp_basis_type), allocatable :: qvszp_bas
   type(basis_cache) :: bcache
   class(exchange_type), allocatable :: exchange
   type(wavefunction_type) :: wfn, wfn_aux
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), allocatable :: lattr(:, :)
   real(wp), allocatable :: overlap(:,:)
   real(wp) :: thr_, energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :), ao_grad(:, :)
   real(wp), parameter :: step = 5.0e-5_wp

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   nspin = size(density, dim=3)

   ! Setup basis set
   allocate(qvszp_bas)
   call make_qvszp_basis(qvszp_bas, mol, error)
   call move_alloc(qvszp_bas, bas)
   if (allocated(error)) return

   call check(error, bas%nao, size(density, 1))
   if (allocated(error)) return

   ! Obtain EEQBC charges adapt the basis set
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
   
   ! Setup potential and wavefunction
   call new_potential(pot, mol, bas, nspin, .false.)
   call pot%reset()
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, nspin,&
      &  0.0_wp, .false.)
   wfn%density(:, :, :) = density
   
   ! Get lattice points and allocate overlap
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   allocate(overlap(bas%nao, bas%nao))

   ! Setup and calculate exchange
   call make_exchange(exchange, mol, bas)

   ! Setup the overlap
   call get_overlap(mol, lattr, cutoff, bas, bcache, overlap)

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         ! Right hand side
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
         if (allocated(error)) return
         call exchange%update(mol, cache)
         call exchange%get_energy_w_overlap(mol, cache, wfn, overlap, er)
         ! Left hand side
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
         if (allocated(error)) return
         call exchange%update(mol, cache)
         call exchange%get_energy_w_overlap(mol, cache, wfn, overlap, el)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   ! Update charges and exchange with derivatives
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .true.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1), wfn_aux%dqatdr(:, :, :, 1), &
      & wfn_aux%dqatdL(:, :, :, 1))
   if (allocated(error)) return
   call exchange%update(mol, cache)

   ! Analytic gradient of the exchange energy
   allocate(ao_grad(bas%nao, bas%nao), source=0.0_wp)
   call exchange%get_gradient_w_overlap(mol, cache, wfn, overlap, &
      & ao_grad, gradient, sigma)

   if (any(abs(gradient - numgrad) > thr_)) then
      call test_failed(error, "Gradient of exchange energy does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', gradient
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', numgrad
      write(*,*) 'Difference:'
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_num_op_grad


subroutine test_num_ao_grad(error, mol, density, qsh, make_exchange, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Density matrix for this structure
   real(wp), intent(in) :: density(:, :, :)

   !> Shell-resolved charges for this structure
   real(wp), intent(in) :: qsh(:)

   !> Factory to create new exchange objects
   procedure(exchange_maker) :: make_exchange

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: mu, nu, nspin
   class(basis_type), allocatable :: bas
   class(qvszp_basis_type), allocatable :: qvszp_bas
   type(basis_cache) :: bcache
   class(exchange_type), allocatable :: exchange
   type(wavefunction_type) :: wfn, wfn_aux
   type(container_cache) :: cache
   type(potential_type) :: pot
   real(wp), allocatable :: lattr(:, :)
   real(wp), allocatable :: overlap(:, :)
   real(wp), allocatable :: ao_grad(:, :), num_ao_grad(:, :)
   real(wp), allocatable :: gradient(:, :)
   real(wp), allocatable :: el(:), er(:)
   real(wp) :: thr_, sigma(3, 3), value
   real(wp), parameter :: step = 5.0e-5_wp

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   nspin = size(density, dim=3)

   ! Setup basis set
   allocate(qvszp_bas)
   call make_qvszp_basis(qvszp_bas, mol, error)
   call move_alloc(qvszp_bas, bas)
   if (allocated(error)) return

   call check(error, bas%nao, size(density, 1))
   if (allocated(error)) return

   ! Obtain EEQBC charges adapt the basis set
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   ! Get lattice points and allocate overlap
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   allocate(overlap(bas%nao, bas%nao))

   ! Setup potential and wavefunction
   call new_potential(pot, mol, bas, nspin, .false.)
   call pot%reset()
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, nspin, 0.0_wp, .false.)
   wfn%density(:, :, :) = density
   wfn%qsh(:, 1) = qsh

   ! Setup and calculate exchange
   call make_exchange(exchange, mol, bas)
   call exchange%update(mol, cache)

   ! Setup the inital overlap
   call get_overlap(mol, lattr, cutoff, bas, bcache, overlap)

   allocate(num_ao_grad(bas%nao, bas%nao), source = 0.0_wp )
   allocate(el(mol%nat), er(mol%nat))

   do mu = 1, bas%nao
      do nu = mu, bas%nao
         er = 0.0_wp 
         el = 0.0_wp
         ! Right hand side - This leads to a double step 
         ! as both the upper and lower triangle are updated
         overlap(mu, nu) = overlap(mu, nu) + step
         if (nu /= mu) then
            overlap(nu, mu) = overlap(nu, mu) + step
         end if
         call pot%reset()
         call exchange%get_energy_w_overlap(mol, cache, wfn, overlap, er)

         ! Left hand side - This leads to a double step 
         ! as both the upper and lower triangle are updated
         overlap(mu, nu) = overlap(mu, nu) - 2.0_wp * step
         if (nu /= mu) then
            overlap(nu, mu) = overlap(nu, mu) - 2.0_wp * step
         end if
         call pot%reset()
         call exchange%get_energy_w_overlap(mol, cache, wfn, overlap, el)

         overlap(mu, nu) = overlap(mu, nu) + step
         if (nu /= mu) then
            overlap(nu, mu) = overlap(nu, mu) + step
         end if

         ! The factor 0.25 accounts for the double step in the overlap
         if (mu == nu) then
            num_ao_grad(mu, nu) = -0.5_wp * (sum(er) - sum(el)) / step
         else
            num_ao_grad(mu, nu) = -0.25_wp * (sum(er) - sum(el)) / step
            num_ao_grad(nu, mu) = num_ao_grad(mu, nu)
         end if
      end do
   end do

   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .true.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1), wfn_aux%dqatdr(:, :, :, 1), &
      & wfn_aux%dqatdL(:, :, :, 1))
   if (allocated(error)) return
   call exchange%update(mol, cache)

   allocate(ao_grad(bas%nao, bas%nao))
   allocate(gradient(3, mol%nat))
   ao_grad(:, :) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call exchange%get_gradient_w_overlap(mol, cache, wfn, overlap, ao_grad, gradient, sigma)

   if (any(abs(ao_grad - num_ao_grad) > thr_)) then
      call test_failed(error, "AO gradient of exchange energy does not match numerical reference")
      write(*,*) 'Analytic AO gradient:'
      print'(3es21.14)', ao_grad
      write(*,*) 'Numeric AO gradient:'
      print'(3es21.14)', num_ao_grad
      write(*,*) 'Difference:'
      print'(3es21.14)', ao_grad - num_ao_grad
      write(*,*) "max error", maxval(abs(ao_grad - num_ao_grad))
   end if

end subroutine test_num_ao_grad


subroutine test_gamma_fock_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(2) = [&
      & -6.66133814775094E-16_wp,  4.44089209850063E-16_wp]

   real(wp), parameter :: ref_g_mulliken(2, 2) = reshape([&
      &  3.30673408241517E-1_wp,  2.23195931032271E-2_wp,  2.23195931032271E-2_wp, &
      &  3.30673408241517E-1_wp], shape(ref_g_mulliken))

   real(wp), parameter :: ref_g_onsfx(1, 1, 2) = reshape([&
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp], shape(ref_g_onsfx))

   real(wp), parameter :: ref_dgdq_onsfx(1, 1, 2) = reshape([&
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp], shape(ref_dgdq_onsfx))

   real(wp), parameter :: ref_g_onsri(1, 2) = reshape([&
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp], shape(ref_g_onsri))

   real(wp), parameter :: ref_dgdq_onsri(1, 2) = reshape([&
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp], shape(ref_dgdq_onsri))

   real(wp), parameter :: ref_g_bocorr(2, 2) = reshape([&
      &  0.00000000000000E+0_wp, -9.86324920000000E-3_wp, -9.86324920000000E-3_wp, &
      &  0.00000000000000E+0_wp], shape(ref_g_bocorr))

   call get_structure(mol, "MB16-43", "H2")
   call test_gamma_generic(error, mol, make_exchange_gxtb, qsh, ref_g_mulliken, &
      & ref_g_onsfx, ref_dgdq_onsfx, ref_g_onsri, ref_dgdq_onsri, ref_g_bocorr, &
      & thr_in=thr1)

end subroutine test_gamma_fock_h2


subroutine test_gamma_fock_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(3) = [&
      &  1.88324089567125E-1_wp,  2.01980267353847E-1_wp, -3.90304356969502E-1_wp]

   real(wp), parameter :: ref_g_mulliken(3, 3) = reshape([&
      &  4.83274044865195E-2_wp,  2.07019966609833E-2_wp,  1.08314047047257E-2_wp, &
      &  2.07019966609833E-2_wp,  2.49298976776369E-1_wp,  9.94675686538323E-2_wp, &
      &  1.08314047047257E-2_wp,  9.94675686538323E-2_wp,  3.30673408241517E-1_wp],&
      & shape(ref_g_mulliken))

   real(wp), parameter :: ref_g_onsfx(2, 2, 2) = reshape([&
      &  0.00000000000000E+0_wp,  6.84282722695823E-3_wp,  6.84282722695823E-3_wp, &
      &  2.23716194462120E-3_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp], shape(ref_g_onsfx))

   real(wp), parameter :: ref_dgdq_onsfx(2, 2, 3) = reshape([&
      & -0.00000000000000E+0_wp, -4.47574875000000E-3_wp, -4.47574875000000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -2.23787437500000E-3_wp, &
      & -2.23787437500000E-3_wp, -1.38420975000000E-3_wp, -0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp],&
      & shape(ref_dgdq_onsfx))

   real(wp), parameter :: ref_g_onsri(2, 2) = reshape([&
      &  0.00000000000000E+0_wp,  3.72860324103533E-4_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp], shape(ref_g_onsri))

   real(wp), parameter :: ref_dgdq_onsri(2, 3) = reshape([&
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -2.30701625000000E-4_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp],&
      & shape(ref_dgdq_onsri))

   real(wp), parameter :: ref_g_bocorr(2, 2) = reshape([&
      &  0.00000000000000E+0_wp, -5.30709137423993E-3_wp, -5.30709137423993E-3_wp, &
      &  0.00000000000000E+0_wp], shape(ref_g_bocorr))

   call get_structure(mol, "MB16-43", "LiH")
   call test_gamma_generic(error, mol, make_exchange_gxtb, qsh, ref_g_mulliken, &
      & ref_g_onsfx, ref_dgdq_onsfx, ref_g_onsri, ref_dgdq_onsri, ref_g_bocorr, &
      & thr_in=thr1)

end subroutine test_gamma_fock_lih


subroutine test_gamma_fock_no(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(4) = [&
      & -3.72959864890934E-1_wp,  4.88258685456765E-1_wp, -2.06922792771332E-1_wp, &
      &  9.16239721052889E-2_wp]

   real(wp), parameter :: ref_g_mulliken(4, 4) = reshape([&
      &  5.03076032500885E-1_wp,  9.10400061909871E-2_wp,  3.94828343507272E-2_wp, &
      &  1.09636905079235E-1_wp,  9.10400061909871E-2_wp,  2.51971568427881E-1_wp, &
      &  1.21685145984265E-1_wp,  2.00826181059149E-1_wp,  3.94828343507272E-2_wp, &
      &  1.21685145984265E-1_wp,  6.27814264178743E-1_wp,  9.37901743049908E-2_wp, &
      &  1.09636905079235E-1_wp,  2.00826181059149E-1_wp,  9.37901743049908E-2_wp, &
      &  2.14291091273156E-1_wp], shape(ref_g_mulliken))

   real(wp), parameter :: ref_g_onsfx(2, 2, 2) = reshape([&
      &  0.00000000000000E+0_wp,  2.33160325416717E-2_wp,  2.33160325416717E-2_wp, &
      &  3.98633121439628E-3_wp,  0.00000000000000E+0_wp,  2.66328504083408E-2_wp, &
      &  2.66328504083408E-2_wp,  5.69637316695407E-3_wp], shape(ref_g_onsfx))

   real(wp), parameter :: ref_dgdq_onsfx(2, 2, 4) = reshape([&
      & -0.00000000000000E+0_wp, -1.19752875000000E-2_wp, -1.19752875000000E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -5.98764375000000E-3_wp, &
      & -5.98764375000000E-3_wp, -2.99741475000000E-3_wp, -0.00000000000000E+0_wp, &
      & -1.34557417500000E-2_wp, -1.34557417500000E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -6.72787087500000E-3_wp, -6.72787087500000E-3_wp, &
      & -3.29926575000000E-3_wp], shape(ref_dgdq_onsfx))

   real(wp), parameter :: ref_g_onsri(2, 2) = reshape([&
      &  0.00000000000000E+0_wp,  6.64388535732714E-4_wp,  0.00000000000000E+0_wp, &
      &  9.49395527825677E-4_wp], shape(ref_g_onsri))

   real(wp), parameter :: ref_dgdq_onsri(2, 4) = reshape([&
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.99569125000000E-4_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -5.49877625000000E-4_wp], shape(ref_dgdq_onsri))

   real(wp), parameter :: ref_g_bocorr(2, 2) = reshape([&
      &  0.00000000000000E+0_wp, -2.09810627867003E-2_wp, -2.09810627867003E-2_wp, &
      &  0.00000000000000E+0_wp], shape(ref_g_bocorr))

   call get_structure(mol, "MB16-43", "NO")
   call test_gamma_generic(error, mol, make_exchange_gxtb, qsh, ref_g_mulliken, &
      & ref_g_onsfx, ref_dgdq_onsfx, ref_g_onsri, ref_dgdq_onsri, ref_g_bocorr, &
      & thr_in=thr1)

end subroutine test_gamma_fock_no


subroutine test_gamma_fock_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(6) = [&
      & -1.65829465017236E-1_wp,  8.50038195478131E-2_wp,  8.08256454270758E-2_wp, &
      & -1.65829465017241E-1_wp,  8.50038195477536E-2_wp,  8.08256454270730E-2_wp]

   real(wp), parameter :: ref_g_mulliken(6, 6) = reshape([&
      &  3.13821024392893E-1_wp,  6.33858386215170E-2_wp,  1.80442929446396E-3_wp, &
      &  2.55301101250079E-2_wp,  9.07292750669926E-2_wp,  4.53480910058135E-3_wp, &
      &  6.33858386215170E-2_wp,  1.95804674140298E-1_wp,  1.26733605145401E-2_wp, &
      &  9.07292750669926E-2_wp,  1.67159421942208E-1_wp,  2.81016705194958E-2_wp, &
      &  1.80442929446396E-3_wp,  1.26733605145401E-2_wp,  3.56845869750545E-3_wp, &
      &  4.53480910058135E-3_wp,  2.81016705194958E-2_wp,  1.53402859416096E-2_wp, &
      &  2.55301101250079E-2_wp,  9.07292750669926E-2_wp,  4.53480910058135E-3_wp, &
      &  3.13821024392893E-1_wp,  6.33858386215170E-2_wp,  1.80442929446396E-3_wp, &
      &  9.07292750669926E-2_wp,  1.67159421942208E-1_wp,  2.81016705194958E-2_wp, &
      &  6.33858386215170E-2_wp,  1.95804674140298E-1_wp,  1.26733605145401E-2_wp, &
      &  4.53480910058135E-3_wp,  2.81016705194958E-2_wp,  1.53402859416096E-2_wp, &
      &  1.80442929446396E-3_wp,  1.26733605145401E-2_wp,  3.56845869750545E-3_wp],&
      & shape(ref_g_mulliken))
      
   real(wp), parameter :: ref_g_onsfx(3, 3, 2) = reshape([&
      &  0.00000000000000E+0_wp,  1.76626048596496E-2_wp,  8.24208462856111E-3_wp, &
      &  1.76626048596496E-2_wp,  4.00725875262220E-3_wp,  8.45766855160350E-3_wp, &
      &  8.24208462856111E-3_wp,  8.45766855160350E-3_wp,  3.09321405262857E-3_wp, &
      &  0.00000000000000E+0_wp,  1.76626048596499E-2_wp,  8.24208462856113E-3_wp, &
      &  1.76626048596499E-2_wp,  4.00725875262234E-3_wp,  8.45766855160364E-3_wp, &
      &  8.24208462856113E-3_wp,  8.45766855160364E-3_wp,  3.09321405262858E-3_wp],&
      & shape(ref_g_onsfx))

   real(wp), parameter :: ref_dgdq_onsfx(3, 3, 6) = reshape([&
      & -0.00000000000000E+0_wp, -9.09735750000000E-3_wp, -4.19699775000000E-3_wp, &
      & -9.09735750000000E-3_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.19699775000000E-3_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -4.54867875000000E-3_wp,  0.00000000000000E+0_wp, &
      & -4.54867875000000E-3_wp, -2.31208725000000E-3_wp, -2.40894225000000E-3_wp, &
      &  0.00000000000000E+0_wp, -2.40894225000000E-3_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.04924943750000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.20447112500000E-3_wp, &
      & -1.04924943750000E-3_wp, -1.20447112500000E-3_wp, -8.69970750000000E-4_wp, &
      & -0.00000000000000E+0_wp, -9.09735750000000E-3_wp, -4.19699775000000E-3_wp, &
      & -9.09735750000000E-3_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.19699775000000E-3_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -4.54867875000000E-3_wp,  0.00000000000000E+0_wp, &
      & -4.54867875000000E-3_wp, -2.31208725000000E-3_wp, -2.40894225000000E-3_wp, &
      &  0.00000000000000E+0_wp, -2.40894225000000E-3_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.04924943750000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.20447112500000E-3_wp, &
      & -1.04924943750000E-3_wp, -1.20447112500000E-3_wp, -8.69970750000000E-4_wp],&
      & shape(ref_dgdq_onsfx))

   real(wp), parameter :: ref_g_onsri(3, 2) = reshape([&
      &  0.00000000000000E+0_wp,  6.67876458770367E-4_wp,  3.09321405262857E-4_wp, &
      &  0.00000000000000E+0_wp,  6.67876458770390E-4_wp,  3.09321405262858E-4_wp],&
      & shape(ref_g_onsri))

   real(wp), parameter :: ref_dgdq_onsri(3, 6) = reshape([&
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -3.85347875000000E-4_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -8.69970750000000E-5_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -3.85347875000000E-4_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -8.69970750000000E-5_wp],&
      & shape(ref_dgdq_onsri))

   real(wp), parameter :: ref_g_bocorr(2, 2) = reshape([&
      &  0.00000000000000E+0_wp, -2.17577953435037E-2_wp, -2.17577953435037E-2_wp, &
      &  0.00000000000000E+0_wp], shape(ref_g_bocorr))

   call get_structure(mol, "MB16-43", "S2")
   call test_gamma_generic(error, mol, make_exchange_gxtb, qsh, ref_g_mulliken, &
      & ref_g_onsfx, ref_dgdq_onsfx, ref_g_onsri, ref_dgdq_onsri, ref_g_bocorr, &
      & thr_in=thr1)

end subroutine test_gamma_fock_s2

subroutine test_gamma_fock_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(13) = [&
      &  8.56114797532967E-1_wp,  1.08081185051881E-1_wp,  8.40158726565170E-1_wp, &
      & -4.26774273578783E-2_wp, -9.50845571950689E-2_wp, -5.37262650096369E-1_wp, &
      &  4.55006307726409E-2_wp, -9.61301031913531E-2_wp, -5.37689488045761E-1_wp, &
      &  4.55254261185908E-2_wp, -9.54730110948796E-2_wp, -5.36576677043751E-1_wp, &
      &  4.55131478317075E-2_wp]

   real(wp), parameter :: ref_g_mulliken(13, 13) = reshape([&
      &  8.47539002683444E-2_wp,  4.11812648459241E-3_wp,  1.92931632001652E-2_wp, &
      &  3.35832134128565E-2_wp,  1.47466662623270E-2_wp,  6.12343120679460E-2_wp, &
      &  5.16950830163081E-2_wp,  1.47439191657366E-2_wp,  6.12213226979248E-2_wp, &
      &  5.16843912586466E-2_wp,  1.47516390928419E-2_wp,  6.12578403314810E-2_wp, &
      &  5.17144472349395E-2_wp,  4.11812648459241E-3_wp,  8.89787525770182E-3_wp, &
      &  2.84128098530754E-2_wp,  3.08288247744108E-2_wp,  1.07035590812117E-2_wp, &
      &  5.49146492666633E-2_wp,  5.42692633529841E-2_wp,  1.07015892136922E-2_wp, &
      &  5.49031933552994E-2_wp,  5.42579615452898E-2_wp,  1.07071247408624E-2_wp, &
      &  5.49353983040604E-2_wp,  5.42897331128294E-2_wp,  1.92931632001652E-2_wp, &
      &  2.84128098530754E-2_wp,  6.71688313670347E-2_wp,  2.11403075481855E-1_wp, &
      &  6.44042224466309E-2_wp,  1.36078617674309E-1_wp,  1.26238864444015E-1_wp, &
      &  6.43904471822949E-2_wp,  1.36044091207986E-1_wp,  1.26207524934737E-1_wp, &
      &  6.44291751819495E-2_wp,  1.36141205354997E-1_wp,  1.26295669936343E-1_wp, &
      &  3.35832134128565E-2_wp,  3.08288247744108E-2_wp,  2.11403075481855E-1_wp, &
      &  2.03519213159664E-1_wp,  8.82988317916133E-2_wp,  1.50900678780445E-1_wp, &
      &  1.46053825227799E-1_wp,  8.82787730666998E-2_wp,  1.50861148514530E-1_wp, &
      &  1.46015958086381E-1_wp,  8.83351762547545E-2_wp,  1.50972346250484E-1_wp, &
      &  1.46122474666890E-1_wp,  1.47466662623270E-2_wp,  1.07035590812117E-2_wp, &
      &  6.44042224466309E-2_wp,  8.82988317916133E-2_wp,  3.54512902023288E-1_wp, &
      &  7.46259472549862E-2_wp,  5.68679384793874E-2_wp,  1.54704474796601E-2_wp, &
      &  5.53087881801101E-2_wp,  4.74023446206159E-2_wp,  1.57086651187799E-2_wp, &
      &  5.60457127325977E-2_wp,  4.80534084458044E-2_wp,  6.12343120679460E-2_wp, &
      &  5.49146492666633E-2_wp,  1.36078617674309E-1_wp,  1.50900678780445E-1_wp, &
      &  7.46259472549862E-2_wp,  2.40252700329031E-1_wp,  3.49918402948991E-1_wp, &
      &  5.53087881801101E-2_wp,  9.91871960132011E-2_wp,  9.53970206357638E-2_wp, &
      &  5.60457127325977E-2_wp,  1.00283098134793E-1_wp,  9.64697539781455E-2_wp, &
      &  5.16950830163081E-2_wp,  5.42692633529841E-2_wp,  1.26238864444015E-1_wp, &
      &  1.46053825227799E-1_wp,  5.68679384793874E-2_wp,  3.49918402948991E-1_wp, &
      &  1.62017342584553E-1_wp,  4.74023446206159E-2_wp,  9.53970206357638E-2_wp, &
      &  9.18858469302032E-2_wp,  4.80534084458044E-2_wp,  9.64697539781455E-2_wp, &
      &  9.29357967695427E-2_wp,  1.47439191657366E-2_wp,  1.07015892136922E-2_wp, &
      &  6.43904471822949E-2_wp,  8.82787730666998E-2_wp,  1.54704474796601E-2_wp, &
      &  5.53087881801101E-2_wp,  4.74023446206159E-2_wp,  3.54512902023288E-1_wp, &
      &  7.46259472549862E-2_wp,  5.68679384793874E-2_wp,  1.52715848232330E-2_wp, &
      &  5.46936513746535E-2_wp,  4.68588206442025E-2_wp,  6.12213226979248E-2_wp, &
      &  5.49031933552994E-2_wp,  1.36044091207986E-1_wp,  1.50861148514530E-1_wp, &
      &  5.53087881801101E-2_wp,  9.91871960132011E-2_wp,  9.53970206357638E-2_wp, &
      &  7.46259472549862E-2_wp,  2.40252700329031E-1_wp,  3.49918402948991E-1_wp, &
      &  5.46936513746535E-2_wp,  9.82740167579058E-2_wp,  9.45029260243677E-2_wp, &
      &  5.16843912586466E-2_wp,  5.42579615452898E-2_wp,  1.26207524934737E-1_wp, &
      &  1.46015958086381E-1_wp,  4.74023446206159E-2_wp,  9.53970206357638E-2_wp, &
      &  9.18858469302032E-2_wp,  5.68679384793874E-2_wp,  3.49918402948991E-1_wp, &
      &  1.62017342584553E-1_wp,  4.68588206442025E-2_wp,  9.45029260243677E-2_wp, &
      &  9.10105576258867E-2_wp,  1.47516390928419E-2_wp,  1.07071247408624E-2_wp, &
      &  6.44291751819495E-2_wp,  8.83351762547545E-2_wp,  1.57086651187799E-2_wp, &
      &  5.60457127325977E-2_wp,  4.80534084458044E-2_wp,  1.52715848232330E-2_wp, &
      &  5.46936513746535E-2_wp,  4.68588206442025E-2_wp,  3.54512902023288E-1_wp, &
      &  7.46259472549862E-2_wp,  5.68679384793874E-2_wp,  6.12578403314810E-2_wp, &
      &  5.49353983040604E-2_wp,  1.36141205354997E-1_wp,  1.50972346250484E-1_wp, &
      &  5.60457127325977E-2_wp,  1.00283098134793E-1_wp,  9.64697539781455E-2_wp, &
      &  5.46936513746535E-2_wp,  9.82740167579058E-2_wp,  9.45029260243677E-2_wp, &
      &  7.46259472549862E-2_wp,  2.40252700329031E-1_wp,  3.49918402948991E-1_wp, &
      &  5.17144472349395E-2_wp,  5.42897331128294E-2_wp,  1.26295669936343E-1_wp, &
      &  1.46122474666890E-1_wp,  4.80534084458044E-2_wp,  9.64697539781455E-2_wp, &
      &  9.29357967695427E-2_wp,  4.68588206442025E-2_wp,  9.45029260243677E-2_wp, &
      &  9.10105576258867E-2_wp,  5.68679384793874E-2_wp,  3.49918402948991E-1_wp, &
      &  1.62017342584553E-1_wp], shape(ref_g_mulliken))

   real(wp), parameter :: ref_g_onsfx(4, 4, 4) = reshape([&
      &  0.00000000000000E+0_wp,  2.41913154786025E-3_wp,  8.07796370420516E-4_wp, &
      &  1.38678968364880E-4_wp,  2.41913154786025E-3_wp,  1.43591763346926E-3_wp, &
      &  7.37778237464658E-4_wp,  1.08989983573585E-4_wp,  8.07796370420516E-4_wp, &
      &  7.37778237464658E-4_wp,  1.64160261943114E-3_wp,  1.49931848987138E-3_wp, &
      &  1.38678968364880E-4_wp,  1.08989983573585E-4_wp,  1.49931848987138E-3_wp, &
      &  3.39056525939677E-3_wp,  0.00000000000000E+0_wp,  2.16666156333045E-2_wp, &
      &  9.12924758867307E-3_wp,  0.00000000000000E+0_wp,  2.16666156333045E-2_wp, &
      &  5.90063394812268E-3_wp,  1.11131843404128E-2_wp,  0.00000000000000E+0_wp, &
      &  9.12924758867307E-3_wp,  1.11131843404128E-2_wp,  3.60140577937541E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  2.16791173698051E-2_wp,  9.13423654780183E-3_wp,  0.00000000000000E+0_wp, &
      &  2.16791173698051E-2_wp,  5.90170322168376E-3_wp,  1.11142939257682E-2_wp, &
      &  0.00000000000000E+0_wp,  9.13423654780183E-3_wp,  1.11142939257682E-2_wp, &
      &  3.60138091123032E-3_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  2.16670671320642E-2_wp,  9.13109718131150E-3_wp, &
      &  0.00000000000000E+0_wp,  2.16670671320642E-2_wp,  5.89891551411318E-3_wp, &
      &  1.11113310214429E-2_wp,  0.00000000000000E+0_wp,  9.13109718131150E-3_wp, &
      &  1.11113310214429E-2_wp,  3.60139322556626E-3_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp], shape(ref_g_onsfx))

   real(wp), parameter :: ref_dgdq_onsfx(4, 4, 13) = reshape([&
      & -0.00000000000000E+0_wp, -2.66416425000000E-3_wp, -1.07415825000000E-3_wp, &
      & -1.43352000000000E-4_wp, -2.66416425000000E-3_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.07415825000000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -1.43352000000000E-4_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.33208212500000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.33208212500000E-3_wp, &
      & -8.39668500000000E-4_wp, -2.37364875000000E-4_wp, -3.07972500000000E-5_wp, &
      &  0.00000000000000E+0_wp, -2.37364875000000E-4_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -3.07972500000000E-5_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -2.68539562500000E-4_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.18682437500000E-4_wp, &
      &  0.00000000000000E+0_wp, -2.68539562500000E-4_wp, -1.18682437500000E-4_wp, &
      & -5.87082375000000E-4_wp, -2.32311750000000E-4_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -2.32311750000000E-4_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -1.79190000000000E-5_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -7.69931250000000E-6_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.16155875000000E-4_wp, &
      & -1.79190000000000E-5_wp, -7.69931250000000E-6_wp, -1.16155875000000E-4_wp, &
      & -4.63482937500000E-4_wp, -0.00000000000000E+0_wp, -9.93017025000000E-3_wp, &
      & -4.80008925000000E-3_wp,  0.00000000000000E+0_wp, -9.93017025000000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.80008925000000E-3_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.96508512500000E-3_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.96508512500000E-3_wp, -2.50510425000000E-3_wp, -2.67731062500000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -2.67731062500000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.20002231250000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -1.33865531250000E-3_wp,  0.00000000000000E+0_wp, -1.20002231250000E-3_wp, &
      & -1.33865531250000E-3_wp, -1.00293600000000E-3_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -0.00000000000000E+0_wp, -9.93017025000000E-3_wp, &
      & -4.80008925000000E-3_wp,  0.00000000000000E+0_wp, -9.93017025000000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.80008925000000E-3_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.96508512500000E-3_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.96508512500000E-3_wp, -2.50510425000000E-3_wp, -2.67731062500000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -2.67731062500000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.20002231250000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -1.33865531250000E-3_wp,  0.00000000000000E+0_wp, -1.20002231250000E-3_wp, &
      & -1.33865531250000E-3_wp, -1.00293600000000E-3_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -0.00000000000000E+0_wp, -9.93017025000000E-3_wp, &
      & -4.80008925000000E-3_wp,  0.00000000000000E+0_wp, -9.93017025000000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.80008925000000E-3_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.96508512500000E-3_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.96508512500000E-3_wp, -2.50510425000000E-3_wp, -2.67731062500000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -2.67731062500000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.20002231250000E-3_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -1.33865531250000E-3_wp,  0.00000000000000E+0_wp, -1.20002231250000E-3_wp, &
      & -1.33865531250000E-3_wp, -1.00293600000000E-3_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp], shape(ref_dgdq_onsfx))

   real(wp), parameter :: ref_g_onsri(4, 4) = reshape([&
      &  0.00000000000000E+0_wp,  2.39319605578211E-4_wp,  1.64160261943114E-4_wp, &
      &  2.42183232814055E-4_wp,  0.00000000000000E+0_wp,  9.83438991353779E-4_wp, &
      &  3.60140577937541E-4_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  9.83617203613960E-4_wp,  3.60138091123032E-4_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  9.83152585685530E-4_wp,  3.60139322556626E-4_wp, &
      &  0.00000000000000E+0_wp], shape(ref_g_onsri))

   real(wp), parameter :: ref_dgdq_onsri(4, 13) = reshape([&
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.39944750000000E-4_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -5.87082375000000E-5_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.31059241071429E-5_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.17517375000000E-4_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.00293600000000E-4_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.17517375000000E-4_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.00293600000000E-4_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.17517375000000E-4_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.00293600000000E-4_wp, &
      &  0.00000000000000E+0_wp], shape(ref_dgdq_onsri))

   real(wp), parameter :: ref_g_bocorr(4, 4) = reshape([&
      &  0.00000000000000E+0_wp, -3.01380798699169E-3_wp, -3.01098405501475E-3_wp, &
      & -3.01892916635643E-3_wp, -3.01380798699169E-3_wp,  0.00000000000000E+0_wp, &
      & -6.53231448608759E-4_wp, -7.04333348409273E-4_wp, -3.01098405501475E-3_wp, &
      & -6.53231448608759E-4_wp,  0.00000000000000E+0_wp, -6.12313934567724E-4_wp, &
      & -3.01892916635643E-3_wp, -7.04333348409273E-4_wp, -6.12313934567724E-4_wp, &
      &  0.00000000000000E+0_wp], shape(ref_g_bocorr))

   call get_structure(mol, "f-block", "CeCl3")
   call test_gamma_generic(error, mol, make_exchange_gxtb, qsh, ref_g_mulliken, &
      & ref_g_onsfx, ref_dgdq_onsfx, ref_g_onsri, ref_dgdq_onsri, ref_g_bocorr, &
      & thr_in=thr1)

end subroutine test_gamma_fock_cecl3


subroutine test_e_fock_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(2) = [&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp]

   real(wp), parameter :: density(2, 2, 1) = reshape([&
      & 5.93683737100003E-01_wp, 5.93683737100003E-01_wp, 5.93683737100003E-01_wp, &
      & 5.93683737100003E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "H2")
   call test_energy_generic(error, mol, density, qsh, make_exchange_gxtb, &
      & -0.186359749872183_wp, thr_in=thr1)

end subroutine test_e_fock_h2


subroutine test_e_fock_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(3) = [&
      & 1.88324720841056E-01_wp, 2.01979331165313E-01_wp,-3.90304052054899E-01_wp]

   real(wp), parameter :: density(5, 5, 1) = reshape([&
      & 7.43135570300144E-02_wp, 0.00000000000000E+00_wp, 1.15037957895324E-01_wp, &
      & 0.00000000000000E+00_wp, 2.77066781440927E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.15037957895324E-01_wp, 0.00000000000000E+00_wp, &
      & 1.78079643683070E-01_wp, 0.00000000000000E+00_wp, 4.28901508842071E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 2.77066781440927E-01_wp, &
      & 0.00000000000000E+00_wp, 4.28901508842071E-01_wp, 0.00000000000000E+00_wp, &
      & 1.03300130482288E+00_wp], shape(density))

   call get_structure(mol, "MB16-43", "LiH")
   call test_energy_generic(error, mol, density, qsh, make_exchange_gxtb, &
      & -0.206091553517411_wp, thr_in=thr1)

end subroutine test_e_fock_lih


subroutine test_e_fock_no(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(4) = [&
      &-3.73400110850439E-01_wp, 4.85422151137329E-01_wp,-2.06420659772590E-01_wp, &
      & 9.43986193854895E-02_wp]

   real(wp), parameter :: density(8, 8, 2) = reshape([&
      & 9.42019520534289E-01_wp,-4.20333969379186E-17_wp,-3.13829020687168E-01_wp, &
      &-9.08969029521383E-17_wp,-1.91142085625176E-01_wp,-1.26214259304181E-16_wp, &
      &-2.54599763899220E-02_wp, 1.39297732856661E-16_wp,-4.20333969379186E-17_wp, &
      & 9.80685944415004E-01_wp, 1.41151208526390E-16_wp,-2.78804784267892E-01_wp, &
      & 4.82906497923563E-17_wp,-2.21785459465347E-01_wp,-1.81892320999620E-17_wp, &
      & 2.26652363982251E-01_wp,-3.13829020687168E-01_wp, 1.41151208526390E-16_wp, &
      & 4.49237030058337E-01_wp,-5.84351340746751E-18_wp,-1.32690321562938E-03_wp, &
      &-1.23697560507224E-16_wp,-4.06976420678055E-01_wp,-6.18707795513060E-17_wp, &
      &-9.08969029521383E-17_wp,-2.78804784267892E-01_wp,-5.84351340746751E-18_wp, &
      & 3.60805147157226E-01_wp,-1.18727758269954E-17_wp, 2.26652363982254E-01_wp, &
      &-5.36317407847742E-17_wp, 2.82142220450980E-01_wp,-1.91142085625176E-01_wp, &
      & 4.82906497923563E-17_wp,-1.32690321562938E-03_wp,-1.18727758269954E-17_wp, &
      & 1.01793801861039E+00_wp, 1.17305253624528E-16_wp, 2.62429582506685E-01_wp, &
      & 5.30251183029614E-17_wp,-1.26214259304181E-16_wp,-2.21785459465347E-01_wp, &
      &-1.23697560507224E-16_wp, 2.26652363982254E-01_wp, 1.17305253624528E-16_wp, &
      & 1.01695424802969E+00_wp,-9.55620890266933E-17_wp,-1.84255425292091E-01_wp, &
      &-2.54599763899220E-02_wp,-1.81892320999620E-17_wp,-4.06976420678055E-01_wp, &
      &-5.36317407847742E-17_wp, 2.62429582506685E-01_wp,-9.55620890266933E-17_wp, &
      & 5.34554575158069E-01_wp, 8.15142917082273E-17_wp, 1.39297732856661E-16_wp, &
      & 2.26652363982251E-01_wp,-6.18707795513060E-17_wp, 2.82142220450980E-01_wp, &
      & 5.30251183029614E-17_wp,-1.84255425292091E-01_wp, 8.15142917082273E-17_wp, &
      & 6.07289829145903E-01_wp, 9.54493646021522E-01_wp, 7.93662775608665E-17_wp, &
      &-3.10425447985195E-01_wp,-2.28608225800888E-17_wp,-2.03716269338494E-01_wp, &
      &-1.39160964040937E-16_wp,-1.78381858562025E-02_wp, 1.41307783093691E-16_wp, &
      & 7.93662775608665E-17_wp, 2.44187320087128E-01_wp,-7.34349602106188E-17_wp, &
      & 2.39137447190470E-05_wp,-3.77198716005244E-17_wp, 3.65846783197429E-01_wp, &
      & 4.41325052793164E-17_wp, 8.34920778010917E-06_wp,-3.10425447985195E-01_wp, &
      &-7.34349602106188E-17_wp, 4.40523292733754E-01_wp, 6.55385937589603E-17_wp, &
      &-4.71240774340272E-03_wp, 4.36977656014520E-17_wp,-4.09806796492085E-01_wp, &
      &-4.18179359368282E-17_wp,-2.28608225800888E-17_wp, 2.39137447190470E-05_wp, &
      & 6.55385937589603E-17_wp, 2.44240488726901E-01_wp,-1.43156675593419E-16_wp, &
      & 8.34920777982584E-06_wp,-7.51658113014306E-17_wp, 3.65865346413942E-01_wp, &
      &-2.03716269338494E-01_wp,-3.77198716005244E-17_wp,-4.71240774340272E-03_wp, &
      &-1.43156675593419E-16_wp, 1.03061284967373E+00_wp, 1.08810042276255E-16_wp, &
      & 2.54769757899732E-01_wp, 1.89391492904298E-16_wp,-1.39160964040937E-16_wp, &
      & 3.65846783197429E-01_wp, 4.36977656014520E-17_wp, 8.34920777982584E-06_wp, &
      & 1.08810042276255E-16_wp, 5.48119654546661E-01_wp, 6.59189420851037E-18_wp, &
      &-2.86536312039953E-05_wp,-1.78381858562025E-02_wp, 4.41325052793164E-17_wp, &
      &-4.09806796492085E-01_wp,-7.51658113014306E-17_wp, 2.54769757899732E-01_wp, &
      & 6.59189420851037E-18_wp, 5.36711347767033E-01_wp, 5.71457519096238E-17_wp, &
      & 1.41307783093691E-16_wp, 8.34920778010917E-06_wp,-4.18179359368282E-17_wp, &
      & 3.65865346413942E-01_wp, 1.89391492904298E-16_wp,-2.86536312039953E-05_wp, &
      & 5.71457519096238E-17_wp, 5.48055947477209E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "NO")
   call test_energy_generic(error, mol, density, qsh, make_exchange_gxtb, &
      & -1.99982865437890_wp, thr_in=thr1)

end subroutine test_e_fock_no


subroutine test_e_fock_n2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(4) = [&
      &-3.82062631895399E-01_wp, 3.82062631845031E-01_wp,-3.82062631895401E-01_wp, &
      & 3.82062631845026E-01_wp]

   real(wp), parameter :: density(8, 8, 1) = reshape([&
      & 1.89799502629230E+00_wp, 8.29256049064486E-17_wp,-6.36539727548007E-01_wp, &
      &-3.22185922461283E-18_wp,-3.58808460007242E-01_wp, 3.14519781711840E-16_wp, &
      &-6.05878687705267E-02_wp, 1.22104063460849E-17_wp, 8.29256049064486E-17_wp, &
      & 7.57491744560831E-01_wp, 2.79922050388198E-16_wp, 3.44290965205769E-18_wp, &
      &-3.70529841988989E-17_wp, 7.57491744560830E-01_wp, 9.79432668805693E-17_wp, &
      &-5.04069129978789E-16_wp,-6.36539727548007E-01_wp, 2.79922050388198E-16_wp, &
      & 1.02865522729829E+00_wp,-3.29194092537802E-17_wp, 6.05878687705266E-02_wp, &
      & 2.42877897867654E-16_wp,-8.13312204080213E-01_wp,-9.14238869535326E-17_wp, &
      &-3.22185922461283E-18_wp, 3.44290965205769E-18_wp,-3.29194092537802E-17_wp, &
      & 7.57491744560831E-01_wp, 2.19607501004388E-16_wp,-5.55111512312578E-17_wp, &
      & 1.72546470413369E-17_wp, 7.57491744560830E-01_wp,-3.58808460007242E-01_wp, &
      &-3.70529841988989E-17_wp, 6.05878687705266E-02_wp, 2.19607501004388E-16_wp, &
      & 1.89799502629230E+00_wp,-4.27209889863763E-16_wp, 6.36539727548011E-01_wp, &
      &-5.68753811725907E-17_wp, 3.14519781711840E-16_wp, 7.57491744560830E-01_wp, &
      & 2.42877897867654E-16_wp,-5.55111512312578E-17_wp,-4.27209889863763E-16_wp, &
      & 7.57491744560829E-01_wp,-5.70717653095793E-17_wp,-6.10622663543836E-16_wp, &
      &-6.05878687705267E-02_wp, 9.79432668805693E-17_wp,-8.13312204080213E-01_wp, &
      & 1.72546470413369E-17_wp, 6.36539727548011E-01_wp,-5.70717653095793E-17_wp, &
      & 1.02865522729830E+00_wp,-1.44135936695812E-17_wp, 1.22104063460849E-17_wp, &
      &-5.04069129978789E-16_wp,-9.14238869535326E-17_wp, 7.57491744560830E-01_wp, &
      &-5.68753811725907E-17_wp,-6.10622663543836E-16_wp,-1.44135936695812E-17_wp, &
      & 7.57491744560828E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "N2")
   call test_energy_generic(error, mol, density, qsh, make_exchange_gxtb, &
      & -1.79390498153404_wp, thr_in=thr1)

end subroutine test_e_fock_n2


subroutine test_e_fock_h2o(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(4) = [&
      & 3.14505633333656E-01_wp, 3.14503683456474E-01_wp,-1.33128918348271E-01_wp, &
      &-4.95880398491699E-01_wp]

   real(wp), parameter :: density(6, 6, 1) = reshape([&
      & 4.00153972503678E-01_wp,-9.25481468955952E-02_wp, 3.73690554621718E-02_wp, &
      & 1.25444506797340E-16_wp,-6.25268361696945E-01_wp, 5.96674376399834E-02_wp, &
      &-9.25481468955952E-02_wp, 4.00158952840405E-01_wp, 3.73744409137897E-02_wp, &
      &-3.33459827262735E-16_wp, 1.51987149503491E-01_wp, 6.09447275536468E-01_wp, &
      & 3.73690554621718E-02_wp, 3.73744409137897E-02_wp, 1.77094632228289E+00_wp, &
      & 1.00394384848864E-16_wp, 3.16552288124044E-01_wp,-4.47523974864359E-01_wp, &
      & 1.25444506797340E-16_wp,-3.33459827262735E-16_wp, 1.00394384848864E-16_wp, &
      & 2.00000000000000E+00_wp,-4.52473080326630E-17_wp, 2.55826164308748E-16_wp, &
      &-6.25268361696945E-01_wp, 1.51987149503491E-01_wp, 3.16552288124044E-01_wp, &
      &-4.52473080326630E-17_wp, 1.05658005895389E+00_wp,-1.93369383729808E-01_wp, &
      & 5.96674376399834E-02_wp, 6.09447275536468E-01_wp,-4.47523974864359E-01_wp, &
      & 2.55826164308748E-16_wp,-1.93369383729808E-01_wp, 1.19318512323793E+00_wp],&
      & shape(density))

   call get_structure(mol, "ICE10", "gas")
   call test_energy_generic(error, mol, density, qsh, make_exchange_gxtb, &
      & -1.24357085549152_wp, thr_in=thr1)

end subroutine test_e_fock_h2o


subroutine test_e_fock_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(6) = [&
      &-1.65829465642048E-01_wp, 8.50037124926217E-02_wp, 8.08257531070440E-02_wp, &
      &-1.65829465642047E-01_wp, 8.50037124926244E-02_wp, 8.08257531070442E-02_wp]

   real(wp), parameter :: density(18, 18, 1) = reshape([&
      & 2.01512275101693E+00_wp, 8.69702020476042E-17_wp,-3.99929220910119E-01_wp, &
      & 1.55773403391567E-16_wp, 1.63576811366064E-04_wp,-7.15317852021623E-18_wp, &
      &-3.51352591933450E-02_wp, 4.51032127453589E-20_wp, 3.55558297044844E-05_wp, &
      &-2.61461507468549E-01_wp,-3.43182042645543E-17_wp, 8.46120361482158E-02_wp, &
      & 1.67783359907526E-16_wp,-2.90703197793973E-05_wp,-6.89667265911756E-18_wp, &
      & 1.18677054674651E-02_wp, 4.56621850648170E-18_wp,-6.31886225739145E-06_wp, &
      & 8.69702020476042E-17_wp, 1.28272149679125E+00_wp,-7.54865599148717E-18_wp, &
      & 5.85799207483420E-01_wp, 2.07694503111621E-18_wp, 1.00951275464803E-02_wp, &
      & 4.98101003909952E-17_wp,-4.32459272824377E-02_wp,-8.03466743686323E-18_wp, &
      &-1.63314702479790E-17_wp, 3.34748212919108E-01_wp, 7.57576586626171E-17_wp, &
      &-5.90366231871004E-01_wp, 2.30216039812600E-19_wp,-8.76093662322534E-02_wp, &
      &-1.39662490020440E-16_wp,-5.29272140403637E-02_wp, 5.43883875213128E-17_wp, &
      &-3.99929220910119E-01_wp,-7.54865599148717E-18_wp, 8.16076846648327E-01_wp, &
      & 6.11285133940764E-17_wp,-2.68312706374975E-05_wp, 6.41243357319597E-18_wp, &
      & 5.39351527038126E-02_wp,-1.27414642211600E-16_wp,-5.83217194169540E-06_wp, &
      &-8.46120361482148E-02_wp, 1.09460039181865E-16_wp,-7.72403994815329E-01_wp, &
      & 3.60887054888469E-17_wp,-1.48774341268614E-07_wp,-9.62026471439360E-21_wp, &
      & 4.74250301651590E-02_wp,-4.22680717602266E-17_wp,-3.23382947650264E-08_wp, &
      & 1.55773403391567E-16_wp, 5.85799207483420E-01_wp, 6.11285133940764E-17_wp, &
      & 1.53738567184553E+00_wp, 3.47097293878672E-18_wp,-4.32459272824378E-02_wp, &
      & 1.47751293249772E-17_wp,-8.70515120970188E-03_wp, 2.37034465968426E-18_wp, &
      &-3.91896272003785E-17_wp,-5.90366231871004E-01_wp, 3.60985123239876E-16_wp, &
      & 7.80986178040750E-02_wp, 2.08615446024900E-19_wp,-5.29272140403639E-02_wp, &
      &-1.37889229695067E-16_wp,-1.10618385584490E-01_wp,-1.41966749806042E-16_wp, &
      & 1.63576811366064E-04_wp, 2.07694503111621E-18_wp,-2.68312706374975E-05_wp, &
      & 3.47097293878672E-18_wp, 1.33450121877225E-08_wp,-2.26220082384044E-19_wp, &
      &-2.54870417630385E-06_wp,-2.32060260310431E-19_wp, 2.90073499286456E-09_wp, &
      &-2.90703197794302E-05_wp,-2.89827131754869E-18_wp, 1.48774341197204E-07_wp, &
      &-2.68991192418201E-18_wp,-2.95700562701623E-09_wp,-1.80963233126845E-19_wp, &
      & 1.42873977881777E-06_wp,-2.73041420985304E-19_wp,-6.42748734559565E-10_wp, &
      &-7.15317852021623E-18_wp, 1.00951275464803E-02_wp, 6.41243357319597E-18_wp, &
      &-4.32459272824378E-02_wp,-2.26220082384044E-19_wp, 6.15075023236154E-03_wp, &
      & 5.73107898983234E-18_wp, 4.53852546312478E-03_wp,-2.59905923013698E-19_wp, &
      &-5.10605669701393E-18_wp, 8.76093662322536E-02_wp,-8.07837585297486E-18_wp, &
      & 5.29272140403642E-02_wp,-1.19981759421529E-21_wp, 1.87463408525633E-04_wp, &
      &-9.45961822144709E-18_wp, 3.32539639889863E-03_wp, 2.06766115803617E-18_wp, &
      &-3.51352591933450E-02_wp, 4.98101003909952E-17_wp, 5.39351527038126E-02_wp, &
      & 1.47751293249772E-17_wp,-2.54870417630385E-06_wp, 5.73107898983234E-18_wp, &
      & 3.73731491668374E-03_wp,-3.41754956218296E-18_wp,-5.53998399309322E-07_wp, &
      & 1.18677054674651E-02_wp, 8.66539900807442E-17_wp,-4.74250301651592E-02_wp, &
      & 4.93091128376413E-17_wp, 1.42873977881271E-06_wp,-2.89559933915559E-18_wp, &
      & 2.76687919541808E-03_wp,-3.62779089330793E-18_wp, 3.10557638597950E-07_wp, &
      & 4.51032127453589E-20_wp,-4.32459272824377E-02_wp,-1.27414642211600E-16_wp, &
      &-8.70515120970188E-03_wp,-2.32060260310431E-19_wp, 4.53852546312478E-03_wp, &
      &-3.41754956218296E-18_wp, 8.12378097099090E-03_wp, 3.61588930516630E-19_wp, &
      & 1.38092961541038E-16_wp, 5.29272140403641E-02_wp, 1.63344576636165E-16_wp, &
      & 1.10618385584490E-01_wp, 5.05397955998414E-21_wp, 3.32539639889862E-03_wp, &
      &-2.05686172502982E-17_wp, 1.63311120018454E-03_wp,-9.04724595991492E-18_wp, &
      & 3.55558297044844E-05_wp,-8.03466743686323E-18_wp,-5.83217194169540E-06_wp, &
      & 2.37034465968426E-18_wp, 2.90073499286456E-09_wp,-2.59905923013698E-19_wp, &
      &-5.53998399309322E-07_wp, 3.61588930516630E-19_wp, 6.30517483271401E-10_wp, &
      &-6.31886225735169E-06_wp,-5.21771524622812E-18_wp, 3.23382948783783E-08_wp, &
      & 5.86797813348608E-18_wp,-6.42748734557216E-10_wp, 4.90163559995298E-19_wp, &
      & 3.10557638591897E-07_wp,-7.58835352038403E-20_wp,-1.39710906195764E-10_wp, &
      &-2.61461507468549E-01_wp,-1.63314702479790E-17_wp,-8.46120361482148E-02_wp, &
      &-3.91896272003785E-17_wp,-2.90703197794302E-05_wp,-5.10605669701393E-18_wp, &
      & 1.18677054674651E-02_wp, 1.38092961541038E-16_wp,-6.31886225735169E-06_wp, &
      & 2.01512275101693E+00_wp,-7.42728956353780E-17_wp, 3.99929220910119E-01_wp, &
      &-5.17302407382590E-17_wp, 1.63576811366029E-04_wp, 2.22760462887619E-18_wp, &
      &-3.51352591933450E-02_wp,-6.73707475021705E-18_wp, 3.55558297043036E-05_wp, &
      &-3.43182042645543E-17_wp, 3.34748212919108E-01_wp, 1.09460039181865E-16_wp, &
      &-5.90366231871004E-01_wp,-2.89827131754869E-18_wp, 8.76093662322536E-02_wp, &
      & 8.66539900807442E-17_wp, 5.29272140403641E-02_wp,-5.21771524622812E-18_wp, &
      &-7.42728956353780E-17_wp, 1.28272149679125E+00_wp,-1.61997535865294E-16_wp, &
      & 5.85799207483420E-01_wp, 1.09347353167198E-20_wp,-1.00951275464804E-02_wp, &
      &-1.38900521802473E-16_wp, 4.32459272824377E-02_wp, 4.98016653761617E-17_wp, &
      & 8.46120361482158E-02_wp, 7.57576586626171E-17_wp,-7.72403994815329E-01_wp, &
      & 3.60985123239876E-16_wp, 1.48774341197204E-07_wp,-8.07837585297486E-18_wp, &
      &-4.74250301651592E-02_wp, 1.63344576636165E-16_wp, 3.23382948783783E-08_wp, &
      & 3.99929220910119E-01_wp,-1.61997535865294E-16_wp, 8.16076846648328E-01_wp, &
      & 2.21926630850298E-16_wp, 2.68312706375597E-05_wp,-7.36781669749971E-18_wp, &
      &-5.39351527038124E-02_wp, 1.07528213863783E-17_wp, 5.83217194156254E-06_wp, &
      & 1.67783359907526E-16_wp,-5.90366231871004E-01_wp, 3.60887054888469E-17_wp, &
      & 7.80986178040750E-02_wp,-2.68991192418201E-18_wp, 5.29272140403642E-02_wp, &
      & 4.93091128376413E-17_wp, 1.10618385584490E-01_wp, 5.86797813348608E-18_wp, &
      &-5.17302407382590E-17_wp, 5.85799207483420E-01_wp, 2.21926630850298E-16_wp, &
      & 1.53738567184553E+00_wp,-6.02339688103840E-20_wp, 4.32459272824376E-02_wp, &
      &-1.45993073254721E-16_wp, 8.70515120970158E-03_wp,-1.47656869356249E-16_wp, &
      &-2.90703197793973E-05_wp, 2.30216039812600E-19_wp,-1.48774341268614E-07_wp, &
      & 2.08615446024900E-19_wp,-2.95700562701623E-09_wp,-1.19981759421529E-21_wp, &
      & 1.42873977881271E-06_wp, 5.05397955998414E-21_wp,-6.42748734557216E-10_wp, &
      & 1.63576811366029E-04_wp, 1.09347353167198E-20_wp, 2.68312706375597E-05_wp, &
      &-6.02339688103840E-20_wp, 1.33450121877156E-08_wp,-1.65913815044373E-20_wp, &
      &-2.54870417630752E-06_wp,-1.74783318955294E-20_wp, 2.90073499285004E-09_wp, &
      &-6.89667265911756E-18_wp,-8.76093662322534E-02_wp,-9.62026471439360E-21_wp, &
      &-5.29272140403639E-02_wp,-1.80963233126845E-19_wp, 1.87463408525633E-04_wp, &
      &-2.89559933915559E-18_wp, 3.32539639889862E-03_wp, 4.90163559995298E-19_wp, &
      & 2.22760462887619E-18_wp,-1.00951275464804E-02_wp,-7.36781669749971E-18_wp, &
      & 4.32459272824376E-02_wp,-1.65913815044373E-20_wp, 6.15075023236151E-03_wp, &
      & 9.22103308043608E-18_wp, 4.53852546312473E-03_wp,-2.44267976183480E-18_wp, &
      & 1.18677054674651E-02_wp,-1.39662490020440E-16_wp, 4.74250301651590E-02_wp, &
      &-1.37889229695067E-16_wp, 1.42873977881777E-06_wp,-9.45961822144709E-18_wp, &
      & 2.76687919541808E-03_wp,-2.05686172502982E-17_wp, 3.10557638591897E-07_wp, &
      &-3.51352591933450E-02_wp,-1.38900521802473E-16_wp,-5.39351527038124E-02_wp, &
      &-1.45993073254721E-16_wp,-2.54870417630752E-06_wp, 9.22103308043608E-18_wp, &
      & 3.73731491668372E-03_wp, 7.37860220244058E-18_wp,-5.53998399300359E-07_wp, &
      & 4.56621850648170E-18_wp,-5.29272140403637E-02_wp,-4.22680717602266E-17_wp, &
      &-1.10618385584490E-01_wp,-2.73041420985304E-19_wp, 3.32539639889863E-03_wp, &
      &-3.62779089330793E-18_wp, 1.63311120018454E-03_wp,-7.58835352038403E-20_wp, &
      &-6.73707475021705E-18_wp, 4.32459272824377E-02_wp, 1.07528213863783E-17_wp, &
      & 8.70515120970158E-03_wp,-1.74783318955294E-20_wp, 4.53852546312473E-03_wp, &
      & 7.37860220244058E-18_wp, 8.12378097099087E-03_wp, 8.58405190407147E-18_wp, &
      &-6.31886225739145E-06_wp, 5.43883875213128E-17_wp,-3.23382947650264E-08_wp, &
      &-1.41966749806042E-16_wp,-6.42748734559565E-10_wp, 2.06766115803617E-18_wp, &
      & 3.10557638597950E-07_wp,-9.04724595991492E-18_wp,-1.39710906195764E-10_wp, &
      & 3.55558297043036E-05_wp, 4.98016653761617E-17_wp, 5.83217194156254E-06_wp, &
      &-1.47656869356249E-16_wp, 2.90073499285004E-09_wp,-2.44267976183480E-18_wp, &
      &-5.53998399300359E-07_wp, 8.58405190407147E-18_wp, 6.30517483265413E-10_wp],&
      & shape(density))

   call get_structure(mol, "MB16-43", "S2")
   call test_energy_generic(error, mol, density, qsh, make_exchange_gxtb, &
      & -1.49369725174462_wp, thr_in=thr1)

end subroutine test_e_fock_s2


subroutine test_e_fock_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(7) = [&
      & 2.93108671070277E-01_wp,-1.70086400930251E-01_wp, 1.76661066979063E-01_wp, &
      &-7.49208342927936E-02_wp,-7.49208342927947E-02_wp,-7.49208342927952E-02_wp, &
      &-7.49208342927936E-02_wp]

   real(wp), parameter :: density(13, 13, 1) = reshape([&
      & 6.67076530840565E-01_wp, 1.37047397649895E-15_wp, 5.78064191406667E-16_wp, &
      &-1.08292025440329E-15_wp,-7.47124145130203E-18_wp, 1.16250868196032E-17_wp, &
      & 1.30914186280549E-17_wp,-1.45818112027175E-16_wp,-1.46951885454401E-18_wp, &
      & 2.51982328297639E-01_wp, 2.51982328297639E-01_wp, 2.51982328297639E-01_wp, &
      & 2.51982328297643E-01_wp, 1.37047397649895E-15_wp, 4.65632859964007E-01_wp, &
      & 8.08916089580512E-17_wp, 4.35134767335844E-16_wp, 4.62236407280884E-16_wp, &
      &-3.89696638482398E-16_wp,-1.34208728185462E-17_wp,-5.58631007293306E-02_wp, &
      & 3.42923026441888E-17_wp, 3.05808254360765E-01_wp,-3.05808254360767E-01_wp, &
      &-3.05808254360767E-01_wp, 3.05808254360766E-01_wp, 5.78064191406667E-16_wp, &
      & 8.08916089580512E-17_wp, 4.65632859964007E-01_wp,-6.93141114274268E-16_wp, &
      &-5.58631007293306E-02_wp,-6.20989199929278E-16_wp, 4.15547922215721E-17_wp, &
      & 5.02288208634189E-16_wp, 1.88740484982810E-18_wp,-3.05808254360767E-01_wp, &
      &-3.05808254360767E-01_wp, 3.05808254360765E-01_wp, 3.05808254360766E-01_wp, &
      &-1.08292025440329E-15_wp, 4.35134767335844E-16_wp,-6.93141114274268E-16_wp, &
      & 4.65632859964007E-01_wp,-3.65092803957423E-16_wp,-5.58631007293305E-02_wp, &
      & 2.39535445984415E-17_wp,-5.13779029809840E-16_wp, 4.04535013407324E-18_wp, &
      & 3.05808254360767E-01_wp,-3.05808254360766E-01_wp, 3.05808254360767E-01_wp, &
      &-3.05808254360766E-01_wp,-7.47124145130203E-18_wp, 4.62236407280884E-16_wp, &
      &-5.58631007293306E-02_wp,-3.65092803957423E-16_wp, 6.70203134576148E-03_wp, &
      & 1.28201004478432E-16_wp,-4.98542895756873E-18_wp,-1.17136813792813E-16_wp, &
      &-2.26436526088623E-19_wp, 3.66885561266806E-02_wp, 3.66885561266806E-02_wp, &
      &-3.66885561266809E-02_wp,-3.66885561266799E-02_wp, 1.16250868196032E-17_wp, &
      &-3.89696638482398E-16_wp,-6.20989199929278E-16_wp,-5.58631007293305E-02_wp, &
      & 1.28201004478432E-16_wp, 6.70203134576145E-03_wp,-2.87376469699906E-18_wp, &
      & 1.01955127023450E-16_wp,-4.85330442620878E-19_wp,-3.66885561266803E-02_wp, &
      & 3.66885561266810E-02_wp,-3.66885561266808E-02_wp, 3.66885561266797E-02_wp, &
      & 1.30914186280549E-17_wp,-1.34208728185462E-17_wp, 4.15547922215721E-17_wp, &
      & 2.39535445984415E-17_wp,-4.98542895756873E-18_wp,-2.87376469699906E-18_wp, &
      & 5.58449254396715E-33_wp, 1.61013458155841E-18_wp,-6.40697937995635E-34_wp, &
      &-1.54288667231892E-17_wp,-2.92637063477538E-17_wp, 5.67825839873141E-17_wp, &
      & 7.69066487781186E-18_wp,-1.45818112027175E-16_wp,-5.58631007293306E-02_wp, &
      & 5.02288208634189E-16_wp,-5.13779029809840E-16_wp,-1.17136813792813E-16_wp, &
      & 1.01955127023450E-16_wp, 1.61013458155841E-18_wp, 6.70203134576147E-03_wp, &
      &-4.11413051261264E-18_wp,-3.66885561266810E-02_wp, 3.66885561266805E-02_wp, &
      & 3.66885561266806E-02_wp,-3.66885561266798E-02_wp,-1.46951885454401E-18_wp, &
      & 3.42923026441888E-17_wp, 1.88740484982810E-18_wp, 4.04535013407324E-18_wp, &
      &-2.26436526088623E-19_wp,-4.85330442620878E-19_wp,-6.40697937995635E-34_wp, &
      &-4.11413051261264E-18_wp, 2.57154647550885E-33_wp, 2.33839056534641E-17_wp, &
      &-2.69732391292976E-17_wp,-1.91804667492846E-17_wp, 2.05494086189032E-17_wp, &
      & 2.51982328297639E-01_wp, 3.05808254360765E-01_wp,-3.05808254360767E-01_wp, &
      & 3.05808254360767E-01_wp, 3.66885561266806E-02_wp,-3.66885561266803E-02_wp, &
      &-1.54288667231892E-17_wp,-3.66885561266810E-02_wp, 2.33839056534641E-17_wp, &
      & 6.97710466990056E-01_wp,-1.05657993237668E-01_wp,-1.05657993237667E-01_wp, &
      &-1.05657993237669E-01_wp, 2.51982328297639E-01_wp,-3.05808254360767E-01_wp, &
      &-3.05808254360767E-01_wp,-3.05808254360766E-01_wp, 3.66885561266806E-02_wp, &
      & 3.66885561266810E-02_wp,-2.92637063477538E-17_wp, 3.66885561266805E-02_wp, &
      &-2.69732391292976E-17_wp,-1.05657993237668E-01_wp, 6.97710466990057E-01_wp, &
      &-1.05657993237668E-01_wp,-1.05657993237669E-01_wp, 2.51982328297639E-01_wp, &
      &-3.05808254360767E-01_wp, 3.05808254360765E-01_wp, 3.05808254360767E-01_wp, &
      &-3.66885561266809E-02_wp,-3.66885561266808E-02_wp, 5.67825839873141E-17_wp, &
      & 3.66885561266806E-02_wp,-1.91804667492846E-17_wp,-1.05657993237667E-01_wp, &
      &-1.05657993237668E-01_wp, 6.97710466990057E-01_wp,-1.05657993237669E-01_wp, &
      & 2.51982328297643E-01_wp, 3.05808254360766E-01_wp, 3.05808254360766E-01_wp, &
      &-3.05808254360766E-01_wp,-3.66885561266799E-02_wp, 3.66885561266797E-02_wp, &
      & 7.69066487781186E-18_wp,-3.66885561266798E-02_wp, 2.05494086189032E-17_wp, &
      &-1.05657993237669E-01_wp,-1.05657993237669E-01_wp,-1.05657993237669E-01_wp, &
      & 6.97710466990055E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "SiH4")
   call test_energy_generic(error, mol, density, qsh, make_exchange_gxtb, &
      & -0.781243149455857_wp, thr_in=thr1)

end subroutine test_e_fock_sih4


subroutine test_e_fock_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(13) = [&
      & 8.56111630751516E-01_wp, 1.08157064201355E-01_wp, 8.39880035447105E-01_wp, &
      &-4.26579824418720E-02_wp,-9.50257569167443E-02_wp,-5.37611462924087E-01_wp, &
      & 4.54632448304181E-02_wp,-9.62011525820314E-02_wp,-5.37310637136599E-01_wp, &
      & 4.55474939339988E-02_wp,-9.54523782951797E-02_wp,-5.36431940022193E-01_wp, &
      & 4.55318410022139E-02_wp]

   real(wp), parameter :: density(43, 43, 2) = reshape([&
      & 6.91508355650380E-03_wp,-4.64128438601378E-04_wp, 7.39852004665911E-04_wp, &
      &-3.08180805700232E-04_wp,-2.70834625071592E-03_wp, 6.73888290097408E-03_wp, &
      &-5.28983297438339E-03_wp, 4.82025419519933E-03_wp, 9.20948535364543E-04_wp, &
      &-2.62666642034588E-03_wp,-1.03344462572843E-02_wp,-2.13948643567722E-03_wp, &
      &-2.86176244725240E-03_wp,-2.19360566596758E-03_wp, 3.65330755839363E-03_wp, &
      & 1.10745388230377E-02_wp,-8.43277221798517E-04_wp, 3.10522969694570E-02_wp, &
      & 2.30387166634654E-02_wp, 2.04828236426981E-02_wp, 1.17529474576173E-03_wp, &
      & 1.12189820768851E-03_wp,-3.37676035843328E-04_wp, 7.56850448846645E-04_wp, &
      &-4.90010925833026E-04_wp,-1.05404313350324E-03_wp,-3.44887507450498E-02_wp, &
      &-1.34351050524750E-02_wp, 2.17894276273119E-02_wp,-1.32317745313018E-03_wp, &
      & 9.23295621268477E-04_wp,-5.90066147083195E-04_wp,-6.04014419841712E-04_wp, &
      &-5.29940796996229E-04_wp,-9.36839279449376E-04_wp, 9.73171638458462E-03_wp, &
      &-1.32361089963587E-02_wp,-4.00527607361329E-02_wp,-7.26493927890594E-04_wp, &
      &-2.84864369040110E-04_wp,-5.84337860555941E-04_wp, 1.06457919323247E-03_wp, &
      & 1.20863238120399E-03_wp,-4.64128438601378E-04_wp, 2.35817379836776E-03_wp, &
      & 3.92587064845560E-04_wp,-1.97054886034146E-04_wp,-2.70852106214767E-03_wp, &
      & 7.57824105301747E-04_wp, 3.58485595994388E-05_wp,-1.71837630862046E-03_wp, &
      &-1.11631208345474E-04_wp, 5.43641884415014E-04_wp, 1.70282647270043E-04_wp, &
      & 3.32527173464662E-03_wp,-3.25193570498183E-03_wp, 6.16981129415339E-03_wp, &
      &-2.47480280545290E-03_wp,-2.24742086894289E-04_wp, 6.95417545032155E-03_wp, &
      &-1.08633249907778E-02_wp,-1.76598589656093E-02_wp,-1.93609624730210E-02_wp, &
      &-8.72287324457150E-04_wp,-5.65631717526762E-04_wp, 1.93375364112228E-04_wp, &
      &-6.13369844793697E-04_wp,-7.09402995502100E-05_wp,-7.79253610457477E-03_wp, &
      &-1.13057960252247E-02_wp,-1.77696559058022E-02_wp, 2.23419903016915E-02_wp, &
      &-8.58013769944745E-04_wp, 6.57139878410130E-04_wp,-1.13106310747154E-04_wp, &
      &-7.25066154483291E-04_wp, 1.72200355417520E-04_wp, 2.47719305207997E-03_wp, &
      & 1.50974537942946E-02_wp, 4.19681624623315E-03_wp, 1.19489455972230E-02_wp, &
      &-4.06772228425090E-04_wp,-1.49187995855161E-04_wp,-1.57324892220013E-05_wp, &
      &-4.42293634056734E-04_wp,-5.24582933230386E-04_wp, 7.39852004665911E-04_wp, &
      & 3.92587064845560E-04_wp, 1.91536935819115E-03_wp, 2.47696444538848E-04_wp, &
      &-1.35580198876781E-03_wp,-1.33623843632451E-03_wp, 2.74713667604232E-04_wp, &
      &-5.62965067558568E-04_wp, 4.27221963605903E-04_wp, 1.26154191937191E-03_wp, &
      & 2.11079460969478E-03_wp,-7.92506392031096E-03_wp,-3.74932686039594E-03_wp, &
      &-9.22902125332738E-03_wp,-1.87861138540621E-03_wp, 2.91424187960775E-03_wp, &
      & 4.48629748476077E-03_wp,-1.94279763315955E-02_wp, 1.08495767742708E-02_wp, &
      &-1.31613586515294E-02_wp,-6.31970749173562E-04_wp,-1.51926912208397E-04_wp, &
      & 4.76567196308843E-04_wp,-6.46992100973259E-05_wp, 2.63809170113323E-04_wp, &
      &-3.70243865248362E-03_wp,-1.60655875718359E-02_wp, 1.10901118217271E-02_wp, &
      & 1.06050858375466E-02_wp,-6.86781164787524E-04_wp,-1.07485745623368E-05_wp, &
      &-5.80144587727174E-04_wp,-3.80504491466814E-05_wp,-1.88723311255307E-04_wp, &
      &-3.49727598971794E-03_wp, 4.62053531076067E-03_wp, 1.19028969178129E-02_wp, &
      &-1.89336042467516E-02_wp,-4.11421365806423E-04_wp,-8.15783274082982E-05_wp, &
      &-5.09349489441788E-04_wp, 4.03617862685355E-05_wp, 4.69679477374292E-04_wp, &
      &-3.08180805700232E-04_wp,-1.97054886034146E-04_wp, 2.47696444538848E-04_wp, &
      & 2.47716089856009E-03_wp,-1.15071806995931E-03_wp,-1.49842966722536E-03_wp, &
      &-1.80625067325603E-06_wp, 1.76900754059816E-03_wp, 2.09637797823798E-03_wp, &
      &-2.83985293835542E-03_wp, 1.88827606958926E-03_wp, 4.80174630401615E-03_wp, &
      &-2.23540115944937E-03_wp, 1.98252298932427E-03_wp, 3.75687033823298E-03_wp, &
      &-1.16227766860251E-03_wp, 4.60019091474426E-03_wp,-1.97342698427398E-02_wp, &
      &-1.21032084322852E-02_wp, 5.11250873715631E-03_wp,-2.83717722384947E-04_wp, &
      &-5.88103188125038E-04_wp, 1.40865283219372E-04_wp,-9.26782496704140E-05_wp, &
      & 7.36601934773917E-04_wp, 5.33025120676840E-03_wp, 2.09425905656062E-02_wp, &
      & 1.10099041175004E-02_wp, 3.15743661071952E-03_wp, 3.35355565772707E-04_wp, &
      &-7.16937802389864E-04_wp, 5.04915814179606E-05_wp, 2.25742379909139E-04_wp, &
      & 6.87321959615046E-04_wp,-8.83698321844673E-03_wp, 1.35558351291241E-02_wp, &
      &-2.03523592200911E-02_wp,-2.27896269355215E-02_wp,-6.80990118806957E-04_wp, &
      &-4.08274548452866E-04_wp,-1.56609917648878E-04_wp, 1.00472108720785E-03_wp, &
      & 6.42113030611923E-04_wp,-2.70834625071592E-03_wp,-2.70852106214767E-03_wp, &
      &-1.35580198876781E-03_wp,-1.15071806995931E-03_wp, 4.40698343138236E-02_wp, &
      & 5.53685140184838E-03_wp, 2.57661510240173E-03_wp, 7.74565338533846E-03_wp, &
      &-4.29167996159675E-04_wp, 8.87283046251010E-04_wp, 7.59817981237967E-03_wp, &
      &-2.90121327491061E-03_wp, 1.65065913001023E-02_wp,-8.56320568075693E-03_wp, &
      &-9.61743773202385E-04_wp, 2.03348322954967E-03_wp, 2.24529453955674E-02_wp, &
      & 6.64065847829613E-02_wp, 8.88285708892675E-02_wp, 6.09778408701114E-04_wp, &
      & 1.00856122335875E-03_wp, 2.43202433037870E-03_wp, 2.80111828722630E-04_wp, &
      & 1.28797404654250E-03_wp,-1.29393581855696E-03_wp,-2.40535817565888E-02_wp, &
      & 7.93344293670733E-02_wp, 7.02584780993291E-02_wp,-1.26549176620155E-02_wp, &
      & 1.71977080519939E-03_wp,-2.33627976151097E-03_wp, 1.31736482134134E-04_wp, &
      & 1.23882021027591E-03_wp, 1.31126864067533E-03_wp,-1.35887092213168E-02_wp, &
      & 7.65773113004499E-02_wp, 3.30174684587050E-02_wp, 6.43422160137272E-02_wp, &
      &-1.00280886934538E-03_wp,-7.03926923327412E-04_wp, 3.69301350068986E-04_wp, &
      &-1.38996067520150E-03_wp,-2.37514633952885E-03_wp, 6.73888290097408E-03_wp, &
      & 7.57824105301747E-04_wp,-1.33623843632451E-03_wp,-1.49842966722536E-03_wp, &
      & 5.53685140184838E-03_wp, 3.15062504394876E-02_wp, 5.47450080184616E-03_wp, &
      &-5.00796852925626E-03_wp,-7.28454124897730E-03_wp, 2.79751133190821E-03_wp, &
      &-1.71486712002216E-02_wp, 7.94980209958558E-03_wp, 1.65027427287929E-03_wp, &
      & 3.93589352195860E-03_wp, 8.76263384418473E-04_wp, 9.74759706412966E-03_wp, &
      & 2.07568903488779E-02_wp, 6.76261657530245E-02_wp,-2.17931462181147E-02_wp, &
      & 7.90366362305845E-02_wp, 2.79402928154148E-03_wp, 6.30313400095883E-04_wp, &
      &-1.78080710395585E-03_wp, 6.08343459325308E-04_wp,-4.94440667285671E-04_wp, &
      & 1.84723481845835E-02_wp,-6.44273893756731E-02_wp, 3.79307091876071E-02_wp, &
      & 6.86738430361502E-02_wp,-2.66785048528304E-03_wp, 3.39078065070956E-04_wp, &
      &-2.13068231597136E-03_wp,-3.19054602405769E-04_wp,-6.05877111680159E-04_wp, &
      &-5.22575499812361E-03_wp, 3.76217969547401E-02_wp,-1.71192760195093E-02_wp, &
      & 4.37862349890877E-02_wp,-7.45923917349099E-04_wp,-5.93743903736391E-04_wp, &
      & 5.99040723188793E-04_wp, 9.36631352457798E-05_wp,-1.11597837195792E-03_wp, &
      &-5.28983297438339E-03_wp, 3.58485595994388E-05_wp, 2.74713667604232E-04_wp, &
      &-1.80625067325603E-06_wp, 2.57661510240173E-03_wp, 5.47450080184616E-03_wp, &
      & 2.93992920871206E-02_wp, 3.63929600632670E-03_wp,-7.16356695186403E-04_wp, &
      & 3.39233717806562E-03_wp,-1.74695534855078E-04_wp, 3.76998134788861E-03_wp, &
      &-1.42151291041326E-05_wp, 5.13883241660874E-03_wp, 6.46716131991330E-04_wp, &
      &-1.71066964711473E-02_wp,-7.55036512350748E-03_wp, 1.56196228580328E-02_wp, &
      &-9.66954644874436E-02_wp, 1.11349817585917E-02_wp, 1.22603152276987E-04_wp, &
      &-2.11297913066481E-03_wp,-1.37122114854505E-03_wp,-1.46955760523574E-03_wp, &
      &-3.20358396726471E-05_wp,-9.72011490896808E-03_wp, 3.77933661329130E-04_wp, &
      & 8.63683672373309E-02_wp,-5.52544574468794E-03_wp, 2.37974672043536E-04_wp, &
      &-1.95351899121100E-03_wp,-1.29213626521191E-03_wp, 1.69539224395771E-03_wp, &
      &-2.70533226346614E-04_wp,-9.91706868410338E-03_wp,-5.51719337384838E-03_wp, &
      & 8.65390101840318E-02_wp, 3.02003084647696E-03_wp, 3.82254327745874E-04_wp, &
      & 9.48495944206912E-04_wp,-1.24471774247931E-03_wp,-2.38179137703504E-03_wp, &
      & 4.90468431720295E-05_wp, 4.82025419519933E-03_wp,-1.71837630862046E-03_wp, &
      &-5.62965067558568E-04_wp, 1.76900754059816E-03_wp, 7.74565338533846E-03_wp, &
      &-5.00796852925626E-03_wp, 3.63929600632670E-03_wp, 3.44559738728420E-02_wp, &
      & 4.73877555844598E-03_wp,-1.08018636242987E-02_wp,-4.09589429070970E-03_wp, &
      &-1.24378156087350E-02_wp, 8.23509093532584E-04_wp,-3.00732828881522E-03_wp, &
      & 1.05004767401701E-02_wp, 5.68415105020803E-03_wp, 1.41430406380897E-02_wp, &
      & 7.87648653119079E-02_wp,-1.31132079427450E-02_wp, 5.05353428353156E-03_wp, &
      & 1.56860526279678E-03_wp, 4.95050490542593E-04_wp,-1.40540921397035E-03_wp, &
      &-4.36565995938913E-06_wp,-1.55186627714891E-03_wp,-1.19396614153255E-02_wp, &
      & 7.87414432104843E-02_wp,-2.99966874007633E-02_wp,-2.09673932343892E-03_wp, &
      & 1.14415800389318E-03_wp,-1.63860594224303E-04_wp, 1.53063969395156E-03_wp, &
      &-6.84469334258254E-04_wp, 1.76370008603692E-03_wp, 2.14468953008398E-02_wp, &
      & 3.28265258373041E-02_wp, 4.80550657744445E-02_wp,-1.04816257994469E-01_wp, &
      &-1.57983395972902E-03_wp,-8.84185331974615E-05_wp,-2.38428284607014E-03_wp, &
      & 4.96040206857101E-04_wp, 2.58029772016527E-03_wp, 9.20948535364543E-04_wp, &
      &-1.11631208345474E-04_wp, 4.27221963605903E-04_wp, 2.09637797823798E-03_wp, &
      &-4.29167996159675E-04_wp,-7.28454124897730E-03_wp,-7.16356695186403E-04_wp, &
      & 4.73877555844598E-03_wp, 4.33658766595577E-02_wp,-3.56786230735714E-03_wp, &
      &-6.98090014751096E-04_wp, 7.94216229028011E-03_wp,-6.12349300091585E-03_wp, &
      &-1.95159298365422E-03_wp, 7.85741657894896E-03_wp,-1.75117606699352E-03_wp, &
      &-8.86026264798512E-03_wp, 3.47893280541423E-02_wp,-3.47126886510459E-02_wp, &
      &-8.92145930930235E-02_wp,-1.31598285372499E-03_wp,-3.50665892391665E-04_wp, &
      &-9.20931256837223E-05_wp,-1.31550509938713E-03_wp,-1.78591283773075E-03_wp, &
      &-8.62500600176080E-03_wp,-3.79155045233419E-02_wp, 3.52606751298176E-02_wp, &
      &-9.70310332800300E-02_wp, 1.09295237553953E-03_wp,-3.65714841904203E-04_wp, &
      &-3.71310081988481E-04_wp, 1.86764564795742E-03_wp,-2.30237597503765E-03_wp, &
      & 2.20223984658526E-02_wp, 7.33725870090370E-02_wp,-7.21093493716756E-02_wp, &
      &-4.68711668322494E-02_wp,-2.11281803772875E-03_wp,-1.44674369270472E-03_wp, &
      & 1.06345967457132E-04_wp, 2.30301477281955E-03_wp, 5.27801603506582E-04_wp, &
      &-2.62666642034588E-03_wp, 5.43641884415014E-04_wp, 1.26154191937191E-03_wp, &
      &-2.83985293835542E-03_wp, 8.87283046251010E-04_wp, 2.79751133190821E-03_wp, &
      & 3.39233717806562E-03_wp,-1.08018636242987E-02_wp,-3.56786230735714E-03_wp, &
      & 2.93613172945072E-02_wp, 2.31394574059330E-02_wp,-6.60986023806925E-02_wp, &
      & 1.25016127644544E-03_wp,-6.24509007128584E-02_wp,-1.03216597964290E-02_wp, &
      & 1.32644245444874E-02_wp, 4.71230456570157E-04_wp,-4.68879647252080E-02_wp, &
      &-2.47262717241474E-03_wp, 4.41570213412796E-02_wp, 1.02316183175296E-03_wp, &
      &-1.82803596260942E-03_wp,-8.23217732713214E-04_wp, 9.93031268178490E-05_wp, &
      & 2.01374750556590E-03_wp, 3.07722003547285E-04_wp,-5.84884648956699E-02_wp, &
      &-1.64182145464190E-02_wp,-5.32623004110966E-02_wp, 5.34718902385395E-04_wp, &
      & 1.82873529738723E-03_wp, 5.62884466208868E-04_wp, 1.47717103400572E-05_wp, &
      &-2.63392884744990E-03_wp,-1.29176440695309E-04_wp, 3.91140770304674E-02_wp, &
      & 5.34660261424884E-02_wp, 5.84896470057643E-02_wp,-1.07874841815666E-03_wp, &
      &-7.53053375859206E-04_wp, 3.71134001508957E-04_wp,-2.11151814934020E-03_wp, &
      &-2.94550580024158E-03_wp,-1.03344462572843E-02_wp, 1.70282647270043E-04_wp, &
      & 2.11079460969478E-03_wp, 1.88827606958926E-03_wp, 7.59817981237967E-03_wp, &
      &-1.71486712002216E-02_wp,-1.74695534855078E-04_wp,-4.09589429070970E-03_wp, &
      &-6.98090014751096E-04_wp, 2.31394574059330E-02_wp, 5.91416085295369E-02_wp, &
      &-1.22441986810700E-01_wp, 9.13592493676379E-04_wp,-1.30860871682892E-01_wp, &
      &-2.79634534925338E-02_wp, 1.48760975395522E-02_wp, 3.28838948018154E-04_wp, &
      &-9.20017100476902E-02_wp,-9.28727236293850E-03_wp,-3.36526330392563E-02_wp, &
      &-1.36605831419777E-03_wp,-4.04469794704718E-03_wp,-5.61547036728571E-04_wp, &
      &-1.79865124144556E-03_wp, 1.74946442106978E-03_wp, 1.11491519164976E-04_wp, &
      & 7.63322292457885E-02_wp, 4.93084007884058E-03_wp,-3.14965757562480E-02_wp, &
      & 2.32955790448110E-03_wp,-9.71790300967990E-04_wp, 2.00789686792316E-03_wp, &
      &-3.47566238648489E-04_wp, 1.90671492659640E-03_wp, 8.63738096668437E-04_wp, &
      & 2.78573450072921E-02_wp, 8.43233407925204E-03_wp, 4.10236255981251E-02_wp, &
      &-7.73918178641165E-04_wp,-1.64698024921482E-03_wp, 1.82838363631060E-03_wp, &
      & 3.86732721439976E-04_wp,-3.09148435784959E-03_wp,-2.13948643567722E-03_wp, &
      & 3.32527173464662E-03_wp,-7.92506392031096E-03_wp, 4.80174630401615E-03_wp, &
      &-2.90121327491061E-03_wp, 7.94980209958558E-03_wp, 3.76998134788861E-03_wp, &
      &-1.24378156087350E-02_wp, 7.94216229028011E-03_wp,-6.60986023806925E-02_wp, &
      &-1.22441986810700E-01_wp, 4.33564582130013E-01_wp, 1.95611504332125E-03_wp, &
      & 4.51515755881632E-01_wp, 8.78258948759680E-02_wp,-9.51362732512129E-02_wp, &
      &-1.25508317134303E-03_wp,-7.63128539187115E-03_wp,-1.11533246612948E-02_wp, &
      &-2.66711934094889E-04_wp,-4.41947865520944E-03_wp, 6.85882517062905E-03_wp, &
      & 7.11134130767397E-03_wp, 3.37358828176490E-03_wp, 1.44852013817612E-03_wp, &
      & 5.29641913037380E-04_wp, 8.80447817429636E-03_wp, 3.91330610474767E-02_wp, &
      & 9.30647517605180E-03_wp, 4.90411824161737E-05_wp,-2.82389512930550E-03_wp, &
      &-2.61623063595918E-03_wp, 3.78836597522281E-03_wp,-7.48604330870602E-04_wp, &
      &-1.70308189354695E-04_wp, 1.13610126762784E-02_wp,-6.50954576516091E-02_wp, &
      & 2.26688923785656E-02_wp, 6.83511061285916E-04_wp, 3.49269145418072E-03_wp, &
      &-3.31410275388691E-03_wp,-2.97731431290906E-03_wp, 4.21575971843549E-03_wp, &
      &-2.86176244725240E-03_wp,-3.25193570498183E-03_wp,-3.74932686039594E-03_wp, &
      &-2.23540115944937E-03_wp, 1.65065913001023E-02_wp, 1.65027427287929E-03_wp, &
      &-1.42151291041326E-05_wp, 8.23509093532584E-04_wp,-6.12349300091585E-03_wp, &
      & 1.25016127644544E-03_wp, 9.13592493676379E-04_wp, 1.95611504332125E-03_wp, &
      & 1.54119080114156E-02_wp, 2.26395335277116E-03_wp, 2.67540292426281E-04_wp, &
      &-4.53779864440760E-03_wp, 7.02122769989784E-04_wp, 5.25285844224189E-02_wp, &
      & 3.34567993236921E-02_wp, 3.58397598243777E-02_wp, 1.95121967815390E-03_wp, &
      & 1.52224078153395E-03_wp,-7.20199271675900E-04_wp, 1.03691084693786E-03_wp, &
      &-7.60692602493730E-04_wp,-3.00213735625271E-04_wp, 4.94121495627030E-02_wp, &
      & 6.09420749900072E-03_wp,-4.30860195133737E-02_wp, 2.32658753451363E-03_wp, &
      &-7.53508690225449E-04_wp, 1.24050701272183E-03_wp, 8.55022975789090E-04_wp, &
      & 1.67234209279600E-04_wp,-2.42765201916911E-04_wp,-2.43254052385456E-02_wp, &
      & 5.40323761942012E-03_wp, 6.11536812838240E-02_wp, 1.67995322470535E-03_wp, &
      & 5.55389209372992E-04_wp, 1.20747771046646E-03_wp,-9.68228710592495E-04_wp, &
      &-1.52225439381485E-03_wp,-2.19360566596758E-03_wp, 6.16981129415339E-03_wp, &
      &-9.22902125332738E-03_wp, 1.98252298932427E-03_wp,-8.56320568075693E-03_wp, &
      & 3.93589352195860E-03_wp, 5.13883241660874E-03_wp,-3.00732828881522E-03_wp, &
      &-1.95159298365422E-03_wp,-6.24509007128584E-02_wp,-1.30860871682892E-01_wp, &
      & 4.51515755881632E-01_wp, 2.26395335277116E-03_wp, 4.94207247301378E-01_wp, &
      & 9.42909644531336E-02_wp,-9.92766824619663E-02_wp,-1.38055776894079E-03_wp, &
      & 2.33078648947951E-03_wp,-3.34565199386503E-02_wp, 7.72600418280543E-03_wp, &
      &-4.36466827514398E-03_wp, 6.92260048146288E-03_wp, 7.01847980989513E-03_wp, &
      & 3.33671398085666E-03_wp, 1.34879934078008E-03_wp, 1.35487620190833E-03_wp, &
      & 8.94664814628892E-03_wp,-6.96974239270687E-02_wp, 1.83044694261429E-03_wp, &
      & 2.45059911759910E-04_wp, 1.87682772385266E-04_wp,-6.56177016050264E-04_wp, &
      & 1.97824728090476E-03_wp,-1.35163183416280E-03_wp,-1.19394373726397E-03_wp, &
      & 2.31845543892078E-02_wp, 2.60925461908099E-02_wp, 2.77813689934099E-02_wp, &
      & 2.20058704998442E-04_wp, 4.37711797391931E-03_wp,-5.12579723062006E-03_wp, &
      &-6.39253710685816E-03_wp, 3.81730219384280E-03_wp, 3.65330755839363E-03_wp, &
      &-2.47480280545290E-03_wp,-1.87861138540621E-03_wp, 3.75687033823298E-03_wp, &
      &-9.61743773202385E-04_wp, 8.76263384418473E-04_wp, 6.46716131991330E-04_wp, &
      & 1.05004767401701E-02_wp, 7.85741657894896E-03_wp,-1.03216597964290E-02_wp, &
      &-2.79634534925338E-02_wp, 8.78258948759680E-02_wp, 2.67540292426281E-04_wp, &
      & 9.42909644531336E-02_wp, 3.69236498248453E-02_wp,-1.38065077020077E-02_wp, &
      &-1.76974667651580E-04_wp,-3.51739862290911E-03_wp,-2.38855436262098E-03_wp, &
      & 7.69776474169230E-02_wp, 1.36643652078677E-03_wp, 1.37587195953166E-03_wp, &
      & 7.48778153956748E-04_wp, 2.13490629757382E-03_wp, 2.25243758392478E-03_wp, &
      &-5.92614759492204E-04_wp, 1.22050640255783E-02_wp,-5.67599294232971E-04_wp, &
      &-5.95871699567429E-02_wp, 1.86239809890842E-03_wp,-4.55777507061880E-04_wp, &
      & 4.50807312848272E-04_wp, 1.66783078270450E-03_wp,-1.27146849958339E-03_wp, &
      & 1.03573655632359E-04_wp, 4.75824479949787E-02_wp,-1.10491431748055E-02_wp, &
      &-6.98072919459784E-02_wp,-2.04606432020437E-03_wp,-4.78902053888293E-06_wp, &
      &-2.27508889527851E-03_wp, 4.87698483462043E-04_wp, 2.67025836021827E-03_wp, &
      & 1.10745388230377E-02_wp,-2.24742086894289E-04_wp, 2.91424187960775E-03_wp, &
      &-1.16227766860251E-03_wp, 2.03348322954967E-03_wp, 9.74759706412966E-03_wp, &
      &-1.71066964711473E-02_wp, 5.68415105020803E-03_wp,-1.75117606699352E-03_wp, &
      & 1.32644245444874E-02_wp, 1.48760975395522E-02_wp,-9.51362732512129E-02_wp, &
      &-4.53779864440760E-03_wp,-9.92766824619663E-02_wp,-1.38065077020077E-02_wp, &
      & 4.32230703620958E-02_wp,-4.67496303155928E-04_wp, 2.48294893113605E-02_wp, &
      & 7.10731854822242E-02_wp, 3.34201329272351E-02_wp, 2.32681980688723E-03_wp, &
      & 5.56817632099659E-04_wp,-1.25016037903597E-03_wp, 9.41199215808539E-04_wp, &
      &-3.02774896627336E-04_wp,-3.55591598093717E-04_wp,-4.93733275076534E-02_wp, &
      &-4.86973240502231E-02_wp, 4.52876588207916E-02_wp,-2.27875919841765E-03_wp, &
      & 2.50121224806569E-03_wp,-5.97769698513096E-05_wp,-2.37571013869545E-03_wp, &
      &-2.62237525358857E-04_wp,-2.86493133067578E-04_wp, 5.67692155285208E-02_wp, &
      &-2.51528185615128E-02_wp,-3.76025089855853E-02_wp,-2.46663081102699E-03_wp, &
      &-2.04537262522129E-03_wp, 2.56519694622059E-04_wp, 2.33412819707485E-03_wp, &
      &-2.89987015622134E-04_wp,-8.43277221798517E-04_wp, 6.95417545032155E-03_wp, &
      & 4.48629748476077E-03_wp, 4.60019091474426E-03_wp, 2.24529453955674E-02_wp, &
      & 2.07568903488779E-02_wp,-7.55036512350748E-03_wp, 1.41430406380897E-02_wp, &
      &-8.86026264798512E-03_wp, 4.71230456570157E-04_wp, 3.28838948018154E-04_wp, &
      &-1.25508317134303E-03_wp, 7.02122769989784E-04_wp,-1.38055776894079E-03_wp, &
      &-1.76974667651580E-04_wp,-4.67496303155928E-04_wp, 9.79337811380466E-01_wp, &
      &-6.74219435004162E-02_wp,-4.24289626423976E-02_wp,-4.55035139131935E-02_wp, &
      &-3.44483338939823E-03_wp,-3.19796194196180E-03_wp, 1.08310954300958E-03_wp, &
      &-2.17887584028734E-03_wp, 1.36316094026119E-03_wp,-1.74330915702325E-03_wp, &
      &-4.58677935890264E-03_wp,-4.64473531872142E-03_wp,-2.29331522657019E-02_wp, &
      & 1.08820023071544E-03_wp, 8.19971309071037E-04_wp, 8.60563453836057E-04_wp, &
      & 3.74866402148292E-04_wp,-2.27285691797216E-03_wp,-4.66588242543135E-04_wp, &
      &-2.18088403427284E-02_wp,-4.84923506617974E-03_wp, 2.45274151312103E-03_wp, &
      & 2.10826592111825E-03_wp, 5.15305178827397E-04_wp, 8.24138845319420E-04_wp, &
      & 7.13070047842620E-04_wp, 1.14992365340433E-03_wp, 3.10522969694570E-02_wp, &
      &-1.08633249907778E-02_wp,-1.94279763315955E-02_wp,-1.97342698427398E-02_wp, &
      & 6.64065847829613E-02_wp, 6.76261657530245E-02_wp, 1.56196228580328E-02_wp, &
      & 7.87648653119079E-02_wp, 3.47893280541423E-02_wp,-4.68879647252080E-02_wp, &
      &-9.20017100476902E-02_wp,-7.63128539187115E-03_wp, 5.25285844224189E-02_wp, &
      & 2.33078648947951E-03_wp,-3.51739862290911E-03_wp, 2.48294893113605E-02_wp, &
      &-6.74219435004162E-02_wp, 8.89616588229133E-01_wp,-1.30912197418645E-02_wp, &
      &-1.38437693771471E-02_wp, 1.53689319173918E-02_wp, 1.30549540811057E-02_wp, &
      &-1.30573252929538E-02_wp,-1.98932965996673E-05_wp,-2.41087735648506E-02_wp, &
      & 5.83917668822536E-03_wp, 4.54852027795121E-03_wp,-1.29884368229359E-02_wp, &
      &-1.20861841270074E-02_wp,-6.93166920400433E-04_wp, 2.31083568203780E-03_wp, &
      & 3.77274788892379E-04_wp,-9.37624772359046E-04_wp,-9.82783410965334E-04_wp, &
      &-1.79306724772635E-02_wp,-4.45043895094737E-02_wp,-4.51857519227536E-03_wp, &
      & 3.23050890198580E-02_wp, 5.49801719243796E-03_wp, 2.23738762999721E-03_wp, &
      & 2.14802109998526E-03_wp, 2.15131149462760E-04_wp,-2.63690796504470E-04_wp, &
      & 2.30387166634654E-02_wp,-1.76598589656093E-02_wp, 1.08495767742708E-02_wp, &
      &-1.21032084322852E-02_wp, 8.88285708892675E-02_wp,-2.17931462181147E-02_wp, &
      &-9.66954644874436E-02_wp,-1.31132079427450E-02_wp,-3.47126886510459E-02_wp, &
      &-2.47262717241474E-03_wp,-9.28727236293850E-03_wp,-1.11533246612948E-02_wp, &
      & 3.34567993236921E-02_wp,-3.34565199386503E-02_wp,-2.38855436262098E-03_wp, &
      & 7.10731854822242E-02_wp,-4.24289626423976E-02_wp,-1.30912197418645E-02_wp, &
      & 9.01706958053382E-01_wp,-9.56238625153878E-03_wp, 1.42821308760973E-03_wp, &
      & 2.13457270793031E-02_wp, 1.22148337982342E-02_wp, 1.49853968405253E-02_wp, &
      &-4.43859381082050E-04_wp, 9.69086766965422E-04_wp,-1.16033600518990E-02_wp, &
      & 2.04795523049283E-02_wp,-1.36569808889512E-02_wp,-7.21222341323190E-04_wp, &
      & 7.26268373831813E-05_wp,-1.52972790544857E-03_wp, 5.08502570020347E-04_wp, &
      &-2.09194354735716E-03_wp, 1.70652047126286E-03_wp,-1.64247670292941E-02_wp, &
      & 2.26152065790688E-02_wp,-9.39650646415531E-03_wp, 5.30549401718718E-04_wp, &
      & 1.21492153053422E-04_wp,-1.14204140349652E-03_wp, 3.16308765273205E-04_wp, &
      & 1.51716666848909E-03_wp, 2.04828236426981E-02_wp,-1.93609624730210E-02_wp, &
      &-1.31613586515294E-02_wp, 5.11250873715631E-03_wp, 6.09778408701114E-04_wp, &
      & 7.90366362305845E-02_wp, 1.11349817585917E-02_wp, 5.05353428353156E-03_wp, &
      &-8.92145930930235E-02_wp, 4.41570213412796E-02_wp,-3.36526330392563E-02_wp, &
      &-2.66711934094889E-04_wp, 3.58397598243777E-02_wp, 7.72600418280543E-03_wp, &
      & 7.69776474169230E-02_wp, 3.34201329272351E-02_wp,-4.55035139131935E-02_wp, &
      &-1.38437693771471E-02_wp,-9.56238625153878E-03_wp, 9.01322100189292E-01_wp, &
      & 2.36536807686828E-02_wp, 3.21986151735749E-04_wp,-8.65332883356591E-03_wp, &
      & 1.29451185710472E-02_wp, 1.70363857376580E-02_wp,-2.30363044077329E-02_wp, &
      & 1.65447430322963E-02_wp, 2.35203484160585E-03_wp,-5.05728174639121E-02_wp, &
      & 5.43932743796690E-03_wp,-9.04895911711816E-05_wp, 2.37278396459574E-03_wp, &
      & 2.62727701317609E-03_wp,-2.77814041002767E-03_wp, 1.34190851188152E-02_wp, &
      & 9.68094770971209E-03_wp,-1.32881963094551E-02_wp,-3.36509581492724E-03_wp, &
      &-2.43229181329884E-03_wp,-1.01164103858773E-03_wp,-7.48267748507518E-04_wp, &
      & 2.35043433033507E-03_wp, 1.60863886138561E-03_wp, 1.17529474576173E-03_wp, &
      &-8.72287324457150E-04_wp,-6.31970749173562E-04_wp,-2.83717722384947E-04_wp, &
      & 1.00856122335875E-03_wp, 2.79402928154148E-03_wp, 1.22603152276987E-04_wp, &
      & 1.56860526279678E-03_wp,-1.31598285372499E-03_wp, 1.02316183175296E-03_wp, &
      &-1.36605831419777E-03_wp,-4.41947865520944E-03_wp, 1.95121967815390E-03_wp, &
      &-4.36466827514398E-03_wp, 1.36643652078677E-03_wp, 2.32681980688723E-03_wp, &
      &-3.44483338939823E-03_wp, 1.53689319173918E-02_wp, 1.42821308760973E-03_wp, &
      & 2.36536807686828E-02_wp, 9.68828848521314E-04_wp, 2.22946911224688E-04_wp, &
      &-5.13301914431355E-04_wp, 3.54712418689991E-04_wp, 1.78411143916457E-05_wp, &
      &-9.90869743129443E-04_wp,-4.65550886867003E-04_wp,-1.08192538127679E-03_wp, &
      &-5.13523542037677E-03_wp, 2.00183991552993E-04_wp, 8.28467511050551E-05_wp, &
      & 1.24197069253257E-04_wp, 6.25152483639254E-05_wp,-1.66773246000592E-04_wp, &
      & 9.05681949044069E-04_wp,-1.85522927119767E-03_wp,-1.54510199955994E-03_wp, &
      &-8.07848283131484E-04_wp, 5.56541903522150E-05_wp,-1.75762882400090E-05_wp, &
      & 6.87763115435970E-05_wp, 1.51322622194485E-04_wp, 3.34775023659391E-05_wp, &
      & 1.12189820768851E-03_wp,-5.65631717526762E-04_wp,-1.51926912208397E-04_wp, &
      &-5.88103188125038E-04_wp, 2.43202433037870E-03_wp, 6.30313400095883E-04_wp, &
      &-2.11297913066481E-03_wp, 4.95050490542593E-04_wp,-3.50665892391665E-04_wp, &
      &-1.82803596260942E-03_wp,-4.04469794704718E-03_wp, 6.85882517062905E-03_wp, &
      & 1.52224078153395E-03_wp, 6.92260048146288E-03_wp, 1.37587195953166E-03_wp, &
      & 5.56817632099659E-04_wp,-3.19796194196180E-03_wp, 1.30549540811057E-02_wp, &
      & 2.13457270793031E-02_wp, 3.21986151735749E-04_wp, 2.22946911224688E-04_wp, &
      & 8.45466246656272E-04_wp, 2.05800937159412E-04_wp, 4.33613368926151E-04_wp, &
      &-3.53468297382543E-04_wp, 6.90995730459461E-04_wp,-2.30536700688280E-03_wp, &
      &-2.21987363481164E-04_wp,-6.65356782173023E-04_wp,-5.73474810553668E-05_wp, &
      & 4.56304467441363E-05_wp,-8.45754544145134E-05_wp, 5.05322897854274E-05_wp, &
      &-1.40434975631984E-04_wp,-3.66882904454643E-04_wp,-3.35446033168041E-03_wp, &
      &-3.70548351221492E-04_wp,-2.78556704937534E-04_wp, 1.62696960532809E-04_wp, &
      & 1.32878219836544E-04_wp,-6.24460321559444E-05_wp,-4.09463185306541E-05_wp, &
      & 1.51524804433381E-04_wp,-3.37676035843328E-04_wp, 1.93375364112228E-04_wp, &
      & 4.76567196308843E-04_wp, 1.40865283219372E-04_wp, 2.80111828722630E-04_wp, &
      &-1.78080710395585E-03_wp,-1.37122114854505E-03_wp,-1.40540921397035E-03_wp, &
      &-9.20931256837223E-05_wp,-8.23217732713214E-04_wp,-5.61547036728571E-04_wp, &
      & 7.11134130767397E-03_wp,-7.20199271675900E-04_wp, 7.01847980989513E-03_wp, &
      & 7.48778153956748E-04_wp,-1.25016037903597E-03_wp, 1.08310954300958E-03_wp, &
      &-1.30573252929538E-02_wp, 1.22148337982342E-02_wp,-8.65332883356591E-03_wp, &
      &-5.13301914431355E-04_wp, 2.05800937159412E-04_wp, 5.61141428283749E-04_wp, &
      & 1.30790693570089E-04_wp, 2.05141020566213E-04_wp, 7.03407939819235E-04_wp, &
      &-5.10132384573942E-04_wp, 1.58086691395969E-03_wp, 1.55347467449728E-03_wp, &
      &-8.38604250224578E-05_wp,-8.01750620859773E-05_wp,-1.15990435048486E-04_wp, &
      & 5.02969184274869E-05_wp, 4.77910989952591E-06_wp, 6.82449872522285E-04_wp, &
      & 1.42540648066691E-03_wp, 1.25335587673295E-03_wp,-6.61457072852844E-04_wp, &
      &-7.14058444774105E-05_wp, 3.86876319922551E-05_wp,-1.39324883856965E-04_wp, &
      &-1.15338713787178E-04_wp, 8.54214155075652E-05_wp, 7.56850448846645E-04_wp, &
      &-6.13369844793697E-04_wp,-6.46992100973259E-05_wp,-9.26782496704140E-05_wp, &
      & 1.28797404654250E-03_wp, 6.08343459325308E-04_wp,-1.46955760523574E-03_wp, &
      &-4.36565995938913E-06_wp,-1.31550509938713E-03_wp, 9.93031268178490E-05_wp, &
      &-1.79865124144556E-03_wp, 3.37358828176490E-03_wp, 1.03691084693786E-03_wp, &
      & 3.33671398085666E-03_wp, 2.13490629757382E-03_wp, 9.41199215808539E-04_wp, &
      &-2.17887584028734E-03_wp,-1.98932965996673E-05_wp, 1.49853968405253E-02_wp, &
      & 1.29451185710472E-02_wp, 3.54712418689991E-04_wp, 4.33613368926151E-04_wp, &
      & 1.30790693570089E-04_wp, 4.82942940074573E-04_wp, 2.51166095202739E-04_wp, &
      &-7.79526274556678E-04_wp,-7.94735706145713E-04_wp, 2.61363715975006E-04_wp, &
      &-3.32328565294179E-03_wp, 1.04443195498850E-04_wp, 6.48943271257431E-07_wp, &
      & 3.48534693272860E-06_wp, 1.06421932313923E-04_wp,-1.50279837156915E-04_wp, &
      & 9.18878731574177E-04_wp, 2.34666955737731E-04_wp,-3.00285320014155E-04_wp, &
      &-1.87675863045336E-03_wp,-4.37472910374419E-05_wp, 1.47388733357481E-05_wp, &
      &-9.08775083684892E-05_wp, 3.16503941397320E-05_wp, 1.29899358588447E-04_wp, &
      &-4.90010925833026E-04_wp,-7.09402995502100E-05_wp, 2.63809170113323E-04_wp, &
      & 7.36601934773917E-04_wp,-1.29393581855696E-03_wp,-4.94440667285671E-04_wp, &
      &-3.20358396726471E-05_wp,-1.55186627714891E-03_wp,-1.78591283773075E-03_wp, &
      & 2.01374750556590E-03_wp, 1.74946442106978E-03_wp, 1.44852013817612E-03_wp, &
      &-7.60692602493730E-04_wp, 1.34879934078008E-03_wp, 2.25243758392478E-03_wp, &
      &-3.02774896627336E-04_wp, 1.36316094026119E-03_wp,-2.41087735648506E-02_wp, &
      &-4.43859381082050E-04_wp, 1.70363857376580E-02_wp, 1.78411143916457E-05_wp, &
      &-3.53468297382543E-04_wp, 2.05141020566213E-04_wp, 2.51166095202739E-04_wp, &
      & 1.00086854659435E-03_wp,-2.26253358021328E-03_wp, 1.21437328753939E-03_wp, &
      & 1.43726580756122E-03_wp,-2.97126147911374E-03_wp, 1.86692869848502E-04_wp, &
      &-1.00760285830176E-04_wp, 4.75139681030006E-05_wp, 1.30570990143511E-04_wp, &
      &-4.85383826114350E-05_wp, 2.19789845074570E-03_wp, 5.36001426961137E-03_wp, &
      &-1.07860118926926E-04_wp,-2.11500611018870E-03_wp,-3.10082495884187E-04_wp, &
      &-1.22377866574111E-04_wp,-1.21743616812446E-04_wp, 4.07538213548353E-05_wp, &
      & 4.44904702132598E-05_wp,-1.05404313350324E-03_wp,-7.79253610457477E-03_wp, &
      &-3.70243865248362E-03_wp, 5.33025120676840E-03_wp,-2.40535817565888E-02_wp, &
      & 1.84723481845835E-02_wp,-9.72011490896808E-03_wp,-1.19396614153255E-02_wp, &
      &-8.62500600176080E-03_wp, 3.07722003547285E-04_wp, 1.11491519164976E-04_wp, &
      & 5.29641913037380E-04_wp,-3.00213735625271E-04_wp, 1.35487620190833E-03_wp, &
      &-5.92614759492204E-04_wp,-3.55591598093717E-04_wp,-1.74330915702325E-03_wp, &
      & 5.83917668822536E-03_wp, 9.69086766965422E-04_wp,-2.30363044077329E-02_wp, &
      &-9.90869743129443E-04_wp, 6.90995730459461E-04_wp, 7.03407939819235E-04_wp, &
      &-7.79526274556678E-04_wp,-2.26253358021328E-03_wp, 9.80374962527289E-01_wp, &
      & 6.95673580178762E-02_wp, 3.39764122184999E-02_wp,-4.74420205416907E-02_wp, &
      & 3.77754358456775E-03_wp,-2.82730146168479E-03_wp, 1.53500502661241E-03_wp, &
      & 1.89164343610649E-03_wp, 1.39975874696665E-03_wp,-3.03459327912099E-03_wp, &
      & 1.59267223668962E-02_wp, 1.87802573636107E-02_wp, 5.08789375835296E-03_wp, &
      &-1.68949070134437E-03_wp,-6.53393508005757E-04_wp,-1.24410308779058E-03_wp, &
      &-1.84667948741942E-03_wp,-3.94588356338352E-04_wp,-3.44887507450498E-02_wp, &
      &-1.13057960252247E-02_wp,-1.60655875718359E-02_wp, 2.09425905656062E-02_wp, &
      & 7.93344293670733E-02_wp,-6.44273893756731E-02_wp, 3.77933661329130E-04_wp, &
      & 7.87414432104843E-02_wp,-3.79155045233419E-02_wp,-5.84884648956699E-02_wp, &
      & 7.63322292457885E-02_wp, 8.80447817429636E-03_wp, 4.94121495627030E-02_wp, &
      & 8.94664814628892E-03_wp, 1.22050640255783E-02_wp,-4.93733275076534E-02_wp, &
      &-4.58677935890264E-03_wp, 4.54852027795121E-03_wp,-1.16033600518990E-02_wp, &
      & 1.65447430322963E-02_wp,-4.65550886867003E-04_wp,-2.30536700688280E-03_wp, &
      &-5.10132384573942E-04_wp,-7.94735706145713E-04_wp, 1.21437328753939E-03_wp, &
      & 6.95673580178762E-02_wp, 8.84180397452405E-01_wp,-5.73003494969123E-03_wp, &
      & 1.37764773885447E-02_wp, 1.70763073164178E-02_wp,-1.40683496133293E-02_wp, &
      & 1.28751303148498E-02_wp, 2.60641678220040E-04_wp, 2.31072589795466E-02_wp, &
      & 9.93166053400878E-03_wp,-2.98030521165622E-03_wp,-3.93048121184033E-02_wp, &
      &-1.30198217131533E-02_wp, 1.41350754924126E-03_wp,-5.48193100399673E-04_wp, &
      & 2.00987188975386E-03_wp, 3.42911090001512E-03_wp, 1.23010649004639E-04_wp, &
      &-1.34351050524750E-02_wp,-1.77696559058022E-02_wp, 1.10901118217271E-02_wp, &
      & 1.10099041175004E-02_wp, 7.02584780993291E-02_wp, 3.79307091876071E-02_wp, &
      & 8.63683672373309E-02_wp,-2.99966874007633E-02_wp, 3.52606751298176E-02_wp, &
      &-1.64182145464190E-02_wp, 4.93084007884058E-03_wp, 3.91330610474767E-02_wp, &
      & 6.09420749900072E-03_wp,-6.96974239270687E-02_wp,-5.67599294232971E-04_wp, &
      &-4.86973240502231E-02_wp,-4.64473531872142E-03_wp,-1.29884368229359E-02_wp, &
      & 2.04795523049283E-02_wp, 2.35203484160585E-03_wp,-1.08192538127679E-03_wp, &
      &-2.21987363481164E-04_wp, 1.58086691395969E-03_wp, 2.61363715975006E-04_wp, &
      & 1.43726580756122E-03_wp, 3.39764122184999E-02_wp,-5.73003494969123E-03_wp, &
      & 9.01059730065419E-01_wp, 8.62885023910562E-03_wp,-3.56783361149083E-04_wp, &
      &-2.30967611245580E-02_wp,-1.61244864531824E-02_wp, 1.68999482428444E-02_wp, &
      & 4.38646820706053E-04_wp, 1.88978542729194E-02_wp,-2.14362575560263E-02_wp, &
      &-1.69515545908515E-02_wp,-3.46246032784671E-02_wp, 1.95194206190891E-03_wp, &
      & 3.15786528053090E-04_wp, 1.90278129874468E-04_wp, 3.95027038964027E-03_wp, &
      & 3.86842771313291E-03_wp, 2.17894276273119E-02_wp, 2.23419903016915E-02_wp, &
      & 1.06050858375466E-02_wp, 3.15743661071952E-03_wp,-1.26549176620155E-02_wp, &
      & 6.86738430361502E-02_wp,-5.52544574468794E-03_wp,-2.09673932343892E-03_wp, &
      &-9.70310332800300E-02_wp,-5.32623004110966E-02_wp,-3.14965757562480E-02_wp, &
      & 9.30647517605180E-03_wp,-4.30860195133737E-02_wp, 1.83044694261429E-03_wp, &
      &-5.95871699567429E-02_wp, 4.52876588207916E-02_wp,-2.29331522657019E-02_wp, &
      &-1.20861841270074E-02_wp,-1.36569808889512E-02_wp,-5.05728174639121E-02_wp, &
      &-5.13523542037677E-03_wp,-6.65356782173023E-04_wp, 1.55347467449728E-03_wp, &
      &-3.32328565294179E-03_wp,-2.97126147911374E-03_wp,-4.74420205416907E-02_wp, &
      & 1.37764773885447E-02_wp, 8.62885023910562E-03_wp, 9.00665204134815E-01_wp, &
      &-2.32504544118005E-02_wp,-5.61086785463042E-04_wp,-1.03203148052035E-02_wp, &
      &-1.38491485560455E-02_wp, 1.79764937901381E-02_wp, 1.33708468639779E-02_wp, &
      & 1.11144973242871E-02_wp,-9.19867803125049E-03_wp,-2.07176142335574E-03_wp, &
      &-1.79485417279020E-03_wp,-9.30705196521984E-04_wp,-4.94262224653313E-04_wp, &
      & 2.70908618935943E-03_wp, 1.56234463093093E-03_wp,-1.32317745313018E-03_wp, &
      &-8.58013769944745E-04_wp,-6.86781164787524E-04_wp, 3.35355565772707E-04_wp, &
      & 1.71977080519939E-03_wp,-2.66785048528304E-03_wp, 2.37974672043536E-04_wp, &
      & 1.14415800389318E-03_wp, 1.09295237553953E-03_wp, 5.34718902385395E-04_wp, &
      & 2.32955790448110E-03_wp, 4.90411824161737E-05_wp, 2.32658753451363E-03_wp, &
      & 2.45059911759910E-04_wp, 1.86239809890842E-03_wp,-2.27875919841765E-03_wp, &
      & 1.08820023071544E-03_wp,-6.93166920400433E-04_wp,-7.21222341323190E-04_wp, &
      & 5.43932743796690E-03_wp, 2.00183991552993E-04_wp,-5.73474810553668E-05_wp, &
      &-8.38604250224578E-05_wp, 1.04443195498850E-04_wp, 1.86692869848502E-04_wp, &
      & 3.77754358456775E-03_wp, 1.70763073164178E-02_wp,-3.56783361149083E-04_wp, &
      &-2.32504544118005E-02_wp, 9.73649928304254E-04_wp,-2.66784572492629E-04_wp, &
      & 5.37753152378892E-04_wp, 3.79172201204915E-04_wp,-2.14001374746345E-05_wp, &
      &-1.08428106403024E-03_wp,-1.81560002174338E-03_wp,-1.02968022752516E-03_wp, &
      & 1.85080347591414E-03_wp, 1.23988028846638E-04_wp, 2.72413933173410E-05_wp, &
      & 9.71753060639905E-05_wp,-1.46388242648612E-05_wp,-7.69753532631084E-05_wp, &
      & 9.23295621268477E-04_wp, 6.57139878410130E-04_wp,-1.07485745623368E-05_wp, &
      &-7.16937802389864E-04_wp,-2.33627976151097E-03_wp, 3.39078065070956E-04_wp, &
      &-1.95351899121100E-03_wp,-1.63860594224303E-04_wp,-3.65714841904203E-04_wp, &
      & 1.82873529738723E-03_wp,-9.71790300967990E-04_wp,-2.82389512930550E-03_wp, &
      &-7.53508690225449E-04_wp, 1.87682772385266E-04_wp,-4.55777507061880E-04_wp, &
      & 2.50121224806569E-03_wp, 8.19971309071037E-04_wp, 2.31083568203780E-03_wp, &
      & 7.26268373831813E-05_wp,-9.04895911711816E-05_wp, 8.28467511050551E-05_wp, &
      & 4.56304467441363E-05_wp,-8.01750620859773E-05_wp, 6.48943271257431E-07_wp, &
      &-1.00760285830176E-04_wp,-2.82730146168479E-03_wp,-1.40683496133293E-02_wp, &
      &-2.30967611245580E-02_wp,-5.61086785463042E-04_wp,-2.66784572492629E-04_wp, &
      & 8.52282998751351E-04_wp, 2.20685799884168E-04_wp,-4.57974388265354E-04_wp, &
      &-3.99086808070311E-04_wp,-1.96684349758697E-03_wp, 3.16068626905044E-03_wp, &
      & 3.82500580393816E-03_wp, 2.13899661384175E-03_wp,-1.32999331985299E-04_wp, &
      &-2.33025828331968E-05_wp,-6.08937262019009E-05_wp,-2.32989245631772E-04_wp, &
      &-1.72467005604405E-04_wp,-5.90066147083195E-04_wp,-1.13106310747154E-04_wp, &
      &-5.80144587727174E-04_wp, 5.04915814179606E-05_wp, 1.31736482134134E-04_wp, &
      &-2.13068231597136E-03_wp,-1.29213626521191E-03_wp, 1.53063969395156E-03_wp, &
      &-3.71310081988481E-04_wp, 5.62884466208868E-04_wp, 2.00789686792316E-03_wp, &
      &-2.61623063595918E-03_wp, 1.24050701272183E-03_wp,-6.56177016050264E-04_wp, &
      & 4.50807312848272E-04_wp,-5.97769698513096E-05_wp, 8.60563453836057E-04_wp, &
      & 3.77274788892379E-04_wp,-1.52972790544857E-03_wp, 2.37278396459574E-03_wp, &
      & 1.24197069253257E-04_wp,-8.45754544145134E-05_wp,-1.15990435048486E-04_wp, &
      & 3.48534693272860E-06_wp, 4.75139681030006E-05_wp, 1.53500502661241E-03_wp, &
      & 1.28751303148498E-02_wp,-1.61244864531824E-02_wp,-1.03203148052035E-02_wp, &
      & 5.37753152378892E-04_wp, 2.20685799884168E-04_wp, 6.11130609792302E-04_wp, &
      &-1.47297049930426E-04_wp, 1.23700973779813E-04_wp,-1.25333148532490E-03_wp, &
      & 2.66385503686812E-04_wp, 1.26880999524027E-04_wp, 2.10248656475350E-03_wp, &
      & 1.31162314378301E-05_wp,-2.06468811675943E-05_wp, 7.42987296409726E-05_wp, &
      &-6.15055301804686E-05_wp,-1.53448604360204E-04_wp,-6.04014419841712E-04_wp, &
      &-7.25066154483291E-04_wp,-3.80504491466814E-05_wp, 2.25742379909139E-04_wp, &
      & 1.23882021027591E-03_wp,-3.19054602405769E-04_wp, 1.69539224395771E-03_wp, &
      &-6.84469334258254E-04_wp, 1.86764564795742E-03_wp, 1.47717103400572E-05_wp, &
      &-3.47566238648489E-04_wp, 3.78836597522281E-03_wp, 8.55022975789090E-04_wp, &
      & 1.97824728090476E-03_wp, 1.66783078270450E-03_wp,-2.37571013869545E-03_wp, &
      & 3.74866402148292E-04_wp,-9.37624772359046E-04_wp, 5.08502570020347E-04_wp, &
      & 2.62727701317609E-03_wp, 6.25152483639254E-05_wp, 5.05322897854274E-05_wp, &
      & 5.02969184274869E-05_wp, 1.06421932313923E-04_wp, 1.30570990143511E-04_wp, &
      & 1.89164343610649E-03_wp, 2.60641678220040E-04_wp, 1.68999482428444E-02_wp, &
      &-1.38491485560455E-02_wp, 3.79172201204915E-04_wp,-4.57974388265354E-04_wp, &
      &-1.47297049930426E-04_wp, 5.67497379548977E-04_wp,-2.66471277992456E-04_wp, &
      &-2.05537354959541E-05_wp,-2.24672138415218E-03_wp,-9.83367039048201E-04_wp, &
      &-9.36370582574761E-04_wp, 1.09878048189182E-04_wp, 6.51750410270745E-05_wp, &
      &-7.84162915168507E-06_wp, 3.18694253551562E-05_wp, 1.13232309414826E-04_wp, &
      &-5.29940796996229E-04_wp, 1.72200355417520E-04_wp,-1.88723311255307E-04_wp, &
      & 6.87321959615046E-04_wp, 1.31126864067533E-03_wp,-6.05877111680159E-04_wp, &
      &-2.70533226346614E-04_wp, 1.76370008603692E-03_wp,-2.30237597503765E-03_wp, &
      &-2.63392884744990E-03_wp, 1.90671492659640E-03_wp,-7.48604330870602E-04_wp, &
      & 1.67234209279600E-04_wp,-1.35163183416280E-03_wp,-1.27146849958339E-03_wp, &
      &-2.62237525358857E-04_wp,-2.27285691797216E-03_wp,-9.82783410965334E-04_wp, &
      &-2.09194354735716E-03_wp,-2.77814041002767E-03_wp,-1.66773246000592E-04_wp, &
      &-1.40434975631984E-04_wp, 4.77910989952591E-06_wp,-1.50279837156915E-04_wp, &
      &-4.85383826114350E-05_wp, 1.39975874696665E-03_wp, 2.31072589795466E-02_wp, &
      & 4.38646820706053E-04_wp, 1.79764937901381E-02_wp,-2.14001374746345E-05_wp, &
      &-3.99086808070311E-04_wp, 1.23700973779813E-04_wp,-2.66471277992456E-04_wp, &
      & 9.79650219582221E-04_wp, 1.37107013798353E-03_wp,-3.44485251776749E-04_wp, &
      &-4.36588225424418E-03_wp,-8.48621751587428E-04_wp, 8.58216130174137E-06_wp, &
      &-6.62860477387958E-05_wp, 1.05721183585390E-04_wp, 2.47637173022898E-04_wp, &
      & 3.18214998614067E-05_wp,-9.36839279449376E-04_wp, 2.47719305207997E-03_wp, &
      &-3.49727598971794E-03_wp,-8.83698321844673E-03_wp,-1.35887092213168E-02_wp, &
      &-5.22575499812361E-03_wp,-9.91706868410338E-03_wp, 2.14468953008398E-02_wp, &
      & 2.20223984658526E-02_wp,-1.29176440695309E-04_wp, 8.63738096668437E-04_wp, &
      &-1.70308189354695E-04_wp,-2.42765201916911E-04_wp,-1.19394373726397E-03_wp, &
      & 1.03573655632359E-04_wp,-2.86493133067578E-04_wp,-4.66588242543135E-04_wp, &
      &-1.79306724772635E-02_wp, 1.70652047126286E-03_wp, 1.34190851188152E-02_wp, &
      & 9.05681949044069E-04_wp,-3.66882904454643E-04_wp, 6.82449872522285E-04_wp, &
      & 9.18878731574177E-04_wp, 2.19789845074570E-03_wp,-3.03459327912099E-03_wp, &
      & 9.93166053400878E-03_wp, 1.88978542729194E-02_wp, 1.33708468639779E-02_wp, &
      &-1.08428106403024E-03_wp,-1.96684349758697E-03_wp,-1.25333148532490E-03_wp, &
      &-2.05537354959541E-05_wp, 1.37107013798353E-03_wp, 9.79735154223244E-01_wp, &
      &-2.24094441813411E-02_wp, 3.40734168483049E-02_wp, 8.19089170158313E-02_wp, &
      & 2.08377777367126E-03_wp, 8.38445627491986E-04_wp, 1.55818874225189E-03_wp, &
      &-3.26712147082133E-03_wp,-3.45181101321062E-03_wp, 9.73171638458462E-03_wp, &
      & 1.50974537942946E-02_wp, 4.62053531076067E-03_wp, 1.35558351291241E-02_wp, &
      & 7.65773113004499E-02_wp, 3.76217969547401E-02_wp,-5.51719337384838E-03_wp, &
      & 3.28265258373041E-02_wp, 7.33725870090370E-02_wp, 3.91140770304674E-02_wp, &
      & 2.78573450072921E-02_wp, 1.13610126762784E-02_wp,-2.43254052385456E-02_wp, &
      & 2.31845543892078E-02_wp, 4.75824479949787E-02_wp, 5.67692155285208E-02_wp, &
      &-2.18088403427284E-02_wp,-4.45043895094737E-02_wp,-1.64247670292941E-02_wp, &
      & 9.68094770971209E-03_wp,-1.85522927119767E-03_wp,-3.35446033168041E-03_wp, &
      & 1.42540648066691E-03_wp, 2.34666955737731E-04_wp, 5.36001426961137E-03_wp, &
      & 1.59267223668962E-02_wp,-2.98030521165622E-03_wp,-2.14362575560263E-02_wp, &
      & 1.11144973242871E-02_wp,-1.81560002174338E-03_wp, 3.16068626905044E-03_wp, &
      & 2.66385503686812E-04_wp,-2.24672138415218E-03_wp,-3.44485251776749E-04_wp, &
      &-2.24094441813411E-02_wp, 9.06985941795760E-01_wp, 4.93082866596486E-03_wp, &
      & 6.81990810666048E-03_wp,-2.76733628649532E-02_wp,-1.32975610219188E-02_wp, &
      &-5.75746263646097E-03_wp,-8.97411598060613E-04_wp,-8.63422164122171E-03_wp, &
      &-1.32361089963587E-02_wp, 4.19681624623315E-03_wp, 1.19028969178129E-02_wp, &
      &-2.03523592200911E-02_wp, 3.30174684587050E-02_wp,-1.71192760195093E-02_wp, &
      & 8.65390101840318E-02_wp, 4.80550657744445E-02_wp,-7.21093493716756E-02_wp, &
      & 5.34660261424884E-02_wp, 8.43233407925204E-03_wp,-6.50954576516091E-02_wp, &
      & 5.40323761942012E-03_wp, 2.60925461908099E-02_wp,-1.10491431748055E-02_wp, &
      &-2.51528185615128E-02_wp,-4.84923506617974E-03_wp,-4.51857519227536E-03_wp, &
      & 2.26152065790688E-02_wp,-1.32881963094551E-02_wp,-1.54510199955994E-03_wp, &
      &-3.70548351221492E-04_wp, 1.25335587673295E-03_wp,-3.00285320014155E-04_wp, &
      &-1.07860118926926E-04_wp, 1.87802573636107E-02_wp,-3.93048121184033E-02_wp, &
      &-1.69515545908515E-02_wp,-9.19867803125049E-03_wp,-1.02968022752516E-03_wp, &
      & 3.82500580393816E-03_wp, 1.26880999524027E-04_wp,-9.83367039048201E-04_wp, &
      &-4.36588225424418E-03_wp, 3.40734168483049E-02_wp, 4.93082866596486E-03_wp, &
      & 9.03728584036943E-01_wp,-1.07443955567044E-02_wp,-5.83128634967848E-04_wp, &
      & 7.58961269677126E-03_wp,-1.49660732658146E-02_wp,-2.62733819940710E-02_wp, &
      &-1.06417180484544E-03_wp,-4.00527607361329E-02_wp, 1.19489455972230E-02_wp, &
      &-1.89336042467516E-02_wp,-2.27896269355215E-02_wp, 6.43422160137272E-02_wp, &
      & 4.37862349890877E-02_wp, 3.02003084647696E-03_wp,-1.04816257994469E-01_wp, &
      &-4.68711668322494E-02_wp, 5.84896470057643E-02_wp, 4.10236255981251E-02_wp, &
      & 2.26688923785656E-02_wp, 6.11536812838240E-02_wp, 2.77813689934099E-02_wp, &
      &-6.98072919459784E-02_wp,-3.76025089855853E-02_wp, 2.45274151312103E-03_wp, &
      & 3.23050890198580E-02_wp,-9.39650646415531E-03_wp,-3.36509581492724E-03_wp, &
      &-8.07848283131484E-04_wp,-2.78556704937534E-04_wp,-6.61457072852844E-04_wp, &
      &-1.87675863045336E-03_wp,-2.11500611018870E-03_wp, 5.08789375835296E-03_wp, &
      &-1.30198217131533E-02_wp,-3.46246032784671E-02_wp,-2.07176142335574E-03_wp, &
      & 1.85080347591414E-03_wp, 2.13899661384175E-03_wp, 2.10248656475350E-03_wp, &
      &-9.36370582574761E-04_wp,-8.48621751587428E-04_wp, 8.19089170158313E-02_wp, &
      & 6.81990810666048E-03_wp,-1.07443955567044E-02_wp, 8.78160106104680E-01_wp, &
      & 8.67164949805118E-03_wp, 7.25413397711342E-04_wp, 1.48329982255582E-02_wp, &
      &-1.41349460904117E-02_wp,-2.62925477940729E-02_wp,-7.26493927890594E-04_wp, &
      &-4.06772228425090E-04_wp,-4.11421365806423E-04_wp,-6.80990118806957E-04_wp, &
      &-1.00280886934538E-03_wp,-7.45923917349099E-04_wp, 3.82254327745874E-04_wp, &
      &-1.57983395972902E-03_wp,-2.11281803772875E-03_wp,-1.07874841815666E-03_wp, &
      &-7.73918178641165E-04_wp, 6.83511061285916E-04_wp, 1.67995322470535E-03_wp, &
      & 2.20058704998442E-04_wp,-2.04606432020437E-03_wp,-2.46663081102699E-03_wp, &
      & 2.10826592111825E-03_wp, 5.49801719243796E-03_wp, 5.30549401718718E-04_wp, &
      &-2.43229181329884E-03_wp, 5.56541903522150E-05_wp, 1.62696960532809E-04_wp, &
      &-7.14058444774105E-05_wp,-4.37472910374419E-05_wp,-3.10082495884187E-04_wp, &
      &-1.68949070134437E-03_wp, 1.41350754924126E-03_wp, 1.95194206190891E-03_wp, &
      &-1.79485417279020E-03_wp, 1.23988028846638E-04_wp,-1.32999331985299E-04_wp, &
      & 1.31162314378301E-05_wp, 1.09878048189182E-04_wp, 8.58216130174137E-06_wp, &
      & 2.08377777367126E-03_wp,-2.76733628649532E-02_wp,-5.83128634967848E-04_wp, &
      & 8.67164949805118E-03_wp, 9.70600754648758E-04_wp, 4.29953409489857E-04_wp, &
      & 3.37914573840815E-04_wp,-1.08472703426642E-04_wp, 7.46795967613791E-06_wp, &
      &-2.84864369040110E-04_wp,-1.49187995855161E-04_wp,-8.15783274082982E-05_wp, &
      &-4.08274548452866E-04_wp,-7.03926923327412E-04_wp,-5.93743903736391E-04_wp, &
      & 9.48495944206912E-04_wp,-8.84185331974615E-05_wp,-1.44674369270472E-03_wp, &
      &-7.53053375859206E-04_wp,-1.64698024921482E-03_wp, 3.49269145418072E-03_wp, &
      & 5.55389209372992E-04_wp, 4.37711797391931E-03_wp,-4.78902053888293E-06_wp, &
      &-2.04537262522129E-03_wp, 5.15305178827397E-04_wp, 2.23738762999721E-03_wp, &
      & 1.21492153053422E-04_wp,-1.01164103858773E-03_wp,-1.75762882400090E-05_wp, &
      & 1.32878219836544E-04_wp, 3.86876319922551E-05_wp, 1.47388733357481E-05_wp, &
      &-1.22377866574111E-04_wp,-6.53393508005757E-04_wp,-5.48193100399673E-04_wp, &
      & 3.15786528053090E-04_wp,-9.30705196521984E-04_wp, 2.72413933173410E-05_wp, &
      &-2.33025828331968E-05_wp,-2.06468811675943E-05_wp, 6.51750410270745E-05_wp, &
      &-6.62860477387958E-05_wp, 8.38445627491986E-04_wp,-1.32975610219188E-02_wp, &
      & 7.58961269677126E-03_wp, 7.25413397711342E-04_wp, 4.29953409489857E-04_wp, &
      & 3.06994419944357E-04_wp,-6.88834984059164E-05_wp,-2.72138730458034E-04_wp, &
      & 1.36166386900313E-04_wp,-5.84337860555941E-04_wp,-1.57324892220013E-05_wp, &
      &-5.09349489441788E-04_wp,-1.56609917648878E-04_wp, 3.69301350068986E-04_wp, &
      & 5.99040723188793E-04_wp,-1.24471774247931E-03_wp,-2.38428284607014E-03_wp, &
      & 1.06345967457132E-04_wp, 3.71134001508957E-04_wp, 1.82838363631060E-03_wp, &
      &-3.31410275388691E-03_wp, 1.20747771046646E-03_wp,-5.12579723062006E-03_wp, &
      &-2.27508889527851E-03_wp, 2.56519694622059E-04_wp, 8.24138845319420E-04_wp, &
      & 2.14802109998526E-03_wp,-1.14204140349652E-03_wp,-7.48267748507518E-04_wp, &
      & 6.87763115435970E-05_wp,-6.24460321559444E-05_wp,-1.39324883856965E-04_wp, &
      &-9.08775083684892E-05_wp,-1.21743616812446E-04_wp,-1.24410308779058E-03_wp, &
      & 2.00987188975386E-03_wp, 1.90278129874468E-04_wp,-4.94262224653313E-04_wp, &
      & 9.71753060639905E-05_wp,-6.08937262019009E-05_wp, 7.42987296409726E-05_wp, &
      &-7.84162915168507E-06_wp, 1.05721183585390E-04_wp, 1.55818874225189E-03_wp, &
      &-5.75746263646097E-03_wp,-1.49660732658146E-02_wp, 1.48329982255582E-02_wp, &
      & 3.37914573840815E-04_wp,-6.88834984059164E-05_wp, 5.89819115781267E-04_wp, &
      & 2.52240655439004E-04_wp,-4.23234699536520E-04_wp, 1.06457919323247E-03_wp, &
      &-4.42293634056734E-04_wp, 4.03617862685355E-05_wp, 1.00472108720785E-03_wp, &
      &-1.38996067520150E-03_wp, 9.36631352457798E-05_wp,-2.38179137703504E-03_wp, &
      & 4.96040206857101E-04_wp, 2.30301477281955E-03_wp,-2.11151814934020E-03_wp, &
      & 3.86732721439976E-04_wp,-2.97731431290906E-03_wp,-9.68228710592495E-04_wp, &
      &-6.39253710685816E-03_wp, 4.87698483462043E-04_wp, 2.33412819707485E-03_wp, &
      & 7.13070047842620E-04_wp, 2.15131149462760E-04_wp, 3.16308765273205E-04_wp, &
      & 2.35043433033507E-03_wp, 1.51322622194485E-04_wp,-4.09463185306541E-05_wp, &
      &-1.15338713787178E-04_wp, 3.16503941397320E-05_wp, 4.07538213548353E-05_wp, &
      &-1.84667948741942E-03_wp, 3.42911090001512E-03_wp, 3.95027038964027E-03_wp, &
      & 2.70908618935943E-03_wp,-1.46388242648612E-05_wp,-2.32989245631772E-04_wp, &
      &-6.15055301804686E-05_wp, 3.18694253551562E-05_wp, 2.47637173022898E-04_wp, &
      &-3.26712147082133E-03_wp,-8.97411598060613E-04_wp,-2.62733819940710E-02_wp, &
      &-1.41349460904117E-02_wp,-1.08472703426642E-04_wp,-2.72138730458034E-04_wp, &
      & 2.52240655439004E-04_wp, 1.08240217357025E-03_wp, 4.36693703044900E-04_wp, &
      & 1.20863238120399E-03_wp,-5.24582933230386E-04_wp, 4.69679477374292E-04_wp, &
      & 6.42113030611923E-04_wp,-2.37514633952885E-03_wp,-1.11597837195792E-03_wp, &
      & 4.90468431720295E-05_wp, 2.58029772016527E-03_wp, 5.27801603506582E-04_wp, &
      &-2.94550580024158E-03_wp,-3.09148435784959E-03_wp, 4.21575971843549E-03_wp, &
      &-1.52225439381485E-03_wp, 3.81730219384280E-03_wp, 2.67025836021827E-03_wp, &
      &-2.89987015622134E-04_wp, 1.14992365340433E-03_wp,-2.63690796504470E-04_wp, &
      & 1.51716666848909E-03_wp, 1.60863886138561E-03_wp, 3.34775023659391E-05_wp, &
      & 1.51524804433381E-04_wp, 8.54214155075652E-05_wp, 1.29899358588447E-04_wp, &
      & 4.44904702132598E-05_wp,-3.94588356338352E-04_wp, 1.23010649004639E-04_wp, &
      & 3.86842771313291E-03_wp, 1.56234463093093E-03_wp,-7.69753532631084E-05_wp, &
      &-1.72467005604405E-04_wp,-1.53448604360204E-04_wp, 1.13232309414826E-04_wp, &
      & 3.18214998614067E-05_wp,-3.45181101321062E-03_wp,-8.63422164122171E-03_wp, &
      &-1.06417180484544E-03_wp,-2.62925477940729E-02_wp, 7.46795967613791E-06_wp, &
      & 1.36166386900313E-04_wp,-4.23234699536520E-04_wp, 4.36693703044900E-04_wp, &
      & 9.39787513840761E-04_wp, 7.12005478656161E-03_wp,-4.51310814188260E-04_wp, &
      & 7.03701001093687E-04_wp,-2.94299487920766E-04_wp,-2.72611976822528E-03_wp, &
      & 6.83963626666539E-03_wp,-5.31303754025144E-03_wp, 4.86750852539054E-03_wp, &
      & 9.51752790468847E-04_wp,-2.43673258893618E-03_wp,-9.06588364218912E-03_wp, &
      &-3.79091162487751E-04_wp,-2.36936520816196E-03_wp,-3.07632771315538E-04_wp, &
      & 3.36923184300490E-03_wp, 8.99317487950987E-03_wp,-1.09160503081918E-03_wp, &
      & 3.20293314794187E-02_wp, 2.32797874206224E-02_wp, 2.11137376206241E-02_wp, &
      & 1.14014043976149E-03_wp, 1.14151243425104E-03_wp,-3.01493542097973E-04_wp, &
      & 7.61823161502331E-04_wp,-4.78022561902225E-04_wp,-1.33276214891788E-03_wp, &
      &-3.51386816605225E-02_wp,-1.39271236840666E-02_wp, 2.23819392409864E-02_wp, &
      &-1.29940447168148E-03_wp, 9.04158731595418E-04_wp,-5.82863124169693E-04_wp, &
      &-5.86304925362903E-04_wp,-5.20787144981459E-04_wp,-1.20082797898505E-03_wp, &
      & 1.01678832885046E-02_wp,-1.38192679518940E-02_wp,-4.08823657751133E-02_wp, &
      &-7.15273250933024E-04_wp,-2.62747725165519E-04_wp,-5.93120928958409E-04_wp, &
      & 1.02934066632207E-03_wp, 1.20857770847103E-03_wp,-4.51310814188260E-04_wp, &
      & 2.28710762544773E-03_wp, 5.41555921947793E-04_wp,-2.60454699305636E-04_wp, &
      &-2.64761084039973E-03_wp, 5.74937304198407E-04_wp,-4.73555100697373E-05_wp, &
      &-1.70452681738136E-03_wp,-2.03891789976642E-04_wp, 1.06857216751498E-03_wp, &
      & 1.23576961126418E-03_wp,-1.23906888668214E-03_wp,-2.62937520452631E-03_wp, &
      & 6.80225612020106E-04_wp,-2.77935808450687E-03_wp, 6.09598880255074E-04_wp, &
      & 7.15239177618886E-03_wp,-1.15189576864124E-02_wp,-1.67088328726513E-02_wp, &
      &-1.96392476573831E-02_wp,-8.08703725341765E-04_wp,-6.27557723472192E-04_wp, &
      & 1.20819065008304E-04_wp,-6.29338096627950E-04_wp,-7.69573321153983E-05_wp, &
      &-8.03102785109213E-03_wp,-1.18519098620367E-02_wp,-1.74731534485127E-02_wp, &
      & 2.21934571114739E-02_wp,-8.42367911860373E-04_wp, 6.54572936555637E-04_wp, &
      &-9.79001755995669E-05_wp,-7.38113927657857E-04_wp, 1.66609509564965E-04_wp, &
      & 2.55368914978415E-03_wp, 1.44150017668853E-02_wp, 4.82113542118391E-03_wp, &
      & 1.15441408104976E-02_wp,-3.91702612956306E-04_wp,-1.85560188597899E-04_wp, &
      & 2.73984706658681E-05_wp,-3.83899731674412E-04_wp,-5.51201269388917E-04_wp, &
      & 7.03701001093687E-04_wp, 5.41555921947793E-04_wp, 1.62759410878750E-03_wp, &
      & 3.60831842873332E-04_wp,-1.63343497710738E-03_wp,-1.08136334407647E-03_wp, &
      & 3.26049513797072E-04_wp,-6.81691203926653E-04_wp, 5.74633576946177E-04_wp, &
      &-1.94158806191590E-04_wp,-5.97682632358859E-04_wp, 1.49160580801637E-03_wp, &
      &-2.93069153287612E-03_wp, 1.04587611852477E-03_wp, 1.49439576750988E-04_wp, &
      & 6.08564494858310E-04_wp, 4.56133941233273E-03_wp,-1.88851883302615E-02_wp, &
      & 8.77992495294718E-03_wp,-1.26538773814726E-02_wp,-7.09842519054558E-04_wp, &
      & 2.78132229357562E-06_wp, 6.06385637342826E-04_wp, 3.06228686912739E-06_wp, &
      & 2.86568553158027E-04_wp,-3.81506946938390E-03_wp,-1.56903929838757E-02_wp, &
      & 9.92777633341959E-03_wp, 1.08044883389166E-02_wp,-6.58362672011616E-04_wp, &
      &-2.66475674848386E-05_wp,-5.84082021986155E-04_wp, 1.60899708568502E-05_wp, &
      &-1.95975352223134E-04_wp,-3.62978936655017E-03_wp, 5.13925313958751E-03_wp, &
      & 1.02205843520703E-02_wp,-1.81593556200167E-02_wp,-3.90295021261056E-04_wp, &
      & 9.37999895527612E-06_wp,-5.75461930289809E-04_wp,-4.74477891207914E-05_wp, &
      & 5.45960716411972E-04_wp,-2.94299487920766E-04_wp,-2.60454699305636E-04_wp, &
      & 3.60831842873332E-04_wp, 2.43634114124241E-03_wp,-1.06882662768228E-03_wp, &
      &-1.65294351161779E-03_wp,-7.33701778288681E-05_wp, 1.79929228447519E-03_wp, &
      & 2.13659131646830E-03_wp,-1.80821431790967E-03_wp, 2.37534010499930E-03_wp, &
      & 6.75003273155266E-04_wp,-1.81229618483608E-03_wp,-1.71562995623376E-03_wp, &
      & 2.35017057420186E-03_wp,-2.96432814230523E-04_wp, 4.73945653321463E-03_wp, &
      &-1.99644814082213E-02_wp,-1.13633091808742E-02_wp, 4.51838948509923E-03_wp, &
      &-2.45727052670871E-04_wp,-6.32157378714537E-04_wp, 8.19701869394340E-05_wp, &
      &-1.19768297628478E-04_wp, 7.02087724447307E-04_wp, 5.48748418783290E-03_wp, &
      & 2.07344255347774E-02_wp, 1.12782220484343E-02_wp, 2.59864135990082E-03_wp, &
      & 3.30597234140415E-04_wp,-6.86166525721677E-04_wp, 6.16125688124352E-05_wp, &
      & 2.02646772793093E-04_wp, 6.65298252647229E-04_wp,-9.07158901302496E-03_wp, &
      & 1.33594273132114E-02_wp,-1.99274198938198E-02_wp,-2.34198218006705E-02_wp, &
      &-6.68613327789153E-04_wp,-4.37632534503615E-04_wp,-1.19585651284702E-04_wp, &
      & 1.01952218604768E-03_wp, 5.93888712235277E-04_wp,-2.72611976822528E-03_wp, &
      &-2.64761084039973E-03_wp,-1.63343497710738E-03_wp,-1.06882662768228E-03_wp, &
      & 4.26628832660554E-02_wp, 5.49747856041775E-03_wp, 2.52276146921499E-03_wp, &
      & 7.36515487172029E-03_wp,-4.06092774742350E-04_wp,-3.22068686629460E-04_wp, &
      & 3.89591208745783E-03_wp, 4.64033606839685E-03_wp, 1.33984180274338E-02_wp, &
      & 7.59055857447997E-04_wp, 8.14807186834400E-04_wp,-2.97423612012955E-04_wp, &
      & 2.23748991423594E-02_wp, 6.59864328579052E-02_wp, 8.69171101021901E-02_wp, &
      & 1.24058369540158E-03_wp, 8.43407994918159E-04_wp, 2.57087657914212E-03_wp, &
      & 4.72280474640535E-04_wp, 1.34942396110424E-03_wp,-1.20043130191039E-03_wp, &
      &-2.39746394216992E-02_wp, 7.91327525317699E-02_wp, 6.88040171310118E-02_wp, &
      &-1.27823598977582E-02_wp, 1.67224810045337E-03_wp,-2.26628072621023E-03_wp, &
      & 1.02410092154246E-04_wp, 1.25422685598679E-03_wp, 1.23413766856317E-03_wp, &
      &-1.35647385583807E-02_wp, 7.56592969150405E-02_wp, 3.15195489468795E-02_wp, &
      & 6.45164569511301E-02_wp,-9.39640522978982E-04_wp,-6.34568358829660E-04_wp, &
      & 3.16393698227943E-04_wp,-1.37698038383402E-03_wp,-2.23650440108797E-03_wp, &
      & 6.83963626666539E-03_wp, 5.74937304198407E-04_wp,-1.08136334407647E-03_wp, &
      &-1.65294351161779E-03_wp, 5.49747856041775E-03_wp, 3.04747470562918E-02_wp, &
      & 5.04578236073762E-03_wp,-4.64580168198365E-03_wp,-7.19056998586050E-03_wp, &
      & 3.26704088425329E-03_wp,-1.16593166664047E-02_wp,-1.76046211725149E-03_wp, &
      & 1.32231457493828E-03_wp,-5.03958332987170E-03_wp,-8.83041595713061E-04_wp, &
      & 9.84459994392474E-03_wp, 2.07361056877234E-02_wp, 6.71002330332418E-02_wp, &
      &-1.96682616947918E-02_wp, 7.82218702145295E-02_wp, 2.84314378706769E-03_wp, &
      & 3.82184571884326E-04_wp,-1.93416448459594E-03_wp, 4.83643649261318E-04_wp, &
      &-5.30311417466847E-04_wp, 1.83841041806474E-02_wp,-6.46154002416890E-02_wp, &
      & 3.76013880919811E-02_wp, 6.76244370380171E-02_wp,-2.56243044267504E-03_wp, &
      & 3.58655471763609E-04_wp,-2.02294067435815E-03_wp,-3.48741034904623E-04_wp, &
      &-5.80912056614872E-04_wp,-5.21088318946672E-03_wp, 3.66481656934027E-02_wp, &
      &-1.62019419351131E-02_wp, 4.22691260020166E-02_wp,-7.23873876527143E-04_wp, &
      &-6.94871573164942E-04_wp, 6.92791471894560E-04_wp, 2.14712344443191E-04_wp, &
      &-1.19295902152296E-03_wp,-5.31303754025144E-03_wp,-4.73555100697373E-05_wp, &
      & 3.26049513797072E-04_wp,-7.33701778288681E-05_wp, 2.52276146921499E-03_wp, &
      & 5.04578236073762E-03_wp, 2.84678467269820E-02_wp, 3.55720928362012E-03_wp, &
      &-7.76387042120417E-04_wp, 3.40221397139322E-03_wp, 1.80491436909120E-03_wp, &
      &-1.74606285771953E-03_wp, 2.05719732236316E-04_wp,-9.55018724841159E-04_wp, &
      &-6.73908004567400E-04_wp,-1.24968367218685E-02_wp,-7.48356444965498E-03_wp, &
      & 1.50312087642997E-02_wp,-9.52279424705729E-02_wp, 1.07162329819618E-02_wp, &
      & 2.74271011527141E-04_wp,-2.28230447764398E-03_wp,-1.54891424625151E-03_wp, &
      &-1.52565429299107E-03_wp,-8.33157167079826E-05_wp,-9.69794623656807E-03_wp, &
      & 5.97324438688014E-04_wp, 8.53965281948876E-02_wp,-6.23157155351273E-03_wp, &
      & 2.39631207747973E-04_wp,-1.86189149033346E-03_wp,-1.22978146596578E-03_wp, &
      & 1.62253168231263E-03_wp,-2.57581016819968E-04_wp,-9.88315126715304E-03_wp, &
      &-6.24835621250338E-03_wp, 8.60672941363428E-02_wp, 3.14814051825769E-03_wp, &
      & 3.77981999348335E-04_wp, 9.34039727567837E-04_wp,-1.19873933534838E-03_wp, &
      &-2.29602091271526E-03_wp, 4.64422684027194E-05_wp, 4.86750852539054E-03_wp, &
      &-1.70452681738136E-03_wp,-6.81691203926653E-04_wp, 1.79929228447519E-03_wp, &
      & 7.36515487172029E-03_wp,-4.64580168198365E-03_wp, 3.55720928362012E-03_wp, &
      & 3.34375782587607E-02_wp, 4.74783493233757E-03_wp,-9.19753668187009E-03_wp, &
      &-4.80342778340768E-03_wp,-5.05146906110698E-03_wp, 6.81728944109617E-04_wp, &
      & 2.02662411163210E-03_wp, 9.45853199801574E-03_wp, 3.79403591761504E-03_wp, &
      & 1.40905795163860E-02_wp, 7.86252686133398E-02_wp,-1.34685700986416E-02_wp, &
      & 5.53966133581774E-03_wp, 1.51709946910195E-03_wp, 4.93883363269233E-04_wp, &
      &-1.33728430469689E-03_wp, 9.51914809403550E-06_wp,-1.49346420012755E-03_wp, &
      &-1.19096389468461E-02_wp, 7.76490812173268E-02_wp,-3.00304876167653E-02_wp, &
      &-2.01300808356742E-03_wp, 1.10201283191560E-03_wp,-1.69175844074889E-04_wp, &
      & 1.46303970828171E-03_wp,-6.46390412013604E-04_wp, 1.68271214865759E-03_wp, &
      & 2.13613175487847E-02_wp, 3.25299748661382E-02_wp, 4.67415065569315E-02_wp, &
      &-1.03864191994931E-01_wp,-1.50087755722131E-03_wp, 3.62418390726754E-05_wp, &
      &-2.39989499316024E-03_wp, 3.71610089264428E-04_wp, 2.60085995965948E-03_wp, &
      & 9.51752790468847E-04_wp,-2.03891789976642E-04_wp, 5.74633576946177E-04_wp, &
      & 2.13659131646830E-03_wp,-4.06092774742350E-04_wp,-7.19056998586050E-03_wp, &
      &-7.76387042120417E-04_wp, 4.74783493233757E-03_wp, 4.18714281584843E-02_wp, &
      &-2.52340854281525E-03_wp, 7.67210126584856E-04_wp, 2.11987439578359E-03_wp, &
      &-4.97769753948644E-03_wp,-5.83490926431050E-03_wp, 5.49862980270049E-03_wp, &
      &-4.23716574381990E-04_wp,-8.82914903560124E-03_wp, 3.37935054072574E-02_wp, &
      &-3.37085582367700E-02_wp,-8.79319158772149E-02_wp,-1.20607118888942E-03_wp, &
      &-4.47456238961509E-04_wp,-1.80831545521759E-04_wp,-1.30457359386324E-03_wp, &
      &-1.72606096100193E-03_wp,-8.62164670054123E-03_wp,-3.71397599691318E-02_wp, &
      & 3.53408962248655E-02_wp,-9.54603867530334E-02_wp, 1.04148511475808E-03_wp, &
      &-3.58853452102961E-04_wp,-3.58903840536546E-04_wp, 1.79457630154229E-03_wp, &
      &-2.19066197265030E-03_wp, 2.19877649688146E-02_wp, 7.20879127329431E-02_wp, &
      &-7.10369548170602E-02_wp,-4.72439787083179E-02_wp,-2.04362405293725E-03_wp, &
      &-1.49241311598283E-03_wp, 1.88746416191365E-04_wp, 2.32194947529299E-03_wp, &
      & 4.24533792527953E-04_wp,-2.43673258893618E-03_wp, 1.06857216751498E-03_wp, &
      &-1.94158806191590E-04_wp,-1.80821431790967E-03_wp,-3.22068686629460E-04_wp, &
      & 3.26704088425329E-03_wp, 3.40221397139322E-03_wp,-9.19753668187009E-03_wp, &
      &-2.52340854281525E-03_wp, 1.18754589864597E-02_wp, 3.12284929361183E-03_wp, &
      &-1.79484074665602E-03_wp, 1.08426219074813E-03_wp, 2.95149250539913E-03_wp, &
      & 9.39535373805293E-04_wp,-8.20550769743827E-04_wp, 2.51208188470779E-05_wp, &
      &-3.43256354112925E-02_wp,-1.00847249630762E-02_wp, 3.53479811853906E-02_wp, &
      & 3.47939694975228E-04_wp,-7.43623597157547E-04_wp, 3.68449955830016E-06_wp, &
      & 3.27464191561612E-04_wp, 1.58744109760933E-03_wp, 5.59683214958685E-04_wp, &
      &-4.32983426569359E-02_wp,-1.56788026891989E-02_wp,-3.91999274242567E-02_wp, &
      & 4.04261927191880E-04_wp, 1.24003896893423E-03_wp, 2.73469257707459E-04_wp, &
      & 2.83112352296633E-04_wp,-2.02065592692274E-03_wp,-8.80011993355233E-04_wp, &
      & 3.13091822106836E-02_wp, 3.67758849915333E-02_wp, 5.10602082770450E-02_wp, &
      &-6.95431065144622E-04_wp,-1.53919581309734E-04_wp,-6.53037749246239E-05_wp, &
      &-2.04327967176310E-03_wp,-1.86505232958969E-03_wp,-9.06588364218912E-03_wp, &
      & 1.23576961126418E-03_wp,-5.97682632358859E-04_wp, 2.37534010499930E-03_wp, &
      & 3.89591208745783E-03_wp,-1.16593166664047E-02_wp, 1.80491436909120E-03_wp, &
      &-4.80342778340768E-03_wp, 7.67210126584856E-04_wp, 3.12284929361183E-03_wp, &
      & 1.53464414780711E-02_wp, 4.96984189818760E-04_wp, 1.15257191840106E-03_wp, &
      & 3.14378953074075E-04_wp,-1.87096506480690E-03_wp,-8.95788403453741E-03_wp, &
      &-9.71288715620663E-04_wp,-6.97162425658716E-02_wp,-2.42710112695174E-02_wp, &
      &-2.43258070526638E-02_wp,-1.88935362844706E-03_wp,-1.84650959078125E-03_wp, &
      & 8.39626505216267E-04_wp,-8.32108659053711E-04_wp, 1.57482674045950E-03_wp, &
      &-6.34511451582411E-04_wp, 6.51239772502499E-02_wp, 1.21424722812014E-04_wp, &
      &-2.51829922520101E-02_wp, 1.91223197619631E-03_wp,-1.05936624101185E-03_wp, &
      & 1.27335936243968E-03_wp, 3.48809823803027E-04_wp, 1.31361041814418E-03_wp, &
      & 5.34522432730452E-05_wp, 2.45086039414661E-02_wp,-2.56042707975578E-03_wp, &
      & 4.06759347016995E-02_wp,-4.26658981254059E-04_wp,-4.73464958472976E-04_wp, &
      & 6.75174793891442E-04_wp,-5.61493875295284E-04_wp,-1.60711065190090E-03_wp, &
      &-3.79091162487751E-04_wp,-1.23906888668214E-03_wp, 1.49160580801637E-03_wp, &
      & 6.75003273155266E-04_wp, 4.64033606839685E-03_wp,-1.76046211725149E-03_wp, &
      &-1.74606285771953E-03_wp,-5.05146906110698E-03_wp, 2.11987439578359E-03_wp, &
      &-1.79484074665602E-03_wp, 4.96984189818760E-04_wp, 5.57587921198552E-03_wp, &
      &-1.29796730481245E-04_wp,-7.71056435516257E-04_wp,-1.85520908913895E-04_wp, &
      &-9.53869418532134E-04_wp, 3.34454098522262E-05_wp,-2.14763237024395E-02_wp, &
      & 4.22704389236198E-02_wp,-1.51686553113530E-02_wp,-7.82696981268825E-04_wp, &
      & 8.26115324539116E-04_wp, 1.18017507170875E-03_wp, 5.55319872616625E-04_wp, &
      & 3.30610206698993E-04_wp,-2.20045263579465E-04_wp,-3.86345949249431E-03_wp, &
      & 4.91636570798305E-02_wp, 3.71609952624270E-04_wp,-1.49962948827029E-04_wp, &
      &-1.27713980885142E-03_wp,-1.02227591396376E-03_wp, 9.54382140678933E-04_wp, &
      &-4.76500934864201E-05_wp, 5.47853281429701E-04_wp,-6.37022577896441E-03_wp, &
      &-1.46041480869821E-02_wp,-4.95678742246128E-03_wp, 1.86282860648305E-04_wp, &
      &-2.24408715978599E-05_wp, 1.12001261850226E-04_wp, 5.82371675898028E-04_wp, &
      & 4.18386404073427E-04_wp,-2.36936520816196E-03_wp,-2.62937520452631E-03_wp, &
      &-2.93069153287612E-03_wp,-1.81229618483608E-03_wp, 1.33984180274338E-02_wp, &
      & 1.32231457493828E-03_wp, 2.05719732236316E-04_wp, 6.81728944109617E-04_wp, &
      &-4.97769753948644E-03_wp, 1.08426219074813E-03_wp, 1.15257191840106E-03_wp, &
      &-1.29796730481245E-04_wp, 9.94185722342727E-03_wp,-7.18121361524183E-05_wp, &
      &-1.96699029239134E-04_wp,-2.73593329213702E-03_wp, 1.26668964141692E-03_wp, &
      & 4.19015268410830E-02_wp, 2.73109384918306E-02_wp, 2.86043258072822E-02_wp, &
      & 1.52385694306641E-03_wp, 1.17176477545813E-03_wp,-5.69308111565666E-04_wp, &
      & 8.01133922692052E-04_wp,-5.92917008182626E-04_wp,-8.92853894712445E-04_wp, &
      & 4.00952224969050E-02_wp, 6.11598468630614E-03_wp,-3.43458365253901E-02_wp, &
      & 1.81169953657323E-03_wp,-6.03449882052736E-04_wp, 9.59524174131419E-04_wp, &
      & 6.67147821678755E-04_wp, 1.39319388056191E-04_wp,-8.42494585548669E-04_wp, &
      &-1.91489697850969E-02_wp, 5.62404928698235E-03_wp, 4.94699954460613E-02_wp, &
      & 1.30383869759727E-03_wp, 4.09728855883975E-04_wp, 9.57204049530645E-04_wp, &
      &-7.51198289646634E-04_wp,-1.22239368337505E-03_wp,-3.07632771315538E-04_wp, &
      & 6.80225612020106E-04_wp, 1.04587611852477E-03_wp,-1.71562995623376E-03_wp, &
      & 7.59055857447997E-04_wp,-5.03958332987170E-03_wp,-9.55018724841159E-04_wp, &
      & 2.02662411163210E-03_wp,-5.83490926431050E-03_wp, 2.95149250539913E-03_wp, &
      & 3.14378953074075E-04_wp,-7.71056435516257E-04_wp,-7.18121361524183E-05_wp, &
      & 6.02070686568459E-03_wp,-1.97514382282722E-04_wp, 3.97736132948771E-04_wp, &
      & 2.34645999821497E-05_wp,-1.51973998328438E-02_wp, 2.87609828406397E-02_wp, &
      &-9.50155837778294E-03_wp,-5.48957450710877E-04_wp, 5.53339536578408E-04_wp, &
      & 8.02603441540097E-04_wp, 3.77924609061101E-04_wp, 1.95211698561655E-04_wp, &
      & 6.69839738764968E-04_wp,-5.69855582943845E-03_wp,-3.24906319103924E-02_wp, &
      &-4.84840539732236E-03_wp,-6.17442835723966E-05_wp, 1.02522766286935E-03_wp, &
      & 4.82528489892734E-04_wp,-5.52941116213947E-04_wp,-4.61446181202303E-04_wp, &
      &-4.40962412742548E-04_wp, 1.09352816143174E-03_wp, 5.78051525646661E-02_wp, &
      &-1.68964519113291E-03_wp,-1.50061773404568E-04_wp, 5.19504171298998E-04_wp, &
      &-1.09043296657485E-03_wp,-1.81994394260434E-03_wp,-3.20958715707679E-05_wp, &
      & 3.36923184300490E-03_wp,-2.77935808450687E-03_wp, 1.49439576750988E-04_wp, &
      & 2.35017057420186E-03_wp, 8.14807186834400E-04_wp,-8.83041595713061E-04_wp, &
      &-6.73908004567400E-04_wp, 9.45853199801574E-03_wp, 5.49862980270049E-03_wp, &
      & 9.39535373805293E-04_wp,-1.87096506480690E-03_wp,-1.85520908913895E-04_wp, &
      &-1.96699029239134E-04_wp,-1.97514382282722E-04_wp, 1.14266944639945E-02_wp, &
      & 3.82848693564779E-03_wp, 4.43424232478341E-04_wp,-4.12247572677397E-03_wp, &
      & 9.67883166460861E-03_wp, 5.65375296258784E-02_wp, 1.60363997930508E-03_wp, &
      & 1.52136956648579E-04_wp,-3.37196377565783E-04_wp, 1.19561506028651E-03_wp, &
      & 1.49002979472088E-03_wp,-9.35609376577390E-04_wp, 7.87102640097848E-03_wp, &
      & 3.89419882648860E-03_wp,-4.71640372580748E-02_wp, 1.35122605955239E-03_wp, &
      &-2.00224351211323E-04_wp, 5.23077778166468E-04_wp, 8.65873372134146E-04_wp, &
      &-8.10925155139404E-04_wp, 1.04270229020181E-03_wp, 3.42456251455934E-02_wp, &
      &-2.45165732131513E-03_wp,-6.15903343214419E-02_wp,-1.62313747856192E-03_wp, &
      &-5.52953471697458E-04_wp,-1.20152057600109E-03_wp, 1.02978264381368E-03_wp, &
      & 1.53238664921058E-03_wp, 8.99317487950987E-03_wp, 6.09598880255074E-04_wp, &
      & 6.08564494858310E-04_wp,-2.96432814230523E-04_wp,-2.97423612012955E-04_wp, &
      & 9.84459994392474E-03_wp,-1.24968367218685E-02_wp, 3.79403591761504E-03_wp, &
      &-4.23716574381990E-04_wp,-8.20550769743827E-04_wp,-8.95788403453741E-03_wp, &
      &-9.53869418532134E-04_wp,-2.73593329213702E-03_wp, 3.97736132948771E-04_wp, &
      & 3.82848693564779E-03_wp, 1.48313519957604E-02_wp, 2.13914454002736E-04_wp, &
      & 2.60340503566586E-02_wp, 4.55835114064785E-02_wp, 3.11676871085774E-02_wp, &
      & 1.24738766350628E-03_wp, 1.53260468361852E-03_wp, 1.46449523729464E-06_wp, &
      & 1.22724085143122E-03_wp,-1.01555733777233E-04_wp, 6.87860581410655E-04_wp, &
      &-3.92830170172995E-02_wp,-4.33396736859480E-02_wp, 3.84108907998124E-02_wp, &
      &-1.77145679916974E-03_wp, 1.76952428228545E-03_wp,-2.92072808132570E-04_wp, &
      &-1.42087409386533E-03_wp,-3.57640032710888E-04_wp, 2.43804604273497E-04_wp, &
      & 4.74548117289432E-02_wp,-2.86052910874011E-02_wp,-2.74882106726929E-02_wp, &
      &-1.80961978070453E-03_wp,-9.88501696530309E-04_wp,-3.96125648447587E-04_wp, &
      & 1.22747267256443E-03_wp, 4.59124605942063E-04_wp,-1.09160503081918E-03_wp, &
      & 7.15239177618886E-03_wp, 4.56133941233273E-03_wp, 4.73945653321463E-03_wp, &
      & 2.23748991423594E-02_wp, 2.07361056877234E-02_wp,-7.48356444965498E-03_wp, &
      & 1.40905795163860E-02_wp,-8.82914903560124E-03_wp, 2.51208188470779E-05_wp, &
      &-9.71288715620663E-04_wp, 3.34454098522262E-05_wp, 1.26668964141692E-03_wp, &
      & 2.34645999821497E-05_wp, 4.43424232478341E-04_wp, 2.13914454002736E-04_wp, &
      & 9.79896488940348E-01_wp,-6.57322929233282E-02_wp,-4.15963736907122E-02_wp, &
      &-4.43149339013263E-02_wp,-3.26278019336732E-03_wp,-3.02728738083426E-03_wp, &
      & 1.02915011544939E-03_wp,-2.06316466007126E-03_wp, 1.29122317026888E-03_wp, &
      &-1.73472742483357E-03_wp,-4.78997962193546E-03_wp,-4.81001363066836E-03_wp, &
      &-2.26440585379319E-02_wp, 1.03307478495128E-03_wp, 8.28947346207256E-04_wp, &
      & 8.36491910921276E-04_wp, 3.56449239284045E-04_wp,-2.26774965634077E-03_wp, &
      &-4.67635439959650E-04_wp,-2.15271189761502E-02_wp,-5.04314644985981E-03_wp, &
      & 2.25472945476261E-03_wp, 2.05817642037157E-03_wp, 5.05903531955649E-04_wp, &
      & 7.91175752116483E-04_wp, 7.16659078003491E-04_wp, 1.18813032100251E-03_wp, &
      & 3.20293314794187E-02_wp,-1.15189576864124E-02_wp,-1.88851883302615E-02_wp, &
      &-1.99644814082213E-02_wp, 6.59864328579052E-02_wp, 6.71002330332418E-02_wp, &
      & 1.50312087642997E-02_wp, 7.86252686133398E-02_wp, 3.37935054072574E-02_wp, &
      &-3.43256354112925E-02_wp,-6.97162425658716E-02_wp,-2.14763237024395E-02_wp, &
      & 4.19015268410830E-02_wp,-1.51973998328438E-02_wp,-4.12247572677397E-03_wp, &
      & 2.60340503566586E-02_wp,-6.57322929233282E-02_wp, 8.97759130716830E-01_wp, &
      &-9.01116304500364E-03_wp,-1.45174628349368E-02_wp, 1.54984343292926E-02_wp, &
      & 1.20456174434359E-02_wp,-1.34358668411596E-02_wp,-3.45682744297207E-04_wp, &
      &-2.38133131656488E-02_wp, 5.97236079805190E-03_wp, 3.58761577868366E-03_wp, &
      &-1.13532555792650E-02_wp,-9.58740704382178E-03_wp,-6.72832129514782E-04_wp, &
      & 2.30563219584726E-03_wp, 4.32903054632106E-04_wp,-1.03486058174053E-03_wp, &
      &-9.31982587104710E-04_wp,-1.75697172800605E-02_wp,-4.80545681158281E-02_wp, &
      &-4.04120334435707E-03_wp, 2.91141309843752E-02_wp, 5.40940931925482E-03_wp, &
      & 2.05121479517989E-03_wp, 2.29185776702489E-03_wp, 4.08333313106208E-04_wp, &
      &-4.29139116977737E-04_wp, 2.32797874206224E-02_wp,-1.67088328726513E-02_wp, &
      & 8.77992495294718E-03_wp,-1.13633091808742E-02_wp, 8.69171101021901E-02_wp, &
      &-1.96682616947918E-02_wp,-9.52279424705729E-02_wp,-1.34685700986416E-02_wp, &
      &-3.37085582367700E-02_wp,-1.00847249630762E-02_wp,-2.42710112695174E-02_wp, &
      & 4.22704389236198E-02_wp, 2.73109384918306E-02_wp, 2.87609828406397E-02_wp, &
      & 9.67883166460861E-03_wp, 4.55835114064785E-02_wp,-4.15963736907122E-02_wp, &
      &-9.01116304500364E-03_wp, 8.98867427012339E-01_wp,-5.85200264439804E-03_wp, &
      &-4.20559005977604E-05_wp, 2.32344411503457E-02_wp, 1.40802881718393E-02_wp, &
      & 1.57179533032114E-02_wp, 4.80630571996852E-05_wp, 1.29675679383981E-03_wp, &
      &-1.01166042844672E-02_wp, 1.64479235418644E-02_wp,-1.13616609265885E-02_wp, &
      &-6.58277296511160E-04_wp,-9.90110157107515E-05_wp,-1.69338125188358E-03_wp, &
      & 8.97862853937416E-04_wp,-2.17183689744087E-03_wp, 1.74301398453031E-03_wp, &
      &-1.23079999997972E-02_wp, 1.53478215818567E-02_wp,-5.76126655951124E-03_wp, &
      & 6.00834559352812E-04_wp, 6.74812882410510E-04_wp,-1.66126653498090E-03_wp, &
      &-2.78697981558906E-04_wp, 2.05970017764324E-03_wp, 2.11137376206241E-02_wp, &
      &-1.96392476573831E-02_wp,-1.26538773814726E-02_wp, 4.51838948509923E-03_wp, &
      & 1.24058369540158E-03_wp, 7.82218702145295E-02_wp, 1.07162329819618E-02_wp, &
      & 5.53966133581774E-03_wp,-8.79319158772149E-02_wp, 3.53479811853906E-02_wp, &
      &-2.43258070526638E-02_wp,-1.51686553113530E-02_wp, 2.86043258072822E-02_wp, &
      &-9.50155837778294E-03_wp, 5.65375296258784E-02_wp, 3.11676871085774E-02_wp, &
      &-4.43149339013263E-02_wp,-1.45174628349368E-02_wp,-5.85200264439804E-03_wp, &
      & 9.09495444174971E-01_wp, 2.36136979255962E-02_wp,-3.50462060520742E-04_wp, &
      &-9.07663574278284E-03_wp, 1.23559391944239E-02_wp, 1.65996363974584E-02_wp, &
      &-2.28367340501708E-02_wp, 1.33600689301378E-02_wp, 1.41674636091019E-03_wp, &
      &-5.45952965985893E-02_wp, 5.34526365153076E-03_wp,-4.98413593713582E-06_wp, &
      & 2.40242825854498E-03_wp, 2.48340276376734E-03_wp,-2.72169924795952E-03_wp, &
      & 1.34726158032522E-02_wp, 1.29098402147179E-02_wp,-9.80559442308669E-03_wp, &
      &-5.05472445563373E-03_wp,-2.40189183443410E-03_wp,-1.16055967334105E-03_wp, &
      &-5.81589715025240E-04_wp, 2.45136476720739E-03_wp, 1.40102216346349E-03_wp, &
      & 1.14014043976149E-03_wp,-8.08703725341765E-04_wp,-7.09842519054558E-04_wp, &
      &-2.45727052670871E-04_wp, 8.43407994918159E-04_wp, 2.84314378706769E-03_wp, &
      & 2.74271011527141E-04_wp, 1.51709946910195E-03_wp,-1.20607118888942E-03_wp, &
      & 3.47939694975228E-04_wp,-1.88935362844706E-03_wp,-7.82696981268825E-04_wp, &
      & 1.52385694306641E-03_wp,-5.48957450710877E-04_wp, 1.60363997930508E-03_wp, &
      & 1.24738766350628E-03_wp,-3.26278019336732E-03_wp, 1.54984343292926E-02_wp, &
      &-4.20559005977604E-05_wp, 2.36136979255962E-02_wp, 9.16664322782240E-04_wp, &
      & 2.19194099720034E-04_wp,-4.79576519455896E-04_wp, 3.36655933977181E-04_wp, &
      & 1.96149013647969E-05_wp,-9.41201078461898E-04_wp,-2.88524317147487E-04_wp, &
      &-1.27997516549154E-03_wp,-5.00411372640482E-03_wp, 1.98057673776146E-04_wp, &
      & 6.86374624017901E-05_wp, 1.12832726753596E-04_wp, 8.04650027810863E-05_wp, &
      &-1.60848353529329E-04_wp, 9.23735814364947E-04_wp,-1.64728607471015E-03_wp, &
      &-1.88915001173799E-03_wp,-5.38651683714022E-04_wp, 6.00370748621287E-05_wp, &
      & 1.06799791194602E-05_wp, 4.33426677929980E-05_wp, 1.15971619292740E-04_wp, &
      & 5.65401791911269E-05_wp, 1.14151243425104E-03_wp,-6.27557723472192E-04_wp, &
      & 2.78132229357562E-06_wp,-6.32157378714537E-04_wp, 2.57087657914212E-03_wp, &
      & 3.82184571884326E-04_wp,-2.28230447764398E-03_wp, 4.93883363269233E-04_wp, &
      &-4.47456238961509E-04_wp,-7.43623597157547E-04_wp,-1.84650959078125E-03_wp, &
      & 8.26115324539116E-04_wp, 1.17176477545813E-03_wp, 5.53339536578408E-04_wp, &
      & 1.52136956648579E-04_wp, 1.53260468361852E-03_wp,-3.02728738083426E-03_wp, &
      & 1.20456174434359E-02_wp, 2.32344411503457E-02_wp,-3.50462060520742E-04_wp, &
      & 2.19194099720034E-04_wp, 7.84936492872693E-04_wp, 1.80915758798787E-04_wp, &
      & 4.08650778681187E-04_wp,-3.38274540502832E-04_wp, 7.13200676898407E-04_wp, &
      &-2.49071523373190E-03_wp, 1.68276351901251E-04_wp,-7.40984775053522E-04_wp, &
      &-5.87717691895947E-05_wp, 5.60579444300148E-05_wp,-6.75347654517853E-05_wp, &
      & 1.75896012197508E-05_wp,-1.27743771739469E-04_wp,-3.13016299573240E-04_wp, &
      &-3.56800317283693E-03_wp, 2.99814085716266E-04_wp,-6.76380607569483E-04_wp, &
      & 1.50521426043470E-04_wp, 8.35092720095804E-05_wp,-1.78734019454987E-05_wp, &
      & 5.07787630967066E-06_wp, 9.85537923192919E-05_wp,-3.01493542097973E-04_wp, &
      & 1.20819065008304E-04_wp, 6.06385637342826E-04_wp, 8.19701869394340E-05_wp, &
      & 4.72280474640535E-04_wp,-1.93416448459594E-03_wp,-1.54891424625151E-03_wp, &
      &-1.33728430469689E-03_wp,-1.80831545521759E-04_wp, 3.68449955830016E-06_wp, &
      & 8.39626505216267E-04_wp, 1.18017507170875E-03_wp,-5.69308111565666E-04_wp, &
      & 8.02603441540097E-04_wp,-3.37196377565783E-04_wp, 1.46449523729464E-06_wp, &
      & 1.02915011544939E-03_wp,-1.34358668411596E-02_wp, 1.40802881718393E-02_wp, &
      &-9.07663574278284E-03_wp,-4.79576519455896E-04_wp, 1.80915758798787E-04_wp, &
      & 5.17997844063907E-04_wp, 1.22196058566877E-04_wp, 1.92582814986119E-04_wp, &
      & 6.94244409194631E-04_wp,-7.43888213748894E-04_wp, 1.89961156823768E-03_wp, &
      & 1.41680520137045E-03_wp,-8.73158437683607E-05_wp,-6.27091019265727E-05_wp, &
      &-9.95063330636391E-05_wp, 1.69169056622561E-05_wp, 1.01932464729247E-05_wp, &
      & 6.86048819619154E-04_wp, 1.12236025836381E-03_wp, 1.84140109929851E-03_wp, &
      &-1.05661045734916E-03_wp,-7.39580970193152E-05_wp,-5.70885132748073E-06_wp, &
      &-9.38004786692939E-05_wp,-6.43822098136181E-05_wp, 4.02430278647076E-05_wp, &
      & 7.61823161502331E-04_wp,-6.29338096627950E-04_wp, 3.06228686912739E-06_wp, &
      &-1.19768297628478E-04_wp, 1.34942396110424E-03_wp, 4.83643649261318E-04_wp, &
      &-1.52565429299107E-03_wp, 9.51914809403550E-06_wp,-1.30457359386324E-03_wp, &
      & 3.27464191561612E-04_wp,-8.32108659053711E-04_wp, 5.55319872616625E-04_wp, &
      & 8.01133922692052E-04_wp, 3.77924609061101E-04_wp, 1.19561506028651E-03_wp, &
      & 1.22724085143122E-03_wp,-2.06316466007126E-03_wp,-3.45682744297207E-04_wp, &
      & 1.57179533032114E-02_wp, 1.23559391944239E-02_wp, 3.36655933977181E-04_wp, &
      & 4.08650778681187E-04_wp, 1.22196058566877E-04_wp, 4.58450887444042E-04_wp, &
      & 2.38862116773585E-04_wp,-7.44842386190068E-04_wp,-8.81280082986540E-04_wp, &
      & 4.46531921382866E-04_wp,-3.31795027746360E-03_wp, 9.74089105250731E-05_wp, &
      & 4.58190789078680E-06_wp, 8.56170943489919E-06_wp, 8.84178093731221E-05_wp, &
      &-1.38692174983584E-04_wp, 9.28641771788832E-04_wp, 1.20019615694106E-04_wp, &
      & 5.49103850892938E-06_wp,-2.03354641721527E-03_wp,-4.17039969912064E-05_wp, &
      &-5.12570566395057E-06_wp,-6.69133776173241E-05_wp, 5.12877242948265E-05_wp, &
      & 1.03773299413043E-04_wp,-4.78022561902225E-04_wp,-7.69573321153983E-05_wp, &
      & 2.86568553158027E-04_wp, 7.02087724447307E-04_wp,-1.20043130191039E-03_wp, &
      &-5.30311417466847E-04_wp,-8.33157167079826E-05_wp,-1.49346420012755E-03_wp, &
      &-1.72606096100193E-03_wp, 1.58744109760933E-03_wp, 1.57482674045950E-03_wp, &
      & 3.30610206698993E-04_wp,-5.92917008182626E-04_wp, 1.95211698561655E-04_wp, &
      & 1.49002979472088E-03_wp,-1.01555733777233E-04_wp, 1.29122317026888E-03_wp, &
      &-2.38133131656488E-02_wp, 4.80630571996852E-05_wp, 1.65996363974584E-02_wp, &
      & 1.96149013647969E-05_wp,-3.38274540502832E-04_wp, 1.92582814986119E-04_wp, &
      & 2.38862116773585E-04_wp, 9.54162780231959E-04_wp,-2.25215415158003E-03_wp, &
      & 1.15317833242442E-03_wp, 1.46641095484580E-03_wp,-2.97128264733445E-03_wp, &
      & 1.74461196351452E-04_wp,-9.53916270802859E-05_wp, 4.56136924682147E-05_wp, &
      & 1.19018582501939E-04_wp,-4.22197524065099E-05_wp, 2.15520688800641E-03_wp, &
      & 5.25039872285271E-03_wp, 1.28694863261561E-05_wp,-2.17948057370663E-03_wp, &
      &-2.95522378422435E-04_wp,-1.25578374403127E-04_wp,-1.08474731496987E-04_wp, &
      & 5.08839831972418E-05_wp, 3.72929823801503E-05_wp,-1.33276214891788E-03_wp, &
      &-8.03102785109213E-03_wp,-3.81506946938390E-03_wp, 5.48748418783290E-03_wp, &
      &-2.39746394216992E-02_wp, 1.83841041806474E-02_wp,-9.69794623656807E-03_wp, &
      &-1.19096389468461E-02_wp,-8.62164670054123E-03_wp, 5.59683214958685E-04_wp, &
      &-6.34511451582411E-04_wp,-2.20045263579465E-04_wp,-8.92853894712445E-04_wp, &
      & 6.69839738764968E-04_wp,-9.35609376577390E-04_wp, 6.87860581410655E-04_wp, &
      &-1.73472742483357E-03_wp, 5.97236079805190E-03_wp, 1.29675679383981E-03_wp, &
      &-2.28367340501708E-02_wp,-9.41201078461898E-04_wp, 7.13200676898407E-04_wp, &
      & 6.94244409194631E-04_wp,-7.44842386190068E-04_wp,-2.25215415158003E-03_wp, &
      & 9.81042082276303E-01_wp, 6.76631500443054E-02_wp, 3.31263519224785E-02_wp, &
      &-4.62192374315659E-02_wp, 3.56676334999497E-03_wp,-2.66831115267742E-03_wp, &
      & 1.45049718606748E-03_wp, 1.78512181474530E-03_wp, 1.32158505057704E-03_wp, &
      &-3.02442989193302E-03_wp, 1.59234254613220E-02_wp, 1.87596427584668E-02_wp, &
      & 4.69726375095024E-03_wp,-1.69035032089387E-03_wp,-6.62367898981521E-04_wp, &
      &-1.25089497431037E-03_wp,-1.79892494783825E-03_wp,-3.50141211609294E-04_wp, &
      &-3.51386816605225E-02_wp,-1.18519098620367E-02_wp,-1.56903929838757E-02_wp, &
      & 2.07344255347774E-02_wp, 7.91327525317699E-02_wp,-6.46154002416890E-02_wp, &
      & 5.97324438688014E-04_wp, 7.76490812173268E-02_wp,-3.71397599691318E-02_wp, &
      &-4.32983426569359E-02_wp, 6.51239772502499E-02_wp,-3.86345949249431E-03_wp, &
      & 4.00952224969050E-02_wp,-5.69855582943845E-03_wp, 7.87102640097848E-03_wp, &
      &-3.92830170172995E-02_wp,-4.78997962193546E-03_wp, 3.58761577868366E-03_wp, &
      &-1.01166042844672E-02_wp, 1.33600689301378E-02_wp,-2.88524317147487E-04_wp, &
      &-2.49071523373190E-03_wp,-7.43888213748894E-04_wp,-8.81280082986540E-04_wp, &
      & 1.15317833242442E-03_wp, 6.76631500443054E-02_wp, 8.93154835465960E-01_wp, &
      &-4.10699248842775E-03_wp, 1.29570205711521E-02_wp, 1.67088075680187E-02_wp, &
      &-1.36859151158571E-02_wp, 1.27380826419643E-02_wp,-3.28701831792318E-05_wp, &
      & 2.27063670207665E-02_wp, 9.65825023407488E-03_wp,-5.19049746898824E-03_wp, &
      &-3.84519952008319E-02_wp,-1.24584238984903E-02_wp, 1.38573985450491E-03_wp, &
      &-6.85719506260365E-04_wp, 2.13681430858576E-03_wp, 3.54363548421033E-03_wp, &
      &-6.38092426844506E-06_wp,-1.39271236840666E-02_wp,-1.74731534485127E-02_wp, &
      & 9.92777633341959E-03_wp, 1.12782220484343E-02_wp, 6.88040171310118E-02_wp, &
      & 3.76013880919811E-02_wp, 8.53965281948876E-02_wp,-3.00304876167653E-02_wp, &
      & 3.53408962248655E-02_wp,-1.56788026891989E-02_wp, 1.21424722812014E-04_wp, &
      & 4.91636570798305E-02_wp, 6.11598468630614E-03_wp,-3.24906319103924E-02_wp, &
      & 3.89419882648860E-03_wp,-4.33396736859480E-02_wp,-4.81001363066836E-03_wp, &
      &-1.13532555792650E-02_wp, 1.64479235418644E-02_wp, 1.41674636091019E-03_wp, &
      &-1.27997516549154E-03_wp, 1.68276351901251E-04_wp, 1.89961156823768E-03_wp, &
      & 4.46531921382866E-04_wp, 1.46641095484580E-03_wp, 3.31263519224785E-02_wp, &
      &-4.10699248842775E-03_wp, 9.08224713552045E-01_wp, 9.12432844483099E-03_wp, &
      &-3.41613266047853E-04_wp,-2.27460155593839E-02_wp,-1.60912578565575E-02_wp, &
      & 1.70167016886020E-02_wp, 3.01226961416096E-04_wp, 1.87142490100715E-02_wp, &
      &-2.32116270340476E-02_wp,-2.18422945282865E-02_wp,-3.26436291902690E-02_wp, &
      & 1.96946187706925E-03_wp, 5.26688653828019E-04_wp,-1.74728068369742E-05_wp, &
      & 3.65712561056565E-03_wp, 4.03104311606952E-03_wp, 2.23819392409864E-02_wp, &
      & 2.21934571114739E-02_wp, 1.08044883389166E-02_wp, 2.59864135990082E-03_wp, &
      &-1.27823598977582E-02_wp, 6.76244370380171E-02_wp,-6.23157155351273E-03_wp, &
      &-2.01300808356742E-03_wp,-9.54603867530334E-02_wp,-3.91999274242567E-02_wp, &
      &-2.51829922520101E-02_wp, 3.71609952624270E-04_wp,-3.43458365253901E-02_wp, &
      &-4.84840539732236E-03_wp,-4.71640372580748E-02_wp, 3.84108907998124E-02_wp, &
      &-2.26440585379319E-02_wp,-9.58740704382178E-03_wp,-1.13616609265885E-02_wp, &
      &-5.45952965985893E-02_wp,-5.00411372640482E-03_wp,-7.40984775053522E-04_wp, &
      & 1.41680520137045E-03_wp,-3.31795027746360E-03_wp,-2.97128264733445E-03_wp, &
      &-4.62192374315659E-02_wp, 1.29570205711521E-02_wp, 9.12432844483099E-03_wp, &
      & 9.09747652808388E-01_wp,-2.27488843952394E-02_wp,-5.32740547571669E-04_wp, &
      &-1.00490960111794E-02_wp,-1.37415944004160E-02_wp, 1.76756499052769E-02_wp, &
      & 1.35400708983845E-02_wp, 8.14153301903184E-03_wp,-1.12402991020031E-02_wp, &
      &-3.83955742874452E-03_wp,-1.74863714612803E-03_wp,-9.80560844842254E-04_wp, &
      &-3.97380098058104E-04_wp, 2.70852077770555E-03_wp, 1.45358935925310E-03_wp, &
      &-1.29940447168148E-03_wp,-8.42367911860373E-04_wp,-6.58362672011616E-04_wp, &
      & 3.30597234140415E-04_wp, 1.67224810045337E-03_wp,-2.56243044267504E-03_wp, &
      & 2.39631207747973E-04_wp, 1.10201283191560E-03_wp, 1.04148511475808E-03_wp, &
      & 4.04261927191880E-04_wp, 1.91223197619631E-03_wp,-1.49962948827029E-04_wp, &
      & 1.81169953657323E-03_wp,-6.17442835723966E-05_wp, 1.35122605955239E-03_wp, &
      &-1.77145679916974E-03_wp, 1.03307478495128E-03_wp,-6.72832129514782E-04_wp, &
      &-6.58277296511160E-04_wp, 5.34526365153076E-03_wp, 1.98057673776146E-04_wp, &
      &-5.87717691895947E-05_wp,-8.73158437683607E-05_wp, 9.74089105250731E-05_wp, &
      & 1.74461196351452E-04_wp, 3.56676334999497E-03_wp, 1.67088075680187E-02_wp, &
      &-3.41613266047853E-04_wp,-2.27488843952394E-02_wp, 9.21471303420809E-04_wp, &
      &-2.50223042796709E-04_wp, 5.11478506725393E-04_wp, 3.56651093117193E-04_wp, &
      &-2.01068073842728E-05_wp,-1.11114151322445E-03_wp,-1.77564009962320E-03_wp, &
      &-9.69235159515710E-04_wp, 1.79530476408968E-03_wp, 1.18044138863789E-04_wp, &
      & 2.24225409830870E-05_wp, 9.58261722177173E-05_wp,-6.29010293656909E-06_wp, &
      &-7.34189767748855E-05_wp, 9.04158731595418E-04_wp, 6.54572936555637E-04_wp, &
      &-2.66475674848386E-05_wp,-6.86166525721677E-04_wp,-2.26628072621023E-03_wp, &
      & 3.58655471763609E-04_wp,-1.86189149033346E-03_wp,-1.69175844074889E-04_wp, &
      &-3.58853452102961E-04_wp, 1.24003896893423E-03_wp,-1.05936624101185E-03_wp, &
      &-1.27713980885142E-03_wp,-6.03449882052736E-04_wp, 1.02522766286935E-03_wp, &
      &-2.00224351211323E-04_wp, 1.76952428228545E-03_wp, 8.28947346207256E-04_wp, &
      & 2.30563219584726E-03_wp,-9.90110157107515E-05_wp,-4.98413593713582E-06_wp, &
      & 6.86374624017901E-05_wp, 5.60579444300148E-05_wp,-6.27091019265727E-05_wp, &
      & 4.58190789078680E-06_wp,-9.53916270802859E-05_wp,-2.66831115267742E-03_wp, &
      &-1.36859151158571E-02_wp,-2.27460155593839E-02_wp,-5.32740547571669E-04_wp, &
      &-2.50223042796709E-04_wp, 8.06895419612611E-04_wp, 2.13209610966298E-04_wp, &
      &-4.34071103274985E-04_wp,-3.77193165303000E-04_wp,-1.92654310539002E-03_wp, &
      & 3.14504809743806E-03_wp, 3.64690912288282E-03_wp, 2.20085092555298E-03_wp, &
      &-1.22517292249389E-04_wp,-1.16285951506442E-05_wp,-6.52491923975451E-05_wp, &
      &-2.32417402929102E-04_wp,-1.56685849054832E-04_wp,-5.82863124169693E-04_wp, &
      &-9.79001755995669E-05_wp,-5.84082021986155E-04_wp, 6.16125688124352E-05_wp, &
      & 1.02410092154246E-04_wp,-2.02294067435815E-03_wp,-1.22978146596578E-03_wp, &
      & 1.46303970828171E-03_wp,-3.58903840536546E-04_wp, 2.73469257707459E-04_wp, &
      & 1.27335936243968E-03_wp,-1.02227591396376E-03_wp, 9.59524174131419E-04_wp, &
      & 4.82528489892734E-04_wp, 5.23077778166468E-04_wp,-2.92072808132570E-04_wp, &
      & 8.36491910921276E-04_wp, 4.32903054632106E-04_wp,-1.69338125188358E-03_wp, &
      & 2.40242825854498E-03_wp, 1.12832726753596E-04_wp,-6.75347654517853E-05_wp, &
      &-9.95063330636391E-05_wp, 8.56170943489919E-06_wp, 4.56136924682147E-05_wp, &
      & 1.45049718606748E-03_wp, 1.27380826419643E-02_wp,-1.60912578565575E-02_wp, &
      &-1.00490960111794E-02_wp, 5.11478506725393E-04_wp, 2.13209610966298E-04_wp, &
      & 5.86625181050966E-04_wp,-1.45137809114396E-04_wp, 1.20314624517607E-04_wp, &
      &-1.26165357059354E-03_wp, 3.43704852993555E-04_wp, 4.41167570974510E-06_wp, &
      & 2.16340103360513E-03_wp, 1.38248755681475E-05_wp,-1.14143379910040E-05_wp, &
      & 6.32898273021533E-05_wp,-6.57957789189623E-05_wp,-1.36875188754333E-04_wp, &
      &-5.86304925362903E-04_wp,-7.38113927657857E-04_wp, 1.60899708568502E-05_wp, &
      & 2.02646772793093E-04_wp, 1.25422685598679E-03_wp,-3.48741034904623E-04_wp, &
      & 1.62253168231263E-03_wp,-6.46390412013604E-04_wp, 1.79457630154229E-03_wp, &
      & 2.83112352296633E-04_wp, 3.48809823803027E-04_wp, 9.54382140678933E-04_wp, &
      & 6.67147821678755E-04_wp,-5.52941116213947E-04_wp, 8.65873372134146E-04_wp, &
      &-1.42087409386533E-03_wp, 3.56449239284045E-04_wp,-1.03486058174053E-03_wp, &
      & 8.97862853937416E-04_wp, 2.48340276376734E-03_wp, 8.04650027810863E-05_wp, &
      & 1.75896012197508E-05_wp, 1.69169056622561E-05_wp, 8.84178093731221E-05_wp, &
      & 1.19018582501939E-04_wp, 1.78512181474530E-03_wp,-3.28701831792318E-05_wp, &
      & 1.70167016886020E-02_wp,-1.37415944004160E-02_wp, 3.56651093117193E-04_wp, &
      &-4.34071103274985E-04_wp,-1.45137809114396E-04_wp, 5.40221254362377E-04_wp, &
      &-2.59624456923779E-04_wp,-4.07410755353966E-05_wp,-2.31771957076137E-03_wp, &
      &-6.97059680648043E-04_wp,-1.10499628597321E-03_wp, 1.02396778243273E-04_wp, &
      & 4.30213896140587E-05_wp, 1.05577585074867E-05_wp, 5.27478632782136E-05_wp, &
      & 9.04761672193272E-05_wp,-5.20787144981459E-04_wp, 1.66609509564965E-04_wp, &
      &-1.95975352223134E-04_wp, 6.65298252647229E-04_wp, 1.23413766856317E-03_wp, &
      &-5.80912056614872E-04_wp,-2.57581016819968E-04_wp, 1.68271214865759E-03_wp, &
      &-2.19066197265030E-03_wp,-2.02065592692274E-03_wp, 1.31361041814418E-03_wp, &
      &-4.76500934864201E-05_wp, 1.39319388056191E-04_wp,-4.61446181202303E-04_wp, &
      &-8.10925155139404E-04_wp,-3.57640032710888E-04_wp,-2.26774965634077E-03_wp, &
      &-9.31982587104710E-04_wp,-2.17183689744087E-03_wp,-2.72169924795952E-03_wp, &
      &-1.60848353529329E-04_wp,-1.27743771739469E-04_wp, 1.01932464729247E-05_wp, &
      &-1.38692174983584E-04_wp,-4.22197524065099E-05_wp, 1.32158505057704E-03_wp, &
      & 2.27063670207665E-02_wp, 3.01226961416096E-04_wp, 1.76756499052769E-02_wp, &
      &-2.01068073842728E-05_wp,-3.77193165303000E-04_wp, 1.20314624517607E-04_wp, &
      &-2.59624456923779E-04_wp, 9.33862912780774E-04_wp, 1.33668788928570E-03_wp, &
      &-3.22448649333838E-04_wp,-4.38999879664258E-03_wp,-8.07257395393632E-04_wp, &
      & 5.60366974513104E-06_wp,-6.18831150062864E-05_wp, 9.70229349338622E-05_wp, &
      & 2.32508819618581E-04_wp, 3.39849124064269E-05_wp,-1.20082797898505E-03_wp, &
      & 2.55368914978415E-03_wp,-3.62978936655017E-03_wp,-9.07158901302496E-03_wp, &
      &-1.35647385583807E-02_wp,-5.21088318946672E-03_wp,-9.88315126715304E-03_wp, &
      & 2.13613175487847E-02_wp, 2.19877649688146E-02_wp,-8.80011993355233E-04_wp, &
      & 5.34522432730452E-05_wp, 5.47853281429701E-04_wp,-8.42494585548669E-04_wp, &
      &-4.40962412742548E-04_wp, 1.04270229020181E-03_wp, 2.43804604273497E-04_wp, &
      &-4.67635439959650E-04_wp,-1.75697172800605E-02_wp, 1.74301398453031E-03_wp, &
      & 1.34726158032522E-02_wp, 9.23735814364947E-04_wp,-3.13016299573240E-04_wp, &
      & 6.86048819619154E-04_wp, 9.28641771788832E-04_wp, 2.15520688800641E-03_wp, &
      &-3.02442989193302E-03_wp, 9.65825023407488E-03_wp, 1.87142490100715E-02_wp, &
      & 1.35400708983845E-02_wp,-1.11114151322445E-03_wp,-1.92654310539002E-03_wp, &
      &-1.26165357059354E-03_wp,-4.07410755353966E-05_wp, 1.33668788928570E-03_wp, &
      & 9.80347237320379E-01_wp,-2.18568961272081E-02_wp, 3.31155807207264E-02_wp, &
      & 7.98627581118252E-02_wp, 1.97136514590489E-03_wp, 7.92991091054886E-04_wp, &
      & 1.47395046707294E-03_wp,-3.09029645057104E-03_wp,-3.26271386811814E-03_wp, &
      & 1.01678832885046E-02_wp, 1.44150017668853E-02_wp, 5.13925313958751E-03_wp, &
      & 1.33594273132114E-02_wp, 7.56592969150405E-02_wp, 3.66481656934027E-02_wp, &
      &-6.24835621250338E-03_wp, 3.25299748661382E-02_wp, 7.20879127329431E-02_wp, &
      & 3.13091822106836E-02_wp, 2.45086039414661E-02_wp,-6.37022577896441E-03_wp, &
      &-1.91489697850969E-02_wp, 1.09352816143174E-03_wp, 3.42456251455934E-02_wp, &
      & 4.74548117289432E-02_wp,-2.15271189761502E-02_wp,-4.80545681158281E-02_wp, &
      &-1.23079999997972E-02_wp, 1.29098402147179E-02_wp,-1.64728607471015E-03_wp, &
      &-3.56800317283693E-03_wp, 1.12236025836381E-03_wp, 1.20019615694106E-04_wp, &
      & 5.25039872285271E-03_wp, 1.59234254613220E-02_wp,-5.19049746898824E-03_wp, &
      &-2.32116270340476E-02_wp, 8.14153301903184E-03_wp,-1.77564009962320E-03_wp, &
      & 3.14504809743806E-03_wp, 3.43704852993555E-04_wp,-2.31771957076137E-03_wp, &
      &-3.22448649333838E-04_wp,-2.18568961272081E-02_wp, 9.14808340354401E-01_wp, &
      & 7.48068306522445E-03_wp, 5.09969621287088E-03_wp,-2.71571357928294E-02_wp, &
      &-1.35284878636377E-02_wp,-5.25811997546465E-03_wp,-4.54836245032895E-04_wp, &
      &-8.92063582671482E-03_wp,-1.38192679518940E-02_wp, 4.82113542118391E-03_wp, &
      & 1.02205843520703E-02_wp,-1.99274198938198E-02_wp, 3.15195489468795E-02_wp, &
      &-1.62019419351131E-02_wp, 8.60672941363428E-02_wp, 4.67415065569315E-02_wp, &
      &-7.10369548170602E-02_wp, 3.67758849915333E-02_wp,-2.56042707975578E-03_wp, &
      &-1.46041480869821E-02_wp, 5.62404928698235E-03_wp, 5.78051525646661E-02_wp, &
      &-2.45165732131513E-03_wp,-2.86052910874011E-02_wp,-5.04314644985981E-03_wp, &
      &-4.04120334435707E-03_wp, 1.53478215818567E-02_wp,-9.80559442308669E-03_wp, &
      &-1.88915001173799E-03_wp, 2.99814085716266E-04_wp, 1.84140109929851E-03_wp, &
      & 5.49103850892938E-06_wp, 1.28694863261561E-05_wp, 1.87596427584668E-02_wp, &
      &-3.84519952008319E-02_wp,-2.18422945282865E-02_wp,-1.12402991020031E-02_wp, &
      &-9.69235159515710E-04_wp, 3.64690912288282E-03_wp, 4.41167570974510E-06_wp, &
      &-6.97059680648043E-04_wp,-4.38999879664258E-03_wp, 3.31155807207264E-02_wp, &
      & 7.48068306522445E-03_wp, 9.07039936210219E-01_wp,-6.65624962293605E-03_wp, &
      &-4.68412144489354E-04_wp, 8.51223172707439E-03_wp,-1.56557272168524E-02_wp, &
      &-2.67310470185593E-02_wp,-5.13493806965360E-05_wp,-4.08823657751133E-02_wp, &
      & 1.15441408104976E-02_wp,-1.81593556200167E-02_wp,-2.34198218006705E-02_wp, &
      & 6.45164569511301E-02_wp, 4.22691260020166E-02_wp, 3.14814051825769E-03_wp, &
      &-1.03864191994931E-01_wp,-4.72439787083179E-02_wp, 5.10602082770450E-02_wp, &
      & 4.06759347016995E-02_wp,-4.95678742246128E-03_wp, 4.94699954460613E-02_wp, &
      &-1.68964519113291E-03_wp,-6.15903343214419E-02_wp,-2.74882106726929E-02_wp, &
      & 2.25472945476261E-03_wp, 2.91141309843752E-02_wp,-5.76126655951124E-03_wp, &
      &-5.05472445563373E-03_wp,-5.38651683714022E-04_wp,-6.76380607569483E-04_wp, &
      &-1.05661045734916E-03_wp,-2.03354641721527E-03_wp,-2.17948057370663E-03_wp, &
      & 4.69726375095024E-03_wp,-1.24584238984903E-02_wp,-3.26436291902690E-02_wp, &
      &-3.83955742874452E-03_wp, 1.79530476408968E-03_wp, 2.20085092555298E-03_wp, &
      & 2.16340103360513E-03_wp,-1.10499628597321E-03_wp,-8.07257395393632E-04_wp, &
      & 7.98627581118252E-02_wp, 5.09969621287088E-03_wp,-6.65624962293605E-03_wp, &
      & 8.85714398346155E-01_wp, 8.42605738546966E-03_wp, 4.42619688345430E-05_wp, &
      & 1.51188717216053E-02_wp,-1.31893721640515E-02_wp,-2.63691965213547E-02_wp, &
      &-7.15273250933024E-04_wp,-3.91702612956306E-04_wp,-3.90295021261056E-04_wp, &
      &-6.68613327789153E-04_wp,-9.39640522978982E-04_wp,-7.23873876527143E-04_wp, &
      & 3.77981999348335E-04_wp,-1.50087755722131E-03_wp,-2.04362405293725E-03_wp, &
      &-6.95431065144622E-04_wp,-4.26658981254059E-04_wp, 1.86282860648305E-04_wp, &
      & 1.30383869759727E-03_wp,-1.50061773404568E-04_wp,-1.62313747856192E-03_wp, &
      &-1.80961978070453E-03_wp, 2.05817642037157E-03_wp, 5.40940931925482E-03_wp, &
      & 6.00834559352812E-04_wp,-2.40189183443410E-03_wp, 6.00370748621287E-05_wp, &
      & 1.50521426043470E-04_wp,-7.39580970193152E-05_wp,-4.17039969912064E-05_wp, &
      &-2.95522378422435E-04_wp,-1.69035032089387E-03_wp, 1.38573985450491E-03_wp, &
      & 1.96946187706925E-03_wp,-1.74863714612803E-03_wp, 1.18044138863789E-04_wp, &
      &-1.22517292249389E-04_wp, 1.38248755681475E-05_wp, 1.02396778243273E-04_wp, &
      & 5.60366974513104E-06_wp, 1.97136514590489E-03_wp,-2.71571357928294E-02_wp, &
      &-4.68412144489354E-04_wp, 8.42605738546966E-03_wp, 9.22078704785179E-04_wp, &
      & 4.11819321249521E-04_wp, 3.19311058760809E-04_wp,-1.03752848569467E-04_wp, &
      & 1.06589860594684E-05_wp,-2.62747725165519E-04_wp,-1.85560188597899E-04_wp, &
      & 9.37999895527612E-06_wp,-4.37632534503615E-04_wp,-6.34568358829660E-04_wp, &
      &-6.94871573164942E-04_wp, 9.34039727567837E-04_wp, 3.62418390726754E-05_wp, &
      &-1.49241311598283E-03_wp,-1.53919581309734E-04_wp,-4.73464958472976E-04_wp, &
      &-2.24408715978599E-05_wp, 4.09728855883975E-04_wp, 5.19504171298998E-04_wp, &
      &-5.52953471697458E-04_wp,-9.88501696530309E-04_wp, 5.05903531955649E-04_wp, &
      & 2.05121479517989E-03_wp, 6.74812882410510E-04_wp,-1.16055967334105E-03_wp, &
      & 1.06799791194602E-05_wp, 8.35092720095804E-05_wp,-5.70885132748073E-06_wp, &
      &-5.12570566395057E-06_wp,-1.25578374403127E-04_wp,-6.62367898981521E-04_wp, &
      &-6.85719506260365E-04_wp, 5.26688653828019E-04_wp,-9.80560844842254E-04_wp, &
      & 2.24225409830870E-05_wp,-1.16285951506442E-05_wp,-1.14143379910040E-05_wp, &
      & 4.30213896140587E-05_wp,-6.18831150062864E-05_wp, 7.92991091054886E-04_wp, &
      &-1.35284878636377E-02_wp, 8.51223172707439E-03_wp, 4.42619688345430E-05_wp, &
      & 4.11819321249521E-04_wp, 2.86722191709918E-04_wp,-6.45162873531794E-05_wp, &
      &-2.51941036738604E-04_wp, 1.26351280783451E-04_wp,-5.93120928958409E-04_wp, &
      & 2.73984706658681E-05_wp,-5.75461930289809E-04_wp,-1.19585651284702E-04_wp, &
      & 3.16393698227943E-04_wp, 6.92791471894560E-04_wp,-1.19873933534838E-03_wp, &
      &-2.39989499316024E-03_wp, 1.88746416191365E-04_wp,-6.53037749246239E-05_wp, &
      & 6.75174793891442E-04_wp, 1.12001261850226E-04_wp, 9.57204049530645E-04_wp, &
      &-1.09043296657485E-03_wp,-1.20152057600109E-03_wp,-3.96125648447587E-04_wp, &
      & 7.91175752116483E-04_wp, 2.29185776702489E-03_wp,-1.66126653498090E-03_wp, &
      &-5.81589715025240E-04_wp, 4.33426677929980E-05_wp,-1.78734019454987E-05_wp, &
      &-9.38004786692939E-05_wp,-6.69133776173241E-05_wp,-1.08474731496987E-04_wp, &
      &-1.25089497431037E-03_wp, 2.13681430858576E-03_wp,-1.74728068369742E-05_wp, &
      &-3.97380098058104E-04_wp, 9.58261722177173E-05_wp,-6.52491923975451E-05_wp, &
      & 6.32898273021533E-05_wp, 1.05577585074867E-05_wp, 9.70229349338622E-05_wp, &
      & 1.47395046707294E-03_wp,-5.25811997546465E-03_wp,-1.56557272168524E-02_wp, &
      & 1.51188717216053E-02_wp, 3.19311058760809E-04_wp,-6.45162873531794E-05_wp, &
      & 5.64230343817690E-04_wp, 2.39346141037040E-04_wp,-4.00869699204397E-04_wp, &
      & 1.02934066632207E-03_wp,-3.83899731674412E-04_wp,-4.74477891207914E-05_wp, &
      & 1.01952218604768E-03_wp,-1.37698038383402E-03_wp, 2.14712344443191E-04_wp, &
      &-2.29602091271526E-03_wp, 3.71610089264428E-04_wp, 2.32194947529299E-03_wp, &
      &-2.04327967176310E-03_wp,-5.61493875295284E-04_wp, 5.82371675898028E-04_wp, &
      &-7.51198289646634E-04_wp,-1.81994394260434E-03_wp, 1.02978264381368E-03_wp, &
      & 1.22747267256443E-03_wp, 7.16659078003491E-04_wp, 4.08333313106208E-04_wp, &
      &-2.78697981558906E-04_wp, 2.45136476720739E-03_wp, 1.15971619292740E-04_wp, &
      & 5.07787630967066E-06_wp,-6.43822098136181E-05_wp, 5.12877242948265E-05_wp, &
      & 5.08839831972418E-05_wp,-1.79892494783825E-03_wp, 3.54363548421033E-03_wp, &
      & 3.65712561056565E-03_wp, 2.70852077770555E-03_wp,-6.29010293656909E-06_wp, &
      &-2.32417402929102E-04_wp,-6.57957789189623E-05_wp, 5.27478632782136E-05_wp, &
      & 2.32508819618581E-04_wp,-3.09029645057104E-03_wp,-4.54836245032895E-04_wp, &
      &-2.67310470185593E-02_wp,-1.31893721640515E-02_wp,-1.03752848569467E-04_wp, &
      &-2.51941036738604E-04_wp, 2.39346141037040E-04_wp, 1.01787218154808E-03_wp, &
      & 4.18069420263225E-04_wp, 1.20857770847103E-03_wp,-5.51201269388917E-04_wp, &
      & 5.45960716411972E-04_wp, 5.93888712235277E-04_wp,-2.23650440108797E-03_wp, &
      &-1.19295902152296E-03_wp, 4.64422684027194E-05_wp, 2.60085995965948E-03_wp, &
      & 4.24533792527953E-04_wp,-1.86505232958969E-03_wp,-1.60711065190090E-03_wp, &
      & 4.18386404073427E-04_wp,-1.22239368337505E-03_wp,-3.20958715707679E-05_wp, &
      & 1.53238664921058E-03_wp, 4.59124605942063E-04_wp, 1.18813032100251E-03_wp, &
      &-4.29139116977737E-04_wp, 2.05970017764324E-03_wp, 1.40102216346349E-03_wp, &
      & 5.65401791911269E-05_wp, 9.85537923192919E-05_wp, 4.02430278647076E-05_wp, &
      & 1.03773299413043E-04_wp, 3.72929823801503E-05_wp,-3.50141211609294E-04_wp, &
      &-6.38092426844506E-06_wp, 4.03104311606952E-03_wp, 1.45358935925310E-03_wp, &
      &-7.34189767748855E-05_wp,-1.56685849054832E-04_wp,-1.36875188754333E-04_wp, &
      & 9.04761672193272E-05_wp, 3.39849124064269E-05_wp,-3.26271386811814E-03_wp, &
      &-8.92063582671482E-03_wp,-5.13493806965360E-05_wp,-2.63691965213547E-02_wp, &
      & 1.06589860594684E-05_wp, 1.26351280783451E-04_wp,-4.00869699204397E-04_wp, &
      & 4.18069420263225E-04_wp, 8.89206033352750E-04_wp], shape(density))

   call get_structure(mol, "f-block", "CeCl3")
   call test_energy_generic(error, mol, density, qsh, make_exchange_gxtb, &
      & -3.45602684277503_wp, thr_in=thr1)

end subroutine test_e_fock_cecl3


subroutine test_e_fock_ce2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(8) = [&
      &-8.71550786234335E-02_wp,-7.18166762937467E-02_wp, 4.08861943259277E-01_wp, &
      &-2.49890188342099E-01_wp,-8.71550786235603E-02_wp,-7.18166762939041E-02_wp, &
      & 4.08861943259043E-01_wp,-2.49890188341574E-01_wp]
   
   real(wp), parameter :: density(32, 32, 1) = reshape([&
      & 6.36026007398742E-01_wp, 1.99623320104828E-17_wp, 6.78407215513985E-03_wp, &
      &-1.84497744291860E-18_wp, 9.65107976232738E-18_wp,-3.47391173958214E-17_wp, &
      & 5.16421460894714E-03_wp, 3.11958303249284E-18_wp,-2.44860247637404E-17_wp, &
      & 7.50889325836972E-19_wp, 6.74062249596675E-19_wp,-7.47257373407983E-18_wp, &
      &-1.03520480376254E-02_wp, 2.50753660983099E-18_wp,-1.95160341757053E-17_wp, &
      & 8.83722844825746E-30_wp, 6.36026007398848E-01_wp, 4.54155077637372E-17_wp, &
      &-6.78407215531017E-03_wp,-8.67695391478395E-18_wp, 1.89566830830589E-18_wp, &
      & 6.32263657180051E-17_wp, 5.16421460896905E-03_wp, 4.11665386035113E-18_wp, &
      &-1.86487216170029E-17_wp,-1.47266519840908E-17_wp, 8.39838980502871E-18_wp, &
      &-6.88288516982198E-17_wp, 1.03520480376186E-02_wp,-3.50607319354894E-18_wp, &
      & 3.61719196879460E-17_wp, 2.18504335499174E-31_wp, 1.99623320104828E-17_wp, &
      & 4.87261607756642E-02_wp, 3.74130208213754E-18_wp, 6.14530855148644E-17_wp, &
      & 7.50730114572648E-18_wp, 1.18968197226871E-01_wp,-9.86551732611997E-18_wp, &
      & 1.56125112837913E-17_wp, 2.35089931305409E-17_wp, 1.59613420193876E-16_wp, &
      &-3.59630392012243E-18_wp, 1.45926928037620E-01_wp,-5.44270283370436E-17_wp, &
      & 3.72122621352501E-17_wp,-2.77321013775153E-17_wp,-3.86489796446206E-20_wp, &
      & 3.80449278353769E-17_wp, 4.87261607756792E-02_wp,-6.49358823473300E-18_wp, &
      & 3.61801788060898E-17_wp,-5.80448036224271E-18_wp,-1.18968197226890E-01_wp, &
      &-3.53811371739650E-17_wp,-4.94396190653390E-17_wp, 1.44775446971247E-17_wp, &
      & 2.23819063966914E-16_wp,-2.42621683411740E-17_wp, 1.45926928037584E-01_wp, &
      & 6.04175268235442E-17_wp, 3.12802712863402E-17_wp,-6.08317383838331E-17_wp, &
      &-3.87667277030004E-20_wp, 6.78407215513985E-03_wp, 3.74130208213754E-18_wp, &
      & 7.23612469785225E-05_wp,-1.72667467338322E-18_wp,-1.15835631094633E-18_wp, &
      & 8.24423022984701E-18_wp, 5.50832892431720E-05_wp,-4.13446974606691E-18_wp, &
      &-7.83164995428246E-17_wp, 8.00924381668769E-21_wp,-1.42745626581918E-18_wp, &
      & 1.04872102637566E-17_wp,-1.10418504941255E-04_wp,-5.08542772380710E-18_wp, &
      &-8.89911141652889E-17_wp, 1.08464729600112E-31_wp, 6.78407215514099E-03_wp, &
      & 4.01279440930580E-18_wp,-7.23612469803392E-05_wp,-1.79954688428111E-18_wp, &
      &-1.24107819697753E-18_wp,-7.94037544109613E-18_wp, 5.50832892434057E-05_wp, &
      & 4.21165393765978E-18_wp,-7.82542368597328E-17_wp,-1.57079535272932E-19_wp, &
      & 1.52422615989048E-18_wp, 9.83276309958546E-18_wp, 1.10418504941183E-04_wp, &
      &-5.14957095434112E-18_wp, 8.91687715689248E-17_wp, 8.44369028099050E-33_wp, &
      &-1.84497744291860E-18_wp, 6.14530855148644E-17_wp,-1.72667467338322E-18_wp, &
      & 4.87261607756647E-02_wp,-3.34403003986037E-18_wp, 1.13416284108215E-16_wp, &
      & 8.23465932535890E-18_wp, 1.18968197226872E-01_wp, 6.77389046461013E-18_wp, &
      & 5.28413089972355E-19_wp,-2.63729759855742E-18_wp, 1.16512675045148E-16_wp, &
      & 8.29262677443047E-18_wp, 1.45926928037621E-01_wp, 2.42484355434170E-17_wp, &
      &-1.04755344235862E-16_wp, 1.32381762691404E-18_wp, 6.09607146462613E-17_wp, &
      & 1.71492294696653E-18_wp, 4.87261607756798E-02_wp,-2.78480522972095E-18_wp, &
      &-9.92670042757072E-17_wp, 1.50120628647457E-18_wp,-1.18968197226890E-01_wp, &
      & 4.25039971179187E-17_wp, 7.03012249162146E-18_wp, 6.26764539114895E-18_wp, &
      & 1.33227792809983E-16_wp,-2.54827460989731E-17_wp, 1.45926928037584E-01_wp, &
      &-1.54452490602129E-17_wp,-1.15024453040306E-16_wp, 9.65107976232738E-18_wp, &
      & 7.50730114572648E-18_wp,-1.15835631094633E-18_wp,-3.34403003986037E-18_wp, &
      & 9.85849302395265E-05_wp, 1.79182358624563E-17_wp, 5.90321317875273E-18_wp, &
      &-8.00717972329512E-18_wp, 6.10092008616651E-03_wp,-1.53474430223198E-19_wp, &
      & 1.12134067473773E-04_wp, 2.26140197094306E-17_wp, 9.63138150357932E-18_wp, &
      &-1.36217502491525E-17_wp, 6.93940730020419E-03_wp,-8.23570112181602E-19_wp, &
      & 1.18022048079148E-17_wp, 7.46102863717653E-18_wp, 2.02291975686579E-18_wp, &
      &-1.88534072005163E-18_wp, 9.85849302395307E-05_wp,-1.73482991883529E-17_wp, &
      & 4.48592584932744E-18_wp, 1.25660253041414E-17_wp, 6.10092008616677E-03_wp, &
      &-2.89719626342310E-19_wp,-1.12134067473766E-04_wp, 2.13296369409181E-17_wp, &
      &-9.96563504787242E-18_wp,-8.47927333238547E-18_wp,-6.93940730020377E-03_wp, &
      &-1.63065940971419E-19_wp,-3.47391173958214E-17_wp, 1.18968197226871E-01_wp, &
      & 8.24423022984701E-18_wp, 1.13416284108215E-16_wp, 1.79182358624563E-17_wp, &
      & 2.90468851354290E-01_wp,-2.47651275898488E-17_wp,-5.36440971760485E-17_wp, &
      & 3.19427628553264E-17_wp, 3.89706895667517E-16_wp,-9.24849599448890E-18_wp, &
      & 3.56290404971976E-01_wp,-1.31528541028460E-16_wp,-2.01842105617732E-17_wp, &
      &-9.66643931501240E-17_wp,-9.43640820409339E-20_wp, 9.41075594169643E-18_wp, &
      & 1.18968197226907E-01_wp,-1.49641218981684E-17_wp, 5.11666666624951E-17_wp, &
      &-1.45833735947789E-17_wp,-2.90468851354336E-01_wp,-8.70632287761936E-17_wp, &
      &-3.12250225675825E-17_wp, 9.89187435175968E-18_wp, 5.46469291265161E-16_wp, &
      &-5.87698368612470E-17_wp, 3.56290404971887E-01_wp, 1.46154745803971E-16_wp, &
      &-3.51816401787238E-17_wp,-1.19570175471162E-16_wp,-9.46515718412620E-20_wp, &
      & 5.16421460894714E-03_wp,-9.86551732611997E-18_wp, 5.50832892431720E-05_wp, &
      & 8.23465932535890E-18_wp, 5.90321317875273E-18_wp,-2.47651275898488E-17_wp, &
      & 4.19308522246379E-05_wp, 2.01673791982782E-17_wp, 3.60271611406113E-16_wp, &
      & 6.09684758336005E-21_wp, 6.63086962479289E-18_wp,-3.00917095362412E-17_wp, &
      &-8.40534775096268E-05_wp, 2.47266891277598E-17_wp, 4.09853647289589E-16_wp, &
      & 5.36581311445888E-33_wp, 5.16421460894800E-03_wp,-9.65885020789207E-18_wp, &
      &-5.50832892445549E-05_wp, 8.17918707691012E-18_wp, 5.84024309787051E-18_wp, &
      & 2.49964298671853E-17_wp, 4.19308522248157E-05_wp,-2.01086245564250E-17_wp, &
      & 3.60319007402261E-16_wp,-1.19573084170256E-19_wp,-6.55720582976770E-18_wp, &
      &-3.05898919742584E-17_wp, 8.40534775095716E-05_wp, 2.46778616093225E-17_wp, &
      &-4.09718409789389E-16_wp,-2.73269110400166E-32_wp, 3.11958303249284E-18_wp, &
      & 1.56125112837913E-17_wp,-4.13446974606691E-18_wp, 1.18968197226872E-01_wp, &
      &-8.00717972329512E-18_wp,-5.36440971760485E-17_wp, 2.01673791982782E-17_wp, &
      & 2.90468851354290E-01_wp, 2.62854382582260E-17_wp, 1.29015608257129E-18_wp, &
      &-6.25999952264570E-18_wp,-1.20498957592446E-16_wp, 2.01229128945001E-17_wp, &
      & 3.56290404971975E-01_wp, 7.02902406782333E-17_wp,-2.55767215295260E-16_wp, &
      & 1.08564092227418E-17_wp, 1.39766376211368E-17_wp, 4.10577711599016E-18_wp, &
      & 1.18968197226908E-01_wp,-6.64179475456893E-18_wp, 8.58831364103808E-17_wp, &
      & 3.72720093683120E-18_wp,-2.90468851354336E-01_wp, 1.13522895665109E-16_wp, &
      & 1.71645166743767E-17_wp, 1.51237377141199E-17_wp,-7.82644961576860E-17_wp, &
      &-6.20937448582820E-17_wp, 3.56290404971886E-01_wp,-4.87966690713799E-17_wp, &
      &-2.80839934798361E-16_wp,-2.44860247637404E-17_wp, 2.35089931305409E-17_wp, &
      &-7.83164995428246E-17_wp, 6.77389046461013E-18_wp, 6.10092008616651E-03_wp, &
      & 3.19427628553264E-17_wp, 3.60271611406113E-16_wp, 2.62854382582260E-17_wp, &
      & 3.77554924544305E-01_wp,-9.49775215934981E-18_wp, 6.93940730020419E-03_wp, &
      & 7.85047107447312E-17_wp, 6.06156795466151E-16_wp,-2.02927490629492E-16_wp, &
      & 4.29444635007022E-01_wp,-5.09665668735307E-17_wp, 1.08636167102537E-16_wp, &
      & 2.06454228447258E-17_wp, 1.31819935414578E-16_wp, 9.70447558908168E-17_wp, &
      & 6.10092008616677E-03_wp, 3.32772011249550E-18_wp, 2.72562904660025E-16_wp, &
      & 2.55838335693153E-16_wp, 3.77554924544321E-01_wp,-1.79292746205137E-17_wp, &
      &-6.93940730020379E-03_wp,-9.79208467595434E-19_wp,-6.26842047839864E-16_wp, &
      & 1.15314259280036E-16_wp,-4.29444635006996E-01_wp,-1.00913219923683E-17_wp, &
      & 7.50889325836972E-19_wp, 1.59613420193876E-16_wp, 8.00924381668769E-21_wp, &
      & 5.28413089972355E-19_wp,-1.53474430223198E-19_wp, 3.89706895667517E-16_wp, &
      & 6.09684758336005E-21_wp, 1.29015608257129E-18_wp,-9.49775215934981E-18_wp, &
      & 5.23094941450829E-31_wp,-1.74567371223344E-19_wp, 4.78016238334604E-16_wp, &
      &-1.22215794349642E-20_wp, 1.58251127786412E-18_wp,-1.08030870326539E-17_wp, &
      & 1.94848126562505E-35_wp, 7.50889325837154E-19_wp, 1.59613420193925E-16_wp, &
      &-8.00924381689913E-21_wp, 5.28413089972434E-19_wp,-1.53474430223249E-19_wp, &
      &-3.89706895667578E-16_wp, 6.09684758330446E-21_wp,-1.29015608257161E-18_wp, &
      &-9.49775215935025E-18_wp, 7.33679203094675E-31_wp, 1.74567371223243E-19_wp, &
      & 4.78016238334485E-16_wp, 1.22215794349762E-20_wp, 1.58251127786369E-18_wp, &
      & 1.08030870326530E-17_wp,-1.12052026646400E-33_wp, 6.74062249596675E-19_wp, &
      &-3.59630392012243E-18_wp,-1.42745626581918E-18_wp,-2.63729759855742E-18_wp, &
      & 1.12134067473773E-04_wp,-9.24849599448890E-18_wp, 6.63086962479289E-18_wp, &
      &-6.25999952264570E-18_wp, 6.93940730020419E-03_wp,-1.74567371223344E-19_wp, &
      & 1.27545346511502E-04_wp,-1.06214851137822E-17_wp, 1.11227818701514E-17_wp, &
      &-1.20009234127102E-17_wp, 7.89313300256411E-03_wp,-9.36758450854265E-19_wp, &
      & 3.12082972633157E-18_wp,-3.64893594596737E-18_wp, 2.41084201803115E-18_wp, &
      &-9.78131575785710E-19_wp, 1.12134067473778E-04_wp, 9.89676259380455E-18_wp, &
      & 5.01879572423427E-18_wp, 1.14453954723415E-17_wp, 6.93940730020449E-03_wp, &
      &-3.29537588045326E-19_wp,-1.27545346511494E-04_wp,-1.20823885568917E-17_wp, &
      &-1.15029739481702E-17_wp,-6.15168405973870E-18_wp,-7.89313300256364E-03_wp, &
      &-1.85477102668095E-19_wp,-7.47257373407983E-18_wp, 1.45926928037620E-01_wp, &
      & 1.04872102637566E-17_wp, 1.16512675045148E-16_wp, 2.26140197094306E-17_wp, &
      & 3.56290404971976E-01_wp,-3.00917095362412E-17_wp,-1.20498957592446E-16_wp, &
      & 7.85047107447312E-17_wp, 4.78016238334604E-16_wp,-1.06214851137822E-17_wp, &
      & 4.37027419922076E-01_wp,-1.61905422978401E-16_wp,-9.16933384556871E-17_wp, &
      &-7.38409255258722E-17_wp,-1.15747409226171E-19_wp, 4.66818603043603E-17_wp, &
      & 1.45926928037665E-01_wp,-1.87298597344789E-17_wp, 4.05535275106961E-17_wp, &
      &-1.72526012441359E-17_wp,-3.56290404972032E-01_wp,-1.06506841820406E-16_wp, &
      & 1.73472347597681E-17_wp, 5.14569913135540E-17_wp, 6.70301700791082E-16_wp, &
      &-7.28101091714507E-17_wp, 4.37027419921967E-01_wp, 1.79845991815207E-16_wp, &
      &-1.09030624000210E-16_wp,-1.91393372967111E-16_wp,-1.16100045513692E-19_wp, &
      &-1.03520480376254E-02_wp,-5.44270283370436E-17_wp,-1.10418504941255E-04_wp, &
      & 8.29262677443047E-18_wp, 9.63138150357932E-18_wp,-1.31528541028460E-16_wp, &
      &-8.40534775096268E-05_wp, 2.01229128945001E-17_wp, 6.06156795466151E-16_wp, &
      &-1.22215794349642E-20_wp, 1.11227818701514E-17_wp,-1.61905422978401E-16_wp, &
      & 1.68491378224602E-04_wp, 2.47043233862882E-17_wp, 6.89329010367911E-16_wp, &
      &-2.43328627619396E-31_wp,-1.03520480376272E-02_wp,-5.48413077769616E-17_wp, &
      & 1.10418504944028E-04_wp, 8.40382497671730E-18_wp, 9.75760966247523E-18_wp, &
      & 1.31064878600120E-16_wp,-8.40534775099834E-05_wp,-2.02406908958029E-17_wp, &
      & 6.06061786705714E-16_wp, 2.39693042421374E-19_wp,-1.12704463656210E-17_wp, &
      &-1.60906779637513E-16_wp,-1.68491378224491E-04_wp, 2.48022017384882E-17_wp, &
      &-6.89600103884371E-16_wp,-3.92090630566152E-32_wp, 2.50753660983099E-18_wp, &
      & 3.72122621352501E-17_wp,-5.08542772380710E-18_wp, 1.45926928037621E-01_wp, &
      &-1.36217502491525E-17_wp,-2.01842105617732E-17_wp, 2.47266891277598E-17_wp, &
      & 3.56290404971975E-01_wp,-2.02927490629492E-16_wp, 1.58251127786412E-18_wp, &
      &-1.20009234127102E-17_wp,-9.16933384556871E-17_wp, 2.47043233862882E-17_wp, &
      & 4.37027419922075E-01_wp,-1.81271777740451E-16_wp,-3.13725221452240E-16_wp, &
      & 1.19975625581186E-17_wp, 3.51920727141113E-17_wp, 5.05023321416519E-18_wp, &
      & 1.45926928037666E-01_wp,-1.19469628797602E-17_wp, 5.94289669462712E-17_wp, &
      & 4.56109222030375E-18_wp,-3.56290404972031E-01_wp,-9.59216313667728E-17_wp, &
      & 2.10540736762259E-17_wp, 2.28732247113919E-17_wp,-4.11883049958516E-17_wp, &
      &-7.61859354340708E-17_wp, 4.37027419921966E-01_wp, 2.07635889748987E-16_wp, &
      &-3.44479532435530E-16_wp,-1.95160341757053E-17_wp,-2.77321013775153E-17_wp, &
      &-8.89911141652889E-17_wp, 2.42484355434170E-17_wp, 6.93940730020419E-03_wp, &
      &-9.66643931501240E-17_wp, 4.09853647289589E-16_wp, 7.02902406782333E-17_wp, &
      & 4.29444635007022E-01_wp,-1.08030870326539E-17_wp, 7.89313300256411E-03_wp, &
      &-7.38409255258722E-17_wp, 6.89329010367911E-16_wp,-1.81271777740451E-16_wp, &
      & 4.88465869592103E-01_wp,-5.79712176578855E-17_wp, 1.31901964586593E-16_wp, &
      &-3.09892298541188E-17_wp, 1.49847858431830E-16_wp, 1.26925784809929E-16_wp, &
      & 6.93940730020448E-03_wp, 1.36782316927511E-16_wp, 3.10090590031171E-16_wp, &
      & 2.50607549961322E-16_wp, 4.29444635007041E-01_wp,-2.03934058194061E-17_wp, &
      &-7.89313300256366E-03_wp,-1.64248811785728E-16_wp,-7.12857165023500E-16_wp, &
      & 1.80707905874698E-16_wp,-4.88465869592074E-01_wp,-1.14782348422061E-17_wp, &
      & 8.83722844825746E-30_wp,-3.86489796446206E-20_wp, 1.08464729600112E-31_wp, &
      &-1.04755344235862E-16_wp,-8.23570112181602E-19_wp,-9.43640820409339E-20_wp, &
      & 5.36581311445888E-33_wp,-2.55767215295260E-16_wp,-5.09665668735307E-17_wp, &
      & 1.94848126562505E-35_wp,-9.36758450854265E-19_wp,-1.15747409226171E-19_wp, &
      &-2.43328627619396E-31_wp,-3.13725221452240E-16_wp,-5.79712176578855E-17_wp, &
      & 2.32091367985052E-31_wp, 8.81243126430586E-30_wp,-3.86489796446311E-20_wp, &
      &-1.15659771140108E-31_wp,-1.04755344235895E-16_wp,-8.23570112181639E-19_wp, &
      & 9.43640820409146E-20_wp, 3.17020579613494E-32_wp, 2.55767215295300E-16_wp, &
      &-5.09665668735330E-17_wp,-1.28711492963744E-32_wp, 9.36758450854204E-19_wp, &
      &-1.15747409226165E-19_wp, 2.83072876256128E-31_wp,-3.13725221452161E-16_wp, &
      & 5.79712176578819E-17_wp, 2.48650912671926E-31_wp, 6.36026007398848E-01_wp, &
      & 3.80449278353769E-17_wp, 6.78407215514099E-03_wp, 1.32381762691404E-18_wp, &
      & 1.18022048079148E-17_wp, 9.41075594169643E-18_wp, 5.16421460894800E-03_wp, &
      & 1.08564092227418E-17_wp, 1.08636167102537E-16_wp, 7.50889325837154E-19_wp, &
      & 3.12082972633157E-18_wp, 4.66818603043603E-17_wp,-1.03520480376272E-02_wp, &
      & 1.19975625581186E-17_wp, 1.31901964586593E-16_wp, 8.81243126430586E-30_wp, &
      & 6.36026007398954E-01_wp, 6.34981035886410E-17_wp,-6.78407215531131E-03_wp, &
      &-5.50815884495145E-18_wp, 4.04679335389210E-18_wp, 1.90764923804851E-17_wp, &
      & 5.16421460896991E-03_wp,-3.62017232989779E-18_wp, 1.14473470249282E-16_wp, &
      &-1.47266519840931E-17_wp, 5.95162232829546E-18_wp,-1.46744176598035E-17_wp, &
      & 1.03520480376204E-02_wp, 5.98395275473543E-18_wp,-1.15246079074340E-16_wp, &
      & 2.07451493237937E-31_wp, 4.54155077637372E-17_wp, 4.87261607756792E-02_wp, &
      & 4.01279440930580E-18_wp, 6.09607146462613E-17_wp, 7.46102863717653E-18_wp, &
      & 1.18968197226907E-01_wp,-9.65885020789207E-18_wp, 1.39766376211368E-17_wp, &
      & 2.06454228447258E-17_wp, 1.59613420193925E-16_wp,-3.64893594596737E-18_wp, &
      & 1.45926928037665E-01_wp,-5.48413077769616E-17_wp, 3.51920727141113E-17_wp, &
      &-3.09892298541188E-17_wp,-3.86489796446311E-20_wp, 6.34981035886410E-17_wp, &
      & 4.87261607756941E-02_wp,-6.76508056190893E-18_wp, 3.56878079377531E-17_wp, &
      &-5.85075287079675E-18_wp,-1.18968197226926E-01_wp,-3.51744700557441E-17_wp, &
      &-4.85722573273506E-17_wp, 1.16139744113067E-17_wp, 2.23819063966983E-16_wp, &
      &-2.42095363153376E-17_wp, 1.45926928037629E-01_wp, 6.08318062634637E-17_wp, &
      & 2.92600818645547E-17_wp,-5.75746099072571E-17_wp,-3.87667277030112E-20_wp, &
      &-6.78407215531017E-03_wp,-6.49358823473300E-18_wp,-7.23612469803392E-05_wp, &
      & 1.71492294696653E-18_wp, 2.02291975686579E-18_wp,-1.49641218981684E-17_wp, &
      &-5.50832892445549E-05_wp, 4.10577711599016E-18_wp, 1.31819935414578E-16_wp, &
      &-8.00924381689913E-21_wp, 2.41084201803115E-18_wp,-1.87298597344789E-17_wp, &
      & 1.10418504944028E-04_wp, 5.05023321416519E-18_wp, 1.49847858431830E-16_wp, &
      &-1.15659771140108E-31_wp,-6.78407215531131E-03_wp,-6.76508056190893E-18_wp, &
      & 7.23612469821558E-05_wp, 1.78779515786626E-18_wp, 2.10564164289910E-18_wp, &
      & 1.46602671094109E-17_wp,-5.50832892447886E-05_wp,-4.18296130758492E-18_wp, &
      & 1.31757672731487E-16_wp, 1.57079535276860E-19_wp,-2.50761191210483E-18_wp, &
      &-1.80754125702893E-17_wp,-1.10418504943955E-04_wp, 5.11437644470088E-18_wp, &
      &-1.50025515835466E-16_wp,-9.84380371715469E-33_wp,-8.67695391478395E-18_wp, &
      & 3.61801788060898E-17_wp,-1.79954688428111E-18_wp, 4.87261607756798E-02_wp, &
      &-1.88534072005163E-18_wp, 5.11666666624951E-17_wp, 8.17918707691012E-18_wp, &
      & 1.18968197226908E-01_wp, 9.70447558908168E-17_wp, 5.28413089972434E-19_wp, &
      &-9.78131575785710E-19_wp, 4.05535275106961E-17_wp, 8.40382497671730E-18_wp, &
      & 1.45926928037666E-01_wp, 1.26925784809929E-16_wp,-1.04755344235895E-16_wp, &
      &-5.50815884495145E-18_wp, 3.56878079377531E-17_wp, 1.78779515786626E-18_wp, &
      & 4.87261607756949E-02_wp,-1.32611590991196E-18_wp,-3.70173868303146E-17_wp, &
      & 1.44573403802346E-18_wp,-1.18968197226927E-01_wp, 1.32774862544140E-16_wp, &
      & 7.03012249162352E-18_wp, 4.60847936837848E-18_wp, 5.72686452748858E-17_wp, &
      &-2.55939443012652E-17_wp, 1.45926928037629E-01_wp,-1.18122598326716E-16_wp, &
      &-1.15024453040342E-16_wp, 1.89566830830589E-18_wp,-5.80448036224271E-18_wp, &
      &-1.24107819697753E-18_wp,-2.78480522972095E-18_wp, 9.85849302395307E-05_wp, &
      &-1.45833735947789E-17_wp, 5.84024309787051E-18_wp,-6.64179475456893E-18_wp, &
      & 6.10092008616677E-03_wp,-1.53474430223249E-19_wp, 1.12134067473778E-04_wp, &
      &-1.72526012441359E-17_wp, 9.75760966247523E-18_wp,-1.19469628797602E-17_wp, &
      & 6.93940730020448E-03_wp,-8.23570112181639E-19_wp, 4.04679335389210E-18_wp, &
      &-5.85075287079675E-18_wp, 2.10564164289910E-18_wp,-1.32611590991196E-18_wp, &
      & 9.85849302395349E-05_wp, 1.51533102688874E-17_wp, 4.42295576844490E-18_wp, &
      & 1.12006403354152E-17_wp, 6.10092008616703E-03_wp,-2.89719626342384E-19_wp, &
      &-1.12134067473771E-04_wp,-1.85369840126384E-17_wp,-1.00918632067683E-17_wp, &
      &-6.80448596299337E-18_wp,-6.93940730020406E-03_wp,-1.63065940971428E-19_wp, &
      & 6.32263657180051E-17_wp,-1.18968197226890E-01_wp,-7.94037544109613E-18_wp, &
      &-9.92670042757072E-17_wp,-1.73482991883529E-17_wp,-2.90468851354336E-01_wp, &
      & 2.49964298671853E-17_wp, 8.58831364103808E-17_wp, 3.32772011249550E-18_wp, &
      &-3.89706895667578E-16_wp, 9.89676259380455E-18_wp,-3.56290404972032E-01_wp, &
      & 1.31064878600120E-16_wp, 5.94289669462712E-17_wp, 1.36782316927511E-16_wp, &
      & 9.43640820409146E-20_wp, 1.90764923804851E-17_wp,-1.18968197226926E-01_wp, &
      & 1.46602671094109E-17_wp,-3.70173868303146E-17_wp, 1.51533102688874E-17_wp, &
      & 2.90468851354381E-01_wp, 8.72945310535409E-17_wp,-1.33490843362475E-19_wp, &
      & 2.53786086160670E-17_wp,-5.46469291265246E-16_wp, 5.81215702619420E-17_wp, &
      &-3.56290404971943E-01_wp,-1.45691083375633E-16_wp, 7.44263965640176E-17_wp, &
      & 7.94522516938117E-17_wp, 9.46515718412435E-20_wp, 5.16421460896905E-03_wp, &
      &-3.53811371739650E-17_wp, 5.50832892434057E-05_wp, 1.50120628647457E-18_wp, &
      & 4.48592584932744E-18_wp,-8.70632287761936E-17_wp, 4.19308522248157E-05_wp, &
      & 3.72720093683120E-18_wp, 2.72562904660025E-16_wp, 6.09684758330446E-21_wp, &
      & 5.01879572423427E-18_wp,-1.06506841820406E-16_wp,-8.40534775099834E-05_wp, &
      & 4.56109222030375E-18_wp, 3.10090590031171E-16_wp, 3.17020579613494E-32_wp, &
      & 5.16421460896991E-03_wp,-3.51744700557441E-17_wp,-5.50832892447886E-05_wp, &
      & 1.44573403802346E-18_wp, 4.42295576844490E-18_wp, 8.72945310535409E-17_wp, &
      & 4.19308522249936E-05_wp,-3.66844629497529E-18_wp, 2.72610300656170E-16_wp, &
      &-1.19573084170877E-19_wp,-4.94513192920884E-18_wp,-1.07005024258406E-16_wp, &
      & 8.40534775099282E-05_wp, 4.51226470187122E-18_wp,-3.09955352530976E-16_wp, &
      &-9.06713103001811E-33_wp, 4.11665386035113E-18_wp,-4.94396190653390E-17_wp, &
      & 4.21165393765978E-18_wp,-1.18968197226890E-01_wp, 1.25660253041414E-17_wp, &
      &-3.12250225675825E-17_wp,-2.01086245564250E-17_wp,-2.90468851354336E-01_wp, &
      & 2.55838335693153E-16_wp,-1.29015608257161E-18_wp, 1.14453954723415E-17_wp, &
      & 1.73472347597681E-17_wp,-2.02406908958029E-17_wp,-3.56290404972031E-01_wp, &
      & 2.50607549961322E-16_wp, 2.55767215295300E-16_wp,-3.62017232989779E-18_wp, &
      &-4.85722573273506E-17_wp,-4.18296130758492E-18_wp,-1.18968197226927E-01_wp, &
      & 1.12006403354152E-17_wp,-1.33490843362475E-19_wp,-3.66844629497529E-18_wp, &
      & 2.90468851354381E-01_wp, 1.68600878286269E-16_wp,-1.71645166743796E-17_wp, &
      &-2.03091336638168E-17_wp,-2.49086271625784E-17_wp, 6.22115228595914E-17_wp, &
      &-3.56290404971942E-01_wp,-2.72101121568159E-16_wp, 2.80839934798405E-16_wp, &
      &-1.86487216170029E-17_wp, 1.44775446971247E-17_wp,-7.82542368597328E-17_wp, &
      & 4.25039971179187E-17_wp, 6.10092008616677E-03_wp, 9.89187435175968E-18_wp, &
      & 3.60319007402261E-16_wp, 1.13522895665109E-16_wp, 3.77554924544321E-01_wp, &
      &-9.49775215935025E-18_wp, 6.93940730020449E-03_wp, 5.14569913135540E-17_wp, &
      & 6.06061786705714E-16_wp,-9.59216313667728E-17_wp, 4.29444635007041E-01_wp, &
      &-5.09665668735330E-17_wp, 1.14473470249282E-16_wp, 1.16139744113067E-17_wp, &
      & 1.31757672731487E-16_wp, 1.32774862544140E-16_wp, 6.10092008616703E-03_wp, &
      & 2.53786086160670E-17_wp, 2.72610300656170E-16_wp, 1.68600878286269E-16_wp, &
      & 3.77554924544337E-01_wp,-1.79292746205145E-17_wp,-6.93940730020409E-03_wp, &
      &-2.80269278987693E-17_wp,-6.26747039079427E-16_wp, 2.22320118542743E-16_wp, &
      &-4.29444635007015E-01_wp,-1.00913219923688E-17_wp,-1.47266519840908E-17_wp, &
      & 2.23819063966914E-16_wp,-1.57079535272932E-19_wp, 7.03012249162146E-18_wp, &
      &-2.89719626342310E-19_wp, 5.46469291265161E-16_wp,-1.19573084170256E-19_wp, &
      & 1.71645166743767E-17_wp,-1.79292746205137E-17_wp, 7.33679203094675E-31_wp, &
      &-3.29537588045326E-19_wp, 6.70301700791082E-16_wp, 2.39693042421374E-19_wp, &
      & 2.10540736762259E-17_wp,-2.03934058194061E-17_wp,-1.28711492963744E-32_wp, &
      &-1.47266519840931E-17_wp, 2.23819063966983E-16_wp, 1.57079535276860E-19_wp, &
      & 7.03012249162352E-18_wp,-2.89719626342384E-19_wp,-5.46469291265246E-16_wp, &
      &-1.19573084170877E-19_wp,-1.71645166743796E-17_wp,-1.79292746205145E-17_wp, &
      & 1.03029864438581E-30_wp, 3.29537588045180E-19_wp, 6.70301700790914E-16_wp, &
      &-2.39693042421191E-19_wp, 2.10540736762206E-17_wp, 2.03934058194045E-17_wp, &
      &-1.62943764584851E-32_wp, 8.39838980502871E-18_wp,-2.42621683411740E-17_wp, &
      & 1.52422615989048E-18_wp, 6.26764539114895E-18_wp,-1.12134067473766E-04_wp, &
      &-5.87698368612470E-17_wp,-6.55720582976770E-18_wp, 1.51237377141199E-17_wp, &
      &-6.93940730020379E-03_wp, 1.74567371223243E-19_wp,-1.27545346511494E-04_wp, &
      &-7.28101091714507E-17_wp,-1.12704463656210E-17_wp, 2.28732247113919E-17_wp, &
      &-7.89313300256366E-03_wp, 9.36758450854204E-19_wp, 5.95162232829546E-18_wp, &
      &-2.42095363153376E-17_wp,-2.50761191210483E-18_wp, 4.60847936837848E-18_wp, &
      &-1.12134067473771E-04_wp, 5.81215702619420E-17_wp,-4.94513192920884E-18_wp, &
      &-2.03091336638168E-17_wp,-6.93940730020409E-03_wp, 3.29537588045180E-19_wp, &
      & 1.27545346511487E-04_wp,-7.13492057283204E-17_wp, 1.16506384436397E-17_wp, &
      & 1.70239853584180E-17_wp, 7.89313300256319E-03_wp, 1.85477102668076E-19_wp, &
      &-6.88288516982198E-17_wp, 1.45926928037584E-01_wp, 9.83276309958546E-18_wp, &
      & 1.33227792809983E-16_wp, 2.13296369409181E-17_wp, 3.56290404971887E-01_wp, &
      &-3.05898919742584E-17_wp,-7.82644961576860E-17_wp,-9.79208467595434E-19_wp, &
      & 4.78016238334485E-16_wp,-1.20823885568917E-17_wp, 4.37027419921967E-01_wp, &
      &-1.60906779637513E-16_wp,-4.11883049958516E-17_wp,-1.64248811785728E-16_wp, &
      &-1.15747409226165E-19_wp,-1.46744176598035E-17_wp, 1.45926928037629E-01_wp, &
      &-1.80754125702893E-17_wp, 5.72686452748858E-17_wp,-1.85369840126384E-17_wp, &
      &-3.56290404971943E-01_wp,-1.07005024258406E-16_wp,-2.49086271625784E-17_wp, &
      &-2.80269278987693E-17_wp, 6.70301700790914E-16_wp,-7.13492057283204E-17_wp, &
      & 4.37027419921857E-01_wp, 1.78847348474315E-16_wp,-5.85255905388085E-17_wp, &
      &-1.00985486707195E-16_wp,-1.16100045513699E-19_wp, 1.03520480376186E-02_wp, &
      & 6.04175268235442E-17_wp, 1.10418504941183E-04_wp,-2.54827460989731E-17_wp, &
      &-9.96563504787242E-18_wp, 1.46154745803971E-16_wp, 8.40534775095716E-05_wp, &
      &-6.20937448582820E-17_wp,-6.26842047839864E-16_wp, 1.22215794349762E-20_wp, &
      &-1.15029739481702E-17_wp, 1.79845991815207E-16_wp,-1.68491378224491E-04_wp, &
      &-7.61859354340708E-17_wp,-7.12857165023500E-16_wp, 2.83072876256128E-31_wp, &
      & 1.03520480376204E-02_wp, 6.08318062634637E-17_wp,-1.10418504943955E-04_wp, &
      &-2.55939443012652E-17_wp,-1.00918632067683E-17_wp,-1.45691083375633E-16_wp, &
      & 8.40534775099282E-05_wp, 6.22115228595914E-17_wp,-6.26747039079427E-16_wp, &
      &-2.39693042421191E-19_wp, 1.16506384436397E-17_wp, 1.78847348474315E-16_wp, &
      & 1.68491378224380E-04_wp,-7.62838137862578E-17_wp, 7.13128258539958E-16_wp, &
      & 8.03366913497250E-32_wp,-3.50607319354894E-18_wp, 3.12802712863402E-17_wp, &
      &-5.14957095434112E-18_wp, 1.45926928037584E-01_wp,-8.47927333238547E-18_wp, &
      &-3.51816401787238E-17_wp, 2.46778616093225E-17_wp, 3.56290404971886E-01_wp, &
      & 1.15314259280036E-16_wp, 1.58251127786369E-18_wp,-6.15168405973870E-18_wp, &
      &-1.09030624000210E-16_wp, 2.48022017384882E-17_wp, 4.37027419921966E-01_wp, &
      & 1.80707905874698E-16_wp,-3.13725221452161E-16_wp, 5.98395275473543E-18_wp, &
      & 2.92600818645547E-17_wp, 5.11437644470088E-18_wp, 1.45926928037629E-01_wp, &
      &-6.80448596299337E-18_wp, 7.44263965640176E-17_wp, 4.51226470187122E-18_wp, &
      &-3.56290404971942E-01_wp, 2.22320118542743E-16_wp, 2.10540736762206E-17_wp, &
      & 1.70239853584180E-17_wp,-5.85255905388085E-17_wp,-7.62838137862578E-17_wp, &
      & 4.37027419921856E-01_wp,-1.54343793866147E-16_wp,-3.44479532435444E-16_wp, &
      & 3.61719196879460E-17_wp,-6.08317383838331E-17_wp, 8.91687715689248E-17_wp, &
      &-1.54452490602129E-17_wp,-6.93940730020377E-03_wp,-1.19570175471162E-16_wp, &
      &-4.09718409789389E-16_wp,-4.87966690713799E-17_wp,-4.29444635006996E-01_wp, &
      & 1.08030870326530E-17_wp,-7.89313300256364E-03_wp,-1.91393372967111E-16_wp, &
      &-6.89600103884371E-16_wp, 2.07635889748987E-16_wp,-4.88465869592074E-01_wp, &
      & 5.79712176578819E-17_wp,-1.15246079074340E-16_wp,-5.75746099072571E-17_wp, &
      &-1.50025515835466E-16_wp,-1.18122598326716E-16_wp,-6.93940730020406E-03_wp, &
      & 7.94522516938117E-17_wp,-3.09955352530976E-16_wp,-2.72101121568159E-16_wp, &
      &-4.29444635007015E-01_wp, 2.03934058194045E-17_wp, 7.89313300256319E-03_wp, &
      &-1.00985486707195E-16_wp, 7.13128258539958E-16_wp,-1.54343793866147E-16_wp, &
      & 4.88465869592044E-01_wp, 1.14782348422054E-17_wp, 2.18504335499174E-31_wp, &
      &-3.87667277030004E-20_wp, 8.44369028099050E-33_wp,-1.15024453040306E-16_wp, &
      &-1.63065940971419E-19_wp,-9.46515718412620E-20_wp,-2.73269110400166E-32_wp, &
      &-2.80839934798361E-16_wp,-1.00913219923683E-17_wp,-1.12052026646400E-33_wp, &
      &-1.85477102668095E-19_wp,-1.16100045513692E-19_wp,-3.92090630566152E-32_wp, &
      &-3.44479532435530E-16_wp,-1.14782348422061E-17_wp, 2.48650912671926E-31_wp, &
      & 2.07451493237937E-31_wp,-3.87667277030112E-20_wp,-9.84380371715469E-33_wp, &
      &-1.15024453040342E-16_wp,-1.63065940971428E-19_wp, 9.46515718412435E-20_wp, &
      &-9.06713103001811E-33_wp, 2.80839934798405E-16_wp,-1.00913219923688E-17_wp, &
      &-1.62943764584851E-32_wp, 1.85477102668076E-19_wp,-1.16100045513699E-19_wp, &
      & 8.03366913497250E-32_wp,-3.44479532435444E-16_wp, 1.14782348422054E-17_wp, &
      & 2.71799965258949E-31_wp], shape(density))

   call get_structure(mol, "f-block", "Ce2")
   call test_energy_generic(error, mol, density, qsh, make_exchange_gxtb, &
      & -0.490545588540915_wp, thr_in=thr1)

end subroutine test_e_fock_ce2


subroutine test_p_fock_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(2) = [&
      & -6.66133814775094E-16_wp,  4.44089209850063E-16_wp]

   real(wp), parameter :: density(2, 2, 1) = reshape([&
      & 5.93683766916992E-1_wp, 5.93683766916992E-1_wp, 5.93683766916992E-1_wp, &
      & 5.93683766916992E-1_wp], shape(density))

   call get_structure(mol, "MB16-43", "H2")
   call test_numpot(error, mol, density, qsh, make_exchange_gxtb, &
      & thr_in=thr1)

end subroutine test_p_fock_h2


subroutine test_p_fock_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(3) = [&
      &  1.88324089567125E-1_wp,  2.01980267353847E-1_wp, -3.90304356969502E-1_wp]

   real(wp), parameter :: density(5, 5, 1) = reshape([&
      & 7.43138968868805E-02_wp, 6.30732585479418E-45_wp, 1.15038033099931E-01_wp, &
      & 0.00000000000000E+00_wp, 2.77067464359636E-01_wp, 6.30732585479418E-45_wp, &
      & 5.35328668056679E-88_wp, 9.76375066914250E-45_wp, 0.00000000000000E+00_wp, &
      & 2.35158544321515E-44_wp, 1.15038033099931E-01_wp, 9.76375066914250E-45_wp, &
      & 1.78079062111965E-01_wp, 0.00000000000000E+00_wp, 4.28900884910328E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 2.77067464359636E-01_wp, &
      & 2.35158544321515E-44_wp, 4.28900884910328E-01_wp, 0.00000000000000E+00_wp, &
      & 1.03300167293785E+00_wp], shape(density))

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpot(error, mol, density, qsh, make_exchange_gxtb, &
      & thr_in=thr1)

end subroutine test_p_fock_lih


subroutine test_p_fock_no(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(4) = [&
      & -3.72959864890934E-1_wp,  4.88258685456765E-1_wp, -2.06922792771332E-1_wp, &
      &  9.16239721052889E-2_wp]

   real(wp), parameter :: density(8, 8, 2) = reshape([&
      & 9.42009046415958E-01_wp,-1.97869914210805E-16_wp,-3.13530611389513E-01_wp, &
      & 1.57008624511250E-16_wp,-1.91210322106419E-01_wp, 1.16548396645701E-16_wp, &
      &-2.58330200259595E-02_wp,-1.75070805017504E-16_wp,-1.97869914210805E-16_wp, &
      & 7.16403612559420E-01_wp, 6.75092136044200E-17_wp, 4.11552712525306E-01_wp, &
      & 1.21088275692219E-16_wp,-3.76653483713363E-03_wp,-1.64250151646619E-17_wp, &
      &-3.38085910709977E-01_wp,-3.13530611389513E-01_wp, 6.75092136044200E-17_wp, &
      & 4.49934652749990E-01_wp,-8.24568998262442E-17_wp,-1.71924645027865E-03_wp, &
      &-1.23858164254808E-17_wp,-4.07076133805151E-01_wp, 1.00925254078120E-16_wp, &
      & 1.57008624511250E-16_wp, 4.11552712525306E-01_wp,-8.24568998262442E-17_wp, &
      & 6.31376987893123E-01_wp,-1.61242008793196E-16_wp,-3.38085910709977E-01_wp, &
      &-6.74819946100851E-18_wp, 6.60818781608317E-02_wp,-1.91210322106419E-01_wp, &
      & 1.21088275692219E-16_wp,-1.71924645027865E-03_wp,-1.61242008793196E-16_wp, &
      & 1.01808686826266E+00_wp, 4.52238546687169E-17_wp, 2.62723878284716E-01_wp, &
      &-1.97376402183765E-17_wp, 1.16548396645701E-16_wp,-3.76653483713363E-03_wp, &
      &-1.23858164254808E-17_wp,-3.38085910709977E-01_wp, 4.52238546687169E-17_wp, &
      & 8.37111178906109E-01_wp, 3.74456682512930E-17_wp, 2.77733761780434E-01_wp, &
      &-2.58330200259595E-02_wp,-1.64250151646619E-17_wp,-4.07076133805151E-01_wp, &
      &-6.74819946100851E-18_wp, 2.62723878284716E-01_wp, 3.74456682512930E-17_wp, &
      & 5.33778495172291E-01_wp, 2.29800571514756E-17_wp,-1.75070805017504E-16_wp, &
      &-3.38085910709977E-01_wp, 1.00925254078120E-16_wp, 6.60818781608317E-02_wp, &
      &-1.97376402183765E-17_wp, 2.77733761780434E-01_wp, 2.29800571514756E-17_wp, &
      & 7.79731495571835E-01_wp, 9.54305653170375E-01_wp,-1.73912460120364E-16_wp, &
      &-3.11307650277130E-01_wp, 1.94985952807570E-16_wp,-2.03773626333088E-01_wp, &
      & 1.87683417049811E-17_wp,-1.70067070427977E-02_wp, 2.26960265829082E-17_wp, &
      &-1.73912460120364E-16_wp, 2.40284367737361E-01_wp, 1.07958284666329E-16_wp, &
      &-1.87742209952433E-05_wp, 1.14092931051935E-16_wp, 3.64460563073953E-01_wp, &
      &-4.86861943671765E-17_wp,-6.78042890246151E-06_wp,-3.11307650277130E-01_wp, &
      & 1.07958284666329E-16_wp, 4.38731609559674E-01_wp,-1.59063207389752E-16_wp, &
      &-4.10640194692124E-03_wp,-3.43859837540706E-17_wp,-4.09548671902707E-01_wp, &
      &-2.10949694797432E-17_wp, 1.94985952807570E-16_wp,-1.87742209952433E-05_wp, &
      &-1.59063207389752E-16_wp, 2.40288246483859E-01_wp,-6.87541757092717E-17_wp, &
      &-6.78042890286945E-06_wp, 6.84810692055402E-17_wp, 3.64461963907834E-01_wp, &
      &-2.03773626333088E-01_wp, 1.14092931051935E-16_wp,-4.10640194692124E-03_wp, &
      &-6.87541757092717E-17_wp, 1.03092149514981E+00_wp,-1.10878908539341E-16_wp, &
      & 2.53665476333654E-01_wp, 5.86026180268729E-17_wp, 1.87683417049811E-17_wp, &
      & 3.64460563073953E-01_wp,-3.43859837540706E-17_wp,-6.78042890286945E-06_wp, &
      &-1.10878908539341E-16_wp, 5.52809588728990E-01_wp, 1.03373764663983E-17_wp, &
      & 2.26235036331901E-05_wp,-1.70067070427977E-02_wp,-4.86861943671765E-17_wp, &
      &-4.09548671902707E-01_wp, 6.84810692055402E-17_wp, 2.53665476333654E-01_wp, &
      & 1.03373764663983E-17_wp, 5.38687781288217E-01_wp,-3.05073307416719E-17_wp, &
      & 2.26960265829082E-17_wp,-6.78042890246151E-06_wp,-2.10949694797432E-17_wp, &
      & 3.64461963907834E-01_wp, 5.86026180268729E-17_wp, 2.26235036331901E-05_wp, &
      &-3.05073307416719E-17_wp, 5.52804914722249E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "NO")
   call test_numpot(error, mol, density, qsh, make_exchange_gxtb, &
      & thr_in=thr1*10)

end subroutine test_p_fock_no


subroutine test_p_fock_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(6) = [&
      & -1.65829465017236E-1_wp,  8.50038195478131E-2_wp,  8.08256454270758E-2_wp, &
      & -1.65829465017241E-1_wp,  8.50038195477536E-2_wp,  8.08256454270730E-2_wp]

   real(wp), parameter :: density(18, 18, 1) = reshape([&
      & 2.01512274410727E+00_wp,-2.98954363036077E-16_wp,-3.99929219814917E-01_wp, &
      &-3.60921385847533E-16_wp, 1.66928202631975E-04_wp, 2.96505862833536E-17_wp, &
      &-3.51352512626524E-02_wp,-2.26041467075836E-17_wp,-1.25121106619603E-05_wp, &
      &-2.61461503009251E-01_wp, 2.00709987824803E-16_wp, 8.46120257485226E-02_wp, &
      & 4.15182726310734E-16_wp,-2.96659169389449E-05_wp, 3.71826627538686E-17_wp, &
      & 1.18677041085955E-02_wp, 5.96815828699915E-17_wp, 2.22361008940176E-06_wp, &
      &-2.98954363036077E-16_wp, 1.45486170073001E+00_wp,-9.00635884019608E-16_wp, &
      & 5.97801199387253E-01_wp, 1.67776446953954E-17_wp,-2.61299185132976E-03_wp, &
      & 1.68298746105631E-16_wp,-4.41320346165751E-02_wp,-2.78389072323402E-18_wp, &
      & 4.34685388862911E-16_wp, 1.61265924256514E-01_wp,-6.14873875231672E-16_wp, &
      &-6.02461837620382E-01_wp, 4.43319518778254E-17_wp,-1.03162491097761E-01_wp, &
      &-2.04398675254949E-16_wp,-5.40117632485737E-02_wp, 1.09406633163272E-16_wp, &
      &-3.99929219814917E-01_wp,-9.00635884019608E-16_wp, 8.16076844347620E-01_wp, &
      &-1.55848848678421E-16_wp,-2.73809977947621E-05_wp, 4.89930425005989E-17_wp, &
      & 5.39351481156965E-02_wp, 2.47586597786769E-16_wp, 2.05234387570902E-06_wp, &
      &-8.46120257485189E-02_wp, 6.66084959945643E-16_wp,-7.72403989719108E-01_wp, &
      & 5.71110633457431E-16_wp,-1.51824865425160E-07_wp, 6.80560093611612E-17_wp, &
      & 4.74250266390748E-02_wp,-1.14353946762662E-17_wp, 1.13800395663470E-08_wp, &
      &-3.60921385847533E-16_wp, 5.97801199387253E-01_wp,-1.55848848678421E-16_wp, &
      & 1.36524527111739E+00_wp, 1.53974719339088E-17_wp,-4.41320346165748E-02_wp, &
      & 1.44159180113409E-16_wp, 4.00284528304191E-03_wp, 4.02705177739493E-18_wp, &
      & 8.69046775550039E-17_wp,-6.02461837620375E-01_wp,-7.97303127069371E-16_wp, &
      & 2.51581030550614E-01_wp,-5.41615397023848E-17_wp,-5.40117632485737E-02_wp, &
      & 4.07523509225128E-18_wp,-9.50655829878735E-02_wp, 9.48838887830889E-17_wp, &
      & 1.66928202631975E-04_wp, 1.67776446953954E-17_wp,-2.73809977947621E-05_wp, &
      & 1.53974719339088E-17_wp, 1.38974440043868E-08_wp,-6.94704531194593E-19_wp, &
      &-2.60092206439580E-06_wp,-6.56220312997718E-19_wp,-1.04168351758665E-09_wp, &
      &-2.96659169388496E-05_wp,-7.73908741086395E-18_wp, 1.51824865387540E-07_wp, &
      &-7.33055884782161E-18_wp,-3.07941415440952E-09_wp,-1.30780677462207E-18_wp, &
      & 1.45801168257227E-06_wp,-1.20266021436374E-18_wp, 2.30817621363979E-10_wp, &
      & 2.96505862833536E-17_wp,-2.61299185132976E-03_wp, 4.89930425005989E-17_wp, &
      &-4.41320346165748E-02_wp,-6.94704531194593E-19_wp, 7.48444822726994E-03_wp, &
      &-1.84444900024127E-17_wp, 4.63153745136131E-03_wp, 2.56316174731192E-19_wp, &
      &-4.56723208050267E-17_wp, 1.03162491097758E-01_wp, 4.43388496380181E-17_wp, &
      & 5.40117632485682E-02_wp, 2.70455166704829E-18_wp, 1.16466225731597E-03_wp, &
      & 1.89516760801754E-17_wp, 3.39354102094470E-03_wp,-1.17177365770997E-17_wp, &
      &-3.51352512626524E-02_wp, 1.68298746105631E-16_wp, 5.39351481156965E-02_wp, &
      & 1.44159180113409E-16_wp,-2.60092206439580E-06_wp,-1.84444900024127E-17_wp, &
      & 3.73731424127323E-03_wp,-7.83051079850619E-18_wp, 1.94952226013385E-07_wp, &
      & 1.18677041085919E-02_wp,-2.17279393388892E-16_wp,-4.74250266390769E-02_wp, &
      &-2.59398305030398E-16_wp, 1.45801168257080E-06_wp,-1.41001874145128E-17_wp, &
      & 2.76687889875095E-03_wp,-1.47129262760623E-17_wp,-1.09285328835862E-07_wp, &
      &-2.26041467075836E-17_wp,-4.41320346165751E-02_wp, 2.47586597786769E-16_wp, &
      & 4.00284528304191E-03_wp,-6.56220312997718E-19_wp, 4.63153745136131E-03_wp, &
      &-7.83051079850619E-18_wp, 6.79013404651283E-03_wp, 6.04554882466166E-19_wp, &
      &-1.74543578973199E-16_wp, 5.40117632485687E-02_wp,-1.45442246694960E-16_wp, &
      & 9.50655829878701E-02_wp,-2.78339838142217E-18_wp, 3.39354102094469E-03_wp, &
      & 3.97034075802654E-17_wp, 6.55936228742013E-04_wp,-1.12393067366785E-17_wp, &
      &-1.25121106619603E-05_wp,-2.78389072323402E-18_wp, 2.05234387570902E-06_wp, &
      & 4.02705177739493E-18_wp,-1.04168351758665E-09_wp, 2.56316174731192E-19_wp, &
      & 1.94952226013385E-07_wp, 6.04554882466166E-19_wp, 7.80794332014712E-11_wp, &
      & 2.22361008999868E-06_wp, 2.61284439827709E-18_wp,-1.13800387946783E-08_wp, &
      & 9.03292745688553E-18_wp, 2.30817621408582E-10_wp, 1.84803242584056E-19_wp, &
      &-1.09285328887014E-07_wp,-1.95198444371785E-19_wp,-1.73009448131981E-11_wp, &
      &-2.61461503009251E-01_wp, 4.34685388862911E-16_wp,-8.46120257485189E-02_wp, &
      & 8.69046775550039E-17_wp,-2.96659169388496E-05_wp,-4.56723208050267E-17_wp, &
      & 1.18677041085919E-02_wp,-1.74543578973199E-16_wp, 2.22361008999868E-06_wp, &
      & 2.01512274410727E+00_wp,-2.95337324834156E-16_wp, 3.99929219814911E-01_wp, &
      &-2.60805282861184E-16_wp, 1.66928202632090E-04_wp,-3.91687455452402E-17_wp, &
      &-3.51352512626490E-02_wp, 3.43718096863094E-18_wp,-1.25121106621160E-05_wp, &
      & 2.00709987824803E-16_wp, 1.61265924256514E-01_wp, 6.66084959945643E-16_wp, &
      &-6.02461837620375E-01_wp,-7.73908741086395E-18_wp, 1.03162491097758E-01_wp, &
      &-2.17279393388892E-16_wp, 5.40117632485687E-02_wp, 2.61284439827709E-18_wp, &
      &-2.95337324834156E-16_wp, 1.45486170073003E+00_wp, 4.61646074650480E-16_wp, &
      & 5.97801199387235E-01_wp, 4.81747951661595E-17_wp, 2.61299185133292E-03_wp, &
      & 2.21800860913915E-16_wp, 4.41320346165780E-02_wp,-1.45847563550088E-16_wp, &
      & 8.46120257485226E-02_wp,-6.14873875231672E-16_wp,-7.72403989719108E-01_wp, &
      &-7.97303127069371E-16_wp, 1.51824865387540E-07_wp, 4.43388496380181E-17_wp, &
      &-4.74250266390769E-02_wp,-1.45442246694960E-16_wp,-1.13800387946783E-08_wp, &
      & 3.99929219814911E-01_wp, 4.61646074650480E-16_wp, 8.16076844347641E-01_wp, &
      & 9.30644287891831E-16_wp, 2.73809977948286E-05_wp, 4.76651490748586E-17_wp, &
      &-5.39351481156958E-02_wp, 9.39184312631712E-17_wp,-2.05234387641953E-06_wp, &
      & 4.15182726310734E-16_wp,-6.02461837620382E-01_wp, 5.71110633457431E-16_wp, &
      & 2.51581030550614E-01_wp,-7.33055884782161E-18_wp, 5.40117632485682E-02_wp, &
      &-2.59398305030398E-16_wp, 9.50655829878701E-02_wp, 9.03292745688553E-18_wp, &
      &-2.60805282861184E-16_wp, 5.97801199387235E-01_wp, 9.30644287891831E-16_wp, &
      & 1.36524527111740E+00_wp,-5.05746530671741E-17_wp, 4.41320346165780E-02_wp, &
      & 3.53371679850913E-16_wp,-4.00284528303917E-03_wp,-1.41955353109642E-16_wp, &
      &-2.96659169389449E-05_wp, 4.43319518778254E-17_wp,-1.51824865425160E-07_wp, &
      &-5.41615397023848E-17_wp,-3.07941415440952E-09_wp, 2.70455166704829E-18_wp, &
      & 1.45801168257080E-06_wp,-2.78339838142217E-18_wp, 2.30817621408582E-10_wp, &
      & 1.66928202632090E-04_wp, 4.81747951661595E-17_wp, 2.73809977948286E-05_wp, &
      &-5.05746530671741E-17_wp, 1.38974440044059E-08_wp,-2.39255576630622E-18_wp, &
      &-2.60092206440109E-06_wp, 3.06620883893745E-18_wp,-1.04168351759322E-09_wp, &
      & 3.71826627538686E-17_wp,-1.03162491097761E-01_wp, 6.80560093611612E-17_wp, &
      &-5.40117632485737E-02_wp,-1.30780677462207E-18_wp, 1.16466225731597E-03_wp, &
      &-1.41001874145128E-17_wp, 3.39354102094469E-03_wp, 1.84803242584056E-19_wp, &
      &-3.91687455452402E-17_wp, 2.61299185133292E-03_wp, 4.76651490748586E-17_wp, &
      & 4.41320346165780E-02_wp,-2.39255576630622E-18_wp, 7.48444822727060E-03_wp, &
      & 1.59370810633493E-17_wp, 4.63153745136205E-03_wp,-9.15386829425696E-18_wp, &
      & 1.18677041085955E-02_wp,-2.04398675254949E-16_wp, 4.74250266390748E-02_wp, &
      & 4.07523509225128E-18_wp, 1.45801168257227E-06_wp, 1.89516760801754E-17_wp, &
      & 2.76687889875095E-03_wp, 3.97034075802654E-17_wp,-1.09285328887014E-07_wp, &
      &-3.51352512626490E-02_wp, 2.21800860913915E-16_wp,-5.39351481156958E-02_wp, &
      & 3.53371679850913E-16_wp,-2.60092206440109E-06_wp, 1.59370810633493E-17_wp, &
      & 3.73731424127310E-03_wp, 1.71392595868698E-18_wp, 1.94952226055385E-07_wp, &
      & 5.96815828699915E-17_wp,-5.40117632485737E-02_wp,-1.14353946762662E-17_wp, &
      &-9.50655829878735E-02_wp,-1.20266021436374E-18_wp, 3.39354102094470E-03_wp, &
      &-1.47129262760623E-17_wp, 6.55936228742013E-04_wp,-1.95198444371785E-19_wp, &
      & 3.43718096863094E-18_wp, 4.41320346165780E-02_wp, 9.39184312631712E-17_wp, &
      &-4.00284528303917E-03_wp, 3.06620883893745E-18_wp, 4.63153745136205E-03_wp, &
      & 1.71392595868698E-18_wp, 6.79013404651345E-03_wp,-8.12591467858445E-18_wp, &
      & 2.22361008940176E-06_wp, 1.09406633163272E-16_wp, 1.13800395663470E-08_wp, &
      & 9.48838887830889E-17_wp, 2.30817621363979E-10_wp,-1.17177365770997E-17_wp, &
      &-1.09285328835862E-07_wp,-1.12393067366785E-17_wp,-1.73009448131981E-11_wp, &
      &-1.25121106621160E-05_wp,-1.45847563550088E-16_wp,-2.05234387641953E-06_wp, &
      &-1.41955353109642E-16_wp,-1.04168351759322E-09_wp,-9.15386829425696E-18_wp, &
      & 1.94952226055385E-07_wp,-8.12591467858445E-18_wp, 7.80794332023474E-11_wp],&
      & shape(density))

   call get_structure(mol, "MB16-43", "S2")
   call test_numpot(error, mol, density, qsh, make_exchange_gxtb, &
      & thr_in=thr1)

end subroutine test_p_fock_s2


subroutine test_p_fock_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(13) = [&
      &  8.56114797532967E-1_wp,  1.08081185051881E-1_wp,  8.40158726565170E-1_wp, &
      & -4.26774273578783E-2_wp, -9.50845571950689E-2_wp, -5.37262650096369E-1_wp, &
      &  4.55006307726409E-2_wp, -9.61301031913531E-2_wp, -5.37689488045761E-1_wp, &
      &  4.55254261185908E-2_wp, -9.54730110948796E-2_wp, -5.36576677043751E-1_wp, &
      &  4.55131478317075E-2_wp]

   real(wp), parameter :: density(43, 43, 2) = reshape([&
      & 6.91475853693089E-03_wp,-4.66439474604429E-04_wp, 7.43124629337322E-04_wp, &
      &-3.09904488020287E-04_wp,-2.72087829479343E-03_wp, 6.74488694982622E-03_wp, &
      &-5.28266720158092E-03_wp, 4.80689218345537E-03_wp, 9.22601160315454E-04_wp, &
      &-2.71966921747329E-03_wp,-1.01894320857934E-02_wp,-2.27063129355292E-03_wp, &
      &-2.82181425186610E-03_wp,-2.14613394982406E-03_wp, 3.61217559598682E-03_wp, &
      & 1.12028617140403E-02_wp,-8.64033744988582E-04_wp, 3.10371374558505E-02_wp, &
      & 2.30333249958885E-02_wp, 2.04795146556320E-02_wp, 1.17566904494576E-03_wp, &
      & 1.12383332750122E-03_wp,-3.37161261858860E-04_wp, 7.55695280205111E-04_wp, &
      &-4.91653601591487E-04_wp,-1.02695166967632E-03_wp,-3.45289149212074E-02_wp, &
      &-1.34121866165070E-02_wp, 2.17842958084033E-02_wp,-1.32239449246865E-03_wp, &
      & 9.23763188396611E-04_wp,-5.87060101043125E-04_wp,-6.06266471529292E-04_wp, &
      &-5.25957119782735E-04_wp,-9.44887800911054E-04_wp, 9.71328579482040E-03_wp, &
      &-1.32264514235801E-02_wp,-4.00506954181809E-02_wp,-7.30975427202173E-04_wp, &
      &-2.86085329171831E-04_wp,-5.83829373776639E-04_wp, 1.06473924707411E-03_wp, &
      & 1.20864032142143E-03_wp,-4.66439474604429E-04_wp, 2.37196189260497E-03_wp, &
      & 3.76291798637178E-04_wp,-1.89524702709497E-04_wp,-2.70383567603876E-03_wp, &
      & 7.51687695834723E-04_wp, 1.57852557800278E-05_wp,-1.70421269428447E-03_wp, &
      &-1.24711040739595E-04_wp, 7.89473512446235E-04_wp,-3.57098061462513E-04_wp, &
      & 4.01348859449634E-03_wp,-3.38878489220265E-03_wp, 6.34879336465873E-03_wp, &
      &-2.30645681571837E-03_wp,-6.70882068372886E-04_wp, 6.96944723684011E-03_wp, &
      &-1.08436093353111E-02_wp,-1.77208852887962E-02_wp,-1.93565751770493E-02_wp, &
      &-8.75401591255748E-04_wp,-5.65296230455721E-04_wp, 1.94865564588492E-04_wp, &
      &-6.09542770974166E-04_wp,-6.86637147630286E-05_wp,-7.77156635669540E-03_wp, &
      &-1.12023415897657E-02_wp,-1.78971328951955E-02_wp, 2.23980508202680E-02_wp, &
      &-8.58573131873999E-04_wp, 6.50886246445002E-04_wp,-1.23589021151666E-04_wp, &
      &-7.14602097206055E-04_wp, 1.61848868703983E-04_wp, 2.47681579670218E-03_wp, &
      & 1.51505155319696E-02_wp, 4.13303227393656E-03_wp, 1.19937863231962E-02_wp, &
      &-4.00974483671832E-04_wp,-1.44661987627843E-04_wp,-2.20905355380273E-05_wp, &
      &-4.45102163842455E-04_wp,-5.21718641131263E-04_wp, 7.43124629337322E-04_wp, &
      & 3.76291798637178E-04_wp, 1.92419514812543E-03_wp, 2.43411943306908E-04_wp, &
      &-1.38558534276645E-03_wp,-1.29700970258115E-03_wp, 3.20274034784215E-04_wp, &
      &-6.00771636389316E-04_wp, 4.53126654291574E-04_wp, 6.89363513188489E-04_wp, &
      & 2.85458472253175E-03_wp,-8.42769400826849E-03_wp,-3.49285659189717E-03_wp, &
      &-8.73740458821759E-03_wp,-2.03180597844227E-03_wp, 3.57585309496864E-03_wp, &
      & 4.50241444890764E-03_wp,-1.94202483147287E-02_wp, 1.08132873675669E-02_wp, &
      &-1.31156713246053E-02_wp,-6.35109734800467E-04_wp,-1.41358099308255E-04_wp, &
      & 4.86112152716323E-04_wp,-6.63799363346240E-05_wp, 2.63148272705645E-04_wp, &
      &-3.69338959899261E-03_wp,-1.62131073853099E-02_wp, 1.12706765571131E-02_wp, &
      & 1.05297118737550E-02_wp,-6.83652635177475E-04_wp,-3.67543219032771E-06_wp, &
      &-5.64274993109804E-04_wp,-4.99490489384980E-05_wp,-1.70869611262570E-04_wp, &
      &-3.50216085251382E-03_wp, 4.54245265618630E-03_wp, 1.19055309005404E-02_wp, &
      &-1.89552860618799E-02_wp,-4.21614601950282E-04_wp,-8.16969303104830E-05_wp, &
      &-5.05453807147833E-04_wp, 3.80106625809516E-05_wp, 4.70419178273609E-04_wp, &
      &-3.09904488020287E-04_wp,-1.89524702709497E-04_wp, 2.43411943306908E-04_wp, &
      & 2.47928402744496E-03_wp,-1.13559737080272E-03_wp,-1.51157026309616E-03_wp, &
      &-1.74075334963782E-05_wp, 1.78516116698550E-03_wp, 2.08500104756471E-03_wp, &
      &-2.62166403396570E-03_wp, 1.57889816861724E-03_wp, 5.04088894827932E-03_wp, &
      &-2.33600703355515E-03_wp, 1.82408462819675E-03_wp, 3.82256307528407E-03_wp, &
      &-1.43702087113228E-03_wp, 4.61273866776237E-03_wp,-1.97237097978807E-02_wp, &
      &-1.20991357049060E-02_wp, 5.08336438918823E-03_wp,-2.83097626905179E-04_wp, &
      &-5.91992885185732E-04_wp, 1.38111314515451E-04_wp,-9.22681999267372E-05_wp, &
      & 7.37417254727631E-04_wp, 5.31578919345921E-03_wp, 2.10073897134008E-02_wp, &
      & 1.09429414011331E-02_wp, 3.20167044114226E-03_wp, 3.34046150738591E-04_wp, &
      &-7.19136030754649E-04_wp, 4.42987781476963E-05_wp, 2.30642331218241E-04_wp, &
      & 6.80365098406377E-04_wp,-8.84877256219645E-03_wp, 1.35788864585981E-02_wp, &
      &-2.03601998000497E-02_wp,-2.27788852327071E-02_wp,-6.78210598851328E-04_wp, &
      &-4.08503618156095E-04_wp,-1.59367877771540E-04_wp, 1.00614012061061E-03_wp, &
      & 6.42501172500220E-04_wp,-2.72087829479343E-03_wp,-2.70383567603876E-03_wp, &
      &-1.38558534276645E-03_wp,-1.13559737080272E-03_wp, 4.40223717129710E-02_wp, &
      & 5.57984043827394E-03_wp, 2.63086726948185E-03_wp, 7.70835864592983E-03_wp, &
      &-3.97661022360440E-04_wp, 2.26637457682227E-04_wp, 7.74500745271912E-03_wp, &
      &-1.79565826935409E-03_wp, 1.67084160854114E-02_wp,-6.61485609393775E-03_wp, &
      &-7.67995123123247E-04_wp, 2.18809239737792E-03_wp, 2.24397055345448E-02_wp, &
      & 6.64549541656991E-02_wp, 8.85871421288148E-02_wp, 6.96854518048087E-04_wp, &
      & 9.90846322158204E-04_wp, 2.46975374062160E-03_wp, 3.12047632632326E-04_wp, &
      & 1.29596421479028E-03_wp,-1.29311438090923E-03_wp,-2.40673203813875E-02_wp, &
      & 7.93390435947543E-02_wp, 7.02903965683193E-02_wp,-1.26672314141576E-02_wp, &
      & 1.72311948156631E-03_wp,-2.33142012668952E-03_wp, 1.41048483943771E-04_wp, &
      & 1.23754087717857E-03_wp, 1.32235237491106E-03_wp,-1.35903661755529E-02_wp, &
      & 7.65937840551499E-02_wp, 3.28686383455122E-02_wp, 6.44067474820589E-02_wp, &
      &-1.00233070467192E-03_wp,-6.88819688392179E-04_wp, 3.53857312562318E-04_wp, &
      &-1.40694498187750E-03_wp,-2.36540545013463E-03_wp, 6.74488694982622E-03_wp, &
      & 7.51687695834723E-04_wp,-1.29700970258115E-03_wp,-1.51157026309616E-03_wp, &
      & 5.57984043827394E-03_wp, 3.14464127057781E-02_wp, 5.42414731300720E-03_wp, &
      &-4.96716029907645E-03_wp,-7.32631149807317E-03_wp, 3.52486701106247E-03_wp, &
      &-1.72853100047814E-02_wp, 6.66597965316477E-03_wp, 1.41280628812100E-03_wp, &
      & 1.69880913558385E-03_wp, 6.38096027475187E-04_wp, 9.57106430361156E-03_wp, &
      & 2.07478739672553E-02_wp, 6.75325136802193E-02_wp,-2.15309610650728E-02_wp, &
      & 7.89233538673288E-02_wp, 2.81701658164297E-03_wp, 5.88757207563303E-04_wp, &
      &-1.82030938334733E-03_wp, 5.97326378135036E-04_wp,-5.00617661061709E-04_wp, &
      & 1.84804474365580E-02_wp,-6.44185602640370E-02_wp, 3.79320344196431E-02_wp, &
      & 6.87221148902145E-02_wp,-2.66981711043064E-03_wp, 3.40407076459113E-04_wp, &
      &-2.12977682745461E-03_wp,-3.25448101994960E-04_wp,-6.09023766967790E-04_wp, &
      &-5.22282280973725E-03_wp, 3.75844989308498E-02_wp,-1.69469206158735E-02_wp, &
      & 4.36814617081137E-02_wp,-7.35276572441748E-04_wp,-6.06690581183132E-04_wp, &
      & 6.06034207881716E-04_wp, 1.10660155867559E-04_wp,-1.13219077295460E-03_wp, &
      &-5.28266720158092E-03_wp, 1.57852557800278E-05_wp, 3.20274034784215E-04_wp, &
      &-1.74075334963782E-05_wp, 2.63086726948185E-03_wp, 5.42414731300720E-03_wp, &
      & 2.93589728175237E-02_wp, 3.67709738831901E-03_wp,-7.47063038600942E-04_wp, &
      & 3.90836704631231E-03_wp, 7.09621953957305E-05_wp, 2.07760552454923E-03_wp, &
      &-1.29405483042985E-04_wp, 2.89718143665865E-03_wp, 3.16800401100101E-04_wp, &
      &-1.69606053549435E-02_wp,-7.54471811196135E-03_wp, 1.55318190585745E-02_wp, &
      &-9.64076933204775E-02_wp, 1.10515854057058E-02_wp, 1.45947727308444E-04_wp, &
      &-2.16103631753385E-03_wp,-1.41298732458153E-03_wp,-1.48077115539182E-03_wp, &
      &-3.68872484716546E-05_wp,-9.73129495503267E-03_wp, 3.62507400512151E-04_wp, &
      & 8.64415228683387E-02_wp,-5.52582904666052E-03_wp, 2.36297575249053E-04_wp, &
      &-1.94458556002731E-03_wp,-1.28072235876436E-03_wp, 1.67826970775164E-03_wp, &
      &-2.64680965054583E-04_wp,-9.91183630949522E-03_wp,-5.57175718616637E-03_wp, &
      & 8.67153633545511E-02_wp, 2.90410287474541E-03_wp, 3.78053496299335E-04_wp, &
      & 9.38622543052269E-04_wp,-1.22901953593348E-03_wp,-2.36862461809407E-03_wp, &
      & 3.42392873702847E-05_wp, 4.80689218345537E-03_wp,-1.70421269428447E-03_wp, &
      &-6.00771636389316E-04_wp, 1.78516116698550E-03_wp, 7.70835864592983E-03_wp, &
      &-4.96716029907645E-03_wp, 3.67709738831901E-03_wp, 3.44193006387385E-02_wp, &
      & 4.76228788479958E-03_wp,-1.12594099044268E-02_wp,-4.22637413712781E-03_wp, &
      &-1.10988128216367E-02_wp, 9.39053014220174E-04_wp,-1.16060082138239E-03_wp, &
      & 1.07624737461691E-02_wp, 5.60381954555915E-03_wp, 1.41316668863593E-02_wp, &
      & 7.88193671551598E-02_wp,-1.33489350096580E-02_wp, 5.13853337547666E-03_wp, &
      & 1.55817983245919E-03_wp, 5.12251036351069E-04_wp,-1.38894336043584E-03_wp, &
      & 4.98434304366330E-06_wp,-1.55155659357601E-03_wp,-1.19435034444843E-02_wp, &
      & 7.87707633106839E-02_wp,-3.00272463633385E-02_wp,-2.09804044213811E-03_wp, &
      & 1.14598547458640E-03_wp,-1.69273880682273E-04_wp, 1.52074117567960E-03_wp, &
      &-6.69143343856180E-04_wp, 1.75734486260661E-03_wp, 2.14421083889658E-02_wp, &
      & 3.28723251997892E-02_wp, 4.78972976104863E-02_wp,-1.04712017981930E-01_wp, &
      &-1.59454356614833E-03_wp,-7.29860323332650E-05_wp,-2.39582887844565E-03_wp, &
      & 4.78040749239259E-04_wp, 2.59547000701250E-03_wp, 9.22601160315454E-04_wp, &
      &-1.24711040739595E-04_wp, 4.53126654291574E-04_wp, 2.08500104756471E-03_wp, &
      &-3.97661022360440E-04_wp,-7.32631149807317E-03_wp,-7.47063038600942E-04_wp, &
      & 4.76228788479958E-03_wp, 4.33364941711876E-02_wp,-3.15071146365869E-03_wp, &
      &-6.67108487341014E-04_wp, 6.99176018564855E-03_wp,-6.22633014911391E-03_wp, &
      &-3.38055775918460E-03_wp, 7.69017399475255E-03_wp,-1.76961128943864E-03_wp, &
      &-8.85374052511688E-03_wp, 3.47268154939585E-02_wp,-3.45265330160837E-02_wp, &
      &-8.92468896301809E-02_wp,-1.30909318967325E-03_wp,-3.70938201880121E-04_wp, &
      &-1.10750429751085E-04_wp,-1.31927484389834E-03_wp,-1.78842876359275E-03_wp, &
      &-8.63307150232595E-03_wp,-3.79182506231813E-02_wp, 3.52704972597514E-02_wp, &
      &-9.70603235678255E-02_wp, 1.09039808896337E-03_wp,-3.59947499874213E-04_wp, &
      &-3.60824594445763E-04_wp, 1.85106237169424E-03_wp,-2.29662810010527E-03_wp, &
      & 2.20191489255149E-02_wp, 7.33216018913339E-02_wp,-7.19861577005277E-02_wp, &
      &-4.69373915752390E-02_wp,-2.10917504939904E-03_wp,-1.45840439450576E-03_wp, &
      & 1.08123286916086E-04_wp, 2.31640389940747E-03_wp, 5.16614467010584E-04_wp, &
      &-2.71966921747329E-03_wp, 7.89473512446235E-04_wp, 6.89363513188489E-04_wp, &
      &-2.62166403396570E-03_wp, 2.26637457682227E-04_wp, 3.52486701106247E-03_wp, &
      & 3.90836704631231E-03_wp,-1.12594099044268E-02_wp,-3.15071146365869E-03_wp, &
      & 2.34151430732665E-02_wp, 1.85415808471368E-02_wp,-4.24734078208696E-02_wp, &
      & 2.35575691060635E-03_wp,-3.33058390377248E-02_wp,-5.57795006945784E-03_wp, &
      & 1.01108773595180E-02_wp, 4.01102854486234E-04_wp,-4.56200941513286E-02_wp, &
      &-6.60382781517692E-03_wp, 4.54256036669546E-02_wp, 8.13365307011402E-04_wp, &
      &-1.45138153891304E-03_wp,-4.61944048289146E-04_wp, 2.51862372221059E-04_wp, &
      & 2.07065310902136E-03_wp, 3.41550940616051E-04_wp,-5.80400261978084E-02_wp, &
      &-1.74917522456008E-02_wp,-5.30554321998480E-02_wp, 5.62332461332891E-04_wp, &
      & 1.78951892335046E-03_wp, 5.49306638133438E-04_wp, 1.18354420616184E-04_wp, &
      &-2.60744135288374E-03_wp,-1.84935397426444E-04_wp, 4.00823788730330E-02_wp, &
      & 5.07256361602297E-02_wp, 6.01728319806439E-02_wp,-1.08074030233043E-03_wp, &
      &-5.41519705143916E-04_wp, 1.74128670059436E-04_wp,-2.34898043484794E-03_wp, &
      &-2.72927683099058E-03_wp,-1.01894320857934E-02_wp,-3.57098061462513E-04_wp, &
      & 2.85458472253175E-03_wp, 1.57889816861724E-03_wp, 7.74500745271912E-03_wp, &
      &-1.72853100047814E-02_wp, 7.09621953957305E-05_wp,-4.22637413712781E-03_wp, &
      &-6.67108487341014E-04_wp, 1.85415808471368E-02_wp, 7.81937061967329E-02_wp, &
      &-1.56952868336072E-01_wp, 4.63674964372455E-03_wp,-1.52201455595406E-01_wp, &
      &-3.58517678579157E-02_wp, 3.07585591489104E-02_wp, 4.52883324935330E-04_wp, &
      &-9.32147819443025E-02_wp,-5.48349175520555E-03_wp,-3.42959429035113E-02_wp, &
      &-1.17137609071555E-03_wp,-4.30791961534006E-03_wp,-8.16503856499543E-04_wp, &
      &-1.99648246102047E-03_wp, 1.66375144799991E-03_wp, 7.55268863773708E-05_wp, &
      & 7.36998183180566E-02_wp, 8.59929205815937E-03_wp,-3.29569302508286E-02_wp, &
      & 2.35295984832753E-03_wp,-7.98432581448504E-04_wp, 2.29759184125395E-03_wp, &
      &-6.70858028737459E-04_wp, 2.19545673457071E-03_wp, 9.16685372636403E-04_wp, &
      & 2.57126867852129E-02_wp, 1.14966740661808E-02_wp, 3.88665260793635E-02_wp, &
      &-9.20983705261979E-04_wp,-1.88651521176194E-03_wp, 2.10285595586614E-03_wp, &
      & 5.97064769855728E-04_wp,-3.32760114530298E-03_wp,-2.27063129355292E-03_wp, &
      & 4.01348859449634E-03_wp,-8.42769400826849E-03_wp, 5.04088894827932E-03_wp, &
      &-1.79565826935409E-03_wp, 6.66597965316477E-03_wp, 2.07760552454923E-03_wp, &
      &-1.10988128216367E-02_wp, 6.99176018564855E-03_wp,-4.24734078208696E-02_wp, &
      &-1.56952868336072E-01_wp, 4.62813275840075E-01_wp,-9.10464660712884E-03_wp, &
      & 4.38047333436907E-01_wp, 9.59547188065274E-02_wp,-1.25382771228993E-01_wp, &
      &-1.38374238063024E-03_wp,-7.41273093046619E-03_wp,-1.10468068242443E-02_wp, &
      &-1.56187513033571E-03_wp,-4.41907081312552E-03_wp, 6.66432321271109E-03_wp, &
      & 6.93367015806001E-03_wp, 3.50242290460517E-03_wp, 1.53391697946383E-03_wp, &
      & 5.43354709829560E-04_wp, 1.40237003801708E-02_wp, 3.28804253498567E-02_wp, &
      & 1.21189122314942E-02_wp,-6.66975468342735E-05_wp,-3.13907873267953E-03_wp, &
      &-3.25713970999458E-03_wp, 4.31017066297599E-03_wp,-1.47887531116780E-03_wp, &
      &-1.84087838429997E-04_wp, 1.41827982064466E-02_wp,-6.62740370090869E-02_wp, &
      & 2.39879765891764E-02_wp, 1.05566543858823E-03_wp, 3.60304872458680E-03_wp, &
      &-3.54139161533717E-03_wp,-2.95929713103497E-03_wp, 4.28984898099735E-03_wp, &
      &-2.82181425186610E-03_wp,-3.38878489220265E-03_wp,-3.49285659189717E-03_wp, &
      &-2.33600703355515E-03_wp, 1.67084160854114E-02_wp, 1.41280628812100E-03_wp, &
      &-1.29405483042985E-04_wp, 9.39053014220174E-04_wp,-6.22633014911391E-03_wp, &
      & 2.35575691060635E-03_wp, 4.63674964372455E-03_wp,-9.10464660712884E-03_wp, &
      & 1.55693316213329E-02_wp,-8.69456072981448E-03_wp,-2.07001602146450E-03_wp, &
      &-1.64641127464159E-03_wp, 7.40412926607099E-04_wp, 5.19987722272070E-02_wp, &
      & 3.50963776224181E-02_wp, 3.53821139878523E-02_wp, 2.03766503438101E-03_wp, &
      & 1.38182723046831E-03_wp,-8.56882698093754E-04_wp, 9.67443679062229E-04_wp, &
      &-7.89383163152358E-04_wp,-3.12016989036287E-04_wp, 4.89673716424762E-02_wp, &
      & 6.85566681167747E-03_wp,-4.33726604714240E-02_wp, 2.32152383555889E-03_wp, &
      &-7.19828442745710E-04_wp, 1.28190858547543E-03_wp, 7.83301665646748E-04_wp, &
      & 1.95937371857659E-04_wp,-2.19901739927516E-04_wp,-2.48692793280075E-02_wp, &
      & 6.56157086454086E-03_wp, 6.03966283304102E-02_wp, 1.66757073509143E-03_wp, &
      & 4.66111281358563E-04_wp, 1.29891259804992E-03_wp,-8.75032383977239E-04_wp, &
      &-1.61471980448458E-03_wp,-2.14613394982406E-03_wp, 6.34879336465873E-03_wp, &
      &-8.73740458821759E-03_wp, 1.82408462819675E-03_wp,-6.61485609393775E-03_wp, &
      & 1.69880913558385E-03_wp, 2.89718143665865E-03_wp,-1.16060082138239E-03_wp, &
      &-3.38055775918460E-03_wp,-3.33058390377248E-02_wp,-1.52201455595406E-01_wp, &
      & 4.38047333436907E-01_wp,-8.69456072981448E-03_wp, 4.36466806654643E-01_wp, &
      & 9.34741416346391E-02_wp,-1.19354104145719E-01_wp,-1.37121127203996E-03_wp, &
      & 5.30583941925366E-04_wp,-2.67897186120732E-02_wp, 4.59737348249017E-03_wp, &
      &-4.02177217100425E-03_wp, 6.15256712974783E-03_wp, 6.28787735430965E-03_wp, &
      & 3.19802362532034E-03_wp, 1.32955044996654E-03_wp, 1.31200798108916E-03_wp, &
      & 1.24744952822560E-02_wp,-7.32813987751117E-02_wp, 3.72566964772602E-03_wp, &
      & 1.02640693375478E-04_wp,-1.47429100083848E-05_wp,-1.17086414454252E-03_wp, &
      & 2.24601344727871E-03_wp,-2.00559116482125E-03_wp,-1.12586025604532E-03_wp, &
      & 2.39815493206550E-02_wp, 2.95058137726838E-02_wp, 2.61728970650068E-02_wp, &
      & 5.32051158914855E-04_wp, 4.14351970344133E-03_wp,-5.00522584021680E-03_wp, &
      &-6.01050825635358E-03_wp, 3.53071480903333E-03_wp, 3.61217559598682E-03_wp, &
      &-2.30645681571837E-03_wp,-2.03180597844227E-03_wp, 3.82256307528407E-03_wp, &
      &-7.67995123123247E-04_wp, 6.38096027475187E-04_wp, 3.16800401100101E-04_wp, &
      & 1.07624737461691E-02_wp, 7.69017399475255E-03_wp,-5.57795006945784E-03_wp, &
      &-3.58517678579157E-02_wp, 9.59547188065274E-02_wp,-2.07001602146450E-03_wp, &
      & 9.34741416346391E-02_wp, 3.90520218749009E-02_wp,-2.06735101088362E-02_wp, &
      &-2.07679312253469E-04_wp,-3.36427527894663E-03_wp,-2.68963264938291E-03_wp, &
      & 7.67642844476174E-02_wp, 1.35665085459366E-03_wp, 1.35542875736036E-03_wp, &
      & 7.30947572329957E-04_wp, 2.17212419187750E-03_wp, 2.27422974116782E-03_wp, &
      &-5.83309848279943E-04_wp, 1.33899406086197E-02_wp,-2.03214366796728E-03_wp, &
      &-5.90050199663154E-02_wp, 1.83884922533798E-03_wp,-5.27133217581021E-04_wp, &
      & 3.10277119126789E-04_wp, 1.78806056491888E-03_wp,-1.42996134356551E-03_wp, &
      & 9.38729390030271E-05_wp, 4.82654826897998E-02_wp,-1.14992339499985E-02_wp, &
      &-6.93910045606182E-02_wp,-1.97145781885766E-03_wp, 3.56893171984317E-05_wp, &
      &-2.34151302881454E-03_wp, 4.72722357649428E-04_wp, 2.70282185999502E-03_wp, &
      & 1.12028617140403E-02_wp,-6.70882068372886E-04_wp, 3.57585309496864E-03_wp, &
      &-1.43702087113228E-03_wp, 2.18809239737792E-03_wp, 9.57106430361156E-03_wp, &
      &-1.69606053549435E-02_wp, 5.60381954555915E-03_wp,-1.76961128943864E-03_wp, &
      & 1.01108773595180E-02_wp, 3.07585591489104E-02_wp,-1.25382771228993E-01_wp, &
      &-1.64641127464159E-03_wp,-1.19354104145719E-01_wp,-2.06735101088362E-02_wp, &
      & 5.64165987850187E-02_wp,-3.60508414739047E-04_wp, 2.36891004261134E-02_wp, &
      & 7.45384160365845E-02_wp, 3.27358666273189E-02_wp, 2.50498630196250E-03_wp, &
      & 3.12765956430376E-04_wp,-1.49064981297606E-03_wp, 7.64346376852005E-04_wp, &
      &-3.81834224017286E-04_wp,-3.89914338131860E-04_wp,-5.16374969915572E-02_wp, &
      &-4.56666343915549E-02_wp, 4.41302193900407E-02_wp,-2.26312525917335E-03_wp, &
      & 2.64354771088890E-03_wp, 1.76499596075038E-04_wp,-2.64863626414929E-03_wp, &
      &-2.83279607122441E-05_wp,-2.41011357525776E-04_wp, 5.49357483690815E-02_wp, &
      &-2.23923554737875E-02_wp,-3.95054115591359E-02_wp,-2.59058745651069E-03_wp, &
      &-2.26084037291830E-03_wp, 4.94952180334675E-04_wp, 2.52712586788310E-03_wp, &
      &-5.02041178157069E-04_wp,-8.64033744988582E-04_wp, 6.96944723684011E-03_wp, &
      & 4.50241444890764E-03_wp, 4.61273866776237E-03_wp, 2.24397055345448E-02_wp, &
      & 2.07478739672553E-02_wp,-7.54471811196135E-03_wp, 1.41316668863593E-02_wp, &
      &-8.85374052511688E-03_wp, 4.01102854486234E-04_wp, 4.52883324935330E-04_wp, &
      &-1.38374238063024E-03_wp, 7.40412926607099E-04_wp,-1.37121127203996E-03_wp, &
      &-2.07679312253469E-04_wp,-3.60508414739047E-04_wp, 9.79385233275976E-01_wp, &
      &-6.73932117684685E-02_wp,-4.24221242518916E-02_wp,-4.54829773946277E-02_wp, &
      &-3.44767227988470E-03_wp,-3.20186908954796E-03_wp, 1.08518088940097E-03_wp, &
      &-2.17906881067377E-03_wp, 1.36803923011790E-03_wp,-1.74429112109240E-03_wp, &
      &-4.60039554543272E-03_wp,-4.62867164242500E-03_wp,-2.29461193380882E-02_wp, &
      & 1.08781582530049E-03_wp, 8.21284940310410E-04_wp, 8.63411726539932E-04_wp, &
      & 3.71234417790668E-04_wp,-2.26885600406139E-03_wp,-4.67052712135014E-04_wp, &
      &-2.18177810703358E-02_wp,-4.85160230143877E-03_wp, 2.44761467322725E-03_wp, &
      & 2.10650456051350E-03_wp, 5.13537734107091E-04_wp, 8.26572953275499E-04_wp, &
      & 7.14219643150120E-04_wp, 1.14929051465742E-03_wp, 3.10371374558505E-02_wp, &
      &-1.08436093353111E-02_wp,-1.94202483147287E-02_wp,-1.97237097978807E-02_wp, &
      & 6.64549541656991E-02_wp, 6.75325136802193E-02_wp, 1.55318190585745E-02_wp, &
      & 7.88193671551598E-02_wp, 3.47268154939585E-02_wp,-4.56200941513286E-02_wp, &
      &-9.32147819443025E-02_wp,-7.41273093046619E-03_wp, 5.19987722272070E-02_wp, &
      & 5.30583941925366E-04_wp,-3.36427527894663E-03_wp, 2.36891004261134E-02_wp, &
      &-6.73932117684685E-02_wp, 8.89614738721909E-01_wp,-1.28950331968545E-02_wp, &
      &-1.39586534563915E-02_wp, 1.53963342226206E-02_wp, 1.29909255553247E-02_wp, &
      &-1.31338339197474E-02_wp,-1.58767803767510E-05_wp,-2.41453153734081E-02_wp, &
      & 5.83630129305774E-03_wp, 4.74494717503171E-03_wp,-1.32063220713228E-02_wp, &
      &-1.19859393613625E-02_wp,-6.98434276117959E-04_wp, 2.29788467998815E-03_wp, &
      & 3.50868101514438E-04_wp,-9.19684784727941E-04_wp,-1.01534697329913E-03_wp, &
      &-1.79249305334630E-02_wp,-4.44215180079381E-02_wp,-4.45142330258693E-03_wp, &
      & 3.22808677054596E-02_wp, 5.51511024939793E-03_wp, 2.23194762526908E-03_wp, &
      & 2.14969208354250E-03_wp, 2.26334563100816E-04_wp,-2.69348432554597E-04_wp, &
      & 2.30333249958885E-02_wp,-1.77208852887962E-02_wp, 1.08132873675669E-02_wp, &
      &-1.20991357049060E-02_wp, 8.85871421288148E-02_wp,-2.15309610650728E-02_wp, &
      &-9.64076933204775E-02_wp,-1.33489350096580E-02_wp,-3.45265330160837E-02_wp, &
      &-6.60382781517692E-03_wp,-5.48349175520555E-03_wp,-1.10468068242443E-02_wp, &
      & 3.50963776224181E-02_wp,-2.67897186120732E-02_wp,-2.68963264938291E-03_wp, &
      & 7.45384160365845E-02_wp,-4.24221242518916E-02_wp,-1.28950331968545E-02_wp, &
      & 9.00965756211378E-01_wp,-9.14415453342455E-03_wp, 1.33841889823395E-03_wp, &
      & 2.15951288119434E-02_wp, 1.24101348943107E-02_wp, 1.49920354705176E-02_wp, &
      &-4.40469099734259E-04_wp, 9.85577803658296E-04_wp,-1.22201335360621E-02_wp, &
      & 2.11341370624279E-02_wp,-1.39791945191471E-02_wp,-7.01721881187569E-04_wp, &
      & 1.08687409942297E-04_wp,-1.44491706267255E-03_wp, 4.54075284812451E-04_wp, &
      &-1.98967421504811E-03_wp, 1.69947916252597E-03_wp,-1.66229817006137E-02_wp, &
      & 2.23035839447183E-02_wp,-9.27699182867206E-03_wp, 4.79918109486991E-04_wp, &
      & 1.40000483980811E-04_wp,-1.14187690596524E-03_wp, 2.76159075773128E-04_wp, &
      & 1.54213587967972E-03_wp, 2.04795146556320E-02_wp,-1.93565751770493E-02_wp, &
      &-1.31156713246053E-02_wp, 5.08336438918823E-03_wp, 6.96854518048087E-04_wp, &
      & 7.89233538673288E-02_wp, 1.10515854057058E-02_wp, 5.13853337547666E-03_wp, &
      &-8.92468896301809E-02_wp, 4.54256036669546E-02_wp,-3.42959429035113E-02_wp, &
      &-1.56187513033571E-03_wp, 3.53821139878523E-02_wp, 4.59737348249017E-03_wp, &
      & 7.67642844476174E-02_wp, 3.27358666273189E-02_wp,-4.54829773946277E-02_wp, &
      &-1.39586534563915E-02_wp,-9.14415453342455E-03_wp, 9.01197413335997E-01_wp, &
      & 2.37424843741826E-02_wp, 2.21729555528611E-04_wp,-8.75069972739857E-03_wp, &
      & 1.29168909226698E-02_wp, 1.70306110002549E-02_wp,-2.30415156877361E-02_wp, &
      & 1.66719242999167E-02_wp, 2.25021243598935E-03_wp,-5.05076972762118E-02_wp, &
      & 5.43204781194369E-03_wp,-9.57154239907275E-05_wp, 2.35591709642219E-03_wp, &
      & 2.63084615326781E-03_wp,-2.80062748714786E-03_wp, 1.34259756367277E-02_wp, &
      & 9.67852731639095E-03_wp,-1.30717775429154E-02_wp,-3.47522554508627E-03_wp, &
      &-2.41917211488525E-03_wp,-1.02717741612734E-03_wp,-7.37519467381799E-04_wp, &
      & 2.37177085329030E-03_wp, 1.59044196104693E-03_wp, 1.17566904494576E-03_wp, &
      &-8.75401591255748E-04_wp,-6.35109734800467E-04_wp,-2.83097626905179E-04_wp, &
      & 9.90846322158204E-04_wp, 2.81701658164297E-03_wp, 1.45947727308444E-04_wp, &
      & 1.55817983245919E-03_wp,-1.30909318967325E-03_wp, 8.13365307011402E-04_wp, &
      &-1.17137609071555E-03_wp,-4.41907081312552E-03_wp, 2.03766503438101E-03_wp, &
      &-4.02177217100425E-03_wp, 1.35665085459366E-03_wp, 2.50498630196250E-03_wp, &
      &-3.44767227988470E-03_wp, 1.53963342226206E-02_wp, 1.33841889823395E-03_wp, &
      & 2.37424843741826E-02_wp, 9.71080523129421E-04_wp, 2.24537932149413E-04_wp, &
      &-5.12935597024076E-04_wp, 3.54826393035252E-04_wp, 1.78175387129317E-05_wp, &
      &-9.92864289498025E-04_wp,-4.93805000442882E-04_wp,-1.05118342010505E-03_wp, &
      &-5.15532684545041E-03_wp, 2.01568567477248E-04_wp, 8.47489533147879E-05_wp, &
      & 1.28915036185613E-04_wp, 5.96716941131621E-05_wp,-1.61436989305845E-04_wp, &
      & 9.05992698758370E-04_wp,-1.86454186027926E-03_wp,-1.56221941949092E-03_wp, &
      &-8.01779151306951E-04_wp, 5.29754230191627E-05_wp,-1.67720095156313E-05_wp, &
      & 6.88550780578720E-05_wp, 1.49459436391235E-04_wp, 3.48745958377767E-05_wp, &
      & 1.12383332750122E-03_wp,-5.65296230455721E-04_wp,-1.41358099308255E-04_wp, &
      &-5.91992885185732E-04_wp, 2.46975374062160E-03_wp, 5.88757207563303E-04_wp, &
      &-2.16103631753385E-03_wp, 5.12251036351069E-04_wp,-3.70938201880121E-04_wp, &
      &-1.45138153891304E-03_wp,-4.30791961534006E-03_wp, 6.66432321271109E-03_wp, &
      & 1.38182723046831E-03_wp, 6.15256712974783E-03_wp, 1.35542875736036E-03_wp, &
      & 3.12765956430376E-04_wp,-3.20186908954796E-03_wp, 1.29909255553247E-02_wp, &
      & 2.15951288119434E-02_wp, 2.21729555528611E-04_wp, 2.24537932149413E-04_wp, &
      & 8.41821697249856E-04_wp, 2.01357871247298E-04_wp, 4.33545604699002E-04_wp, &
      &-3.54326607641147E-04_wp, 6.92419916188416E-04_wp,-2.26880794835676E-03_wp, &
      &-2.60472690045053E-04_wp,-6.40001479018286E-04_wp,-5.97084418734319E-05_wp, &
      & 4.30725140318215E-05_wp,-9.13448849020790E-05_wp, 5.35251386103987E-05_wp, &
      &-1.48760158877089E-04_wp,-3.65628894303910E-04_wp,-3.34642237401371E-03_wp, &
      &-3.20723614997892E-04_wp,-3.02716907724124E-04_wp, 1.66642160459605E-04_wp, &
      & 1.29518345211162E-04_wp,-6.06645445272479E-05_wp,-3.58157141419361E-05_wp, &
      & 1.47884993198466E-04_wp,-3.37161261858860E-04_wp, 1.94865564588492E-04_wp, &
      & 4.86112152716323E-04_wp, 1.38111314515451E-04_wp, 3.12047632632326E-04_wp, &
      &-1.82030938334733E-03_wp,-1.41298732458153E-03_wp,-1.38894336043584E-03_wp, &
      &-1.10750429751085E-04_wp,-4.61944048289146E-04_wp,-8.16503856499543E-04_wp, &
      & 6.93367015806001E-03_wp,-8.56882698093754E-04_wp, 6.28787735430965E-03_wp, &
      & 7.30947572329957E-04_wp,-1.49064981297606E-03_wp, 1.08518088940097E-03_wp, &
      &-1.31338339197474E-02_wp, 1.24101348943107E-02_wp,-8.75069972739857E-03_wp, &
      &-5.12935597024076E-04_wp, 2.01357871247298E-04_wp, 5.57704342370171E-04_wp, &
      & 1.30088685793062E-04_wp, 2.05807133619415E-04_wp, 7.04869867156499E-04_wp, &
      &-4.71965410425160E-04_wp, 1.54322398343414E-03_wp, 1.57912740212179E-03_wp, &
      &-8.60224070514466E-05_wp,-8.26444834400528E-05_wp,-1.22481720955816E-04_wp, &
      & 5.33348705143874E-05_wp,-3.16133131285991E-06_wp, 6.83482385414134E-04_wp, &
      & 1.43421135623885E-03_wp, 1.29954197473103E-03_wp,-6.84034941761651E-04_wp, &
      &-6.77514813138002E-05_wp, 3.57114046695395E-05_wp,-1.37991308366822E-04_wp, &
      &-1.10652382975503E-04_wp, 8.18569440148499E-05_wp, 7.55695280205111E-04_wp, &
      &-6.09542770974166E-04_wp,-6.63799363346240E-05_wp,-9.22681999267372E-05_wp, &
      & 1.29596421479028E-03_wp, 5.97326378135036E-04_wp,-1.48077115539182E-03_wp, &
      & 4.98434304366330E-06_wp,-1.31927484389834E-03_wp, 2.51862372221059E-04_wp, &
      &-1.99648246102047E-03_wp, 3.50242290460517E-03_wp, 9.67443679062229E-04_wp, &
      & 3.19802362532034E-03_wp, 2.17212419187750E-03_wp, 7.64346376852005E-04_wp, &
      &-2.17906881067377E-03_wp,-1.58767803767510E-05_wp, 1.49920354705176E-02_wp, &
      & 1.29168909226698E-02_wp, 3.54826393035252E-04_wp, 4.33545604699002E-04_wp, &
      & 1.30088685793062E-04_wp, 4.82820109679038E-04_wp, 2.51108863693032E-04_wp, &
      &-7.79172865408991E-04_wp,-7.66705658343985E-04_wp, 2.26017817015977E-04_wp, &
      &-3.30674201547611E-03_wp, 1.03484025350007E-04_wp,-1.06917681599775E-06_wp, &
      &-1.97825617009888E-07_wp, 1.08987968326159E-04_wp,-1.54513083585034E-04_wp, &
      & 9.18499763897787E-04_wp, 2.49064579616363E-04_wp,-3.01515131459455E-04_wp, &
      &-1.87144089307596E-03_wp,-4.15809966344438E-05_wp, 1.50166180255650E-05_wp, &
      &-9.17945277861075E-05_wp, 3.21446802815042E-05_wp, 1.29877997002024E-04_wp, &
      &-4.91653601591487E-04_wp,-6.86637147630286E-05_wp, 2.63148272705645E-04_wp, &
      & 7.37417254727631E-04_wp,-1.29311438090923E-03_wp,-5.00617661061709E-04_wp, &
      &-3.68872484716546E-05_wp,-1.55155659357601E-03_wp,-1.78842876359275E-03_wp, &
      & 2.07065310902136E-03_wp, 1.66375144799991E-03_wp, 1.53391697946383E-03_wp, &
      &-7.89383163152358E-04_wp, 1.32955044996654E-03_wp, 2.27422974116782E-03_wp, &
      &-3.81834224017286E-04_wp, 1.36803923011790E-03_wp,-2.41453153734081E-02_wp, &
      &-4.40469099734259E-04_wp, 1.70306110002549E-02_wp, 1.78175387129317E-05_wp, &
      &-3.54326607641147E-04_wp, 2.05807133619415E-04_wp, 2.51108863693032E-04_wp, &
      & 1.00296300268016E-03_wp,-2.26339054665822E-03_wp, 1.22673110868481E-03_wp, &
      & 1.42194057769481E-03_wp,-2.96462717477230E-03_wp, 1.86342724132444E-04_wp, &
      &-1.01495495974295E-04_wp, 4.60297941198317E-05_wp, 1.31730529442160E-04_wp, &
      &-5.02008128431916E-05_wp, 2.19806549831237E-03_wp, 5.36800952916668E-03_wp, &
      &-1.11865135945953E-04_wp,-2.11155401443314E-03_wp,-3.09423947555407E-04_wp, &
      &-1.21929197698885E-04_wp,-1.22852682598501E-04_wp, 4.05393279044262E-05_wp, &
      & 4.47136744943953E-05_wp,-1.02695166967632E-03_wp,-7.77156635669540E-03_wp, &
      &-3.69338959899261E-03_wp, 5.31578919345921E-03_wp,-2.40673203813875E-02_wp, &
      & 1.84804474365580E-02_wp,-9.73129495503267E-03_wp,-1.19435034444843E-02_wp, &
      &-8.63307150232595E-03_wp, 3.41550940616051E-04_wp, 7.55268863773708E-05_wp, &
      & 5.43354709829560E-04_wp,-3.12016989036287E-04_wp, 1.31200798108916E-03_wp, &
      &-5.83309848279943E-04_wp,-3.89914338131860E-04_wp,-1.74429112109240E-03_wp, &
      & 5.83630129305774E-03_wp, 9.85577803658296E-04_wp,-2.30415156877361E-02_wp, &
      &-9.92864289498025E-04_wp, 6.92419916188416E-04_wp, 7.04869867156499E-04_wp, &
      &-7.79172865408991E-04_wp,-2.26339054665822E-03_wp, 9.80319653876714E-01_wp, &
      & 6.96148081109739E-02_wp, 3.39853186201137E-02_wp,-4.74581076068337E-02_wp, &
      & 3.77653932652171E-03_wp,-2.82491840284833E-03_wp, 1.53601953767167E-03_wp, &
      & 1.88938374714941E-03_wp, 1.39873724097286E-03_wp,-3.03607735046960E-03_wp, &
      & 1.59310863120097E-02_wp, 1.87885921811317E-02_wp, 5.08467535689675E-03_wp, &
      &-1.69023816854923E-03_wp,-6.52194452940651E-04_wp,-1.24575925065985E-03_wp, &
      &-1.84891013913080E-03_wp,-3.95093353419281E-04_wp,-3.45289149212074E-02_wp, &
      &-1.12023415897657E-02_wp,-1.62131073853099E-02_wp, 2.10073897134008E-02_wp, &
      & 7.93390435947543E-02_wp,-6.44185602640370E-02_wp, 3.62507400512151E-04_wp, &
      & 7.87707633106839E-02_wp,-3.79182506231813E-02_wp,-5.80400261978084E-02_wp, &
      & 7.36998183180566E-02_wp, 1.40237003801708E-02_wp, 4.89673716424762E-02_wp, &
      & 1.24744952822560E-02_wp, 1.33899406086197E-02_wp,-5.16374969915572E-02_wp, &
      &-4.60039554543272E-03_wp, 4.74494717503171E-03_wp,-1.22201335360621E-02_wp, &
      & 1.66719242999167E-02_wp,-4.93805000442882E-04_wp,-2.26880794835676E-03_wp, &
      &-4.71965410425160E-04_wp,-7.66705658343985E-04_wp, 1.22673110868481E-03_wp, &
      & 6.96148081109739E-02_wp, 8.84483522477083E-01_wp,-6.27087331785375E-03_wp, &
      & 1.39890422776069E-02_wp, 1.70641049794487E-02_wp,-1.41030186631501E-02_wp, &
      & 1.27499204873507E-02_wp, 3.86236907823217E-04_wp, 2.29766788920185E-02_wp, &
      & 9.93178929004829E-03_wp,-2.68014119423897E-03_wp,-3.97830522407406E-02_wp, &
      &-1.26975405205949E-02_wp, 1.43642412751397E-03_wp,-5.15623780968630E-04_wp, &
      & 1.97029419195767E-03_wp, 3.39979714316987E-03_wp, 1.57372281668683E-04_wp, &
      &-1.34121866165070E-02_wp,-1.78971328951955E-02_wp, 1.12706765571131E-02_wp, &
      & 1.09429414011331E-02_wp, 7.02903965683193E-02_wp, 3.79320344196431E-02_wp, &
      & 8.64415228683387E-02_wp,-3.00272463633385E-02_wp, 3.52704972597514E-02_wp, &
      &-1.74917522456008E-02_wp, 8.59929205815937E-03_wp, 3.28804253498567E-02_wp, &
      & 6.85566681167747E-03_wp,-7.32813987751117E-02_wp,-2.03214366796728E-03_wp, &
      &-4.56666343915549E-02_wp,-4.62867164242500E-03_wp,-1.32063220713228E-02_wp, &
      & 2.11341370624279E-02_wp, 2.25021243598935E-03_wp,-1.05118342010505E-03_wp, &
      &-2.60472690045053E-04_wp, 1.54322398343414E-03_wp, 2.26017817015977E-04_wp, &
      & 1.42194057769481E-03_wp, 3.39853186201137E-02_wp,-6.27087331785375E-03_wp, &
      & 9.01725438066334E-01_wp, 8.33818371870478E-03_wp,-3.41304494414273E-04_wp, &
      &-2.29885786612176E-02_wp,-1.59420543751911E-02_wp, 1.67323852681780E-02_wp, &
      & 5.90543543263457E-04_wp, 1.89070265565559E-02_wp,-2.18204966861355E-02_wp, &
      &-1.64116449632517E-02_wp,-3.50049447182234E-02_wp, 1.92048599024827E-03_wp, &
      & 2.73294147837168E-04_wp, 2.40805559719466E-04_wp, 3.98704788703455E-03_wp, &
      & 3.83007723351944E-03_wp, 2.17842958084033E-02_wp, 2.23980508202680E-02_wp, &
      & 1.05297118737550E-02_wp, 3.20167044114226E-03_wp,-1.26672314141576E-02_wp, &
      & 6.87221148902145E-02_wp,-5.52582904666052E-03_wp,-2.09804044213811E-03_wp, &
      &-9.70603235678255E-02_wp,-5.30554321998480E-02_wp,-3.29569302508286E-02_wp, &
      & 1.21189122314942E-02_wp,-4.33726604714240E-02_wp, 3.72566964772602E-03_wp, &
      &-5.90050199663154E-02_wp, 4.41302193900407E-02_wp,-2.29461193380882E-02_wp, &
      &-1.19859393613625E-02_wp,-1.39791945191471E-02_wp,-5.05076972762118E-02_wp, &
      &-5.15532684545041E-03_wp,-6.40001479018286E-04_wp, 1.57912740212179E-03_wp, &
      &-3.30674201547611E-03_wp,-2.96462717477230E-03_wp,-4.74581076068337E-02_wp, &
      & 1.39890422776069E-02_wp, 8.33818371870478E-03_wp, 9.00709592745117E-01_wp, &
      &-2.32277295829378E-02_wp,-5.91310561324898E-04_wp,-1.03765797970804E-02_wp, &
      &-1.37570262673697E-02_wp, 1.79170712797512E-02_wp, 1.33699745251597E-02_wp, &
      & 1.12806043159290E-02_wp,-9.44493581999449E-03_wp,-1.90246054656196E-03_wp, &
      &-1.78419681380516E-03_wp,-9.11286541924727E-04_wp,-5.17440346509960E-04_wp, &
      & 2.69307620534643E-03_wp, 1.58235914556952E-03_wp,-1.32239449246865E-03_wp, &
      &-8.58573131873999E-04_wp,-6.83652635177475E-04_wp, 3.34046150738591E-04_wp, &
      & 1.72311948156631E-03_wp,-2.66981711043064E-03_wp, 2.36297575249053E-04_wp, &
      & 1.14598547458640E-03_wp, 1.09039808896337E-03_wp, 5.62332461332891E-04_wp, &
      & 2.35295984832753E-03_wp,-6.66975468342735E-05_wp, 2.32152383555889E-03_wp, &
      & 1.02640693375478E-04_wp, 1.83884922533798E-03_wp,-2.26312525917335E-03_wp, &
      & 1.08781582530049E-03_wp,-6.98434276117959E-04_wp,-7.01721881187569E-04_wp, &
      & 5.43204781194369E-03_wp, 2.01568567477248E-04_wp,-5.97084418734319E-05_wp, &
      &-8.60224070514466E-05_wp, 1.03484025350007E-04_wp, 1.86342724132444E-04_wp, &
      & 3.77653932652171E-03_wp, 1.70641049794487E-02_wp,-3.41304494414273E-04_wp, &
      &-2.32277295829378E-02_wp, 9.71954797250523E-04_wp,-2.66463312336675E-04_wp, &
      & 5.36546546854730E-04_wp, 3.78261896706480E-04_wp,-2.20377842347737E-05_wp, &
      &-1.08352826422085E-03_wp,-1.82074210885778E-03_wp,-1.01687330466814E-03_wp, &
      & 1.84168810430206E-03_wp, 1.24242288971551E-04_wp, 2.60909585424519E-05_wp, &
      & 9.82423307899891E-05_wp,-1.33251711439347E-05_wp,-7.80121783922660E-05_wp, &
      & 9.23763188396611E-04_wp, 6.50886246445002E-04_wp,-3.67543219032771E-06_wp, &
      &-7.19136030754649E-04_wp,-2.33142012668952E-03_wp, 3.40407076459113E-04_wp, &
      &-1.94458556002731E-03_wp,-1.69273880682273E-04_wp,-3.59947499874213E-04_wp, &
      & 1.78951892335046E-03_wp,-7.98432581448504E-04_wp,-3.13907873267953E-03_wp, &
      &-7.19828442745710E-04_wp,-1.47429100083848E-05_wp,-5.27133217581021E-04_wp, &
      & 2.64354771088890E-03_wp, 8.21284940310410E-04_wp, 2.29788467998815E-03_wp, &
      & 1.08687409942297E-04_wp,-9.57154239907275E-05_wp, 8.47489533147879E-05_wp, &
      & 4.30725140318215E-05_wp,-8.26444834400528E-05_wp,-1.06917681599775E-06_wp, &
      &-1.01495495974295E-04_wp,-2.82491840284833E-03_wp,-1.41030186631501E-02_wp, &
      &-2.29885786612176E-02_wp,-5.91310561324898E-04_wp,-2.66463312336675E-04_wp, &
      & 8.50416540479365E-04_wp, 2.20268366240989E-04_wp,-4.57895934641153E-04_wp, &
      &-3.97935121573486E-04_wp,-1.96526810804084E-03_wp, 3.13903002693610E-03_wp, &
      & 3.85187428687074E-03_wp, 2.11734626341689E-03_wp,-1.34264453636114E-04_wp, &
      &-2.51838077169850E-05_wp,-5.85408576231827E-05_wp,-2.31098244614289E-04_wp, &
      &-1.74472472563144E-04_wp,-5.87060101043125E-04_wp,-1.23589021151666E-04_wp, &
      &-5.64274993109804E-04_wp, 4.42987781476963E-05_wp, 1.41048483943771E-04_wp, &
      &-2.12977682745461E-03_wp,-1.28072235876436E-03_wp, 1.52074117567960E-03_wp, &
      &-3.60824594445763E-04_wp, 5.49306638133438E-04_wp, 2.29759184125395E-03_wp, &
      &-3.25713970999458E-03_wp, 1.28190858547543E-03_wp,-1.17086414454252E-03_wp, &
      & 3.10277119126789E-04_wp, 1.76499596075038E-04_wp, 8.63411726539932E-04_wp, &
      & 3.50868101514438E-04_wp,-1.44491706267255E-03_wp, 2.35591709642219E-03_wp, &
      & 1.28915036185613E-04_wp,-9.13448849020790E-05_wp,-1.22481720955816E-04_wp, &
      &-1.97825617009888E-07_wp, 4.60297941198317E-05_wp, 1.53601953767167E-03_wp, &
      & 1.27499204873507E-02_wp,-1.59420543751911E-02_wp,-1.03765797970804E-02_wp, &
      & 5.36546546854730E-04_wp, 2.20268366240989E-04_wp, 6.09265712052873E-04_wp, &
      &-1.47519943592550E-04_wp, 1.22238728620312E-04_wp,-1.25095952083955E-03_wp, &
      & 2.27846314199926E-04_wp, 1.88274865285337E-04_wp, 2.05809775538211E-03_wp, &
      & 1.14606708088033E-05_wp,-2.53947944807633E-05_wp, 7.94171880800760E-05_wp, &
      &-5.67477088132154E-05_wp,-1.57956030126034E-04_wp,-6.06266471529292E-04_wp, &
      &-7.14602097206055E-04_wp,-4.99490489384980E-05_wp, 2.30642331218241E-04_wp, &
      & 1.23754087717857E-03_wp,-3.25448101994960E-04_wp, 1.67826970775164E-03_wp, &
      &-6.69143343856180E-04_wp, 1.85106237169424E-03_wp, 1.18354420616184E-04_wp, &
      &-6.70858028737459E-04_wp, 4.31017066297599E-03_wp, 7.83301665646748E-04_wp, &
      & 2.24601344727871E-03_wp, 1.78806056491888E-03_wp,-2.64863626414929E-03_wp, &
      & 3.71234417790668E-04_wp,-9.19684784727941E-04_wp, 4.54075284812451E-04_wp, &
      & 2.63084615326781E-03_wp, 5.96716941131621E-05_wp, 5.35251386103987E-05_wp, &
      & 5.33348705143874E-05_wp, 1.08987968326159E-04_wp, 1.31730529442160E-04_wp, &
      & 1.88938374714941E-03_wp, 3.86236907823217E-04_wp, 1.67323852681780E-02_wp, &
      &-1.37570262673697E-02_wp, 3.78261896706480E-04_wp,-4.57895934641153E-04_wp, &
      &-1.47519943592550E-04_wp, 5.67070951884331E-04_wp,-2.66365841732221E-04_wp, &
      &-2.17685295256138E-05_wp,-2.21004109845814E-03_wp,-1.02849359838044E-03_wp, &
      &-9.02442359029762E-04_wp, 1.12455065949355E-04_wp, 6.84526285750054E-05_wp, &
      &-1.17897739759140E-05_wp, 2.91415636443471E-05_wp, 1.16342256034267E-04_wp, &
      &-5.25957119782735E-04_wp, 1.61848868703983E-04_wp,-1.70869611262570E-04_wp, &
      & 6.80365098406377E-04_wp, 1.32235237491106E-03_wp,-6.09023766967790E-04_wp, &
      &-2.64680965054583E-04_wp, 1.75734486260661E-03_wp,-2.29662810010527E-03_wp, &
      &-2.60744135288374E-03_wp, 2.19545673457071E-03_wp,-1.47887531116780E-03_wp, &
      & 1.95937371857659E-04_wp,-2.00559116482125E-03_wp,-1.42996134356551E-03_wp, &
      &-2.83279607122441E-05_wp,-2.26885600406139E-03_wp,-1.01534697329913E-03_wp, &
      &-1.98967421504811E-03_wp,-2.80062748714786E-03_wp,-1.61436989305845E-04_wp, &
      &-1.48760158877089E-04_wp,-3.16133131285991E-06_wp,-1.54513083585034E-04_wp, &
      &-5.02008128431916E-05_wp, 1.39873724097286E-03_wp, 2.29766788920185E-02_wp, &
      & 5.90543543263457E-04_wp, 1.79170712797512E-02_wp,-2.20377842347737E-05_wp, &
      &-3.97935121573486E-04_wp, 1.22238728620312E-04_wp,-2.66365841732221E-04_wp, &
      & 9.76918208351756E-04_wp, 1.37310447814570E-03_wp,-3.85050874145553E-04_wp, &
      &-4.29155366733989E-03_wp,-8.98951630102113E-04_wp, 7.11466472221123E-06_wp, &
      &-7.22055400537447E-05_wp, 1.11585039120817E-04_wp, 2.53643008087482E-04_wp, &
      & 2.63123070323647E-05_wp,-9.44887800911054E-04_wp, 2.47681579670218E-03_wp, &
      &-3.50216085251382E-03_wp,-8.84877256219645E-03_wp,-1.35903661755529E-02_wp, &
      &-5.22282280973725E-03_wp,-9.91183630949522E-03_wp, 2.14421083889658E-02_wp, &
      & 2.20191489255149E-02_wp,-1.84935397426444E-04_wp, 9.16685372636403E-04_wp, &
      &-1.84087838429997E-04_wp,-2.19901739927516E-04_wp,-1.12586025604532E-03_wp, &
      & 9.38729390030271E-05_wp,-2.41011357525776E-04_wp,-4.67052712135014E-04_wp, &
      &-1.79249305334630E-02_wp, 1.69947916252597E-03_wp, 1.34259756367277E-02_wp, &
      & 9.05992698758370E-04_wp,-3.65628894303910E-04_wp, 6.83482385414134E-04_wp, &
      & 9.18499763897787E-04_wp, 2.19806549831237E-03_wp,-3.03607735046960E-03_wp, &
      & 9.93178929004829E-03_wp, 1.89070265565559E-02_wp, 1.33699745251597E-02_wp, &
      &-1.08352826422085E-03_wp,-1.96526810804084E-03_wp,-1.25095952083955E-03_wp, &
      &-2.17685295256138E-05_wp, 1.37310447814570E-03_wp, 9.79754791069693E-01_wp, &
      &-2.23989504876099E-02_wp, 3.40636152248508E-02_wp, 8.18979621687480E-02_wp, &
      & 2.09012487006980E-03_wp, 8.40503347739840E-04_wp, 1.56038857676967E-03_wp, &
      &-3.26760910351179E-03_wp,-3.45211447044840E-03_wp, 9.71328579482040E-03_wp, &
      & 1.51505155319696E-02_wp, 4.54245265618630E-03_wp, 1.35788864585981E-02_wp, &
      & 7.65937840551499E-02_wp, 3.75844989308498E-02_wp,-5.57175718616637E-03_wp, &
      & 3.28723251997892E-02_wp, 7.33216018913339E-02_wp, 4.00823788730330E-02_wp, &
      & 2.57126867852129E-02_wp, 1.41827982064466E-02_wp,-2.48692793280075E-02_wp, &
      & 2.39815493206550E-02_wp, 4.82654826897998E-02_wp, 5.49357483690815E-02_wp, &
      &-2.18177810703358E-02_wp,-4.44215180079381E-02_wp,-1.66229817006137E-02_wp, &
      & 9.67852731639095E-03_wp,-1.86454186027926E-03_wp,-3.34642237401371E-03_wp, &
      & 1.43421135623885E-03_wp, 2.49064579616363E-04_wp, 5.36800952916668E-03_wp, &
      & 1.59310863120097E-02_wp,-2.68014119423897E-03_wp,-2.18204966861355E-02_wp, &
      & 1.12806043159290E-02_wp,-1.82074210885778E-03_wp, 3.13903002693610E-03_wp, &
      & 2.27846314199926E-04_wp,-2.21004109845814E-03_wp,-3.85050874145553E-04_wp, &
      &-2.23989504876099E-02_wp, 9.07216790441305E-01_wp, 4.71099439338643E-03_wp, &
      & 6.98669050585301E-03_wp,-2.76452181520743E-02_wp,-1.32617256121615E-02_wp, &
      &-5.84858535888869E-03_wp,-9.24621410537155E-04_wp,-8.64305835007620E-03_wp, &
      &-1.32264514235801E-02_wp, 4.13303227393656E-03_wp, 1.19055309005404E-02_wp, &
      &-2.03601998000497E-02_wp, 3.28686383455122E-02_wp,-1.69469206158735E-02_wp, &
      & 8.67153633545511E-02_wp, 4.78972976104863E-02_wp,-7.19861577005277E-02_wp, &
      & 5.07256361602297E-02_wp, 1.14966740661808E-02_wp,-6.62740370090869E-02_wp, &
      & 6.56157086454086E-03_wp, 2.95058137726838E-02_wp,-1.14992339499985E-02_wp, &
      &-2.23923554737875E-02_wp,-4.85160230143877E-03_wp,-4.45142330258693E-03_wp, &
      & 2.23035839447183E-02_wp,-1.30717775429154E-02_wp,-1.56221941949092E-03_wp, &
      &-3.20723614997892E-04_wp, 1.29954197473103E-03_wp,-3.01515131459455E-04_wp, &
      &-1.11865135945953E-04_wp, 1.87885921811317E-02_wp,-3.97830522407406E-02_wp, &
      &-1.64116449632517E-02_wp,-9.44493581999449E-03_wp,-1.01687330466814E-03_wp, &
      & 3.85187428687074E-03_wp, 1.88274865285337E-04_wp,-1.02849359838044E-03_wp, &
      &-4.29155366733989E-03_wp, 3.40636152248508E-02_wp, 4.71099439338643E-03_wp, &
      & 9.03653999747146E-01_wp,-1.07477493790883E-02_wp,-6.79035487613655E-04_wp, &
      & 7.65027028289048E-03_wp,-1.49468871725830E-02_wp,-2.63446016222323E-02_wp, &
      &-1.04417300489744E-03_wp,-4.00506954181809E-02_wp, 1.19937863231962E-02_wp, &
      &-1.89552860618799E-02_wp,-2.27788852327071E-02_wp, 6.44067474820589E-02_wp, &
      & 4.36814617081137E-02_wp, 2.90410287474541E-03_wp,-1.04712017981930E-01_wp, &
      &-4.69373915752390E-02_wp, 6.01728319806439E-02_wp, 3.88665260793635E-02_wp, &
      & 2.39879765891764E-02_wp, 6.03966283304102E-02_wp, 2.61728970650068E-02_wp, &
      &-6.93910045606182E-02_wp,-3.95054115591359E-02_wp, 2.44761467322725E-03_wp, &
      & 3.22808677054596E-02_wp,-9.27699182867206E-03_wp,-3.47522554508627E-03_wp, &
      &-8.01779151306951E-04_wp,-3.02716907724124E-04_wp,-6.84034941761651E-04_wp, &
      &-1.87144089307596E-03_wp,-2.11155401443314E-03_wp, 5.08467535689675E-03_wp, &
      &-1.26975405205949E-02_wp,-3.50049447182234E-02_wp,-1.90246054656196E-03_wp, &
      & 1.84168810430206E-03_wp, 2.11734626341689E-03_wp, 2.05809775538211E-03_wp, &
      &-9.02442359029762E-04_wp,-8.98951630102113E-04_wp, 8.18979621687480E-02_wp, &
      & 6.98669050585301E-03_wp,-1.07477493790883E-02_wp, 8.78220908101580E-01_wp, &
      & 8.78198962890330E-03_wp, 7.29973779501288E-04_wp, 1.48135170964987E-02_wp, &
      &-1.41265088250698E-02_wp,-2.63165868327967E-02_wp,-7.30975427202173E-04_wp, &
      &-4.00974483671832E-04_wp,-4.21614601950282E-04_wp,-6.78210598851328E-04_wp, &
      &-1.00233070467192E-03_wp,-7.35276572441748E-04_wp, 3.78053496299335E-04_wp, &
      &-1.59454356614833E-03_wp,-2.10917504939904E-03_wp,-1.08074030233043E-03_wp, &
      &-9.20983705261979E-04_wp, 1.05566543858823E-03_wp, 1.66757073509143E-03_wp, &
      & 5.32051158914855E-04_wp,-1.97145781885766E-03_wp,-2.59058745651069E-03_wp, &
      & 2.10650456051350E-03_wp, 5.51511024939793E-03_wp, 4.79918109486991E-04_wp, &
      &-2.41917211488525E-03_wp, 5.29754230191627E-05_wp, 1.66642160459605E-04_wp, &
      &-6.77514813138002E-05_wp,-4.15809966344438E-05_wp,-3.09423947555407E-04_wp, &
      &-1.69023816854923E-03_wp, 1.43642412751397E-03_wp, 1.92048599024827E-03_wp, &
      &-1.78419681380516E-03_wp, 1.24242288971551E-04_wp,-1.34264453636114E-04_wp, &
      & 1.14606708088033E-05_wp, 1.12455065949355E-04_wp, 7.11466472221123E-06_wp, &
      & 2.09012487006980E-03_wp,-2.76452181520743E-02_wp,-6.79035487613655E-04_wp, &
      & 8.78198962890330E-03_wp, 9.72829001385098E-04_wp, 4.31613119157574E-04_wp, &
      & 3.39329532587728E-04_wp,-1.10362577227696E-04_wp, 8.20860388583884E-06_wp, &
      &-2.86085329171831E-04_wp,-1.44661987627843E-04_wp,-8.16969303104830E-05_wp, &
      &-4.08503618156095E-04_wp,-6.88819688392179E-04_wp,-6.06690581183132E-04_wp, &
      & 9.38622543052269E-04_wp,-7.29860323332650E-05_wp,-1.45840439450576E-03_wp, &
      &-5.41519705143916E-04_wp,-1.88651521176194E-03_wp, 3.60304872458680E-03_wp, &
      & 4.66111281358563E-04_wp, 4.14351970344133E-03_wp, 3.56893171984317E-05_wp, &
      &-2.26084037291830E-03_wp, 5.13537734107091E-04_wp, 2.23194762526908E-03_wp, &
      & 1.40000483980811E-04_wp,-1.02717741612734E-03_wp,-1.67720095156313E-05_wp, &
      & 1.29518345211162E-04_wp, 3.57114046695395E-05_wp, 1.50166180255650E-05_wp, &
      &-1.21929197698885E-04_wp,-6.52194452940651E-04_wp,-5.15623780968630E-04_wp, &
      & 2.73294147837168E-04_wp,-9.11286541924727E-04_wp, 2.60909585424519E-05_wp, &
      &-2.51838077169850E-05_wp,-2.53947944807633E-05_wp, 6.84526285750054E-05_wp, &
      &-7.22055400537447E-05_wp, 8.40503347739840E-04_wp,-1.32617256121615E-02_wp, &
      & 7.65027028289048E-03_wp, 7.29973779501288E-04_wp, 4.31613119157574E-04_wp, &
      & 3.07330023331940E-04_wp,-6.96308160686061E-05_wp,-2.72899832828417E-04_wp, &
      & 1.35850391213342E-04_wp,-5.83829373776639E-04_wp,-2.20905355380273E-05_wp, &
      &-5.05453807147833E-04_wp,-1.59367877771540E-04_wp, 3.53857312562318E-04_wp, &
      & 6.06034207881716E-04_wp,-1.22901953593348E-03_wp,-2.39582887844565E-03_wp, &
      & 1.08123286916086E-04_wp, 1.74128670059436E-04_wp, 2.10285595586614E-03_wp, &
      &-3.54139161533717E-03_wp, 1.29891259804992E-03_wp,-5.00522584021680E-03_wp, &
      &-2.34151302881454E-03_wp, 4.94952180334675E-04_wp, 8.26572953275499E-04_wp, &
      & 2.14969208354250E-03_wp,-1.14187690596524E-03_wp,-7.37519467381799E-04_wp, &
      & 6.88550780578720E-05_wp,-6.06645445272479E-05_wp,-1.37991308366822E-04_wp, &
      &-9.17945277861075E-05_wp,-1.22852682598501E-04_wp,-1.24575925065985E-03_wp, &
      & 1.97029419195767E-03_wp, 2.40805559719466E-04_wp,-5.17440346509960E-04_wp, &
      & 9.82423307899891E-05_wp,-5.85408576231827E-05_wp, 7.94171880800760E-05_wp, &
      &-1.17897739759140E-05_wp, 1.11585039120817E-04_wp, 1.56038857676967E-03_wp, &
      &-5.84858535888869E-03_wp,-1.49468871725830E-02_wp, 1.48135170964987E-02_wp, &
      & 3.39329532587728E-04_wp,-6.96308160686061E-05_wp, 5.91760718215047E-04_wp, &
      & 2.52912255874616E-04_wp,-4.23264327989667E-04_wp, 1.06473924707411E-03_wp, &
      &-4.45102163842455E-04_wp, 3.80106625809516E-05_wp, 1.00614012061061E-03_wp, &
      &-1.40694498187750E-03_wp, 1.10660155867559E-04_wp,-2.36862461809407E-03_wp, &
      & 4.78040749239259E-04_wp, 2.31640389940747E-03_wp,-2.34898043484794E-03_wp, &
      & 5.97064769855728E-04_wp,-2.95929713103497E-03_wp,-8.75032383977239E-04_wp, &
      &-6.01050825635358E-03_wp, 4.72722357649428E-04_wp, 2.52712586788310E-03_wp, &
      & 7.14219643150120E-04_wp, 2.26334563100816E-04_wp, 2.76159075773128E-04_wp, &
      & 2.37177085329030E-03_wp, 1.49459436391235E-04_wp,-3.58157141419361E-05_wp, &
      &-1.10652382975503E-04_wp, 3.21446802815042E-05_wp, 4.05393279044262E-05_wp, &
      &-1.84891013913080E-03_wp, 3.39979714316987E-03_wp, 3.98704788703455E-03_wp, &
      & 2.69307620534643E-03_wp,-1.33251711439347E-05_wp,-2.31098244614289E-04_wp, &
      &-5.67477088132154E-05_wp, 2.91415636443471E-05_wp, 2.53643008087482E-04_wp, &
      &-3.26760910351179E-03_wp,-9.24621410537155E-04_wp,-2.63446016222323E-02_wp, &
      &-1.41265088250698E-02_wp,-1.10362577227696E-04_wp,-2.72899832828417E-04_wp, &
      & 2.52912255874616E-04_wp, 1.08310028976880E-03_wp, 4.38482863598217E-04_wp, &
      & 1.20864032142143E-03_wp,-5.21718641131263E-04_wp, 4.70419178273609E-04_wp, &
      & 6.42501172500220E-04_wp,-2.36540545013463E-03_wp,-1.13219077295460E-03_wp, &
      & 3.42392873702847E-05_wp, 2.59547000701250E-03_wp, 5.16614467010584E-04_wp, &
      &-2.72927683099058E-03_wp,-3.32760114530298E-03_wp, 4.28984898099735E-03_wp, &
      &-1.61471980448458E-03_wp, 3.53071480903333E-03_wp, 2.70282185999502E-03_wp, &
      &-5.02041178157069E-04_wp, 1.14929051465742E-03_wp,-2.69348432554597E-04_wp, &
      & 1.54213587967972E-03_wp, 1.59044196104693E-03_wp, 3.48745958377767E-05_wp, &
      & 1.47884993198466E-04_wp, 8.18569440148499E-05_wp, 1.29877997002024E-04_wp, &
      & 4.47136744943953E-05_wp,-3.95093353419281E-04_wp, 1.57372281668683E-04_wp, &
      & 3.83007723351944E-03_wp, 1.58235914556952E-03_wp,-7.80121783922660E-05_wp, &
      &-1.74472472563144E-04_wp,-1.57956030126034E-04_wp, 1.16342256034267E-04_wp, &
      & 2.63123070323647E-05_wp,-3.45211447044840E-03_wp,-8.64305835007620E-03_wp, &
      &-1.04417300489744E-03_wp,-2.63165868327967E-02_wp, 8.20860388583884E-06_wp, &
      & 1.35850391213342E-04_wp,-4.23264327989667E-04_wp, 4.38482863598217E-04_wp, &
      & 9.40813877255382E-04_wp, 7.11952239175541E-03_wp,-4.48759433563640E-04_wp, &
      & 7.04758808977037E-04_wp,-2.94624945675146E-04_wp,-2.74166581356118E-03_wp, &
      & 6.84153465016012E-03_wp,-5.31279761822284E-03_wp, 4.85658416213751E-03_wp, &
      & 9.50019121811072E-04_wp,-2.43299151738347E-03_wp,-9.06679526463687E-03_wp, &
      &-3.80607531989747E-04_wp,-2.37874651841234E-03_wp,-3.06390494210974E-04_wp, &
      & 3.36294433121599E-03_wp, 8.99503672589565E-03_wp,-1.11688855965280E-03_wp, &
      & 3.20046292353626E-02_wp, 2.32606695048069E-02_wp, 2.10966188523550E-02_wp, &
      & 1.14006433275942E-03_wp, 1.14147229303161E-03_wp,-3.02272998088171E-04_wp, &
      & 7.61032562708509E-04_wp,-4.79144450852202E-04_wp,-1.30082860861911E-03_wp, &
      &-3.51696354973918E-02_wp,-1.39407590739072E-02_wp, 2.24022500427104E-02_wp, &
      &-1.30001097214204E-03_wp, 9.03676498870888E-04_wp,-5.83171540616180E-04_wp, &
      &-5.86376985234874E-04_wp,-5.20172683811180E-04_wp,-1.20825923162132E-03_wp, &
      & 1.01661826271289E-02_wp,-1.38165545028930E-02_wp,-4.08742499801316E-02_wp, &
      &-7.18067691148890E-04_wp,-2.63641797653816E-04_wp,-5.93966433083350E-04_wp, &
      & 1.02979125404260E-03_wp, 1.20885516832022E-03_wp,-4.48759433563640E-04_wp, &
      & 2.28687129134323E-03_wp, 5.41396449947113E-04_wp,-2.59430166947036E-04_wp, &
      &-2.64919936937197E-03_wp, 5.77530559995243E-04_wp,-5.26976013453336E-05_wp, &
      &-1.70360657293970E-03_wp,-2.12815610708007E-04_wp, 1.06766405478897E-03_wp, &
      & 1.23297231618397E-03_wp,-1.24027234557057E-03_wp,-2.62975652749672E-03_wp, &
      & 6.81819221360484E-04_wp,-2.77984913399787E-03_wp, 6.15223438264382E-04_wp, &
      & 7.17489050802568E-03_wp,-1.15284447884864E-02_wp,-1.66995151518811E-02_wp, &
      &-1.96285390992294E-02_wp,-8.09390428732436E-04_wp,-6.27932511020789E-04_wp, &
      & 1.21450766756593E-04_wp,-6.28894926244619E-04_wp,-7.62484937972927E-05_wp, &
      &-8.00507020681381E-03_wp,-1.18379012605945E-02_wp,-1.74855485990918E-02_wp, &
      & 2.22060923985063E-02_wp,-8.41884382085804E-04_wp, 6.54018323036504E-04_wp, &
      &-9.79735922629530E-05_wp,-7.38014920537504E-04_wp, 1.67408952638742E-04_wp, &
      & 2.55413860729610E-03_wp, 1.44077180987701E-02_wp, 4.81975145073823E-03_wp, &
      & 1.15402494798051E-02_wp,-3.91472860045293E-04_wp,-1.85447400076972E-04_wp, &
      & 2.69372961745411E-05_wp,-3.83960703078049E-04_wp,-5.52056921452457E-04_wp, &
      & 7.04758808977037E-04_wp, 5.41396449947113E-04_wp, 1.62716519681099E-03_wp, &
      & 3.61628176171863E-04_wp,-1.63409127414158E-03_wp,-1.07250498654332E-03_wp, &
      & 3.31847658399270E-04_wp,-6.83598102194624E-04_wp, 5.75709000782809E-04_wp, &
      &-1.94547758751706E-04_wp,-6.02151183303941E-04_wp, 1.49179954780895E-03_wp, &
      &-2.93152544359318E-03_wp, 1.04205333697131E-03_wp, 1.48481072358356E-04_wp, &
      & 6.08528200623643E-04_wp, 4.57690583399948E-03_wp,-1.88737263402908E-02_wp, &
      & 8.75664181962437E-03_wp,-1.26486079338156E-02_wp,-7.10043290946883E-04_wp, &
      & 2.89620397186540E-06_wp, 6.06131016453156E-04_wp, 2.82293389951878E-06_wp, &
      & 2.86985578991552E-04_wp,-3.80262208043235E-03_wp,-1.56986199671817E-02_wp, &
      & 9.95289435444545E-03_wp, 1.08109379474327E-02_wp,-6.58426823850392E-04_wp, &
      &-2.72580344706362E-05_wp,-5.84203711930900E-04_wp, 1.67615675300633E-05_wp, &
      &-1.95672888720847E-04_wp,-3.63428057490289E-03_wp, 5.13544987508727E-03_wp, &
      & 1.02100793239385E-02_wp,-1.81555166175954E-02_wp,-3.91264029795281E-04_wp, &
      & 9.88365504774939E-06_wp,-5.75792938733305E-04_wp,-4.71112289724211E-05_wp, &
      & 5.46257757174961E-04_wp,-2.94624945675146E-04_wp,-2.59430166947036E-04_wp, &
      & 3.61628176171863E-04_wp, 2.43622668168911E-03_wp,-1.06123973390391E-03_wp, &
      &-1.65324946777580E-03_wp,-7.25296306812267E-05_wp, 1.80304070856627E-03_wp, &
      & 2.13517632377509E-03_wp,-1.81460080347582E-03_wp, 2.37543992326638E-03_wp, &
      & 6.76725101117473E-04_wp,-1.81030095889503E-03_wp,-1.71651898481580E-03_wp, &
      & 2.34703941229104E-03_wp,-2.97757430264917E-04_wp, 4.75605525498004E-03_wp, &
      &-1.99507322125151E-02_wp,-1.13575027794183E-02_wp, 4.49366063601879E-03_wp, &
      &-2.46172687351467E-04_wp,-6.32133131700029E-04_wp, 8.27214390433664E-05_wp, &
      &-1.20209486448796E-04_wp, 7.01975281929019E-04_wp, 5.46887348152095E-03_wp, &
      & 2.07465098916424E-02_wp, 1.12870324819581E-02_wp, 2.62360032992966E-03_wp, &
      & 3.30170019630296E-04_wp,-6.85792113938563E-04_wp, 6.13175091245107E-05_wp, &
      & 2.02457181360983E-04_wp, 6.65644753860528E-04_wp,-9.08181098067992E-03_wp, &
      & 1.33513181820185E-02_wp,-1.99268786228131E-02_wp,-2.34217202320817E-02_wp, &
      &-6.69813653944916E-04_wp,-4.38712775661301E-04_wp,-1.20006019567996E-04_wp, &
      & 1.02053886193382E-03_wp, 5.94007597988431E-04_wp,-2.74166581356118E-03_wp, &
      &-2.64919936937197E-03_wp,-1.63409127414158E-03_wp,-1.06123973390391E-03_wp, &
      & 4.26660636479625E-02_wp, 5.48867639267580E-03_wp, 2.54176254686327E-03_wp, &
      & 7.36832296027571E-03_wp,-4.03975654473346E-04_wp,-3.26942203748600E-04_wp, &
      & 3.91589546659811E-03_wp, 4.64342996080713E-03_wp, 1.33992841815513E-02_wp, &
      & 7.49785230432728E-04_wp, 8.17240065738855E-04_wp,-3.19238355445102E-04_wp, &
      & 2.23646244137191E-02_wp, 6.59548065857688E-02_wp, 8.68754738499332E-02_wp, &
      & 1.25040456174465E-03_wp, 8.43040832431947E-04_wp, 2.57139294590380E-03_wp, &
      & 4.69721320177037E-04_wp, 1.34881283254919E-03_wp,-1.20147515573047E-03_wp, &
      &-2.39915293849376E-02_wp, 7.91814684776186E-02_wp, 6.88501216429309E-02_wp, &
      &-1.27801136787665E-02_wp, 1.67318715111030E-03_wp,-2.26525736225413E-03_wp, &
      & 1.03312589059841E-04_wp, 1.25511859423494E-03_wp, 1.23392237280277E-03_wp, &
      &-1.35620731019103E-02_wp, 7.56618991833905E-02_wp, 3.15084854540213E-02_wp, &
      & 6.45053437127255E-02_wp,-9.38359904068336E-04_wp,-6.34043613156057E-04_wp, &
      & 3.14862311605364E-04_wp,-1.37787803489053E-03_wp,-2.24225399100506E-03_wp, &
      & 6.84153465016012E-03_wp, 5.77530559995243E-04_wp,-1.07250498654332E-03_wp, &
      &-1.65324946777580E-03_wp, 5.48867639267580E-03_wp, 3.04783663710076E-02_wp, &
      & 5.04606442095931E-03_wp,-4.66371208138162E-03_wp,-7.19052804954136E-03_wp, &
      & 3.26384647666834E-03_wp,-1.16648070070522E-02_wp,-1.74839240093770E-03_wp, &
      & 1.30572623509032E-03_wp,-5.03897420550643E-03_wp,-8.96763548844018E-04_wp, &
      & 9.84759792548253E-03_wp, 2.07224740193056E-02_wp, 6.70688097195501E-02_wp, &
      &-1.96459370841132E-02_wp, 7.81809041418510E-02_wp, 2.84420597616094E-03_wp, &
      & 3.81213248123848E-04_wp,-1.93389141426929E-03_wp, 4.82420477979836E-04_wp, &
      &-5.32810176688890E-04_wp, 1.83968711176338E-02_wp,-6.46585082511956E-02_wp, &
      & 3.76360397992007E-02_wp, 6.76726646968229E-02_wp,-2.56333530435437E-03_wp, &
      & 3.57609160919946E-04_wp,-2.02312098196746E-03_wp,-3.47311396600001E-04_wp, &
      &-5.79280405069388E-04_wp,-5.21056540703202E-03_wp, 3.66511396790250E-02_wp, &
      &-1.61959501492248E-02_wp, 4.22610487435240E-02_wp,-7.23150166635905E-04_wp, &
      &-6.96300015777383E-04_wp, 6.92000783727210E-04_wp, 2.14568368478059E-04_wp, &
      &-1.19565402301583E-03_wp,-5.31279761822284E-03_wp,-5.26976013453336E-05_wp, &
      & 3.31847658399270E-04_wp,-7.25296306812267E-05_wp, 2.54176254686327E-03_wp, &
      & 5.04606442095931E-03_wp, 2.84679587556288E-02_wp, 3.54894557191667E-03_wp, &
      &-7.72046463325575E-04_wp, 3.39812061949490E-03_wp, 1.80295732067456E-03_wp, &
      &-1.72699146503989E-03_wp, 2.08068962045653E-04_wp,-9.57025303578533E-04_wp, &
      &-6.72245961471300E-04_wp,-1.25000483023922E-02_wp,-7.48190813025679E-03_wp, &
      & 1.50194621101878E-02_wp,-9.51746257217116E-02_wp, 1.07105823451813E-02_wp, &
      & 2.74276798196795E-04_wp,-2.28477526111468E-03_wp,-1.54670546994734E-03_wp, &
      &-1.52509696471363E-03_wp,-8.34604782021031E-05_wp,-9.70631584445533E-03_wp, &
      & 5.89106335179422E-04_wp, 8.54634077370166E-02_wp,-6.23156571113277E-03_wp, &
      & 2.39237960359928E-04_wp,-1.86157192226547E-03_wp,-1.22914115729581E-03_wp, &
      & 1.62380339070335E-03_wp,-2.57899109827770E-04_wp,-9.88014923090812E-03_wp, &
      &-6.25144814937098E-03_wp, 8.60612099656798E-02_wp, 3.14099372304354E-03_wp, &
      & 3.78532516914578E-04_wp, 9.38872653962262E-04_wp,-1.19990637124721E-03_wp, &
      &-2.29785457537120E-03_wp, 4.74724238068227E-05_wp, 4.85658416213751E-03_wp, &
      &-1.70360657293970E-03_wp,-6.83598102194624E-04_wp, 1.80304070856627E-03_wp, &
      & 7.36832296027571E-03_wp,-4.66371208138162E-03_wp, 3.54894557191667E-03_wp, &
      & 3.34393781938834E-02_wp, 4.73561694058205E-03_wp,-9.20066557038393E-03_wp, &
      &-4.78340485526953E-03_wp,-5.05261443662532E-03_wp, 6.84791172569323E-04_wp, &
      & 2.03099825523925E-03_wp, 9.45917682510677E-03_wp, 3.78879388099240E-03_wp, &
      & 1.40822121394520E-02_wp, 7.85903417261112E-02_wp,-1.34573967461192E-02_wp, &
      & 5.53985700081892E-03_wp, 1.51665855913107E-03_wp, 4.92005695112074E-04_wp, &
      &-1.33814279600227E-03_wp, 9.60493988868284E-06_wp,-1.49543254126587E-03_wp, &
      &-1.19162210963841E-02_wp, 7.77031024200781E-02_wp,-3.00536616447836E-02_wp, &
      &-2.00732274112082E-03_wp, 1.10329669950774E-03_wp,-1.68670894927570E-04_wp, &
      & 1.46237940335649E-03_wp,-6.46986175422616E-04_wp, 1.68257298559226E-03_wp, &
      & 2.13600011676250E-02_wp, 3.25281161052200E-02_wp, 4.67381503378004E-02_wp, &
      &-1.03858523700384E-01_wp,-1.50740974288056E-03_wp, 3.82651994861127E-05_wp, &
      &-2.40313487407472E-03_wp, 3.71631530998101E-04_wp, 2.60157115829106E-03_wp, &
      & 9.50019121811072E-04_wp,-2.12815610708007E-04_wp, 5.75709000782809E-04_wp, &
      & 2.13517632377509E-03_wp,-4.03975654473346E-04_wp,-7.19052804954136E-03_wp, &
      &-7.72046463325575E-04_wp, 4.73561694058205E-03_wp, 4.18770813681395E-02_wp, &
      &-2.50290440454828E-03_wp, 7.65727418568397E-04_wp, 2.12624345393224E-03_wp, &
      &-4.97027877918030E-03_wp,-5.83472669844520E-03_wp, 5.51732011957694E-03_wp, &
      &-4.27826077944579E-04_wp,-8.82403162605848E-03_wp, 3.37743473393132E-02_wp, &
      &-3.36885096134573E-02_wp,-8.78793958800221E-02_wp,-1.20767003906758E-03_wp, &
      &-4.48556039119347E-04_wp,-1.81086224544769E-04_wp,-1.30225149564277E-03_wp, &
      &-1.72564923333505E-03_wp,-8.62878718054140E-03_wp,-3.71725125852140E-02_wp, &
      & 3.53696185434452E-02_wp,-9.55348684789663E-02_wp, 1.04090454647008E-03_wp, &
      &-3.59221400041849E-04_wp,-3.57354903215663E-04_wp, 1.79429741173201E-03_wp, &
      &-2.19262341202439E-03_wp, 2.19832369216965E-02_wp, 7.20780368313875E-02_wp, &
      &-7.10320337351810E-02_wp,-4.72446672147925E-02_wp,-2.04798561893996E-03_wp, &
      &-1.49699282601155E-03_wp, 1.86593826628628E-04_wp, 2.32423304833064E-03_wp, &
      & 4.22160666294196E-04_wp,-2.43299151738347E-03_wp, 1.06766405478897E-03_wp, &
      &-1.94547758751706E-04_wp,-1.81460080347582E-03_wp,-3.26942203748600E-04_wp, &
      & 3.26384647666834E-03_wp, 3.39812061949490E-03_wp,-9.20066557038393E-03_wp, &
      &-2.50290440454828E-03_wp, 1.18831170923295E-02_wp, 3.11417561356810E-03_wp, &
      &-1.79633667273998E-03_wp, 1.08428680718719E-03_wp, 2.95504975538360E-03_wp, &
      & 9.40583767874230E-04_wp,-8.18360385983953E-04_wp, 2.70527565876776E-05_wp, &
      &-3.42955794705007E-02_wp,-1.00765000114553E-02_wp, 3.53152003767956E-02_wp, &
      & 3.49311376473168E-04_wp,-7.42699380839901E-04_wp, 4.39431782533191E-06_wp, &
      & 3.27030474840371E-04_wp, 1.58749721679993E-03_wp, 5.59939284194209E-04_wp, &
      &-4.33600687024431E-02_wp,-1.57001763492994E-02_wp,-3.92643295028350E-02_wp, &
      & 4.03441445098017E-04_wp, 1.24040844923198E-03_wp, 2.73735765525016E-04_wp, &
      & 2.82584056419365E-04_wp,-2.02260087963467E-03_wp,-8.83664108598077E-04_wp, &
      & 3.13114544459522E-02_wp, 3.67817089890598E-02_wp, 5.10642880072331E-02_wp, &
      &-6.93371177426215E-04_wp,-1.52128659470300E-04_wp,-6.60971351199133E-05_wp, &
      &-2.04580162271769E-03_wp,-1.86857557022708E-03_wp,-9.06679526463687E-03_wp, &
      & 1.23297231618397E-03_wp,-6.02151183303941E-04_wp, 2.37543992326638E-03_wp, &
      & 3.91589546659811E-03_wp,-1.16648070070522E-02_wp, 1.80295732067456E-03_wp, &
      &-4.78340485526953E-03_wp, 7.65727418568397E-04_wp, 3.11417561356810E-03_wp, &
      & 1.53506396343196E-02_wp, 4.93458019576049E-04_wp, 1.16823409155536E-03_wp, &
      & 3.13052003010560E-04_wp,-1.86532900379958E-03_wp,-8.96110846963001E-03_wp, &
      &-9.72183490348629E-04_wp,-6.96677285675238E-02_wp,-2.42640812741514E-02_wp, &
      &-2.43105967018234E-02_wp,-1.88909058022195E-03_wp,-1.84620872265765E-03_wp, &
      & 8.40605401675080E-04_wp,-8.31609991216725E-04_wp, 1.57654392909966E-03_wp, &
      &-6.32238839106823E-04_wp, 6.52028457208073E-02_wp, 1.09811119079166E-04_wp, &
      &-2.52064987261363E-02_wp, 1.91383514382521E-03_wp,-1.05856014409741E-03_wp, &
      & 1.27412333131700E-03_wp, 3.48350771197910E-04_wp, 1.31355340280529E-03_wp, &
      & 5.25854327632037E-05_wp, 2.45048495855122E-02_wp,-2.56655009577307E-03_wp, &
      & 4.06754362288397E-02_wp,-4.24473085272663E-04_wp,-4.73872134127600E-04_wp, &
      & 6.75497357337969E-04_wp,-5.61485634191755E-04_wp,-1.60968883271744E-03_wp, &
      &-3.80607531989747E-04_wp,-1.24027234557057E-03_wp, 1.49179954780895E-03_wp, &
      & 6.76725101117473E-04_wp, 4.64342996080713E-03_wp,-1.74839240093770E-03_wp, &
      &-1.72699146503989E-03_wp,-5.05261443662532E-03_wp, 2.12624345393224E-03_wp, &
      &-1.79633667273998E-03_wp, 4.93458019576049E-04_wp, 5.57899823592994E-03_wp, &
      &-1.29897186515243E-04_wp,-7.80032096532300E-04_wp,-1.84892158652279E-04_wp, &
      &-9.62092557670708E-04_wp, 3.41667470649450E-05_wp,-2.14598915207260E-02_wp, &
      & 4.22313704080252E-02_wp,-1.51574111619879E-02_wp,-7.82803317697228E-04_wp, &
      & 8.27111347824794E-04_wp, 1.17930469715137E-03_wp, 5.54949096390952E-04_wp, &
      & 3.31242618197280E-04_wp,-2.19959469457463E-04_wp,-3.87427809377845E-03_wp, &
      & 4.92359767360659E-02_wp, 3.76026077550443E-04_wp,-1.50246094084492E-04_wp, &
      &-1.27777504135482E-03_wp,-1.02227177755345E-03_wp, 9.55873747394603E-04_wp, &
      &-4.76745213909973E-05_wp, 5.49263377585196E-04_wp,-6.37450174955805E-03_wp, &
      &-1.46067658017428E-02_wp,-4.95983018633171E-03_wp, 1.86494551430705E-04_wp, &
      &-2.30162533491128E-05_wp, 1.12528083885423E-04_wp, 5.83275218807406E-04_wp, &
      & 4.19280315359715E-04_wp,-2.37874651841234E-03_wp,-2.62975652749672E-03_wp, &
      &-2.93152544359318E-03_wp,-1.81030095889503E-03_wp, 1.33992841815513E-02_wp, &
      & 1.30572623509032E-03_wp, 2.08068962045653E-04_wp, 6.84791172569323E-04_wp, &
      &-4.97027877918030E-03_wp, 1.08428680718719E-03_wp, 1.16823409155536E-03_wp, &
      &-1.29897186515243E-04_wp, 9.94428221843365E-03_wp,-7.19846587009332E-05_wp, &
      &-1.93960408498675E-04_wp,-2.74856880867014E-03_wp, 1.27031368525090E-03_wp, &
      & 4.18745285287946E-02_wp, 2.72966013888811E-02_wp, 2.85915868985803E-02_wp, &
      & 1.52421575467949E-03_wp, 1.17133429597765E-03_wp,-5.70259377186291E-04_wp, &
      & 8.00501008899362E-04_wp,-5.93648368390723E-04_wp,-8.89167773593586E-04_wp, &
      & 4.01388520677411E-02_wp, 6.11538319289382E-03_wp,-3.43827040288937E-02_wp, &
      & 1.81283369635852E-03_wp,-6.03244529432564E-04_wp, 9.60330421289763E-04_wp, &
      & 6.66961474578803E-04_wp, 1.38916670200151E-04_wp,-8.40001714014691E-04_wp, &
      &-1.91461291146198E-02_wp, 5.62002262438386E-03_wp, 4.94693462764905E-02_wp, &
      & 1.30712345342734E-03_wp, 4.10116747608054E-04_wp, 9.58776573157139E-04_wp, &
      &-7.52025031398366E-04_wp,-1.22301380113740E-03_wp,-3.06390494210974E-04_wp, &
      & 6.81819221360484E-04_wp, 1.04205333697131E-03_wp,-1.71651898481580E-03_wp, &
      & 7.49785230432728E-04_wp,-5.03897420550643E-03_wp,-9.57025303578533E-04_wp, &
      & 2.03099825523925E-03_wp,-5.83472669844520E-03_wp, 2.95504975538360E-03_wp, &
      & 3.13052003010560E-04_wp,-7.80032096532300E-04_wp,-7.19846587009332E-05_wp, &
      & 6.02218778760986E-03_wp,-1.96369388860250E-04_wp, 4.00523412010053E-04_wp, &
      & 2.43415200592200E-05_wp,-1.51845049947363E-02_wp, 2.87327242622429E-02_wp, &
      &-9.49249765705375E-03_wp,-5.48713980681319E-04_wp, 5.54022846553796E-04_wp, &
      & 8.01687607586127E-04_wp, 3.77777202539428E-04_wp, 1.95523826272871E-04_wp, &
      & 6.68905535755028E-04_wp,-5.70282419272425E-03_wp,-3.25394854692952E-02_wp, &
      &-4.86038808135031E-03_wp,-6.16005794354335E-05_wp, 1.02574808923692E-03_wp, &
      & 4.82756606241590E-04_wp,-5.53868614428592E-04_wp,-4.61704646589813E-04_wp, &
      &-4.42109441645743E-04_wp, 1.09535528037961E-03_wp, 5.78076943435148E-02_wp, &
      &-1.68594875769831E-03_wp,-1.50326172865105E-04_wp, 5.22664015461809E-04_wp, &
      &-1.09131629766252E-03_wp,-1.82188558275370E-03_wp,-3.25003079162797E-05_wp, &
      & 3.36294433121599E-03_wp,-2.77984913399787E-03_wp, 1.48481072358356E-04_wp, &
      & 2.34703941229104E-03_wp, 8.17240065738855E-04_wp,-8.96763548844018E-04_wp, &
      &-6.72245961471300E-04_wp, 9.45917682510677E-03_wp, 5.51732011957694E-03_wp, &
      & 9.40583767874230E-04_wp,-1.86532900379958E-03_wp,-1.84892158652279E-04_wp, &
      &-1.93960408498675E-04_wp,-1.96369388860250E-04_wp, 1.14272521505618E-02_wp, &
      & 3.81891246340949E-03_wp, 4.46942267787996E-04_wp,-4.11175285560145E-03_wp, &
      & 9.67588448002621E-03_wp, 5.64940155335164E-02_wp, 1.60503540214095E-03_wp, &
      & 1.52655440518585E-04_wp,-3.37131360832134E-04_wp, 1.19451943984326E-03_wp, &
      & 1.48935194145954E-03_wp,-9.31069490873247E-04_wp, 7.87169218444709E-03_wp, &
      & 3.89347577592347E-03_wp,-4.72243223724468E-02_wp, 1.35155308784679E-03_wp, &
      &-2.00017008301805E-04_wp, 5.23748289028198E-04_wp, 8.65609593211044E-04_wp, &
      &-8.12160921508423E-04_wp, 1.03952145445938E-03_wp, 3.42404939671757E-02_wp, &
      &-2.44882129383753E-03_wp,-6.15913706401558E-02_wp,-1.62689399094736E-03_wp, &
      &-5.53268037298817E-04_wp,-1.20341681703315E-03_wp, 1.03007344334789E-03_wp, &
      & 1.53200056222922E-03_wp, 8.99503672589565E-03_wp, 6.15223438264382E-04_wp, &
      & 6.08528200623643E-04_wp,-2.97757430264917E-04_wp,-3.19238355445102E-04_wp, &
      & 9.84759792548253E-03_wp,-1.25000483023922E-02_wp, 3.78879388099240E-03_wp, &
      &-4.27826077944579E-04_wp,-8.18360385983953E-04_wp,-8.96110846963001E-03_wp, &
      &-9.62092557670708E-04_wp,-2.74856880867014E-03_wp, 4.00523412010053E-04_wp, &
      & 3.81891246340949E-03_wp, 1.48383546641339E-02_wp, 2.14569885078231E-04_wp, &
      & 2.60228009071645E-02_wp, 4.55480551778791E-02_wp, 3.11414586934253E-02_wp, &
      & 1.24763488181786E-03_wp, 1.53320635726685E-03_wp, 7.48314601602185E-08_wp, &
      & 1.22592371170185E-03_wp,-1.03037430243343E-04_wp, 6.85864626230428E-04_wp, &
      &-3.93198699151506E-02_wp,-4.33959140311928E-02_wp, 3.84544626969289E-02_wp, &
      &-1.77253159605633E-03_wp, 1.76953567734680E-03_wp,-2.92690161508324E-04_wp, &
      &-1.42187034413839E-03_wp,-3.56619680821674E-04_wp, 2.40006331976448E-04_wp, &
      & 4.74589429392366E-02_wp,-2.86025601155137E-02_wp,-2.74869037414883E-02_wp, &
      &-1.81297319011831E-03_wp,-9.90759348859535E-04_wp,-3.97806973349569E-04_wp, &
      & 1.22820555419500E-03_wp, 4.56941288207854E-04_wp,-1.11688855965280E-03_wp, &
      & 7.17489050802568E-03_wp, 4.57690583399948E-03_wp, 4.75605525498004E-03_wp, &
      & 2.23646244137191E-02_wp, 2.07224740193056E-02_wp,-7.48190813025679E-03_wp, &
      & 1.40822121394520E-02_wp,-8.82403162605848E-03_wp, 2.70527565876776E-05_wp, &
      &-9.72183490348629E-04_wp, 3.41667470649450E-05_wp, 1.27031368525090E-03_wp, &
      & 2.43415200592200E-05_wp, 4.46942267787996E-04_wp, 2.14569885078231E-04_wp, &
      & 9.79953185022657E-01_wp,-6.56913535165536E-02_wp,-4.15746208633436E-02_wp, &
      &-4.42906788242026E-02_wp,-3.26497777863993E-03_wp,-3.02957050163202E-03_wp, &
      & 1.03160173383662E-03_wp,-2.06224126752325E-03_wp, 1.29590186401136E-03_wp, &
      &-1.73578590097926E-03_wp,-4.78755981358135E-03_wp,-4.80964039210276E-03_wp, &
      &-2.26550292609627E-02_wp, 1.03244385712790E-03_wp, 8.28645456652537E-04_wp, &
      & 8.36183698869665E-04_wp, 3.56016852022877E-04_wp,-2.26721397739007E-03_wp, &
      &-4.67936049592108E-04_wp,-2.15316812457990E-02_wp,-5.04588128389885E-03_wp, &
      & 2.25098483451193E-03_wp, 2.05812072051282E-03_wp, 5.05766248283295E-04_wp, &
      & 7.91454294590884E-04_wp, 7.16816433450587E-04_wp, 1.18902455312498E-03_wp, &
      & 3.20046292353626E-02_wp,-1.15284447884864E-02_wp,-1.88737263402908E-02_wp, &
      &-1.99507322125151E-02_wp, 6.59548065857688E-02_wp, 6.70688097195501E-02_wp, &
      & 1.50194621101878E-02_wp, 7.85903417261112E-02_wp, 3.37743473393132E-02_wp, &
      &-3.42955794705007E-02_wp,-6.96677285675238E-02_wp,-2.14598915207260E-02_wp, &
      & 4.18745285287946E-02_wp,-1.51845049947363E-02_wp,-4.11175285560145E-03_wp, &
      & 2.60228009071645E-02_wp,-6.56913535165536E-02_wp, 8.97832343104443E-01_wp, &
      &-9.01190569007232E-03_wp,-1.45144645636867E-02_wp, 1.55005079704011E-02_wp, &
      & 1.20367234816803E-02_wp,-1.34556376466504E-02_wp,-3.44948217512164E-04_wp, &
      &-2.38477695614008E-02_wp, 5.97412493903298E-03_wp, 3.58652467089654E-03_wp, &
      &-1.13563354350246E-02_wp,-9.59168996535191E-03_wp,-6.71551375241214E-04_wp, &
      & 2.30409327479312E-03_wp, 4.33495391693635E-04_wp,-1.03517926424556E-03_wp, &
      &-9.32392910000516E-04_wp,-1.75636665246333E-02_wp,-4.80394095409376E-02_wp, &
      &-4.04081526462394E-03_wp, 2.91067062790983E-02_wp, 5.41082241160368E-03_wp, &
      & 2.05070858833928E-03_wp, 2.29339909979128E-03_wp, 4.07757821937252E-04_wp, &
      &-4.27488088447255E-04_wp, 2.32606695048069E-02_wp,-1.66995151518811E-02_wp, &
      & 8.75664181962437E-03_wp,-1.13575027794183E-02_wp, 8.68754738499332E-02_wp, &
      &-1.96459370841132E-02_wp,-9.51746257217116E-02_wp,-1.34573967461192E-02_wp, &
      &-3.36885096134573E-02_wp,-1.00765000114553E-02_wp,-2.42640812741514E-02_wp, &
      & 4.22313704080252E-02_wp, 2.72966013888811E-02_wp, 2.87327242622429E-02_wp, &
      & 9.67588448002621E-03_wp, 4.55480551778791E-02_wp,-4.15746208633436E-02_wp, &
      &-9.01190569007232E-03_wp, 8.98949522434251E-01_wp,-5.85249818311784E-03_wp, &
      &-4.20396123671864E-05_wp, 2.32696964109612E-02_wp, 1.40699101827079E-02_wp, &
      & 1.57206108208992E-02_wp, 4.82043997276390E-05_wp, 1.29832452985145E-03_wp, &
      &-1.01189621165404E-02_wp, 1.64501682042580E-02_wp,-1.13646125614409E-02_wp, &
      &-6.58322439410420E-04_wp,-9.95021813203161E-05_wp,-1.69154717257233E-03_wp, &
      & 8.98861322776528E-04_wp,-2.17156562704799E-03_wp, 1.74725131414713E-03_wp, &
      &-1.23060391666262E-02_wp, 1.53478173145738E-02_wp,-5.75999068005323E-03_wp, &
      & 6.00827861785183E-04_wp, 6.76325363903052E-04_wp,-1.66014282991117E-03_wp, &
      &-2.79625792410964E-04_wp, 2.06014446932532E-03_wp, 2.10966188523550E-02_wp, &
      &-1.96285390992294E-02_wp,-1.26486079338156E-02_wp, 4.49366063601879E-03_wp, &
      & 1.25040456174465E-03_wp, 7.81809041418510E-02_wp, 1.07105823451813E-02_wp, &
      & 5.53985700081892E-03_wp,-8.78793958800221E-02_wp, 3.53152003767956E-02_wp, &
      &-2.43105967018234E-02_wp,-1.51574111619879E-02_wp, 2.85915868985803E-02_wp, &
      &-9.49249765705375E-03_wp, 5.64940155335164E-02_wp, 3.11414586934253E-02_wp, &
      &-4.42906788242026E-02_wp,-1.45144645636867E-02_wp,-5.85249818311784E-03_wp, &
      & 9.09568901510486E-01_wp, 2.36486322778065E-02_wp,-3.50001068673827E-04_wp, &
      &-9.07751827589900E-03_wp, 1.23464513582564E-02_wp, 1.66020669497395E-02_wp, &
      &-2.28312252496091E-02_wp, 1.33634035581096E-02_wp, 1.41764609986090E-03_wp, &
      &-5.45979916231744E-02_wp, 5.34400854366388E-03_wp,-4.96044569641485E-06_wp, &
      & 2.40201841492479E-03_wp, 2.48198884031917E-03_wp,-2.72121814524756E-03_wp, &
      & 1.34763066806291E-02_wp, 1.29062067504193E-02_wp,-9.80190748893888E-03_wp, &
      &-5.04907733624908E-03_wp,-2.40085902556881E-03_wp,-1.16119924081096E-03_wp, &
      &-5.81093006974993E-04_wp, 2.45027490200436E-03_wp, 1.39966952488807E-03_wp, &
      & 1.14006433275942E-03_wp,-8.09390428732436E-04_wp,-7.10043290946883E-04_wp, &
      &-2.46172687351467E-04_wp, 8.43040832431947E-04_wp, 2.84420597616094E-03_wp, &
      & 2.74276798196795E-04_wp, 1.51665855913107E-03_wp,-1.20767003906758E-03_wp, &
      & 3.49311376473168E-04_wp,-1.88909058022195E-03_wp,-7.82803317697228E-04_wp, &
      & 1.52421575467949E-03_wp,-5.48713980681319E-04_wp, 1.60503540214095E-03_wp, &
      & 1.24763488181786E-03_wp,-3.26497777863993E-03_wp, 1.55005079704011E-02_wp, &
      &-4.20396123671864E-05_wp, 2.36486322778065E-02_wp, 9.18511119490182E-04_wp, &
      & 2.19078907284031E-04_wp,-4.80302880422083E-04_wp, 3.36883951800672E-04_wp, &
      & 1.96443072180799E-05_wp,-9.42305899810999E-04_wp,-2.89537305206217E-04_wp, &
      &-1.28022878279255E-03_wp,-5.00721188868593E-03_wp, 1.98157365462255E-04_wp, &
      & 6.86048439806152E-05_wp, 1.12897189556701E-04_wp, 8.04397669199020E-05_wp, &
      &-1.60993997929043E-04_wp, 9.24348836836794E-04_wp,-1.64629006596400E-03_wp, &
      &-1.88908493067740E-03_wp,-5.38795740743565E-04_wp, 5.99810932075002E-05_wp, &
      & 1.05424830110665E-05_wp, 4.34016120608695E-05_wp, 1.16049685662807E-04_wp, &
      & 5.66604383383426E-05_wp, 1.14147229303161E-03_wp,-6.27932511020789E-04_wp, &
      & 2.89620397186540E-06_wp,-6.32133131700029E-04_wp, 2.57139294590380E-03_wp, &
      & 3.81213248123848E-04_wp,-2.28477526111468E-03_wp, 4.92005695112074E-04_wp, &
      &-4.48556039119347E-04_wp,-7.42699380839901E-04_wp,-1.84620872265765E-03_wp, &
      & 8.27111347824794E-04_wp, 1.17133429597765E-03_wp, 5.54022846553796E-04_wp, &
      & 1.52655440518585E-04_wp, 1.53320635726685E-03_wp,-3.02957050163202E-03_wp, &
      & 1.20367234816803E-02_wp, 2.32696964109612E-02_wp,-3.50001068673827E-04_wp, &
      & 2.19078907284031E-04_wp, 7.86474390598049E-04_wp, 1.81027823094485E-04_wp, &
      & 4.09324704116455E-04_wp,-3.38483160385626E-04_wp, 7.13341494118537E-04_wp, &
      &-2.49338759770228E-03_wp, 1.68735866975632E-04_wp,-7.41309576023266E-04_wp, &
      &-5.88299190152080E-05_wp, 5.59832212371451E-05_wp,-6.75477375787617E-05_wp, &
      & 1.76396292110500E-05_wp,-1.27813548615544E-04_wp,-3.12515170166504E-04_wp, &
      &-3.56717280611389E-03_wp, 3.00321001380816E-04_wp,-6.77069550327497E-04_wp, &
      & 1.50498405466455E-04_wp, 8.35670166748990E-05_wp,-1.78507090937967E-05_wp, &
      & 5.03854201102447E-06_wp, 9.88365521261211E-05_wp,-3.02272998088171E-04_wp, &
      & 1.21450766756593E-04_wp, 6.06131016453156E-04_wp, 8.27214390433664E-05_wp, &
      & 4.69721320177037E-04_wp,-1.93389141426929E-03_wp,-1.54670546994734E-03_wp, &
      &-1.33814279600227E-03_wp,-1.81086224544769E-04_wp, 4.39431782533191E-06_wp, &
      & 8.40605401675080E-04_wp, 1.17930469715137E-03_wp,-5.70259377186291E-04_wp, &
      & 8.01687607586127E-04_wp,-3.37131360832134E-04_wp, 7.48314601602185E-08_wp, &
      & 1.03160173383662E-03_wp,-1.34556376466504E-02_wp, 1.40699101827079E-02_wp, &
      &-9.07751827589900E-03_wp,-4.80302880422083E-04_wp, 1.81027823094485E-04_wp, &
      & 5.18247582373672E-04_wp, 1.22109199673396E-04_wp, 1.93563736986049E-04_wp, &
      & 6.94512792358934E-04_wp,-7.43041114473191E-04_wp, 1.90087784391754E-03_wp, &
      & 1.41762342890851E-03_wp,-8.72768819206869E-05_wp,-6.27592023471956E-05_wp, &
      &-9.94397095236798E-05_wp, 1.69999823758290E-05_wp, 1.02867249071332E-05_wp, &
      & 6.86263434466547E-04_wp, 1.12323534619302E-03_wp, 1.84122040819053E-03_wp, &
      &-1.05685564533451E-03_wp,-7.41528292047348E-05_wp,-5.64959811539522E-06_wp, &
      &-9.38849407025415E-05_wp,-6.44076461307482E-05_wp, 4.02013059080795E-05_wp, &
      & 7.61032562708509E-04_wp,-6.28894926244619E-04_wp, 2.82293389951878E-06_wp, &
      &-1.20209486448796E-04_wp, 1.34881283254919E-03_wp, 4.82420477979836E-04_wp, &
      &-1.52509696471363E-03_wp, 9.60493988868284E-06_wp,-1.30225149564277E-03_wp, &
      & 3.27030474840371E-04_wp,-8.31609991216725E-04_wp, 5.54949096390952E-04_wp, &
      & 8.00501008899362E-04_wp, 3.77777202539428E-04_wp, 1.19451943984326E-03_wp, &
      & 1.22592371170185E-03_wp,-2.06224126752325E-03_wp,-3.44948217512164E-04_wp, &
      & 1.57206108208992E-02_wp, 1.23464513582564E-02_wp, 3.36883951800672E-04_wp, &
      & 4.09324704116455E-04_wp, 1.22109199673396E-04_wp, 4.58251142523724E-04_wp, &
      & 2.38692760788248E-04_wp,-7.45009963602636E-04_wp,-8.81995629091509E-04_wp, &
      & 4.45232108062044E-04_wp,-3.31802835554287E-03_wp, 9.72796042636198E-05_wp, &
      & 4.59931850098990E-06_wp, 8.59895681224011E-06_wp, 8.83186095660427E-05_wp, &
      &-1.38666929038775E-04_wp, 9.28623086750347E-04_wp, 1.19596404451573E-04_wp, &
      & 5.87522199738515E-06_wp,-2.03335470252928E-03_wp,-4.17552021346422E-05_wp, &
      &-5.09405249892863E-06_wp,-6.69138499163273E-05_wp, 5.12364190514972E-05_wp, &
      & 1.03787374695023E-04_wp,-4.79144450852202E-04_wp,-7.62484937972927E-05_wp, &
      & 2.86985578991552E-04_wp, 7.01975281929019E-04_wp,-1.20147515573047E-03_wp, &
      &-5.32810176688890E-04_wp,-8.34604782021031E-05_wp,-1.49543254126587E-03_wp, &
      &-1.72564923333505E-03_wp, 1.58749721679993E-03_wp, 1.57654392909966E-03_wp, &
      & 3.31242618197280E-04_wp,-5.93648368390723E-04_wp, 1.95523826272871E-04_wp, &
      & 1.48935194145954E-03_wp,-1.03037430243343E-04_wp, 1.29590186401136E-03_wp, &
      &-2.38477695614008E-02_wp, 4.82043997276390E-05_wp, 1.66020669497395E-02_wp, &
      & 1.96443072180799E-05_wp,-3.38483160385626E-04_wp, 1.93563736986049E-04_wp, &
      & 2.38692760788248E-04_wp, 9.55996160124030E-04_wp,-2.25331678846465E-03_wp, &
      & 1.15282575115821E-03_wp, 1.46723806272367E-03_wp,-2.97238116407435E-03_wp, &
      & 1.74422156274703E-04_wp,-9.53909049759476E-05_wp, 4.56054718637186E-05_wp, &
      & 1.19023687658626E-04_wp,-4.22374644343754E-05_wp, 2.15549509066774E-03_wp, &
      & 5.25101958893504E-03_wp, 1.29756344050034E-05_wp,-2.17978866795582E-03_wp, &
      &-2.95874145259878E-04_wp,-1.25680696991265E-04_wp,-1.08710889283027E-04_wp, &
      & 5.08755648891538E-05_wp, 3.70646570214242E-05_wp,-1.30082860861911E-03_wp, &
      &-8.00507020681381E-03_wp,-3.80262208043235E-03_wp, 5.46887348152095E-03_wp, &
      &-2.39915293849376E-02_wp, 1.83968711176338E-02_wp,-9.70631584445533E-03_wp, &
      &-1.19162210963841E-02_wp,-8.62878718054140E-03_wp, 5.59939284194209E-04_wp, &
      &-6.32238839106823E-04_wp,-2.19959469457463E-04_wp,-8.89167773593586E-04_wp, &
      & 6.68905535755028E-04_wp,-9.31069490873247E-04_wp, 6.85864626230428E-04_wp, &
      &-1.73578590097926E-03_wp, 5.97412493903298E-03_wp, 1.29832452985145E-03_wp, &
      &-2.28312252496091E-02_wp,-9.42305899810999E-04_wp, 7.13341494118537E-04_wp, &
      & 6.94512792358934E-04_wp,-7.45009963602636E-04_wp,-2.25331678846465E-03_wp, &
      & 9.80974511856630E-01_wp, 6.77159511266487E-02_wp, 3.31547266836877E-02_wp, &
      &-4.62531679479594E-02_wp, 3.56741815719038E-03_wp,-2.66711332001602E-03_wp, &
      & 1.45224189985198E-03_wp, 1.78354903062118E-03_wp, 1.32108871221236E-03_wp, &
      &-3.02557120515213E-03_wp, 1.59240225965453E-02_wp, 1.87604065635365E-02_wp, &
      & 4.69837569301331E-03_wp,-1.69113411385159E-03_wp,-6.61974620886249E-04_wp, &
      &-1.25193134896850E-03_wp,-1.80059796256980E-03_wp,-3.51490801006137E-04_wp, &
      &-3.51696354973918E-02_wp,-1.18379012605945E-02_wp,-1.56986199671817E-02_wp, &
      & 2.07465098916424E-02_wp, 7.91814684776186E-02_wp,-6.46585082511956E-02_wp, &
      & 5.89106335179422E-04_wp, 7.77031024200781E-02_wp,-3.71725125852140E-02_wp, &
      &-4.33600687024431E-02_wp, 6.52028457208073E-02_wp,-3.87427809377845E-03_wp, &
      & 4.01388520677411E-02_wp,-5.70282419272425E-03_wp, 7.87169218444709E-03_wp, &
      &-3.93198699151506E-02_wp,-4.78755981358135E-03_wp, 3.58652467089654E-03_wp, &
      &-1.01189621165404E-02_wp, 1.33634035581096E-02_wp,-2.89537305206217E-04_wp, &
      &-2.49338759770228E-03_wp,-7.43041114473191E-04_wp,-8.81995629091509E-04_wp, &
      & 1.15282575115821E-03_wp, 6.77159511266487E-02_wp, 8.93041610503400E-01_wp, &
      &-4.10359825689074E-03_wp, 1.29601523129676E-02_wp, 1.67095113530985E-02_wp, &
      &-1.36647926268533E-02_wp, 1.27251154417880E-02_wp,-3.25733892434158E-05_wp, &
      & 2.26849602998059E-02_wp, 9.66633023541495E-03_wp,-5.19297868799894E-03_wp, &
      &-3.84633749324139E-02_wp,-1.24613040121662E-02_wp, 1.38676485413549E-03_wp, &
      &-6.88153886369830E-04_wp, 2.13804102298866E-03_wp, 3.54546470887354E-03_wp, &
      &-6.75694770549573E-06_wp,-1.39407590739072E-02_wp,-1.74855485990918E-02_wp, &
      & 9.95289435444545E-03_wp, 1.12870324819581E-02_wp, 6.88501216429309E-02_wp, &
      & 3.76360397992007E-02_wp, 8.54634077370166E-02_wp,-3.00536616447836E-02_wp, &
      & 3.53696185434452E-02_wp,-1.57001763492994E-02_wp, 1.09811119079166E-04_wp, &
      & 4.92359767360659E-02_wp, 6.11538319289382E-03_wp,-3.25394854692952E-02_wp, &
      & 3.89347577592347E-03_wp,-4.33959140311928E-02_wp,-4.80964039210276E-03_wp, &
      &-1.13563354350246E-02_wp, 1.64501682042580E-02_wp, 1.41764609986090E-03_wp, &
      &-1.28022878279255E-03_wp, 1.68735866975632E-04_wp, 1.90087784391754E-03_wp, &
      & 4.45232108062044E-04_wp, 1.46723806272367E-03_wp, 3.31547266836877E-02_wp, &
      &-4.10359825689074E-03_wp, 9.08117769603523E-01_wp, 9.12873283600010E-03_wp, &
      &-3.41667105883292E-04_wp,-2.27240427085510E-02_wp,-1.60656013699116E-02_wp, &
      & 1.70173883505167E-02_wp, 3.01970480128171E-04_wp, 1.87221098337107E-02_wp, &
      &-2.32165734653989E-02_wp,-2.18459102876561E-02_wp,-3.26536302873989E-02_wp, &
      & 1.96920444748044E-03_wp, 5.26310641005301E-04_wp,-1.60834820235921E-05_wp, &
      & 3.65950247565864E-03_wp, 4.03502323498709E-03_wp, 2.24022500427104E-02_wp, &
      & 2.22060923985063E-02_wp, 1.08109379474327E-02_wp, 2.62360032992966E-03_wp, &
      &-1.27801136787665E-02_wp, 6.76726646968229E-02_wp,-6.23156571113277E-03_wp, &
      &-2.00732274112082E-03_wp,-9.55348684789663E-02_wp,-3.92643295028350E-02_wp, &
      &-2.52064987261363E-02_wp, 3.76026077550443E-04_wp,-3.43827040288937E-02_wp, &
      &-4.86038808135031E-03_wp,-4.72243223724468E-02_wp, 3.84544626969289E-02_wp, &
      &-2.26550292609627E-02_wp,-9.59168996535191E-03_wp,-1.13646125614409E-02_wp, &
      &-5.45979916231744E-02_wp,-5.00721188868593E-03_wp,-7.41309576023266E-04_wp, &
      & 1.41762342890851E-03_wp,-3.31802835554287E-03_wp,-2.97238116407435E-03_wp, &
      &-4.62531679479594E-02_wp, 1.29601523129676E-02_wp, 9.12873283600010E-03_wp, &
      & 9.09645248969603E-01_wp,-2.27265262407807E-02_wp,-5.32299199303875E-04_wp, &
      &-1.00495621815055E-02_wp,-1.37190880213870E-02_wp, 1.76759686356936E-02_wp, &
      & 1.35404942972656E-02_wp, 8.14988231623426E-03_wp,-1.12416390562618E-02_wp, &
      &-3.84377541807347E-03_wp,-1.75012611402529E-03_wp,-9.81189201589943E-04_wp, &
      &-3.98266659773776E-04_wp, 2.71148530315031E-03_wp, 1.45446242844221E-03_wp, &
      &-1.30001097214204E-03_wp,-8.41884382085804E-04_wp,-6.58426823850392E-04_wp, &
      & 3.30170019630296E-04_wp, 1.67318715111030E-03_wp,-2.56333530435437E-03_wp, &
      & 2.39237960359928E-04_wp, 1.10329669950774E-03_wp, 1.04090454647008E-03_wp, &
      & 4.03441445098017E-04_wp, 1.91383514382521E-03_wp,-1.50246094084492E-04_wp, &
      & 1.81283369635852E-03_wp,-6.16005794354335E-05_wp, 1.35155308784679E-03_wp, &
      &-1.77253159605633E-03_wp, 1.03244385712790E-03_wp,-6.71551375241214E-04_wp, &
      &-6.58322439410420E-04_wp, 5.34400854366388E-03_wp, 1.98157365462255E-04_wp, &
      &-5.88299190152080E-05_wp,-8.72768819206869E-05_wp, 9.72796042636198E-05_wp, &
      & 1.74422156274703E-04_wp, 3.56741815719038E-03_wp, 1.67095113530985E-02_wp, &
      &-3.41667105883292E-04_wp,-2.27265262407807E-02_wp, 9.20468008197577E-04_wp, &
      &-2.49874089186843E-04_wp, 5.11060765003503E-04_wp, 3.55781691261206E-04_wp, &
      &-2.00686966101574E-05_wp,-1.11069365429552E-03_wp,-1.77595705940930E-03_wp, &
      &-9.69766591163740E-04_wp, 1.79487344769609E-03_wp, 1.18227059103663E-04_wp, &
      & 2.23628169815625E-05_wp, 9.59540695118637E-05_wp,-6.27125251357207E-06_wp, &
      &-7.33810769916223E-05_wp, 9.03676498870888E-04_wp, 6.54018323036504E-04_wp, &
      &-2.72580344706362E-05_wp,-6.85792113938563E-04_wp,-2.26525736225413E-03_wp, &
      & 3.57609160919946E-04_wp,-1.86157192226547E-03_wp,-1.68670894927570E-04_wp, &
      &-3.59221400041849E-04_wp, 1.24040844923198E-03_wp,-1.05856014409741E-03_wp, &
      &-1.27777504135482E-03_wp,-6.03244529432564E-04_wp, 1.02574808923692E-03_wp, &
      &-2.00017008301805E-04_wp, 1.76953567734680E-03_wp, 8.28645456652537E-04_wp, &
      & 2.30409327479312E-03_wp,-9.95021813203161E-05_wp,-4.96044569641485E-06_wp, &
      & 6.86048439806152E-05_wp, 5.59832212371451E-05_wp,-6.27592023471956E-05_wp, &
      & 4.59931850098990E-06_wp,-9.53909049759476E-05_wp,-2.66711332001602E-03_wp, &
      &-1.36647926268533E-02_wp,-2.27240427085510E-02_wp,-5.32299199303875E-04_wp, &
      &-2.49874089186843E-04_wp, 8.05224008513408E-04_wp, 2.12701081422194E-04_wp, &
      &-4.33732331900178E-04_wp,-3.76373327444631E-04_wp,-1.92622491079814E-03_wp, &
      & 3.14448950646059E-03_wp, 3.64600111382421E-03_wp, 2.20063557970657E-03_wp, &
      &-1.22466148179089E-04_wp,-1.14817581188101E-05_wp,-6.53524280632195E-05_wp, &
      &-2.32450086452929E-04_wp,-1.56880573246909E-04_wp,-5.83171540616180E-04_wp, &
      &-9.79735922629530E-05_wp,-5.84203711930900E-04_wp, 6.13175091245107E-05_wp, &
      & 1.03312589059841E-04_wp,-2.02312098196746E-03_wp,-1.22914115729581E-03_wp, &
      & 1.46237940335649E-03_wp,-3.57354903215663E-04_wp, 2.73735765525016E-04_wp, &
      & 1.27412333131700E-03_wp,-1.02227177755345E-03_wp, 9.60330421289763E-04_wp, &
      & 4.82756606241590E-04_wp, 5.23748289028198E-04_wp,-2.92690161508324E-04_wp, &
      & 8.36183698869665E-04_wp, 4.33495391693635E-04_wp,-1.69154717257233E-03_wp, &
      & 2.40201841492479E-03_wp, 1.12897189556701E-04_wp,-6.75477375787617E-05_wp, &
      &-9.94397095236798E-05_wp, 8.59895681224011E-06_wp, 4.56054718637186E-05_wp, &
      & 1.45224189985198E-03_wp, 1.27251154417880E-02_wp,-1.60656013699116E-02_wp, &
      &-1.00495621815055E-02_wp, 5.11060765003503E-04_wp, 2.12701081422194E-04_wp, &
      & 5.85429402823081E-04_wp,-1.44929109080006E-04_wp, 1.19686671170196E-04_wp, &
      &-1.26135270204623E-03_wp, 3.43298532897914E-04_wp, 3.83817970991369E-06_wp, &
      & 2.16308617823652E-03_wp, 1.39871480504802E-05_wp,-1.14181893240379E-05_wp, &
      & 6.32963494893167E-05_wp,-6.58156831461232E-05_wp,-1.36918692981879E-04_wp, &
      &-5.86376985234874E-04_wp,-7.38014920537504E-04_wp, 1.67615675300633E-05_wp, &
      & 2.02457181360983E-04_wp, 1.25511859423494E-03_wp,-3.47311396600001E-04_wp, &
      & 1.62380339070335E-03_wp,-6.46986175422616E-04_wp, 1.79429741173201E-03_wp, &
      & 2.82584056419365E-04_wp, 3.48350771197910E-04_wp, 9.55873747394603E-04_wp, &
      & 6.66961474578803E-04_wp,-5.53868614428592E-04_wp, 8.65609593211044E-04_wp, &
      &-1.42187034413839E-03_wp, 3.56016852022877E-04_wp,-1.03517926424556E-03_wp, &
      & 8.98861322776528E-04_wp, 2.48198884031917E-03_wp, 8.04397669199020E-05_wp, &
      & 1.76396292110500E-05_wp, 1.69999823758290E-05_wp, 8.83186095660427E-05_wp, &
      & 1.19023687658626E-04_wp, 1.78354903062118E-03_wp,-3.25733892434158E-05_wp, &
      & 1.70173883505167E-02_wp,-1.37190880213870E-02_wp, 3.55781691261206E-04_wp, &
      &-4.33732331900178E-04_wp,-1.44929109080006E-04_wp, 5.39623595018980E-04_wp, &
      &-2.59213800285308E-04_wp,-4.02492245771770E-05_wp,-2.31816357724684E-03_wp, &
      &-6.97563672918741E-04_wp,-1.10536488762744E-03_wp, 1.02416470270342E-04_wp, &
      & 4.30064827780300E-05_wp, 1.06480846519699E-05_wp, 5.28437067393580E-05_wp, &
      & 9.07005439053409E-05_wp,-5.20172683811180E-04_wp, 1.67408952638742E-04_wp, &
      &-1.95672888720847E-04_wp, 6.65644753860528E-04_wp, 1.23392237280277E-03_wp, &
      &-5.79280405069388E-04_wp,-2.57899109827770E-04_wp, 1.68257298559226E-03_wp, &
      &-2.19262341202439E-03_wp,-2.02260087963467E-03_wp, 1.31355340280529E-03_wp, &
      &-4.76745213909973E-05_wp, 1.38916670200151E-04_wp,-4.61704646589813E-04_wp, &
      &-8.12160921508423E-04_wp,-3.56619680821674E-04_wp,-2.26721397739007E-03_wp, &
      &-9.32392910000516E-04_wp,-2.17156562704799E-03_wp,-2.72121814524756E-03_wp, &
      &-1.60993997929043E-04_wp,-1.27813548615544E-04_wp, 1.02867249071332E-05_wp, &
      &-1.38666929038775E-04_wp,-4.22374644343754E-05_wp, 1.32108871221236E-03_wp, &
      & 2.26849602998059E-02_wp, 3.01970480128171E-04_wp, 1.76759686356936E-02_wp, &
      &-2.00686966101574E-05_wp,-3.76373327444631E-04_wp, 1.19686671170196E-04_wp, &
      &-2.59213800285308E-04_wp, 9.32898777817163E-04_wp, 1.33666846700357E-03_wp, &
      &-3.22194483116882E-04_wp,-4.38986218965597E-03_wp,-8.07406407670481E-04_wp, &
      & 5.55356868898791E-06_wp,-6.21201449864920E-05_wp, 9.70425331312843E-05_wp, &
      & 2.32661759561976E-04_wp, 3.40426800813873E-05_wp,-1.20825923162132E-03_wp, &
      & 2.55413860729610E-03_wp,-3.63428057490289E-03_wp,-9.08181098067992E-03_wp, &
      &-1.35620731019103E-02_wp,-5.21056540703202E-03_wp,-9.88014923090812E-03_wp, &
      & 2.13600011676250E-02_wp, 2.19832369216965E-02_wp,-8.83664108598077E-04_wp, &
      & 5.25854327632037E-05_wp, 5.49263377585196E-04_wp,-8.40001714014691E-04_wp, &
      &-4.42109441645743E-04_wp, 1.03952145445938E-03_wp, 2.40006331976448E-04_wp, &
      &-4.67936049592108E-04_wp,-1.75636665246333E-02_wp, 1.74725131414713E-03_wp, &
      & 1.34763066806291E-02_wp, 9.24348836836794E-04_wp,-3.12515170166504E-04_wp, &
      & 6.86263434466547E-04_wp, 9.28623086750347E-04_wp, 2.15549509066774E-03_wp, &
      &-3.02557120515213E-03_wp, 9.66633023541495E-03_wp, 1.87221098337107E-02_wp, &
      & 1.35404942972656E-02_wp,-1.11069365429552E-03_wp,-1.92622491079814E-03_wp, &
      &-1.26135270204623E-03_wp,-4.02492245771770E-05_wp, 1.33666846700357E-03_wp, &
      & 9.80364656316394E-01_wp,-2.18477518470214E-02_wp, 3.31124366005084E-02_wp, &
      & 7.98537612496862E-02_wp, 1.97816816411409E-03_wp, 7.95287884239187E-04_wp, &
      & 1.47614735664940E-03_wp,-3.09112108744084E-03_wp,-3.26315136987645E-03_wp, &
      & 1.01661826271289E-02_wp, 1.44077180987701E-02_wp, 5.13544987508727E-03_wp, &
      & 1.33513181820185E-02_wp, 7.56618991833905E-02_wp, 3.66511396790250E-02_wp, &
      &-6.25144814937098E-03_wp, 3.25281161052200E-02_wp, 7.20780368313875E-02_wp, &
      & 3.13114544459522E-02_wp, 2.45048495855122E-02_wp,-6.37450174955805E-03_wp, &
      &-1.91461291146198E-02_wp, 1.09535528037961E-03_wp, 3.42404939671757E-02_wp, &
      & 4.74589429392366E-02_wp,-2.15316812457990E-02_wp,-4.80394095409376E-02_wp, &
      &-1.23060391666262E-02_wp, 1.29062067504193E-02_wp,-1.64629006596400E-03_wp, &
      &-3.56717280611389E-03_wp, 1.12323534619302E-03_wp, 1.19596404451573E-04_wp, &
      & 5.25101958893504E-03_wp, 1.59240225965453E-02_wp,-5.19297868799894E-03_wp, &
      &-2.32165734653989E-02_wp, 8.14988231623426E-03_wp,-1.77595705940930E-03_wp, &
      & 3.14448950646059E-03_wp, 3.43298532897914E-04_wp,-2.31816357724684E-03_wp, &
      &-3.22194483116882E-04_wp,-2.18477518470214E-02_wp, 9.14810407256861E-01_wp, &
      & 7.48033986208798E-03_wp, 5.09953559406378E-03_wp,-2.71827543681093E-02_wp, &
      &-1.35396965598036E-02_wp,-5.28749471868861E-03_wp,-4.55353025074698E-04_wp, &
      &-8.97158346345382E-03_wp,-1.38165545028930E-02_wp, 4.81975145073823E-03_wp, &
      & 1.02100793239385E-02_wp,-1.99268786228131E-02_wp, 3.15084854540213E-02_wp, &
      &-1.61959501492248E-02_wp, 8.60612099656798E-02_wp, 4.67381503378004E-02_wp, &
      &-7.10320337351810E-02_wp, 3.67817089890598E-02_wp,-2.56655009577307E-03_wp, &
      &-1.46067658017428E-02_wp, 5.62002262438386E-03_wp, 5.78076943435148E-02_wp, &
      &-2.44882129383753E-03_wp,-2.86025601155137E-02_wp,-5.04588128389885E-03_wp, &
      &-4.04081526462394E-03_wp, 1.53478173145738E-02_wp,-9.80190748893888E-03_wp, &
      &-1.88908493067740E-03_wp, 3.00321001380816E-04_wp, 1.84122040819053E-03_wp, &
      & 5.87522199738515E-06_wp, 1.29756344050034E-05_wp, 1.87604065635365E-02_wp, &
      &-3.84633749324139E-02_wp,-2.18459102876561E-02_wp,-1.12416390562618E-02_wp, &
      &-9.69766591163740E-04_wp, 3.64600111382421E-03_wp, 3.83817970991369E-06_wp, &
      &-6.97563672918741E-04_wp,-4.38986218965597E-03_wp, 3.31124366005084E-02_wp, &
      & 7.48033986208798E-03_wp, 9.07045088096441E-01_wp,-6.65538482288344E-03_wp, &
      &-4.69539950187602E-04_wp, 8.56203889791274E-03_wp,-1.56700021141420E-02_wp, &
      &-2.67566253419244E-02_wp,-5.14238323186605E-05_wp,-4.08742499801316E-02_wp, &
      & 1.15402494798051E-02_wp,-1.81555166175954E-02_wp,-2.34217202320817E-02_wp, &
      & 6.45053437127255E-02_wp, 4.22610487435240E-02_wp, 3.14099372304354E-03_wp, &
      &-1.03858523700384E-01_wp,-4.72446672147925E-02_wp, 5.10642880072331E-02_wp, &
      & 4.06754362288397E-02_wp,-4.95983018633171E-03_wp, 4.94693462764905E-02_wp, &
      &-1.68594875769831E-03_wp,-6.15913706401558E-02_wp,-2.74869037414883E-02_wp, &
      & 2.25098483451193E-03_wp, 2.91067062790983E-02_wp,-5.75999068005323E-03_wp, &
      &-5.04907733624908E-03_wp,-5.38795740743565E-04_wp,-6.77069550327497E-04_wp, &
      &-1.05685564533451E-03_wp,-2.03335470252928E-03_wp,-2.17978866795582E-03_wp, &
      & 4.69837569301331E-03_wp,-1.24613040121662E-02_wp,-3.26536302873989E-02_wp, &
      &-3.84377541807347E-03_wp, 1.79487344769609E-03_wp, 2.20063557970657E-03_wp, &
      & 2.16308617823652E-03_wp,-1.10536488762744E-03_wp,-8.07406407670481E-04_wp, &
      & 7.98537612496862E-02_wp, 5.09953559406378E-03_wp,-6.65538482288344E-03_wp, &
      & 8.85720560609372E-01_wp, 8.47405040928117E-03_wp, 4.32609508448040E-05_wp, &
      & 1.51332117219224E-02_wp,-1.32002681370741E-02_wp,-2.63946101633402E-02_wp, &
      &-7.18067691148890E-04_wp,-3.91472860045293E-04_wp,-3.91264029795281E-04_wp, &
      &-6.69813653944916E-04_wp,-9.38359904068336E-04_wp,-7.23150166635905E-04_wp, &
      & 3.78532516914578E-04_wp,-1.50740974288056E-03_wp,-2.04798561893996E-03_wp, &
      &-6.93371177426215E-04_wp,-4.24473085272663E-04_wp, 1.86494551430705E-04_wp, &
      & 1.30712345342734E-03_wp,-1.50326172865105E-04_wp,-1.62689399094736E-03_wp, &
      &-1.81297319011831E-03_wp, 2.05812072051282E-03_wp, 5.41082241160368E-03_wp, &
      & 6.00827861785183E-04_wp,-2.40085902556881E-03_wp, 5.99810932075002E-05_wp, &
      & 1.50498405466455E-04_wp,-7.41528292047348E-05_wp,-4.17552021346422E-05_wp, &
      &-2.95874145259878E-04_wp,-1.69113411385159E-03_wp, 1.38676485413549E-03_wp, &
      & 1.96920444748044E-03_wp,-1.75012611402529E-03_wp, 1.18227059103663E-04_wp, &
      &-1.22466148179089E-04_wp, 1.39871480504802E-05_wp, 1.02416470270342E-04_wp, &
      & 5.55356868898791E-06_wp, 1.97816816411409E-03_wp,-2.71827543681093E-02_wp, &
      &-4.69539950187602E-04_wp, 8.47405040928117E-03_wp, 9.24530237513042E-04_wp, &
      & 4.12508019867131E-04_wp, 3.21310649550534E-04_wp,-1.04518495185366E-04_wp, &
      & 1.07469239917818E-05_wp,-2.63641797653816E-04_wp,-1.85447400076972E-04_wp, &
      & 9.88365504774939E-06_wp,-4.38712775661301E-04_wp,-6.34043613156057E-04_wp, &
      &-6.96300015777383E-04_wp, 9.38872653962262E-04_wp, 3.82651994861127E-05_wp, &
      &-1.49699282601155E-03_wp,-1.52128659470300E-04_wp,-4.73872134127600E-04_wp, &
      &-2.30162533491128E-05_wp, 4.10116747608054E-04_wp, 5.22664015461809E-04_wp, &
      &-5.53268037298817E-04_wp,-9.90759348859535E-04_wp, 5.05766248283295E-04_wp, &
      & 2.05070858833928E-03_wp, 6.76325363903052E-04_wp,-1.16119924081096E-03_wp, &
      & 1.05424830110665E-05_wp, 8.35670166748990E-05_wp,-5.64959811539522E-06_wp, &
      &-5.09405249892863E-06_wp,-1.25680696991265E-04_wp,-6.61974620886249E-04_wp, &
      &-6.88153886369830E-04_wp, 5.26310641005301E-04_wp,-9.81189201589943E-04_wp, &
      & 2.23628169815625E-05_wp,-1.14817581188101E-05_wp,-1.14181893240379E-05_wp, &
      & 4.30064827780300E-05_wp,-6.21201449864920E-05_wp, 7.95287884239187E-04_wp, &
      &-1.35396965598036E-02_wp, 8.56203889791274E-03_wp, 4.32609508448040E-05_wp, &
      & 4.12508019867131E-04_wp, 2.88003393911521E-04_wp,-6.50182407583703E-05_wp, &
      &-2.53634831707529E-04_wp, 1.27233371965883E-04_wp,-5.93966433083350E-04_wp, &
      & 2.69372961745411E-05_wp,-5.75792938733305E-04_wp,-1.20006019567996E-04_wp, &
      & 3.14862311605364E-04_wp, 6.92000783727210E-04_wp,-1.19990637124721E-03_wp, &
      &-2.40313487407472E-03_wp, 1.86593826628628E-04_wp,-6.60971351199133E-05_wp, &
      & 6.75497357337969E-04_wp, 1.12528083885423E-04_wp, 9.58776573157139E-04_wp, &
      &-1.09131629766252E-03_wp,-1.20341681703315E-03_wp,-3.97806973349569E-04_wp, &
      & 7.91454294590884E-04_wp, 2.29339909979128E-03_wp,-1.66014282991117E-03_wp, &
      &-5.81093006974993E-04_wp, 4.34016120608695E-05_wp,-1.78507090937967E-05_wp, &
      &-9.38849407025415E-05_wp,-6.69138499163273E-05_wp,-1.08710889283027E-04_wp, &
      &-1.25193134896850E-03_wp, 2.13804102298866E-03_wp,-1.60834820235921E-05_wp, &
      &-3.98266659773776E-04_wp, 9.59540695118637E-05_wp,-6.53524280632195E-05_wp, &
      & 6.32963494893167E-05_wp, 1.06480846519699E-05_wp, 9.70425331312843E-05_wp, &
      & 1.47614735664940E-03_wp,-5.28749471868861E-03_wp,-1.56700021141420E-02_wp, &
      & 1.51332117219224E-02_wp, 3.21310649550534E-04_wp,-6.50182407583703E-05_wp, &
      & 5.65537994689567E-04_wp, 2.39816673197796E-04_wp,-4.01145908943532E-04_wp, &
      & 1.02979125404260E-03_wp,-3.83960703078049E-04_wp,-4.71112289724211E-05_wp, &
      & 1.02053886193382E-03_wp,-1.37787803489053E-03_wp, 2.14568368478059E-04_wp, &
      &-2.29785457537120E-03_wp, 3.71631530998101E-04_wp, 2.32423304833064E-03_wp, &
      &-2.04580162271769E-03_wp,-5.61485634191755E-04_wp, 5.83275218807406E-04_wp, &
      &-7.52025031398366E-04_wp,-1.82188558275370E-03_wp, 1.03007344334789E-03_wp, &
      & 1.22820555419500E-03_wp, 7.16816433450587E-04_wp, 4.07757821937252E-04_wp, &
      &-2.79625792410964E-04_wp, 2.45027490200436E-03_wp, 1.16049685662807E-04_wp, &
      & 5.03854201102447E-06_wp,-6.44076461307482E-05_wp, 5.12364190514972E-05_wp, &
      & 5.08755648891538E-05_wp,-1.80059796256980E-03_wp, 3.54546470887354E-03_wp, &
      & 3.65950247565864E-03_wp, 2.71148530315031E-03_wp,-6.27125251357207E-06_wp, &
      &-2.32450086452929E-04_wp,-6.58156831461232E-05_wp, 5.28437067393580E-05_wp, &
      & 2.32661759561976E-04_wp,-3.09112108744084E-03_wp,-4.55353025074698E-04_wp, &
      &-2.67566253419244E-02_wp,-1.32002681370741E-02_wp,-1.04518495185366E-04_wp, &
      &-2.53634831707529E-04_wp, 2.39816673197796E-04_wp, 1.01972520181169E-03_wp, &
      & 4.18800952996832E-04_wp, 1.20885516832022E-03_wp,-5.52056921452457E-04_wp, &
      & 5.46257757174961E-04_wp, 5.94007597988431E-04_wp,-2.24225399100506E-03_wp, &
      &-1.19565402301583E-03_wp, 4.74724238068227E-05_wp, 2.60157115829106E-03_wp, &
      & 4.22160666294196E-04_wp,-1.86857557022708E-03_wp,-1.60968883271744E-03_wp, &
      & 4.19280315359715E-04_wp,-1.22301380113740E-03_wp,-3.25003079162797E-05_wp, &
      & 1.53200056222922E-03_wp, 4.56941288207854E-04_wp, 1.18902455312498E-03_wp, &
      &-4.27488088447255E-04_wp, 2.06014446932532E-03_wp, 1.39966952488807E-03_wp, &
      & 5.66604383383426E-05_wp, 9.88365521261211E-05_wp, 4.02013059080795E-05_wp, &
      & 1.03787374695023E-04_wp, 3.70646570214242E-05_wp,-3.51490801006137E-04_wp, &
      &-6.75694770549573E-06_wp, 4.03502323498709E-03_wp, 1.45446242844221E-03_wp, &
      &-7.33810769916223E-05_wp,-1.56880573246909E-04_wp,-1.36918692981879E-04_wp, &
      & 9.07005439053409E-05_wp, 3.40426800813873E-05_wp,-3.26315136987645E-03_wp, &
      &-8.97158346345382E-03_wp,-5.14238323186605E-05_wp,-2.63946101633402E-02_wp, &
      & 1.07469239917818E-05_wp, 1.27233371965883E-04_wp,-4.01145908943532E-04_wp, &
      & 4.18800952996832E-04_wp, 8.91699965008798E-04_wp], shape(density))

   call get_structure(mol, "f-block", "CeCl3")
   call test_numpot(error, mol, density, qsh, make_exchange_gxtb, thr_in=thr1*10)

end subroutine test_p_fock_cecl3


subroutine test_p_fock_ce2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qsh(8) = [&
      & -5.81684067068857E-1_wp, -3.32402447916258E-1_wp,  1.04173946069987E+0_wp, &
      & -1.27652922694715E-1_wp, -5.81684070494655E-1_wp, -3.32402447074014E-1_wp, &
      &  1.04173945811682E+0_wp, -1.27652963568184E-1_wp]
   
   real(wp), parameter :: density(32, 32, 1) = reshape([&
      & 1.34115685561296E+00_wp, 2.95040209552549E-17_wp,-3.27055576896394E-01_wp, &
      & 1.43604125002084E-17_wp,-3.41890335403403E-17_wp, 7.27095612532062E-19_wp, &
      &-1.06435520703580E-01_wp, 2.73735199807599E-18_wp, 4.96021353460711E-17_wp, &
      & 1.72433300486755E-17_wp, 3.10850086749740E-17_wp, 4.13920779621525E-17_wp, &
      & 5.98906667478317E-02_wp, 1.29262174241255E-17_wp, 4.44083909509455E-17_wp, &
      &-5.01528789376767E-27_wp,-1.33570401925884E-01_wp,-3.38149550350916E-16_wp, &
      &-4.25155301003653E-01_wp,-1.28989578591008E-17_wp,-5.17357619198806E-17_wp, &
      &-8.39801861579282E-17_wp, 1.81108201463904E-01_wp,-9.72496982710943E-18_wp, &
      & 7.96062373381183E-17_wp,-3.18260052314977E-17_wp, 1.39123540364297E-16_wp, &
      & 1.33855770960173E-17_wp, 8.92061338128334E-02_wp, 1.44260904550692E-17_wp, &
      &-2.01125186749760E-16_wp, 9.17410126990192E-28_wp, 2.95040209552549E-17_wp, &
      & 2.77140754922035E-03_wp,-1.14156519725075E-17_wp,-3.85456657566366E-18_wp, &
      &-1.56769994717750E-17_wp, 2.46293500505109E-02_wp,-1.12490353247610E-17_wp, &
      &-1.30104260698261E-17_wp, 6.43443161089344E-19_wp,-1.87671993779360E-16_wp, &
      &-2.44141767132226E-17_wp, 4.19084840085586E-02_wp,-2.52382312583193E-17_wp, &
      &-7.21140625583997E-18_wp,-1.07163774483076E-17_wp, 1.97655006366489E-21_wp, &
      &-1.82787972147955E-17_wp, 2.77140751053707E-03_wp,-1.36420649358352E-17_wp, &
      &-4.04163230174865E-18_wp,-1.33561836867554E-17_wp,-2.46293501150133E-02_wp, &
      & 1.43755402599569E-17_wp, 6.93889390390723E-18_wp,-2.10392103991529E-18_wp, &
      &-1.77025310061228E-16_wp, 2.89803562889256E-17_wp, 4.19084846595035E-02_wp, &
      & 3.12830876064525E-17_wp,-3.17171050218737E-17_wp, 2.06275819441853E-17_wp, &
      &-2.50285987685287E-21_wp,-3.27055576896394E-01_wp,-1.14156519725075E-17_wp, &
      & 1.95823894853419E-01_wp,-1.15221139896263E-18_wp, 1.89167741369865E-17_wp, &
      &-3.77018592081629E-17_wp, 7.63664572202967E-02_wp, 4.64417497426684E-18_wp, &
      &-3.42277748132774E-17_wp,-1.13721493819988E-17_wp, 3.84564801954096E-18_wp, &
      & 5.40539925312216E-18_wp,-3.92154752318439E-02_wp, 1.05106929051690E-18_wp, &
      &-3.32753192075161E-17_wp, 3.43293026134509E-27_wp, 4.25155300877227E-01_wp, &
      &-2.68102246439161E-16_wp, 1.87854646548846E-01_wp, 1.45063455819008E-17_wp, &
      & 1.54136010738722E-17_wp,-7.89082478289947E-17_wp,-7.03003339238226E-02_wp, &
      & 1.97389772270935E-18_wp,-2.93432673339620E-17_wp, 6.33107175306980E-18_wp, &
      &-6.45162606238304E-17_wp, 1.51647595042051E-17_wp,-3.68340001692192E-02_wp, &
      & 2.54424582926323E-18_wp, 9.41050782841875E-17_wp,-3.31686667557185E-28_wp, &
      & 1.43604125002084E-17_wp,-3.85456657566366E-18_wp,-1.15221139896263E-18_wp, &
      & 2.77140754922034E-03_wp, 1.16292082144050E-18_wp,-2.16635673099271E-17_wp, &
      & 6.95835778407141E-19_wp, 2.46293500505109E-02_wp,-1.84280626035159E-20_wp, &
      &-4.25800029707754E-20_wp, 4.27771513008680E-18_wp,-3.44426121919881E-17_wp, &
      & 2.78304841488652E-18_wp, 4.19084840085585E-02_wp, 1.60468047089764E-18_wp, &
      & 1.73226265634137E-16_wp, 4.54019653199472E-18_wp,-1.27851740368253E-18_wp, &
      &-3.05601911633501E-18_wp, 2.77140751053706E-03_wp, 3.10068186843677E-18_wp, &
      & 2.49944922548076E-17_wp, 2.71256570939184E-18_wp,-2.46293501150133E-02_wp, &
      & 1.22079844037651E-18_wp,-9.93988121204300E-19_wp,-3.60039053674948E-18_wp, &
      &-5.29495672833062E-17_wp,-3.58333013113701E-18_wp, 4.19084846595034E-02_wp, &
      & 8.79472678323498E-19_wp, 1.64987512999437E-16_wp,-3.41890335403403E-17_wp, &
      &-1.56769994717750E-17_wp, 1.89167741369865E-17_wp, 1.16292082144050E-18_wp, &
      & 9.10037888950928E-32_wp,-1.39320652366933E-16_wp, 7.30816180766086E-18_wp, &
      & 1.03348148850728E-17_wp,-6.92923611264855E-33_wp, 1.06049198889524E-30_wp, &
      & 1.40147514418647E-31_wp,-2.37063394681843E-16_wp,-3.76995335331624E-18_wp, &
      & 1.75853777486710E-17_wp, 5.81146219295055E-32_wp, 7.26769488701110E-32_wp, &
      & 3.91883209360840E-17_wp,-1.56769992529556E-17_wp, 1.85106510366749E-17_wp, &
      & 1.16292080520847E-18_wp, 7.84267347127136E-32_wp, 1.39320652731803E-16_wp, &
      &-6.99902444140602E-18_wp,-1.03348149121388E-17_wp, 9.47901645394113E-33_wp, &
      & 1.00164156833194E-30_wp,-1.71778478447135E-31_wp,-2.37063398364039E-16_wp, &
      &-3.64859041169069E-18_wp, 1.75853780218166E-17_wp,-1.07080703791633E-31_wp, &
      & 6.92451933067693E-32_wp, 7.27095612532062E-19_wp, 2.46293500505109E-02_wp, &
      &-3.77018592081629E-17_wp,-2.16635673099271E-17_wp,-1.39320652366933E-16_wp, &
      & 2.18879711171045E-01_wp,-7.92251697232656E-17_wp,-2.41100800849787E-18_wp, &
      & 5.71824481626132E-18_wp,-1.66783092972721E-15_wp,-2.16967477278523E-16_wp, &
      & 3.72438446674291E-01_wp,-2.35964040632013E-16_wp, 1.28397204413363E-16_wp, &
      &-9.52358708563163E-17_wp, 1.75654942652798E-20_wp,-1.36451437807371E-16_wp, &
      & 2.46293497067347E-02_wp,-3.83581862621285E-17_wp,-2.27945384644813E-17_wp, &
      &-1.18695687125686E-16_wp,-2.18879711744274E-01_wp, 9.24488963624558E-17_wp, &
      &-5.55111512312578E-17_wp,-1.86974332899130E-17_wp,-1.57321442330803E-15_wp, &
      & 2.57546869939526E-16_wp, 3.72438452459204E-01_wp, 2.60621297069682E-16_wp, &
      &-9.34514695314122E-17_wp, 1.83316212926483E-16_wp,-2.22427813078371E-20_wp, &
      &-1.06435520703580E-01_wp,-1.12490353247610E-17_wp, 7.63664572202967E-02_wp, &
      & 6.95835778407141E-19_wp, 7.30816180766086E-18_wp,-7.92251697232656E-17_wp, &
      & 3.03415147246508E-02_wp, 9.33216879207763E-18_wp,-1.35488181338000E-17_wp, &
      &-4.48132144278951E-18_wp, 2.49567157860423E-18_wp,-1.03530147053401E-16_wp, &
      &-1.54418973064381E-02_wp, 1.31231813173097E-17_wp,-1.32730563661014E-17_wp, &
      & 1.35787992059083E-27_wp, 1.81108199345750E-01_wp,-1.32496562140540E-16_wp, &
      & 7.03003330580197E-02_wp, 6.77288398570413E-18_wp, 5.32072401654417E-18_wp, &
      & 2.63679918315433E-17_wp,-2.57240327880902E-02_wp,-6.64332971457817E-18_wp, &
      &-1.06306472855005E-17_wp, 1.90464148457209E-18_wp,-2.43267135760346E-17_wp, &
      &-1.00035099535184E-16_wp,-1.36291380722532E-02_wp, 1.38115301449283E-17_wp, &
      & 3.55314728153116E-17_wp,-1.19650353049022E-28_wp, 2.73735199807599E-18_wp, &
      &-1.30104260698261E-17_wp, 4.64417497426684E-18_wp, 2.46293500505109E-02_wp, &
      & 1.03348148850728E-17_wp,-2.41100800849787E-18_wp, 9.33216879207763E-18_wp, &
      & 2.18879711171045E-01_wp,-1.63769202671905E-19_wp,-3.78406199628771E-19_wp, &
      & 3.80158246248963E-17_wp, 1.53505234060622E-17_wp, 2.24574725183172E-17_wp, &
      & 3.72438446674291E-01_wp, 1.42607091649422E-17_wp, 1.53945251951349E-15_wp, &
      & 1.21988114168630E-19_wp, 1.00254756044942E-17_wp, 1.13799341050549E-18_wp, &
      & 2.46293497067347E-02_wp, 2.75555932415949E-17_wp, 2.85084061341327E-17_wp, &
      & 1.07483602006170E-17_wp,-2.18879711744274E-01_wp, 1.08491701762195E-17_wp, &
      &-8.83351905066482E-18_wp,-3.19964773398613E-17_wp,-1.36262959273222E-16_wp, &
      &-3.81283996669708E-17_wp, 3.72438452459204E-01_wp, 7.81582645987379E-18_wp, &
      & 1.46623516731399E-15_wp, 4.96021353460711E-17_wp, 6.43443161089344E-19_wp, &
      &-3.42277748132774E-17_wp,-1.84280626035159E-20_wp,-6.92923611264855E-33_wp, &
      & 5.71824481626132E-18_wp,-1.35488181338000E-17_wp,-1.63769202671905E-19_wp, &
      & 6.20409866477369E-33_wp,-4.15675230592345E-32_wp,-6.72578100587264E-33_wp, &
      & 9.72997545399451E-18_wp, 6.90774688247140E-18_wp,-2.78664235848457E-19_wp, &
      & 3.42367704410442E-33_wp,-1.15138340899647E-33_wp,-7.97975116639943E-17_wp, &
      & 6.43443152108223E-19_wp,-3.17748173780746E-17_wp,-1.84280623462999E-20_wp, &
      &-5.56835826553294E-33_wp,-5.71824483123694E-18_wp, 1.16816478879813E-17_wp, &
      & 1.63769203100801E-19_wp, 4.34114355013797E-33_wp,-4.19980981347647E-32_wp, &
      & 1.77305803405670E-32_wp, 9.72997560512568E-18_wp, 6.17472181633275E-18_wp, &
      &-2.78664240176825E-19_wp,-1.12469707088764E-32_wp,-1.09764104060678E-33_wp, &
      & 1.72433300486755E-17_wp,-1.87671993779360E-16_wp,-1.13721493819988E-17_wp, &
      &-4.25800029707754E-20_wp, 1.06049198889524E-30_wp,-1.66783092972721E-15_wp, &
      &-4.48132144278951E-18_wp,-3.78406199628771E-19_wp,-4.15675230592345E-32_wp, &
      & 1.27092891975898E-29_wp, 1.65288837181326E-30_wp,-2.83792571481203E-15_wp, &
      & 2.28971263478392E-18_wp,-6.43883420934220E-19_wp, 7.27615544040963E-31_wp, &
      &-2.79530089298938E-33_wp,-2.59592601471085E-17_wp,-1.87671991159836E-16_wp, &
      &-1.06641003950090E-17_wp,-4.25800023764352E-20_wp, 9.03557987996926E-31_wp, &
      & 1.66783093409513E-15_wp, 3.94236060996094E-18_wp, 3.78406200620225E-19_wp, &
      & 1.44089616294457E-31_wp, 1.19873549172850E-29_wp,-1.95873601060336E-30_wp, &
      &-2.83792575889221E-15_wp, 2.07812411870623E-18_wp,-6.43883430933571E-19_wp, &
      &-1.40222409804720E-30_wp,-2.36538725322292E-33_wp, 3.10850086749740E-17_wp, &
      &-2.44141767132226E-17_wp, 3.84564801954096E-18_wp, 4.27771513008680E-18_wp, &
      & 1.40147514418647E-31_wp,-2.16967477278523E-16_wp, 2.49567157860423E-18_wp, &
      & 3.80158246248963E-17_wp,-6.72578100587264E-33_wp, 1.65288837181326E-30_wp, &
      & 2.23519972446692E-31_wp,-3.69184652995568E-16_wp,-1.03459788273466E-18_wp, &
      & 6.46864645269670E-17_wp, 9.57003491444447E-32_wp, 2.67360302817163E-31_wp, &
      & 3.55511390567388E-17_wp,-2.44141763724498E-17_wp,-1.56761885797389E-18_wp, &
      & 4.27771507037849E-18_wp, 1.21520915923955E-31_wp, 2.16967477846744E-16_wp, &
      & 1.62486129743607E-18_wp,-3.80158247244567E-17_wp, 2.12859207917019E-32_wp, &
      & 1.55705742722027E-30_wp,-2.60640833892739E-31_wp,-3.69184658729941E-16_wp, &
      & 5.83066485252031E-19_wp, 6.46864655317119E-17_wp,-1.80583133059065E-31_wp, &
      & 2.54683107179597E-31_wp, 4.13920779621525E-17_wp, 4.19084840085586E-02_wp, &
      & 5.40539925312216E-18_wp,-3.44426121919881E-17_wp,-2.37063394681843E-16_wp, &
      & 3.72438446674291E-01_wp,-1.03530147053401E-16_wp, 1.53505234060622E-17_wp, &
      & 9.72997545399451E-18_wp,-2.83792571481203E-15_wp,-3.69184652995568E-16_wp, &
      & 6.33728890718255E-01_wp,-4.16540388168571E-16_wp, 2.57766699401504E-16_wp, &
      &-1.62050194692015E-16_wp, 2.98888616229728E-20_wp, 3.22089028957633E-17_wp, &
      & 4.19084834236004E-02_wp,-2.04514622825458E-17_wp,-3.67026106494832E-17_wp, &
      &-2.01968638863391E-16_wp,-3.72438447649679E-01_wp, 1.44863137958540E-16_wp, &
      &-1.11022302462516E-16_wp,-3.18149314709649E-17_wp,-2.67692941007473E-15_wp, &
      & 4.38233199746618E-16_wp, 6.33728900561671E-01_wp, 4.35825887214287E-16_wp, &
      &-1.36729465914179E-16_wp, 3.11924779264715E-16_wp,-3.78475779031596E-20_wp, &
      & 5.98906667478317E-02_wp,-2.52382312583193E-17_wp,-3.92154752318439E-02_wp, &
      & 2.78304841488652E-18_wp,-3.76995335331624E-18_wp,-2.35964040632013E-16_wp, &
      &-1.54418973064381E-02_wp, 2.24574725183172E-17_wp, 6.90774688247140E-18_wp, &
      & 2.28971263478392E-18_wp,-1.03459788273466E-18_wp,-4.16540388168571E-16_wp, &
      & 7.89276986890700E-03_wp, 3.96073062894617E-17_wp, 6.74241441322119E-18_wp, &
      &-6.92374106977196E-28_wp,-8.92061298959614E-02_wp, 3.17807754267225E-17_wp, &
      &-3.68339985833876E-02_wp,-3.44909746440281E-19_wp,-2.90343249620542E-18_wp, &
      & 2.61276534738854E-16_wp, 1.36291376558167E-02_wp,-2.38114673055874E-17_wp, &
      & 5.66049568925336E-18_wp,-1.11800158227415E-18_wp, 1.26987491254611E-17_wp, &
      &-4.18412246384521E-16_wp, 7.18110568545843E-03_wp, 3.92801248949416E-17_wp, &
      &-1.85354424939333E-17_wp, 6.40177057547808E-29_wp, 1.29262174241255E-17_wp, &
      &-7.21140625583997E-18_wp, 1.05106929051690E-18_wp, 4.19084840085585E-02_wp, &
      & 1.75853777486710E-17_wp, 1.28397204413363E-16_wp, 1.31231813173097E-17_wp, &
      & 3.72438446674291E-01_wp,-2.78664235848457E-19_wp,-6.43883420934220E-19_wp, &
      & 6.46864645269670E-17_wp, 2.57766699401504E-16_wp, 3.96073062894617E-17_wp, &
      & 6.33728890718256E-01_wp, 2.42655490609378E-17_wp, 2.61948127594329E-15_wp, &
      &-1.69694570077344E-17_wp, 3.23861724163708E-17_wp,-4.19123223859288E-18_wp, &
      & 4.19084834236003E-02_wp, 4.68876822304845E-17_wp,-8.20267603302476E-17_wp, &
      & 2.04943051581071E-17_wp,-3.72438447649679E-01_wp, 1.84605876283281E-17_wp, &
      &-1.50308226207710E-17_wp,-5.44441431129022E-17_wp,-7.26596145091637E-18_wp, &
      &-6.36998672946740E-17_wp, 6.33728900561672E-01_wp, 1.32991507098458E-17_wp, &
      & 2.49489706127629E-15_wp, 4.44083909509455E-17_wp,-1.07163774483076E-17_wp, &
      &-3.32753192075161E-17_wp, 1.60468047089764E-18_wp, 5.81146219295055E-32_wp, &
      &-9.52358708563163E-17_wp,-1.32730563661014E-17_wp, 1.42607091649422E-17_wp, &
      & 3.42367704410442E-33_wp, 7.27615544040963E-31_wp, 9.57003491444447E-32_wp, &
      &-1.62050194692015E-16_wp, 6.74241441322119E-18_wp, 2.42655490609378E-17_wp, &
      & 4.81779831110240E-32_wp, 1.00292583863982E-31_wp,-8.03426297896837E-17_wp, &
      &-1.07163772987285E-17_wp,-3.03561584380044E-17_wp, 1.60468044849953E-18_wp, &
      & 5.11865097234841E-32_wp, 9.52358711057317E-17_wp, 1.10510160861017E-17_wp, &
      &-1.42607092022899E-17_wp, 1.33985466103712E-32_wp, 6.83162195812944E-31_wp, &
      &-1.03622630068020E-31_wp,-1.62050197209065E-16_wp, 5.87007225823383E-18_wp, &
      & 2.42655494378433E-17_wp,-9.46260833251232E-32_wp, 9.55395613748457E-32_wp, &
      &-5.01528789376767E-27_wp, 1.97655006366489E-21_wp, 3.43293026134509E-27_wp, &
      & 1.73226265634137E-16_wp, 7.26769488701110E-32_wp, 1.75654942652798E-20_wp, &
      & 1.35787992059083E-27_wp, 1.53945251951349E-15_wp,-1.15138340899647E-33_wp, &
      &-2.79530089298938E-33_wp, 2.67360302817163E-31_wp, 2.98888616229728E-20_wp, &
      &-6.92374106977196E-28_wp, 2.61948127594329E-15_wp, 1.00292583863982E-31_wp, &
      & 1.08274725303022E-29_wp, 7.97402628879861E-27_wp, 1.97655003624271E-21_wp, &
      & 3.19255058891742E-27_wp, 1.73226263216247E-16_wp, 1.93797964494555E-31_wp, &
      &-1.75654943110708E-20_wp,-1.17478250486975E-27_wp,-1.53945252354520E-15_wp, &
      & 7.63042579108507E-32_wp,-6.22552788228898E-32_wp,-2.25021010193546E-31_wp, &
      & 2.98888620861064E-20_wp,-6.20974135858138E-28_wp, 2.61948131663047E-15_wp, &
      & 5.49859724849833E-32_wp, 1.03125109680276E-29_wp,-1.33570401925884E-01_wp, &
      &-1.82787972147955E-17_wp, 4.25155300877227E-01_wp, 4.54019653199472E-18_wp, &
      & 3.91883209360840E-17_wp,-1.36451437807371E-16_wp, 1.81108199345750E-01_wp, &
      & 1.21988114168630E-19_wp,-7.97975116639943E-17_wp,-2.59592601471085E-17_wp, &
      & 3.55511390567388E-17_wp, 3.22089028957633E-17_wp,-8.92061298959614E-02_wp, &
      &-1.69694570077344E-17_wp,-8.03426297896837E-17_wp, 7.97402628879861E-27_wp, &
      & 1.34115685929163E+00_wp,-1.15311752418817E-15_wp, 3.27055575235494E-01_wp, &
      & 3.77336724555385E-17_wp, 1.46139257079928E-17_wp,-3.18343013343537E-16_wp, &
      &-1.06435521305301E-01_wp, 1.71950507898743E-17_wp,-4.15164874344940E-17_wp, &
      &-1.66725294795057E-18_wp,-1.17320282908270E-16_wp, 4.49073590815545E-17_wp, &
      &-5.98906689093605E-02_wp,-1.08312508326572E-17_wp, 1.72434692616746E-16_wp, &
      &-4.56666095776420E-28_wp,-3.38149550350916E-16_wp, 2.77140751053707E-03_wp, &
      &-2.68102246439161E-16_wp,-1.27851740368253E-18_wp,-1.56769992529556E-17_wp, &
      & 2.46293497067347E-02_wp,-1.32496562140540E-16_wp, 1.00254756044942E-17_wp, &
      & 6.43443152108223E-19_wp,-1.87671991159836E-16_wp,-2.44141763724498E-17_wp, &
      & 4.19084834236004E-02_wp, 3.17807754267225E-17_wp, 3.23861724163708E-17_wp, &
      &-1.07163772987285E-17_wp, 1.97655003624271E-21_wp,-1.15311752418817E-15_wp, &
      & 2.77140747185379E-03_wp,-1.48271906593600E-16_wp,-1.46091150545397E-18_wp, &
      &-1.33561835003299E-17_wp,-2.46293497712371E-02_wp, 4.27145019848108E-17_wp, &
      &-1.64798730217797E-17_wp,-2.10392101054878E-18_wp,-1.77025307590310E-16_wp, &
      & 2.89803558844181E-17_wp, 4.19084840745453E-02_wp, 5.18274735896228E-17_wp, &
      & 7.71213850291132E-18_wp, 2.06275816562656E-17_wp,-2.50285984176088E-21_wp, &
      &-4.25155301003653E-01_wp,-1.36420649358352E-17_wp, 1.87854646548846E-01_wp, &
      &-3.05601911633501E-18_wp, 1.85106510366749E-17_wp,-3.83581862621285E-17_wp, &
      & 7.03003330580197E-02_wp, 1.13799341050549E-18_wp,-3.17748173780746E-17_wp, &
      &-1.06641003950090E-17_wp,-1.56761885797389E-18_wp,-2.04514622825458E-17_wp, &
      &-3.68339985833876E-02_wp,-4.19123223859288E-18_wp,-3.03561584380044E-17_wp, &
      & 3.19255058891742E-27_wp, 3.27055575235494E-01_wp,-1.48271906593600E-16_wp, &
      & 1.95823894368358E-01_wp, 1.21204737067721E-17_wp, 1.84292376512702E-17_wp, &
      &-3.45429873126161E-17_wp,-7.63664579683641E-02_wp, 4.64095645358721E-18_wp, &
      &-3.24375335106879E-17_wp, 9.05193648475045E-18_wp,-6.62874741202060E-17_wp, &
      &-9.44853510406184E-18_wp,-3.92154768353488E-02_wp,-3.31854425541767E-18_wp, &
      & 9.64357846109246E-17_wp,-3.69137687524134E-28_wp,-1.28989578591008E-17_wp, &
      &-4.04163230174865E-18_wp, 1.45063455819008E-17_wp, 2.77140751053706E-03_wp, &
      & 1.16292080520847E-18_wp,-2.27945384644813E-17_wp, 6.77288398570413E-18_wp, &
      & 2.46293497067347E-02_wp,-1.84280623462999E-20_wp,-4.25800023764352E-20_wp, &
      & 4.27771507037849E-18_wp,-3.67026106494832E-17_wp,-3.44909746440281E-19_wp, &
      & 4.19084834236003E-02_wp, 1.60468044849953E-18_wp, 1.73226263216247E-16_wp, &
      & 3.77336724555385E-17_wp,-1.46091150545397E-18_wp, 1.21204737067721E-17_wp, &
      & 2.77140747185378E-03_wp, 3.10068182515749E-18_wp, 2.63140382607757E-17_wp, &
      &-2.99753945246292E-18_wp,-2.46293497712371E-02_wp, 1.22079842333662E-18_wp, &
      &-9.93988107330212E-19_wp,-3.60039048649527E-18_wp,-5.36431779274770E-17_wp, &
      &-6.56723150896588E-18_wp, 4.19084840745452E-02_wp, 8.79472666047834E-19_wp, &
      & 1.64987510696543E-16_wp,-5.17357619198806E-17_wp,-1.33561836867554E-17_wp, &
      & 1.54136010738722E-17_wp, 3.10068186843677E-18_wp, 7.84267347127136E-32_wp, &
      &-1.18695687125686E-16_wp, 5.32072401654417E-18_wp, 2.75555932415949E-17_wp, &
      &-5.56835826553294E-33_wp, 9.03557987996926E-31_wp, 1.21520915923955E-31_wp, &
      &-2.01968638863391E-16_wp,-2.90343249620542E-18_wp, 4.68876822304845E-17_wp, &
      & 5.11865097234841E-32_wp, 1.93797964494555E-31_wp, 1.46139257079928E-17_wp, &
      &-1.33561835003299E-17_wp, 1.84292376512702E-17_wp, 3.10068182515749E-18_wp, &
      & 6.98993880473903E-32_wp, 1.18695687436541E-16_wp,-7.61620111599659E-18_wp, &
      &-2.75555933137609E-17_wp, 8.19504835423912E-33_wp, 8.53215352189635E-31_wp, &
      &-1.49796527583809E-31_wp,-2.01968642000476E-16_wp,-3.80460525829567E-18_wp, &
      & 4.68876829587692E-17_wp,-8.95816636479575E-32_wp, 1.84601943108622E-31_wp, &
      &-8.39801861579282E-17_wp,-2.46293501150133E-02_wp,-7.89082478289947E-17_wp, &
      & 2.49944922548076E-17_wp, 1.39320652731803E-16_wp,-2.18879711744274E-01_wp, &
      & 2.63679918315433E-17_wp, 2.85084061341327E-17_wp,-5.71824483123694E-18_wp, &
      & 1.66783093409513E-15_wp, 2.16967477846744E-16_wp,-3.72438447649679E-01_wp, &
      & 2.61276534738854E-16_wp,-8.20267603302476E-17_wp, 9.52358711057317E-17_wp, &
      &-1.75654943110708E-20_wp,-3.18343013343537E-16_wp,-2.46293497712371E-02_wp, &
      &-3.45429873126161E-17_wp, 2.63140382607757E-17_wp, 1.18695687436541E-16_wp, &
      & 2.18879712317504E-01_wp,-7.28625878311845E-17_wp, 2.10980540486193E-17_wp, &
      & 1.86974333388802E-17_wp, 1.57321442742816E-15_wp,-2.57546870614022E-16_wp, &
      &-3.72438453434591E-01_wp,-2.48370488048218E-16_wp, 1.36648660865405E-16_wp, &
      &-1.83316213406574E-16_wp, 2.22427813662622E-20_wp, 1.81108201463904E-01_wp, &
      & 1.43755402599569E-17_wp,-7.03003339238226E-02_wp, 2.71256570939184E-18_wp, &
      &-6.99902444140602E-18_wp, 9.24488963624558E-17_wp,-2.57240327880902E-02_wp, &
      & 1.07483602006170E-17_wp, 1.16816478879813E-17_wp, 3.94236060996094E-18_wp, &
      & 1.62486129743607E-18_wp, 1.44863137958540E-16_wp, 1.36291376558167E-02_wp, &
      & 2.04943051581071E-17_wp, 1.10510160861017E-17_wp,-1.17478250486975E-27_wp, &
      &-1.06435521305301E-01_wp, 4.27145019848108E-17_wp,-7.63664579683641E-02_wp, &
      &-2.99753945246292E-18_wp,-7.61620111599659E-18_wp,-7.28625878311845E-17_wp, &
      & 3.03415153930987E-02_wp,-1.27984661650380E-17_wp, 1.29859766419532E-17_wp, &
      &-3.97574070488552E-18_wp, 2.56749466204716E-17_wp, 1.40421498114251E-16_wp, &
      & 1.54418981293995E-02_wp, 2.02782674158168E-17_wp,-3.73055872161586E-17_wp, &
      & 1.48274390135123E-28_wp,-9.72496982710943E-18_wp, 6.93889390390723E-18_wp, &
      & 1.97389772270935E-18_wp,-2.46293501150133E-02_wp,-1.03348149121388E-17_wp, &
      &-5.55111512312578E-17_wp,-6.64332971457817E-18_wp,-2.18879711744274E-01_wp, &
      & 1.63769203100801E-19_wp, 3.78406200620225E-19_wp,-3.80158247244567E-17_wp, &
      &-1.11022302462516E-16_wp,-2.38114673055874E-17_wp,-3.72438447649679E-01_wp, &
      &-1.42607092022899E-17_wp,-1.53945252354520E-15_wp, 1.71950507898743E-17_wp, &
      &-1.64798730217797E-17_wp, 4.64095645358721E-18_wp,-2.46293497712371E-02_wp, &
      &-2.75555933137609E-17_wp, 2.10980540486193E-17_wp,-1.27984661650380E-17_wp, &
      & 2.18879712317503E-01_wp,-1.08491702046326E-17_wp, 8.83351907379951E-18_wp, &
      & 3.19964774236576E-17_wp, 4.28968317636652E-17_wp, 3.70251628283329E-17_wp, &
      &-3.72438453434591E-01_wp,-7.81582648034288E-18_wp,-1.46623517115395E-15_wp, &
      & 7.96062373381183E-17_wp,-2.10392103991529E-18_wp,-2.93432673339620E-17_wp, &
      & 1.22079844037651E-18_wp, 9.47901645394113E-33_wp,-1.86974332899130E-17_wp, &
      &-1.06306472855005E-17_wp, 1.08491701762195E-17_wp, 4.34114355013797E-33_wp, &
      & 1.44089616294457E-31_wp, 2.12859207917019E-32_wp,-3.18149314709649E-17_wp, &
      & 5.66049568925336E-18_wp, 1.84605876283281E-17_wp, 1.33985466103712E-32_wp, &
      & 7.63042579108507E-32_wp,-4.15164874344940E-17_wp,-2.10392101054878E-18_wp, &
      &-3.24375335106879E-17_wp, 1.22079842333662E-18_wp, 8.19504835423912E-33_wp, &
      & 1.86974333388802E-17_wp, 1.29859766419532E-17_wp,-1.08491702046326E-17_wp, &
      & 7.70971343637110E-33_wp, 1.32184635135000E-31_wp,-1.27114553952805E-32_wp, &
      &-3.18149319651314E-17_wp, 6.58516571499225E-18_wp, 1.84605879150679E-17_wp, &
      &-3.10652176789749E-32_wp, 7.26785075070276E-32_wp,-3.18260052314977E-17_wp, &
      &-1.77025310061228E-16_wp, 6.33107175306980E-18_wp,-9.93988121204300E-19_wp, &
      & 1.00164156833194E-30_wp,-1.57321442330803E-15_wp, 1.90464148457209E-18_wp, &
      &-8.83351905066482E-18_wp,-4.19980981347647E-32_wp, 1.19873549172850E-29_wp, &
      & 1.55705742722027E-30_wp,-2.67692941007473E-15_wp,-1.11800158227415E-18_wp, &
      &-1.50308226207710E-17_wp, 6.83162195812944E-31_wp,-6.22552788228898E-32_wp, &
      &-1.66725294795057E-18_wp,-1.77025307590310E-16_wp, 9.05193648475045E-18_wp, &
      &-9.93988107330212E-19_wp, 8.53215352189635E-31_wp, 1.57321442742816E-15_wp, &
      &-3.97574070488552E-18_wp, 8.83351907379951E-18_wp, 1.32184635135000E-31_wp, &
      & 1.13087266179515E-29_wp,-1.85277078485464E-30_wp,-2.67692945165423E-15_wp, &
      &-1.93108662221892E-18_wp,-1.50308228542362E-17_wp,-1.31369693028273E-30_wp, &
      &-5.90142574190684E-32_wp, 1.39123540364297E-16_wp, 2.89803562889256E-17_wp, &
      &-6.45162606238304E-17_wp,-3.60039053674948E-18_wp,-1.71778478447135E-31_wp, &
      & 2.57546869939526E-16_wp,-2.43267135760346E-17_wp,-3.19964773398613E-17_wp, &
      & 1.77305803405670E-32_wp,-1.95873601060336E-30_wp,-2.60640833892739E-31_wp, &
      & 4.38233199746618E-16_wp, 1.26987491254611E-17_wp,-5.44441431129022E-17_wp, &
      &-1.03622630068020E-31_wp,-2.25021010193546E-31_wp,-1.17320282908270E-16_wp, &
      & 2.89803558844181E-17_wp,-6.62874741202060E-17_wp,-3.60039048649527E-18_wp, &
      &-1.49796527583809E-31_wp,-2.57546870614022E-16_wp, 2.56749466204716E-17_wp, &
      & 3.19964774236576E-17_wp,-1.27114553952805E-32_wp,-1.85277078485464E-30_wp, &
      & 3.30215928631260E-31_wp, 4.38233206553490E-16_wp, 1.32280472591487E-17_wp, &
      &-5.44441439585579E-17_wp, 1.81819591139079E-31_wp,-2.14364723955075E-31_wp, &
      & 1.33855770960173E-17_wp, 4.19084846595035E-02_wp, 1.51647595042051E-17_wp, &
      &-5.29495672833062E-17_wp,-2.37063398364039E-16_wp, 3.72438452459204E-01_wp, &
      &-1.00035099535184E-16_wp,-1.36262959273222E-16_wp, 9.72997560512568E-18_wp, &
      &-2.83792575889221E-15_wp,-3.69184658729941E-16_wp, 6.33728900561671E-01_wp, &
      &-4.18412246384521E-16_wp,-7.26596145091637E-18_wp,-1.62050197209065E-16_wp, &
      & 2.98888620861064E-20_wp, 4.49073590815545E-17_wp, 4.19084840745453E-02_wp, &
      &-9.44853510406184E-18_wp,-5.36431779274770E-17_wp,-2.01968642000476E-16_wp, &
      &-3.72438453434591E-01_wp, 1.40421498114251E-16_wp, 4.28968317636652E-17_wp, &
      &-3.18149319651314E-17_wp,-2.67692945165423E-15_wp, 4.38233206553490E-16_wp, &
      & 6.33728910405086E-01_wp, 4.33582422825006E-16_wp,-3.97493879491129E-16_wp, &
      & 3.11924784109697E-16_wp,-3.78475784920443E-20_wp, 8.92061338128334E-02_wp, &
      & 3.12830876064525E-17_wp,-3.68340001692192E-02_wp,-3.58333013113701E-18_wp, &
      &-3.64859041169069E-18_wp, 2.60621297069682E-16_wp,-1.36291380722532E-02_wp, &
      &-3.81283996669708E-17_wp, 6.17472181633275E-18_wp, 2.07812411870623E-18_wp, &
      & 5.83066485252031E-19_wp, 4.35825887214287E-16_wp, 7.18110568545843E-03_wp, &
      &-6.36998672946740E-17_wp, 5.87007225823383E-18_wp,-6.20974135858138E-28_wp, &
      &-5.98906689093605E-02_wp, 5.18274735896228E-17_wp,-3.92154768353488E-02_wp, &
      &-6.56723150896588E-18_wp,-3.80460525829567E-18_wp,-2.48370488048218E-16_wp, &
      & 1.54418981293995E-02_wp, 3.70251628283329E-17_wp, 6.58516571499225E-18_wp, &
      &-1.93108662221892E-18_wp, 1.32280472591487E-17_wp, 4.33582422825006E-16_wp, &
      & 7.89277053508258E-03_wp,-6.38416276888399E-17_wp,-1.92319359732299E-17_wp, &
      & 7.47963975966199E-29_wp, 1.44260904550692E-17_wp,-3.17171050218737E-17_wp, &
      & 2.54424582926323E-18_wp, 4.19084846595034E-02_wp, 1.75853780218166E-17_wp, &
      &-9.34514695314122E-17_wp, 1.38115301449283E-17_wp, 3.72438452459204E-01_wp, &
      &-2.78664240176825E-19_wp,-6.43883430933571E-19_wp, 6.46864655317119E-17_wp, &
      &-1.36729465914179E-16_wp, 3.92801248949416E-17_wp, 6.33728900561672E-01_wp, &
      & 2.42655494378433E-17_wp, 2.61948131663047E-15_wp,-1.08312508326572E-17_wp, &
      & 7.71213850291132E-18_wp,-3.31854425541767E-18_wp, 4.19084840745452E-02_wp, &
      & 4.68876829587692E-17_wp, 1.36648660865405E-16_wp, 2.02782674158168E-17_wp, &
      &-3.72438453434591E-01_wp, 1.84605879150679E-17_wp,-1.50308228542362E-17_wp, &
      &-5.44441439585579E-17_wp,-3.97493879491129E-16_wp,-6.38416276888399E-17_wp, &
      & 6.33728910405087E-01_wp, 1.32991509164151E-17_wp, 2.49489710002836E-15_wp, &
      &-2.01125186749760E-16_wp, 2.06275819441853E-17_wp, 9.41050782841875E-17_wp, &
      & 8.79472678323498E-19_wp,-1.07080703791633E-31_wp, 1.83316212926483E-16_wp, &
      & 3.55314728153116E-17_wp, 7.81582645987379E-18_wp,-1.12469707088764E-32_wp, &
      &-1.40222409804720E-30_wp,-1.80583133059065E-31_wp, 3.11924779264715E-16_wp, &
      &-1.85354424939333E-17_wp, 1.32991507098458E-17_wp,-9.46260833251232E-32_wp, &
      & 5.49859724849833E-32_wp, 1.72434692616746E-16_wp, 2.06275816562656E-17_wp, &
      & 9.64357846109246E-17_wp, 8.79472666047834E-19_wp,-8.95816636479575E-32_wp, &
      &-1.83316213406574E-16_wp,-3.73055872161586E-17_wp,-7.81582648034288E-18_wp, &
      &-3.10652176789749E-32_wp,-1.31369693028273E-30_wp, 1.81819591139079E-31_wp, &
      & 3.11924784109697E-16_wp,-1.92319359732299E-17_wp, 1.32991509164151E-17_wp, &
      & 2.01463782277123E-31_wp, 5.23381637033705E-32_wp, 9.17410126990192E-28_wp, &
      &-2.50285987685287E-21_wp,-3.31686667557185E-28_wp, 1.64987512999437E-16_wp, &
      & 6.92451933067693E-32_wp,-2.22427813078371E-20_wp,-1.19650353049022E-28_wp, &
      & 1.46623516731399E-15_wp,-1.09764104060678E-33_wp,-2.36538725322292E-33_wp, &
      & 2.54683107179597E-31_wp,-3.78475779031596E-20_wp, 6.40177057547808E-29_wp, &
      & 2.49489706127629E-15_wp, 9.55395613748457E-32_wp, 1.03125109680276E-29_wp, &
      &-4.56666095776420E-28_wp,-2.50285984176088E-21_wp,-3.69137687524134E-28_wp, &
      & 1.64987510696543E-16_wp, 1.84601943108622E-31_wp, 2.22427813662622E-20_wp, &
      & 1.48274390135123E-28_wp,-1.46623517115395E-15_wp, 7.26785075070276E-32_wp, &
      &-5.90142574190684E-32_wp,-2.14364723955075E-31_wp,-3.78475784920443E-20_wp, &
      & 7.47963975966199E-29_wp, 2.49489710002836E-15_wp, 5.23381637033705E-32_wp, &
      & 9.82204131603147E-30_wp], shape(density))

   call get_structure(mol, "f-block", "Ce2")
   call test_numpot(error, mol, density, qsh, make_exchange_gxtb, &
      & thr_in=thr1)

end subroutine test_p_fock_ce2


subroutine test_op_g_fock_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: density(2, 2, 1) = reshape([&
      & 5.93683766916992E-1_wp, 5.93683766916992E-1_wp, 5.93683766916992E-1_wp, &
      & 5.93683766916992E-1_wp], shape(density))

   call get_structure(mol, "MB16-43", "H2")
   call test_num_op_grad(error, mol, density, make_exchange_gxtb, thr_in=thr1)

end subroutine test_op_g_fock_h2


subroutine test_op_g_fock_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: density(5, 5, 1) = reshape([&
      & 7.43138968868805E-02_wp, 6.30732585479418E-45_wp, 1.15038033099931E-01_wp, &
      & 0.00000000000000E+00_wp, 2.77067464359636E-01_wp, 6.30732585479418E-45_wp, &
      & 5.35328668056679E-88_wp, 9.76375066914250E-45_wp, 0.00000000000000E+00_wp, &
      & 2.35158544321515E-44_wp, 1.15038033099931E-01_wp, 9.76375066914250E-45_wp, &
      & 1.78079062111965E-01_wp, 0.00000000000000E+00_wp, 4.28900884910328E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 2.77067464359636E-01_wp, &
      & 2.35158544321515E-44_wp, 4.28900884910328E-01_wp, 0.00000000000000E+00_wp, &
      & 1.03300167293785E+00_wp], shape(density))

   call get_structure(mol, "MB16-43", "LiH")
   call test_num_op_grad(error, mol, density, make_exchange_gxtb, thr_in=thr1)

end subroutine test_op_g_fock_lih


subroutine test_op_g_fock_no(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: density(8, 8, 2) = reshape([&
      & 9.42009046415958E-01_wp,-1.97869914210805E-16_wp,-3.13530611389513E-01_wp, &
      & 1.57008624511250E-16_wp,-1.91210322106419E-01_wp, 1.16548396645701E-16_wp, &
      &-2.58330200259595E-02_wp,-1.75070805017504E-16_wp,-1.97869914210805E-16_wp, &
      & 7.16403612559420E-01_wp, 6.75092136044200E-17_wp, 4.11552712525306E-01_wp, &
      & 1.21088275692219E-16_wp,-3.76653483713363E-03_wp,-1.64250151646619E-17_wp, &
      &-3.38085910709977E-01_wp,-3.13530611389513E-01_wp, 6.75092136044200E-17_wp, &
      & 4.49934652749990E-01_wp,-8.24568998262442E-17_wp,-1.71924645027865E-03_wp, &
      &-1.23858164254808E-17_wp,-4.07076133805151E-01_wp, 1.00925254078120E-16_wp, &
      & 1.57008624511250E-16_wp, 4.11552712525306E-01_wp,-8.24568998262442E-17_wp, &
      & 6.31376987893123E-01_wp,-1.61242008793196E-16_wp,-3.38085910709977E-01_wp, &
      &-6.74819946100851E-18_wp, 6.60818781608317E-02_wp,-1.91210322106419E-01_wp, &
      & 1.21088275692219E-16_wp,-1.71924645027865E-03_wp,-1.61242008793196E-16_wp, &
      & 1.01808686826266E+00_wp, 4.52238546687169E-17_wp, 2.62723878284716E-01_wp, &
      &-1.97376402183765E-17_wp, 1.16548396645701E-16_wp,-3.76653483713363E-03_wp, &
      &-1.23858164254808E-17_wp,-3.38085910709977E-01_wp, 4.52238546687169E-17_wp, &
      & 8.37111178906109E-01_wp, 3.74456682512930E-17_wp, 2.77733761780434E-01_wp, &
      &-2.58330200259595E-02_wp,-1.64250151646619E-17_wp,-4.07076133805151E-01_wp, &
      &-6.74819946100851E-18_wp, 2.62723878284716E-01_wp, 3.74456682512930E-17_wp, &
      & 5.33778495172291E-01_wp, 2.29800571514756E-17_wp,-1.75070805017504E-16_wp, &
      &-3.38085910709977E-01_wp, 1.00925254078120E-16_wp, 6.60818781608317E-02_wp, &
      &-1.97376402183765E-17_wp, 2.77733761780434E-01_wp, 2.29800571514756E-17_wp, &
      & 7.79731495571835E-01_wp, 9.54305653170375E-01_wp,-1.73912460120364E-16_wp, &
      &-3.11307650277130E-01_wp, 1.94985952807570E-16_wp,-2.03773626333088E-01_wp, &
      & 1.87683417049811E-17_wp,-1.70067070427977E-02_wp, 2.26960265829082E-17_wp, &
      &-1.73912460120364E-16_wp, 2.40284367737361E-01_wp, 1.07958284666329E-16_wp, &
      &-1.87742209952433E-05_wp, 1.14092931051935E-16_wp, 3.64460563073953E-01_wp, &
      &-4.86861943671765E-17_wp,-6.78042890246151E-06_wp,-3.11307650277130E-01_wp, &
      & 1.07958284666329E-16_wp, 4.38731609559674E-01_wp,-1.59063207389752E-16_wp, &
      &-4.10640194692124E-03_wp,-3.43859837540706E-17_wp,-4.09548671902707E-01_wp, &
      &-2.10949694797432E-17_wp, 1.94985952807570E-16_wp,-1.87742209952433E-05_wp, &
      &-1.59063207389752E-16_wp, 2.40288246483859E-01_wp,-6.87541757092717E-17_wp, &
      &-6.78042890286945E-06_wp, 6.84810692055402E-17_wp, 3.64461963907834E-01_wp, &
      &-2.03773626333088E-01_wp, 1.14092931051935E-16_wp,-4.10640194692124E-03_wp, &
      &-6.87541757092717E-17_wp, 1.03092149514981E+00_wp,-1.10878908539341E-16_wp, &
      & 2.53665476333654E-01_wp, 5.86026180268729E-17_wp, 1.87683417049811E-17_wp, &
      & 3.64460563073953E-01_wp,-3.43859837540706E-17_wp,-6.78042890286945E-06_wp, &
      &-1.10878908539341E-16_wp, 5.52809588728990E-01_wp, 1.03373764663983E-17_wp, &
      & 2.26235036331901E-05_wp,-1.70067070427977E-02_wp,-4.86861943671765E-17_wp, &
      &-4.09548671902707E-01_wp, 6.84810692055402E-17_wp, 2.53665476333654E-01_wp, &
      & 1.03373764663983E-17_wp, 5.38687781288217E-01_wp,-3.05073307416719E-17_wp, &
      & 2.26960265829082E-17_wp,-6.78042890246151E-06_wp,-2.10949694797432E-17_wp, &
      & 3.64461963907834E-01_wp, 5.86026180268729E-17_wp, 2.26235036331901E-05_wp, &
      &-3.05073307416719E-17_wp, 5.52804914722249E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "NO")
   call test_num_op_grad(error, mol, density, make_exchange_gxtb, thr_in=thr1)

end subroutine test_op_g_fock_no


subroutine test_op_g_fock_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: density(18, 18, 1) = reshape([&
      & 2.01512274410727E+00_wp,-2.98954363036077E-16_wp,-3.99929219814917E-01_wp, &
      &-3.60921385847533E-16_wp, 1.66928202631975E-04_wp, 2.96505862833536E-17_wp, &
      &-3.51352512626524E-02_wp,-2.26041467075836E-17_wp,-1.25121106619603E-05_wp, &
      &-2.61461503009251E-01_wp, 2.00709987824803E-16_wp, 8.46120257485226E-02_wp, &
      & 4.15182726310734E-16_wp,-2.96659169389449E-05_wp, 3.71826627538686E-17_wp, &
      & 1.18677041085955E-02_wp, 5.96815828699915E-17_wp, 2.22361008940176E-06_wp, &
      &-2.98954363036077E-16_wp, 1.45486170073001E+00_wp,-9.00635884019608E-16_wp, &
      & 5.97801199387253E-01_wp, 1.67776446953954E-17_wp,-2.61299185132976E-03_wp, &
      & 1.68298746105631E-16_wp,-4.41320346165751E-02_wp,-2.78389072323402E-18_wp, &
      & 4.34685388862911E-16_wp, 1.61265924256514E-01_wp,-6.14873875231672E-16_wp, &
      &-6.02461837620382E-01_wp, 4.43319518778254E-17_wp,-1.03162491097761E-01_wp, &
      &-2.04398675254949E-16_wp,-5.40117632485737E-02_wp, 1.09406633163272E-16_wp, &
      &-3.99929219814917E-01_wp,-9.00635884019608E-16_wp, 8.16076844347620E-01_wp, &
      &-1.55848848678421E-16_wp,-2.73809977947621E-05_wp, 4.89930425005989E-17_wp, &
      & 5.39351481156965E-02_wp, 2.47586597786769E-16_wp, 2.05234387570902E-06_wp, &
      &-8.46120257485189E-02_wp, 6.66084959945643E-16_wp,-7.72403989719108E-01_wp, &
      & 5.71110633457431E-16_wp,-1.51824865425160E-07_wp, 6.80560093611612E-17_wp, &
      & 4.74250266390748E-02_wp,-1.14353946762662E-17_wp, 1.13800395663470E-08_wp, &
      &-3.60921385847533E-16_wp, 5.97801199387253E-01_wp,-1.55848848678421E-16_wp, &
      & 1.36524527111739E+00_wp, 1.53974719339088E-17_wp,-4.41320346165748E-02_wp, &
      & 1.44159180113409E-16_wp, 4.00284528304191E-03_wp, 4.02705177739493E-18_wp, &
      & 8.69046775550039E-17_wp,-6.02461837620375E-01_wp,-7.97303127069371E-16_wp, &
      & 2.51581030550614E-01_wp,-5.41615397023848E-17_wp,-5.40117632485737E-02_wp, &
      & 4.07523509225128E-18_wp,-9.50655829878735E-02_wp, 9.48838887830889E-17_wp, &
      & 1.66928202631975E-04_wp, 1.67776446953954E-17_wp,-2.73809977947621E-05_wp, &
      & 1.53974719339088E-17_wp, 1.38974440043868E-08_wp,-6.94704531194593E-19_wp, &
      &-2.60092206439580E-06_wp,-6.56220312997718E-19_wp,-1.04168351758665E-09_wp, &
      &-2.96659169388496E-05_wp,-7.73908741086395E-18_wp, 1.51824865387540E-07_wp, &
      &-7.33055884782161E-18_wp,-3.07941415440952E-09_wp,-1.30780677462207E-18_wp, &
      & 1.45801168257227E-06_wp,-1.20266021436374E-18_wp, 2.30817621363979E-10_wp, &
      & 2.96505862833536E-17_wp,-2.61299185132976E-03_wp, 4.89930425005989E-17_wp, &
      &-4.41320346165748E-02_wp,-6.94704531194593E-19_wp, 7.48444822726994E-03_wp, &
      &-1.84444900024127E-17_wp, 4.63153745136131E-03_wp, 2.56316174731192E-19_wp, &
      &-4.56723208050267E-17_wp, 1.03162491097758E-01_wp, 4.43388496380181E-17_wp, &
      & 5.40117632485682E-02_wp, 2.70455166704829E-18_wp, 1.16466225731597E-03_wp, &
      & 1.89516760801754E-17_wp, 3.39354102094470E-03_wp,-1.17177365770997E-17_wp, &
      &-3.51352512626524E-02_wp, 1.68298746105631E-16_wp, 5.39351481156965E-02_wp, &
      & 1.44159180113409E-16_wp,-2.60092206439580E-06_wp,-1.84444900024127E-17_wp, &
      & 3.73731424127323E-03_wp,-7.83051079850619E-18_wp, 1.94952226013385E-07_wp, &
      & 1.18677041085919E-02_wp,-2.17279393388892E-16_wp,-4.74250266390769E-02_wp, &
      &-2.59398305030398E-16_wp, 1.45801168257080E-06_wp,-1.41001874145128E-17_wp, &
      & 2.76687889875095E-03_wp,-1.47129262760623E-17_wp,-1.09285328835862E-07_wp, &
      &-2.26041467075836E-17_wp,-4.41320346165751E-02_wp, 2.47586597786769E-16_wp, &
      & 4.00284528304191E-03_wp,-6.56220312997718E-19_wp, 4.63153745136131E-03_wp, &
      &-7.83051079850619E-18_wp, 6.79013404651283E-03_wp, 6.04554882466166E-19_wp, &
      &-1.74543578973199E-16_wp, 5.40117632485687E-02_wp,-1.45442246694960E-16_wp, &
      & 9.50655829878701E-02_wp,-2.78339838142217E-18_wp, 3.39354102094469E-03_wp, &
      & 3.97034075802654E-17_wp, 6.55936228742013E-04_wp,-1.12393067366785E-17_wp, &
      &-1.25121106619603E-05_wp,-2.78389072323402E-18_wp, 2.05234387570902E-06_wp, &
      & 4.02705177739493E-18_wp,-1.04168351758665E-09_wp, 2.56316174731192E-19_wp, &
      & 1.94952226013385E-07_wp, 6.04554882466166E-19_wp, 7.80794332014712E-11_wp, &
      & 2.22361008999868E-06_wp, 2.61284439827709E-18_wp,-1.13800387946783E-08_wp, &
      & 9.03292745688553E-18_wp, 2.30817621408582E-10_wp, 1.84803242584056E-19_wp, &
      &-1.09285328887014E-07_wp,-1.95198444371785E-19_wp,-1.73009448131981E-11_wp, &
      &-2.61461503009251E-01_wp, 4.34685388862911E-16_wp,-8.46120257485189E-02_wp, &
      & 8.69046775550039E-17_wp,-2.96659169388496E-05_wp,-4.56723208050267E-17_wp, &
      & 1.18677041085919E-02_wp,-1.74543578973199E-16_wp, 2.22361008999868E-06_wp, &
      & 2.01512274410727E+00_wp,-2.95337324834156E-16_wp, 3.99929219814911E-01_wp, &
      &-2.60805282861184E-16_wp, 1.66928202632090E-04_wp,-3.91687455452402E-17_wp, &
      &-3.51352512626490E-02_wp, 3.43718096863094E-18_wp,-1.25121106621160E-05_wp, &
      & 2.00709987824803E-16_wp, 1.61265924256514E-01_wp, 6.66084959945643E-16_wp, &
      &-6.02461837620375E-01_wp,-7.73908741086395E-18_wp, 1.03162491097758E-01_wp, &
      &-2.17279393388892E-16_wp, 5.40117632485687E-02_wp, 2.61284439827709E-18_wp, &
      &-2.95337324834156E-16_wp, 1.45486170073003E+00_wp, 4.61646074650480E-16_wp, &
      & 5.97801199387235E-01_wp, 4.81747951661595E-17_wp, 2.61299185133292E-03_wp, &
      & 2.21800860913915E-16_wp, 4.41320346165780E-02_wp,-1.45847563550088E-16_wp, &
      & 8.46120257485226E-02_wp,-6.14873875231672E-16_wp,-7.72403989719108E-01_wp, &
      &-7.97303127069371E-16_wp, 1.51824865387540E-07_wp, 4.43388496380181E-17_wp, &
      &-4.74250266390769E-02_wp,-1.45442246694960E-16_wp,-1.13800387946783E-08_wp, &
      & 3.99929219814911E-01_wp, 4.61646074650480E-16_wp, 8.16076844347641E-01_wp, &
      & 9.30644287891831E-16_wp, 2.73809977948286E-05_wp, 4.76651490748586E-17_wp, &
      &-5.39351481156958E-02_wp, 9.39184312631712E-17_wp,-2.05234387641953E-06_wp, &
      & 4.15182726310734E-16_wp,-6.02461837620382E-01_wp, 5.71110633457431E-16_wp, &
      & 2.51581030550614E-01_wp,-7.33055884782161E-18_wp, 5.40117632485682E-02_wp, &
      &-2.59398305030398E-16_wp, 9.50655829878701E-02_wp, 9.03292745688553E-18_wp, &
      &-2.60805282861184E-16_wp, 5.97801199387235E-01_wp, 9.30644287891831E-16_wp, &
      & 1.36524527111740E+00_wp,-5.05746530671741E-17_wp, 4.41320346165780E-02_wp, &
      & 3.53371679850913E-16_wp,-4.00284528303917E-03_wp,-1.41955353109642E-16_wp, &
      &-2.96659169389449E-05_wp, 4.43319518778254E-17_wp,-1.51824865425160E-07_wp, &
      &-5.41615397023848E-17_wp,-3.07941415440952E-09_wp, 2.70455166704829E-18_wp, &
      & 1.45801168257080E-06_wp,-2.78339838142217E-18_wp, 2.30817621408582E-10_wp, &
      & 1.66928202632090E-04_wp, 4.81747951661595E-17_wp, 2.73809977948286E-05_wp, &
      &-5.05746530671741E-17_wp, 1.38974440044059E-08_wp,-2.39255576630622E-18_wp, &
      &-2.60092206440109E-06_wp, 3.06620883893745E-18_wp,-1.04168351759322E-09_wp, &
      & 3.71826627538686E-17_wp,-1.03162491097761E-01_wp, 6.80560093611612E-17_wp, &
      &-5.40117632485737E-02_wp,-1.30780677462207E-18_wp, 1.16466225731597E-03_wp, &
      &-1.41001874145128E-17_wp, 3.39354102094469E-03_wp, 1.84803242584056E-19_wp, &
      &-3.91687455452402E-17_wp, 2.61299185133292E-03_wp, 4.76651490748586E-17_wp, &
      & 4.41320346165780E-02_wp,-2.39255576630622E-18_wp, 7.48444822727060E-03_wp, &
      & 1.59370810633493E-17_wp, 4.63153745136205E-03_wp,-9.15386829425696E-18_wp, &
      & 1.18677041085955E-02_wp,-2.04398675254949E-16_wp, 4.74250266390748E-02_wp, &
      & 4.07523509225128E-18_wp, 1.45801168257227E-06_wp, 1.89516760801754E-17_wp, &
      & 2.76687889875095E-03_wp, 3.97034075802654E-17_wp,-1.09285328887014E-07_wp, &
      &-3.51352512626490E-02_wp, 2.21800860913915E-16_wp,-5.39351481156958E-02_wp, &
      & 3.53371679850913E-16_wp,-2.60092206440109E-06_wp, 1.59370810633493E-17_wp, &
      & 3.73731424127310E-03_wp, 1.71392595868698E-18_wp, 1.94952226055385E-07_wp, &
      & 5.96815828699915E-17_wp,-5.40117632485737E-02_wp,-1.14353946762662E-17_wp, &
      &-9.50655829878735E-02_wp,-1.20266021436374E-18_wp, 3.39354102094470E-03_wp, &
      &-1.47129262760623E-17_wp, 6.55936228742013E-04_wp,-1.95198444371785E-19_wp, &
      & 3.43718096863094E-18_wp, 4.41320346165780E-02_wp, 9.39184312631712E-17_wp, &
      &-4.00284528303917E-03_wp, 3.06620883893745E-18_wp, 4.63153745136205E-03_wp, &
      & 1.71392595868698E-18_wp, 6.79013404651345E-03_wp,-8.12591467858445E-18_wp, &
      & 2.22361008940176E-06_wp, 1.09406633163272E-16_wp, 1.13800395663470E-08_wp, &
      & 9.48838887830889E-17_wp, 2.30817621363979E-10_wp,-1.17177365770997E-17_wp, &
      &-1.09285328835862E-07_wp,-1.12393067366785E-17_wp,-1.73009448131981E-11_wp, &
      &-1.25121106621160E-05_wp,-1.45847563550088E-16_wp,-2.05234387641953E-06_wp, &
      &-1.41955353109642E-16_wp,-1.04168351759322E-09_wp,-9.15386829425696E-18_wp, &
      & 1.94952226055385E-07_wp,-8.12591467858445E-18_wp, 7.80794332023474E-11_wp],&
      & shape(density))

   call get_structure(mol, "MB16-43", "S2")
   call test_num_op_grad(error, mol, density, make_exchange_gxtb, thr_in=thr1)

end subroutine test_op_g_fock_s2


subroutine test_op_g_fock_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: density(43, 43, 2) = reshape([&
      & 6.91475853693089E-03_wp,-4.66439474604429E-04_wp, 7.43124629337322E-04_wp, &
      &-3.09904488020287E-04_wp,-2.72087829479343E-03_wp, 6.74488694982622E-03_wp, &
      &-5.28266720158092E-03_wp, 4.80689218345537E-03_wp, 9.22601160315454E-04_wp, &
      &-2.71966921747329E-03_wp,-1.01894320857934E-02_wp,-2.27063129355292E-03_wp, &
      &-2.82181425186610E-03_wp,-2.14613394982406E-03_wp, 3.61217559598682E-03_wp, &
      & 1.12028617140403E-02_wp,-8.64033744988582E-04_wp, 3.10371374558505E-02_wp, &
      & 2.30333249958885E-02_wp, 2.04795146556320E-02_wp, 1.17566904494576E-03_wp, &
      & 1.12383332750122E-03_wp,-3.37161261858860E-04_wp, 7.55695280205111E-04_wp, &
      &-4.91653601591487E-04_wp,-1.02695166967632E-03_wp,-3.45289149212074E-02_wp, &
      &-1.34121866165070E-02_wp, 2.17842958084033E-02_wp,-1.32239449246865E-03_wp, &
      & 9.23763188396611E-04_wp,-5.87060101043125E-04_wp,-6.06266471529292E-04_wp, &
      &-5.25957119782735E-04_wp,-9.44887800911054E-04_wp, 9.71328579482040E-03_wp, &
      &-1.32264514235801E-02_wp,-4.00506954181809E-02_wp,-7.30975427202173E-04_wp, &
      &-2.86085329171831E-04_wp,-5.83829373776639E-04_wp, 1.06473924707411E-03_wp, &
      & 1.20864032142143E-03_wp,-4.66439474604429E-04_wp, 2.37196189260497E-03_wp, &
      & 3.76291798637178E-04_wp,-1.89524702709497E-04_wp,-2.70383567603876E-03_wp, &
      & 7.51687695834723E-04_wp, 1.57852557800278E-05_wp,-1.70421269428447E-03_wp, &
      &-1.24711040739595E-04_wp, 7.89473512446235E-04_wp,-3.57098061462513E-04_wp, &
      & 4.01348859449634E-03_wp,-3.38878489220265E-03_wp, 6.34879336465873E-03_wp, &
      &-2.30645681571837E-03_wp,-6.70882068372886E-04_wp, 6.96944723684011E-03_wp, &
      &-1.08436093353111E-02_wp,-1.77208852887962E-02_wp,-1.93565751770493E-02_wp, &
      &-8.75401591255748E-04_wp,-5.65296230455721E-04_wp, 1.94865564588492E-04_wp, &
      &-6.09542770974166E-04_wp,-6.86637147630286E-05_wp,-7.77156635669540E-03_wp, &
      &-1.12023415897657E-02_wp,-1.78971328951955E-02_wp, 2.23980508202680E-02_wp, &
      &-8.58573131873999E-04_wp, 6.50886246445002E-04_wp,-1.23589021151666E-04_wp, &
      &-7.14602097206055E-04_wp, 1.61848868703983E-04_wp, 2.47681579670218E-03_wp, &
      & 1.51505155319696E-02_wp, 4.13303227393656E-03_wp, 1.19937863231962E-02_wp, &
      &-4.00974483671832E-04_wp,-1.44661987627843E-04_wp,-2.20905355380273E-05_wp, &
      &-4.45102163842455E-04_wp,-5.21718641131263E-04_wp, 7.43124629337322E-04_wp, &
      & 3.76291798637178E-04_wp, 1.92419514812543E-03_wp, 2.43411943306908E-04_wp, &
      &-1.38558534276645E-03_wp,-1.29700970258115E-03_wp, 3.20274034784215E-04_wp, &
      &-6.00771636389316E-04_wp, 4.53126654291574E-04_wp, 6.89363513188489E-04_wp, &
      & 2.85458472253175E-03_wp,-8.42769400826849E-03_wp,-3.49285659189717E-03_wp, &
      &-8.73740458821759E-03_wp,-2.03180597844227E-03_wp, 3.57585309496864E-03_wp, &
      & 4.50241444890764E-03_wp,-1.94202483147287E-02_wp, 1.08132873675669E-02_wp, &
      &-1.31156713246053E-02_wp,-6.35109734800467E-04_wp,-1.41358099308255E-04_wp, &
      & 4.86112152716323E-04_wp,-6.63799363346240E-05_wp, 2.63148272705645E-04_wp, &
      &-3.69338959899261E-03_wp,-1.62131073853099E-02_wp, 1.12706765571131E-02_wp, &
      & 1.05297118737550E-02_wp,-6.83652635177475E-04_wp,-3.67543219032771E-06_wp, &
      &-5.64274993109804E-04_wp,-4.99490489384980E-05_wp,-1.70869611262570E-04_wp, &
      &-3.50216085251382E-03_wp, 4.54245265618630E-03_wp, 1.19055309005404E-02_wp, &
      &-1.89552860618799E-02_wp,-4.21614601950282E-04_wp,-8.16969303104830E-05_wp, &
      &-5.05453807147833E-04_wp, 3.80106625809516E-05_wp, 4.70419178273609E-04_wp, &
      &-3.09904488020287E-04_wp,-1.89524702709497E-04_wp, 2.43411943306908E-04_wp, &
      & 2.47928402744496E-03_wp,-1.13559737080272E-03_wp,-1.51157026309616E-03_wp, &
      &-1.74075334963782E-05_wp, 1.78516116698550E-03_wp, 2.08500104756471E-03_wp, &
      &-2.62166403396570E-03_wp, 1.57889816861724E-03_wp, 5.04088894827932E-03_wp, &
      &-2.33600703355515E-03_wp, 1.82408462819675E-03_wp, 3.82256307528407E-03_wp, &
      &-1.43702087113228E-03_wp, 4.61273866776237E-03_wp,-1.97237097978807E-02_wp, &
      &-1.20991357049060E-02_wp, 5.08336438918823E-03_wp,-2.83097626905179E-04_wp, &
      &-5.91992885185732E-04_wp, 1.38111314515451E-04_wp,-9.22681999267372E-05_wp, &
      & 7.37417254727631E-04_wp, 5.31578919345921E-03_wp, 2.10073897134008E-02_wp, &
      & 1.09429414011331E-02_wp, 3.20167044114226E-03_wp, 3.34046150738591E-04_wp, &
      &-7.19136030754649E-04_wp, 4.42987781476963E-05_wp, 2.30642331218241E-04_wp, &
      & 6.80365098406377E-04_wp,-8.84877256219645E-03_wp, 1.35788864585981E-02_wp, &
      &-2.03601998000497E-02_wp,-2.27788852327071E-02_wp,-6.78210598851328E-04_wp, &
      &-4.08503618156095E-04_wp,-1.59367877771540E-04_wp, 1.00614012061061E-03_wp, &
      & 6.42501172500220E-04_wp,-2.72087829479343E-03_wp,-2.70383567603876E-03_wp, &
      &-1.38558534276645E-03_wp,-1.13559737080272E-03_wp, 4.40223717129710E-02_wp, &
      & 5.57984043827394E-03_wp, 2.63086726948185E-03_wp, 7.70835864592983E-03_wp, &
      &-3.97661022360440E-04_wp, 2.26637457682227E-04_wp, 7.74500745271912E-03_wp, &
      &-1.79565826935409E-03_wp, 1.67084160854114E-02_wp,-6.61485609393775E-03_wp, &
      &-7.67995123123247E-04_wp, 2.18809239737792E-03_wp, 2.24397055345448E-02_wp, &
      & 6.64549541656991E-02_wp, 8.85871421288148E-02_wp, 6.96854518048087E-04_wp, &
      & 9.90846322158204E-04_wp, 2.46975374062160E-03_wp, 3.12047632632326E-04_wp, &
      & 1.29596421479028E-03_wp,-1.29311438090923E-03_wp,-2.40673203813875E-02_wp, &
      & 7.93390435947543E-02_wp, 7.02903965683193E-02_wp,-1.26672314141576E-02_wp, &
      & 1.72311948156631E-03_wp,-2.33142012668952E-03_wp, 1.41048483943771E-04_wp, &
      & 1.23754087717857E-03_wp, 1.32235237491106E-03_wp,-1.35903661755529E-02_wp, &
      & 7.65937840551499E-02_wp, 3.28686383455122E-02_wp, 6.44067474820589E-02_wp, &
      &-1.00233070467192E-03_wp,-6.88819688392179E-04_wp, 3.53857312562318E-04_wp, &
      &-1.40694498187750E-03_wp,-2.36540545013463E-03_wp, 6.74488694982622E-03_wp, &
      & 7.51687695834723E-04_wp,-1.29700970258115E-03_wp,-1.51157026309616E-03_wp, &
      & 5.57984043827394E-03_wp, 3.14464127057781E-02_wp, 5.42414731300720E-03_wp, &
      &-4.96716029907645E-03_wp,-7.32631149807317E-03_wp, 3.52486701106247E-03_wp, &
      &-1.72853100047814E-02_wp, 6.66597965316477E-03_wp, 1.41280628812100E-03_wp, &
      & 1.69880913558385E-03_wp, 6.38096027475187E-04_wp, 9.57106430361156E-03_wp, &
      & 2.07478739672553E-02_wp, 6.75325136802193E-02_wp,-2.15309610650728E-02_wp, &
      & 7.89233538673288E-02_wp, 2.81701658164297E-03_wp, 5.88757207563303E-04_wp, &
      &-1.82030938334733E-03_wp, 5.97326378135036E-04_wp,-5.00617661061709E-04_wp, &
      & 1.84804474365580E-02_wp,-6.44185602640370E-02_wp, 3.79320344196431E-02_wp, &
      & 6.87221148902145E-02_wp,-2.66981711043064E-03_wp, 3.40407076459113E-04_wp, &
      &-2.12977682745461E-03_wp,-3.25448101994960E-04_wp,-6.09023766967790E-04_wp, &
      &-5.22282280973725E-03_wp, 3.75844989308498E-02_wp,-1.69469206158735E-02_wp, &
      & 4.36814617081137E-02_wp,-7.35276572441748E-04_wp,-6.06690581183132E-04_wp, &
      & 6.06034207881716E-04_wp, 1.10660155867559E-04_wp,-1.13219077295460E-03_wp, &
      &-5.28266720158092E-03_wp, 1.57852557800278E-05_wp, 3.20274034784215E-04_wp, &
      &-1.74075334963782E-05_wp, 2.63086726948185E-03_wp, 5.42414731300720E-03_wp, &
      & 2.93589728175237E-02_wp, 3.67709738831901E-03_wp,-7.47063038600942E-04_wp, &
      & 3.90836704631231E-03_wp, 7.09621953957305E-05_wp, 2.07760552454923E-03_wp, &
      &-1.29405483042985E-04_wp, 2.89718143665865E-03_wp, 3.16800401100101E-04_wp, &
      &-1.69606053549435E-02_wp,-7.54471811196135E-03_wp, 1.55318190585745E-02_wp, &
      &-9.64076933204775E-02_wp, 1.10515854057058E-02_wp, 1.45947727308444E-04_wp, &
      &-2.16103631753385E-03_wp,-1.41298732458153E-03_wp,-1.48077115539182E-03_wp, &
      &-3.68872484716546E-05_wp,-9.73129495503267E-03_wp, 3.62507400512151E-04_wp, &
      & 8.64415228683387E-02_wp,-5.52582904666052E-03_wp, 2.36297575249053E-04_wp, &
      &-1.94458556002731E-03_wp,-1.28072235876436E-03_wp, 1.67826970775164E-03_wp, &
      &-2.64680965054583E-04_wp,-9.91183630949522E-03_wp,-5.57175718616637E-03_wp, &
      & 8.67153633545511E-02_wp, 2.90410287474541E-03_wp, 3.78053496299335E-04_wp, &
      & 9.38622543052269E-04_wp,-1.22901953593348E-03_wp,-2.36862461809407E-03_wp, &
      & 3.42392873702847E-05_wp, 4.80689218345537E-03_wp,-1.70421269428447E-03_wp, &
      &-6.00771636389316E-04_wp, 1.78516116698550E-03_wp, 7.70835864592983E-03_wp, &
      &-4.96716029907645E-03_wp, 3.67709738831901E-03_wp, 3.44193006387385E-02_wp, &
      & 4.76228788479958E-03_wp,-1.12594099044268E-02_wp,-4.22637413712781E-03_wp, &
      &-1.10988128216367E-02_wp, 9.39053014220174E-04_wp,-1.16060082138239E-03_wp, &
      & 1.07624737461691E-02_wp, 5.60381954555915E-03_wp, 1.41316668863593E-02_wp, &
      & 7.88193671551598E-02_wp,-1.33489350096580E-02_wp, 5.13853337547666E-03_wp, &
      & 1.55817983245919E-03_wp, 5.12251036351069E-04_wp,-1.38894336043584E-03_wp, &
      & 4.98434304366330E-06_wp,-1.55155659357601E-03_wp,-1.19435034444843E-02_wp, &
      & 7.87707633106839E-02_wp,-3.00272463633385E-02_wp,-2.09804044213811E-03_wp, &
      & 1.14598547458640E-03_wp,-1.69273880682273E-04_wp, 1.52074117567960E-03_wp, &
      &-6.69143343856180E-04_wp, 1.75734486260661E-03_wp, 2.14421083889658E-02_wp, &
      & 3.28723251997892E-02_wp, 4.78972976104863E-02_wp,-1.04712017981930E-01_wp, &
      &-1.59454356614833E-03_wp,-7.29860323332650E-05_wp,-2.39582887844565E-03_wp, &
      & 4.78040749239259E-04_wp, 2.59547000701250E-03_wp, 9.22601160315454E-04_wp, &
      &-1.24711040739595E-04_wp, 4.53126654291574E-04_wp, 2.08500104756471E-03_wp, &
      &-3.97661022360440E-04_wp,-7.32631149807317E-03_wp,-7.47063038600942E-04_wp, &
      & 4.76228788479958E-03_wp, 4.33364941711876E-02_wp,-3.15071146365869E-03_wp, &
      &-6.67108487341014E-04_wp, 6.99176018564855E-03_wp,-6.22633014911391E-03_wp, &
      &-3.38055775918460E-03_wp, 7.69017399475255E-03_wp,-1.76961128943864E-03_wp, &
      &-8.85374052511688E-03_wp, 3.47268154939585E-02_wp,-3.45265330160837E-02_wp, &
      &-8.92468896301809E-02_wp,-1.30909318967325E-03_wp,-3.70938201880121E-04_wp, &
      &-1.10750429751085E-04_wp,-1.31927484389834E-03_wp,-1.78842876359275E-03_wp, &
      &-8.63307150232595E-03_wp,-3.79182506231813E-02_wp, 3.52704972597514E-02_wp, &
      &-9.70603235678255E-02_wp, 1.09039808896337E-03_wp,-3.59947499874213E-04_wp, &
      &-3.60824594445763E-04_wp, 1.85106237169424E-03_wp,-2.29662810010527E-03_wp, &
      & 2.20191489255149E-02_wp, 7.33216018913339E-02_wp,-7.19861577005277E-02_wp, &
      &-4.69373915752390E-02_wp,-2.10917504939904E-03_wp,-1.45840439450576E-03_wp, &
      & 1.08123286916086E-04_wp, 2.31640389940747E-03_wp, 5.16614467010584E-04_wp, &
      &-2.71966921747329E-03_wp, 7.89473512446235E-04_wp, 6.89363513188489E-04_wp, &
      &-2.62166403396570E-03_wp, 2.26637457682227E-04_wp, 3.52486701106247E-03_wp, &
      & 3.90836704631231E-03_wp,-1.12594099044268E-02_wp,-3.15071146365869E-03_wp, &
      & 2.34151430732665E-02_wp, 1.85415808471368E-02_wp,-4.24734078208696E-02_wp, &
      & 2.35575691060635E-03_wp,-3.33058390377248E-02_wp,-5.57795006945784E-03_wp, &
      & 1.01108773595180E-02_wp, 4.01102854486234E-04_wp,-4.56200941513286E-02_wp, &
      &-6.60382781517692E-03_wp, 4.54256036669546E-02_wp, 8.13365307011402E-04_wp, &
      &-1.45138153891304E-03_wp,-4.61944048289146E-04_wp, 2.51862372221059E-04_wp, &
      & 2.07065310902136E-03_wp, 3.41550940616051E-04_wp,-5.80400261978084E-02_wp, &
      &-1.74917522456008E-02_wp,-5.30554321998480E-02_wp, 5.62332461332891E-04_wp, &
      & 1.78951892335046E-03_wp, 5.49306638133438E-04_wp, 1.18354420616184E-04_wp, &
      &-2.60744135288374E-03_wp,-1.84935397426444E-04_wp, 4.00823788730330E-02_wp, &
      & 5.07256361602297E-02_wp, 6.01728319806439E-02_wp,-1.08074030233043E-03_wp, &
      &-5.41519705143916E-04_wp, 1.74128670059436E-04_wp,-2.34898043484794E-03_wp, &
      &-2.72927683099058E-03_wp,-1.01894320857934E-02_wp,-3.57098061462513E-04_wp, &
      & 2.85458472253175E-03_wp, 1.57889816861724E-03_wp, 7.74500745271912E-03_wp, &
      &-1.72853100047814E-02_wp, 7.09621953957305E-05_wp,-4.22637413712781E-03_wp, &
      &-6.67108487341014E-04_wp, 1.85415808471368E-02_wp, 7.81937061967329E-02_wp, &
      &-1.56952868336072E-01_wp, 4.63674964372455E-03_wp,-1.52201455595406E-01_wp, &
      &-3.58517678579157E-02_wp, 3.07585591489104E-02_wp, 4.52883324935330E-04_wp, &
      &-9.32147819443025E-02_wp,-5.48349175520555E-03_wp,-3.42959429035113E-02_wp, &
      &-1.17137609071555E-03_wp,-4.30791961534006E-03_wp,-8.16503856499543E-04_wp, &
      &-1.99648246102047E-03_wp, 1.66375144799991E-03_wp, 7.55268863773708E-05_wp, &
      & 7.36998183180566E-02_wp, 8.59929205815937E-03_wp,-3.29569302508286E-02_wp, &
      & 2.35295984832753E-03_wp,-7.98432581448504E-04_wp, 2.29759184125395E-03_wp, &
      &-6.70858028737459E-04_wp, 2.19545673457071E-03_wp, 9.16685372636403E-04_wp, &
      & 2.57126867852129E-02_wp, 1.14966740661808E-02_wp, 3.88665260793635E-02_wp, &
      &-9.20983705261979E-04_wp,-1.88651521176194E-03_wp, 2.10285595586614E-03_wp, &
      & 5.97064769855728E-04_wp,-3.32760114530298E-03_wp,-2.27063129355292E-03_wp, &
      & 4.01348859449634E-03_wp,-8.42769400826849E-03_wp, 5.04088894827932E-03_wp, &
      &-1.79565826935409E-03_wp, 6.66597965316477E-03_wp, 2.07760552454923E-03_wp, &
      &-1.10988128216367E-02_wp, 6.99176018564855E-03_wp,-4.24734078208696E-02_wp, &
      &-1.56952868336072E-01_wp, 4.62813275840075E-01_wp,-9.10464660712884E-03_wp, &
      & 4.38047333436907E-01_wp, 9.59547188065274E-02_wp,-1.25382771228993E-01_wp, &
      &-1.38374238063024E-03_wp,-7.41273093046619E-03_wp,-1.10468068242443E-02_wp, &
      &-1.56187513033571E-03_wp,-4.41907081312552E-03_wp, 6.66432321271109E-03_wp, &
      & 6.93367015806001E-03_wp, 3.50242290460517E-03_wp, 1.53391697946383E-03_wp, &
      & 5.43354709829560E-04_wp, 1.40237003801708E-02_wp, 3.28804253498567E-02_wp, &
      & 1.21189122314942E-02_wp,-6.66975468342735E-05_wp,-3.13907873267953E-03_wp, &
      &-3.25713970999458E-03_wp, 4.31017066297599E-03_wp,-1.47887531116780E-03_wp, &
      &-1.84087838429997E-04_wp, 1.41827982064466E-02_wp,-6.62740370090869E-02_wp, &
      & 2.39879765891764E-02_wp, 1.05566543858823E-03_wp, 3.60304872458680E-03_wp, &
      &-3.54139161533717E-03_wp,-2.95929713103497E-03_wp, 4.28984898099735E-03_wp, &
      &-2.82181425186610E-03_wp,-3.38878489220265E-03_wp,-3.49285659189717E-03_wp, &
      &-2.33600703355515E-03_wp, 1.67084160854114E-02_wp, 1.41280628812100E-03_wp, &
      &-1.29405483042985E-04_wp, 9.39053014220174E-04_wp,-6.22633014911391E-03_wp, &
      & 2.35575691060635E-03_wp, 4.63674964372455E-03_wp,-9.10464660712884E-03_wp, &
      & 1.55693316213329E-02_wp,-8.69456072981448E-03_wp,-2.07001602146450E-03_wp, &
      &-1.64641127464159E-03_wp, 7.40412926607099E-04_wp, 5.19987722272070E-02_wp, &
      & 3.50963776224181E-02_wp, 3.53821139878523E-02_wp, 2.03766503438101E-03_wp, &
      & 1.38182723046831E-03_wp,-8.56882698093754E-04_wp, 9.67443679062229E-04_wp, &
      &-7.89383163152358E-04_wp,-3.12016989036287E-04_wp, 4.89673716424762E-02_wp, &
      & 6.85566681167747E-03_wp,-4.33726604714240E-02_wp, 2.32152383555889E-03_wp, &
      &-7.19828442745710E-04_wp, 1.28190858547543E-03_wp, 7.83301665646748E-04_wp, &
      & 1.95937371857659E-04_wp,-2.19901739927516E-04_wp,-2.48692793280075E-02_wp, &
      & 6.56157086454086E-03_wp, 6.03966283304102E-02_wp, 1.66757073509143E-03_wp, &
      & 4.66111281358563E-04_wp, 1.29891259804992E-03_wp,-8.75032383977239E-04_wp, &
      &-1.61471980448458E-03_wp,-2.14613394982406E-03_wp, 6.34879336465873E-03_wp, &
      &-8.73740458821759E-03_wp, 1.82408462819675E-03_wp,-6.61485609393775E-03_wp, &
      & 1.69880913558385E-03_wp, 2.89718143665865E-03_wp,-1.16060082138239E-03_wp, &
      &-3.38055775918460E-03_wp,-3.33058390377248E-02_wp,-1.52201455595406E-01_wp, &
      & 4.38047333436907E-01_wp,-8.69456072981448E-03_wp, 4.36466806654643E-01_wp, &
      & 9.34741416346391E-02_wp,-1.19354104145719E-01_wp,-1.37121127203996E-03_wp, &
      & 5.30583941925366E-04_wp,-2.67897186120732E-02_wp, 4.59737348249017E-03_wp, &
      &-4.02177217100425E-03_wp, 6.15256712974783E-03_wp, 6.28787735430965E-03_wp, &
      & 3.19802362532034E-03_wp, 1.32955044996654E-03_wp, 1.31200798108916E-03_wp, &
      & 1.24744952822560E-02_wp,-7.32813987751117E-02_wp, 3.72566964772602E-03_wp, &
      & 1.02640693375478E-04_wp,-1.47429100083848E-05_wp,-1.17086414454252E-03_wp, &
      & 2.24601344727871E-03_wp,-2.00559116482125E-03_wp,-1.12586025604532E-03_wp, &
      & 2.39815493206550E-02_wp, 2.95058137726838E-02_wp, 2.61728970650068E-02_wp, &
      & 5.32051158914855E-04_wp, 4.14351970344133E-03_wp,-5.00522584021680E-03_wp, &
      &-6.01050825635358E-03_wp, 3.53071480903333E-03_wp, 3.61217559598682E-03_wp, &
      &-2.30645681571837E-03_wp,-2.03180597844227E-03_wp, 3.82256307528407E-03_wp, &
      &-7.67995123123247E-04_wp, 6.38096027475187E-04_wp, 3.16800401100101E-04_wp, &
      & 1.07624737461691E-02_wp, 7.69017399475255E-03_wp,-5.57795006945784E-03_wp, &
      &-3.58517678579157E-02_wp, 9.59547188065274E-02_wp,-2.07001602146450E-03_wp, &
      & 9.34741416346391E-02_wp, 3.90520218749009E-02_wp,-2.06735101088362E-02_wp, &
      &-2.07679312253469E-04_wp,-3.36427527894663E-03_wp,-2.68963264938291E-03_wp, &
      & 7.67642844476174E-02_wp, 1.35665085459366E-03_wp, 1.35542875736036E-03_wp, &
      & 7.30947572329957E-04_wp, 2.17212419187750E-03_wp, 2.27422974116782E-03_wp, &
      &-5.83309848279943E-04_wp, 1.33899406086197E-02_wp,-2.03214366796728E-03_wp, &
      &-5.90050199663154E-02_wp, 1.83884922533798E-03_wp,-5.27133217581021E-04_wp, &
      & 3.10277119126789E-04_wp, 1.78806056491888E-03_wp,-1.42996134356551E-03_wp, &
      & 9.38729390030271E-05_wp, 4.82654826897998E-02_wp,-1.14992339499985E-02_wp, &
      &-6.93910045606182E-02_wp,-1.97145781885766E-03_wp, 3.56893171984317E-05_wp, &
      &-2.34151302881454E-03_wp, 4.72722357649428E-04_wp, 2.70282185999502E-03_wp, &
      & 1.12028617140403E-02_wp,-6.70882068372886E-04_wp, 3.57585309496864E-03_wp, &
      &-1.43702087113228E-03_wp, 2.18809239737792E-03_wp, 9.57106430361156E-03_wp, &
      &-1.69606053549435E-02_wp, 5.60381954555915E-03_wp,-1.76961128943864E-03_wp, &
      & 1.01108773595180E-02_wp, 3.07585591489104E-02_wp,-1.25382771228993E-01_wp, &
      &-1.64641127464159E-03_wp,-1.19354104145719E-01_wp,-2.06735101088362E-02_wp, &
      & 5.64165987850187E-02_wp,-3.60508414739047E-04_wp, 2.36891004261134E-02_wp, &
      & 7.45384160365845E-02_wp, 3.27358666273189E-02_wp, 2.50498630196250E-03_wp, &
      & 3.12765956430376E-04_wp,-1.49064981297606E-03_wp, 7.64346376852005E-04_wp, &
      &-3.81834224017286E-04_wp,-3.89914338131860E-04_wp,-5.16374969915572E-02_wp, &
      &-4.56666343915549E-02_wp, 4.41302193900407E-02_wp,-2.26312525917335E-03_wp, &
      & 2.64354771088890E-03_wp, 1.76499596075038E-04_wp,-2.64863626414929E-03_wp, &
      &-2.83279607122441E-05_wp,-2.41011357525776E-04_wp, 5.49357483690815E-02_wp, &
      &-2.23923554737875E-02_wp,-3.95054115591359E-02_wp,-2.59058745651069E-03_wp, &
      &-2.26084037291830E-03_wp, 4.94952180334675E-04_wp, 2.52712586788310E-03_wp, &
      &-5.02041178157069E-04_wp,-8.64033744988582E-04_wp, 6.96944723684011E-03_wp, &
      & 4.50241444890764E-03_wp, 4.61273866776237E-03_wp, 2.24397055345448E-02_wp, &
      & 2.07478739672553E-02_wp,-7.54471811196135E-03_wp, 1.41316668863593E-02_wp, &
      &-8.85374052511688E-03_wp, 4.01102854486234E-04_wp, 4.52883324935330E-04_wp, &
      &-1.38374238063024E-03_wp, 7.40412926607099E-04_wp,-1.37121127203996E-03_wp, &
      &-2.07679312253469E-04_wp,-3.60508414739047E-04_wp, 9.79385233275976E-01_wp, &
      &-6.73932117684685E-02_wp,-4.24221242518916E-02_wp,-4.54829773946277E-02_wp, &
      &-3.44767227988470E-03_wp,-3.20186908954796E-03_wp, 1.08518088940097E-03_wp, &
      &-2.17906881067377E-03_wp, 1.36803923011790E-03_wp,-1.74429112109240E-03_wp, &
      &-4.60039554543272E-03_wp,-4.62867164242500E-03_wp,-2.29461193380882E-02_wp, &
      & 1.08781582530049E-03_wp, 8.21284940310410E-04_wp, 8.63411726539932E-04_wp, &
      & 3.71234417790668E-04_wp,-2.26885600406139E-03_wp,-4.67052712135014E-04_wp, &
      &-2.18177810703358E-02_wp,-4.85160230143877E-03_wp, 2.44761467322725E-03_wp, &
      & 2.10650456051350E-03_wp, 5.13537734107091E-04_wp, 8.26572953275499E-04_wp, &
      & 7.14219643150120E-04_wp, 1.14929051465742E-03_wp, 3.10371374558505E-02_wp, &
      &-1.08436093353111E-02_wp,-1.94202483147287E-02_wp,-1.97237097978807E-02_wp, &
      & 6.64549541656991E-02_wp, 6.75325136802193E-02_wp, 1.55318190585745E-02_wp, &
      & 7.88193671551598E-02_wp, 3.47268154939585E-02_wp,-4.56200941513286E-02_wp, &
      &-9.32147819443025E-02_wp,-7.41273093046619E-03_wp, 5.19987722272070E-02_wp, &
      & 5.30583941925366E-04_wp,-3.36427527894663E-03_wp, 2.36891004261134E-02_wp, &
      &-6.73932117684685E-02_wp, 8.89614738721909E-01_wp,-1.28950331968545E-02_wp, &
      &-1.39586534563915E-02_wp, 1.53963342226206E-02_wp, 1.29909255553247E-02_wp, &
      &-1.31338339197474E-02_wp,-1.58767803767510E-05_wp,-2.41453153734081E-02_wp, &
      & 5.83630129305774E-03_wp, 4.74494717503171E-03_wp,-1.32063220713228E-02_wp, &
      &-1.19859393613625E-02_wp,-6.98434276117959E-04_wp, 2.29788467998815E-03_wp, &
      & 3.50868101514438E-04_wp,-9.19684784727941E-04_wp,-1.01534697329913E-03_wp, &
      &-1.79249305334630E-02_wp,-4.44215180079381E-02_wp,-4.45142330258693E-03_wp, &
      & 3.22808677054596E-02_wp, 5.51511024939793E-03_wp, 2.23194762526908E-03_wp, &
      & 2.14969208354250E-03_wp, 2.26334563100816E-04_wp,-2.69348432554597E-04_wp, &
      & 2.30333249958885E-02_wp,-1.77208852887962E-02_wp, 1.08132873675669E-02_wp, &
      &-1.20991357049060E-02_wp, 8.85871421288148E-02_wp,-2.15309610650728E-02_wp, &
      &-9.64076933204775E-02_wp,-1.33489350096580E-02_wp,-3.45265330160837E-02_wp, &
      &-6.60382781517692E-03_wp,-5.48349175520555E-03_wp,-1.10468068242443E-02_wp, &
      & 3.50963776224181E-02_wp,-2.67897186120732E-02_wp,-2.68963264938291E-03_wp, &
      & 7.45384160365845E-02_wp,-4.24221242518916E-02_wp,-1.28950331968545E-02_wp, &
      & 9.00965756211378E-01_wp,-9.14415453342455E-03_wp, 1.33841889823395E-03_wp, &
      & 2.15951288119434E-02_wp, 1.24101348943107E-02_wp, 1.49920354705176E-02_wp, &
      &-4.40469099734259E-04_wp, 9.85577803658296E-04_wp,-1.22201335360621E-02_wp, &
      & 2.11341370624279E-02_wp,-1.39791945191471E-02_wp,-7.01721881187569E-04_wp, &
      & 1.08687409942297E-04_wp,-1.44491706267255E-03_wp, 4.54075284812451E-04_wp, &
      &-1.98967421504811E-03_wp, 1.69947916252597E-03_wp,-1.66229817006137E-02_wp, &
      & 2.23035839447183E-02_wp,-9.27699182867206E-03_wp, 4.79918109486991E-04_wp, &
      & 1.40000483980811E-04_wp,-1.14187690596524E-03_wp, 2.76159075773128E-04_wp, &
      & 1.54213587967972E-03_wp, 2.04795146556320E-02_wp,-1.93565751770493E-02_wp, &
      &-1.31156713246053E-02_wp, 5.08336438918823E-03_wp, 6.96854518048087E-04_wp, &
      & 7.89233538673288E-02_wp, 1.10515854057058E-02_wp, 5.13853337547666E-03_wp, &
      &-8.92468896301809E-02_wp, 4.54256036669546E-02_wp,-3.42959429035113E-02_wp, &
      &-1.56187513033571E-03_wp, 3.53821139878523E-02_wp, 4.59737348249017E-03_wp, &
      & 7.67642844476174E-02_wp, 3.27358666273189E-02_wp,-4.54829773946277E-02_wp, &
      &-1.39586534563915E-02_wp,-9.14415453342455E-03_wp, 9.01197413335997E-01_wp, &
      & 2.37424843741826E-02_wp, 2.21729555528611E-04_wp,-8.75069972739857E-03_wp, &
      & 1.29168909226698E-02_wp, 1.70306110002549E-02_wp,-2.30415156877361E-02_wp, &
      & 1.66719242999167E-02_wp, 2.25021243598935E-03_wp,-5.05076972762118E-02_wp, &
      & 5.43204781194369E-03_wp,-9.57154239907275E-05_wp, 2.35591709642219E-03_wp, &
      & 2.63084615326781E-03_wp,-2.80062748714786E-03_wp, 1.34259756367277E-02_wp, &
      & 9.67852731639095E-03_wp,-1.30717775429154E-02_wp,-3.47522554508627E-03_wp, &
      &-2.41917211488525E-03_wp,-1.02717741612734E-03_wp,-7.37519467381799E-04_wp, &
      & 2.37177085329030E-03_wp, 1.59044196104693E-03_wp, 1.17566904494576E-03_wp, &
      &-8.75401591255748E-04_wp,-6.35109734800467E-04_wp,-2.83097626905179E-04_wp, &
      & 9.90846322158204E-04_wp, 2.81701658164297E-03_wp, 1.45947727308444E-04_wp, &
      & 1.55817983245919E-03_wp,-1.30909318967325E-03_wp, 8.13365307011402E-04_wp, &
      &-1.17137609071555E-03_wp,-4.41907081312552E-03_wp, 2.03766503438101E-03_wp, &
      &-4.02177217100425E-03_wp, 1.35665085459366E-03_wp, 2.50498630196250E-03_wp, &
      &-3.44767227988470E-03_wp, 1.53963342226206E-02_wp, 1.33841889823395E-03_wp, &
      & 2.37424843741826E-02_wp, 9.71080523129421E-04_wp, 2.24537932149413E-04_wp, &
      &-5.12935597024076E-04_wp, 3.54826393035252E-04_wp, 1.78175387129317E-05_wp, &
      &-9.92864289498025E-04_wp,-4.93805000442882E-04_wp,-1.05118342010505E-03_wp, &
      &-5.15532684545041E-03_wp, 2.01568567477248E-04_wp, 8.47489533147879E-05_wp, &
      & 1.28915036185613E-04_wp, 5.96716941131621E-05_wp,-1.61436989305845E-04_wp, &
      & 9.05992698758370E-04_wp,-1.86454186027926E-03_wp,-1.56221941949092E-03_wp, &
      &-8.01779151306951E-04_wp, 5.29754230191627E-05_wp,-1.67720095156313E-05_wp, &
      & 6.88550780578720E-05_wp, 1.49459436391235E-04_wp, 3.48745958377767E-05_wp, &
      & 1.12383332750122E-03_wp,-5.65296230455721E-04_wp,-1.41358099308255E-04_wp, &
      &-5.91992885185732E-04_wp, 2.46975374062160E-03_wp, 5.88757207563303E-04_wp, &
      &-2.16103631753385E-03_wp, 5.12251036351069E-04_wp,-3.70938201880121E-04_wp, &
      &-1.45138153891304E-03_wp,-4.30791961534006E-03_wp, 6.66432321271109E-03_wp, &
      & 1.38182723046831E-03_wp, 6.15256712974783E-03_wp, 1.35542875736036E-03_wp, &
      & 3.12765956430376E-04_wp,-3.20186908954796E-03_wp, 1.29909255553247E-02_wp, &
      & 2.15951288119434E-02_wp, 2.21729555528611E-04_wp, 2.24537932149413E-04_wp, &
      & 8.41821697249856E-04_wp, 2.01357871247298E-04_wp, 4.33545604699002E-04_wp, &
      &-3.54326607641147E-04_wp, 6.92419916188416E-04_wp,-2.26880794835676E-03_wp, &
      &-2.60472690045053E-04_wp,-6.40001479018286E-04_wp,-5.97084418734319E-05_wp, &
      & 4.30725140318215E-05_wp,-9.13448849020790E-05_wp, 5.35251386103987E-05_wp, &
      &-1.48760158877089E-04_wp,-3.65628894303910E-04_wp,-3.34642237401371E-03_wp, &
      &-3.20723614997892E-04_wp,-3.02716907724124E-04_wp, 1.66642160459605E-04_wp, &
      & 1.29518345211162E-04_wp,-6.06645445272479E-05_wp,-3.58157141419361E-05_wp, &
      & 1.47884993198466E-04_wp,-3.37161261858860E-04_wp, 1.94865564588492E-04_wp, &
      & 4.86112152716323E-04_wp, 1.38111314515451E-04_wp, 3.12047632632326E-04_wp, &
      &-1.82030938334733E-03_wp,-1.41298732458153E-03_wp,-1.38894336043584E-03_wp, &
      &-1.10750429751085E-04_wp,-4.61944048289146E-04_wp,-8.16503856499543E-04_wp, &
      & 6.93367015806001E-03_wp,-8.56882698093754E-04_wp, 6.28787735430965E-03_wp, &
      & 7.30947572329957E-04_wp,-1.49064981297606E-03_wp, 1.08518088940097E-03_wp, &
      &-1.31338339197474E-02_wp, 1.24101348943107E-02_wp,-8.75069972739857E-03_wp, &
      &-5.12935597024076E-04_wp, 2.01357871247298E-04_wp, 5.57704342370171E-04_wp, &
      & 1.30088685793062E-04_wp, 2.05807133619415E-04_wp, 7.04869867156499E-04_wp, &
      &-4.71965410425160E-04_wp, 1.54322398343414E-03_wp, 1.57912740212179E-03_wp, &
      &-8.60224070514466E-05_wp,-8.26444834400528E-05_wp,-1.22481720955816E-04_wp, &
      & 5.33348705143874E-05_wp,-3.16133131285991E-06_wp, 6.83482385414134E-04_wp, &
      & 1.43421135623885E-03_wp, 1.29954197473103E-03_wp,-6.84034941761651E-04_wp, &
      &-6.77514813138002E-05_wp, 3.57114046695395E-05_wp,-1.37991308366822E-04_wp, &
      &-1.10652382975503E-04_wp, 8.18569440148499E-05_wp, 7.55695280205111E-04_wp, &
      &-6.09542770974166E-04_wp,-6.63799363346240E-05_wp,-9.22681999267372E-05_wp, &
      & 1.29596421479028E-03_wp, 5.97326378135036E-04_wp,-1.48077115539182E-03_wp, &
      & 4.98434304366330E-06_wp,-1.31927484389834E-03_wp, 2.51862372221059E-04_wp, &
      &-1.99648246102047E-03_wp, 3.50242290460517E-03_wp, 9.67443679062229E-04_wp, &
      & 3.19802362532034E-03_wp, 2.17212419187750E-03_wp, 7.64346376852005E-04_wp, &
      &-2.17906881067377E-03_wp,-1.58767803767510E-05_wp, 1.49920354705176E-02_wp, &
      & 1.29168909226698E-02_wp, 3.54826393035252E-04_wp, 4.33545604699002E-04_wp, &
      & 1.30088685793062E-04_wp, 4.82820109679038E-04_wp, 2.51108863693032E-04_wp, &
      &-7.79172865408991E-04_wp,-7.66705658343985E-04_wp, 2.26017817015977E-04_wp, &
      &-3.30674201547611E-03_wp, 1.03484025350007E-04_wp,-1.06917681599775E-06_wp, &
      &-1.97825617009888E-07_wp, 1.08987968326159E-04_wp,-1.54513083585034E-04_wp, &
      & 9.18499763897787E-04_wp, 2.49064579616363E-04_wp,-3.01515131459455E-04_wp, &
      &-1.87144089307596E-03_wp,-4.15809966344438E-05_wp, 1.50166180255650E-05_wp, &
      &-9.17945277861075E-05_wp, 3.21446802815042E-05_wp, 1.29877997002024E-04_wp, &
      &-4.91653601591487E-04_wp,-6.86637147630286E-05_wp, 2.63148272705645E-04_wp, &
      & 7.37417254727631E-04_wp,-1.29311438090923E-03_wp,-5.00617661061709E-04_wp, &
      &-3.68872484716546E-05_wp,-1.55155659357601E-03_wp,-1.78842876359275E-03_wp, &
      & 2.07065310902136E-03_wp, 1.66375144799991E-03_wp, 1.53391697946383E-03_wp, &
      &-7.89383163152358E-04_wp, 1.32955044996654E-03_wp, 2.27422974116782E-03_wp, &
      &-3.81834224017286E-04_wp, 1.36803923011790E-03_wp,-2.41453153734081E-02_wp, &
      &-4.40469099734259E-04_wp, 1.70306110002549E-02_wp, 1.78175387129317E-05_wp, &
      &-3.54326607641147E-04_wp, 2.05807133619415E-04_wp, 2.51108863693032E-04_wp, &
      & 1.00296300268016E-03_wp,-2.26339054665822E-03_wp, 1.22673110868481E-03_wp, &
      & 1.42194057769481E-03_wp,-2.96462717477230E-03_wp, 1.86342724132444E-04_wp, &
      &-1.01495495974295E-04_wp, 4.60297941198317E-05_wp, 1.31730529442160E-04_wp, &
      &-5.02008128431916E-05_wp, 2.19806549831237E-03_wp, 5.36800952916668E-03_wp, &
      &-1.11865135945953E-04_wp,-2.11155401443314E-03_wp,-3.09423947555407E-04_wp, &
      &-1.21929197698885E-04_wp,-1.22852682598501E-04_wp, 4.05393279044262E-05_wp, &
      & 4.47136744943953E-05_wp,-1.02695166967632E-03_wp,-7.77156635669540E-03_wp, &
      &-3.69338959899261E-03_wp, 5.31578919345921E-03_wp,-2.40673203813875E-02_wp, &
      & 1.84804474365580E-02_wp,-9.73129495503267E-03_wp,-1.19435034444843E-02_wp, &
      &-8.63307150232595E-03_wp, 3.41550940616051E-04_wp, 7.55268863773708E-05_wp, &
      & 5.43354709829560E-04_wp,-3.12016989036287E-04_wp, 1.31200798108916E-03_wp, &
      &-5.83309848279943E-04_wp,-3.89914338131860E-04_wp,-1.74429112109240E-03_wp, &
      & 5.83630129305774E-03_wp, 9.85577803658296E-04_wp,-2.30415156877361E-02_wp, &
      &-9.92864289498025E-04_wp, 6.92419916188416E-04_wp, 7.04869867156499E-04_wp, &
      &-7.79172865408991E-04_wp,-2.26339054665822E-03_wp, 9.80319653876714E-01_wp, &
      & 6.96148081109739E-02_wp, 3.39853186201137E-02_wp,-4.74581076068337E-02_wp, &
      & 3.77653932652171E-03_wp,-2.82491840284833E-03_wp, 1.53601953767167E-03_wp, &
      & 1.88938374714941E-03_wp, 1.39873724097286E-03_wp,-3.03607735046960E-03_wp, &
      & 1.59310863120097E-02_wp, 1.87885921811317E-02_wp, 5.08467535689675E-03_wp, &
      &-1.69023816854923E-03_wp,-6.52194452940651E-04_wp,-1.24575925065985E-03_wp, &
      &-1.84891013913080E-03_wp,-3.95093353419281E-04_wp,-3.45289149212074E-02_wp, &
      &-1.12023415897657E-02_wp,-1.62131073853099E-02_wp, 2.10073897134008E-02_wp, &
      & 7.93390435947543E-02_wp,-6.44185602640370E-02_wp, 3.62507400512151E-04_wp, &
      & 7.87707633106839E-02_wp,-3.79182506231813E-02_wp,-5.80400261978084E-02_wp, &
      & 7.36998183180566E-02_wp, 1.40237003801708E-02_wp, 4.89673716424762E-02_wp, &
      & 1.24744952822560E-02_wp, 1.33899406086197E-02_wp,-5.16374969915572E-02_wp, &
      &-4.60039554543272E-03_wp, 4.74494717503171E-03_wp,-1.22201335360621E-02_wp, &
      & 1.66719242999167E-02_wp,-4.93805000442882E-04_wp,-2.26880794835676E-03_wp, &
      &-4.71965410425160E-04_wp,-7.66705658343985E-04_wp, 1.22673110868481E-03_wp, &
      & 6.96148081109739E-02_wp, 8.84483522477083E-01_wp,-6.27087331785375E-03_wp, &
      & 1.39890422776069E-02_wp, 1.70641049794487E-02_wp,-1.41030186631501E-02_wp, &
      & 1.27499204873507E-02_wp, 3.86236907823217E-04_wp, 2.29766788920185E-02_wp, &
      & 9.93178929004829E-03_wp,-2.68014119423897E-03_wp,-3.97830522407406E-02_wp, &
      &-1.26975405205949E-02_wp, 1.43642412751397E-03_wp,-5.15623780968630E-04_wp, &
      & 1.97029419195767E-03_wp, 3.39979714316987E-03_wp, 1.57372281668683E-04_wp, &
      &-1.34121866165070E-02_wp,-1.78971328951955E-02_wp, 1.12706765571131E-02_wp, &
      & 1.09429414011331E-02_wp, 7.02903965683193E-02_wp, 3.79320344196431E-02_wp, &
      & 8.64415228683387E-02_wp,-3.00272463633385E-02_wp, 3.52704972597514E-02_wp, &
      &-1.74917522456008E-02_wp, 8.59929205815937E-03_wp, 3.28804253498567E-02_wp, &
      & 6.85566681167747E-03_wp,-7.32813987751117E-02_wp,-2.03214366796728E-03_wp, &
      &-4.56666343915549E-02_wp,-4.62867164242500E-03_wp,-1.32063220713228E-02_wp, &
      & 2.11341370624279E-02_wp, 2.25021243598935E-03_wp,-1.05118342010505E-03_wp, &
      &-2.60472690045053E-04_wp, 1.54322398343414E-03_wp, 2.26017817015977E-04_wp, &
      & 1.42194057769481E-03_wp, 3.39853186201137E-02_wp,-6.27087331785375E-03_wp, &
      & 9.01725438066334E-01_wp, 8.33818371870478E-03_wp,-3.41304494414273E-04_wp, &
      &-2.29885786612176E-02_wp,-1.59420543751911E-02_wp, 1.67323852681780E-02_wp, &
      & 5.90543543263457E-04_wp, 1.89070265565559E-02_wp,-2.18204966861355E-02_wp, &
      &-1.64116449632517E-02_wp,-3.50049447182234E-02_wp, 1.92048599024827E-03_wp, &
      & 2.73294147837168E-04_wp, 2.40805559719466E-04_wp, 3.98704788703455E-03_wp, &
      & 3.83007723351944E-03_wp, 2.17842958084033E-02_wp, 2.23980508202680E-02_wp, &
      & 1.05297118737550E-02_wp, 3.20167044114226E-03_wp,-1.26672314141576E-02_wp, &
      & 6.87221148902145E-02_wp,-5.52582904666052E-03_wp,-2.09804044213811E-03_wp, &
      &-9.70603235678255E-02_wp,-5.30554321998480E-02_wp,-3.29569302508286E-02_wp, &
      & 1.21189122314942E-02_wp,-4.33726604714240E-02_wp, 3.72566964772602E-03_wp, &
      &-5.90050199663154E-02_wp, 4.41302193900407E-02_wp,-2.29461193380882E-02_wp, &
      &-1.19859393613625E-02_wp,-1.39791945191471E-02_wp,-5.05076972762118E-02_wp, &
      &-5.15532684545041E-03_wp,-6.40001479018286E-04_wp, 1.57912740212179E-03_wp, &
      &-3.30674201547611E-03_wp,-2.96462717477230E-03_wp,-4.74581076068337E-02_wp, &
      & 1.39890422776069E-02_wp, 8.33818371870478E-03_wp, 9.00709592745117E-01_wp, &
      &-2.32277295829378E-02_wp,-5.91310561324898E-04_wp,-1.03765797970804E-02_wp, &
      &-1.37570262673697E-02_wp, 1.79170712797512E-02_wp, 1.33699745251597E-02_wp, &
      & 1.12806043159290E-02_wp,-9.44493581999449E-03_wp,-1.90246054656196E-03_wp, &
      &-1.78419681380516E-03_wp,-9.11286541924727E-04_wp,-5.17440346509960E-04_wp, &
      & 2.69307620534643E-03_wp, 1.58235914556952E-03_wp,-1.32239449246865E-03_wp, &
      &-8.58573131873999E-04_wp,-6.83652635177475E-04_wp, 3.34046150738591E-04_wp, &
      & 1.72311948156631E-03_wp,-2.66981711043064E-03_wp, 2.36297575249053E-04_wp, &
      & 1.14598547458640E-03_wp, 1.09039808896337E-03_wp, 5.62332461332891E-04_wp, &
      & 2.35295984832753E-03_wp,-6.66975468342735E-05_wp, 2.32152383555889E-03_wp, &
      & 1.02640693375478E-04_wp, 1.83884922533798E-03_wp,-2.26312525917335E-03_wp, &
      & 1.08781582530049E-03_wp,-6.98434276117959E-04_wp,-7.01721881187569E-04_wp, &
      & 5.43204781194369E-03_wp, 2.01568567477248E-04_wp,-5.97084418734319E-05_wp, &
      &-8.60224070514466E-05_wp, 1.03484025350007E-04_wp, 1.86342724132444E-04_wp, &
      & 3.77653932652171E-03_wp, 1.70641049794487E-02_wp,-3.41304494414273E-04_wp, &
      &-2.32277295829378E-02_wp, 9.71954797250523E-04_wp,-2.66463312336675E-04_wp, &
      & 5.36546546854730E-04_wp, 3.78261896706480E-04_wp,-2.20377842347737E-05_wp, &
      &-1.08352826422085E-03_wp,-1.82074210885778E-03_wp,-1.01687330466814E-03_wp, &
      & 1.84168810430206E-03_wp, 1.24242288971551E-04_wp, 2.60909585424519E-05_wp, &
      & 9.82423307899891E-05_wp,-1.33251711439347E-05_wp,-7.80121783922660E-05_wp, &
      & 9.23763188396611E-04_wp, 6.50886246445002E-04_wp,-3.67543219032771E-06_wp, &
      &-7.19136030754649E-04_wp,-2.33142012668952E-03_wp, 3.40407076459113E-04_wp, &
      &-1.94458556002731E-03_wp,-1.69273880682273E-04_wp,-3.59947499874213E-04_wp, &
      & 1.78951892335046E-03_wp,-7.98432581448504E-04_wp,-3.13907873267953E-03_wp, &
      &-7.19828442745710E-04_wp,-1.47429100083848E-05_wp,-5.27133217581021E-04_wp, &
      & 2.64354771088890E-03_wp, 8.21284940310410E-04_wp, 2.29788467998815E-03_wp, &
      & 1.08687409942297E-04_wp,-9.57154239907275E-05_wp, 8.47489533147879E-05_wp, &
      & 4.30725140318215E-05_wp,-8.26444834400528E-05_wp,-1.06917681599775E-06_wp, &
      &-1.01495495974295E-04_wp,-2.82491840284833E-03_wp,-1.41030186631501E-02_wp, &
      &-2.29885786612176E-02_wp,-5.91310561324898E-04_wp,-2.66463312336675E-04_wp, &
      & 8.50416540479365E-04_wp, 2.20268366240989E-04_wp,-4.57895934641153E-04_wp, &
      &-3.97935121573486E-04_wp,-1.96526810804084E-03_wp, 3.13903002693610E-03_wp, &
      & 3.85187428687074E-03_wp, 2.11734626341689E-03_wp,-1.34264453636114E-04_wp, &
      &-2.51838077169850E-05_wp,-5.85408576231827E-05_wp,-2.31098244614289E-04_wp, &
      &-1.74472472563144E-04_wp,-5.87060101043125E-04_wp,-1.23589021151666E-04_wp, &
      &-5.64274993109804E-04_wp, 4.42987781476963E-05_wp, 1.41048483943771E-04_wp, &
      &-2.12977682745461E-03_wp,-1.28072235876436E-03_wp, 1.52074117567960E-03_wp, &
      &-3.60824594445763E-04_wp, 5.49306638133438E-04_wp, 2.29759184125395E-03_wp, &
      &-3.25713970999458E-03_wp, 1.28190858547543E-03_wp,-1.17086414454252E-03_wp, &
      & 3.10277119126789E-04_wp, 1.76499596075038E-04_wp, 8.63411726539932E-04_wp, &
      & 3.50868101514438E-04_wp,-1.44491706267255E-03_wp, 2.35591709642219E-03_wp, &
      & 1.28915036185613E-04_wp,-9.13448849020790E-05_wp,-1.22481720955816E-04_wp, &
      &-1.97825617009888E-07_wp, 4.60297941198317E-05_wp, 1.53601953767167E-03_wp, &
      & 1.27499204873507E-02_wp,-1.59420543751911E-02_wp,-1.03765797970804E-02_wp, &
      & 5.36546546854730E-04_wp, 2.20268366240989E-04_wp, 6.09265712052873E-04_wp, &
      &-1.47519943592550E-04_wp, 1.22238728620312E-04_wp,-1.25095952083955E-03_wp, &
      & 2.27846314199926E-04_wp, 1.88274865285337E-04_wp, 2.05809775538211E-03_wp, &
      & 1.14606708088033E-05_wp,-2.53947944807633E-05_wp, 7.94171880800760E-05_wp, &
      &-5.67477088132154E-05_wp,-1.57956030126034E-04_wp,-6.06266471529292E-04_wp, &
      &-7.14602097206055E-04_wp,-4.99490489384980E-05_wp, 2.30642331218241E-04_wp, &
      & 1.23754087717857E-03_wp,-3.25448101994960E-04_wp, 1.67826970775164E-03_wp, &
      &-6.69143343856180E-04_wp, 1.85106237169424E-03_wp, 1.18354420616184E-04_wp, &
      &-6.70858028737459E-04_wp, 4.31017066297599E-03_wp, 7.83301665646748E-04_wp, &
      & 2.24601344727871E-03_wp, 1.78806056491888E-03_wp,-2.64863626414929E-03_wp, &
      & 3.71234417790668E-04_wp,-9.19684784727941E-04_wp, 4.54075284812451E-04_wp, &
      & 2.63084615326781E-03_wp, 5.96716941131621E-05_wp, 5.35251386103987E-05_wp, &
      & 5.33348705143874E-05_wp, 1.08987968326159E-04_wp, 1.31730529442160E-04_wp, &
      & 1.88938374714941E-03_wp, 3.86236907823217E-04_wp, 1.67323852681780E-02_wp, &
      &-1.37570262673697E-02_wp, 3.78261896706480E-04_wp,-4.57895934641153E-04_wp, &
      &-1.47519943592550E-04_wp, 5.67070951884331E-04_wp,-2.66365841732221E-04_wp, &
      &-2.17685295256138E-05_wp,-2.21004109845814E-03_wp,-1.02849359838044E-03_wp, &
      &-9.02442359029762E-04_wp, 1.12455065949355E-04_wp, 6.84526285750054E-05_wp, &
      &-1.17897739759140E-05_wp, 2.91415636443471E-05_wp, 1.16342256034267E-04_wp, &
      &-5.25957119782735E-04_wp, 1.61848868703983E-04_wp,-1.70869611262570E-04_wp, &
      & 6.80365098406377E-04_wp, 1.32235237491106E-03_wp,-6.09023766967790E-04_wp, &
      &-2.64680965054583E-04_wp, 1.75734486260661E-03_wp,-2.29662810010527E-03_wp, &
      &-2.60744135288374E-03_wp, 2.19545673457071E-03_wp,-1.47887531116780E-03_wp, &
      & 1.95937371857659E-04_wp,-2.00559116482125E-03_wp,-1.42996134356551E-03_wp, &
      &-2.83279607122441E-05_wp,-2.26885600406139E-03_wp,-1.01534697329913E-03_wp, &
      &-1.98967421504811E-03_wp,-2.80062748714786E-03_wp,-1.61436989305845E-04_wp, &
      &-1.48760158877089E-04_wp,-3.16133131285991E-06_wp,-1.54513083585034E-04_wp, &
      &-5.02008128431916E-05_wp, 1.39873724097286E-03_wp, 2.29766788920185E-02_wp, &
      & 5.90543543263457E-04_wp, 1.79170712797512E-02_wp,-2.20377842347737E-05_wp, &
      &-3.97935121573486E-04_wp, 1.22238728620312E-04_wp,-2.66365841732221E-04_wp, &
      & 9.76918208351756E-04_wp, 1.37310447814570E-03_wp,-3.85050874145553E-04_wp, &
      &-4.29155366733989E-03_wp,-8.98951630102113E-04_wp, 7.11466472221123E-06_wp, &
      &-7.22055400537447E-05_wp, 1.11585039120817E-04_wp, 2.53643008087482E-04_wp, &
      & 2.63123070323647E-05_wp,-9.44887800911054E-04_wp, 2.47681579670218E-03_wp, &
      &-3.50216085251382E-03_wp,-8.84877256219645E-03_wp,-1.35903661755529E-02_wp, &
      &-5.22282280973725E-03_wp,-9.91183630949522E-03_wp, 2.14421083889658E-02_wp, &
      & 2.20191489255149E-02_wp,-1.84935397426444E-04_wp, 9.16685372636403E-04_wp, &
      &-1.84087838429997E-04_wp,-2.19901739927516E-04_wp,-1.12586025604532E-03_wp, &
      & 9.38729390030271E-05_wp,-2.41011357525776E-04_wp,-4.67052712135014E-04_wp, &
      &-1.79249305334630E-02_wp, 1.69947916252597E-03_wp, 1.34259756367277E-02_wp, &
      & 9.05992698758370E-04_wp,-3.65628894303910E-04_wp, 6.83482385414134E-04_wp, &
      & 9.18499763897787E-04_wp, 2.19806549831237E-03_wp,-3.03607735046960E-03_wp, &
      & 9.93178929004829E-03_wp, 1.89070265565559E-02_wp, 1.33699745251597E-02_wp, &
      &-1.08352826422085E-03_wp,-1.96526810804084E-03_wp,-1.25095952083955E-03_wp, &
      &-2.17685295256138E-05_wp, 1.37310447814570E-03_wp, 9.79754791069693E-01_wp, &
      &-2.23989504876099E-02_wp, 3.40636152248508E-02_wp, 8.18979621687480E-02_wp, &
      & 2.09012487006980E-03_wp, 8.40503347739840E-04_wp, 1.56038857676967E-03_wp, &
      &-3.26760910351179E-03_wp,-3.45211447044840E-03_wp, 9.71328579482040E-03_wp, &
      & 1.51505155319696E-02_wp, 4.54245265618630E-03_wp, 1.35788864585981E-02_wp, &
      & 7.65937840551499E-02_wp, 3.75844989308498E-02_wp,-5.57175718616637E-03_wp, &
      & 3.28723251997892E-02_wp, 7.33216018913339E-02_wp, 4.00823788730330E-02_wp, &
      & 2.57126867852129E-02_wp, 1.41827982064466E-02_wp,-2.48692793280075E-02_wp, &
      & 2.39815493206550E-02_wp, 4.82654826897998E-02_wp, 5.49357483690815E-02_wp, &
      &-2.18177810703358E-02_wp,-4.44215180079381E-02_wp,-1.66229817006137E-02_wp, &
      & 9.67852731639095E-03_wp,-1.86454186027926E-03_wp,-3.34642237401371E-03_wp, &
      & 1.43421135623885E-03_wp, 2.49064579616363E-04_wp, 5.36800952916668E-03_wp, &
      & 1.59310863120097E-02_wp,-2.68014119423897E-03_wp,-2.18204966861355E-02_wp, &
      & 1.12806043159290E-02_wp,-1.82074210885778E-03_wp, 3.13903002693610E-03_wp, &
      & 2.27846314199926E-04_wp,-2.21004109845814E-03_wp,-3.85050874145553E-04_wp, &
      &-2.23989504876099E-02_wp, 9.07216790441305E-01_wp, 4.71099439338643E-03_wp, &
      & 6.98669050585301E-03_wp,-2.76452181520743E-02_wp,-1.32617256121615E-02_wp, &
      &-5.84858535888869E-03_wp,-9.24621410537155E-04_wp,-8.64305835007620E-03_wp, &
      &-1.32264514235801E-02_wp, 4.13303227393656E-03_wp, 1.19055309005404E-02_wp, &
      &-2.03601998000497E-02_wp, 3.28686383455122E-02_wp,-1.69469206158735E-02_wp, &
      & 8.67153633545511E-02_wp, 4.78972976104863E-02_wp,-7.19861577005277E-02_wp, &
      & 5.07256361602297E-02_wp, 1.14966740661808E-02_wp,-6.62740370090869E-02_wp, &
      & 6.56157086454086E-03_wp, 2.95058137726838E-02_wp,-1.14992339499985E-02_wp, &
      &-2.23923554737875E-02_wp,-4.85160230143877E-03_wp,-4.45142330258693E-03_wp, &
      & 2.23035839447183E-02_wp,-1.30717775429154E-02_wp,-1.56221941949092E-03_wp, &
      &-3.20723614997892E-04_wp, 1.29954197473103E-03_wp,-3.01515131459455E-04_wp, &
      &-1.11865135945953E-04_wp, 1.87885921811317E-02_wp,-3.97830522407406E-02_wp, &
      &-1.64116449632517E-02_wp,-9.44493581999449E-03_wp,-1.01687330466814E-03_wp, &
      & 3.85187428687074E-03_wp, 1.88274865285337E-04_wp,-1.02849359838044E-03_wp, &
      &-4.29155366733989E-03_wp, 3.40636152248508E-02_wp, 4.71099439338643E-03_wp, &
      & 9.03653999747146E-01_wp,-1.07477493790883E-02_wp,-6.79035487613655E-04_wp, &
      & 7.65027028289048E-03_wp,-1.49468871725830E-02_wp,-2.63446016222323E-02_wp, &
      &-1.04417300489744E-03_wp,-4.00506954181809E-02_wp, 1.19937863231962E-02_wp, &
      &-1.89552860618799E-02_wp,-2.27788852327071E-02_wp, 6.44067474820589E-02_wp, &
      & 4.36814617081137E-02_wp, 2.90410287474541E-03_wp,-1.04712017981930E-01_wp, &
      &-4.69373915752390E-02_wp, 6.01728319806439E-02_wp, 3.88665260793635E-02_wp, &
      & 2.39879765891764E-02_wp, 6.03966283304102E-02_wp, 2.61728970650068E-02_wp, &
      &-6.93910045606182E-02_wp,-3.95054115591359E-02_wp, 2.44761467322725E-03_wp, &
      & 3.22808677054596E-02_wp,-9.27699182867206E-03_wp,-3.47522554508627E-03_wp, &
      &-8.01779151306951E-04_wp,-3.02716907724124E-04_wp,-6.84034941761651E-04_wp, &
      &-1.87144089307596E-03_wp,-2.11155401443314E-03_wp, 5.08467535689675E-03_wp, &
      &-1.26975405205949E-02_wp,-3.50049447182234E-02_wp,-1.90246054656196E-03_wp, &
      & 1.84168810430206E-03_wp, 2.11734626341689E-03_wp, 2.05809775538211E-03_wp, &
      &-9.02442359029762E-04_wp,-8.98951630102113E-04_wp, 8.18979621687480E-02_wp, &
      & 6.98669050585301E-03_wp,-1.07477493790883E-02_wp, 8.78220908101580E-01_wp, &
      & 8.78198962890330E-03_wp, 7.29973779501288E-04_wp, 1.48135170964987E-02_wp, &
      &-1.41265088250698E-02_wp,-2.63165868327967E-02_wp,-7.30975427202173E-04_wp, &
      &-4.00974483671832E-04_wp,-4.21614601950282E-04_wp,-6.78210598851328E-04_wp, &
      &-1.00233070467192E-03_wp,-7.35276572441748E-04_wp, 3.78053496299335E-04_wp, &
      &-1.59454356614833E-03_wp,-2.10917504939904E-03_wp,-1.08074030233043E-03_wp, &
      &-9.20983705261979E-04_wp, 1.05566543858823E-03_wp, 1.66757073509143E-03_wp, &
      & 5.32051158914855E-04_wp,-1.97145781885766E-03_wp,-2.59058745651069E-03_wp, &
      & 2.10650456051350E-03_wp, 5.51511024939793E-03_wp, 4.79918109486991E-04_wp, &
      &-2.41917211488525E-03_wp, 5.29754230191627E-05_wp, 1.66642160459605E-04_wp, &
      &-6.77514813138002E-05_wp,-4.15809966344438E-05_wp,-3.09423947555407E-04_wp, &
      &-1.69023816854923E-03_wp, 1.43642412751397E-03_wp, 1.92048599024827E-03_wp, &
      &-1.78419681380516E-03_wp, 1.24242288971551E-04_wp,-1.34264453636114E-04_wp, &
      & 1.14606708088033E-05_wp, 1.12455065949355E-04_wp, 7.11466472221123E-06_wp, &
      & 2.09012487006980E-03_wp,-2.76452181520743E-02_wp,-6.79035487613655E-04_wp, &
      & 8.78198962890330E-03_wp, 9.72829001385098E-04_wp, 4.31613119157574E-04_wp, &
      & 3.39329532587728E-04_wp,-1.10362577227696E-04_wp, 8.20860388583884E-06_wp, &
      &-2.86085329171831E-04_wp,-1.44661987627843E-04_wp,-8.16969303104830E-05_wp, &
      &-4.08503618156095E-04_wp,-6.88819688392179E-04_wp,-6.06690581183132E-04_wp, &
      & 9.38622543052269E-04_wp,-7.29860323332650E-05_wp,-1.45840439450576E-03_wp, &
      &-5.41519705143916E-04_wp,-1.88651521176194E-03_wp, 3.60304872458680E-03_wp, &
      & 4.66111281358563E-04_wp, 4.14351970344133E-03_wp, 3.56893171984317E-05_wp, &
      &-2.26084037291830E-03_wp, 5.13537734107091E-04_wp, 2.23194762526908E-03_wp, &
      & 1.40000483980811E-04_wp,-1.02717741612734E-03_wp,-1.67720095156313E-05_wp, &
      & 1.29518345211162E-04_wp, 3.57114046695395E-05_wp, 1.50166180255650E-05_wp, &
      &-1.21929197698885E-04_wp,-6.52194452940651E-04_wp,-5.15623780968630E-04_wp, &
      & 2.73294147837168E-04_wp,-9.11286541924727E-04_wp, 2.60909585424519E-05_wp, &
      &-2.51838077169850E-05_wp,-2.53947944807633E-05_wp, 6.84526285750054E-05_wp, &
      &-7.22055400537447E-05_wp, 8.40503347739840E-04_wp,-1.32617256121615E-02_wp, &
      & 7.65027028289048E-03_wp, 7.29973779501288E-04_wp, 4.31613119157574E-04_wp, &
      & 3.07330023331940E-04_wp,-6.96308160686061E-05_wp,-2.72899832828417E-04_wp, &
      & 1.35850391213342E-04_wp,-5.83829373776639E-04_wp,-2.20905355380273E-05_wp, &
      &-5.05453807147833E-04_wp,-1.59367877771540E-04_wp, 3.53857312562318E-04_wp, &
      & 6.06034207881716E-04_wp,-1.22901953593348E-03_wp,-2.39582887844565E-03_wp, &
      & 1.08123286916086E-04_wp, 1.74128670059436E-04_wp, 2.10285595586614E-03_wp, &
      &-3.54139161533717E-03_wp, 1.29891259804992E-03_wp,-5.00522584021680E-03_wp, &
      &-2.34151302881454E-03_wp, 4.94952180334675E-04_wp, 8.26572953275499E-04_wp, &
      & 2.14969208354250E-03_wp,-1.14187690596524E-03_wp,-7.37519467381799E-04_wp, &
      & 6.88550780578720E-05_wp,-6.06645445272479E-05_wp,-1.37991308366822E-04_wp, &
      &-9.17945277861075E-05_wp,-1.22852682598501E-04_wp,-1.24575925065985E-03_wp, &
      & 1.97029419195767E-03_wp, 2.40805559719466E-04_wp,-5.17440346509960E-04_wp, &
      & 9.82423307899891E-05_wp,-5.85408576231827E-05_wp, 7.94171880800760E-05_wp, &
      &-1.17897739759140E-05_wp, 1.11585039120817E-04_wp, 1.56038857676967E-03_wp, &
      &-5.84858535888869E-03_wp,-1.49468871725830E-02_wp, 1.48135170964987E-02_wp, &
      & 3.39329532587728E-04_wp,-6.96308160686061E-05_wp, 5.91760718215047E-04_wp, &
      & 2.52912255874616E-04_wp,-4.23264327989667E-04_wp, 1.06473924707411E-03_wp, &
      &-4.45102163842455E-04_wp, 3.80106625809516E-05_wp, 1.00614012061061E-03_wp, &
      &-1.40694498187750E-03_wp, 1.10660155867559E-04_wp,-2.36862461809407E-03_wp, &
      & 4.78040749239259E-04_wp, 2.31640389940747E-03_wp,-2.34898043484794E-03_wp, &
      & 5.97064769855728E-04_wp,-2.95929713103497E-03_wp,-8.75032383977239E-04_wp, &
      &-6.01050825635358E-03_wp, 4.72722357649428E-04_wp, 2.52712586788310E-03_wp, &
      & 7.14219643150120E-04_wp, 2.26334563100816E-04_wp, 2.76159075773128E-04_wp, &
      & 2.37177085329030E-03_wp, 1.49459436391235E-04_wp,-3.58157141419361E-05_wp, &
      &-1.10652382975503E-04_wp, 3.21446802815042E-05_wp, 4.05393279044262E-05_wp, &
      &-1.84891013913080E-03_wp, 3.39979714316987E-03_wp, 3.98704788703455E-03_wp, &
      & 2.69307620534643E-03_wp,-1.33251711439347E-05_wp,-2.31098244614289E-04_wp, &
      &-5.67477088132154E-05_wp, 2.91415636443471E-05_wp, 2.53643008087482E-04_wp, &
      &-3.26760910351179E-03_wp,-9.24621410537155E-04_wp,-2.63446016222323E-02_wp, &
      &-1.41265088250698E-02_wp,-1.10362577227696E-04_wp,-2.72899832828417E-04_wp, &
      & 2.52912255874616E-04_wp, 1.08310028976880E-03_wp, 4.38482863598217E-04_wp, &
      & 1.20864032142143E-03_wp,-5.21718641131263E-04_wp, 4.70419178273609E-04_wp, &
      & 6.42501172500220E-04_wp,-2.36540545013463E-03_wp,-1.13219077295460E-03_wp, &
      & 3.42392873702847E-05_wp, 2.59547000701250E-03_wp, 5.16614467010584E-04_wp, &
      &-2.72927683099058E-03_wp,-3.32760114530298E-03_wp, 4.28984898099735E-03_wp, &
      &-1.61471980448458E-03_wp, 3.53071480903333E-03_wp, 2.70282185999502E-03_wp, &
      &-5.02041178157069E-04_wp, 1.14929051465742E-03_wp,-2.69348432554597E-04_wp, &
      & 1.54213587967972E-03_wp, 1.59044196104693E-03_wp, 3.48745958377767E-05_wp, &
      & 1.47884993198466E-04_wp, 8.18569440148499E-05_wp, 1.29877997002024E-04_wp, &
      & 4.47136744943953E-05_wp,-3.95093353419281E-04_wp, 1.57372281668683E-04_wp, &
      & 3.83007723351944E-03_wp, 1.58235914556952E-03_wp,-7.80121783922660E-05_wp, &
      &-1.74472472563144E-04_wp,-1.57956030126034E-04_wp, 1.16342256034267E-04_wp, &
      & 2.63123070323647E-05_wp,-3.45211447044840E-03_wp,-8.64305835007620E-03_wp, &
      &-1.04417300489744E-03_wp,-2.63165868327967E-02_wp, 8.20860388583884E-06_wp, &
      & 1.35850391213342E-04_wp,-4.23264327989667E-04_wp, 4.38482863598217E-04_wp, &
      & 9.40813877255382E-04_wp, 7.11952239175541E-03_wp,-4.48759433563640E-04_wp, &
      & 7.04758808977037E-04_wp,-2.94624945675146E-04_wp,-2.74166581356118E-03_wp, &
      & 6.84153465016012E-03_wp,-5.31279761822284E-03_wp, 4.85658416213751E-03_wp, &
      & 9.50019121811072E-04_wp,-2.43299151738347E-03_wp,-9.06679526463687E-03_wp, &
      &-3.80607531989747E-04_wp,-2.37874651841234E-03_wp,-3.06390494210974E-04_wp, &
      & 3.36294433121599E-03_wp, 8.99503672589565E-03_wp,-1.11688855965280E-03_wp, &
      & 3.20046292353626E-02_wp, 2.32606695048069E-02_wp, 2.10966188523550E-02_wp, &
      & 1.14006433275942E-03_wp, 1.14147229303161E-03_wp,-3.02272998088171E-04_wp, &
      & 7.61032562708509E-04_wp,-4.79144450852202E-04_wp,-1.30082860861911E-03_wp, &
      &-3.51696354973918E-02_wp,-1.39407590739072E-02_wp, 2.24022500427104E-02_wp, &
      &-1.30001097214204E-03_wp, 9.03676498870888E-04_wp,-5.83171540616180E-04_wp, &
      &-5.86376985234874E-04_wp,-5.20172683811180E-04_wp,-1.20825923162132E-03_wp, &
      & 1.01661826271289E-02_wp,-1.38165545028930E-02_wp,-4.08742499801316E-02_wp, &
      &-7.18067691148890E-04_wp,-2.63641797653816E-04_wp,-5.93966433083350E-04_wp, &
      & 1.02979125404260E-03_wp, 1.20885516832022E-03_wp,-4.48759433563640E-04_wp, &
      & 2.28687129134323E-03_wp, 5.41396449947113E-04_wp,-2.59430166947036E-04_wp, &
      &-2.64919936937197E-03_wp, 5.77530559995243E-04_wp,-5.26976013453336E-05_wp, &
      &-1.70360657293970E-03_wp,-2.12815610708007E-04_wp, 1.06766405478897E-03_wp, &
      & 1.23297231618397E-03_wp,-1.24027234557057E-03_wp,-2.62975652749672E-03_wp, &
      & 6.81819221360484E-04_wp,-2.77984913399787E-03_wp, 6.15223438264382E-04_wp, &
      & 7.17489050802568E-03_wp,-1.15284447884864E-02_wp,-1.66995151518811E-02_wp, &
      &-1.96285390992294E-02_wp,-8.09390428732436E-04_wp,-6.27932511020789E-04_wp, &
      & 1.21450766756593E-04_wp,-6.28894926244619E-04_wp,-7.62484937972927E-05_wp, &
      &-8.00507020681381E-03_wp,-1.18379012605945E-02_wp,-1.74855485990918E-02_wp, &
      & 2.22060923985063E-02_wp,-8.41884382085804E-04_wp, 6.54018323036504E-04_wp, &
      &-9.79735922629530E-05_wp,-7.38014920537504E-04_wp, 1.67408952638742E-04_wp, &
      & 2.55413860729610E-03_wp, 1.44077180987701E-02_wp, 4.81975145073823E-03_wp, &
      & 1.15402494798051E-02_wp,-3.91472860045293E-04_wp,-1.85447400076972E-04_wp, &
      & 2.69372961745411E-05_wp,-3.83960703078049E-04_wp,-5.52056921452457E-04_wp, &
      & 7.04758808977037E-04_wp, 5.41396449947113E-04_wp, 1.62716519681099E-03_wp, &
      & 3.61628176171863E-04_wp,-1.63409127414158E-03_wp,-1.07250498654332E-03_wp, &
      & 3.31847658399270E-04_wp,-6.83598102194624E-04_wp, 5.75709000782809E-04_wp, &
      &-1.94547758751706E-04_wp,-6.02151183303941E-04_wp, 1.49179954780895E-03_wp, &
      &-2.93152544359318E-03_wp, 1.04205333697131E-03_wp, 1.48481072358356E-04_wp, &
      & 6.08528200623643E-04_wp, 4.57690583399948E-03_wp,-1.88737263402908E-02_wp, &
      & 8.75664181962437E-03_wp,-1.26486079338156E-02_wp,-7.10043290946883E-04_wp, &
      & 2.89620397186540E-06_wp, 6.06131016453156E-04_wp, 2.82293389951878E-06_wp, &
      & 2.86985578991552E-04_wp,-3.80262208043235E-03_wp,-1.56986199671817E-02_wp, &
      & 9.95289435444545E-03_wp, 1.08109379474327E-02_wp,-6.58426823850392E-04_wp, &
      &-2.72580344706362E-05_wp,-5.84203711930900E-04_wp, 1.67615675300633E-05_wp, &
      &-1.95672888720847E-04_wp,-3.63428057490289E-03_wp, 5.13544987508727E-03_wp, &
      & 1.02100793239385E-02_wp,-1.81555166175954E-02_wp,-3.91264029795281E-04_wp, &
      & 9.88365504774939E-06_wp,-5.75792938733305E-04_wp,-4.71112289724211E-05_wp, &
      & 5.46257757174961E-04_wp,-2.94624945675146E-04_wp,-2.59430166947036E-04_wp, &
      & 3.61628176171863E-04_wp, 2.43622668168911E-03_wp,-1.06123973390391E-03_wp, &
      &-1.65324946777580E-03_wp,-7.25296306812267E-05_wp, 1.80304070856627E-03_wp, &
      & 2.13517632377509E-03_wp,-1.81460080347582E-03_wp, 2.37543992326638E-03_wp, &
      & 6.76725101117473E-04_wp,-1.81030095889503E-03_wp,-1.71651898481580E-03_wp, &
      & 2.34703941229104E-03_wp,-2.97757430264917E-04_wp, 4.75605525498004E-03_wp, &
      &-1.99507322125151E-02_wp,-1.13575027794183E-02_wp, 4.49366063601879E-03_wp, &
      &-2.46172687351467E-04_wp,-6.32133131700029E-04_wp, 8.27214390433664E-05_wp, &
      &-1.20209486448796E-04_wp, 7.01975281929019E-04_wp, 5.46887348152095E-03_wp, &
      & 2.07465098916424E-02_wp, 1.12870324819581E-02_wp, 2.62360032992966E-03_wp, &
      & 3.30170019630296E-04_wp,-6.85792113938563E-04_wp, 6.13175091245107E-05_wp, &
      & 2.02457181360983E-04_wp, 6.65644753860528E-04_wp,-9.08181098067992E-03_wp, &
      & 1.33513181820185E-02_wp,-1.99268786228131E-02_wp,-2.34217202320817E-02_wp, &
      &-6.69813653944916E-04_wp,-4.38712775661301E-04_wp,-1.20006019567996E-04_wp, &
      & 1.02053886193382E-03_wp, 5.94007597988431E-04_wp,-2.74166581356118E-03_wp, &
      &-2.64919936937197E-03_wp,-1.63409127414158E-03_wp,-1.06123973390391E-03_wp, &
      & 4.26660636479625E-02_wp, 5.48867639267580E-03_wp, 2.54176254686327E-03_wp, &
      & 7.36832296027571E-03_wp,-4.03975654473346E-04_wp,-3.26942203748600E-04_wp, &
      & 3.91589546659811E-03_wp, 4.64342996080713E-03_wp, 1.33992841815513E-02_wp, &
      & 7.49785230432728E-04_wp, 8.17240065738855E-04_wp,-3.19238355445102E-04_wp, &
      & 2.23646244137191E-02_wp, 6.59548065857688E-02_wp, 8.68754738499332E-02_wp, &
      & 1.25040456174465E-03_wp, 8.43040832431947E-04_wp, 2.57139294590380E-03_wp, &
      & 4.69721320177037E-04_wp, 1.34881283254919E-03_wp,-1.20147515573047E-03_wp, &
      &-2.39915293849376E-02_wp, 7.91814684776186E-02_wp, 6.88501216429309E-02_wp, &
      &-1.27801136787665E-02_wp, 1.67318715111030E-03_wp,-2.26525736225413E-03_wp, &
      & 1.03312589059841E-04_wp, 1.25511859423494E-03_wp, 1.23392237280277E-03_wp, &
      &-1.35620731019103E-02_wp, 7.56618991833905E-02_wp, 3.15084854540213E-02_wp, &
      & 6.45053437127255E-02_wp,-9.38359904068336E-04_wp,-6.34043613156057E-04_wp, &
      & 3.14862311605364E-04_wp,-1.37787803489053E-03_wp,-2.24225399100506E-03_wp, &
      & 6.84153465016012E-03_wp, 5.77530559995243E-04_wp,-1.07250498654332E-03_wp, &
      &-1.65324946777580E-03_wp, 5.48867639267580E-03_wp, 3.04783663710076E-02_wp, &
      & 5.04606442095931E-03_wp,-4.66371208138162E-03_wp,-7.19052804954136E-03_wp, &
      & 3.26384647666834E-03_wp,-1.16648070070522E-02_wp,-1.74839240093770E-03_wp, &
      & 1.30572623509032E-03_wp,-5.03897420550643E-03_wp,-8.96763548844018E-04_wp, &
      & 9.84759792548253E-03_wp, 2.07224740193056E-02_wp, 6.70688097195501E-02_wp, &
      &-1.96459370841132E-02_wp, 7.81809041418510E-02_wp, 2.84420597616094E-03_wp, &
      & 3.81213248123848E-04_wp,-1.93389141426929E-03_wp, 4.82420477979836E-04_wp, &
      &-5.32810176688890E-04_wp, 1.83968711176338E-02_wp,-6.46585082511956E-02_wp, &
      & 3.76360397992007E-02_wp, 6.76726646968229E-02_wp,-2.56333530435437E-03_wp, &
      & 3.57609160919946E-04_wp,-2.02312098196746E-03_wp,-3.47311396600001E-04_wp, &
      &-5.79280405069388E-04_wp,-5.21056540703202E-03_wp, 3.66511396790250E-02_wp, &
      &-1.61959501492248E-02_wp, 4.22610487435240E-02_wp,-7.23150166635905E-04_wp, &
      &-6.96300015777383E-04_wp, 6.92000783727210E-04_wp, 2.14568368478059E-04_wp, &
      &-1.19565402301583E-03_wp,-5.31279761822284E-03_wp,-5.26976013453336E-05_wp, &
      & 3.31847658399270E-04_wp,-7.25296306812267E-05_wp, 2.54176254686327E-03_wp, &
      & 5.04606442095931E-03_wp, 2.84679587556288E-02_wp, 3.54894557191667E-03_wp, &
      &-7.72046463325575E-04_wp, 3.39812061949490E-03_wp, 1.80295732067456E-03_wp, &
      &-1.72699146503989E-03_wp, 2.08068962045653E-04_wp,-9.57025303578533E-04_wp, &
      &-6.72245961471300E-04_wp,-1.25000483023922E-02_wp,-7.48190813025679E-03_wp, &
      & 1.50194621101878E-02_wp,-9.51746257217116E-02_wp, 1.07105823451813E-02_wp, &
      & 2.74276798196795E-04_wp,-2.28477526111468E-03_wp,-1.54670546994734E-03_wp, &
      &-1.52509696471363E-03_wp,-8.34604782021031E-05_wp,-9.70631584445533E-03_wp, &
      & 5.89106335179422E-04_wp, 8.54634077370166E-02_wp,-6.23156571113277E-03_wp, &
      & 2.39237960359928E-04_wp,-1.86157192226547E-03_wp,-1.22914115729581E-03_wp, &
      & 1.62380339070335E-03_wp,-2.57899109827770E-04_wp,-9.88014923090812E-03_wp, &
      &-6.25144814937098E-03_wp, 8.60612099656798E-02_wp, 3.14099372304354E-03_wp, &
      & 3.78532516914578E-04_wp, 9.38872653962262E-04_wp,-1.19990637124721E-03_wp, &
      &-2.29785457537120E-03_wp, 4.74724238068227E-05_wp, 4.85658416213751E-03_wp, &
      &-1.70360657293970E-03_wp,-6.83598102194624E-04_wp, 1.80304070856627E-03_wp, &
      & 7.36832296027571E-03_wp,-4.66371208138162E-03_wp, 3.54894557191667E-03_wp, &
      & 3.34393781938834E-02_wp, 4.73561694058205E-03_wp,-9.20066557038393E-03_wp, &
      &-4.78340485526953E-03_wp,-5.05261443662532E-03_wp, 6.84791172569323E-04_wp, &
      & 2.03099825523925E-03_wp, 9.45917682510677E-03_wp, 3.78879388099240E-03_wp, &
      & 1.40822121394520E-02_wp, 7.85903417261112E-02_wp,-1.34573967461192E-02_wp, &
      & 5.53985700081892E-03_wp, 1.51665855913107E-03_wp, 4.92005695112074E-04_wp, &
      &-1.33814279600227E-03_wp, 9.60493988868284E-06_wp,-1.49543254126587E-03_wp, &
      &-1.19162210963841E-02_wp, 7.77031024200781E-02_wp,-3.00536616447836E-02_wp, &
      &-2.00732274112082E-03_wp, 1.10329669950774E-03_wp,-1.68670894927570E-04_wp, &
      & 1.46237940335649E-03_wp,-6.46986175422616E-04_wp, 1.68257298559226E-03_wp, &
      & 2.13600011676250E-02_wp, 3.25281161052200E-02_wp, 4.67381503378004E-02_wp, &
      &-1.03858523700384E-01_wp,-1.50740974288056E-03_wp, 3.82651994861127E-05_wp, &
      &-2.40313487407472E-03_wp, 3.71631530998101E-04_wp, 2.60157115829106E-03_wp, &
      & 9.50019121811072E-04_wp,-2.12815610708007E-04_wp, 5.75709000782809E-04_wp, &
      & 2.13517632377509E-03_wp,-4.03975654473346E-04_wp,-7.19052804954136E-03_wp, &
      &-7.72046463325575E-04_wp, 4.73561694058205E-03_wp, 4.18770813681395E-02_wp, &
      &-2.50290440454828E-03_wp, 7.65727418568397E-04_wp, 2.12624345393224E-03_wp, &
      &-4.97027877918030E-03_wp,-5.83472669844520E-03_wp, 5.51732011957694E-03_wp, &
      &-4.27826077944579E-04_wp,-8.82403162605848E-03_wp, 3.37743473393132E-02_wp, &
      &-3.36885096134573E-02_wp,-8.78793958800221E-02_wp,-1.20767003906758E-03_wp, &
      &-4.48556039119347E-04_wp,-1.81086224544769E-04_wp,-1.30225149564277E-03_wp, &
      &-1.72564923333505E-03_wp,-8.62878718054140E-03_wp,-3.71725125852140E-02_wp, &
      & 3.53696185434452E-02_wp,-9.55348684789663E-02_wp, 1.04090454647008E-03_wp, &
      &-3.59221400041849E-04_wp,-3.57354903215663E-04_wp, 1.79429741173201E-03_wp, &
      &-2.19262341202439E-03_wp, 2.19832369216965E-02_wp, 7.20780368313875E-02_wp, &
      &-7.10320337351810E-02_wp,-4.72446672147925E-02_wp,-2.04798561893996E-03_wp, &
      &-1.49699282601155E-03_wp, 1.86593826628628E-04_wp, 2.32423304833064E-03_wp, &
      & 4.22160666294196E-04_wp,-2.43299151738347E-03_wp, 1.06766405478897E-03_wp, &
      &-1.94547758751706E-04_wp,-1.81460080347582E-03_wp,-3.26942203748600E-04_wp, &
      & 3.26384647666834E-03_wp, 3.39812061949490E-03_wp,-9.20066557038393E-03_wp, &
      &-2.50290440454828E-03_wp, 1.18831170923295E-02_wp, 3.11417561356810E-03_wp, &
      &-1.79633667273998E-03_wp, 1.08428680718719E-03_wp, 2.95504975538360E-03_wp, &
      & 9.40583767874230E-04_wp,-8.18360385983953E-04_wp, 2.70527565876776E-05_wp, &
      &-3.42955794705007E-02_wp,-1.00765000114553E-02_wp, 3.53152003767956E-02_wp, &
      & 3.49311376473168E-04_wp,-7.42699380839901E-04_wp, 4.39431782533191E-06_wp, &
      & 3.27030474840371E-04_wp, 1.58749721679993E-03_wp, 5.59939284194209E-04_wp, &
      &-4.33600687024431E-02_wp,-1.57001763492994E-02_wp,-3.92643295028350E-02_wp, &
      & 4.03441445098017E-04_wp, 1.24040844923198E-03_wp, 2.73735765525016E-04_wp, &
      & 2.82584056419365E-04_wp,-2.02260087963467E-03_wp,-8.83664108598077E-04_wp, &
      & 3.13114544459522E-02_wp, 3.67817089890598E-02_wp, 5.10642880072331E-02_wp, &
      &-6.93371177426215E-04_wp,-1.52128659470300E-04_wp,-6.60971351199133E-05_wp, &
      &-2.04580162271769E-03_wp,-1.86857557022708E-03_wp,-9.06679526463687E-03_wp, &
      & 1.23297231618397E-03_wp,-6.02151183303941E-04_wp, 2.37543992326638E-03_wp, &
      & 3.91589546659811E-03_wp,-1.16648070070522E-02_wp, 1.80295732067456E-03_wp, &
      &-4.78340485526953E-03_wp, 7.65727418568397E-04_wp, 3.11417561356810E-03_wp, &
      & 1.53506396343196E-02_wp, 4.93458019576049E-04_wp, 1.16823409155536E-03_wp, &
      & 3.13052003010560E-04_wp,-1.86532900379958E-03_wp,-8.96110846963001E-03_wp, &
      &-9.72183490348629E-04_wp,-6.96677285675238E-02_wp,-2.42640812741514E-02_wp, &
      &-2.43105967018234E-02_wp,-1.88909058022195E-03_wp,-1.84620872265765E-03_wp, &
      & 8.40605401675080E-04_wp,-8.31609991216725E-04_wp, 1.57654392909966E-03_wp, &
      &-6.32238839106823E-04_wp, 6.52028457208073E-02_wp, 1.09811119079166E-04_wp, &
      &-2.52064987261363E-02_wp, 1.91383514382521E-03_wp,-1.05856014409741E-03_wp, &
      & 1.27412333131700E-03_wp, 3.48350771197910E-04_wp, 1.31355340280529E-03_wp, &
      & 5.25854327632037E-05_wp, 2.45048495855122E-02_wp,-2.56655009577307E-03_wp, &
      & 4.06754362288397E-02_wp,-4.24473085272663E-04_wp,-4.73872134127600E-04_wp, &
      & 6.75497357337969E-04_wp,-5.61485634191755E-04_wp,-1.60968883271744E-03_wp, &
      &-3.80607531989747E-04_wp,-1.24027234557057E-03_wp, 1.49179954780895E-03_wp, &
      & 6.76725101117473E-04_wp, 4.64342996080713E-03_wp,-1.74839240093770E-03_wp, &
      &-1.72699146503989E-03_wp,-5.05261443662532E-03_wp, 2.12624345393224E-03_wp, &
      &-1.79633667273998E-03_wp, 4.93458019576049E-04_wp, 5.57899823592994E-03_wp, &
      &-1.29897186515243E-04_wp,-7.80032096532300E-04_wp,-1.84892158652279E-04_wp, &
      &-9.62092557670708E-04_wp, 3.41667470649450E-05_wp,-2.14598915207260E-02_wp, &
      & 4.22313704080252E-02_wp,-1.51574111619879E-02_wp,-7.82803317697228E-04_wp, &
      & 8.27111347824794E-04_wp, 1.17930469715137E-03_wp, 5.54949096390952E-04_wp, &
      & 3.31242618197280E-04_wp,-2.19959469457463E-04_wp,-3.87427809377845E-03_wp, &
      & 4.92359767360659E-02_wp, 3.76026077550443E-04_wp,-1.50246094084492E-04_wp, &
      &-1.27777504135482E-03_wp,-1.02227177755345E-03_wp, 9.55873747394603E-04_wp, &
      &-4.76745213909973E-05_wp, 5.49263377585196E-04_wp,-6.37450174955805E-03_wp, &
      &-1.46067658017428E-02_wp,-4.95983018633171E-03_wp, 1.86494551430705E-04_wp, &
      &-2.30162533491128E-05_wp, 1.12528083885423E-04_wp, 5.83275218807406E-04_wp, &
      & 4.19280315359715E-04_wp,-2.37874651841234E-03_wp,-2.62975652749672E-03_wp, &
      &-2.93152544359318E-03_wp,-1.81030095889503E-03_wp, 1.33992841815513E-02_wp, &
      & 1.30572623509032E-03_wp, 2.08068962045653E-04_wp, 6.84791172569323E-04_wp, &
      &-4.97027877918030E-03_wp, 1.08428680718719E-03_wp, 1.16823409155536E-03_wp, &
      &-1.29897186515243E-04_wp, 9.94428221843365E-03_wp,-7.19846587009332E-05_wp, &
      &-1.93960408498675E-04_wp,-2.74856880867014E-03_wp, 1.27031368525090E-03_wp, &
      & 4.18745285287946E-02_wp, 2.72966013888811E-02_wp, 2.85915868985803E-02_wp, &
      & 1.52421575467949E-03_wp, 1.17133429597765E-03_wp,-5.70259377186291E-04_wp, &
      & 8.00501008899362E-04_wp,-5.93648368390723E-04_wp,-8.89167773593586E-04_wp, &
      & 4.01388520677411E-02_wp, 6.11538319289382E-03_wp,-3.43827040288937E-02_wp, &
      & 1.81283369635852E-03_wp,-6.03244529432564E-04_wp, 9.60330421289763E-04_wp, &
      & 6.66961474578803E-04_wp, 1.38916670200151E-04_wp,-8.40001714014691E-04_wp, &
      &-1.91461291146198E-02_wp, 5.62002262438386E-03_wp, 4.94693462764905E-02_wp, &
      & 1.30712345342734E-03_wp, 4.10116747608054E-04_wp, 9.58776573157139E-04_wp, &
      &-7.52025031398366E-04_wp,-1.22301380113740E-03_wp,-3.06390494210974E-04_wp, &
      & 6.81819221360484E-04_wp, 1.04205333697131E-03_wp,-1.71651898481580E-03_wp, &
      & 7.49785230432728E-04_wp,-5.03897420550643E-03_wp,-9.57025303578533E-04_wp, &
      & 2.03099825523925E-03_wp,-5.83472669844520E-03_wp, 2.95504975538360E-03_wp, &
      & 3.13052003010560E-04_wp,-7.80032096532300E-04_wp,-7.19846587009332E-05_wp, &
      & 6.02218778760986E-03_wp,-1.96369388860250E-04_wp, 4.00523412010053E-04_wp, &
      & 2.43415200592200E-05_wp,-1.51845049947363E-02_wp, 2.87327242622429E-02_wp, &
      &-9.49249765705375E-03_wp,-5.48713980681319E-04_wp, 5.54022846553796E-04_wp, &
      & 8.01687607586127E-04_wp, 3.77777202539428E-04_wp, 1.95523826272871E-04_wp, &
      & 6.68905535755028E-04_wp,-5.70282419272425E-03_wp,-3.25394854692952E-02_wp, &
      &-4.86038808135031E-03_wp,-6.16005794354335E-05_wp, 1.02574808923692E-03_wp, &
      & 4.82756606241590E-04_wp,-5.53868614428592E-04_wp,-4.61704646589813E-04_wp, &
      &-4.42109441645743E-04_wp, 1.09535528037961E-03_wp, 5.78076943435148E-02_wp, &
      &-1.68594875769831E-03_wp,-1.50326172865105E-04_wp, 5.22664015461809E-04_wp, &
      &-1.09131629766252E-03_wp,-1.82188558275370E-03_wp,-3.25003079162797E-05_wp, &
      & 3.36294433121599E-03_wp,-2.77984913399787E-03_wp, 1.48481072358356E-04_wp, &
      & 2.34703941229104E-03_wp, 8.17240065738855E-04_wp,-8.96763548844018E-04_wp, &
      &-6.72245961471300E-04_wp, 9.45917682510677E-03_wp, 5.51732011957694E-03_wp, &
      & 9.40583767874230E-04_wp,-1.86532900379958E-03_wp,-1.84892158652279E-04_wp, &
      &-1.93960408498675E-04_wp,-1.96369388860250E-04_wp, 1.14272521505618E-02_wp, &
      & 3.81891246340949E-03_wp, 4.46942267787996E-04_wp,-4.11175285560145E-03_wp, &
      & 9.67588448002621E-03_wp, 5.64940155335164E-02_wp, 1.60503540214095E-03_wp, &
      & 1.52655440518585E-04_wp,-3.37131360832134E-04_wp, 1.19451943984326E-03_wp, &
      & 1.48935194145954E-03_wp,-9.31069490873247E-04_wp, 7.87169218444709E-03_wp, &
      & 3.89347577592347E-03_wp,-4.72243223724468E-02_wp, 1.35155308784679E-03_wp, &
      &-2.00017008301805E-04_wp, 5.23748289028198E-04_wp, 8.65609593211044E-04_wp, &
      &-8.12160921508423E-04_wp, 1.03952145445938E-03_wp, 3.42404939671757E-02_wp, &
      &-2.44882129383753E-03_wp,-6.15913706401558E-02_wp,-1.62689399094736E-03_wp, &
      &-5.53268037298817E-04_wp,-1.20341681703315E-03_wp, 1.03007344334789E-03_wp, &
      & 1.53200056222922E-03_wp, 8.99503672589565E-03_wp, 6.15223438264382E-04_wp, &
      & 6.08528200623643E-04_wp,-2.97757430264917E-04_wp,-3.19238355445102E-04_wp, &
      & 9.84759792548253E-03_wp,-1.25000483023922E-02_wp, 3.78879388099240E-03_wp, &
      &-4.27826077944579E-04_wp,-8.18360385983953E-04_wp,-8.96110846963001E-03_wp, &
      &-9.62092557670708E-04_wp,-2.74856880867014E-03_wp, 4.00523412010053E-04_wp, &
      & 3.81891246340949E-03_wp, 1.48383546641339E-02_wp, 2.14569885078231E-04_wp, &
      & 2.60228009071645E-02_wp, 4.55480551778791E-02_wp, 3.11414586934253E-02_wp, &
      & 1.24763488181786E-03_wp, 1.53320635726685E-03_wp, 7.48314601602185E-08_wp, &
      & 1.22592371170185E-03_wp,-1.03037430243343E-04_wp, 6.85864626230428E-04_wp, &
      &-3.93198699151506E-02_wp,-4.33959140311928E-02_wp, 3.84544626969289E-02_wp, &
      &-1.77253159605633E-03_wp, 1.76953567734680E-03_wp,-2.92690161508324E-04_wp, &
      &-1.42187034413839E-03_wp,-3.56619680821674E-04_wp, 2.40006331976448E-04_wp, &
      & 4.74589429392366E-02_wp,-2.86025601155137E-02_wp,-2.74869037414883E-02_wp, &
      &-1.81297319011831E-03_wp,-9.90759348859535E-04_wp,-3.97806973349569E-04_wp, &
      & 1.22820555419500E-03_wp, 4.56941288207854E-04_wp,-1.11688855965280E-03_wp, &
      & 7.17489050802568E-03_wp, 4.57690583399948E-03_wp, 4.75605525498004E-03_wp, &
      & 2.23646244137191E-02_wp, 2.07224740193056E-02_wp,-7.48190813025679E-03_wp, &
      & 1.40822121394520E-02_wp,-8.82403162605848E-03_wp, 2.70527565876776E-05_wp, &
      &-9.72183490348629E-04_wp, 3.41667470649450E-05_wp, 1.27031368525090E-03_wp, &
      & 2.43415200592200E-05_wp, 4.46942267787996E-04_wp, 2.14569885078231E-04_wp, &
      & 9.79953185022657E-01_wp,-6.56913535165536E-02_wp,-4.15746208633436E-02_wp, &
      &-4.42906788242026E-02_wp,-3.26497777863993E-03_wp,-3.02957050163202E-03_wp, &
      & 1.03160173383662E-03_wp,-2.06224126752325E-03_wp, 1.29590186401136E-03_wp, &
      &-1.73578590097926E-03_wp,-4.78755981358135E-03_wp,-4.80964039210276E-03_wp, &
      &-2.26550292609627E-02_wp, 1.03244385712790E-03_wp, 8.28645456652537E-04_wp, &
      & 8.36183698869665E-04_wp, 3.56016852022877E-04_wp,-2.26721397739007E-03_wp, &
      &-4.67936049592108E-04_wp,-2.15316812457990E-02_wp,-5.04588128389885E-03_wp, &
      & 2.25098483451193E-03_wp, 2.05812072051282E-03_wp, 5.05766248283295E-04_wp, &
      & 7.91454294590884E-04_wp, 7.16816433450587E-04_wp, 1.18902455312498E-03_wp, &
      & 3.20046292353626E-02_wp,-1.15284447884864E-02_wp,-1.88737263402908E-02_wp, &
      &-1.99507322125151E-02_wp, 6.59548065857688E-02_wp, 6.70688097195501E-02_wp, &
      & 1.50194621101878E-02_wp, 7.85903417261112E-02_wp, 3.37743473393132E-02_wp, &
      &-3.42955794705007E-02_wp,-6.96677285675238E-02_wp,-2.14598915207260E-02_wp, &
      & 4.18745285287946E-02_wp,-1.51845049947363E-02_wp,-4.11175285560145E-03_wp, &
      & 2.60228009071645E-02_wp,-6.56913535165536E-02_wp, 8.97832343104443E-01_wp, &
      &-9.01190569007232E-03_wp,-1.45144645636867E-02_wp, 1.55005079704011E-02_wp, &
      & 1.20367234816803E-02_wp,-1.34556376466504E-02_wp,-3.44948217512164E-04_wp, &
      &-2.38477695614008E-02_wp, 5.97412493903298E-03_wp, 3.58652467089654E-03_wp, &
      &-1.13563354350246E-02_wp,-9.59168996535191E-03_wp,-6.71551375241214E-04_wp, &
      & 2.30409327479312E-03_wp, 4.33495391693635E-04_wp,-1.03517926424556E-03_wp, &
      &-9.32392910000516E-04_wp,-1.75636665246333E-02_wp,-4.80394095409376E-02_wp, &
      &-4.04081526462394E-03_wp, 2.91067062790983E-02_wp, 5.41082241160368E-03_wp, &
      & 2.05070858833928E-03_wp, 2.29339909979128E-03_wp, 4.07757821937252E-04_wp, &
      &-4.27488088447255E-04_wp, 2.32606695048069E-02_wp,-1.66995151518811E-02_wp, &
      & 8.75664181962437E-03_wp,-1.13575027794183E-02_wp, 8.68754738499332E-02_wp, &
      &-1.96459370841132E-02_wp,-9.51746257217116E-02_wp,-1.34573967461192E-02_wp, &
      &-3.36885096134573E-02_wp,-1.00765000114553E-02_wp,-2.42640812741514E-02_wp, &
      & 4.22313704080252E-02_wp, 2.72966013888811E-02_wp, 2.87327242622429E-02_wp, &
      & 9.67588448002621E-03_wp, 4.55480551778791E-02_wp,-4.15746208633436E-02_wp, &
      &-9.01190569007232E-03_wp, 8.98949522434251E-01_wp,-5.85249818311784E-03_wp, &
      &-4.20396123671864E-05_wp, 2.32696964109612E-02_wp, 1.40699101827079E-02_wp, &
      & 1.57206108208992E-02_wp, 4.82043997276390E-05_wp, 1.29832452985145E-03_wp, &
      &-1.01189621165404E-02_wp, 1.64501682042580E-02_wp,-1.13646125614409E-02_wp, &
      &-6.58322439410420E-04_wp,-9.95021813203161E-05_wp,-1.69154717257233E-03_wp, &
      & 8.98861322776528E-04_wp,-2.17156562704799E-03_wp, 1.74725131414713E-03_wp, &
      &-1.23060391666262E-02_wp, 1.53478173145738E-02_wp,-5.75999068005323E-03_wp, &
      & 6.00827861785183E-04_wp, 6.76325363903052E-04_wp,-1.66014282991117E-03_wp, &
      &-2.79625792410964E-04_wp, 2.06014446932532E-03_wp, 2.10966188523550E-02_wp, &
      &-1.96285390992294E-02_wp,-1.26486079338156E-02_wp, 4.49366063601879E-03_wp, &
      & 1.25040456174465E-03_wp, 7.81809041418510E-02_wp, 1.07105823451813E-02_wp, &
      & 5.53985700081892E-03_wp,-8.78793958800221E-02_wp, 3.53152003767956E-02_wp, &
      &-2.43105967018234E-02_wp,-1.51574111619879E-02_wp, 2.85915868985803E-02_wp, &
      &-9.49249765705375E-03_wp, 5.64940155335164E-02_wp, 3.11414586934253E-02_wp, &
      &-4.42906788242026E-02_wp,-1.45144645636867E-02_wp,-5.85249818311784E-03_wp, &
      & 9.09568901510486E-01_wp, 2.36486322778065E-02_wp,-3.50001068673827E-04_wp, &
      &-9.07751827589900E-03_wp, 1.23464513582564E-02_wp, 1.66020669497395E-02_wp, &
      &-2.28312252496091E-02_wp, 1.33634035581096E-02_wp, 1.41764609986090E-03_wp, &
      &-5.45979916231744E-02_wp, 5.34400854366388E-03_wp,-4.96044569641485E-06_wp, &
      & 2.40201841492479E-03_wp, 2.48198884031917E-03_wp,-2.72121814524756E-03_wp, &
      & 1.34763066806291E-02_wp, 1.29062067504193E-02_wp,-9.80190748893888E-03_wp, &
      &-5.04907733624908E-03_wp,-2.40085902556881E-03_wp,-1.16119924081096E-03_wp, &
      &-5.81093006974993E-04_wp, 2.45027490200436E-03_wp, 1.39966952488807E-03_wp, &
      & 1.14006433275942E-03_wp,-8.09390428732436E-04_wp,-7.10043290946883E-04_wp, &
      &-2.46172687351467E-04_wp, 8.43040832431947E-04_wp, 2.84420597616094E-03_wp, &
      & 2.74276798196795E-04_wp, 1.51665855913107E-03_wp,-1.20767003906758E-03_wp, &
      & 3.49311376473168E-04_wp,-1.88909058022195E-03_wp,-7.82803317697228E-04_wp, &
      & 1.52421575467949E-03_wp,-5.48713980681319E-04_wp, 1.60503540214095E-03_wp, &
      & 1.24763488181786E-03_wp,-3.26497777863993E-03_wp, 1.55005079704011E-02_wp, &
      &-4.20396123671864E-05_wp, 2.36486322778065E-02_wp, 9.18511119490182E-04_wp, &
      & 2.19078907284031E-04_wp,-4.80302880422083E-04_wp, 3.36883951800672E-04_wp, &
      & 1.96443072180799E-05_wp,-9.42305899810999E-04_wp,-2.89537305206217E-04_wp, &
      &-1.28022878279255E-03_wp,-5.00721188868593E-03_wp, 1.98157365462255E-04_wp, &
      & 6.86048439806152E-05_wp, 1.12897189556701E-04_wp, 8.04397669199020E-05_wp, &
      &-1.60993997929043E-04_wp, 9.24348836836794E-04_wp,-1.64629006596400E-03_wp, &
      &-1.88908493067740E-03_wp,-5.38795740743565E-04_wp, 5.99810932075002E-05_wp, &
      & 1.05424830110665E-05_wp, 4.34016120608695E-05_wp, 1.16049685662807E-04_wp, &
      & 5.66604383383426E-05_wp, 1.14147229303161E-03_wp,-6.27932511020789E-04_wp, &
      & 2.89620397186540E-06_wp,-6.32133131700029E-04_wp, 2.57139294590380E-03_wp, &
      & 3.81213248123848E-04_wp,-2.28477526111468E-03_wp, 4.92005695112074E-04_wp, &
      &-4.48556039119347E-04_wp,-7.42699380839901E-04_wp,-1.84620872265765E-03_wp, &
      & 8.27111347824794E-04_wp, 1.17133429597765E-03_wp, 5.54022846553796E-04_wp, &
      & 1.52655440518585E-04_wp, 1.53320635726685E-03_wp,-3.02957050163202E-03_wp, &
      & 1.20367234816803E-02_wp, 2.32696964109612E-02_wp,-3.50001068673827E-04_wp, &
      & 2.19078907284031E-04_wp, 7.86474390598049E-04_wp, 1.81027823094485E-04_wp, &
      & 4.09324704116455E-04_wp,-3.38483160385626E-04_wp, 7.13341494118537E-04_wp, &
      &-2.49338759770228E-03_wp, 1.68735866975632E-04_wp,-7.41309576023266E-04_wp, &
      &-5.88299190152080E-05_wp, 5.59832212371451E-05_wp,-6.75477375787617E-05_wp, &
      & 1.76396292110500E-05_wp,-1.27813548615544E-04_wp,-3.12515170166504E-04_wp, &
      &-3.56717280611389E-03_wp, 3.00321001380816E-04_wp,-6.77069550327497E-04_wp, &
      & 1.50498405466455E-04_wp, 8.35670166748990E-05_wp,-1.78507090937967E-05_wp, &
      & 5.03854201102447E-06_wp, 9.88365521261211E-05_wp,-3.02272998088171E-04_wp, &
      & 1.21450766756593E-04_wp, 6.06131016453156E-04_wp, 8.27214390433664E-05_wp, &
      & 4.69721320177037E-04_wp,-1.93389141426929E-03_wp,-1.54670546994734E-03_wp, &
      &-1.33814279600227E-03_wp,-1.81086224544769E-04_wp, 4.39431782533191E-06_wp, &
      & 8.40605401675080E-04_wp, 1.17930469715137E-03_wp,-5.70259377186291E-04_wp, &
      & 8.01687607586127E-04_wp,-3.37131360832134E-04_wp, 7.48314601602185E-08_wp, &
      & 1.03160173383662E-03_wp,-1.34556376466504E-02_wp, 1.40699101827079E-02_wp, &
      &-9.07751827589900E-03_wp,-4.80302880422083E-04_wp, 1.81027823094485E-04_wp, &
      & 5.18247582373672E-04_wp, 1.22109199673396E-04_wp, 1.93563736986049E-04_wp, &
      & 6.94512792358934E-04_wp,-7.43041114473191E-04_wp, 1.90087784391754E-03_wp, &
      & 1.41762342890851E-03_wp,-8.72768819206869E-05_wp,-6.27592023471956E-05_wp, &
      &-9.94397095236798E-05_wp, 1.69999823758290E-05_wp, 1.02867249071332E-05_wp, &
      & 6.86263434466547E-04_wp, 1.12323534619302E-03_wp, 1.84122040819053E-03_wp, &
      &-1.05685564533451E-03_wp,-7.41528292047348E-05_wp,-5.64959811539522E-06_wp, &
      &-9.38849407025415E-05_wp,-6.44076461307482E-05_wp, 4.02013059080795E-05_wp, &
      & 7.61032562708509E-04_wp,-6.28894926244619E-04_wp, 2.82293389951878E-06_wp, &
      &-1.20209486448796E-04_wp, 1.34881283254919E-03_wp, 4.82420477979836E-04_wp, &
      &-1.52509696471363E-03_wp, 9.60493988868284E-06_wp,-1.30225149564277E-03_wp, &
      & 3.27030474840371E-04_wp,-8.31609991216725E-04_wp, 5.54949096390952E-04_wp, &
      & 8.00501008899362E-04_wp, 3.77777202539428E-04_wp, 1.19451943984326E-03_wp, &
      & 1.22592371170185E-03_wp,-2.06224126752325E-03_wp,-3.44948217512164E-04_wp, &
      & 1.57206108208992E-02_wp, 1.23464513582564E-02_wp, 3.36883951800672E-04_wp, &
      & 4.09324704116455E-04_wp, 1.22109199673396E-04_wp, 4.58251142523724E-04_wp, &
      & 2.38692760788248E-04_wp,-7.45009963602636E-04_wp,-8.81995629091509E-04_wp, &
      & 4.45232108062044E-04_wp,-3.31802835554287E-03_wp, 9.72796042636198E-05_wp, &
      & 4.59931850098990E-06_wp, 8.59895681224011E-06_wp, 8.83186095660427E-05_wp, &
      &-1.38666929038775E-04_wp, 9.28623086750347E-04_wp, 1.19596404451573E-04_wp, &
      & 5.87522199738515E-06_wp,-2.03335470252928E-03_wp,-4.17552021346422E-05_wp, &
      &-5.09405249892863E-06_wp,-6.69138499163273E-05_wp, 5.12364190514972E-05_wp, &
      & 1.03787374695023E-04_wp,-4.79144450852202E-04_wp,-7.62484937972927E-05_wp, &
      & 2.86985578991552E-04_wp, 7.01975281929019E-04_wp,-1.20147515573047E-03_wp, &
      &-5.32810176688890E-04_wp,-8.34604782021031E-05_wp,-1.49543254126587E-03_wp, &
      &-1.72564923333505E-03_wp, 1.58749721679993E-03_wp, 1.57654392909966E-03_wp, &
      & 3.31242618197280E-04_wp,-5.93648368390723E-04_wp, 1.95523826272871E-04_wp, &
      & 1.48935194145954E-03_wp,-1.03037430243343E-04_wp, 1.29590186401136E-03_wp, &
      &-2.38477695614008E-02_wp, 4.82043997276390E-05_wp, 1.66020669497395E-02_wp, &
      & 1.96443072180799E-05_wp,-3.38483160385626E-04_wp, 1.93563736986049E-04_wp, &
      & 2.38692760788248E-04_wp, 9.55996160124030E-04_wp,-2.25331678846465E-03_wp, &
      & 1.15282575115821E-03_wp, 1.46723806272367E-03_wp,-2.97238116407435E-03_wp, &
      & 1.74422156274703E-04_wp,-9.53909049759476E-05_wp, 4.56054718637186E-05_wp, &
      & 1.19023687658626E-04_wp,-4.22374644343754E-05_wp, 2.15549509066774E-03_wp, &
      & 5.25101958893504E-03_wp, 1.29756344050034E-05_wp,-2.17978866795582E-03_wp, &
      &-2.95874145259878E-04_wp,-1.25680696991265E-04_wp,-1.08710889283027E-04_wp, &
      & 5.08755648891538E-05_wp, 3.70646570214242E-05_wp,-1.30082860861911E-03_wp, &
      &-8.00507020681381E-03_wp,-3.80262208043235E-03_wp, 5.46887348152095E-03_wp, &
      &-2.39915293849376E-02_wp, 1.83968711176338E-02_wp,-9.70631584445533E-03_wp, &
      &-1.19162210963841E-02_wp,-8.62878718054140E-03_wp, 5.59939284194209E-04_wp, &
      &-6.32238839106823E-04_wp,-2.19959469457463E-04_wp,-8.89167773593586E-04_wp, &
      & 6.68905535755028E-04_wp,-9.31069490873247E-04_wp, 6.85864626230428E-04_wp, &
      &-1.73578590097926E-03_wp, 5.97412493903298E-03_wp, 1.29832452985145E-03_wp, &
      &-2.28312252496091E-02_wp,-9.42305899810999E-04_wp, 7.13341494118537E-04_wp, &
      & 6.94512792358934E-04_wp,-7.45009963602636E-04_wp,-2.25331678846465E-03_wp, &
      & 9.80974511856630E-01_wp, 6.77159511266487E-02_wp, 3.31547266836877E-02_wp, &
      &-4.62531679479594E-02_wp, 3.56741815719038E-03_wp,-2.66711332001602E-03_wp, &
      & 1.45224189985198E-03_wp, 1.78354903062118E-03_wp, 1.32108871221236E-03_wp, &
      &-3.02557120515213E-03_wp, 1.59240225965453E-02_wp, 1.87604065635365E-02_wp, &
      & 4.69837569301331E-03_wp,-1.69113411385159E-03_wp,-6.61974620886249E-04_wp, &
      &-1.25193134896850E-03_wp,-1.80059796256980E-03_wp,-3.51490801006137E-04_wp, &
      &-3.51696354973918E-02_wp,-1.18379012605945E-02_wp,-1.56986199671817E-02_wp, &
      & 2.07465098916424E-02_wp, 7.91814684776186E-02_wp,-6.46585082511956E-02_wp, &
      & 5.89106335179422E-04_wp, 7.77031024200781E-02_wp,-3.71725125852140E-02_wp, &
      &-4.33600687024431E-02_wp, 6.52028457208073E-02_wp,-3.87427809377845E-03_wp, &
      & 4.01388520677411E-02_wp,-5.70282419272425E-03_wp, 7.87169218444709E-03_wp, &
      &-3.93198699151506E-02_wp,-4.78755981358135E-03_wp, 3.58652467089654E-03_wp, &
      &-1.01189621165404E-02_wp, 1.33634035581096E-02_wp,-2.89537305206217E-04_wp, &
      &-2.49338759770228E-03_wp,-7.43041114473191E-04_wp,-8.81995629091509E-04_wp, &
      & 1.15282575115821E-03_wp, 6.77159511266487E-02_wp, 8.93041610503400E-01_wp, &
      &-4.10359825689074E-03_wp, 1.29601523129676E-02_wp, 1.67095113530985E-02_wp, &
      &-1.36647926268533E-02_wp, 1.27251154417880E-02_wp,-3.25733892434158E-05_wp, &
      & 2.26849602998059E-02_wp, 9.66633023541495E-03_wp,-5.19297868799894E-03_wp, &
      &-3.84633749324139E-02_wp,-1.24613040121662E-02_wp, 1.38676485413549E-03_wp, &
      &-6.88153886369830E-04_wp, 2.13804102298866E-03_wp, 3.54546470887354E-03_wp, &
      &-6.75694770549573E-06_wp,-1.39407590739072E-02_wp,-1.74855485990918E-02_wp, &
      & 9.95289435444545E-03_wp, 1.12870324819581E-02_wp, 6.88501216429309E-02_wp, &
      & 3.76360397992007E-02_wp, 8.54634077370166E-02_wp,-3.00536616447836E-02_wp, &
      & 3.53696185434452E-02_wp,-1.57001763492994E-02_wp, 1.09811119079166E-04_wp, &
      & 4.92359767360659E-02_wp, 6.11538319289382E-03_wp,-3.25394854692952E-02_wp, &
      & 3.89347577592347E-03_wp,-4.33959140311928E-02_wp,-4.80964039210276E-03_wp, &
      &-1.13563354350246E-02_wp, 1.64501682042580E-02_wp, 1.41764609986090E-03_wp, &
      &-1.28022878279255E-03_wp, 1.68735866975632E-04_wp, 1.90087784391754E-03_wp, &
      & 4.45232108062044E-04_wp, 1.46723806272367E-03_wp, 3.31547266836877E-02_wp, &
      &-4.10359825689074E-03_wp, 9.08117769603523E-01_wp, 9.12873283600010E-03_wp, &
      &-3.41667105883292E-04_wp,-2.27240427085510E-02_wp,-1.60656013699116E-02_wp, &
      & 1.70173883505167E-02_wp, 3.01970480128171E-04_wp, 1.87221098337107E-02_wp, &
      &-2.32165734653989E-02_wp,-2.18459102876561E-02_wp,-3.26536302873989E-02_wp, &
      & 1.96920444748044E-03_wp, 5.26310641005301E-04_wp,-1.60834820235921E-05_wp, &
      & 3.65950247565864E-03_wp, 4.03502323498709E-03_wp, 2.24022500427104E-02_wp, &
      & 2.22060923985063E-02_wp, 1.08109379474327E-02_wp, 2.62360032992966E-03_wp, &
      &-1.27801136787665E-02_wp, 6.76726646968229E-02_wp,-6.23156571113277E-03_wp, &
      &-2.00732274112082E-03_wp,-9.55348684789663E-02_wp,-3.92643295028350E-02_wp, &
      &-2.52064987261363E-02_wp, 3.76026077550443E-04_wp,-3.43827040288937E-02_wp, &
      &-4.86038808135031E-03_wp,-4.72243223724468E-02_wp, 3.84544626969289E-02_wp, &
      &-2.26550292609627E-02_wp,-9.59168996535191E-03_wp,-1.13646125614409E-02_wp, &
      &-5.45979916231744E-02_wp,-5.00721188868593E-03_wp,-7.41309576023266E-04_wp, &
      & 1.41762342890851E-03_wp,-3.31802835554287E-03_wp,-2.97238116407435E-03_wp, &
      &-4.62531679479594E-02_wp, 1.29601523129676E-02_wp, 9.12873283600010E-03_wp, &
      & 9.09645248969603E-01_wp,-2.27265262407807E-02_wp,-5.32299199303875E-04_wp, &
      &-1.00495621815055E-02_wp,-1.37190880213870E-02_wp, 1.76759686356936E-02_wp, &
      & 1.35404942972656E-02_wp, 8.14988231623426E-03_wp,-1.12416390562618E-02_wp, &
      &-3.84377541807347E-03_wp,-1.75012611402529E-03_wp,-9.81189201589943E-04_wp, &
      &-3.98266659773776E-04_wp, 2.71148530315031E-03_wp, 1.45446242844221E-03_wp, &
      &-1.30001097214204E-03_wp,-8.41884382085804E-04_wp,-6.58426823850392E-04_wp, &
      & 3.30170019630296E-04_wp, 1.67318715111030E-03_wp,-2.56333530435437E-03_wp, &
      & 2.39237960359928E-04_wp, 1.10329669950774E-03_wp, 1.04090454647008E-03_wp, &
      & 4.03441445098017E-04_wp, 1.91383514382521E-03_wp,-1.50246094084492E-04_wp, &
      & 1.81283369635852E-03_wp,-6.16005794354335E-05_wp, 1.35155308784679E-03_wp, &
      &-1.77253159605633E-03_wp, 1.03244385712790E-03_wp,-6.71551375241214E-04_wp, &
      &-6.58322439410420E-04_wp, 5.34400854366388E-03_wp, 1.98157365462255E-04_wp, &
      &-5.88299190152080E-05_wp,-8.72768819206869E-05_wp, 9.72796042636198E-05_wp, &
      & 1.74422156274703E-04_wp, 3.56741815719038E-03_wp, 1.67095113530985E-02_wp, &
      &-3.41667105883292E-04_wp,-2.27265262407807E-02_wp, 9.20468008197577E-04_wp, &
      &-2.49874089186843E-04_wp, 5.11060765003503E-04_wp, 3.55781691261206E-04_wp, &
      &-2.00686966101574E-05_wp,-1.11069365429552E-03_wp,-1.77595705940930E-03_wp, &
      &-9.69766591163740E-04_wp, 1.79487344769609E-03_wp, 1.18227059103663E-04_wp, &
      & 2.23628169815625E-05_wp, 9.59540695118637E-05_wp,-6.27125251357207E-06_wp, &
      &-7.33810769916223E-05_wp, 9.03676498870888E-04_wp, 6.54018323036504E-04_wp, &
      &-2.72580344706362E-05_wp,-6.85792113938563E-04_wp,-2.26525736225413E-03_wp, &
      & 3.57609160919946E-04_wp,-1.86157192226547E-03_wp,-1.68670894927570E-04_wp, &
      &-3.59221400041849E-04_wp, 1.24040844923198E-03_wp,-1.05856014409741E-03_wp, &
      &-1.27777504135482E-03_wp,-6.03244529432564E-04_wp, 1.02574808923692E-03_wp, &
      &-2.00017008301805E-04_wp, 1.76953567734680E-03_wp, 8.28645456652537E-04_wp, &
      & 2.30409327479312E-03_wp,-9.95021813203161E-05_wp,-4.96044569641485E-06_wp, &
      & 6.86048439806152E-05_wp, 5.59832212371451E-05_wp,-6.27592023471956E-05_wp, &
      & 4.59931850098990E-06_wp,-9.53909049759476E-05_wp,-2.66711332001602E-03_wp, &
      &-1.36647926268533E-02_wp,-2.27240427085510E-02_wp,-5.32299199303875E-04_wp, &
      &-2.49874089186843E-04_wp, 8.05224008513408E-04_wp, 2.12701081422194E-04_wp, &
      &-4.33732331900178E-04_wp,-3.76373327444631E-04_wp,-1.92622491079814E-03_wp, &
      & 3.14448950646059E-03_wp, 3.64600111382421E-03_wp, 2.20063557970657E-03_wp, &
      &-1.22466148179089E-04_wp,-1.14817581188101E-05_wp,-6.53524280632195E-05_wp, &
      &-2.32450086452929E-04_wp,-1.56880573246909E-04_wp,-5.83171540616180E-04_wp, &
      &-9.79735922629530E-05_wp,-5.84203711930900E-04_wp, 6.13175091245107E-05_wp, &
      & 1.03312589059841E-04_wp,-2.02312098196746E-03_wp,-1.22914115729581E-03_wp, &
      & 1.46237940335649E-03_wp,-3.57354903215663E-04_wp, 2.73735765525016E-04_wp, &
      & 1.27412333131700E-03_wp,-1.02227177755345E-03_wp, 9.60330421289763E-04_wp, &
      & 4.82756606241590E-04_wp, 5.23748289028198E-04_wp,-2.92690161508324E-04_wp, &
      & 8.36183698869665E-04_wp, 4.33495391693635E-04_wp,-1.69154717257233E-03_wp, &
      & 2.40201841492479E-03_wp, 1.12897189556701E-04_wp,-6.75477375787617E-05_wp, &
      &-9.94397095236798E-05_wp, 8.59895681224011E-06_wp, 4.56054718637186E-05_wp, &
      & 1.45224189985198E-03_wp, 1.27251154417880E-02_wp,-1.60656013699116E-02_wp, &
      &-1.00495621815055E-02_wp, 5.11060765003503E-04_wp, 2.12701081422194E-04_wp, &
      & 5.85429402823081E-04_wp,-1.44929109080006E-04_wp, 1.19686671170196E-04_wp, &
      &-1.26135270204623E-03_wp, 3.43298532897914E-04_wp, 3.83817970991369E-06_wp, &
      & 2.16308617823652E-03_wp, 1.39871480504802E-05_wp,-1.14181893240379E-05_wp, &
      & 6.32963494893167E-05_wp,-6.58156831461232E-05_wp,-1.36918692981879E-04_wp, &
      &-5.86376985234874E-04_wp,-7.38014920537504E-04_wp, 1.67615675300633E-05_wp, &
      & 2.02457181360983E-04_wp, 1.25511859423494E-03_wp,-3.47311396600001E-04_wp, &
      & 1.62380339070335E-03_wp,-6.46986175422616E-04_wp, 1.79429741173201E-03_wp, &
      & 2.82584056419365E-04_wp, 3.48350771197910E-04_wp, 9.55873747394603E-04_wp, &
      & 6.66961474578803E-04_wp,-5.53868614428592E-04_wp, 8.65609593211044E-04_wp, &
      &-1.42187034413839E-03_wp, 3.56016852022877E-04_wp,-1.03517926424556E-03_wp, &
      & 8.98861322776528E-04_wp, 2.48198884031917E-03_wp, 8.04397669199020E-05_wp, &
      & 1.76396292110500E-05_wp, 1.69999823758290E-05_wp, 8.83186095660427E-05_wp, &
      & 1.19023687658626E-04_wp, 1.78354903062118E-03_wp,-3.25733892434158E-05_wp, &
      & 1.70173883505167E-02_wp,-1.37190880213870E-02_wp, 3.55781691261206E-04_wp, &
      &-4.33732331900178E-04_wp,-1.44929109080006E-04_wp, 5.39623595018980E-04_wp, &
      &-2.59213800285308E-04_wp,-4.02492245771770E-05_wp,-2.31816357724684E-03_wp, &
      &-6.97563672918741E-04_wp,-1.10536488762744E-03_wp, 1.02416470270342E-04_wp, &
      & 4.30064827780300E-05_wp, 1.06480846519699E-05_wp, 5.28437067393580E-05_wp, &
      & 9.07005439053409E-05_wp,-5.20172683811180E-04_wp, 1.67408952638742E-04_wp, &
      &-1.95672888720847E-04_wp, 6.65644753860528E-04_wp, 1.23392237280277E-03_wp, &
      &-5.79280405069388E-04_wp,-2.57899109827770E-04_wp, 1.68257298559226E-03_wp, &
      &-2.19262341202439E-03_wp,-2.02260087963467E-03_wp, 1.31355340280529E-03_wp, &
      &-4.76745213909973E-05_wp, 1.38916670200151E-04_wp,-4.61704646589813E-04_wp, &
      &-8.12160921508423E-04_wp,-3.56619680821674E-04_wp,-2.26721397739007E-03_wp, &
      &-9.32392910000516E-04_wp,-2.17156562704799E-03_wp,-2.72121814524756E-03_wp, &
      &-1.60993997929043E-04_wp,-1.27813548615544E-04_wp, 1.02867249071332E-05_wp, &
      &-1.38666929038775E-04_wp,-4.22374644343754E-05_wp, 1.32108871221236E-03_wp, &
      & 2.26849602998059E-02_wp, 3.01970480128171E-04_wp, 1.76759686356936E-02_wp, &
      &-2.00686966101574E-05_wp,-3.76373327444631E-04_wp, 1.19686671170196E-04_wp, &
      &-2.59213800285308E-04_wp, 9.32898777817163E-04_wp, 1.33666846700357E-03_wp, &
      &-3.22194483116882E-04_wp,-4.38986218965597E-03_wp,-8.07406407670481E-04_wp, &
      & 5.55356868898791E-06_wp,-6.21201449864920E-05_wp, 9.70425331312843E-05_wp, &
      & 2.32661759561976E-04_wp, 3.40426800813873E-05_wp,-1.20825923162132E-03_wp, &
      & 2.55413860729610E-03_wp,-3.63428057490289E-03_wp,-9.08181098067992E-03_wp, &
      &-1.35620731019103E-02_wp,-5.21056540703202E-03_wp,-9.88014923090812E-03_wp, &
      & 2.13600011676250E-02_wp, 2.19832369216965E-02_wp,-8.83664108598077E-04_wp, &
      & 5.25854327632037E-05_wp, 5.49263377585196E-04_wp,-8.40001714014691E-04_wp, &
      &-4.42109441645743E-04_wp, 1.03952145445938E-03_wp, 2.40006331976448E-04_wp, &
      &-4.67936049592108E-04_wp,-1.75636665246333E-02_wp, 1.74725131414713E-03_wp, &
      & 1.34763066806291E-02_wp, 9.24348836836794E-04_wp,-3.12515170166504E-04_wp, &
      & 6.86263434466547E-04_wp, 9.28623086750347E-04_wp, 2.15549509066774E-03_wp, &
      &-3.02557120515213E-03_wp, 9.66633023541495E-03_wp, 1.87221098337107E-02_wp, &
      & 1.35404942972656E-02_wp,-1.11069365429552E-03_wp,-1.92622491079814E-03_wp, &
      &-1.26135270204623E-03_wp,-4.02492245771770E-05_wp, 1.33666846700357E-03_wp, &
      & 9.80364656316394E-01_wp,-2.18477518470214E-02_wp, 3.31124366005084E-02_wp, &
      & 7.98537612496862E-02_wp, 1.97816816411409E-03_wp, 7.95287884239187E-04_wp, &
      & 1.47614735664940E-03_wp,-3.09112108744084E-03_wp,-3.26315136987645E-03_wp, &
      & 1.01661826271289E-02_wp, 1.44077180987701E-02_wp, 5.13544987508727E-03_wp, &
      & 1.33513181820185E-02_wp, 7.56618991833905E-02_wp, 3.66511396790250E-02_wp, &
      &-6.25144814937098E-03_wp, 3.25281161052200E-02_wp, 7.20780368313875E-02_wp, &
      & 3.13114544459522E-02_wp, 2.45048495855122E-02_wp,-6.37450174955805E-03_wp, &
      &-1.91461291146198E-02_wp, 1.09535528037961E-03_wp, 3.42404939671757E-02_wp, &
      & 4.74589429392366E-02_wp,-2.15316812457990E-02_wp,-4.80394095409376E-02_wp, &
      &-1.23060391666262E-02_wp, 1.29062067504193E-02_wp,-1.64629006596400E-03_wp, &
      &-3.56717280611389E-03_wp, 1.12323534619302E-03_wp, 1.19596404451573E-04_wp, &
      & 5.25101958893504E-03_wp, 1.59240225965453E-02_wp,-5.19297868799894E-03_wp, &
      &-2.32165734653989E-02_wp, 8.14988231623426E-03_wp,-1.77595705940930E-03_wp, &
      & 3.14448950646059E-03_wp, 3.43298532897914E-04_wp,-2.31816357724684E-03_wp, &
      &-3.22194483116882E-04_wp,-2.18477518470214E-02_wp, 9.14810407256861E-01_wp, &
      & 7.48033986208798E-03_wp, 5.09953559406378E-03_wp,-2.71827543681093E-02_wp, &
      &-1.35396965598036E-02_wp,-5.28749471868861E-03_wp,-4.55353025074698E-04_wp, &
      &-8.97158346345382E-03_wp,-1.38165545028930E-02_wp, 4.81975145073823E-03_wp, &
      & 1.02100793239385E-02_wp,-1.99268786228131E-02_wp, 3.15084854540213E-02_wp, &
      &-1.61959501492248E-02_wp, 8.60612099656798E-02_wp, 4.67381503378004E-02_wp, &
      &-7.10320337351810E-02_wp, 3.67817089890598E-02_wp,-2.56655009577307E-03_wp, &
      &-1.46067658017428E-02_wp, 5.62002262438386E-03_wp, 5.78076943435148E-02_wp, &
      &-2.44882129383753E-03_wp,-2.86025601155137E-02_wp,-5.04588128389885E-03_wp, &
      &-4.04081526462394E-03_wp, 1.53478173145738E-02_wp,-9.80190748893888E-03_wp, &
      &-1.88908493067740E-03_wp, 3.00321001380816E-04_wp, 1.84122040819053E-03_wp, &
      & 5.87522199738515E-06_wp, 1.29756344050034E-05_wp, 1.87604065635365E-02_wp, &
      &-3.84633749324139E-02_wp,-2.18459102876561E-02_wp,-1.12416390562618E-02_wp, &
      &-9.69766591163740E-04_wp, 3.64600111382421E-03_wp, 3.83817970991369E-06_wp, &
      &-6.97563672918741E-04_wp,-4.38986218965597E-03_wp, 3.31124366005084E-02_wp, &
      & 7.48033986208798E-03_wp, 9.07045088096441E-01_wp,-6.65538482288344E-03_wp, &
      &-4.69539950187602E-04_wp, 8.56203889791274E-03_wp,-1.56700021141420E-02_wp, &
      &-2.67566253419244E-02_wp,-5.14238323186605E-05_wp,-4.08742499801316E-02_wp, &
      & 1.15402494798051E-02_wp,-1.81555166175954E-02_wp,-2.34217202320817E-02_wp, &
      & 6.45053437127255E-02_wp, 4.22610487435240E-02_wp, 3.14099372304354E-03_wp, &
      &-1.03858523700384E-01_wp,-4.72446672147925E-02_wp, 5.10642880072331E-02_wp, &
      & 4.06754362288397E-02_wp,-4.95983018633171E-03_wp, 4.94693462764905E-02_wp, &
      &-1.68594875769831E-03_wp,-6.15913706401558E-02_wp,-2.74869037414883E-02_wp, &
      & 2.25098483451193E-03_wp, 2.91067062790983E-02_wp,-5.75999068005323E-03_wp, &
      &-5.04907733624908E-03_wp,-5.38795740743565E-04_wp,-6.77069550327497E-04_wp, &
      &-1.05685564533451E-03_wp,-2.03335470252928E-03_wp,-2.17978866795582E-03_wp, &
      & 4.69837569301331E-03_wp,-1.24613040121662E-02_wp,-3.26536302873989E-02_wp, &
      &-3.84377541807347E-03_wp, 1.79487344769609E-03_wp, 2.20063557970657E-03_wp, &
      & 2.16308617823652E-03_wp,-1.10536488762744E-03_wp,-8.07406407670481E-04_wp, &
      & 7.98537612496862E-02_wp, 5.09953559406378E-03_wp,-6.65538482288344E-03_wp, &
      & 8.85720560609372E-01_wp, 8.47405040928117E-03_wp, 4.32609508448040E-05_wp, &
      & 1.51332117219224E-02_wp,-1.32002681370741E-02_wp,-2.63946101633402E-02_wp, &
      &-7.18067691148890E-04_wp,-3.91472860045293E-04_wp,-3.91264029795281E-04_wp, &
      &-6.69813653944916E-04_wp,-9.38359904068336E-04_wp,-7.23150166635905E-04_wp, &
      & 3.78532516914578E-04_wp,-1.50740974288056E-03_wp,-2.04798561893996E-03_wp, &
      &-6.93371177426215E-04_wp,-4.24473085272663E-04_wp, 1.86494551430705E-04_wp, &
      & 1.30712345342734E-03_wp,-1.50326172865105E-04_wp,-1.62689399094736E-03_wp, &
      &-1.81297319011831E-03_wp, 2.05812072051282E-03_wp, 5.41082241160368E-03_wp, &
      & 6.00827861785183E-04_wp,-2.40085902556881E-03_wp, 5.99810932075002E-05_wp, &
      & 1.50498405466455E-04_wp,-7.41528292047348E-05_wp,-4.17552021346422E-05_wp, &
      &-2.95874145259878E-04_wp,-1.69113411385159E-03_wp, 1.38676485413549E-03_wp, &
      & 1.96920444748044E-03_wp,-1.75012611402529E-03_wp, 1.18227059103663E-04_wp, &
      &-1.22466148179089E-04_wp, 1.39871480504802E-05_wp, 1.02416470270342E-04_wp, &
      & 5.55356868898791E-06_wp, 1.97816816411409E-03_wp,-2.71827543681093E-02_wp, &
      &-4.69539950187602E-04_wp, 8.47405040928117E-03_wp, 9.24530237513042E-04_wp, &
      & 4.12508019867131E-04_wp, 3.21310649550534E-04_wp,-1.04518495185366E-04_wp, &
      & 1.07469239917818E-05_wp,-2.63641797653816E-04_wp,-1.85447400076972E-04_wp, &
      & 9.88365504774939E-06_wp,-4.38712775661301E-04_wp,-6.34043613156057E-04_wp, &
      &-6.96300015777383E-04_wp, 9.38872653962262E-04_wp, 3.82651994861127E-05_wp, &
      &-1.49699282601155E-03_wp,-1.52128659470300E-04_wp,-4.73872134127600E-04_wp, &
      &-2.30162533491128E-05_wp, 4.10116747608054E-04_wp, 5.22664015461809E-04_wp, &
      &-5.53268037298817E-04_wp,-9.90759348859535E-04_wp, 5.05766248283295E-04_wp, &
      & 2.05070858833928E-03_wp, 6.76325363903052E-04_wp,-1.16119924081096E-03_wp, &
      & 1.05424830110665E-05_wp, 8.35670166748990E-05_wp,-5.64959811539522E-06_wp, &
      &-5.09405249892863E-06_wp,-1.25680696991265E-04_wp,-6.61974620886249E-04_wp, &
      &-6.88153886369830E-04_wp, 5.26310641005301E-04_wp,-9.81189201589943E-04_wp, &
      & 2.23628169815625E-05_wp,-1.14817581188101E-05_wp,-1.14181893240379E-05_wp, &
      & 4.30064827780300E-05_wp,-6.21201449864920E-05_wp, 7.95287884239187E-04_wp, &
      &-1.35396965598036E-02_wp, 8.56203889791274E-03_wp, 4.32609508448040E-05_wp, &
      & 4.12508019867131E-04_wp, 2.88003393911521E-04_wp,-6.50182407583703E-05_wp, &
      &-2.53634831707529E-04_wp, 1.27233371965883E-04_wp,-5.93966433083350E-04_wp, &
      & 2.69372961745411E-05_wp,-5.75792938733305E-04_wp,-1.20006019567996E-04_wp, &
      & 3.14862311605364E-04_wp, 6.92000783727210E-04_wp,-1.19990637124721E-03_wp, &
      &-2.40313487407472E-03_wp, 1.86593826628628E-04_wp,-6.60971351199133E-05_wp, &
      & 6.75497357337969E-04_wp, 1.12528083885423E-04_wp, 9.58776573157139E-04_wp, &
      &-1.09131629766252E-03_wp,-1.20341681703315E-03_wp,-3.97806973349569E-04_wp, &
      & 7.91454294590884E-04_wp, 2.29339909979128E-03_wp,-1.66014282991117E-03_wp, &
      &-5.81093006974993E-04_wp, 4.34016120608695E-05_wp,-1.78507090937967E-05_wp, &
      &-9.38849407025415E-05_wp,-6.69138499163273E-05_wp,-1.08710889283027E-04_wp, &
      &-1.25193134896850E-03_wp, 2.13804102298866E-03_wp,-1.60834820235921E-05_wp, &
      &-3.98266659773776E-04_wp, 9.59540695118637E-05_wp,-6.53524280632195E-05_wp, &
      & 6.32963494893167E-05_wp, 1.06480846519699E-05_wp, 9.70425331312843E-05_wp, &
      & 1.47614735664940E-03_wp,-5.28749471868861E-03_wp,-1.56700021141420E-02_wp, &
      & 1.51332117219224E-02_wp, 3.21310649550534E-04_wp,-6.50182407583703E-05_wp, &
      & 5.65537994689567E-04_wp, 2.39816673197796E-04_wp,-4.01145908943532E-04_wp, &
      & 1.02979125404260E-03_wp,-3.83960703078049E-04_wp,-4.71112289724211E-05_wp, &
      & 1.02053886193382E-03_wp,-1.37787803489053E-03_wp, 2.14568368478059E-04_wp, &
      &-2.29785457537120E-03_wp, 3.71631530998101E-04_wp, 2.32423304833064E-03_wp, &
      &-2.04580162271769E-03_wp,-5.61485634191755E-04_wp, 5.83275218807406E-04_wp, &
      &-7.52025031398366E-04_wp,-1.82188558275370E-03_wp, 1.03007344334789E-03_wp, &
      & 1.22820555419500E-03_wp, 7.16816433450587E-04_wp, 4.07757821937252E-04_wp, &
      &-2.79625792410964E-04_wp, 2.45027490200436E-03_wp, 1.16049685662807E-04_wp, &
      & 5.03854201102447E-06_wp,-6.44076461307482E-05_wp, 5.12364190514972E-05_wp, &
      & 5.08755648891538E-05_wp,-1.80059796256980E-03_wp, 3.54546470887354E-03_wp, &
      & 3.65950247565864E-03_wp, 2.71148530315031E-03_wp,-6.27125251357207E-06_wp, &
      &-2.32450086452929E-04_wp,-6.58156831461232E-05_wp, 5.28437067393580E-05_wp, &
      & 2.32661759561976E-04_wp,-3.09112108744084E-03_wp,-4.55353025074698E-04_wp, &
      &-2.67566253419244E-02_wp,-1.32002681370741E-02_wp,-1.04518495185366E-04_wp, &
      &-2.53634831707529E-04_wp, 2.39816673197796E-04_wp, 1.01972520181169E-03_wp, &
      & 4.18800952996832E-04_wp, 1.20885516832022E-03_wp,-5.52056921452457E-04_wp, &
      & 5.46257757174961E-04_wp, 5.94007597988431E-04_wp,-2.24225399100506E-03_wp, &
      &-1.19565402301583E-03_wp, 4.74724238068227E-05_wp, 2.60157115829106E-03_wp, &
      & 4.22160666294196E-04_wp,-1.86857557022708E-03_wp,-1.60968883271744E-03_wp, &
      & 4.19280315359715E-04_wp,-1.22301380113740E-03_wp,-3.25003079162797E-05_wp, &
      & 1.53200056222922E-03_wp, 4.56941288207854E-04_wp, 1.18902455312498E-03_wp, &
      &-4.27488088447255E-04_wp, 2.06014446932532E-03_wp, 1.39966952488807E-03_wp, &
      & 5.66604383383426E-05_wp, 9.88365521261211E-05_wp, 4.02013059080795E-05_wp, &
      & 1.03787374695023E-04_wp, 3.70646570214242E-05_wp,-3.51490801006137E-04_wp, &
      &-6.75694770549573E-06_wp, 4.03502323498709E-03_wp, 1.45446242844221E-03_wp, &
      &-7.33810769916223E-05_wp,-1.56880573246909E-04_wp,-1.36918692981879E-04_wp, &
      & 9.07005439053409E-05_wp, 3.40426800813873E-05_wp,-3.26315136987645E-03_wp, &
      &-8.97158346345382E-03_wp,-5.14238323186605E-05_wp,-2.63946101633402E-02_wp, &
      & 1.07469239917818E-05_wp, 1.27233371965883E-04_wp,-4.01145908943532E-04_wp, &
      & 4.18800952996832E-04_wp, 8.91699965008798E-04_wp], shape(density))

   call get_structure(mol, "f-block", "CeCl3")
   call test_num_op_grad(error, mol, density, make_exchange_gxtb, thr_in=thr1)

end subroutine test_op_g_fock_cecl3


!> Test analytic vs numerical gradient of exchange energy for pyrazine (PBC, op grad)
subroutine test_op_g_fock_pyrazine(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(56,56,1) = reshape([ &
   &  4.996323064703260E-01_wp, -6.038045829436379E-02_wp, -3.865880949483817E-03_wp, &
   & -3.469278033590685E-03_wp,  3.674180808123230E-02_wp,  4.244475807793323E-04_wp, &
   & -4.719659729581716E-03_wp, -2.043891119358536E-03_wp,  3.215652341046649E-01_wp, &
   &  3.949845594623018E-01_wp,  1.877266064541833E-01_wp,  2.805496228740723E-01_wp, &
   & -1.226809308310684E-01_wp, -6.235703051699102E-02_wp, -2.468880734654038E-02_wp, &
   & -1.208993709601304E-01_wp,  5.498706220187589E-03_wp, -2.682407928818319E-03_wp, &
   &  2.693644321560950E-04_wp,  1.303015842936295E-03_wp, -4.287832154916927E-03_wp, &
   &  2.044003485932207E-04_wp,  8.751385748047206E-03_wp,  9.014283332274178E-03_wp, &
   &  1.315713396064210E-02_wp,  1.421731639133684E-02_wp,  1.982936508356516E-03_wp, &
   & -1.638356792386353E-02_wp,  4.415729811646429E-02_wp,  1.087878338388485E-01_wp, &
   &  6.085361576972199E-02_wp,  4.337274631722113E-02_wp,  7.676830732914519E-03_wp, &
   &  6.010314084803126E-03_wp,  6.699758088043848E-03_wp, -1.649545855874165E-03_wp, &
   &  3.086428269557153E-04_wp,  5.982866467384993E-03_wp, -3.083173614074472E-03_wp, &
   &  6.032315381140116E-03_wp, -2.018616341381384E-01_wp, -1.348005264422091E-01_wp, &
   & -5.574264986967414E-02_wp, -1.305267507991367E-01_wp,  1.041002211525260E-02_wp, &
   &  7.413909982224848E-02_wp,  2.866334386143702E-02_wp,  1.049620885912693E-01_wp, &
   &  3.468787035654383E-03_wp, -1.556266923650076E-03_wp, -4.697567182304204E-03_wp, &
   &  5.252068851896829E-03_wp, -1.851743395625369E-03_wp,  8.190947240366564E-03_wp, &
   & -2.024491870468237E-03_wp, -3.771975464662083E-03_wp, -6.038045829436379E-02_wp, &
   &  4.996318510375576E-01_wp, -3.466346084102220E-03_wp, -3.866530058951323E-03_wp, &
   &  4.244353408670472E-04_wp,  3.674180706459462E-02_wp, -2.043925312956940E-03_wp, &
   & -4.717513073777494E-03_wp, -1.226814369553450E-01_wp, -6.235649514184920E-02_wp, &
   & -2.468847162697405E-02_wp,  1.208993706003673E-01_wp,  3.215654129343700E-01_wp, &
   &  3.949847425377017E-01_wp,  1.877266874435544E-01_wp, -2.805489103522327E-01_wp, &
   & -4.287430237917707E-03_wp,  2.041263539903958E-04_wp,  8.751701679502094E-03_wp, &
   & -9.013011234728530E-03_wp,  5.499073985967830E-03_wp, -2.682798239092303E-03_wp, &
   &  2.695262479973700E-04_wp, -1.304167939683527E-03_wp,  4.415770787953931E-02_wp, &
   &  1.087881602220137E-01_wp,  6.085384709963981E-02_wp, -4.337273634581573E-02_wp, &
   &  1.315687805394647E-02_wp,  1.421742136850592E-02_wp,  1.982987914676789E-03_wp, &
   &  1.638378434837623E-02_wp,  3.095258120151857E-04_wp,  5.983459493418869E-03_wp, &
   & -3.083347611846924E-03_wp, -6.032919064557218E-03_wp,  7.675707547637729E-03_wp, &
   &  6.008621613144670E-03_wp,  6.700480961870181E-03_wp,  1.649505164758265E-03_wp, &
   &  1.041055766967904E-02_wp,  7.413920473761294E-02_wp,  2.866332331412033E-02_wp, &
   & -1.049613757799108E-01_wp, -2.018623595147256E-01_wp, -1.348005538473235E-01_wp, &
   & -5.574271573258759E-02_wp,  1.305276789102963E-01_wp, -1.852728792443303E-03_wp, &
   &  8.191399118302628E-03_wp, -2.024632422877527E-03_wp,  3.770544286479827E-03_wp, &
   &  3.468244026829145E-03_wp, -1.556663367068126E-03_wp, -4.697511095964893E-03_wp, &
   & -5.251923220272459E-03_wp, -3.865880949483817E-03_wp, -3.466346084102220E-03_wp, &
   &  4.974125368885180E-01_wp, -5.965944466231984E-02_wp, -4.745397789703586E-03_wp, &
   & -2.084498589645663E-03_wp,  3.614820674096364E-02_wp, -4.438882658831626E-04_wp, &
   &  5.478947603167341E-03_wp,  2.662859497750907E-03_wp,  3.560082383852504E-04_wp, &
   &  1.268974885733104E-03_wp, -4.723925904439215E-03_wp, -3.780206736377426E-04_wp, &
   &  8.752624480373896E-03_wp,  9.316120032689994E-03_wp,  3.222353137213799E-01_wp, &
   & -3.955576278448117E-01_wp,  1.878861062451542E-01_wp,  2.787480693862803E-01_wp, &
   & -1.220082606284287E-01_wp,  6.229036665474397E-02_wp, -2.472051104696680E-02_wp, &
   & -1.207796158775903E-01_wp,  7.535473075301793E-03_wp, -6.085886474843062E-03_wp, &
   &  6.870536941104054E-03_wp, -1.590870018319618E-03_wp,  3.428998781745260E-04_wp, &
   & -6.004927066623619E-03_wp, -3.043269941348812E-03_wp,  5.958736679177928E-03_wp, &
   &  1.364727763523009E-02_wp, -1.431013175074783E-02_wp,  2.075142813212527E-03_wp, &
   & -1.604065601170488E-02_wp,  4.499654897229581E-02_wp, -1.089461510977031E-01_wp, &
   &  6.121347798932544E-02_wp,  4.338681583913085E-02_wp,  3.475680838830167E-03_wp, &
   &  1.592836372275454E-03_wp, -4.736727799549166E-03_wp,  5.294951166937761E-03_wp, &
   & -1.756646422878684E-03_wp, -8.182581808537367E-03_wp, -1.750276924826151E-03_wp, &
   & -3.250093390684439E-03_wp, -2.025727736049031E-01_wp,  1.363116666354803E-01_wp, &
   & -5.649549802151434E-02_wp, -1.294885328690167E-01_wp,  9.232026474174537E-03_wp, &
   & -7.303372313506439E-02_wp,  2.808423687111002E-02_wp,  1.037044247528985E-01_wp, &
   & -3.469278033590682E-03_wp, -3.866530058951323E-03_wp, -5.965944466231984E-02_wp, &
   &  4.974211788155059E-01_wp, -2.084534616294795E-03_wp, -4.747543047387469E-03_wp, &
   & -4.486579619528315E-04_wp,  3.614820761444790E-02_wp, -4.727944164665414E-03_wp, &
   & -3.786005378653440E-04_wp,  8.753238377034695E-03_wp, -9.319573280149776E-03_wp, &
   &  5.479481708318162E-03_wp,  2.662153218083750E-03_wp,  3.568278630241579E-04_wp, &
   & -1.267766929892562E-03_wp, -1.220031725743797E-01_wp,  6.229101417531866E-02_wp, &
   & -2.471988022804309E-02_wp,  1.207793448348777E-01_wp,  3.222375528530796E-01_wp, &
   & -3.955613432847502E-01_wp,  1.878888234682954E-01_wp, -2.787558694019131E-01_wp, &
   &  3.436186553753724E-04_wp, -6.005550526819929E-03_wp, -3.042063894358351E-03_wp, &
   & -5.958499739953396E-03_wp,  7.533427596972730E-03_wp, -6.085629091190239E-03_wp, &
   &  6.871785977987184E-03_wp,  1.589292442514503E-03_wp,  4.498941623563778E-02_wp, &
   & -1.089524752117987E-01_wp,  6.121655806829838E-02_wp, -4.338592338238546E-02_wp, &
   &  1.364786448470975E-02_wp, -1.431166066642352E-02_wp,  2.075332588360217E-03_wp, &
   &  1.604230038290055E-02_wp, -1.757649574346473E-03_wp, -8.182952073746532E-03_wp, &
   & -1.752859939197428E-03_wp,  3.253758271617169E-03_wp,  3.473784430573348E-03_wp, &
   &  1.592236996815013E-03_wp, -4.736230410173405E-03_wp, -5.293442289541150E-03_wp, &
   &  9.226336221603473E-03_wp, -7.303191575335419E-02_wp,  2.808355212843569E-02_wp, &
   & -1.037087748751767E-01_wp, -2.025836973893619E-01_wp,  1.363173580525894E-01_wp, &
   & -5.649750140866367E-02_wp,  1.294895097260828E-01_wp,  3.674180808123230E-02_wp, &
   &  4.244353408670472E-04_wp, -4.745397789703586E-03_wp, -2.084534616294799E-03_wp, &
   &  5.012942247855861E-01_wp, -6.081444943148732E-02_wp, -3.898846109092222E-03_wp, &
   & -3.605864269072001E-03_wp,  1.348077095700679E-02_wp, -1.422486198769357E-02_wp, &
   & -2.177761207993463E-03_wp,  1.631366894119213E-02_wp,  4.533861073404212E-02_wp, &
   & -1.095742599472761E-01_wp, -6.138654825322425E-02_wp, -4.390566297627865E-02_wp, &
   &  7.625658910081461E-03_wp, -5.986317947460313E-03_wp, -6.725184651742494E-03_wp, &
   &  1.586802840367607E-03_wp,  3.046710348517129E-04_wp, -5.982650890830495E-03_wp, &
   &  3.142870872497287E-03_wp, -6.012592350652365E-03_wp,  3.214390604969148E-01_wp, &
   & -3.961037152394013E-01_wp, -1.882539383786511E-01_wp, -2.789144572943701E-01_wp, &
   & -1.221179632877208E-01_wp,  6.230722519523442E-02_wp,  2.484861103711767E-02_wp, &
   &  1.200929772037539E-01_wp,  5.575998655428776E-03_wp,  2.765339311704464E-03_wp, &
   & -2.818126004903321E-04_wp, -1.336184314911525E-03_wp, -4.353406478498622E-03_wp, &
   & -2.126723791504598E-04_wp, -8.766400524998663E-03_wp, -9.087342946319639E-03_wp, &
   &  9.637106365162448E-03_wp, -7.369382888401302E-02_wp, -2.841649893839564E-02_wp, &
   & -1.047361603082406E-01_wp, -2.024245396594210E-01_wp,  1.349196958026813E-01_wp, &
   &  5.590846316495357E-02_wp,  1.309283953591511E-01_wp, -1.826507296104563E-03_wp, &
   & -8.232201340172253E-03_wp,  1.976190866293467E-03_wp,  3.755665465210820E-03_wp, &
   &  3.529801573999167E-03_wp,  1.581255206164669E-03_wp,  4.713819655831289E-03_wp, &
   & -5.331683711214077E-03_wp,  4.244475807793288E-04_wp,  3.674180706459462E-02_wp, &
   & -2.084498589645663E-03_wp, -4.747543047387472E-03_wp, -6.081444943148732E-02_wp, &
   &  5.012946780451271E-01_wp, -3.608821299372046E-03_wp, -3.898189592416251E-03_wp, &
   &  4.533819980979545E-02_wp, -1.095739336851609E-01_wp, -6.138631582149122E-02_wp, &
   &  4.390566938028945E-02_wp,  1.348102084384013E-02_wp, -1.422475405315354E-02_wp, &
   & -2.177708235331760E-03_wp, -1.631345209173193E-02_wp,  3.037854712508248E-04_wp, &
   & -5.982059491727628E-03_wp,  3.142698949422845E-03_wp,  6.011986702729309E-03_wp, &
   &  7.626778774808675E-03_wp, -5.988008281389998E-03_wp, -6.724463977621850E-03_wp, &
   & -1.586837277104500E-03_wp, -1.221174549502065E-01_wp,  6.230776427109196E-02_wp, &
   &  2.484894940597390E-02_wp, -1.200929750944426E-01_wp,  3.214388832129277E-01_wp, &
   & -3.961035306007012E-01_wp, -1.882538539354670E-01_wp,  2.789151759543163E-01_wp, &
   & -4.353801046840801E-03_wp, -2.129408121150662E-04_wp, -8.766086966417072E-03_wp, &
   &  9.088616371803124E-03_wp,  5.575618806522659E-03_wp,  2.764939982309179E-03_wp, &
   & -2.816460996223850E-04_wp,  1.335025177614989E-03_wp, -2.024238088144413E-01_wp, &
   &  1.349196812143939E-01_wp,  5.590840217582008E-02_wp, -1.309274647645181E-01_wp, &
   &  9.636569924340352E-03_wp, -7.369372293629420E-02_wp, -2.841652071952432E-02_wp, &
   &  1.047368757642663E-01_wp,  3.530348390705637E-03_wp,  1.580861854842841E-03_wp, &
   &  4.713874274092306E-03_wp,  5.331831361425191E-03_wp, -1.825513260223358E-03_wp, &
   & -8.231748457816464E-03_wp,  1.976050064171017E-03_wp, -3.757103059626676E-03_wp, &
   & -4.719659729581716E-03_wp, -2.043925312956940E-03_wp,  3.614820674096364E-02_wp, &
   & -4.486579619528315E-04_wp, -3.898846109092222E-03_wp, -3.608821299372046E-03_wp, &
   &  4.973625140154048E-01_wp, -5.965491465948145E-02_wp,  7.456698243138926E-03_wp, &
   &  5.929519221641824E-03_wp, -6.834988463389063E-03_wp,  1.630850450541114E-03_wp, &
   &  2.602638508272352E-04_wp,  5.994359611709595E-03_wp,  3.016288542567534E-03_wp, &
   & -6.034873875218512E-03_wp,  1.364409319729182E-02_wp,  1.431458658200893E-02_wp, &
   & -2.069210923895054E-03_wp,  1.604383983791018E-02_wp,  4.500960801472729E-02_wp, &
   &  1.089361360416053E-01_wp, -6.122161828353037E-02_wp, -4.339848056135150E-02_wp, &
   &  5.530047037939930E-03_wp, -2.676768520472069E-03_wp, -3.256242685531711E-04_wp, &
   & -1.149768177310682E-03_wp, -4.647040632469871E-03_wp,  2.699674398540623E-04_wp, &
   & -8.879099182546374E-03_wp, -9.434213232497941E-03_wp,  3.222313839815492E-01_wp, &
   &  3.955559494071835E-01_wp, -1.878748571960719E-01_wp, -2.787259594675615E-01_wp, &
   & -1.220033309745911E-01_wp, -6.228592308352851E-02_wp,  2.472409049003296E-02_wp, &
   &  1.207716734700453E-01_wp, -1.699130993047991E-03_wp,  8.535738856329883E-03_wp, &
   &  1.905012656259403E-03_wp,  3.456342577874180E-03_wp,  3.495014311407543E-03_wp, &
   & -1.567395487485505E-03_wp,  4.746948076908401E-03_wp, -5.296451585602324E-03_wp, &
   &  9.228907270398393E-03_wp,  7.303210790879353E-02_wp, -2.808387708646047E-02_wp, &
   & -1.037016869246249E-01_wp, -2.025846191385777E-01_wp, -1.363381950886083E-01_wp, &
   &  5.647603482094272E-02_wp,  1.294866496059647E-01_wp, -2.043891119358536E-03_wp, &
   & -4.717513073777491E-03_wp, -4.438882658831592E-04_wp,  3.614820761444790E-02_wp, &
   & -3.605864269071998E-03_wp, -3.898189592416248E-03_wp, -5.965491465948145E-02_wp, &
   &  4.973538962562015E-01_wp,  2.595589979218433E-04_wp,  5.993735420793160E-03_wp, &
   &  3.017483840519584E-03_wp,  6.035105234422779E-03_wp,  7.458746853310369E-03_wp, &
   &  5.929788372724423E-03_wp, -6.833731434090051E-03_wp, -1.632421248836141E-03_wp, &
   &  4.501673538816142E-02_wp,  1.089298169099016E-01_wp, -6.121853700392915E-02_wp, &
   &  4.339936991984415E-02_wp,  1.364350653474842E-02_wp,  1.431305862501090E-02_wp, &
   & -2.069022477564587E-03_wp, -1.604219428636813E-02_wp, -4.642999498868949E-03_wp, &
   &  2.694244657963935E-04_wp, -8.878478256508787E-03_wp,  9.430726882578795E-03_wp, &
   &  5.529488195520937E-03_wp, -2.677493390270879E-03_wp, -3.248024266876738E-04_wp, &
   &  1.150998895518944E-03_wp, -1.220084238782400E-01_wp, -6.228527908579715E-02_wp, &
   &  2.472472171481012E-02_wp, -1.207719456799169E-01_wp,  3.222291498982033E-01_wp, &
   &  3.955522415388684E-01_wp, -1.878721411889985E-01_wp,  2.787181662448878E-01_wp, &
   &  3.496911942846001E-03_wp, -1.567996147383244E-03_wp,  4.747446067989234E-03_wp, &
   &  5.297970446505873E-03_wp, -1.698159717415534E-03_wp,  8.535313896577645E-03_wp, &
   &  1.902395363242915E-03_wp, -3.452610340164516E-03_wp, -2.025736978208567E-01_wp, &
   & -1.363325014206645E-01_wp,  5.647403804236660E-02_wp, -1.294856726660747E-01_wp, &
   &  9.234596091040548E-03_wp,  7.303391508409832E-02_wp, -2.808456425118116E-02_wp, &
   &  1.036973409355553E-01_wp,  3.215652341046649E-01_wp, -1.226814369553450E-01_wp, &
   &  5.478947603167341E-03_wp, -4.727944164665414E-03_wp,  1.348077095700679E-02_wp, &
   &  4.533819980979545E-02_wp,  7.456698243138926E-03_wp,  2.595589979218433E-04_wp, &
   &  7.652548804759316E-01_wp,  9.011113009691714E-02_wp,  5.816373824724456E-02_wp, &
   & -3.810986804827858E-03_wp,  9.977579253307377E-02_wp, -4.969041233017572E-02_wp, &
   & -7.850717761760965E-03_wp,  2.680674855097192E-01_wp, -1.365031801248799E-02_wp, &
   & -8.385831999992309E-04_wp, -6.387612270062180E-03_wp,  3.271310642066478E-03_wp, &
   & -1.271757437904245E-02_wp,  2.538353402630409E-02_wp,  6.618028957324348E-04_wp, &
   &  1.741922551530340E-02_wp, -1.176725613967957E-01_wp, -6.266923413405520E-02_wp, &
   & -3.827659664859119E-02_wp, -4.061416741639562E-02_wp, -5.676640709098192E-02_wp, &
   & -3.928322745863190E-02_wp, -1.682168622246193E-02_wp, -4.446532505059938E-03_wp, &
   &  1.081665546945541E-02_wp, -2.382727643844985E-03_wp,  1.353447597982711E-02_wp, &
   & -1.943350706591306E-03_wp, -2.525532319575198E-02_wp, -7.140781347833789E-03_wp, &
   & -4.077573607576199E-03_wp,  4.528041766467475E-03_wp, -2.551968871121886E-02_wp, &
   &  1.557099978822604E-01_wp,  8.986905531255818E-02_wp, -1.796270662500462E-01_wp, &
   &  2.845113839599321E-02_wp, -7.332375822996767E-03_wp,  2.419272907748681E-03_wp, &
   & -1.739403599517548E-01_wp,  7.804943062986344E-03_wp, -4.325931302248842E-03_wp, &
   & -1.164661813543127E-02_wp,  8.205582248790529E-03_wp, -9.863558958445527E-02_wp, &
   &  1.864218687887590E-02_wp,  1.685404835755008E-02_wp,  1.014979472509048E-01_wp, &
   &  3.949845594623018E-01_wp, -6.235649514184920E-02_wp,  2.662859497750907E-03_wp, &
   & -3.786005378653440E-04_wp, -1.422486198769357E-02_wp, -1.095739336851609E-01_wp, &
   &  5.929519221641824E-03_wp,  5.993735420793161E-03_wp,  9.011113009691714E-02_wp, &
   &  6.981933186750613E-01_wp, -1.909107182270493E-03_wp,  5.766345658030282E-02_wp, &
   & -4.969048169587328E-02_wp,  9.605938118212901E-02_wp, -1.728686037915520E-01_wp, &
   & -1.415489149731182E-02_wp,  7.455032270066268E-04_wp,  2.016497841158198E-04_wp, &
   & -1.944142628172348E-03_wp, -3.904162851399948E-03_wp, -2.499457886108047E-02_wp, &
   &  2.600264462723397E-02_wp,  2.927867900513289E-03_wp,  2.504845943564731E-02_wp, &
   &  6.294328682808437E-02_wp, -3.237302861018303E-02_wp,  1.404742416438065E-01_wp, &
   &  7.468996494821596E-02_wp,  3.917135434981721E-02_wp,  1.164844574162803E-01_wp, &
   &  1.120779559391549E-01_wp, -1.520330401052979E-03_wp, -2.156391535833864E-03_wp, &
   & -5.517827719244788E-03_wp,  2.294197713089968E-02_wp,  5.398441203417207E-03_wp, &
   & -7.248945552398283E-03_wp,  4.733761954882748E-03_wp,  4.686288305285726E-04_wp, &
   &  8.189091634808692E-03_wp, -1.310796772062899E-01_wp, -2.736718176449321E-01_wp, &
   & -3.712663813199765E-01_wp,  2.912542799115598E-01_wp,  7.242104618418099E-02_wp, &
   &  3.615075610362662E-02_wp,  2.884311599600671E-02_wp,  1.927387162567067E-02_wp, &
   & -2.538903847321170E-03_wp, -8.372682940089578E-04_wp, -8.593498681167266E-03_wp, &
   &  5.747215247280658E-03_wp, -5.981308331902151E-02_wp,  1.301470143766196E-02_wp, &
   &  8.957318822404307E-03_wp,  6.321787876977525E-02_wp,  1.877266064541833E-01_wp, &
   & -2.468847162697405E-02_wp,  3.560082383852522E-04_wp,  8.753238377034695E-03_wp, &
   & -2.177761207993461E-03_wp, -6.138631582149122E-02_wp, -6.834988463389063E-03_wp, &
   &  3.017483840519586E-03_wp,  5.816373824724456E-02_wp, -1.909107182270493E-03_wp, &
   &  7.413061096584009E-01_wp,  3.903881899085133E-02_wp, -7.850688989647627E-03_wp, &
   & -1.728684968762369E-01_wp,  4.370992964277925E-01_wp, -1.624189519282423E-02_wp, &
   & -6.396904583426172E-03_wp,  1.901842921261877E-03_wp,  3.744928412373609E-03_wp, &
   &  3.557615641350124E-03_wp,  6.110983579137473E-04_wp, -2.670956385410298E-03_wp, &
   &  1.515430962282247E-02_wp,  7.957655760353795E-04_wp,  3.871829409834249E-02_wp, &
   &  1.403038387581643E-01_wp, -3.076280567036424E-01_wp,  2.311774941603729E-02_wp, &
   &  1.691130145172703E-02_wp,  1.118345868652891E-01_wp, -8.605644667181840E-02_wp, &
   &  4.687733005224171E-03_wp, -1.321534025919795E-02_wp, -2.269316060642431E-02_wp, &
   &  1.871304601260290E-02_wp,  1.461370864808542E-02_wp,  3.994431557940251E-03_wp, &
   & -4.770325680931548E-04_wp, -1.070308403456159E-03_wp, -3.022108110398367E-03_wp, &
   & -8.243877639499057E-02_wp, -3.621129115947413E-01_wp,  3.583325159333834E-01_wp, &
   &  1.489142714551477E-01_wp, -3.813695409354407E-03_wp,  3.259378511984293E-02_wp, &
   & -1.363227504153741E-02_wp, -3.386000962854872E-03_wp,  1.300248247308219E-03_wp, &
   & -5.641419553982402E-03_wp, -1.665767637263845E-03_wp,  8.506711270493465E-03_wp, &
   & -2.300906688145581E-02_wp,  8.242867184823468E-03_wp,  6.843528572579309E-03_wp, &
   &  2.284068750337667E-02_wp,  2.805496228740723E-01_wp,  1.208993706003673E-01_wp, &
   &  1.268974885733111E-03_wp, -9.319573280149776E-03_wp,  1.631366894119213E-02_wp, &
   &  4.390566938028945E-02_wp,  1.630850450541114E-03_wp,  6.035105234422779E-03_wp, &
   & -3.810986804827856E-03_wp,  5.766345658030282E-02_wp,  3.903881899085133E-02_wp, &
   &  6.446570157793475E-01_wp, -2.680667678174304E-01_wp,  1.415469007851684E-02_wp, &
   &  1.624169878959325E-02_wp, -5.068422655790606E-01_wp,  3.289415666046548E-03_wp, &
   &  3.915339536663696E-03_wp,  3.568652924204524E-03_wp,  1.615845114342901E-03_wp, &
   & -1.735561385890670E-02_wp,  2.520168181443818E-02_wp, -7.062358149822860E-04_wp, &
   &  1.896662531270825E-02_wp,  3.999767452904321E-02_wp,  7.385684969425246E-02_wp, &
   &  2.301656662585119E-02_wp, -4.690265114230360E-02_wp, -4.568530638601414E-03_wp, &
   &  1.260822103582827E-03_wp, -4.673202978761417E-03_wp,  7.706322331090823E-03_wp, &
   &  1.979348048503108E-03_wp, -5.455464238544920E-03_wp,  1.479408949786869E-02_wp, &
   & -4.704712498599953E-03_wp,  4.434271568564158E-03_wp,  8.155084996599015E-03_wp, &
   &  3.054875561499033E-03_wp,  1.357188479253855E-02_wp,  4.380665652887400E-02_wp, &
   &  2.164122060185590E-01_wp,  9.992231915273242E-02_wp, -1.860131997232384E-01_wp, &
   & -5.800049750378886E-03_wp,  5.536611255510542E-02_wp,  1.517387819873661E-02_wp, &
   &  2.147301316072957E-01_wp, -2.533085154898232E-02_wp,  2.977890074213256E-03_wp, &
   & -3.429125807313961E-03_wp, -6.449076849290553E-03_wp, -7.088130909614150E-02_wp, &
   &  1.683742740209905E-02_wp,  2.034317339994986E-03_wp,  6.327359645025230E-02_wp, &
   & -1.226809308310684E-01_wp,  3.215654129343700E-01_wp, -4.723925904439215E-03_wp, &
   &  5.479481708318162E-03_wp,  4.533861073404212E-02_wp,  1.348102084384013E-02_wp, &
   &  2.602638508272352E-04_wp,  7.458746853310369E-03_wp,  9.977579253307377E-02_wp, &
   & -4.969048169587328E-02_wp, -7.850688989647627E-03_wp, -2.680667678174304E-01_wp, &
   &  7.652542147071683E-01_wp,  9.011089921769284E-02_wp,  5.816369977169853E-02_wp, &
   &  3.809658207773107E-03_wp, -1.271661267336451E-02_wp,  2.538177417988710E-02_wp, &
   &  6.626828046074034E-04_wp, -1.741763285657855E-02_wp, -1.365207447542371E-02_wp, &
   & -8.376918796911199E-04_wp, -6.388006078369474E-03_wp, -3.271139003729714E-03_wp, &
   & -5.676648055065017E-02_wp, -3.928338545144464E-02_wp, -1.682148564197242E-02_wp, &
   &  4.446549979852882E-03_wp, -1.176725614866071E-01_wp, -6.266887430989296E-02_wp, &
   & -3.827636804163989E-02_wp,  4.061440131228599E-02_wp, -2.525710611760287E-02_wp, &
   & -7.142392630901574E-03_wp, -4.076781715761658E-03_wp, -4.526852586364580E-03_wp, &
   &  1.081588729290596E-02_wp, -2.384170357324099E-03_wp,  1.353505807534879E-02_wp, &
   &  1.942893511781486E-03_wp,  2.845024133171710E-02_wp, -7.332213225838278E-03_wp, &
   &  2.419222393111799E-03_wp,  1.739388997798245E-01_wp, -2.552041991790559E-02_wp, &
   &  1.557102188120999E-01_wp,  8.986925833650755E-02_wp,  1.796276237353349E-01_wp, &
   & -9.863723434119841E-02_wp,  1.864332326969586E-02_wp,  1.685365810430683E-02_wp, &
   & -1.014984336142544E-01_wp,  7.805716983696948E-03_wp, -4.325915303674373E-03_wp, &
   & -1.164665903359848E-02_wp, -8.205281927793753E-03_wp, -6.235703051699102E-02_wp, &
   &  3.949847425377017E-01_wp, -3.780206736377426E-04_wp,  2.662153218083750E-03_wp, &
   & -1.095742599472761E-01_wp, -1.422475405315354E-02_wp,  5.994359611709595E-03_wp, &
   &  5.929788372724423E-03_wp, -4.969041233017572E-02_wp,  9.605938118212901E-02_wp, &
   & -1.728684968762369E-01_wp,  1.415469007851684E-02_wp,  9.011089921769284E-02_wp, &
   &  6.981944691435855E-01_wp, -1.907787141053764E-03_wp, -5.766302883292199E-02_wp, &
   & -2.499181127332156E-02_wp,  2.600013710838952E-02_wp,  2.929189598464992E-03_wp, &
   & -2.504603203560259E-02_wp,  7.450574745177401E-04_wp,  2.020636007006103E-04_wp, &
   & -1.944265591397268E-03_wp,  3.903117843780120E-03_wp,  3.917149979601663E-02_wp, &
   &  1.164858179226037E-01_wp,  1.120760805058622E-01_wp,  1.519884169794081E-03_wp, &
   &  6.294365115732777E-02_wp, -3.237303057718443E-02_wp,  1.404741559814788E-01_wp, &
   & -7.468997191983635E-02_wp, -7.248937186017703E-03_wp,  4.733589079461940E-03_wp, &
   &  4.686365515341698E-04_wp, -8.189599121879210E-03_wp, -2.155701853227797E-03_wp, &
   & -5.518611163450458E-03_wp,  2.294243295262437E-02_wp, -5.397665312309056E-03_wp, &
   &  7.242096930076826E-02_wp,  3.615078167594412E-02_wp,  2.884443599340384E-02_wp, &
   & -1.927460744094261E-02_wp, -1.310793590904125E-01_wp, -2.736720740635467E-01_wp, &
   & -3.712668920587331E-01_wp, -2.912548161985849E-01_wp, -5.981538476871626E-02_wp, &
   &  1.301553465986978E-02_wp,  8.956889561815496E-03_wp, -6.321907174909422E-02_wp, &
   & -2.539141675934529E-03_wp, -8.379272217330506E-04_wp, -8.593342563929714E-03_wp, &
   & -5.746306227340858E-03_wp, -2.468880734654038E-02_wp,  1.877266874435544E-01_wp, &
   &  8.752624480373897E-03_wp,  3.568278630241579E-04_wp, -6.138654825322425E-02_wp, &
   & -2.177708235331758E-03_wp,  3.016288542567534E-03_wp, -6.833731434090051E-03_wp, &
   & -7.850717761760965E-03_wp, -1.728686037915520E-01_wp,  4.370992964277925E-01_wp, &
   &  1.624169878959326E-02_wp,  5.816369977169853E-02_wp, -1.907787141053764E-03_wp, &
   &  7.413051723049228E-01_wp, -3.903853820096768E-02_wp,  6.100966513710596E-04_wp, &
   & -2.668441239887111E-03_wp,  1.515304925386562E-02_wp, -7.964942445096990E-04_wp, &
   & -6.396597997644129E-03_wp,  1.901965081265378E-03_wp,  3.744905056558747E-03_wp, &
   & -3.557342388154321E-03_wp,  1.691109495800374E-02_wp,  1.118327127908419E-01_wp, &
   & -8.605189328272377E-02_wp, -4.687748761840149E-03_wp,  3.871852186723950E-02_wp, &
   &  1.403039356289937E-01_wp, -3.076280678749312E-01_wp, -2.311781510475472E-02_wp, &
   &  3.993925220480317E-03_wp, -4.772392440508207E-04_wp, -1.070229159507171E-03_wp, &
   &  3.021112433352000E-03_wp, -1.321253458196261E-02_wp, -2.268954798545994E-02_wp, &
   &  1.871131041721295E-02_wp, -1.461153704388195E-02_wp, -3.813587344754937E-03_wp, &
   &  3.259515165671301E-02_wp, -1.363462500872640E-02_wp,  3.385748079643317E-03_wp, &
   & -8.243836686100103E-02_wp, -3.621132033005950E-01_wp,  3.583332212693929E-01_wp, &
   & -1.489146544218456E-01_wp, -2.301012739452538E-02_wp,  8.240915671879087E-03_wp, &
   &  6.844292821216313E-03_wp, -2.284112114653753E-02_wp,  1.300202442455607E-03_wp, &
   & -5.640644943804287E-03_wp, -1.666141520371214E-03_wp, -8.506322112167258E-03_wp, &
   & -1.208993709601304E-01_wp, -2.805489103522327E-01_wp,  9.316120032689994E-03_wp, &
   & -1.267766929892562E-03_wp, -4.390566297627865E-02_wp, -1.631345209173193E-02_wp, &
   & -6.034873875218512E-03_wp, -1.632421248836141E-03_wp,  2.680674855097192E-01_wp, &
   & -1.415489149731182E-02_wp, -1.624189519282423E-02_wp, -5.068422655790606E-01_wp, &
   &  3.809658207773107E-03_wp, -5.766302883292199E-02_wp, -3.903853820096768E-02_wp, &
   &  6.446557665702752E-01_wp,  1.735381890505609E-02_wp, -2.519998642506840E-02_wp, &
   &  7.053852889883534E-04_wp,  1.896466816081822E-02_wp, -3.289468198094815E-03_wp, &
   & -3.915625267464379E-03_wp, -3.568450626889205E-03_wp,  1.616932918705391E-03_wp, &
   &  4.568541832258806E-03_wp, -1.260377757965680E-03_wp,  4.673222407782869E-03_wp, &
   &  7.705877872698733E-03_wp, -3.999744105521979E-02_wp, -7.385684388230018E-02_wp, &
   & -2.301650449419141E-02_wp, -4.690265128373323E-02_wp, -4.435351612508362E-03_wp, &
   & -8.155828201260798E-03_wp, -3.054655989632594E-03_wp,  1.357272418942471E-02_wp, &
   & -1.979007091270055E-03_wp,  5.456380640725584E-03_wp, -1.479433721130869E-02_wp, &
   & -4.705408422851859E-03_wp,  5.799744495423170E-03_wp, -5.536604633567437E-02_wp, &
   & -1.517379452047215E-02_wp,  2.147294740426981E-01_wp, -4.380651862943320E-02_wp, &
   & -2.164121130792328E-01_wp, -9.992231574193339E-02_wp, -1.860132438716162E-01_wp, &
   &  7.088307240209218E-02_wp, -1.683846527489821E-02_wp, -2.034015039736788E-03_wp, &
   &  6.327477665737241E-02_wp,  2.533132611978335E-02_wp, -2.977223798203978E-03_wp, &
   &  3.428901092250826E-03_wp, -6.449877852665542E-03_wp,  5.498706220187589E-03_wp, &
   & -4.287430237917706E-03_wp,  3.222353137213799E-01_wp, -1.220031725743797E-01_wp, &
   &  7.625658910081461E-03_wp,  3.037854712508248E-04_wp,  1.364409319729182E-02_wp, &
   &  4.501673538816142E-02_wp, -1.365031801248799E-02_wp,  7.455032270066266E-04_wp, &
   & -6.396904583426172E-03_wp,  3.289415666046548E-03_wp, -1.271661267336451E-02_wp, &
   & -2.499181127332156E-02_wp,  6.100966513710596E-04_wp,  1.735381890505609E-02_wp, &
   &  7.682438507271151E-01_wp, -9.102557409115047E-02_wp,  5.863868124164489E-02_wp, &
   & -1.828197450038696E-03_wp,  9.698989293743471E-02_wp,  5.001905546644120E-02_wp, &
   & -7.790567700853807E-03_wp,  2.659189700936845E-01_wp,  1.099936125805206E-02_wp, &
   &  2.222100732782892E-03_wp,  1.347503569396740E-02_wp, -2.125415548578480E-03_wp, &
   & -2.532533083378794E-02_wp,  7.285778277825599E-03_wp, -4.053994783348412E-03_wp, &
   &  4.523362039197296E-03_wp, -1.184597708358510E-01_wp,  6.289158309803693E-02_wp, &
   & -3.864036874936527E-02_wp, -4.087198470480553E-02_wp, -5.726299452013267E-02_wp, &
   &  3.921417238547908E-02_wp, -1.688741882898376E-02_wp, -4.553650550635379E-03_wp, &
   &  7.619742304759612E-03_wp,  4.223219909973934E-03_wp, -1.167476765732309E-02_wp, &
   &  7.880149076208952E-03_wp, -9.845927225922103E-02_wp, -1.850913627906792E-02_wp, &
   &  1.681748597117476E-02_wp,  1.011023440756921E-01_wp, -2.369189222685061E-02_wp, &
   & -1.558062066062249E-01_wp,  8.988384688691844E-02_wp, -1.817969489626665E-01_wp, &
   &  2.878518138124739E-02_wp,  7.063319593985584E-03_wp,  2.684049542240108E-03_wp, &
   & -1.741285000906176E-01_wp, -2.682407928818319E-03_wp,  2.041263539903967E-04_wp, &
   & -3.955576278448117E-01_wp,  6.229101417531866E-02_wp, -5.986317947460313E-03_wp, &
   & -5.982059491727627E-03_wp,  1.431458658200893E-02_wp,  1.089298169099016E-01_wp, &
   & -8.385831999992309E-04_wp,  2.016497841158197E-04_wp,  1.901842921261877E-03_wp, &
   &  3.915339536663696E-03_wp,  2.538177417988710E-02_wp,  2.600013710838952E-02_wp, &
   & -2.668441239887111E-03_wp, -2.519998642506840E-02_wp, -9.102557409115047E-02_wp, &
   &  6.974283951825311E-01_wp,  1.374898864819328E-03_wp, -5.770459098428200E-02_wp, &
   &  5.002495036194888E-02_wp,  9.243869095911050E-02_wp,  1.711520292633035E-01_wp, &
   &  1.398625602132661E-02_wp,  2.230742411596002E-03_wp, -5.490598181416147E-03_wp, &
   & -2.292847571360686E-02_wp, -5.427468452794968E-03_wp,  7.213630957379532E-03_wp, &
   &  4.745391244676303E-03_wp, -5.064257343623238E-04_wp, -8.222919454152797E-03_wp, &
   & -6.287832318848177E-02_wp, -3.255395128101228E-02_wp, -1.408816241791581E-01_wp, &
   & -7.471691933004269E-02_wp, -3.922119341459263E-02_wp,  1.160414267458867E-01_wp, &
   & -1.122392498472225E-01_wp,  1.352397637141819E-03_wp,  2.590228232093147E-03_wp, &
   & -7.729410252126835E-04_wp,  8.659799042456904E-03_wp, -5.802103710602418E-03_wp, &
   &  5.936510222427418E-02_wp,  1.282726726055198E-02_wp, -8.754346717420952E-03_wp, &
   & -6.279365091767158E-02_wp,  1.327559948026534E-01_wp, -2.716716161450161E-01_wp, &
   &  3.731497359428295E-01_wp, -2.914955018238549E-01_wp, -7.213726217635565E-02_wp, &
   &  3.599017371484808E-02_wp, -2.878545675566218E-02_wp, -1.934774622854778E-02_wp, &
   &  2.693644321560950E-04_wp,  8.751701679502096E-03_wp,  1.878861062451542E-01_wp, &
   & -2.471988022804309E-02_wp, -6.725184651742494E-03_wp,  3.142698949422847E-03_wp, &
   & -2.069210923895054E-03_wp, -6.121853700392915E-02_wp, -6.387612270062179E-03_wp, &
   & -1.944142628172348E-03_wp,  3.744928412373609E-03_wp,  3.568652924204524E-03_wp, &
   &  6.626828046074034E-04_wp,  2.929189598464992E-03_wp,  1.515304925386562E-02_wp, &
   &  7.053852889883500E-04_wp,  5.863868124164489E-02_wp,  1.374898864819328E-03_wp, &
   &  7.396511576803894E-01_wp,  3.917359720529495E-02_wp, -7.792833946951926E-03_wp, &
   &  1.711515263277960E-01_wp,  4.308136998618690E-01_wp, -1.622826245827758E-02_wp, &
   & -1.348557840054358E-02_wp,  2.291881278196172E-02_wp,  1.886138231619269E-02_wp, &
   &  1.478720529204197E-02_wp,  3.960966670450329E-03_wp,  4.981560665251116E-04_wp, &
   & -1.200450866205804E-03_wp, -2.988513172102313E-03_wp,  3.863424072780783E-02_wp, &
   & -1.408760185979790E-01_wp, -3.088581278052090E-01_wp,  2.324793300914090E-02_wp, &
   &  1.690736399359133E-02_wp, -1.122364450025662E-01_wp, -8.652690893444472E-02_wp, &
   &  4.619687899269517E-03_wp,  1.291179652397112E-03_wp,  5.740911507820033E-03_wp, &
   & -1.714098005125333E-03_wp,  8.529977331087440E-03_wp, -2.278501721523940E-02_wp, &
   & -8.212466815202671E-03_wp,  6.820742093190971E-03_wp,  2.266829586314910E-02_wp, &
   & -8.322877082234913E-02_wp,  3.640456085580128E-01_wp,  3.635064035974106E-01_wp, &
   &  1.490313042946019E-01_wp, -4.040031653327232E-03_wp, -3.253383232104642E-02_wp, &
   & -1.360134300980547E-02_wp, -3.378510275726214E-03_wp,  1.303015842936295E-03_wp, &
   & -9.013011234728523E-03_wp,  2.787480693862803E-01_wp,  1.207793448348777E-01_wp, &
   &  1.586802840367607E-03_wp,  6.011986702729309E-03_wp,  1.604383983791018E-02_wp, &
   &  4.339936991984415E-02_wp,  3.271310642066478E-03_wp, -3.904162851399948E-03_wp, &
   &  3.557615641350124E-03_wp,  1.615845114342901E-03_wp, -1.741763285657855E-02_wp, &
   & -2.504603203560259E-02_wp, -7.964942445096990E-04_wp,  1.896466816081822E-02_wp, &
   & -1.828197450038696E-03_wp, -5.770459098428200E-02_wp,  3.917359720529495E-02_wp, &
   &  6.469361463689909E-01_wp, -2.659344401302082E-01_wp, -1.397627535382276E-02_wp, &
   &  1.622320187810272E-02_wp, -5.094810692800513E-01_wp,  2.152427869093051E-03_wp, &
   &  5.455485845626132E-03_wp,  1.481425697550504E-02_wp, -4.801974069113951E-03_wp, &
   &  4.496214070664370E-03_wp, -8.247967093940546E-03_wp,  3.094075457696282E-03_wp, &
   &  1.365520300857091E-02_wp,  4.085432781020713E-02_wp, -7.471553613655965E-02_wp, &
   &  2.325760370067884E-02_wp, -4.756254639370319E-02_wp, -4.548300829255204E-03_wp, &
   & -1.343524553838545E-03_wp, -4.603471265628236E-03_wp,  7.858653910191471E-03_wp, &
   & -2.530417551779322E-02_wp, -2.943663228912056E-03_wp, -3.463083877359813E-03_wp, &
   & -6.333394339464310E-03_wp, -7.072222891782816E-02_wp, -1.668584366919810E-02_wp, &
   &  1.876645488428163E-03_wp,  6.306833010156396E-02_wp,  4.287940936913998E-02_wp, &
   & -2.163725859997901E-01_wp,  1.000300771004840E-01_wp, -1.847694129517921E-01_wp, &
   & -7.216306434360623E-03_wp, -5.498940836255696E-02_wp,  1.489760162623975E-02_wp, &
   &  2.154876430367316E-01_wp, -4.287832154916927E-03_wp,  5.499073985967830E-03_wp, &
   & -1.220082606284287E-01_wp,  3.222375528530796E-01_wp,  3.046710348517129E-04_wp, &
   &  7.626778774808675E-03_wp,  4.500960801472729E-02_wp,  1.364350653474841E-02_wp, &
   & -1.271757437904245E-02_wp, -2.499457886108047E-02_wp,  6.110983579137468E-04_wp, &
   & -1.735561385890670E-02_wp, -1.365207447542371E-02_wp,  7.450574745177401E-04_wp, &
   & -6.396597997644129E-03_wp, -3.289468198094813E-03_wp,  9.698989293743471E-02_wp, &
   &  5.002495036194888E-02_wp, -7.792833946951926E-03_wp, -2.659344401302082E-01_wp, &
   &  7.682608314274988E-01_wp, -9.103581616419490E-02_wp,  5.864436321839947E-02_wp, &
   &  1.832440298511231E-03_wp, -2.532713202170716E-02_wp,  7.285763756416801E-03_wp, &
   & -4.053482619620960E-03_wp, -4.524457640091498E-03_wp,  1.100015613576399E-02_wp, &
   &  2.222798067329443E-03_wp,  1.347788452917695E-02_wp,  2.125772940411393E-03_wp, &
   & -5.725877587258311E-02_wp,  3.921484588210891E-02_wp, -1.688683597876312E-02_wp, &
   &  4.553102888578676E-03_wp, -1.184597738112737E-01_wp,  6.289372614147358E-02_wp, &
   & -3.864159531118652E-02_wp,  4.086846697200497E-02_wp, -9.847144852560728E-02_wp, &
   & -1.851123192800146E-02_wp,  1.681962644965352E-02_wp, -1.011141159736219E-01_wp, &
   &  7.618980820543752E-03_wp,  4.223400132357042E-03_wp, -1.167579672914539E-02_wp, &
   & -7.880669598994257E-03_wp,  2.878693085465274E-02_wp,  7.064998248426048E-03_wp, &
   &  2.683597700459387E-03_wp,  1.741398116942549E-01_wp, -2.370309164556275E-02_wp, &
   & -1.557979903549532E-01_wp,  8.988113142672710E-02_wp,  1.817932770756827E-01_wp, &
   &  2.044003485932207E-04_wp, -2.682798239092303E-03_wp,  6.229036665474397E-02_wp, &
   & -3.955613432847502E-01_wp, -5.982650890830495E-03_wp, -5.988008281389998E-03_wp, &
   &  1.089361360416053E-01_wp,  1.431305862501090E-02_wp,  2.538353402630409E-02_wp, &
   &  2.600264462723397E-02_wp, -2.670956385410298E-03_wp,  2.520168181443818E-02_wp, &
   & -8.376918796911196E-04_wp,  2.020636007006104E-04_wp,  1.901965081265378E-03_wp, &
   & -3.915625267464379E-03_wp,  5.001905546644120E-02_wp,  9.243869095911050E-02_wp, &
   &  1.711515263277960E-01_wp, -1.397627535382276E-02_wp, -9.103581616419490E-02_wp, &
   &  6.974359477282280E-01_wp,  1.367490431307298E-03_wp,  5.770281101804980E-02_wp, &
   &  7.215243499261652E-03_wp,  4.745218750400342E-03_wp, -5.066319636012327E-04_wp, &
   &  8.223668426745046E-03_wp,  2.229291255920331E-03_wp, -5.489802225053643E-03_wp, &
   & -2.293214173344136E-02_wp,  5.426543864703204E-03_wp, -3.922186295684765E-02_wp, &
   &  1.160506387841913E-01_wp, -1.122402208931718E-01_wp, -1.354358795737720E-03_wp, &
   & -6.287618022372808E-02_wp, -3.255395287938539E-02_wp, -1.408814002950515E-01_wp, &
   &  7.471791946172623E-02_wp,  5.937310100667977E-02_wp,  1.282979193187429E-02_wp, &
   & -8.756560241211521E-03_wp,  6.280096415787730E-02_wp,  2.591875878241655E-03_wp, &
   & -7.736922218217283E-04_wp,  8.662091590280946E-03_wp,  5.801361411158843E-03_wp, &
   & -7.213721278339691E-02_wp,  3.598800667611273E-02_wp, -2.878459476131501E-02_wp, &
   &  1.934151879091920E-02_wp,  1.327614522838112E-01_wp, -2.716840289418065E-01_wp, &
   &  3.731547362524139E-01_wp,  2.914978418533170E-01_wp,  8.751385748047206E-03_wp, &
   &  2.695262479973700E-04_wp, -2.472051104696680E-02_wp,  1.878888234682954E-01_wp, &
   &  3.142870872497287E-03_wp, -6.724463977621848E-03_wp, -6.122161828353037E-02_wp, &
   & -2.069022477564587E-03_wp,  6.618028957324348E-04_wp,  2.927867900513289E-03_wp, &
   &  1.515430962282247E-02_wp, -7.062358149822860E-04_wp, -6.388006078369474E-03_wp, &
   & -1.944265591397268E-03_wp,  3.744905056558747E-03_wp, -3.568450626889205E-03_wp, &
   & -7.790567700853807E-03_wp,  1.711520292633035E-01_wp,  4.308136998618690E-01_wp, &
   &  1.622320187810272E-02_wp,  5.864436321839947E-02_wp,  1.367490431307298E-03_wp, &
   &  7.396492985257906E-01_wp, -3.917253432483862E-02_wp,  3.960174899363565E-03_wp, &
   &  4.981634997720502E-04_wp, -1.200371732536248E-03_wp,  2.988292434367271E-03_wp, &
   & -1.348499348927481E-02_wp,  2.291834816195670E-02_wp,  1.886314429570616E-02_wp, &
   & -1.478695575275324E-02_wp,  1.690678139094201E-02_wp, -1.122374164166262E-01_wp, &
   & -8.651941146892095E-02_wp, -4.618229847397556E-03_wp,  3.863301465131361E-02_wp, &
   & -1.408762422576576E-01_wp, -3.088581282385118E-01_wp, -2.324896861371008E-02_wp, &
   & -2.278812893083784E-02_wp, -8.213618357573981E-03_wp,  6.821899693253411E-03_wp, &
   & -2.267115234887514E-02_wp,  1.290951888940569E-03_wp,  5.741259563876600E-03_wp, &
   & -1.715226971518527E-03_wp, -8.530080048747328E-03_wp, -4.038779749340381E-03_wp, &
   & -3.253315639511126E-02_wp, -1.360180780112587E-02_wp,  3.382014459363425E-03_wp, &
   & -8.323103634066595E-02_wp,  3.640522108596773E-01_wp,  3.635041436227353E-01_wp, &
   & -1.490334705749429E-01_wp,  9.014283332274178E-03_wp, -1.304167939683534E-03_wp, &
   & -1.207796158775903E-01_wp, -2.787558694019131E-01_wp, -6.012592350652365E-03_wp, &
   & -1.586837277104500E-03_wp, -4.339848056135150E-02_wp, -1.604219428636813E-02_wp, &
   &  1.741922551530340E-02_wp,  2.504845943564731E-02_wp,  7.957655760353795E-04_wp, &
   &  1.896662531270825E-02_wp, -3.271139003729715E-03_wp,  3.903117843780118E-03_wp, &
   & -3.557342388154318E-03_wp,  1.616932918705391E-03_wp,  2.659189700936845E-01_wp, &
   &  1.398625602132661E-02_wp, -1.622826245827758E-02_wp, -5.094810692800513E-01_wp, &
   &  1.832440298511231E-03_wp,  5.770281101804980E-02_wp, -3.917253432483862E-02_wp, &
   &  6.469341614537668E-01_wp, -4.495001945148239E-03_wp,  8.248479686085566E-03_wp, &
   & -3.093065411874579E-03_wp,  1.365604986778353E-02_wp, -2.152912146361097E-03_wp, &
   & -5.456266676215703E-03_wp, -1.481646099455012E-02_wp, -4.801281492113077E-03_wp, &
   &  4.547750849219394E-03_wp,  1.345483753142764E-03_wp,  4.602013044171604E-03_wp, &
   &  7.856583401107661E-03_wp, -4.085784348490161E-02_wp,  7.471453562286141E-02_wp, &
   & -2.325656785513532E-02_wp, -4.756254716429453E-02_wp,  7.073249475298077E-02_wp, &
   &  1.668755103203241E-02_wp, -1.878607371143480E-03_wp,  6.307787433385099E-02_wp, &
   &  2.530055882102593E-02_wp,  2.942679438949034E-03_wp,  3.464590689409829E-03_wp, &
   & -6.329093098400240E-03_wp,  7.215211206464081E-03_wp,  5.498834317754229E-02_wp, &
   & -1.489667751705982E-02_wp,  2.154832128200224E-01_wp, -4.287278436401225E-02_wp, &
   &  2.163694037783475E-01_wp, -1.000289298114973E-01_wp, -1.847699740988629E-01_wp, &
   &  1.315713396064210E-02_wp,  4.415770787953931E-02_wp,  7.535473075301793E-03_wp, &
   &  3.436186553753724E-04_wp,  3.214390604969148E-01_wp, -1.221174549502065E-01_wp, &
   &  5.530047037939930E-03_wp, -4.642999498868950E-03_wp, -1.176725613967957E-01_wp, &
   &  6.294328682808437E-02_wp,  3.871829409834249E-02_wp,  3.999767452904320E-02_wp, &
   & -5.676648055065017E-02_wp,  3.917149979601663E-02_wp,  1.691109495800374E-02_wp, &
   &  4.568541832258806E-03_wp,  1.099936125805206E-02_wp,  2.230742411596002E-03_wp, &
   & -1.348557840054358E-02_wp,  2.152427869093053E-03_wp, -2.532713202170716E-02_wp, &
   &  7.215243499261652E-03_wp,  3.960174899363565E-03_wp, -4.495001945148239E-03_wp, &
   &  7.633978057163138E-01_wp, -8.856601185970754E-02_wp, -5.769065171547247E-02_wp, &
   &  5.077798452434275E-03_wp,  9.773141860581004E-02_wp,  4.896967825493526E-02_wp, &
   &  7.349243681255485E-03_wp, -2.672299925715932E-01_wp, -1.390797301235984E-02_wp, &
   &  6.881743823780077E-04_wp,  6.423569787478840E-03_wp, -3.205742315796407E-03_wp, &
   & -1.229475779570280E-02_wp, -2.510285036163237E-02_wp, -7.394156897273328E-04_wp, &
   & -1.717362995183904E-02_wp,  2.929510346456513E-02_wp,  7.434275886919354E-03_wp, &
   & -2.443516507676706E-03_wp,  1.756037405016918E-01_wp, -2.418957168368160E-02_wp, &
   & -1.562434616964579E-01_wp, -9.035563306592274E-02_wp,  1.833298769665542E-01_wp, &
   & -9.884426245318094E-02_wp, -1.863273215096028E-02_wp, -1.681206142369586E-02_wp, &
   & -1.015205171001342E-01_wp,  7.727620842577443E-03_wp,  4.280987002069590E-03_wp, &
   &  1.158991091617597E-02_wp, -7.916934585009305E-03_wp,  1.421731639133684E-02_wp, &
   &  1.087881602220137E-01_wp, -6.085886474843062E-03_wp, -6.005550526819929E-03_wp, &
   & -3.961037152394013E-01_wp,  6.230776427109196E-02_wp, -2.676768520472069E-03_wp, &
   &  2.694244657963927E-04_wp, -6.266923413405520E-02_wp, -3.237302861018303E-02_wp, &
   &  1.403038387581643E-01_wp,  7.385684969425246E-02_wp, -3.928338545144464E-02_wp, &
   &  1.164858179226037E-01_wp,  1.118327127908419E-01_wp, -1.260377757965678E-03_wp, &
   &  2.222100732782892E-03_wp, -5.490598181416147E-03_wp,  2.291881278196172E-02_wp, &
   &  5.455485845626132E-03_wp,  7.285763756416801E-03_wp,  4.745218750400342E-03_wp, &
   &  4.981634997720500E-04_wp,  8.248479686085566E-03_wp, -8.856601185970754E-02_wp, &
   &  6.974547505826840E-01_wp, -9.365952876686562E-04_wp,  5.801189027210991E-02_wp, &
   &  4.896959824805414E-02_wp,  9.315785434449748E-02_wp, -1.689608278007086E-01_wp, &
   & -1.413706344437824E-02_wp, -7.943164815645875E-04_wp,  1.928757151626588E-04_wp, &
   & -1.972490835491246E-03_wp, -3.952082793087780E-03_wp,  2.532226254985261E-02_wp, &
   &  2.629472700284261E-02_wp,  2.804395208647878E-03_wp,  2.524346575856477E-02_wp, &
   & -7.215882064800733E-02_wp,  3.627696595311906E-02_wp,  2.875119908564139E-02_wp, &
   &  1.957479661351896E-02_wp,  1.309512863422543E-01_wp, -2.690443581569927E-01_wp, &
   & -3.743573887854515E-01_wp,  2.931831809220679E-01_wp,  5.975943572505005E-02_wp, &
   &  1.297006349097501E-02_wp,  8.960267273041791E-03_wp,  6.330195250473385E-02_wp, &
   &  2.589722239178300E-03_wp, -7.760709388883756E-04_wp, -8.590194921635699E-03_wp, &
   &  5.794296157964703E-03_wp,  1.982936508356516E-03_wp,  6.085384709963981E-02_wp, &
   &  6.870536941104054E-03_wp, -3.042063894358351E-03_wp, -1.882539383786511E-01_wp, &
   &  2.484894940597390E-02_wp, -3.256242685531693E-04_wp, -8.878478256508785E-03_wp, &
   & -3.827659664859119E-02_wp,  1.404742416438065E-01_wp, -3.076280567036424E-01_wp, &
   &  2.301656662585120E-02_wp, -1.682148564197242E-02_wp,  1.120760805058622E-01_wp, &
   & -8.605189328272377E-02_wp,  4.673222407782869E-03_wp,  1.347503569396740E-02_wp, &
   & -2.292847571360686E-02_wp,  1.886138231619269E-02_wp,  1.481425697550504E-02_wp, &
   & -4.053482619620960E-03_wp, -5.066319636012325E-04_wp, -1.200371732536248E-03_wp, &
   & -3.093065411874579E-03_wp, -5.769065171547247E-02_wp, -9.365952876686562E-04_wp, &
   &  7.392888419899679E-01_wp,  3.939144224208752E-02_wp,  7.349269093404890E-03_wp, &
   & -1.689609335868921E-01_wp,  4.273886601932486E-01_wp, -1.632080191159903E-02_wp, &
   &  6.446084723645063E-03_wp,  1.940044737722060E-03_wp,  3.821405390004228E-03_wp, &
   &  3.647333896716159E-03_wp, -6.423393741742480E-04_wp, -2.859958918711297E-03_wp, &
   &  1.524260432956563E-02_wp,  7.288003370454037E-04_wp,  3.618941053847149E-03_wp, &
   &  3.252597635589729E-02_wp, -1.327504017809231E-02_wp, -3.220333341428742E-03_wp, &
   &  8.258060529965380E-02_wp, -3.653724543429441E-01_wp,  3.680801231838481E-01_wp, &
   &  1.498598638686168E-01_wp,  2.288650062627777E-02_wp,  8.306702077661978E-03_wp, &
   &  6.864908474993436E-03_wp,  2.284663513307648E-02_wp, -1.244046434831659E-03_wp, &
   & -5.722674596992679E-03_wp, -1.672278219204951E-03_wp,  8.473545712606514E-03_wp, &
   & -1.638356792386353E-02_wp, -4.337273634581573E-02_wp, -1.590870018319618E-03_wp, &
   & -5.958499739953396E-03_wp, -2.789144572943701E-01_wp, -1.200929750944426E-01_wp, &
   & -1.149768177310682E-03_wp,  9.430726882578802E-03_wp, -4.061416741639562E-02_wp, &
   &  7.468996494821596E-02_wp,  2.311774941603729E-02_wp, -4.690265114230360E-02_wp, &
   &  4.446549979852882E-03_wp,  1.519884169794081E-03_wp, -4.687748761840145E-03_wp, &
   &  7.705877872698733E-03_wp, -2.125415548578480E-03_wp, -5.427468452794968E-03_wp, &
   &  1.478720529204197E-02_wp, -4.801974069113951E-03_wp, -4.524457640091498E-03_wp, &
   &  8.223668426745046E-03_wp,  2.988292434367275E-03_wp,  1.365604986778353E-02_wp, &
   &  5.077798452434277E-03_wp,  5.801189027210991E-02_wp,  3.939144224208752E-02_wp, &
   &  6.447250354712125E-01_wp,  2.672307169271428E-01_wp,  1.413726774321543E-02_wp, &
   &  1.632100214893510E-02_wp, -5.077523260433411E-01_wp, -3.464413185729152E-03_wp, &
   &  3.775733927542388E-03_wp,  3.561560280737855E-03_wp,  1.755662350481010E-03_wp, &
   &  1.717986995224270E-02_wp,  2.509731414543290E-02_wp, -7.003288826061063E-04_wp, &
   &  1.884198707659332E-02_wp,  6.534370866653386E-03_wp,  5.564409057226768E-02_wp, &
   &  1.531544562200643E-02_wp,  2.161460773258060E-01_wp, -4.447215084892581E-02_wp, &
   &  2.168683175933099E-01_wp,  1.004116692618975E-01_wp, -1.852658196507851E-01_wp, &
   &  7.107247800333402E-02_wp,  1.689252113263069E-02_wp,  2.089735104077940E-03_wp, &
   &  6.345913253029471E-02_wp,  2.527593529481101E-02_wp,  2.964544962364357E-03_wp, &
   & -3.421246430889103E-03_wp, -6.356288505501850E-03_wp,  4.415729811646429E-02_wp, &
   &  1.315687805394647E-02_wp,  3.428998781745260E-04_wp,  7.533427596972730E-03_wp, &
   & -1.221179632877208E-01_wp,  3.214388832129277E-01_wp, -4.647040632469871E-03_wp, &
   &  5.529488195520937E-03_wp, -5.676640709098192E-02_wp,  3.917135434981721E-02_wp, &
   &  1.691130145172703E-02_wp, -4.568530638601416E-03_wp, -1.176725614866071E-01_wp, &
   &  6.294365115732777E-02_wp,  3.871852186723949E-02_wp, -3.999744105521978E-02_wp, &
   & -2.532533083378794E-02_wp,  7.213630957379532E-03_wp,  3.960966670450329E-03_wp, &
   &  4.496214070664369E-03_wp,  1.100015613576399E-02_wp,  2.229291255920331E-03_wp, &
   & -1.348499348927481E-02_wp, -2.152912146361095E-03_wp,  9.773141860581004E-02_wp, &
   &  4.896959824805414E-02_wp,  7.349269093404890E-03_wp,  2.672307169271428E-01_wp, &
   &  7.633985219222896E-01_wp, -8.856623880928000E-02_wp, -5.769068496469183E-02_wp, &
   & -5.079137464534792E-03_wp, -1.229567321301728E-02_wp, -2.510457513442781E-02_wp, &
   & -7.385526035989136E-04_wp,  1.717517770848850E-02_wp, -1.390618957521301E-02_wp, &
   &  6.890824967592879E-04_wp,  6.423166689272697E-03_wp,  3.205947159947167E-03_wp, &
   & -2.418884425484044E-02_wp, -1.562432385840267E-01_wp, -9.035542567922990E-02_wp, &
   & -1.833293027304498E-01_wp,  2.929601302874185E-02_wp,  7.434446578429053E-03_wp, &
   & -2.443563686568657E-03_wp, -1.756052107039139E-01_wp,  7.726829452574939E-03_wp, &
   &  4.281002320965354E-03_wp,  1.158987035355509E-02_wp,  7.917211034894438E-03_wp, &
   & -9.884265042117255E-02_wp, -1.863159696218155E-02_wp, -1.681244984664289E-02_wp, &
   &  1.015200503060320E-01_wp,  1.087878338388485E-01_wp,  1.421742136850592E-02_wp, &
   & -6.004927066623618E-03_wp, -6.085629091190239E-03_wp,  6.230722519523442E-02_wp, &
   & -3.961035306007012E-01_wp,  2.699674398540623E-04_wp, -2.677493390270879E-03_wp, &
   & -3.928322745863190E-02_wp,  1.164844574162803E-01_wp,  1.118345868652891E-01_wp, &
   &  1.260822103582827E-03_wp, -6.266887430989296E-02_wp, -3.237303057718443E-02_wp, &
   &  1.403039356289937E-01_wp, -7.385684388230018E-02_wp,  7.285778277825599E-03_wp, &
   &  4.745391244676303E-03_wp,  4.981560665251116E-04_wp, -8.247967093940545E-03_wp, &
   &  2.222798067329443E-03_wp, -5.489802225053643E-03_wp,  2.291834816195670E-02_wp, &
   & -5.456266676215703E-03_wp,  4.896967825493526E-02_wp,  9.315785434449748E-02_wp, &
   & -1.689609335868921E-01_wp,  1.413726774321543E-02_wp, -8.856623880928000E-02_wp, &
   &  6.974535742708103E-01_wp, -9.378587107186976E-04_wp, -5.801231522156750E-02_wp, &
   &  2.532508239159577E-02_wp,  2.629728211238974E-02_wp,  2.803051313766203E-03_wp, &
   & -2.524593097806645E-02_wp, -7.947613217261478E-04_wp,  1.924630227055868E-04_wp, &
   & -1.972368471019810E-03_wp,  3.953144503594116E-03_wp,  1.309516163820550E-01_wp, &
   & -2.690441096373638E-01_wp, -3.743568770468897E-01_wp, -2.931826377618207E-01_wp, &
   & -7.215889877554434E-02_wp,  3.627695352670001E-02_wp,  2.874985981178474E-02_wp, &
   & -1.957406411145724E-02_wp,  2.589483077764269E-03_wp, -7.754048407027505E-04_wp, &
   & -8.590353289713753E-03_wp, -5.795216302978648E-03_wp,  5.975712114865252E-02_wp, &
   &  1.296921634983868E-02_wp,  8.960702273653885E-03_wp, -6.330075475231250E-02_wp, &
   &  6.085361576972199E-02_wp,  1.982987914676789E-03_wp, -3.043269941348813E-03_wp, &
   &  6.871785977987184E-03_wp,  2.484861103711767E-02_wp, -1.882538539354670E-01_wp, &
   & -8.879099182546374E-03_wp, -3.248024266876738E-04_wp, -1.682168622246193E-02_wp, &
   &  1.120779559391549E-01_wp, -8.605644667181840E-02_wp, -4.673202978761417E-03_wp, &
   & -3.827636804163989E-02_wp,  1.404741559814788E-01_wp, -3.076280678749312E-01_wp, &
   & -2.301650449419141E-02_wp, -4.053994783348412E-03_wp, -5.064257343623236E-04_wp, &
   & -1.200450866205804E-03_wp,  3.094075457696282E-03_wp,  1.347788452917695E-02_wp, &
   & -2.293214173344136E-02_wp,  1.886314429570616E-02_wp, -1.481646099455012E-02_wp, &
   &  7.349243681255486E-03_wp, -1.689608278007086E-01_wp,  4.273886601932486E-01_wp, &
   &  1.632100214893510E-02_wp, -5.769068496469183E-02_wp, -9.378587107186976E-04_wp, &
   &  7.392896501429503E-01_wp, -3.939172089031184E-02_wp, -6.433652749597172E-04_wp, &
   & -2.862518162679253E-03_wp,  1.524388513770821E-02_wp, -7.280562139118651E-04_wp, &
   &  6.446404476308520E-03_wp,  1.939924747404366E-03_wp,  3.821428149507959E-03_wp, &
   & -3.647612696143962E-03_wp,  8.258101038071362E-02_wp, -3.653721654972368E-01_wp, &
   &  3.680794146991755E-01_wp, -1.498594746533048E-01_wp,  3.619049329913439E-03_wp, &
   &  3.252458845540963E-02_wp, -1.327264103522864E-02_wp,  3.220587350774420E-03_wp, &
   & -1.244094321036358E-03_wp, -5.723459428157306E-03_wp, -1.671899435629965E-03_wp, &
   & -8.473940779935289E-03_wp,  2.288543553170518E-02_wp,  8.308680173787811E-03_wp, &
   &  6.864133072952403E-03_wp, -2.284620243063615E-02_wp,  4.337274631722113E-02_wp, &
   &  1.638378434837624E-02_wp,  5.958736679177928E-03_wp,  1.589292442514503E-03_wp, &
   &  1.200929772037539E-01_wp,  2.789151759543163E-01_wp, -9.434213232497941E-03_wp, &
   &  1.150998895518944E-03_wp, -4.446532505059938E-03_wp, -1.520330401052981E-03_wp, &
   &  4.687733005224171E-03_wp,  7.706322331090823E-03_wp,  4.061440131228599E-02_wp, &
   & -7.468997191983635E-02_wp, -2.311781510475472E-02_wp, -4.690265128373323E-02_wp, &
   &  4.523362039197296E-03_wp, -8.222919454152799E-03_wp, -2.988513172102317E-03_wp, &
   &  1.365520300857091E-02_wp,  2.125772940411394E-03_wp,  5.426543864703204E-03_wp, &
   & -1.478695575275324E-02_wp, -4.801281492113077E-03_wp, -2.672299925715932E-01_wp, &
   & -1.413706344437824E-02_wp, -1.632080191159903E-02_wp, -5.077523260433411E-01_wp, &
   & -5.079137464534792E-03_wp, -5.801231522156750E-02_wp, -3.939172089031184E-02_wp, &
   &  6.447263223074291E-01_wp, -1.718163176067052E-02_wp, -2.509898920227886E-02_wp, &
   &  7.011706174227667E-04_wp,  1.884393418196202E-02_wp,  3.464331170947334E-03_wp, &
   & -3.775467895130504E-03_wp, -3.561752932146639E-03_wp,  1.754560513036185E-03_wp, &
   &  4.447228593521083E-02_wp, -2.168684192675906E-01_wp, -1.004116743940034E-01_wp, &
   & -1.852657794787091E-01_wp, -6.534687037887052E-03_wp, -5.564415487573214E-02_wp, &
   & -1.531553131254082E-02_wp,  2.161467408179067E-01_wp, -2.527544829182120E-02_wp, &
   & -2.965209619472379E-03_wp,  3.421469452061979E-03_wp, -6.355489319584323E-03_wp, &
   & -7.107072288070147E-02_wp, -1.689148415298210E-02_wp, -2.090036052235476E-03_wp, &
   &  6.345794212690893E-02_wp,  7.676830732914519E-03_wp,  3.095258120151857E-04_wp, &
   &  1.364727763523009E-02_wp,  4.498941623563778E-02_wp,  5.575998655428776E-03_wp, &
   & -4.353801046840802E-03_wp,  3.222313839815492E-01_wp, -1.220084238782400E-01_wp, &
   &  1.081665546945541E-02_wp, -2.156391535833864E-03_wp, -1.321534025919795E-02_wp, &
   &  1.979348048503108E-03_wp, -2.525710611760287E-02_wp, -7.248937186017703E-03_wp, &
   &  3.993925220480317E-03_wp, -4.435351612508361E-03_wp, -1.184597708358510E-01_wp, &
   & -6.287832318848177E-02_wp,  3.863424072780783E-02_wp,  4.085432781020713E-02_wp, &
   & -5.725877587258311E-02_wp, -3.922186295684765E-02_wp,  1.690678139094201E-02_wp, &
   &  4.547750849219394E-03_wp, -1.390797301235984E-02_wp, -7.943164815645875E-04_wp, &
   &  6.446084723645063E-03_wp, -3.464413185729152E-03_wp, -1.229567321301728E-02_wp, &
   &  2.532508239159577E-02_wp, -6.433652749597172E-04_wp, -1.718163176067052E-02_wp, &
   &  7.683483967062564E-01_wp,  9.109679294113550E-02_wp, -5.864543153957989E-02_wp, &
   &  1.796960995305847E-03_wp,  9.693675898919289E-02_wp, -5.005208321108362E-02_wp, &
   &  7.791440829249158E-03_wp, -2.659276333325515E-01_wp, -9.864018496849172E-02_wp, &
   &  1.848916476504280E-02_wp, -1.681058788278299E-02_wp, -1.013484377614760E-01_wp, &
   &  7.667271490821401E-03_wp, -4.298131697609277E-03_wp,  1.162637117966820E-02_wp, &
   & -8.019162958230901E-03_wp,  2.878275285056010E-02_wp, -7.069186144772829E-03_wp, &
   & -2.697581873754085E-03_wp,  1.741109383746041E-01_wp, -2.370761223236972E-02_wp, &
   &  1.558022075373159E-01_wp, -8.986673754216588E-02_wp,  1.817765936749479E-01_wp, &
   &  6.010314084803126E-03_wp,  5.983459493418869E-03_wp, -1.431013175074783E-02_wp, &
   & -1.089524752117987E-01_wp,  2.765339311704463E-03_wp, -2.129408121150670E-04_wp, &
   &  3.955559494071835E-01_wp, -6.228527908579715E-02_wp, -2.382727643844985E-03_wp, &
   & -5.517827719244788E-03_wp, -2.269316060642431E-02_wp, -5.455464238544920E-03_wp, &
   & -7.142392630901574E-03_wp,  4.733589079461940E-03_wp, -4.772392440508207E-04_wp, &
   & -8.155828201260798E-03_wp,  6.289158309803693E-02_wp, -3.255395128101228E-02_wp, &
   & -1.408760185979790E-01_wp, -7.471553613655965E-02_wp,  3.921484588210891E-02_wp, &
   &  1.160506387841913E-01_wp, -1.122374164166262E-01_wp,  1.345483753142766E-03_wp, &
   &  6.881743823780077E-04_wp,  1.928757151626587E-04_wp,  1.940044737722060E-03_wp, &
   &  3.775733927542388E-03_wp, -2.510457513442781E-02_wp,  2.629728211238974E-02_wp, &
   & -2.862518162679253E-03_wp, -2.509898920227886E-02_wp,  9.109679294113550E-02_wp, &
   &  6.974579433298894E-01_wp,  1.355158928859740E-03_wp, -5.773476009805139E-02_wp, &
   & -5.004618954245851E-02_wp,  9.241924934243366E-02_wp,  1.711383766257452E-01_wp, &
   &  1.399467092466735E-02_wp, -5.949677131429966E-02_wp,  1.291050525765418E-02_wp, &
   & -8.827125108553768E-03_wp, -6.289379435518315E-02_wp, -2.603789280439722E-03_wp, &
   & -8.625490550460051E-04_wp,  8.683327553712171E-03_wp, -5.810647525999135E-03_wp, &
   &  7.212321202108825E-02_wp,  3.599899354309610E-02_wp, -2.877566873736249E-02_wp, &
   & -1.936516509696889E-02_wp, -1.327733209628317E-01_wp, -2.716838217555189E-01_wp, &
   &  3.731613743080805E-01_wp, -2.914932250600488E-01_wp,  6.699758088043848E-03_wp, &
   & -3.083347611846924E-03_wp,  2.075142813212527E-03_wp,  6.121655806829838E-02_wp, &
   & -2.818126004903321E-04_wp, -8.766086966417072E-03_wp, -1.878748571960719E-01_wp, &
   &  2.472472171481012E-02_wp,  1.353447597982711E-02_wp,  2.294197713089968E-02_wp, &
   &  1.871304601260290E-02_wp,  1.479408949786869E-02_wp, -4.076781715761658E-03_wp, &
   &  4.686365515341696E-04_wp, -1.070229159507171E-03_wp, -3.054655989632594E-03_wp, &
   & -3.864036874936527E-02_wp, -1.408816241791581E-01_wp, -3.088581278052090E-01_wp, &
   &  2.325760370067884E-02_wp, -1.688683597876312E-02_wp, -1.122402208931718E-01_wp, &
   & -8.651941146892095E-02_wp,  4.602013044171604E-03_wp,  6.423569787478840E-03_wp, &
   & -1.972490835491246E-03_wp,  3.821405390004228E-03_wp,  3.561560280737851E-03_wp, &
   & -7.385526035989136E-04_wp,  2.803051313766203E-03_wp,  1.524388513770821E-02_wp, &
   &  7.011706174227701E-04_wp, -5.864543153957989E-02_wp,  1.355158928859740E-03_wp, &
   &  7.396015068586778E-01_wp,  3.916392139818056E-02_wp,  7.789176470773139E-03_wp, &
   &  1.711388811988714E-01_wp,  4.308173452428330E-01_wp, -1.622263693389168E-02_wp, &
   &  2.281957939324207E-02_wp, -8.100167831144842E-03_wp,  6.640300786843723E-03_wp, &
   &  2.262282856065698E-02_wp, -1.254378777132380E-03_wp,  5.768955684406493E-03_wp, &
   & -1.628748226245409E-03_wp,  8.543827026121172E-03_wp,  4.031007117573190E-03_wp, &
   & -3.251104180394714E-02_wp, -1.356440831792615E-02_wp, -3.398848896955086E-03_wp, &
   &  8.321663167120150E-02_wp,  3.640443328413350E-01_wp,  3.635082154629348E-01_wp, &
   &  1.490410602015968E-01_wp, -1.649545855874165E-03_wp, -6.032919064557225E-03_wp, &
   & -1.604065601170488E-02_wp, -4.338592338238546E-02_wp, -1.336184314911525E-03_wp, &
   &  9.088616371803117E-03_wp, -2.787259594675615E-01_wp, -1.207719456799169E-01_wp, &
   & -1.943350706591308E-03_wp,  5.398441203417209E-03_wp,  1.461370864808542E-02_wp, &
   & -4.704712498599953E-03_wp, -4.526852586364580E-03_wp, -8.189599121879212E-03_wp, &
   &  3.021112433352000E-03_wp,  1.357272418942471E-02_wp, -4.087198470480553E-02_wp, &
   & -7.471691933004269E-02_wp,  2.324793300914090E-02_wp, -4.756254639370320E-02_wp, &
   &  4.553102888578674E-03_wp, -1.354358795737720E-03_wp, -4.618229847397556E-03_wp, &
   &  7.856583401107661E-03_wp, -3.205742315796407E-03_wp, -3.952082793087780E-03_wp, &
   &  3.647333896716162E-03_wp,  1.755662350481010E-03_wp,  1.717517770848850E-02_wp, &
   & -2.524593097806645E-02_wp, -7.280562139118651E-04_wp,  1.884393418196202E-02_wp, &
   &  1.796960995305845E-03_wp, -5.773476009805139E-02_wp,  3.916392139818056E-02_wp, &
   &  6.469514847578206E-01_wp,  2.659121715620078E-01_wp, -1.400464649369516E-02_wp, &
   &  1.622769772571203E-02_wp, -5.095060355609149E-01_wp,  7.075341151345149E-02_wp, &
   & -1.667256867967071E-02_wp,  1.851969658863674E-03_wp,  6.308156713610985E-02_wp, &
   &  2.533952928739046E-02_wp, -2.947921825878810E-03_wp, -3.456728106606126E-03_wp, &
   & -6.405374650658768E-03_wp,  7.221437239959924E-03_wp, -5.498611974593316E-02_wp, &
   &  1.491645318736809E-02_wp,  2.154920445364152E-01_wp, -4.287043025926325E-02_wp, &
   & -2.163632746459495E-01_wp,  1.000384377202653E-01_wp, -1.847774919476955E-01_wp, &
   &  3.086428269557153E-04_wp,  7.675707547637730E-03_wp,  4.499654897229581E-02_wp, &
   &  1.364786448470975E-02_wp, -4.353406478498622E-03_wp,  5.575618806522659E-03_wp, &
   & -1.220033309745911E-01_wp,  3.222291498982033E-01_wp, -2.525532319575198E-02_wp, &
   & -7.248945552398283E-03_wp,  3.994431557940251E-03_wp,  4.434271568564158E-03_wp, &
   &  1.081588729290596E-02_wp, -2.155701853227797E-03_wp, -1.321253458196261E-02_wp, &
   & -1.979007091270054E-03_wp, -5.726299452013267E-02_wp, -3.922119341459263E-02_wp, &
   &  1.690736399359133E-02_wp, -4.548300829255204E-03_wp, -1.184597738112737E-01_wp, &
   & -6.287618022372808E-02_wp,  3.863301465131361E-02_wp, -4.085784348490161E-02_wp, &
   & -1.229475779570280E-02_wp,  2.532226254985261E-02_wp, -6.423393741742480E-04_wp, &
   &  1.717986995224270E-02_wp, -1.390618957521301E-02_wp, -7.947613217261478E-04_wp, &
   &  6.446404476308520E-03_wp,  3.464331170947334E-03_wp,  9.693675898919289E-02_wp, &
   & -5.004618954245851E-02_wp,  7.789176470773139E-03_wp,  2.659121715620078E-01_wp, &
   &  7.683313973806455E-01_wp,  9.108653280498506E-02_wp, -5.863974522645114E-02_wp, &
   & -1.792736686722264E-03_wp,  7.668033705483055E-03_wp, -4.297938855882854E-03_wp, &
   &  1.162534510098400E-02_wp,  8.018640094390198E-03_wp, -9.862799290665338E-02_wp, &
   &  1.848709597519913E-02_wp, -1.680843299756937E-02_wp,  1.013366365612862E-01_wp, &
   & -2.369640534435551E-02_wp,  1.558104232052441E-01_wp, -8.986945779671203E-02_wp, &
   & -1.817802613867497E-01_wp,  2.878100748007888E-02_wp, -7.067508205037884E-03_wp, &
   & -2.698031152407583E-03_wp, -1.740996340442926E-01_wp,  5.982866467384993E-03_wp, &
   &  6.008621613144669E-03_wp, -1.089461510977031E-01_wp, -1.431166066642352E-02_wp, &
   & -2.126723791504598E-04_wp,  2.764939982309179E-03_wp, -6.228592308352851E-02_wp, &
   &  3.955522415388684E-01_wp, -7.140781347833789E-03_wp,  4.733761954882749E-03_wp, &
   & -4.770325680931546E-04_wp,  8.155084996599015E-03_wp, -2.384170357324099E-03_wp, &
   & -5.518611163450458E-03_wp, -2.268954798545994E-02_wp,  5.456380640725584E-03_wp, &
   &  3.921417238547908E-02_wp,  1.160414267458867E-01_wp, -1.122364450025662E-01_wp, &
   & -1.343524553838545E-03_wp,  6.289372614147358E-02_wp, -3.255395287938539E-02_wp, &
   & -1.408762422576576E-01_wp,  7.471453562286141E-02_wp, -2.510285036163237E-02_wp, &
   &  2.629472700284261E-02_wp, -2.859958918711297E-03_wp,  2.509731414543290E-02_wp, &
   &  6.890824967592879E-04_wp,  1.924630227055869E-04_wp,  1.939924747404366E-03_wp, &
   & -3.775467895130504E-03_wp, -5.005208321108362E-02_wp,  9.241924934243366E-02_wp, &
   &  1.711388811988714E-01_wp, -1.400464649369516E-02_wp,  9.108653280498506E-02_wp, &
   &  6.974503762779322E-01_wp,  1.362571248361449E-03_wp,  5.773652696911138E-02_wp, &
   & -2.602141625735618E-03_wp, -8.617948050810009E-04_wp,  8.681035179764049E-03_wp, &
   &  5.811381672063432E-03_wp, -5.948875302707791E-02_wp,  1.290799380288620E-02_wp, &
   & -8.824902780438812E-03_wp,  6.288644488652945E-02_wp, -1.327678594452733E-01_wp, &
   & -2.716714091968739E-01_wp,  3.731563710172496E-01_wp,  2.914908853714635E-01_wp, &
   &  7.212326273495168E-02_wp,  3.600116055821483E-02_wp, -2.877652838059358E-02_wp, &
   &  1.937138678356154E-02_wp, -3.083173614074472E-03_wp,  6.700480961870179E-03_wp, &
   &  6.121347798932544E-02_wp,  2.075332588360217E-03_wp, -8.766400524998665E-03_wp, &
   & -2.816460996223867E-04_wp,  2.472409049003296E-02_wp, -1.878721411889985E-01_wp, &
   & -4.077573607576199E-03_wp,  4.686288305285728E-04_wp, -1.070308403456159E-03_wp, &
   &  3.054875561499033E-03_wp,  1.353505807534879E-02_wp,  2.294243295262437E-02_wp, &
   &  1.871131041721295E-02_wp, -1.479433721130869E-02_wp, -1.688741882898376E-02_wp, &
   & -1.122392498472225E-01_wp, -8.652690893444472E-02_wp, -4.603471265628236E-03_wp, &
   & -3.864159531118652E-02_wp, -1.408814002950515E-01_wp, -3.088581282385118E-01_wp, &
   & -2.325656785513532E-02_wp, -7.394156897273328E-04_wp,  2.804395208647878E-03_wp, &
   &  1.524260432956563E-02_wp, -7.003288826061063E-04_wp,  6.423166689272697E-03_wp, &
   & -1.972368471019810E-03_wp,  3.821428149507959E-03_wp, -3.561752932146639E-03_wp, &
   &  7.791440829249158E-03_wp,  1.711383766257452E-01_wp,  4.308173452428330E-01_wp, &
   &  1.622769772571203E-02_wp, -5.863974522645114E-02_wp,  1.362571248361450E-03_wp, &
   &  7.396033689595294E-01_wp, -3.916498227122598E-02_wp, -1.254607342720801E-03_wp, &
   &  5.768605983220164E-03_wp, -1.627619035255788E-03_wp, -8.543720900047503E-03_wp, &
   &  2.281645445543275E-02_wp, -8.099021776958312E-03_wp,  6.639140177634144E-03_wp, &
   & -2.261994847862543E-02_wp,  8.321436681991920E-02_wp,  3.640377316426116E-01_wp, &
   &  3.635104767603448E-01_wp, -1.490388921797315E-01_wp,  4.032258295112501E-03_wp, &
   & -3.251171672723583E-02_wp, -1.356394305104294E-02_wp,  3.395346439081716E-03_wp, &
   &  6.032315381140116E-03_wp,  1.649505164758265E-03_wp,  4.338681583913085E-02_wp, &
   &  1.604230038290055E-02_wp, -9.087342946319639E-03_wp,  1.335025177614989E-03_wp, &
   &  1.207716734700453E-01_wp,  2.787181662448878E-01_wp,  4.528041766467476E-03_wp, &
   &  8.189091634808692E-03_wp, -3.022108110398367E-03_wp,  1.357188479253855E-02_wp, &
   &  1.942893511781486E-03_wp, -5.397665312309055E-03_wp, -1.461153704388195E-02_wp, &
   & -4.705408422851873E-03_wp, -4.553650550635379E-03_wp,  1.352397637141819E-03_wp, &
   &  4.619687899269514E-03_wp,  7.858653910191471E-03_wp,  4.086846697200497E-02_wp, &
   &  7.471791946172623E-02_wp, -2.324896861371008E-02_wp, -4.756254716429453E-02_wp, &
   & -1.717362995183904E-02_wp,  2.524346575856478E-02_wp,  7.288003370454037E-04_wp, &
   &  1.884198707659332E-02_wp,  3.205947159947167E-03_wp,  3.953144503594116E-03_wp, &
   & -3.647612696143962E-03_wp,  1.754560513036185E-03_wp, -2.659276333325515E-01_wp, &
   &  1.399467092466735E-02_wp, -1.622263693389169E-02_wp, -5.095060355609149E-01_wp, &
   & -1.792736686722264E-03_wp,  5.773652696911138E-02_wp, -3.916498227122598E-02_wp, &
   &  6.469534635538421E-01_wp, -2.534314551616326E-02_wp,  2.948904725756175E-03_wp, &
   &  3.455221081550927E-03_wp, -6.409682799313698E-03_wp, -7.074313263645753E-02_wp, &
   &  1.667086952269183E-02_wp, -1.850003447867720E-03_wp,  6.307199782805303E-02_wp, &
   &  4.287705717319071E-02_wp,  2.163664561045142E-01_wp, -1.000395884419529E-01_wp, &
   & -1.847769327371833E-01_wp, -7.222533434345535E-03_wp,  5.498718414079060E-02_wp, &
   & -1.491737454609560E-02_wp,  2.154964720833920E-01_wp, -2.018616341381384E-01_wp, &
   &  1.041055766967904E-02_wp,  3.475680838830167E-03_wp, -1.757649574346473E-03_wp, &
   &  9.637106365162448E-03_wp, -2.024238088144413E-01_wp, -1.699130993047991E-03_wp, &
   &  3.496911942846015E-03_wp, -2.551968871121886E-02_wp, -1.310796772062899E-01_wp, &
   & -8.243877639499057E-02_wp,  4.380665652887400E-02_wp,  2.845024133171710E-02_wp, &
   &  7.242096930076826E-02_wp, -3.813587344754937E-03_wp,  5.799744495423170E-03_wp, &
   &  7.619742304759612E-03_wp,  2.590228232093148E-03_wp,  1.291179652397115E-03_wp, &
   & -2.530417551779322E-02_wp, -9.847144852560728E-02_wp,  5.937310100667977E-02_wp, &
   & -2.278812893083784E-02_wp,  7.073249475298077E-02_wp,  2.929510346456513E-02_wp, &
   & -7.215882064800733E-02_wp,  3.618941053847149E-03_wp,  6.534370866653386E-03_wp, &
   & -2.418884425484044E-02_wp,  1.309516163820550E-01_wp,  8.258101038071362E-02_wp, &
   &  4.447228593521083E-02_wp, -9.864018496849172E-02_wp, -5.949677131429966E-02_wp, &
   &  2.281957939324207E-02_wp,  7.075341151345149E-02_wp,  7.668033705483055E-03_wp, &
   & -2.602141625735618E-03_wp, -1.254607342720801E-03_wp, -2.534314551616325E-02_wp, &
   &  1.456572036789029E+00_wp, -8.726255237097133E-04_wp, -3.544596542823982E-05_wp, &
   &  6.699204649868689E-01_wp, -1.508894608281372E-02_wp, -4.253747956213052E-04_wp, &
   & -1.851098518103374E-04_wp, -1.041810958169407E-01_wp, -1.699249422782070E-02_wp, &
   & -6.047954117260467E-07_wp,  1.943524763295232E-05_wp,  1.805611650153299E-02_wp, &
   &  6.582235634799814E-02_wp,  4.623487213329100E-05_wp,  1.377026171748543E-05_wp, &
   & -3.880333508754497E-02_wp, -1.348005264422091E-01_wp,  7.413920473761294E-02_wp, &
   &  1.592836372275454E-03_wp, -8.182952073746532E-03_wp, -7.369382888401302E-02_wp, &
   &  1.349196812143939E-01_wp,  8.535738856329883E-03_wp, -1.567996147383244E-03_wp, &
   &  1.557099978822604E-01_wp, -2.736718176449321E-01_wp, -3.621129115947413E-01_wp, &
   &  2.164122060185590E-01_wp, -7.332213225838278E-03_wp,  3.615078167594412E-02_wp, &
   &  3.259515165671301E-02_wp, -5.536604633567437E-02_wp,  4.223219909973934E-03_wp, &
   & -7.729410252126835E-04_wp,  5.740911507820033E-03_wp, -2.943663228912056E-03_wp, &
   & -1.851123192800146E-02_wp,  1.282979193187429E-02_wp, -8.213618357573981E-03_wp, &
   &  1.668755103203241E-02_wp,  7.434275886919354E-03_wp,  3.627696595311906E-02_wp, &
   &  3.252597635589729E-02_wp,  5.564409057226768E-02_wp, -1.562432385840267E-01_wp, &
   & -2.690441096373638E-01_wp, -3.653721654972368E-01_wp, -2.168684192675906E-01_wp, &
   &  1.848916476504280E-02_wp,  1.291050525765418E-02_wp, -8.100167831144842E-03_wp, &
   & -1.667256867967071E-02_wp, -4.297938855882854E-03_wp, -8.617948050810009E-04_wp, &
   &  5.768605983220164E-03_wp,  2.948904725756175E-03_wp, -8.726255237097133E-04_wp, &
   &  8.285199551995719E-01_wp, -7.150214824724886E-02_wp, -6.156410684641352E-04_wp, &
   & -4.253707314880794E-04_wp, -1.343572025099251E-01_wp,  9.138812138953431E-02_wp, &
   & -6.093398779526438E-04_wp, -7.752754778230711E-05_wp,  8.699073221492716E-03_wp, &
   &  3.607080227520981E-03_wp, -5.829185296378039E-05_wp, -1.468009595234313E-04_wp, &
   &  2.171478713889929E-02_wp, -3.962764959114419E-03_wp,  1.576945435411921E-04_wp, &
   & -5.574264986967414E-02_wp,  2.866332331412033E-02_wp, -4.736727799549166E-03_wp, &
   & -1.752859939197428E-03_wp, -2.841649893839564E-02_wp,  5.590840217582008E-02_wp, &
   &  1.905012656259403E-03_wp,  4.747446067989234E-03_wp,  8.986905531255818E-02_wp, &
   & -3.712663813199765E-01_wp,  3.583325159333834E-01_wp,  9.992231915273242E-02_wp, &
   &  2.419222393111799E-03_wp,  2.884443599340384E-02_wp, -1.363462500872640E-02_wp, &
   & -1.517379452047215E-02_wp, -1.167476765732309E-02_wp,  8.659799042456904E-03_wp, &
   & -1.714098005125333E-03_wp, -3.463083877359813E-03_wp,  1.681962644965352E-02_wp, &
   & -8.756560241211521E-03_wp,  6.821899693253411E-03_wp, -1.878607371143480E-03_wp, &
   & -2.443516507676706E-03_wp,  2.875119908564139E-02_wp, -1.327504017809231E-02_wp, &
   &  1.531544562200643E-02_wp, -9.035542567922990E-02_wp, -3.743568770468897E-01_wp, &
   &  3.680794146991755E-01_wp, -1.004116743940034E-01_wp, -1.681058788278299E-02_wp, &
   & -8.827125108553768E-03_wp,  6.640300786843723E-03_wp,  1.851969658863674E-03_wp, &
   &  1.162534510098400E-02_wp,  8.681035179764049E-03_wp, -1.627619035255788E-03_wp, &
   &  3.455221081550927E-03_wp, -3.544596542823982E-05_wp, -7.150214824724886E-02_wp, &
   &  1.003158080018918E+00_wp, -4.383252872756901E-04_wp, -1.851109241314647E-04_wp, &
   &  9.138802097100168E-02_wp, -2.922802212004755E-01_wp, -2.490858318982935E-04_wp, &
   &  1.441291545281137E-06_wp, -3.583186603423463E-03_wp,  1.299395560340380E-03_wp, &
   & -1.749536999720342E-05_wp, -9.732663513955692E-05_wp,  3.958493355851037E-03_wp, &
   &  8.312175296145786E-03_wp,  7.085480222764424E-05_wp, -1.305267507991367E-01_wp, &
   & -1.049613757799108E-01_wp,  5.294951166937761E-03_wp,  3.253758271617169E-03_wp, &
   & -1.047361603082406E-01_wp, -1.309274647645181E-01_wp,  3.456342577874180E-03_wp, &
   &  5.297970446505873E-03_wp, -1.796270662500462E-01_wp,  2.912542799115598E-01_wp, &
   &  1.489142714551477E-01_wp, -1.860131997232384E-01_wp,  1.739388997798245E-01_wp, &
   & -1.927460744094261E-02_wp,  3.385748079643321E-03_wp,  2.147294740426981E-01_wp, &
   &  7.880149076208948E-03_wp, -5.802103710602418E-03_wp,  8.529977331087443E-03_wp, &
   & -6.333394339464338E-03_wp, -1.011141159736219E-01_wp,  6.280096415787730E-02_wp, &
   & -2.267115234887514E-02_wp,  6.307787433385099E-02_wp,  1.756037405016918E-01_wp, &
   &  1.957479661351896E-02_wp, -3.220333341428738E-03_wp,  2.161460773258060E-01_wp, &
   & -1.833293027304498E-01_wp, -2.931826377618207E-01_wp, -1.498594746533048E-01_wp, &
   & -1.852657794787091E-01_wp, -1.013484377614760E-01_wp, -6.289379435518315E-02_wp, &
   &  2.262282856065699E-02_wp,  6.308156713610985E-02_wp,  8.018640094390198E-03_wp, &
   &  5.811381672063434E-03_wp, -8.543720900047503E-03_wp, -6.409682799313726E-03_wp, &
   &  6.699204649868689E-01_wp, -6.156410684641352E-04_wp, -4.383252872756901E-04_wp, &
   &  1.452755860951080E+00_wp,  1.041805008618904E-01_wp,  6.093382860972202E-04_wp, &
   &  2.490832151561429E-04_wp,  6.898185782588656E-02_wp,  1.850068833011317E-02_wp, &
   & -1.183635594894300E-06_wp,  3.129922706115531E-05_wp,  4.746355917106823E-02_wp, &
   &  3.934915349339888E-02_wp,  5.519002938136542E-05_wp,  3.374468824312119E-05_wp, &
   &  9.318507444294744E-05_wp,  1.041002211525260E-02_wp, -2.018623595147256E-01_wp, &
   & -1.756646422878698E-03_wp,  3.473784430573348E-03_wp, -2.024245396594210E-01_wp, &
   &  9.636569924340338E-03_wp,  3.495014311407543E-03_wp, -1.698159717415534E-03_wp, &
   &  2.845113839599321E-02_wp,  7.242104618418099E-02_wp, -3.813695409354407E-03_wp, &
   & -5.800049750378886E-03_wp, -2.552041991790559E-02_wp, -1.310793590904125E-01_wp, &
   & -8.243836686100103E-02_wp, -4.380651862943320E-02_wp, -9.845927225922103E-02_wp, &
   &  5.936510222427419E-02_wp, -2.278501721523940E-02_wp, -7.072222891782816E-02_wp, &
   &  7.618980820543752E-03_wp,  2.591875878241655E-03_wp,  1.290951888940569E-03_wp, &
   &  2.530055882102593E-02_wp, -2.418957168368160E-02_wp,  1.309512863422543E-01_wp, &
   &  8.258060529965380E-02_wp, -4.447215084892579E-02_wp,  2.929601302874184E-02_wp, &
   & -7.215889877554434E-02_wp,  3.619049329913439E-03_wp, -6.534687037887052E-03_wp, &
   &  7.667271490821401E-03_wp, -2.603789280439722E-03_wp, -1.254378777132380E-03_wp, &
   &  2.533952928739046E-02_wp, -9.862799290665338E-02_wp, -5.948875302707791E-02_wp, &
   &  2.281645445543275E-02_wp, -7.074313263645753E-02_wp, -1.508894608281372E-02_wp, &
   & -4.253707314880794E-04_wp, -1.851109241314647E-04_wp,  1.041805008618904E-01_wp, &
   &  1.456573096809524E+00_wp, -8.726282058784825E-04_wp, -3.544522920170469E-05_wp, &
   & -6.699231995473770E-01_wp,  6.581076592144455E-02_wp,  4.623390617771426E-05_wp, &
   &  1.377154371745163E-05_wp,  3.879861021810282E-02_wp, -1.699238442421719E-02_wp, &
   & -6.023940012383251E-07_wp,  1.943480597946303E-05_wp, -1.805123526243574E-02_wp, &
   &  7.413909982224848E-02_wp, -1.348005538473235E-01_wp, -8.182581808537367E-03_wp, &
   &  1.592236996815013E-03_wp,  1.349196958026813E-01_wp, -7.369372293629420E-02_wp, &
   & -1.567395487485505E-03_wp,  8.535313896577645E-03_wp, -7.332375822996767E-03_wp, &
   &  3.615075610362662E-02_wp,  3.259378511984293E-02_wp,  5.536611255510542E-02_wp, &
   &  1.557102188120999E-01_wp, -2.736720740635467E-01_wp, -3.621132033005950E-01_wp, &
   & -2.164121130792328E-01_wp, -1.850913627906792E-02_wp,  1.282726726055198E-02_wp, &
   & -8.212466815202671E-03_wp, -1.668584366919810E-02_wp,  4.223400132357042E-03_wp, &
   & -7.736922218217283E-04_wp,  5.741259563876600E-03_wp,  2.942679438949034E-03_wp, &
   & -1.562434616964579E-01_wp, -2.690443581569927E-01_wp, -3.653724543429441E-01_wp, &
   &  2.168683175933099E-01_wp,  7.434446578429053E-03_wp,  3.627695352670001E-02_wp, &
   &  3.252458845540963E-02_wp, -5.564415487573214E-02_wp, -4.298131697609277E-03_wp, &
   & -8.625490550460051E-04_wp,  5.768955684406493E-03_wp, -2.947921825878810E-03_wp, &
   &  1.848709597519913E-02_wp,  1.290799380288620E-02_wp, -8.099021776958312E-03_wp, &
   &  1.667086952269183E-02_wp, -4.253747956213052E-04_wp, -1.343572025099251E-01_wp, &
   &  9.138802097100168E-02_wp,  6.093382860972202E-04_wp, -8.726282058784825E-04_wp, &
   &  8.285169171274349E-01_wp, -7.150004766304523E-02_wp,  6.156319550196720E-04_wp, &
   & -1.467970976631559E-04_wp,  2.171627047748486E-02_wp, -3.962973474598682E-03_wp, &
   & -1.576948134907236E-04_wp, -7.752189093408663E-05_wp,  8.698957812115696E-03_wp, &
   &  3.607148540381250E-03_wp,  5.829149479418008E-05_wp,  2.866334386143702E-02_wp, &
   & -5.574271573258759E-02_wp, -1.750276924826151E-03_wp, -4.736230410173405E-03_wp, &
   &  5.590846316495357E-02_wp, -2.841652071952432E-02_wp,  4.746948076908401E-03_wp, &
   &  1.902395363242915E-03_wp,  2.419272907748681E-03_wp,  2.884311599600671E-02_wp, &
   & -1.363227504153741E-02_wp,  1.517387819873661E-02_wp,  8.986925833650755E-02_wp, &
   & -3.712668920587331E-01_wp,  3.583332212693929E-01_wp, -9.992231574193339E-02_wp, &
   &  1.681748597117476E-02_wp, -8.754346717420952E-03_wp,  6.820742093190971E-03_wp, &
   &  1.876645488428163E-03_wp, -1.167579672914539E-02_wp,  8.662091590280946E-03_wp, &
   & -1.715226971518527E-03_wp,  3.464590689409829E-03_wp, -9.035563306592274E-02_wp, &
   & -3.743573887854515E-01_wp,  3.680801231838481E-01_wp,  1.004116692618975E-01_wp, &
   & -2.443563686568657E-03_wp,  2.874985981178474E-02_wp, -1.327264103522864E-02_wp, &
   & -1.531553131254082E-02_wp,  1.162637117966820E-02_wp,  8.683327553712171E-03_wp, &
   & -1.628748226245409E-03_wp, -3.456728106606126E-03_wp, -1.680843299756937E-02_wp, &
   & -8.824902780438812E-03_wp,  6.639140177634144E-03_wp, -1.850003447867720E-03_wp, &
   & -1.851098518103374E-04_wp,  9.138812138953431E-02_wp, -2.922802212004755E-01_wp, &
   &  2.490832151561429E-04_wp, -3.544522920170469E-05_wp, -7.150004766304523E-02_wp, &
   &  1.003152337660315E+00_wp,  4.383147823305224E-04_wp, -9.732282011006750E-05_wp, &
   &  3.958805902420505E-03_wp,  8.312188661074398E-03_wp, -7.085321583330165E-05_wp, &
   &  1.443433154872220E-06_wp, -3.583826023239487E-03_wp,  1.299658324821738E-03_wp, &
   &  1.749630712685038E-05_wp,  1.049620885912693E-01_wp,  1.305276789102963E-01_wp, &
   & -3.250093390684439E-03_wp, -5.293442289541164E-03_wp,  1.309283953591511E-01_wp, &
   &  1.047368757642663E-01_wp, -5.296451585602324E-03_wp, -3.452610340164516E-03_wp, &
   & -1.739403599517548E-01_wp,  1.927387162567067E-02_wp, -3.386000962854876E-03_wp, &
   &  2.147301316072956E-01_wp,  1.796276237353349E-01_wp, -2.912548161985849E-01_wp, &
   & -1.489146544218456E-01_wp, -1.860132438716162E-01_wp,  1.011023440756921E-01_wp, &
   & -6.279365091767158E-02_wp,  2.266829586314910E-02_wp,  6.306833010156396E-02_wp, &
   & -7.880669598994257E-03_wp,  5.801361411158841E-03_wp, -8.530080048747328E-03_wp, &
   & -6.329093098400240E-03_wp,  1.833298769665542E-01_wp,  2.931831809220679E-01_wp, &
   &  1.498598638686168E-01_wp, -1.852658196507851E-01_wp, -1.756052107039139E-01_wp, &
   & -1.957406411145725E-02_wp,  3.220587350774420E-03_wp,  2.161467408179066E-01_wp, &
   & -8.019162958230901E-03_wp, -5.810647525999135E-03_wp,  8.543827026121172E-03_wp, &
   & -6.405374650658741E-03_wp,  1.013366365612862E-01_wp,  6.288644488652945E-02_wp, &
   & -2.261994847862543E-02_wp,  6.307199782805303E-02_wp, -1.041810958169407E-01_wp, &
   & -6.093398779526438E-04_wp, -2.490858318982935E-04_wp,  6.898185782588662E-02_wp, &
   & -6.699231995473770E-01_wp,  6.156319550196720E-04_wp,  4.383147823305224E-04_wp, &
   &  1.452756989081952E+00_wp, -3.934016995923928E-02_wp, -5.519275426865360E-05_wp, &
   & -3.374646025114020E-05_wp,  9.670663156857717E-05_wp, -1.849931819948261E-02_wp, &
   &  1.180397284356355E-06_wp, -3.129870738785489E-05_wp,  4.745860498902796E-02_wp, &
   &  3.468787035654369E-03_wp, -1.852728792443303E-03_wp, -2.025727736049031E-01_wp, &
   &  9.226336221603473E-03_wp, -1.826507296104563E-03_wp,  3.530348390705637E-03_wp, &
   &  9.228907270398379E-03_wp, -2.025736978208567E-01_wp,  7.804943062986344E-03_wp, &
   & -2.538903847321168E-03_wp,  1.300248247308219E-03_wp, -2.533085154898232E-02_wp, &
   & -9.863723434119841E-02_wp, -5.981538476871626E-02_wp, -2.301012739452538E-02_wp, &
   &  7.088307240209220E-02_wp, -2.369189222685061E-02_wp,  1.327559948026534E-01_wp, &
   & -8.322877082234913E-02_wp,  4.287940936913996E-02_wp,  2.878693085465274E-02_wp, &
   & -7.213721278339691E-02_wp, -4.038779749340381E-03_wp,  7.215211206464095E-03_wp, &
   & -9.884426245318094E-02_wp,  5.975943572505005E-02_wp,  2.288650062627778E-02_wp, &
   &  7.107247800333402E-02_wp,  7.726829452574935E-03_wp,  2.589483077764267E-03_wp, &
   & -1.244094321036355E-03_wp, -2.527544829182120E-02_wp,  2.878275285056010E-02_wp, &
   &  7.212321202108825E-02_wp,  4.031007117573194E-03_wp,  7.221437239959924E-03_wp, &
   & -2.369640534435551E-02_wp, -1.327678594452733E-01_wp,  8.321436681991920E-02_wp, &
   &  4.287705717319069E-02_wp, -1.699249422782070E-02_wp, -7.752754778230711E-05_wp, &
   &  1.441291545281137E-06_wp,  1.850068833011317E-02_wp,  6.581076592144458E-02_wp, &
   & -1.467970976631559E-04_wp, -9.732282011006750E-05_wp, -3.934016995923928E-02_wp, &
   &  1.452035474479663E+00_wp,  8.005824500055966E-06_wp,  3.307496337576725E-06_wp, &
   &  6.702493915081253E-01_wp, -1.502911934037832E-02_wp,  2.563373224162916E-06_wp, &
   &  8.779585124951143E-06_wp, -1.044178478991612E-01_wp, -1.556266923650076E-03_wp, &
   &  8.191399118302628E-03_wp,  1.363116666354803E-01_wp, -7.303191575335419E-02_wp, &
   & -8.232201340172253E-03_wp,  1.580861854842841E-03_wp,  7.303210790879353E-02_wp, &
   & -1.363325014206645E-01_wp, -4.325931302248842E-03_wp, -8.372682940089578E-04_wp, &
   & -5.641419553982402E-03_wp,  2.977890074213256E-03_wp,  1.864332326969586E-02_wp, &
   &  1.301553465986978E-02_wp,  8.240915671879087E-03_wp, -1.683846527489821E-02_wp, &
   & -1.558062066062249E-01_wp, -2.716716161450161E-01_wp,  3.640456085580128E-01_wp, &
   & -2.163725859997901E-01_wp,  7.064998248426048E-03_wp,  3.598800667611273E-02_wp, &
   & -3.253315639511126E-02_wp,  5.498834317754229E-02_wp, -1.863273215096028E-02_wp, &
   &  1.297006349097501E-02_wp,  8.306702077661978E-03_wp,  1.689252113263069E-02_wp, &
   &  4.281002320965354E-03_wp, -7.754048407027505E-04_wp, -5.723459428157306E-03_wp, &
   & -2.965209619472379E-03_wp, -7.069186144772829E-03_wp,  3.599899354309610E-02_wp, &
   & -3.251104180394714E-02_wp, -5.498611974593316E-02_wp,  1.558104232052441E-01_wp, &
   & -2.716714091968739E-01_wp,  3.640377316426116E-01_wp,  2.163664561045142E-01_wp, &
   & -6.047954117260462E-07_wp,  8.699073221492716E-03_wp, -3.583186603423463E-03_wp, &
   & -1.183635594894300E-06_wp,  4.623390617771426E-05_wp,  2.171627047748486E-02_wp, &
   &  3.958805902420505E-03_wp, -5.519275426865360E-05_wp,  8.005824500055966E-06_wp, &
   &  8.275253657015017E-01_wp,  7.250286969534346E-02_wp,  7.680658951093900E-07_wp, &
   &  2.562851309669045E-06_wp, -1.339197102470138E-01_wp, -9.021164454194042E-02_wp, &
   & -2.243100947738236E-06_wp, -4.697567182304204E-03_wp, -2.024632422877527E-03_wp, &
   & -5.649549802151434E-02_wp,  2.808355212843569E-02_wp,  1.976190866293467E-03_wp, &
   &  4.713874274092306E-03_wp, -2.808387708646047E-02_wp,  5.647403804236660E-02_wp, &
   & -1.164661813543127E-02_wp, -8.593498681167266E-03_wp, -1.665767637263845E-03_wp, &
   & -3.429125807313961E-03_wp,  1.685365810430683E-02_wp,  8.956889561815496E-03_wp, &
   &  6.844292821216313E-03_wp, -2.034015039736788E-03_wp,  8.988384688691844E-02_wp, &
   &  3.731497359428295E-01_wp,  3.635064035974106E-01_wp,  1.000300771004840E-01_wp, &
   &  2.683597700459387E-03_wp, -2.878459476131501E-02_wp, -1.360180780112587E-02_wp, &
   & -1.489667751705982E-02_wp, -1.681206142369586E-02_wp,  8.960267273041791E-03_wp, &
   &  6.864908474993436E-03_wp,  2.089735104077940E-03_wp,  1.158987035355509E-02_wp, &
   & -8.590353289713753E-03_wp, -1.671899435629965E-03_wp,  3.421469452061979E-03_wp, &
   & -2.697581873754085E-03_wp, -2.877566873736249E-02_wp, -1.356440831792615E-02_wp, &
   &  1.491645318736809E-02_wp, -8.986945779671203E-02_wp,  3.731563710172496E-01_wp, &
   &  3.635104767603448E-01_wp, -1.000395884419529E-01_wp,  1.943524763295232E-05_wp, &
   &  3.607080227520981E-03_wp,  1.299395560340380E-03_wp,  3.129922706115531E-05_wp, &
   &  1.377154371745163E-05_wp, -3.962973474598682E-03_wp,  8.312188661074398E-03_wp, &
   & -3.374646025114020E-05_wp,  3.307496337576725E-06_wp,  7.250286969534346E-02_wp, &
   &  1.004321753334088E+00_wp,  1.734075911784468E-05_wp,  8.780440738912821E-06_wp, &
   & -9.021154108160916E-02_wp, -2.896920308217301E-01_wp, -7.604913427095871E-06_wp, &
   &  5.252068851896829E-03_wp,  3.770544286479827E-03_wp, -1.294885328690167E-01_wp, &
   & -1.037087748751767E-01_wp,  3.755665465210820E-03_wp,  5.331831361425191E-03_wp, &
   & -1.037016869246249E-01_wp, -1.294856726660747E-01_wp,  8.205582248790529E-03_wp, &
   &  5.747215247280658E-03_wp,  8.506711270493468E-03_wp, -6.449076849290553E-03_wp, &
   & -1.014984336142544E-01_wp, -6.321907174909422E-02_wp, -2.284112114653753E-02_wp, &
   &  6.327477665737241E-02_wp, -1.817969489626665E-01_wp, -2.914955018238549E-01_wp, &
   &  1.490313042946019E-01_wp, -1.847694129517921E-01_wp,  1.741398116942549E-01_wp, &
   &  1.934151879091920E-02_wp,  3.382014459363425E-03_wp,  2.154832128200224E-01_wp, &
   & -1.015205171001342E-01_wp,  6.330195250473385E-02_wp,  2.284663513307648E-02_wp, &
   &  6.345913253029471E-02_wp,  7.917211034894435E-03_wp, -5.795216302978646E-03_wp, &
   & -8.473940779935293E-03_wp, -6.355489319584323E-03_wp,  1.741109383746041E-01_wp, &
   & -1.936516509696889E-02_wp, -3.398848896955086E-03_wp,  2.154920445364152E-01_wp, &
   & -1.817802613867497E-01_wp,  2.914908853714635E-01_wp, -1.490388921797315E-01_wp, &
   & -1.847769327371833E-01_wp,  1.805611650153296E-02_wp, -5.829185296378038E-05_wp, &
   & -1.749536999720342E-05_wp,  4.746355917106823E-02_wp,  3.879861021810280E-02_wp, &
   & -1.576948134907236E-04_wp, -7.085321583330165E-05_wp,  9.670663156857717E-05_wp, &
   &  6.702493915081253E-01_wp,  7.680658951093892E-07_wp,  1.734075911784468E-05_wp, &
   &  1.454209265605659E+00_wp,  1.044159377731140E-01_wp,  2.243144933113974E-06_wp, &
   &  7.603611023070860E-06_wp,  6.855067163425388E-02_wp, -1.851743395625369E-03_wp, &
   &  3.468244026829145E-03_wp,  9.232026474174537E-03_wp, -2.025836973893619E-01_wp, &
   &  3.529801573999167E-03_wp, -1.825513260223371E-03_wp, -2.025846191385778E-01_wp, &
   &  9.234596091040548E-03_wp, -9.863558958445527E-02_wp, -5.981308331902150E-02_wp, &
   & -2.300906688145581E-02_wp, -7.088130909614150E-02_wp,  7.805716983696948E-03_wp, &
   & -2.539141675934529E-03_wp,  1.300202442455607E-03_wp,  2.533132611978335E-02_wp, &
   &  2.878518138124739E-02_wp, -7.213726217635565E-02_wp, -4.040031653327232E-03_wp, &
   & -7.216306434360636E-03_wp, -2.370309164556275E-02_wp,  1.327614522838112E-01_wp, &
   & -8.323103634066595E-02_wp, -4.287278436401225E-02_wp,  7.727620842577441E-03_wp, &
   &  2.589722239178298E-03_wp, -1.244046434831659E-03_wp,  2.527593529481102E-02_wp, &
   & -9.884265042117255E-02_wp,  5.975712114865252E-02_wp,  2.288543553170518E-02_wp, &
   & -7.107072288070147E-02_wp, -2.370761223236972E-02_wp, -1.327733209628317E-01_wp, &
   &  8.321663167120150E-02_wp, -4.287043025926325E-02_wp,  2.878100748007888E-02_wp, &
   &  7.212326273495168E-02_wp,  4.032258295112501E-03_wp, -7.222533434345535E-03_wp, &
   &  6.582235634799817E-02_wp, -1.468009595234313E-04_wp, -9.732663513955692E-05_wp, &
   &  3.934915349339888E-02_wp, -1.699238442421716E-02_wp, -7.752189093408663E-05_wp, &
   &  1.443433154872220E-06_wp, -1.849931819948264E-02_wp, -1.502911934037832E-02_wp, &
   &  2.562851309669046E-06_wp,  8.780440738912821E-06_wp,  1.044159377731140E-01_wp, &
   &  1.452052682128937E+00_wp,  8.007481338732793E-06_wp,  3.307038913590702E-06_wp, &
   & -6.702534475659446E-01_wp,  8.190947240366564E-03_wp, -1.556663367068126E-03_wp, &
   & -7.303372313506439E-02_wp,  1.363173580525894E-01_wp,  1.581255206164669E-03_wp, &
   & -8.231748457816464E-03_wp, -1.363381950886083E-01_wp,  7.303391508409832E-02_wp, &
   &  1.864218687887590E-02_wp,  1.301470143766196E-02_wp,  8.242867184823468E-03_wp, &
   &  1.683742740209905E-02_wp, -4.325915303674373E-03_wp, -8.379272217330506E-04_wp, &
   & -5.640644943804287E-03_wp, -2.977223798203978E-03_wp,  7.063319593985584E-03_wp, &
   &  3.599017371484808E-02_wp, -3.253383232104642E-02_wp, -5.498940836255696E-02_wp, &
   & -1.557979903549532E-01_wp, -2.716840289418065E-01_wp,  3.640522108596773E-01_wp, &
   &  2.163694037783475E-01_wp,  4.280987002069590E-03_wp, -7.760709388883756E-04_wp, &
   & -5.722674596992679E-03_wp,  2.964544962364357E-03_wp, -1.863159696218155E-02_wp, &
   &  1.296921634983868E-02_wp,  8.308680173787811E-03_wp, -1.689148415298210E-02_wp, &
   &  1.558022075373159E-01_wp, -2.716838217555189E-01_wp,  3.640443328413350E-01_wp, &
   & -2.163632746459495E-01_wp, -7.067508205037884E-03_wp,  3.600116055821483E-02_wp, &
   & -3.251171672723583E-02_wp,  5.498718414079060E-02_wp,  4.623487213329100E-05_wp, &
   &  2.171478713889929E-02_wp,  3.958493355851037E-03_wp,  5.519002938136542E-05_wp, &
   & -6.023940012383255E-07_wp,  8.698957812115696E-03_wp, -3.583826023239487E-03_wp, &
   &  1.180397284356355E-06_wp,  2.563373224162915E-06_wp, -1.339197102470138E-01_wp, &
   & -9.021154108160916E-02_wp,  2.243144933113974E-06_wp,  8.007481338732793E-06_wp, &
   &  8.275371387806604E-01_wp,  7.249821911786589E-02_wp, -7.685815583641724E-07_wp, &
   & -2.024491870468237E-03_wp, -4.697511095964893E-03_wp,  2.808423687111002E-02_wp, &
   & -5.649750140866367E-02_wp,  4.713819655831289E-03_wp,  1.976050064171017E-03_wp, &
   &  5.647603482094272E-02_wp, -2.808456425118116E-02_wp,  1.685404835755008E-02_wp, &
   &  8.957318822404307E-03_wp,  6.843528572579309E-03_wp,  2.034317339994986E-03_wp, &
   & -1.164665903359848E-02_wp, -8.593342563929714E-03_wp, -1.666141520371214E-03_wp, &
   &  3.428901092250826E-03_wp,  2.684049542240108E-03_wp, -2.878545675566218E-02_wp, &
   & -1.360134300980547E-02_wp,  1.489760162623975E-02_wp,  8.988113142672710E-02_wp, &
   &  3.731547362524139E-01_wp,  3.635041436227353E-01_wp, -1.000289298114973E-01_wp, &
   &  1.158991091617597E-02_wp, -8.590194921635699E-03_wp, -1.672278219204951E-03_wp, &
   & -3.421246430889103E-03_wp, -1.681244984664289E-02_wp,  8.960702273653885E-03_wp, &
   &  6.864133072952403E-03_wp, -2.090036052235476E-03_wp, -8.986673754216588E-02_wp, &
   &  3.731613743080805E-01_wp,  3.635082154629348E-01_wp,  1.000384377202653E-01_wp, &
   & -2.698031152407583E-03_wp, -2.877652838059358E-02_wp, -1.356394305104294E-02_wp, &
   & -1.491737454609560E-02_wp,  1.377026171748543E-05_wp, -3.962764959114419E-03_wp, &
   &  8.312175296145786E-03_wp,  3.374468824312119E-05_wp,  1.943480597946303E-05_wp, &
   &  3.607148540381250E-03_wp,  1.299658324821738E-03_wp, -3.129870738785489E-05_wp, &
   &  8.779585124951143E-06_wp, -9.021164454194042E-02_wp, -2.896920308217301E-01_wp, &
   &  7.603611023070860E-06_wp,  3.307038913590702E-06_wp,  7.249821911786589E-02_wp, &
   &  1.004322548291383E+00_wp, -1.734031483932472E-05_wp, -3.771975464662083E-03_wp, &
   & -5.251923220272459E-03_wp,  1.037044247528985E-01_wp,  1.294895097260828E-01_wp, &
   & -5.331683711214077E-03_wp, -3.757103059626676E-03_wp,  1.294866496059647E-01_wp, &
   &  1.036973409355553E-01_wp,  1.014979472509048E-01_wp,  6.321787876977525E-02_wp, &
   &  2.284068750337667E-02_wp,  6.327359645025230E-02_wp, -8.205281927793753E-03_wp, &
   & -5.746306227340856E-03_wp, -8.506322112167262E-03_wp, -6.449877852665542E-03_wp, &
   & -1.741285000906176E-01_wp, -1.934774622854778E-02_wp, -3.378510275726217E-03_wp, &
   &  2.154876430367316E-01_wp,  1.817932770756827E-01_wp,  2.914978418533170E-01_wp, &
   & -1.490334705749429E-01_wp, -1.847699740988629E-01_wp, -7.916934585009305E-03_wp, &
   &  5.794296157964703E-03_wp,  8.473545712606514E-03_wp, -6.356288505501823E-03_wp, &
   &  1.015200503060320E-01_wp, -6.330075475231250E-02_wp, -2.284620243063615E-02_wp, &
   &  6.345794212690893E-02_wp,  1.817765936749479E-01_wp, -2.914932250600488E-01_wp, &
   &  1.490410602015968E-01_wp, -1.847774919476955E-01_wp, -1.740996340442926E-01_wp, &
   &  1.937138678356154E-02_wp,  3.395346439081719E-03_wp,  2.154964720833920E-01_wp, &
   & -3.880333508754497E-02_wp,  1.576945435411920E-04_wp,  7.085480222764424E-05_wp, &
   &  9.318507444294744E-05_wp, -1.805123526243571E-02_wp,  5.829149479418008E-05_wp, &
   &  1.749630712685040E-05_wp,  4.745860498902790E-02_wp, -1.044178478991612E-01_wp, &
   & -2.243100947738235E-06_wp, -7.604913427095871E-06_wp,  6.855067163425388E-02_wp, &
   & -6.702534475659446E-01_wp, -7.685815583641716E-07_wp, -1.734031483932472E-05_wp, &
   &  1.454213654512245E+00_wp], &
   &  shape(density))

   call get_structure(mol, "X23", "pyrazine")

   call test_num_op_grad(error, mol, density, make_exchange_gxtb, thr_in=thr2)

end subroutine test_op_g_fock_pyrazine


!> Test analytic vs numerical gradient of exchange energy for oxacb (PBC, op grad)
subroutine test_op_g_fock_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(52,52,1) = reshape([ &
   &  3.168846521055535E-01_wp,  3.250755204046418E-03_wp,  4.690629305014816E-02_wp, &
   &  1.319149483089800E-03_wp, -1.196224587481389E-02_wp,  3.242940457173099E-02_wp, &
   & -5.871962715475485E-02_wp, -9.800345301311790E-02_wp, -6.279131220983131E-03_wp, &
   & -4.601741626533677E-04_wp,  7.192266070701130E-04_wp,  6.857108714237196E-04_wp, &
   &  7.434613310779223E-02_wp,  1.286760939909034E-02_wp, -1.631890529088685E-02_wp, &
   &  1.298721004409743E-01_wp, -6.679887117572364E-03_wp,  1.072389224409596E-03_wp, &
   & -1.355789093345902E-03_wp,  2.361020252331918E-03_wp, -3.832625659357759E-02_wp, &
   &  4.198805679545000E-02_wp, -6.413881380697377E-02_wp,  5.753380864868125E-02_wp, &
   & -3.040323340029377E-04_wp, -4.493523318597993E-04_wp,  8.203910821059849E-03_wp, &
   & -8.552249264904208E-04_wp,  4.911426046033805E-02_wp,  8.458327097158366E-03_wp, &
   & -7.555890109641587E-03_wp, -1.550649509124217E-01_wp,  1.364518771015944E-04_wp, &
   &  1.031567456603362E-02_wp, -4.661771912720552E-03_wp,  1.352730910145476E-03_wp, &
   &  7.580502462269115E-02_wp,  3.423914070860962E-02_wp, -4.526923976551495E-02_wp, &
   &  5.290300153629187E-01_wp,  1.454694527520587E-03_wp, -5.230543099499541E-03_wp, &
   & -5.475539339841637E-03_wp,  5.029086448775175E-03_wp, -5.873750796922662E-02_wp, &
   &  2.997882420757660E-02_wp, -5.591177901303645E-02_wp, -2.873784305744250E-02_wp, &
   &  5.401793605232210E-03_wp,  1.165434529800888E-02_wp, -4.088414634391401E-03_wp, &
   & -5.238704930788645E-03_wp,  3.250755204046418E-03_wp,  3.170972776290459E-01_wp, &
   &  1.258457034745813E-03_wp,  4.712529740150160E-02_wp, -6.299521649534997E-03_wp, &
   & -4.148464495008240E-04_wp, -6.908893418827178E-04_wp, -6.310477581421732E-04_wp, &
   & -1.168795641010645E-02_wp,  3.239728633724738E-02_wp,  5.872652032400529E-02_wp, &
   &  9.837219950312542E-02_wp, -6.703513719039950E-03_wp,  1.021237784893391E-03_wp, &
   &  1.303333505416583E-03_wp, -2.390897664126387E-03_wp,  7.444902492001537E-02_wp, &
   &  1.281459086295719E-02_wp,  1.625674221823976E-02_wp, -1.302057363170864E-01_wp, &
   & -2.910217311767796E-04_wp, -4.251158967129966E-04_wp, -8.169084722126843E-03_wp, &
   &  8.450451233251027E-04_wp, -3.851718429301771E-02_wp,  4.207907422154082E-02_wp, &
   &  6.421113390698539E-02_wp, -5.744332214139890E-02_wp,  1.605240717571781E-04_wp, &
   &  1.027345385064401E-02_wp,  4.621800031846540E-03_wp, -1.358714909350574E-03_wp, &
   &  4.913915356091953E-02_wp,  8.212285660166099E-03_wp,  7.164192470386322E-03_wp, &
   &  1.548485620523062E-01_wp,  1.437809962620081E-03_wp, -5.220645740682702E-03_wp, &
   &  5.464860841707854E-03_wp, -5.022889917614801E-03_wp,  7.579206001018166E-02_wp, &
   &  3.427650757146209E-02_wp,  4.529702101290851E-02_wp, -5.291873073042418E-01_wp, &
   &  5.398541212814784E-03_wp,  1.170724185580760E-02_wp,  4.068259828654927E-03_wp, &
   &  5.306119173344457E-03_wp, -5.866460440253323E-02_wp,  2.995383714867841E-02_wp, &
   &  5.588706921176288E-02_wp,  2.859858066898498E-02_wp,  4.690629305014816E-02_wp, &
   &  1.258457034745813E-03_wp,  3.168959771888755E-01_wp,  3.272390413607148E-03_wp, &
   &  7.433829043550980E-02_wp, -1.286423908354723E-02_wp,  1.632162713121671E-02_wp, &
   & -1.298589941021538E-01_wp, -6.661092935395364E-03_wp, -1.044545316910896E-03_wp, &
   &  1.328230984236084E-03_wp, -2.359384643701530E-03_wp, -1.196966000154187E-02_wp, &
   & -3.242534561842155E-02_wp,  5.872483092002760E-02_wp,  9.798681880857191E-02_wp, &
   & -6.251871088207496E-03_wp,  4.174789300134936E-04_wp, -7.021794444146690E-04_wp, &
   & -6.985629288825479E-04_wp,  4.911014170955832E-02_wp, -8.470644344694443E-03_wp, &
   &  7.545304111834568E-03_wp,  1.550495190069728E-01_wp,  1.132634881668910E-04_wp, &
   & -1.029895438504152E-02_wp,  4.659438803525548E-03_wp, -1.339644572136041E-03_wp, &
   & -3.832627968016269E-02_wp, -4.199367163358746E-02_wp,  6.411417849114696E-02_wp, &
   & -5.752151081594777E-02_wp, -2.969041615684436E-04_wp,  4.515508073346490E-04_wp, &
   & -8.192072030219997E-03_wp,  8.878425530952222E-04_wp, -5.873365961043855E-02_wp, &
   & -2.997885454140316E-02_wp,  5.589274797182139E-02_wp,  2.873288137202622E-02_wp, &
   &  5.408112689239400E-03_wp, -1.169935146473086E-02_wp,  4.068025141748031E-03_wp, &
   &  5.322386258012152E-03_wp,  7.579801364543501E-02_wp, -3.425062092156357E-02_wp, &
   &  4.526598713611546E-02_wp, -5.290432544172488E-01_wp,  1.442174296498074E-03_wp, &
   &  5.218762433237576E-03_wp,  5.494918943424746E-03_wp, -4.977729138402987E-03_wp, &
   &  1.319149483089800E-03_wp,  4.712529740150160E-02_wp,  3.272390413607148E-03_wp, &
   &  3.167239546732092E-01_wp, -6.736872707004743E-03_wp, -1.016447668114869E-03_wp, &
   & -1.374585876305031E-03_wp,  2.362501271069201E-03_wp,  7.093392456161390E-02_wp, &
   & -1.209850983039165E-02_wp, -1.515898622269251E-02_wp,  1.257751136631665E-01_wp, &
   & -6.250505483431419E-03_wp,  5.026181659209923E-04_wp,  7.919192080460943E-04_wp, &
   &  6.612697887665494E-04_wp, -1.272759551324487E-02_wp, -3.216708991641609E-02_wp, &
   & -5.856052840859258E-02_wp, -9.623347293437912E-02_wp,  2.078487962056577E-04_wp, &
   & -1.030952412861568E-02_wp, -4.693807740685106E-03_wp,  1.293707498310058E-03_wp, &
   &  5.203536420446843E-02_wp, -9.234309186277682E-03_wp, -9.221765000020827E-03_wp, &
   & -1.560169469960605E-01_wp, -4.101483768953129E-04_wp,  9.553232181870224E-04_wp, &
   &  8.567822520753905E-03_wp, -6.397531873150973E-04_wp, -3.794187581816817E-02_wp, &
   & -4.210580981934167E-02_wp, -6.422550781484228E-02_wp,  5.660925794823445E-02_wp, &
   &  5.475717228827796E-03_wp, -1.180731862149381E-02_wp, -4.014856647449936E-03_wp, &
   & -5.311396711876894E-03_wp, -5.839639244371739E-02_wp, -3.020781403039341E-02_wp, &
   & -5.626759640739265E-02_wp, -2.956141066515159E-02_wp,  1.510806635622226E-03_wp, &
   &  4.971150257396431E-03_wp, -5.627628235781873E-03_wp,  5.089000567812710E-03_wp, &
   &  7.647175198668643E-02_wp, -3.340646505662848E-02_wp, -4.523805977454380E-02_wp, &
   &  5.286168946426076E-01_wp, -1.196224587481389E-02_wp, -6.299521649534996E-03_wp, &
   &  7.433829043550980E-02_wp, -6.736872707004743E-03_wp,  8.429326196000994E-01_wp, &
   &  5.454186006526193E-03_wp, -1.317394141270070E-02_wp, -1.380164659261662E-01_wp, &
   & -1.039840771817273E-03_wp,  9.441813035990095E-04_wp, -7.025956212604888E-04_wp, &
   & -4.875661919418878E-04_wp,  1.891309702292605E-01_wp,  8.367566084559982E-03_wp, &
   & -6.144229616403341E-03_wp,  2.753621765898573E-01_wp, -7.895148297830563E-04_wp, &
   &  2.255651556425698E-03_wp, -4.236227136247043E-03_wp,  1.053909754895402E-03_wp, &
   & -1.657257723234933E-01_wp, -1.708980155892887E-01_wp,  2.642395600556753E-01_wp, &
   & -1.789389989696000E-01_wp, -2.318060720164255E-02_wp, -1.544549398688108E-02_wp, &
   & -1.220172254632664E-02_wp, -1.000250708768312E-02_wp, -8.138703378017074E-02_wp, &
   &  1.534008046510517E-02_wp, -3.036148229123299E-02_wp, -3.127961791029916E-01_wp, &
   & -2.595677925547819E-02_wp,  1.313400746515209E-02_wp,  1.774030040201211E-02_wp, &
   &  2.455335723719948E-03_wp, -1.321201998006772E-01_wp,  1.488951503221468E-01_wp, &
   & -2.510248276258650E-01_wp, -1.534607830666470E-01_wp, -1.460308024330004E-02_wp, &
   &  2.145605836418894E-03_wp,  1.982681150420149E-03_wp,  2.367182058618314E-04_wp, &
   & -1.491414911419830E-01_wp,  2.672404872968473E-02_wp, -6.106109832140092E-02_wp, &
   & -2.170756894938684E-01_wp, -1.264404032600393E-02_wp, -5.711514899287208E-03_wp, &
   & -3.087023996236248E-03_wp,  2.538143068475399E-03_wp,  3.242940457173099E-02_wp, &
   & -4.148464495008240E-04_wp, -1.286423908354723E-02_wp, -1.016447668114869E-03_wp, &
   &  5.454186006526193E-03_wp,  4.900975988926149E-01_wp, -5.586784445635188E-03_wp, &
   &  4.986232267582763E-03_wp,  9.499478900720762E-04_wp,  5.167756375487324E-03_wp, &
   & -3.588284470306161E-03_wp, -3.182467059190740E-03_wp, -8.349427620557459E-03_wp, &
   &  3.010369337275243E-02_wp,  3.143850179021387E-02_wp, -1.188768059898162E-02_wp, &
   & -2.352905761865181E-03_wp,  8.293513906870725E-03_wp, -5.786828146520061E-03_wp, &
   &  8.788780580587822E-04_wp, -1.803512460547285E-02_wp,  2.893336223420668E-01_wp, &
   &  4.456579031710494E-01_wp, -1.957747719980853E-01_wp,  3.297270138654383E-03_wp, &
   &  3.233388733204503E-03_wp,  1.591282222122241E-03_wp, -2.811112246012069E-03_wp, &
   &  7.217962268249751E-03_wp, -1.465429531912288E-02_wp, -1.038968907927841E-02_wp, &
   &  3.264253319671220E-02_wp, -6.197049531325822E-03_wp,  3.576015982288077E-03_wp, &
   &  5.300441441243784E-03_wp, -9.017402911145283E-03_wp, -2.412750523444410E-03_wp, &
   &  1.499577451346859E-01_wp,  3.926320642675794E-01_wp,  1.033017651324187E-01_wp, &
   & -1.256618474760451E-02_wp,  1.419760878664407E-02_wp, -1.513744763688585E-03_wp, &
   &  8.101542765591399E-03_wp,  1.767956891538433E-03_wp, -1.683065106280155E-02_wp, &
   & -8.700242824565430E-03_wp, -4.025253508741049E-02_wp,  1.000643989248778E-02_wp, &
   &  1.982971406572112E-02_wp, -1.893954077802771E-03_wp, -3.014670166330845E-03_wp, &
   & -5.871962715475485E-02_wp, -6.908893418827178E-04_wp,  1.632162713121671E-02_wp, &
   & -1.374585876305031E-03_wp, -1.317394141270070E-02_wp, -5.586784445635188E-03_wp, &
   &  4.961051989029899E-01_wp, -4.159445251992173E-03_wp,  6.940628533841650E-04_wp, &
   &  3.581781097085035E-03_wp, -1.759937428925466E-03_wp, -6.500358384336613E-03_wp, &
   &  6.156068262093434E-03_wp,  3.143822674633701E-02_wp, -2.005312615417137E-03_wp, &
   &  4.403932312142212E-03_wp, -4.302988406396988E-03_wp,  5.876240630642964E-03_wp, &
   & -2.416968585163443E-03_wp,  1.120519030515370E-03_wp,  2.702426709506782E-02_wp, &
   &  4.451127413192278E-01_wp, -1.070069421144721E-01_wp,  2.849107342135790E-01_wp, &
   &  8.499489701407577E-03_wp,  1.023902985142095E-02_wp, -9.798591834190845E-03_wp, &
   & -4.745814818565420E-03_wp,  7.517433013925840E-03_wp, -1.738383260717502E-02_wp, &
   &  1.319690499959618E-02_wp, -4.064972177795283E-02_wp, -6.635914471490414E-03_wp, &
   &  1.637587219690404E-02_wp, -1.155542345146213E-02_wp,  3.226255954096729E-04_wp, &
   &  1.066750419360065E-02_wp,  3.907479332333729E-01_wp, -2.335759858806779E-01_wp, &
   & -1.745947607698172E-01_wp, -6.462818904493103E-03_wp,  1.281852403078568E-02_wp, &
   & -4.627416287810797E-04_wp,  3.411384254096665E-03_wp, -8.267152240045386E-03_wp, &
   & -2.559729433992322E-02_wp,  5.062772494027674E-03_wp,  6.058348225014799E-02_wp, &
   &  7.642773678844308E-03_wp,  2.531309065105236E-02_wp, -7.459720442475663E-03_wp, &
   &  4.013830365771267E-03_wp, -9.800345301311790E-02_wp, -6.310477581421741E-04_wp, &
   & -1.298589941021538E-01_wp,  2.362501271069201E-03_wp, -1.380164659261662E-01_wp, &
   &  4.986232267582763E-03_wp, -4.159445251992173E-03_wp,  6.337177280410977E-01_wp, &
   &  4.893621576621127E-04_wp,  3.207817638247836E-03_wp, -6.516583487466641E-03_wp, &
   &  1.251090620324090E-02_wp, -2.753678417295266E-01_wp, -1.190960129238832E-02_wp, &
   &  4.394774512976690E-03_wp, -4.023182634145880E-01_wp,  1.055441718594615E-03_wp, &
   & -9.481975499901815E-04_wp,  1.155982015342645E-03_wp, -1.027700775080411E-02_wp, &
   &  4.241317133803125E-02_wp, -2.119840694276171E-01_wp,  3.057780897872059E-01_wp, &
   & -1.438377687257794E-02_wp,  2.319570021325564E-02_wp,  1.049207368818320E-02_wp, &
   &  5.002988700501795E-03_wp, -1.759662411672122E-03_wp,  8.482048844191342E-02_wp, &
   &  3.388121259468538E-03_wp,  1.836313813427145E-02_wp,  2.974977636596358E-01_wp, &
   &  2.142340779274882E-02_wp, -2.440427953757240E-02_wp,  7.502270778977416E-03_wp, &
   & -3.559797879290408E-03_wp,  6.772742090462043E-02_wp,  1.487991423992202E-01_wp, &
   & -2.451460387601280E-01_wp, -5.172003090114057E-03_wp,  5.290270182577418E-03_wp, &
   & -2.483563415016234E-03_wp,  9.466882832942922E-03_wp,  6.822251226400594E-04_wp, &
   &  1.703210708227039E-01_wp, -4.213592394024056E-02_wp,  6.307130221910287E-02_wp, &
   &  1.939267990732981E-01_wp,  5.294122943726928E-03_wp, -9.600813313738004E-03_wp, &
   & -4.616744038439015E-03_wp,  8.115573032660257E-03_wp, -6.279131220983131E-03_wp, &
   & -1.168795641010645E-02_wp, -6.661092935395363E-03_wp,  7.093392456161390E-02_wp, &
   & -1.039840771817277E-03_wp,  9.499478900720760E-04_wp,  6.940628533841650E-04_wp, &
   &  4.893621576621127E-04_wp,  8.440838334803196E-01_wp,  4.905098956976956E-03_wp, &
   &  1.357208373914456E-02_wp,  1.393150004901561E-01_wp, -7.672996513895941E-04_wp, &
   &  2.247057077481621E-03_wp,  4.209175505347790E-03_wp, -1.058023478311760E-03_wp, &
   &  1.893107319151598E-01_wp,  8.285327706836283E-03_wp,  6.217682393742921E-03_wp, &
   & -2.752109323970900E-01_wp, -2.319576000115207E-02_wp, -1.543709149218754E-02_wp, &
   &  1.217820086890074E-02_wp,  9.942851359656751E-03_wp, -1.653780775395164E-01_wp, &
   & -1.709439446336162E-01_wp, -2.641444438990908E-01_wp,  1.784668723847591E-01_wp, &
   & -2.592953444016959E-02_wp,  1.302343834342591E-02_wp, -1.779073162236685E-02_wp, &
   & -2.491510203795122E-03_wp, -8.145251752221179E-02_wp,  1.563599543080730E-02_wp, &
   &  3.077920992193937E-02_wp,  3.132114535967224E-01_wp, -1.460952102501112E-02_wp, &
   &  2.147256620212746E-03_wp, -1.966606484542853E-03_wp, -2.359994158908571E-04_wp, &
   & -1.323594660617543E-01_wp,  1.487433493467488E-01_wp,  2.512430203608930E-01_wp, &
   &  1.533997236662861E-01_wp, -1.266821157991563E-02_wp, -5.673665778168770E-03_wp, &
   &  3.101807402292274E-03_wp, -2.473754141221987E-03_wp, -1.480200394334560E-01_wp, &
   &  2.682219595221672E-02_wp,  6.117560019213610E-02_wp,  2.169113743768337E-01_wp, &
   & -4.601741626533677E-04_wp,  3.239728633724738E-02_wp, -1.044545316910896E-03_wp, &
   & -1.209850983039165E-02_wp,  9.441813035990097E-04_wp,  5.167756375487324E-03_wp, &
   &  3.581781097085035E-03_wp,  3.207817638247836E-03_wp,  4.905098956976956E-03_wp, &
   &  4.897814171966381E-01_wp,  5.734599801011469E-03_wp, -5.054332082185083E-03_wp, &
   & -2.307349711758816E-03_wp,  8.249685263131983E-03_wp,  5.820539241638587E-03_wp, &
   & -8.881627805059455E-04_wp, -8.528644720827837E-03_wp,  2.963016928987688E-02_wp, &
   & -3.126916799135237E-02_wp,  1.208270934726472E-02_wp,  3.258950771850526E-03_wp, &
   &  3.195851458648564E-03_wp, -1.526228509487418E-03_wp,  2.870855741916439E-03_wp, &
   & -1.781838984708946E-02_wp,  2.891932065592874E-01_wp, -4.457175475793600E-01_wp, &
   &  1.958397655878534E-01_wp, -6.247785722943803E-03_wp,  3.551522894860388E-03_wp, &
   & -5.251071076191872E-03_wp,  9.054567357963395E-03_wp,  7.266027835159628E-03_wp, &
   & -1.466010502540641E-02_wp,  1.037978655013428E-02_wp, -3.279230862883733E-02_wp, &
   & -1.254445802445833E-02_wp,  1.418234526306885E-02_wp,  1.525759156952472E-03_wp, &
   & -8.086516903264363E-03_wp, -2.121388674350320E-03_wp,  1.500688127451139E-01_wp, &
   & -3.926375660752298E-01_wp, -1.032282700263639E-01_wp,  1.001834442550433E-02_wp, &
   &  1.979827992580431E-02_wp,  1.853235953651607E-03_wp,  2.987968882560821E-03_wp, &
   &  1.700385822401850E-03_wp, -1.661051130919642E-02_wp,  8.539519630381751E-03_wp, &
   &  3.987033032478870E-02_wp,  7.192266070701130E-04_wp,  5.872652032400529E-02_wp, &
   &  1.328230984236084E-03_wp, -1.515898622269251E-02_wp, -7.025956212604888E-04_wp, &
   & -3.588284470306161E-03_wp, -1.759937428925466E-03_wp, -6.516583487466641E-03_wp, &
   &  1.357208373914456E-02_wp,  5.734599801011469E-03_wp,  4.957485235927724E-01_wp, &
   & -4.020065693571041E-03_wp,  4.253023278789240E-03_wp, -5.805874148151157E-03_wp, &
   & -2.402379390926253E-03_wp,  1.140437619749657E-03_wp, -6.409071717277027E-03_wp, &
   & -3.129140864602862E-02_wp, -2.352591384986751E-03_wp,  4.550801438672616E-03_wp, &
   & -8.447627089361608E-03_wp, -1.024918492357502E-02_wp, -9.832928613073454E-03_wp, &
   & -4.755797903129810E-03_wp, -2.704219484225604E-02_wp, -4.452416252597802E-01_wp, &
   & -1.070319465842709E-01_wp,  2.847711408940242E-01_wp,  6.666155912529718E-03_wp, &
   & -1.633220298493662E-02_wp, -1.158231976716216E-02_wp,  3.173781797194465E-04_wp, &
   & -7.536133854468864E-03_wp,  1.741201503734447E-02_wp,  1.319859421964997E-02_wp, &
   & -4.053426550408897E-02_wp,  6.443821487081486E-03_wp, -1.282593655401865E-02_wp, &
   & -4.838100087988706E-04_wp,  3.409247581638760E-03_wp, -1.063982081247404E-02_wp, &
   & -3.905781596073190E-01_wp, -2.338639896244302E-01_wp, -1.745248313432676E-01_wp, &
   & -7.608286765451494E-03_wp, -2.533577437102717E-02_wp, -7.433462681310299E-03_wp, &
   &  4.005414443716922E-03_wp,  8.116145381064886E-03_wp,  2.544418327570327E-02_wp, &
   &  5.209091343999873E-03_wp,  6.041076413473551E-02_wp,  6.857108714237196E-04_wp, &
   &  9.837219950312542E-02_wp, -2.359384643701530E-03_wp,  1.257751136631665E-01_wp, &
   & -4.875661919418912E-04_wp, -3.182467059190740E-03_wp, -6.500358384336613E-03_wp, &
   &  1.251090620324090E-02_wp,  1.393150004901561E-01_wp, -5.054332082185083E-03_wp, &
   & -4.020065693571041E-03_wp,  6.349376627643307E-01_wp, -1.055599311096726E-03_wp, &
   &  8.319828994270612E-04_wp,  1.072738527777089E-03_wp, -1.027305973761432E-02_wp, &
   &  2.754122263302456E-01_wp,  1.173123880273071E-02_wp,  4.375154962221900E-03_wp, &
   & -4.029130250229577E-01_wp, -2.320237897976008E-02_wp, -1.048011620248967E-02_wp, &
   &  5.013979106768630E-03_wp, -1.771176088237056E-03_wp, -4.236455577071013E-02_wp, &
   &  2.121421887890352E-01_wp,  3.057313329878171E-01_wp, -1.457122936094603E-02_wp, &
   & -2.139127442612924E-02_wp,  2.444029883229076E-02_wp,  7.502693923921305E-03_wp, &
   & -3.596138573388669E-03_wp, -8.475244399355290E-02_wp, -3.181080000113450E-03_wp, &
   &  1.863860377306688E-02_wp,  2.979406044790362E-01_wp, -5.286287347104757E-03_wp, &
   &  2.474080389321256E-03_wp,  9.464848782751618E-03_wp,  6.742819824984703E-04_wp, &
   & -6.784656854515582E-02_wp, -1.485893963175227E-01_wp, -2.452245826391641E-01_wp, &
   & -5.214816827278326E-03_wp, -5.256371142047481E-03_wp,  9.637880193722827E-03_wp, &
   & -4.589717515932207E-03_wp,  8.060262539910920E-03_wp, -1.691793086208223E-01_wp, &
   &  4.207039178782482E-02_wp,  6.303054044855709E-02_wp,  1.943019231543751E-01_wp, &
   &  7.434613310779223E-02_wp, -6.703513719039949E-03_wp, -1.196966000154187E-02_wp, &
   & -6.250505483431420E-03_wp,  1.891309702292605E-01_wp, -8.349427620557459E-03_wp, &
   &  6.156068262093434E-03_wp, -2.753678417295266E-01_wp, -7.672996513895941E-04_wp, &
   & -2.307349711758816E-03_wp,  4.253023278789240E-03_wp, -1.055599311096722E-03_wp, &
   &  8.429243436740793E-01_wp, -5.413530946956013E-03_wp,  1.314471515825148E-02_wp, &
   &  1.380122855737993E-01_wp, -1.061303902389656E-03_wp, -9.265554151292729E-04_wp, &
   &  6.888225060027338E-04_wp,  4.656231092177304E-04_wp, -8.139473464426875E-02_wp, &
   & -1.533919345367875E-02_wp,  3.032257098052348E-02_wp,  3.127931244323070E-01_wp, &
   & -2.599185858576010E-02_wp, -1.308363624836873E-02_wp, -1.773278463004659E-02_wp, &
   & -2.470383346368663E-03_wp, -1.657269212141140E-01_wp,  1.708975073338503E-01_wp, &
   & -2.642221834488391E-01_wp,  1.789397839834518E-01_wp, -2.316632131853606E-02_wp, &
   &  1.542197400806148E-02_wp,  1.221152261208723E-02_wp,  9.981223456534682E-03_wp, &
   & -1.491410307302132E-01_wp, -2.671754073244886E-02_wp,  6.101614619024202E-02_wp, &
   &  2.170839367149032E-01_wp, -1.267048579151314E-02_wp,  5.701640309369226E-03_wp, &
   &  3.116726806188214E-03_wp, -2.476142530256957E-03_wp, -1.321078178110328E-01_wp, &
   & -1.488842647684905E-01_wp,  2.510445353237979E-01_wp,  1.534613979273354E-01_wp, &
   & -1.460755041629342E-02_wp, -2.136201982272271E-03_wp, -1.997947093469191E-03_wp, &
   & -2.491047885149165E-04_wp,  1.286760939909034E-02_wp,  1.021237784893391E-03_wp, &
   & -3.242534561842155E-02_wp,  5.026181659209923E-04_wp,  8.367566084559982E-03_wp, &
   &  3.010369337275243E-02_wp,  3.143822674633701E-02_wp, -1.190960129238832E-02_wp, &
   &  2.247057077481621E-03_wp,  8.249685263131983E-03_wp, -5.805874148151157E-03_wp, &
   &  8.319828994270617E-04_wp, -5.413530946956012E-03_wp,  4.900871807697212E-01_wp, &
   & -5.596756913160784E-03_wp,  4.995209703286868E-03_wp, -9.858064083849093E-04_wp, &
   &  5.198954990384540E-03_wp, -3.532797347700027E-03_wp, -3.108967316539873E-03_wp, &
   & -7.222939499670370E-03_wp, -1.463712866321349E-02_wp, -1.037421588469870E-02_wp, &
   &  3.265421205689111E-02_wp,  6.249373377971345E-03_wp,  3.491074201379736E-03_wp, &
   &  5.271139416584750E-03_wp, -8.917424258505349E-03_wp,  1.803693884379085E-02_wp, &
   &  2.893247956503804E-01_wp,  4.456639329301224E-01_wp, -1.958230124990306E-01_wp, &
   & -3.247202360336833E-03_wp,  3.297768879419099E-03_wp,  1.593243135905389E-03_wp, &
   & -2.842032462672788E-03_wp, -1.774536610663534E-03_wp, -1.683422150244912E-02_wp, &
   & -8.703708739375222E-03_wp, -4.024801588740388E-02_wp, -1.000069826473252E-02_wp, &
   &  1.987396831224482E-02_wp, -1.826680563705499E-03_wp, -2.969229187384274E-03_wp, &
   &  2.410365605460931E-03_wp,  1.499476098159285E-01_wp,  3.926132423729201E-01_wp, &
   &  1.032751244669945E-01_wp,  1.246165939538359E-02_wp,  1.420831397488669E-02_wp, &
   & -1.507519814669938E-03_wp,  7.879338749816585E-03_wp, -1.631890529088685E-02_wp, &
   &  1.303333505416583E-03_wp,  5.872483092002760E-02_wp,  7.919192080460943E-04_wp, &
   & -6.144229616403341E-03_wp,  3.143850179021387E-02_wp, -2.005312615417137E-03_wp, &
   &  4.394774512976690E-03_wp,  4.209175505347790E-03_wp,  5.820539241638587E-03_wp, &
   & -2.402379390926253E-03_wp,  1.072738527777089E-03_wp,  1.314471515825148E-02_wp, &
   & -5.596756913160784E-03_wp,  4.960974267090863E-01_wp, -4.141781656964817E-03_wp, &
   & -6.743778300490629E-04_wp,  3.638971667978895E-03_wp, -1.720320164099396E-03_wp, &
   & -6.438766157513710E-03_wp, -7.515659630910805E-03_wp, -1.737030673046555E-02_wp, &
   &  1.321063980480507E-02_wp, -4.063934982645971E-02_wp,  6.668506526406026E-03_wp, &
   &  1.629544532065390E-02_wp, -1.153080802371634E-02_wp,  4.005588825677863E-04_wp, &
   & -2.701109625402425E-02_wp,  4.451157004621552E-01_wp, -1.069518229891273E-01_wp, &
   &  2.849175611106504E-01_wp, -8.452122821130040E-03_wp,  1.027468115563575E-02_wp, &
   & -9.831466240989489E-03_wp, -4.756019537487924E-03_wp,  8.265497770558025E-03_wp, &
   & -2.559730740604082E-02_wp,  5.061542100640776E-03_wp,  6.058992671227835E-02_wp, &
   & -7.583582685916045E-03_wp,  2.536263019625774E-02_wp, -7.410634617386340E-03_wp, &
   &  4.014880863428783E-03_wp, -1.066620604694135E-02_wp,  3.907223175034765E-01_wp, &
   & -2.336135451174082E-01_wp, -1.745923449463822E-01_wp,  6.352354891197232E-03_wp, &
   &  1.283705618323999E-02_wp, -4.543834522675481E-04_wp,  3.266875793296565E-03_wp, &
   &  1.298721004409743E-01_wp, -2.390897664126386E-03_wp,  9.798681880857191E-02_wp, &
   &  6.612697887665494E-04_wp,  2.753621765898573E-01_wp, -1.188768059898162E-02_wp, &
   &  4.403932312142212E-03_wp, -4.023182634145880E-01_wp, -1.058023478311763E-03_wp, &
   & -8.881627805059455E-04_wp,  1.140437619749657E-03_wp, -1.027305973761432E-02_wp, &
   &  1.380122855737993E-01_wp,  4.995209703286868E-03_wp, -4.141781656964817E-03_wp, &
   &  6.337001704949426E-01_wp, -4.720398103786848E-04_wp,  3.147109815353715E-03_wp, &
   & -6.459458576537677E-03_wp,  1.250007039652695E-02_wp, -8.481778745057701E-02_wp, &
   &  3.384726228089784E-03_wp,  1.832842557297677E-02_wp,  2.974753120073315E-01_wp, &
   & -2.144652430875094E-02_wp, -2.441402512761884E-02_wp,  7.503367099585823E-03_wp, &
   & -3.553397090138383E-03_wp, -4.241743832695642E-02_wp, -2.120264316741839E-01_wp, &
   &  3.057743526797048E-01_wp, -1.436800564547772E-02_wp, -2.318548203365187E-02_wp, &
   &  1.053919970004938E-02_wp,  4.979510908235431E-03_wp, -1.753212643164456E-03_wp, &
   & -1.703238437187407E-01_wp, -4.212739390044036E-02_wp,  6.303042784252559E-02_wp, &
   &  1.939408725198874E-01_wp, -5.268995145825184E-03_wp, -9.590785761891103E-03_wp, &
   & -4.626907472414760E-03_wp,  8.044404888566635E-03_wp, -6.771802335365747E-02_wp, &
   &  1.487610897931907E-01_wp, -2.451443254802192E-01_wp, -5.173852886849756E-03_wp, &
   & -5.332602643252117E-03_wp, -2.446411327017794E-03_wp,  9.479628340630961E-03_wp, &
   &  6.297357615509258E-04_wp, -6.679887117572363E-03_wp,  7.444902492001537E-02_wp, &
   & -6.251871088207495E-03_wp, -1.272759551324487E-02_wp, -7.895148297830563E-04_wp, &
   & -2.352905761865181E-03_wp, -4.302988406396988E-03_wp,  1.055441718594612E-03_wp, &
   &  1.893107319151598E-01_wp, -8.528644720827839E-03_wp, -6.409071717277027E-03_wp, &
   &  2.754122263302456E-01_wp, -1.061303902389652E-03_wp, -9.858064083849093E-04_wp, &
   & -6.743778300490630E-04_wp, -4.720398103786813E-04_wp,  8.407661590793284E-01_wp, &
   & -5.677688409346515E-03_wp, -1.355150160609650E-02_wp, -1.364764058780396E-01_wp, &
   & -2.595350954541548E-02_wp, -1.324295687161667E-02_wp,  1.775204009372663E-02_wp, &
   &  2.564204627362940E-03_wp, -8.086244878243089E-02_wp, -1.527232052331534E-02_wp, &
   & -3.035138164936899E-02_wp, -3.135596086621796E-01_wp, -2.318708098773159E-02_wp, &
   &  1.548467799425112E-02_wp, -1.218955416095540E-02_wp, -9.994284547877375E-03_wp, &
   & -1.658753423440057E-01_wp,  1.713553084114756E-01_wp,  2.647943143328361E-01_wp, &
   & -1.785910201043079E-01_wp, -1.260646194170170E-02_wp,  5.592320105044995E-03_wp, &
   & -3.104053425941847E-03_wp,  2.438377742502516E-03_wp, -1.491780784715833E-01_wp, &
   & -2.663163457527537E-02_wp, -6.121419307829655E-02_wp, -2.169361995806018E-01_wp, &
   & -1.457201130417061E-02_wp, -2.027423210222607E-03_wp,  1.980517257124989E-03_wp, &
   &  2.242175637279242E-04_wp, -1.310292759259390E-01_wp, -1.497798779426327E-01_wp, &
   & -2.524053952636748E-01_wp, -1.496382287741379E-01_wp,  1.072389224409596E-03_wp, &
   &  1.281459086295719E-02_wp,  4.174789300134937E-04_wp, -3.216708991641609E-02_wp, &
   &  2.255651556425698E-03_wp,  8.293513906870725E-03_wp,  5.876240630642964E-03_wp, &
   & -9.481975499901815E-04_wp,  8.285327706836283E-03_wp,  2.963016928987688E-02_wp, &
   & -3.129140864602862E-02_wp,  1.173123880273071E-02_wp, -9.265554151292727E-04_wp, &
   &  5.198954990384540E-03_wp,  3.638971667978895E-03_wp,  3.147109815353715E-03_wp, &
   & -5.677688409346515E-03_wp,  4.904617372594265E-01_wp,  4.905872042167136E-03_wp, &
   & -4.693617839185362E-03_wp,  6.413622381211287E-03_wp,  3.683762713926030E-03_wp, &
   & -5.325594795411811E-03_wp,  8.982476159279553E-03_wp, -7.284095090466867E-03_wp, &
   & -1.508384427528231E-02_wp,  1.055241187815559E-02_wp, -3.263312566941264E-02_wp, &
   & -3.408673525016986E-03_wp,  3.368833142021846E-03_wp, -1.663002589937036E-03_wp, &
   &  2.746552100602511E-03_wp,  1.795448754187384E-02_wp,  2.894669651689678E-01_wp, &
   & -4.456721790724243E-01_wp,  1.955411855967179E-01_wp, -1.000012702991848E-02_wp, &
   &  1.970151246932410E-02_wp,  1.760256337875836E-03_wp,  3.012932969958350E-03_wp, &
   & -1.687030687679118E-03_wp, -1.700436400697725E-02_wp,  8.892644060682691E-03_wp, &
   &  4.026448777005662E-02_wp,  1.253150227241959E-02_wp,  1.406281643489339E-02_wp, &
   &  1.440319008298704E-03_wp, -8.079919789563846E-03_wp,  2.288466861541040E-03_wp, &
   &  1.501352108339646E-01_wp, -3.931103870497424E-01_wp, -1.029729431720197E-01_wp, &
   & -1.355789093345902E-03_wp,  1.625674221823976E-02_wp, -7.021794444146690E-04_wp, &
   & -5.856052840859258E-02_wp, -4.236227136247043E-03_wp, -5.786828146520061E-03_wp, &
   & -2.416968585163443E-03_wp,  1.155982015342645E-03_wp,  6.217682393742921E-03_wp, &
   & -3.126916799135237E-02_wp, -2.352591384986751E-03_wp,  4.375154962221900E-03_wp, &
   &  6.888225060027338E-04_wp, -3.532797347700027E-03_wp, -1.720320164099396E-03_wp, &
   & -6.459458576537677E-03_wp, -1.355150160609650E-02_wp,  4.905872042167136E-03_wp, &
   &  4.958265611318840E-01_wp, -3.628005973533152E-03_wp, -6.718447216272541E-03_wp, &
   & -1.645504776088864E-02_wp, -1.149157679255069E-02_wp,  3.584188089654364E-04_wp, &
   &  7.415848264184872E-03_wp,  1.754113370741720E-02_wp,  1.289431821764276E-02_wp, &
   & -4.081635698877482E-02_wp,  8.511548438474442E-03_wp, -1.031655271816231E-02_wp, &
   & -9.727541059417626E-03_wp, -4.674478667837700E-03_wp,  2.684783123649586E-02_wp, &
   & -4.451203244006532E-01_wp, -1.068394515688850E-01_wp,  2.846339733760712E-01_wp, &
   &  7.628098106089499E-03_wp, -2.527738044445453E-02_wp, -7.398974178769744E-03_wp, &
   &  3.985294758549605E-03_wp, -8.159842431310724E-03_wp,  2.564908723643024E-02_wp, &
   &  4.937990448470266E-03_wp,  6.064273938241257E-02_wp, -6.439218353030748E-03_wp, &
   & -1.275554334950709E-02_wp, -4.352368225150864E-04_wp,  3.404030448830548E-03_wp, &
   &  1.048516538151298E-02_wp, -3.912198384430868E-01_wp, -2.339238248737927E-01_wp, &
   & -1.738028570516824E-01_wp,  2.361020252331918E-03_wp, -1.302057363170864E-01_wp, &
   & -6.985629288825479E-04_wp, -9.623347293437912E-02_wp,  1.053909754895402E-03_wp, &
   &  8.788780580587822E-04_wp,  1.120519030515370E-03_wp, -1.027700775080411E-02_wp, &
   & -2.752109323970900E-01_wp,  1.208270934726472E-02_wp,  4.550801438672616E-03_wp, &
   & -4.029130250229577E-01_wp,  4.656231092177304E-04_wp, -3.108967316539873E-03_wp, &
   & -6.438766157513710E-03_wp,  1.250007039652695E-02_wp, -1.364764058780396E-01_wp, &
   & -4.693617839185362E-03_wp, -3.628005973533152E-03_wp,  6.342868118918844E-01_wp, &
   &  2.138547919130964E-02_wp,  2.444185879871064E-02_wp,  7.485209089154482E-03_wp, &
   & -3.597661001285501E-03_wp,  8.399271700388317E-02_wp, -3.474362913185432E-03_wp, &
   &  1.852404496066776E-02_wp,  2.989390691340027E-01_wp,  2.312396198293607E-02_wp, &
   & -1.042978442278318E-02_wp,  5.007225909869392E-03_wp, -1.762373847291776E-03_wp, &
   &  4.204746720219069E-02_wp,  2.121026114011793E-01_wp,  3.061856776281458E-01_wp, &
   & -1.332693544380374E-02_wp,  5.250219250263501E-03_wp,  9.621043622666827E-03_wp, &
   & -4.610957486082562E-03_wp,  8.077821380194819E-03_wp,  1.706126424970112E-01_wp, &
   &  4.204645221277876E-02_wp,  6.342946382614620E-02_wp,  1.939857159631502E-01_wp, &
   &  5.312023266393023E-03_wp,  2.457308742080725E-03_wp,  9.423759394375709E-03_wp, &
   &  6.756893186891516E-04_wp,  6.511143463338417E-02_wp, -1.483828137330170E-01_wp, &
   & -2.444355309205807E-01_wp, -6.722946300012508E-03_wp, -3.832625659357759E-02_wp, &
   & -2.910217311767796E-04_wp,  4.911014170955832E-02_wp,  2.078487962056579E-04_wp, &
   & -1.657257723234933E-01_wp, -1.803512460547285E-02_wp,  2.702426709506782E-02_wp, &
   &  4.241317133803125E-02_wp, -2.319576000115207E-02_wp,  3.258950771850526E-03_wp, &
   & -8.447627089361608E-03_wp, -2.320237897976008E-02_wp, -8.139473464426875E-02_wp, &
   & -7.222939499670370E-03_wp, -7.515659630910805E-03_wp, -8.481778745057701E-02_wp, &
   & -2.595350954541548E-02_wp,  6.413622381211287E-03_wp, -6.718447216272541E-03_wp, &
   &  2.138547919130964E-02_wp,  1.926848483643393E+00_wp,  2.523216717240165E-01_wp, &
   & -3.744497202212876E-01_wp,  2.079336600638556E-01_wp, -5.777826636612601E-03_wp, &
   &  2.153730360168114E-02_wp, -1.457198913841331E-02_wp,  1.263999023166349E-03_wp, &
   &  1.056851822576655E-02_wp, -2.409212357888751E-03_wp, -3.129192927315042E-02_wp, &
   &  6.532425843818973E-02_wp, -2.350629245771730E-02_wp,  5.056989348941861E-02_wp, &
   & -3.076405058003849E-02_wp, -4.543925211347843E-03_wp, -5.717935275029361E-03_wp, &
   & -5.390317807562760E-02_wp,  7.309749637482553E-02_wp, -9.714300020215932E-02_wp, &
   & -1.573785119430015E-02_wp,  3.262800478449652E-02_wp, -6.973219408918997E-03_wp, &
   &  1.143549362508695E-02_wp, -1.300895199971275E-01_wp,  1.099100123152549E-02_wp, &
   & -2.864002799454562E-03_wp,  2.243791573011786E-01_wp, -2.023228360698578E-02_wp, &
   &  1.666726096580421E-02_wp, -1.215357414037543E-02_wp,  5.164434699717015E-03_wp, &
   &  4.198805679545000E-02_wp, -4.251158967129975E-04_wp, -8.470644344694444E-03_wp, &
   & -1.030952412861568E-02_wp, -1.708980155892887E-01_wp,  2.893336223420668E-01_wp, &
   &  4.451127413192278E-01_wp, -2.119840694276171E-01_wp, -1.543709149218754E-02_wp, &
   &  3.195851458648564E-03_wp, -1.024918492357502E-02_wp, -1.048011620248967E-02_wp, &
   & -1.533919345367875E-02_wp, -1.463712866321349E-02_wp, -1.737030673046555E-02_wp, &
   &  3.384726228089781E-03_wp, -1.324295687161667E-02_wp,  3.683762713926030E-03_wp, &
   & -1.645504776088864E-02_wp,  2.444185879871064E-02_wp,  2.523216717240165E-01_wp, &
   &  1.369474203578749E+00_wp,  6.794320592905591E-02_wp, -1.906212965051925E-01_wp, &
   &  2.157697344921147E-02_wp,  7.978963797049222E-02_wp, -1.371073448303978E-03_wp, &
   &  2.087773084448162E-02_wp,  2.409764949665100E-03_wp, -6.684545227570542E-02_wp, &
   & -3.909244019541337E-02_wp, -1.790214766124290E-02_wp, -5.060168129579690E-02_wp, &
   &  5.269332533863068E-02_wp,  1.027592095432595E-02_wp,  1.450273323327533E-03_wp, &
   &  6.005522562545247E-02_wp, -2.513313120542579E-01_wp, -3.861543712368549E-01_wp, &
   &  4.879030348393663E-02_wp, -1.852974369621770E-02_wp,  1.232878797800032E-02_wp, &
   &  6.188808398671801E-03_wp,  1.223242374577508E-02_wp,  1.429155078064441E-02_wp, &
   & -1.083191611256763E-01_wp, -4.031613472964876E-02_wp, -2.774561432925469E-02_wp, &
   &  1.715256929930816E-02_wp,  8.203746951009575E-02_wp, -4.044006799140985E-03_wp, &
   & -9.652428849881690E-03_wp, -6.413881380697377E-02_wp, -8.169084722126841E-03_wp, &
   &  7.545304111834568E-03_wp, -4.693807740685106E-03_wp,  2.642395600556753E-01_wp, &
   &  4.456579031710494E-01_wp, -1.070069421144721E-01_wp,  3.057780897872059E-01_wp, &
   &  1.217820086890074E-02_wp, -1.526228509487418E-03_wp, -9.832928613073454E-03_wp, &
   &  5.013979106768630E-03_wp,  3.032257098052349E-02_wp, -1.037421588469870E-02_wp, &
   &  1.321063980480507E-02_wp,  1.832842557297677E-02_wp,  1.775204009372663E-02_wp, &
   & -5.325594795411811E-03_wp, -1.149157679255069E-02_wp,  7.485209089154482E-03_wp, &
   & -3.744497202212876E-01_wp,  6.794320592905591E-02_wp,  1.320567713468630E+00_wp, &
   &  3.129704518993774E-01_wp,  1.462018052720512E-02_wp,  1.358868318415217E-03_wp, &
   &  1.978877150861749E-02_wp,  1.550653626890516E-02_wp,  3.130074965933824E-02_wp, &
   & -3.908937699891116E-02_wp,  9.695864526332573E-03_wp,  3.720983979051862E-02_wp, &
   & -3.079888151813678E-02_wp, -1.012824449158330E-02_wp,  5.213976022273503E-02_wp, &
   & -2.086277882728174E-02_wp, -9.901469545667278E-02_wp, -3.957205218392055E-01_wp, &
   &  1.265239701221103E-01_wp, -5.861095255407232E-02_wp, -1.029243857275900E-02_wp, &
   & -1.401008290642191E-02_wp,  3.775030908808071E-02_wp,  1.963260925995177E-02_wp, &
   &  6.933316499867500E-03_wp, -4.201691464225484E-02_wp, -7.274251842956297E-02_wp, &
   &  3.770254238684999E-02_wp,  1.865495047479866E-03_wp, -8.635557223496085E-03_wp, &
   &  1.338003834431531E-02_wp, -6.908828442012717E-03_wp,  5.753380864868125E-02_wp, &
   &  8.450451233251027E-04_wp,  1.550495190069728E-01_wp,  1.293707498310058E-03_wp, &
   & -1.789389989696000E-01_wp, -1.957747719980853E-01_wp,  2.849107342135790E-01_wp, &
   & -1.438377687257794E-02_wp,  9.942851359656737E-03_wp,  2.870855741916439E-03_wp, &
   & -4.755797903129810E-03_wp, -1.771176088237056E-03_wp,  3.127931244323070E-01_wp, &
   &  3.265421205689111E-02_wp, -4.063934982645971E-02_wp,  2.974753120073315E-01_wp, &
   &  2.564204627362926E-03_wp,  8.982476159279553E-03_wp,  3.584188089654364E-04_wp, &
   & -3.597661001285488E-03_wp,  2.079336600638556E-01_wp, -1.906212965051925E-01_wp, &
   &  3.129704518993774E-01_wp,  1.626991060067219E+00_wp, -1.232789978303808E-03_wp, &
   & -2.093591850124991E-02_wp,  1.554724019872558E-02_wp,  1.210343840871203E-02_wp, &
   & -6.533595999018363E-02_wp, -1.791338449551878E-02_wp,  3.718892145420852E-02_wp, &
   & -1.041156167570835E-01_wp, -4.518157760969749E-03_wp, -1.532375845760117E-03_wp, &
   & -2.085969784006983E-02_wp, -4.433865201468268E-03_wp, -1.855302611096530E-01_wp, &
   & -8.374764020602153E-02_wp,  1.376598521306487E-01_wp,  1.101342225395888E-03_wp, &
   & -1.555877947121680E-02_wp, -8.035277723205586E-03_wp,  2.886706449458729E-02_wp, &
   &  7.146330771728562E-03_wp, -1.820906817729780E-01_wp,  3.857194786084731E-02_wp, &
   & -5.465485333034970E-02_wp,  3.483625814784412E-01_wp, -2.462283074712185E-02_wp, &
   & -9.868078060312466E-03_wp, -4.459001831734477E-04_wp, -3.127256758552299E-03_wp, &
   & -3.040323340029377E-04_wp, -3.851718429301771E-02_wp,  1.132634881668910E-04_wp, &
   &  5.203536420446843E-02_wp, -2.318060720164255E-02_wp,  3.297270138654383E-03_wp, &
   &  8.499489701407577E-03_wp,  2.319570021325564E-02_wp, -1.653780775395164E-01_wp, &
   & -1.781838984708946E-02_wp, -2.704219484225604E-02_wp, -4.236455577071013E-02_wp, &
   & -2.599185858576010E-02_wp,  6.249373377971345E-03_wp,  6.668506526406026E-03_wp, &
   & -2.144652430875094E-02_wp, -8.086244878243089E-02_wp, -7.284095090466867E-03_wp, &
   &  7.415848264184872E-03_wp,  8.399271700388317E-02_wp, -5.777826636612601E-03_wp, &
   &  2.157697344921147E-02_wp,  1.462018052720513E-02_wp, -1.232789978303808E-03_wp, &
   &  1.925190137555957E+00_wp,  2.527509683417552E-01_wp,  3.750063822716700E-01_wp, &
   & -2.067220167151342E-01_wp, -2.354356998774298E-02_wp,  5.058114623212935E-02_wp, &
   &  3.076438051088457E-02_wp,  4.456717413048608E-03_wp,  1.041404345108276E-02_wp, &
   & -2.521274400848204E-03_wp,  3.119459860731957E-02_wp, -6.512169792188434E-02_wp, &
   & -1.578594610730463E-02_wp,  3.269703607432843E-02_wp,  6.964522965623769E-03_wp, &
   & -1.142018525754264E-02_wp, -5.688597060383519E-03_wp, -5.384009938581511E-02_wp, &
   & -7.288773528300484E-02_wp,  9.743519498628025E-02_wp, -2.021590405397043E-02_wp, &
   &  1.672461362425607E-02_wp,  1.223352579956685E-02_wp, -5.276732047483599E-03_wp, &
   & -1.318657238440606E-01_wp,  1.087088887548559E-02_wp,  2.919030379864971E-03_wp, &
   & -2.260728331184097E-01_wp, -4.493523318597993E-04_wp,  4.207907422154082E-02_wp, &
   & -1.029895438504152E-02_wp, -9.234309186277682E-03_wp, -1.544549398688108E-02_wp, &
   &  3.233388733204503E-03_wp,  1.023902985142095E-02_wp,  1.049207368818320E-02_wp, &
   & -1.709439446336162E-01_wp,  2.891932065592874E-01_wp, -4.452416252597802E-01_wp, &
   &  2.121421887890352E-01_wp, -1.308363624836873E-02_wp,  3.491074201379736E-03_wp, &
   &  1.629544532065390E-02_wp, -2.441402512761885E-02_wp, -1.527232052331534E-02_wp, &
   & -1.508384427528231E-02_wp,  1.754113370741720E-02_wp, -3.474362913185432E-03_wp, &
   &  2.153730360168114E-02_wp,  7.978963797049222E-02_wp,  1.358868318415217E-03_wp, &
   & -2.093591850124991E-02_wp,  2.527509683417552E-01_wp,  1.369469894449481E+00_wp, &
   & -6.826659247602748E-02_wp,  1.903041744449289E-01_wp, -5.057698564821952E-02_wp, &
   &  5.275107182007691E-02_wp, -1.012976688131399E-02_wp, -1.534249792763384E-03_wp, &
   &  2.367165320956746E-03_wp, -6.641376893124294E-02_wp,  3.887425976509527E-02_wp, &
   &  1.810891325053174E-02_wp, -1.849640619419218E-02_wp,  1.237647302717487E-02_wp, &
   & -6.178139001288267E-03_wp, -1.225637791523566E-02_wp,  5.992976277236097E-02_wp, &
   & -2.511369516436848E-01_wp,  3.859972450964884E-01_wp, -4.872401455665089E-02_wp, &
   &  1.710229838831284E-02_wp,  8.217955684208357E-02_wp,  4.089208890809513E-03_wp, &
   &  9.652999950917271E-03_wp,  1.484411960351194E-02_wp, -1.082103747147362E-01_wp, &
   &  4.020646820081600E-02_wp,  2.841422012598827E-02_wp,  8.203910821059849E-03_wp, &
   &  6.421113390698539E-02_wp,  4.659438803525548E-03_wp, -9.221765000020827E-03_wp, &
   & -1.220172254632664E-02_wp,  1.591282222122241E-03_wp, -9.798591834190845E-03_wp, &
   &  5.002988700501795E-03_wp, -2.641444438990908E-01_wp, -4.457175475793600E-01_wp, &
   & -1.070319465842709E-01_wp,  3.057313329878171E-01_wp, -1.773278463004659E-02_wp, &
   &  5.271139416584750E-03_wp, -1.153080802371634E-02_wp,  7.503367099585816E-03_wp, &
   & -3.035138164936899E-02_wp,  1.055241187815559E-02_wp,  1.289431821764276E-02_wp, &
   &  1.852404496066776E-02_wp, -1.457198913841331E-02_wp, -1.371073448303985E-03_wp, &
   &  1.978877150861749E-02_wp,  1.554724019872558E-02_wp,  3.750063822716700E-01_wp, &
   & -6.826659247602748E-02_wp,  1.320864734011985E+00_wp,  3.122766792850088E-01_wp, &
   &  3.077026689983809E-02_wp,  1.012784983664287E-02_wp,  5.210857233250323E-02_wp, &
   & -2.083376288009661E-02_wp, -3.134950464921349E-02_wp,  3.884260085549837E-02_wp, &
   &  9.898472207251173E-03_wp,  3.721641606431228E-02_wp,  1.029362979057584E-02_wp, &
   &  1.395801224594882E-02_wp,  3.774051610790029E-02_wp,  1.960703470239653E-02_wp, &
   &  9.891979102940525E-02_wp,  3.953897623034650E-01_wp,  1.265492571105951E-01_wp, &
   & -5.860987383820147E-02_wp, -1.836785378448062E-03_wp,  8.619070045387887E-03_wp, &
   &  1.333887966099229E-02_wp, -6.941334314587246E-03_wp, -5.799877605709401E-03_wp, &
   &  4.187856249259576E-02_wp, -7.279888050845371E-02_wp,  3.860847673142320E-02_wp, &
   & -8.552249264904208E-04_wp, -5.744332214139890E-02_wp, -1.339644572136044E-03_wp, &
   & -1.560169469960605E-01_wp, -1.000250708768312E-02_wp, -2.811112246012069E-03_wp, &
   & -4.745814818565420E-03_wp, -1.759662411672108E-03_wp,  1.784668723847591E-01_wp, &
   &  1.958397655878534E-01_wp,  2.847711408940242E-01_wp, -1.457122936094604E-02_wp, &
   & -2.470383346368663E-03_wp, -8.917424258505349E-03_wp,  4.005588825677863E-04_wp, &
   & -3.553397090138397E-03_wp, -3.135596086621796E-01_wp, -3.263312566941264E-02_wp, &
   & -4.081635698877482E-02_wp,  2.989390691340027E-01_wp,  1.263999023166349E-03_wp, &
   &  2.087773084448162E-02_wp,  1.550653626890516E-02_wp,  1.210343840871198E-02_wp, &
   & -2.067220167151342E-01_wp,  1.903041744449289E-01_wp,  3.122766792850088E-01_wp, &
   &  1.626483199093528E+00_wp,  4.493699358442295E-03_wp,  1.594533242877436E-03_wp, &
   & -2.079105029582395E-02_wp, -4.407595175006052E-03_wp,  6.516354303963282E-02_wp, &
   &  1.762804123837060E-02_wp,  3.684016865220782E-02_wp, -1.045161900059688E-01_wp, &
   &  1.561297077817547E-02_wp,  7.973861302618956E-03_wp,  2.884145038176758E-02_wp, &
   &  7.138503141274094E-03_wp,  1.853806435677593E-01_wp,  8.358022480031102E-02_wp, &
   &  1.375238069651874E-01_wp,  5.478855100479740E-04_wp,  2.465131547229450E-02_wp, &
   &  9.788974676214313E-03_wp, -5.403053145017767E-04_wp, -3.141782913342259E-03_wp, &
   &  1.827837373453676E-01_wp, -3.848499253204109E-02_wp, -5.478951737105753E-02_wp, &
   &  3.490443238856631E-01_wp,  4.911426046033805E-02_wp,  1.605240717571781E-04_wp, &
   & -3.832627968016269E-02_wp, -4.101483768953129E-04_wp, -8.138703378017074E-02_wp, &
   &  7.217962268249751E-03_wp,  7.517433013925840E-03_wp,  8.482048844191342E-02_wp, &
   & -2.592953444016959E-02_wp, -6.247785722943803E-03_wp,  6.666155912529718E-03_wp, &
   & -2.139127442612924E-02_wp, -1.657269212141140E-01_wp,  1.803693884379085E-02_wp, &
   & -2.701109625402425E-02_wp, -4.241743832695642E-02_wp, -2.318708098773159E-02_wp, &
   & -3.408673525016986E-03_wp,  8.511548438474442E-03_wp,  2.312396198293607E-02_wp, &
   &  1.056851822576655E-02_wp,  2.409764949665099E-03_wp,  3.130074965933824E-02_wp, &
   & -6.533595999018363E-02_wp, -2.354356998774298E-02_wp, -5.057698564821952E-02_wp, &
   &  3.077026689983809E-02_wp,  4.493699358442292E-03_wp,  1.926833644485508E+00_wp, &
   & -2.523513143277224E-01_wp,  3.744524585342033E-01_wp, -2.079574806303811E-01_wp, &
   & -5.750832580695996E-03_wp, -2.148327938043032E-02_wp,  1.457785541285807E-02_wp, &
   & -1.219422583879158E-03_wp, -1.300885372020220E-01_wp, -1.099180630490448E-02_wp, &
   &  2.871982032181086E-03_wp, -2.243806020850800E-01_wp, -2.020289511613963E-02_wp, &
   & -1.664764211717569E-02_wp,  1.215162326489200E-02_wp, -5.254592228574507E-03_wp, &
   & -5.714568455347840E-03_wp,  5.391273638120669E-02_wp, -7.309206912203756E-02_wp, &
   &  9.713180632562089E-02_wp, -1.567940047099500E-02_wp, -3.261839753809036E-02_wp, &
   &  6.969332880562385E-03_wp, -1.141390689429487E-02_wp,  8.458327097158366E-03_wp, &
   &  1.027345385064401E-02_wp, -4.199367163358746E-02_wp,  9.553232181870224E-04_wp, &
   &  1.534008046510517E-02_wp, -1.465429531912288E-02_wp, -1.738383260717502E-02_wp, &
   &  3.388121259468538E-03_wp,  1.302343834342590E-02_wp,  3.551522894860388E-03_wp, &
   & -1.633220298493662E-02_wp,  2.444029883229076E-02_wp,  1.708975073338503E-01_wp, &
   &  2.893247956503804E-01_wp,  4.451157004621552E-01_wp, -2.120264316741839E-01_wp, &
   &  1.548467799425111E-02_wp,  3.368833142021845E-03_wp, -1.031655271816231E-02_wp, &
   & -1.042978442278318E-02_wp, -2.409212357888751E-03_wp, -6.684545227570542E-02_wp, &
   & -3.908937699891116E-02_wp, -1.791338449551878E-02_wp,  5.058114623212935E-02_wp, &
   &  5.275107182007691E-02_wp,  1.012784983664288E-02_wp,  1.594533242877449E-03_wp, &
   & -2.523513143277224E-01_wp,  1.369435097980728E+00_wp,  6.791475599129615E-02_wp, &
   & -1.906039112332074E-01_wp, -2.157597236931997E-02_wp,  7.969536622478790E-02_wp, &
   & -1.312748481039426E-03_wp,  2.088049510944025E-02_wp, -1.428839084912438E-02_wp, &
   & -1.082981672896379E-01_wp, -4.031454033397728E-02_wp, -2.774211684120608E-02_wp, &
   & -1.705124725543574E-02_wp,  8.207495225072624E-02_wp, -4.037645668190044E-03_wp, &
   & -9.690935545305562E-03_wp, -6.004869827815557E-02_wp, -2.513101784696622E-01_wp, &
   & -3.861292384178510E-01_wp,  4.879087861986919E-02_wp,  1.826470446091687E-02_wp, &
   &  1.235733857969077E-02_wp,  6.309438731553277E-03_wp,  1.206818938600337E-02_wp, &
   & -7.555890109641587E-03_wp,  4.621800031846540E-03_wp,  6.411417849114696E-02_wp, &
   &  8.567822520753905E-03_wp, -3.036148229123299E-02_wp, -1.038968907927841E-02_wp, &
   &  1.319690499959618E-02_wp,  1.836313813427145E-02_wp, -1.779073162236685E-02_wp, &
   & -5.251071076191872E-03_wp, -1.158231976716216E-02_wp,  7.502693923921305E-03_wp, &
   & -2.642221834488391E-01_wp,  4.456639329301224E-01_wp, -1.069518229891273E-01_wp, &
   &  3.057743526797048E-01_wp, -1.218955416095540E-02_wp, -1.663002589937036E-03_wp, &
   & -9.727541059417626E-03_wp,  5.007225909869392E-03_wp, -3.129192927315042E-02_wp, &
   & -3.909244019541337E-02_wp,  9.695864526332573E-03_wp,  3.718892145420852E-02_wp, &
   &  3.076438051088457E-02_wp, -1.012976688131399E-02_wp,  5.210857233250324E-02_wp, &
   & -2.079105029582395E-02_wp,  3.744524585342033E-01_wp,  6.791475599129615E-02_wp, &
   &  1.320615389099381E+00_wp,  3.129550954820600E-01_wp, -1.460280125534990E-02_wp, &
   &  1.385141888962091E-03_wp,  1.978366320534735E-02_wp,  1.552445640266042E-02_wp, &
   & -6.923747915993997E-03_wp, -4.200759717868002E-02_wp, -7.274181263876665E-02_wp, &
   &  3.768763400950512E-02_wp, -1.829700293572264E-03_wp, -8.665056645390731E-03_wp, &
   &  1.337301262745814E-02_wp, -6.957459169783020E-03_wp,  9.904109033335100E-02_wp, &
   & -3.957237012371094E-01_wp,  1.265482254234212E-01_wp, -5.857169747439377E-02_wp, &
   &  1.012750122886469E-02_wp, -1.395813206873921E-02_wp,  3.775597739913987E-02_wp, &
   &  1.952764309612394E-02_wp, -1.550649509124217E-01_wp, -1.358714909350574E-03_wp, &
   & -5.752151081594777E-02_wp, -6.397531873150938E-04_wp, -3.127961791029916E-01_wp, &
   &  3.264253319671220E-02_wp, -4.064972177795283E-02_wp,  2.974977636596358E-01_wp, &
   & -2.491510203795122E-03_wp,  9.054567357963395E-03_wp,  3.173781797194460E-04_wp, &
   & -3.596138573388682E-03_wp,  1.789397839834518E-01_wp, -1.958230124990306E-01_wp, &
   &  2.849175611106504E-01_wp, -1.436800564547772E-02_wp, -9.994284547877375E-03_wp, &
   &  2.746552100602511E-03_wp, -4.674478667837700E-03_wp, -1.762373847291776E-03_wp, &
   &  6.532425843818973E-02_wp, -1.790214766124291E-02_wp,  3.720983979051862E-02_wp, &
   & -1.041156167570835E-01_wp,  4.456717413048608E-03_wp, -1.534249792763398E-03_wp, &
   & -2.083376288009661E-02_wp, -4.407595175006052E-03_wp, -2.079574806303811E-01_wp, &
   & -1.906039112332074E-01_wp,  3.129550954820600E-01_wp,  1.626945363839795E+00_wp, &
   &  1.251441046804934E-03_wp, -2.086821921748747E-02_wp,  1.553401346564973E-02_wp, &
   &  1.209593579929269E-02_wp,  1.820968810900535E-01_wp,  3.858234803709026E-02_wp, &
   & -5.463397867212968E-02_wp,  3.483291098658333E-01_wp,  2.462892753347599E-02_wp, &
   & -9.904551241728204E-03_wp, -4.696345324793072E-04_wp, -3.155407564017262E-03_wp, &
   &  1.855366474611164E-01_wp, -8.373652887747721E-02_wp,  1.377131350035825E-01_wp, &
   &  1.105617959962876E-03_wp,  1.547051707812786E-02_wp, -8.020935835705831E-03_wp, &
   &  2.883807099290106E-02_wp,  7.037763979237069E-03_wp,  1.364518771015942E-04_wp, &
   &  4.913915356091953E-02_wp, -2.969041615684438E-04_wp, -3.794187581816817E-02_wp, &
   & -2.595677925547819E-02_wp, -6.197049531325822E-03_wp, -6.635914471490414E-03_wp, &
   &  2.142340779274882E-02_wp, -8.145251752221179E-02_wp,  7.266027835159628E-03_wp, &
   & -7.536133854468864E-03_wp, -8.475244399355290E-02_wp, -2.316632131853606E-02_wp, &
   & -3.247202360336833E-03_wp, -8.452122821130040E-03_wp, -2.318548203365187E-02_wp, &
   & -1.658753423440057E-01_wp,  1.795448754187384E-02_wp,  2.684783123649586E-02_wp, &
   &  4.204746720219069E-02_wp, -2.350629245771730E-02_wp, -5.060168129579690E-02_wp, &
   & -3.079888151813678E-02_wp, -4.518157760969753E-03_wp,  1.041404345108276E-02_wp, &
   &  2.367165320956746E-03_wp, -3.134950464921349E-02_wp,  6.516354303963282E-02_wp, &
   & -5.750832580695996E-03_wp, -2.157597236931998E-02_wp, -1.460280125534990E-02_wp, &
   &  1.251441046804934E-03_wp,  1.927020018257394E+00_wp, -2.523565268061509E-01_wp, &
   & -3.744924338344182E-01_wp,  2.072922152721698E-01_wp, -2.018650391307496E-02_wp, &
   & -1.666377900170123E-02_wp, -1.219942853830065E-02_wp,  5.272614688697778E-03_wp, &
   & -1.300685891544499E-01_wp, -1.102890178333774E-02_wp, -2.882560110566906E-03_wp, &
   &  2.243981124728023E-01_wp, -1.575723640445514E-02_wp, -3.268078629293569E-02_wp, &
   & -6.985144224365554E-03_wp,  1.142916324618238E-02_wp, -5.716655186579303E-03_wp, &
   &  5.392176570274911E-02_wp,  7.310876805353318E-02_wp, -9.781969897944398E-02_wp, &
   &  1.031567456603362E-02_wp,  8.212285660166100E-03_wp,  4.515508073346490E-04_wp, &
   & -4.210580981934167E-02_wp,  1.313400746515208E-02_wp,  3.576015982288077E-03_wp, &
   &  1.637587219690404E-02_wp, -2.440427953757240E-02_wp,  1.563599543080730E-02_wp, &
   & -1.466010502540641E-02_wp,  1.741201503734447E-02_wp, -3.181080000113450E-03_wp, &
   &  1.542197400806148E-02_wp,  3.297768879419099E-03_wp,  1.027468115563575E-02_wp, &
   &  1.053919970004938E-02_wp,  1.713553084114756E-01_wp,  2.894669651689678E-01_wp, &
   & -4.451203244006532E-01_wp,  2.121026114011793E-01_wp,  5.056989348941861E-02_wp, &
   &  5.269332533863068E-02_wp, -1.012824449158331E-02_wp, -1.532375845760117E-03_wp, &
   & -2.521274400848204E-03_wp, -6.641376893124294E-02_wp,  3.884260085549838E-02_wp, &
   &  1.762804123837061E-02_wp, -2.148327938043032E-02_wp,  7.969536622478790E-02_wp, &
   &  1.385141888962091E-03_wp, -2.086821921748747E-02_wp, -2.523565268061509E-01_wp, &
   &  1.368906368249732E+00_wp, -6.811835580875319E-02_wp,  1.903316662837480E-01_wp, &
   & -1.714111923267463E-02_wp,  8.217847282120570E-02_wp,  4.112494385897093E-03_wp, &
   &  9.651216831666802E-03_wp, -1.422754669816838E-02_wp, -1.081169841849596E-01_wp, &
   &  4.004869888694660E-02_wp,  2.756351257020310E-02_wp,  1.852396766443213E-02_wp, &
   &  1.238235927609199E-02_wp, -6.151559272887804E-03_wp, -1.225161878104451E-02_wp, &
   & -6.011355816262258E-02_wp, -2.515057885779451E-01_wp,  3.867521900533673E-01_wp, &
   & -4.899040059800329E-02_wp, -4.661771912720550E-03_wp,  7.164192470386323E-03_wp, &
   & -8.192072030219997E-03_wp, -6.422550781484228E-02_wp,  1.774030040201211E-02_wp, &
   &  5.300441441243784E-03_wp, -1.155542345146213E-02_wp,  7.502270778977416E-03_wp, &
   &  3.077920992193937E-02_wp,  1.037978655013428E-02_wp,  1.319859421964997E-02_wp, &
   &  1.863860377306688E-02_wp,  1.221152261208723E-02_wp,  1.593243135905389E-03_wp, &
   & -9.831466240989489E-03_wp,  4.979510908235438E-03_wp,  2.647943143328361E-01_wp, &
   & -4.456721790724243E-01_wp, -1.068394515688850E-01_wp,  3.061856776281458E-01_wp, &
   & -3.076405058003849E-02_wp,  1.027592095432595E-02_wp,  5.213976022273503E-02_wp, &
   & -2.085969784006980E-02_wp,  3.119459860731957E-02_wp,  3.887425976509527E-02_wp, &
   &  9.898472207251166E-03_wp,  3.684016865220782E-02_wp,  1.457785541285807E-02_wp, &
   & -1.312748481039419E-03_wp,  1.978366320534735E-02_wp,  1.553401346564973E-02_wp, &
   & -3.744924338344182E-01_wp, -6.811835580875318E-02_wp,  1.319829158812635E+00_wp, &
   &  3.123637026649840E-01_wp,  1.853550601840084E-03_wp,  8.664935894429557E-03_wp, &
   &  1.334342994033325E-02_wp, -6.929175554918170E-03_wp,  7.017335992862174E-03_wp, &
   &  4.180137185454618E-02_wp, -7.275537057343315E-02_wp,  3.732851220087915E-02_wp, &
   & -1.033697188617647E-02_wp,  1.394012851303677E-02_wp,  3.771933936637650E-02_wp, &
   &  1.963186604467970E-02_wp, -9.913262964485561E-02_wp,  3.963914766797381E-01_wp, &
   &  1.269458223872552E-01_wp, -5.914315248664746E-02_wp,  1.352730910145476E-03_wp, &
   &  1.548485620523062E-01_wp,  8.878425530952187E-04_wp,  5.660925794823445E-02_wp, &
   &  2.455335723719934E-03_wp, -9.017402911145283E-03_wp,  3.226255954096729E-04_wp, &
   & -3.559797879290408E-03_wp,  3.132114535967224E-01_wp, -3.279230862883733E-02_wp, &
   & -4.053426550408897E-02_wp,  2.979406044790362E-01_wp,  9.981223456534682E-03_wp, &
   & -2.842032462672789E-03_wp, -4.756019537487924E-03_wp, -1.753212643164470E-03_wp, &
   & -1.785910201043079E-01_wp,  1.955411855967179E-01_wp,  2.846339733760712E-01_wp, &
   & -1.332693544380373E-02_wp, -4.543925211347843E-03_wp,  1.450273323327533E-03_wp, &
   & -2.086277882728174E-02_wp, -4.433865201468268E-03_wp, -6.512169792188434E-02_wp, &
   &  1.810891325053175E-02_wp,  3.721641606431228E-02_wp, -1.045161900059688E-01_wp, &
   & -1.219422583879158E-03_wp,  2.088049510944025E-02_wp,  1.552445640266045E-02_wp, &
   &  1.209593579929269E-02_wp,  2.072922152721698E-01_wp,  1.903316662837480E-01_wp, &
   &  3.123637026649839E-01_wp,  1.627538602459993E+00_wp, -2.462085710491223E-02_wp, &
   &  9.861755228827623E-03_wp, -4.799783325223905E-04_wp, -3.170363180926034E-03_wp, &
   & -1.820692703569220E-01_wp, -3.843769378159749E-02_wp, -5.476937202521867E-02_wp, &
   &  3.482421834244690E-01_wp, -1.558443463552646E-02_wp,  8.070091651500730E-03_wp, &
   &  2.884833644800450E-02_wp,  7.139320580811570E-03_wp, -1.850282949772885E-01_wp, &
   &  8.381254486714156E-02_wp,  1.377989518292175E-01_wp,  2.035320217682310E-03_wp, &
   &  7.580502462269115E-02_wp,  1.437809962620081E-03_wp, -5.873365961043855E-02_wp, &
   &  5.475717228827796E-03_wp, -1.321201998006772E-01_wp, -2.412750523444410E-03_wp, &
   &  1.066750419360065E-02_wp,  6.772742090462043E-02_wp, -1.460952102501112E-02_wp, &
   & -1.254445802445833E-02_wp,  6.443821487081486E-03_wp, -5.286287347104757E-03_wp, &
   & -1.491410307302132E-01_wp, -1.774536610663534E-03_wp,  8.265497770558025E-03_wp, &
   & -1.703238437187407E-01_wp, -1.260646194170170E-02_wp, -1.000012702991848E-02_wp, &
   &  7.628098106089499E-03_wp,  5.250219250263501E-03_wp, -5.717935275029361E-03_wp, &
   &  6.005522562545247E-02_wp, -9.901469545667280E-02_wp, -1.855302611096530E-01_wp, &
   & -1.578594610730463E-02_wp, -1.849640619419218E-02_wp,  1.029362979057585E-02_wp, &
   &  1.561297077817547E-02_wp, -1.300885372020220E-01_wp, -1.428839084912438E-02_wp, &
   & -6.923747915993997E-03_wp,  1.820968810900535E-01_wp, -2.018650391307496E-02_wp, &
   & -1.714111923267464E-02_wp,  1.853550601840084E-03_wp, -2.462085710491225E-02_wp, &
   &  1.834500576467604E+00_wp, -2.328380235053190E-01_wp,  3.700742976438847E-01_wp, &
   & -1.645813292093071E-01_wp, -1.482033905084581E-02_wp,  2.965917328149744E-03_wp, &
   & -2.277210531366966E-04_wp,  2.065108905999737E-02_wp,  5.520562654251384E-02_wp, &
   & -1.290129365990534E-02_wp,  8.686895758138347E-02_wp,  4.638743655401732E-02_wp, &
   & -2.770509482107263E-02_wp, -3.450542723335484E-02_wp,  9.280670468506228E-03_wp, &
   & -1.063485190409389E-02_wp,  3.423914070860962E-02_wp, -5.220645740682702E-03_wp, &
   & -2.997885454140316E-02_wp, -1.180731862149381E-02_wp,  1.488951503221468E-01_wp, &
   &  1.499577451346859E-01_wp,  3.907479332333729E-01_wp,  1.487991423992202E-01_wp, &
   &  2.147256620212746E-03_wp,  1.418234526306885E-02_wp, -1.282593655401865E-02_wp, &
   &  2.474080389321252E-03_wp, -2.671754073244886E-02_wp, -1.683422150244912E-02_wp, &
   & -2.559730740604082E-02_wp, -4.212739390044037E-02_wp,  5.592320105044995E-03_wp, &
   &  1.970151246932410E-02_wp, -2.527738044445453E-02_wp,  9.621043622666831E-03_wp, &
   & -5.390317807562760E-02_wp, -2.513313120542579E-01_wp, -3.957205218392055E-01_wp, &
   & -8.374764020602153E-02_wp,  3.269703607432843E-02_wp,  1.237647302717487E-02_wp, &
   &  1.395801224594882E-02_wp,  7.973861302618956E-03_wp, -1.099180630490448E-02_wp, &
   & -1.082981672896379E-01_wp, -4.200759717868002E-02_wp,  3.858234803709028E-02_wp, &
   & -1.666377900170123E-02_wp,  8.217847282120570E-02_wp,  8.664935894429550E-03_wp, &
   &  9.861755228827609E-03_wp, -2.328380235053190E-01_wp,  1.575497944172972E+00_wp, &
   &  1.839537781589221E-01_wp,  1.362040843363563E-01_wp,  2.944664462031690E-03_wp, &
   &  3.323176378487962E-02_wp,  3.788817515003366E-03_wp, -6.942485255779902E-03_wp, &
   &  1.291021010402296E-02_wp, -3.077184722091095E-02_wp,  9.352556651385235E-03_wp, &
   &  1.910357442093814E-02_wp,  3.463671161763840E-02_wp, -6.595425127400040E-03_wp, &
   &  1.505601503706575E-02_wp, -5.994551409174813E-04_wp, -4.526923976551495E-02_wp, &
   &  5.464860841707856E-03_wp,  5.589274797182139E-02_wp, -4.014856647449936E-03_wp, &
   & -2.510248276258650E-01_wp,  3.926320642675794E-01_wp, -2.335759858806779E-01_wp, &
   & -2.451460387601280E-01_wp, -1.966606484542853E-03_wp,  1.525759156952471E-03_wp, &
   & -4.838100087988706E-04_wp,  9.464848782751618E-03_wp,  6.101614619024202E-02_wp, &
   & -8.703708739375223E-03_wp,  5.061542100640776E-03_wp,  6.303042784252560E-02_wp, &
   & -3.104053425941847E-03_wp,  1.760256337875835E-03_wp, -7.398974178769744E-03_wp, &
   & -4.610957486082562E-03_wp,  7.309749637482553E-02_wp, -3.861543712368549E-01_wp, &
   &  1.265239701221103E-01_wp,  1.376598521306487E-01_wp,  6.964522965623769E-03_wp, &
   & -6.178139001288267E-03_wp,  3.774051610790029E-02_wp,  2.884145038176758E-02_wp, &
   &  2.871982032181086E-03_wp, -4.031454033397728E-02_wp, -7.274181263876665E-02_wp, &
   & -5.463397867212968E-02_wp, -1.219942853830065E-02_wp,  4.112494385897093E-03_wp, &
   &  1.334342994033325E-02_wp, -4.799783325223905E-04_wp,  3.700742976438847E-01_wp, &
   &  1.839537781589221E-01_wp,  1.374558780674070E+00_wp, -2.251410320954913E-01_wp, &
   &  2.341240897192308E-04_wp, -3.831522060671021E-03_wp,  6.668809435068521E-03_wp, &
   &  4.415191249561461E-03_wp, -8.684431372327756E-02_wp,  9.349134536507570E-03_wp, &
   & -9.924441515200284E-03_wp, -8.271099172811339E-03_wp,  9.219977888707088E-03_wp, &
   & -1.510790172846692E-02_wp,  4.870582257507241E-02_wp, -3.017551899557763E-02_wp, &
   &  5.290300153629187E-01_wp, -5.022889917614801E-03_wp,  2.873288137202622E-02_wp, &
   & -5.311396711876894E-03_wp, -1.534607830666470E-01_wp,  1.033017651324187E-01_wp, &
   & -1.745947607698172E-01_wp, -5.172003090114060E-03_wp, -2.359994158908536E-04_wp, &
   & -8.086516903264363E-03_wp,  3.409247581638760E-03_wp,  6.742819824984703E-04_wp, &
   &  2.170839367149032E-01_wp, -4.024801588740388E-02_wp,  6.058992671227835E-02_wp, &
   &  1.939408725198874E-01_wp,  2.438377742502520E-03_wp,  3.012932969958350E-03_wp, &
   &  3.985294758549605E-03_wp,  8.077821380194812E-03_wp, -9.714300020215932E-02_wp, &
   &  4.879030348393663E-02_wp, -5.861095255407231E-02_wp,  1.101342225395888E-03_wp, &
   & -1.142018525754264E-02_wp, -1.225637791523566E-02_wp,  1.960703470239653E-02_wp, &
   &  7.138503141274080E-03_wp, -2.243806020850800E-01_wp, -2.774211684120608E-02_wp, &
   &  3.768763400950512E-02_wp,  3.483291098658333E-01_wp,  5.272614688697777E-03_wp, &
   &  9.651216831666802E-03_wp, -6.929175554918170E-03_wp, -3.170363180926034E-03_wp, &
   & -1.645813292093071E-01_wp,  1.362040843363563E-01_wp, -2.251410320954913E-01_wp, &
   &  1.274331705001594E+00_wp, -2.063632970626589E-02_wp,  6.914349002436846E-03_wp, &
   &  4.422524409965087E-03_wp,  3.470498377624705E-03_wp, -4.638404157466271E-02_wp, &
   &  1.909759554937811E-02_wp, -8.277297221126335E-03_wp, -3.474715691841671E-02_wp, &
   & -1.051823147338514E-02_wp,  6.739784781324351E-04_wp, -3.017793667473354E-02_wp, &
   &  7.107277939760917E-03_wp,  1.454694527520587E-03_wp,  7.579206001018166E-02_wp, &
   &  5.408112689239400E-03_wp, -5.839639244371739E-02_wp, -1.460308024330004E-02_wp, &
   & -1.256618474760451E-02_wp, -6.462818904493103E-03_wp,  5.290270182577418E-03_wp, &
   & -1.323594660617543E-01_wp, -2.121388674350319E-03_wp, -1.063982081247404E-02_wp, &
   & -6.784656854515582E-02_wp, -1.267048579151313E-02_wp, -1.000069826473252E-02_wp, &
   & -7.583582685916045E-03_wp, -5.268995145825184E-03_wp, -1.491780784715833E-01_wp, &
   & -1.687030687679118E-03_wp, -8.159842431310724E-03_wp,  1.706126424970112E-01_wp, &
   & -1.573785119430015E-02_wp, -1.852974369621771E-02_wp, -1.029243857275900E-02_wp, &
   & -1.555877947121680E-02_wp, -5.688597060383519E-03_wp,  5.992976277236097E-02_wp, &
   &  9.891979102940525E-02_wp,  1.853806435677593E-01_wp, -2.020289511613963E-02_wp, &
   & -1.705124725543574E-02_wp, -1.829700293572267E-03_wp,  2.462892753347599E-02_wp, &
   & -1.300685891544499E-01_wp, -1.422754669816838E-02_wp,  7.017335992862177E-03_wp, &
   & -1.820692703569220E-01_wp, -1.482033905084581E-02_wp,  2.944664462031687E-03_wp, &
   &  2.341240897192343E-04_wp, -2.063632970626589E-02_wp,  1.834803540998220E+00_wp, &
   & -2.325650932054020E-01_wp, -3.700591009079603E-01_wp,  1.649016559026844E-01_wp, &
   & -2.762248232292545E-02_wp, -3.455460835604964E-02_wp, -9.232478319507905E-03_wp, &
   &  1.052209109133401E-02_wp,  5.491827716436429E-02_wp, -1.283899963125250E-02_wp, &
   & -8.668269080210475E-02_wp, -4.604827592579480E-02_wp, -5.230543099499541E-03_wp, &
   &  3.427650757146209E-02_wp, -1.169935146473086E-02_wp, -3.020781403039341E-02_wp, &
   &  2.145605836418894E-03_wp,  1.419760878664407E-02_wp,  1.281852403078568E-02_wp, &
   & -2.483563415016234E-03_wp,  1.487433493467488E-01_wp,  1.500688127451139E-01_wp, &
   & -3.905781596073190E-01_wp, -1.485893963175227E-01_wp,  5.701640309369222E-03_wp, &
   &  1.987396831224482E-02_wp,  2.536263019625774E-02_wp, -9.590785761891107E-03_wp, &
   & -2.663163457527537E-02_wp, -1.700436400697725E-02_wp,  2.564908723643024E-02_wp, &
   &  4.204645221277876E-02_wp,  3.262800478449652E-02_wp,  1.232878797800032E-02_wp, &
   & -1.401008290642191E-02_wp, -8.035277723205586E-03_wp, -5.384009938581511E-02_wp, &
   & -2.511369516436848E-01_wp,  3.953897623034650E-01_wp,  8.358022480031102E-02_wp, &
   & -1.664764211717569E-02_wp,  8.207495225072622E-02_wp, -8.665056645390731E-03_wp, &
   & -9.904551241728204E-03_wp, -1.102890178333774E-02_wp, -1.081169841849596E-01_wp, &
   &  4.180137185454618E-02_wp, -3.843769378159749E-02_wp,  2.965917328149740E-03_wp, &
   &  3.323176378487962E-02_wp, -3.831522060671021E-03_wp,  6.914349002436850E-03_wp, &
   & -2.325650932054020E-01_wp,  1.576019052337825E+00_wp, -1.839774554074385E-01_wp, &
   & -1.360141528537016E-01_wp,  3.455938062385694E-02_wp, -6.618405144708205E-03_wp, &
   & -1.510475926221091E-02_wp,  6.641684508544819E-04_wp,  1.299078842167502E-02_wp, &
   & -3.074703975247960E-02_wp, -9.303734205697342E-03_wp, -1.855138818545762E-02_wp, &
   & -5.475539339841637E-03_wp,  4.529702101290851E-02_wp,  4.068025141748031E-03_wp, &
   & -5.626759640739265E-02_wp,  1.982681150420149E-03_wp, -1.513744763688584E-03_wp, &
   & -4.627416287810797E-04_wp,  9.466882832942922E-03_wp,  2.512430203608930E-01_wp, &
   & -3.926375660752298E-01_wp, -2.338639896244302E-01_wp, -2.452245826391641E-01_wp, &
   &  3.116726806188214E-03_wp, -1.826680563705499E-03_wp, -7.410634617386340E-03_wp, &
   & -4.626907472414760E-03_wp, -6.121419307829654E-02_wp,  8.892644060682691E-03_wp, &
   &  4.937990448470266E-03_wp,  6.342946382614620E-02_wp, -6.973219408918997E-03_wp, &
   &  6.188808398671794E-03_wp,  3.775030908808071E-02_wp,  2.886706449458729E-02_wp, &
   & -7.288773528300484E-02_wp,  3.859972450964884E-01_wp,  1.265492571105951E-01_wp, &
   &  1.375238069651874E-01_wp,  1.215162326489200E-02_wp, -4.037645668190044E-03_wp, &
   &  1.337301262745814E-02_wp, -4.696345324793072E-04_wp, -2.882560110566906E-03_wp, &
   &  4.004869888694660E-02_wp, -7.275537057343313E-02_wp, -5.476937202521867E-02_wp, &
   & -2.277210531366966E-04_wp,  3.788817515003366E-03_wp,  6.668809435068521E-03_wp, &
   &  4.422524409965094E-03_wp, -3.700591009079603E-01_wp, -1.839774554074385E-01_wp, &
   &  1.374236415345905E+00_wp, -2.251858314954866E-01_wp, -9.236323039855960E-03_wp, &
   &  1.507883069385210E-02_wp,  4.867441834435678E-02_wp, -3.019206203668821E-02_wp, &
   &  8.704559606665746E-02_wp, -9.398686908794328E-03_wp, -9.857900553131523E-03_wp, &
   & -7.705620697797476E-03_wp,  5.029086448775175E-03_wp, -5.291873073042418E-01_wp, &
   &  5.322386258012152E-03_wp, -2.956141066515159E-02_wp,  2.367182058618279E-04_wp, &
   &  8.101542765591399E-03_wp,  3.411384254096665E-03_wp,  6.822251226400594E-04_wp, &
   &  1.533997236662861E-01_wp, -1.032282700263639E-01_wp, -1.745248313432676E-01_wp, &
   & -5.214816827278326E-03_wp, -2.476142530256957E-03_wp, -2.969229187384275E-03_wp, &
   &  4.014880863428783E-03_wp,  8.044404888566635E-03_wp, -2.169361995806018E-01_wp, &
   &  4.026448777005662E-02_wp,  6.064273938241257E-02_wp,  1.939857159631502E-01_wp, &
   &  1.143549362508695E-02_wp,  1.223242374577508E-02_wp,  1.963260925995177E-02_wp, &
   &  7.146330771728576E-03_wp,  9.743519498628025E-02_wp, -4.872401455665089E-02_wp, &
   & -5.860987383820147E-02_wp,  5.478855100479740E-04_wp, -5.254592228574507E-03_wp, &
   & -9.690935545305555E-03_wp, -6.957459169783013E-03_wp, -3.155407564017262E-03_wp, &
   &  2.243981124728023E-01_wp,  2.756351257020310E-02_wp,  3.732851220087915E-02_wp, &
   &  3.482421834244690E-01_wp,  2.065108905999738E-02_wp, -6.942485255779899E-03_wp, &
   &  4.415191249561461E-03_wp,  3.470498377624705E-03_wp,  1.649016559026844E-01_wp, &
   & -1.360141528537016E-01_wp, -2.251858314954866E-01_wp,  1.274144072399136E+00_wp, &
   &  1.051870076218589E-02_wp, -6.703544177567994E-04_wp, -3.021956250439821E-02_wp, &
   &  7.184346264968118E-03_wp,  4.700936232345972E-02_wp, -1.907704158349379E-02_wp, &
   & -8.244670798347911E-03_wp, -3.335941760308349E-02_wp, -5.873750796922662E-02_wp, &
   &  5.398541212814784E-03_wp,  7.579801364543501E-02_wp,  1.510806635622226E-03_wp, &
   & -1.491414911419830E-01_wp,  1.767956891538433E-03_wp, -8.267152240045386E-03_wp, &
   &  1.703210708227039E-01_wp, -1.266821157991562E-02_wp,  1.001834442550433E-02_wp, &
   & -7.608286765451494E-03_wp, -5.256371142047481E-03_wp, -1.321078178110328E-01_wp, &
   &  2.410365605460931E-03_wp, -1.066620604694135E-02_wp, -6.771802335365747E-02_wp, &
   & -1.457201130417061E-02_wp,  1.253150227241959E-02_wp, -6.439218353030748E-03_wp, &
   &  5.312023266393023E-03_wp, -1.300895199971275E-01_wp,  1.429155078064440E-02_wp, &
   &  6.933316499867500E-03_wp, -1.820906817729780E-01_wp, -2.021590405397043E-02_wp, &
   &  1.710229838831284E-02_wp, -1.836785378448062E-03_wp,  2.465131547229449E-02_wp, &
   & -5.714568455347841E-03_wp, -6.004869827815557E-02_wp,  9.904109033335100E-02_wp, &
   &  1.855366474611164E-01_wp, -1.575723640445514E-02_wp,  1.852396766443213E-02_wp, &
   & -1.033697188617646E-02_wp, -1.558443463552646E-02_wp,  5.520562654251384E-02_wp, &
   &  1.291021010402296E-02_wp, -8.684431372327756E-02_wp, -4.638404157466271E-02_wp, &
   & -2.762248232292545E-02_wp,  3.455938062385694E-02_wp, -9.236323039855960E-03_wp, &
   &  1.051870076218589E-02_wp,  1.834510243202380E+00_wp,  2.328151835234852E-01_wp, &
   & -3.700575059576854E-01_wp,  1.646065678500633E-01_wp, -1.488508890941887E-02_wp, &
   & -2.930551134459511E-03_wp,  2.414861299464589E-04_wp, -2.072585863659663E-02_wp, &
   &  2.997882420757660E-02_wp,  1.170724185580760E-02_wp, -3.425062092156357E-02_wp, &
   &  4.971150257396431E-03_wp,  2.672404872968473E-02_wp, -1.683065106280155E-02_wp, &
   & -2.559729433992322E-02_wp, -4.213592394024056E-02_wp, -5.673665778168770E-03_wp, &
   &  1.979827992580431E-02_wp, -2.533577437102717E-02_wp,  9.637880193722827E-03_wp, &
   & -1.488842647684905E-01_wp,  1.499476098159285E-01_wp,  3.907223175034765E-01_wp, &
   &  1.487610897931907E-01_wp, -2.027423210222611E-03_wp,  1.406281643489339E-02_wp, &
   & -1.275554334950709E-02_wp,  2.457308742080725E-03_wp,  1.099100123152548E-02_wp, &
   & -1.083191611256763E-01_wp, -4.201691464225484E-02_wp,  3.857194786084731E-02_wp, &
   &  1.672461362425607E-02_wp,  8.217955684208357E-02_wp,  8.619070045387887E-03_wp, &
   &  9.788974676214313E-03_wp,  5.391273638120669E-02_wp, -2.513101784696622E-01_wp, &
   & -3.957237012371094E-01_wp, -8.373652887747721E-02_wp, -3.268078629293569E-02_wp, &
   &  1.238235927609200E-02_wp,  1.394012851303676E-02_wp,  8.070091651500730E-03_wp, &
   & -1.290129365990533E-02_wp, -3.077184722091095E-02_wp,  9.349134536507570E-03_wp, &
   &  1.909759554937811E-02_wp, -3.455460835604964E-02_wp, -6.618405144708205E-03_wp, &
   &  1.507883069385210E-02_wp, -6.703544177567994E-04_wp,  2.328151835234852E-01_wp, &
   &  1.575586046478149E+00_wp,  1.839754782159632E-01_wp,  1.361947607761956E-01_wp, &
   & -2.835842105689206E-03_wp,  3.330132930893714E-02_wp,  3.839979065064121E-03_wp, &
   & -6.632617820337198E-03_wp, -5.591177901303645E-02_wp,  4.068259828654927E-03_wp, &
   &  4.526598713611546E-02_wp, -5.627628235781873E-03_wp, -6.106109832140092E-02_wp, &
   & -8.700242824565430E-03_wp,  5.062772494027674E-03_wp,  6.307130221910287E-02_wp, &
   &  3.101807402292274E-03_wp,  1.853235953651607E-03_wp, -7.433462681310299E-03_wp, &
   & -4.589717515932207E-03_wp,  2.510445353237979E-01_wp,  3.926132423729201E-01_wp, &
   & -2.336135451174082E-01_wp, -2.451443254802192E-01_wp,  1.980517257124989E-03_wp, &
   &  1.440319008298704E-03_wp, -4.352368225150865E-04_wp,  9.423759394375709E-03_wp, &
   & -2.864002799454561E-03_wp, -4.031613472964875E-02_wp, -7.274251842956297E-02_wp, &
   & -5.465485333034970E-02_wp,  1.223352579956685E-02_wp,  4.089208890809513E-03_wp, &
   &  1.333887966099229E-02_wp, -5.403053145017489E-04_wp, -7.309206912203756E-02_wp, &
   & -3.861292384178510E-01_wp,  1.265482254234212E-01_wp,  1.377131350035825E-01_wp, &
   & -6.985144224365554E-03_wp, -6.151559272887810E-03_wp,  3.771933936637650E-02_wp, &
   &  2.884833644800447E-02_wp,  8.686895758138347E-02_wp,  9.352556651385235E-03_wp, &
   & -9.924441515200284E-03_wp, -8.277297221126328E-03_wp, -9.232478319507902E-03_wp, &
   & -1.510475926221091E-02_wp,  4.867441834435678E-02_wp, -3.021956250439822E-02_wp, &
   & -3.700575059576854E-01_wp,  1.839754782159632E-01_wp,  1.374534701408154E+00_wp, &
   & -2.251451830774975E-01_wp, -1.404128029104043E-04_wp, -3.815997945520934E-03_wp, &
   &  6.658412764855771E-03_wp,  4.558105366480492E-03_wp, -2.873784305744250E-02_wp, &
   &  5.306119173344457E-03_wp, -5.290432544172488E-01_wp,  5.089000567812710E-03_wp, &
   & -2.170756894938684E-01_wp, -4.025253508741049E-02_wp,  6.058348225014799E-02_wp, &
   &  1.939267990732981E-01_wp, -2.473754141221987E-03_wp,  2.987968882560821E-03_wp, &
   &  4.005414443716922E-03_wp,  8.060262539910920E-03_wp,  1.534613979273354E-01_wp, &
   &  1.032751244669945E-01_wp, -1.745923449463822E-01_wp, -5.173852886849753E-03_wp, &
   &  2.242175637279242E-04_wp, -8.079919789563848E-03_wp,  3.404030448830548E-03_wp, &
   &  6.756893186891516E-04_wp,  2.243791573011786E-01_wp, -2.774561432925469E-02_wp, &
   &  3.770254238684999E-02_wp,  3.483625814784412E-01_wp, -5.276732047483599E-03_wp, &
   &  9.652999950917278E-03_wp, -6.941334314587246E-03_wp, -3.141782913342245E-03_wp, &
   &  9.713180632562089E-02_wp,  4.879087861986919E-02_wp, -5.857169747439377E-02_wp, &
   &  1.105617959962876E-03_wp,  1.142916324618238E-02_wp, -1.225161878104451E-02_wp, &
   &  1.963186604467970E-02_wp,  7.139320580811570E-03_wp,  4.638743655401732E-02_wp, &
   &  1.910357442093815E-02_wp, -8.271099172811339E-03_wp, -3.474715691841671E-02_wp, &
   &  1.052209109133401E-02_wp,  6.641684508544854E-04_wp, -3.019206203668821E-02_wp, &
   &  7.184346264968118E-03_wp,  1.646065678500633E-01_wp,  1.361947607761956E-01_wp, &
   & -2.251451830774975E-01_wp,  1.274329978592291E+00_wp,  2.060787791617670E-02_wp, &
   &  6.886547060421595E-03_wp,  4.422047794420313E-03_wp,  3.483173458245455E-03_wp, &
   &  5.401793605232210E-03_wp, -5.866460440253323E-02_wp,  1.442174296498074E-03_wp, &
   &  7.647175198668643E-02_wp, -1.264404032600393E-02_wp,  1.000643989248778E-02_wp, &
   &  7.642773678844308E-03_wp,  5.294122943726928E-03_wp, -1.480200394334560E-01_wp, &
   &  1.700385822401850E-03_wp,  8.116145381064886E-03_wp, -1.691793086208223E-01_wp, &
   & -1.460755041629342E-02_wp,  1.246165939538359E-02_wp,  6.352354891197232E-03_wp, &
   & -5.332602643252117E-03_wp, -1.310292759259390E-01_wp,  2.288466861541040E-03_wp, &
   &  1.048516538151298E-02_wp,  6.511143463338417E-02_wp, -2.023228360698578E-02_wp, &
   &  1.715256929930816E-02_wp,  1.865495047479869E-03_wp, -2.462283074712185E-02_wp, &
   & -1.318657238440606E-01_wp,  1.484411960351194E-02_wp, -5.799877605709401E-03_wp, &
   &  1.827837373453675E-01_wp, -1.567940047099500E-02_wp,  1.826470446091687E-02_wp, &
   &  1.012750122886468E-02_wp,  1.547051707812786E-02_wp, -5.716655186579303E-03_wp, &
   & -6.011355816262258E-02_wp, -9.913262964485561E-02_wp, -1.850282949772885E-01_wp, &
   & -2.770509482107263E-02_wp,  3.463671161763840E-02_wp,  9.219977888707088E-03_wp, &
   & -1.051823147338514E-02_wp,  5.491827716436429E-02_wp,  1.299078842167502E-02_wp, &
   &  8.704559606665746E-02_wp,  4.700936232345972E-02_wp, -1.488508890941887E-02_wp, &
   & -2.835842105689206E-03_wp, -1.404128029104078E-04_wp,  2.060787791617670E-02_wp, &
   &  1.834411652882381E+00_wp,  2.325373376080316E-01_wp,  3.702278907330752E-01_wp, &
   & -1.660517967600043E-01_wp,  1.165434529800888E-02_wp,  2.995383714867841E-02_wp, &
   &  5.218762433237576E-03_wp, -3.340646505662848E-02_wp, -5.711514899287208E-03_wp, &
   &  1.982971406572112E-02_wp,  2.531309065105236E-02_wp, -9.600813313738004E-03_wp, &
   &  2.682219595221672E-02_wp, -1.661051130919642E-02_wp,  2.544418327570327E-02_wp, &
   &  4.207039178782483E-02_wp, -2.136201982272271E-03_wp,  1.420831397488669E-02_wp, &
   &  1.283705618323999E-02_wp, -2.446411327017794E-03_wp, -1.497798779426327E-01_wp, &
   &  1.501352108339646E-01_wp, -3.912198384430868E-01_wp, -1.483828137330170E-01_wp, &
   &  1.666726096580421E-02_wp,  8.203746951009575E-02_wp, -8.635557223496085E-03_wp, &
   & -9.868078060312466E-03_wp,  1.087088887548559E-02_wp, -1.082103747147362E-01_wp, &
   &  4.187856249259576E-02_wp, -3.848499253204109E-02_wp, -3.261839753809036E-02_wp, &
   &  1.235733857969077E-02_wp, -1.395813206873921E-02_wp, -8.020935835705831E-03_wp, &
   &  5.392176570274911E-02_wp, -2.515057885779451E-01_wp,  3.963914766797381E-01_wp, &
   &  8.381254486714157E-02_wp, -3.450542723335484E-02_wp, -6.595425127400040E-03_wp, &
   & -1.510790172846692E-02_wp,  6.739784781324316E-04_wp, -1.283899963125250E-02_wp, &
   & -3.074703975247960E-02_wp, -9.398686908794328E-03_wp, -1.907704158349379E-02_wp, &
   & -2.930551134459515E-03_wp,  3.330132930893714E-02_wp, -3.815997945520934E-03_wp, &
   &  6.886547060421598E-03_wp,  2.325373376080316E-01_wp,  1.574752228399172E+00_wp, &
   & -1.842076159964946E-01_wp, -1.359585205036384E-01_wp, -4.088414634391401E-03_wp, &
   &  5.588706921176288E-02_wp,  5.494918943424746E-03_wp, -4.523805977454380E-02_wp, &
   & -3.087023996236248E-03_wp, -1.893954077802771E-03_wp, -7.459720442475663E-03_wp, &
   & -4.616744038439015E-03_wp,  6.117560019213610E-02_wp,  8.539519630381751E-03_wp, &
   &  5.209091343999873E-03_wp,  6.303054044855710E-02_wp, -1.997947093469191E-03_wp, &
   & -1.507519814669938E-03_wp, -4.543834522675482E-04_wp,  9.479628340630968E-03_wp, &
   & -2.524053952636748E-01_wp, -3.931103870497424E-01_wp, -2.339238248737927E-01_wp, &
   & -2.444355309205807E-01_wp, -1.215357414037543E-02_wp, -4.044006799140985E-03_wp, &
   &  1.338003834431531E-02_wp, -4.459001831734477E-04_wp,  2.919030379864969E-03_wp, &
   &  4.020646820081600E-02_wp, -7.279888050845371E-02_wp, -5.478951737105753E-02_wp, &
   &  6.969332880562385E-03_wp,  6.309438731553277E-03_wp,  3.775597739913987E-02_wp, &
   &  2.883807099290103E-02_wp,  7.310876805353318E-02_wp,  3.867521900533673E-01_wp, &
   &  1.269458223872552E-01_wp,  1.377989518292176E-01_wp,  9.280670468506228E-03_wp, &
   &  1.505601503706575E-02_wp,  4.870582257507241E-02_wp, -3.017793667473354E-02_wp, &
   & -8.668269080210475E-02_wp, -9.303734205697335E-03_wp, -9.857900553131536E-03_wp, &
   & -8.244670798347911E-03_wp,  2.414861299464555E-04_wp,  3.839979065064121E-03_wp, &
   &  6.658412764855771E-03_wp,  4.422047794420313E-03_wp,  3.702278907330751E-01_wp, &
   & -1.842076159964946E-01_wp,  1.373453579338022E+00_wp, -2.240452486999068E-01_wp, &
   & -5.238704930788645E-03_wp,  2.859858066898498E-02_wp, -4.977729138402988E-03_wp, &
   &  5.286168946426076E-01_wp,  2.538143068475399E-03_wp, -3.014670166330845E-03_wp, &
   &  4.013830365771267E-03_wp,  8.115573032660264E-03_wp,  2.169113743768337E-01_wp, &
   &  3.987033032478870E-02_wp,  6.041076413473551E-02_wp,  1.943019231543751E-01_wp, &
   & -2.491047885149200E-04_wp,  7.879338749816587E-03_wp,  3.266875793296566E-03_wp, &
   &  6.297357615509189E-04_wp, -1.496382287741379E-01_wp, -1.029729431720197E-01_wp, &
   & -1.738028570516824E-01_wp, -6.722946300012508E-03_wp,  5.164434699717015E-03_wp, &
   & -9.652428849881690E-03_wp, -6.908828442012724E-03_wp, -3.127256758552313E-03_wp, &
   & -2.260728331184097E-01_wp,  2.841422012598828E-02_wp,  3.860847673142320E-02_wp, &
   &  3.490443238856631E-01_wp, -1.141390689429487E-02_wp,  1.206818938600336E-02_wp, &
   &  1.952764309612394E-02_wp,  7.037763979237069E-03_wp, -9.781969897944398E-02_wp, &
   & -4.899040059800329E-02_wp, -5.914315248664746E-02_wp,  2.035320217682296E-03_wp, &
   & -1.063485190409389E-02_wp, -5.994551409174813E-04_wp, -3.017551899557763E-02_wp, &
   &  7.107277939760924E-03_wp, -4.604827592579479E-02_wp, -1.855138818545762E-02_wp, &
   & -7.705620697797476E-03_wp, -3.335941760308349E-02_wp, -2.072585863659662E-02_wp, &
   & -6.632617820337198E-03_wp,  4.558105366480485E-03_wp,  3.483173458245455E-03_wp, &
   & -1.660517967600043E-01_wp, -1.359585205036384E-01_wp, -2.240452486999068E-01_wp, &
   &  1.276983551217846E+00_wp], &
   &  shape(density))

   call get_structure(mol, "X23", "oxacb")

   call test_num_op_grad(error, mol, density, make_exchange_gxtb, thr_in=thr2)

end subroutine test_op_g_fock_oxacb



subroutine test_op_g_fock_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(48,48,1) = reshape([ &
   &  5.926107199394257E-01_wp,  6.313752291162211E-04_wp,  6.313752291159144E-04_wp, &
   &  6.313752291162483E-04_wp, -1.002169991707603E-02_wp, -1.652680392332771E-05_wp, &
   & -2.796388824430549E-05_wp, -2.111433219217073E-05_wp, -1.002169991707614E-02_wp, &
   & -2.111433219181250E-05_wp, -1.652680392330776E-05_wp, -2.796388824409472E-05_wp, &
   & -1.002169991707632E-02_wp, -2.796388824422173E-05_wp, -2.111433219224735E-05_wp, &
   & -1.652680392333367E-05_wp,  2.675463048516357E-01_wp, -1.541417593578825E-01_wp, &
   & -1.541417593578830E-01_wp, -1.541417593578828E-01_wp, -8.075824133965537E-03_wp, &
   &  4.405532199375228E-03_wp,  1.010542468054945E-03_wp, -6.422706566373671E-03_wp, &
   & -8.067347093899502E-03_wp,  6.427870946056665E-03_wp, -4.329167530505221E-03_wp, &
   & -1.068428092528660E-03_wp, -8.067347093899963E-03_wp, -1.068428092528640E-03_wp, &
   &  6.427870946056199E-03_wp, -4.329167530505186E-03_wp,  2.667631089770241E-01_wp, &
   &  1.540440853286230E-01_wp,  1.540440853286227E-01_wp,  1.540440853286232E-01_wp, &
   & -8.067347093899700E-03_wp, -4.329167530504967E-03_wp, -1.068428092529037E-03_wp, &
   &  6.427870946056195E-03_wp, -8.075824133965659E-03_wp, -6.422706566372709E-03_wp, &
   &  4.405532199374985E-03_wp,  1.010542468055552E-03_wp, -8.075824133965872E-03_wp, &
   &  1.010542468055235E-03_wp, -6.422706566373787E-03_wp,  4.405532199374922E-03_wp, &
   &  6.313752291162211E-04_wp,  3.762713946849557E-01_wp,  4.065079782192155E-02_wp, &
   &  4.065079782192278E-02_wp, -1.616106489918758E-05_wp,  1.725572476348868E-03_wp, &
   & -3.724237379732360E-03_wp, -2.195886734634207E-03_wp, -1.883335093316854E-05_wp, &
   &  4.332437613238558E-03_wp, -2.195319400570014E-03_wp,  1.484721906251033E-03_wp, &
   &  2.694325003484910E-05_wp, -1.426698129811196E-04_wp, -1.484219597730527E-03_wp, &
   &  3.723791929503827E-03_wp,  1.320963754698337E-01_wp,  1.668205732358428E-01_wp, &
   & -2.739412762089972E-01_wp, -2.740564962399666E-01_wp,  1.399680745356982E-03_wp, &
   & -7.587990071245018E-05_wp, -1.161248023338480E-02_wp, -1.473786647645935E-02_wp, &
   &  5.766356950312447E-03_wp,  6.198485081399732E-03_wp, -3.638159143167712E-04_wp, &
   &  9.636383016019847E-03_wp, -4.320725957862492E-03_wp,  4.689580279485156E-03_wp, &
   &  3.825186283070631E-03_wp,  5.481470771605821E-04_wp, -1.315626420540229E-01_wp, &
   &  1.638287135548553E-01_wp, -2.731058070716672E-01_wp, -2.732219358542302E-01_wp, &
   & -1.395889570934521E-03_wp, -1.046208735740102E-04_wp, -1.162636921253166E-02_wp, &
   & -1.467070544573311E-02_wp, -5.818301210943058E-03_wp,  6.107986945215566E-03_wp, &
   & -4.721368135901439E-04_wp,  9.594368787516452E-03_wp,  4.299900967136565E-03_wp, &
   &  4.812990261177597E-03_wp,  3.932117466623117E-03_wp,  6.423804436812781E-04_wp, &
   &  6.313752291159144E-04_wp,  4.065079782192155E-02_wp,  3.762713946849568E-01_wp, &
   &  4.065079782192142E-02_wp,  2.694325003506606E-05_wp,  3.723791929504425E-03_wp, &
   & -1.426698129804384E-04_wp, -1.484219597730451E-03_wp, -1.616106489880607E-05_wp, &
   & -2.195886734633821E-03_wp,  1.725572476349207E-03_wp, -3.724237379731516E-03_wp, &
   & -1.883335093320207E-05_wp,  1.484721906252062E-03_wp,  4.332437613239408E-03_wp, &
   & -2.195319400570727E-03_wp,  1.320963754698336E-01_wp, -2.740564962399750E-01_wp, &
   &  1.668205732358601E-01_wp, -2.739412762090060E-01_wp,  4.299900967137431E-03_wp, &
   &  6.423804436743620E-04_wp,  4.812990261162023E-03_wp,  3.932117466614690E-03_wp, &
   & -1.395889570933940E-03_wp, -1.467070544572615E-02_wp, -1.046208735897015E-04_wp, &
   & -1.162636921254039E-02_wp,  5.766356950311810E-03_wp,  9.636383016012297E-03_wp, &
   &  6.198485081384104E-03_wp, -3.638159143073545E-04_wp, -1.315626420540227E-01_wp, &
   & -2.732219358542230E-01_wp,  1.638287135548399E-01_wp, -2.731058070716580E-01_wp, &
   & -4.320725957862151E-03_wp,  5.481470771673907E-04_wp,  4.689580279501103E-03_wp, &
   &  3.825186283078443E-03_wp,  1.399680745356025E-03_wp, -1.473786647646734E-02_wp, &
   & -7.587990069659981E-05_wp, -1.161248023337612E-02_wp, -5.818301210942681E-03_wp, &
   &  9.594368787522746E-03_wp,  6.107986945231997E-03_wp, -4.721368135991234E-04_wp, &
   &  6.313752291162483E-04_wp,  4.065079782192278E-02_wp,  4.065079782192142E-02_wp, &
   &  3.762713946849561E-01_wp, -1.883335093282772E-05_wp, -2.195319400570493E-03_wp, &
   &  1.484721906250675E-03_wp,  4.332437613238209E-03_wp,  2.694325003497794E-05_wp, &
   & -1.484219597729876E-03_wp,  3.723791929503831E-03_wp, -1.426698129813187E-04_wp, &
   & -1.616106489922502E-05_wp, -3.724237379732451E-03_wp, -2.195886734634651E-03_wp, &
   &  1.725572476349354E-03_wp,  1.320963754698339E-01_wp, -2.739412762089891E-01_wp, &
   & -2.740564962399838E-01_wp,  1.668205732358522E-01_wp, -5.818301210943532E-03_wp, &
   & -4.721368135830033E-04_wp,  9.594368787532023E-03_wp,  6.107986945223884E-03_wp, &
   & -4.320725957862965E-03_wp,  3.825186283063490E-03_wp,  5.481470771758441E-04_wp, &
   &  4.689580279494033E-03_wp, -1.395889570933077E-03_wp, -1.162636921252262E-02_wp, &
   & -1.467070544571731E-02_wp, -1.046208735832610E-04_wp, -1.315626420540235E-01_wp, &
   & -2.731058070716741E-01_wp, -2.732219358542142E-01_wp,  1.638287135548470E-01_wp, &
   &  5.766356950312010E-03_wp, -3.638159143236981E-04_wp,  9.636383016004037E-03_wp, &
   &  6.198485081391500E-03_wp,  4.299900967137845E-03_wp,  3.932117466631558E-03_wp, &
   &  6.423804436659595E-04_wp,  4.812990261169226E-03_wp,  1.399680745356946E-03_wp, &
   & -1.161248023339158E-02_wp, -1.473786647647505E-02_wp, -7.587990070342471E-05_wp, &
   & -1.002169991707603E-02_wp, -1.616106489918758E-05_wp,  2.694325003506606E-05_wp, &
   & -1.883335093282772E-05_wp,  5.926110700186242E-01_wp,  6.306062341021410E-04_wp, &
   & -6.296531093957627E-04_wp,  6.326277694937785E-04_wp, -1.002133098859540E-02_wp, &
   & -9.710340367719686E-06_wp,  3.903237184080666E-05_wp,  1.602074751050377E-05_wp, &
   & -1.002133098859570E-02_wp, -3.831356674489186E-05_wp, -1.424885720255911E-05_wp, &
   & -7.951709963053610E-06_wp, -8.079011306016456E-03_wp,  4.402759099971587E-03_wp, &
   & -1.013469187784766E-03_wp, -6.422687048651490E-03_wp,  2.675493330426835E-01_wp, &
   & -1.541473648586917E-01_wp,  1.541541864861006E-01_wp, -1.541018959659995E-01_wp, &
   & -8.034954641436455E-03_wp,  1.010497659766758E-03_wp,  6.420796560402954E-03_wp, &
   &  4.352317846913126E-03_wp, -8.036016332142348E-03_wp, -6.418480558179469E-03_wp, &
   & -4.350822085788245E-03_wp,  1.011694615660657E-03_wp, -8.065118430879478E-03_wp, &
   & -4.332027630450989E-03_wp,  1.067083504962701E-03_wp,  6.427992704465525E-03_wp, &
   &  2.667583325664288E-01_wp,  1.540383064982062E-01_wp, -1.540210211521305E-01_wp, &
   &  1.540730265314308E-01_wp, -8.109002085582512E-03_wp, -1.069687978375465E-03_wp, &
   & -6.429653425925694E-03_wp, -4.381969175444946E-03_wp, -8.107086159785016E-03_wp, &
   &  6.432110408948821E-03_wp,  4.383983009643901E-03_wp, -1.067922889876838E-03_wp, &
   & -1.652680392332771E-05_wp,  1.725572476348868E-03_wp,  3.723791929504425E-03_wp, &
   & -2.195319400570493E-03_wp,  6.306062341021410E-04_wp,  3.762791743573559E-01_wp, &
   & -4.064974927694869E-02_wp,  4.064981454521582E-02_wp, -7.951709963547733E-06_wp, &
   & -1.433240133743578E-04_wp,  1.484524434316980E-03_wp,  3.724452388100439E-03_wp, &
   &  3.903237184117638E-05_wp,  4.332380099949831E-03_wp,  2.195736794044228E-03_wp, &
   &  1.484524434316362E-03_wp,  1.396339312649658E-03_wp, -7.314752331259542E-05_wp, &
   &  1.160904696779823E-02_wp, -1.473751930385213E-02_wp,  1.320967967791354E-01_wp, &
   &  1.668162937630294E-01_wp,  2.739658646940500E-01_wp, -2.740603628281728E-01_wp, &
   &  4.297024085104940E-03_wp,  4.771615287614561E-03_wp, -3.898182431279804E-03_wp, &
   &  5.656963574224485E-04_wp, -5.821796072573855E-03_wp,  6.115206261215075E-03_wp, &
   &  3.791956525489121E-04_wp,  9.609460701613550E-03_wp, -1.398601382637568E-03_wp, &
   & -1.073486531672740E-04_wp,  1.162793322873198E-02_wp, -1.467051088107924E-02_wp, &
   & -1.315674126781334E-01_wp,  1.638070991959940E-01_wp,  2.730885160458183E-01_wp, &
   & -2.732252164136403E-01_wp, -4.321933415546995E-03_wp,  4.731900665515745E-03_wp, &
   & -3.856068523530349E-03_wp,  6.239611840168801E-04_wp,  5.761822710439690E-03_wp, &
   &  6.190050571836645E-03_wp,  4.561397155203686E-04_wp,  9.617562110800296E-03_wp, &
   & -2.796388824430549E-05_wp, -3.724237379732360E-03_wp, -1.426698129804384E-04_wp, &
   &  1.484721906250675E-03_wp, -6.296531093957627E-04_wp, -4.064974927694869E-02_wp, &
   &  3.762764645270570E-01_wp, -4.064843059989548E-02_wp, -3.831356674483379E-05_wp, &
   & -1.484800648992754E-03_wp,  4.332380099949475E-03_wp,  2.195840583569926E-03_wp, &
   &  1.602074751052610E-05_wp,  2.195840583569932E-03_wp,  1.725382374465028E-03_wp, &
   &  3.724452388099907E-03_wp, -4.292494327558612E-03_wp, -6.368124605281580E-04_wp, &
   &  4.816052244620153E-03_wp, -3.928962612510762E-03_wp, -1.320952064701817E-01_wp, &
   &  2.740838536132255E-01_wp,  1.668245955234447E-01_wp,  2.739411961597586E-01_wp, &
   &  5.822898098022933E-03_wp, -9.609452785748099E-03_wp,  6.116284732840767E-03_wp, &
   &  3.790844803209448E-04_wp, -1.395618826533400E-03_wp,  1.469071831614684E-02_wp, &
   & -9.349014039826982E-05_wp,  1.159285730895873E-02_wp,  4.326116442786803E-03_wp, &
   & -5.525493108820648E-04_wp,  4.687392010207347E-03_wp, -3.826572590914473E-03_wp, &
   &  1.315600071952982E-01_wp,  2.732019480397966E-01_wp,  1.638352610319876E-01_wp, &
   &  2.730949693464216E-01_wp, -5.759563306425127E-03_wp, -9.618332711403532E-03_wp, &
   &  6.190453672740990E-03_wp,  4.559846368499422E-04_wp,  1.400622608534982E-03_wp, &
   &  1.471706657317676E-02_wp, -8.709375104890387E-05_wp,  1.164521447807223E-02_wp, &
   & -2.111433219217073E-05_wp, -2.195886734634207E-03_wp, -1.484219597730451E-03_wp, &
   &  4.332437613238209E-03_wp,  6.326277694937785E-04_wp,  4.064981454521582E-02_wp, &
   & -4.064843059989548E-02_wp,  3.762765192929706E-01_wp, -1.424885720291436E-05_wp, &
   & -3.723983998564238E-03_wp,  2.195736794043728E-03_wp,  1.725382374464937E-03_wp, &
   & -9.710340367449303E-06_wp, -1.484800648991972E-03_wp, -3.723983998564540E-03_wp, &
   & -1.433240133749422E-04_wp, -5.811201835825922E-03_wp, -4.700665322186146E-04_wp, &
   & -9.592080771953320E-03_wp,  6.110681135735346E-03_wp,  1.320870341341682E-01_wp, &
   & -2.739393300672709E-01_wp,  2.740531769205482E-01_wp,  1.668759851262432E-01_wp, &
   &  1.398990802389267E-03_wp, -1.159011711667311E-02_wp,  1.469243854516845E-02_wp, &
   & -9.070800688119617E-05_wp,  4.298267508628735E-03_wp,  3.899519562816855E-03_wp, &
   & -5.699489028350755E-04_wp,  4.773621692205526E-03_wp,  5.771644567295561E-03_wp, &
   & -3.649895203934780E-04_wp, -9.636430266575064E-03_wp,  6.197209703178194E-03_wp, &
   & -1.315679157401549E-01_wp, -2.731148504997652E-01_wp,  2.732140653180960E-01_wp, &
   &  1.637840674194031E-01_wp, -1.394706983281260E-03_wp, -1.164734926221247E-02_wp, &
   &  1.471548099200957E-02_wp, -8.954771676624464E-05_wp, -4.322027981023274E-03_wp, &
   &  3.857120906489046E-03_wp, -6.205072471903762E-04_wp,  4.728892895504842E-03_wp, &
   & -1.002169991707614E-02_wp, -1.883335093316854E-05_wp, -1.616106489880607E-05_wp, &
   &  2.694325003497794E-05_wp, -1.002133098859540E-02_wp, -7.951709963547733E-06_wp, &
   & -3.831356674483379E-05_wp, -1.424885720291436E-05_wp,  5.926110700186238E-01_wp, &
   &  6.326277694937601E-04_wp,  6.306062341016666E-04_wp, -6.296531093957749E-04_wp, &
   & -1.002133098859575E-02_wp,  1.602074751069214E-05_wp, -9.710340366988948E-06_wp, &
   &  3.903237184092060E-05_wp, -8.079011306016029E-03_wp, -6.422687048651994E-03_wp, &
   &  4.402759099972186E-03_wp, -1.013469187783903E-03_wp, -8.107086159784718E-03_wp, &
   & -1.067922889876811E-03_wp,  6.432110408948320E-03_wp,  4.383983009644042E-03_wp, &
   &  2.667583325664288E-01_wp,  1.540730265314310E-01_wp,  1.540383064982057E-01_wp, &
   & -1.540210211521316E-01_wp, -8.034954641436892E-03_wp,  4.352317846912497E-03_wp, &
   &  1.010497659766600E-03_wp,  6.420796560402731E-03_wp, -8.065118430879696E-03_wp, &
   &  6.427992704465975E-03_wp, -4.332027630450825E-03_wp,  1.067083504962251E-03_wp, &
   & -8.036016332142520E-03_wp,  1.011694615660650E-03_wp, -6.418480558180584E-03_wp, &
   & -4.350822085788758E-03_wp,  2.675493330426832E-01_wp, -1.541018959659998E-01_wp, &
   & -1.541473648586913E-01_wp,  1.541541864861002E-01_wp, -8.109002085581485E-03_wp, &
   & -4.381969175445271E-03_wp, -1.069687978374858E-03_wp, -6.429653425926268E-03_wp, &
   & -2.111433219181250E-05_wp,  4.332437613238558E-03_wp, -2.195886734633821E-03_wp, &
   & -1.484219597729876E-03_wp, -9.710340367719686E-06_wp, -1.433240133743578E-04_wp, &
   & -1.484800648992754E-03_wp, -3.723983998564238E-03_wp,  6.326277694937601E-04_wp, &
   &  3.762765192929699E-01_wp,  4.064981454521582E-02_wp, -4.064843059989595E-02_wp, &
   & -1.424885720248168E-05_wp,  1.725382374465140E-03_wp, -3.723983998564704E-03_wp, &
   &  2.195736794043624E-03_wp, -5.811201835826739E-03_wp,  6.110681135737091E-03_wp, &
   & -4.700665322182961E-04_wp, -9.592080771956809E-03_wp, -4.322027981022819E-03_wp, &
   &  4.728892895507215E-03_wp,  3.857120906488954E-03_wp, -6.205072471940197E-04_wp, &
   & -1.315679157401556E-01_wp,  1.637840674194009E-01_wp, -2.731148504997642E-01_wp, &
   &  2.732140653180922E-01_wp,  1.398990802389900E-03_wp, -9.070800687766887E-05_wp, &
   & -1.159011711667271E-02_wp,  1.469243854517151E-02_wp,  5.771644567295471E-03_wp, &
   &  6.197209703175011E-03_wp, -3.649895203933488E-04_wp, -9.636430266573199E-03_wp, &
   &  4.298267508628347E-03_wp,  4.773621692202064E-03_wp,  3.899519562816786E-03_wp, &
   & -5.699489028334303E-04_wp,  1.320870341341678E-01_wp,  1.668759851262443E-01_wp, &
   & -2.739393300672721E-01_wp,  2.740531769205501E-01_wp, -1.394706983282141E-03_wp, &
   & -8.954771676848478E-05_wp, -1.164734926221190E-02_wp,  1.471548099200719E-02_wp, &
   & -1.652680392330776E-05_wp, -2.195319400570014E-03_wp,  1.725572476349207E-03_wp, &
   &  3.723791929503831E-03_wp,  3.903237184080666E-05_wp,  1.484524434316980E-03_wp, &
   &  4.332380099949475E-03_wp,  2.195736794043728E-03_wp,  6.306062341016666E-04_wp, &
   &  4.064981454521582E-02_wp,  3.762791743573567E-01_wp, -4.064974927694905E-02_wp, &
   & -7.951709962637613E-06_wp,  3.724452388100046E-03_wp, -1.433240133747481E-04_wp, &
   &  1.484524434316568E-03_wp,  1.396339312649446E-03_wp, -1.473751930385433E-02_wp, &
   & -7.314752333028112E-05_wp,  1.160904696781784E-02_wp,  5.761822710439412E-03_wp, &
   &  9.617562110797586E-03_wp,  6.190050571855734E-03_wp,  4.561397155401570E-04_wp, &
   & -1.315674126781348E-01_wp, -2.732252164136370E-01_wp,  1.638070991960129E-01_wp, &
   &  2.730885160458391E-01_wp,  4.297024085104610E-03_wp,  5.656963574184763E-04_wp, &
   &  4.771615287632452E-03_wp, -3.898182431300511E-03_wp, -1.398601382637399E-03_wp, &
   & -1.467051088107671E-02_wp, -1.073486531488078E-04_wp,  1.162793322871141E-02_wp, &
   & -5.821796072573463E-03_wp,  9.609460701616227E-03_wp,  6.115206261197411E-03_wp, &
   &  3.791956525273337E-04_wp,  1.320967967791354E-01_wp, -2.740603628281753E-01_wp, &
   &  1.668162937630118E-01_wp,  2.739658646940301E-01_wp, -4.321933415546939E-03_wp, &
   &  6.239611840200175E-04_wp,  4.731900665498968E-03_wp, -3.856068523510137E-03_wp, &
   & -2.796388824409472E-05_wp,  1.484721906251033E-03_wp, -3.724237379731516E-03_wp, &
   & -1.426698129813187E-04_wp,  1.602074751050377E-05_wp,  3.724452388100439E-03_wp, &
   &  2.195840583569926E-03_wp,  1.725382374464937E-03_wp, -6.296531093957749E-04_wp, &
   & -4.064843059989595E-02_wp, -4.064974927694905E-02_wp,  3.762764645270575E-01_wp, &
   & -3.831356674485109E-05_wp,  2.195840583570273E-03_wp, -1.484800648992161E-03_wp, &
   &  4.332380099949160E-03_wp, -4.292494327558464E-03_wp, -3.928962612509327E-03_wp, &
   & -6.368124605437831E-04_wp,  4.816052244635683E-03_wp,  1.400622608535111E-03_wp, &
   &  1.164521447807355E-02_wp,  1.471706657319379E-02_wp, -8.709375103392857E-05_wp, &
   &  1.315600071952971E-01_wp,  2.730949693464221E-01_wp,  2.732019480398148E-01_wp, &
   &  1.638352610320041E-01_wp,  5.822898098022540E-03_wp,  3.790844803206732E-04_wp, &
   & -9.609452785730773E-03_wp,  6.116284732824143E-03_wp,  4.326116442787243E-03_wp, &
   & -3.826572590915450E-03_wp, -5.525493108651967E-04_wp,  4.687392010190264E-03_wp, &
   & -1.395618826532265E-03_wp,  1.159285730895686E-02_wp,  1.469071831612966E-02_wp, &
   & -9.349014041595554E-05_wp, -1.320952064701819E-01_wp,  2.739411961597603E-01_wp, &
   &  2.740838536132076E-01_wp,  1.668245955234277E-01_wp, -5.759563306424172E-03_wp, &
   &  4.559846368483941E-04_wp, -9.618332711418983E-03_wp,  6.190453672757687E-03_wp, &
   & -1.002169991707632E-02_wp,  2.694325003484910E-05_wp, -1.883335093320207E-05_wp, &
   & -1.616106489922502E-05_wp, -1.002133098859570E-02_wp,  3.903237184117638E-05_wp, &
   &  1.602074751052610E-05_wp, -9.710340367449303E-06_wp, -1.002133098859575E-02_wp, &
   & -1.424885720248168E-05_wp, -7.951709962637613E-06_wp, -3.831356674485109E-05_wp, &
   &  5.926110700186241E-01_wp, -6.296531093953800E-04_wp,  6.326277694940046E-04_wp, &
   &  6.306062341014848E-04_wp, -8.079011306016704E-03_wp, -1.013469187784047E-03_wp, &
   & -6.422687048652284E-03_wp,  4.402759099972208E-03_wp, -8.109002085582346E-03_wp, &
   & -6.429653425925691E-03_wp, -4.381969175444765E-03_wp, -1.069687978375247E-03_wp, &
   & -8.036016332142737E-03_wp, -4.350822085788690E-03_wp,  1.011694615661138E-03_wp, &
   & -6.418480558179483E-03_wp,  2.667583325664297E-01_wp, -1.540210211521312E-01_wp, &
   &  1.540730265314304E-01_wp,  1.540383064982052E-01_wp, -8.065118430879795E-03_wp, &
   &  1.067083504962477E-03_wp,  6.427992704465863E-03_wp, -4.332027630450433E-03_wp, &
   & -8.034954641436611E-03_wp,  6.420796560403153E-03_wp,  4.352317846913310E-03_wp, &
   &  1.010497659766653E-03_wp, -8.107086159785748E-03_wp,  4.383983009644266E-03_wp, &
   & -1.067922889876910E-03_wp,  6.432110408949499E-03_wp,  2.675493330426822E-01_wp, &
   &  1.541541864861017E-01_wp, -1.541018959659994E-01_wp, -1.541473648586924E-01_wp, &
   & -2.796388824422173E-05_wp, -1.426698129811196E-04_wp,  1.484721906252062E-03_wp, &
   & -3.724237379732451E-03_wp, -3.831356674489186E-05_wp,  4.332380099949831E-03_wp, &
   &  2.195840583569932E-03_wp, -1.484800648991972E-03_wp,  1.602074751069214E-05_wp, &
   &  1.725382374465140E-03_wp,  3.724452388100046E-03_wp,  2.195840583570273E-03_wp, &
   & -6.296531093953800E-04_wp,  3.762764645270585E-01_wp, -4.064843059989492E-02_wp, &
   & -4.064974927694884E-02_wp, -4.292494327558534E-03_wp,  4.816052244616950E-03_wp, &
   & -3.928962612489824E-03_wp, -6.368124605451239E-04_wp, -5.759563306424233E-03_wp, &
   &  6.190453672740583E-03_wp,  4.559846368317824E-04_wp, -9.618332711420860E-03_wp, &
   & -1.395618826532905E-03_wp, -9.349014039769449E-05_wp,  1.159285730894220E-02_wp, &
   &  1.469071831613034E-02_wp,  1.315600071952972E-01_wp,  1.638352610319860E-01_wp, &
   &  2.730949693464056E-01_wp,  2.732019480398136E-01_wp,  4.326116442787101E-03_wp, &
   &  4.687392010209177E-03_wp, -3.826572590932637E-03_wp, -5.525493108660352E-04_wp, &
   &  5.822898098023783E-03_wp,  6.116284732840437E-03_wp,  3.790844803380568E-04_wp, &
   & -9.609452785732511E-03_wp,  1.400622608533285E-03_wp, -8.709375105211670E-05_wp, &
   &  1.164521447808818E-02_wp,  1.471706657319282E-02_wp, -1.320952064701815E-01_wp, &
   &  1.668245955234454E-01_wp,  2.739411961597784E-01_wp,  2.740838536132084E-01_wp, &
   & -2.111433219224735E-05_wp, -1.484219597730527E-03_wp,  4.332437613239408E-03_wp, &
   & -2.195886734634651E-03_wp, -1.424885720255911E-05_wp,  2.195736794044228E-03_wp, &
   &  1.725382374465028E-03_wp, -3.723983998564540E-03_wp, -9.710340366988948E-06_wp, &
   & -3.723983998564704E-03_wp, -1.433240133747481E-04_wp, -1.484800648992161E-03_wp, &
   &  6.326277694940046E-04_wp, -4.064843059989492E-02_wp,  3.762765192929720E-01_wp, &
   &  4.064981454521510E-02_wp, -5.811201835826484E-03_wp, -9.592080771958292E-03_wp, &
   &  6.110681135738235E-03_wp, -4.700665322178458E-04_wp, -1.394706983281508E-03_wp, &
   &  1.471548099200723E-02_wp, -8.954771676644244E-05_wp, -1.164734926221099E-02_wp, &
   &  4.298267508628526E-03_wp, -5.699489028320571E-04_wp,  4.773621692205643E-03_wp, &
   &  3.899519562819931E-03_wp, -1.315679157401564E-01_wp,  2.732140653180914E-01_wp, &
   &  1.637840674194038E-01_wp, -2.731148504997675E-01_wp,  5.771644567295396E-03_wp, &
   & -9.636430266571006E-03_wp,  6.197209703176661E-03_wp, -3.649895203967040E-04_wp, &
   &  1.398990802390154E-03_wp,  1.469243854517028E-02_wp, -9.070800688030699E-05_wp, &
   & -1.159011711667622E-02_wp, -4.322027981024076E-03_wp, -6.205072471946606E-04_wp, &
   &  4.728892895503906E-03_wp,  3.857120906485082E-03_wp,  1.320870341341683E-01_wp, &
   &  2.740531769205512E-01_wp,  1.668759851262438E-01_wp, -2.739393300672687E-01_wp, &
   & -1.652680392333367E-05_wp,  3.723791929503827E-03_wp, -2.195319400570727E-03_wp, &
   &  1.725572476349354E-03_wp, -7.951709963053610E-06_wp,  1.484524434316362E-03_wp, &
   &  3.724452388099907E-03_wp, -1.433240133749422E-04_wp,  3.903237184092060E-05_wp, &
   &  2.195736794043624E-03_wp,  1.484524434316568E-03_wp,  4.332380099949160E-03_wp, &
   &  6.306062341014848E-04_wp, -4.064974927694884E-02_wp,  4.064981454521510E-02_wp, &
   &  3.762791743573564E-01_wp,  1.396339312649570E-03_wp,  1.160904696779991E-02_wp, &
   & -1.473751930383382E-02_wp, -7.314752333317925E-05_wp, -4.321933415546646E-03_wp, &
   & -3.856068523527721E-03_wp,  6.239611839995426E-04_wp,  4.731900665495690E-03_wp, &
   & -5.821796072573427E-03_wp,  3.791956525460162E-04_wp,  9.609460701595565E-03_wp, &
   &  6.115206261195275E-03_wp, -1.315674126781342E-01_wp,  2.730885160458212E-01_wp, &
   & -2.732252164136577E-01_wp,  1.638070991960162E-01_wp, -1.398601382637849E-03_wp, &
   &  1.162793322873002E-02_wp, -1.467051088109777E-02_wp, -1.073486531463678E-04_wp, &
   &  4.297024085104993E-03_wp, -3.898182431282156E-03_wp,  5.656963574396877E-04_wp, &
   &  4.771615287635194E-03_wp,  5.761822710439592E-03_wp,  4.561397155227758E-04_wp, &
   &  9.617562110819180E-03_wp,  6.190050571858821E-03_wp,  1.320967967791345E-01_wp, &
   &  2.739658646940491E-01_wp, -2.740603628281547E-01_wp,  1.668162937630085E-01_wp, &
   &  2.675463048516357E-01_wp,  1.320963754698337E-01_wp,  1.320963754698336E-01_wp, &
   &  1.320963754698339E-01_wp, -8.079011306016456E-03_wp,  1.396339312649658E-03_wp, &
   & -4.292494327558612E-03_wp, -5.811201835825922E-03_wp, -8.079011306016029E-03_wp, &
   & -5.811201835826739E-03_wp,  1.396339312649446E-03_wp, -4.292494327558464E-03_wp, &
   & -8.079011306016704E-03_wp, -4.292494327558534E-03_wp, -5.811201835826484E-03_wp, &
   &  1.396339312649570E-03_wp,  1.306770741509492E+00_wp,  3.963681467213998E-01_wp, &
   &  3.963681467214000E-01_wp,  3.963681467213994E-01_wp, -4.029901788217475E-02_wp, &
   &  1.065946765565237E-04_wp,  3.642714200684565E-02_wp,  1.815444896577345E-02_wp, &
   & -2.255582463889032E-02_wp,  2.306123297462149E-02_wp, -1.724157098669024E-02_wp, &
   &  7.848904838977764E-03_wp, -2.255582463889024E-02_wp,  7.848904838978604E-03_wp, &
   &  2.306123297462161E-02_wp, -1.724157098669097E-02_wp, -2.161239237817555E-02_wp, &
   & -4.594919082002166E-02_wp, -4.594919082002197E-02_wp, -4.594919082002181E-02_wp, &
   & -2.255582463888970E-02_wp, -1.724157098669069E-02_wp,  7.848904838977870E-03_wp, &
   &  2.306123297462158E-02_wp, -4.029901788217441E-02_wp,  1.815444896577326E-02_wp, &
   &  1.065946765560373E-04_wp,  3.642714200684513E-02_wp, -4.029901788217424E-02_wp, &
   &  3.642714200684585E-02_wp,  1.815444896577343E-02_wp,  1.065946765559029E-04_wp, &
   & -1.541417593578825E-01_wp,  1.668205732358428E-01_wp, -2.740564962399750E-01_wp, &
   & -2.739412762089891E-01_wp,  4.402759099971587E-03_wp, -7.314752331259542E-05_wp, &
   & -6.368124605281580E-04_wp, -4.700665322186146E-04_wp, -6.422687048651994E-03_wp, &
   &  6.110681135737091E-03_wp, -1.473751930385433E-02_wp, -3.928962612509327E-03_wp, &
   & -1.013469187784047E-03_wp,  4.816052244616950E-03_wp, -9.592080771958292E-03_wp, &
   &  1.160904696779991E-02_wp,  3.963681467213998E-01_wp,  1.528839532883521E+00_wp, &
   & -7.960022803225419E-02_wp, -7.960022799236466E-02_wp,  1.075621700610065E-04_wp, &
   & -1.566472012310031E-02_wp,  3.968099543396018E-03_wp,  4.284465464028109E-03_wp, &
   & -2.311055743440250E-02_wp,  9.087948809263599E-03_wp, -6.384750134713867E-03_wp, &
   & -5.634928713373913E-03_wp,  7.829529797438364E-03_wp,  5.282500701841642E-02_wp, &
   &  5.498291996347014E-03_wp, -7.273193053826593E-03_wp,  4.637311361959880E-02_wp, &
   & -2.496902479358869E-01_wp,  2.009882403790397E-01_wp,  2.009852442730616E-01_wp, &
   &  1.720890802055751E-02_wp,  2.187942866920409E-03_wp,  7.354800366739529E-03_wp, &
   & -6.320786338396098E-03_wp,  1.814755668476659E-02_wp,  7.906486372857148E-02_wp, &
   &  4.285025385193222E-03_wp, -2.243131359845996E-03_wp, -3.642800403924312E-02_wp, &
   &  1.397885422753201E-02_wp,  2.245130836457984E-03_wp, -3.969989770460119E-03_wp, &
   & -1.541417593578830E-01_wp, -2.739412762089972E-01_wp,  1.668205732358601E-01_wp, &
   & -2.740564962399838E-01_wp, -1.013469187784766E-03_wp,  1.160904696779823E-02_wp, &
   &  4.816052244620153E-03_wp, -9.592080771953320E-03_wp,  4.402759099972186E-03_wp, &
   & -4.700665322182961E-04_wp, -7.314752333028112E-05_wp, -6.368124605437831E-04_wp, &
   & -6.422687048652284E-03_wp, -3.928962612489824E-03_wp,  6.110681135738235E-03_wp, &
   & -1.473751930383382E-02_wp,  3.963681467214000E-01_wp, -7.960022803225419E-02_wp, &
   &  1.528839532883653E+00_wp, -7.960022799249776E-02_wp, -3.642800403924346E-02_wp, &
   & -3.969989810322441E-03_wp,  1.397885422743508E-02_wp,  2.245130876227394E-03_wp, &
   &  1.720890802055159E-02_wp, -6.320786298267878E-03_wp,  2.187942866783865E-03_wp, &
   &  7.354800406725454E-03_wp, -2.311055743439637E-02_wp, -5.634928753465253E-03_wp, &
   &  9.087948809156803E-03_wp, -6.384750174708173E-03_wp,  4.637311361959871E-02_wp, &
   &  2.009852443131717E-01_wp, -2.496902479359525E-01_wp,  2.009882403389958E-01_wp, &
   &  7.829529797437762E-03_wp, -7.273193013697443E-03_wp,  5.282500701851719E-02_wp, &
   &  5.498291956321469E-03_wp,  1.075621700668561E-04_wp,  4.284465424162789E-03_wp, &
   & -1.566472012303417E-02_wp,  3.968099503592115E-03_wp,  1.814755668476133E-02_wp, &
   & -2.243131319949097E-03_wp,  7.906486372866436E-02_wp,  4.285025424987640E-03_wp, &
   & -1.541417593578828E-01_wp, -2.740564962399666E-01_wp, -2.739412762090060E-01_wp, &
   &  1.668205732358522E-01_wp, -6.422687048651490E-03_wp, -1.473751930385213E-02_wp, &
   & -3.928962612510762E-03_wp,  6.110681135735346E-03_wp, -1.013469187783903E-03_wp, &
   & -9.592080771956809E-03_wp,  1.160904696781784E-02_wp,  4.816052244635683E-03_wp, &
   &  4.402759099972208E-03_wp, -6.368124605451239E-04_wp, -4.700665322178458E-04_wp, &
   & -7.314752333317925E-05_wp,  3.963681467213994E-01_wp, -7.960022799236464E-02_wp, &
   & -7.960022799249775E-02_wp,  1.528839532843764E+00_wp,  1.814755668476572E-02_wp, &
   &  4.285025425055142E-03_wp, -2.243131359751072E-03_wp,  7.906486368880035E-02_wp, &
   &  7.829529797443335E-03_wp,  5.498291956216789E-03_wp, -7.273193053690188E-03_wp, &
   &  5.282500697842979E-02_wp,  1.720890802055125E-02_wp,  7.354800406829370E-03_wp, &
   & -6.320786338290391E-03_wp,  2.187942906913876E-03_wp,  4.637311361959900E-02_wp, &
   &  2.009882403389308E-01_wp,  2.009852442731270E-01_wp, -2.496902478958431E-01_wp, &
   & -2.311055743440203E-02_wp, -6.384750174843759E-03_wp, -5.634928713476720E-03_wp, &
   &  9.087948849287764E-03_wp, -3.642800403924789E-02_wp,  2.245130876323137E-03_wp, &
   & -3.969989770527052E-03_wp,  1.397885426733542E-02_wp,  1.075621700657132E-04_wp, &
   &  3.968099503497802E-03_wp,  4.284465463934503E-03_wp, -1.566472016289619E-02_wp, &
   & -8.075824133965537E-03_wp,  1.399680745356982E-03_wp,  4.299900967137431E-03_wp, &
   & -5.818301210943532E-03_wp,  2.675493330426835E-01_wp,  1.320967967791354E-01_wp, &
   & -1.320952064701817E-01_wp,  1.320870341341682E-01_wp, -8.107086159784718E-03_wp, &
   & -4.322027981022819E-03_wp,  5.761822710439412E-03_wp,  1.400622608535111E-03_wp, &
   & -8.109002085582346E-03_wp, -5.759563306424233E-03_wp, -1.394706983281508E-03_wp, &
   & -4.321933415546646E-03_wp, -4.029901788217475E-02_wp,  1.075621700610065E-04_wp, &
   & -3.642800403924346E-02_wp,  1.814755668476572E-02_wp,  1.306735148166990E+00_wp, &
   &  3.963998439273382E-01_wp, -3.964024246028324E-01_wp,  3.963273261897348E-01_wp, &
   & -4.028816596427909E-02_wp,  3.637482439148959E-02_wp, -1.806419269042958E-02_wp, &
   &  1.409580674174422E-04_wp, -4.028018006299157E-02_wp,  1.806138589744576E-02_wp, &
   & -1.458206527468416E-04_wp,  3.637650743496148E-02_wp, -2.255350359608576E-02_wp, &
   & -1.724206130738177E-02_wp, -7.848932803848118E-03_wp,  2.306398633753205E-02_wp, &
   & -2.161446192056909E-02_wp, -4.595454260063544E-02_wp,  4.596223619374454E-02_wp, &
   & -4.593460516438750E-02_wp, -2.252660675277589E-02_wp,  7.895767063900935E-03_wp, &
   & -2.313180601722654E-02_wp, -1.722157274331605E-02_wp, -2.252660675277604E-02_wp, &
   &  2.312902995055617E-02_wp,  1.721793885756783E-02_wp,  7.886773669638863E-03_wp, &
   &  4.405532199375228E-03_wp, -7.587990071245019E-05_wp,  6.423804436743620E-04_wp, &
   & -4.721368135830033E-04_wp, -1.541473648586917E-01_wp,  1.668162937630294E-01_wp, &
   &  2.740838536132255E-01_wp, -2.739393300672709E-01_wp, -1.067922889876811E-03_wp, &
   &  4.728892895507215E-03_wp,  9.617562110797586E-03_wp,  1.164521447807355E-02_wp, &
   & -6.429653425925691E-03_wp,  6.190453672740583E-03_wp,  1.471548099200723E-02_wp, &
   & -3.856068523527721E-03_wp,  1.065946765565237E-04_wp, -1.566472012310031E-02_wp, &
   & -3.969989810322448E-03_wp,  4.285025425055142E-03_wp,  3.963998439273382E-01_wp, &
   &  1.528796196898468E+00_wp,  7.959146314444271E-02_wp, -7.956353956197099E-02_wp, &
   & -3.638671220221724E-02_wp,  1.400670222150991E-02_wp, -2.207357045166773E-03_wp, &
   & -3.908689675460611E-03_wp,  1.812729918007535E-02_wp,  7.900309190998850E-02_wp, &
   & -4.198164018491705E-03_wp, -2.352970220666648E-03_wp,  1.720871391508201E-02_wp, &
   &  2.188039210862483E-03_wp, -7.357997821639504E-03_wp, -6.322190661943047E-03_wp, &
   &  4.637336383843259E-02_wp, -2.496907228624170E-01_wp, -2.009986631370533E-01_wp, &
   &  2.009757664580810E-01_wp,  7.886773669638350E-03_wp,  5.284121649500558E-02_wp, &
   & -5.618509734235139E-03_wp, -7.338105170410281E-03_wp, -2.313180601722668E-02_wp, &
   &  9.071651020712057E-03_wp,  6.321471759956468E-03_wp, -5.618509694401669E-03_wp, &
   &  1.010542468054945E-03_wp, -1.161248023338480E-02_wp,  4.812990261162023E-03_wp, &
   &  9.594368787532023E-03_wp,  1.541541864861006E-01_wp,  2.739658646940500E-01_wp, &
   &  1.668245955234447E-01_wp,  2.740531769205482E-01_wp,  6.432110408948320E-03_wp, &
   &  3.857120906488954E-03_wp,  6.190050571855734E-03_wp,  1.471706657319379E-02_wp, &
   & -4.381969175444765E-03_wp,  4.559846368317824E-04_wp, -8.954771676644244E-05_wp, &
   &  6.239611839995426E-04_wp,  3.642714200684565E-02_wp,  3.968099543396025E-03_wp, &
   &  1.397885422743508E-02_wp, -2.243131359751072E-03_wp, -3.964024246028324E-01_wp, &
   &  7.959146314444271E-02_wp,  1.528791104278062E+00_wp,  7.956794282497384E-02_wp, &
   & -1.813635682470457E-02_wp,  2.359242460137127E-03_wp,  7.900993954231281E-02_wp, &
   & -4.196085293246457E-03_wp, -1.329010635212524E-04_wp, -4.262042716463819E-03_wp, &
   & -1.565305775643328E-02_wp, -3.982699889496608E-03_wp, -7.830740345404170E-03_wp, &
   &  7.274386773369348E-03_wp,  5.282159533484578E-02_wp, -5.498856495633420E-03_wp, &
   & -4.636275770260743E-02_wp, -2.009911020334240E-01_wp, -2.497047114216899E-01_wp, &
   & -2.009642477465611E-01_wp,  2.312902995055177E-02_wp,  5.610420661574793E-03_wp, &
   &  9.071651020684517E-03_wp,  6.320222201022616E-03_wp, -1.722157274331205E-02_wp, &
   &  6.320222161172298E-03_wp,  2.170407961605759E-03_wp, -7.338105170382359E-03_wp, &
   & -6.422706566373671E-03_wp, -1.473786647645935E-02_wp,  3.932117466614690E-03_wp, &
   &  6.107986945223884E-03_wp, -1.541018959659995E-01_wp, -2.740603628281728E-01_wp, &
   &  2.739411961597586E-01_wp,  1.668759851262432E-01_wp,  4.383983009644042E-03_wp, &
   & -6.205072471940197E-04_wp,  4.561397155401570E-04_wp, -8.709375103392856E-05_wp, &
   & -1.069687978375247E-03_wp, -9.618332711420860E-03_wp, -1.164734926221099E-02_wp, &
   &  4.731900665495690E-03_wp,  1.815444896577345E-02_wp,  4.284465464028123E-03_wp, &
   &  2.245130876227394E-03_wp,  7.906486368880038E-02_wp,  3.963273261897348E-01_wp, &
   & -7.956353956197099E-02_wp,  7.956794282497384E-02_wp,  1.528836580818669E+00_wp, &
   &  1.288467804743353E-04_wp,  3.983391719479368E-03_wp, -4.261639573518425E-03_wp, &
   & -1.565680964673249E-02_wp, -3.638778655468135E-02_wp,  2.206137474764766E-03_wp, &
   &  3.908482794874726E-03_wp,  1.400580831898560E-02_wp, -2.311309351552541E-02_wp, &
   & -6.387442542612268E-03_wp,  5.637726485026015E-03_wp,  9.096835064801717E-03_wp, &
   &  4.639021136181751E-02_wp,  2.010139205395234E-01_wp, -2.010118728260463E-01_wp, &
   & -2.497045686259713E-01_wp,  1.721793885756373E-02_wp,  7.334109154952306E-03_wp, &
   &  6.321471799766568E-03_wp,  2.170408001480023E-03_wp,  7.895767063907195E-03_wp, &
   &  5.610420621700363E-03_wp,  7.334109154897767E-03_wp,  5.284121645519599E-02_wp, &
   & -8.067347093899502E-03_wp,  5.766356950312447E-03_wp, -1.395889570933940E-03_wp, &
   & -4.320725957862965E-03_wp, -8.034954641436455E-03_wp,  4.297024085104940E-03_wp, &
   &  5.822898098022933E-03_wp,  1.398990802389267E-03_wp,  2.667583325664288E-01_wp, &
   & -1.315679157401556E-01_wp, -1.315674126781348E-01_wp,  1.315600071952971E-01_wp, &
   & -8.036016332142737E-03_wp, -1.395618826532905E-03_wp,  4.298267508628526E-03_wp, &
   & -5.821796072573427E-03_wp, -2.255582463889032E-02_wp, -2.311055743440250E-02_wp, &
   &  1.720890802055159E-02_wp,  7.829529797443335E-03_wp, -4.028816596427909E-02_wp, &
   & -3.638671220221724E-02_wp, -1.813635682470457E-02_wp,  1.288467804743353E-04_wp, &
   &  1.309781562657337E+00_wp, -3.957889034317806E-01_wp, -3.957438745213444E-01_wp, &
   &  3.957138937754857E-01_wp, -2.258339700932242E-02_wp,  1.722989401278175E-02_wp, &
   & -7.786666353756151E-03_wp,  2.304357958047222E-02_wp, -4.026948104248459E-02_wp, &
   & -1.804462980363553E-02_wp, -1.672620921187963E-04_wp,  3.633494139675414E-02_wp, &
   & -2.258339700932225E-02_wp, -7.790642179568902E-03_wp, -2.304500118481929E-02_wp, &
   & -1.723099295689275E-02_wp, -2.161446192056871E-02_wp,  4.639021136181701E-02_wp, &
   &  4.637336383842627E-02_wp, -4.636275770261205E-02_wp, -4.028018006299187E-02_wp, &
   & -1.329010635151454E-04_wp, -3.638778655468153E-02_wp,  1.812729918008017E-02_wp, &
   &  6.427870946056665E-03_wp,  6.198485081399732E-03_wp, -1.467070544572615E-02_wp, &
   &  3.825186283063490E-03_wp,  1.010497659766758E-03_wp,  4.771615287614561E-03_wp, &
   & -9.609452785748099E-03_wp, -1.159011711667311E-02_wp,  1.540730265314310E-01_wp, &
   &  1.637840674194009E-01_wp, -2.732252164136370E-01_wp,  2.730949693464221E-01_wp, &
   & -4.350822085788690E-03_wp, -9.349014039769450E-05_wp, -5.699489028320571E-04_wp, &
   &  3.791956525460162E-04_wp,  2.306123297462149E-02_wp,  9.087948809263599E-03_wp, &
   & -6.320786298267878E-03_wp,  5.498291956216789E-03_wp,  3.637482439148959E-02_wp, &
   &  1.400670222150991E-02_wp,  2.359242460137127E-03_wp,  3.983391719479354E-03_wp, &
   & -3.957889034317806E-01_wp,  1.531462777264657E+00_wp, -8.098799896831285E-02_wp, &
   &  8.097836946446238E-02_wp, -1.723099295688689E-02_wp,  2.204196856083765E-03_wp, &
   &  7.294577257796801E-03_wp,  6.385321090147958E-03_wp, -1.804291761830162E-02_wp, &
   &  7.895023855106011E-02_wp,  4.174445770797119E-03_wp,  2.320905169525905E-03_wp, &
   & -7.786666353760748E-03_wp,  5.280602081260962E-02_wp,  5.521675939285944E-03_wp, &
   &  7.294577298167279E-03_wp, -4.593460516438819E-02_wp, -2.497045686660115E-01_wp, &
   &  2.009757664981821E-01_wp, -2.009642477464963E-01_wp, -1.458206527459482E-04_wp, &
   & -1.565305775649867E-02_wp,  3.908482834914621E-03_wp, -4.198164058591497E-03_wp, &
   & -4.329167530505221E-03_wp, -3.638159143167712E-04_wp, -1.046208735897015E-04_wp, &
   &  5.481470771758441E-04_wp,  6.420796560402954E-03_wp, -3.898182431279804E-03_wp, &
   &  6.116284732840767E-03_wp,  1.469243854516845E-02_wp,  1.540383064982057E-01_wp, &
   & -2.731148504997642E-01_wp,  1.638070991960129E-01_wp,  2.732019480398148E-01_wp, &
   &  1.011694615661138E-03_wp,  1.159285730894220E-02_wp,  4.773621692205643E-03_wp, &
   &  9.609460701595565E-03_wp, -1.724157098669024E-02_wp, -6.384750134713867E-03_wp, &
   &  2.187942866783865E-03_wp, -7.273193053690175E-03_wp, -1.806419269042958E-02_wp, &
   & -2.207357045166780E-03_wp,  7.900993954231281E-02_wp, -4.261639573518411E-03_wp, &
   & -3.957438745213444E-01_wp, -8.098799896831285E-02_wp,  1.531492082597841E+00_wp, &
   &  8.096005244473173E-02_wp, -7.790642179574978E-03_wp, -7.294004108261191E-03_wp, &
   &  5.280602081271573E-02_wp, -5.518768767943244E-03_wp, -1.677518521963276E-04_wp, &
   &  4.174210182478100E-03_wp, -1.564584338666537E-02_wp, -3.921119179288918E-03_wp, &
   &  2.304357958047743E-02_wp, -5.518768808176283E-03_wp,  9.113463583086644E-03_wp, &
   &  6.385321090041446E-03_wp, -4.595454260064046E-02_wp,  2.010139205796296E-01_wp, &
   & -2.496907228624863E-01_wp, -2.009911019933822E-01_wp,  3.637650743496574E-02_wp, &
   & -3.982699929539577E-03_wp,  1.400580827887927E-02_wp, -2.352970220597939E-03_wp, &
   & -1.068428092528660E-03_wp,  9.636383016019847E-03_wp, -1.162636921254039E-02_wp, &
   &  4.689580279494033E-03_wp,  4.352317846913126E-03_wp,  5.656963574224485E-04_wp, &
   &  3.790844803209448E-04_wp, -9.070800688119619E-05_wp, -1.540210211521316E-01_wp, &
   &  2.732140653180922E-01_wp,  2.730885160458391E-01_wp,  1.638352610320041E-01_wp, &
   & -6.418480558179483E-03_wp,  1.469071831613034E-02_wp,  3.899519562819931E-03_wp, &
   &  6.115206261195275E-03_wp,  7.848904838977764E-03_wp, -5.634928713373899E-03_wp, &
   &  7.354800406725454E-03_wp,  5.282500697842982E-02_wp,  1.409580674174422E-04_wp, &
   & -3.908689675460625E-03_wp, -4.196085293246443E-03_wp, -1.565680964673249E-02_wp, &
   &  3.957138937754857E-01_wp,  8.097836946446238E-02_wp,  8.096005244473173E-02_wp, &
   &  1.531507932996908E+00_wp, -2.304500118482008E-02_wp,  6.386786583804857E-03_wp, &
   &  5.521675899028591E-03_wp,  9.113463623417861E-03_wp, -3.633703846910721E-02_wp, &
   & -2.318659345945565E-03_wp,  3.922443497281367E-03_wp,  1.403428243334412E-02_wp, &
   &  1.722989401278325E-02_wp, -7.294004148590563E-03_wp,  6.386786583732998E-03_wp, &
   &  2.204196896340604E-03_wp,  4.596223619374088E-02_wp, -2.010118728259834E-01_wp, &
   & -2.009986630970265E-01_wp, -2.497047113815913E-01_wp,  1.806138589744916E-02_wp, &
   & -4.262042756561418E-03_wp,  2.206137474699568E-03_wp,  7.900309186996354E-02_wp, &
   & -8.067347093899963E-03_wp, -4.320725957862492E-03_wp,  5.766356950311810E-03_wp, &
   & -1.395889570933077E-03_wp, -8.036016332142348E-03_wp, -5.821796072573855E-03_wp, &
   & -1.395618826533400E-03_wp,  4.298267508628735E-03_wp, -8.034954641436892E-03_wp, &
   &  1.398990802389900E-03_wp,  4.297024085104610E-03_wp,  5.822898098022540E-03_wp, &
   &  2.667583325664297E-01_wp,  1.315600071952972E-01_wp, -1.315679157401564E-01_wp, &
   & -1.315674126781342E-01_wp, -2.255582463889024E-02_wp,  7.829529797438364E-03_wp, &
   & -2.311055743439637E-02_wp,  1.720890802055125E-02_wp, -4.028018006299157E-02_wp, &
   &  1.812729918007535E-02_wp, -1.329010635212524E-04_wp, -3.638778655468135E-02_wp, &
   & -2.258339700932242E-02_wp, -1.723099295688689E-02_wp, -7.790642179574978E-03_wp, &
   & -2.304500118482008E-02_wp,  1.309781562657334E+00_wp,  3.957138937754797E-01_wp, &
   & -3.957889034317847E-01_wp, -3.957438745213429E-01_wp, -4.026948104248455E-02_wp, &
   &  3.633494139675881E-02_wp, -1.804462980364132E-02_wp, -1.672620921185480E-04_wp, &
   & -2.258339700932276E-02_wp,  2.304357958047660E-02_wp,  1.722989401278798E-02_wp, &
   & -7.786666353754837E-03_wp, -4.028816596427713E-02_wp,  1.288467804703043E-04_wp, &
   & -3.638671220221046E-02_wp, -1.813635682470510E-02_wp, -2.161446192056909E-02_wp, &
   & -4.636275770260764E-02_wp,  4.639021136182270E-02_wp,  4.637336383842686E-02_wp, &
   & -1.068428092528640E-03_wp,  4.689580279485156E-03_wp,  9.636383016012297E-03_wp, &
   & -1.162636921252262E-02_wp, -6.418480558179469E-03_wp,  6.115206261215075E-03_wp, &
   &  1.469071831614684E-02_wp,  3.899519562816855E-03_wp,  4.352317846912497E-03_wp, &
   & -9.070800687766887E-05_wp,  5.656963574184763E-04_wp,  3.790844803206732E-04_wp, &
   & -1.540210211521312E-01_wp,  1.638352610319860E-01_wp,  2.732140653180914E-01_wp, &
   &  2.730885160458212E-01_wp,  7.848904838978604E-03_wp,  5.282500701841642E-02_wp, &
   & -5.634928753465253E-03_wp,  7.354800406829357E-03_wp,  1.806138589744576E-02_wp, &
   &  7.900309190998850E-02_wp, -4.262042716463826E-03_wp,  2.206137474764766E-03_wp, &
   &  1.722989401278175E-02_wp,  2.204196856083765E-03_wp, -7.294004108261191E-03_wp, &
   &  6.386786583804843E-03_wp,  3.957138937754797E-01_wp,  1.531507933037127E+00_wp, &
   &  8.097836950479281E-02_wp,  8.096005244462763E-02_wp, -3.633703846910621E-02_wp, &
   &  1.403428239306602E-02_wp, -2.318659305634754E-03_wp,  3.922443497246825E-03_wp, &
   & -2.304500118481468E-02_wp,  9.113463583192212E-03_wp,  6.386786543514121E-03_wp, &
   &  5.521675898956899E-03_wp,  1.409580674179586E-04_wp, -1.565680960672970E-02_wp, &
   & -3.908689715521114E-03_wp, -4.196085293310017E-03_wp,  4.596223619374492E-02_wp, &
   & -2.497047114216257E-01_wp, -2.010118728660495E-01_wp, -2.009986630969919E-01_wp, &
   &  6.427870946056199E-03_wp,  3.825186283070631E-03_wp,  6.198485081384104E-03_wp, &
   & -1.467070544571731E-02_wp, -4.350822085788245E-03_wp,  3.791956525489121E-04_wp, &
   & -9.349014039826982E-05_wp, -5.699489028350755E-04_wp,  1.010497659766600E-03_wp, &
   & -1.159011711667271E-02_wp,  4.771615287632452E-03_wp, -9.609452785730773E-03_wp, &
   &  1.540730265314304E-01_wp,  2.730949693464056E-01_wp,  1.637840674194038E-01_wp, &
   & -2.732252164136577E-01_wp,  2.306123297462161E-02_wp,  5.498291996347021E-03_wp, &
   &  9.087948809156809E-03_wp, -6.320786338290391E-03_wp, -1.458206527468416E-04_wp, &
   & -4.198164018491705E-03_wp, -1.565305775643328E-02_wp,  3.908482794874726E-03_wp, &
   & -7.786666353756151E-03_wp,  7.294577257796801E-03_wp,  5.280602081271572E-02_wp, &
   &  5.521675899028591E-03_wp, -3.957889034317847E-01_wp,  8.097836950479281E-02_wp, &
   &  1.531462777264731E+00_wp, -8.098799892804963E-02_wp, -1.804291761830099E-02_wp, &
   &  2.320905129175342E-03_wp,  7.895023855109590E-02_wp,  4.174445811111219E-03_wp, &
   & -1.723099295688687E-02_wp,  6.385321049778216E-03_wp,  2.204196856010261E-03_wp, &
   &  7.294577298092214E-03_wp,  3.637482439148552E-02_wp,  3.983391759584891E-03_wp, &
   &  1.400670222147266E-02_wp,  2.359242500209136E-03_wp, -4.593460516438379E-02_wp, &
   & -2.009642477866347E-01_wp, -2.497045686660787E-01_wp,  2.009757664581194E-01_wp, &
   & -4.329167530505186E-03_wp,  5.481470771605821E-04_wp, -3.638159143073545E-04_wp, &
   & -1.046208735832610E-04_wp,  1.011694615660657E-03_wp,  9.609460701613550E-03_wp, &
   &  1.159285730895873E-02_wp,  4.773621692205526E-03_wp,  6.420796560402731E-03_wp, &
   &  1.469243854517151E-02_wp, -3.898182431300511E-03_wp,  6.116284732824143E-03_wp, &
   &  1.540383064982052E-01_wp,  2.732019480398136E-01_wp, -2.731148504997675E-01_wp, &
   &  1.638070991960162E-01_wp, -1.724157098669097E-02_wp, -7.273193053826593E-03_wp, &
   & -6.384750174708173E-03_wp,  2.187942906913848E-03_wp,  3.637650743496148E-02_wp, &
   & -2.352970220666648E-03_wp, -3.982699889496621E-03_wp,  1.400580831898560E-02_wp, &
   &  2.304357958047222E-02_wp,  6.385321090147944E-03_wp, -5.518768767943230E-03_wp, &
   &  9.113463623417861E-03_wp, -3.957438745213429E-01_wp,  8.096005244462763E-02_wp, &
   & -8.098799892804963E-02_wp,  1.531492082557467E+00_wp, -1.677518521958641E-04_wp, &
   & -3.921119179222415E-03_wp,  4.174210222760599E-03_wp, -1.564584342701444E-02_wp, &
   & -7.790642179569497E-03_wp, -5.518768767805285E-03_wp, -7.294004148487201E-03_wp, &
   &  5.280602077234545E-02_wp, -1.806419269042532E-02_wp, -4.261639573614626E-03_wp, &
   & -2.207357085199105E-03_wp,  7.900993950217548E-02_wp, -4.595454260063910E-02_wp, &
   & -2.009911019932855E-01_wp,  2.010139205396199E-01_wp, -2.496907228223858E-01_wp, &
   &  2.667631089770241E-01_wp, -1.315626420540229E-01_wp, -1.315626420540227E-01_wp, &
   & -1.315626420540235E-01_wp, -8.065118430879478E-03_wp, -1.398601382637568E-03_wp, &
   &  4.326116442786803E-03_wp,  5.771644567295561E-03_wp, -8.065118430879696E-03_wp, &
   &  5.771644567295471E-03_wp, -1.398601382637399E-03_wp,  4.326116442787243E-03_wp, &
   & -8.065118430879795E-03_wp,  4.326116442787101E-03_wp,  5.771644567295396E-03_wp, &
   & -1.398601382637849E-03_wp, -2.161239237817555E-02_wp,  4.637311361959880E-02_wp, &
   &  4.637311361959871E-02_wp,  4.637311361959900E-02_wp, -2.255350359608576E-02_wp, &
   &  1.720871391508201E-02_wp, -7.830740345404170E-03_wp, -2.311309351552541E-02_wp, &
   & -4.026948104248459E-02_wp, -1.804291761830162E-02_wp, -1.677518521963276E-04_wp, &
   & -3.633703846910721E-02_wp, -4.026948104248455E-02_wp, -3.633703846910621E-02_wp, &
   & -1.804291761830099E-02_wp, -1.677518521958641E-04_wp,  1.309775452113892E+00_wp, &
   & -3.957521367355804E-01_wp, -3.957521367355803E-01_wp, -3.957521367355801E-01_wp, &
   & -4.026948104248463E-02_wp, -1.677518521964074E-04_wp, -3.633703846910725E-02_wp, &
   & -1.804291761830118E-02_wp, -2.255350359608572E-02_wp, -2.311309351552587E-02_wp, &
   &  1.720871391508191E-02_wp, -7.830740345404210E-03_wp, -2.255350359608517E-02_wp, &
   & -7.830740345404680E-03_wp, -2.311309351552493E-02_wp,  1.720871391508208E-02_wp, &
   &  1.540440853286230E-01_wp,  1.638287135548553E-01_wp, -2.732219358542230E-01_wp, &
   & -2.731058070716741E-01_wp, -4.332027630450989E-03_wp, -1.073486531672740E-04_wp, &
   & -5.525493108820648E-04_wp, -3.649895203934780E-04_wp,  6.427992704465975E-03_wp, &
   &  6.197209703175011E-03_wp, -1.467051088107671E-02_wp, -3.826572590915450E-03_wp, &
   &  1.067083504962477E-03_wp,  4.687392010209177E-03_wp, -9.636430266571006E-03_wp, &
   &  1.162793322873002E-02_wp, -4.594919082002166E-02_wp, -2.496902479358869E-01_wp, &
   &  2.009852443131717E-01_wp,  2.009882403389308E-01_wp, -1.724206130738177E-02_wp, &
   &  2.188039210862476E-03_wp,  7.274386773369355E-03_wp, -6.387442542612268E-03_wp, &
   & -1.804462980363553E-02_wp,  7.895023855106011E-02_wp,  4.174210182478100E-03_wp, &
   & -2.318659345945551E-03_wp,  3.633494139675881E-02_wp,  1.403428239306602E-02_wp, &
   &  2.320905129175349E-03_wp, -3.921119179222415E-03_wp, -3.957521367355804E-01_wp, &
   &  1.531465936200819E+00_wp, -8.096125854899061E-02_wp, -8.096125850866288E-02_wp, &
   & -1.672620921236748E-04_wp, -1.564584338673149E-02_wp,  3.922443537557163E-03_wp, &
   &  4.174445811146080E-03_wp,  2.306398633753208E-02_wp,  9.096835024743219E-03_wp, &
   & -6.322190621863628E-03_wp, -5.498856495607149E-03_wp, -7.848932803848064E-03_wp, &
   &  5.282159533481898E-02_wp,  5.637726525084957E-03_wp, -7.357997861720061E-03_wp, &
   &  1.540440853286227E-01_wp, -2.731058070716672E-01_wp,  1.638287135548399E-01_wp, &
   & -2.732219358542142E-01_wp,  1.067083504962701E-03_wp,  1.162793322873198E-02_wp, &
   &  4.687392010207347E-03_wp, -9.636430266575064E-03_wp, -4.332027630450825E-03_wp, &
   & -3.649895203933488E-04_wp, -1.073486531488078E-04_wp, -5.525493108651967E-04_wp, &
   &  6.427992704465863E-03_wp, -3.826572590932637E-03_wp,  6.197209703176661E-03_wp, &
   & -1.467051088109777E-02_wp, -4.594919082002197E-02_wp,  2.009882403790397E-01_wp, &
   & -2.496902479359525E-01_wp,  2.009852442731270E-01_wp, -7.848932803848118E-03_wp, &
   & -7.357997821639504E-03_wp,  5.282159533484578E-02_wp,  5.637726485026001E-03_wp, &
   & -1.672620921187963E-04_wp,  4.174445770797119E-03_wp, -1.564584338666537E-02_wp, &
   &  3.922443497281367E-03_wp, -1.804462980364132E-02_wp, -2.318659305634747E-03_wp, &
   &  7.895023855109590E-02_wp,  4.174210222760613E-03_wp, -3.957521367355803E-01_wp, &
   & -8.096125854899063E-02_wp,  1.531465936200814E+00_wp, -8.096125850865754E-02_wp, &
   &  3.633494139675919E-02_wp, -3.921119219572250E-03_wp,  1.403428239303365E-02_wp, &
   &  2.320905169490267E-03_wp, -1.724206130738499E-02_wp, -6.387442502527312E-03_wp, &
   &  2.188039210866057E-03_wp,  7.274386813460626E-03_wp,  2.306398633753603E-02_wp, &
   & -5.498856535724551E-03_wp,  9.096835024717372E-03_wp, -6.322190661946420E-03_wp, &
   &  1.540440853286232E-01_wp, -2.732219358542302E-01_wp, -2.731058070716580E-01_wp, &
   &  1.638287135548470E-01_wp,  6.427992704465525E-03_wp, -1.467051088107924E-02_wp, &
   & -3.826572590914473E-03_wp,  6.197209703178194E-03_wp,  1.067083504962251E-03_wp, &
   & -9.636430266573199E-03_wp,  1.162793322871141E-02_wp,  4.687392010190264E-03_wp, &
   & -4.332027630450433E-03_wp, -5.525493108660352E-04_wp, -3.649895203967040E-04_wp, &
   & -1.073486531463678E-04_wp, -4.594919082002181E-02_wp,  2.009852442730616E-01_wp, &
   &  2.009882403389958E-01_wp, -2.496902478958431E-01_wp,  2.306398633753205E-02_wp, &
   & -6.322190661943033E-03_wp, -5.498856495633420E-03_wp,  9.096835064801745E-03_wp, &
   &  3.633494139675414E-02_wp,  2.320905169525905E-03_wp, -3.921119179288904E-03_wp, &
   &  1.403428243334412E-02_wp, -1.672620921185480E-04_wp,  3.922443497246825E-03_wp, &
   &  4.174445811111205E-03_wp, -1.564584342701444E-02_wp, -3.957521367355801E-01_wp, &
   & -8.096125850866288E-02_wp, -8.096125850865754E-02_wp,  1.531465936160484E+00_wp, &
   & -1.804462980363579E-02_wp,  4.174210222827296E-03_wp, -2.318659345911564E-03_wp, &
   &  7.895023851074685E-02_wp, -7.848932803843849E-03_wp,  5.637726485000633E-03_wp, &
   & -7.357997861721893E-03_wp,  5.282159529472830E-02_wp, -1.724206130738620E-02_wp, &
   &  7.274386813486883E-03_wp, -6.387442542586635E-03_wp,  2.188039250945628E-03_wp, &
   & -8.067347093899700E-03_wp, -1.395889570934521E-03_wp, -4.320725957862151E-03_wp, &
   &  5.766356950312010E-03_wp,  2.667583325664288E-01_wp, -1.315674126781334E-01_wp, &
   &  1.315600071952982E-01_wp, -1.315679157401549E-01_wp, -8.036016332142520E-03_wp, &
   &  4.298267508628347E-03_wp, -5.821796072573463E-03_wp, -1.395618826532265E-03_wp, &
   & -8.034954641436611E-03_wp,  5.822898098023783E-03_wp,  1.398990802390154E-03_wp, &
   &  4.297024085104993E-03_wp, -2.255582463888970E-02_wp,  1.720890802055751E-02_wp, &
   &  7.829529797437762E-03_wp, -2.311055743440203E-02_wp, -2.161446192056909E-02_wp, &
   &  4.637336383843259E-02_wp, -4.636275770260743E-02_wp,  4.639021136181751E-02_wp, &
   & -2.258339700932225E-02_wp, -7.786666353760748E-03_wp,  2.304357958047743E-02_wp, &
   &  1.722989401278325E-02_wp, -2.258339700932276E-02_wp, -2.304500118481468E-02_wp, &
   & -1.723099295688687E-02_wp, -7.790642179569497E-03_wp, -4.026948104248463E-02_wp, &
   & -1.672620921236748E-04_wp,  3.633494139675919E-02_wp, -1.804462980363579E-02_wp, &
   &  1.309781562657336E+00_wp, -3.957438745213488E-01_wp,  3.957138937754792E-01_wp, &
   & -3.957889034317807E-01_wp, -4.028018006299169E-02_wp, -3.638778655467666E-02_wp, &
   &  1.812729918007519E-02_wp, -1.329010635157967E-04_wp, -4.028816596427742E-02_wp, &
   & -1.813635682471038E-02_wp,  1.288467804690999E-04_wp, -3.638671220221605E-02_wp, &
   & -4.329167530504967E-03_wp, -1.046208735740102E-04_wp,  5.481470771673907E-04_wp, &
   & -3.638159143236981E-04_wp,  1.540383064982062E-01_wp,  1.638070991959940E-01_wp, &
   &  2.732019480397966E-01_wp, -2.731148504997652E-01_wp,  1.011694615660650E-03_wp, &
   &  4.773621692202064E-03_wp,  9.609460701616227E-03_wp,  1.159285730895686E-02_wp, &
   &  6.420796560403153E-03_wp,  6.116284732840437E-03_wp,  1.469243854517028E-02_wp, &
   & -3.898182431282156E-03_wp, -1.724157098669069E-02_wp,  2.187942866920416E-03_wp, &
   & -7.273193013697436E-03_wp, -6.384750174843759E-03_wp, -4.595454260063544E-02_wp, &
   & -2.496907228624170E-01_wp, -2.009911020334240E-01_wp,  2.010139205395234E-01_wp, &
   & -7.790642179568902E-03_wp,  5.280602081260962E-02_wp, -5.518768808176290E-03_wp, &
   & -7.294004148590563E-03_wp,  2.304357958047660E-02_wp,  9.113463583192212E-03_wp, &
   &  6.385321049778216E-03_wp, -5.518768767805285E-03_wp, -1.677518521964074E-04_wp, &
   & -1.564584338673149E-02_wp, -3.921119219572250E-03_wp,  4.174210222827296E-03_wp, &
   & -3.957438745213488E-01_wp,  1.531492082597702E+00_wp,  8.096005248495813E-02_wp, &
   & -8.098799892794255E-02_wp,  3.637650743496145E-02_wp,  1.400580827897723E-02_wp, &
   & -2.352970180565267E-03_wp, -3.982699889400448E-03_wp, -1.806419269042981E-02_wp, &
   &  7.900993954221643E-02_wp, -4.261639533508860E-03_wp, -2.207357085267078E-03_wp, &
   & -1.068428092529037E-03_wp, -1.162636921253166E-02_wp,  4.689580279501103E-03_wp, &
   &  9.636383016004037E-03_wp, -1.540210211521305E-01_wp,  2.730885160458183E-01_wp, &
   &  1.638352610319876E-01_wp,  2.732140653180960E-01_wp, -6.418480558180584E-03_wp, &
   &  3.899519562816786E-03_wp,  6.115206261197411E-03_wp,  1.469071831612966E-02_wp, &
   &  4.352317846913310E-03_wp,  3.790844803380568E-04_wp, -9.070800688030699E-05_wp, &
   &  5.656963574396877E-04_wp,  7.848904838977870E-03_wp,  7.354800366739529E-03_wp, &
   &  5.282500701851719E-02_wp, -5.634928713476733E-03_wp,  4.596223619374454E-02_wp, &
   & -2.009986631370533E-01_wp, -2.497047114216899E-01_wp, -2.010118728260463E-01_wp, &
   & -2.304500118481929E-02_wp,  5.521675939285937E-03_wp,  9.113463583086644E-03_wp, &
   &  6.386786583732984E-03_wp,  1.722989401278798E-02_wp,  6.386786543514121E-03_wp, &
   &  2.204196856010254E-03_wp, -7.294004148487201E-03_wp, -3.633703846910725E-02_wp, &
   &  3.922443537557163E-03_wp,  1.403428239303364E-02_wp, -2.318659345911564E-03_wp, &
   &  3.957138937754792E-01_wp,  8.096005248495812E-02_wp,  1.531507933037199E+00_wp, &
   &  8.097836946453590E-02_wp,  1.806138589745004E-02_wp,  2.206137434698857E-03_wp, &
   &  7.900309191002428E-02_wp, -4.262042756497733E-03_wp,  1.409580674128876E-04_wp, &
   & -4.196085253212599E-03_wp, -1.565680960666847E-02_wp, -3.908689675496249E-03_wp, &
   &  6.427870946056195E-03_wp, -1.467070544573311E-02_wp,  3.825186283078443E-03_wp, &
   &  6.198485081391500E-03_wp,  1.540730265314308E-01_wp, -2.732252164136403E-01_wp, &
   &  2.730949693464216E-01_wp,  1.637840674194031E-01_wp, -4.350822085788758E-03_wp, &
   & -5.699489028334303E-04_wp,  3.791956525273337E-04_wp, -9.349014041595554E-05_wp, &
   &  1.010497659766653E-03_wp, -9.609452785732511E-03_wp, -1.159011711667622E-02_wp, &
   &  4.771615287635194E-03_wp,  2.306123297462158E-02_wp, -6.320786338396098E-03_wp, &
   &  5.498291956321469E-03_wp,  9.087948849287764E-03_wp, -4.593460516438750E-02_wp, &
   &  2.009757664580810E-01_wp, -2.009642477465611E-01_wp, -2.497045686259713E-01_wp, &
   & -1.723099295689275E-02_wp,  7.294577298167279E-03_wp,  6.385321090041446E-03_wp, &
   &  2.204196896340577E-03_wp, -7.786666353754837E-03_wp,  5.521675898956899E-03_wp, &
   &  7.294577298092200E-03_wp,  5.280602077234545E-02_wp, -1.804291761830118E-02_wp, &
   &  4.174445811146080E-03_wp,  2.320905169490253E-03_wp,  7.895023851074685E-02_wp, &
   & -3.957889034317807E-01_wp, -8.098799892794255E-02_wp,  8.097836946453590E-02_wp, &
   &  1.531462777224363E+00_wp, -1.458206527413296E-04_wp,  3.908482794810486E-03_wp, &
   & -4.198164058553874E-03_wp, -1.565305779657139E-02_wp,  3.637482439148515E-02_wp, &
   &  2.359242500275166E-03_wp,  3.983391719543428E-03_wp,  1.400670226157377E-02_wp, &
   & -8.075824133965659E-03_wp, -5.818301210943058E-03_wp,  1.399680745356025E-03_wp, &
   &  4.299900967137845E-03_wp, -8.109002085582512E-03_wp, -4.321933415546995E-03_wp, &
   & -5.759563306425127E-03_wp, -1.394706983281260E-03_wp,  2.675493330426832E-01_wp, &
   &  1.320870341341678E-01_wp,  1.320967967791354E-01_wp, -1.320952064701819E-01_wp, &
   & -8.107086159785748E-03_wp,  1.400622608533285E-03_wp, -4.322027981024076E-03_wp, &
   &  5.761822710439592E-03_wp, -4.029901788217441E-02_wp,  1.814755668476659E-02_wp, &
   &  1.075621700668561E-04_wp, -3.642800403924789E-02_wp, -2.252660675277589E-02_wp, &
   &  7.886773669638350E-03_wp,  2.312902995055177E-02_wp,  1.721793885756373E-02_wp, &
   & -2.161446192056871E-02_wp, -4.593460516438819E-02_wp, -4.595454260064046E-02_wp, &
   &  4.596223619374088E-02_wp, -4.028816596427713E-02_wp,  1.409580674179586E-04_wp, &
   &  3.637482439148552E-02_wp, -1.806419269042532E-02_wp, -2.255350359608572E-02_wp, &
   &  2.306398633753208E-02_wp, -1.724206130738499E-02_wp, -7.848932803843849E-03_wp, &
   & -4.028018006299169E-02_wp,  3.637650743496145E-02_wp,  1.806138589745004E-02_wp, &
   & -1.458206527413296E-04_wp,  1.306735148166990E+00_wp,  3.963273261897356E-01_wp, &
   &  3.963998439273415E-01_wp, -3.964024246028278E-01_wp, -2.252660675277623E-02_wp, &
   & -1.722157274331661E-02_wp,  7.895767063907448E-03_wp, -2.313180601723036E-02_wp, &
   & -6.422706566372709E-03_wp,  6.107986945215566E-03_wp, -1.473786647646734E-02_wp, &
   &  3.932117466631558E-03_wp, -1.069687978375465E-03_wp,  4.731900665515745E-03_wp, &
   & -9.618332711403532E-03_wp, -1.164734926221247E-02_wp, -1.541018959659998E-01_wp, &
   &  1.668759851262443E-01_wp, -2.740603628281753E-01_wp,  2.739411961597603E-01_wp, &
   &  4.383983009644266E-03_wp, -8.709375105211670E-05_wp, -6.205072471946606E-04_wp, &
   &  4.561397155227758E-04_wp,  1.815444896577326E-02_wp,  7.906486372857147E-02_wp, &
   &  4.284465424162789E-03_wp,  2.245130876323137E-03_wp,  7.895767063900935E-03_wp, &
   &  5.284121649500558E-02_wp,  5.610420661574786E-03_wp,  7.334109154952320E-03_wp, &
   &  4.639021136181701E-02_wp, -2.497045686660115E-01_wp,  2.010139205796296E-01_wp, &
   & -2.010118728259834E-01_wp,  1.288467804703043E-04_wp, -1.565680960672970E-02_wp, &
   &  3.983391759584891E-03_wp, -4.261639573614626E-03_wp, -2.311309351552587E-02_wp, &
   &  9.096835024743219E-03_wp, -6.387442502527312E-03_wp,  5.637726485000619E-03_wp, &
   & -3.638778655467666E-02_wp,  1.400580827897723E-02_wp,  2.206137434698850E-03_wp, &
   &  3.908482794810472E-03_wp,  3.963273261897356E-01_wp,  1.528836580858457E+00_wp, &
   & -7.956353960180840E-02_wp,  7.956794282491923E-02_wp,  1.721793885756801E-02_wp, &
   &  2.170407961661583E-03_wp,  7.334109115111315E-03_wp,  6.321471799794559E-03_wp, &
   &  4.405532199374985E-03_wp, -4.721368135901439E-04_wp, -7.587990069659981E-05_wp, &
   &  6.423804436659595E-04_wp, -6.429653425925694E-03_wp, -3.856068523530349E-03_wp, &
   &  6.190453672740990E-03_wp,  1.471548099200957E-02_wp, -1.541473648586913E-01_wp, &
   & -2.739393300672721E-01_wp,  1.668162937630118E-01_wp,  2.740838536132076E-01_wp, &
   & -1.067922889876910E-03_wp,  1.164521447808818E-02_wp,  4.728892895503906E-03_wp, &
   &  9.617562110819180E-03_wp,  1.065946765560373E-04_wp,  4.285025385193215E-03_wp, &
   & -1.566472012303418E-02_wp, -3.969989770527052E-03_wp, -2.313180601722654E-02_wp, &
   & -5.618509734235139E-03_wp,  9.071651020684517E-03_wp,  6.321471799766581E-03_wp, &
   &  4.637336383842627E-02_wp,  2.009757664981821E-01_wp, -2.496907228624863E-01_wp, &
   & -2.009986630970265E-01_wp, -3.638671220221046E-02_wp, -3.908689715521121E-03_wp, &
   &  1.400670222147266E-02_wp, -2.207357085199119E-03_wp,  1.720871391508191E-02_wp, &
   & -6.322190621863621E-03_wp,  2.188039210866057E-03_wp, -7.357997861721907E-03_wp, &
   &  1.812729918007519E-02_wp, -2.352970180565260E-03_wp,  7.900309191002428E-02_wp, &
   & -4.198164058553888E-03_wp,  3.963998439273415E-01_wp, -7.956353960180841E-02_wp, &
   &  1.528796196898468E+00_wp,  7.959146310460047E-02_wp,  7.886773669634445E-03_wp, &
   & -7.338105130539736E-03_wp,  5.284121649503365E-02_wp, -5.618509694401475E-03_wp, &
   &  1.010542468055552E-03_wp,  9.594368787516452E-03_wp, -1.161248023337612E-02_wp, &
   &  4.812990261169226E-03_wp, -4.381969175444946E-03_wp,  6.239611840168801E-04_wp, &
   &  4.559846368499422E-04_wp, -8.954771676624464E-05_wp,  1.541541864861002E-01_wp, &
   &  2.740531769205501E-01_wp,  2.739658646940301E-01_wp,  1.668245955234277E-01_wp, &
   &  6.432110408949499E-03_wp,  1.471706657319282E-02_wp,  3.857120906485082E-03_wp, &
   &  6.190050571858821E-03_wp,  3.642714200684513E-02_wp, -2.243131359846010E-03_wp, &
   &  3.968099503592101E-03_wp,  1.397885426733542E-02_wp, -1.722157274331605E-02_wp, &
   & -7.338105170410281E-03_wp,  6.320222201022616E-03_wp,  2.170408001480051E-03_wp, &
   & -4.636275770261205E-02_wp, -2.009642477464962E-01_wp, -2.009911019933821E-01_wp, &
   & -2.497047113815913E-01_wp, -1.813635682470510E-02_wp, -4.196085293310031E-03_wp, &
   &  2.359242500209136E-03_wp,  7.900993950217546E-02_wp, -7.830740345404210E-03_wp, &
   & -5.498856495607135E-03_wp,  7.274386813460626E-03_wp,  5.282159529472827E-02_wp, &
   & -1.329010635157967E-04_wp, -3.982699889400448E-03_wp, -4.262042756497747E-03_wp, &
   & -1.565305779657139E-02_wp, -3.964024246028278E-01_wp,  7.956794282491922E-02_wp, &
   &  7.959146310460047E-02_wp,  1.528791104238155E+00_wp,  2.312902995055146E-02_wp, &
   &  6.320222201079459E-03_wp,  5.610420621755291E-03_wp,  9.071651060554770E-03_wp, &
   & -8.075824133965872E-03_wp,  4.299900967136565E-03_wp, -5.818301210942681E-03_wp, &
   &  1.399680745356946E-03_wp, -8.107086159785016E-03_wp,  5.761822710439690E-03_wp, &
   &  1.400622608534982E-03_wp, -4.322027981023274E-03_wp, -8.109002085581485E-03_wp, &
   & -1.394706983282141E-03_wp, -4.321933415546939E-03_wp, -5.759563306424172E-03_wp, &
   &  2.675493330426822E-01_wp, -1.320952064701815E-01_wp,  1.320870341341683E-01_wp, &
   &  1.320967967791345E-01_wp, -4.029901788217424E-02_wp, -3.642800403924312E-02_wp, &
   &  1.814755668476133E-02_wp,  1.075621700657132E-04_wp, -2.252660675277604E-02_wp, &
   & -2.313180601722668E-02_wp, -1.722157274331205E-02_wp,  7.895767063907195E-03_wp, &
   & -4.028018006299187E-02_wp, -1.458206527459482E-04_wp,  3.637650743496574E-02_wp, &
   &  1.806138589744916E-02_wp, -2.161446192056909E-02_wp,  4.596223619374492E-02_wp, &
   & -4.593460516438379E-02_wp, -4.595454260063910E-02_wp, -2.255350359608517E-02_wp, &
   & -7.848932803848064E-03_wp,  2.306398633753603E-02_wp, -1.724206130738620E-02_wp, &
   & -4.028816596427742E-02_wp, -1.806419269042981E-02_wp,  1.409580674128876E-04_wp, &
   &  3.637482439148515E-02_wp, -2.252660675277623E-02_wp,  1.721793885756801E-02_wp, &
   &  7.886773669634445E-03_wp,  2.312902995055146E-02_wp,  1.306735148166993E+00_wp, &
   & -3.964024246028335E-01_wp,  3.963273261897298E-01_wp,  3.963998439273430E-01_wp, &
   &  1.010542468055235E-03_wp,  4.812990261177597E-03_wp,  9.594368787522746E-03_wp, &
   & -1.161248023339158E-02_wp,  6.432110408948821E-03_wp,  6.190050571836645E-03_wp, &
   &  1.471706657317676E-02_wp,  3.857120906489046E-03_wp, -4.381969175445271E-03_wp, &
   & -8.954771676848478E-05_wp,  6.239611840200175E-04_wp,  4.559846368483941E-04_wp, &
   &  1.541541864861017E-01_wp,  1.668245955234454E-01_wp,  2.740531769205512E-01_wp, &
   &  2.739658646940491E-01_wp,  3.642714200684585E-02_wp,  1.397885422753201E-02_wp, &
   & -2.243131319949097E-03_wp,  3.968099503497816E-03_wp,  2.312902995055617E-02_wp, &
   &  9.071651020712050E-03_wp,  6.320222161172305E-03_wp,  5.610420621700363E-03_wp, &
   & -1.329010635151454E-04_wp, -1.565305775649867E-02_wp, -3.982699929539577E-03_wp, &
   & -4.262042756561404E-03_wp, -4.636275770260764E-02_wp, -2.497047114216257E-01_wp, &
   & -2.009642477866347E-01_wp, -2.009911019932855E-01_wp, -7.830740345404680E-03_wp, &
   &  5.282159533481898E-02_wp, -5.498856535724544E-03_wp,  7.274386813486869E-03_wp, &
   & -1.813635682471038E-02_wp,  7.900993954221643E-02_wp, -4.196085253212599E-03_wp, &
   &  2.359242500275166E-03_wp, -1.722157274331661E-02_wp,  2.170407961661576E-03_wp, &
   & -7.338105130539743E-03_wp,  6.320222201079445E-03_wp, -3.964024246028335E-01_wp, &
   &  1.528791104278004E+00_wp,  7.956794286479173E-02_wp,  7.959146310457306E-02_wp, &
   & -6.422706566373787E-03_wp,  3.932117466623117E-03_wp,  6.107986945231997E-03_wp, &
   & -1.473786647647505E-02_wp,  4.383983009643901E-03_wp,  4.561397155203686E-04_wp, &
   & -8.709375104890387E-05_wp, -6.205072471903762E-04_wp, -1.069687978374858E-03_wp, &
   & -1.164734926221190E-02_wp,  4.731900665498968E-03_wp, -9.618332711418983E-03_wp, &
   & -1.541018959659994E-01_wp,  2.739411961597784E-01_wp,  1.668759851262438E-01_wp, &
   & -2.740603628281547E-01_wp,  1.815444896577343E-02_wp,  2.245130836457984E-03_wp, &
   &  7.906486372866436E-02_wp,  4.284465463934503E-03_wp,  1.721793885756783E-02_wp, &
   &  6.321471759956475E-03_wp,  2.170407961605753E-03_wp,  7.334109154897767E-03_wp, &
   & -3.638778655468153E-02_wp,  3.908482834914621E-03_wp,  1.400580827887927E-02_wp, &
   &  2.206137474699568E-03_wp,  4.639021136182270E-02_wp, -2.010118728660495E-01_wp, &
   & -2.497045686660787E-01_wp,  2.010139205396199E-01_wp, -2.311309351552493E-02_wp, &
   &  5.637726525084957E-03_wp,  9.096835024717365E-03_wp, -6.387442542586635E-03_wp, &
   &  1.288467804690999E-04_wp, -4.261639533508867E-03_wp, -1.565680960666847E-02_wp, &
   &  3.983391719543428E-03_wp,  7.895767063907448E-03_wp,  7.334109115111315E-03_wp, &
   &  5.284121649503366E-02_wp,  5.610420621755305E-03_wp,  3.963273261897298E-01_wp, &
   &  7.956794286479173E-02_wp,  1.528836580858510E+00_wp, -7.956353956199855E-02_wp, &
   &  4.405532199374922E-03_wp,  6.423804436812781E-04_wp, -4.721368135991234E-04_wp, &
   & -7.587990070342471E-05_wp, -1.067922889876838E-03_wp,  9.617562110800296E-03_wp, &
   &  1.164521447807223E-02_wp,  4.728892895504842E-03_wp, -6.429653425926268E-03_wp, &
   &  1.471548099200719E-02_wp, -3.856068523510137E-03_wp,  6.190453672757687E-03_wp, &
   & -1.541473648586924E-01_wp,  2.740838536132084E-01_wp, -2.739393300672687E-01_wp, &
   &  1.668162937630085E-01_wp,  1.065946765559029E-04_wp, -3.969989770460106E-03_wp, &
   &  4.285025424987654E-03_wp, -1.566472016289619E-02_wp,  7.886773669638863E-03_wp, &
   & -5.618509694401669E-03_wp, -7.338105170382359E-03_wp,  5.284121645519596E-02_wp, &
   &  1.812729918008017E-02_wp, -4.198164058591497E-03_wp, -2.352970220597952E-03_wp, &
   &  7.900309186996354E-02_wp,  4.637336383842686E-02_wp, -2.009986630969918E-01_wp, &
   &  2.009757664581194E-01_wp, -2.496907228223858E-01_wp,  1.720871391508208E-02_wp, &
   & -7.357997861720075E-03_wp, -6.322190661946434E-03_wp,  2.188039250945628E-03_wp, &
   & -3.638671220221605E-02_wp, -2.207357085267078E-03_wp, -3.908689675496235E-03_wp, &
   &  1.400670226157380E-02_wp, -2.313180601723036E-02_wp,  6.321471799794559E-03_wp, &
   & -5.618509694401475E-03_wp,  9.071651060554770E-03_wp,  3.963998439273430E-01_wp, &
   &  7.959146310457305E-02_wp, -7.956353956199853E-02_wp,  1.528796196858635E+00_wp],&
   &  shape(density))

   call get_structure(mol, "X23", "CO2")

   call test_num_op_grad(error, mol, density, make_exchange_gxtb, thr_in=thr2)

end subroutine test_op_g_fock_co2


subroutine test_bvk_exchange_kernel(error)
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: num(2) = [1, 1]
   integer, parameter :: num_sc(6) = [1, 1, 1, 1, 1, 1]
   integer, parameter :: nmesh(3) = [3, 1, 1]
   integer, parameter :: nmesh_large(3) = [9, 9, 1]
   real(wp), parameter :: xyz(3, 2) = reshape([ &
      & 0.40_wp, 0.30_wp, 0.20_wp, &
      & 1.35_wp, 0.35_wp, 0.25_wp], [3, 2])/autoaa
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      & 3.20_wp, 0.00_wp, 0.00_wp, &
      & 0.15_wp, 20.0_wp, 0.00_wp, &
      & 0.05_wp, 0.10_wp, 21.0_wp], [3, 3])/autoaa
   logical, parameter :: periodic(3) = [.true., .true., .true.]

   type(structure_type) :: mol, mol_probe, mol_sc
   class(qvszp_basis_type), allocatable :: qbas, qbas_sc
   class(basis_type), allocatable :: bas, bas_sc
   class(exchange_type), allocatable :: exchange, exchange_sc
   type(basis_cache) :: bcache, bcache_sc
   type(container_cache) :: cache, cache_cp2k, cache_sc
   type(exchange_cache), pointer :: ptr, ptr_cp2k, ptr_sc
   type(exchange_bvk_kernel) :: kernel, kernel_large, kernel_one
   type(cp2k_exchange_stream_type) :: stream, stream_large, stream_order, &
      & stream_permuted, stream_reverse, stream_reverse_missing
   type(xtb_calculator) :: calc_cp2k
   type(wavefunction_type) :: wfn, wfn_aux, wfn_sc, wfn_aux_sc
   integer, parameter :: korder(3) = [3, 1, 2]
   integer :: acell, bcell, bvk_builds, i, iat, icell, ik, ir, ish, j, jat, jsh
   integer :: first, ii, izp, jao, jj, jzp, last, nao, ncell, ni, nj, nsh
   real(wp) :: angle, deriv_fd, deriv_ref, energy_large, &
      & energy_large_oracle, energy_mesh, energy_m, energy_p, energy_order, &
      & energy_stream, energy_stream_order, energy_twist, onsite_probe, order_err
   real(wp) :: energy_sc, err, kfrac(3, 3), kfrac_bad(3, 3), kfrac_order(3, 3), &
      & kfrac_twist(3, 3), lattice_sc(3, 3), step
   real(wp) :: weights(3), weights_bad(3), weights_order(3), xyz_sc(3, 6), shift(3)
   real(wp) :: coord_direction(3, 2), identity(3, 3), lattice0(3, 3), &
      & strain_direction(3, 3), transform(3, 3), xyz0(3, 2)
   complex(wp) :: phase
   real(wp), allocatable :: density_sc(:, :, :), overlap_sc(:, :)
   real(wp), allocatable :: ao_grad_sc(:, :), bocorr_grad_r(:, :, :), &
      & bocorr_grad_large(:, :, :), bocorr_grad_sc(:, :), &
      & gradient_cp2k(:, :), gradient_large(:, :), gradient_large_oracle(:, :), &
      & gradient_mesh(:, :), gradient_sc(:, :), gradient_stream_probe(:, :), &
      & gradient_stream_reverse(:, :), &
      & kernel_bocorr_backup(:, :, :), kernel_bocorr_direction(:, :, :), &
      & kernel_mulliken_backup(:, :, :), &
      & kernel_mulliken_direction(:, :, :), kfrac_large(:, :), &
      & mulliken_grad_large(:, :, :), mulliken_grad_r(:, :, :), &
      & mulliken_grad_sc(:, :), p_r(:, :, :), s_r(:, :, :), &
      & sigma_cp2k(:, :), sigma_large(:, :), sigma_large_oracle(:, :), &
      & sigma_mesh(:, :), sigma_sc(:, :), sigma_stream_probe(:, :), &
      & sigma_stream_reverse(:, :), vsh_cp2k(:), vsh_large(:), &
      & vsh_large_oracle(:), vsh_mesh(:), vsh_order(:), vsh_ref(:), &
      & vsh_stream(:), vsh_stream_order(:), vsh_twist(:), weights_large(:)
   complex(wp), allocatable :: density_mesh(:, :, :, :), &
      & density_m(:, :, :, :), density_p(:, :, :, :), &
      & density_large(:, :, :, :), density_order(:, :, :, :), &
      & direction(:, :, :, :), fock_cp2k(:, :, :, :), &
      & fock_large(:, :, :, :), fock_large_oracle(:, :, :, :), &
      & fock_mesh(:, :, :, :), &
      & fock_stream(:, :, :, :), fock_stream_order(:, :, :, :), &
      & fock_order(:, :, :, :), &
      & fock_m(:, :, :, :), fock_p(:, :, :, :), fock_ref(:, :, :, :), &
      & fock_twist(:, :, :, :), overlap_direction(:, :, :), &
      & overlap_grad_cp2k(:, :, :), overlap_grad_large(:, :, :), &
      & overlap_grad_large_oracle(:, :, :), overlap_grad_mesh(:, :, :), &
      & overlap_grad_stream_probe(:, :), overlap_grad_stream_reverse(:, :, :), &
      & overlap_large(:, :, :), overlap_m(:, :, :), overlap_mesh(:, :, :), &
      & overlap_order(:, :, :), overlap_p(:, :, :)

   call new(mol, num, xyz, lattice=lattice, periodic=periodic)
   allocate(qbas)
   call make_qvszp_basis(qbas, mol, error)
   if (allocated(error)) return
   call move_alloc(qbas, bas)
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
   call make_exchange_gxtb(exchange, mol, bas)
   call exchange%update(mol, cache)
   call view(cache, ptr)

   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix(mol, [1, 1, 1], kernel_one)
      call fock%get_bvk_Kmatrix(mol, nmesh, kernel)
   class default
      call test_failed(error, "Unexpected exchange implementation")
      return
   end select

   err = max(maxval(abs(kernel_one%g_mulliken_r(:, :, 1) &
      & -ptr%g_mulliken)), maxval(abs(kernel_one%g_bocorr_r(:, :, 1) &
      & -ptr%g_bocorr)))

   lattice_sc = lattice
   lattice_sc(:, 1) = real(nmesh(1), wp)*lattice(:, 1)
   do icell = 1, nmesh(1)
      shift = real(icell-1, wp)*lattice(:, 1)
      xyz_sc(:, 2*icell-1:2*icell) = xyz + spread(shift, 2, 2)
   end do
   call new(mol_sc, num_sc, xyz_sc, lattice=lattice_sc, periodic=periodic)
   allocate(qbas_sc)
   call make_qvszp_basis(qbas_sc, mol_sc, error)
   if (allocated(error)) return
   call move_alloc(qbas_sc, bas_sc)
   call new_wavefunction(wfn_aux_sc, mol_sc%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol_sc, error, wfn_aux_sc%qat(:, 1))
   if (allocated(error)) return
   call bas_sc%update(mol_sc, bcache_sc, .false., wfn_aux=wfn_aux_sc)
   call make_exchange_gxtb(exchange_sc, mol_sc, bas_sc)
   call exchange_sc%update(mol_sc, cache_sc)
   call view(cache_sc, ptr_sc)

   nsh = bas%nsh
   if (bas_sc%nsh /= nmesh(1)*nsh) then
      call test_failed(error, "Unexpected replicated basis layout")
      return
   end if
   do icell = 1, product(nmesh)
      first = (icell-1)*nsh + 1
      last = icell*nsh
      err = max(err, maxval(abs(kernel%g_mulliken_r(:, :, icell) &
         & -ptr_sc%g_mulliken(1:nsh, first:last))), &
         & maxval(abs(kernel%g_bocorr_r(:, :, icell) &
         & -ptr_sc%g_bocorr(1:2, 2*icell-1:2*icell))))
   end do
   if (any(kernel%reps /= reshape([0, 0, 0, 1, 0, 0, 2, 0, 0], &
      & shape(kernel%reps)))) err = max(err, 1.0_wp)

   call check(error, err, 0.0_wp, thr=1.0e-12_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Image-resolved exchange kernels do not reproduce the BvK supercell")
      return
   end if

   ! Construct deliberately non-symmetric real image blocks.  R=2 is the
   ! transpose of R=1, so their Fourier transforms are genuinely complex
   ! Hermitian matrices and the assembled supercell is real symmetric.
   nao = bas%nao
   ncell = product(nmesh)
   if (bas_sc%nao /= ncell*nao) then
      call test_failed(error, "Unexpected replicated AO layout")
      return
   end if
   allocate(p_r(nao, nao, ncell), s_r(nao, nao, ncell), source=0.0_wp)
   do i = 1, nao
      s_r(i, i, 1) = 1.0_wp + 0.04_wp*real(i, wp)
      p_r(i, i, 1) = 0.55_wp + 0.07_wp*real(i, wp)
      do j = 1, i-1
         s_r(j, i, 1) = 0.03_wp/real(i+j, wp)
         s_r(i, j, 1) = s_r(j, i, 1)
         p_r(j, i, 1) = 0.04_wp*cos(real(2*i+3*j, wp))
         p_r(i, j, 1) = p_r(j, i, 1)
      end do
      do j = 1, nao
         s_r(i, j, 2) = 0.012_wp*sin(real(3*i+5*j, wp))
         p_r(i, j, 2) = 0.018_wp*cos(real(7*i+2*j, wp))
      end do
   end do
   s_r(:, :, 3) = transpose(s_r(:, :, 2))
   p_r(:, :, 3) = transpose(p_r(:, :, 2))

   allocate(overlap_sc(ncell*nao, ncell*nao), &
      & density_sc(ncell*nao, ncell*nao, 1), source=0.0_wp)
   do acell = 0, ncell-1
      do bcell = 0, ncell-1
         ir = modulo(bcell-acell, ncell) + 1
         overlap_sc(acell*nao+1:(acell+1)*nao, &
            & bcell*nao+1:(bcell+1)*nao) = s_r(:, :, ir)
         density_sc(acell*nao+1:(acell+1)*nao, &
            & bcell*nao+1:(bcell+1)*nao, 1) = p_r(:, :, ir)
      end do
   end do

   kfrac = 0.0_wp
   weights = 1.0_wp/real(ncell, wp)
   do ik = 1, ncell
      kfrac(1, ik) = real(ik-1, wp)/real(ncell, wp)
   end do
   allocate(overlap_mesh(nao, nao, ncell), &
      & density_mesh(nao, nao, 1, ncell), &
      & fock_mesh(nao, nao, 1, ncell), &
      & fock_ref(nao, nao, 1, ncell), source=(0.0_wp, 0.0_wp))
   allocate(vsh_mesh(nsh), vsh_ref(nsh), source=0.0_wp)
   do ik = 1, ncell
      do ir = 1, ncell
         angle = 2.0_wp*pi*kfrac(1, ik)*real(ir-1, wp)
         phase = exp(cmplx(0.0_wp, angle, wp))
         overlap_mesh(:, :, ik) = overlap_mesh(:, :, ik) &
            & + phase*s_r(:, :, ir)
         density_mesh(:, :, 1, ik) = density_mesh(:, :, 1, ik) &
            & + phase*p_r(:, :, ir)
      end do
   end do

   call new_wavefunction(wfn, mol%nat, nsh, nao, 1, 0.0_wp, .false.)
   do i = 1, nsh
      wfn%qsh(i, 1) = 0.08_wp*sin(real(3*i, wp))
   end do
   call exchange%get_onsite_Kmatrix(mol, wfn, ptr)

   call new_wavefunction(wfn_sc, mol_sc%nat, bas_sc%nsh, bas_sc%nao, &
      & 1, 0.0_wp, .false.)
   do icell = 1, ncell
      first = (icell-1)*nsh + 1
      last = icell*nsh
      wfn_sc%qsh(first:last, 1) = wfn%qsh(:, 1)
   end do
   call exchange_sc%get_onsite_Kmatrix(mol_sc, wfn_sc, ptr_sc)
   allocate(ptr_sc%prev_F(bas_sc%nao, bas_sc%nao, 1), &
      & ptr_sc%prev_vsh(bas_sc%nsh, 1))
   call exchange_sc%get_KFock(mol_sc, ptr_sc, density_sc, overlap_sc)

   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_mesh, vsh_mesh, energy_mesh)
   class default
      call test_failed(error, "Unexpected exchange implementation")
      return
   end select

   do ik = 1, ncell
      do ir = 1, ncell
         first = (ir-1)*nao + 1
         last = ir*nao
         angle = 2.0_wp*pi*kfrac(1, ik)*real(ir-1, wp)
         phase = exp(cmplx(0.0_wp, angle, wp))
         fock_ref(:, :, 1, ik) = fock_ref(:, :, 1, ik) &
            & + phase*ptr_sc%prev_F(1:nao, first:last, 1)
      end do
   end do
   energy_sc = 0.5_wp*sum(density_sc*ptr_sc%prev_F)/real(ncell, wp)
   vsh_ref = 0.0_wp
   do icell = 1, ncell
      first = (icell-1)*nsh + 1
      last = icell*nsh
      vsh_ref = vsh_ref + ptr_sc%prev_vsh(first:last, 1)/real(ncell, wp)
   end do

   err = max(maxval(abs(fock_mesh-fock_ref)), abs(energy_mesh-energy_sc), &
      & maxval(abs(vsh_mesh-vsh_ref)))
   do icell = 1, ncell
      first = (icell-1)*nsh + 1
      last = icell*nsh
      err = max(err, maxval(abs(vsh_mesh-ptr_sc%prev_vsh(first:last, 1))))
   end do

   ! A perturbation at one non-Gamma block must change the other blocks via
   ! the image kernels, while the BZ energy derivative remains the Fock map.
   allocate(direction(nao, nao, 1, ncell), &
      & density_m(nao, nao, 1, ncell), density_p(nao, nao, 1, ncell), &
      & fock_m(nao, nao, 1, ncell), fock_p(nao, nao, 1, ncell), &
      & source=(0.0_wp, 0.0_wp))
   do i = 1, nao
      direction(i, i, 1, 2) = cmplx(0.01_wp*sin(real(5*i, wp)), 0.0_wp, wp)
      do j = 1, i-1
         direction(j, i, 1, 2) = cmplx( &
            & 0.01_wp*cos(real(3*i+j, wp)), &
            & 0.01_wp*sin(real(i+4*j, wp)), wp)
         direction(i, j, 1, 2) = conjg(direction(j, i, 1, 2))
      end do
   end do
   step = 1.0e-5_wp
   density_m = density_mesh-step*direction
   density_p = density_mesh+step*direction
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac, weights, &
         & density_m, overlap_mesh, fock_m, vsh_mesh, energy_m)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac, weights, &
         & density_p, overlap_mesh, fock_p, vsh_mesh, energy_p)
   end select
   deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
   deriv_ref = 0.0_wp
   do ik = 1, ncell
      deriv_ref = deriv_ref + weights(ik)*real(sum( &
         & conjg(direction(:, :, :, ik))*fock_mesh(:, :, :, ik)), wp)
   end do
   err = max(err, abs(deriv_fd-deriv_ref))
   if (maxval(abs(fock_p(:, :, :, 1)-fock_m(:, :, :, 1))) &
      & < 1.0e-12_wp) err = max(err, 1.0_wp)

   ! Reverse the complete-mesh functional with respect to independently
   ! Hermitian S(k) blocks.  The public response is unweighted, so the BZ
   ! contraction below contains exactly one factor of w_k.
   allocate(overlap_grad_mesh(nao, nao, ncell), &
      & overlap_direction(nao, nao, ncell), overlap_m(nao, nao, ncell), &
      & overlap_p(nao, nao, ncell), source=(0.0_wp, 0.0_wp))
   allocate(mulliken_grad_r(nao, nao, ncell), &
      & bocorr_grad_r(nao, nao, ncell), source=0.0_wp)
   do ik = 1, ncell
      do i = 1, nao
         overlap_direction(i, i, ik) = cmplx( &
            & 0.01_wp*sin(real(5*i+7*ik, wp)), 0.0_wp, wp)
         do j = 1, i-1
            overlap_direction(j, i, ik) = cmplx( &
               & 0.01_wp*cos(real(3*i+j+2*ik, wp)), &
               & 0.01_wp*sin(real(i+4*j+3*ik, wp)), wp)
            overlap_direction(i, j, ik) = &
               & conjg(overlap_direction(j, i, ik))
         end do
      end do
   end do
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_KGrad_kmesh(mol, ptr, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, overlap_grad_mesh, &
         & mulliken_grad_r, bocorr_grad_r)
      overlap_m = overlap_mesh-step*overlap_direction
      overlap_p = overlap_mesh+step*overlap_direction
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac, weights, &
         & density_mesh, overlap_m, fock_m, vsh_mesh, energy_m)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac, weights, &
         & density_mesh, overlap_p, fock_p, vsh_mesh, energy_p)
   end select
   deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
   deriv_ref = 0.0_wp
   do ik = 1, ncell
      deriv_ref = deriv_ref + weights(ik)*real(sum( &
         & conjg(overlap_grad_mesh(:, :, ik)) &
         & *overlap_direction(:, :, ik)), wp)
      err = max(err, maxval(abs(overlap_grad_mesh(:, :, ik) &
         & -conjg(transpose(overlap_grad_mesh(:, :, ik))))))
   end do
   err = max(err, abs(deriv_fd-deriv_ref))

   ! Perturb every oriented image kernel independently.  In particular no
   ! R/-R or row/column symmetrization is allowed in this contraction.
   allocate(kernel_mulliken_backup(nsh, nsh, ncell), &
      & kernel_mulliken_direction(nsh, nsh, ncell), &
      & kernel_bocorr_backup(mol%nat, mol%nat, ncell), &
      & kernel_bocorr_direction(mol%nat, mol%nat, ncell))
   kernel_mulliken_backup = kernel%g_mulliken_r
   kernel_bocorr_backup = kernel%g_bocorr_r
   do icell = 1, ncell
      do ish = 1, nsh
         do jsh = 1, nsh
            kernel_mulliken_direction(ish, jsh, icell) = &
               & 0.01_wp*sin(real(11*ish+7*jsh+5*icell, wp))
         end do
      end do
      do iat = 1, mol%nat
         do jat = 1, mol%nat
            kernel_bocorr_direction(iat, jat, icell) = &
               & 0.01_wp*cos(real(3*iat+13*jat+2*icell, wp))
         end do
      end do
   end do
   kernel%g_mulliken_r = kernel_mulliken_backup &
      & - step*kernel_mulliken_direction
   kernel%g_bocorr_r = kernel_bocorr_backup-step*kernel_bocorr_direction
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_m, vsh_mesh, energy_m)
   end select
   kernel%g_mulliken_r = kernel_mulliken_backup &
      & + step*kernel_mulliken_direction
   kernel%g_bocorr_r = kernel_bocorr_backup+step*kernel_bocorr_direction
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_p, vsh_mesh, energy_p)
   end select
   kernel%g_mulliken_r = kernel_mulliken_backup
   kernel%g_bocorr_r = kernel_bocorr_backup
   deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
   deriv_ref = 0.0_wp
   do icell = 1, ncell
      do ish = 1, nsh
         ii = bas%iao_sh(ish)
         ni = bas%nao_sh(ish)
         do jsh = 1, nsh
            jj = bas%iao_sh(jsh)
            nj = bas%nao_sh(jsh)
            deriv_ref = deriv_ref &
               & + kernel_mulliken_direction(ish, jsh, icell) &
               & *sum(mulliken_grad_r(ii+1:ii+ni, jj+1:jj+nj, icell))
         end do
      end do
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, mol%nat
            jzp = mol%id(jat)
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(bas%ish_at(iat)+ish)
               ni = bas%nao_sh(bas%ish_at(iat)+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(bas%ish_at(jat)+jsh)
                  nj = bas%nao_sh(bas%ish_at(jat)+jsh)
                  deriv_ref = deriv_ref &
                     & + kernel_bocorr_direction(iat, jat, icell) &
                     & *sum(bocorr_grad_r(ii+1:ii+ni, jj+1:jj+nj, icell))
               end do
            end do
         end do
      end do
   end do
   err = max(err, abs(deriv_fd-deriv_ref))

   ! Contract the oriented image responses with the very same BvK
   ! Wigner-Seitz literals used in the forward kernel construction.
   allocate(gradient_mesh(3, mol%nat), sigma_mesh(3, 3), source=0.0_wp)
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix_derivs(mol, kernel, mulliken_grad_r, &
         & bocorr_grad_r, gradient_mesh, sigma_mesh)
   end select
   do iat = 1, mol%nat
      do i = 1, 3
         coord_direction(i, iat) = 0.03_wp*sin(real(7*i+11*iat, wp))
      end do
   end do
   xyz0 = mol%xyz
   mol%xyz = xyz0-step*coord_direction
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix(mol, nmesh, kernel)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_m, vsh_mesh, energy_m)
   end select
   mol%xyz = xyz0+step*coord_direction
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix(mol, nmesh, kernel)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_p, vsh_mesh, energy_p)
   end select
   deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
   deriv_ref = sum(gradient_mesh*coord_direction)
   err = max(err, abs(deriv_fd-deriv_ref))

   identity = 0.0_wp
   do i = 1, 3
      identity(i, i) = 1.0_wp
   end do
   strain_direction = reshape([ &
      & 0.07_wp, 0.02_wp, -0.01_wp, &
      & 0.02_wp, -0.04_wp, 0.03_wp, &
      & -0.01_wp, 0.03_wp, 0.05_wp], [3, 3])
   lattice0 = mol%lattice
   transform = identity-step*strain_direction
   mol%xyz = matmul(transform, xyz0)
   mol%lattice = matmul(transform, lattice0)
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix(mol, nmesh, kernel)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_m, vsh_mesh, energy_m)
   end select
   transform = identity+step*strain_direction
   mol%xyz = matmul(transform, xyz0)
   mol%lattice = matmul(transform, lattice0)
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix(mol, nmesh, kernel)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac, weights, &
         & density_mesh, overlap_mesh, fock_p, vsh_mesh, energy_p)
   end select
   deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
   deriv_ref = sum(sigma_mesh*strain_direction)
   err = max(err, abs(deriv_fd-deriv_ref))
   mol%xyz = xyz0
   mol%lattice = lattice0
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_bvk_Kmatrix(mol, nmesh, kernel)
   end select

   ! The real Gamma exchange gradient of the explicit 3x1x1 BvK supercell is
   ! an independent oracle for the image-kernel geometry contraction.
   allocate(mulliken_grad_sc(bas_sc%nao, bas_sc%nao), &
      & bocorr_grad_sc(bas_sc%nao, bas_sc%nao), &
      & ao_grad_sc(bas_sc%nao, bas_sc%nao), &
      & gradient_sc(3, mol_sc%nat), sigma_sc(3, 3), source=0.0_wp)
   select type (fock_sc => exchange_sc)
   type is (exchange_fock)
      call fock_sc%get_KGrad(mol_sc, ptr_sc, density_sc, overlap_sc, &
         & mulliken_grad_sc, bocorr_grad_sc, ao_grad_sc)
      call fock_sc%get_mulliken_derivs(mol_sc, ptr_sc, mulliken_grad_sc, &
         & gradient_sc, sigma_sc)
      call fock_sc%get_bocorr_derivs(mol_sc, ptr_sc, bocorr_grad_sc, &
         & gradient_sc, sigma_sc)
   end select
   do icell = 0, ncell-1
      err = max(err, maxval(abs(gradient_mesh &
         & -gradient_sc(:, 2*icell+1:2*icell+2))))
   end do
   err = max(err, maxval(abs(sigma_mesh-sigma_sc/real(ncell, wp))))

   ! A common Monkhorst-Pack twist must retain the variational Fock map, and
   ! neither the result nor the energy may depend on CP2K's k-point ordering.
   allocate(density_order(nao, nao, 1, ncell), &
      & fock_order(nao, nao, 1, ncell), fock_twist(nao, nao, 1, ncell), &
      & overlap_order(nao, nao, ncell), source=(0.0_wp, 0.0_wp))
   allocate(vsh_order(nsh), vsh_twist(nsh), source=0.0_wp)
   do ik = 1, ncell
      kfrac_twist(:, ik) = kfrac(:, ik) + [0.137_wp, -0.073_wp, 0.041_wp]
   end do
   do ik = 1, ncell
      kfrac_order(:, ik) = kfrac_twist(:, korder(ik))
      weights_order(ik) = weights(korder(ik))
      density_order(:, :, :, ik) = density_mesh(:, :, :, korder(ik))
      overlap_order(:, :, ik) = overlap_mesh(:, :, korder(ik))
   end do
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac_twist, weights, &
         & density_mesh, overlap_mesh, fock_twist, vsh_twist, energy_twist)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac_order, weights_order, &
         & density_order, overlap_order, fock_order, vsh_order, energy_order)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac_twist, weights, &
         & density_m, overlap_mesh, fock_m, vsh_mesh, energy_m)
      call fock%get_KFock_kmesh(mol, ptr, kernel, kfrac_twist, weights, &
         & density_p, overlap_mesh, fock_p, vsh_mesh, energy_p)
   end select
   deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
   deriv_ref = 0.0_wp
   order_err = 0.0_wp
   do ik = 1, ncell
      deriv_ref = deriv_ref + weights(ik)*real(sum( &
         & conjg(direction(:, :, :, ik))*fock_twist(:, :, :, ik)), wp)
      order_err = max(order_err, maxval(abs(fock_order(:, :, :, ik) &
         & -fock_twist(:, :, :, korder(ik)))))
   end do
   err = max(err, abs(deriv_fd-deriv_ref), abs(energy_order-energy_twist), &
      & maxval(abs(vsh_order-vsh_twist)), order_err)

   ! Exercise the public CP2K whole-mesh path and use its persistent plan as
   ! the cache-hit/invalidation oracle.  The first call builds the plan, an
   ! identical second call must reuse it, and a reordered mesh must rebuild
   ! the ordering descriptor without changing the physical result.
   allocate(calc_cp2k%exchange, source=exchange)
   allocate(fock_cp2k(nao, nao, 1, ncell), vsh_cp2k(nsh))
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, wfn, nmesh, &
      & kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   call view(cache_cp2k, ptr_cp2k)
   if (.not.ptr_cp2k%bvk_matches(mol, nmesh, kfrac, weights)) &
      & err = max(err, 1.0_wp)
   if (any(ptr_cp2k%bvk_input_to_grid /= [1, 2, 3])) &
      & err = max(err, 1.0_wp)
   if (any(ptr_cp2k%bvk_grid_to_input /= [1, 2, 3])) &
      & err = max(err, 1.0_wp)
   err = max(err, maxval(abs(fock_cp2k-fock_ref)), &
      & maxval(abs(vsh_cp2k-vsh_ref)), abs(energy_mesh-energy_sc))

   ! The memory-reduced stream must reject duplicate and missing input
   ! blocks without corrupting its transaction.  Completing that transaction
   ! afterwards must reproduce the unchanged whole-mesh oracle exactly.
   allocate(fock_stream(nao, nao, 1, ncell), &
      & fock_stream_order(nao, nao, 1, ncell), &
      & vsh_stream(nsh), vsh_stream_order(nsh), &
      & gradient_stream_probe(3, mol%nat), sigma_stream_probe(3, 3))
   call cp2k_exchange_stream_begin(stream, calc_cp2k, mol, cache_cp2k, wfn, &
      & nmesh, kfrac, weights, error)
   if (allocated(error)) return
   if (cp2k_exchange_stream_has_full_mesh_storage(stream)) &
      & err = max(err, 1.0_wp)
   call cp2k_exchange_stream_push(stream, 1, density_mesh(:, :, :, 1), &
      & overlap_mesh(:, :, 1), error)
   if (allocated(error)) return
   call cp2k_exchange_stream_push(stream, 1, density_mesh(:, :, :, 1), &
      & overlap_mesh(:, :, 1), error)
   if (allocated(error)) then
      if (trim(error%message) /= "duplicate g-xTB exchange stream block") then
         deallocate(error)
         call test_failed(error, "Unexpected duplicate stream-block diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Duplicate exchange stream block was accepted")
      return
   end if
   call cp2k_exchange_stream_push(stream, 2, density_mesh(:, :, :, 2), &
      & overlap_mesh(:, :, 2), error)
   if (allocated(error)) return
   call cp2k_exchange_stream_apply(stream, calc_cp2k, mol, cache_cp2k, &
      & vsh_stream, energy_stream, error)
   if (allocated(error)) then
      if (trim(error%message) /= "g-xTB exchange stream has missing blocks") then
         deallocate(error)
         call test_failed(error, "Unexpected missing stream-block diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Incomplete exchange stream was applied")
      return
   end if
   call cp2k_exchange_stream_push(stream, 3, density_mesh(:, :, :, 3), &
      & overlap_mesh(:, :, 3), error)
   if (allocated(error)) return
   onsite_probe = ptr_cp2k%g_onsfx(1, 1, 1)
   ptr_cp2k%g_onsfx(1, 1, 1) = onsite_probe + 1.0e-8_wp
   call cp2k_exchange_stream_apply(stream, calc_cp2k, mol, cache_cp2k, &
      & vsh_stream, energy_stream, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream state was invalidated before apply") then
         deallocate(error)
         call test_failed(error, "Unexpected invalidated stream diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Invalidated exchange stream was applied")
      return
   end if
   ptr_cp2k%g_onsfx(1, 1, 1) = onsite_probe
   call cp2k_exchange_stream_apply(stream, calc_cp2k, mol, cache_cp2k, &
      & vsh_stream, energy_stream, error)
   if (allocated(error)) return
   if (cp2k_exchange_stream_has_full_mesh_storage(stream)) &
      & err = max(err, 1.0_wp)
   do ik = 1, ncell
      call cp2k_exchange_stream_get(stream, ik, &
         & fock_stream(:, :, :, ik), error, overlap_mesh(:, :, ik))
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_end(stream, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(fock_stream-fock_cp2k)), &
      & maxval(abs(vsh_stream-vsh_cp2k)), abs(energy_stream-energy_mesh))

   ! Arrival order is independent of the physical k-point ordering.  Push the
   ! same indexed blocks in a nontrivial order and compare the full response.
   call cp2k_exchange_stream_begin(stream_order, calc_cp2k, mol, cache_cp2k, &
      & wfn, nmesh, kfrac, weights, error)
   if (allocated(error)) return
   if (cp2k_exchange_stream_has_full_mesh_storage(stream_order)) &
      & err = max(err, 1.0_wp)
   do j = 1, ncell
      ik = korder(j)
      call cp2k_exchange_stream_push(stream_order, ik, &
         & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_order, calc_cp2k, mol, cache_cp2k, &
      & vsh_stream_order, energy_stream_order, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_reverse_apply(stream_order, calc_cp2k, mol, &
      & cache_cp2k, gradient_stream_probe, sigma_stream_probe, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream reverse requires oracle mode") then
         deallocate(error)
         call test_failed(error, "Unexpected reduced reverse diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Reduced exchange stream accepted reverse apply")
      return
   end if
   do ik = 1, ncell
      call cp2k_exchange_stream_get(stream_order, ik, &
         & fock_stream_order(:, :, :, ik), error, overlap_mesh(:, :, ik))
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_end(stream_order, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(fock_stream_order-fock_stream)), &
      & maxval(abs(vsh_stream_order-vsh_stream)), &
      & abs(energy_stream_order-energy_stream))

   ! A common twist and a permutation of the physical k-point list must use
   ! the cached representatives without changing the reduced-stream result.
   call cp2k_exchange_stream_begin(stream_permuted, calc_cp2k, mol, &
      & cache_cp2k, wfn, nmesh, kfrac_order, weights_order, error, &
      & mode=cp2k_exchange_stream_reduced)
   if (allocated(error)) return
   if (cp2k_exchange_stream_has_full_mesh_storage(stream_permuted)) &
      & err = max(err, 1.0_wp)
   do ik = 1, ncell
      call cp2k_exchange_stream_push(stream_permuted, ik, &
         & density_order(:, :, :, ik), overlap_order(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_permuted, calc_cp2k, mol, &
      & cache_cp2k, vsh_stream_order, energy_stream_order, error)
   if (allocated(error)) return
   do ik = 1, ncell
      call cp2k_exchange_stream_get(stream_permuted, ik, &
         & fock_stream_order(:, :, :, ik), error, overlap_order(:, :, ik))
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_end(stream_permuted, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(fock_stream_order-fock_order)), &
      & maxval(abs(vsh_stream_order-vsh_order)), &
      & abs(energy_stream_order-energy_order))

   ! The cached reverse path must reproduce the independently evaluated
   ! overlap response, atomic forces, and homogeneous-strain derivative.
   allocate(overlap_grad_cp2k(nao, nao, ncell), &
      & gradient_cp2k(3, mol%nat), sigma_cp2k(3, 3))
   call cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, nmesh, &
      & kfrac, weights, density_mesh, overlap_mesh, overlap_grad_cp2k, &
      & gradient_cp2k, sigma_cp2k, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(overlap_grad_cp2k-overlap_grad_mesh)), &
      & maxval(abs(gradient_cp2k-gradient_mesh)), &
      & maxval(abs(sigma_cp2k-sigma_mesh)))

   ! Exercise the separate reverse-stream phase without changing the forward
   ! transaction tests above.  A stale onsite state must be rejected before
   ! the unchanged whole-mesh gradient oracle is entered.
   allocate(overlap_grad_stream_reverse(nao, nao, ncell), &
      & overlap_grad_stream_probe(nao, nao), &
      & gradient_stream_reverse(3, mol%nat), &
      & sigma_stream_reverse(3, 3))
   call cp2k_exchange_stream_begin(stream_reverse, calc_cp2k, mol, &
      & cache_cp2k, wfn, nmesh, kfrac, weights, error, &
      & mode=cp2k_exchange_stream_oracle)
   if (allocated(error)) return
   if (.not.cp2k_exchange_stream_has_full_mesh_storage(stream_reverse)) &
      & err = max(err, 1.0_wp)
   do ik = 1, ncell
      call cp2k_exchange_stream_push(stream_reverse, ik, &
         & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_reverse, calc_cp2k, mol, &
      & cache_cp2k, vsh_stream, energy_stream, error)
   if (allocated(error)) return
   do ik = 1, ncell
      call cp2k_exchange_stream_get(stream_reverse, ik, &
         & fock_stream_order(:, :, :, ik), error)
      if (allocated(error)) return
   end do
   err = max(err, maxval(abs(fock_stream_order-fock_cp2k)), &
      & maxval(abs(vsh_stream-vsh_cp2k)), abs(energy_stream-energy_mesh))

   onsite_probe = ptr_cp2k%g_onsfx(1, 1, 1)
   ptr_cp2k%g_onsfx(1, 1, 1) = onsite_probe + 1.0e-8_wp
   call cp2k_exchange_stream_reverse_apply(stream_reverse, calc_cp2k, mol, &
      & cache_cp2k, gradient_stream_probe, sigma_stream_probe, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream state was invalidated before reverse apply") then
         deallocate(error)
         call test_failed(error, "Unexpected stale reverse-stream diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Stale exchange reverse stream was applied")
      return
   end if
   ptr_cp2k%g_onsfx(1, 1, 1) = onsite_probe

   call cp2k_exchange_stream_reverse_apply(stream_reverse, calc_cp2k, mol, &
      & cache_cp2k, gradient_stream_reverse, sigma_stream_reverse, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_reverse_apply(stream_reverse, calc_cp2k, mol, &
      & cache_cp2k, gradient_stream_probe, sigma_stream_probe, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream reverse was already applied") then
         deallocate(error)
         call test_failed(error, "Unexpected repeated reverse-apply diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Exchange reverse stream was applied twice")
      return
   end if

   call cp2k_exchange_stream_reverse_get(stream_reverse, 1, &
      & overlap_grad_stream_reverse(:, :, 1), error)
   if (allocated(error)) return
   call cp2k_exchange_stream_reverse_get(stream_reverse, 1, &
      & overlap_grad_stream_probe, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "duplicate g-xTB exchange stream reverse pull") then
         deallocate(error)
         call test_failed(error, "Unexpected duplicate reverse-pull diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Exchange overlap adjoint was pulled twice")
      return
   end if
   do ik = 2, ncell
      call cp2k_exchange_stream_reverse_get(stream_reverse, ik, &
         & overlap_grad_stream_reverse(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_end(stream_reverse, error)
   if (allocated(error)) return
   err = max(err, &
      & maxval(abs(overlap_grad_stream_reverse-overlap_grad_cp2k)), &
      & maxval(abs(overlap_grad_stream_reverse-overlap_grad_mesh)), &
      & maxval(abs(gradient_stream_reverse-gradient_cp2k)), &
      & maxval(abs(gradient_stream_reverse-gradient_mesh)), &
      & maxval(abs(sigma_stream_reverse-sigma_cp2k)), &
      & maxval(abs(sigma_stream_reverse-sigma_mesh)))

   ! Once reverse evaluation was requested, closing the stream is legal only
   ! after every overlap-adjoint block was pulled exactly once.
   call cp2k_exchange_stream_begin(stream_reverse_missing, calc_cp2k, mol, &
      & cache_cp2k, wfn, nmesh, kfrac, weights, error, &
      & mode=cp2k_exchange_stream_oracle)
   if (allocated(error)) return
   do j = 1, ncell
      ik = korder(j)
      call cp2k_exchange_stream_push(stream_reverse_missing, ik, &
         & density_mesh(:, :, :, ik), overlap_mesh(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_reverse_missing, calc_cp2k, mol, &
      & cache_cp2k, vsh_stream, energy_stream, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_reverse_apply(stream_reverse_missing, calc_cp2k, &
      & mol, cache_cp2k, gradient_stream_probe, sigma_stream_probe, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_reverse_get(stream_reverse_missing, 1, &
      & overlap_grad_stream_probe, error)
   if (allocated(error)) return
   call cp2k_exchange_stream_end(stream_reverse_missing, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange stream ended with missing reverse pulls") then
         deallocate(error)
         call test_failed(error, "Unexpected missing reverse-pull diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Exchange stream with missing reverse pulls ended")
      return
   end if

   bvk_builds = ptr_cp2k%bvk_plan_builds
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, wfn, nmesh, &
      & kfrac, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (ptr_cp2k%bvk_plan_builds /= bvk_builds) err = max(err, 1.0_wp)

   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, wfn, nmesh, &
      & kfrac_twist, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (ptr_cp2k%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)
   err = max(err, maxval(abs(fock_cp2k-fock_twist)), &
      & maxval(abs(vsh_cp2k-vsh_twist)), abs(energy_mesh-energy_twist))

   bvk_builds = ptr_cp2k%bvk_plan_builds
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (ptr_cp2k%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)
   if (any(ptr_cp2k%bvk_input_to_grid /= korder)) err = max(err, 1.0_wp)
   if (any(ptr_cp2k%bvk_grid_to_input /= [2, 3, 1])) &
      & err = max(err, 1.0_wp)
   err = max(err, maxval(abs(fock_cp2k-fock_order)), &
      & maxval(abs(vsh_cp2k-vsh_order)), abs(energy_mesh-energy_order))

   ! Every geometry, cell, mesh, and model dependency must trigger an actual
   ! rebuild.  The AO inputs can remain fixed here because this test probes
   ! plan ownership and invalidation rather than a self-consistent state.
   bvk_builds = ptr_cp2k%bvk_plan_builds
   mol_probe = mol
   mol_probe%xyz(1, 1) = mol_probe%xyz(1, 1) + 1.0e-4_wp
   call cp2k_exchange_kmesh(calc_cp2k, mol_probe, cache_cp2k, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (ptr_cp2k%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)
   bvk_builds = ptr_cp2k%bvk_plan_builds
   mol_probe = mol
   mol_probe%lattice(1, 1) = mol_probe%lattice(1, 1) + 1.0e-4_wp
   call cp2k_exchange_kmesh(calc_cp2k, mol_probe, cache_cp2k, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (ptr_cp2k%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)

   bvk_builds = ptr_cp2k%bvk_plan_builds
   kfrac_bad = 0.0_wp
   do ik = 1, ncell
      kfrac_bad(2, ik) = real(ik-1, wp)/real(ncell, wp)
   end do
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, wfn, [1, 3, 1], &
      & kfrac_bad, weights, density_mesh, overlap_mesh, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (ptr_cp2k%bvk_plan_builds /= bvk_builds+1 &
      & .or. any(ptr_cp2k%bvk_kernel%nmesh /= [1, 3, 1])) &
      & err = max(err, 1.0_wp)

   bvk_builds = ptr_cp2k%bvk_plan_builds
   select type (fock_cp2k_model => calc_cp2k%exchange)
   type is (exchange_fock)
      fock_cp2k_model%corr_exp = fock_cp2k_model%corr_exp + 1.0e-6_wp
   end select
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (ptr_cp2k%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)
   select type (fock_cp2k_model => calc_cp2k%exchange)
   type is (exchange_fock)
      fock_cp2k_model%corr_exp = fock_cp2k_model%corr_exp - 1.0e-6_wp
   end select

   ! Reject incomplete and duplicate grids.  A boundary rejection leaves a
   ! valid plan untouched; failed internal construction invalidates it, and
   ! the next valid call must recover without incrementing on the failures.
   bvk_builds = ptr_cp2k%bvk_plan_builds
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, wfn, [4, 1, 1], &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB exchange requires a complete regular k mesh") then
         deallocate(error)
         call test_failed(error, "Unexpected incomplete-grid diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Incomplete whole-mesh exchange grid was accepted")
      return
   end if
   if (ptr_cp2k%bvk_plan_builds /= bvk_builds &
      & .or. .not.ptr_cp2k%bvk_plan_valid) err = max(err, 1.0_wp)

   kfrac_bad = kfrac_order
   kfrac_bad(:, 2) = kfrac_bad(:, 1)
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, wfn, nmesh, &
      & kfrac_bad, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) then
      if (trim(error%message) /= &
         & "g-xTB k-point mesh contains a duplicate grid point") then
         deallocate(error)
         call test_failed(error, "Unexpected duplicate-grid diagnostic")
         return
      end if
      deallocate(error)
   else
      call test_failed(error, "Duplicate whole-mesh exchange grid point was accepted")
      return
   end if
   if (ptr_cp2k%bvk_plan_builds /= bvk_builds &
      & .or. ptr_cp2k%bvk_plan_valid) err = max(err, 1.0_wp)
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, wfn, nmesh, &
      & kfrac_order, weights_order, density_order, overlap_order, fock_cp2k, &
      & vsh_cp2k, energy_mesh, error)
   if (allocated(error)) return
   if (ptr_cp2k%bvk_plan_builds /= bvk_builds+1 &
      & .or. .not.ptr_cp2k%bvk_plan_valid) err = max(err, 1.0_wp)

   ! Exercise the production O(Nk) grid proof beyond the small-mesh dense
   ! orthogonality oracle.  The cached 9x9x1 forward and reverse paths are
   ! compared with the unchanged dense-transform implementation.
   allocate(kfrac_large(3, product(nmesh_large)), &
      & weights_large(product(nmesh_large)))
   kfrac_large = 0.0_wp
   weights_large = 1.0_wp/real(product(nmesh_large), wp)
   do ik = 1, product(nmesh_large)
      kfrac_large(1, ik) = real(modulo(ik-1, nmesh_large(1)), wp) &
         & /real(nmesh_large(1), wp)
      kfrac_large(2, ik) = real((ik-1)/nmesh_large(1), wp) &
         & /real(nmesh_large(2), wp)
   end do
   allocate(density_large(nao, nao, 1, product(nmesh_large)), &
      & overlap_large(nao, nao, product(nmesh_large)), &
      & fock_large(nao, nao, 1, product(nmesh_large)), &
      & fock_large_oracle(nao, nao, 1, product(nmesh_large)), &
      & vsh_large(nsh), vsh_large_oracle(nsh))
   do ik = 1, product(nmesh_large)
      density_large(:, :, :, ik) = density_mesh(:, :, :, 1)
      overlap_large(:, :, ik) = overlap_mesh(:, :, 1)
   end do
   bvk_builds = ptr_cp2k%bvk_plan_builds
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, wfn, nmesh_large, &
      & kfrac_large, weights_large, density_large, overlap_large, fock_large, &
      & vsh_large, energy_large, error)
   if (allocated(error)) return
   if (ptr_cp2k%bvk_plan_builds /= bvk_builds+1) err = max(err, 1.0_wp)
   kernel_large = ptr_cp2k%bvk_kernel
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_KFock_kmesh(mol, ptr, kernel_large, kfrac_large, &
         & weights_large, density_large, overlap_large, fock_large_oracle, &
         & vsh_large_oracle, energy_large_oracle)
   end select
   err = max(err, maxval(abs(fock_large-fock_large_oracle)), &
      & maxval(abs(vsh_large-vsh_large_oracle)), &
      & abs(energy_large-energy_large_oracle))

   ! The production grid proof and reduced image accumulation are also tested
   ! beyond the small-mesh dense phase oracle, with reverse arrival order.
   call cp2k_exchange_stream_begin(stream_large, calc_cp2k, mol, cache_cp2k, &
      & wfn, nmesh_large, kfrac_large, weights_large, error)
   if (allocated(error)) return
   if (cp2k_exchange_stream_has_full_mesh_storage(stream_large)) &
      & err = max(err, 1.0_wp)
   do j = 1, product(nmesh_large)
      ik = product(nmesh_large)-j+1
      call cp2k_exchange_stream_push(stream_large, ik, &
         & density_large(:, :, :, ik), overlap_large(:, :, ik), error)
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_apply(stream_large, calc_cp2k, mol, cache_cp2k, &
      & vsh_large, energy_large, error)
   if (allocated(error)) return
   if (cp2k_exchange_stream_has_full_mesh_storage(stream_large)) &
      & err = max(err, 1.0_wp)
   do ik = 1, product(nmesh_large)
      call cp2k_exchange_stream_get(stream_large, ik, &
         & fock_large(:, :, :, ik), error, overlap_large(:, :, ik))
      if (allocated(error)) return
   end do
   call cp2k_exchange_stream_end(stream_large, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(fock_large-fock_large_oracle)), &
      & maxval(abs(vsh_large-vsh_large_oracle)), &
      & abs(energy_large-energy_large_oracle))

   allocate(overlap_grad_large(nao, nao, product(nmesh_large)), &
      & overlap_grad_large_oracle(nao, nao, product(nmesh_large)), &
      & mulliken_grad_large(nao, nao, product(nmesh_large)), &
      & bocorr_grad_large(nao, nao, product(nmesh_large)), &
      & gradient_large(3, mol%nat), gradient_large_oracle(3, mol%nat), &
      & sigma_large(3, 3), sigma_large_oracle(3, 3))
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_KGrad_kmesh(mol, ptr, kernel_large, kfrac_large, &
         & weights_large, density_large, overlap_large, &
         & overlap_grad_large_oracle, mulliken_grad_large, bocorr_grad_large)
      gradient_large_oracle = 0.0_wp
      sigma_large_oracle = 0.0_wp
      call fock%get_bvk_Kmatrix_derivs(mol, kernel_large, &
         & mulliken_grad_large, bocorr_grad_large, gradient_large_oracle, &
         & sigma_large_oracle)
   end select
   call cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, &
      & nmesh_large, kfrac_large, weights_large, density_large, &
      & overlap_large, overlap_grad_large, gradient_large, sigma_large, error)
   if (allocated(error)) return
   err = max(err, maxval(abs(overlap_grad_large-overlap_grad_large_oracle)), &
      & maxval(abs(gradient_large-gradient_large_oracle)), &
      & maxval(abs(sigma_large-sigma_large_oracle)))

   bvk_builds = ptr_cp2k%bvk_plan_builds
   call cp2k_exchange_kmesh(calc_cp2k, mol, cache_cp2k, wfn, nmesh_large, &
      & kfrac_large, weights_large, density_large, overlap_large, fock_large, &
      & vsh_large, energy_large, error)
   if (allocated(error)) return
   if (ptr_cp2k%bvk_plan_builds /= bvk_builds) err = max(err, 1.0_wp)
   err = max(err, maxval(abs(fock_large-fock_large_oracle)), &
      & maxval(abs(vsh_large-vsh_large_oracle)), &
      & abs(energy_large-energy_large_oracle))

   call check(error, err, 0.0_wp, thr=1.0e-9_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Whole-mesh exchange does not reproduce the explicit BvK supercell")
   end if
end subroutine test_bvk_exchange_kernel


subroutine test_kpoint_fock_complex(error)
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(basis_type), allocatable :: bas
   class(qvszp_basis_type), allocatable :: qvszp_bas
   class(exchange_type), allocatable :: exchange
   type(xtb_calculator) :: calc_cp2k
   type(basis_cache) :: bcache
   type(container_cache) :: cache, cache_cp2k
   type(exchange_cache), pointer :: ptr
   type(exchange_bvk_kernel) :: kernel_mesh
   type(wavefunction_type) :: wfn, wfn_aux
   integer, parameter :: num(3) = [1, 1, 8]
   real(wp), parameter :: xyz(3, 3) = reshape([ &
      & 0.20_wp, 0.30_wp, 0.40_wp, &
      & 0.96_wp, 0.30_wp, 0.90_wp, &
      & 0.20_wp, 0.30_wp, 1.40_wp], [3, 3])/autoaa
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      & 7.80_wp, 0.20_wp, 0.10_wp, &
      & 0.40_wp, 8.30_wp, 0.30_wp, &
      & 0.10_wp, 0.50_wp, 9.10_wp], [3, 3])/autoaa
   logical, parameter :: periodic(3) = [.true., .true., .true.]
   integer :: i, iao, iat, ii, is, ish, izp, j, jao, jat, jj, js, jsh, jzp, nao
   real(wp) :: deriv_fd, deriv_ref, energy, energy_m, energy_p, energy_ref
   real(wp) :: err, q0, step, theta
   real(wp) :: coord_direction(3, 3), identity(3, 3), lattice0(3, 3), &
      & strain_direction(3, 3), transform(3, 3), xyz0(3, 3)
   real(wp), allocatable :: lattr(:, :), overlap(:, :), vsh(:), vsh_m(:), &
      & vsh_p(:), vsh_tr(:), vsh_u(:), ao_grad_gamma(:, :), &
      & mulliken_grad(:, :), mulliken_grad_gamma(:, :), bocorr_grad(:, :), &
      & bocorr_grad_gamma(:, :), bocorr_grad_mesh_r(:, :, :), &
      & gradient_cp_k(:, :), gradient_cp_mesh(:, :), gradient_gamma(:, :), &
      & gradient_k(:, :), gradient_mesh(:, :), &
      & kernel_backup(:, :), kernel_direction(:, :), sigma_gamma(:, :), sigma_k(:, :)
   real(wp), allocatable :: mulliken_grad_mesh_r(:, :, :), &
      & sigma_cp_k(:, :), sigma_cp_mesh(:, :), sigma_mesh(:, :), vsh_mesh(:)
   complex(wp), allocatable :: density_k(:, :, :), density_m(:, :, :), &
      & density_p(:, :, :), density_u(:, :, :), direction(:, :, :), &
      & density_uks(:, :, :), &
      & fock_k(:, :, :), fock_m(:, :, :), fock_p(:, :, :), fock_tr(:, :, :), &
      & fock_u(:, :, :), fock_uks(:, :, :), overlap_grad(:, :), &
      & overlap_grad_cp_k(:, :), &
      & overlap_grad_uks(:, :), overlap_k(:, :), overlap_m(:, :), &
      & overlap_p(:, :), overlap_u(:, :), unitary(:)
   complex(wp), allocatable :: density_mesh(:, :, :, :), &
      & fock_mesh(:, :, :, :), overlap_grad_cp_mesh(:, :, :), &
      & overlap_grad_mesh(:, :, :), &
      & overlap_mesh(:, :, :)

   call new(mol, num, xyz, lattice=lattice, periodic=periodic)
   allocate(qvszp_bas)
   call make_qvszp_basis(qvszp_bas, mol, error)
   if (allocated(error)) return
   call move_alloc(qvszp_bas, bas)

   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   nao = bas%nao
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   allocate(overlap(nao, nao))
   call get_overlap(mol, lattr, cutoff, bas, bcache, overlap)
   call new_wavefunction(wfn, mol%nat, bas%nsh, nao, 1, 0.0_wp, .false.)
   do i = 1, nao
      do j = 1, nao
         wfn%density(i, j, 1) = 0.2_wp/(1.0_wp + abs(i-j))
      end do
   end do
   wfn%qsh = 0.0_wp

   call make_exchange_gxtb(exchange, mol, bas)
   call exchange%update(mol, cache)
   call view(cache, ptr)
   call exchange%get_onsite_Kmatrix(mol, wfn, ptr)
   allocate(ptr%prev_F(nao, nao, 1), ptr%prev_vsh(bas%nsh, 1))
   call exchange%get_KFock(mol, ptr, wfn%density, overlap)
   allocate(calc_cp2k%exchange, source=exchange)
   call cp2k_prepare_exchange(calc_cp2k, mol, cache_cp2k, wfn, error)
   if (allocated(error)) return

   allocate(density_k(nao, nao, 1), overlap_k(nao, nao), &
      & fock_k(nao, nao, 1), vsh(bas%nsh))
   density_k = cmplx(wfn%density, 0.0_wp, wp)
   overlap_k = cmplx(overlap, 0.0_wp, wp)
   select type (fock => exchange)
   type is (exchange_fock)
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_k, &
         & fock_k, vsh, energy)
      energy_ref = 0.5_wp*sum(wfn%density*ptr%prev_F)
      err = max(maxval(abs(real(fock_k, wp)-ptr%prev_F)), &
         & maxval(abs(aimag(fock_k))), maxval(abs(vsh-ptr%prev_vsh(:, 1))), &
         & abs(energy-energy_ref))
      call check(error, err, 0.0_wp, thr=1.0e-11_wp)
      if (allocated(error)) then
         call test_failed(error, "Complex exchange does not reproduce Gamma")
         return
      end if

      call fock%get_bvk_Kmatrix(mol, [1, 1, 1], kernel_mesh)
      allocate(density_mesh(nao, nao, 1, 1), overlap_mesh(nao, nao, 1), &
         & fock_mesh(nao, nao, 1, 1), vsh_mesh(bas%nsh))
      density_mesh(:, :, :, 1) = density_k
      overlap_mesh(:, :, 1) = overlap_k
      call fock%get_KFock_kmesh(mol, ptr, kernel_mesh, &
         & reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), [1.0_wp], &
         & density_mesh, overlap_mesh, fock_mesh, vsh_mesh, energy_ref)
      err = max(maxval(abs(fock_mesh(:, :, :, 1)-fock_k)), &
         & maxval(abs(vsh_mesh-vsh)), abs(energy_ref-energy))
      call check(error, err, 0.0_wp, thr=1.0e-11_wp)
      if (allocated(error)) then
         call test_failed(error, "Whole-mesh exchange does not reproduce Gamma")
         return
      end if

      ! At Gamma the direct overlap derivative and tblite's historical
      ! energy-weighted-density correction must obey G_S = -W_exchange.
      ! CP2K's independent-pair Pulay kernel supplies its own factor of two,
      ! so this relation fixes the sign used when folding the response.
      allocate(ao_grad_gamma(nao, nao), mulliken_grad(nao, nao), &
         & mulliken_grad_gamma(nao, nao), bocorr_grad(nao, nao), &
         & bocorr_grad_gamma(nao, nao), gradient_gamma(3, mol%nat), &
         & gradient_k(3, mol%nat), sigma_gamma(3, 3), sigma_k(3, 3), &
         & source=0.0_wp)
      allocate(overlap_grad(nao, nao), source=(0.0_wp, 0.0_wp))
      call fock%get_KGrad(mol, ptr, wfn%density, overlap, &
         & mulliken_grad_gamma, bocorr_grad_gamma, ao_grad_gamma)
      call fock%get_KGrad_kpoint(mol, ptr, density_k, overlap_k, &
         & mulliken_grad, bocorr_grad, overlap_grad)

      call fock%get_mulliken_derivs(mol, ptr, mulliken_grad_gamma, &
         & gradient_gamma, sigma_gamma)
      call fock%get_bocorr_derivs(mol, ptr, bocorr_grad_gamma, &
         & gradient_gamma, sigma_gamma)
      call fock%get_mulliken_derivs_direct(mol, ptr, mulliken_grad, &
         & gradient_k, sigma_k)
      call fock%get_bocorr_derivs_direct(mol, ptr, bocorr_grad, &
         & gradient_k, sigma_k)
      err = max(maxval(abs(overlap_grad+ao_grad_gamma)), &
         & maxval(abs(gradient_k-gradient_gamma)), &
         & maxval(abs(sigma_k-sigma_gamma)))

      ! The whole-mesh reverse mode must reduce exactly to the independently
      ! validated single-k/Gamma implementation at N=1.
      allocate(overlap_grad_mesh(nao, nao, 1), &
         & mulliken_grad_mesh_r(nao, nao, 1), &
         & bocorr_grad_mesh_r(nao, nao, 1), gradient_mesh(3, mol%nat), &
         & sigma_mesh(3, 3))
      call fock%get_KGrad_kmesh(mol, ptr, kernel_mesh, &
         & reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), [1.0_wp], &
         & density_mesh, overlap_mesh, overlap_grad_mesh, &
         & mulliken_grad_mesh_r, bocorr_grad_mesh_r)
      gradient_mesh = 0.0_wp
      sigma_mesh = 0.0_wp
      call fock%get_bvk_Kmatrix_derivs(mol, kernel_mesh, &
         & mulliken_grad_mesh_r, bocorr_grad_mesh_r, gradient_mesh, sigma_mesh)
      err = max(err, maxval(abs(overlap_grad_mesh(:, :, 1)-overlap_grad)), &
         & maxval(abs(gradient_mesh-gradient_k)), &
         & maxval(abs(sigma_mesh-sigma_k)))

      ! Exercise the public CP2K ABI as well: both wrappers must expose the
      ! same N=1 functional, including overlap sign and strain convention.
      allocate(overlap_grad_cp_k(nao, nao), &
         & overlap_grad_cp_mesh(nao, nao, 1), &
         & gradient_cp_k(3, mol%nat), gradient_cp_mesh(3, mol%nat), &
         & sigma_cp_k(3, 3), sigma_cp_mesh(3, 3))
      call cp2k_exchange_kpoint_gradient(calc_cp2k, mol, cache_cp2k, &
         & density_k, overlap_k, overlap_grad_cp_k, gradient_cp_k, &
         & sigma_cp_k, error)
      if (allocated(error)) return
      call cp2k_exchange_kmesh_gradient(calc_cp2k, mol, cache_cp2k, &
         & [1, 1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), &
         & [1.0_wp], density_mesh, overlap_mesh, overlap_grad_cp_mesh, &
         & gradient_cp_mesh, sigma_cp_mesh, error)
      if (allocated(error)) return
      err = max(err, maxval(abs(overlap_grad_cp_k-overlap_grad)), &
         & maxval(abs(overlap_grad_cp_mesh(:, :, 1)-overlap_grad_cp_k)), &
         & maxval(abs(gradient_cp_mesh-gradient_cp_k)), &
         & maxval(abs(sigma_cp_mesh-sigma_cp_k)))
      call check(error, err, 0.0_wp, thr=1.0e-10_wp)
      if (allocated(error)) then
         call test_failed(error, "Complex exchange gradient does not reproduce Gamma")
         return
      end if

      do i = 1, nao
         do j = 1, i-1
            overlap_k(j, i) = overlap_k(j, i) &
               & + cmplx(0.0_wp, 1.0e-3_wp*sin(real(3*i+5*j, wp)), wp)
            overlap_k(i, j) = conjg(overlap_k(j, i))
            density_k(j, i, 1) = density_k(j, i, 1) &
               & + cmplx(0.0_wp, 2.0e-3_wp*cos(real(7*i+2*j, wp)), wp)
            density_k(i, j, 1) = conjg(density_k(j, i, 1))
         end do
      end do
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_k, &
         & fock_k, vsh, energy)
      err = maxval(abs(fock_k(:, :, 1)-conjg(transpose(fock_k(:, :, 1)))))

      allocate(fock_tr(nao, nao, 1), vsh_tr(bas%nsh))
      call fock%get_KFock_kpoint(mol, ptr, conjg(density_k), &
         & conjg(overlap_k), fock_tr, vsh_tr, energy_ref)
      err = max(err, maxval(abs(fock_tr-conjg(fock_k))), &
         & maxval(abs(vsh_tr-vsh)), abs(energy_ref-energy))

      allocate(unitary(nao), density_u(nao, nao, 1), overlap_u(nao, nao), &
         & fock_u(nao, nao, 1), vsh_u(bas%nsh))
      do i = 1, nao
         theta = 0.071_wp*real(i, wp)
         unitary(i) = exp(cmplx(0.0_wp, theta, wp))
      end do
      do i = 1, nao
         do j = 1, nao
            overlap_u(i, j) = conjg(unitary(i))*overlap_k(i, j)*unitary(j)
            density_u(i, j, 1) = conjg(unitary(i)) &
               & * density_k(i, j, 1)*unitary(j)
         end do
      end do
      call fock%get_KFock_kpoint(mol, ptr, density_u, overlap_u, &
         & fock_u, vsh_u, energy_ref)
      do i = 1, nao
         do j = 1, nao
            err = max(err, abs(fock_u(i, j, 1) &
               & - conjg(unitary(i))*fock_k(i, j, 1)*unitary(j)))
         end do
      end do
      err = max(err, maxval(abs(vsh_u-vsh)), abs(energy_ref-energy))

      allocate(direction(nao, nao, 1), density_m(nao, nao, 1), &
         & density_p(nao, nao, 1), fock_m(nao, nao, 1), &
         & fock_p(nao, nao, 1), source=(0.0_wp, 0.0_wp))
      allocate(vsh_m(bas%nsh), vsh_p(bas%nsh))
      do i = 1, nao
         direction(i, i, 1) = cmplx(1.0e-2_wp*sin(real(11*i, wp)), &
            & 0.0_wp, wp)
         do j = 1, i-1
            direction(j, i, 1) = cmplx( &
               & 1.0e-2_wp*cos(real(3*i+j, wp)), &
               & 1.0e-2_wp*sin(real(i+4*j, wp)), wp)
            direction(i, j, 1) = conjg(direction(j, i, 1))
         end do
      end do
      step = 1.0e-5_wp
      density_m = density_k-step*direction
      density_p = density_k+step*direction
      call fock%get_KFock_kpoint(mol, ptr, density_m, overlap_k, &
         & fock_m, vsh_m, energy_m)
      call fock%get_KFock_kpoint(mol, ptr, density_p, overlap_k, &
         & fock_p, vsh_p, energy_p)
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = real(sum(conjg(direction)*fock_k), wp)
      err = max(err, abs(deriv_fd-deriv_ref))

      ! Reverse the complex exchange functional with respect to S(k).  The
      ! overlap response is the missing Pulay ingredient of a k-point force;
      ! validate it independently of any real-space integral convention.
      allocate(overlap_m(nao, nao), overlap_p(nao, nao))
      call fock%get_KGrad_kpoint(mol, ptr, density_k, overlap_k, &
         & mulliken_grad, bocorr_grad, overlap_grad)
      overlap_m = overlap_k-step*direction(:, :, 1)
      overlap_p = overlap_k+step*direction(:, :, 1)
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_m, &
         & fock_m, vsh_m, energy_m)
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_p, &
         & fock_p, vsh_p, energy_p)
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = real(sum(conjg(overlap_grad)*direction(:, :, 1)), wp)
      err = max(err, abs(deriv_fd-deriv_ref))
      err = max(err, maxval(abs(overlap_grad &
         & - conjg(transpose(overlap_grad)))))

      ! The geometry-dependent Mulliken kernel is real and symmetric.  Check
      ! the AO-resolved reverse response before it is contracted with the
      ! pairwise kernel derivative in get_mulliken_derivs.
      allocate(kernel_backup(fock%nsh, fock%nsh), &
         & kernel_direction(fock%nsh, fock%nsh), source=0.0_wp)
      do ish = 1, fock%nsh
         kernel_direction(ish, ish) = 1.0e-2_wp*sin(real(13*ish, wp))
         do jsh = 1, ish-1
            kernel_direction(jsh, ish) = &
               & 1.0e-2_wp*cos(real(5*ish+7*jsh, wp))
            kernel_direction(ish, jsh) = kernel_direction(jsh, ish)
         end do
      end do
      kernel_backup = ptr%g_mulliken
      ptr%g_mulliken = kernel_backup-step*kernel_direction
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_k, &
         & fock_m, vsh_m, energy_m)
      ptr%g_mulliken = kernel_backup+step*kernel_direction
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_k, &
         & fock_p, vsh_p, energy_p)
      ptr%g_mulliken = kernel_backup
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = 0.0_wp
      do ish = 1, fock%nsh
         ii = fock%iao_sh(ish)
         do jsh = 1, ish
            jj = fock%iao_sh(jsh)
            do iao = 1, fock%nao_sh(ish)
               do jao = 1, fock%nao_sh(jsh)
                  deriv_ref = deriv_ref + kernel_direction(jsh, ish) &
                     & * mulliken_grad(jj+jao, ii+iao)
               end do
            end do
         end do
      end do
      err = max(err, abs(deriv_fd-deriv_ref))

      ! Repeat the same directional check for the atom-block bond-order
      ! correction kernel.  As above, the response stores one derivative for
      ! each independent symmetric block, so only one triangle is contracted.
      deallocate(kernel_backup, kernel_direction)
      allocate(kernel_backup(mol%nat, mol%nat), &
         & kernel_direction(mol%nat, mol%nat), source=0.0_wp)
      do iat = 1, mol%nat
         do jat = 1, iat-1
            kernel_direction(jat, iat) = &
               & 1.0e-2_wp*sin(real(3*iat+11*jat, wp))
            kernel_direction(iat, jat) = kernel_direction(jat, iat)
         end do
      end do
      kernel_backup = ptr%g_bocorr
      ptr%g_bocorr = kernel_backup-step*kernel_direction
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_k, &
         & fock_m, vsh_m, energy_m)
      ptr%g_bocorr = kernel_backup+step*kernel_direction
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_k, &
         & fock_p, vsh_p, energy_p)
      ptr%g_bocorr = kernel_backup
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = 0.0_wp
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = fock%ish_at(iat)
         do jat = 1, iat-1
            jzp = mol%id(jat)
            js = fock%ish_at(jat)
            do ish = 1, fock%nsh_id(izp)
               ii = fock%iao_sh(is+ish)
               do jsh = 1, fock%nsh_id(jzp)
                  jj = fock%iao_sh(js+jsh)
                  do iao = 1, fock%nao_sh(is+ish)
                     do jao = 1, fock%nao_sh(js+jsh)
                        deriv_ref = deriv_ref + kernel_direction(jat, iat) &
                           & * bocorr_grad(jj+jao, ii+iao)
                     end do
                  end do
               end do
            end do
         end do
      end do
      err = max(err, abs(deriv_fd-deriv_ref))

      ! Contract the direct kernel responses with their real-space derivatives
      ! and validate the complete explicit exchange geometry response for a
      ! genuinely complex Hermitian k-point density.  The overlap is held
      ! fixed here, so the separate overlap-response check above remains the
      ! sole Pulay contribution.
      gradient_k = 0.0_wp
      sigma_k = 0.0_wp
      call fock%get_mulliken_derivs_direct(mol, ptr, mulliken_grad, &
         & gradient_k, sigma_k)
      call fock%get_bocorr_derivs_direct(mol, ptr, bocorr_grad, &
         & gradient_k, sigma_k)

      do iat = 1, mol%nat
         do i = 1, 3
            coord_direction(i, iat) = 0.03_wp*sin(real(7*i+11*iat, wp))
         end do
      end do
      xyz0 = mol%xyz
      mol%xyz = xyz0-step*coord_direction
      call ptr%update(mol)
      call fock%get_mulliken_Kmatrix(mol, ptr)
      call fock%get_bocorr_Kmatrix(mol, ptr)
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_k, &
         & fock_m, vsh_m, energy_m)
      mol%xyz = xyz0+step*coord_direction
      call ptr%update(mol)
      call fock%get_mulliken_Kmatrix(mol, ptr)
      call fock%get_bocorr_Kmatrix(mol, ptr)
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_k, &
         & fock_p, vsh_p, energy_p)
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = sum(gradient_k*coord_direction)
      err = max(err, abs(deriv_fd-deriv_ref))

      ! Apply a symmetric homogeneous strain to both Cartesian coordinates and
      ! lattice vectors.  This also covers diagonal self-image kernels, whose
      ! historical derivative representation carries a factor of two.
      identity = 0.0_wp
      do i = 1, 3
         identity(i, i) = 1.0_wp
      end do
      strain_direction = reshape([ &
         & 0.07_wp, 0.02_wp, -0.01_wp, &
         & 0.02_wp, -0.04_wp, 0.03_wp, &
         & -0.01_wp, 0.03_wp, 0.05_wp], [3, 3])
      lattice0 = mol%lattice
      mol%xyz = xyz0
      transform = identity-step*strain_direction
      mol%xyz = matmul(transform, xyz0)
      mol%lattice = matmul(transform, lattice0)
      call ptr%update(mol)
      call fock%get_mulliken_Kmatrix(mol, ptr)
      call fock%get_bocorr_Kmatrix(mol, ptr)
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_k, &
         & fock_m, vsh_m, energy_m)
      transform = identity+step*strain_direction
      mol%xyz = matmul(transform, xyz0)
      mol%lattice = matmul(transform, lattice0)
      call ptr%update(mol)
      call fock%get_mulliken_Kmatrix(mol, ptr)
      call fock%get_bocorr_Kmatrix(mol, ptr)
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_k, &
         & fock_p, vsh_p, energy_p)
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = sum(sigma_k*strain_direction)
      err = max(err, abs(deriv_fd-deriv_ref))

      mol%xyz = xyz0
      mol%lattice = lattice0
      call ptr%update(mol)
      call fock%get_mulliken_Kmatrix(mol, ptr)
      call fock%get_bocorr_Kmatrix(mol, ptr)

      ! The shell potential is the derivative of the same complex-k energy
      ! with respect to the charge-dependent onsite exchange kernels.
      q0 = wfn%qsh(1, 1)
      wfn%qsh(1, 1) = q0-step
      call exchange%get_onsite_Kmatrix(mol, wfn, ptr)
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_k, &
         & fock_m, vsh_m, energy_m)
      wfn%qsh(1, 1) = q0+step
      call exchange%get_onsite_Kmatrix(mol, wfn, ptr)
      call fock%get_KFock_kpoint(mol, ptr, density_k, overlap_k, &
         & fock_p, vsh_p, energy_p)
      err = max(err, abs((energy_p-energy_m)/(2.0_wp*step)-vsh(1)))
      wfn%qsh(1, 1) = q0
      call exchange%get_onsite_Kmatrix(mol, wfn, ptr)

      ! Exercise the unrestricted spin factor and accumulation of both spin
      ! channels in the common overlap response.
      allocate(density_uks(nao, nao, 2), fock_uks(nao, nao, 2), &
         & overlap_grad_uks(nao, nao))
      density_uks(:, :, 1) = 0.61_wp*density_k(:, :, 1)
      density_uks(:, :, 2) = 0.37_wp*density_k(:, :, 1)
      call fock%get_KFock_kpoint(mol, ptr, density_uks, overlap_k, &
         & fock_uks, vsh, energy)
      call fock%get_KGrad_kpoint(mol, ptr, density_uks, overlap_k, &
         & mulliken_grad, bocorr_grad, overlap_grad_uks)
      overlap_m = overlap_k-step*direction(:, :, 1)
      overlap_p = overlap_k+step*direction(:, :, 1)
      call fock%get_KFock_kpoint(mol, ptr, density_uks, overlap_m, &
         & fock_uks, vsh_m, energy_m)
      call fock%get_KFock_kpoint(mol, ptr, density_uks, overlap_p, &
         & fock_uks, vsh_p, energy_p)
      deriv_fd = (energy_p-energy_m)/(2.0_wp*step)
      deriv_ref = real(sum(conjg(overlap_grad_uks) &
         & * direction(:, :, 1)), wp)
      err = max(err, abs(deriv_fd-deriv_ref))
   class default
      call test_failed(error, "Unexpected exchange implementation")
      return
   end select

   call check(error, err, 0.0_wp, thr=1.0e-9_wp)
   if (allocated(error)) then
      call test_failed(error, &
         & "Complex exchange violates Hermiticity, covariance, or variationality")
   end if
end subroutine test_kpoint_fock_complex


end module test_exchange
