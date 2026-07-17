program benchmark_separable_dft
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use tblite_exchange_transform, only : bvk_separable_k_to_r, &
      & bvk_separable_r_to_k
   implicit none

   integer, parameter :: nrow = 64, ncase = 6, repeat = 20
   integer, parameter :: meshes(3, ncase) = reshape([ &
      & 216, 1, 1, 18, 12, 1, 6, 6, 6, &
      & 512, 1, 1, 32, 16, 1, 8, 8, 8], [3, ncase])
   integer :: c, g(3), icell, igrid, ik, j, nmesh(3), nk, rep, root
   integer :: count0, count1, rate
   integer, allocatable :: input_to_grid(:), reps(:, :)
   real(wp), parameter :: twist(3) = [0.037_wp, 0.021_wp, 0.013_wp]
   real(wp) :: angle, dense_time, err, sep_time
   complex(wp) :: phase
   complex(wp), allocatable :: dense_k(:, :), dense_r(:, :), &
      & input_k(:, :), input_r(:, :), phase_forward(:, :), &
      & phase_inverse(:, :), phase_roots(:, :), sep_k(:, :), sep_r(:, :), &
      & twist_phase(:)

   call system_clock(count_rate=rate)
   write(*, '(a)') 'n1 n2 n3 nk rows dense_s separable_dft_s speedup max_error'
   do c = 1, ncase
      nmesh = meshes(:, c)
      nk = product(nmesh)
      allocate(input_to_grid(nk), reps(3, nk), input_k(nrow, nk), &
         & input_r(nrow, nk), dense_k(nrow, nk), dense_r(nrow, nk), &
         & sep_k(nrow, nk), sep_r(nrow, nk), phase_forward(nk, nk), &
         & phase_inverse(nk, nk), phase_roots(maxval(nmesh), 3), &
         & twist_phase(nk))
      do ik = 1, nk
         input_to_grid(ik) = 1+modulo(17*(ik-1)+5, nk)
         do j = 1, nrow
            input_k(j, ik) = cmplx(sin(real(7*j+11*ik, wp)), &
               & cos(real(13*j+5*ik, wp)), wp)
            input_r(j, ik) = cmplx(cos(real(3*j+17*ik, wp)), &
               & sin(real(19*j+2*ik, wp)), wp)
         end do
      end do
      do icell = 1, nk
         igrid = nk-icell
         reps(1, icell) = modulo(igrid, nmesh(1))
         reps(2, icell) = modulo(igrid/nmesh(1), nmesh(2))
         reps(3, icell) = igrid/(nmesh(1)*nmesh(2))
         angle = 2.0_wp*pi*dot_product(twist, real(reps(:, icell), wp))
         twist_phase(icell) = exp(cmplx(0.0_wp, angle, wp))
      end do
      phase_roots = (1.0_wp, 0.0_wp)
      do j = 1, 3
         do root = 0, nmesh(j)-1
            angle = 2.0_wp*pi*real(root, wp)/real(nmesh(j), wp)
            phase_roots(root+1, j) = exp(cmplx(0.0_wp, angle, wp))
         end do
      end do
      do ik = 1, nk
         igrid = input_to_grid(ik)-1
         g(1) = modulo(igrid, nmesh(1))
         g(2) = modulo(igrid/nmesh(1), nmesh(2))
         g(3) = igrid/(nmesh(1)*nmesh(2))
         do icell = 1, nk
            angle = 2.0_wp*pi*dot_product(twist+real(g, wp)/real(nmesh, wp), &
               & real(reps(:, icell), wp))
            phase = exp(cmplx(0.0_wp, angle, wp))
            phase_forward(icell, ik) = phase
            phase_inverse(ik, icell) = conjg(phase)/real(nk, wp)
         end do
      end do

      call system_clock(count0)
      do rep = 1, repeat
         dense_r = matmul(input_k, phase_inverse)
         dense_k = matmul(input_r, phase_forward)
      end do
      call system_clock(count1)
      dense_time = real(count1-count0, wp)/real(rate*repeat, wp)

      call system_clock(count0)
      do rep = 1, repeat
         call bvk_separable_k_to_r(input_k, sep_r, nmesh, input_to_grid, &
            & reps, phase_roots, twist_phase, .true.)
         call bvk_separable_r_to_k(input_r, sep_k, nmesh, input_to_grid, &
            & reps, phase_roots, twist_phase)
      end do
      call system_clock(count1)
      sep_time = real(count1-count0, wp)/real(rate*repeat, wp)
      err = max(maxval(abs(sep_r-dense_r)), maxval(abs(sep_k-dense_k)))
      write(*, '(3(i4,1x),i6,1x,i5,1x,f10.6,1x,f15.6,1x,f8.3,1x,es12.4)') &
         & nmesh, nk, nrow, dense_time, sep_time, dense_time/sep_time, err
      deallocate(input_to_grid, reps, input_k, input_r, dense_k, dense_r, &
         & sep_k, sep_r, phase_forward, phase_inverse, phase_roots, twist_phase)
   end do
end program benchmark_separable_dft
