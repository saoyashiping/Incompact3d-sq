!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_interp

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use variables, only : nx, ny, nz
  use param, only : xlx, yly, zlz
  use fiber_types, only : fiber_active, fiber_nl, fiber_x, fiber_uinterp, fiber_uexact, fiber_uerror, &
       fiber_sumw, fiber_interp_max_error, interp_test_case
  use fiber_delta, only : delta_kernel_3d

  implicit none

contains

  subroutine interp_euler_to_lagrangian_velocity(uxe, uye, uze, xg, yg, zg)

    real(mytype), intent(in), dimension(:,:,:) :: uxe, uye, uze
    real(mytype), intent(in), dimension(:) :: xg, yg, zg

    integer :: i, j, k, l
    real(mytype) :: hx, hy, hz, weight, sumw
    real(mytype) :: rx, ry, rz

    if (.not.fiber_active) return

    hx = xg(2) - xg(1)
    if (size(yg) > 1) then
      hy = yg(2) - yg(1)
    else
      hy = 1._mytype
    endif
    hz = zg(2) - zg(1)

    if (.not.allocated(fiber_uinterp)) allocate(fiber_uinterp(3, fiber_nl))
    if (.not.allocated(fiber_sumw)) allocate(fiber_sumw(fiber_nl))
    fiber_uinterp = 0._mytype
    fiber_sumw = 0._mytype

    do l = 1, fiber_nl
      sumw = 0._mytype

      do k = 1, size(zg)
        rz = fiber_x(3,l) - zg(k)
        if (abs(rz) > 2._mytype * hz) cycle

        do j = 1, size(yg)
          ry = fiber_x(2,l) - yg(j)
          if (abs(ry) > 2._mytype * hy) cycle

          do i = 1, size(xg)
            rx = fiber_x(1,l) - xg(i)
            if (abs(rx) > 2._mytype * hx) cycle

            weight = delta_kernel_3d(rx, ry, rz, hx, hy, hz) * hx * hy * hz
            sumw = sumw + weight
            fiber_uinterp(1,l) = fiber_uinterp(1,l) + uxe(i,j,k) * weight
            fiber_uinterp(2,l) = fiber_uinterp(2,l) + uye(i,j,k) * weight
            fiber_uinterp(3,l) = fiber_uinterp(3,l) + uze(i,j,k) * weight
          enddo
        enddo
      enddo

      fiber_sumw(l) = sumw
      if (sumw > 0._mytype) then
        fiber_uinterp(:,l) = fiber_uinterp(:,l) / sumw
      endif
    enddo

  end subroutine interp_euler_to_lagrangian_velocity

  subroutine run_fiber_interp_operator_test()

    integer :: i, j, k, l
    real(mytype), allocatable, dimension(:) :: xg, yg, zg
    real(mytype), allocatable, dimension(:,:,:) :: uxe, uye, uze
    real(mytype) :: hx, hy, hz
    real(mytype) :: c1, c2, c3
    real(mytype) :: a1, b1, c_1, d1
    real(mytype) :: a2, b2, c_2, d2
    real(mytype) :: a3, b3, c_3, d3

    if (.not.fiber_active) then
      if (nrank == 0) write(*,*) 'Error: interp_test_active requires fiber_active = true.'
      stop
    endif

    allocate(xg(nx), yg(ny), zg(nz))
    allocate(uxe(nx,ny,nz), uye(nx,ny,nz), uze(nx,ny,nz))

    hx = xlx / real(nx, mytype)
    if (ny > 1) then
      hy = yly / real(ny - 1, mytype)
    else
      hy = 1._mytype
    endif
    hz = zlz / real(nz, mytype)

    do i = 1, nx
      xg(i) = real(i - 1, mytype) * hx
    enddo
    do j = 1, ny
      yg(j) = real(j - 1, mytype) * hy
    enddo
    do k = 1, nz
      zg(k) = real(k - 1, mytype) * hz
    enddo

    if (interp_test_case == 1) then
      c1 = 1.25_mytype
      c2 = -0.5_mytype
      c3 = 0.75_mytype

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            uxe(i,j,k) = c1
            uye(i,j,k) = c2
            uze(i,j,k) = c3
          enddo
        enddo
      enddo

      call interp_euler_to_lagrangian_velocity(uxe, uye, uze, xg, yg, zg)

      if (.not.allocated(fiber_uexact)) allocate(fiber_uexact(3, fiber_nl))
      if (.not.allocated(fiber_uerror)) allocate(fiber_uerror(3, fiber_nl))

      do l = 1, fiber_nl
        fiber_uexact(:,l) = (/c1, c2, c3/)
      enddo

    else if (interp_test_case == 2) then
      a1 = 0.3_mytype; b1 = -0.2_mytype; c_1 = 0.1_mytype; d1 = 0.5_mytype
      a2 = -0.1_mytype; b2 = 0.25_mytype; c_2 = 0.05_mytype; d2 = -0.3_mytype
      a3 = 0.2_mytype; b3 = 0.1_mytype; c_3 = -0.15_mytype; d3 = 0.1_mytype

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            uxe(i,j,k) = a1 * xg(i) + b1 * yg(j) + c_1 * zg(k) + d1
            uye(i,j,k) = a2 * xg(i) + b2 * yg(j) + c_2 * zg(k) + d2
            uze(i,j,k) = a3 * xg(i) + b3 * yg(j) + c_3 * zg(k) + d3
          enddo
        enddo
      enddo

      call interp_euler_to_lagrangian_velocity(uxe, uye, uze, xg, yg, zg)

      if (.not.allocated(fiber_uexact)) allocate(fiber_uexact(3, fiber_nl))
      if (.not.allocated(fiber_uerror)) allocate(fiber_uerror(3, fiber_nl))

      do l = 1, fiber_nl
        fiber_uexact(1,l) = a1 * fiber_x(1,l) + b1 * fiber_x(2,l) + c_1 * fiber_x(3,l) + d1
        fiber_uexact(2,l) = a2 * fiber_x(1,l) + b2 * fiber_x(2,l) + c_2 * fiber_x(3,l) + d2
        fiber_uexact(3,l) = a3 * fiber_x(1,l) + b3 * fiber_x(2,l) + c_3 * fiber_x(3,l) + d3
      enddo

    else
      if (nrank == 0) write(*,*) 'Error: interp_test_case must be 1 (constant) or 2 (linear).'
      stop
    endif

    fiber_uerror = fiber_uinterp - fiber_uexact
    fiber_interp_max_error = maxval(abs(fiber_uerror))

    if (nrank == 0) then
      if (interp_test_case == 1) write(*,'(A,ES12.4)') 'Fiber interpolation constant-field max abs error: ', fiber_interp_max_error
      if (interp_test_case == 2) write(*,'(A,ES12.4)') 'Fiber interpolation linear-field max abs error  : ', fiber_interp_max_error
      write(*,'(A,ES12.4)') 'Fiber interpolation sumw_min                   : ', minval(fiber_sumw)
      write(*,'(A,ES12.4)') 'Fiber interpolation sumw_max                   : ', maxval(fiber_sumw)
    endif

    deallocate(xg, yg, zg, uxe, uye, uze)

  end subroutine run_fiber_interp_operator_test

end module fiber_interp
