!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_spread

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use variables, only : nx, ny, nz
  use param, only : xlx, yly, zlz
  use fiber_types, only : fiber_active, fiber_nl, fiber_x, fiber_quad_w, fiber_test_force, spread_test_case
  use fiber_delta, only : delta_kernel_3d

  implicit none

contains

  pure real(mytype) function minimum_image(delta, period_length)

    real(mytype), intent(in) :: delta, period_length

    minimum_image = delta
    if (period_length > 0._mytype) then
      minimum_image = delta - period_length * nint(delta / period_length)
    endif

  end function minimum_image

  subroutine run_fiber_spread_conservation_test(lag_total, eul_total, abs_error, rel_error, sumw_min, sumw_max)

    real(mytype), intent(out), dimension(3) :: lag_total, eul_total, abs_error, rel_error
    real(mytype), intent(out) :: sumw_min, sumw_max

    real(mytype), allocatable, dimension(:,:,:) :: fx, fy, fz
    real(mytype), allocatable, dimension(:) :: xg, yg, zg, sumw_l

    integer :: i, j, k, l
    real(mytype) :: hx, hy, hz, dV, weight
    real(mytype) :: rx, ry, rz

    if (.not.fiber_active) then
      if (nrank == 0) write(*,*) 'Error: spread_test_active requires fiber_active = true.'
      stop
    endif

    if (.not.allocated(fiber_quad_w)) then
      if (nrank == 0) write(*,*) 'Error: fiber_quad_w is not allocated.'
      stop
    endif

    if (.not.allocated(fiber_test_force)) then
      allocate(fiber_test_force(3, fiber_nl))
    endif

    fiber_test_force = 0._mytype
    do l = 1, fiber_nl
      fiber_test_force(:,l) = (/1._mytype, 0.5_mytype, -0.25_mytype/)
    enddo

    allocate(xg(nx), yg(ny), zg(nz))
    allocate(fx(nx,ny,nz), fy(nx,ny,nz), fz(nx,ny,nz))
    allocate(sumw_l(fiber_nl))

    hx = xlx / real(nx, mytype)
    if (ny > 1) then
      hy = yly / real(ny - 1, mytype)
    else
      hy = 1._mytype
    endif
    hz = zlz / real(nz, mytype)
    dV = hx * hy * hz

    do i = 1, nx
      xg(i) = real(i - 1, mytype) * hx
    enddo
    do j = 1, ny
      yg(j) = real(j - 1, mytype) * hy
    enddo
    do k = 1, nz
      zg(k) = real(k - 1, mytype) * hz
    enddo

    fx = 0._mytype
    fy = 0._mytype
    fz = 0._mytype
    sumw_l = 0._mytype

    do l = 1, fiber_nl
      do k = 1, nz
        rz = minimum_image(zg(k) - fiber_x(3,l), zlz)
        if (abs(rz) > 2._mytype * hz) cycle

        do j = 1, ny
          ry = yg(j) - fiber_x(2,l)
          if (abs(ry) > 2._mytype * hy) cycle

          do i = 1, nx
            rx = minimum_image(xg(i) - fiber_x(1,l), xlx)
            if (abs(rx) > 2._mytype * hx) cycle

            weight = delta_kernel_3d(rx, ry, rz, hx, hy, hz)

            fx(i,j,k) = fx(i,j,k) + fiber_test_force(1,l) * fiber_quad_w(l) * weight
            fy(i,j,k) = fy(i,j,k) + fiber_test_force(2,l) * fiber_quad_w(l) * weight
            fz(i,j,k) = fz(i,j,k) + fiber_test_force(3,l) * fiber_quad_w(l) * weight

            sumw_l(l) = sumw_l(l) + weight * dV
          enddo
        enddo
      enddo
    enddo

    lag_total(1) = sum(fiber_test_force(1,:) * fiber_quad_w(:))
    lag_total(2) = sum(fiber_test_force(2,:) * fiber_quad_w(:))
    lag_total(3) = sum(fiber_test_force(3,:) * fiber_quad_w(:))

    eul_total(1) = sum(fx) * dV
    eul_total(2) = sum(fy) * dV
    eul_total(3) = sum(fz) * dV

    abs_error = abs(eul_total - lag_total)

    rel_error = 0._mytype
    do i = 1, 3
      if (abs(lag_total(i)) > 0._mytype) then
        rel_error(i) = abs_error(i) / abs(lag_total(i))
      endif
    enddo

    sumw_min = minval(sumw_l)
    sumw_max = maxval(sumw_l)

    if (nrank == 0) then
      write(*,'(A,I6)') 'Fiber spread conservation test case: ', spread_test_case
      write(*,'(A,3ES12.4)') 'Lagrangian total force          : ', lag_total
      write(*,'(A,3ES12.4)') 'Euler integrated spread force   : ', eul_total
      write(*,'(A,3ES12.4)') 'Absolute conservation error     : ', abs_error
      write(*,'(A,3ES12.4)') 'Relative conservation error     : ', rel_error
      write(*,'(A,2ES12.4)') 'Spread kernel sumw min/max      : ', sumw_min, sumw_max
    endif

    deallocate(xg, yg, zg, fx, fy, fz, sumw_l)

  end subroutine run_fiber_spread_conservation_test

end module fiber_spread
