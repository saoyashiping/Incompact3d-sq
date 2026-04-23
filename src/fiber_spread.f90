!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_spread

  use decomp_2d_constants, only : mytype
  use decomp_2d, only : xstart
  use decomp_2d_mpi, only : nrank
  use mpi, only : MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MIN, MPI_MAX, MPI_COMM_WORLD, MPI_IN_PLACE
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

  pure real(mytype) function local_dy(yp_vec, j, jmin, jmax)

    real(mytype), intent(in), dimension(:) :: yp_vec
    integer, intent(in) :: j, jmin, jmax

    if (j <= jmin) then
      local_dy = yp_vec(jmin+1) - yp_vec(jmin)
    else if (j >= jmax) then
      local_dy = yp_vec(jmax) - yp_vec(jmax-1)
    else
      local_dy = 0.5_mytype * (yp_vec(j+1) - yp_vec(j-1))
    endif

  end function local_dy

  subroutine spread_lagrangian_force_to_euler(fiber_force_lag, fx, fy, fz, eul_total, sumw_min, sumw_max, &
       xg_in, yg_in, zg_in, hy_loc_min_out, hy_loc_max_out)

    real(mytype), intent(in), dimension(3, fiber_nl) :: fiber_force_lag
    real(mytype), intent(out), dimension(:,:,:) :: fx, fy, fz
    real(mytype), intent(out), dimension(3) :: eul_total
    real(mytype), intent(out) :: sumw_min, sumw_max
    real(mytype), intent(in), optional, dimension(:) :: xg_in, yg_in, zg_in
    real(mytype), intent(out), optional :: hy_loc_min_out, hy_loc_max_out

    real(mytype), allocatable, dimension(:) :: xg, yg, zg, sumw_l, sumw_l_global
    integer :: i, j, k, l, ierr
    integer :: nx_loc, ny_loc, nz_loc
    integer :: ig, jg, kg, shift_x, shift_y, shift_z
    integer :: jmin_g, jmax_g
    real(mytype) :: hx, hy, hz, hy_loc, dVloc, weight
    real(mytype) :: rx, ry, rz
    real(mytype) :: hy_loc_min, hy_loc_max
    logical :: use_input_grid

    use_input_grid = present(xg_in) .and. present(yg_in) .and. present(zg_in)

    nx_loc = size(fx,1)
    ny_loc = size(fx,2)
    nz_loc = size(fx,3)

    allocate(xg(nx_loc), yg(ny_loc), zg(nz_loc), sumw_l(fiber_nl), sumw_l_global(fiber_nl))

    if (use_input_grid) then
      shift_x = xstart(1) - lbound(fx,1)
      shift_y = xstart(2) - lbound(fx,2)
      shift_z = xstart(3) - lbound(fx,3)
      do i = 1, nx_loc
        ig = i + shift_x
        xg(i) = xg_in(ig)
      enddo
      do j = 1, ny_loc
        jg = j + shift_y
        yg(j) = yg_in(jg)
      enddo
      do k = 1, nz_loc
        kg = k + shift_z
        zg(k) = zg_in(kg)
      enddo
      hx = xg(2) - xg(1)
      hz = zg(2) - zg(1)
      hy = -1._mytype
      jmin_g = 1
      jmax_g = size(yg_in)
    else
      hx = xlx / real(nx, mytype)
      if (ny > 1) then
        hy = yly / real(ny - 1, mytype)
      else
        hy = 1._mytype
      endif
      hz = zlz / real(nz, mytype)

      do i = 1, nx_loc
        xg(i) = real(i - 1, mytype) * hx
      enddo
      do j = 1, ny_loc
        yg(j) = real(j - 1, mytype) * hy
      enddo
      do k = 1, nz_loc
        zg(k) = real(k - 1, mytype) * hz
      enddo
      jmin_g = 1
      jmax_g = ny_loc
    endif

    fx = 0._mytype
    fy = 0._mytype
    fz = 0._mytype
    sumw_l = 0._mytype
    hy_loc_min = huge(1._mytype)
    hy_loc_max = 0._mytype

    do l = 1, fiber_nl
      do k = 1, nz_loc
        rz = minimum_image(zg(k) - fiber_x(3,l), zlz)
        if (abs(rz) > 2._mytype * hz) cycle
        do j = 1, ny_loc
          if (use_input_grid) then
            jg = j + shift_y
            hy_loc = local_dy(yg_in, jg, jmin_g, jmax_g)
          else
            hy_loc = hy
          endif
          hy_loc_min = min(hy_loc_min, hy_loc)
          hy_loc_max = max(hy_loc_max, hy_loc)

          ry = yg(j) - fiber_x(2,l)
          if (abs(ry) > 2._mytype * hy_loc) cycle

          do i = 1, nx_loc
            rx = minimum_image(xg(i) - fiber_x(1,l), xlx)
            if (abs(rx) > 2._mytype * hx) cycle
            weight = delta_kernel_3d(rx, ry, rz, hx, hy_loc, hz)
            fx(i,j,k) = fx(i,j,k) + fiber_force_lag(1,l) * fiber_quad_w(l) * weight
            fy(i,j,k) = fy(i,j,k) + fiber_force_lag(2,l) * fiber_quad_w(l) * weight
            fz(i,j,k) = fz(i,j,k) + fiber_force_lag(3,l) * fiber_quad_w(l) * weight

            dVloc = hx * hy_loc * hz
            sumw_l(l) = sumw_l(l) + weight * dVloc
          enddo
        enddo
      enddo
    enddo

    eul_total = 0._mytype
    do k = 1, nz_loc
      do j = 1, ny_loc
        if (use_input_grid) then
          jg = j + shift_y
          hy_loc = local_dy(yg_in, jg, jmin_g, jmax_g)
        else
          hy_loc = hy
        endif
        dVloc = hx * hy_loc * hz
        do i = 1, nx_loc
          eul_total(1) = eul_total(1) + fx(i,j,k) * dVloc
          eul_total(2) = eul_total(2) + fy(i,j,k) * dVloc
          eul_total(3) = eul_total(3) + fz(i,j,k) * dVloc
        enddo
      enddo
    enddo

    call MPI_ALLREDUCE(sumw_l, sumw_l_global, fiber_nl, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, hy_loc_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, hy_loc_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    sumw_min = minval(sumw_l_global)
    sumw_max = maxval(sumw_l_global)

    if (present(hy_loc_min_out)) hy_loc_min_out = hy_loc_min
    if (present(hy_loc_max_out)) hy_loc_max_out = hy_loc_max

    deallocate(xg, yg, zg, sumw_l, sumw_l_global)

  end subroutine spread_lagrangian_force_to_euler

  subroutine run_fiber_spread_conservation_test(lag_total, eul_total, abs_error, rel_error, sumw_min, sumw_max)

    real(mytype), intent(out), dimension(3) :: lag_total, eul_total, abs_error, rel_error
    real(mytype), intent(out) :: sumw_min, sumw_max

    real(mytype), allocatable, dimension(:,:,:) :: fx, fy, fz
    integer :: i, l

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

    allocate(fx(nx,ny,nz), fy(nx,ny,nz), fz(nx,ny,nz))

    call spread_lagrangian_force_to_euler(fiber_test_force, fx, fy, fz, eul_total, sumw_min, sumw_max)

    lag_total(1) = sum(fiber_test_force(1,:) * fiber_quad_w(:))
    lag_total(2) = sum(fiber_test_force(2,:) * fiber_quad_w(:))
    lag_total(3) = sum(fiber_test_force(3,:) * fiber_quad_w(:))

    abs_error = abs(eul_total - lag_total)

    rel_error = 0._mytype
    do i = 1, 3
      if (abs(lag_total(i)) > 0._mytype) then
        rel_error(i) = abs_error(i) / abs(lag_total(i))
      endif
    enddo

    if (nrank == 0) then
      write(*,'(A,I6)') 'Fiber spread conservation test case: ', spread_test_case
      write(*,'(A,3ES12.4)') 'Lagrangian total force          : ', lag_total
      write(*,'(A,3ES12.4)') 'Euler integrated spread force   : ', eul_total
      write(*,'(A,3ES12.4)') 'Absolute conservation error     : ', abs_error
      write(*,'(A,3ES12.4)') 'Relative conservation error     : ', rel_error
      write(*,'(A,2ES12.4)') 'Spread kernel sumw min/max      : ', sumw_min, sumw_max
    endif

    deallocate(fx, fy, fz)

  end subroutine run_fiber_spread_conservation_test

end module fiber_spread
