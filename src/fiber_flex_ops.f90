!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_flex_ops

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use fiber_types, only : fiber_active, fiber_flexible_active, fiber_flex_initialized, fiber_nl, fiber_length, &
       fiber_ds, fiber_s_ref, fiber_x, fiber_xs, fiber_xss, fiber_xssss, fiber_flex_operator_case
  use fiber_io, only : write_fiber_flex_operator_points, write_fiber_flex_operator_summary

  implicit none

contains

  subroutine compute_flexible_spatial_operators(x_in, xs, xss, xssss, bc_moment_left, bc_moment_right, &
       bc_shear_left, bc_shear_right)

    ! Step 3.2 only: spatial operators + free-end boundary verification layer.
    ! No flexible dynamics/tension/inextensibility/IBM coupling is implemented here.
    !
    ! Free-end closure uses ghost-point elimination from:
    !   X''(0)=0, X'''(0)=0, X''(L)=0, X'''(L)=0
    ! with central formulas at the boundary to eliminate ghost values.

    real(mytype), intent(in) :: x_in(3, fiber_nl)
    real(mytype), intent(out) :: xs(3, fiber_nl), xss(3, fiber_nl), xssss(3, fiber_nl)
    real(mytype), intent(out) :: bc_moment_left(3), bc_moment_right(3)
    real(mytype), intent(out) :: bc_shear_left(3), bc_shear_right(3)

    real(mytype) :: xext(3, fiber_nl + 4)
    real(mytype) :: ds, ds2, ds3, ds4
    integer :: i

    ds = fiber_ds
    ds2 = ds * ds
    ds3 = ds2 * ds
    ds4 = ds2 * ds2

    xext = 0._mytype
    do i = 1, fiber_nl
      xext(:, i + 2) = x_in(:, i)
    enddo

    xext(:, 2) = 2._mytype * xext(:, 3) - xext(:, 4)
    xext(:, 1) = 4._mytype * xext(:, 3) - 4._mytype * xext(:, 4) + xext(:, 5)
    xext(:, fiber_nl + 3) = 2._mytype * xext(:, fiber_nl + 2) - xext(:, fiber_nl + 1)
    xext(:, fiber_nl + 4) = 4._mytype * xext(:, fiber_nl + 2) - 4._mytype * xext(:, fiber_nl + 1) + xext(:, fiber_nl)

    do i = 1, fiber_nl
      xs(:, i) = (xext(:, i + 3) - xext(:, i + 1)) / (2._mytype * ds)
      xss(:, i) = (xext(:, i + 3) - 2._mytype * xext(:, i + 2) + xext(:, i + 1)) / ds2
      xssss(:, i) = (xext(:, i) - 4._mytype * xext(:, i + 1) + 6._mytype * xext(:, i + 2) - &
           4._mytype * xext(:, i + 3) + xext(:, i + 4)) / ds4
    enddo

    bc_moment_left = xss(:, 1)
    bc_moment_right = xss(:, fiber_nl)
    bc_shear_left = (xext(:, 1) - 2._mytype * xext(:, 2) + 2._mytype * xext(:, 4) - xext(:, 5)) / (2._mytype * ds3)
    bc_shear_right = (xext(:, fiber_nl) - 2._mytype * xext(:, fiber_nl + 1) + &
         2._mytype * xext(:, fiber_nl + 3) - xext(:, fiber_nl + 4)) / (2._mytype * ds3)

  end subroutine compute_flexible_spatial_operators

  subroutine run_flexible_operator_test()

    real(mytype), parameter :: a1 = 1._mytype, a2 = -0.5_mytype, a3 = 0.25_mytype
    real(mytype), allocatable :: x_exact(:,:), xs_exact(:,:), xss_exact(:,:), xssss_exact(:,:)
    real(mytype) :: bc_moment_left(3), bc_moment_right(3), bc_shear_left(3), bc_shear_right(3)
    real(mytype) :: err_xs_max, err_xss_max, err_xssss_max
    real(mytype) :: bc_moment_left_max, bc_moment_right_max, bc_shear_left_max, bc_shear_right_max
    real(mytype) :: smin, s, q, invL, invL2, invL4
    real(mytype) :: g, g1, g2, g4
    integer :: l

    if (.not.fiber_active) then
      if (nrank == 0) write(*,*) 'Error: fiber_flex_operator_test_active requires fiber_active = true.'
      stop
    endif
    if (.not.fiber_flexible_active) then
      if (nrank == 0) write(*,*) 'Error: fiber_flex_operator_test_active requires fiber_flexible_active = true.'
      stop
    endif
    if (.not.fiber_flex_initialized) then
      if (nrank == 0) write(*,*) 'Error: run_flexible_operator_test requires flexible state initialization first.'
      stop
    endif
    if (fiber_nl < 5) then
      if (nrank == 0) write(*,*) 'Error: run_flexible_operator_test requires fiber_nl >= 5.'
      stop
    endif

    if (allocated(x_exact)) deallocate(x_exact)
    if (allocated(xs_exact)) deallocate(xs_exact)
    if (allocated(xss_exact)) deallocate(xss_exact)
    if (allocated(xssss_exact)) deallocate(xssss_exact)
    allocate(x_exact(3, fiber_nl), xs_exact(3, fiber_nl), xss_exact(3, fiber_nl), xssss_exact(3, fiber_nl))

    if (allocated(fiber_xs)) deallocate(fiber_xs)
    if (allocated(fiber_xss)) deallocate(fiber_xss)
    if (allocated(fiber_xssss)) deallocate(fiber_xssss)
    allocate(fiber_xs(3, fiber_nl), fiber_xss(3, fiber_nl), fiber_xssss(3, fiber_nl))

    smin = -0.5_mytype * fiber_length
    invL = 1._mytype / fiber_length
    invL2 = invL * invL
    invL4 = invL2 * invL2

    do l = 1, fiber_nl
      s = fiber_s_ref(l)
      q = (s - smin) * invL

      g = q**4 - 4._mytype*q**5 + 6._mytype*q**6 - 4._mytype*q**7 + q**8
      g1 = 4._mytype*q**3 - 20._mytype*q**4 + 36._mytype*q**5 - 28._mytype*q**6 + 8._mytype*q**7
      g2 = 12._mytype*q**2 - 80._mytype*q**3 + 180._mytype*q**4 - 168._mytype*q**5 + 56._mytype*q**6
      g4 = 24._mytype - 480._mytype*q + 2160._mytype*q**2 - 3360._mytype*q**3 + 1680._mytype*q**4

      x_exact(1,l) = a1 * g
      x_exact(2,l) = a2 * g
      x_exact(3,l) = a3 * g

      xs_exact(1,l) = a1 * g1 * invL
      xs_exact(2,l) = a2 * g1 * invL
      xs_exact(3,l) = a3 * g1 * invL

      xss_exact(1,l) = a1 * g2 * invL2
      xss_exact(2,l) = a2 * g2 * invL2
      xss_exact(3,l) = a3 * g2 * invL2

      xssss_exact(1,l) = a1 * g4 * invL4
      xssss_exact(2,l) = a2 * g4 * invL4
      xssss_exact(3,l) = a3 * g4 * invL4
    enddo

    fiber_x = x_exact

    call compute_flexible_spatial_operators(fiber_x, fiber_xs, fiber_xss, fiber_xssss, bc_moment_left, bc_moment_right, &
         bc_shear_left, bc_shear_right)

    err_xs_max = maxval(abs(fiber_xs - xs_exact))
    err_xss_max = maxval(abs(fiber_xss - xss_exact))
    err_xssss_max = maxval(abs(fiber_xssss - xssss_exact))
    bc_moment_left_max = maxval(abs(bc_moment_left))
    bc_moment_right_max = maxval(abs(bc_moment_right))
    bc_shear_left_max = maxval(abs(bc_shear_left))
    bc_shear_right_max = maxval(abs(bc_shear_right))

    call write_fiber_flex_operator_points(x_exact, fiber_xs, xs_exact, fiber_xss, xss_exact, fiber_xssss, xssss_exact)
    call write_fiber_flex_operator_summary(fiber_flex_operator_case, err_xs_max, err_xss_max, err_xssss_max, &
         bc_moment_left_max, bc_moment_right_max, bc_shear_left_max, bc_shear_right_max)

    deallocate(x_exact, xs_exact, xss_exact, xssss_exact)

  end subroutine run_flexible_operator_test

end module fiber_flex_ops
