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

  subroutine estimate_boundary_d2_independent(x_in, d2_left_indep, d2_right_indep)

    ! Independent boundary second-derivative estimator (D2) from physical nodes only.
    ! No ghost values, closure constraints, or bending operators are used here.
    ! This uses a higher-order one-sided stencil to improve stability/credibility
    ! of boundary-error decay under grid refinement.

    real(mytype), intent(in) :: x_in(3, fiber_nl)
    real(mytype), intent(out) :: d2_left_indep(3), d2_right_indep(3)
    real(mytype) :: ds2
    integer :: i

    if (fiber_nl < 6) then
      if (nrank == 0) write(*,*) 'Error: estimate_boundary_d2_independent requires fiber_nl >= 6.'
      stop
    endif

    ds2 = fiber_ds * fiber_ds
    do i = 1, 3
      d2_left_indep(i) = (45._mytype * x_in(i,1) - 154._mytype * x_in(i,2) + 214._mytype * x_in(i,3) - &
           156._mytype * x_in(i,4) + 61._mytype * x_in(i,5) - 10._mytype * x_in(i,6)) / (12._mytype * ds2)
      d2_right_indep(i) = (45._mytype * x_in(i,fiber_nl) - 154._mytype * x_in(i,fiber_nl-1) + &
           214._mytype * x_in(i,fiber_nl-2) - 156._mytype * x_in(i,fiber_nl-3) + &
           61._mytype * x_in(i,fiber_nl-4) - 10._mytype * x_in(i,fiber_nl-5)) / (12._mytype * ds2)
    enddo

  end subroutine estimate_boundary_d2_independent

  subroutine estimate_boundary_d3_independent(x_in, d3_left_indep, d3_right_indep)

    ! Independent boundary shear estimate from physical nodes only.
    ! Uses second-order one-sided stencils and does not use ghost values.

    real(mytype), intent(in) :: x_in(3, fiber_nl)
    real(mytype), intent(out) :: d3_left_indep(3), d3_right_indep(3)
    real(mytype) :: ds3
    integer :: i

    ds3 = fiber_ds * fiber_ds * fiber_ds
    do i = 1, 3
      d3_left_indep(i) = (-5._mytype * x_in(i,1) + 18._mytype * x_in(i,2) - 24._mytype * x_in(i,3) + &
           14._mytype * x_in(i,4) - 3._mytype * x_in(i,5)) / (2._mytype * ds3)
      d3_right_indep(i) = (5._mytype * x_in(i,fiber_nl) - 18._mytype * x_in(i,fiber_nl-1) + &
           24._mytype * x_in(i,fiber_nl-2) - 14._mytype * x_in(i,fiber_nl-3) + 3._mytype * x_in(i,fiber_nl-4)) / &
           (2._mytype * ds3)
    enddo

  end subroutine estimate_boundary_d3_independent

  subroutine construct_free_end_ghost_scalar(u, uext)

    real(mytype), intent(in) :: u(fiber_nl)
    real(mytype), intent(out) :: uext(fiber_nl + 4)
    integer :: i

    uext = 0._mytype
    do i = 1, fiber_nl
      uext(i + 2) = u(i)
    enddo

    uext(2) = 2._mytype * uext(3) - uext(4)
    uext(1) = 4._mytype * uext(3) - 4._mytype * uext(4) + uext(5)
    uext(fiber_nl + 3) = 2._mytype * uext(fiber_nl + 2) - uext(fiber_nl + 1)
    uext(fiber_nl + 4) = 4._mytype * uext(fiber_nl + 2) - 4._mytype * uext(fiber_nl + 1) + uext(fiber_nl)

  end subroutine construct_free_end_ghost_scalar

  subroutine apply_free_end_d2_scalar(u, d2u)

    real(mytype), intent(in) :: u(fiber_nl)
    real(mytype), intent(out) :: d2u(fiber_nl)
    real(mytype) :: uext(fiber_nl + 4), ds2
    integer :: i

    call construct_free_end_ghost_scalar(u, uext)
    ds2 = fiber_ds * fiber_ds

    do i = 1, fiber_nl
      d2u(i) = (uext(i + 3) - 2._mytype * uext(i + 2) + uext(i + 1)) / ds2
    enddo

  end subroutine apply_free_end_d2_scalar

  subroutine apply_free_end_d4_scalar(u, d4u)

    real(mytype), intent(in) :: u(fiber_nl)
    real(mytype), intent(out) :: d4u(fiber_nl)
    real(mytype) :: uext(fiber_nl + 4), ds4
    integer :: i

    call construct_free_end_ghost_scalar(u, uext)
    ds4 = fiber_ds**4

    do i = 1, fiber_nl
      d4u(i) = (uext(i) - 4._mytype * uext(i + 1) + 6._mytype * uext(i + 2) - 4._mytype * uext(i + 3) + uext(i + 4)) / ds4
    enddo

  end subroutine apply_free_end_d4_scalar

  subroutine apply_free_end_bending_operator_scalar(u, bu, gamma_b)

    ! Step 3.3 bending-only operator semantics:
    !   B(u) = -d_ss( gamma * d_ss u )
    ! For constant gamma_b, this is equivalent to:
    !   B(u) = -gamma_b * d4u

    real(mytype), intent(in) :: u(fiber_nl), gamma_b
    real(mytype), intent(out) :: bu(fiber_nl)
    real(mytype) :: d4u(fiber_nl)

    call apply_free_end_d4_scalar(u, d4u)
    bu = -gamma_b * d4u

  end subroutine apply_free_end_bending_operator_scalar

  subroutine compute_bending_force_free_end(x_in, gamma_b, fb_out)

    ! Structural bending force operator in vector form:
    !   F_b = - d_ss( gamma * d_ss X )
    ! For constant gamma_b in this Step 3.2 verification layer:
    !   F_b = -gamma_b * X_ssss
    ! This routine computes the force via the bending operator, component by component.

    real(mytype), intent(in) :: x_in(3, fiber_nl), gamma_b
    real(mytype), intent(out) :: fb_out(3, fiber_nl)
    integer :: i

    do i = 1, 3
      call apply_free_end_bending_operator_scalar(x_in(i,:), fb_out(i,:), gamma_b)
    enddo

  end subroutine compute_bending_force_free_end

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

    real(mytype) :: xext(3, fiber_nl + 4), uext(fiber_nl + 4)
    real(mytype) :: ds, ds3
    real(mytype) :: d2tmp(fiber_nl), d4tmp(fiber_nl)
    integer :: i

    ds = fiber_ds
    ds3 = ds * ds * ds

    do i = 1, 3
      call construct_free_end_ghost_scalar(x_in(i,:), uext)
      xext(i,:) = uext
      call apply_free_end_d2_scalar(x_in(i,:), d2tmp)
      call apply_free_end_d4_scalar(x_in(i,:), d4tmp)
      xss(i,:) = d2tmp
      xssss(i,:) = d4tmp
      xs(i,:) = (uext(4:fiber_nl+3) - uext(2:fiber_nl+1)) / (2._mytype * ds)
      bc_shear_left(i) = (uext(1) - 2._mytype * uext(2) + 2._mytype * uext(4) - uext(5)) / (2._mytype * ds3)
      bc_shear_right(i) = (uext(fiber_nl) - 2._mytype * uext(fiber_nl + 1) + 2._mytype * uext(fiber_nl + 3) - &
           uext(fiber_nl + 4)) / (2._mytype * ds3)
    enddo

    bc_moment_left = xss(:,1)
    bc_moment_right = xss(:,fiber_nl)

  end subroutine compute_flexible_spatial_operators

  subroutine run_flexible_operator_test()

    real(mytype), parameter :: a1 = 1._mytype, a2 = -0.5_mytype, a3 = 0.25_mytype
    real(mytype), parameter :: gamma_test = 1._mytype
    real(mytype), allocatable :: x_exact(:,:), xs_exact(:,:), xss_exact(:,:), xsss_exact(:,:), xssss_exact(:,:)
    real(mytype), allocatable :: fb_exact(:,:), fb_num(:,:)
    real(mytype) :: bc_moment_left(3), bc_moment_right(3), bc_shear_left(3), bc_shear_right(3)
    real(mytype) :: d2_left_indep(3), d2_right_indep(3), d3_left_indep(3), d3_right_indep(3)
    real(mytype) :: err_xs_max, err_xss_max, err_xssss_max, err_fb_max
    real(mytype) :: bc_moment_left_max, bc_moment_right_max, bc_shear_left_max, bc_shear_right_max
    real(mytype) :: bc_d2_left_indep_err_max, bc_d2_right_indep_err_max
    real(mytype) :: bc_d3_left_indep_err_max, bc_d3_right_indep_err_max
    real(mytype) :: smin, s, q, invL, invL2, invL4
    real(mytype) :: invL3, g, g1, g2, g3, g4
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
    if (allocated(xsss_exact)) deallocate(xsss_exact)
    if (allocated(xssss_exact)) deallocate(xssss_exact)
    if (allocated(fb_exact)) deallocate(fb_exact)
    if (allocated(fb_num)) deallocate(fb_num)
    allocate(x_exact(3, fiber_nl), xs_exact(3, fiber_nl), xss_exact(3, fiber_nl), xsss_exact(3, fiber_nl), &
         xssss_exact(3, fiber_nl))
    allocate(fb_exact(3, fiber_nl), fb_num(3, fiber_nl))

    if (allocated(fiber_xs)) deallocate(fiber_xs)
    if (allocated(fiber_xss)) deallocate(fiber_xss)
    if (allocated(fiber_xssss)) deallocate(fiber_xssss)
    allocate(fiber_xs(3, fiber_nl), fiber_xss(3, fiber_nl), fiber_xssss(3, fiber_nl))

    smin = -0.5_mytype * fiber_length
    invL = 1._mytype / fiber_length
    invL2 = invL * invL
    invL3 = invL2 * invL
    invL4 = invL2 * invL2

    do l = 1, fiber_nl
      s = fiber_s_ref(l)
      q = (s - smin) * invL

      select case (fiber_flex_operator_case)
      case (1)
        g = q**4 - 4._mytype*q**5 + 6._mytype*q**6 - 4._mytype*q**7 + q**8
        g1 = 4._mytype*q**3 - 20._mytype*q**4 + 36._mytype*q**5 - 28._mytype*q**6 + 8._mytype*q**7
        g2 = 12._mytype*q**2 - 80._mytype*q**3 + 180._mytype*q**4 - 168._mytype*q**5 + 56._mytype*q**6
        g3 = 24._mytype*q - 240._mytype*q**2 + 720._mytype*q**3 - 840._mytype*q**4 + 336._mytype*q**5
        g4 = 24._mytype - 480._mytype*q + 2160._mytype*q**2 - 3360._mytype*q**3 + 1680._mytype*q**4
      case (2)
        ! Distinct free-end manufactured basis for case 2:
        ! g(q)=q^4(1-q)^4(2q-1), which keeps g''=g'''=0 at both ends but
        ! changes xss/xssss/Fb exact distributions versus case 1.
        g = -q**4 + 6._mytype*q**5 - 14._mytype*q**6 + 16._mytype*q**7 - 9._mytype*q**8 + 2._mytype*q**9
        g1 = -4._mytype*q**3 + 30._mytype*q**4 - 84._mytype*q**5 + 112._mytype*q**6 - 72._mytype*q**7 + 18._mytype*q**8
        g2 = -12._mytype*q**2 + 120._mytype*q**3 - 420._mytype*q**4 + 672._mytype*q**5 - 504._mytype*q**6 + 144._mytype*q**7
        g3 = -24._mytype*q + 360._mytype*q**2 - 1680._mytype*q**3 + 3360._mytype*q**4 - 3024._mytype*q**5 + 1008._mytype*q**6
        g4 = -24._mytype + 720._mytype*q - 5040._mytype*q**2 + 13440._mytype*q**3 - 15120._mytype*q**4 + 6048._mytype*q**5
      case default
        if (nrank == 0) write(*,*) 'Error: fiber_flex_operator_case must be 1 or 2.'
        stop
      end select

      x_exact(1,l) = a1 * g
      x_exact(2,l) = a2 * g
      x_exact(3,l) = a3 * g

      xs_exact(1,l) = a1 * g1 * invL
      xs_exact(2,l) = a2 * g1 * invL
      xs_exact(3,l) = a3 * g1 * invL

      xss_exact(1,l) = a1 * g2 * invL2
      xss_exact(2,l) = a2 * g2 * invL2
      xss_exact(3,l) = a3 * g2 * invL2

      xsss_exact(1,l) = a1 * g3 * invL3
      xsss_exact(2,l) = a2 * g3 * invL3
      xsss_exact(3,l) = a3 * g3 * invL3

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
    fb_exact = -gamma_test * xssss_exact
    call compute_bending_force_free_end(x_exact, gamma_test, fb_num)
    err_fb_max = maxval(abs(fb_num - fb_exact))
    bc_moment_left_max = maxval(abs(bc_moment_left))
    bc_moment_right_max = maxval(abs(bc_moment_right))
    bc_shear_left_max = maxval(abs(bc_shear_left))
    bc_shear_right_max = maxval(abs(bc_shear_right))
    call estimate_boundary_d2_independent(fiber_x, d2_left_indep, d2_right_indep)
    call estimate_boundary_d3_independent(fiber_x, d3_left_indep, d3_right_indep)
    bc_d2_left_indep_err_max = maxval(abs(d2_left_indep - xss_exact(:,1)))
    bc_d2_right_indep_err_max = maxval(abs(d2_right_indep - xss_exact(:,fiber_nl)))
    bc_d3_left_indep_err_max = maxval(abs(d3_left_indep - xsss_exact(:,1)))
    bc_d3_right_indep_err_max = maxval(abs(d3_right_indep - xsss_exact(:,fiber_nl)))

    call write_fiber_flex_operator_points(x_exact, fiber_xs, xs_exact, fiber_xss, xss_exact, fiber_xssss, xssss_exact)
    call write_fiber_flex_operator_summary(fiber_flex_operator_case, err_xs_max, err_xss_max, err_xssss_max, &
         gamma_test, err_fb_max, &
         bc_moment_left_max, bc_moment_right_max, bc_shear_left_max, bc_shear_right_max, &
         bc_d2_left_indep_err_max, bc_d2_right_indep_err_max, bc_d3_left_indep_err_max, bc_d3_right_indep_err_max)

    deallocate(x_exact, xs_exact, xss_exact, xsss_exact, xssss_exact, fb_exact, fb_num)

  end subroutine run_flexible_operator_test

end module fiber_flex_ops
