!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_flex_constraint

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use fiber_types, only : fiber_active, fiber_flexible_active, fiber_flex_initialized, fiber_nl, fiber_length, fiber_ds, &
       fiber_s_ref, fiber_x, fiber_flex_constraint_test_active, fiber_flex_constraint_case, fiber_flex_constraint_nsteps, &
       fiber_flex_constraint_output_interval, fiber_flex_constraint_dt, fiber_flex_constraint_force_amp
  use fiber_types, only : fiber_bending_gamma
  use fiber_flex_ops, only : apply_free_end_bending_operator_scalar
  use fiber_io, only : write_fiber_flex_constraint_initial_points, write_fiber_flex_constraint_final_points, &
       write_fiber_flex_constraint_tension_last, write_fiber_flex_constraint_series, write_fiber_flex_constraint_summary

  implicit none

contains

  subroutine compute_centerline_velocity(xn, xnm1, dt_c, un)
    real(mytype), intent(in) :: xn(3,fiber_nl), xnm1(3,fiber_nl), dt_c
    real(mytype), intent(out) :: un(3,fiber_nl)
    un = (xn - xnm1) / dt_c
  end subroutine compute_centerline_velocity

  subroutine compute_constraint_error(x_in, max_seg_err, max_inext_err)
    real(mytype), intent(in) :: x_in(3,fiber_nl)
    real(mytype), intent(out) :: max_seg_err, max_inext_err
    real(mytype) :: dsvec(3), seglen, g
    integer :: i
    max_seg_err = 0._mytype
    max_inext_err = 0._mytype
    do i = 1, fiber_nl - 1
      dsvec = (x_in(:,i+1) - x_in(:,i)) / fiber_ds
      seglen = sqrt(sum((x_in(:,i+1)-x_in(:,i))**2))
      g = dot_product(dsvec, dsvec) - 1._mytype
      max_seg_err = max(max_seg_err, abs(seglen - fiber_ds))
      max_inext_err = max(max_inext_err, abs(g))
    enddo
  end subroutine compute_constraint_error

  subroutine apply_tension_force_from_halfgrid(t_half, x_ref, ftension)
    real(mytype), intent(in) :: t_half(fiber_nl-1), x_ref(3,fiber_nl)
    real(mytype), intent(out) :: ftension(3,fiber_nl)
    real(mytype) :: xs_half(3,fiber_nl-1), flux(3,fiber_nl-1)
    integer :: i

    do i = 1, fiber_nl-1
      xs_half(:,i) = (x_ref(:,i+1) - x_ref(:,i)) / fiber_ds
      flux(:,i) = t_half(i) * xs_half(:,i)
    enddo

    ftension = 0._mytype
    ftension(:,1) = flux(:,1) / fiber_ds
    do i = 2, fiber_nl-1
      ftension(:,i) = (flux(:,i) - flux(:,i-1)) / fiber_ds
    enddo
    ftension(:,fiber_nl) = -flux(:,fiber_nl-1) / fiber_ds
  end subroutine apply_tension_force_from_halfgrid

  subroutine apply_tension_poisson_operator(t_half_unknown, x_ref, lhs_t)
    real(mytype), intent(in) :: t_half_unknown(fiber_nl-1), x_ref(3,fiber_nl)
    real(mytype), intent(out) :: lhs_t(fiber_nl-1)
    real(mytype) :: ft(3,fiber_nl), xs_half(3,fiber_nl-1), dft(3,fiber_nl-1)
    integer :: i

    call apply_tension_force_from_halfgrid(t_half_unknown, x_ref, ft)
    do i = 1, fiber_nl-1
      xs_half(:,i) = (x_ref(:,i+1) - x_ref(:,i)) / fiber_ds
      dft(:,i) = (ft(:,i+1) - ft(:,i)) / fiber_ds
      lhs_t(i) = 2._mytype * dot_product(xs_half(:,i), dft(:,i))
    enddo
  end subroutine apply_tension_poisson_operator

  subroutine solve_dense_small(ain, b, x, n)
    integer, intent(in) :: n
    real(mytype), intent(in) :: ain(n,n), b(n)
    real(mytype), intent(out) :: x(n)
    real(mytype) :: a(n,n), rhs(n), rowtmp(n), piv, fac, tmp
    integer :: i, k, p
    a = ain; rhs = b; x = 0._mytype
    do k = 1, n
      p = k; piv = abs(a(k,k))
      do i = k+1, n
        if (abs(a(i,k)) > piv) then
          piv = abs(a(i,k)); p = i
        endif
      enddo
      if (piv <= 1.0e-14_mytype) then
        if (nrank == 0) write(*,*) 'Error: singular matrix in solve_dense_small.'
        stop
      endif
      if (p /= k) then
        rowtmp = a(k,:); a(k,:) = a(p,:); a(p,:) = rowtmp
        tmp = rhs(k); rhs(k) = rhs(p); rhs(p) = tmp
      endif
      do i = k+1, n
        fac = a(i,k) / a(k,k)
        a(i,k:n) = a(i,k:n) - fac * a(k,k:n)
        rhs(i) = rhs(i) - fac * rhs(k)
      enddo
    enddo
    do i = n, 1, -1
      if (i < n) then
        x(i) = (rhs(i) - sum(a(i,i+1:n)*x(i+1:n))) / a(i,i)
      else
        x(i) = rhs(i) / a(i,i)
      endif
    enddo
  end subroutine solve_dense_small

  subroutine solve_constraint_tension(xn, xnm1, x_ref, fb_ref, fext_ref, dt_c, t_half)
    real(mytype), intent(in) :: xn(3,fiber_nl), xnm1(3,fiber_nl), x_ref(3,fiber_nl)
    real(mytype), intent(in) :: fb_ref(3,fiber_nl), fext_ref(3,fiber_nl), dt_c
    real(mytype), intent(out) :: t_half(fiber_nl-1)
    real(mytype) :: amat(fiber_nl-1,fiber_nl-1), rhs(fiber_nl-1), evec(fiber_nl-1), col(fiber_nl-1)
    real(mytype) :: un(3,fiber_nl), xs_half(3,fiber_nl-1), dun(3,fiber_nl-1), drift, velterm, forceterm
    real(mytype) :: kappa_drift
    integer :: i, j

    call compute_centerline_velocity(xn, xnm1, dt_c, un)
    do i = 1, fiber_nl-1
      xs_half(:,i) = (x_ref(:,i+1) - x_ref(:,i)) / fiber_ds
      dun(:,i) = (un(:,i+1) - un(:,i)) / fiber_ds
    enddo

    do j = 1, fiber_nl-1
      evec = 0._mytype; evec(j) = 1._mytype
      call apply_tension_poisson_operator(evec, x_ref, col)
      amat(:,j) = col
    enddo

    kappa_drift = 2._mytype / (dt_c * dt_c)
    do i = 1, fiber_nl-1
      drift = dot_product(xs_half(:,i), xs_half(:,i)) - 1._mytype
      velterm = dot_product(dun(:,i), dun(:,i))
      forceterm = 2._mytype * dot_product(xs_half(:,i), (fb_ref(:,i+1)+fext_ref(:,i+1)-fb_ref(:,i)-fext_ref(:,i))/fiber_ds)
      rhs(i) = -kappa_drift * drift - velterm - forceterm
    enddo

    call solve_dense_small(amat, rhs, t_half, fiber_nl-1)
  end subroutine solve_constraint_tension

  subroutine advance_constrained_structure_step(xnm1, xn, xnp1, t_half, fext_n, dt_c)
    real(mytype), intent(in) :: xnm1(3,fiber_nl), xn(3,fiber_nl), fext_n(3,fiber_nl), dt_c
    real(mytype), intent(out) :: xnp1(3,fiber_nl), t_half(fiber_nl-1)
    real(mytype) :: x_ref(3,fiber_nl), fb_ref(3,fiber_nl), ftension(3,fiber_nl)
    integer :: c

    x_ref = 2._mytype * xn - xnm1
    do c = 1, 3
      call apply_free_end_bending_operator_scalar(x_ref(c,:), fb_ref(c,:), fiber_bending_gamma)
    enddo
    call solve_constraint_tension(xn, xnm1, x_ref, fb_ref, fext_n, dt_c, t_half)
    call apply_tension_force_from_halfgrid(t_half, x_ref, ftension)

    xnp1 = 2._mytype * xn - xnm1 + dt_c * dt_c * (ftension + fb_ref + fext_n)
  end subroutine advance_constrained_structure_step

  subroutine initialize_constraint_test_state(case_id, xnm1, xn)
    integer, intent(in) :: case_id
    real(mytype), intent(out) :: xnm1(3,fiber_nl), xn(3,fiber_nl)
    real(mytype) :: vel0(3), amp, smin, q, g
    integer :: l

    select case (case_id)
    case (1)
      xn = fiber_x
      xnm1 = xn
    case (2)
      vel0 = (/0.05_mytype, 0._mytype, 0._mytype/)
      xn = fiber_x
      xnm1 = xn - fiber_flex_constraint_dt * spread(vel0, 2, fiber_nl)
    case (3)
      smin = -0.5_mytype * fiber_length
      amp = 0.05_mytype * fiber_length
      do l = 1, fiber_nl
        q = (fiber_s_ref(l) - smin) / fiber_length
        g = q**4 * (1._mytype - q)**4
        xn(1,l) = fiber_s_ref(l)
        xn(2,l) = amp * g
        xn(3,l) = 0._mytype
      enddo
      xnm1 = xn
    case default
      if (nrank == 0) write(*,*) 'Error: fiber_flex_constraint_case must be 1,2,3.'
      stop
    end select
  end subroutine initialize_constraint_test_state

  subroutine run_flexible_constraint_test()
    real(mytype), allocatable :: xnm1(:,:), xn(:,:), xnp1(:,:), xinit(:,:), fext(:,:), t_half(:)
    real(mytype) :: max_seg_err, max_inext_err, max_abs_tension
    real(mytype) :: max_seg_global, max_inext_global, max_t_global, tnow
    real(mytype) :: case_metric
    integer :: step, outint, i

    if (.not.fiber_active .or. .not.fiber_flexible_active .or. .not.fiber_flex_initialized .or. &
         .not.fiber_flex_constraint_test_active) then
      if (nrank == 0) write(*,*) 'Error: invalid state for flexible constraint test.'
      stop
    endif

    allocate(xnm1(3,fiber_nl), xn(3,fiber_nl), xnp1(3,fiber_nl), xinit(3,fiber_nl), fext(3,fiber_nl), t_half(fiber_nl-1))
    call initialize_constraint_test_state(fiber_flex_constraint_case, xnm1, xn)
    xinit = xn
    call write_fiber_flex_constraint_initial_points(xinit)

    max_seg_global = 0._mytype
    max_inext_global = 0._mytype
    max_t_global = 0._mytype
    outint = max(1, fiber_flex_constraint_output_interval)
    call write_fiber_flex_constraint_series(0, 0._mytype, 0._mytype, 0._mytype, 0._mytype, .true.)

    do step = 1, fiber_flex_constraint_nsteps
      fext = 0._mytype
      if (fiber_flex_constraint_case == 3) then
        do i = 1, fiber_nl
          fext(2,i) = fiber_flex_constraint_force_amp * sin(2._mytype * acos(-1._mytype) * real(step,mytype) / &
               real(max(1,fiber_flex_constraint_nsteps),mytype))
        enddo
      endif

      call advance_constrained_structure_step(xnm1, xn, xnp1, t_half, fext, fiber_flex_constraint_dt)
      call compute_constraint_error(xnp1, max_seg_err, max_inext_err)
      max_abs_tension = maxval(abs(t_half))
      max_seg_global = max(max_seg_global, max_seg_err)
      max_inext_global = max(max_inext_global, max_inext_err)
      max_t_global = max(max_t_global, max_abs_tension)

      tnow = real(step,mytype) * fiber_flex_constraint_dt
      if (mod(step,outint) == 0 .or. step == fiber_flex_constraint_nsteps) then
        call write_fiber_flex_constraint_series(step, tnow, max_seg_err, max_inext_err, max_abs_tension, .false.)
      endif

      xnm1 = xn
      xn = xnp1
    enddo

    call write_fiber_flex_constraint_final_points(xn)
    call write_fiber_flex_constraint_tension_last(t_half)

    case_metric = 0._mytype
    select case (fiber_flex_constraint_case)
    case (1)
      case_metric = maxval(abs(xn - xinit))
    case (2)
      case_metric = maxval(abs((xn(:,2:fiber_nl)-xn(:,1:fiber_nl-1)) - (xinit(:,2:fiber_nl)-xinit(:,1:fiber_nl-1))))
    case (3)
      case_metric = sqrt(sum((xn - xinit)**2) / real(3*fiber_nl,mytype))
    end select

    call write_fiber_flex_constraint_summary(fiber_flex_constraint_case, fiber_flex_constraint_dt, &
         fiber_flex_constraint_nsteps, max_seg_global, max_inext_global, max_t_global, case_metric)

    deallocate(xnm1, xn, xnp1, xinit, fext, t_half)
  end subroutine run_flexible_constraint_test

end module fiber_flex_constraint
