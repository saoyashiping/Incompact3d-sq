!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_flex_bending

  ! Step 3.3 scope definition:
  ! This module is a bending-only implicit intermediate layer used to
  ! verify the free-end bending operator and the numerically stiff bending term
  ! in isolation.
  !
  ! It does not include structural inertia, tension/inextensibility constraints,
  ! or fluid-structure coupling terms, and therefore is not the final
  ! structural time-integration solver.

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use fiber_types, only : fiber_active, fiber_flexible_active, fiber_flex_initialized, fiber_nl, fiber_length, fiber_ds, &
       fiber_s_ref, fiber_x, fiber_flex_bending_test_active, fiber_flex_bending_case, fiber_flex_bending_nsteps, &
       fiber_flex_bending_output_interval, fiber_flex_bending_dt, fiber_bending_gamma, &
       fiber_flex_bending_linear_solver, fiber_flex_bending_cg_tol, fiber_flex_bending_cg_maxit, &
       fiber_flex_bending_iter_tol_effective, fiber_flex_bending_iter_maxit_effective
  use fiber_flex_ops, only : apply_free_end_d2_scalar, apply_free_end_d4_scalar, apply_free_end_bending_operator_scalar
  use fiber_io, only : write_fiber_flex_bending_initial_points, write_fiber_flex_bending_final_points, &
       write_fiber_flex_bending_series, write_fiber_flex_bending_summary

  implicit none

contains

  subroutine assemble_free_end_bending_matrix(bmat, gamma_b)
    real(mytype), intent(out) :: bmat(fiber_nl, fiber_nl)
    real(mytype), intent(in) :: gamma_b
    real(mytype) :: evec(fiber_nl), bcol(fiber_nl)
    integer :: j

    bmat = 0._mytype
    do j = 1, fiber_nl
      evec = 0._mytype
      evec(j) = 1._mytype
      call apply_free_end_bending_operator_scalar(evec, bcol, gamma_b)
      bmat(:,j) = bcol
    enddo
  end subroutine assemble_free_end_bending_matrix

  subroutine build_bending_backward_euler_matrix(dt_b, gamma_b, amat)
    ! Step 3.3 bending-only backward Euler:
    ! continuous target: x_t = B(x), with B(x) = -gamma*d4(x)
    ! backward Euler: (x^{n+1} - x^n)/dt_b = B(x^{n+1})
    ! therefore: (I - dt_b * B) x^{n+1} = x^n
    ! since B = -gamma*d4, this is (I + dt_b * gamma * d4) x^{n+1} = x^n.
    real(mytype), intent(in) :: dt_b, gamma_b
    real(mytype), intent(out) :: amat(fiber_nl, fiber_nl)
    real(mytype) :: bmat(fiber_nl, fiber_nl)
    integer :: i

    call assemble_free_end_bending_matrix(bmat, gamma_b)
    amat = -dt_b * bmat
    do i = 1, fiber_nl
      amat(i,i) = amat(i,i) + 1._mytype
    enddo
  end subroutine build_bending_backward_euler_matrix

  subroutine solve_linear_system_dense(ain, b, x)
    real(mytype), intent(in) :: ain(fiber_nl, fiber_nl), b(fiber_nl)
    real(mytype), intent(out) :: x(fiber_nl)
    real(mytype) :: a(fiber_nl, fiber_nl), rhs(fiber_nl), rowtmp(fiber_nl), piv, fac, tmp
    integer :: i, j, k, p

    a = ain
    rhs = b
    x = 0._mytype

    do k = 1, fiber_nl
      p = k
      piv = abs(a(k,k))
      do i = k+1, fiber_nl
        if (abs(a(i,k)) > piv) then
          piv = abs(a(i,k))
          p = i
        endif
      enddo
      if (piv <= 1.0e-14_mytype) then
        if (nrank == 0) write(*,*) 'Error: singular matrix in solve_linear_system_dense.'
        stop
      endif
      if (p /= k) then
        rowtmp = a(k,:); a(k,:) = a(p,:); a(p,:) = rowtmp
        tmp = rhs(k); rhs(k) = rhs(p); rhs(p) = tmp
      endif
      do i = k+1, fiber_nl
        fac = a(i,k) / a(k,k)
        a(i,k:fiber_nl) = a(i,k:fiber_nl) - fac * a(k,k:fiber_nl)
        rhs(i) = rhs(i) - fac * rhs(k)
      enddo
    enddo

    do i = fiber_nl, 1, -1
      if (i < fiber_nl) then
        x(i) = (rhs(i) - sum(a(i,i+1:fiber_nl) * x(i+1:fiber_nl))) / a(i,i)
      else
        x(i) = rhs(i) / a(i,i)
      endif
    enddo
  end subroutine solve_linear_system_dense

  subroutine apply_bending_backward_euler_operator_scalar(v_in, av_out, dt_b, gamma_b)
    ! Operator action for Step 3.3 bending-only backward Euler:
    !   A v = (I + dt_b * gamma_b * D4) v
    ! This is the matrix-free primary path and reuses the validated D4 operator.
    real(mytype), intent(in) :: v_in(fiber_nl), dt_b, gamma_b
    real(mytype), intent(out) :: av_out(fiber_nl)
    real(mytype) :: d4v(fiber_nl)

    call apply_free_end_d4_scalar(v_in, d4v)
    av_out = v_in + dt_b * gamma_b * d4v
  end subroutine apply_bending_backward_euler_operator_scalar

  subroutine build_bending_operator_diagonal(diagA, dt_b, gamma_b)
    ! Build Jacobi diagonal directly from operator action to stay
    ! consistent with the validated matrix-free A*v implementation.
    real(mytype), intent(out) :: diagA(fiber_nl)
    real(mytype), intent(in) :: dt_b, gamma_b
    real(mytype) :: ej(fiber_nl), aej(fiber_nl)
    integer :: j

    do j = 1, fiber_nl
      ej = 0._mytype
      ej(j) = 1._mytype
      call apply_bending_backward_euler_operator_scalar(ej, aej, dt_b, gamma_b)
      diagA(j) = aej(j)
    enddo
  end subroutine build_bending_operator_diagonal

  subroutine apply_jacobi_preconditioner(r_in, z_out, diagA)
    ! Simple diagonal (Jacobi) preconditioner: z = M^{-1} r.
    real(mytype), intent(in) :: r_in(fiber_nl), diagA(fiber_nl)
    real(mytype), intent(out) :: z_out(fiber_nl)
    real(mytype), parameter :: diag_eps = 1.0e-30_mytype
    integer :: i

    do i = 1, fiber_nl
      if (abs(diagA(i)) <= diag_eps) then
        if (nrank == 0) write(*,*) 'Error: Jacobi preconditioner has near-zero diagonal in Step 3.3 at i=', i
        stop
      endif
      z_out(i) = r_in(i) / diagA(i)
    enddo
  end subroutine apply_jacobi_preconditioner

  subroutine solve_bending_scalar_bicgstab(rhs, x, dt_b, gamma_b, tol, maxit, it_used, final_residual, converged, &
       breakdown_detected, restart_count, failure_mode)
    ! Primary operator-based linear solver for Step 3.3.
    ! Uses a right-preconditioned BiCGSTAB primary path with simple
    ! Jacobi preconditioning:
    !   A M^{-1} y = b,   x = M^{-1} y
    ! Residuals/shadow residuals remain in the physical residual space.
    real(mytype), intent(in) :: rhs(fiber_nl), dt_b, gamma_b, tol
    integer, intent(in) :: maxit
    real(mytype), intent(inout) :: x(fiber_nl)
    real(mytype), intent(out) :: final_residual
    integer, intent(out) :: it_used, restart_count, failure_mode
    logical, intent(out) :: converged, breakdown_detected
    real(mytype) :: r(fiber_nl), rhat(fiber_nl), p(fiber_nl), v(fiber_nl), s(fiber_nl), t(fiber_nl)
    real(mytype) :: phat(fiber_nl), shat(fiber_nl), diagA(fiber_nl)
    real(mytype) :: rho, rho_old, alpha, omega, beta, rhat_v
    real(mytype) :: rhs_norm, target, tt, eps_break, old_residual
    real(mytype) :: ax(fiber_nl), one_minus_eps, best_residual
    real(mytype) :: x_best(fiber_nl)
    integer :: k, max_restarts, stagnation_count

    ! Warm start: keep incoming x as initial guess rather than forcing x=0.
    call apply_bending_backward_euler_operator_scalar(x, ax, dt_b, gamma_b)
    r = rhs - ax
    rhs_norm = sqrt(max(dot_product(rhs, rhs), 0._mytype))
    target = tol * max(rhs_norm, 1._mytype)
    final_residual = sqrt(max(dot_product(r, r), 0._mytype))
    it_used = 0
    restart_count = 0
    breakdown_detected = .false.
    failure_mode = 0
    x_best = x
    best_residual = final_residual
    converged = (final_residual <= target)
    if (converged) return

    eps_break = 1.0e-30_mytype
    one_minus_eps = 0.999_mytype
    max_restarts = 2
    stagnation_count = 0
    call build_bending_operator_diagonal(diagA, dt_b, gamma_b)

    do
      rhat = r
      p = 0._mytype
      v = 0._mytype
      rho_old = 1._mytype
      alpha = 1._mytype
      omega = 1._mytype

      do k = 1, maxit
        rho = dot_product(rhat, r)
        if (abs(rho) <= eps_break) then
          breakdown_detected = .true.
          failure_mode = 2
          final_residual = sqrt(max(dot_product(r, r), 0._mytype))
          if (final_residual <= target) then
            converged = .true.
            x = x_best
            return
          endif
          exit
        endif

        if (it_used == 0 .or. k == 1) then
          p = r
        else
          if (abs(omega) <= eps_break) then
            breakdown_detected = .true.
            failure_mode = 4
            final_residual = sqrt(max(dot_product(r, r), 0._mytype))
            if (final_residual <= target) then
              converged = .true.
              x = x_best
              return
            endif
            exit
          endif
          beta = (rho / rho_old) * (alpha / omega)
          p = r + beta * (p - omega * v)
        endif

        call apply_jacobi_preconditioner(p, phat, diagA)
        call apply_bending_backward_euler_operator_scalar(phat, v, dt_b, gamma_b)
        rhat_v = dot_product(rhat, v)
        if (abs(rhat_v) <= eps_break) then
          breakdown_detected = .true.
          failure_mode = 3
          final_residual = sqrt(max(dot_product(r, r), 0._mytype))
          if (final_residual <= target) then
            converged = .true.
            x = x_best
            return
          endif
          exit
        endif
        alpha = rho / rhat_v
        s = r - alpha * v

        old_residual = final_residual
        it_used = it_used + 1
        final_residual = sqrt(max(dot_product(s, s), 0._mytype))
        if (final_residual <= target) then
          x = x + alpha * phat
          x_best = x
          converged = .true.
          return
        endif

        call apply_jacobi_preconditioner(s, shat, diagA)
        call apply_bending_backward_euler_operator_scalar(shat, t, dt_b, gamma_b)
        tt = dot_product(t, t)
        if (tt <= eps_break) then
          breakdown_detected = .true.
          failure_mode = 4
          final_residual = sqrt(max(dot_product(s, s), 0._mytype))
          if (final_residual <= target) then
            x = x + alpha * phat
            x_best = x
            converged = .true.
            return
          endif
          exit
        endif
        omega = dot_product(t, s) / tt
        if (abs(omega) <= eps_break) then
          breakdown_detected = .true.
          failure_mode = 4
          final_residual = sqrt(max(dot_product(s, s), 0._mytype))
          if (final_residual <= target) then
            x = x + alpha * phat
            x_best = x
            converged = .true.
            return
          endif
          exit
        endif
        x = x + alpha * phat + omega * shat
        r = s - omega * t
        final_residual = sqrt(max(dot_product(r, r), 0._mytype))
        if (final_residual < best_residual) then
          best_residual = final_residual
          x_best = x
        endif
        if (final_residual <= target) then
          converged = .true.
          return
        endif

        if (final_residual >= one_minus_eps * old_residual) then
          stagnation_count = stagnation_count + 1
        else
          stagnation_count = 0
        endif
        if (stagnation_count >= 8) then
          failure_mode = 5
          exit
        endif
        if (it_used >= maxit) exit
        rho_old = rho
      enddo

      if (final_residual <= target) then
        converged = .true.
        return
      endif
      if (it_used >= maxit) exit
      if (restart_count >= max_restarts) exit

      ! Limited restart from current best iterate for robustness.
      restart_count = restart_count + 1
      x = x_best
      call apply_bending_backward_euler_operator_scalar(x, ax, dt_b, gamma_b)
      r = rhs - ax
      final_residual = sqrt(max(dot_product(r, r), 0._mytype))
      if (final_residual <= target) then
        converged = .true.
        return
      endif
      stagnation_count = 0
    enddo
    if (failure_mode == 0) failure_mode = 1
  end subroutine solve_bending_scalar_bicgstab

  subroutine solve_free_end_implicit_bending_rhs_scalar(rhs, x, alpha, gamma_b, it_used, final_residual, converged, &
       breakdown_detected, restart_count, failure_mode)
    ! Generic free-end implicit bending solve:
    !   (I + alpha * gamma_b * D4) x = rhs
    ! Reuses the validated Step 3.3 operator/solver route.
    real(mytype), intent(in) :: rhs(fiber_nl), alpha, gamma_b
    real(mytype), intent(inout) :: x(fiber_nl)
    real(mytype), intent(out) :: final_residual
    integer, intent(out) :: it_used, restart_count, failure_mode
    logical, intent(out) :: converged, breakdown_detected
    real(mytype) :: amat(fiber_nl, fiber_nl)

    if (fiber_flex_bending_linear_solver == 1) then
      call solve_bending_scalar_bicgstab(rhs, x, alpha, gamma_b, fiber_flex_bending_iter_tol_effective, &
           fiber_flex_bending_iter_maxit_effective, it_used, final_residual, converged, breakdown_detected, restart_count, failure_mode)
    else
      call build_bending_backward_euler_matrix(alpha, gamma_b, amat)
      call solve_linear_system_dense(amat, rhs, x)
      it_used = fiber_nl
      final_residual = 0._mytype
      converged = .true.
      breakdown_detected = .false.
      restart_count = 0
      failure_mode = 0
    endif
  end subroutine solve_free_end_implicit_bending_rhs_scalar

  subroutine advance_bending_backward_euler(x_old, x_new, dt_b, gamma_b, max_linear_iterations, max_linear_residual, &
       max_linear_restarts, linear_breakdown_detected, linear_failure_mode)
    ! Bending-only backward-Euler kernel:
    ! advances the isolated implicit bending operator and is not a full
    ! constrained structural dynamics step.
    real(mytype), intent(in) :: x_old(3, fiber_nl), dt_b, gamma_b
    real(mytype), intent(out) :: x_new(3, fiber_nl)
    integer, intent(out) :: max_linear_iterations
    real(mytype), intent(out) :: max_linear_residual
    integer, intent(out) :: max_linear_restarts
    integer, intent(out) :: linear_failure_mode
    logical, intent(out) :: linear_breakdown_detected
    real(mytype) :: amat(fiber_nl, fiber_nl)
    integer :: it_used, restart_count, failure_mode
    real(mytype) :: final_residual
    logical :: converged, breakdown_detected
    integer :: c
    character(len=64) :: failure_mode_label

    max_linear_iterations = 0
    max_linear_residual = 0._mytype
    max_linear_restarts = 0
    linear_failure_mode = 0
    linear_breakdown_detected = .false.
    if (fiber_flex_bending_linear_solver == 1) then
      ! Primary iterative path (BiCGSTAB); use merged effective iterative controls.
      do c = 1, 3
        x_new(c,:) = x_old(c,:)
        call solve_bending_scalar_bicgstab(x_old(c,:), x_new(c,:), dt_b, gamma_b, fiber_flex_bending_iter_tol_effective, &
             fiber_flex_bending_iter_maxit_effective, it_used, final_residual, converged, breakdown_detected, restart_count, failure_mode)
        max_linear_iterations = max(max_linear_iterations, it_used)
        max_linear_residual = max(max_linear_residual, final_residual)
        max_linear_restarts = max(max_linear_restarts, restart_count)
        linear_failure_mode = max(linear_failure_mode, failure_mode)
        linear_breakdown_detected = linear_breakdown_detected .or. breakdown_detected
        if (.not.converged) then
          select case (failure_mode)
          case (1); failure_mode_label = 'maxit'
          case (2); failure_mode_label = 'breakdown_rho'
          case (3); failure_mode_label = 'breakdown_alpha_den'
          case (4); failure_mode_label = 'breakdown_omega'
          case (5); failure_mode_label = 'stagnation'
          case default; failure_mode_label = 'unknown'
          end select
          if (nrank == 0) write(*,*) 'Error: BiCGSTAB primary path failed after restart limit due to ', &
               trim(failure_mode_label), ' in Step 3.3.'
          stop
        endif
      enddo
    else
      ! Dense route retained as reference/fallback implementation.
      call build_bending_backward_euler_matrix(dt_b, gamma_b, amat)
      do c = 1, 3
        call solve_linear_system_dense(amat, x_old(c,:), x_new(c,:))
      enddo
      max_linear_iterations = fiber_nl
      max_linear_residual = 0._mytype
      max_linear_restarts = 0
      linear_failure_mode = 0
      linear_breakdown_detected = .false.
    endif
  end subroutine advance_bending_backward_euler

  subroutine compute_bending_energy(x_in, gamma_b, ebend)
    ! Bending-only diagnostic energy used in Step 3.3 verification.
    ! This excludes tension contribution and full structural-energy terms.
    real(mytype), intent(in) :: x_in(3, fiber_nl), gamma_b
    real(mytype), intent(out) :: ebend
    real(mytype) :: d2(3, fiber_nl), norm2(fiber_nl)
    integer :: i

    do i = 1, 3
      call apply_free_end_d2_scalar(x_in(i,:), d2(i,:))
    enddo
    norm2 = d2(1,:)**2 + d2(2,:)**2 + d2(3,:)**2
    ebend = 0.5_mytype * gamma_b * fiber_ds * (0.5_mytype * norm2(1) + sum(norm2(2:fiber_nl-1)) + 0.5_mytype * norm2(fiber_nl))
  end subroutine compute_bending_energy

  subroutine initialize_bending_test_shape(case_id, x0)
    integer, intent(in) :: case_id
    real(mytype), intent(out) :: x0(3, fiber_nl)
    real(mytype) :: smin, q, g, amp
    integer :: l

    select case (case_id)
    case (1)
      x0 = fiber_x
    case (2)
      smin = -0.5_mytype * fiber_length
      amp = 0.1_mytype * fiber_length
      do l = 1, fiber_nl
        q = (fiber_s_ref(l) - smin) / fiber_length
        g = q**4 * (1._mytype - q)**4
        x0(1,l) = fiber_s_ref(l)
        x0(2,l) = amp * g
        x0(3,l) = 0._mytype
      enddo
    case default
      if (nrank == 0) write(*,*) 'Error: fiber_flex_bending_case must be 1 or 2.'
      stop
    end select
  end subroutine initialize_bending_test_shape

  subroutine run_flexible_bending_test()
    ! Step 3.3 intermediate-kernel verification driver.
    ! This validates bending-only implicit relaxation behavior and should
    ! not be interpreted as complete structural dynamics validation.
    real(mytype), allocatable :: xold(:,:), xnew(:,:), xinit(:,:), energy(:), maxupd(:)
    real(mytype) :: e0, ecur, tnow, max_update, straight_err, final_disp
    logical :: energy_monotone
    integer :: step, outint, linear_iters_step, max_linear_iterations_used
    integer :: linear_restarts_step, max_linear_restarts_used
    integer :: linear_failure_mode_step, final_linear_failure_mode
    real(mytype) :: linear_residual_step, max_final_linear_residual
    logical :: linear_breakdown_step, linear_breakdown_detected_any

    if (.not.fiber_active .or. .not.fiber_flexible_active .or. .not.fiber_flex_initialized .or. &
         .not.fiber_flex_bending_test_active) then
      if (nrank == 0) write(*,*) 'Error: invalid state for run_flexible_bending_test.'
      stop
    endif

    allocate(xold(3,fiber_nl), xnew(3,fiber_nl), xinit(3,fiber_nl), &
         energy(0:fiber_flex_bending_nsteps), maxupd(0:fiber_flex_bending_nsteps))

    call initialize_bending_test_shape(fiber_flex_bending_case, xold)
    xinit = xold
    call write_fiber_flex_bending_initial_points(xold)
    call compute_bending_energy(xold, fiber_bending_gamma, e0)
    energy(0) = e0
    maxupd(0) = 0._mytype
    max_linear_iterations_used = 0
    max_linear_restarts_used = 0
    final_linear_failure_mode = 0
    max_final_linear_residual = 0._mytype
    linear_breakdown_detected_any = .false.

    outint = max(1, fiber_flex_bending_output_interval)
    call write_fiber_flex_bending_series(0, 0._mytype, energy(0), maxupd(0), 0, 0._mytype, 0, .false., 0, .true.)
    do step = 1, fiber_flex_bending_nsteps
      call advance_bending_backward_euler(xold, xnew, fiber_flex_bending_dt, fiber_bending_gamma, &
           linear_iters_step, linear_residual_step, linear_restarts_step, linear_breakdown_step, linear_failure_mode_step)
      max_update = maxval(abs(xnew - xold))
      call compute_bending_energy(xnew, fiber_bending_gamma, ecur)
      energy(step) = ecur
      maxupd(step) = max_update
      max_linear_iterations_used = max(max_linear_iterations_used, linear_iters_step)
      max_linear_restarts_used = max(max_linear_restarts_used, linear_restarts_step)
      max_final_linear_residual = max(max_final_linear_residual, linear_residual_step)
      final_linear_failure_mode = linear_failure_mode_step
      linear_breakdown_detected_any = linear_breakdown_detected_any .or. linear_breakdown_step
      tnow = real(step, mytype) * fiber_flex_bending_dt
      if (mod(step, outint) == 0 .or. step == fiber_flex_bending_nsteps) then
        call write_fiber_flex_bending_series(step, tnow, ecur, max_update, linear_iters_step, linear_residual_step, &
             linear_restarts_step, linear_breakdown_step, linear_failure_mode_step, .false.)
      endif
      xold = xnew
    enddo

    call write_fiber_flex_bending_final_points(xold)

    energy_monotone = .true.
    do step = 1, fiber_flex_bending_nsteps
      if (energy(step) > energy(step-1) + 1.0e-12_mytype) energy_monotone = .false.
    enddo

    straight_err = 0._mytype
    final_disp = 0._mytype
    if (fiber_flex_bending_case == 1) then
      straight_err = maxval(abs(xold - xinit))
    else if (fiber_flex_bending_case == 2) then
      final_disp = sqrt(sum((xold - xinit)**2) / real(3*fiber_nl, mytype))
    endif

    call write_fiber_flex_bending_summary(fiber_flex_bending_case, fiber_flex_bending_dt, fiber_bending_gamma, &
         fiber_flex_bending_nsteps, energy(0), energy(fiber_flex_bending_nsteps), energy_monotone, &
         maxupd(fiber_flex_bending_nsteps), straight_err, final_disp, fiber_flex_bending_linear_solver, &
         fiber_flex_bending_iter_tol_effective, fiber_flex_bending_iter_maxit_effective, &
         max_linear_iterations_used, max_final_linear_residual, &
         max_linear_restarts_used, linear_breakdown_detected_any, final_linear_failure_mode)

    deallocate(xold, xnew, xinit, energy, maxupd)
  end subroutine run_flexible_bending_test

end module fiber_flex_bending
