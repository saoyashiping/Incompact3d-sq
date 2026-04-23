!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_flex_constraint

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use fiber_types, only : fiber_active, fiber_flexible_active, fiber_flex_initialized, fiber_nl, fiber_length, fiber_ds, &
       fiber_s_ref, fiber_x, fiber_flex_constraint_test_active, fiber_flex_constraint_case, fiber_flex_constraint_nsteps, &
       fiber_flex_constraint_output_interval, fiber_flex_constraint_dt, fiber_flex_constraint_force_amp, &
       fiber_structure_rho_tilde
  use fiber_types, only : fiber_bending_gamma
  use fiber_flex_bending, only : solve_free_end_implicit_bending_rhs_scalar
  use fiber_flex_ops, only : apply_free_end_bending_operator_scalar, apply_free_end_d2_scalar
  use fiber_io, only : write_fiber_flex_constraint_initial_points, write_fiber_flex_constraint_final_points, &
       write_fiber_flex_constraint_tension_last, write_fiber_flex_constraint_series, write_fiber_flex_constraint_summary

  implicit none

contains

  ! Step 3.4 refinement:
  ! Huang-style staggered-Lagrangian semantics (node X_i, half-node T_{i+1/2})
  ! for standalone constraint+tension layer only (no fluid coupling).

  subroutine lag_diff_node_to_half_centered(x_node, x_half)
    ! Node -> half centered difference:
    !   x_half(i+1/2) = (x_node(i+1) - x_node(i)) / ds.
    real(mytype), intent(in) :: x_node(3,fiber_nl)
    real(mytype), intent(out) :: x_half(3,fiber_nl-1)
    integer :: i
    do i = 1, fiber_nl-1
      x_half(:,i) = (x_node(:,i+1) - x_node(:,i)) / fiber_ds
    enddo
  end subroutine lag_diff_node_to_half_centered

  subroutine lag_diff_half_to_node_forward(q_half, dq_node)
    ! Half -> node forward difference:
    !   dq_node(i) = (q_half(i+1/2) - q_half(i-1/2)) / ds
    ! with explicit free-end closure q_half(1/2)=q_half(N+1/2)=0.
    real(mytype), intent(in) :: q_half(3,fiber_nl-1)
    real(mytype), intent(out) :: dq_node(3,fiber_nl)
    integer :: i
    dq_node = 0._mytype
    dq_node(:,1) = q_half(:,1) / fiber_ds
    do i = 2, fiber_nl-1
      dq_node(:,i) = (q_half(:,i) - q_half(:,i-1)) / fiber_ds
    enddo
    dq_node(:,fiber_nl) = -q_half(:,fiber_nl-1) / fiber_ds
  end subroutine lag_diff_half_to_node_forward

  subroutine lag_diff_half_to_node_backward(q_half, dq_node)
    ! Half -> node backward difference:
    !   dq_node(i) = (q_half(i+1/2) - q_half(i-1/2)) / ds
    ! same discrete node result under this uniform 1D staggered grid,
    ! but kept as a distinct helper to keep D_s^+/D_s^- roles explicit.
    real(mytype), intent(in) :: q_half(3,fiber_nl-1)
    real(mytype), intent(out) :: dq_node(3,fiber_nl)
    call lag_diff_half_to_node_forward(q_half, dq_node)
  end subroutine lag_diff_half_to_node_backward

  subroutine lag_ds_center_node_to_half(x_node, xs_half)
    ! D_s^0 on node-centered X_i to half-grid (i+1/2):
    !   (D_s^0 X)_{i+1/2} = (X_{i+1} - X_i)/ds.
    real(mytype), intent(in) :: x_node(3,fiber_nl)
    real(mytype), intent(out) :: xs_half(3,fiber_nl-1)
    call lag_diff_node_to_half_centered(x_node, xs_half)
  end subroutine lag_ds_center_node_to_half

  subroutine lag_ds_plus_node(x_node, x_plus)
    ! Huang-style D_s^+ on node-centered field, represented on half-grid:
    !   (D_s^+ X)_{i+1/2} = (X_{i+1} - X_i)/ds.
    real(mytype), intent(in) :: x_node(3,fiber_nl)
    real(mytype), intent(out) :: x_plus(3,fiber_nl-1)
    call lag_diff_node_to_half_centered(x_node, x_plus)
  end subroutine lag_ds_plus_node

  subroutine lag_ds_minus_node(x_node, x_minus)
    ! Huang-style D_s^- expressed on half-grid.
    ! Under this uniform single-fiber staggered layout, it shares the same
    ! nearest-neighbor node->half stencil as D_s^+, but is kept separate
    ! to preserve explicit operator semantics and reviewability.
    real(mytype), intent(in) :: x_node(3,fiber_nl)
    real(mytype), intent(out) :: x_minus(3,fiber_nl-1)
    call lag_diff_node_to_half_centered(x_node, x_minus)
  end subroutine lag_ds_minus_node

  subroutine lag_dss_node_free_end(x_node, xss_node)
    ! D_s^+ D_s^- on nodes under free-end closure.
    ! This uses the validated Step 3.2 D2 operator and is semantically
    ! separate from bending-force operators.
    real(mytype), intent(in) :: x_node(3,fiber_nl)
    real(mytype), intent(out) :: xss_node(3,fiber_nl)
    integer :: c
    do c = 1, 3
      call apply_free_end_d2_scalar(x_node(c,:), xss_node(c,:))
    enddo
  end subroutine lag_dss_node_free_end

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

  subroutine compute_tension_term_huang_style(t_half, x_ref, ftension)
    ! Huang-style tension force term:
    !   d_s ( T d_s X )
    ! with T stored at half-grid internal unknowns and explicit zero-end BC:
    !   T_{1/2} = 0, T_{N+1/2} = 0.
    real(mytype), intent(in) :: t_half(fiber_nl-1), x_ref(3,fiber_nl)
    real(mytype), intent(out) :: ftension(3,fiber_nl)
    real(mytype) :: xs_half(3,fiber_nl-1), flux_half(3,fiber_nl-1)
    integer :: i

    call lag_ds_center_node_to_half(x_ref, xs_half)
    do i = 1, fiber_nl-1
      flux_half(:,i) = t_half(i) * xs_half(:,i)
    enddo
    call lag_diff_half_to_node_forward(flux_half, ftension)
  end subroutine compute_tension_term_huang_style

  subroutine apply_tension_poisson_operator(t_half_unknown, x_ref, lhs_t)
    real(mytype), intent(in) :: t_half_unknown(fiber_nl-1), x_ref(3,fiber_nl)
    real(mytype), intent(out) :: lhs_t(fiber_nl-1)
    real(mytype) :: ft(3,fiber_nl), xs_half(3,fiber_nl-1), dft(3,fiber_nl-1)
    integer :: i

    call compute_tension_term_huang_style(t_half_unknown, x_ref, ft)
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

  subroutine build_tension_rhs_huang_style(xn, xnm1, x_ref, fb_ref, fext_ref, dt_c, rhs_t, drift_term, vel_term, force_term)
    real(mytype), intent(in) :: xn(3,fiber_nl), xnm1(3,fiber_nl), x_ref(3,fiber_nl)
    real(mytype), intent(in) :: fb_ref(3,fiber_nl), fext_ref(3,fiber_nl), dt_c
    real(mytype), intent(out) :: rhs_t(fiber_nl-1), drift_term(fiber_nl-1), vel_term(fiber_nl-1), force_term(fiber_nl-1)
    real(mytype) :: un(3,fiber_nl), xs_plus(3,fiber_nl-1), u_s_half(3,fiber_nl-1), dforce(3,fiber_nl-1)
    real(mytype) :: drift, velterm, forceterm, kappa_drift, inv_rho
    integer :: i

    call compute_centerline_velocity(xn, xnm1, dt_c, un)
    call lag_ds_plus_node(x_ref, xs_plus)
    call lag_ds_minus_node(un, u_s_half)

    ! Drift control follows Huang-style differentiated inextensibility:
    !   |D_s X_ref|^2 - 1 -> damped through a dt-scaled correction term.
    ! This is a practical discrete realization of Huang (2007)-style
    ! constraint-time-differentiation drift suppression.
    kappa_drift = 2._mytype / (dt_c * dt_c)
    inv_rho = 1._mytype / fiber_structure_rho_tilde
    do i = 1, fiber_nl-1
      ! Force term uses bending force evaluated at X_ref = 2X^n - X^{n-1},
      ! plus structural external force only (no fluid-structure forcing here).
      dforce(:,i) = (fb_ref(:,i+1) + fext_ref(:,i+1) - fb_ref(:,i) - fext_ref(:,i)) / fiber_ds
      drift = dot_product(xs_plus(:,i), xs_plus(:,i)) - 1._mytype
      ! Velocity term corresponds to discrete D_s U · D_s U, not an empirical add-on.
      velterm = dot_product(u_s_half(:,i), u_s_half(:,i))
      ! Force contribution is 2 D_s X_ref · D_s(F_b(X_ref)+F_ext).
      forceterm = 2._mytype * dot_product(xs_plus(:,i), dforce(:,i))
      drift_term(i) = -kappa_drift * drift
      vel_term(i) = -velterm * inv_rho
      force_term(i) = -forceterm * inv_rho
      rhs_t(i) = drift_term(i) + vel_term(i) + force_term(i)
    enddo
  end subroutine build_tension_rhs_huang_style

  subroutine assemble_tension_matrix_huang_style(x_ref, amat)
    ! Internal unknowns are T_{i+1/2}, i=1..N-1, with explicit end BC:
    ! T_{1/2}=0 and T_{N+1/2}=0 (not part of unknown vector).
    real(mytype), intent(in) :: x_ref(3,fiber_nl)
    real(mytype), intent(out) :: amat(fiber_nl-1,fiber_nl-1)
    real(mytype) :: evec(fiber_nl-1), col(fiber_nl-1)
    integer :: j
    do j = 1, fiber_nl-1
      evec = 0._mytype; evec(j) = 1._mytype
      call apply_tension_poisson_operator(evec, x_ref, col)
      amat(:,j) = col
    enddo
  end subroutine assemble_tension_matrix_huang_style

  subroutine solve_constraint_tension_huang_style(xn, xnm1, x_ref, fb_ref, fext_ref, dt_c, t_half, &
       max_abs_drift_term, max_abs_vel_term, max_abs_force_term)
    real(mytype), intent(in) :: xn(3,fiber_nl), xnm1(3,fiber_nl), x_ref(3,fiber_nl)
    real(mytype), intent(in) :: fb_ref(3,fiber_nl), fext_ref(3,fiber_nl), dt_c
    real(mytype), intent(out) :: t_half(fiber_nl-1)
    real(mytype), intent(out) :: max_abs_drift_term, max_abs_vel_term, max_abs_force_term
    real(mytype) :: amat(fiber_nl-1,fiber_nl-1), rhs(fiber_nl-1)
    real(mytype) :: drift_term(fiber_nl-1), vel_term(fiber_nl-1), force_term(fiber_nl-1)

    call assemble_tension_matrix_huang_style(x_ref, amat)
    call build_tension_rhs_huang_style(xn, xnm1, x_ref, fb_ref, fext_ref, dt_c, rhs, drift_term, vel_term, force_term)

    call solve_dense_small(amat, rhs, t_half, fiber_nl-1)
    max_abs_drift_term = maxval(abs(drift_term))
    max_abs_vel_term = maxval(abs(vel_term))
    max_abs_force_term = maxval(abs(force_term))
  end subroutine solve_constraint_tension_huang_style

  subroutine solve_constrained_step_with_implicit_bending(rhs_vec, x_out, dt_c, rho_tilde, gamma_b, implicit_update_norm)
    ! Solve (I + alpha*gamma_b*D4) X^{n+1} = rhs with alpha = dt_c^2/rho_tilde
    ! by reusing the validated Step 3.3 implicit bending kernel route.
    real(mytype), intent(in) :: rhs_vec(3,fiber_nl), dt_c, rho_tilde, gamma_b
    real(mytype), intent(out) :: x_out(3,fiber_nl), implicit_update_norm
    real(mytype) :: alpha, final_residual
    integer :: c, it_used, restart_count, failure_mode
    logical :: converged, breakdown_detected

    alpha = (dt_c * dt_c) / rho_tilde
    do c = 1, 3
      x_out(c,:) = rhs_vec(c,:)
      call solve_free_end_implicit_bending_rhs_scalar(rhs_vec(c,:), x_out(c,:), alpha, gamma_b, it_used, final_residual, &
           converged, breakdown_detected, restart_count, failure_mode)
      if (.not.converged) then
        if (nrank == 0) write(*,*) 'Error: implicit bending solve failed in Step 3.4.'
        stop
      endif
    enddo
    implicit_update_norm = sqrt(sum((x_out - rhs_vec)**2) / real(3*fiber_nl, mytype))
  end subroutine solve_constrained_step_with_implicit_bending

  subroutine build_constraint_correction_rhs(x_prov, dt_c, rhs_corr)
    ! Post-bending constraint correction RHS:
    ! driven only by current inextensibility violation of provisional X.
    ! This is separate from the main Huang-style time-advance RHS.
    real(mytype), intent(in) :: x_prov(3,fiber_nl), dt_c
    real(mytype), intent(out) :: rhs_corr(fiber_nl-1)
    real(mytype) :: xs_half(3,fiber_nl-1), kappa_corr
    integer :: i

    call lag_ds_center_node_to_half(x_prov, xs_half)
    kappa_corr = 2._mytype / (dt_c * dt_c)
    do i = 1, fiber_nl-1
      rhs_corr(i) = -kappa_corr * (dot_product(xs_half(:,i), xs_half(:,i)) - 1._mytype)
    enddo
  end subroutine build_constraint_correction_rhs

  subroutine apply_post_bending_constraint_correction(x_prov, x_corr, dt_c, rho_tilde, correction_scale, &
       post_bending_correction_norm, post_bending_max_inext_err_after_correction)
    ! This is a post-bending constraint correction used within the Step 3.4
    ! structure-only solver. It is not a full monolithic constrained implicit
    ! structure solve. The goal is to restore consistency between the
    ! provisional implicit-bending configuration and inextensibility.
    real(mytype), intent(in) :: x_prov(3,fiber_nl), dt_c, rho_tilde, correction_scale
    real(mytype), intent(out) :: x_corr(3,fiber_nl), post_bending_correction_norm, post_bending_max_inext_err_after_correction
    real(mytype) :: rhs_corr(fiber_nl-1), amat_corr(fiber_nl-1,fiber_nl-1), t_corr(fiber_nl-1)
    real(mytype) :: f_corr(3,fiber_nl), max_seg_after

    call build_constraint_correction_rhs(x_prov, dt_c, rhs_corr)
    call assemble_tension_matrix_huang_style(x_prov, amat_corr)
    call solve_dense_small(amat_corr, rhs_corr, t_corr, fiber_nl-1)
    call compute_tension_term_huang_style(t_corr, x_prov, f_corr)

    x_corr = x_prov + correction_scale * (dt_c * dt_c / rho_tilde) * f_corr
    post_bending_correction_norm = sqrt(sum((x_corr - x_prov)**2) / real(3*fiber_nl, mytype))
    call compute_constraint_error(x_corr, max_seg_after, post_bending_max_inext_err_after_correction)
  end subroutine apply_post_bending_constraint_correction

  subroutine advance_constrained_structure_step(xnm1, xn, xnp1, t_half, fext_n, dt_c, &
       max_abs_drift_term, max_abs_vel_term, max_abs_force_term, implicit_bending_update_norm, &
       post_bending_correction_norm, post_bending_max_inext_err_after_correction)
    real(mytype), intent(in) :: xnm1(3,fiber_nl), xn(3,fiber_nl), fext_n(3,fiber_nl), dt_c
    real(mytype), intent(out) :: xnp1(3,fiber_nl), t_half(fiber_nl-1)
    real(mytype), intent(out) :: max_abs_drift_term, max_abs_vel_term, max_abs_force_term, implicit_bending_update_norm
    real(mytype), intent(out) :: post_bending_correction_norm, post_bending_max_inext_err_after_correction
    real(mytype) :: x_ref(3,fiber_nl), fb_ref(3,fiber_nl), ftension(3,fiber_nl), rhs_final(3,fiber_nl), x_prov(3,fiber_nl)
    real(mytype), parameter :: post_bending_correction_scale = 1.0_mytype
    integer :: c

    ! Step A: predictor configuration.
    x_ref = 2._mytype * xn - xnm1
    ! Step B: evaluate predicted bending force for tension RHS/diagnostics.
    do c = 1, 3
      call apply_free_end_bending_operator_scalar(x_ref(c,:), fb_ref(c,:), fiber_bending_gamma)
    enddo
    call solve_constraint_tension_huang_style(xn, xnm1, x_ref, fb_ref, fext_n, dt_c, t_half, &
         max_abs_drift_term, max_abs_vel_term, max_abs_force_term)
    call compute_tension_term_huang_style(t_half, x_ref, ftension)

    ! Step C: final RHS excludes explicit bending push; bending is implicit in Step D.
    rhs_final = x_ref + (dt_c * dt_c / fiber_structure_rho_tilde) * (ftension + fext_n)
    ! Step D: semi-implicit constrained update with 3.3 kernel:
    !   (I + (dt^2/rho_tilde) * gamma * D4) X^{n+1} = rhs_final
    call solve_constrained_step_with_implicit_bending(rhs_final, x_prov, dt_c, fiber_structure_rho_tilde, &
         fiber_bending_gamma, implicit_bending_update_norm)
    call apply_post_bending_constraint_correction(x_prov, xnp1, dt_c, fiber_structure_rho_tilde, post_bending_correction_scale, &
         post_bending_correction_norm, post_bending_max_inext_err_after_correction)
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
      ! Constraint capability test: small-amplitude free-end-consistent shape,
      ! zero initial velocity (xnm1 = xn), no rigid-body-style forcing.
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
    real(mytype) :: max_abs_drift_term_step, max_abs_vel_term_step, max_abs_force_term_step
    real(mytype) :: max_abs_drift_term_global, max_abs_vel_term_global, max_abs_force_term_global
    real(mytype) :: implicit_bending_update_norm_step, implicit_bending_update_norm_global
    real(mytype) :: post_bending_correction_norm_step, post_bending_correction_norm_global
    real(mytype) :: post_bending_max_inext_err_after_correction_step, post_bending_max_inext_err_after_correction_global
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
    max_abs_drift_term_global = 0._mytype
    max_abs_vel_term_global = 0._mytype
    max_abs_force_term_global = 0._mytype
    implicit_bending_update_norm_global = 0._mytype
    post_bending_correction_norm_global = 0._mytype
    post_bending_max_inext_err_after_correction_global = 0._mytype
    outint = max(1, fiber_flex_constraint_output_interval)
    call write_fiber_flex_constraint_series(0, 0._mytype, 0._mytype, 0._mytype, 0._mytype, &
         0._mytype, 0._mytype, 0._mytype, 0._mytype, 0._mytype, 0._mytype, .true.)

    do step = 1, fiber_flex_constraint_nsteps
      fext = 0._mytype
      if (fiber_flex_constraint_case == 3) then
        ! Keep forcing zero or very small and spatially localized so case 3
        ! emphasizes constrained structural response rather than rigid forcing.
        do i = 1, fiber_nl
          fext(2,i) = 0.1_mytype * fiber_flex_constraint_force_amp * exp(-((real(i-1,mytype)/real(max(1,fiber_nl-1),mytype) - 0.5_mytype)**2) / 0.02_mytype)
        enddo
      endif

      call advance_constrained_structure_step(xnm1, xn, xnp1, t_half, fext, fiber_flex_constraint_dt, &
           max_abs_drift_term_step, max_abs_vel_term_step, max_abs_force_term_step, implicit_bending_update_norm_step, &
           post_bending_correction_norm_step, post_bending_max_inext_err_after_correction_step)
      call compute_constraint_error(xnp1, max_seg_err, max_inext_err)
      max_abs_tension = maxval(abs(t_half))
      max_seg_global = max(max_seg_global, max_seg_err)
      max_inext_global = max(max_inext_global, max_inext_err)
      max_t_global = max(max_t_global, max_abs_tension)
      max_abs_drift_term_global = max(max_abs_drift_term_global, max_abs_drift_term_step)
      max_abs_vel_term_global = max(max_abs_vel_term_global, max_abs_vel_term_step)
      max_abs_force_term_global = max(max_abs_force_term_global, max_abs_force_term_step)
      implicit_bending_update_norm_global = max(implicit_bending_update_norm_global, implicit_bending_update_norm_step)
      post_bending_correction_norm_global = max(post_bending_correction_norm_global, post_bending_correction_norm_step)
      post_bending_max_inext_err_after_correction_global = max(post_bending_max_inext_err_after_correction_global, &
           post_bending_max_inext_err_after_correction_step)

      tnow = real(step,mytype) * fiber_flex_constraint_dt
      if (mod(step,outint) == 0 .or. step == fiber_flex_constraint_nsteps) then
        call write_fiber_flex_constraint_series(step, tnow, max_seg_err, max_inext_err, max_abs_tension, &
             max_abs_drift_term_step, max_abs_vel_term_step, max_abs_force_term_step, &
             implicit_bending_update_norm_step, post_bending_correction_norm_step, &
             post_bending_max_inext_err_after_correction_step, .false.)
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
         fiber_flex_constraint_nsteps, max_seg_global, max_inext_global, max_t_global, case_metric, &
         max_abs_drift_term_global, max_abs_vel_term_global, max_abs_force_term_global, &
         implicit_bending_update_norm_global, post_bending_correction_norm_global, &
         post_bending_max_inext_err_after_correction_global, 1.0_mytype)

    deallocate(xnm1, xn, xnp1, xinit, fext, t_half)
  end subroutine run_flexible_constraint_test

end module fiber_flex_constraint
