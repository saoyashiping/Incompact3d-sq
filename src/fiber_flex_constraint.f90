!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_flex_constraint

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use fiber_types, only : fiber_active, fiber_flexible_active, fiber_flex_initialized, fiber_nl, fiber_length, fiber_ds, &
       fiber_center, fiber_direction, &
       fiber_s_ref, fiber_x, fiber_flex_constraint_test_active, fiber_flex_constraint_case, fiber_flex_constraint_nsteps, &
       fiber_flex_constraint_output_interval, fiber_flex_constraint_dt, fiber_flex_constraint_force_amp, &
       fiber_flex_structure_test_active, fiber_flex_structure_case, fiber_flex_structure_nsteps, &
       fiber_flex_structure_output_interval, fiber_flex_structure_dt, fiber_flex_structure_force_mode, &
       fiber_flex_structure_force_amp, fiber_flex_structure_force_omega, fiber_flex_structure_force_direction, &
       fiber_flex_structure_initial_shape_amp, fiber_flex_structure_adaptive_substep_active, &
       fiber_flex_structure_max_substep_splits, &
       fiber_structure_rho_tilde, fiber_flex_constraint_outer_maxit, fiber_flex_constraint_outer_tol_x, &
       fiber_flex_constraint_outer_tol_g, fiber_flex_constraint_outer_tol_rx_rel, &
       fiber_flex_constraint_line_search_active, fiber_flex_constraint_line_search_beta, &
       fiber_flex_constraint_line_search_max_backtracks, fiber_flex_constraint_tension_warm_start_active, &
       fiber_tension_half_prevstep
  use fiber_types, only : fiber_bending_gamma
  use fiber_flex_bending, only : solve_free_end_implicit_bending_rhs_scalar
  use fiber_flex_ops, only : apply_free_end_bending_operator_scalar, apply_free_end_d2_scalar
  use fiber_io, only : write_fiber_flex_constraint_initial_points, write_fiber_flex_constraint_final_points, &
       write_fiber_flex_constraint_tension_last, write_fiber_flex_constraint_series, write_fiber_flex_constraint_summary, &
       write_fiber_flex_structure_series, write_fiber_flex_structure_summary

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

  pure function max_abs_array_2d(a) result(vmax)
    real(mytype), intent(in) :: a(:,:)
    real(mytype) :: vmax
    vmax = maxval(abs(a))
  end function max_abs_array_2d

  subroutine check_implicit_projection_residual(solution, rhs, alpha, gamma_b, projection_rhs_negligible_tol, projection_rel_tol, &
       res_norm, rel_res, projection_negligible, projection_ok)
    real(mytype), intent(in) :: solution(3,fiber_nl), rhs(3,fiber_nl), alpha, gamma_b
    real(mytype), intent(in) :: projection_rhs_negligible_tol, projection_rel_tol
    real(mytype), intent(out) :: res_norm, rel_res
    logical, intent(out) :: projection_negligible, projection_ok
    real(mytype) :: lhs(3,fiber_nl), d4sol(3,fiber_nl), rhs_norm, tiny_norm
    integer :: c

    do c = 1, 3
      call apply_free_end_bending_operator_scalar(solution(c,:), d4sol(c,:))
    enddo
    lhs = solution + alpha * gamma_b * d4sol
    rhs_norm = max_abs_array_2d(rhs)
    res_norm = max_abs_array_2d(lhs - rhs)
    tiny_norm = 1.0e-30_mytype
    rel_res = res_norm / max(rhs_norm, tiny_norm)
    projection_negligible = rhs_norm <= projection_rhs_negligible_tol
    if (projection_negligible) then
      projection_ok = .true.
    else
      projection_ok = rel_res <= projection_rel_tol
    endif
  end subroutine check_implicit_projection_residual

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
    ! Deprecated 3.4 predictor-corrector / post-correction reference path.
    ! Not used by the default coupled final-time solver.
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

  subroutine build_coupled_structure_residual(x_trial, t_half_trial, xn, xnm1, rx, rg, fext_in, dt_c, rho_tilde, gamma_b)
    ! Coupled final-time residual for Step 3.4 structure-only solve:
    ! unknowns are node X and half-grid internal tension T (zero-end BC).
    real(mytype), intent(in) :: x_trial(3,fiber_nl), t_half_trial(fiber_nl-1), xn(3,fiber_nl), xnm1(3,fiber_nl)
    real(mytype), intent(in) :: fext_in(3,fiber_nl), dt_c, rho_tilde, gamma_b
    real(mytype), intent(out) :: rx(3,fiber_nl), rg(fiber_nl-1)
    real(mytype) :: ftension(3,fiber_nl), fb_trial(3,fiber_nl), xs_half(3,fiber_nl-1)
    integer :: c, i

    call compute_tension_term_huang_style(t_half_trial, x_trial, ftension)
    do c = 1, 3
      call apply_free_end_bending_operator_scalar(x_trial(c,:), fb_trial(c,:), gamma_b)
    enddo
    rx = (rho_tilde/(dt_c*dt_c)) * (x_trial - 2._mytype * xn + xnm1) - ftension - fb_trial - fext_in

    call lag_ds_center_node_to_half(x_trial, xs_half)
    do i = 1, fiber_nl-1
      rg(i) = dot_product(xs_half(:,i), xs_half(:,i)) - 1._mytype
    enddo
  end subroutine build_coupled_structure_residual

  subroutine compute_scaled_momentum_residual_reference(x_trial, t_half_trial, xn, xnm1, fext_in, dt_c, rho_tilde, gamma_b, rx_scale)
    ! Deprecated for Step 3.5 default path.
    ! Kept only for older 3.4 reference routines.
    real(mytype), intent(in) :: x_trial(3,fiber_nl), t_half_trial(fiber_nl-1), xn(3,fiber_nl), xnm1(3,fiber_nl)
    real(mytype), intent(in) :: fext_in(3,fiber_nl), dt_c, rho_tilde, gamma_b
    real(mytype), intent(out) :: rx_scale
    real(mytype) :: ftension(3,fiber_nl), fb_trial(3,fiber_nl), inertia_ref(3,fiber_nl)
    real(mytype) :: inertia_ref_max, ftension_max, fb_max, fext_max
    integer :: c

    inertia_ref = (rho_tilde/(dt_c*dt_c)) * (2._mytype * xn - xnm1)
    inertia_ref_max = maxval(abs(inertia_ref))
    call compute_tension_term_huang_style(t_half_trial, x_trial, ftension)
    ftension_max = maxval(abs(ftension))
    do c = 1, 3
      call apply_free_end_bending_operator_scalar(x_trial(c,:), fb_trial(c,:), gamma_b)
    enddo
    fb_max = maxval(abs(fb_trial))
    fext_max = maxval(abs(fext_in))
    ! this is a scaled momentum residual used for coupled outer-iteration convergence assessment
    rx_scale = max(1._mytype, max(inertia_ref_max, max(ftension_max, max(fb_max, fext_max))))
  end subroutine compute_scaled_momentum_residual_reference

  subroutine compute_momentum_residual_scale_final_time(x_trial, t_half_trial, xn, xnm1, fext_in, dt_c, rho_tilde, gamma_b, rx_scale)
    ! momentum residual scale is translation-invariant and based on final-time residual term magnitudes
    real(mytype), intent(in) :: x_trial(3,fiber_nl), t_half_trial(fiber_nl-1), xn(3,fiber_nl), xnm1(3,fiber_nl)
    real(mytype), intent(in) :: fext_in(3,fiber_nl), dt_c, rho_tilde, gamma_b
    real(mytype), intent(out) :: rx_scale
    real(mytype) :: ftension(3,fiber_nl), fb_trial(3,fiber_nl), inertia_term(3,fiber_nl)
    real(mytype) :: inertia_max, ftension_max, fb_max, fext_max
    integer :: c

    inertia_term = (rho_tilde/(dt_c*dt_c)) * (x_trial - 2._mytype * xn + xnm1)
    inertia_max = maxval(abs(inertia_term))
    call compute_tension_term_huang_style(t_half_trial, x_trial, ftension)
    ftension_max = maxval(abs(ftension))
    do c = 1, 3
      call apply_free_end_bending_operator_scalar(x_trial(c,:), fb_trial(c,:), gamma_b)
    enddo
    fb_max = maxval(abs(fb_trial))
    fext_max = maxval(abs(fext_in))
    rx_scale = max(1._mytype, max(inertia_max, max(ftension_max, max(fb_max, fext_max))))
  end subroutine compute_momentum_residual_scale_final_time

  subroutine pack_structure_unknown(x, t_half, uvec)
    ! Unknown layout:
    !   1..(3*fiber_nl): X(c,l) in component-major order
    !   (3*fiber_nl+1)..(4*fiber_nl-1): T_half(i), i=1..fiber_nl-1
    real(mytype), intent(in) :: x(3,fiber_nl), t_half(fiber_nl-1)
    real(mytype), intent(out) :: uvec(4*fiber_nl-1)
    integer :: l, c, k
    k = 0
    do c = 1, 3
      do l = 1, fiber_nl
        k = k + 1
        uvec(k) = x(c,l)
      enddo
    enddo
    do l = 1, fiber_nl-1
      k = k + 1
      uvec(k) = t_half(l)
    enddo
  end subroutine pack_structure_unknown

  subroutine unpack_structure_unknown(uvec, x, t_half)
    real(mytype), intent(in) :: uvec(4*fiber_nl-1)
    real(mytype), intent(out) :: x(3,fiber_nl), t_half(fiber_nl-1)
    integer :: l, c, k
    k = 0
    do c = 1, 3
      do l = 1, fiber_nl
        k = k + 1
        x(c,l) = uvec(k)
      enddo
    enddo
    do l = 1, fiber_nl-1
      k = k + 1
      t_half(l) = uvec(k)
    enddo
  end subroutine unpack_structure_unknown

  subroutine build_structure_newton_residual_vector(uvec, xnm1, xn, fext_in, dt_s, rho_tilde, gamma_b, rvec)
    real(mytype), intent(in) :: uvec(4*fiber_nl-1), xnm1(3,fiber_nl), xn(3,fiber_nl), fext_in(3,fiber_nl)
    real(mytype), intent(in) :: dt_s, rho_tilde, gamma_b
    real(mytype), intent(out) :: rvec(4*fiber_nl-1)
    real(mytype) :: x_trial(3,fiber_nl), t_trial(fiber_nl-1), rx(3,fiber_nl), rg(fiber_nl-1)
    integer :: l, c, k
    call unpack_structure_unknown(uvec, x_trial, t_trial)
    call build_coupled_structure_residual(x_trial, t_trial, xn, xnm1, rx, rg, fext_in, dt_s, rho_tilde, gamma_b)
    k = 0
    do c = 1, 3
      do l = 1, fiber_nl
        k = k + 1
        rvec(k) = rx(c,l)
      enddo
    enddo
    do l = 1, fiber_nl-1
      k = k + 1
      rvec(k) = rg(l)
    enddo
  end subroutine build_structure_newton_residual_vector

  subroutine build_structure_newton_scaled_residual_vector(uvec, xnm1, xn, fext_in, dt_s, rho_tilde, gamma_b, rvec_scaled, &
       rx_scale, rg_scale, rvec_unscaled)
    ! scaled residual is used for Newton/KKT solve conditioning; unscaled residual is retained for diagnostics
    real(mytype), intent(in) :: uvec(4*fiber_nl-1), xnm1(3,fiber_nl), xn(3,fiber_nl), fext_in(3,fiber_nl)
    real(mytype), intent(in) :: dt_s, rho_tilde, gamma_b
    real(mytype), intent(out) :: rvec_scaled(4*fiber_nl-1), rx_scale, rg_scale
    real(mytype), intent(out), optional :: rvec_unscaled(4*fiber_nl-1)
    real(mytype) :: x_trial(3,fiber_nl), t_trial(fiber_nl-1), rx(3,fiber_nl), rg(fiber_nl-1), rraw(4*fiber_nl-1)
    integer :: l, c, k

    call unpack_structure_unknown(uvec, x_trial, t_trial)
    call build_coupled_structure_residual(x_trial, t_trial, xn, xnm1, rx, rg, fext_in, dt_s, rho_tilde, gamma_b)
    call compute_momentum_residual_scale_final_time(x_trial, t_trial, xn, xnm1, fext_in, dt_s, rho_tilde, gamma_b, rx_scale)
    rg_scale = 1._mytype
    k = 0
    do c = 1, 3
      do l = 1, fiber_nl
        k = k + 1
        rraw(k) = rx(c,l)
        rvec_scaled(k) = rraw(k) / rx_scale
      enddo
    enddo
    do l = 1, fiber_nl-1
      k = k + 1
      rraw(k) = rg(l)
      rvec_scaled(k) = rraw(k) / rg_scale
    enddo
    if (present(rvec_unscaled)) rvec_unscaled = rraw
  end subroutine build_structure_newton_scaled_residual_vector

  subroutine build_newton_unknown_scaling(uvec, fext_in, gamma_b, u_scale, x_scale, t_scale)
    real(mytype), intent(in) :: uvec(4*fiber_nl-1), fext_in(3,fiber_nl), gamma_b
    real(mytype), intent(out) :: u_scale(4*fiber_nl-1), x_scale, t_scale
    real(mytype) :: max_abs_tension, max_abs_fext, tiny_s
    integer :: j

    x_scale = max(fiber_length, 1._mytype)
    max_abs_tension = maxval(abs(uvec(3*fiber_nl+1:4*fiber_nl-1)))
    max_abs_fext = maxval(abs(fext_in))
    tiny_s = max(fiber_ds*fiber_ds, 1.0e-14_mytype)
    t_scale = max(1._mytype, max(max_abs_tension, max(max_abs_fext * fiber_length, gamma_b / tiny_s)))
    do j = 1, 3*fiber_nl
      u_scale(j) = x_scale
    enddo
    do j = 3*fiber_nl+1, 4*fiber_nl-1
      u_scale(j) = t_scale
    enddo
  end subroutine build_newton_unknown_scaling

  subroutine solve_dense_general_pivoted(ain, b, x, n, info, min_abs_pivot, max_abs_pivot)
    integer, intent(in) :: n
    real(mytype), intent(in) :: ain(n,n), b(n)
    real(mytype), intent(out) :: x(n)
    integer, intent(out) :: info
    real(mytype), intent(out) :: min_abs_pivot, max_abs_pivot
    real(mytype) :: a(n,n), rhs(n), rowtmp(n), piv, fac, tmp, pivot_tol
    integer :: i, k, p

    a = ain; rhs = b; x = 0._mytype
    info = 0
    min_abs_pivot = huge(1._mytype)
    max_abs_pivot = 0._mytype
    pivot_tol = 1.0e-14_mytype
    do k = 1, n
      p = k; piv = abs(a(k,k))
      do i = k+1, n
        if (abs(a(i,k)) > piv) then
          piv = abs(a(i,k)); p = i
        endif
      enddo
      min_abs_pivot = min(min_abs_pivot, piv)
      max_abs_pivot = max(max_abs_pivot, piv)
      if (piv <= pivot_tol) then
        info = 1
        return
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
      if (abs(a(i,i)) <= pivot_tol) then
        info = 2
        return
      endif
      if (i < n) then
        x(i) = (rhs(i) - sum(a(i,i+1:n)*x(i+1:n))) / a(i,i)
      else
        x(i) = rhs(i) / a(i,i)
      endif
    enddo
  end subroutine solve_dense_general_pivoted

  subroutine solve_newton_linear_system_with_lm_fallback(jac, rvec, u_scale, du, method, lm_used, mu_final, &
       min_abs_pivot, max_abs_pivot, failure_mode)
    real(mytype), intent(in) :: jac(4*fiber_nl-1,4*fiber_nl-1), rvec(4*fiber_nl-1), u_scale(4*fiber_nl-1)
    real(mytype), intent(out) :: du(4*fiber_nl-1), mu_final, min_abs_pivot, max_abs_pivot
    character(len=*), intent(out) :: method, failure_mode
    logical, intent(out) :: lm_used
    real(mytype) :: a(4*fiber_nl-1,4*fiber_nl-1), dz(4*fiber_nl-1), rhs(4*fiber_nl-1)
    real(mytype) :: ata(4*fiber_nl-1,4*fiber_nl-1), ata0(4*fiber_nl-1,4*fiber_nl-1), atr(4*fiber_nl-1), mu
    integer :: n, j, it_lm, info
    real(mytype) :: minp, maxp, minp_lm, maxp_lm

    n = 4*fiber_nl-1
    do j = 1, n
      a(:,j) = jac(:,j) * u_scale(j)
    enddo
    rhs = -rvec
    call solve_dense_general_pivoted(a, rhs, dz, n, info, minp, maxp)
    min_abs_pivot = minp
    max_abs_pivot = maxp
    if (info == 0) then
      du = u_scale * dz
      method = 'pivoted_lu'
      lm_used = .false.
      mu_final = 0._mytype
      failure_mode = 'none'
      return
    endif

    ata0 = matmul(transpose(a), a)
    atr = -matmul(transpose(a), rvec)
    mu = 1.0e-12_mytype
    lm_used = .true.
    do it_lm = 1, 8
      ata = ata0
      do j = 1, n
        ata(j,j) = ata(j,j) + mu
      enddo
      call solve_dense_general_pivoted(ata, atr, dz, n, info, minp_lm, maxp_lm)
      if (info == 0) then
        du = u_scale * dz
        method = 'lm_fallback'
        mu_final = mu
        min_abs_pivot = min(min_abs_pivot, minp_lm)
        max_abs_pivot = max(max_abs_pivot, maxp_lm)
        failure_mode = 'lu_near_singular_recovered_by_lm'
        return
      endif
      mu = 10._mytype * mu
    enddo
    du = 0._mytype
    method = 'failed'
    mu_final = mu
    failure_mode = 'lu_and_lm_failed'
  end subroutine solve_newton_linear_system_with_lm_fallback

  subroutine compute_half_segment_constraint_residual(x_in, g_out, max_abs_g, max_seg_err, max_inext_err)
    real(mytype), intent(in) :: x_in(3,fiber_nl)
    real(mytype), intent(out) :: g_out(fiber_nl-1), max_abs_g, max_seg_err, max_inext_err
    real(mytype) :: dsvec(3), seglen
    integer :: i
    max_abs_g = 0._mytype
    max_seg_err = 0._mytype
    max_inext_err = 0._mytype
    do i = 1, fiber_nl-1
      dsvec = (x_in(:,i+1)-x_in(:,i))/fiber_ds
      g_out(i) = dot_product(dsvec, dsvec) - 1._mytype
      seglen = sqrt(sum((x_in(:,i+1)-x_in(:,i))**2))
      max_abs_g = max(max_abs_g, abs(g_out(i)))
      max_seg_err = max(max_seg_err, abs(seglen - fiber_ds))
      max_inext_err = max(max_inext_err, abs(g_out(i)))
    enddo
  end subroutine compute_half_segment_constraint_residual

  subroutine solve_position_given_tension_schur(xnm1, xn, x_geom, t_half, fext_in, dt_s, rho_tilde, gamma_b, x_new, solve_ok)
    ! implicit bending kernel absorbs the high-order stiffness and is reused from Step 3.3
    real(mytype), intent(in) :: xnm1(3,fiber_nl), xn(3,fiber_nl), x_geom(3,fiber_nl), t_half(fiber_nl-1)
    real(mytype), intent(in) :: fext_in(3,fiber_nl), dt_s, rho_tilde, gamma_b
    real(mytype), intent(out) :: x_new(3,fiber_nl)
    logical, intent(out) :: solve_ok
    real(mytype) :: ftension(3,fiber_nl), rhs(3,fiber_nl), alpha, dxmax
    integer :: c, it_used, restart_count, failure_mode
    logical :: converged, breakdown_detected

    alpha = dt_s * dt_s / rho_tilde
    solve_ok = .true.
    call compute_tension_term_huang_style(t_half, x_geom, ftension)
    rhs = 2._mytype * xn - xnm1 + alpha * (fext_in + ftension)
    do c = 1, 3
      x_new(c,:) = x_geom(c,:)
      call solve_free_end_implicit_bending_rhs_scalar(rhs(c,:), x_new(c,:), alpha, gamma_b, it_used, dxmax, &
           converged, breakdown_detected, restart_count, failure_mode)
      if (.not.converged) then
        solve_ok = .false.
        return
      endif
    enddo
  end subroutine solve_position_given_tension_schur

  subroutine build_tension_schur_matrix(x_geom, dt_s, rho_tilde, gamma_b, smat, schur_matrix_ok, &
       schur_matrix_failure_mode, schur_matrix_failed_column)
    real(mytype), intent(in) :: x_geom(3,fiber_nl), dt_s, rho_tilde, gamma_b
    real(mytype), intent(out) :: smat(fiber_nl-1,fiber_nl-1)
    logical, intent(out) :: schur_matrix_ok
    character(len=*), intent(out) :: schur_matrix_failure_mode
    integer, intent(out) :: schur_matrix_failed_column
    real(mytype) :: et(fiber_nl-1), delta_t(fiber_nl-1), delta_ft(3,fiber_nl), rhs_delta(3,fiber_nl), delta_x(3,fiber_nl)
    real(mytype) :: dsx(3,fiber_nl-1), dsdx(3,fiber_nl-1), alpha
    integer :: j, i, c, it_used, restart_count, failure_mode
    real(mytype) :: dxmax
    logical :: converged, breakdown_detected

    alpha = dt_s * dt_s / rho_tilde
    schur_matrix_ok = .true.
    schur_matrix_failure_mode = 'none'
    schur_matrix_failed_column = 0
    call lag_ds_center_node_to_half(x_geom, dsx)
    do j = 1, fiber_nl-1
      delta_t = 0._mytype
      delta_t(j) = 1._mytype
      call compute_tension_term_huang_style(delta_t, x_geom, delta_ft)
      rhs_delta = alpha * delta_ft
      do c = 1, 3
        delta_x(c,:) = 0._mytype
        call solve_free_end_implicit_bending_rhs_scalar(rhs_delta(c,:), delta_x(c,:), alpha, gamma_b, it_used, dxmax, &
             converged, breakdown_detected, restart_count, failure_mode)
        if (.not.converged) then
          schur_matrix_ok = .false.
          schur_matrix_failure_mode = 'bending_basis_solve_failed'
          schur_matrix_failed_column = j
          return
        endif
      enddo
      call lag_ds_center_node_to_half(delta_x, dsdx)
      do i = 1, fiber_nl-1
        smat(i,j) = 2._mytype * dot_product(dsx(:,i), dsdx(:,i))
      enddo
    enddo
  end subroutine build_tension_schur_matrix

  subroutine project_momentum_residual_through_bending(rx, dt_s, rho_tilde, gamma_b, w_out, projection_ok, failure_mode, &
       schur_rx_norm, schur_alpha_rx_norm, schur_w_proj_norm, schur_w_projection_residual_norm, &
       schur_w_projection_relative_residual, schur_w_projection_negligible)
    real(mytype), intent(in) :: rx(3,fiber_nl), dt_s, rho_tilde, gamma_b
    real(mytype), intent(out) :: w_out(3,fiber_nl)
    logical, intent(out) :: projection_ok
    character(len=*), intent(out) :: failure_mode
    real(mytype), intent(out) :: schur_rx_norm, schur_alpha_rx_norm, schur_w_proj_norm
    real(mytype), intent(out) :: schur_w_projection_residual_norm, schur_w_projection_relative_residual
    logical, intent(out) :: schur_w_projection_negligible
    real(mytype) :: alpha, rhs_w(3,fiber_nl), dxmax
    real(mytype) :: projection_rhs_negligible_tol, projection_rel_tol
    integer :: c, it_used, restart_count, failure_code
    logical :: converged, breakdown_detected

    alpha = dt_s * dt_s / rho_tilde
    schur_rx_norm = max_abs_array_2d(rx)
    rhs_w = alpha * rx
    schur_alpha_rx_norm = max_abs_array_2d(rhs_w)
    schur_w_proj_norm = 0._mytype
    schur_w_projection_residual_norm = 0._mytype
    schur_w_projection_relative_residual = 0._mytype
    schur_w_projection_negligible = .false.
    projection_ok = .true.
    failure_mode = 'none'
    internal_update_failure = .false.
    if (schur_rx_norm > 0._mytype .and. schur_alpha_rx_norm <= 0._mytype) then
      projection_ok = .false.
      internal_update_failure = .true.
      failure_mode = 'momentum_projection_rhs_zero_internal_failure'
      return
    endif
    do c = 1, 3
      w_out(c,:) = 0._mytype
      call solve_free_end_implicit_bending_rhs_scalar(rhs_w(c,:), w_out(c,:), alpha, gamma_b, it_used, dxmax, &
           converged, breakdown_detected, restart_count, failure_code)
      if (.not.converged) then
        projection_ok = .false.
        failure_mode = 'momentum_projection_bending_solve_failed'
        return
      endif
    enddo
    schur_w_proj_norm = max_abs_array_2d(w_out)
    projection_rhs_negligible_tol = 1000._mytype * epsilon(1._mytype) * max(1._mytype, fiber_length)
    projection_rel_tol = 1.0e-8_mytype
    call check_implicit_projection_residual(w_out, rhs_w, alpha, gamma_b, projection_rhs_negligible_tol, projection_rel_tol, &
         schur_w_projection_residual_norm, schur_w_projection_relative_residual, schur_w_projection_negligible, projection_ok)
    if (.not.projection_ok) failure_mode = 'momentum_projection_residual_check_failed'
  end subroutine project_momentum_residual_through_bending

  subroutine compute_constraint_linearized_projection(x_geom, vec_node, cvec)
    real(mytype), intent(in) :: x_geom(3,fiber_nl), vec_node(3,fiber_nl)
    real(mytype), intent(out) :: cvec(fiber_nl-1)
    real(mytype) :: dsx(3,fiber_nl-1), dsvec(3,fiber_nl-1)
    integer :: i
    call lag_ds_center_node_to_half(x_geom, dsx)
    call lag_ds_center_node_to_half(vec_node, dsvec)
    do i = 1, fiber_nl-1
      cvec(i) = 2._mytype * dot_product(dsx(:,i), dsvec(:,i))
    enddo
  end subroutine compute_constraint_linearized_projection

  subroutine apply_momentum_consistent_schur_update(x_geom, delta_t, dt_s, rho_tilde, gamma_b, zcorr, update_ok, &
       failure_mode, schur_delta_t_norm, schur_delta_ft_norm, schur_alpha_delta_ft_norm, schur_zcorr_norm, &
       schur_z_projection_residual_norm, schur_z_projection_relative_residual, schur_z_projection_negligible)
    real(mytype), intent(in) :: x_geom(3,fiber_nl), delta_t(fiber_nl-1), dt_s, rho_tilde, gamma_b
    real(mytype), intent(out) :: zcorr(3,fiber_nl)
    logical, intent(out) :: update_ok
    character(len=*), intent(out) :: failure_mode
    real(mytype), intent(out) :: schur_delta_t_norm, schur_delta_ft_norm, schur_alpha_delta_ft_norm, schur_zcorr_norm
    real(mytype), intent(out) :: schur_z_projection_residual_norm, schur_z_projection_relative_residual
    logical, intent(out) :: schur_z_projection_negligible
    real(mytype) :: delta_ft(3,fiber_nl), rhs_z(3,fiber_nl), alpha, dxmax
    real(mytype) :: projection_rhs_negligible_tol, projection_rel_tol
    integer :: c, it_used, restart_count, failure_code
    logical :: converged, breakdown_detected

    alpha = dt_s * dt_s / rho_tilde
    schur_delta_t_norm = maxval(abs(delta_t))
    call compute_tension_term_huang_style(delta_t, x_geom, delta_ft)
    schur_delta_ft_norm = max_abs_array_2d(delta_ft)
    rhs_z = alpha * delta_ft
    schur_alpha_delta_ft_norm = max_abs_array_2d(rhs_z)
    schur_zcorr_norm = 0._mytype
    schur_z_projection_residual_norm = 0._mytype
    schur_z_projection_relative_residual = 0._mytype
    schur_z_projection_negligible = .false.
    update_ok = .true.
    failure_mode = 'none'
    internal_update_failure = .false.
    if (schur_delta_t_norm > 0._mytype .and. schur_delta_ft_norm <= 0._mytype) then
      update_ok = .false.
      internal_update_failure = .true.
      failure_mode = 'tension_response_force_zero_internal_failure'
      return
    endif
    do c = 1, 3
      zcorr(c,:) = 0._mytype
      call solve_free_end_implicit_bending_rhs_scalar(rhs_z(c,:), zcorr(c,:), alpha, gamma_b, it_used, dxmax, &
           converged, breakdown_detected, restart_count, failure_code)
      if (.not.converged) then
        update_ok = .false.
        failure_mode = 'tension_response_bending_solve_failed'
        return
      endif
    enddo
    schur_zcorr_norm = max_abs_array_2d(zcorr)
    projection_rhs_negligible_tol = 1000._mytype * epsilon(1._mytype) * max(1._mytype, fiber_length)
    projection_rel_tol = 1.0e-8_mytype
    call check_implicit_projection_residual(zcorr, rhs_z, alpha, gamma_b, projection_rhs_negligible_tol, projection_rel_tol, &
         schur_z_projection_residual_norm, schur_z_projection_relative_residual, schur_z_projection_negligible, update_ok)
    if (.not.update_ok) failure_mode = 'tension_response_residual_check_failed'
  end subroutine apply_momentum_consistent_schur_update

  subroutine solve_schur_tension_correction(smat, gvec, delta_t, method, reg_used, mu_final, min_abs_pivot, max_abs_pivot, &
       failure_mode)
    real(mytype), intent(in) :: smat(fiber_nl-1,fiber_nl-1), gvec(fiber_nl-1)
    real(mytype), intent(out) :: delta_t(fiber_nl-1), mu_final, min_abs_pivot, max_abs_pivot
    character(len=*), intent(out) :: method, failure_mode
    logical, intent(out) :: reg_used
    real(mytype) :: rhs(fiber_nl-1), dz(fiber_nl-1), sts(fiber_nl-1,fiber_nl-1), sts0(fiber_nl-1,fiber_nl-1), stg(fiber_nl-1), mu
    real(mytype) :: minp, maxp, minp2, maxp2
    integer :: info, i, it

    rhs = -gvec
    call solve_dense_general_pivoted(smat, rhs, dz, fiber_nl-1, info, minp, maxp)
    min_abs_pivot = minp
    max_abs_pivot = maxp
    if (info == 0) then
      delta_t = dz
      method = 'pivoted_lu'
      reg_used = .false.
      mu_final = 0._mytype
      failure_mode = 'none'
      return
    endif
    sts0 = matmul(transpose(smat), smat)
    stg = -matmul(transpose(smat), gvec)
    mu = 1.0e-12_mytype
    reg_used = .true.
    do it = 1, 8
      sts = sts0
      do i = 1, fiber_nl-1
        sts(i,i) = sts(i,i) + mu
      enddo
      call solve_dense_general_pivoted(sts, stg, dz, fiber_nl-1, info, minp2, maxp2)
      if (info == 0) then
        delta_t = dz
        method = 'lm_fallback'
        mu_final = mu
        min_abs_pivot = min(min_abs_pivot, minp2)
        max_abs_pivot = max(max_abs_pivot, maxp2)
        failure_mode = 'lu_near_singular_recovered_by_lm'
        return
      endif
      mu = 10._mytype * mu
    enddo
    delta_t = 0._mytype
    method = 'failed'
    mu_final = mu
    failure_mode = 'lu_and_lm_failed'
  end subroutine solve_schur_tension_correction

  subroutine build_structure_newton_fd_jacobian(uvec, u_scale, xnm1, xn, fext_in, dt_s, rho_tilde, gamma_b, jac, &
       eps_base, hmin, hmax)
    ! dense finite-difference Jacobian is used only for Step 3.5 academic validation;
    ! later stages may replace it with analytic / matrix-free Newton-Krylov.
    real(mytype), intent(in) :: uvec(4*fiber_nl-1), u_scale(4*fiber_nl-1), xnm1(3,fiber_nl), xn(3,fiber_nl), fext_in(3,fiber_nl)
    real(mytype), intent(in) :: dt_s, rho_tilde, gamma_b
    real(mytype), intent(out) :: jac(4*fiber_nl-1,4*fiber_nl-1)
    real(mytype), intent(out) :: eps_base, hmin, hmax
    real(mytype) :: up(4*fiber_nl-1), um(4*fiber_nl-1), rp(4*fiber_nl-1), rm(4*fiber_nl-1), epsj, eps0
    real(mytype) :: rx_scale_dummy, rg_scale_dummy
    integer :: j
    eps0 = sqrt(epsilon(1._mytype))
    eps_base = eps0
    hmin = huge(1._mytype)
    hmax = 0._mytype
    do j = 1, 4*fiber_nl-1
      epsj = eps0 * max(1._mytype, max(abs(uvec(j)), u_scale(j)))
      hmin = min(hmin, epsj)
      hmax = max(hmax, epsj)
      up = uvec
      um = uvec
      up(j) = up(j) + epsj
      um(j) = um(j) - epsj
      call build_structure_newton_scaled_residual_vector(up, xnm1, xn, fext_in, dt_s, rho_tilde, gamma_b, rp, &
           rx_scale_dummy, rg_scale_dummy)
      call build_structure_newton_scaled_residual_vector(um, xnm1, xn, fext_in, dt_s, rho_tilde, gamma_b, rm, &
           rx_scale_dummy, rg_scale_dummy)
      jac(:,j) = (rp - rm) / (2._mytype * epsj)
    enddo
  end subroutine build_structure_newton_fd_jacobian

  subroutine solve_structure_dynamics_newton_step(xnm1, xn, xnp1, t_half, fext_in, dt_s, rho_tilde, gamma_b, &
       newton_converged_strict, newton_converged_effective, newton_iterations_used, newton_final_residual_x, &
       newton_final_residual_g, newton_final_residual_norm, newton_final_update_norm, newton_accepted_lambda, &
       newton_convergence_mode, newton_linear_solver_method, newton_lm_regularization_used, newton_lm_mu_final, &
       newton_min_abs_pivot, newton_max_abs_pivot, newton_linear_failure_mode, fd_jacobian_eps_base, &
       fd_jacobian_min_perturb, fd_jacobian_max_perturb, newton_unknown_scaling_active, newton_update_norm_scaled)
    real(mytype), intent(in) :: xnm1(3,fiber_nl), xn(3,fiber_nl), fext_in(3,fiber_nl), dt_s, rho_tilde, gamma_b
    real(mytype), intent(out) :: xnp1(3,fiber_nl), t_half(fiber_nl-1)
    logical, intent(out) :: newton_converged_strict, newton_converged_effective
    integer, intent(out) :: newton_iterations_used
    real(mytype), intent(out) :: newton_final_residual_x, newton_final_residual_g, newton_final_residual_norm
    real(mytype), intent(out) :: newton_final_update_norm, newton_accepted_lambda
    character(len=*), intent(out) :: newton_convergence_mode
    character(len=*), intent(out) :: newton_linear_solver_method, newton_linear_failure_mode
    logical, intent(out) :: newton_lm_regularization_used, newton_unknown_scaling_active
    real(mytype), intent(out) :: newton_lm_mu_final, newton_min_abs_pivot, newton_max_abs_pivot
    real(mytype), intent(out) :: fd_jacobian_eps_base, fd_jacobian_min_perturb, fd_jacobian_max_perturb
    real(mytype), intent(out) :: newton_update_norm_scaled
    real(mytype) :: u(4*fiber_nl-1), r(4*fiber_nl-1), jmat(4*fiber_nl-1,4*fiber_nl-1), du(4*fiber_nl-1)
    real(mytype) :: utrial(4*fiber_nl-1), rtrial(4*fiber_nl-1), r_unscaled(4*fiber_nl-1), lambda, beta_ls
    real(mytype) :: phi, phi_trial, rx_scale, rg_scale
    real(mytype) :: xk(3,fiber_nl), tk(fiber_nl-1), xtrial(3,fiber_nl), ttrial(fiber_nl-1), xpred(3,fiber_nl)
    real(mytype) :: rx(3,fiber_nl), rg(fiber_nl-1), rx_rel, du_inf, alpha, du_scaled_inf
    real(mytype) :: u_scale(4*fiber_nl-1), x_scale_dummy, t_scale_dummy
    real(mytype) :: phi_hist(3), du_hist(3)
    integer :: it, ibt, maxit, maxbt
    logical :: accepted, plateau, geometry_regular, final_nan_inf, final_geometry_regular

    alpha = dt_s * dt_s / rho_tilde
    xpred = 2._mytype * xn - xnm1 + alpha * fext_in
    if (fiber_flex_constraint_tension_warm_start_active .and. allocated(fiber_tension_half_prevstep) .and. &
         size(fiber_tension_half_prevstep) == fiber_nl-1) then
      tk = fiber_tension_half_prevstep
    else
      tk = 0._mytype
    endif
    xk = xpred
    call pack_structure_unknown(xk, tk, u)

    maxit = fiber_flex_constraint_outer_maxit
    maxbt = fiber_flex_constraint_line_search_max_backtracks
    beta_ls = fiber_flex_constraint_line_search_beta
    phi_hist = huge(1._mytype)
    du_hist = huge(1._mytype)
    newton_convergence_mode = 'failed'
    newton_converged_strict = .false.
    newton_converged_effective = .false.
    newton_accepted_lambda = 0._mytype
    newton_final_update_norm = 0._mytype
    newton_update_norm_scaled = 0._mytype
    newton_iterations_used = 0
    newton_linear_solver_method = 'failed'
    newton_lm_regularization_used = .false.
    newton_lm_mu_final = 0._mytype
    newton_min_abs_pivot = huge(1._mytype)
    newton_max_abs_pivot = 0._mytype
    newton_linear_failure_mode = 'none'
    fd_jacobian_eps_base = sqrt(epsilon(1._mytype))
    fd_jacobian_min_perturb = huge(1._mytype)
    fd_jacobian_max_perturb = 0._mytype
    newton_unknown_scaling_active = .true.

    do it = 1, maxit
      call build_structure_newton_scaled_residual_vector(u, xnm1, xn, fext_in, dt_s, rho_tilde, gamma_b, r, rx_scale, rg_scale, r_unscaled)
      call build_newton_unknown_scaling(u, fext_in, gamma_b, u_scale, x_scale_dummy, t_scale_dummy)
      call unpack_structure_unknown(u, xk, tk)
      call build_coupled_structure_residual(xk, tk, xn, xnm1, rx, rg, fext_in, dt_s, rho_tilde, gamma_b)
      newton_final_residual_x = maxval(abs(rx))
      newton_final_residual_g = maxval(abs(rg))
      rx_rel = newton_final_residual_x / rx_scale
      phi = maxval(abs(r))
      newton_final_residual_norm = phi
      newton_iterations_used = it
      if (newton_final_residual_g <= fiber_flex_constraint_outer_tol_g .and. rx_rel <= fiber_flex_constraint_outer_tol_rx_rel .and. &
           newton_final_update_norm <= fiber_flex_constraint_outer_tol_x) then
        newton_converged_strict = .true.
        newton_converged_effective = .true.
        newton_convergence_mode = 'strict'
        exit
      endif

      call build_structure_newton_fd_jacobian(u, u_scale, xnm1, xn, fext_in, dt_s, rho_tilde, gamma_b, jmat, &
           fd_jacobian_eps_base, fd_jacobian_min_perturb, fd_jacobian_max_perturb)
      call solve_newton_linear_system_with_lm_fallback(jmat, r, u_scale, du, newton_linear_solver_method, &
           newton_lm_regularization_used, newton_lm_mu_final, newton_min_abs_pivot, newton_max_abs_pivot, &
           newton_linear_failure_mode)
      if (trim(newton_linear_solver_method) == 'failed') then
        newton_convergence_mode = 'failed'
        exit
      endif
      du_inf = maxval(abs(du))
      du_scaled_inf = maxval(abs(du / max(u_scale, 1.0e-30_mytype)))
      if (fiber_flex_constraint_line_search_active) then
        lambda = 1._mytype
        accepted = .false.
        do ibt = 0, maxbt
          utrial = u + lambda * du
          call build_structure_newton_scaled_residual_vector(utrial, xnm1, xn, fext_in, dt_s, rho_tilde, gamma_b, rtrial, &
               rx_scale, rg_scale)
          call unpack_structure_unknown(utrial, xtrial, ttrial)
          call build_coupled_structure_residual(xtrial, ttrial, xn, xnm1, rx, rg, fext_in, dt_s, rho_tilde, gamma_b)
          phi_trial = maxval(abs(rtrial))
          if (phi_trial < phi) then
            accepted = .true.
            exit
          endif
          lambda = beta_ls * lambda
        enddo
      else
        lambda = 1._mytype
        utrial = u + du
        accepted = .true.
      endif

      phi_hist = cshift(phi_hist, -1); phi_hist(3) = phi
      du_hist = cshift(du_hist, -1); du_hist(3) = du_inf
      geometry_regular = (all(xk == xk) .and. maxval(abs(xk)) < 1.0e30_mytype)
      plateau = (it >= 3 .and. abs(phi_hist(3)-phi_hist(2)) <= 1.0e-2_mytype * max(phi_hist(2), 1.0e-14_mytype) .and. &
           abs(phi_hist(2)-phi_hist(1)) <= 1.0e-2_mytype * max(phi_hist(1), 1.0e-14_mytype))
      if (.not.accepted) then
        if (newton_final_residual_g <= fiber_flex_constraint_outer_tol_g .and. rx_rel <= fiber_flex_constraint_outer_tol_rx_rel .and. &
             (du_inf <= 10._mytype * fiber_flex_constraint_outer_tol_x .or. plateau) .and. geometry_regular) then
          newton_converged_strict = .false.
          newton_converged_effective = .true.
          newton_convergence_mode = 'plateau_accepted'
          exit
        else
          newton_converged_strict = .false.
          newton_converged_effective = .false.
          newton_convergence_mode = 'failed'
          exit
        endif
      endif

      u = utrial
      newton_accepted_lambda = lambda
      newton_final_update_norm = lambda * du_inf
      newton_update_norm_scaled = lambda * du_scaled_inf
    enddo

    call unpack_structure_unknown(u, xnp1, t_half)
    ! final post-loop residual recheck prevents stale convergence labels after the last accepted Newton step
    call build_structure_newton_scaled_residual_vector(u, xnm1, xn, fext_in, dt_s, rho_tilde, gamma_b, r, rx_scale, rg_scale)
    call build_coupled_structure_residual(xnp1, t_half, xn, xnm1, rx, rg, fext_in, dt_s, rho_tilde, gamma_b)
    newton_final_residual_x = maxval(abs(rx))
    newton_final_residual_g = maxval(abs(rg))
    rx_rel = newton_final_residual_x / max(rx_scale, 1.0e-30_mytype)
    newton_final_residual_norm = maxval(abs(r))
    newton_update_norm_scaled = max(newton_update_norm_scaled, 0._mytype)
    final_nan_inf = (.not.all(xnp1 == xnp1)) .or. (.not.all(t_half == t_half)) .or. (.not.all(r == r))
    final_geometry_regular = (.not.final_nan_inf) .and. maxval(abs(xnp1)) < 1.0e30_mytype
    plateau = (newton_iterations_used >= 3 .and. abs(phi_hist(3)-phi_hist(2)) <= 1.0e-2_mytype * max(phi_hist(2), 1.0e-14_mytype) .and. &
         abs(phi_hist(2)-phi_hist(1)) <= 1.0e-2_mytype * max(phi_hist(1), 1.0e-14_mytype))
    if (newton_final_residual_g <= fiber_flex_constraint_outer_tol_g .and. rx_rel <= fiber_flex_constraint_outer_tol_rx_rel .and. &
         newton_final_update_norm <= fiber_flex_constraint_outer_tol_x .and. final_geometry_regular) then
      newton_convergence_mode = 'strict'
    else if (newton_final_residual_g <= fiber_flex_constraint_outer_tol_g .and. rx_rel <= fiber_flex_constraint_outer_tol_rx_rel .and. &
         (newton_final_update_norm <= 10._mytype * fiber_flex_constraint_outer_tol_x .or. plateau) .and. final_geometry_regular) then
      newton_convergence_mode = 'plateau_accepted'
    else
      newton_convergence_mode = 'failed'
    endif
    if (trim(newton_convergence_mode) == 'strict') then
      newton_converged_strict = .true.
      newton_converged_effective = .true.
    else if (trim(newton_convergence_mode) == 'plateau_accepted') then
      newton_converged_strict = .false.
      newton_converged_effective = .true.
    else
      newton_converged_strict = .false.
      newton_converged_effective = .false.
    endif
    if (newton_converged_effective) then
      if (.not.allocated(fiber_tension_half_prevstep)) allocate(fiber_tension_half_prevstep(fiber_nl-1))
      fiber_tension_half_prevstep = t_half
    endif
  end subroutine solve_structure_dynamics_newton_step

  subroutine solve_structure_dynamics_schur_step(xnm1, xn, fext_in, dt_s, rho_tilde, gamma_b, t_half_initial, xnp1, t_half_out, &
       schur_solver_converged, schur_iterations_used, schur_final_constraint_residual, schur_final_momentum_residual, &
       schur_final_rx_rel, schur_final_update_norm_x, schur_final_update_norm_t, schur_accepted_lambda, &
       schur_convergence_mode, schur_linear_solver_method, schur_regularization_used, schur_regularization_mu_final, &
       schur_min_abs_pivot, schur_max_abs_pivot, schur_linear_failure_mode, schur_position_solve_converged, &
       position_solve_failure_mode, schur_matrix_ok, schur_matrix_failure_mode, schur_matrix_failed_column, &
       schur_rhs_momentum_projection_active, schur_rhs_momentum_projection_ok, schur_rhs_failure_mode, &
       schur_momentum_consistent_rhs_active, schur_combined_merit_initial, schur_combined_merit_final, schur_delta_x_norm, &
       schur_delta_t_norm, schur_nonzero_update_accepted, schur_pre_state_effective_ok, schur_rx_norm, &
       schur_alpha_rx_norm, schur_w_proj_norm, schur_delta_ft_norm, schur_alpha_delta_ft_norm, schur_zcorr_norm, &
       schur_x_trial_update_norm, schur_internal_update_failure, schur_w_projection_residual_norm, &
       schur_w_projection_relative_residual, schur_w_projection_negligible, schur_z_projection_residual_norm, &
       schur_z_projection_relative_residual, schur_z_projection_negligible, schur_zero_update_warning, &
       schur_projection_failure_mode, schur_tension_response_ok)
    real(mytype), intent(in) :: xnm1(3,fiber_nl), xn(3,fiber_nl), fext_in(3,fiber_nl), dt_s, rho_tilde, gamma_b
    real(mytype), intent(in) :: t_half_initial(fiber_nl-1)
    real(mytype), intent(out) :: xnp1(3,fiber_nl), t_half_out(fiber_nl-1)
    logical, intent(out) :: schur_solver_converged, schur_regularization_used
    integer, intent(out) :: schur_iterations_used
    real(mytype), intent(out) :: schur_final_constraint_residual, schur_final_momentum_residual, schur_final_rx_rel
    real(mytype), intent(out) :: schur_final_update_norm_x, schur_final_update_norm_t, schur_accepted_lambda
    real(mytype), intent(out) :: schur_regularization_mu_final, schur_min_abs_pivot, schur_max_abs_pivot
    character(len=*), intent(out) :: schur_convergence_mode, schur_linear_solver_method, schur_linear_failure_mode
    logical, intent(out) :: schur_position_solve_converged
    character(len=*), intent(out) :: position_solve_failure_mode
    logical, intent(out) :: schur_matrix_ok
    character(len=*), intent(out) :: schur_matrix_failure_mode
    integer, intent(out) :: schur_matrix_failed_column
    logical, intent(out) :: schur_rhs_momentum_projection_active, schur_rhs_momentum_projection_ok
    logical, intent(out) :: schur_momentum_consistent_rhs_active, schur_nonzero_update_accepted, schur_pre_state_effective_ok
    logical, intent(out) :: schur_internal_update_failure
    logical, intent(out) :: schur_w_projection_negligible, schur_z_projection_negligible, schur_zero_update_warning
    logical, intent(out) :: schur_tension_response_ok
    character(len=*), intent(out) :: schur_rhs_failure_mode
    character(len=*), intent(out) :: schur_projection_failure_mode
    real(mytype), intent(out) :: schur_combined_merit_initial, schur_combined_merit_final, schur_delta_x_norm, schur_delta_t_norm
    real(mytype), intent(out) :: schur_rx_norm, schur_alpha_rx_norm, schur_w_proj_norm
    real(mytype), intent(out) :: schur_delta_ft_norm, schur_alpha_delta_ft_norm, schur_zcorr_norm, schur_x_trial_update_norm
    real(mytype), intent(out) :: schur_w_projection_residual_norm, schur_w_projection_relative_residual
    real(mytype), intent(out) :: schur_z_projection_residual_norm, schur_z_projection_relative_residual
    real(mytype) :: xk(3,fiber_nl), xtrial(3,fiber_nl), tk(fiber_nl-1), ttrial(fiber_nl-1), delta_t(fiber_nl-1)
    real(mytype) :: gk(fiber_nl-1), gtrial(fiber_nl-1), smat(fiber_nl-1,fiber_nl-1), lambda, phi_current, phi_trial
    real(mytype) :: max_abs_g, max_seg_err, max_inext_err, max_abs_g_trial, rx_scale, alpha, rx_rel_trial
    real(mytype) :: rhs_t(fiber_nl-1), rhs_solver(fiber_nl-1), cw(fiber_nl-1), wproj(3,fiber_nl), zcorr(3,fiber_nl), delta_x(3,fiber_nl)
    real(mytype) :: rx(3,fiber_nl), rg(fiber_nl-1), beta_ls
    real(mytype) :: delta_x_floor
    integer :: it, ibt, maxit, maxbt
    logical :: accepted, solve_ok, geometry_sane, geometry_finite, plateau_detected, strict_success, plateau_success
    logical :: update_ok, pre_state_effective_ok, accepted_nonzero_update
    logical :: rx_ok_current, rx_ok_trial, rx_ok_pre, rx_ok_plateau
    real(mytype) :: plateau_g_tol, plateau_rx_rel_tol, g_scale, rx_scale_success, tiny_merit, update_floor, update_norm, phi_prev
    real(mytype) :: rx_tol_guard

    alpha = dt_s * dt_s / rho_tilde
    xk = 2._mytype * xn - xnm1 + alpha * fext_in
    tk = t_half_initial
    maxit = fiber_flex_constraint_outer_maxit
    maxbt = fiber_flex_constraint_line_search_max_backtracks
    beta_ls = fiber_flex_constraint_line_search_beta
    plateau_g_tol = 10._mytype * fiber_flex_constraint_outer_tol_g
    plateau_rx_rel_tol = 10._mytype * fiber_flex_constraint_outer_tol_rx_rel
    schur_solver_converged = .false.
    schur_convergence_mode = 'failed'
    schur_accepted_lambda = 0._mytype
    schur_final_update_norm_x = 0._mytype
    schur_final_update_norm_t = 0._mytype
    schur_iterations_used = 0
    schur_position_solve_converged = .true.
    position_solve_failure_mode = 'none'
    schur_matrix_ok = .true.
    schur_matrix_failure_mode = 'none'
    schur_matrix_failed_column = 0
    schur_linear_solver_method = 'failed'
    schur_linear_failure_mode = 'none'
    schur_regularization_used = .false.
    schur_regularization_mu_final = 0._mytype
    schur_min_abs_pivot = huge(1._mytype)
    schur_max_abs_pivot = 0._mytype
    schur_rhs_momentum_projection_active = .true.
    schur_rhs_momentum_projection_ok = .true.
    schur_rhs_failure_mode = 'none'
    schur_projection_failure_mode = 'none'
    schur_momentum_consistent_rhs_active = .true.
    schur_combined_merit_initial = huge(1._mytype)
    schur_combined_merit_final = huge(1._mytype)
    schur_delta_x_norm = 0._mytype
    schur_delta_t_norm = 0._mytype
    schur_rx_norm = 0._mytype
    schur_alpha_rx_norm = 0._mytype
    schur_w_proj_norm = 0._mytype
    schur_delta_ft_norm = 0._mytype
    schur_alpha_delta_ft_norm = 0._mytype
    schur_zcorr_norm = 0._mytype
    schur_x_trial_update_norm = 0._mytype
    schur_w_projection_residual_norm = 0._mytype
    schur_w_projection_relative_residual = 0._mytype
    schur_w_projection_negligible = .false.
    schur_z_projection_residual_norm = 0._mytype
    schur_z_projection_relative_residual = 0._mytype
    schur_z_projection_negligible = .false.
    schur_zero_update_warning = .false.
    schur_tension_response_ok = .true.
    schur_internal_update_failure = .false.
    schur_nonzero_update_accepted = .false.
    schur_pre_state_effective_ok = .false.
    tiny_merit = 1.0e-30_mytype
    g_scale = max(fiber_flex_constraint_outer_tol_g, tiny_merit)
    rx_scale_success = max(fiber_flex_constraint_outer_tol_rx_rel, tiny_merit)
    update_floor = 10._mytype * epsilon(1._mytype)
    delta_x_floor = 100._mytype * epsilon(1._mytype) * max(1._mytype, fiber_length)
    rx_tol_guard = 1._mytype + 1.0e-5_mytype
    phi_prev = huge(1._mytype)

    do it = 1, maxit
      solve_ok = .true.
      call compute_half_segment_constraint_residual(xk, gk, max_abs_g, max_seg_err, max_inext_err)
      call build_coupled_structure_residual(xk, tk, xn, xnm1, rx, rg, fext_in, dt_s, rho_tilde, gamma_b)
      call compute_momentum_residual_scale_final_time(xk, tk, xn, xnm1, fext_in, dt_s, rho_tilde, gamma_b, rx_scale)
      schur_final_momentum_residual = maxval(abs(rx))
      schur_final_rx_rel = schur_final_momentum_residual / max(rx_scale, 1.0e-30_mytype)
      rx_ok_current = schur_final_rx_rel <= fiber_flex_constraint_outer_tol_rx_rel * rx_tol_guard
      geometry_finite = all(xk == xk) .and. all(tk == tk) .and. all(rx == rx) .and. all(gk == gk)
      geometry_sane = max_seg_err <= 1.0e-2_mytype .and. max_inext_err <= 1.0e-2_mytype
      schur_iterations_used = it
      strict_success = solve_ok .and. geometry_finite .and. geometry_sane .and. &
           (max_abs_g <= fiber_flex_constraint_outer_tol_g) .and. &
           rx_ok_current
      phi_current = max(max_abs_g / g_scale, schur_final_rx_rel / rx_scale_success)
      if (it == 1) schur_combined_merit_initial = phi_current
      schur_combined_merit_final = phi_current
      if (strict_success) then
        schur_solver_converged = .true.
        schur_convergence_mode = 'strict'
        exit
      endif
      call project_momentum_residual_through_bending(rx, dt_s, rho_tilde, gamma_b, wproj, update_ok, schur_rhs_failure_mode, &
           schur_rx_norm, schur_alpha_rx_norm, schur_w_proj_norm, schur_w_projection_residual_norm, &
           schur_w_projection_relative_residual, schur_w_projection_negligible)
      schur_rhs_momentum_projection_ok = update_ok
      if (.not.update_ok) then
        schur_position_solve_converged = .false.
        position_solve_failure_mode = 'momentum_projection_failed'
        schur_projection_failure_mode = 'w_projection_residual_failed'
        schur_convergence_mode = 'failed'
        exit
      endif
      call compute_constraint_linearized_projection(xk, wproj, cw)
      rhs_t = -gk + cw
      rhs_solver = -rhs_t
      call build_tension_schur_matrix(xk, dt_s, rho_tilde, gamma_b, smat, schur_matrix_ok, schur_matrix_failure_mode, &
           schur_matrix_failed_column)
      if (.not.schur_matrix_ok) then
        schur_solver_converged = .false.
        schur_convergence_mode = 'failed'
        schur_linear_solver_method = 'failed'
        schur_linear_failure_mode = trim(schur_matrix_failure_mode)
        exit
      endif
      call solve_schur_tension_correction(smat, rhs_solver, delta_t, schur_linear_solver_method, schur_regularization_used, &
           schur_regularization_mu_final, schur_min_abs_pivot, schur_max_abs_pivot, schur_linear_failure_mode)
      if (trim(schur_linear_solver_method) == 'failed') then
        schur_convergence_mode = 'failed'
        exit
      endif
      call apply_momentum_consistent_schur_update(xk, delta_t, dt_s, rho_tilde, gamma_b, zcorr, update_ok, &
           schur_rhs_failure_mode, schur_delta_t_norm, schur_delta_ft_norm, schur_alpha_delta_ft_norm, schur_zcorr_norm, &
           schur_z_projection_residual_norm, schur_z_projection_relative_residual, schur_z_projection_negligible)
      schur_tension_response_ok = update_ok
      schur_rhs_momentum_projection_ok = schur_rhs_momentum_projection_ok .and. update_ok
      if (.not.update_ok) then
        schur_position_solve_converged = .false.
        position_solve_failure_mode = 'tension_response_failed'
        schur_projection_failure_mode = 'z_projection_residual_failed'
        schur_convergence_mode = 'failed'
        exit
      endif
      delta_x = -wproj + zcorr
      schur_delta_x_norm = max_abs_array_2d(delta_x)
      schur_zero_update_warning = schur_zero_update_warning .or. (schur_delta_x_norm <= delta_x_floor)
      lambda = 1._mytype
      accepted = .false.
      do ibt = 0, maxbt
        ttrial = tk + lambda * delta_t
        xtrial = xk + lambda * delta_x
        schur_x_trial_update_norm = max_abs_array_2d(xtrial - xk)
        call compute_half_segment_constraint_residual(xtrial, gtrial, max_abs_g_trial, max_seg_err, max_inext_err)
        call build_coupled_structure_residual(xtrial, ttrial, xn, xnm1, rx, rg, fext_in, dt_s, rho_tilde, gamma_b)
        call compute_momentum_residual_scale_final_time(xtrial, ttrial, xn, xnm1, fext_in, dt_s, rho_tilde, gamma_b, rx_scale)
        rx_rel_trial = maxval(abs(rx)) / max(rx_scale, 1.0e-30_mytype)
        rx_ok_trial = rx_rel_trial <= fiber_flex_constraint_outer_tol_rx_rel * rx_tol_guard
        geometry_finite = all(xtrial == xtrial) .and. all(ttrial == ttrial) .and. all(rx == rx) .and. all(gtrial == gtrial)
        geometry_sane = max_seg_err <= 1.0e-2_mytype .and. max_inext_err <= 1.0e-2_mytype
        phi_trial = max(max_abs_g_trial / g_scale, rx_rel_trial / rx_scale_success)
        if ((phi_trial < phi_current .or. (max_abs_g_trial <= fiber_flex_constraint_outer_tol_g .and. &
             rx_ok_trial)) .and. geometry_finite .and. geometry_sane) then
          accepted = .true.
          exit
        endif
        lambda = beta_ls * lambda
      enddo
      if (trim(schur_convergence_mode) == 'failed' .and. schur_internal_update_failure) exit
      update_norm = lambda * max(schur_delta_x_norm, schur_delta_t_norm)
      accepted_nonzero_update = accepted .and. (lambda > 0._mytype) .and. (schur_delta_x_norm > delta_x_floor) .and. &
           (update_norm > update_floor) .and. (.not.schur_internal_update_failure)
      schur_nonzero_update_accepted = schur_nonzero_update_accepted .or. accepted_nonzero_update
      rx_ok_pre = schur_final_rx_rel <= plateau_rx_rel_tol * rx_tol_guard
      pre_state_effective_ok = solve_ok .and. geometry_finite .and. geometry_sane .and. &
           (max_abs_g <= plateau_g_tol) .and. rx_ok_pre
      schur_pre_state_effective_ok = schur_pre_state_effective_ok .or. pre_state_effective_ok
      plateau_detected = (it >= 2) .and. (abs(phi_current - phi_prev) <= 1.0e-2_mytype * max(phi_prev, tiny_merit))
      if (.not.accepted) then
        strict_success = solve_ok .and. geometry_finite .and. geometry_sane .and. &
             (max_abs_g <= fiber_flex_constraint_outer_tol_g) .and. &
             rx_ok_current
        rx_ok_plateau = schur_final_rx_rel <= plateau_rx_rel_tol * rx_tol_guard
        plateau_success = (.not.schur_internal_update_failure) .and. solve_ok .and. geometry_finite .and. geometry_sane .and. &
             (max_abs_g <= plateau_g_tol) .and. rx_ok_plateau .and. &
             ( (pre_state_effective_ok .and. plateau_detected) .or. accepted_nonzero_update )
        if (strict_success) then
          schur_solver_converged = .true.
          schur_convergence_mode = 'strict'
        else if (plateau_success) then
          schur_solver_converged = .true.
          schur_convergence_mode = 'plateau_accepted'
        else
          schur_solver_converged = .false.
          schur_convergence_mode = 'failed'
        endif
        exit
      endif
      schur_accepted_lambda = lambda
      schur_final_update_norm_t = lambda * maxval(abs(delta_t))
      schur_final_update_norm_x = lambda * maxval(abs(delta_x))
      tk = ttrial
      xk = xtrial
      call compute_half_segment_constraint_residual(xk, gk, max_abs_g, max_seg_err, max_inext_err)
      call build_coupled_structure_residual(xk, tk, xn, xnm1, rx, rg, fext_in, dt_s, rho_tilde, gamma_b)
      call compute_momentum_residual_scale_final_time(xk, tk, xn, xnm1, fext_in, dt_s, rho_tilde, gamma_b, rx_scale)
      schur_final_momentum_residual = maxval(abs(rx))
      schur_final_constraint_residual = maxval(abs(gk))
      schur_final_rx_rel = schur_final_momentum_residual / max(rx_scale, 1.0e-30_mytype)
      schur_combined_merit_final = max(schur_final_constraint_residual / g_scale, schur_final_rx_rel / rx_scale_success)
      phi_prev = phi_current
    enddo

    xnp1 = xk
    t_half_out = tk
    call build_coupled_structure_residual(xnp1, t_half_out, xn, xnm1, rx, rg, fext_in, dt_s, rho_tilde, gamma_b)
    call compute_momentum_residual_scale_final_time(xnp1, t_half_out, xn, xnm1, fext_in, dt_s, rho_tilde, gamma_b, rx_scale)
    schur_final_momentum_residual = maxval(abs(rx))
    schur_final_constraint_residual = maxval(abs(rg))
    schur_final_rx_rel = schur_final_momentum_residual / max(rx_scale, 1.0e-30_mytype)
    schur_combined_merit_final = max(schur_final_constraint_residual / g_scale, schur_final_rx_rel / rx_scale_success)
  end subroutine solve_structure_dynamics_schur_step

  subroutine advance_structure_dynamics_with_retry(case_id, t_macro_start, xnm1, xn, xnp1, t_half, dt_macro, rho_tilde, gamma_b, &
       step_success, nsub_used, rejection_count, retry_count, min_substep_dt, schur_solver_converged, schur_iterations_used, &
       schur_final_constraint_residual, schur_final_momentum_residual, schur_final_rx_rel, schur_final_update_norm_x, &
       schur_final_update_norm_t, schur_accepted_lambda, schur_convergence_mode, schur_linear_solver_method, &
       schur_regularization_used, schur_regularization_mu_final, schur_min_abs_pivot, schur_max_abs_pivot, &
       schur_linear_failure_mode, schur_position_solve_converged, position_solve_failure_mode, bending_kernel_failure_count, &
       schur_matrix_ok, schur_matrix_failure_mode, schur_matrix_failed_column, structure_local_tension_warm_start_active, &
       schur_rhs_momentum_projection_active, schur_rhs_momentum_projection_ok, schur_rhs_failure_mode, &
       schur_momentum_consistent_rhs_active, schur_combined_merit_initial, schur_combined_merit_final, schur_delta_x_norm, &
       schur_delta_t_norm, schur_nonzero_update_accepted, schur_pre_state_effective_ok, schur_rx_norm, &
       schur_alpha_rx_norm, schur_w_proj_norm, schur_delta_ft_norm, schur_alpha_delta_ft_norm, schur_zcorr_norm, &
       schur_x_trial_update_norm, schur_internal_update_failure, schur_w_projection_residual_norm, &
       schur_w_projection_relative_residual, schur_w_projection_negligible, schur_z_projection_residual_norm, &
       schur_z_projection_relative_residual, schur_z_projection_negligible, schur_zero_update_warning, &
       schur_projection_failure_mode, schur_tension_response_ok)
    integer, intent(in) :: case_id
    real(mytype), intent(in) :: t_macro_start
    real(mytype), intent(in) :: xnm1(3,fiber_nl), xn(3,fiber_nl), dt_macro, rho_tilde, gamma_b
    real(mytype), intent(out) :: xnp1(3,fiber_nl), t_half(fiber_nl-1), min_substep_dt
    logical, intent(out) :: step_success, schur_solver_converged, schur_regularization_used
    integer, intent(out) :: nsub_used, rejection_count, retry_count, schur_iterations_used
    real(mytype), intent(out) :: schur_final_constraint_residual, schur_final_momentum_residual, schur_final_rx_rel
    real(mytype), intent(out) :: schur_final_update_norm_x, schur_final_update_norm_t, schur_accepted_lambda
    real(mytype), intent(out) :: schur_regularization_mu_final, schur_min_abs_pivot, schur_max_abs_pivot
    character(len=*), intent(out) :: schur_convergence_mode, schur_linear_solver_method, schur_linear_failure_mode
    logical, intent(out) :: schur_position_solve_converged
    character(len=*), intent(out) :: position_solve_failure_mode
    integer, intent(out) :: bending_kernel_failure_count
    logical, intent(out) :: schur_matrix_ok, structure_local_tension_warm_start_active
    character(len=*), intent(out) :: schur_matrix_failure_mode
    integer, intent(out) :: schur_matrix_failed_column
    logical, intent(out) :: schur_rhs_momentum_projection_active, schur_rhs_momentum_projection_ok
    logical, intent(out) :: schur_momentum_consistent_rhs_active, schur_nonzero_update_accepted, schur_pre_state_effective_ok
    logical, intent(out) :: schur_internal_update_failure
    logical, intent(out) :: schur_w_projection_negligible, schur_z_projection_negligible, schur_zero_update_warning
    logical, intent(out) :: schur_tension_response_ok
    character(len=*), intent(out) :: schur_rhs_failure_mode
    character(len=*), intent(out) :: schur_projection_failure_mode
    real(mytype), intent(out) :: schur_combined_merit_initial, schur_combined_merit_final, schur_delta_x_norm, schur_delta_t_norm
    real(mytype), intent(out) :: schur_rx_norm, schur_alpha_rx_norm, schur_w_proj_norm
    real(mytype), intent(out) :: schur_delta_ft_norm, schur_alpha_delta_ft_norm, schur_zcorr_norm, schur_x_trial_update_norm
    real(mytype), intent(out) :: schur_w_projection_residual_norm, schur_w_projection_relative_residual
    real(mytype), intent(out) :: schur_z_projection_residual_norm, schur_z_projection_relative_residual
    real(mytype) :: xprev(3,fiber_nl), xcur(3,fiber_nl), xnext(3,fiber_nl), veln(3,fiber_nl), dt_sub
    real(mytype) :: fext_sub(3,fiber_nl), t_sub
    real(mytype) :: t_prev_local(fiber_nl-1), t_sub_out(fiber_nl-1), tension_prev_backup(fiber_nl-1)
    integer :: split, nsub, isub
    logical :: ok_local

    step_success = .false.
    rejection_count = 0
    retry_count = 0
    min_substep_dt = dt_macro
    bending_kernel_failure_count = 0
    schur_matrix_ok = .true.
    schur_matrix_failure_mode = 'none'
    schur_matrix_failed_column = 0
    schur_rhs_momentum_projection_active = .true.
    schur_rhs_momentum_projection_ok = .true.
    schur_rhs_failure_mode = 'none'
    schur_projection_failure_mode = 'none'
    schur_momentum_consistent_rhs_active = .true.
    schur_combined_merit_initial = huge(1._mytype)
    schur_combined_merit_final = huge(1._mytype)
    schur_delta_x_norm = 0._mytype
    schur_delta_t_norm = 0._mytype
    schur_rx_norm = 0._mytype
    schur_alpha_rx_norm = 0._mytype
    schur_w_proj_norm = 0._mytype
    schur_delta_ft_norm = 0._mytype
    schur_alpha_delta_ft_norm = 0._mytype
    schur_zcorr_norm = 0._mytype
    schur_x_trial_update_norm = 0._mytype
    schur_w_projection_residual_norm = 0._mytype
    schur_w_projection_relative_residual = 0._mytype
    schur_w_projection_negligible = .false.
    schur_z_projection_residual_norm = 0._mytype
    schur_z_projection_relative_residual = 0._mytype
    schur_z_projection_negligible = .false.
    schur_zero_update_warning = .false.
    schur_tension_response_ok = .true.
    schur_internal_update_failure = .false.
    schur_nonzero_update_accepted = .false.
    schur_pre_state_effective_ok = .false.
    structure_local_tension_warm_start_active = .true.
    if (.not.allocated(fiber_tension_half_prevstep)) allocate(fiber_tension_half_prevstep(fiber_nl-1))
    tension_prev_backup = fiber_tension_half_prevstep
    do split = 0, merge(fiber_flex_structure_max_substep_splits, 0, fiber_flex_structure_adaptive_substep_active)
      nsub = 2**split
      dt_sub = dt_macro / real(nsub, mytype)
      min_substep_dt = min(min_substep_dt, dt_sub)
      veln = (xn - xnm1) / dt_macro
      xprev = xn - dt_sub * veln
      xcur = xn
      t_prev_local = tension_prev_backup
      ok_local = .true.
      do isub = 1, nsub
        t_sub = t_macro_start + real(isub,mytype) * dt_sub
        call compute_prescribed_structure_force(case_id, t_sub, xcur, fext_sub)
        call solve_structure_dynamics_schur_step(xprev, xcur, fext_sub, dt_sub, rho_tilde, gamma_b, t_prev_local, xnext, t_sub_out, &
             schur_solver_converged, schur_iterations_used, schur_final_constraint_residual, schur_final_momentum_residual, &
             schur_final_rx_rel, schur_final_update_norm_x, schur_final_update_norm_t, schur_accepted_lambda, &
             schur_convergence_mode, schur_linear_solver_method, schur_regularization_used, schur_regularization_mu_final, &
             schur_min_abs_pivot, schur_max_abs_pivot, schur_linear_failure_mode, schur_position_solve_converged, &
             position_solve_failure_mode, schur_matrix_ok, schur_matrix_failure_mode, schur_matrix_failed_column, &
             schur_rhs_momentum_projection_active, schur_rhs_momentum_projection_ok, schur_rhs_failure_mode, &
             schur_momentum_consistent_rhs_active, schur_combined_merit_initial, schur_combined_merit_final, schur_delta_x_norm, &
             schur_delta_t_norm, schur_nonzero_update_accepted, schur_pre_state_effective_ok, schur_rx_norm, &
             schur_alpha_rx_norm, schur_w_proj_norm, schur_delta_ft_norm, schur_alpha_delta_ft_norm, schur_zcorr_norm, &
             schur_x_trial_update_norm, schur_internal_update_failure, schur_w_projection_residual_norm, &
             schur_w_projection_relative_residual, schur_w_projection_negligible, schur_z_projection_residual_norm, &
             schur_z_projection_relative_residual, schur_z_projection_negligible, schur_zero_update_warning, &
             schur_projection_failure_mode, schur_tension_response_ok)
        if (.not.schur_solver_converged) then
          if (.not.schur_position_solve_converged) bending_kernel_failure_count = bending_kernel_failure_count + 1
          ok_local = .false.
          exit
        endif
        t_prev_local = t_sub_out
        xprev = xcur
        xcur = xnext
      enddo
      if (ok_local) then
        step_success = .true.
        xnp1 = xcur
        t_half = t_prev_local
        nsub_used = nsub
        fiber_tension_half_prevstep = t_prev_local
        return
      endif
      rejection_count = rejection_count + 1
      retry_count = retry_count + nsub
    enddo
    fiber_tension_half_prevstep = tension_prev_backup
    t_half = tension_prev_backup
    xnp1 = xn
    nsub_used = 2**4
  end subroutine advance_structure_dynamics_with_retry

  subroutine solve_coupled_constrained_structure_step(xnm1, xn, xnp1, t_half, fext_in, dt_c, rho_tilde, gamma_b, &
       max_abs_drift_term, max_abs_vel_term, max_abs_force_term, implicit_bending_update_norm, &
       max_outer_iterations_used, final_coupled_residual_x, final_coupled_residual_g, coupled_solver_converged, &
       coupled_solver_converged_strict, coupled_solver_converged_effective, &
       max_abs_constraint_residual_during_outer, max_abs_momentum_residual_during_outer, &
       max_abs_delta_x_during_outer, max_abs_delta_t_during_outer, accepted_line_search_lambda, &
       final_coupled_residual_x_rel, final_coupled_convergence_mode, plateau_detected_final)
    ! Final-time coupled structure-only constrained implicit solve (Step 3.4):
    ! solves X^{n+1} and T^{n+1/2} consistently at the same time layer.
    real(mytype), intent(in) :: xnm1(3,fiber_nl), xn(3,fiber_nl), fext_in(3,fiber_nl), dt_c, rho_tilde, gamma_b
    real(mytype), intent(out) :: xnp1(3,fiber_nl), t_half(fiber_nl-1)
    real(mytype), intent(out) :: max_abs_drift_term, max_abs_vel_term, max_abs_force_term, implicit_bending_update_norm
    integer, intent(out) :: max_outer_iterations_used
    real(mytype), intent(out) :: final_coupled_residual_x, final_coupled_residual_g
    logical, intent(out) :: coupled_solver_converged, coupled_solver_converged_strict, coupled_solver_converged_effective
    real(mytype), intent(out) :: max_abs_constraint_residual_during_outer, max_abs_momentum_residual_during_outer
    real(mytype), intent(out) :: max_abs_delta_x_during_outer, max_abs_delta_t_during_outer, accepted_line_search_lambda
    real(mytype), intent(out) :: final_coupled_residual_x_rel
    character(len=*), intent(out) :: final_coupled_convergence_mode
    logical, intent(out) :: plateau_detected_final
    real(mytype) :: xk(3,fiber_nl), xnew(3,fiber_nl), rx(3,fiber_nl), rg(fiber_nl-1), rhs_corr(fiber_nl-1)
    real(mytype) :: amat_corr(fiber_nl-1,fiber_nl-1), delta_t(fiber_nl-1), ftension(3,fiber_nl), rhs_eff(3,fiber_nl)
    real(mytype) :: alpha, dxmax, dtmax, lambda_ls, phi_current, phi_trial, max_abs_rx, max_abs_rg
    real(mytype) :: rx_scale, rx_rel, prev_dxmax, rel_improve1, rel_improve2
    real(mytype) :: x_base(3,fiber_nl), t_base(fiber_nl-1), t_trial(fiber_nl-1), x_trial(3,fiber_nl), delta_x(3,fiber_nl)
    integer :: coupled_max_outer, ls_bt, ls_bt_max
    logical :: line_search_active, accepted_step, has_prevstep_tension, plateau_detected, line_search_limited, geometry_regular
    real(mytype) :: rx_rel_hist(3), dx_hist(3), lambda_hist(3), phi_hist(3)
    real(mytype) :: rhs_tmp(fiber_nl-1), drift_tmp(fiber_nl-1), vel_tmp(fiber_nl-1), force_tmp(fiber_nl-1)
    real(mytype) :: fb_tmp(3,fiber_nl)
    integer :: it, c, it_used, restart_count, failure_mode
    logical :: converged, breakdown_detected

    xk = 2._mytype * xn - xnm1
    if (fiber_flex_constraint_tension_warm_start_active .and. allocated(fiber_tension_half_prevstep)) then
      has_prevstep_tension = (size(fiber_tension_half_prevstep) == fiber_nl-1)
    else
      has_prevstep_tension = .false.
    endif
    if (has_prevstep_tension) then
      ! half-grid tension is warm-started from the previous time step to improve robustness of the coupled final-time solve
      t_half = fiber_tension_half_prevstep
    else
      t_half = 0._mytype
    endif
    alpha = dt_c * dt_c / rho_tilde
    coupled_max_outer = fiber_flex_constraint_outer_maxit
    line_search_active = fiber_flex_constraint_line_search_active
    ls_bt_max = fiber_flex_constraint_line_search_max_backtracks
    coupled_solver_converged = .false.
    coupled_solver_converged_strict = .false.
    coupled_solver_converged_effective = .false.
    max_abs_constraint_residual_during_outer = 0._mytype
    max_abs_momentum_residual_during_outer = 0._mytype
    max_abs_delta_x_during_outer = 0._mytype
    max_abs_delta_t_during_outer = 0._mytype
    accepted_line_search_lambda = 0._mytype
    final_coupled_residual_x_rel = huge(1._mytype)
    final_coupled_convergence_mode = 'failed'
    plateau_detected_final = .false.
    rx_rel_hist = huge(1._mytype)
    dx_hist = huge(1._mytype)
    lambda_hist = 1._mytype
    phi_hist = huge(1._mytype)
    prev_dxmax = huge(1._mytype)
    implicit_bending_update_norm = 0._mytype
    max_outer_iterations_used = 0

    do it = 1, coupled_max_outer
      call build_coupled_structure_residual(xk, t_half, xn, xnm1, rx, rg, fext_in, dt_c, rho_tilde, gamma_b)
      final_coupled_residual_x = maxval(abs(rx))
      final_coupled_residual_g = maxval(abs(rg))
      max_abs_rx = final_coupled_residual_x
      max_abs_rg = final_coupled_residual_g
      phi_current = max(max_abs_rx, max_abs_rg)
      call compute_momentum_residual_scale_final_time(xk, t_half, xn, xnm1, fext_in, dt_c, rho_tilde, gamma_b, rx_scale)
      rx_rel = max_abs_rx / rx_scale
      final_coupled_residual_x_rel = rx_rel
      max_abs_momentum_residual_during_outer = max(max_abs_momentum_residual_during_outer, max_abs_rx)
      max_abs_constraint_residual_during_outer = max(max_abs_constraint_residual_during_outer, max_abs_rg)
      max_outer_iterations_used = it
      rx_rel_hist = cshift(rx_rel_hist, -1)
      rx_rel_hist(3) = rx_rel
      dx_hist = cshift(dx_hist, -1)
      dx_hist(3) = prev_dxmax
      lambda_hist = cshift(lambda_hist, -1)
      lambda_hist(3) = accepted_line_search_lambda
      phi_hist = cshift(phi_hist, -1)
      phi_hist(3) = phi_current
      rel_improve1 = abs(rx_rel_hist(3) - rx_rel_hist(2)) / max(rx_rel_hist(2), 1.0e-14_mytype)
      rel_improve2 = abs(rx_rel_hist(2) - rx_rel_hist(1)) / max(rx_rel_hist(1), 1.0e-14_mytype)
      line_search_limited = (lambda_hist(2) < 1._mytype .and. lambda_hist(3) < 1._mytype)
      geometry_regular = (all(xk == xk) .and. maxval(abs(xk)) < 1.0e30_mytype)
      plateau_detected = (it >= 3 .and. rel_improve1 <= 1.0e-2_mytype .and. rel_improve2 <= 1.0e-2_mytype .and. &
           abs(phi_hist(3)-phi_hist(2)) <= 1.0e-2_mytype * max(phi_hist(2), 1.0e-14_mytype) .and. &
           (dx_hist(3) <= 10._mytype * fiber_flex_constraint_outer_tol_x .or. line_search_limited))
      plateau_detected_final = plateau_detected
      ! convergence is declared when inextensibility is satisfied and the scaled momentum residual is sufficiently small;
      ! exact machine-zero absolute momentum residual is not required for the structure-only coupled DAE solve
      if (max_abs_rg <= fiber_flex_constraint_outer_tol_g .and. rx_rel <= fiber_flex_constraint_outer_tol_rx_rel .and. &
           prev_dxmax <= fiber_flex_constraint_outer_tol_x) then
        coupled_solver_converged = .true.
        coupled_solver_converged_strict = .true.
        coupled_solver_converged_effective = .true.
        final_coupled_convergence_mode = 'strict'
        exit
      endif
      if (max_abs_rg <= fiber_flex_constraint_outer_tol_g .and. rx_rel <= fiber_flex_constraint_outer_tol_rx_rel .and. &
           plateau_detected .and. geometry_regular) then
        coupled_solver_converged = .true.
        coupled_solver_converged_strict = .false.
        coupled_solver_converged_effective = .true.
        final_coupled_convergence_mode = 'plateau_accepted'
        exit
      endif

      rhs_corr = -2._mytype/(dt_c*dt_c) * rg
      call assemble_tension_matrix_huang_style(xk, amat_corr)
      call solve_dense_small(amat_corr, rhs_corr, delta_t, fiber_nl-1)
      dtmax = maxval(abs(delta_t))
      max_abs_delta_t_during_outer = max(max_abs_delta_t_during_outer, dtmax)
      t_base = t_half
      t_trial = t_base + delta_t

      call compute_tension_term_huang_style(t_trial, xk, ftension)
      rhs_eff = 2._mytype * xn - xnm1 + alpha * (ftension + fext_in)
      do c = 1, 3
        xnew(c,:) = xk(c,:)
        call solve_free_end_implicit_bending_rhs_scalar(rhs_eff(c,:), xnew(c,:), alpha, gamma_b, it_used, dxmax, &
             converged, breakdown_detected, restart_count, failure_mode)
        if (.not.converged) then
          if (nrank == 0) write(*,*) 'Error: coupled X-update implicit solve failed in Step 3.4.'
          stop
        endif
      enddo
      x_base = xk
      delta_x = xnew - x_base
      dxmax = maxval(abs(delta_x))

      if (.not.line_search_active) then
        lambda_ls = 1._mytype
        accepted_step = .true.
      else
        lambda_ls = 1._mytype
        accepted_step = .false.
        do ls_bt = 0, ls_bt_max
          x_trial = x_base + lambda_ls * delta_x
          t_half = t_base + lambda_ls * delta_t
          call build_coupled_structure_residual(x_trial, t_half, xn, xnm1, rx, rg, fext_in, dt_c, rho_tilde, gamma_b)
          phi_trial = max(maxval(abs(rx)), maxval(abs(rg)))
          if (phi_trial < phi_current) then
            accepted_step = .true.
            exit
          endif
          lambda_ls = fiber_flex_constraint_line_search_beta * lambda_ls
        enddo
      endif

      if (.not.accepted_step) then
        ! Line-search exhaustion does not automatically imply failure:
        ! first check whether we are already on an acceptable coupled plateau.
        geometry_regular = (all(xk == xk) .and. maxval(abs(xk)) < 1.0e30_mytype)
        if (max_abs_rg <= fiber_flex_constraint_outer_tol_g .and. rx_rel <= fiber_flex_constraint_outer_tol_rx_rel .and. &
             (prev_dxmax <= 10._mytype * fiber_flex_constraint_outer_tol_x .or. plateau_detected) .and. geometry_regular) then
          coupled_solver_converged = .true.
          coupled_solver_converged_strict = .false.
          coupled_solver_converged_effective = .true.
          final_coupled_residual_x = max_abs_rx
          final_coupled_residual_g = max_abs_rg
          final_coupled_residual_x_rel = rx_rel
          final_coupled_convergence_mode = 'plateau_accepted'
          plateau_detected_final = .true.
          exit
        endif
        coupled_solver_converged = .false.
        coupled_solver_converged_strict = .false.
        coupled_solver_converged_effective = .false.
        final_coupled_residual_x = max_abs_rx
        final_coupled_residual_g = max_abs_rg
        final_coupled_residual_x_rel = rx_rel
        final_coupled_convergence_mode = 'failed'
        if (nrank == 0) write(*,*) 'Error: coupled line search exhausted without acceptable plateau in Step 3.4.'
        exit
      endif

      accepted_line_search_lambda = lambda_ls
      xk = x_base + lambda_ls * delta_x
      t_half = t_base + lambda_ls * delta_t
      dxmax = maxval(abs(lambda_ls * delta_x))
      dtmax = maxval(abs(lambda_ls * delta_t))
      max_abs_delta_x_during_outer = max(max_abs_delta_x_during_outer, dxmax)
      max_abs_delta_t_during_outer = max(max_abs_delta_t_during_outer, dtmax)
      prev_dxmax = dxmax
      implicit_bending_update_norm = max(implicit_bending_update_norm, dxmax)
    enddo

    if (trim(final_coupled_convergence_mode) == 'strict') then
      coupled_solver_converged_strict = .true.
      coupled_solver_converged_effective = .true.
    else if (trim(final_coupled_convergence_mode) == 'plateau_accepted') then
      coupled_solver_converged_strict = .false.
      coupled_solver_converged_effective = .true.
    else
      coupled_solver_converged_strict = .false.
      coupled_solver_converged_effective = .false.
    endif
    coupled_solver_converged = coupled_solver_converged_effective
    xnp1 = xk
    if (coupled_solver_converged) then
      if (.not.allocated(fiber_tension_half_prevstep)) allocate(fiber_tension_half_prevstep(fiber_nl-1))
      fiber_tension_half_prevstep = t_half
    endif
    do c = 1, 3
      call apply_free_end_bending_operator_scalar(xk(c,:), fb_tmp(c,:), gamma_b)
    enddo
    call build_tension_rhs_huang_style(xn, xnm1, xk, fb_tmp, fext_in, dt_c, rhs_tmp, drift_tmp, vel_tmp, force_tmp)
    max_abs_drift_term = maxval(abs(drift_tmp))
    max_abs_vel_term = maxval(abs(vel_tmp))
    max_abs_force_term = maxval(abs(force_tmp))
  end subroutine solve_coupled_constrained_structure_step

  subroutine advance_constrained_structure_step(xnm1, xn, xnp1, t_half, fext_n, dt_c, &
       max_abs_drift_term, max_abs_vel_term, max_abs_force_term, implicit_bending_update_norm, &
       post_bending_correction_norm, post_bending_max_inext_err_after_correction, max_outer_iterations_used, &
       final_coupled_residual_x, final_coupled_residual_g, coupled_solver_converged, &
       coupled_solver_converged_strict, coupled_solver_converged_effective, &
       max_abs_constraint_residual_during_outer, max_abs_momentum_residual_during_outer, &
       max_abs_delta_x_during_outer, max_abs_delta_t_during_outer, accepted_line_search_lambda, &
       final_coupled_residual_x_rel, final_coupled_convergence_mode, plateau_detected)
    real(mytype), intent(in) :: xnm1(3,fiber_nl), xn(3,fiber_nl), fext_n(3,fiber_nl), dt_c
    real(mytype), intent(out) :: xnp1(3,fiber_nl), t_half(fiber_nl-1)
    real(mytype), intent(out) :: max_abs_drift_term, max_abs_vel_term, max_abs_force_term, implicit_bending_update_norm
    real(mytype), intent(out) :: post_bending_correction_norm, post_bending_max_inext_err_after_correction
    integer, intent(out) :: max_outer_iterations_used
    real(mytype), intent(out) :: final_coupled_residual_x, final_coupled_residual_g
    logical, intent(out) :: coupled_solver_converged, coupled_solver_converged_strict, coupled_solver_converged_effective
    real(mytype), intent(out) :: max_abs_constraint_residual_during_outer, max_abs_momentum_residual_during_outer
    real(mytype), intent(out) :: max_abs_delta_x_during_outer, max_abs_delta_t_during_outer, accepted_line_search_lambda
    real(mytype), intent(out) :: final_coupled_residual_x_rel
    character(len=*), intent(out) :: final_coupled_convergence_mode
    logical, intent(out) :: plateau_detected

    ! Default path now uses final-time coupled solve.
    ! Deprecated 3.4 predictor-corrector / post-correction reference path is intentionally isolated.
    call solve_coupled_constrained_structure_step(xnm1, xn, xnp1, t_half, fext_n, dt_c, fiber_structure_rho_tilde, &
         fiber_bending_gamma, max_abs_drift_term, max_abs_vel_term, max_abs_force_term, implicit_bending_update_norm, &
         max_outer_iterations_used, final_coupled_residual_x, final_coupled_residual_g, coupled_solver_converged, &
         coupled_solver_converged_strict, coupled_solver_converged_effective, &
         max_abs_constraint_residual_during_outer, max_abs_momentum_residual_during_outer, &
         max_abs_delta_x_during_outer, max_abs_delta_t_during_outer, accepted_line_search_lambda, &
         final_coupled_residual_x_rel, final_coupled_convergence_mode, plateau_detected)
    post_bending_correction_norm = 0._mytype
    post_bending_max_inext_err_after_correction = 0._mytype
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

  subroutine compute_prescribed_structure_force(case_id, time_now, x_state, fext_out)
    ! prescribed structure-only forcing, not from IBM
    integer, intent(in) :: case_id
    real(mytype), intent(in) :: time_now, x_state(3,fiber_nl)
    real(mytype), intent(out) :: fext_out(3,fiber_nl)
    real(mytype) :: dir(3), norm_dir, q, gshape, amp_now
    integer :: l

    fext_out = 0._mytype
    dir = fiber_flex_structure_force_direction
    norm_dir = sqrt(sum(dir**2))
    if (norm_dir > 0._mytype) dir = dir / norm_dir
    select case (fiber_flex_structure_force_mode)
    case (0)
      amp_now = 0._mytype
    case (1)
      amp_now = fiber_flex_structure_force_amp * sin(fiber_flex_structure_force_omega * time_now)
    case (2)
      amp_now = fiber_flex_structure_force_amp
    case default
      amp_now = 0._mytype
    end select
    if (case_id == 3) amp_now = 5._mytype * amp_now
    do l = 1, fiber_nl
      q = (fiber_s_ref(l) - fiber_s_ref(1)) / max(fiber_length, 1.0e-12_mytype)
      gshape = q**4 * (1._mytype - q)**4
      fext_out(:,l) = amp_now * gshape * dir
    enddo
  end subroutine compute_prescribed_structure_force

  subroutine reparameterize_curve_to_uniform_arclength(x_in, x_out)
    ! Deprecated uniform-on-current-curve-length reparameterization.
    ! Keeps equal spacing along the *current* curve length and does not
    ! enforce projection to target fiber_length.
    real(mytype), intent(in) :: x_in(3,fiber_nl)
    real(mytype), intent(out) :: x_out(3,fiber_nl)
    real(mytype) :: seglen(fiber_nl-1), sacc(fiber_nl), starget, total_len, w
    integer :: i, j

    sacc(1) = 0._mytype
    do i = 1, fiber_nl-1
      seglen(i) = sqrt(sum((x_in(:,i+1)-x_in(:,i))**2))
      sacc(i+1) = sacc(i) + seglen(i)
    enddo
    total_len = max(sacc(fiber_nl), 1.0e-14_mytype)
    x_out(:,1) = x_in(:,1)
    x_out(:,fiber_nl) = x_in(:,fiber_nl)
    do i = 2, fiber_nl-1
      starget = (real(i-1,mytype)/real(fiber_nl-1,mytype)) * total_len
      j = 1
      do while (j < fiber_nl-1 .and. sacc(j+1) < starget)
        j = j + 1
      enddo
      if (sacc(j+1) > sacc(j)) then
        w = (starget - sacc(j)) / (sacc(j+1) - sacc(j))
      else
        w = 0._mytype
      endif
      x_out(:,i) = (1._mytype - w) * x_in(:,j) + w * x_in(:,j+1)
    enddo
  end subroutine reparameterize_curve_to_uniform_arclength

  subroutine build_local_frame_from_tangent(et, en)
    real(mytype), intent(in) :: et(3)
    real(mytype), intent(out) :: en(3)
    real(mytype) :: eref(3), proj

    eref = (/0._mytype, 1._mytype, 0._mytype/)
    proj = dot_product(eref, et)
    en = eref - proj * et
    if (sqrt(sum(en**2)) < 1.0e-10_mytype) then
      eref = (/0._mytype, 0._mytype, 1._mytype/)
      proj = dot_product(eref, et)
      en = eref - proj * et
    endif
    if (sqrt(sum(en**2)) < 1.0e-10_mytype) then
      eref = (/1._mytype, 0._mytype, 0._mytype/)
      proj = dot_product(eref, et)
      en = eref - proj * et
    endif
    if (sqrt(sum(en**2)) < 1.0e-12_mytype) then
      if (nrank == 0) write(*,*) 'Error: cannot build transverse basis for inextensible initial projection.'
      stop
    endif
    en = en / sqrt(sum(en**2))
  end subroutine build_local_frame_from_tangent

  subroutine project_initial_curve_to_target_arclength(x_shape, x_out)
    ! initial amplitude controls the tangent-angle perturbation used to build an inextensible initial curve
    real(mytype), intent(in) :: x_shape(3,fiber_nl)
    real(mytype), intent(out) :: x_out(3,fiber_nl)
    real(mytype) :: et(3), en(3), delta(3,fiber_nl), yshape(fiber_nl), ys(fiber_nl), theta(fiber_nl)
    real(mytype) :: ds_target, s_mid, center_now(3), shift(3)
    integer :: i

    et = fiber_direction
    if (sqrt(sum(et**2)) <= 1.0e-12_mytype) then
      if (nrank == 0) write(*,*) 'Error: fiber_direction norm is zero in initial inextensible projection.'
      stop
    endif
    et = et / sqrt(sum(et**2))
    call build_local_frame_from_tangent(et, en)

    delta = x_shape - fiber_x
    do i = 1, fiber_nl
      yshape(i) = dot_product(delta(:,i), en)
    enddo
    ys(1) = (yshape(2) - yshape(1)) / fiber_ds
    do i = 2, fiber_nl-1
      ys(i) = (yshape(i+1) - yshape(i-1)) / (2._mytype * fiber_ds)
    enddo
    ys(fiber_nl) = (yshape(fiber_nl) - yshape(fiber_nl-1)) / fiber_ds
    theta = atan(ys)

    ds_target = fiber_length / real(max(1,fiber_nl-1), mytype)
    x_out(:,1) = 0._mytype
    do i = 1, fiber_nl-1
      s_mid = 0.5_mytype * (theta(i) + theta(i+1))
      x_out(:,i+1) = x_out(:,i) + ds_target * (cos(s_mid) * et + sin(s_mid) * en)
    enddo
    center_now = 0._mytype
    do i = 1, fiber_nl
      center_now = center_now + x_out(:,i)
    enddo
    center_now = center_now / real(fiber_nl, mytype)
    shift = fiber_center - center_now
    do i = 1, fiber_nl
      x_out(:,i) = x_out(:,i) + shift
    enddo
  end subroutine project_initial_curve_to_target_arclength

  subroutine compute_free_end_boundary_residuals_independent(x_in, max_d2_independent, max_d3_independent)
    ! Independent one-sided physical-node boundary estimator for free-end d2/d3 residuals.
    ! Used for diagnostics only; does not alter solver operators.
    real(mytype), intent(in) :: x_in(3,fiber_nl)
    real(mytype), intent(out) :: max_d2_independent, max_d3_independent
    real(mytype) :: d2_left(3), d2_right(3), d3_left(3), d3_right(3), invds2, invds3

    invds2 = 1._mytype / (fiber_ds * fiber_ds)
    invds3 = 1._mytype / (fiber_ds * fiber_ds * fiber_ds)
    d2_left = (2._mytype*x_in(:,1) - 5._mytype*x_in(:,2) + 4._mytype*x_in(:,3) - x_in(:,4)) * invds2
    d2_right = (2._mytype*x_in(:,fiber_nl) - 5._mytype*x_in(:,fiber_nl-1) + 4._mytype*x_in(:,fiber_nl-2) - &
         x_in(:,fiber_nl-3)) * invds2
    d3_left = (-5._mytype*x_in(:,1) + 18._mytype*x_in(:,2) - 24._mytype*x_in(:,3) + 14._mytype*x_in(:,4) - &
         3._mytype*x_in(:,5)) * (0.5_mytype*invds3)
    d3_right = (5._mytype*x_in(:,fiber_nl) - 18._mytype*x_in(:,fiber_nl-1) + 24._mytype*x_in(:,fiber_nl-2) - &
         14._mytype*x_in(:,fiber_nl-3) + 3._mytype*x_in(:,fiber_nl-4)) * (0.5_mytype*invds3)
    max_d2_independent = max(maxval(abs(d2_left)), maxval(abs(d2_right)))
    max_d3_independent = max(maxval(abs(d3_left)), maxval(abs(d3_right)))
  end subroutine compute_free_end_boundary_residuals_independent

  subroutine initialize_structure_dynamics_test_state(case_id, xnm1, xn, initial_shape_projection_active, &
       initial_max_seg_err, initial_max_inext_err, initial_projection_semantics)
    integer, intent(in) :: case_id
    real(mytype), intent(out) :: xnm1(3,fiber_nl), xn(3,fiber_nl)
    logical, intent(out) :: initial_shape_projection_active
    real(mytype), intent(out) :: initial_max_seg_err, initial_max_inext_err
    character(len=*), intent(out) :: initial_projection_semantics
    real(mytype) :: q, gshape, amp0
    real(mytype) :: xproj(3,fiber_nl)
    integer :: l

    xn = fiber_x
    amp0 = fiber_flex_structure_initial_shape_amp
    initial_shape_projection_active = .false.
    initial_projection_semantics = 'none'
    select case (case_id)
    case (0)
      do l = 1, fiber_nl
        q = (fiber_s_ref(l) - fiber_s_ref(1)) / max(fiber_length, 1.0e-12_mytype)
        gshape = q**4 * (1._mytype - q)**4
        xn(2,l) = xn(2,l) + amp0 * gshape
      enddo
      call project_initial_curve_to_target_arclength(xn, xproj)
      xn = xproj
      initial_shape_projection_active = .true.
      initial_projection_semantics = 'target_length_tangent_angle_inextensible_projection'
      xnm1 = xn
    case (1)
      xnm1 = xn
    case (2)
      do l = 1, fiber_nl
        q = (fiber_s_ref(l) - fiber_s_ref(1)) / max(fiber_length, 1.0e-12_mytype)
        gshape = q**4 * (1._mytype - q)**4
        xn(2,l) = xn(2,l) + 1.5_mytype * amp0 * gshape
      enddo
      call project_initial_curve_to_target_arclength(xn, xproj)
      xn = xproj
      initial_shape_projection_active = .true.
      initial_projection_semantics = 'target_length_tangent_angle_inextensible_projection'
      xnm1 = xn
    case (3)
      ! optional stronger harmonic-load validation case; keeps prescribed-force semantics
      xnm1 = xn
      if (fiber_flex_structure_force_mode == 0) then
        if (nrank == 0) write(*,*) 'Warning: structure case 3 is most informative with non-zero force mode.'
      endif
    case default
      if (nrank == 0) write(*,*) 'Error: fiber_flex_structure_case must be 0,1,2,3.'
      stop
    end select
    call compute_constraint_error(xn, initial_max_seg_err, initial_max_inext_err)
  end subroutine initialize_structure_dynamics_test_state

  subroutine compute_structure_dynamics_metrics(xn, xnp1, fext_now, dt_s, kinetic_energy, bending_energy, external_power, &
       max_end_bc_d2_residual, max_end_bc_d3_residual, max_end_bc_d2_residual_independent, &
       max_end_bc_d3_residual_independent)
    real(mytype), intent(in) :: xn(3,fiber_nl), xnp1(3,fiber_nl), fext_now(3,fiber_nl), dt_s
    real(mytype), intent(out) :: kinetic_energy, bending_energy, external_power, max_end_bc_d2_residual, max_end_bc_d3_residual
    real(mytype), intent(out) :: max_end_bc_d2_residual_independent, max_end_bc_d3_residual_independent
    real(mytype) :: u(3,fiber_nl), xss(3,fiber_nl), d3_left(3), d3_right(3)
    integer :: c
    u = (xnp1 - xn) / dt_s
    kinetic_energy = 0.5_mytype * fiber_structure_rho_tilde * fiber_ds * sum(u**2)
    do c = 1, 3
      call apply_free_end_d2_scalar(xnp1(c,:), xss(c,:))
    enddo
    bending_energy = 0.5_mytype * fiber_bending_gamma * fiber_ds * sum(xss**2)
    external_power = fiber_ds * sum(fext_now * u)
    max_end_bc_d2_residual = max(maxval(abs(xss(:,1))), maxval(abs(xss(:,fiber_nl))))
    d3_left = (xss(:,2) - xss(:,1)) / fiber_ds
    d3_right = (xss(:,fiber_nl) - xss(:,fiber_nl-1)) / fiber_ds
    max_end_bc_d3_residual = max(maxval(abs(d3_left)), maxval(abs(d3_right)))
    if (fiber_nl >= 5) then
      call compute_free_end_boundary_residuals_independent(xnp1, max_end_bc_d2_residual_independent, max_end_bc_d3_residual_independent)
    else
      max_end_bc_d2_residual_independent = max_end_bc_d2_residual
      max_end_bc_d3_residual_independent = max_end_bc_d3_residual
    endif
  end subroutine compute_structure_dynamics_metrics

  subroutine run_flexible_structure_dynamics_test()
    real(mytype), allocatable :: xnm1(:,:), xn(:,:), xnp1(:,:), xinit(:,:), fext(:,:), t_half(:)
    real(mytype) :: tnow, max_seg_err, max_inext_err, max_seg_err_global, max_inext_err_global
    real(mytype) :: max_abs_drift_term, max_abs_vel_term, max_abs_force_term, implicit_bending_update_norm
    real(mytype) :: post_bending_correction_norm, post_bending_max_inext_err_after_correction
    real(mytype) :: final_coupled_residual_x, final_coupled_residual_g, final_coupled_residual_x_rel
    real(mytype) :: max_abs_constraint_residual_during_outer, max_abs_momentum_residual_during_outer
    real(mytype) :: max_abs_delta_x_during_outer, max_abs_delta_t_during_outer, accepted_line_search_lambda
    real(mytype) :: kinetic_energy, bending_energy, external_power
    real(mytype) :: max_abs_fext, max_abs_fext_global, max_abs_fext_final, external_power_final
    real(mytype) :: max_end_bc_d2_residual, max_end_bc_d3_residual
    real(mytype) :: max_end_bc_d2_residual_independent, max_end_bc_d3_residual_independent
    real(mytype) :: kinetic_energy_final, bending_energy_final, final_displacement_norm
    real(mytype) :: max_end_bc_d2_residual_global, max_end_bc_d3_residual_global
    real(mytype) :: max_end_bc_d2_residual_independent_global, max_end_bc_d3_residual_independent_global
    real(mytype) :: max_forced_predictor_norm_global, max_newton_accepted_update_norm_global, rx_scale_step
    real(mytype) :: initial_max_seg_err, initial_max_inext_err
    logical :: initial_shape_projection_active
    character(len=96) :: initial_projection_semantics
    integer :: step, outint, max_outer_iterations_used
    logical :: coupled_solver_converged, coupled_solver_converged_strict, coupled_solver_converged_effective
    logical :: plateau_detected
    character(len=32) :: final_coupled_convergence_mode
    real(mytype) :: newton_final_residual_norm, newton_final_update_norm, newton_accepted_lambda
    integer :: newton_iterations_used
    logical :: step_success, structure_step_failed_any, structure_adaptive_substepping_active
    integer :: structure_step_rejection_count, structure_substep_retry_count, structure_max_substeps_used, nsub_used
    real(mytype) :: structure_min_substep_dt, max_schur_accepted_update_norm_global
    real(mytype) :: structure_failed_time
    integer :: structure_failed_step_index
    real(mytype) :: schur_final_constraint_residual, schur_final_momentum_residual, schur_final_rx_rel
    real(mytype) :: schur_final_update_norm_x, schur_final_update_norm_t, schur_accepted_lambda
    real(mytype) :: schur_regularization_mu_final, schur_min_abs_pivot, schur_max_abs_pivot
    logical :: schur_solver_converged, schur_regularization_used
    integer :: schur_iterations_used
    character(len=32) :: schur_convergence_mode, schur_linear_solver_method, schur_linear_failure_mode
    logical :: schur_position_solve_converged
    character(len=32) :: position_solve_failure_mode
    integer :: bending_kernel_failure_count
    logical :: schur_matrix_ok, structure_local_tension_warm_start_active
    character(len=64) :: schur_matrix_failure_mode
    integer :: schur_matrix_failed_column
    logical :: schur_rhs_momentum_projection_active, schur_rhs_momentum_projection_ok
    logical :: schur_momentum_consistent_rhs_active, schur_nonzero_update_accepted, schur_pre_state_effective_ok
    logical :: schur_internal_update_failure
    logical :: schur_w_projection_negligible, schur_z_projection_negligible, schur_zero_update_warning
    logical :: schur_tension_response_ok
    character(len=64) :: schur_rhs_failure_mode
    character(len=64) :: schur_projection_failure_mode
    real(mytype) :: schur_combined_merit_initial, schur_combined_merit_final, schur_delta_x_norm, schur_delta_t_norm
    real(mytype) :: schur_rx_norm, schur_alpha_rx_norm, schur_w_proj_norm
    real(mytype) :: schur_delta_ft_norm, schur_alpha_delta_ft_norm, schur_zcorr_norm, schur_x_trial_update_norm
    real(mytype) :: schur_w_projection_residual_norm, schur_w_projection_relative_residual
    real(mytype) :: schur_z_projection_residual_norm, schur_z_projection_relative_residual

    if (.not.fiber_active .or. .not.fiber_flexible_active .or. .not.fiber_flex_initialized .or. &
         .not.fiber_flex_structure_test_active) then
      if (nrank == 0) write(*,*) 'Error: invalid state for flexible structure dynamics test.'
      stop
    endif

    allocate(xnm1(3,fiber_nl), xn(3,fiber_nl), xnp1(3,fiber_nl), xinit(3,fiber_nl), fext(3,fiber_nl), t_half(fiber_nl-1))
    if (allocated(fiber_tension_half_prevstep)) fiber_tension_half_prevstep = 0._mytype
    call initialize_structure_dynamics_test_state(fiber_flex_structure_case, xnm1, xn, initial_shape_projection_active, &
         initial_max_seg_err, initial_max_inext_err, initial_projection_semantics)
    xinit = xn
    outint = max(1, fiber_flex_structure_output_interval)
    max_seg_err_global = 0._mytype
    max_inext_err_global = 0._mytype
    max_end_bc_d2_residual_global = 0._mytype
    max_end_bc_d3_residual_global = 0._mytype
    max_end_bc_d2_residual_independent_global = 0._mytype
    max_end_bc_d3_residual_independent_global = 0._mytype
    max_forced_predictor_norm_global = 0._mytype
    max_newton_accepted_update_norm_global = 0._mytype
    coupled_solver_converged_strict = .false.
    coupled_solver_converged_effective = .false.
    kinetic_energy_final = 0._mytype
    bending_energy_final = 0._mytype
    final_coupled_residual_x = 0._mytype
    final_coupled_residual_g = 0._mytype
    final_coupled_residual_x_rel = 0._mytype
    final_coupled_convergence_mode = 'failed'
    max_abs_fext_global = 0._mytype
    max_abs_fext_final = 0._mytype
    external_power_final = 0._mytype
    structure_adaptive_substepping_active = .true.
    structure_step_rejection_count = 0
    structure_substep_retry_count = 0
    structure_max_substeps_used = 1
    structure_min_substep_dt = fiber_flex_structure_dt
    structure_step_failed_any = .false.
    structure_failed_step_index = 0
    structure_failed_time = 0._mytype
    max_schur_accepted_update_norm_global = 0._mytype

    call write_fiber_flex_structure_series(0, 0._mytype, 0._mytype, 0._mytype, 0._mytype, 0._mytype, 0._mytype, &
         0._mytype, 0._mytype, 0._mytype, 0._mytype, 0, 0._mytype, 0._mytype, 0._mytype, 'failed', .true.)

    do step = 1, fiber_flex_structure_nsteps
      tnow = real(step,mytype) * fiber_flex_structure_dt
      call compute_prescribed_structure_force(fiber_flex_structure_case, tnow, xn, fext)
      max_abs_fext = maxval(abs(fext))
      max_abs_fext_global = max(max_abs_fext_global, max_abs_fext)
      max_abs_fext_final = max_abs_fext
      max_forced_predictor_norm_global = max(max_forced_predictor_norm_global, &
           (fiber_flex_structure_dt*fiber_flex_structure_dt/max(fiber_structure_rho_tilde,1.0e-12_mytype))*max_abs_fext)
      call advance_structure_dynamics_with_retry(fiber_flex_structure_case, tnow - fiber_flex_structure_dt, xnm1, xn, xnp1, t_half, &
           fiber_flex_structure_dt, fiber_structure_rho_tilde, &
           fiber_bending_gamma, step_success, nsub_used, structure_step_rejection_count, structure_substep_retry_count, &
           structure_min_substep_dt, schur_solver_converged, schur_iterations_used, schur_final_constraint_residual, &
           schur_final_momentum_residual, schur_final_rx_rel, schur_final_update_norm_x, schur_final_update_norm_t, &
           schur_accepted_lambda, schur_convergence_mode, schur_linear_solver_method, schur_regularization_used, &
           schur_regularization_mu_final, schur_min_abs_pivot, schur_max_abs_pivot, schur_linear_failure_mode, &
           schur_position_solve_converged, position_solve_failure_mode, bending_kernel_failure_count, schur_matrix_ok, &
           schur_matrix_failure_mode, schur_matrix_failed_column, structure_local_tension_warm_start_active, &
           schur_rhs_momentum_projection_active, schur_rhs_momentum_projection_ok, schur_rhs_failure_mode, &
           schur_momentum_consistent_rhs_active, schur_combined_merit_initial, schur_combined_merit_final, schur_delta_x_norm, &
           schur_delta_t_norm, schur_nonzero_update_accepted, schur_pre_state_effective_ok, schur_rx_norm, &
           schur_alpha_rx_norm, schur_w_proj_norm, schur_delta_ft_norm, schur_alpha_delta_ft_norm, schur_zcorr_norm, &
           schur_x_trial_update_norm, schur_internal_update_failure, schur_w_projection_residual_norm, &
           schur_w_projection_relative_residual, schur_w_projection_negligible, schur_z_projection_residual_norm, &
           schur_z_projection_relative_residual, schur_z_projection_negligible, schur_zero_update_warning, &
           schur_projection_failure_mode, schur_tension_response_ok)
      structure_max_substeps_used = max(structure_max_substeps_used, nsub_used)
      if (.not.step_success) then
        structure_step_failed_any = .true.
        structure_failed_step_index = step
        structure_failed_time = tnow
      endif
      max_schur_accepted_update_norm_global = max(max_schur_accepted_update_norm_global, schur_final_update_norm_x)
      max_newton_accepted_update_norm_global = max(max_newton_accepted_update_norm_global, schur_final_update_norm_x)
      coupled_solver_converged = step_success
      coupled_solver_converged_effective = step_success
      coupled_solver_converged_strict = step_success .and. trim(schur_convergence_mode) == 'strict'
      final_coupled_convergence_mode = schur_convergence_mode
      final_coupled_residual_x = schur_final_momentum_residual
      final_coupled_residual_g = schur_final_constraint_residual
      final_coupled_residual_x_rel = schur_final_rx_rel
      newton_iterations_used = schur_iterations_used
      newton_final_residual_norm = max(schur_final_constraint_residual, schur_final_rx_rel)
      newton_final_update_norm = schur_final_update_norm_x
      newton_accepted_lambda = schur_accepted_lambda
      call compute_momentum_residual_scale_final_time(xnp1, t_half, xn, xnm1, fext, fiber_flex_structure_dt, &
           fiber_structure_rho_tilde, fiber_bending_gamma, rx_scale_step)
      max_outer_iterations_used = newton_iterations_used
      accepted_line_search_lambda = newton_accepted_lambda
      max_abs_momentum_residual_during_outer = final_coupled_residual_x
      max_abs_constraint_residual_during_outer = final_coupled_residual_g
      call compute_constraint_error(xnp1, max_seg_err, max_inext_err)
      call compute_structure_dynamics_metrics(xn, xnp1, fext, fiber_flex_structure_dt, kinetic_energy, bending_energy, external_power, &
           max_end_bc_d2_residual, max_end_bc_d3_residual, max_end_bc_d2_residual_independent, &
           max_end_bc_d3_residual_independent)
      max_seg_err_global = max(max_seg_err_global, max_seg_err)
      max_inext_err_global = max(max_inext_err_global, max_inext_err)
      max_end_bc_d2_residual_global = max(max_end_bc_d2_residual_global, max_end_bc_d2_residual)
      max_end_bc_d3_residual_global = max(max_end_bc_d3_residual_global, max_end_bc_d3_residual)
      max_end_bc_d2_residual_independent_global = max(max_end_bc_d2_residual_independent_global, max_end_bc_d2_residual_independent)
      max_end_bc_d3_residual_independent_global = max(max_end_bc_d3_residual_independent_global, max_end_bc_d3_residual_independent)
      kinetic_energy_final = kinetic_energy
      bending_energy_final = bending_energy
      external_power_final = external_power
      if (mod(step,outint) == 0 .or. step == fiber_flex_structure_nsteps) then
        call write_fiber_flex_structure_series(step, tnow, kinetic_energy, bending_energy, external_power, max_abs_fext, &
             max_seg_err, max_inext_err, max_abs_momentum_residual_during_outer, max_abs_constraint_residual_during_outer, &
             final_coupled_residual_x_rel, max_outer_iterations_used, accepted_line_search_lambda, &
             newton_final_update_norm, newton_final_residual_norm, final_coupled_convergence_mode, .false.)
      endif
      if (.not.step_success) exit
      xnm1 = xn
      xn = xnp1
    enddo

    final_displacement_norm = sqrt(sum((xn - xinit)**2) / real(3*fiber_nl,mytype))
    if (trim(final_coupled_convergence_mode) == 'strict') then
      coupled_solver_converged_strict = .true.
      coupled_solver_converged_effective = .true.
    else if (trim(final_coupled_convergence_mode) == 'plateau_accepted') then
      coupled_solver_converged_strict = .false.
      coupled_solver_converged_effective = .true.
    else
      coupled_solver_converged_strict = .false.
      coupled_solver_converged_effective = .false.
    endif
    call write_fiber_flex_structure_summary(fiber_flex_structure_case, fiber_flex_structure_dt, fiber_flex_structure_nsteps, &
         coupled_solver_converged_strict, coupled_solver_converged_effective, final_coupled_convergence_mode, &
         final_coupled_residual_x, final_coupled_residual_x_rel, final_coupled_residual_g, kinetic_energy_final, &
         bending_energy_final, max_seg_err_global, max_inext_err_global, max_end_bc_d2_residual_global, &
         max_end_bc_d3_residual_global, max_end_bc_d2_residual_independent_global, max_end_bc_d3_residual_independent_global, &
         final_displacement_norm, max_abs_fext_global, max_abs_fext_final, external_power_final, &
         max_forced_predictor_norm_global, max_newton_accepted_update_norm_global, initial_shape_projection_active, &
         initial_max_seg_err, initial_max_inext_err, trim(initial_projection_semantics), &
         newton_iterations_used, newton_final_residual_norm, &
         newton_final_update_norm, newton_accepted_lambda, schur_linear_solver_method, &
         schur_regularization_used, schur_regularization_mu_final, schur_min_abs_pivot, schur_max_abs_pivot, &
         schur_linear_failure_mode, sqrt(epsilon(1._mytype)), 0._mytype, 0._mytype, &
         .true., schur_final_update_norm_x, 'schur_tension_implicit_bending', .true., schur_solver_converged, &
         schur_iterations_used, schur_final_constraint_residual, schur_final_momentum_residual, schur_final_rx_rel, &
         schur_accepted_lambda, schur_convergence_mode, structure_adaptive_substepping_active, &
         structure_step_rejection_count, structure_substep_retry_count, structure_max_substeps_used, structure_min_substep_dt, &
         structure_step_failed_any, max_schur_accepted_update_norm_global, schur_position_solve_converged, &
         trim(position_solve_failure_mode), bending_kernel_failure_count, .true., schur_matrix_ok, &
         trim(schur_matrix_failure_mode), schur_matrix_failed_column, structure_failed_step_index, structure_failed_time, &
         structure_local_tension_warm_start_active, schur_rhs_momentum_projection_active, schur_rhs_momentum_projection_ok, &
         trim(schur_rhs_failure_mode), schur_momentum_consistent_rhs_active, schur_combined_merit_initial, &
         schur_combined_merit_final, schur_delta_x_norm, schur_delta_t_norm, schur_nonzero_update_accepted, &
         schur_pre_state_effective_ok, schur_rx_norm, schur_alpha_rx_norm, schur_w_proj_norm, schur_delta_ft_norm, &
         schur_alpha_delta_ft_norm, schur_zcorr_norm, schur_x_trial_update_norm, schur_internal_update_failure, &
         schur_w_projection_residual_norm, schur_w_projection_relative_residual, schur_w_projection_negligible, &
         schur_z_projection_residual_norm, schur_z_projection_relative_residual, schur_z_projection_negligible, &
         schur_zero_update_warning, trim(schur_projection_failure_mode), schur_tension_response_ok)

    deallocate(xnm1, xn, xnp1, xinit, fext, t_half)
  end subroutine run_flexible_structure_dynamics_test

  subroutine run_flexible_constraint_test()
    real(mytype), allocatable :: xnm1(:,:), xn(:,:), xnp1(:,:), xinit(:,:), fext(:,:), t_half(:)
    real(mytype) :: max_seg_err, max_inext_err, max_abs_tension
    real(mytype) :: max_seg_global, max_inext_global, max_t_global, tnow
    real(mytype) :: max_abs_drift_term_step, max_abs_vel_term_step, max_abs_force_term_step
    real(mytype) :: max_abs_drift_term_global, max_abs_vel_term_global, max_abs_force_term_global
    real(mytype) :: implicit_bending_update_norm_step, implicit_bending_update_norm_global
    real(mytype) :: post_bending_correction_norm_step, post_bending_correction_norm_global
    real(mytype) :: post_bending_max_inext_err_after_correction_step, post_bending_max_inext_err_after_correction_global
    real(mytype) :: final_coupled_residual_x_step, final_coupled_residual_g_step
    real(mytype) :: final_coupled_residual_x_rel_step, final_coupled_residual_x_rel_global
    real(mytype) :: max_abs_constraint_residual_during_outer_step, max_abs_momentum_residual_during_outer_step
    real(mytype) :: max_abs_constraint_residual_during_outer_global, max_abs_momentum_residual_during_outer_global
    real(mytype) :: max_abs_delta_x_during_outer_step, max_abs_delta_t_during_outer_step, accepted_line_search_lambda_step
    real(mytype) :: max_abs_delta_x_during_outer_global, max_abs_delta_t_during_outer_global
    real(mytype) :: accepted_line_search_lambda_global
    logical :: plateau_detected_step, plateau_detected_global
    real(mytype) :: final_coupled_residual_x_global, final_coupled_residual_g_global
    character(len=32) :: final_coupled_convergence_mode_step, final_coupled_convergence_mode_global
    integer :: max_outer_iterations_used_step, max_outer_iterations_used_global
    logical :: coupled_solver_converged_step, coupled_solver_converged_all
    logical :: coupled_solver_converged_strict_step, coupled_solver_converged_effective_step
    real(mytype) :: case_metric
    integer :: step, outint, i

    if (.not.fiber_active .or. .not.fiber_flexible_active .or. .not.fiber_flex_initialized .or. &
         .not.fiber_flex_constraint_test_active) then
      if (nrank == 0) write(*,*) 'Error: invalid state for flexible constraint test.'
      stop
    endif

    allocate(xnm1(3,fiber_nl), xn(3,fiber_nl), xnp1(3,fiber_nl), xinit(3,fiber_nl), fext(3,fiber_nl), t_half(fiber_nl-1))
    if (allocated(fiber_tension_half_prevstep)) fiber_tension_half_prevstep = 0._mytype
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
    max_abs_constraint_residual_during_outer_global = 0._mytype
    max_abs_momentum_residual_during_outer_global = 0._mytype
    max_abs_delta_x_during_outer_global = 0._mytype
    max_abs_delta_t_during_outer_global = 0._mytype
    accepted_line_search_lambda_global = 0._mytype
    final_coupled_residual_x_global = 0._mytype
    final_coupled_residual_g_global = 0._mytype
    final_coupled_residual_x_rel_global = 0._mytype
    plateau_detected_global = .false.
    final_coupled_convergence_mode_global = 'strict'
    max_outer_iterations_used_global = 0
    coupled_solver_converged_all = .true.
    outint = max(1, fiber_flex_constraint_output_interval)
    call write_fiber_flex_constraint_series(0, 0._mytype, 0._mytype, 0._mytype, 0._mytype, &
         0._mytype, 0._mytype, 0._mytype, 0._mytype, 0._mytype, 0._mytype, &
         0, 0._mytype, 0._mytype, 0._mytype, 0._mytype, 0._mytype, 0._mytype, .false., .true.)

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
           post_bending_correction_norm_step, post_bending_max_inext_err_after_correction_step, &
           max_outer_iterations_used_step, final_coupled_residual_x_step, final_coupled_residual_g_step, &
           coupled_solver_converged_step, coupled_solver_converged_strict_step, coupled_solver_converged_effective_step, &
           max_abs_constraint_residual_during_outer_step, max_abs_momentum_residual_during_outer_step, &
           max_abs_delta_x_during_outer_step, max_abs_delta_t_during_outer_step, accepted_line_search_lambda_step, &
           final_coupled_residual_x_rel_step, final_coupled_convergence_mode_step, plateau_detected_step)
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
      max_abs_constraint_residual_during_outer_global = max(max_abs_constraint_residual_during_outer_global, &
           max_abs_constraint_residual_during_outer_step)
      max_abs_momentum_residual_during_outer_global = max(max_abs_momentum_residual_during_outer_global, &
           max_abs_momentum_residual_during_outer_step)
      max_abs_delta_x_during_outer_global = max(max_abs_delta_x_during_outer_global, max_abs_delta_x_during_outer_step)
      max_abs_delta_t_during_outer_global = max(max_abs_delta_t_during_outer_global, max_abs_delta_t_during_outer_step)
      accepted_line_search_lambda_global = max(accepted_line_search_lambda_global, accepted_line_search_lambda_step)
      final_coupled_residual_x_global = max(final_coupled_residual_x_global, final_coupled_residual_x_step)
      final_coupled_residual_g_global = max(final_coupled_residual_g_global, final_coupled_residual_g_step)
      final_coupled_residual_x_rel_global = max(final_coupled_residual_x_rel_global, final_coupled_residual_x_rel_step)
      plateau_detected_global = plateau_detected_global .or. plateau_detected_step
      if (final_coupled_convergence_mode_step == 'failed') final_coupled_convergence_mode_global = 'failed'
      if (final_coupled_convergence_mode_step == 'plateau_accepted' .and. final_coupled_convergence_mode_global /= 'failed') then
        final_coupled_convergence_mode_global = 'plateau_accepted'
      endif
      max_outer_iterations_used_global = max(max_outer_iterations_used_global, max_outer_iterations_used_step)
      coupled_solver_converged_all = coupled_solver_converged_all .and. coupled_solver_converged_step

      tnow = real(step,mytype) * fiber_flex_constraint_dt
      if (mod(step,outint) == 0 .or. step == fiber_flex_constraint_nsteps) then
        call write_fiber_flex_constraint_series(step, tnow, max_seg_err, max_inext_err, max_abs_tension, &
             max_abs_drift_term_step, max_abs_vel_term_step, max_abs_force_term_step, &
             implicit_bending_update_norm_step, post_bending_correction_norm_step, &
             post_bending_max_inext_err_after_correction_step, max_outer_iterations_used_step, &
             final_coupled_residual_x_step, final_coupled_residual_g_step, max_abs_delta_x_during_outer_step, &
             max_abs_delta_t_during_outer_step, accepted_line_search_lambda_step, final_coupled_residual_x_rel_step, &
             plateau_detected_step, .false.)
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
         post_bending_max_inext_err_after_correction_global, 1.0_mytype, max_outer_iterations_used_global, &
         final_coupled_residual_x_global, final_coupled_residual_g_global, coupled_solver_converged_all, &
         max_abs_constraint_residual_during_outer_global, max_abs_momentum_residual_during_outer_global, &
         max_abs_delta_x_during_outer_global, max_abs_delta_t_during_outer_global, accepted_line_search_lambda_global, &
         final_coupled_residual_x_rel_global, final_coupled_convergence_mode_global, plateau_detected_global)

    deallocate(xnm1, xn, xnp1, xinit, fext, t_half)
  end subroutine run_flexible_constraint_test

end module fiber_flex_constraint
