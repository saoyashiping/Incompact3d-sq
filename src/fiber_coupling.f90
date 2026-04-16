!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_coupling

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank, nproc
  use mpi, only : MPI_COMM_WORLD, MPI_ABORT, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_SUM
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use param, only : ifirst, xlx, yly, zlz, gdt
  use variables, only : xp, yp, zp
  use fiber_types, only : fiber_active, fiber_nl, rigid_coupling_test_active, rigid_free_test_active, &
       rigid_two_way_test_active, ibm_beta, coupling_ramp_steps, rigid_output_interval, rigid_two_way_output_interval, &
       rigid_two_way_force_relaxation, rigid_two_way_velocity_relaxation, rigid_two_way_omega_relaxation, &
       fiber_x, fiber_xc, fiber_uc, fiber_p, fiber_omega, fiber_force_total, fiber_torque_total, fiber_xdot, &
       fiber_uinterp, fiber_slip, fiber_coupling_force, fiber_quad_w, fiber_euler_force_x, fiber_euler_force_y, &
       fiber_euler_force_z, rigid_motion_case, fiber_mass, fiber_inertia_perp
  use fiber_rigid_motion, only : update_rigid_motion, compute_rigid_spacing_error
  use fiber_rigid_free, only : wrap_periodic, periodic_delta, cross_product, compute_rigid_free_spacing_error, &
       update_rigid_free_xdot, update_rigid_free_geometry
  use fiber_interp, only : run_fiber_interp_solver_readonly
  use fiber_spread, only : spread_lagrangian_force_to_euler
  use fiber_io, only : write_fiber_rigid_coupling_points, write_fiber_rigid_coupling_summary, &
       write_fiber_rigid_two_way_points, write_fiber_rigid_two_way_summary, write_fiber_rigid_two_way_momentum

  implicit none

  real(mytype), allocatable, save :: fiber_coupling_force_prev(:,:)
  real(mytype), allocatable, save :: rigid_two_way_force_prev(:,:)
  real(mytype), parameter :: coupling_force_relaxation = 0.25_mytype

contains

  subroutine run_rigid_coupling_step(uxe, uye, uze, time, itime, isubstep, nsubsteps)

    real(mytype), intent(in), dimension(:,:,:) :: uxe, uye, uze
    real(mytype), intent(in) :: time
    integer, intent(in) :: itime
    integer, intent(in) :: isubstep, nsubsteps

    real(mytype) :: slip_max, slip_rms, spacing_error_max
    real(mytype) :: lag_total(3), eul_total(3), abs_force_balance(3)
    real(mytype) :: sumw_min, sumw_max, spread_hy_loc_min, spread_hy_loc_max
    real(mytype) :: u_interp_max_norm, xdot_max_norm, coupling_force_max_norm, euler_force_max_norm
    real(mytype) :: beta_eff, ramp_factor, ramp_linear, relax_alpha
    logical :: failed_flag
    logical :: is_finite
    character(len=64) :: fail_quantity
    integer :: failure_code
    integer :: npts, coupling_step
    integer :: ierr
    logical :: output_now

    if (.not.fiber_active) return
    if (.not.rigid_coupling_test_active .and. .not.rigid_free_test_active) return

    if (nproc /= 1) then
      if (nrank == 0) write(*,*) 'Error: rigid_coupling_test_active currently supports only nproc = 1.'
      stop
    endif

    call update_rigid_motion(time)
    call run_fiber_interp_solver_readonly(uxe, uye, uze, itime)

    if (.not.allocated(fiber_slip)) allocate(fiber_slip(3, fiber_nl))
    if (.not.allocated(fiber_coupling_force)) allocate(fiber_coupling_force(3, fiber_nl))
    fiber_slip = 0._mytype
    fiber_coupling_force = 0._mytype
    relax_alpha = max(0._mytype, min(1._mytype, coupling_force_relaxation))

    if (.not.allocated(fiber_coupling_force_prev)) then
      allocate(fiber_coupling_force_prev(3, fiber_nl))
      fiber_coupling_force_prev = 0._mytype
    else if (size(fiber_coupling_force_prev,2) /= fiber_nl) then
      deallocate(fiber_coupling_force_prev)
      allocate(fiber_coupling_force_prev(3, fiber_nl))
      fiber_coupling_force_prev = 0._mytype
    endif

    coupling_step = max(1, itime - ifirst + 1)
    if (coupling_ramp_steps > 0) then
      ramp_linear = real(coupling_step, mytype) / real(coupling_ramp_steps, mytype)
      ramp_linear = min(1._mytype, ramp_linear)
      ramp_factor = ramp_linear * ramp_linear * (3._mytype - 2._mytype * ramp_linear)
    else
      ramp_factor = 1._mytype
    endif
    beta_eff = ibm_beta * ramp_factor
    fiber_slip = fiber_uinterp - fiber_xdot
    fiber_coupling_force = beta_eff * fiber_slip
    fiber_coupling_force = relax_alpha * fiber_coupling_force + (1._mytype - relax_alpha) * fiber_coupling_force_prev
    fiber_coupling_force_prev = fiber_coupling_force

    if (.not.allocated(fiber_euler_force_x)) allocate(fiber_euler_force_x(size(uxe,1), size(uxe,2), size(uxe,3)))
    if (.not.allocated(fiber_euler_force_y)) allocate(fiber_euler_force_y(size(uxe,1), size(uxe,2), size(uxe,3)))
    if (.not.allocated(fiber_euler_force_z)) allocate(fiber_euler_force_z(size(uxe,1), size(uxe,2), size(uxe,3)))
    fiber_euler_force_x = 0._mytype
    fiber_euler_force_y = 0._mytype
    fiber_euler_force_z = 0._mytype

    call spread_lagrangian_force_to_euler(fiber_coupling_force, fiber_euler_force_x, fiber_euler_force_y, fiber_euler_force_z, &
         eul_total, sumw_min, sumw_max, xp, yp, zp, spread_hy_loc_min, spread_hy_loc_max)

    lag_total(1) = sum(fiber_coupling_force(1,:) * fiber_quad_w(:))
    lag_total(2) = sum(fiber_coupling_force(2,:) * fiber_quad_w(:))
    lag_total(3) = sum(fiber_coupling_force(3,:) * fiber_quad_w(:))

    abs_force_balance = abs(eul_total - lag_total)

    slip_max = maxval(abs(fiber_slip))
    npts = 3 * fiber_nl
    slip_rms = sqrt(sum(fiber_slip**2) / real(npts, mytype))
    u_interp_max_norm = maxval(sqrt(fiber_uinterp(1,:)**2 + fiber_uinterp(2,:)**2 + fiber_uinterp(3,:)**2))
    xdot_max_norm = maxval(sqrt(fiber_xdot(1,:)**2 + fiber_xdot(2,:)**2 + fiber_xdot(3,:)**2))
    coupling_force_max_norm = maxval(sqrt(fiber_coupling_force(1,:)**2 + fiber_coupling_force(2,:)**2 + &
         fiber_coupling_force(3,:)**2))
    euler_force_max_norm = maxval(sqrt(fiber_euler_force_x**2 + fiber_euler_force_y**2 + fiber_euler_force_z**2))

    call compute_rigid_spacing_error(spacing_error_max)

    failed_flag = .false.
    fail_quantity = 'none'
    failure_code = 0

    is_finite = all(ieee_is_finite(fiber_slip))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      fail_quantity = 'fiber_slip'
      failure_code = 1
    endif
    is_finite = all(ieee_is_finite(fiber_coupling_force))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      fail_quantity = 'fiber_coupling_force'
      failure_code = 2
    endif
    is_finite = all(ieee_is_finite(fiber_euler_force_x)) .and. all(ieee_is_finite(fiber_euler_force_y)) .and. &
         all(ieee_is_finite(fiber_euler_force_z))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      fail_quantity = 'fiber_euler_force_xyz'
      failure_code = 3
    endif
    is_finite = all(ieee_is_finite((/slip_max, slip_rms, spacing_error_max, u_interp_max_norm, xdot_max_norm, &
         coupling_force_max_norm, euler_force_max_norm/))) .and. &
         all(ieee_is_finite(lag_total)) .and. all(ieee_is_finite(eul_total)) .and. all(ieee_is_finite(abs_force_balance))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      fail_quantity = 'rigid_summary_scalars'
      failure_code = 4
    endif

    output_now = (isubstep == nsubsteps) .and. ((itime == ifirst) .or. (mod(itime, max(1, rigid_output_interval)) == 0))
    if (failed_flag) output_now = .true.

    if (nrank == 0 .and. itime == ifirst) then
      write(*,'(A,I10,A,I10)') 'DEBUG_RAMP ifirst=', ifirst, ' itime=', itime
      write(*,'(A,I10)') 'DEBUG_RAMP coupling_ramp_steps=', coupling_ramp_steps
      write(*,'(A,I10)') 'DEBUG_RAMP coupling_step=', coupling_step
      write(*,'(A,ES12.4)') 'DEBUG_RAMP ramp_factor=', ramp_factor
      write(*,'(A,ES12.4)') 'DEBUG_RAMP ibm_beta=', ibm_beta
      write(*,'(A,ES12.4)') 'DEBUG_RAMP beta_eff=', beta_eff
    endif

    if (output_now) then
      if (nrank == 0) then
        write(*,'(A,3ES12.4)') 'Rigid coupling beta_in/beta_eff/ramp_factor  : ', ibm_beta, beta_eff, ramp_factor
        write(*,'(A,2ES12.4)') 'Rigid coupling spread sumw min/max            : ', sumw_min, sumw_max
        write(*,'(A,2ES12.4)') 'Rigid coupling spread hy_loc min/max          : ', spread_hy_loc_min, spread_hy_loc_max
        write(*,'(A,4ES12.4)') 'Rigid coupling max norms (U, Xdot, Ffs, Feul): ', &
             u_interp_max_norm, xdot_max_norm, coupling_force_max_norm, euler_force_max_norm
      endif
      call write_fiber_rigid_coupling_points(itime)
      call write_fiber_rigid_coupling_summary(itime, time, slip_max, slip_rms, lag_total, eul_total, &
           abs_force_balance, spacing_error_max, u_interp_max_norm, xdot_max_norm, coupling_force_max_norm, &
           euler_force_max_norm, ibm_beta, beta_eff, ramp_factor, failed_flag, failure_code)
    endif

    if (failed_flag) then
      if (nrank == 0) then
        write(*,'(A)') 'RIGID COUPLING TEST FAILED'
        write(*,'(A,I6)') 'failure_code                                   : ', failure_code
        write(*,'(A,I10)') 'Failure at itime                               : ', itime
        write(*,'(A,ES12.4)') 'Failure time                                   : ', time
        write(*,'(A,I6)') 'Failure motion case                            : ', rigid_motion_case
        write(*,'(A,A)') 'first_nonfinite_quantity                       : ', trim(fail_quantity)
        write(*,'(A,ES12.4)') 'beta_input                                     : ', ibm_beta
        write(*,'(A,ES12.4)') 'beta_eff                                       : ', beta_eff
        write(*,'(A,ES12.4)') 'ramp_factor                                    : ', ramp_factor
      endif
      call MPI_ABORT(MPI_COMM_WORLD, 911, ierr)
    endif

  end subroutine run_rigid_coupling_step

  subroutine run_rigid_two_way_step(uxe, uye, uze, time, itime, isubstep, nsubsteps)

    real(mytype), intent(in), dimension(:,:,:) :: uxe, uye, uze
    real(mytype), intent(in) :: time
    integer, intent(in) :: itime
    integer, intent(in) :: isubstep, nsubsteps

    real(mytype) :: slip_max, slip_rms, spacing_error_max, p_norm_error
    real(mytype) :: lag_total_local(3), lag_total(3), eul_total_local(3), eul_total(3), abs_force_balance(3)
    real(mytype) :: sumw_min, sumw_max, spread_hy_loc_min, spread_hy_loc_max
    real(mytype) :: beta_eff, ramp_factor, ramp_linear, coupling_step, dt_stage
    real(mytype) :: r(3), torque_local(3), pnorm
    real(mytype) :: lag_impulse(3), euler_impulse(3)
    real(mytype) :: relax_alpha, alpha_u, alpha_omega
    real(mytype) :: uc_norm, omega_norm, xc_abs_max, lag_force_norm
    real(mytype) :: xc_old(3), uc_old(3), omega_old(3), p_old(3)
    real(mytype) :: xc_pred(3), uc_pred_raw(3), uc_pred_used(3)
    real(mytype) :: omega_pred_raw(3), omega_pred_used(3), p_pred(3)
    real(mytype) :: xc_new_raw(3), xc_used(3), uc_new_raw(3), uc_used(3)
    real(mytype) :: omega_new_raw(3), omega_used(3), p_new_raw(3), p_used(3)
    real(mytype) :: lag_total_pred_local(3), lag_total_pred(3), lag_total_corr_local(3), lag_total_corr(3)
    real(mytype) :: torque_pred_local(3), torque_pred(3), torque_corr_local(3), torque_corr(3)
    real(mytype) :: torque_pred_perp(3), omega_dot_pred(3)
    real(mytype) :: torque_used(3), torque_used_perp(3), omega_dot_used(3)
    real(mytype), allocatable :: x_backup(:,:), xdot_backup(:,:), uinterp_backup(:,:), slip_backup(:,:)
    real(mytype), allocatable :: x_pred(:,:), xdot_pred(:,:), uinterp_pred(:,:), slip_pred(:,:), force_pred(:,:)
    real(mytype), allocatable :: force_corr(:,:)
    real(mytype), parameter :: fail_uc_limit = 1.0e3_mytype
    real(mytype), parameter :: fail_omega_limit = 1.0e5_mytype
    real(mytype), parameter :: fail_xc_limit_factor = 10._mytype
    real(mytype), parameter :: fail_lag_force_limit = 1.0e8_mytype
    real(mytype), parameter :: fail_slip_limit = 1.0e5_mytype
    real(mytype), parameter :: fail_spacing_limit = 0.25_mytype
    real(mytype), parameter :: fail_pnorm_error_limit = 0.25_mytype
    integer :: l, npts, ierr
    integer :: failure_code
    logical :: output_now
    logical :: failed_flag, is_finite
    character(len=64) :: fail_stage, fail_quantity

    if (.not.fiber_active) return
    if (.not.rigid_two_way_test_active) return

    failed_flag = .false.
    failure_code = 0
    fail_stage = 'none'
    fail_quantity = 'none'

    call update_rigid_free_geometry()
    call update_rigid_free_xdot()
    call run_fiber_interp_solver_readonly(uxe, uye, uze, itime)

    if (.not.allocated(fiber_slip)) allocate(fiber_slip(3, fiber_nl))
    if (.not.allocated(fiber_coupling_force)) allocate(fiber_coupling_force(3, fiber_nl))
    fiber_slip = 0._mytype
    fiber_coupling_force = 0._mytype
    alpha_u = max(0._mytype, min(1._mytype, rigid_two_way_velocity_relaxation))
    alpha_omega = max(0._mytype, min(1._mytype, rigid_two_way_omega_relaxation))
    relax_alpha = max(0._mytype, min(1._mytype, rigid_two_way_force_relaxation))

    if (.not.allocated(rigid_two_way_force_prev)) then
      allocate(rigid_two_way_force_prev(3, fiber_nl))
      rigid_two_way_force_prev = 0._mytype
    else if (size(rigid_two_way_force_prev,2) /= fiber_nl) then
      deallocate(rigid_two_way_force_prev)
      allocate(rigid_two_way_force_prev(3, fiber_nl))
      rigid_two_way_force_prev = 0._mytype
    endif

    coupling_step = real(max(1, itime - ifirst + 1), mytype)
    dt_stage = gdt(isubstep)
    if (coupling_ramp_steps > 0) then
      ramp_linear = coupling_step / real(coupling_ramp_steps, mytype)
      ramp_linear = min(1._mytype, ramp_linear)
      ramp_factor = ramp_linear * ramp_linear * (3._mytype - 2._mytype * ramp_linear)
    else
      ramp_factor = 1._mytype
    endif
    beta_eff = ibm_beta * ramp_factor

    allocate(x_backup(3, fiber_nl), xdot_backup(3, fiber_nl), uinterp_backup(3, fiber_nl), slip_backup(3, fiber_nl), &
         x_pred(3, fiber_nl), xdot_pred(3, fiber_nl), uinterp_pred(3, fiber_nl), &
         slip_pred(3, fiber_nl), force_pred(3, fiber_nl), force_corr(3, fiber_nl))

    x_backup = fiber_x
    xdot_backup = fiber_xdot
    uinterp_backup = fiber_uinterp

    xc_old = fiber_xc
    uc_old = fiber_uc
    omega_old = fiber_omega
    p_old = fiber_p

    ! Predictor pass: use current state
    slip_pred = uinterp_backup - xdot_backup
    force_pred = beta_eff * slip_pred
    force_pred = relax_alpha * force_pred + (1._mytype - relax_alpha) * rigid_two_way_force_prev

    is_finite = all(ieee_is_finite(uinterp_backup))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 101
      fail_stage = 'post_interp'
      fail_quantity = 'fiber_uinterp'
    endif
    is_finite = all(ieee_is_finite(slip_pred))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 102
      fail_stage = 'post_interp'
      fail_quantity = 'slip_pred'
    endif

    lag_total_pred_local = 0._mytype
    torque_pred_local = 0._mytype
    do l = 1, fiber_nl
      lag_total_pred_local = lag_total_pred_local + force_pred(:,l) * fiber_quad_w(l)
      r(1) = periodic_delta(x_backup(1,l) - xc_old(1), xlx)
      r(2) = x_backup(2,l) - xc_old(2)
      r(3) = periodic_delta(x_backup(3,l) - xc_old(3), zlz)
      torque_pred_local = torque_pred_local + cross_product(r, force_pred(:,l)) * fiber_quad_w(l)
    enddo
    call MPI_ALLREDUCE(lag_total_pred_local, lag_total_pred, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(torque_pred_local, torque_pred, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    torque_pred_perp = torque_pred - dot_product(torque_pred, p_old) * p_old
    omega_dot_pred = torque_pred_perp / fiber_inertia_perp
    uc_pred_raw = uc_old + dt_stage * (lag_total_pred / fiber_mass)
    uc_pred_used = alpha_u * uc_pred_raw + (1._mytype - alpha_u) * uc_old
    xc_pred = xc_old + dt_stage * uc_pred_used
    xc_pred(1) = wrap_periodic(xc_pred(1), xlx)
    xc_pred(3) = wrap_periodic(xc_pred(3), zlz)
    omega_pred_raw = omega_old + dt_stage * omega_dot_pred
    omega_pred_used = alpha_omega * omega_pred_raw + (1._mytype - alpha_omega) * omega_old
    p_pred = p_old + dt_stage * cross_product(omega_pred_used, p_old)
    pnorm = sqrt(sum(p_pred**2))
    if (pnorm > 0._mytype) p_pred = p_pred / pnorm

    ! Corrector pass: use predicted state
    fiber_xc = xc_pred
    fiber_uc = uc_pred_used
    fiber_omega = omega_pred_used
    fiber_p = p_pred
    call update_rigid_free_geometry()
    call update_rigid_free_xdot()
    call run_fiber_interp_solver_readonly(uxe, uye, uze, itime)

    x_pred = fiber_x
    xdot_pred = fiber_xdot
    uinterp_pred = fiber_uinterp
    slip_backup = uinterp_pred - xdot_pred
    force_corr = beta_eff * slip_backup
    force_corr = relax_alpha * force_corr + (1._mytype - relax_alpha) * force_pred

    if (.not.allocated(fiber_euler_force_x)) allocate(fiber_euler_force_x(size(uxe,1), size(uxe,2), size(uxe,3)))
    if (.not.allocated(fiber_euler_force_y)) allocate(fiber_euler_force_y(size(uxe,1), size(uxe,2), size(uxe,3)))
    if (.not.allocated(fiber_euler_force_z)) allocate(fiber_euler_force_z(size(uxe,1), size(uxe,2), size(uxe,3)))
    fiber_euler_force_x = 0._mytype
    fiber_euler_force_y = 0._mytype
    fiber_euler_force_z = 0._mytype

    fiber_coupling_force = 0.5_mytype * (force_pred + force_corr)
    rigid_two_way_force_prev = fiber_coupling_force

    call spread_lagrangian_force_to_euler(fiber_coupling_force, fiber_euler_force_x, fiber_euler_force_y, fiber_euler_force_z, &
         eul_total_local, sumw_min, sumw_max, xp, yp, zp, spread_hy_loc_min, spread_hy_loc_max)

    lag_total_corr_local = 0._mytype
    torque_corr_local = 0._mytype
    do l = 1, fiber_nl
      lag_total_corr_local = lag_total_corr_local + force_corr(:,l) * fiber_quad_w(l)
      r(1) = periodic_delta(x_pred(1,l) - xc_pred(1), xlx)
      r(2) = x_pred(2,l) - xc_pred(2)
      r(3) = periodic_delta(x_pred(3,l) - xc_pred(3), zlz)
      torque_corr_local = torque_corr_local + cross_product(r, force_corr(:,l)) * fiber_quad_w(l)
    enddo
    call MPI_ALLREDUCE(lag_total_corr_local, lag_total_corr, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(torque_corr_local, torque_corr, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    lag_total = 0.5_mytype * (lag_total_pred + lag_total_corr)
    torque_used = 0.5_mytype * (torque_pred + torque_corr)
    fiber_force_total = lag_total
    fiber_torque_total = torque_used
    call MPI_ALLREDUCE(eul_total_local, eul_total, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    abs_force_balance = abs(eul_total - fiber_force_total)

    is_finite = all(ieee_is_finite(fiber_coupling_force))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 103
      fail_stage = 'post_force_build'
      fail_quantity = 'force_used'
    endif
    is_finite = all(ieee_is_finite(lag_total))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 104
      fail_stage = 'post_force_build'
      fail_quantity = 'lag_total_force_used'
    endif
    is_finite = all(ieee_is_finite(fiber_torque_total))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 105
      fail_stage = 'post_force_build'
      fail_quantity = 'lag_total_torque'
    endif

    torque_used_perp = fiber_torque_total - dot_product(fiber_torque_total, p_old) * p_old
    omega_dot_used = torque_used_perp / fiber_inertia_perp
    uc_new_raw = uc_old + dt_stage * (fiber_force_total / fiber_mass)
    uc_used = alpha_u * uc_new_raw + (1._mytype - alpha_u) * uc_old
    xc_new_raw = xc_old + dt_stage * uc_used
    xc_used = xc_new_raw
    xc_used(1) = wrap_periodic(xc_used(1), xlx)
    xc_used(3) = wrap_periodic(xc_used(3), zlz)
    omega_new_raw = omega_old + dt_stage * omega_dot_used
    omega_used = alpha_omega * omega_new_raw + (1._mytype - alpha_omega) * omega_old
    p_new_raw = p_old + dt_stage * cross_product(omega_used, p_old)
    pnorm = sqrt(sum(p_new_raw**2))
    p_used = p_new_raw
    if (pnorm > 0._mytype) p_used = p_new_raw / pnorm

    fiber_uc = uc_used
    fiber_xc = xc_used
    fiber_omega = omega_used
    fiber_p = p_used
    p_norm_error = abs(sqrt(sum(fiber_p**2)) - 1._mytype)

    is_finite = all(ieee_is_finite(fiber_uc))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 106
      fail_stage = 'post_rigid_update'
      fail_quantity = 'fiber_uc'
    endif
    is_finite = all(ieee_is_finite(fiber_xc))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 107
      fail_stage = 'post_rigid_update'
      fail_quantity = 'fiber_xc'
    endif
    is_finite = all(ieee_is_finite(fiber_omega))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 108
      fail_stage = 'post_rigid_update'
      fail_quantity = 'fiber_omega'
    endif
    is_finite = all(ieee_is_finite(fiber_p))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 109
      fail_stage = 'post_rigid_update'
      fail_quantity = 'fiber_p'
    endif

    call update_rigid_free_geometry()
    call update_rigid_free_xdot()
    call run_fiber_interp_solver_readonly(uxe, uye, uze, itime)
    fiber_slip = fiber_uinterp - fiber_xdot
    call compute_rigid_free_spacing_error(spacing_error_max)

    slip_max = maxval(abs(fiber_slip))
    npts = 3 * fiber_nl
    slip_rms = sqrt(sum(fiber_slip**2) / real(npts, mytype))
    lag_impulse = lag_total * dt_stage
    euler_impulse = eul_total * dt_stage

    is_finite = all(ieee_is_finite(fiber_x))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 110
      fail_stage = 'post_geometry'
      fail_quantity = 'fiber_x'
    endif
    is_finite = ieee_is_finite(spacing_error_max) .and. ieee_is_finite(p_norm_error)
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 111
      fail_stage = 'post_geometry'
      fail_quantity = 'spacing_or_pnorm_error'
    endif

    uc_norm = sqrt(fiber_uc(1)**2 + fiber_uc(2)**2 + fiber_uc(3)**2)
    omega_norm = sqrt(fiber_omega(1)**2 + fiber_omega(2)**2 + fiber_omega(3)**2)
    xc_abs_max = maxval(abs(fiber_xc))
    lag_force_norm = sqrt(lag_total(1)**2 + lag_total(2)**2 + lag_total(3)**2)
    if (uc_norm > fail_uc_limit .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 112
      fail_stage = 'threshold'
      fail_quantity = '|fiber_uc|'
    endif
    if (omega_norm > fail_omega_limit .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 113
      fail_stage = 'threshold'
      fail_quantity = '|fiber_omega|'
    endif
    if (xc_abs_max > fail_xc_limit_factor * max(xlx, max(yly, zlz)) .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 114
      fail_stage = 'threshold'
      fail_quantity = '|fiber_xc|'
    endif
    if (lag_force_norm > fail_lag_force_limit .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 115
      fail_stage = 'threshold'
      fail_quantity = '|lag_total_force|'
    endif
    if (slip_max > fail_slip_limit .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 116
      fail_stage = 'threshold'
      fail_quantity = 'slip_max'
    endif
    if (spacing_error_max > fail_spacing_limit .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 117
      fail_stage = 'threshold'
      fail_quantity = 'spacing_error_max'
    endif
    if (p_norm_error > fail_pnorm_error_limit .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 118
      fail_stage = 'threshold'
      fail_quantity = 'p_norm_error'
    endif

    output_now = (isubstep == nsubsteps) .and. &
         ((itime == ifirst) .or. (mod(itime, max(1, rigid_two_way_output_interval)) == 0))
    if (failed_flag) output_now = .true.
    if (output_now) then
      call write_fiber_rigid_two_way_points(itime)
      call write_fiber_rigid_two_way_summary(itime, time, slip_max, slip_rms, spacing_error_max, p_norm_error, &
           lag_total, eul_total)
      call write_fiber_rigid_two_way_momentum(itime, time, lag_total, eul_total, abs_force_balance, dt_stage, &
           lag_impulse, euler_impulse)
    endif

    if (failed_flag) then
      if (nrank == 0) then
        write(*,'(A)') 'RIGID TWO-WAY TEST FAILED (FAIL-FAST)'
        write(*,'(A,I6)') 'failure_code                                   : ', failure_code
        write(*,'(A,I10)') 'Failure at itime                               : ', itime
        write(*,'(A,ES12.4)') 'Failure time                                   : ', time
        write(*,'(A,A)') 'failure_stage                                  : ', trim(fail_stage)
        write(*,'(A,A)') 'failure_quantity                               : ', trim(fail_quantity)
      endif
      call MPI_ABORT(MPI_COMM_WORLD, 913, ierr)
    endif

    deallocate(x_backup, xdot_backup, uinterp_backup, slip_backup, x_pred, xdot_pred, uinterp_pred, &
         slip_pred, force_pred, force_corr)

  end subroutine run_rigid_two_way_step

  subroutine add_rigid_coupling_rhs(dux, duy, duz)

    real(mytype), intent(inout), dimension(:,:,:) :: dux, duy, duz

    if (.not.rigid_coupling_test_active .and. .not.rigid_free_test_active .and. .not.rigid_two_way_test_active) return
    if (.not.allocated(fiber_euler_force_x)) return

    dux = dux + fiber_euler_force_x
    duy = duy + fiber_euler_force_y
    duz = duz + fiber_euler_force_z

  end subroutine add_rigid_coupling_rhs

end module fiber_coupling
