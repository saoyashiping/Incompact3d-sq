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
       rigid_two_way_subiterations, rigid_two_way_subiter_slip_tol, rigid_two_way_subiter_verbose, &
       rigid_two_way_startup_fit, rigid_two_way_initialized, &
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

  subroutine startup_fit_two_way_state(uxe, uye, uze, itime)

    real(mytype), intent(in), dimension(:,:,:) :: uxe, uye, uze
    integer, intent(in) :: itime

    real(mytype) :: normal_mat(6,6), rhs(6), sol(6)
    real(mytype) :: row1(6), row2(6), row3(6), r(3), wq, diag_scale, wsum
    integer :: l, i, j
    logical :: solved

    call update_rigid_free_geometry()
    call update_rigid_free_xdot()
    call run_fiber_interp_solver_readonly(uxe, uye, uze, itime)

    normal_mat = 0._mytype
    rhs = 0._mytype
    wsum = 0._mytype

    do l = 1, fiber_nl
      r(1) = periodic_delta(fiber_x(1,l) - fiber_xc(1), xlx)
      r(2) = fiber_x(2,l) - fiber_xc(2)
      r(3) = periodic_delta(fiber_x(3,l) - fiber_xc(3), zlz)

      if (allocated(fiber_quad_w)) then
        wq = max(0._mytype, fiber_quad_w(l))
      else
        wq = 1._mytype
      endif
      wsum = wsum + wq

      row1 = (/1._mytype, 0._mytype, 0._mytype, 0._mytype, r(3), -r(2)/)
      row2 = (/0._mytype, 1._mytype, 0._mytype, -r(3), 0._mytype, r(1)/)
      row3 = (/0._mytype, 0._mytype, 1._mytype, r(2), -r(1), 0._mytype/)

      do i = 1, 6
        do j = 1, 6
          normal_mat(i,j) = normal_mat(i,j) + wq * row1(i) * row1(j)
          normal_mat(i,j) = normal_mat(i,j) + wq * row2(i) * row2(j)
          normal_mat(i,j) = normal_mat(i,j) + wq * row3(i) * row3(j)
        enddo
      enddo
      rhs = rhs + wq * row1 * fiber_uinterp(1,l) + wq * row2 * fiber_uinterp(2,l) + wq * row3 * fiber_uinterp(3,l)
    enddo

    if (wsum <= 0._mytype) wsum = real(fiber_nl, mytype)
    diag_scale = max(1._mytype, maxval(abs(normal_mat)))
    do i = 1, 6
      normal_mat(i,i) = normal_mat(i,i) + 1.0e-12_mytype * diag_scale
    enddo

    call solve_linear_6x6(normal_mat, rhs, sol, solved)

    if (solved) then
      fiber_uc = sol(1:3)
      fiber_omega = sol(4:6)
    else
      fiber_uc = sum(fiber_uinterp, dim=2) / real(fiber_nl, mytype)
      fiber_omega = 0._mytype
    endif

    call update_rigid_free_xdot()

  end subroutine startup_fit_two_way_state

  subroutine solve_linear_6x6(a_in, b_in, x_out, solved)

    real(mytype), intent(in) :: a_in(6,6), b_in(6)
    real(mytype), intent(out) :: x_out(6)
    logical, intent(out) :: solved

    real(mytype) :: a(6,6), b(6), pivot_row(6), pivot_val, factor, tmp
    integer :: i, j, k, p

    a = a_in
    b = b_in
    solved = .true.

    do k = 1, 6
      p = k
      pivot_val = abs(a(k,k))
      do i = k+1, 6
        if (abs(a(i,k)) > pivot_val) then
          pivot_val = abs(a(i,k))
          p = i
        endif
      enddo

      if (pivot_val <= 1.0e-14_mytype) then
        solved = .false.
        x_out = 0._mytype
        return
      endif

      if (p /= k) then
        pivot_row = a(k,:)
        a(k,:) = a(p,:)
        a(p,:) = pivot_row
        tmp = b(k)
        b(k) = b(p)
        b(p) = tmp
      endif

      do i = k+1, 6
        factor = a(i,k) / a(k,k)
        a(i,k:6) = a(i,k:6) - factor * a(k,k:6)
        b(i) = b(i) - factor * b(k)
      enddo
    enddo

    x_out = 0._mytype
    do i = 6, 1, -1
      x_out(i) = (b(i) - sum(a(i,i+1:6) * x_out(i+1:6))) / a(i,i)
    enddo

  end subroutine solve_linear_6x6

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
    real(mytype) :: relax_alpha, alpha_u, alpha_omega, slip_max_iter, slip_rms_iter
    real(mytype) :: uc_norm, omega_norm, xc_abs_max, lag_force_norm
    real(mytype) :: torque_used_perp(3), omega_dot_used(3)
    real(mytype) :: uc_new_raw(3), xc_new_raw(3), omega_new_raw(3), p_new_raw(3)
    real(mytype) :: force_seed_norm
    real(mytype), allocatable :: force_seed(:,:), uxe_eff(:,:,:), uye_eff(:,:,:), uze_eff(:,:,:)
    real(mytype), parameter :: fail_uc_limit = 1.0e3_mytype
    real(mytype), parameter :: fail_omega_limit = 1.0e5_mytype
    real(mytype), parameter :: fail_xc_limit_factor = 10._mytype
    real(mytype), parameter :: fail_lag_force_limit = 1.0e8_mytype
    real(mytype), parameter :: fail_slip_limit = 1.0e5_mytype
    real(mytype), parameter :: fail_spacing_limit = 0.25_mytype
    real(mytype), parameter :: fail_pnorm_error_limit = 0.25_mytype
    integer :: l, npts, ierr, isubiter, nsubiter
    integer :: failure_code
    logical :: output_now
    logical :: failed_flag, is_finite
    character(len=64) :: fail_stage, fail_quantity

    if (.not.fiber_active) return
    if (.not.rigid_two_way_test_active) return

    if (.not.rigid_two_way_initialized) then
      if (rigid_two_way_startup_fit) then
        call startup_fit_two_way_state(uxe, uye, uze, itime)
      endif
      rigid_two_way_initialized = .true.
    endif

    failed_flag = .false.
    failure_code = 0
    fail_stage = 'none'
    fail_quantity = 'none'

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

    nsubiter = max(1, rigid_two_way_subiterations)
    allocate(force_seed(3, fiber_nl), uxe_eff(size(uxe,1), size(uxe,2), size(uxe,3)), &
         uye_eff(size(uye,1), size(uye,2), size(uye,3)), uze_eff(size(uze,1), size(uze,2), size(uze,3)))
    force_seed = rigid_two_way_force_prev
    uxe_eff = uxe
    uye_eff = uye
    uze_eff = uze

    do isubiter = 1, nsubiter

      ! Stage A (prepare forcing): rigid geometry/interpolation/slip/force/spread
      call update_rigid_free_geometry()
      call update_rigid_free_xdot()
      call run_fiber_interp_solver_readonly(uxe_eff, uye_eff, uze_eff, itime)

      fiber_slip = fiber_uinterp - fiber_xdot
      fiber_coupling_force = beta_eff * fiber_slip
      fiber_coupling_force = relax_alpha * fiber_coupling_force + (1._mytype - relax_alpha) * force_seed

      if (.not.allocated(fiber_euler_force_x)) allocate(fiber_euler_force_x(size(uxe,1), size(uxe,2), size(uxe,3)))
      if (.not.allocated(fiber_euler_force_y)) allocate(fiber_euler_force_y(size(uxe,1), size(uxe,2), size(uxe,3)))
      if (.not.allocated(fiber_euler_force_z)) allocate(fiber_euler_force_z(size(uxe,1), size(uxe,2), size(uxe,3)))
      fiber_euler_force_x = 0._mytype
      fiber_euler_force_y = 0._mytype
      fiber_euler_force_z = 0._mytype

      call spread_lagrangian_force_to_euler(fiber_coupling_force, fiber_euler_force_x, fiber_euler_force_y, fiber_euler_force_z, &
           eul_total_local, sumw_min, sumw_max, xp, yp, zp, spread_hy_loc_min, spread_hy_loc_max)

      lag_total_local = 0._mytype
      torque_local = 0._mytype
      do l = 1, fiber_nl
        lag_total_local = lag_total_local + fiber_coupling_force(:,l) * fiber_quad_w(l)
        r(1) = periodic_delta(fiber_x(1,l) - fiber_xc(1), xlx)
        r(2) = fiber_x(2,l) - fiber_xc(2)
        r(3) = periodic_delta(fiber_x(3,l) - fiber_xc(3), zlz)
        torque_local = torque_local + cross_product(r, fiber_coupling_force(:,l)) * fiber_quad_w(l)
      enddo
      call MPI_ALLREDUCE(lag_total_local, lag_total, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(torque_local, fiber_torque_total, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(eul_total_local, eul_total, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      fiber_force_total = lag_total

      is_finite = all(ieee_is_finite(fiber_uinterp)) .and. all(ieee_is_finite(fiber_slip)) .and. all(ieee_is_finite(fiber_coupling_force))
      if (.not.is_finite .and. .not.failed_flag) then
        failed_flag = .true.
        failure_code = 201
        fail_stage = 'prepare_forcing'
        fail_quantity = 'interp_slip_force'
      endif

      ! Stage B (finalize update): rigid-body update from updated force/torque
      torque_used_perp = fiber_torque_total - dot_product(fiber_torque_total, fiber_p) * fiber_p
      omega_dot_used = torque_used_perp / fiber_inertia_perp
      uc_new_raw = fiber_uc + dt_stage * (fiber_force_total / fiber_mass)
      fiber_uc = alpha_u * uc_new_raw + (1._mytype - alpha_u) * fiber_uc
      xc_new_raw = fiber_xc + dt_stage * fiber_uc
      fiber_xc = xc_new_raw
      fiber_xc(1) = wrap_periodic(fiber_xc(1), xlx)
      fiber_xc(3) = wrap_periodic(fiber_xc(3), zlz)
      omega_new_raw = fiber_omega + dt_stage * omega_dot_used
      fiber_omega = alpha_omega * omega_new_raw + (1._mytype - alpha_omega) * fiber_omega
      p_new_raw = fiber_p + dt_stage * cross_product(fiber_omega, fiber_p)
      pnorm = sqrt(sum(p_new_raw**2))
      if (pnorm > 0._mytype) then
        fiber_p = p_new_raw / pnorm
      else
        fiber_p = p_new_raw
      endif

      call update_rigid_free_geometry()
      call update_rigid_free_xdot()
      call run_fiber_interp_solver_readonly(uxe_eff, uye_eff, uze_eff, itime)
      fiber_slip = fiber_uinterp - fiber_xdot
      slip_max_iter = maxval(abs(fiber_slip))
      npts = 3 * fiber_nl
      slip_rms_iter = sqrt(sum(fiber_slip**2) / real(npts, mytype))

      force_seed = fiber_coupling_force
      force_seed_norm = sqrt(sum(force_seed**2))
      if (rigid_two_way_subiter_verbose .and. nrank == 0) then
        write(*,'(A,I3,A,ES12.4,A,ES12.4,A,ES12.4)') 'two_way_subiter=', isubiter, ' slip_max=', slip_max_iter, &
             ' slip_rms=', slip_rms_iter, ' force_seed_norm=', force_seed_norm
      endif

      if (slip_max_iter <= max(0._mytype, rigid_two_way_subiter_slip_tol)) exit

      if (isubiter < nsubiter) then
        uxe_eff = uxe + dt_stage * fiber_euler_force_x
        uye_eff = uye + dt_stage * fiber_euler_force_y
        uze_eff = uze + dt_stage * fiber_euler_force_z
      endif

    enddo

    rigid_two_way_force_prev = fiber_coupling_force

    call update_rigid_free_geometry()
    call update_rigid_free_xdot()
    call run_fiber_interp_solver_readonly(uxe, uye, uze, itime)
    fiber_slip = fiber_uinterp - fiber_xdot
    call compute_rigid_free_spacing_error(spacing_error_max)

    slip_max = maxval(abs(fiber_slip))
    npts = 3 * fiber_nl
    slip_rms = sqrt(sum(fiber_slip**2) / real(npts, mytype))
    p_norm_error = abs(sqrt(sum(fiber_p**2)) - 1._mytype)
    abs_force_balance = abs(eul_total - fiber_force_total)
    lag_impulse = fiber_force_total * dt_stage
    euler_impulse = eul_total * dt_stage

    uc_norm = sqrt(fiber_uc(1)**2 + fiber_uc(2)**2 + fiber_uc(3)**2)
    omega_norm = sqrt(fiber_omega(1)**2 + fiber_omega(2)**2 + fiber_omega(3)**2)
    xc_abs_max = maxval(abs(fiber_xc))
    lag_force_norm = sqrt(fiber_force_total(1)**2 + fiber_force_total(2)**2 + fiber_force_total(3)**2)

    is_finite = all(ieee_is_finite(fiber_uc)) .and. all(ieee_is_finite(fiber_xc)) .and. all(ieee_is_finite(fiber_omega)) .and. &
         all(ieee_is_finite(fiber_p)) .and. all(ieee_is_finite(fiber_x))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 202
      fail_stage = 'post_rigid_update'
      fail_quantity = 'non_finite_state'
    endif

    if ((.not.ieee_is_finite(spacing_error_max) .or. .not.ieee_is_finite(p_norm_error)) .and. .not.failed_flag) then
      failed_flag = .true.
      failure_code = 203
      fail_stage = 'post_geometry'
      fail_quantity = 'spacing_or_pnorm_error'
    endif

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
           fiber_force_total, eul_total)
      call write_fiber_rigid_two_way_momentum(itime, time, fiber_force_total, eul_total, abs_force_balance, dt_stage, &
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

    deallocate(force_seed, uxe_eff, uye_eff, uze_eff)

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
