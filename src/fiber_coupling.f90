!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_coupling

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank, nproc
  use mpi, only : MPI_COMM_WORLD, MPI_ABORT, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_SUM
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use param, only : ifirst, xlx, zlz, gdt
  use variables, only : xp, yp, zp
  use fiber_types, only : fiber_active, fiber_nl, rigid_coupling_test_active, rigid_free_test_active, &
       rigid_two_way_test_active, ibm_beta, coupling_ramp_steps, rigid_output_interval, rigid_two_way_output_interval, &
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
    real(mytype) :: lag_total_local(3), lag_total(3), eul_total(3), abs_force_balance(3)
    real(mytype) :: sumw_min, sumw_max, spread_hy_loc_min, spread_hy_loc_max
    real(mytype) :: beta_eff, ramp_factor, coupling_step, dt_stage
    real(mytype) :: r(3), torque_local(3), omega_dot(3), torque_perp(3), pnorm
    real(mytype) :: lag_impulse(3), euler_impulse(3)
    integer :: l, npts, ierr
    logical :: output_now

    if (.not.fiber_active) return
    if (.not.rigid_two_way_test_active) return

    call update_rigid_free_geometry()
    call update_rigid_free_xdot()
    call run_fiber_interp_solver_readonly(uxe, uye, uze, itime)

    if (.not.allocated(fiber_slip)) allocate(fiber_slip(3, fiber_nl))
    if (.not.allocated(fiber_coupling_force)) allocate(fiber_coupling_force(3, fiber_nl))
    fiber_slip = 0._mytype
    fiber_coupling_force = 0._mytype

    coupling_step = real(max(1, itime - ifirst + 1), mytype)
    dt_stage = gdt(isubstep)
    if (coupling_ramp_steps > 0) then
      ramp_factor = coupling_step / real(coupling_ramp_steps, mytype)
      ramp_factor = min(1._mytype, ramp_factor)
    else
      ramp_factor = 1._mytype
    endif
    beta_eff = ibm_beta * ramp_factor

    fiber_slip = fiber_uinterp - fiber_xdot
    fiber_coupling_force = beta_eff * fiber_slip

    if (.not.allocated(fiber_euler_force_x)) allocate(fiber_euler_force_x(size(uxe,1), size(uxe,2), size(uxe,3)))
    if (.not.allocated(fiber_euler_force_y)) allocate(fiber_euler_force_y(size(uxe,1), size(uxe,2), size(uxe,3)))
    if (.not.allocated(fiber_euler_force_z)) allocate(fiber_euler_force_z(size(uxe,1), size(uxe,2), size(uxe,3)))
    fiber_euler_force_x = 0._mytype
    fiber_euler_force_y = 0._mytype
    fiber_euler_force_z = 0._mytype

    call spread_lagrangian_force_to_euler(fiber_coupling_force, fiber_euler_force_x, fiber_euler_force_y, fiber_euler_force_z, &
         eul_total, sumw_min, sumw_max, xp, yp, zp, spread_hy_loc_min, spread_hy_loc_max)

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
    fiber_force_total = lag_total
    abs_force_balance = abs(eul_total - lag_total)

    torque_perp = fiber_torque_total - dot_product(fiber_torque_total, fiber_p) * fiber_p
    omega_dot = torque_perp / fiber_inertia_perp
    fiber_uc = fiber_uc + dt_stage * (fiber_force_total / fiber_mass)
    fiber_xc = fiber_xc + dt_stage * fiber_uc
    fiber_xc(1) = wrap_periodic(fiber_xc(1), xlx)
    fiber_xc(3) = wrap_periodic(fiber_xc(3), zlz)
    fiber_omega = fiber_omega + dt_stage * omega_dot
    fiber_p = fiber_p + dt_stage * cross_product(fiber_omega, fiber_p)
    pnorm = sqrt(sum(fiber_p**2))
    if (pnorm > 0._mytype) fiber_p = fiber_p / pnorm
    p_norm_error = abs(sqrt(sum(fiber_p**2)) - 1._mytype)

    call update_rigid_free_geometry()
    call update_rigid_free_xdot()
    call compute_rigid_free_spacing_error(spacing_error_max)

    slip_max = maxval(abs(fiber_slip))
    npts = 3 * fiber_nl
    slip_rms = sqrt(sum(fiber_slip**2) / real(npts, mytype))
    lag_impulse = lag_total * dt_stage
    euler_impulse = eul_total * dt_stage

    output_now = (isubstep == nsubsteps) .and. &
         ((itime == ifirst) .or. (mod(itime, max(1, rigid_two_way_output_interval)) == 0))
    if (output_now) then
      call write_fiber_rigid_two_way_points(itime)
      call write_fiber_rigid_two_way_summary(itime, time, slip_max, slip_rms, spacing_error_max, p_norm_error, &
           lag_total, eul_total)
      call write_fiber_rigid_two_way_momentum(itime, time, lag_total, eul_total, abs_force_balance, dt_stage, &
           lag_impulse, euler_impulse)
    endif

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
