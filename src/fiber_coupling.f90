!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_coupling

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank, nproc
  use mpi, only : MPI_COMM_WORLD, MPI_ABORT
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use param, only : ifirst
  use variables, only : xp, yp, zp
  use fiber_types, only : fiber_active, fiber_nl, rigid_coupling_test_active, ibm_beta, &
       coupling_ramp_steps, rigid_output_interval, fiber_xdot, fiber_uinterp, fiber_slip, fiber_coupling_force, &
       fiber_quad_w, fiber_euler_force_x, fiber_euler_force_y, fiber_euler_force_z, rigid_motion_case
  use fiber_rigid_motion, only : update_rigid_motion, compute_rigid_spacing_error
  use fiber_interp, only : run_fiber_interp_solver_readonly
  use fiber_spread, only : spread_lagrangian_force_to_euler
  use fiber_io, only : write_fiber_rigid_coupling_points, write_fiber_rigid_coupling_summary

  implicit none

contains

  subroutine run_rigid_coupling_step(uxe, uye, uze, time, itime)

    real(mytype), intent(in), dimension(:,:,:) :: uxe, uye, uze
    real(mytype), intent(in) :: time
    integer, intent(in) :: itime

    real(mytype) :: slip_max, slip_rms, spacing_error_max
    real(mytype) :: lag_total(3), eul_total(3), abs_force_balance(3)
    real(mytype) :: sumw_min, sumw_max, spread_hy_loc_min, spread_hy_loc_max
    real(mytype) :: u_interp_max_norm, xdot_max_norm, coupling_force_max_norm, euler_force_max_norm
    real(mytype) :: beta_eff, ramp_factor
    logical :: failed_flag
    logical :: is_finite
    character(len=64) :: fail_quantity
    integer :: failure_code
    integer :: npts, coupling_step
    integer :: ierr
    logical :: output_now

    if (.not.fiber_active) return
    if (.not.rigid_coupling_test_active) return

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

    coupling_step = max(1, itime - ifirst + 1)
    if (coupling_ramp_steps > 0) then
      ramp_factor = real(coupling_step, mytype) / real(coupling_ramp_steps, mytype)
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

    output_now = (itime == ifirst) .or. (mod(itime, max(1, rigid_output_interval)) == 0)
    if (failed_flag) output_now = .true.

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

  subroutine add_rigid_coupling_rhs(dux, duy, duz)

    real(mytype), intent(inout), dimension(:,:,:) :: dux, duy, duz

    if (.not.rigid_coupling_test_active) return
    if (.not.allocated(fiber_euler_force_x)) return

    dux = dux + fiber_euler_force_x
    duy = duy + fiber_euler_force_y
    duz = duz + fiber_euler_force_z

  end subroutine add_rigid_coupling_rhs

end module fiber_coupling
