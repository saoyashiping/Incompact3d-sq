!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_rigid_free

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use mpi, only : MPI_COMM_WORLD, MPI_ABORT, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_SUM
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use param, only : ifirst, dt, xlx, zlz
  use variables, only : xp, yp, zp
  use fiber_types, only : fiber_active, fiber_nl, rigid_free_test_active, rigid_free_case, ibm_beta, &
       coupling_ramp_steps, free_output_interval, fiber_mass, fiber_inertia_perp, fiber_x, fiber_s_ref, &
       fiber_xc, fiber_uc, fiber_p, fiber_omega, fiber_force_total, fiber_torque_total, fiber_xdot, fiber_uinterp, &
       fiber_slip, fiber_coupling_force, fiber_quad_w, fiber_euler_force_x, fiber_euler_force_y, fiber_euler_force_z
  use fiber_interp, only : run_fiber_interp_solver_readonly
  use fiber_spread, only : spread_lagrangian_force_to_euler
  use fiber_io, only : write_fiber_rigid_free_points, write_fiber_rigid_free_summary

  implicit none

contains

  pure real(mytype) function wrap_periodic(value, period_length)

    real(mytype), intent(in) :: value, period_length

    wrap_periodic = value
    if (period_length > 0._mytype) then
      wrap_periodic = modulo(value, period_length)
      if (wrap_periodic < 0._mytype) wrap_periodic = wrap_periodic + period_length
    endif

  end function wrap_periodic

  pure function cross_product(a, b) result(c)

    real(mytype), intent(in) :: a(3), b(3)
    real(mytype) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

  end function cross_product

  subroutine init_rigid_free_state()

    real(mytype) :: p_norm

    if (.not.fiber_active) return
    if (.not.rigid_free_test_active) return

    if (.not.allocated(fiber_s_ref)) then
      if (nrank == 0) write(*,*) 'Error: fiber_s_ref is not allocated in init_rigid_free_state.'
      stop
    endif

    if (fiber_mass <= 0._mytype) then
      if (nrank == 0) write(*,*) 'Error: fiber_mass must be > 0 for rigid free test.'
      stop
    endif

    if (fiber_inertia_perp <= 0._mytype) then
      if (nrank == 0) write(*,*) 'Error: fiber_inertia_perp must be > 0 for rigid free test.'
      stop
    endif

    p_norm = sqrt(sum(fiber_p**2))
    if (p_norm <= 0._mytype) then
      if (nrank == 0) write(*,*) 'Error: fiber_p must be non-zero for rigid free test.'
      stop
    endif
    fiber_p = fiber_p / p_norm

    call update_rigid_free_geometry()

  end subroutine init_rigid_free_state

  subroutine update_rigid_free_geometry()

    integer :: l

    do l = 1, fiber_nl
      fiber_x(1,l) = wrap_periodic(fiber_xc(1) + fiber_s_ref(l) * fiber_p(1), xlx)
      fiber_x(2,l) = fiber_xc(2) + fiber_s_ref(l) * fiber_p(2)
      fiber_x(3,l) = wrap_periodic(fiber_xc(3) + fiber_s_ref(l) * fiber_p(3), zlz)
    enddo

  end subroutine update_rigid_free_geometry

  subroutine update_rigid_free_xdot()

    integer :: l
    real(mytype) :: r(3)

    if (.not.allocated(fiber_xdot)) allocate(fiber_xdot(3, fiber_nl))

    do l = 1, fiber_nl
      r = fiber_x(:,l) - fiber_xc
      fiber_xdot(:,l) = fiber_uc + cross_product(fiber_omega, r)
    enddo

  end subroutine update_rigid_free_xdot

  subroutine compute_rigid_free_spacing_error(spacing_error_max)

    real(mytype), intent(out) :: spacing_error_max

    integer :: l
    real(mytype) :: ds_ref, ds_now

    spacing_error_max = 0._mytype

    if (fiber_nl < 2) return
    ds_ref = abs(fiber_s_ref(2) - fiber_s_ref(1))

    do l = 1, fiber_nl - 1
      ds_now = sqrt(sum((fiber_x(:,l+1) - fiber_x(:,l))**2))
      spacing_error_max = max(spacing_error_max, abs(ds_now - ds_ref))
    enddo

  end subroutine compute_rigid_free_spacing_error

  subroutine run_rigid_free_step(uxe, uye, uze, time, itime)

    real(mytype), intent(in), dimension(:,:,:) :: uxe, uye, uze
    real(mytype), intent(in) :: time
    integer, intent(in) :: itime

    real(mytype) :: slip_max, slip_rms, spacing_error_max, p_norm_error
    real(mytype) :: lag_total_local(3), lag_total(3), eul_total(3), abs_force_balance(3)
    real(mytype) :: sumw_min, sumw_max, spread_hy_loc_min, spread_hy_loc_max
    real(mytype) :: u_interp_max_norm, xdot_max_norm, coupling_force_max_norm, euler_force_max_norm
    real(mytype) :: beta_eff, ramp_factor, coupling_step
    real(mytype) :: r(3), torque_local(3), omega_dot(3), torque_perp(3), pnorm
    real(mytype) :: omega_norm
    logical :: failed_flag
    logical :: is_finite
    character(len=64) :: fail_quantity
    integer :: failure_code
    integer :: npts
    integer :: ierr
    logical :: output_now

    if (.not.fiber_active) return
    if (.not.rigid_free_test_active) return

    call update_rigid_free_geometry()
    call update_rigid_free_xdot()
    call run_fiber_interp_solver_readonly(uxe, uye, uze, itime)

    if (.not.allocated(fiber_slip)) allocate(fiber_slip(3, fiber_nl))
    if (.not.allocated(fiber_coupling_force)) allocate(fiber_coupling_force(3, fiber_nl))
    fiber_slip = 0._mytype
    fiber_coupling_force = 0._mytype

    coupling_step = real(max(1, itime - ifirst + 1), mytype)
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
    do npts = 1, fiber_nl
      lag_total_local = lag_total_local + fiber_coupling_force(:,npts) * fiber_quad_w(npts)
      r = fiber_x(:,npts) - fiber_xc
      torque_local = torque_local + cross_product(r, fiber_coupling_force(:,npts)) * fiber_quad_w(npts)
    enddo

    call MPI_ALLREDUCE(lag_total_local, lag_total, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(torque_local, fiber_torque_total, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    fiber_force_total = lag_total

    abs_force_balance = abs(eul_total - lag_total)

    torque_perp = fiber_torque_total - dot_product(fiber_torque_total, fiber_p) * fiber_p
    omega_dot = torque_perp / fiber_inertia_perp
    fiber_uc = fiber_uc + dt * (fiber_force_total / fiber_mass)
    fiber_xc = fiber_xc + dt * fiber_uc
    fiber_omega = fiber_omega + dt * omega_dot
    fiber_p = fiber_p + dt * cross_product(fiber_omega, fiber_p)
    pnorm = sqrt(sum(fiber_p**2))
    if (pnorm > 0._mytype) fiber_p = fiber_p / pnorm
    p_norm_error = abs(sqrt(sum(fiber_p**2)) - 1._mytype)

    call update_rigid_free_geometry()
    call update_rigid_free_xdot()
    call compute_rigid_free_spacing_error(spacing_error_max)

    slip_max = maxval(abs(fiber_slip))
    npts = 3 * fiber_nl
    slip_rms = sqrt(sum(fiber_slip**2) / real(npts, mytype))
    u_interp_max_norm = maxval(sqrt(fiber_uinterp(1,:)**2 + fiber_uinterp(2,:)**2 + fiber_uinterp(3,:)**2))
    xdot_max_norm = maxval(sqrt(fiber_xdot(1,:)**2 + fiber_xdot(2,:)**2 + fiber_xdot(3,:)**2))
    coupling_force_max_norm = maxval(sqrt(fiber_coupling_force(1,:)**2 + fiber_coupling_force(2,:)**2 + &
         fiber_coupling_force(3,:)**2))
    euler_force_max_norm = maxval(sqrt(fiber_euler_force_x**2 + fiber_euler_force_y**2 + fiber_euler_force_z**2))
    omega_norm = sqrt(sum(fiber_omega**2))

    failed_flag = .false.
    fail_quantity = 'none'
    failure_code = 0

    is_finite = all(ieee_is_finite(fiber_uc))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      fail_quantity = 'fiber_uc'
      failure_code = 1
    endif
    is_finite = all(ieee_is_finite(fiber_omega))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      fail_quantity = 'fiber_omega'
      failure_code = 2
    endif
    is_finite = all(ieee_is_finite(fiber_xc))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      fail_quantity = 'fiber_xc'
      failure_code = 3
    endif
    is_finite = all(ieee_is_finite(fiber_p))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      fail_quantity = 'fiber_p'
      failure_code = 4
    endif
    is_finite = all(ieee_is_finite(fiber_force_total))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      fail_quantity = 'fiber_force_total'
      failure_code = 5
    endif
    is_finite = all(ieee_is_finite(fiber_torque_total))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      fail_quantity = 'fiber_torque_total'
      failure_code = 6
    endif
    is_finite = all(ieee_is_finite((/slip_max, slip_rms, spacing_error_max, p_norm_error, u_interp_max_norm, xdot_max_norm, &
         coupling_force_max_norm, euler_force_max_norm, omega_norm/)))
    if (.not.is_finite .and. .not.failed_flag) then
      failed_flag = .true.
      fail_quantity = 'rigid_free_summary_scalars'
      failure_code = 7
    endif

    output_now = (itime == ifirst) .or. (mod(itime, max(1, free_output_interval)) == 0)
    if (failed_flag) output_now = .true.

    if (output_now) then
      call write_fiber_rigid_free_points(itime)
      call write_fiber_rigid_free_summary(itime, time, fiber_xc, fiber_uc, fiber_p, fiber_omega, fiber_force_total, &
           fiber_torque_total, slip_max, slip_rms, spacing_error_max, p_norm_error, failed_flag, failure_code)
      if (nrank == 0) then
        write(*,'(A,3ES12.4)') 'Rigid free total force                      : ', fiber_force_total
        write(*,'(A,3ES12.4)') 'Rigid free total torque                     : ', fiber_torque_total
      endif
    endif

    if (failed_flag) then
      if (nrank == 0) then
        write(*,'(A)') 'RIGID FREE TEST FAILED'
        write(*,'(A,I6)') 'failure_code                                   : ', failure_code
        write(*,'(A,I10)') 'Failure at itime                               : ', itime
        write(*,'(A,ES12.4)') 'Failure time                                   : ', time
        write(*,'(A,I6)') 'Failure free case                              : ', rigid_free_case
        write(*,'(A,A)') 'first_nonfinite_quantity                       : ', trim(fail_quantity)
      endif
      call MPI_ABORT(MPI_COMM_WORLD, 912, ierr)
    endif

  end subroutine run_rigid_free_step

end module fiber_rigid_free
