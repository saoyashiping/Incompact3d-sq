!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_coupling

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank, nproc
  use variables, only : xp, yp, zp
  use fiber_types, only : fiber_active, fiber_nl, rigid_coupling_test_active, ibm_beta, &
       rigid_output_interval, fiber_xdot, fiber_uinterp, fiber_slip, fiber_coupling_force, &
       fiber_quad_w, fiber_euler_force_x, fiber_euler_force_y, fiber_euler_force_z
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
    integer :: npts

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

    fiber_slip = fiber_uinterp - fiber_xdot
    fiber_coupling_force = ibm_beta * fiber_slip

    if (.not.allocated(fiber_euler_force_x)) allocate(fiber_euler_force_x(size(uxe,1), size(uxe,2), size(uxe,3)))
    if (.not.allocated(fiber_euler_force_y)) allocate(fiber_euler_force_y(size(uxe,1), size(uxe,2), size(uxe,3)))
    if (.not.allocated(fiber_euler_force_z)) allocate(fiber_euler_force_z(size(uxe,1), size(uxe,2), size(uxe,3)))

    call spread_lagrangian_force_to_euler(fiber_coupling_force, fiber_euler_force_x, fiber_euler_force_y, fiber_euler_force_z, &
         eul_total, sumw_min, sumw_max, xp, yp, zp, spread_hy_loc_min, spread_hy_loc_max)

    lag_total(1) = sum(fiber_coupling_force(1,:) * fiber_quad_w(:))
    lag_total(2) = sum(fiber_coupling_force(2,:) * fiber_quad_w(:))
    lag_total(3) = sum(fiber_coupling_force(3,:) * fiber_quad_w(:))

    abs_force_balance = abs(eul_total - lag_total)

    slip_max = maxval(abs(fiber_slip))
    npts = 3 * fiber_nl
    slip_rms = sqrt(sum(fiber_slip**2) / real(npts, mytype))

    call compute_rigid_spacing_error(spacing_error_max)

    if (mod(itime, max(1, rigid_output_interval)) == 0) then
      if (nrank == 0) then
        write(*,'(A,2ES12.4)') 'Rigid coupling spread sumw min/max            : ', sumw_min, sumw_max
        write(*,'(A,2ES12.4)') 'Rigid coupling spread hy_loc min/max          : ', spread_hy_loc_min, spread_hy_loc_max
      endif
      call write_fiber_rigid_coupling_points(itime)
      call write_fiber_rigid_coupling_summary(itime, time, slip_max, slip_rms, lag_total, eul_total, &
           abs_force_balance, spacing_error_max)
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
