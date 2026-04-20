!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_io

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use param, only : xlx, zlz
  use fiber_types, only : fiber_active, fiber_nl, fiber_x, fiber_uinterp, fiber_uexact, fiber_uerror, &
       fiber_sumw, fiber_test_force, fiber_quad_w, spread_test_case, fiber_xdot, fiber_slip, fiber_coupling_force, &
       fiber_xc, fiber_uc, fiber_p, fiber_omega, fiber_force_total, fiber_torque_total, &
       rigid_two_way_min_wall_gap, fiber_s_ref, fiber_x_old, fiber_x_nm1, fiber_flexible_active, &
       fiber_flex_initialized, fiber_length, fiber_ds, fiber_center, fiber_direction, &
       fiber_flex_operator_test_active, fiber_flex_bending_test_active

  implicit none

contains

  pure real(mytype) function periodic_delta_local(delta, period_length)

    real(mytype), intent(in) :: delta, period_length

    periodic_delta_local = delta
    if (period_length > 0._mytype) then
      if (periodic_delta_local >  0.5_mytype * period_length) periodic_delta_local = periodic_delta_local - period_length
      if (periodic_delta_local < -0.5_mytype * period_length) periodic_delta_local = periodic_delta_local + period_length
    endif

  end function periodic_delta_local

  subroutine write_fiber_points(filename)

    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile, l

    if (.not.fiber_active) return
    if (nrank /= 0) return

    if (.not.allocated(fiber_x)) then
      write(*,*) 'Error: fiber_x is not allocated in write_fiber_points.'
      stop
    endif

    output_file = 'fiber_points.dat'
    if (present(filename)) output_file = filename

    open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
    write(ifile,'(A)') 'index x y z'
    do l = 1, fiber_nl
      write(ifile,'(I8,1X,ES24.16,1X,ES24.16,1X,ES24.16)') l, fiber_x(1,l), fiber_x(2,l), fiber_x(3,l)
    enddo
    close(ifile)

  end subroutine write_fiber_points

  subroutine write_fiber_interp(filename)

    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile, l

    if (.not.fiber_active) return
    if (nrank /= 0) return

    if (.not.allocated(fiber_x)) then
      write(*,*) 'Error: fiber_x is not allocated in write_fiber_interp.'
      stop
    endif

    if (.not.allocated(fiber_uinterp)) then
      write(*,*) 'Error: fiber_uinterp is not allocated in write_fiber_interp.'
      stop
    endif

    output_file = 'fiber_interp.dat'
    if (present(filename)) output_file = filename

    open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
    write(ifile,'(A)') 'index x y z u_interp v_interp w_interp u_exact v_exact w_exact err_u err_v err_w sumw'
    do l = 1, fiber_nl
      write(ifile,'(I8,1X,13(ES24.16,1X))') l, fiber_x(1,l), fiber_x(2,l), fiber_x(3,l), &
           fiber_uinterp(1,l), fiber_uinterp(2,l), fiber_uinterp(3,l), &
           fiber_uexact(1,l), fiber_uexact(2,l), fiber_uexact(3,l), &
           fiber_uerror(1,l), fiber_uerror(2,l), fiber_uerror(3,l), fiber_sumw(l)
    enddo
    close(ifile)

  end subroutine write_fiber_interp


  subroutine write_fiber_interp_solver(itime, filename)

    integer, intent(in) :: itime
    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile, l

    if (.not.fiber_active) return
    if (nrank /= 0) return

    if (.not.allocated(fiber_x)) then
      write(*,*) 'Error: fiber_x is not allocated in write_fiber_interp_solver.'
      stop
    endif

    if (.not.allocated(fiber_uinterp)) then
      write(*,*) 'Error: fiber_uinterp is not allocated in write_fiber_interp_solver.'
      stop
    endif

    output_file = 'fiber_interp_solver.dat'
    if (present(filename)) output_file = filename

    open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
    write(ifile,'(A)') 'itime index x y z u_interp v_interp w_interp sumw'
    do l = 1, fiber_nl
      write(ifile,'(I10,1X,I8,1X,8(ES24.16,1X))') itime, l, fiber_x(1,l), fiber_x(2,l), fiber_x(3,l), &
           fiber_uinterp(1,l), fiber_uinterp(2,l), fiber_uinterp(3,l), fiber_sumw(l)
    enddo
    close(ifile)

  end subroutine write_fiber_interp_solver


  subroutine write_fiber_spread_lagrangian(filename)

    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile, l

    if (.not.fiber_active) return
    if (nrank /= 0) return

    output_file = 'fiber_spread_lagrangian.dat'
    if (present(filename)) output_file = filename

    open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
    write(ifile,'(A)') 'index x y z Fx Fy Fz quadrature_weight'
    do l = 1, fiber_nl
      write(ifile,'(I8,1X,7(ES24.16,1X))') l, fiber_x(1,l), fiber_x(2,l), fiber_x(3,l), &
           fiber_test_force(1,l), fiber_test_force(2,l), fiber_test_force(3,l), fiber_quad_w(l)
    enddo
    close(ifile)

  end subroutine write_fiber_spread_lagrangian

  subroutine write_fiber_spread_summary(lag_total, eul_total, abs_error, rel_error, sumw_min, sumw_max, filename)

    real(mytype), intent(in), dimension(3) :: lag_total, eul_total, abs_error, rel_error
    real(mytype), intent(in) :: sumw_min, sumw_max
    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile

    if (nrank /= 0) return

    output_file = 'fiber_spread_summary.dat'
    if (present(filename)) output_file = filename

    open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
    write(ifile,'(A)') 'test_case lag_Fx lag_Fy lag_Fz eul_Fx eul_Fy eul_Fz abs_err_Fx abs_err_Fy abs_err_Fz ' // &
         'rel_err_Fx rel_err_Fy rel_err_Fz spread_sumw_min spread_sumw_max'
    write(ifile,'(I8,1X,14(ES24.16,1X))') spread_test_case, lag_total(1), lag_total(2), lag_total(3), &
         eul_total(1), eul_total(2), eul_total(3), abs_error(1), abs_error(2), abs_error(3), &
         rel_error(1), rel_error(2), rel_error(3), sumw_min, sumw_max
    close(ifile)

  end subroutine write_fiber_spread_summary

  subroutine write_fiber_rigid_coupling_points(itime, filename)

    integer, intent(in) :: itime
    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile, l
    logical :: file_exists

    if (.not.fiber_active) return
    if (nrank /= 0) return

    output_file = 'fiber_rigid_coupling_points.dat'
    if (present(filename)) output_file = filename

    inquire(file=trim(output_file), exist=file_exists)
    if (file_exists) then
      open(newunit=ifile, file=trim(output_file), status='old', action='write', position='append', form='formatted')
    else
      open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
      write(ifile,'(A)') 'itime index x y z xdot_x xdot_y xdot_z u_interp v_interp w_interp slip_x slip_y slip_z Ffs_x Ffs_y Ffs_z'
    endif

    do l = 1, fiber_nl
      write(ifile,'(I10,1X,I8,1X,15(ES24.16,1X))') itime, l, fiber_x(1,l), fiber_x(2,l), fiber_x(3,l), &
           fiber_xdot(1,l), fiber_xdot(2,l), fiber_xdot(3,l), &
           fiber_uinterp(1,l), fiber_uinterp(2,l), fiber_uinterp(3,l), &
           fiber_slip(1,l), fiber_slip(2,l), fiber_slip(3,l), &
           fiber_coupling_force(1,l), fiber_coupling_force(2,l), fiber_coupling_force(3,l)
    enddo
    close(ifile)

  end subroutine write_fiber_rigid_coupling_points

  subroutine write_fiber_rigid_coupling_summary(itime, time, slip_max, slip_rms, lag_total, eul_total, &
       abs_force_balance, spacing_error_max, u_interp_max_norm, xdot_max_norm, coupling_force_max_norm, &
       euler_force_max_norm, beta_input, beta_eff, ramp_factor, failed_flag, failure_code, filename)

    integer, intent(in) :: itime
    real(mytype), intent(in) :: time, slip_max, slip_rms, spacing_error_max
    real(mytype), intent(in) :: u_interp_max_norm, xdot_max_norm, coupling_force_max_norm, euler_force_max_norm
    real(mytype), intent(in) :: beta_input, beta_eff, ramp_factor
    real(mytype), intent(in), dimension(3) :: lag_total, eul_total, abs_force_balance
    logical, intent(in) :: failed_flag
    integer, intent(in) :: failure_code
    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile
    logical :: file_exists

    if (nrank /= 0) return

    output_file = 'fiber_rigid_coupling_summary.dat'
    if (present(filename)) output_file = filename

    inquire(file=trim(output_file), exist=file_exists)
    if (file_exists) then
      open(newunit=ifile, file=trim(output_file), status='old', action='write', position='append', form='formatted')
    else
      open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
      write(ifile,'(A)') 'itime time slip_max slip_rms lag_total_Fx lag_total_Fy lag_total_Fz euler_total_Fx ' // &
           'euler_total_Fy euler_total_Fz abs_force_balance_Fx abs_force_balance_Fy abs_force_balance_Fz spacing_error_max ' // &
           'u_interp_max_norm xdot_max_norm coupling_force_max_norm euler_force_max_norm beta_input beta_eff ' // &
           'ramp_factor failed_flag failure_code'
    endif

    write(ifile,'(I10,1X,21(ES24.16,1X),I8)') itime, time, slip_max, slip_rms, &
         lag_total(1), lag_total(2), lag_total(3), eul_total(1), eul_total(2), eul_total(3), &
         abs_force_balance(1), abs_force_balance(2), abs_force_balance(3), spacing_error_max, &
         u_interp_max_norm, xdot_max_norm, coupling_force_max_norm, euler_force_max_norm, &
         beta_input, beta_eff, ramp_factor, merge(1._mytype, 0._mytype, failed_flag), failure_code
    close(ifile)

  end subroutine write_fiber_rigid_coupling_summary

  subroutine write_fiber_rigid_free_points(itime, filename)

    integer, intent(in) :: itime
    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile, l
    logical :: file_exists

    if (.not.fiber_active) return
    if (nrank /= 0) return

    output_file = 'fiber_rigid_free_points.dat'
    if (present(filename)) output_file = filename

    inquire(file=trim(output_file), exist=file_exists)
    if (file_exists) then
      open(newunit=ifile, file=trim(output_file), status='old', action='write', position='append', form='formatted')
    else
      open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
      write(ifile,'(A)') 'itime index x y z u_interp v_interp w_interp xdot_x xdot_y xdot_z ' // &
           'slip_x slip_y slip_z Ffs_x Ffs_y Ffs_z'
    endif

    do l = 1, fiber_nl
      write(ifile,'(I10,1X,I8,1X,15(ES24.16,1X))') itime, l, fiber_x(1,l), fiber_x(2,l), fiber_x(3,l), &
           fiber_uinterp(1,l), fiber_uinterp(2,l), fiber_uinterp(3,l), &
           fiber_xdot(1,l), fiber_xdot(2,l), fiber_xdot(3,l), &
           fiber_slip(1,l), fiber_slip(2,l), fiber_slip(3,l), &
           fiber_coupling_force(1,l), fiber_coupling_force(2,l), fiber_coupling_force(3,l)
    enddo
    close(ifile)

  end subroutine write_fiber_rigid_free_points

  subroutine write_fiber_rigid_free_summary(itime, time, xc, uc, p, omega, total_force, total_torque, &
       slip_max, slip_rms, spacing_error_max, p_norm_error, failed_flag, failure_code, filename)

    integer, intent(in) :: itime
    real(mytype), intent(in) :: time, slip_max, slip_rms, spacing_error_max, p_norm_error
    real(mytype), intent(in), dimension(3) :: xc, uc, p, omega, total_force, total_torque
    logical, intent(in) :: failed_flag
    integer, intent(in) :: failure_code
    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile
    logical :: file_exists

    if (nrank /= 0) return

    output_file = 'fiber_rigid_free_summary.dat'
    if (present(filename)) output_file = filename

    inquire(file=trim(output_file), exist=file_exists)
    if (file_exists) then
      open(newunit=ifile, file=trim(output_file), status='old', action='write', position='append', form='formatted')
    else
      open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
      write(ifile,'(A)') 'itime time xc yc zc uc_x uc_y uc_z p_x p_y p_z omega_x omega_y omega_z ' // &
           'total_Fx total_Fy total_Fz total_Mx total_My total_Mz slip_max slip_rms spacing_error_max p_norm_error ' // &
           'failed_flag failure_code'
    endif

    write(ifile,'(I10,1X,24(ES24.16,1X),I8)') itime, time, xc(1), xc(2), xc(3), uc(1), uc(2), uc(3), &
         p(1), p(2), p(3), omega(1), omega(2), omega(3), total_force(1), total_force(2), total_force(3), &
         total_torque(1), total_torque(2), total_torque(3), slip_max, slip_rms, spacing_error_max, p_norm_error, &
         merge(1._mytype, 0._mytype, failed_flag), failure_code
    close(ifile)

  end subroutine write_fiber_rigid_free_summary

  subroutine write_fiber_rigid_two_way_points(itime, filename)

    integer, intent(in) :: itime
    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile, l
    logical :: file_exists

    if (.not.fiber_active) return
    if (nrank /= 0) return

    output_file = 'rigid_two_way_points.dat'
    if (present(filename)) output_file = filename

    inquire(file=trim(output_file), exist=file_exists)
    if (file_exists) then
      open(newunit=ifile, file=trim(output_file), status='old', action='write', position='append', form='formatted')
    else
      open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
      write(ifile,'(A)') 'itime index x y z xdot_x xdot_y xdot_z u_interp v_interp w_interp slip_x slip_y slip_z Ffs_x Ffs_y Ffs_z'
    endif

    do l = 1, fiber_nl
      write(ifile,'(I10,1X,I8,1X,15(ES24.16,1X))') itime, l, fiber_x(1,l), fiber_x(2,l), fiber_x(3,l), &
           fiber_xdot(1,l), fiber_xdot(2,l), fiber_xdot(3,l), &
           fiber_uinterp(1,l), fiber_uinterp(2,l), fiber_uinterp(3,l), &
           fiber_slip(1,l), fiber_slip(2,l), fiber_slip(3,l), &
           fiber_coupling_force(1,l), fiber_coupling_force(2,l), fiber_coupling_force(3,l)
    enddo
    close(ifile)

  end subroutine write_fiber_rigid_two_way_points

  subroutine write_fiber_rigid_two_way_summary(itime, time, slip_max, slip_rms, spacing_error_max, p_norm_error, &
       lag_total, eul_total, filename)

    integer, intent(in) :: itime
    real(mytype), intent(in) :: time, slip_max, slip_rms, spacing_error_max, p_norm_error
    real(mytype), intent(in) :: lag_total(3), eul_total(3)
    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile
    logical :: file_exists

    if (nrank /= 0) return

    output_file = 'rigid_two_way_summary.dat'
    if (present(filename)) output_file = filename

    inquire(file=trim(output_file), exist=file_exists)
    if (file_exists) then
      open(newunit=ifile, file=trim(output_file), status='old', action='write', position='append', form='formatted')
    else
      open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
      write(ifile,'(A)') 'itime time xc yc zc uc_x uc_y uc_z p_x p_y p_z omega_x omega_y omega_z ' // &
           'total_Fx total_Fy total_Fz total_Mx total_My total_Mz lag_total_Fx lag_total_Fy lag_total_Fz ' // &
           'eul_total_Fx eul_total_Fy eul_total_Fz slip_max slip_rms spacing_error_max p_norm_error'
    endif

    write(ifile,'(I10,1X,30(ES24.16,1X))') itime, time, fiber_xc(1), fiber_xc(2), fiber_xc(3), &
         fiber_uc(1), fiber_uc(2), fiber_uc(3), fiber_p(1), fiber_p(2), fiber_p(3), fiber_omega(1), fiber_omega(2), &
         fiber_omega(3), fiber_force_total(1), fiber_force_total(2), fiber_force_total(3), fiber_torque_total(1), &
         fiber_torque_total(2), fiber_torque_total(3), lag_total(1), lag_total(2), lag_total(3), eul_total(1), &
         eul_total(2), eul_total(3), slip_max, slip_rms, spacing_error_max, p_norm_error
    close(ifile)

  end subroutine write_fiber_rigid_two_way_summary

  subroutine write_fiber_rigid_two_way_momentum(itime, time, lag_total, eul_total, abs_force_balance, dt_stage, &
       lag_impulse, euler_impulse, filename)

    integer, intent(in) :: itime
    real(mytype), intent(in) :: time, dt_stage
    real(mytype), intent(in) :: lag_total(3), eul_total(3), abs_force_balance(3), lag_impulse(3), euler_impulse(3)
    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile
    logical :: file_exists

    if (nrank /= 0) return

    output_file = 'rigid_two_way_momentum_check.dat'
    if (present(filename)) output_file = filename

    inquire(file=trim(output_file), exist=file_exists)
    if (file_exists) then
      open(newunit=ifile, file=trim(output_file), status='old', action='write', position='append', form='formatted')
    else
      open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
      write(ifile,'(A)') 'itime time lag_total_force_x lag_total_force_y lag_total_force_z ' // &
           'euler_total_force_x euler_total_force_y euler_total_force_z ' // &
           'abs_force_balance_x abs_force_balance_y abs_force_balance_z dt_stage ' // &
           'lag_impulse_x lag_impulse_y lag_impulse_z ' // &
           'euler_feedback_impulse_x euler_feedback_impulse_y euler_feedback_impulse_z'
    endif

    write(ifile,'(I10,1X,18(ES24.16,1X))') itime, time, lag_total(1), lag_total(2), lag_total(3), eul_total(1), &
         eul_total(2), eul_total(3), abs_force_balance(1), abs_force_balance(2), abs_force_balance(3), dt_stage, &
         lag_impulse(1), lag_impulse(2), lag_impulse(3), euler_impulse(1), euler_impulse(2), euler_impulse(3)
    close(ifile)

  end subroutine write_fiber_rigid_two_way_momentum

  subroutine write_fiber_rigid_two_way_turb_series(itime, time, slip_max, slip_rms, min_wall_gap, filename)

    integer, intent(in) :: itime
    real(mytype), intent(in) :: time, slip_max, slip_rms, min_wall_gap
    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile
    logical :: file_exists

    if (nrank /= 0) return

    output_file = 'rigid_two_way_turb_series.dat'
    if (present(filename)) output_file = filename

    inquire(file=trim(output_file), exist=file_exists)
    if (file_exists) then
      open(newunit=ifile, file=trim(output_file), status='old', action='write', position='append', form='formatted')
    else
      open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
      write(ifile,'(A)') 'itime time xc yc zc uc_x uc_y uc_z p_x p_y p_z omega_x omega_y omega_z ' // &
           'slip_max slip_rms min_wall_gap min_wall_gap_limit'
    endif

    write(ifile,'(I10,1X,17(ES24.16,1X))') itime, time, fiber_xc(1), fiber_xc(2), fiber_xc(3), &
         fiber_uc(1), fiber_uc(2), fiber_uc(3), fiber_p(1), fiber_p(2), fiber_p(3), &
         fiber_omega(1), fiber_omega(2), fiber_omega(3), &
         slip_max, slip_rms, min_wall_gap, rigid_two_way_min_wall_gap
    close(ifile)

  end subroutine write_fiber_rigid_two_way_turb_series

  subroutine write_fiber_flex_init_points(filename)

    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile, l

    if (.not.fiber_active) return
    if (.not.fiber_flexible_active) return
    if (nrank /= 0) return

    if (.not.allocated(fiber_x)) then
      write(*,*) 'Error: fiber_x is not allocated in write_fiber_flex_init_points.'
      stop
    endif
    if (.not.allocated(fiber_s_ref)) then
      write(*,*) 'Error: fiber_s_ref is not allocated in write_fiber_flex_init_points.'
      stop
    endif
    if (.not.allocated(fiber_x_old)) then
      write(*,*) 'Error: fiber_x_old is not allocated in write_fiber_flex_init_points.'
      stop
    endif
    if (.not.allocated(fiber_x_nm1)) then
      write(*,*) 'Error: fiber_x_nm1 is not allocated in write_fiber_flex_init_points.'
      stop
    endif

    output_file = 'fiber_flex_init_points.dat'
    if (present(filename)) output_file = filename

    open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
    write(ifile,'(A)') '# Flexible fiber init points (Step 3.1 data layer only)'
    write(ifile,'(A)') '# l s_ref x y z x_old y_old z_old x_nm1 y_nm1 z_nm1'
    do l = 1, fiber_nl
      write(ifile,'(I8,1X,10(ES24.16,1X))') l, fiber_s_ref(l), &
           fiber_x(1,l), fiber_x(2,l), fiber_x(3,l), &
           fiber_x_old(1,l), fiber_x_old(2,l), fiber_x_old(3,l), &
           fiber_x_nm1(1,l), fiber_x_nm1(2,l), fiber_x_nm1(3,l)
    enddo
    close(ifile)

  end subroutine write_fiber_flex_init_points

  subroutine write_fiber_flex_init_summary(filename)

    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile, l
    real(mytype) :: dx, dy, dz, seg_len, spacing_error_max, end_to_end_length

    if (.not.fiber_active) return
    if (.not.fiber_flexible_active) return
    if (nrank /= 0) return

    if (.not.allocated(fiber_x)) then
      write(*,*) 'Error: fiber_x is not allocated in write_fiber_flex_init_summary.'
      stop
    endif

    spacing_error_max = 0._mytype
    do l = 1, fiber_nl - 1
      dx = periodic_delta_local(fiber_x(1,l+1) - fiber_x(1,l), xlx)
      dy = fiber_x(2,l+1) - fiber_x(2,l)
      dz = periodic_delta_local(fiber_x(3,l+1) - fiber_x(3,l), zlz)
      seg_len = sqrt(dx*dx + dy*dy + dz*dz)
      spacing_error_max = max(spacing_error_max, abs(seg_len - fiber_ds))
    enddo

    dx = periodic_delta_local(fiber_x(1,fiber_nl) - fiber_x(1,1), xlx)
    dy = fiber_x(2,fiber_nl) - fiber_x(2,1)
    dz = periodic_delta_local(fiber_x(3,fiber_nl) - fiber_x(3,1), zlz)
    end_to_end_length = sqrt(dx*dx + dy*dy + dz*dz)

    output_file = 'fiber_flex_init_summary.dat'
    if (present(filename)) output_file = filename

    open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
    write(ifile,'(A)') '# Flexible fiber init summary (Step 3.1 data layer only)'
    write(ifile,'(A,I8)') 'fiber_nl ', fiber_nl
    write(ifile,'(A,ES24.16)') 'fiber_length ', fiber_length
    write(ifile,'(A,ES24.16)') 'fiber_ds ', fiber_ds
    write(ifile,'(A,3(ES24.16,1X))') 'fiber_center ', fiber_center(1), fiber_center(2), fiber_center(3)
    write(ifile,'(A,3(ES24.16,1X))') 'fiber_direction ', fiber_direction(1), fiber_direction(2), fiber_direction(3)
    write(ifile,'(A,ES24.16)') 'end_to_end_length ', end_to_end_length
    write(ifile,'(A,ES24.16)') 'spacing_error_max ', spacing_error_max
    write(ifile,'(A,L1)') 'fiber_flex_initialized ', fiber_flex_initialized
    close(ifile)

  end subroutine write_fiber_flex_init_summary

  subroutine write_fiber_flex_operator_points(x_exact, xs_num, xs_exact, xss_num, xss_exact, xssss_num, xssss_exact, filename)

    real(mytype), intent(in) :: x_exact(3, fiber_nl), xs_num(3, fiber_nl), xs_exact(3, fiber_nl)
    real(mytype), intent(in) :: xss_num(3, fiber_nl), xss_exact(3, fiber_nl)
    real(mytype), intent(in) :: xssss_num(3, fiber_nl), xssss_exact(3, fiber_nl)
    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile, l

    if (.not.fiber_active) return
    if (.not.fiber_flexible_active) return
    if (.not.fiber_flex_operator_test_active) return
    if (nrank /= 0) return

    output_file = 'fiber_flex_operator_points.dat'
    if (present(filename)) output_file = filename

    open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
    write(ifile,'(A)') '# l s_ref x_exact_x x_exact_y x_exact_z ' // &
         'xs_num_x xs_num_y xs_num_z xs_exact_x xs_exact_y xs_exact_z ' // &
         'xss_num_x xss_num_y xss_num_z xss_exact_x xss_exact_y xss_exact_z ' // &
         'xssss_num_x xssss_num_y xssss_num_z xssss_exact_x xssss_exact_y xssss_exact_z'
    do l = 1, fiber_nl
      write(ifile,'(I8,1X,22(ES24.16,1X))') l, fiber_s_ref(l), &
           x_exact(1,l), x_exact(2,l), x_exact(3,l), &
           xs_num(1,l), xs_num(2,l), xs_num(3,l), xs_exact(1,l), xs_exact(2,l), xs_exact(3,l), &
           xss_num(1,l), xss_num(2,l), xss_num(3,l), xss_exact(1,l), xss_exact(2,l), xss_exact(3,l), &
           xssss_num(1,l), xssss_num(2,l), xssss_num(3,l), xssss_exact(1,l), xssss_exact(2,l), xssss_exact(3,l)
    enddo
    close(ifile)

  end subroutine write_fiber_flex_operator_points

  subroutine write_fiber_flex_operator_summary(flex_case, err_xs_max, err_xss_max, err_xssss_max, &
       bc_moment_left_max, bc_moment_right_max, bc_shear_left_max, bc_shear_right_max, filename)

    integer, intent(in) :: flex_case
    real(mytype), intent(in) :: err_xs_max, err_xss_max, err_xssss_max
    real(mytype), intent(in) :: bc_moment_left_max, bc_moment_right_max, bc_shear_left_max, bc_shear_right_max
    character(len=*), intent(in), optional :: filename

    character(len=256) :: output_file
    integer :: ifile

    if (.not.fiber_active) return
    if (.not.fiber_flexible_active) return
    if (.not.fiber_flex_operator_test_active) return
    if (nrank /= 0) return

    output_file = 'fiber_flex_operator_summary.dat'
    if (present(filename)) output_file = filename

    open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
    write(ifile,'(A)') '# Flexible operator summary (Step 3.2 spatial operator verification only)'
    write(ifile,'(A,I8)') 'fiber_nl ', fiber_nl
    write(ifile,'(A,ES24.16)') 'fiber_ds ', fiber_ds
    write(ifile,'(A,I8)') 'fiber_flex_operator_case ', flex_case
    write(ifile,'(A,ES24.16)') 'err_xs_max ', err_xs_max
    write(ifile,'(A,ES24.16)') 'err_xss_max ', err_xss_max
    write(ifile,'(A,ES24.16)') 'err_xssss_max ', err_xssss_max
    write(ifile,'(A,ES24.16)') 'bc_moment_left_max ', bc_moment_left_max
    write(ifile,'(A,ES24.16)') 'bc_moment_right_max ', bc_moment_right_max
    write(ifile,'(A,ES24.16)') 'bc_shear_left_max ', bc_shear_left_max
    write(ifile,'(A,ES24.16)') 'bc_shear_right_max ', bc_shear_right_max
    close(ifile)

  end subroutine write_fiber_flex_operator_summary

  subroutine write_fiber_flex_bending_initial_points(x0, filename)
    real(mytype), intent(in) :: x0(3, fiber_nl)
    character(len=*), intent(in), optional :: filename
    character(len=256) :: output_file
    integer :: ifile, l

    if (nrank /= 0) return
    if (.not.fiber_flex_bending_test_active) return
    output_file = 'fiber_flex_bending_initial_points.dat'
    if (present(filename)) output_file = filename
    open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
    write(ifile,'(A)') '# l s_ref x0 y0 z0'
    do l = 1, fiber_nl
      write(ifile,'(I8,1X,4(ES24.16,1X))') l, fiber_s_ref(l), x0(1,l), x0(2,l), x0(3,l)
    enddo
    close(ifile)
  end subroutine write_fiber_flex_bending_initial_points

  subroutine write_fiber_flex_bending_final_points(xf, filename)
    real(mytype), intent(in) :: xf(3, fiber_nl)
    character(len=*), intent(in), optional :: filename
    character(len=256) :: output_file
    integer :: ifile, l

    if (nrank /= 0) return
    if (.not.fiber_flex_bending_test_active) return
    output_file = 'fiber_flex_bending_final_points.dat'
    if (present(filename)) output_file = filename
    open(newunit=ifile, file=trim(output_file), status='replace', action='write', form='formatted')
    write(ifile,'(A)') '# l s_ref x_final y_final z_final'
    do l = 1, fiber_nl
      write(ifile,'(I8,1X,4(ES24.16,1X))') l, fiber_s_ref(l), xf(1,l), xf(2,l), xf(3,l)
    enddo
    close(ifile)
  end subroutine write_fiber_flex_bending_final_points

  subroutine write_fiber_flex_bending_series(step, time, bending_energy, max_update_norm, overwrite)
    integer, intent(in) :: step
    real(mytype), intent(in) :: time, bending_energy, max_update_norm
    logical, intent(in) :: overwrite
    integer :: ifile

    if (nrank /= 0) return
    if (.not.fiber_flex_bending_test_active) return
    if (overwrite) then
      open(newunit=ifile, file='fiber_flex_bending_series.dat', status='replace', action='write', form='formatted')
      write(ifile,'(A)') 'step time bending_energy max_update_norm'
    else
      open(newunit=ifile, file='fiber_flex_bending_series.dat', status='old', action='write', position='append', form='formatted')
    endif
    write(ifile,'(I10,1X,3(ES24.16,1X))') step, time, bending_energy, max_update_norm
    close(ifile)
  end subroutine write_fiber_flex_bending_series

  subroutine write_fiber_flex_bending_summary(case_id, dt_b, gamma_b, nsteps, e_init, e_final, energy_monotone, &
       max_update_last_step, straight_preservation_error_max, final_displacement_norm)
    integer, intent(in) :: case_id, nsteps
    real(mytype), intent(in) :: dt_b, gamma_b, e_init, e_final, max_update_last_step
    real(mytype), intent(in) :: straight_preservation_error_max, final_displacement_norm
    logical, intent(in) :: energy_monotone
    integer :: ifile

    if (nrank /= 0) return
    if (.not.fiber_flex_bending_test_active) return
    open(newunit=ifile, file='fiber_flex_bending_summary.dat', status='replace', action='write', form='formatted')
    write(ifile,'(A)') '# Step 3.3 bending-only implicit summary'
    write(ifile,'(A,I8)') 'fiber_nl ', fiber_nl
    write(ifile,'(A,ES24.16)') 'fiber_ds ', fiber_ds
    write(ifile,'(A,I8)') 'fiber_flex_bending_case ', case_id
    write(ifile,'(A,ES24.16)') 'fiber_flex_bending_dt ', dt_b
    write(ifile,'(A,ES24.16)') 'fiber_bending_gamma ', gamma_b
    write(ifile,'(A,I8)') 'fiber_flex_bending_nsteps ', nsteps
    write(ifile,'(A,ES24.16)') 'initial_bending_energy ', e_init
    write(ifile,'(A,ES24.16)') 'final_bending_energy ', e_final
    write(ifile,'(A,L1)') 'energy_monotone_decay ', energy_monotone
    write(ifile,'(A,ES24.16)') 'max_update_last_step ', max_update_last_step
    if (case_id == 1) then
      write(ifile,'(A,ES24.16)') 'straight_preservation_error_max ', straight_preservation_error_max
    else if (case_id == 2) then
      write(ifile,'(A,ES24.16)') 'final_displacement_norm ', final_displacement_norm
    endif
    close(ifile)
  end subroutine write_fiber_flex_bending_summary

end module fiber_io
