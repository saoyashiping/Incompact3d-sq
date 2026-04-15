!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_io

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use fiber_types, only : fiber_active, fiber_nl, fiber_x, fiber_uinterp, fiber_uexact, fiber_uerror, &
       fiber_sumw, fiber_test_force, fiber_quad_w, spread_test_case, fiber_xdot, fiber_slip, fiber_coupling_force, &
       fiber_xc, fiber_uc, fiber_p, fiber_omega, fiber_force_total, fiber_torque_total

  implicit none

contains

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
      write(ifile,'(A)') 'itime time lag_total_Fx lag_total_Fy lag_total_Fz eul_total_Fx eul_total_Fy eul_total_Fz ' // &
           'abs_force_balance_Fx abs_force_balance_Fy abs_force_balance_Fz dt_stage lag_impulse_x lag_impulse_y lag_impulse_z ' // &
           'fluid_momentum_delta_x fluid_momentum_delta_y fluid_momentum_delta_z'
    endif

    write(ifile,'(I10,1X,18(ES24.16,1X))') itime, time, lag_total(1), lag_total(2), lag_total(3), eul_total(1), &
         eul_total(2), eul_total(3), abs_force_balance(1), abs_force_balance(2), abs_force_balance(3), dt_stage, &
         lag_impulse(1), lag_impulse(2), lag_impulse(3), euler_impulse(1), euler_impulse(2), euler_impulse(3)
    close(ifile)

  end subroutine write_fiber_rigid_two_way_momentum

end module fiber_io
