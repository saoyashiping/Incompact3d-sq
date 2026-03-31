!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_io

  use decomp_2d_mpi, only : nrank
  use fiber_types, only : fiber_active, fiber_nl, fiber_x, fiber_uinterp, fiber_uexact, fiber_uerror, &
       fiber_sumw, fiber_test_force, fiber_quad_w, spread_test_case

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
    write(ifile,'(A)') 'test_case lag_Fx lag_Fy lag_Fz eul_Fx eul_Fy eul_Fz abs_err_Fx abs_err_Fy abs_err_Fz rel_err_Fx rel_err_Fy rel_err_Fz spread_sumw_min spread_sumw_max'
    write(ifile,'(I8,1X,14(ES24.16,1X))') spread_test_case, lag_total(1), lag_total(2), lag_total(3), &
         eul_total(1), eul_total(2), eul_total(3), abs_error(1), abs_error(2), abs_error(3), &
         rel_error(1), rel_error(2), rel_error(3), sumw_min, sumw_max
    close(ifile)

  end subroutine write_fiber_spread_summary

end module fiber_io
