!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_io

  use decomp_2d_mpi, only : nrank
  use fiber_types, only : fiber_active, fiber_nl, fiber_x

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

end module fiber_io
