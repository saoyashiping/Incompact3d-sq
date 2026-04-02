!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_init

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use fiber_types, only : fiber_active, fiber_nl, fiber_length, fiber_center, fiber_direction, fiber_x, fiber_quad_w

  implicit none

contains

  subroutine init_fiber()

    integer :: l
    real(mytype) :: dir_norm, ds, s

    if (.not.fiber_active) return

    if (fiber_nl < 2) then
      if (nrank == 0) write(*,*) 'Error: fiber_nl must be >= 2 when fiber_active is true.'
      stop
    endif

    if (fiber_length <= 0._mytype) then
      if (nrank == 0) write(*,*) 'Error: fiber_length must be > 0 when fiber_active is true.'
      stop
    endif

    dir_norm = sqrt(sum(fiber_direction**2))
    if (dir_norm <= 0._mytype) then
      if (nrank == 0) write(*,*) 'Error: fiber_direction must be a non-zero vector when fiber_active is true.'
      stop
    endif

    fiber_direction = fiber_direction / dir_norm

    if (allocated(fiber_x)) deallocate(fiber_x)
    allocate(fiber_x(3, fiber_nl))

    ds = fiber_length / real(fiber_nl - 1, mytype)

    if (allocated(fiber_quad_w)) deallocate(fiber_quad_w)
    allocate(fiber_quad_w(fiber_nl))
    fiber_quad_w = ds
    fiber_quad_w(1) = 0.5_mytype * ds
    fiber_quad_w(fiber_nl) = 0.5_mytype * ds

    do l = 1, fiber_nl
      s = -0.5_mytype * fiber_length + real(l - 1, mytype) * ds
      fiber_x(:, l) = fiber_center + s * fiber_direction
    enddo

  end subroutine init_fiber

end module fiber_init
