!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_types

  use decomp_2d_constants, only : mytype

  implicit none

  logical :: fiber_active
  integer :: fiber_nl
  real(mytype) :: fiber_length
  real(mytype), dimension(3) :: fiber_center
  real(mytype), dimension(3) :: fiber_direction
  real(mytype), allocatable, dimension(:,:) :: fiber_x

contains

  subroutine fiber_set_defaults()

    fiber_active = .false.
    fiber_nl = 2
    fiber_length = 1._mytype
    fiber_center = 0._mytype
    fiber_direction = (/1._mytype, 0._mytype, 0._mytype/)

    if (allocated(fiber_x)) deallocate(fiber_x)

  end subroutine fiber_set_defaults

end module fiber_types
