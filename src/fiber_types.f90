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
  real(mytype), allocatable, dimension(:,:) :: fiber_uinterp
  real(mytype), allocatable, dimension(:,:) :: fiber_uexact
  real(mytype), allocatable, dimension(:,:) :: fiber_uerror
  real(mytype), allocatable, dimension(:) :: fiber_sumw
  real(mytype), allocatable, dimension(:) :: fiber_quad_w
  real(mytype), allocatable, dimension(:,:) :: fiber_test_force
  real(mytype) :: fiber_interp_max_error
  logical :: interp_test_active
  integer :: interp_test_case
  logical :: interp_solver_test_active
  integer :: interp_solver_output_step
  logical :: spread_test_active
  integer :: spread_test_case

contains

  subroutine fiber_set_defaults()

    fiber_active = .false.
    fiber_nl = 2
    fiber_length = 1._mytype
    fiber_center = 0._mytype
    fiber_direction = (/1._mytype, 0._mytype, 0._mytype/)

    interp_test_active = .false.
    interp_test_case = 1
    fiber_interp_max_error = 0._mytype
    interp_solver_test_active = .false.
    interp_solver_output_step = 1
    spread_test_active = .false.
    spread_test_case = 1

    if (allocated(fiber_x)) deallocate(fiber_x)
    if (allocated(fiber_uinterp)) deallocate(fiber_uinterp)
    if (allocated(fiber_uexact)) deallocate(fiber_uexact)
    if (allocated(fiber_uerror)) deallocate(fiber_uerror)
    if (allocated(fiber_sumw)) deallocate(fiber_sumw)
    if (allocated(fiber_quad_w)) deallocate(fiber_quad_w)
    if (allocated(fiber_test_force)) deallocate(fiber_test_force)

  end subroutine fiber_set_defaults

end module fiber_types
