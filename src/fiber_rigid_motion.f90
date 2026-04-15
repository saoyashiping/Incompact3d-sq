!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fiber_rigid_motion

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : nrank
  use param, only : xlx, zlz
  use fiber_types, only : fiber_active, fiber_nl, fiber_x, fiber_x_ref, fiber_xdot, &
       rigid_motion_case, rigid_translation_velocity

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

  subroutine init_rigid_motion_reference()

    if (.not.fiber_active) return

    if (.not.allocated(fiber_x_ref)) allocate(fiber_x_ref(3, fiber_nl))
    if (.not.allocated(fiber_xdot)) allocate(fiber_xdot(3, fiber_nl))

    fiber_x_ref = fiber_x
    fiber_xdot = 0._mytype

  end subroutine init_rigid_motion_reference

  subroutine update_rigid_motion(time)

    real(mytype), intent(in) :: time

    integer :: l

    if (.not.fiber_active) return

    if (.not.allocated(fiber_x_ref)) then
      if (nrank == 0) write(*,*) 'Error: fiber_x_ref is not allocated in update_rigid_motion.'
      stop
    endif

    if (.not.allocated(fiber_xdot)) allocate(fiber_xdot(3, fiber_nl))

    select case (rigid_motion_case)
    case (1)
      fiber_x = fiber_x_ref
      fiber_xdot = 0._mytype
    case (2)
      do l = 1, fiber_nl
        fiber_x(1,l) = wrap_periodic(fiber_x_ref(1,l) + rigid_translation_velocity(1) * time, xlx)
        fiber_x(2,l) = fiber_x_ref(2,l) + rigid_translation_velocity(2) * time
        fiber_x(3,l) = wrap_periodic(fiber_x_ref(3,l) + rigid_translation_velocity(3) * time, zlz)
      enddo
      do l = 1, fiber_nl
        fiber_xdot(:,l) = rigid_translation_velocity
      enddo
    case default
      if (nrank == 0) write(*,*) 'Error: rigid_motion_case must be 1 (stationary) or 2 (constant translation).'
      stop
    end select

  end subroutine update_rigid_motion

  subroutine compute_rigid_spacing_error(spacing_error_max)

    real(mytype), intent(out) :: spacing_error_max

    integer :: l
    real(mytype) :: ds_ref, ds_now

    spacing_error_max = 0._mytype

    if (.not.fiber_active) return
    if (fiber_nl < 2) return

    ds_ref = sqrt(sum((fiber_x_ref(:,2) - fiber_x_ref(:,1))**2))

    do l = 1, fiber_nl - 1
      ds_now = sqrt(sum((fiber_x(:,l+1) - fiber_x(:,l))**2))
      spacing_error_max = max(spacing_error_max, abs(ds_now - ds_ref))
    enddo

  end subroutine compute_rigid_spacing_error

end module fiber_rigid_motion
