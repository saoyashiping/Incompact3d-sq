module fibre_prod_wall_contact
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_wall_geometry, only : fibre_prod_wall_geometry_type, fibre_prod_wall_geometry_validate, &
                                       fibre_prod_wall_gap_point
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_wall_contact_state_type
  public :: fibre_prod_wall_contact_evaluate
  public :: fibre_prod_wall_contact_force_candidate
  public :: fibre_prod_wall_contact_status_string

  type :: fibre_prod_wall_contact_state_type
    logical :: contact_active = .false.
    logical :: near_wall_warning = .false.
    logical :: penetration_detected = .false.
    integer :: nearest_wall = 0
    real(dp) :: min_gap = huge(1.0_dp)
    real(dp) :: max_penetration = 0.0_dp
    real(dp) :: normal_force_norm = 0.0_dp
  end type fibre_prod_wall_contact_state_type

contains

  subroutine fibre_prod_wall_contact_evaluate(geometry, x, state, status)
    type(fibre_prod_wall_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: x(:, :)
    type(fibre_prod_wall_contact_state_type), intent(out) :: state
    integer, intent(out) :: status
    integer :: inode
    real(dp) :: lower_gap
    real(dp) :: upper_gap
    real(dp) :: min_gap

    call reset_contact_state(state)
    status = validate_contact_inputs(geometry, x)
    if (status /= 0) return
    do inode = 1, size(x, 1)
      call fibre_prod_wall_gap_point(geometry, x(inode, 2), lower_gap, upper_gap, min_gap, status)
      if (status /= 0) return
      if (min_gap < state%min_gap) then
        state%min_gap = min_gap
        if (lower_gap <= upper_gap) then
          state%nearest_wall = -1
        else
          state%nearest_wall = 1
        end if
      end if
      if (min_gap < 0.0_dp) then
        state%penetration_detected = .true.
        state%contact_active = .true.
        state%max_penetration = max(state%max_penetration, -min_gap)
      else if (min_gap <= geometry%warning_gap) then
        state%near_wall_warning = .true.
      end if
    end do
  end subroutine fibre_prod_wall_contact_evaluate

  subroutine fibre_prod_wall_contact_force_candidate(geometry, x, v, k_wall, c_wall, force, state, status)
    type(fibre_prod_wall_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in) :: v(:, :)
    real(dp), intent(in) :: k_wall
    real(dp), intent(in) :: c_wall
    real(dp), intent(out) :: force(:, :)
    type(fibre_prod_wall_contact_state_type), intent(out) :: state
    integer, intent(out) :: status
    integer :: inode
    real(dp) :: lower_gap
    real(dp) :: upper_gap
    real(dp) :: min_gap
    real(dp) :: penetration
    real(dp) :: damping

    force = 0.0_dp
    call fibre_prod_wall_contact_evaluate(geometry, x, state, status)
    if (status /= 0) return
    if (size(v, 1) /= size(x, 1) .or. size(v, 2) /= 3 .or. size(force, 1) /= size(x, 1) .or. &
        size(force, 2) /= 3) then
      status = 20
    else if (.not. all(ieee_is_finite(v))) then
      status = 21
    else if (.not. ieee_is_finite(k_wall) .or. k_wall < 0.0_dp) then
      status = 22
    else if (.not. ieee_is_finite(c_wall) .or. c_wall < 0.0_dp) then
      status = 23
    else
      do inode = 1, size(x, 1)
        call fibre_prod_wall_gap_point(geometry, x(inode, 2), lower_gap, upper_gap, min_gap, status)
        if (status /= 0) return
        if (lower_gap < 0.0_dp .and. lower_gap <= upper_gap) then
          penetration = -lower_gap
          damping = c_wall * max(0.0_dp, -v(inode, 2))
          force(inode, 2) = force(inode, 2) + k_wall * penetration + damping
        else if (upper_gap < 0.0_dp) then
          penetration = -upper_gap
          damping = c_wall * max(0.0_dp, v(inode, 2))
          force(inode, 2) = force(inode, 2) - k_wall * penetration - damping
        end if
      end do
      if (all(ieee_is_finite(force))) then
        state%normal_force_norm = sqrt(sum(force(:, 2) * force(:, 2)))
      else
        status = 24
      end if
    end if
  end subroutine fibre_prod_wall_contact_force_candidate

  pure function fibre_prod_wall_contact_status_string(status) result(message)
    integer, intent(in) :: status
    character(len=96) :: message

    select case (status)
    case (0)
      message = 'success'
    case (10)
      message = 'invalid_wall_point'
    case (20)
      message = 'shape_mismatch'
    case (21)
      message = 'non_finite_velocity'
    case (22)
      message = 'invalid_k_wall'
    case (23)
      message = 'invalid_c_wall'
    case default
      message = 'wall_contact_failure'
    end select
  end function fibre_prod_wall_contact_status_string

  pure integer function validate_contact_inputs(geometry, x) result(status)
    type(fibre_prod_wall_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: x(:, :)

    status = fibre_prod_wall_geometry_validate(geometry)
    if (status == 0 .and. (size(x, 1) < 1 .or. size(x, 2) /= 3)) status = 10
    if (status == 0 .and. .not. all(ieee_is_finite(x))) status = 11
  end function validate_contact_inputs

  pure subroutine reset_contact_state(state)
    type(fibre_prod_wall_contact_state_type), intent(out) :: state

    state%contact_active = .false.
    state%near_wall_warning = .false.
    state%penetration_detected = .false.
    state%nearest_wall = 0
    state%min_gap = huge(1.0_dp)
    state%max_penetration = 0.0_dp
    state%normal_force_norm = 0.0_dp
  end subroutine reset_contact_state

end module fibre_prod_wall_contact
