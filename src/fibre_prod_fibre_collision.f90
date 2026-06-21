module fibre_prod_fibre_collision
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_collision_geometry, only : fibre_prod_distance_point_point, &
                                            fibre_prod_distance_segment_segment, &
                                            fibre_prod_collision_effective_gap, &
                                            fibre_prod_collision_geometry_validate
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_collision_pair_type
  public :: fibre_prod_collision_state_type
  public :: fibre_prod_collision_evaluate_two_fibres
  public :: fibre_prod_collision_force_candidate_two_fibres
  public :: fibre_prod_collision_status_string

  type :: fibre_prod_collision_pair_type
    integer :: fibre_i = 0
    integer :: fibre_j = 0
    integer :: node_i = 0
    integer :: node_j = 0
    integer :: segment_i = 0
    integer :: segment_j = 0
    real(dp) :: distance = huge(1.0_dp)
    real(dp) :: gap = huge(1.0_dp)
    real(dp) :: penetration = 0.0_dp
    real(dp) :: normal(3) = [1.0_dp, 0.0_dp, 0.0_dp]
    logical :: near_contact = .false.
    logical :: overlap_detected = .false.
  end type fibre_prod_collision_pair_type

  type :: fibre_prod_collision_state_type
    logical :: initialized = .false.
    logical :: near_contact_warning = .false.
    logical :: overlap_detected = .false.
    integer :: npairs_checked = 0
    integer :: npairs_near = 0
    integer :: npairs_overlap = 0
    real(dp) :: min_gap = huge(1.0_dp)
    real(dp) :: max_penetration = 0.0_dp
    real(dp) :: force_norm = 0.0_dp
    type(fibre_prod_collision_pair_type) :: closest_pair
  end type fibre_prod_collision_state_type

contains

  subroutine fibre_prod_collision_evaluate_two_fibres(x_i, x_j, radius_i, radius_j, warning_gap, &
                                                      state, status)
    real(dp), intent(in) :: x_i(:, :)
    real(dp), intent(in) :: x_j(:, :)
    real(dp), intent(in) :: radius_i
    real(dp), intent(in) :: radius_j
    real(dp), intent(in) :: warning_gap
    type(fibre_prod_collision_state_type), intent(out) :: state
    integer, intent(out) :: status
    type(fibre_prod_collision_pair_type) :: pair
    integer :: i
    integer :: j
    real(dp) :: distance
    real(dp) :: gap

    call reset_state(state)
    status = validate_inputs(x_i, x_j, radius_i, radius_j, warning_gap)
    if (status /= 0) return
    do i = 1, size(x_i, 1)
      do j = 1, size(x_j, 1)
        distance = fibre_prod_distance_point_point(x_i(i, :), x_j(j, :), status)
        if (status /= 0) return
        gap = fibre_prod_collision_effective_gap(distance, radius_i, radius_j, status)
        if (status /= 0) return
        call fill_pair(pair, i, j, 0, 0, distance, gap, x_i(i, :) - x_j(j, :), warning_gap)
        call accumulate_pair(state, pair)
      end do
    end do
    do i = 1, size(x_i, 1) - 1
      do j = 1, size(x_j, 1) - 1
        distance = fibre_prod_distance_segment_segment(x_i(i, :), x_i(i + 1, :), &
                                                       x_j(j, :), x_j(j + 1, :), status)
        if (status /= 0) return
        gap = fibre_prod_collision_effective_gap(distance, radius_i, radius_j, status)
        if (status /= 0) return
        call fill_pair(pair, i, j, i, j, distance, gap, segment_normal(x_i, x_j, i, j), warning_gap)
        call accumulate_pair(state, pair)
      end do
    end do
    state%initialized = .true.
  end subroutine fibre_prod_collision_evaluate_two_fibres

  subroutine fibre_prod_collision_force_candidate_two_fibres(x_i, v_i, x_j, v_j, radius_i, radius_j, &
                                                             warning_gap, k_collision, c_collision, &
                                                             force_i, force_j, state, status)
    real(dp), intent(in) :: x_i(:, :)
    real(dp), intent(in) :: v_i(:, :)
    real(dp), intent(in) :: x_j(:, :)
    real(dp), intent(in) :: v_j(:, :)
    real(dp), intent(in) :: radius_i
    real(dp), intent(in) :: radius_j
    real(dp), intent(in) :: warning_gap
    real(dp), intent(in) :: k_collision
    real(dp), intent(in) :: c_collision
    real(dp), intent(out) :: force_i(:, :)
    real(dp), intent(out) :: force_j(:, :)
    type(fibre_prod_collision_state_type), intent(out) :: state
    integer, intent(out) :: status
    real(dp) :: relative_velocity(3)
    real(dp) :: closing_speed
    real(dp) :: magnitude
    real(dp) :: force_vector(3)
    integer :: ni
    integer :: nj

    force_i = 0.0_dp
    force_j = 0.0_dp
    call fibre_prod_collision_evaluate_two_fibres(x_i, x_j, radius_i, radius_j, warning_gap, state, status)
    if (status /= 0) return
    if (.not. same_shape(x_i, v_i) .or. .not. same_shape(x_j, v_j) .or. &
        .not. same_shape(x_i, force_i) .or. .not. same_shape(x_j, force_j)) then
      status = 20
    else if (.not. all(ieee_is_finite(v_i)) .or. .not. all(ieee_is_finite(v_j))) then
      status = 21
    else if (.not. ieee_is_finite(k_collision) .or. k_collision < 0.0_dp) then
      status = 22
    else if (.not. ieee_is_finite(c_collision) .or. c_collision < 0.0_dp) then
      status = 23
    else if (state%overlap_detected) then
      ni = max(1, state%closest_pair%node_i)
      nj = max(1, state%closest_pair%node_j)
      relative_velocity = v_i(ni, :) - v_j(nj, :)
      closing_speed = max(0.0_dp, -sum(relative_velocity * state%closest_pair%normal))
      magnitude = k_collision * state%closest_pair%penetration + c_collision * closing_speed
      force_vector = magnitude * state%closest_pair%normal
      force_i(ni, :) = force_i(ni, :) + force_vector
      force_j(nj, :) = force_j(nj, :) - force_vector
    end if
    if (status == 0) then
      if (all(ieee_is_finite(force_i)) .and. all(ieee_is_finite(force_j))) then
        state%force_norm = sqrt(sum(force_i * force_i) + sum(force_j * force_j))
      else
        status = 24
      end if
    end if
  end subroutine fibre_prod_collision_force_candidate_two_fibres

  pure function fibre_prod_collision_status_string(status) result(message)
    integer, intent(in) :: status
    character(len=96) :: message

    select case (status)
    case (0)
      message = 'success'
    case (10)
      message = 'invalid_collision_geometry'
    case (11)
      message = 'invalid_state_shape'
    case (12)
      message = 'non_finite_coordinate'
    case (20)
      message = 'force_or_velocity_shape_mismatch'
    case (21)
      message = 'non_finite_velocity'
    case (22)
      message = 'invalid_collision_stiffness'
    case (23)
      message = 'invalid_collision_damping'
    case default
      message = 'fibre_collision_failure'
    end select
  end function fibre_prod_collision_status_string

  pure integer function validate_inputs(x_i, x_j, radius_i, radius_j, warning_gap) result(status)
    real(dp), intent(in) :: x_i(:, :)
    real(dp), intent(in) :: x_j(:, :)
    real(dp), intent(in) :: radius_i
    real(dp), intent(in) :: radius_j
    real(dp), intent(in) :: warning_gap

    status = fibre_prod_collision_geometry_validate(radius_i, radius_j, warning_gap)
    if (status /= 0) then
      status = 10
    else if (size(x_i, 1) < 2 .or. size(x_j, 1) < 2 .or. size(x_i, 2) /= 3 .or. size(x_j, 2) /= 3) then
      status = 11
    else if (.not. all(ieee_is_finite(x_i)) .or. .not. all(ieee_is_finite(x_j))) then
      status = 12
    end if
  end function validate_inputs

  pure logical function same_shape(a, b) result(matches)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in) :: b(:, :)

    matches = size(a, 1) == size(b, 1) .and. size(a, 2) == size(b, 2) .and. size(a, 2) == 3
  end function same_shape

  pure subroutine reset_state(state)
    type(fibre_prod_collision_state_type), intent(out) :: state

    state%initialized = .false.
    state%near_contact_warning = .false.
    state%overlap_detected = .false.
    state%npairs_checked = 0
    state%npairs_near = 0
    state%npairs_overlap = 0
    state%min_gap = huge(1.0_dp)
    state%max_penetration = 0.0_dp
    state%force_norm = 0.0_dp
    state%closest_pair = fibre_prod_collision_pair_type()
  end subroutine reset_state

  pure subroutine fill_pair(pair, node_i, node_j, segment_i, segment_j, distance, gap, raw_normal, warning_gap)
    type(fibre_prod_collision_pair_type), intent(out) :: pair
    integer, intent(in) :: node_i
    integer, intent(in) :: node_j
    integer, intent(in) :: segment_i
    integer, intent(in) :: segment_j
    real(dp), intent(in) :: distance
    real(dp), intent(in) :: gap
    real(dp), intent(in) :: raw_normal(3)
    real(dp), intent(in) :: warning_gap
    real(dp) :: norm_value

    pair%fibre_i = 1
    pair%fibre_j = 2
    pair%node_i = node_i
    pair%node_j = node_j
    pair%segment_i = segment_i
    pair%segment_j = segment_j
    pair%distance = distance
    pair%gap = gap
    pair%penetration = max(0.0_dp, -gap)
    norm_value = sqrt(sum(raw_normal * raw_normal))
    if (norm_value > epsilon(1.0_dp)) then
      pair%normal = raw_normal / norm_value
    else
      pair%normal = [1.0_dp, 0.0_dp, 0.0_dp]
    end if
    pair%near_contact = gap >= 0.0_dp .and. gap <= warning_gap
    pair%overlap_detected = gap < 0.0_dp
  end subroutine fill_pair

  pure subroutine accumulate_pair(state, pair)
    type(fibre_prod_collision_state_type), intent(inout) :: state
    type(fibre_prod_collision_pair_type), intent(in) :: pair

    state%npairs_checked = state%npairs_checked + 1
    if (pair%near_contact) then
      state%near_contact_warning = .true.
      state%npairs_near = state%npairs_near + 1
    end if
    if (pair%overlap_detected) then
      state%overlap_detected = .true.
      state%npairs_overlap = state%npairs_overlap + 1
      state%max_penetration = max(state%max_penetration, pair%penetration)
    end if
    if (pair%gap < state%min_gap) then
      state%min_gap = pair%gap
      state%closest_pair = pair
    end if
  end subroutine accumulate_pair

  pure function segment_normal(x_i, x_j, segment_i, segment_j) result(normal)
    real(dp), intent(in) :: x_i(:, :)
    real(dp), intent(in) :: x_j(:, :)
    integer, intent(in) :: segment_i
    integer, intent(in) :: segment_j
    real(dp) :: normal(3)
    real(dp) :: midpoint_i(3)
    real(dp) :: midpoint_j(3)

    midpoint_i = 0.5_dp * (x_i(segment_i, :) + x_i(segment_i + 1, :))
    midpoint_j = 0.5_dp * (x_j(segment_j, :) + x_j(segment_j + 1, :))
    normal = midpoint_i - midpoint_j
  end function segment_normal

end module fibre_prod_fibre_collision
