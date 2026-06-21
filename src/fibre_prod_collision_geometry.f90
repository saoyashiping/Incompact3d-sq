module fibre_prod_collision_geometry
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_distance_point_point
  public :: fibre_prod_distance_point_segment
  public :: fibre_prod_distance_segment_segment
  public :: fibre_prod_collision_effective_gap
  public :: fibre_prod_collision_geometry_validate

contains

  real(dp) function fibre_prod_distance_point_point(p, q, status) result(distance)
    real(dp), intent(in) :: p(3)
    real(dp), intent(in) :: q(3)
    integer, intent(out) :: status
    real(dp) :: d(3)

    status = 0
    distance = huge(1.0_dp)
    if (.not. all(ieee_is_finite(p)) .or. .not. all(ieee_is_finite(q))) then
      status = 1
    else
      d = p - q
      distance = sqrt(sum(d * d))
    end if
  end function fibre_prod_distance_point_point

  real(dp) function fibre_prod_distance_point_segment(p, a, b, status) result(distance)
    real(dp), intent(in) :: p(3)
    real(dp), intent(in) :: a(3)
    real(dp), intent(in) :: b(3)
    integer, intent(out) :: status
    real(dp) :: ab(3)
    real(dp) :: closest(3)
    real(dp) :: denom
    real(dp) :: t

    status = 0
    distance = huge(1.0_dp)
    if (.not. all(ieee_is_finite(p)) .or. .not. all(ieee_is_finite(a)) .or. &
        .not. all(ieee_is_finite(b))) then
      status = 1
      return
    end if
    ab = b - a
    denom = sum(ab * ab)
    if (denom <= epsilon(1.0_dp)) then
      distance = fibre_prod_distance_point_point(p, a, status)
    else
      t = max(0.0_dp, min(1.0_dp, sum((p - a) * ab) / denom))
      closest = a + t * ab
      distance = fibre_prod_distance_point_point(p, closest, status)
    end if
  end function fibre_prod_distance_point_segment

  real(dp) function fibre_prod_distance_segment_segment(a0, a1, b0, b1, status) result(distance)
    real(dp), intent(in) :: a0(3)
    real(dp), intent(in) :: a1(3)
    real(dp), intent(in) :: b0(3)
    real(dp), intent(in) :: b1(3)
    integer, intent(out) :: status
    real(dp) :: u(3)
    real(dp) :: v(3)
    real(dp) :: w(3)
    real(dp) :: a
    real(dp) :: b
    real(dp) :: c
    real(dp) :: d
    real(dp) :: e
    real(dp) :: denom
    real(dp) :: s
    real(dp) :: t
    real(dp) :: p_closest(3)
    real(dp) :: q_closest(3)

    status = 0
    distance = huge(1.0_dp)
    if (.not. all(ieee_is_finite(a0)) .or. .not. all(ieee_is_finite(a1)) .or. &
        .not. all(ieee_is_finite(b0)) .or. .not. all(ieee_is_finite(b1))) then
      status = 1
      return
    end if
    u = a1 - a0
    v = b1 - b0
    w = a0 - b0
    a = sum(u * u)
    b = sum(u * v)
    c = sum(v * v)
    d = sum(u * w)
    e = sum(v * w)
    if (a <= epsilon(1.0_dp) .and. c <= epsilon(1.0_dp)) then
      distance = fibre_prod_distance_point_point(a0, b0, status)
      return
    else if (a <= epsilon(1.0_dp)) then
      distance = fibre_prod_distance_point_segment(a0, b0, b1, status)
      return
    else if (c <= epsilon(1.0_dp)) then
      distance = fibre_prod_distance_point_segment(b0, a0, a1, status)
      return
    end if
    denom = a * c - b * b
    if (abs(denom) <= epsilon(1.0_dp)) then
      s = 0.0_dp
    else
      s = max(0.0_dp, min(1.0_dp, (b * e - c * d) / denom))
    end if
    t = max(0.0_dp, min(1.0_dp, (b * s + e) / c))
    s = max(0.0_dp, min(1.0_dp, (b * t - d) / a))
    p_closest = a0 + s * u
    q_closest = b0 + t * v
    distance = fibre_prod_distance_point_point(p_closest, q_closest, status)
  end function fibre_prod_distance_segment_segment

  real(dp) function fibre_prod_collision_effective_gap(distance, radius_i, radius_j, status) result(gap)
    real(dp), intent(in) :: distance
    real(dp), intent(in) :: radius_i
    real(dp), intent(in) :: radius_j
    integer, intent(out) :: status

    status = fibre_prod_collision_geometry_validate(radius_i, radius_j, 0.0_dp)
    if (status == 0 .and. (.not. ieee_is_finite(distance) .or. distance < 0.0_dp)) status = 4
    if (status == 0) then
      gap = distance - radius_i - radius_j
    else
      gap = huge(1.0_dp)
    end if
  end function fibre_prod_collision_effective_gap

  pure integer function fibre_prod_collision_geometry_validate(radius_i, radius_j, warning_gap) result(status)
    real(dp), intent(in) :: radius_i
    real(dp), intent(in) :: radius_j
    real(dp), intent(in) :: warning_gap

    status = 0
    if (.not. ieee_is_finite(radius_i) .or. radius_i < 0.0_dp) status = 1
    if (status == 0 .and. (.not. ieee_is_finite(radius_j) .or. radius_j < 0.0_dp)) status = 2
    if (status == 0 .and. (.not. ieee_is_finite(warning_gap) .or. warning_gap < 0.0_dp)) status = 3
  end function fibre_prod_collision_geometry_validate

end module fibre_prod_collision_geometry
