module fibre_prod_wall_geometry
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_wall_geometry_type
  public :: fibre_prod_wall_geometry_init
  public :: fibre_prod_wall_geometry_validate
  public :: fibre_prod_wall_gap_point
  public :: fibre_prod_wall_gap_segment_sampled
  public :: fibre_prod_wall_destroy

  type :: fibre_prod_wall_geometry_type
    real(dp) :: y_lower = 0.0_dp
    real(dp) :: y_upper = 0.0_dp
    real(dp) :: fibre_radius = 0.0_dp
    real(dp) :: warning_gap = 0.0_dp
    logical :: initialized = .false.
  end type fibre_prod_wall_geometry_type

contains

  subroutine fibre_prod_wall_geometry_init(geometry, y_lower, y_upper, fibre_radius, warning_gap, status)
    type(fibre_prod_wall_geometry_type), intent(inout) :: geometry
    real(dp), intent(in) :: y_lower
    real(dp), intent(in) :: y_upper
    real(dp), intent(in) :: fibre_radius
    real(dp), intent(in) :: warning_gap
    integer, intent(out) :: status

    geometry%y_lower = y_lower
    geometry%y_upper = y_upper
    geometry%fibre_radius = fibre_radius
    geometry%warning_gap = warning_gap
    geometry%initialized = .true.
    status = fibre_prod_wall_geometry_validate(geometry)
    if (status /= 0) call fibre_prod_wall_destroy(geometry)
  end subroutine fibre_prod_wall_geometry_init

  pure integer function fibre_prod_wall_geometry_validate(geometry) result(status)
    type(fibre_prod_wall_geometry_type), intent(in) :: geometry

    status = 0
    if (.not. geometry%initialized) status = 1
    if (status == 0 .and. (.not. ieee_is_finite(geometry%y_lower))) status = 2
    if (status == 0 .and. (.not. ieee_is_finite(geometry%y_upper))) status = 3
    if (status == 0 .and. geometry%y_upper <= geometry%y_lower) status = 4
    if (status == 0 .and. (.not. ieee_is_finite(geometry%fibre_radius) .or. &
        geometry%fibre_radius < 0.0_dp)) status = 5
    if (status == 0 .and. (.not. ieee_is_finite(geometry%warning_gap) .or. &
        geometry%warning_gap < 0.0_dp)) status = 6
  end function fibre_prod_wall_geometry_validate

  pure subroutine fibre_prod_wall_gap_point(geometry, y_node, lower_gap, upper_gap, min_gap, status)
    type(fibre_prod_wall_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: y_node
    real(dp), intent(out) :: lower_gap
    real(dp), intent(out) :: upper_gap
    real(dp), intent(out) :: min_gap
    integer, intent(out) :: status

    lower_gap = huge(1.0_dp)
    upper_gap = huge(1.0_dp)
    min_gap = huge(1.0_dp)
    status = fibre_prod_wall_geometry_validate(geometry)
    if (status /= 0) return
    if (.not. ieee_is_finite(y_node)) then
      status = 10
      return
    end if
    lower_gap = y_node - geometry%y_lower - geometry%fibre_radius
    upper_gap = geometry%y_upper - y_node - geometry%fibre_radius
    min_gap = min(lower_gap, upper_gap)
  end subroutine fibre_prod_wall_gap_point

  pure subroutine fibre_prod_wall_gap_segment_sampled(geometry, y0, y1, nsample, min_gap, status)
    type(fibre_prod_wall_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: y0
    real(dp), intent(in) :: y1
    integer, intent(in) :: nsample
    real(dp), intent(out) :: min_gap
    integer, intent(out) :: status
    integer :: isample
    real(dp) :: alpha
    real(dp) :: y_node
    real(dp) :: lower_gap
    real(dp) :: upper_gap
    real(dp) :: point_gap

    min_gap = huge(1.0_dp)
    status = fibre_prod_wall_geometry_validate(geometry)
    if (status /= 0) return
    if (nsample < 2 .or. .not. ieee_is_finite(y0) .or. .not. ieee_is_finite(y1)) then
      status = 11
      return
    end if
    do isample = 1, nsample
      alpha = real(isample - 1, dp) / real(nsample - 1, dp)
      y_node = (1.0_dp - alpha) * y0 + alpha * y1
      call fibre_prod_wall_gap_point(geometry, y_node, lower_gap, upper_gap, point_gap, status)
      if (status /= 0) return
      min_gap = min(min_gap, point_gap)
    end do
  end subroutine fibre_prod_wall_gap_segment_sampled

  subroutine fibre_prod_wall_destroy(geometry)
    type(fibre_prod_wall_geometry_type), intent(inout) :: geometry

    geometry%y_lower = 0.0_dp
    geometry%y_upper = 0.0_dp
    geometry%fibre_radius = 0.0_dp
    geometry%warning_gap = 0.0_dp
    geometry%initialized = .false.
  end subroutine fibre_prod_wall_destroy

end module fibre_prod_wall_geometry
