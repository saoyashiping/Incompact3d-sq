module fibre_prod_ibm_interpolation
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, &
                                      fibre_prod_grid_is_initialized, &
                                      fibre_prod_grid_compute_min_spacing, &
                                      fibre_prod_grid_compute_max_spacing
  use fibre_prod_ibm_delta, only : dp, fibre_prod_delta_weight_3d
  implicit none
  private

  public :: fibre_prod_ibm_interpolate_scalar
  public :: fibre_prod_ibm_interpolate_velocity

contains

  subroutine fibre_prod_ibm_interpolate_scalar(grid, scalar_field, point, value, status, weight_sum)
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(in) :: scalar_field(:, :, :)
    real(dp), intent(in) :: point(3)
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    real(dp), intent(out), optional :: weight_sum
    real(dp) :: wrapped_point(3)
    real(dp) :: local_weight_sum
    real(dp) :: weighted_value
    real(dp) :: hx
    real(dp) :: hy
    real(dp) :: hz
    real(dp) :: dx
    real(dp) :: dy
    real(dp) :: dz
    real(dp) :: weight
    integer :: i
    integer :: j
    integer :: k

    value = 0.0_dp
    local_weight_sum = 0.0_dp
    weighted_value = 0.0_dp
    status = 0

    if (.not. fibre_prod_grid_is_initialized(grid)) then
      status = 1
    else if (.not. all(ieee_is_finite(point))) then
      status = 2
    else if (.not. scalar_shape_matches_grid(grid, scalar_field)) then
      status = 3
    else
      wrapped_point = point
      call wrap_or_fail(grid%x, grid%periodic_x, wrapped_point(1), status)
      if (status == 0) call wrap_or_fail(grid%y, grid%periodic_y, wrapped_point(2), status)
      if (status == 0) call wrap_or_fail(grid%z, grid%periodic_z, wrapped_point(3), status)

      if (status == 0) then
        hx = max(fibre_prod_grid_compute_max_spacing(grid), fibre_prod_grid_compute_min_spacing(grid))
        hy = hx
        hz = hx
        do k = 1, grid%nz_local
          do j = 1, grid%ny_local
            do i = 1, grid%nx_local
              dx = periodic_distance(wrapped_point(1), grid%x(i), grid%x, grid%periodic_x)
              dy = periodic_distance(wrapped_point(2), grid%y(j), grid%y, grid%periodic_y)
              dz = periodic_distance(wrapped_point(3), grid%z(k), grid%z, grid%periodic_z)
              weight = fibre_prod_delta_weight_3d(dx, dy, dz, hx, hy, hz)
              local_weight_sum = local_weight_sum + weight
              weighted_value = weighted_value + weight * scalar_field(i, j, k)
            end do
          end do
        end do

        if (local_weight_sum > 0.0_dp .and. ieee_is_finite(local_weight_sum)) then
          value = weighted_value / local_weight_sum
        else
          status = 4
        end if
      end if
    end if

    if (present(weight_sum)) weight_sum = local_weight_sum
  end subroutine fibre_prod_ibm_interpolate_scalar

  subroutine fibre_prod_ibm_interpolate_velocity(grid, velocity_field, point, velocity, status, weight_sum)
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(in) :: velocity_field(:, :, :, :)
    real(dp), intent(in) :: point(3)
    real(dp), intent(out) :: velocity(3)
    integer, intent(out) :: status
    real(dp), intent(out), optional :: weight_sum
    real(dp) :: component_weight_sum
    integer :: component

    velocity = 0.0_dp
    status = 0
    component_weight_sum = 0.0_dp
    if (size(velocity_field, 4) /= 3) then
      status = 5
    else
      do component = 1, 3
        call fibre_prod_ibm_interpolate_scalar(grid, velocity_field(:, :, :, component), point, &
                                               velocity(component), status, component_weight_sum)
        if (status /= 0) exit
      end do
    end if

    if (present(weight_sum)) weight_sum = component_weight_sum
  end subroutine fibre_prod_ibm_interpolate_velocity

  pure logical function scalar_shape_matches_grid(grid, scalar_field) result(matches)
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(in) :: scalar_field(:, :, :)

    matches = size(scalar_field, 1) == grid%nx_local .and. &
              size(scalar_field, 2) == grid%ny_local .and. &
              size(scalar_field, 3) == grid%nz_local
  end function scalar_shape_matches_grid

  subroutine wrap_or_fail(coord, is_periodic, value, status)
    real(dp), intent(in) :: coord(:)
    logical, intent(in) :: is_periodic
    real(dp), intent(inout) :: value
    integer, intent(inout) :: status
    real(dp) :: period

    if (status /= 0) return
    if (value >= coord(1) .and. value <= coord(size(coord))) return

    if (.not. is_periodic) then
      status = 6
      return
    end if

    period = coord(size(coord)) - coord(1)
    if (period <= 0.0_dp) then
      status = 7
      return
    end if

    do while (value < coord(1))
      value = value + period
    end do
    do while (value > coord(size(coord)))
      value = value - period
    end do
  end subroutine wrap_or_fail

  pure real(dp) function periodic_distance(point_coord, grid_coord, coord, is_periodic) result(distance)
    real(dp), intent(in) :: point_coord
    real(dp), intent(in) :: grid_coord
    real(dp), intent(in) :: coord(:)
    logical, intent(in) :: is_periodic
    real(dp) :: period

    distance = point_coord - grid_coord
    if (is_periodic) then
      period = coord(size(coord)) - coord(1)
      if (period > 0.0_dp) distance = distance - nint(distance / period) * period
    end if
  end function periodic_distance

end module fibre_prod_ibm_interpolation
