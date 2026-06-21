module fibre_prod_ibm_spreading
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_is_initialized, &
                                      fibre_prod_grid_compute_min_spacing, &
                                      fibre_prod_grid_compute_max_spacing
  use fibre_prod_ibm_delta, only : dp, fibre_prod_delta_weight_3d
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, &
                                          fibre_prod_force_buffer_is_finite
  implicit none
  private

  public :: fibre_prod_spread_point_force
  public :: fibre_prod_spread_multiple_point_forces
  public :: fibre_prod_spreading_weight_sum
  public :: fibre_prod_spreading_status_string

contains

  subroutine fibre_prod_spread_point_force(grid, buffer, point, force, status, weight_sum)
    type(fibre_prod_grid_type), intent(in) :: grid
    type(fibre_prod_force_buffer_type), intent(inout) :: buffer
    real(dp), intent(in) :: point(3)
    real(dp), intent(in) :: force(3)
    integer, intent(out) :: status
    real(dp), intent(out), optional :: weight_sum
    real(dp) :: wrapped_point(3)
    real(dp) :: local_weight_sum
    real(dp) :: hx
    real(dp) :: hy
    real(dp) :: hz
    real(dp) :: dx
    real(dp) :: dy
    real(dp) :: dz
    real(dp) :: weight
    real(dp) :: density_scale
    integer :: i
    integer :: j
    integer :: k

    status = 0
    local_weight_sum = 0.0_dp
    if (.not. fibre_prod_grid_is_initialized(grid)) then
      status = 1
    else if (.not. fibre_prod_force_buffer_is_finite(buffer)) then
      status = 2
    else if (.not. buffer_shape_matches_grid(buffer, grid)) then
      status = 3
    else if (.not. all(ieee_is_finite(point)) .or. .not. all(ieee_is_finite(force))) then
      status = 4
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
              weight = fibre_prod_delta_weight_3d(dx, dy, dz, hx, hy, hz) * grid%cell_volume(i, j, k)
              local_weight_sum = local_weight_sum + weight
            end do
          end do
        end do
        if (local_weight_sum <= 0.0_dp .or. .not. ieee_is_finite(local_weight_sum)) then
          status = 5
        else
          do k = 1, grid%nz_local
            do j = 1, grid%ny_local
              do i = 1, grid%nx_local
                dx = periodic_distance(wrapped_point(1), grid%x(i), grid%x, grid%periodic_x)
                dy = periodic_distance(wrapped_point(2), grid%y(j), grid%y, grid%periodic_y)
                dz = periodic_distance(wrapped_point(3), grid%z(k), grid%z, grid%periodic_z)
                weight = fibre_prod_delta_weight_3d(dx, dy, dz, hx, hy, hz) * grid%cell_volume(i, j, k)
                density_scale = weight / local_weight_sum / grid%cell_volume(i, j, k)
                buffer%fx(i, j, k) = buffer%fx(i, j, k) + force(1) * density_scale
                buffer%fy(i, j, k) = buffer%fy(i, j, k) + force(2) * density_scale
                buffer%fz(i, j, k) = buffer%fz(i, j, k) + force(3) * density_scale
              end do
            end do
          end do
        end if
      end if
    end if
    if (present(weight_sum)) weight_sum = local_weight_sum
  end subroutine fibre_prod_spread_point_force

  subroutine fibre_prod_spread_multiple_point_forces(grid, buffer, points, forces, status, min_weight_sum)
    type(fibre_prod_grid_type), intent(in) :: grid
    type(fibre_prod_force_buffer_type), intent(inout) :: buffer
    real(dp), intent(in) :: points(:, :)
    real(dp), intent(in) :: forces(:, :)
    integer, intent(out) :: status
    real(dp), intent(out), optional :: min_weight_sum
    real(dp) :: current_weight_sum
    real(dp) :: local_min_weight_sum
    integer :: lpoint

    status = 0
    local_min_weight_sum = huge(1.0_dp)
    if (size(points, 2) /= 3 .or. size(forces, 2) /= 3 .or. size(points, 1) /= size(forces, 1)) then
      status = 10
    else
      do lpoint = 1, size(points, 1)
        call fibre_prod_spread_point_force(grid, buffer, points(lpoint, :), forces(lpoint, :), &
                                           status, current_weight_sum)
        if (status /= 0) exit
        local_min_weight_sum = min(local_min_weight_sum, current_weight_sum)
      end do
    end if
    if (present(min_weight_sum)) min_weight_sum = local_min_weight_sum
  end subroutine fibre_prod_spread_multiple_point_forces

  subroutine fibre_prod_spreading_weight_sum(grid, point, weight_sum, status)
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(in) :: point(3)
    real(dp), intent(out) :: weight_sum
    integer, intent(out) :: status
    weight_sum = 0.0_dp
    status = 99
    call compute_weight_sum_only(grid, point, weight_sum, status)
  end subroutine fibre_prod_spreading_weight_sum

  pure function fibre_prod_spreading_status_string(status) result(message)
    integer, intent(in) :: status
    character(len=96) :: message

    select case (status)
    case (0)
      message = 'success'
    case (1)
      message = 'grid_not_initialized'
    case (2)
      message = 'force_buffer_not_finite'
    case (3)
      message = 'force_buffer_shape_mismatch'
    case (4)
      message = 'non_finite_input'
    case (5)
      message = 'non_positive_weight_sum'
    case (6)
      message = 'non_periodic_out_of_domain'
    case default
      message = 'unknown_spreading_status'
    end select
  end function fibre_prod_spreading_status_string

  subroutine compute_weight_sum_only(grid, point, weight_sum, status)
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(in) :: point(3)
    real(dp), intent(out) :: weight_sum
    integer, intent(out) :: status
    real(dp) :: wrapped_point(3)
    real(dp) :: hx
    real(dp) :: dx
    real(dp) :: dy
    real(dp) :: dz
    integer :: i
    integer :: j
    integer :: k

    weight_sum = 0.0_dp
    status = 0
    if (.not. fibre_prod_grid_is_initialized(grid)) then
      status = 1
    else if (.not. all(ieee_is_finite(point))) then
      status = 4
    else
      wrapped_point = point
      call wrap_or_fail(grid%x, grid%periodic_x, wrapped_point(1), status)
      if (status == 0) call wrap_or_fail(grid%y, grid%periodic_y, wrapped_point(2), status)
      if (status == 0) call wrap_or_fail(grid%z, grid%periodic_z, wrapped_point(3), status)
      if (status == 0) then
        hx = max(fibre_prod_grid_compute_max_spacing(grid), fibre_prod_grid_compute_min_spacing(grid))
        do k = 1, grid%nz_local
          do j = 1, grid%ny_local
            do i = 1, grid%nx_local
              dx = periodic_distance(wrapped_point(1), grid%x(i), grid%x, grid%periodic_x)
              dy = periodic_distance(wrapped_point(2), grid%y(j), grid%y, grid%periodic_y)
              dz = periodic_distance(wrapped_point(3), grid%z(k), grid%z, grid%periodic_z)
              weight_sum = weight_sum + fibre_prod_delta_weight_3d(dx, dy, dz, hx, hx, hx) * &
                                         grid%cell_volume(i, j, k)
            end do
          end do
        end do
        if (weight_sum <= 0.0_dp .or. .not. ieee_is_finite(weight_sum)) status = 5
      end if
    end if
  end subroutine compute_weight_sum_only

  pure logical function buffer_shape_matches_grid(buffer, grid) result(matches)
    type(fibre_prod_force_buffer_type), intent(in) :: buffer
    type(fibre_prod_grid_type), intent(in) :: grid

    matches = buffer%nx_local == grid%nx_local .and. buffer%ny_local == grid%ny_local .and. &
              buffer%nz_local == grid%nz_local
  end function buffer_shape_matches_grid

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

end module fibre_prod_ibm_spreading
