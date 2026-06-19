program fibre_prod_ibm_spreading_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, &
                                      fibre_prod_grid_destroy
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, &
                                          fibre_prod_force_buffer_allocate, &
                                          fibre_prod_force_buffer_reset_to_zero, &
                                          fibre_prod_force_buffer_destroy, &
                                          fibre_prod_force_buffer_is_finite, &
                                          fibre_prod_force_buffer_total_force
  use fibre_prod_ibm_spreading, only : fibre_prod_spread_point_force, &
                                       fibre_prod_spread_multiple_point_forces, &
                                       fibre_prod_spreading_weight_sum
  implicit none

  integer, parameter :: dp = real64
  type(fibre_prod_grid_type) :: grid
  type(fibre_prod_force_buffer_type) :: buffer
  real(dp) :: x_coord(8)
  real(dp) :: y_coord(6)
  real(dp) :: z_coord(7)
  real(dp) :: point(3)
  real(dp) :: force(3)
  real(dp) :: total_force(3)
  real(dp) :: points(2, 3)
  real(dp) :: forces(2, 3)
  real(dp) :: weight_sum
  real(dp) :: power_lagrangian
  real(dp) :: power_eulerian
  real(dp) :: velocity(3)
  integer :: status

  x_coord = [0.0_dp, 0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp, 1.0_dp, 1.2_dp, 1.4_dp]
  y_coord = [0.0_dp, 0.03_dp, 0.12_dp, 0.30_dp, 0.62_dp, 1.0_dp]
  z_coord = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp, 1.25_dp, 1.5_dp]

  call fibre_prod_grid_init_from_coordinates(grid, x_coord, y_coord, z_coord, &
                                             1, 8, 1, 6, 1, 7, &
                                             .true., .false., .true., status)
  if (status /= 0) error stop 1

  call fibre_prod_force_buffer_allocate(buffer, grid, status)
  if (status /= 0 .or. .not. buffer%allocated) error stop 2
  call fibre_prod_force_buffer_reset_to_zero(buffer, status)
  if (status /= 0 .or. any(buffer%fx /= 0.0_dp) .or. any(buffer%fy /= 0.0_dp) .or. &
      any(buffer%fz /= 0.0_dp)) error stop 3

  point = [0.35_dp, 0.20_dp, 0.65_dp]
  force = [1.0_dp, -2.0_dp, 0.5_dp]
  call fibre_prod_spread_point_force(grid, buffer, point, force, status, weight_sum)
  if (status /= 0 .or. weight_sum <= 0.0_dp) error stop 4
  call fibre_prod_force_buffer_total_force(buffer, grid, total_force, status)
  if (status /= 0 .or. maxval(abs(total_force - force)) > 1.0e-12_dp) error stop 5

  call fibre_prod_force_buffer_reset_to_zero(buffer, status)
  points(1, :) = [0.35_dp, 0.20_dp, 0.65_dp]
  points(2, :) = [0.85_dp, 0.40_dp, 1.10_dp]
  forces(1, :) = [1.0_dp, -2.0_dp, 0.5_dp]
  forces(2, :) = [-0.25_dp, 0.75_dp, 1.25_dp]
  call fibre_prod_spread_multiple_point_forces(grid, buffer, points, forces, status, weight_sum)
  if (status /= 0 .or. weight_sum <= 0.0_dp) error stop 6
  call fibre_prod_force_buffer_total_force(buffer, grid, total_force, status)
  if (status /= 0 .or. maxval(abs(total_force - sum(forces, dim=1))) > 1.0e-12_dp) error stop 7

  call fibre_prod_force_buffer_reset_to_zero(buffer, status)
  point = [-0.05_dp, 0.20_dp, 0.65_dp]
  call fibre_prod_spread_point_force(grid, buffer, point, force, status, weight_sum)
  if (status /= 0) error stop 8

  call fibre_prod_force_buffer_reset_to_zero(buffer, status)
  point = [0.35_dp, 0.20_dp, 1.55_dp]
  call fibre_prod_spread_point_force(grid, buffer, point, force, status, weight_sum)
  if (status /= 0) error stop 9

  point = [0.35_dp, -0.01_dp, 0.65_dp]
  call fibre_prod_spread_point_force(grid, buffer, point, force, status, weight_sum)
  if (status == 0) error stop 10

  point = [0.35_dp, 1.01_dp, 0.65_dp]
  call fibre_prod_spread_point_force(grid, buffer, point, force, status, weight_sum)
  if (status == 0) error stop 11

  point = [ieee_value(0.0_dp, ieee_quiet_nan), 0.20_dp, 0.65_dp]
  call fibre_prod_spread_point_force(grid, buffer, point, force, status, weight_sum)
  if (status == 0) error stop 12

  point = [0.35_dp, 0.20_dp, 0.65_dp]
  call fibre_prod_spreading_weight_sum(grid, point, weight_sum, status)
  if (status /= 0 .or. weight_sum <= 0.0_dp) error stop 13
  if (.not. fibre_prod_force_buffer_is_finite(buffer)) error stop 14

  call fibre_prod_force_buffer_reset_to_zero(buffer, status)
  call fibre_prod_force_buffer_total_force(buffer, grid, total_force, status)
  if (status /= 0 .or. maxval(abs(total_force)) > 1.0e-12_dp) error stop 15

  force = [1.25_dp, -0.50_dp, 0.75_dp]
  velocity = [2.0_dp, -1.0_dp, 0.5_dp]
  point = [0.35_dp, 0.20_dp, 0.65_dp]
  call fibre_prod_spread_point_force(grid, buffer, point, force, status, weight_sum)
  if (status /= 0) error stop 16
  power_lagrangian = sum(force * velocity)
  power_eulerian = sum((buffer%fx * velocity(1) + buffer%fy * velocity(2) + buffer%fz * velocity(3)) * &
                       grid%cell_volume)
  if (abs(power_lagrangian - power_eulerian) > 1.0e-12_dp) error stop 17

  call fibre_prod_force_buffer_destroy(buffer)
  if (buffer%allocated .or. allocated(buffer%fx) .or. allocated(buffer%fy) .or. allocated(buffer%fz)) error stop 18
  call fibre_prod_grid_destroy(grid)
  if (grid%initialized .or. allocated(grid%x) .or. allocated(grid%y) .or. allocated(grid%z) .or. &
      allocated(grid%cell_volume)) error stop 19

  print *, 'R6_FIBRE_PROD_IBM_SPREADING_CHECK PASS'
end program fibre_prod_ibm_spreading_check
