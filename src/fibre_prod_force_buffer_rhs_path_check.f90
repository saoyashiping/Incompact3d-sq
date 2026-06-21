program fibre_prod_force_buffer_rhs_path_check
  use, intrinsic :: iso_fortran_env, only : real64
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, &
                                      fibre_prod_grid_destroy
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, &
                                          fibre_prod_force_buffer_allocate, &
                                          fibre_prod_force_buffer_reset_to_zero, &
                                          fibre_prod_force_buffer_destroy, &
                                          fibre_prod_force_buffer_total_force
  use fibre_prod_ibm_spreading, only : fibre_prod_spread_point_force
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, fibre_prod_runtime_config_default
  use fibre_prod_main_hook, only : fibre_prod_main_hook_init, fibre_prod_main_hook_apply_force_buffer
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: n = 7
  real(dp), parameter :: tol = 1.0e-10_dp
  type(fibre_prod_grid_type) :: grid
  type(fibre_prod_force_buffer_type) :: buffer
  type(fibre_prod_force_buffer_type) :: invalid_buffer
  type(fibre_prod_runtime_config_type) :: lambda0_config
  real(dp) :: x(n)
  real(dp) :: y(n)
  real(dp) :: z(n)
  real(dp) :: point(3)
  real(dp) :: force(3)
  real(dp) :: total_force(3)
  real(dp) :: rhs_total_force(3)
  real(dp) :: rhs_x(n, n, n)
  real(dp) :: rhs_y(n, n, n)
  real(dp) :: rhs_z(n, n, n)
  real(dp) :: baseline_x(n, n, n)
  real(dp) :: baseline_y(n, n, n)
  real(dp) :: baseline_z(n, n, n)
  real(dp) :: lambda_fsi
  real(dp) :: penalty_beta
  real(dp) :: force_scale
  real(dp) :: weight_sum
  integer :: status
  integer :: i

  do i = 1, n
    x(i) = real(i - 1, dp) / real(n - 1, dp)
    y(i) = real(i - 1, dp) / real(n - 1, dp)
    z(i) = real(i - 1, dp) / real(n - 1, dp)
  end do

  call fibre_prod_grid_init_from_coordinates(grid, x, y, z, 1, n, 1, n, 1, n, .true., .true., .true., status)
  if (status /= 0) error stop 1
  call fibre_prod_force_buffer_allocate(buffer, grid, status)
  if (status /= 0) error stop 2
  call fibre_prod_force_buffer_reset_to_zero(buffer, status)
  if (status /= 0) error stop 3

  point = [0.5_dp, 0.5_dp, 0.5_dp]
  force = [1.0_dp, -0.5_dp, 0.25_dp]
  call fibre_prod_spread_point_force(grid, buffer, point, force, status, weight_sum)
  if (status /= 0) error stop 4
  if (weight_sum <= 0.0_dp) error stop 5
  call fibre_prod_force_buffer_total_force(buffer, grid, total_force, status)
  if (status /= 0) error stop 6
  if (maxval(abs(total_force - force)) > tol) error stop 7

  rhs_x = 0.0_dp
  rhs_y = 0.0_dp
  rhs_z = 0.0_dp
  baseline_x = rhs_x
  baseline_y = rhs_y
  baseline_z = rhs_z
  call fibre_prod_main_hook_init(status)
  if (status /= 0) error stop 8
  call read_env_real('FIBRE_PROD_LAMBDA', lambda_fsi, status)
  if (status /= 0) error stop 9
  call read_env_real('FIBRE_PROD_PENALTY_BETA', penalty_beta, status)
  if (status /= 0) error stop 10
  force_scale = lambda_fsi * penalty_beta
  call fibre_prod_main_hook_apply_force_buffer(rhs_x, rhs_y, rhs_z, buffer, status)
  if (status /= 0) error stop 11
  if (maxval(abs(rhs_x - force_scale * buffer%fx)) > tol) error stop 12
  if (maxval(abs(rhs_y - force_scale * buffer%fy)) > tol) error stop 13
  if (maxval(abs(rhs_z - force_scale * buffer%fz)) > tol) error stop 14
  if (maxval(abs(rhs_x)) <= 0.0_dp .and. maxval(abs(rhs_y)) <= 0.0_dp .and. maxval(abs(rhs_z)) <= 0.0_dp) error stop 15
  if (is_uniform_rhs(rhs_x, rhs_y, rhs_z)) error stop 16
  call rhs_integral(grid, rhs_x, rhs_y, rhs_z, rhs_total_force)
  if (maxval(abs(rhs_total_force - force_scale * total_force)) > tol) error stop 17

  call fibre_prod_runtime_config_default(lambda0_config)
  lambda0_config%enabled = .true.
  lambda0_config%lambda_fsi = 0.0_dp
  lambda0_config%penalty_beta = penalty_beta
  rhs_x = baseline_x
  rhs_y = baseline_y
  rhs_z = baseline_z
  call fibre_prod_main_hook_init(status, lambda0_config)
  if (status /= 0) error stop 18
  call fibre_prod_main_hook_apply_force_buffer(rhs_x, rhs_y, rhs_z, buffer, status)
  if (status /= 0) error stop 19
  if (maxval(abs(rhs_x - baseline_x)) > 0.0_dp) error stop 20
  if (maxval(abs(rhs_y - baseline_y)) > 0.0_dp) error stop 21
  if (maxval(abs(rhs_z - baseline_z)) > 0.0_dp) error stop 22

  rhs_x = baseline_x
  rhs_y = baseline_y
  rhs_z = baseline_z
  call fibre_prod_main_hook_init(status)
  if (status /= 0) error stop 23
  call fibre_prod_main_hook_apply_force_buffer(rhs_x, rhs_y, rhs_z, invalid_buffer, status)
  if (status == 0) error stop 24
  if (maxval(abs(rhs_x - baseline_x)) > 0.0_dp) error stop 25
  if (maxval(abs(rhs_y - baseline_y)) > 0.0_dp) error stop 26
  if (maxval(abs(rhs_z - baseline_z)) > 0.0_dp) error stop 27

  call fibre_prod_force_buffer_destroy(buffer)
  call fibre_prod_grid_destroy(grid)
  print *, 'P0_2_FORCE_BUFFER_TO_RHS_PATH_CHECK PASS'
contains

  subroutine read_env_real(name, value, status)
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    character(len=128) :: raw
    integer :: length

    value = 0.0_dp
    call get_environment_variable(name, raw, length=length, status=status)
    if (status /= 0 .or. length <= 0) then
      status = 1
      return
    end if
    read(raw(1:min(length, len(raw))), *, iostat=status) value
  end subroutine read_env_real

  subroutine rhs_integral(grid, rhs_x, rhs_y, rhs_z, total)
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(in) :: rhs_x(:, :, :)
    real(dp), intent(in) :: rhs_y(:, :, :)
    real(dp), intent(in) :: rhs_z(:, :, :)
    real(dp), intent(out) :: total(3)

    total(1) = sum(rhs_x * grid%cell_volume)
    total(2) = sum(rhs_y * grid%cell_volume)
    total(3) = sum(rhs_z * grid%cell_volume)
  end subroutine rhs_integral

  logical function is_uniform_rhs(rhs_x, rhs_y, rhs_z) result(is_uniform)
    real(dp), intent(in) :: rhs_x(:, :, :)
    real(dp), intent(in) :: rhs_y(:, :, :)
    real(dp), intent(in) :: rhs_z(:, :, :)

    is_uniform = all(rhs_x == rhs_x(1, 1, 1)) .and. all(rhs_y == rhs_y(1, 1, 1)) .and. &
                 all(rhs_z == rhs_z(1, 1, 1)) .and. rhs_x(1, 1, 1) == rhs_y(1, 1, 1) .and. &
                 rhs_y(1, 1, 1) == rhs_z(1, 1, 1)
  end function is_uniform_rhs
end program fibre_prod_force_buffer_rhs_path_check
