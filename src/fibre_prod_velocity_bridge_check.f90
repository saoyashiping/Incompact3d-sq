program fibre_prod_velocity_bridge_check
  use, intrinsic :: iso_fortran_env, only : real64
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, fibre_prod_grid_destroy
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, fibre_prod_runtime_config_default
  use fibre_prod_main_hook, only : fibre_prod_main_hook_init
  use fibre_prod_runtime_bridge, only : fibre_prod_runtime_bridge_type, &
                                        fibre_prod_runtime_bridge_init_from_rhs, &
                                        fibre_prod_runtime_bridge_apply_lambda0_noop, &
                                        fibre_prod_runtime_bridge_finalize
  use fibre_prod_velocity_bridge, only : fibre_prod_velocity_sample_type, &
                                         fibre_prod_velocity_sample_allocate, &
                                         fibre_prod_velocity_sample_set_points, &
                                         fibre_prod_velocity_bridge_sample, &
                                         fibre_prod_velocity_sample_finalize
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: n = 9
  real(dp), parameter :: tol = 1.0e-8_dp
  type(fibre_prod_grid_type) :: grid
  type(fibre_prod_velocity_sample_type) :: sample
  type(fibre_prod_runtime_bridge_type) :: bridge
  type(fibre_prod_runtime_config_type) :: config
  real(dp) :: x(n), y(n), z(n)
  real(dp) :: ux(n, n, n), uy(n, n, n), uz(n, n, n)
  real(dp) :: ux0(n, n, n), uy0(n, n, n), uz0(n, n, n)
  real(dp) :: rhs_x(n, n, n), rhs_y(n, n, n), rhs_z(n, n, n)
  real(dp) :: rhs_x0(n, n, n), rhs_y0(n, n, n), rhs_z0(n, n, n)
  real(dp) :: points(3, 3)
  real(dp) :: bad_points(1, 3)
  real(dp) :: expected(3)
  integer :: status
  integer :: ipoint

  call initialize_grid_coordinates(x, y, z)
  call fibre_prod_grid_init_from_coordinates(grid, x, y, z, 1, n, 1, n, 1, n, .false., .false., .false., status)
  if (status /= 0) error stop 1
  call initialize_velocity_field(grid, ux, uy, uz)
  ux0 = ux
  uy0 = uy
  uz0 = uz
  call initialize_rhs(rhs_x, rhs_y, rhs_z)
  rhs_x0 = rhs_x
  rhs_y0 = rhs_y
  rhs_z0 = rhs_z

  call fibre_prod_velocity_sample_allocate(sample, 3, status)
  if (status /= 0) error stop 2
  points(1, :) = [0.5_dp, 0.5_dp, 0.5_dp]
  points(2, :) = [0.375_dp, 0.375_dp, 0.375_dp]
  points(3, :) = [0.625_dp, 0.5_dp, 0.375_dp]
  call fibre_prod_velocity_sample_set_points(sample, points, status)
  if (status /= 0) error stop 3
  call fibre_prod_velocity_bridge_sample(grid, ux, uy, uz, sample, status)
  if (status /= 0) error stop 4
  do ipoint = 1, sample%npoint
    expected = analytic_velocity(sample%x(ipoint, :))
    if (maxval(abs(sample%u(ipoint, :) - expected)) > tol) error stop 5
  end do
  if (.not. same_field(ux, ux0) .or. .not. same_field(uy, uy0) .or. .not. same_field(uz, uz0)) error stop 6
  if (.not. same_field(rhs_x, rhs_x0) .or. .not. same_field(rhs_y, rhs_y0) .or. &
      .not. same_field(rhs_z, rhs_z0)) error stop 7

  call fibre_prod_runtime_bridge_init_from_rhs(rhs_x, rhs_y, rhs_z, bridge, status)
  if (status /= 0) error stop 8
  call fibre_prod_runtime_config_default(config)
  config%enabled = .true.
  config%lambda_fsi = 0.0_dp
  config%penalty_beta = 2.0_dp
  call fibre_prod_main_hook_init(status, config)
  if (status /= 0) error stop 9
  call fibre_prod_runtime_bridge_apply_lambda0_noop(rhs_x, rhs_y, rhs_z, bridge, status)
  if (status /= 0) error stop 10
  call fibre_prod_velocity_bridge_sample(grid, ux, uy, uz, sample, status)
  if (status /= 0) error stop 11
  if (.not. same_field(ux, ux0) .or. .not. same_field(uy, uy0) .or. .not. same_field(uz, uz0)) error stop 12
  if (.not. same_field(rhs_x, rhs_x0) .or. .not. same_field(rhs_y, rhs_y0) .or. &
      .not. same_field(rhs_z, rhs_z0)) error stop 13

  config%lambda_fsi = 1.0e-3_dp
  call fibre_prod_main_hook_init(status, config)
  if (status /= 0) error stop 14
  call fibre_prod_velocity_bridge_sample(grid, ux, uy, uz, sample, status)
  if (status /= 0) error stop 15
  if (.not. same_field(rhs_x, rhs_x0) .or. .not. same_field(rhs_y, rhs_y0) .or. &
      .not. same_field(rhs_z, rhs_z0)) error stop 16

  bad_points(1, :) = [1.25_dp, 0.5_dp, 0.5_dp]
  call fibre_prod_velocity_sample_finalize(sample)
  call fibre_prod_velocity_sample_allocate(sample, 1, status)
  if (status /= 0) error stop 17
  call fibre_prod_velocity_sample_set_points(sample, bad_points, status)
  if (status /= 0) error stop 18
  call fibre_prod_velocity_bridge_sample(grid, ux, uy, uz, sample, status)
  if (status == 0) error stop 19
  if (.not. same_field(ux, ux0) .or. .not. same_field(uy, uy0) .or. .not. same_field(uz, uz0)) error stop 20
  if (.not. same_field(rhs_x, rhs_x0) .or. .not. same_field(rhs_y, rhs_y0) .or. &
      .not. same_field(rhs_z, rhs_z0)) error stop 21

  call fibre_prod_velocity_sample_finalize(sample)
  call fibre_prod_runtime_bridge_finalize(bridge)
  call fibre_prod_grid_destroy(grid)
  print *, 'P0_4_VELOCITY_BRIDGE_CHECK PASS'
contains

  subroutine initialize_grid_coordinates(x, y, z)
    real(dp), intent(out) :: x(:), y(:), z(:)
    integer :: i

    do i = 1, size(x)
      x(i) = real(i - 1, dp) / real(size(x) - 1, dp)
      y(i) = real(i - 1, dp) / real(size(y) - 1, dp)
      z(i) = real(i - 1, dp) / real(size(z) - 1, dp)
    end do
  end subroutine initialize_grid_coordinates

  subroutine initialize_velocity_field(grid, ux, uy, uz)
    type(fibre_prod_grid_type), intent(in) :: grid
    real(dp), intent(out) :: ux(:, :, :), uy(:, :, :), uz(:, :, :)
    integer :: i, j, k
    real(dp) :: point(3)
    real(dp) :: velocity(3)

    do k = 1, size(ux, 3)
      do j = 1, size(ux, 2)
        do i = 1, size(ux, 1)
          point = [grid%x(i), grid%y(j), grid%z(k)]
          velocity = analytic_velocity(point)
          ux(i, j, k) = velocity(1)
          uy(i, j, k) = velocity(2)
          uz(i, j, k) = velocity(3)
        end do
      end do
    end do
  end subroutine initialize_velocity_field

  subroutine initialize_rhs(rhs_x, rhs_y, rhs_z)
    real(dp), intent(out) :: rhs_x(:, :, :), rhs_y(:, :, :), rhs_z(:, :, :)
    integer :: i, j, k

    do k = 1, size(rhs_x, 3)
      do j = 1, size(rhs_x, 2)
        do i = 1, size(rhs_x, 1)
          rhs_x(i, j, k) = real(i + j + k, dp)
          rhs_y(i, j, k) = -real(i + 2*j + k, dp)
          rhs_z(i, j, k) = real(i - j + 2*k, dp)
        end do
      end do
    end do
  end subroutine initialize_rhs

  pure function analytic_velocity(point) result(velocity)
    real(dp), intent(in) :: point(3)
    real(dp) :: velocity(3)

    velocity(1) = 1.0_dp + point(1) + 2.0_dp * point(2) + 3.0_dp * point(3)
    velocity(2) = -2.0_dp + 2.0_dp * point(1) - point(2) + 0.5_dp * point(3)
    velocity(3) = 0.25_dp + 0.5_dp * point(1) + point(2) - point(3)
  end function analytic_velocity

  logical function same_field(lhs, rhs) result(matches)
    real(dp), intent(in) :: lhs(:, :, :)
    real(dp), intent(in) :: rhs(:, :, :)

    matches = all(lhs == rhs)
  end function same_field
end program fibre_prod_velocity_bridge_check
