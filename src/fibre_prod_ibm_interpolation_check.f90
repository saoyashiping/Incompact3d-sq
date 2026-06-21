program fibre_prod_ibm_interpolation_check
  use fibre_prod_grid_adapter, only : dp, fibre_prod_grid_type, &
                                      fibre_prod_grid_init_from_coordinates, &
                                      fibre_prod_grid_destroy
  use fibre_prod_ibm_delta, only : fibre_prod_delta_kernel_1d
  use fibre_prod_ibm_interpolation, only : fibre_prod_ibm_interpolate_scalar, &
                                           fibre_prod_ibm_interpolate_velocity
  implicit none

  type(fibre_prod_grid_type) :: grid
  real(dp) :: x_coord(8)
  real(dp) :: y_coord(6)
  real(dp) :: z_coord(7)
  real(dp), allocatable :: scalar_field(:, :, :)
  real(dp), allocatable :: velocity_field(:, :, :, :)
  real(dp) :: point(3)
  real(dp) :: value
  real(dp) :: expected
  real(dp) :: velocity(3)
  real(dp) :: weight_sum
  integer :: status
  integer :: i
  integer :: j
  integer :: k

  x_coord = [0.0_dp, 0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp, 1.0_dp, 1.2_dp, 1.4_dp]
  y_coord = [0.0_dp, 0.03_dp, 0.12_dp, 0.30_dp, 0.62_dp, 1.0_dp]
  z_coord = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp, 1.25_dp, 1.5_dp]

  call fibre_prod_grid_init_from_coordinates(grid, x_coord, y_coord, z_coord, &
                                             1, 8, 1, 6, 1, 7, &
                                             .true., .false., .true., status)
  if (status /= 0) then
    print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK FAIL: grid init', status
    error stop 1
  end if

  if (fibre_prod_delta_kernel_1d(0.0_dp, 0.2_dp) <= 0.0_dp) then
    print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK FAIL: delta center not positive'
    error stop 2
  end if
  if (fibre_prod_delta_kernel_1d(0.41_dp, 0.2_dp) /= 0.0_dp) then
    print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK FAIL: delta support not compact'
    error stop 3
  end if

  allocate(scalar_field(grid%nx_local, grid%ny_local, grid%nz_local))
  allocate(velocity_field(grid%nx_local, grid%ny_local, grid%nz_local, 3))

  scalar_field = 4.25_dp
  velocity_field(:, :, :, 1) = 1.5_dp
  velocity_field(:, :, :, 2) = -2.0_dp
  velocity_field(:, :, :, 3) = 0.75_dp
  point = [0.35_dp, 0.20_dp, 0.65_dp]

  call fibre_prod_ibm_interpolate_scalar(grid, scalar_field, point, value, status, weight_sum)
  if (status /= 0 .or. abs(value - 4.25_dp) > 1.0e-12_dp .or. weight_sum <= 0.0_dp) then
    print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK FAIL: constant scalar', status, value, weight_sum
    error stop 4
  end if

  call fibre_prod_ibm_interpolate_velocity(grid, velocity_field, point, velocity, status, weight_sum)
  if (status /= 0 .or. maxval(abs(velocity - [1.5_dp, -2.0_dp, 0.75_dp])) > 1.0e-12_dp) then
    print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK FAIL: constant velocity', status, velocity
    error stop 5
  end if

  do k = 1, grid%nz_local
    do j = 1, grid%ny_local
      do i = 1, grid%nx_local
        scalar_field(i, j, k) = 2.0_dp + grid%x(i) + 0.5_dp * grid%y(j) - 0.25_dp * grid%z(k)
      end do
    end do
  end do
  expected = 2.0_dp + point(1) + 0.5_dp * point(2) - 0.25_dp * point(3)
  call fibre_prod_ibm_interpolate_scalar(grid, scalar_field, point, value, status, weight_sum)
  if (status /= 0 .or. abs(value - expected) > 0.25_dp) then
    print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK FAIL: affine scalar', status, value, expected
    error stop 6
  end if

  scalar_field = 7.0_dp
  point = [-0.05_dp, 0.20_dp, 0.65_dp]
  call fibre_prod_ibm_interpolate_scalar(grid, scalar_field, point, value, status, weight_sum)
  if (status /= 0 .or. abs(value - 7.0_dp) > 1.0e-12_dp) then
    print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK FAIL: x periodic wrap', status, value
    error stop 7
  end if

  point = [0.35_dp, 0.20_dp, 1.55_dp]
  call fibre_prod_ibm_interpolate_scalar(grid, scalar_field, point, value, status, weight_sum)
  if (status /= 0 .or. abs(value - 7.0_dp) > 1.0e-12_dp) then
    print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK FAIL: z periodic wrap', status, value
    error stop 8
  end if

  point = [0.35_dp, -0.01_dp, 0.65_dp]
  call fibre_prod_ibm_interpolate_scalar(grid, scalar_field, point, value, status, weight_sum)
  if (status == 0) then
    print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK FAIL: lower y out-of-domain accepted'
    error stop 9
  end if

  point = [0.35_dp, 1.01_dp, 0.65_dp]
  call fibre_prod_ibm_interpolate_scalar(grid, scalar_field, point, value, status, weight_sum)
  if (status == 0) then
    print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK FAIL: upper y out-of-domain accepted'
    error stop 10
  end if

  point = [0.35_dp, 0.20_dp, 0.65_dp]
  call fibre_prod_ibm_interpolate_scalar(grid, scalar_field, point, value, status, weight_sum)
  if (status /= 0 .or. weight_sum <= 0.0_dp) then
    print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK FAIL: internal point status/weight', status, weight_sum
    error stop 11
  end if

  call fibre_prod_grid_destroy(grid)
  if (grid%initialized .or. allocated(grid%x) .or. allocated(grid%y) .or. allocated(grid%z) .or. &
      allocated(grid%dx) .or. allocated(grid%dy) .or. allocated(grid%dz) .or. &
      allocated(grid%cell_volume)) then
    print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK FAIL: grid destroy'
    error stop 12
  end if

  deallocate(scalar_field)
  deallocate(velocity_field)
  print *, 'R4_FIBRE_PROD_IBM_INTERPOLATION_CHECK PASS'
end program fibre_prod_ibm_interpolation_check
