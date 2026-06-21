program fibre_prod_grid_adapter_check
  use fibre_prod_grid_adapter, only : dp, fibre_prod_grid_type, &
                                      fibre_prod_grid_init_from_coordinates, &
                                      fibre_prod_grid_validate, &
                                      fibre_prod_grid_destroy, &
                                      fibre_prod_grid_is_initialized, &
                                      fibre_prod_grid_compute_min_spacing, &
                                      fibre_prod_grid_compute_max_spacing, &
                                      fibre_prod_grid_compute_total_local_volume, &
                                      fibre_prod_grid_find_cell, &
                                      fibre_prod_grid_bounds_string
  implicit none

  type(fibre_prod_grid_type) :: grid
  real(dp) :: x_coord(6)
  real(dp) :: y_coord(5)
  real(dp) :: z_coord(4)
  real(dp) :: point(3)
  real(dp) :: min_spacing
  real(dp) :: max_spacing
  real(dp) :: total_volume
  integer :: status
  integer :: i_cell
  integer :: j_cell
  integer :: k_cell
  character(len=256) :: summary

  x_coord = [0.0_dp, 0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp, 1.0_dp]
  y_coord = [0.0_dp, 0.04_dp, 0.16_dp, 0.49_dp, 1.0_dp]
  z_coord = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp]

  call fibre_prod_grid_init_from_coordinates(grid, x_coord, y_coord, z_coord, &
                                             2, 5, 1, 4, 2, 4, &
                                             .true., .false., .true., status)
  if (status /= 0) then
    print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK FAIL: init status', status
    error stop 1
  end if

  if (.not. fibre_prod_grid_is_initialized(grid)) then
    print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK FAIL: grid not initialized'
    error stop 2
  end if

  if (.not. fibre_prod_grid_validate(grid)) then
    print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK FAIL: grid validation failed'
    error stop 3
  end if

  if (grid%nx_local /= 4 .or. grid%ny_local /= 4 .or. grid%nz_local /= 3) then
    print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK FAIL: local sizes mismatch'
    error stop 4
  end if

  if (.not. all(grid%dx > 0.0_dp) .or. .not. all(grid%dy > 0.0_dp) .or. &
      .not. all(grid%dz > 0.0_dp)) then
    print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK FAIL: non-positive spacing'
    error stop 5
  end if

  if (.not. all(grid%cell_volume > 0.0_dp)) then
    print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK FAIL: non-positive cell volume'
    error stop 6
  end if

  min_spacing = fibre_prod_grid_compute_min_spacing(grid)
  max_spacing = fibre_prod_grid_compute_max_spacing(grid)
  if (min_spacing <= 0.0_dp .or. max_spacing < min_spacing) then
    print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK FAIL: spacing extrema invalid', min_spacing, max_spacing
    error stop 7
  end if

  total_volume = fibre_prod_grid_compute_total_local_volume(grid)
  if (total_volume <= 0.0_dp) then
    print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK FAIL: total local volume invalid', total_volume
    error stop 8
  end if

  point = [0.5_dp, 0.20_dp, 0.60_dp]
  call fibre_prod_grid_find_cell(grid, point, i_cell, j_cell, k_cell, status)
  if (status /= 0 .or. i_cell < grid%istart .or. i_cell >= grid%iend .or. &
      j_cell < grid%jstart .or. j_cell >= grid%jend .or. &
      k_cell < grid%kstart .or. k_cell >= grid%kend) then
    print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK FAIL: internal point cell lookup failed', &
             status, i_cell, j_cell, k_cell
    error stop 9
  end if

  point = [-1.0_dp, 0.20_dp, 0.60_dp]
  call fibre_prod_grid_find_cell(grid, point, i_cell, j_cell, k_cell, status)
  if (status == 0) then
    print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK FAIL: out-of-domain point accepted'
    error stop 10
  end if

  summary = fibre_prod_grid_bounds_string(grid)
  if (len_trim(summary) == 0) then
    print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK FAIL: empty bounds summary'
    error stop 11
  end if

  call fibre_prod_grid_destroy(grid)
  if (grid%initialized .or. allocated(grid%x) .or. allocated(grid%y) .or. allocated(grid%z) .or. &
      allocated(grid%dx) .or. allocated(grid%dy) .or. allocated(grid%dz) .or. &
      allocated(grid%cell_volume)) then
    print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK FAIL: destroy did not deallocate arrays'
    error stop 12
  end if

  print *, 'R3_FIBRE_PROD_GRID_ADAPTER_CHECK PASS'
end program fibre_prod_grid_adapter_check
