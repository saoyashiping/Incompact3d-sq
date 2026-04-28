program fibre_ibm_spreading_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_uniform_ibm_grid, destroy_ibm_grid
  use fibre_ibm_grid, only : destroy_lagrangian_points
  use fibre_ibm_spreading, only : clear_eulerian_force_density, spread_lag_force_to_eulerian
  use fibre_ibm_spreading, only : compute_eulerian_total_force, compute_lagrangian_total_force
  use fibre_ibm_spreading, only : compute_force_density_norms

  implicit none

  type(ibm_grid_t) :: grid
  type(ibm_lagrangian_points_t) :: lag

  integer, parameter :: nlag = 5
  integer :: i, j, k, l, io_unit
  integer :: single_lag_nonzero_grid_count
  integer :: spreading_status

  real(mytype), allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
  real(mytype), allocatable :: fx_a(:,:,:), fy_a(:,:,:), fz_a(:,:,:)
  real(mytype), allocatable :: fx_b(:,:,:), fy_b(:,:,:), fz_b(:,:,:)
  real(mytype), allocatable :: fx_ab(:,:,:), fy_ab(:,:,:), fz_ab(:,:,:)
  real(mytype), allocatable :: force_a(:,:), force_b(:,:)

  real(mytype) :: lag_total_force(3), eul_total_force(3)
  real(mytype) :: max_abs_fd, l2_fd
  real(mytype) :: ds_l, ref_norm

  real(mytype) :: zero_force_grid_norm, zero_force_total_force_norm
  real(mytype) :: single_lag_total_force_error_norm, single_lag_relative_force_error
  real(mytype) :: single_lag_max_abs_force_density
  real(mytype) :: multi_lag_total_force_error_norm, multi_lag_relative_force_error
  real(mytype) :: superposition_max_error, superposition_l2_error
  real(mytype) :: periodic_lag_total_force_error_norm, periodic_lag_relative_force_error
  real(mytype) :: clear_grid_force_norm_after_clear

  call init_uniform_ibm_grid(grid, nx=16, ny=12, nz=10, &
                             xmin=0._mytype, xmax=2._mytype, &
                             ymin=-1._mytype, ymax=1._mytype, &
                             zmin=0._mytype, zmax=1._mytype, &
                             periodic_x=.true., periodic_y=.false., periodic_z=.true.)

  allocate(lag%x(3, nlag), lag%v(3, nlag), lag%force(3, nlag), lag%weight(nlag))
  lag%nl = nlag
  lag%x(:, 1) = [0.73_mytype,  0.17_mytype, 0.41_mytype]
  lag%x(:, 2) = [1.22_mytype, -0.25_mytype, 0.62_mytype]
  lag%x(:, 3) = [1.00_mytype,  0.00_mytype, 0.50_mytype]
  lag%x(:, 4) = [1.85_mytype, -0.10_mytype, 0.12_mytype]
  lag%x(:, 5) = [0.05_mytype,  0.00_mytype, 0.95_mytype]
  lag%v = 0._mytype
  lag%force = 0._mytype
  ds_l = 0.03125_mytype
  lag%weight = ds_l

  allocate(fx(grid%nx, grid%ny, grid%nz), fy(grid%nx, grid%ny, grid%nz), fz(grid%nx, grid%ny, grid%nz))
  allocate(fx_a(grid%nx, grid%ny, grid%nz), fy_a(grid%nx, grid%ny, grid%nz), fz_a(grid%nx, grid%ny, grid%nz))
  allocate(fx_b(grid%nx, grid%ny, grid%nz), fy_b(grid%nx, grid%ny, grid%nz), fz_b(grid%nx, grid%ny, grid%nz))
  allocate(fx_ab(grid%nx, grid%ny, grid%nz), fy_ab(grid%nx, grid%ny, grid%nz), fz_ab(grid%nx, grid%ny, grid%nz))
  allocate(force_a(3, nlag), force_b(3, nlag))

  call clear_eulerian_force_density(fx, fy, fz)
  lag%force = 0._mytype
  call spread_lag_force_to_eulerian(grid, lag, fx, fy, fz)
  call compute_force_density_norms(fx, fy, fz, max_abs_fd, l2_fd)
  zero_force_grid_norm = l2_fd
  call compute_eulerian_total_force(grid, fx, fy, fz, eul_total_force)
  zero_force_total_force_norm = sqrt(sum(eul_total_force**2))

  call clear_eulerian_force_density(fx, fy, fz)
  lag%force = 0._mytype
  lag%force(:, 1) = [0.2_mytype, -0.1_mytype, 0.05_mytype]
  call spread_lag_force_to_eulerian(grid, lag, fx, fy, fz)
  call compute_lagrangian_total_force(lag, lag_total_force)
  call compute_eulerian_total_force(grid, fx, fy, fz, eul_total_force)
  single_lag_total_force_error_norm = sqrt(sum((eul_total_force - lag_total_force)**2))
  ref_norm = max(sqrt(sum(lag_total_force**2)), 1.0e-30_mytype)
  single_lag_relative_force_error = single_lag_total_force_error_norm / ref_norm

  call compute_force_density_norms(fx, fy, fz, single_lag_max_abs_force_density, l2_fd)
  single_lag_nonzero_grid_count = 0
  do k = 1, grid%nz
    do j = 1, grid%ny
      do i = 1, grid%nx
        if (abs(fx(i, j, k)) > 1.0e-30_mytype .or. abs(fy(i, j, k)) > 1.0e-30_mytype .or. abs(fz(i, j, k)) > 1.0e-30_mytype) then
          single_lag_nonzero_grid_count = single_lag_nonzero_grid_count + 1
        end if
      end do
    end do
  end do

  call clear_eulerian_force_density(fx, fy, fz)
  do l = 1, nlag
    lag%force(1, l) = 0.1_mytype + 0.02_mytype * real(l, mytype)
    lag%force(2, l) = -0.05_mytype + 0.01_mytype * real(l, mytype)
    lag%force(3, l) = 0.03_mytype - 0.005_mytype * real(l, mytype)
  end do
  call spread_lag_force_to_eulerian(grid, lag, fx, fy, fz)
  call compute_lagrangian_total_force(lag, lag_total_force)
  call compute_eulerian_total_force(grid, fx, fy, fz, eul_total_force)
  multi_lag_total_force_error_norm = sqrt(sum((eul_total_force - lag_total_force)**2))
  ref_norm = max(sqrt(sum(lag_total_force**2)), 1.0e-30_mytype)
  multi_lag_relative_force_error = multi_lag_total_force_error_norm / ref_norm

  do l = 1, nlag
    force_a(1, l) = 0.05_mytype + 0.01_mytype * real(l, mytype)
    force_a(2, l) = -0.02_mytype * real(l, mytype)
    force_a(3, l) = 0.01_mytype

    force_b(1, l) = -0.03_mytype + 0.005_mytype * real(l, mytype)
    force_b(2, l) = 0.04_mytype - 0.003_mytype * real(l, mytype)
    force_b(3, l) = -0.02_mytype + 0.002_mytype * real(l, mytype)
  end do

  call clear_eulerian_force_density(fx_a, fy_a, fz_a)
  lag%force = force_a
  call spread_lag_force_to_eulerian(grid, lag, fx_a, fy_a, fz_a)

  call clear_eulerian_force_density(fx_b, fy_b, fz_b)
  lag%force = force_b
  call spread_lag_force_to_eulerian(grid, lag, fx_b, fy_b, fz_b)

  call clear_eulerian_force_density(fx_ab, fy_ab, fz_ab)
  lag%force = force_a + force_b
  call spread_lag_force_to_eulerian(grid, lag, fx_ab, fy_ab, fz_ab)

  superposition_max_error = maxval(abs(fx_ab - fx_a - fx_b))
  superposition_max_error = max(superposition_max_error, maxval(abs(fy_ab - fy_a - fy_b)))
  superposition_max_error = max(superposition_max_error, maxval(abs(fz_ab - fz_a - fz_b)))

  superposition_l2_error = sqrt(sum((fx_ab - fx_a - fx_b)**2 + (fy_ab - fy_a - fy_b)**2 + (fz_ab - fz_a - fz_b)**2))

  call clear_eulerian_force_density(fx, fy, fz)
  lag%force = 0._mytype
  lag%force(:, 5) = [0.13_mytype, 0.07_mytype, -0.02_mytype]
  call spread_lag_force_to_eulerian(grid, lag, fx, fy, fz)
  call compute_lagrangian_total_force(lag, lag_total_force)
  call compute_eulerian_total_force(grid, fx, fy, fz, eul_total_force)
  periodic_lag_total_force_error_norm = sqrt(sum((eul_total_force - lag_total_force)**2))
  ref_norm = max(sqrt(sum(lag_total_force**2)), 1.0e-30_mytype)
  periodic_lag_relative_force_error = periodic_lag_total_force_error_norm / ref_norm

  fx = 1._mytype
  fy = -2._mytype
  fz = 0.5_mytype
  call clear_eulerian_force_density(fx, fy, fz)
  call compute_force_density_norms(fx, fy, fz, max_abs_fd, l2_fd)
  clear_grid_force_norm_after_clear = l2_fd

  spreading_status = 1
  if (ieee_is_nan(zero_force_grid_norm) .or. ieee_is_nan(zero_force_total_force_norm) .or. &
      ieee_is_nan(single_lag_total_force_error_norm) .or. ieee_is_nan(single_lag_relative_force_error) .or. &
      ieee_is_nan(single_lag_max_abs_force_density) .or. ieee_is_nan(multi_lag_total_force_error_norm) .or. &
      ieee_is_nan(multi_lag_relative_force_error) .or. ieee_is_nan(superposition_max_error) .or. &
      ieee_is_nan(superposition_l2_error) .or. ieee_is_nan(periodic_lag_total_force_error_norm) .or. &
      ieee_is_nan(periodic_lag_relative_force_error) .or. ieee_is_nan(clear_grid_force_norm_after_clear)) then
    spreading_status = 0
  end if

  open(newunit=io_unit, file='stage3_outputs/fibre_ibm_spreading_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,ES24.16)') 'zero_force_grid_norm = ', zero_force_grid_norm
  write(io_unit, '(A,ES24.16)') 'zero_force_total_force_norm = ', zero_force_total_force_norm

  write(io_unit, '(A,ES24.16)') 'single_lag_total_force_error_norm = ', single_lag_total_force_error_norm
  write(io_unit, '(A,ES24.16)') 'single_lag_relative_force_error = ', single_lag_relative_force_error
  write(io_unit, '(A,I0)') 'single_lag_nonzero_grid_count = ', single_lag_nonzero_grid_count
  write(io_unit, '(A,ES24.16)') 'single_lag_max_abs_force_density = ', single_lag_max_abs_force_density

  write(io_unit, '(A,ES24.16)') 'multi_lag_total_force_error_norm = ', multi_lag_total_force_error_norm
  write(io_unit, '(A,ES24.16)') 'multi_lag_relative_force_error = ', multi_lag_relative_force_error

  write(io_unit, '(A,ES24.16)') 'superposition_max_error = ', superposition_max_error
  write(io_unit, '(A,ES24.16)') 'superposition_l2_error = ', superposition_l2_error

  write(io_unit, '(A,ES24.16)') 'periodic_lag_total_force_error_norm = ', periodic_lag_total_force_error_norm
  write(io_unit, '(A,ES24.16)') 'periodic_lag_relative_force_error = ', periodic_lag_relative_force_error

  write(io_unit, '(A,ES24.16)') 'clear_grid_force_norm_after_clear = ', clear_grid_force_norm_after_clear

  write(io_unit, '(A,I0)') 'spreading_status = ', spreading_status

  close(io_unit)

  deallocate(fx, fy, fz, fx_a, fy_a, fz_a, fx_b, fy_b, fz_b, fx_ab, fy_ab, fz_ab, force_a, force_b)
  call destroy_lagrangian_points(lag)
  call destroy_ibm_grid(grid)

end program fibre_ibm_spreading_check
