program fibre_ibm_power_consistency_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_uniform_ibm_grid, destroy_ibm_grid
  use fibre_ibm_grid, only : destroy_lagrangian_points
  use fibre_ibm_interpolation, only : interpolate_vector_to_lag
  use fibre_ibm_spreading, only : clear_eulerian_force_density, spread_lag_force_to_eulerian
  use fibre_ibm_power_diagnostics, only : compute_eulerian_power, compute_lagrangian_power
  use fibre_ibm_power_diagnostics, only : compute_power_consistency_error

  implicit none

  type(ibm_grid_t) :: grid
  type(ibm_lagrangian_points_t) :: lag

  integer, parameter :: nlag = 5
  integer :: i, j, k, l, io_unit
  integer :: power_consistency_status

  real(mytype), allocatable :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
  real(mytype), allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
  real(mytype), allocatable :: u_lag(:,:)

  real(mytype) :: lx, lz, pi
  real(mytype) :: x, y, z

  real(mytype) :: zero_force_power_eulerian, zero_force_power_lagrangian, zero_force_power_abs_error
  real(mytype) :: unused_rel_error
  real(mytype) :: constant_velocity_power_eulerian, constant_velocity_power_lagrangian
  real(mytype) :: constant_velocity_power_abs_error, constant_velocity_power_rel_error
  real(mytype) :: nonuniform_power_eulerian, nonuniform_power_lagrangian
  real(mytype) :: nonuniform_power_abs_error, nonuniform_power_rel_error
  real(mytype) :: periodic_power_eulerian, periodic_power_lagrangian
  real(mytype) :: periodic_power_abs_error, periodic_power_rel_error
  real(mytype) :: deterministic_power_eulerian, deterministic_power_lagrangian
  real(mytype) :: deterministic_power_abs_error, deterministic_power_rel_error

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
  lag%weight = 0.03125_mytype

  allocate(ux(grid%nx, grid%ny, grid%nz), uy(grid%nx, grid%ny, grid%nz), uz(grid%nx, grid%ny, grid%nz))
  allocate(fx(grid%nx, grid%ny, grid%nz), fy(grid%nx, grid%ny, grid%nz), fz(grid%nx, grid%ny, grid%nz))
  allocate(u_lag(3, nlag))

  ! Test A: zero force power
  do k = 1, grid%nz
    z = grid%z(k)
    do j = 1, grid%ny
      y = grid%y(j)
      do i = 1, grid%nx
        x = grid%x(i)
        ux(i, j, k) = 1.0_mytype + 0.2_mytype * x - 0.1_mytype * y + 0.05_mytype * z
        uy(i, j, k) = -0.3_mytype + 0.4_mytype * y
        uz(i, j, k) = 0.7_mytype - 0.2_mytype * z + 0.1_mytype * x
      end do
    end do
  end do
  lag%force = 0._mytype
  call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
  call clear_eulerian_force_density(fx, fy, fz)
  call spread_lag_force_to_eulerian(grid, lag, fx, fy, fz)
  call compute_eulerian_power(grid, ux, uy, uz, fx, fy, fz, zero_force_power_eulerian)
  call compute_lagrangian_power(lag, u_lag, zero_force_power_lagrangian)
  call compute_power_consistency_error(zero_force_power_eulerian, zero_force_power_lagrangian, &
                                       zero_force_power_abs_error, unused_rel_error)

  ! Test B: constant velocity + multi-point force
  ux = 0.2_mytype
  uy = -0.1_mytype
  uz = 0.05_mytype
  do l = 1, nlag
    lag%force(1, l) = 0.1_mytype + 0.02_mytype * real(l, mytype)
    lag%force(2, l) = -0.05_mytype + 0.01_mytype * real(l, mytype)
    lag%force(3, l) = 0.03_mytype - 0.005_mytype * real(l, mytype)
  end do
  call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
  call clear_eulerian_force_density(fx, fy, fz)
  call spread_lag_force_to_eulerian(grid, lag, fx, fy, fz)
  call compute_eulerian_power(grid, ux, uy, uz, fx, fy, fz, constant_velocity_power_eulerian)
  call compute_lagrangian_power(lag, u_lag, constant_velocity_power_lagrangian)
  call compute_power_consistency_error(constant_velocity_power_eulerian, constant_velocity_power_lagrangian, &
                                       constant_velocity_power_abs_error, constant_velocity_power_rel_error)

  ! Test C: nonuniform velocity + multi-point force
  do k = 1, grid%nz
    z = grid%z(k)
    do j = 1, grid%ny
      y = grid%y(j)
      do i = 1, grid%nx
        x = grid%x(i)
        ux(i, j, k) = 1.0_mytype + 0.2_mytype * x - 0.1_mytype * y + 0.05_mytype * z
        uy(i, j, k) = -0.3_mytype + 0.4_mytype * y
        uz(i, j, k) = 0.7_mytype - 0.2_mytype * z + 0.1_mytype * x
      end do
    end do
  end do
  call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
  call clear_eulerian_force_density(fx, fy, fz)
  call spread_lag_force_to_eulerian(grid, lag, fx, fy, fz)
  call compute_eulerian_power(grid, ux, uy, uz, fx, fy, fz, nonuniform_power_eulerian)
  call compute_lagrangian_power(lag, u_lag, nonuniform_power_lagrangian)
  call compute_power_consistency_error(nonuniform_power_eulerian, nonuniform_power_lagrangian, &
                                       nonuniform_power_abs_error, nonuniform_power_rel_error)

  ! Test D: periodic-near-boundary lag point
  pi = acos(-1._mytype)
  lx = grid%xmax - grid%xmin
  lz = grid%zmax - grid%zmin
  do k = 1, grid%nz
    z = grid%z(k)
    do j = 1, grid%ny
      y = grid%y(j)
      do i = 1, grid%nx
        x = grid%x(i)
        ux(i, j, k) = sin(2._mytype * pi * (x - grid%xmin) / lx) + 0.1_mytype * cos(2._mytype * pi * (z - grid%zmin) / lz)
        uy(i, j, k) = 0.2_mytype * y
        uz(i, j, k) = cos(2._mytype * pi * (z - grid%zmin) / lz)
      end do
    end do
  end do
  lag%force = 0._mytype
  lag%force(:, 5) = [0.13_mytype, 0.07_mytype, -0.02_mytype]
  lag%weight = 0._mytype
  lag%weight(5) = 0.03125_mytype
  call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
  call clear_eulerian_force_density(fx, fy, fz)
  call spread_lag_force_to_eulerian(grid, lag, fx, fy, fz)
  call compute_eulerian_power(grid, ux, uy, uz, fx, fy, fz, periodic_power_eulerian)
  call compute_lagrangian_power(lag, u_lag, periodic_power_lagrangian)
  call compute_power_consistency_error(periodic_power_eulerian, periodic_power_lagrangian, &
                                       periodic_power_abs_error, periodic_power_rel_error)

  ! Test E: deterministic pseudo-random-like fields/forces
  do k = 1, grid%nz
    do j = 1, grid%ny
      do i = 1, grid%nx
        ux(i, j, k) = sin(0.17_mytype * real(i, mytype) + 0.31_mytype * real(j, mytype) + 0.23_mytype * real(k, mytype))
        uy(i, j, k) = cos(0.11_mytype * real(i, mytype) - 0.19_mytype * real(j, mytype) + 0.29_mytype * real(k, mytype))
        uz(i, j, k) = sin(0.07_mytype * real(i, mytype)) * cos(0.13_mytype * real(k, mytype)) + 0.05_mytype * real(j, mytype)
      end do
    end do
  end do
  do l = 1, nlag
    lag%force(1, l) = sin(0.4_mytype * real(l, mytype))
    lag%force(2, l) = cos(0.3_mytype * real(l, mytype))
    lag%force(3, l) = 0.2_mytype * sin(0.7_mytype * real(l, mytype))
  end do
  lag%weight = 0.03125_mytype
  call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
  call clear_eulerian_force_density(fx, fy, fz)
  call spread_lag_force_to_eulerian(grid, lag, fx, fy, fz)
  call compute_eulerian_power(grid, ux, uy, uz, fx, fy, fz, deterministic_power_eulerian)
  call compute_lagrangian_power(lag, u_lag, deterministic_power_lagrangian)
  call compute_power_consistency_error(deterministic_power_eulerian, deterministic_power_lagrangian, &
                                       deterministic_power_abs_error, deterministic_power_rel_error)

  power_consistency_status = 1
  if (ieee_is_nan(zero_force_power_eulerian) .or. ieee_is_nan(zero_force_power_lagrangian) .or. &
      ieee_is_nan(zero_force_power_abs_error) .or. ieee_is_nan(constant_velocity_power_eulerian) .or. &
      ieee_is_nan(constant_velocity_power_lagrangian) .or. ieee_is_nan(constant_velocity_power_abs_error) .or. &
      ieee_is_nan(constant_velocity_power_rel_error) .or. ieee_is_nan(nonuniform_power_eulerian) .or. &
      ieee_is_nan(nonuniform_power_lagrangian) .or. ieee_is_nan(nonuniform_power_abs_error) .or. &
      ieee_is_nan(nonuniform_power_rel_error) .or. ieee_is_nan(periodic_power_eulerian) .or. &
      ieee_is_nan(periodic_power_lagrangian) .or. ieee_is_nan(periodic_power_abs_error) .or. &
      ieee_is_nan(periodic_power_rel_error) .or. ieee_is_nan(deterministic_power_eulerian) .or. &
      ieee_is_nan(deterministic_power_lagrangian) .or. ieee_is_nan(deterministic_power_abs_error) .or. &
      ieee_is_nan(deterministic_power_rel_error)) then
    power_consistency_status = 0
  end if

  open(newunit=io_unit, file='stage3_outputs/fibre_ibm_power_consistency_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,ES24.16)') 'zero_force_power_eulerian = ', zero_force_power_eulerian
  write(io_unit, '(A,ES24.16)') 'zero_force_power_lagrangian = ', zero_force_power_lagrangian
  write(io_unit, '(A,ES24.16)') 'zero_force_power_abs_error = ', zero_force_power_abs_error

  write(io_unit, '(A,ES24.16)') 'constant_velocity_power_eulerian = ', constant_velocity_power_eulerian
  write(io_unit, '(A,ES24.16)') 'constant_velocity_power_lagrangian = ', constant_velocity_power_lagrangian
  write(io_unit, '(A,ES24.16)') 'constant_velocity_power_abs_error = ', constant_velocity_power_abs_error
  write(io_unit, '(A,ES24.16)') 'constant_velocity_power_rel_error = ', constant_velocity_power_rel_error

  write(io_unit, '(A,ES24.16)') 'nonuniform_power_eulerian = ', nonuniform_power_eulerian
  write(io_unit, '(A,ES24.16)') 'nonuniform_power_lagrangian = ', nonuniform_power_lagrangian
  write(io_unit, '(A,ES24.16)') 'nonuniform_power_abs_error = ', nonuniform_power_abs_error
  write(io_unit, '(A,ES24.16)') 'nonuniform_power_rel_error = ', nonuniform_power_rel_error

  write(io_unit, '(A,ES24.16)') 'periodic_power_eulerian = ', periodic_power_eulerian
  write(io_unit, '(A,ES24.16)') 'periodic_power_lagrangian = ', periodic_power_lagrangian
  write(io_unit, '(A,ES24.16)') 'periodic_power_abs_error = ', periodic_power_abs_error
  write(io_unit, '(A,ES24.16)') 'periodic_power_rel_error = ', periodic_power_rel_error

  write(io_unit, '(A,ES24.16)') 'deterministic_power_eulerian = ', deterministic_power_eulerian
  write(io_unit, '(A,ES24.16)') 'deterministic_power_lagrangian = ', deterministic_power_lagrangian
  write(io_unit, '(A,ES24.16)') 'deterministic_power_abs_error = ', deterministic_power_abs_error
  write(io_unit, '(A,ES24.16)') 'deterministic_power_rel_error = ', deterministic_power_rel_error

  write(io_unit, '(A,I0)') 'power_consistency_status = ', power_consistency_status

  close(io_unit)

  deallocate(ux, uy, uz, fx, fy, fz, u_lag)
  call destroy_lagrangian_points(lag)
  call destroy_ibm_grid(grid)

end program fibre_ibm_power_consistency_check
