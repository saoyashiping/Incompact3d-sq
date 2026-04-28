program fibre_ibm_interpolation_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_uniform_ibm_grid, destroy_ibm_grid
  use fibre_ibm_grid, only : destroy_lagrangian_points
  use fibre_ibm_interpolation, only : interpolate_scalar_to_lag, interpolate_vector_to_lag
  use fibre_ibm_interpolation, only : compute_interpolation_weight_sum

  implicit none

  type(ibm_grid_t) :: grid
  type(ibm_lagrangian_points_t) :: lag

  integer :: i, j, k, l, io_unit
  integer, parameter :: nlag = 5
  integer, parameter :: n_inner = 3
  integer :: interpolation_status

  real(mytype), allocatable :: scalar(:,:,:), ux(:,:,:), uy(:,:,:), uz(:,:,:)
  real(mytype), allocatable :: q_lag(:), wsum(:), q_exact(:)
  real(mytype), allocatable :: q_periodic_grid(:,:,:), q_periodic_lag(:), q_periodic_exact(:)
  real(mytype), allocatable :: u_lag(:,:), u_exact(:,:)

  real(mytype) :: q0
  real(mytype) :: a0, ax, ay, az
  real(mytype) :: scalar_constant_max_error, scalar_constant_weight_sum_max_error
  real(mytype) :: vector_constant_max_error
  real(mytype) :: scalar_linear_inner_max_error, vector_linear_inner_max_error
  real(mytype) :: periodic_wrap_constant_error
  real(mytype) :: periodic_wrap_periodic_field_error_p1, periodic_wrap_periodic_field_error_p2
  real(mytype) :: interpolation_weight_sum_min, interpolation_weight_sum_max
  real(mytype) :: lx, lz, pi

  call init_uniform_ibm_grid(grid, nx=16, ny=12, nz=10, &
                             xmin=0._mytype, xmax=2._mytype, &
                             ymin=-1._mytype, ymax=1._mytype, &
                             zmin=0._mytype, zmax=1._mytype, &
                             periodic_x=.true., periodic_y=.false., periodic_z=.true.)

  allocate(lag%x(3, nlag), lag%v(3, nlag), lag%force(3, nlag), lag%weight(nlag))
  lag%nl = nlag
  lag%x(:, 1) = [0.73_mytype,  0.17_mytype, 0.41_mytype]
  lag%x(:, 2) = [1.22_mytype, -0.25_mytype, 0.62_mytype]
  lag%x(:, 3) = [0.15_mytype,  0.33_mytype, 0.80_mytype]
  lag%x(:, 4) = [1.85_mytype, -0.10_mytype, 0.12_mytype]
  lag%x(:, 5) = [0.05_mytype,  0.00_mytype, 0.95_mytype]
  lag%v = 0._mytype
  lag%force = 0._mytype
  lag%weight = 1._mytype

  allocate(scalar(grid%nx, grid%ny, grid%nz), ux(grid%nx, grid%ny, grid%nz), uy(grid%nx, grid%ny, grid%nz), uz(grid%nx, grid%ny, grid%nz))
  allocate(q_lag(nlag), wsum(nlag), q_exact(nlag), u_lag(3, nlag), u_exact(3, nlag))
  allocate(q_periodic_grid(grid%nx, grid%ny, grid%nz), q_periodic_lag(nlag), q_periodic_exact(nlag))

  q0 = 3.141592653589793_mytype

  scalar = q0
  call interpolate_scalar_to_lag(grid, scalar, lag, q_lag)
  call compute_interpolation_weight_sum(grid, lag, wsum)

  scalar_constant_max_error = maxval(abs(q_lag - q0))
  scalar_constant_weight_sum_max_error = maxval(abs(wsum - 1._mytype))

  ux = 0.2_mytype
  uy = -0.1_mytype
  uz = 0.05_mytype
  call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)

  vector_constant_max_error = 0._mytype
  do l = 1, nlag
    vector_constant_max_error = max(vector_constant_max_error, abs(u_lag(1, l) - 0.2_mytype))
    vector_constant_max_error = max(vector_constant_max_error, abs(u_lag(2, l) + 0.1_mytype))
    vector_constant_max_error = max(vector_constant_max_error, abs(u_lag(3, l) - 0.05_mytype))
  end do

  a0 = 1.25_mytype
  ax = 0.7_mytype
  ay = -0.4_mytype
  az = 0.3_mytype

  do k = 1, grid%nz
    do j = 1, grid%ny
      do i = 1, grid%nx
        scalar(i, j, k) = a0 + ax * grid%x(i) + ay * grid%y(j) + az * grid%z(k)
      end do
    end do
  end do
  call interpolate_scalar_to_lag(grid, scalar, lag, q_lag)

  scalar_linear_inner_max_error = 0._mytype
  do l = 1, n_inner
    q_exact(l) = a0 + ax * lag%x(1, l) + ay * lag%x(2, l) + az * lag%x(3, l)
    scalar_linear_inner_max_error = max(scalar_linear_inner_max_error, abs(q_lag(l) - q_exact(l)))
  end do

  do k = 1, grid%nz
    do j = 1, grid%ny
      do i = 1, grid%nx
        ux(i, j, k) = 1.0_mytype + 0.2_mytype * grid%x(i) - 0.1_mytype * grid%y(j) + 0.05_mytype * grid%z(k)
        uy(i, j, k) = -0.3_mytype + 0.4_mytype * grid%y(j)
        uz(i, j, k) = 0.7_mytype - 0.2_mytype * grid%z(k) + 0.1_mytype * grid%x(i)
      end do
    end do
  end do
  call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)

  vector_linear_inner_max_error = 0._mytype
  do l = 1, n_inner
    u_exact(1, l) = 1.0_mytype + 0.2_mytype * lag%x(1, l) - 0.1_mytype * lag%x(2, l) + 0.05_mytype * lag%x(3, l)
    u_exact(2, l) = -0.3_mytype + 0.4_mytype * lag%x(2, l)
    u_exact(3, l) = 0.7_mytype - 0.2_mytype * lag%x(3, l) + 0.1_mytype * lag%x(1, l)
    vector_linear_inner_max_error = max(vector_linear_inner_max_error, abs(u_lag(1, l) - u_exact(1, l)))
    vector_linear_inner_max_error = max(vector_linear_inner_max_error, abs(u_lag(2, l) - u_exact(2, l)))
    vector_linear_inner_max_error = max(vector_linear_inner_max_error, abs(u_lag(3, l) - u_exact(3, l)))
  end do

  scalar = q0
  call interpolate_scalar_to_lag(grid, scalar, lag, q_lag)
  periodic_wrap_constant_error = abs(q_lag(5) - q0)

  lx = grid%xmax - grid%xmin
  lz = grid%zmax - grid%zmin
  pi = acos(-1.0_mytype)
  do k = 1, grid%nz
    do j = 1, grid%ny
      do i = 1, grid%nx
        q_periodic_grid(i, j, k) = sin(2._mytype * pi * grid%x(i) / lx) + cos(2._mytype * pi * grid%z(k) / lz) + 0.2_mytype * grid%y(j)
      end do
    end do
  end do

  call interpolate_scalar_to_lag(grid, q_periodic_grid, lag, q_periodic_lag)
  do l = 1, nlag
    q_periodic_exact(l) = sin(2._mytype * pi * lag%x(1, l) / lx) + cos(2._mytype * pi * lag%x(3, l) / lz) + 0.2_mytype * lag%x(2, l)
  end do
  periodic_wrap_periodic_field_error_p1 = abs(q_periodic_lag(5) - q_periodic_exact(5))
  periodic_wrap_periodic_field_error_p2 = abs(q_periodic_lag(4) - q_periodic_exact(4))

  interpolation_weight_sum_min = minval(wsum)
  interpolation_weight_sum_max = maxval(wsum)

  interpolation_status = 1
  if (ieee_is_nan(scalar_constant_max_error) .or. ieee_is_nan(scalar_constant_weight_sum_max_error) .or. &
      ieee_is_nan(vector_constant_max_error) .or. ieee_is_nan(scalar_linear_inner_max_error) .or. &
      ieee_is_nan(vector_linear_inner_max_error) .or. ieee_is_nan(periodic_wrap_constant_error) .or. &
      ieee_is_nan(periodic_wrap_periodic_field_error_p1) .or. ieee_is_nan(periodic_wrap_periodic_field_error_p2) .or. &
      ieee_is_nan(interpolation_weight_sum_min) .or. ieee_is_nan(interpolation_weight_sum_max)) then
    interpolation_status = 0
  end if

  open(newunit=io_unit, file='stage3_outputs/fibre_ibm_interpolation_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,ES24.16)') 'scalar_constant_max_error = ', scalar_constant_max_error
  write(io_unit, '(A,ES24.16)') 'scalar_constant_weight_sum_max_error = ', scalar_constant_weight_sum_max_error
  write(io_unit, '(A,ES24.16)') 'vector_constant_max_error = ', vector_constant_max_error
  write(io_unit, '(A,ES24.16)') 'scalar_linear_inner_max_error = ', scalar_linear_inner_max_error
  write(io_unit, '(A,ES24.16)') 'vector_linear_inner_max_error = ', vector_linear_inner_max_error
  write(io_unit, '(A,ES24.16)') 'periodic_wrap_constant_error = ', periodic_wrap_constant_error
  write(io_unit, '(A,ES24.16)') 'periodic_wrap_periodic_field_error_p1 = ', periodic_wrap_periodic_field_error_p1
  write(io_unit, '(A,ES24.16)') 'periodic_wrap_periodic_field_error_p2 = ', periodic_wrap_periodic_field_error_p2
  write(io_unit, '(A,ES24.16)') 'interpolation_weight_sum_min = ', interpolation_weight_sum_min
  write(io_unit, '(A,ES24.16)') 'interpolation_weight_sum_max = ', interpolation_weight_sum_max
  write(io_unit, '(A,I0)') 'interpolation_status = ', interpolation_status

  close(io_unit)

  deallocate(scalar, ux, uy, uz, q_lag, wsum, q_exact, q_periodic_grid, q_periodic_lag, q_periodic_exact, u_lag, u_exact)
  call destroy_lagrangian_points(lag)
  call destroy_ibm_grid(grid)

end program fibre_ibm_interpolation_check
