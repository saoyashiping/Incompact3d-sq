program fibre_ibm_boundary_safety_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_uniform_ibm_grid, destroy_ibm_grid
  use fibre_ibm_grid, only : destroy_lagrangian_points
  use fibre_ibm_boundary_safety, only : ibm_point_safety_t
  use fibre_ibm_boundary_safety, only : init_point_safety, destroy_point_safety
  use fibre_ibm_boundary_safety, only : check_ibm_point_boundary_safety, summarize_ibm_boundary_safety
  use fibre_ibm_boundary_safety, only : assert_no_unsafe_ibm_points

  implicit none

  type(ibm_grid_t) :: grid
  type(ibm_lagrangian_points_t) :: lag
  type(ibm_point_safety_t) :: safety

  integer :: io_unit
  integer :: interior_safe_count, interior_periodic_wrap_count, interior_unsafe_count, interior_outside_count
  integer :: periodic_safe_count, periodic_wrap_count, periodic_unsafe_count, periodic_outside_count
  integer :: wall_safe_count, wall_periodic_wrap_count, wall_unsafe_count, wall_outside_count
  integer :: outside_safe_count, outside_periodic_wrap_count, outside_unsafe_count, outside_outside_count
  integer :: mixed_safe_count, mixed_periodic_wrap_count, mixed_unsafe_count, mixed_outside_count
  integer :: mixed_total_unsafe_count
  integer :: boundary_safety_policy_status, boundary_safety_status

  real(mytype) :: wall_min_wall_distance_min, wall_min_wall_distance_max

  call init_uniform_ibm_grid(grid, nx=16, ny=12, nz=10, &
                             xmin=0._mytype, xmax=2._mytype, &
                             ymin=-1._mytype, ymax=1._mytype, &
                             zmin=0._mytype, zmax=1._mytype, &
                             periodic_x=.true., periodic_y=.false., periodic_z=.true.)

  call init_point_safety(safety, 4)

  ! Test A: safe interior points
  call set_lag_points(lag, [0.73_mytype, 0.17_mytype, 0.41_mytype], &
                           [1.22_mytype, -0.25_mytype, 0.62_mytype], &
                           [1.00_mytype, 0.00_mytype, 0.50_mytype])
  call check_ibm_point_boundary_safety(grid, lag, safety)
  call summarize_ibm_boundary_safety(safety, interior_safe_count, interior_periodic_wrap_count, &
                                    interior_unsafe_count, interior_outside_count)

  ! Test B: periodic-near-boundary points
  call set_lag_points(lag, [0.05_mytype, 0.00_mytype, 0.50_mytype], &
                           [1.95_mytype, 0.00_mytype, 0.50_mytype], &
                           [1.00_mytype, 0.00_mytype, 0.95_mytype], &
                           [1.00_mytype, 0.00_mytype, 0.05_mytype])
  call check_ibm_point_boundary_safety(grid, lag, safety)
  call summarize_ibm_boundary_safety(safety, periodic_safe_count, periodic_wrap_count, &
                                    periodic_unsafe_count, periodic_outside_count)

  ! Test C: near y-wall unsafe points
  call set_lag_points(lag, [1.00_mytype, -0.95_mytype, 0.50_mytype], &
                           [1.00_mytype, 0.95_mytype, 0.50_mytype], &
                           [1.00_mytype, -0.80_mytype, 0.50_mytype], &
                           [1.00_mytype, 0.80_mytype, 0.50_mytype])
  call check_ibm_point_boundary_safety(grid, lag, safety)
  call summarize_ibm_boundary_safety(safety, wall_safe_count, wall_periodic_wrap_count, wall_unsafe_count, wall_outside_count)
  wall_min_wall_distance_min = minval(safety%min_wall_distance)
  wall_min_wall_distance_max = maxval(safety%min_wall_distance)

  ! Test D: outside nonperiodic y domain
  call set_lag_points(lag, [1.00_mytype, -1.10_mytype, 0.50_mytype], &
                           [1.00_mytype, 1.10_mytype, 0.50_mytype])
  call check_ibm_point_boundary_safety(grid, lag, safety)
  call summarize_ibm_boundary_safety(safety, outside_safe_count, outside_periodic_wrap_count, &
                                    outside_unsafe_count, outside_outside_count)

  ! Test E: mixed-point classification
  call set_lag_points(lag, [1.00_mytype, 0.00_mytype, 0.50_mytype], &
                           [0.05_mytype, 0.00_mytype, 0.95_mytype], &
                           [1.00_mytype, 0.95_mytype, 0.50_mytype], &
                           [1.00_mytype, 1.10_mytype, 0.50_mytype])
  call check_ibm_point_boundary_safety(grid, lag, safety)
  call summarize_ibm_boundary_safety(safety, mixed_safe_count, mixed_periodic_wrap_count, &
                                    mixed_unsafe_count, mixed_outside_count)
  call assert_no_unsafe_ibm_points(safety, mixed_total_unsafe_count)

  ! Stage 3.8 only detects unsafe points. It does not apply wall-corrected kernels or silently
  ! renormalize truncated kernels. Later coupling stages must check unsafe_count before
  ! interpolation/spreading near walls.
  boundary_safety_policy_status = 1

  boundary_safety_status = 1
  if (ieee_is_nan(wall_min_wall_distance_min) .or. ieee_is_nan(wall_min_wall_distance_max)) then
    boundary_safety_status = 0
  end if

  open(newunit=io_unit, file='stage3_outputs/fibre_ibm_boundary_safety_check.dat', &
       status='replace', action='write', form='formatted')

  write(io_unit, '(A,I0)') 'interior_safe_count = ', interior_safe_count
  write(io_unit, '(A,I0)') 'interior_periodic_wrap_count = ', interior_periodic_wrap_count
  write(io_unit, '(A,I0)') 'interior_unsafe_count = ', interior_unsafe_count
  write(io_unit, '(A,I0)') 'interior_outside_count = ', interior_outside_count

  write(io_unit, '(A,I0)') 'periodic_safe_count = ', periodic_safe_count
  write(io_unit, '(A,I0)') 'periodic_wrap_count = ', periodic_wrap_count
  write(io_unit, '(A,I0)') 'periodic_unsafe_count = ', periodic_unsafe_count
  write(io_unit, '(A,I0)') 'periodic_outside_count = ', periodic_outside_count

  write(io_unit, '(A,I0)') 'wall_safe_count = ', wall_safe_count
  write(io_unit, '(A,I0)') 'wall_periodic_wrap_count = ', wall_periodic_wrap_count
  write(io_unit, '(A,I0)') 'wall_unsafe_count = ', wall_unsafe_count
  write(io_unit, '(A,I0)') 'wall_outside_count = ', wall_outside_count
  write(io_unit, '(A,ES24.16)') 'wall_min_wall_distance_min = ', wall_min_wall_distance_min
  write(io_unit, '(A,ES24.16)') 'wall_min_wall_distance_max = ', wall_min_wall_distance_max

  write(io_unit, '(A,I0)') 'outside_safe_count = ', outside_safe_count
  write(io_unit, '(A,I0)') 'outside_periodic_wrap_count = ', outside_periodic_wrap_count
  write(io_unit, '(A,I0)') 'outside_unsafe_count = ', outside_unsafe_count
  write(io_unit, '(A,I0)') 'outside_outside_count = ', outside_outside_count

  write(io_unit, '(A,I0)') 'mixed_safe_count = ', mixed_safe_count
  write(io_unit, '(A,I0)') 'mixed_periodic_wrap_count = ', mixed_periodic_wrap_count
  write(io_unit, '(A,I0)') 'mixed_unsafe_count = ', mixed_unsafe_count
  write(io_unit, '(A,I0)') 'mixed_outside_count = ', mixed_outside_count
  write(io_unit, '(A,I0)') 'mixed_total_unsafe_count = ', mixed_total_unsafe_count

  write(io_unit, '(A,I0)') 'boundary_safety_policy_status = ', boundary_safety_policy_status
  write(io_unit, '(A,I0)') 'boundary_safety_status = ', boundary_safety_status

  close(io_unit)

  call destroy_point_safety(safety)
  call destroy_lagrangian_points(lag)
  call destroy_ibm_grid(grid)

contains

  subroutine set_lag_points(lag_in, p1, p2, p3, p4)
    type(ibm_lagrangian_points_t), intent(inout) :: lag_in
    real(mytype), intent(in) :: p1(3)
    real(mytype), intent(in), optional :: p2(3), p3(3), p4(3)

    integer :: n

    n = 1
    if (present(p2)) n = n + 1
    if (present(p3)) n = n + 1
    if (present(p4)) n = n + 1

    call destroy_lagrangian_points(lag_in)
    lag_in%nl = n
    allocate(lag_in%x(3, n), lag_in%v(3, n), lag_in%force(3, n), lag_in%weight(n))

    lag_in%x(:, 1) = p1
    if (present(p2)) lag_in%x(:, 2) = p2
    if (present(p3)) lag_in%x(:, 3) = p3
    if (present(p4)) lag_in%x(:, 4) = p4

    lag_in%v = 0._mytype
    lag_in%force = 0._mytype
    lag_in%weight = 0.03125_mytype
  end subroutine set_lag_points

end program fibre_ibm_boundary_safety_check
