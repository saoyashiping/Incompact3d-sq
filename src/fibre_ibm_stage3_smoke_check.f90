program fibre_ibm_stage3_smoke_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use fibre_parameters, only : mytype
  use fibre_types, only : fibre_t, fibre_allocate
  use fibre_geometry, only : fibre_end_to_end_distance
  use fibre_external_force, only : set_fibre_external_force
  use fibre_structure_integrator, only : advance_fibre_structure_freefree
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_uniform_ibm_grid, destroy_ibm_grid
  use fibre_ibm_grid, only : destroy_lagrangian_points
  use fibre_ibm_interpolation, only : interpolate_vector_to_lag
  use fibre_ibm_feedback, only : compute_ibm_feedback_forces
  use fibre_ibm_power_diagnostics, only : compute_eulerian_power, compute_lagrangian_power
  use fibre_ibm_force_buffer, only : ibm_force_buffer_t, allocate_ibm_force_buffer
  use fibre_ibm_force_buffer, only : destroy_ibm_force_buffer, clear_ibm_force_buffer
  use fibre_ibm_force_buffer, only : accumulate_lag_force_to_buffer, compute_ibm_force_buffer_total_force
  use fibre_ibm_force_buffer, only : compute_ibm_force_buffer_norms
  use fibre_ibm_boundary_safety, only : ibm_point_safety_t, init_point_safety, destroy_point_safety
  use fibre_ibm_boundary_safety, only : check_ibm_point_boundary_safety, summarize_ibm_boundary_safety
  use fibre_ibm_boundary_safety, only : assert_no_unsafe_ibm_points

  implicit none

  type(ibm_grid_t) :: grid
  type(fibre_t) :: fibre
  type(ibm_lagrangian_points_t) :: lag
  type(ibm_force_buffer_t) :: buffer
  type(ibm_point_safety_t) :: safety

  integer, parameter :: nl = 33
  integer, parameter :: nsteps = 20
  integer, parameter :: nsteps_nonuniform = 5
  integer :: i, j, k, l, step, io_unit, ierr
  integer :: safe_count, periodic_wrap_count, unsafe_near_wall_count, outside_count, unsafe_count
  integer :: uniform_smoke_direction_check, uniform_smoke_solver_failure_count, uniform_smoke_nan_detected
  integer :: nonuniform_smoke_nan_detected, nonuniform_smoke_unsafe_count
  integer :: stage3_smoke_status

  real(mytype), parameter :: fibre_length = 1.0_mytype
  real(mytype), parameter :: rho_tilde = 1.0_mytype
  real(mytype), parameter :: gamma = 1.0_mytype
  real(mytype), parameter :: beta_drag = 10.0_mytype
  real(mytype), parameter :: dt = 1.0e-5_mytype

  real(mytype), allocatable :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
  real(mytype), allocatable :: u_lag(:,:), force_on_structure(:,:), force_on_fluid(:,:)
  real(mytype), allocatable :: x_initial(:,:)

  real(mytype) :: lag_total_force(3), buffer_total_force(3)
  real(mytype) :: x, y, z
  real(mytype) :: tension_residual, relative_tension_residual

  real(mytype) :: smoke_force_conservation_error, smoke_force_conservation_relative_error, smoke_force_buffer_max_abs
  real(mytype) :: smoke_power_eulerian, smoke_power_lagrangian, smoke_power_abs_error, smoke_power_relative_error
  real(mytype) :: smoke_power_recomputed_abs_error, smoke_power_error_consistency_check

  real(mytype) :: smoke_boundary_safe_count, smoke_boundary_periodic_wrap_count
  real(mytype) :: smoke_boundary_unsafe_count, smoke_boundary_outside_count, smoke_boundary_total_unsafe_count

  real(mytype) :: uniform_smoke_final_center_velocity_x, uniform_smoke_length_error, uniform_smoke_shape_error_max
  real(mytype) :: zero_slip_smoke_f_ext_norm, zero_slip_smoke_force_buffer_norm

  real(mytype) :: nonuniform_smoke_f_ext_norm, nonuniform_smoke_force_buffer_norm
  real(mytype) :: nonuniform_smoke_force_conservation_relative_error, nonuniform_smoke_power_relative_error
  real(mytype) :: nonuniform_smoke_power_eulerian, nonuniform_smoke_power_lagrangian
  real(mytype) :: nonuniform_smoke_power_recomputed_abs_error, nonuniform_smoke_power_error_consistency_check
  real(mytype) :: nonuniform_smoke_length_error
  real(mytype) :: tmp_buffer_l2_uniform, tmp_buffer_l2_nonuniform

  call init_uniform_ibm_grid(grid, nx=16, ny=12, nz=10, &
                             xmin=0._mytype, xmax=2._mytype, &
                             ymin=-1._mytype, ymax=1._mytype, &
                             zmin=0._mytype, zmax=1._mytype, &
                             periodic_x=.true., periodic_y=.false., periodic_z=.true.)

  call fibre_allocate(fibre, nl)
  fibre%length = fibre_length
  fibre%ds = fibre%length / real(fibre%nl - 1, mytype)
  fibre%rho_tilde = rho_tilde
  fibre%gamma = gamma

  allocate(ux(grid%nx, grid%ny, grid%nz), uy(grid%nx, grid%ny, grid%nz), uz(grid%nx, grid%ny, grid%nz))
  allocate(u_lag(3, fibre%nl), force_on_structure(3, fibre%nl), force_on_fluid(3, fibre%nl))
  allocate(x_initial(3, fibre%nl))

  call allocate_ibm_force_buffer(buffer, grid)
  call init_point_safety(safety, fibre%nl)

  ! Test A: boundary safety pre-check
  call reset_straight_fibre(fibre)
  call sync_lag_from_fibre(fibre, lag)
  call check_ibm_point_boundary_safety(grid, lag, safety)
  call summarize_ibm_boundary_safety(safety, safe_count, periodic_wrap_count, unsafe_near_wall_count, outside_count)
  call assert_no_unsafe_ibm_points(safety, unsafe_count)
  smoke_boundary_safe_count = real(safe_count, mytype)
  smoke_boundary_periodic_wrap_count = real(periodic_wrap_count, mytype)
  smoke_boundary_unsafe_count = real(unsafe_near_wall_count, mytype)
  smoke_boundary_outside_count = real(outside_count, mytype)
  smoke_boundary_total_unsafe_count = real(unsafe_count, mytype)

  ! Test B/C/D: uniform-flow one-way smoke chain
  call reset_straight_fibre(fibre)
  x_initial = fibre%x
  ux = 0.2_mytype
  uy = 0._mytype
  uz = 0._mytype
  uniform_smoke_solver_failure_count = 0
  uniform_smoke_nan_detected = 0

  do step = 1, nsteps
    call sync_lag_from_fibre(fibre, lag)
    call check_ibm_point_boundary_safety(grid, lag, safety)
    call assert_no_unsafe_ibm_points(safety, unsafe_count)
    if (unsafe_count > 0) then
      uniform_smoke_solver_failure_count = uniform_smoke_solver_failure_count + 1
      exit
    end if

    call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
    call compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)
    call set_fibre_external_force(fibre, force_on_structure)
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
    if (ierr /= 0) uniform_smoke_solver_failure_count = uniform_smoke_solver_failure_count + 1

    call clear_ibm_force_buffer(buffer)
    lag%force = force_on_fluid
    call accumulate_lag_force_to_buffer(grid, lag, buffer)

    lag_total_force = 0._mytype
    do l = 1, lag%nl
      lag_total_force = lag_total_force + force_on_fluid(:, l) * lag%weight(l)
    end do
    call compute_ibm_force_buffer_total_force(grid, buffer, buffer_total_force)
    smoke_force_conservation_error = sqrt(sum((buffer_total_force - lag_total_force)**2))
    smoke_force_conservation_relative_error = smoke_force_conservation_error / max(sqrt(sum(lag_total_force**2)), 1.0e-30_mytype)

    call compute_eulerian_power(grid, ux, uy, uz, buffer%fx, buffer%fy, buffer%fz, smoke_power_eulerian)
    call compute_lagrangian_power(lag, u_lag, smoke_power_lagrangian)
    smoke_power_abs_error = abs(smoke_power_eulerian - smoke_power_lagrangian)
    smoke_power_relative_error = smoke_power_abs_error / max(abs(smoke_power_lagrangian), 1.0e-300_mytype)
    smoke_power_recomputed_abs_error = abs(smoke_power_eulerian - smoke_power_lagrangian)
    smoke_power_error_consistency_check = abs(smoke_power_recomputed_abs_error - smoke_power_abs_error)

    if (ieee_is_nan(smoke_force_conservation_error) .or. ieee_is_nan(smoke_power_abs_error)) then
      uniform_smoke_nan_detected = 1
    end if
  end do

  uniform_smoke_final_center_velocity_x = sum(fibre%v(1, :)) / real(fibre%nl, mytype)
  uniform_smoke_direction_check = 0
  if (uniform_smoke_final_center_velocity_x > 0._mytype) uniform_smoke_direction_check = 1
  uniform_smoke_length_error = abs(fibre_end_to_end_distance(fibre) - fibre%length)
  uniform_smoke_shape_error_max = maxval(abs(fibre%x(2:3, :) - x_initial(2:3, :)))
  call compute_ibm_force_buffer_norms(buffer, smoke_force_buffer_max_abs, tmp_buffer_l2_uniform)

  ! Test E: zero-slip smoke
  call reset_straight_fibre(fibre)
  ux = 0.2_mytype
  uy = -0.1_mytype
  uz = 0.05_mytype
  do l = 1, fibre%nl
    fibre%v(:, l) = [0.2_mytype, -0.1_mytype, 0.05_mytype]
    fibre%x_old(:, l) = fibre%x(:, l) - fibre%v(:, l) * dt
  end do
  call sync_lag_from_fibre(fibre, lag)
  call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
  call compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)
  call set_fibre_external_force(fibre, force_on_structure)
  zero_slip_smoke_f_ext_norm = sqrt(sum(fibre%f_ext**2))
  call clear_ibm_force_buffer(buffer)
  lag%force = force_on_fluid
  call accumulate_lag_force_to_buffer(grid, lag, buffer)
  zero_slip_smoke_force_buffer_norm = sqrt(sum(buffer%fx**2 + buffer%fy**2 + buffer%fz**2))

  ! Test F: nonuniform-flow short smoke
  call reset_straight_fibre(fibre)
  x_initial = fibre%x
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

  nonuniform_smoke_nan_detected = 0
  nonuniform_smoke_unsafe_count = 0
  do step = 1, nsteps_nonuniform
    call sync_lag_from_fibre(fibre, lag)
    call check_ibm_point_boundary_safety(grid, lag, safety)
    call assert_no_unsafe_ibm_points(safety, unsafe_count)
    nonuniform_smoke_unsafe_count = max(nonuniform_smoke_unsafe_count, unsafe_count)

    call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
    call compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)
    call set_fibre_external_force(fibre, force_on_structure)
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)

    call clear_ibm_force_buffer(buffer)
    lag%force = force_on_fluid
    call accumulate_lag_force_to_buffer(grid, lag, buffer)

    lag_total_force = 0._mytype
    do l = 1, lag%nl
      lag_total_force = lag_total_force + force_on_fluid(:, l) * lag%weight(l)
    end do
    call compute_ibm_force_buffer_total_force(grid, buffer, buffer_total_force)
    nonuniform_smoke_force_conservation_relative_error = &
      sqrt(sum((buffer_total_force - lag_total_force)**2)) / max(sqrt(sum(lag_total_force**2)), 1.0e-30_mytype)

    call compute_eulerian_power(grid, ux, uy, uz, buffer%fx, buffer%fy, buffer%fz, nonuniform_smoke_power_eulerian)
    call compute_lagrangian_power(lag, u_lag, nonuniform_smoke_power_lagrangian)
    nonuniform_smoke_power_recomputed_abs_error = abs(nonuniform_smoke_power_eulerian - nonuniform_smoke_power_lagrangian)
    nonuniform_smoke_power_relative_error = nonuniform_smoke_power_recomputed_abs_error / &
      max(abs(nonuniform_smoke_power_lagrangian), 1.0e-300_mytype)
    nonuniform_smoke_power_error_consistency_check = abs(nonuniform_smoke_power_recomputed_abs_error - &
      abs(nonuniform_smoke_power_eulerian - nonuniform_smoke_power_lagrangian))

    if (ieee_is_nan(nonuniform_smoke_force_conservation_relative_error) .or. &
        ieee_is_nan(nonuniform_smoke_power_relative_error) .or. &
        ieee_is_nan(nonuniform_smoke_power_recomputed_abs_error) .or. &
        ieee_is_nan(nonuniform_smoke_power_error_consistency_check)) then
      nonuniform_smoke_nan_detected = 1
    end if
  end do

  nonuniform_smoke_f_ext_norm = sqrt(sum(fibre%f_ext**2))
  call compute_ibm_force_buffer_norms(buffer, nonuniform_smoke_force_buffer_norm, tmp_buffer_l2_nonuniform)
  nonuniform_smoke_length_error = abs(fibre_end_to_end_distance(fibre) - fibre%length)

  stage3_smoke_status = 1
  if (uniform_smoke_nan_detected /= 0 .or. nonuniform_smoke_nan_detected /= 0) stage3_smoke_status = 0

  open(newunit=io_unit, file='stage3_outputs/fibre_ibm_stage3_smoke_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,I0)') 'smoke_boundary_safe_count = ', int(smoke_boundary_safe_count)
  write(io_unit, '(A,I0)') 'smoke_boundary_periodic_wrap_count = ', int(smoke_boundary_periodic_wrap_count)
  write(io_unit, '(A,I0)') 'smoke_boundary_unsafe_count = ', int(smoke_boundary_unsafe_count)
  write(io_unit, '(A,I0)') 'smoke_boundary_outside_count = ', int(smoke_boundary_outside_count)
  write(io_unit, '(A,I0)') 'smoke_boundary_total_unsafe_count = ', int(smoke_boundary_total_unsafe_count)

  write(io_unit, '(A,ES24.16)') 'uniform_smoke_final_center_velocity_x = ', uniform_smoke_final_center_velocity_x
  write(io_unit, '(A,I0)') 'uniform_smoke_direction_check = ', uniform_smoke_direction_check
  write(io_unit, '(A,ES24.16)') 'uniform_smoke_length_error = ', uniform_smoke_length_error
  write(io_unit, '(A,ES24.16)') 'uniform_smoke_shape_error_max = ', uniform_smoke_shape_error_max
  write(io_unit, '(A,I0)') 'uniform_smoke_solver_failure_count = ', uniform_smoke_solver_failure_count
  write(io_unit, '(A,I0)') 'uniform_smoke_nan_detected = ', uniform_smoke_nan_detected

  write(io_unit, '(A,ES24.16)') 'smoke_force_conservation_error = ', smoke_force_conservation_error
  write(io_unit, '(A,ES24.16)') 'smoke_force_conservation_relative_error = ', smoke_force_conservation_relative_error
  write(io_unit, '(A,ES24.16)') 'smoke_force_buffer_max_abs = ', smoke_force_buffer_max_abs

  write(io_unit, '(A,ES24.16)') 'smoke_power_eulerian = ', smoke_power_eulerian
  write(io_unit, '(A,ES24.16)') 'smoke_power_lagrangian = ', smoke_power_lagrangian
  write(io_unit, '(A,ES24.16)') 'smoke_power_abs_error = ', smoke_power_abs_error
  write(io_unit, '(A,ES24.16)') 'smoke_power_relative_error = ', smoke_power_relative_error
  write(io_unit, '(A,ES24.16)') 'smoke_power_recomputed_abs_error = ', smoke_power_recomputed_abs_error
  write(io_unit, '(A,ES24.16)') 'smoke_power_error_consistency_check = ', smoke_power_error_consistency_check

  write(io_unit, '(A,ES24.16)') 'zero_slip_smoke_f_ext_norm = ', zero_slip_smoke_f_ext_norm
  write(io_unit, '(A,ES24.16)') 'zero_slip_smoke_force_buffer_norm = ', zero_slip_smoke_force_buffer_norm

  write(io_unit, '(A,ES24.16)') 'nonuniform_smoke_f_ext_norm = ', nonuniform_smoke_f_ext_norm
  write(io_unit, '(A,ES24.16)') 'nonuniform_smoke_force_buffer_norm = ', nonuniform_smoke_force_buffer_norm
  write(io_unit, '(A,ES24.16)') 'nonuniform_smoke_force_conservation_relative_error = ', &
                             nonuniform_smoke_force_conservation_relative_error
  write(io_unit, '(A,ES24.16)') 'nonuniform_smoke_power_relative_error = ', nonuniform_smoke_power_relative_error
  write(io_unit, '(A,ES24.16)') 'nonuniform_smoke_power_recomputed_abs_error = ', nonuniform_smoke_power_recomputed_abs_error
  write(io_unit, '(A,ES24.16)') 'nonuniform_smoke_power_error_consistency_check = ', nonuniform_smoke_power_error_consistency_check
  write(io_unit, '(A,ES24.16)') 'nonuniform_smoke_length_error = ', nonuniform_smoke_length_error
  write(io_unit, '(A,I0)') 'nonuniform_smoke_nan_detected = ', nonuniform_smoke_nan_detected
  write(io_unit, '(A,I0)') 'nonuniform_smoke_unsafe_count = ', nonuniform_smoke_unsafe_count

  write(io_unit, '(A,I0)') 'stage3_smoke_status = ', stage3_smoke_status

  close(io_unit)

  call destroy_point_safety(safety)
  call destroy_ibm_force_buffer(buffer)
  deallocate(ux, uy, uz, u_lag, force_on_structure, force_on_fluid, x_initial)
  call destroy_lagrangian_points(lag)
  call destroy_ibm_grid(grid)

contains

  subroutine reset_straight_fibre(f)
    type(fibre_t), intent(inout) :: f
    integer :: p
    real(mytype) :: s

    do p = 1, f%nl
      s = real(p - 1, mytype) / real(f%nl - 1, mytype)
      f%x(1, p) = 0.5_mytype + s
      f%x(2, p) = 0._mytype
      f%x(3, p) = 0.5_mytype
    end do
    f%x_old = f%x
    f%v = 0._mytype
    f%f_ext = 0._mytype
    f%tension = 0._mytype
  end subroutine reset_straight_fibre

  subroutine sync_lag_from_fibre(f, lgr)
    type(fibre_t), intent(in) :: f
    type(ibm_lagrangian_points_t), intent(inout) :: lgr

    call destroy_lagrangian_points(lgr)
    lgr%nl = f%nl
    allocate(lgr%x(3, lgr%nl), lgr%v(3, lgr%nl), lgr%force(3, lgr%nl), lgr%weight(lgr%nl))
    lgr%x = f%x
    lgr%v = f%v
    lgr%force = 0._mytype
    lgr%weight = f%ds
    if (lgr%nl >= 1) lgr%weight(1) = 0.5_mytype * f%ds
    if (lgr%nl >= 2) lgr%weight(lgr%nl) = 0.5_mytype * f%ds
  end subroutine sync_lag_from_fibre

end program fibre_ibm_stage3_smoke_check
