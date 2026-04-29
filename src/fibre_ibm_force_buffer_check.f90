program fibre_ibm_force_buffer_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_uniform_ibm_grid, destroy_ibm_grid
  use fibre_ibm_grid, only : destroy_lagrangian_points
  use fibre_ibm_interpolation, only : interpolate_vector_to_lag
  use fibre_ibm_feedback, only : compute_ibm_feedback_forces
  use fibre_ibm_force_buffer, only : ibm_force_buffer_t
  use fibre_ibm_force_buffer, only : allocate_ibm_force_buffer, destroy_ibm_force_buffer
  use fibre_ibm_force_buffer, only : clear_ibm_force_buffer, accumulate_lag_force_to_buffer
  use fibre_ibm_force_buffer, only : compute_ibm_force_buffer_total_force, compute_ibm_force_buffer_norms

  implicit none

  type(ibm_grid_t) :: grid
  type(ibm_lagrangian_points_t) :: lag
  type(ibm_force_buffer_t) :: buffer, buffer_a, buffer_b, buffer_ab

  integer, parameter :: nlag = 5
  integer :: i, j, k, l, io_unit
  integer :: buffer_allocated_flag, force_buffer_status

  real(mytype), parameter :: beta_drag = 10.0_mytype

  real(mytype), allocatable :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
  real(mytype), allocatable :: u_lag(:,:), force_on_structure(:,:), force_on_fluid(:,:)
  real(mytype), allocatable :: force_a(:,:), force_b(:,:)

  real(mytype) :: lag_total_force(3), buffer_total_force(3)
  real(mytype) :: total_a(3), total_b(3), total_ab(3)
  real(mytype) :: ref_norm, x, y, z

  real(mytype) :: buffer_initial_max_abs, buffer_initial_l2_norm
  real(mytype) :: single_accumulate_total_force_error, single_accumulate_relative_force_error
  real(mytype) :: single_accumulate_buffer_max_abs, single_accumulate_buffer_l2_norm
  real(mytype) :: clear_after_accumulate_max_abs, clear_after_accumulate_l2_norm
  real(mytype) :: double_accumulate_total_force_error, double_accumulate_relative_force_error
  real(mytype) :: buffer_superposition_max_error, buffer_superposition_l2_error
  real(mytype) :: feedback_buffer_total_force_error, feedback_buffer_relative_force_error
  real(mytype) :: feedback_buffer_max_abs, tmp_l2_norm

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
  allocate(u_lag(3, nlag), force_on_structure(3, nlag), force_on_fluid(3, nlag))
  allocate(force_a(3, nlag), force_b(3, nlag))

  ! Test A: allocation and zero init
  call allocate_ibm_force_buffer(buffer, grid)
  buffer_allocated_flag = 0
  if (buffer%is_allocated) buffer_allocated_flag = 1
  call compute_ibm_force_buffer_norms(buffer, buffer_initial_max_abs, buffer_initial_l2_norm)

  ! Test B: single accumulate conservation
  do l = 1, nlag
    lag%force(1, l) = -0.1_mytype - 0.02_mytype * real(l, mytype)
    lag%force(2, l) = 0.05_mytype - 0.01_mytype * real(l, mytype)
    lag%force(3, l) = -0.03_mytype + 0.005_mytype * real(l, mytype)
  end do
  call clear_ibm_force_buffer(buffer)
  call accumulate_lag_force_to_buffer(grid, lag, buffer)
  call compute_ibm_force_buffer_total_force(grid, buffer, buffer_total_force)
  lag_total_force = 0._mytype
  do l = 1, nlag
    lag_total_force = lag_total_force + lag%force(:, l) * lag%weight(l)
  end do
  single_accumulate_total_force_error = sqrt(sum((buffer_total_force - lag_total_force)**2))
  ref_norm = max(sqrt(sum(lag_total_force**2)), 1.0e-30_mytype)
  single_accumulate_relative_force_error = single_accumulate_total_force_error / ref_norm
  call compute_ibm_force_buffer_norms(buffer, single_accumulate_buffer_max_abs, single_accumulate_buffer_l2_norm)

  ! Test C: clear after accumulate
  call clear_ibm_force_buffer(buffer)
  call compute_ibm_force_buffer_norms(buffer, clear_after_accumulate_max_abs, clear_after_accumulate_l2_norm)

  ! Test D/E: double accumulate and array superposition
  do l = 1, nlag
    force_a(1, l) = -0.04_mytype - 0.01_mytype * real(l, mytype)
    force_a(2, l) =  0.03_mytype - 0.005_mytype * real(l, mytype)
    force_a(3, l) =  0.02_mytype + 0.002_mytype * real(l, mytype)

    force_b(1, l) = 0.06_mytype - 0.008_mytype * real(l, mytype)
    force_b(2, l) = -0.01_mytype + 0.004_mytype * real(l, mytype)
    force_b(3, l) = -0.03_mytype + 0.003_mytype * real(l, mytype)
  end do

  call allocate_ibm_force_buffer(buffer_a, grid)
  call allocate_ibm_force_buffer(buffer_b, grid)
  call allocate_ibm_force_buffer(buffer_ab, grid)

  lag%force = force_a
  call accumulate_lag_force_to_buffer(grid, lag, buffer_a)
  call compute_ibm_force_buffer_total_force(grid, buffer_a, total_a)

  lag%force = force_b
  call accumulate_lag_force_to_buffer(grid, lag, buffer_b)
  call compute_ibm_force_buffer_total_force(grid, buffer_b, total_b)

  lag%force = force_a
  call accumulate_lag_force_to_buffer(grid, lag, buffer_ab)
  lag%force = force_b
  call accumulate_lag_force_to_buffer(grid, lag, buffer_ab)
  call compute_ibm_force_buffer_total_force(grid, buffer_ab, total_ab)

  double_accumulate_total_force_error = sqrt(sum((total_ab - total_a - total_b)**2))
  ref_norm = max(sqrt(sum((total_a + total_b)**2)), 1.0e-30_mytype)
  double_accumulate_relative_force_error = double_accumulate_total_force_error / ref_norm

  buffer_superposition_max_error = maxval(abs(buffer_ab%fx - buffer_a%fx - buffer_b%fx))
  buffer_superposition_max_error = max(buffer_superposition_max_error, maxval(abs(buffer_ab%fy - buffer_a%fy - buffer_b%fy)))
  buffer_superposition_max_error = max(buffer_superposition_max_error, maxval(abs(buffer_ab%fz - buffer_a%fz - buffer_b%fz)))
  buffer_superposition_l2_error = sqrt(sum((buffer_ab%fx - buffer_a%fx - buffer_b%fx)**2 + &
                                           (buffer_ab%fy - buffer_a%fy - buffer_b%fy)**2 + &
                                           (buffer_ab%fz - buffer_a%fz - buffer_b%fz)**2))

  ! Test F: feedback force_on_fluid to buffer
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
  do l = 1, nlag
    lag%v(1, l) = 0.02_mytype * real(l, mytype)
    lag%v(2, l) = -0.01_mytype * real(l, mytype)
    lag%v(3, l) = 0.005_mytype * real(l, mytype)
  end do

  call interpolate_vector_to_lag(grid, ux, uy, uz, lag, u_lag)
  call compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)
  lag%force = force_on_fluid

  call clear_ibm_force_buffer(buffer)
  call accumulate_lag_force_to_buffer(grid, lag, buffer)
  call compute_ibm_force_buffer_total_force(grid, buffer, buffer_total_force)

  lag_total_force = 0._mytype
  do l = 1, nlag
    lag_total_force = lag_total_force + lag%force(:, l) * lag%weight(l)
  end do

  feedback_buffer_total_force_error = sqrt(sum((buffer_total_force - lag_total_force)**2))
  ref_norm = max(sqrt(sum(lag_total_force**2)), 1.0e-30_mytype)
  feedback_buffer_relative_force_error = feedback_buffer_total_force_error / ref_norm
  call compute_ibm_force_buffer_norms(buffer, feedback_buffer_max_abs, tmp_l2_norm)

  force_buffer_status = 1
  if (ieee_is_nan(buffer_initial_max_abs) .or. ieee_is_nan(buffer_initial_l2_norm) .or. &
      ieee_is_nan(single_accumulate_total_force_error) .or. ieee_is_nan(single_accumulate_relative_force_error) .or. &
      ieee_is_nan(single_accumulate_buffer_max_abs) .or. ieee_is_nan(single_accumulate_buffer_l2_norm) .or. &
      ieee_is_nan(clear_after_accumulate_max_abs) .or. ieee_is_nan(clear_after_accumulate_l2_norm) .or. &
      ieee_is_nan(double_accumulate_total_force_error) .or. ieee_is_nan(double_accumulate_relative_force_error) .or. &
      ieee_is_nan(buffer_superposition_max_error) .or. ieee_is_nan(buffer_superposition_l2_error) .or. &
      ieee_is_nan(feedback_buffer_total_force_error) .or. ieee_is_nan(feedback_buffer_relative_force_error) .or. &
      ieee_is_nan(feedback_buffer_max_abs)) then
    force_buffer_status = 0
  end if

  open(newunit=io_unit, file='stage3_outputs/fibre_ibm_force_buffer_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,I0)') 'buffer_allocated_flag = ', buffer_allocated_flag
  write(io_unit, '(A,I0)') 'buffer_nx = ', buffer%nx
  write(io_unit, '(A,I0)') 'buffer_ny = ', buffer%ny
  write(io_unit, '(A,I0)') 'buffer_nz = ', buffer%nz
  write(io_unit, '(A,ES24.16)') 'buffer_initial_max_abs = ', buffer_initial_max_abs
  write(io_unit, '(A,ES24.16)') 'buffer_initial_l2_norm = ', buffer_initial_l2_norm

  write(io_unit, '(A,ES24.16)') 'single_accumulate_total_force_error = ', single_accumulate_total_force_error
  write(io_unit, '(A,ES24.16)') 'single_accumulate_relative_force_error = ', single_accumulate_relative_force_error
  write(io_unit, '(A,ES24.16)') 'single_accumulate_buffer_max_abs = ', single_accumulate_buffer_max_abs
  write(io_unit, '(A,ES24.16)') 'single_accumulate_buffer_l2_norm = ', single_accumulate_buffer_l2_norm

  write(io_unit, '(A,ES24.16)') 'clear_after_accumulate_max_abs = ', clear_after_accumulate_max_abs
  write(io_unit, '(A,ES24.16)') 'clear_after_accumulate_l2_norm = ', clear_after_accumulate_l2_norm

  write(io_unit, '(A,ES24.16)') 'double_accumulate_total_force_error = ', double_accumulate_total_force_error
  write(io_unit, '(A,ES24.16)') 'double_accumulate_relative_force_error = ', double_accumulate_relative_force_error

  write(io_unit, '(A,ES24.16)') 'buffer_superposition_max_error = ', buffer_superposition_max_error
  write(io_unit, '(A,ES24.16)') 'buffer_superposition_l2_error = ', buffer_superposition_l2_error

  write(io_unit, '(A,ES24.16)') 'feedback_buffer_total_force_error = ', feedback_buffer_total_force_error
  write(io_unit, '(A,ES24.16)') 'feedback_buffer_relative_force_error = ', feedback_buffer_relative_force_error
  write(io_unit, '(A,ES24.16)') 'feedback_buffer_max_abs = ', feedback_buffer_max_abs

  write(io_unit, '(A,I0)') 'force_buffer_status = ', force_buffer_status

  close(io_unit)

  call destroy_ibm_force_buffer(buffer)
  call destroy_ibm_force_buffer(buffer_a)
  call destroy_ibm_force_buffer(buffer_b)
  call destroy_ibm_force_buffer(buffer_ab)
  deallocate(ux, uy, uz, u_lag, force_on_structure, force_on_fluid, force_a, force_b)
  call destroy_lagrangian_points(lag)
  call destroy_ibm_grid(grid)

end program fibre_ibm_force_buffer_check
