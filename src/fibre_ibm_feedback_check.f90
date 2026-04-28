program fibre_ibm_feedback_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_uniform_ibm_grid, destroy_ibm_grid
  use fibre_ibm_grid, only : destroy_lagrangian_points
  use fibre_ibm_interpolation, only : interpolate_vector_to_lag
  use fibre_ibm_spreading, only : clear_eulerian_force_density, spread_lag_force_to_eulerian
  use fibre_ibm_spreading, only : compute_eulerian_total_force
  use fibre_ibm_feedback, only : compute_ibm_feedback_forces, compute_feedback_slip_norms
  use fibre_ibm_feedback, only : compute_feedback_pair_power, compute_feedback_expected_dissipation

  implicit none

  type(ibm_grid_t) :: grid
  type(ibm_lagrangian_points_t) :: lag

  integer, parameter :: nlag = 5
  integer :: i, j, k, l, io_unit
  integer :: feedback_status

  real(mytype), parameter :: beta_drag = 10.0_mytype

  real(mytype), allocatable :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
  real(mytype), allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
  real(mytype), allocatable :: u_lag(:,:)
  real(mytype), allocatable :: force_on_structure(:,:), force_on_fluid(:,:)

  real(mytype) :: x, y, z
  real(mytype) :: tmp_rms_slip
  real(mytype) :: ref_norm

  real(mytype) :: total_structure_force(3), total_fluid_force(3), total_eulerian_force(3)
  real(mytype) :: total_lagrangian_fluid_force(3)

  real(mytype) :: zero_slip_force_structure_norm, zero_slip_force_fluid_norm
  real(mytype) :: zero_slip_max_slip_norm, zero_slip_power_total

  real(mytype) :: constant_slip_structure_force_error, constant_slip_fluid_force_error
  real(mytype) :: constant_slip_action_reaction_error

  real(mytype) :: action_reaction_total_force_error

  real(mytype) :: feedback_power_structure, feedback_power_fluid, feedback_power_total
  real(mytype) :: feedback_expected_dissipation, feedback_dissipation_identity_error

  real(mytype) :: feedback_spread_total_force_error, feedback_spread_relative_force_error

  real(mytype) :: interpolated_zero_slip_force_norm

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
  allocate(u_lag(3, nlag), force_on_structure(3, nlag), force_on_fluid(3, nlag))

  ! Test A: zero slip
  do l = 1, nlag
    u_lag(1, l) = 0.12_mytype + 0.01_mytype * real(l, mytype)
    u_lag(2, l) = -0.08_mytype + 0.02_mytype * real(l, mytype)
    u_lag(3, l) = 0.03_mytype * sin(0.4_mytype * real(l, mytype))
  end do
  lag%v = u_lag
  call compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)
  zero_slip_force_structure_norm = sqrt(sum(force_on_structure**2))
  zero_slip_force_fluid_norm = sqrt(sum(force_on_fluid**2))
  call compute_feedback_slip_norms(lag, u_lag, zero_slip_max_slip_norm, tmp_rms_slip)
  call compute_feedback_pair_power(lag, u_lag, force_on_structure, force_on_fluid, &
                                   feedback_power_structure, feedback_power_fluid, zero_slip_power_total)

  ! Test B: constant slip sign and magnitude
  lag%v = 0._mytype
  do l = 1, nlag
    u_lag(:, l) = [0.2_mytype, -0.1_mytype, 0.05_mytype]
  end do
  call compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)

  constant_slip_structure_force_error = 0._mytype
  constant_slip_fluid_force_error = 0._mytype
  do l = 1, nlag
    constant_slip_structure_force_error = max(constant_slip_structure_force_error, &
      maxval(abs(force_on_structure(:, l) - beta_drag * [0.2_mytype, -0.1_mytype, 0.05_mytype])))
    constant_slip_fluid_force_error = max(constant_slip_fluid_force_error, &
      maxval(abs(force_on_fluid(:, l) + beta_drag * [0.2_mytype, -0.1_mytype, 0.05_mytype])))
  end do
  constant_slip_action_reaction_error = maxval(abs(force_on_structure + force_on_fluid))

  ! Test C: action-reaction total force
  do l = 1, nlag
    u_lag(1, l) = 0.1_mytype + 0.02_mytype * real(l, mytype)
    u_lag(2, l) = -0.2_mytype + 0.01_mytype * real(l, mytype)
    u_lag(3, l) = 0.05_mytype * sin(0.3_mytype * real(l, mytype))

    lag%v(1, l) = -0.03_mytype + 0.01_mytype * real(l, mytype)
    lag%v(2, l) = 0.04_mytype * cos(0.2_mytype * real(l, mytype))
    lag%v(3, l) = 0.02_mytype
  end do
  call compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)

  total_structure_force = 0._mytype
  total_fluid_force = 0._mytype
  do l = 1, nlag
    total_structure_force = total_structure_force + force_on_structure(:, l) * lag%weight(l)
    total_fluid_force = total_fluid_force + force_on_fluid(:, l) * lag%weight(l)
  end do
  action_reaction_total_force_error = sqrt(sum((total_structure_force + total_fluid_force)**2))

  ! Test D: relative-slip dissipation power
  call compute_feedback_pair_power(lag, u_lag, force_on_structure, force_on_fluid, &
                                   feedback_power_structure, feedback_power_fluid, feedback_power_total)
  call compute_feedback_expected_dissipation(lag, u_lag, beta_drag, feedback_expected_dissipation)
  feedback_dissipation_identity_error = abs(feedback_power_total + feedback_expected_dissipation)

  ! Test E: force_on_fluid spreading total-force conservation
  lag%force = force_on_fluid
  call clear_eulerian_force_density(fx, fy, fz)
  call spread_lag_force_to_eulerian(grid, lag, fx, fy, fz)
  call compute_eulerian_total_force(grid, fx, fy, fz, total_eulerian_force)

  total_lagrangian_fluid_force = 0._mytype
  do l = 1, nlag
    total_lagrangian_fluid_force = total_lagrangian_fluid_force + force_on_fluid(:, l) * lag%weight(l)
  end do

  feedback_spread_total_force_error = sqrt(sum((total_eulerian_force - total_lagrangian_fluid_force)**2))
  ref_norm = max(sqrt(sum(total_lagrangian_fluid_force**2)), 1.0e-30_mytype)
  feedback_spread_relative_force_error = feedback_spread_total_force_error / ref_norm

  ! Test F: interpolated field linked zero slip
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
  lag%v = u_lag
  call compute_ibm_feedback_forces(lag, u_lag, beta_drag, force_on_structure, force_on_fluid)
  interpolated_zero_slip_force_norm = sqrt(sum(force_on_structure**2 + force_on_fluid**2))

  feedback_status = 1
  if (ieee_is_nan(zero_slip_force_structure_norm) .or. ieee_is_nan(zero_slip_force_fluid_norm) .or. &
      ieee_is_nan(zero_slip_max_slip_norm) .or. ieee_is_nan(zero_slip_power_total) .or. &
      ieee_is_nan(constant_slip_structure_force_error) .or. ieee_is_nan(constant_slip_fluid_force_error) .or. &
      ieee_is_nan(constant_slip_action_reaction_error) .or. ieee_is_nan(action_reaction_total_force_error) .or. &
      ieee_is_nan(feedback_power_structure) .or. ieee_is_nan(feedback_power_fluid) .or. &
      ieee_is_nan(feedback_power_total) .or. ieee_is_nan(feedback_expected_dissipation) .or. &
      ieee_is_nan(feedback_dissipation_identity_error) .or. ieee_is_nan(feedback_spread_total_force_error) .or. &
      ieee_is_nan(feedback_spread_relative_force_error) .or. ieee_is_nan(interpolated_zero_slip_force_norm)) then
    feedback_status = 0
  end if

  open(newunit=io_unit, file='stage3_outputs/fibre_ibm_feedback_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,ES24.16)') 'zero_slip_force_structure_norm = ', zero_slip_force_structure_norm
  write(io_unit, '(A,ES24.16)') 'zero_slip_force_fluid_norm = ', zero_slip_force_fluid_norm
  write(io_unit, '(A,ES24.16)') 'zero_slip_max_slip_norm = ', zero_slip_max_slip_norm
  write(io_unit, '(A,ES24.16)') 'zero_slip_power_total = ', zero_slip_power_total

  write(io_unit, '(A,ES24.16)') 'constant_slip_structure_force_error = ', constant_slip_structure_force_error
  write(io_unit, '(A,ES24.16)') 'constant_slip_fluid_force_error = ', constant_slip_fluid_force_error
  write(io_unit, '(A,ES24.16)') 'constant_slip_action_reaction_error = ', constant_slip_action_reaction_error

  write(io_unit, '(A,ES24.16)') 'action_reaction_total_force_error = ', action_reaction_total_force_error

  write(io_unit, '(A,ES24.16)') 'feedback_power_structure = ', feedback_power_structure
  write(io_unit, '(A,ES24.16)') 'feedback_power_fluid = ', feedback_power_fluid
  write(io_unit, '(A,ES24.16)') 'feedback_power_total = ', feedback_power_total
  write(io_unit, '(A,ES24.16)') 'feedback_expected_dissipation = ', feedback_expected_dissipation
  write(io_unit, '(A,ES24.16)') 'feedback_dissipation_identity_error = ', feedback_dissipation_identity_error

  write(io_unit, '(A,ES24.16)') 'feedback_spread_total_force_error = ', feedback_spread_total_force_error
  write(io_unit, '(A,ES24.16)') 'feedback_spread_relative_force_error = ', feedback_spread_relative_force_error

  write(io_unit, '(A,ES24.16)') 'interpolated_zero_slip_force_norm = ', interpolated_zero_slip_force_norm

  write(io_unit, '(A,I0)') 'feedback_status = ', feedback_status

  close(io_unit)

  deallocate(ux, uy, uz, fx, fy, fz, u_lag, force_on_structure, force_on_fluid)
  call destroy_lagrangian_points(lag)
  call destroy_ibm_grid(grid)

end program fibre_ibm_feedback_check
