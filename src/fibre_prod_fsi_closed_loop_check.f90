program fibre_prod_fsi_closed_loop_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, &
                                      fibre_prod_grid_destroy
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, &
                               fibre_prod_state_init_straight, fibre_prod_state_destroy, &
                               fibre_prod_state_all_finite
  use fibre_prod_fsi_config, only : fibre_prod_fsi_config_type, fibre_prod_fsi_config_default, &
                                    fibre_prod_fsi_config_validate
  use fibre_prod_fluid_surrogate, only : fibre_prod_fluid_surrogate_type, &
                                         fibre_prod_fluid_surrogate_allocate, &
                                         fibre_prod_fluid_surrogate_reset_constant, &
                                         fibre_prod_fluid_surrogate_destroy, &
                                         fibre_prod_fluid_surrogate_is_finite
  use fibre_prod_fsi_coupling, only : fibre_prod_fsi_closed_loop_step
  use fibre_prod_fsi_diagnostics, only : fibre_prod_fsi_force_norm, &
                                         fibre_prod_fsi_action_reaction_residual, &
                                         fibre_prod_fsi_structure_kinetic_energy, &
                                         fibre_prod_fsi_fluid_kinetic_energy
  implicit none

  integer, parameter :: dp = real64
  type(fibre_prod_grid_type) :: grid
  type(fibre_prod_state_type) :: state
  type(fibre_prod_fluid_surrogate_type) :: fluid
  type(fibre_prod_fsi_config_type) :: config
  real(dp) :: x_coord(8)
  real(dp) :: y_coord(6)
  real(dp) :: z_coord(7)
  real(dp) :: total_fibre_force(3)
  real(dp) :: reaction_integral(3)
  real(dp), allocatable :: x_before(:, :, :)
  real(dp), allocatable :: v_before(:, :, :)
  real(dp), allocatable :: u_before(:, :, :)
  real(dp) :: action_residual
  real(dp) :: structure_energy
  real(dp) :: fluid_energy
  integer :: status
  integer :: istep

  call setup_grid_state_fluid(grid, state, fluid, status)
  if (status /= 0) error stop 1

  config = fibre_prod_fsi_config_default()
  config%fsi_enabled = .false.
  config%lambda_fsi = 0.0_dp
  config%penalty_beta = 2.0_dp
  config%dt = 1.0e-3_dp
  config%rho_tilde = 1.0_dp
  config%gamma = 0.0_dp
  config%ds = 0.05_dp
  config%nstep = 1
  if (fibre_prod_fsi_config_validate(config) /= 0) error stop 2

  allocate(x_before(state%nfibre, state%nnode, 3), v_before(state%nfibre, state%nnode, 3))
  allocate(u_before(fluid%nx_local, fluid%ny_local, fluid%nz_local))
  x_before = state%x
  v_before = state%v
  u_before = fluid%u
  call fibre_prod_fsi_closed_loop_step(grid, fluid, state, config, status, total_fibre_force, reaction_integral)
  if (status /= 0) error stop 3
  if (fibre_prod_fsi_force_norm(total_fibre_force) > 1.0e-14_dp) error stop 4
  if (fibre_prod_fsi_force_norm(reaction_integral) > 1.0e-14_dp) error stop 5
  if (maxval(abs(state%x - x_before)) > 1.0e-12_dp) error stop 6
  if (maxval(abs(state%v - v_before)) > 1.0e-12_dp) error stop 7
  if (maxval(abs(fluid%u - u_before)) > 1.0e-12_dp) error stop 8
  if (.not. fibre_prod_state_all_finite(state) .or. .not. fibre_prod_fluid_surrogate_is_finite(fluid)) error stop 9

  call reset_state_fluid(state, fluid, status)
  if (status /= 0) error stop 10
  config%fsi_enabled = .true.
  config%lambda_fsi = 0.1_dp
  config%nstep = 3
  call fibre_prod_fsi_closed_loop_step(grid, fluid, state, config, status, total_fibre_force, reaction_integral)
  if (status /= 0) error stop 11
  if (fibre_prod_fsi_force_norm(total_fibre_force) <= 0.0_dp) error stop 12
  if (state%v(1, 3, 1) <= 0.0_dp) error stop 13
  action_residual = fibre_prod_fsi_action_reaction_residual(total_fibre_force, reaction_integral)
  if (action_residual > 1.0e-10_dp) error stop 14
  if (fluid%u(1, 1, 1) >= 1.0_dp) error stop 15
  do istep = 2, config%nstep
    call fibre_prod_fsi_closed_loop_step(grid, fluid, state, config, status, total_fibre_force, reaction_integral)
    if (status /= 0) error stop 16
    if (.not. fibre_prod_state_all_finite(state) .or. .not. fibre_prod_fluid_surrogate_is_finite(fluid)) error stop 17
  end do
  structure_energy = fibre_prod_fsi_structure_kinetic_energy(state, config%rho_tilde)
  fluid_energy = fibre_prod_fsi_fluid_kinetic_energy(fluid, grid)
  if (structure_energy < 0.0_dp .or. fluid_energy < 0.0_dp) error stop 18
  if (maxval(abs(state%v)) > 1.0_dp) error stop 19

  config%lambda_fsi = -0.1_dp
  if (fibre_prod_fsi_config_validate(config) == 0) error stop 20
  config%lambda_fsi = 0.1_dp
  config%dt = -1.0e-3_dp
  if (fibre_prod_fsi_config_validate(config) == 0) error stop 21
  config%dt = 1.0e-3_dp
  state%x(1, 3, 2) = -0.1_dp
  call fibre_prod_fsi_closed_loop_step(grid, fluid, state, config, status, total_fibre_force, reaction_integral)
  if (status == 0) error stop 22
  call reset_state_fluid(state, fluid, status)
  fluid%u(1, 1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
  call fibre_prod_fsi_closed_loop_step(grid, fluid, state, config, status, total_fibre_force, reaction_integral)
  if (status == 0) error stop 23

  call fibre_prod_state_destroy(state)
  call fibre_prod_fluid_surrogate_destroy(fluid)
  call fibre_prod_grid_destroy(grid)
  print *, 'R7_FIBRE_PROD_FSI_CLOSED_LOOP_CHECK PASS'

contains

  subroutine setup_grid_state_fluid(grid, state, fluid, status)
    type(fibre_prod_grid_type), intent(out) :: grid
    type(fibre_prod_state_type), intent(out) :: state
    type(fibre_prod_fluid_surrogate_type), intent(out) :: fluid
    integer, intent(out) :: status
    real(dp) :: origin(3)
    real(dp) :: direction(3)

    x_coord = [0.0_dp, 0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp, 1.0_dp, 1.2_dp, 1.4_dp]
    y_coord = [0.0_dp, 0.03_dp, 0.12_dp, 0.30_dp, 0.62_dp, 1.0_dp]
    z_coord = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp, 1.25_dp, 1.5_dp]
    call fibre_prod_grid_init_from_coordinates(grid, x_coord, y_coord, z_coord, &
                                               1, 8, 1, 6, 1, 7, &
                                               .true., .false., .true., status)
    if (status /= 0) return
    call fibre_prod_state_allocate(state, 1, 5, status)
    if (status /= 0) return
    origin = [0.30_dp, 0.20_dp, 0.50_dp]
    direction = [1.0_dp, 0.0_dp, 0.0_dp]
    call fibre_prod_state_init_straight(state, origin, direction, 0.05_dp, stat=status)
    if (status /= 0) return
    call fibre_prod_fluid_surrogate_allocate(fluid, grid, status)
    if (status /= 0) return
    call fibre_prod_fluid_surrogate_reset_constant(fluid, [1.0_dp, 0.0_dp, 0.0_dp], status)
  end subroutine setup_grid_state_fluid

  subroutine reset_state_fluid(state, fluid, status)
    type(fibre_prod_state_type), intent(inout) :: state
    type(fibre_prod_fluid_surrogate_type), intent(inout) :: fluid
    integer, intent(out) :: status
    real(dp) :: origin(3)
    real(dp) :: direction(3)

    origin = [0.30_dp, 0.20_dp, 0.50_dp]
    direction = [1.0_dp, 0.0_dp, 0.0_dp]
    call fibre_prod_state_init_straight(state, origin, direction, 0.05_dp, stat=status)
    if (status /= 0) return
    call fibre_prod_fluid_surrogate_reset_constant(fluid, [1.0_dp, 0.0_dp, 0.0_dp], status)
  end subroutine reset_state_fluid

end program fibre_prod_fsi_closed_loop_check
