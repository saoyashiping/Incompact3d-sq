program fibre_prod_reaction_spreading_buffer_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, fibre_prod_grid_destroy
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, fibre_prod_force_buffer_allocate, &
                                          fibre_prod_force_buffer_destroy
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, fibre_prod_state_destroy, &
                              fibre_prod_state_attach_sampled_velocity, fibre_prod_state_get_structure_coordinates, &
                              fibre_prod_state_get_structure_velocity_or_zero, fibre_prod_state_get_structure_input_force
  use fibre_prod_hydro_input_candidate, only : fibre_prod_hydro_input_candidate_type, &
                                               fibre_prod_hydro_input_candidate_init, &
                                               fibre_prod_hydro_input_candidate_compute, &
                                               fibre_prod_hydro_input_candidate_attach_to_state, &
                                               fibre_prod_hydro_input_candidate_finalize
  use fibre_prod_structure_input_handoff, only : fibre_prod_structure_input_handoff_type, &
                                                 fibre_prod_structure_input_handoff_init, &
                                                 fibre_prod_structure_input_handoff_from_candidate, &
                                                 fibre_prod_structure_input_handoff_attach_to_state, &
                                                 fibre_prod_structure_input_handoff_finalize
  use fibre_prod_structure_dry_step, only : fibre_prod_structure_dry_step_type, &
                                            fibre_prod_structure_dry_step_init, &
                                            fibre_prod_structure_dry_step_predict, &
                                            fibre_prod_structure_dry_step_finalize
  use fibre_prod_structure_commit_gate, only : fibre_prod_structure_commit_gate_type, &
                                               fibre_prod_structure_commit_gate_init, &
                                               fibre_prod_structure_commit_gate_set_enabled, &
                                               fibre_prod_structure_commit_gate_evaluate, &
                                               fibre_prod_structure_commit_gate_commit_to_state, &
                                               fibre_prod_structure_commit_gate_finalize
  use fibre_prod_reaction_force_candidate, only : fibre_prod_reaction_force_candidate_type, &
                                                  fibre_prod_reaction_force_candidate_init, &
                                                  fibre_prod_reaction_force_candidate_from_structure_input, &
                                                  fibre_prod_reaction_force_candidate_finalize
  use fibre_prod_reaction_spreading_buffer, only : fibre_prod_reaction_spreading_buffer_type, &
                                                   fibre_prod_reaction_spreading_buffer_init, &
                                                   fibre_prod_reaction_spreading_buffer_apply, &
                                                   fibre_prod_reaction_spreading_buffer_check_finite_bounded, &
                                                   fibre_prod_reaction_spreading_buffer_finalize
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, fibre_prod_runtime_config_default
  use fibre_prod_main_hook, only : fibre_prod_main_hook_init
  use fibre_prod_runtime_bridge, only : fibre_prod_runtime_bridge_type, fibre_prod_runtime_bridge_init_from_rhs, &
                                        fibre_prod_runtime_bridge_apply_lambda0_noop, fibre_prod_runtime_bridge_finalize
  implicit none

  integer, parameter :: dp = real64
  integer :: status
  type(fibre_prod_grid_type) :: grid
  type(fibre_prod_force_buffer_type) :: buffer
  type(fibre_prod_force_buffer_type) :: unallocated_buffer
  type(fibre_prod_state_type) :: state
  type(fibre_prod_hydro_input_candidate_type) :: hydro
  type(fibre_prod_structure_input_handoff_type) :: handoff
  type(fibre_prod_structure_dry_step_type) :: dry
  type(fibre_prod_structure_commit_gate_type) :: gate
  type(fibre_prod_reaction_force_candidate_type) :: react
  type(fibre_prod_reaction_force_candidate_type) :: uninit_react
  type(fibre_prod_reaction_spreading_buffer_type) :: spread
  type(fibre_prod_runtime_bridge_type) :: bridge
  type(fibre_prod_runtime_config_type) :: config
  real(dp) :: coord(9)
  real(dp) :: structure_input(3,3), x_before(3,3), structure_u_before(3,3)
  real(dp) :: sampled_before(3,3), hydro_before(3,3), input_before(3,3)
  real(dp) :: rhs_x(3,3,3), rhs_y(3,3,3), rhs_z(3,3,3)
  real(dp) :: rhs_x0(3,3,3), rhs_y0(3,3,3), rhs_z0(3,3,3)
  real(dp) :: ux(3,3,3), uy(3,3,3), uz(3,3,3), ux0(3,3,3), uy0(3,3,3), uz0(3,3,3)
  real(dp) :: bad_points(2,3), bad_force(2,3), nonfinite_force(3,3)
  integer :: i

  do i = 1, 9
    coord(i) = real(i - 1, dp) / 8.0_dp
  end do
  call fibre_prod_grid_init_from_coordinates(grid, coord, coord, coord, 1, 9, 1, 9, 1, 9, .false., .false., .false., status)
  call require(status == 0, 'grid setup failed')
  call setup_pipeline(state, hydro, handoff, dry, gate, structure_input, status)
  call require(status == 0, 'structure-side setup failed')
  call fibre_prod_state_get_structure_coordinates(state, x_before, status)
  call require(status == 0, 'coordinate read failed')
  call fibre_prod_state_get_structure_velocity_or_zero(state, structure_u_before, status)
  call require(status == 0, 'structure velocity read failed')
  sampled_before = state%sampled_u
  hydro_before = state%hydro_force_candidate
  input_before = state%structure_input_force

  call fibre_prod_reaction_force_candidate_init(react, 3, status)
  call require(status == 0, 'reaction candidate init failed')
  call fibre_prod_reaction_force_candidate_from_structure_input(react, structure_input, status)
  call require(status == 0, 'reaction candidate compute failed')
  call require(maxval(abs(react%reaction_force + structure_input)) <= 0.0_dp, 'reaction sign convention failed')
  call require(react%has_reaction_force, 'reaction candidate flag missing')
  call require(abs(react%max_abs_reaction_force - maxval(abs(react%reaction_force))) <= 0.0_dp, 'max diagnostic failed')
  call require(abs(react%sum_abs_reaction_force - sum(abs(react%reaction_force))) <= 0.0_dp, 'sum diagnostic failed')
  call require(maxval(abs(react%net_reaction_force - sum(react%reaction_force, dim=1))) <= 0.0_dp, &
               'net diagnostic failed')

  call fibre_prod_force_buffer_allocate(buffer, grid, status)
  call require(status == 0, 'force buffer allocation failed')
  call fibre_prod_reaction_spreading_buffer_init(spread, grid%nx_local, grid%ny_local, grid%nz_local, 3, status)
  call require(status == 0, 'spreading buffer init failed')
  call fibre_prod_reaction_spreading_buffer_apply(spread, grid, x_before, react%reaction_force, buffer, status)
  call require(status == 0 .and. spread%has_spread, 'reaction spreading failed')
  call require(maxval(abs(buffer%fx)) + maxval(abs(buffer%fy)) + maxval(abs(buffer%fz)) > 0.0_dp, &
               'Eulerian force buffer is zero')
  call require(spread%max_abs_force_buffer > 0.0_dp .and. spread%sum_abs_force_buffer > 0.0_dp, &
               'force buffer diagnostics missing')
  call fibre_prod_reaction_spreading_buffer_check_finite_bounded(spread, 1.0e6_dp, status)
  call require(status == 0, 'finite bounded check failed')
  call require(maxval(abs(spread%conservation_error)) <= 1.0e-8_dp, 'spreading conservation failed')

  call require(all(state%x(1,:,:) == x_before), 'coordinates changed')
  call require(all(state%structure_u == structure_u_before), 'structure velocity changed')
  call require(all(state%sampled_u == sampled_before), 'sampled velocity changed')
  call require(all(state%hydro_force_candidate == hydro_before), 'hydro candidate changed')
  call require(all(state%structure_input_force == input_before), 'structure input changed')

  ux = reshape([(real(i, dp), i=1,27)], shape(ux))
  uy = ux + 5.0_dp
  uz = ux - 5.0_dp
  ux0 = ux; uy0 = uy; uz0 = uz
  rhs_x = ux * 0.1_dp; rhs_y = uy * 0.1_dp; rhs_z = uz * 0.1_dp
  rhs_x0 = rhs_x; rhs_y0 = rhs_y; rhs_z0 = rhs_z
  call require(all(ux == ux0) .and. all(uy == uy0) .and. all(uz == uz0), 'velocity fields changed')
  call require(all(rhs_x == rhs_x0) .and. all(rhs_y == rhs_y0) .and. all(rhs_z == rhs_z0), &
               'RHS changed before runtime bridge')

  call fibre_prod_runtime_bridge_init_from_rhs(rhs_x, rhs_y, rhs_z, bridge, status)
  call require(status == 0, 'runtime bridge init failed')
  call fibre_prod_runtime_config_default(config)
  config%enabled = .true.; config%lambda_fsi = 0.0_dp; config%penalty_beta = 2.0_dp
  call fibre_prod_main_hook_init(status, config)
  call require(status == 0, 'lambda0 main hook init failed')
  call fibre_prod_runtime_bridge_apply_lambda0_noop(rhs_x, rhs_y, rhs_z, bridge, status)
  call require(status == 0, 'lambda0 runtime path failed')
  call require(all(rhs_x == rhs_x0) .and. all(rhs_y == rhs_y0) .and. all(rhs_z == rhs_z0), 'lambda0 RHS changed')
  config%lambda_fsi = 1.0e-3_dp
  call fibre_prod_main_hook_init(status, config)
  call require(status == 0, 'lambda positive main hook init failed')
  call fibre_prod_reaction_spreading_buffer_apply(spread, grid, x_before, react%reaction_force, buffer, status)
  call require(status == 0, 'lambda positive candidate spreading failed')
  call require(all(rhs_x == rhs_x0) .and. all(rhs_y == rhs_y0) .and. all(rhs_z == rhs_z0), &
               'lambda positive without RHS gate changed RHS')

  call fibre_prod_reaction_force_candidate_from_structure_input(uninit_react, structure_input, status)
  call require(status /= 0, 'uninitialized candidate guard failed')
  bad_points = 0.0_dp
  call fibre_prod_reaction_spreading_buffer_apply(spread, grid, bad_points, react%reaction_force, buffer, status)
  call require(status /= 0, 'bad fibre point shape guard failed')
  bad_force = 0.0_dp
  call fibre_prod_reaction_spreading_buffer_apply(spread, grid, x_before, bad_force, buffer, status)
  call require(status /= 0, 'bad reaction force shape guard failed')
  call fibre_prod_reaction_spreading_buffer_apply(spread, grid, x_before, react%reaction_force, unallocated_buffer, status)
  call require(status /= 0, 'unallocated force buffer guard failed')
  nonfinite_force = react%reaction_force
  nonfinite_force(1,1) = ieee_value(0.0_dp, ieee_quiet_nan)
  call fibre_prod_reaction_spreading_buffer_apply(spread, grid, x_before, nonfinite_force, buffer, status)
  call require(status /= 0, 'nonfinite reaction force guard failed')

  print *, 'P0_10_DIAGNOSTIC: reaction spreading buffer only, no RHS feedback in P0_10'
  print *, 'P0_10_REACTION_SPREADING_BUFFER_CHECK PASS'

  call fibre_prod_runtime_bridge_finalize(bridge)
  call fibre_prod_reaction_spreading_buffer_finalize(spread)
  call fibre_prod_reaction_force_candidate_finalize(react)
  call fibre_prod_force_buffer_destroy(buffer)
  call fibre_prod_grid_destroy(grid)
  call fibre_prod_structure_commit_gate_finalize(gate)
  call fibre_prod_structure_dry_step_finalize(dry)
  call fibre_prod_structure_input_handoff_finalize(handoff)
  call fibre_prod_hydro_input_candidate_finalize(hydro)
  call fibre_prod_state_destroy(state)
contains
  subroutine setup_pipeline(state, hydro, handoff, dry, gate, structure_input, status)
    type(fibre_prod_state_type), intent(inout) :: state
    type(fibre_prod_hydro_input_candidate_type), intent(inout) :: hydro
    type(fibre_prod_structure_input_handoff_type), intent(inout) :: handoff
    type(fibre_prod_structure_dry_step_type), intent(inout) :: dry
    type(fibre_prod_structure_commit_gate_type), intent(inout) :: gate
    real(dp), intent(out) :: structure_input(3,3)
    integer, intent(out) :: status
    real(dp) :: points(3,3), sampled_u(3,3), structure_u(3,3), state_x(3,3)

    points(1,:) = [0.25_dp, 0.5_dp, 0.5_dp]
    points(2,:) = [0.50_dp, 0.5_dp, 0.5_dp]
    points(3,:) = [0.75_dp, 0.5_dp, 0.5_dp]
    sampled_u(1,:) = [1.0_dp, 0.5_dp, -0.25_dp]
    sampled_u(2,:) = [1.5_dp, 0.25_dp, -0.50_dp]
    sampled_u(3,:) = [2.0_dp, 0.0_dp, -0.75_dp]
    call fibre_prod_state_allocate(state, 1, 3, status)
    if (status /= 0) return
    state%x(1,:,:) = points
    call fibre_prod_state_attach_sampled_velocity(state, sampled_u, status)
    if (status /= 0) return
    call fibre_prod_state_get_structure_velocity_or_zero(state, structure_u, status)
    if (status /= 0) return
    call fibre_prod_hydro_input_candidate_init(hydro, 3, 2.5_dp, status)
    if (status /= 0) return
    call fibre_prod_hydro_input_candidate_compute(hydro, state%sampled_u, structure_u, status)
    if (status /= 0) return
    call fibre_prod_hydro_input_candidate_attach_to_state(hydro, state, status)
    if (status /= 0) return
    call fibre_prod_structure_input_handoff_init(handoff, 3, status)
    if (status /= 0) return
    call fibre_prod_structure_input_handoff_from_candidate(handoff, state%hydro_force_candidate, status)
    if (status /= 0) return
    call fibre_prod_structure_input_handoff_attach_to_state(handoff, state, status)
    if (status /= 0) return
    call fibre_prod_state_get_structure_coordinates(state, state_x, status)
    if (status /= 0) return
    call fibre_prod_structure_dry_step_init(dry, 3, 1.0e-4_dp, 2.0_dp, status)
    if (status /= 0) return
    call fibre_prod_structure_dry_step_predict(dry, state_x, structure_u, state%structure_input_force, status)
    if (status /= 0) return
    call fibre_prod_structure_commit_gate_init(gate, 3, 1.0e-2_dp, status)
    if (status /= 0) return
    call fibre_prod_structure_commit_gate_set_enabled(gate, .true., status)
    if (status /= 0) return
    call fibre_prod_structure_commit_gate_evaluate(gate, dry%x_trial, dry%u_trial, dry%dx_trial, status)
    if (status /= 0) return
    call fibre_prod_structure_commit_gate_commit_to_state(gate, state, dry%x_trial, dry%u_trial, status)
    if (status /= 0) return
    call fibre_prod_state_get_structure_input_force(state, structure_input, status)
  end subroutine setup_pipeline

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      print *, 'P0_10_REACTION_SPREADING_BUFFER_CHECK FAIL: ', trim(message)
      error stop 1
    end if
  end subroutine require
end program fibre_prod_reaction_spreading_buffer_check
