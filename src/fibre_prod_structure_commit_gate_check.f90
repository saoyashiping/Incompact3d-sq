program fibre_prod_structure_commit_gate_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, fibre_prod_grid_destroy
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, fibre_prod_force_buffer_allocate, &
                                          fibre_prod_force_buffer_destroy
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, fibre_prod_state_destroy, &
                              fibre_prod_state_attach_sampled_velocity, fibre_prod_state_get_structure_coordinates, &
                              fibre_prod_state_get_structure_velocity_or_zero
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
                                               fibre_prod_structure_commit_gate_finalize, &
                                               fibre_prod_commit_reject_disabled, fibre_prod_commit_reject_nonfinite, &
                                               fibre_prod_commit_reject_bound, fibre_prod_commit_reject_shape
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, fibre_prod_runtime_config_default
  use fibre_prod_main_hook, only : fibre_prod_main_hook_init
  use fibre_prod_runtime_bridge, only : fibre_prod_runtime_bridge_type, fibre_prod_runtime_bridge_init_from_rhs, &
                                        fibre_prod_runtime_bridge_apply_lambda0_noop, fibre_prod_runtime_bridge_finalize
  implicit none

  integer, parameter :: dp = real64
  integer :: status
  type(fibre_prod_state_type) :: state
  type(fibre_prod_hydro_input_candidate_type) :: candidate
  type(fibre_prod_structure_input_handoff_type) :: handoff
  type(fibre_prod_structure_dry_step_type) :: dry
  type(fibre_prod_structure_commit_gate_type) :: gate
  type(fibre_prod_grid_type) :: grid
  type(fibre_prod_force_buffer_type) :: force_buffer
  type(fibre_prod_runtime_bridge_type) :: bridge
  type(fibre_prod_runtime_config_type) :: config
  real(dp) :: x0(3,3), sampled_u0(3,3), hydro0(3,3), input0(3,3), state_u0(3,3)
  real(dp) :: rhs_x(2,2,2), rhs_y(2,2,2), rhs_z(2,2,2)
  real(dp) :: rhs_x0(2,2,2), rhs_y0(2,2,2), rhs_z0(2,2,2)
  real(dp) :: ux(2,2,2), uy(2,2,2), uz(2,2,2), ux0(2,2,2), uy0(2,2,2), uz0(2,2,2)
  real(dp) :: coord(2), bad_x(2,3), nan_x(3,3), huge_dx(3,3)
  integer :: i

  call setup_state(state, candidate, handoff, dry, x0, sampled_u0, hydro0, input0, state_u0, status)
  call require(status == 0, 'deterministic setup failed')
  call fibre_prod_structure_commit_gate_init(gate, 3, 1.0e-2_dp, status)
  call require(status == 0 .and. gate%initialized, 'gate init failed')
  call require(.not. gate%gate_enabled, 'gate must default disabled')

  call fibre_prod_structure_commit_gate_evaluate(gate, dry%x_trial, dry%u_trial, dry%dx_trial, status)
  call require(status == 0 .and. gate%rejected .and. gate%reject_code == fibre_prod_commit_reject_disabled, &
               'disabled gate should reject without fatal status')
  call fibre_prod_structure_commit_gate_commit_to_state(gate, state, dry%x_trial, dry%u_trial, status)
  call require(status == 0, 'disabled no-commit should be nonfatal')
  call require(all(state%x(1,:,:) == x0), 'disabled gate changed coordinates')
  call require(.not. state%has_structure_u, 'disabled gate attached structure velocity')

  call fibre_prod_structure_commit_gate_set_enabled(gate, .true., status)
  call require(status == 0, 'enable gate failed')
  call fibre_prod_structure_commit_gate_evaluate(gate, dry%x_trial, dry%u_trial, dry%dx_trial, status)
  call require(status == 0 .and. gate%accepted, 'enabled bounded trial should be accepted')
  call fibre_prod_structure_commit_gate_commit_to_state(gate, state, dry%x_trial, dry%u_trial, status)
  call require(status == 0 .and. gate%committed, 'accepted trial commit failed')
  call require(maxval(abs(state%x(1,:,:) - dry%x_trial)) <= 0.0_dp, 'coordinates not committed')
  call require(state%has_structure_u .and. maxval(abs(state%structure_u - dry%u_trial)) <= 0.0_dp, &
               'structure velocity not committed')
  call require(all(state%sampled_u == sampled_u0), 'sampled velocity changed')
  call require(all(state%hydro_force_candidate == hydro0), 'hydro candidate changed')
  call require(all(state%structure_input_force == input0), 'structure input changed')

  state%x(1,:,:) = x0
  call fibre_prod_structure_commit_gate_finalize(gate)
  call fibre_prod_structure_commit_gate_init(gate, 3, 1.0e-2_dp, status)
  call fibre_prod_structure_commit_gate_set_enabled(gate, .true., status)
  nan_x = dry%x_trial
  nan_x(1,1) = ieee_value(0.0_dp, ieee_quiet_nan)
  call fibre_prod_structure_commit_gate_evaluate(gate, nan_x, dry%u_trial, dry%dx_trial, status)
  call require(status == fibre_prod_commit_reject_nonfinite .and. gate%rejected, 'nonfinite trial not rejected')
  call fibre_prod_structure_commit_gate_commit_to_state(gate, state, nan_x, dry%u_trial, status)
  call require(all(state%x(1,:,:) == x0), 'nonfinite rejection changed coordinates')

  huge_dx = dry%dx_trial
  huge_dx(1,1) = 2.0e-2_dp
  call fibre_prod_structure_commit_gate_evaluate(gate, dry%x_trial, dry%u_trial, huge_dx, status)
  call require(status == fibre_prod_commit_reject_bound .and. gate%rejected, 'oversized trial not rejected')
  call fibre_prod_structure_commit_gate_commit_to_state(gate, state, dry%x_trial, dry%u_trial, status)
  call require(all(state%x(1,:,:) == x0), 'boundedness rejection changed coordinates')

  bad_x = 0.0_dp
  call fibre_prod_structure_commit_gate_evaluate(gate, bad_x, dry%u_trial, dry%dx_trial, status)
  call require(status == fibre_prod_commit_reject_shape .and. gate%rejected, 'shape mismatch not rejected')

  ux = reshape([(real(i, dp), i=1,8)], shape(ux))
  uy = ux + 10.0_dp
  uz = ux - 10.0_dp
  ux0 = ux; uy0 = uy; uz0 = uz
  rhs_x = ux * 0.25_dp; rhs_y = uy * 0.25_dp; rhs_z = uz * 0.25_dp
  rhs_x0 = rhs_x; rhs_y0 = rhs_y; rhs_z0 = rhs_z
  coord = [0.0_dp, 1.0_dp]
  call fibre_prod_grid_init_from_coordinates(grid, coord, coord, coord, 1, 2, 1, 2, 1, 2, .false., .false., .false., status)
  call require(status == 0, 'grid init failed')
  call fibre_prod_force_buffer_allocate(force_buffer, grid, status)
  call require(status == 0, 'force buffer allocation failed')
  call fibre_prod_structure_commit_gate_evaluate(gate, dry%x_trial, dry%u_trial, dry%dx_trial, status)
  call require(status == 0 .and. gate%accepted, 'lambda path gate acceptance failed')
  call fibre_prod_runtime_bridge_init_from_rhs(rhs_x, rhs_y, rhs_z, bridge, status)
  call require(status == 0, 'lambda0 runtime bridge init failed')
  call fibre_prod_runtime_config_default(config)
  config%enabled = .true.; config%lambda_fsi = 0.0_dp; config%penalty_beta = 2.0_dp
  call fibre_prod_main_hook_init(status, config)
  call require(status == 0, 'lambda0 main hook init failed')
  call fibre_prod_runtime_bridge_apply_lambda0_noop(rhs_x, rhs_y, rhs_z, bridge, status)
  call require(status == 0, 'lambda0 runtime bridge path failed')
  call require(all(ux == ux0) .and. all(uy == uy0) .and. all(uz == uz0), 'velocity fields changed')
  call require(all(rhs_x == rhs_x0) .and. all(rhs_y == rhs_y0) .and. all(rhs_z == rhs_z0), 'RHS fields changed')
  call require(all(force_buffer%fx == 0.0_dp) .and. all(force_buffer%fy == 0.0_dp) .and. &
               all(force_buffer%fz == 0.0_dp), 'force buffer changed')

  print *, 'P0_9_DIAGNOSTIC: structure dry-step commit only, no back-coupling payload, no RHS feedback'
  print *, 'P0_9_STRUCTURE_COMMIT_GATE_CHECK PASS'

  call fibre_prod_runtime_bridge_finalize(bridge)
  call fibre_prod_force_buffer_destroy(force_buffer)
  call fibre_prod_grid_destroy(grid)
  call fibre_prod_structure_commit_gate_finalize(gate)
  call fibre_prod_structure_dry_step_finalize(dry)
  call fibre_prod_structure_input_handoff_finalize(handoff)
  call fibre_prod_hydro_input_candidate_finalize(candidate)
  call fibre_prod_state_destroy(state)
contains
  subroutine setup_state(state, candidate, handoff, dry, x0, sampled_u0, hydro0, input0, state_u0, status)
    type(fibre_prod_state_type), intent(inout) :: state
    type(fibre_prod_hydro_input_candidate_type), intent(inout) :: candidate
    type(fibre_prod_structure_input_handoff_type), intent(inout) :: handoff
    type(fibre_prod_structure_dry_step_type), intent(inout) :: dry
    real(dp), intent(out) :: x0(3,3), sampled_u0(3,3), hydro0(3,3), input0(3,3), state_u0(3,3)
    integer, intent(out) :: status

    x0(1,:) = [0.25_dp, 0.5_dp, 0.5_dp]
    x0(2,:) = [0.50_dp, 0.5_dp, 0.5_dp]
    x0(3,:) = [0.75_dp, 0.5_dp, 0.5_dp]
    sampled_u0(1,:) = [1.0_dp, 0.5_dp, -0.25_dp]
    sampled_u0(2,:) = [1.5_dp, 0.25_dp, -0.50_dp]
    sampled_u0(3,:) = [2.0_dp, 0.0_dp, -0.75_dp]
    call fibre_prod_state_allocate(state, 1, 3, status)
    if (status /= 0) return
    state%x(1,:,:) = x0
    call fibre_prod_state_attach_sampled_velocity(state, sampled_u0, status)
    if (status /= 0) return
    call fibre_prod_state_get_structure_velocity_or_zero(state, state_u0, status)
    if (status /= 0) return
    call fibre_prod_hydro_input_candidate_init(candidate, 3, 2.5_dp, status)
    if (status /= 0) return
    call fibre_prod_hydro_input_candidate_compute(candidate, state%sampled_u, state_u0, status)
    if (status /= 0) return
    call fibre_prod_hydro_input_candidate_attach_to_state(candidate, state, status)
    if (status /= 0) return
    hydro0 = state%hydro_force_candidate
    call fibre_prod_structure_input_handoff_init(handoff, 3, status)
    if (status /= 0) return
    call fibre_prod_structure_input_handoff_from_candidate(handoff, state%hydro_force_candidate, status)
    if (status /= 0) return
    call fibre_prod_structure_input_handoff_attach_to_state(handoff, state, status)
    if (status /= 0) return
    input0 = state%structure_input_force
    call fibre_prod_state_get_structure_coordinates(state, x0, status)
    if (status /= 0) return
    call fibre_prod_structure_dry_step_init(dry, 3, 1.0e-4_dp, 2.0_dp, status)
    if (status /= 0) return
    call fibre_prod_structure_dry_step_predict(dry, x0, state_u0, state%structure_input_force, status)
  end subroutine setup_state

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      print *, 'P0_9_STRUCTURE_COMMIT_GATE_CHECK FAIL: ', trim(message)
      error stop 1
    end if
  end subroutine require
end program fibre_prod_structure_commit_gate_check
