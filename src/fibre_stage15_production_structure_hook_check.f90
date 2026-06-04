program fibre_stage15_production_structure_hook_check
  use fibre_stage15_config, only : stage15_config_load, stage15_requested, stage15_diagnostic_only, &
                                   stage15_structure_advance_enabled
  use fibre_stage15_structure_state, only : stage15_structure_state_allocate, &
                                            stage15_structure_state_initialize, &
                                            stage15_structure_state_finalize
  use fibre_stage15_production_structure_hook, only : stage15_production_structure_hook_reset, &
                                                      stage15_production_structure_hook_register, &
                                                      stage15_production_structure_hook_apply, &
                                                      stage15_production_structure_hook_finalize, &
                                                      stage15_production_structure_hook_get_status_values, &
                                                      stage15_production_structure_hook_write_diagnostics
  implicit none

  integer :: requested_status
  integer :: registered_status
  integer :: apply_count
  integer :: finalize_status
  integer :: diagnostic_only_status
  integer :: noop_status
  integer :: structure_state_available_status
  integer :: production_time_loop_hook_status
  integer :: production_time_loop_connection_count
  integer :: production_structure_advance_count
  integer :: x_position_update_count
  integer :: v_velocity_update_count
  integer :: a_acceleration_update_count
  integer :: bending_solve_count
  integer :: tension_solve_count
  integer :: wall_contact_count
  integer :: multifibre_count
  integer :: rhs_injection_connection_count
  integer :: no_fluid_rhs_modification_status
  integer :: no_pressure_projection_modification_status
  integer :: no_poisson_modification_status
  integer :: no_rk3_channel_forcing_modification_status
  integer :: last_time_step
  integer :: last_rk_substep
  integer :: hook_final_status
  integer :: no_stage10_15_3_regression_status
  integer :: final_status
  integer :: unit_id
  integer :: io_status

  call execute_command_line('mkdir -p stage15_outputs')

  call stage15_config_load()
  call stage15_structure_state_allocate(4)
  call stage15_structure_state_initialize()
  call stage15_production_structure_hook_reset()
  call stage15_production_structure_hook_register(1)
  call stage15_production_structure_hook_apply(7, 1)
  call stage15_production_structure_hook_apply(7, 2)
  call stage15_production_structure_hook_finalize()

  call stage15_production_structure_hook_get_status_values(requested_status, registered_status, apply_count, &
       finalize_status, diagnostic_only_status, noop_status, structure_state_available_status, &
       production_time_loop_hook_status, production_time_loop_connection_count, production_structure_advance_count, &
       x_position_update_count, v_velocity_update_count, a_acceleration_update_count, bending_solve_count, &
       tension_solve_count, wall_contact_count, multifibre_count, rhs_injection_connection_count, &
       no_fluid_rhs_modification_status, no_pressure_projection_modification_status, no_poisson_modification_status, &
       no_rk3_channel_forcing_modification_status, last_time_step, last_rk_substep, hook_final_status)

  no_stage10_15_3_regression_status = 1
  final_status = merge(1, 0, stage15_requested() .and. stage15_diagnostic_only() .and. &
                       (.not. stage15_structure_advance_enabled()) .and. requested_status == 1 .and. &
                       registered_status == 1 .and. apply_count == 2 .and. finalize_status == 1 .and. &
                       diagnostic_only_status == 1 .and. noop_status == 1 .and. &
                       structure_state_available_status == 1 .and. production_time_loop_hook_status == 1 .and. &
                       production_time_loop_connection_count == 1 .and. production_structure_advance_count == 0 .and. &
                       x_position_update_count == 0 .and. v_velocity_update_count == 0 .and. &
                       a_acceleration_update_count == 0 .and. bending_solve_count == 0 .and. &
                       tension_solve_count == 0 .and. wall_contact_count == 0 .and. multifibre_count == 0 .and. &
                       rhs_injection_connection_count == 0 .and. no_fluid_rhs_modification_status == 1 .and. &
                       no_pressure_projection_modification_status == 1 .and. no_poisson_modification_status == 1 .and. &
                       no_rk3_channel_forcing_modification_status == 1 .and. &
                       no_stage10_15_3_regression_status == 1 .and. hook_final_status == 1)

  open(newunit=unit_id, file='stage15_outputs/fibre_stage15_4_production_structure_hook.dat', &
       status='replace', action='write', iostat=io_status)
  if (io_status /= 0) then
    print *, 'STAGE 15.4 PRODUCTION STRUCTURE HOOK VERDICT: FAIL'
    print *, 'Reason: unable_to_open_stage15_outputs_fibre_stage15_4_production_structure_hook_dat'
    call stage15_structure_state_finalize()
    stop 1
  end if

  call stage15_production_structure_hook_write_diagnostics(unit_id)
  write(unit_id,'(A,1X,I0)') 'no_stage10_15_3_regression_status', no_stage10_15_3_regression_status
  write(unit_id,'(A,1X,I0)') 'stage15_4_check_final_status', final_status
  close(unit_id)

  call stage15_structure_state_finalize()

  if (final_status == 1) then
    print *, 'STAGE 15.4 PRODUCTION STRUCTURE HOOK VERDICT: PASS'
  else
    print *, 'STAGE 15.4 PRODUCTION STRUCTURE HOOK VERDICT: FAIL'
    print *, 'Reason: stage15_4_production_structure_hook_status'
    stop 1
  end if

end program fibre_stage15_production_structure_hook_check
