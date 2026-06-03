program fibre_stage15_config_check
  use fibre_stage15_config, only : stage15_config_load, &
                                   stage15_requested, &
                                   stage15_structure_advance_enabled, &
                                   stage15_diagnostic_only, &
                                   stage15_require_stage14_closed, &
                                   stage15_get_config_status, &
                                   stage15_get_status_values, &
                                   stage15_write_config_diagnostics
  implicit none

  integer :: requested_flag
  integer :: disabled_by_default_status
  integer :: request_parse_status
  integer :: structure_advance_requested_flag
  integer :: structure_advance_enabled_flag
  integer :: structure_advance_parse_status
  integer :: structure_advance_disabled_status
  integer :: structure_advance_blocked_status
  integer :: diagnostic_only_flag
  integer :: diagnostic_only_default_status
  integer :: diagnostic_only_parse_status
  integer :: diagnostic_only_enforced_status
  integer :: require_stage14_closed_flag
  integer :: require_stage14_closed_default_status
  integer :: require_stage14_closed_parse_status
  integer :: no_structure_state_allocation_status
  integer :: no_structure_advance_status
  integer :: no_bending_solve_status
  integer :: no_tension_solve_status
  integer :: no_fibre_position_update_status
  integer :: no_fibre_velocity_update_status
  integer :: no_wall_contact_status
  integer :: no_multifibre_logic_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: no_rhs_modification_status
  integer :: no_pressure_modification_status
  integer :: no_projection_modification_status
  integer :: no_poisson_modification_status
  integer :: no_rk3_modification_status
  integer :: no_channel_forcing_modification_status
  integer :: no_production_ibm_forcing_status
  integer :: no_stage11_14_reimplementation_status
  integer :: no_production_behavior_change_status
  integer :: config_status
  integer :: unit_id
  integer :: io_status
  integer :: final_status

  call execute_command_line('mkdir -p stage15_outputs')

  call stage15_config_load()
  call stage15_get_status_values(requested_flag, &
                                 disabled_by_default_status, &
                                 request_parse_status, &
                                 structure_advance_requested_flag, &
                                 structure_advance_enabled_flag, &
                                 structure_advance_parse_status, &
                                 structure_advance_disabled_status, &
                                 structure_advance_blocked_status, &
                                 diagnostic_only_flag, &
                                 diagnostic_only_default_status, &
                                 diagnostic_only_parse_status, &
                                 diagnostic_only_enforced_status, &
                                 require_stage14_closed_flag, &
                                 require_stage14_closed_default_status, &
                                 require_stage14_closed_parse_status, &
                                 no_structure_state_allocation_status, &
                                 no_structure_advance_status, &
                                 no_bending_solve_status, &
                                 no_tension_solve_status, &
                                 no_fibre_position_update_status, &
                                 no_fibre_velocity_update_status, &
                                 no_wall_contact_status, &
                                 no_multifibre_logic_status, &
                                 no_fluid_field_access_status, &
                                 no_fluid_field_modification_status, &
                                 no_rhs_modification_status, &
                                 no_pressure_modification_status, &
                                 no_projection_modification_status, &
                                 no_poisson_modification_status, &
                                 no_rk3_modification_status, &
                                 no_channel_forcing_modification_status, &
                                 no_production_ibm_forcing_status, &
                                 no_stage11_14_reimplementation_status, &
                                 no_production_behavior_change_status, &
                                 config_status)

  final_status = config_status
  if (stage15_structure_advance_enabled()) final_status = 0
  if (.not. stage15_diagnostic_only()) final_status = 0
  if (stage15_get_config_status() /= 1) final_status = 0

  open(newunit=unit_id, file='stage15_outputs/fibre_stage15_0_config.dat', &
       status='replace', action='write', iostat=io_status)
  if (io_status /= 0) then
    print *, 'STAGE 15.0 CONFIG VERDICT: FAIL'
    print *, 'Reason: unable_to_open_stage15_outputs_fibre_stage15_0_config_dat'
    stop 1
  end if

  call stage15_write_config_diagnostics(unit_id)
  write(unit_id, '(A,1X,I0)') 'stage15_0_requested_value', logical_to_int_local(stage15_requested())
  write(unit_id, '(A,1X,I0)') 'stage15_0_structure_advance_enabled_value', &
                               logical_to_int_local(stage15_structure_advance_enabled())
  write(unit_id, '(A,1X,I0)') 'stage15_0_diagnostic_only_value', logical_to_int_local(stage15_diagnostic_only())
  write(unit_id, '(A,1X,I0)') 'stage15_0_require_stage14_closed_value', &
                               logical_to_int_local(stage15_require_stage14_closed())
  write(unit_id, '(A,1X,I0)') 'stage15_0_check_status', final_status
  close(unit_id)

  if (final_status == 1) then
    print *, 'STAGE 15.0 CONFIG VERDICT: PASS'
  else
    print *, 'STAGE 15.0 CONFIG VERDICT: FAIL'
    print *, 'Reason: stage15_0_config_status'
    stop 1
  end if

contains

  integer function logical_to_int_local(value)
    logical, intent(in) :: value
    if (value) then
      logical_to_int_local = 1
    else
      logical_to_int_local = 0
    end if
  end function logical_to_int_local

end program fibre_stage15_config_check
