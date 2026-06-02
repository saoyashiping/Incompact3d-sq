program fibre_stage14_config_check
  use decomp_2d_constants, only : mytype
  use fibre_stage14_config
  implicit none

  integer :: requested_flag
  integer :: rhs_injection_enabled_flag
  integer :: disabled_by_default_status
  integer :: gain_parse_status
  integer :: gain_default_zero_status
  integer :: gain_finite_status
  integer :: safe_fallback_status
  integer :: max_steps_parse_status
  integer :: require_stage13_status
  integer :: require_stage13_parse_status
  integer :: diagnostic_only_status
  integer :: diagnostic_only_forced_status
  integer :: no_rhs_buffer_allocation_status
  integer :: no_rhs_injection_status
  integer :: no_production_hook_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: no_pressure_modification_status
  integer :: no_projection_modification_status
  integer :: no_rk3_modification_status
  integer :: no_channel_forcing_modification_status
  integer :: no_production_ibm_forcing_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: config_status
  integer :: unit_id
  integer :: io_status
  integer :: final_status
  real(mytype) :: injection_gain
  integer :: max_steps
  integer :: require_stage13_int
  integer :: diagnostic_only_int

  call execute_command_line('mkdir -p stage14_outputs')

  call stage14_config_load()
  call stage14_get_status_values(requested_flag, &
                                 rhs_injection_enabled_flag, &
                                 disabled_by_default_status, &
                                 gain_parse_status, &
                                 gain_default_zero_status, &
                                 gain_finite_status, &
                                 safe_fallback_status, &
                                 max_steps_parse_status, &
                                 require_stage13_status, &
                                 require_stage13_parse_status, &
                                 diagnostic_only_status, &
                                 diagnostic_only_forced_status, &
                                 no_rhs_buffer_allocation_status, &
                                 no_rhs_injection_status, &
                                 no_production_hook_status, &
                                 no_fluid_field_access_status, &
                                 no_fluid_field_modification_status, &
                                 no_pressure_modification_status, &
                                 no_projection_modification_status, &
                                 no_rk3_modification_status, &
                                 no_channel_forcing_modification_status, &
                                 no_production_ibm_forcing_status, &
                                 no_feedback_application_status, &
                                 no_twoway_force_status, &
                                 no_structure_advance_status, &
                                 config_status)

  injection_gain = stage14_get_injection_gain()
  max_steps = stage14_get_max_steps()
  require_stage13_int = logical_to_int_local(stage14_require_stage13())
  diagnostic_only_int = 1
  final_status = config_status

  open(newunit=unit_id, file='stage14_outputs/fibre_stage14_0_config.dat', &
       status='replace', action='write', iostat=io_status)
  if (io_status /= 0) then
    final_status = 0
    print *, 'STAGE 14.0 CONFIG VERDICT: FAIL'
    print *, 'Reason: unable_to_open_stage14_outputs_fibre_stage14_0_config_dat'
    stop 1
  end if

  write(unit_id, '(A,1X,I0)') 'stage14_0_requested_flag', requested_flag
  write(unit_id, '(A,1X,I0)') 'stage14_0_rhs_injection_enabled_flag', rhs_injection_enabled_flag
  write(unit_id, '(A,1X,I0)') 'stage14_0_disabled_by_default_status', disabled_by_default_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_gain_parse_status', gain_parse_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_gain_default_zero_status', gain_default_zero_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_gain_finite_status', gain_finite_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_safe_fallback_status', safe_fallback_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_max_steps_parse_status', max_steps_parse_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_require_stage13_status', require_stage13_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_require_stage13_parse_status', require_stage13_parse_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_diagnostic_only_status', diagnostic_only_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_diagnostic_only_forced_status', diagnostic_only_forced_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_rhs_buffer_allocation_status', no_rhs_buffer_allocation_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_rhs_injection_status', no_rhs_injection_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_production_hook_status', no_production_hook_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_fluid_field_access_status', no_fluid_field_access_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_fluid_field_modification_status', no_fluid_field_modification_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_pressure_modification_status', no_pressure_modification_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_projection_modification_status', no_projection_modification_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_rk3_modification_status', no_rk3_modification_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_channel_forcing_modification_status', no_channel_forcing_modification_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_production_ibm_forcing_status', no_production_ibm_forcing_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_feedback_application_status', no_feedback_application_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_twoway_force_status', no_twoway_force_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_no_structure_advance_status', no_structure_advance_status
  write(unit_id, '(A,1X,I0)') 'stage14_0_config_status', config_status
  write(unit_id, '(A,1X,ES24.16)') 'stage14_0_injection_gain', injection_gain
  write(unit_id, '(A,1X,I0)') 'stage14_0_max_steps', max_steps
  write(unit_id, '(A,1X,I0)') 'stage14_0_require_stage13_value', require_stage13_int
  write(unit_id, '(A,1X,I0)') 'stage14_0_diagnostic_only_value', diagnostic_only_int
  close(unit_id)

  if (final_status == 1) then
    print *, 'STAGE 14.0 CONFIG VERDICT: PASS'
  else
    print *, 'STAGE 14.0 CONFIG VERDICT: FAIL'
    print *, 'Reason: stage14_0_config_status'
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

end program fibre_stage14_config_check
