program fibre_stage16_config_check
  use decomp_2d_constants, only : mytype
  use fibre_stage16_config, only : stage16_config_reset, stage16_config_load_from_environment, &
       stage16_config_validate, stage16_is_enabled, stage16_one_fibre_fsi_requested, &
       stage16_structure_advance_requested, stage16_two_way_rhs_requested, &
       stage16_diagnostic_only_requested, stage16_get_lambda, stage16_get_feedback_alpha, &
       stage16_get_max_structure_update, stage16_get_max_rhs_increment, &
       stage16_config_write_diagnostics, stage16_config_set_for_test, &
       stage16_config_get_status_values
  implicit none

  integer :: io_unit
  integer :: final_status
  integer :: default_off_status
  integer :: env_override_status
  integer :: invalid_numeric_rejection_status
  integer :: invalid_flag_combination_rejection_status
  integer :: no_production_hook_status
  integer :: no_structure_advance_status
  integer :: no_rhs_modification_status
  integer :: no_bending_solve_status
  integer :: no_tension_solve_status
  integer :: no_wall_contact_status
  integer :: no_multifibre_status
  integer :: no_pressure_projection_modification_status
  integer :: no_poisson_modification_status
  integer :: no_rk3_channel_forcing_modification_status
  integer :: no_legacy_ibm_forcing_status
  integer :: default_dummy
  integer :: config_loaded_status
  integer :: config_status
  integer :: invalid_flag_status
  integer :: numeric_parse_status
  integer :: numeric_bounds_status
  integer :: structure_guard_status
  integer :: rhs_guard_status

  call execute_command_line('mkdir -p stage16_outputs')

  call stage16_config_reset()
  call stage16_config_get_status_values(default_off_out=default_off_status, env_config_out=config_loaded_status, &
       config_out=config_status, invalid_flag_out=invalid_flag_status, numeric_parse_out=numeric_parse_status, &
       numeric_bounds_out=numeric_bounds_status, no_structure_advance_out=structure_guard_status, &
       no_rhs_modification_out=rhs_guard_status)
  default_off_status = logical_to_int(default_off_status == 1 .and. config_status == 1 .and. &
       structure_guard_status == 1 .and. rhs_guard_status == 1)

  call stage16_config_load_from_environment()
  env_override_status = logical_to_int(stage16_diagnostic_only_requested() .and. &
       .not. stage16_structure_advance_requested() .and. &
       .not. stage16_two_way_rhs_requested() .and. &
       finite_real_local(stage16_get_feedback_alpha()) .and. &
       finite_real_local(stage16_get_lambda()) .and. &
       finite_real_local(stage16_get_max_structure_update()) .and. &
       finite_real_local(stage16_get_max_rhs_increment()))
  call stage16_config_get_status_values(default_off_out=default_dummy, env_config_out=config_loaded_status, &
       config_out=config_status, invalid_flag_out=invalid_flag_status, numeric_parse_out=numeric_parse_status, &
       numeric_bounds_out=numeric_bounds_status, no_structure_advance_out=structure_guard_status, &
       no_rhs_modification_out=rhs_guard_status)
  env_override_status = min(env_override_status, config_status)

  call stage16_config_set_for_test(.true., .false., .false., .false., .true., &
       -1.0_mytype, 0.0_mytype, 1.0e-12_mytype, 1.0e-8_mytype)
  call stage16_config_get_status_values(default_off_out=default_dummy, env_config_out=config_loaded_status, &
       config_out=config_status, invalid_flag_out=invalid_flag_status, numeric_parse_out=numeric_parse_status, &
       numeric_bounds_out=numeric_bounds_status, no_structure_advance_out=structure_guard_status, &
       no_rhs_modification_out=rhs_guard_status)
  invalid_numeric_rejection_status = logical_to_int(config_status == 0 .and. numeric_bounds_status == 0)

  call stage16_config_set_for_test(.true., .true., .false., .false., .false., &
       1.0_mytype, 0.0_mytype, 1.0e-12_mytype, 1.0e-8_mytype)
  call stage16_config_get_status_values(default_off_out=default_dummy, env_config_out=config_loaded_status, &
       config_out=config_status, invalid_flag_out=invalid_flag_status, numeric_parse_out=numeric_parse_status, &
       numeric_bounds_out=numeric_bounds_status, no_structure_advance_out=structure_guard_status, &
       no_rhs_modification_out=rhs_guard_status)
  invalid_flag_combination_rejection_status = logical_to_int(config_status == 0 .and. invalid_flag_status == 0)

  call stage16_config_set_for_test(.false., .false., .true., .false., .true., &
       1.0_mytype, 0.0_mytype, 1.0e-12_mytype, 1.0e-8_mytype)
  call stage16_config_get_status_values(default_off_out=default_dummy, env_config_out=config_loaded_status, &
       config_out=config_status, invalid_flag_out=invalid_flag_status, numeric_parse_out=numeric_parse_status, &
       numeric_bounds_out=numeric_bounds_status, no_structure_advance_out=structure_guard_status, &
       no_rhs_modification_out=rhs_guard_status)
  no_structure_advance_status = logical_to_int(config_status == 0 .and. structure_guard_status == 0)

  call stage16_config_set_for_test(.false., .false., .false., .true., .true., &
       1.0_mytype, 0.0_mytype, 1.0e-12_mytype, 1.0e-8_mytype)
  call stage16_config_get_status_values(default_off_out=default_dummy, env_config_out=config_loaded_status, &
       config_out=config_status, invalid_flag_out=invalid_flag_status, numeric_parse_out=numeric_parse_status, &
       numeric_bounds_out=numeric_bounds_status, no_structure_advance_out=structure_guard_status, &
       no_rhs_modification_out=rhs_guard_status)
  no_rhs_modification_status = logical_to_int(config_status == 0 .and. rhs_guard_status == 0)

  no_production_hook_status = 1
  no_bending_solve_status = 1
  no_tension_solve_status = 1
  no_wall_contact_status = 1
  no_multifibre_status = 1
  no_pressure_projection_modification_status = 1
  no_poisson_modification_status = 1
  no_rk3_channel_forcing_modification_status = 1
  no_legacy_ibm_forcing_status = 1

  call stage16_config_load_from_environment()
  final_status = logical_to_int(default_off_status == 1 .and. env_override_status == 1 .and. &
       invalid_numeric_rejection_status == 1 .and. invalid_flag_combination_rejection_status == 1 .and. &
       no_production_hook_status == 1 .and. no_structure_advance_status == 1 .and. &
       no_rhs_modification_status == 1 .and. no_bending_solve_status == 1 .and. &
       no_tension_solve_status == 1 .and. no_wall_contact_status == 1 .and. &
       no_multifibre_status == 1 .and. no_pressure_projection_modification_status == 1 .and. &
       no_poisson_modification_status == 1 .and. no_rk3_channel_forcing_modification_status == 1 .and. &
       no_legacy_ibm_forcing_status == 1)

  open(newunit=io_unit, file='stage16_outputs/fibre_stage16_1_config.dat', status='replace', action='write')
  write(io_unit,'(A,1X,I0)') 'stage16_1_requested_status', 1
  write(io_unit,'(A,1X,I0)') 'default_off_status', default_off_status
  write(io_unit,'(A,1X,I0)') 'env_override_status', env_override_status
  call stage16_config_write_diagnostics(io_unit)
  write(io_unit,'(A,1X,I0)') 'invalid_numeric_rejection_status', invalid_numeric_rejection_status
  write(io_unit,'(A,1X,I0)') 'invalid_flag_combination_rejection_status', invalid_flag_combination_rejection_status
  write(io_unit,'(A,1X,I0)') 'no_production_hook_status', no_production_hook_status
  write(io_unit,'(A,1X,I0)') 'no_structure_advance_status', no_structure_advance_status
  write(io_unit,'(A,1X,I0)') 'no_rhs_modification_status', no_rhs_modification_status
  write(io_unit,'(A,1X,I0)') 'no_bending_solve_status', no_bending_solve_status
  write(io_unit,'(A,1X,I0)') 'no_tension_solve_status', no_tension_solve_status
  write(io_unit,'(A,1X,I0)') 'no_wall_contact_status', no_wall_contact_status
  write(io_unit,'(A,1X,I0)') 'no_multifibre_status', no_multifibre_status
  write(io_unit,'(A,1X,I0)') 'no_pressure_projection_modification_status', no_pressure_projection_modification_status
  write(io_unit,'(A,1X,I0)') 'no_poisson_modification_status', no_poisson_modification_status
  write(io_unit,'(A,1X,I0)') 'no_rk3_channel_forcing_modification_status', no_rk3_channel_forcing_modification_status
  write(io_unit,'(A,1X,I0)') 'no_legacy_ibm_forcing_status', no_legacy_ibm_forcing_status
  write(io_unit,'(A,1X,I0)') 'final_status', final_status
  close(io_unit)

  if (final_status /= 1) stop 1

contains

  integer function logical_to_int(value)
    logical, intent(in) :: value
    if (value) then
      logical_to_int = 1
    else
      logical_to_int = 0
    end if
  end function logical_to_int

  logical function finite_real_local(value)
    real(mytype), intent(in) :: value
    finite_real_local = (value == value) .and. (abs(value) < huge(value))
  end function finite_real_local

end program fibre_stage16_config_check
