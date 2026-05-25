program fibre_stage11_config_check
  use fibre_stage11_config, only : stage11_config_load, stage11_get_status_values
  implicit none

  integer :: requested_flag, readonly_mode_status, disabled_by_default_status
  integer :: no_lagrangian_state_status, no_velocity_sampling_status
  integer :: no_fluid_field_access_status, no_fluid_field_modification_status
  integer :: no_rhs_injection_status, no_ibm_spreading_status
  integer :: no_feedback_force_status, no_twoway_force_status
  integer :: no_structure_advance_status, config_status
  integer :: verdict, unit_dat

  call stage11_config_load()
  call stage11_get_status_values(requested_flag, readonly_mode_status, disabled_by_default_status, &
       no_lagrangian_state_status, no_velocity_sampling_status, no_fluid_field_access_status, &
       no_fluid_field_modification_status, no_rhs_injection_status, no_ibm_spreading_status, &
       no_feedback_force_status, no_twoway_force_status, no_structure_advance_status, config_status)

  call execute_command_line('mkdir -p stage11_outputs')
  open(newunit=unit_dat, file='stage11_outputs/fibre_stage11_0_config.dat', status='replace', action='write')
  write(unit_dat,*) 'stage11_0_requested_flag ', requested_flag
  write(unit_dat,*) 'stage11_0_readonly_mode_status ', readonly_mode_status
  write(unit_dat,*) 'stage11_0_disabled_by_default_status ', disabled_by_default_status
  write(unit_dat,*) 'stage11_0_no_lagrangian_state_status ', no_lagrangian_state_status
  write(unit_dat,*) 'stage11_0_no_velocity_sampling_status ', no_velocity_sampling_status
  write(unit_dat,*) 'stage11_0_no_fluid_field_access_status ', no_fluid_field_access_status
  write(unit_dat,*) 'stage11_0_no_fluid_field_modification_status ', no_fluid_field_modification_status
  write(unit_dat,*) 'stage11_0_no_rhs_injection_status ', no_rhs_injection_status
  write(unit_dat,*) 'stage11_0_no_ibm_spreading_status ', no_ibm_spreading_status
  write(unit_dat,*) 'stage11_0_no_feedback_force_status ', no_feedback_force_status
  write(unit_dat,*) 'stage11_0_no_twoway_force_status ', no_twoway_force_status
  write(unit_dat,*) 'stage11_0_no_structure_advance_status ', no_structure_advance_status
  write(unit_dat,*) 'stage11_0_config_status ', config_status
  close(unit_dat)

  verdict = 1
  if (requested_flag /= 1) verdict = 0
  if (readonly_mode_status /= 1) verdict = 0
  if (no_lagrangian_state_status /= 1) verdict = 0
  if (no_velocity_sampling_status /= 1) verdict = 0
  if (no_fluid_field_access_status /= 1) verdict = 0
  if (no_fluid_field_modification_status /= 1) verdict = 0
  if (no_rhs_injection_status /= 1) verdict = 0
  if (no_ibm_spreading_status /= 1) verdict = 0
  if (no_feedback_force_status /= 1) verdict = 0
  if (no_twoway_force_status /= 1) verdict = 0
  if (no_structure_advance_status /= 1) verdict = 0
  if (config_status /= 1) verdict = 0

  if (verdict == 1) then
    write(*,'(A)') 'STAGE 11.0 CONFIG VERDICT: PASS'
  else
    write(*,'(A)') 'STAGE 11.0 CONFIG VERDICT: FAIL'
  endif
end program fibre_stage11_config_check
