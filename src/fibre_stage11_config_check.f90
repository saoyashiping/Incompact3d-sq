program fibre_stage11_config_check
  use fibre_stage11_config, only: stage11_config_load, stage11_get_status_values
  implicit none

  integer :: requested_flag
  integer :: readonly_mode_status
  integer :: disabled_by_default_status
  integer :: no_lagrangian_state_status
  integer :: no_velocity_sampling_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: no_rhs_injection_status
  integer :: no_ibm_spreading_status
  integer :: no_feedback_force_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: config_status
  integer :: io_unit
  integer :: ios
  logical :: pass

  call execute_command_line('mkdir -p stage11_outputs', exitstat=ios)

  call stage11_config_load()
  call stage11_get_status_values(requested_flag, readonly_mode_status, disabled_by_default_status, &
                                 no_lagrangian_state_status, no_velocity_sampling_status, &
                                 no_fluid_field_access_status, no_fluid_field_modification_status, &
                                 no_rhs_injection_status, no_ibm_spreading_status, no_feedback_force_status, &
                                 no_twoway_force_status, no_structure_advance_status, config_status)

  open(newunit=io_unit, file='stage11_outputs/fibre_stage11_0_config.dat', status='replace', action='write', iostat=ios)
  if (ios /= 0) then
    print *, 'STAGE 11.0 CONFIG VERDICT: FAIL'
    print *, 'Reason: unable to open stage11_outputs/fibre_stage11_0_config.dat'
    stop 1
  end if

  write(io_unit, '(A,1X,I0)') 'stage11_0_requested_flag', requested_flag
  write(io_unit, '(A,1X,I0)') 'stage11_0_readonly_mode_status', readonly_mode_status
  write(io_unit, '(A,1X,I0)') 'stage11_0_disabled_by_default_status', disabled_by_default_status
  write(io_unit, '(A,1X,I0)') 'stage11_0_no_lagrangian_state_status', no_lagrangian_state_status
  write(io_unit, '(A,1X,I0)') 'stage11_0_no_velocity_sampling_status', no_velocity_sampling_status
  write(io_unit, '(A,1X,I0)') 'stage11_0_no_fluid_field_access_status', no_fluid_field_access_status
  write(io_unit, '(A,1X,I0)') 'stage11_0_no_fluid_field_modification_status', no_fluid_field_modification_status
  write(io_unit, '(A,1X,I0)') 'stage11_0_no_rhs_injection_status', no_rhs_injection_status
  write(io_unit, '(A,1X,I0)') 'stage11_0_no_ibm_spreading_status', no_ibm_spreading_status
  write(io_unit, '(A,1X,I0)') 'stage11_0_no_feedback_force_status', no_feedback_force_status
  write(io_unit, '(A,1X,I0)') 'stage11_0_no_twoway_force_status', no_twoway_force_status
  write(io_unit, '(A,1X,I0)') 'stage11_0_no_structure_advance_status', no_structure_advance_status
  write(io_unit, '(A,1X,I0)') 'stage11_0_config_status', config_status
  close(io_unit)

  pass = .true.
  if (requested_flag /= 1) then
    print *, 'Reason: stage11_0_requested_flag /= 1'
    pass = .false.
  end if
  if (readonly_mode_status /= 1) then
    print *, 'Reason: stage11_0_readonly_mode_status /= 1'
    pass = .false.
  end if
  if (no_lagrangian_state_status /= 1) then
    print *, 'Reason: stage11_0_no_lagrangian_state_status /= 1'
    pass = .false.
  end if
  if (no_velocity_sampling_status /= 1) then
    print *, 'Reason: stage11_0_no_velocity_sampling_status /= 1'
    pass = .false.
  end if
  if (no_fluid_field_access_status /= 1) then
    print *, 'Reason: stage11_0_no_fluid_field_access_status /= 1'
    pass = .false.
  end if
  if (no_fluid_field_modification_status /= 1) then
    print *, 'Reason: stage11_0_no_fluid_field_modification_status /= 1'
    pass = .false.
  end if
  if (no_rhs_injection_status /= 1) then
    print *, 'Reason: stage11_0_no_rhs_injection_status /= 1'
    pass = .false.
  end if
  if (no_ibm_spreading_status /= 1) then
    print *, 'Reason: stage11_0_no_ibm_spreading_status /= 1'
    pass = .false.
  end if
  if (no_feedback_force_status /= 1) then
    print *, 'Reason: stage11_0_no_feedback_force_status /= 1'
    pass = .false.
  end if
  if (no_twoway_force_status /= 1) then
    print *, 'Reason: stage11_0_no_twoway_force_status /= 1'
    pass = .false.
  end if
  if (no_structure_advance_status /= 1) then
    print *, 'Reason: stage11_0_no_structure_advance_status /= 1'
    pass = .false.
  end if
  if (config_status /= 1) then
    print *, 'Reason: stage11_0_config_status /= 1'
    pass = .false.
  end if

  if (pass) then
    print *, 'STAGE 11.0 CONFIG VERDICT: PASS'
  else
    print *, 'STAGE 11.0 CONFIG VERDICT: FAIL'
    stop 1
  end if

end program fibre_stage11_config_check
