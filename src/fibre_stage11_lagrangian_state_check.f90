program fibre_stage11_lagrangian_state_check
  use fibre_stage11_config, only : stage11_config_load, stage11_requested, stage11_readonly_mode, stage11_get_max_points
  use fibre_stage11_lagrangian_state, only : stage11_lagrangian_state_init, stage11_lagrangian_state_finalize, &
                                             stage11_lagrangian_state_get_status_values, stage11_lagrangian_state_write_diagnostics
  implicit none

  integer :: requested_flag
  integer :: readonly_mode_status
  integer :: allocated_status
  integer :: point_count_status
  integer :: coordinates_finite_status
  integer :: velocity_placeholder_status
  integer :: valid_flag_status
  integer :: no_fluid_field_access_status
  integer :: no_velocity_sampling_status
  integer :: no_fluid_field_modification_status
  integer :: no_rhs_injection_status
  integer :: no_ibm_spreading_status
  integer :: no_feedback_force_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: lagrangian_state_status
  integer :: n_points
  integer :: io_unit
  integer :: ios
  logical :: pass

  call execute_command_line('mkdir -p stage11_outputs', exitstat=ios)

  call stage11_config_load()
  if (stage11_requested()) then
    requested_flag = 1
  else
    requested_flag = 0
  end if
  if (stage11_readonly_mode()) then
    readonly_mode_status = 1
  else
    readonly_mode_status = 0
  end if

  n_points = stage11_get_max_points()
  if (n_points <= 0) n_points = 8

  call stage11_lagrangian_state_init(n_points)
  call stage11_lagrangian_state_get_status_values(allocated_status, point_count_status, coordinates_finite_status, &
                                                   velocity_placeholder_status, valid_flag_status, &
                                                   no_fluid_field_access_status, no_velocity_sampling_status, &
                                                   no_fluid_field_modification_status, no_rhs_injection_status, &
                                                   no_ibm_spreading_status, no_feedback_force_status, &
                                                   no_twoway_force_status, no_structure_advance_status, lagrangian_state_status)

  open(newunit=io_unit, file='stage11_outputs/fibre_stage11_1_lagrangian_state.dat', status='replace', action='write', iostat=ios)
  if (ios /= 0) then
    print *, 'STAGE 11.1 LAGRANGIAN STATE VERDICT: FAIL'
    print *, 'Reason: unable to open stage11_outputs/fibre_stage11_1_lagrangian_state.dat'
    stop 1
  end if

  write(io_unit, '(A,1X,I0)') 'stage11_1_requested_flag', requested_flag
  write(io_unit, '(A,1X,I0)') 'stage11_1_readonly_mode_status', readonly_mode_status
  call stage11_lagrangian_state_write_diagnostics(io_unit)
  close(io_unit)

  pass = .true.
  if (requested_flag /= 1) then
    print *, 'Reason: stage11_1_requested_flag /= 1'
    pass = .false.
  end if
  if (readonly_mode_status /= 1) then
    print *, 'Reason: stage11_1_readonly_mode_status /= 1'
    pass = .false.
  end if
  if (allocated_status /= 1) then
    print *, 'Reason: stage11_1_allocated_status /= 1'
    pass = .false.
  end if
  if (point_count_status /= 1) then
    print *, 'Reason: stage11_1_point_count_status /= 1'
    pass = .false.
  end if
  if (coordinates_finite_status /= 1) then
    print *, 'Reason: stage11_1_coordinates_finite_status /= 1'
    pass = .false.
  end if
  if (velocity_placeholder_status /= 1) then
    print *, 'Reason: stage11_1_velocity_placeholder_status /= 1'
    pass = .false.
  end if
  if (valid_flag_status /= 1) then
    print *, 'Reason: stage11_1_valid_flag_status /= 1'
    pass = .false.
  end if
  if (no_fluid_field_access_status /= 1) then
    print *, 'Reason: stage11_1_no_fluid_field_access_status /= 1'
    pass = .false.
  end if
  if (no_velocity_sampling_status /= 1) then
    print *, 'Reason: stage11_1_no_velocity_sampling_status /= 1'
    pass = .false.
  end if
  if (no_fluid_field_modification_status /= 1) then
    print *, 'Reason: stage11_1_no_fluid_field_modification_status /= 1'
    pass = .false.
  end if
  if (no_rhs_injection_status /= 1) then
    print *, 'Reason: stage11_1_no_rhs_injection_status /= 1'
    pass = .false.
  end if
  if (no_ibm_spreading_status /= 1) then
    print *, 'Reason: stage11_1_no_ibm_spreading_status /= 1'
    pass = .false.
  end if
  if (no_feedback_force_status /= 1) then
    print *, 'Reason: stage11_1_no_feedback_force_status /= 1'
    pass = .false.
  end if
  if (no_twoway_force_status /= 1) then
    print *, 'Reason: stage11_1_no_twoway_force_status /= 1'
    pass = .false.
  end if
  if (no_structure_advance_status /= 1) then
    print *, 'Reason: stage11_1_no_structure_advance_status /= 1'
    pass = .false.
  end if
  if (lagrangian_state_status /= 1) then
    print *, 'Reason: stage11_1_lagrangian_state_status /= 1'
    pass = .false.
  end if

  if (pass) then
    print *, 'STAGE 11.1 LAGRANGIAN STATE VERDICT: PASS'
  else
    print *, 'STAGE 11.1 LAGRANGIAN STATE VERDICT: FAIL'
    stop 1
  end if

  call stage11_lagrangian_state_finalize()
end program fibre_stage11_lagrangian_state_check
