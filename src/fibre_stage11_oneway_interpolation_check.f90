program fibre_stage11_oneway_interpolation_check
  use fibre_stage11_config, only : stage11_config_load, stage11_requested, stage11_readonly_mode, stage11_get_max_points
  use fibre_stage11_lagrangian_state, only : stage11_lagrangian_state_init, stage11_lagrangian_state_finalize, stage11_lagrangian_state_is_allocated
  use fibre_stage11_grid_metadata, only : stage11_grid_metadata_init_uniform, stage11_grid_metadata_finalize, stage11_grid_metadata_is_initialized
  use fibre_stage11_oneway_interpolation, only : stage11_oneway_interpolation_init, stage11_oneway_interpolation_finalize, &
                                                 stage11_oneway_interpolation_is_initialized, stage11_oneway_interpolation_prepare, &
                                                 stage11_oneway_interpolation_sample_velocity, stage11_oneway_interpolation_get_status_values, &
                                                 stage11_oneway_interpolation_write_diagnostics
  implicit none
  integer :: requested_flag, readonly_mode_status
  integer :: lagrangian_state_available_status, grid_metadata_available_status
  integer :: interp_initialized_status, interface_available_status, prepare_called_status
  integer :: sample_interface_called_status, sample_not_performed_status
  integer :: lagrangian_state_input_status, grid_metadata_input_status
  integer :: velocity_placeholder_unchanged_status
  integer :: no_fluid_field_access_status, no_velocity_sampling_status, no_fluid_field_modification_status
  integer :: no_rhs_injection_status, no_ibm_spreading_status, no_feedback_force_status
  integer :: no_twoway_force_status, no_structure_advance_status, oneway_interpolation_api_status
  integer :: n_points, io_unit, ios
  logical :: pass

  call execute_command_line('mkdir -p stage11_outputs', exitstat=ios)
  call stage11_config_load()
  requested_flag = 0; if (stage11_requested()) requested_flag = 1
  readonly_mode_status = 0; if (stage11_readonly_mode()) readonly_mode_status = 1
  n_points = stage11_get_max_points(); if (n_points <= 0) n_points = 8

  call stage11_lagrangian_state_init(n_points)
  lagrangian_state_available_status = 0
  if (stage11_lagrangian_state_is_allocated()) lagrangian_state_available_status = 1

  call stage11_grid_metadata_init_uniform()
  grid_metadata_available_status = 0
  if (stage11_grid_metadata_is_initialized()) grid_metadata_available_status = 1

  call stage11_oneway_interpolation_init()
  call stage11_oneway_interpolation_prepare(.true., .true.)
  call stage11_oneway_interpolation_sample_velocity()

  call stage11_oneway_interpolation_get_status_values(interp_initialized_status, interface_available_status, prepare_called_status, &
                                                      sample_interface_called_status, sample_not_performed_status, &
                                                      lagrangian_state_input_status, grid_metadata_input_status, &
                                                      velocity_placeholder_unchanged_status, no_fluid_field_access_status, &
                                                      no_velocity_sampling_status, no_fluid_field_modification_status, &
                                                      no_rhs_injection_status, no_ibm_spreading_status, no_feedback_force_status, &
                                                      no_twoway_force_status, no_structure_advance_status, oneway_interpolation_api_status)

  open(newunit=io_unit,file='stage11_outputs/fibre_stage11_3_oneway_interpolation_api.dat',status='replace',action='write',iostat=ios)
  if (ios /= 0) then
    print *, 'STAGE 11.3 ONEWAY INTERPOLATION API VERDICT: FAIL'
    print *, 'Reason: unable to open stage11_outputs/fibre_stage11_3_oneway_interpolation_api.dat'
    stop 1
  end if

  write(io_unit,'(A,1X,I0)') 'stage11_3_requested_flag', requested_flag
  write(io_unit,'(A,1X,I0)') 'stage11_3_readonly_mode_status', readonly_mode_status
  write(io_unit,'(A,1X,I0)') 'stage11_3_lagrangian_state_available_status', lagrangian_state_available_status
  write(io_unit,'(A,1X,I0)') 'stage11_3_grid_metadata_available_status', grid_metadata_available_status
  call stage11_oneway_interpolation_write_diagnostics(io_unit)
  close(io_unit)

  pass = .true.
  if (requested_flag /= 1) then; print *, 'Reason: stage11_3_requested_flag /= 1'; pass=.false.; end if
  if (readonly_mode_status /= 1) then; print *, 'Reason: stage11_3_readonly_mode_status /= 1'; pass=.false.; end if
  if (lagrangian_state_available_status /= 1) then; print *, 'Reason: stage11_3_lagrangian_state_available_status /= 1'; pass=.false.; end if
  if (grid_metadata_available_status /= 1) then; print *, 'Reason: stage11_3_grid_metadata_available_status /= 1'; pass=.false.; end if
  if (interp_initialized_status /= 1) then; print *, 'Reason: stage11_3_interpolation_initialized_status /= 1'; pass=.false.; end if
  if (interface_available_status /= 1) then; print *, 'Reason: stage11_3_interface_available_status /= 1'; pass=.false.; end if
  if (prepare_called_status /= 1) then; print *, 'Reason: stage11_3_prepare_called_status /= 1'; pass=.false.; end if
  if (sample_interface_called_status /= 1) then; print *, 'Reason: stage11_3_sample_interface_called_status /= 1'; pass=.false.; end if
  if (sample_not_performed_status /= 1) then; print *, 'Reason: stage11_3_sample_not_performed_status /= 1'; pass=.false.; end if
  if (lagrangian_state_input_status /= 1) then; print *, 'Reason: stage11_3_lagrangian_state_input_status /= 1'; pass=.false.; end if
  if (grid_metadata_input_status /= 1) then; print *, 'Reason: stage11_3_grid_metadata_input_status /= 1'; pass=.false.; end if
  if (velocity_placeholder_unchanged_status /= 1) then; print *, 'Reason: stage11_3_velocity_placeholder_unchanged_status /= 1'; pass=.false.; end if
  if (no_fluid_field_access_status /= 1) then; print *, 'Reason: stage11_3_no_fluid_field_access_status /= 1'; pass=.false.; end if
  if (no_velocity_sampling_status /= 1) then; print *, 'Reason: stage11_3_no_velocity_sampling_status /= 1'; pass=.false.; end if
  if (no_fluid_field_modification_status /= 1) then; print *, 'Reason: stage11_3_no_fluid_field_modification_status /= 1'; pass=.false.; end if
  if (no_rhs_injection_status /= 1) then; print *, 'Reason: stage11_3_no_rhs_injection_status /= 1'; pass=.false.; end if
  if (no_ibm_spreading_status /= 1) then; print *, 'Reason: stage11_3_no_ibm_spreading_status /= 1'; pass=.false.; end if
  if (no_feedback_force_status /= 1) then; print *, 'Reason: stage11_3_no_feedback_force_status /= 1'; pass=.false.; end if
  if (no_twoway_force_status /= 1) then; print *, 'Reason: stage11_3_no_twoway_force_status /= 1'; pass=.false.; end if
  if (no_structure_advance_status /= 1) then; print *, 'Reason: stage11_3_no_structure_advance_status /= 1'; pass=.false.; end if
  if (oneway_interpolation_api_status /= 1) then; print *, 'Reason: stage11_3_oneway_interpolation_api_status /= 1'; pass=.false.; end if

  if (pass) then
    print *, 'STAGE 11.3 ONEWAY INTERPOLATION API VERDICT: PASS'
  else
    print *, 'STAGE 11.3 ONEWAY INTERPOLATION API VERDICT: FAIL'
    stop 1
  end if

  call stage11_oneway_interpolation_finalize()
  call stage11_grid_metadata_finalize()
  call stage11_lagrangian_state_finalize()
end program fibre_stage11_oneway_interpolation_check
