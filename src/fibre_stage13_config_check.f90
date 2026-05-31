program fibre_stage13_config_check
  use fibre_stage13_config, only : stage13_config_load, stage13_get_max_points, &
                                  stage13_get_max_eulerian_points, stage13_get_normalization_mode, &
                                  stage13_get_status_values
  implicit none

  integer :: requested_flag
  integer :: readonly_status
  integer :: spreading_readonly_status
  integer :: disabled_default_status
  integer :: max_points_status
  integer :: max_eulerian_points_status
  integer :: normalization_status
  integer :: no_force_density_allocation_status
  integer :: no_spreading_status
  integer :: no_rhs_injection_status
  integer :: no_ibm_spreading_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: config_status
  integer :: max_points
  integer :: max_eulerian_points
  integer :: normalization_mode
  integer :: unit_id
  integer :: io_status
  integer :: verdict_status

  call stage13_config_load()
  call stage13_get_status_values(requested_flag, readonly_status, spreading_readonly_status, disabled_default_status, &
                                 max_points_status, max_eulerian_points_status, normalization_status, &
                                 no_force_density_allocation_status, no_spreading_status, no_rhs_injection_status, &
                                 no_ibm_spreading_status, no_feedback_application_status, no_twoway_force_status, &
                                 no_structure_advance_status, no_fluid_field_access_status, &
                                 no_fluid_field_modification_status, config_status)
  max_points = stage13_get_max_points()
  max_eulerian_points = stage13_get_max_eulerian_points()
  normalization_mode = stage13_get_normalization_mode()

  open(newunit=unit_id, file='stage13_outputs/fibre_stage13_1_config.dat', status='replace', action='write', &
       iostat=io_status)
  if (io_status /= 0) then
    write(*,'(A)') 'STAGE 13.1 CONFIG VERDICT: FAIL'
    write(*,'(A)') 'Reason: unable_to_open_stage13_outputs_fibre_stage13_1_config_dat'
    stop 1
  end if

  write(unit_id,'(A,1X,I0)') 'stage13_1_requested_flag', requested_flag
  write(unit_id,'(A,1X,I0)') 'stage13_1_readonly_mode_status', readonly_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_spreading_readonly_status', spreading_readonly_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_disabled_by_default_status', disabled_default_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_max_points_status', max_points_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_max_eulerian_points_status', max_eulerian_points_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_normalization_mode_status', normalization_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_no_force_density_allocation_status', no_force_density_allocation_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_no_spreading_status', no_spreading_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_no_rhs_injection_status', no_rhs_injection_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_no_ibm_spreading_status', no_ibm_spreading_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_no_feedback_application_status', no_feedback_application_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_no_twoway_force_status', no_twoway_force_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_no_structure_advance_status', no_structure_advance_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_no_fluid_field_access_status', no_fluid_field_access_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_no_fluid_field_modification_status', no_fluid_field_modification_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_config_status', config_status
  write(unit_id,'(A,1X,I0)') 'stage13_1_max_points', max_points
  write(unit_id,'(A,1X,I0)') 'stage13_1_max_eulerian_points', max_eulerian_points
  write(unit_id,'(A,1X,I0)') 'stage13_1_normalization_mode_code', normalization_mode
  close(unit_id)

  verdict_status = 1
  if (readonly_status /= 1) verdict_status = 0
  if (spreading_readonly_status /= 1) verdict_status = 0
  if (max_points_status /= 1) verdict_status = 0
  if (max_eulerian_points_status /= 1) verdict_status = 0
  if (normalization_status /= 1) verdict_status = 0
  if (no_force_density_allocation_status /= 1) verdict_status = 0
  if (no_spreading_status /= 1) verdict_status = 0
  if (no_rhs_injection_status /= 1) verdict_status = 0
  if (no_ibm_spreading_status /= 1) verdict_status = 0
  if (no_feedback_application_status /= 1) verdict_status = 0
  if (no_twoway_force_status /= 1) verdict_status = 0
  if (no_structure_advance_status /= 1) verdict_status = 0
  if (no_fluid_field_access_status /= 1) verdict_status = 0
  if (no_fluid_field_modification_status /= 1) verdict_status = 0
  if (config_status /= 1) verdict_status = 0

  if (verdict_status == 1) then
    write(*,'(A)') 'STAGE 13.1 CONFIG VERDICT: PASS'
  else
    write(*,'(A)') 'STAGE 13.1 CONFIG VERDICT: FAIL'
    write(*,'(A)') 'Reason: stage13_1_config_status_failure'
    stop 1
  end if
end program fibre_stage13_config_check
