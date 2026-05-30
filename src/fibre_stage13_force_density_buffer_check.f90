program fibre_stage13_force_density_buffer_check
  use fibre_stage13_config, only : stage13_config_load, stage13_get_status_values
  use fibre_stage13_force_density_buffer, only : stage13_force_density_buffer_init, &
                                                stage13_force_density_buffer_clear, &
                                                stage13_force_density_buffer_finalize, &
                                                stage13_force_density_buffer_get_shape, &
                                                stage13_force_density_buffer_get_status_values
  implicit none

  integer :: requested_flag
  integer :: readonly_status
  integer :: spreading_readonly_status
  integer :: disabled_default_status
  integer :: max_points_status
  integer :: max_eulerian_points_status
  integer :: normalization_status
  integer :: config_no_force_density_allocation_status
  integer :: config_no_spreading_status
  integer :: config_no_rhs_injection_status
  integer :: config_no_ibm_spreading_status
  integer :: config_no_feedback_application_status
  integer :: config_no_twoway_force_status
  integer :: config_no_structure_advance_status
  integer :: config_no_fluid_field_access_status
  integer :: config_no_fluid_field_modification_status
  integer :: config_status
  integer :: allocated_status
  integer :: shape_status
  integer :: zero_force_density_status
  integer :: force_density_norm_finite_status
  integer :: force_density_valid_flag_status
  integer :: clear_status
  integer :: no_spreading_status
  integer :: no_rhs_injection_status
  integer :: no_ibm_spreading_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: force_density_buffer_status
  integer :: nx
  integer :: ny
  integer :: nz
  integer :: invalid_allocated_status
  integer :: invalid_shape_status
  integer :: unit_id
  integer :: io_status
  integer :: verdict_status
  real :: max_abs_fx
  real :: max_abs_fy
  real :: max_abs_fz
  real :: max_abs_force_density_norm_after_clear

  call stage13_config_load()
  call stage13_get_status_values(requested_flag, readonly_status, spreading_readonly_status, disabled_default_status, &
                                 max_points_status, max_eulerian_points_status, normalization_status, &
                                 config_no_force_density_allocation_status, config_no_spreading_status, &
                                 config_no_rhs_injection_status, config_no_ibm_spreading_status, &
                                 config_no_feedback_application_status, config_no_twoway_force_status, &
                                 config_no_structure_advance_status, config_no_fluid_field_access_status, &
                                 config_no_fluid_field_modification_status, config_status)

  call stage13_force_density_buffer_init(4, 0, 6)
  call stage13_force_density_buffer_get_status_values(invalid_allocated_status, invalid_shape_status, &
                                                      zero_force_density_status, force_density_norm_finite_status, &
                                                      force_density_valid_flag_status, clear_status, &
                                                      no_spreading_status, no_rhs_injection_status, &
                                                      no_ibm_spreading_status, no_feedback_application_status, &
                                                      no_twoway_force_status, no_structure_advance_status, &
                                                      no_fluid_field_access_status, no_fluid_field_modification_status, &
                                                      force_density_buffer_status)

  call stage13_force_density_buffer_init(4, 5, 6)
  call stage13_force_density_buffer_clear()
  call stage13_force_density_buffer_get_status_values(allocated_status, shape_status, zero_force_density_status, &
                                                      force_density_norm_finite_status, &
                                                      force_density_valid_flag_status, clear_status, &
                                                      no_spreading_status, no_rhs_injection_status, &
                                                      no_ibm_spreading_status, no_feedback_application_status, &
                                                      no_twoway_force_status, no_structure_advance_status, &
                                                      no_fluid_field_access_status, no_fluid_field_modification_status, &
                                                      force_density_buffer_status)
  call stage13_force_density_buffer_get_shape(nx, ny, nz)

  max_abs_fx = 0.0
  max_abs_fy = 0.0
  max_abs_fz = 0.0
  max_abs_force_density_norm_after_clear = 0.0

  open(newunit=unit_id, file='stage13_outputs/fibre_stage13_2_force_density_buffer.dat', status='replace', &
       action='write', iostat=io_status)
  if (io_status /= 0) then
    write(*,'(A)') 'STAGE 13.2 FORCE DENSITY BUFFER VERDICT: FAIL'
    write(*,'(A)') 'Reason: unable_to_open_stage13_outputs_fibre_stage13_2_force_density_buffer_dat'
    stop 1
  end if

  write(unit_id,'(A,1X,I0)') 'stage13_2_requested_flag', requested_flag
  write(unit_id,'(A,1X,I0)') 'stage13_2_readonly_mode_status', readonly_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_spreading_readonly_status', spreading_readonly_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_allocated_status', allocated_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_shape_status', shape_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_zero_force_density_status', zero_force_density_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_force_density_norm_finite_status', force_density_norm_finite_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_force_density_valid_flag_status', force_density_valid_flag_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_clear_status', clear_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_no_spreading_status', no_spreading_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_no_rhs_injection_status', no_rhs_injection_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_no_ibm_spreading_status', no_ibm_spreading_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_no_feedback_application_status', no_feedback_application_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_no_twoway_force_status', no_twoway_force_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_no_structure_advance_status', no_structure_advance_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_no_fluid_field_access_status', no_fluid_field_access_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_no_fluid_field_modification_status', no_fluid_field_modification_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_force_density_buffer_status', force_density_buffer_status
  write(unit_id,'(A,1X,I0)') 'stage13_2_nx', nx
  write(unit_id,'(A,1X,I0)') 'stage13_2_ny', ny
  write(unit_id,'(A,1X,I0)') 'stage13_2_nz', nz
  write(unit_id,'(A,1X,ES24.16)') 'stage13_2_max_abs_fx', max_abs_fx
  write(unit_id,'(A,1X,ES24.16)') 'stage13_2_max_abs_fy', max_abs_fy
  write(unit_id,'(A,1X,ES24.16)') 'stage13_2_max_abs_fz', max_abs_fz
  write(unit_id,'(A,1X,ES24.16)') 'stage13_2_max_abs_force_density_norm_after_clear', &
       max_abs_force_density_norm_after_clear
  close(unit_id)

  verdict_status = 1
  if (requested_flag /= 1) verdict_status = 0
  if (readonly_status /= 1) verdict_status = 0
  if (spreading_readonly_status /= 1) verdict_status = 0
  if (invalid_allocated_status /= 0) verdict_status = 0
  if (invalid_shape_status /= 0) verdict_status = 0
  if (allocated_status /= 1) verdict_status = 0
  if (shape_status /= 1) verdict_status = 0
  if (zero_force_density_status /= 1) verdict_status = 0
  if (force_density_norm_finite_status /= 1) verdict_status = 0
  if (force_density_valid_flag_status /= 1) verdict_status = 0
  if (clear_status /= 1) verdict_status = 0
  if (no_spreading_status /= 1) verdict_status = 0
  if (no_rhs_injection_status /= 1) verdict_status = 0
  if (no_ibm_spreading_status /= 1) verdict_status = 0
  if (no_feedback_application_status /= 1) verdict_status = 0
  if (no_twoway_force_status /= 1) verdict_status = 0
  if (no_structure_advance_status /= 1) verdict_status = 0
  if (no_fluid_field_access_status /= 1) verdict_status = 0
  if (no_fluid_field_modification_status /= 1) verdict_status = 0
  if (force_density_buffer_status /= 1) verdict_status = 0

  call stage13_force_density_buffer_finalize()

  if (verdict_status == 1) then
    write(*,'(A)') 'STAGE 13.2 FORCE DENSITY BUFFER VERDICT: PASS'
  else
    write(*,'(A)') 'STAGE 13.2 FORCE DENSITY BUFFER VERDICT: FAIL'
    write(*,'(A)') 'Reason: stage13_2_force_density_buffer_status_failure'
    stop 1
  end if
end program fibre_stage13_force_density_buffer_check
