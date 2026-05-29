program fibre_stage12_config_check
  use fibre_stage12_config, only: stage12_config_load, stage12_get_status_values
  implicit none

  integer :: requested_flag
  integer :: readonly_mode_status
  integer :: disabled_by_default_status
  integer :: feedback_gain_status
  integer :: force_sign_status
  integer :: no_force_computation_status
  integer :: no_eulerian_force_density_status
  integer :: no_rhs_injection_status
  integer :: no_ibm_spreading_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: config_status
  integer :: io_unit
  integer :: ios
  logical :: pass

  call execute_command_line('mkdir -p stage12_outputs', exitstat=ios)

  call stage12_config_load()
  call stage12_get_status_values(requested_flag, readonly_mode_status, disabled_by_default_status, &
                                 feedback_gain_status, force_sign_status, no_force_computation_status, &
                                 no_eulerian_force_density_status, no_rhs_injection_status, &
                                 no_ibm_spreading_status, no_feedback_application_status, &
                                 no_twoway_force_status, no_structure_advance_status, &
                                 no_fluid_field_access_status, no_fluid_field_modification_status, &
                                 config_status)

  open(newunit=io_unit, file='stage12_outputs/fibre_stage12_0_config.dat', status='replace', action='write', iostat=ios)
  if (ios /= 0) then
    print *, 'STAGE 12.0 CONFIG VERDICT: FAIL'
    print *, 'Reason: unable to open stage12_outputs/fibre_stage12_0_config.dat'
    stop 1
  end if

  write(io_unit, '(A,1X,I0)') 'stage12_0_requested_flag', requested_flag
  write(io_unit, '(A,1X,I0)') 'stage12_0_readonly_mode_status', readonly_mode_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_disabled_by_default_status', disabled_by_default_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_feedback_gain_status', feedback_gain_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_force_sign_status', force_sign_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_no_force_computation_status', no_force_computation_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_no_eulerian_force_density_status', no_eulerian_force_density_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_no_rhs_injection_status', no_rhs_injection_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_no_ibm_spreading_status', no_ibm_spreading_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_no_feedback_application_status', no_feedback_application_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_no_twoway_force_status', no_twoway_force_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_no_structure_advance_status', no_structure_advance_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_no_fluid_field_access_status', no_fluid_field_access_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_no_fluid_field_modification_status', no_fluid_field_modification_status
  write(io_unit, '(A,1X,I0)') 'stage12_0_config_status', config_status
  close(io_unit)

  pass = .true.
  call require_status('stage12_0_requested_flag', requested_flag, 1, pass)
  call require_status('stage12_0_readonly_mode_status', readonly_mode_status, 1, pass)
  call require_status('stage12_0_feedback_gain_status', feedback_gain_status, 1, pass)
  call require_status('stage12_0_force_sign_status', force_sign_status, 1, pass)
  call require_status('stage12_0_no_force_computation_status', no_force_computation_status, 1, pass)
  call require_status('stage12_0_no_eulerian_force_density_status', no_eulerian_force_density_status, 1, pass)
  call require_status('stage12_0_no_rhs_injection_status', no_rhs_injection_status, 1, pass)
  call require_status('stage12_0_no_ibm_spreading_status', no_ibm_spreading_status, 1, pass)
  call require_status('stage12_0_no_feedback_application_status', no_feedback_application_status, 1, pass)
  call require_status('stage12_0_no_twoway_force_status', no_twoway_force_status, 1, pass)
  call require_status('stage12_0_no_structure_advance_status', no_structure_advance_status, 1, pass)
  call require_status('stage12_0_no_fluid_field_access_status', no_fluid_field_access_status, 1, pass)
  call require_status('stage12_0_no_fluid_field_modification_status', no_fluid_field_modification_status, 1, pass)
  call require_status('stage12_0_config_status', config_status, 1, pass)

  if (pass) then
    print *, 'STAGE 12.0 CONFIG VERDICT: PASS'
  else
    print *, 'STAGE 12.0 CONFIG VERDICT: FAIL'
    stop 1
  end if

contains

  subroutine require_status(name, actual_value, expected_value, pass)
    character(len=*), intent(in) :: name
    integer, intent(in) :: actual_value
    integer, intent(in) :: expected_value
    logical, intent(inout) :: pass

    if (actual_value /= expected_value) then
      print *, 'Reason: ', trim(name), ' /= ', expected_value
      pass = .false.
    end if
  end subroutine require_status

end program fibre_stage12_config_check
