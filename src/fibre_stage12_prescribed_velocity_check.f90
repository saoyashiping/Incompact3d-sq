program fibre_stage12_prescribed_velocity_check
  use fibre_stage12_config, only: stage12_config_load, stage12_get_max_points, stage12_get_status_values
  use fibre_stage12_prescribed_velocity, only: stage12_prescribed_velocity_init, &
                                               stage12_prescribed_velocity_set_constant, &
                                               stage12_prescribed_velocity_clear, &
                                               stage12_prescribed_velocity_finalize, &
                                               stage12_prescribed_velocity_get_max_abs_velocity, &
                                               stage12_prescribed_velocity_get_status_values
  implicit none

  integer, parameter :: mytype = kind(1.0d0)

  integer :: requested_flag
  integer :: readonly_mode_status
  integer :: disabled_by_default_status
  integer :: feedback_gain_status
  integer :: force_sign_status
  integer :: config_no_force_computation_status
  integer :: config_no_eulerian_force_density_status
  integer :: config_no_rhs_injection_status
  integer :: config_no_ibm_spreading_status
  integer :: config_no_feedback_application_status
  integer :: config_no_twoway_force_status
  integer :: config_no_structure_advance_status
  integer :: config_no_fluid_field_access_status
  integer :: config_no_fluid_field_modification_status
  integer :: config_status
  integer :: allocated_status
  integer :: point_count_status
  integer :: zero_velocity_status
  integer :: constant_velocity_status
  integer :: velocity_norm_finite_status
  integer :: velocity_valid_flag_status
  integer :: clear_status
  integer :: no_force_computation_status
  integer :: no_eulerian_force_density_status
  integer :: no_rhs_injection_status
  integer :: no_ibm_spreading_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: prescribed_velocity_status
  integer :: n_points
  integer :: io_unit
  integer :: ios
  logical :: pass
  real(mytype) :: vx_const
  real(mytype) :: vy_const
  real(mytype) :: vz_const
  real(mytype) :: constant_velocity_norm
  real(mytype) :: max_abs_velocity_after_clear

  call execute_command_line('mkdir -p stage12_outputs', exitstat=ios)

  call stage12_config_load()
  call stage12_get_status_values(requested_flag, readonly_mode_status, disabled_by_default_status, &
                                 feedback_gain_status, force_sign_status, config_no_force_computation_status, &
                                 config_no_eulerian_force_density_status, config_no_rhs_injection_status, &
                                 config_no_ibm_spreading_status, config_no_feedback_application_status, &
                                 config_no_twoway_force_status, config_no_structure_advance_status, &
                                 config_no_fluid_field_access_status, config_no_fluid_field_modification_status, &
                                 config_status)

  n_points = stage12_get_max_points()
  if (n_points <= 0) n_points = 8

  vx_const = 0.25_mytype
  vy_const = -0.10_mytype
  vz_const = 0.05_mytype
  constant_velocity_norm = sqrt(vx_const * vx_const + vy_const * vy_const + vz_const * vz_const)

  call stage12_prescribed_velocity_init(n_points)
  call stage12_prescribed_velocity_get_status_values(allocated_status, point_count_status, zero_velocity_status, &
                                                     constant_velocity_status, velocity_norm_finite_status, &
                                                     velocity_valid_flag_status, clear_status, &
                                                     no_force_computation_status, no_eulerian_force_density_status, &
                                                     no_rhs_injection_status, no_ibm_spreading_status, &
                                                     no_feedback_application_status, no_twoway_force_status, &
                                                     no_structure_advance_status, no_fluid_field_access_status, &
                                                     no_fluid_field_modification_status, prescribed_velocity_status)

  call stage12_prescribed_velocity_set_constant(vx_const, vy_const, vz_const)
  call stage12_prescribed_velocity_get_status_values(allocated_status, point_count_status, zero_velocity_status, &
                                                     constant_velocity_status, velocity_norm_finite_status, &
                                                     velocity_valid_flag_status, clear_status, &
                                                     no_force_computation_status, no_eulerian_force_density_status, &
                                                     no_rhs_injection_status, no_ibm_spreading_status, &
                                                     no_feedback_application_status, no_twoway_force_status, &
                                                     no_structure_advance_status, no_fluid_field_access_status, &
                                                     no_fluid_field_modification_status, prescribed_velocity_status)

  call stage12_prescribed_velocity_clear()
  call stage12_prescribed_velocity_get_status_values(allocated_status, point_count_status, zero_velocity_status, &
                                                     constant_velocity_status, velocity_norm_finite_status, &
                                                     velocity_valid_flag_status, clear_status, &
                                                     no_force_computation_status, no_eulerian_force_density_status, &
                                                     no_rhs_injection_status, no_ibm_spreading_status, &
                                                     no_feedback_application_status, no_twoway_force_status, &
                                                     no_structure_advance_status, no_fluid_field_access_status, &
                                                     no_fluid_field_modification_status, prescribed_velocity_status)
  max_abs_velocity_after_clear = stage12_prescribed_velocity_get_max_abs_velocity()

  open(newunit=io_unit, file='stage12_outputs/fibre_stage12_2_prescribed_velocity.dat', status='replace', &
       action='write', iostat=ios)
  if (ios /= 0) then
    print *, 'STAGE 12.2 PRESCRIBED VELOCITY VERDICT: FAIL'
    print *, 'Reason: unable to open stage12_outputs/fibre_stage12_2_prescribed_velocity.dat'
    call stage12_prescribed_velocity_finalize()
    stop 1
  end if

  write(io_unit, '(A,1X,I0)') 'stage12_2_requested_flag', requested_flag
  write(io_unit, '(A,1X,I0)') 'stage12_2_readonly_mode_status', readonly_mode_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_allocated_status', allocated_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_point_count_status', point_count_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_zero_velocity_status', zero_velocity_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_constant_velocity_status', constant_velocity_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_velocity_norm_finite_status', velocity_norm_finite_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_velocity_valid_flag_status', velocity_valid_flag_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_clear_status', clear_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_no_force_computation_status', no_force_computation_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_no_eulerian_force_density_status', no_eulerian_force_density_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_no_rhs_injection_status', no_rhs_injection_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_no_ibm_spreading_status', no_ibm_spreading_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_no_feedback_application_status', no_feedback_application_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_no_twoway_force_status', no_twoway_force_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_no_structure_advance_status', no_structure_advance_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_no_fluid_field_access_status', no_fluid_field_access_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_no_fluid_field_modification_status', no_fluid_field_modification_status
  write(io_unit, '(A,1X,I0)') 'stage12_2_prescribed_velocity_status', prescribed_velocity_status
  write(io_unit, '(A,1X,ES24.16)') 'stage12_2_constant_vx', vx_const
  write(io_unit, '(A,1X,ES24.16)') 'stage12_2_constant_vy', vy_const
  write(io_unit, '(A,1X,ES24.16)') 'stage12_2_constant_vz', vz_const
  write(io_unit, '(A,1X,ES24.16)') 'stage12_2_constant_velocity_norm', constant_velocity_norm
  write(io_unit, '(A,1X,ES24.16)') 'stage12_2_max_abs_velocity_after_clear', max_abs_velocity_after_clear
  close(io_unit)

  pass = .true.
  call require_status('stage12_2_requested_flag', requested_flag, 1, pass)
  call require_status('stage12_2_readonly_mode_status', readonly_mode_status, 1, pass)
  call require_status('stage12_2_allocated_status', allocated_status, 1, pass)
  call require_status('stage12_2_point_count_status', point_count_status, 1, pass)
  call require_status('stage12_2_zero_velocity_status', zero_velocity_status, 1, pass)
  call require_status('stage12_2_constant_velocity_status', constant_velocity_status, 1, pass)
  call require_status('stage12_2_velocity_norm_finite_status', velocity_norm_finite_status, 1, pass)
  call require_status('stage12_2_velocity_valid_flag_status', velocity_valid_flag_status, 1, pass)
  call require_status('stage12_2_clear_status', clear_status, 1, pass)
  call require_status('stage12_2_no_force_computation_status', no_force_computation_status, 1, pass)
  call require_status('stage12_2_no_eulerian_force_density_status', no_eulerian_force_density_status, 1, pass)
  call require_status('stage12_2_no_rhs_injection_status', no_rhs_injection_status, 1, pass)
  call require_status('stage12_2_no_ibm_spreading_status', no_ibm_spreading_status, 1, pass)
  call require_status('stage12_2_no_feedback_application_status', no_feedback_application_status, 1, pass)
  call require_status('stage12_2_no_twoway_force_status', no_twoway_force_status, 1, pass)
  call require_status('stage12_2_no_structure_advance_status', no_structure_advance_status, 1, pass)
  call require_status('stage12_2_no_fluid_field_access_status', no_fluid_field_access_status, 1, pass)
  call require_status('stage12_2_no_fluid_field_modification_status', no_fluid_field_modification_status, 1, pass)
  call require_status('stage12_2_prescribed_velocity_status', prescribed_velocity_status, 1, pass)

  if (abs(max_abs_velocity_after_clear) > 0.0_mytype) then
    print *, 'Reason: stage12_2_max_abs_velocity_after_clear /= 0'
    pass = .false.
  end if

  call stage12_prescribed_velocity_finalize()

  if (pass) then
    print *, 'STAGE 12.2 PRESCRIBED VELOCITY VERDICT: PASS'
  else
    print *, 'STAGE 12.2 PRESCRIBED VELOCITY VERDICT: FAIL'
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

end program fibre_stage12_prescribed_velocity_check
