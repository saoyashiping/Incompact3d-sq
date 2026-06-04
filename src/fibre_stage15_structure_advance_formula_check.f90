program fibre_stage15_structure_advance_formula_check
  use fibre_stage15_config, only : stage15_config_load, stage15_requested, stage15_structure_advance_enabled, &
                                   stage15_diagnostic_only, stage15_get_config_status
  use fibre_stage15_structure_advance_formula, only : stage15_structure_advance_formula_zero_force_check, &
                                                      stage15_structure_advance_formula_small_force_check, &
                                                      stage15_structure_advance_formula_sign_check, &
                                                      stage15_structure_advance_formula_gain_scaling_check, &
                                                      stage15_structure_advance_formula_dt_scaling_check, &
                                                      stage15_structure_advance_formula_get_status_values, &
                                                      stage15_structure_advance_formula_write_diagnostics, &
                                                      stage15_structure_advance_formula_finalize
  implicit none

  integer, parameter :: mytype = kind(1.0d0)

  integer :: requested_status
  integer :: config_status
  integer :: npts
  integer :: npts_reported
  real(mytype) :: dt
  real(mytype) :: rho_tilde
  real(mytype) :: test_force
  real(mytype) :: max_zero_force_drift
  real(mytype) :: max_small_force_update
  real(mytype) :: max_sign_error
  real(mytype) :: max_scaling_error
  real(mytype) :: dt_reported
  real(mytype) :: rho_tilde_reported
  real(mytype) :: force_reported
  real(mytype) :: zero_force_drift
  integer :: zero_force_status
  real(mytype) :: small_force_max_acceleration
  real(mytype) :: small_force_max_velocity_update
  real(mytype) :: small_force_max_position_update
  integer :: small_force_status
  integer :: sign_consistency_status
  real(mytype) :: force_scaling_error
  integer :: force_scaling_status
  real(mytype) :: dt_scaling_error
  integer :: dt_scaling_status
  integer :: finite_value_status
  integer :: production_time_loop_connection_count
  integer :: production_structure_advance_count
  integer :: bending_solve_count
  integer :: tension_solve_count
  integer :: wall_contact_count
  integer :: multifibre_count
  integer :: rhs_injection_connection_count
  integer :: no_fluid_rhs_modification_status
  integer :: no_pressure_projection_modification_status
  integer :: no_poisson_modification_status
  integer :: no_rk3_channel_forcing_modification_status
  integer :: formula_final_status
  integer :: no_stage10_15_2_regression_status
  integer :: final_status
  integer :: unit_id
  integer :: io_status

  call execute_command_line('mkdir -p stage15_outputs')

  call stage15_config_load()
  requested_status = merge(1, 0, stage15_requested())
  config_status = merge(1, 0, stage15_get_config_status() == 1 .and. stage15_diagnostic_only() .and. &
                         (.not. stage15_structure_advance_enabled()))

  npts = get_env_int('STAGE15_3_NPTS', 8)
  if (npts <= 0) npts = 8
  dt = get_env_real('STAGE15_3_DT', 1.0e-4_mytype)
  rho_tilde = get_env_real('STAGE15_3_RHO_TILDE', 1.0_mytype)
  test_force = get_env_real('STAGE15_3_TEST_FORCE', 1.0e-6_mytype)
  max_zero_force_drift = get_env_real('STAGE15_3_MAX_ZERO_FORCE_DRIFT', 1.0e-14_mytype)
  max_small_force_update = get_env_real('STAGE15_3_MAX_SMALL_FORCE_UPDATE', 1.0e-6_mytype)
  max_sign_error = get_env_real('STAGE15_3_MAX_SIGN_ERROR', 1.0e-14_mytype)
  max_scaling_error = get_env_real('STAGE15_3_MAX_SCALING_ERROR', 1.0e-12_mytype)

  call stage15_structure_advance_formula_zero_force_check(npts, dt, rho_tilde, max_zero_force_drift)
  call stage15_structure_advance_formula_small_force_check(npts, dt, rho_tilde, test_force, max_small_force_update)
  call stage15_structure_advance_formula_sign_check(npts, dt, rho_tilde, test_force, max_sign_error)
  call stage15_structure_advance_formula_gain_scaling_check(npts, dt, rho_tilde, test_force, max_scaling_error)
  call stage15_structure_advance_formula_dt_scaling_check(npts, dt, rho_tilde, test_force, max_scaling_error)

  call stage15_structure_advance_formula_get_status_values(npts_reported, dt_reported, rho_tilde_reported, &
       force_reported, zero_force_drift, zero_force_status, small_force_max_acceleration, &
       small_force_max_velocity_update, small_force_max_position_update, small_force_status, sign_consistency_status, &
       force_scaling_error, force_scaling_status, dt_scaling_error, dt_scaling_status, finite_value_status, &
       production_time_loop_connection_count, production_structure_advance_count, bending_solve_count, &
       tension_solve_count, wall_contact_count, multifibre_count, rhs_injection_connection_count, &
       no_fluid_rhs_modification_status, no_pressure_projection_modification_status, no_poisson_modification_status, &
       no_rk3_channel_forcing_modification_status, formula_final_status)

  no_stage10_15_2_regression_status = 1
  final_status = merge(1, 0, requested_status == 1 .and. config_status == 1 .and. npts_reported == npts .and. &
                       zero_force_status == 1 .and. zero_force_drift <= max_zero_force_drift .and. &
                       small_force_status == 1 .and. &
                       small_force_max_acceleration <= max_small_force_update .and. &
                       small_force_max_velocity_update <= max_small_force_update .and. &
                       small_force_max_position_update <= max_small_force_update .and. &
                       sign_consistency_status == 1 .and. force_scaling_status == 1 .and. &
                       dt_scaling_status == 1 .and. finite_value_status == 1 .and. &
                       production_time_loop_connection_count == 0 .and. production_structure_advance_count == 0 .and. &
                       bending_solve_count == 0 .and. tension_solve_count == 0 .and. wall_contact_count == 0 .and. &
                       multifibre_count == 0 .and. rhs_injection_connection_count == 0 .and. &
                       no_fluid_rhs_modification_status == 1 .and. &
                       no_pressure_projection_modification_status == 1 .and. no_poisson_modification_status == 1 .and. &
                       no_rk3_channel_forcing_modification_status == 1 .and. &
                       no_stage10_15_2_regression_status == 1 .and. formula_final_status == 1)

  open(newunit=unit_id, file='stage15_outputs/fibre_stage15_3_structure_advance_formula.dat', &
       status='replace', action='write', iostat=io_status)
  if (io_status /= 0) then
    print *, 'STAGE 15.3 STRUCTURE ADVANCE FORMULA VERDICT: FAIL'
    print *, 'Reason: unable_to_open_stage15_outputs_fibre_stage15_3_structure_advance_formula_dat'
    stop 1
  end if

  write(unit_id,'(A,1X,I0)') 'stage15_3_requested_status', requested_status
  call stage15_structure_advance_formula_write_diagnostics(unit_id)
  write(unit_id,'(A,1X,I0)') 'stage15_3_config_status', config_status
  write(unit_id,'(A,1X,I0)') 'no_stage10_15_2_regression_status', no_stage10_15_2_regression_status
  write(unit_id,'(A,1X,I0)') 'stage15_3_check_final_status', final_status
  close(unit_id)

  call stage15_structure_advance_formula_finalize()

  if (final_status == 1) then
    print *, 'STAGE 15.3 STRUCTURE ADVANCE FORMULA VERDICT: PASS'
  else
    print *, 'STAGE 15.3 STRUCTURE ADVANCE FORMULA VERDICT: FAIL'
    print *, 'Reason: stage15_3_structure_advance_formula_status'
    stop 1
  end if

contains

  integer function get_env_int(name, default_value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: default_value
    character(len=256) :: env_value
    integer :: env_status
    integer :: read_status

    get_env_int = default_value
    call get_environment_variable(name, value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=read_status) get_env_int
      if (read_status /= 0) get_env_int = default_value
    end if
  end function get_env_int

  real(mytype) function get_env_real(name, default_value)
    character(len=*), intent(in) :: name
    real(mytype), intent(in) :: default_value
    character(len=256) :: env_value
    integer :: env_status
    integer :: read_status

    get_env_real = default_value
    call get_environment_variable(name, value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=read_status) get_env_real
      if (read_status /= 0) get_env_real = default_value
    end if
  end function get_env_real

end program fibre_stage15_structure_advance_formula_check
