program fibre_stage15_controlled_structure_step_check
  use fibre_stage15_controlled_structure_step, only : &
       stage15_controlled_structure_step_apply, &
       stage15_controlled_structure_step_get_status_values, &
       stage15_controlled_structure_step_write_diagnostics
  implicit none

  integer, parameter :: mytype = kind(1.0d0)
  integer :: npts
  integer :: np
  integer :: stage15_enable
  integer :: controlled_enable
  integer :: structure_enable
  integer :: diagnostic_only
  real(mytype) :: dt
  real(mytype) :: rho_tilde
  real(mytype) :: test_force
  real(mytype) :: max_position_update_allowed
  real(mytype) :: max_velocity_update_allowed
  real(mytype) :: max_acceleration_allowed
  real(mytype) :: max_zero_force_drift_allowed
  integer :: requested_status
  integer :: controlled_status
  integer :: structure_status
  integer :: diagnostic_status
  integer :: np_status
  integer :: npts_status
  real(mytype) :: dt_status
  real(mytype) :: rho_status
  real(mytype) :: force_status
  real(mytype) :: zero_force_drift
  integer :: zero_force_status
  integer :: small_force_status
  real(mytype) :: max_acceleration
  real(mytype) :: max_velocity_update
  real(mytype) :: max_position_update
  integer :: forced_component_nonzero_status
  integer :: sign_consistency_status
  integer :: bounded_update_status
  integer :: controlled_update_count
  integer :: production_full_structure_advance_count
  integer :: bending_solve_count
  integer :: tension_solve_count
  integer :: wall_contact_count
  integer :: multifibre_count
  integer :: rhs_injection_connection_count
  integer :: no_fluid_rhs_modification_status
  integer :: no_pressure_projection_modification_status
  integer :: no_poisson_modification_status
  integer :: no_rk3_channel_forcing_modification_status
  integer :: final_status
  integer :: check_status
  integer :: unit_id

  npts = get_env_int('STAGE15_6_NPTS', 8)
  np = get_env_int('STAGE15_6_NP', 1)
  stage15_enable = get_env_int('STAGE15_6_ENABLE', 1)
  controlled_enable = get_env_int('STAGE15_6_CONTROLLED_STEP_ENABLE', 1)
  structure_enable = get_env_int('STAGE15_6_STRUCTURE_ADVANCE_ENABLE', 1)
  diagnostic_only = get_env_int('STAGE15_6_DIAGNOSTIC_ONLY', 1)
  dt = get_env_real('STAGE15_6_DT', 1.0e-4_mytype)
  rho_tilde = get_env_real('STAGE15_6_RHO_TILDE', 1.0_mytype)
  test_force = get_env_real('STAGE15_6_TEST_FORCE', 1.0e-6_mytype)
  max_position_update_allowed = get_env_real('STAGE15_6_MAX_POSITION_UPDATE', 1.0e-12_mytype)
  max_velocity_update_allowed = get_env_real('STAGE15_6_MAX_VELOCITY_UPDATE', 1.0e-9_mytype)
  max_acceleration_allowed = get_env_real('STAGE15_6_MAX_ACCELERATION', 1.0e-5_mytype)
  max_zero_force_drift_allowed = get_env_real('STAGE15_6_MAX_ZERO_FORCE_DRIFT', 1.0e-14_mytype)

  call execute_command_line('mkdir -p stage15_outputs')
  call stage15_controlled_structure_step_apply(npts, np, dt, rho_tilde, test_force, &
                                               max_position_update_allowed, max_velocity_update_allowed, &
                                               max_acceleration_allowed, max_zero_force_drift_allowed, &
                                               stage15_enable, controlled_enable, structure_enable, diagnostic_only)
  call stage15_controlled_structure_step_get_status_values(requested_status, controlled_status, structure_status, &
       diagnostic_status, np_status, npts_status, dt_status, rho_status, force_status, zero_force_drift, &
       zero_force_status, small_force_status, max_acceleration, max_velocity_update, max_position_update, &
       forced_component_nonzero_status, sign_consistency_status, bounded_update_status, controlled_update_count, &
       production_full_structure_advance_count, bending_solve_count, tension_solve_count, wall_contact_count, &
       multifibre_count, rhs_injection_connection_count, no_fluid_rhs_modification_status, &
       no_pressure_projection_modification_status, no_poisson_modification_status, &
       no_rk3_channel_forcing_modification_status, final_status)

  open(newunit=unit_id, file='stage15_outputs/fibre_stage15_6_controlled_structure_step_np1.dat', &
       status='replace', action='write')
  call stage15_controlled_structure_step_write_diagnostics(unit_id)
  write(unit_id, '(A,1X,I0)') 'stage15_6_check_final_status', final_status
  close(unit_id)

  check_status = 1
  if (requested_status /= 1) check_status = 0
  if (controlled_status /= 1) check_status = 0
  if (structure_status /= 1) check_status = 0
  if (diagnostic_status /= 1) check_status = 0
  if (np_status /= 1) check_status = 0
  if (npts_status /= npts) check_status = 0
  if (zero_force_status /= 1) check_status = 0
  if (small_force_status /= 1) check_status = 0
  if (forced_component_nonzero_status /= 1) check_status = 0
  if (sign_consistency_status /= 1) check_status = 0
  if (bounded_update_status /= 1) check_status = 0
  if (controlled_update_count /= 1) check_status = 0
  if (production_full_structure_advance_count /= 0) check_status = 0
  if (bending_solve_count /= 0) check_status = 0
  if (tension_solve_count /= 0) check_status = 0
  if (wall_contact_count /= 0) check_status = 0
  if (multifibre_count /= 0) check_status = 0
  if (rhs_injection_connection_count /= 0) check_status = 0
  if (no_fluid_rhs_modification_status /= 1) check_status = 0
  if (no_pressure_projection_modification_status /= 1) check_status = 0
  if (no_poisson_modification_status /= 1) check_status = 0
  if (no_rk3_channel_forcing_modification_status /= 1) check_status = 0
  if (final_status /= 1) check_status = 0

  if (check_status == 1) then
    write(*, '(A)') 'STAGE 15.6 CONTROLLED STRUCTURE STEP CHECK: PASS'
    stop 0
  end if

  write(*, '(A)') 'STAGE 15.6 CONTROLLED STRUCTURE STEP CHECK: FAIL'
  stop 1

contains

  integer function get_env_int(name, default_value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: default_value
    character(len=128) :: value
    integer :: env_status
    integer :: read_status

    get_env_int = default_value
    call get_environment_variable(name, value=value, status=env_status)
    if (env_status == 0) then
      read(value, *, iostat=read_status) get_env_int
      if (read_status /= 0) get_env_int = default_value
    end if
  end function get_env_int

  real(mytype) function get_env_real(name, default_value)
    character(len=*), intent(in) :: name
    real(mytype), intent(in) :: default_value
    character(len=128) :: value
    integer :: env_status
    integer :: read_status

    get_env_real = default_value
    call get_environment_variable(name, value=value, status=env_status)
    if (env_status == 0) then
      read(value, *, iostat=read_status) get_env_real
      if (read_status /= 0) get_env_real = default_value
    end if
  end function get_env_real

end program fibre_stage15_controlled_structure_step_check
