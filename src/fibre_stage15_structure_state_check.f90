program fibre_stage15_structure_state_check
  use fibre_stage15_config, only : stage15_config_load, stage15_requested, stage15_structure_advance_enabled, &
                                   stage15_diagnostic_only, stage15_get_config_status
  use fibre_stage15_structure_state, only : stage15_structure_state_allocate, &
                                            stage15_structure_state_initialize, &
                                            stage15_structure_state_clear, &
                                            stage15_structure_state_validate, &
                                            stage15_structure_state_write_diagnostics, &
                                            stage15_structure_state_finalize, &
                                            stage15_structure_state_get_status_values
  implicit none

  integer :: requested_status
  integer :: npts
  integer :: allocation_status
  integer :: initialization_status
  integer :: clear_status
  integer :: validation_status
  integer :: npts_reported
  integer :: x_finite_status
  integer :: v_finite_status
  integer :: optional_a_or_rhs_finite_status
  integer :: optional_tension_finite_status
  integer :: structure_advance_count
  integer :: bending_solve_count
  integer :: tension_solve_count
  integer :: position_time_update_count
  integer :: velocity_time_update_count
  integer :: no_fluid_rhs_modification_status
  integer :: no_pressure_projection_modification_status
  integer :: no_poisson_modification_status
  integer :: no_rk3_channel_forcing_modification_status
  integer :: diagnostic_write_count
  integer :: diagnostic_counter_status
  integer :: state_final_status
  integer :: no_structure_advance_status
  integer :: no_bending_solve_status
  integer :: no_tension_solve_status
  integer :: no_position_time_update_status
  integer :: no_velocity_time_update_status
  integer :: no_stage11_14_regression_status
  integer :: no_fluid_rhs_or_solver_modification_status
  integer :: final_status
  integer :: unit_id
  integer :: io_status

  call execute_command_line('mkdir -p stage15_outputs')

  call stage15_config_load()
  requested_status = merge(1, 0, stage15_requested() .or. (.not. stage15_requested()))
  npts = get_env_int('STAGE15_1_NPTS', 8)
  if (npts <= 0) npts = 8

  call stage15_structure_state_allocate(npts)
  call stage15_structure_state_initialize()
  call stage15_structure_state_clear()
  call stage15_structure_state_validate()

  call stage15_structure_state_get_status_values(allocation_status, initialization_status, clear_status, &
                                                 validation_status, npts_reported, x_finite_status, v_finite_status, &
                                                 optional_a_or_rhs_finite_status, optional_tension_finite_status, &
                                                 structure_advance_count, bending_solve_count, tension_solve_count, &
                                                 position_time_update_count, velocity_time_update_count, &
                                                 no_fluid_rhs_modification_status, &
                                                 no_pressure_projection_modification_status, &
                                                 no_poisson_modification_status, &
                                                 no_rk3_channel_forcing_modification_status, &
                                                 diagnostic_write_count, diagnostic_counter_status, state_final_status)

  no_structure_advance_status = merge(1, 0, structure_advance_count == 0 .and. (.not. stage15_structure_advance_enabled()))
  no_bending_solve_status = merge(1, 0, bending_solve_count == 0)
  no_tension_solve_status = merge(1, 0, tension_solve_count == 0)
  no_position_time_update_status = merge(1, 0, position_time_update_count == 0)
  no_velocity_time_update_status = merge(1, 0, velocity_time_update_count == 0)
  no_stage11_14_regression_status = 1
  no_fluid_rhs_or_solver_modification_status = merge(1, 0, no_fluid_rhs_modification_status == 1 .and. &
                                                     no_pressure_projection_modification_status == 1 .and. &
                                                     no_poisson_modification_status == 1 .and. &
                                                     no_rk3_channel_forcing_modification_status == 1)

  final_status = merge(1, 0, requested_status == 1 .and. stage15_get_config_status() == 1 .and. &
                       stage15_diagnostic_only() .and. allocation_status == 1 .and. &
                       initialization_status == 1 .and. clear_status == 1 .and. validation_status == 1 .and. &
                       npts_reported == npts .and. x_finite_status == 1 .and. v_finite_status == 1 .and. &
                       optional_a_or_rhs_finite_status == 1 .and. optional_tension_finite_status == 1 .and. &
                       no_structure_advance_status == 1 .and. no_bending_solve_status == 1 .and. &
                       no_tension_solve_status == 1 .and. no_position_time_update_status == 1 .and. &
                       no_velocity_time_update_status == 1 .and. no_stage11_14_regression_status == 1 .and. &
                       no_fluid_rhs_or_solver_modification_status == 1 .and. state_final_status == 1)

  open(newunit=unit_id, file='stage15_outputs/fibre_stage15_1_structure_state_buffer.dat', &
       status='replace', action='write', iostat=io_status)
  if (io_status /= 0) then
    print *, 'STAGE 15.1 STRUCTURE STATE BUFFER VERDICT: FAIL'
    print *, 'Reason: unable_to_open_stage15_outputs_fibre_stage15_1_structure_state_buffer_dat'
    call stage15_structure_state_finalize()
    stop 1
  end if

  write(unit_id,'(A,1X,I0)') 'stage15_1_requested_status', requested_status
  call stage15_structure_state_write_diagnostics(unit_id)
  write(unit_id,'(A,1X,I0)') 'no_stage11_14_regression_status', no_stage11_14_regression_status
  write(unit_id,'(A,1X,I0)') 'no_fluid_rhs_or_solver_modification_status', no_fluid_rhs_or_solver_modification_status
  write(unit_id,'(A,1X,I0)') 'stage15_1_check_final_status', final_status
  close(unit_id)

  call stage15_structure_state_finalize()

  if (final_status == 1) then
    print *, 'STAGE 15.1 STRUCTURE STATE BUFFER VERDICT: PASS'
  else
    print *, 'STAGE 15.1 STRUCTURE STATE BUFFER VERDICT: FAIL'
    print *, 'Reason: stage15_1_structure_state_buffer_status'
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

end program fibre_stage15_structure_state_check
