program fibre_stage16_one_fibre_case_check
  use decomp_2d_constants, only : mytype
  use fibre_stage16_one_fibre_case, only : stage16_one_fibre_case_load_from_environment, &
       stage16_one_fibre_case_write_diagnostics, stage16_one_fibre_case_set_for_test, &
       stage16_one_fibre_case_get_status_values
  implicit none

  integer :: io_unit
  integer :: final_status
  integer :: case_status
  integer :: one_fibre_count_status
  integer :: npts_valid_status
  integer :: point_spacing_status
  integer :: geometry_status
  integer :: velocity_status
  integer :: acceleration_status
  integer :: contact_status
  integer :: multifibre_status
  integer :: invalid_npts_rejection_status
  integer :: invalid_geometry_rejection_status
  integer :: invalid_velocity_rejection_status
  integer :: invalid_acceleration_rejection_status
  integer :: wall_contact_rejection_status
  integer :: multifibre_rejection_status
  integer :: no_production_hook_status
  integer :: no_structure_advance_status
  integer :: no_rhs_modification_status
  integer :: no_bending_solve_status
  integer :: no_tension_solve_status
  integer :: no_wall_contact_status
  integer :: no_multifibre_status
  integer :: no_pressure_projection_modification_status
  integer :: no_poisson_modification_status
  integer :: no_rk3_channel_forcing_modification_status
  integer :: no_legacy_ibm_forcing_status

  call execute_command_line('mkdir -p stage16_outputs')

  call stage16_one_fibre_case_load_from_environment()
  call stage16_one_fibre_case_get_status_values(case_status, one_fibre_count_status, npts_valid_status, &
       point_spacing_status, geometry_status, velocity_status, acceleration_status, contact_status, &
       multifibre_status, final_status)

  call stage16_one_fibre_case_set_for_test(1, 1, .false., .false., .false., .false.)
  call stage16_one_fibre_case_get_status_values(case_out=case_status, count_out=one_fibre_count_status, &
       npts_out=npts_valid_status, spacing_out=point_spacing_status, geometry_out=geometry_status, &
       velocity_out=velocity_status, acceleration_out=acceleration_status, contact_out=contact_status, &
       multifibre_out=multifibre_status, final_out=final_status)
  invalid_npts_rejection_status = logical_to_int(final_status == 0 .and. npts_valid_status == 0)

  call stage16_one_fibre_case_set_for_test(8, 1, .true., .false., .false., .false.)
  call stage16_one_fibre_case_get_status_values(case_out=case_status, count_out=one_fibre_count_status, &
       npts_out=npts_valid_status, spacing_out=point_spacing_status, geometry_out=geometry_status, &
       velocity_out=velocity_status, acceleration_out=acceleration_status, contact_out=contact_status, &
       multifibre_out=multifibre_status, final_out=final_status)
  invalid_geometry_rejection_status = logical_to_int(final_status == 0 .and. geometry_status == 0)

  call stage16_one_fibre_case_set_for_test(8, 1, .false., .true., .false., .false.)
  call stage16_one_fibre_case_get_status_values(case_out=case_status, count_out=one_fibre_count_status, &
       npts_out=npts_valid_status, spacing_out=point_spacing_status, geometry_out=geometry_status, &
       velocity_out=velocity_status, acceleration_out=acceleration_status, contact_out=contact_status, &
       multifibre_out=multifibre_status, final_out=final_status)
  invalid_velocity_rejection_status = logical_to_int(final_status == 0 .and. velocity_status == 0)

  call stage16_one_fibre_case_set_for_test(8, 1, .false., .false., .true., .false.)
  call stage16_one_fibre_case_get_status_values(case_out=case_status, count_out=one_fibre_count_status, &
       npts_out=npts_valid_status, spacing_out=point_spacing_status, geometry_out=geometry_status, &
       velocity_out=velocity_status, acceleration_out=acceleration_status, contact_out=contact_status, &
       multifibre_out=multifibre_status, final_out=final_status)
  invalid_acceleration_rejection_status = logical_to_int(final_status == 0 .and. acceleration_status == 0)

  call stage16_one_fibre_case_set_for_test(8, 1, .false., .false., .false., .true.)
  call stage16_one_fibre_case_get_status_values(case_out=case_status, count_out=one_fibre_count_status, &
       npts_out=npts_valid_status, spacing_out=point_spacing_status, geometry_out=geometry_status, &
       velocity_out=velocity_status, acceleration_out=acceleration_status, contact_out=contact_status, &
       multifibre_out=multifibre_status, final_out=final_status)
  wall_contact_rejection_status = logical_to_int(final_status == 0 .and. contact_status == 0)

  call stage16_one_fibre_case_set_for_test(8, 2, .false., .false., .false., .false.)
  call stage16_one_fibre_case_get_status_values(case_out=case_status, count_out=one_fibre_count_status, &
       npts_out=npts_valid_status, spacing_out=point_spacing_status, geometry_out=geometry_status, &
       velocity_out=velocity_status, acceleration_out=acceleration_status, contact_out=contact_status, &
       multifibre_out=multifibre_status, final_out=final_status)
  multifibre_rejection_status = logical_to_int(final_status == 0 .and. multifibre_status == 0)

  call stage16_one_fibre_case_load_from_environment()
  call stage16_one_fibre_case_get_status_values(case_status, one_fibre_count_status, npts_valid_status, &
       point_spacing_status, geometry_status, velocity_status, acceleration_status, contact_status, &
       multifibre_status, final_status)

  no_production_hook_status = 1
  no_structure_advance_status = 1
  no_rhs_modification_status = 1
  no_bending_solve_status = 1
  no_tension_solve_status = 1
  no_wall_contact_status = 1
  no_multifibre_status = 1
  no_pressure_projection_modification_status = 1
  no_poisson_modification_status = 1
  no_rk3_channel_forcing_modification_status = 1
  no_legacy_ibm_forcing_status = 1

  final_status = logical_to_int(final_status == 1 .and. invalid_npts_rejection_status == 1 .and. &
       invalid_geometry_rejection_status == 1 .and. invalid_velocity_rejection_status == 1 .and. &
       invalid_acceleration_rejection_status == 1 .and. wall_contact_rejection_status == 1 .and. &
       multifibre_rejection_status == 1)

  if (rank0_write_allowed()) then
    open(newunit=io_unit, file='stage16_outputs/fibre_stage16_2_one_fibre_case_definition.dat', &
         status='replace', action='write')
    call stage16_one_fibre_case_write_diagnostics(io_unit)
    write(io_unit,'(A,1X,I0)') 'invalid_npts_rejection_status', invalid_npts_rejection_status
    write(io_unit,'(A,1X,I0)') 'invalid_geometry_rejection_status', invalid_geometry_rejection_status
    write(io_unit,'(A,1X,I0)') 'invalid_velocity_rejection_status', invalid_velocity_rejection_status
    write(io_unit,'(A,1X,I0)') 'invalid_acceleration_rejection_status', invalid_acceleration_rejection_status
    write(io_unit,'(A,1X,I0)') 'wall_contact_rejection_status', wall_contact_rejection_status
    write(io_unit,'(A,1X,I0)') 'multifibre_rejection_status', multifibre_rejection_status
    write(io_unit,'(A,1X,I0)') 'no_production_hook_status', no_production_hook_status
    write(io_unit,'(A,1X,I0)') 'no_structure_advance_status', no_structure_advance_status
    write(io_unit,'(A,1X,I0)') 'no_rhs_modification_status', no_rhs_modification_status
    write(io_unit,'(A,1X,I0)') 'no_bending_solve_status', no_bending_solve_status
    write(io_unit,'(A,1X,I0)') 'no_tension_solve_status', no_tension_solve_status
    write(io_unit,'(A,1X,I0)') 'no_wall_contact_status', no_wall_contact_status
    write(io_unit,'(A,1X,I0)') 'no_multifibre_status', no_multifibre_status
    write(io_unit,'(A,1X,I0)') 'no_pressure_projection_modification_status', no_pressure_projection_modification_status
    write(io_unit,'(A,1X,I0)') 'no_poisson_modification_status', no_poisson_modification_status
    write(io_unit,'(A,1X,I0)') 'no_rk3_channel_forcing_modification_status', no_rk3_channel_forcing_modification_status
    write(io_unit,'(A,1X,I0)') 'no_legacy_ibm_forcing_status', no_legacy_ibm_forcing_status
    write(io_unit,'(A,1X,I0)') 'final_status', final_status
    close(io_unit)
  end if

  if (final_status /= 1) stop 1

contains

  integer function logical_to_int(value)
    logical, intent(in) :: value
    if (value) then
      logical_to_int = 1
    else
      logical_to_int = 0
    end if
  end function logical_to_int

  logical function rank0_write_allowed()
    character(len=32) :: env_value
    integer :: env_status
    integer :: rank_value
    integer :: read_status

    rank0_write_allowed = .true.
    call get_environment_variable('OMPI_COMM_WORLD_RANK', value=env_value, status=env_status)
    if (env_status /= 0) call get_environment_variable('PMI_RANK', value=env_value, status=env_status)
    if (env_status /= 0) call get_environment_variable('MPI_LOCALRANKID', value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=read_status) rank_value
      if (read_status == 0) rank0_write_allowed = (rank_value == 0)
    end if
  end function rank0_write_allowed

end program fibre_stage16_one_fibre_case_check
