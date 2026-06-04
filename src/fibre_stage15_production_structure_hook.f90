module fibre_stage15_production_structure_hook
  use fibre_stage15_config, only : stage15_config_load, stage15_requested, &
                                   stage15_structure_advance_enabled, stage15_diagnostic_only
  use fibre_stage15_structure_state, only : stage15_structure_state_is_allocated
  implicit none
  private

  integer :: stage15_4_requested_status = 0
  integer :: hook_registered_status = 0
  integer :: hook_apply_count = 0
  integer :: hook_finalize_status = 0
  integer :: diagnostic_only_status = 0
  integer :: noop_status = 0
  integer :: structure_state_available_status = 0
  integer :: production_time_loop_hook_status = 0
  integer :: production_time_loop_connection_count = 0
  integer :: production_structure_advance_count = 0
  integer :: x_position_update_count = 0
  integer :: v_velocity_update_count = 0
  integer :: a_acceleration_update_count = 0
  integer :: bending_solve_count = 0
  integer :: tension_solve_count = 0
  integer :: wall_contact_count = 0
  integer :: multifibre_count = 0
  integer :: rhs_injection_connection_count = 0
  integer :: no_fluid_rhs_modification_status = 1
  integer :: no_pressure_projection_modification_status = 1
  integer :: no_poisson_modification_status = 1
  integer :: no_rk3_channel_forcing_modification_status = 1
  integer :: last_time_step = -1
  integer :: last_rk_substep = -1
  integer :: final_status = 0

  public :: stage15_production_structure_hook_reset
  public :: stage15_production_structure_hook_register
  public :: stage15_production_structure_hook_apply
  public :: stage15_production_structure_hook_finalize
  public :: stage15_production_structure_hook_get_status_values
  public :: stage15_production_structure_hook_write_diagnostics

contains

  subroutine stage15_production_structure_hook_reset()
    call stage15_config_load()
    stage15_4_requested_status = merge(1, 0, stage15_requested())
    hook_registered_status = 0
    hook_apply_count = 0
    hook_finalize_status = 0
    diagnostic_only_status = merge(1, 0, stage15_diagnostic_only())
    noop_status = merge(1, 0, stage15_diagnostic_only() .and. (.not. stage15_structure_advance_enabled()))
    structure_state_available_status = merge(1, 0, stage15_structure_state_is_allocated())
    production_time_loop_hook_status = 0
    production_time_loop_connection_count = 0
    production_structure_advance_count = 0
    x_position_update_count = 0
    v_velocity_update_count = 0
    a_acceleration_update_count = 0
    bending_solve_count = 0
    tension_solve_count = 0
    wall_contact_count = 0
    multifibre_count = 0
    rhs_injection_connection_count = 0
    no_fluid_rhs_modification_status = 1
    no_pressure_projection_modification_status = 1
    no_poisson_modification_status = 1
    no_rk3_channel_forcing_modification_status = 1
    last_time_step = -1
    last_rk_substep = -1
    call update_final_status()
  end subroutine stage15_production_structure_hook_reset

  subroutine stage15_production_structure_hook_register(production_time_loop_context)
    integer, intent(in) :: production_time_loop_context

    if (stage15_4_requested_status /= 1) call stage15_production_structure_hook_reset()
    if (stage15_4_requested_status == 1 .and. noop_status == 1) then
      hook_registered_status = 1
      production_time_loop_hook_status = merge(1, 0, production_time_loop_context == 1)
      production_time_loop_connection_count = merge(1, 0, production_time_loop_context == 1)
    else
      hook_registered_status = 0
    end if
    structure_state_available_status = merge(1, 0, stage15_structure_state_is_allocated())
    call update_final_status()
  end subroutine stage15_production_structure_hook_register

  subroutine stage15_production_structure_hook_apply(time_step, rk_substep)
    integer, intent(in) :: time_step
    integer, intent(in) :: rk_substep

    if (hook_registered_status == 1 .and. noop_status == 1) then
      hook_apply_count = hook_apply_count + 1
      last_time_step = time_step
      last_rk_substep = rk_substep
      structure_state_available_status = merge(1, 0, stage15_structure_state_is_allocated())
    end if
    call update_final_status()
  end subroutine stage15_production_structure_hook_apply

  subroutine stage15_production_structure_hook_finalize()
    hook_finalize_status = 1
    call update_final_status()
    call write_default_diagnostics()
  end subroutine stage15_production_structure_hook_finalize

  subroutine stage15_production_structure_hook_get_status_values(requested_out, registered_out, apply_count_out, &
                                                                 finalize_out, diagnostic_only_out, noop_out, &
                                                                 state_available_out, time_loop_hook_out, &
                                                                 time_loop_count_out, structure_count_out, &
                                                                 x_update_count_out, v_update_count_out, &
                                                                 a_update_count_out, bending_count_out, &
                                                                 tension_count_out, wall_count_out, &
                                                                 multifibre_count_out, rhs_count_out, &
                                                                 no_fluid_rhs_out, no_pressure_projection_out, &
                                                                 no_poisson_out, no_rk3_channel_out, &
                                                                 last_time_step_out, last_rk_substep_out, final_out)
    integer, intent(out) :: requested_out
    integer, intent(out) :: registered_out
    integer, intent(out) :: apply_count_out
    integer, intent(out) :: finalize_out
    integer, intent(out) :: diagnostic_only_out
    integer, intent(out) :: noop_out
    integer, intent(out) :: state_available_out
    integer, intent(out) :: time_loop_hook_out
    integer, intent(out) :: time_loop_count_out
    integer, intent(out) :: structure_count_out
    integer, intent(out) :: x_update_count_out
    integer, intent(out) :: v_update_count_out
    integer, intent(out) :: a_update_count_out
    integer, intent(out) :: bending_count_out
    integer, intent(out) :: tension_count_out
    integer, intent(out) :: wall_count_out
    integer, intent(out) :: multifibre_count_out
    integer, intent(out) :: rhs_count_out
    integer, intent(out) :: no_fluid_rhs_out
    integer, intent(out) :: no_pressure_projection_out
    integer, intent(out) :: no_poisson_out
    integer, intent(out) :: no_rk3_channel_out
    integer, intent(out) :: last_time_step_out
    integer, intent(out) :: last_rk_substep_out
    integer, intent(out) :: final_out

    call update_final_status()
    requested_out = stage15_4_requested_status
    registered_out = hook_registered_status
    apply_count_out = hook_apply_count
    finalize_out = hook_finalize_status
    diagnostic_only_out = diagnostic_only_status
    noop_out = noop_status
    state_available_out = structure_state_available_status
    time_loop_hook_out = production_time_loop_hook_status
    time_loop_count_out = production_time_loop_connection_count
    structure_count_out = production_structure_advance_count
    x_update_count_out = x_position_update_count
    v_update_count_out = v_velocity_update_count
    a_update_count_out = a_acceleration_update_count
    bending_count_out = bending_solve_count
    tension_count_out = tension_solve_count
    wall_count_out = wall_contact_count
    multifibre_count_out = multifibre_count
    rhs_count_out = rhs_injection_connection_count
    no_fluid_rhs_out = no_fluid_rhs_modification_status
    no_pressure_projection_out = no_pressure_projection_modification_status
    no_poisson_out = no_poisson_modification_status
    no_rk3_channel_out = no_rk3_channel_forcing_modification_status
    last_time_step_out = last_time_step
    last_rk_substep_out = last_rk_substep
    final_out = final_status
  end subroutine stage15_production_structure_hook_get_status_values

  subroutine stage15_production_structure_hook_write_diagnostics(unit_id)
    integer, intent(in) :: unit_id

    call update_final_status()
    write(unit_id,'(A,1X,I0)') 'stage15_4_requested_status', stage15_4_requested_status
    write(unit_id,'(A,1X,I0)') 'hook_registered_status', hook_registered_status
    write(unit_id,'(A,1X,I0)') 'hook_apply_count', hook_apply_count
    write(unit_id,'(A,1X,I0)') 'hook_finalize_status', hook_finalize_status
    write(unit_id,'(A,1X,I0)') 'diagnostic_only_status', diagnostic_only_status
    write(unit_id,'(A,1X,I0)') 'noop_status', noop_status
    write(unit_id,'(A,1X,I0)') 'structure_state_available_status', structure_state_available_status
    write(unit_id,'(A,1X,I0)') 'production_time_loop_hook_status', production_time_loop_hook_status
    write(unit_id,'(A,1X,I0)') 'production_time_loop_connection_count', production_time_loop_connection_count
    write(unit_id,'(A,1X,I0)') 'production_structure_advance_count', production_structure_advance_count
    write(unit_id,'(A,1X,I0)') 'x_position_update_count', x_position_update_count
    write(unit_id,'(A,1X,I0)') 'v_velocity_update_count', v_velocity_update_count
    write(unit_id,'(A,1X,I0)') 'a_acceleration_update_count', a_acceleration_update_count
    write(unit_id,'(A,1X,I0)') 'bending_solve_count', bending_solve_count
    write(unit_id,'(A,1X,I0)') 'tension_solve_count', tension_solve_count
    write(unit_id,'(A,1X,I0)') 'wall_contact_count', wall_contact_count
    write(unit_id,'(A,1X,I0)') 'multifibre_count', multifibre_count
    write(unit_id,'(A,1X,I0)') 'rhs_injection_connection_count', rhs_injection_connection_count
    write(unit_id,'(A,1X,I0)') 'no_fluid_rhs_modification_status', no_fluid_rhs_modification_status
    write(unit_id,'(A,1X,I0)') 'no_pressure_projection_modification_status', no_pressure_projection_modification_status
    write(unit_id,'(A,1X,I0)') 'no_poisson_modification_status', no_poisson_modification_status
    write(unit_id,'(A,1X,I0)') 'no_rk3_channel_forcing_modification_status', no_rk3_channel_forcing_modification_status
    write(unit_id,'(A,1X,I0)') 'last_time_step', last_time_step
    write(unit_id,'(A,1X,I0)') 'last_rk_substep', last_rk_substep
    write(unit_id,'(A,1X,I0)') 'final_status', final_status
  end subroutine stage15_production_structure_hook_write_diagnostics

  subroutine update_final_status()
    final_status = merge(1, 0, stage15_4_requested_status == 1 .and. hook_registered_status == 1 .and. &
                         diagnostic_only_status == 1 .and. noop_status == 1 .and. &
                         production_time_loop_hook_status == 1 .and. production_time_loop_connection_count >= 1 .and. &
                         production_structure_advance_count == 0 .and. x_position_update_count == 0 .and. &
                         v_velocity_update_count == 0 .and. a_acceleration_update_count == 0 .and. &
                         bending_solve_count == 0 .and. tension_solve_count == 0 .and. wall_contact_count == 0 .and. &
                         multifibre_count == 0 .and. rhs_injection_connection_count == 0 .and. &
                         no_fluid_rhs_modification_status == 1 .and. &
                         no_pressure_projection_modification_status == 1 .and. &
                         no_poisson_modification_status == 1 .and. &
                         no_rk3_channel_forcing_modification_status == 1)
  end subroutine update_final_status

  subroutine write_default_diagnostics()
    integer :: unit_id
    integer :: io_status

    if (.not. rank0_write_allowed()) return
    call execute_command_line('mkdir -p stage15_outputs')
    open(newunit=unit_id, file='stage15_outputs/fibre_stage15_4_production_structure_hook.dat', &
         status='replace', action='write', iostat=io_status)
    if (io_status /= 0) return
    call stage15_production_structure_hook_write_diagnostics(unit_id)
    close(unit_id)
  end subroutine write_default_diagnostics

  logical function rank0_write_allowed()
    character(len=32) :: value
    integer :: status
    integer :: ios
    integer :: rank_value

    rank0_write_allowed = .true.
    call get_environment_variable('OMPI_COMM_WORLD_RANK', value=value, status=status)
    if (status /= 0) call get_environment_variable('PMI_RANK', value=value, status=status)
    if (status /= 0) call get_environment_variable('MPI_RANK', value=value, status=status)
    if (status == 0) then
      read(value, *, iostat=ios) rank_value
      if (ios == 0) rank0_write_allowed = (rank_value == 0)
    end if
  end function rank0_write_allowed

end module fibre_stage15_production_structure_hook
