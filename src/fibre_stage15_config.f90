module fibre_stage15_config
  implicit none

  private

  public :: stage15_config_load
  public :: stage15_requested
  public :: stage15_structure_advance_enabled
  public :: stage15_diagnostic_only
  public :: stage15_require_stage14_closed
  public :: stage15_get_config_status
  public :: stage15_get_status_values
  public :: stage15_write_config_diagnostics

  logical :: config_loaded = .false.
  logical :: requested_value = .false.
  logical :: structure_advance_requested_value = .false.
  logical :: structure_advance_enabled_value = .false.
  logical :: diagnostic_only_value = .true.
  logical :: require_stage14_closed_value = .true.

  integer :: requested_flag = 0
  integer :: disabled_by_default_status = 1
  integer :: request_parse_status = 1
  integer :: structure_advance_requested_flag = 0
  integer :: structure_advance_enabled_flag = 0
  integer :: structure_advance_parse_status = 1
  integer :: structure_advance_disabled_status = 1
  integer :: structure_advance_blocked_status = 1
  integer :: diagnostic_only_flag = 1
  integer :: diagnostic_only_default_status = 1
  integer :: diagnostic_only_parse_status = 1
  integer :: diagnostic_only_enforced_status = 1
  integer :: require_stage14_closed_flag = 1
  integer :: require_stage14_closed_default_status = 1
  integer :: require_stage14_closed_parse_status = 1
  integer :: no_structure_state_allocation_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_bending_solve_status = 1
  integer :: no_tension_solve_status = 1
  integer :: no_fibre_position_update_status = 1
  integer :: no_fibre_velocity_update_status = 1
  integer :: no_wall_contact_status = 1
  integer :: no_multifibre_logic_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: no_rhs_modification_status = 1
  integer :: no_pressure_modification_status = 1
  integer :: no_projection_modification_status = 1
  integer :: no_poisson_modification_status = 1
  integer :: no_rk3_modification_status = 1
  integer :: no_channel_forcing_modification_status = 1
  integer :: no_production_ibm_forcing_status = 1
  integer :: no_stage11_14_reimplementation_status = 1
  integer :: no_production_behavior_change_status = 1
  integer :: config_status = 1

contains

  subroutine stage15_config_load()
    character(len=256) :: env_value
    integer :: env_status
    logical :: parsed_bool
    integer :: parsed_status

    config_loaded = .true.

    requested_value = .false.
    structure_advance_requested_value = .false.
    structure_advance_enabled_value = .false.
    diagnostic_only_value = .true.
    require_stage14_closed_value = .true.

    requested_flag = 0
    disabled_by_default_status = 1
    request_parse_status = 1
    structure_advance_requested_flag = 0
    structure_advance_enabled_flag = 0
    structure_advance_parse_status = 1
    structure_advance_disabled_status = 1
    structure_advance_blocked_status = 1
    diagnostic_only_flag = 1
    diagnostic_only_default_status = 1
    diagnostic_only_parse_status = 1
    diagnostic_only_enforced_status = 1
    require_stage14_closed_flag = 1
    require_stage14_closed_default_status = 1
    require_stage14_closed_parse_status = 1
    no_structure_state_allocation_status = 1
    no_structure_advance_status = 1
    no_bending_solve_status = 1
    no_tension_solve_status = 1
    no_fibre_position_update_status = 1
    no_fibre_velocity_update_status = 1
    no_wall_contact_status = 1
    no_multifibre_logic_status = 1
    no_fluid_field_access_status = 1
    no_fluid_field_modification_status = 1
    no_rhs_modification_status = 1
    no_pressure_modification_status = 1
    no_projection_modification_status = 1
    no_poisson_modification_status = 1
    no_rk3_modification_status = 1
    no_channel_forcing_modification_status = 1
    no_production_ibm_forcing_status = 1
    no_stage11_14_reimplementation_status = 1
    no_production_behavior_change_status = 1
    config_status = 1

    call get_environment_variable('X3D_STAGE15_ENABLE', value=env_value, status=env_status)
    if (env_status == 0) then
      call parse_bool(env_value, parsed_bool, parsed_status)
      request_parse_status = parsed_status
      if (parsed_status == 1) then
        requested_value = parsed_bool
      else
        requested_value = .false.
      end if
    end if
    requested_flag = logical_to_int(requested_value)
    disabled_by_default_status = logical_to_int(.not. requested_value)

    call get_environment_variable('X3D_STAGE15_STRUCTURE_ADVANCE_ENABLE', value=env_value, status=env_status)
    if (env_status == 0) then
      call parse_bool(env_value, parsed_bool, parsed_status)
      structure_advance_parse_status = parsed_status
      if (parsed_status == 1) then
        structure_advance_requested_value = parsed_bool
      else
        structure_advance_requested_value = .false.
      end if
    end if
    structure_advance_requested_flag = logical_to_int(structure_advance_requested_value)
    structure_advance_enabled_value = .false.
    structure_advance_enabled_flag = 0
    structure_advance_disabled_status = 1
    structure_advance_blocked_status = 1

    call get_environment_variable('X3D_STAGE15_DIAGNOSTIC_ONLY', value=env_value, status=env_status)
    if (env_status == 0) then
      call parse_bool(env_value, parsed_bool, parsed_status)
      diagnostic_only_parse_status = parsed_status
    end if
    diagnostic_only_value = .true.
    diagnostic_only_flag = 1
    diagnostic_only_default_status = 1
    diagnostic_only_enforced_status = 1

    call get_environment_variable('X3D_STAGE15_REQUIRE_STAGE14_CLOSED', value=env_value, status=env_status)
    if (env_status == 0) then
      call parse_bool(env_value, parsed_bool, parsed_status)
      require_stage14_closed_parse_status = parsed_status
      if (parsed_status == 1) then
        require_stage14_closed_value = parsed_bool
      else
        require_stage14_closed_value = .true.
      end if
    end if
    require_stage14_closed_flag = logical_to_int(require_stage14_closed_value)
    require_stage14_closed_default_status = logical_to_int(require_stage14_closed_value)

    config_status = logical_to_int(request_parse_status == 1 .and. &
                                   structure_advance_parse_status == 1 .and. &
                                   structure_advance_enabled_flag == 0 .and. &
                                   structure_advance_disabled_status == 1 .and. &
                                   diagnostic_only_parse_status == 1 .and. &
                                   diagnostic_only_flag == 1 .and. &
                                   diagnostic_only_enforced_status == 1 .and. &
                                   require_stage14_closed_parse_status == 1 .and. &
                                   no_structure_state_allocation_status == 1 .and. &
                                   no_structure_advance_status == 1 .and. &
                                   no_bending_solve_status == 1 .and. &
                                   no_tension_solve_status == 1 .and. &
                                   no_fibre_position_update_status == 1 .and. &
                                   no_fibre_velocity_update_status == 1 .and. &
                                   no_wall_contact_status == 1 .and. &
                                   no_multifibre_logic_status == 1 .and. &
                                   no_fluid_field_access_status == 1 .and. &
                                   no_fluid_field_modification_status == 1 .and. &
                                   no_rhs_modification_status == 1 .and. &
                                   no_pressure_modification_status == 1 .and. &
                                   no_projection_modification_status == 1 .and. &
                                   no_poisson_modification_status == 1 .and. &
                                   no_rk3_modification_status == 1 .and. &
                                   no_channel_forcing_modification_status == 1 .and. &
                                   no_production_ibm_forcing_status == 1 .and. &
                                   no_stage11_14_reimplementation_status == 1 .and. &
                                   no_production_behavior_change_status == 1)
  end subroutine stage15_config_load

  logical function stage15_requested()
    if (.not. config_loaded) call stage15_config_load()
    stage15_requested = requested_value
  end function stage15_requested

  logical function stage15_structure_advance_enabled()
    if (.not. config_loaded) call stage15_config_load()
    stage15_structure_advance_enabled = structure_advance_enabled_value
  end function stage15_structure_advance_enabled

  logical function stage15_diagnostic_only()
    if (.not. config_loaded) call stage15_config_load()
    stage15_diagnostic_only = diagnostic_only_value
  end function stage15_diagnostic_only

  logical function stage15_require_stage14_closed()
    if (.not. config_loaded) call stage15_config_load()
    stage15_require_stage14_closed = require_stage14_closed_value
  end function stage15_require_stage14_closed

  integer function stage15_get_config_status()
    if (.not. config_loaded) call stage15_config_load()
    stage15_get_config_status = config_status
  end function stage15_get_config_status

  subroutine stage15_get_status_values(out_requested_flag, &
                                       out_disabled_by_default_status, &
                                       out_request_parse_status, &
                                       out_structure_advance_requested_flag, &
                                       out_structure_advance_enabled_flag, &
                                       out_structure_advance_parse_status, &
                                       out_structure_advance_disabled_status, &
                                       out_structure_advance_blocked_status, &
                                       out_diagnostic_only_flag, &
                                       out_diagnostic_only_default_status, &
                                       out_diagnostic_only_parse_status, &
                                       out_diagnostic_only_enforced_status, &
                                       out_require_stage14_closed_flag, &
                                       out_require_stage14_closed_default_status, &
                                       out_require_stage14_closed_parse_status, &
                                       out_no_structure_state_allocation_status, &
                                       out_no_structure_advance_status, &
                                       out_no_bending_solve_status, &
                                       out_no_tension_solve_status, &
                                       out_no_fibre_position_update_status, &
                                       out_no_fibre_velocity_update_status, &
                                       out_no_wall_contact_status, &
                                       out_no_multifibre_logic_status, &
                                       out_no_fluid_field_access_status, &
                                       out_no_fluid_field_modification_status, &
                                       out_no_rhs_modification_status, &
                                       out_no_pressure_modification_status, &
                                       out_no_projection_modification_status, &
                                       out_no_poisson_modification_status, &
                                       out_no_rk3_modification_status, &
                                       out_no_channel_forcing_modification_status, &
                                       out_no_production_ibm_forcing_status, &
                                       out_no_stage11_14_reimplementation_status, &
                                       out_no_production_behavior_change_status, &
                                       out_config_status)
    integer, intent(out) :: out_requested_flag
    integer, intent(out) :: out_disabled_by_default_status
    integer, intent(out) :: out_request_parse_status
    integer, intent(out) :: out_structure_advance_requested_flag
    integer, intent(out) :: out_structure_advance_enabled_flag
    integer, intent(out) :: out_structure_advance_parse_status
    integer, intent(out) :: out_structure_advance_disabled_status
    integer, intent(out) :: out_structure_advance_blocked_status
    integer, intent(out) :: out_diagnostic_only_flag
    integer, intent(out) :: out_diagnostic_only_default_status
    integer, intent(out) :: out_diagnostic_only_parse_status
    integer, intent(out) :: out_diagnostic_only_enforced_status
    integer, intent(out) :: out_require_stage14_closed_flag
    integer, intent(out) :: out_require_stage14_closed_default_status
    integer, intent(out) :: out_require_stage14_closed_parse_status
    integer, intent(out) :: out_no_structure_state_allocation_status
    integer, intent(out) :: out_no_structure_advance_status
    integer, intent(out) :: out_no_bending_solve_status
    integer, intent(out) :: out_no_tension_solve_status
    integer, intent(out) :: out_no_fibre_position_update_status
    integer, intent(out) :: out_no_fibre_velocity_update_status
    integer, intent(out) :: out_no_wall_contact_status
    integer, intent(out) :: out_no_multifibre_logic_status
    integer, intent(out) :: out_no_fluid_field_access_status
    integer, intent(out) :: out_no_fluid_field_modification_status
    integer, intent(out) :: out_no_rhs_modification_status
    integer, intent(out) :: out_no_pressure_modification_status
    integer, intent(out) :: out_no_projection_modification_status
    integer, intent(out) :: out_no_poisson_modification_status
    integer, intent(out) :: out_no_rk3_modification_status
    integer, intent(out) :: out_no_channel_forcing_modification_status
    integer, intent(out) :: out_no_production_ibm_forcing_status
    integer, intent(out) :: out_no_stage11_14_reimplementation_status
    integer, intent(out) :: out_no_production_behavior_change_status
    integer, intent(out) :: out_config_status

    if (.not. config_loaded) call stage15_config_load()

    out_requested_flag = requested_flag
    out_disabled_by_default_status = disabled_by_default_status
    out_request_parse_status = request_parse_status
    out_structure_advance_requested_flag = structure_advance_requested_flag
    out_structure_advance_enabled_flag = structure_advance_enabled_flag
    out_structure_advance_parse_status = structure_advance_parse_status
    out_structure_advance_disabled_status = structure_advance_disabled_status
    out_structure_advance_blocked_status = structure_advance_blocked_status
    out_diagnostic_only_flag = diagnostic_only_flag
    out_diagnostic_only_default_status = diagnostic_only_default_status
    out_diagnostic_only_parse_status = diagnostic_only_parse_status
    out_diagnostic_only_enforced_status = diagnostic_only_enforced_status
    out_require_stage14_closed_flag = require_stage14_closed_flag
    out_require_stage14_closed_default_status = require_stage14_closed_default_status
    out_require_stage14_closed_parse_status = require_stage14_closed_parse_status
    out_no_structure_state_allocation_status = no_structure_state_allocation_status
    out_no_structure_advance_status = no_structure_advance_status
    out_no_bending_solve_status = no_bending_solve_status
    out_no_tension_solve_status = no_tension_solve_status
    out_no_fibre_position_update_status = no_fibre_position_update_status
    out_no_fibre_velocity_update_status = no_fibre_velocity_update_status
    out_no_wall_contact_status = no_wall_contact_status
    out_no_multifibre_logic_status = no_multifibre_logic_status
    out_no_fluid_field_access_status = no_fluid_field_access_status
    out_no_fluid_field_modification_status = no_fluid_field_modification_status
    out_no_rhs_modification_status = no_rhs_modification_status
    out_no_pressure_modification_status = no_pressure_modification_status
    out_no_projection_modification_status = no_projection_modification_status
    out_no_poisson_modification_status = no_poisson_modification_status
    out_no_rk3_modification_status = no_rk3_modification_status
    out_no_channel_forcing_modification_status = no_channel_forcing_modification_status
    out_no_production_ibm_forcing_status = no_production_ibm_forcing_status
    out_no_stage11_14_reimplementation_status = no_stage11_14_reimplementation_status
    out_no_production_behavior_change_status = no_production_behavior_change_status
    out_config_status = config_status
  end subroutine stage15_get_status_values

  subroutine stage15_write_config_diagnostics(unit_id)
    integer, intent(in) :: unit_id

    if (.not. config_loaded) call stage15_config_load()

    write(unit_id, '(A,1X,I0)') 'stage15_0_requested_flag', requested_flag
    write(unit_id, '(A,1X,I0)') 'stage15_0_disabled_by_default_status', disabled_by_default_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_request_parse_status', request_parse_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_structure_advance_requested_flag', structure_advance_requested_flag
    write(unit_id, '(A,1X,I0)') 'stage15_0_structure_advance_enabled_flag', structure_advance_enabled_flag
    write(unit_id, '(A,1X,I0)') 'stage15_0_structure_advance_parse_status', structure_advance_parse_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_structure_advance_disabled_status', structure_advance_disabled_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_structure_advance_blocked_status', structure_advance_blocked_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_diagnostic_only_flag', diagnostic_only_flag
    write(unit_id, '(A,1X,I0)') 'stage15_0_diagnostic_only_default_status', diagnostic_only_default_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_diagnostic_only_parse_status', diagnostic_only_parse_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_diagnostic_only_enforced_status', diagnostic_only_enforced_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_require_stage14_closed_flag', require_stage14_closed_flag
    write(unit_id, '(A,1X,I0)') 'stage15_0_require_stage14_closed_default_status', require_stage14_closed_default_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_require_stage14_closed_parse_status', require_stage14_closed_parse_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_structure_state_allocation_status', no_structure_state_allocation_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_structure_advance_status', no_structure_advance_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_bending_solve_status', no_bending_solve_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_tension_solve_status', no_tension_solve_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_fibre_position_update_status', no_fibre_position_update_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_fibre_velocity_update_status', no_fibre_velocity_update_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_wall_contact_status', no_wall_contact_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_multifibre_logic_status', no_multifibre_logic_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_fluid_field_access_status', no_fluid_field_access_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_fluid_field_modification_status', no_fluid_field_modification_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_rhs_modification_status', no_rhs_modification_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_pressure_modification_status', no_pressure_modification_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_projection_modification_status', no_projection_modification_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_poisson_modification_status', no_poisson_modification_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_rk3_modification_status', no_rk3_modification_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_channel_forcing_modification_status', no_channel_forcing_modification_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_production_ibm_forcing_status', no_production_ibm_forcing_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_stage11_14_reimplementation_status', no_stage11_14_reimplementation_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_no_production_behavior_change_status', no_production_behavior_change_status
    write(unit_id, '(A,1X,I0)') 'stage15_0_config_status', config_status
  end subroutine stage15_write_config_diagnostics

  integer function logical_to_int(value)
    logical, intent(in) :: value
    if (value) then
      logical_to_int = 1
    else
      logical_to_int = 0
    end if
  end function logical_to_int

  subroutine parse_bool(text, value, parse_status)
    character(len=*), intent(in) :: text
    logical, intent(out) :: value
    integer, intent(out) :: parse_status
    character(len=256) :: lowered

    lowered = to_lower(trim(adjustl(text)))
    select case (trim(lowered))
    case ('1', 'true', 'yes', 'on')
      value = .true.
      parse_status = 1
    case ('0', 'false', 'no', 'off')
      value = .false.
      parse_status = 1
    case default
      value = .false.
      parse_status = 0
    end select
  end subroutine parse_bool

  function to_lower(text) result(lowered)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lowered
    integer :: i
    integer :: code

    lowered = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lowered(i:i) = achar(code + iachar('a') - iachar('A'))
      end if
    end do
  end function to_lower

end module fibre_stage15_config
