module fibre_stage14_config
  use decomp_2d_constants, only : mytype
  implicit none

  private

  public :: stage14_config_load
  public :: stage14_requested
  public :: stage14_rhs_injection_enabled
  public :: stage14_get_injection_gain
  public :: stage14_get_max_steps
  public :: stage14_require_stage13
  public :: stage14_get_status_values

  integer, parameter :: default_max_steps = 3

  logical :: config_loaded = .false.
  logical :: requested_value = .false.
  logical :: rhs_injection_enabled_value = .false.
  real(mytype) :: injection_gain_value = 0.0_mytype
  integer :: max_steps_value = default_max_steps
  logical :: require_stage13_value = .true.
  logical :: diagnostic_only_value = .true.

  integer :: requested_flag = 0
  integer :: rhs_injection_enabled_flag = 0
  integer :: disabled_by_default_status = 1
  integer :: gain_parse_status = 1
  integer :: gain_default_zero_status = 1
  integer :: gain_finite_status = 1
  integer :: safe_fallback_status = 1
  integer :: max_steps_parse_status = 1
  integer :: require_stage13_status = 1
  integer :: require_stage13_parse_status = 1
  integer :: diagnostic_only_status = 1
  integer :: diagnostic_only_forced_status = 1
  integer :: no_rhs_buffer_allocation_status = 1
  integer :: no_rhs_injection_status = 1
  integer :: no_production_hook_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: no_pressure_modification_status = 1
  integer :: no_projection_modification_status = 1
  integer :: no_rk3_modification_status = 1
  integer :: no_channel_forcing_modification_status = 1
  integer :: no_production_ibm_forcing_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: config_status = 1

contains

  subroutine stage14_config_load()
    character(len=256) :: env_value
    integer :: env_status
    integer :: parsed_int
    integer :: parse_iostat
    logical :: parsed_bool
    integer :: bool_parse_status
    real(mytype) :: parsed_gain

    config_loaded = .true.

    requested_value = .false.
    rhs_injection_enabled_value = .false.
    injection_gain_value = 0.0_mytype
    max_steps_value = default_max_steps
    require_stage13_value = .true.
    diagnostic_only_value = .true.

    requested_flag = 0
    rhs_injection_enabled_flag = 0
    disabled_by_default_status = 1
    gain_parse_status = 1
    gain_default_zero_status = 1
    gain_finite_status = 1
    safe_fallback_status = 1
    max_steps_parse_status = 1
    require_stage13_status = 1
    require_stage13_parse_status = 1
    diagnostic_only_status = 1
    diagnostic_only_forced_status = 1
    no_rhs_buffer_allocation_status = 1
    no_rhs_injection_status = 1
    no_production_hook_status = 1
    no_fluid_field_access_status = 1
    no_fluid_field_modification_status = 1
    no_pressure_modification_status = 1
    no_projection_modification_status = 1
    no_rk3_modification_status = 1
    no_channel_forcing_modification_status = 1
    no_production_ibm_forcing_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    config_status = 1

    call get_environment_variable('X3D_STAGE14_RHS_INJECTION', value=env_value, status=env_status)
    if (env_status == 0) then
      if (trim(adjustl(env_value)) == '1') then
        requested_value = .true.
        rhs_injection_enabled_value = .true.
      end if
    end if
    requested_flag = logical_to_int(requested_value)
    rhs_injection_enabled_flag = logical_to_int(rhs_injection_enabled_value)

    call get_environment_variable('X3D_STAGE14_INJECTION_GAIN', value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=parse_iostat) parsed_gain
      if (parse_iostat == 0 .and. finite_real(parsed_gain)) then
        injection_gain_value = parsed_gain
        gain_parse_status = 1
      else
        injection_gain_value = 0.0_mytype
        gain_parse_status = 0
        safe_fallback_status = 1
      end if
    else
      injection_gain_value = 0.0_mytype
      gain_parse_status = 1
    end if
    gain_finite_status = logical_to_int(finite_real(injection_gain_value))
    gain_default_zero_status = logical_to_int(injection_gain_value == 0.0_mytype)

    call get_environment_variable('X3D_STAGE14_MAX_STEPS', value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=parse_iostat) parsed_int
      if (parse_iostat == 0 .and. parsed_int > 0) then
        max_steps_value = parsed_int
        max_steps_parse_status = 1
      else
        max_steps_value = default_max_steps
        max_steps_parse_status = 0
        safe_fallback_status = 1
      end if
    else
      max_steps_value = default_max_steps
      max_steps_parse_status = 1
    end if

    call get_environment_variable('X3D_STAGE14_REQUIRE_STAGE13', value=env_value, status=env_status)
    if (env_status == 0) then
      call parse_bool(env_value, parsed_bool, bool_parse_status)
      if (bool_parse_status == 1) then
        require_stage13_value = parsed_bool
        require_stage13_parse_status = 1
      else
        require_stage13_value = .true.
        require_stage13_parse_status = 0
        safe_fallback_status = 1
      end if
    else
      require_stage13_value = .true.
      require_stage13_parse_status = 1
    end if
    require_stage13_status = logical_to_int(require_stage13_value)

    call get_environment_variable('X3D_STAGE14_DIAGNOSTIC_ONLY', value=env_value, status=env_status)
    if (env_status == 0) then
      call parse_bool(env_value, parsed_bool, bool_parse_status)
      if (bool_parse_status == 0) then
        safe_fallback_status = 1
      end if
    end if
    diagnostic_only_value = .true.
    diagnostic_only_status = 1
    diagnostic_only_forced_status = 1

    config_status = logical_to_int(gain_finite_status == 1 .and. &
                                   max_steps_value > 0 .and. &
                                   require_stage13_status == 1 .and. &
                                   diagnostic_only_status == 1 .and. &
                                   no_rhs_buffer_allocation_status == 1 .and. &
                                   no_rhs_injection_status == 1 .and. &
                                   no_production_hook_status == 1 .and. &
                                   no_fluid_field_access_status == 1 .and. &
                                   no_fluid_field_modification_status == 1 .and. &
                                   no_pressure_modification_status == 1 .and. &
                                   no_projection_modification_status == 1 .and. &
                                   no_rk3_modification_status == 1 .and. &
                                   no_channel_forcing_modification_status == 1 .and. &
                                   no_production_ibm_forcing_status == 1 .and. &
                                   no_feedback_application_status == 1 .and. &
                                   no_twoway_force_status == 1 .and. &
                                   no_structure_advance_status == 1)
  end subroutine stage14_config_load

  logical function stage14_requested()
    if (.not. config_loaded) call stage14_config_load()
    stage14_requested = requested_value
  end function stage14_requested

  logical function stage14_rhs_injection_enabled()
    if (.not. config_loaded) call stage14_config_load()
    stage14_rhs_injection_enabled = rhs_injection_enabled_value
  end function stage14_rhs_injection_enabled

  real(mytype) function stage14_get_injection_gain()
    if (.not. config_loaded) call stage14_config_load()
    stage14_get_injection_gain = injection_gain_value
  end function stage14_get_injection_gain

  integer function stage14_get_max_steps()
    if (.not. config_loaded) call stage14_config_load()
    stage14_get_max_steps = max_steps_value
  end function stage14_get_max_steps

  logical function stage14_require_stage13()
    if (.not. config_loaded) call stage14_config_load()
    stage14_require_stage13 = require_stage13_value
  end function stage14_require_stage13

  subroutine stage14_get_status_values(requested_flag_out, &
                                       rhs_injection_enabled_flag_out, &
                                       disabled_by_default_status_out, &
                                       gain_parse_status_out, &
                                       gain_default_zero_status_out, &
                                       gain_finite_status_out, &
                                       safe_fallback_status_out, &
                                       max_steps_parse_status_out, &
                                       require_stage13_status_out, &
                                       require_stage13_parse_status_out, &
                                       diagnostic_only_status_out, &
                                       diagnostic_only_forced_status_out, &
                                       no_rhs_buffer_allocation_status_out, &
                                       no_rhs_injection_status_out, &
                                       no_production_hook_status_out, &
                                       no_fluid_field_access_status_out, &
                                       no_fluid_field_modification_status_out, &
                                       no_pressure_modification_status_out, &
                                       no_projection_modification_status_out, &
                                       no_rk3_modification_status_out, &
                                       no_channel_forcing_modification_status_out, &
                                       no_production_ibm_forcing_status_out, &
                                       no_feedback_application_status_out, &
                                       no_twoway_force_status_out, &
                                       no_structure_advance_status_out, &
                                       config_status_out)
    integer, intent(out) :: requested_flag_out
    integer, intent(out) :: rhs_injection_enabled_flag_out
    integer, intent(out) :: disabled_by_default_status_out
    integer, intent(out) :: gain_parse_status_out
    integer, intent(out) :: gain_default_zero_status_out
    integer, intent(out) :: gain_finite_status_out
    integer, intent(out) :: safe_fallback_status_out
    integer, intent(out) :: max_steps_parse_status_out
    integer, intent(out) :: require_stage13_status_out
    integer, intent(out) :: require_stage13_parse_status_out
    integer, intent(out) :: diagnostic_only_status_out
    integer, intent(out) :: diagnostic_only_forced_status_out
    integer, intent(out) :: no_rhs_buffer_allocation_status_out
    integer, intent(out) :: no_rhs_injection_status_out
    integer, intent(out) :: no_production_hook_status_out
    integer, intent(out) :: no_fluid_field_access_status_out
    integer, intent(out) :: no_fluid_field_modification_status_out
    integer, intent(out) :: no_pressure_modification_status_out
    integer, intent(out) :: no_projection_modification_status_out
    integer, intent(out) :: no_rk3_modification_status_out
    integer, intent(out) :: no_channel_forcing_modification_status_out
    integer, intent(out) :: no_production_ibm_forcing_status_out
    integer, intent(out) :: no_feedback_application_status_out
    integer, intent(out) :: no_twoway_force_status_out
    integer, intent(out) :: no_structure_advance_status_out
    integer, intent(out) :: config_status_out

    if (.not. config_loaded) call stage14_config_load()

    requested_flag_out = requested_flag
    rhs_injection_enabled_flag_out = rhs_injection_enabled_flag
    disabled_by_default_status_out = disabled_by_default_status
    gain_parse_status_out = gain_parse_status
    gain_default_zero_status_out = gain_default_zero_status
    gain_finite_status_out = gain_finite_status
    safe_fallback_status_out = safe_fallback_status
    max_steps_parse_status_out = max_steps_parse_status
    require_stage13_status_out = require_stage13_status
    require_stage13_parse_status_out = require_stage13_parse_status
    diagnostic_only_status_out = diagnostic_only_status
    diagnostic_only_forced_status_out = diagnostic_only_forced_status
    no_rhs_buffer_allocation_status_out = no_rhs_buffer_allocation_status
    no_rhs_injection_status_out = no_rhs_injection_status
    no_production_hook_status_out = no_production_hook_status
    no_fluid_field_access_status_out = no_fluid_field_access_status
    no_fluid_field_modification_status_out = no_fluid_field_modification_status
    no_pressure_modification_status_out = no_pressure_modification_status
    no_projection_modification_status_out = no_projection_modification_status
    no_rk3_modification_status_out = no_rk3_modification_status
    no_channel_forcing_modification_status_out = no_channel_forcing_modification_status
    no_production_ibm_forcing_status_out = no_production_ibm_forcing_status
    no_feedback_application_status_out = no_feedback_application_status
    no_twoway_force_status_out = no_twoway_force_status
    no_structure_advance_status_out = no_structure_advance_status
    config_status_out = config_status
  end subroutine stage14_get_status_values

  integer function logical_to_int(value)
    logical, intent(in) :: value
    if (value) then
      logical_to_int = 1
    else
      logical_to_int = 0
    end if
  end function logical_to_int

  logical function finite_real(value)
    real(mytype), intent(in) :: value
    finite_real = (value == value) .and. (abs(value) < huge(value))
  end function finite_real

  subroutine parse_bool(text, value, parse_status)
    character(len=*), intent(in) :: text
    logical, intent(out) :: value
    integer, intent(out) :: parse_status
    character(len=256) :: lowered

    lowered = to_lower(trim(adjustl(text)))
    select case (trim(lowered))
    case ('1', 'true', 'yes')
      value = .true.
      parse_status = 1
    case ('0', 'false', 'no')
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

end module fibre_stage14_config
