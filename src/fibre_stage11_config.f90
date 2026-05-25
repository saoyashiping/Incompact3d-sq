module fibre_stage11_config
  implicit none
  private

  logical, save :: loaded = .false.
  logical, save :: requested = .false.
  logical, save :: readonly_mode = .true.
  integer, save :: max_points = 8
  integer, save :: max_steps = 3

  integer, save :: requested_flag = 0
  integer, save :: readonly_mode_status = 1
  integer, save :: disabled_by_default_status = 1
  integer, save :: no_lagrangian_state_status = 1
  integer, save :: no_velocity_sampling_status = 1
  integer, save :: no_fluid_field_access_status = 1
  integer, save :: no_fluid_field_modification_status = 1
  integer, save :: no_rhs_injection_status = 1
  integer, save :: no_ibm_spreading_status = 1
  integer, save :: no_feedback_force_status = 1
  integer, save :: no_twoway_force_status = 1
  integer, save :: no_structure_advance_status = 1
  integer, save :: config_status = 1

  public :: stage11_config_load
  public :: stage11_requested
  public :: stage11_readonly_mode
  public :: stage11_get_max_points
  public :: stage11_get_max_steps
  public :: stage11_get_status_values

contains

  subroutine stage11_config_load()
    character(len=64) :: env_value
    integer :: env_status, io_status
    integer :: parsed_points, parsed_steps

    requested = .false.
    readonly_mode = .true.
    max_points = 8
    max_steps = 3
    config_status = 1

    call get_environment_variable('X3D_STAGE11_ONEWAY_HOOK', value=env_value, status=env_status)
    if (env_status == 0) then
      requested = (trim(adjustl(env_value)) == '1')
    endif

    call get_environment_variable('X3D_STAGE11_FORCE_READONLY', value=env_value, status=env_status)
    if (env_status == 0) then
      readonly_mode = (trim(adjustl(env_value)) == '1')
    else
      readonly_mode = .true.
    endif

    call get_environment_variable('X3D_STAGE11_MAX_POINTS', value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=io_status) parsed_points
      if (io_status == 0 .and. parsed_points >= 0) then
        max_points = parsed_points
      else
        config_status = 0
      endif
    endif

    call get_environment_variable('X3D_STAGE11_MAX_STEPS', value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=io_status) parsed_steps
      if (io_status == 0 .and. parsed_steps >= 1) then
        max_steps = parsed_steps
      else
        config_status = 0
      endif
    endif

    requested_flag = merge(1, 0, requested)
    readonly_mode_status = merge(1, 0, readonly_mode)
    disabled_by_default_status = merge(1, 0, .not. requested)

    loaded = .true.
  end subroutine stage11_config_load

  logical function stage11_requested()
    if (.not. loaded) call stage11_config_load()
    stage11_requested = requested
  end function stage11_requested

  logical function stage11_readonly_mode()
    if (.not. loaded) call stage11_config_load()
    stage11_readonly_mode = readonly_mode
  end function stage11_readonly_mode

  integer function stage11_get_max_points()
    if (.not. loaded) call stage11_config_load()
    stage11_get_max_points = max_points
  end function stage11_get_max_points

  integer function stage11_get_max_steps()
    if (.not. loaded) call stage11_config_load()
    stage11_get_max_steps = max_steps
  end function stage11_get_max_steps

  subroutine stage11_get_status_values(o_requested_flag, o_readonly_mode_status, o_disabled_by_default_status, &
                                       o_no_lagrangian_state_status, o_no_velocity_sampling_status, &
                                       o_no_fluid_field_access_status, o_no_fluid_field_modification_status, &
                                       o_no_rhs_injection_status, o_no_ibm_spreading_status, &
                                       o_no_feedback_force_status, o_no_twoway_force_status, &
                                       o_no_structure_advance_status, o_config_status)
    integer, intent(out) :: o_requested_flag, o_readonly_mode_status, o_disabled_by_default_status
    integer, intent(out) :: o_no_lagrangian_state_status, o_no_velocity_sampling_status
    integer, intent(out) :: o_no_fluid_field_access_status, o_no_fluid_field_modification_status
    integer, intent(out) :: o_no_rhs_injection_status, o_no_ibm_spreading_status
    integer, intent(out) :: o_no_feedback_force_status, o_no_twoway_force_status
    integer, intent(out) :: o_no_structure_advance_status, o_config_status

    if (.not. loaded) call stage11_config_load()

    o_requested_flag = requested_flag
    o_readonly_mode_status = readonly_mode_status
    o_disabled_by_default_status = disabled_by_default_status
    o_no_lagrangian_state_status = no_lagrangian_state_status
    o_no_velocity_sampling_status = no_velocity_sampling_status
    o_no_fluid_field_access_status = no_fluid_field_access_status
    o_no_fluid_field_modification_status = no_fluid_field_modification_status
    o_no_rhs_injection_status = no_rhs_injection_status
    o_no_ibm_spreading_status = no_ibm_spreading_status
    o_no_feedback_force_status = no_feedback_force_status
    o_no_twoway_force_status = no_twoway_force_status
    o_no_structure_advance_status = no_structure_advance_status
    o_config_status = config_status
  end subroutine stage11_get_status_values

end module fibre_stage11_config
