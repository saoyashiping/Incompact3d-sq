module fibre_stage11_config
  implicit none
  private

  integer, parameter :: default_max_points = 8
  integer, parameter :: default_max_steps = 3

  logical :: loaded = .false.
  logical :: requested = .false.
  logical :: readonly_mode = .true.
  integer :: max_points = default_max_points
  integer :: max_steps = default_max_steps

  integer :: requested_flag = 0
  integer :: readonly_mode_status = 1
  integer :: disabled_by_default_status = 1
  integer :: no_lagrangian_state_status = 1
  integer :: no_velocity_sampling_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_force_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: config_status = 1

  public :: stage11_config_load
  public :: stage11_requested
  public :: stage11_readonly_mode
  public :: stage11_get_max_points
  public :: stage11_get_max_steps
  public :: stage11_get_status_values

contains

  subroutine stage11_config_load()
    character(len=64) :: value
    integer :: status

    requested = .false.
    readonly_mode = .true.
    max_points = default_max_points
    max_steps = default_max_steps

    requested_flag = 0
    readonly_mode_status = 1
    disabled_by_default_status = 1
    no_lagrangian_state_status = 1
    no_velocity_sampling_status = 1
    no_fluid_field_access_status = 1
    no_fluid_field_modification_status = 1
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_force_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    config_status = 1

    call get_environment_variable('X3D_STAGE11_ONEWAY_HOOK', value=value, status=status)
    if (status == 0) then
      if (trim(adjustl(value)) == '1') then
        requested = .true.
        requested_flag = 1
        disabled_by_default_status = 0
      else
        config_status = 0
      end if
    end if

    call get_environment_variable('X3D_STAGE11_FORCE_READONLY', value=value, status=status)
    if (status == 0) then
      if (trim(adjustl(value)) /= '1') then
        config_status = 0
      end if
    end if

    call parse_positive_int_env('X3D_STAGE11_MAX_POINTS', default_max_points, max_points, config_status)
    call parse_positive_int_env('X3D_STAGE11_MAX_STEPS', default_max_steps, max_steps, config_status)

    readonly_mode = .true.
    readonly_mode_status = 1

    loaded = .true.
  end subroutine stage11_config_load

  subroutine parse_positive_int_env(name, default_value, parsed_value, status_accum)
    character(len=*), intent(in) :: name
    integer, intent(in) :: default_value
    integer, intent(out) :: parsed_value
    integer, intent(inout) :: status_accum
    character(len=64) :: value
    integer :: status
    integer :: ios
    integer :: tmp

    parsed_value = default_value
    call get_environment_variable(name, value=value, status=status)
    if (status == 0) then
      read(value, *, iostat=ios) tmp
      if (ios == 0 .and. tmp > 0) then
        parsed_value = tmp
      else
        status_accum = 0
      end if
    end if
  end subroutine parse_positive_int_env

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

  subroutine stage11_get_status_values(out_requested_flag, out_readonly_mode_status, &
                                       out_disabled_by_default_status, out_no_lagrangian_state_status, &
                                       out_no_velocity_sampling_status, out_no_fluid_field_access_status, &
                                       out_no_fluid_field_modification_status, out_no_rhs_injection_status, &
                                       out_no_ibm_spreading_status, out_no_feedback_force_status, &
                                       out_no_twoway_force_status, out_no_structure_advance_status, out_config_status)
    integer, intent(out) :: out_requested_flag
    integer, intent(out) :: out_readonly_mode_status
    integer, intent(out) :: out_disabled_by_default_status
    integer, intent(out) :: out_no_lagrangian_state_status
    integer, intent(out) :: out_no_velocity_sampling_status
    integer, intent(out) :: out_no_fluid_field_access_status
    integer, intent(out) :: out_no_fluid_field_modification_status
    integer, intent(out) :: out_no_rhs_injection_status
    integer, intent(out) :: out_no_ibm_spreading_status
    integer, intent(out) :: out_no_feedback_force_status
    integer, intent(out) :: out_no_twoway_force_status
    integer, intent(out) :: out_no_structure_advance_status
    integer, intent(out) :: out_config_status

    if (.not. loaded) call stage11_config_load()

    out_requested_flag = requested_flag
    out_readonly_mode_status = readonly_mode_status
    out_disabled_by_default_status = disabled_by_default_status
    out_no_lagrangian_state_status = no_lagrangian_state_status
    out_no_velocity_sampling_status = no_velocity_sampling_status
    out_no_fluid_field_access_status = no_fluid_field_access_status
    out_no_fluid_field_modification_status = no_fluid_field_modification_status
    out_no_rhs_injection_status = no_rhs_injection_status
    out_no_ibm_spreading_status = no_ibm_spreading_status
    out_no_feedback_force_status = no_feedback_force_status
    out_no_twoway_force_status = no_twoway_force_status
    out_no_structure_advance_status = no_structure_advance_status
    out_config_status = config_status
  end subroutine stage11_get_status_values

end module fibre_stage11_config
