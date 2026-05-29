module fibre_stage12_config
  implicit none
  private

  integer, parameter :: default_max_points = 8
  real(kind=kind(1.0d0)), parameter :: default_feedback_gain = 1.0d0
  integer, parameter :: default_force_sign = 1

  logical :: loaded = .false.
  logical :: requested = .false.
  logical :: readonly_mode = .true.
  real(kind=kind(1.0d0)) :: feedback_gain = default_feedback_gain
  integer :: force_sign = default_force_sign
  integer :: max_points = default_max_points

  integer :: requested_flag = 0
  integer :: readonly_mode_status = 1
  integer :: disabled_by_default_status = 1
  integer :: feedback_gain_status = 1
  integer :: force_sign_status = 1
  integer :: no_force_computation_status = 1
  integer :: no_eulerian_force_density_status = 1
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: config_status = 1

  public :: stage12_config_load
  public :: stage12_requested
  public :: stage12_readonly_mode
  public :: stage12_get_feedback_gain
  public :: stage12_get_force_sign
  public :: stage12_get_max_points
  public :: stage12_get_status_values

contains

  subroutine stage12_config_load()
    character(len=64) :: value
    integer :: status

    requested = .false.
    readonly_mode = .true.
    feedback_gain = default_feedback_gain
    force_sign = default_force_sign
    max_points = default_max_points

    requested_flag = 0
    readonly_mode_status = 1
    disabled_by_default_status = 1
    feedback_gain_status = 1
    force_sign_status = 1
    no_force_computation_status = 1
    no_eulerian_force_density_status = 1
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    no_fluid_field_access_status = 1
    no_fluid_field_modification_status = 1
    config_status = 1

    call get_environment_variable('X3D_STAGE12_FEEDBACK_CANDIDATE', value=value, status=status)
    if (status == 0) then
      if (trim(adjustl(value)) == '1') then
        requested = .true.
        requested_flag = 1
        disabled_by_default_status = 0
      else
        config_status = 0
      end if
    end if

    call get_environment_variable('X3D_STAGE12_FORCE_READONLY', value=value, status=status)
    if (status == 0) then
      if (trim(adjustl(value)) /= '1') then
        config_status = 0
      end if
    end if

    call parse_positive_real_env('X3D_STAGE12_FEEDBACK_GAIN', default_feedback_gain, feedback_gain, feedback_gain_status)
    call parse_force_sign_env('X3D_STAGE12_FORCE_SIGN', default_force_sign, force_sign, force_sign_status)
    call parse_positive_int_env('X3D_STAGE12_MAX_POINTS', default_max_points, max_points, config_status)

    if (feedback_gain_status /= 1) config_status = 0
    if (force_sign_status /= 1) config_status = 0

    readonly_mode = .true.
    readonly_mode_status = 1

    loaded = .true.
  end subroutine stage12_config_load

  subroutine parse_positive_real_env(name, default_value, parsed_value, parse_status)
    character(len=*), intent(in) :: name
    real(kind=kind(1.0d0)), intent(in) :: default_value
    real(kind=kind(1.0d0)), intent(out) :: parsed_value
    integer, intent(out) :: parse_status
    character(len=64) :: value
    integer :: status
    integer :: ios
    real(kind=kind(1.0d0)) :: tmp

    parsed_value = default_value
    parse_status = 1
    call get_environment_variable(name, value=value, status=status)
    if (status == 0) then
      read(value, *, iostat=ios) tmp
      if (ios == 0 .and. tmp > 0.0d0) then
        parsed_value = tmp
      else
        parsed_value = default_value
        parse_status = 0
      end if
    end if
  end subroutine parse_positive_real_env

  subroutine parse_force_sign_env(name, default_value, parsed_value, parse_status)
    character(len=*), intent(in) :: name
    integer, intent(in) :: default_value
    integer, intent(out) :: parsed_value
    integer, intent(out) :: parse_status
    character(len=64) :: value
    character(len=64) :: trimmed_value
    integer :: status

    parsed_value = default_value
    parse_status = 1
    call get_environment_variable(name, value=value, status=status)
    if (status == 0) then
      trimmed_value = trim(adjustl(value))
      if (trimmed_value == '+1' .or. trimmed_value == '1') then
        parsed_value = 1
      else if (trimmed_value == '-1') then
        parsed_value = -1
      else
        parsed_value = default_value
        parse_status = 0
      end if
    end if
  end subroutine parse_force_sign_env

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
        parsed_value = default_value
        status_accum = 0
      end if
    end if
  end subroutine parse_positive_int_env

  logical function stage12_requested()
    if (.not. loaded) call stage12_config_load()
    stage12_requested = requested
  end function stage12_requested

  logical function stage12_readonly_mode()
    if (.not. loaded) call stage12_config_load()
    stage12_readonly_mode = readonly_mode
  end function stage12_readonly_mode

  real(kind=kind(1.0d0)) function stage12_get_feedback_gain()
    if (.not. loaded) call stage12_config_load()
    stage12_get_feedback_gain = feedback_gain
  end function stage12_get_feedback_gain

  integer function stage12_get_force_sign()
    if (.not. loaded) call stage12_config_load()
    stage12_get_force_sign = force_sign
  end function stage12_get_force_sign

  integer function stage12_get_max_points()
    if (.not. loaded) call stage12_config_load()
    stage12_get_max_points = max_points
  end function stage12_get_max_points

  subroutine stage12_get_status_values(out_requested_flag, out_readonly_mode_status, &
                                       out_disabled_by_default_status, out_feedback_gain_status, &
                                       out_force_sign_status, out_no_force_computation_status, &
                                       out_no_eulerian_force_density_status, out_no_rhs_injection_status, &
                                       out_no_ibm_spreading_status, out_no_feedback_application_status, &
                                       out_no_twoway_force_status, out_no_structure_advance_status, &
                                       out_no_fluid_field_access_status, out_no_fluid_field_modification_status, &
                                       out_config_status)
    integer, intent(out) :: out_requested_flag
    integer, intent(out) :: out_readonly_mode_status
    integer, intent(out) :: out_disabled_by_default_status
    integer, intent(out) :: out_feedback_gain_status
    integer, intent(out) :: out_force_sign_status
    integer, intent(out) :: out_no_force_computation_status
    integer, intent(out) :: out_no_eulerian_force_density_status
    integer, intent(out) :: out_no_rhs_injection_status
    integer, intent(out) :: out_no_ibm_spreading_status
    integer, intent(out) :: out_no_feedback_application_status
    integer, intent(out) :: out_no_twoway_force_status
    integer, intent(out) :: out_no_structure_advance_status
    integer, intent(out) :: out_no_fluid_field_access_status
    integer, intent(out) :: out_no_fluid_field_modification_status
    integer, intent(out) :: out_config_status

    if (.not. loaded) call stage12_config_load()

    out_requested_flag = requested_flag
    out_readonly_mode_status = readonly_mode_status
    out_disabled_by_default_status = disabled_by_default_status
    out_feedback_gain_status = feedback_gain_status
    out_force_sign_status = force_sign_status
    out_no_force_computation_status = no_force_computation_status
    out_no_eulerian_force_density_status = no_eulerian_force_density_status
    out_no_rhs_injection_status = no_rhs_injection_status
    out_no_ibm_spreading_status = no_ibm_spreading_status
    out_no_feedback_application_status = no_feedback_application_status
    out_no_twoway_force_status = no_twoway_force_status
    out_no_structure_advance_status = no_structure_advance_status
    out_no_fluid_field_access_status = no_fluid_field_access_status
    out_no_fluid_field_modification_status = no_fluid_field_modification_status
    out_config_status = config_status
  end subroutine stage12_get_status_values

end module fibre_stage12_config
