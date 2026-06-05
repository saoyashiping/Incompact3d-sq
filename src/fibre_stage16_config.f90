module fibre_stage16_config
  use decomp_2d_constants, only : mytype
  implicit none

  private

  public :: stage16_config_reset
  public :: stage16_config_load_from_environment
  public :: stage16_config_validate
  public :: stage16_is_enabled
  public :: stage16_one_fibre_fsi_requested
  public :: stage16_structure_advance_requested
  public :: stage16_two_way_rhs_requested
  public :: stage16_diagnostic_only_requested
  public :: stage16_get_lambda
  public :: stage16_get_feedback_alpha
  public :: stage16_get_max_structure_update
  public :: stage16_get_max_rhs_increment
  public :: stage16_config_write_diagnostics
  public :: stage16_config_set_for_test
  public :: stage16_config_get_status_values

  logical :: config_loaded = .false.
  logical :: stage16_enable_value = .false.
  logical :: one_fibre_fsi_enable_value = .false.
  logical :: structure_advance_enable_value = .false.
  logical :: two_way_rhs_enable_value = .false.
  logical :: diagnostic_only_value = .true.
  logical :: require_stage15_closed_value = .true.

  real(mytype) :: feedback_alpha_value = 1.0_mytype
  real(mytype) :: lambda_value = 0.0_mytype
  real(mytype) :: max_structure_update_value = 1.0e-12_mytype
  real(mytype) :: max_rhs_increment_value = 1.0e-8_mytype

  integer :: master_enable_status = 0
  integer :: one_fibre_fsi_enable_status = 0
  integer :: structure_advance_enable_status = 0
  integer :: two_way_rhs_enable_status = 0
  integer :: diagnostic_only_status = 1
  integer :: require_stage15_closed_status = 1
  integer :: default_off_status = 1
  integer :: bool_parse_status = 1
  integer :: numeric_parse_status = 1
  integer :: numeric_bounds_status = 1
  integer :: diagnostic_only_guard_status = 1
  integer :: invalid_flag_combination_status = 1
  integer :: no_production_hook_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_rhs_modification_status = 1
  integer :: no_bending_solve_status = 1
  integer :: no_tension_solve_status = 1
  integer :: no_wall_contact_status = 1
  integer :: no_multifibre_status = 1
  integer :: no_pressure_projection_modification_status = 1
  integer :: no_poisson_modification_status = 1
  integer :: no_rk3_channel_forcing_modification_status = 1
  integer :: no_legacy_ibm_forcing_status = 1
  integer :: config_status = 1

contains

  subroutine stage16_config_reset()
    config_loaded = .false.
    stage16_enable_value = .false.
    one_fibre_fsi_enable_value = .false.
    structure_advance_enable_value = .false.
    two_way_rhs_enable_value = .false.
    diagnostic_only_value = .true.
    require_stage15_closed_value = .true.
    feedback_alpha_value = 1.0_mytype
    lambda_value = 0.0_mytype
    max_structure_update_value = 1.0e-12_mytype
    max_rhs_increment_value = 1.0e-8_mytype

    master_enable_status = 0
    one_fibre_fsi_enable_status = 0
    structure_advance_enable_status = 0
    two_way_rhs_enable_status = 0
    diagnostic_only_status = 1
    require_stage15_closed_status = 1
    default_off_status = 1
    bool_parse_status = 1
    numeric_parse_status = 1
    numeric_bounds_status = 1
    diagnostic_only_guard_status = 1
    invalid_flag_combination_status = 1
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
    config_status = 1
  end subroutine stage16_config_reset

  subroutine stage16_config_load_from_environment()
    character(len=256) :: env_value
    integer :: env_status
    logical :: parsed_bool
    integer :: parsed_status
    real(mytype) :: parsed_real

    call stage16_config_reset()
    config_loaded = .true.

    call load_bool_env('STAGE16_1_ENABLE', stage16_enable_value, parsed_status)
    bool_parse_status = min(bool_parse_status, parsed_status)
    call load_bool_env('STAGE16_1_ONE_FIBRE_FSI_ENABLE', one_fibre_fsi_enable_value, parsed_status)
    bool_parse_status = min(bool_parse_status, parsed_status)
    call load_bool_env('STAGE16_1_STRUCTURE_ADVANCE_ENABLE', structure_advance_enable_value, parsed_status)
    bool_parse_status = min(bool_parse_status, parsed_status)
    call load_bool_env('STAGE16_1_TWO_WAY_RHS_ENABLE', two_way_rhs_enable_value, parsed_status)
    bool_parse_status = min(bool_parse_status, parsed_status)
    call load_bool_env('STAGE16_1_DIAGNOSTIC_ONLY', diagnostic_only_value, parsed_status)
    bool_parse_status = min(bool_parse_status, parsed_status)
    call load_bool_env('STAGE16_1_REQUIRE_STAGE15_CLOSED', require_stage15_closed_value, parsed_status)
    bool_parse_status = min(bool_parse_status, parsed_status)

    call get_environment_variable('STAGE16_1_FEEDBACK_ALPHA', value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=parsed_status) parsed_real
      if (parsed_status == 0 .and. finite_real(parsed_real)) then
        feedback_alpha_value = parsed_real
      else
        numeric_parse_status = 0
      end if
    end if

    call get_environment_variable('STAGE16_1_LAMBDA', value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=parsed_status) parsed_real
      if (parsed_status == 0 .and. finite_real(parsed_real)) then
        lambda_value = parsed_real
      else
        numeric_parse_status = 0
      end if
    end if

    call get_environment_variable('STAGE16_1_MAX_STRUCTURE_UPDATE', value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=parsed_status) parsed_real
      if (parsed_status == 0 .and. finite_real(parsed_real)) then
        max_structure_update_value = parsed_real
      else
        numeric_parse_status = 0
      end if
    end if

    call get_environment_variable('STAGE16_1_MAX_RHS_INCREMENT', value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=parsed_status) parsed_real
      if (parsed_status == 0 .and. finite_real(parsed_real)) then
        max_rhs_increment_value = parsed_real
      else
        numeric_parse_status = 0
      end if
    end if

    call stage16_config_validate()
  contains
    subroutine load_bool_env(name, value, parse_status)
      character(len=*), intent(in) :: name
      logical, intent(inout) :: value
      integer, intent(out) :: parse_status
      character(len=256) :: local_env
      integer :: local_status
      logical :: local_bool

      parse_status = 1
      call get_environment_variable(name, value=local_env, status=local_status)
      if (local_status == 0) then
        call parse_bool(local_env, local_bool, parse_status)
        if (parse_status == 1) value = local_bool
      end if
    end subroutine load_bool_env
  end subroutine stage16_config_load_from_environment

  subroutine stage16_config_set_for_test(master_enable, one_fibre_fsi_enable, structure_advance_enable, &
       two_way_rhs_enable, diagnostic_only, feedback_alpha, lambda_in, max_structure_update, max_rhs_increment)
    logical, intent(in) :: master_enable
    logical, intent(in) :: one_fibre_fsi_enable
    logical, intent(in) :: structure_advance_enable
    logical, intent(in) :: two_way_rhs_enable
    logical, intent(in) :: diagnostic_only
    real(mytype), intent(in) :: feedback_alpha
    real(mytype), intent(in) :: lambda_in
    real(mytype), intent(in) :: max_structure_update
    real(mytype), intent(in) :: max_rhs_increment

    call stage16_config_reset()
    config_loaded = .true.
    stage16_enable_value = master_enable
    one_fibre_fsi_enable_value = one_fibre_fsi_enable
    structure_advance_enable_value = structure_advance_enable
    two_way_rhs_enable_value = two_way_rhs_enable
    diagnostic_only_value = diagnostic_only
    feedback_alpha_value = feedback_alpha
    lambda_value = lambda_in
    max_structure_update_value = max_structure_update
    max_rhs_increment_value = max_rhs_increment
    call stage16_config_validate()
  end subroutine stage16_config_set_for_test

  subroutine stage16_config_validate()
    default_off_status = logical_to_int(.not. stage16_enable_value .and. .not. one_fibre_fsi_enable_value .and. &
         .not. structure_advance_enable_value .and. .not. two_way_rhs_enable_value)
    master_enable_status = logical_to_int(stage16_enable_value)
    one_fibre_fsi_enable_status = logical_to_int(one_fibre_fsi_enable_value)
    structure_advance_enable_status = logical_to_int(structure_advance_enable_value)
    two_way_rhs_enable_status = logical_to_int(two_way_rhs_enable_value)
    diagnostic_only_status = logical_to_int(diagnostic_only_value)
    require_stage15_closed_status = logical_to_int(require_stage15_closed_value)

    numeric_bounds_status = logical_to_int(finite_real(feedback_alpha_value) .and. finite_real(lambda_value) .and. &
         finite_real(max_structure_update_value) .and. finite_real(max_rhs_increment_value) .and. &
         feedback_alpha_value >= 0.0_mytype .and. lambda_value >= 0.0_mytype .and. &
         max_structure_update_value >= 0.0_mytype .and. max_rhs_increment_value >= 0.0_mytype)

    diagnostic_only_guard_status = logical_to_int(diagnostic_only_value)
    invalid_flag_combination_status = logical_to_int(.not. ((one_fibre_fsi_enable_value .and. .not. diagnostic_only_value) .or. &
         structure_advance_enable_value .or. two_way_rhs_enable_value))

    no_production_hook_status = 1
    no_structure_advance_status = logical_to_int(.not. structure_advance_enable_value)
    no_rhs_modification_status = logical_to_int(.not. two_way_rhs_enable_value)
    no_bending_solve_status = 1
    no_tension_solve_status = 1
    no_wall_contact_status = 1
    no_multifibre_status = 1
    no_pressure_projection_modification_status = 1
    no_poisson_modification_status = 1
    no_rk3_channel_forcing_modification_status = 1
    no_legacy_ibm_forcing_status = 1

    config_status = logical_to_int(bool_parse_status == 1 .and. numeric_parse_status == 1 .and. &
         numeric_bounds_status == 1 .and. diagnostic_only_guard_status == 1 .and. &
         invalid_flag_combination_status == 1)
  end subroutine stage16_config_validate

  logical function stage16_is_enabled()
    if (.not. config_loaded) call stage16_config_load_from_environment()
    stage16_is_enabled = stage16_enable_value
  end function stage16_is_enabled

  logical function stage16_one_fibre_fsi_requested()
    if (.not. config_loaded) call stage16_config_load_from_environment()
    stage16_one_fibre_fsi_requested = one_fibre_fsi_enable_value
  end function stage16_one_fibre_fsi_requested

  logical function stage16_structure_advance_requested()
    if (.not. config_loaded) call stage16_config_load_from_environment()
    stage16_structure_advance_requested = structure_advance_enable_value
  end function stage16_structure_advance_requested

  logical function stage16_two_way_rhs_requested()
    if (.not. config_loaded) call stage16_config_load_from_environment()
    stage16_two_way_rhs_requested = two_way_rhs_enable_value
  end function stage16_two_way_rhs_requested

  logical function stage16_diagnostic_only_requested()
    if (.not. config_loaded) call stage16_config_load_from_environment()
    stage16_diagnostic_only_requested = diagnostic_only_value
  end function stage16_diagnostic_only_requested

  real(mytype) function stage16_get_lambda()
    if (.not. config_loaded) call stage16_config_load_from_environment()
    stage16_get_lambda = lambda_value
  end function stage16_get_lambda

  real(mytype) function stage16_get_feedback_alpha()
    if (.not. config_loaded) call stage16_config_load_from_environment()
    stage16_get_feedback_alpha = feedback_alpha_value
  end function stage16_get_feedback_alpha

  real(mytype) function stage16_get_max_structure_update()
    if (.not. config_loaded) call stage16_config_load_from_environment()
    stage16_get_max_structure_update = max_structure_update_value
  end function stage16_get_max_structure_update

  real(mytype) function stage16_get_max_rhs_increment()
    if (.not. config_loaded) call stage16_config_load_from_environment()
    stage16_get_max_rhs_increment = max_rhs_increment_value
  end function stage16_get_max_rhs_increment

  subroutine stage16_config_get_status_values(default_off_out, env_config_out, config_out, invalid_flag_out, &
       numeric_parse_out, numeric_bounds_out, no_structure_advance_out, no_rhs_modification_out)
    integer, intent(out) :: default_off_out
    integer, intent(out) :: env_config_out
    integer, intent(out) :: config_out
    integer, intent(out) :: invalid_flag_out
    integer, intent(out) :: numeric_parse_out
    integer, intent(out) :: numeric_bounds_out
    integer, intent(out) :: no_structure_advance_out
    integer, intent(out) :: no_rhs_modification_out

    default_off_out = default_off_status
    env_config_out = logical_to_int(config_loaded)
    config_out = config_status
    invalid_flag_out = invalid_flag_combination_status
    numeric_parse_out = numeric_parse_status
    numeric_bounds_out = numeric_bounds_status
    no_structure_advance_out = no_structure_advance_status
    no_rhs_modification_out = no_rhs_modification_status
  end subroutine stage16_config_get_status_values

  subroutine stage16_config_write_diagnostics(unit)
    integer, intent(in) :: unit

    write(unit,'(A,1X,I0)') 'master_enable_status', master_enable_status
    write(unit,'(A,1X,I0)') 'one_fibre_fsi_enable_status', one_fibre_fsi_enable_status
    write(unit,'(A,1X,I0)') 'structure_advance_enable_status', structure_advance_enable_status
    write(unit,'(A,1X,I0)') 'two_way_rhs_enable_status', two_way_rhs_enable_status
    write(unit,'(A,1X,I0)') 'diagnostic_only_status', diagnostic_only_status
    write(unit,'(A,1X,I0)') 'require_stage15_closed_status', require_stage15_closed_status
    write(unit,'(A,1X,ES24.16E3)') 'feedback_alpha', feedback_alpha_value
    write(unit,'(A,1X,ES24.16E3)') 'lambda_value', lambda_value
    write(unit,'(A,1X,ES24.16E3)') 'max_structure_update', max_structure_update_value
    write(unit,'(A,1X,ES24.16E3)') 'max_rhs_increment', max_rhs_increment_value
    write(unit,'(A,1X,I0)') 'bool_parse_status', bool_parse_status
    write(unit,'(A,1X,I0)') 'numeric_parse_status', numeric_parse_status
    write(unit,'(A,1X,I0)') 'numeric_bounds_status', numeric_bounds_status
    write(unit,'(A,1X,I0)') 'diagnostic_only_guard_status', diagnostic_only_guard_status
    write(unit,'(A,1X,I0)') 'invalid_flag_combination_status', invalid_flag_combination_status
    write(unit,'(A,1X,I0)') 'no_production_hook_status', no_production_hook_status
    write(unit,'(A,1X,I0)') 'no_structure_advance_status', no_structure_advance_status
    write(unit,'(A,1X,I0)') 'no_rhs_modification_status', no_rhs_modification_status
    write(unit,'(A,1X,I0)') 'no_bending_solve_status', no_bending_solve_status
    write(unit,'(A,1X,I0)') 'no_tension_solve_status', no_tension_solve_status
    write(unit,'(A,1X,I0)') 'no_wall_contact_status', no_wall_contact_status
    write(unit,'(A,1X,I0)') 'no_multifibre_status', no_multifibre_status
    write(unit,'(A,1X,I0)') 'no_pressure_projection_modification_status', no_pressure_projection_modification_status
    write(unit,'(A,1X,I0)') 'no_poisson_modification_status', no_poisson_modification_status
    write(unit,'(A,1X,I0)') 'no_rk3_channel_forcing_modification_status', no_rk3_channel_forcing_modification_status
    write(unit,'(A,1X,I0)') 'no_legacy_ibm_forcing_status', no_legacy_ibm_forcing_status
    write(unit,'(A,1X,I0)') 'config_status', config_status
  end subroutine stage16_config_write_diagnostics

  subroutine parse_bool(value, parsed_bool, parse_status)
    character(len=*), intent(in) :: value
    logical, intent(out) :: parsed_bool
    integer, intent(out) :: parse_status
    character(len=256) :: lowered

    lowered = lower_string(trim(adjustl(value)))
    select case (trim(lowered))
    case ('1', 'true', 't', 'yes', 'y', 'on')
      parsed_bool = .true.
      parse_status = 1
    case ('0', 'false', 'f', 'no', 'n', 'off')
      parsed_bool = .false.
      parse_status = 1
    case default
      parsed_bool = .false.
      parse_status = 0
    end select
  end subroutine parse_bool

  character(len=256) function lower_string(value)
    character(len=*), intent(in) :: value
    integer :: i
    integer :: code

    lower_string = ' '
    do i = 1, min(len_trim(value), len(lower_string))
      code = iachar(value(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower_string(i:i) = achar(code + 32)
      else
        lower_string(i:i) = value(i:i)
      end if
    end do
  end function lower_string

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

end module fibre_stage16_config
