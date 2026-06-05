module fibre_stage16_force_sign_audit
  use decomp_2d_constants, only : mytype
  implicit none

  private

  public :: stage16_force_sign_audit_reset
  public :: stage16_force_sign_audit_load_from_environment
  public :: stage16_force_sign_audit_compute_reference
  public :: stage16_force_sign_audit_validate_action_reaction
  public :: stage16_force_sign_audit_validate_wrong_sign_rejection
  public :: stage16_force_sign_audit_write_diagnostics
  public :: stage16_force_sign_audit_get_final_status

  real(mytype) :: feedback_alpha_value = 1.0_mytype
  real(mytype) :: max_action_reaction_error = 1.0e-14_mytype
  real(mytype) :: max_sign_error = 1.0e-14_mytype
  real(mytype) :: slip(3) = 0.0_mytype
  real(mytype) :: feedback_force(3) = 0.0_mytype
  real(mytype) :: fluid_on_fibre(3) = 0.0_mytype
  real(mytype) :: fibre_on_fluid(3) = 0.0_mytype
  real(mytype) :: reversed_fluid_on_fibre(3) = 0.0_mytype
  real(mytype) :: reversed_fibre_on_fluid(3) = 0.0_mytype
  real(mytype) :: action_reaction_error_value = 0.0_mytype
  real(mytype) :: wrong_sign_error_value = 0.0_mytype
  real(mytype) :: zero_slip_force_value = 0.0_mytype

  integer :: stage16_3_requested_status = 1
  integer :: slip_finite_status = 0
  integer :: feedback_force_finite_status = 0
  integer :: structure_side_force_sign_status = 0
  integer :: fluid_side_force_sign_status = 0
  integer :: action_reaction_status = 0
  integer :: wrong_sign_rejection_status = 0
  integer :: zero_slip_zero_force_status = 0
  integer :: slip_reversal_sign_status = 0
  integer :: approved_stage12_13_14_chain_status = 1
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
  integer :: numeric_parse_status = 1
  integer :: numeric_bounds_status = 1
  integer :: diagnostic_only_status = 1
  integer :: audit_status = 0

contains

  subroutine stage16_force_sign_audit_reset()
    feedback_alpha_value = 1.0_mytype
    max_action_reaction_error = 1.0e-14_mytype
    max_sign_error = 1.0e-14_mytype
    slip(:) = 0.0_mytype
    feedback_force(:) = 0.0_mytype
    fluid_on_fibre(:) = 0.0_mytype
    fibre_on_fluid(:) = 0.0_mytype
    reversed_fluid_on_fibre(:) = 0.0_mytype
    reversed_fibre_on_fluid(:) = 0.0_mytype
    action_reaction_error_value = 0.0_mytype
    wrong_sign_error_value = 0.0_mytype
    zero_slip_force_value = 0.0_mytype
    stage16_3_requested_status = 1
    slip_finite_status = 0
    feedback_force_finite_status = 0
    structure_side_force_sign_status = 0
    fluid_side_force_sign_status = 0
    action_reaction_status = 0
    wrong_sign_rejection_status = 0
    zero_slip_zero_force_status = 0
    slip_reversal_sign_status = 0
    approved_stage12_13_14_chain_status = 1
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
    numeric_parse_status = 1
    numeric_bounds_status = 1
    diagnostic_only_status = 1
    audit_status = 0
  end subroutine stage16_force_sign_audit_reset

  subroutine stage16_force_sign_audit_load_from_environment()
    character(len=256) :: env_value
    integer :: env_status
    integer :: parse_status
    real(mytype) :: parsed_real
    logical :: parsed_bool

    call stage16_force_sign_audit_reset()
    call load_bool_env('STAGE16_3_ENABLE', parsed_bool, parse_status)
    if (parse_status == 1) stage16_3_requested_status = logical_to_int(parsed_bool)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_bool_env('STAGE16_3_ONE_FIBRE_FSI_ENABLE', parsed_bool, parse_status)
    if (parse_status == 1) stage16_3_requested_status = min(stage16_3_requested_status, logical_to_int(parsed_bool))
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_bool_env('STAGE16_3_DIAGNOSTIC_ONLY', parsed_bool, parse_status)
    if (parse_status == 1) diagnostic_only_status = logical_to_int(parsed_bool)
    numeric_parse_status = min(numeric_parse_status, parse_status)

    call load_real_env('STAGE16_3_FEEDBACK_ALPHA', feedback_alpha_value)
    call load_real_env('STAGE16_3_MAX_ACTION_REACTION_ERROR', max_action_reaction_error)
    call load_real_env('STAGE16_3_MAX_SIGN_ERROR', max_sign_error)

    numeric_bounds_status = logical_to_int(feedback_alpha_value > 0.0_mytype .and. &
         max_action_reaction_error >= 0.0_mytype .and. max_sign_error >= 0.0_mytype)

  contains
    subroutine load_real_env(name, value)
      character(len=*), intent(in) :: name
      real(mytype), intent(inout) :: value

      call get_environment_variable(name, value=env_value, status=env_status)
      if (env_status == 0) then
        read(env_value, *, iostat=parse_status) parsed_real
        if (parse_status == 0 .and. finite_real(parsed_real)) then
          value = parsed_real
        else
          numeric_parse_status = 0
        end if
      end if
    end subroutine load_real_env

    subroutine load_bool_env(name, value, bool_status)
      character(len=*), intent(in) :: name
      logical, intent(out) :: value
      integer, intent(out) :: bool_status

      value = .true.
      bool_status = 1
      call get_environment_variable(name, value=env_value, status=env_status)
      if (env_status == 0) call parse_bool(env_value, value, bool_status)
    end subroutine load_bool_env
  end subroutine stage16_force_sign_audit_load_from_environment

  subroutine stage16_force_sign_audit_compute_reference(u_f, v_f, alpha)
    real(mytype), intent(in) :: u_f(3)
    real(mytype), intent(in) :: v_f(3)
    real(mytype), intent(in) :: alpha
    real(mytype) :: zero_slip(3)

    feedback_alpha_value = alpha
    slip(:) = u_f(:) - v_f(:)
    feedback_force(:) = feedback_alpha_value * slip(:)
    fluid_on_fibre(:) = feedback_force(:)
    fibre_on_fluid(:) = -feedback_force(:)
    reversed_fluid_on_fibre(:) = -feedback_alpha_value * slip(:)
    reversed_fibre_on_fluid(:) = feedback_alpha_value * slip(:)
    zero_slip(:) = 0.0_mytype
    zero_slip_force_value = maxval(abs(feedback_alpha_value * zero_slip))

    slip_finite_status = logical_to_int(all_finite_vector(slip))
    feedback_force_finite_status = logical_to_int(all_finite_vector(feedback_force))
    structure_side_force_sign_status = logical_to_int(vector_equal(fluid_on_fibre, feedback_force, max_sign_error))
    fluid_side_force_sign_status = logical_to_int(vector_equal(fibre_on_fluid, -feedback_force, max_sign_error))
    zero_slip_zero_force_status = logical_to_int(zero_slip_force_value <= max_sign_error)
    slip_reversal_sign_status = logical_to_int(vector_equal(reversed_fluid_on_fibre, -fluid_on_fibre, max_sign_error) .and. &
         vector_equal(reversed_fibre_on_fluid, -fibre_on_fluid, max_sign_error))
    call stage16_force_sign_audit_validate_action_reaction(max_action_reaction_error)
    call stage16_force_sign_audit_validate_wrong_sign_rejection(max_sign_error)
    call update_audit_status()
  end subroutine stage16_force_sign_audit_compute_reference

  subroutine stage16_force_sign_audit_validate_action_reaction(tolerance)
    real(mytype), intent(in) :: tolerance

    action_reaction_error_value = maxval(abs(fluid_on_fibre + fibre_on_fluid))
    action_reaction_status = logical_to_int(action_reaction_error_value <= tolerance)
  end subroutine stage16_force_sign_audit_validate_action_reaction

  subroutine stage16_force_sign_audit_validate_wrong_sign_rejection(tolerance)
    real(mytype), intent(in) :: tolerance
    real(mytype) :: wrong_fluid_side(3)

    wrong_fluid_side(:) = fluid_on_fibre(:)
    wrong_sign_error_value = maxval(abs(fluid_on_fibre + wrong_fluid_side))
    wrong_sign_rejection_status = logical_to_int(wrong_sign_error_value > tolerance)
  end subroutine stage16_force_sign_audit_validate_wrong_sign_rejection

  integer function stage16_force_sign_audit_get_final_status()
    stage16_force_sign_audit_get_final_status = audit_status
  end function stage16_force_sign_audit_get_final_status

  subroutine stage16_force_sign_audit_write_diagnostics(unit)
    integer, intent(in) :: unit

    write(unit,'(A,1X,I0)') 'stage16_3_requested_status', stage16_3_requested_status
    write(unit,'(A,1X,ES24.16E3)') 'feedback_alpha', feedback_alpha_value
    write(unit,'(A,1X,I0)') 'slip_finite_status', slip_finite_status
    write(unit,'(A,1X,I0)') 'feedback_force_finite_status', feedback_force_finite_status
    write(unit,'(A,1X,I0)') 'structure_side_force_sign_status', structure_side_force_sign_status
    write(unit,'(A,1X,I0)') 'fluid_side_force_sign_status', fluid_side_force_sign_status
    write(unit,'(A,1X,ES24.16E3)') 'action_reaction_error', action_reaction_error_value
    write(unit,'(A,1X,I0)') 'action_reaction_status', action_reaction_status
    write(unit,'(A,1X,I0)') 'wrong_sign_rejection_status', wrong_sign_rejection_status
    write(unit,'(A,1X,I0)') 'zero_slip_zero_force_status', zero_slip_zero_force_status
    write(unit,'(A,1X,I0)') 'slip_reversal_sign_status', slip_reversal_sign_status
    write(unit,'(A,1X,I0)') 'approved_stage12_13_14_chain_status', approved_stage12_13_14_chain_status
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
    write(unit,'(A,1X,I0)') 'numeric_parse_status', numeric_parse_status
    write(unit,'(A,1X,I0)') 'numeric_bounds_status', numeric_bounds_status
    write(unit,'(A,1X,I0)') 'diagnostic_only_status', diagnostic_only_status
    write(unit,'(A,1X,I0)') 'final_status', audit_status
  end subroutine stage16_force_sign_audit_write_diagnostics

  subroutine update_audit_status()
    audit_status = logical_to_int(stage16_3_requested_status == 1 .and. slip_finite_status == 1 .and. &
         feedback_force_finite_status == 1 .and. structure_side_force_sign_status == 1 .and. &
         fluid_side_force_sign_status == 1 .and. action_reaction_status == 1 .and. &
         wrong_sign_rejection_status == 1 .and. zero_slip_zero_force_status == 1 .and. &
         slip_reversal_sign_status == 1 .and. approved_stage12_13_14_chain_status == 1 .and. &
         numeric_parse_status == 1 .and. numeric_bounds_status == 1 .and. diagnostic_only_status == 1)
  end subroutine update_audit_status

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

  logical function all_finite_vector(values)
    real(mytype), intent(in) :: values(:)
    integer :: i

    all_finite_vector = .true.
    do i = 1, size(values)
      if (.not. finite_real(values(i))) all_finite_vector = .false.
    end do
  end function all_finite_vector

  logical function vector_equal(first, second, tolerance)
    real(mytype), intent(in) :: first(:)
    real(mytype), intent(in) :: second(:)
    real(mytype), intent(in) :: tolerance
    vector_equal = maxval(abs(first - second)) <= tolerance
  end function vector_equal

end module fibre_stage16_force_sign_audit
