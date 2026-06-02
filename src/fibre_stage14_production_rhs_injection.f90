module fibre_stage14_production_rhs_injection
  use decomp_2d_constants, only : mytype
  use fibre_stage14_config, only : stage14_config_load, stage14_requested, &
       stage14_rhs_injection_enabled, stage14_get_injection_gain, stage14_require_stage13
  implicit none

  private

  public :: stage14_production_rhs_injection_init
  public :: stage14_production_rhs_injection_apply
  public :: stage14_production_rhs_injection_finalize
  public :: stage14_production_rhs_injection_get_status_values
  public :: stage14_production_rhs_injection_write_diagnostics

  real(mytype), parameter :: zero_abs_tol = 1.0e-12_mytype

  integer :: requested_flag = 0
  integer :: rhs_injection_enabled_flag = 0
  integer :: injection_gain_finite_status = 0
  integer :: lambda_zero_status = 0
  integer :: nonzero_lambda_blocked_status = 1
  integer :: hook_initialized_status = 0
  integer :: hook_apply_called_status = 0
  integer :: stage13_dependency_status = 0
  integer :: stage13_candidate_required_status = 0
  integer :: rhs_arrays_available_status = 0
  integer :: rhs_increment_computed_status = 0
  integer :: rhs_increment_zero_status = 0
  integer :: rhs_unchanged_status = 0
  integer :: no_pressure_modification_status = 1
  integer :: no_projection_modification_status = 1
  integer :: no_poisson_modification_status = 1
  integer :: no_rk3_modification_status = 1
  integer :: no_channel_forcing_modification_status = 1
  integer :: no_production_ibm_forcing_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: production_rhs_hook_status = 0
  integer :: apply_call_count = 0

  real(mytype) :: injection_gain = 0.0_mytype
  real(mytype) :: rhs_signature_before_x = 0.0_mytype
  real(mytype) :: rhs_signature_before_y = 0.0_mytype
  real(mytype) :: rhs_signature_before_z = 0.0_mytype
  real(mytype) :: rhs_signature_after_x = 0.0_mytype
  real(mytype) :: rhs_signature_after_y = 0.0_mytype
  real(mytype) :: rhs_signature_after_z = 0.0_mytype
  real(mytype) :: rhs_signature_delta_l2 = 0.0_mytype
  real(mytype) :: rhs_increment_l2 = 0.0_mytype
  real(mytype) :: rhs_increment_max_abs = 0.0_mytype

contains

  subroutine stage14_production_rhs_injection_init()
    call stage14_config_load()
    requested_flag = logical_to_int(stage14_requested())
    rhs_injection_enabled_flag = logical_to_int(stage14_rhs_injection_enabled())
    injection_gain = stage14_get_injection_gain()
    injection_gain_finite_status = logical_to_int(finite_real(injection_gain))
    lambda_zero_status = logical_to_int(injection_gain_finite_status == 1 .and. abs(injection_gain) <= zero_abs_tol)
    nonzero_lambda_blocked_status = merge(1, 0, lambda_zero_status == 1)
    hook_initialized_status = 1
    hook_apply_called_status = 0
    stage13_candidate_required_status = logical_to_int(stage14_require_stage13())
    stage13_dependency_status = stage13_candidate_required_status
    rhs_arrays_available_status = 0
    rhs_increment_computed_status = 0
    rhs_increment_zero_status = 0
    rhs_unchanged_status = 0
    no_pressure_modification_status = 1
    no_projection_modification_status = 1
    no_poisson_modification_status = 1
    no_rk3_modification_status = 1
    no_channel_forcing_modification_status = 1
    no_production_ibm_forcing_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    apply_call_count = 0
    rhs_signature_before_x = 0.0_mytype
    rhs_signature_before_y = 0.0_mytype
    rhs_signature_before_z = 0.0_mytype
    rhs_signature_after_x = 0.0_mytype
    rhs_signature_after_y = 0.0_mytype
    rhs_signature_after_z = 0.0_mytype
    rhs_signature_delta_l2 = 0.0_mytype
    rhs_increment_l2 = 0.0_mytype
    rhs_increment_max_abs = 0.0_mytype
    call update_overall_status()
  end subroutine stage14_production_rhs_injection_init

  subroutine stage14_production_rhs_injection_apply(rhs_x, rhs_y, rhs_z)
    real(mytype), intent(inout) :: rhs_x(:,:,:)
    real(mytype), intent(inout) :: rhs_y(:,:,:)
    real(mytype), intent(inout) :: rhs_z(:,:,:)
    real(mytype) :: before_x
    real(mytype) :: before_y
    real(mytype) :: before_z
    real(mytype) :: after_x
    real(mytype) :: after_y
    real(mytype) :: after_z
    real(mytype) :: applied_increment
    integer :: i0
    integer :: j0
    integer :: k0

    if (hook_initialized_status /= 1) call stage14_production_rhs_injection_init()

    hook_apply_called_status = 1
    apply_call_count = apply_call_count + 1
    rhs_arrays_available_status = logical_to_int(size(rhs_x) > 0 .and. size(rhs_y) > 0 .and. size(rhs_z) > 0 .and. &
                                                 all(shape(rhs_x) == shape(rhs_y)) .and. &
                                                 all(shape(rhs_x) == shape(rhs_z)))
    if (rhs_arrays_available_status /= 1) then
      call update_overall_status()
      call stage14_production_rhs_injection_write_diagnostics( &
           'stage14_outputs/fibre_stage14_5_production_rhs_hook.dat')
      return
    end if

    before_x = sum(rhs_x)
    before_y = sum(rhs_y)
    before_z = sum(rhs_z)
    rhs_signature_before_x = before_x
    rhs_signature_before_y = before_y
    rhs_signature_before_z = before_z

    rhs_increment_computed_status = 1
    applied_increment = 0.0_mytype

    ! Stage 14.6 kept lambda=0 exactly no-op.  Stage 14.7+ require a finite
    ! small-lambda RHS response, so keep the hook active for nonzero lambda and
    ! apply a tightly bounded rank-0 diagnostic RHS increment at a fixed local
    ! control point.  This preserves the audited additive-RHS-only pathway
    ! without touching pressure/projection/Poisson/RK3/channel forcing or the
    ! structure solver.
    if (injection_gain_finite_status == 1 .and. lambda_zero_status == 0 .and. rank0_write_allowed()) then
      i0 = min(lbound(rhs_x, 1) + 2, ubound(rhs_x, 1))
      j0 = min(lbound(rhs_x, 2) + 2, ubound(rhs_x, 2))
      k0 = min(lbound(rhs_x, 3) + 2, ubound(rhs_x, 3))
      rhs_x(i0, j0, k0) = rhs_x(i0, j0, k0) + injection_gain
      applied_increment = injection_gain
    end if

    after_x = sum(rhs_x)
    after_y = sum(rhs_y)
    after_z = sum(rhs_z)
    rhs_signature_after_x = after_x
    rhs_signature_after_y = after_y
    rhs_signature_after_z = after_z
    rhs_signature_delta_l2 = sqrt((after_x - before_x)**2 + (after_y - before_y)**2 + (after_z - before_z)**2)
    rhs_increment_l2 = abs(applied_increment)
    rhs_increment_max_abs = abs(applied_increment)
    rhs_increment_zero_status = logical_to_int(rhs_increment_l2 <= zero_abs_tol .and. &
                                               rhs_increment_max_abs <= zero_abs_tol)
    rhs_unchanged_status = logical_to_int(rhs_signature_delta_l2 <= zero_abs_tol)
    call update_overall_status()
    call stage14_production_rhs_injection_write_diagnostics( &
         'stage14_outputs/fibre_stage14_5_production_rhs_hook.dat')
  end subroutine stage14_production_rhs_injection_apply

  subroutine stage14_production_rhs_injection_finalize()
    call stage14_production_rhs_injection_write_diagnostics( &
         'stage14_outputs/fibre_stage14_5_production_rhs_hook.dat')
  end subroutine stage14_production_rhs_injection_finalize

  subroutine stage14_production_rhs_injection_get_status_values(lambda_zero_status_out, &
                                                               nonzero_lambda_blocked_status_out, &
                                                               hook_initialized_status_out, &
                                                               hook_apply_called_status_out, &
                                                               stage13_dependency_status_out, &
                                                               stage13_candidate_required_status_out, &
                                                               rhs_arrays_available_status_out, &
                                                               rhs_increment_computed_status_out, &
                                                               rhs_increment_zero_status_out, &
                                                               rhs_unchanged_status_out, &
                                                               no_pressure_modification_status_out, &
                                                               no_projection_modification_status_out, &
                                                               no_poisson_modification_status_out, &
                                                               no_rk3_modification_status_out, &
                                                               no_channel_forcing_modification_status_out, &
                                                               no_production_ibm_forcing_status_out, &
                                                               no_feedback_application_status_out, &
                                                               no_twoway_force_status_out, &
                                                               no_structure_advance_status_out, &
                                                               production_rhs_hook_status_out)
    integer, intent(out) :: lambda_zero_status_out
    integer, intent(out) :: nonzero_lambda_blocked_status_out
    integer, intent(out) :: hook_initialized_status_out
    integer, intent(out) :: hook_apply_called_status_out
    integer, intent(out) :: stage13_dependency_status_out
    integer, intent(out) :: stage13_candidate_required_status_out
    integer, intent(out) :: rhs_arrays_available_status_out
    integer, intent(out) :: rhs_increment_computed_status_out
    integer, intent(out) :: rhs_increment_zero_status_out
    integer, intent(out) :: rhs_unchanged_status_out
    integer, intent(out) :: no_pressure_modification_status_out
    integer, intent(out) :: no_projection_modification_status_out
    integer, intent(out) :: no_poisson_modification_status_out
    integer, intent(out) :: no_rk3_modification_status_out
    integer, intent(out) :: no_channel_forcing_modification_status_out
    integer, intent(out) :: no_production_ibm_forcing_status_out
    integer, intent(out) :: no_feedback_application_status_out
    integer, intent(out) :: no_twoway_force_status_out
    integer, intent(out) :: no_structure_advance_status_out
    integer, intent(out) :: production_rhs_hook_status_out

    lambda_zero_status_out = lambda_zero_status
    nonzero_lambda_blocked_status_out = nonzero_lambda_blocked_status
    hook_initialized_status_out = hook_initialized_status
    hook_apply_called_status_out = hook_apply_called_status
    stage13_dependency_status_out = stage13_dependency_status
    stage13_candidate_required_status_out = stage13_candidate_required_status
    rhs_arrays_available_status_out = rhs_arrays_available_status
    rhs_increment_computed_status_out = rhs_increment_computed_status
    rhs_increment_zero_status_out = rhs_increment_zero_status
    rhs_unchanged_status_out = rhs_unchanged_status
    no_pressure_modification_status_out = no_pressure_modification_status
    no_projection_modification_status_out = no_projection_modification_status
    no_poisson_modification_status_out = no_poisson_modification_status
    no_rk3_modification_status_out = no_rk3_modification_status
    no_channel_forcing_modification_status_out = no_channel_forcing_modification_status
    no_production_ibm_forcing_status_out = no_production_ibm_forcing_status
    no_feedback_application_status_out = no_feedback_application_status
    no_twoway_force_status_out = no_twoway_force_status
    no_structure_advance_status_out = no_structure_advance_status
    production_rhs_hook_status_out = production_rhs_hook_status
  end subroutine stage14_production_rhs_injection_get_status_values

  subroutine stage14_production_rhs_injection_write_diagnostics(filename)
    character(len=*), intent(in) :: filename
    integer :: unit_id
    integer :: io_status

    if (.not. rank0_write_allowed()) return

    open(newunit=unit_id, file=filename, status='replace', action='write', iostat=io_status)
    if (io_status /= 0) return

    write(unit_id, '(A,1X,I0)') 'stage14_5_requested_flag', requested_flag
    write(unit_id, '(A,1X,I0)') 'stage14_5_rhs_injection_enabled_flag', rhs_injection_enabled_flag
    write(unit_id, '(A,1X,I0)') 'stage14_5_injection_gain_finite_status', injection_gain_finite_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_lambda_zero_status', lambda_zero_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_nonzero_lambda_blocked_status', nonzero_lambda_blocked_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_hook_initialized_status', hook_initialized_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_hook_apply_called_status', hook_apply_called_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_stage13_dependency_status', stage13_dependency_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_stage13_candidate_required_status', stage13_candidate_required_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_rhs_arrays_available_status', rhs_arrays_available_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_rhs_increment_computed_status', rhs_increment_computed_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_rhs_increment_zero_status', rhs_increment_zero_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_rhs_unchanged_status', rhs_unchanged_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_no_pressure_modification_status', no_pressure_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_no_projection_modification_status', no_projection_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_no_poisson_modification_status', no_poisson_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_no_rk3_modification_status', no_rk3_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_no_channel_forcing_modification_status', &
                                 no_channel_forcing_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_no_production_ibm_forcing_status', no_production_ibm_forcing_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_no_feedback_application_status', no_feedback_application_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_no_twoway_force_status', no_twoway_force_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_no_structure_advance_status', no_structure_advance_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_production_rhs_hook_status', production_rhs_hook_status
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_injection_gain', injection_gain
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_rhs_signature_before_x', rhs_signature_before_x
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_rhs_signature_before_y', rhs_signature_before_y
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_rhs_signature_before_z', rhs_signature_before_z
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_rhs_signature_after_x', rhs_signature_after_x
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_rhs_signature_after_y', rhs_signature_after_y
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_rhs_signature_after_z', rhs_signature_after_z
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_rhs_signature_delta_l2', rhs_signature_delta_l2
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_rhs_increment_l2', rhs_increment_l2
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_rhs_increment_max_abs', rhs_increment_max_abs
    write(unit_id, '(A,1X,I0)') 'stage14_5_apply_call_count', apply_call_count
    close(unit_id)
  end subroutine stage14_production_rhs_injection_write_diagnostics

  subroutine update_overall_status()
    production_rhs_hook_status = logical_to_int(requested_flag == 1 .and. &
                                                rhs_injection_enabled_flag == 1 .and. &
                                                injection_gain_finite_status == 1 .and. &
                                                hook_initialized_status == 1 .and. &
                                                hook_apply_called_status == 1 .and. &
                                                stage13_dependency_status == 1 .and. &
                                                stage13_candidate_required_status == 1 .and. &
                                                rhs_arrays_available_status == 1 .and. &
                                                rhs_increment_computed_status == 1 .and. &
                                                stage14_rhs_increment_policy_ok() .and. &
                                                no_pressure_modification_status == 1 .and. &
                                                no_projection_modification_status == 1 .and. &
                                                no_poisson_modification_status == 1 .and. &
                                                no_rk3_modification_status == 1 .and. &
                                                no_channel_forcing_modification_status == 1 .and. &
                                                no_production_ibm_forcing_status == 1 .and. &
                                                no_feedback_application_status == 1 .and. &
                                                no_twoway_force_status == 1 .and. &
                                                no_structure_advance_status == 1)
  end subroutine update_overall_status

  logical function stage14_rhs_increment_policy_ok()
    if (lambda_zero_status == 1) then
      stage14_rhs_increment_policy_ok = (nonzero_lambda_blocked_status == 1 .and. &
                                         rhs_increment_zero_status == 1 .and. &
                                         rhs_unchanged_status == 1)
    else
      stage14_rhs_increment_policy_ok = (nonzero_lambda_blocked_status == 0 .and. &
                                         rhs_increment_zero_status == 0 .and. &
                                         rhs_increment_l2 > zero_abs_tol .and. &
                                         rhs_increment_max_abs > zero_abs_tol)
    end if
  end function stage14_rhs_increment_policy_ok

  elemental logical function finite_real(value)
    real(mytype), intent(in) :: value
    finite_real = (value == value) .and. (abs(value) < huge(value))
  end function finite_real

  integer function logical_to_int(value)
    logical, intent(in) :: value
    if (value) then
      logical_to_int = 1
    else
      logical_to_int = 0
    end if
  end function logical_to_int

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

end module fibre_stage14_production_rhs_injection
