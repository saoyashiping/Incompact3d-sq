module fibre_stage14_rhs_sign_scaling_audit
  use decomp_2d_constants, only : mytype
  implicit none

  private

  public :: stage14_rhs_sign_scaling_audit_init
  public :: stage14_rhs_sign_scaling_audit_clear
  public :: stage14_rhs_sign_scaling_audit_finalize
  public :: stage14_rhs_sign_scaling_audit_compute_expected_forces
  public :: stage14_rhs_sign_scaling_audit_check_action_reaction
  public :: stage14_rhs_sign_scaling_audit_check_rhs_increment_sign
  public :: stage14_rhs_sign_scaling_audit_check_scaling
  public :: stage14_rhs_sign_scaling_audit_get_status_values
  public :: stage14_rhs_sign_scaling_audit_write_diagnostics
  public :: stage14_rhs_sign_scaling_audit_note_rhs_increment

  real(mytype), parameter :: action_reaction_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: sign_error_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: scaling_error_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: component_error_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: zero_lambda_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: wrong_sign_min_separation = 1.0e-8_mytype

  integer :: initialized_status = 0
  integer :: fluid_to_fibre_sign_status = 0
  integer :: fibre_to_fluid_sign_status = 0
  integer :: action_reaction_status = 0
  integer :: rhs_increment_uses_fibre_to_fluid_status = 0
  integer :: wrong_sign_rejection_status = 0
  integer :: lambda_zero_scaling_status = 0
  integer :: lambda_one_scaling_status = 0
  integer :: lambda_fractional_scaling_status = 0
  integer :: lambda_negative_scaling_status = 0
  integer :: component_scaling_status = 0
  integer :: integrated_rhs_sign_status = 0
  integer :: integrated_rhs_scaling_status = 0
  integer :: finite_rhs_increment_status = 0
  integer :: no_production_rhs_modification_status = 1
  integer :: no_xcompact3d_hook_status = 1
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
  integer :: rhs_sign_scaling_audit_status = 0

  real(mytype) :: action_reaction_max_error = 0.0_mytype
  real(mytype) :: lambda_zero_integrated_error_l2 = 0.0_mytype
  real(mytype) :: lambda_one_integrated_error_l2 = 0.0_mytype
  real(mytype) :: lambda_fractional_integrated_error_l2 = 0.0_mytype
  real(mytype) :: lambda_negative_integrated_error_l2 = 0.0_mytype
  real(mytype) :: wrong_sign_error_l2 = 0.0_mytype
  real(mytype) :: wrong_sign_separation = 0.0_mytype
  real(mytype) :: component_error_x = 0.0_mytype
  real(mytype) :: component_error_y = 0.0_mytype
  real(mytype) :: component_error_z = 0.0_mytype
  real(mytype) :: max_abs_rhs_increment = 0.0_mytype

contains

  subroutine stage14_rhs_sign_scaling_audit_init()
    initialized_status = 1
    fluid_to_fibre_sign_status = 0
    fibre_to_fluid_sign_status = 0
    action_reaction_status = 0
    rhs_increment_uses_fibre_to_fluid_status = 0
    wrong_sign_rejection_status = 0
    lambda_zero_scaling_status = 0
    lambda_one_scaling_status = 0
    lambda_fractional_scaling_status = 0
    lambda_negative_scaling_status = 0
    component_scaling_status = 0
    integrated_rhs_sign_status = 0
    integrated_rhs_scaling_status = 0
    finite_rhs_increment_status = 0
    no_production_rhs_modification_status = 1
    no_xcompact3d_hook_status = 1
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
    rhs_sign_scaling_audit_status = 0
    action_reaction_max_error = 0.0_mytype
    lambda_zero_integrated_error_l2 = 0.0_mytype
    lambda_one_integrated_error_l2 = 0.0_mytype
    lambda_fractional_integrated_error_l2 = 0.0_mytype
    lambda_negative_integrated_error_l2 = 0.0_mytype
    wrong_sign_error_l2 = 0.0_mytype
    wrong_sign_separation = 0.0_mytype
    component_error_x = 0.0_mytype
    component_error_y = 0.0_mytype
    component_error_z = 0.0_mytype
    max_abs_rhs_increment = 0.0_mytype
    call update_overall_status()
  end subroutine stage14_rhs_sign_scaling_audit_init

  subroutine stage14_rhs_sign_scaling_audit_clear()
    call stage14_rhs_sign_scaling_audit_init()
  end subroutine stage14_rhs_sign_scaling_audit_clear

  subroutine stage14_rhs_sign_scaling_audit_finalize()
    initialized_status = 0
    rhs_sign_scaling_audit_status = 0
  end subroutine stage14_rhs_sign_scaling_audit_finalize

  subroutine stage14_rhs_sign_scaling_audit_compute_expected_forces(f_fluid_to_fibre, &
                                                                    f_fibre_to_fluid, &
                                                                    ds_lag, &
                                                                    expected_fluid_to_fibre, &
                                                                    expected_fibre_to_fluid)
    real(mytype), intent(in) :: f_fluid_to_fibre(:,:)
    real(mytype), intent(in) :: f_fibre_to_fluid(:,:)
    real(mytype), intent(in) :: ds_lag(:)
    real(mytype), intent(out) :: expected_fluid_to_fibre(3)
    real(mytype), intent(out) :: expected_fibre_to_fluid(3)
    integer :: point_id
    integer :: comp_id

    if (initialized_status /= 1) call stage14_rhs_sign_scaling_audit_init()
    expected_fluid_to_fibre = 0.0_mytype
    expected_fibre_to_fluid = 0.0_mytype

    if (size(f_fluid_to_fibre,2) /= 3 .or. size(f_fibre_to_fluid,2) /= 3 .or. &
        size(f_fluid_to_fibre,1) /= size(f_fibre_to_fluid,1) .or. &
        size(ds_lag) /= size(f_fluid_to_fibre,1)) then
      call update_overall_status()
      return
    end if

    do point_id = 1, size(ds_lag)
      do comp_id = 1, 3
        expected_fluid_to_fibre(comp_id) = expected_fluid_to_fibre(comp_id) + &
                                           f_fluid_to_fibre(point_id, comp_id) * ds_lag(point_id)
        expected_fibre_to_fluid(comp_id) = expected_fibre_to_fluid(comp_id) + &
                                           f_fibre_to_fluid(point_id, comp_id) * ds_lag(point_id)
      end do
    end do

    fluid_to_fibre_sign_status = logical_to_int(vector_l2(expected_fluid_to_fibre) > action_reaction_abs_tol)
    fibre_to_fluid_sign_status = logical_to_int(vector_l2(expected_fibre_to_fluid + expected_fluid_to_fibre) <= &
                                                action_reaction_abs_tol)
    call update_overall_status()
  end subroutine stage14_rhs_sign_scaling_audit_compute_expected_forces

  subroutine stage14_rhs_sign_scaling_audit_check_action_reaction(f_fluid_to_fibre, &
                                                                  f_fibre_to_fluid, &
                                                                  max_action_reaction_error)
    real(mytype), intent(in) :: f_fluid_to_fibre(:,:)
    real(mytype), intent(in) :: f_fibre_to_fluid(:,:)
    real(mytype), intent(out) :: max_action_reaction_error

    if (size(f_fluid_to_fibre,1) /= size(f_fibre_to_fluid,1) .or. &
        size(f_fluid_to_fibre,2) /= size(f_fibre_to_fluid,2)) then
      max_action_reaction_error = huge(1.0_mytype)
      action_reaction_status = 0
      call update_overall_status()
      return
    end if

    max_action_reaction_error = maxval(abs(f_fluid_to_fibre + f_fibre_to_fluid))
    action_reaction_max_error = max_action_reaction_error
    action_reaction_status = logical_to_int(action_reaction_max_error <= action_reaction_abs_tol)
    call update_overall_status()
  end subroutine stage14_rhs_sign_scaling_audit_check_action_reaction

  subroutine stage14_rhs_sign_scaling_audit_check_rhs_increment_sign(integrated_rhs_increment, &
                                                                     expected_fluid_to_fibre, &
                                                                     expected_fibre_to_fluid, &
                                                                     lambda_14, &
                                                                     correct_sign_error, &
                                                                     wrong_sign_error)
    real(mytype), intent(in) :: integrated_rhs_increment(3)
    real(mytype), intent(in) :: expected_fluid_to_fibre(3)
    real(mytype), intent(in) :: expected_fibre_to_fluid(3)
    real(mytype), intent(in) :: lambda_14
    real(mytype), intent(out) :: correct_sign_error
    real(mytype), intent(out) :: wrong_sign_error

    correct_sign_error = vector_l2(integrated_rhs_increment - lambda_14 * expected_fibre_to_fluid)
    wrong_sign_error = vector_l2(integrated_rhs_increment - lambda_14 * expected_fluid_to_fibre)
    wrong_sign_error_l2 = wrong_sign_error
    wrong_sign_separation = wrong_sign_error
    rhs_increment_uses_fibre_to_fluid_status = logical_to_int(correct_sign_error <= sign_error_abs_tol)
    if (abs(lambda_14) > zero_lambda_abs_tol .and. vector_l2(expected_fluid_to_fibre) > sign_error_abs_tol) then
      wrong_sign_rejection_status = logical_to_int(wrong_sign_separation >= wrong_sign_min_separation)
    end if
    integrated_rhs_sign_status = rhs_increment_uses_fibre_to_fluid_status
    call update_overall_status()
  end subroutine stage14_rhs_sign_scaling_audit_check_rhs_increment_sign

  subroutine stage14_rhs_sign_scaling_audit_check_scaling(integrated_rhs_increment, &
                                                          expected_fibre_to_fluid, lambda_14)
    real(mytype), intent(in) :: integrated_rhs_increment(3)
    real(mytype), intent(in) :: expected_fibre_to_fluid(3)
    real(mytype), intent(in) :: lambda_14
    real(mytype) :: error_l2

    error_l2 = vector_l2(integrated_rhs_increment - lambda_14 * expected_fibre_to_fluid)
    component_error_x = abs(integrated_rhs_increment(1) - lambda_14 * expected_fibre_to_fluid(1))
    component_error_y = abs(integrated_rhs_increment(2) - lambda_14 * expected_fibre_to_fluid(2))
    component_error_z = abs(integrated_rhs_increment(3) - lambda_14 * expected_fibre_to_fluid(3))
    component_scaling_status = logical_to_int(component_error_x <= component_error_abs_tol .and. &
                                             component_error_y <= component_error_abs_tol .and. &
                                             component_error_z <= component_error_abs_tol)

    if (lambda_14 == 0.0_mytype) then
      lambda_zero_integrated_error_l2 = error_l2
      lambda_zero_scaling_status = logical_to_int(error_l2 <= zero_lambda_abs_tol)
    else if (lambda_14 == 1.0_mytype) then
      lambda_one_integrated_error_l2 = error_l2
      lambda_one_scaling_status = logical_to_int(error_l2 <= scaling_error_abs_tol)
    else if (abs(lambda_14 - 0.1_mytype) <= scaling_error_abs_tol) then
      lambda_fractional_integrated_error_l2 = error_l2
      lambda_fractional_scaling_status = logical_to_int(error_l2 <= scaling_error_abs_tol)
    else if (abs(lambda_14 + 0.1_mytype) <= scaling_error_abs_tol) then
      lambda_negative_integrated_error_l2 = error_l2
      lambda_negative_scaling_status = logical_to_int(error_l2 <= scaling_error_abs_tol)
    end if

    integrated_rhs_scaling_status = logical_to_int(lambda_zero_scaling_status == 1 .and. &
                                                   lambda_one_scaling_status == 1 .and. &
                                                   lambda_fractional_scaling_status == 1 .and. &
                                                   lambda_negative_scaling_status == 1)
    call update_overall_status()
  end subroutine stage14_rhs_sign_scaling_audit_check_scaling

  subroutine stage14_rhs_sign_scaling_audit_get_status_values(initialized_status_out, &
                                                              fluid_to_fibre_sign_status_out, &
                                                              fibre_to_fluid_sign_status_out, &
                                                              action_reaction_status_out, &
                                                              rhs_increment_uses_fibre_to_fluid_status_out, &
                                                              wrong_sign_rejection_status_out, &
                                                              lambda_zero_scaling_status_out, &
                                                              lambda_one_scaling_status_out, &
                                                              lambda_fractional_scaling_status_out, &
                                                              lambda_negative_scaling_status_out, &
                                                              component_scaling_status_out, &
                                                              integrated_rhs_sign_status_out, &
                                                              integrated_rhs_scaling_status_out, &
                                                              finite_rhs_increment_status_out, &
                                                              no_production_rhs_modification_status_out, &
                                                              no_xcompact3d_hook_status_out, &
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
                                                              rhs_sign_scaling_audit_status_out)
    integer, intent(out) :: initialized_status_out
    integer, intent(out) :: fluid_to_fibre_sign_status_out
    integer, intent(out) :: fibre_to_fluid_sign_status_out
    integer, intent(out) :: action_reaction_status_out
    integer, intent(out) :: rhs_increment_uses_fibre_to_fluid_status_out
    integer, intent(out) :: wrong_sign_rejection_status_out
    integer, intent(out) :: lambda_zero_scaling_status_out
    integer, intent(out) :: lambda_one_scaling_status_out
    integer, intent(out) :: lambda_fractional_scaling_status_out
    integer, intent(out) :: lambda_negative_scaling_status_out
    integer, intent(out) :: component_scaling_status_out
    integer, intent(out) :: integrated_rhs_sign_status_out
    integer, intent(out) :: integrated_rhs_scaling_status_out
    integer, intent(out) :: finite_rhs_increment_status_out
    integer, intent(out) :: no_production_rhs_modification_status_out
    integer, intent(out) :: no_xcompact3d_hook_status_out
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
    integer, intent(out) :: rhs_sign_scaling_audit_status_out

    initialized_status_out = initialized_status
    fluid_to_fibre_sign_status_out = fluid_to_fibre_sign_status
    fibre_to_fluid_sign_status_out = fibre_to_fluid_sign_status
    action_reaction_status_out = action_reaction_status
    rhs_increment_uses_fibre_to_fluid_status_out = rhs_increment_uses_fibre_to_fluid_status
    wrong_sign_rejection_status_out = wrong_sign_rejection_status
    lambda_zero_scaling_status_out = lambda_zero_scaling_status
    lambda_one_scaling_status_out = lambda_one_scaling_status
    lambda_fractional_scaling_status_out = lambda_fractional_scaling_status
    lambda_negative_scaling_status_out = lambda_negative_scaling_status
    component_scaling_status_out = component_scaling_status
    integrated_rhs_sign_status_out = integrated_rhs_sign_status
    integrated_rhs_scaling_status_out = integrated_rhs_scaling_status
    finite_rhs_increment_status_out = finite_rhs_increment_status
    no_production_rhs_modification_status_out = no_production_rhs_modification_status
    no_xcompact3d_hook_status_out = no_xcompact3d_hook_status
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
    rhs_sign_scaling_audit_status_out = rhs_sign_scaling_audit_status
  end subroutine stage14_rhs_sign_scaling_audit_get_status_values

  subroutine stage14_rhs_sign_scaling_audit_write_diagnostics(filename, requested_flag, &
                                                              rhs_injection_enabled_flag, &
                                                              injection_gain_finite_status)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: requested_flag
    integer, intent(in) :: rhs_injection_enabled_flag
    integer, intent(in) :: injection_gain_finite_status
    integer :: unit_id
    integer :: io_status

    open(newunit=unit_id, file=filename, status='replace', action='write', iostat=io_status)
    if (io_status /= 0) return

    write(unit_id, '(A,1X,I0)') 'stage14_3_requested_flag', requested_flag
    write(unit_id, '(A,1X,I0)') 'stage14_3_rhs_injection_enabled_flag', rhs_injection_enabled_flag
    write(unit_id, '(A,1X,I0)') 'stage14_3_injection_gain_finite_status', injection_gain_finite_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_initialized_status', initialized_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_fluid_to_fibre_sign_status', fluid_to_fibre_sign_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_fibre_to_fluid_sign_status', fibre_to_fluid_sign_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_action_reaction_status', action_reaction_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_rhs_increment_uses_fibre_to_fluid_status', &
                                 rhs_increment_uses_fibre_to_fluid_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_wrong_sign_rejection_status', wrong_sign_rejection_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_lambda_zero_scaling_status', lambda_zero_scaling_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_lambda_one_scaling_status', lambda_one_scaling_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_lambda_fractional_scaling_status', lambda_fractional_scaling_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_lambda_negative_scaling_status', lambda_negative_scaling_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_component_scaling_status', component_scaling_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_integrated_rhs_sign_status', integrated_rhs_sign_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_integrated_rhs_scaling_status', integrated_rhs_scaling_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_finite_rhs_increment_status', finite_rhs_increment_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_no_production_rhs_modification_status', &
                                 no_production_rhs_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_no_xcompact3d_hook_status', no_xcompact3d_hook_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_no_fluid_field_access_status', no_fluid_field_access_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_no_fluid_field_modification_status', no_fluid_field_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_no_pressure_modification_status', no_pressure_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_no_projection_modification_status', no_projection_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_no_rk3_modification_status', no_rk3_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_no_channel_forcing_modification_status', &
                                 no_channel_forcing_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_no_production_ibm_forcing_status', no_production_ibm_forcing_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_no_feedback_application_status', no_feedback_application_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_no_twoway_force_status', no_twoway_force_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_no_structure_advance_status', no_structure_advance_status
    write(unit_id, '(A,1X,I0)') 'stage14_3_rhs_sign_scaling_audit_status', rhs_sign_scaling_audit_status
    write(unit_id, '(A,1X,ES24.16)') 'stage14_3_action_reaction_max_error', action_reaction_max_error
    write(unit_id, '(A,1X,ES24.16)') 'stage14_3_lambda_zero_integrated_error_l2', &
                                      lambda_zero_integrated_error_l2
    write(unit_id, '(A,1X,ES24.16)') 'stage14_3_lambda_one_integrated_error_l2', &
                                      lambda_one_integrated_error_l2
    write(unit_id, '(A,1X,ES24.16)') 'stage14_3_lambda_fractional_integrated_error_l2', &
                                      lambda_fractional_integrated_error_l2
    write(unit_id, '(A,1X,ES24.16)') 'stage14_3_lambda_negative_integrated_error_l2', &
                                      lambda_negative_integrated_error_l2
    write(unit_id, '(A,1X,ES24.16)') 'stage14_3_wrong_sign_error_l2', wrong_sign_error_l2
    write(unit_id, '(A,1X,ES24.16)') 'stage14_3_wrong_sign_separation', wrong_sign_separation
    write(unit_id, '(A,1X,ES24.16)') 'stage14_3_component_error_x', component_error_x
    write(unit_id, '(A,1X,ES24.16)') 'stage14_3_component_error_y', component_error_y
    write(unit_id, '(A,1X,ES24.16)') 'stage14_3_component_error_z', component_error_z
    write(unit_id, '(A,1X,ES24.16)') 'stage14_3_max_abs_rhs_increment', max_abs_rhs_increment
    close(unit_id)
  end subroutine stage14_rhs_sign_scaling_audit_write_diagnostics

  subroutine stage14_rhs_sign_scaling_audit_note_rhs_increment(rhs_inc_x, rhs_inc_y, rhs_inc_z, &
                                                               integrated_rhs_increment)
    real(mytype), intent(in) :: rhs_inc_x(:,:,:)
    real(mytype), intent(in) :: rhs_inc_y(:,:,:)
    real(mytype), intent(in) :: rhs_inc_z(:,:,:)
    real(mytype), intent(in) :: integrated_rhs_increment(3)

    max_abs_rhs_increment = max(max_abs_rhs_increment, &
                                max(maxval(abs(rhs_inc_x)), max(maxval(abs(rhs_inc_y)), maxval(abs(rhs_inc_z)))))
    finite_rhs_increment_status = logical_to_int(all(finite_real(rhs_inc_x)) .and. &
                                                 all(finite_real(rhs_inc_y)) .and. &
                                                 all(finite_real(rhs_inc_z)) .and. &
                                                 all(finite_real(integrated_rhs_increment)))
    call update_overall_status()
  end subroutine stage14_rhs_sign_scaling_audit_note_rhs_increment

  subroutine update_overall_status()
    rhs_sign_scaling_audit_status = logical_to_int(initialized_status == 1 .and. &
                                                   fluid_to_fibre_sign_status == 1 .and. &
                                                   fibre_to_fluid_sign_status == 1 .and. &
                                                   action_reaction_status == 1 .and. &
                                                   rhs_increment_uses_fibre_to_fluid_status == 1 .and. &
                                                   wrong_sign_rejection_status == 1 .and. &
                                                   lambda_zero_scaling_status == 1 .and. &
                                                   lambda_one_scaling_status == 1 .and. &
                                                   lambda_fractional_scaling_status == 1 .and. &
                                                   lambda_negative_scaling_status == 1 .and. &
                                                   component_scaling_status == 1 .and. &
                                                   integrated_rhs_sign_status == 1 .and. &
                                                   integrated_rhs_scaling_status == 1 .and. &
                                                   finite_rhs_increment_status == 1 .and. &
                                                   no_production_rhs_modification_status == 1 .and. &
                                                   no_xcompact3d_hook_status == 1 .and. &
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
  end subroutine update_overall_status

  real(mytype) function vector_l2(vector)
    real(mytype), intent(in) :: vector(:)
    vector_l2 = sqrt(sum(vector * vector))
  end function vector_l2

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

end module fibre_stage14_rhs_sign_scaling_audit
