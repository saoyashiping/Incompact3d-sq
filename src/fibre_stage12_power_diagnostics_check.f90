program fibre_stage12_power_diagnostics_check
  use fibre_stage12_config, only: stage12_config_load, stage12_requested, stage12_readonly_mode, &
                                  stage12_get_max_points
  use fibre_stage12_feedback_formula, only: stage12_feedback_formula_init, &
                                            stage12_feedback_formula_compute_controlled, &
                                            stage12_feedback_formula_finalize
  use fibre_stage12_sign_convention_audit, only: stage12_sign_convention_audit_init, &
                                                 stage12_sign_convention_compute_pair, &
                                                 stage12_sign_convention_audit_finalize
  use fibre_stage12_power_diagnostics, only: stage12_power_diagnostics_init, &
                                             stage12_power_diagnostics_compute, &
                                             stage12_power_diagnostics_finalize, &
                                             stage12_power_diagnostics_get_status_values, &
                                             stage12_power_diagnostics_record_statuses
  implicit none

  integer, parameter :: mytype = kind(1.0d0)
  real(mytype), parameter :: zero_power_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: power_formula_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: pair_consistency_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: gain_scaling_abs_tol = 1.0e-12_mytype

  integer :: n_points
  integer :: alloc_status
  integer :: ios
  integer :: io_unit
  integer :: requested_flag
  integer :: readonly_mode_status
  integer :: initialized_status
  integer :: zero_slip_power_status
  integer :: positive_slip_power_status
  integer :: force_norm_finite_status
  integer :: slip_power_finite_status
  integer :: structure_power_finite_status
  integer :: fluid_power_finite_status
  integer :: pair_power_consistency_status
  integer :: gain_scaling_power_status
  integer :: action_reaction_power_status
  integer :: finite_diagnostics_status
  integer :: no_eulerian_force_density_status
  integer :: no_rhs_injection_status
  integer :: no_ibm_spreading_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: power_diagnostics_status
  real(mytype) :: zero_slip_power_max_error
  real(mytype) :: total_force_l2
  real(mytype) :: p_slip
  real(mytype) :: p_structure
  real(mytype) :: p_fluid
  real(mytype) :: p_pair
  real(mytype) :: pair_power_consistency_error
  real(mytype) :: gain_scaling_power_error
  real(mytype) :: force_l2_gain_scaling_error
  real(mytype) :: positive_slip_power_error
  real(mytype) :: action_reaction_power_error
  real(mytype), allocatable :: u_f(:,:)
  real(mytype), allocatable :: v_f(:,:)
  real(mytype), allocatable :: f_fluid_to_fibre(:,:)
  real(mytype), allocatable :: f_fibre_to_fluid(:,:)
  real(mytype), allocatable :: force_norm(:)
  real(mytype), allocatable :: formula_force(:,:)
  real(mytype), allocatable :: formula_norm(:)

  call execute_command_line('mkdir -p stage12_outputs', exitstat=ios)
  call stage12_config_load()
  requested_flag = logical_status(stage12_requested())
  readonly_mode_status = logical_status(stage12_readonly_mode())
  n_points = stage12_get_max_points()
  if (n_points <= 0) n_points = 8

  allocate(u_f(n_points, 3), v_f(n_points, 3), f_fluid_to_fibre(n_points, 3), &
           f_fibre_to_fluid(n_points, 3), force_norm(n_points), formula_force(n_points, 3), &
           formula_norm(n_points), stat=alloc_status)
  if (alloc_status /= 0) then
    call write_failed_allocation(requested_flag, readonly_mode_status)
    print *, 'STAGE 12.5 POWER DIAGNOSTICS VERDICT: FAIL'
    print *, 'Reason: allocation failed for power diagnostics arrays.'
    stop 1
  end if

  call stage12_feedback_formula_init()
  call stage12_sign_convention_audit_init()
  call stage12_power_diagnostics_init()

  zero_slip_power_max_error = run_zero_slip_test(n_points, u_f, v_f, f_fluid_to_fibre, f_fibre_to_fluid, &
                                                 force_norm, total_force_l2, p_slip, p_structure, p_fluid, p_pair)
  call run_positive_slip_test(n_points, u_f, v_f, f_fluid_to_fibre, f_fibre_to_fluid, force_norm, &
                              formula_force, formula_norm, total_force_l2, p_slip, p_structure, p_fluid, &
                              p_pair, positive_slip_power_error, pair_power_consistency_error, &
                              action_reaction_power_error)
  call run_mixed_slip_test(n_points, u_f, v_f, f_fluid_to_fibre, f_fibre_to_fluid, force_norm, &
                           total_force_l2, p_slip, p_structure, p_fluid, p_pair, &
                           pair_power_consistency_error, action_reaction_power_error)
  call run_gain_scaling_test(n_points, u_f, v_f, f_fluid_to_fibre, f_fibre_to_fluid, force_norm, &
                             gain_scaling_power_error, force_l2_gain_scaling_error)

  zero_slip_power_status = status_from_error(zero_slip_power_max_error, zero_power_abs_tol)
  positive_slip_power_status = status_from_error(positive_slip_power_error, power_formula_abs_tol)
  pair_power_consistency_status = status_from_error(pair_power_consistency_error, pair_consistency_abs_tol)
  gain_scaling_power_status = status_from_error(gain_scaling_power_error, gain_scaling_abs_tol)
  if (force_l2_gain_scaling_error > gain_scaling_abs_tol) gain_scaling_power_status = 0
  action_reaction_power_status = status_from_error(action_reaction_power_error, pair_consistency_abs_tol)

  call stage12_power_diagnostics_record_statuses(zero_slip_power_status, positive_slip_power_status, &
                                                 pair_power_consistency_status, gain_scaling_power_status, &
                                                 action_reaction_power_status, 1)

  call stage12_power_diagnostics_get_status_values(initialized_status, zero_slip_power_status, &
                                                   positive_slip_power_status, force_norm_finite_status, &
                                                   slip_power_finite_status, structure_power_finite_status, &
                                                   fluid_power_finite_status, pair_power_consistency_status, &
                                                   gain_scaling_power_status, action_reaction_power_status, &
                                                   finite_diagnostics_status, no_eulerian_force_density_status, &
                                                   no_rhs_injection_status, no_ibm_spreading_status, &
                                                   no_feedback_application_status, no_twoway_force_status, &
                                                   no_structure_advance_status, no_fluid_field_access_status, &
                                                   no_fluid_field_modification_status, power_diagnostics_status)

  open(newunit=io_unit, file='stage12_outputs/fibre_stage12_5_power_diagnostics.dat', status='replace', &
       action='write', iostat=ios)
  if (ios /= 0) then
    print *, 'STAGE 12.5 POWER DIAGNOSTICS VERDICT: FAIL'
    print *, 'Reason: could not open stage12_outputs/fibre_stage12_5_power_diagnostics.dat.'
    stop 1
  end if

  call write_int(io_unit, 'stage12_5_requested_flag', requested_flag)
  call write_int(io_unit, 'stage12_5_readonly_mode_status', readonly_mode_status)
  call write_int(io_unit, 'stage12_5_initialized_status', initialized_status)
  call write_int(io_unit, 'stage12_5_zero_slip_power_status', zero_slip_power_status)
  call write_int(io_unit, 'stage12_5_positive_slip_power_status', positive_slip_power_status)
  call write_int(io_unit, 'stage12_5_force_norm_finite_status', force_norm_finite_status)
  call write_int(io_unit, 'stage12_5_slip_power_finite_status', slip_power_finite_status)
  call write_int(io_unit, 'stage12_5_structure_power_finite_status', structure_power_finite_status)
  call write_int(io_unit, 'stage12_5_fluid_power_finite_status', fluid_power_finite_status)
  call write_int(io_unit, 'stage12_5_pair_power_consistency_status', pair_power_consistency_status)
  call write_int(io_unit, 'stage12_5_gain_scaling_power_status', gain_scaling_power_status)
  call write_int(io_unit, 'stage12_5_action_reaction_power_status', action_reaction_power_status)
  call write_int(io_unit, 'stage12_5_finite_diagnostics_status', finite_diagnostics_status)
  call write_int(io_unit, 'stage12_5_no_eulerian_force_density_status', no_eulerian_force_density_status)
  call write_int(io_unit, 'stage12_5_no_rhs_injection_status', no_rhs_injection_status)
  call write_int(io_unit, 'stage12_5_no_ibm_spreading_status', no_ibm_spreading_status)
  call write_int(io_unit, 'stage12_5_no_feedback_application_status', no_feedback_application_status)
  call write_int(io_unit, 'stage12_5_no_twoway_force_status', no_twoway_force_status)
  call write_int(io_unit, 'stage12_5_no_structure_advance_status', no_structure_advance_status)
  call write_int(io_unit, 'stage12_5_no_fluid_field_access_status', no_fluid_field_access_status)
  call write_int(io_unit, 'stage12_5_no_fluid_field_modification_status', no_fluid_field_modification_status)
  call write_int(io_unit, 'stage12_5_power_diagnostics_status', power_diagnostics_status)
  call write_real(io_unit, 'stage12_5_zero_slip_power_max_error', zero_slip_power_max_error)
  call write_real(io_unit, 'stage12_5_total_force_l2', total_force_l2)
  call write_real(io_unit, 'stage12_5_p_slip', p_slip)
  call write_real(io_unit, 'stage12_5_p_structure', p_structure)
  call write_real(io_unit, 'stage12_5_p_fluid', p_fluid)
  call write_real(io_unit, 'stage12_5_p_pair', p_pair)
  call write_real(io_unit, 'stage12_5_pair_power_consistency_error', pair_power_consistency_error)
  call write_real(io_unit, 'stage12_5_gain_scaling_power_error', gain_scaling_power_error)
  call write_real(io_unit, 'stage12_5_force_l2_gain_scaling_error', force_l2_gain_scaling_error)
  close(io_unit)

  if (power_diagnostics_status == 1 .and. requested_flag == 1 .and. readonly_mode_status == 1) then
    print *, 'STAGE 12.5 POWER DIAGNOSTICS VERDICT: PASS'
  else
    print *, 'STAGE 12.5 POWER DIAGNOSTICS VERDICT: FAIL'
    if (requested_flag /= 1) print *, 'Reason: Stage 12 feedback candidate was not requested.'
    if (readonly_mode_status /= 1) print *, 'Reason: Stage 12 readonly mode was not enforced.'
    if (zero_slip_power_status /= 1) print *, 'Reason: zero-slip power diagnostic failed.'
    if (positive_slip_power_status /= 1) print *, 'Reason: positive-slip power diagnostic failed.'
    if (pair_power_consistency_status /= 1) print *, 'Reason: pair-power consistency failed.'
    if (gain_scaling_power_status /= 1) print *, 'Reason: gain-scaling power diagnostic failed.'
    if (action_reaction_power_status /= 1) print *, 'Reason: action-reaction power diagnostic failed.'
    if (finite_diagnostics_status /= 1) print *, 'Reason: non-finite power diagnostic value detected.'
    stop 1
  end if

  call stage12_power_diagnostics_finalize()
  call stage12_sign_convention_audit_finalize()
  call stage12_feedback_formula_finalize()

contains

  integer function logical_status(value)
    logical, intent(in) :: value
    if (value) then
      logical_status = 1
    else
      logical_status = 0
    end if
  end function logical_status

  integer function status_from_error(error_value, tolerance)
    real(mytype), intent(in) :: error_value
    real(mytype), intent(in) :: tolerance
    if (error_value <= tolerance) then
      status_from_error = 1
    else
      status_from_error = 0
    end if
  end function status_from_error

  subroutine write_int(io_unit_in, key, value)
    integer, intent(in) :: io_unit_in
    character(len=*), intent(in) :: key
    integer, intent(in) :: value
    write(io_unit_in, '(A,1X,I0)') trim(key), value
  end subroutine write_int

  subroutine write_real(io_unit_in, key, value)
    integer, intent(in) :: io_unit_in
    character(len=*), intent(in) :: key
    real(mytype), intent(in) :: value
    write(io_unit_in, '(A,1X,ES24.16)') trim(key), value
  end subroutine write_real

  subroutine write_failed_allocation(requested_status, readonly_status)
    integer, intent(in) :: requested_status
    integer, intent(in) :: readonly_status
    integer :: failed_unit
    integer :: failed_ios

    open(newunit=failed_unit, file='stage12_outputs/fibre_stage12_5_power_diagnostics.dat', &
         status='replace', action='write', iostat=failed_ios)
    if (failed_ios /= 0) return
    call write_int(failed_unit, 'stage12_5_requested_flag', requested_status)
    call write_int(failed_unit, 'stage12_5_readonly_mode_status', readonly_status)
    call write_int(failed_unit, 'stage12_5_initialized_status', 0)
    call write_int(failed_unit, 'stage12_5_power_diagnostics_status', 0)
    close(failed_unit)
  end subroutine write_failed_allocation

  real(mytype) function run_zero_slip_test(n_local, u_local, v_local, fluid_force, fibre_force, norm_local, &
                                           l2_value, slip_power, structure_power, fluid_power, pair_power)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: fluid_force(:,:)
    real(mytype), intent(inout) :: fibre_force(:,:)
    real(mytype), intent(inout) :: norm_local(:)
    real(mytype), intent(out) :: l2_value
    real(mytype), intent(out) :: slip_power
    real(mytype), intent(out) :: structure_power
    real(mytype), intent(out) :: fluid_power
    real(mytype), intent(out) :: pair_power
    integer :: k

    do k = 1, n_local
      u_local(k, :) = (/0.05_mytype * real(k, mytype), -0.10_mytype, 0.20_mytype/)
    end do
    v_local(:, :) = u_local(:, :)
    fluid_force(:, :) = 0.0_mytype
    fibre_force(:, :) = 0.0_mytype
    call stage12_power_diagnostics_compute(u_local, v_local, fluid_force, fibre_force, norm_local, l2_value, &
                                           slip_power, structure_power, fluid_power, pair_power)
    run_zero_slip_test = max(maxval(abs(norm_local(:))), abs(l2_value))
    run_zero_slip_test = max(run_zero_slip_test, abs(slip_power))
    run_zero_slip_test = max(run_zero_slip_test, abs(structure_power))
    run_zero_slip_test = max(run_zero_slip_test, abs(fluid_power))
    run_zero_slip_test = max(run_zero_slip_test, abs(pair_power))
  end function run_zero_slip_test

  subroutine run_positive_slip_test(n_local, u_local, v_local, fluid_force, fibre_force, norm_local, &
                                    formula_force, formula_norm, l2_value, slip_power, structure_power, &
                                    fluid_power, pair_power, positive_error, pair_error, action_error)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: fluid_force(:,:)
    real(mytype), intent(inout) :: fibre_force(:,:)
    real(mytype), intent(inout) :: norm_local(:)
    real(mytype), intent(inout) :: formula_force(:,:)
    real(mytype), intent(inout) :: formula_norm(:)
    real(mytype), intent(out) :: l2_value
    real(mytype), intent(out) :: slip_power
    real(mytype), intent(out) :: structure_power
    real(mytype), intent(out) :: fluid_power
    real(mytype), intent(out) :: pair_power
    real(mytype), intent(out) :: positive_error
    real(mytype), intent(out) :: pair_error
    real(mytype), intent(out) :: action_error
    real(mytype) :: slip(3)
    real(mytype) :: expected_slip_power

    slip = (/0.10_mytype, 0.20_mytype, 0.30_mytype/)
    call fill_constant_slip(n_local, slip, u_local, v_local)
    call stage12_feedback_formula_compute_controlled(u_local, v_local, 2.0_mytype, 1, formula_force, formula_norm)
    call stage12_sign_convention_compute_pair(u_local, v_local, 2.0_mytype, fluid_force, fibre_force, norm_local)
    call stage12_power_diagnostics_compute(u_local, v_local, fluid_force, fibre_force, norm_local, l2_value, &
                                           slip_power, structure_power, fluid_power, pair_power)
    expected_slip_power = 2.0_mytype * real(n_local, mytype) * sum(slip(:) * slip(:))
    positive_error = abs(slip_power - expected_slip_power)
    positive_error = max(positive_error, maxval(abs(formula_force(:, :) - fluid_force(:, :))))
    positive_error = max(positive_error, maxval(abs(formula_norm(:) - norm_local(:))))
    pair_error = abs(pair_power + slip_power)
    action_error = maxval(abs(fluid_force(:, :) + fibre_force(:, :)))
  end subroutine run_positive_slip_test

  subroutine run_mixed_slip_test(n_local, u_local, v_local, fluid_force, fibre_force, norm_local, &
                                 l2_value, slip_power, structure_power, fluid_power, pair_power, &
                                 pair_error, action_error)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: fluid_force(:,:)
    real(mytype), intent(inout) :: fibre_force(:,:)
    real(mytype), intent(inout) :: norm_local(:)
    real(mytype), intent(out) :: l2_value
    real(mytype), intent(out) :: slip_power
    real(mytype), intent(out) :: structure_power
    real(mytype), intent(out) :: fluid_power
    real(mytype), intent(out) :: pair_power
    real(mytype), intent(inout) :: pair_error
    real(mytype), intent(inout) :: action_error
    integer :: k

    do k = 1, n_local
      v_local(k, :) = (/0.15_mytype, -0.05_mytype, 0.25_mytype/)
      u_local(k, 1) = v_local(k, 1) + 0.02_mytype * real(k, mytype)
      u_local(k, 2) = v_local(k, 2) - 0.03_mytype * real(k, mytype)
      u_local(k, 3) = v_local(k, 3) + 0.01_mytype * real(k, mytype)
    end do
    call stage12_sign_convention_compute_pair(u_local, v_local, 1.5_mytype, fluid_force, fibre_force, norm_local)
    call stage12_power_diagnostics_compute(u_local, v_local, fluid_force, fibre_force, norm_local, l2_value, &
                                           slip_power, structure_power, fluid_power, pair_power)
    pair_error = max(pair_error, abs(pair_power + slip_power))
    action_error = max(action_error, maxval(abs(fluid_force(:, :) + fibre_force(:, :))))
  end subroutine run_mixed_slip_test

  subroutine run_gain_scaling_test(n_local, u_local, v_local, fluid_force, fibre_force, norm_local, &
                                   gain_error, l2_gain_error)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: fluid_force(:,:)
    real(mytype), intent(inout) :: fibre_force(:,:)
    real(mytype), intent(inout) :: norm_local(:)
    real(mytype), intent(out) :: gain_error
    real(mytype), intent(out) :: l2_gain_error
    real(mytype) :: p_slip_one
    real(mytype) :: p_slip_two
    real(mytype) :: p_structure_local
    real(mytype) :: p_fluid_local
    real(mytype) :: p_pair_local
    real(mytype) :: l2_one
    real(mytype) :: l2_two
    real(mytype) :: slip(3)

    slip = (/0.07_mytype, -0.04_mytype, 0.02_mytype/)
    call fill_constant_slip(n_local, slip, u_local, v_local)
    call stage12_sign_convention_compute_pair(u_local, v_local, 1.0_mytype, fluid_force, fibre_force, norm_local)
    call stage12_power_diagnostics_compute(u_local, v_local, fluid_force, fibre_force, norm_local, l2_one, &
                                           p_slip_one, p_structure_local, p_fluid_local, p_pair_local)
    call stage12_sign_convention_compute_pair(u_local, v_local, 2.0_mytype, fluid_force, fibre_force, norm_local)
    call stage12_power_diagnostics_compute(u_local, v_local, fluid_force, fibre_force, norm_local, l2_two, &
                                           p_slip_two, p_structure_local, p_fluid_local, p_pair_local)
    gain_error = abs(p_slip_two - 2.0_mytype * p_slip_one)
    l2_gain_error = abs(l2_two - 2.0_mytype * l2_one)
  end subroutine run_gain_scaling_test

  subroutine fill_constant_slip(n_local, slip, u_local, v_local)
    integer, intent(in) :: n_local
    real(mytype), intent(in) :: slip(3)
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    integer :: k

    do k = 1, n_local
      v_local(k, :) = (/0.20_mytype, -0.10_mytype, 0.05_mytype/)
      u_local(k, :) = v_local(k, :) + slip(:)
    end do
  end subroutine fill_constant_slip

end program fibre_stage12_power_diagnostics_check
