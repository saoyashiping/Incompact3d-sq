module fibre_stage15_controlled_structure_step
  implicit none
  private

  integer, parameter :: mytype = kind(1.0d0)

  integer :: requested_status = 0
  integer :: controlled_step_enabled_status = 0
  integer :: structure_advance_enable_status = 0
  integer :: diagnostic_only_status = 0
  integer :: np_value = 0
  integer :: n_points = 0
  real(mytype) :: dt_value = 0.0_mytype
  real(mytype) :: rho_tilde_value = 0.0_mytype
  real(mytype) :: test_force_magnitude = 0.0_mytype
  real(mytype) :: zero_force_drift = huge(1.0_mytype)
  integer :: zero_force_status = 0
  integer :: small_force_status = 0
  real(mytype) :: max_acceleration = huge(1.0_mytype)
  real(mytype) :: max_velocity_update = huge(1.0_mytype)
  real(mytype) :: max_position_update = huge(1.0_mytype)
  integer :: forced_component_nonzero_status = 0
  integer :: sign_consistency_status = 0
  integer :: bounded_update_status = 0
  integer :: finite_value_status = 0
  integer :: controlled_update_count = 0
  integer :: production_full_structure_advance_count = 0
  integer :: bending_solve_count = 0
  integer :: tension_solve_count = 0
  integer :: wall_contact_count = 0
  integer :: multifibre_count = 0
  integer :: rhs_injection_connection_count = 0
  integer :: no_fluid_rhs_modification_status = 1
  integer :: no_pressure_projection_modification_status = 1
  integer :: no_poisson_modification_status = 1
  integer :: no_rk3_channel_forcing_modification_status = 1
  integer :: final_status = 0

  real(mytype), allocatable :: x_f(:,:)
  real(mytype), allocatable :: v_f(:,:)
  real(mytype), allocatable :: a_f(:,:)

  public :: stage15_controlled_structure_step_reset
  public :: stage15_controlled_structure_step_apply
  public :: stage15_controlled_structure_step_validate
  public :: stage15_controlled_structure_step_get_status_values
  public :: stage15_controlled_structure_step_write_diagnostics

contains

  subroutine stage15_controlled_structure_step_reset()
    if (allocated(x_f)) deallocate(x_f)
    if (allocated(v_f)) deallocate(v_f)
    if (allocated(a_f)) deallocate(a_f)
    requested_status = 0
    controlled_step_enabled_status = 0
    structure_advance_enable_status = 0
    diagnostic_only_status = 0
    np_value = 0
    n_points = 0
    dt_value = 0.0_mytype
    rho_tilde_value = 0.0_mytype
    test_force_magnitude = 0.0_mytype
    zero_force_drift = huge(1.0_mytype)
    zero_force_status = 0
    small_force_status = 0
    max_acceleration = huge(1.0_mytype)
    max_velocity_update = huge(1.0_mytype)
    max_position_update = huge(1.0_mytype)
    forced_component_nonzero_status = 0
    sign_consistency_status = 0
    bounded_update_status = 0
    finite_value_status = 0
    controlled_update_count = 0
    production_full_structure_advance_count = 0
    bending_solve_count = 0
    tension_solve_count = 0
    wall_contact_count = 0
    multifibre_count = 0
    rhs_injection_connection_count = 0
    no_fluid_rhs_modification_status = 1
    no_pressure_projection_modification_status = 1
    no_poisson_modification_status = 1
    no_rk3_channel_forcing_modification_status = 1
    final_status = 0
  end subroutine stage15_controlled_structure_step_reset

  subroutine stage15_controlled_structure_step_apply(npts, np, dt, rho_tilde, force_magnitude, max_pos_update, &
                                                     max_vel_update, max_accel, max_zero_drift, stage15_enable, &
                                                     controlled_enable, structure_enable, diagnostic_only)
    integer, intent(in) :: npts
    integer, intent(in) :: np
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: rho_tilde
    real(mytype), intent(in) :: force_magnitude
    real(mytype), intent(in) :: max_pos_update
    real(mytype), intent(in) :: max_vel_update
    real(mytype), intent(in) :: max_accel
    real(mytype), intent(in) :: max_zero_drift
    integer, intent(in) :: stage15_enable
    integer, intent(in) :: controlled_enable
    integer, intent(in) :: structure_enable
    integer, intent(in) :: diagnostic_only
    real(mytype), allocatable :: x_old(:,:)
    real(mytype), allocatable :: v_old(:,:)
    real(mytype), allocatable :: force_test(:,:)
    real(mytype), allocatable :: x_zero(:,:)
    real(mytype), allocatable :: v_zero(:,:)
    real(mytype), allocatable :: a_zero(:,:)
    integer :: alloc_status

    call stage15_controlled_structure_step_reset()
    requested_status = merge(1, 0, stage15_enable == 1)
    controlled_step_enabled_status = merge(1, 0, controlled_enable == 1)
    structure_advance_enable_status = merge(1, 0, structure_enable == 1)
    diagnostic_only_status = merge(1, 0, diagnostic_only == 1)
    np_value = np
    n_points = npts
    dt_value = dt
    rho_tilde_value = rho_tilde
    test_force_magnitude = force_magnitude

    if (.not. controls_are_valid(npts, np, dt, rho_tilde, force_magnitude, max_pos_update, &
                                 max_vel_update, max_accel, max_zero_drift)) then
      call stage15_controlled_structure_step_validate()
      return
    end if

    allocate(x_old(npts, 3), v_old(npts, 3), force_test(npts, 3), x_zero(npts, 3), &
             v_zero(npts, 3), a_zero(npts, 3), x_f(npts, 3), v_f(npts, 3), a_f(npts, 3), stat=alloc_status)
    if (alloc_status /= 0) then
      call stage15_controlled_structure_step_validate()
      return
    end if

    call initialize_reference_state(x_old, v_old)
    force_test(:, :) = 0.0_mytype
    call compute_stage15_3_candidate(x_old, v_old, force_test, dt, rho_tilde, a_zero, v_zero, x_zero)
    zero_force_drift = max(maxval(abs(a_zero)), max(maxval(abs(v_zero - v_old)), maxval(abs(x_zero - x_old))))
    zero_force_status = merge(1, 0, zero_force_drift <= max_zero_drift .and. all_candidate_finite(a_zero, v_zero, x_zero))

    force_test(:, :) = 0.0_mytype
    force_test(:, 1) = force_magnitude
    call compute_stage15_3_candidate(x_old, v_old, force_test, dt, rho_tilde, a_f, v_f, x_f)

    max_acceleration = maxval(abs(a_f))
    max_velocity_update = maxval(abs(v_f - v_old))
    max_position_update = maxval(abs(x_f - x_old))
    finite_value_status = merge(1, 0, all_candidate_finite(a_f, v_f, x_f))
    forced_component_nonzero_status = merge(1, 0, abs(force_magnitude) > 0.0_mytype .and. &
                                            maxval(abs(v_f(:, 1) - v_old(:, 1))) > 0.0_mytype .and. &
                                            maxval(abs(x_f(:, 1) - x_old(:, 1))) > 0.0_mytype)
    sign_consistency_status = sign_check(force_magnitude, a_f(:, 1), v_f(:, 1) - v_old(:, 1), x_f(:, 1) - x_old(:, 1))
    bounded_update_status = merge(1, 0, max_acceleration <= max_accel .and. &
                                  max_velocity_update <= max_vel_update .and. &
                                  max_position_update <= max_pos_update)
    small_force_status = merge(1, 0, finite_value_status == 1 .and. forced_component_nonzero_status == 1 .and. &
                               sign_consistency_status == 1 .and. bounded_update_status == 1)
    controlled_update_count = merge(1, 0, small_force_status == 1)
    call stage15_controlled_structure_step_validate()
  end subroutine stage15_controlled_structure_step_apply

  subroutine stage15_controlled_structure_step_validate()
    final_status = merge(1, 0, requested_status == 1 .and. controlled_step_enabled_status == 1 .and. &
                         structure_advance_enable_status == 1 .and. diagnostic_only_status == 1 .and. &
                         np_value == 1 .and. zero_force_status == 1 .and. small_force_status == 1 .and. &
                         finite_value_status == 1 .and. controlled_update_count == 1 .and. &
                         production_full_structure_advance_count == 0 .and. bending_solve_count == 0 .and. &
                         tension_solve_count == 0 .and. wall_contact_count == 0 .and. multifibre_count == 0 .and. &
                         rhs_injection_connection_count == 0 .and. no_fluid_rhs_modification_status == 1 .and. &
                         no_pressure_projection_modification_status == 1 .and. no_poisson_modification_status == 1 .and. &
                         no_rk3_channel_forcing_modification_status == 1)
  end subroutine stage15_controlled_structure_step_validate

  subroutine stage15_controlled_structure_step_get_status_values(requested_out, controlled_out, structure_enable_out, &
                                                                 diagnostic_only_out, np_out, npts_out, dt_out, rho_out, &
                                                                 force_out, zero_drift_out, zero_status_out, small_status_out, &
                                                                 max_accel_out, max_vel_out, max_pos_out, forced_nonzero_out, &
                                                                 sign_status_out, bounded_out, update_count_out, &
                                                                 full_advance_out, bending_out, tension_out, wall_out, &
                                                                 multifibre_out, rhs_out, no_fluid_out, no_pressure_out, &
                                                                 no_poisson_out, no_rk3_out, final_out)
    integer, intent(out) :: requested_out
    integer, intent(out) :: controlled_out
    integer, intent(out) :: structure_enable_out
    integer, intent(out) :: diagnostic_only_out
    integer, intent(out) :: np_out
    integer, intent(out) :: npts_out
    real(mytype), intent(out) :: dt_out
    real(mytype), intent(out) :: rho_out
    real(mytype), intent(out) :: force_out
    real(mytype), intent(out) :: zero_drift_out
    integer, intent(out) :: zero_status_out
    integer, intent(out) :: small_status_out
    real(mytype), intent(out) :: max_accel_out
    real(mytype), intent(out) :: max_vel_out
    real(mytype), intent(out) :: max_pos_out
    integer, intent(out) :: forced_nonzero_out
    integer, intent(out) :: sign_status_out
    integer, intent(out) :: bounded_out
    integer, intent(out) :: update_count_out
    integer, intent(out) :: full_advance_out
    integer, intent(out) :: bending_out
    integer, intent(out) :: tension_out
    integer, intent(out) :: wall_out
    integer, intent(out) :: multifibre_out
    integer, intent(out) :: rhs_out
    integer, intent(out) :: no_fluid_out
    integer, intent(out) :: no_pressure_out
    integer, intent(out) :: no_poisson_out
    integer, intent(out) :: no_rk3_out
    integer, intent(out) :: final_out

    requested_out = requested_status
    controlled_out = controlled_step_enabled_status
    structure_enable_out = structure_advance_enable_status
    diagnostic_only_out = diagnostic_only_status
    np_out = np_value
    npts_out = n_points
    dt_out = dt_value
    rho_out = rho_tilde_value
    force_out = test_force_magnitude
    zero_drift_out = zero_force_drift
    zero_status_out = zero_force_status
    small_status_out = small_force_status
    max_accel_out = max_acceleration
    max_vel_out = max_velocity_update
    max_pos_out = max_position_update
    forced_nonzero_out = forced_component_nonzero_status
    sign_status_out = sign_consistency_status
    bounded_out = bounded_update_status
    update_count_out = controlled_update_count
    full_advance_out = production_full_structure_advance_count
    bending_out = bending_solve_count
    tension_out = tension_solve_count
    wall_out = wall_contact_count
    multifibre_out = multifibre_count
    rhs_out = rhs_injection_connection_count
    no_fluid_out = no_fluid_rhs_modification_status
    no_pressure_out = no_pressure_projection_modification_status
    no_poisson_out = no_poisson_modification_status
    no_rk3_out = no_rk3_channel_forcing_modification_status
    final_out = final_status
  end subroutine stage15_controlled_structure_step_get_status_values

  subroutine stage15_controlled_structure_step_write_diagnostics(unit_id)
    integer, intent(in) :: unit_id

    write(unit_id, '(A,1X,I0)') 'stage15_6_requested_status', requested_status
    write(unit_id, '(A,1X,I0)') 'controlled_step_enabled_status', controlled_step_enabled_status
    write(unit_id, '(A,1X,I0)') 'structure_advance_enable_status', structure_advance_enable_status
    write(unit_id, '(A,1X,I0)') 'diagnostic_only_status', diagnostic_only_status
    write(unit_id, '(A,1X,I0)') 'np', np_value
    write(unit_id, '(A,1X,I0)') 'npts', n_points
    write(unit_id, '(A,1X,ES24.16)') 'dt', dt_value
    write(unit_id, '(A,1X,ES24.16)') 'rho_tilde', rho_tilde_value
    write(unit_id, '(A,1X,ES24.16)') 'test_force_magnitude', test_force_magnitude
    write(unit_id, '(A,1X,ES24.16)') 'zero_force_drift', zero_force_drift
    write(unit_id, '(A,1X,I0)') 'zero_force_status', zero_force_status
    write(unit_id, '(A,1X,I0)') 'small_force_status', small_force_status
    write(unit_id, '(A,1X,ES24.16)') 'max_acceleration', max_acceleration
    write(unit_id, '(A,1X,ES24.16)') 'max_velocity_update', max_velocity_update
    write(unit_id, '(A,1X,ES24.16)') 'max_position_update', max_position_update
    write(unit_id, '(A,1X,I0)') 'forced_component_nonzero_status', forced_component_nonzero_status
    write(unit_id, '(A,1X,I0)') 'sign_consistency_status', sign_consistency_status
    write(unit_id, '(A,1X,I0)') 'bounded_update_status', bounded_update_status
    write(unit_id, '(A,1X,I0)') 'controlled_update_count', controlled_update_count
    write(unit_id, '(A,1X,I0)') 'production_full_structure_advance_count', production_full_structure_advance_count
    write(unit_id, '(A,1X,I0)') 'bending_solve_count', bending_solve_count
    write(unit_id, '(A,1X,I0)') 'tension_solve_count', tension_solve_count
    write(unit_id, '(A,1X,I0)') 'wall_contact_count', wall_contact_count
    write(unit_id, '(A,1X,I0)') 'multifibre_count', multifibre_count
    write(unit_id, '(A,1X,I0)') 'rhs_injection_connection_count', rhs_injection_connection_count
    write(unit_id, '(A,1X,I0)') 'no_fluid_rhs_modification_status', no_fluid_rhs_modification_status
    write(unit_id, '(A,1X,I0)') 'no_pressure_projection_modification_status', no_pressure_projection_modification_status
    write(unit_id, '(A,1X,I0)') 'no_poisson_modification_status', no_poisson_modification_status
    write(unit_id, '(A,1X,I0)') 'no_rk3_channel_forcing_modification_status', no_rk3_channel_forcing_modification_status
    write(unit_id, '(A,1X,I0)') 'final_status', final_status
  end subroutine stage15_controlled_structure_step_write_diagnostics

  logical function controls_are_valid(npts, np, dt, rho_tilde, force_magnitude, max_pos_update, &
                                      max_vel_update, max_accel, max_zero_drift)
    integer, intent(in) :: npts
    integer, intent(in) :: np
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: rho_tilde
    real(mytype), intent(in) :: force_magnitude
    real(mytype), intent(in) :: max_pos_update
    real(mytype), intent(in) :: max_vel_update
    real(mytype), intent(in) :: max_accel
    real(mytype), intent(in) :: max_zero_drift

    controls_are_valid = requested_status == 1 .and. controlled_step_enabled_status == 1 .and. &
                         structure_advance_enable_status == 1 .and. diagnostic_only_status == 1 .and. &
                         np == 1 .and. npts > 0 .and. is_finite(dt) .and. is_finite(rho_tilde) .and. &
                         is_finite(force_magnitude) .and. is_finite(max_pos_update) .and. &
                         is_finite(max_vel_update) .and. is_finite(max_accel) .and. &
                         is_finite(max_zero_drift) .and. dt > 0.0_mytype .and. rho_tilde > 0.0_mytype .and. &
                         abs(force_magnitude) > 0.0_mytype .and. max_pos_update >= 0.0_mytype .and. &
                         max_vel_update >= 0.0_mytype .and. max_accel >= 0.0_mytype .and. max_zero_drift >= 0.0_mytype
  end function controls_are_valid

  subroutine initialize_reference_state(x_ref, v_ref)
    real(mytype), intent(out) :: x_ref(:,:)
    real(mytype), intent(out) :: v_ref(:,:)
    integer :: i
    integer :: local_npts

    local_npts = size(x_ref, 1)
    x_ref(:, :) = 0.0_mytype
    v_ref(:, :) = 0.0_mytype
    do i = 1, local_npts
      if (local_npts > 1) then
        x_ref(i, 1) = real(i - 1, mytype) / real(local_npts - 1, mytype)
      else
        x_ref(i, 1) = 0.0_mytype
      end if
    end do
  end subroutine initialize_reference_state

  subroutine compute_stage15_3_candidate(x_old, v_old, force_test, dt, rho_tilde, a_cand, v_new, x_new)
    real(mytype), intent(in) :: x_old(:,:)
    real(mytype), intent(in) :: v_old(:,:)
    real(mytype), intent(in) :: force_test(:,:)
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: rho_tilde
    real(mytype), intent(out) :: a_cand(:,:)
    real(mytype), intent(out) :: v_new(:,:)
    real(mytype), intent(out) :: x_new(:,:)

    ! Stage 15.6 intentionally reuses the Stage 15.3 checked explicit candidate formula:
    ! A_f_cand = F_test / rho_tilde; V_f_new = V_f_old + dt*A_f_cand; X_f_new = X_f_old + dt*V_f_new.
    a_cand(:, :) = force_test(:, :) / rho_tilde
    v_new(:, :) = v_old(:, :) + dt * a_cand(:, :)
    x_new(:, :) = x_old(:, :) + dt * v_new(:, :)
  end subroutine compute_stage15_3_candidate

  integer function sign_check(force_magnitude, a_component, v_delta_component, x_delta_component)
    real(mytype), intent(in) :: force_magnitude
    real(mytype), intent(in) :: a_component(:)
    real(mytype), intent(in) :: v_delta_component(:)
    real(mytype), intent(in) :: x_delta_component(:)

    if (force_magnitude > 0.0_mytype) then
      sign_check = merge(1, 0, minval(a_component) >= 0.0_mytype .and. &
                         minval(v_delta_component) >= 0.0_mytype .and. minval(x_delta_component) >= 0.0_mytype)
    else if (force_magnitude < 0.0_mytype) then
      sign_check = merge(1, 0, maxval(a_component) <= 0.0_mytype .and. &
                         maxval(v_delta_component) <= 0.0_mytype .and. maxval(x_delta_component) <= 0.0_mytype)
    else
      sign_check = 0
    end if
  end function sign_check

  logical function all_candidate_finite(a_cand, v_new, x_new)
    real(mytype), intent(in) :: a_cand(:,:)
    real(mytype), intent(in) :: v_new(:,:)
    real(mytype), intent(in) :: x_new(:,:)

    all_candidate_finite = all_finite_rank2(a_cand) .and. all_finite_rank2(v_new) .and. all_finite_rank2(x_new)
  end function all_candidate_finite

  logical function all_finite_rank2(values)
    real(mytype), intent(in) :: values(:,:)

    all_finite_rank2 = all(values == values) .and. all(abs(values) < huge(1.0_mytype))
  end function all_finite_rank2

  logical function is_finite(value)
    real(mytype), intent(in) :: value

    is_finite = (value == value) .and. (abs(value) < huge(1.0_mytype))
  end function is_finite

end module fibre_stage15_controlled_structure_step
