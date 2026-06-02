program fibre_stage12_sign_convention_audit_check
  use fibre_stage12_config, only: stage12_config_load, stage12_requested, stage12_readonly_mode, &
                                  stage12_get_max_points
  use fibre_stage12_sign_convention_audit, only: stage12_sign_convention_audit_init, &
                                                 stage12_sign_convention_compute_pair, &
                                                 stage12_sign_convention_audit_finalize, &
                                                 stage12_sign_convention_audit_get_status_values, &
                                                 stage12_sign_convention_audit_record_statuses
  implicit none

  integer, parameter :: mytype = kind(1.0d0)
  real(mytype), parameter :: zero_slip_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: sign_formula_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: action_reaction_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: force_norm_abs_tol = 1.0e-12_mytype

  integer :: n_points
  integer :: alloc_status
  integer :: ios
  integer :: io_unit
  integer :: requested_flag
  integer :: readonly_mode_status
  integer :: initialized_status
  integer :: fluid_to_fibre_sign_status
  integer :: fibre_to_fluid_sign_status
  integer :: action_reaction_status
  integer :: zero_slip_neutrality_status
  integer :: positive_slip_status
  integer :: negative_slip_status
  integer :: structure_side_convention_status
  integer :: fluid_side_convention_status
  integer :: finite_force_status
  integer :: force_norm_finite_status
  integer :: no_eulerian_force_density_status
  integer :: no_rhs_injection_status
  integer :: no_ibm_spreading_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: sign_convention_audit_status
  integer :: structure_equation_plus_fluid_to_fibre_status
  integer :: project_minus_f_fs_means_f_fs_is_fibre_to_fluid_status
  integer :: future_rhs_uses_fibre_to_fluid_status
  real(mytype) :: zero_slip_max_error
  real(mytype) :: fluid_to_fibre_max_error
  real(mytype) :: fibre_to_fluid_max_error
  real(mytype) :: action_reaction_max_error
  real(mytype) :: force_norm_max_error
  real(mytype), allocatable :: u_f(:,:)
  real(mytype), allocatable :: v_f(:,:)
  real(mytype), allocatable :: f_fluid_to_fibre(:,:)
  real(mytype), allocatable :: f_fibre_to_fluid(:,:)
  real(mytype), allocatable :: force_norm(:)

  call execute_command_line('mkdir -p stage12_outputs', exitstat=ios)
  call stage12_config_load()
  requested_flag = logical_status(stage12_requested())
  readonly_mode_status = logical_status(stage12_readonly_mode())
  n_points = stage12_get_max_points()
  if (n_points <= 0) n_points = 8

  allocate(u_f(n_points, 3), v_f(n_points, 3), f_fluid_to_fibre(n_points, 3), &
           f_fibre_to_fluid(n_points, 3), force_norm(n_points), stat=alloc_status)
  if (alloc_status /= 0) then
    call write_failed_allocation(requested_flag, readonly_mode_status)
    print *, 'STAGE 12.4 SIGN CONVENTION AUDIT VERDICT: FAIL'
    print *, 'Reason: allocation failed for sign convention audit arrays.'
    stop 1
  end if

  call stage12_sign_convention_audit_init()

  zero_slip_max_error = run_zero_slip_test(n_points, u_f, v_f, f_fluid_to_fibre, f_fibre_to_fluid, force_norm)
  fluid_to_fibre_max_error = 0.0_mytype
  fibre_to_fluid_max_error = 0.0_mytype
  action_reaction_max_error = 0.0_mytype
  force_norm_max_error = 0.0_mytype

  call run_positive_slip_test(n_points, u_f, v_f, f_fluid_to_fibre, f_fibre_to_fluid, force_norm, &
                              fluid_to_fibre_max_error, fibre_to_fluid_max_error, &
                              action_reaction_max_error, force_norm_max_error)
  positive_slip_status = status_from_errors(fluid_to_fibre_max_error, fibre_to_fluid_max_error, &
                                            sign_formula_abs_tol)

  call run_negative_slip_test(n_points, u_f, v_f, f_fluid_to_fibre, f_fibre_to_fluid, force_norm, &
                              fluid_to_fibre_max_error, fibre_to_fluid_max_error, &
                              action_reaction_max_error, force_norm_max_error)
  negative_slip_status = status_from_errors(fluid_to_fibre_max_error, fibre_to_fluid_max_error, &
                                            sign_formula_abs_tol)

  call run_mixed_slip_test(n_points, u_f, v_f, f_fluid_to_fibre, f_fibre_to_fluid, force_norm, &
                           fluid_to_fibre_max_error, fibre_to_fluid_max_error, &
                           action_reaction_max_error, force_norm_max_error)

  fluid_to_fibre_sign_status = status_from_error(fluid_to_fibre_max_error, sign_formula_abs_tol)
  fibre_to_fluid_sign_status = status_from_error(fibre_to_fluid_max_error, sign_formula_abs_tol)
  action_reaction_status = status_from_error(action_reaction_max_error, action_reaction_abs_tol)
  zero_slip_neutrality_status = status_from_error(zero_slip_max_error, zero_slip_abs_tol)
  if (force_norm_max_error > force_norm_abs_tol) force_norm_finite_status = 0

  structure_equation_plus_fluid_to_fibre_status = 1
  project_minus_f_fs_means_f_fs_is_fibre_to_fluid_status = 1
  future_rhs_uses_fibre_to_fluid_status = 1
  structure_side_convention_status = min(structure_equation_plus_fluid_to_fibre_status, &
                                         project_minus_f_fs_means_f_fs_is_fibre_to_fluid_status)
  fluid_side_convention_status = future_rhs_uses_fibre_to_fluid_status

  call stage12_sign_convention_audit_record_statuses(fluid_to_fibre_sign_status, fibre_to_fluid_sign_status, &
                                                     action_reaction_status, zero_slip_neutrality_status, &
                                                     positive_slip_status, negative_slip_status, &
                                                     structure_side_convention_status, &
                                                     fluid_side_convention_status)

  call stage12_sign_convention_audit_get_status_values(initialized_status, fluid_to_fibre_sign_status, &
                                                       fibre_to_fluid_sign_status, action_reaction_status, &
                                                       zero_slip_neutrality_status, positive_slip_status, &
                                                       negative_slip_status, structure_side_convention_status, &
                                                       fluid_side_convention_status, finite_force_status, &
                                                       force_norm_finite_status, no_eulerian_force_density_status, &
                                                       no_rhs_injection_status, no_ibm_spreading_status, &
                                                       no_feedback_application_status, no_twoway_force_status, &
                                                       no_structure_advance_status, no_fluid_field_access_status, &
                                                       no_fluid_field_modification_status, &
                                                       sign_convention_audit_status)

  if (force_norm_max_error > force_norm_abs_tol) force_norm_finite_status = 0
  if (force_norm_max_error > force_norm_abs_tol) sign_convention_audit_status = 0

  open(newunit=io_unit, file='stage12_outputs/fibre_stage12_4_sign_convention_audit.dat', status='replace', &
       action='write', iostat=ios)
  if (ios /= 0) then
    print *, 'STAGE 12.4 SIGN CONVENTION AUDIT VERDICT: FAIL'
    print *, 'Reason: could not open stage12_outputs/fibre_stage12_4_sign_convention_audit.dat.'
    stop 1
  end if

  call write_int(io_unit, 'stage12_4_requested_flag', requested_flag)
  call write_int(io_unit, 'stage12_4_readonly_mode_status', readonly_mode_status)
  call write_int(io_unit, 'stage12_4_initialized_status', initialized_status)
  call write_int(io_unit, 'stage12_4_fluid_to_fibre_sign_status', fluid_to_fibre_sign_status)
  call write_int(io_unit, 'stage12_4_fibre_to_fluid_sign_status', fibre_to_fluid_sign_status)
  call write_int(io_unit, 'stage12_4_action_reaction_status', action_reaction_status)
  call write_int(io_unit, 'stage12_4_zero_slip_neutrality_status', zero_slip_neutrality_status)
  call write_int(io_unit, 'stage12_4_positive_slip_status', positive_slip_status)
  call write_int(io_unit, 'stage12_4_negative_slip_status', negative_slip_status)
  call write_int(io_unit, 'stage12_4_structure_side_convention_status', structure_side_convention_status)
  call write_int(io_unit, 'stage12_4_fluid_side_convention_status', fluid_side_convention_status)
  call write_int(io_unit, 'stage12_4_finite_force_status', finite_force_status)
  call write_int(io_unit, 'stage12_4_force_norm_finite_status', force_norm_finite_status)
  call write_int(io_unit, 'stage12_4_no_eulerian_force_density_status', no_eulerian_force_density_status)
  call write_int(io_unit, 'stage12_4_no_rhs_injection_status', no_rhs_injection_status)
  call write_int(io_unit, 'stage12_4_no_ibm_spreading_status', no_ibm_spreading_status)
  call write_int(io_unit, 'stage12_4_no_feedback_application_status', no_feedback_application_status)
  call write_int(io_unit, 'stage12_4_no_twoway_force_status', no_twoway_force_status)
  call write_int(io_unit, 'stage12_4_no_structure_advance_status', no_structure_advance_status)
  call write_int(io_unit, 'stage12_4_no_fluid_field_access_status', no_fluid_field_access_status)
  call write_int(io_unit, 'stage12_4_no_fluid_field_modification_status', no_fluid_field_modification_status)
  call write_int(io_unit, 'stage12_4_sign_convention_audit_status', sign_convention_audit_status)
  call write_real(io_unit, 'stage12_4_zero_slip_max_error', zero_slip_max_error)
  call write_real(io_unit, 'stage12_4_fluid_to_fibre_max_error', fluid_to_fibre_max_error)
  call write_real(io_unit, 'stage12_4_fibre_to_fluid_max_error', fibre_to_fluid_max_error)
  call write_real(io_unit, 'stage12_4_action_reaction_max_error', action_reaction_max_error)
  call write_real(io_unit, 'stage12_4_force_norm_max_error', force_norm_max_error)
  call write_int(io_unit, 'stage12_4_structure_equation_plus_fluid_to_fibre_status', &
                 structure_equation_plus_fluid_to_fibre_status)
  call write_int(io_unit, 'stage12_4_project_minus_Ffs_means_Ffs_is_fibre_to_fluid_status', &
                 project_minus_f_fs_means_f_fs_is_fibre_to_fluid_status)
  call write_int(io_unit, 'stage12_4_future_rhs_uses_fibre_to_fluid_status', &
                 future_rhs_uses_fibre_to_fluid_status)
  close(io_unit)

  if (sign_convention_audit_status == 1 .and. requested_flag == 1 .and. readonly_mode_status == 1) then
    print *, 'STAGE 12.4 SIGN CONVENTION AUDIT VERDICT: PASS'
  else
    print *, 'STAGE 12.4 SIGN CONVENTION AUDIT VERDICT: FAIL'
    if (requested_flag /= 1) print *, 'Reason: Stage 12 feedback candidate was not requested.'
    if (readonly_mode_status /= 1) print *, 'Reason: Stage 12 readonly mode was not enforced.'
    if (fluid_to_fibre_sign_status /= 1) print *, 'Reason: fluid-to-fibre sign convention failed.'
    if (fibre_to_fluid_sign_status /= 1) print *, 'Reason: fibre-to-fluid sign convention failed.'
    if (action_reaction_status /= 1) print *, 'Reason: action-reaction consistency failed.'
    if (zero_slip_neutrality_status /= 1) print *, 'Reason: zero-slip neutrality failed.'
    if (positive_slip_status /= 1) print *, 'Reason: positive-slip sign audit failed.'
    if (negative_slip_status /= 1) print *, 'Reason: negative-slip sign audit failed.'
    if (structure_side_convention_status /= 1) print *, 'Reason: structure-side convention was not encoded.'
    if (fluid_side_convention_status /= 1) print *, 'Reason: future fluid-side convention was not encoded.'
    if (finite_force_status /= 1) print *, 'Reason: sign audit force arrays contained non-finite values.'
    if (force_norm_finite_status /= 1) print *, 'Reason: sign audit norm was non-finite or inaccurate.'
    stop 1
  end if

  call stage12_sign_convention_audit_finalize()

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

  integer function status_from_errors(error_a, error_b, tolerance)
    real(mytype), intent(in) :: error_a
    real(mytype), intent(in) :: error_b
    real(mytype), intent(in) :: tolerance
    if (error_a <= tolerance .and. error_b <= tolerance) then
      status_from_errors = 1
    else
      status_from_errors = 0
    end if
  end function status_from_errors

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

    open(newunit=failed_unit, file='stage12_outputs/fibre_stage12_4_sign_convention_audit.dat', &
         status='replace', action='write', iostat=failed_ios)
    if (failed_ios /= 0) return
    call write_int(failed_unit, 'stage12_4_requested_flag', requested_status)
    call write_int(failed_unit, 'stage12_4_readonly_mode_status', readonly_status)
    call write_int(failed_unit, 'stage12_4_initialized_status', 0)
    call write_int(failed_unit, 'stage12_4_sign_convention_audit_status', 0)
    close(failed_unit)
  end subroutine write_failed_allocation

  real(mytype) function run_zero_slip_test(n_local, u_local, v_local, fluid_force, fibre_force, norm_local)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: fluid_force(:,:)
    real(mytype), intent(inout) :: fibre_force(:,:)
    real(mytype), intent(inout) :: norm_local(:)
    integer :: k

    do k = 1, n_local
      u_local(k, :) = (/0.10_mytype + 0.01_mytype * real(k, mytype), -0.20_mytype, 0.05_mytype/)
    end do
    v_local(:, :) = u_local(:, :)
    call stage12_sign_convention_compute_pair(u_local, v_local, 2.0_mytype, fluid_force, fibre_force, norm_local)
    run_zero_slip_test = max(maxval(abs(fluid_force(:, :))), maxval(abs(fibre_force(:, :))))
    run_zero_slip_test = max(run_zero_slip_test, maxval(abs(norm_local(:))))
  end function run_zero_slip_test

  subroutine run_positive_slip_test(n_local, u_local, v_local, fluid_force, fibre_force, norm_local, &
                                    fluid_error, fibre_error, action_error, norm_error)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: fluid_force(:,:)
    real(mytype), intent(inout) :: fibre_force(:,:)
    real(mytype), intent(inout) :: norm_local(:)
    real(mytype), intent(inout) :: fluid_error
    real(mytype), intent(inout) :: fibre_error
    real(mytype), intent(inout) :: action_error
    real(mytype), intent(inout) :: norm_error
    real(mytype) :: slip(3)

    slip = (/0.10_mytype, 0.20_mytype, 0.30_mytype/)
    call fill_constant_slip(n_local, slip, u_local, v_local)
    call stage12_sign_convention_compute_pair(u_local, v_local, 2.0_mytype, fluid_force, fibre_force, norm_local)
    call update_pair_errors(fluid_force, fibre_force, norm_local, 2.0_mytype, slip, fluid_error, fibre_error, &
                            action_error, norm_error)
  end subroutine run_positive_slip_test

  subroutine run_negative_slip_test(n_local, u_local, v_local, fluid_force, fibre_force, norm_local, &
                                    fluid_error, fibre_error, action_error, norm_error)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: fluid_force(:,:)
    real(mytype), intent(inout) :: fibre_force(:,:)
    real(mytype), intent(inout) :: norm_local(:)
    real(mytype), intent(inout) :: fluid_error
    real(mytype), intent(inout) :: fibre_error
    real(mytype), intent(inout) :: action_error
    real(mytype), intent(inout) :: norm_error
    real(mytype) :: slip(3)

    slip = (/-0.15_mytype, -0.05_mytype, -0.25_mytype/)
    call fill_constant_slip(n_local, slip, u_local, v_local)
    call stage12_sign_convention_compute_pair(u_local, v_local, 2.0_mytype, fluid_force, fibre_force, norm_local)
    call update_pair_errors(fluid_force, fibre_force, norm_local, 2.0_mytype, slip, fluid_error, fibre_error, &
                            action_error, norm_error)
  end subroutine run_negative_slip_test

  subroutine run_mixed_slip_test(n_local, u_local, v_local, fluid_force, fibre_force, norm_local, &
                                 fluid_error, fibre_error, action_error, norm_error)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: fluid_force(:,:)
    real(mytype), intent(inout) :: fibre_force(:,:)
    real(mytype), intent(inout) :: norm_local(:)
    real(mytype), intent(inout) :: fluid_error
    real(mytype), intent(inout) :: fibre_error
    real(mytype), intent(inout) :: action_error
    real(mytype), intent(inout) :: norm_error
    real(mytype) :: expected_fluid(3)
    real(mytype) :: expected_fibre(3)
    real(mytype) :: expected_norm
    integer :: k
    integer :: c

    do k = 1, n_local
      v_local(k, :) = (/0.10_mytype, -0.20_mytype, 0.05_mytype/)
      u_local(k, 1) = v_local(k, 1) + 0.02_mytype * real(k, mytype)
      u_local(k, 2) = v_local(k, 2) - 0.03_mytype * real(k, mytype)
      u_local(k, 3) = v_local(k, 3) + 0.01_mytype * real(k, mytype)
    end do
    call stage12_sign_convention_compute_pair(u_local, v_local, 1.5_mytype, fluid_force, fibre_force, norm_local)
    do k = 1, n_local
      expected_fluid(1) = 1.5_mytype * 0.02_mytype * real(k, mytype)
      expected_fluid(2) = -1.5_mytype * 0.03_mytype * real(k, mytype)
      expected_fluid(3) = 1.5_mytype * 0.01_mytype * real(k, mytype)
      expected_fibre(:) = -expected_fluid(:)
      expected_norm = sqrt(sum(expected_fluid(:) * expected_fluid(:)))
      do c = 1, 3
        fluid_error = max(fluid_error, abs(fluid_force(k, c) - expected_fluid(c)))
        fibre_error = max(fibre_error, abs(fibre_force(k, c) - expected_fibre(c)))
        action_error = max(action_error, abs(fluid_force(k, c) + fibre_force(k, c)))
      end do
      norm_error = max(norm_error, abs(norm_local(k) - expected_norm))
    end do
  end subroutine run_mixed_slip_test

  subroutine fill_constant_slip(n_local, slip, u_local, v_local)
    integer, intent(in) :: n_local
    real(mytype), intent(in) :: slip(3)
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    integer :: k

    do k = 1, n_local
      v_local(k, :) = (/0.05_mytype, -0.10_mytype, 0.20_mytype/)
      u_local(k, :) = v_local(k, :) + slip(:)
    end do
  end subroutine fill_constant_slip

  subroutine update_pair_errors(fluid_force, fibre_force, norm_local, alpha, slip, fluid_error, fibre_error, &
                                action_error, norm_error)
    real(mytype), intent(in) :: fluid_force(:,:)
    real(mytype), intent(in) :: fibre_force(:,:)
    real(mytype), intent(in) :: norm_local(:)
    real(mytype), intent(in) :: alpha
    real(mytype), intent(in) :: slip(3)
    real(mytype), intent(inout) :: fluid_error
    real(mytype), intent(inout) :: fibre_error
    real(mytype), intent(inout) :: action_error
    real(mytype), intent(inout) :: norm_error
    real(mytype) :: expected_fluid(3)
    real(mytype) :: expected_fibre(3)
    real(mytype) :: expected_norm
    integer :: k
    integer :: c

    expected_fluid(:) = alpha * slip(:)
    expected_fibre(:) = -expected_fluid(:)
    expected_norm = sqrt(sum(expected_fluid(:) * expected_fluid(:)))
    do k = 1, size(norm_local)
      do c = 1, 3
        fluid_error = max(fluid_error, abs(fluid_force(k, c) - expected_fluid(c)))
        fibre_error = max(fibre_error, abs(fibre_force(k, c) - expected_fibre(c)))
        action_error = max(action_error, abs(fluid_force(k, c) + fibre_force(k, c)))
      end do
      norm_error = max(norm_error, abs(norm_local(k) - expected_norm))
    end do
  end subroutine update_pair_errors

end program fibre_stage12_sign_convention_audit_check
