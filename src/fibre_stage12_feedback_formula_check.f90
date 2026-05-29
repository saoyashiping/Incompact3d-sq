program fibre_stage12_feedback_formula_check
  use fibre_stage12_config, only: stage12_config_load, stage12_requested, stage12_readonly_mode, &
                                  stage12_get_max_points
  use fibre_stage12_force_buffer, only: stage12_force_buffer_init, stage12_force_buffer_finalize
  use fibre_stage12_prescribed_velocity, only: stage12_prescribed_velocity_init, &
                                               stage12_prescribed_velocity_finalize
  use fibre_stage12_feedback_formula, only: stage12_feedback_formula_init, &
                                            stage12_feedback_formula_compute_controlled, &
                                            stage12_feedback_formula_finalize, &
                                            stage12_feedback_formula_get_status_values, &
                                            stage12_feedback_formula_record_test_statuses
  implicit none

  integer, parameter :: mytype = kind(1.0d0)
  real(mytype), parameter :: zero_slip_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: formula_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: gain_scaling_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: force_norm_abs_tol = 1.0e-12_mytype

  integer :: n_points
  integer :: alloc_status
  integer :: ios
  integer :: io_unit
  integer :: requested_flag
  integer :: readonly_mode_status
  integer :: initialized_status
  integer :: zero_slip_status
  integer :: positive_slip_sign_status
  integer :: negative_slip_sign_status
  integer :: force_sign_minus_status
  integer :: multicomponent_slip_status
  integer :: gain_scaling_status
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
  integer :: feedback_formula_status
  real(mytype) :: zero_slip_max_error
  real(mytype) :: positive_slip_max_error
  real(mytype) :: negative_slip_max_error
  real(mytype) :: force_sign_minus_max_error
  real(mytype) :: multicomponent_max_error
  real(mytype) :: gain_scaling_max_error
  real(mytype) :: force_norm_max_error
  real(mytype), allocatable :: u_f(:,:)
  real(mytype), allocatable :: v_f(:,:)
  real(mytype), allocatable :: f_fs_cand(:,:)
  real(mytype), allocatable :: force_norm(:)
  real(mytype), allocatable :: force_alpha1(:,:)
  real(mytype), allocatable :: force_alpha2(:,:)
  real(mytype), allocatable :: norm_alpha1(:)
  real(mytype), allocatable :: norm_alpha2(:)

  call execute_command_line('mkdir -p stage12_outputs', exitstat=ios)
  call stage12_config_load()
  requested_flag = logical_status(stage12_requested())
  readonly_mode_status = logical_status(stage12_readonly_mode())
  n_points = stage12_get_max_points()
  if (n_points <= 0) n_points = 8

  allocate(u_f(n_points, 3), v_f(n_points, 3), f_fs_cand(n_points, 3), force_norm(n_points), &
           force_alpha1(n_points, 3), force_alpha2(n_points, 3), norm_alpha1(n_points), norm_alpha2(n_points), &
           stat=alloc_status)
  if (alloc_status /= 0) then
    call write_failed_allocation(requested_flag, readonly_mode_status)
    print *, 'STAGE 12.3 FEEDBACK FORMULA VERDICT: FAIL'
    print *, 'Reason: allocation failed for controlled feedback formula arrays.'
    stop 1
  end if

  call stage12_force_buffer_init(n_points)
  call stage12_prescribed_velocity_init(n_points)
  call stage12_feedback_formula_init()

  force_norm_max_error = 0.0_mytype
  zero_slip_max_error = run_zero_slip_test(n_points, u_f, v_f, f_fs_cand, force_norm)
  positive_slip_max_error = run_positive_slip_test(n_points, u_f, v_f, f_fs_cand, force_norm)
  negative_slip_max_error = run_negative_slip_test(n_points, u_f, v_f, f_fs_cand, force_norm)
  force_sign_minus_max_error = run_force_sign_minus_test(n_points, u_f, v_f, f_fs_cand, force_norm)
  multicomponent_max_error = run_multicomponent_test(n_points, u_f, v_f, f_fs_cand, force_norm)
  gain_scaling_max_error = run_gain_scaling_test(n_points, u_f, v_f, force_alpha1, force_alpha2, &
                                                 norm_alpha1, norm_alpha2)

  zero_slip_status = status_from_error(zero_slip_max_error, zero_slip_abs_tol)
  positive_slip_sign_status = status_from_error(positive_slip_max_error, formula_abs_tol)
  negative_slip_sign_status = status_from_error(negative_slip_max_error, formula_abs_tol)
  force_sign_minus_status = status_from_error(force_sign_minus_max_error, formula_abs_tol)
  multicomponent_slip_status = status_from_error(multicomponent_max_error, formula_abs_tol)
  gain_scaling_status = status_from_error(gain_scaling_max_error, gain_scaling_abs_tol)

  call stage12_feedback_formula_record_test_statuses(zero_slip_status, positive_slip_sign_status, &
                                                     negative_slip_sign_status, force_sign_minus_status, &
                                                     multicomponent_slip_status, gain_scaling_status)

  call stage12_feedback_formula_get_status_values(initialized_status, zero_slip_status, positive_slip_sign_status, &
                                                  negative_slip_sign_status, force_sign_minus_status, &
                                                  multicomponent_slip_status, gain_scaling_status, &
                                                  finite_force_status, force_norm_finite_status, &
                                                  no_eulerian_force_density_status, no_rhs_injection_status, &
                                                  no_ibm_spreading_status, no_feedback_application_status, &
                                                  no_twoway_force_status, no_structure_advance_status, &
                                                  no_fluid_field_access_status, no_fluid_field_modification_status, &
                                                  feedback_formula_status)

  if (force_norm_max_error > force_norm_abs_tol) force_norm_finite_status = 0
  if (force_norm_max_error > force_norm_abs_tol) feedback_formula_status = 0

  open(newunit=io_unit, file='stage12_outputs/fibre_stage12_3_feedback_formula.dat', status='replace', &
       action='write', iostat=ios)
  if (ios /= 0) then
    print *, 'STAGE 12.3 FEEDBACK FORMULA VERDICT: FAIL'
    print *, 'Reason: could not open stage12_outputs/fibre_stage12_3_feedback_formula.dat.'
    stop 1
  end if

  call write_int(io_unit, 'stage12_3_requested_flag', requested_flag)
  call write_int(io_unit, 'stage12_3_readonly_mode_status', readonly_mode_status)
  call write_int(io_unit, 'stage12_3_initialized_status', initialized_status)
  call write_int(io_unit, 'stage12_3_zero_slip_status', zero_slip_status)
  call write_int(io_unit, 'stage12_3_positive_slip_sign_status', positive_slip_sign_status)
  call write_int(io_unit, 'stage12_3_negative_slip_sign_status', negative_slip_sign_status)
  call write_int(io_unit, 'stage12_3_force_sign_minus_status', force_sign_minus_status)
  call write_int(io_unit, 'stage12_3_multicomponent_slip_status', multicomponent_slip_status)
  call write_int(io_unit, 'stage12_3_gain_scaling_status', gain_scaling_status)
  call write_int(io_unit, 'stage12_3_finite_force_status', finite_force_status)
  call write_int(io_unit, 'stage12_3_force_norm_finite_status', force_norm_finite_status)
  call write_int(io_unit, 'stage12_3_no_eulerian_force_density_status', no_eulerian_force_density_status)
  call write_int(io_unit, 'stage12_3_no_rhs_injection_status', no_rhs_injection_status)
  call write_int(io_unit, 'stage12_3_no_ibm_spreading_status', no_ibm_spreading_status)
  call write_int(io_unit, 'stage12_3_no_feedback_application_status', no_feedback_application_status)
  call write_int(io_unit, 'stage12_3_no_twoway_force_status', no_twoway_force_status)
  call write_int(io_unit, 'stage12_3_no_structure_advance_status', no_structure_advance_status)
  call write_int(io_unit, 'stage12_3_no_fluid_field_access_status', no_fluid_field_access_status)
  call write_int(io_unit, 'stage12_3_no_fluid_field_modification_status', no_fluid_field_modification_status)
  call write_int(io_unit, 'stage12_3_feedback_formula_status', feedback_formula_status)
  call write_real(io_unit, 'stage12_3_zero_slip_max_error', zero_slip_max_error)
  call write_real(io_unit, 'stage12_3_positive_slip_max_error', positive_slip_max_error)
  call write_real(io_unit, 'stage12_3_negative_slip_max_error', negative_slip_max_error)
  call write_real(io_unit, 'stage12_3_force_sign_minus_max_error', force_sign_minus_max_error)
  call write_real(io_unit, 'stage12_3_multicomponent_max_error', multicomponent_max_error)
  call write_real(io_unit, 'stage12_3_gain_scaling_max_error', gain_scaling_max_error)
  call write_real(io_unit, 'stage12_3_force_norm_max_error', force_norm_max_error)
  close(io_unit)

  if (feedback_formula_status == 1 .and. requested_flag == 1 .and. readonly_mode_status == 1) then
    print *, 'STAGE 12.3 FEEDBACK FORMULA VERDICT: PASS'
  else
    print *, 'STAGE 12.3 FEEDBACK FORMULA VERDICT: FAIL'
    if (requested_flag /= 1) print *, 'Reason: Stage 12 feedback candidate was not requested.'
    if (readonly_mode_status /= 1) print *, 'Reason: Stage 12 readonly mode was not enforced.'
    if (feedback_formula_status /= 1) print *, 'Reason: controlled feedback formula status failed.'
    if (zero_slip_status /= 1) print *, 'Reason: zero-slip formula check failed.'
    if (positive_slip_sign_status /= 1) print *, 'Reason: positive-slip formula check failed.'
    if (negative_slip_sign_status /= 1) print *, 'Reason: negative-slip formula check failed.'
    if (force_sign_minus_status /= 1) print *, 'Reason: force_sign=-1 formula check failed.'
    if (multicomponent_slip_status /= 1) print *, 'Reason: multicomponent formula check failed.'
    if (gain_scaling_status /= 1) print *, 'Reason: gain-scaling formula check failed.'
    if (finite_force_status /= 1) print *, 'Reason: force candidate contained non-finite values.'
    if (force_norm_finite_status /= 1) print *, 'Reason: force norm was non-finite or inaccurate.'
    stop 1
  end if

  call stage12_feedback_formula_finalize()
  call stage12_prescribed_velocity_finalize()
  call stage12_force_buffer_finalize()

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

    open(newunit=failed_unit, file='stage12_outputs/fibre_stage12_3_feedback_formula.dat', status='replace', &
         action='write', iostat=failed_ios)
    if (failed_ios /= 0) return
    call write_int(failed_unit, 'stage12_3_requested_flag', requested_status)
    call write_int(failed_unit, 'stage12_3_readonly_mode_status', readonly_status)
    call write_int(failed_unit, 'stage12_3_initialized_status', 0)
    call write_int(failed_unit, 'stage12_3_feedback_formula_status', 0)
    close(failed_unit)
  end subroutine write_failed_allocation

  real(mytype) function run_zero_slip_test(n_local, u_local, v_local, f_local, norm_local)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: f_local(:,:)
    real(mytype), intent(inout) :: norm_local(:)
    integer :: k

    do k = 1, n_local
      u_local(k, 1) = 0.25_mytype + real(k, mytype) * 0.01_mytype
      u_local(k, 2) = -0.10_mytype
      u_local(k, 3) = 0.05_mytype
    end do
    v_local(:, :) = u_local(:, :)
    call stage12_feedback_formula_compute_controlled(u_local, v_local, 2.0_mytype, 1, f_local, norm_local)
    force_norm_max_error = max(force_norm_max_error, compute_force_norm_error(f_local, norm_local))
    run_zero_slip_test = max(maxval(abs(f_local(:, :))), maxval(abs(norm_local(:))))
  end function run_zero_slip_test

  real(mytype) function run_positive_slip_test(n_local, u_local, v_local, f_local, norm_local)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: f_local(:,:)
    real(mytype), intent(inout) :: norm_local(:)
    real(mytype) :: expected(3)
    integer :: k

    expected = (/0.20_mytype, 0.40_mytype, 0.60_mytype/)
    do k = 1, n_local
      v_local(k, :) = 0.0_mytype
      u_local(k, :) = expected(:) / 2.0_mytype
    end do
    call stage12_feedback_formula_compute_controlled(u_local, v_local, 2.0_mytype, 1, f_local, norm_local)
    force_norm_max_error = max(force_norm_max_error, compute_force_norm_error(f_local, norm_local))
    run_positive_slip_test = max_component_error(f_local, expected)
    run_positive_slip_test = max(run_positive_slip_test, force_norm_max_error)
  end function run_positive_slip_test

  real(mytype) function run_negative_slip_test(n_local, u_local, v_local, f_local, norm_local)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: f_local(:,:)
    real(mytype) :: expected(3)
    real(mytype), intent(inout) :: norm_local(:)
    integer :: k

    expected = (/-0.30_mytype, -0.10_mytype, -0.50_mytype/)
    do k = 1, n_local
      v_local(k, :) = 0.0_mytype
      u_local(k, :) = expected(:) / 2.0_mytype
    end do
    call stage12_feedback_formula_compute_controlled(u_local, v_local, 2.0_mytype, 1, f_local, norm_local)
    force_norm_max_error = max(force_norm_max_error, compute_force_norm_error(f_local, norm_local))
    run_negative_slip_test = max_component_error(f_local, expected)
    run_negative_slip_test = max(run_negative_slip_test, force_norm_max_error)
  end function run_negative_slip_test

  real(mytype) function run_force_sign_minus_test(n_local, u_local, v_local, f_local, norm_local)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: f_local(:,:)
    real(mytype), intent(inout) :: norm_local(:)
    real(mytype) :: slip(3)
    real(mytype) :: expected(3)
    integer :: k

    slip = (/0.20_mytype, -0.15_mytype, 0.05_mytype/)
    expected = -2.0_mytype * slip(:)
    do k = 1, n_local
      v_local(k, :) = 0.0_mytype
      u_local(k, :) = slip(:)
    end do
    call stage12_feedback_formula_compute_controlled(u_local, v_local, 2.0_mytype, -1, f_local, norm_local)
    force_norm_max_error = max(force_norm_max_error, compute_force_norm_error(f_local, norm_local))
    run_force_sign_minus_test = max_component_error(f_local, expected)
    run_force_sign_minus_test = max(run_force_sign_minus_test, force_norm_max_error)
  end function run_force_sign_minus_test

  real(mytype) function run_multicomponent_test(n_local, u_local, v_local, f_local, norm_local)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: f_local(:,:)
    real(mytype), intent(inout) :: norm_local(:)
    real(mytype) :: expected(3)
    real(mytype) :: slip(3)
    real(mytype) :: max_error
    integer :: k
    integer :: c

    max_error = 0.0_mytype
    do k = 1, n_local
      slip(1) = 0.05_mytype * real(k, mytype)
      slip(2) = -0.03_mytype * real(k, mytype)
      slip(3) = 0.02_mytype * real(k, mytype)
      v_local(k, :) = (/0.10_mytype, 0.20_mytype, -0.15_mytype/)
      u_local(k, :) = v_local(k, :) + slip(:)
    end do
    call stage12_feedback_formula_compute_controlled(u_local, v_local, 1.5_mytype, 1, f_local, norm_local)
    do k = 1, n_local
      expected(1) = 1.5_mytype * 0.05_mytype * real(k, mytype)
      expected(2) = -1.5_mytype * 0.03_mytype * real(k, mytype)
      expected(3) = 1.5_mytype * 0.02_mytype * real(k, mytype)
      do c = 1, 3
        max_error = max(max_error, abs(f_local(k, c) - expected(c)))
      end do
    end do
    force_norm_max_error = max(force_norm_max_error, compute_force_norm_error(f_local, norm_local))
    run_multicomponent_test = max(max_error, force_norm_max_error)
  end function run_multicomponent_test

  real(mytype) function run_gain_scaling_test(n_local, u_local, v_local, force_one, force_two, norm_one, norm_two)
    integer, intent(in) :: n_local
    real(mytype), intent(inout) :: u_local(:,:)
    real(mytype), intent(inout) :: v_local(:,:)
    real(mytype), intent(inout) :: force_one(:,:)
    real(mytype), intent(inout) :: force_two(:,:)
    real(mytype), intent(inout) :: norm_one(:)
    real(mytype), intent(inout) :: norm_two(:)
    integer :: k

    do k = 1, n_local
      v_local(k, :) = (/0.10_mytype, -0.10_mytype, 0.05_mytype/)
      u_local(k, :) = v_local(k, :) + (/0.03_mytype, -0.04_mytype, 0.07_mytype/)
    end do
    call stage12_feedback_formula_compute_controlled(u_local, v_local, 1.0_mytype, 1, force_one, norm_one)
    call stage12_feedback_formula_compute_controlled(u_local, v_local, 2.0_mytype, 1, force_two, norm_two)
    force_norm_max_error = max(force_norm_max_error, compute_force_norm_error(force_one, norm_one))
    force_norm_max_error = max(force_norm_max_error, compute_force_norm_error(force_two, norm_two))
    run_gain_scaling_test = maxval(abs(force_two(:, :) - 2.0_mytype * force_one(:, :)))
    run_gain_scaling_test = max(run_gain_scaling_test, maxval(abs(norm_two(:) - 2.0_mytype * norm_one(:))))
  end function run_gain_scaling_test

  real(mytype) function max_component_error(f_local, expected)
    real(mytype), intent(in) :: f_local(:,:)
    real(mytype), intent(in) :: expected(3)
    real(mytype) :: max_error
    integer :: k
    integer :: c

    max_error = 0.0_mytype
    do k = 1, size(f_local, 1)
      do c = 1, 3
        max_error = max(max_error, abs(f_local(k, c) - expected(c)))
      end do
    end do
    max_component_error = max_error
  end function max_component_error

  real(mytype) function compute_force_norm_error(f_local, norm_local)
    real(mytype), intent(in) :: f_local(:,:)
    real(mytype), intent(in) :: norm_local(:)
    real(mytype) :: expected_norm
    real(mytype) :: max_error
    integer :: k

    max_error = 0.0_mytype
    do k = 1, size(norm_local)
      expected_norm = sqrt(f_local(k, 1) * f_local(k, 1) + f_local(k, 2) * f_local(k, 2) + &
                           f_local(k, 3) * f_local(k, 3))
      max_error = max(max_error, abs(norm_local(k) - expected_norm))
    end do
    compute_force_norm_error = max_error
  end function compute_force_norm_error

end program fibre_stage12_feedback_formula_check
