program fibre_stage14_rk_timing_audit_check
  use decomp_2d_constants, only : mytype
  use fibre_stage14_config
  use fibre_stage14_rhs_accumulator
  use fibre_stage14_rk_timing_audit
  implicit none

  integer, parameter :: nx = 4
  integer, parameter :: ny = 5
  integer, parameter :: nz = 6
  integer, parameter :: expected_substeps = 3
  real(mytype), parameter :: zero_abs_tol = 1.0e-12_mytype

  real(mytype), allocatable :: fx_cand(:,:,:)
  real(mytype), allocatable :: fy_cand(:,:,:)
  real(mytype), allocatable :: fz_cand(:,:,:)
  real(mytype), allocatable :: rhs_inc_x(:,:,:)
  real(mytype), allocatable :: rhs_inc_y(:,:,:)
  real(mytype), allocatable :: rhs_inc_z(:,:,:)
  real(mytype) :: lambda_zero_max_abs_rhs_increment
  integer :: substep_index
  integer :: i
  integer :: j
  integer :: k
  integer :: requested_flag
  integer :: rhs_injection_enabled_flag
  integer :: injection_gain_finite_status
  integer :: duplicate_detected
  integer :: skipped_detected
  integer :: wrong_order_detected
  integer :: invalid_event_rejected
  integer :: invalid_substep_rejected
  integer :: initialized_status
  integer :: rk_substep_count_status
  integer :: stage13_candidate_before_increment_status
  integer :: increment_before_rhs_addition_status
  integer :: rhs_addition_before_projection_status
  integer :: once_per_substep_status
  integer :: no_duplicate_injection_status
  integer :: no_skipped_substep_status
  integer :: duplicate_detection_status
  integer :: skipped_substep_detection_status
  integer :: wrong_order_detection_status
  integer :: invalid_event_rejection_status
  integer :: invalid_substep_rejection_status
  integer :: valid_event_id_status
  integer :: valid_substep_index_status
  integer :: lambda_zero_increment_status
  integer :: timing_diagnostics_finite_status
  integer :: no_production_rhs_modification_status
  integer :: no_xcompact3d_hook_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: no_pressure_modification_status
  integer :: no_projection_modification_status
  integer :: no_poisson_modification_status
  integer :: no_rk3_modification_status
  integer :: no_channel_forcing_modification_status
  integer :: no_production_ibm_forcing_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: rk_timing_audit_status

  call execute_command_line('mkdir -p stage14_outputs')
  call stage14_config_load()
  requested_flag = logical_to_int_local(stage14_requested())
  rhs_injection_enabled_flag = logical_to_int_local(stage14_rhs_injection_enabled())
  injection_gain_finite_status = logical_to_int_local(finite_real_local(stage14_get_injection_gain()))

  duplicate_detected = run_duplicate_detection_case()
  skipped_detected = run_skipped_substep_detection_case()
  wrong_order_detected = run_wrong_order_detection_case()
  invalid_event_rejected = run_invalid_event_rejection_case()
  invalid_substep_rejected = run_invalid_substep_rejection_case()

  call stage14_rk_timing_audit_init()
  call record_valid_sequence()
  call stage14_rk_timing_audit_check_sequence()

  allocate(fx_cand(nx,ny,nz), fy_cand(nx,ny,nz), fz_cand(nx,ny,nz))
  allocate(rhs_inc_x(nx,ny,nz), rhs_inc_y(nx,ny,nz), rhs_inc_z(nx,ny,nz))
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        fx_cand(i,j,k) = 0.01_mytype * real(i, mytype) - 0.02_mytype * real(j, mytype)
        fy_cand(i,j,k) = -0.03_mytype * real(j, mytype) + 0.01_mytype * real(k, mytype)
        fz_cand(i,j,k) = 0.02_mytype * real(k, mytype) - 0.01_mytype * real(i, mytype)
      end do
    end do
  end do
  call stage14_rhs_accumulator_init(nx, ny, nz)
  call stage14_rhs_accumulator_compute_from_force_density(fx_cand, fy_cand, fz_cand, 0.0_mytype)
  rhs_inc_x = 0.0_mytype * fx_cand
  rhs_inc_y = 0.0_mytype * fy_cand
  rhs_inc_z = 0.0_mytype * fz_cand
  lambda_zero_max_abs_rhs_increment = max(maxval(abs(rhs_inc_x)), max(maxval(abs(rhs_inc_y)), maxval(abs(rhs_inc_z))))
  call stage14_rk_timing_audit_note_lambda_zero_increment(lambda_zero_max_abs_rhs_increment)

  call stage14_rk_timing_audit_set_negative_detection_statuses(duplicate_detected, skipped_detected, &
                                                               wrong_order_detected, invalid_event_rejected, &
                                                               invalid_substep_rejected)

  call stage14_rk_timing_audit_get_status_values(initialized_status, &
                                                 rk_substep_count_status, &
                                                 stage13_candidate_before_increment_status, &
                                                 increment_before_rhs_addition_status, &
                                                 rhs_addition_before_projection_status, &
                                                 once_per_substep_status, &
                                                 no_duplicate_injection_status, &
                                                 no_skipped_substep_status, &
                                                 duplicate_detection_status, &
                                                 skipped_substep_detection_status, &
                                                 wrong_order_detection_status, &
                                                 invalid_event_rejection_status, &
                                                 invalid_substep_rejection_status, &
                                                 valid_event_id_status, &
                                                 valid_substep_index_status, &
                                                 lambda_zero_increment_status, &
                                                 timing_diagnostics_finite_status, &
                                                 no_production_rhs_modification_status, &
                                                 no_xcompact3d_hook_status, &
                                                 no_fluid_field_access_status, &
                                                 no_fluid_field_modification_status, &
                                                 no_pressure_modification_status, &
                                                 no_projection_modification_status, &
                                                 no_poisson_modification_status, &
                                                 no_rk3_modification_status, &
                                                 no_channel_forcing_modification_status, &
                                                 no_production_ibm_forcing_status, &
                                                 no_feedback_application_status, &
                                                 no_twoway_force_status, &
                                                 no_structure_advance_status, &
                                                 rk_timing_audit_status)

  call stage14_rk_timing_audit_write_diagnostics('stage14_outputs/fibre_stage14_4_rk_timing_audit.dat', &
                                                 requested_flag, rhs_injection_enabled_flag, &
                                                 injection_gain_finite_status)

  if (rk_timing_audit_status == 1 .and. requested_flag == 1 .and. &
      rhs_injection_enabled_flag == 1 .and. injection_gain_finite_status == 1 .and. &
      lambda_zero_max_abs_rhs_increment <= zero_abs_tol) then
    print *, 'STAGE 14.4 RK TIMING AUDIT VERDICT: PASS'
  else
    print *, 'STAGE 14.4 RK TIMING AUDIT VERDICT: FAIL'
    if (requested_flag /= 1) print *, 'Reason: stage14_4_requested_flag'
    if (rhs_injection_enabled_flag /= 1) print *, 'Reason: stage14_4_rhs_injection_enabled_flag'
    if (injection_gain_finite_status /= 1) print *, 'Reason: stage14_4_injection_gain_finite_status'
    if (initialized_status /= 1) print *, 'Reason: stage14_4_initialized_status'
    if (rk_substep_count_status /= 1) print *, 'Reason: stage14_4_rk_substep_count_status'
    if (stage13_candidate_before_increment_status /= 1) print *, 'Reason: stage14_4_stage13_candidate_before_increment_status'
    if (increment_before_rhs_addition_status /= 1) print *, 'Reason: stage14_4_increment_before_rhs_addition_status'
    if (rhs_addition_before_projection_status /= 1) print *, 'Reason: stage14_4_rhs_addition_before_projection_status'
    if (once_per_substep_status /= 1) print *, 'Reason: stage14_4_once_per_substep_status'
    if (no_duplicate_injection_status /= 1) print *, 'Reason: stage14_4_no_duplicate_injection_status'
    if (no_skipped_substep_status /= 1) print *, 'Reason: stage14_4_no_skipped_substep_status'
    if (duplicate_detection_status /= 1) print *, 'Reason: stage14_4_duplicate_detection_status'
    if (skipped_substep_detection_status /= 1) print *, 'Reason: stage14_4_skipped_substep_detection_status'
    if (wrong_order_detection_status /= 1) print *, 'Reason: stage14_4_wrong_order_detection_status'
    if (invalid_event_rejection_status /= 1) print *, 'Reason: stage14_4_invalid_event_rejection_status'
    if (invalid_substep_rejection_status /= 1) print *, 'Reason: stage14_4_invalid_substep_rejection_status'
    if (lambda_zero_increment_status /= 1) print *, 'Reason: stage14_4_lambda_zero_increment_status'
    if (timing_diagnostics_finite_status /= 1) print *, 'Reason: stage14_4_timing_diagnostics_finite_status'
    if (no_production_rhs_modification_status /= 1) print *, 'Reason: stage14_4_no_production_rhs_modification_status'
    if (no_xcompact3d_hook_status /= 1) print *, 'Reason: stage14_4_no_xcompact3d_hook_status'
    if (no_fluid_field_access_status /= 1) print *, 'Reason: stage14_4_no_fluid_field_access_status'
    if (no_fluid_field_modification_status /= 1) print *, 'Reason: stage14_4_no_fluid_field_modification_status'
    if (no_pressure_modification_status /= 1) print *, 'Reason: stage14_4_no_pressure_modification_status'
    if (no_projection_modification_status /= 1) print *, 'Reason: stage14_4_no_projection_modification_status'
    if (no_poisson_modification_status /= 1) print *, 'Reason: stage14_4_no_poisson_modification_status'
    if (no_rk3_modification_status /= 1) print *, 'Reason: stage14_4_no_rk3_modification_status'
    if (no_channel_forcing_modification_status /= 1) print *, 'Reason: stage14_4_no_channel_forcing_modification_status'
    if (no_production_ibm_forcing_status /= 1) print *, 'Reason: stage14_4_no_production_ibm_forcing_status'
    if (no_feedback_application_status /= 1) print *, 'Reason: stage14_4_no_feedback_application_status'
    if (no_twoway_force_status /= 1) print *, 'Reason: stage14_4_no_twoway_force_status'
    if (no_structure_advance_status /= 1) print *, 'Reason: stage14_4_no_structure_advance_status'
    if (rk_timing_audit_status /= 1) print *, 'Reason: stage14_4_rk_timing_audit_status'
    if (lambda_zero_max_abs_rhs_increment > zero_abs_tol) print *, 'Reason: stage14_4_lambda_zero_max_abs_rhs_increment'
    stop 1
  end if

  call stage14_rhs_accumulator_finalize()
  call stage14_rk_timing_audit_finalize()
  deallocate(fx_cand, fy_cand, fz_cand)
  deallocate(rhs_inc_x, rhs_inc_y, rhs_inc_z)

contains

  subroutine record_valid_sequence()
    integer :: substep
    do substep = 1, expected_substeps
      call stage14_rk_timing_audit_record_event(substep, EVENT_STAGE13_CANDIDATE_READY)
      call stage14_rk_timing_audit_record_event(substep, EVENT_STAGE14_INCREMENT_READY)
      call stage14_rk_timing_audit_record_event(substep, EVENT_RHS_ADDITION_PLANNED)
      call stage14_rk_timing_audit_record_event(substep, EVENT_PROJECTION_PLANNED)
    end do
  end subroutine record_valid_sequence

  integer function run_duplicate_detection_case()
    integer :: dup_initialized
    integer :: dup_rk_count
    integer :: dup_stage13_before_increment
    integer :: dup_increment_before_rhs
    integer :: dup_rhs_before_projection
    integer :: dup_once
    integer :: dup_no_duplicate
    integer :: dup_no_skipped
    integer :: ignored_status(24)

    call stage14_rk_timing_audit_init()
    call record_valid_sequence()
    call stage14_rk_timing_audit_record_event(2, EVENT_RHS_ADDITION_PLANNED)
    call stage14_rk_timing_audit_check_sequence()
    call get_subset_status(dup_initialized, dup_rk_count, dup_stage13_before_increment, &
                           dup_increment_before_rhs, dup_rhs_before_projection, dup_once, &
                           dup_no_duplicate, dup_no_skipped, ignored_status(1), ignored_status(2))
    run_duplicate_detection_case = logical_to_int_local(dup_no_duplicate == 0)
  end function run_duplicate_detection_case

  integer function run_skipped_substep_detection_case()
    integer :: skip_initialized
    integer :: skip_rk_count
    integer :: skip_stage13_before_increment
    integer :: skip_increment_before_rhs
    integer :: skip_rhs_before_projection
    integer :: skip_once
    integer :: skip_no_duplicate
    integer :: skip_no_skipped
    integer :: ignored_status(24)

    call stage14_rk_timing_audit_init()
    call stage14_rk_timing_audit_record_event(1, EVENT_STAGE13_CANDIDATE_READY)
    call stage14_rk_timing_audit_record_event(1, EVENT_STAGE14_INCREMENT_READY)
    call stage14_rk_timing_audit_record_event(1, EVENT_RHS_ADDITION_PLANNED)
    call stage14_rk_timing_audit_record_event(1, EVENT_PROJECTION_PLANNED)
    call stage14_rk_timing_audit_record_event(2, EVENT_STAGE13_CANDIDATE_READY)
    call stage14_rk_timing_audit_record_event(2, EVENT_STAGE14_INCREMENT_READY)
    call stage14_rk_timing_audit_record_event(2, EVENT_PROJECTION_PLANNED)
    call stage14_rk_timing_audit_record_event(3, EVENT_STAGE13_CANDIDATE_READY)
    call stage14_rk_timing_audit_record_event(3, EVENT_STAGE14_INCREMENT_READY)
    call stage14_rk_timing_audit_record_event(3, EVENT_RHS_ADDITION_PLANNED)
    call stage14_rk_timing_audit_record_event(3, EVENT_PROJECTION_PLANNED)
    call stage14_rk_timing_audit_check_sequence()
    call get_subset_status(skip_initialized, skip_rk_count, skip_stage13_before_increment, &
                           skip_increment_before_rhs, skip_rhs_before_projection, skip_once, &
                           skip_no_duplicate, skip_no_skipped, ignored_status(1), ignored_status(2))
    run_skipped_substep_detection_case = logical_to_int_local(skip_no_skipped == 0)
  end function run_skipped_substep_detection_case

  integer function run_wrong_order_detection_case()
    integer :: wrong_initialized
    integer :: wrong_rk_count
    integer :: wrong_stage13_before_increment
    integer :: wrong_increment_before_rhs
    integer :: wrong_rhs_before_projection
    integer :: wrong_once
    integer :: wrong_no_duplicate
    integer :: wrong_no_skipped
    integer :: ignored_status(24)

    call stage14_rk_timing_audit_init()
    call stage14_rk_timing_audit_record_event(1, EVENT_STAGE13_CANDIDATE_READY)
    call stage14_rk_timing_audit_record_event(1, EVENT_STAGE14_INCREMENT_READY)
    call stage14_rk_timing_audit_record_event(1, EVENT_PROJECTION_PLANNED)
    call stage14_rk_timing_audit_record_event(1, EVENT_RHS_ADDITION_PLANNED)
    do substep_index = 2, expected_substeps
      call stage14_rk_timing_audit_record_event(substep_index, EVENT_STAGE13_CANDIDATE_READY)
      call stage14_rk_timing_audit_record_event(substep_index, EVENT_STAGE14_INCREMENT_READY)
      call stage14_rk_timing_audit_record_event(substep_index, EVENT_RHS_ADDITION_PLANNED)
      call stage14_rk_timing_audit_record_event(substep_index, EVENT_PROJECTION_PLANNED)
    end do
    call stage14_rk_timing_audit_check_sequence()
    call get_subset_status(wrong_initialized, wrong_rk_count, wrong_stage13_before_increment, &
                           wrong_increment_before_rhs, wrong_rhs_before_projection, wrong_once, &
                           wrong_no_duplicate, wrong_no_skipped, ignored_status(1), ignored_status(2))
    run_wrong_order_detection_case = logical_to_int_local(wrong_rhs_before_projection == 0)
  end function run_wrong_order_detection_case

  integer function run_invalid_event_rejection_case()
    integer :: invalid_event_id
    integer :: invalid_substep_id
    integer :: ignored_status(24)

    call stage14_rk_timing_audit_init()
    call stage14_rk_timing_audit_record_event(1, 99)
    call get_subset_status(ignored_status(1), ignored_status(2), ignored_status(3), ignored_status(4), &
                           ignored_status(5), ignored_status(6), ignored_status(7), ignored_status(8), &
                           invalid_event_id, invalid_substep_id)
    run_invalid_event_rejection_case = logical_to_int_local(invalid_event_id == 0)
  end function run_invalid_event_rejection_case

  integer function run_invalid_substep_rejection_case()
    integer :: invalid_event_id
    integer :: invalid_substep_id
    integer :: ignored_status(24)

    call stage14_rk_timing_audit_init()
    call stage14_rk_timing_audit_record_event(0, EVENT_STAGE13_CANDIDATE_READY)
    call get_subset_status(ignored_status(1), ignored_status(2), ignored_status(3), ignored_status(4), &
                           ignored_status(5), ignored_status(6), ignored_status(7), ignored_status(8), &
                           invalid_event_id, invalid_substep_id)
    run_invalid_substep_rejection_case = logical_to_int_local(invalid_substep_id == 0)
  end function run_invalid_substep_rejection_case

  subroutine get_subset_status(initialized_out, rk_count_out, stage13_before_increment_out, increment_before_rhs_out, &
                               rhs_before_projection_out, once_out, no_duplicate_out, no_skipped_out, &
                               valid_event_out, valid_substep_out)
    integer, intent(out) :: initialized_out
    integer, intent(out) :: rk_count_out
    integer, intent(out) :: stage13_before_increment_out
    integer, intent(out) :: increment_before_rhs_out
    integer, intent(out) :: rhs_before_projection_out
    integer, intent(out) :: once_out
    integer, intent(out) :: no_duplicate_out
    integer, intent(out) :: no_skipped_out
    integer, intent(out) :: valid_event_out
    integer, intent(out) :: valid_substep_out
    integer :: duplicate_detection_out
    integer :: skipped_detection_out
    integer :: wrong_order_detection_out
    integer :: invalid_event_rejection_out
    integer :: invalid_substep_rejection_out
    integer :: lambda_zero_out
    integer :: finite_out
    integer :: no_rhs_out
    integer :: no_hook_out
    integer :: no_fluid_access_out
    integer :: no_fluid_mod_out
    integer :: no_pressure_out
    integer :: no_projection_out
    integer :: no_poisson_out
    integer :: no_rk3_out
    integer :: no_channel_out
    integer :: no_ibm_out
    integer :: no_feedback_out
    integer :: no_twoway_out
    integer :: no_structure_out
    integer :: status_out

    call stage14_rk_timing_audit_get_status_values(initialized_out, rk_count_out, &
                                                   stage13_before_increment_out, increment_before_rhs_out, &
                                                   rhs_before_projection_out, once_out, no_duplicate_out, &
                                                   no_skipped_out, duplicate_detection_out, skipped_detection_out, &
                                                   wrong_order_detection_out, invalid_event_rejection_out, &
                                                   invalid_substep_rejection_out, valid_event_out, valid_substep_out, &
                                                   lambda_zero_out, finite_out, no_rhs_out, no_hook_out, &
                                                   no_fluid_access_out, no_fluid_mod_out, no_pressure_out, &
                                                   no_projection_out, no_poisson_out, no_rk3_out, no_channel_out, &
                                                   no_ibm_out, no_feedback_out, no_twoway_out, no_structure_out, &
                                                   status_out)
  end subroutine get_subset_status

  logical function finite_real_local(value)
    real(mytype), intent(in) :: value
    finite_real_local = (value == value) .and. (abs(value) < huge(value))
  end function finite_real_local

  integer function logical_to_int_local(value)
    logical, intent(in) :: value
    if (value) then
      logical_to_int_local = 1
    else
      logical_to_int_local = 0
    end if
  end function logical_to_int_local

end program fibre_stage14_rk_timing_audit_check
