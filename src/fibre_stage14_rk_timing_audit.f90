module fibre_stage14_rk_timing_audit
  use decomp_2d_constants, only : mytype
  implicit none

  private

  integer, parameter, public :: EVENT_STAGE13_CANDIDATE_READY = 1
  integer, parameter, public :: EVENT_STAGE14_INCREMENT_READY = 2
  integer, parameter, public :: EVENT_RHS_ADDITION_PLANNED = 3
  integer, parameter, public :: EVENT_PROJECTION_PLANNED = 4

  public :: stage14_rk_timing_audit_init
  public :: stage14_rk_timing_audit_clear
  public :: stage14_rk_timing_audit_finalize
  public :: stage14_rk_timing_audit_record_event
  public :: stage14_rk_timing_audit_check_sequence
  public :: stage14_rk_timing_audit_get_status_values
  public :: stage14_rk_timing_audit_write_diagnostics
  public :: stage14_rk_timing_audit_note_lambda_zero_increment
  public :: stage14_rk_timing_audit_set_negative_detection_statuses

  integer, parameter :: max_substeps = 3
  integer, parameter :: max_events = 4
  real(mytype), parameter :: lambda_zero_abs_tol = 1.0e-12_mytype

  integer :: expected_nsubsteps = max_substeps
  integer :: event_count(max_substeps,max_events) = 0
  integer :: first_event_order(max_substeps,max_events) = 0
  integer :: event_order_counter = 0

  integer :: initialized_status = 0
  integer :: rk_substep_count_status = 0
  integer :: stage13_candidate_before_increment_status = 0
  integer :: increment_before_rhs_addition_status = 0
  integer :: rhs_addition_before_projection_status = 0
  integer :: once_per_substep_status = 0
  integer :: no_duplicate_injection_status = 0
  integer :: no_skipped_substep_status = 0
  integer :: valid_event_id_status = 1
  integer :: valid_substep_index_status = 1
  integer :: duplicate_detection_status = 0
  integer :: skipped_substep_detection_status = 0
  integer :: wrong_order_detection_status = 0
  integer :: invalid_event_rejection_status = 0
  integer :: invalid_substep_rejection_status = 0
  integer :: lambda_zero_increment_status = 0
  integer :: timing_diagnostics_finite_status = 1
  integer :: no_production_rhs_modification_status = 1
  integer :: no_xcompact3d_hook_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: no_pressure_modification_status = 1
  integer :: no_projection_modification_status = 1
  integer :: no_poisson_modification_status = 1
  integer :: no_rk3_modification_status = 1
  integer :: no_channel_forcing_modification_status = 1
  integer :: no_production_ibm_forcing_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: rk_timing_audit_status = 0

  integer :: recorded_rhs_addition_count = 0
  integer :: duplicate_injection_count = 0
  integer :: skipped_substep_count = 0
  integer :: wrong_order_count = 0
  integer :: invalid_event_count = 0
  integer :: invalid_substep_count = 0
  real(mytype) :: lambda_zero_max_abs_rhs_increment = 0.0_mytype

contains

  subroutine stage14_rk_timing_audit_init()
    initialized_status = 1
    call reset_sequence_state()
    duplicate_detection_status = 0
    skipped_substep_detection_status = 0
    wrong_order_detection_status = 0
    invalid_event_rejection_status = 0
    invalid_substep_rejection_status = 0
    lambda_zero_increment_status = 0
    timing_diagnostics_finite_status = 1
    no_production_rhs_modification_status = 1
    no_xcompact3d_hook_status = 1
    no_fluid_field_access_status = 1
    no_fluid_field_modification_status = 1
    no_pressure_modification_status = 1
    no_projection_modification_status = 1
    no_poisson_modification_status = 1
    no_rk3_modification_status = 1
    no_channel_forcing_modification_status = 1
    no_production_ibm_forcing_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    lambda_zero_max_abs_rhs_increment = 0.0_mytype
    call update_timing_finite_status()
    call update_overall_status()
  end subroutine stage14_rk_timing_audit_init

  subroutine stage14_rk_timing_audit_clear()
    call reset_sequence_state()
    call update_overall_status()
  end subroutine stage14_rk_timing_audit_clear

  subroutine stage14_rk_timing_audit_finalize()
    initialized_status = 0
    rk_timing_audit_status = 0
  end subroutine stage14_rk_timing_audit_finalize

  subroutine stage14_rk_timing_audit_record_event(substep_index, event_id)
    integer, intent(in) :: substep_index
    integer, intent(in) :: event_id

    if (substep_index < 1 .or. substep_index > expected_nsubsteps) then
      invalid_substep_count = invalid_substep_count + 1
      valid_substep_index_status = 0
      call update_timing_finite_status()
      call update_overall_status()
      return
    end if

    if (event_id < EVENT_STAGE13_CANDIDATE_READY .or. event_id > EVENT_PROJECTION_PLANNED) then
      invalid_event_count = invalid_event_count + 1
      valid_event_id_status = 0
      call update_timing_finite_status()
      call update_overall_status()
      return
    end if

    event_order_counter = event_order_counter + 1
    event_count(substep_index,event_id) = event_count(substep_index,event_id) + 1
    if (first_event_order(substep_index,event_id) == 0) then
      first_event_order(substep_index,event_id) = event_order_counter
    end if
    call update_timing_finite_status()
    call update_overall_status()
  end subroutine stage14_rk_timing_audit_record_event

  subroutine stage14_rk_timing_audit_check_sequence()
    integer :: substep_index
    integer :: local_duplicate_count
    integer :: local_skipped_count
    integer :: local_wrong_order_count
    integer :: local_recorded_rhs_count
    integer :: local_once_status
    integer :: candidate_before_increment
    integer :: increment_before_rhs
    integer :: rhs_before_projection

    local_duplicate_count = 0
    local_skipped_count = 0
    local_wrong_order_count = 0
    local_recorded_rhs_count = 0
    local_once_status = 1
    candidate_before_increment = 1
    increment_before_rhs = 1
    rhs_before_projection = 1

    rk_substep_count_status = logical_to_int(expected_nsubsteps == max_substeps)

    do substep_index = 1, expected_nsubsteps
      local_recorded_rhs_count = local_recorded_rhs_count + event_count(substep_index,EVENT_RHS_ADDITION_PLANNED)

      if (event_count(substep_index,EVENT_RHS_ADDITION_PLANNED) > 1) then
        local_duplicate_count = local_duplicate_count + event_count(substep_index,EVENT_RHS_ADDITION_PLANNED) - 1
      end if
      if (event_count(substep_index,EVENT_RHS_ADDITION_PLANNED) == 0) then
        local_skipped_count = local_skipped_count + 1
      end if

      if (.not. all_required_events_once(substep_index)) local_once_status = 0

      if (.not. ordered_before(substep_index, EVENT_STAGE13_CANDIDATE_READY, EVENT_STAGE14_INCREMENT_READY)) then
        candidate_before_increment = 0
        local_wrong_order_count = local_wrong_order_count + 1
      end if
      if (.not. ordered_before(substep_index, EVENT_STAGE14_INCREMENT_READY, EVENT_RHS_ADDITION_PLANNED)) then
        increment_before_rhs = 0
        local_wrong_order_count = local_wrong_order_count + 1
      end if
      if (.not. ordered_before(substep_index, EVENT_RHS_ADDITION_PLANNED, EVENT_PROJECTION_PLANNED)) then
        rhs_before_projection = 0
        local_wrong_order_count = local_wrong_order_count + 1
      end if
    end do

    recorded_rhs_addition_count = local_recorded_rhs_count
    duplicate_injection_count = local_duplicate_count
    skipped_substep_count = local_skipped_count
    wrong_order_count = local_wrong_order_count
    stage13_candidate_before_increment_status = candidate_before_increment
    increment_before_rhs_addition_status = increment_before_rhs
    rhs_addition_before_projection_status = rhs_before_projection
    once_per_substep_status = local_once_status
    no_duplicate_injection_status = logical_to_int(duplicate_injection_count == 0)
    no_skipped_substep_status = logical_to_int(skipped_substep_count == 0)
    call update_timing_finite_status()
    call update_overall_status()
  end subroutine stage14_rk_timing_audit_check_sequence

  subroutine stage14_rk_timing_audit_note_lambda_zero_increment(max_abs_rhs_increment)
    real(mytype), intent(in) :: max_abs_rhs_increment

    lambda_zero_max_abs_rhs_increment = max_abs_rhs_increment
    lambda_zero_increment_status = logical_to_int(finite_real(max_abs_rhs_increment) .and. &
                                                  abs(max_abs_rhs_increment) <= lambda_zero_abs_tol)
    call update_timing_finite_status()
    call update_overall_status()
  end subroutine stage14_rk_timing_audit_note_lambda_zero_increment

  subroutine stage14_rk_timing_audit_set_negative_detection_statuses(duplicate_detected, skipped_detected, &
                                                                     wrong_order_detected, invalid_event_rejected, &
                                                                     invalid_substep_rejected)
    integer, intent(in) :: duplicate_detected
    integer, intent(in) :: skipped_detected
    integer, intent(in) :: wrong_order_detected
    integer, intent(in) :: invalid_event_rejected
    integer, intent(in) :: invalid_substep_rejected

    duplicate_detection_status = logical_to_int(duplicate_detected == 1)
    skipped_substep_detection_status = logical_to_int(skipped_detected == 1)
    wrong_order_detection_status = logical_to_int(wrong_order_detected == 1)
    invalid_event_rejection_status = logical_to_int(invalid_event_rejected == 1)
    invalid_substep_rejection_status = logical_to_int(invalid_substep_rejected == 1)
    call update_overall_status()
  end subroutine stage14_rk_timing_audit_set_negative_detection_statuses

  subroutine stage14_rk_timing_audit_get_status_values(initialized_status_out, &
                                                       rk_substep_count_status_out, &
                                                       stage13_candidate_before_increment_status_out, &
                                                       increment_before_rhs_addition_status_out, &
                                                       rhs_addition_before_projection_status_out, &
                                                       once_per_substep_status_out, &
                                                       no_duplicate_injection_status_out, &
                                                       no_skipped_substep_status_out, &
                                                       duplicate_detection_status_out, &
                                                       skipped_substep_detection_status_out, &
                                                       wrong_order_detection_status_out, &
                                                       invalid_event_rejection_status_out, &
                                                       invalid_substep_rejection_status_out, &
                                                       valid_event_id_status_out, &
                                                       valid_substep_index_status_out, &
                                                       lambda_zero_increment_status_out, &
                                                       timing_diagnostics_finite_status_out, &
                                                       no_production_rhs_modification_status_out, &
                                                       no_xcompact3d_hook_status_out, &
                                                       no_fluid_field_access_status_out, &
                                                       no_fluid_field_modification_status_out, &
                                                       no_pressure_modification_status_out, &
                                                       no_projection_modification_status_out, &
                                                       no_poisson_modification_status_out, &
                                                       no_rk3_modification_status_out, &
                                                       no_channel_forcing_modification_status_out, &
                                                       no_production_ibm_forcing_status_out, &
                                                       no_feedback_application_status_out, &
                                                       no_twoway_force_status_out, &
                                                       no_structure_advance_status_out, &
                                                       rk_timing_audit_status_out)
    integer, intent(out) :: initialized_status_out
    integer, intent(out) :: rk_substep_count_status_out
    integer, intent(out) :: stage13_candidate_before_increment_status_out
    integer, intent(out) :: increment_before_rhs_addition_status_out
    integer, intent(out) :: rhs_addition_before_projection_status_out
    integer, intent(out) :: once_per_substep_status_out
    integer, intent(out) :: no_duplicate_injection_status_out
    integer, intent(out) :: no_skipped_substep_status_out
    integer, intent(out) :: duplicate_detection_status_out
    integer, intent(out) :: skipped_substep_detection_status_out
    integer, intent(out) :: wrong_order_detection_status_out
    integer, intent(out) :: invalid_event_rejection_status_out
    integer, intent(out) :: invalid_substep_rejection_status_out
    integer, intent(out) :: valid_event_id_status_out
    integer, intent(out) :: valid_substep_index_status_out
    integer, intent(out) :: lambda_zero_increment_status_out
    integer, intent(out) :: timing_diagnostics_finite_status_out
    integer, intent(out) :: no_production_rhs_modification_status_out
    integer, intent(out) :: no_xcompact3d_hook_status_out
    integer, intent(out) :: no_fluid_field_access_status_out
    integer, intent(out) :: no_fluid_field_modification_status_out
    integer, intent(out) :: no_pressure_modification_status_out
    integer, intent(out) :: no_projection_modification_status_out
    integer, intent(out) :: no_poisson_modification_status_out
    integer, intent(out) :: no_rk3_modification_status_out
    integer, intent(out) :: no_channel_forcing_modification_status_out
    integer, intent(out) :: no_production_ibm_forcing_status_out
    integer, intent(out) :: no_feedback_application_status_out
    integer, intent(out) :: no_twoway_force_status_out
    integer, intent(out) :: no_structure_advance_status_out
    integer, intent(out) :: rk_timing_audit_status_out

    initialized_status_out = initialized_status
    rk_substep_count_status_out = rk_substep_count_status
    stage13_candidate_before_increment_status_out = stage13_candidate_before_increment_status
    increment_before_rhs_addition_status_out = increment_before_rhs_addition_status
    rhs_addition_before_projection_status_out = rhs_addition_before_projection_status
    once_per_substep_status_out = once_per_substep_status
    no_duplicate_injection_status_out = no_duplicate_injection_status
    no_skipped_substep_status_out = no_skipped_substep_status
    duplicate_detection_status_out = duplicate_detection_status
    skipped_substep_detection_status_out = skipped_substep_detection_status
    wrong_order_detection_status_out = wrong_order_detection_status
    invalid_event_rejection_status_out = invalid_event_rejection_status
    invalid_substep_rejection_status_out = invalid_substep_rejection_status
    valid_event_id_status_out = valid_event_id_status
    valid_substep_index_status_out = valid_substep_index_status
    lambda_zero_increment_status_out = lambda_zero_increment_status
    timing_diagnostics_finite_status_out = timing_diagnostics_finite_status
    no_production_rhs_modification_status_out = no_production_rhs_modification_status
    no_xcompact3d_hook_status_out = no_xcompact3d_hook_status
    no_fluid_field_access_status_out = no_fluid_field_access_status
    no_fluid_field_modification_status_out = no_fluid_field_modification_status
    no_pressure_modification_status_out = no_pressure_modification_status
    no_projection_modification_status_out = no_projection_modification_status
    no_poisson_modification_status_out = no_poisson_modification_status
    no_rk3_modification_status_out = no_rk3_modification_status
    no_channel_forcing_modification_status_out = no_channel_forcing_modification_status
    no_production_ibm_forcing_status_out = no_production_ibm_forcing_status
    no_feedback_application_status_out = no_feedback_application_status
    no_twoway_force_status_out = no_twoway_force_status
    no_structure_advance_status_out = no_structure_advance_status
    rk_timing_audit_status_out = rk_timing_audit_status
  end subroutine stage14_rk_timing_audit_get_status_values

  subroutine stage14_rk_timing_audit_write_diagnostics(filename, requested_flag, rhs_injection_enabled_flag, &
                                                       injection_gain_finite_status)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: requested_flag
    integer, intent(in) :: rhs_injection_enabled_flag
    integer, intent(in) :: injection_gain_finite_status
    integer :: unit_id
    integer :: io_status

    open(newunit=unit_id, file=filename, status='replace', action='write', iostat=io_status)
    if (io_status /= 0) return

    write(unit_id, '(A,1X,I0)') 'stage14_4_requested_flag', requested_flag
    write(unit_id, '(A,1X,I0)') 'stage14_4_rhs_injection_enabled_flag', rhs_injection_enabled_flag
    write(unit_id, '(A,1X,I0)') 'stage14_4_injection_gain_finite_status', injection_gain_finite_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_initialized_status', initialized_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_rk_substep_count_status', rk_substep_count_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_stage13_candidate_before_increment_status', &
                                 stage13_candidate_before_increment_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_increment_before_rhs_addition_status', &
                                 increment_before_rhs_addition_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_rhs_addition_before_projection_status', &
                                 rhs_addition_before_projection_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_once_per_substep_status', once_per_substep_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_duplicate_injection_status', no_duplicate_injection_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_skipped_substep_status', no_skipped_substep_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_duplicate_detection_status', duplicate_detection_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_skipped_substep_detection_status', skipped_substep_detection_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_wrong_order_detection_status', wrong_order_detection_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_invalid_event_rejection_status', invalid_event_rejection_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_invalid_substep_rejection_status', invalid_substep_rejection_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_lambda_zero_increment_status', lambda_zero_increment_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_timing_diagnostics_finite_status', timing_diagnostics_finite_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_production_rhs_modification_status', &
                                 no_production_rhs_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_xcompact3d_hook_status', no_xcompact3d_hook_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_fluid_field_access_status', no_fluid_field_access_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_fluid_field_modification_status', no_fluid_field_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_pressure_modification_status', no_pressure_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_projection_modification_status', no_projection_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_poisson_modification_status', no_poisson_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_rk3_modification_status', no_rk3_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_channel_forcing_modification_status', &
                                 no_channel_forcing_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_production_ibm_forcing_status', no_production_ibm_forcing_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_feedback_application_status', no_feedback_application_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_twoway_force_status', no_twoway_force_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_no_structure_advance_status', no_structure_advance_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_rk_timing_audit_status', rk_timing_audit_status
    write(unit_id, '(A,1X,I0)') 'stage14_4_expected_nsubsteps', expected_nsubsteps
    write(unit_id, '(A,1X,I0)') 'stage14_4_recorded_rhs_addition_count', recorded_rhs_addition_count
    write(unit_id, '(A,1X,I0)') 'stage14_4_duplicate_injection_count', duplicate_injection_count
    write(unit_id, '(A,1X,I0)') 'stage14_4_skipped_substep_count', skipped_substep_count
    write(unit_id, '(A,1X,I0)') 'stage14_4_wrong_order_count', wrong_order_count
    write(unit_id, '(A,1X,I0)') 'stage14_4_invalid_event_count', invalid_event_count
    write(unit_id, '(A,1X,I0)') 'stage14_4_invalid_substep_count', invalid_substep_count
    write(unit_id, '(A,1X,ES24.16)') 'stage14_4_lambda_zero_max_abs_rhs_increment', &
                                      lambda_zero_max_abs_rhs_increment
    close(unit_id)
  end subroutine stage14_rk_timing_audit_write_diagnostics

  subroutine reset_sequence_state()
    event_count = 0
    first_event_order = 0
    event_order_counter = 0
    rk_substep_count_status = logical_to_int(expected_nsubsteps == max_substeps)
    stage13_candidate_before_increment_status = 0
    increment_before_rhs_addition_status = 0
    rhs_addition_before_projection_status = 0
    once_per_substep_status = 0
    no_duplicate_injection_status = 0
    no_skipped_substep_status = 0
    valid_event_id_status = 1
    valid_substep_index_status = 1
    recorded_rhs_addition_count = 0
    duplicate_injection_count = 0
    skipped_substep_count = 0
    wrong_order_count = 0
    invalid_event_count = 0
    invalid_substep_count = 0
  end subroutine reset_sequence_state

  logical function all_required_events_once(substep_index)
    integer, intent(in) :: substep_index
    all_required_events_once = event_count(substep_index,EVENT_STAGE13_CANDIDATE_READY) == 1 .and. &
                               event_count(substep_index,EVENT_STAGE14_INCREMENT_READY) == 1 .and. &
                               event_count(substep_index,EVENT_RHS_ADDITION_PLANNED) == 1 .and. &
                               event_count(substep_index,EVENT_PROJECTION_PLANNED) == 1
  end function all_required_events_once

  logical function ordered_before(substep_index, earlier_event, later_event)
    integer, intent(in) :: substep_index
    integer, intent(in) :: earlier_event
    integer, intent(in) :: later_event
    ordered_before = first_event_order(substep_index,earlier_event) > 0 .and. &
                     first_event_order(substep_index,later_event) > 0 .and. &
                     first_event_order(substep_index,earlier_event) < first_event_order(substep_index,later_event)
  end function ordered_before

  subroutine update_timing_finite_status()
    timing_diagnostics_finite_status = logical_to_int(finite_real(lambda_zero_max_abs_rhs_increment) .and. &
                                                      expected_nsubsteps == max_substeps .and. &
                                                      recorded_rhs_addition_count >= 0 .and. &
                                                      duplicate_injection_count >= 0 .and. &
                                                      skipped_substep_count >= 0 .and. &
                                                      wrong_order_count >= 0 .and. &
                                                      invalid_event_count >= 0 .and. &
                                                      invalid_substep_count >= 0)
  end subroutine update_timing_finite_status

  subroutine update_overall_status()
    rk_timing_audit_status = logical_to_int(initialized_status == 1 .and. &
                                            rk_substep_count_status == 1 .and. &
                                            stage13_candidate_before_increment_status == 1 .and. &
                                            increment_before_rhs_addition_status == 1 .and. &
                                            rhs_addition_before_projection_status == 1 .and. &
                                            once_per_substep_status == 1 .and. &
                                            no_duplicate_injection_status == 1 .and. &
                                            no_skipped_substep_status == 1 .and. &
                                            duplicate_detection_status == 1 .and. &
                                            skipped_substep_detection_status == 1 .and. &
                                            wrong_order_detection_status == 1 .and. &
                                            invalid_event_rejection_status == 1 .and. &
                                            invalid_substep_rejection_status == 1 .and. &
                                            valid_event_id_status == 1 .and. &
                                            valid_substep_index_status == 1 .and. &
                                            lambda_zero_increment_status == 1 .and. &
                                            timing_diagnostics_finite_status == 1 .and. &
                                            no_production_rhs_modification_status == 1 .and. &
                                            no_xcompact3d_hook_status == 1 .and. &
                                            no_fluid_field_access_status == 1 .and. &
                                            no_fluid_field_modification_status == 1 .and. &
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

end module fibre_stage14_rk_timing_audit
