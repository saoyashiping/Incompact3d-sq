program fibre_stage12_production_feedback_candidate_check
  use decomp_2d_constants, only : mytype
  use fibre_stage12_config, only : stage12_config_load
  use fibre_stage12_production_feedback_candidate, only : &
       stage12_production_feedback_candidate_init, &
       stage12_production_feedback_candidate_sample, &
       stage12_production_feedback_candidate_finalize, &
       stage12_production_feedback_candidate_get_status_values
  implicit none

  integer, parameter :: nx = 8
  integer, parameter :: ny = 9
  integer, parameter :: nz = 8
  real(mytype), parameter :: unchanged_tol = 0.0_mytype

  real(mytype) :: ux(nx, ny, nz)
  real(mytype) :: uy(nx, ny, nz)
  real(mytype) :: uz(nx, ny, nz)
  real(mytype) :: ux_before(nx, ny, nz)
  real(mytype) :: uy_before(nx, ny, nz)
  real(mytype) :: uz_before(nx, ny, nz)
  real(mytype) :: signature_before
  real(mytype) :: signature_after
  real(mytype) :: field_delta
  integer :: requested_flag
  integer :: readonly_mode_status
  integer :: hook_initialized_status
  integer :: hook_sample_called_status
  integer :: sampled_velocity_available_status
  integer :: force_candidate_computed_status
  integer :: force_candidate_finite_status
  integer :: force_norm_finite_status
  integer :: power_diagnostics_finite_status
  integer :: action_reaction_status
  integer :: pair_power_consistency_status
  integer :: field_modified_status
  integer :: rhs_modified_status
  integer :: no_eulerian_force_density_status
  integer :: no_rhs_injection_status
  integer :: no_ibm_spreading_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: production_feedback_candidate_status
  integer :: field_unchanged_status
  integer :: io_unit
  integer :: ios

  call execute_command_line('mkdir -p stage12_outputs', exitstat=ios)
  call stage12_config_load()
  call fill_velocity_fields(ux, uy, uz)
  ux_before(:, :, :) = ux(:, :, :)
  uy_before(:, :, :) = uy(:, :, :)
  uz_before(:, :, :) = uz(:, :, :)
  signature_before = field_signature(ux, uy, uz)

  call stage12_production_feedback_candidate_init()
  call stage12_production_feedback_candidate_sample(ux, uy, uz)
  call stage12_production_feedback_candidate_finalize()

  signature_after = field_signature(ux, uy, uz)
  field_delta = maxval(abs(ux(:, :, :) - ux_before(:, :, :)))
  field_delta = max(field_delta, maxval(abs(uy(:, :, :) - uy_before(:, :, :))))
  field_delta = max(field_delta, maxval(abs(uz(:, :, :) - uz_before(:, :, :))))
  field_delta = max(field_delta, abs(signature_after - signature_before))
  field_unchanged_status = merge(1, 0, field_delta <= unchanged_tol)

  call stage12_production_feedback_candidate_get_status_values(requested_flag, readonly_mode_status, &
       hook_initialized_status, hook_sample_called_status, sampled_velocity_available_status, &
       force_candidate_computed_status, force_candidate_finite_status, force_norm_finite_status, &
       power_diagnostics_finite_status, action_reaction_status, pair_power_consistency_status, &
       field_modified_status, rhs_modified_status, no_eulerian_force_density_status, no_rhs_injection_status, &
       no_ibm_spreading_status, no_feedback_application_status, no_twoway_force_status, &
       no_structure_advance_status, production_feedback_candidate_status)

  if (field_unchanged_status /= 1) production_feedback_candidate_status = 0
  if (field_modified_status /= 0) production_feedback_candidate_status = 0

  open(newunit=io_unit, file='stage12_outputs/fibre_stage12_6_production_feedback_candidate_check.dat', &
       status='replace', action='write', iostat=ios)
  if (ios /= 0) then
    print *, 'STAGE 12.6 PRODUCTION FEEDBACK CANDIDATE CHECK VERDICT: FAIL'
    print *, 'Reason: could not open Stage 12.6 standalone check diagnostics file.'
    stop 1
  end if

  call write_int(io_unit, 'stage12_6_check_requested_flag', requested_flag)
  call write_int(io_unit, 'stage12_6_check_readonly_mode_status', readonly_mode_status)
  call write_int(io_unit, 'stage12_6_check_hook_initialized_status', hook_initialized_status)
  call write_int(io_unit, 'stage12_6_check_hook_sample_called_status', hook_sample_called_status)
  call write_int(io_unit, 'stage12_6_check_sampled_velocity_available_status', sampled_velocity_available_status)
  call write_int(io_unit, 'stage12_6_check_force_candidate_computed_status', force_candidate_computed_status)
  call write_int(io_unit, 'stage12_6_check_force_candidate_finite_status', force_candidate_finite_status)
  call write_int(io_unit, 'stage12_6_check_force_norm_finite_status', force_norm_finite_status)
  call write_int(io_unit, 'stage12_6_check_power_diagnostics_finite_status', power_diagnostics_finite_status)
  call write_int(io_unit, 'stage12_6_check_action_reaction_status', action_reaction_status)
  call write_int(io_unit, 'stage12_6_check_pair_power_consistency_status', pair_power_consistency_status)
  call write_int(io_unit, 'stage12_6_check_field_unchanged_status', field_unchanged_status)
  call write_int(io_unit, 'stage12_6_check_rhs_modified_status', rhs_modified_status)
  call write_int(io_unit, 'stage12_6_check_no_eulerian_force_density_status', no_eulerian_force_density_status)
  call write_int(io_unit, 'stage12_6_check_no_rhs_injection_status', no_rhs_injection_status)
  call write_int(io_unit, 'stage12_6_check_no_ibm_spreading_status', no_ibm_spreading_status)
  call write_int(io_unit, 'stage12_6_check_no_feedback_application_status', no_feedback_application_status)
  call write_int(io_unit, 'stage12_6_check_no_twoway_force_status', no_twoway_force_status)
  call write_int(io_unit, 'stage12_6_check_no_structure_advance_status', no_structure_advance_status)
  call write_int(io_unit, 'stage12_6_check_production_feedback_candidate_status', &
                 production_feedback_candidate_status)
  call write_real(io_unit, 'stage12_6_check_field_delta', field_delta)
  call write_real(io_unit, 'stage12_6_check_signature_before', signature_before)
  call write_real(io_unit, 'stage12_6_check_signature_after', signature_after)
  close(io_unit)

  if (requested_flag == 1 .and. readonly_mode_status == 1 .and. hook_initialized_status == 1 .and. &
      hook_sample_called_status == 1 .and. sampled_velocity_available_status == 1 .and. &
      force_candidate_computed_status == 1 .and. force_candidate_finite_status == 1 .and. &
      force_norm_finite_status == 1 .and. power_diagnostics_finite_status == 1 .and. &
      action_reaction_status == 1 .and. pair_power_consistency_status == 1 .and. &
      field_unchanged_status == 1 .and. rhs_modified_status == 0 .and. &
      no_eulerian_force_density_status == 1 .and. no_rhs_injection_status == 1 .and. &
      no_ibm_spreading_status == 1 .and. no_feedback_application_status == 1 .and. &
      no_twoway_force_status == 1 .and. no_structure_advance_status == 1 .and. &
      production_feedback_candidate_status == 1) then
    print *, 'STAGE 12.6 PRODUCTION FEEDBACK CANDIDATE CHECK VERDICT: PASS'
  else
    print *, 'STAGE 12.6 PRODUCTION FEEDBACK CANDIDATE CHECK VERDICT: FAIL'
    if (requested_flag /= 1) print *, 'Reason: Stage 12 feedback candidate was not requested.'
    if (readonly_mode_status /= 1) print *, 'Reason: Stage 12 readonly mode was not enforced.'
    if (hook_initialized_status /= 1) print *, 'Reason: production candidate hook was not initialized.'
    if (hook_sample_called_status /= 1) print *, 'Reason: production candidate sample was not called.'
    if (sampled_velocity_available_status /= 1) print *, 'Reason: sampled velocity was unavailable or non-finite.'
    if (force_candidate_computed_status /= 1) print *, 'Reason: force candidate was not computed.'
    if (force_candidate_finite_status /= 1) print *, 'Reason: force candidate was non-finite.'
    if (force_norm_finite_status /= 1) print *, 'Reason: force norm was non-finite.'
    if (power_diagnostics_finite_status /= 1) print *, 'Reason: power diagnostics were non-finite.'
    if (action_reaction_status /= 1) print *, 'Reason: action-reaction diagnostic failed.'
    if (pair_power_consistency_status /= 1) print *, 'Reason: pair-power consistency failed.'
    if (field_unchanged_status /= 1) print *, 'Reason: synthetic production velocity arrays changed.'
    if (rhs_modified_status /= 0) print *, 'Reason: RHS modification status was nonzero.'
    if (production_feedback_candidate_status /= 1) print *, 'Reason: production feedback candidate status failed.'
    stop 1
  end if

contains

  subroutine fill_velocity_fields(ux_field, uy_field, uz_field)
    real(mytype), intent(out) :: ux_field(:, :, :)
    real(mytype), intent(out) :: uy_field(:, :, :)
    real(mytype), intent(out) :: uz_field(:, :, :)
    integer :: i
    integer :: j
    integer :: k

    do k = 1, size(ux_field, 3)
      do j = 1, size(ux_field, 2)
        do i = 1, size(ux_field, 1)
          ux_field(i, j, k) = 0.10_mytype + 0.01_mytype * real(i, mytype) + 0.001_mytype * real(k, mytype)
          uy_field(i, j, k) = -0.05_mytype + 0.02_mytype * real(j, mytype)
          uz_field(i, j, k) = 0.03_mytype + 0.004_mytype * real(i + j + k, mytype)
        end do
      end do
    end do
  end subroutine fill_velocity_fields

  real(mytype) function field_signature(ux_field, uy_field, uz_field)
    real(mytype), intent(in) :: ux_field(:, :, :)
    real(mytype), intent(in) :: uy_field(:, :, :)
    real(mytype), intent(in) :: uz_field(:, :, :)

    field_signature = sum(ux_field(:, :, :)) + 2.0_mytype * sum(uy_field(:, :, :)) + &
                      3.0_mytype * sum(uz_field(:, :, :))
  end function field_signature

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
    write(io_unit_in, '(A,1X,ES24.16E3)') trim(key), value
  end subroutine write_real

end program fibre_stage12_production_feedback_candidate_check
