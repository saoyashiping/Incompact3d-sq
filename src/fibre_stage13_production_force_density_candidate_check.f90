program fibre_stage13_production_force_density_candidate_check
  use decomp_2d_constants, only : mytype
  use fibre_stage13_config, only : stage13_config_load
  use fibre_stage13_production_force_density_candidate, only : &
       stage13_production_force_density_candidate_init, &
       stage13_production_force_density_candidate_sample, &
       stage13_production_force_density_candidate_finalize, &
       stage13_production_force_density_candidate_get_status_values
  implicit none

  integer, parameter :: nx = 8
  integer, parameter :: ny = 9
  integer, parameter :: nz = 8
  real(mytype), parameter :: signature_abs_tol = 1.0e-12_mytype

  real(mytype) :: ux(nx,ny,nz)
  real(mytype) :: uy(nx,ny,nz)
  real(mytype) :: uz(nx,ny,nz)
  real(mytype) :: signature_before(3)
  real(mytype) :: signature_after(3)
  real(mytype) :: signature_error
  integer :: i
  integer :: j
  integer :: k
  integer :: requested_flag
  integer :: readonly_mode_status
  integer :: spreading_readonly_status
  integer :: hook_initialized_status
  integer :: hook_sample_called_status
  integer :: sampled_velocity_available_status
  integer :: force_density_candidate_computed_status
  integer :: force_density_candidate_finite_status
  integer :: force_density_norm_finite_status
  integer :: integrated_force_finite_status
  integer :: integrated_force_conservation_status
  integer :: spreading_input_sign_status
  integer :: wrong_sign_rejection_status
  integer :: field_modified_status
  integer :: rhs_modified_status
  integer :: no_rhs_injection_status
  integer :: no_production_ibm_forcing_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: production_force_density_candidate_status
  integer :: field_unchanged_status
  integer :: verdict_status
  integer :: io_unit
  integer :: ios

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        ux(i,j,k) = 0.10_mytype * real(i, mytype) + 0.01_mytype * real(j, mytype) + &
                   0.001_mytype * real(k, mytype)
        uy(i,j,k) = -0.05_mytype * real(i, mytype) + 0.02_mytype * real(j, mytype) - &
                   0.002_mytype * real(k, mytype)
        uz(i,j,k) = 0.03_mytype * real(i, mytype) - 0.01_mytype * real(j, mytype) + &
                   0.004_mytype * real(k, mytype)
      end do
    end do
  end do

  call compute_signature(ux, uy, uz, signature_before)
  call execute_command_line('mkdir -p stage13_outputs')
  call stage13_config_load()
  call stage13_production_force_density_candidate_init()
  call stage13_production_force_density_candidate_sample(ux, uy, uz)
  call stage13_production_force_density_candidate_finalize()
  call compute_signature(ux, uy, uz, signature_after)
  signature_error = maxval(abs(signature_after - signature_before))
  field_unchanged_status = merge(1, 0, signature_error <= signature_abs_tol)

  call stage13_production_force_density_candidate_get_status_values(requested_flag, readonly_mode_status, &
       spreading_readonly_status, hook_initialized_status, hook_sample_called_status, sampled_velocity_available_status, &
       force_density_candidate_computed_status, force_density_candidate_finite_status, &
       force_density_norm_finite_status, integrated_force_finite_status, integrated_force_conservation_status, &
       spreading_input_sign_status, wrong_sign_rejection_status, field_modified_status, rhs_modified_status, &
       no_rhs_injection_status, no_production_ibm_forcing_status, no_feedback_application_status, &
       no_twoway_force_status, no_structure_advance_status, production_force_density_candidate_status)

  verdict_status = 1
  if (requested_flag /= 1) verdict_status = 0
  if (readonly_mode_status /= 1) verdict_status = 0
  if (spreading_readonly_status /= 1) verdict_status = 0
  if (hook_initialized_status /= 1) verdict_status = 0
  if (hook_sample_called_status /= 1) verdict_status = 0
  if (sampled_velocity_available_status /= 1) verdict_status = 0
  if (force_density_candidate_computed_status /= 1) verdict_status = 0
  if (force_density_candidate_finite_status /= 1) verdict_status = 0
  if (force_density_norm_finite_status /= 1) verdict_status = 0
  if (integrated_force_finite_status /= 1) verdict_status = 0
  if (integrated_force_conservation_status /= 1) verdict_status = 0
  if (spreading_input_sign_status /= 1) verdict_status = 0
  if (wrong_sign_rejection_status /= 1) verdict_status = 0
  if (field_unchanged_status /= 1) verdict_status = 0
  if (field_modified_status /= 0) verdict_status = 0
  if (rhs_modified_status /= 0) verdict_status = 0
  if (no_rhs_injection_status /= 1) verdict_status = 0
  if (no_production_ibm_forcing_status /= 1) verdict_status = 0
  if (no_feedback_application_status /= 1) verdict_status = 0
  if (no_twoway_force_status /= 1) verdict_status = 0
  if (no_structure_advance_status /= 1) verdict_status = 0
  if (production_force_density_candidate_status /= 1) verdict_status = 0

  open(newunit=io_unit, file='stage13_outputs/fibre_stage13_6_production_force_density_candidate_check.dat', &
       status='replace', action='write', iostat=ios)
  if (ios == 0) then
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_requested_flag', requested_flag
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_readonly_mode_status', readonly_mode_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_spreading_readonly_status', spreading_readonly_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_hook_initialized_status', hook_initialized_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_hook_sample_called_status', hook_sample_called_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_sampled_velocity_available_status', &
         sampled_velocity_available_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_force_density_candidate_computed_status', &
         force_density_candidate_computed_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_force_density_candidate_finite_status', &
         force_density_candidate_finite_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_force_density_norm_finite_status', &
         force_density_norm_finite_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_integrated_force_finite_status', integrated_force_finite_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_integrated_force_conservation_status', &
         integrated_force_conservation_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_spreading_input_sign_status', spreading_input_sign_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_wrong_sign_rejection_status', wrong_sign_rejection_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_field_unchanged_status', field_unchanged_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_rhs_modified_status', rhs_modified_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_no_rhs_injection_status', no_rhs_injection_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_no_production_ibm_forcing_status', &
         no_production_ibm_forcing_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_no_feedback_application_status', no_feedback_application_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_no_twoway_force_status', no_twoway_force_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_no_structure_advance_status', no_structure_advance_status
    write(io_unit,'(A,1X,I0)') 'stage13_6_check_production_force_density_candidate_status', &
         production_force_density_candidate_status
    write(io_unit,'(A,1X,ES24.16E3)') 'stage13_6_check_signature_error', signature_error
    close(io_unit)
  end if

  if (verdict_status == 1) then
    write(*,'(A)') 'STAGE 13.6 PRODUCTION FORCE DENSITY CANDIDATE CHECK VERDICT: PASS'
  else
    write(*,'(A)') 'STAGE 13.6 PRODUCTION FORCE DENSITY CANDIDATE CHECK VERDICT: FAIL'
    call print_failures()
    stop 1
  end if

contains

  subroutine compute_signature(ux_in, uy_in, uz_in, signature)
    real(mytype), intent(in) :: ux_in(:,:,:)
    real(mytype), intent(in) :: uy_in(:,:,:)
    real(mytype), intent(in) :: uz_in(:,:,:)
    real(mytype), intent(out) :: signature(3)

    signature(1) = sum(ux_in)
    signature(2) = sum(uy_in)
    signature(3) = sum(uz_in)
  end subroutine compute_signature

  subroutine print_failures()
    if (requested_flag /= 1) write(*,'(A)') 'FAIL: requested flag is not enabled'
    if (readonly_mode_status /= 1) write(*,'(A)') 'FAIL: readonly mode is not enabled'
    if (spreading_readonly_status /= 1) write(*,'(A)') 'FAIL: spreading readonly mode is not enabled'
    if (hook_initialized_status /= 1) write(*,'(A)') 'FAIL: hook was not initialized'
    if (hook_sample_called_status /= 1) write(*,'(A)') 'FAIL: sample hook was not called'
    if (sampled_velocity_available_status /= 1) write(*,'(A)') 'FAIL: sampled velocity was unavailable or nonfinite'
    if (force_density_candidate_computed_status /= 1) write(*,'(A)') 'FAIL: force-density candidate was not computed'
    if (force_density_candidate_finite_status /= 1) write(*,'(A)') 'FAIL: force-density candidate was nonfinite'
    if (force_density_norm_finite_status /= 1) write(*,'(A)') 'FAIL: force-density norm was nonfinite'
    if (integrated_force_finite_status /= 1) write(*,'(A)') 'FAIL: integrated force was nonfinite'
    if (integrated_force_conservation_status /= 1) write(*,'(A)') 'FAIL: integrated force was not conserved'
    if (spreading_input_sign_status /= 1) write(*,'(A)') 'FAIL: spreading input sign was not fibre-to-fluid'
    if (wrong_sign_rejection_status /= 1) write(*,'(A)') 'FAIL: wrong-sign force was not rejected'
    if (field_unchanged_status /= 1) write(*,'(A)') 'FAIL: input velocity field signature changed'
    if (field_modified_status /= 0) write(*,'(A)') 'FAIL: module reported a production field modification'
    if (rhs_modified_status /= 0) write(*,'(A)') 'FAIL: module reported RHS modification'
    if (no_rhs_injection_status /= 1) write(*,'(A)') 'FAIL: RHS injection guard failed'
    if (no_production_ibm_forcing_status /= 1) write(*,'(A)') 'FAIL: production IBM forcing guard failed'
    if (no_feedback_application_status /= 1) write(*,'(A)') 'FAIL: feedback-application guard failed'
    if (no_twoway_force_status /= 1) write(*,'(A)') 'FAIL: two-way-force guard failed'
    if (no_structure_advance_status /= 1) write(*,'(A)') 'FAIL: structure-advance guard failed'
    if (production_force_density_candidate_status /= 1) write(*,'(A)') 'FAIL: aggregate Stage 13.6 status failed'
  end subroutine print_failures

end program fibre_stage13_production_force_density_candidate_check
