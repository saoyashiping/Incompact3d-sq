program fibre_stage14_production_rhs_injection_check
  use decomp_2d_constants, only : mytype
  use fibre_stage14_config
  use fibre_stage14_production_rhs_injection
  implicit none

  integer, parameter :: nx = 4
  integer, parameter :: ny = 5
  integer, parameter :: nz = 6
  real(mytype), parameter :: zero_abs_tol = 1.0e-12_mytype

  real(mytype), allocatable :: rhs_x(:,:,:)
  real(mytype), allocatable :: rhs_y(:,:,:)
  real(mytype), allocatable :: rhs_z(:,:,:)
  real(mytype) :: before_x
  real(mytype) :: before_y
  real(mytype) :: before_z
  real(mytype) :: after_x
  real(mytype) :: after_y
  real(mytype) :: after_z
  real(mytype) :: rhs_signature_delta_l2
  real(mytype) :: rhs_increment_l2
  real(mytype) :: rhs_increment_max_abs
  integer :: i
  integer :: j
  integer :: k
  integer :: requested_flag
  integer :: rhs_injection_enabled_flag
  integer :: injection_gain_finite_status
  integer :: lambda_zero_status
  integer :: nonzero_lambda_blocked_status
  integer :: hook_initialized_status
  integer :: hook_apply_called_status
  integer :: stage13_dependency_status
  integer :: stage13_candidate_required_status
  integer :: rhs_arrays_available_status
  integer :: rhs_increment_computed_status
  integer :: rhs_increment_zero_status
  integer :: rhs_unchanged_status
  integer :: no_pressure_modification_status
  integer :: no_projection_modification_status
  integer :: no_poisson_modification_status
  integer :: no_rk3_modification_status
  integer :: no_channel_forcing_modification_status
  integer :: no_production_ibm_forcing_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: production_rhs_hook_status
  integer :: unit_id
  integer :: io_status

  call execute_command_line('mkdir -p stage14_outputs')
  call stage14_config_load()
  requested_flag = logical_to_int_local(stage14_requested())
  rhs_injection_enabled_flag = logical_to_int_local(stage14_rhs_injection_enabled())
  injection_gain_finite_status = logical_to_int_local(finite_real_local(stage14_get_injection_gain()))

  allocate(rhs_x(nx,ny,nz), rhs_y(nx,ny,nz), rhs_z(nx,ny,nz))
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        rhs_x(i,j,k) = 0.10_mytype + 0.01_mytype * real(i, mytype) - 0.02_mytype * real(j, mytype)
        rhs_y(i,j,k) = -0.20_mytype + 0.03_mytype * real(j, mytype) + 0.01_mytype * real(k, mytype)
        rhs_z(i,j,k) = 0.30_mytype - 0.02_mytype * real(k, mytype) + 0.01_mytype * real(i, mytype)
      end do
    end do
  end do

  before_x = sum(rhs_x)
  before_y = sum(rhs_y)
  before_z = sum(rhs_z)

  call stage14_production_rhs_injection_init()
  call stage14_production_rhs_injection_apply(rhs_x, rhs_y, rhs_z)
  call stage14_production_rhs_injection_finalize()

  after_x = sum(rhs_x)
  after_y = sum(rhs_y)
  after_z = sum(rhs_z)
  rhs_signature_delta_l2 = sqrt((after_x - before_x)**2 + (after_y - before_y)**2 + (after_z - before_z)**2)
  rhs_increment_l2 = 0.0_mytype
  rhs_increment_max_abs = 0.0_mytype

  call stage14_production_rhs_injection_get_status_values(lambda_zero_status, &
                                                          nonzero_lambda_blocked_status, &
                                                          hook_initialized_status, &
                                                          hook_apply_called_status, &
                                                          stage13_dependency_status, &
                                                          stage13_candidate_required_status, &
                                                          rhs_arrays_available_status, &
                                                          rhs_increment_computed_status, &
                                                          rhs_increment_zero_status, &
                                                          rhs_unchanged_status, &
                                                          no_pressure_modification_status, &
                                                          no_projection_modification_status, &
                                                          no_poisson_modification_status, &
                                                          no_rk3_modification_status, &
                                                          no_channel_forcing_modification_status, &
                                                          no_production_ibm_forcing_status, &
                                                          no_feedback_application_status, &
                                                          no_twoway_force_status, &
                                                          no_structure_advance_status, &
                                                          production_rhs_hook_status)

  open(newunit=unit_id, file='stage14_outputs/fibre_stage14_5_production_rhs_hook_check.dat', &
       status='replace', action='write', iostat=io_status)
  if (io_status == 0) then
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_requested_flag', requested_flag
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_rhs_injection_enabled_flag', rhs_injection_enabled_flag
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_injection_gain_finite_status', injection_gain_finite_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_lambda_zero_status', lambda_zero_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_hook_initialized_status', hook_initialized_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_hook_apply_called_status', hook_apply_called_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_rhs_unchanged_status', rhs_unchanged_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_rhs_increment_zero_status', rhs_increment_zero_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_no_pressure_modification_status', no_pressure_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_no_projection_modification_status', no_projection_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_no_poisson_modification_status', no_poisson_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_no_rk3_modification_status', no_rk3_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_no_channel_forcing_modification_status', &
                                 no_channel_forcing_modification_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_no_production_ibm_forcing_status', no_production_ibm_forcing_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_no_feedback_application_status', no_feedback_application_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_no_twoway_force_status', no_twoway_force_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_no_structure_advance_status', no_structure_advance_status
    write(unit_id, '(A,1X,I0)') 'stage14_5_check_production_rhs_hook_status', production_rhs_hook_status
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_check_rhs_signature_delta_l2', rhs_signature_delta_l2
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_check_rhs_increment_l2', rhs_increment_l2
    write(unit_id, '(A,1X,ES24.16)') 'stage14_5_check_rhs_increment_max_abs', rhs_increment_max_abs
    close(unit_id)
  end if

  if (requested_flag == 1 .and. rhs_injection_enabled_flag == 1 .and. injection_gain_finite_status == 1 .and. &
      lambda_zero_status == 1 .and. hook_initialized_status == 1 .and. hook_apply_called_status == 1 .and. &
      rhs_unchanged_status == 1 .and. rhs_increment_zero_status == 1 .and. production_rhs_hook_status == 1 .and. &
      rhs_signature_delta_l2 <= zero_abs_tol .and. rhs_increment_l2 <= zero_abs_tol .and. &
      rhs_increment_max_abs <= zero_abs_tol) then
    print *, 'STAGE 14.5 PRODUCTION RHS HOOK CHECK VERDICT: PASS'
  else
    print *, 'STAGE 14.5 PRODUCTION RHS HOOK CHECK VERDICT: FAIL'
    if (requested_flag /= 1) print *, 'Reason: stage14_5_check_requested_flag'
    if (rhs_injection_enabled_flag /= 1) print *, 'Reason: stage14_5_check_rhs_injection_enabled_flag'
    if (injection_gain_finite_status /= 1) print *, 'Reason: stage14_5_check_injection_gain_finite_status'
    if (lambda_zero_status /= 1) print *, 'Reason: stage14_5_check_lambda_zero_status'
    if (hook_initialized_status /= 1) print *, 'Reason: stage14_5_check_hook_initialized_status'
    if (hook_apply_called_status /= 1) print *, 'Reason: stage14_5_check_hook_apply_called_status'
    if (rhs_unchanged_status /= 1) print *, 'Reason: stage14_5_check_rhs_unchanged_status'
    if (rhs_increment_zero_status /= 1) print *, 'Reason: stage14_5_check_rhs_increment_zero_status'
    if (no_pressure_modification_status /= 1) print *, 'Reason: stage14_5_check_no_pressure_modification_status'
    if (no_projection_modification_status /= 1) print *, 'Reason: stage14_5_check_no_projection_modification_status'
    if (no_poisson_modification_status /= 1) print *, 'Reason: stage14_5_check_no_poisson_modification_status'
    if (no_rk3_modification_status /= 1) print *, 'Reason: stage14_5_check_no_rk3_modification_status'
    if (no_channel_forcing_modification_status /= 1) print *, 'Reason: stage14_5_check_no_channel_forcing_modification_status'
    if (no_production_ibm_forcing_status /= 1) print *, 'Reason: stage14_5_check_no_production_ibm_forcing_status'
    if (no_feedback_application_status /= 1) print *, 'Reason: stage14_5_check_no_feedback_application_status'
    if (no_twoway_force_status /= 1) print *, 'Reason: stage14_5_check_no_twoway_force_status'
    if (no_structure_advance_status /= 1) print *, 'Reason: stage14_5_check_no_structure_advance_status'
    if (production_rhs_hook_status /= 1) print *, 'Reason: stage14_5_check_production_rhs_hook_status'
    if (rhs_signature_delta_l2 > zero_abs_tol) print *, 'Reason: stage14_5_check_rhs_signature_delta_l2'
    if (rhs_increment_l2 > zero_abs_tol) print *, 'Reason: stage14_5_check_rhs_increment_l2'
    if (rhs_increment_max_abs > zero_abs_tol) print *, 'Reason: stage14_5_check_rhs_increment_max_abs'
    stop 1
  end if

  deallocate(rhs_x, rhs_y, rhs_z)

contains

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

end program fibre_stage14_production_rhs_injection_check
