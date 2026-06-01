program fibre_stage14_rhs_accumulator_check
  use decomp_2d_constants, only : mytype
  use fibre_stage14_config
  use fibre_stage14_rhs_accumulator
  implicit none

  integer, parameter :: nx = 4
  integer, parameter :: ny = 5
  integer, parameter :: nz = 6
  real(mytype), allocatable :: fx_cand(:,:,:)
  real(mytype), allocatable :: fy_cand(:,:,:)
  real(mytype), allocatable :: fz_cand(:,:,:)
  integer :: i
  integer :: j
  integer :: k
  integer :: allocated_status
  integer :: shape_status
  integer :: zero_initialization_status
  integer :: clear_status
  integer :: lambda_zero_status
  integer :: lambda_one_scaling_status
  integer :: lambda_fractional_scaling_status
  integer :: component_scaling_status
  integer :: finite_accumulator_status
  integer :: rhs_increment_norm_finite_status
  integer :: rhs_increment_valid_flag_status
  integer :: no_rhs_addition_status
  integer :: no_production_hook_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: no_pressure_modification_status
  integer :: no_projection_modification_status
  integer :: no_rk3_modification_status
  integer :: no_channel_forcing_modification_status
  integer :: no_production_ibm_forcing_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: rhs_accumulator_status
  integer :: requested_flag
  integer :: rhs_injection_enabled_flag
  integer :: injection_gain_finite_status
  logical :: invalid_shape_rejected

  call execute_command_line('mkdir -p stage14_outputs')

  call stage14_config_load()
  requested_flag = logical_to_int_local(stage14_requested())
  rhs_injection_enabled_flag = logical_to_int_local(stage14_rhs_injection_enabled())
  injection_gain_finite_status = logical_to_int_local(finite_real_local(stage14_get_injection_gain()))

  call stage14_rhs_accumulator_init(0, ny, nz)
  invalid_shape_rejected = .not. stage14_rhs_accumulator_is_allocated()
  call stage14_rhs_accumulator_finalize()
  if (.not. invalid_shape_rejected) then
    print *, 'STAGE 14.1 RHS ACCUMULATOR VERDICT: FAIL'
    print *, 'Reason: invalid_shape_was_allocated'
    stop 1
  end if

  allocate(fx_cand(nx,ny,nz), fy_cand(nx,ny,nz), fz_cand(nx,ny,nz))
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        fx_cand(i,j,k) = 0.25_mytype * real(i, mytype) + &
                         0.01_mytype * real(j, mytype) + &
                         0.001_mytype * real(k, mytype)
        fy_cand(i,j,k) = -0.15_mytype * real(i, mytype) + &
                         0.02_mytype * real(j, mytype) - &
                         0.002_mytype * real(k, mytype)
        fz_cand(i,j,k) = 0.05_mytype * real(i + j + k, mytype)
      end do
    end do
  end do

  call stage14_rhs_accumulator_init(nx, ny, nz)
  call stage14_rhs_accumulator_compute_from_force_density(fx_cand, fy_cand, fz_cand, 0.0_mytype)
  call stage14_rhs_accumulator_compute_from_force_density(fx_cand, fy_cand, fz_cand, 1.0_mytype)
  call stage14_rhs_accumulator_compute_from_force_density(fx_cand, fy_cand, fz_cand, 0.1_mytype)
  call stage14_rhs_accumulator_clear()

  call stage14_rhs_accumulator_get_status_values(allocated_status, &
                                                 shape_status, &
                                                 zero_initialization_status, &
                                                 clear_status, &
                                                 lambda_zero_status, &
                                                 lambda_one_scaling_status, &
                                                 lambda_fractional_scaling_status, &
                                                 component_scaling_status, &
                                                 finite_accumulator_status, &
                                                 rhs_increment_norm_finite_status, &
                                                 rhs_increment_valid_flag_status, &
                                                 no_rhs_addition_status, &
                                                 no_production_hook_status, &
                                                 no_fluid_field_access_status, &
                                                 no_fluid_field_modification_status, &
                                                 no_pressure_modification_status, &
                                                 no_projection_modification_status, &
                                                 no_rk3_modification_status, &
                                                 no_channel_forcing_modification_status, &
                                                 no_production_ibm_forcing_status, &
                                                 no_feedback_application_status, &
                                                 no_twoway_force_status, &
                                                 no_structure_advance_status, &
                                                 rhs_accumulator_status)

  call stage14_rhs_accumulator_write_diagnostics( &
       'stage14_outputs/fibre_stage14_1_rhs_accumulator.dat', &
       requested_flag, rhs_injection_enabled_flag, injection_gain_finite_status)

  if (rhs_accumulator_status == 1 .and. requested_flag == 1 .and. &
      rhs_injection_enabled_flag == 1 .and. injection_gain_finite_status == 1) then
    print *, 'STAGE 14.1 RHS ACCUMULATOR VERDICT: PASS'
  else
    print *, 'STAGE 14.1 RHS ACCUMULATOR VERDICT: FAIL'
    if (requested_flag /= 1) print *, 'Reason: stage14_1_requested_flag'
    if (rhs_injection_enabled_flag /= 1) print *, 'Reason: stage14_1_rhs_injection_enabled_flag'
    if (injection_gain_finite_status /= 1) print *, 'Reason: stage14_1_injection_gain_finite_status'
    if (allocated_status /= 1) print *, 'Reason: stage14_1_allocated_status'
    if (shape_status /= 1) print *, 'Reason: stage14_1_shape_status'
    if (zero_initialization_status /= 1) print *, 'Reason: stage14_1_zero_initialization_status'
    if (clear_status /= 1) print *, 'Reason: stage14_1_clear_status'
    if (lambda_zero_status /= 1) print *, 'Reason: stage14_1_lambda_zero_status'
    if (lambda_one_scaling_status /= 1) print *, 'Reason: stage14_1_lambda_one_scaling_status'
    if (lambda_fractional_scaling_status /= 1) print *, 'Reason: stage14_1_lambda_fractional_scaling_status'
    if (component_scaling_status /= 1) print *, 'Reason: stage14_1_component_scaling_status'
    if (finite_accumulator_status /= 1) print *, 'Reason: stage14_1_finite_accumulator_status'
    if (rhs_increment_norm_finite_status /= 1) print *, 'Reason: stage14_1_rhs_increment_norm_finite_status'
    if (rhs_increment_valid_flag_status /= 1) print *, 'Reason: stage14_1_rhs_increment_valid_flag_status'
    if (no_rhs_addition_status /= 1) print *, 'Reason: stage14_1_no_rhs_addition_status'
    if (no_production_hook_status /= 1) print *, 'Reason: stage14_1_no_production_hook_status'
    if (no_fluid_field_access_status /= 1) print *, 'Reason: stage14_1_no_fluid_field_access_status'
    if (no_fluid_field_modification_status /= 1) print *, 'Reason: stage14_1_no_fluid_field_modification_status'
    if (no_pressure_modification_status /= 1) print *, 'Reason: stage14_1_no_pressure_modification_status'
    if (no_projection_modification_status /= 1) print *, 'Reason: stage14_1_no_projection_modification_status'
    if (no_rk3_modification_status /= 1) print *, 'Reason: stage14_1_no_rk3_modification_status'
    if (no_channel_forcing_modification_status /= 1) print *, 'Reason: stage14_1_no_channel_forcing_modification_status'
    if (no_production_ibm_forcing_status /= 1) print *, 'Reason: stage14_1_no_production_ibm_forcing_status'
    if (no_feedback_application_status /= 1) print *, 'Reason: stage14_1_no_feedback_application_status'
    if (no_twoway_force_status /= 1) print *, 'Reason: stage14_1_no_twoway_force_status'
    if (no_structure_advance_status /= 1) print *, 'Reason: stage14_1_no_structure_advance_status'
    if (rhs_accumulator_status /= 1) print *, 'Reason: stage14_1_rhs_accumulator_status'
    stop 1
  end if

  call stage14_rhs_accumulator_finalize()
  deallocate(fx_cand, fy_cand, fz_cand)

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

end program fibre_stage14_rhs_accumulator_check
