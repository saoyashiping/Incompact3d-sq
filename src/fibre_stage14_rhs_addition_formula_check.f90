program fibre_stage14_rhs_addition_formula_check
  use decomp_2d_constants, only : mytype
  use fibre_stage14_config
  use fibre_stage14_rhs_accumulator
  use fibre_stage14_rhs_addition_formula
  implicit none

  integer, parameter :: nx = 4
  integer, parameter :: ny = 5
  integer, parameter :: nz = 6
  real(mytype), allocatable :: fx_cand(:,:,:)
  real(mytype), allocatable :: fy_cand(:,:,:)
  real(mytype), allocatable :: fz_cand(:,:,:)
  real(mytype), allocatable :: rhs_old_x(:,:,:)
  real(mytype), allocatable :: rhs_old_y(:,:,:)
  real(mytype), allocatable :: rhs_old_z(:,:,:)
  real(mytype), allocatable :: rhs_old_x_ref(:,:,:)
  real(mytype), allocatable :: rhs_old_y_ref(:,:,:)
  real(mytype), allocatable :: rhs_old_z_ref(:,:,:)
  real(mytype), allocatable :: rhs_inc_x(:,:,:)
  real(mytype), allocatable :: rhs_inc_y(:,:,:)
  real(mytype), allocatable :: rhs_inc_z(:,:,:)
  real(mytype), allocatable :: rhs_new_x(:,:,:)
  real(mytype), allocatable :: rhs_new_y(:,:,:)
  real(mytype), allocatable :: rhs_new_z(:,:,:)
  real(mytype), allocatable :: rhs_bad(:,:,:)
  integer :: i
  integer :: j
  integer :: k
  integer :: initialized_status
  integer :: shape_status
  integer :: lambda_zero_invariance_status
  integer :: lambda_one_addition_status
  integer :: lambda_fractional_addition_status
  integer :: component_addition_status
  integer :: additive_not_overwrite_status
  integer :: rhs_old_preserved_status
  integer :: finite_rhs_old_status
  integer :: finite_rhs_increment_status
  integer :: finite_rhs_new_status
  integer :: no_production_rhs_modification_status
  integer :: no_xcompact3d_hook_status
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
  integer :: rhs_addition_formula_status
  integer :: requested_flag
  integer :: rhs_injection_enabled_flag
  integer :: injection_gain_finite_status
  real(mytype) :: preservation_error
  real(mytype) :: tolerance

  tolerance = 1.0e-12_mytype
  call execute_command_line('mkdir -p stage14_outputs')

  call stage14_config_load()
  requested_flag = logical_to_int_local(stage14_requested())
  rhs_injection_enabled_flag = logical_to_int_local(stage14_rhs_injection_enabled())
  injection_gain_finite_status = logical_to_int_local(finite_real_local(stage14_get_injection_gain()))

  allocate(fx_cand(nx,ny,nz), fy_cand(nx,ny,nz), fz_cand(nx,ny,nz))
  allocate(rhs_old_x(nx,ny,nz), rhs_old_y(nx,ny,nz), rhs_old_z(nx,ny,nz))
  allocate(rhs_old_x_ref(nx,ny,nz), rhs_old_y_ref(nx,ny,nz), rhs_old_z_ref(nx,ny,nz))
  allocate(rhs_inc_x(nx,ny,nz), rhs_inc_y(nx,ny,nz), rhs_inc_z(nx,ny,nz))
  allocate(rhs_new_x(nx,ny,nz), rhs_new_y(nx,ny,nz), rhs_new_z(nx,ny,nz))
  allocate(rhs_bad(nx+1,ny,nz))

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
        rhs_old_x(i,j,k) = -0.30_mytype + 0.07_mytype * real(i, mytype) - &
                           0.013_mytype * real(j, mytype) + 0.002_mytype * real(k, mytype)
        rhs_old_y(i,j,k) = 0.41_mytype - 0.021_mytype * real(i, mytype) + &
                           0.03_mytype * real(j, mytype) - 0.004_mytype * real(k, mytype)
        rhs_old_z(i,j,k) = -0.18_mytype + 0.011_mytype * real(i*j, mytype) - &
                           0.006_mytype * real(k, mytype)
      end do
    end do
  end do
  rhs_old_x_ref = rhs_old_x
  rhs_old_y_ref = rhs_old_y
  rhs_old_z_ref = rhs_old_z
  rhs_bad = 0.0_mytype

  call stage14_rhs_accumulator_init(nx, ny, nz)
  call stage14_rhs_addition_formula_init()

  call stage14_rhs_addition_formula_apply_controlled(rhs_old_x, rhs_old_y, rhs_old_z, &
                                                     rhs_inc_x, rhs_inc_y, rhs_inc_z, &
                                                     rhs_bad, rhs_new_y, rhs_new_z)
  call stage14_rhs_addition_formula_clear()

  call run_case(0.0_mytype)
  call run_case(1.0_mytype)
  call run_case(0.1_mytype)

  preservation_error = max(maxval(abs(rhs_old_x - rhs_old_x_ref)), &
                           max(maxval(abs(rhs_old_y - rhs_old_y_ref)), &
                               maxval(abs(rhs_old_z - rhs_old_z_ref))))

  call stage14_rhs_addition_formula_get_status_values(initialized_status, &
                                                      shape_status, &
                                                      lambda_zero_invariance_status, &
                                                      lambda_one_addition_status, &
                                                      lambda_fractional_addition_status, &
                                                      component_addition_status, &
                                                      additive_not_overwrite_status, &
                                                      rhs_old_preserved_status, &
                                                      finite_rhs_old_status, &
                                                      finite_rhs_increment_status, &
                                                      finite_rhs_new_status, &
                                                      no_production_rhs_modification_status, &
                                                      no_xcompact3d_hook_status, &
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
                                                      rhs_addition_formula_status)

  call stage14_rhs_addition_formula_write_diagnostics( &
       'stage14_outputs/fibre_stage14_2_rhs_addition_formula.dat', &
       requested_flag, rhs_injection_enabled_flag, injection_gain_finite_status, nx, ny, nz)

  if (rhs_addition_formula_status == 1 .and. requested_flag == 1 .and. &
      rhs_injection_enabled_flag == 1 .and. injection_gain_finite_status == 1 .and. &
      preservation_error <= tolerance) then
    print *, 'STAGE 14.2 RHS ADDITION FORMULA VERDICT: PASS'
  else
    print *, 'STAGE 14.2 RHS ADDITION FORMULA VERDICT: FAIL'
    if (requested_flag /= 1) print *, 'Reason: stage14_2_requested_flag'
    if (rhs_injection_enabled_flag /= 1) print *, 'Reason: stage14_2_rhs_injection_enabled_flag'
    if (injection_gain_finite_status /= 1) print *, 'Reason: stage14_2_injection_gain_finite_status'
    if (initialized_status /= 1) print *, 'Reason: stage14_2_initialized_status'
    if (shape_status /= 1) print *, 'Reason: stage14_2_shape_status'
    if (lambda_zero_invariance_status /= 1) print *, 'Reason: stage14_2_lambda_zero_invariance_status'
    if (lambda_one_addition_status /= 1) print *, 'Reason: stage14_2_lambda_one_addition_status'
    if (lambda_fractional_addition_status /= 1) print *, 'Reason: stage14_2_lambda_fractional_addition_status'
    if (component_addition_status /= 1) print *, 'Reason: stage14_2_component_addition_status'
    if (additive_not_overwrite_status /= 1) print *, 'Reason: stage14_2_additive_not_overwrite_status'
    if (rhs_old_preserved_status /= 1) print *, 'Reason: stage14_2_rhs_old_preserved_status'
    if (finite_rhs_old_status /= 1) print *, 'Reason: stage14_2_finite_rhs_old_status'
    if (finite_rhs_increment_status /= 1) print *, 'Reason: stage14_2_finite_rhs_increment_status'
    if (finite_rhs_new_status /= 1) print *, 'Reason: stage14_2_finite_rhs_new_status'
    if (no_production_rhs_modification_status /= 1) print *, 'Reason: stage14_2_no_production_rhs_modification_status'
    if (no_xcompact3d_hook_status /= 1) print *, 'Reason: stage14_2_no_xcompact3d_hook_status'
    if (no_fluid_field_access_status /= 1) print *, 'Reason: stage14_2_no_fluid_field_access_status'
    if (no_fluid_field_modification_status /= 1) print *, 'Reason: stage14_2_no_fluid_field_modification_status'
    if (no_pressure_modification_status /= 1) print *, 'Reason: stage14_2_no_pressure_modification_status'
    if (no_projection_modification_status /= 1) print *, 'Reason: stage14_2_no_projection_modification_status'
    if (no_rk3_modification_status /= 1) print *, 'Reason: stage14_2_no_rk3_modification_status'
    if (no_channel_forcing_modification_status /= 1) print *, 'Reason: stage14_2_no_channel_forcing_modification_status'
    if (no_production_ibm_forcing_status /= 1) print *, 'Reason: stage14_2_no_production_ibm_forcing_status'
    if (no_feedback_application_status /= 1) print *, 'Reason: stage14_2_no_feedback_application_status'
    if (no_twoway_force_status /= 1) print *, 'Reason: stage14_2_no_twoway_force_status'
    if (no_structure_advance_status /= 1) print *, 'Reason: stage14_2_no_structure_advance_status'
    if (rhs_addition_formula_status /= 1) print *, 'Reason: stage14_2_rhs_addition_formula_status'
    if (preservation_error > tolerance) print *, 'Reason: stage14_2_rhs_old_preservation_error'
    stop 1
  end if

  call stage14_rhs_accumulator_finalize()
  call stage14_rhs_addition_formula_finalize()
  deallocate(fx_cand, fy_cand, fz_cand)
  deallocate(rhs_old_x, rhs_old_y, rhs_old_z)
  deallocate(rhs_old_x_ref, rhs_old_y_ref, rhs_old_z_ref)
  deallocate(rhs_inc_x, rhs_inc_y, rhs_inc_z)
  deallocate(rhs_new_x, rhs_new_y, rhs_new_z, rhs_bad)

contains

  subroutine run_case(lambda_value)
    real(mytype), intent(in) :: lambda_value

    call stage14_rhs_accumulator_compute_from_force_density(fx_cand, fy_cand, fz_cand, lambda_value)
    rhs_inc_x = lambda_value * fx_cand
    rhs_inc_y = lambda_value * fy_cand
    rhs_inc_z = lambda_value * fz_cand
    call stage14_rhs_addition_formula_apply_controlled(rhs_old_x, rhs_old_y, rhs_old_z, &
                                                       rhs_inc_x, rhs_inc_y, rhs_inc_z, &
                                                       rhs_new_x, rhs_new_y, rhs_new_z)
  end subroutine run_case

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

end program fibre_stage14_rhs_addition_formula_check
