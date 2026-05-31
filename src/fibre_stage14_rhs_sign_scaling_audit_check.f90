program fibre_stage14_rhs_sign_scaling_audit_check
  use decomp_2d_constants, only : mytype
  use fibre_stage14_config
  use fibre_stage14_rhs_accumulator
  use fibre_stage14_rhs_addition_formula
  use fibre_stage14_rhs_sign_scaling_audit
  implicit none

  integer, parameter :: nx = 4
  integer, parameter :: ny = 5
  integer, parameter :: nz = 6
  real(mytype), parameter :: local_abs_tol = 1.0e-12_mytype
  real(mytype), allocatable :: volume_eul(:,:,:)
  real(mytype), allocatable :: fx_cand(:,:,:)
  real(mytype), allocatable :: fy_cand(:,:,:)
  real(mytype), allocatable :: fz_cand(:,:,:)
  real(mytype), allocatable :: rhs_inc_x(:,:,:)
  real(mytype), allocatable :: rhs_inc_y(:,:,:)
  real(mytype), allocatable :: rhs_inc_z(:,:,:)
  real(mytype), allocatable :: rhs_old_x(:,:,:)
  real(mytype), allocatable :: rhs_old_y(:,:,:)
  real(mytype), allocatable :: rhs_old_z(:,:,:)
  real(mytype), allocatable :: rhs_new_x(:,:,:)
  real(mytype), allocatable :: rhs_new_y(:,:,:)
  real(mytype), allocatable :: rhs_new_z(:,:,:)
  real(mytype) :: f_fluid_to_fibre(2,3)
  real(mytype) :: f_fibre_to_fluid(2,3)
  real(mytype) :: ds_lag(2)
  real(mytype) :: expected_fluid_to_fibre(3)
  real(mytype) :: expected_fibre_to_fluid(3)
  real(mytype) :: integrated_rhs_increment(3)
  real(mytype) :: correct_sign_error
  real(mytype) :: wrong_sign_error
  real(mytype) :: action_reaction_error
  real(mytype) :: additive_formula_error
  real(mytype) :: total_volume
  integer :: i
  integer :: j
  integer :: k
  integer :: requested_flag
  integer :: rhs_injection_enabled_flag
  integer :: injection_gain_finite_status
  integer :: initialized_status
  integer :: fluid_to_fibre_sign_status
  integer :: fibre_to_fluid_sign_status
  integer :: action_reaction_status
  integer :: rhs_increment_uses_fibre_to_fluid_status
  integer :: wrong_sign_rejection_status
  integer :: lambda_zero_scaling_status
  integer :: lambda_one_scaling_status
  integer :: lambda_fractional_scaling_status
  integer :: lambda_negative_scaling_status
  integer :: component_scaling_status
  integer :: integrated_rhs_sign_status
  integer :: integrated_rhs_scaling_status
  integer :: finite_rhs_increment_status
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
  integer :: rhs_sign_scaling_audit_status

  call execute_command_line('mkdir -p stage14_outputs')
  call stage14_config_load()
  requested_flag = logical_to_int_local(stage14_requested())
  rhs_injection_enabled_flag = logical_to_int_local(stage14_rhs_injection_enabled())
  injection_gain_finite_status = logical_to_int_local(finite_real_local(stage14_get_injection_gain()))

  allocate(volume_eul(nx,ny,nz), fx_cand(nx,ny,nz), fy_cand(nx,ny,nz), fz_cand(nx,ny,nz))
  allocate(rhs_inc_x(nx,ny,nz), rhs_inc_y(nx,ny,nz), rhs_inc_z(nx,ny,nz))
  allocate(rhs_old_x(nx,ny,nz), rhs_old_y(nx,ny,nz), rhs_old_z(nx,ny,nz))
  allocate(rhs_new_x(nx,ny,nz), rhs_new_y(nx,ny,nz), rhs_new_z(nx,ny,nz))

  f_fluid_to_fibre(1, :) = (/ 1.25_mytype, -0.50_mytype,  0.75_mytype /)
  f_fluid_to_fibre(2, :) = (/ 0.20_mytype,  0.40_mytype, -0.10_mytype /)
  f_fibre_to_fluid = -f_fluid_to_fibre
  ds_lag = (/ 1.0_mytype, 0.5_mytype /)

  total_volume = 0.0_mytype
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        volume_eul(i,j,k) = 1.0_mytype + 0.01_mytype * real(i, mytype) + &
                            0.02_mytype * real(j, mytype) + 0.03_mytype * real(k, mytype)
        rhs_old_x(i,j,k) = 0.10_mytype + 0.001_mytype * real(i, mytype)
        rhs_old_y(i,j,k) = -0.20_mytype + 0.002_mytype * real(j, mytype)
        rhs_old_z(i,j,k) = 0.30_mytype - 0.003_mytype * real(k, mytype)
        total_volume = total_volume + volume_eul(i,j,k)
      end do
    end do
  end do

  call stage14_rhs_sign_scaling_audit_init()
  call stage14_rhs_sign_scaling_audit_compute_expected_forces(f_fluid_to_fibre, f_fibre_to_fluid, &
                                                              ds_lag, expected_fluid_to_fibre, &
                                                              expected_fibre_to_fluid)
  call stage14_rhs_sign_scaling_audit_check_action_reaction(f_fluid_to_fibre, f_fibre_to_fluid, &
                                                            action_reaction_error)

  fx_cand = expected_fibre_to_fluid(1) / total_volume
  fy_cand = expected_fibre_to_fluid(2) / total_volume
  fz_cand = expected_fibre_to_fluid(3) / total_volume

  additive_formula_error = 0.0_mytype

  call stage14_rhs_accumulator_init(nx, ny, nz)
  call stage14_rhs_addition_formula_init()

  call run_lambda_case(0.0_mytype)
  call run_lambda_case(1.0_mytype)
  call run_lambda_case(0.1_mytype)
  call run_lambda_case(-0.1_mytype)

  call stage14_rhs_sign_scaling_audit_get_status_values(initialized_status, &
                                                        fluid_to_fibre_sign_status, &
                                                        fibre_to_fluid_sign_status, &
                                                        action_reaction_status, &
                                                        rhs_increment_uses_fibre_to_fluid_status, &
                                                        wrong_sign_rejection_status, &
                                                        lambda_zero_scaling_status, &
                                                        lambda_one_scaling_status, &
                                                        lambda_fractional_scaling_status, &
                                                        lambda_negative_scaling_status, &
                                                        component_scaling_status, &
                                                        integrated_rhs_sign_status, &
                                                        integrated_rhs_scaling_status, &
                                                        finite_rhs_increment_status, &
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
                                                        rhs_sign_scaling_audit_status)

  call stage14_rhs_sign_scaling_audit_write_diagnostics( &
       'stage14_outputs/fibre_stage14_3_rhs_sign_scaling_audit.dat', &
       requested_flag, rhs_injection_enabled_flag, injection_gain_finite_status)

  if (rhs_sign_scaling_audit_status == 1 .and. requested_flag == 1 .and. &
      rhs_injection_enabled_flag == 1 .and. injection_gain_finite_status == 1 .and. &
      additive_formula_error <= local_abs_tol) then
    print *, 'STAGE 14.3 RHS SIGN SCALING AUDIT VERDICT: PASS'
  else
    print *, 'STAGE 14.3 RHS SIGN SCALING AUDIT VERDICT: FAIL'
    if (requested_flag /= 1) print *, 'Reason: stage14_3_requested_flag'
    if (rhs_injection_enabled_flag /= 1) print *, 'Reason: stage14_3_rhs_injection_enabled_flag'
    if (injection_gain_finite_status /= 1) print *, 'Reason: stage14_3_injection_gain_finite_status'
    if (initialized_status /= 1) print *, 'Reason: stage14_3_initialized_status'
    if (fluid_to_fibre_sign_status /= 1) print *, 'Reason: stage14_3_fluid_to_fibre_sign_status'
    if (fibre_to_fluid_sign_status /= 1) print *, 'Reason: stage14_3_fibre_to_fluid_sign_status'
    if (action_reaction_status /= 1) print *, 'Reason: stage14_3_action_reaction_status'
    if (rhs_increment_uses_fibre_to_fluid_status /= 1) print *, 'Reason: stage14_3_rhs_increment_uses_fibre_to_fluid_status'
    if (wrong_sign_rejection_status /= 1) print *, 'Reason: stage14_3_wrong_sign_rejection_status'
    if (lambda_zero_scaling_status /= 1) print *, 'Reason: stage14_3_lambda_zero_scaling_status'
    if (lambda_one_scaling_status /= 1) print *, 'Reason: stage14_3_lambda_one_scaling_status'
    if (lambda_fractional_scaling_status /= 1) print *, 'Reason: stage14_3_lambda_fractional_scaling_status'
    if (lambda_negative_scaling_status /= 1) print *, 'Reason: stage14_3_lambda_negative_scaling_status'
    if (component_scaling_status /= 1) print *, 'Reason: stage14_3_component_scaling_status'
    if (integrated_rhs_sign_status /= 1) print *, 'Reason: stage14_3_integrated_rhs_sign_status'
    if (integrated_rhs_scaling_status /= 1) print *, 'Reason: stage14_3_integrated_rhs_scaling_status'
    if (finite_rhs_increment_status /= 1) print *, 'Reason: stage14_3_finite_rhs_increment_status'
    if (additive_formula_error > local_abs_tol) print *, 'Reason: stage14_3_additive_formula_consistency'
    if (no_production_rhs_modification_status /= 1) print *, 'Reason: stage14_3_no_production_rhs_modification_status'
    if (no_xcompact3d_hook_status /= 1) print *, 'Reason: stage14_3_no_xcompact3d_hook_status'
    if (no_fluid_field_access_status /= 1) print *, 'Reason: stage14_3_no_fluid_field_access_status'
    if (no_fluid_field_modification_status /= 1) print *, 'Reason: stage14_3_no_fluid_field_modification_status'
    if (no_pressure_modification_status /= 1) print *, 'Reason: stage14_3_no_pressure_modification_status'
    if (no_projection_modification_status /= 1) print *, 'Reason: stage14_3_no_projection_modification_status'
    if (no_rk3_modification_status /= 1) print *, 'Reason: stage14_3_no_rk3_modification_status'
    if (no_channel_forcing_modification_status /= 1) print *, 'Reason: stage14_3_no_channel_forcing_modification_status'
    if (no_production_ibm_forcing_status /= 1) print *, 'Reason: stage14_3_no_production_ibm_forcing_status'
    if (no_feedback_application_status /= 1) print *, 'Reason: stage14_3_no_feedback_application_status'
    if (no_twoway_force_status /= 1) print *, 'Reason: stage14_3_no_twoway_force_status'
    if (no_structure_advance_status /= 1) print *, 'Reason: stage14_3_no_structure_advance_status'
    if (rhs_sign_scaling_audit_status /= 1) print *, 'Reason: stage14_3_rhs_sign_scaling_audit_status'
    stop 1
  end if

  call stage14_rhs_accumulator_finalize()
  call stage14_rhs_addition_formula_finalize()
  call stage14_rhs_sign_scaling_audit_finalize()
  deallocate(volume_eul, fx_cand, fy_cand, fz_cand)
  deallocate(rhs_inc_x, rhs_inc_y, rhs_inc_z)
  deallocate(rhs_old_x, rhs_old_y, rhs_old_z)
  deallocate(rhs_new_x, rhs_new_y, rhs_new_z)

contains

  subroutine run_lambda_case(lambda_value)
    real(mytype), intent(in) :: lambda_value

    call stage14_rhs_accumulator_compute_from_force_density(fx_cand, fy_cand, fz_cand, lambda_value)
    rhs_inc_x = lambda_value * fx_cand
    rhs_inc_y = lambda_value * fy_cand
    rhs_inc_z = lambda_value * fz_cand
    integrated_rhs_increment(1) = sum(rhs_inc_x * volume_eul)
    integrated_rhs_increment(2) = sum(rhs_inc_y * volume_eul)
    integrated_rhs_increment(3) = sum(rhs_inc_z * volume_eul)
    call stage14_rhs_sign_scaling_audit_note_rhs_increment(rhs_inc_x, rhs_inc_y, rhs_inc_z, &
                                                           integrated_rhs_increment)
    call stage14_rhs_sign_scaling_audit_check_scaling(integrated_rhs_increment, expected_fibre_to_fluid, &
                                                       lambda_value)
    if (lambda_value == 1.0_mytype) then
      call stage14_rhs_sign_scaling_audit_check_rhs_increment_sign(integrated_rhs_increment, &
                                                                   expected_fluid_to_fibre, &
                                                                   expected_fibre_to_fluid, &
                                                                   lambda_value, &
                                                                   correct_sign_error, wrong_sign_error)
      call stage14_rhs_addition_formula_apply_controlled(rhs_old_x, rhs_old_y, rhs_old_z, &
                                                         rhs_inc_x, rhs_inc_y, rhs_inc_z, &
                                                         rhs_new_x, rhs_new_y, rhs_new_z)
      additive_formula_error = maxval(abs(rhs_new_x - rhs_old_x - rhs_inc_x))
      additive_formula_error = max(additive_formula_error, maxval(abs(rhs_new_y - rhs_old_y - rhs_inc_y)))
      additive_formula_error = max(additive_formula_error, maxval(abs(rhs_new_z - rhs_old_z - rhs_inc_z)))
    end if
  end subroutine run_lambda_case

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

end program fibre_stage14_rhs_sign_scaling_audit_check
