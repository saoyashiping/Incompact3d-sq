program fibre_stage13_conservation_sign_audit_check
  use decomp_2d_constants, only : mytype
  use fibre_stage13_config, only : stage13_config_load, stage13_get_status_values
  use fibre_stage13_force_density_buffer, only : stage13_force_density_buffer_init, &
                                                stage13_force_density_buffer_finalize
  use fibre_stage13_spreading_kernel, only : stage13_spreading_kernel_init, &
                                            stage13_spreading_kernel_clear, &
                                            stage13_spreading_kernel_spread_controlled, &
                                            stage13_spreading_kernel_finalize
  use fibre_stage13_volume_normalization_audit, only : stage13_volume_normalization_audit_init, &
                                                       stage13_volume_normalization_audit_check_volumes, &
                                                       stage13_volume_normalization_audit_integrate_force_density, &
                                                       stage13_volume_normalization_audit_finalize
  use fibre_stage13_conservation_sign_audit, only : stage13_conservation_sign_audit_init, &
                                                    stage13_conservation_sign_audit_clear, &
                                                    stage13_conservation_sign_audit_compute_expected_forces, &
                                                    stage13_conservation_sign_audit_check_action_reaction, &
                                                    stage13_conservation_sign_audit_check_integrated_force_sign, &
                                                    stage13_conservation_sign_audit_record_statuses, &
                                                    stage13_conservation_sign_audit_get_status_values, &
                                                    stage13_conservation_sign_audit_finalize
  implicit none

  integer, parameter :: nx = 6
  integer, parameter :: ny = 7
  integer, parameter :: nz = 6
  real(mytype), parameter :: action_reaction_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: force_conservation_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: component_conservation_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: wrong_sign_min_separation = 1.0e-8_mytype
  real(mytype), parameter :: zero_after_clear_abs_tol = 1.0e-12_mytype

  integer :: requested_flag
  integer :: readonly_status
  integer :: spreading_readonly_status
  integer :: disabled_default_status
  integer :: max_points_status
  integer :: max_eulerian_points_status
  integer :: normalization_status
  integer :: config_no_force_density_allocation_status
  integer :: config_no_spreading_status
  integer :: config_no_rhs_injection_status
  integer :: config_no_ibm_spreading_status
  integer :: config_no_feedback_application_status
  integer :: config_no_twoway_force_status
  integer :: config_no_structure_advance_status
  integer :: config_no_fluid_field_access_status
  integer :: config_no_fluid_field_modification_status
  integer :: config_status
  integer :: initialized_status
  integer :: fluid_to_fibre_sign_status
  integer :: fibre_to_fluid_sign_status
  integer :: action_reaction_status
  integer :: spreading_input_sign_status
  integer :: integrated_force_fibre_to_fluid_status
  integer :: wrong_sign_rejection_status
  integer :: component_sign_conservation_status
  integer :: multipoint_sign_conservation_status
  integer :: nonuniform_volume_sign_conservation_status
  integer :: boundary_sign_conservation_status
  integer :: finite_force_density_status
  integer :: clear_status
  integer :: no_rhs_injection_status
  integer :: no_ibm_spreading_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: audit_status
  integer :: verdict_status
  integer :: unit_id
  integer :: io_status
  integer :: i
  integer :: j
  integer :: k
  real(mytype) :: x_eul(nx)
  real(mytype) :: y_eul(ny)
  real(mytype) :: z_eul(nz)
  real(mytype) :: volume_eul(nx,ny,nz)
  real(mytype) :: fx_cand(nx,ny,nz)
  real(mytype) :: fy_cand(nx,ny,nz)
  real(mytype) :: fz_cand(nx,ny,nz)
  real(mytype) :: force_density_norm(nx,ny,nz)
  real(mytype) :: x_lag2(2)
  real(mytype) :: y_lag2(2)
  real(mytype) :: z_lag2(2)
  real(mytype) :: ds_lag2(2)
  real(mytype) :: f_fluid_to_fibre2(2,3)
  real(mytype) :: f_fibre_to_fluid2(2,3)
  real(mytype) :: x_lag1(1)
  real(mytype) :: y_lag1(1)
  real(mytype) :: z_lag1(1)
  real(mytype) :: ds_lag1(1)
  real(mytype) :: f_fluid_to_fibre1(1,3)
  real(mytype) :: f_fibre_to_fluid1(1,3)
  real(mytype) :: expected_fluid_to_fibre(3)
  real(mytype) :: expected_fibre_to_fluid(3)
  real(mytype) :: expected_boundary_fluid(3)
  real(mytype) :: expected_boundary_fibre(3)
  real(mytype) :: integrated_force(3)
  real(mytype) :: min_volume
  real(mytype) :: max_volume
  real(mytype) :: action_reaction_max_error
  real(mytype) :: fibre_to_fluid_error
  real(mytype) :: wrong_sign_error
  real(mytype) :: integrated_force_fibre_to_fluid_error_l2
  real(mytype) :: wrong_sign_error_l2
  real(mytype) :: wrong_sign_separation
  real(mytype) :: component_error_x
  real(mytype) :: component_error_y
  real(mytype) :: component_error_z
  real(mytype) :: nonuniform_sign_error_l2
  real(mytype) :: boundary_sign_error_l2
  real(mytype) :: max_abs_force_density
  real(mytype) :: max_abs_force_density_norm_after_clear

  do i = 1, nx
    x_eul(i) = real(i - 1, mytype)
  end do
  do j = 1, ny
    y_eul(j) = real(j - 1, mytype)
  end do
  do k = 1, nz
    z_eul(k) = real(k - 1, mytype)
  end do
  volume_eul(:,:,:) = 1.0_mytype
  fx_cand(:,:,:) = 0.0_mytype
  fy_cand(:,:,:) = 0.0_mytype
  fz_cand(:,:,:) = 0.0_mytype
  force_density_norm(:,:,:) = 0.0_mytype
  max_abs_force_density = 0.0_mytype
  max_abs_force_density_norm_after_clear = 0.0_mytype
  finite_force_density_status = 1

  call stage13_config_load()
  call stage13_get_status_values(requested_flag, readonly_status, spreading_readonly_status, disabled_default_status, &
                                 max_points_status, max_eulerian_points_status, normalization_status, &
                                 config_no_force_density_allocation_status, config_no_spreading_status, &
                                 config_no_rhs_injection_status, config_no_ibm_spreading_status, &
                                 config_no_feedback_application_status, config_no_twoway_force_status, &
                                 config_no_structure_advance_status, config_no_fluid_field_access_status, &
                                 config_no_fluid_field_modification_status, config_status)

  call stage13_force_density_buffer_init(nx, ny, nz)
  call stage13_spreading_kernel_init()
  call stage13_volume_normalization_audit_init()
  call stage13_conservation_sign_audit_init()

  x_lag2 = (/ 2.25_mytype, 3.25_mytype /)
  y_lag2 = (/ 2.50_mytype, 4.25_mytype /)
  z_lag2 = (/ 1.25_mytype, 2.75_mytype /)
  ds_lag2 = (/ 1.0_mytype, 0.5_mytype /)
  f_fluid_to_fibre2(1,1) = 1.25_mytype
  f_fluid_to_fibre2(1,2) = -0.50_mytype
  f_fluid_to_fibre2(1,3) = 0.75_mytype
  f_fluid_to_fibre2(2,1) = 0.20_mytype
  f_fluid_to_fibre2(2,2) = 0.40_mytype
  f_fluid_to_fibre2(2,3) = -0.10_mytype
  f_fibre_to_fluid2(:,:) = -f_fluid_to_fibre2(:,:)

  call stage13_conservation_sign_audit_compute_expected_forces(f_fluid_to_fibre2, f_fibre_to_fluid2, ds_lag2, &
                                                               expected_fluid_to_fibre, expected_fibre_to_fluid)
  call stage13_conservation_sign_audit_check_action_reaction(f_fluid_to_fibre2, f_fibre_to_fluid2, &
                                                             action_reaction_max_error)

  call run_sign_case(x_lag2, y_lag2, z_lag2, ds_lag2, f_fibre_to_fluid2, volume_eul, expected_fibre_to_fluid, &
                     integrated_force)
  call stage13_conservation_sign_audit_check_integrated_force_sign(integrated_force, expected_fluid_to_fibre, &
                                                                   expected_fibre_to_fluid, fibre_to_fluid_error, &
                                                                   wrong_sign_error)
  integrated_force_fibre_to_fluid_error_l2 = fibre_to_fluid_error
  wrong_sign_error_l2 = wrong_sign_error
  wrong_sign_separation = wrong_sign_error_l2 - integrated_force_fibre_to_fluid_error_l2
  component_error_x = abs(integrated_force(1) - expected_fibre_to_fluid(1))
  component_error_y = abs(integrated_force(2) - expected_fibre_to_fluid(2))
  component_error_z = abs(integrated_force(3) - expected_fibre_to_fluid(3))

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        volume_eul(i,j,k) = 1.0_mytype + 0.01_mytype * real(i, mytype) + &
                            0.02_mytype * real(j, mytype) + 0.03_mytype * real(k, mytype)
      end do
    end do
  end do
  call run_sign_case(x_lag2, y_lag2, z_lag2, ds_lag2, f_fibre_to_fluid2, volume_eul, expected_fibre_to_fluid, &
                     integrated_force)
  nonuniform_sign_error_l2 = vector_l2(integrated_force - expected_fibre_to_fluid)

  x_lag1 = (/ 0.25_mytype /)
  y_lag1 = (/ 1.25_mytype /)
  z_lag1 = (/ 0.25_mytype /)
  ds_lag1 = (/ 1.20_mytype /)
  f_fluid_to_fibre1(1,1) = 0.30_mytype
  f_fluid_to_fibre1(1,2) = -0.40_mytype
  f_fluid_to_fibre1(1,3) = 0.50_mytype
  f_fibre_to_fluid1(:,:) = -f_fluid_to_fibre1(:,:)
  call stage13_conservation_sign_audit_compute_expected_forces(f_fluid_to_fibre1, f_fibre_to_fluid1, ds_lag1, &
                                                               expected_boundary_fluid, expected_boundary_fibre)
  call run_sign_case(x_lag1, y_lag1, z_lag1, ds_lag1, f_fibre_to_fluid1, volume_eul, expected_boundary_fibre, &
                     integrated_force)
  boundary_sign_error_l2 = vector_l2(integrated_force - expected_boundary_fibre)

  call stage13_spreading_kernel_clear(fx_cand, fy_cand, fz_cand, force_density_norm)
  max_abs_force_density_norm_after_clear = maxval(abs(force_density_norm))
  clear_status = merge(1, 0, max_abs_force_density_norm_after_clear <= zero_after_clear_abs_tol)
  finite_force_density_status = min(finite_force_density_status, &
       merge(1, 0, all_finite_rank3(fx_cand) .and. all_finite_rank3(fy_cand) .and. &
                   all_finite_rank3(fz_cand) .and. all_finite_rank3(force_density_norm)))

  call stage13_conservation_sign_audit_clear(clear_status == 1)
  call stage13_conservation_sign_audit_record_statuses(1, &
       merge(1, 0, component_error_x <= component_conservation_abs_tol .and. &
                   component_error_y <= component_conservation_abs_tol .and. &
                   component_error_z <= component_conservation_abs_tol), &
       merge(1, 0, integrated_force_fibre_to_fluid_error_l2 <= force_conservation_abs_tol), &
       merge(1, 0, nonuniform_sign_error_l2 <= force_conservation_abs_tol), &
       merge(1, 0, boundary_sign_error_l2 <= force_conservation_abs_tol), &
       finite_force_density_status, clear_status)
  call stage13_conservation_sign_audit_get_status_values(initialized_status, fluid_to_fibre_sign_status, &
                                                         fibre_to_fluid_sign_status, action_reaction_status, &
                                                         spreading_input_sign_status, &
                                                         integrated_force_fibre_to_fluid_status, &
                                                         wrong_sign_rejection_status, &
                                                         component_sign_conservation_status, &
                                                         multipoint_sign_conservation_status, &
                                                         nonuniform_volume_sign_conservation_status, &
                                                         boundary_sign_conservation_status, &
                                                         finite_force_density_status, clear_status, &
                                                         no_rhs_injection_status, no_ibm_spreading_status, &
                                                         no_feedback_application_status, no_twoway_force_status, &
                                                         no_structure_advance_status, no_fluid_field_access_status, &
                                                         no_fluid_field_modification_status, audit_status)

  verdict_status = 1
  if (requested_flag /= 1) verdict_status = 0
  if (readonly_status /= 1) verdict_status = 0
  if (spreading_readonly_status /= 1) verdict_status = 0
  if (initialized_status /= 1) verdict_status = 0
  if (fluid_to_fibre_sign_status /= 1) verdict_status = 0
  if (fibre_to_fluid_sign_status /= 1) verdict_status = 0
  if (action_reaction_status /= 1) verdict_status = 0
  if (spreading_input_sign_status /= 1) verdict_status = 0
  if (integrated_force_fibre_to_fluid_status /= 1) verdict_status = 0
  if (wrong_sign_rejection_status /= 1) verdict_status = 0
  if (component_sign_conservation_status /= 1) verdict_status = 0
  if (multipoint_sign_conservation_status /= 1) verdict_status = 0
  if (nonuniform_volume_sign_conservation_status /= 1) verdict_status = 0
  if (boundary_sign_conservation_status /= 1) verdict_status = 0
  if (finite_force_density_status /= 1) verdict_status = 0
  if (clear_status /= 1) verdict_status = 0
  if (no_rhs_injection_status /= 1) verdict_status = 0
  if (no_ibm_spreading_status /= 1) verdict_status = 0
  if (no_feedback_application_status /= 1) verdict_status = 0
  if (no_twoway_force_status /= 1) verdict_status = 0
  if (no_structure_advance_status /= 1) verdict_status = 0
  if (no_fluid_field_access_status /= 1) verdict_status = 0
  if (no_fluid_field_modification_status /= 1) verdict_status = 0
  if (audit_status /= 1) verdict_status = 0
  if (action_reaction_max_error > action_reaction_abs_tol) verdict_status = 0
  if (integrated_force_fibre_to_fluid_error_l2 > force_conservation_abs_tol) verdict_status = 0
  if (component_error_x > component_conservation_abs_tol) verdict_status = 0
  if (component_error_y > component_conservation_abs_tol) verdict_status = 0
  if (component_error_z > component_conservation_abs_tol) verdict_status = 0
  if (nonuniform_sign_error_l2 > force_conservation_abs_tol) verdict_status = 0
  if (boundary_sign_error_l2 > force_conservation_abs_tol) verdict_status = 0
  if (wrong_sign_separation < wrong_sign_min_separation) verdict_status = 0
  if (max_abs_force_density_norm_after_clear > zero_after_clear_abs_tol) verdict_status = 0

  open(newunit=unit_id, file='stage13_outputs/fibre_stage13_5_conservation_sign_audit.dat', status='replace', &
       action='write', iostat=io_status)
  if (io_status /= 0) then
    write(*,'(A)') 'STAGE 13.5 CONSERVATION SIGN AUDIT VERDICT: FAIL'
    write(*,'(A)') 'Reason: unable_to_open_stage13_outputs_fibre_stage13_5_conservation_sign_audit_dat'
    stop 1
  end if

  write(unit_id,'(A,1X,I0)') 'stage13_5_requested_flag', requested_flag
  write(unit_id,'(A,1X,I0)') 'stage13_5_readonly_mode_status', readonly_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_spreading_readonly_status', spreading_readonly_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_initialized_status', initialized_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_fluid_to_fibre_sign_status', fluid_to_fibre_sign_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_fibre_to_fluid_sign_status', fibre_to_fluid_sign_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_action_reaction_status', action_reaction_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_spreading_input_sign_status', spreading_input_sign_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_integrated_force_fibre_to_fluid_status', &
       integrated_force_fibre_to_fluid_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_wrong_sign_rejection_status', wrong_sign_rejection_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_component_sign_conservation_status', component_sign_conservation_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_multipoint_sign_conservation_status', multipoint_sign_conservation_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_nonuniform_volume_sign_conservation_status', &
       nonuniform_volume_sign_conservation_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_boundary_sign_conservation_status', boundary_sign_conservation_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_finite_force_density_status', finite_force_density_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_clear_status', clear_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_no_rhs_injection_status', no_rhs_injection_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_no_ibm_spreading_status', no_ibm_spreading_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_no_feedback_application_status', no_feedback_application_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_no_twoway_force_status', no_twoway_force_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_no_structure_advance_status', no_structure_advance_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_no_fluid_field_access_status', no_fluid_field_access_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_no_fluid_field_modification_status', no_fluid_field_modification_status
  write(unit_id,'(A,1X,I0)') 'stage13_5_conservation_sign_audit_status', audit_status
  write(unit_id,'(A,1X,ES24.16)') 'stage13_5_action_reaction_max_error', action_reaction_max_error
  write(unit_id,'(A,1X,ES24.16)') 'stage13_5_integrated_force_fibre_to_fluid_error_l2', &
       integrated_force_fibre_to_fluid_error_l2
  write(unit_id,'(A,1X,ES24.16)') 'stage13_5_wrong_sign_error_l2', wrong_sign_error_l2
  write(unit_id,'(A,1X,ES24.16)') 'stage13_5_wrong_sign_separation', wrong_sign_separation
  write(unit_id,'(A,1X,ES24.16)') 'stage13_5_component_error_x', component_error_x
  write(unit_id,'(A,1X,ES24.16)') 'stage13_5_component_error_y', component_error_y
  write(unit_id,'(A,1X,ES24.16)') 'stage13_5_component_error_z', component_error_z
  write(unit_id,'(A,1X,ES24.16)') 'stage13_5_nonuniform_sign_error_l2', nonuniform_sign_error_l2
  write(unit_id,'(A,1X,ES24.16)') 'stage13_5_boundary_sign_error_l2', boundary_sign_error_l2
  write(unit_id,'(A,1X,ES24.16)') 'stage13_5_max_abs_force_density', max_abs_force_density
  write(unit_id,'(A,1X,ES24.16)') 'stage13_5_max_abs_force_density_norm_after_clear', &
       max_abs_force_density_norm_after_clear
  close(unit_id)

  call stage13_conservation_sign_audit_finalize()
  call stage13_volume_normalization_audit_finalize()
  call stage13_spreading_kernel_finalize()
  call stage13_force_density_buffer_finalize()

  if (verdict_status == 1) then
    write(*,'(A)') 'STAGE 13.5 CONSERVATION SIGN AUDIT VERDICT: PASS'
  else
    write(*,'(A)') 'STAGE 13.5 CONSERVATION SIGN AUDIT VERDICT: FAIL'
    write(*,'(A,1X,ES24.16)') 'Reason: stage13_5_integrated_force_fibre_to_fluid_error_l2', &
         integrated_force_fibre_to_fluid_error_l2
    write(*,'(A,1X,ES24.16)') 'Reason: stage13_5_wrong_sign_separation', wrong_sign_separation
    stop 1
  end if

contains

  subroutine run_sign_case(x_lag, y_lag, z_lag, ds_lag, f_spread, volume_case, expected, integrated)
    real(mytype), intent(in) :: x_lag(:)
    real(mytype), intent(in) :: y_lag(:)
    real(mytype), intent(in) :: z_lag(:)
    real(mytype), intent(in) :: ds_lag(:)
    real(mytype), intent(in) :: f_spread(:,:)
    real(mytype), intent(in) :: volume_case(:,:,:)
    real(mytype), intent(in) :: expected(3)
    real(mytype), intent(out) :: integrated(3)
    real(mytype) :: local_min_volume
    real(mytype) :: local_max_volume

    call stage13_spreading_kernel_clear(fx_cand, fy_cand, fz_cand, force_density_norm)
    call stage13_volume_normalization_audit_check_volumes(volume_case, local_min_volume, local_max_volume)
    call stage13_spreading_kernel_spread_controlled(x_lag, y_lag, z_lag, ds_lag, f_spread, x_eul, y_eul, z_eul, &
                                                    volume_case, fx_cand, fy_cand, fz_cand, force_density_norm)
    call stage13_volume_normalization_audit_integrate_force_density(fx_cand, fy_cand, fz_cand, volume_case, &
                                                                    integrated)
    max_abs_force_density = max(max_abs_force_density, maxval(abs(fx_cand)), maxval(abs(fy_cand)), &
                                maxval(abs(fz_cand)))
    finite_force_density_status = min(finite_force_density_status, &
         merge(1, 0, all_finite_rank3(fx_cand) .and. all_finite_rank3(fy_cand) .and. &
                     all_finite_rank3(fz_cand) .and. all_finite_rank3(force_density_norm)))
  end subroutine run_sign_case

  real(mytype) function vector_l2(values)
    real(mytype), intent(in) :: values(3)

    vector_l2 = sqrt(sum(values * values))
  end function vector_l2

  logical function all_finite_rank3(values)
    real(mytype), intent(in) :: values(:,:,:)

    all_finite_rank3 = all(values == values) .and. all(abs(values) < huge(1.0_mytype))
  end function all_finite_rank3

end program fibre_stage13_conservation_sign_audit_check
