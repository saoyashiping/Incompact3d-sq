program fibre_stage13_spreading_kernel_check
  use decomp_2d_constants, only : mytype
  use fibre_stage13_config, only : stage13_config_load, stage13_get_status_values
  use fibre_stage13_force_density_buffer, only : stage13_force_density_buffer_init, &
                                                stage13_force_density_buffer_finalize
  use fibre_stage13_spreading_kernel, only : stage13_spreading_kernel_init, &
                                            stage13_spreading_kernel_clear, &
                                            stage13_spreading_kernel_spread_controlled, &
                                            stage13_spreading_kernel_get_status_values, &
                                            stage13_spreading_kernel_get_diagnostics, &
                                            stage13_spreading_kernel_finalize
  implicit none

  integer, parameter :: nx = 6
  integer, parameter :: ny = 7
  integer, parameter :: nz = 6
  real(mytype), parameter :: weight_sum_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: force_conservation_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: zero_force_abs_tol = 1.0e-12_mytype

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
  integer :: zero_force_spreading_status
  integer :: single_point_spreading_status
  integer :: compact_support_status
  integer :: weight_normalization_status
  integer :: nonnegative_weight_status
  integer :: component_spreading_status
  integer :: boundary_safe_status
  integer :: integrated_force_conservation_status
  integer :: finite_force_density_status
  integer :: force_density_norm_finite_status
  integer :: clear_status
  integer :: no_rhs_injection_status
  integer :: no_ibm_spreading_status
  integer :: no_feedback_application_status
  integer :: no_twoway_force_status
  integer :: no_structure_advance_status
  integer :: no_fluid_field_access_status
  integer :: no_fluid_field_modification_status
  integer :: kernel_status_tmp
  integer :: clear_initialized_tmp
  integer :: clear_zero_tmp
  integer :: clear_single_tmp
  integer :: clear_compact_tmp
  integer :: clear_weight_norm_tmp
  integer :: clear_nonnegative_tmp
  integer :: clear_component_tmp
  integer :: clear_boundary_tmp
  integer :: clear_integrated_tmp
  integer :: clear_finite_tmp
  integer :: clear_norm_finite_tmp
  integer :: clear_no_rhs_tmp
  integer :: clear_no_ibm_tmp
  integer :: clear_no_feedback_tmp
  integer :: clear_no_twoway_tmp
  integer :: clear_no_structure_tmp
  integer :: clear_no_fluid_access_tmp
  integer :: clear_no_fluid_mod_tmp
  integer :: final_kernel_status
  integer :: support_count_tmp
  integer :: support_count_after_clear_tmp
  integer :: unit_id
  integer :: io_status
  integer :: verdict_status
  integer :: i
  integer :: j
  integer :: k
  real(mytype) :: zero_force_max_abs
  real(mytype) :: zero_max_tmp
  real(mytype) :: weight_error_tmp
  real(mytype) :: min_weight_tmp
  real(mytype) :: error_x_tmp
  real(mytype) :: error_y_tmp
  real(mytype) :: error_z_tmp
  real(mytype) :: error_l2_tmp
  real(mytype) :: max_force_tmp
  real(mytype) :: norm_after_clear_tmp
  real(mytype) :: max_weight_sum_error
  real(mytype) :: min_weight
  real(mytype) :: integrated_force_error_x
  real(mytype) :: integrated_force_error_y
  real(mytype) :: integrated_force_error_z
  real(mytype) :: integrated_force_error_l2
  real(mytype) :: max_abs_force_density
  real(mytype) :: max_abs_force_density_norm_after_clear
  real(mytype) :: x_eul(nx)
  real(mytype) :: y_eul(ny)
  real(mytype) :: z_eul(nz)
  real(mytype) :: volume_eul(nx,ny,nz)
  real(mytype) :: fx_cand(nx,ny,nz)
  real(mytype) :: fy_cand(nx,ny,nz)
  real(mytype) :: fz_cand(nx,ny,nz)
  real(mytype) :: force_density_norm(nx,ny,nz)
  real(mytype) :: x_lag1(1)
  real(mytype) :: y_lag1(1)
  real(mytype) :: z_lag1(1)
  real(mytype) :: ds_lag1(1)
  real(mytype) :: f_lag1(1,3)
  real(mytype) :: x_lag2(2)
  real(mytype) :: y_lag2(2)
  real(mytype) :: z_lag2(2)
  real(mytype) :: ds_lag2(2)
  real(mytype) :: f_lag2(2,3)
  real(mytype) :: x_lag3(3)
  real(mytype) :: y_lag3(3)
  real(mytype) :: z_lag3(3)
  real(mytype) :: ds_lag3(3)
  real(mytype) :: f_lag3(3,3)

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

  initialized_status = 1
  zero_force_spreading_status = 0
  single_point_spreading_status = 0
  compact_support_status = 1
  weight_normalization_status = 1
  nonnegative_weight_status = 1
  component_spreading_status = 1
  boundary_safe_status = 1
  integrated_force_conservation_status = 1
  finite_force_density_status = 1
  force_density_norm_finite_status = 1
  clear_status = 0
  no_rhs_injection_status = 1
  no_ibm_spreading_status = 1
  no_feedback_application_status = 1
  no_twoway_force_status = 1
  no_structure_advance_status = 1
  no_fluid_field_access_status = 1
  no_fluid_field_modification_status = 1
  zero_force_max_abs = 0.0_mytype
  max_weight_sum_error = 0.0_mytype
  min_weight = huge(1.0_mytype)
  integrated_force_error_x = 0.0_mytype
  integrated_force_error_y = 0.0_mytype
  integrated_force_error_z = 0.0_mytype
  integrated_force_error_l2 = 0.0_mytype
  max_abs_force_density = 0.0_mytype
  max_abs_force_density_norm_after_clear = 0.0_mytype
  support_count_tmp = 0

  x_lag2 = (/ 2.25_mytype, 3.25_mytype /)
  y_lag2 = (/ 2.50_mytype, 3.50_mytype /)
  z_lag2 = (/ 1.25_mytype, 2.25_mytype /)
  ds_lag2 = (/ 1.0_mytype, 0.5_mytype /)
  f_lag2(:,:) = 0.0_mytype
  call stage13_spreading_kernel_clear(fx_cand, fy_cand, fz_cand, force_density_norm)
  call stage13_spreading_kernel_spread_controlled(x_lag2, y_lag2, z_lag2, ds_lag2, f_lag2, &
                                                  x_eul, y_eul, z_eul, volume_eul, &
                                                  fx_cand, fy_cand, fz_cand, force_density_norm)
  call collect_kernel_status(zero_force_max_abs, max_weight_sum_error, min_weight, integrated_force_error_x, &
                             integrated_force_error_y, integrated_force_error_z, integrated_force_error_l2, &
                             max_abs_force_density, max_abs_force_density_norm_after_clear, support_count_tmp)

  x_lag1 = (/ 2.25_mytype /)
  y_lag1 = (/ 3.25_mytype /)
  z_lag1 = (/ 2.25_mytype /)
  ds_lag1 = (/ 1.0_mytype /)
  f_lag1(1,1) = 1.25_mytype
  f_lag1(1,2) = -0.50_mytype
  f_lag1(1,3) = 0.75_mytype
  call stage13_spreading_kernel_clear(fx_cand, fy_cand, fz_cand, force_density_norm)
  call stage13_spreading_kernel_spread_controlled(x_lag1, y_lag1, z_lag1, ds_lag1, f_lag1, &
                                                  x_eul, y_eul, z_eul, volume_eul, &
                                                  fx_cand, fy_cand, fz_cand, force_density_norm)
  call collect_kernel_status(zero_force_max_abs, max_weight_sum_error, min_weight, integrated_force_error_x, &
                             integrated_force_error_y, integrated_force_error_z, integrated_force_error_l2, &
                             max_abs_force_density, max_abs_force_density_norm_after_clear, support_count_tmp)

  x_lag1 = (/ 0.25_mytype /)
  y_lag1 = (/ 1.25_mytype /)
  z_lag1 = (/ 0.25_mytype /)
  ds_lag1 = (/ 1.0_mytype /)
  f_lag1(1,1) = 0.25_mytype
  f_lag1(1,2) = 0.50_mytype
  f_lag1(1,3) = -0.75_mytype
  call stage13_spreading_kernel_clear(fx_cand, fy_cand, fz_cand, force_density_norm)
  call stage13_spreading_kernel_spread_controlled(x_lag1, y_lag1, z_lag1, ds_lag1, f_lag1, &
                                                  x_eul, y_eul, z_eul, volume_eul, &
                                                  fx_cand, fy_cand, fz_cand, force_density_norm)
  call collect_kernel_status(zero_force_max_abs, max_weight_sum_error, min_weight, integrated_force_error_x, &
                             integrated_force_error_y, integrated_force_error_z, integrated_force_error_l2, &
                             max_abs_force_density, max_abs_force_density_norm_after_clear, support_count_tmp)

  x_lag3 = (/ 1.25_mytype, 3.50_mytype, 0.25_mytype /)
  y_lag3 = (/ 2.25_mytype, 4.50_mytype, 1.25_mytype /)
  z_lag3 = (/ 1.25_mytype, 2.50_mytype, 0.25_mytype /)
  ds_lag3 = (/ 0.50_mytype, 1.25_mytype, 0.75_mytype /)
  f_lag3(1,1) = 0.75_mytype
  f_lag3(1,2) = -0.25_mytype
  f_lag3(1,3) = 0.50_mytype
  f_lag3(2,1) = -0.60_mytype
  f_lag3(2,2) = 0.40_mytype
  f_lag3(2,3) = 0.20_mytype
  f_lag3(3,1) = 0.10_mytype
  f_lag3(3,2) = 0.30_mytype
  f_lag3(3,3) = -0.20_mytype
  call stage13_spreading_kernel_clear(fx_cand, fy_cand, fz_cand, force_density_norm)
  call stage13_spreading_kernel_spread_controlled(x_lag3, y_lag3, z_lag3, ds_lag3, f_lag3, &
                                                  x_eul, y_eul, z_eul, volume_eul, &
                                                  fx_cand, fy_cand, fz_cand, force_density_norm)
  call collect_kernel_status(zero_force_max_abs, max_weight_sum_error, min_weight, integrated_force_error_x, &
                             integrated_force_error_y, integrated_force_error_z, integrated_force_error_l2, &
                             max_abs_force_density, max_abs_force_density_norm_after_clear, support_count_tmp)

  call stage13_spreading_kernel_clear(fx_cand, fy_cand, fz_cand, force_density_norm)
  call stage13_spreading_kernel_get_diagnostics(zero_max_tmp, support_count_after_clear_tmp, weight_error_tmp, min_weight_tmp, &
                                                error_x_tmp, error_y_tmp, error_z_tmp, error_l2_tmp, &
                                                max_force_tmp, norm_after_clear_tmp)
  max_abs_force_density_norm_after_clear = norm_after_clear_tmp
  call stage13_spreading_kernel_get_status_values(clear_initialized_tmp, clear_zero_tmp, clear_single_tmp, &
                                                  clear_compact_tmp, clear_weight_norm_tmp, &
                                                  clear_nonnegative_tmp, clear_component_tmp, &
                                                  clear_boundary_tmp, clear_integrated_tmp, clear_finite_tmp, &
                                                  clear_norm_finite_tmp, clear_status, clear_no_rhs_tmp, &
                                                  clear_no_ibm_tmp, clear_no_feedback_tmp, clear_no_twoway_tmp, &
                                                  clear_no_structure_tmp, clear_no_fluid_access_tmp, &
                                                  clear_no_fluid_mod_tmp, kernel_status_tmp)
  initialized_status = min(initialized_status, clear_initialized_tmp)
  no_rhs_injection_status = min(no_rhs_injection_status, clear_no_rhs_tmp)
  no_ibm_spreading_status = min(no_ibm_spreading_status, clear_no_ibm_tmp)
  no_feedback_application_status = min(no_feedback_application_status, clear_no_feedback_tmp)
  no_twoway_force_status = min(no_twoway_force_status, clear_no_twoway_tmp)
  no_structure_advance_status = min(no_structure_advance_status, clear_no_structure_tmp)
  no_fluid_field_access_status = min(no_fluid_field_access_status, clear_no_fluid_access_tmp)
  no_fluid_field_modification_status = min(no_fluid_field_modification_status, clear_no_fluid_mod_tmp)
  finite_force_density_status = min(finite_force_density_status, clear_finite_tmp)
  force_density_norm_finite_status = min(force_density_norm_finite_status, clear_norm_finite_tmp)

  if (min_weight == huge(1.0_mytype)) min_weight = 0.0_mytype
  final_kernel_status = 1
  if (requested_flag /= 1) final_kernel_status = 0
  if (readonly_status /= 1) final_kernel_status = 0
  if (spreading_readonly_status /= 1) final_kernel_status = 0
  if (initialized_status /= 1) final_kernel_status = 0
  if (zero_force_spreading_status /= 1) final_kernel_status = 0
  if (single_point_spreading_status /= 1) final_kernel_status = 0
  if (compact_support_status /= 1) final_kernel_status = 0
  if (weight_normalization_status /= 1) final_kernel_status = 0
  if (nonnegative_weight_status /= 1) final_kernel_status = 0
  if (component_spreading_status /= 1) final_kernel_status = 0
  if (boundary_safe_status /= 1) final_kernel_status = 0
  if (integrated_force_conservation_status /= 1) final_kernel_status = 0
  if (finite_force_density_status /= 1) final_kernel_status = 0
  if (force_density_norm_finite_status /= 1) final_kernel_status = 0
  if (clear_status /= 1) final_kernel_status = 0
  if (no_rhs_injection_status /= 1) final_kernel_status = 0
  if (no_ibm_spreading_status /= 1) final_kernel_status = 0
  if (no_feedback_application_status /= 1) final_kernel_status = 0
  if (no_twoway_force_status /= 1) final_kernel_status = 0
  if (no_structure_advance_status /= 1) final_kernel_status = 0
  if (no_fluid_field_access_status /= 1) final_kernel_status = 0
  if (no_fluid_field_modification_status /= 1) final_kernel_status = 0
  if (zero_force_max_abs > zero_force_abs_tol) final_kernel_status = 0
  if (max_weight_sum_error > weight_sum_abs_tol) final_kernel_status = 0
  if (integrated_force_error_l2 > force_conservation_abs_tol) final_kernel_status = 0
  if (max_abs_force_density_norm_after_clear > zero_force_abs_tol) final_kernel_status = 0
  if (support_count_tmp > 8) final_kernel_status = 0

  open(newunit=unit_id, file='stage13_outputs/fibre_stage13_3_spreading_kernel.dat', status='replace', &
       action='write', iostat=io_status)
  if (io_status /= 0) then
    write(*,'(A)') 'STAGE 13.3 SPREADING KERNEL VERDICT: FAIL'
    write(*,'(A)') 'Reason: unable_to_open_stage13_outputs_fibre_stage13_3_spreading_kernel_dat'
    stop 1
  end if

  write(unit_id,'(A,1X,I0)') 'stage13_3_requested_flag', requested_flag
  write(unit_id,'(A,1X,I0)') 'stage13_3_readonly_mode_status', readonly_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_spreading_readonly_status', spreading_readonly_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_initialized_status', initialized_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_zero_force_spreading_status', zero_force_spreading_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_single_point_spreading_status', single_point_spreading_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_compact_support_status', compact_support_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_weight_normalization_status', weight_normalization_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_nonnegative_weight_status', nonnegative_weight_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_component_spreading_status', component_spreading_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_boundary_safe_status', boundary_safe_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_integrated_force_conservation_status', &
       integrated_force_conservation_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_finite_force_density_status', finite_force_density_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_force_density_norm_finite_status', force_density_norm_finite_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_clear_status', clear_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_no_rhs_injection_status', no_rhs_injection_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_no_ibm_spreading_status', no_ibm_spreading_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_no_feedback_application_status', no_feedback_application_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_no_twoway_force_status', no_twoway_force_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_no_structure_advance_status', no_structure_advance_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_no_fluid_field_access_status', no_fluid_field_access_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_no_fluid_field_modification_status', no_fluid_field_modification_status
  write(unit_id,'(A,1X,I0)') 'stage13_3_spreading_kernel_status', final_kernel_status
  write(unit_id,'(A,1X,ES24.16)') 'stage13_3_zero_force_max_abs', zero_force_max_abs
  write(unit_id,'(A,1X,I0)') 'stage13_3_single_point_support_count', support_count_tmp
  write(unit_id,'(A,1X,ES24.16)') 'stage13_3_max_weight_sum_error', max_weight_sum_error
  write(unit_id,'(A,1X,ES24.16)') 'stage13_3_min_weight', min_weight
  write(unit_id,'(A,1X,ES24.16)') 'stage13_3_integrated_force_error_x', integrated_force_error_x
  write(unit_id,'(A,1X,ES24.16)') 'stage13_3_integrated_force_error_y', integrated_force_error_y
  write(unit_id,'(A,1X,ES24.16)') 'stage13_3_integrated_force_error_z', integrated_force_error_z
  write(unit_id,'(A,1X,ES24.16)') 'stage13_3_integrated_force_error_l2', integrated_force_error_l2
  write(unit_id,'(A,1X,ES24.16)') 'stage13_3_max_abs_force_density', max_abs_force_density
  write(unit_id,'(A,1X,ES24.16)') 'stage13_3_max_abs_force_density_norm_after_clear', &
       max_abs_force_density_norm_after_clear
  close(unit_id)

  verdict_status = final_kernel_status
  call stage13_spreading_kernel_finalize()
  call stage13_force_density_buffer_finalize()

  if (verdict_status == 1) then
    write(*,'(A)') 'STAGE 13.3 SPREADING KERNEL VERDICT: PASS'
  else
    write(*,'(A)') 'STAGE 13.3 SPREADING KERNEL VERDICT: FAIL'
    write(*,'(A,1X,ES24.16)') 'Reason: stage13_3_spreading_kernel_status_failure max_weight_sum_error', &
         max_weight_sum_error
    write(*,'(A,1X,ES24.16)') 'Reason: stage13_3_spreading_kernel_status_failure integrated_force_error_l2', &
         integrated_force_error_l2
    stop 1
  end if

contains

  subroutine collect_kernel_status(zero_max, weight_error, weight_min, error_x, error_y, error_z, error_l2, &
                                   max_force, norm_after_clear, support_count)
    real(mytype), intent(inout) :: zero_max
    real(mytype), intent(inout) :: weight_error
    real(mytype), intent(inout) :: weight_min
    real(mytype), intent(inout) :: error_x
    real(mytype), intent(inout) :: error_y
    real(mytype), intent(inout) :: error_z
    real(mytype), intent(inout) :: error_l2
    real(mytype), intent(inout) :: max_force
    real(mytype), intent(inout) :: norm_after_clear
    integer, intent(inout) :: support_count
    integer :: initialized_tmp
    integer :: zero_tmp_status
    integer :: single_tmp_status
    integer :: compact_tmp_status
    integer :: weight_norm_tmp_status
    integer :: nonnegative_tmp_status
    integer :: component_tmp_status
    integer :: boundary_tmp_status
    integer :: integrated_tmp_status
    integer :: finite_tmp_status
    integer :: norm_finite_tmp_status
    integer :: clear_tmp_status
    integer :: no_rhs_tmp_status
    integer :: no_ibm_tmp_status
    integer :: no_feedback_tmp_status
    integer :: no_twoway_tmp_status
    integer :: no_structure_tmp_status
    integer :: no_fluid_access_tmp_status
    integer :: no_fluid_mod_tmp_status
    integer :: kernel_tmp_status
    integer :: support_tmp_local
    real(mytype) :: zero_tmp
    real(mytype) :: weight_error_tmp_local
    real(mytype) :: min_weight_tmp_local
    real(mytype) :: error_x_tmp_local
    real(mytype) :: error_y_tmp_local
    real(mytype) :: error_z_tmp_local
    real(mytype) :: error_l2_tmp_local
    real(mytype) :: max_force_tmp_local
    real(mytype) :: norm_after_clear_tmp_local

    call stage13_spreading_kernel_get_status_values(initialized_tmp, zero_tmp_status, single_tmp_status, &
                                                    compact_tmp_status, weight_norm_tmp_status, &
                                                    nonnegative_tmp_status, component_tmp_status, &
                                                    boundary_tmp_status, integrated_tmp_status, &
                                                    finite_tmp_status, norm_finite_tmp_status, clear_tmp_status, &
                                                    no_rhs_tmp_status, no_ibm_tmp_status, no_feedback_tmp_status, &
                                                    no_twoway_tmp_status, no_structure_tmp_status, &
                                                    no_fluid_access_tmp_status, no_fluid_mod_tmp_status, &
                                                    kernel_tmp_status)
    call stage13_spreading_kernel_get_diagnostics(zero_tmp, support_tmp_local, weight_error_tmp_local, &
                                                  min_weight_tmp_local, error_x_tmp_local, error_y_tmp_local, &
                                                  error_z_tmp_local, error_l2_tmp_local, max_force_tmp_local, &
                                                  norm_after_clear_tmp_local)

    initialized_status = min(initialized_status, initialized_tmp)
    zero_force_spreading_status = max(zero_force_spreading_status, zero_tmp_status)
    single_point_spreading_status = max(single_point_spreading_status, single_tmp_status)
    compact_support_status = min(compact_support_status, compact_tmp_status)
    weight_normalization_status = min(weight_normalization_status, weight_norm_tmp_status)
    nonnegative_weight_status = min(nonnegative_weight_status, nonnegative_tmp_status)
    component_spreading_status = min(component_spreading_status, component_tmp_status)
    boundary_safe_status = min(boundary_safe_status, boundary_tmp_status)
    integrated_force_conservation_status = min(integrated_force_conservation_status, integrated_tmp_status)
    finite_force_density_status = min(finite_force_density_status, finite_tmp_status)
    force_density_norm_finite_status = min(force_density_norm_finite_status, norm_finite_tmp_status)
    no_rhs_injection_status = min(no_rhs_injection_status, no_rhs_tmp_status)
    no_ibm_spreading_status = min(no_ibm_spreading_status, no_ibm_tmp_status)
    no_feedback_application_status = min(no_feedback_application_status, no_feedback_tmp_status)
    no_twoway_force_status = min(no_twoway_force_status, no_twoway_tmp_status)
    no_structure_advance_status = min(no_structure_advance_status, no_structure_tmp_status)
    no_fluid_field_access_status = min(no_fluid_field_access_status, no_fluid_access_tmp_status)
    no_fluid_field_modification_status = min(no_fluid_field_modification_status, no_fluid_mod_tmp_status)
    if (zero_tmp_status == 1) zero_max = max(zero_max, zero_tmp)
    weight_error = max(weight_error, weight_error_tmp_local)
    weight_min = min(weight_min, min_weight_tmp_local)
    error_x = error_x_tmp_local
    error_y = error_y_tmp_local
    error_z = error_z_tmp_local
    error_l2 = error_l2_tmp_local
    max_force = max(max_force, max_force_tmp_local)
    norm_after_clear = norm_after_clear_tmp_local
    support_count = max(support_count, support_tmp_local)
  end subroutine collect_kernel_status

end program fibre_stage13_spreading_kernel_check
