program fibre_stage11_grid_metadata_check
  use decomp_2d_constants, only : mytype
  use fibre_stage11_config, only : stage11_config_load, stage11_requested, stage11_readonly_mode, stage11_get_max_points
  use fibre_stage11_lagrangian_state, only : stage11_lagrangian_state_init, stage11_lagrangian_state_finalize, &
                                             stage11_lagrangian_state_is_allocated
  use fibre_stage11_grid_metadata, only : stage11_grid_metadata_init_uniform, stage11_grid_metadata_finalize, &
                                          stage11_grid_metadata_is_initialized, stage11_grid_metadata_point_in_domain, &
                                          stage11_grid_metadata_get_status_values, stage11_grid_metadata_write_diagnostics
  implicit none

  integer :: requested_flag, readonly_mode_status, lagrangian_state_available_status
  integer :: grid_initialized_status, global_size_status, local_bounds_status
  integer :: extents_finite_status, spacing_positive_status, periodic_flags_status
  integer :: staggered_layout_status, nonuniform_y_policy_status, domain_safety_status
  integer :: no_fluid_field_access_status, no_velocity_sampling_status, no_fluid_field_modification_status
  integer :: no_rhs_injection_status, no_ibm_spreading_status, no_feedback_force_status
  integer :: no_twoway_force_status, no_structure_advance_status, grid_metadata_status
  integer :: n_points, k, io_unit, ios
  logical :: pass
  real(mytype) :: xk, yk, zk, denom

  call execute_command_line('mkdir -p stage11_outputs', exitstat=ios)
  call stage11_config_load()

  requested_flag = 0
  if (stage11_requested()) requested_flag = 1
  readonly_mode_status = 0
  if (stage11_readonly_mode()) readonly_mode_status = 1

  n_points = stage11_get_max_points()
  if (n_points <= 0) n_points = 8

  call stage11_lagrangian_state_init(n_points)
  lagrangian_state_available_status = 0
  if (stage11_lagrangian_state_is_allocated()) lagrangian_state_available_status = 1

  call stage11_grid_metadata_init_uniform()

  call stage11_grid_metadata_get_status_values(grid_initialized_status, global_size_status, local_bounds_status, &
                                               extents_finite_status, spacing_positive_status, periodic_flags_status, &
                                               staggered_layout_status, nonuniform_y_policy_status, domain_safety_status, &
                                               no_fluid_field_access_status, no_velocity_sampling_status, &
                                               no_fluid_field_modification_status, no_rhs_injection_status, &
                                               no_ibm_spreading_status, no_feedback_force_status, no_twoway_force_status, &
                                               no_structure_advance_status, grid_metadata_status)

  if (domain_safety_status == 1) then
    denom = real(max(1, n_points - 1), mytype)
    do k = 1, n_points
      xk = real(k - 1, mytype) / denom
      yk = 0.5_mytype
      zk = 0.5_mytype
      if (.not. stage11_grid_metadata_point_in_domain(xk, yk, zk)) then
        domain_safety_status = 0
      end if
    end do
  end if

  if (domain_safety_status == 0) grid_metadata_status = 0

  open(newunit=io_unit,file='stage11_outputs/fibre_stage11_2_grid_metadata.dat',status='replace',action='write',iostat=ios)
  if (ios /= 0) then
    print *, 'STAGE 11.2 GRID METADATA VERDICT: FAIL'
    print *, 'Reason: unable to open stage11_outputs/fibre_stage11_2_grid_metadata.dat'
    stop 1
  end if

  write(io_unit,'(A,1X,I0)') 'stage11_2_requested_flag', requested_flag
  write(io_unit,'(A,1X,I0)') 'stage11_2_readonly_mode_status', readonly_mode_status
  write(io_unit,'(A,1X,I0)') 'stage11_2_lagrangian_state_available_status', lagrangian_state_available_status
  call stage11_grid_metadata_write_diagnostics(io_unit)
  close(io_unit)

  pass = .true.
  if (requested_flag /= 1) then; print *, 'Reason: stage11_2_requested_flag /= 1'; pass=.false.; end if
  if (readonly_mode_status /= 1) then; print *, 'Reason: stage11_2_readonly_mode_status /= 1'; pass=.false.; end if
  if (lagrangian_state_available_status /= 1) then; print *, 'Reason: stage11_2_lagrangian_state_available_status /= 1'; pass=.false.; end if
  if (grid_initialized_status /= 1) then; print *, 'Reason: stage11_2_grid_initialized_status /= 1'; pass=.false.; end if
  if (global_size_status /= 1) then; print *, 'Reason: stage11_2_global_size_status /= 1'; pass=.false.; end if
  if (local_bounds_status /= 1) then; print *, 'Reason: stage11_2_local_bounds_status /= 1'; pass=.false.; end if
  if (extents_finite_status /= 1) then; print *, 'Reason: stage11_2_extents_finite_status /= 1'; pass=.false.; end if
  if (spacing_positive_status /= 1) then; print *, 'Reason: stage11_2_spacing_positive_status /= 1'; pass=.false.; end if
  if (periodic_flags_status /= 1) then; print *, 'Reason: stage11_2_periodic_flags_status /= 1'; pass=.false.; end if
  if (staggered_layout_status /= 1) then; print *, 'Reason: stage11_2_staggered_layout_status /= 1'; pass=.false.; end if
  if (nonuniform_y_policy_status /= 1) then; print *, 'Reason: stage11_2_nonuniform_y_policy_status /= 1'; pass=.false.; end if
  if (domain_safety_status /= 1) then; print *, 'Reason: stage11_2_domain_safety_status /= 1'; pass=.false.; end if
  if (no_fluid_field_access_status /= 1) then; print *, 'Reason: stage11_2_no_fluid_field_access_status /= 1'; pass=.false.; end if
  if (no_velocity_sampling_status /= 1) then; print *, 'Reason: stage11_2_no_velocity_sampling_status /= 1'; pass=.false.; end if
  if (no_fluid_field_modification_status /= 1) then; print *, 'Reason: stage11_2_no_fluid_field_modification_status /= 1'; pass=.false.; end if
  if (no_rhs_injection_status /= 1) then; print *, 'Reason: stage11_2_no_rhs_injection_status /= 1'; pass=.false.; end if
  if (no_ibm_spreading_status /= 1) then; print *, 'Reason: stage11_2_no_ibm_spreading_status /= 1'; pass=.false.; end if
  if (no_feedback_force_status /= 1) then; print *, 'Reason: stage11_2_no_feedback_force_status /= 1'; pass=.false.; end if
  if (no_twoway_force_status /= 1) then; print *, 'Reason: stage11_2_no_twoway_force_status /= 1'; pass=.false.; end if
  if (no_structure_advance_status /= 1) then; print *, 'Reason: stage11_2_no_structure_advance_status /= 1'; pass=.false.; end if
  if (grid_metadata_status /= 1) then; print *, 'Reason: stage11_2_grid_metadata_status /= 1'; pass=.false.; end if

  if (pass) then
    print *, 'STAGE 11.2 GRID METADATA VERDICT: PASS'
  else
    print *, 'STAGE 11.2 GRID METADATA VERDICT: FAIL'
    stop 1
  end if

  call stage11_grid_metadata_finalize()
  call stage11_lagrangian_state_finalize()
end program fibre_stage11_grid_metadata_check
