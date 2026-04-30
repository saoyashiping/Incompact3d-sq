program fibre_stage4_oneway_skeleton_check

  use fibre_types, only : fibre_t
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_uniform_ibm_grid, init_lagrangian_points_from_fibre
  use fibre_ibm_force_buffer, only : ibm_force_buffer_t, allocate_ibm_force_buffer, clear_ibm_force_buffer, compute_ibm_force_buffer_norms
  use fibre_stage4_config, only : stage4_oneway_config_t, init_stage4_oneway_config, validate_stage4_oneway_config
  use fibre_stage4_diagnostics, only : stage4_safety_flags

  implicit none

  type(stage4_oneway_config_t) :: config
  type(ibm_grid_t) :: grid
  type(fibre_t) :: fibre
  type(ibm_lagrangian_points_t) :: lag
  type(ibm_force_buffer_t) :: force_buffer

  integer :: config_status, rhs_disabled_flag, diagnostics_enabled_flag
  integer :: i, unit_id
  integer :: interpolation_called, feedback_called, structure_advance_called
  integer :: spreading_called, fluid_rhs_modified, skeleton_status
  real(mytype) :: lag_weight_sum, expected_fibre_length
  real(mytype) :: force_buffer_initial_max_abs, force_buffer_initial_l2_norm

  call init_stage4_oneway_config(config)
  call validate_stage4_oneway_config(config, config_status)
  call stage4_safety_flags(config, rhs_disabled_flag, diagnostics_enabled_flag)

  call init_uniform_ibm_grid(grid, 16, 12, 10, 0._mytype, 2._mytype, -1._mytype, 1._mytype, 0._mytype, 1._mytype, .true., .false., .true.)

  call fibre_init_straight_free_free(fibre, 33, 1._mytype, 1._mytype, 1._mytype)
  do i = 1, fibre%nl
    fibre%x(1, i) = 0.5_mytype + (real(i - 1, mytype) / real(fibre%nl - 1, mytype))
    fibre%x(2, i) = 0._mytype
    fibre%x(3, i) = 0.5_mytype
  end do

  call init_lagrangian_points_from_fibre(lag, fibre)

  call allocate_ibm_force_buffer(force_buffer, grid)
  call clear_ibm_force_buffer(force_buffer)
  call compute_ibm_force_buffer_norms(force_buffer, force_buffer_initial_max_abs, force_buffer_initial_l2_norm)

  lag_weight_sum = sum(lag%weight)
  expected_fibre_length = fibre%length

  interpolation_called = 0
  feedback_called = 0
  structure_advance_called = 0
  spreading_called = 0
  fluid_rhs_modified = 0

  if (config_status == 1 .and. rhs_disabled_flag == 1 .and. diagnostics_enabled_flag == 1) then
    skeleton_status = 1
  else
    skeleton_status = 0
  end if

  open(newunit=unit_id, file='stage4_outputs/fibre_stage4_oneway_skeleton_check.dat', status='replace', action='write', form='formatted')
  write(unit_id,'(A,1X,I0)') 'stage4_config_status', config_status
  write(unit_id,'(A,1X,I0)') 'stage4_enable_oneway', merge(1, 0, config%enable_stage4_oneway)
  write(unit_id,'(A,1X,I0)') 'stage4_enable_ibm_diagnostics', merge(1, 0, config%enable_ibm_diagnostics)
  write(unit_id,'(A,1X,I0)') 'stage4_apply_ibm_to_fluid_rhs', merge(1, 0, config%apply_ibm_to_fluid_rhs)
  write(unit_id,'(A,1X,I0)') 'stage4_rhs_disabled_flag', rhs_disabled_flag
  write(unit_id,'(A,1X,I0)') 'stage4_coupling_mode', config%coupling_mode
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_beta_drag', config%beta_drag
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_structure_dt', config%structure_dt

  write(unit_id,'(A,1X,I0)') 'stage4_grid_nx', grid%nx
  write(unit_id,'(A,1X,I0)') 'stage4_grid_ny', grid%ny
  write(unit_id,'(A,1X,I0)') 'stage4_grid_nz', grid%nz
  write(unit_id,'(A,1X,I0)') 'stage4_grid_periodic_x', merge(1, 0, grid%periodic_x)
  write(unit_id,'(A,1X,I0)') 'stage4_grid_periodic_y', merge(1, 0, grid%periodic_y)
  write(unit_id,'(A,1X,I0)') 'stage4_grid_periodic_z', merge(1, 0, grid%periodic_z)

  write(unit_id,'(A,1X,I0)') 'stage4_fibre_nl', fibre%nl
  write(unit_id,'(A,1X,I0)') 'stage4_lag_nl', lag%nl
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_lag_weight_sum', lag_weight_sum
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_expected_fibre_length', expected_fibre_length

  write(unit_id,'(A,1X,I0)') 'stage4_force_buffer_allocated', merge(1, 0, force_buffer%is_allocated)
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_force_buffer_initial_max_abs', force_buffer_initial_max_abs
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_force_buffer_initial_l2_norm', force_buffer_initial_l2_norm

  write(unit_id,'(A,1X,I0)') 'stage4_interpolation_called', interpolation_called
  write(unit_id,'(A,1X,I0)') 'stage4_feedback_called', feedback_called
  write(unit_id,'(A,1X,I0)') 'stage4_structure_advance_called', structure_advance_called
  write(unit_id,'(A,1X,I0)') 'stage4_spreading_called', spreading_called
  write(unit_id,'(A,1X,I0)') 'stage4_fluid_rhs_modified', fluid_rhs_modified

  write(unit_id,'(A,1X,I0)') 'stage4_skeleton_status', skeleton_status
  close(unit_id)

end program fibre_stage4_oneway_skeleton_check
