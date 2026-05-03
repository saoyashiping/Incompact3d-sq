program fibre_stage4_grid_adapter_check

  use fibre_parameters, only : mytype
  use fibre_stage4_config, only : stage4_oneway_config_t, init_stage4_oneway_config
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t, init_stage4_grid_adapter_from_arrays, destroy_stage4_grid_adapter, stage4_adapter_rhs_disabled_flag

  implicit none

  type(stage4_grid_adapter_t) :: adapter_uniform, adapter_nonuniform, adapter_layout0, adapter_layout1, adapter_layout2
  type(stage4_oneway_config_t) :: config

  integer, parameter :: nx = 16, ny = 12, nz = 10
  real(mytype) :: x(nx), y_uniform(ny), y_nonuniform(ny), z(nz)
  real(mytype) :: eta
  integer :: i, unit_id
  integer :: layout_staggered_requires_component_coordinates
  integer :: adapter_stage4_rhs_disabled_flag_value
  integer :: uniform_adapter_to_ibm_grid_possible, nonuniform_adapter_to_ibm_grid_possible
  integer :: stage4_grid_adapter_status

  do i = 1, nx
    x(i) = 0._mytype + (real(i, mytype) - 0.5_mytype) * (2._mytype / real(nx, mytype))
  end do
  do i = 1, ny
    y_uniform(i) = -1._mytype + (real(i, mytype) - 0.5_mytype) * (2._mytype / real(ny, mytype))
    eta = (real(i, mytype) - 0.5_mytype) / real(ny, mytype)
    y_nonuniform(i) = -1._mytype + 2._mytype * eta * eta
  end do
  do i = 1, nz
    z(i) = 0._mytype + (real(i, mytype) - 0.5_mytype) * (1._mytype / real(nz, mytype))
  end do

  call init_stage4_grid_adapter_from_arrays(adapter_uniform, x, y_uniform, z, .true., .false., .true., 1)
  call init_stage4_grid_adapter_from_arrays(adapter_nonuniform, x, y_nonuniform, z, .true., .false., .true., 1)

  call init_stage4_grid_adapter_from_arrays(adapter_layout0, x, y_uniform, z, .true., .false., .true., 0)
  call init_stage4_grid_adapter_from_arrays(adapter_layout1, x, y_uniform, z, .true., .false., .true., 1)
  call init_stage4_grid_adapter_from_arrays(adapter_layout2, x, y_uniform, z, .true., .false., .true., 2)

  layout_staggered_requires_component_coordinates = 1

  call init_stage4_oneway_config(config)
  call stage4_adapter_rhs_disabled_flag(config%apply_ibm_to_fluid_rhs, adapter_stage4_rhs_disabled_flag_value)

  uniform_adapter_to_ibm_grid_possible = merge(1, 0, adapter_uniform%uniform_ibm_compatible)
  nonuniform_adapter_to_ibm_grid_possible = merge(1, 0, adapter_nonuniform%uniform_ibm_compatible)

  if (uniform_adapter_to_ibm_grid_possible == 1 .and. nonuniform_adapter_to_ibm_grid_possible == 0) then
    stage4_grid_adapter_status = 1
  else
    stage4_grid_adapter_status = 0
  end if

  open(newunit=unit_id, file='stage4_outputs/fibre_stage4_grid_adapter_check.dat', status='replace', action='write', form='formatted')

  write(unit_id,'(A,1X,I0)') 'uniform_adapter_nx', adapter_uniform%nx
  write(unit_id,'(A,1X,I0)') 'uniform_adapter_ny', adapter_uniform%ny
  write(unit_id,'(A,1X,I0)') 'uniform_adapter_nz', adapter_uniform%nz
  write(unit_id,'(A,1X,I0)') 'uniform_adapter_periodic_x', merge(1, 0, adapter_uniform%periodic_x)
  write(unit_id,'(A,1X,I0)') 'uniform_adapter_periodic_y', merge(1, 0, adapter_uniform%periodic_y)
  write(unit_id,'(A,1X,I0)') 'uniform_adapter_periodic_z', merge(1, 0, adapter_uniform%periodic_z)
  write(unit_id,'(A,1X,I0)') 'uniform_adapter_uniform_x', merge(1, 0, adapter_uniform%uniform_x)
  write(unit_id,'(A,1X,I0)') 'uniform_adapter_uniform_y', merge(1, 0, adapter_uniform%uniform_y)
  write(unit_id,'(A,1X,I0)') 'uniform_adapter_uniform_z', merge(1, 0, adapter_uniform%uniform_z)
  write(unit_id,'(A,1X,I0)') 'uniform_adapter_uniform_ibm_compatible', merge(1, 0, adapter_uniform%uniform_ibm_compatible)
  write(unit_id,'(A,1X,I0)') 'uniform_adapter_velocity_layout_mode', adapter_uniform%velocity_layout_mode

  write(unit_id,'(A,1X,I0)') 'nonuniform_y_adapter_uniform_x', merge(1, 0, adapter_nonuniform%uniform_x)
  write(unit_id,'(A,1X,I0)') 'nonuniform_y_adapter_uniform_y', merge(1, 0, adapter_nonuniform%uniform_y)
  write(unit_id,'(A,1X,I0)') 'nonuniform_y_adapter_uniform_z', merge(1, 0, adapter_nonuniform%uniform_z)
  write(unit_id,'(A,1X,I0)') 'nonuniform_y_adapter_uniform_ibm_compatible', merge(1, 0, adapter_nonuniform%uniform_ibm_compatible)
  write(unit_id,'(A,1X,ES24.16E3)') 'nonuniform_y_adapter_dy_min', adapter_nonuniform%dy_min
  write(unit_id,'(A,1X,ES24.16E3)') 'nonuniform_y_adapter_dy_max', adapter_nonuniform%dy_max

  write(unit_id,'(A,1X,I0)') 'layout_unknown_mode', adapter_layout0%velocity_layout_mode
  write(unit_id,'(A,1X,I0)') 'layout_collocated_mode', adapter_layout1%velocity_layout_mode
  write(unit_id,'(A,1X,I0)') 'layout_staggered_mode', adapter_layout2%velocity_layout_mode
  write(unit_id,'(A,1X,I0)') 'layout_staggered_requires_component_coordinates', layout_staggered_requires_component_coordinates

  write(unit_id,'(A,1X,I0)') 'adapter_stage4_apply_ibm_to_fluid_rhs', merge(1, 0, config%apply_ibm_to_fluid_rhs)
  write(unit_id,'(A,1X,I0)') 'adapter_stage4_rhs_disabled_flag', adapter_stage4_rhs_disabled_flag_value

  write(unit_id,'(A,1X,I0)') 'uniform_adapter_to_ibm_grid_possible', uniform_adapter_to_ibm_grid_possible
  write(unit_id,'(A,1X,I0)') 'nonuniform_adapter_to_ibm_grid_possible', nonuniform_adapter_to_ibm_grid_possible

  write(unit_id,'(A,1X,I0)') 'stage4_grid_adapter_status', stage4_grid_adapter_status
  close(unit_id)

  call destroy_stage4_grid_adapter(adapter_uniform)
  call destroy_stage4_grid_adapter(adapter_nonuniform)
  call destroy_stage4_grid_adapter(adapter_layout0)
  call destroy_stage4_grid_adapter(adapter_layout1)
  call destroy_stage4_grid_adapter(adapter_layout2)

end program fibre_stage4_grid_adapter_check
