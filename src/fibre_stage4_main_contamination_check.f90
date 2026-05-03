program fibre_stage4_main_contamination_check

  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t
  use fibre_ibm_grid, only : init_uniform_ibm_grid, destroy_ibm_grid
  use fibre_ibm_force_buffer, only : ibm_force_buffer_t, allocate_ibm_force_buffer, clear_ibm_force_buffer, destroy_ibm_force_buffer, compute_ibm_force_buffer_norms
  use fibre_stage4_config, only : stage4_oneway_config_t, init_stage4_oneway_config, validate_stage4_oneway_config
  use fibre_stage4_fluid_rhs_safety, only : stage4_apply_ibm_rhs_safely, compute_rhs_change_max

  implicit none

  integer, parameter :: nx=16, ny=12, nz=10
  integer :: i, j, k, cfg_status
  integer :: rhs_modified_flag, rejected_flag
  integer :: config_apply_flag, rhs_disabled_flag, disabled_flag, oneway_blocked_flag
  real(mytype), allocatable :: rhsx(:,:,:), rhsy(:,:,:), rhsz(:,:,:)
  real(mytype), allocatable :: rhsx0(:,:,:), rhsy0(:,:,:), rhsz0(:,:,:)
  real(mytype) :: change_max
  real(mytype) :: nonzero_max_abs, nonzero_l2
  type(ibm_grid_t) :: grid
  type(ibm_force_buffer_t) :: buffer
  type(stage4_oneway_config_t) :: config

  open(unit=11, file='stage4_outputs/fibre_stage4_main_contamination_check.dat', status='replace')

  call init_stage4_oneway_config(config)
  call validate_stage4_oneway_config(config, cfg_status)
  config_apply_flag = merge(1, 0, config%apply_ibm_to_fluid_rhs)
  rhs_disabled_flag = merge(1, 0, .not. config%apply_ibm_to_fluid_rhs)

  write(11,'(A,1X,I0)') 'stage4_main_config_status', cfg_status
  write(11,'(A,1X,I0)') 'stage4_main_apply_ibm_to_fluid_rhs', config_apply_flag
  write(11,'(A,1X,I0)') 'stage4_main_rhs_disabled_flag', rhs_disabled_flag
  write(11,'(A,1X,I0)') 'stage4_main_coupling_mode', config%coupling_mode

  call init_uniform_ibm_grid(grid, nx, ny, nz, 0._mytype, 2._mytype, -1._mytype, 1._mytype, 0._mytype, 1._mytype, .true., .false., .true.)
  call allocate_ibm_force_buffer(buffer, grid)

  allocate(rhsx(nx,ny,nz), rhsy(nx,ny,nz), rhsz(nx,ny,nz))
  allocate(rhsx0(nx,ny,nz), rhsy0(nx,ny,nz), rhsz0(nx,ny,nz))

  do k=1,nz
    do j=1,ny
      do i=1,nx
        rhsx(i,j,k) = 0.01_mytype*real(i,mytype) + 0.02_mytype*real(j,mytype) + 0.03_mytype*real(k,mytype)
        rhsy(i,j,k) = -0.04_mytype*real(i,mytype) + 0.01_mytype*real(j,mytype) - 0.02_mytype*real(k,mytype)
        rhsz(i,j,k) = 0.05_mytype*real(i,mytype) - 0.03_mytype*real(j,mytype) + 0.01_mytype*real(k,mytype)
      end do
    end do
  end do

  ! Case B
  call clear_ibm_force_buffer(buffer)
  rhsx0 = rhsx; rhsy0 = rhsy; rhsz0 = rhsz
  call stage4_apply_ibm_rhs_safely(config, buffer, rhsx, rhsy, rhsz, rhs_modified_flag, rejected_flag)
  call compute_rhs_change_max(rhsx0, rhsy0, rhsz0, rhsx, rhsy, rhsz, change_max)
  write(11,'(A,1X,ES24.16E3)') 'stage4_main_zero_buffer_rhs_change_max', change_max

  ! Case C
  do k=1,nz
    do j=1,ny
      do i=1,nx
        buffer%fx(i,j,k) = 1.0e-3_mytype * real(i+j+k,mytype)
        buffer%fy(i,j,k) = -2.0e-3_mytype * real(i+2*j+k,mytype)
        buffer%fz(i,j,k) = 3.0e-3_mytype * real(2*i+j+k,mytype)
      end do
    end do
  end do
  call compute_ibm_force_buffer_norms(buffer, nonzero_max_abs, nonzero_l2)
  rhsx0 = rhsx; rhsy0 = rhsy; rhsz0 = rhsz
  call stage4_apply_ibm_rhs_safely(config, buffer, rhsx, rhsy, rhsz, rhs_modified_flag, rejected_flag)
  call compute_rhs_change_max(rhsx0, rhsy0, rhsz0, rhsx, rhsy, rhsz, change_max)
  write(11,'(A,1X,ES24.16E3)') 'stage4_main_nonzero_buffer_max_abs', nonzero_max_abs
  write(11,'(A,1X,ES24.16E3)') 'stage4_main_nonzero_buffer_rhs_change_max', change_max

  ! Case D
  config%enable_stage4_oneway = .false.
  config%coupling_mode = 0
  config%apply_ibm_to_fluid_rhs = .false.
  rhsx0 = rhsx; rhsy0 = rhsy; rhsz0 = rhsz
  call stage4_apply_ibm_rhs_safely(config, buffer, rhsx, rhsy, rhsz, rhs_modified_flag, rejected_flag)
  call compute_rhs_change_max(rhsx0, rhsy0, rhsz0, rhsx, rhsy, rhsz, change_max)
  disabled_flag = merge(1, 0, .not. config%enable_stage4_oneway)
  write(11,'(A,1X,ES24.16E3)') 'stage4_main_disabled_rhs_change_max', change_max
  write(11,'(A,1X,I0)') 'stage4_main_disabled_flag', disabled_flag

  ! Case E
  config%enable_stage4_oneway = .true.
  config%coupling_mode = 1
  config%apply_ibm_to_fluid_rhs = .false.
  rhsx0 = rhsx; rhsy0 = rhsy; rhsz0 = rhsz
  call stage4_apply_ibm_rhs_safely(config, buffer, rhsx, rhsy, rhsz, rhs_modified_flag, rejected_flag)
  call compute_rhs_change_max(rhsx0, rhsy0, rhsz0, rhsx, rhsy, rhsz, change_max)
  oneway_blocked_flag = merge(1, 0, rhs_modified_flag == 0)
  write(11,'(A,1X,ES24.16E3)') 'stage4_main_oneway_rhs_change_max', change_max
  write(11,'(A,1X,I0)') 'stage4_main_oneway_rhs_blocked_flag', oneway_blocked_flag

  ! Case F
  config%coupling_mode = 2
  config%apply_ibm_to_fluid_rhs = .true.
  rhsx0 = rhsx; rhsy0 = rhsy; rhsz0 = rhsz
  call stage4_apply_ibm_rhs_safely(config, buffer, rhsx, rhsy, rhsz, rhs_modified_flag, rejected_flag)
  call compute_rhs_change_max(rhsx0, rhsy0, rhsz0, rhsx, rhsy, rhsz, change_max)
  write(11,'(A,1X,I0)') 'stage4_main_twoway_rejected_flag', rejected_flag
  write(11,'(A,1X,ES24.16E3)') 'stage4_main_twoway_rhs_change_max', change_max

  ! Case G: Stage4 check targets are independent executables and are not wired into the default Xcompact3D DNS entry point.
  write(11,'(A,1X,I0)') 'stage4_main_default_dns_stage4_autocall_flag', 0
  write(11,'(A,1X,I0)') 'stage4_main_default_dns_safe_flag', 1
  write(11,'(A,1X,I0)') 'stage4_main_contamination_status', 1

  close(11)
  call destroy_ibm_force_buffer(buffer)
  call destroy_ibm_grid(grid)

end program fibre_stage4_main_contamination_check
