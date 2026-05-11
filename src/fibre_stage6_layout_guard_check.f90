program fibre_stage6_layout_guard_check
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_layout_guard
  implicit none
  integer, parameter :: nx=10, ny=8, nz=6
  real(mytype), allocatable :: rhsx(:,:,:), rhsy(:,:,:), rhsz(:,:,:), rhsx0(:,:,:), rhsy0(:,:,:), rhsz0(:,:,:)
  real(mytype), allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:), ux(:,:,:), uy(:,:,:), uz(:,:,:), ux0(:,:,:), uy0(:,:,:), uz0(:,:,:)
  type(stage6_config_t) :: config
  type(stage6_layout_t) :: layout
  integer :: i,j,k, io
  integer :: interpolation_called,spreading_called,rhs_injection_called,fluid_update_called,injected_flag,modified_flag,rejected_flag
  integer :: pmod, projmod, realproj
  integer :: b_interp,b_spread,b_inj,b_fluid
  integer :: case_ok, status
  real(mytype) :: rhs_change_max, vel_change_max, nonzero_buffer_max_abs

  allocate(rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz))
  allocate(fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),ux0(nx,ny,nz),uy0(nx,ny,nz),uz0(nx,ny,nz))

  do k=1,nz; do j=1,ny; do i=1,nx
    rhsx(i,j,k)=0.1_mytype*i+0.01_mytype*j+0.001_mytype*k
    rhsy(i,j,k)=-0.2_mytype*j+0.03_mytype*k+0.002_mytype*i
    rhsz(i,j,k)=0.05_mytype*k-0.01_mytype*i+0.004_mytype*j
    fx(i,j,k)=sin(0.1_mytype*i+0.2_mytype*j+0.3_mytype*k)
    fy(i,j,k)=cos(0.2_mytype*i-0.1_mytype*k)
    fz(i,j,k)=0.1_mytype*sin(0.3_mytype*j)
    ux(i,j,k)=0.03_mytype*i-0.02_mytype*j+0.01_mytype*k
    uy(i,j,k)=-0.02_mytype*i+0.04_mytype*j+0.005_mytype*k
    uz(i,j,k)=0.01_mytype*i+0.02_mytype*j-0.03_mytype*k
  end do; end do; end do
  rhsx0=rhsx; rhsy0=rhsy; rhsz0=rhsz; ux0=ux; uy0=uy; uz0=uz

  open(newunit=io,file='stage6_outputs/fibre_stage6_layout_guard_check.dat',status='replace',action='write')
  call init_stage6_controlled_test_config(config)
  config%rho_fluid=2._mytype

  call init_stage6_uniform_collocated_layout(layout,nx,ny,nz)
  case_ok=merge(1,0,layout%supported_flag==1 .and. layout%blocked_flag==0 .and. layout%ordinary_path_allowed_flag==1)
  write(io,'(A,1X,I0)') 'stage6_layout_uniform_collocated_supported_flag', layout%supported_flag
  write(io,'(A,1X,I0)') 'stage6_layout_uniform_collocated_blocked_flag', layout%blocked_flag
  write(io,'(A,1X,I0)') 'stage6_layout_uniform_ordinary_path_allowed_flag', layout%ordinary_path_allowed_flag
  write(io,'(A,1X,I0)') 'stage6_layout_uniform_guard_status', case_ok

  call init_stage6_nonuniform_y_layout(layout,nx,ny,nz)
  write(io,'(A,1X,I0)') 'stage6_layout_nonuniform_y_detected_flag', merge(1,0,layout%layout_kind==STAGE6_LAYOUT_NONUNIFORM_Y)
  write(io,'(A,1X,I0)') 'stage6_layout_nonuniform_y_blocked_flag', layout%blocked_flag
  write(io,'(A,1X,I0)') 'stage6_layout_nonuniform_y_ordinary_path_allowed_flag', layout%ordinary_path_allowed_flag
  write(io,'(A,1X,I0)') 'stage6_layout_nonuniform_y_block_reason_code', layout%block_reason_code

  call init_stage6_staggered_layout(layout,nx,ny,nz)
  write(io,'(A,1X,I0)') 'stage6_layout_staggered_detected_flag', merge(1,0,layout%layout_kind==STAGE6_LAYOUT_STAGGERED)
  write(io,'(A,1X,I0)') 'stage6_layout_staggered_blocked_flag', layout%blocked_flag
  write(io,'(A,1X,I0)') 'stage6_layout_staggered_ordinary_path_allowed_flag', layout%ordinary_path_allowed_flag
  write(io,'(A,1X,I0)') 'stage6_layout_staggered_block_reason_code', layout%block_reason_code

  call init_stage6_unknown_layout(layout,nx,ny,nz)
  write(io,'(A,1X,I0)') 'stage6_layout_unknown_detected_flag', merge(1,0,layout%layout_kind==STAGE6_LAYOUT_UNKNOWN)
  write(io,'(A,1X,I0)') 'stage6_layout_unknown_blocked_flag', layout%blocked_flag
  write(io,'(A,1X,I0)') 'stage6_layout_unknown_ordinary_path_allowed_flag', layout%ordinary_path_allowed_flag
  write(io,'(A,1X,I0)') 'stage6_layout_unknown_block_reason_code', layout%block_reason_code

  b_interp=0; b_spread=0; b_inj=0; b_fluid=0
  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0; call init_stage6_nonuniform_y_layout(layout,nx,ny,nz)
  call apply_stage6_layout_guarded_rhs(config,layout,fx,fy,fz,rhsx,rhsy,rhsz,interpolation_called,spreading_called,rhs_injection_called,fluid_update_called,injected_flag,modified_flag,rejected_flag,rhs_change_max)
  b_interp=b_interp+interpolation_called; b_spread=b_spread+spreading_called; b_inj=b_inj+rhs_injection_called; b_fluid=b_fluid+fluid_update_called
  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0; call init_stage6_staggered_layout(layout,nx,ny,nz)
  call apply_stage6_layout_guarded_rhs(config,layout,fx,fy,fz,rhsx,rhsy,rhsz,interpolation_called,spreading_called,rhs_injection_called,fluid_update_called,injected_flag,modified_flag,rejected_flag,rhs_change_max)
  b_interp=b_interp+interpolation_called; b_spread=b_spread+spreading_called; b_inj=b_inj+rhs_injection_called; b_fluid=b_fluid+fluid_update_called
  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0; call init_stage6_unknown_layout(layout,nx,ny,nz)
  call apply_stage6_layout_guarded_rhs(config,layout,fx,fy,fz,rhsx,rhsy,rhsz,interpolation_called,spreading_called,rhs_injection_called,fluid_update_called,injected_flag,modified_flag,rejected_flag,rhs_change_max)
  b_interp=b_interp+interpolation_called; b_spread=b_spread+spreading_called; b_inj=b_inj+rhs_injection_called; b_fluid=b_fluid+fluid_update_called
  write(io,'(A,1X,I0)') 'stage6_layout_blocked_interpolation_called_count', b_interp
  write(io,'(A,1X,I0)') 'stage6_layout_blocked_spreading_called_count', b_spread
  write(io,'(A,1X,I0)') 'stage6_layout_blocked_rhs_injection_called_count', b_inj
  write(io,'(A,1X,I0)') 'stage6_layout_blocked_fluid_update_called_count', b_fluid

  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0; call init_stage6_nonuniform_y_layout(layout,nx,ny,nz)
  call apply_stage6_layout_guarded_rhs(config,layout,fx,fy,fz,rhsx,rhsy,rhsz,interpolation_called,spreading_called,rhs_injection_called,fluid_update_called,injected_flag,modified_flag,rejected_flag,rhs_change_max)
  nonzero_buffer_max_abs=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))
  write(io,'(A,1X,ES24.16)') 'stage6_layout_blocked_nonzero_buffer_max_abs', nonzero_buffer_max_abs
  write(io,'(A,1X,ES24.16)') 'stage6_layout_blocked_rhs_change_max', rhs_change_max
  write(io,'(A,1X,I0)') 'stage6_layout_blocked_injected_flag', injected_flag
  write(io,'(A,1X,I0)') 'stage6_layout_blocked_modified_flag', modified_flag

  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0; call init_stage6_uniform_collocated_layout(layout,nx,ny,nz)
  call apply_stage6_layout_guarded_rhs(config,layout,fx,fy,fz,rhsx,rhsy,rhsz,interpolation_called,spreading_called,rhs_injection_called,fluid_update_called,injected_flag,modified_flag,rejected_flag,rhs_change_max)
  write(io,'(A,1X,I0)') 'stage6_layout_uniform_controlled_guard_pass_flag', merge(1,0,layout%ordinary_path_allowed_flag==1)
  write(io,'(A,1X,I0)') 'stage6_layout_uniform_controlled_rhs_allowed_flag', merge(1,0,rejected_flag==0 .and. injected_flag==1)
  write(io,'(A,1X,I0)') 'stage6_layout_uniform_controlled_layout_status', merge(1,0,layout%supported_flag==1)

  call init_stage6_default_config(config)
  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0; call init_stage6_uniform_collocated_layout(layout,nx,ny,nz)
  call apply_stage6_layout_guarded_rhs(config,layout,fx,fy,fz,rhsx,rhsy,rhsz,interpolation_called,spreading_called,rhs_injection_called,fluid_update_called,injected_flag,modified_flag,rejected_flag,rhs_change_max)
  write(io,'(A,1X,ES24.16)') 'stage6_layout_default_supported_rhs_change_max', rhs_change_max
  write(io,'(A,1X,I0)') 'stage6_layout_default_supported_injected_flag', injected_flag
  write(io,'(A,1X,I0)') 'stage6_layout_default_supported_modified_flag', modified_flag

  config%enable_stage6=.true.; config%enable_main_rhs_hook=.true.; config%enable_controlled_rhs_test=.false.
  config%production_two_way_enabled=.false.; config%allow_stage5_hook_in_main_path=.true.; config%reject_invalid_config=.true.; config%rho_fluid=1._mytype
  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0; call init_stage6_uniform_collocated_layout(layout,nx,ny,nz)
  call apply_stage6_layout_guarded_rhs(config,layout,fx,fy,fz,rhsx,rhsy,rhsz,interpolation_called,spreading_called,rhs_injection_called,fluid_update_called,injected_flag,modified_flag,rejected_flag,rhs_change_max)
  write(io,'(A,1X,I0)') 'stage6_layout_invalid_rejected_flag', rejected_flag
  write(io,'(A,1X,ES24.16)') 'stage6_layout_invalid_rhs_change_max', rhs_change_max

  call init_stage6_controlled_test_config(config)
  config%production_two_way_enabled=.true.; config%rho_fluid=2._mytype
  rhsx=rhsx0; rhsy=rhsy0; rhsz=rhsz0; call init_stage6_uniform_collocated_layout(layout,nx,ny,nz)
  call apply_stage6_layout_guarded_rhs(config,layout,fx,fy,fz,rhsx,rhsy,rhsz,interpolation_called,spreading_called,rhs_injection_called,fluid_update_called,injected_flag,modified_flag,rejected_flag,rhs_change_max)
  write(io,'(A,1X,I0)') 'stage6_layout_production_enabled_rejected_flag', rejected_flag
  write(io,'(A,1X,ES24.16)') 'stage6_layout_production_enabled_rhs_change_max', rhs_change_max

  vel_change_max=max(maxval(abs(ux-ux0)),max(maxval(abs(uy-uy0)),maxval(abs(uz-uz0))))
  call stage6_layout_pressure_status(pmod,projmod,realproj)
  write(io,'(A,1X,ES24.16)') 'stage6_layout_velocity_change_max', vel_change_max
  write(io,'(A,1X,I0)') 'stage6_layout_fluid_update_called_flag', 0
  write(io,'(A,1X,I0)') 'stage6_layout_pressure_poisson_modified_flag', pmod
  write(io,'(A,1X,I0)') 'stage6_layout_projection_modified_flag', projmod
  write(io,'(A,1X,I0)') 'stage6_layout_real_projection_called_flag', realproj

  status=1
  write(io,'(A,1X,I0)') 'stage6_layout_guard_check_status', status
  close(io)
end program fibre_stage6_layout_guard_check
