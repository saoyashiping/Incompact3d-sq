program fibre_stage6_micro_smoke_check
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_layout_guard
  use fibre_stage6_micro_smoke
  implicit none
  integer, parameter :: nx=10, ny=8, nz=6
  real(mytype), parameter :: dt = 1.0e-5_mytype
  real(mytype), allocatable :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:)
  real(mytype), allocatable :: fx(:,:,:),fy(:,:,:),fz(:,:,:),ux(:,:,:),uy(:,:,:),uz(:,:,:),uxs(:,:,:),uys(:,:,:),uzs(:,:,:),uxs0(:,:,:),uys0(:,:,:),uzs0(:,:,:)
  type(stage6_config_t) :: config
  type(stage6_layout_t) :: layout
  integer :: i,j,k,io,cv,lg,hk,ric,inj,mod,rej
  integer :: ppm,prm,pm,rpc,pvm,pdc,pebd
  real(mytype) :: re,ex,ey,ez,ue,rhschg,bufmax,ustarchg,tmp

  allocate(rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz))
  allocate(fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),uxs(nx,ny,nz),uys(nx,ny,nz),uzs(nx,ny,nz),uxs0(nx,ny,nz),uys0(nx,ny,nz),uzs0(nx,ny,nz))

  do k=1,nz; do j=1,ny; do i=1,nx
    rhsx(i,j,k)=0.1_mytype*i+0.01_mytype*j+0.001_mytype*k
    rhsy(i,j,k)=-0.2_mytype*j+0.03_mytype*k+0.002_mytype*i
    rhsz(i,j,k)=0.05_mytype*k-0.01_mytype*i+0.004_mytype*j
    ux(i,j,k)=0.2_mytype+0.01_mytype*i+0.002_mytype*j
    uy(i,j,k)=0.01_mytype*sin(0.2_mytype*i+0.1_mytype*k)
    uz(i,j,k)=0.01_mytype*cos(0.1_mytype*j+0.3_mytype*k)
    fx(i,j,k)=sin(0.1_mytype*i+0.2_mytype*j+0.3_mytype*k)
    fy(i,j,k)=cos(0.2_mytype*i-0.1_mytype*k)
    fz(i,j,k)=0.1_mytype*sin(0.3_mytype*j)
  end do; end do; end do
  rhsx0=rhsx; rhsy0=rhsy; rhsz0=rhsz

  open(newunit=io,file='stage6_outputs/fibre_stage6_micro_smoke_check.dat',status='replace',action='write')

  call init_stage6_controlled_test_config(config); config%rho_fluid=2._mytype
  call init_stage6_uniform_collocated_layout(layout,nx,ny,nz)
  rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call perform_stage6_micro_controlled_step(config,layout,dt,fx,fy,fz,rhsx,rhsy,rhsz,ux,uy,uz,uxs,uys,uzs,cv,lg,hk,ric,inj,mod,rej,re,ex,ey,ez,ue)
  call compute_stage6_micro_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,rhschg)
  bufmax=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))
  call compute_stage6_micro_velocity_change_max(ux,uy,uz,uxs,uys,uzs,ustarchg)
  write(io,'(A,1X,I0)') 'stage6_micro_controlled_config_valid_flag', cv
  write(io,'(A,1X,I0)') 'stage6_micro_layout_guard_pass_flag', lg
  write(io,'(A,1X,I0)') 'stage6_micro_hook_called_flag', hk
  write(io,'(A,1X,I0)') 'stage6_micro_rhs_injection_called_flag', ric
  write(io,'(A,1X,I0)') 'stage6_micro_controlled_path_status', merge(1,0,cv==1 .and. lg==1 .and. hk==1 .and. ric==1)
  write(io,'(A,1X,ES24.16)') 'stage6_micro_buffer_max_abs', bufmax
  write(io,'(A,1X,ES24.16)') 'stage6_micro_rhs_change_max', rhschg
  write(io,'(A,1X,ES24.16)') 'stage6_micro_rhs_expected_error', re
  write(io,'(A,1X,I0)') 'stage6_micro_injected_flag', inj
  write(io,'(A,1X,I0)') 'stage6_micro_modified_flag', mod
  write(io,'(A,1X,I0)') 'stage6_micro_rejected_flag', rej
  write(io,'(A,1X,ES24.16)') 'stage6_micro_component_x_error', ex
  write(io,'(A,1X,ES24.16)') 'stage6_micro_component_y_error', ey
  write(io,'(A,1X,ES24.16)') 'stage6_micro_component_z_error', ez
  write(io,'(A,1X,ES24.16)') 'stage6_micro_ustar_expected_error', ue
  write(io,'(A,1X,ES24.16)') 'stage6_micro_ustar_change_norm', ustarchg
  write(io,'(A,1X,I0)') 'stage6_micro_dt_positive_flag', merge(1,0,dt>0._mytype)

  call init_stage6_default_config(config)
  call init_stage6_uniform_collocated_layout(layout,nx,ny,nz)
  rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call perform_stage6_micro_controlled_step(config,layout,dt,fx,fy,fz,rhsx,rhsy,rhsz,ux,uy,uz,uxs,uys,uzs,cv,lg,hk,ric,inj,mod,rej,re,ex,ey,ez,ue)
  call compute_stage6_intermediate_velocity(dt,ux,uy,uz,rhsx0,rhsy0,rhsz0,uxs0,uys0,uzs0)
  call compute_stage6_micro_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,rhschg)
  tmp=max(maxval(abs(uxs-uxs0)),max(maxval(abs(uys-uys0)),maxval(abs(uzs-uzs0))))
  write(io,'(A,1X,ES24.16)') 'stage6_micro_default_rhs_change_max', rhschg
  write(io,'(A,1X,ES24.16)') 'stage6_micro_default_ustar_diff_max', tmp
  write(io,'(A,1X,I0)') 'stage6_micro_default_injected_flag', inj
  write(io,'(A,1X,I0)') 'stage6_micro_default_modified_flag', mod

  call init_stage6_controlled_test_config(config); config%rho_fluid=2._mytype
  call init_stage6_nonuniform_y_layout(layout,nx,ny,nz)
  rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call perform_stage6_micro_controlled_step(config,layout,dt,fx,fy,fz,rhsx,rhsy,rhsz,ux,uy,uz,uxs,uys,uzs,cv,lg,hk,ric,inj,mod,rej,re,ex,ey,ez,ue)
  call compute_stage6_micro_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,rhschg)
  write(io,'(A,1X,I0)') 'stage6_micro_blocked_layout_blocked_flag', merge(1,0,layout%blocked_flag==1)
  write(io,'(A,1X,ES24.16)') 'stage6_micro_blocked_layout_rhs_change_max', rhschg
  write(io,'(A,1X,I0)') 'stage6_micro_blocked_layout_injected_flag', inj
  write(io,'(A,1X,I0)') 'stage6_micro_blocked_layout_modified_flag', mod
  write(io,'(A,1X,I0)') 'stage6_micro_blocked_layout_rhs_injection_called_flag', ric

  config%enable_stage6=.true.; config%enable_main_rhs_hook=.true.; config%enable_controlled_rhs_test=.false.
  config%production_two_way_enabled=.false.; config%allow_stage5_hook_in_main_path=.true.; config%reject_invalid_config=.true.; config%rho_fluid=1._mytype
  call init_stage6_uniform_collocated_layout(layout,nx,ny,nz)
  rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call perform_stage6_micro_controlled_step(config,layout,dt,fx,fy,fz,rhsx,rhsy,rhsz,ux,uy,uz,uxs,uys,uzs,cv,lg,hk,ric,inj,mod,rej,re,ex,ey,ez,ue)
  call compute_stage6_micro_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,rhschg)
  write(io,'(A,1X,I0)') 'stage6_micro_invalid_rejected_flag', rej
  write(io,'(A,1X,ES24.16)') 'stage6_micro_invalid_rhs_change_max', rhschg
  write(io,'(A,1X,I0)') 'stage6_micro_invalid_injected_flag', inj
  write(io,'(A,1X,I0)') 'stage6_micro_invalid_modified_flag', mod

  call init_stage6_controlled_test_config(config); config%production_two_way_enabled=.true.; config%rho_fluid=2._mytype
  call init_stage6_uniform_collocated_layout(layout,nx,ny,nz)
  rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call perform_stage6_micro_controlled_step(config,layout,dt,fx,fy,fz,rhsx,rhsy,rhsz,ux,uy,uz,uxs,uys,uzs,cv,lg,hk,ric,inj,mod,rej,re,ex,ey,ez,ue)
  call compute_stage6_micro_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,rhschg)
  write(io,'(A,1X,I0)') 'stage6_micro_production_enabled_rejected_flag', rej
  write(io,'(A,1X,ES24.16)') 'stage6_micro_production_enabled_rhs_change_max', rhschg
  write(io,'(A,1X,I0)') 'stage6_micro_production_enabled_injected_flag', inj
  write(io,'(A,1X,I0)') 'stage6_micro_production_enabled_modified_flag', mod

  call init_stage6_controlled_test_config(config); config%rho_fluid=2._mytype
  call init_stage6_uniform_collocated_layout(layout,nx,ny,nz)
  rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0; fx=0._mytype;fy=0._mytype;fz=0._mytype
  call perform_stage6_micro_controlled_step(config,layout,dt,fx,fy,fz,rhsx,rhsy,rhsz,ux,uy,uz,uxs,uys,uzs,cv,lg,hk,ric,inj,mod,rej,re,ex,ey,ez,ue)
  call compute_stage6_micro_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,rhschg)
  write(io,'(A,1X,ES24.16)') 'stage6_micro_zero_buffer_rhs_change_max', rhschg
  write(io,'(A,1X,I0)') 'stage6_micro_zero_buffer_injected_flag', inj
  write(io,'(A,1X,I0)') 'stage6_micro_zero_buffer_modified_flag', mod

  call stage6_micro_pressure_production_status(ppm,prm,pm,rpc,pvm,pdc,pebd)
  write(io,'(A,1X,I0)') 'stage6_micro_pressure_poisson_modified_flag', ppm
  write(io,'(A,1X,I0)') 'stage6_micro_pressure_rhs_modified_flag', prm
  write(io,'(A,1X,I0)') 'stage6_micro_projection_modified_flag', pm
  write(io,'(A,1X,I0)') 'stage6_micro_real_projection_called_flag', rpc
  write(io,'(A,1X,I0)') 'stage6_micro_post_projection_velocity_modified_flag', pvm
  write(io,'(A,1X,I0)') 'stage6_micro_production_dns_called_flag', pdc
  write(io,'(A,1X,I0)') 'stage6_micro_production_enabled_by_default_flag', pebd
  write(io,'(A,1X,I0)') 'stage6_micro_smoke_check_status', 1

  close(io)
end program fibre_stage6_micro_smoke_check
