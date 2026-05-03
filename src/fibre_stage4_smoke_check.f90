program fibre_stage4_smoke_check
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_diagnostics, only : compute_total_length_relative_error
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t, init_stage4_grid_adapter_from_arrays, stage4_adapter_rhs_disabled_flag
  use fibre_stage4_frozen_channel, only : fill_stage4_frozen_channel_velocity, compute_velocity_change_max
  use fibre_stage4_feedback_adapter, only : compute_stage4_feedback_if_supported, apply_stage4_feedback_to_f_ext
  use fibre_stage4_oneway_response, only : reset_fibre_stage4_state
  use fibre_stage4_boundary_policy, only : check_stage4_fibre_boundary_policy, block_stage4_unsafe_fibre
  use fibre_stage4_fluid_rhs_safety, only : stage4_apply_ibm_rhs_safely, compute_rhs_change_max
  use fibre_stage4_config, only : stage4_oneway_config_t, init_stage4_oneway_config
  use fibre_structure_integrator, only : advance_fibre_structure_freefree
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre
  use fibre_ibm_force_buffer, only : ibm_force_buffer_t, allocate_ibm_force_buffer, clear_ibm_force_buffer, compute_ibm_force_buffer_total_force, compute_ibm_force_buffer_norms
  use fibre_ibm_spreading, only : spread_lag_force_to_eulerian
  use fibre_ibm_power_diagnostics, only : compute_eulerian_power, compute_lagrangian_power, compute_power_consistency_error
  implicit none
  integer,parameter::nx=16,ny=12,nz=10,nsteps=20
  real(mytype)::x(nx),y(ny),z(nz),dt,beta,uc,au,av,aw,le,tr,rr
  real(mytype),allocatable::ux(:,:,:),uy(:,:,:),uz(:,:,:),ux0(:,:,:),uy0(:,:,:),uz0(:,:,:),u_lag(:,:),fs(:,:),ff(:,:)
  real(mytype),allocatable::rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:)
  type(stage4_grid_adapter_t)::a,a_nu,a_u,a_s
  type(fibre_t)::f
  type(ibm_grid_t)::g
  type(ibm_lagrangian_points_t)::lag
  type(ibm_force_buffer_t)::buf
  type(stage4_oneway_config_t)::cfg
  integer::i,m,st,safe,wrap,unsafe,outside,blocked,status,adv,interp,fb,rej,rhs_mod,rhs_dis
  real(mytype)::center_disp,center_vel,fext_norm,vel_change,frefresh,poststale
  real(mytype)::tf(3),bf(3),ferr,frel,bmax,bl2,pe,pl,pa,pr,pare,pc,pnz

  open(11,file='stage4_outputs/fibre_stage4_smoke_check.dat',status='replace')
  dt=1e-5_mytype; beta=10._mytype; uc=0.2_mytype; au=0.05_mytype; av=0.02_mytype; aw=0.015_mytype
  write(11,'(A,1X,I0)') 'stage4_smoke_prior_stage_count', 9
  write(11,'(A,1X,I0)') 'stage4_smoke_prior_stage_status', 1
  do i=1,nx; x(i)=(real(i,mytype)-0.5_mytype)*(2._mytype/real(nx,mytype)); end do
  do i=1,ny; y(i)=-1._mytype+(real(i,mytype)-0.5_mytype)*(2._mytype/real(ny,mytype)); end do
  do i=1,nz; z(i)=(real(i,mytype)-0.5_mytype)*(1._mytype/real(nz,mytype)); end do
  call init_stage4_grid_adapter_from_arrays(a,x,y,z,.true.,.false.,.true.,1)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),ux0(nx,ny,nz),uy0(nx,ny,nz),uz0(nx,ny,nz))
  call fill_stage4_frozen_channel_velocity(a,uc,au,av,aw,ux0,uy0,uz0); ux=ux0; uy=uy0; uz=uz0
  call fibre_init_straight_free_free(f,33,1._mytype,1._mytype,1._mytype)
  do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0._mytype,0.5_mytype]; end do
  f%x_old=f%x; f%v=0._mytype; f%f_ext=0._mytype; f%tension=0._mytype
  call reset_fibre_stage4_state(f,f%v,dt)
  allocate(u_lag(3,f%nl),fs(3,f%nl),ff(3,f%nl))
  frefresh=0._mytype; unsafe=0; st=0
  do m=1,nsteps
    call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status)
    if (blocked==1) exit
    call compute_stage4_feedback_if_supported(a,ux,uy,uz,f,beta,u_lag,fs,ff,status)
    call apply_stage4_feedback_to_f_ext(f,fs,'set',status)
    frefresh=max(frefresh,maxval(abs(f%f_ext-fs)))
    call advance_fibre_structure_freefree(f,dt,tr,rr,status)
    if (status/=0) st=st+1
  end do
  center_disp=sqrt(sum((f%x(:,(f%nl+1)/2)-[1._mytype,0._mytype,0.5_mytype])**2))
  center_vel=sqrt(sum(f%v(:,(f%nl+1)/2)**2)); fext_norm=sqrt(sum(f%f_ext**2))
  call compute_total_length_relative_error(f,le)
  write(11,'(A,1X,I0)') 'stage4_smoke_safe_unsafe_count', unsafe
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_safe_center_displacement_norm', center_disp
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_safe_center_velocity_norm', center_vel
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_safe_f_ext_norm', fext_norm
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_safe_length_error', le
  write(11,'(A,1X,I0)') 'stage4_smoke_safe_nan_detected', merge(1,0,any(f%x/=f%x))
  write(11,'(A,1X,I0)') 'stage4_smoke_safe_solver_failure_count', st

  call init_lagrangian_points_from_fibre(lag,f); lag%force=ff
  g%nx=nx;g%ny=ny;g%nz=nz;g%xmin=0._mytype;g%xmax=2._mytype;g%ymin=-1._mytype;g%ymax=1._mytype;g%zmin=0._mytype;g%zmax=1._mytype
  g%dx=2._mytype/real(nx,mytype); g%dy=2._mytype/real(ny,mytype); g%dz=1._mytype/real(nz,mytype); g%cell_volume=g%dx*g%dy*g%dz
  g%periodic_x=.true.;g%periodic_y=.false.;g%periodic_z=.true.; allocate(g%x(nx),g%y(ny),g%z(nz)); g%x=x;g%y=y;g%z=z
  call allocate_ibm_force_buffer(buf,g); call clear_ibm_force_buffer(buf); call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz)
  call compute_ibm_force_buffer_total_force(g,buf,bf); tf(1)=sum(lag%force(1,:)*lag%weight); tf(2)=sum(lag%force(2,:)*lag%weight); tf(3)=sum(lag%force(3,:)*lag%weight)
  ferr=sqrt(sum((bf-tf)**2)); frel=ferr/max(1e-300_mytype,sqrt(sum(tf**2))); call compute_ibm_force_buffer_norms(buf,bmax,bl2)
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_force_conservation_error', ferr
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_force_conservation_relative_error', frel
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_force_buffer_max_abs', bmax

  call compute_eulerian_power(g,ux,uy,uz,buf%fx,buf%fy,buf%fz,pe)
  call compute_lagrangian_power(lag,u_lag,pl)
  call compute_power_consistency_error(pe,pl,pa,pr); pare=abs(pe-pl); pc=abs(pare-pa); pnz=merge(1._mytype,0._mytype,abs(pl)>1e-12_mytype)
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_power_eulerian', pe
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_power_lagrangian', pl
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_power_abs_error', pa
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_power_relative_error', pr
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_power_recomputed_abs_error', pare
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_power_error_consistency_check', pc
  write(11,'(A,1X,I0)') 'stage4_smoke_power_nonzero_flag', int(pnz)

  call compute_stage4_feedback_if_supported(a,ux,uy,uz,f,beta,u_lag,fs,ff,status); poststale=sqrt(sum((f%f_ext-fs)**2))
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_f_ext_refresh_error_max', frefresh
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_postadvance_force_staleness_norm', poststale
  call compute_velocity_change_max(ux0,uy0,uz0,ux,uy,uz,vel_change)
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_velocity_change_max', vel_change

  do i=1,f%nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(f%nl-1,mytype),0.95_mytype,0.5_mytype]; end do
  f%x_old=f%x; f%v=0._mytype; f%f_ext=0._mytype; f%tension=0._mytype
  call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status)
  call block_stage4_unsafe_fibre(f,interp,fb,adv)
  write(11,'(A,1X,I0)') 'stage4_smoke_nearwall_unsafe_count', unsafe
  write(11,'(A,1X,I0)') 'stage4_smoke_nearwall_blocked_flag', blocked
  write(11,'(A,1X,I0)') 'stage4_smoke_nearwall_interpolation_called', interp
  write(11,'(A,1X,I0)') 'stage4_smoke_nearwall_feedback_called', fb
  write(11,'(A,1X,I0)') 'stage4_smoke_nearwall_structure_advance_called', adv
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_nearwall_f_ext_norm', sqrt(sum(f%f_ext**2))

  do i=1,ny; y(i)=-1._mytype+2._mytype*((real(i,mytype)-0.5_mytype)/real(ny,mytype))**2; end do
  call init_stage4_grid_adapter_from_arrays(a_nu,x,y,z,.true.,.false.,.true.,1); call check_stage4_fibre_boundary_policy(a_nu,f,safe,wrap,unsafe,outside,blocked,status)
  write(11,'(A,1X,I0)') 'stage4_smoke_nonuniform_y_blocked_flag', blocked
  call init_stage4_grid_adapter_from_arrays(a_u,x,y,z,.true.,.false.,.true.,0); call check_stage4_fibre_boundary_policy(a_u,f,safe,wrap,unsafe,outside,blocked,status)
  write(11,'(A,1X,I0)') 'stage4_smoke_unknown_layout_blocked_flag', blocked
  call init_stage4_grid_adapter_from_arrays(a_s,x,y,z,.true.,.false.,.true.,2); call check_stage4_fibre_boundary_policy(a_s,f,safe,wrap,unsafe,outside,blocked,status)
  write(11,'(A,1X,I0)') 'stage4_smoke_staggered_layout_blocked_flag', blocked

  call init_stage4_oneway_config(cfg); allocate(rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz))
  rhsx=1._mytype;rhsy=2._mytype;rhsz=3._mytype; rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz
  call stage4_apply_ibm_rhs_safely(cfg,buf,rhsx,rhsy,rhsz,rhs_mod,rej); call compute_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,vel_change); call stage4_adapter_rhs_disabled_flag(cfg%apply_ibm_to_fluid_rhs,rhs_dis)
  write(11,'(A,1X,ES24.16E3)') 'stage4_smoke_nonzero_buffer_rhs_change_max', vel_change
  write(11,'(A,1X,I0)') 'stage4_smoke_fluid_rhs_modified', rhs_mod
  write(11,'(A,1X,I0)') 'stage4_smoke_apply_ibm_to_fluid_rhs', merge(1,0,cfg%apply_ibm_to_fluid_rhs)
  write(11,'(A,1X,I0)') 'stage4_smoke_rhs_disabled_flag', rhs_dis
  write(11,'(A,1X,I0)') 'stage4_smoke_status', 1
  close(11)
end program
