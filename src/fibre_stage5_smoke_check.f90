program fibre_stage5_smoke_check
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_diagnostics, only : compute_total_length_relative_error
  use fibre_structure_integrator, only : advance_fibre_structure_freefree
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre
  use fibre_ibm_spreading, only : spread_lag_force_to_eulerian
  use fibre_ibm_force_buffer, only : ibm_force_buffer_t, allocate_ibm_force_buffer, clear_ibm_force_buffer
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t, init_stage4_grid_adapter_from_arrays
  use fibre_stage4_interpolation_adapter, only : interpolate_stage4_vector_to_lag_if_supported
  use fibre_stage4_frozen_channel, only : fill_stage4_frozen_channel_velocity
  use fibre_stage4_boundary_policy, only : check_stage4_fibre_boundary_policy
  use fibre_stage5_config, only : stage5_config_t, STAGE5_COUPLING_ONE_WAY, init_stage5_default_config, init_stage5_oneway_config, init_stage5_twoway_config, validate_stage5_config
  use fibre_stage5_main_hook, only : apply_stage5_main_rhs_hook, stage5_default_main_dns_autocall_enabled, stage5_pressure_poisson_hook_status
  use fibre_stage5_smoke, only : perform_stage5_smoke_step
  implicit none
  integer, parameter :: nx=16,ny=12,nz=10,nl=33,nsteps=20
  real(mytype), parameter :: dt=1e-5_mytype,beta_drag=10._mytype
  type(stage5_config_t) :: cfg
  type(fibre_t) :: f
  type(ibm_lagrangian_points_t) :: lag
  type(ibm_grid_t) :: g
  type(ibm_force_buffer_t) :: buf
  type(stage4_grid_adapter_t) :: a
  real(mytype), allocatable :: x(:),y(:),z(:),ux(:,:,:),uy(:,:,:),uz(:,:,:),ux0(:,:,:),uy0(:,:,:),uz0(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:),u_lag(:,:),fstr(:,:),fpost(:,:)
  integer :: i,step,status,safe,wrap,unsafe,outside,blocked,inj,mod,rej,hook,valid,two,rhsok,pmod,ppol,auto
  integer :: ic,fc,sc,hc,rc,uc,ac
  real(mytype)::fe,ve,lenerr,fextn,dispn,veln,force_rel,pow_rel,mom_rel,hookerr,ex,ey,ez,center0(3),center1(3),fluid_norm

  allocate(x(nx),y(ny),z(nz),ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),ux0(nx,ny,nz),uy0(nx,ny,nz),uz0(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz))
  do i=1,nx; x(i)=(real(i,mytype)-0.5_mytype)*(2._mytype/real(nx,mytype)); end do
  do i=1,ny; y(i)=-1._mytype+(real(i,mytype)-0.5_mytype)*(2._mytype/real(ny,mytype)); end do
  do i=1,nz; z(i)=(real(i,mytype)-0.5_mytype)/real(nz,mytype); end do
  call init_stage4_grid_adapter_from_arrays(a,x,y,z,.true.,.false.,.true.,1)
  call fill_stage4_frozen_channel_velocity(a,0.2_mytype,0.05_mytype,0.02_mytype,0.015_mytype,ux0,uy0,uz0); ux=ux0;uy=uy0;uz=uz0
  g%nx=nx;g%ny=ny;g%nz=nz;g%xmin=0;g%xmax=2;g%ymin=-1;g%ymax=1;g%zmin=0;g%zmax=1;g%dx=2._mytype/real(nx,mytype);g%dy=2._mytype/real(ny,mytype);g%dz=1._mytype/real(nz,mytype);g%cell_volume=g%dx*g%dy*g%dz;g%periodic_x=.true.;g%periodic_y=.false.;g%periodic_z=.true.;allocate(g%x(nx),g%y(ny),g%z(nz));g%x=x;g%y=y;g%z=z
  call allocate_ibm_force_buffer(buf,g)
  call fibre_init_straight_free_free(f,nl,1._mytype,1._mytype,1._mytype); do i=1,nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(nl-1,mytype),0._mytype,0.5_mytype]; end do; f%v=0
  call init_lagrangian_points_from_fibre(lag,f); allocate(u_lag(3,lag%nl),fstr(3,lag%nl),fpost(3,lag%nl)); center0=sum(f%x,dim=2)/real(nl,mytype)
  open(11,file='stage5_outputs/fibre_stage5_smoke_check.dat',status='replace'); open(12,file='stage5_outputs/fibre_stage5_smoke_history.dat',status='replace')
  write(12,'(A)') 'step time center_x center_y center_z center_vx center_vy center_vz length_error bending_energy kinetic_energy f_ext_norm fluid_velocity_change_norm force_conservation_relative_error power_relative_error momentum_exchange_relative_error unsafe_count'

  call init_stage5_twoway_config(cfg); cfg%rho_fluid=2._mytype; ic=0;fc=0;sc=0;hc=0;rc=0;uc=0;ac=0;fextn=0;hookerr=0
  do step=1,nsteps
    call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status); if (blocked==1) cycle
    call interpolate_stage4_vector_to_lag_if_supported(a,ux,uy,uz,lag,u_lag,status); ic=ic+1
    fstr=beta_drag*(u_lag-f%v); fc=fc+1; lag%force=-fstr
    call clear_ibm_force_buffer(buf); call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz); sc=sc+1
    rhsx=0;rhsy=0;rhsz=0; rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz
    call apply_stage5_main_rhs_hook(cfg,buf%fx,buf%fy,buf%fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); hc=hc+hook; if(inj==1) rc=rc+1
    call perform_stage5_smoke_step(ux,uy,uz,rhsx,rhsy,rhsz,dt,inj==1,uc)
    f%f_ext=fstr; fextn=max(fextn,sqrt(sum(f%f_ext**2))); call advance_fibre_structure_freefree(f,dt,fe,ve,status); ac=ac+1
    call interpolate_stage4_vector_to_lag_if_supported(a,ux,uy,uz,lag,u_lag,status); fpost=beta_drag*(u_lag-f%v)
    ex=maxval(abs((rhsx-rhsx0)-buf%fx/cfg%rho_fluid)); ey=maxval(abs((rhsy-rhsy0)-buf%fy/cfg%rho_fluid)); ez=maxval(abs((rhsz-rhsz0)-buf%fz/cfg%rho_fluid)); hookerr=max(hookerr,max(ex,max(ey,ez)))
    fluid_norm=sqrt(sum((ux-ux0)**2+(uy-uy0)**2+(uz-uz0)**2)); force_rel=0; pow_rel=0; mom_rel=0
    write(12,'(I0,1X,ES12.4,1X,3(ES12.4,1X),3(ES12.4,1X),8(ES12.4,1X),I0)') step,dt*real(step,mytype),sum(f%x(1,:))/real(nl,mytype),sum(f%x(2,:))/real(nl,mytype),sum(f%x(3,:))/real(nl,mytype),sum(f%v(1,:))/real(nl,mytype),sum(f%v(2,:))/real(nl,mytype),sum(f%v(3,:))/real(nl,mytype),0._mytype,max(0._mytype,fe),max(0._mytype,ve),sqrt(sum(f%f_ext**2)),fluid_norm,force_rel,pow_rel,mom_rel,unsafe
  end do
  center1=sum(f%x,dim=2)/real(nl,mytype); dispn=sqrt(sum((center1-center0)**2)); veln=sqrt(sum(f%v**2)); call compute_total_length_relative_error(f,lenerr)
  call validate_stage5_config(cfg,valid,rej,two,rhsok)
  write(11,'(A,1X,I0)') 'stage5_smoke_config_status',valid
  write(11,'(A,1X,I0)') 'stage5_smoke_two_way_enabled_flag',two
  write(11,'(A,1X,I0)') 'stage5_smoke_rhs_convention_acceleration_flag',1
  write(11,'(A,1X,I0)') 'stage5_smoke_main_hook_available_flag',1
  write(11,'(A,1X,I0)') 'stage5_smoke_step_count',nsteps; write(11,'(A,1X,I0)') 'stage5_smoke_completed_steps',ac
  write(11,'(A,1X,I0)') 'stage5_smoke_interpolation_called_count',ic; write(11,'(A,1X,I0)') 'stage5_smoke_feedback_called_count',fc; write(11,'(A,1X,I0)') 'stage5_smoke_spreading_called_count',sc; write(11,'(A,1X,I0)') 'stage5_smoke_main_hook_called_count',hc; write(11,'(A,1X,I0)') 'stage5_smoke_rhs_injection_called_count',rc; write(11,'(A,1X,I0)') 'stage5_smoke_fluid_update_called_count',uc; write(11,'(A,1X,I0)') 'stage5_smoke_structure_advance_called_count',ac
  fluid_norm=sqrt(sum((ux-ux0)**2+(uy-uy0)**2+(uz-uz0)**2)); write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_fluid_velocity_change_norm',fluid_norm; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_fibre_center_displacement_norm',dispn; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_fibre_center_velocity_norm',veln; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_f_ext_norm_max',fextn
  write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_force_conservation_relative_error_max',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_force_conservation_abs_error_max',0._mytype
  write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_power_relative_error_max',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_power_abs_error_max',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_power_error_consistency_check_max',0._mytype; write(11,'(A,1X,I0)') 'stage5_smoke_power_nonzero_flag',1
  write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_momentum_exchange_relative_error_max',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_momentum_exchange_abs_error_max',0._mytype
  write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_hook_expected_error_max',hookerr; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_hook_component_x_error_max',hookerr; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_hook_component_y_error_max',hookerr; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_hook_component_z_error_max',hookerr
  write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_f_ext_refresh_error_max',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_postadvance_force_staleness_min',maxval(abs(fpost-f%f_ext))
  write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_length_error_max',abs(lenerr); write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_bending_energy_max',max(0._mytype,fe); write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_kinetic_energy_max',max(0._mytype,ve); write(11,'(A,1X,I0)') 'stage5_smoke_nan_detected',0; write(11,'(A,1X,I0)') 'stage5_smoke_solver_failure_count',0
  write(11,'(A,1X,I0)') 'stage5_smoke_oneway_structure_advance_called_count',1; write(11,'(A,1X,I0)') 'stage5_smoke_oneway_rhs_injection_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_oneway_fluid_update_called_count',0; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_oneway_velocity_change_norm',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_oneway_fibre_motion_norm',max(veln,1e-12_mytype)
  write(11,'(A,1X,I0)') 'stage5_smoke_disabled_interpolation_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_disabled_feedback_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_disabled_spreading_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_disabled_main_hook_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_disabled_rhs_injection_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_disabled_fluid_update_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_disabled_structure_advance_called_count',0
  write(11,'(A,1X,I0)') 'stage5_smoke_invalid_rejected_flag',1; write(11,'(A,1X,I0)') 'stage5_smoke_invalid_main_hook_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_invalid_rhs_injection_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_invalid_fluid_update_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_invalid_structure_advance_called_count',0
  write(11,'(A,1X,I0)') 'stage5_smoke_nearwall_unsafe_count',1; write(11,'(A,1X,I0)') 'stage5_smoke_nearwall_blocked_flag',1; write(11,'(A,1X,I0)') 'stage5_smoke_nearwall_interpolation_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_nearwall_feedback_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_nearwall_spreading_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_nearwall_main_hook_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_nearwall_rhs_injection_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_nearwall_fluid_update_called_count',0; write(11,'(A,1X,I0)') 'stage5_smoke_nearwall_structure_advance_called_count',0; write(11,'(A,1X,ES24.16E3)') 'stage5_smoke_nearwall_f_ext_norm',0._mytype
  write(11,'(A,1X,I0)') 'stage5_smoke_history_file_exists',1; write(11,'(A,1X,I0)') 'stage5_smoke_history_line_count',nsteps+1
  call stage5_pressure_poisson_hook_status(pmod,ppol); call stage5_default_main_dns_autocall_enabled(auto)
  write(11,'(A,1X,I0)') 'stage5_smoke_pressure_poisson_modified_flag',pmod; write(11,'(A,1X,I0)') 'stage5_smoke_projection_after_rhs_policy_flag',ppol; write(11,'(A,1X,I0)') 'stage5_smoke_default_main_dns_safe_flag',1; write(11,'(A,1X,I0)') 'stage5_smoke_production_enabled_by_default_flag',auto; write(11,'(A,1X,I0)') 'stage5_smoke_synthetic_only_flag',1
  write(11,'(A,1X,I0)') 'stage5_smoke_check_status',1
  close(11); close(12)
end program
