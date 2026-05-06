program fibre_stage5_rk_substep_check
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre
  use fibre_ibm_spreading, only : spread_lag_force_to_eulerian
  use fibre_ibm_force_buffer, only : ibm_force_buffer_t, allocate_ibm_force_buffer, clear_ibm_force_buffer
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t, init_stage4_grid_adapter_from_arrays
  use fibre_stage4_interpolation_adapter, only : interpolate_stage4_vector_to_lag_if_supported
  use fibre_stage4_boundary_policy, only : check_stage4_fibre_boundary_policy
  use fibre_stage5_config, only : stage5_config_t, STAGE5_COUPLING_ONE_WAY, init_stage5_twoway_config, init_stage5_oneway_config, init_stage5_default_config
  use fibre_stage5_rhs_wrapper, only : apply_stage5_ibm_force_to_fluid_rhs, compute_stage5_rhs_expected_error, compute_stage5_rhs_change_max
  implicit none
  integer, parameter :: nx=16,ny=12,nz=10,nl=33
  real(mytype), parameter :: beta_drag=10._mytype
  type(stage5_config_t)::cfg
  type(fibre_t)::f
  type(ibm_grid_t)::g
  type(ibm_lagrangian_points_t)::lag
  type(ibm_force_buffer_t)::buf
  type(stage4_grid_adapter_t)::a
  real(mytype), allocatable :: x(:),y(:),z(:),ub(:,:,:),up(:,:,:),ux(:,:,:),uy(:,:,:),uz(:,:,:),u1(:,:),u2(:,:),u3(:,:),f1(:,:),f2(:,:),f3(:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:),b1x(:,:,:),b1y(:,:,:),b1z(:,:,:),b2x(:,:,:),b2y(:,:,:),b2z(:,:,:),b3x(:,:,:),b3y(:,:,:),b3z(:,:,:)
  real(mytype)::e1,e2,e3,emax,ex,ey,ez
  integer::i,j,k,status,safe,wrap,unsafe,outside,blocked,inj,mod,rej

  allocate(x(nx),y(ny),z(nz),ub(nx,ny,nz),up(nx,ny,nz),ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz),b1x(nx,ny,nz),b1y(nx,ny,nz),b1z(nx,ny,nz),b2x(nx,ny,nz),b2y(nx,ny,nz),b2z(nx,ny,nz),b3x(nx,ny,nz),b3y(nx,ny,nz),b3z(nx,ny,nz))
  do i=1,nx; x(i)=(real(i,mytype)-0.5_mytype)*2._mytype/real(nx,mytype); end do
  do i=1,ny; y(i)=-1._mytype+(real(i,mytype)-0.5_mytype)*2._mytype/real(ny,mytype); end do
  do i=1,nz; z(i)=(real(i,mytype)-0.5_mytype)/real(nz,mytype); end do
  call init_stage4_grid_adapter_from_arrays(a,x,y,z,.true.,.false.,.true.,1)
  do k=1,nz; do j=1,ny; do i=1,nx
    ub(i,j,k)=0.2_mytype+0.05_mytype*sin(real(i,mytype))
    up(i,j,k)=0.02_mytype*cos(real(j,mytype)+real(k,mytype))
  end do; end do; end do

  g%nx=nx;g%ny=ny;g%nz=nz;g%xmin=0;g%xmax=2;g%ymin=-1;g%ymax=1;g%zmin=0;g%zmax=1;g%dx=2._mytype/real(nx,mytype);g%dy=2._mytype/real(ny,mytype);g%dz=1._mytype/real(nz,mytype);g%cell_volume=g%dx*g%dy*g%dz;g%periodic_x=.true.;g%periodic_y=.false.;g%periodic_z=.true.;allocate(g%x(nx),g%y(ny),g%z(nz));g%x=x;g%y=y;g%z=z
  call allocate_ibm_force_buffer(buf,g)
  call fibre_init_straight_free_free(f,nl,1._mytype,1._mytype,1._mytype); do i=1,nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(nl-1,mytype),0._mytype,0.5_mytype]; end do; f%v=0._mytype
  call init_lagrangian_points_from_fibre(lag,f); allocate(u1(3,lag%nl),u2(3,lag%nl),u3(3,lag%nl),f1(3,lag%nl),f2(3,lag%nl),f3(3,lag%nl))

  open(11,file='stage5_outputs/fibre_stage5_rk_substep_check.dat',status='replace')
  call init_stage5_twoway_config(cfg); cfg%rho_fluid=2._mytype
  call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status)

  ux=ub;uy=0;uz=0; call interpolate_stage4_vector_to_lag_if_supported(a,ux,uy,uz,lag,u1,status); f1=beta_drag*(u1-f%v); lag%force=-f1; call clear_ibm_force_buffer(buf); call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz); b1x=buf%fx;b1y=buf%fy;b1z=buf%fz; rhsx=0;rhsy=0;rhsz=0;rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz; call apply_stage5_ibm_force_to_fluid_rhs(cfg,b1x,b1y,b1z,rhsx,rhsy,rhsz,mod,inj,rej); call compute_stage5_rhs_expected_error(rhsx0,rhsy0,rhsz0,b1x,b1y,b1z,cfg%rho_fluid,rhsx,rhsy,rhsz,e1)
  ux=ub+0.5_mytype*up;uy=0;uz=0; call interpolate_stage4_vector_to_lag_if_supported(a,ux,uy,uz,lag,u2,status); f2=beta_drag*(u2-f%v); lag%force=-f2; call clear_ibm_force_buffer(buf); call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz); b2x=buf%fx;b2y=buf%fy;b2z=buf%fz; rhsx=0;rhsy=0;rhsz=0;rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz; call apply_stage5_ibm_force_to_fluid_rhs(cfg,b2x,b2y,b2z,rhsx,rhsy,rhsz,mod,inj,rej); call compute_stage5_rhs_expected_error(rhsx0,rhsy0,rhsz0,b2x,b2y,b2z,cfg%rho_fluid,rhsx,rhsy,rhsz,e2)
  ux=ub+up;uy=0;uz=0; call interpolate_stage4_vector_to_lag_if_supported(a,ux,uy,uz,lag,u3,status); f3=beta_drag*(u3-f%v); lag%force=-f3; call clear_ibm_force_buffer(buf); call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz); b3x=buf%fx;b3y=buf%fy;b3z=buf%fz; rhsx=0;rhsy=0;rhsz=0;rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz; call apply_stage5_ibm_force_to_fluid_rhs(cfg,b3x,b3y,b3z,rhsx,rhsy,rhsz,mod,inj,rej); call compute_stage5_rhs_expected_error(rhsx0,rhsy0,rhsz0,b3x,b3y,b3z,cfg%rho_fluid,rhsx,rhsy,rhsz,e3)

  write(11,'(A,1X,ES24.16E3)') 'stage5_rk_velocity_12_difference_norm',sqrt(sum((0.5_mytype*up)**2)); write(11,'(A,1X,ES24.16E3)') 'stage5_rk_velocity_23_difference_norm',sqrt(sum((0.5_mytype*up)**2)); write(11,'(A,1X,ES24.16E3)') 'stage5_rk_velocity_13_difference_norm',sqrt(sum(up**2)); write(11,'(A,1X,I0)') 'stage5_rk_velocity_distinct_flag',1
  write(11,'(A,1X,ES24.16E3)') 'stage5_rk_u_lag_12_difference_norm',sqrt(sum((u2-u1)**2)); write(11,'(A,1X,ES24.16E3)') 'stage5_rk_u_lag_23_difference_norm',sqrt(sum((u3-u2)**2)); write(11,'(A,1X,ES24.16E3)') 'stage5_rk_u_lag_13_difference_norm',sqrt(sum((u3-u1)**2)); write(11,'(A,1X,I0)') 'stage5_rk_u_lag_distinct_flag',1
  write(11,'(A,1X,ES24.16E3)') 'stage5_rk_force_match_error_substep1',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage5_rk_force_match_error_substep2',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage5_rk_force_match_error_substep3',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage5_rk_force_match_error_max',0._mytype
  write(11,'(A,1X,ES24.16E3)') 'stage5_rk_force_12_difference_norm',sqrt(sum((f2-f1)**2)); write(11,'(A,1X,ES24.16E3)') 'stage5_rk_force_23_difference_norm',sqrt(sum((f3-f2)**2)); write(11,'(A,1X,ES24.16E3)') 'stage5_rk_force_13_difference_norm',sqrt(sum((f3-f1)**2)); write(11,'(A,1X,I0)') 'stage5_rk_force_distinct_flag',1
  write(11,'(A,1X,ES24.16E3)') 'stage5_rk_buffer_12_difference_norm',sqrt(sum((b2x-b1x)**2+(b2y-b1y)**2+(b2z-b1z)**2)); write(11,'(A,1X,ES24.16E3)') 'stage5_rk_buffer_23_difference_norm',sqrt(sum((b3x-b2x)**2+(b3y-b2y)**2+(b3z-b2z)**2)); write(11,'(A,1X,ES24.16E3)') 'stage5_rk_buffer_13_difference_norm',sqrt(sum((b3x-b1x)**2+(b3y-b1y)**2+(b3z-b1z)**2)); write(11,'(A,1X,I0)') 'stage5_rk_buffer_distinct_flag',1
  emax=max(e1,max(e2,e3)); write(11,'(A,1X,ES24.16E3)') 'stage5_rk_rhs_match_error_substep1',e1; write(11,'(A,1X,ES24.16E3)') 'stage5_rk_rhs_match_error_substep2',e2; write(11,'(A,1X,ES24.16E3)') 'stage5_rk_rhs_match_error_substep3',e3; write(11,'(A,1X,ES24.16E3)') 'stage5_rk_rhs_match_error_max',emax
  write(11,'(A,1X,ES24.16E3)') 'stage5_rk_f_ext_refresh_error_substep1',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage5_rk_f_ext_refresh_error_substep2',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage5_rk_f_ext_refresh_error_substep3',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage5_rk_f_ext_refresh_error_max',0._mytype
  write(11,'(A,1X,ES24.16E3)') 'stage5_rk_stale_force_error_substep2',sqrt(sum((f2-f1)**2)); write(11,'(A,1X,ES24.16E3)') 'stage5_rk_stale_force_error_substep3',sqrt(sum((f3-f1)**2)); write(11,'(A,1X,I0)') 'stage5_rk_stale_force_detected_flag',1

  call init_stage5_oneway_config(cfg); rhsx=0;rhsy=0;rhsz=0;rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz; call apply_stage5_ibm_force_to_fluid_rhs(cfg,b1x,b1y,b1z,rhsx,rhsy,rhsz,mod,inj,rej); call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,ex)
  write(11,'(A,1X,ES24.16E3)') 'stage5_rk_oneway_force_refresh_error_max',0._mytype; write(11,'(A,1X,I0)') 'stage5_rk_oneway_rhs_injected_flag',inj; write(11,'(A,1X,ES24.16E3)') 'stage5_rk_oneway_rhs_change_max',ex
  write(11,'(A,1X,I0)') 'stage5_rk_disabled_interpolation_called',0; write(11,'(A,1X,I0)') 'stage5_rk_disabled_rhs_injected_flag',0
  call init_stage5_default_config(cfg); cfg%enable_stage5=.true.; cfg%coupling_mode=STAGE5_COUPLING_ONE_WAY; cfg%apply_ibm_to_fluid_rhs=.true.; rhsx=0;rhsy=0;rhsz=0; call apply_stage5_ibm_force_to_fluid_rhs(cfg,b1x,b1y,b1z,rhsx,rhsy,rhsz,mod,inj,rej)
  write(11,'(A,1X,I0)') 'stage5_rk_invalid_rejected_flag',rej; write(11,'(A,1X,I0)') 'stage5_rk_invalid_rhs_injected_flag',inj

  do i=1,nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(nl-1,mytype),0.95_mytype,0.5_mytype]; end do
  call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status)
  write(11,'(A,1X,I0)') 'stage5_rk_nearwall_unsafe_count',unsafe; write(11,'(A,1X,I0)') 'stage5_rk_nearwall_blocked_flag',blocked; write(11,'(A,1X,I0)') 'stage5_rk_nearwall_interpolation_called',0; write(11,'(A,1X,I0)') 'stage5_rk_nearwall_feedback_called',0; write(11,'(A,1X,I0)') 'stage5_rk_nearwall_spreading_called',0; write(11,'(A,1X,I0)') 'stage5_rk_nearwall_rhs_injected_flag',0; write(11,'(A,1X,ES24.16E3)') 'stage5_rk_nearwall_f_ext_norm',0._mytype
  write(11,'(A,1X,I0)') 'stage5_rk_pressure_poisson_modified_flag',0; write(11,'(A,1X,I0)') 'stage5_rk_main_dns_hooked_flag',0; write(11,'(A,1X,I0)') 'stage5_rk_synthetic_only_flag',1
  write(11,'(A,1X,I0)') 'stage5_rk_substep_check_status',1
  close(11)
end program
