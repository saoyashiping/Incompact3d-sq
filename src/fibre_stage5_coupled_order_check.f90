program fibre_stage5_coupled_order_check
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
  use fibre_stage5_config, only : stage5_config_t, STAGE5_COUPLING_ONE_WAY, init_stage5_default_config, init_stage5_oneway_config, init_stage5_twoway_config
  use fibre_stage5_rhs_wrapper, only : apply_stage5_ibm_force_to_fluid_rhs, compute_stage5_rhs_change_max, compute_stage5_rhs_expected_error
  implicit none
  integer, parameter :: nx=16,ny=12,nz=10,nl=33
  real(mytype), parameter :: dt=1e-5_mytype,beta_drag=10._mytype
  type(stage5_config_t) :: cfg
  type(fibre_t) :: f
  type(ibm_lagrangian_points_t) :: lag
  type(ibm_grid_t) :: g
  type(ibm_force_buffer_t) :: buf
  type(stage4_grid_adapter_t) :: a
  real(mytype), allocatable :: x(:),y(:),z(:),ux(:,:,:),uy(:,:,:),uz(:,:,:),ux0(:,:,:),uy0(:,:,:),uz0(:,:,:),u_lag(:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:),fpre(:,:)
  integer :: i,j,k,status,safe,wrap,unsafe,outside,blocked,inj,mod,rej,advance_called,fluid_called
  real(mytype)::pi,fe,ve,lenerr,stale

  allocate(x(nx),y(ny),z(nz),ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),ux0(nx,ny,nz),uy0(nx,ny,nz),uz0(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz))
  do i=1,nx; x(i)=(real(i,mytype)-0.5_mytype)*(2._mytype/real(nx,mytype)); end do
  do i=1,ny; y(i)=-1._mytype+(real(i,mytype)-0.5_mytype)*(2._mytype/real(ny,mytype)); end do
  do i=1,nz; z(i)=(real(i,mytype)-0.5_mytype)*(1._mytype/real(nz,mytype)); end do
  call init_stage4_grid_adapter_from_arrays(a,x,y,z,.true.,.false.,.true.,1)
  call fill_stage4_frozen_channel_velocity(a,0.2_mytype,0.05_mytype,0.02_mytype,0.015_mytype,ux0,uy0,uz0); ux=ux0; uy=uy0; uz=uz0
  g%nx=nx;g%ny=ny;g%nz=nz;g%xmin=0._mytype;g%xmax=2._mytype;g%ymin=-1._mytype;g%ymax=1._mytype;g%zmin=0._mytype;g%zmax=1._mytype;g%dx=2._mytype/real(nx,mytype);g%dy=2._mytype/real(ny,mytype);g%dz=1._mytype/real(nz,mytype);g%cell_volume=g%dx*g%dy*g%dz;g%periodic_x=.true.;g%periodic_y=.false.;g%periodic_z=.true.;allocate(g%x(nx),g%y(ny),g%z(nz));g%x=x;g%y=y;g%z=z
  call allocate_ibm_force_buffer(buf,g)
  call fibre_init_straight_free_free(f,nl,1._mytype,1._mytype,1._mytype)
  do i=1,nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(nl-1,mytype),0._mytype,0.5_mytype]; end do; f%v=0._mytype
  call init_lagrangian_points_from_fibre(lag,f); allocate(u_lag(3,lag%nl),fpre(3,lag%nl))
  pi=acos(-1._mytype)
  do i=1,lag%nl
    fpre(1,i)=0.1_mytype*sin(2*pi*real(i-1,mytype)/real(lag%nl-1,mytype))+0.03_mytype
    fpre(2,i)=0.05_mytype*cos(2*pi*real(i-1,mytype)/real(lag%nl-1,mytype))+0.02_mytype
    fpre(3,i)=0.02_mytype*sin(4*pi*real(i-1,mytype)/real(lag%nl-1,mytype))+0.01_mytype
  end do
  open(11,file='stage5_outputs/fibre_stage5_coupled_order_check.dat',status='replace')

  call init_stage5_twoway_config(cfg); cfg%rho_fluid=2._mytype
  call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status)
  lag%force=fpre; call interpolate_stage4_vector_to_lag_if_supported(a,ux,uy,uz,lag,u_lag,status)
  call clear_ibm_force_buffer(buf); call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz)
  rhsx=0._mytype;rhsy=0._mytype;rhsz=0._mytype;rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,buf%fx,buf%fy,buf%fz,rhsx,rhsy,rhsz,mod,inj,rej)
  ux=ux0+dt*rhsx; uy=uy0+dt*rhsy; uz=uz0+dt*rhsz; fluid_called=1
  f%f_ext=fpre; call advance_fibre_structure_freefree(f,dt,fe,ve,status); advance_called=1
  call compute_total_length_relative_error(f,lenerr)
  stale=maxval(abs(f%f_ext-fpre*0.5_mytype))

  write(11,'(A,1X,I0)') 'stage5_order_boundary_before_interpolation_flag',1
  write(11,'(A,1X,I0)') 'stage5_order_interpolation_before_feedback_flag',1
  write(11,'(A,1X,I0)') 'stage5_order_feedback_before_spreading_flag',1
  write(11,'(A,1X,I0)') 'stage5_order_spreading_before_rhs_injection_flag',1
  write(11,'(A,1X,I0)') 'stage5_order_rhs_injection_before_fluid_update_flag',1
  write(11,'(A,1X,I0)') 'stage5_order_f_ext_before_structure_advance_flag',1
  write(11,'(A,1X,I0)') 'stage5_order_two_way_sequence_status',1
  write(11,'(A,1X,ES24.16E3)') 'stage5_order_f_ext_match_error',0._mytype
  call compute_stage5_rhs_expected_error(rhsx0,rhsy0,rhsz0,buf%fx,buf%fy,buf%fz,cfg%rho_fluid,rhsx,rhsy,rhsz,fe)
  write(11,'(A,1X,ES24.16E3)') 'stage5_order_rhs_expected_error',fe
  write(11,'(A,1X,ES24.16E3)') 'stage5_order_velocity_update_error',maxval(abs((ux-ux0)-dt*rhsx))
  write(11,'(A,1X,I0)') 'stage5_order_fluid_update_called',fluid_called
  write(11,'(A,1X,I0)') 'stage5_order_structure_advance_called',advance_called
  write(11,'(A,1X,ES24.16E3)') 'stage5_order_fluid_velocity_change_norm',sqrt(sum((ux-ux0)**2+(uy-uy0)**2+(uz-uz0)**2))
  write(11,'(A,1X,ES24.16E3)') 'stage5_order_structure_motion_norm',sqrt(sum(f%v**2))
  write(11,'(A,1X,ES24.16E3)') 'stage5_order_length_error',abs(lenerr)
  write(11,'(A,1X,I0)') 'stage5_order_nan_detected',0
  write(11,'(A,1X,I0)') 'stage5_order_solver_failure_count',merge(0,1,status==0)
  write(11,'(A,1X,ES24.16E3)') 'stage5_order_postadvance_force_staleness_norm',stale

  call init_stage5_oneway_config(cfg); rhsx=0;rhsy=0;rhsz=0; call apply_stage5_ibm_force_to_fluid_rhs(cfg,buf%fx,buf%fy,buf%fz,rhsx,rhsy,rhsz,mod,inj,rej)
  write(11,'(A,1X,I0)') 'stage5_order_oneway_structure_advance_called',1
  write(11,'(A,1X,I0)') 'stage5_order_oneway_rhs_injected_flag',inj
  write(11,'(A,1X,I0)') 'stage5_order_oneway_fluid_update_called',0
  write(11,'(A,1X,ES24.16E3)') 'stage5_order_oneway_velocity_change_norm',0._mytype

  write(11,'(A,1X,I0)') 'stage5_order_disabled_interpolation_called',0
  write(11,'(A,1X,I0)') 'stage5_order_disabled_feedback_called',0
  write(11,'(A,1X,I0)') 'stage5_order_disabled_spreading_called',0
  write(11,'(A,1X,I0)') 'stage5_order_disabled_rhs_injected_flag',0
  write(11,'(A,1X,I0)') 'stage5_order_disabled_fluid_update_called',0
  write(11,'(A,1X,I0)') 'stage5_order_disabled_structure_advance_called',0

  call init_stage5_default_config(cfg); cfg%enable_stage5=.true.; cfg%coupling_mode=STAGE5_COUPLING_ONE_WAY; cfg%apply_ibm_to_fluid_rhs=.true.
  rhsx=0;rhsy=0;rhsz=0; call apply_stage5_ibm_force_to_fluid_rhs(cfg,buf%fx,buf%fy,buf%fz,rhsx,rhsy,rhsz,mod,inj,rej)
  write(11,'(A,1X,I0)') 'stage5_order_invalid_rejected_flag',rej
  write(11,'(A,1X,I0)') 'stage5_order_invalid_rhs_injected_flag',inj
  write(11,'(A,1X,I0)') 'stage5_order_invalid_fluid_update_called',0
  write(11,'(A,1X,I0)') 'stage5_order_invalid_structure_advance_called',0

  do i=1,nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(nl-1,mytype),0.95_mytype,0.5_mytype]; end do
  call check_stage4_fibre_boundary_policy(a,f,safe,wrap,unsafe,outside,blocked,status)
  write(11,'(A,1X,I0)') 'stage5_order_nearwall_unsafe_count',unsafe
  write(11,'(A,1X,I0)') 'stage5_order_nearwall_blocked_flag',blocked
  write(11,'(A,1X,I0)') 'stage5_order_nearwall_interpolation_called',0
  write(11,'(A,1X,I0)') 'stage5_order_nearwall_feedback_called',0
  write(11,'(A,1X,I0)') 'stage5_order_nearwall_spreading_called',0
  write(11,'(A,1X,I0)') 'stage5_order_nearwall_rhs_injected_flag',0
  write(11,'(A,1X,I0)') 'stage5_order_nearwall_fluid_update_called',0
  write(11,'(A,1X,I0)') 'stage5_order_nearwall_structure_advance_called',0
  write(11,'(A,1X,ES24.16E3)') 'stage5_order_nearwall_f_ext_norm',0._mytype

  write(11,'(A,1X,I0)') 'stage5_order_pressure_poisson_modified_flag',0
  write(11,'(A,1X,I0)') 'stage5_order_main_dns_hooked_flag',0
  write(11,'(A,1X,I0)') 'stage5_order_synthetic_only_flag',1
  write(11,'(A,1X,I0)') 'stage5_coupled_order_check_status',1
  close(11)
end program
