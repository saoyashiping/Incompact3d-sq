program fibre_stage5_momentum_exchange_check
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre
  use fibre_ibm_spreading, only : spread_lag_force_to_eulerian
  use fibre_ibm_force_buffer, only : ibm_force_buffer_t, allocate_ibm_force_buffer, clear_ibm_force_buffer
  use fibre_stage5_config, only : stage5_config_t, STAGE5_COUPLING_ONE_WAY, init_stage5_default_config, init_stage5_oneway_config, init_stage5_twoway_config
  use fibre_stage5_rhs_wrapper, only : apply_stage5_ibm_force_to_fluid_rhs, compute_stage5_rhs_change_max
  use fibre_stage5_momentum_exchange, only : compute_stage5_fluid_momentum_change, compute_stage5_eulerian_force_impulse, compute_stage5_lagrangian_force_impulse
  implicit none
  integer, parameter :: nx=16,ny=12,nz=10,nl=33
  real(mytype), parameter :: dt=1e-5_mytype
  type(stage5_config_t) :: cfg
  type(fibre_t) :: f
  type(ibm_lagrangian_points_t) :: lag
  type(ibm_grid_t) :: g
  type(ibm_force_buffer_t) :: buf
  real(mytype), allocatable :: ux0(:,:,:),uy0(:,:,:),uz0(:,:,:),ux1(:,:,:),uy1(:,:,:),uz1(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),fx(:,:,:),fy(:,:,:),fz(:,:,:),rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:),fstruct(:,:)
  real(mytype) :: dpf(3),dpe(3),dps(3),norm0,bmax,chg,acterr,ferr,frel,ex,ey,ez,vchg,zerr
  integer :: i,j,k,inj,mod,rej,status
  real(mytype) :: pi

  allocate(ux0(nx,ny,nz),uy0(nx,ny,nz),uz0(nx,ny,nz),ux1(nx,ny,nz),uy1(nx,ny,nz),uz1(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz))
  g%nx=nx;g%ny=ny;g%nz=nz;g%dx=2._mytype/real(nx,mytype);g%dy=2._mytype/real(ny,mytype);g%dz=1._mytype/real(nz,mytype);g%cell_volume=g%dx*g%dy*g%dz
  g%xmin=0._mytype;g%xmax=2._mytype;g%ymin=-1._mytype;g%ymax=1._mytype;g%zmin=0._mytype;g%zmax=1._mytype;g%periodic_x=.true.;g%periodic_y=.false.;g%periodic_z=.true.; allocate(g%x(nx),g%y(ny),g%z(nz)); do i=1,nx; g%x(i)=(real(i,mytype)-0.5_mytype)*g%dx; end do; do j=1,ny; g%y(j)=g%ymin+(real(j,mytype)-0.5_mytype)*g%dy; end do; do k=1,nz; g%z(k)=(real(k,mytype)-0.5_mytype)*g%dz; end do
  call allocate_ibm_force_buffer(buf,g)

  do k=1,nz; do j=1,ny; do i=1,nx
    ux0(i,j,k)=0.1_mytype*real(i,mytype)+0.01_mytype*real(j,mytype)
    uy0(i,j,k)=-0.2_mytype*real(j,mytype)+0.03_mytype*real(k,mytype)
    uz0(i,j,k)=0.05_mytype*real(k,mytype)-0.01_mytype*real(i,mytype)
  end do; end do; end do
  ux1=ux0; uy1=uy0; uz1=uz0

  call fibre_init_straight_free_free(f,nl,1._mytype,1._mytype,1._mytype)
  do i=1,nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(nl-1,mytype),0._mytype,0.5_mytype]; end do
  call init_lagrangian_points_from_fibre(lag,f); allocate(fstruct(3,lag%nl))

  open(11,file='stage5_outputs/fibre_stage5_momentum_exchange_check.dat',status='replace')

  lag%force=0._mytype; fstruct=0._mytype; call clear_ibm_force_buffer(buf); call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz)
  call compute_stage5_eulerian_force_impulse(dt,g%cell_volume,buf%fx,buf%fy,buf%fz,dpe); call compute_stage5_lagrangian_force_impulse(dt,fstruct,lag%weight,dps)
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_zero_fluid_impulse_norm', sqrt(sum(dpe**2))
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_zero_structure_impulse_norm', sqrt(sum(dps**2))
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_zero_exchange_error', sqrt(sum((dpe+dps)**2))

  pi=acos(-1._mytype)
  do i=1,lag%nl
    lag%force(1,i)=0.1_mytype*sin(2._mytype*pi*real(i-1,mytype)/real(lag%nl-1,mytype))+0.03_mytype
    lag%force(2,i)=0.05_mytype*cos(2._mytype*pi*real(i-1,mytype)/real(lag%nl-1,mytype))+0.02_mytype
    lag%force(3,i)=0.02_mytype*sin(4._mytype*pi*real(i-1,mytype)/real(lag%nl-1,mytype))+0.01_mytype
  end do
  fstruct=-lag%force
  call clear_ibm_force_buffer(buf); call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz); fx=buf%fx;fy=buf%fy;fz=buf%fz
  bmax=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))
  rhsx=0._mytype; rhsy=0._mytype; rhsz=0._mytype; rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz
  call init_stage5_twoway_config(cfg); cfg%rho_fluid=2._mytype
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,mod,inj,rej)
  call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg)
  ux1=ux0+dt*rhsx; uy1=uy0+dt*rhsy; uz1=uz0+dt*rhsz
  call compute_stage5_fluid_momentum_change(cfg%rho_fluid,g%cell_volume,ux0,uy0,uz0,ux1,uy1,uz1,dpf)
  call compute_stage5_eulerian_force_impulse(dt,g%cell_volume,fx,fy,fz,dpe)
  call compute_stage5_lagrangian_force_impulse(dt,fstruct,lag%weight,dps)
  acterr=maxval(abs(fstruct+lag%force))
  ferr=sqrt(sum((dpf-dpe)**2)); frel=ferr/max(sqrt(sum(dpe**2)),tiny(1._mytype))
  zerr=sqrt(sum((dpf+dps)**2)); norm0=max(sqrt(sum(dps**2)),tiny(1._mytype))
  ex=abs(dpf(1)+dps(1)); ey=abs(dpf(2)+dps(2)); ez=abs(dpf(3)+dps(3))
  vchg=sqrt(sum((ux1-ux0)**2 + (uy1-uy0)**2 + (uz1-uz0)**2))

  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_nonzero_force_buffer_max_abs', bmax
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_nonzero_rhs_change_max', chg
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_fluid_velocity_change_norm', vchg
  write(11,'(A,1X,I0)') 'stage5_momentum_rhs_injected_flag', inj
  write(11,'(A,1X,I0)') 'stage5_momentum_rhs_modified_flag', mod
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_action_reaction_error', acterr
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_fluid_impulse_from_velocity_norm', sqrt(sum(dpf**2))
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_eulerian_force_impulse_norm', sqrt(sum(dpe**2))
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_fluid_impulse_error', ferr
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_fluid_impulse_relative_error', frel
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_structure_impulse_norm', sqrt(sum(dps**2))
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_exchange_abs_error', zerr
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_exchange_relative_error', zerr/norm0
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_component_x_error', ex
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_component_y_error', ey
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_component_z_error', ez

  call init_stage5_default_config(cfg); rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0; call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,mod,inj,rej)
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_disabled_velocity_change_norm', sqrt(sum((rhsx-rhsx0)**2 + (rhsy-rhsy0)**2 + (rhsz-rhsz0)**2))
  write(11,'(A,1X,I0)') 'stage5_momentum_disabled_injected_flag', inj
  call init_stage5_oneway_config(cfg); rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0; call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,mod,inj,rej)
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_oneway_velocity_change_norm', sqrt(sum((rhsx-rhsx0)**2 + (rhsy-rhsy0)**2 + (rhsz-rhsz0)**2))
  write(11,'(A,1X,I0)') 'stage5_momentum_oneway_injected_flag', inj

  call init_stage5_default_config(cfg); cfg%enable_stage5=.true.; cfg%coupling_mode=STAGE5_COUPLING_ONE_WAY; cfg%apply_ibm_to_fluid_rhs=.true.
  rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0; call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,mod,inj,rej)
  write(11,'(A,1X,I0)') 'stage5_momentum_invalid_rejected_flag', rej
  write(11,'(A,1X,ES24.16E3)') 'stage5_momentum_invalid_velocity_change_norm', sqrt(sum((rhsx-rhsx0)**2 + (rhsy-rhsy0)**2 + (rhsz-rhsz0)**2))
  write(11,'(A,1X,I0)') 'stage5_momentum_invalid_injected_flag', inj

  write(11,'(A,1X,I0)') 'stage5_momentum_pressure_poisson_modified_flag', 0
  write(11,'(A,1X,I0)') 'stage5_momentum_main_dns_hooked_flag', 0
  write(11,'(A,1X,I0)') 'stage5_momentum_synthetic_only_flag', 1
  status=1
  write(11,'(A,1X,I0)') 'stage5_momentum_exchange_check_status', status
  close(11)
end program
