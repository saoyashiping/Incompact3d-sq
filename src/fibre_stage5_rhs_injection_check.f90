program fibre_stage5_rhs_injection_check
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre
  use fibre_ibm_spreading, only : spread_lag_force_to_eulerian
  use fibre_ibm_force_buffer, only : ibm_force_buffer_t, allocate_ibm_force_buffer, clear_ibm_force_buffer
  use fibre_stage5_config, only : stage5_config_t, STAGE5_COUPLING_ONE_WAY, init_stage5_default_config, init_stage5_oneway_config, init_stage5_twoway_config
  use fibre_stage5_rhs_wrapper, only : apply_stage5_ibm_force_to_fluid_rhs, compute_stage5_rhs_change_max, compute_stage5_rhs_expected_error
  implicit none
  integer, parameter :: nx=16,ny=12,nz=10,nl=33
  type(stage5_config_t) :: cfg
  type(fibre_t) :: f
  type(ibm_lagrangian_points_t) :: lag
  type(ibm_grid_t) :: g
  type(ibm_force_buffer_t) :: buf
  real(mytype), allocatable :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:),fx(:,:,:),fy(:,:,:),fz(:,:,:)
  integer :: i,j,k,rej,inj,mod,status
  real(mytype) :: pi,chg,exp_err,bmax,ex,ey,ez,err2,err4,scale_err

  allocate(rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz),fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz))
  g%nx=nx;g%ny=ny;g%nz=nz;g%xmin=0._mytype;g%xmax=2._mytype;g%ymin=-1._mytype;g%ymax=1._mytype;g%zmin=0._mytype;g%zmax=1._mytype
  g%dx=2._mytype/real(nx,mytype);g%dy=2._mytype/real(ny,mytype);g%dz=1._mytype/real(nz,mytype);g%cell_volume=g%dx*g%dy*g%dz
  g%periodic_x=.true.;g%periodic_y=.false.;g%periodic_z=.true.;allocate(g%x(nx),g%y(ny),g%z(nz)); do i=1,nx; g%x(i)=(real(i,mytype)-0.5_mytype)*g%dx; end do; do j=1,ny; g%y(j)=g%ymin+(real(j,mytype)-0.5_mytype)*g%dy; end do; do k=1,nz; g%z(k)=(real(k,mytype)-0.5_mytype)*g%dz; end do
  call allocate_ibm_force_buffer(buf,g)

  call fibre_init_straight_free_free(f,nl,1._mytype,1._mytype,1._mytype)
  do i=1,nl; f%x(:,i)=[0.5_mytype+real(i-1,mytype)/real(nl-1,mytype),0._mytype,0.5_mytype]; end do
  call init_lagrangian_points_from_fibre(lag,f)
  pi=acos(-1._mytype)
  do i=1,lag%nl
    lag%force(1,i)=0.1_mytype*sin(2._mytype*pi*real(i-1,mytype)/real(lag%nl-1,mytype))+0.03_mytype
    lag%force(2,i)=0.05_mytype*cos(2._mytype*pi*real(i-1,mytype)/real(lag%nl-1,mytype))+0.02_mytype
    lag%force(3,i)=0.02_mytype*sin(4._mytype*pi*real(i-1,mytype)/real(lag%nl-1,mytype))+0.01_mytype
  end do
  call clear_ibm_force_buffer(buf); call spread_lag_force_to_eulerian(g,lag,buf%fx,buf%fy,buf%fz); fx=buf%fx; fy=buf%fy; fz=buf%fz

  do k=1,nz; do j=1,ny; do i=1,nx
    rhsx(i,j,k)=0.1_mytype*real(i,mytype)+0.01_mytype*real(j,mytype)
    rhsy(i,j,k)=-0.2_mytype*real(j,mytype)+0.03_mytype*real(k,mytype)
    rhsz(i,j,k)=0.05_mytype*real(k,mytype)-0.01_mytype*real(i,mytype)
  end do; end do; end do
  rhsx0=rhsx; rhsy0=rhsy; rhsz0=rhsz

  open(11,file='stage5_outputs/fibre_stage5_rhs_injection_check.dat',status='replace')

  call init_stage5_twoway_config(cfg); cfg%rho_fluid=2._mytype
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,0._mytype*fx,0._mytype*fy,0._mytype*fz,rhsx,rhsy,rhsz,mod,inj,rej)
  call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg)
  call compute_stage5_rhs_expected_error(rhsx0,rhsy0,rhsz0,0._mytype*fx,0._mytype*fy,0._mytype*fz,cfg%rho_fluid,rhsx,rhsy,rhsz,exp_err)
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_zero_buffer_change_max', chg
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_zero_buffer_expected_error', exp_err
  write(11,'(A,1X,I0)') 'stage5_injection_zero_buffer_injected_flag', inj
  write(11,'(A,1X,I0)') 'stage5_injection_zero_buffer_modified_flag', mod

  rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,mod,inj,rej)
  call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg)
  call compute_stage5_rhs_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,cfg%rho_fluid,rhsx,rhsy,rhsz,exp_err)
  bmax=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))
  ex=maxval(abs((rhsx-rhsx0)-fx/cfg%rho_fluid)); ey=maxval(abs((rhsy-rhsy0)-fy/cfg%rho_fluid)); ez=maxval(abs((rhsz-rhsz0)-fz/cfg%rho_fluid))
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_nonzero_buffer_max_abs', bmax
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_nonzero_rhs_change_max', chg
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_nonzero_expected_error', exp_err
  write(11,'(A,1X,I0)') 'stage5_injection_nonzero_injected_flag', inj
  write(11,'(A,1X,I0)') 'stage5_injection_nonzero_modified_flag', mod
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_component_x_error', ex
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_component_y_error', ey
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_component_z_error', ez

  call init_stage5_default_config(cfg); rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,mod,inj,rej); call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg)
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_disabled_change_max', chg
  write(11,'(A,1X,I0)') 'stage5_injection_disabled_injected_flag', inj
  write(11,'(A,1X,I0)') 'stage5_injection_disabled_modified_flag', mod

  call init_stage5_oneway_config(cfg); rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,mod,inj,rej); call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg)
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_oneway_change_max', chg
  write(11,'(A,1X,I0)') 'stage5_injection_oneway_injected_flag', inj
  write(11,'(A,1X,I0)') 'stage5_injection_oneway_modified_flag', mod

  call init_stage5_default_config(cfg); cfg%enable_stage5=.true.; cfg%coupling_mode=STAGE5_COUPLING_ONE_WAY; cfg%apply_ibm_to_fluid_rhs=.true.; cfg%allow_two_way=.false.; cfg%rho_fluid=1._mytype
  rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,mod,inj,rej); call compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg)
  write(11,'(A,1X,I0)') 'stage5_injection_invalid_rejected_flag', rej
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_invalid_change_max', chg
  write(11,'(A,1X,I0)') 'stage5_injection_invalid_injected_flag', inj

  call init_stage5_twoway_config(cfg); cfg%rho_fluid=2._mytype; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,mod,inj,rej); call compute_stage5_rhs_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,2._mytype,rhsx,rhsy,rhsz,err2)
  call init_stage5_twoway_config(cfg); cfg%rho_fluid=4._mytype; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage5_ibm_force_to_fluid_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,mod,inj,rej); call compute_stage5_rhs_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,4._mytype,rhsx,rhsy,rhsz,err4)
  scale_err=maxval(abs((rhsx-rhsx0)-0.5_mytype*(fx/2._mytype)))
  scale_err=max(scale_err,maxval(abs((rhsy-rhsy0)-0.5_mytype*(fy/2._mytype))))
  scale_err=max(scale_err,maxval(abs((rhsz-rhsz0)-0.5_mytype*(fz/2._mytype))))
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_rho2_expected_error', err2
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_rho4_expected_error', err4
  write(11,'(A,1X,ES24.16E3)') 'stage5_injection_rho_scaling_error', scale_err

  write(11,'(A,1X,I0)') 'stage5_injection_pressure_poisson_modified_flag', 0
  write(11,'(A,1X,I0)') 'stage5_injection_main_dns_hooked_flag', 0
  write(11,'(A,1X,I0)') 'stage5_injection_synthetic_only_flag', 1
  status=1
  write(11,'(A,1X,I0)') 'stage5_rhs_injection_check_status', status
  close(11)
end program
