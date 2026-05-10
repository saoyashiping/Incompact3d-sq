program fibre_stage6_projection_order_check
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_projection_order
  use fibre_stage6_controlled_rhs_hook, only : compute_stage6_rhs_change_max, compute_stage6_rhs_expected_error
  implicit none
  integer, parameter :: nx=10,ny=8,nz=6
  real(mytype), parameter :: dt=1e-5_mytype
  type(stage6_config_t) :: cfg
  real(mytype) :: rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz)
  real(mytype) :: rhsx2(nx,ny,nz),rhsy2(nx,ny,nz),rhsz2(nx,ny,nz),ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),uxs(nx,ny,nz),uys(nx,ny,nz),uzs(nx,ny,nz)
  real(mytype) :: fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),prhs(nx,ny,nz),prhs0(nx,ny,nz)
  real(mytype) :: chg,err,ex,ey,ez,ustarerr,ustarchg,bmax,prhsdiff
  integer :: i,j,k,hook,mod,inj,rej,bf,af,ppd,pol,pvm,pff,pmod,rpc
  open(11,file='stage6_outputs/fibre_stage6_projection_order_check.dat',status='replace')
  do k=1,nz;do j=1,ny;do i=1,nx
    rhsx(i,j,k)=0.1_mytype*i+0.01_mytype*j+0.001_mytype*k; rhsy(i,j,k)=-0.2_mytype*j+0.03_mytype*k+0.002_mytype*i; rhsz(i,j,k)=0.05_mytype*k-0.01_mytype*i+0.004_mytype*j
    ux(i,j,k)=0.2_mytype+0.01_mytype*i+0.002_mytype*j; uy(i,j,k)=0.01_mytype*sin(0.2_mytype*i+0.1_mytype*k); uz(i,j,k)=0.01_mytype*cos(0.1_mytype*j+0.3_mytype*k)
    prhs(i,j,k)=0.001_mytype*cos(0.2_mytype*i-0.1_mytype*j+0.05_mytype*k)
    fx(i,j,k)=sin(0.1_mytype*i+0.2_mytype*j+0.3_mytype*k); fy(i,j,k)=cos(0.2_mytype*i-0.1_mytype*k); fz(i,j,k)=0.1_mytype*sin(0.3_mytype*j)
  end do;end do;end do
  rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz; prhs0=prhs
  call init_stage6_controlled_test_config(cfg); cfg%rho_fluid=2._mytype
  call apply_stage6_projection_order_rhs(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej,bf,af,ppd)
  call compute_stage6_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg); call compute_stage6_rhs_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,cfg%rho_fluid,rhsx,rhsy,rhsz,err,ex,ey,ez)
  call compute_stage6_intermediate_velocity(dt,ux,uy,uz,rhsx,rhsy,rhsz,uxs,uys,uzs); call compute_stage6_ustar_expected_error(dt,ux,uy,uz,rhsx,rhsy,rhsz,uxs,uys,uzs,ustarerr)
  ustarchg=sqrt(sum((uxs-ux)**2+(uys-uy)**2+(uzs-uz)**2)); bmax=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))
  write(11,'(A,1X,I0)') 'stage6_projection_rhs_before_projection_flag',bf
  write(11,'(A,1X,I0)') 'stage6_projection_rhs_after_projection_flag',af
  call stage6_projection_policy_status(pol,pvm,pff,pmod,rpc)
  write(11,'(A,1X,I0)') 'stage6_projection_after_rhs_policy_flag',pol
  write(11,'(A,1X,I0)') 'stage6_projection_order_status',1
  write(11,'(A,1X,I0)') 'stage6_projection_momentum_rhs_modified_flag',mod
  write(11,'(A,1X,I0)') 'stage6_projection_pressure_rhs_modified_flag',0
  write(11,'(A,1X,I0)') 'stage6_projection_pressure_poisson_direct_modify_flag',ppd
  write(11,'(A,1X,ES24.16E3)') 'stage6_projection_rhs_expected_error',err
  write(11,'(A,1X,ES24.16E3)') 'stage6_projection_component_x_error',ex
  write(11,'(A,1X,ES24.16E3)') 'stage6_projection_component_y_error',ey
  write(11,'(A,1X,ES24.16E3)') 'stage6_projection_component_z_error',ez
  write(11,'(A,1X,ES24.16E3)') 'stage6_projection_ustar_expected_error',ustarerr
  write(11,'(A,1X,ES24.16E3)') 'stage6_projection_ustar_change_norm',ustarchg
  write(11,'(A,1X,I0)') 'stage6_projection_dt_positive_flag',merge(1,0,dt>0._mytype)
  write(11,'(A,1X,I0)') 'stage6_projection_post_projection_velocity_modified_flag',pvm
  write(11,'(A,1X,I0)') 'stage6_projection_post_projection_direct_forcing_forbidden_flag',pff

  call init_stage6_default_config(cfg); rhsx2=rhsx0;rhsy2=rhsy0;rhsz2=rhsz0
  call apply_stage6_projection_order_rhs(cfg,fx,fy,fz,rhsx2,rhsy2,rhsz2,hook,mod,inj,rej,bf,af,ppd); call compute_stage6_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx2,rhsy2,rhsz2,chg)
  call compute_stage6_intermediate_velocity(dt,ux,uy,uz,rhsx0,rhsy0,rhsz0,uxs,uys,uzs); call compute_stage6_intermediate_velocity(dt,ux,uy,uz,rhsx2,rhsy2,rhsz2,rhsx,rhsy,rhsz); call compute_stage6_ustar_expected_error(1._mytype,uxs,uys,uzs,rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,ustarerr)
  write(11,'(A,1X,ES24.16E3)') 'stage6_projection_default_rhs_change_max',chg
  write(11,'(A,1X,ES24.16E3)') 'stage6_projection_default_ustar_change_max',max(maxval(abs(rhsx-uxs)),max(maxval(abs(rhsy-uys)),maxval(abs(rhsz-uzs))))
  write(11,'(A,1X,I0)') 'stage6_projection_default_injected_flag',inj

  call init_stage6_default_config(cfg); cfg%enable_stage6=.true.; cfg%enable_main_rhs_hook=.true.; cfg%allow_stage5_hook_in_main_path=.true.; cfg%rho_fluid=1._mytype; rhsx2=rhsx0;rhsy2=rhsy0;rhsz2=rhsz0
  call apply_stage6_projection_order_rhs(cfg,fx,fy,fz,rhsx2,rhsy2,rhsz2,hook,mod,inj,rej,bf,af,ppd); call compute_stage6_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx2,rhsy2,rhsz2,chg)
  write(11,'(A,1X,I0)') 'stage6_projection_invalid_rejected_flag',rej; write(11,'(A,1X,ES24.16E3)') 'stage6_projection_invalid_rhs_change_max',chg; write(11,'(A,1X,I0)') 'stage6_projection_invalid_injected_flag',inj

  call init_stage6_default_config(cfg); cfg%enable_stage6=.true.; cfg%enable_main_rhs_hook=.true.; cfg%production_two_way_enabled=.true.; cfg%allow_stage5_hook_in_main_path=.true.; cfg%rho_fluid=1._mytype; rhsx2=rhsx0;rhsy2=rhsy0;rhsz2=rhsz0
  call apply_stage6_projection_order_rhs(cfg,fx,fy,fz,rhsx2,rhsy2,rhsz2,hook,mod,inj,rej,bf,af,ppd); call compute_stage6_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx2,rhsy2,rhsz2,chg)
  write(11,'(A,1X,I0)') 'stage6_projection_production_enabled_rejected_flag',rej; write(11,'(A,1X,ES24.16E3)') 'stage6_projection_production_enabled_rhs_change_max',chg; write(11,'(A,1X,I0)') 'stage6_projection_production_enabled_injected_flag',inj

  prhsdiff=maxval(abs(prhs-prhs0))
  write(11,'(A,1X,I0)') 'stage6_projection_pressure_poisson_modified_flag',0
  write(11,'(A,1X,ES24.16E3)') 'stage6_projection_pressure_rhs_diff_max',prhsdiff
  write(11,'(A,1X,I0)') 'stage6_projection_projection_modified_flag',pmod
  write(11,'(A,1X,I0)') 'stage6_projection_real_projection_called_flag',rpc
  write(11,'(A,1X,I0)') 'stage6_projection_order_check_status',1
  close(11)
end program
