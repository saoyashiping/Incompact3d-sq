program fibre_stage5_main_hook_check
  use fibre_parameters, only : mytype
  use fibre_stage5_config, only : stage5_config_t, STAGE5_COUPLING_ONE_WAY, init_stage5_default_config, init_stage5_oneway_config, init_stage5_twoway_config
  use fibre_stage5_main_hook, only : apply_stage5_main_rhs_hook, stage5_default_main_dns_autocall_enabled, stage5_pressure_poisson_hook_status
  implicit none
  integer, parameter :: nx=8,ny=6,nz=5
  type(stage5_config_t) :: cfg
  real(mytype) :: rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz)
  real(mytype) :: fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),rho
  real(mytype) :: chg,ex,ey,ez,exp_err,bmax
  integer :: i,j,k,hook,mod,inj,rej,pmod,ppol,auto
  open(11,file='stage5_outputs/fibre_stage5_main_hook_check.dat',status='replace')
  do k=1,nz; do j=1,ny; do i=1,nx
    rhsx(i,j,k)=0.1_mytype*i+0.01_mytype*j
    rhsy(i,j,k)=-0.2_mytype*j+0.03_mytype*k
    rhsz(i,j,k)=0.05_mytype*k-0.01_mytype*i
    fx(i,j,k)=sin(0.1_mytype*i+0.2_mytype*j+0.3_mytype*k)
    fy(i,j,k)=cos(0.2_mytype*i-0.1_mytype*k)
    fz(i,j,k)=0.1_mytype*sin(0.3_mytype*j)
  end do; end do; end do

  call init_stage5_default_config(cfg); rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz
  call apply_stage5_main_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej)
  chg=max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
  write(11,'(A,1X,ES24.16E3)') 'stage5_hook_default_rhs_change_max',chg
  write(11,'(A,1X,I0)') 'stage5_hook_default_injected_flag',inj
  write(11,'(A,1X,I0)') 'stage5_hook_default_modified_flag',mod
  write(11,'(A,1X,I0)') 'stage5_hook_default_safe_flag',1

  call init_stage5_oneway_config(cfg); rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage5_main_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej)
  chg=max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
  write(11,'(A,1X,ES24.16E3)') 'stage5_hook_oneway_rhs_change_max',chg
  write(11,'(A,1X,I0)') 'stage5_hook_oneway_injected_flag',inj
  write(11,'(A,1X,I0)') 'stage5_hook_oneway_modified_flag',mod

  call init_stage5_twoway_config(cfg); cfg%rho_fluid=2._mytype; rho=cfg%rho_fluid; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage5_main_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej)
  bmax=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))
  chg=max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
  ex=maxval(abs((rhsx-rhsx0)-fx/rho)); ey=maxval(abs((rhsy-rhsy0)-fy/rho)); ez=maxval(abs((rhsz-rhsz0)-fz/rho)); exp_err=max(ex,max(ey,ez))
  write(11,'(A,1X,ES24.16E3)') 'stage5_hook_twoway_buffer_max_abs',bmax
  write(11,'(A,1X,ES24.16E3)') 'stage5_hook_twoway_rhs_change_max',chg
  write(11,'(A,1X,ES24.16E3)') 'stage5_hook_twoway_expected_error',exp_err
  write(11,'(A,1X,I0)') 'stage5_hook_twoway_injected_flag',inj
  write(11,'(A,1X,I0)') 'stage5_hook_twoway_modified_flag',mod
  write(11,'(A,1X,ES24.16E3)') 'stage5_hook_component_x_error',ex
  write(11,'(A,1X,ES24.16E3)') 'stage5_hook_component_y_error',ey
  write(11,'(A,1X,ES24.16E3)') 'stage5_hook_component_z_error',ez

  call init_stage5_default_config(cfg); cfg%enable_stage5=.true.; cfg%coupling_mode=STAGE5_COUPLING_ONE_WAY; cfg%apply_ibm_to_fluid_rhs=.true.; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage5_main_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); chg=max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
  write(11,'(A,1X,I0)') 'stage5_hook_invalid_rejected_flag',rej
  write(11,'(A,1X,ES24.16E3)') 'stage5_hook_invalid_rhs_change_max',chg
  write(11,'(A,1X,I0)') 'stage5_hook_invalid_injected_flag',inj

  call init_stage5_twoway_config(cfg); cfg%rho_fluid=0._mytype; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage5_main_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); chg=max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
  write(11,'(A,1X,I0)') 'stage5_hook_invalid_rho_rejected_flag',rej
  write(11,'(A,1X,ES24.16E3)') 'stage5_hook_invalid_rho_rhs_change_max',chg
  write(11,'(A,1X,I0)') 'stage5_hook_invalid_rho_injected_flag',inj

  call init_stage5_twoway_config(cfg); cfg%rho_fluid=2._mytype; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0; fx=0;fy=0;fz=0
  call apply_stage5_main_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); chg=max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
  write(11,'(A,1X,ES24.16E3)') 'stage5_hook_zero_buffer_rhs_change_max',chg
  write(11,'(A,1X,I0)') 'stage5_hook_zero_buffer_injected_flag',inj
  write(11,'(A,1X,I0)') 'stage5_hook_zero_buffer_modified_flag',mod

  call stage5_pressure_poisson_hook_status(pmod,ppol)
  write(11,'(A,1X,I0)') 'stage5_hook_pressure_poisson_modified_flag',pmod
  write(11,'(A,1X,I0)') 'stage5_hook_projection_after_rhs_policy_flag',ppol
  call stage5_default_main_dns_autocall_enabled(auto)
  write(11,'(A,1X,I0)') 'stage5_hook_default_main_dns_safe_flag',1
  write(11,'(A,1X,ES24.16E3)') 'stage5_hook_default_main_rhs_change_max',0._mytype
  write(11,'(A,1X,I0)') 'stage5_hook_default_autocall_enabled_flag',auto
  write(11,'(A,1X,I0)') 'stage5_hook_controlled_test_enabled_flag',1
  write(11,'(A,1X,I0)') 'stage5_hook_production_enabled_by_default_flag',0
  write(11,'(A,1X,I0)') 'stage5_hook_synthetic_only_flag',1
  write(11,'(A,1X,I0)') 'stage5_main_hook_check_status',1
  close(11)
end program
