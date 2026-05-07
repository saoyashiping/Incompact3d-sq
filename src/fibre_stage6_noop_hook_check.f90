program fibre_stage6_noop_hook_check
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_main_noop_hook
  implicit none
  type(stage6_config_t) :: cfg
  real(mytype) :: rhsx(8,6,5),rhsy(8,6,5),rhsz(8,6,5),rhsx0(8,6,5),rhsy0(8,6,5),rhsz0(8,6,5)
  real(mytype) :: fx(8,6,5),fy(8,6,5),fz(8,6,5),ux(8,6,5),uy(8,6,5),uz(8,6,5),ux0(8,6,5),uy0(8,6,5),uz0(8,6,5)
  real(mytype) :: chg,vchg
  integer :: i,j,k,hook,mod,inj,rej,pp,pm,pr,safe,prod,ctl
  open(11,file='stage6_outputs/fibre_stage6_noop_hook_check.dat',status='replace')
  do k=1,5;do j=1,6;do i=1,8
    rhsx(i,j,k)=0.1_mytype*i+0.01_mytype*j; rhsy(i,j,k)=-0.2_mytype*j+0.03_mytype*k; rhsz(i,j,k)=0.05_mytype*k-0.01_mytype*i
    fx(i,j,k)=sin(0.1_mytype*i+0.2_mytype*j+0.3_mytype*k); fy(i,j,k)=cos(0.2_mytype*i-0.1_mytype*k); fz(i,j,k)=0.1_mytype*sin(0.3_mytype*j)
    ux(i,j,k)=0.2_mytype+0.01_mytype*i; uy(i,j,k)=0.03_mytype*j; uz(i,j,k)=0.05_mytype*k
  end do;end do;end do
  ux0=ux;uy0=uy;uz0=uz

  call init_stage6_default_config(cfg); rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz
  call apply_stage6_main_noop_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej)
  chg=max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
  write(11,'(A,1X,I0)') 'stage6_noop_default_hook_called_flag',hook
  write(11,'(A,1X,ES24.16E3)') 'stage6_noop_default_rhs_change_max',chg
  write(11,'(A,1X,I0)') 'stage6_noop_default_rhs_modified_flag',mod
  write(11,'(A,1X,I0)') 'stage6_noop_default_injected_flag',inj
  write(11,'(A,1X,I0)') 'stage6_noop_default_rejected_flag',rej

  call init_stage6_default_config(cfg); rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage6_main_noop_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej)
  chg=max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
  write(11,'(A,1X,ES24.16E3)') 'stage6_noop_oneway_rhs_change_max',chg
  write(11,'(A,1X,I0)') 'stage6_noop_oneway_rhs_modified_flag',mod
  write(11,'(A,1X,I0)') 'stage6_noop_oneway_injected_flag',inj

  call init_stage6_default_config(cfg); cfg%enable_stage6=.true.; cfg%enable_main_rhs_hook=.true.; cfg%allow_stage5_hook_in_main_path=.true.; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage6_main_noop_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); chg=max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
  write(11,'(A,1X,I0)') 'stage6_noop_invalid_rejected_flag',rej
  write(11,'(A,1X,ES24.16E3)') 'stage6_noop_invalid_rhs_change_max',chg
  write(11,'(A,1X,I0)') 'stage6_noop_invalid_injected_flag',inj
  write(11,'(A,1X,I0)') 'stage6_noop_invalid_rhs_modified_flag',mod

  call init_stage6_default_config(cfg); cfg%enable_stage6=.true.; cfg%enable_main_rhs_hook=.true.; cfg%production_two_way_enabled=.true.; cfg%allow_stage5_hook_in_main_path=.true.; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage6_main_noop_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); chg=max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
  write(11,'(A,1X,I0)') 'stage6_noop_production_enabled_rejected_flag',rej
  write(11,'(A,1X,ES24.16E3)') 'stage6_noop_production_enabled_rhs_change_max',chg
  write(11,'(A,1X,I0)') 'stage6_noop_production_enabled_injected_flag',inj

  call init_stage6_controlled_test_config(cfg); cfg%rho_fluid=0._mytype; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage6_main_noop_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); chg=max(maxval(abs(rhsx-rhsx0)),max(maxval(abs(rhsy-rhsy0)),maxval(abs(rhsz-rhsz0))))
  write(11,'(A,1X,I0)') 'stage6_noop_invalid_rho_rejected_flag',rej
  write(11,'(A,1X,ES24.16E3)') 'stage6_noop_invalid_rho_rhs_change_max',chg
  write(11,'(A,1X,I0)') 'stage6_noop_invalid_rho_injected_flag',inj

  vchg=max(maxval(abs(ux-ux0)),max(maxval(abs(uy-uy0)),maxval(abs(uz-uz0))))
  write(11,'(A,1X,ES24.16E3)') 'stage6_noop_velocity_change_max',vchg
  write(11,'(A,1X,I0)') 'stage6_noop_fluid_update_called_flag',0

  call stage6_noop_hook_pressure_status(pp,pm,pr)
  write(11,'(A,1X,I0)') 'stage6_noop_pressure_poisson_modified_flag',pp
  write(11,'(A,1X,I0)') 'stage6_noop_projection_modified_flag',pm
  write(11,'(A,1X,I0)') 'stage6_noop_pressure_rhs_modified_flag',pr

  call stage6_noop_hook_production_safety_status(safe,prod,ctl)
  write(11,'(A,1X,I0)') 'stage6_noop_default_main_dns_safe_flag',safe
  write(11,'(A,1X,I0)') 'stage6_noop_production_enabled_by_default_flag',prod
  write(11,'(A,1X,I0)') 'stage6_noop_controlled_test_only_flag',ctl
  write(11,'(A,1X,I0)') 'stage6_noop_hook_check_status',1
  close(11)
end program
