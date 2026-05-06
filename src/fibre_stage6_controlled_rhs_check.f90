program fibre_stage6_controlled_rhs_check
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_controlled_rhs_hook
  implicit none
  type(stage6_config_t) :: cfg
  real(mytype) :: rhsx(8,6,5),rhsy(8,6,5),rhsz(8,6,5),rhsx0(8,6,5),rhsy0(8,6,5),rhsz0(8,6,5),rhsx2(8,6,5),rhsy2(8,6,5),rhsz2(8,6,5)
  real(mytype) :: fx(8,6,5),fy(8,6,5),fz(8,6,5),ux(8,6,5),uy(8,6,5),uz(8,6,5),ux0(8,6,5),uy0(8,6,5),uz0(8,6,5)
  real(mytype) :: chg,err,ex,ey,ez,bmax,e2,e4,sc,vchg
  integer :: i,j,k,hook,mod,inj,rej,pp,pm,pr
  open(11,file='stage6_outputs/fibre_stage6_controlled_rhs_check.dat',status='replace')
  do k=1,5;do j=1,6;do i=1,8
    rhsx(i,j,k)=0.1_mytype*i+0.01_mytype*j; rhsy(i,j,k)=-0.2_mytype*j+0.03_mytype*k; rhsz(i,j,k)=0.05_mytype*k-0.01_mytype*i
    fx(i,j,k)=sin(0.1_mytype*i+0.2_mytype*j+0.3_mytype*k); fy(i,j,k)=cos(0.2_mytype*i-0.1_mytype*k); fz(i,j,k)=0.1_mytype*sin(0.3_mytype*j)
    ux(i,j,k)=0.1_mytype*i; uy(i,j,k)=0.2_mytype*j; uz(i,j,k)=0.3_mytype*k
  end do;end do;end do
  ux0=ux;uy0=uy;uz0=uz; rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz
  call init_stage6_controlled_test_config(cfg); cfg%rho_fluid=2._mytype
  call apply_stage6_controlled_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej)
  call compute_stage6_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg); call compute_stage6_rhs_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,cfg%rho_fluid,rhsx,rhsy,rhsz,err,ex,ey,ez)
  bmax=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))
  write(11,'(A,1X,I0)') 'stage6_controlled_hook_called_flag',hook
  write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_buffer_max_abs',bmax
  write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_rhs_change_max',chg
  write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_expected_error',err
  write(11,'(A,1X,I0)') 'stage6_controlled_injected_flag',inj
  write(11,'(A,1X,I0)') 'stage6_controlled_modified_flag',mod
  write(11,'(A,1X,I0)') 'stage6_controlled_rejected_flag',rej
  write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_component_x_error',ex
  write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_component_y_error',ey
  write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_component_z_error',ez

  rhsx2=rhsx0;rhsy2=rhsy0;rhsz2=rhsz0; cfg%rho_fluid=2._mytype; call apply_stage6_controlled_rhs_hook(cfg,fx,fy,fz,rhsx2,rhsy2,rhsz2,hook,mod,inj,rej); call compute_stage6_rhs_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,2._mytype,rhsx2,rhsy2,rhsz2,e2,ex,ey,ez)
  rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0; cfg%rho_fluid=4._mytype; call apply_stage6_controlled_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); call compute_stage6_rhs_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,4._mytype,rhsx,rhsy,rhsz,e4,ex,ey,ez)
  sc=max(maxval(abs((rhsx-rhsx0)-0.5_mytype*(rhsx2-rhsx0))),max(maxval(abs((rhsy-rhsy0)-0.5_mytype*(rhsy2-rhsy0))),maxval(abs((rhsz-rhsz0)-0.5_mytype*(rhsz2-rhsz0)))))
  write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_rho2_expected_error',e2; write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_rho4_expected_error',e4; write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_rho_scaling_error',sc

  cfg%rho_fluid=2._mytype; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0; fx=0;fy=0;fz=0; call apply_stage6_controlled_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); call compute_stage6_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg)
  write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_zero_buffer_rhs_change_max',chg; write(11,'(A,1X,I0)') 'stage6_controlled_zero_buffer_injected_flag',inj; write(11,'(A,1X,I0)') 'stage6_controlled_zero_buffer_modified_flag',mod

  do k=1,5;do j=1,6;do i=1,8; fx(i,j,k)=sin(0.1_mytype*i+0.2_mytype*j+0.3_mytype*k); fy(i,j,k)=cos(0.2_mytype*i-0.1_mytype*k); fz(i,j,k)=0.1_mytype*sin(0.3_mytype*j); end do;end do;end do
  call init_stage6_default_config(cfg); rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0; call apply_stage6_controlled_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); call compute_stage6_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg)
  write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_default_rhs_change_max',chg; write(11,'(A,1X,I0)') 'stage6_controlled_default_injected_flag',inj; write(11,'(A,1X,I0)') 'stage6_controlled_default_modified_flag',mod

  call init_stage6_default_config(cfg); cfg%enable_stage6=.true.; cfg%enable_main_rhs_hook=.true.; cfg%allow_stage5_hook_in_main_path=.true.; cfg%rho_fluid=1._mytype; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage6_controlled_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); call compute_stage6_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg)
  write(11,'(A,1X,I0)') 'stage6_controlled_invalid_rejected_flag',rej; write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_invalid_rhs_change_max',chg; write(11,'(A,1X,I0)') 'stage6_controlled_invalid_injected_flag',inj

  call init_stage6_default_config(cfg); cfg%enable_stage6=.true.; cfg%enable_main_rhs_hook=.true.; cfg%production_two_way_enabled=.true.; cfg%allow_stage5_hook_in_main_path=.true.; cfg%rho_fluid=1._mytype; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage6_controlled_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); call compute_stage6_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg)
  write(11,'(A,1X,I0)') 'stage6_controlled_production_enabled_rejected_flag',rej; write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_production_enabled_rhs_change_max',chg; write(11,'(A,1X,I0)') 'stage6_controlled_production_enabled_injected_flag',inj

  call init_stage6_controlled_test_config(cfg); cfg%rho_fluid=0._mytype; rhsx=rhsx0;rhsy=rhsy0;rhsz=rhsz0
  call apply_stage6_controlled_rhs_hook(cfg,fx,fy,fz,rhsx,rhsy,rhsz,hook,mod,inj,rej); call compute_stage6_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,chg)
  write(11,'(A,1X,I0)') 'stage6_controlled_invalid_rho_rejected_flag',rej; write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_invalid_rho_rhs_change_max',chg; write(11,'(A,1X,I0)') 'stage6_controlled_invalid_rho_injected_flag',inj

  vchg=max(maxval(abs(ux-ux0)),max(maxval(abs(uy-uy0)),maxval(abs(uz-uz0))))
  write(11,'(A,1X,ES24.16E3)') 'stage6_controlled_velocity_change_max',vchg; write(11,'(A,1X,I0)') 'stage6_controlled_fluid_update_called_flag',0
  call stage6_controlled_pressure_status(pp,pm,pr)
  write(11,'(A,1X,I0)') 'stage6_controlled_pressure_poisson_modified_flag',pp; write(11,'(A,1X,I0)') 'stage6_controlled_projection_modified_flag',pm; write(11,'(A,1X,I0)') 'stage6_controlled_pressure_rhs_modified_flag',pr
  write(11,'(A,1X,I0)') 'stage6_controlled_rhs_hook_check_status',1
  close(11)
end program
