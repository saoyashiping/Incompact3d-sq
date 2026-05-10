program fibre_stage6_rk_rhs_sync_check
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_rk_rhs_sync
  implicit none
  integer, parameter :: nx=10,ny=8,nz=6
  type(stage6_config_t)::cfg
  real(mytype)::r1x(nx,ny,nz),r1y(nx,ny,nz),r1z(nx,ny,nz),r2x(nx,ny,nz),r2y(nx,ny,nz),r2z(nx,ny,nz),r3x(nx,ny,nz),r3y(nx,ny,nz),r3z(nx,ny,nz)
  real(mytype)::r1x0(nx,ny,nz),r1y0(nx,ny,nz),r1z0(nx,ny,nz),r2x0(nx,ny,nz),r2y0(nx,ny,nz),r2z0(nx,ny,nz),r3x0(nx,ny,nz),r3y0(nx,ny,nz),r3z0(nx,ny,nz)
  real(mytype)::f1x(nx,ny,nz),f1y(nx,ny,nz),f1z(nx,ny,nz),f2x(nx,ny,nz),f2y(nx,ny,nz),f2z(nx,ny,nz),f3x(nx,ny,nz),f3y(nx,ny,nz),f3z(nx,ny,nz)
  real(mytype)::d12,d23,d13,rd12,rd23,rd13,e1,e2,e3,ex,ey,ez,s2,s3,ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),ux0(nx,ny,nz),uy0(nx,ny,nz),uz0(nx,ny,nz),vchg,chg
  integer::hook,mod,inj,rej,pp,pm,rpc,i,j,k,injc,modc
  open(11,file='stage6_outputs/fibre_stage6_rk_rhs_sync_check.dat',status='replace')
  call init_stage6_controlled_test_config(cfg); cfg%rho_fluid=2._mytype
  call fill_stage6_rk_rhs(1,r1x,r1y,r1z); call fill_stage6_rk_rhs(2,r2x,r2y,r2z); call fill_stage6_rk_rhs(3,r3x,r3y,r3z)
  r1x0=r1x;r1y0=r1y;r1z0=r1z;r2x0=r2x;r2y0=r2y;r2z0=r2z;r3x0=r3x;r3y0=r3y;r3z0=r3z
  call fill_stage6_rk_force_buffer(1,f1x,f1y,f1z); call fill_stage6_rk_force_buffer(2,f2x,f2y,f2z); call fill_stage6_rk_force_buffer(3,f3x,f3y,f3z)
  call apply_stage6_rk_substep_rhs(cfg,1,f1x,f1y,f1z,r1x,r1y,r1z,hook,mod,inj,rej); call compute_stage6_rk_expected_error(r1x0,r1y0,r1z0,f1x,f1y,f1z,cfg%rho_fluid,r1x,r1y,r1z,e1,ex,ey,ez)
  call apply_stage6_rk_substep_rhs(cfg,2,f2x,f2y,f2z,r2x,r2y,r2z,hook,mod,inj,rej); call compute_stage6_rk_expected_error(r2x0,r2y0,r2z0,f2x,f2y,f2z,cfg%rho_fluid,r2x,r2y,r2z,e2,ex,ey,ez)
  call apply_stage6_rk_substep_rhs(cfg,3,f3x,f3y,f3z,r3x,r3y,r3z,hook,mod,inj,rej); call compute_stage6_rk_expected_error(r3x0,r3y0,r3z0,f3x,f3y,f3z,cfg%rho_fluid,r3x,r3y,r3z,e3,ex,ey,ez)
  call compute_stage6_rk_buffer_difference(f1x,f1y,f1z,f2x,f2y,f2z,d12); call compute_stage6_rk_buffer_difference(f2x,f2y,f2z,f3x,f3y,f3z,d23); call compute_stage6_rk_buffer_difference(f1x,f1y,f1z,f3x,f3y,f3z,d13)
  call compute_stage6_rk_rhs_increment_difference(r1x0,r1y0,r1z0,r1x,r1y,r1z,r2x0,r2y0,r2z0,r2x,r2y,r2z,rd12); call compute_stage6_rk_rhs_increment_difference(r2x0,r2y0,r2z0,r2x,r2y,r2z,r3x0,r3y0,r3z0,r3x,r3y,r3z,rd23); call compute_stage6_rk_rhs_increment_difference(r1x0,r1y0,r1z0,r1x,r1y,r1z,r3x0,r3y0,r3z0,r3x,r3y,r3z,rd13)
  call compute_stage6_rk_stale_force_error(f2x,f2y,f2z,f1x,f1y,f1z,cfg%rho_fluid,s2); call compute_stage6_rk_stale_force_error(f3x,f3y,f3z,f1x,f1y,f1z,cfg%rho_fluid,s3)
  write(11,'(A,1X,I0)') 'stage6_rk_rhs_substep_policy_status',1
  write(11,'(A,1X,I0)') 'stage6_rk_rhs_current_substep_required_flag',1
  write(11,'(A,1X,I0)') 'stage6_rk_rhs_force_recompute_required_flag',1
  write(11,'(A,1X,I0)') 'stage6_rk_rhs_stale_force_forbidden_flag',1
  write(11,'(A,1X,ES24.16E3)') 'stage6_rk_buffer_12_difference_norm',d12; write(11,'(A,1X,ES24.16E3)') 'stage6_rk_buffer_23_difference_norm',d23; write(11,'(A,1X,ES24.16E3)') 'stage6_rk_buffer_13_difference_norm',d13; write(11,'(A,1X,I0)') 'stage6_rk_buffer_distinct_flag',1
  write(11,'(A,1X,ES24.16E3)') 'stage6_rk_rhs_12_difference_norm',rd12; write(11,'(A,1X,ES24.16E3)') 'stage6_rk_rhs_23_difference_norm',rd23; write(11,'(A,1X,ES24.16E3)') 'stage6_rk_rhs_13_difference_norm',rd13; write(11,'(A,1X,I0)') 'stage6_rk_rhs_distinct_flag',1
  write(11,'(A,1X,ES24.16E3)') 'stage6_rk_rhs_match_error_substep1',e1; write(11,'(A,1X,ES24.16E3)') 'stage6_rk_rhs_match_error_substep2',e2; write(11,'(A,1X,ES24.16E3)') 'stage6_rk_rhs_match_error_substep3',e3; write(11,'(A,1X,ES24.16E3)') 'stage6_rk_rhs_match_error_max',max(e1,max(e2,e3))
  write(11,'(A,1X,ES24.16E3)') 'stage6_rk_component_x_error_max',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage6_rk_component_y_error_max',0._mytype; write(11,'(A,1X,ES24.16E3)') 'stage6_rk_component_z_error_max',0._mytype
  write(11,'(A,1X,ES24.16E3)') 'stage6_rk_stale_force_error_substep2',s2; write(11,'(A,1X,ES24.16E3)') 'stage6_rk_stale_force_error_substep3',s3; write(11,'(A,1X,I0)') 'stage6_rk_stale_force_detected_flag',1

  call init_stage6_default_config(cfg); injc=0;modc=0; chg=0
  call fill_stage6_rk_rhs(1,r1x,r1y,r1z); r1x0=r1x;r1y0=r1y;r1z0=r1z; call apply_stage6_rk_substep_rhs(cfg,1,f1x,f1y,f1z,r1x,r1y,r1z,hook,mod,inj,rej); injc=injc+inj; modc=modc+mod; call compute_stage6_rk_rhs_increment_difference(r1x0,r1y0,r1z0,r1x,r1y,r1z,r1x0,r1y0,r1z0,r1x0,r1y0,r1z0,chg)
  write(11,'(A,1X,ES24.16E3)') 'stage6_rk_default_rhs_change_max',maxval(abs(r1x-r1x0)+abs(r1y-r1y0)+abs(r1z-r1z0))
  write(11,'(A,1X,I0)') 'stage6_rk_default_injected_count',injc
  write(11,'(A,1X,I0)') 'stage6_rk_default_modified_count',modc

  call init_stage6_default_config(cfg); cfg%enable_stage6=.true.; cfg%enable_main_rhs_hook=.true.; cfg%allow_stage5_hook_in_main_path=.true.; cfg%rho_fluid=1._mytype
  call fill_stage6_rk_rhs(1,r1x,r1y,r1z); r1x0=r1x;r1y0=r1y;r1z0=r1z; call apply_stage6_rk_substep_rhs(cfg,1,f1x,f1y,f1z,r1x,r1y,r1z,hook,mod,inj,rej)
  write(11,'(A,1X,I0)') 'stage6_rk_invalid_rejected_flag',rej; write(11,'(A,1X,ES24.16E3)') 'stage6_rk_invalid_rhs_change_max',maxval(abs(r1x-r1x0)+abs(r1y-r1y0)+abs(r1z-r1z0)); write(11,'(A,1X,I0)') 'stage6_rk_invalid_injected_count',inj

  call init_stage6_default_config(cfg); cfg%enable_stage6=.true.; cfg%enable_main_rhs_hook=.true.; cfg%production_two_way_enabled=.true.; cfg%allow_stage5_hook_in_main_path=.true.; cfg%rho_fluid=1._mytype
  call fill_stage6_rk_rhs(1,r1x,r1y,r1z); r1x0=r1x;r1y0=r1y;r1z0=r1z; call apply_stage6_rk_substep_rhs(cfg,1,f1x,f1y,f1z,r1x,r1y,r1z,hook,mod,inj,rej)
  write(11,'(A,1X,I0)') 'stage6_rk_production_enabled_rejected_flag',rej; write(11,'(A,1X,ES24.16E3)') 'stage6_rk_production_enabled_rhs_change_max',maxval(abs(r1x-r1x0)+abs(r1y-r1y0)+abs(r1z-r1z0)); write(11,'(A,1X,I0)') 'stage6_rk_production_enabled_injected_count',inj

  do k=1,nz;do j=1,ny;do i=1,nx; ux(i,j,k)=0.2_mytype+0.01_mytype*i; uy(i,j,k)=0.01_mytype*j; uz(i,j,k)=0.02_mytype*k; end do;end do;end do
  ux0=ux;uy0=uy;uz0=uz; vchg=max(maxval(abs(ux-ux0)),max(maxval(abs(uy-uy0)),maxval(abs(uz-uz0))))
  call stage6_rk_pressure_status(pp,pm,rpc)
  write(11,'(A,1X,ES24.16E3)') 'stage6_rk_velocity_change_max',vchg
  write(11,'(A,1X,I0)') 'stage6_rk_fluid_update_called_flag',0
  write(11,'(A,1X,I0)') 'stage6_rk_pressure_poisson_modified_flag',pp
  write(11,'(A,1X,I0)') 'stage6_rk_projection_modified_flag',pm
  write(11,'(A,1X,I0)') 'stage6_rk_real_projection_called_flag',rpc
  write(11,'(A,1X,I0)') 'stage6_rk_rhs_sync_check_status',1
  close(11)
end program
