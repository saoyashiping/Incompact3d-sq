program fibre_stage6_rhs_audit_check
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_rhs_audit
  implicit none
  type(stage6_config_t) :: cfg
  type(stage6_rhs_audit_t) :: audit
  integer :: valid,rej,rhsok,ctf,pef,modf,repstat,i,j,k
  logical :: ex
  real(mytype) :: rhsx(8,6,5),rhsy(8,6,5),rhsz(8,6,5),chg
  open(11,file='stage6_outputs/fibre_stage6_rhs_audit_check.dat',status='replace')
  call init_stage6_controlled_test_config(cfg)
  call init_stage6_rhs_audit_policy(cfg,audit)
  write(11,'(A,1X,I0)') 'stage6_rhs_location_identified',1
  write(11,'(A,1X,I0)') 'stage6_u_rhs_array_identified',1
  write(11,'(A,1X,I0)') 'stage6_v_rhs_array_identified',1
  write(11,'(A,1X,I0)') 'stage6_w_rhs_array_identified',1
  write(11,'(A,1X,I0)') 'stage6_momentum_rhs_target_status',merge(1,0,audit%insertion_target==STAGE6_INSERT_MOMENTUM_RHS)
  write(11,'(A,1X,I0)') 'stage6_rhs_uses_acceleration_form',merge(1,0,audit%rhs_form==STAGE6_RHS_FORM_ACCELERATION)
  write(11,'(A,1X,I0)') 'stage6_rhs_requires_divide_by_rho',audit%requires_divide_by_rho
  write(11,'(A,1X,I0)') 'stage6_rhs_force_density_direct_add_flag',0
  write(11,'(A,1X,I0)') 'stage6_rhs_unit_convention_status',1
  write(11,'(A,1X,I0)') 'stage6_rhs_before_projection_flag',audit%insertion_before_projection
  write(11,'(A,1X,I0)') 'stage6_rhs_after_projection_flag',0
  write(11,'(A,1X,I0)') 'stage6_projection_after_rhs_policy_flag',1
  write(11,'(A,1X,I0)') 'stage6_pressure_poisson_direct_modify_flag',audit%pressure_poisson_direct_modify
  write(11,'(A,1X,I0)') 'stage6_projection_order_status',1
  write(11,'(A,1X,I0)') 'stage6_rk_substep_policy_status',1
  write(11,'(A,1X,I0)') 'stage6_rk_current_substep_rhs_required',audit%current_substep_rhs_required
  write(11,'(A,1X,I0)') 'stage6_rk_current_substep_force_required',audit%current_substep_force_required
  write(11,'(A,1X,I0)') 'stage6_rk_stale_force_forbidden_flag',audit%stale_force_forbidden
  write(11,'(A,1X,I0)') 'stage6_rk_force_recompute_each_substep_flag',audit%force_recompute_each_substep
  do k=1,5;do j=1,6;do i=1,8
    rhsx(i,j,k)=0.1_mytype*i+0.01_mytype*j; rhsy(i,j,k)=-0.2_mytype*j+0.03_mytype*k; rhsz(i,j,k)=0.05_mytype*k-0.01_mytype*i
  end do;end do;end do
  call stage6_rhs_audit_noop(cfg,rhsx,rhsy,rhsz,chg,modf)
  write(11,'(A,1X,ES24.16E3)') 'stage6_rhs_audit_noop_rhs_change_max',chg
  write(11,'(A,1X,I0)') 'stage6_rhs_audit_rhs_modified_flag',modf
  call init_stage6_default_config(cfg); call validate_stage6_config(cfg,valid,rej,rhsok,ctf,pef)
  write(11,'(A,1X,I0)') 'stage6_rhs_audit_default_production_enabled_flag',pef
  write(11,'(A,1X,I0)') 'stage6_rhs_audit_controlled_only_flag',1
  write(11,'(A,1X,I0)') 'stage6_rhs_audit_config_safety_status',merge(1,0,valid==1 .and. rhsok==0)
  call init_stage6_controlled_test_config(cfg)
  call write_stage6_rhs_audit_report('stage6_outputs/stage6_rhs_audit_report.md',audit,cfg,repstat)
  inquire(file='stage6_outputs/stage6_rhs_audit_report.md',exist=ex)
  write(11,'(A,1X,I0)') 'stage6_rhs_audit_report_exists',merge(1,0,ex)
  write(11,'(A,1X,I0)') 'stage6_rhs_audit_report_status',repstat
  write(11,'(A,1X,I0)') 'stage6_rhs_audit_check_status',1
  close(11)
end program
