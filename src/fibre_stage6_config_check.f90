program fibre_stage6_config_check
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  implicit none
  type(stage6_config_t) :: cfg
  integer :: valid,rej,rhsok,ctf,pef,modf,i,j,k
  logical :: ex
  real(mytype) :: chg
  real(mytype) :: rhsx(8,6,5),rhsy(8,6,5),rhsz(8,6,5)
  open(11,file='stage6_outputs/fibre_stage6_config_check.dat',status='replace')
  inquire(file='stage5_checks/STAGE5_CLOSED.md',exist=ex)
  write(11,'(A,1X,I0)') 'stage6_stage5_closed_marker_status',merge(1,0,ex)
  write(11,'(A,1X,I0)') 'stage6_stage5_dependency_status',merge(1,0,ex)

  call init_stage6_default_config(cfg); call validate_stage6_config(cfg,valid,rej,rhsok,ctf,pef)
  write(11,'(A,1X,I0)') 'stage6_default_enable_stage6',merge(1,0,cfg%enable_stage6)
  write(11,'(A,1X,I0)') 'stage6_default_main_rhs_hook_enabled',merge(1,0,cfg%enable_main_rhs_hook)
  write(11,'(A,1X,I0)') 'stage6_default_controlled_test_enabled',merge(1,0,cfg%enable_controlled_rhs_test)
  write(11,'(A,1X,I0)') 'stage6_default_production_two_way_enabled',merge(1,0,cfg%production_two_way_enabled)
  write(11,'(A,1X,I0)') 'stage6_default_allow_stage5_hook_in_main_path',merge(1,0,cfg%allow_stage5_hook_in_main_path)
  write(11,'(A,1X,I0)') 'stage6_default_valid_flag',valid; write(11,'(A,1X,I0)') 'stage6_default_rejected_flag',rej; write(11,'(A,1X,I0)') 'stage6_default_rhs_allowed_flag',rhsok; write(11,'(A,1X,I0)') 'stage6_default_config_status',cfg%config_status

  call init_stage6_controlled_test_config(cfg); call validate_stage6_config(cfg,valid,rej,rhsok,ctf,pef)
  write(11,'(A,1X,I0)') 'stage6_controlled_enable_stage6',merge(1,0,cfg%enable_stage6)
  write(11,'(A,1X,I0)') 'stage6_controlled_main_rhs_hook_enabled',merge(1,0,cfg%enable_main_rhs_hook)
  write(11,'(A,1X,I0)') 'stage6_controlled_test_enabled',merge(1,0,cfg%enable_controlled_rhs_test)
  write(11,'(A,1X,I0)') 'stage6_controlled_production_two_way_enabled',merge(1,0,cfg%production_two_way_enabled)
  write(11,'(A,1X,I0)') 'stage6_controlled_allow_stage5_hook_in_main_path',merge(1,0,cfg%allow_stage5_hook_in_main_path)
  write(11,'(A,1X,I0)') 'stage6_controlled_valid_flag',valid; write(11,'(A,1X,I0)') 'stage6_controlled_rejected_flag',rej; write(11,'(A,1X,I0)') 'stage6_controlled_rhs_allowed_flag',rhsok; write(11,'(A,1X,I0)') 'stage6_controlled_config_status',cfg%config_status

  call init_stage6_default_config(cfg); cfg%enable_stage6=.true.; cfg%enable_main_rhs_hook=.true.; cfg%production_two_way_enabled=.true.; cfg%allow_stage5_hook_in_main_path=.true.
  call validate_stage6_config(cfg,valid,rej,rhsok,ctf,pef)
  write(11,'(A,1X,I0)') 'stage6_invalid_production_enabled_rejected_flag',rej; write(11,'(A,1X,I0)') 'stage6_invalid_production_enabled_valid_flag',valid; write(11,'(A,1X,I0)') 'stage6_invalid_production_enabled_rhs_allowed_flag',rhsok

  call init_stage6_default_config(cfg); cfg%enable_stage6=.true.; cfg%enable_main_rhs_hook=.true.; cfg%allow_stage5_hook_in_main_path=.true.
  call validate_stage6_config(cfg,valid,rej,rhsok,ctf,pef)
  write(11,'(A,1X,I0)') 'stage6_invalid_hook_without_controlled_test_rejected_flag',rej; write(11,'(A,1X,I0)') 'stage6_invalid_hook_without_controlled_test_valid_flag',valid; write(11,'(A,1X,I0)') 'stage6_invalid_hook_without_controlled_test_rhs_allowed_flag',rhsok

  call init_stage6_default_config(cfg); cfg%enable_stage6=.true.; cfg%enable_main_rhs_hook=.true.; cfg%enable_controlled_rhs_test=.true.
  call validate_stage6_config(cfg,valid,rej,rhsok,ctf,pef)
  write(11,'(A,1X,I0)') 'stage6_invalid_controlled_without_allow_hook_rejected_flag',rej; write(11,'(A,1X,I0)') 'stage6_invalid_controlled_without_allow_hook_valid_flag',valid

  call init_stage6_controlled_test_config(cfg); cfg%rho_fluid=0._mytype; call validate_stage6_config(cfg,valid,rej,rhsok,ctf,pef)
  write(11,'(A,1X,I0)') 'stage6_invalid_rho_rejected_flag',rej; write(11,'(A,1X,I0)') 'stage6_invalid_rho_valid_flag',valid

  do k=1,5;do j=1,6;do i=1,8
    rhsx(i,j,k)=0.1_mytype*i+0.01_mytype*j; rhsy(i,j,k)=-0.2_mytype*j+0.03_mytype*k; rhsz(i,j,k)=0.05_mytype*k-0.01_mytype*i
  end do;end do;end do
  call stage6_config_noop_rhs_guard(cfg,rhsx,rhsy,rhsz,chg,modf)
  write(11,'(A,1X,ES24.16E3)') 'stage6_config_noop_rhs_change_max',chg
  write(11,'(A,1X,I0)') 'stage6_config_noop_rhs_modified_flag',modf
  write(11,'(A,1X,I0)') 'stage6_config_check_status',1
  close(11)
end program
