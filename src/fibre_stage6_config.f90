module fibre_stage6_config
  use fibre_parameters, only : mytype
  implicit none
  type stage6_config_t
    logical :: enable_stage6
    logical :: enable_main_rhs_hook
    logical :: enable_controlled_rhs_test
    logical :: production_two_way_enabled
    logical :: allow_stage5_hook_in_main_path
    logical :: reject_invalid_config
    real(mytype) :: rho_fluid
    integer :: config_status
  end type stage6_config_t
contains
  subroutine init_stage6_default_config(config)
    type(stage6_config_t), intent(out) :: config
    config%enable_stage6=.false.; config%enable_main_rhs_hook=.false.; config%enable_controlled_rhs_test=.false.
    config%production_two_way_enabled=.false.; config%allow_stage5_hook_in_main_path=.false.; config%reject_invalid_config=.true.
    config%rho_fluid=1._mytype; config%config_status=1
  end subroutine
  subroutine init_stage6_controlled_test_config(config)
    type(stage6_config_t), intent(out) :: config
    config%enable_stage6=.true.; config%enable_main_rhs_hook=.true.; config%enable_controlled_rhs_test=.true.
    config%production_two_way_enabled=.false.; config%allow_stage5_hook_in_main_path=.true.; config%reject_invalid_config=.true.
    config%rho_fluid=1._mytype; config%config_status=1
  end subroutine
  subroutine validate_stage6_config(config, valid, rejected_flag, rhs_allowed_flag, controlled_test_flag, production_enabled_flag)
    type(stage6_config_t), intent(in) :: config
    integer, intent(out) :: valid,rejected_flag,rhs_allowed_flag,controlled_test_flag,production_enabled_flag
    valid=1;rejected_flag=0;rhs_allowed_flag=0;controlled_test_flag=0;production_enabled_flag=0
    if (config%rho_fluid<=0._mytype) then
      valid=0
    else if (config%production_two_way_enabled) then
      valid=0; production_enabled_flag=1
    else if ((.not.config%enable_stage6) .and. (config%enable_main_rhs_hook .or. config%enable_controlled_rhs_test .or. config%allow_stage5_hook_in_main_path)) then
      valid=0
    else if ((.not.config%enable_stage6) .and. (.not.config%enable_main_rhs_hook) .and. (.not.config%enable_controlled_rhs_test) .and. (.not.config%allow_stage5_hook_in_main_path)) then
      valid=1
    else if (config%enable_main_rhs_hook .and. (.not.config%enable_controlled_rhs_test)) then
      valid=0
    else if (config%enable_controlled_rhs_test .and. (.not.config%allow_stage5_hook_in_main_path)) then
      valid=0
    else if (config%enable_stage6 .and. config%enable_main_rhs_hook .and. config%enable_controlled_rhs_test .and. config%allow_stage5_hook_in_main_path) then
      valid=1; rhs_allowed_flag=1; controlled_test_flag=1
    else
      valid=0
    end if
    if (valid==0) then
      rejected_flag=1; rhs_allowed_flag=0; controlled_test_flag=0
    end if
  end subroutine
  subroutine stage6_config_noop_rhs_guard(config, rhsx, rhsy, rhsz, rhs_change_max, rhs_modified_flag)
    type(stage6_config_t), intent(in) :: config
    real(mytype), intent(inout) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    real(mytype), intent(out) :: rhs_change_max
    integer, intent(out) :: rhs_modified_flag
    real(mytype), allocatable :: x0(:,:,:),y0(:,:,:),z0(:,:,:)
    allocate(x0(size(rhsx,1),size(rhsx,2),size(rhsx,3)),y0(size(rhsy,1),size(rhsy,2),size(rhsy,3)),z0(size(rhsz,1),size(rhsz,2),size(rhsz,3)))
    x0=rhsx; y0=rhsy; z0=rhsz
    rhs_change_max=max(maxval(abs(rhsx-x0)),max(maxval(abs(rhsy-y0)),maxval(abs(rhsz-z0))))
    rhs_modified_flag=0
    deallocate(x0,y0,z0)
  end subroutine
end module fibre_stage6_config
