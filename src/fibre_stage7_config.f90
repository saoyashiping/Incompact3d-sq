module fibre_stage7_config
  use fibre_parameters, only : mytype
  implicit none
  type stage7_config_t
    logical :: enable_stage7
    logical :: enable_nonuniform_y_ibm
    logical :: enable_component_specific_layout
    logical :: enable_real_layout_interpolation
    logical :: enable_real_layout_spreading
    logical :: enable_real_layout_rhs_candidate
    logical :: enable_controlled_real_layout_test
    logical :: production_dns_enabled
    logical :: production_two_way_enabled
    logical :: reject_invalid_config
    real(mytype) :: rho_fluid
    integer :: config_status
  end type
contains
  subroutine init_stage7_default_config(config)
    type(stage7_config_t), intent(out) :: config
    config%enable_stage7=.false.; config%enable_nonuniform_y_ibm=.false.; config%enable_component_specific_layout=.false.
    config%enable_real_layout_interpolation=.false.; config%enable_real_layout_spreading=.false.; config%enable_real_layout_rhs_candidate=.false.
    config%enable_controlled_real_layout_test=.false.; config%production_dns_enabled=.false.; config%production_two_way_enabled=.false.
    config%reject_invalid_config=.true.; config%rho_fluid=1._mytype; config%config_status=1
  end subroutine
  subroutine init_stage7_controlled_real_layout_config(config)
    type(stage7_config_t), intent(out) :: config
    config%enable_stage7=.true.; config%enable_nonuniform_y_ibm=.true.; config%enable_component_specific_layout=.true.
    config%enable_real_layout_interpolation=.true.; config%enable_real_layout_spreading=.true.; config%enable_real_layout_rhs_candidate=.true.
    config%enable_controlled_real_layout_test=.true.; config%production_dns_enabled=.false.; config%production_two_way_enabled=.false.
    config%reject_invalid_config=.true.; config%rho_fluid=1._mytype; config%config_status=1
  end subroutine
  subroutine validate_stage7_config(config, valid, rejected_flag, rhs_allowed_flag, controlled_test_flag, production_dns_flag, production_twoway_flag)
    type(stage7_config_t), intent(in) :: config
    integer, intent(out) :: valid,rejected_flag,rhs_allowed_flag,controlled_test_flag,production_dns_flag,production_twoway_flag
    valid=1; rejected_flag=0; rhs_allowed_flag=0; controlled_test_flag=0; production_dns_flag=0; production_twoway_flag=0
    if (config%rho_fluid<=0._mytype) then
      valid=0
    else if (config%production_dns_enabled) then
      valid=0; production_dns_flag=1
    else if (config%production_two_way_enabled) then
      valid=0; production_twoway_flag=1
    else if ((.not.config%enable_stage7) .and. (config%enable_nonuniform_y_ibm .or. config%enable_component_specific_layout .or. &
             config%enable_real_layout_interpolation .or. config%enable_real_layout_spreading .or. config%enable_real_layout_rhs_candidate .or. &
             config%enable_controlled_real_layout_test .or. config%production_dns_enabled .or. config%production_two_way_enabled)) then
      valid=0
    else if ((.not.config%enable_controlled_real_layout_test) .and. (config%enable_nonuniform_y_ibm .or. config%enable_component_specific_layout .or. &
             config%enable_real_layout_interpolation .or. config%enable_real_layout_spreading .or. config%enable_real_layout_rhs_candidate)) then
      valid=0
    else if (config%enable_real_layout_rhs_candidate .and. ((.not.config%enable_real_layout_interpolation) .or. (.not.config%enable_real_layout_spreading))) then
      valid=0
    else if (config%enable_real_layout_interpolation .neqv. config%enable_real_layout_spreading) then
      valid=0
    else if ((.not.config%enable_stage7) .and. (.not.config%enable_nonuniform_y_ibm) .and. (.not.config%enable_component_specific_layout) .and. &
             (.not.config%enable_real_layout_interpolation) .and. (.not.config%enable_real_layout_spreading) .and. (.not.config%enable_real_layout_rhs_candidate) .and. &
             (.not.config%enable_controlled_real_layout_test) .and. (.not.config%production_dns_enabled) .and. (.not.config%production_two_way_enabled)) then
      valid=1; rhs_allowed_flag=0; controlled_test_flag=0
    else if (config%enable_stage7 .and. config%enable_nonuniform_y_ibm .and. config%enable_component_specific_layout .and. config%enable_real_layout_interpolation .and. &
             config%enable_real_layout_spreading .and. config%enable_real_layout_rhs_candidate .and. config%enable_controlled_real_layout_test .and. &
             (.not.config%production_dns_enabled) .and. (.not.config%production_two_way_enabled)) then
      valid=1; rhs_allowed_flag=1; controlled_test_flag=1
    else
      valid=0
    end if
    if (valid==0) then
      rejected_flag=1; rhs_allowed_flag=0; controlled_test_flag=0
    end if
  end subroutine
  subroutine stage7_config_noop_rhs_guard(config, rhsx, rhsy, rhsz, rhs_change_max, rhs_modified_flag)
    type(stage7_config_t), intent(in) :: config
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
  subroutine stage7_config_pressure_status(pressure_poisson_modified_flag, projection_modified_flag, real_projection_called_flag, production_dns_called_flag, fluid_update_called_flag, fibre_advance_called_flag)
    integer, intent(out) :: pressure_poisson_modified_flag, projection_modified_flag, real_projection_called_flag, production_dns_called_flag, fluid_update_called_flag, fibre_advance_called_flag
    pressure_poisson_modified_flag=0; projection_modified_flag=0; real_projection_called_flag=0
    production_dns_called_flag=0; fluid_update_called_flag=0; fibre_advance_called_flag=0
  end subroutine
end module fibre_stage7_config
