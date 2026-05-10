module fibre_stage6_main_noop_hook
  use fibre_parameters, only : mytype
  use fibre_stage6_config, only : stage6_config_t, validate_stage6_config
  use fibre_stage6_rhs_audit
  implicit none
contains
  subroutine apply_stage6_main_noop_hook(config,fx,fy,fz,rhsx,rhsy,rhsz,hook_called_flag,rhs_modified_flag,injected_flag,rejected_flag)
    type(stage6_config_t), intent(in) :: config
    real(mytype), intent(in) :: fx(:,:,:),fy(:,:,:),fz(:,:,:)
    real(mytype), intent(inout) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    integer, intent(out) :: hook_called_flag,rhs_modified_flag,injected_flag,rejected_flag
    integer :: valid,rhs_allowed_flag,controlled_test_flag,production_enabled_flag
    hook_called_flag=1; rhs_modified_flag=0; injected_flag=0
    call validate_stage6_config(config,valid,rejected_flag,rhs_allowed_flag,controlled_test_flag,production_enabled_flag)
    if (rejected_flag==1 .or. valid==0) return
  end subroutine
  subroutine stage6_noop_hook_pressure_status(pressure_poisson_modified_flag, projection_modified_flag, pressure_rhs_modified_flag)
    integer, intent(out) :: pressure_poisson_modified_flag, projection_modified_flag, pressure_rhs_modified_flag
    pressure_poisson_modified_flag=0; projection_modified_flag=0; pressure_rhs_modified_flag=0
  end subroutine
  subroutine stage6_noop_hook_production_safety_status(default_main_dns_safe_flag, production_enabled_by_default_flag, controlled_test_only_flag)
    integer, intent(out) :: default_main_dns_safe_flag, production_enabled_by_default_flag, controlled_test_only_flag
    default_main_dns_safe_flag=1; production_enabled_by_default_flag=0; controlled_test_only_flag=1
  end subroutine
end module
