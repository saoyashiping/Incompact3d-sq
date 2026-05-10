module fibre_stage5_main_hook
  use fibre_parameters, only : mytype
  use fibre_stage5_config, only : stage5_config_t
  use fibre_stage5_rhs_wrapper, only : apply_stage5_ibm_force_to_fluid_rhs
  implicit none
contains
  subroutine apply_stage5_main_rhs_hook(config, fx, fy, fz, rhsx, rhsy, rhsz, hook_called_flag, rhs_modified_flag, injected_flag, rejected_flag)
    type(stage5_config_t), intent(in) :: config
    real(mytype), intent(in) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
    real(mytype), intent(inout) :: rhsx(:,:,:), rhsy(:,:,:), rhsz(:,:,:)
    integer, intent(out) :: hook_called_flag, rhs_modified_flag, injected_flag, rejected_flag
    hook_called_flag = 1
    call apply_stage5_ibm_force_to_fluid_rhs(config, fx, fy, fz, rhsx, rhsy, rhsz, rhs_modified_flag, injected_flag, rejected_flag)
  end subroutine

  subroutine stage5_default_main_dns_autocall_enabled(flag)
    integer, intent(out) :: flag
    flag = 0
  end subroutine

  subroutine stage5_pressure_poisson_hook_status(modified_flag, projection_after_rhs_policy_flag)
    integer, intent(out) :: modified_flag, projection_after_rhs_policy_flag
    modified_flag = 0
    projection_after_rhs_policy_flag = 1
  end subroutine
end module fibre_stage5_main_hook
