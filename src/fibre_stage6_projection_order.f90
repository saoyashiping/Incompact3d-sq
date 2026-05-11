module fibre_stage6_projection_order
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_rhs_audit
  use fibre_stage6_controlled_rhs_hook
  implicit none
contains
  subroutine apply_stage6_projection_order_rhs(config,fx,fy,fz,rhsx,rhsy,rhsz,hook_called_flag,rhs_modified_flag,injected_flag,rejected_flag,rhs_before_projection_flag,rhs_after_projection_flag,pressure_poisson_direct_modify_flag)
    type(stage6_config_t), intent(in) :: config
    real(mytype), intent(in) :: fx(:,:,:),fy(:,:,:),fz(:,:,:)
    real(mytype), intent(inout) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    integer, intent(out) :: hook_called_flag,rhs_modified_flag,injected_flag,rejected_flag,rhs_before_projection_flag,rhs_after_projection_flag,pressure_poisson_direct_modify_flag
    rhs_before_projection_flag=1; rhs_after_projection_flag=0; pressure_poisson_direct_modify_flag=0
    call apply_stage6_controlled_rhs_hook(config,fx,fy,fz,rhsx,rhsy,rhsz,hook_called_flag,rhs_modified_flag,injected_flag,rejected_flag)
  end subroutine
  subroutine compute_stage6_intermediate_velocity(dt,ux0,uy0,uz0,rhsx,rhsy,rhsz,uxstar,uystar,uzstar)
    real(mytype), intent(in) :: dt,ux0(:,:,:),uy0(:,:,:),uz0(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    real(mytype), intent(out) :: uxstar(:,:,:),uystar(:,:,:),uzstar(:,:,:)
    uxstar=ux0+dt*rhsx; uystar=uy0+dt*rhsy; uzstar=uz0+dt*rhsz
  end subroutine
  subroutine compute_stage6_ustar_expected_error(dt,ux0,uy0,uz0,rhsx,rhsy,rhsz,uxstar,uystar,uzstar,err)
    real(mytype), intent(in) :: dt,ux0(:,:,:),uy0(:,:,:),uz0(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),uxstar(:,:,:),uystar(:,:,:),uzstar(:,:,:)
    real(mytype), intent(out) :: err
    err=max(maxval(abs(uxstar-(ux0+dt*rhsx))),max(maxval(abs(uystar-(uy0+dt*rhsy))),maxval(abs(uzstar-(uz0+dt*rhsz)))))
  end subroutine
  subroutine stage6_projection_policy_status(projection_after_rhs_policy_flag,post_projection_velocity_modified_flag,post_projection_direct_forcing_forbidden_flag,projection_modified_flag,real_projection_called_flag)
    integer, intent(out) :: projection_after_rhs_policy_flag,post_projection_velocity_modified_flag,post_projection_direct_forcing_forbidden_flag,projection_modified_flag,real_projection_called_flag
    projection_after_rhs_policy_flag=1; post_projection_velocity_modified_flag=0; post_projection_direct_forcing_forbidden_flag=1; projection_modified_flag=0; real_projection_called_flag=0
  end subroutine
end module
