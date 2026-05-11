module fibre_stage6_controlled_rhs_hook
  use fibre_parameters, only : mytype
  use fibre_stage6_config, only : stage6_config_t, validate_stage6_config
  use fibre_stage6_rhs_audit
  implicit none
contains
  subroutine apply_stage6_controlled_rhs_hook(config,fx,fy,fz,rhsx,rhsy,rhsz,hook_called_flag,rhs_modified_flag,injected_flag,rejected_flag)
    type(stage6_config_t), intent(in) :: config
    real(mytype), intent(in) :: fx(:,:,:),fy(:,:,:),fz(:,:,:)
    real(mytype), intent(inout) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    integer, intent(out) :: hook_called_flag,rhs_modified_flag,injected_flag,rejected_flag
    integer :: valid,rhs_allowed_flag,controlled_test_flag,production_enabled_flag
    hook_called_flag=1; rhs_modified_flag=0; injected_flag=0
    if (any(shape(rhsx)/=shape(fx)) .or. any(shape(rhsy)/=shape(fy)) .or. any(shape(rhsz)/=shape(fz))) then
      rejected_flag=1; return
    end if
    call validate_stage6_config(config,valid,rejected_flag,rhs_allowed_flag,controlled_test_flag,production_enabled_flag)
    if (rejected_flag==1 .or. valid==0) return
    if (rhs_allowed_flag==0 .or. controlled_test_flag==0) return
    if (config%rho_fluid<=0._mytype) then; rejected_flag=1; return; end if
    rhsx=rhsx+fx/config%rho_fluid; rhsy=rhsy+fy/config%rho_fluid; rhsz=rhsz+fz/config%rho_fluid
    injected_flag=1
    if (max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))>0._mytype) rhs_modified_flag=1
  end subroutine
  subroutine compute_stage6_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx1,rhsy1,rhsz1,change_max)
    real(mytype), intent(in) :: rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:),rhsx1(:,:,:),rhsy1(:,:,:),rhsz1(:,:,:)
    real(mytype), intent(out) :: change_max
    change_max=max(maxval(abs(rhsx1-rhsx0)),max(maxval(abs(rhsy1-rhsy0)),maxval(abs(rhsz1-rhsz0))))
  end subroutine
  subroutine compute_stage6_rhs_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,rho,rhsx1,rhsy1,rhsz1,err,ex,ey,ez)
    real(mytype), intent(in) :: rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:),fx(:,:,:),fy(:,:,:),fz(:,:,:),rho,rhsx1(:,:,:),rhsy1(:,:,:),rhsz1(:,:,:)
    real(mytype), intent(out) :: err,ex,ey,ez
    ex=maxval(abs((rhsx1-rhsx0)-fx/rho)); ey=maxval(abs((rhsy1-rhsy0)-fy/rho)); ez=maxval(abs((rhsz1-rhsz0)-fz/rho)); err=max(ex,max(ey,ez))
  end subroutine
  subroutine stage6_controlled_pressure_status(pressure_poisson_modified_flag, projection_modified_flag, pressure_rhs_modified_flag)
    integer, intent(out) :: pressure_poisson_modified_flag, projection_modified_flag, pressure_rhs_modified_flag
    pressure_poisson_modified_flag=0; projection_modified_flag=0; pressure_rhs_modified_flag=0
  end subroutine
end module
