module fibre_stage5_rhs_wrapper
  use fibre_parameters, only : mytype
  use fibre_stage5_config, only : stage5_config_t, validate_stage5_config
  use fibre_stage5_rhs_audit, only : stage5_rhs_audit_t
  implicit none

contains

  subroutine apply_stage5_ibm_force_to_fluid_rhs(config, fx, fy, fz, rhsx, rhsy, rhsz, rhs_modified_flag, injected_flag, rejected_flag)
    type(stage5_config_t), intent(in) :: config
    real(mytype), intent(in) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
    real(mytype), intent(inout) :: rhsx(:,:,:), rhsy(:,:,:), rhsz(:,:,:)
    integer, intent(out) :: rhs_modified_flag, injected_flag, rejected_flag
    integer :: valid, two_way_enabled_flag, rhs_allowed_flag
    real(mytype) :: buffer_max_abs

    if (any(shape(rhsx) /= shape(fx)) .or. any(shape(rhsy) /= shape(fy)) .or. any(shape(rhsz) /= shape(fz))) then
      error stop 'stage5 rhs wrapper shape mismatch'
    end if

    call validate_stage5_config(config, valid, rejected_flag, two_way_enabled_flag, rhs_allowed_flag)
    rhs_modified_flag = 0
    injected_flag = 0

    if (rejected_flag==1 .or. valid==0) return
    if (rhs_allowed_flag==0) return
    if (config%rho_fluid<=0._mytype) then
      rejected_flag = 1
      return
    end if

    rhsx = rhsx + fx / config%rho_fluid
    rhsy = rhsy + fy / config%rho_fluid
    rhsz = rhsz + fz / config%rho_fluid
    injected_flag = 1

    buffer_max_abs = max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))
    if (buffer_max_abs > 0._mytype) rhs_modified_flag = 1
  end subroutine

  subroutine compute_stage5_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx1,rhsy1,rhsz1,change_max)
    real(mytype), intent(in) :: rhsx0(:,:,:), rhsy0(:,:,:), rhsz0(:,:,:)
    real(mytype), intent(in) :: rhsx1(:,:,:), rhsy1(:,:,:), rhsz1(:,:,:)
    real(mytype), intent(out) :: change_max
    change_max = max(maxval(abs(rhsx1-rhsx0)),max(maxval(abs(rhsy1-rhsy0)),maxval(abs(rhsz1-rhsz0))))
  end subroutine

  subroutine compute_stage5_rhs_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,rho,rhsx1,rhsy1,rhsz1,error_max)
    real(mytype), intent(in) :: rhsx0(:,:,:), rhsy0(:,:,:), rhsz0(:,:,:)
    real(mytype), intent(in) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
    real(mytype), intent(in) :: rho
    real(mytype), intent(in) :: rhsx1(:,:,:), rhsy1(:,:,:), rhsz1(:,:,:)
    real(mytype), intent(out) :: error_max
    real(mytype), allocatable :: ex(:,:,:), ey(:,:,:), ez(:,:,:)

    allocate(ex(size(rhsx0,1),size(rhsx0,2),size(rhsx0,3)))
    allocate(ey(size(rhsy0,1),size(rhsy0,2),size(rhsy0,3)))
    allocate(ez(size(rhsz0,1),size(rhsz0,2),size(rhsz0,3)))
    ex = rhsx0 + fx/rho
    ey = rhsy0 + fy/rho
    ez = rhsz0 + fz/rho
    error_max = max(maxval(abs(rhsx1-ex)),max(maxval(abs(rhsy1-ey)),maxval(abs(rhsz1-ez))))
    deallocate(ex,ey,ez)
  end subroutine

end module fibre_stage5_rhs_wrapper
