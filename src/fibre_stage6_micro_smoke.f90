module fibre_stage6_micro_smoke
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_rhs_audit
  use fibre_stage6_controlled_rhs_hook
  use fibre_stage6_projection_order, only : compute_stage6_intermediate_velocity
  use fibre_stage6_layout_guard
  implicit none
contains
  subroutine perform_stage6_micro_controlled_step(config, layout, dt, fx,fy,fz, rhsx,rhsy,rhsz, ux,uy,uz, &
                                                  uxstar,uystar,uzstar, &
                                                  config_valid_flag, layout_guard_pass_flag, hook_called_flag, &
                                                  rhs_injection_called_flag, injected_flag, modified_flag, rejected_flag, &
                                                  rhs_expected_error, ex,ey,ez, ustar_expected_error)
    type(stage6_config_t), intent(in) :: config
    type(stage6_layout_t), intent(inout) :: layout
    real(mytype), intent(in) :: dt, fx(:,:,:),fy(:,:,:),fz(:,:,:), ux(:,:,:),uy(:,:,:),uz(:,:,:)
    real(mytype), intent(inout) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    real(mytype), intent(out) :: uxstar(:,:,:),uystar(:,:,:),uzstar(:,:,:)
    integer, intent(out) :: config_valid_flag, layout_guard_pass_flag, hook_called_flag
    integer, intent(out) :: rhs_injection_called_flag, injected_flag, modified_flag, rejected_flag
    real(mytype), intent(out) :: rhs_expected_error, ex,ey,ez, ustar_expected_error
    integer :: valid, rhs_allowed_flag, controlled_test_flag, production_enabled_flag
    real(mytype), allocatable :: rhsx0(:,:,:), rhsy0(:,:,:), rhsz0(:,:,:)

    allocate(rhsx0(size(rhsx,1),size(rhsx,2),size(rhsx,3)))
    allocate(rhsy0(size(rhsy,1),size(rhsy,2),size(rhsy,3)))
    allocate(rhsz0(size(rhsz,1),size(rhsz,2),size(rhsz,3)))
    rhsx0 = rhsx; rhsy0 = rhsy; rhsz0 = rhsz

    call validate_stage6_config(config,valid,rejected_flag,rhs_allowed_flag,controlled_test_flag,production_enabled_flag)
    config_valid_flag = valid
    call validate_stage6_layout_guard(layout)
    layout_guard_pass_flag = layout%ordinary_path_allowed_flag

    hook_called_flag = 0; rhs_injection_called_flag = 0; injected_flag = 0; modified_flag = 0
    rhs_expected_error = 0._mytype; ex = 0._mytype; ey = 0._mytype; ez = 0._mytype

    if (valid == 0) then
      rejected_flag = 1
      call compute_stage6_intermediate_velocity(dt,ux,uy,uz,rhsx,rhsy,rhsz,uxstar,uystar,uzstar)
      ustar_expected_error = 0._mytype
      deallocate(rhsx0,rhsy0,rhsz0)
      return
    end if

    if (layout%ordinary_path_allowed_flag == 0) then
      rejected_flag = 1
      call compute_stage6_intermediate_velocity(dt,ux,uy,uz,rhsx,rhsy,rhsz,uxstar,uystar,uzstar)
      ustar_expected_error = 0._mytype
      deallocate(rhsx0,rhsy0,rhsz0)
      return
    end if

    call apply_stage6_controlled_rhs_hook(config,fx,fy,fz,rhsx,rhsy,rhsz,hook_called_flag,modified_flag,injected_flag,rejected_flag)
    rhs_injection_called_flag = injected_flag
    call compute_stage6_intermediate_velocity(dt,ux,uy,uz,rhsx,rhsy,rhsz,uxstar,uystar,uzstar)
    call compute_stage6_rhs_expected_error(rhsx0,rhsy0,rhsz0,fx,fy,fz,config%rho_fluid,rhsx,rhsy,rhsz,rhs_expected_error,ex,ey,ez)
    ustar_expected_error = max(maxval(abs(uxstar-(ux+dt*rhsx))),max(maxval(abs(uystar-(uy+dt*rhsy))),maxval(abs(uzstar-(uz+dt*rhsz)))))
    deallocate(rhsx0,rhsy0,rhsz0)
  end subroutine

  subroutine compute_stage6_micro_rhs_change_max(rhsx0,rhsy0,rhsz0,rhsx1,rhsy1,rhsz1,change_max)
    real(mytype), intent(in) :: rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:),rhsx1(:,:,:),rhsy1(:,:,:),rhsz1(:,:,:)
    real(mytype), intent(out) :: change_max
    change_max=max(maxval(abs(rhsx1-rhsx0)),max(maxval(abs(rhsy1-rhsy0)),maxval(abs(rhsz1-rhsz0))))
  end subroutine

  subroutine compute_stage6_micro_velocity_change_max(ux0,uy0,uz0,ux1,uy1,uz1,change_max)
    real(mytype), intent(in) :: ux0(:,:,:),uy0(:,:,:),uz0(:,:,:),ux1(:,:,:),uy1(:,:,:),uz1(:,:,:)
    real(mytype), intent(out) :: change_max
    change_max=max(maxval(abs(ux1-ux0)),max(maxval(abs(uy1-uy0)),maxval(abs(uz1-uz0))))
  end subroutine

  subroutine stage6_micro_pressure_production_status(pressure_poisson_modified_flag, pressure_rhs_modified_flag, &
                                                     projection_modified_flag, real_projection_called_flag, &
                                                     post_projection_velocity_modified_flag, production_dns_called_flag, &
                                                     production_enabled_by_default_flag)
    integer, intent(out) :: pressure_poisson_modified_flag, pressure_rhs_modified_flag
    integer, intent(out) :: projection_modified_flag, real_projection_called_flag
    integer, intent(out) :: post_projection_velocity_modified_flag, production_dns_called_flag, production_enabled_by_default_flag
    pressure_poisson_modified_flag = 0
    pressure_rhs_modified_flag = 0
    projection_modified_flag = 0
    real_projection_called_flag = 0
    post_projection_velocity_modified_flag = 0
    production_dns_called_flag = 0
    production_enabled_by_default_flag = 0
  end subroutine
end module fibre_stage6_micro_smoke
