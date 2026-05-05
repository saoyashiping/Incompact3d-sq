module fibre_stage5_spreading_candidate
  use fibre_parameters, only : mytype
  use fibre_ibm_types
  use fibre_ibm_grid
  use fibre_ibm_spreading
  use fibre_ibm_force_buffer
  use fibre_ibm_power_diagnostics
  use fibre_ibm_interpolation
  use fibre_ibm_boundary_safety
  use fibre_stage4_grid_adapter
  use fibre_stage4_interpolation_adapter
  use fibre_stage4_boundary_policy
  use fibre_stage5_config
  use fibre_stage5_rhs_wrapper
  implicit none
contains
  subroutine compute_stage5_rhs_candidate_from_buffer(config, fx, fy, fz, ax, ay, az, status)
    type(stage5_config_t), intent(in) :: config
    real(mytype), intent(in) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
    real(mytype), intent(out) :: ax(:,:,:), ay(:,:,:), az(:,:,:)
    integer, intent(out) :: status
    if (config%rho_fluid <= 0._mytype) then
      status = 0; ax = 0._mytype; ay = 0._mytype; az = 0._mytype
    else
      ax = fx / config%rho_fluid
      ay = fy / config%rho_fluid
      az = fz / config%rho_fluid
      status = 1
    end if
  end subroutine

  subroutine compute_stage5_rhs_candidate_expected_error(config, fx, fy, fz, ax, ay, az, err, ex, ey, ez)
    type(stage5_config_t), intent(in) :: config
    real(mytype), intent(in) :: fx(:,:,:), fy(:,:,:), fz(:,:,:), ax(:,:,:), ay(:,:,:), az(:,:,:)
    real(mytype), intent(out) :: err, ex, ey, ez
    if (config%rho_fluid <= 0._mytype) then
      ex=0._mytype; ey=0._mytype; ez=0._mytype; err=0._mytype
    else
      ex = maxval(abs(ax - fx/config%rho_fluid))
      ey = maxval(abs(ay - fy/config%rho_fluid))
      ez = maxval(abs(az - fz/config%rho_fluid))
      err = max(ex,max(ey,ez))
    end if
  end subroutine

  subroutine compute_stage5_real_rhs_noop_change(rhsx0,rhsy0,rhsz0,rhsx1,rhsy1,rhsz1,change_max)
    real(mytype), intent(in) :: rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:),rhsx1(:,:,:),rhsy1(:,:,:),rhsz1(:,:,:)
    real(mytype), intent(out) :: change_max
    change_max = max(maxval(abs(rhsx1-rhsx0)),max(maxval(abs(rhsy1-rhsy0)),maxval(abs(rhsz1-rhsz0))))
  end subroutine
end module
