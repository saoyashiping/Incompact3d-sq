module fibre_stage8_config
  use fibre_parameters, only: mytype
  implicit none

  type stage8_config_t
    integer :: enable_stage8
    integer :: controlled_test_enabled
    integer :: production_dns_enabled
    integer :: production_two_way_enabled

    integer :: enable_runtime_grid_bridge
    integer :: enable_lagrangian_state
    integer :: enable_velocity_to_fibre
    integer :: enable_feedback_candidate
    integer :: enable_force_density_candidate
    integer :: enable_rhs_candidate

    real(mytype) :: rho_fluid
  end type stage8_config_t

contains

  subroutine init_stage8_default_config(cfg)
    type(stage8_config_t), intent(out) :: cfg

    cfg%enable_stage8 = 0
    cfg%controlled_test_enabled = 0
    cfg%production_dns_enabled = 0
    cfg%production_two_way_enabled = 0

    cfg%enable_runtime_grid_bridge = 0
    cfg%enable_lagrangian_state = 0
    cfg%enable_velocity_to_fibre = 0
    cfg%enable_feedback_candidate = 0
    cfg%enable_force_density_candidate = 0
    cfg%enable_rhs_candidate = 0

    cfg%rho_fluid = 1.0_mytype
  end subroutine init_stage8_default_config

  subroutine init_stage8_controlled_test_config(cfg)
    type(stage8_config_t), intent(out) :: cfg

    cfg%enable_stage8 = 1
    cfg%controlled_test_enabled = 1
    cfg%production_dns_enabled = 0
    cfg%production_two_way_enabled = 0

    cfg%enable_runtime_grid_bridge = 1
    cfg%enable_lagrangian_state = 1
    cfg%enable_velocity_to_fibre = 1
    cfg%enable_feedback_candidate = 0
    cfg%enable_force_density_candidate = 0
    cfg%enable_rhs_candidate = 0

    cfg%rho_fluid = 1.0_mytype
  end subroutine init_stage8_controlled_test_config

  subroutine validate_stage8_config(cfg, valid, rejected, rhs_allowed)
    type(stage8_config_t), intent(in) :: cfg
    integer, intent(out) :: valid, rejected, rhs_allowed

    valid = 0
    rejected = 0
    rhs_allowed = 0

    if (cfg%rho_fluid <= 0._mytype) then
      rejected = 1
      return
    end if

    if (cfg%production_dns_enabled == 1) then
      rejected = 1
      return
    end if

    if (cfg%production_two_way_enabled == 1) then
      rejected = 1
      return
    end if

    if (cfg%enable_stage8 == 0) then
      valid = 1
      rhs_allowed = 0
      return
    end if

    if (cfg%controlled_test_enabled == 0) then
      rejected = 1
      return
    end if

    valid = 1
    if (cfg%enable_rhs_candidate == 1 .and. cfg%controlled_test_enabled == 1 .and. &
        cfg%production_dns_enabled == 0 .and. cfg%production_two_way_enabled == 0) then
      rhs_allowed = 1
    end if
  end subroutine validate_stage8_config

end module fibre_stage8_config
