module fibre_stage4_config

  use fibre_parameters, only : mytype

  implicit none
  private

  type, public :: stage4_oneway_config_t
    logical :: enable_stage4_oneway
    logical :: enable_ibm_diagnostics
    logical :: apply_ibm_to_fluid_rhs
    integer :: coupling_mode
    real(mytype) :: beta_drag
    real(mytype) :: structure_dt
  end type stage4_oneway_config_t

  public :: init_stage4_oneway_config
  public :: validate_stage4_oneway_config

contains

  subroutine init_stage4_oneway_config(config)
    type(stage4_oneway_config_t), intent(out) :: config

    config%enable_stage4_oneway = .true.
    config%enable_ibm_diagnostics = .true.
    config%apply_ibm_to_fluid_rhs = .false.
    config%coupling_mode = 1
    config%beta_drag = 10.0_mytype
    config%structure_dt = 1.0e-5_mytype
  end subroutine init_stage4_oneway_config

  subroutine validate_stage4_oneway_config(config, status)
    type(stage4_oneway_config_t), intent(in) :: config
    integer, intent(out) :: status

    status = 1

    if (.not. config%enable_stage4_oneway) status = 0
    if (.not. config%enable_ibm_diagnostics) status = 0
    if (config%apply_ibm_to_fluid_rhs) status = 0
    if (config%coupling_mode /= 1) status = 0
    if (config%beta_drag <= 0._mytype) status = 0
    if (config%structure_dt <= 0._mytype) status = 0
  end subroutine validate_stage4_oneway_config

end module fibre_stage4_config
