module fibre_stage4_diagnostics

  use fibre_parameters, only : mytype
  use fibre_stage4_config, only : stage4_oneway_config_t

  implicit none
  private

  public :: stage4_safety_flags
  public :: write_stage4_config_diagnostics

contains

  subroutine stage4_safety_flags(config, rhs_disabled_flag, diagnostics_enabled_flag)
    type(stage4_oneway_config_t), intent(in) :: config
    integer, intent(out) :: rhs_disabled_flag
    integer, intent(out) :: diagnostics_enabled_flag

    if (.not. config%apply_ibm_to_fluid_rhs) then
      rhs_disabled_flag = 1
    else
      rhs_disabled_flag = 0
    end if

    if (config%enable_ibm_diagnostics) then
      diagnostics_enabled_flag = 1
    else
      diagnostics_enabled_flag = 0
    end if
  end subroutine stage4_safety_flags

  subroutine write_stage4_config_diagnostics(config, filename)
    type(stage4_oneway_config_t), intent(in) :: config
    character(len=*), intent(in) :: filename

    integer :: unit_id

    open(newunit=unit_id, file=trim(filename), status='replace', action='write', form='formatted')
    write(unit_id,'(A,1X,I0)') 'stage4_enable_oneway', merge(1, 0, config%enable_stage4_oneway)
    write(unit_id,'(A,1X,I0)') 'stage4_enable_ibm_diagnostics', merge(1, 0, config%enable_ibm_diagnostics)
    write(unit_id,'(A,1X,I0)') 'stage4_apply_ibm_to_fluid_rhs', merge(1, 0, config%apply_ibm_to_fluid_rhs)
    write(unit_id,'(A,1X,I0)') 'stage4_coupling_mode', config%coupling_mode
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_beta_drag', config%beta_drag
    write(unit_id,'(A,1X,ES24.16E3)') 'stage4_structure_dt', config%structure_dt
    close(unit_id)
  end subroutine write_stage4_config_diagnostics

end module fibre_stage4_diagnostics
