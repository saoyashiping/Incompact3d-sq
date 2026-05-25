module fibre_stage10_noop_hook
  use fibre_stage10_config, only : stage10_config_load, stage10_requested, stage10_noop_mode
  implicit none
  private

  integer, save :: initialized_status = 0
  integer, save :: pre_step_status = 0
  integer, save :: pre_rhs_status = 0
  integer, save :: post_projection_status = 0
  integer, save :: post_step_status = 0
  integer, save :: finalize_status = 0
  integer, save :: requested_flag = 0
  integer, save :: noop_mode_status = 1
  integer, save :: no_fibre_state_status = 1
  integer, save :: no_force_status = 1
  integer, save :: no_rhs_injection_status = 1
  integer, save :: no_ibm_call_status = 1
  integer, save :: no_structure_advance_status = 1

  public :: stage10_hook_init
  public :: stage10_hook_pre_step
  public :: stage10_hook_pre_rhs
  public :: stage10_hook_post_projection
  public :: stage10_hook_post_step
  public :: stage10_hook_finalize
  public :: stage10_hook_get_status_values

contains

  subroutine stage10_hook_init()
    call stage10_config_load()

    if (stage10_requested()) then
      requested_flag = 1
    else
      requested_flag = 0
    endif

    if (stage10_noop_mode()) then
      noop_mode_status = 1
    else
      noop_mode_status = 0
    endif

    initialized_status = 1
  end subroutine stage10_hook_init

  subroutine stage10_hook_pre_step()
    pre_step_status = 1
  end subroutine stage10_hook_pre_step

  subroutine stage10_hook_pre_rhs()
    pre_rhs_status = 1
  end subroutine stage10_hook_pre_rhs

  subroutine stage10_hook_post_projection()
    post_projection_status = 1
  end subroutine stage10_hook_post_projection

  subroutine stage10_hook_post_step()
    post_step_status = 1
  end subroutine stage10_hook_post_step

  subroutine stage10_hook_finalize()
    finalize_status = 1
  end subroutine stage10_hook_finalize

  subroutine stage10_hook_get_status_values(requested, noop_mode, init_status, pre_step, pre_rhs, &
                                            post_projection, post_step, finalized, no_fibre_state, &
                                            no_force, no_rhs_injection, no_ibm_call, no_structure_advance)
    integer, intent(out) :: requested, noop_mode, init_status, pre_step, pre_rhs
    integer, intent(out) :: post_projection, post_step, finalized
    integer, intent(out) :: no_fibre_state, no_force, no_rhs_injection
    integer, intent(out) :: no_ibm_call, no_structure_advance

    requested = requested_flag
    noop_mode = noop_mode_status
    init_status = initialized_status
    pre_step = pre_step_status
    pre_rhs = pre_rhs_status
    post_projection = post_projection_status
    post_step = post_step_status
    finalized = finalize_status
    no_fibre_state = no_fibre_state_status
    no_force = no_force_status
    no_rhs_injection = no_rhs_injection_status
    no_ibm_call = no_ibm_call_status
    no_structure_advance = no_structure_advance_status
  end subroutine stage10_hook_get_status_values

end module fibre_stage10_noop_hook
