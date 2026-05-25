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
  integer, save :: field_modified_status = 0

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
    call stage10_write_main_noop_hook_dat_if_requested()
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

  subroutine stage10_write_main_noop_hook_dat_if_requested()
    character(len=64) :: value
    integer :: env_status
    integer :: unit_id
    integer :: main_status

    call get_environment_variable('X3D_STAGE10_3_MAIN_NOOP_HOOK', value=value, status=env_status)
    if (env_status /= 0) return
    if (trim(value) /= '1') return
    if (.not. stage10_rank0_from_environment()) return

    main_status = 0
    if (requested_flag == 1 .and. noop_mode_status == 1 .and. &
        initialized_status == 1 .and. pre_step_status == 1 .and. pre_rhs_status == 1 .and. &
        post_projection_status == 1 .and. post_step_status == 1 .and. finalize_status == 1 .and. &
        no_fibre_state_status == 1 .and. no_force_status == 1 .and. no_rhs_injection_status == 1 .and. &
        no_ibm_call_status == 1 .and. no_structure_advance_status == 1 .and. field_modified_status == 0) then
      main_status = 1
    endif

    open(newunit=unit_id, file='stage10_outputs/fibre_stage10_3_main_noop_hook.dat', status='replace', action='write', form='formatted')
    write(unit_id,'(A,1X,I0)') 'stage10_3_requested_flag', requested_flag
    write(unit_id,'(A,1X,I0)') 'stage10_3_noop_mode_status', noop_mode_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_hook_init_status', initialized_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_hook_pre_step_status', pre_step_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_hook_pre_rhs_status', pre_rhs_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_hook_post_projection_status', post_projection_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_hook_post_step_status', post_step_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_hook_finalize_status', finalize_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_no_fibre_state_status', no_fibre_state_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_no_force_status', no_force_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_no_rhs_injection_status', no_rhs_injection_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_no_ibm_call_status', no_ibm_call_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_no_structure_advance_status', no_structure_advance_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_field_modified_status', field_modified_status
    write(unit_id,'(A,1X,I0)') 'stage10_3_main_noop_hook_status', main_status
    close(unit_id)
  end subroutine stage10_write_main_noop_hook_dat_if_requested

  logical function stage10_rank0_from_environment()
    character(len=64) :: value
    integer :: env_status
    integer :: rank_value
    integer :: read_status

    stage10_rank0_from_environment = .true.

    call get_environment_variable('OMPI_COMM_WORLD_RANK', value=value, status=env_status)
    if (env_status /= 0) call get_environment_variable('PMI_RANK', value=value, status=env_status)
    if (env_status /= 0) call get_environment_variable('MPI_RANKID', value=value, status=env_status)
    if (env_status /= 0) return

    read(value, *, iostat=read_status) rank_value
    if (read_status /= 0) return
    stage10_rank0_from_environment = (rank_value == 0)
  end function stage10_rank0_from_environment

end module fibre_stage10_noop_hook
