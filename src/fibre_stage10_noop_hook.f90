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
    if (stage10_main_noop_diag_requested()) then
      if (stage10_rank_is_writer()) call stage10_write_main_noop_diag()
    endif
  end subroutine stage10_hook_finalize

  logical function stage10_main_noop_diag_requested()
    character(len=64) :: value
    integer :: status

    value = ''
    call get_environment_variable('X3D_STAGE10_3_MAIN_NOOP_HOOK', value=value, status=status)
    stage10_main_noop_diag_requested = .false.
    if (status == 0) then
      select case (trim(adjustl(value)))
      case ('1', 'true', 'TRUE', 'yes', 'YES', 'on', 'ON')
        stage10_main_noop_diag_requested = .true.
      end select
    endif
  end function stage10_main_noop_diag_requested

  logical function stage10_rank_is_writer()
    character(len=64) :: value
    integer :: status, ios, rank_value

    stage10_rank_is_writer = .true.

    value = ''
    call get_environment_variable('OMPI_COMM_WORLD_RANK', value=value, status=status)
    if (status /= 0) call get_environment_variable('PMI_RANK', value=value, status=status)
    if (status /= 0) call get_environment_variable('PMIX_RANK', value=value, status=status)
    if (status /= 0) call get_environment_variable('SLURM_PROCID', value=value, status=status)

    if (status == 0) then
      read(value, *, iostat=ios) rank_value
      if (ios == 0) stage10_rank_is_writer = (rank_value == 0)
    endif
  end function stage10_rank_is_writer

  subroutine stage10_write_main_noop_diag()
    integer :: unit_dat
    integer :: main_status

    main_status = 1
    if (requested_flag /= 1) main_status = 0
    if (noop_mode_status /= 1) main_status = 0
    if (initialized_status /= 1) main_status = 0
    if (pre_step_status /= 1) main_status = 0
    if (pre_rhs_status /= 1) main_status = 0
    if (post_projection_status /= 1) main_status = 0
    if (post_step_status /= 1) main_status = 0
    if (finalize_status /= 1) main_status = 0
    if (no_fibre_state_status /= 1) main_status = 0
    if (no_force_status /= 1) main_status = 0
    if (no_rhs_injection_status /= 1) main_status = 0
    if (no_ibm_call_status /= 1) main_status = 0
    if (no_structure_advance_status /= 1) main_status = 0

    call execute_command_line('mkdir -p stage10_outputs')
    open(newunit=unit_dat, file='stage10_outputs/fibre_stage10_3_main_noop_hook.dat', &
         status='replace', action='write')
    write(unit_dat,'(A,I0)') 'stage10_3_requested_flag ', requested_flag
    write(unit_dat,'(A,I0)') 'stage10_3_noop_mode_status ', noop_mode_status
    write(unit_dat,'(A,I0)') 'stage10_3_hook_init_status ', initialized_status
    write(unit_dat,'(A,I0)') 'stage10_3_hook_pre_step_status ', pre_step_status
    write(unit_dat,'(A,I0)') 'stage10_3_hook_pre_rhs_status ', pre_rhs_status
    write(unit_dat,'(A,I0)') 'stage10_3_hook_post_projection_status ', post_projection_status
    write(unit_dat,'(A,I0)') 'stage10_3_hook_post_step_status ', post_step_status
    write(unit_dat,'(A,I0)') 'stage10_3_hook_finalize_status ', finalize_status
    write(unit_dat,'(A,I0)') 'stage10_3_no_fibre_state_status ', no_fibre_state_status
    write(unit_dat,'(A,I0)') 'stage10_3_no_force_status ', no_force_status
    write(unit_dat,'(A,I0)') 'stage10_3_no_rhs_injection_status ', no_rhs_injection_status
    write(unit_dat,'(A,I0)') 'stage10_3_no_ibm_call_status ', no_ibm_call_status
    write(unit_dat,'(A,I0)') 'stage10_3_no_structure_advance_status ', no_structure_advance_status
    write(unit_dat,'(A,I0)') 'stage10_3_field_modified_status ', 0
    write(unit_dat,'(A,I0)') 'stage10_3_main_noop_hook_status ', main_status
    close(unit_dat)
  end subroutine stage10_write_main_noop_diag

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
