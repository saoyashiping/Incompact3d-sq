module fibre_stage10_noop_hook
  use mpi
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
    call stage10_hook_write_stage10_3_dat()
  end subroutine stage10_hook_finalize


  logical function stage10_hook_stage10_3_diag_requested() result(requested)
    character(len=64) :: env_value
    integer :: env_status

    requested = .false.
    env_value = ''
    call get_environment_variable('X3D_STAGE10_3_MAIN_NOOP_HOOK', value=env_value, status=env_status)
    if (env_status==0) requested = (trim(adjustl(env_value))=='1')
  end function stage10_hook_stage10_3_diag_requested

  subroutine stage10_hook_write_stage10_3_dat()
    integer :: u, ierr, mpi_rank
    logical :: mpi_is_initialized
    integer :: main_noop_hook_status

    if (.not. stage10_hook_stage10_3_diag_requested()) return

    mpi_rank = 0
    mpi_is_initialized = .false.
    ierr = 0
    call MPI_Initialized(mpi_is_initialized, ierr)
    if (ierr==0 .and. mpi_is_initialized) then
      call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
      if (ierr/=0) mpi_rank = 0
    endif
    if (mpi_rank/=0) return

    main_noop_hook_status = min(requested_flag, min(noop_mode_status, min(initialized_status, &
      min(pre_step_status, min(pre_rhs_status, min(post_projection_status, &
      min(post_step_status, min(finalize_status, min(no_fibre_state_status, &
      min(no_force_status, min(no_rhs_injection_status, min(no_ibm_call_status, &
      min(no_structure_advance_status, merge(1,0,field_modified_status==0))))))))))))))

    call execute_command_line('mkdir -p stage10_outputs', wait=.true.)
    open(newunit=u, file='stage10_outputs/fibre_stage10_3_main_noop_hook.dat', status='replace', action='write')
    write(u,*) 'stage10_3_requested_flag ', requested_flag
    write(u,*) 'stage10_3_noop_mode_status ', noop_mode_status
    write(u,*) 'stage10_3_hook_init_status ', initialized_status
    write(u,*) 'stage10_3_hook_pre_step_status ', pre_step_status
    write(u,*) 'stage10_3_hook_pre_rhs_status ', pre_rhs_status
    write(u,*) 'stage10_3_hook_post_projection_status ', post_projection_status
    write(u,*) 'stage10_3_hook_post_step_status ', post_step_status
    write(u,*) 'stage10_3_hook_finalize_status ', finalize_status
    write(u,*) 'stage10_3_no_fibre_state_status ', no_fibre_state_status
    write(u,*) 'stage10_3_no_force_status ', no_force_status
    write(u,*) 'stage10_3_no_rhs_injection_status ', no_rhs_injection_status
    write(u,*) 'stage10_3_no_ibm_call_status ', no_ibm_call_status
    write(u,*) 'stage10_3_no_structure_advance_status ', no_structure_advance_status
    write(u,*) 'stage10_3_field_modified_status ', field_modified_status
    write(u,*) 'stage10_3_main_noop_hook_status ', main_noop_hook_status
    close(u)
  end subroutine stage10_hook_write_stage10_3_dat

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
