program fibre_stage16_closed_loop_dryrun_check
  use fibre_stage16_closed_loop_dryrun, only : stage16_closed_loop_dryrun_load_from_environment, &
       stage16_closed_loop_dryrun_initialize, stage16_closed_loop_dryrun_apply_one_step, &
       stage16_closed_loop_dryrun_write_diagnostics, stage16_closed_loop_dryrun_get_final_status
  implicit none

  integer :: io_unit
  integer :: final_status

  call execute_command_line('mkdir -p stage16_outputs')
  call stage16_closed_loop_dryrun_load_from_environment()
  call stage16_closed_loop_dryrun_initialize()
  call stage16_closed_loop_dryrun_apply_one_step()
  final_status = stage16_closed_loop_dryrun_get_final_status()

  if (rank0_write_allowed()) then
    open(newunit=io_unit, file='stage16_outputs/fibre_stage16_5_closed_loop_dryrun_np1.dat', &
         status='replace', action='write')
    call stage16_closed_loop_dryrun_write_diagnostics(io_unit)
    close(io_unit)
  end if

  if (final_status /= 1) stop 1

contains

  logical function rank0_write_allowed()
    character(len=32) :: rank_env
    integer :: rank_status
    integer :: rank_value
    integer :: parse_status

    rank0_write_allowed = .true.
    call get_environment_variable('OMPI_COMM_WORLD_RANK', value=rank_env, status=rank_status)
    if (rank_status /= 0) call get_environment_variable('PMI_RANK', value=rank_env, status=rank_status)
    if (rank_status /= 0) call get_environment_variable('MPI_LOCALRANKID', value=rank_env, status=rank_status)
    if (rank_status == 0) then
      read(rank_env, *, iostat=parse_status) rank_value
      if (parse_status == 0) rank0_write_allowed = (rank_value == 0)
    end if
  end function rank0_write_allowed

end program fibre_stage16_closed_loop_dryrun_check
