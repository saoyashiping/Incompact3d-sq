program fibre_stage16_structure_force_input_check
  use decomp_2d_constants, only : mytype
  use fibre_stage16_structure_force_input, only : stage16_structure_force_input_load_from_environment, &
       stage16_structure_force_input_reset, stage16_structure_force_input_set_from_stage12_candidate, &
       stage16_structure_force_input_get_force, stage16_structure_force_input_write_diagnostics, &
       stage16_structure_force_input_get_final_status
  implicit none

  integer :: io_unit
  integer :: final_status
  integer :: wrong_sign_rejected_status
  integer :: reset_clear_status
  integer :: readback_status
  real(mytype) :: u_f(3)
  real(mytype) :: v_f(3)
  real(mytype) :: force_readback(3)
  real(mytype) :: alpha
  character(len=256) :: env_value
  integer :: env_status
  integer :: read_status

  call execute_command_line('mkdir -p stage16_outputs')
  call stage16_structure_force_input_load_from_environment()

  alpha = 1.0_mytype
  call get_environment_variable('STAGE16_4_FEEDBACK_ALPHA', value=env_value, status=env_status)
  if (env_status == 0) then
    read(env_value, *, iostat=read_status) alpha
    if (read_status /= 0) alpha = 1.0_mytype
  end if

  u_f(:) = (/ 2.0e-7_mytype, -1.0e-7_mytype, 1.0e-7_mytype /)
  v_f(:) = (/ 1.0e-7_mytype, -2.0e-7_mytype, 0.0_mytype /)

  call stage16_structure_force_input_reset()
  call stage16_structure_force_input_get_force(1, force_readback)
  reset_clear_status = logical_to_int(maxval(abs(force_readback)) == 0.0_mytype)

  call stage16_structure_force_input_load_from_environment()
  call stage16_structure_force_input_set_from_stage12_candidate(u_f, v_f, alpha, .true.)
  wrong_sign_rejected_status = logical_to_int(stage16_structure_force_input_get_final_status() == 0)

  call stage16_structure_force_input_load_from_environment()
  call stage16_structure_force_input_set_from_stage12_candidate(u_f, v_f, alpha, .false.)
  call stage16_structure_force_input_get_force(1, force_readback)
  readback_status = logical_to_int(maxval(abs(force_readback - alpha * (u_f - v_f))) <= 1.0e-14_mytype)
  final_status = logical_to_int(stage16_structure_force_input_get_final_status() == 1 .and. &
       wrong_sign_rejected_status == 1 .and. reset_clear_status == 1 .and. readback_status == 1)

  if (rank0_write_allowed()) then
    open(newunit=io_unit, file='stage16_outputs/fibre_stage16_4_structure_force_input.dat', &
         status='replace', action='write')
    call stage16_structure_force_input_write_diagnostics(io_unit)
    write(io_unit,'(A,1X,I0)') 'force_input_reset_status', reset_clear_status
    write(io_unit,'(A,1X,I0)') 'force_input_readback_status', readback_status
    write(io_unit,'(A,1X,I0)') 'wrong_sign_rejection_status', wrong_sign_rejected_status
    write(io_unit,'(A,1X,I0)') 'final_status', final_status
    close(io_unit)
  end if

  if (final_status /= 1) stop 1

contains

  integer function logical_to_int(value)
    logical, intent(in) :: value
    if (value) then
      logical_to_int = 1
    else
      logical_to_int = 0
    end if
  end function logical_to_int

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

end program fibre_stage16_structure_force_input_check
