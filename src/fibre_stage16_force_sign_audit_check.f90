program fibre_stage16_force_sign_audit_check
  use decomp_2d_constants, only : mytype
  use fibre_stage16_force_sign_audit, only : stage16_force_sign_audit_load_from_environment, &
       stage16_force_sign_audit_compute_reference, stage16_force_sign_audit_write_diagnostics, &
       stage16_force_sign_audit_get_final_status
  implicit none

  integer :: io_unit
  integer :: final_status
  real(mytype) :: u_f(3)
  real(mytype) :: v_f(3)
  real(mytype) :: alpha
  character(len=256) :: env_value
  integer :: env_status
  integer :: read_status

  call execute_command_line('mkdir -p stage16_outputs')
  call stage16_force_sign_audit_load_from_environment()

  alpha = 1.0_mytype
  call get_environment_variable('STAGE16_3_FEEDBACK_ALPHA', value=env_value, status=env_status)
  if (env_status == 0) then
    read(env_value, *, iostat=read_status) alpha
    if (read_status /= 0) alpha = 1.0_mytype
  end if

  u_f(:) = (/ 1.25_mytype, -0.50_mytype, 0.75_mytype /)
  v_f(:) = (/ 0.25_mytype, 0.25_mytype, -0.25_mytype /)
  call stage16_force_sign_audit_compute_reference(u_f, v_f, alpha)
  final_status = stage16_force_sign_audit_get_final_status()

  if (rank0_write_allowed()) then
    open(newunit=io_unit, file='stage16_outputs/fibre_stage16_3_force_sign_audit.dat', &
         status='replace', action='write')
    call stage16_force_sign_audit_write_diagnostics(io_unit)
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

end program fibre_stage16_force_sign_audit_check
