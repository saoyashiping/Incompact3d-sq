program fibre_stage10_config_check
  use fibre_stage10_config, only : stage10_config_load, stage10_get_max_steps, stage10_get_status_values
  implicit none

  integer :: requested_flag
  integer :: noop_mode_status
  integer :: disabled_by_default_status
  integer :: no_fibre_state_status
  integer :: no_force_status
  integer :: no_rhs_injection_status
  integer :: config_status
  integer :: max_steps
  integer :: unit

  call stage10_config_load()
  max_steps=stage10_get_max_steps()
  call stage10_get_status_values(requested_flag, noop_mode_status, disabled_by_default_status, no_fibre_state_status, no_force_status, no_rhs_injection_status)

  config_status=1
  if (noop_mode_status/=1) config_status=0
  if (no_fibre_state_status/=1) config_status=0
  if (no_force_status/=1) config_status=0
  if (no_rhs_injection_status/=1) config_status=0
  if (max_steps<=0) config_status=0

  call execute_command_line('mkdir -p stage10_outputs')
  open(newunit=unit,file='stage10_outputs/fibre_stage10_0_config.dat',status='replace',action='write')
  write(unit,*) 'stage10_0_requested_flag ',requested_flag
  write(unit,*) 'stage10_0_noop_mode_status ',noop_mode_status
  write(unit,*) 'stage10_0_disabled_by_default_status ',disabled_by_default_status
  write(unit,*) 'stage10_0_no_fibre_state_status ',no_fibre_state_status
  write(unit,*) 'stage10_0_no_force_status ',no_force_status
  write(unit,*) 'stage10_0_no_rhs_injection_status ',no_rhs_injection_status
  write(unit,*) 'stage10_0_config_status ',config_status
  close(unit)

  if (config_status==1) then
    write(*,'(A)') 'STAGE 10.0 CONFIG VERDICT: PASS'
  else
    write(*,'(A)') 'STAGE 10.0 CONFIG VERDICT: FAIL'
    if (noop_mode_status/=1) write(*,'(A)') 'reason: noop mode disabled'
    if (no_fibre_state_status/=1) write(*,'(A)') 'reason: fibre state active'
    if (no_force_status/=1) write(*,'(A)') 'reason: force status active'
    if (no_rhs_injection_status/=1) write(*,'(A)') 'reason: rhs injection active'
  endif
end program fibre_stage10_config_check
