program fibre_stage10_noop_hook_check
  use fibre_stage10_noop_hook, only : stage10_hook_init, stage10_hook_pre_step, stage10_hook_pre_rhs, &
                                       stage10_hook_post_projection, stage10_hook_post_step, &
                                       stage10_hook_finalize, stage10_hook_get_status_values
  implicit none

  integer :: requested_flag, noop_mode_status
  integer :: init_status, pre_step_status, pre_rhs_status
  integer :: post_projection_status, post_step_status, finalize_status
  integer :: no_fibre_state_status, no_force_status, no_rhs_injection_status
  integer :: no_ibm_call_status, no_structure_advance_status
  integer :: hook_build_status
  integer :: unit_dat

  call stage10_hook_init()
  call stage10_hook_pre_step()
  call stage10_hook_pre_rhs()
  call stage10_hook_post_projection()
  call stage10_hook_post_step()
  call stage10_hook_finalize()

  call stage10_hook_get_status_values(requested_flag, noop_mode_status, init_status, pre_step_status, &
                                      pre_rhs_status, post_projection_status, post_step_status, finalize_status, &
                                      no_fibre_state_status, no_force_status, no_rhs_injection_status, &
                                      no_ibm_call_status, no_structure_advance_status)

  hook_build_status = 1
  if (requested_flag /= 1) hook_build_status = 0
  if (noop_mode_status /= 1) hook_build_status = 0
  if (init_status /= 1 .or. pre_step_status /= 1 .or. pre_rhs_status /= 1) hook_build_status = 0
  if (post_projection_status /= 1 .or. post_step_status /= 1 .or. finalize_status /= 1) hook_build_status = 0
  if (no_fibre_state_status /= 1 .or. no_force_status /= 1) hook_build_status = 0
  if (no_rhs_injection_status /= 1 .or. no_ibm_call_status /= 1 .or. no_structure_advance_status /= 1) then
    hook_build_status = 0
  endif

  call execute_command_line('mkdir -p stage10_outputs')
  open(newunit=unit_dat, file='stage10_outputs/fibre_stage10_1_noop_hook.dat', status='replace', action='write')
  write(unit_dat,*) 'stage10_1_requested_flag ', requested_flag
  write(unit_dat,*) 'stage10_1_noop_mode_status ', noop_mode_status
  write(unit_dat,*) 'stage10_1_hook_init_status ', init_status
  write(unit_dat,*) 'stage10_1_hook_pre_step_status ', pre_step_status
  write(unit_dat,*) 'stage10_1_hook_pre_rhs_status ', pre_rhs_status
  write(unit_dat,*) 'stage10_1_hook_post_projection_status ', post_projection_status
  write(unit_dat,*) 'stage10_1_hook_post_step_status ', post_step_status
  write(unit_dat,*) 'stage10_1_hook_finalize_status ', finalize_status
  write(unit_dat,*) 'stage10_1_no_fibre_state_status ', no_fibre_state_status
  write(unit_dat,*) 'stage10_1_no_force_status ', no_force_status
  write(unit_dat,*) 'stage10_1_no_rhs_injection_status ', no_rhs_injection_status
  write(unit_dat,*) 'stage10_1_no_ibm_call_status ', no_ibm_call_status
  write(unit_dat,*) 'stage10_1_no_structure_advance_status ', no_structure_advance_status
  write(unit_dat,*) 'stage10_1_hook_build_status ', hook_build_status
  close(unit_dat)

  if (hook_build_status == 1) then
    write(*,'(A)') 'STAGE 10.1 HOOK BUILD VERDICT: PASS'
  else
    write(*,'(A)') 'STAGE 10.1 HOOK BUILD VERDICT: FAIL'
  endif
end program fibre_stage10_noop_hook_check
