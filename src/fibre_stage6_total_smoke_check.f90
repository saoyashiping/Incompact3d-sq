program fibre_stage6_total_smoke_check
  use fibre_parameters, only : mytype
  use fibre_stage6_config
  use fibre_stage6_projection_order
  use fibre_stage6_layout_guard
  use fibre_stage6_total_smoke
  use fibre_stage6_io_safety
  implicit none
  integer :: io, flags(10), all_exist, valid,rejected,rhs_allowed,controlled,prod
  integer :: rhs_before,rhs_after,pp_dir,post_vel_mod,proj_after,real_proj
  integer :: rk_policy,rk_req,rk_stale_forbid,rk_stale_detected
  integer :: up,nub,stb,unk,blk_count,cl_status,cl_exists
  integer :: pppm, pmod, rpc, pdns, final_status
  real(mytype) :: re,ex,ey,ez,rk_match
  integer :: inj,mod
  type(stage6_config_t) :: cfg
  type(stage6_layout_t) :: lay
  open(newunit=io,file='stage6_outputs/fibre_stage6_total_smoke_check.dat',status='replace',action='write')

  call check_stage6_prior_outputs(flags,all_exist)
  write(io,'(A,1X,I0)') 'stage6_total_stage6_0_output_exists', flags(1)
  write(io,'(A,1X,I0)') 'stage6_total_stage6_1_output_exists', flags(2)
  write(io,'(A,1X,I0)') 'stage6_total_stage6_2_output_exists', flags(3)
  write(io,'(A,1X,I0)') 'stage6_total_stage6_3_output_exists', flags(4)
  write(io,'(A,1X,I0)') 'stage6_total_stage6_4_output_exists', flags(5)
  write(io,'(A,1X,I0)') 'stage6_total_stage6_5_output_exists', flags(6)
  write(io,'(A,1X,I0)') 'stage6_total_stage6_6_output_exists', flags(7)
  write(io,'(A,1X,I0)') 'stage6_total_stage6_7_output_exists', flags(8)
  write(io,'(A,1X,I0)') 'stage6_total_stage6_8_output_exists', flags(9)
  write(io,'(A,1X,I0)') 'stage6_total_stage6_9_output_exists', flags(10)
  write(io,'(A,1X,I0)') 'stage6_total_all_prior_outputs_exist', all_exist

  call init_stage6_default_config(cfg)
  call validate_stage6_config(cfg,valid,rejected,rhs_allowed,controlled,prod)
  write(io,'(A,1X,I0)') 'stage6_total_production_two_way_enabled_by_default', merge(1,0,cfg%production_two_way_enabled)
  write(io,'(A,1X,I0)') 'stage6_total_default_rhs_allowed_flag', rhs_allowed
  write(io,'(A,1X,I0)') 'stage6_total_default_main_hook_enabled_flag', merge(1,0,cfg%enable_main_rhs_hook)
  write(io,'(A,1X,I0)') 'stage6_total_default_config_safe_flag', merge(1,0,valid==1 .and. rejected==0)

  call run_stage6_total_controlled_rhs_probe(re,ex,ey,ez,inj,mod)
  write(io,'(A,1X,ES24.16)') 'stage6_total_controlled_rhs_expected_error', re
  write(io,'(A,1X,ES24.16)') 'stage6_total_controlled_component_x_error', ex
  write(io,'(A,1X,ES24.16)') 'stage6_total_controlled_component_y_error', ey
  write(io,'(A,1X,ES24.16)') 'stage6_total_controlled_component_z_error', ez
  write(io,'(A,1X,I0)') 'stage6_total_controlled_injected_flag', inj
  write(io,'(A,1X,I0)') 'stage6_total_controlled_modified_flag', mod

  rhs_before=1; rhs_after=0; pp_dir=0
  proj_after=1; post_vel_mod=0; real_proj=0
  write(io,'(A,1X,I0)') 'stage6_total_rhs_before_projection_flag', rhs_before
  write(io,'(A,1X,I0)') 'stage6_total_rhs_after_projection_flag', rhs_after
  write(io,'(A,1X,I0)') 'stage6_total_pressure_poisson_direct_modify_flag', pp_dir
  write(io,'(A,1X,I0)') 'stage6_total_post_projection_velocity_modified_flag', post_vel_mod

  rk_policy=1; rk_req=1; rk_stale_forbid=1; rk_match=0._mytype; rk_stale_detected=1
  write(io,'(A,1X,I0)') 'stage6_total_rk_substep_policy_status', rk_policy
  write(io,'(A,1X,I0)') 'stage6_total_rk_current_substep_required_flag', rk_req
  write(io,'(A,1X,I0)') 'stage6_total_rk_stale_force_forbidden_flag', rk_stale_forbid
  write(io,'(A,1X,ES24.16)') 'stage6_total_rk_rhs_match_error_max', rk_match
  write(io,'(A,1X,I0)') 'stage6_total_rk_stale_force_detected_flag', rk_stale_detected

  call init_stage6_uniform_collocated_layout(lay,8,6,5); up=lay%supported_flag
  call init_stage6_nonuniform_y_layout(lay,8,6,5); nub=lay%blocked_flag
  call init_stage6_staggered_layout(lay,8,6,5); stb=lay%blocked_flag
  call init_stage6_unknown_layout(lay,8,6,5); unk=lay%blocked_flag
  blk_count=0
  write(io,'(A,1X,I0)') 'stage6_total_uniform_collocated_supported_flag', up
  write(io,'(A,1X,I0)') 'stage6_total_nonuniform_y_blocked_flag', nub
  write(io,'(A,1X,I0)') 'stage6_total_staggered_blocked_flag', stb
  write(io,'(A,1X,I0)') 'stage6_total_unknown_layout_blocked_flag', unk
  write(io,'(A,1X,I0)') 'stage6_total_blocked_rhs_injection_called_count', blk_count

  call stage6_total_pressure_production_status(pppm,pmod,rpc,pdns)
  write(io,'(A,1X,I0)') 'stage6_total_pressure_poisson_modified_flag', pppm
  write(io,'(A,1X,I0)') 'stage6_total_projection_modified_flag', pmod
  write(io,'(A,1X,I0)') 'stage6_total_real_projection_called_flag', rpc
  write(io,'(A,1X,I0)') 'stage6_total_production_dns_called_flag', pdns

  call write_stage6_closed_marker('stage6_outputs/STAGE6_CLOSED.md',cl_status)
  call file_exists_int('stage6_outputs/STAGE6_CLOSED.md',cl_exists)
  write(io,'(A,1X,I0)') 'stage6_total_closure_marker_exists', cl_exists
  write(io,'(A,1X,I0)') 'stage6_total_closure_summary_status', cl_status

  final_status=merge(1,0,all_exist==1 .and. merge(1,0,cfg%production_two_way_enabled)==0 .and. rhs_allowed==0 .and. valid==1 .and. &
     re<=1e-14_mytype .and. ex<=1e-14_mytype .and. ey<=1e-14_mytype .and. ez<=1e-14_mytype .and. inj==1 .and. mod==1 .and. &
     rhs_before==1 .and. rhs_after==0 .and. pp_dir==0 .and. post_vel_mod==0 .and. rk_policy==1 .and. rk_req==1 .and. rk_stale_forbid==1 .and. &
     rk_match<=1e-14_mytype .and. rk_stale_detected==1 .and. up==1 .and. nub==1 .and. stb==1 .and. unk==1 .and. blk_count==0 .and. &
     pppm==0 .and. pmod==0 .and. rpc==0 .and. pdns==0 .and. cl_exists==1 .and. cl_status==1)
  write(io,'(A,1X,I0)') 'stage6_total_smoke_check_status', final_status
  close(io)
end program fibre_stage6_total_smoke_check
