program fibre_stage8_total_smoke_check
use fibre_parameters, only: mytype
implicit none
integer::io,final_status
integer::s7m,s7o,s7s,s7c,dep7
integer::e80,e81,e82,e83,e84,e85,e86,e87,e88,e89,all_outputs
integer::s80,s81,s82,s83,s84,s85,s86,s87,s88,s89,all_status
integer::cfg_default_off,cfg_controlled,cfg_invalid_prod,cfg_summary
integer::ctrl_rhs_allowed_raw,default_rhs_allowed_raw
integer::ctrl_enable_status,ctrl_test_enabled_status,ctrl_valid_status,ctrl_not_rejected_status,ctrl_rhs_allowed_status,ctrl_config_status
integer::num_grid,num_state,num_v,num_fb,num_ow,num_tw,num_rk,num_bd,num_int,num_summary
integer::cons_f,cons_comp,cons_ar,cons_pair,cons_fp,cons_fd,cons_rho,cons_nodbl,cons_summary
integer::bd_safe,bd_near,bd_outy,bd_invcoord,bd_invlayout,bd_per,bd_mix,bd_summary
integer::rk_nsub,rk_ctr,rk_ord,rk_cons,rk_blk,rk_summary
integer::rhs_untouched,pp_untouched,proj_untouched,dns_not_called,fluid_not_called,fibre_not_called,noside
integer::marker_written,marker_status
real(mytype)::x
integer::tmpi
call ensure_dir('stage8_outputs')
call rm_file('stage8_outputs/STAGE8_CLOSED.md')
call file_exists_int('stage7_outputs/STAGE7_CLOSED.md',s7m)
call file_exists_int('stage7_outputs/fibre_stage7_total_smoke_check.dat',s7o)
call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_smoke_check_status',s7s)
call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_closed_marker_status',s7c)
dep7=merge(1,0,s7m*s7o*s7s*s7c==1)
call file_nonempty_int('stage8_outputs/fibre_stage8_config_check.dat',e80)
call file_nonempty_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat',e81)
call file_nonempty_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat',e82)
call file_nonempty_int('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat',e83)
call file_nonempty_int('stage8_outputs/fibre_stage8_feedback_candidate_check.dat',e84)
call file_nonempty_int('stage8_outputs/fibre_stage8_oneway_forcing_check.dat',e85)
call file_nonempty_int('stage8_outputs/fibre_stage8_twoway_force_density_check.dat',e86)
call file_nonempty_int('stage8_outputs/fibre_stage8_rk_sync_check.dat',e87)
call file_nonempty_int('stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat',e88)
call file_nonempty_int('stage8_outputs/fibre_stage8_integrated_audit_check.dat',e89)
all_outputs=merge(1,0,e80*e81*e82*e83*e84*e85*e86*e87*e88*e89==1)
call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_check_status',s80)
call get_int('stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat','stage8_runtime_grid_bridge_check_status',s81)
call get_int('stage8_outputs/fibre_stage8_lagrangian_state_check.dat','stage8_lagrangian_state_check_status',s82)
call get_int('stage8_outputs/fibre_stage8_velocity_to_fibre_check.dat','stage8_velocity_to_fibre_check_status',s83)
call get_int('stage8_outputs/fibre_stage8_feedback_candidate_check.dat','stage8_feedback_candidate_check_status',s84)
call get_int('stage8_outputs/fibre_stage8_oneway_forcing_check.dat','stage8_oneway_forcing_check_status',s85)
call get_int('stage8_outputs/fibre_stage8_twoway_force_density_check.dat','stage8_twoway_force_density_check_status',s86)
call get_int('stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_sync_check_status',s87)
call get_int('stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_safe_workflow_check_status',s88)
call get_int('stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_audit_check_status',s89)
all_status=merge(1,0,s80*s81*s82*s83*s84*s85*s86*s87*s88*s89==1)
call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_default_enable_stage8',cfg_default_off); cfg_default_off=merge(1,0,cfg_default_off==0)
call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_controlled_enable_stage8',tmpi); ctrl_enable_status=merge(1,0,tmpi==1)
call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_controlled_controlled_test_enabled',tmpi); ctrl_test_enabled_status=merge(1,0,tmpi==1)
call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_controlled_valid_flag',tmpi); ctrl_valid_status=merge(1,0,tmpi==1)
call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_controlled_rejected_flag',tmpi); ctrl_not_rejected_status=merge(1,0,tmpi==0)
call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_controlled_rhs_allowed_flag',ctrl_rhs_allowed_raw)
call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_default_rhs_allowed_flag',default_rhs_allowed_raw)
ctrl_rhs_allowed_status=merge(1,0,ctrl_rhs_allowed_raw==1)
call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_controlled_config_status',tmpi); ctrl_config_status=merge(1,0,tmpi==1)
cfg_controlled=merge(1,0,ctrl_enable_status*ctrl_test_enabled_status*ctrl_valid_status*ctrl_not_rejected_status*ctrl_rhs_allowed_status*ctrl_config_status==1)
call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_invalid_production_dns_rejected_flag',cfg_invalid_prod)
call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_invalid_production_twoway_rejected_flag',tmpi); cfg_invalid_prod=merge(1,0,cfg_invalid_prod==1.and.tmpi==1)
cfg_summary=merge(1,0,cfg_default_off==1.and.cfg_controlled==1.and.cfg_invalid_prod==1)
num_grid=s81; num_state=s82; num_v=s83; num_fb=s84; num_ow=s85; num_tw=s86; num_rk=s87; num_bd=s88; num_int=s89
num_summary=merge(1,0,num_grid*num_state*num_v*num_fb*num_ow*num_tw*num_rk*num_bd*num_int==1)
call get_int('stage8_outputs/fibre_stage8_twoway_force_density_check.dat','stage8_twoway_force_conservation_status',cons_f)
call get_int('stage8_outputs/fibre_stage8_twoway_force_density_check.dat','stage8_twoway_component_conservation_status',cons_comp)
call get_int('stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_action_reaction_status',cons_ar)
call get_int('stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_pair_power_status',cons_pair)
call get_int('stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_fluid_power_consistency_status',cons_fp)
call get_int('stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_force_density_convention_status',cons_fd)
call get_int('stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_rho_convention_status',cons_rho)
call get_int('stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_double_division_detected_flag',cons_nodbl); cons_nodbl=merge(1,0,cons_nodbl==0)
cons_summary=merge(1,0,cons_f*cons_comp*cons_ar*cons_pair*cons_fp*cons_fd*cons_rho*cons_nodbl==1)
call get_int('stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_safe_workflow_status',bd_safe)
call get_int('stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_nearwall_status',bd_near)
call get_int('stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_outside_y_status',bd_outy)
call get_int('stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_invalid_coord_status',bd_invcoord)
call get_int('stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_invalid_layout_status',bd_invlayout)
call get_int('stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_periodic_status',bd_per)
call get_int('stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_mixed_status',bd_mix)
bd_summary=merge(1,0,bd_safe*bd_near*bd_outy*bd_invcoord*bd_invlayout*bd_per*bd_mix==1)
call get_int('stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_nsubstep',rk_nsub); rk_nsub=merge(1,0,rk_nsub==3)
call get_int('stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_substep_valid_count',rk_ctr); call get_int('stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_substep_rejected_count',tmpi); rk_ctr=merge(1,0,rk_ctr==3.and.tmpi==0)
call get_int('stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_event_order_status',rk_ord)
call get_int('stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_force_conservation_status',rk_cons)
call get_int('stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_blocked_substep_status',rk_blk)
rk_summary=merge(1,0,rk_nsub*rk_ctr*rk_ord*rk_cons*rk_blk==1)
rhs_untouched=all_zero5('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_noop_rhs_modified_flag', &
'stage8_outputs/fibre_stage8_twoway_force_density_check.dat','stage8_twoway_rhs_modified_flag', &
'stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_rhs_modified_flag', &
'stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_rhs_modified_flag', &
'stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_production_rhs_modified_flag')
pp_untouched=all_zero5('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_pressure_poisson_modified_flag', &
'stage8_outputs/fibre_stage8_twoway_force_density_check.dat','stage8_twoway_pressure_poisson_modified_flag', &
'stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_pressure_poisson_modified_flag', &
'stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_pressure_poisson_modified_flag', &
'stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_pressure_poisson_modified_flag')
proj_untouched=all_zero5('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_projection_modified_flag', &
'stage8_outputs/fibre_stage8_twoway_force_density_check.dat','stage8_twoway_projection_modified_flag', &
'stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_projection_modified_flag', &
'stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_projection_modified_flag', &
'stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_projection_modified_flag')
dns_not_called=all_zero5('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_production_dns_called_flag', &
'stage8_outputs/fibre_stage8_twoway_force_density_check.dat','stage8_twoway_production_dns_called_flag', &
'stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_production_dns_called_flag', &
'stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_production_dns_called_flag', &
'stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_production_dns_called_flag')
fluid_not_called=all_zero5('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_fluid_update_called_flag', &
'stage8_outputs/fibre_stage8_twoway_force_density_check.dat','stage8_twoway_fluid_update_called_flag', &
'stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_fluid_update_called_flag', &
'stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_fluid_update_called_flag', &
'stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_fluid_update_called_flag')
fibre_not_called=all_zero5('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_fibre_advance_called_flag', &
'stage8_outputs/fibre_stage8_twoway_force_density_check.dat','stage8_twoway_fibre_advance_called_flag', &
'stage8_outputs/fibre_stage8_rk_sync_check.dat','stage8_rk_fibre_advance_called_flag', &
'stage8_outputs/fibre_stage8_boundary_safe_workflow_check.dat','stage8_boundary_fibre_advance_called_flag', &
'stage8_outputs/fibre_stage8_integrated_audit_check.dat','stage8_integrated_fibre_advance_called_flag')
noside=merge(1,0,rhs_untouched*pp_untouched*proj_untouched*dns_not_called*fluid_not_called*fibre_not_called==1)
final_status=merge(1,0,dep7*all_outputs*all_status*cfg_summary*num_summary*cons_summary*bd_summary*rk_summary*noside==1)
if(final_status==1) then
  call write_marker('stage8_outputs/STAGE8_CLOSED.md')
  marker_written=1; marker_status=1
else
  call rm_file('stage8_outputs/STAGE8_CLOSED.md')
  marker_written=0; marker_status=0
endif
open(newunit=io,file='stage8_outputs/fibre_stage8_total_smoke_check.dat',status='replace',action='write')
write(io,'(A,1X,I0)') 'stage8_total_stage7_closed_marker_exists',s7m
write(io,'(A,1X,I0)') 'stage8_total_stage7_total_smoke_output_exists',s7o
write(io,'(A,1X,I0)') 'stage8_total_stage7_total_smoke_status',s7s
write(io,'(A,1X,I0)') 'stage8_total_stage7_closed_marker_status',s7c
write(io,'(A,1X,I0)') 'stage8_total_stage7_dependency_status',dep7
write(io,'(A,1X,I0)') 'stage8_total_stage8_0_output_exists',e80; write(io,'(A,1X,I0)') 'stage8_total_stage8_1_output_exists',e81; write(io,'(A,1X,I0)') 'stage8_total_stage8_2_output_exists',e82; write(io,'(A,1X,I0)') 'stage8_total_stage8_3_output_exists',e83; write(io,'(A,1X,I0)') 'stage8_total_stage8_4_output_exists',e84; write(io,'(A,1X,I0)') 'stage8_total_stage8_5_output_exists',e85; write(io,'(A,1X,I0)') 'stage8_total_stage8_6_output_exists',e86; write(io,'(A,1X,I0)') 'stage8_total_stage8_7_output_exists',e87; write(io,'(A,1X,I0)') 'stage8_total_stage8_8_output_exists',e88; write(io,'(A,1X,I0)') 'stage8_total_stage8_9_output_exists',e89; write(io,'(A,1X,I0)') 'stage8_total_all_stage8_outputs_exist',all_outputs
write(io,'(A,1X,I0)') 'stage8_total_stage8_0_status',s80; write(io,'(A,1X,I0)') 'stage8_total_stage8_1_status',s81; write(io,'(A,1X,I0)') 'stage8_total_stage8_2_status',s82; write(io,'(A,1X,I0)') 'stage8_total_stage8_3_status',s83; write(io,'(A,1X,I0)') 'stage8_total_stage8_4_status',s84; write(io,'(A,1X,I0)') 'stage8_total_stage8_5_status',s85; write(io,'(A,1X,I0)') 'stage8_total_stage8_6_status',s86; write(io,'(A,1X,I0)') 'stage8_total_stage8_7_status',s87; write(io,'(A,1X,I0)') 'stage8_total_stage8_8_status',s88; write(io,'(A,1X,I0)') 'stage8_total_stage8_9_status',s89; write(io,'(A,1X,I0)') 'stage8_total_all_stage8_status',all_status
write(io,'(A,1X,I0)') 'stage8_total_default_production_disabled_status',cfg_default_off; write(io,'(A,1X,I0)') 'stage8_total_controlled_path_status',cfg_controlled; write(io,'(A,1X,I0)') 'stage8_total_invalid_production_rejection_status',cfg_invalid_prod; write(io,'(A,1X,I0)') 'stage8_total_config_summary_status',cfg_summary
write(io,'(A,1X,I0)') 'stage8_total_controlled_enable_status',ctrl_enable_status
write(io,'(A,1X,I0)') 'stage8_total_controlled_test_enabled_status',ctrl_test_enabled_status
write(io,'(A,1X,I0)') 'stage8_total_controlled_valid_status',ctrl_valid_status
write(io,'(A,1X,I0)') 'stage8_total_controlled_not_rejected_status',ctrl_not_rejected_status
write(io,'(A,1X,I0)') 'stage8_total_controlled_rhs_allowed_status',ctrl_rhs_allowed_status
write(io,'(A,1X,I0)') 'stage8_total_controlled_config_status',ctrl_config_status
write(io,'(A,1X,I0)') 'stage8_total_grid_numeric_status',num_grid; write(io,'(A,1X,I0)') 'stage8_total_lagrangian_state_numeric_status',num_state; write(io,'(A,1X,I0)') 'stage8_total_velocity_interpolation_numeric_status',num_v; write(io,'(A,1X,I0)') 'stage8_total_feedback_numeric_status',num_fb; write(io,'(A,1X,I0)') 'stage8_total_oneway_numeric_status',num_ow; write(io,'(A,1X,I0)') 'stage8_total_twoway_numeric_status',num_tw; write(io,'(A,1X,I0)') 'stage8_total_rk_sync_numeric_status',num_rk; write(io,'(A,1X,I0)') 'stage8_total_boundary_safe_numeric_status',num_bd; write(io,'(A,1X,I0)') 'stage8_total_integrated_numeric_status',num_int; write(io,'(A,1X,I0)') 'stage8_total_numeric_summary_status',num_summary
write(io,'(A,1X,I0)') 'stage8_total_force_conservation_status',cons_f; write(io,'(A,1X,I0)') 'stage8_total_component_conservation_status',cons_comp; write(io,'(A,1X,I0)') 'stage8_total_action_reaction_status',cons_ar; write(io,'(A,1X,I0)') 'stage8_total_pair_power_status',cons_pair; write(io,'(A,1X,I0)') 'stage8_total_fluid_power_consistency_status',cons_fp; write(io,'(A,1X,I0)') 'stage8_total_force_density_convention_status',cons_fd; write(io,'(A,1X,I0)') 'stage8_total_rho_convention_status',cons_rho; write(io,'(A,1X,I0)') 'stage8_total_no_double_rho_division_status',cons_nodbl; write(io,'(A,1X,I0)') 'stage8_total_conservation_power_rho_summary_status',cons_summary
write(io,'(A,1X,I0)') 'stage8_total_boundary_safe_workflow_status',bd_safe; write(io,'(A,1X,I0)') 'stage8_total_nearwall_blocked_status',bd_near; write(io,'(A,1X,I0)') 'stage8_total_outside_y_blocked_status',bd_outy; write(io,'(A,1X,I0)') 'stage8_total_invalid_coord_rejection_status',bd_invcoord; write(io,'(A,1X,I0)') 'stage8_total_invalid_layout_rejection_status',bd_invlayout; write(io,'(A,1X,I0)') 'stage8_total_periodic_boundary_status',bd_per; write(io,'(A,1X,I0)') 'stage8_total_mixed_safe_blocked_status',bd_mix; write(io,'(A,1X,I0)') 'stage8_total_boundary_summary_status',bd_summary
write(io,'(A,1X,I0)') 'stage8_total_rk_substep_count_status',rk_nsub; write(io,'(A,1X,I0)') 'stage8_total_rk_counter_status',rk_ctr; write(io,'(A,1X,I0)') 'stage8_total_rk_event_order_status',rk_ord; write(io,'(A,1X,I0)') 'stage8_total_rk_conservation_status',rk_cons; write(io,'(A,1X,I0)') 'stage8_total_rk_blocked_substep_status',rk_blk; write(io,'(A,1X,I0)') 'stage8_total_rk_summary_status',rk_summary
write(io,'(A,1X,I0)') 'stage8_total_rhs_untouched_status',rhs_untouched; write(io,'(A,1X,I0)') 'stage8_total_pressure_poisson_untouched_status',pp_untouched; write(io,'(A,1X,I0)') 'stage8_total_projection_untouched_status',proj_untouched; write(io,'(A,1X,I0)') 'stage8_total_production_dns_not_called_status',dns_not_called; write(io,'(A,1X,I0)') 'stage8_total_fluid_update_not_called_status',fluid_not_called; write(io,'(A,1X,I0)') 'stage8_total_fibre_advance_not_called_status',fibre_not_called; write(io,'(A,1X,I0)') 'stage8_total_no_production_side_effect_status',noside
write(io,'(A,1X,I0)') 'stage8_total_closed_marker_written_flag',marker_written
write(io,'(A,1X,I0)') 'stage8_total_closed_marker_status',marker_status
write(io,'(A,1X,I0)') 'stage8_total_smoke_check_status',final_status
close(io)
contains
subroutine ensure_dir(p); character(len=*),intent(in)::p; call execute_command_line('mkdir -p '//trim(p)); end
subroutine rm_file(path); character(len=*),intent(in)::path; logical::ex; inquire(file=path,exist=ex); if(ex) call execute_command_line('rm -f '//trim(path)); end
subroutine file_exists_int(path,flag); character(len=*),intent(in)::path; integer,intent(out)::flag; logical::ex; inquire(file=path,exist=ex); flag=merge(1,0,ex); end
subroutine file_nonempty_int(path,flag); character(len=*),intent(in)::path; integer,intent(out)::flag; logical::ex; integer(kind=8)::sz; inquire(file=path,exist=ex,size=sz); flag=merge(1,0,ex.and.sz>0); end
subroutine get_int(path,key,val)
character(len=*),intent(in)::path,key; integer,intent(out)::val
integer::u,ios; character(len=256)::k; real(mytype)::xx; logical::ex
val=0; inquire(file=path,exist=ex); if(.not.ex)return
open(newunit=u,file=path,status='old',action='read',iostat=ios); if(ios/=0)return
do; read(u,*,iostat=ios) k,xx; if(ios/=0)exit; if(trim(k)==trim(key)) then; val=nint(xx); exit; endif; enddo; close(u)
end
integer function all_zero5(f1,k1,f2,k2,f3,k3,f4,k4,f5,k5)
character(len=*),intent(in)::f1,k1,f2,k2,f3,k3,f4,k4,f5,k5
integer::a,b,c,d,e
call get_int(f1,k1,a); call get_int(f2,k2,b); call get_int(f3,k3,c); call get_int(f4,k4,d); call get_int(f5,k5,e)
all_zero5=merge(1,0,a==0.and.b==0.and.c==0.and.d==0.and.e==0)
end function
subroutine write_marker(path)
character(len=*),intent(in)::path; integer::u
open(newunit=u,file=path,status='replace',action='write')
write(u,'(A)') 'STAGE 8 CLOSED'
write(u,'(A)') 'all Stage 8 checks passed'
write(u,'(A)') 'no production RHS/projection/DNS touched'
write(u,'(A)') 'Stage 9 may begin only after this marker and status are both present.'
close(u)
end
end program
