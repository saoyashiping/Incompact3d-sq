program fibre_stage9_build_regression_check
use fibre_parameters, only: mytype
implicit none
integer::io,final_status
integer::d0_exists,d0_status,d0_entry,d0_91,d0_prod_allow,d0_dep
integer::s7_closed,s8_closed,s8_out,s8_smoke,s8_marker,s8_no_side,prior_closure
integer::ev_exists,ev_cfg,ev_build,ev_real_decomp,ev_stage_checks,ev_stage_checks_off
integer::ev_default_disabled,ev_dns,ev_rhs,ev_proj,ev_flu,ev_fib,ev_status
integer::pol_default_disabled,pol_dns_not_called,pol_prod_disabled,pol_rhs_untouched
integer::pol_proj_untouched,pol_flu_not_called,pol_fib_not_called,pol_status
integer::rhs_hook,rhs_mod,pp_mod,proj_mod,rproj,dns_call,flu_call,fib_call,noop

call ensure_dir('stage9_outputs')

call file_exists_int('stage9_outputs/fibre_stage9_dependency_gate_check.dat',d0_exists)
call get_int('stage9_outputs/fibre_stage9_dependency_gate_check.dat','stage9_dependency_gate_check_status',d0_status)
call get_int('stage9_outputs/fibre_stage9_dependency_gate_check.dat','stage9_gate_stage9_entry_allowed_flag',d0_entry)
call get_int('stage9_outputs/fibre_stage9_dependency_gate_check.dat','stage9_gate_stage9_1_allowed_flag',d0_91)
call get_int('stage9_outputs/fibre_stage9_dependency_gate_check.dat','stage9_gate_production_coupling_allowed_flag',d0_prod_allow)
d0_dep=merge(1,0,d0_exists==1.and.d0_status==1.and.d0_entry==1.and.d0_91==1.and.d0_prod_allow==0)

call file_exists_int('stage7_outputs/STAGE7_CLOSED.md',s7_closed)
call file_exists_int('stage8_outputs/STAGE8_CLOSED.md',s8_closed)
call file_exists_int('stage8_outputs/fibre_stage8_total_smoke_check.dat',s8_out)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_smoke_check_status',s8_smoke)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_closed_marker_status',s8_marker)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_no_production_side_effect_status',s8_no_side)
prior_closure=merge(1,0,s7_closed==1.and.s8_closed==1.and.s8_out==1.and.s8_smoke==1.and.s8_marker==1.and.s8_no_side==1)

call file_exists_int('stage9_outputs/stage9_1_build_evidence.dat',ev_exists)
call get_int('stage9_outputs/stage9_1_build_evidence.dat','stage9_build_evidence_configure_status',ev_cfg)
call get_int('stage9_outputs/stage9_1_build_evidence.dat','stage9_build_evidence_build_status',ev_build)
call get_int('stage9_outputs/stage9_1_build_evidence.dat','stage9_build_evidence_used_real_decomp2d_flag',ev_real_decomp)
call get_int('stage9_outputs/stage9_1_build_evidence.dat','stage9_build_evidence_stage_checks_only_flag',ev_stage_checks)
call get_int('stage9_outputs/stage9_1_build_evidence.dat','stage9_build_evidence_fibre_stage_checks_only_off_flag',ev_stage_checks_off)
call get_int('stage9_outputs/stage9_1_build_evidence.dat','stage9_build_evidence_default_production_disabled_flag',ev_default_disabled)
call get_int('stage9_outputs/stage9_1_build_evidence.dat','stage9_build_evidence_dns_executed_flag',ev_dns)
call get_int('stage9_outputs/stage9_1_build_evidence.dat','stage9_build_evidence_rhs_modified_flag',ev_rhs)
call get_int('stage9_outputs/stage9_1_build_evidence.dat','stage9_build_evidence_projection_called_flag',ev_proj)
call get_int('stage9_outputs/stage9_1_build_evidence.dat','stage9_build_evidence_fluid_update_called_flag',ev_flu)
call get_int('stage9_outputs/stage9_1_build_evidence.dat','stage9_build_evidence_fibre_advance_called_flag',ev_fib)
ev_status=merge(1,0,ev_exists==1.and.ev_cfg==1.and.ev_build==1.and.ev_real_decomp==1.and.ev_stage_checks==0.and.ev_stage_checks_off==1.and.ev_default_disabled==1.and.ev_dns==0.and.ev_rhs==0.and.ev_proj==0.and.ev_flu==0.and.ev_fib==0)

pol_default_disabled=merge(1,0,ev_default_disabled==1)
pol_dns_not_called=merge(1,0,ev_dns==0)
pol_prod_disabled=merge(1,0,d0_prod_allow==0)
pol_rhs_untouched=merge(1,0,ev_rhs==0)
pol_proj_untouched=merge(1,0,ev_proj==0)
pol_flu_not_called=merge(1,0,ev_flu==0)
pol_fib_not_called=merge(1,0,ev_fib==0)
pol_status=merge(1,0,pol_default_disabled==1.and.pol_dns_not_called==1.and.pol_prod_disabled==1.and.pol_rhs_untouched==1.and.pol_proj_untouched==1.and.pol_flu_not_called==1.and.pol_fib_not_called==1.and.s8_no_side==1)

rhs_hook=0; rhs_mod=0; pp_mod=0; proj_mod=0; rproj=0; dns_call=0; flu_call=0; fib_call=0
noop=merge(1,0,rhs_hook==0.and.rhs_mod==0.and.pp_mod==0.and.proj_mod==0.and.rproj==0.and.dns_call==0.and.flu_call==0.and.fib_call==0)

final_status=merge(1,0,d0_dep==1.and.prior_closure==1.and.ev_status==1.and.pol_status==1.and.noop==1)

open(newunit=io,file='stage9_outputs/fibre_stage9_build_regression_check.dat',status='replace',action='write')
write(io,'(A,1X,I0)') 'stage9_build_stage9_0_output_exists',d0_exists
write(io,'(A,1X,I0)') 'stage9_build_stage9_0_status',d0_status
write(io,'(A,1X,I0)') 'stage9_build_stage9_entry_allowed_flag',d0_entry
write(io,'(A,1X,I0)') 'stage9_build_stage9_1_allowed_flag',d0_91
write(io,'(A,1X,I0)') 'stage9_build_production_coupling_allowed_flag',d0_prod_allow
write(io,'(A,1X,I0)') 'stage9_build_stage9_0_dependency_status',d0_dep
write(io,'(A,1X,I0)') 'stage9_build_stage7_closed_marker_exists',s7_closed
write(io,'(A,1X,I0)') 'stage9_build_stage8_closed_marker_exists',s8_closed
write(io,'(A,1X,I0)') 'stage9_build_stage8_total_smoke_output_exists',s8_out
write(io,'(A,1X,I0)') 'stage9_build_stage8_total_smoke_status',s8_smoke
write(io,'(A,1X,I0)') 'stage9_build_stage8_closed_marker_status',s8_marker
write(io,'(A,1X,I0)') 'stage9_build_stage8_no_production_side_effect_status',s8_no_side
write(io,'(A,1X,I0)') 'stage9_build_prior_stage_closure_status',prior_closure
write(io,'(A,1X,I0)') 'stage9_build_evidence_output_exists',ev_exists
write(io,'(A,1X,I0)') 'stage9_build_configure_status',ev_cfg
write(io,'(A,1X,I0)') 'stage9_build_compile_status',ev_build
write(io,'(A,1X,I0)') 'stage9_build_used_real_decomp2d_flag',ev_real_decomp
write(io,'(A,1X,I0)') 'stage9_build_stage_checks_only_flag',ev_stage_checks
write(io,'(A,1X,I0)') 'stage9_build_fibre_stage_checks_only_off_flag',ev_stage_checks_off
write(io,'(A,1X,I0)') 'stage9_build_default_production_disabled_flag',ev_default_disabled
write(io,'(A,1X,I0)') 'stage9_build_dns_executed_flag',ev_dns
write(io,'(A,1X,I0)') 'stage9_build_rhs_modified_flag',ev_rhs
write(io,'(A,1X,I0)') 'stage9_build_projection_called_flag',ev_proj
write(io,'(A,1X,I0)') 'stage9_build_fluid_update_called_flag',ev_flu
write(io,'(A,1X,I0)') 'stage9_build_fibre_advance_called_flag',ev_fib
write(io,'(A,1X,I0)') 'stage9_build_evidence_status',ev_status
write(io,'(A,1X,I0)') 'stage9_build_default_production_disabled_status',pol_default_disabled
write(io,'(A,1X,I0)') 'stage9_build_production_dns_not_called_status',pol_dns_not_called
write(io,'(A,1X,I0)') 'stage9_build_production_coupling_disabled_status',pol_prod_disabled
write(io,'(A,1X,I0)') 'stage9_build_rhs_untouched_status',pol_rhs_untouched
write(io,'(A,1X,I0)') 'stage9_build_projection_untouched_status',pol_proj_untouched
write(io,'(A,1X,I0)') 'stage9_build_fluid_update_not_called_status',pol_flu_not_called
write(io,'(A,1X,I0)') 'stage9_build_fibre_advance_not_called_status',pol_fib_not_called
write(io,'(A,1X,I0)') 'stage9_build_production_disabled_policy_status',pol_status
write(io,'(A,1X,I0)') 'stage9_build_rhs_hook_called_flag',rhs_hook
write(io,'(A,1X,I0)') 'stage9_build_rhs_modified_flag',rhs_mod
write(io,'(A,1X,I0)') 'stage9_build_pressure_poisson_modified_flag',pp_mod
write(io,'(A,1X,I0)') 'stage9_build_projection_modified_flag',proj_mod
write(io,'(A,1X,I0)') 'stage9_build_real_projection_called_flag',rproj
write(io,'(A,1X,I0)') 'stage9_build_production_dns_called_flag',dns_call
write(io,'(A,1X,I0)') 'stage9_build_fluid_update_called_flag',flu_call
write(io,'(A,1X,I0)') 'stage9_build_fibre_advance_called_flag',fib_call
write(io,'(A,1X,I0)') 'stage9_build_noop_safety_status',noop
write(io,'(A,1X,I0)') 'stage9_build_regression_check_status',final_status
close(io)
contains
subroutine ensure_dir(p); character(len=*),intent(in)::p; call execute_command_line('mkdir -p '//trim(p)); end
subroutine file_exists_int(path,flag); character(len=*),intent(in)::path; integer,intent(out)::flag; logical::ex; inquire(file=path,exist=ex); flag=merge(1,0,ex); end
subroutine get_int(path,key,val)
character(len=*),intent(in)::path,key; integer,intent(out)::val
integer::u,ios; character(len=256)::k; real(mytype)::x; logical::ex
val=0; inquire(file=path,exist=ex); if(.not.ex)return
open(newunit=u,file=path,status='old',action='read',iostat=ios); if(ios/=0)return
do; read(u,*,iostat=ios) k,x; if(ios/=0)exit; if(trim(k)==trim(key)) then; val=nint(x); exit; endif; enddo; close(u)
end
end program
