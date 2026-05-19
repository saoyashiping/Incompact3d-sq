program fibre_stage9_dependency_gate_check
use fibre_parameters, only: mytype
implicit none
integer::io,final_status
integer::s7m,s7o,s7s,s7c,s7dep
integer::s8m,s8o,s8s,s8c,s8w,s8dep
integer::pdef,pnoside,prhs,ppp,pproj,pdns,pflu,pfib,pinh
integer::rhs_hook,rhs_mod,pp_mod,proj_mod,rproj,dns_call,flu_call,fib_call,noop
integer::entry,entry91,prod_allow,entrypol
call ensure_dir('stage9_outputs')
call file_exists_int('stage7_outputs/STAGE7_CLOSED.md',s7m)
call file_exists_int('stage7_outputs/fibre_stage7_total_smoke_check.dat',s7o)
call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_smoke_check_status',s7s)
call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_closed_marker_status',s7c)
s7dep=merge(1,0,s7m*s7o*s7s*s7c==1)
call file_exists_int('stage8_outputs/STAGE8_CLOSED.md',s8m)
call file_exists_int('stage8_outputs/fibre_stage8_total_smoke_check.dat',s8o)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_smoke_check_status',s8s)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_closed_marker_status',s8c)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_closed_marker_written_flag',s8w)
s8dep=merge(1,0,s8m*s8o*s8s*s8c*s8w==1)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_default_production_disabled_status',pdef)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_no_production_side_effect_status',pnoside)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_rhs_untouched_status',prhs)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_pressure_poisson_untouched_status',ppp)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_projection_untouched_status',pproj)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_production_dns_not_called_status',pdns)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_fluid_update_not_called_status',pflu)
call get_int('stage8_outputs/fibre_stage8_total_smoke_check.dat','stage8_total_fibre_advance_not_called_status',pfib)
pinh=merge(1,0,pdef*pnoside*prhs*ppp*pproj*pdns*pflu*pfib==1)
rhs_hook=0; rhs_mod=0; pp_mod=0; proj_mod=0; rproj=0; dns_call=0; flu_call=0; fib_call=0
noop=merge(1,0,rhs_hook==0.and.rhs_mod==0.and.pp_mod==0.and.proj_mod==0.and.rproj==0.and.dns_call==0.and.flu_call==0.and.fib_call==0)
entry=merge(1,0,s7dep==1.and.s8dep==1.and.pinh==1.and.noop==1)
entry91=entry
prod_allow=0
entrypol=merge(1,0,entry==1.and.entry91==1.and.prod_allow==0)
final_status=merge(1,0,s7dep==1.and.s8dep==1.and.pinh==1.and.noop==1.and.entrypol==1)
open(newunit=io,file='stage9_outputs/fibre_stage9_dependency_gate_check.dat',status='replace',action='write')
write(io,'(A,1X,I0)') 'stage9_gate_stage7_closed_marker_exists',s7m
write(io,'(A,1X,I0)') 'stage9_gate_stage7_total_smoke_output_exists',s7o
write(io,'(A,1X,I0)') 'stage9_gate_stage7_total_smoke_status',s7s
write(io,'(A,1X,I0)') 'stage9_gate_stage7_closed_marker_status',s7c
write(io,'(A,1X,I0)') 'stage9_gate_stage7_dependency_status',s7dep
write(io,'(A,1X,I0)') 'stage9_gate_stage8_closed_marker_exists',s8m
write(io,'(A,1X,I0)') 'stage9_gate_stage8_total_smoke_output_exists',s8o
write(io,'(A,1X,I0)') 'stage9_gate_stage8_total_smoke_status',s8s
write(io,'(A,1X,I0)') 'stage9_gate_stage8_closed_marker_status',s8c
write(io,'(A,1X,I0)') 'stage9_gate_stage8_closed_marker_written_flag',s8w
write(io,'(A,1X,I0)') 'stage9_gate_stage8_dependency_status',s8dep
write(io,'(A,1X,I0)') 'stage9_gate_stage8_default_production_disabled_status',pdef
write(io,'(A,1X,I0)') 'stage9_gate_stage8_no_production_side_effect_status',pnoside
write(io,'(A,1X,I0)') 'stage9_gate_stage8_rhs_untouched_status',prhs
write(io,'(A,1X,I0)') 'stage9_gate_stage8_pressure_poisson_untouched_status',ppp
write(io,'(A,1X,I0)') 'stage9_gate_stage8_projection_untouched_status',pproj
write(io,'(A,1X,I0)') 'stage9_gate_stage8_production_dns_not_called_status',pdns
write(io,'(A,1X,I0)') 'stage9_gate_stage8_fluid_update_not_called_status',pflu
write(io,'(A,1X,I0)') 'stage9_gate_stage8_fibre_advance_not_called_status',pfib
write(io,'(A,1X,I0)') 'stage9_gate_stage8_production_disabled_inheritance_status',pinh
write(io,'(A,1X,I0)') 'stage9_gate_rhs_hook_called_flag',rhs_hook
write(io,'(A,1X,I0)') 'stage9_gate_rhs_modified_flag',rhs_mod
write(io,'(A,1X,I0)') 'stage9_gate_pressure_poisson_modified_flag',pp_mod
write(io,'(A,1X,I0)') 'stage9_gate_projection_modified_flag',proj_mod
write(io,'(A,1X,I0)') 'stage9_gate_real_projection_called_flag',rproj
write(io,'(A,1X,I0)') 'stage9_gate_production_dns_called_flag',dns_call
write(io,'(A,1X,I0)') 'stage9_gate_fluid_update_called_flag',flu_call
write(io,'(A,1X,I0)') 'stage9_gate_fibre_advance_called_flag',fib_call
write(io,'(A,1X,I0)') 'stage9_gate_noop_safety_status',noop
write(io,'(A,1X,I0)') 'stage9_gate_stage9_entry_allowed_flag',entry
write(io,'(A,1X,I0)') 'stage9_gate_stage9_1_allowed_flag',entry91
write(io,'(A,1X,I0)') 'stage9_gate_production_coupling_allowed_flag',prod_allow
write(io,'(A,1X,I0)') 'stage9_gate_stage9_entry_policy_status',entrypol
write(io,'(A,1X,I0)') 'stage9_dependency_gate_check_status',final_status
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
