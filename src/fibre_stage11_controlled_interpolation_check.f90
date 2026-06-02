program fibre_stage11_controlled_interpolation_check
use fibre_stage11_config, only: stage11_config_load,stage11_requested,stage11_readonly_mode,stage11_get_max_points
use fibre_stage11_lagrangian_state, only: stage11_lagrangian_state_init,stage11_lagrangian_state_finalize,stage11_lagrangian_state_is_allocated
use fibre_stage11_grid_metadata, only: stage11_grid_metadata_init_uniform,stage11_grid_metadata_finalize,stage11_grid_metadata_is_initialized
use fibre_stage11_controlled_interpolation, only: stage11_controlled_interpolation_init,stage11_controlled_interpolation_finalize,stage11_controlled_interpolation_get_status_values,stage11_controlled_interpolation_write_diagnostics
implicit none
integer::rf,ro,lsa,gma,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,np,io,ios
logical::pass
call execute_command_line('mkdir -p stage11_outputs', exitstat=ios)
call stage11_config_load(); rf=merge(1,0,stage11_requested()); ro=merge(1,0,stage11_readonly_mode())
np=stage11_get_max_points(); if(np<=0)np=8
call stage11_lagrangian_state_init(np); lsa=merge(1,0,stage11_lagrangian_state_is_allocated())
call stage11_grid_metadata_init_uniform(); gma=merge(1,0,stage11_grid_metadata_is_initialized())
call stage11_controlled_interpolation_init()
call stage11_controlled_interpolation_get_status_values(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)
open(newunit=io,file='stage11_outputs/fibre_stage11_4_controlled_interpolation.dat',status='replace',action='write',iostat=ios)
if(ios/=0) then; print *,'STAGE 11.4 CONTROLLED INTERPOLATION VERDICT: FAIL'; stop 1; end if
write(io,'(A,1X,I0)') 'stage11_4_requested_flag',rf
write(io,'(A,1X,I0)') 'stage11_4_readonly_mode_status',ro
write(io,'(A,1X,I0)') 'stage11_4_lagrangian_state_available_status',lsa
write(io,'(A,1X,I0)') 'stage11_4_grid_metadata_available_status',gma
call stage11_controlled_interpolation_write_diagnostics(io)
close(io)
pass=.true.
if(rf/=1.or.ro/=1.or.lsa/=1.or.gma/=1) pass=.false.
if(a/=1.or.b/=1.or.c/=1.or.d/=1.or.e/=1.or.f/=1.or.g/=1.or.h/=1.or.i/=1.or.j/=1.or.k/=1.or.l/=1.or.m/=1.or.n/=1.or.o/=1.or.p/=1.or.q/=1) pass=.false.
if(pass) then; print *,'STAGE 11.4 CONTROLLED INTERPOLATION VERDICT: PASS'; else; print *,'STAGE 11.4 CONTROLLED INTERPOLATION VERDICT: FAIL'; stop 1; end if
call stage11_controlled_interpolation_finalize(); call stage11_grid_metadata_finalize(); call stage11_lagrangian_state_finalize()
end program
