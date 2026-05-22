module fibre_stage9_7_stats_visu_io_smoke
  use mpi
  use, intrinsic :: ieee_arithmetic
  use decomp_2d_mpi, only : nrank
  use decomp_2d_constants, only : mytype
  use param, only : itype, itype_channel
  implicit none
  private
  logical, save :: enabled=.false.
  integer, save :: requested_max_steps=0, completed_steps=0, finalise_reached=0
  integer, save :: require_stats=1, require_visu=1, require_coarse_io=1
  integer, save :: stats_path=0, stats_output=0, stats_finite=1
  integer, save :: visu_path=0, visu_output=0, visu_finite=1
  integer, save :: coarse_desc=0, io_open=0, io_write=0, io_close=0
  integer, save :: expected_stats_files=0, expected_visu_files=0, output_nonempty=0, output_field_finite=1
  public :: stage9_7_requested, stage9_7_get_max_steps, stage9_7_get_requirements, stage9_7_begin
  public :: stage9_7_record_stats_path, stage9_7_record_visu_path, stage9_7_record_coarse_io_path
  public :: stage9_7_record_output_file_status, stage9_7_record_field_finite_status, stage9_7_after_completed_step
  public :: stage9_7_finalise_mark, stage9_7_should_stop, stage9_7_progress_note, stage9_7_final_audit
contains
logical function stage9_7_requested() result(v)
character(len=64)::e;integer::st; v=.false.;e=''; call get_environment_variable('X3D_STAGE9_7_STATS_VISU_IO_SMOKE', value=e, status=st); if(st==0) v=(trim(adjustl(e))=='1')
end function
integer function stage9_7_get_max_steps(d) result(v)
integer,intent(in)::d; character(len=64)::e; integer::st,ios,t; v=d; e=''; call get_environment_variable('X3D_STAGE9_7_MAX_STEPS', value=e, status=st); if(st==0) then; read(e,*,iostat=ios) t; if(ios==0.and.t>0) v=t; endif
end function
subroutine stage9_7_get_requirements(rstats,rvisu,rcoarse)
integer,intent(out)::rstats,rvisu,rcoarse; character(len=64)::e; integer::st,ios,t
rstats=1; rvisu=1; rcoarse=1
e=''; call get_environment_variable('X3D_STAGE9_7_REQUIRE_STATS', value=e, status=st); if(st==0) then; read(e,*,iostat=ios) t; if(ios==0) rstats=merge(1,0,t/=0); endif
e=''; call get_environment_variable('X3D_STAGE9_7_REQUIRE_VISU', value=e, status=st); if(st==0) then; read(e,*,iostat=ios) t; if(ios==0) rvisu=merge(1,0,t/=0); endif
e=''; call get_environment_variable('X3D_STAGE9_7_REQUIRE_COARSE_IO', value=e, status=st); if(st==0) then; read(e,*,iostat=ios) t; if(ios==0) rcoarse=merge(1,0,t/=0); endif
end subroutine
subroutine stage9_7_begin(en,mx,rs,rv,rc)
logical,intent(in)::en; integer,intent(in)::mx,rs,rv,rc
enabled=en; requested_max_steps=mx; require_stats=rs; require_visu=rv; require_coarse_io=rc; completed_steps=0; finalise_reached=0
stats_path=0;stats_output=0;stats_finite=1;visu_path=0;visu_output=0;visu_finite=1;coarse_desc=0;io_open=0;io_write=0;io_close=0
expected_stats_files=0;expected_visu_files=0;output_nonempty=0;output_field_finite=1
end subroutine
subroutine stage9_7_record_stats_path(output_called,finite_ok)
integer,intent(in)::output_called,finite_ok; if(.not.enabled) return; stats_path=1; stats_output=max(stats_output,output_called); stats_finite=min(stats_finite,finite_ok)
end subroutine
subroutine stage9_7_record_visu_path(output_called,finite_ok)
integer,intent(in)::output_called,finite_ok; if(.not.enabled) return; visu_path=1; visu_output=max(visu_output,output_called); visu_finite=min(visu_finite,finite_ok)
end subroutine
subroutine stage9_7_record_coarse_io_path(desc_ok,open_ok,write_ok,close_ok)
integer,intent(in)::desc_ok,open_ok,write_ok,close_ok; if(.not.enabled) return; coarse_desc=max(coarse_desc,desc_ok); io_open=max(io_open,open_ok); io_write=max(io_write,write_ok); io_close=max(io_close,close_ok)
end subroutine
subroutine stage9_7_record_output_file_status(sf,vf,ne)
integer,intent(in)::sf,vf,ne; if(.not.enabled) return; expected_stats_files=max(expected_stats_files,sf); expected_visu_files=max(expected_visu_files,vf); output_nonempty=max(output_nonempty,ne); if (vf==1 .and. ne==1) visu_output=max(visu_output,1); if (sf==1 .and. ne==1) stats_output=max(stats_output,1)
end subroutine
subroutine stage9_7_record_field_finite_status(ux,uy,uz,pp,divu)
real(mytype),intent(in)::ux(:,:,:),uy(:,:,:),uz(:,:,:),pp(:,:,:),divu(:,:,:)
integer::loc,glob,ierr; if(.not.enabled) return
loc=merge(1,0,ieee_is_finite(maxval(abs(ux))).and.ieee_is_finite(maxval(abs(uy))).and.ieee_is_finite(maxval(abs(uz))).and.ieee_is_finite(maxval(abs(pp))).and.ieee_is_finite(maxval(abs(divu))))
call MPI_Allreduce(loc,glob,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr); output_field_finite=min(output_field_finite,glob)
end subroutine
subroutine stage9_7_after_completed_step(); if(enabled) completed_steps=completed_steps+1; end
subroutine stage9_7_finalise_mark(); finalise_reached=1; end
subroutine stage9_7_progress_note()
  if (enabled .and. nrank==0) then
    write(*,'(A,I0,A,I0)') '[STAGE9.7] completed outer step ',completed_steps,' / ',requested_max_steps
  endif
end subroutine
logical function stage9_7_should_stop() result(v)
  v = enabled .and. (completed_steps>=requested_max_steps)
end function
subroutine stage9_7_final_audit()
integer::u,s_case,s_noc,s_time,s_stats,s_visu,s_coarse,s_nonan,s_final,i
s_case=merge(1,0,itype==itype_channel); s_noc=1; s_time=merge(1,0,completed_steps>=1.and.completed_steps<=requested_max_steps)
s_stats=merge(1,0,(require_stats==0).or.(stats_path==1.and.stats_output==1.and.stats_finite==1.and.expected_stats_files==1))
s_visu=merge(1,0,(require_visu==0).or.(visu_path==1.and.visu_finite==1.and.expected_visu_files==1.and.output_nonempty==1.and.io_write==1))
s_coarse=merge(1,0,(require_coarse_io==0).or.(coarse_desc==1.and.io_open==1.and.io_write==1.and.io_close==1))
s_nonan=merge(1,0,stats_finite==1.and.visu_finite==1.and.output_field_finite==1.and.output_nonempty==1)
s_final=min(s_case,min(s_noc,min(s_time,min(s_stats,min(s_visu,min(s_coarse,min(output_field_finite,min(s_nonan,finalise_reached))))))))
call MPI_Allreduce(MPI_IN_PLACE,s_final,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,i)
if(nrank==0) then
open(newunit=u,file='stage9_outputs/fibre_stage9_7_stats_visu_io_smoke.dat',status='replace',action='write')
write(u,*) 'stage9_7_requested_flag ',merge(1,0,enabled);write(u,*) 'stage9_7_channel_case_status ',s_case;write(u,*) 'stage9_7_no_fibre_coupling_status ',s_noc
write(u,*) 'stage9_7_requested_max_steps ',requested_max_steps;write(u,*) 'stage9_7_completed_steps ',completed_steps;write(u,*) 'stage9_7_time_advance_status ',s_time
write(u,*) 'stage9_7_stats_path_executed_status ',stats_path;write(u,*) 'stage9_7_stats_output_status ',stats_output;write(u,*) 'stage9_7_stats_finite_status ',stats_finite
write(u,*) 'stage9_7_visu_path_executed_status ',visu_path;write(u,*) 'stage9_7_visu_output_status ',visu_output;write(u,*) 'stage9_7_visu_field_finite_status ',visu_finite
write(u,*) 'stage9_7_coarse_mesh_descriptor_status ',coarse_desc;write(u,*) 'stage9_7_decomp_io_open_status ',io_open;write(u,*) 'stage9_7_decomp_io_write_status ',io_write;write(u,*) 'stage9_7_decomp_io_close_status ',io_close
write(u,*) 'stage9_7_expected_stats_files_status ',expected_stats_files;write(u,*) 'stage9_7_expected_visu_files_status ',expected_visu_files;write(u,*) 'stage9_7_output_files_nonempty_status ',output_nonempty
write(u,*) 'stage9_7_output_field_finite_status ',output_field_finite;write(u,*) 'stage9_7_no_nan_inf_status ',s_nonan;write(u,*) 'stage9_7_finalise_reached_status ',finalise_reached
write(u,*) 'stage9_7_stats_visu_io_smoke_status ',s_final
close(u)
if(s_final==1) then; write(*,'(A)') 'STAGE 9.7 STATS VISU IO SMOKE VERDICT: PASS'; else; write(*,'(A)') 'STAGE 9.7 STATS VISU IO SMOKE VERDICT: FAIL'; endif
endif
end subroutine
end module
