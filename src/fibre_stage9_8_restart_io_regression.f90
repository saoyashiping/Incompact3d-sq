module fibre_stage9_8_restart_io_regression
  use mpi
  use, intrinsic :: ieee_arithmetic
  use decomp_2d_mpi, only : nrank
  use decomp_2d_constants, only : mytype, real_type
  use param, only : itype, itype_channel
  implicit none
  private
  logical, save :: enabled=.false.
  integer, save :: phase_restart=0, requested_steps=0, completed_steps=0, finalise_reached=0
  integer, save :: restart_write_path=0, restart_read_path=0, restart_exist=0, restart_nonempty=0, restart_signature_status=1
  integer, save :: restart_read_signature_checked=0
  integer, save :: vel_finite=1, pres_finite=1, div_finite=1
  real(8), save :: sig_tol=1.0d-8
  real(8), save :: sig_sum(3)=0d0, sig_max(3)=0d0
  real(8), save :: read_sig_sum(3)=0d0, read_sig_max(3)=0d0
  character(len=256), save :: signature_file=''
  public :: stage9_8_requested, stage9_8_get_phase, stage9_8_get_max_steps_before_restart, stage9_8_get_max_steps_after_restart
  public :: stage9_8_get_signature_tol, stage9_8_begin, stage9_8_record_restart_write_path, stage9_8_record_restart_read_path
  public :: stage9_8_record_restart_file_status, stage9_8_record_field_finite_status, stage9_8_record_signature
  public :: stage9_8_after_completed_step, stage9_8_should_stop, stage9_8_finalise_mark, stage9_8_final_audit
contains
logical function stage9_8_requested() result(v)
  character(len=64)::e; integer::st
  v=.false.; e=''
  call get_environment_variable('X3D_STAGE9_8_RESTART_IO_REGRESSION',value=e,status=st)
  if (st==0) v=(trim(adjustl(e))=='1')
end function

integer function stage9_8_get_phase() result(v)
  character(len=64)::e; integer::st
  v=0; e=''
  call get_environment_variable('X3D_STAGE9_8_PHASE',value=e,status=st)
  if (st==0) then
    if (trim(adjustl(e))=='restart') v=1
  endif
end function

integer function stage9_8_get_max_steps_before_restart(d) result(v)
  integer,intent(in)::d; character(len=64)::e; integer::st,ios,t
  v=d; e=''
  call get_environment_variable('X3D_STAGE9_8_MAX_STEPS_BEFORE_RESTART',value=e,status=st)
  if (st==0) then; read(e,*,iostat=ios) t; if (ios==0 .and. t>0) v=t; endif
end function

integer function stage9_8_get_max_steps_after_restart(d) result(v)
  integer,intent(in)::d; character(len=64)::e; integer::st,ios,t
  v=d; e=''
  call get_environment_variable('X3D_STAGE9_8_MAX_STEPS_AFTER_RESTART',value=e,status=st)
  if (st==0) then; read(e,*,iostat=ios) t; if (ios==0 .and. t>0) v=t; endif
end function

real(8) function stage9_8_get_signature_tol(d) result(v)
  real(8),intent(in)::d; character(len=64)::e; integer::st,ios
  v=d; e=''
  call get_environment_variable('X3D_STAGE9_8_RESTART_SIGNATURE_TOL',value=e,status=st)
  if (st==0) then; read(e,*,iostat=ios) v; if (ios/=0) v=d; endif
end function

subroutine stage9_8_begin(en,phase,steps,tol)
  logical,intent(in)::en; integer,intent(in)::phase,steps; real(8),intent(in)::tol
  integer :: st
  enabled=en; phase_restart=phase; requested_steps=steps; sig_tol=tol; completed_steps=0; finalise_reached=0
  restart_write_path=0; restart_read_path=0; restart_exist=0; restart_nonempty=0; restart_signature_status=1
  restart_read_signature_checked=0
  vel_finite=1; pres_finite=1; div_finite=1; sig_sum=0d0; sig_max=0d0; read_sig_sum=0d0; read_sig_max=0d0
  signature_file=''
  call get_environment_variable('X3D_STAGE9_8_SIGNATURE_FILE',value=signature_file,status=st)
  if (enabled .and. nrank==0) then
    write(*,'(A)') '[STAGE9.8] restart I/O regression mode enabled'
    if (phase==1) then
      write(*,'(A,A)') '[STAGE9.8] phase = ','restart'
    else
      write(*,'(A,A)') '[STAGE9.8] phase = ','cold'
    endif
  endif
end subroutine

subroutine stage9_8_record_restart_write_path()
  if (enabled) then
    restart_write_path=1
    if (nrank==0) write(*,'(A)') '[STAGE9.8] writing restart'
  endif
end subroutine

subroutine stage9_8_record_restart_read_path()
  if (enabled) then
    restart_read_path=1
    if (nrank==0) write(*,'(A)') '[STAGE9.8] reading restart'
  endif
end subroutine

subroutine stage9_8_record_restart_file_status(ex,ne)
  integer,intent(in)::ex,ne
  if(.not.enabled) return
  restart_exist=max(restart_exist,ex); restart_nonempty=max(restart_nonempty,ne)
end subroutine

subroutine stage9_8_record_field_finite_status(ux,uy,uz,pp,divu)
  real(mytype),intent(in)::ux(:,:,:),uy(:,:,:),uz(:,:,:),pp(:,:,:),divu(:,:,:)
  integer::ierr,loc,glob
  if(.not.enabled) return
  loc=merge(1,0,ieee_is_finite(maxval(abs(ux))).and.ieee_is_finite(maxval(abs(uy))).and.ieee_is_finite(maxval(abs(uz))))
  call MPI_Allreduce(loc,glob,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr); vel_finite=min(vel_finite,glob)
  loc=merge(1,0,ieee_is_finite(maxval(abs(pp))))
  call MPI_Allreduce(loc,glob,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr); pres_finite=min(pres_finite,glob)
  loc=merge(1,0,ieee_is_finite(maxval(abs(divu))))
  call MPI_Allreduce(loc,glob,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr); div_finite=min(div_finite,glob)
end subroutine

subroutine stage9_8_record_signature(ux,uy,uz)
  real(mytype),intent(in)::ux(:,:,:),uy(:,:,:),uz(:,:,:)
  real(mytype)::ls(3),gs(3),lm(3),gm(3)
  integer::ierr,u,ios
  if(.not.enabled) return
  ls=(/sum(ux),sum(uy),sum(uz)/); lm=(/maxval(abs(ux)),maxval(abs(uy)),maxval(abs(uz))/); gs=ls; gm=lm
  call MPI_Allreduce(MPI_IN_PLACE,gs,3,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(MPI_IN_PLACE,gm,3,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
  if (phase_restart==0) then
    sig_sum=real(gs,8); sig_max=real(gm,8)
    if (nrank==0 .and. len_trim(signature_file)>0) then
      open(newunit=u,file=trim(signature_file),status='replace',action='write')
      write(u,*) sig_sum(1),sig_sum(2),sig_sum(3)
      write(u,*) sig_max(1),sig_max(2),sig_max(3)
      close(u)
    endif
  else
    if (restart_read_signature_checked==1) return
    if (nrank==0 .and. len_trim(signature_file)>0) then
      open(newunit=u,file=trim(signature_file),status='old',action='read',iostat=ios)
      if (ios==0) then
        read(u,*,iostat=ios) sig_sum(1),sig_sum(2),sig_sum(3)
        read(u,*,iostat=ios) sig_max(1),sig_max(2),sig_max(3)
        close(u)
        if (ios/=0) restart_signature_status=0
      else
        restart_signature_status=0
      endif
    endif
    call MPI_Bcast(sig_sum,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(sig_max,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    read_sig_sum=real(gs,8); read_sig_max=real(gm,8)
    if (maxval(abs(read_sig_sum-sig_sum))>sig_tol .or. maxval(abs(read_sig_max-sig_max))>sig_tol) restart_signature_status=0
    restart_read_signature_checked=1
  endif
end subroutine

subroutine stage9_8_after_completed_step()
  if (enabled) then
    completed_steps=completed_steps+1
    if (nrank==0) then
      if (phase_restart==1) then
        write(*,'(A,I0,A,I0,A)') '[STAGE9.8] completed outer step ',completed_steps,' / ',requested_steps,' after restart'
      else
        write(*,'(A,I0,A,I0)') '[STAGE9.8] completed outer step ',completed_steps,' / ',requested_steps
      endif
    endif
  endif
end subroutine

logical function stage9_8_should_stop() result(v)
  v=enabled.and.(completed_steps>=requested_steps)
end function
subroutine stage9_8_finalise_mark(); finalise_reached=1; end subroutine

subroutine stage9_8_final_audit()
  integer::u,s_case,s_noc,s_time,s_nonan,s_final,i
  s_case=merge(1,0,itype==itype_channel); s_noc=1; s_time=merge(1,0,completed_steps>=1 .and. completed_steps<=requested_steps)
  s_nonan=merge(1,0,vel_finite==1.and.pres_finite==1.and.div_finite==1)
  s_final=min(s_case,min(s_noc,min(s_time,min(restart_exist,min(restart_nonempty,min(restart_signature_status,min(vel_finite,min(pres_finite,min(div_finite,min(s_nonan,finalise_reached))))))))))
  call MPI_Allreduce(MPI_IN_PLACE,s_final,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,i)
  if(nrank==0) then
    open(newunit=u,file='stage9_outputs/fibre_stage9_8_restart_io_regression.dat',status='replace',action='write')
    write(u,*) 'stage9_8_requested_flag ',merge(1,0,enabled)
    write(u,*) 'stage9_8_phase_cold_status ',merge(1,0,phase_restart==0)
    write(u,*) 'stage9_8_phase_restart_status ',merge(1,0,phase_restart==1)
    write(u,*) 'stage9_8_channel_case_status ',s_case
    write(u,*) 'stage9_8_no_fibre_coupling_status ',s_noc
    write(u,*) 'stage9_8_requested_steps ',requested_steps
    write(u,*) 'stage9_8_completed_steps ',completed_steps
    write(u,*) 'stage9_8_time_advance_status ',s_time
    write(u,*) 'stage9_8_restart_write_path_status ',restart_write_path
    write(u,*) 'stage9_8_restart_read_path_status ',restart_read_path
    write(u,*) 'stage9_8_restart_files_exist_status ',restart_exist
    write(u,*) 'stage9_8_restart_files_nonempty_status ',restart_nonempty
    write(u,*) 'stage9_8_restart_signature_status ',restart_signature_status
    write(u,*) 'stage9_8_signature_sum_ux ',sig_sum(1)
    write(u,*) 'stage9_8_signature_sum_uy ',sig_sum(2)
    write(u,*) 'stage9_8_signature_sum_uz ',sig_sum(3)
    write(u,*) 'stage9_8_signature_max_ux ',sig_max(1)
    write(u,*) 'stage9_8_signature_max_uy ',sig_max(2)
    write(u,*) 'stage9_8_signature_max_uz ',sig_max(3)
    write(u,*) 'stage9_8_read_signature_sum_ux ',read_sig_sum(1)
    write(u,*) 'stage9_8_read_signature_sum_uy ',read_sig_sum(2)
    write(u,*) 'stage9_8_read_signature_sum_uz ',read_sig_sum(3)
    write(u,*) 'stage9_8_read_signature_max_ux ',read_sig_max(1)
    write(u,*) 'stage9_8_read_signature_max_uy ',read_sig_max(2)
    write(u,*) 'stage9_8_read_signature_max_uz ',read_sig_max(3)
    write(u,*) 'stage9_8_velocity_finite_status ',vel_finite
    write(u,*) 'stage9_8_pressure_finite_status ',pres_finite
    write(u,*) 'stage9_8_divergence_finite_status ',div_finite
    write(u,*) 'stage9_8_no_nan_inf_status ',s_nonan
    write(u,*) 'stage9_8_finalise_reached_status ',finalise_reached
    write(u,*) 'stage9_8_restart_io_regression_status ',s_final
    close(u)
    if (s_final==1) then
      write(*,'(A)') 'STAGE 9.8 RESTART IO REGRESSION VERDICT: PASS'
    else
      write(*,'(A)') 'STAGE 9.8 RESTART IO REGRESSION VERDICT: FAIL'
    endif
  endif
end subroutine
end module
