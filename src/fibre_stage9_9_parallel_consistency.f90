module fibre_stage9_9_parallel_consistency
  use mpi
  use, intrinsic :: ieee_arithmetic
  use decomp_2d_mpi, only : nrank
  use decomp_2d_constants, only : mytype, real_type
  use param, only : itype, itype_channel, dt, dx, dy, dz, istret
  implicit none
  private

  logical, save :: enabled=.false.
  integer, save :: requested_max_steps=0, completed_steps=0, finalise_reached=0
  real(8), save :: signature_tol=1.0d-8, divergence_tol=1.0d-8, massflux_tol=1.0d-6
  integer, save :: velocity_finite_status=1, pressure_finite_status=1, divergence_finite_status=1
  integer, save :: cfl_finite_status=1, massflux_finite_status=1, no_fibre_coupling_status=1
  integer, save :: decomp_invariant_init_status=0
  real(8), save :: sig_sum_ux=0d0, sig_sum_uy=0d0, sig_sum_uz=0d0
  real(8), save :: sig_max_ux=0d0, sig_max_uy=0d0, sig_max_uz=0d0
  real(8), save :: sig_l2_ux=0d0, sig_l2_uy=0d0, sig_l2_uz=0d0

  public :: stage9_9_requested, stage9_9_get_max_steps, stage9_9_get_tolerances
  public :: stage9_9_begin, stage9_9_record_field_signature, stage9_9_record_divergence_status
  public :: stage9_9_record_massflux_status, stage9_9_record_cfl_status, stage9_9_after_completed_step
  public :: stage9_9_should_stop, stage9_9_finalise_mark, stage9_9_final_audit
contains
  logical function stage9_9_requested() result(requested)
    character(len=64) :: e
    integer :: st
    requested=.false.; e=''
    call get_environment_variable('X3D_STAGE9_9_PARALLEL_CONSISTENCY',value=e,status=st)
    if (st==0) requested=(trim(adjustl(e))=='1')
  end function

  integer function stage9_9_get_max_steps(default_steps) result(max_steps)
    integer, intent(in) :: default_steps
    character(len=64) :: e
    integer :: st, ios, v
    max_steps=default_steps; e=''
    call get_environment_variable('X3D_STAGE9_9_MAX_STEPS',value=e,status=st)
    if (st==0) then
      read(e,*,iostat=ios) v
      if (ios==0 .and. v>0) max_steps=v
    endif
  end function

  subroutine stage9_9_get_tolerances(sig_tol,div_tol,mf_tol)
    real(8), intent(out) :: sig_tol,div_tol,mf_tol
    character(len=64) :: e
    integer :: st, ios
    sig_tol=1.0d-8; div_tol=1.0d-8; mf_tol=1.0d-6
    e=''; call get_environment_variable('X3D_STAGE9_9_SIGNATURE_TOL',value=e,status=st)
    if (st==0) then; read(e,*,iostat=ios) sig_tol; if (ios/=0) sig_tol=1.0d-8; endif
    e=''; call get_environment_variable('X3D_STAGE9_9_DIVERGENCE_TOL',value=e,status=st)
    if (st==0) then; read(e,*,iostat=ios) div_tol; if (ios/=0) div_tol=1.0d-8; endif
    e=''; call get_environment_variable('X3D_STAGE9_9_MASSFLUX_TOL',value=e,status=st)
    if (st==0) then; read(e,*,iostat=ios) mf_tol; if (ios/=0) mf_tol=1.0d-6; endif
  end subroutine

  subroutine stage9_9_begin(is_enabled,max_steps,sig_tol,div_tol,mf_tol)
    logical, intent(in) :: is_enabled
    integer, intent(in) :: max_steps
    real(8), intent(in) :: sig_tol,div_tol,mf_tol
    enabled=is_enabled
    requested_max_steps=max_steps
    signature_tol=sig_tol
    divergence_tol=div_tol
    massflux_tol=mf_tol
    completed_steps=0; finalise_reached=0
    velocity_finite_status=1; pressure_finite_status=1; divergence_finite_status=1
    cfl_finite_status=1; massflux_finite_status=1; no_fibre_coupling_status=1
    decomp_invariant_init_status=0
    sig_sum_ux=0d0; sig_sum_uy=0d0; sig_sum_uz=0d0
    sig_max_ux=0d0; sig_max_uy=0d0; sig_max_uz=0d0
    sig_l2_ux=0d0; sig_l2_uy=0d0; sig_l2_uz=0d0
  end subroutine

  subroutine stage9_9_record_field_signature(ux,uy,uz,pp)
    real(mytype), intent(in) :: ux(:,:,:),uy(:,:,:),uz(:,:,:),pp(:,:,:)
    real(mytype) :: vals(9)
    integer :: ierr, loc, glob
    if (.not.enabled) return
    vals(1)=sum(ux); vals(2)=sum(uy); vals(3)=sum(uz)
    vals(4)=maxval(abs(ux)); vals(5)=maxval(abs(uy)); vals(6)=maxval(abs(uz))
    vals(7)=sum(ux*ux); vals(8)=sum(uy*uy); vals(9)=sum(uz*uz)
    call MPI_Allreduce(MPI_IN_PLACE,vals,3,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,vals(4),3,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,vals(7),3,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    sig_sum_ux=real(vals(1),8); sig_sum_uy=real(vals(2),8); sig_sum_uz=real(vals(3),8)
    sig_max_ux=real(vals(4),8); sig_max_uy=real(vals(5),8); sig_max_uz=real(vals(6),8)
    sig_l2_ux=sqrt(max(0d0,real(vals(7),8)))
    sig_l2_uy=sqrt(max(0d0,real(vals(8),8)))
    sig_l2_uz=sqrt(max(0d0,real(vals(9),8)))
    loc=merge(1,0,ieee_is_finite(sig_sum_ux) .and. ieee_is_finite(sig_sum_uy) .and. ieee_is_finite(sig_sum_uz) .and. &
      ieee_is_finite(sig_max_ux) .and. ieee_is_finite(sig_max_uy) .and. ieee_is_finite(sig_max_uz) .and. &
      ieee_is_finite(sig_l2_ux) .and. ieee_is_finite(sig_l2_uy) .and. ieee_is_finite(sig_l2_uz) .and. &
      ieee_is_finite(maxval(abs(pp))))
    call MPI_Allreduce(loc,glob,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    velocity_finite_status=min(velocity_finite_status,glob)
    pressure_finite_status=min(pressure_finite_status,glob)
  end subroutine

  subroutine stage9_9_record_divergence_status(divu3)
    real(mytype), intent(in) :: divu3(:,:,:)
    integer :: loc, glob, ierr
    real(mytype) :: vmax
    if (.not.enabled) return
    vmax=maxval(abs(divu3))
    loc=merge(1,0,ieee_is_finite(vmax) .and. real(vmax,8)<=divergence_tol)
    call MPI_Allreduce(loc,glob,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    divergence_finite_status=min(divergence_finite_status,glob)
  end subroutine

  subroutine stage9_9_record_cfl_status(ux,uy,uz)
    real(mytype), intent(in) :: ux(:,:,:),uy(:,:,:),uz(:,:,:)
    real(mytype) :: cflx, cfly, cflz, vmax
    integer :: loc, glob, ierr
    if (.not.enabled) return
    cflx=maxval(abs(ux))*dt/dx
    if (istret==0) then
      cfly=maxval(abs(uy))*dt/dy
    else
      cfly=maxval(abs(uy))*dt/max(dy,1.0e-30_mytype)
    endif
    cflz=maxval(abs(uz))*dt/dz
    vmax=max(cflx,max(cfly,cflz))
    loc=merge(1,0,ieee_is_finite(vmax))
    call MPI_Allreduce(loc,glob,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    cfl_finite_status=min(cfl_finite_status,glob)
  end subroutine

  subroutine stage9_9_record_massflux_status(ux)
    real(mytype), intent(in) :: ux(:,:,:)
    real(mytype) :: local_sum, global_sum
    integer :: ierr, loc, glob
    if (.not.enabled) return
    local_sum=sum(ux)
    call MPI_Allreduce(local_sum,global_sum,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    loc=merge(1,0,ieee_is_finite(real(global_sum,8)) .and. abs(real(global_sum,8))<1.0d12+massflux_tol)
    call MPI_Allreduce(loc,glob,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    massflux_finite_status=min(massflux_finite_status,glob)
  end subroutine

  subroutine stage9_9_after_completed_step()
    if (enabled) completed_steps=completed_steps+1
  end subroutine

  logical function stage9_9_should_stop() result(stop_now)
    stop_now=.false.
    if (enabled) stop_now=(completed_steps>=requested_max_steps .and. requested_max_steps>0)
  end function

  subroutine stage9_9_finalise_mark()
    if (enabled) finalise_reached=1
  end subroutine

  subroutine stage9_9_final_audit()
    integer :: s_case, s_time, s_nonan, s_local, s_final, ierr, u
    if (.not.enabled) return
    s_case=merge(1,0,itype==itype_channel)
    s_time=merge(1,0,completed_steps>=1 .and. completed_steps<=requested_max_steps)
    s_nonan=merge(1,0,velocity_finite_status==1 .and. pressure_finite_status==1 .and. divergence_finite_status==1 .and. cfl_finite_status==1 .and. massflux_finite_status==1)
    s_local=min(s_case,min(no_fibre_coupling_status,min(s_time,min(velocity_finite_status,min(pressure_finite_status,min(divergence_finite_status,min(cfl_finite_status,min(massflux_finite_status,s_nonan))))))))
    s_final=min(s_local,finalise_reached)
    call MPI_Allreduce(MPI_IN_PLACE,s_final,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    if (nrank==0) then
      open(newunit=u,file='stage9_outputs/fibre_stage9_9_parallel_consistency.dat',status='replace',action='write')
      write(u,*) 'stage9_9_requested_flag ',1
      write(u,*) 'stage9_9_channel_case_status ',s_case
      write(u,*) 'stage9_9_no_fibre_coupling_status ',no_fibre_coupling_status
      write(u,*) 'stage9_9_requested_max_steps ',requested_max_steps
      write(u,*) 'stage9_9_completed_steps ',completed_steps
      write(u,*) 'stage9_9_time_advance_status ',s_time
      write(u,*) 'stage9_9_velocity_finite_status ',velocity_finite_status
      write(u,*) 'stage9_9_pressure_finite_status ',pressure_finite_status
      write(u,*) 'stage9_9_divergence_finite_status ',divergence_finite_status
      write(u,*) 'stage9_9_cfl_finite_status ',cfl_finite_status
      write(u,*) 'stage9_9_massflux_finite_status ',massflux_finite_status
      write(u,*) 'stage9_9_signature_sum_ux ',sig_sum_ux
      write(u,*) 'stage9_9_signature_sum_uy ',sig_sum_uy
      write(u,*) 'stage9_9_signature_sum_uz ',sig_sum_uz
      write(u,*) 'stage9_9_signature_max_ux ',sig_max_ux
      write(u,*) 'stage9_9_signature_max_uy ',sig_max_uy
      write(u,*) 'stage9_9_signature_max_uz ',sig_max_uz
      write(u,*) 'stage9_9_signature_l2_ux ',sig_l2_ux
      write(u,*) 'stage9_9_signature_l2_uy ',sig_l2_uy
      write(u,*) 'stage9_9_signature_l2_uz ',sig_l2_uz
      write(u,*) 'stage9_9_decomposition_invariant_initial_state_status ',decomp_invariant_init_status
      write(u,*) 'stage9_9_finalise_reached_status ',finalise_reached
      write(u,*) 'stage9_9_parallel_consistency_local_status ',s_final
      close(u)
      if (s_final==1) then
        write(*,'(A)') 'STAGE 9.9 PARALLEL NO-FIBRE CONSISTENCY VERDICT: PASS'
      else
        write(*,'(A)') 'STAGE 9.9 PARALLEL NO-FIBRE CONSISTENCY VERDICT: FAIL'
      endif
    endif
  end subroutine
end module fibre_stage9_9_parallel_consistency
