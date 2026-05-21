module fibre_stage9_5_projection_regression
  use mpi
  use, intrinsic :: ieee_arithmetic
  use decomp_2d_mpi, only : nrank
  use decomp_2d_constants, only : mytype
  use param, only : itype, itype_channel
  implicit none
  private
  integer, parameter :: max_track_steps=64

  logical, save :: proj_enabled=.false.
  integer, save :: requested_max_steps=0, completed_steps=0, projection_samples=0, finalise_reached=0
  real(8), save :: div_max_tol=1.0d-8, div_mean_tol=1.0d-9
  real(8), save :: before_max(max_track_steps)=0d0, before_mean(max_track_steps)=0d0
  real(8), save :: after_max(max_track_steps)=0d0, after_mean(max_track_steps)=0d0
  integer, save :: pressure_finite_status=1, velocity_finite_status=1

  public :: stage9_5_projection_requested, stage9_5_get_max_steps, stage9_5_get_divergence_tolerances
  public :: stage9_5_begin, stage9_5_record_divergence_before_projection, stage9_5_record_divergence_after_projection, stage9_5_record_projection_divergence_pair
  public :: stage9_5_record_pressure_finite_status, stage9_5_after_completed_step, stage9_5_finalise_mark, stage9_5_final_audit
contains
  logical function stage9_5_projection_requested() result(requested)
    character(len=64) :: e; integer :: st
    requested=.false.; e=''
    call get_environment_variable('X3D_STAGE9_5_PROJECTION_REGRESSION',value=e,status=st)
    if (st==0) requested=(trim(adjustl(e))=='1')
  end function

  integer function stage9_5_get_max_steps(default_steps) result(max_steps)
    integer,intent(in)::default_steps
    character(len=64)::e; integer::st,ios,v
    max_steps=default_steps; e=''
    call get_environment_variable('X3D_STAGE9_5_MAX_STEPS',value=e,status=st)
    if (st==0) then; read(e,*,iostat=ios) v; if (ios==0 .and. v>0) max_steps=v; endif
  end function

  subroutine stage9_5_get_divergence_tolerances(max_tol,mean_tol)
    real(8),intent(out)::max_tol,mean_tol
    character(len=64)::e; integer::st,ios
    max_tol=1.0d-8; mean_tol=1.0d-9
    e=''; call get_environment_variable('X3D_STAGE9_5_DIV_MAX_TOL',value=e,status=st)
    if (st==0) then; read(e,*,iostat=ios) max_tol; if (ios/=0) max_tol=1.0d-8; endif
    e=''; call get_environment_variable('X3D_STAGE9_5_DIV_MEAN_TOL',value=e,status=st)
    if (st==0) then; read(e,*,iostat=ios) mean_tol; if (ios/=0) mean_tol=1.0d-9; endif
  end subroutine

  subroutine stage9_5_begin(enabled,max_steps,max_tol,mean_tol)
    logical,intent(in)::enabled; integer,intent(in)::max_steps; real(8),intent(in)::max_tol,mean_tol
    proj_enabled=enabled; requested_max_steps=max_steps; div_max_tol=max_tol; div_mean_tol=mean_tol
    completed_steps=0; projection_samples=0; finalise_reached=0; pressure_finite_status=1; velocity_finite_status=1
    before_max=0d0; before_mean=0d0; after_max=0d0; after_mean=0d0
  end subroutine

  subroutine reduce_div_stats(div,maxv,meanv)
    real(8),intent(in) :: div(:,:,:)
    real(8),intent(out)::maxv,meanv
    real(8)::lmax,lsum,gsum; integer(kind=8)::ln,gn; integer::ierr
    lmax=maxval(abs(div)); lsum=sum(abs(div)); ln=size(div,kind=8)
    call MPI_Allreduce(lmax,maxv,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(lsum,gsum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(ln,gn,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
    if (gn>0_8) then; meanv=gsum/dble(gn); else; meanv=0d0; endif
  end subroutine

  subroutine stage9_5_record_divergence_before_projection(div)
    real(8),intent(in)::div(:,:,:)
    if (.not.proj_enabled) return
    if (projection_samples>=max_track_steps) return
    projection_samples=projection_samples+1
    call reduce_div_stats(div,before_max(projection_samples),before_mean(projection_samples))
  end subroutine

  subroutine stage9_5_record_divergence_after_projection(div)
    real(8),intent(in)::div(:,:,:)
    if (.not.proj_enabled) return
    if (projection_samples<=0 .or. projection_samples>max_track_steps) return
    call reduce_div_stats(div,after_max(projection_samples),after_mean(projection_samples))
  end subroutine


  subroutine stage9_5_record_projection_divergence_pair(nlock, div_max, div_mean)
    integer, intent(in) :: nlock
    real(mytype), intent(in) :: div_max, div_mean
    if (.not.proj_enabled) return
    if (nlock == 1) then
      if (projection_samples>=max_track_steps) return
      projection_samples = projection_samples + 1
      before_max(projection_samples) = real(div_max,8)
      before_mean(projection_samples) = real(div_mean,8)
    else if (nlock == 2) then
      if (projection_samples<=0 .or. projection_samples>max_track_steps) return
      after_max(projection_samples) = real(div_max,8)
      after_mean(projection_samples) = real(div_mean,8)
    end if
  end subroutine
  subroutine stage9_5_record_pressure_finite_status(pp,px,py,pz,ux,uy,uz)
    real(8),intent(in)::pp(:,:,:,:),px(:,:,:),py(:,:,:),pz(:,:,:),ux(:,:,:),uy(:,:,:),uz(:,:,:)
    integer :: locp,locv,gp,gv,ierr
    if (.not.proj_enabled) return
    locp=merge(1,0,ieee_is_finite(maxval(abs(pp))).and.ieee_is_finite(maxval(abs(px))).and.ieee_is_finite(maxval(abs(py))).and.ieee_is_finite(maxval(abs(pz))))
    locv=merge(1,0,ieee_is_finite(maxval(abs(ux))).and.ieee_is_finite(maxval(abs(uy))).and.ieee_is_finite(maxval(abs(uz))))
    call MPI_Allreduce(locp,gp,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(locv,gv,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    pressure_finite_status=min(pressure_finite_status,gp)
    velocity_finite_status=min(velocity_finite_status,gv)
  end subroutine

  subroutine stage9_5_after_completed_step()
    if (proj_enabled) completed_steps=completed_steps+1
  end subroutine
  subroutine stage9_5_finalise_mark(); finalise_reached=1; end subroutine

  subroutine stage9_5_final_audit()
    integer::i,u,s_case,s_nocouple,s_time,s_path,s_bfin,s_afin,s_reduc,s_amax,s_amean,s_nonan,s_final,rank
    character(len=128) :: key
    real(8)::tinyv,rrm,rra
    rank=nrank; tinyv=1.0d-300
    s_case=merge(1,0,itype==itype_channel); s_nocouple=1; s_time=merge(1,0,completed_steps>=1 .and. completed_steps<=requested_max_steps)
    s_path=merge(1,0,projection_samples>=1)
    s_bfin=1; s_afin=1; s_reduc=1; s_amax=1; s_amean=1
    do i=1,projection_samples
      if (.not.ieee_is_finite(before_max(i)) .or. .not.ieee_is_finite(before_mean(i))) s_bfin=0
      if (.not.ieee_is_finite(after_max(i)) .or. .not.ieee_is_finite(after_mean(i))) s_afin=0
      rra=after_max(i)/max(before_max(i),tinyv); rrm=after_mean(i)/max(before_mean(i),tinyv)
      if (rra>1.d0 .or. rrm>1.d0) s_reduc=0
      if (after_max(i)>div_max_tol) s_amax=0
      if (after_mean(i)>div_mean_tol) s_amean=0
    end do
    s_nonan=merge(1,0,pressure_finite_status==1 .and. velocity_finite_status==1 .and. s_bfin==1 .and. s_afin==1)
    s_final=min(s_case,min(s_nocouple,min(s_time,min(s_path,min(s_bfin,min(s_afin,min(s_reduc,min(s_amax,min(s_amean,min(pressure_finite_status,min(velocity_finite_status,min(s_nonan,finalise_reached))))))))))))
    call MPI_Allreduce(MPI_IN_PLACE,s_final,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,i)
    if (rank==0) then
      open(newunit=u,file='stage9_outputs/fibre_stage9_5_projection_regression.dat',status='replace',action='write')
      write(u,*) 'stage9_5_projection_requested_flag ',merge(1,0,proj_enabled)
      write(u,*) 'stage9_5_channel_case_status ',s_case
      write(u,*) 'stage9_5_no_fibre_coupling_status ',s_nocouple
      write(u,*) 'stage9_5_requested_max_steps ',requested_max_steps
      write(u,*) 'stage9_5_completed_steps ',completed_steps
      write(u,*) 'stage9_5_div_max_tolerance ',div_max_tol
      write(u,*) 'stage9_5_div_mean_tolerance ',div_mean_tol
      write(u,*) 'stage9_5_time_advance_status ',s_time
      write(u,*) 'stage9_5_projection_path_executed_status ',s_path
      write(u,*) 'stage9_5_divergence_before_finite_status ',s_bfin
      write(u,*) 'stage9_5_divergence_after_finite_status ',s_afin
      write(u,*) 'stage9_5_divergence_reduction_status ',s_reduc
      write(u,*) 'stage9_5_divergence_after_max_status ',s_amax
      write(u,*) 'stage9_5_divergence_after_mean_status ',s_amean
      write(u,*) 'stage9_5_pressure_finite_status ',pressure_finite_status
      write(u,*) 'stage9_5_velocity_finite_status ',velocity_finite_status
      write(u,*) 'stage9_5_no_nan_inf_status ',s_nonan
      write(u,*) 'stage9_5_finalise_reached_status ',finalise_reached
      write(u,*) 'stage9_5_projection_regression_status ',s_final
      do i=1,projection_samples
        rra=after_max(i)/max(before_max(i),tinyv); rrm=after_mean(i)/max(before_mean(i),tinyv)
        write(key,'(A,I0,A)') 'stage9_5_step_',i,'_div_before_max'
        write(u,*) trim(key),' ',before_max(i)
        write(key,'(A,I0,A)') 'stage9_5_step_',i,'_div_before_mean'
        write(u,*) trim(key),' ',before_mean(i)
        write(key,'(A,I0,A)') 'stage9_5_step_',i,'_div_after_max'
        write(u,*) trim(key),' ',after_max(i)
        write(key,'(A,I0,A)') 'stage9_5_step_',i,'_div_after_mean'
        write(u,*) trim(key),' ',after_mean(i)
        write(key,'(A,I0,A)') 'stage9_5_step_',i,'_div_reduction_ratio_max'
        write(u,*) trim(key),' ',rra
        write(key,'(A,I0,A)') 'stage9_5_step_',i,'_div_reduction_ratio_mean'
        write(u,*) trim(key),' ',rrm
      end do
      close(u)
      if (s_final==1) then
        write(*,'(A)') 'STAGE 9.5 PRESSURE PROJECTION REGRESSION VERDICT: PASS'
      else
        write(*,'(A)') 'STAGE 9.5 PRESSURE PROJECTION REGRESSION VERDICT: FAIL'
      endif
    endif
  end subroutine
end module fibre_stage9_5_projection_regression
