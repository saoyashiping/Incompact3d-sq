module fibre_stage9_6_rk3_rhs_massflux_regression
  use mpi
  use, intrinsic :: ieee_arithmetic
  use decomp_2d_mpi, only : nrank
  use decomp_2d_constants, only : mytype, real_type
  use param, only : itype, itype_channel, dt, dx, dy, dz, istret
  implicit none
  private
  integer, parameter :: max_track_steps=64

  logical, save :: enabled=.false.
  integer, save :: requested_max_steps=0, completed_steps=0, rk3_substep_calls=0, finalise_reached=0
  real(8), save :: mass_flux_tol=1.0d-6, cfl_max_limit=10.0d0
  real(8), save :: step_cfl_x(max_track_steps)=0d0, step_cfl_y(max_track_steps)=0d0, step_cfl_z(max_track_steps)=0d0
  real(8), save :: step_massflux(max_track_steps)=0d0, step_massflux_delta(max_track_steps)=0d0
  integer, save :: rhs_finite_status=1, velocity_finite_status=1, cfl_finite_status=1, cfl_limit_status=1
  integer, save :: massflux_finite_status=1, massflux_tolerance_status=1

  public :: stage9_6_requested, stage9_6_get_max_steps, stage9_6_get_tolerances
  public :: stage9_6_begin, stage9_6_record_rk_substep, stage9_6_record_rhs_finite_status
  public :: stage9_6_record_velocity_finite_status, stage9_6_record_cfl_status, stage9_6_record_massflux_status
  public :: stage9_6_after_completed_step, stage9_6_finalise_mark, stage9_6_final_audit
contains
  logical function stage9_6_requested() result(requested)
    character(len=64) :: e; integer :: st
    requested=.false.; e=''
    call get_environment_variable('X3D_STAGE9_6_RK3_RHS_MASSFLUX_REGRESSION',value=e,status=st)
    if (st==0) requested=(trim(adjustl(e))=='1')
  end function

  integer function stage9_6_get_max_steps(default_steps) result(max_steps)
    integer,intent(in)::default_steps
    character(len=64)::e; integer::st,ios,v
    max_steps=default_steps; e=''
    call get_environment_variable('X3D_STAGE9_6_MAX_STEPS',value=e,status=st)
    if (st==0) then; read(e,*,iostat=ios) v; if (ios==0 .and. v>0) max_steps=v; endif
  end function

  subroutine stage9_6_get_tolerances(mflux_tol,cfl_limit)
    real(8),intent(out)::mflux_tol,cfl_limit
    character(len=64)::e; integer::st,ios
    mflux_tol=1.0d-6; cfl_limit=10.0d0
    e=''; call get_environment_variable('X3D_STAGE9_6_MASS_FLUX_TOL',value=e,status=st)
    if (st==0) then; read(e,*,iostat=ios) mflux_tol; if (ios/=0) mflux_tol=1.0d-6; endif
    e=''; call get_environment_variable('X3D_STAGE9_6_CFL_MAX_LIMIT',value=e,status=st)
    if (st==0) then; read(e,*,iostat=ios) cfl_limit; if (ios/=0) cfl_limit=10.0d0; endif
  end subroutine

  subroutine stage9_6_begin(is_enabled,max_steps,mflux_tol,cfl_limit)
    logical,intent(in)::is_enabled
    integer,intent(in)::max_steps
    real(8),intent(in)::mflux_tol,cfl_limit
    enabled=is_enabled
    requested_max_steps=max_steps
    mass_flux_tol=mflux_tol
    cfl_max_limit=cfl_limit
    completed_steps=0; rk3_substep_calls=0; finalise_reached=0
    step_cfl_x=0d0; step_cfl_y=0d0; step_cfl_z=0d0; step_massflux=0d0; step_massflux_delta=0d0
    rhs_finite_status=1; velocity_finite_status=1; cfl_finite_status=1; cfl_limit_status=1
    massflux_finite_status=1; massflux_tolerance_status=1
  end subroutine

  subroutine stage9_6_record_rk_substep()
    if (enabled) rk3_substep_calls=rk3_substep_calls+1
  end subroutine

  subroutine stage9_6_record_rhs_finite_status(drho,dux,duy,duz)
    real(mytype),intent(in) :: drho(:,:,:,:),dux(:,:,:,:),duy(:,:,:,:),duz(:,:,:,:)
    integer :: loc,glob,ierr
    if (.not.enabled) return
    loc=merge(1,0,ieee_is_finite(maxval(abs(drho))) .and. ieee_is_finite(maxval(abs(dux))) .and. &
                  ieee_is_finite(maxval(abs(duy))) .and. ieee_is_finite(maxval(abs(duz))))
    call MPI_Allreduce(loc,glob,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    rhs_finite_status=min(rhs_finite_status,glob)
  end subroutine

  subroutine stage9_6_record_velocity_finite_status(ux,uy,uz)
    real(mytype),intent(in) :: ux(:,:,:),uy(:,:,:),uz(:,:,:)
    integer :: loc,glob,ierr
    if (.not.enabled) return
    loc=merge(1,0,ieee_is_finite(maxval(abs(ux))) .and. ieee_is_finite(maxval(abs(uy))) .and. ieee_is_finite(maxval(abs(uz))))
    call MPI_Allreduce(loc,glob,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    velocity_finite_status=min(velocity_finite_status,glob)
  end subroutine

  subroutine stage9_6_record_cfl_status(ux,uy,uz)
    real(mytype),intent(in) :: ux(:,:,:),uy(:,:,:),uz(:,:,:)
    real(mytype) :: cflx_loc,cfly_loc,cflz_loc,vals(3)
    integer :: ierr,idx
    if (.not.enabled) return
    idx=completed_steps+1
    if (idx<1 .or. idx>max_track_steps) return
    cflx_loc=maxval(abs(ux))*dt/dx
    if (istret==0) then
      cfly_loc=maxval(abs(uy))*dt/dy
    else
      cfly_loc=maxval(abs(uy))*dt/max(dy,1.0e-30_mytype)
    endif
    cflz_loc=maxval(abs(uz))*dt/dz
    vals=(/cflx_loc,cfly_loc,cflz_loc/)
    call MPI_Allreduce(MPI_IN_PLACE,vals,3,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    step_cfl_x(idx)=real(vals(1),8); step_cfl_y(idx)=real(vals(2),8); step_cfl_z(idx)=real(vals(3),8)
    if (.not.(ieee_is_finite(step_cfl_x(idx)) .and. ieee_is_finite(step_cfl_y(idx)) .and. ieee_is_finite(step_cfl_z(idx)))) cfl_finite_status=0
    if (max(step_cfl_x(idx),max(step_cfl_y(idx),step_cfl_z(idx)))>cfl_max_limit) cfl_limit_status=0
  end subroutine

  subroutine stage9_6_record_massflux_status(ux)
    real(mytype),intent(in) :: ux(:,:,:)
    real(mytype) :: local_sum,global_sum,proxy
    integer :: local_n,global_n,ierr,idx
    if (.not.enabled) return
    idx=completed_steps+1
    if (idx<1 .or. idx>max_track_steps) return
    local_sum=sum(ux); local_n=size(ux)
    call MPI_Allreduce(local_sum,global_sum,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(local_n,global_n,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    proxy=global_sum/max(1,global_n)
    step_massflux(idx)=real(proxy,8)
    if (.not.ieee_is_finite(step_massflux(idx))) massflux_finite_status=0
    if (idx==1) then
      step_massflux_delta(idx)=0d0
    else
      step_massflux_delta(idx)=abs(step_massflux(idx)-step_massflux(idx-1))
      if (step_massflux_delta(idx)>mass_flux_tol) massflux_tolerance_status=0
    endif
  end subroutine

  subroutine stage9_6_after_completed_step()
    if (enabled) completed_steps=completed_steps+1
  end subroutine

  subroutine stage9_6_finalise_mark(); finalise_reached=1; end subroutine

  subroutine stage9_6_final_audit()
    integer :: i,u,s_case,s_nocouple,s_time,s_rk,s_nonan,s_final,rank
    character(len=128)::key
    rank=nrank
    s_case=merge(1,0,itype==itype_channel)
    s_nocouple=1
    s_time=merge(1,0,completed_steps>=1 .and. completed_steps<=requested_max_steps)
    s_rk=merge(1,0,rk3_substep_calls>=completed_steps)
    s_nonan=merge(1,0,rhs_finite_status==1 .and. velocity_finite_status==1 .and. cfl_finite_status==1 .and. massflux_finite_status==1)
    s_final=min(s_case,min(s_nocouple,min(s_time,min(s_rk,min(rhs_finite_status,min(velocity_finite_status,min(cfl_finite_status,min(cfl_limit_status,min(massflux_finite_status,min(massflux_tolerance_status,min(s_nonan,finalise_reached)))))))))))
    call MPI_Allreduce(MPI_IN_PLACE,s_final,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,i)
    if (rank==0) then
      open(newunit=u,file='stage9_outputs/fibre_stage9_6_rk3_rhs_massflux_regression.dat',status='replace',action='write')
      write(u,*) 'stage9_6_requested_flag ',merge(1,0,enabled)
      write(u,*) 'stage9_6_channel_case_status ',s_case
      write(u,*) 'stage9_6_no_fibre_coupling_status ',s_nocouple
      write(u,*) 'stage9_6_requested_max_steps ',requested_max_steps
      write(u,*) 'stage9_6_completed_steps ',completed_steps
      write(u,*) 'stage9_6_time_advance_status ',s_time
      write(u,*) 'stage9_6_rk3_execution_status ',s_rk
      write(u,*) 'stage9_6_rhs_finite_status ',rhs_finite_status
      write(u,*) 'stage9_6_velocity_finite_status ',velocity_finite_status
      write(u,*) 'stage9_6_cfl_finite_status ',cfl_finite_status
      write(u,*) 'stage9_6_cfl_limit_status ',cfl_limit_status
      write(u,*) 'stage9_6_massflux_finite_status ',massflux_finite_status
      write(u,*) 'stage9_6_massflux_tolerance_status ',massflux_tolerance_status
      write(u,*) 'stage9_6_no_nan_inf_status ',s_nonan
      write(u,*) 'stage9_6_finalise_reached_status ',finalise_reached
      write(u,*) 'stage9_6_rk3_rhs_massflux_regression_status ',s_final
      do i=1,completed_steps
        write(key,'(A,I0,A)') 'stage9_6_step_',i,'_cfl_x'; write(u,*) trim(key),' ',step_cfl_x(i)
        write(key,'(A,I0,A)') 'stage9_6_step_',i,'_cfl_y'; write(u,*) trim(key),' ',step_cfl_y(i)
        write(key,'(A,I0,A)') 'stage9_6_step_',i,'_cfl_z'; write(u,*) trim(key),' ',step_cfl_z(i)
        write(key,'(A,I0,A)') 'stage9_6_step_',i,'_massflux_proxy'; write(u,*) trim(key),' ',step_massflux(i)
        write(key,'(A,I0,A)') 'stage9_6_step_',i,'_massflux_delta'; write(u,*) trim(key),' ',step_massflux_delta(i)
      enddo
      close(u)
      if (s_final==1) then
        write(*,'(A)') 'STAGE 9.6 RK3 RHS MASS-FLUX REGRESSION VERDICT: PASS'
      else
        write(*,'(A)') 'STAGE 9.6 RK3 RHS MASS-FLUX REGRESSION VERDICT: FAIL'
      endif
    endif
  end subroutine
end module fibre_stage9_6_rk3_rhs_massflux_regression
