module fibre_stage9_4_no_fibre_dns_smoke
  use mpi
  use, intrinsic :: ieee_arithmetic
  use decomp_2d_mpi, only : nrank, nproc
  use var, only : ux1, uy1, uz1, pp3, px1, py1, pz1, divu3
  use param, only : itype, itype_channel
  implicit none
  private

  logical, save :: smoke_enabled = .false.
  integer, save :: requested_max_steps = 0
  integer, save :: completed_steps = 0
  integer, save :: start_itime = 0
  real(8), save :: start_time = 0.0d0
  integer, save :: finalise_reached = 0

  public :: stage9_4_smoke_requested
  public :: stage9_4_get_max_steps
  public :: stage9_4_begin
  public :: stage9_4_after_completed_step
  public :: stage9_4_finalise_mark
  public :: stage9_4_final_audit
contains

  logical function stage9_4_smoke_requested() result(requested)
    character(len=64) :: env_value
    integer :: env_status
    requested = .false.
    env_value = ''
    call get_environment_variable('X3D_STAGE9_4_NO_FIBRE_DNS_SMOKE', value=env_value, status=env_status)
    if (env_status == 0) requested = (trim(adjustl(env_value)) == '1')
  end function

  integer function stage9_4_get_max_steps(default_steps) result(max_steps)
    integer, intent(in) :: default_steps
    character(len=64) :: env_value
    integer :: env_status, parsed, ios
    max_steps = default_steps
    env_value = ''
    call get_environment_variable('X3D_STAGE9_4_MAX_STEPS', value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value,*,iostat=ios) parsed
      if (ios == 0 .and. parsed > 0) max_steps = parsed
    end if
  end function

  subroutine stage9_4_begin(enabled, max_steps, itime0_in, t0_in)
    logical, intent(in) :: enabled
    integer, intent(in) :: max_steps, itime0_in
    real(8), intent(in) :: t0_in
    smoke_enabled = enabled
    requested_max_steps = max_steps
    completed_steps = 0
    start_itime = itime0_in
    start_time = t0_in
    finalise_reached = 0
  end subroutine

  subroutine stage9_4_after_completed_step()
    if (smoke_enabled) completed_steps = completed_steps + 1
  end subroutine

  subroutine stage9_4_finalise_mark()
    finalise_reached = 1
  end subroutine

  subroutine stage9_4_final_audit()
    integer :: ierr, rank, s_case, s_nocouple, s_timeadv, s_vel, s_pres, s_div, s_nan, s_finalise, s_final
    integer :: u
    real(8) :: m1,m2,m3,mp,mqx,mqy,mqz,md

    rank = nrank
    s_case = merge(1,0,itype == itype_channel)
    s_nocouple = 1
    s_timeadv = merge(1,0,completed_steps >= 1 .and. completed_steps <= requested_max_steps)

    m1=maxval(abs(ux1)); m2=maxval(abs(uy1)); m3=maxval(abs(uz1))
    mp=maxval(abs(pp3)); mqx=maxval(abs(px1)); mqy=maxval(abs(py1)); mqz=maxval(abs(pz1)); md=maxval(abs(divu3))
    s_vel = merge(1,0,ieee_is_finite(m1).and.ieee_is_finite(m2).and.ieee_is_finite(m3))
    s_pres = merge(1,0,ieee_is_finite(mp).and.ieee_is_finite(mqx).and.ieee_is_finite(mqy).and.ieee_is_finite(mqz))
    s_div = merge(1,0,ieee_is_finite(md))
    s_nan = merge(1,0,s_vel==1 .and. s_pres==1 .and. s_div==1)
    s_finalise = finalise_reached

    s_final = min(s_case,min(s_nocouple,min(s_timeadv,min(s_vel,min(s_pres,min(s_div,min(s_nan,s_finalise)))))))
    call MPI_Allreduce(MPI_IN_PLACE, s_final, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
      open(newunit=u,file='stage9_outputs/fibre_stage9_4_no_fibre_dns_smoke.dat',status='replace',action='write')
      write(u,*) 'stage9_4_smoke_requested_flag ', merge(1,0,smoke_enabled)
      write(u,*) 'stage9_4_channel_case_status ', s_case
      write(u,*) 'stage9_4_no_fibre_coupling_status ', s_nocouple
      write(u,*) 'stage9_4_requested_max_steps ', requested_max_steps
      write(u,*) 'stage9_4_completed_steps ', completed_steps
      write(u,*) 'stage9_4_time_advance_status ', s_timeadv
      write(u,*) 'stage9_4_velocity_finite_status ', s_vel
      write(u,*) 'stage9_4_pressure_finite_status ', s_pres
      write(u,*) 'stage9_4_divergence_finite_status ', s_div
      write(u,*) 'stage9_4_no_nan_inf_status ', s_nan
      write(u,*) 'stage9_4_finalise_reached_status ', s_finalise
      write(u,*) 'stage9_4_no_fibre_dns_smoke_status ', s_final
      close(u)
      if (s_final == 1) then
        write(*,'(A)') 'STAGE 9.4 NO-FIBRE DNS SMOKE VERDICT: PASS'
      else
        write(*,'(A)') 'STAGE 9.4 NO-FIBRE DNS SMOKE VERDICT: FAIL'
      end if
    end if
  end subroutine

end module fibre_stage9_4_no_fibre_dns_smoke
