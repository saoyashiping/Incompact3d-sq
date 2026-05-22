module fibre_stage9_3_channel_init_dryrun
  use mpi
  use, intrinsic :: ieee_arithmetic
  use decomp_2d
  use decomp_2d_mpi, only : nrank
  use var, only : ux1,uy1,uz1,rho1,pp3,px1,py1,pz1,divu3,ph1,ph2,ph3,ph4,phG
  use param, only : itype, itype_channel
  use variables, only : nx,ny,nz,p_row,p_col
  implicit none
  private
  public :: stage9_3_dryrun_requested, stage9_3_channel_init_dryrun_audit
contains

  logical function stage9_3_dryrun_requested() result(requested)
    implicit none
    character(len=64) :: env_value
    integer :: env_status

    requested = .false.
    env_value = ''
    call get_environment_variable('X3D_STAGE9_3_CHANNEL_INIT_DRYRUN', value=env_value, status=env_status)
    if (env_status == 0) then
      requested = trim(adjustl(env_value)) == '1'
    end if
  end function stage9_3_dryrun_requested

  subroutine stage9_3_channel_init_dryrun_audit()
    implicit none
    integer :: ierr, rank
    integer :: s_dry,s_case,s_grid,s_decomp,s_stat,s_visu,s_ph1,s_ph2,s_ph3,s_ph4,s_phG,s_alloc,s_vel,s_pres,s_div,s_notime,s_nocoupling,s_final
    real(8) :: m1,m2,m3,mp,mqx,mqy,mqz,md
    integer :: u

    rank = nrank
    s_dry=1; s_case=merge(1,0,itype==itype_channel); s_grid=merge(1,0,nx>0.and.ny>0.and.nz>0.and.p_row>=0.and.p_col>=0)
    s_decomp=1
    if (.not.(all(xsize>=0).and.all(ysize>=0).and.all(zsize>=0))) s_decomp=0
    if (.not.all((xsize==0).or.(xstart<=xend))) s_decomp=0
    if (.not.all((ysize==0).or.(ystart<=yend))) s_decomp=0
    if (.not.all((zsize==0).or.(zstart<=zend))) s_decomp=0
    if (.not.all((xsize==0).or.(xsize==xend-xstart+1))) s_decomp=0
    s_stat=merge(1,0,all(xstS>=1).and.all(xenS>=xstS-1))
    s_visu=merge(1,0,all(xstV>=1).and.all(xenV>=xstV-1))
    s_ph1=merge(1,0,all(ph1%xst<=ph1%xen).and.all(ph1%yst<=ph1%yen).and.all(ph1%zst<=ph1%zen))
    s_ph2=merge(1,0,all(ph2%xst<=ph2%xen).and.all(ph2%yst<=ph2%yen).and.all(ph2%zst<=ph2%zen))
    s_ph3=merge(1,0,all(ph3%xst<=ph3%xen).and.all(ph3%yst<=ph3%yen).and.all(ph3%zst<=ph3%zen))
    s_ph4=merge(1,0,all(ph4%xst<=ph4%xen).and.all(ph4%yst<=ph4%yen).and.all(ph4%zst<=ph4%zen))
    s_phG=merge(1,0,all(phG%xst<=phG%xen).and.all(phG%yst<=phG%yen).and.all(phG%zst<=phG%zen))
    s_alloc=merge(1,0,allocated(ux1).and.allocated(uy1).and.allocated(uz1).and.allocated(rho1).and.allocated(pp3).and.allocated(px1).and.allocated(py1).and.allocated(pz1).and.allocated(divu3))
    m1=maxval(abs(ux1)); m2=maxval(abs(uy1)); m3=maxval(abs(uz1))
    mp=maxval(abs(pp3)); mqx=maxval(abs(px1)); mqy=maxval(abs(py1)); mqz=maxval(abs(pz1)); md=maxval(abs(divu3))
    s_vel=merge(1,0,ieee_is_finite(m1).and.ieee_is_finite(m2).and.ieee_is_finite(m3))
    s_pres=merge(1,0,ieee_is_finite(mp).and.ieee_is_finite(mqx).and.ieee_is_finite(mqy).and.ieee_is_finite(mqz))
    s_div=merge(1,0,allocated(divu3).and.ieee_is_finite(md))
    s_notime=1; s_nocoupling=1
    s_final=min(s_dry,min(s_case,min(s_grid,min(s_decomp,min(s_stat,min(s_visu,min(s_ph1,min(s_ph2,min(s_ph3,min(s_ph4,min(s_phG,min(s_alloc,min(s_vel,min(s_pres,min(s_div,min(s_notime,s_nocoupling))))))))))))))))
    call MPI_Allreduce(MPI_IN_PLACE,s_final,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
    if (rank==0) then
      open(newunit=u,file='stage9_outputs/fibre_stage9_3_channel_init_dryrun.dat',status='replace',action='write')
      write(u,*) 'stage9_3_dryrun_requested_flag ',s_dry
      write(u,*) 'stage9_3_channel_case_status ',s_case
      write(u,*) 'stage9_3_global_grid_status ',s_grid
      write(u,*) 'stage9_3_local_decomposition_status ',s_decomp
      write(u,*) 'stage9_3_coarse_stat_bounds_status ',s_stat
      write(u,*) 'stage9_3_visu_bounds_status ',s_visu
      write(u,*) 'stage9_3_ph1_descriptor_status ',s_ph1
      write(u,*) 'stage9_3_ph2_descriptor_status ',s_ph2
      write(u,*) 'stage9_3_ph3_descriptor_status ',s_ph3
      write(u,*) 'stage9_3_ph4_descriptor_status ',s_ph4
      write(u,*) 'stage9_3_phG_descriptor_status ',s_phG
      write(u,*) 'stage9_3_core_arrays_allocated_status ',s_alloc
      write(u,*) 'stage9_3_velocity_finite_status ',s_vel
      write(u,*) 'stage9_3_pressure_finite_status ',s_pres
      write(u,*) 'stage9_3_divergence_finite_status ',s_div
      write(u,*) 'stage9_3_no_time_advance_status ',s_notime
      write(u,*) 'stage9_3_no_fibre_coupling_status ',s_nocoupling
      write(u,*) 'stage9_3_channel_init_dryrun_status ',s_final
      close(u)
      if (s_final==1) then
        write(*,'(A)') 'STAGE 9.3 CHANNEL INIT DRY-RUN VERDICT: PASS'
      else
        write(*,'(A)') 'STAGE 9.3 CHANNEL INIT DRY-RUN VERDICT: FAIL'
      endif
    endif
  end subroutine stage9_3_channel_init_dryrun_audit

end module fibre_stage9_3_channel_init_dryrun
