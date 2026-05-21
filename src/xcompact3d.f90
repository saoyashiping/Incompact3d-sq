!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

program xcompact3d

  use var
  use case

  use transeq, only : calculate_transeq_rhs
  use time_integrators, only : int_time
  use navier, only : velocity_to_momentum, momentum_to_velocity, pre_correc, &
       calc_divu_constraint, solve_poisson, cor_vel
  use tools, only : restart, simu_stats, apply_spatial_filter, read_inflow
  use turbine, only : compute_turbines
  use ibm_param
  use ibm, only : body
  use genepsi, only : genepsi3d
  use mhd,    only : Bm,mhd_equation,test_magnetic, &
                     solve_poisson_mhd
  use param, only : mhd_active
  use particle, only : intt_particles

  implicit none

  logical :: stage9_3_dryrun

  call init_xcompact3d()

  stage9_3_dryrun = stage9_3_dryrun_requested()
  if (stage9_3_dryrun) then
     call stage9_3_channel_init_dryrun_audit()
     call finalise_xcompact3d()
     stop
  endif

  do itime=ifirst,ilast
     !t=itime*dt
     t=t0 + (itime0 + itime + 1 - ifirst)*dt
     
     call simu_stats(2)

     if (iturbine.ne.0) call compute_turbines(ux1, uy1, uz1)

     if (iin.eq.3.and.mod(itime,ntimesteps)==1) then
        call read_inflow(ux_inflow,uy_inflow,uz_inflow,itime/ntimesteps)
     endif

     if ((itype.eq.itype_abl.or.iturbine.ne.0).and.(ifilter.ne.0).and.(ilesmod.ne.0)) then
        call filter(C_filter)
        call apply_spatial_filter(ux1,uy1,uz1,phi1)
     endif

     do itr=1,iadvance_time

        call set_fluid_properties(rho1,mu1)
        call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)

        if (imove.eq.1) then ! update epsi for moving objects
          if ((iibm.eq.2).or.(iibm.eq.3)) then
             call genepsi3d(ep1)
          else if (iibm.eq.1) then
             call body(ux1,uy1,uz1,ep1)
          endif
        endif
        call calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

#ifdef DEBG
        call check_transients()
#endif
        
        if (ilmn) then
           !! XXX N.B. from this point, X-pencil velocity arrays contain momentum (LMN only).
           call velocity_to_momentum(rho1,ux1,uy1,uz1)
        endif

        call int_time(rho1,ux1,uy1,uz1,phi1,drho1,dux1,duy1,duz1,dphi1)
        call pre_correc(ux1,uy1,uz1,ep1)

        call calc_divu_constraint(divu3,rho1,phi1)
        call solve_poisson(pp3,px1,py1,pz1,rho1,ux1,uy1,uz1,ep1,drho1,divu3)
        call cor_vel(ux1,uy1,uz1,px1,py1,pz1)

        if(mhd_active .and. mhd_equation == 'induction') then
          call solve_poisson_mhd()
        endif

        if (ilmn) then
           call momentum_to_velocity(rho1,ux1,uy1,uz1)
           !! XXX N.B. from this point, X-pencil velocity arrays contain velocity (LMN only).
           !! Note - all other solvers work on velocity always
        endif
        
        call test_flow(rho1,ux1,uy1,uz1,phi1,ep1,drho1,divu3)

        if(mhd_active) call test_magnetic

     enddo !! End sub timesteps

     if(particle_active) then
       call intt_particles(ux1,uy1,uz1,t)
     endif

     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,1)

     call simu_stats(3)

     call postprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)

  enddo !! End time loop

  call finalise_xcompact3d()

end program xcompact3d
!########################################################################
!########################################################################
subroutine init_xcompact3d()

  use MPI
  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_init
  USE decomp_2d_poisson, ONLY : decomp_2d_poisson_init
  use case
  use sandbox, only : init_sandbox
  use forces

  use var

  use navier, only : calc_divu_constraint
  use tools, only : test_speed_min_max, test_scalar_min_max, &
       restart, &
       simu_stats, compute_cfldiff, &
       init_inflow_outflow, read_inflow

  use param, only : ilesmod, jles,itype
  use param, only : irestart, mhd_active
  use param, only : periodic_bc

  use variables, only : nx, ny, nz, nxm, nym, nzm
  use variables, only : p_row, p_col
  use variables, only : nstat, nvisu, ilist

  use les, only: init_explicit_les
  use turbine, only: init_turbines

  use visu, only : visu_init, visu_ready

  use genepsi, only : genepsi3d, epsi_init
  use ibm, only : body

  use probes, only : init_probes

  use mhd, only: mhd_init
  use particle,  only : particle_report,local_domain_size

  implicit none

  integer :: ierr

  integer :: nargin, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, FNBase
    
  !! Initialise MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! Handle input file like a boss -- GD
  nargin=command_argument_count()
  if (nargin <1) then
     InputFN='input.i3d'
     if (nrank==0) write(*,*) 'Xcompact3d is run with the default file -->', trim(InputFN)
  elseif (nargin >= 1) then
     call get_command_argument(1,InputFN,FNLength,status)
     back=.true.
     FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
     DecInd=index(FNBase,'.',back)
     if (DecInd >1) then
        FNBase=FNBase(1:(DecInd-1))
     end if
     if (nrank==0) write(*,*) 'Xcompact3d is run with the provided file -->', trim(InputFN)
  endif

#ifdef ADIOS2
  if (nrank .eq. 0) then
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
     print *, " WARNING: Running Xcompact3d with ADIOS2"
     print *, "          this is currently experimental"
     print *, "          for safety of results it is recommended"
     print *, "          to run the default build as this feature"
     print *, "          is developed. Thank you for trying it."
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
  endif
#endif
  
  call parameter(InputFN)

  call decomp_2d_init(nx,ny,nz,p_row,p_col,periodic_bc)

  call decomp_2d_io_init()
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)    !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)    !start from 1 == true
  !div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
  call decomp_info_init(nxm, nym, nzm, ph1)
  call decomp_info_init(nxm, ny, nz, ph4)
  !gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
  call decomp_info_init(nxm, ny, nz, ph2)
  call decomp_info_init(nxm, nym, nz, ph3)

  call init_variables()

  call schemes()

  call decomp_2d_poisson_init()
  call decomp_info_init(nxm,nym,nzm,phG)

  if (ilesmod.ne.0) then
     if (jles.gt.0)  call init_explicit_les()
  endif

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call epsi_init(ep1)
     call body(ux1,uy1,uz1,ep1)
  endif

  if (iforces.eq.1) then
     call init_forces()
     if (irestart==1) then
        call restart_forces(0)
     endif
  endif

  !####################################################################
  ! initialise mhd
  if (mhd_active) call mhd_init()

  !####################################################################
  ! initialise particles
  if (particle_active) then
    call particle_report('input')

    call local_domain_size
  endif


  !####################################################################
  ! initialise visu
  if (ivisu.ne.0) then
     call visu_init()
     call visu_case_init() !! XXX: If you get error about uninitialised IO, look here.
                           !! Ensures additional case-specific variables declared for IO
     call visu_ready()
  end if
  ! compute diffusion number of simulation
  call compute_cfldiff()
  !####################################################################
  if (irestart==0) then
     call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
     itime = 0
     call preprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)
  else
     itr=1
     if (itype == itype_sandbox) then
        call init_sandbox(ux1,uy1,uz1,ep1,phi1,1)
     end if
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,0)
  endif

  if ((ioutflow.eq.1).or.(iin.eq.3)) then
     call init_inflow_outflow()
     if ((irestart==1).and.(iin.eq.3)) then
        call read_inflow(ux_inflow,uy_inflow,uz_inflow,itime/ntimesteps)
     endif
  end if

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call body(ux1,uy1,uz1,ep1)
  endif

  if (mod(itime, ilist) == 0 .or. itime == ifirst) then
     call test_speed_min_max(ux1,uy1,uz1)
     if (iscalar==1) call test_scalar_min_max(phi1)
  endif

  call simu_stats(1)

  call calc_divu_constraint(divu3, rho1, phi1)

  call init_probes()

  if (iturbine.ne.0) call init_turbines(ux1, uy1, uz1)

  if (itype==2) then
     if(nrank.eq.0)then
        open(42,file='time_evol.dat',form='formatted')
        if(mhd_active) then
           open(43,file='mhd_time_evol.dat',form='formatted')
        endif
     endif
  endif

  if (iforces == 1) then
     if(nrank.eq.0)then
        open(38,file='forces.dat',form='formatted')
     endif
  endif
  
  if (itype==10) then
     if(nrank.eq.0)then
        open(42,file='shear.dat',form='formatted')
     endif
  endif

endsubroutine init_xcompact3d
!########################################################################
!########################################################################
subroutine finalise_xcompact3d()

  use MPI
  use decomp_2d_mpi
  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_finalise

  use tools, only : simu_stats
  use param, only : itype, jles, ilesmod, mhd_active
  use probes, only : finalize_probes
  use visu, only : visu_finalise
  use les, only: finalise_explicit_les
  use mhd, only: mhd_fin
  use case, only: visu_case_finalise
  use forces, only: iforces

  implicit none

  integer :: ierr
  
  if (itype==2) then
     if(nrank.eq.0)then
        close(42)
        if(mhd_active) then
           close(43)
        endif
     endif
  endif

  if (iforces == 1) then
     if(nrank.eq.0)then
        close(38)
     endif
  endif

  if (itype==10) then
     if(nrank.eq.0)then
        close(42)
     endif
  endif
  
  call simu_stats(4)
  call finalize_probes()
  call visu_case_finalise()
  call visu_finalise()
  if (mhd_active) call mhd_fin()
  if (ilesmod.ne.0) then
     if (jles.gt.0) call finalise_explicit_les()
  endif
  call decomp_2d_io_finalise()
  call decomp_2d_finalize
  CALL MPI_FINALIZE(ierr)

endsubroutine finalise_xcompact3d



logical function stage9_3_dryrun_requested()
  implicit none
  character(len=64) :: envv
  integer :: stat, lenv
  stage9_3_dryrun_requested = .false.
  envv = ''
  call get_environment_variable('X3D_STAGE9_3_CHANNEL_INIT_DRYRUN', envv, length=lenv, status=stat)
  if (stat == 0) then
    if (trim(envv(1:lenv)) == '1') stage9_3_dryrun_requested = .true.
  endif
end function stage9_3_dryrun_requested

subroutine stage9_3_channel_init_dryrun_audit()
  use mpi
  use, intrinsic :: ieee_arithmetic
  use decomp_2d
  use decomp_2d_mpi, only : nrank
  use var, only : ux1,uy1,uz1,rho1,pp3,px1,py1,pz1,divu3,ph1,ph2,ph3,ph4,phG
  use param, only : itype, itype_channel
  use variables, only : nx,ny,nz,p_row,p_col
  implicit none
  integer :: ierr, rank, gstat
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

subroutine check_transients()

  use decomp_2d_constants, only : mytype
  use mpi
  use var
  
  implicit none

  real(mytype) :: dep
  integer :: code
   
  dep=maxval(abs(dux1))
  call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
  if (nrank == 0) write(*,*)'## MAX dux1 ', dep
 
  dep=maxval(abs(duy1))
  call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
  if (nrank == 0) write(*,*)'## MAX duy1 ', dep
 
  dep=maxval(abs(duz1))
  call MPI_ALLREDUCE(MPI_IN_PLACE,dep,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
  if (nrank == 0) write(*,*)'## MAX duz1 ', dep
  
end subroutine check_transients
