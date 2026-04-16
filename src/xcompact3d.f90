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
  use fiber_types, only : fiber_active, interp_solver_test_active, interp_solver_output_step, &
       rigid_coupling_test_active, rigid_free_test_active, rigid_kinematics_test_active, rigid_two_way_test_active, &
       rigid_kinematics_standalone, rigid_two_way_subiterations
  use fiber_io, only : write_fiber_interp_solver
  use fiber_interp, only : run_fiber_interp_solver_readonly
  use fiber_coupling, only : run_rigid_coupling_step, rigid_two_way_prepare_forcing, rigid_two_way_finalize_update, &
       rigid_two_way_commit_state
  use fiber_rigid_free, only : run_rigid_free_step
  use fiber_rigid_kinematics, only : rigid_kinematics_step

  implicit none

  logical :: solver_interp_done, two_way_converged
  integer :: two_way_isubiter, two_way_nsubiter
  real(mytype), allocatable :: ux_base(:,:,:), uy_base(:,:,:), uz_base(:,:,:)
  real(mytype), allocatable :: ux_work(:,:,:), uy_work(:,:,:), uz_work(:,:,:)
  real(mytype), allocatable :: rho_base(:,:,:,:), phi_base(:,:,:,:)

  solver_interp_done = .false.

  call init_xcompact3d()

  do itime=ifirst,ilast
     !t=itime*dt
     t=t0 + (itime0 + itime + 1 - ifirst)*dt
     
     call simu_stats(2)

     if (fiber_active .and. rigid_kinematics_test_active .and. rigid_kinematics_standalone) then
        do itr=1,iadvance_time
           call rigid_kinematics_step(ux1, uy1, uz1, t, itime, itr)
        enddo
        call simu_stats(3)
        cycle
     endif

     if (iturbine.ne.0) call compute_turbines(ux1, uy1, uz1)

     if (iin.eq.3.and.mod(itime,ntimesteps)==1) then
        call read_inflow(ux_inflow,uy_inflow,uz_inflow,itime/ntimesteps)
     endif

     if ((itype.eq.itype_abl.or.iturbine.ne.0).and.(ifilter.ne.0).and.(ilesmod.ne.0)) then
        call filter(C_filter)
        call apply_spatial_filter(ux1,uy1,uz1,phi1)
     endif

     do itr=1,iadvance_time

        if (fiber_active .and. rigid_two_way_test_active) then
           two_way_nsubiter = max(1, rigid_two_way_subiterations)
           two_way_converged = .false.
           if (.not.allocated(ux_base)) then
              allocate(ux_base(size(ux1,1), size(ux1,2), size(ux1,3)))
              allocate(uy_base(size(uy1,1), size(uy1,2), size(uy1,3)))
              allocate(uz_base(size(uz1,1), size(uz1,2), size(uz1,3)))
              allocate(ux_work(size(ux1,1), size(ux1,2), size(ux1,3)))
              allocate(uy_work(size(uy1,1), size(uy1,2), size(uy1,3)))
              allocate(uz_work(size(uz1,1), size(uz1,2), size(uz1,3)))
              allocate(rho_base(size(rho1,1), size(rho1,2), size(rho1,3), size(rho1,4)))
              allocate(phi_base(size(phi1,1), size(phi1,2), size(phi1,3), size(phi1,4)))
           endif
           ux_base = ux1
           uy_base = uy1
           uz_base = uz1
           rho_base = rho1
           phi_base = phi1

           do two_way_isubiter = 1, two_way_nsubiter
              ux_work = ux_base
              uy_work = uy_base
              uz_work = uz_base
              rho1 = rho_base
              phi1 = phi_base
              call set_fluid_properties(rho1,mu1)
              call boundary_conditions(rho1,ux_work,uy_work,uz_work,phi1,ep1)

              if (imove.eq.1) then ! update epsi for moving objects
                if ((iibm.eq.2).or.(iibm.eq.3)) then
                   call genepsi3d(ep1)
                else if (iibm.eq.1) then
                   call body(ux_work,uy_work,uz_work,ep1)
                endif
              endif

              call rigid_two_way_prepare_forcing(ux_work, uy_work, uz_work, t, itime, itr, iadvance_time, &
                   two_way_isubiter, two_way_nsubiter)

              call calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux_work,uy_work,uz_work,ep1,phi1,divu3)

#ifdef DEBG
              call check_transients()
#endif

              if (ilmn) then
                 !! XXX N.B. from this point, X-pencil velocity arrays contain momentum (LMN only).
                 call velocity_to_momentum(rho1,ux_work,uy_work,uz_work)
              endif

              call int_time(rho1,ux_work,uy_work,uz_work,phi1,drho1,dux1,duy1,duz1,dphi1)
              call pre_correc(ux_work,uy_work,uz_work,ep1)

              call calc_divu_constraint(divu3,rho1,phi1)
              call solve_poisson(pp3,px1,py1,pz1,rho1,ux_work,uy_work,uz_work,ep1,drho1,divu3)
              call cor_vel(ux_work,uy_work,uz_work,px1,py1,pz1)

              if(mhd_active .and. mhd_equation == 'induction') then
                call solve_poisson_mhd()
              endif

              if (ilmn) then
                 call momentum_to_velocity(rho1,ux_work,uy_work,uz_work)
                 !! XXX N.B. from this point, X-pencil velocity arrays contain velocity (LMN only).
                 !! Note - all other solvers work on velocity always
              endif

              call rigid_two_way_finalize_update(ux_work, uy_work, uz_work, t, itime, itr, iadvance_time, &
                   two_way_isubiter, two_way_nsubiter, two_way_converged)

              call test_flow(rho1,ux_work,uy_work,uz_work,phi1,ep1,drho1,divu3)

              if(mhd_active) call test_magnetic

              if (two_way_converged) exit
           enddo
           call rigid_two_way_commit_state(ux1, uy1, uz1, ux_work, uy_work, uz_work, t, itime, itr, iadvance_time, &
                min(two_way_isubiter, two_way_nsubiter), two_way_nsubiter)
        else
           call set_fluid_properties(rho1,mu1)
           call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)

           if (imove.eq.1) then ! update epsi for moving objects
             if ((iibm.eq.2).or.(iibm.eq.3)) then
                call genepsi3d(ep1)
             else if (iibm.eq.1) then
                call body(ux1,uy1,uz1,ep1)
             endif
           endif

           if (fiber_active .and. rigid_coupling_test_active) then
              call run_rigid_coupling_step(ux1, uy1, uz1, t, itime, itr, iadvance_time)
           else if (fiber_active .and. rigid_free_test_active) then
              call run_rigid_free_step(ux1, uy1, uz1, t, itime, itr, iadvance_time)
           else if (fiber_active .and. rigid_kinematics_test_active) then
              call rigid_kinematics_step(ux1, uy1, uz1, t, itime, itr)
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
        endif

     enddo !! End sub timesteps

     if(particle_active) then
       call intt_particles(ux1,uy1,uz1,t)
     endif

     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,1)

     call simu_stats(3)

     call postprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)

     if (fiber_active .and. interp_solver_test_active .and. .not.solver_interp_done) then
        if (itime >= interp_solver_output_step) then
           call run_fiber_interp_solver_readonly(ux1, uy1, uz1, itime)
           call write_fiber_interp_solver(itime)
           solver_interp_done = .true.
        endif
     endif

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
  use fiber_types, only : fiber_active, interp_test_active, interp_solver_test_active, &
       spread_test_active, rigid_coupling_test_active, rigid_free_test_active, rigid_kinematics_test_active, &
       rigid_two_way_test_active, &
       rigid_kinematics_standalone
  use fiber_init, only : init_fiber
  use fiber_rigid_motion, only : init_rigid_motion_reference
  use fiber_rigid_free, only : init_rigid_free_state
  use fiber_io, only : write_fiber_points, write_fiber_interp, &
       write_fiber_spread_lagrangian, write_fiber_spread_summary
  use fiber_interp, only : run_fiber_interp_operator_test
  use fiber_spread, only : run_fiber_spread_conservation_test

  implicit none

  integer :: ierr
  real(mytype) :: lag_total(3), eul_total(3), abs_error(3), rel_error(3)
  real(mytype) :: spread_sumw_min, spread_sumw_max

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

  if (interp_test_active .and. interp_solver_test_active) then
     if (nrank == 0) write(*,*) "Error: interp_test_active and interp_solver_test_active cannot both be true."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_coupling_test_active .and. (interp_test_active .or. interp_solver_test_active .or. spread_test_active)) then
     if (nrank == 0) write(*,*) "Error: rigid_coupling_test_active cannot be combined with other fiber test modes."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_free_test_active .and. (interp_test_active .or. interp_solver_test_active .or. spread_test_active)) then
     if (nrank == 0) write(*,*) "Error: rigid_free_test_active cannot be combined with other fiber test modes."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_kinematics_test_active .and. (interp_test_active .or. interp_solver_test_active .or. spread_test_active)) then
     if (nrank == 0) write(*,*) "Error: rigid_kinematics_test_active cannot be combined with interpolation or spread test modes."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_two_way_test_active .and. (interp_test_active .or. interp_solver_test_active .or. spread_test_active)) then
     if (nrank == 0) write(*,*) "Error: rigid_two_way_test_active cannot be combined with interpolation or spread test modes."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_coupling_test_active .and. rigid_free_test_active) then
     if (nrank == 0) write(*,*) "Error: rigid_coupling_test_active and rigid_free_test_active cannot both be true."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_two_way_test_active .and. &
       (rigid_coupling_test_active .or. rigid_free_test_active .or. rigid_kinematics_test_active)) then
     if (nrank == 0) write(*,*) "Error: rigid_two_way_test_active cannot be combined with other rigid fiber test modes."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_kinematics_test_active .and. (rigid_coupling_test_active .or. rigid_free_test_active)) then
     if (nrank == 0) then
        write(*,*) "Error: rigid_kinematics_test_active cannot be combined with rigid_coupling_test_active."
     endif
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_coupling_test_active .and. .not.fiber_active) then
     if (nrank == 0) write(*,*) "Error: rigid_coupling_test_active requires fiber_active = true."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_free_test_active .and. .not.fiber_active) then
     if (nrank == 0) write(*,*) "Error: rigid_free_test_active requires fiber_active = true."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_kinematics_test_active .and. .not.fiber_active) then
     if (nrank == 0) write(*,*) "Error: rigid_kinematics_test_active requires fiber_active = true."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_two_way_test_active .and. .not.fiber_active) then
     if (nrank == 0) write(*,*) "Error: rigid_two_way_test_active requires fiber_active = true."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_kinematics_standalone .and. .not.rigid_kinematics_test_active) then
     if (nrank == 0) write(*,*) "Error: rigid_kinematics_standalone requires rigid_kinematics_test_active = true."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (rigid_kinematics_standalone .and. (interp_test_active .or. interp_solver_test_active .or. spread_test_active .or. &
       rigid_coupling_test_active .or. rigid_free_test_active .or. rigid_two_way_test_active)) then
     if (nrank == 0) write(*,*) "Error: rigid_kinematics_standalone cannot be combined with other fiber test modes."
     call MPI_FINALIZE(ierr)
     stop
  endif


  if (spread_test_active .and. (interp_test_active .or. interp_solver_test_active)) then
     if (nrank == 0) write(*,*) "Error: spread_test_active cannot be combined with interpolation test modes."
     call MPI_FINALIZE(ierr)
     stop
  endif

  if (fiber_active) then
     call init_fiber()
     if (rigid_coupling_test_active) call init_rigid_motion_reference()
     if (rigid_free_test_active .or. rigid_two_way_test_active) call init_rigid_free_state()
     call write_fiber_points()

     if (interp_test_active) then
        call run_fiber_interp_operator_test()
        call write_fiber_interp()
        if (nrank == 0) write(*,*) "Interpolation operator test complete. Exiting before solver initialization."
        call MPI_FINALIZE(ierr)
        stop
     endif

     if (spread_test_active) then
        call run_fiber_spread_conservation_test(lag_total, eul_total, abs_error, rel_error, spread_sumw_min, spread_sumw_max)
        call write_fiber_spread_lagrangian()
        call write_fiber_spread_summary(lag_total, eul_total, abs_error, rel_error, spread_sumw_min, spread_sumw_max)
        if (nrank == 0) write(*,*) "Spreading conservation test complete. Exiting before solver initialization."
        call MPI_FINALIZE(ierr)
        stop
     endif
  endif

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
  real(mytype) :: lag_total(3), eul_total(3), abs_error(3), rel_error(3)
  real(mytype) :: spread_sumw_min, spread_sumw_max
  
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
