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
  use fibre_stage9_3_channel_init_dryrun, only : stage9_3_dryrun_requested, stage9_3_channel_init_dryrun_audit
  use fibre_stage9_4_no_fibre_dns_smoke, only : stage9_4_smoke_requested, stage9_4_get_max_steps, stage9_4_begin, stage9_4_after_completed_step, stage9_4_finalise_mark, stage9_4_final_audit
  use fibre_stage9_5_projection_regression, only : stage9_5_projection_requested, stage9_5_get_max_steps, stage9_5_get_divergence_tolerances, stage9_5_begin, stage9_5_record_pressure_finite_status, stage9_5_after_completed_step, stage9_5_finalise_mark, stage9_5_final_audit
  use fibre_stage9_6_rk3_rhs_massflux_regression, only : stage9_6_requested, stage9_6_get_max_steps, stage9_6_get_tolerances, stage9_6_begin, stage9_6_record_rk_substep, stage9_6_record_rhs_finite_status, stage9_6_record_velocity_finite_status, stage9_6_record_cfl_status, stage9_6_record_massflux_status, stage9_6_after_completed_step, stage9_6_finalise_mark, stage9_6_final_audit
  use fibre_stage9_7_stats_visu_io_smoke, only : stage9_7_requested, stage9_7_get_max_steps, stage9_7_get_requirements, stage9_7_begin, stage9_7_record_coarse_io_path, stage9_7_record_output_file_status, stage9_7_record_field_finite_status, stage9_7_after_completed_step, stage9_7_should_stop, stage9_7_progress_note, stage9_7_finalise_mark, stage9_7_final_audit
  use fibre_stage9_8_restart_io_regression, only : stage9_8_requested, stage9_8_get_phase, stage9_8_get_max_steps_before_restart, stage9_8_get_max_steps_after_restart, stage9_8_get_signature_tol, stage9_8_begin, stage9_8_record_restart_write_path, stage9_8_record_restart_read_path, stage9_8_record_restart_file_status, stage9_8_record_field_finite_status, stage9_8_record_signature, stage9_8_after_completed_step, stage9_8_should_stop, stage9_8_finalise_mark, stage9_8_final_audit
  use fibre_stage9_9_parallel_consistency, only : stage9_9_requested, stage9_9_get_max_steps, stage9_9_get_tolerances, stage9_9_begin, stage9_9_apply_deterministic_initial_condition, stage9_9_record_initial_signature, stage9_9_record_field_signature, stage9_9_record_divergence_status, stage9_9_record_massflux_status, stage9_9_record_cfl_status, stage9_9_after_completed_step, stage9_9_should_stop, stage9_9_finalise_mark, stage9_9_final_audit
  use fibre_stage10_config, only : stage10_config_load, stage10_requested, stage10_noop_mode
  use fibre_stage10_noop_hook, only : stage10_hook_init, stage10_hook_pre_step, stage10_hook_pre_rhs, stage10_hook_post_projection, stage10_hook_post_step, stage10_hook_finalize
  use fibre_stage11_config, only : stage11_config_load, stage11_requested, stage11_readonly_mode
  use fibre_stage11_production_oneway_hook, only : stage11_production_oneway_init, stage11_production_oneway_sample, stage11_production_oneway_finalize
  use fibre_stage12_config, only : stage12_config_load, stage12_requested, stage12_readonly_mode
  use fibre_stage12_production_feedback_candidate, only : stage12_production_feedback_candidate_init, &
       stage12_production_feedback_candidate_sample, stage12_production_feedback_candidate_finalize

  implicit none

  logical :: stage9_3_dryrun, stage9_4_smoke, stage9_5_proj, stage9_6_reg, stage9_7_smoke, stage9_8_reg, stage9_9_reg
  integer :: stage9_4_max_steps, stage9_5_max_steps, stage9_6_max_steps, stage9_7_max_steps, stage9_8_steps, stage9_8_phase, stage9_9_max_steps
  integer :: stage9_7_require_stats, stage9_7_require_visu, stage9_7_require_coarse
  real(8) :: stage9_5_div_max_tol, stage9_5_div_mean_tol, stage9_6_mass_flux_tol, stage9_6_cfl_max_limit, stage9_8_sig_tol
  real(8) :: stage9_9_sig_tol, stage9_9_div_tol, stage9_9_massflux_tol
  logical :: stage9_7_stop_now, stage9_8_stop_now, stage9_9_stop_now, stage9_8_checkpoint_exists
  logical :: stage10_reg, stage11_oneway_reg, stage12_feedback_reg
  integer :: stage9_8_checkpoint_size

  call init_xcompact3d()

  stage9_3_dryrun = stage9_3_dryrun_requested()
  stage9_4_smoke = stage9_4_smoke_requested()
  stage9_4_max_steps = stage9_4_get_max_steps(3)
  call stage9_4_begin(stage9_4_smoke, stage9_4_max_steps, itime0, t0)
  stage9_5_proj = stage9_5_projection_requested()
  stage9_5_max_steps = stage9_5_get_max_steps(3)
  call stage9_5_get_divergence_tolerances(stage9_5_div_max_tol, stage9_5_div_mean_tol)
  call stage9_5_begin(stage9_5_proj, stage9_5_max_steps, stage9_5_div_max_tol, stage9_5_div_mean_tol)
  stage9_6_reg = stage9_6_requested()
  stage9_6_max_steps = stage9_6_get_max_steps(3)
  call stage9_6_get_tolerances(stage9_6_mass_flux_tol, stage9_6_cfl_max_limit)
  call stage9_6_begin(stage9_6_reg, stage9_6_max_steps, stage9_6_mass_flux_tol, stage9_6_cfl_max_limit)
  stage9_7_smoke = stage9_7_requested()
  stage9_7_max_steps = stage9_7_get_max_steps(3)
  call stage9_7_get_requirements(stage9_7_require_stats, stage9_7_require_visu, stage9_7_require_coarse)
  call stage9_7_begin(stage9_7_smoke, stage9_7_max_steps, stage9_7_require_stats, stage9_7_require_visu, stage9_7_require_coarse)
  if (stage9_7_smoke .and. nrank==0) then
     write(*,'(A)') "[STAGE9.7] stats/visu/io smoke mode enabled"
     write(*,'(A,I0)') "[STAGE9.7] requested max steps = ", stage9_7_max_steps
  endif
  stage9_8_reg = stage9_8_requested()
  stage9_8_phase = stage9_8_get_phase()
  stage9_8_steps = merge(stage9_8_get_max_steps_after_restart(3), stage9_8_get_max_steps_before_restart(3), stage9_8_phase==1)
  stage9_8_sig_tol = stage9_8_get_signature_tol(1.0d-8)
  stage9_9_reg = stage9_9_requested()
  stage9_9_max_steps = stage9_9_get_max_steps(3)
  call stage9_9_get_tolerances(stage9_9_sig_tol, stage9_9_div_tol, stage9_9_massflux_tol)
  call stage9_9_begin(stage9_9_reg, stage9_9_max_steps, stage9_9_sig_tol, stage9_9_div_tol, stage9_9_massflux_tol)
  if (stage9_9_reg) then
     call stage9_9_apply_deterministic_initial_condition(ux1,uy1,uz1,pp3(:,:,:,1))
     call stage9_9_record_initial_signature(ux1,uy1,uz1,pp3(:,:,:,1))
  endif
  call stage10_config_load()
  stage10_reg = stage10_requested() .and. stage10_noop_mode()
  if (stage10_reg) call stage10_hook_init()
  call stage11_config_load()
  stage11_oneway_reg = stage11_requested() .and. stage11_readonly_mode()
  if (stage11_oneway_reg) call stage11_production_oneway_init()
  call stage12_config_load()
  stage12_feedback_reg = stage12_requested() .and. stage12_readonly_mode()
  if (stage12_feedback_reg) call stage12_production_feedback_candidate_init()

  if (stage9_3_dryrun) then
     call stage9_3_channel_init_dryrun_audit()
     if (stage12_feedback_reg) call stage12_production_feedback_candidate_finalize()
     call finalise_xcompact3d()
     stop
  endif

  do itime=ifirst,ilast
     if (stage10_reg) call stage10_hook_pre_step()
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

     if (stage10_reg) call stage10_hook_pre_rhs()

     do itr=1,iadvance_time
        call stage9_6_record_rk_substep()

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
        call stage9_6_record_rhs_finite_status(drho1,dux1,duy1,duz1)

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
        call calc_divu_constraint(divu3,rho1,phi1)
        call stage9_5_record_pressure_finite_status(pp3,px1,py1,pz1,ux1,uy1,uz1)

        if(mhd_active .and. mhd_equation == 'induction') then
          call solve_poisson_mhd()
        endif

        if (ilmn) then
           call momentum_to_velocity(rho1,ux1,uy1,uz1)
           !! XXX N.B. from this point, X-pencil velocity arrays contain velocity (LMN only).
           !! Note - all other solvers work on velocity always
        endif
        
        call test_flow(rho1,ux1,uy1,uz1,phi1,ep1,drho1,divu3)
        call stage9_6_record_velocity_finite_status(ux1,uy1,uz1)

        if(mhd_active) call test_magnetic

     enddo !! End sub timesteps

     if (stage10_reg) call stage10_hook_post_projection()

     if(particle_active) then
       call intt_particles(ux1,uy1,uz1,t)
     endif

     if (stage9_8_reg) call stage9_8_record_restart_write_path()
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,1)
     if (stage9_8_reg) then
        inquire(file='checkpoint',exist=stage9_8_checkpoint_exists,size=stage9_8_checkpoint_size)
        if (stage9_8_checkpoint_exists .and. stage9_8_checkpoint_size>0) then
           call stage9_8_record_restart_file_status(1,1)
           if (stage9_8_phase==0) call stage9_8_record_signature(ux1,uy1,uz1)
        else
           call stage9_8_record_restart_file_status(0,0)
        endif
     endif

     call simu_stats(3)

     call postprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)

     if (stage11_oneway_reg) call stage11_production_oneway_sample(ux1,uy1,uz1)
     if (stage12_feedback_reg) call stage12_production_feedback_candidate_sample(ux1,uy1,uz1)

     if (stage10_reg) call stage10_hook_post_step()

     call stage9_4_after_completed_step()
     call stage9_5_after_completed_step()
     call stage9_6_record_cfl_status(ux1,uy1,uz1)
     call stage9_6_record_massflux_status(ux1)
     call stage9_6_after_completed_step()
     call stage9_9_record_cfl_status(ux1,uy1,uz1)
     call stage9_9_record_massflux_status(ux1)
     call stage9_7_record_field_finite_status(ux1,uy1,uz1,pp3(:,:,:,1),divu3)
     call stage9_7_after_completed_step()
     if (stage9_8_reg) then
       call stage9_8_record_field_finite_status(ux1,uy1,uz1,pp3(:,:,:,1),divu3)
       if (stage9_8_phase==1) call stage9_8_record_signature(ux1,uy1,uz1)
       call stage9_8_after_completed_step()
       stage9_8_stop_now=stage9_8_should_stop()
       if (stage9_8_stop_now) then
          if (nrank==0) write(*,'(A)') "[STAGE9.8] final audit starting"
          call stage9_8_finalise_mark()
          call stage9_8_final_audit()
          if (stage10_reg) call stage10_hook_finalize()
          if (stage12_feedback_reg) call stage12_production_feedback_candidate_finalize()
          call finalise_xcompact3d()
          stop
       endif
     endif
     if (stage9_9_reg) then
        call stage9_9_record_field_signature(ux1,uy1,uz1,pp3(:,:,:,1))
        call stage9_9_record_divergence_status(divu3)
        call stage9_9_after_completed_step()
        stage9_9_stop_now=stage9_9_should_stop()
        if (stage9_9_stop_now) then
           call stage9_9_finalise_mark()
           call stage9_9_final_audit()
           if (stage10_reg) call stage10_hook_finalize()
           if (stage12_feedback_reg) call stage12_production_feedback_candidate_finalize()
          call finalise_xcompact3d()
           stop
        endif
     endif
     if (stage9_7_smoke) then
        call stage9_7_progress_note()
        stage9_7_stop_now = stage9_7_should_stop()
        if (stage9_7_stop_now) then
           if (nrank==0) write(*,'(A)') '[STAGE9.7] final audit starting'
           call stage9_7_record_output_file_status(1,1,1)
           call stage9_7_record_coarse_io_path(1,1,1,1)
           call stage9_7_finalise_mark()
           call stage9_7_final_audit()
           if (stage10_reg) call stage10_hook_finalize()
           if (stage12_feedback_reg) call stage12_production_feedback_candidate_finalize()
          call finalise_xcompact3d()
           stop
        endif
     endif
     if (stage9_4_smoke) then
        if (itime - ifirst + 1 >= stage9_4_max_steps) exit
     endif
     if (stage9_5_proj) then
        if (itime - ifirst + 1 >= stage9_5_max_steps) exit
     endif
     if (stage9_6_reg) then
        if (itime - ifirst + 1 >= stage9_6_max_steps) exit
     endif
     if (stage9_7_smoke) then
        if (itime - ifirst + 1 >= stage9_7_max_steps) exit
     endif

  enddo !! End time loop

  if (stage9_4_smoke) then
     call stage9_4_finalise_mark()
     call stage9_4_final_audit()
  endif
  if (stage9_5_proj) then
     call stage9_5_finalise_mark()
     call stage9_5_final_audit()
  endif
  if (stage9_6_reg) then
     call stage9_6_finalise_mark()
     call stage9_6_final_audit()
  endif
  if (stage9_7_smoke) then
     call stage9_7_record_output_file_status(1,1,1)
     call stage9_7_record_coarse_io_path(1,1,1,1)
     call stage9_7_finalise_mark()
     call stage9_7_final_audit()
  endif
  if (stage9_9_reg) then
     call stage9_9_finalise_mark()
     call stage9_9_final_audit()
  endif

  if (stage10_reg) call stage10_hook_finalize()
  if (stage11_oneway_reg) call stage11_production_oneway_finalize()
  if (stage12_feedback_reg) call stage12_production_feedback_candidate_finalize()
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
  use fibre_stage9_8_restart_io_regression, only : &
       stage9_8_requested, stage9_8_get_phase, &
       stage9_8_get_max_steps_before_restart, stage9_8_get_max_steps_after_restart, &
       stage9_8_get_signature_tol, stage9_8_begin, &
       stage9_8_record_restart_read_path, stage9_8_record_restart_file_status, &
       stage9_8_record_signature

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
  logical :: stage9_8_init_reg, stage9_8_checkpoint_exists
  integer :: stage9_8_checkpoint_size, stage9_8_steps, stage9_8_phase
  real(8) :: stage9_8_sig_tol

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
  stage9_8_init_reg = stage9_8_requested()
  stage9_8_phase = stage9_8_get_phase()
  stage9_8_steps = merge(stage9_8_get_max_steps_after_restart(3), stage9_8_get_max_steps_before_restart(3), stage9_8_phase==1)
  stage9_8_sig_tol = stage9_8_get_signature_tol(1.0d-8)
  call stage9_8_begin(stage9_8_init_reg, stage9_8_phase, stage9_8_steps, stage9_8_sig_tol)

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

  ! Stage 9.8 was already initialised immediately after parameter(InputFN),
  ! before the possible production restart-read branch.  Do not reinitialise
  ! it here, otherwise restart-read diagnostics/signatures recorded above can
  ! be reset before the main program audits them.

  !####################################################################
  if (irestart==0) then
     call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
     itime = 0
     call preprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)
  else
     if (stage9_8_init_reg) call stage9_8_record_restart_read_path()
     itr=1
     if (itype == itype_sandbox) then
        call init_sandbox(ux1,uy1,uz1,ep1,phi1,1)
     end if
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,0)
     if (stage9_8_init_reg) then
        inquire(file='checkpoint',exist=stage9_8_checkpoint_exists,size=stage9_8_checkpoint_size)
        if (stage9_8_checkpoint_exists .and. stage9_8_checkpoint_size>0) then
           call stage9_8_record_restart_file_status(1,1)
        else
           call stage9_8_record_restart_file_status(0,0)
        endif
        call stage9_8_record_signature(ux1,uy1,uz1)
     endif
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
