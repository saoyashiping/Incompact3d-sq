module fibre_stage8_boundary_safe_workflow
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata, only: stage7_channel_grid_t
  use fibre_stage7_boundary_safety, only: stage7_boundary_safety_result_t, classify_stage7_boundary_point
  use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t
  use fibre_stage8_lagrangian_state, only: stage8_lagrangian_state_t
  use fibre_stage8_velocity_to_fibre, only: clear_stage8_fluid_velocity_lag, interpolate_stage8_velocity_to_state
  use fibre_stage8_feedback_candidate, only: clear_stage8_slip_and_feedback, assemble_stage8_slip_velocity, assemble_stage8_linear_feedback_candidate
  use fibre_stage8_twoway_force_density, only: build_stage8_twoway_force_density_candidate
  implicit none
  type stage8_boundary_workflow_status_t
    integer :: boundary_classification_called_flag=0,velocity_interpolation_called_flag=0,slip_assembly_called_flag=0,feedback_candidate_called_flag=0,force_density_candidate_called_flag=0
    integer :: safe_count=0,blocked_count=0,unsafe_count=0,valid_count=0,rejected_flag=0
    integer :: velocity_write_allowed_flag=0,slip_write_allowed_flag=0,feedback_write_allowed_flag=0,force_density_write_allowed_flag=0
    integer :: rhs_hook_called_flag=0,rhs_modified_flag=0,pressure_poisson_called_flag=0,projection_called_flag=0,production_dns_called_flag=0,fluid_update_called_flag=0,fibre_advance_called_flag=0
    integer :: workflow_status=0
  end type
contains
  subroutine init_stage8_boundary_workflow_status(st)
    type(stage8_boundary_workflow_status_t),intent(out)::st; st=stage8_boundary_workflow_status_t()
  end
  subroutine apply_stage8_boundary_safe_coupling_workflow(grid,layout,ux,uy,uz,beta,state,fx,fy,fz,st,valid,rejected)
    type(stage7_channel_grid_t),intent(in)::grid; type(stage7_velocity_layout_t),intent(in)::layout
    real(mytype),intent(in)::ux(:,:,:),uy(:,:,:),uz(:,:,:),beta; type(stage8_lagrangian_state_t),intent(inout)::state
    real(mytype),intent(inout)::fx(:,:,:),fy(:,:,:),fz(:,:,:); type(stage8_boundary_workflow_status_t),intent(inout)::st
    integer,intent(out)::valid,rejected; type(stage7_boundary_safety_result_t)::res; integer::l,vc,bc,uc,fv,fr,fvc,fbc,fuc
    call init_stage8_boundary_workflow_status(st); fx=0; fy=0; fz=0; valid=0; rejected=1
    if(state%allocated_flag/=1) then; st%rejected_flag=1; return; endif
    st%boundary_classification_called_flag=1
    do l=1,state%nlag
      call classify_stage7_boundary_point(grid,layout,state%x(1,l),state%x(2,l),state%x(3,l),1,res)
      state%point_valid_flag(l)=res%safe_flag; state%point_blocked_flag(l)=res%blocked_flag; state%point_unsafe_flag(l)=res%unsafe_flag; state%point_status_code(l)=res%status_code
    end do
    st%safe_count=sum(state%point_valid_flag); st%blocked_count=sum(state%point_blocked_flag); st%unsafe_count=sum(state%point_unsafe_flag)
    call clear_stage8_fluid_velocity_lag(state); call clear_stage8_slip_and_feedback(state)
    call interpolate_stage8_velocity_to_state(grid,layout,ux,uy,uz,state,vc,bc,uc); st%velocity_interpolation_called_flag=1
    call assemble_stage8_slip_velocity(state,vc,bc,uc); st%slip_assembly_called_flag=1
    call assemble_stage8_linear_feedback_candidate(state,beta,fv,fr,fvc,fbc,fuc); st%feedback_candidate_called_flag=1
    call build_stage8_twoway_force_density_candidate(grid,layout,state,fx,fy,fz,fv,fr,fvc,fbc,fuc); st%force_density_candidate_called_flag=1
    st%valid_count=fvc; st%rejected_flag=fr
    st%velocity_write_allowed_flag=merge(1,0,max_blocked_abs(state%u_fluid_lag,state)<=1e-14_mytype)
    st%slip_write_allowed_flag=merge(1,0,max_blocked_abs(state%slip,state)<=1e-14_mytype)
    st%feedback_write_allowed_flag=merge(1,0,max(max_blocked_abs(state%force_structure,state),max_blocked_abs(state%force_fluid,state))<=1e-14_mytype)
    st%force_density_write_allowed_flag=merge(1,0,st%blocked_count==0 .or. maxval(abs(fx))+maxval(abs(fy))+maxval(abs(fz))<=1e-14_mytype)
    if(st%safe_count>0.and.st%blocked_count==0.and.st%unsafe_count==0.and.fv==1.and.fr==0) then; valid=1; rejected=0; endif
    if(st%blocked_count>0.or.st%unsafe_count>0) then; valid=0; rejected=1; endif
    st%workflow_status=merge(1,0,st%boundary_classification_called_flag==1.and.st%velocity_interpolation_called_flag==1.and.st%slip_assembly_called_flag==1.and.st%feedback_candidate_called_flag==1.and.st%force_density_candidate_called_flag==1)
  contains
    function max_blocked_abs(a,state) result(err)
      real(mytype),intent(in)::a(:,:); type(stage8_lagrangian_state_t),intent(in)::state; real(mytype)::err; integer::ll
      err=0; do ll=1,state%nlag; if(state%point_blocked_flag(ll)==1.or.state%point_unsafe_flag(ll)==1) err=max(err,maxval(abs(a(:,ll)))); end do
    end function
  end
  subroutine compute_stage8_blocked_state_write_error(state,velocity_err,slip_err,force_err)
    type(stage8_lagrangian_state_t),intent(in)::state; real(mytype),intent(out)::velocity_err,slip_err,force_err; integer::l
    velocity_err=0; slip_err=0; force_err=0
    do l=1,state%nlag
      if(state%point_blocked_flag(l)==1.or.state%point_unsafe_flag(l)==1) then
        velocity_err=max(velocity_err,maxval(abs(state%u_fluid_lag(:,l)))); slip_err=max(slip_err,maxval(abs(state%slip(:,l))))
        force_err=max(force_err,max(maxval(abs(state%force_structure(:,l))),maxval(abs(state%force_fluid(:,l)))))
      end if
    end do
  end
  subroutine compute_stage8_force_density_buffer_norm(fx,fy,fz,norm_max)
    real(mytype),intent(in)::fx(:,:,:),fy(:,:,:),fz(:,:,:); real(mytype),intent(out)::norm_max
    norm_max=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))
  end
end module
