module fibre_stage8_rk_sync
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata, only: stage7_channel_grid_t
  use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t
  use fibre_stage8_lagrangian_state, only: stage8_lagrangian_state_t
  use fibre_stage8_velocity_to_fibre, only: clear_stage8_fluid_velocity_lag, interpolate_stage8_velocity_to_state
  use fibre_stage8_feedback_candidate, only: clear_stage8_slip_and_feedback, assemble_stage8_slip_velocity, assemble_stage8_linear_feedback_candidate
  use fibre_stage8_twoway_force_density, only: build_stage8_twoway_force_density_candidate
  implicit none
  type stage8_rk_sync_status_t
    integer :: nsubstep=0,velocity_interpolation_count=0,slip_assembly_count=0,feedback_candidate_count=0,force_density_candidate_count=0,force_buffer_clear_count=0
    integer :: stale_velocity_detected_flag=0,stale_slip_detected_flag=0,stale_feedback_detected_flag=0,stale_force_density_detected_flag=0
    integer :: rhs_hook_called_flag=0,rhs_modified_flag=0,pressure_poisson_called_flag=0,projection_called_flag=0,production_dns_called_flag=0,fluid_update_called_flag=0,fibre_advance_called_flag=0
    integer :: clear_buffer_event_index=0,velocity_event_index=0,slip_event_index=0,feedback_event_index=0,force_density_event_index=0,hypothetical_rhs_gate_event_index=0,hypothetical_projection_gate_event_index=0
    integer :: event_order_status=0,rk_sync_status=0
    integer :: last_valid_count=0,last_blocked_count=0,last_unsafe_count=0
    integer :: last_valid_flag=0,last_rejected_flag=0
  end type
contains
  subroutine init_stage8_rk_sync_status(st,nsubstep)
    type(stage8_rk_sync_status_t),intent(out)::st; integer,intent(in)::nsubstep; st=stage8_rk_sync_status_t(); st%nsubstep=nsubstep
  end
  subroutine clear_stage8_force_density_buffer(fx,fy,fz)
    real(mytype),intent(inout)::fx(:,:,:),fy(:,:,:),fz(:,:,:); fx=0; fy=0; fz=0
  end
  subroutine compute_stage8_force_density_signature(grid,fx,fy,fz,sig)
    type(stage7_channel_grid_t),intent(in)::grid; real(mytype),intent(in)::fx(:,:,:),fy(:,:,:),fz(:,:,:); real(mytype),intent(out)::sig
    integer::i,j,k; sig=0
    do k=1,size(fx,3); do j=1,size(fx,2); do i=1,size(fx,1)
      sig=sig+(abs(fx(i,j,k))+2*abs(fy(i,j,k))+3*abs(fz(i,j,k)))*grid%volume_y(j)
    end do; end do; end do
  end
  subroutine apply_stage8_controlled_rk_substep_coupling(grid,layout,ux,uy,uz,beta,state,fx,fy,fz,st,substep_id,valid,rejected)
    type(stage7_channel_grid_t),intent(in)::grid; type(stage7_velocity_layout_t),intent(in)::layout
    real(mytype),intent(in)::ux(:,:,:),uy(:,:,:),uz(:,:,:),beta
    type(stage8_lagrangian_state_t),intent(inout)::state; real(mytype),intent(inout)::fx(:,:,:),fy(:,:,:),fz(:,:,:)
    type(stage8_rk_sync_status_t),intent(inout)::st; integer,intent(in)::substep_id; integer,intent(out)::valid,rejected
    integer::vc,bc,uc,fv,fr,fvc,fbc,fuc,e
    e=0; valid=0; rejected=1
    call clear_stage8_force_density_buffer(fx,fy,fz); st%force_buffer_clear_count=st%force_buffer_clear_count+1; e=e+1; st%clear_buffer_event_index=e
    call clear_stage8_fluid_velocity_lag(state); call clear_stage8_slip_and_feedback(state)
    call interpolate_stage8_velocity_to_state(grid,layout,ux,uy,uz,state,vc,bc,uc); st%velocity_interpolation_count=st%velocity_interpolation_count+1; e=e+1; st%velocity_event_index=e
    call assemble_stage8_slip_velocity(state,vc,bc,uc); st%slip_assembly_count=st%slip_assembly_count+1; e=e+1; st%slip_event_index=e
    call assemble_stage8_linear_feedback_candidate(state,beta,fv,fr,fvc,fbc,fuc); st%feedback_candidate_count=st%feedback_candidate_count+1; e=e+1; st%feedback_event_index=e
    call build_stage8_twoway_force_density_candidate(grid,layout,state,fx,fy,fz,fv,fr,fvc,fbc,fuc); st%force_density_candidate_count=st%force_density_candidate_count+1; e=e+1; st%force_density_event_index=e
    st%last_valid_count=fvc; st%last_blocked_count=fbc; st%last_unsafe_count=fuc
    st%last_valid_flag=fv; st%last_rejected_flag=fr
    e=e+1; st%hypothetical_rhs_gate_event_index=e; e=e+1; st%hypothetical_projection_gate_event_index=e
    st%event_order_status=merge(1,0,st%clear_buffer_event_index<st%velocity_event_index.and.st%velocity_event_index<st%slip_event_index.and.st%slip_event_index<st%feedback_event_index.and.st%feedback_event_index<st%force_density_event_index.and.st%force_density_event_index<st%hypothetical_rhs_gate_event_index.and.st%hypothetical_rhs_gate_event_index<st%hypothetical_projection_gate_event_index)
    if(fv==1.and.fr==0.and.fvc>0.and.fbc==0.and.fuc==0.and.st%rhs_hook_called_flag==0.and.st%projection_called_flag==0.and.st%production_dns_called_flag==0) then; valid=1; rejected=0; endif
  end
end module
