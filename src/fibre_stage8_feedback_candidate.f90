module fibre_stage8_feedback_candidate
  use, intrinsic :: ieee_arithmetic
  use fibre_parameters, only: mytype
  use fibre_stage8_lagrangian_state, only: stage8_lagrangian_state_t
  implicit none
  private
  public :: clear_stage8_slip_and_feedback
  public :: assemble_stage8_slip_velocity
  public :: assemble_stage8_linear_feedback_candidate
  public :: compute_stage8_feedback_action_reaction_error
  public :: compute_stage8_feedback_sign_metrics
contains
  subroutine clear_stage8_slip_and_feedback(state)
    type(stage8_lagrangian_state_t), intent(inout) :: state
    if(state%allocated_flag==1) then
      state%slip=0; state%force_structure=0; state%force_fluid=0
    end if
  end subroutine

  subroutine assemble_stage8_slip_velocity(state, valid_count, blocked_count, unsafe_count)
    type(stage8_lagrangian_state_t), intent(inout) :: state
    integer, intent(out) :: valid_count, blocked_count, unsafe_count
    integer :: l
    valid_count=0; blocked_count=0; unsafe_count=0
    if(state%allocated_flag/=1) return
    do l=1,state%nlag
      if(state%point_valid_flag(l)==1 .and. state%point_blocked_flag(l)==0 .and. state%point_unsafe_flag(l)==0) then
        state%slip(:,l)=state%u_fluid_lag(:,l)-state%v_fibre(:,l); valid_count=valid_count+1
      else
        state%slip(:,l)=0
        if(state%point_blocked_flag(l)==1) blocked_count=blocked_count+1
        if(state%point_unsafe_flag(l)==1) unsafe_count=unsafe_count+1
      end if
    end do
  end subroutine

  subroutine assemble_stage8_linear_feedback_candidate(state, beta, valid, rejected, valid_count, blocked_count, unsafe_count)
    type(stage8_lagrangian_state_t), intent(inout) :: state
    real(mytype), intent(in) :: beta
    integer, intent(out) :: valid,rejected,valid_count,blocked_count,unsafe_count
    integer :: l
    valid=0; rejected=1; valid_count=0; blocked_count=0; unsafe_count=0
    if(state%allocated_flag/=1) return
    if(beta<=0._mytype .or. .not.ieee_is_finite(beta)) return
    do l=1,state%nlag
      if(state%point_valid_flag(l)==1 .and. state%point_blocked_flag(l)==0 .and. state%point_unsafe_flag(l)==0) then
        state%force_structure(:,l)=beta*state%slip(:,l)
        state%force_fluid(:,l)=-state%force_structure(:,l)
        valid_count=valid_count+1
      else
        state%force_structure(:,l)=0; state%force_fluid(:,l)=0
        if(state%point_blocked_flag(l)==1) blocked_count=blocked_count+1
        if(state%point_unsafe_flag(l)==1) unsafe_count=unsafe_count+1
      end if
    end do
    valid=1; rejected=0
  end subroutine

  subroutine compute_stage8_feedback_action_reaction_error(state, err_max)
    type(stage8_lagrangian_state_t), intent(in) :: state
    real(mytype), intent(out) :: err_max
    integer :: l
    err_max=0
    if(state%allocated_flag/=1) return
    do l=1,state%nlag
      err_max=max(err_max,maxval(abs(state%force_structure(:,l)+state%force_fluid(:,l))))
    end do
  end subroutine

  subroutine compute_stage8_feedback_sign_metrics(state, structure_dot_slip, fluid_dot_slip, total_pair_power)
    type(stage8_lagrangian_state_t), intent(in) :: state
    real(mytype), intent(out) :: structure_dot_slip,fluid_dot_slip,total_pair_power
    integer :: l
    structure_dot_slip=0; fluid_dot_slip=0; total_pair_power=0
    if(state%allocated_flag/=1) return
    do l=1,state%nlag
      structure_dot_slip=structure_dot_slip+dot_product(state%force_structure(:,l),state%slip(:,l))*state%ds(l)
      fluid_dot_slip=fluid_dot_slip+dot_product(state%force_fluid(:,l),state%slip(:,l))*state%ds(l)
      total_pair_power=total_pair_power+(dot_product(state%force_structure(:,l),state%v_fibre(:,l))+dot_product(state%force_fluid(:,l),state%u_fluid_lag(:,l)))*state%ds(l)
    end do
  end subroutine
end module
