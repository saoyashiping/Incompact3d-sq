module fibre_stage8_oneway_forcing
  use, intrinsic :: ieee_arithmetic
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata, only: stage7_channel_grid_t
  use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t
  use fibre_stage8_lagrangian_state, only: stage8_lagrangian_state_t
  use fibre_stage8_velocity_to_fibre, only: clear_stage8_fluid_velocity_lag,interpolate_stage8_velocity_to_state
  use fibre_stage8_feedback_candidate, only: clear_stage8_slip_and_feedback,assemble_stage8_slip_velocity,assemble_stage8_linear_feedback_candidate
  implicit none
  private
  public :: apply_stage8_oneway_fluid_to_fibre_forcing
  public :: compute_stage8_oneway_force_error
  public :: compute_stage8_oneway_structure_power
contains
  subroutine apply_stage8_oneway_fluid_to_fibre_forcing(grid, layout, ux, uy, uz, beta, state, valid, rejected, valid_count, blocked_count, unsafe_count)
    type(stage7_channel_grid_t), intent(in) :: grid
    type(stage7_velocity_layout_t), intent(in) :: layout
    real(mytype), intent(in) :: ux(:,:,:),uy(:,:,:),uz(:,:,:),beta
    type(stage8_lagrangian_state_t), intent(inout) :: state
    integer, intent(out) :: valid,rejected,valid_count,blocked_count,unsafe_count
    integer :: vc,bc,uc,fv,fr,fvc,fbc,fuc
    valid=0; rejected=1; valid_count=0; blocked_count=0; unsafe_count=0
    if(state%allocated_flag/=1) return
    if(beta<=0._mytype .or. .not.ieee_is_finite(beta)) return
    call clear_stage8_fluid_velocity_lag(state)
    call clear_stage8_slip_and_feedback(state)
    call interpolate_stage8_velocity_to_state(grid,layout,ux,uy,uz,state,vc,bc,uc)
    call assemble_stage8_slip_velocity(state,vc,bc,uc)
    call assemble_stage8_linear_feedback_candidate(state,beta,fv,fr,fvc,fbc,fuc)
    valid_count=fvc; blocked_count=fbc; unsafe_count=fuc
    if(fv==1 .and. fr==0 .and. valid_count>0 .and. blocked_count==0 .and. unsafe_count==0) then
      valid=1; rejected=0
    else
      valid=0; rejected=1
    end if
  end subroutine

  subroutine compute_stage8_oneway_force_error(state, expected_force, err_max)
    type(stage8_lagrangian_state_t), intent(in) :: state
    real(mytype), intent(in) :: expected_force(:,:)
    real(mytype), intent(out) :: err_max
    err_max=huge(1._mytype)
    if(state%allocated_flag/=1) return
    err_max=maxval(abs(state%force_structure-expected_force))
  end subroutine

  subroutine compute_stage8_oneway_structure_power(state,power)
    type(stage8_lagrangian_state_t), intent(in) :: state
    real(mytype), intent(out) :: power
    integer :: l
    power=0
    if(state%allocated_flag/=1) return
    do l=1,state%nlag
      power=power+dot_product(state%force_structure(:,l),state%v_fibre(:,l))*state%ds(l)
    end do
  end subroutine
end module
