module fibre_stage8_velocity_to_fibre
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata, only: stage7_channel_grid_t
  use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t, interpolate_stage7_velocity
  use fibre_stage7_boundary_safety, only: stage7_boundary_safety_result_t, classify_stage7_boundary_point
  use fibre_stage8_lagrangian_state, only: stage8_lagrangian_state_t
  implicit none
  private
  public :: interpolate_stage8_velocity_to_state
  public :: clear_stage8_fluid_velocity_lag
  public :: compute_stage8_velocity_state_error
contains
subroutine interpolate_stage8_velocity_to_state(grid, layout, ux, uy, uz, state, valid_count, blocked_count, unsafe_count)
  type(stage7_channel_grid_t), intent(in) :: grid
  type(stage7_velocity_layout_t), intent(in) :: layout
  real(mytype), intent(in) :: ux(:,:,:),uy(:,:,:),uz(:,:,:)
  type(stage8_lagrangian_state_t), intent(inout) :: state
  integer, intent(out) :: valid_count, blocked_count, unsafe_count
  type(stage7_boundary_safety_result_t) :: res
  real(mytype) :: ulag(3)
  integer :: l,v,u
  valid_count=0; blocked_count=0; unsafe_count=0
  if(state%allocated_flag/=1 .or. state%nlag<2) return
  do l=1,state%nlag
    call classify_stage7_boundary_point(grid,layout,state%x(1,l),state%x(2,l),state%x(3,l),1,res)
    state%point_valid_flag(l)=res%safe_flag
    state%point_blocked_flag(l)=res%blocked_flag
    state%point_unsafe_flag(l)=res%unsafe_flag
    state%point_status_code(l)=res%status_code
    if(res%safe_flag==1) then
      call interpolate_stage7_velocity(grid,layout,ux,uy,uz,state%x(1,l),state%x(2,l),state%x(3,l),ulag,v,u)
      if(v==1 .and. u==0) then
        state%u_fluid_lag(:,l)=ulag; valid_count=valid_count+1
      else
        state%u_fluid_lag(:,l)=0; unsafe_count=unsafe_count+1
      end if
    else
      state%u_fluid_lag(:,l)=0
      if(res%blocked_flag==1) blocked_count=blocked_count+1
      unsafe_count=unsafe_count+max(0,res%unsafe_flag)
    end if
  end do
end subroutine
subroutine clear_stage8_fluid_velocity_lag(state)
  type(stage8_lagrangian_state_t), intent(inout) :: state
  if(state%allocated_flag==1 .and. allocated(state%u_fluid_lag)) state%u_fluid_lag=0
end subroutine
subroutine compute_stage8_velocity_state_error(state, expected, err_max)
  type(stage8_lagrangian_state_t), intent(in) :: state
  real(mytype), intent(in) :: expected(:,:)
  real(mytype), intent(out) :: err_max
  err_max=huge(1._mytype)
  if(state%allocated_flag/=1) return
  err_max=maxval(abs(state%u_fluid_lag-expected))
end subroutine
end module
