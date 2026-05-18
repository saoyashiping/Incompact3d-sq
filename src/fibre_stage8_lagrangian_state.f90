module fibre_stage8_lagrangian_state
  use, intrinsic :: ieee_arithmetic
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata, only: stage7_channel_grid_t
  use fibre_stage7_boundary_safety, only: stage7_boundary_safety_result_t, classify_stage7_boundary_point
  use fibre_stage7_velocity_interpolation, only: stage7_velocity_layout_t
  implicit none
  private
  public :: stage8_lagrangian_state_t
  public :: init_stage8_lagrangian_state, allocate_stage8_lagrangian_state
  public :: clear_stage8_lagrangian_state, destroy_stage8_lagrangian_state
  public :: build_stage8_straight_fibre_state, validate_stage8_lagrangian_state
  public :: compute_stage8_lagrangian_total_length, compute_stage8_lagrangian_segment_error

  type stage8_lagrangian_state_t
    integer :: nlag
    integer :: allocated_flag
    integer :: valid_flag
    integer :: rejected_flag
    real(mytype), allocatable :: x(:,:)
    real(mytype), allocatable :: v_fibre(:,:)
    real(mytype), allocatable :: u_fluid_lag(:,:)
    real(mytype), allocatable :: slip(:,:)
    real(mytype), allocatable :: force_structure(:,:)
    real(mytype), allocatable :: force_fluid(:,:)
    real(mytype), allocatable :: ds(:)
    integer, allocatable :: point_valid_flag(:)
    integer, allocatable :: point_blocked_flag(:)
    integer, allocatable :: point_unsafe_flag(:)
    integer, allocatable :: point_status_code(:)
  end type
contains
subroutine init_stage8_lagrangian_state(state)
  type(stage8_lagrangian_state_t), intent(inout) :: state
  state%nlag=0; state%allocated_flag=0; state%valid_flag=0; state%rejected_flag=0
end subroutine
subroutine allocate_stage8_lagrangian_state(state,nlag,valid,rejected)
  type(stage8_lagrangian_state_t), intent(inout) :: state
  integer, intent(in) :: nlag
  integer, intent(out) :: valid,rejected
  call destroy_stage8_lagrangian_state(state)
  if(nlag<2) then; valid=0; rejected=1; state%valid_flag=0; state%rejected_flag=1; return; end if
  allocate(state%x(3,nlag),state%v_fibre(3,nlag),state%u_fluid_lag(3,nlag),state%slip(3,nlag))
  allocate(state%force_structure(3,nlag),state%force_fluid(3,nlag),state%ds(nlag))
  allocate(state%point_valid_flag(nlag),state%point_blocked_flag(nlag),state%point_unsafe_flag(nlag),state%point_status_code(nlag))
  state%x=0; state%v_fibre=0; state%u_fluid_lag=0; state%slip=0; state%force_structure=0; state%force_fluid=0; state%ds=0
  state%point_valid_flag=0; state%point_blocked_flag=0; state%point_unsafe_flag=0; state%point_status_code=0
  state%nlag=nlag; state%allocated_flag=1; state%valid_flag=1; state%rejected_flag=0; valid=1; rejected=0
end subroutine
subroutine clear_stage8_lagrangian_state(state)
  type(stage8_lagrangian_state_t), intent(inout) :: state
  if(state%allocated_flag==1) then
    state%x=0; state%v_fibre=0; state%u_fluid_lag=0; state%slip=0; state%force_structure=0; state%force_fluid=0; state%ds=0
    state%point_valid_flag=0; state%point_blocked_flag=0; state%point_unsafe_flag=0; state%point_status_code=0
    state%valid_flag=0; state%rejected_flag=0
  end if
end subroutine
subroutine destroy_stage8_lagrangian_state(state)
  type(stage8_lagrangian_state_t), intent(inout) :: state
  if(allocated(state%x)) deallocate(state%x)
  if(allocated(state%v_fibre)) deallocate(state%v_fibre)
  if(allocated(state%u_fluid_lag)) deallocate(state%u_fluid_lag)
  if(allocated(state%slip)) deallocate(state%slip)
  if(allocated(state%force_structure)) deallocate(state%force_structure)
  if(allocated(state%force_fluid)) deallocate(state%force_fluid)
  if(allocated(state%ds)) deallocate(state%ds)
  if(allocated(state%point_valid_flag)) deallocate(state%point_valid_flag)
  if(allocated(state%point_blocked_flag)) deallocate(state%point_blocked_flag)
  if(allocated(state%point_unsafe_flag)) deallocate(state%point_unsafe_flag)
  if(allocated(state%point_status_code)) deallocate(state%point_status_code)
  call init_stage8_lagrangian_state(state)
end subroutine
subroutine build_stage8_straight_fibre_state(state,grid,x0,y0,z0,length,direction,valid,rejected)
  type(stage8_lagrangian_state_t), intent(inout) :: state
  type(stage7_channel_grid_t), intent(in) :: grid
  real(mytype), intent(in) :: x0,y0,z0,length,direction(3)
  integer, intent(out) :: valid,rejected
  real(mytype) :: dnorm,dl,s
  real(mytype) :: du(3)
  integer :: l
  valid=0; rejected=1; state%valid_flag=0; state%rejected_flag=1
  if(state%allocated_flag/=1 .or. state%nlag<2) return
  if(length<=0 .or. .not.ieee_is_finite(length)) return
  dnorm=sqrt(sum(direction**2)); if(dnorm<=0 .or. .not.ieee_is_finite(dnorm)) return
  du=direction/dnorm; dl=length/real(state%nlag-1,mytype)
  do l=1,state%nlag
    s=-0.5_mytype*length + real(l-1,mytype)*dl
    state%x(:,l)=[x0,y0,z0] + s*du
    if(.not.all(ieee_is_finite(state%x(:,l)))) return
  end do
  state%ds=dl; state%ds(1)=0.5_mytype*dl; state%ds(state%nlag)=0.5_mytype*dl
  if(any(state%ds<=0) .or. .not.all(ieee_is_finite(state%ds))) return
  state%valid_flag=1; state%rejected_flag=0; valid=1; rejected=0
end subroutine
subroutine validate_stage8_lagrangian_state(state,grid,layout,valid,rejected,safe_count,blocked_count,unsafe_count)
  type(stage8_lagrangian_state_t), intent(inout) :: state
  type(stage7_channel_grid_t), intent(in) :: grid
  type(stage7_velocity_layout_t), intent(in) :: layout
  integer, intent(out) :: valid,rejected,safe_count,blocked_count,unsafe_count
  type(stage7_boundary_safety_result_t) :: res
  integer :: l
  valid=0; rejected=1; safe_count=0; blocked_count=0; unsafe_count=0
  if(state%allocated_flag/=1 .or. state%nlag<2) return
  if(.not.all(ieee_is_finite(state%x)) .or. .not.all(ieee_is_finite(state%ds)) .or. any(state%ds<=0)) return
  do l=1,state%nlag
    call classify_stage7_boundary_point(grid,layout,state%x(1,l),state%x(2,l),state%x(3,l),1,res)
    state%point_valid_flag(l)=res%safe_flag
    state%point_blocked_flag(l)=res%blocked_flag
    state%point_unsafe_flag(l)=merge(1,0,res%safe_flag==0)
    state%point_status_code(l)=res%status_code
  end do
  safe_count=sum(state%point_valid_flag); blocked_count=sum(state%point_blocked_flag); unsafe_count=sum(state%point_unsafe_flag)
  if(safe_count==state%nlag .and. blocked_count==0 .and. unsafe_count==0) then
    valid=1; rejected=0; state%valid_flag=1; state%rejected_flag=0
  else
    valid=0; rejected=1; state%valid_flag=0; state%rejected_flag=1
  end if
end subroutine
subroutine compute_stage8_lagrangian_total_length(state,total_length)
  type(stage8_lagrangian_state_t), intent(in) :: state
  real(mytype), intent(out) :: total_length
  total_length=0; if(state%allocated_flag==1 .and. allocated(state%ds)) total_length=sum(state%ds)
end subroutine
subroutine compute_stage8_lagrangian_segment_error(state,expected_dl,err_max)
  type(stage8_lagrangian_state_t), intent(in) :: state
  real(mytype), intent(in) :: expected_dl
  real(mytype), intent(out) :: err_max
  integer :: l
  real(mytype) :: seg
  err_max=huge(1.0_mytype)
  if(state%allocated_flag/=1 .or. state%nlag<2) return
  err_max=0
  do l=1,state%nlag-1
    seg=sqrt(sum((state%x(:,l+1)-state%x(:,l))**2))
    err_max=max(err_max,abs(seg-expected_dl))
  end do
end subroutine
end module
