module fibre_stage4_analytic_channel
  use fibre_parameters, only : mytype
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t
  implicit none
  private
  public :: fill_stage4_poiseuille_velocity, fill_stage4_constant_velocity
contains
  subroutine fill_stage4_poiseuille_velocity(adapter, uc, ux, uy, uz)
    type(stage4_grid_adapter_t), intent(in) :: adapter
    real(mytype), intent(in) :: uc
    real(mytype), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    integer :: i,j,k
    real(mytype) :: yc,h
    yc = 0.5_mytype*(adapter%ymin+adapter%ymax)
    h = 0.5_mytype*(adapter%ymax-adapter%ymin)
    if (h <= 0._mytype) error stop 'fill_stage4_poiseuille_velocity: h must be > 0'
    do k=1,adapter%nz; do j=1,adapter%ny; do i=1,adapter%nx
      ux(i,j,k)=uc*(1._mytype-((adapter%y(j)-yc)/h)**2)
      uy(i,j,k)=0._mytype; uz(i,j,k)=0._mytype
    end do; end do; end do
  end subroutine
  subroutine fill_stage4_constant_velocity(adapter, u0, ux, uy, uz)
    type(stage4_grid_adapter_t), intent(in) :: adapter
    real(mytype), intent(in) :: u0(3)
    real(mytype), intent(out) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    ux=u0(1); uy=u0(2); uz=u0(3)
  end subroutine
end module
