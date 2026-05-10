module fibre_stage5_momentum_exchange
  use fibre_parameters, only : mytype
  implicit none
contains
  subroutine compute_stage5_fluid_momentum_change(rho, dvol, ux0,uy0,uz0, ux1,uy1,uz1, dp)
    real(mytype), intent(in) :: rho, dvol
    real(mytype), intent(in) :: ux0(:,:,:),uy0(:,:,:),uz0(:,:,:),ux1(:,:,:),uy1(:,:,:),uz1(:,:,:)
    real(mytype), intent(out) :: dp(3)
    dp(1)=rho*sum((ux1-ux0)*dvol); dp(2)=rho*sum((uy1-uy0)*dvol); dp(3)=rho*sum((uz1-uz0)*dvol)
  end subroutine
  subroutine compute_stage5_eulerian_force_impulse(dt, dvol, fx,fy,fz, impulse)
    real(mytype), intent(in) :: dt,dvol
    real(mytype), intent(in) :: fx(:,:,:),fy(:,:,:),fz(:,:,:)
    real(mytype), intent(out) :: impulse(3)
    impulse(1)=dt*sum(fx*dvol); impulse(2)=dt*sum(fy*dvol); impulse(3)=dt*sum(fz*dvol)
  end subroutine
  subroutine compute_stage5_lagrangian_force_impulse(dt, lag_force, lag_weight, impulse)
    real(mytype), intent(in) :: dt, lag_force(:,:), lag_weight(:)
    real(mytype), intent(out) :: impulse(3)
    integer :: l
    impulse=0._mytype
    do l=1,size(lag_weight)
      impulse(:)=impulse(:)+dt*lag_force(:,l)*lag_weight(l)
    end do
  end subroutine
end module
