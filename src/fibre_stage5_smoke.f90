module fibre_stage5_smoke
  use fibre_parameters, only : mytype
  implicit none
contains
  subroutine perform_stage5_smoke_step(ux,uy,uz,rhsx,rhsy,rhsz,dt,do_update,fluid_update_called_count)
    real(mytype), intent(inout) :: ux(:,:,:),uy(:,:,:),uz(:,:,:)
    real(mytype), intent(in) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),dt
    logical, intent(in) :: do_update
    integer, intent(inout) :: fluid_update_called_count
    if (.not. do_update) return
    ux=ux+dt*rhsx; uy=uy+dt*rhsy; uz=uz+dt*rhsz
    fluid_update_called_count=fluid_update_called_count+1
  end subroutine
end module fibre_stage5_smoke
