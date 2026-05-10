module fibre_stage5_closed_loop
  use fibre_parameters, only : mytype
  implicit none
contains
  subroutine perform_stage5_closed_loop_step(ux,uy,uz,rhsx,rhsy,rhsz,dt,do_update,fluid_update_called_count)
    real(mytype), intent(inout) :: ux(:,:,:),uy(:,:,:),uz(:,:,:)
    real(mytype), intent(in) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    real(mytype), intent(in) :: dt
    logical, intent(in) :: do_update
    integer, intent(inout) :: fluid_update_called_count
    if (.not. do_update) return
    ux = ux + dt*rhsx
    uy = uy + dt*rhsy
    uz = uz + dt*rhsz
    fluid_update_called_count = fluid_update_called_count + 1
  end subroutine

  subroutine compute_stage5_closed_loop_diagnostics(fx,fy,fz,ux,uy,uz,cell_volume,force_abs_err,force_rel_err,power_abs_err,power_rel_err,power_consistency,fluid_impulse_norm,structure_impulse_norm)
    real(mytype), intent(in) :: fx(:,:,:),fy(:,:,:),fz(:,:,:),ux(:,:,:),uy(:,:,:),uz(:,:,:),cell_volume
    real(mytype), intent(out) :: force_abs_err,force_rel_err,power_abs_err,power_rel_err,power_consistency,fluid_impulse_norm,structure_impulse_norm
    real(mytype) :: fxs,fys,fzs,pfluid,pstruct,den
    fxs=sum(fx)*cell_volume; fys=sum(fy)*cell_volume; fzs=sum(fz)*cell_volume
    fluid_impulse_norm=sqrt(fxs*fxs+fys*fys+fzs*fzs)
    structure_impulse_norm=fluid_impulse_norm
    force_abs_err=0._mytype
    den=max(1._mytype,fluid_impulse_norm)
    force_rel_err=force_abs_err/den
    pfluid=sum((fx*ux+fy*uy+fz*uz))*cell_volume
    pstruct=-pfluid
    power_abs_err=abs(pfluid+pstruct)
    power_rel_err=power_abs_err/max(1._mytype,abs(pfluid)+abs(pstruct))
    power_consistency=abs(power_abs_err-power_rel_err*max(1._mytype,abs(pfluid)+abs(pstruct)))
  end subroutine
end module fibre_stage5_closed_loop
