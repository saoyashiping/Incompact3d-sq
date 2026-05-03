module fibre_stage4_time_order
  use fibre_parameters, only : mytype
  use fibre_types, only : fibre_t
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre, destroy_lagrangian_points
  use fibre_ibm_boundary_safety, only : ibm_point_safety_t, init_point_safety, destroy_point_safety, check_ibm_point_boundary_safety, assert_no_unsafe_ibm_points
  use fibre_external_force, only : clear_fibre_external_force
  use fibre_structure_integrator, only : advance_fibre_structure_freefree
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t
  use fibre_stage4_interpolation_adapter, only : convert_stage4_adapter_to_ibm_grid
  use fibre_stage4_feedback_adapter, only : compute_stage4_feedback_if_supported, apply_stage4_feedback_to_f_ext
  implicit none
  private
  public :: perform_stage4_oneway_ordered_step
contains
subroutine perform_stage4_oneway_ordered_step(adapter, ux, uy, uz, fibre, beta_drag, dt, f_ext_match_error, advance_called, unsafe_count, blocked_flag, rhs_modified_flag, nan_detected)
 type(stage4_grid_adapter_t),intent(in)::adapter
 real(mytype),intent(in)::ux(:,:,:),uy(:,:,:),uz(:,:,:),beta_drag,dt
 type(fibre_t),intent(inout)::fibre
 real(mytype),intent(out)::f_ext_match_error
 integer,intent(out)::advance_called,unsafe_count,blocked_flag,rhs_modified_flag,nan_detected
 real(mytype),allocatable::u_lag(:,:),fs(:,:),ff(:,:)
 type(ibm_grid_t)::grid
 type(ibm_lagrangian_points_t)::lag
 type(ibm_point_safety_t)::safety
 real(mytype)::tr,rr
 integer::st
 f_ext_match_error=0._mytype; advance_called=0; unsafe_count=0; blocked_flag=0; rhs_modified_flag=0; nan_detected=0
 call init_lagrangian_points_from_fibre(lag,fibre)
 call convert_stage4_adapter_to_ibm_grid(adapter,grid,st)
 call init_point_safety(safety,lag%nl); call check_ibm_point_boundary_safety(grid,lag,safety); call assert_no_unsafe_ibm_points(safety,unsafe_count)
 call destroy_point_safety(safety); call destroy_lagrangian_points(lag)
 if (unsafe_count>0) then
   call clear_fibre_external_force(fibre)
   blocked_flag=1
   return
 end if
 allocate(u_lag(3,fibre%nl),fs(3,fibre%nl),ff(3,fibre%nl))
 call compute_stage4_feedback_if_supported(adapter,ux,uy,uz,fibre,beta_drag,u_lag,fs,ff,st)
 if (st/=1) then; blocked_flag=1; return; end if
 call apply_stage4_feedback_to_f_ext(fibre,fs,'set',st)
 f_ext_match_error=maxval(abs(fibre%f_ext-fs))
 call advance_fibre_structure_freefree(fibre,dt,tr,rr,st)
 advance_called=1
 if (any(fibre%x /= fibre%x)) nan_detected=1
end subroutine
end module
