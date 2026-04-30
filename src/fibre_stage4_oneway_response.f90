module fibre_stage4_oneway_response
  use fibre_parameters, only : mytype
  use fibre_types, only : fibre_t
  use fibre_diagnostics, only : compute_total_length_relative_error
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t, ibm_point_safety_t
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre, destroy_lagrangian_points
  use fibre_ibm_boundary_safety, only : init_point_safety, destroy_point_safety, check_ibm_point_boundary_safety, assert_no_unsafe_ibm_points
  use fibre_structure_integrator, only : advance_fibre_structure_freefree
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t
  use fibre_stage4_interpolation_adapter, only : convert_stage4_adapter_to_ibm_grid
  use fibre_stage4_feedback_adapter, only : compute_stage4_feedback_if_supported, apply_stage4_feedback_to_f_ext
  implicit none
  private
  public :: advance_fibre_oneway_stage4
contains
subroutine advance_fibre_oneway_stage4(adapter, ux, uy, uz, fibre, beta_drag, dt, nsteps, max_length_error, max_f_ext_norm, final_center_velocity, final_center_displacement, solver_failure_count, nan_detected, unsafe_count_max)
 type(stage4_grid_adapter_t),intent(in)::adapter
 real(mytype),intent(in)::ux(:,:,:),uy(:,:,:),uz(:,:,:),beta_drag,dt
 integer,intent(in)::nsteps
 type(fibre_t),intent(inout)::fibre
 real(mytype),intent(out)::max_length_error,max_f_ext_norm,final_center_velocity(3),final_center_displacement(3)
 integer,intent(out)::solver_failure_count,nan_detected,unsafe_count_max
 real(mytype),allocatable::u_lag(:,:),fs(:,:),ff(:,:)
 type(ibm_grid_t)::grid
 type(ibm_lagrangian_points_t)::lag
 type(ibm_point_safety_t)::safety
  real(mytype)::tr,rr,len_err
 integer::st,unsafe,ierr,cidx
 real(mytype)::center0(3)
 cidx=(fibre%nl+1)/2; center0=fibre%x(:,cidx)
 allocate(u_lag(3,fibre%nl),fs(3,fibre%nl),ff(3,fibre%nl))
 max_length_error=0._mytype;max_f_ext_norm=0._mytype;solver_failure_count=0;nan_detected=0;unsafe_count_max=0
 call convert_stage4_adapter_to_ibm_grid(adapter,grid,st)
 do st=1,nsteps
   call init_lagrangian_points_from_fibre(lag,fibre); call init_point_safety(safety,lag%nl); call check_ibm_point_boundary_safety(grid,lag,safety); call assert_no_unsafe_ibm_points(safety,unsafe)
   unsafe_count_max=max(unsafe_count_max,unsafe)
   call destroy_point_safety(safety); call destroy_lagrangian_points(lag)
   if (unsafe>0) exit
   call compute_stage4_feedback_if_supported(adapter,ux,uy,uz,fibre,beta_drag,u_lag,fs,ff,ierr)
   if (ierr/=1) then; solver_failure_count=solver_failure_count+1; exit; end if
   call apply_stage4_feedback_to_f_ext(fibre,fs,'set',ierr)
   max_f_ext_norm=max(max_f_ext_norm,sqrt(sum(fibre%f_ext**2)))
   call advance_fibre_structure_freefree(fibre,dt,tr,rr,ierr)
   if (ierr/=0) solver_failure_count=solver_failure_count+1
   call compute_total_length_relative_error(fibre, len_err)
   max_length_error=max(max_length_error,abs(len_err))
   if (any(fibre%x /= fibre%x)) nan_detected=1
 end do
 final_center_velocity=fibre%v(:,cidx)
 final_center_displacement=fibre%x(:,cidx)-center0
end subroutine
end module
