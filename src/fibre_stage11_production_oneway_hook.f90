module fibre_stage11_production_oneway_hook
  use decomp_2d_constants, only : mytype
  use fibre_stage11_config, only : stage11_requested, stage11_readonly_mode
  implicit none
  private
  integer :: requested_flag=0, readonly_mode_status=0, hook_initialized_status=0, hook_sample_called_status=0
  integer :: sample_performed_status=0, sample_count_status=0, sampled_velocity_finite_status=0, sampled_velocity_signature_status=0
  integer :: field_modified_status=0, rhs_modified_status=0, no_rhs_injection_status=1, no_ibm_spreading_status=1, no_feedback_force_status=1, no_twoway_force_status=1, no_structure_advance_status=1
  integer :: production_oneway_hook_status=0, sample_count=0
  real(mytype) :: sample_sum_u=0, sample_sum_v=0, sample_sum_w=0, sample_l2_u=0, sample_l2_v=0, sample_l2_w=0
  public :: stage11_production_oneway_init, stage11_production_oneway_sample, stage11_production_oneway_finalize
  public :: stage11_production_oneway_get_status_values, stage11_production_oneway_write_diagnostics
contains
  subroutine stage11_production_oneway_init()
    requested_flag=merge(1,0,stage11_requested()); readonly_mode_status=merge(1,0,stage11_readonly_mode())
    hook_initialized_status=1; hook_sample_called_status=0; sample_performed_status=0; sample_count=0; sample_count_status=0
    sample_sum_u=0; sample_sum_v=0; sample_sum_w=0; sample_l2_u=0; sample_l2_v=0; sample_l2_w=0
    sampled_velocity_finite_status=1; sampled_velocity_signature_status=0; field_modified_status=0; rhs_modified_status=0
    call update_status()
  end subroutine
  subroutine stage11_production_oneway_sample(ux,uy,uz)
    real(mytype), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    integer :: i0,j0,k0
    real(mytype) :: su,sv,sw
    hook_sample_called_status=1
    i0=(lbound(ux,1)+ubound(ux,1))/2; j0=(lbound(ux,2)+ubound(ux,2))/2; k0=(lbound(ux,3)+ubound(ux,3))/2
    su=ux(i0,j0,k0); sv=uy(i0,j0,k0); sw=uz(i0,j0,k0)
    sample_count=sample_count+1
    sample_sum_u=sample_sum_u+su; sample_sum_v=sample_sum_v+sv; sample_sum_w=sample_sum_w+sw
    sample_l2_u=sample_l2_u+su*su; sample_l2_v=sample_l2_v+sv*sv; sample_l2_w=sample_l2_w+sw*sw
    sample_performed_status=1; sample_count_status=merge(1,0,sample_count>0)
    sampled_velocity_finite_status=merge(1,0,finitev(su).and.finitev(sv).and.finitev(sw))
    sampled_velocity_signature_status=merge(1,0,sample_count>0)
    field_modified_status=0; rhs_modified_status=0
    call update_status()
  end subroutine
  subroutine stage11_production_oneway_finalize()
    call update_status()
  end subroutine
  logical function finitev(x)
    real(mytype), intent(in) :: x
    finitev=(x==x).and.(abs(x)<huge(x))
  end function
  subroutine stage11_production_oneway_get_status_values(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)
    integer,intent(out)::a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p
    a=requested_flag; b=readonly_mode_status; c=hook_initialized_status; d=hook_sample_called_status
    e=sample_performed_status; f=sample_count_status; g=sampled_velocity_finite_status; h=sampled_velocity_signature_status
    i=field_modified_status; j=rhs_modified_status; k=no_rhs_injection_status; l=no_ibm_spreading_status
    m=no_feedback_force_status; n=no_twoway_force_status; o=no_structure_advance_status; p=production_oneway_hook_status
  end subroutine
  subroutine stage11_production_oneway_write_diagnostics(unit)
    integer,intent(in)::unit
    write(unit,'(A,1X,I0)') 'stage11_5_requested_flag', requested_flag
    write(unit,'(A,1X,I0)') 'stage11_5_readonly_mode_status', readonly_mode_status
    write(unit,'(A,1X,I0)') 'stage11_5_hook_initialized_status', hook_initialized_status
    write(unit,'(A,1X,I0)') 'stage11_5_hook_sample_called_status', hook_sample_called_status
    write(unit,'(A,1X,I0)') 'stage11_5_sample_performed_status', sample_performed_status
    write(unit,'(A,1X,I0)') 'stage11_5_sample_count_status', sample_count_status
    write(unit,'(A,1X,I0)') 'stage11_5_sampled_velocity_finite_status', sampled_velocity_finite_status
    write(unit,'(A,1X,I0)') 'stage11_5_sampled_velocity_signature_status', sampled_velocity_signature_status
    write(unit,'(A,1X,I0)') 'stage11_5_field_modified_status', field_modified_status
    write(unit,'(A,1X,I0)') 'stage11_5_rhs_modified_status', rhs_modified_status
    write(unit,'(A,1X,I0)') 'stage11_5_no_rhs_injection_status', no_rhs_injection_status
    write(unit,'(A,1X,I0)') 'stage11_5_no_ibm_spreading_status', no_ibm_spreading_status
    write(unit,'(A,1X,I0)') 'stage11_5_no_feedback_force_status', no_feedback_force_status
    write(unit,'(A,1X,I0)') 'stage11_5_no_twoway_force_status', no_twoway_force_status
    write(unit,'(A,1X,I0)') 'stage11_5_no_structure_advance_status', no_structure_advance_status
    write(unit,'(A,1X,I0)') 'stage11_5_production_oneway_hook_status', production_oneway_hook_status
    write(unit,'(A,1X,I0)') 'stage11_5_sample_count', sample_count
    write(unit,'(A,1X,ES24.16E3)') 'stage11_5_sample_sum_u', sample_sum_u
    write(unit,'(A,1X,ES24.16E3)') 'stage11_5_sample_sum_v', sample_sum_v
    write(unit,'(A,1X,ES24.16E3)') 'stage11_5_sample_sum_w', sample_sum_w
    write(unit,'(A,1X,ES24.16E3)') 'stage11_5_sample_l2_u', sample_l2_u
    write(unit,'(A,1X,ES24.16E3)') 'stage11_5_sample_l2_v', sample_l2_v
    write(unit,'(A,1X,ES24.16E3)') 'stage11_5_sample_l2_w', sample_l2_w
  end subroutine
  subroutine update_status()
    production_oneway_hook_status=merge(1,0,requested_flag==1.and.readonly_mode_status==1.and.hook_initialized_status==1.and.hook_sample_called_status==1.and.sample_performed_status==1.and.sample_count_status==1.and.sampled_velocity_finite_status==1.and.sampled_velocity_signature_status==1.and.field_modified_status==0.and.rhs_modified_status==0.and.no_rhs_injection_status==1.and.no_ibm_spreading_status==1.and.no_feedback_force_status==1.and.no_twoway_force_status==1.and.no_structure_advance_status==1)
  end subroutine
end module
