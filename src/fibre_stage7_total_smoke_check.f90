program fibre_stage7_total_smoke_check
  use fibre_parameters, only : mytype
  implicit none
  integer :: io, status6_smoke, status6_prior, final_status
  integer :: stage6_closed_exists, stage6_out_exists, stage6_dep
  integer :: out_exists(10), st(10), all_out_exists, all_st_ok
  integer :: default_prod_disabled, controlled_rhs_ok, force_density_ok, no_rho_div_ok, prod_reject_ok, safety_ok
  integer :: pp_untouched, proj_untouched, dns_not_called, fluid_not_called, fibre_not_called, noop_ok
  integer :: grid_num_ok, interp_num_ok, vel_num_ok, spread_num_ok, power_num_ok, adapter_num_ok, rhs_num_ok, numeric_ok
  integer :: marker_written, marker_status
  integer :: i
  character(len=*), parameter :: stage6_marker='stage6_outputs/STAGE6_CLOSED.md'
  character(len=*), parameter :: stage6_dat='stage6_outputs/fibre_stage6_total_smoke_check.dat'
  character(len=*), parameter :: stage7_files(10) = [ character(len=72) :: &
    'stage7_outputs/fibre_stage7_config_check.dat', &
    'stage7_outputs/fibre_stage7_grid_metadata_check.dat', &
    'stage7_outputs/fibre_stage7_scalar_interpolation_check.dat', &
    'stage7_outputs/fibre_stage7_scalar_interpolation_robustness_check.dat', &
    'stage7_outputs/fibre_stage7_velocity_interpolation_check.dat', &
    'stage7_outputs/fibre_stage7_force_spreading_check.dat', &
    'stage7_outputs/fibre_stage7_power_consistency_check.dat', &
    'stage7_outputs/fibre_stage7_boundary_safety_check.dat', &
    'stage7_outputs/fibre_stage7_channel_grid_adapter_check.dat', &
    'stage7_outputs/fibre_stage7_rhs_candidate_check.dat' ]

  call ensure_dir('stage7_outputs')
  call file_exists_int(stage6_marker, stage6_closed_exists)
  call file_exists_int(stage6_dat, stage6_out_exists)
  call get_int_value(stage6_dat,'stage6_total_smoke_check_status',status6_smoke)
  call get_int_value(stage6_dat,'stage6_total_all_prior_outputs_exist',status6_prior)
  stage6_dep = merge(1,0,stage6_closed_exists==1 .and. stage6_out_exists==1 .and. status6_smoke==1 .and. status6_prior==1)

  do i=1,10
    call file_exists_int(stage7_files(i), out_exists(i))
  end do
  all_out_exists = merge(1,0,all(out_exists==1))

  call get_int_value(stage7_files(1),'stage7_config_check_status',st(1))
  call get_int_value(stage7_files(2),'stage7_grid_metadata_check_status',st(2))
  call get_int_value(stage7_files(3),'stage7_scalar_interpolation_check_status',st(3))
  call get_int_value(stage7_files(4),'stage7_scalar_interpolation_robustness_check_status',st(4))
  call get_int_value(stage7_files(5),'stage7_velocity_interpolation_check_status',st(5))
  call get_int_value(stage7_files(6),'stage7_force_spreading_check_status',st(6))
  call get_int_value(stage7_files(7),'stage7_power_consistency_check_status',st(7))
  call get_int_value(stage7_files(8),'stage7_boundary_safety_check_status',st(8))
  call get_int_value(stage7_files(9),'stage7_channel_grid_adapter_check_status',st(9))
  call get_int_value(stage7_files(10),'stage7_rhs_candidate_check_status',st(10))
  all_st_ok = merge(1,0,all(st==1))

  default_prod_disabled = merge(1,0,get_flag(stage7_files(1),'stage7_default_production_dns_enabled')==0 .and. &
    get_flag(stage7_files(1),'stage7_default_production_two_way_enabled')==0)
  controlled_rhs_ok = merge(1,0,get_flag(stage7_files(1),'stage7_default_rhs_allowed_flag')==0 .and. &
    get_flag(stage7_files(1),'stage7_controlled_rhs_allowed_flag')==1 .and. get_flag(stage7_files(10),'stage7_rhs_controlled_hook_called_flag')==1 .and. &
    get_flag(stage7_files(10),'stage7_rhs_controlled_injected_flag')==1)
  force_density_ok = merge(1,0,get_flag(stage7_files(6),'stage7_spread_force_density_convention_flag')==1 .and. &
    get_flag(stage7_files(7),'stage7_power_force_density_convention_flag')==1)
  no_rho_div_ok = merge(1,0,get_flag(stage7_files(6),'stage7_spread_no_rho_division_flag')==1 .and. &
    get_flag(stage7_files(7),'stage7_power_no_rho_division_flag')==1 .and. get_flag(stage7_files(7),'stage7_power_stage6_rhs_hook_called_flag')==0 .and. &
    get_flag(stage7_files(10),'stage7_rhs_double_division_detected_flag')==0 .and. get_flag(stage7_files(10),'stage7_rhs_default_injected_flag')==0 .and. &
    get_flag(stage7_files(10),'stage7_rhs_default_modified_flag')==0 .and. get_flag(stage7_files(10),'stage7_rhs_production_injected_flag')==0)
  prod_reject_ok = merge(1,0,get_flag(stage7_files(1),'stage7_invalid_production_dns_rejected_flag')==1 .and. &
    get_flag(stage7_files(1),'stage7_invalid_production_twoway_rejected_flag')==1 .and. get_flag(stage7_files(10),'stage7_rhs_production_dns_rejected_flag')==1 .and. &
    get_flag(stage7_files(10),'stage7_rhs_production_twoway_rejected_flag')==1)
  safety_ok = merge(1,0,default_prod_disabled==1 .and. controlled_rhs_ok==1 .and. force_density_ok==1 .and. no_rho_div_ok==1 .and. prod_reject_ok==1)

  pp_untouched = merge(1,0,all([get_flag(stage7_files(1),'stage7_config_pressure_poisson_modified_flag'),get_flag(stage7_files(2),'stage7_grid_pressure_poisson_modified_flag'), &
    get_flag(stage7_files(3),'stage7_interp_pressure_poisson_modified_flag'),get_flag(stage7_files(5),'stage7_vel_pressure_poisson_modified_flag'), &
    get_flag(stage7_files(6),'stage7_spread_pressure_poisson_modified_flag'),get_flag(stage7_files(7),'stage7_power_pressure_poisson_modified_flag'), &
    get_flag(stage7_files(8),'stage7_boundary_pressure_poisson_modified_flag'),get_flag(stage7_files(9),'stage7_adapter_pressure_poisson_modified_flag'), &
    get_flag(stage7_files(10),'stage7_rhs_pressure_poisson_modified_flag')]==0))
  proj_untouched = merge(1,0,all([get_flag(stage7_files(1),'stage7_config_projection_modified_flag'),get_flag(stage7_files(2),'stage7_grid_projection_modified_flag'), &
    get_flag(stage7_files(3),'stage7_interp_projection_modified_flag'),get_flag(stage7_files(5),'stage7_vel_projection_modified_flag'), &
    get_flag(stage7_files(6),'stage7_spread_projection_modified_flag'),get_flag(stage7_files(7),'stage7_power_projection_modified_flag'), &
    get_flag(stage7_files(8),'stage7_boundary_projection_modified_flag'),get_flag(stage7_files(9),'stage7_adapter_projection_modified_flag'), &
    get_flag(stage7_files(10),'stage7_rhs_projection_modified_flag'),get_flag(stage7_files(1),'stage7_config_real_projection_called_flag'), &
    get_flag(stage7_files(10),'stage7_rhs_real_projection_called_flag')]==0))
  dns_not_called = merge(1,0,all([get_flag(stage7_files(1),'stage7_config_production_dns_called_flag'),get_flag(stage7_files(10),'stage7_rhs_production_dns_called_flag')]==0))
  fluid_not_called = merge(1,0,all([get_flag(stage7_files(1),'stage7_config_fluid_update_called_flag'),get_flag(stage7_files(10),'stage7_rhs_fluid_update_called_flag')]==0))
  fibre_not_called = merge(1,0,all([get_flag(stage7_files(1),'stage7_config_fibre_advance_called_flag'),get_flag(stage7_files(10),'stage7_rhs_fibre_advance_called_flag')]==0))
  noop_ok = merge(1,0,pp_untouched==1 .and. proj_untouched==1 .and. dns_not_called==1 .and. fluid_not_called==1 .and. fibre_not_called==1)

  grid_num_ok = merge(1,0,get_real(stage7_files(2),'stage7_grid_total_volume_error')<=1e-12_mytype .and. &
    get_real(stage7_files(2),'stage7_grid_dy_face_consistency_error_max')<=1e-12_mytype .and. get_real(stage7_files(2),'stage7_grid_volume_formula_error_max')<=1e-12_mytype)
  interp_num_ok = merge(1,0,get_real(stage7_files(3),'stage7_interp_weight_sum_error_max')<=1e-12_mytype .and. get_real(stage7_files(3),'stage7_interp_constant_error_max')<=1e-12_mytype .and. &
    get_real(stage7_files(3),'stage7_interp_linear_error_max')<=1e-11_mytype .and. get_real(stage7_files(3),'stage7_interp_poiseuille_error_max')<=1e-11_mytype .and. &
    get_real(stage7_files(3),'stage7_interp_periodic_wrap_error_max')<=1e-12_mytype .and. get_real(stage7_files(4),'stage7_interp_robust_weight_sum_error_max')<=1e-12_mytype .and. &
    get_real(stage7_files(4),'stage7_interp_robust_cubic_y_error_max')<=1e-11_mytype .and. get_real(stage7_files(4),'stage7_interp_robust_periodic_shift_error_max')<=1e-12_mytype)
  vel_num_ok = merge(1,0,get_real(stage7_files(5),'stage7_vel_constant_error_max')<=1e-12_mytype .and. get_real(stage7_files(5),'stage7_vel_linear_error_max')<=1e-11_mytype .and. &
    get_real(stage7_files(5),'stage7_vel_poiseuille_error_max')<=1e-11_mytype .and. get_real(stage7_files(5),'stage7_vel_periodic_x_shift_error_max')<=1e-12_mytype .and. &
    get_real(stage7_files(5),'stage7_vel_periodic_z_shift_error_max')<=1e-12_mytype)
  spread_num_ok = merge(1,0,get_real(stage7_files(6),'stage7_spread_single_force_relative_error')<=1e-12_mytype .and. &
    get_real(stage7_files(6),'stage7_spread_multipoint_force_relative_error')<=1e-12_mytype .and. get_real(stage7_files(6),'stage7_spread_volume_scaling_error_max')<=1e-12_mytype)
  power_num_ok = merge(1,0,get_real(stage7_files(7),'stage7_power_single_relative_error')<=1e-12_mytype .and. &
    get_real(stage7_files(7),'stage7_power_multipoint_relative_error')<=1e-12_mytype .and. get_real(stage7_files(7),'stage7_power_periodic_relative_error')<=1e-12_mytype)
  adapter_num_ok = merge(1,0,get_real(stage7_files(9),'stage7_adapter_metadata_match_error_max')<=1e-12_mytype .and. get_real(stage7_files(9),'stage7_adapter_volume_formula_error_max')<=1e-12_mytype)
  rhs_num_ok = merge(1,0,get_real(stage7_files(10),'stage7_rhs_controlled_increment_error_max')<=1e-12_mytype .and. get_real(stage7_files(10),'stage7_rhs_rho_scaling_error_max')<=1e-12_mytype .and. &
    get_real(stage7_files(10),'stage7_rhs_candidate_force_conservation_error')<=1e-12_mytype)
  numeric_ok = merge(1,0,grid_num_ok==1 .and. interp_num_ok==1 .and. vel_num_ok==1 .and. spread_num_ok==1 .and. power_num_ok==1 .and. adapter_num_ok==1 .and. rhs_num_ok==1)

  final_status = merge(1,0,stage6_dep==1 .and. all_out_exists==1 .and. all_st_ok==1 .and. safety_ok==1 .and. noop_ok==1 .and. numeric_ok==1)
  if (final_status==1) then
    call write_stage7_closed_marker('stage7_outputs/STAGE7_CLOSED.md', marker_status)
    call file_exists_int('stage7_outputs/STAGE7_CLOSED.md', marker_written)
  else
    marker_status = 0
    marker_written = 0
    call remove_file_if_exists('stage7_outputs/STAGE7_CLOSED.md')
  end if

  open(newunit=io,file='stage7_outputs/fibre_stage7_total_smoke_check.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage7_total_stage6_closed_marker_exists', stage6_closed_exists
  write(io,'(A,1X,I0)') 'stage7_total_stage6_total_smoke_output_exists', stage6_out_exists
  write(io,'(A,1X,I0)') 'stage7_total_stage6_total_smoke_status', status6_smoke
  write(io,'(A,1X,I0)') 'stage7_total_stage6_all_prior_outputs_status', status6_prior
  write(io,'(A,1X,I0)') 'stage7_total_stage6_dependency_status', stage6_dep
  do i=1,10
    write(io,'(A,I0,A,1X,I0)') 'stage7_total_stage7_', i-1, '_output_exists', out_exists(i)
  end do
  write(io,'(A,1X,I0)') 'stage7_total_all_stage7_outputs_exist', all_out_exists
  do i=1,10
    write(io,'(A,I0,A,1X,I0)') 'stage7_total_stage7_', i-1, '_status', st(i)
  end do
  write(io,'(A,1X,I0)') 'stage7_total_all_stage7_status', all_st_ok
  write(io,'(A,1X,I0)') 'stage7_total_default_production_disabled_status', default_prod_disabled
  write(io,'(A,1X,I0)') 'stage7_total_controlled_rhs_path_status', controlled_rhs_ok
  write(io,'(A,1X,I0)') 'stage7_total_force_density_convention_status', force_density_ok
  write(io,'(A,1X,I0)') 'stage7_total_no_double_rho_division_status', no_rho_div_ok
  write(io,'(A,1X,I0)') 'stage7_total_production_rejection_status', prod_reject_ok
  write(io,'(A,1X,I0)') 'stage7_total_safety_summary_status', safety_ok
  write(io,'(A,1X,I0)') 'stage7_total_pressure_poisson_untouched_status', pp_untouched
  write(io,'(A,1X,I0)') 'stage7_total_projection_untouched_status', proj_untouched
  write(io,'(A,1X,I0)') 'stage7_total_production_dns_not_called_status', dns_not_called
  write(io,'(A,1X,I0)') 'stage7_total_fluid_update_not_called_status', fluid_not_called
  write(io,'(A,1X,I0)') 'stage7_total_fibre_advance_not_called_status', fibre_not_called
  write(io,'(A,1X,I0)') 'stage7_total_noop_safety_summary_status', noop_ok
  write(io,'(A,1X,I0)') 'stage7_total_grid_numeric_status', grid_num_ok
  write(io,'(A,1X,I0)') 'stage7_total_interpolation_numeric_status', interp_num_ok
  write(io,'(A,1X,I0)') 'stage7_total_velocity_numeric_status', vel_num_ok
  write(io,'(A,1X,I0)') 'stage7_total_spreading_numeric_status', spread_num_ok
  write(io,'(A,1X,I0)') 'stage7_total_power_numeric_status', power_num_ok
  write(io,'(A,1X,I0)') 'stage7_total_adapter_numeric_status', adapter_num_ok
  write(io,'(A,1X,I0)') 'stage7_total_rhs_numeric_status', rhs_num_ok
  write(io,'(A,1X,I0)') 'stage7_total_numeric_summary_status', numeric_ok
  write(io,'(A,1X,I0)') 'stage7_total_closed_marker_written_flag', marker_written
  write(io,'(A,1X,I0)') 'stage7_total_closed_marker_status', marker_status
  write(io,'(A,1X,I0)') 'stage7_total_smoke_check_status', final_status
  close(io)
contains
  subroutine ensure_dir(path)
    character(len=*), intent(in) :: path
    call execute_command_line('mkdir -p ' // trim(path))
  end subroutine
  subroutine remove_file_if_exists(path)
    character(len=*), intent(in) :: path
    logical :: ex
    inquire(file=path, exist=ex)
    if (ex) call execute_command_line('rm -f ' // trim(path))
  end subroutine
  subroutine file_exists_int(path, flag)
    character(len=*), intent(in) :: path
    integer, intent(out) :: flag
    logical :: ex
    inquire(file=path, exist=ex); flag=merge(1,0,ex)
  end subroutine
  subroutine get_int_value(path,key,val)
    character(len=*), intent(in) :: path,key
    integer, intent(out) :: val
    real(mytype) :: tmp
    integer :: found
    call get_value(path,key,tmp,found); val=0; if(found==1) val=nint(tmp)
  end subroutine
  function get_flag(path,key) result(v)
    character(len=*), intent(in) :: path,key
    integer :: v
    call get_int_value(path,key,v)
  end function
  function get_real(path,key) result(v)
    character(len=*), intent(in) :: path,key
    real(mytype) :: v
    integer :: found
    call get_value(path,key,v,found)
    if(found/=1) v=huge(1._mytype)
  end function
  subroutine get_value(path,key,val,found)
    character(len=*), intent(in) :: path,key
    real(mytype), intent(out) :: val
    integer, intent(out) :: found
    integer :: u,ios
    character(len=256) :: k
    real(mytype) :: x
    logical :: ex
    found=0; val=0._mytype
    inquire(file=path,exist=ex); if(.not.ex) return
    open(newunit=u,file=path,status='old',action='read',iostat=ios); if(ios/=0)return
    do
      read(u,*,iostat=ios) k, x
      if(ios/=0) exit
      if(trim(k)==trim(key)) then
        val=x; found=1; exit
      end if
    end do
    close(u)
  end subroutine
  subroutine write_stage7_closed_marker(path,status)
    character(len=*), intent(in) :: path
    integer, intent(out) :: status
    integer :: u,ios
    open(newunit=u,file=path,status='replace',action='write',iostat=ios)
    if(ios/=0) then; status=0; return; end if
    write(u,'(A)') '# Stage 7 Closed'
    write(u,'(A)') ''
    write(u,'(A)') 'Stage 7.0 through Stage 7.9 have passed.'
    write(u,'(A)') 'Stage 7 real-layout nonuniform-y IBM foundation is closed.'
    write(u,'(A)') 'Production DNS and production two-way coupling remain disabled by default.'
    write(u,'(A)') 'RHS candidate is enabled only through controlled test configuration.'
    write(u,'(A)') 'Pressure/projection/DNS production paths were not modified by Stage 7 checks.'
    close(u); status=1
  end subroutine
end program fibre_stage7_total_smoke_check
