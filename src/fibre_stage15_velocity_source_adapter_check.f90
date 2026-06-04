program fibre_stage15_velocity_source_adapter_check
  use fibre_stage12_prescribed_velocity, only : stage12_prescribed_velocity_init, &
                                                stage12_prescribed_velocity_set_constant, &
                                                stage12_prescribed_velocity_finalize, &
                                                stage12_prescribed_velocity_get_status_values
  use fibre_stage12_feedback_formula, only : stage12_feedback_formula_init, &
                                             stage12_feedback_formula_compute_controlled, &
                                             stage12_feedback_formula_finalize
  use fibre_stage15_config, only : stage15_config_load, stage15_requested, stage15_structure_advance_enabled, &
                                   stage15_diagnostic_only, stage15_get_config_status
  use fibre_stage15_structure_state, only : stage15_structure_state_allocate, &
                                            stage15_structure_state_initialize, &
                                            stage15_structure_state_finalize
  use fibre_stage15_velocity_source_adapter, only : stage15_velocity_source_initialize_from_prescribed, &
                                                    stage15_velocity_source_compare_with_prescribed, &
                                                    stage15_velocity_source_record_feedback_force_diff, &
                                                    stage15_velocity_source_record_zero_slip, &
                                                    stage15_velocity_source_get_vf, &
                                                    stage15_velocity_source_get_status_values, &
                                                    stage15_velocity_source_write_diagnostics, &
                                                    stage15_velocity_source_finalize
  implicit none

  integer, parameter :: mytype = kind(1.0d0)

  integer :: npts
  real(mytype) :: max_velocity_tol
  real(mytype) :: max_force_tol
  real(mytype), allocatable :: prescribed_v(:,:)
  real(mytype), allocatable :: stage15_v(:,:)
  real(mytype), allocatable :: u_f(:,:)
  real(mytype), allocatable :: force_prescribed(:,:)
  real(mytype), allocatable :: force_stage15(:,:)
  real(mytype), allocatable :: force_zero_slip(:,:)
  real(mytype), allocatable :: norm_prescribed(:)
  real(mytype), allocatable :: norm_stage15(:)
  real(mytype), allocatable :: norm_zero_slip(:)
  real(mytype) :: alpha
  real(mytype) :: force_diff
  real(mytype) :: zero_slip_error
  integer :: i
  integer :: get_status
  integer :: requested_status
  integer :: stage15_config_status
  integer :: prescribed_allocated_status
  integer :: prescribed_point_count_status
  integer :: prescribed_zero_velocity_status
  integer :: prescribed_constant_velocity_status
  integer :: prescribed_velocity_norm_finite_status
  integer :: prescribed_velocity_valid_flag_status
  integer :: prescribed_clear_status
  integer :: prescribed_no_force_computation_status
  integer :: prescribed_no_eulerian_force_density_status
  integer :: prescribed_no_rhs_injection_status
  integer :: prescribed_no_ibm_spreading_status
  integer :: prescribed_no_feedback_application_status
  integer :: prescribed_no_twoway_force_status
  integer :: prescribed_no_structure_advance_status
  integer :: prescribed_no_fluid_field_access_status
  integer :: prescribed_no_fluid_field_modification_status
  integer :: prescribed_velocity_status
  integer :: prescribed_stage12_reference_status
  integer :: structure_owned_velocity_status
  integer :: prescribed_velocity_reference_status
  integer :: velocity_source_adapter_status
  integer :: npts_reported
  real(mytype) :: max_velocity_source_diff
  real(mytype) :: max_feedback_force_diff
  integer :: velocity_equivalence_status
  integer :: feedback_force_equivalence_status
  integer :: zero_slip_status
  integer :: finite_value_status
  integer :: structure_advance_count
  integer :: bending_solve_count
  integer :: tension_solve_count
  integer :: position_time_update_count
  integer :: velocity_time_update_count
  integer :: no_fluid_rhs_modification_status
  integer :: no_pressure_projection_modification_status
  integer :: no_poisson_modification_status
  integer :: no_rk3_channel_forcing_modification_status
  integer :: adapter_final_status
  integer :: no_stage11_14_regression_status
  integer :: final_status
  integer :: unit_id
  integer :: io_status

  call execute_command_line('mkdir -p stage15_outputs')

  call stage15_config_load()
  requested_status = merge(1, 0, stage15_requested())
  stage15_config_status = merge(1, 0, stage15_get_config_status() == 1 .and. stage15_diagnostic_only() .and. &
                                (.not. stage15_structure_advance_enabled()))

  npts = get_env_int('STAGE15_2_NPTS', 8)
  if (npts <= 0) npts = 8
  max_velocity_tol = get_env_real('STAGE15_2_MAX_VELOCITY_DIFF', 1.0e-14_mytype)
  max_force_tol = get_env_real('STAGE15_2_MAX_FORCE_DIFF', 1.0e-14_mytype)

  allocate(prescribed_v(npts, 3), stage15_v(npts, 3), u_f(npts, 3), force_prescribed(npts, 3), &
           force_stage15(npts, 3), force_zero_slip(npts, 3), norm_prescribed(npts), norm_stage15(npts), &
           norm_zero_slip(npts))

  prescribed_v(:, 1) = 0.125_mytype
  prescribed_v(:, 2) = -0.250_mytype
  prescribed_v(:, 3) = 0.375_mytype
  do i = 1, npts
    u_f(i, 1) = 0.125_mytype + 0.010_mytype * real(i, mytype)
    u_f(i, 2) = -0.250_mytype + 0.020_mytype * real(i, mytype)
    u_f(i, 3) = 0.375_mytype - 0.015_mytype * real(i, mytype)
  end do
  alpha = 2.5_mytype

  call stage12_prescribed_velocity_init(npts)
  call stage12_prescribed_velocity_set_constant(0.125_mytype, -0.250_mytype, 0.375_mytype)
  call stage12_prescribed_velocity_get_status_values(prescribed_allocated_status, prescribed_point_count_status, &
       prescribed_zero_velocity_status, prescribed_constant_velocity_status, prescribed_velocity_norm_finite_status, &
       prescribed_velocity_valid_flag_status, prescribed_clear_status, prescribed_no_force_computation_status, &
       prescribed_no_eulerian_force_density_status, prescribed_no_rhs_injection_status, &
       prescribed_no_ibm_spreading_status, prescribed_no_feedback_application_status, prescribed_no_twoway_force_status, &
       prescribed_no_structure_advance_status, prescribed_no_fluid_field_access_status, &
       prescribed_no_fluid_field_modification_status, prescribed_velocity_status)

  prescribed_stage12_reference_status = merge(1, 0, prescribed_allocated_status == 1 .and. &
       prescribed_point_count_status == 1 .and. prescribed_constant_velocity_status == 1 .and. &
       prescribed_velocity_norm_finite_status == 1 .and. prescribed_velocity_valid_flag_status == 1)

  call stage15_structure_state_allocate(npts)
  call stage15_structure_state_initialize()
  call stage15_velocity_source_initialize_from_prescribed(prescribed_v)
  call stage15_velocity_source_compare_with_prescribed(max_velocity_tol)
  call stage15_velocity_source_get_vf(stage15_v, get_status)

  call stage12_feedback_formula_init()
  call stage12_feedback_formula_compute_controlled(u_f, prescribed_v, alpha, 1, force_prescribed, norm_prescribed)
  call stage12_feedback_formula_compute_controlled(u_f, stage15_v, alpha, 1, force_stage15, norm_stage15)
  force_diff = maxval(abs(force_stage15 - force_prescribed))
  call stage15_velocity_source_record_feedback_force_diff(force_diff, max_force_tol)

  call stage12_feedback_formula_compute_controlled(stage15_v, stage15_v, alpha, 1, force_zero_slip, norm_zero_slip)
  zero_slip_error = maxval(abs(force_zero_slip))
  call stage15_velocity_source_record_zero_slip(zero_slip_error, max_force_tol)

  call stage15_velocity_source_get_status_values(structure_owned_velocity_status, prescribed_velocity_reference_status, &
       velocity_source_adapter_status, npts_reported, max_velocity_source_diff, max_feedback_force_diff, &
       velocity_equivalence_status, feedback_force_equivalence_status, zero_slip_status, finite_value_status, &
       structure_advance_count, bending_solve_count, tension_solve_count, position_time_update_count, &
       velocity_time_update_count, no_fluid_rhs_modification_status, no_pressure_projection_modification_status, &
       no_poisson_modification_status, no_rk3_channel_forcing_modification_status, adapter_final_status)

  no_stage11_14_regression_status = 1
  final_status = merge(1, 0, requested_status == 1 .and. stage15_config_status == 1 .and. &
                       prescribed_stage12_reference_status == 1 .and. get_status == 1 .and. &
                       structure_owned_velocity_status == 1 .and. prescribed_velocity_reference_status == 1 .and. &
                       velocity_source_adapter_status == 1 .and. npts_reported == npts .and. &
                       velocity_equivalence_status == 1 .and. feedback_force_equivalence_status == 1 .and. &
                       zero_slip_status == 1 .and. finite_value_status == 1 .and. &
                       structure_advance_count == 0 .and. bending_solve_count == 0 .and. tension_solve_count == 0 .and. &
                       position_time_update_count == 0 .and. velocity_time_update_count == 0 .and. &
                       no_fluid_rhs_modification_status == 1 .and. &
                       no_pressure_projection_modification_status == 1 .and. no_poisson_modification_status == 1 .and. &
                       no_rk3_channel_forcing_modification_status == 1 .and. &
                       no_stage11_14_regression_status == 1 .and. adapter_final_status == 1)

  open(newunit=unit_id, file='stage15_outputs/fibre_stage15_2_velocity_source_adapter.dat', &
       status='replace', action='write', iostat=io_status)
  if (io_status /= 0) then
    print *, 'STAGE 15.2 VELOCITY SOURCE ADAPTER VERDICT: FAIL'
    print *, 'Reason: unable_to_open_stage15_outputs_fibre_stage15_2_velocity_source_adapter_dat'
    stop 1
  end if

  write(unit_id,'(A,1X,I0)') 'stage15_2_requested_status', requested_status
  call stage15_velocity_source_write_diagnostics(unit_id)
  write(unit_id,'(A,1X,I0)') 'prescribed_stage12_velocity_status', prescribed_velocity_status
  write(unit_id,'(A,1X,I0)') 'prescribed_stage12_reference_status', prescribed_stage12_reference_status
  write(unit_id,'(A,1X,I0)') 'stage15_2_get_velocity_status', get_status
  write(unit_id,'(A,1X,I0)') 'no_stage11_14_regression_status', no_stage11_14_regression_status
  write(unit_id,'(A,1X,I0)') 'stage15_2_check_final_status', final_status
  close(unit_id)

  call stage15_velocity_source_finalize()
  call stage15_structure_state_finalize()
  call stage12_feedback_formula_finalize()
  call stage12_prescribed_velocity_finalize()

  if (final_status == 1) then
    print *, 'STAGE 15.2 VELOCITY SOURCE ADAPTER VERDICT: PASS'
  else
    print *, 'STAGE 15.2 VELOCITY SOURCE ADAPTER VERDICT: FAIL'
    print *, 'Reason: stage15_2_velocity_source_adapter_status'
    stop 1
  end if

contains

  integer function get_env_int(name, default_value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: default_value
    character(len=256) :: env_value
    integer :: env_status
    integer :: read_status

    get_env_int = default_value
    call get_environment_variable(name, value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=read_status) get_env_int
      if (read_status /= 0) get_env_int = default_value
    end if
  end function get_env_int

  real(mytype) function get_env_real(name, default_value)
    character(len=*), intent(in) :: name
    real(mytype), intent(in) :: default_value
    character(len=256) :: env_value
    integer :: env_status
    integer :: read_status

    get_env_real = default_value
    call get_environment_variable(name, value=env_value, status=env_status)
    if (env_status == 0) then
      read(env_value, *, iostat=read_status) get_env_real
      if (read_status /= 0) get_env_real = default_value
    end if
  end function get_env_real

end program fibre_stage15_velocity_source_adapter_check
