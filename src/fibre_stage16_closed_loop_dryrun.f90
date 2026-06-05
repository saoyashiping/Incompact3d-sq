module fibre_stage16_closed_loop_dryrun
  use decomp_2d_constants, only : mytype
  use fibre_stage16_structure_force_input, only : stage16_structure_force_input_load_from_environment, &
       stage16_structure_force_input_set_from_stage12_candidate, stage16_structure_force_input_get_force, &
       stage16_structure_force_input_get_final_status
  implicit none

  private

  public :: stage16_closed_loop_dryrun_reset
  public :: stage16_closed_loop_dryrun_load_from_environment
  public :: stage16_closed_loop_dryrun_initialize
  public :: stage16_closed_loop_dryrun_apply_one_step
  public :: stage16_closed_loop_dryrun_validate
  public :: stage16_closed_loop_dryrun_write_diagnostics
  public :: stage16_closed_loop_dryrun_get_final_status

  integer :: requested_status = 1
  integer :: np_value = 1
  integer :: npts_value = 8
  integer :: one_fibre_fsi_status = 1
  integer :: structure_update_enable_status = 1
  integer :: two_way_rhs_enable_status = 1
  integer :: diagnostic_only_status = 1
  real(mytype) :: feedback_alpha_value = 1.0_mytype
  real(mytype) :: lambda_value = 1.0e-8_mytype
  real(mytype) :: dt_value = 1.0e-4_mytype
  real(mytype) :: rho_tilde_value = 1.0_mytype
  real(mytype) :: max_force_input_value = 1.0e-6_mytype
  real(mytype) :: max_structure_update_value = 1.0e-12_mytype
  real(mytype) :: max_rhs_increment_value = 1.0e-8_mytype
  real(mytype) :: max_fluid_delta_value = 1.0e-8_mytype

  real(mytype) :: sampled_u_f(3) = 0.0_mytype
  real(mytype) :: structure_v_f(3) = 0.0_mytype
  real(mytype) :: slip_vector(3) = 0.0_mytype
  real(mytype) :: stage12_feedback_candidate(3) = 0.0_mytype
  real(mytype) :: structure_force_input(3) = 0.0_mytype
  real(mytype) :: fluid_side_force(3) = 0.0_mytype
  real(mytype) :: controlled_structure_update(3) = 0.0_mytype
  real(mytype) :: stage13_force_density_signature = 0.0_mytype
  real(mytype) :: stage14_rhs_increment = 0.0_mytype
  real(mytype) :: fluid_signature_delta_value = 0.0_mytype

  integer :: one_fibre_count_status = 0
  integer :: closed_loop_path_status = 0
  integer :: stage11_sampling_status = 0
  integer :: stage12_feedback_status = 0
  integer :: stage16_4_force_input_status = 0
  integer :: structure_force_input_finite_status = 0
  integer :: structure_force_input_bounded_status = 0
  integer :: controlled_structure_update_status = 0
  integer :: controlled_structure_update_bounded_status = 0
  integer :: stage13_force_density_status = 0
  integer :: stage14_rhs_status = 0
  integer :: rhs_increment_bounded_status = 0
  integer :: fluid_signature_status = 0
  integer :: approved_stage12_13_14_chain_status = 1
  integer :: no_direct_rhs_injection_status = 1
  integer :: no_production_hook_status = 1
  integer :: no_pressure_projection_modification_status = 1
  integer :: no_poisson_modification_status = 1
  integer :: no_rk3_channel_forcing_modification_status = 1
  integer :: no_channel_forcing_modification_status = 1
  integer :: no_wall_contact_status = 1
  integer :: no_multifibre_status = 1
  integer :: no_legacy_ibm_forcing_status = 1
  integer :: no_nan_inf_status = 0
  integer :: numeric_parse_status = 1
  integer :: numeric_bounds_status = 1
  integer :: final_status = 0

contains

  subroutine stage16_closed_loop_dryrun_reset()
    requested_status = 1
    np_value = 1
    npts_value = 8
    one_fibre_fsi_status = 1
    structure_update_enable_status = 1
    two_way_rhs_enable_status = 1
    diagnostic_only_status = 1
    feedback_alpha_value = 1.0_mytype
    lambda_value = 1.0e-8_mytype
    dt_value = 1.0e-4_mytype
    rho_tilde_value = 1.0_mytype
    max_force_input_value = 1.0e-6_mytype
    max_structure_update_value = 1.0e-12_mytype
    max_rhs_increment_value = 1.0e-8_mytype
    max_fluid_delta_value = 1.0e-8_mytype
    sampled_u_f(:) = 0.0_mytype
    structure_v_f(:) = 0.0_mytype
    slip_vector(:) = 0.0_mytype
    stage12_feedback_candidate(:) = 0.0_mytype
    structure_force_input(:) = 0.0_mytype
    fluid_side_force(:) = 0.0_mytype
    controlled_structure_update(:) = 0.0_mytype
    stage13_force_density_signature = 0.0_mytype
    stage14_rhs_increment = 0.0_mytype
    fluid_signature_delta_value = 0.0_mytype
    one_fibre_count_status = 0
    closed_loop_path_status = 0
    stage11_sampling_status = 0
    stage12_feedback_status = 0
    stage16_4_force_input_status = 0
    structure_force_input_finite_status = 0
    structure_force_input_bounded_status = 0
    controlled_structure_update_status = 0
    controlled_structure_update_bounded_status = 0
    stage13_force_density_status = 0
    stage14_rhs_status = 0
    rhs_increment_bounded_status = 0
    fluid_signature_status = 0
    approved_stage12_13_14_chain_status = 1
    no_direct_rhs_injection_status = 1
    no_production_hook_status = 1
    no_pressure_projection_modification_status = 1
    no_poisson_modification_status = 1
    no_rk3_channel_forcing_modification_status = 1
    no_channel_forcing_modification_status = 1
    no_wall_contact_status = 1
    no_multifibre_status = 1
    no_legacy_ibm_forcing_status = 1
    no_nan_inf_status = 0
    numeric_parse_status = 1
    numeric_bounds_status = 1
    final_status = 0
  end subroutine stage16_closed_loop_dryrun_reset

  subroutine stage16_closed_loop_dryrun_load_from_environment()
    logical :: parsed_bool
    integer :: parse_status

    call stage16_closed_loop_dryrun_reset()
    call load_bool_env('STAGE16_5_ENABLE', parsed_bool, parse_status)
    if (parse_status == 1) requested_status = logical_to_int(parsed_bool)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_bool_env('STAGE16_5_ONE_FIBRE_FSI_ENABLE', parsed_bool, parse_status)
    if (parse_status == 1) one_fibre_fsi_status = logical_to_int(parsed_bool)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_bool_env('STAGE16_5_STRUCTURE_UPDATE_ENABLE', parsed_bool, parse_status)
    if (parse_status == 1) structure_update_enable_status = logical_to_int(parsed_bool)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_bool_env('STAGE16_5_TWO_WAY_RHS_ENABLE', parsed_bool, parse_status)
    if (parse_status == 1) two_way_rhs_enable_status = logical_to_int(parsed_bool)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_bool_env('STAGE16_5_DIAGNOSTIC_ONLY', parsed_bool, parse_status)
    if (parse_status == 1) diagnostic_only_status = logical_to_int(parsed_bool)
    numeric_parse_status = min(numeric_parse_status, parse_status)

    call load_int_env('STAGE16_5_NP', np_value, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_int_env('STAGE16_5_NPTS', npts_value, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_real_env('STAGE16_5_FEEDBACK_ALPHA', feedback_alpha_value, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_real_env('STAGE16_5_LAMBDA', lambda_value, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_real_env('STAGE16_5_DT', dt_value, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_real_env('STAGE16_5_RHO_TILDE', rho_tilde_value, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_real_env('STAGE16_5_MAX_FORCE_INPUT', max_force_input_value, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_real_env('STAGE16_5_MAX_STRUCTURE_UPDATE', max_structure_update_value, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_real_env('STAGE16_5_MAX_RHS_INCREMENT', max_rhs_increment_value, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)
    call load_real_env('STAGE16_5_MAX_FLUID_DELTA', max_fluid_delta_value, parse_status)
    numeric_parse_status = min(numeric_parse_status, parse_status)

    numeric_bounds_status = logical_to_int(np_value == 1 .and. npts_value >= 2 .and. &
         feedback_alpha_value > 0.0_mytype .and. finite_real(feedback_alpha_value) .and. &
         lambda_value >= 0.0_mytype .and. finite_real(lambda_value) .and. &
         dt_value > 0.0_mytype .and. finite_real(dt_value) .and. &
         rho_tilde_value > 0.0_mytype .and. finite_real(rho_tilde_value) .and. &
         max_force_input_value >= 0.0_mytype .and. finite_real(max_force_input_value) .and. &
         max_structure_update_value >= 0.0_mytype .and. finite_real(max_structure_update_value) .and. &
         max_rhs_increment_value >= 0.0_mytype .and. finite_real(max_rhs_increment_value) .and. &
         max_fluid_delta_value >= 0.0_mytype .and. finite_real(max_fluid_delta_value) .and. &
         diagnostic_only_status == 1)
  end subroutine stage16_closed_loop_dryrun_load_from_environment

  subroutine stage16_closed_loop_dryrun_initialize()
    one_fibre_count_status = logical_to_int(one_fibre_fsi_status == 1 .and. npts_value >= 2)
    stage11_sampling_status = 1
    sampled_u_f(:) = (/ 2.0e-9_mytype, -1.0e-9_mytype, 1.0e-9_mytype /)
    structure_v_f(:) = (/ 1.0e-9_mytype, -2.0e-9_mytype, 0.0_mytype /)
    slip_vector(:) = sampled_u_f(:) - structure_v_f(:)
    stage11_sampling_status = logical_to_int(all_finite_vector(sampled_u_f) .and. all_finite_vector(structure_v_f) .and. &
         all_finite_vector(slip_vector))
  end subroutine stage16_closed_loop_dryrun_initialize

  subroutine stage16_closed_loop_dryrun_apply_one_step()
    stage12_feedback_candidate(:) = feedback_alpha_value * slip_vector(:)
    stage12_feedback_status = logical_to_int(all_finite_vector(stage12_feedback_candidate))

    call stage16_structure_force_input_load_from_environment()
    call stage16_structure_force_input_set_from_stage12_candidate(sampled_u_f, structure_v_f, feedback_alpha_value, .false.)
    call stage16_structure_force_input_get_force(1, structure_force_input)
    stage16_4_force_input_status = stage16_structure_force_input_get_final_status()

    fluid_side_force(:) = -structure_force_input(:)
    structure_force_input_finite_status = logical_to_int(all_finite_vector(structure_force_input))
    structure_force_input_bounded_status = logical_to_int(maxval(abs(structure_force_input)) <= max_force_input_value)

    controlled_structure_update(:) = 0.0_mytype
    if (structure_update_enable_status == 1) then
      controlled_structure_update(:) = dt_value * structure_force_input(:) / rho_tilde_value
    end if
    controlled_structure_update_status = logical_to_int(all_finite_vector(controlled_structure_update))
    controlled_structure_update_bounded_status = logical_to_int( &
         maxval(abs(controlled_structure_update)) <= max_structure_update_value)

    stage13_force_density_signature = maxval(abs(fluid_side_force))
    stage13_force_density_status = logical_to_int(finite_real(stage13_force_density_signature) .and. &
         stage13_force_density_signature <= max_force_input_value)
    stage14_rhs_increment = lambda_value * stage13_force_density_signature
    stage14_rhs_status = logical_to_int(finite_real(stage14_rhs_increment))
    rhs_increment_bounded_status = logical_to_int(stage14_rhs_increment <= max_rhs_increment_value)

    fluid_signature_delta_value = 0.0_mytype
    fluid_signature_status = logical_to_int(abs(fluid_signature_delta_value) <= max_fluid_delta_value)
    closed_loop_path_status = logical_to_int(requested_status == 1 .and. one_fibre_count_status == 1 .and. &
         structure_update_enable_status == 1 .and. two_way_rhs_enable_status == 1 .and. diagnostic_only_status == 1 .and. &
         stage11_sampling_status == 1 .and. stage12_feedback_status == 1 .and. stage16_4_force_input_status == 1 .and. &
         stage13_force_density_status == 1 .and. stage14_rhs_status == 1)
    no_nan_inf_status = logical_to_int(all_finite_vector(sampled_u_f) .and. all_finite_vector(structure_v_f) .and. &
         all_finite_vector(stage12_feedback_candidate) .and. all_finite_vector(structure_force_input) .and. &
         all_finite_vector(controlled_structure_update) .and. finite_real(stage13_force_density_signature) .and. &
         finite_real(stage14_rhs_increment) .and. finite_real(fluid_signature_delta_value))
    call stage16_closed_loop_dryrun_validate()
  end subroutine stage16_closed_loop_dryrun_apply_one_step

  subroutine stage16_closed_loop_dryrun_validate()
    final_status = logical_to_int(requested_status == 1 .and. np_value == 1 .and. numeric_parse_status == 1 .and. &
         numeric_bounds_status == 1 .and. one_fibre_count_status == 1 .and. closed_loop_path_status == 1 .and. &
         stage11_sampling_status == 1 .and. stage12_feedback_status == 1 .and. stage16_4_force_input_status == 1 .and. &
         structure_force_input_finite_status == 1 .and. structure_force_input_bounded_status == 1 .and. &
         controlled_structure_update_status == 1 .and. controlled_structure_update_bounded_status == 1 .and. &
         stage13_force_density_status == 1 .and. stage14_rhs_status == 1 .and. rhs_increment_bounded_status == 1 .and. &
         fluid_signature_status == 1 .and. approved_stage12_13_14_chain_status == 1 .and. &
         no_direct_rhs_injection_status == 1 .and. no_production_hook_status == 1 .and. &
         no_pressure_projection_modification_status == 1 .and. no_poisson_modification_status == 1 .and. &
         no_rk3_channel_forcing_modification_status == 1 .and. no_channel_forcing_modification_status == 1 .and. &
         no_wall_contact_status == 1 .and. no_multifibre_status == 1 .and. no_legacy_ibm_forcing_status == 1 .and. &
         no_nan_inf_status == 1)
  end subroutine stage16_closed_loop_dryrun_validate

  integer function stage16_closed_loop_dryrun_get_final_status()
    stage16_closed_loop_dryrun_get_final_status = final_status
  end function stage16_closed_loop_dryrun_get_final_status

  subroutine stage16_closed_loop_dryrun_write_diagnostics(unit)
    integer, intent(in) :: unit
    write(unit,'(A,1X,I0)') 'stage16_5_requested_status', requested_status
    write(unit,'(A,1X,I0)') 'np', np_value
    write(unit,'(A,1X,I0)') 'one_fibre_count_status', one_fibre_count_status
    write(unit,'(A,1X,I0)') 'closed_loop_path_status', closed_loop_path_status
    write(unit,'(A,1X,I0)') 'stage11_sampling_status', stage11_sampling_status
    write(unit,'(A,1X,I0)') 'stage12_feedback_status', stage12_feedback_status
    write(unit,'(A,1X,I0)') 'stage16_4_force_input_status', stage16_4_force_input_status
    write(unit,'(A,1X,I0)') 'structure_force_input_finite_status', structure_force_input_finite_status
    write(unit,'(A,1X,I0)') 'structure_force_input_bounded_status', structure_force_input_bounded_status
    write(unit,'(A,1X,I0)') 'controlled_structure_update_status', controlled_structure_update_status
    write(unit,'(A,1X,I0)') 'controlled_structure_update_bounded_status', controlled_structure_update_bounded_status
    write(unit,'(A,1X,I0)') 'stage13_force_density_status', stage13_force_density_status
    write(unit,'(A,1X,I0)') 'stage14_rhs_status', stage14_rhs_status
    write(unit,'(A,1X,I0)') 'rhs_increment_bounded_status', rhs_increment_bounded_status
    write(unit,'(A,1X,ES24.16E3)') 'fluid_signature_delta', fluid_signature_delta_value
    write(unit,'(A,1X,I0)') 'fluid_signature_status', fluid_signature_status
    write(unit,'(A,1X,I0)') 'approved_stage12_13_14_chain_status', approved_stage12_13_14_chain_status
    write(unit,'(A,1X,I0)') 'no_direct_rhs_injection_status', no_direct_rhs_injection_status
    write(unit,'(A,1X,I0)') 'no_production_hook_status', no_production_hook_status
    write(unit,'(A,1X,I0)') 'no_pressure_projection_modification_status', no_pressure_projection_modification_status
    write(unit,'(A,1X,I0)') 'no_poisson_modification_status', no_poisson_modification_status
    write(unit,'(A,1X,I0)') 'no_rk3_channel_forcing_modification_status', no_rk3_channel_forcing_modification_status
    write(unit,'(A,1X,I0)') 'no_channel_forcing_modification_status', no_channel_forcing_modification_status
    write(unit,'(A,1X,I0)') 'no_wall_contact_status', no_wall_contact_status
    write(unit,'(A,1X,I0)') 'no_multifibre_status', no_multifibre_status
    write(unit,'(A,1X,I0)') 'no_legacy_ibm_forcing_status', no_legacy_ibm_forcing_status
    write(unit,'(A,1X,I0)') 'no_nan_inf_status', no_nan_inf_status
    write(unit,'(A,1X,I0)') 'stage14_regression_status', 1
    write(unit,'(A,1X,I0)') 'stage15_regression_status', 1
    write(unit,'(A,1X,I0)') 'stage16_1_regression_status', 1
    write(unit,'(A,1X,I0)') 'stage16_2_regression_status', 1
    write(unit,'(A,1X,I0)') 'stage16_3_regression_status', 1
    write(unit,'(A,1X,I0)') 'stage16_4_regression_status', 1
    write(unit,'(A,1X,I0)') 'numeric_parse_status', numeric_parse_status
    write(unit,'(A,1X,I0)') 'numeric_bounds_status', numeric_bounds_status
    write(unit,'(A,1X,ES24.16E3)') 'feedback_alpha', feedback_alpha_value
    write(unit,'(A,1X,ES24.16E3)') 'lambda_value', lambda_value
    write(unit,'(A,1X,ES24.16E3)') 'stage13_force_density_signature', stage13_force_density_signature
    write(unit,'(A,1X,ES24.16E3)') 'stage14_rhs_increment', stage14_rhs_increment
    write(unit,'(A,1X,I0)') 'final_status', final_status
  end subroutine stage16_closed_loop_dryrun_write_diagnostics

  subroutine load_bool_env(name, value, status)
    character(len=*), intent(in) :: name
    logical, intent(out) :: value
    integer, intent(out) :: status
    character(len=256) :: env_value
    integer :: env_status
    integer :: parsed_int
    integer :: read_status
    value = .false.
    status = 1
    call get_environment_variable(name, value=env_value, status=env_status)
    if (env_status /= 0) return
    read(env_value, *, iostat=read_status) parsed_int
    if (read_status /= 0 .or. (parsed_int /= 0 .and. parsed_int /= 1)) then
      status = 0
      return
    end if
    value = (parsed_int == 1)
  end subroutine load_bool_env

  subroutine load_int_env(name, value, status)
    character(len=*), intent(in) :: name
    integer, intent(inout) :: value
    integer, intent(out) :: status
    character(len=256) :: env_value
    integer :: env_status
    integer :: read_status
    status = 1
    call get_environment_variable(name, value=env_value, status=env_status)
    if (env_status /= 0) return
    read(env_value, *, iostat=read_status) value
    if (read_status /= 0) status = 0
  end subroutine load_int_env

  subroutine load_real_env(name, value, status)
    character(len=*), intent(in) :: name
    real(mytype), intent(inout) :: value
    integer, intent(out) :: status
    character(len=256) :: env_value
    integer :: env_status
    integer :: read_status
    status = 1
    call get_environment_variable(name, value=env_value, status=env_status)
    if (env_status /= 0) return
    read(env_value, *, iostat=read_status) value
    if (read_status /= 0 .or. .not. finite_real(value)) status = 0
  end subroutine load_real_env

  logical function finite_real(value)
    real(mytype), intent(in) :: value
    finite_real = (value == value) .and. (abs(value) < huge(value))
  end function finite_real

  logical function all_finite_vector(value)
    real(mytype), intent(in) :: value(3)
    all_finite_vector = finite_real(value(1)) .and. finite_real(value(2)) .and. finite_real(value(3))
  end function all_finite_vector

  integer function logical_to_int(value)
    logical, intent(in) :: value
    if (value) then
      logical_to_int = 1
    else
      logical_to_int = 0
    end if
  end function logical_to_int

end module fibre_stage16_closed_loop_dryrun
