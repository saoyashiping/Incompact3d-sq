module fibre_stage11_oneway_interpolation
  implicit none
  private

  logical :: initialized = .false.
  integer :: prepare_called_status = 0
  integer :: sample_interface_called_status = 0
  integer :: sample_performed_status = 0
  integer :: interface_available_status = 1
  integer :: lagrangian_state_input_status = 0
  integer :: grid_metadata_input_status = 0
  integer :: velocity_placeholder_unchanged_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_velocity_sampling_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_force_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: interpolation_api_status = 0

  public :: stage11_oneway_interpolation_init
  public :: stage11_oneway_interpolation_finalize
  public :: stage11_oneway_interpolation_is_initialized
  public :: stage11_oneway_interpolation_prepare
  public :: stage11_oneway_interpolation_sample_velocity
  public :: stage11_oneway_interpolation_get_status_values
  public :: stage11_oneway_interpolation_write_diagnostics

contains

  subroutine stage11_oneway_interpolation_init()
    initialized = .true.
    prepare_called_status = 0
    sample_interface_called_status = 0
    sample_performed_status = 0
    interface_available_status = 1
    lagrangian_state_input_status = 0
    grid_metadata_input_status = 0
    velocity_placeholder_unchanged_status = 1
    no_fluid_field_access_status = 1
    no_velocity_sampling_status = 1
    no_fluid_field_modification_status = 1
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_force_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    call update_status()
  end subroutine

  subroutine stage11_oneway_interpolation_finalize()
    initialized = .false.
    prepare_called_status = 0
    sample_interface_called_status = 0
    sample_performed_status = 0
    interface_available_status = 1
    lagrangian_state_input_status = 0
    grid_metadata_input_status = 0
    velocity_placeholder_unchanged_status = 1
    no_fluid_field_access_status = 1
    no_velocity_sampling_status = 1
    no_fluid_field_modification_status = 1
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_force_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    interpolation_api_status = 0
  end subroutine

  logical function stage11_oneway_interpolation_is_initialized()
    stage11_oneway_interpolation_is_initialized = initialized
  end function

  subroutine stage11_oneway_interpolation_prepare(lagrangian_state_ready, grid_metadata_ready)
    logical, intent(in) :: lagrangian_state_ready, grid_metadata_ready
    prepare_called_status = 1
    lagrangian_state_input_status = 0
    grid_metadata_input_status = 0
    if (lagrangian_state_ready) lagrangian_state_input_status = 1
    if (grid_metadata_ready) grid_metadata_input_status = 1
    call update_status()
  end subroutine

  subroutine stage11_oneway_interpolation_sample_velocity()
    sample_interface_called_status = 1
    sample_performed_status = 0
    no_velocity_sampling_status = 1
    no_fluid_field_access_status = 1
    velocity_placeholder_unchanged_status = 1
    call update_status()
  end subroutine

  subroutine stage11_oneway_interpolation_write_diagnostics(unit)
    integer, intent(in) :: unit
    write(unit,'(A,1X,I0)') 'stage11_3_interpolation_initialized_status', merge(1,0,initialized)
    write(unit,'(A,1X,I0)') 'stage11_3_interface_available_status', interface_available_status
    write(unit,'(A,1X,I0)') 'stage11_3_prepare_called_status', prepare_called_status
    write(unit,'(A,1X,I0)') 'stage11_3_sample_interface_called_status', sample_interface_called_status
    write(unit,'(A,1X,I0)') 'stage11_3_sample_not_performed_status', merge(1,0,sample_performed_status==0)
    write(unit,'(A,1X,I0)') 'stage11_3_lagrangian_state_input_status', lagrangian_state_input_status
    write(unit,'(A,1X,I0)') 'stage11_3_grid_metadata_input_status', grid_metadata_input_status
    write(unit,'(A,1X,I0)') 'stage11_3_velocity_placeholder_unchanged_status', velocity_placeholder_unchanged_status
    write(unit,'(A,1X,I0)') 'stage11_3_no_fluid_field_access_status', no_fluid_field_access_status
    write(unit,'(A,1X,I0)') 'stage11_3_no_velocity_sampling_status', no_velocity_sampling_status
    write(unit,'(A,1X,I0)') 'stage11_3_no_fluid_field_modification_status', no_fluid_field_modification_status
    write(unit,'(A,1X,I0)') 'stage11_3_no_rhs_injection_status', no_rhs_injection_status
    write(unit,'(A,1X,I0)') 'stage11_3_no_ibm_spreading_status', no_ibm_spreading_status
    write(unit,'(A,1X,I0)') 'stage11_3_no_feedback_force_status', no_feedback_force_status
    write(unit,'(A,1X,I0)') 'stage11_3_no_twoway_force_status', no_twoway_force_status
    write(unit,'(A,1X,I0)') 'stage11_3_no_structure_advance_status', no_structure_advance_status
    write(unit,'(A,1X,I0)') 'stage11_3_oneway_interpolation_api_status', interpolation_api_status
  end subroutine

  subroutine stage11_oneway_interpolation_get_status_values(initialized_status, interface_status, prepare_status, sample_called_status, &
                                                            sample_not_performed_status, lagrangian_input_status, grid_input_status, &
                                                            velocity_unchanged_status, no_fluid_access, no_velocity_sampling, &
                                                            no_fluid_modification, no_rhs_injection, no_ibm_spreading, &
                                                            no_feedback_force, no_twoway_force, no_structure_advance, api_status)
    integer, intent(out) :: initialized_status, interface_status, prepare_status, sample_called_status
    integer, intent(out) :: sample_not_performed_status, lagrangian_input_status, grid_input_status
    integer, intent(out) :: velocity_unchanged_status, no_fluid_access, no_velocity_sampling
    integer, intent(out) :: no_fluid_modification, no_rhs_injection, no_ibm_spreading
    integer, intent(out) :: no_feedback_force, no_twoway_force, no_structure_advance, api_status

    initialized_status = merge(1,0,initialized)
    interface_status = interface_available_status
    prepare_status = prepare_called_status
    sample_called_status = sample_interface_called_status
    sample_not_performed_status = merge(1,0,sample_performed_status == 0)
    lagrangian_input_status = lagrangian_state_input_status
    grid_input_status = grid_metadata_input_status
    velocity_unchanged_status = velocity_placeholder_unchanged_status
    no_fluid_access = no_fluid_field_access_status
    no_velocity_sampling = no_velocity_sampling_status
    no_fluid_modification = no_fluid_field_modification_status
    no_rhs_injection = no_rhs_injection_status
    no_ibm_spreading = no_ibm_spreading_status
    no_feedback_force = no_feedback_force_status
    no_twoway_force = no_twoway_force_status
    no_structure_advance = no_structure_advance_status
    api_status = interpolation_api_status
  end subroutine

  subroutine update_status()
    interpolation_api_status = 1
    if (.not. initialized) interpolation_api_status = 0
    if (interface_available_status /= 1) interpolation_api_status = 0
    if (prepare_called_status /= 1) interpolation_api_status = 0
    if (sample_interface_called_status /= 1) interpolation_api_status = 0
    if (sample_performed_status /= 0) interpolation_api_status = 0
    if (lagrangian_state_input_status /= 1) interpolation_api_status = 0
    if (grid_metadata_input_status /= 1) interpolation_api_status = 0
    if (velocity_placeholder_unchanged_status /= 1) interpolation_api_status = 0
    if (no_fluid_field_access_status /= 1) interpolation_api_status = 0
    if (no_velocity_sampling_status /= 1) interpolation_api_status = 0
    if (no_fluid_field_modification_status /= 1) interpolation_api_status = 0
    if (no_rhs_injection_status /= 1) interpolation_api_status = 0
    if (no_ibm_spreading_status /= 1) interpolation_api_status = 0
    if (no_feedback_force_status /= 1) interpolation_api_status = 0
    if (no_twoway_force_status /= 1) interpolation_api_status = 0
    if (no_structure_advance_status /= 1) interpolation_api_status = 0
  end subroutine

end module fibre_stage11_oneway_interpolation
