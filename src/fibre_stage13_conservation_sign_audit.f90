module fibre_stage13_conservation_sign_audit
  use decomp_2d_constants, only : mytype
  implicit none
  private

  real(mytype), parameter :: action_reaction_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: force_conservation_abs_tol = 1.0e-12_mytype
  real(mytype), parameter :: wrong_sign_min_separation = 1.0e-8_mytype

  integer :: initialized_status = 0
  integer :: fluid_to_fibre_sign_status = 0
  integer :: fibre_to_fluid_sign_status = 0
  integer :: action_reaction_status = 0
  integer :: spreading_input_sign_status = 0
  integer :: integrated_force_fibre_to_fluid_status = 0
  integer :: wrong_sign_rejection_status = 0
  integer :: component_sign_conservation_status = 0
  integer :: multipoint_sign_conservation_status = 0
  integer :: nonuniform_volume_sign_conservation_status = 0
  integer :: boundary_sign_conservation_status = 0
  integer :: finite_force_density_status = 0
  integer :: clear_status = 0
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: conservation_sign_audit_status = 0

  real(mytype) :: action_reaction_max_error = 0.0_mytype
  real(mytype) :: integrated_force_fibre_to_fluid_error_l2 = 0.0_mytype
  real(mytype) :: wrong_sign_error_l2 = 0.0_mytype
  real(mytype) :: wrong_sign_separation = 0.0_mytype

  public :: stage13_conservation_sign_audit_init
  public :: stage13_conservation_sign_audit_clear
  public :: stage13_conservation_sign_audit_compute_expected_forces
  public :: stage13_conservation_sign_audit_check_action_reaction
  public :: stage13_conservation_sign_audit_check_integrated_force_sign
  public :: stage13_conservation_sign_audit_record_statuses
  public :: stage13_conservation_sign_audit_get_status_values
  public :: stage13_conservation_sign_audit_get_diagnostics
  public :: stage13_conservation_sign_audit_write_diagnostics
  public :: stage13_conservation_sign_audit_finalize

contains

  subroutine stage13_conservation_sign_audit_init()
    initialized_status = 1
    fluid_to_fibre_sign_status = 0
    fibre_to_fluid_sign_status = 0
    action_reaction_status = 0
    spreading_input_sign_status = 0
    integrated_force_fibre_to_fluid_status = 0
    wrong_sign_rejection_status = 0
    component_sign_conservation_status = 0
    multipoint_sign_conservation_status = 0
    nonuniform_volume_sign_conservation_status = 0
    boundary_sign_conservation_status = 0
    finite_force_density_status = 0
    clear_status = 0
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    no_fluid_field_access_status = 1
    no_fluid_field_modification_status = 1
    action_reaction_max_error = 0.0_mytype
    integrated_force_fibre_to_fluid_error_l2 = 0.0_mytype
    wrong_sign_error_l2 = 0.0_mytype
    wrong_sign_separation = 0.0_mytype
    call update_status()
  end subroutine stage13_conservation_sign_audit_init

  subroutine stage13_conservation_sign_audit_clear(clear_ok)
    logical, intent(in), optional :: clear_ok

    if (present(clear_ok)) then
      clear_status = merge(1, 0, clear_ok)
    else
      clear_status = 1
    end if
    call update_status()
  end subroutine stage13_conservation_sign_audit_clear

  subroutine stage13_conservation_sign_audit_compute_expected_forces(f_fluid_to_fibre, f_fibre_to_fluid, ds_lag, &
                                                                    expected_fluid_to_fibre, &
                                                                    expected_fibre_to_fluid)
    real(mytype), intent(in) :: f_fluid_to_fibre(:,:)
    real(mytype), intent(in) :: f_fibre_to_fluid(:,:)
    real(mytype), intent(in) :: ds_lag(:)
    real(mytype), intent(out) :: expected_fluid_to_fibre(3)
    real(mytype), intent(out) :: expected_fibre_to_fluid(3)
    integer :: p

    expected_fluid_to_fibre(:) = 0.0_mytype
    expected_fibre_to_fluid(:) = 0.0_mytype
    if (.not. force_shapes_valid(f_fluid_to_fibre, f_fibre_to_fluid, ds_lag)) then
      fluid_to_fibre_sign_status = 0
      fibre_to_fluid_sign_status = 0
      call update_status()
      return
    end if

    do p = 1, size(ds_lag)
      expected_fluid_to_fibre(1) = expected_fluid_to_fibre(1) + f_fluid_to_fibre(p,1) * ds_lag(p)
      expected_fluid_to_fibre(2) = expected_fluid_to_fibre(2) + f_fluid_to_fibre(p,2) * ds_lag(p)
      expected_fluid_to_fibre(3) = expected_fluid_to_fibre(3) + f_fluid_to_fibre(p,3) * ds_lag(p)
      expected_fibre_to_fluid(1) = expected_fibre_to_fluid(1) + f_fibre_to_fluid(p,1) * ds_lag(p)
      expected_fibre_to_fluid(2) = expected_fibre_to_fluid(2) + f_fibre_to_fluid(p,2) * ds_lag(p)
      expected_fibre_to_fluid(3) = expected_fibre_to_fluid(3) + f_fibre_to_fluid(p,3) * ds_lag(p)
    end do

    fluid_to_fibre_sign_status = merge(1, 0, all_finite_vector(expected_fluid_to_fibre))
    fibre_to_fluid_sign_status = merge(1, 0, vector_l2(expected_fibre_to_fluid + expected_fluid_to_fibre) <= &
                                             force_conservation_abs_tol)
    call update_status()
  end subroutine stage13_conservation_sign_audit_compute_expected_forces

  subroutine stage13_conservation_sign_audit_check_action_reaction(f_fluid_to_fibre, f_fibre_to_fluid, &
                                                                  max_action_reaction_error)
    real(mytype), intent(in) :: f_fluid_to_fibre(:,:)
    real(mytype), intent(in) :: f_fibre_to_fluid(:,:)
    real(mytype), intent(out) :: max_action_reaction_error

    if (.not. force_shapes_valid_no_ds(f_fluid_to_fibre, f_fibre_to_fluid)) then
      max_action_reaction_error = huge(1.0_mytype)
      action_reaction_status = 0
      call update_status()
      return
    end if

    max_action_reaction_error = maxval(abs(f_fluid_to_fibre(:,1:3) + f_fibre_to_fluid(:,1:3)))
    action_reaction_max_error = max_action_reaction_error
    action_reaction_status = merge(1, 0, max_action_reaction_error <= action_reaction_abs_tol)
    call update_status()
  end subroutine stage13_conservation_sign_audit_check_action_reaction

  subroutine stage13_conservation_sign_audit_check_integrated_force_sign(integrated_eulerian_force, &
                                                                        expected_fluid_to_fibre, &
                                                                        expected_fibre_to_fluid, &
                                                                        fibre_to_fluid_error, wrong_sign_error)
    real(mytype), intent(in) :: integrated_eulerian_force(3)
    real(mytype), intent(in) :: expected_fluid_to_fibre(3)
    real(mytype), intent(in) :: expected_fibre_to_fluid(3)
    real(mytype), intent(out) :: fibre_to_fluid_error
    real(mytype), intent(out) :: wrong_sign_error

    fibre_to_fluid_error = vector_l2(integrated_eulerian_force - expected_fibre_to_fluid)
    wrong_sign_error = vector_l2(integrated_eulerian_force - expected_fluid_to_fibre)
    integrated_force_fibre_to_fluid_error_l2 = fibre_to_fluid_error
    wrong_sign_error_l2 = wrong_sign_error
    wrong_sign_separation = wrong_sign_error - fibre_to_fluid_error
    integrated_force_fibre_to_fluid_status = merge(1, 0, fibre_to_fluid_error <= force_conservation_abs_tol)
    wrong_sign_rejection_status = merge(1, 0, integrated_force_fibre_to_fluid_status == 1 .and. &
                                             wrong_sign_separation >= wrong_sign_min_separation)
    call update_status()
  end subroutine stage13_conservation_sign_audit_check_integrated_force_sign

  subroutine stage13_conservation_sign_audit_record_statuses(spreading_input_status, component_status, &
                                                            multipoint_status, nonuniform_status, boundary_status, &
                                                            finite_status, clear_status_in)
    integer, intent(in) :: spreading_input_status
    integer, intent(in) :: component_status
    integer, intent(in) :: multipoint_status
    integer, intent(in) :: nonuniform_status
    integer, intent(in) :: boundary_status
    integer, intent(in) :: finite_status
    integer, intent(in) :: clear_status_in

    spreading_input_sign_status = spreading_input_status
    component_sign_conservation_status = component_status
    multipoint_sign_conservation_status = multipoint_status
    nonuniform_volume_sign_conservation_status = nonuniform_status
    boundary_sign_conservation_status = boundary_status
    finite_force_density_status = finite_status
    clear_status = clear_status_in
    call update_status()
  end subroutine stage13_conservation_sign_audit_record_statuses

  subroutine stage13_conservation_sign_audit_get_status_values(initialized_out, fluid_to_fibre_out, &
                                                              fibre_to_fluid_out, action_reaction_out, &
                                                              spreading_input_out, integrated_force_out, &
                                                              wrong_sign_out, component_out, multipoint_out, &
                                                              nonuniform_out, boundary_out, finite_force_out, &
                                                              clear_out, no_rhs_out, no_ibm_out, no_feedback_out, &
                                                              no_twoway_out, no_structure_out, no_fluid_access_out, &
                                                              no_fluid_modification_out, audit_status_out)
    integer, intent(out) :: initialized_out
    integer, intent(out) :: fluid_to_fibre_out
    integer, intent(out) :: fibre_to_fluid_out
    integer, intent(out) :: action_reaction_out
    integer, intent(out) :: spreading_input_out
    integer, intent(out) :: integrated_force_out
    integer, intent(out) :: wrong_sign_out
    integer, intent(out) :: component_out
    integer, intent(out) :: multipoint_out
    integer, intent(out) :: nonuniform_out
    integer, intent(out) :: boundary_out
    integer, intent(out) :: finite_force_out
    integer, intent(out) :: clear_out
    integer, intent(out) :: no_rhs_out
    integer, intent(out) :: no_ibm_out
    integer, intent(out) :: no_feedback_out
    integer, intent(out) :: no_twoway_out
    integer, intent(out) :: no_structure_out
    integer, intent(out) :: no_fluid_access_out
    integer, intent(out) :: no_fluid_modification_out
    integer, intent(out) :: audit_status_out

    call update_status()
    initialized_out = initialized_status
    fluid_to_fibre_out = fluid_to_fibre_sign_status
    fibre_to_fluid_out = fibre_to_fluid_sign_status
    action_reaction_out = action_reaction_status
    spreading_input_out = spreading_input_sign_status
    integrated_force_out = integrated_force_fibre_to_fluid_status
    wrong_sign_out = wrong_sign_rejection_status
    component_out = component_sign_conservation_status
    multipoint_out = multipoint_sign_conservation_status
    nonuniform_out = nonuniform_volume_sign_conservation_status
    boundary_out = boundary_sign_conservation_status
    finite_force_out = finite_force_density_status
    clear_out = clear_status
    no_rhs_out = no_rhs_injection_status
    no_ibm_out = no_ibm_spreading_status
    no_feedback_out = no_feedback_application_status
    no_twoway_out = no_twoway_force_status
    no_structure_out = no_structure_advance_status
    no_fluid_access_out = no_fluid_field_access_status
    no_fluid_modification_out = no_fluid_field_modification_status
    audit_status_out = conservation_sign_audit_status
  end subroutine stage13_conservation_sign_audit_get_status_values

  subroutine stage13_conservation_sign_audit_get_diagnostics(action_error_out, fibre_to_fluid_error_out, &
                                                            wrong_sign_error_out, wrong_sign_separation_out)
    real(mytype), intent(out) :: action_error_out
    real(mytype), intent(out) :: fibre_to_fluid_error_out
    real(mytype), intent(out) :: wrong_sign_error_out
    real(mytype), intent(out) :: wrong_sign_separation_out

    action_error_out = action_reaction_max_error
    fibre_to_fluid_error_out = integrated_force_fibre_to_fluid_error_l2
    wrong_sign_error_out = wrong_sign_error_l2
    wrong_sign_separation_out = wrong_sign_separation
  end subroutine stage13_conservation_sign_audit_get_diagnostics

  subroutine stage13_conservation_sign_audit_write_diagnostics()
    integer :: unit_id
    integer :: io_status

    open(newunit=unit_id, file='stage13_outputs/fibre_stage13_conservation_sign_audit.dat', &
         status='replace', action='write', iostat=io_status)
    if (io_status /= 0) return
    call update_status()
    write(unit_id,'(A,1X,I0)') 'stage13_conservation_sign_action_reaction_status', action_reaction_status
    write(unit_id,'(A,1X,I0)') 'stage13_conservation_sign_wrong_sign_rejection_status', wrong_sign_rejection_status
    write(unit_id,'(A,1X,I0)') 'stage13_conservation_sign_audit_status', conservation_sign_audit_status
    close(unit_id)
  end subroutine stage13_conservation_sign_audit_write_diagnostics

  subroutine stage13_conservation_sign_audit_finalize()
    initialized_status = 0
    call update_status()
  end subroutine stage13_conservation_sign_audit_finalize

  subroutine update_status()
    if (initialized_status == 1 .and. fluid_to_fibre_sign_status == 1 .and. fibre_to_fluid_sign_status == 1 .and. &
        action_reaction_status == 1 .and. spreading_input_sign_status == 1 .and. &
        integrated_force_fibre_to_fluid_status == 1 .and. wrong_sign_rejection_status == 1 .and. &
        component_sign_conservation_status == 1 .and. multipoint_sign_conservation_status == 1 .and. &
        nonuniform_volume_sign_conservation_status == 1 .and. boundary_sign_conservation_status == 1 .and. &
        finite_force_density_status == 1 .and. clear_status == 1 .and. no_rhs_injection_status == 1 .and. &
        no_ibm_spreading_status == 1 .and. no_feedback_application_status == 1 .and. &
        no_twoway_force_status == 1 .and. no_structure_advance_status == 1 .and. &
        no_fluid_field_access_status == 1 .and. no_fluid_field_modification_status == 1) then
      conservation_sign_audit_status = 1
    else
      conservation_sign_audit_status = 0
    end if
  end subroutine update_status

  logical function force_shapes_valid(f_fluid_to_fibre, f_fibre_to_fluid, ds_lag)
    real(mytype), intent(in) :: f_fluid_to_fibre(:,:)
    real(mytype), intent(in) :: f_fibre_to_fluid(:,:)
    real(mytype), intent(in) :: ds_lag(:)

    force_shapes_valid = force_shapes_valid_no_ds(f_fluid_to_fibre, f_fibre_to_fluid) .and. &
                         size(ds_lag) == size(f_fluid_to_fibre, 1)
  end function force_shapes_valid

  logical function force_shapes_valid_no_ds(f_fluid_to_fibre, f_fibre_to_fluid)
    real(mytype), intent(in) :: f_fluid_to_fibre(:,:)
    real(mytype), intent(in) :: f_fibre_to_fluid(:,:)

    force_shapes_valid_no_ds = size(f_fluid_to_fibre, 1) == size(f_fibre_to_fluid, 1) .and. &
                               size(f_fluid_to_fibre, 2) >= 3 .and. size(f_fibre_to_fluid, 2) >= 3
  end function force_shapes_valid_no_ds

  real(mytype) function vector_l2(values)
    real(mytype), intent(in) :: values(3)

    vector_l2 = sqrt(sum(values * values))
  end function vector_l2

  logical function all_finite_vector(values)
    real(mytype), intent(in) :: values(3)

    all_finite_vector = all(values == values) .and. all(abs(values) < huge(1.0_mytype))
  end function all_finite_vector

end module fibre_stage13_conservation_sign_audit
