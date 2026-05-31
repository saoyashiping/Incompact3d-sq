module fibre_stage13_volume_normalization_audit
  use decomp_2d_constants, only : mytype
  implicit none
  private

  integer :: initialized_status = 0
  integer :: positive_volume_status = 0
  integer :: finite_volume_status = 0
  integer :: invalid_zero_volume_rejection_status = 0
  integer :: invalid_negative_volume_rejection_status = 0
  integer :: uniform_volume_conservation_status = 0
  integer :: nonuniform_volume_conservation_status = 0
  integer :: density_times_volume_integral_status = 0
  integer :: component_volume_normalization_status = 0
  integer :: boundary_volume_normalization_status = 0
  integer :: finite_force_density_status = 0
  integer :: clear_status = 0
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: volume_normalization_audit_status = 0

  real(mytype) :: last_min_volume = 0.0_mytype
  real(mytype) :: last_max_volume = 0.0_mytype

  public :: stage13_volume_normalization_audit_init
  public :: stage13_volume_normalization_audit_clear
  public :: stage13_volume_normalization_audit_check_volumes
  public :: stage13_volume_normalization_audit_integrate_force_density
  public :: stage13_volume_normalization_audit_record_statuses
  public :: stage13_volume_normalization_audit_get_status_values
  public :: stage13_volume_normalization_audit_write_diagnostics
  public :: stage13_volume_normalization_audit_finalize

contains

  subroutine stage13_volume_normalization_audit_init()
    initialized_status = 1
    positive_volume_status = 0
    finite_volume_status = 0
    invalid_zero_volume_rejection_status = 0
    invalid_negative_volume_rejection_status = 0
    uniform_volume_conservation_status = 0
    nonuniform_volume_conservation_status = 0
    density_times_volume_integral_status = 0
    component_volume_normalization_status = 0
    boundary_volume_normalization_status = 0
    finite_force_density_status = 0
    clear_status = 0
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    no_fluid_field_access_status = 1
    no_fluid_field_modification_status = 1
    last_min_volume = 0.0_mytype
    last_max_volume = 0.0_mytype
    call update_status()
  end subroutine stage13_volume_normalization_audit_init

  subroutine stage13_volume_normalization_audit_clear(clear_ok)
    logical, intent(in), optional :: clear_ok

    if (present(clear_ok)) then
      clear_status = merge(1, 0, clear_ok)
    else
      clear_status = 1
    end if
    call update_status()
  end subroutine stage13_volume_normalization_audit_clear

  subroutine stage13_volume_normalization_audit_check_volumes(volume_eul, min_volume, max_volume)
    real(mytype), intent(in) :: volume_eul(:,:,:)
    real(mytype), intent(out) :: min_volume
    real(mytype), intent(out) :: max_volume
    logical :: finite_ok
    logical :: positive_ok

    min_volume = minval(volume_eul)
    max_volume = maxval(volume_eul)
    finite_ok = all_finite_rank3(volume_eul)
    positive_ok = finite_ok .and. all(volume_eul > 0.0_mytype)

    finite_volume_status = merge(1, 0, finite_ok)
    positive_volume_status = merge(1, 0, positive_ok)
    if (finite_ok .and. any(volume_eul == 0.0_mytype)) invalid_zero_volume_rejection_status = 1
    if (finite_ok .and. any(volume_eul < 0.0_mytype)) invalid_negative_volume_rejection_status = 1
    last_min_volume = min_volume
    last_max_volume = max_volume
    call update_status()
  end subroutine stage13_volume_normalization_audit_check_volumes

  subroutine stage13_volume_normalization_audit_integrate_force_density(fx_cand, fy_cand, fz_cand, volume_eul, &
                                                                       integrated_force)
    real(mytype), intent(in) :: fx_cand(:,:,:)
    real(mytype), intent(in) :: fy_cand(:,:,:)
    real(mytype), intent(in) :: fz_cand(:,:,:)
    real(mytype), intent(in) :: volume_eul(:,:,:)
    real(mytype), intent(out) :: integrated_force(3)

    if (.not. shapes_match(fx_cand, fy_cand, fz_cand, volume_eul)) then
      integrated_force(:) = 0.0_mytype
      finite_force_density_status = 0
      call update_status()
      return
    end if

    integrated_force(1) = sum(fx_cand * volume_eul)
    integrated_force(2) = sum(fy_cand * volume_eul)
    integrated_force(3) = sum(fz_cand * volume_eul)
    finite_force_density_status = merge(1, 0, all_finite_rank3(fx_cand) .and. all_finite_rank3(fy_cand) .and. &
                                             all_finite_rank3(fz_cand))
    call update_status()
  end subroutine stage13_volume_normalization_audit_integrate_force_density

  subroutine stage13_volume_normalization_audit_record_statuses(uniform_status, nonuniform_status, integral_status, &
                                                               component_status, boundary_status, finite_status, &
                                                               clear_status_in)
    integer, intent(in) :: uniform_status
    integer, intent(in) :: nonuniform_status
    integer, intent(in) :: integral_status
    integer, intent(in) :: component_status
    integer, intent(in) :: boundary_status
    integer, intent(in) :: finite_status
    integer, intent(in) :: clear_status_in

    uniform_volume_conservation_status = uniform_status
    nonuniform_volume_conservation_status = nonuniform_status
    density_times_volume_integral_status = integral_status
    component_volume_normalization_status = component_status
    boundary_volume_normalization_status = boundary_status
    finite_force_density_status = min(finite_force_density_status, finite_status)
    clear_status = clear_status_in
    call update_status()
  end subroutine stage13_volume_normalization_audit_record_statuses

  subroutine stage13_volume_normalization_audit_get_status_values(initialized_out, positive_out, finite_volume_out, &
                                                                  invalid_zero_out, invalid_negative_out, &
                                                                  uniform_out, nonuniform_out, density_integral_out, &
                                                                  component_out, boundary_out, finite_force_out, &
                                                                  clear_out, no_rhs_out, no_ibm_out, no_feedback_out, &
                                                                  no_twoway_out, no_structure_out, no_fluid_access_out, &
                                                                  no_fluid_modification_out, audit_status_out)
    integer, intent(out) :: initialized_out
    integer, intent(out) :: positive_out
    integer, intent(out) :: finite_volume_out
    integer, intent(out) :: invalid_zero_out
    integer, intent(out) :: invalid_negative_out
    integer, intent(out) :: uniform_out
    integer, intent(out) :: nonuniform_out
    integer, intent(out) :: density_integral_out
    integer, intent(out) :: component_out
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
    positive_out = positive_volume_status
    finite_volume_out = finite_volume_status
    invalid_zero_out = invalid_zero_volume_rejection_status
    invalid_negative_out = invalid_negative_volume_rejection_status
    uniform_out = uniform_volume_conservation_status
    nonuniform_out = nonuniform_volume_conservation_status
    density_integral_out = density_times_volume_integral_status
    component_out = component_volume_normalization_status
    boundary_out = boundary_volume_normalization_status
    finite_force_out = finite_force_density_status
    clear_out = clear_status
    no_rhs_out = no_rhs_injection_status
    no_ibm_out = no_ibm_spreading_status
    no_feedback_out = no_feedback_application_status
    no_twoway_out = no_twoway_force_status
    no_structure_out = no_structure_advance_status
    no_fluid_access_out = no_fluid_field_access_status
    no_fluid_modification_out = no_fluid_field_modification_status
    audit_status_out = volume_normalization_audit_status
  end subroutine stage13_volume_normalization_audit_get_status_values

  subroutine stage13_volume_normalization_audit_write_diagnostics()
    integer :: unit_id
    integer :: io_status

    open(newunit=unit_id, file='stage13_outputs/fibre_stage13_volume_normalization_audit.dat', &
         status='replace', action='write', iostat=io_status)
    if (io_status /= 0) return
    call update_status()
    write(unit_id,'(A,1X,I0)') 'stage13_volume_normalization_positive_volume_status', positive_volume_status
    write(unit_id,'(A,1X,I0)') 'stage13_volume_normalization_finite_volume_status', finite_volume_status
    write(unit_id,'(A,1X,I0)') 'stage13_volume_normalization_audit_status', volume_normalization_audit_status
    write(unit_id,'(A,1X,ES24.16)') 'stage13_volume_normalization_min_volume', last_min_volume
    write(unit_id,'(A,1X,ES24.16)') 'stage13_volume_normalization_max_volume', last_max_volume
    close(unit_id)
  end subroutine stage13_volume_normalization_audit_write_diagnostics

  subroutine stage13_volume_normalization_audit_finalize()
    initialized_status = 0
    call update_status()
  end subroutine stage13_volume_normalization_audit_finalize

  subroutine update_status()
    if (initialized_status == 1 .and. positive_volume_status == 1 .and. finite_volume_status == 1 .and. &
        invalid_zero_volume_rejection_status == 1 .and. invalid_negative_volume_rejection_status == 1 .and. &
        uniform_volume_conservation_status == 1 .and. nonuniform_volume_conservation_status == 1 .and. &
        density_times_volume_integral_status == 1 .and. component_volume_normalization_status == 1 .and. &
        boundary_volume_normalization_status == 1 .and. finite_force_density_status == 1 .and. &
        clear_status == 1 .and. no_rhs_injection_status == 1 .and. no_ibm_spreading_status == 1 .and. &
        no_feedback_application_status == 1 .and. no_twoway_force_status == 1 .and. &
        no_structure_advance_status == 1 .and. no_fluid_field_access_status == 1 .and. &
        no_fluid_field_modification_status == 1) then
      volume_normalization_audit_status = 1
    else
      volume_normalization_audit_status = 0
    end if
  end subroutine update_status

  logical function shapes_match(fx_cand, fy_cand, fz_cand, volume_eul)
    real(mytype), intent(in) :: fx_cand(:,:,:)
    real(mytype), intent(in) :: fy_cand(:,:,:)
    real(mytype), intent(in) :: fz_cand(:,:,:)
    real(mytype), intent(in) :: volume_eul(:,:,:)

    shapes_match = all(shape(fx_cand) == shape(volume_eul)) .and. all(shape(fy_cand) == shape(volume_eul)) .and. &
                   all(shape(fz_cand) == shape(volume_eul))
  end function shapes_match

  logical function all_finite_rank3(values)
    real(mytype), intent(in) :: values(:,:,:)

    all_finite_rank3 = all(values == values) .and. all(abs(values) < huge(1.0_mytype))
  end function all_finite_rank3

end module fibre_stage13_volume_normalization_audit
