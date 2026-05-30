module fibre_stage12_sign_convention_audit
  implicit none
  private

  integer, parameter :: mytype = kind(1.0d0)

  integer :: initialized_status = 0
  integer :: fluid_to_fibre_sign_status = 0
  integer :: fibre_to_fluid_sign_status = 0
  integer :: action_reaction_status = 0
  integer :: zero_slip_neutrality_status = 0
  integer :: positive_slip_status = 0
  integer :: negative_slip_status = 0
  integer :: structure_side_convention_status = 0
  integer :: fluid_side_convention_status = 0
  integer :: finite_force_status = 0
  integer :: force_norm_finite_status = 0
  integer :: no_eulerian_force_density_status = 1
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: sign_convention_audit_status = 0

  public :: stage12_sign_convention_audit_init
  public :: stage12_sign_convention_compute_pair
  public :: stage12_sign_convention_audit_clear
  public :: stage12_sign_convention_audit_finalize
  public :: stage12_sign_convention_audit_get_status_values
  public :: stage12_sign_convention_audit_write_diagnostics
  public :: stage12_sign_convention_audit_record_statuses

contains

  subroutine stage12_sign_convention_audit_init()
    initialized_status = 1
    fluid_to_fibre_sign_status = 0
    fibre_to_fluid_sign_status = 0
    action_reaction_status = 0
    zero_slip_neutrality_status = 0
    positive_slip_status = 0
    negative_slip_status = 0
    structure_side_convention_status = 0
    fluid_side_convention_status = 0
    finite_force_status = 0
    force_norm_finite_status = 0
    no_eulerian_force_density_status = 1
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    no_fluid_field_access_status = 1
    no_fluid_field_modification_status = 1
    call update_sign_convention_audit_status()
  end subroutine stage12_sign_convention_audit_init

  subroutine stage12_sign_convention_compute_pair(u_f, v_f, alpha, f_fluid_to_fibre, f_fibre_to_fluid, force_norm)
    real(mytype), intent(in) :: u_f(:,:)
    real(mytype), intent(in) :: v_f(:,:)
    real(mytype), intent(in) :: alpha
    real(mytype), intent(out) :: f_fluid_to_fibre(:,:)
    real(mytype), intent(out) :: f_fibre_to_fluid(:,:)
    real(mytype), intent(out) :: force_norm(:)
    integer :: n_points
    integer :: k
    integer :: c
    real(mytype) :: norm_sum

    finite_force_status = 0
    force_norm_finite_status = 0
    n_points = size(u_f, 1)
    if (size(u_f, 2) /= 3 .or. size(v_f, 1) /= n_points .or. size(v_f, 2) /= 3 .or. &
        size(f_fluid_to_fibre, 1) /= n_points .or. size(f_fluid_to_fibre, 2) /= 3 .or. &
        size(f_fibre_to_fluid, 1) /= n_points .or. size(f_fibre_to_fluid, 2) /= 3 .or. &
        size(force_norm) /= n_points) then
      f_fluid_to_fibre(:, :) = 0.0_mytype
      f_fibre_to_fluid(:, :) = 0.0_mytype
      force_norm(:) = 0.0_mytype
      call update_sign_convention_audit_status()
      return
    end if

    do k = 1, n_points
      norm_sum = 0.0_mytype
      do c = 1, 3
        f_fluid_to_fibre(k, c) = alpha * (u_f(k, c) - v_f(k, c))
        f_fibre_to_fluid(k, c) = -f_fluid_to_fibre(k, c)
        norm_sum = norm_sum + f_fluid_to_fibre(k, c) * f_fluid_to_fibre(k, c)
      end do
      force_norm(k) = sqrt(norm_sum)
    end do

    call update_finite_statuses(f_fluid_to_fibre, f_fibre_to_fluid, force_norm)
    call update_sign_convention_audit_status()
  end subroutine stage12_sign_convention_compute_pair

  subroutine stage12_sign_convention_audit_record_statuses(fluid_to_fibre_sign, fibre_to_fluid_sign, &
                                                           action_reaction, zero_slip_neutrality, &
                                                           positive_slip, negative_slip, &
                                                           structure_side_convention, fluid_side_convention)
    integer, intent(in) :: fluid_to_fibre_sign
    integer, intent(in) :: fibre_to_fluid_sign
    integer, intent(in) :: action_reaction
    integer, intent(in) :: zero_slip_neutrality
    integer, intent(in) :: positive_slip
    integer, intent(in) :: negative_slip
    integer, intent(in) :: structure_side_convention
    integer, intent(in) :: fluid_side_convention

    fluid_to_fibre_sign_status = fluid_to_fibre_sign
    fibre_to_fluid_sign_status = fibre_to_fluid_sign
    action_reaction_status = action_reaction
    zero_slip_neutrality_status = zero_slip_neutrality
    positive_slip_status = positive_slip
    negative_slip_status = negative_slip
    structure_side_convention_status = structure_side_convention
    fluid_side_convention_status = fluid_side_convention
    call update_sign_convention_audit_status()
  end subroutine stage12_sign_convention_audit_record_statuses

  subroutine stage12_sign_convention_audit_clear()
    fluid_to_fibre_sign_status = 0
    fibre_to_fluid_sign_status = 0
    action_reaction_status = 0
    zero_slip_neutrality_status = 0
    positive_slip_status = 0
    negative_slip_status = 0
    structure_side_convention_status = 0
    fluid_side_convention_status = 0
    finite_force_status = 0
    force_norm_finite_status = 0
    call update_sign_convention_audit_status()
  end subroutine stage12_sign_convention_audit_clear

  subroutine stage12_sign_convention_audit_finalize()
    initialized_status = 0
    call stage12_sign_convention_audit_clear()
  end subroutine stage12_sign_convention_audit_finalize

  subroutine stage12_sign_convention_audit_get_status_values(out_initialized_status, &
                                                             out_fluid_to_fibre_sign_status, &
                                                             out_fibre_to_fluid_sign_status, &
                                                             out_action_reaction_status, &
                                                             out_zero_slip_neutrality_status, &
                                                             out_positive_slip_status, &
                                                             out_negative_slip_status, &
                                                             out_structure_side_convention_status, &
                                                             out_fluid_side_convention_status, &
                                                             out_finite_force_status, &
                                                             out_force_norm_finite_status, &
                                                             out_no_eulerian_force_density_status, &
                                                             out_no_rhs_injection_status, &
                                                             out_no_ibm_spreading_status, &
                                                             out_no_feedback_application_status, &
                                                             out_no_twoway_force_status, &
                                                             out_no_structure_advance_status, &
                                                             out_no_fluid_field_access_status, &
                                                             out_no_fluid_field_modification_status, &
                                                             out_sign_convention_audit_status)
    integer, intent(out) :: out_initialized_status
    integer, intent(out) :: out_fluid_to_fibre_sign_status
    integer, intent(out) :: out_fibre_to_fluid_sign_status
    integer, intent(out) :: out_action_reaction_status
    integer, intent(out) :: out_zero_slip_neutrality_status
    integer, intent(out) :: out_positive_slip_status
    integer, intent(out) :: out_negative_slip_status
    integer, intent(out) :: out_structure_side_convention_status
    integer, intent(out) :: out_fluid_side_convention_status
    integer, intent(out) :: out_finite_force_status
    integer, intent(out) :: out_force_norm_finite_status
    integer, intent(out) :: out_no_eulerian_force_density_status
    integer, intent(out) :: out_no_rhs_injection_status
    integer, intent(out) :: out_no_ibm_spreading_status
    integer, intent(out) :: out_no_feedback_application_status
    integer, intent(out) :: out_no_twoway_force_status
    integer, intent(out) :: out_no_structure_advance_status
    integer, intent(out) :: out_no_fluid_field_access_status
    integer, intent(out) :: out_no_fluid_field_modification_status
    integer, intent(out) :: out_sign_convention_audit_status

    call update_sign_convention_audit_status()

    out_initialized_status = initialized_status
    out_fluid_to_fibre_sign_status = fluid_to_fibre_sign_status
    out_fibre_to_fluid_sign_status = fibre_to_fluid_sign_status
    out_action_reaction_status = action_reaction_status
    out_zero_slip_neutrality_status = zero_slip_neutrality_status
    out_positive_slip_status = positive_slip_status
    out_negative_slip_status = negative_slip_status
    out_structure_side_convention_status = structure_side_convention_status
    out_fluid_side_convention_status = fluid_side_convention_status
    out_finite_force_status = finite_force_status
    out_force_norm_finite_status = force_norm_finite_status
    out_no_eulerian_force_density_status = no_eulerian_force_density_status
    out_no_rhs_injection_status = no_rhs_injection_status
    out_no_ibm_spreading_status = no_ibm_spreading_status
    out_no_feedback_application_status = no_feedback_application_status
    out_no_twoway_force_status = no_twoway_force_status
    out_no_structure_advance_status = no_structure_advance_status
    out_no_fluid_field_access_status = no_fluid_field_access_status
    out_no_fluid_field_modification_status = no_fluid_field_modification_status
    out_sign_convention_audit_status = sign_convention_audit_status
  end subroutine stage12_sign_convention_audit_get_status_values

  subroutine stage12_sign_convention_audit_write_diagnostics(filename)
    character(len=*), intent(in), optional :: filename
    character(len=256) :: output_file
    integer :: io_unit
    integer :: ios

    if (present(filename)) then
      output_file = filename
    else
      output_file = 'stage12_outputs/fibre_stage12_sign_convention_audit.dat'
    end if

    call update_sign_convention_audit_status()
    open(newunit=io_unit, file=trim(output_file), status='replace', action='write', iostat=ios)
    if (ios /= 0) return
    write(io_unit, '(A,1X,I0)') 'stage12_sign_convention_initialized_status', initialized_status
    write(io_unit, '(A,1X,I0)') 'stage12_sign_convention_fluid_to_fibre_sign_status', fluid_to_fibre_sign_status
    write(io_unit, '(A,1X,I0)') 'stage12_sign_convention_fibre_to_fluid_sign_status', fibre_to_fluid_sign_status
    write(io_unit, '(A,1X,I0)') 'stage12_sign_convention_action_reaction_status', action_reaction_status
    write(io_unit, '(A,1X,I0)') 'stage12_sign_convention_status', sign_convention_audit_status
    close(io_unit)
  end subroutine stage12_sign_convention_audit_write_diagnostics

  subroutine update_finite_statuses(f_fluid_to_fibre, f_fibre_to_fluid, force_norm)
    real(mytype), intent(in) :: f_fluid_to_fibre(:,:)
    real(mytype), intent(in) :: f_fibre_to_fluid(:,:)
    real(mytype), intent(in) :: force_norm(:)
    integer :: i
    integer :: j

    finite_force_status = 1
    do i = 1, size(f_fluid_to_fibre, 1)
      do j = 1, size(f_fluid_to_fibre, 2)
        if (.not. is_finite(f_fluid_to_fibre(i, j))) finite_force_status = 0
        if (.not. is_finite(f_fibre_to_fluid(i, j))) finite_force_status = 0
      end do
    end do

    force_norm_finite_status = 1
    do i = 1, size(force_norm)
      if (.not. is_finite(force_norm(i))) force_norm_finite_status = 0
    end do
  end subroutine update_finite_statuses

  subroutine update_sign_convention_audit_status()
    if (initialized_status == 1 .and. fluid_to_fibre_sign_status == 1 .and. &
        fibre_to_fluid_sign_status == 1 .and. action_reaction_status == 1 .and. &
        zero_slip_neutrality_status == 1 .and. positive_slip_status == 1 .and. &
        negative_slip_status == 1 .and. structure_side_convention_status == 1 .and. &
        fluid_side_convention_status == 1 .and. finite_force_status == 1 .and. &
        force_norm_finite_status == 1 .and. no_eulerian_force_density_status == 1 .and. &
        no_rhs_injection_status == 1 .and. no_ibm_spreading_status == 1 .and. &
        no_feedback_application_status == 1 .and. no_twoway_force_status == 1 .and. &
        no_structure_advance_status == 1 .and. no_fluid_field_access_status == 1 .and. &
        no_fluid_field_modification_status == 1) then
      sign_convention_audit_status = 1
    else
      sign_convention_audit_status = 0
    end if
  end subroutine update_sign_convention_audit_status

  logical function is_finite(value)
    real(mytype), intent(in) :: value

    is_finite = (value == value) .and. (abs(value) < huge(value))
  end function is_finite

end module fibre_stage12_sign_convention_audit
