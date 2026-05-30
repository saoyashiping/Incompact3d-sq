module fibre_stage12_power_diagnostics
  implicit none
  private

  integer, parameter :: mytype = kind(1.0d0)

  integer :: initialized_status = 0
  integer :: zero_slip_power_status = 0
  integer :: positive_slip_power_status = 0
  integer :: force_norm_finite_status = 0
  integer :: slip_power_finite_status = 0
  integer :: structure_power_finite_status = 0
  integer :: fluid_power_finite_status = 0
  integer :: pair_power_consistency_status = 0
  integer :: gain_scaling_power_status = 0
  integer :: action_reaction_power_status = 0
  integer :: finite_diagnostics_status = 0
  integer :: no_eulerian_force_density_status = 1
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: power_diagnostics_status = 0

  public :: stage12_power_diagnostics_init
  public :: stage12_power_diagnostics_compute
  public :: stage12_power_diagnostics_clear
  public :: stage12_power_diagnostics_finalize
  public :: stage12_power_diagnostics_get_status_values
  public :: stage12_power_diagnostics_write_diagnostics
  public :: stage12_power_diagnostics_record_statuses

contains

  subroutine stage12_power_diagnostics_init()
    initialized_status = 1
    zero_slip_power_status = 0
    positive_slip_power_status = 0
    force_norm_finite_status = 0
    slip_power_finite_status = 0
    structure_power_finite_status = 0
    fluid_power_finite_status = 0
    pair_power_consistency_status = 0
    gain_scaling_power_status = 0
    action_reaction_power_status = 0
    finite_diagnostics_status = 0
    no_eulerian_force_density_status = 1
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_application_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    no_fluid_field_access_status = 1
    no_fluid_field_modification_status = 1
    call update_power_diagnostics_status()
  end subroutine stage12_power_diagnostics_init

  subroutine stage12_power_diagnostics_compute(u_f, v_f, f_fluid_to_fibre, f_fibre_to_fluid, force_norm, &
                                               total_force_l2, p_slip, p_structure, p_fluid, p_pair)
    real(mytype), intent(in) :: u_f(:,:)
    real(mytype), intent(in) :: v_f(:,:)
    real(mytype), intent(in) :: f_fluid_to_fibre(:,:)
    real(mytype), intent(in) :: f_fibre_to_fluid(:,:)
    real(mytype), intent(out) :: force_norm(:)
    real(mytype), intent(out) :: total_force_l2
    real(mytype), intent(out) :: p_slip
    real(mytype), intent(out) :: p_structure
    real(mytype), intent(out) :: p_fluid
    real(mytype), intent(out) :: p_pair
    integer :: n_points
    integer :: k
    integer :: c
    real(mytype) :: force_sq_sum
    real(mytype) :: point_force_sq
    real(mytype) :: slip_component

    call reset_finite_statuses()
    total_force_l2 = 0.0_mytype
    p_slip = 0.0_mytype
    p_structure = 0.0_mytype
    p_fluid = 0.0_mytype
    p_pair = 0.0_mytype
    force_norm(:) = 0.0_mytype

    n_points = size(u_f, 1)
    if (size(u_f, 2) /= 3 .or. size(v_f, 1) /= n_points .or. size(v_f, 2) /= 3 .or. &
        size(f_fluid_to_fibre, 1) /= n_points .or. size(f_fluid_to_fibre, 2) /= 3 .or. &
        size(f_fibre_to_fluid, 1) /= n_points .or. size(f_fibre_to_fluid, 2) /= 3 .or. &
        size(force_norm) /= n_points) then
      call update_power_diagnostics_status()
      return
    end if

    force_sq_sum = 0.0_mytype
    do k = 1, n_points
      point_force_sq = 0.0_mytype
      do c = 1, 3
        slip_component = u_f(k, c) - v_f(k, c)
        point_force_sq = point_force_sq + f_fluid_to_fibre(k, c) * f_fluid_to_fibre(k, c)
        p_slip = p_slip + f_fluid_to_fibre(k, c) * slip_component
        p_structure = p_structure + f_fluid_to_fibre(k, c) * v_f(k, c)
        p_fluid = p_fluid + f_fibre_to_fluid(k, c) * u_f(k, c)
      end do
      force_norm(k) = sqrt(point_force_sq)
      force_sq_sum = force_sq_sum + point_force_sq
    end do
    total_force_l2 = sqrt(force_sq_sum)
    p_pair = p_structure + p_fluid

    call update_finite_statuses(force_norm, total_force_l2, p_slip, p_structure, p_fluid, p_pair)
    call update_power_diagnostics_status()
  end subroutine stage12_power_diagnostics_compute

  subroutine stage12_power_diagnostics_record_statuses(zero_slip_power, positive_slip_power, &
                                                       pair_power_consistency, gain_scaling_power, &
                                                       action_reaction_power, finite_diagnostics)
    integer, intent(in) :: zero_slip_power
    integer, intent(in) :: positive_slip_power
    integer, intent(in) :: pair_power_consistency
    integer, intent(in) :: gain_scaling_power
    integer, intent(in) :: action_reaction_power
    integer, intent(in) :: finite_diagnostics

    zero_slip_power_status = zero_slip_power
    positive_slip_power_status = positive_slip_power
    pair_power_consistency_status = pair_power_consistency
    gain_scaling_power_status = gain_scaling_power
    action_reaction_power_status = action_reaction_power
    finite_diagnostics_status = finite_diagnostics
    call update_power_diagnostics_status()
  end subroutine stage12_power_diagnostics_record_statuses

  subroutine stage12_power_diagnostics_clear()
    zero_slip_power_status = 0
    positive_slip_power_status = 0
    force_norm_finite_status = 0
    slip_power_finite_status = 0
    structure_power_finite_status = 0
    fluid_power_finite_status = 0
    pair_power_consistency_status = 0
    gain_scaling_power_status = 0
    action_reaction_power_status = 0
    finite_diagnostics_status = 0
    call update_power_diagnostics_status()
  end subroutine stage12_power_diagnostics_clear

  subroutine stage12_power_diagnostics_finalize()
    initialized_status = 0
    call stage12_power_diagnostics_clear()
  end subroutine stage12_power_diagnostics_finalize

  subroutine stage12_power_diagnostics_get_status_values(out_initialized_status, out_zero_slip_power_status, &
                                                         out_positive_slip_power_status, &
                                                         out_force_norm_finite_status, &
                                                         out_slip_power_finite_status, &
                                                         out_structure_power_finite_status, &
                                                         out_fluid_power_finite_status, &
                                                         out_pair_power_consistency_status, &
                                                         out_gain_scaling_power_status, &
                                                         out_action_reaction_power_status, &
                                                         out_finite_diagnostics_status, &
                                                         out_no_eulerian_force_density_status, &
                                                         out_no_rhs_injection_status, &
                                                         out_no_ibm_spreading_status, &
                                                         out_no_feedback_application_status, &
                                                         out_no_twoway_force_status, &
                                                         out_no_structure_advance_status, &
                                                         out_no_fluid_field_access_status, &
                                                         out_no_fluid_field_modification_status, &
                                                         out_power_diagnostics_status)
    integer, intent(out) :: out_initialized_status
    integer, intent(out) :: out_zero_slip_power_status
    integer, intent(out) :: out_positive_slip_power_status
    integer, intent(out) :: out_force_norm_finite_status
    integer, intent(out) :: out_slip_power_finite_status
    integer, intent(out) :: out_structure_power_finite_status
    integer, intent(out) :: out_fluid_power_finite_status
    integer, intent(out) :: out_pair_power_consistency_status
    integer, intent(out) :: out_gain_scaling_power_status
    integer, intent(out) :: out_action_reaction_power_status
    integer, intent(out) :: out_finite_diagnostics_status
    integer, intent(out) :: out_no_eulerian_force_density_status
    integer, intent(out) :: out_no_rhs_injection_status
    integer, intent(out) :: out_no_ibm_spreading_status
    integer, intent(out) :: out_no_feedback_application_status
    integer, intent(out) :: out_no_twoway_force_status
    integer, intent(out) :: out_no_structure_advance_status
    integer, intent(out) :: out_no_fluid_field_access_status
    integer, intent(out) :: out_no_fluid_field_modification_status
    integer, intent(out) :: out_power_diagnostics_status

    call update_power_diagnostics_status()

    out_initialized_status = initialized_status
    out_zero_slip_power_status = zero_slip_power_status
    out_positive_slip_power_status = positive_slip_power_status
    out_force_norm_finite_status = force_norm_finite_status
    out_slip_power_finite_status = slip_power_finite_status
    out_structure_power_finite_status = structure_power_finite_status
    out_fluid_power_finite_status = fluid_power_finite_status
    out_pair_power_consistency_status = pair_power_consistency_status
    out_gain_scaling_power_status = gain_scaling_power_status
    out_action_reaction_power_status = action_reaction_power_status
    out_finite_diagnostics_status = finite_diagnostics_status
    out_no_eulerian_force_density_status = no_eulerian_force_density_status
    out_no_rhs_injection_status = no_rhs_injection_status
    out_no_ibm_spreading_status = no_ibm_spreading_status
    out_no_feedback_application_status = no_feedback_application_status
    out_no_twoway_force_status = no_twoway_force_status
    out_no_structure_advance_status = no_structure_advance_status
    out_no_fluid_field_access_status = no_fluid_field_access_status
    out_no_fluid_field_modification_status = no_fluid_field_modification_status
    out_power_diagnostics_status = power_diagnostics_status
  end subroutine stage12_power_diagnostics_get_status_values

  subroutine stage12_power_diagnostics_write_diagnostics(filename)
    character(len=*), intent(in), optional :: filename
    character(len=256) :: output_file
    integer :: io_unit
    integer :: ios

    if (present(filename)) then
      output_file = filename
    else
      output_file = 'stage12_outputs/fibre_stage12_power_diagnostics.dat'
    end if

    call update_power_diagnostics_status()
    open(newunit=io_unit, file=trim(output_file), status='replace', action='write', iostat=ios)
    if (ios /= 0) return
    write(io_unit, '(A,1X,I0)') 'stage12_power_diagnostics_initialized_status', initialized_status
    write(io_unit, '(A,1X,I0)') 'stage12_power_diagnostics_zero_slip_power_status', zero_slip_power_status
    write(io_unit, '(A,1X,I0)') 'stage12_power_diagnostics_positive_slip_power_status', positive_slip_power_status
    write(io_unit, '(A,1X,I0)') 'stage12_power_diagnostics_status', power_diagnostics_status
    close(io_unit)
  end subroutine stage12_power_diagnostics_write_diagnostics

  subroutine reset_finite_statuses()
    force_norm_finite_status = 0
    slip_power_finite_status = 0
    structure_power_finite_status = 0
    fluid_power_finite_status = 0
    finite_diagnostics_status = 0
  end subroutine reset_finite_statuses

  subroutine update_finite_statuses(force_norm, total_force_l2, p_slip, p_structure, p_fluid, p_pair)
    real(mytype), intent(in) :: force_norm(:)
    real(mytype), intent(in) :: total_force_l2
    real(mytype), intent(in) :: p_slip
    real(mytype), intent(in) :: p_structure
    real(mytype), intent(in) :: p_fluid
    real(mytype), intent(in) :: p_pair
    integer :: i

    force_norm_finite_status = 1
    do i = 1, size(force_norm)
      if (.not. is_finite(force_norm(i))) force_norm_finite_status = 0
    end do
    if (.not. is_finite(total_force_l2)) force_norm_finite_status = 0

    slip_power_finite_status = merge(1, 0, is_finite(p_slip))
    structure_power_finite_status = merge(1, 0, is_finite(p_structure))
    fluid_power_finite_status = merge(1, 0, is_finite(p_fluid))

    if (force_norm_finite_status == 1 .and. slip_power_finite_status == 1 .and. &
        structure_power_finite_status == 1 .and. fluid_power_finite_status == 1 .and. is_finite(p_pair)) then
      finite_diagnostics_status = 1
    else
      finite_diagnostics_status = 0
    end if
  end subroutine update_finite_statuses

  subroutine update_power_diagnostics_status()
    if (initialized_status == 1 .and. zero_slip_power_status == 1 .and. positive_slip_power_status == 1 .and. &
        force_norm_finite_status == 1 .and. slip_power_finite_status == 1 .and. &
        structure_power_finite_status == 1 .and. fluid_power_finite_status == 1 .and. &
        pair_power_consistency_status == 1 .and. gain_scaling_power_status == 1 .and. &
        action_reaction_power_status == 1 .and. finite_diagnostics_status == 1 .and. &
        no_eulerian_force_density_status == 1 .and. no_rhs_injection_status == 1 .and. &
        no_ibm_spreading_status == 1 .and. no_feedback_application_status == 1 .and. &
        no_twoway_force_status == 1 .and. no_structure_advance_status == 1 .and. &
        no_fluid_field_access_status == 1 .and. no_fluid_field_modification_status == 1) then
      power_diagnostics_status = 1
    else
      power_diagnostics_status = 0
    end if
  end subroutine update_power_diagnostics_status

  logical function is_finite(value)
    real(mytype), intent(in) :: value

    is_finite = (value == value) .and. (abs(value) < huge(value))
  end function is_finite

end module fibre_stage12_power_diagnostics
