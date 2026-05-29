module fibre_stage12_feedback_formula
  implicit none
  private

  integer, parameter :: mytype = kind(1.0d0)

  integer :: initialized_status = 0
  integer :: zero_slip_status = 0
  integer :: positive_slip_sign_status = 0
  integer :: negative_slip_sign_status = 0
  integer :: force_sign_minus_status = 0
  integer :: multicomponent_slip_status = 0
  integer :: gain_scaling_status = 0
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
  integer :: feedback_formula_status = 0

  public :: stage12_feedback_formula_init
  public :: stage12_feedback_formula_compute_controlled
  public :: stage12_feedback_formula_clear
  public :: stage12_feedback_formula_finalize
  public :: stage12_feedback_formula_get_status_values
  public :: stage12_feedback_formula_write_diagnostics
  public :: stage12_feedback_formula_record_test_statuses

contains

  subroutine stage12_feedback_formula_init()
    initialized_status = 1
    zero_slip_status = 0
    positive_slip_sign_status = 0
    negative_slip_sign_status = 0
    force_sign_minus_status = 0
    multicomponent_slip_status = 0
    gain_scaling_status = 0
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
    call update_feedback_formula_status()
  end subroutine stage12_feedback_formula_init

  subroutine stage12_feedback_formula_compute_controlled(u_f, v_f, alpha, force_sign, f_fs_cand, force_norm)
    real(mytype), intent(in) :: u_f(:,:)
    real(mytype), intent(in) :: v_f(:,:)
    real(mytype), intent(in) :: alpha
    integer, intent(in) :: force_sign
    real(mytype), intent(out) :: f_fs_cand(:,:)
    real(mytype), intent(out) :: force_norm(:)
    integer :: n_points
    integer :: k
    integer :: c
    real(mytype) :: sign_value
    real(mytype) :: norm_sum

    finite_force_status = 0
    force_norm_finite_status = 0
    if (force_sign /= 1 .and. force_sign /= -1) then
      f_fs_cand(:, :) = 0.0_mytype
      force_norm(:) = 0.0_mytype
      call update_feedback_formula_status()
      return
    end if

    n_points = size(u_f, 1)
    if (size(u_f, 2) /= 3 .or. size(v_f, 1) /= n_points .or. size(v_f, 2) /= 3 .or. &
        size(f_fs_cand, 1) /= n_points .or. size(f_fs_cand, 2) /= 3 .or. size(force_norm) /= n_points) then
      f_fs_cand(:, :) = 0.0_mytype
      force_norm(:) = 0.0_mytype
      call update_feedback_formula_status()
      return
    end if

    sign_value = real(force_sign, mytype)
    do k = 1, n_points
      norm_sum = 0.0_mytype
      do c = 1, 3
        f_fs_cand(k, c) = sign_value * alpha * (u_f(k, c) - v_f(k, c))
        norm_sum = norm_sum + f_fs_cand(k, c) * f_fs_cand(k, c)
      end do
      force_norm(k) = sqrt(norm_sum)
    end do

    call update_finite_statuses(f_fs_cand, force_norm)
    call update_feedback_formula_status()
  end subroutine stage12_feedback_formula_compute_controlled

  subroutine stage12_feedback_formula_record_test_statuses(zero_slip, positive_slip, negative_slip, &
                                                           force_sign_minus, multicomponent_slip, gain_scaling)
    integer, intent(in) :: zero_slip
    integer, intent(in) :: positive_slip
    integer, intent(in) :: negative_slip
    integer, intent(in) :: force_sign_minus
    integer, intent(in) :: multicomponent_slip
    integer, intent(in) :: gain_scaling

    zero_slip_status = zero_slip
    positive_slip_sign_status = positive_slip
    negative_slip_sign_status = negative_slip
    force_sign_minus_status = force_sign_minus
    multicomponent_slip_status = multicomponent_slip
    gain_scaling_status = gain_scaling
    call update_feedback_formula_status()
  end subroutine stage12_feedback_formula_record_test_statuses

  subroutine stage12_feedback_formula_clear()
    zero_slip_status = 0
    positive_slip_sign_status = 0
    negative_slip_sign_status = 0
    force_sign_minus_status = 0
    multicomponent_slip_status = 0
    gain_scaling_status = 0
    finite_force_status = 0
    force_norm_finite_status = 0
    call update_feedback_formula_status()
  end subroutine stage12_feedback_formula_clear

  subroutine stage12_feedback_formula_finalize()
    initialized_status = 0
    call stage12_feedback_formula_clear()
  end subroutine stage12_feedback_formula_finalize

  subroutine stage12_feedback_formula_get_status_values(out_initialized_status, out_zero_slip_status, &
                                                        out_positive_slip_sign_status, &
                                                        out_negative_slip_sign_status, &
                                                        out_force_sign_minus_status, &
                                                        out_multicomponent_slip_status, out_gain_scaling_status, &
                                                        out_finite_force_status, out_force_norm_finite_status, &
                                                        out_no_eulerian_force_density_status, &
                                                        out_no_rhs_injection_status, out_no_ibm_spreading_status, &
                                                        out_no_feedback_application_status, out_no_twoway_force_status, &
                                                        out_no_structure_advance_status, out_no_fluid_field_access_status, &
                                                        out_no_fluid_field_modification_status, out_feedback_formula_status)
    integer, intent(out) :: out_initialized_status
    integer, intent(out) :: out_zero_slip_status
    integer, intent(out) :: out_positive_slip_sign_status
    integer, intent(out) :: out_negative_slip_sign_status
    integer, intent(out) :: out_force_sign_minus_status
    integer, intent(out) :: out_multicomponent_slip_status
    integer, intent(out) :: out_gain_scaling_status
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
    integer, intent(out) :: out_feedback_formula_status

    call update_feedback_formula_status()

    out_initialized_status = initialized_status
    out_zero_slip_status = zero_slip_status
    out_positive_slip_sign_status = positive_slip_sign_status
    out_negative_slip_sign_status = negative_slip_sign_status
    out_force_sign_minus_status = force_sign_minus_status
    out_multicomponent_slip_status = multicomponent_slip_status
    out_gain_scaling_status = gain_scaling_status
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
    out_feedback_formula_status = feedback_formula_status
  end subroutine stage12_feedback_formula_get_status_values

  subroutine stage12_feedback_formula_write_diagnostics(filename)
    character(len=*), intent(in), optional :: filename
    character(len=256) :: output_file
    integer :: io_unit
    integer :: ios

    if (present(filename)) then
      output_file = filename
    else
      output_file = 'stage12_outputs/fibre_stage12_feedback_formula.dat'
    end if

    call update_feedback_formula_status()
    open(newunit=io_unit, file=trim(output_file), status='replace', action='write', iostat=ios)
    if (ios /= 0) return
    write(io_unit, '(A,1X,I0)') 'stage12_feedback_formula_initialized_status', initialized_status
    write(io_unit, '(A,1X,I0)') 'stage12_feedback_formula_zero_slip_status', zero_slip_status
    write(io_unit, '(A,1X,I0)') 'stage12_feedback_formula_positive_slip_sign_status', positive_slip_sign_status
    write(io_unit, '(A,1X,I0)') 'stage12_feedback_formula_negative_slip_sign_status', negative_slip_sign_status
    write(io_unit, '(A,1X,I0)') 'stage12_feedback_formula_force_sign_minus_status', force_sign_minus_status
    write(io_unit, '(A,1X,I0)') 'stage12_feedback_formula_multicomponent_slip_status', multicomponent_slip_status
    write(io_unit, '(A,1X,I0)') 'stage12_feedback_formula_gain_scaling_status', gain_scaling_status
    write(io_unit, '(A,1X,I0)') 'stage12_feedback_formula_status', feedback_formula_status
    close(io_unit)
  end subroutine stage12_feedback_formula_write_diagnostics

  subroutine update_finite_statuses(f_fs_cand, force_norm)
    real(mytype), intent(in) :: f_fs_cand(:,:)
    real(mytype), intent(in) :: force_norm(:)
    integer :: i
    integer :: j

    finite_force_status = 1
    do i = 1, size(f_fs_cand, 1)
      do j = 1, size(f_fs_cand, 2)
        if (.not. is_finite(f_fs_cand(i, j))) finite_force_status = 0
      end do
    end do

    force_norm_finite_status = 1
    do i = 1, size(force_norm)
      if (.not. is_finite(force_norm(i))) force_norm_finite_status = 0
    end do
  end subroutine update_finite_statuses

  subroutine update_feedback_formula_status()
    if (initialized_status == 1 .and. zero_slip_status == 1 .and. positive_slip_sign_status == 1 .and. &
        negative_slip_sign_status == 1 .and. force_sign_minus_status == 1 .and. &
        multicomponent_slip_status == 1 .and. gain_scaling_status == 1 .and. finite_force_status == 1 .and. &
        force_norm_finite_status == 1 .and. no_eulerian_force_density_status == 1 .and. &
        no_rhs_injection_status == 1 .and. no_ibm_spreading_status == 1 .and. &
        no_feedback_application_status == 1 .and. no_twoway_force_status == 1 .and. &
        no_structure_advance_status == 1 .and. no_fluid_field_access_status == 1 .and. &
        no_fluid_field_modification_status == 1) then
      feedback_formula_status = 1
    else
      feedback_formula_status = 0
    end if
  end subroutine update_feedback_formula_status

  logical function is_finite(value)
    real(mytype), intent(in) :: value

    is_finite = (value == value) .and. (abs(value) < huge(value))
  end function is_finite

end module fibre_stage12_feedback_formula
