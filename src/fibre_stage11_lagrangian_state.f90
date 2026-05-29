module fibre_stage11_lagrangian_state
  use decomp_2d_constants, only : mytype
  implicit none
  private

  integer :: n_points = 0
  real(mytype), allocatable :: x_f(:,:)
  real(mytype), allocatable :: u_f(:,:)
  integer, allocatable :: valid_point_flag(:)

  integer :: allocated_status = 0
  integer :: point_count_status = 1
  integer :: coordinates_finite_status = 1
  integer :: velocity_placeholder_status = 1
  integer :: valid_flag_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_velocity_sampling_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_force_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: state_status = 1

  public :: stage11_lagrangian_state_init
  public :: stage11_lagrangian_state_finalize
  public :: stage11_lagrangian_state_is_allocated
  public :: stage11_lagrangian_state_get_count
  public :: stage11_lagrangian_state_get_status_values
  public :: stage11_lagrangian_state_write_diagnostics

contains

  subroutine stage11_lagrangian_state_init(input_points)
    integer, intent(in) :: input_points
    integer :: k
    real(mytype) :: denom

    call stage11_lagrangian_state_finalize()

    n_points = max(0, input_points)
    point_count_status = 1
    if (n_points <= 0) then
      point_count_status = 0
      state_status = 0
      return
    end if

    allocate(x_f(n_points,3))
    allocate(u_f(n_points,3))
    allocate(valid_point_flag(n_points))

    allocated_status = 1
    if (.not. stage11_lagrangian_state_is_allocated()) then
      allocated_status = 0
      state_status = 0
      return
    end if

    denom = real(max(1, n_points - 1), mytype)
    do k = 1, n_points
      x_f(k,1) = real(k - 1, mytype) / denom
      x_f(k,2) = 0.5_mytype
      x_f(k,3) = 0.5_mytype
    end do

    u_f(:,:) = 0.0_mytype
    valid_point_flag(:) = 1

    call update_statuses()
  end subroutine stage11_lagrangian_state_init

  subroutine stage11_lagrangian_state_finalize()
    if (allocated(x_f)) deallocate(x_f)
    if (allocated(u_f)) deallocate(u_f)
    if (allocated(valid_point_flag)) deallocate(valid_point_flag)

    n_points = 0
    allocated_status = 0
    point_count_status = 1
    coordinates_finite_status = 1
    velocity_placeholder_status = 1
    valid_flag_status = 1
    no_fluid_field_access_status = 1
    no_velocity_sampling_status = 1
    no_fluid_field_modification_status = 1
    no_rhs_injection_status = 1
    no_ibm_spreading_status = 1
    no_feedback_force_status = 1
    no_twoway_force_status = 1
    no_structure_advance_status = 1
    state_status = 1
  end subroutine stage11_lagrangian_state_finalize

  logical function stage11_lagrangian_state_is_allocated()
    stage11_lagrangian_state_is_allocated = allocated(x_f) .and. allocated(u_f) .and. allocated(valid_point_flag)
  end function stage11_lagrangian_state_is_allocated

  integer function stage11_lagrangian_state_get_count()
    stage11_lagrangian_state_get_count = n_points
  end function stage11_lagrangian_state_get_count

  subroutine stage11_lagrangian_state_write_diagnostics(unit)
    integer, intent(in) :: unit
    write(unit, '(A,1X,I0)') 'stage11_1_allocated_status', allocated_status
    write(unit, '(A,1X,I0)') 'stage11_1_point_count_status', point_count_status
    write(unit, '(A,1X,I0)') 'stage11_1_coordinates_finite_status', coordinates_finite_status
    write(unit, '(A,1X,I0)') 'stage11_1_velocity_placeholder_status', velocity_placeholder_status
    write(unit, '(A,1X,I0)') 'stage11_1_valid_flag_status', valid_flag_status
    write(unit, '(A,1X,I0)') 'stage11_1_no_fluid_field_access_status', no_fluid_field_access_status
    write(unit, '(A,1X,I0)') 'stage11_1_no_velocity_sampling_status', no_velocity_sampling_status
    write(unit, '(A,1X,I0)') 'stage11_1_no_fluid_field_modification_status', no_fluid_field_modification_status
    write(unit, '(A,1X,I0)') 'stage11_1_no_rhs_injection_status', no_rhs_injection_status
    write(unit, '(A,1X,I0)') 'stage11_1_no_ibm_spreading_status', no_ibm_spreading_status
    write(unit, '(A,1X,I0)') 'stage11_1_no_feedback_force_status', no_feedback_force_status
    write(unit, '(A,1X,I0)') 'stage11_1_no_twoway_force_status', no_twoway_force_status
    write(unit, '(A,1X,I0)') 'stage11_1_no_structure_advance_status', no_structure_advance_status
    write(unit, '(A,1X,I0)') 'stage11_1_lagrangian_state_status', state_status
  end subroutine stage11_lagrangian_state_write_diagnostics

  subroutine stage11_lagrangian_state_get_status_values(out_allocated_status, out_point_count_status, &
                                                        out_coordinates_finite_status, out_velocity_placeholder_status, &
                                                        out_valid_flag_status, out_no_fluid_field_access_status, &
                                                        out_no_velocity_sampling_status, out_no_fluid_field_modification_status, &
                                                        out_no_rhs_injection_status, out_no_ibm_spreading_status, &
                                                        out_no_feedback_force_status, out_no_twoway_force_status, &
                                                        out_no_structure_advance_status, out_state_status)
    integer, intent(out) :: out_allocated_status
    integer, intent(out) :: out_point_count_status
    integer, intent(out) :: out_coordinates_finite_status
    integer, intent(out) :: out_velocity_placeholder_status
    integer, intent(out) :: out_valid_flag_status
    integer, intent(out) :: out_no_fluid_field_access_status
    integer, intent(out) :: out_no_velocity_sampling_status
    integer, intent(out) :: out_no_fluid_field_modification_status
    integer, intent(out) :: out_no_rhs_injection_status
    integer, intent(out) :: out_no_ibm_spreading_status
    integer, intent(out) :: out_no_feedback_force_status
    integer, intent(out) :: out_no_twoway_force_status
    integer, intent(out) :: out_no_structure_advance_status
    integer, intent(out) :: out_state_status

    out_allocated_status = allocated_status
    out_point_count_status = point_count_status
    out_coordinates_finite_status = coordinates_finite_status
    out_velocity_placeholder_status = velocity_placeholder_status
    out_valid_flag_status = valid_flag_status
    out_no_fluid_field_access_status = no_fluid_field_access_status
    out_no_velocity_sampling_status = no_velocity_sampling_status
    out_no_fluid_field_modification_status = no_fluid_field_modification_status
    out_no_rhs_injection_status = no_rhs_injection_status
    out_no_ibm_spreading_status = no_ibm_spreading_status
    out_no_feedback_force_status = no_feedback_force_status
    out_no_twoway_force_status = no_twoway_force_status
    out_no_structure_advance_status = no_structure_advance_status
    out_state_status = state_status
  end subroutine stage11_lagrangian_state_get_status_values

  subroutine update_statuses()
    integer :: k, j

    coordinates_finite_status = 1
    velocity_placeholder_status = 1
    valid_flag_status = 1

    do k = 1, n_points
      do j = 1, 3
        if (.not. finite_value(x_f(k,j))) coordinates_finite_status = 0
      end do
    end do

    do k = 1, n_points
      do j = 1, 3
        if (u_f(k,j) /= 0.0_mytype) velocity_placeholder_status = 0
      end do
      if (valid_point_flag(k) /= 1) valid_flag_status = 0
    end do

    state_status = 1
    if (allocated_status /= 1) state_status = 0
    if (point_count_status /= 1) state_status = 0
    if (coordinates_finite_status /= 1) state_status = 0
    if (velocity_placeholder_status /= 1) state_status = 0
    if (valid_flag_status /= 1) state_status = 0
    if (no_fluid_field_access_status /= 1) state_status = 0
    if (no_velocity_sampling_status /= 1) state_status = 0
    if (no_fluid_field_modification_status /= 1) state_status = 0
    if (no_rhs_injection_status /= 1) state_status = 0
    if (no_ibm_spreading_status /= 1) state_status = 0
    if (no_feedback_force_status /= 1) state_status = 0
    if (no_twoway_force_status /= 1) state_status = 0
    if (no_structure_advance_status /= 1) state_status = 0
  end subroutine update_statuses

  logical function finite_value(value)
    real(mytype), intent(in) :: value
    finite_value = (value == value) .and. (abs(value) < huge(value))
  end function finite_value

end module fibre_stage11_lagrangian_state
