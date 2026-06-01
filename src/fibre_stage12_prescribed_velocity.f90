module fibre_stage12_prescribed_velocity
  implicit none
  private

  integer, parameter :: mytype = kind(1.0d0)

  integer :: n_points = 0
  real(mytype), allocatable :: v_f(:,:)
  real(mytype), allocatable :: velocity_norm(:)
  integer, allocatable :: velocity_valid_flag(:)

  integer :: allocated_status = 0
  integer :: point_count_status = 0
  integer :: zero_velocity_status = 0
  integer :: constant_velocity_status = 0
  integer :: velocity_norm_finite_status = 0
  integer :: velocity_valid_flag_status = 0
  integer :: clear_status = 0
  integer :: no_force_computation_status = 1
  integer :: no_eulerian_force_density_status = 1
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: prescribed_velocity_status = 0

  public :: stage12_prescribed_velocity_init
  public :: stage12_prescribed_velocity_set_zero
  public :: stage12_prescribed_velocity_set_constant
  public :: stage12_prescribed_velocity_clear
  public :: stage12_prescribed_velocity_finalize
  public :: stage12_prescribed_velocity_is_allocated
  public :: stage12_prescribed_velocity_get_count
  public :: stage12_prescribed_velocity_get_status_values
  public :: stage12_prescribed_velocity_write_diagnostics
  public :: stage12_prescribed_velocity_get_max_abs_velocity

contains

  subroutine stage12_prescribed_velocity_init(requested_points)
    integer, intent(in) :: requested_points
    integer :: alloc_status

    call stage12_prescribed_velocity_finalize()

    if (requested_points <= 0) then
      call update_status_values()
      return
    end if

    n_points = requested_points
    allocate(v_f(n_points, 3), stat=alloc_status)
    if (alloc_status /= 0) then
      n_points = 0
      call update_status_values()
      return
    end if

    allocate(velocity_norm(n_points), stat=alloc_status)
    if (alloc_status /= 0) then
      call stage12_prescribed_velocity_finalize()
      return
    end if

    allocate(velocity_valid_flag(n_points), stat=alloc_status)
    if (alloc_status /= 0) then
      call stage12_prescribed_velocity_finalize()
      return
    end if

    v_f(:, :) = 0.0_mytype
    velocity_norm(:) = 0.0_mytype
    velocity_valid_flag(:) = 1
    clear_status = 0
    constant_velocity_status = 0
    call update_status_values()
  end subroutine stage12_prescribed_velocity_init

  subroutine stage12_prescribed_velocity_set_zero()
    if (stage12_prescribed_velocity_is_allocated()) then
      v_f(:, :) = 0.0_mytype
      velocity_norm(:) = 0.0_mytype
      velocity_valid_flag(:) = 1
    end if
    call update_status_values()
  end subroutine stage12_prescribed_velocity_set_zero

  subroutine stage12_prescribed_velocity_set_constant(vx_const, vy_const, vz_const)
    real(mytype), intent(in) :: vx_const
    real(mytype), intent(in) :: vy_const
    real(mytype), intent(in) :: vz_const
    real(mytype) :: norm_value

    if (stage12_prescribed_velocity_is_allocated()) then
      norm_value = sqrt(vx_const * vx_const + vy_const * vy_const + vz_const * vz_const)
      v_f(:, 1) = vx_const
      v_f(:, 2) = vy_const
      v_f(:, 3) = vz_const
      velocity_norm(:) = norm_value
      velocity_valid_flag(:) = 1
      constant_velocity_status = 1
      if (.not. is_finite(vx_const) .or. .not. is_finite(vy_const) .or. .not. is_finite(vz_const)) then
        constant_velocity_status = 0
      end if
    else
      constant_velocity_status = 0
    end if
    call update_status_values()
  end subroutine stage12_prescribed_velocity_set_constant

  subroutine stage12_prescribed_velocity_clear()
    if (stage12_prescribed_velocity_is_allocated()) then
      v_f(:, :) = 0.0_mytype
      velocity_norm(:) = 0.0_mytype
      velocity_valid_flag(:) = 1
      clear_status = 1
    else
      clear_status = 0
    end if
    call update_status_values()
  end subroutine stage12_prescribed_velocity_clear

  subroutine stage12_prescribed_velocity_finalize()
    if (allocated(v_f)) deallocate(v_f)
    if (allocated(velocity_norm)) deallocate(velocity_norm)
    if (allocated(velocity_valid_flag)) deallocate(velocity_valid_flag)
    n_points = 0
    allocated_status = 0
    point_count_status = 0
    zero_velocity_status = 0
    constant_velocity_status = 0
    velocity_norm_finite_status = 0
    velocity_valid_flag_status = 0
    clear_status = 0
    prescribed_velocity_status = 0
  end subroutine stage12_prescribed_velocity_finalize

  logical function stage12_prescribed_velocity_is_allocated()
    stage12_prescribed_velocity_is_allocated = allocated(v_f) .and. allocated(velocity_norm) .and. &
                                               allocated(velocity_valid_flag)
  end function stage12_prescribed_velocity_is_allocated

  integer function stage12_prescribed_velocity_get_count()
    stage12_prescribed_velocity_get_count = n_points
  end function stage12_prescribed_velocity_get_count

  real(mytype) function stage12_prescribed_velocity_get_max_abs_velocity()
    if (allocated(v_f)) then
      stage12_prescribed_velocity_get_max_abs_velocity = maxval(abs(v_f))
    else
      stage12_prescribed_velocity_get_max_abs_velocity = 0.0_mytype
    end if
  end function stage12_prescribed_velocity_get_max_abs_velocity

  subroutine stage12_prescribed_velocity_get_status_values(out_allocated_status, out_point_count_status, &
                                                           out_zero_velocity_status, out_constant_velocity_status, &
                                                           out_velocity_norm_finite_status, &
                                                           out_velocity_valid_flag_status, out_clear_status, &
                                                           out_no_force_computation_status, &
                                                           out_no_eulerian_force_density_status, &
                                                           out_no_rhs_injection_status, &
                                                           out_no_ibm_spreading_status, &
                                                           out_no_feedback_application_status, &
                                                           out_no_twoway_force_status, &
                                                           out_no_structure_advance_status, &
                                                           out_no_fluid_field_access_status, &
                                                           out_no_fluid_field_modification_status, &
                                                           out_prescribed_velocity_status)
    integer, intent(out) :: out_allocated_status
    integer, intent(out) :: out_point_count_status
    integer, intent(out) :: out_zero_velocity_status
    integer, intent(out) :: out_constant_velocity_status
    integer, intent(out) :: out_velocity_norm_finite_status
    integer, intent(out) :: out_velocity_valid_flag_status
    integer, intent(out) :: out_clear_status
    integer, intent(out) :: out_no_force_computation_status
    integer, intent(out) :: out_no_eulerian_force_density_status
    integer, intent(out) :: out_no_rhs_injection_status
    integer, intent(out) :: out_no_ibm_spreading_status
    integer, intent(out) :: out_no_feedback_application_status
    integer, intent(out) :: out_no_twoway_force_status
    integer, intent(out) :: out_no_structure_advance_status
    integer, intent(out) :: out_no_fluid_field_access_status
    integer, intent(out) :: out_no_fluid_field_modification_status
    integer, intent(out) :: out_prescribed_velocity_status

    call update_status_values()

    out_allocated_status = allocated_status
    out_point_count_status = point_count_status
    out_zero_velocity_status = zero_velocity_status
    out_constant_velocity_status = constant_velocity_status
    out_velocity_norm_finite_status = velocity_norm_finite_status
    out_velocity_valid_flag_status = velocity_valid_flag_status
    out_clear_status = clear_status
    out_no_force_computation_status = no_force_computation_status
    out_no_eulerian_force_density_status = no_eulerian_force_density_status
    out_no_rhs_injection_status = no_rhs_injection_status
    out_no_ibm_spreading_status = no_ibm_spreading_status
    out_no_feedback_application_status = no_feedback_application_status
    out_no_twoway_force_status = no_twoway_force_status
    out_no_structure_advance_status = no_structure_advance_status
    out_no_fluid_field_access_status = no_fluid_field_access_status
    out_no_fluid_field_modification_status = no_fluid_field_modification_status
    out_prescribed_velocity_status = prescribed_velocity_status
  end subroutine stage12_prescribed_velocity_get_status_values

  subroutine stage12_prescribed_velocity_write_diagnostics(filename)
    character(len=*), intent(in), optional :: filename
    character(len=256) :: output_file
    integer :: io_unit
    integer :: ios

    if (present(filename)) then
      output_file = filename
    else
      output_file = 'stage12_outputs/fibre_stage12_prescribed_velocity.dat'
    end if

    call update_status_values()
    open(newunit=io_unit, file=trim(output_file), status='replace', action='write', iostat=ios)
    if (ios /= 0) return
    write(io_unit, '(A,1X,I0)') 'stage12_prescribed_velocity_allocated_status', allocated_status
    write(io_unit, '(A,1X,I0)') 'stage12_prescribed_velocity_point_count_status', point_count_status
    write(io_unit, '(A,1X,I0)') 'stage12_prescribed_velocity_zero_status', zero_velocity_status
    write(io_unit, '(A,1X,I0)') 'stage12_prescribed_velocity_constant_status', constant_velocity_status
    write(io_unit, '(A,1X,I0)') 'stage12_prescribed_velocity_norm_finite_status', velocity_norm_finite_status
    write(io_unit, '(A,1X,I0)') 'stage12_prescribed_velocity_valid_flag_status', velocity_valid_flag_status
    write(io_unit, '(A,1X,I0)') 'stage12_prescribed_velocity_clear_status', clear_status
    write(io_unit, '(A,1X,I0)') 'stage12_prescribed_velocity_status', prescribed_velocity_status
    close(io_unit)
  end subroutine stage12_prescribed_velocity_write_diagnostics

  subroutine update_status_values()
    if (stage12_prescribed_velocity_is_allocated()) then
      allocated_status = 1
      if (n_points > 0 .and. size(v_f, 1) == n_points .and. size(v_f, 2) == 3 .and. &
          size(velocity_norm) == n_points .and. size(velocity_valid_flag) == n_points) then
        point_count_status = 1
      else
        point_count_status = 0
      end if
      call update_zero_velocity_status()
      call update_velocity_norm_finite_status()
      call update_velocity_valid_flag_status()
    else
      allocated_status = 0
      point_count_status = 0
      zero_velocity_status = 0
      velocity_norm_finite_status = 0
      velocity_valid_flag_status = 0
    end if

    if (allocated_status == 1 .and. point_count_status == 1 .and. zero_velocity_status == 1 .and. &
        constant_velocity_status == 1 .and. velocity_norm_finite_status == 1 .and. &
        velocity_valid_flag_status == 1 .and. clear_status == 1 .and. &
        no_force_computation_status == 1 .and. no_eulerian_force_density_status == 1 .and. &
        no_rhs_injection_status == 1 .and. no_ibm_spreading_status == 1 .and. &
        no_feedback_application_status == 1 .and. no_twoway_force_status == 1 .and. &
        no_structure_advance_status == 1 .and. no_fluid_field_access_status == 1 .and. &
        no_fluid_field_modification_status == 1) then
      prescribed_velocity_status = 1
    else
      prescribed_velocity_status = 0
    end if
  end subroutine update_status_values

  subroutine update_zero_velocity_status()
    integer :: i
    integer :: j

    zero_velocity_status = 1
    do i = 1, n_points
      do j = 1, 3
        if (.not. is_finite(v_f(i, j)) .or. abs(v_f(i, j)) > 0.0_mytype) zero_velocity_status = 0
      end do
      if (.not. is_finite(velocity_norm(i)) .or. abs(velocity_norm(i)) > 0.0_mytype) zero_velocity_status = 0
    end do
  end subroutine update_zero_velocity_status

  subroutine update_velocity_norm_finite_status()
    integer :: i

    velocity_norm_finite_status = 1
    do i = 1, n_points
      if (.not. is_finite(velocity_norm(i))) velocity_norm_finite_status = 0
    end do
  end subroutine update_velocity_norm_finite_status

  subroutine update_velocity_valid_flag_status()
    integer :: i

    velocity_valid_flag_status = 1
    do i = 1, n_points
      if (velocity_valid_flag(i) /= 0 .and. velocity_valid_flag(i) /= 1) velocity_valid_flag_status = 0
    end do
  end subroutine update_velocity_valid_flag_status

  logical function is_finite(value)
    real(mytype), intent(in) :: value

    is_finite = (value == value) .and. (abs(value) < huge(value))
  end function is_finite

end module fibre_stage12_prescribed_velocity
