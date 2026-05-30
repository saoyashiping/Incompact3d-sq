module fibre_stage12_force_buffer
  implicit none
  private

  integer, parameter :: mytype = kind(1.0d0)

  integer :: n_points = 0
  real(mytype), allocatable :: f_fs_cand(:,:)
  real(mytype), allocatable :: force_norm(:)
  integer, allocatable :: force_valid_flag(:)

  integer :: allocated_status = 0
  integer :: point_count_status = 0
  integer :: zero_force_status = 0
  integer :: force_norm_finite_status = 0
  integer :: force_valid_flag_status = 0
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
  integer :: force_buffer_status = 0

  public :: stage12_force_buffer_init
  public :: stage12_force_buffer_clear
  public :: stage12_force_buffer_finalize
  public :: stage12_force_buffer_is_allocated
  public :: stage12_force_buffer_get_count
  public :: stage12_force_buffer_get_status_values
  public :: stage12_force_buffer_write_diagnostics

contains

  subroutine stage12_force_buffer_init(requested_points)
    integer, intent(in) :: requested_points
    integer :: alloc_status

    call stage12_force_buffer_finalize()

    if (requested_points <= 0) then
      call update_status_values()
      return
    end if

    n_points = requested_points
    allocate(f_fs_cand(n_points, 3), stat=alloc_status)
    if (alloc_status /= 0) then
      n_points = 0
      call update_status_values()
      return
    end if

    allocate(force_norm(n_points), stat=alloc_status)
    if (alloc_status /= 0) then
      call stage12_force_buffer_finalize()
      return
    end if

    allocate(force_valid_flag(n_points), stat=alloc_status)
    if (alloc_status /= 0) then
      call stage12_force_buffer_finalize()
      return
    end if

    f_fs_cand(:, :) = 0.0_mytype
    force_norm(:) = 0.0_mytype
    force_valid_flag(:) = 1
    clear_status = 0
    call update_status_values()
  end subroutine stage12_force_buffer_init

  subroutine stage12_force_buffer_clear()
    if (stage12_force_buffer_is_allocated()) then
      f_fs_cand(:, :) = 0.0_mytype
      force_norm(:) = 0.0_mytype
      force_valid_flag(:) = 1
      clear_status = 1
    else
      clear_status = 0
    end if
    call update_status_values()
  end subroutine stage12_force_buffer_clear

  subroutine stage12_force_buffer_finalize()
    if (allocated(f_fs_cand)) deallocate(f_fs_cand)
    if (allocated(force_norm)) deallocate(force_norm)
    if (allocated(force_valid_flag)) deallocate(force_valid_flag)
    n_points = 0
    allocated_status = 0
    point_count_status = 0
    zero_force_status = 0
    force_norm_finite_status = 0
    force_valid_flag_status = 0
    clear_status = 0
    force_buffer_status = 0
  end subroutine stage12_force_buffer_finalize

  logical function stage12_force_buffer_is_allocated()
    stage12_force_buffer_is_allocated = allocated(f_fs_cand) .and. allocated(force_norm) .and. allocated(force_valid_flag)
  end function stage12_force_buffer_is_allocated

  integer function stage12_force_buffer_get_count()
    stage12_force_buffer_get_count = n_points
  end function stage12_force_buffer_get_count

  subroutine stage12_force_buffer_get_status_values(out_allocated_status, out_point_count_status, &
                                                    out_zero_force_status, out_force_norm_finite_status, &
                                                    out_force_valid_flag_status, out_clear_status, &
                                                    out_no_force_computation_status, &
                                                    out_no_eulerian_force_density_status, &
                                                    out_no_rhs_injection_status, out_no_ibm_spreading_status, &
                                                    out_no_feedback_application_status, out_no_twoway_force_status, &
                                                    out_no_structure_advance_status, out_no_fluid_field_access_status, &
                                                    out_no_fluid_field_modification_status, out_force_buffer_status)
    integer, intent(out) :: out_allocated_status
    integer, intent(out) :: out_point_count_status
    integer, intent(out) :: out_zero_force_status
    integer, intent(out) :: out_force_norm_finite_status
    integer, intent(out) :: out_force_valid_flag_status
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
    integer, intent(out) :: out_force_buffer_status

    call update_status_values()

    out_allocated_status = allocated_status
    out_point_count_status = point_count_status
    out_zero_force_status = zero_force_status
    out_force_norm_finite_status = force_norm_finite_status
    out_force_valid_flag_status = force_valid_flag_status
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
    out_force_buffer_status = force_buffer_status
  end subroutine stage12_force_buffer_get_status_values

  subroutine stage12_force_buffer_write_diagnostics(filename)
    character(len=*), intent(in), optional :: filename
    character(len=256) :: output_file
    integer :: io_unit
    integer :: ios

    if (present(filename)) then
      output_file = filename
    else
      output_file = 'stage12_outputs/fibre_stage12_force_buffer.dat'
    end if

    call update_status_values()
    open(newunit=io_unit, file=trim(output_file), status='replace', action='write', iostat=ios)
    if (ios /= 0) return
    write(io_unit, '(A,1X,I0)') 'stage12_force_buffer_allocated_status', allocated_status
    write(io_unit, '(A,1X,I0)') 'stage12_force_buffer_point_count_status', point_count_status
    write(io_unit, '(A,1X,I0)') 'stage12_force_buffer_zero_force_status', zero_force_status
    write(io_unit, '(A,1X,I0)') 'stage12_force_buffer_force_norm_finite_status', force_norm_finite_status
    write(io_unit, '(A,1X,I0)') 'stage12_force_buffer_force_valid_flag_status', force_valid_flag_status
    write(io_unit, '(A,1X,I0)') 'stage12_force_buffer_clear_status', clear_status
    write(io_unit, '(A,1X,I0)') 'stage12_force_buffer_status', force_buffer_status
    close(io_unit)
  end subroutine stage12_force_buffer_write_diagnostics

  subroutine update_status_values()
    if (stage12_force_buffer_is_allocated()) then
      allocated_status = 1
      if (n_points > 0 .and. size(f_fs_cand, 1) == n_points .and. size(f_fs_cand, 2) == 3 .and. &
          size(force_norm) == n_points .and. size(force_valid_flag) == n_points) then
        point_count_status = 1
      else
        point_count_status = 0
      end if
      call update_zero_force_status()
      call update_force_norm_finite_status()
      call update_force_valid_flag_status()
    else
      allocated_status = 0
      point_count_status = 0
      zero_force_status = 0
      force_norm_finite_status = 0
      force_valid_flag_status = 0
    end if

    if (allocated_status == 1 .and. point_count_status == 1 .and. zero_force_status == 1 .and. &
        force_norm_finite_status == 1 .and. force_valid_flag_status == 1 .and. clear_status == 1 .and. &
        no_force_computation_status == 1 .and. no_eulerian_force_density_status == 1 .and. &
        no_rhs_injection_status == 1 .and. no_ibm_spreading_status == 1 .and. &
        no_feedback_application_status == 1 .and. no_twoway_force_status == 1 .and. &
        no_structure_advance_status == 1 .and. no_fluid_field_access_status == 1 .and. &
        no_fluid_field_modification_status == 1) then
      force_buffer_status = 1
    else
      force_buffer_status = 0
    end if
  end subroutine update_status_values

  subroutine update_zero_force_status()
    integer :: i
    integer :: j

    zero_force_status = 1
    do i = 1, n_points
      do j = 1, 3
        if (.not. is_finite(f_fs_cand(i, j)) .or. abs(f_fs_cand(i, j)) > 0.0_mytype) zero_force_status = 0
      end do
      if (.not. is_finite(force_norm(i)) .or. abs(force_norm(i)) > 0.0_mytype) zero_force_status = 0
    end do
  end subroutine update_zero_force_status

  subroutine update_force_norm_finite_status()
    integer :: i

    force_norm_finite_status = 1
    do i = 1, n_points
      if (.not. is_finite(force_norm(i))) force_norm_finite_status = 0
    end do
  end subroutine update_force_norm_finite_status

  subroutine update_force_valid_flag_status()
    integer :: i

    force_valid_flag_status = 1
    do i = 1, n_points
      if (force_valid_flag(i) /= 0 .and. force_valid_flag(i) /= 1) force_valid_flag_status = 0
    end do
  end subroutine update_force_valid_flag_status

  logical function is_finite(value)
    real(mytype), intent(in) :: value

    is_finite = (value == value) .and. (abs(value) < huge(value))
  end function is_finite

end module fibre_stage12_force_buffer
