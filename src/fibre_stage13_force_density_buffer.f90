module fibre_stage13_force_density_buffer
  use decomp_2d_constants, only : mytype
  implicit none
  private

  integer :: nx_buf = 0
  integer :: ny_buf = 0
  integer :: nz_buf = 0
  real(mytype), allocatable :: fx_cand(:,:,:)
  real(mytype), allocatable :: fy_cand(:,:,:)
  real(mytype), allocatable :: fz_cand(:,:,:)
  real(mytype), allocatable :: force_density_norm(:,:,:)
  integer, allocatable :: force_density_valid_flag(:,:,:)

  integer :: allocated_status = 0
  integer :: shape_status = 0
  integer :: zero_force_density_status = 1
  integer :: force_density_norm_finite_status = 1
  integer :: force_density_valid_flag_status = 1
  integer :: clear_status = 0
  integer :: no_spreading_status = 1
  integer :: no_rhs_injection_status = 1
  integer :: no_ibm_spreading_status = 1
  integer :: no_feedback_application_status = 1
  integer :: no_twoway_force_status = 1
  integer :: no_structure_advance_status = 1
  integer :: no_fluid_field_access_status = 1
  integer :: no_fluid_field_modification_status = 1
  integer :: force_density_buffer_status = 0

  public :: stage13_force_density_buffer_init
  public :: stage13_force_density_buffer_clear
  public :: stage13_force_density_buffer_finalize
  public :: stage13_force_density_buffer_is_allocated
  public :: stage13_force_density_buffer_get_shape
  public :: stage13_force_density_buffer_get_status_values
  public :: stage13_force_density_buffer_write_diagnostics

contains

  subroutine stage13_force_density_buffer_init(nx, ny, nz)
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    integer :: alloc_status

    if (nx <= 0 .or. ny <= 0 .or. nz <= 0) then
      call stage13_force_density_buffer_finalize()
      shape_status = 0
      allocated_status = 0
      clear_status = 0
      call update_status()
      return
    end if

    call stage13_force_density_buffer_finalize()
    allocate(fx_cand(nx,ny,nz), fy_cand(nx,ny,nz), fz_cand(nx,ny,nz), &
             force_density_norm(nx,ny,nz), force_density_valid_flag(nx,ny,nz), stat=alloc_status)
    if (alloc_status /= 0) then
      call stage13_force_density_buffer_finalize()
      shape_status = 0
      allocated_status = 0
      clear_status = 0
      call update_status()
      return
    end if

    nx_buf = nx
    ny_buf = ny
    nz_buf = nz
    shape_status = 1
    allocated_status = 1
    fx_cand(:,:,:) = 0.0_mytype
    fy_cand(:,:,:) = 0.0_mytype
    fz_cand(:,:,:) = 0.0_mytype
    force_density_norm(:,:,:) = 0.0_mytype
    force_density_valid_flag(:,:,:) = 1
    clear_status = 1
    call update_status()
  end subroutine stage13_force_density_buffer_init

  subroutine stage13_force_density_buffer_clear()
    if (.not. stage13_force_density_buffer_is_allocated()) then
      clear_status = 0
      call update_status()
      return
    end if

    fx_cand(:,:,:) = 0.0_mytype
    fy_cand(:,:,:) = 0.0_mytype
    fz_cand(:,:,:) = 0.0_mytype
    force_density_norm(:,:,:) = 0.0_mytype
    force_density_valid_flag(:,:,:) = 1
    clear_status = 1
    call update_status()
  end subroutine stage13_force_density_buffer_clear

  subroutine stage13_force_density_buffer_finalize()
    if (allocated(fx_cand)) deallocate(fx_cand)
    if (allocated(fy_cand)) deallocate(fy_cand)
    if (allocated(fz_cand)) deallocate(fz_cand)
    if (allocated(force_density_norm)) deallocate(force_density_norm)
    if (allocated(force_density_valid_flag)) deallocate(force_density_valid_flag)
    nx_buf = 0
    ny_buf = 0
    nz_buf = 0
    allocated_status = 0
    shape_status = 0
    zero_force_density_status = 1
    force_density_norm_finite_status = 1
    force_density_valid_flag_status = 1
    clear_status = 0
    call update_status()
  end subroutine stage13_force_density_buffer_finalize

  logical function stage13_force_density_buffer_is_allocated()
    stage13_force_density_buffer_is_allocated = allocated(fx_cand) .and. allocated(fy_cand) .and. &
                                                allocated(fz_cand) .and. allocated(force_density_norm) .and. &
                                                allocated(force_density_valid_flag)
  end function stage13_force_density_buffer_is_allocated

  subroutine stage13_force_density_buffer_get_shape(nx, ny, nz)
    integer, intent(out) :: nx
    integer, intent(out) :: ny
    integer, intent(out) :: nz

    nx = nx_buf
    ny = ny_buf
    nz = nz_buf
  end subroutine stage13_force_density_buffer_get_shape

  subroutine stage13_force_density_buffer_get_status_values(allocated_out, shape_out, zero_out, norm_finite_out, &
                                                           valid_flag_out, clear_out, no_spreading_out, &
                                                           no_rhs_injection_out, no_ibm_spreading_out, &
                                                           no_feedback_application_out, no_twoway_force_out, &
                                                           no_structure_advance_out, no_fluid_field_access_out, &
                                                           no_fluid_field_modification_out, buffer_status_out)
    integer, intent(out) :: allocated_out
    integer, intent(out) :: shape_out
    integer, intent(out) :: zero_out
    integer, intent(out) :: norm_finite_out
    integer, intent(out) :: valid_flag_out
    integer, intent(out) :: clear_out
    integer, intent(out) :: no_spreading_out
    integer, intent(out) :: no_rhs_injection_out
    integer, intent(out) :: no_ibm_spreading_out
    integer, intent(out) :: no_feedback_application_out
    integer, intent(out) :: no_twoway_force_out
    integer, intent(out) :: no_structure_advance_out
    integer, intent(out) :: no_fluid_field_access_out
    integer, intent(out) :: no_fluid_field_modification_out
    integer, intent(out) :: buffer_status_out

    call update_status()
    allocated_out = allocated_status
    shape_out = shape_status
    zero_out = zero_force_density_status
    norm_finite_out = force_density_norm_finite_status
    valid_flag_out = force_density_valid_flag_status
    clear_out = clear_status
    no_spreading_out = no_spreading_status
    no_rhs_injection_out = no_rhs_injection_status
    no_ibm_spreading_out = no_ibm_spreading_status
    no_feedback_application_out = no_feedback_application_status
    no_twoway_force_out = no_twoway_force_status
    no_structure_advance_out = no_structure_advance_status
    no_fluid_field_access_out = no_fluid_field_access_status
    no_fluid_field_modification_out = no_fluid_field_modification_status
    buffer_status_out = force_density_buffer_status
  end subroutine stage13_force_density_buffer_get_status_values

  subroutine stage13_force_density_buffer_write_diagnostics()
    integer :: unit_id
    integer :: io_status

    open(newunit=unit_id, file='stage13_outputs/fibre_stage13_force_density_buffer.dat', &
         status='replace', action='write', iostat=io_status)
    if (io_status /= 0) return
    call update_status()
    write(unit_id,'(A,1X,I0)') 'stage13_force_density_buffer_allocated_status', allocated_status
    write(unit_id,'(A,1X,I0)') 'stage13_force_density_buffer_shape_status', shape_status
    write(unit_id,'(A,1X,I0)') 'stage13_force_density_buffer_zero_force_density_status', zero_force_density_status
    write(unit_id,'(A,1X,I0)') 'stage13_force_density_buffer_status', force_density_buffer_status
    close(unit_id)
  end subroutine stage13_force_density_buffer_write_diagnostics

  subroutine update_status()
    allocated_status = merge(1, 0, stage13_force_density_buffer_is_allocated())
    if (allocated_status == 1) then
      shape_status = merge(1, 0, nx_buf > 0 .and. ny_buf > 0 .and. nz_buf > 0)
      zero_force_density_status = merge(1, 0, all_zero_force_density())
      force_density_norm_finite_status = merge(1, 0, all_finite_rank3(force_density_norm))
      force_density_valid_flag_status = merge(1, 0, all_valid_flags())
    end if

    if (allocated_status == 1 .and. shape_status == 1 .and. zero_force_density_status == 1 .and. &
        force_density_norm_finite_status == 1 .and. force_density_valid_flag_status == 1 .and. &
        clear_status == 1 .and. no_spreading_status == 1 .and. no_rhs_injection_status == 1 .and. &
        no_ibm_spreading_status == 1 .and. no_feedback_application_status == 1 .and. &
        no_twoway_force_status == 1 .and. no_structure_advance_status == 1 .and. &
        no_fluid_field_access_status == 1 .and. no_fluid_field_modification_status == 1) then
      force_density_buffer_status = 1
    else
      force_density_buffer_status = 0
    end if
  end subroutine update_status

  logical function all_zero_force_density()
    all_zero_force_density = all(abs(fx_cand) == 0.0_mytype) .and. all(abs(fy_cand) == 0.0_mytype) .and. &
                             all(abs(fz_cand) == 0.0_mytype) .and. &
                             all(abs(force_density_norm) == 0.0_mytype)
  end function all_zero_force_density

  logical function all_valid_flags()
    all_valid_flags = all(force_density_valid_flag == 0 .or. force_density_valid_flag == 1)
  end function all_valid_flags

  logical function all_finite_rank3(values)
    real(mytype), intent(in) :: values(:,:,:)

    all_finite_rank3 = all(values == values) .and. all(abs(values) < huge(1.0_mytype))
  end function all_finite_rank3

end module fibre_stage13_force_density_buffer
