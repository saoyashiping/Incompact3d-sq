module fibre_diagnostics

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_geometry, only : fibre_segment_length, fibre_second_derivative_internal
  use fibre_geometry, only : fibre_end_to_end_distance, fibre_center_of_mass

  implicit none
  private

  type, public :: fibre_diagnostics_t
    real(mytype) :: max_segment_length_error = 0._mytype
    real(mytype) :: min_segment_length = 0._mytype
    real(mytype) :: max_segment_length = 0._mytype
    real(mytype) :: end_to_end_distance = 0._mytype
    real(mytype) :: center_of_mass(3) = 0._mytype
    real(mytype) :: max_curvature = 0._mytype
  end type fibre_diagnostics_t

  public :: fibre_compute_diagnostics
  public :: fibre_write_diagnostics

contains

  subroutine fibre_compute_diagnostics(fibre, diag)
    type(fibre_t), intent(in) :: fibre
    type(fibre_diagnostics_t), intent(out) :: diag

    real(mytype), allocatable :: seg_length(:)
    real(mytype), allocatable :: x_ss(:,:)
    integer :: i

    allocate(seg_length(fibre%nl - 1))
    allocate(x_ss(3, fibre%nl))

    call fibre_segment_length(fibre, seg_length)
    call fibre_second_derivative_internal(fibre, x_ss)

    diag%max_segment_length_error = maxval(abs(seg_length / fibre%ds - 1._mytype))
    diag%min_segment_length = minval(seg_length)
    diag%max_segment_length = maxval(seg_length)
    diag%end_to_end_distance = fibre_end_to_end_distance(fibre)
    diag%center_of_mass = fibre_center_of_mass(fibre)

    diag%max_curvature = 0._mytype
    do i = 2, fibre%nl - 1
      diag%max_curvature = max(diag%max_curvature, sqrt(sum(x_ss(:, i)**2)))
    end do

    deallocate(seg_length)
    deallocate(x_ss)
  end subroutine fibre_compute_diagnostics

  subroutine fibre_write_diagnostics(filename, fibre, diag)
    character(len=*), intent(in) :: filename
    type(fibre_t), intent(in) :: fibre
    type(fibre_diagnostics_t), intent(in) :: diag

    integer :: io_unit

    open(newunit=io_unit, file=filename, status='replace', action='write', form='formatted')
    write(io_unit, '(A)') '# free-free fibre geometry diagnostics'
    write(io_unit, '(A,I0)') 'nl = ', fibre%nl
    write(io_unit, '(A,ES24.16)') 'length = ', fibre%length
    write(io_unit, '(A,ES24.16)') 'ds = ', fibre%ds
    write(io_unit, '(A,ES24.16)') 'max_segment_length_error = ', diag%max_segment_length_error
    write(io_unit, '(A,ES24.16)') 'min_segment_length = ', diag%min_segment_length
    write(io_unit, '(A,ES24.16)') 'max_segment_length = ', diag%max_segment_length
    write(io_unit, '(A,ES24.16)') 'end_to_end_distance = ', diag%end_to_end_distance
    write(io_unit, '(A,3(ES24.16,1X))') 'center_of_mass = ', diag%center_of_mass
    write(io_unit, '(A,ES24.16)') 'max_curvature = ', diag%max_curvature
    close(io_unit)
  end subroutine fibre_write_diagnostics

end module fibre_diagnostics
