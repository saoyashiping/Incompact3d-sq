module fibre_geometry

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t

  implicit none
  private

  public :: fibre_segment_tangent
  public :: fibre_segment_length
  public :: fibre_second_derivative_internal
  public :: fibre_end_to_end_distance
  public :: fibre_center_of_mass

contains

  subroutine fibre_segment_tangent(fibre, tangent)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: tangent(3, fibre%nl - 1)

    integer :: i

    do i = 1, fibre%nl - 1
      tangent(:, i) = (fibre%x(:, i + 1) - fibre%x(:, i)) / fibre%ds
    end do
  end subroutine fibre_segment_tangent

  subroutine fibre_segment_length(fibre, seg_length)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: seg_length(fibre%nl - 1)

    integer :: i

    do i = 1, fibre%nl - 1
      seg_length(i) = sqrt(sum((fibre%x(:, i + 1) - fibre%x(:, i))**2))
    end do
  end subroutine fibre_segment_length

  subroutine fibre_second_derivative_internal(fibre, x_ss)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: x_ss(3, fibre%nl)

    integer :: i

    x_ss = 0._mytype

    do i = 2, fibre%nl - 1
      x_ss(:, i) = (fibre%x(:, i + 1) - 2._mytype * fibre%x(:, i) + fibre%x(:, i - 1)) / (fibre%ds**2)
    end do
  end subroutine fibre_second_derivative_internal

  function fibre_end_to_end_distance(fibre) result(distance)
    type(fibre_t), intent(in) :: fibre
    real(mytype) :: distance

    distance = sqrt(sum((fibre%x(:, fibre%nl) - fibre%x(:, 1))**2))
  end function fibre_end_to_end_distance

  function fibre_center_of_mass(fibre) result(x_cm)
    type(fibre_t), intent(in) :: fibre
    real(mytype) :: x_cm(3)
    real(mytype) :: integral(3)

    if (fibre%length <= 0._mytype) then
      x_cm = 0._mytype
      return
    end if

    integral = 0.5_mytype * fibre%x(:, 1) * fibre%ds
    if (fibre%nl > 2) then
      integral = integral + sum(fibre%x(:, 2:fibre%nl - 1), dim=2) * fibre%ds
    end if
    integral = integral + 0.5_mytype * fibre%x(:, fibre%nl) * fibre%ds

    x_cm = integral / fibre%length
  end function fibre_center_of_mass

end module fibre_geometry
