module fibre_bending_operator

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_boundary_freefree, only : fibre_compute_freefree_ghost_points

  implicit none
  private

  public :: compute_fibre_d4_freefree
  public :: compute_fibre_bending_force_freefree

contains

  subroutine compute_fibre_d4_freefree(fibre, d4x)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: d4x(3, fibre%nl)

    integer :: i
    real(mytype) :: x_g_l1(3), x_g_l2(3), x_g_r1(3), x_g_r2(3)
    real(mytype) :: ds4

    if (fibre%nl < 5) then
      error stop 'compute_fibre_d4_freefree: nl must be >= 5'
    end if

    call fibre_compute_freefree_ghost_points(fibre, x_g_l1, x_g_l2, x_g_r1, x_g_r2)

    ds4 = fibre%ds**4

    d4x(:, 1) = (fibre%x(:, 3) - 4._mytype * fibre%x(:, 2) + 6._mytype * fibre%x(:, 1) - 4._mytype * x_g_l1 + x_g_l2) / ds4
    d4x(:, 2) = (fibre%x(:, 4) - 4._mytype * fibre%x(:, 3) + 6._mytype * fibre%x(:, 2) - 4._mytype * fibre%x(:, 1) + x_g_l1) / ds4

    do i = 3, fibre%nl - 2
      d4x(:, i) = (fibre%x(:, i + 2) - 4._mytype * fibre%x(:, i + 1) + 6._mytype * fibre%x(:, i) - 4._mytype * fibre%x(:, i - 1) + &
                   fibre%x(:, i - 2)) / ds4
    end do

    d4x(:, fibre%nl - 1) = (x_g_r1 - 4._mytype * fibre%x(:, fibre%nl) + &
                            6._mytype * fibre%x(:, fibre%nl - 1) - 4._mytype * fibre%x(:, fibre%nl - 2) + &
                            fibre%x(:, fibre%nl - 3)) / ds4

    d4x(:, fibre%nl) = (x_g_r2 - 4._mytype * x_g_r1 + 6._mytype * fibre%x(:, fibre%nl) - 4._mytype * fibre%x(:, fibre%nl - 1) + &
                        fibre%x(:, fibre%nl - 2)) / ds4
  end subroutine compute_fibre_d4_freefree

  subroutine compute_fibre_bending_force_freefree(fibre, fb)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: fb(3, fibre%nl)

    real(mytype) :: d4x(3, fibre%nl)

    call compute_fibre_d4_freefree(fibre, d4x)
    fb = -fibre%gamma * d4x
  end subroutine compute_fibre_bending_force_freefree

end module fibre_bending_operator
