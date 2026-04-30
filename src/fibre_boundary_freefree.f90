module fibre_boundary_freefree

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t

  implicit none
  private

  type, public :: fibre_boundary_residual_t
    real(mytype) :: left_d2_norm = 0._mytype
    real(mytype) :: left_d3_norm = 0._mytype
    real(mytype) :: right_d2_norm = 0._mytype
    real(mytype) :: right_d3_norm = 0._mytype
    real(mytype) :: max_freefree_boundary_residual = 0._mytype
  end type fibre_boundary_residual_t

  public :: fibre_compute_freefree_ghost_points
  public :: fibre_compute_freefree_boundary_residual

contains

  subroutine fibre_compute_freefree_ghost_points(fibre, x_g_l1, x_g_l2, x_g_r1, x_g_r2)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: x_g_l1(3), x_g_l2(3), x_g_r1(3), x_g_r2(3)

    if (fibre%nl < 3) then
      error stop 'fibre_compute_freefree_ghost_points: nl must be >= 3'
    end if

    x_g_l1 = 2._mytype * fibre%x(:, 1) - fibre%x(:, 2)
    x_g_l2 = 4._mytype * fibre%x(:, 1) - 4._mytype * fibre%x(:, 2) + fibre%x(:, 3)

    x_g_r1 = 2._mytype * fibre%x(:, fibre%nl) - fibre%x(:, fibre%nl - 1)
    x_g_r2 = fibre%x(:, fibre%nl - 2) - 4._mytype * fibre%x(:, fibre%nl - 1) + 4._mytype * fibre%x(:, fibre%nl)
  end subroutine fibre_compute_freefree_ghost_points

  subroutine fibre_compute_freefree_boundary_residual(fibre, residual)
    type(fibre_t), intent(in) :: fibre
    type(fibre_boundary_residual_t), intent(out) :: residual

    real(mytype) :: x_g_l1(3), x_g_l2(3), x_g_r1(3), x_g_r2(3)
    real(mytype) :: d2_left(3), d3_left(3), d2_right(3), d3_right(3)

    call fibre_compute_freefree_ghost_points(fibre, x_g_l1, x_g_l2, x_g_r1, x_g_r2)

    d2_left = (x_g_l1 - 2._mytype * fibre%x(:, 1) + fibre%x(:, 2)) / (fibre%ds**2)
    d3_left = (x_g_l2 - 2._mytype * x_g_l1 + 2._mytype * fibre%x(:, 2) - fibre%x(:, 3)) / (2._mytype * fibre%ds**3)

    d2_right = (fibre%x(:, fibre%nl - 1) - 2._mytype * fibre%x(:, fibre%nl) + x_g_r1) / (fibre%ds**2)
    d3_right = (fibre%x(:, fibre%nl - 2) - 2._mytype * fibre%x(:, fibre%nl - 1) + 2._mytype * x_g_r1 - x_g_r2) / &
               (2._mytype * fibre%ds**3)

    residual%left_d2_norm = sqrt(sum(d2_left**2))
    residual%left_d3_norm = sqrt(sum(d3_left**2))
    residual%right_d2_norm = sqrt(sum(d2_right**2))
    residual%right_d3_norm = sqrt(sum(d3_right**2))

    residual%max_freefree_boundary_residual = max(residual%left_d2_norm, residual%left_d3_norm, residual%right_d2_norm, &
                                                  residual%right_d3_norm)
  end subroutine fibre_compute_freefree_boundary_residual

end module fibre_boundary_freefree
