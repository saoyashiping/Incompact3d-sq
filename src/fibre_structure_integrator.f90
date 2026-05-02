module fibre_structure_integrator

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_bending_operator, only : compute_fibre_d4_freefree
  use fibre_implicit_bending, only : build_freefree_d4_matrix
  use fibre_tension_solver, only : solve_tension_freefree

  implicit none
  private

  public :: build_tension_stiffness_matrix_freefree
  public :: build_structure_advance_matrix_freefree
  public :: advance_fibre_structure_freefree

contains

  subroutine build_tension_stiffness_matrix_freefree(fibre, tension, dtmat)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(in) :: tension(fibre%nl - 1)
    real(mytype), intent(out) :: dtmat(fibre%nl, fibre%nl)

    integer :: i, nseg
    real(mytype) :: inv_ds2

    nseg = fibre%nl - 1
    inv_ds2 = 1._mytype / (fibre%ds**2)

    dtmat = 0._mytype

    dtmat(1, 1) = -tension(1) * inv_ds2
    dtmat(1, 2) = tension(1) * inv_ds2

    do i = 2, fibre%nl - 1
      dtmat(i, i - 1) = tension(i - 1) * inv_ds2
      dtmat(i, i) = -(tension(i - 1) + tension(i)) * inv_ds2
      dtmat(i, i + 1) = tension(i) * inv_ds2
    end do

    dtmat(fibre%nl, fibre%nl - 1) = tension(nseg) * inv_ds2
    dtmat(fibre%nl, fibre%nl) = -tension(nseg) * inv_ds2
  end subroutine build_tension_stiffness_matrix_freefree

  subroutine build_structure_advance_matrix_freefree(fibre, dt, tension, amat)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: tension(fibre%nl - 1)
    real(mytype), intent(out) :: amat(fibre%nl, fibre%nl)

    real(mytype) :: dtmat(fibre%nl, fibre%nl)
    real(mytype) :: d4mat(fibre%nl, fibre%nl)
    integer :: i

    if (fibre%rho_tilde <= 0._mytype) then
      error stop 'build_structure_advance_matrix_freefree: rho_tilde must be > 0'
    end if

    call build_tension_stiffness_matrix_freefree(fibre, tension, dtmat)
    call build_freefree_d4_matrix(fibre, d4mat)

    amat = -dt * dt * dtmat + dt * dt * (fibre%gamma / fibre%rho_tilde) * d4mat
    do i = 1, fibre%nl
      amat(i, i) = amat(i, i) + 1._mytype
    end do
  end subroutine build_structure_advance_matrix_freefree

  subroutine advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
    type(fibre_t), intent(inout) :: fibre
    real(mytype), intent(in) :: dt
    real(mytype), intent(out) :: tension_residual
    real(mytype), intent(out) :: relative_tension_residual
    integer, intent(out) :: ierr

    real(mytype) :: x_n(3, fibre%nl), x_nm1(3, fibre%nl), x_star(3, fibre%nl)
    real(mytype) :: x_new(3, fibre%nl), rhs(3, fibre%nl)
    real(mytype) :: d4x_star(3, fibre%nl)
    real(mytype) :: a_non_tension(3, fibre%nl)
    real(mytype) :: tension(fibre%nl - 1)
    real(mytype) :: amat(fibre%nl, fibre%nl)
    real(mytype) :: sol(fibre%nl)
    type(fibre_t) :: fibre_tmp
    integer :: comp, solve_status

    if (fibre%rho_tilde <= 0._mytype) then
      error stop 'advance_fibre_structure_freefree: rho_tilde must be > 0'
    end if

    x_n = fibre%x
    x_nm1 = fibre%x_old
    x_star = 2._mytype * x_n - x_nm1

    fibre_tmp = fibre
    fibre_tmp%x = x_star
    call compute_fibre_d4_freefree(fibre_tmp, d4x_star)
    a_non_tension = -(fibre%gamma / fibre%rho_tilde) * d4x_star + fibre%f_ext / fibre%rho_tilde

    call solve_tension_freefree(fibre, dt, a_non_tension, tension, tension_residual, relative_tension_residual, solve_status)
    if (solve_status /= 0) then
      ierr = solve_status
      return
    end if

    call build_structure_advance_matrix_freefree(fibre, dt, tension, amat)

    rhs = 2._mytype * x_n - x_nm1 + dt * dt * (fibre%f_ext / fibre%rho_tilde)

    ierr = 0
    do comp = 1, 3
      call solve_dense_linear_system(amat, rhs(comp, :), sol, solve_status)
      if (solve_status /= 0) then
        ierr = solve_status
        return
      end if
      x_new(comp, :) = sol
    end do

    fibre%x_old = x_n
    fibre%x = x_new
    fibre%v = (fibre%x - fibre%x_old) / dt
    fibre%tension = tension
  end subroutine advance_fibre_structure_freefree

  subroutine solve_dense_linear_system(ain, bin, xout, ierr)
    real(mytype), intent(in) :: ain(:,:), bin(:)
    real(mytype), intent(out) :: xout(:)
    integer, intent(out) :: ierr

    real(mytype), allocatable :: a(:,:), b(:), temp_row(:)
    real(mytype) :: pivot_val, factor, eps
    integer :: n, i, j, k, pivot_row

    n = size(bin)
    allocate(a(n, n), b(n), temp_row(n))

    a = ain
    b = bin

    eps = epsilon(1.0_mytype)
    ierr = 0

    do k = 1, n - 1
      pivot_row = k
      pivot_val = abs(a(k, k))

      do i = k + 1, n
        if (abs(a(i, k)) > pivot_val) then
          pivot_val = abs(a(i, k))
          pivot_row = i
        end if
      end do

      if (pivot_val <= eps) then
        ierr = 1
        deallocate(a, b, temp_row)
        return
      end if

      if (pivot_row /= k) then
        temp_row = a(k, :)
        a(k, :) = a(pivot_row, :)
        a(pivot_row, :) = temp_row

        factor = b(k)
        b(k) = b(pivot_row)
        b(pivot_row) = factor
      end if

      do i = k + 1, n
        factor = a(i, k) / a(k, k)
        a(i, k) = 0._mytype
        do j = k + 1, n
          a(i, j) = a(i, j) - factor * a(k, j)
        end do
        b(i) = b(i) - factor * b(k)
      end do
    end do

    if (abs(a(n, n)) <= eps) then
      ierr = 1
      deallocate(a, b, temp_row)
      return
    end if

    xout(n) = b(n) / a(n, n)
    do i = n - 1, 1, -1
      xout(i) = b(i)
      do j = i + 1, n
        xout(i) = xout(i) - a(i, j) * xout(j)
      end do
      if (abs(a(i, i)) <= eps) then
        ierr = 1
        deallocate(a, b, temp_row)
        return
      end if
      xout(i) = xout(i) / a(i, i)
    end do

    deallocate(a, b, temp_row)
  end subroutine solve_dense_linear_system

end module fibre_structure_integrator
