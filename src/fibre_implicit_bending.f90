module fibre_implicit_bending

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t

  implicit none
  private

  public :: build_freefree_d4_matrix
  public :: build_implicit_bending_matrix
  public :: implicit_bending_relax_step

contains

  subroutine build_freefree_d4_matrix(fibre, d4mat)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(out) :: d4mat(fibre%nl, fibre%nl)

    integer :: i
    real(mytype) :: ds4

    if (fibre%nl < 5) then
      error stop 'build_freefree_d4_matrix: nl must be >= 5'
    end if

    d4mat = 0._mytype
    ds4 = fibre%ds**4

    d4mat(1, 1) = 2._mytype
    d4mat(1, 2) = -4._mytype
    d4mat(1, 3) = 2._mytype

    d4mat(2, 1) = -2._mytype
    d4mat(2, 2) = 5._mytype
    d4mat(2, 3) = -4._mytype
    d4mat(2, 4) = 1._mytype

    do i = 3, fibre%nl - 2
      d4mat(i, i - 2) = 1._mytype
      d4mat(i, i - 1) = -4._mytype
      d4mat(i, i) = 6._mytype
      d4mat(i, i + 1) = -4._mytype
      d4mat(i, i + 2) = 1._mytype
    end do

    d4mat(fibre%nl - 1, fibre%nl - 3) = 1._mytype
    d4mat(fibre%nl - 1, fibre%nl - 2) = -4._mytype
    d4mat(fibre%nl - 1, fibre%nl - 1) = 5._mytype
    d4mat(fibre%nl - 1, fibre%nl) = -2._mytype

    d4mat(fibre%nl, fibre%nl - 2) = 2._mytype
    d4mat(fibre%nl, fibre%nl - 1) = -4._mytype
    d4mat(fibre%nl, fibre%nl) = 2._mytype

    d4mat = d4mat / ds4
  end subroutine build_freefree_d4_matrix

  subroutine build_implicit_bending_matrix(fibre, dt, amat)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(in) :: dt
    real(mytype), intent(out) :: amat(fibre%nl, fibre%nl)

    real(mytype) :: d4mat(fibre%nl, fibre%nl)
    integer :: i

    call build_freefree_d4_matrix(fibre, d4mat)

    amat = dt * fibre%gamma * d4mat
    do i = 1, fibre%nl
      amat(i, i) = amat(i, i) + 1._mytype
    end do
  end subroutine build_implicit_bending_matrix

  subroutine implicit_bending_relax_step(fibre, dt, ierr)
    type(fibre_t), intent(inout) :: fibre
    real(mytype), intent(in) :: dt
    integer, intent(out) :: ierr

    real(mytype) :: amat(fibre%nl, fibre%nl)
    real(mytype) :: rhs(fibre%nl), sol(fibre%nl)
    real(mytype) :: x_old(3, fibre%nl)
    integer :: comp, solve_status

    call build_implicit_bending_matrix(fibre, dt, amat)

    x_old = fibre%x
    ierr = 0

    do comp = 1, 3
      rhs = x_old(comp, :)
      call solve_dense_linear_system(amat, rhs, sol, solve_status)
      if (solve_status /= 0) then
        ierr = solve_status
        return
      end if
      fibre%x(comp, :) = sol
    end do
  end subroutine implicit_bending_relax_step

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

end module fibre_implicit_bending
