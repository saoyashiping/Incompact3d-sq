module fibre_tension_solver

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t

  implicit none
  private

  public :: assemble_tension_system_freefree
  public :: solve_tension_freefree
  public :: compute_tension_acceleration_freefree
  public :: apply_tension_operator_freefree

contains

  subroutine assemble_tension_system_freefree(fibre, dt, a_non_tension, amat, rhs)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: a_non_tension(3, fibre%nl)
    real(mytype), intent(out) :: amat(fibre%nl - 1, fibre%nl - 1)
    real(mytype), intent(out) :: rhs(fibre%nl - 1)

    integer :: nseg, j
    real(mytype) :: x_star(3, fibre%nl)
    real(mytype) :: t_n(3, fibre%nl - 1), t_old(3, fibre%nl - 1), t_star(3, fibre%nl - 1)
    real(mytype) :: vs(3, fibre%nl - 1)
    real(mytype) :: constraint_correction, non_tension_gradient_term
    real(mytype), allocatable :: tau_basis(:), op_val(:)

    nseg = fibre%nl - 1
    if (size(a_non_tension, 1) /= 3) error stop 'assemble_tension_system_freefree: a_non_tension first dimension must be 3'
    x_star = 2._mytype * fibre%x - fibre%x_old

    do j = 1, nseg
      t_n(:, j) = (fibre%x(:, j + 1) - fibre%x(:, j)) / fibre%ds
      t_old(:, j) = (fibre%x_old(:, j + 1) - fibre%x_old(:, j)) / fibre%ds
      t_star(:, j) = (x_star(:, j + 1) - x_star(:, j)) / fibre%ds
      vs(:, j) = (fibre%v(:, j + 1) - fibre%v(:, j)) / fibre%ds

      constraint_correction = (1._mytype - 2._mytype * dot_product(t_n(:, j), t_n(:, j)) + &
                               dot_product(t_old(:, j), t_old(:, j))) / (2._mytype * dt * dt)

      non_tension_gradient_term = dot_product(t_star(:, j), (a_non_tension(:, j + 1) - a_non_tension(:, j)) / fibre%ds)

      rhs(j) = constraint_correction - dot_product(vs(:, j), vs(:, j)) - non_tension_gradient_term
    end do

    amat = 0._mytype
    allocate(tau_basis(nseg), op_val(nseg))

    do j = 1, nseg
      tau_basis = 0._mytype
      tau_basis(j) = 1._mytype
      call apply_tension_operator_from_tstar(t_star, fibre%ds, tau_basis, op_val)
      amat(:, j) = op_val
    end do

    deallocate(tau_basis, op_val)
  end subroutine assemble_tension_system_freefree

  subroutine solve_tension_freefree(fibre, dt, a_non_tension, tension, residual_norm, relative_residual, ierr)
    type(fibre_t), intent(inout) :: fibre
    real(mytype), intent(in) :: dt
    real(mytype), intent(in) :: a_non_tension(3, fibre%nl)
    real(mytype), intent(out) :: tension(fibre%nl - 1)
    real(mytype), intent(out) :: residual_norm
    real(mytype), intent(out) :: relative_residual
    integer, intent(out) :: ierr

    real(mytype) :: amat(fibre%nl - 1, fibre%nl - 1)
    real(mytype) :: rhs(fibre%nl - 1), residual(fibre%nl - 1)
    real(mytype) :: rhs_norm

    call assemble_tension_system_freefree(fibre, dt, a_non_tension, amat, rhs)
    call solve_dense_linear_system(amat, rhs, tension, ierr)

    if (ierr /= 0) then
      residual_norm = huge(1._mytype)
      relative_residual = huge(1._mytype)
      return
    end if

    residual = matmul(amat, tension) - rhs
    residual_norm = sqrt(sum(residual**2))
    rhs_norm = sqrt(sum(rhs**2))
    relative_residual = residual_norm / max(rhs_norm, tiny(1._mytype))

    fibre%tension = tension
  end subroutine solve_tension_freefree

  subroutine compute_tension_acceleration_freefree(fibre, tension, a_tension)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(in) :: tension(fibre%nl - 1)
    real(mytype), intent(out) :: a_tension(3, fibre%nl)

    integer :: j
    real(mytype) :: x_star(3, fibre%nl)
    real(mytype) :: t_star(3, fibre%nl - 1)

    if (size(a_tension, 1) /= 3) error stop 'compute_tension_acceleration_freefree: a_tension first dimension must be 3'
    x_star = 2._mytype * fibre%x - fibre%x_old

    do j = 1, fibre%nl - 1
      t_star(:, j) = (x_star(:, j + 1) - x_star(:, j)) / fibre%ds
    end do

    call compute_tension_acceleration_from_tstar(t_star, fibre%ds, tension, a_tension)
  end subroutine compute_tension_acceleration_freefree

  subroutine apply_tension_operator_freefree(fibre, tension, op_val)
    type(fibre_t), intent(in) :: fibre
    real(mytype), intent(in) :: tension(fibre%nl - 1)
    real(mytype), intent(out) :: op_val(fibre%nl - 1)

    integer :: j
    real(mytype) :: x_star(3, fibre%nl)
    real(mytype) :: t_star(3, fibre%nl - 1)

    x_star = 2._mytype * fibre%x - fibre%x_old

    do j = 1, fibre%nl - 1
      t_star(:, j) = (x_star(:, j + 1) - x_star(:, j)) / fibre%ds
    end do

    call apply_tension_operator_from_tstar(t_star, fibre%ds, tension, op_val)
  end subroutine apply_tension_operator_freefree

  subroutine apply_tension_operator_from_tstar(t_star, ds, tension, op_val)
    real(mytype), intent(in) :: t_star(:, :)
    real(mytype), intent(in) :: ds
    real(mytype), intent(in) :: tension(:)
    real(mytype), intent(out) :: op_val(:)

    real(mytype) :: a_tension(3, size(tension) + 1)
    integer :: j, nseg

    if (size(t_star, 1) /= 3) error stop 'apply_tension_operator_from_tstar: t_star first dimension must be 3'
    nseg = size(t_star, 2)
    if (size(tension) /= nseg) error stop 'apply_tension_operator_from_tstar: tension size mismatch'
    if (size(op_val) /= nseg) error stop 'apply_tension_operator_from_tstar: op_val size mismatch'

    call compute_tension_acceleration_from_tstar(t_star, ds, tension, a_tension)

    do j = 1, size(tension)
      op_val(j) = dot_product(t_star(:, j), (a_tension(:, j + 1) - a_tension(:, j)) / ds)
    end do
  end subroutine apply_tension_operator_from_tstar

  subroutine compute_tension_acceleration_from_tstar(t_star, ds, tension, a_tension)
    real(mytype), intent(in) :: t_star(:, :)
    real(mytype), intent(in) :: ds
    real(mytype), intent(in) :: tension(:)
    real(mytype), intent(out) :: a_tension(:, :)

    integer :: i, nseg

    if (size(t_star, 1) /= 3) error stop 'compute_tension_acceleration_from_tstar: t_star first dimension must be 3'
    nseg = size(t_star, 2)
    if (size(tension) /= nseg) error stop 'compute_tension_acceleration_from_tstar: tension size mismatch'
    if (size(a_tension, 1) /= 3) error stop 'compute_tension_acceleration_from_tstar: a_tension first dimension must be 3'
    if (size(a_tension, 2) /= nseg + 1) error stop 'compute_tension_acceleration_from_tstar: a_tension second dimension mismatch'

    a_tension = 0._mytype

    a_tension(:, 1) = tension(1) * t_star(:, 1) / ds

    do i = 2, nseg
      a_tension(:, i) = (tension(i) * t_star(:, i) - tension(i - 1) * t_star(:, i - 1)) / ds
    end do

    a_tension(:, nseg + 1) = -tension(nseg) * t_star(:, nseg) / ds
  end subroutine compute_tension_acceleration_from_tstar

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

end module fibre_tension_solver
