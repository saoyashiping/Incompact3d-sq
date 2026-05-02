program fibre_tension_solver_check

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_parameters, only : fibre_init_straight_free_free
  use fibre_bending_operator, only : compute_fibre_bending_force_freefree
  use fibre_tension_solver, only : assemble_tension_system_freefree
  use fibre_tension_solver, only : solve_tension_freefree
  use fibre_tension_solver, only : apply_tension_operator_freefree

  implicit none

  type(fibre_t) :: fibre

  integer :: i, nseg, io_unit, ierr
  real(mytype) :: dt, pi, s, eps, amp

  real(mytype), allocatable :: a_non_tension(:,:), amat(:,:), rhs(:), tension(:), y1(:), y2(:), tau_test(:)
  real(mytype), allocatable :: fb(:,:)

  real(mytype) :: straight_rhs_norm, straight_tension_norm, straight_residual_norm, straight_relative_residual
  real(mytype) :: translation_rhs_norm, translation_tension_norm, translation_residual_norm
  real(mytype) :: dummy_relative_residual
  real(mytype) :: stretch_rhs_norm, stretch_tension_norm, stretch_residual_norm, stretch_relative_residual
  real(mytype) :: curved_rhs_norm, curved_tension_norm, curved_residual_norm, curved_relative_residual
  real(mytype) :: curved_endpoint_change_norm
  real(mytype) :: tension_matrix_apply_maxdiff
  real(mytype) :: left_before(3), left_after(3), right_before(3), right_after(3)

  dt = 1.0e-3_mytype

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  nseg = fibre%nl - 1
  allocate(a_non_tension(3, fibre%nl), amat(nseg, nseg), rhs(nseg), tension(nseg))

  a_non_tension = 0._mytype
  call assemble_tension_system_freefree(fibre, dt, a_non_tension, amat, rhs)
  call solve_tension_freefree(fibre, dt, a_non_tension, tension, straight_residual_norm, straight_relative_residual, ierr)
  if (ierr /= 0) error stop 'solve_tension_freefree failed in straight test'

  straight_rhs_norm = sqrt(sum(rhs**2))
  straight_tension_norm = sqrt(sum(tension**2))

  do i = 1, fibre%nl
    fibre%v(:, i) = [0.2_mytype, -0.1_mytype, 0.05_mytype]
    fibre%x_old(:, i) = fibre%x(:, i) - dt * fibre%v(:, i)
  end do

  a_non_tension = 0._mytype
  call assemble_tension_system_freefree(fibre, dt, a_non_tension, amat, rhs)
  call solve_tension_freefree(fibre, dt, a_non_tension, tension, translation_residual_norm, dummy_relative_residual, ierr)
  if (ierr /= 0) error stop 'solve_tension_freefree failed in translation test'

  translation_rhs_norm = sqrt(sum(rhs**2))
  translation_tension_norm = sqrt(sum(tension**2))

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  eps = 1.0e-3_mytype
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = (1._mytype + eps) * s
    fibre%x(2, i) = 0._mytype
    fibre%x(3, i) = 0._mytype
  end do
  fibre%x_old = fibre%x
  fibre%v = 0._mytype

  a_non_tension = 0._mytype
  call assemble_tension_system_freefree(fibre, dt, a_non_tension, amat, rhs)
  call solve_tension_freefree(fibre, dt, a_non_tension, tension, stretch_residual_norm, stretch_relative_residual, ierr)
  if (ierr /= 0) error stop 'solve_tension_freefree failed in stretch test'

  stretch_rhs_norm = sqrt(sum(rhs**2))
  stretch_tension_norm = sqrt(sum(tension**2))

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  amp = 0.01_mytype
  pi = acos(-1.0_mytype)
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = s
    fibre%x(2, i) = amp * sin(pi * s / fibre%length)
    fibre%x(3, i) = 0._mytype
  end do
  fibre%x_old = fibre%x
  fibre%v = 0._mytype

  left_before = fibre%x(:, 1)
  right_before = fibre%x(:, fibre%nl)

  deallocate(a_non_tension)
  allocate(a_non_tension(3, fibre%nl), fb(3, fibre%nl))
  call compute_fibre_bending_force_freefree(fibre, fb)
  a_non_tension = fb / fibre%rho_tilde

  nseg = fibre%nl - 1
  deallocate(amat, rhs, tension)
  allocate(amat(nseg, nseg), rhs(nseg), tension(nseg))

  call assemble_tension_system_freefree(fibre, dt, a_non_tension, amat, rhs)
  call solve_tension_freefree(fibre, dt, a_non_tension, tension, curved_residual_norm, curved_relative_residual, ierr)
  if (ierr /= 0) error stop 'solve_tension_freefree failed in curved test'

  curved_rhs_norm = sqrt(sum(rhs**2))
  curved_tension_norm = sqrt(sum(tension**2))

  left_after = fibre%x(:, 1)
  right_after = fibre%x(:, fibre%nl)
  curved_endpoint_change_norm = sqrt(sum((left_after - left_before)**2) + sum((right_after - right_before)**2))

  allocate(tau_test(nseg), y1(nseg), y2(nseg))
  do i = 1, nseg
    tau_test(i) = sin(pi * real(i, mytype) / real(nseg + 1, mytype))
  end do

  y1 = matmul(amat, tau_test)
  call apply_tension_operator_freefree(fibre, tau_test, y2)
  tension_matrix_apply_maxdiff = maxval(abs(y1 - y2))

  open(newunit=io_unit, file='stage2_outputs/fibre_tension_solver_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,ES24.16)') 'straight_rhs_norm = ', straight_rhs_norm
  write(io_unit, '(A,ES24.16)') 'straight_tension_norm = ', straight_tension_norm
  write(io_unit, '(A,ES24.16)') 'straight_residual_norm = ', straight_residual_norm
  write(io_unit, '(A,ES24.16)') 'straight_relative_residual = ', straight_relative_residual

  write(io_unit, '(A,ES24.16)') 'translation_rhs_norm = ', translation_rhs_norm
  write(io_unit, '(A,ES24.16)') 'translation_tension_norm = ', translation_tension_norm
  write(io_unit, '(A,ES24.16)') 'translation_residual_norm = ', translation_residual_norm

  write(io_unit, '(A,ES24.16)') 'stretch_rhs_norm = ', stretch_rhs_norm
  write(io_unit, '(A,ES24.16)') 'stretch_tension_norm = ', stretch_tension_norm
  write(io_unit, '(A,ES24.16)') 'stretch_residual_norm = ', stretch_residual_norm
  write(io_unit, '(A,ES24.16)') 'stretch_relative_residual = ', stretch_relative_residual

  write(io_unit, '(A,ES24.16)') 'curved_rhs_norm = ', curved_rhs_norm
  write(io_unit, '(A,ES24.16)') 'curved_tension_norm = ', curved_tension_norm
  write(io_unit, '(A,ES24.16)') 'curved_residual_norm = ', curved_residual_norm
  write(io_unit, '(A,ES24.16)') 'curved_relative_residual = ', curved_relative_residual
  write(io_unit, '(A,ES24.16)') 'curved_endpoint_change_norm = ', curved_endpoint_change_norm

  write(io_unit, '(A,ES24.16)') 'tension_matrix_apply_maxdiff = ', tension_matrix_apply_maxdiff

  close(io_unit)

end program fibre_tension_solver_check
