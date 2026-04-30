program fibre_implicit_bending_check

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_parameters, only : fibre_init_straight_free_free
  use fibre_bending_operator, only : compute_fibre_d4_freefree
  use fibre_boundary_freefree, only : fibre_boundary_residual_t
  use fibre_boundary_freefree, only : fibre_compute_freefree_boundary_residual
  use fibre_diagnostics, only : compute_bending_energy
  use fibre_diagnostics, only : compute_curvature_diagnostics
  use fibre_diagnostics, only : compute_total_length_relative_error
  use fibre_implicit_bending, only : build_freefree_d4_matrix
  use fibre_implicit_bending, only : implicit_bending_relax_step

  implicit none

  type(fibre_t) :: fibre
  type(fibre_boundary_residual_t) :: boundary_res

  integer :: i, step, io_unit, ierr, nsteps
  integer :: sine_energy_increase_count, large_dt_energy_increase_count
  integer :: large_dt_solver_failure_count, endpoint_fixed_constraint_detected

  real(mytype) :: pi, s, amp, dt
  real(mytype) :: straight_preservation_error_max
  real(mytype) :: straight_final_bending_energy, straight_final_max_curvature
  real(mytype) :: sine_initial_bending_energy, sine_final_bending_energy
  real(mytype) :: sine_max_energy_increase
  real(mytype) :: sine_initial_max_curvature, sine_final_max_curvature
  real(mytype) :: sine_final_total_length_relative_error
  real(mytype) :: sine_final_freefree_boundary_residual
  real(mytype) :: large_dt_initial_bending_energy, large_dt_final_bending_energy
  real(mytype) :: d4_matrix_operator_maxdiff

  real(mytype), allocatable :: d4mat(:,:), d4op(:,:), d4mat_vec(:)
  real(mytype), allocatable :: x_initial(:,:)
  real(mytype) :: energy_before, energy_after, energy_diff
  real(mytype) :: max_curv_dummy, rms_curv_dummy
  real(mytype) :: left_before(3), right_before(3), left_after(3), right_after(3)

  pi = acos(-1.0_mytype)

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = s
    fibre%x(2, i) = 0.01_mytype * sin(2._mytype * pi * s / fibre%length)
    fibre%x(3, i) = 0.003_mytype * sin(3._mytype * pi * s / fibre%length)
  end do

  allocate(d4mat(fibre%nl, fibre%nl), d4op(3, fibre%nl), d4mat_vec(fibre%nl))
  call build_freefree_d4_matrix(fibre, d4mat)
  call compute_fibre_d4_freefree(fibre, d4op)

  d4_matrix_operator_maxdiff = 0._mytype
  do i = 1, 3
    d4mat_vec = matmul(d4mat, fibre%x(i, :))
    d4_matrix_operator_maxdiff = max(d4_matrix_operator_maxdiff, maxval(abs(d4mat_vec - d4op(i, :))))
  end do
  deallocate(d4mat, d4op, d4mat_vec)

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  allocate(x_initial(3, fibre%nl))
  x_initial = fibre%x

  dt = 1.0e-3_mytype
  nsteps = 10
  do step = 1, nsteps
    call implicit_bending_relax_step(fibre, dt, ierr)
    if (ierr /= 0) error stop 'implicit_bending_relax_step failed in straight test'
  end do

  straight_preservation_error_max = 0._mytype
  do i = 1, fibre%nl
    straight_preservation_error_max = max(straight_preservation_error_max, sqrt(sum((fibre%x(:, i) - x_initial(:, i))**2)))
  end do

  call compute_bending_energy(fibre, straight_final_bending_energy)
  call compute_curvature_diagnostics(fibre, straight_final_max_curvature, rms_curv_dummy)
  deallocate(x_initial)

  call fibre_init_straight_free_free(fibre, nl=65, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  amp = 0.01_mytype
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = s
    fibre%x(2, i) = amp * sin(pi * s / fibre%length)
    fibre%x(3, i) = 0._mytype
  end do

  call compute_bending_energy(fibre, sine_initial_bending_energy)
  call compute_curvature_diagnostics(fibre, sine_initial_max_curvature, rms_curv_dummy)

  dt = 1.0e-4_mytype
  nsteps = 50
  sine_energy_increase_count = 0
  sine_max_energy_increase = 0._mytype

  do step = 1, nsteps
    call compute_bending_energy(fibre, energy_before)
    call implicit_bending_relax_step(fibre, dt, ierr)
    if (ierr /= 0) error stop 'implicit_bending_relax_step failed in sine small-dt test'
    call compute_bending_energy(fibre, energy_after)

    energy_diff = energy_after - energy_before
    if (energy_diff > 0._mytype) then
      sine_energy_increase_count = sine_energy_increase_count + 1
      sine_max_energy_increase = max(sine_max_energy_increase, energy_diff)
    end if
  end do

  call compute_bending_energy(fibre, sine_final_bending_energy)
  call compute_curvature_diagnostics(fibre, sine_final_max_curvature, rms_curv_dummy)
  call compute_total_length_relative_error(fibre, sine_final_total_length_relative_error)
  call fibre_compute_freefree_boundary_residual(fibre, boundary_res)
  sine_final_freefree_boundary_residual = boundary_res%max_freefree_boundary_residual

  call fibre_init_straight_free_free(fibre, nl=65, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  amp = 0.01_mytype
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = s
    fibre%x(2, i) = amp * sin(pi * s / fibre%length)
    fibre%x(3, i) = 0._mytype
  end do

  call compute_bending_energy(fibre, large_dt_initial_bending_energy)

  dt = 1.0e-2_mytype
  nsteps = 10
  large_dt_energy_increase_count = 0
  large_dt_solver_failure_count = 0

  do step = 1, nsteps
    call compute_bending_energy(fibre, energy_before)
    call implicit_bending_relax_step(fibre, dt, ierr)
    if (ierr /= 0) then
      large_dt_solver_failure_count = large_dt_solver_failure_count + 1
      cycle
    end if

    call compute_bending_energy(fibre, energy_after)
    if (energy_after - energy_before > 0._mytype) then
      large_dt_energy_increase_count = large_dt_energy_increase_count + 1
    end if
  end do

  call compute_bending_energy(fibre, large_dt_final_bending_energy)

  call fibre_init_straight_free_free(fibre, nl=65, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = s
    fibre%x(2, i) = amp * sin(pi * s / fibre%length)
    fibre%x(3, i) = 0._mytype
  end do
  fibre%x(2, 1) = fibre%x(2, 1) + 1.0e-3_mytype
  fibre%x(2, fibre%nl) = fibre%x(2, fibre%nl) - 5.0e-4_mytype

  left_before = fibre%x(:, 1)
  right_before = fibre%x(:, fibre%nl)
  call implicit_bending_relax_step(fibre, 1.0e-4_mytype, ierr)
  if (ierr /= 0) error stop 'implicit_bending_relax_step failed in endpoint test'
  left_after = fibre%x(:, 1)
  right_after = fibre%x(:, fibre%nl)

  endpoint_fixed_constraint_detected = 0
  if (sqrt(sum((left_after - left_before)**2)) <= 1.0e-15_mytype .and. &
      sqrt(sum((right_after - right_before)**2)) <= 1.0e-15_mytype) then
    endpoint_fixed_constraint_detected = 1
  end if

  open(newunit=io_unit, file='stage2_outputs/fibre_implicit_bending_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,ES24.16)') 'd4_matrix_operator_maxdiff = ', d4_matrix_operator_maxdiff

  write(io_unit, '(A,ES24.16)') 'straight_preservation_error_max = ', straight_preservation_error_max
  write(io_unit, '(A,ES24.16)') 'straight_final_bending_energy = ', straight_final_bending_energy
  write(io_unit, '(A,ES24.16)') 'straight_final_max_curvature = ', straight_final_max_curvature

  write(io_unit, '(A,ES24.16)') 'sine_initial_bending_energy = ', sine_initial_bending_energy
  write(io_unit, '(A,ES24.16)') 'sine_final_bending_energy = ', sine_final_bending_energy
  write(io_unit, '(A,I0)') 'sine_energy_increase_count = ', sine_energy_increase_count
  write(io_unit, '(A,ES24.16)') 'sine_max_energy_increase = ', sine_max_energy_increase
  write(io_unit, '(A,ES24.16)') 'sine_initial_max_curvature = ', sine_initial_max_curvature
  write(io_unit, '(A,ES24.16)') 'sine_final_max_curvature = ', sine_final_max_curvature
  write(io_unit, '(A,ES24.16)') 'sine_final_total_length_relative_error = ', sine_final_total_length_relative_error
  write(io_unit, '(A,ES24.16)') 'sine_final_freefree_boundary_residual = ', sine_final_freefree_boundary_residual

  write(io_unit, '(A,ES24.16)') 'large_dt_initial_bending_energy = ', large_dt_initial_bending_energy
  write(io_unit, '(A,ES24.16)') 'large_dt_final_bending_energy = ', large_dt_final_bending_energy
  write(io_unit, '(A,I0)') 'large_dt_energy_increase_count = ', large_dt_energy_increase_count
  write(io_unit, '(A,I0)') 'large_dt_solver_failure_count = ', large_dt_solver_failure_count

  write(io_unit, '(A,I0)') 'endpoint_fixed_constraint_detected = ', endpoint_fixed_constraint_detected

  close(io_unit)

end program fibre_implicit_bending_check
