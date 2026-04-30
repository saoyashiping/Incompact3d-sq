program fibre_freefree_boundary_check

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t
  use fibre_parameters, only : fibre_init_straight_free_free
  use fibre_geometry, only : fibre_center_of_mass
  use fibre_boundary_freefree, only : fibre_boundary_residual_t
  use fibre_boundary_freefree, only : fibre_compute_freefree_ghost_points
  use fibre_boundary_freefree, only : fibre_compute_freefree_boundary_residual

  implicit none

  type(fibre_t) :: fibre
  type(fibre_boundary_residual_t) :: residual

  integer :: i, io_unit
  real(mytype) :: s, pi
  real(mytype) :: x_left_before(3), x_right_before(3)
  real(mytype) :: x_left_after(3), x_right_after(3)
  real(mytype) :: x_cm(3)
  real(mytype) :: left_endpoint_change_norm, right_endpoint_change_norm
  real(mytype) :: x_g_l1(3), x_g_l2(3), x_g_r1(3), x_g_r2(3)

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)

  pi = acos(-1.0_mytype)

  do i = 1, fibre%nl
    s = real(i - 1, mytype) * fibre%ds
    fibre%x(1, i) = s
    fibre%x(2, i) = 0.01_mytype * sin(pi * s / fibre%length)
    fibre%x(3, i) = 0.005_mytype * sin(2._mytype * pi * s / fibre%length)
  end do

  x_left_before = fibre%x(:, 1)
  x_right_before = fibre%x(:, fibre%nl)

  call fibre_compute_freefree_ghost_points(fibre, x_g_l1, x_g_l2, x_g_r1, x_g_r2)
  call fibre_compute_freefree_boundary_residual(fibre, residual)

  x_left_after = fibre%x(:, 1)
  x_right_after = fibre%x(:, fibre%nl)

  left_endpoint_change_norm = sqrt(sum((x_left_after - x_left_before)**2))
  right_endpoint_change_norm = sqrt(sum((x_right_after - x_right_before)**2))

  x_cm = fibre_center_of_mass(fibre)

  open(newunit=io_unit, file='stage2_outputs/fibre_freefree_boundary_check.dat', status='replace', action='write', form='formatted')
  write(io_unit, '(A,I0)') 'nl = ', fibre%nl
  write(io_unit, '(A,ES24.16)') 'length = ', fibre%length
  write(io_unit, '(A,ES24.16)') 'ds = ', fibre%ds
  write(io_unit, '(A,ES24.16)') 'rho_tilde = ', fibre%rho_tilde
  write(io_unit, '(A,ES24.16)') 'gamma = ', fibre%gamma
  write(io_unit, '(A,ES24.16)') 'center_of_mass_x = ', x_cm(1)
  write(io_unit, '(A,ES24.16)') 'center_of_mass_y = ', x_cm(2)
  write(io_unit, '(A,ES24.16)') 'center_of_mass_z = ', x_cm(3)
  write(io_unit, '(A,ES24.16)') 'left_d2_norm = ', residual%left_d2_norm
  write(io_unit, '(A,ES24.16)') 'left_d3_norm = ', residual%left_d3_norm
  write(io_unit, '(A,ES24.16)') 'right_d2_norm = ', residual%right_d2_norm
  write(io_unit, '(A,ES24.16)') 'right_d3_norm = ', residual%right_d3_norm
  write(io_unit, '(A,ES24.16)') 'max_freefree_boundary_residual = ', residual%max_freefree_boundary_residual
  write(io_unit, '(A,ES24.16)') 'left_endpoint_change_norm = ', left_endpoint_change_norm
  write(io_unit, '(A,ES24.16)') 'right_endpoint_change_norm = ', right_endpoint_change_norm
  write(io_unit, '(A,I0)') 'tension_size = ', size(fibre%tension)
  close(io_unit)

end program fibre_freefree_boundary_check
