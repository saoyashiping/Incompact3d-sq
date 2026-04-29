program fibre_ibm_structure_coupling_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use fibre_parameters, only : mytype
  use fibre_types, only : fibre_t, fibre_allocate
  use fibre_geometry, only : fibre_center_of_mass, fibre_end_to_end_distance
  use fibre_external_force, only : clear_fibre_external_force
  use fibre_structure_integrator, only : advance_fibre_structure_freefree
  use fibre_ibm_types, only : ibm_grid_t
  use fibre_ibm_grid, only : init_uniform_ibm_grid, destroy_ibm_grid
  use fibre_ibm_structure_coupling, only : apply_ibm_structure_force_to_fibre
  use fibre_ibm_structure_coupling, only : compute_and_apply_ibm_structure_force

  implicit none

  type(ibm_grid_t) :: grid
  type(fibre_t) :: fibre

  integer, parameter :: nl = 33
  integer, parameter :: nsteps = 20
  integer :: i, j, k, l, step, io_unit, ierr
  integer :: uniform_drag_direction_check, reverse_drag_direction_check
  integer :: uniform_drag_solver_failure_count
  integer :: ibm_structure_coupling_status

  real(mytype), parameter :: fibre_length = 1.0_mytype
  real(mytype), parameter :: rho_tilde = 1.0_mytype
  real(mytype), parameter :: gamma = 1.0_mytype
  real(mytype), parameter :: beta_drag = 10.0_mytype
  real(mytype), parameter :: dt = 1.0e-5_mytype

  real(mytype), allocatable :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
  real(mytype), allocatable :: u_lag(:,:), force_on_structure(:,:), force_on_fluid(:,:)
  real(mytype), allocatable :: x_initial(:,:), ftmp(:,:)

  real(mytype) :: x, y, z
  real(mytype) :: t_final, expected_vx, tension_residual, relative_tension_residual
  real(mytype) :: cm(3)

  real(mytype) :: zero_slip_f_ext_norm, zero_slip_force_structure_norm, zero_slip_force_fluid_norm
  real(mytype) :: uniform_drag_final_center_velocity_x, uniform_drag_expected_center_velocity_x
  real(mytype) :: uniform_drag_center_velocity_error, uniform_drag_shape_error_max, uniform_drag_length_error
  real(mytype) :: reverse_drag_final_center_velocity_x
  real(mytype) :: force_set_error, force_add_error, force_clear_error
  real(mytype) :: nonuniform_force_variation_norm, nonuniform_f_ext_norm, nonuniform_force_matches_f_ext_error

  call init_uniform_ibm_grid(grid, nx=16, ny=12, nz=10, &
                             xmin=0._mytype, xmax=2._mytype, &
                             ymin=-1._mytype, ymax=1._mytype, &
                             zmin=0._mytype, zmax=1._mytype, &
                             periodic_x=.true., periodic_y=.false., periodic_z=.true.)

  call fibre_allocate(fibre, nl)
  fibre%length = fibre_length
  fibre%ds = fibre%length / real(fibre%nl - 1, mytype)
  fibre%rho_tilde = rho_tilde
  fibre%gamma = gamma

  allocate(ux(grid%nx, grid%ny, grid%nz), uy(grid%nx, grid%ny, grid%nz), uz(grid%nx, grid%ny, grid%nz))
  allocate(u_lag(3, fibre%nl), force_on_structure(3, fibre%nl), force_on_fluid(3, fibre%nl))
  allocate(x_initial(3, fibre%nl), ftmp(3, fibre%nl))

  call reset_straight_fibre(fibre)

  ! Test A: zero slip
  ux = 0.2_mytype
  uy = -0.1_mytype
  uz = 0.05_mytype
  do l = 1, fibre%nl
    fibre%v(:, l) = [0.2_mytype, -0.1_mytype, 0.05_mytype]
  end do
  call compute_and_apply_ibm_structure_force(grid, ux, uy, uz, fibre, beta_drag, 'set', u_lag, force_on_structure, force_on_fluid)
  zero_slip_f_ext_norm = sqrt(sum(fibre%f_ext**2))
  zero_slip_force_structure_norm = sqrt(sum(force_on_structure**2))
  zero_slip_force_fluid_norm = sqrt(sum(force_on_fluid**2))

  ! Test B: uniform drag on stationary fibre
  call reset_straight_fibre(fibre)
  x_initial = fibre%x
  call clear_fibre_external_force(fibre)
  ux = 0.2_mytype
  uy = 0._mytype
  uz = 0._mytype
  uniform_drag_solver_failure_count = 0
  do step = 1, nsteps
    call compute_and_apply_ibm_structure_force(grid, ux, uy, uz, fibre, beta_drag, 'set', u_lag, force_on_structure, force_on_fluid)
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
    if (ierr /= 0) uniform_drag_solver_failure_count = uniform_drag_solver_failure_count + 1
  end do

  cm = fibre_center_of_mass(fibre)
  uniform_drag_final_center_velocity_x = sum(fibre%v(1, :)) / real(fibre%nl, mytype)
  t_final = dt * real(nsteps, mytype)
  expected_vx = 0.2_mytype * (1._mytype - exp(-beta_drag * t_final / fibre%rho_tilde))
  uniform_drag_expected_center_velocity_x = expected_vx
  uniform_drag_center_velocity_error = abs(uniform_drag_final_center_velocity_x - uniform_drag_expected_center_velocity_x)
  uniform_drag_direction_check = 0
  if (uniform_drag_final_center_velocity_x > 0._mytype) uniform_drag_direction_check = 1

  uniform_drag_shape_error_max = maxval(abs((fibre%x(2:3, :) - x_initial(2:3, :))))
  uniform_drag_length_error = abs(fibre_end_to_end_distance(fibre) - fibre%length)

  ! Test C: reverse uniform drag direction
  call reset_straight_fibre(fibre)
  ux = -0.2_mytype
  uy = 0._mytype
  uz = 0._mytype
  do step = 1, nsteps
    call compute_and_apply_ibm_structure_force(grid, ux, uy, uz, fibre, beta_drag, 'set', u_lag, force_on_structure, force_on_fluid)
    call advance_fibre_structure_freefree(fibre, dt, tension_residual, relative_tension_residual, ierr)
  end do
  reverse_drag_final_center_velocity_x = sum(fibre%v(1, :)) / real(fibre%nl, mytype)
  reverse_drag_direction_check = 0
  if (reverse_drag_final_center_velocity_x < 0._mytype) reverse_drag_direction_check = 1

  ! Test D: set/add/clear modes
  ftmp(1, :) = 1.0_mytype
  ftmp(2, :) = 2.0_mytype
  ftmp(3, :) = 3.0_mytype
  call apply_ibm_structure_force_to_fibre(fibre, ftmp, 'set')
  force_set_error = maxval(abs(fibre%f_ext - ftmp))

  force_on_structure(1, :) = -0.5_mytype
  force_on_structure(2, :) = 0.25_mytype
  force_on_structure(3, :) = 0.75_mytype
  call apply_ibm_structure_force_to_fibre(fibre, force_on_structure, 'add')
  force_add_error = maxval(abs(fibre%f_ext - (ftmp + force_on_structure)))

  call apply_ibm_structure_force_to_fibre(fibre, force_on_structure, 'clear')
  force_clear_error = maxval(abs(fibre%f_ext))

  ! Test E: nonuniform flow gives nonuniform structure force
  call reset_straight_fibre(fibre)
  ux = 0._mytype
  uy = 0._mytype
  uz = 0._mytype
  do k = 1, grid%nz
    z = grid%z(k)
    do j = 1, grid%ny
      y = grid%y(j)
      do i = 1, grid%nx
        x = grid%x(i)
        ux(i, j, k) = 1.0_mytype + 0.2_mytype * x - 0.1_mytype * y + 0.05_mytype * z
        uy(i, j, k) = -0.3_mytype + 0.4_mytype * y
        uz(i, j, k) = 0.7_mytype - 0.2_mytype * z + 0.1_mytype * x
      end do
    end do
  end do
  call compute_and_apply_ibm_structure_force(grid, ux, uy, uz, fibre, beta_drag, 'set', u_lag, force_on_structure, force_on_fluid)
  nonuniform_f_ext_norm = sqrt(sum(fibre%f_ext**2))
  nonuniform_force_variation_norm = sqrt(sum((force_on_structure(:, 2:fibre%nl) - force_on_structure(:, 1:fibre%nl-1))**2))
  nonuniform_force_matches_f_ext_error = maxval(abs(force_on_structure - fibre%f_ext))

  ibm_structure_coupling_status = 1
  if (ieee_is_nan(zero_slip_f_ext_norm) .or. ieee_is_nan(zero_slip_force_structure_norm) .or. &
      ieee_is_nan(zero_slip_force_fluid_norm) .or. ieee_is_nan(uniform_drag_final_center_velocity_x) .or. &
      ieee_is_nan(uniform_drag_expected_center_velocity_x) .or. ieee_is_nan(uniform_drag_center_velocity_error) .or. &
      ieee_is_nan(uniform_drag_shape_error_max) .or. ieee_is_nan(uniform_drag_length_error) .or. &
      ieee_is_nan(reverse_drag_final_center_velocity_x) .or. ieee_is_nan(force_set_error) .or. &
      ieee_is_nan(force_add_error) .or. ieee_is_nan(force_clear_error) .or. ieee_is_nan(nonuniform_force_variation_norm) .or. &
      ieee_is_nan(nonuniform_f_ext_norm) .or. ieee_is_nan(nonuniform_force_matches_f_ext_error)) then
    ibm_structure_coupling_status = 0
  end if

  open(newunit=io_unit, file='stage3_outputs/fibre_ibm_structure_coupling_check.dat', &
       status='replace', action='write', form='formatted')

  write(io_unit, '(A,ES24.16)') 'zero_slip_f_ext_norm = ', zero_slip_f_ext_norm
  write(io_unit, '(A,ES24.16)') 'zero_slip_force_structure_norm = ', zero_slip_force_structure_norm
  write(io_unit, '(A,ES24.16)') 'zero_slip_force_fluid_norm = ', zero_slip_force_fluid_norm

  write(io_unit, '(A,ES24.16)') 'uniform_drag_final_center_velocity_x = ', uniform_drag_final_center_velocity_x
  write(io_unit, '(A,ES24.16)') 'uniform_drag_expected_center_velocity_x = ', uniform_drag_expected_center_velocity_x
  write(io_unit, '(A,ES24.16)') 'uniform_drag_center_velocity_error = ', uniform_drag_center_velocity_error
  write(io_unit, '(A,I0)') 'uniform_drag_direction_check = ', uniform_drag_direction_check
  write(io_unit, '(A,ES24.16)') 'uniform_drag_shape_error_max = ', uniform_drag_shape_error_max
  write(io_unit, '(A,ES24.16)') 'uniform_drag_length_error = ', uniform_drag_length_error
  write(io_unit, '(A,I0)') 'uniform_drag_solver_failure_count = ', uniform_drag_solver_failure_count

  write(io_unit, '(A,ES24.16)') 'reverse_drag_final_center_velocity_x = ', reverse_drag_final_center_velocity_x
  write(io_unit, '(A,I0)') 'reverse_drag_direction_check = ', reverse_drag_direction_check

  write(io_unit, '(A,ES24.16)') 'force_set_error = ', force_set_error
  write(io_unit, '(A,ES24.16)') 'force_add_error = ', force_add_error
  write(io_unit, '(A,ES24.16)') 'force_clear_error = ', force_clear_error

  write(io_unit, '(A,ES24.16)') 'nonuniform_force_variation_norm = ', nonuniform_force_variation_norm
  write(io_unit, '(A,ES24.16)') 'nonuniform_f_ext_norm = ', nonuniform_f_ext_norm
  write(io_unit, '(A,ES24.16)') 'nonuniform_force_matches_f_ext_error = ', nonuniform_force_matches_f_ext_error

  write(io_unit, '(A,I0)') 'ibm_structure_coupling_status = ', ibm_structure_coupling_status

  close(io_unit)

  deallocate(ux, uy, uz, u_lag, force_on_structure, force_on_fluid, x_initial, ftmp)
  call destroy_ibm_grid(grid)

contains

  subroutine reset_straight_fibre(f)
    type(fibre_t), intent(inout) :: f

    integer :: p
    real(mytype) :: s

    do p = 1, f%nl
      s = real(p - 1, mytype) / real(f%nl - 1, mytype)
      f%x(1, p) = 0.5_mytype + s
      f%x(2, p) = 0._mytype
      f%x(3, p) = 0.5_mytype
    end do
    f%x_old = f%x
    f%v = 0._mytype
    f%tension = 0._mytype
    f%f_ext = 0._mytype
  end subroutine reset_straight_fibre

end program fibre_ibm_structure_coupling_check
