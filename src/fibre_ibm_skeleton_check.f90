program fibre_ibm_skeleton_check

  use fibre_parameters, only : mytype
  use fibre_types, only : fibre_t
  use fibre_parameters, only : fibre_init_straight_free_free
  use fibre_ibm_types, only : ibm_grid_t, ibm_lagrangian_points_t, ibm_stencil_t
  use fibre_ibm_grid, only : init_uniform_ibm_grid, destroy_ibm_grid
  use fibre_ibm_grid, only : init_lagrangian_points_from_fibre, destroy_lagrangian_points
  use fibre_ibm_stencil, only : init_empty_ibm_stencil, destroy_ibm_stencil
  use fibre_ibm_diagnostics, only : compute_lagrangian_weight_sum
  use fibre_ibm_diagnostics, only : compute_grid_basic_diagnostics, check_lagrangian_points_inside_box
  use fibre_ibm_diagnostics, only : compute_stencil_total_count

  implicit none

  type(ibm_grid_t) :: grid
  type(ibm_lagrangian_points_t) :: lag
  type(ibm_stencil_t) :: stencil
  type(fibre_t) :: fibre

  integer :: i, io_unit
  integer :: grid_nx, grid_ny, grid_nz
  integer :: grid_periodic_x, grid_periodic_y, grid_periodic_z
  integer :: lag_inside_count, lag_outside_count
  integer :: stencil_total_count, ibm_skeleton_status
  real(mytype) :: grid_dx, grid_dy, grid_dz, grid_cell_volume
  real(mytype) :: lag_weight_sum, expected_fibre_length

  call init_uniform_ibm_grid(grid, nx=16, ny=12, nz=10, &
                             xmin=0._mytype, xmax=2._mytype, &
                             ymin=-1._mytype, ymax=1._mytype, &
                             zmin=0._mytype, zmax=1._mytype, &
                             periodic_x=.true., periodic_y=.false., periodic_z=.true.)

  call fibre_init_straight_free_free(fibre, nl=33, length=1.0_mytype, rho_tilde=1.0_mytype, gamma=1.0_mytype)
  do i = 1, fibre%nl
    fibre%x(1, i) = 0.5_mytype + real(i - 1, mytype) * fibre%ds
    fibre%x(2, i) = 0._mytype
    fibre%x(3, i) = 0.5_mytype
  end do
  fibre%x_old = fibre%x
  fibre%v = 0._mytype

  call init_lagrangian_points_from_fibre(lag, fibre)
  call init_empty_ibm_stencil(stencil, nl=lag%nl, max_points_per_lag=0)

  call compute_grid_basic_diagnostics(grid, grid_nx, grid_ny, grid_nz, grid_dx, grid_dy, grid_dz, grid_cell_volume, &
                                      grid_periodic_x, grid_periodic_y, grid_periodic_z)
  call compute_lagrangian_weight_sum(lag, lag_weight_sum)
  call check_lagrangian_points_inside_box(grid, lag, lag_inside_count, lag_outside_count)

  call compute_stencil_total_count(stencil, stencil_total_count)
  expected_fibre_length = fibre%length
  ibm_skeleton_status = 1

  open(newunit=io_unit, file='stage3_outputs/fibre_ibm_skeleton_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,I0)') 'grid_nx = ', grid_nx
  write(io_unit, '(A,I0)') 'grid_ny = ', grid_ny
  write(io_unit, '(A,I0)') 'grid_nz = ', grid_nz
  write(io_unit, '(A,ES24.16)') 'grid_dx = ', grid_dx
  write(io_unit, '(A,ES24.16)') 'grid_dy = ', grid_dy
  write(io_unit, '(A,ES24.16)') 'grid_dz = ', grid_dz
  write(io_unit, '(A,ES24.16)') 'grid_cell_volume = ', grid_cell_volume
  write(io_unit, '(A,I0)') 'grid_periodic_x = ', grid_periodic_x
  write(io_unit, '(A,I0)') 'grid_periodic_y = ', grid_periodic_y
  write(io_unit, '(A,I0)') 'grid_periodic_z = ', grid_periodic_z

  write(io_unit, '(A,I0)') 'lag_nl = ', lag%nl
  write(io_unit, '(A,ES24.16)') 'lag_weight_sum = ', lag_weight_sum
  write(io_unit, '(A,ES24.16)') 'expected_fibre_length = ', expected_fibre_length
  write(io_unit, '(A,I0)') 'lag_inside_count = ', lag_inside_count
  write(io_unit, '(A,I0)') 'lag_outside_count = ', lag_outside_count

  write(io_unit, '(A,I0)') 'stencil_nl = ', stencil%nl
  write(io_unit, '(A,I0)') 'stencil_max_points_per_lag = ', stencil%max_points_per_lag
  write(io_unit, '(A,I0)') 'stencil_total_count = ', stencil_total_count

  write(io_unit, '(A,I0)') 'ibm_skeleton_status = ', ibm_skeleton_status

  close(io_unit)

  call destroy_ibm_stencil(stencil)
  call destroy_lagrangian_points(lag)
  call destroy_ibm_grid(grid)

end program fibre_ibm_skeleton_check
