program fibre_ibm_delta_check

  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t
  use fibre_ibm_grid, only : init_uniform_ibm_grid, destroy_ibm_grid
  use fibre_ibm_delta, only : peskin_phi_4pt, delta_weight_1d
  use fibre_ibm_delta, only : periodic_distance_1d, delta_weight_3d

  implicit none

  type(ibm_grid_t) :: grid

  integer :: i, j, k, nx1d, io_unit
  real(mytype) :: xmin1d, xmax1d, dx1d, xlag1d, xi
  real(mytype) :: xlag3d(3), q0, w
  real(mytype) :: dpx, dpy, dpz
  integer :: delta_kernel_status

  real(mytype) :: phi_0, phi_05, phi_10, phi_15, phi_20
  real(mytype) :: weight_sum_1d, first_moment_1d
  real(mytype) :: weight_sum_3d, first_moment_x_3d, first_moment_y_3d, first_moment_z_3d
  real(mytype) :: constant_field_value, constant_field_interpolated, constant_field_error
  real(mytype) :: periodic_wrap_distance_x, nonperiodic_distance_x
  real(mytype) :: phi_outside_positive, phi_outside_negative

  phi_0 = peskin_phi_4pt(0._mytype)
  phi_05 = peskin_phi_4pt(0.5_mytype)
  phi_10 = peskin_phi_4pt(1.0_mytype)
  phi_15 = peskin_phi_4pt(1.5_mytype)
  phi_20 = peskin_phi_4pt(2.0_mytype)

  nx1d = 64
  xmin1d = 0._mytype
  xmax1d = 1._mytype
  dx1d = (xmax1d - xmin1d) / real(nx1d, mytype)
  xlag1d = 0.37123_mytype

  weight_sum_1d = 0._mytype
  first_moment_1d = 0._mytype
  do i = 1, nx1d
    xi = xmin1d + (real(i, mytype) - 0.5_mytype) * dx1d
    w = delta_weight_1d(xi, xlag1d, dx1d)
    weight_sum_1d = weight_sum_1d + w
    first_moment_1d = first_moment_1d + (xi - xlag1d) * w
  end do

  call init_uniform_ibm_grid(grid, nx=16, ny=12, nz=10, &
                             xmin=0._mytype, xmax=2._mytype, &
                             ymin=-1._mytype, ymax=1._mytype, &
                             zmin=0._mytype, zmax=1._mytype, &
                             periodic_x=.true., periodic_y=.false., periodic_z=.true.)

  xlag3d = [0.73_mytype, 0.17_mytype, 0.41_mytype]

  weight_sum_3d = 0._mytype
  first_moment_x_3d = 0._mytype
  first_moment_y_3d = 0._mytype
  first_moment_z_3d = 0._mytype

  do k = 1, grid%nz
    do j = 1, grid%ny
      do i = 1, grid%nx
        w = delta_weight_3d(grid, xlag3d, i, j, k)
        dpx = periodic_distance_1d(grid%x(i), xlag3d(1), grid%xmin, grid%xmax, grid%periodic_x)
        dpy = periodic_distance_1d(grid%y(j), xlag3d(2), grid%ymin, grid%ymax, grid%periodic_y)
        dpz = periodic_distance_1d(grid%z(k), xlag3d(3), grid%zmin, grid%zmax, grid%periodic_z)

        weight_sum_3d = weight_sum_3d + w
        first_moment_x_3d = first_moment_x_3d + dpx * w
        first_moment_y_3d = first_moment_y_3d + dpy * w
        first_moment_z_3d = first_moment_z_3d + dpz * w
      end do
    end do
  end do

  q0 = 3.141592653589793_mytype
  constant_field_value = q0
  constant_field_interpolated = q0 * weight_sum_3d
  constant_field_error = abs(constant_field_interpolated - constant_field_value)

  periodic_wrap_distance_x = periodic_distance_1d(0.05_mytype, 1.95_mytype, 0._mytype, 2._mytype, .true.)
  nonperiodic_distance_x = periodic_distance_1d(0.05_mytype, 1.95_mytype, 0._mytype, 2._mytype, .false.)

  phi_outside_positive = peskin_phi_4pt(2.1_mytype)
  phi_outside_negative = peskin_phi_4pt(-2.1_mytype)

  delta_kernel_status = 1
  if (ieee_is_nan(phi_0) .or. ieee_is_nan(phi_05) .or. ieee_is_nan(phi_10) .or. ieee_is_nan(phi_15) .or. &
      ieee_is_nan(phi_20) .or. ieee_is_nan(weight_sum_1d) .or. ieee_is_nan(first_moment_1d) .or. &
      ieee_is_nan(weight_sum_3d) .or. ieee_is_nan(first_moment_x_3d) .or. ieee_is_nan(first_moment_y_3d) .or. &
      ieee_is_nan(first_moment_z_3d) .or. ieee_is_nan(constant_field_interpolated) .or. ieee_is_nan(constant_field_error) .or. &
      phi_0 < 0._mytype .or. phi_05 < 0._mytype .or. phi_10 < 0._mytype .or. phi_15 < 0._mytype .or. &
      phi_20 < 0._mytype .or. phi_outside_positive < 0._mytype .or. phi_outside_negative < 0._mytype) then
    delta_kernel_status = 0
  end if

  open(newunit=io_unit, file='stage3_outputs/fibre_ibm_delta_check.dat', status='replace', action='write', form='formatted')

  write(io_unit, '(A,ES24.16)') 'phi_0 = ', phi_0
  write(io_unit, '(A,ES24.16)') 'phi_05 = ', phi_05
  write(io_unit, '(A,ES24.16)') 'phi_10 = ', phi_10
  write(io_unit, '(A,ES24.16)') 'phi_15 = ', phi_15
  write(io_unit, '(A,ES24.16)') 'phi_20 = ', phi_20

  write(io_unit, '(A,ES24.16)') 'weight_sum_1d = ', weight_sum_1d
  write(io_unit, '(A,ES24.16)') 'first_moment_1d = ', first_moment_1d

  write(io_unit, '(A,ES24.16)') 'weight_sum_3d = ', weight_sum_3d
  write(io_unit, '(A,ES24.16)') 'first_moment_x_3d = ', first_moment_x_3d
  write(io_unit, '(A,ES24.16)') 'first_moment_y_3d = ', first_moment_y_3d
  write(io_unit, '(A,ES24.16)') 'first_moment_z_3d = ', first_moment_z_3d

  write(io_unit, '(A,ES24.16)') 'constant_field_value = ', constant_field_value
  write(io_unit, '(A,ES24.16)') 'constant_field_interpolated = ', constant_field_interpolated
  write(io_unit, '(A,ES24.16)') 'constant_field_error = ', constant_field_error

  write(io_unit, '(A,ES24.16)') 'periodic_wrap_distance_x = ', periodic_wrap_distance_x
  write(io_unit, '(A,ES24.16)') 'nonperiodic_distance_x = ', nonperiodic_distance_x

  write(io_unit, '(A,ES24.16)') 'phi_outside_positive = ', phi_outside_positive
  write(io_unit, '(A,ES24.16)') 'phi_outside_negative = ', phi_outside_negative

  write(io_unit, '(A,I0)') 'delta_kernel_status = ', delta_kernel_status

  close(io_unit)

  call destroy_ibm_grid(grid)

end program fibre_ibm_delta_check
