module fibre_ibm_delta

  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_grid_t

  implicit none
  private

  public :: peskin_phi_4pt
  public :: delta_weight_1d
  public :: periodic_distance_1d
  public :: delta_weight_1d_periodic
  public :: delta_weight_3d

contains

  pure function peskin_phi_4pt(r) result(phi)
    real(mytype), intent(in) :: r
    real(mytype) :: phi
    real(mytype) :: a, arg

    a = abs(r)

    if (a < 1._mytype) then
      arg = max(1._mytype + 4._mytype * a - 4._mytype * r * r, 0._mytype)
      phi = 0.125_mytype * (3._mytype - 2._mytype * a + sqrt(arg))
    else if (a < 2._mytype) then
      arg = max(-7._mytype + 12._mytype * a - 4._mytype * r * r, 0._mytype)
      phi = 0.125_mytype * (5._mytype - 2._mytype * a - sqrt(arg))
    else
      phi = 0._mytype
    end if
  end function peskin_phi_4pt

  pure function delta_weight_1d(coord, xlag, dx) result(weight)
    real(mytype), intent(in) :: coord, xlag, dx
    real(mytype) :: weight

    if (dx <= 0._mytype) error stop 'delta_weight_1d: dx must be > 0'

    weight = peskin_phi_4pt((coord - xlag) / dx)
  end function delta_weight_1d

  pure function periodic_distance_1d(coord, xlag, xmin, xmax, periodic) result(d)
    real(mytype), intent(in) :: coord, xlag, xmin, xmax
    logical, intent(in) :: periodic
    real(mytype) :: d, l

    d = coord - xlag

    if (periodic) then
      l = xmax - xmin
      if (l <= 0._mytype) error stop 'periodic_distance_1d: xmax must be > xmin when periodic'
      if (d > 0.5_mytype * l) d = d - l
      if (d < -0.5_mytype * l) d = d + l
    end if
  end function periodic_distance_1d

  pure function delta_weight_1d_periodic(coord, xlag, dx, xmin, xmax, periodic) result(weight)
    real(mytype), intent(in) :: coord, xlag, dx, xmin, xmax
    logical, intent(in) :: periodic
    real(mytype) :: weight
    real(mytype) :: d

    if (dx <= 0._mytype) error stop 'delta_weight_1d_periodic: dx must be > 0'

    d = periodic_distance_1d(coord, xlag, xmin, xmax, periodic)
    weight = peskin_phi_4pt(d / dx)
  end function delta_weight_1d_periodic

  pure function delta_weight_3d(grid, xlag, i, j, k) result(weight)
    type(ibm_grid_t), intent(in) :: grid
    real(mytype), intent(in) :: xlag(3)
    integer, intent(in) :: i, j, k
    real(mytype) :: weight
    real(mytype) :: wx, wy, wz

    if (i < 1 .or. i > grid%nx) error stop 'delta_weight_3d: i out of range'
    if (j < 1 .or. j > grid%ny) error stop 'delta_weight_3d: j out of range'
    if (k < 1 .or. k > grid%nz) error stop 'delta_weight_3d: k out of range'

    wx = delta_weight_1d_periodic(grid%x(i), xlag(1), grid%dx, grid%xmin, grid%xmax, grid%periodic_x)
    wy = delta_weight_1d_periodic(grid%y(j), xlag(2), grid%dy, grid%ymin, grid%ymax, grid%periodic_y)
    wz = delta_weight_1d_periodic(grid%z(k), xlag(3), grid%dz, grid%zmin, grid%zmax, grid%periodic_z)

    weight = wx * wy * wz
  end function delta_weight_3d

end module fibre_ibm_delta
