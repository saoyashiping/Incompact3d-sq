module fibre_ibm_types

  use fibre_parameters, only : mytype

  implicit none
  private

  type, public :: ibm_grid_t
    integer :: nx = 0, ny = 0, nz = 0
    real(mytype) :: xmin = 0._mytype, xmax = 0._mytype
    real(mytype) :: ymin = 0._mytype, ymax = 0._mytype
    real(mytype) :: zmin = 0._mytype, zmax = 0._mytype
    real(mytype) :: dx = 0._mytype, dy = 0._mytype, dz = 0._mytype
    real(mytype) :: cell_volume = 0._mytype
    logical :: periodic_x = .false.
    logical :: periodic_y = .false.
    logical :: periodic_z = .false.
    real(mytype), allocatable :: x(:)
    real(mytype), allocatable :: y(:)
    real(mytype), allocatable :: z(:)
  end type ibm_grid_t

  type, public :: ibm_lagrangian_points_t
    integer :: nl = 0
    real(mytype), allocatable :: x(:,:)
    real(mytype), allocatable :: v(:,:)
    real(mytype), allocatable :: force(:,:)
    real(mytype), allocatable :: weight(:)
  end type ibm_lagrangian_points_t

  type, public :: ibm_stencil_t
    integer :: nl = 0
    integer :: max_points_per_lag = 0
    integer, allocatable :: count(:)
    integer, allocatable :: i_idx(:,:)
    integer, allocatable :: j_idx(:,:)
    integer, allocatable :: k_idx(:,:)
    real(mytype), allocatable :: weight(:,:)
  end type ibm_stencil_t

end module fibre_ibm_types
