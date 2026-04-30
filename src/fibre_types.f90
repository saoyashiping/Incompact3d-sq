module fibre_types

  use decomp_2d_constants, only : mytype

  implicit none
  private

  type, public :: fibre_t
    integer :: nl = 0
    real(mytype) :: length = 0._mytype
    real(mytype) :: ds = 0._mytype
    real(mytype) :: rho_tilde = 0._mytype
    real(mytype) :: gamma = 0._mytype
    real(mytype), allocatable :: x(:,:)
    real(mytype), allocatable :: x_old(:,:)
    real(mytype), allocatable :: v(:,:)
    ! Stage 2.1 only stores segment/interface tension placeholders T_{i+1/2}.
    ! The final constrained tension equation and exact boundary tension
    ! treatment will be implemented in Stage 2.5.
    real(mytype), allocatable :: tension(:)
    real(mytype), allocatable :: f_ext(:,:)
  end type fibre_t

  public :: fibre_allocate

contains

  subroutine fibre_allocate(fibre, nl)
    type(fibre_t), intent(inout) :: fibre
    integer, intent(in) :: nl

    fibre%nl = nl

    if (allocated(fibre%x)) deallocate(fibre%x)
    if (allocated(fibre%x_old)) deallocate(fibre%x_old)
    if (allocated(fibre%v)) deallocate(fibre%v)
    if (allocated(fibre%tension)) deallocate(fibre%tension)
    if (allocated(fibre%f_ext)) deallocate(fibre%f_ext)

    allocate(fibre%x(3, nl))
    allocate(fibre%x_old(3, nl))
    allocate(fibre%v(3, nl))
    allocate(fibre%tension(nl - 1))
    allocate(fibre%f_ext(3, nl))

    fibre%f_ext = 0._mytype
  end subroutine fibre_allocate

end module fibre_types
