module fibre_parameters

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t, fibre_allocate

  implicit none
  private

  public :: fibre_init_straight_free_free

contains

  subroutine fibre_init_straight_free_free(fibre, nl, length, rho_tilde, gamma)
    type(fibre_t), intent(inout) :: fibre
    integer, intent(in) :: nl
    real(mytype), intent(in) :: length
    real(mytype), intent(in) :: rho_tilde
    real(mytype), intent(in) :: gamma

    integer :: i

    if (nl < 2) then
      error stop 'fibre_init_straight_free_free: nl must be >= 2'
    end if

    call fibre_allocate(fibre, nl)

    fibre%length = length
    fibre%rho_tilde = rho_tilde
    fibre%gamma = gamma
    fibre%ds = length / real(nl - 1, mytype)

    do i = 1, nl
      fibre%x(1, i) = real(i - 1, mytype) * fibre%ds
      fibre%x(2, i) = 0._mytype
      fibre%x(3, i) = 0._mytype
    end do

    fibre%x_old = fibre%x
    fibre%v = 0._mytype
    fibre%tension = 0._mytype
    fibre%f_ext = 0._mytype
  end subroutine fibre_init_straight_free_free

end module fibre_parameters
