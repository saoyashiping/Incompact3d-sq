module fibre_external_force

  use decomp_2d_constants, only : mytype
  use fibre_types, only : fibre_t

  implicit none
  private

  public :: clear_fibre_external_force
  public :: set_fibre_external_force
  public :: add_fibre_external_force

contains

  subroutine clear_fibre_external_force(fibre)
    type(fibre_t), intent(inout) :: fibre

    if (.not. allocated(fibre%f_ext)) then
      error stop 'clear_fibre_external_force: fibre%f_ext is not allocated'
    end if

    fibre%f_ext = 0._mytype
  end subroutine clear_fibre_external_force

  subroutine set_fibre_external_force(fibre, f_ext_in)
    type(fibre_t), intent(inout) :: fibre
    real(mytype), intent(in) :: f_ext_in(3, fibre%nl)

    if (.not. allocated(fibre%f_ext)) then
      error stop 'set_fibre_external_force: fibre%f_ext is not allocated'
    end if

    fibre%f_ext = f_ext_in
  end subroutine set_fibre_external_force

  subroutine add_fibre_external_force(fibre, f_add)
    type(fibre_t), intent(inout) :: fibre
    real(mytype), intent(in) :: f_add(3, fibre%nl)

    if (.not. allocated(fibre%f_ext)) then
      error stop 'add_fibre_external_force: fibre%f_ext is not allocated'
    end if

    fibre%f_ext = fibre%f_ext + f_add
  end subroutine add_fibre_external_force

end module fibre_external_force
