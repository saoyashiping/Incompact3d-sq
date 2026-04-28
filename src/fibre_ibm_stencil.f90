module fibre_ibm_stencil

  use fibre_parameters, only : mytype
  use fibre_ibm_types, only : ibm_stencil_t

  implicit none
  private

  public :: init_empty_ibm_stencil
  public :: destroy_ibm_stencil

contains

  subroutine init_empty_ibm_stencil(stencil, nl, max_points_per_lag)
    type(ibm_stencil_t), intent(inout) :: stencil
    integer, intent(in) :: nl, max_points_per_lag

    if (nl <= 0) error stop 'init_empty_ibm_stencil: nl must be > 0'
    if (max_points_per_lag < 0) error stop 'init_empty_ibm_stencil: max_points_per_lag must be >= 0'

    call destroy_ibm_stencil(stencil)

    stencil%nl = nl
    stencil%max_points_per_lag = max_points_per_lag

    allocate(stencil%count(nl))
    stencil%count = 0

    if (max_points_per_lag > 0) then
      allocate(stencil%i_idx(max_points_per_lag, nl))
      allocate(stencil%j_idx(max_points_per_lag, nl))
      allocate(stencil%k_idx(max_points_per_lag, nl))
      allocate(stencil%weight(max_points_per_lag, nl))

      stencil%i_idx = 0
      stencil%j_idx = 0
      stencil%k_idx = 0
      stencil%weight = 0._mytype
    end if
  end subroutine init_empty_ibm_stencil

  subroutine destroy_ibm_stencil(stencil)
    type(ibm_stencil_t), intent(inout) :: stencil

    if (allocated(stencil%count)) deallocate(stencil%count)
    if (allocated(stencil%i_idx)) deallocate(stencil%i_idx)
    if (allocated(stencil%j_idx)) deallocate(stencil%j_idx)
    if (allocated(stencil%k_idx)) deallocate(stencil%k_idx)
    if (allocated(stencil%weight)) deallocate(stencil%weight)

    stencil%nl = 0
    stencil%max_points_per_lag = 0
  end subroutine destroy_ibm_stencil

end module fibre_ibm_stencil
