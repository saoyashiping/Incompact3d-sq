module fibre_stage5_coupled_order
  use fibre_parameters, only : mytype
  implicit none
contains
  subroutine perform_stage5_synthetic_two_way_ordered_step(dummy)
    integer, intent(out) :: dummy
    dummy = 1
  end subroutine
end module
