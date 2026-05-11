module fibre_stage5_rk_substep_sync
  use fibre_parameters, only : mytype
  implicit none
contains
  subroutine perform_stage5_synthetic_rk_substep(flag)
    integer, intent(out) :: flag
    flag = 1
  end subroutine
end module
