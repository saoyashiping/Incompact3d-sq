module fibre_prod_main_diagnostics
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_rhs_signature_type
  public :: fibre_prod_main_diagnostics_type
  public :: fibre_prod_main_signature
  public :: fibre_prod_main_diagnostics_reset
  public :: fibre_prod_main_diagnostics_record
  public :: fibre_prod_main_diagnostics_finite
  public :: fibre_prod_main_diagnostics_write

  type :: fibre_prod_rhs_signature_type
    real(dp) :: sum_value = 0.0_dp
    real(dp) :: max_abs = 0.0_dp
    logical :: finite = .true.
  end type fibre_prod_rhs_signature_type

  type :: fibre_prod_main_diagnostics_type
    integer :: hook_calls = 0
    integer :: modified_cells = 0
    real(dp) :: lambda_fsi = 0.0_dp
    logical :: enabled = .false.
    logical :: no_contamination = .true.
    logical :: small_lambda_response = .false.
    type(fibre_prod_rhs_signature_type) :: before
    type(fibre_prod_rhs_signature_type) :: after
  end type fibre_prod_main_diagnostics_type

contains

  function fibre_prod_main_signature(rhs_x, rhs_y, rhs_z) result(signature)
    real(dp), intent(in) :: rhs_x(:, :, :)
    real(dp), intent(in) :: rhs_y(:, :, :)
    real(dp), intent(in) :: rhs_z(:, :, :)
    type(fibre_prod_rhs_signature_type) :: signature

    signature%finite = all(ieee_is_finite(rhs_x)) .and. all(ieee_is_finite(rhs_y)) .and. &
                       all(ieee_is_finite(rhs_z))
    if (signature%finite) then
      signature%sum_value = sum(rhs_x) + sum(rhs_y) + sum(rhs_z)
      signature%max_abs = max(maxval(abs(rhs_x)), max(maxval(abs(rhs_y)), maxval(abs(rhs_z))))
    else
      signature%sum_value = huge(1.0_dp)
      signature%max_abs = huge(1.0_dp)
    end if
  end function fibre_prod_main_signature

  subroutine fibre_prod_main_diagnostics_reset(diag)
    type(fibre_prod_main_diagnostics_type), intent(out) :: diag
    diag = fibre_prod_main_diagnostics_type()
  end subroutine fibre_prod_main_diagnostics_reset

  subroutine fibre_prod_main_diagnostics_record(diag, enabled, lambda_fsi, before, after, modified_cells)
    type(fibre_prod_main_diagnostics_type), intent(inout) :: diag
    logical, intent(in) :: enabled
    real(dp), intent(in) :: lambda_fsi
    type(fibre_prod_rhs_signature_type), intent(in) :: before
    type(fibre_prod_rhs_signature_type), intent(in) :: after
    integer, intent(in) :: modified_cells
    real(dp) :: delta_sum
    real(dp) :: delta_max

    diag%hook_calls = diag%hook_calls + 1
    diag%enabled = enabled
    diag%lambda_fsi = lambda_fsi
    diag%before = before
    diag%after = after
    diag%modified_cells = diag%modified_cells + modified_cells
    delta_sum = abs(after%sum_value - before%sum_value)
    delta_max = abs(after%max_abs - before%max_abs)
    if (enabled .and. lambda_fsi == 0.0_dp) diag%no_contamination = delta_sum <= 0.0_dp .and. delta_max <= 0.0_dp
    if (enabled .and. lambda_fsi > 0.0_dp) diag%small_lambda_response = modified_cells > 0 .and. &
                                                   after%finite .and. (delta_sum > 0.0_dp .or. delta_max > 0.0_dp)
  end subroutine fibre_prod_main_diagnostics_record

  logical function fibre_prod_main_diagnostics_finite(diag) result(is_finite)
    type(fibre_prod_main_diagnostics_type), intent(in) :: diag

    is_finite = ieee_is_finite(diag%lambda_fsi) .and. diag%before%finite .and. diag%after%finite .and. &
                ieee_is_finite(diag%before%sum_value) .and. ieee_is_finite(diag%after%sum_value) .and. &
                ieee_is_finite(diag%before%max_abs) .and. ieee_is_finite(diag%after%max_abs)
  end function fibre_prod_main_diagnostics_finite

  subroutine fibre_prod_main_diagnostics_write(diag, path, status)
    type(fibre_prod_main_diagnostics_type), intent(in) :: diag
    character(len=*), intent(in) :: path
    integer, intent(out) :: status
    integer :: unit_id

    status = 0
    open(newunit=unit_id, file=trim(path), status='replace', action='write', iostat=status)
    if (status /= 0) return
    write(unit_id,'(A,I0)') 'hook_calls=', diag%hook_calls
    write(unit_id,'(A,L1)') 'enabled=', diag%enabled
    write(unit_id,'(A,ES24.16)') 'lambda_fsi=', diag%lambda_fsi
    write(unit_id,'(A,I0)') 'modified_cells=', diag%modified_cells
    write(unit_id,'(A,L1)') 'before_finite=', diag%before%finite
    write(unit_id,'(A,L1)') 'after_finite=', diag%after%finite
    write(unit_id,'(A,L1)') 'no_contamination=', diag%no_contamination
    write(unit_id,'(A,L1)') 'small_lambda_response=', diag%small_lambda_response
    close(unit_id)
  end subroutine fibre_prod_main_diagnostics_write

end module fibre_prod_main_diagnostics
