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
  public :: fibre_prod_main_audit_write

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
    integer :: last_status = 0
    integer :: failed_calls = 0
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

  subroutine fibre_prod_main_diagnostics_record(diag, enabled, lambda_fsi, before, after, modified_cells, status_code)
    type(fibre_prod_main_diagnostics_type), intent(inout) :: diag
    logical, intent(in) :: enabled
    real(dp), intent(in) :: lambda_fsi
    type(fibre_prod_rhs_signature_type), intent(in) :: before
    type(fibre_prod_rhs_signature_type), intent(in) :: after
    integer, intent(in) :: modified_cells
    integer, intent(in), optional :: status_code
    real(dp) :: delta_sum
    real(dp) :: delta_max
    integer :: local_status

    local_status = 0
    if (present(status_code)) local_status = status_code

    diag%hook_calls = diag%hook_calls + 1
    diag%last_status = local_status
    if (local_status /= 0) diag%failed_calls = diag%failed_calls + 1
    diag%enabled = enabled
    diag%lambda_fsi = lambda_fsi
    diag%before = before
    diag%after = after
    diag%modified_cells = diag%modified_cells + modified_cells
    delta_sum = abs(after%sum_value - before%sum_value)
    delta_max = abs(after%max_abs - before%max_abs)

    if (enabled .and. lambda_fsi == 0.0_dp) then
      diag%no_contamination = diag%no_contamination .and. &
                               before%finite .and. after%finite .and. &
                               modified_cells == 0 .and. &
                               delta_sum <= 0.0_dp .and. delta_max <= 0.0_dp
    end if

    if (enabled .and. lambda_fsi > 0.0_dp) then
      diag%small_lambda_response = diag%small_lambda_response .or. &
                                   (modified_cells > 0 .and. after%finite .and. &
                                   (delta_sum > 0.0_dp .or. delta_max > 0.0_dp))
    end if
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
    write(unit_id,'(A,I0)') 'last_status=', diag%last_status
    write(unit_id,'(A,I0)') 'failed_calls=', diag%failed_calls
    write(unit_id,'(A,L1)') 'before_finite=', diag%before%finite
    write(unit_id,'(A,L1)') 'after_finite=', diag%after%finite
    write(unit_id,'(A,ES24.16)') 'before_sum=', diag%before%sum_value
    write(unit_id,'(A,ES24.16)') 'after_sum=', diag%after%sum_value
    write(unit_id,'(A,ES24.16)') 'before_max_abs=', diag%before%max_abs
    write(unit_id,'(A,ES24.16)') 'after_max_abs=', diag%after%max_abs
    write(unit_id,'(A,L1)') 'no_contamination=', diag%no_contamination
    write(unit_id,'(A,L1)') 'small_lambda_response=', diag%small_lambda_response
    close(unit_id)
  end subroutine fibre_prod_main_diagnostics_write

  subroutine fibre_prod_main_audit_write(diag, path, mode, status)
    type(fibre_prod_main_diagnostics_type), intent(in) :: diag
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: mode
    integer, intent(out) :: status
    integer :: unit_id
    logical :: finite_ok
    logical :: pass_ok
    real(dp) :: delta_sum
    real(dp) :: delta_max

    finite_ok = fibre_prod_main_diagnostics_finite(diag)
    delta_sum = abs(diag%after%sum_value - diag%before%sum_value)
    delta_max = abs(diag%after%max_abs - diag%before%max_abs)
    pass_ok = .false.
    if (trim(mode) == 'lambda0') then
      pass_ok = diag%hook_calls > 0 .and. diag%enabled .and. diag%lambda_fsi == 0.0_dp .and. &
                diag%modified_cells == 0 .and. diag%failed_calls == 0 .and. diag%last_status == 0 .and. &
                diag%no_contamination .and. finite_ok
    else if (trim(mode) == 'smalllambda') then
      pass_ok = diag%hook_calls > 0 .and. diag%enabled .and. diag%lambda_fsi > 0.0_dp .and. &
                diag%modified_cells > 0 .and. diag%failed_calls == 0 .and. diag%last_status == 0 .and. &
                diag%small_lambda_response .and. finite_ok
    end if

    status = 0
    open(newunit=unit_id, file=trim(path), status='replace', action='write', iostat=status)
    if (status /= 0) return
    if (trim(mode) == 'lambda0') then
      write(unit_id,'(A)') '# R10 Lambda=0 No-Contamination Audit'
    else
      write(unit_id,'(A)') '# R10 Small-Lambda Response Audit'
    end if
    write(unit_id,'(A)') ''
    if (pass_ok) then
      write(unit_id,'(A)') 'Result: PASS'
    else
      write(unit_id,'(A)') 'Result: FAIL'
    end if
    write(unit_id,'(A,I0)') 'hook_calls=', diag%hook_calls
    write(unit_id,'(A,L1)') 'enabled=', diag%enabled
    write(unit_id,'(A,ES24.16)') 'lambda_fsi=', diag%lambda_fsi
    write(unit_id,'(A,I0)') 'modified_cells=', diag%modified_cells
    write(unit_id,'(A,I0)') 'last_status=', diag%last_status
    write(unit_id,'(A,I0)') 'failed_calls=', diag%failed_calls
    write(unit_id,'(A,L1)') 'finite=', finite_ok
    write(unit_id,'(A,L1)') 'no_contamination=', diag%no_contamination
    write(unit_id,'(A,L1)') 'small_lambda_response=', diag%small_lambda_response
    write(unit_id,'(A,ES24.16)') 'delta_sum=', delta_sum
    write(unit_id,'(A,ES24.16)') 'delta_max_abs=', delta_max
    write(unit_id,'(A)') ''
    if (trim(mode) == 'lambda0') then
      write(unit_id,'(A)') 'Meaning: lambda=0 hook was called and did not modify the RHS buffer.'
    else
      write(unit_id,'(A)') 'Meaning: small lambda hook was called and produced a finite controlled RHS response.'
    end if
    close(unit_id)
  end subroutine fibre_prod_main_audit_write

end module fibre_prod_main_diagnostics
