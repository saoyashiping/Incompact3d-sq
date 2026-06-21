module fibre_prod_main_hook
  use, intrinsic :: iso_fortran_env, only : real64
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, &
                                        fibre_prod_runtime_config_load_env, &
                                        fibre_prod_runtime_config_validate
  use fibre_prod_main_diagnostics, only : fibre_prod_main_diagnostics_type, &
                                          fibre_prod_rhs_signature_type, &
                                          fibre_prod_main_diagnostics_reset, &
                                          fibre_prod_main_diagnostics_record, &
                                          fibre_prod_main_diagnostics_write
  use fibre_prod_rhs_adapter, only : fibre_prod_rhs_adapter_apply
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_main_hook_init
  public :: fibre_prod_main_hook_apply
  public :: fibre_prod_main_hook_finalize
  public :: fibre_prod_main_hook_get_diagnostics

  type(fibre_prod_runtime_config_type), save :: runtime_config
  type(fibre_prod_main_diagnostics_type), save :: runtime_diagnostics
  logical, save :: hook_initialized = .false.

contains

  subroutine fibre_prod_main_hook_init(status)
    integer, intent(out) :: status

    call fibre_prod_runtime_config_load_env(runtime_config, status)
    call fibre_prod_main_diagnostics_reset(runtime_diagnostics)
    hook_initialized = status == 0
  end subroutine fibre_prod_main_hook_init

  subroutine fibre_prod_main_hook_apply(rhs_x, rhs_y, rhs_z, status)
    real(dp), intent(inout) :: rhs_x(:, :, :)
    real(dp), intent(inout) :: rhs_y(:, :, :)
    real(dp), intent(inout) :: rhs_z(:, :, :)
    integer, intent(out) :: status
    type(fibre_prod_rhs_signature_type) :: before
    type(fibre_prod_rhs_signature_type) :: after
    integer :: modified_cells

    if (.not. hook_initialized) then
      status = 1
      return
    end if
    status = fibre_prod_runtime_config_validate(runtime_config)
    if (status /= 0) return
    call fibre_prod_rhs_adapter_apply(runtime_config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status)
    if (status /= 0) return
    call fibre_prod_main_diagnostics_record(runtime_diagnostics, runtime_config%enabled, &
                                            runtime_config%lambda_fsi, before, after, modified_cells)
  end subroutine fibre_prod_main_hook_apply

  subroutine fibre_prod_main_hook_finalize(status)
    integer, intent(out) :: status

    status = 0
    if (.not. hook_initialized) return
    if (runtime_config%diagnostics_enabled) then
      call fibre_prod_main_diagnostics_write(runtime_diagnostics, &
           'production_recovery/R10_evidence/R10_MAIN_HOOK_DIAGNOSTICS.txt', status)
    end if
  end subroutine fibre_prod_main_hook_finalize

  subroutine fibre_prod_main_hook_get_diagnostics(diag)
    type(fibre_prod_main_diagnostics_type), intent(out) :: diag

    diag = runtime_diagnostics
  end subroutine fibre_prod_main_hook_get_diagnostics

end module fibre_prod_main_hook
