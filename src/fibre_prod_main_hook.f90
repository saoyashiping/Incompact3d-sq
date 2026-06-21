module fibre_prod_main_hook
  use, intrinsic :: iso_fortran_env, only : real64
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, &
                                        fibre_prod_runtime_config_load_env, &
                                        fibre_prod_runtime_config_validate
  use fibre_prod_main_diagnostics, only : fibre_prod_main_diagnostics_type, &
                                          fibre_prod_rhs_signature_type, &
                                          fibre_prod_main_signature, &
                                          fibre_prod_main_diagnostics_reset, &
                                          fibre_prod_main_diagnostics_record, &
                                          fibre_prod_main_diagnostics_write, &
                                          fibre_prod_main_audit_write
  use fibre_prod_rhs_adapter, only : fibre_prod_rhs_adapter_apply
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, fibre_prod_force_buffer_is_finite
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_main_hook_init
  public :: fibre_prod_main_hook_apply
  public :: fibre_prod_main_hook_apply_force_buffer
  public :: fibre_prod_main_hook_runtime_bridge_ready
  public :: fibre_prod_main_hook_get_runtime_config_status
  public :: fibre_prod_main_hook_finalize
  public :: fibre_prod_main_hook_get_diagnostics

  type(fibre_prod_runtime_config_type), save :: runtime_config
  type(fibre_prod_main_diagnostics_type), save :: runtime_diagnostics
  logical, save :: hook_initialized = .false.

contains

  subroutine fibre_prod_main_hook_init(status, config_override)
    integer, intent(out) :: status
    type(fibre_prod_runtime_config_type), intent(in), optional :: config_override

    if (present(config_override)) then
      runtime_config = config_override
      status = fibre_prod_runtime_config_validate(runtime_config)
    else
      call fibre_prod_runtime_config_load_env(runtime_config, status)
    end if
    call fibre_prod_main_diagnostics_reset(runtime_diagnostics)
    hook_initialized = status == 0
  end subroutine fibre_prod_main_hook_init

  subroutine fibre_prod_main_hook_apply(rhs_x, rhs_y, rhs_z, status, force_x, force_y, force_z)
    real(dp), intent(inout) :: rhs_x(:, :, :)
    real(dp), intent(inout) :: rhs_y(:, :, :)
    real(dp), intent(inout) :: rhs_z(:, :, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: force_x(:, :, :)
    real(dp), intent(in), optional :: force_y(:, :, :)
    real(dp), intent(in), optional :: force_z(:, :, :)
    type(fibre_prod_rhs_signature_type) :: before
    type(fibre_prod_rhs_signature_type) :: after
    integer :: modified_cells

    if (.not. hook_initialized) then
      status = 1
      return
    end if
    status = fibre_prod_runtime_config_validate(runtime_config)
    if (status /= 0) return
    call fibre_prod_rhs_adapter_apply(runtime_config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status, &
                                      force_x, force_y, force_z)
    call fibre_prod_main_diagnostics_record(runtime_diagnostics, runtime_config%enabled, &
                                            runtime_config%lambda_fsi, before, after, modified_cells, status)
  end subroutine fibre_prod_main_hook_apply


  subroutine fibre_prod_main_hook_apply_force_buffer(rhs_x, rhs_y, rhs_z, force_buffer, status)
    real(dp), intent(inout) :: rhs_x(:, :, :)
    real(dp), intent(inout) :: rhs_y(:, :, :)
    real(dp), intent(inout) :: rhs_z(:, :, :)
    type(fibre_prod_force_buffer_type), intent(in) :: force_buffer
    integer, intent(out) :: status
    type(fibre_prod_rhs_signature_type) :: before
    type(fibre_prod_rhs_signature_type) :: after
    integer :: modified_cells

    before = fibre_prod_main_signature(rhs_x, rhs_y, rhs_z)
    after = before
    modified_cells = 0
    if (.not. hook_initialized) then
      status = 1
      return
    end if
    status = fibre_prod_runtime_config_validate(runtime_config)
    if (status /= 0) return
    if (size(rhs_y, 1) /= size(rhs_x, 1) .or. size(rhs_y, 2) /= size(rhs_x, 2) .or. &
        size(rhs_y, 3) /= size(rhs_x, 3) .or. size(rhs_z, 1) /= size(rhs_x, 1) .or. &
        size(rhs_z, 2) /= size(rhs_x, 2) .or. size(rhs_z, 3) /= size(rhs_x, 3)) then
      status = 23
      call fibre_prod_main_diagnostics_record(runtime_diagnostics, runtime_config%enabled, &
                                              runtime_config%lambda_fsi, before, after, modified_cells, status)
      return
    end if
    if (.not. force_buffer%allocated .or. .not. allocated(force_buffer%fx) .or. &
        .not. allocated(force_buffer%fy) .or. .not. allocated(force_buffer%fz) .or. &
        .not. fibre_prod_force_buffer_is_finite(force_buffer)) then
      status = 21
      call fibre_prod_main_diagnostics_record(runtime_diagnostics, runtime_config%enabled, &
                                              runtime_config%lambda_fsi, before, after, modified_cells, status)
      return
    end if
    if (size(force_buffer%fx, 1) /= size(rhs_x, 1) .or. size(force_buffer%fx, 2) /= size(rhs_x, 2) .or. &
        size(force_buffer%fx, 3) /= size(rhs_x, 3) .or. size(force_buffer%fy, 1) /= size(rhs_y, 1) .or. &
        size(force_buffer%fy, 2) /= size(rhs_y, 2) .or. size(force_buffer%fy, 3) /= size(rhs_y, 3) .or. &
        size(force_buffer%fz, 1) /= size(rhs_z, 1) .or. size(force_buffer%fz, 2) /= size(rhs_z, 2) .or. &
        size(force_buffer%fz, 3) /= size(rhs_z, 3)) then
      status = 22
      call fibre_prod_main_diagnostics_record(runtime_diagnostics, runtime_config%enabled, &
                                              runtime_config%lambda_fsi, before, after, modified_cells, status)
      return
    end if

    call fibre_prod_rhs_adapter_apply(runtime_config, rhs_x, rhs_y, rhs_z, before, after, modified_cells, status, &
                                      force_buffer%fx, force_buffer%fy, force_buffer%fz)
    call fibre_prod_main_diagnostics_record(runtime_diagnostics, runtime_config%enabled, &
                                            runtime_config%lambda_fsi, before, after, modified_cells, status)
  end subroutine fibre_prod_main_hook_apply_force_buffer


  logical function fibre_prod_main_hook_runtime_bridge_ready() result(is_ready)

    is_ready = hook_initialized
  end function fibre_prod_main_hook_runtime_bridge_ready

  subroutine fibre_prod_main_hook_get_runtime_config_status(enabled, lambda_fsi, status)
    logical, intent(out) :: enabled
    real(dp), intent(out) :: lambda_fsi
    integer, intent(out) :: status

    enabled = runtime_config%enabled
    lambda_fsi = runtime_config%lambda_fsi
    if (.not. hook_initialized) then
      status = 1
    else
      status = fibre_prod_runtime_config_validate(runtime_config)
    end if
  end subroutine fibre_prod_main_hook_get_runtime_config_status

  subroutine fibre_prod_main_hook_finalize(status)
    integer, intent(out) :: status
    character(len=512) :: diagnostics_file
    character(len=512) :: audit_file
    character(len=512) :: mode_file
    integer :: write_status

    status = 0
    if (.not. hook_initialized) return
    if (.not. runtime_config%diagnostics_enabled) return

    diagnostics_file = join_path(runtime_config%diagnostics_dir, 'R10_MAIN_HOOK_DIAGNOSTICS.txt')
    call fibre_prod_main_diagnostics_write(runtime_diagnostics, diagnostics_file, write_status)
    if (write_status /= 0 .and. status == 0) status = write_status

    if (runtime_config%enabled .and. runtime_config%lambda_fsi == 0.0_dp) then
      audit_file = join_path(runtime_config%diagnostics_dir, 'R10_LAMBDA0_NO_CONTAMINATION_AUDIT.txt')
      mode_file = join_path(runtime_config%diagnostics_dir, 'R10_MAIN_HOOK_DIAGNOSTICS_lambda0.txt')
      call fibre_prod_main_diagnostics_write(runtime_diagnostics, mode_file, write_status)
      if (write_status /= 0 .and. status == 0) status = write_status
      call fibre_prod_main_audit_write(runtime_diagnostics, audit_file, 'lambda0', write_status)
      if (write_status /= 0 .and. status == 0) status = write_status
    else if (runtime_config%enabled .and. runtime_config%lambda_fsi > 0.0_dp) then
      audit_file = join_path(runtime_config%diagnostics_dir, 'R10_SMALL_LAMBDA_RESPONSE_AUDIT.txt')
      mode_file = join_path(runtime_config%diagnostics_dir, 'R10_MAIN_HOOK_DIAGNOSTICS_smalllambda.txt')
      call fibre_prod_main_diagnostics_write(runtime_diagnostics, mode_file, write_status)
      if (write_status /= 0 .and. status == 0) status = write_status
      call fibre_prod_main_audit_write(runtime_diagnostics, audit_file, 'smalllambda', write_status)
      if (write_status /= 0 .and. status == 0) status = write_status
    end if
  end subroutine fibre_prod_main_hook_finalize

  subroutine fibre_prod_main_hook_get_diagnostics(diag)
    type(fibre_prod_main_diagnostics_type), intent(out) :: diag

    diag = runtime_diagnostics
  end subroutine fibre_prod_main_hook_get_diagnostics

  pure function join_path(directory, filename) result(path)
    character(len=*), intent(in) :: directory
    character(len=*), intent(in) :: filename
    character(len=512) :: path
    integer :: n

    path = ''
    n = len_trim(directory)
    if (n <= 0 .or. trim(directory) == '.') then
      path = trim(filename)
    else if (directory(n:n) == '/') then
      path = trim(directory)//trim(filename)
    else
      path = trim(directory)//'/'//trim(filename)
    end if
  end function join_path

end module fibre_prod_main_hook
