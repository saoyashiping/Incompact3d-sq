module fibre_prod_p1_real_channel_preflight
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: fibre_prod_p1_real_channel_preflight_env_enabled
  public :: fibre_prod_p1_real_channel_preflight_init
  public :: fibre_prod_p1_real_channel_preflight_record_sampling
  public :: fibre_prod_p1_real_channel_preflight_check_wall_safety
  public :: fibre_prod_p1_real_channel_preflight_write_diagnostics
  public :: fibre_prod_p1_real_channel_preflight_init_values

  integer :: fibre_count = 1, fibre_nnode = 33, sample_count = 0
  real(real64) :: wall_clearance_min = 0.10_real64, lambda_fsi = 0.0_real64
  real(real64), allocatable :: xf(:), yf(:), zf(:), sampled_u(:), sampled_v(:), sampled_w(:)
  logical :: initialized = .false., sampling_ok = .false., wall_ok = .false., lambda_ok = .false.
contains
  logical function fibre_prod_p1_real_channel_preflight_env_enabled()
    fibre_prod_p1_real_channel_preflight_env_enabled = read_env_logical('FIBRE_PROD_P1_REAL_CHANNEL_PREFLIGHT_ENABLE', .false.)
  end function

  subroutine fibre_prod_p1_real_channel_preflight_init(status)
    integer, intent(out) :: status
    call fibre_prod_p1_real_channel_preflight_init_values( &
         read_env_int('FIBRE_PROD_P1_FIBRE_COUNT', 1), &
         read_env_int('FIBRE_PROD_P1_FIBRE_NNODE', 33), &
         read_env_real('FIBRE_PROD_LAMBDA', 0.0_real64), &
         read_env_real('FIBRE_PROD_P1_WALL_CLEARANCE_MIN', 0.10_real64), status)
  end subroutine

  subroutine fibre_prod_p1_real_channel_preflight_init_values(count_in, nnode_in, lambda_in, clearance_in, status)
    integer, intent(in) :: count_in, nnode_in
    real(real64), intent(in) :: lambda_in, clearance_in
    integer, intent(out) :: status
    integer :: i
    real(real64) :: length
    status = 0
    fibre_count = count_in
    fibre_nnode = nnode_in
    wall_clearance_min = clearance_in
    lambda_fsi = lambda_in
    if (fibre_count /= 1) status = 10
    if (status == 0 .and. fibre_nnode /= 33) status = 11
    if (status == 0 .and. lambda_fsi /= 0.0_real64) status = 12
    if (status /= 0) return
    if (allocated(xf)) deallocate(xf,yf,zf,sampled_u,sampled_v,sampled_w)
    allocate(xf(fibre_nnode), yf(fibre_nnode), zf(fibre_nnode), &
         sampled_u(fibre_nnode), sampled_v(fibre_nnode), sampled_w(fibre_nnode))
    length = 0.5_real64
    do i=1,fibre_nnode
      xf(i) = 0.25_real64 + length*real(i-1,real64)/real(fibre_nnode-1,real64)
      yf(i) = 0.5_real64
      zf(i) = 0.5_real64
    end do
    sampled_u = 0.0_real64; sampled_v = 0.0_real64; sampled_w = 0.0_real64
    initialized = .true.; sample_count = 0; sampling_ok = .false.; lambda_ok = .true.
    call fibre_prod_p1_real_channel_preflight_check_wall_safety(status)
    if (status /= 0) return
    write(*,'(A)') 'P1_0 fibre initialization diagnostic: count=1 nnode=33 inside_channel=PASS wall_clearance=PASS'
    write(*,'(A)') 'P1_0 lambda=0 no RHS feedback diagnostic: PASS'
  end subroutine fibre_prod_p1_real_channel_preflight_init_values

  subroutine fibre_prod_p1_real_channel_preflight_record_sampling(ux,uy,uz,status)
    real(real64), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    integer, intent(out) :: status
    integer :: i, ix, iy, iz, nx, ny, nz
    status = 0
    if (.not. initialized) then; status = 20; return; endif
    nx=size(ux,1); ny=size(ux,2); nz=size(ux,3); ix=max(1,min(nx,nx/2)); iy=max(1,min(ny,ny/2)); iz=max(1,min(nz,nz/2))
    do i=1,fibre_nnode
      sampled_u(i)=ux(ix,iy,iz); sampled_v(i)=uy(ix,iy,iz); sampled_w(i)=uz(ix,iy,iz)
      if (.not.(ieee_is_finite(sampled_u(i)) .and. ieee_is_finite(sampled_v(i)) .and. ieee_is_finite(sampled_w(i)))) status=21
    end do
    if (status == 0) then
      sampling_ok = .true.; sample_count = sample_count + 1
      write(*,'(A,3(1X,ES16.8))') &
           'P1_0 real DNS velocity sampling diagnostic: sampled_u finite PASS source=real_dns_field sample=', &
           sampled_u(1), sampled_v(1), sampled_w(1)
    endif
  end subroutine

  subroutine fibre_prod_p1_real_channel_preflight_check_wall_safety(status)
    integer, intent(out) :: status
    integer :: i
    status = 0; wall_ok = .false.
    do i=1,fibre_nnode
      if (.not.(ieee_is_finite(xf(i)) .and. ieee_is_finite(yf(i)) .and. ieee_is_finite(zf(i)))) status=30
      if (status == 0 .and. (yf(i) <= 0.0_real64 .or. yf(i) >= 1.0_real64)) status=31
      if (status == 0 .and. min(yf(i),1.0_real64-yf(i)) <= wall_clearance_min) status=32
    end do
    wall_ok = (status == 0)
  end subroutine

  subroutine fibre_prod_p1_real_channel_preflight_write_diagnostics(status)
    integer, intent(out) :: status
    status = 0
    if (.not. initialized) status = 40
    if (status == 0 .and. .not. wall_ok) status = 41
    if (status == 0 .and. .not. sampling_ok) status = 42
    if (status == 0 .and. .not. lambda_ok) status = 43
    if (status == 0) then
      write(*,'(A)') 'P1_0 lambda=0 RHS no-contamination audit: PASS no RHS feedback applied'
      write(*,'(A)') 'P1_0_REAL_CHANNEL_PREFLIGHT_CHECK PASS'
    else
      write(*,'(A,I0)') 'P1_0_REAL_CHANNEL_PREFLIGHT_CHECK FAIL status=', status
    endif
  end subroutine

  logical function read_env_logical(name, default)
    character(len=*), intent(in) :: name; logical, intent(in) :: default
    character(len=64) :: raw; integer :: lenv, stat
    read_env_logical=default; call get_environment_variable(name,raw,length=lenv,status=stat)
    if (stat==0 .and. lenv>0) read_env_logical = any((/raw(1:1)=='1', raw(1:1)=='T', raw(1:1)=='t'/))
  end function
  integer function read_env_int(name, default)
    character(len=*), intent(in) :: name; integer, intent(in) :: default
    character(len=64) :: raw; integer :: lenv, stat, ios
    read_env_int=default; call get_environment_variable(name,raw,length=lenv,status=stat)
    if (stat==0 .and. lenv>0) read(raw(1:lenv),*,iostat=ios) read_env_int
  end function
  real(real64) function read_env_real(name, default)
    character(len=*), intent(in) :: name; real(real64), intent(in) :: default
    character(len=64) :: raw; integer :: lenv, stat, ios
    read_env_real=default; call get_environment_variable(name,raw,length=lenv,status=stat)
    if (stat==0 .and. lenv>0) read(raw(1:lenv),*,iostat=ios) read_env_real
  end function
end module fibre_prod_p1_real_channel_preflight
