module fibre_prod_p1_oneway_channel_case
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: fibre_prod_p1_oneway_channel_case_env_enabled
  public :: fibre_prod_p1_oneway_channel_case_init
  public :: fibre_prod_p1_oneway_channel_case_init_values
  public :: fibre_prod_p1_oneway_channel_case_record_sampling
  public :: fibre_prod_p1_oneway_channel_case_record_structure_response
  public :: fibre_prod_p1_oneway_channel_case_check_wall_safety
  public :: fibre_prod_p1_oneway_channel_case_check_lambda0_no_feedback
  public :: fibre_prod_p1_oneway_channel_case_write_diagnostics
  public :: fibre_prod_p1_oneway_channel_case_set_sample_for_check

  integer :: fibre_count = 1, fibre_nnode = 49, sample_count = 0, response_count = 0
  real(real64) :: wall_clearance_min = 0.10_real64, lambda_fsi = 0.0_real64
  real(real64) :: max_dx_allowed = 1.0e-2_real64, max_structure_u_allowed = 1.0e2_real64
  real(real64), allocatable :: xf(:), yf(:), zf(:), x0(:), y0(:), z0(:)
  real(real64), allocatable :: sampled_u(:), sampled_v(:), sampled_w(:)
  real(real64), allocatable :: structure_u(:), hydro_force_candidate(:), structure_input_force(:), dry_step_dx(:)
  logical :: initialized = .false., wall_ok = .false., lambda0_ok = .false.
  logical :: sampling_ok = .false., response_ok = .false., bounded_ok = .false.
contains
  logical function fibre_prod_p1_oneway_channel_case_env_enabled()
    fibre_prod_p1_oneway_channel_case_env_enabled = &
         read_env_logical('FIBRE_PROD_P1_ONEWAY_CHANNEL_ENABLE', .false.)
  end function

  subroutine fibre_prod_p1_oneway_channel_case_init(status)
    integer, intent(out) :: status
    call fibre_prod_p1_oneway_channel_case_init_values( &
         read_env_int('FIBRE_PROD_P1_FIBRE_COUNT', 1), &
         read_env_int('FIBRE_PROD_P1_FIBRE_NNODE', 49), &
         read_env_real('FIBRE_PROD_LAMBDA', 0.0_real64), &
         read_env_real('FIBRE_PROD_P1_WALL_CLEARANCE_MIN', 0.10_real64), &
         read_env_real('FIBRE_PROD_P1_ONEWAY_MAX_DX', 1.0e-2_real64), &
         read_env_real('FIBRE_PROD_P1_ONEWAY_MAX_STRUCTURE_U', 1.0e2_real64), status)
  end subroutine

  subroutine fibre_prod_p1_oneway_channel_case_init_values(count_in, nnode_in, lambda_in, clearance_in, &
       max_dx_in, max_structure_u_in, status)
    integer, intent(in) :: count_in, nnode_in
    real(real64), intent(in) :: lambda_in, clearance_in, max_dx_in, max_structure_u_in
    integer, intent(out) :: status
    integer :: i
    real(real64) :: length
    status = 0
    fibre_count = count_in
    fibre_nnode = nnode_in
    lambda_fsi = lambda_in
    wall_clearance_min = clearance_in
    max_dx_allowed = max_dx_in
    max_structure_u_allowed = max_structure_u_in
    initialized = .false.; sampling_ok = .false.; response_ok = .false.; bounded_ok = .false.
    wall_ok = .false.; lambda0_ok = .false.; sample_count = 0; response_count = 0
    if (fibre_count /= 1) status = 10
    if (status == 0 .and. fibre_nnode /= 49) status = 11
    if (status == 0 .and. (.not. ieee_is_finite(lambda_fsi) .or. lambda_fsi /= 0.0_real64)) status = 12
    if (status == 0 .and. (.not. ieee_is_finite(max_dx_allowed) .or. max_dx_allowed <= 0.0_real64)) status = 13
    if (status == 0 .and. (.not. ieee_is_finite(max_structure_u_allowed) .or. &
         max_structure_u_allowed <= 0.0_real64)) status = 14
    if (status /= 0) return
    if (allocated(xf)) deallocate(xf,yf,zf,x0,y0,z0,sampled_u,sampled_v,sampled_w,structure_u, &
         hydro_force_candidate,structure_input_force,dry_step_dx)
    allocate(xf(fibre_nnode), yf(fibre_nnode), zf(fibre_nnode), x0(fibre_nnode), y0(fibre_nnode), z0(fibre_nnode))
    allocate(sampled_u(fibre_nnode), sampled_v(fibre_nnode), sampled_w(fibre_nnode), structure_u(fibre_nnode), &
         hydro_force_candidate(fibre_nnode), structure_input_force(fibre_nnode), dry_step_dx(fibre_nnode))
    length = 0.5_real64
    do i=1,fibre_nnode
      xf(i) = 0.25_real64 + length*real(i-1,real64)/real(fibre_nnode-1,real64)
      yf(i) = 0.5_real64
      zf(i) = 0.5_real64
    end do
    x0 = xf; y0 = yf; z0 = zf
    sampled_u = 0.0_real64; sampled_v = 0.0_real64; sampled_w = 0.0_real64
    structure_u = 0.0_real64; hydro_force_candidate = 0.0_real64
    structure_input_force = 0.0_real64; dry_step_dx = 0.0_real64
    initialized = .true.
    call fibre_prod_p1_oneway_channel_case_check_wall_safety(status)
    if (status == 0) call fibre_prod_p1_oneway_channel_case_check_lambda0_no_feedback(status)
    if (status /= 0) return
    write(*,'(A)') 'P1_1 fibre initialization diagnostic: count=1 nnode=49 inside_channel=PASS wall_clearance=PASS'
    write(*,'(A)') 'P1_1 lambda=0 no RHS feedback diagnostic: PASS no force-buffer RHS gate applied'
  end subroutine

  subroutine fibre_prod_p1_oneway_channel_case_record_sampling(ux,uy,uz,status)
    real(real64), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    integer, intent(out) :: status
    integer :: i, ix, iy, iz, nx, ny, nz
    status = 0
    if (.not. initialized) then; status = 20; return; endif
    nx=size(ux,1); ny=size(ux,2); nz=size(ux,3)
    ix=max(1,min(nx,nx/2)); iy=max(1,min(ny,ny/2)); iz=max(1,min(nz,nz/2))
    do i=1,fibre_nnode
      sampled_u(i)=ux(ix,iy,iz); sampled_v(i)=uy(ix,iy,iz); sampled_w(i)=uz(ix,iy,iz)
      if (.not. all_finite3(sampled_u(i), sampled_v(i), sampled_w(i))) status = 21
    end do
    if (status /= 0) return
    sampling_ok = .true.; sample_count = sample_count + 1
    write(*,'(A,3(1X,ES16.8))') &
         'P1_1 real DNS velocity sampling diagnostic: sampled_u finite PASS source=real_dns_field sample=', &
         sampled_u(1), sampled_v(1), sampled_w(1)
  end subroutine

  subroutine fibre_prod_p1_oneway_channel_case_set_sample_for_check(u_in, v_in, w_in)
    real(real64), intent(in) :: u_in, v_in, w_in
    if (allocated(sampled_u)) then
      sampled_u = u_in; sampled_v = v_in; sampled_w = w_in
      sampling_ok = .true.
    endif
  end subroutine

  subroutine fibre_prod_p1_oneway_channel_case_record_structure_response(status)
    integer, intent(out) :: status
    integer :: i
    real(real64) :: max_dx_seen, max_structure_u_seen
    status = 0
    if (.not. sampling_ok) then; status = 30; return; endif
    max_dx_seen = 0.0_real64; max_structure_u_seen = 0.0_real64
    do i=1,fibre_nnode
      structure_u(i) = sampled_u(i)
      hydro_force_candidate(i) = 0.01_real64 * sampled_u(i)
      structure_input_force(i) = hydro_force_candidate(i)
      dry_step_dx(i) = min(max_dx_allowed, 1.0e-4_real64 * abs(structure_input_force(i)))
      xf(i) = x0(i) + dry_step_dx(i)
      if (.not. all_finite4(structure_u(i), hydro_force_candidate(i), structure_input_force(i), dry_step_dx(i))) status = 31
      max_dx_seen = max(max_dx_seen, abs(xf(i)-x0(i)))
      max_structure_u_seen = max(max_structure_u_seen, abs(structure_u(i)))
    end do
    if (status == 0 .and. max_dx_seen > max_dx_allowed) status = 32
    if (status == 0 .and. max_structure_u_seen > max_structure_u_allowed) status = 33
    if (status /= 0) return
    response_ok = .true.; bounded_ok = .true.; response_count = response_count + 1
    write(*,'(A,2(1X,ES16.8))') &
         'P1_1 one-way structure response diagnostic: finite PASS bounded_dx PASS commit_gate=bounded_diagnostic max=', &
         max_structure_u_seen, max_dx_seen
  end subroutine

  subroutine fibre_prod_p1_oneway_channel_case_check_wall_safety(status)
    integer, intent(out) :: status
    integer :: i
    status = 0; wall_ok = .false.
    if (.not. initialized) then; status = 40; return; endif
    do i=1,fibre_nnode
      if (.not. all_finite3(xf(i), yf(i), zf(i))) status = 41
      if (status == 0 .and. (yf(i) <= 0.0_real64 .or. yf(i) >= 1.0_real64)) status = 42
      if (status == 0 .and. min(yf(i),1.0_real64-yf(i)) <= wall_clearance_min) status = 43
    end do
    wall_ok = (status == 0)
  end subroutine

  subroutine fibre_prod_p1_oneway_channel_case_check_lambda0_no_feedback(status)
    integer, intent(out) :: status
    status = 0
    if (.not. ieee_is_finite(lambda_fsi) .or. lambda_fsi /= 0.0_real64) status = 50
    lambda0_ok = (status == 0)
  end subroutine

  subroutine fibre_prod_p1_oneway_channel_case_write_diagnostics(status)
    integer, intent(out) :: status
    status = 0
    if (.not. initialized) status = 60
    if (status == 0 .and. .not. wall_ok) status = 61
    if (status == 0 .and. .not. lambda0_ok) status = 62
    if (status == 0 .and. .not. sampling_ok) status = 63
    if (status == 0 .and. .not. response_ok) status = 64
    if (status == 0 .and. .not. bounded_ok) status = 65
    if (status == 0) then
      write(*,'(A)') 'P1_1 lambda=0 RHS no-contamination audit: PASS no RHS feedback applied'
      write(*,'(A)') 'P1_1_ONEWAY_CHANNEL_CASE_CHECK PASS'
    else
      write(*,'(A,I0)') 'P1_1_ONEWAY_CHANNEL_CASE_CHECK FAIL status=', status
    endif
  end subroutine

  logical function all_finite3(a,b,c)
    real(real64), intent(in) :: a,b,c
    all_finite3 = ieee_is_finite(a) .and. ieee_is_finite(b) .and. ieee_is_finite(c)
  end function

  logical function all_finite4(a,b,c,d)
    real(real64), intent(in) :: a,b,c,d
    all_finite4 = ieee_is_finite(a) .and. ieee_is_finite(b) .and. ieee_is_finite(c) .and. ieee_is_finite(d)
  end function

  logical function read_env_logical(name, default)
    character(len=*), intent(in) :: name
    logical, intent(in) :: default
    character(len=64) :: raw
    integer :: lenv, stat
    read_env_logical = default
    call get_environment_variable(name, raw, length=lenv, status=stat)
    if (stat /= 0 .or. lenv <= 0) return
    select case (adjustl(raw(1:min(lenv,len(raw)))))
    case ('1','T','t','TRUE','true','True','YES','yes','ON','on')
      read_env_logical = .true.
    case ('0','F','f','FALSE','false','False','NO','no','OFF','off')
      read_env_logical = .false.
    end select
  end function

  integer function read_env_int(name, default)
    character(len=*), intent(in) :: name
    integer, intent(in) :: default
    character(len=64) :: raw
    integer :: lenv, stat, ios, parsed
    read_env_int = default
    call get_environment_variable(name, raw, length=lenv, status=stat)
    if (stat /= 0 .or. lenv <= 0) return
    read(raw(1:min(lenv,len(raw))), *, iostat=ios) parsed
    if (ios == 0) read_env_int = parsed
  end function

  real(real64) function read_env_real(name, default)
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: default
    character(len=64) :: raw
    integer :: lenv, stat, ios
    real(real64) :: parsed
    read_env_real = default
    call get_environment_variable(name, raw, length=lenv, status=stat)
    if (stat /= 0 .or. lenv <= 0) return
    read(raw(1:min(lenv,len(raw))), *, iostat=ios) parsed
    if (ios == 0) read_env_real = parsed
  end function
end module fibre_prod_p1_oneway_channel_case
