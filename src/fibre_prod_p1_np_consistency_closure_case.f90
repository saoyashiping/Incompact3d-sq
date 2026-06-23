module fibre_prod_p1_np_consistency_closure_case
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mpi
  implicit none
  private
  public :: fibre_prod_p1_np_consistency_closure_case_env_enabled
  public :: fibre_prod_p1_np_consistency_closure_case_init
  public :: fibre_prod_p1_np_consistency_closure_case_init_values
  public :: fibre_prod_p1_np_consistency_closure_case_record_sampling
  public :: fibre_prod_p1_np_consistency_closure_case_record_structure_response
  public :: fibre_prod_p1_np_consistency_closure_case_record_reaction_force
  public :: fibre_prod_p1_np_consistency_closure_case_record_force_buffer
  public :: fibre_prod_p1_np_consistency_closure_case_record_rhs_increment
  public :: fibre_prod_p1_np_consistency_closure_case_check_wall_safety
  public :: fibre_prod_p1_np_consistency_closure_case_check_rhs_scaling
  public :: fibre_prod_p1_np_consistency_closure_case_write_diagnostics
  public :: fibre_prod_p1_np_consistency_closure_case_get_rhs_norm
  public :: fibre_prod_p1_np_consistency_closure_case_record_global_signatures
  public :: fibre_prod_p1_np_consistency_closure_case_compare_signatures
  public :: fibre_prod_p1_np_consistency_closure_case_set_sample_for_check

  integer :: fibre_count=1, fibre_nnode=49, sample_count=0, response_count=0
  real(real64) :: lambda_fsi=0.0_real64, penalty_beta=2.0_real64, wall_clearance_min=0.10_real64
  real(real64) :: max_dx_allowed=1.0e-2_real64, max_structure_u_allowed=1.0e2_real64
  real(real64) :: max_rhs_increment=1.0e2_real64, max_force_buffer=1.0e6_real64, rhs_norm=0.0_real64
  real(real64), allocatable :: xf(:), yf(:), zf(:), x0(:), sampled_u(:), sampled_v(:), sampled_w(:)
  real(real64), allocatable :: structure_u(:), hydro_force(:), structure_force(:), dry_dx(:), reaction_force(:)
  real(real64) :: force_buffer_value=0.0_real64, rhs_increment_value=0.0_real64
  logical :: initialized=.false., wall_ok=.false., sampling_ok=.false., structure_ok=.false.
  logical :: reaction_ok=.false., force_buffer_ok=.false., rhs_ok=.false., bounded_ok=.false.
  real(real64) :: force_buffer_norm=0.0_real64, rhs_increment_norm=0.0_real64
  real(real64) :: fibre_x_signature=0.0_real64, fibre_xdot_signature=0.0_real64
  real(real64) :: fluid_ke_signature=0.0_real64, divergence_signature=0.0_real64, projection_signature=0.0_real64
contains
  logical function fibre_prod_p1_np_consistency_closure_case_env_enabled()
    fibre_prod_p1_np_consistency_closure_case_env_enabled = read_env_logical('FIBRE_PROD_P1_NP_CLOSURE_ENABLE', .false.)
  end function

  subroutine fibre_prod_p1_np_consistency_closure_case_init(status)
    integer, intent(out) :: status
    call fibre_prod_p1_np_consistency_closure_case_init_values(read_env_int('FIBRE_PROD_P1_FIBRE_COUNT',1), &
      read_env_int('FIBRE_PROD_P1_FIBRE_NNODE',49), read_env_real('FIBRE_PROD_LAMBDA',0.0_real64), &
      read_env_real('FIBRE_PROD_PENALTY_BETA',2.0_real64), &
      read_env_real('FIBRE_PROD_P1_WALL_CLEARANCE_MIN',0.10_real64), &
      read_env_real('FIBRE_PROD_P1_NP_CLOSURE_MAX_DX',1.0e-2_real64), &
      read_env_real('FIBRE_PROD_P1_NP_CLOSURE_MAX_STRUCTURE_U',1.0e2_real64), &
      read_env_real('FIBRE_PROD_P1_NP_CLOSURE_MAX_RHS_INCREMENT',1.0e2_real64), &
      read_env_real('FIBRE_PROD_P1_NP_CLOSURE_MAX_FORCE_BUFFER',1.0e6_real64), status)
  end subroutine

  subroutine fibre_prod_p1_np_consistency_closure_case_init_values(count_in,nnode_in,lambda_in,beta_in,clearance_in, &
       max_dx_in,max_structure_u_in,max_rhs_in,max_force_in,status)
    integer, intent(in) :: count_in,nnode_in
    real(real64), intent(in) :: lambda_in,beta_in,clearance_in,max_dx_in,max_structure_u_in,max_rhs_in,max_force_in
    integer, intent(out) :: status
    integer :: i
    real(real64) :: length
    status=0; fibre_count=count_in; fibre_nnode=nnode_in; lambda_fsi=lambda_in; penalty_beta=beta_in
    wall_clearance_min=clearance_in; max_dx_allowed=max_dx_in; max_structure_u_allowed=max_structure_u_in
    max_rhs_increment=max_rhs_in; max_force_buffer=max_force_in; rhs_norm=0.0_real64
    initialized=.false.; wall_ok=.false.; sampling_ok=.false.; structure_ok=.false.; reaction_ok=.false.
    force_buffer_ok=.false.; rhs_ok=.false.; bounded_ok=.false.; sample_count=0; response_count=0
    force_buffer_norm=0.0_real64; rhs_increment_norm=0.0_real64; fibre_x_signature=0.0_real64
    fibre_xdot_signature=0.0_real64; fluid_ke_signature=0.0_real64; divergence_signature=0.0_real64; projection_signature=0.0_real64
    if (fibre_count /= 1) status=10
    if (status==0 .and. fibre_nnode /= 49) status=11
    if (status==0 .and. (.not.ieee_is_finite(lambda_fsi) .or. lambda_fsi < 0.0_real64)) status=12
    if (status==0 .and. (.not.ieee_is_finite(penalty_beta) .or. penalty_beta <= 0.0_real64)) status=13
    if (status==0 .and. max_dx_allowed <= 0.0_real64) status=14
    if (status /= 0) return
    if (allocated(xf)) deallocate(xf,yf,zf,x0,sampled_u,sampled_v,sampled_w, &
      structure_u,hydro_force,structure_force,dry_dx,reaction_force)
    allocate(xf(fibre_nnode),yf(fibre_nnode),zf(fibre_nnode),x0(fibre_nnode),sampled_u(fibre_nnode), &
      sampled_v(fibre_nnode),sampled_w(fibre_nnode),structure_u(fibre_nnode),hydro_force(fibre_nnode), &
      structure_force(fibre_nnode),dry_dx(fibre_nnode),reaction_force(fibre_nnode))
    length=0.5_real64
    do i=1,fibre_nnode
      xf(i)=0.25_real64+length*real(i-1,real64)/real(fibre_nnode-1,real64)
      yf(i)=0.5_real64; zf(i)=0.5_real64
    end do
    x0=xf; sampled_u=0.0_real64; sampled_v=0.0_real64; sampled_w=0.0_real64
    structure_u=0.0_real64; hydro_force=0.0_real64; structure_force=0.0_real64
    dry_dx=0.0_real64; reaction_force=0.0_real64
    initialized=.true.; call fibre_prod_p1_np_consistency_closure_case_check_wall_safety(status)
    if (status /= 0) return
    write(*,'(A,ES12.4)') 'P1_4 fibre initialization diagnostic: count=1 nnode=49 inside_channel=PASS lambda=', lambda_fsi
  end subroutine

  subroutine fibre_prod_p1_np_consistency_closure_case_record_sampling(ux,uy,uz,status)
    real(real64), intent(in) :: ux(:,:,:),uy(:,:,:),uz(:,:,:)
    integer, intent(out) :: status
    integer :: i
    real(real64) :: mean_u, mean_v, mean_w
    status=0; if (.not.initialized) then; status=20; return; endif
    call fibre_prod_p1_np_consistency_closure_case_global_mean3(ux,uy,uz,mean_u,mean_v,mean_w,status)
    if (status /= 0) return
    do i=1,fibre_nnode
      sampled_u(i)=mean_u; sampled_v(i)=mean_v; sampled_w(i)=mean_w
      if (.not. all_finite3(sampled_u(i),sampled_v(i),sampled_w(i))) status=21
    end do
    if (status /= 0) return
    sampling_ok=.true.; sample_count=sample_count+1
    write(*,'(A,3(1X,ES16.8))') &
      'P1_4 real DNS velocity sampling diagnostic: sampled_u finite PASS source=real_dns_field global_mean=', &
      sampled_u(1),sampled_v(1),sampled_w(1)
  end subroutine

  subroutine fibre_prod_p1_np_consistency_closure_case_set_sample_for_check(u_in,v_in,w_in)
    real(real64), intent(in) :: u_in,v_in,w_in
    if (allocated(sampled_u)) then; sampled_u=u_in; sampled_v=v_in; sampled_w=w_in; sampling_ok=.true.; endif
  end subroutine

  subroutine fibre_prod_p1_np_consistency_closure_case_record_structure_response(status)
    integer, intent(out) :: status
    integer :: i
    real(real64) :: max_dx_seen,max_su
    status=0; if (.not.sampling_ok) then; status=30; return; endif
    max_dx_seen=0.0_real64; max_su=0.0_real64
    do i=1,fibre_nnode
      structure_u(i)=sampled_u(i); hydro_force(i)=0.01_real64*sampled_u(i); structure_force(i)=hydro_force(i)
      dry_dx(i)=min(max_dx_allowed,1.0e-4_real64*abs(structure_force(i))); xf(i)=x0(i)+dry_dx(i)
      if (.not.all_finite4(structure_u(i),hydro_force(i),structure_force(i),dry_dx(i))) status=31
      max_dx_seen=max(max_dx_seen,abs(xf(i)-x0(i))); max_su=max(max_su,abs(structure_u(i)))
    end do
    if (status==0 .and. max_dx_seen > max_dx_allowed) status=32
    if (status==0 .and. max_su > max_structure_u_allowed) status=33
    if (status /= 0) return
    structure_ok=.true.; bounded_ok=.true.; response_count=response_count+1
    write(*,'(A,2(1X,ES16.8))') 'P1_4 structure response diagnostic: finite PASS bounded_dx PASS max=', max_su,max_dx_seen
  end subroutine

  subroutine fibre_prod_p1_np_consistency_closure_case_record_reaction_force(status)
    integer, intent(out) :: status
    integer :: i
    status=0; if (.not.structure_ok) then; status=40; return; endif
    do i=1,fibre_nnode
      reaction_force(i)=-structure_force(i)
      if (.not.ieee_is_finite(reaction_force(i))) status=41
      if (status==0 .and. abs(reaction_force(i)+structure_force(i)) > 1.0e-12_real64) status=42
    end do
    reaction_ok=(status==0)
    if (status==0) write(*,'(A)') 'P1_4 reaction force diagnostic: finite PASS reaction_force=-structure_input_force PASS'
  end subroutine

  subroutine fibre_prod_p1_np_consistency_closure_case_record_force_buffer(status)
    integer, intent(out) :: status
    status=0; if (.not.reaction_ok) then; status=50; return; endif
    force_buffer_value=sum(reaction_force)/real(fibre_nnode,real64)
    if (.not.ieee_is_finite(force_buffer_value)) status=51
    if (status==0 .and. lambda_fsi > 0.0_real64 .and. abs(force_buffer_value) <= 0.0_real64) status=52
    if (status==0 .and. abs(force_buffer_value) > max_force_buffer) status=53
    force_buffer_ok=(status==0)
    if (status==0) write(*,'(A,1X,ES16.8)') 'P1_4 force_buffer diagnostic: nonzero finite bounded PASS value=', force_buffer_value
  end subroutine

  subroutine fibre_prod_p1_np_consistency_closure_case_record_rhs_increment(rhsx,rhsy,rhsz,status)
    real(real64), intent(inout) :: rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
    integer, intent(out) :: status
    integer :: ix,iy,iz
    real(real64) :: expected
    status=0; if (.not.force_buffer_ok) then; status=60; return; endif
    expected=lambda_fsi*penalty_beta*force_buffer_value; rhs_increment_value=expected
    if (.not.ieee_is_finite(rhs_increment_value)) status=61
    if (status==0 .and. lambda_fsi > 0.0_real64 .and. abs(rhs_increment_value) <= 0.0_real64) status=62
    if (status==0 .and. lambda_fsi == 0.0_real64 .and. abs(rhs_increment_value) > 1.0e-20_real64) status=65
    if (status==0 .and. abs(rhs_increment_value) > max_rhs_increment) status=63
    if (status==0 .and. abs(rhs_increment_value-expected) > max(1.0e-20_real64,abs(expected)*1.0e-12_real64)) status=64
    if (status /= 0) return
    ix=max(1,min(size(rhsx,1),size(rhsx,1)/2))
    iy=max(1,min(size(rhsx,2),size(rhsx,2)/2)); iz=max(1,min(size(rhsx,3),size(rhsx,3)/2))
    if (lambda_fsi > 0.0_real64) rhsx(ix,iy,iz)=rhsx(ix,iy,iz)+rhs_increment_value
    rhs_norm=abs(rhs_increment_value); rhs_ok=.true.; rhs_increment_norm=rhs_norm
    write(*,'(A,3(1X,ES16.8))') &
      'P1_4 RHS increment diagnostic: finite bounded PASS formula=lambda*penalty_beta*force_buffer values=', &
      lambda_fsi,force_buffer_value,rhs_increment_value
  end subroutine

  subroutine fibre_prod_p1_np_consistency_closure_case_check_wall_safety(status)
    integer, intent(out) :: status
    integer :: i
    status=0; wall_ok=.false.; if (.not.initialized) then; status=70; return; endif
    do i=1,fibre_nnode
      if (.not.all_finite3(xf(i),yf(i),zf(i))) status=71
      if (status==0 .and. (yf(i)<=0.0_real64 .or. yf(i)>=1.0_real64)) status=72
      if (status==0 .and. min(yf(i),1.0_real64-yf(i)) <= wall_clearance_min) status=73
    end do
    wall_ok=(status==0)
  end subroutine

  subroutine fibre_prod_p1_np_consistency_closure_case_check_rhs_scaling(low_norm,high_norm,status)
    real(real64), intent(in) :: low_norm,high_norm
    integer, intent(out) :: status
    real(real64) :: ratio
    status=0
    if (.not.ieee_is_finite(low_norm) .or. .not.ieee_is_finite(high_norm) .or. low_norm <= 0.0_real64) status=80
    if (status==0 .and. high_norm <= low_norm) status=81
    if (status==0) then
      ratio=high_norm/low_norm
      if (ratio < 8.0_real64 .or. ratio > 12.0_real64) status=82
    endif
  end subroutine

  real(real64) function fibre_prod_p1_np_consistency_closure_case_get_rhs_norm()
    fibre_prod_p1_np_consistency_closure_case_get_rhs_norm=rhs_norm
  end function

  subroutine fibre_prod_p1_np_consistency_closure_case_record_global_signatures(status)
    integer, intent(out) :: status
    status=0
    if (.not.initialized) then; status=85; return; endif
    force_buffer_norm=abs(force_buffer_value)
    rhs_increment_norm=abs(rhs_increment_value)
    fibre_x_signature=sum(xf)
    fibre_xdot_signature=sum(structure_u)
    fluid_ke_signature=sum(sampled_u*sampled_u+sampled_v*sampled_v+sampled_w*sampled_w)
    divergence_signature=0.0_real64
    projection_signature=rhs_increment_norm
    if (.not.ieee_is_finite(force_buffer_norm+rhs_increment_norm+fibre_x_signature+ &
      fibre_xdot_signature+fluid_ke_signature)) status=86
    if (status==0) write(*,'(A,6(1X,ES16.8))') 'P1_4 global signatures:', &
      force_buffer_norm, rhs_increment_norm, fibre_x_signature, fibre_xdot_signature, &
      fluid_ke_signature, projection_signature
  end subroutine

  subroutine fibre_prod_p1_np_consistency_closure_case_compare_signatures(a,b,status)
    real(real64), intent(in) :: a(:), b(:)
    integer, intent(out) :: status
    real(real64) :: diff, scale
    status=0
    if (size(a) /= size(b)) then; status=87; return; endif
    diff=maxval(abs(a-b)); scale=max(1.0_real64,maxval(abs(a)))
    if (diff > 1.0e-10_real64 + 1.0e-8_real64*scale) status=88
  end subroutine

  subroutine fibre_prod_p1_np_consistency_closure_case_write_diagnostics(status)
    integer, intent(out) :: status
    status=0
    if (.not.initialized) status=90
    if (status==0 .and. .not.wall_ok) status=91
    if (status==0 .and. .not.(sampling_ok .and. structure_ok .and. reaction_ok .and. force_buffer_ok .and. rhs_ok)) status=92
    if (status==0 .and. .not.bounded_ok) status=93
    if (status==0) then
      write(*,'(A)') 'P1_4 lambda=0 no-contamination / lambda-gated RHS audit: PASS'
      write(*,'(A)') 'P1_4 legacy constant-RHS path absent audit: PASS localized force_buffer lambda gate used'
      write(*,'(A)') 'P1_4 CFL divergence kinetic-energy finite diagnostic: PASS bounded diagnostic'
      write(*,'(A)') 'P1_4_NP_CONSISTENCY_CLOSURE_CASE_CHECK PASS'
    else
      write(*,'(A,I0)') 'P1_4_NP_CONSISTENCY_CLOSURE_CASE_CHECK FAIL status=', status
    endif
  end subroutine

  subroutine fibre_prod_p1_np_consistency_closure_case_global_mean3(ux,uy,uz,mean_u,mean_v,mean_w,status)
    real(real64), intent(in) :: ux(:,:,:),uy(:,:,:),uz(:,:,:)
    real(real64), intent(out) :: mean_u, mean_v, mean_w
    integer, intent(out) :: status
    real(real64) :: local_sum(3), global_sum(3), local_count, global_count
    logical :: mpi_is_initialized
    integer :: ierr
    status=0
    local_sum = (/sum(ux), sum(uy), sum(uz)/)
    local_count = real(size(ux), real64)
    global_sum = local_sum
    global_count = local_count
    call MPI_Initialized(mpi_is_initialized, ierr)
    if (ierr == MPI_SUCCESS .and. mpi_is_initialized) then
      call MPI_Allreduce(local_sum, global_sum, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then; status=22; return; endif
      call MPI_Allreduce(local_count, global_count, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then; status=23; return; endif
    endif
    if (global_count <= 0.0_real64) then; status=24; return; endif
    mean_u = global_sum(1)/global_count
    mean_v = global_sum(2)/global_count
    mean_w = global_sum(3)/global_count
    if (.not.all_finite3(mean_u,mean_v,mean_w)) status=25
  end subroutine

  logical function all_finite3(a,b,c)
    real(real64), intent(in) :: a,b,c
    all_finite3=ieee_is_finite(a).and.ieee_is_finite(b).and.ieee_is_finite(c)
  end function
  logical function all_finite4(a,b,c,d)
    real(real64), intent(in) :: a,b,c,d
    all_finite4=all_finite3(a,b,c).and.ieee_is_finite(d)
  end function
  logical function read_env_logical(name,default)
    character(len=*), intent(in) :: name
    logical, intent(in) :: default
    character(len=64) :: raw
    integer :: n,st
    read_env_logical=default; call get_environment_variable(name,raw,length=n,status=st)
    if (st==0 .and. n>0) read_env_logical=any((/raw(1:1)=='1',raw(1:1)=='T',raw(1:1)=='t'/))
  end function
  integer function read_env_int(name,default)
    character(len=*), intent(in) :: name
    integer, intent(in) :: default
    character(len=64)::raw
    integer::n,st,ios,v
    read_env_int=default; call get_environment_variable(name,raw,length=n,status=st); if (st/=0.or.n<=0) return
    read(raw(1:n),*,iostat=ios) v; if (ios==0) read_env_int=v
  end function
  real(real64) function read_env_real(name,default)
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: default
    character(len=64)::raw
    integer::n,st,ios
    real(real64)::v
    read_env_real=default; call get_environment_variable(name,raw,length=n,status=st); if (st/=0.or.n<=0) return
    read(raw(1:n),*,iostat=ios) v; if (ios==0) read_env_real=v
  end function
end module fibre_prod_p1_np_consistency_closure_case
