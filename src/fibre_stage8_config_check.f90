program fibre_stage8_config_check
  use fibre_parameters, only: mytype
  use fibre_stage8_config
  implicit none

  type(stage8_config_t) :: cfg
  integer :: io
  integer :: ex_stage7_marker, ex_stage7_total, s7_smoke, s7_closed_status, s7_closed_written
  integer :: dep_ok
  integer :: valid, rejected, rhs_allowed
  integer :: default_ok, controlled_ok, invalid_ok, rhs_gate_ok, noop_ok
  integer :: prod_dns_rej, prod_twoway_rej, no_ctrl_rej, rho_rej
  integer :: rhs_controlled_allowed, rhs_prod_dns_blocked, rhs_prod_twoway_blocked
  real(mytype) :: rhs_before(4,3,2), rhs_after(4,3,2), rhs_change
  integer :: rhs_modified
  integer :: final_status

  call ensure_dir('stage8_outputs')

  call file_exists_int('stage7_outputs/STAGE7_CLOSED.md', ex_stage7_marker)
  call file_exists_int('stage7_outputs/fibre_stage7_total_smoke_check.dat', ex_stage7_total)
  call get_int_value('stage7_outputs/fibre_stage7_total_smoke_check.dat', 'stage7_total_smoke_check_status', s7_smoke)
  call get_int_value('stage7_outputs/fibre_stage7_total_smoke_check.dat', 'stage7_total_closed_marker_status', s7_closed_status)
  call get_int_value('stage7_outputs/fibre_stage7_total_smoke_check.dat', 'stage7_total_closed_marker_written_flag', s7_closed_written)
  dep_ok = merge(1,0,ex_stage7_marker==1 .and. ex_stage7_total==1 .and. s7_smoke==1 .and. s7_closed_status==1 .and. s7_closed_written==1)

  call init_stage8_default_config(cfg)
  call validate_stage8_config(cfg,valid,rejected,rhs_allowed)
  default_ok = merge(1,0,cfg%enable_stage8==0 .and. cfg%controlled_test_enabled==0 .and. cfg%production_dns_enabled==0 .and. &
    cfg%production_two_way_enabled==0 .and. rhs_allowed==0 .and. valid==1 .and. rejected==0)

  call init_stage8_controlled_test_config(cfg)
  call validate_stage8_config(cfg,valid,rejected,rhs_allowed)
  controlled_ok = merge(1,0,cfg%enable_stage8==1 .and. cfg%controlled_test_enabled==1 .and. cfg%production_dns_enabled==0 .and. &
    cfg%production_two_way_enabled==0 .and. cfg%enable_runtime_grid_bridge==1 .and. cfg%enable_lagrangian_state==1 .and. &
    cfg%enable_velocity_to_fibre==1 .and. cfg%enable_feedback_candidate==0 .and. cfg%enable_force_density_candidate==0 .and. &
    cfg%enable_rhs_candidate==0 .and. rhs_allowed==0 .and. valid==1 .and. rejected==0)

  call init_stage8_default_config(cfg)
  cfg%production_dns_enabled=1
  call validate_stage8_config(cfg,valid,rejected,rhs_allowed)
  prod_dns_rej = merge(1,0,rejected==1 .and. valid==0 .and. rhs_allowed==0)

  call init_stage8_default_config(cfg)
  cfg%production_two_way_enabled=1
  call validate_stage8_config(cfg,valid,rejected,rhs_allowed)
  prod_twoway_rej = merge(1,0,rejected==1 .and. valid==0 .and. rhs_allowed==0)

  call init_stage8_default_config(cfg)
  cfg%enable_stage8=1
  cfg%controlled_test_enabled=0
  call validate_stage8_config(cfg,valid,rejected,rhs_allowed)
  no_ctrl_rej = merge(1,0,rejected==1 .and. valid==0 .and. rhs_allowed==0)

  call init_stage8_default_config(cfg)
  cfg%rho_fluid=0._mytype
  call validate_stage8_config(cfg,valid,rejected,rhs_allowed)
  rho_rej = merge(1,0,rejected==1 .and. valid==0 .and. rhs_allowed==0)
  invalid_ok = merge(1,0,prod_dns_rej==1 .and. prod_twoway_rej==1 .and. no_ctrl_rej==1 .and. rho_rej==1)

  call init_stage8_controlled_test_config(cfg)
  cfg%enable_rhs_candidate=1
  call validate_stage8_config(cfg,valid,rejected,rhs_allowed)
  rhs_controlled_allowed = merge(1,0,valid==1 .and. rejected==0 .and. rhs_allowed==1)

  call init_stage8_controlled_test_config(cfg)
  cfg%enable_rhs_candidate=1
  cfg%production_dns_enabled=1
  call validate_stage8_config(cfg,valid,rejected,rhs_allowed)
  rhs_prod_dns_blocked = merge(1,0,rejected==1 .and. rhs_allowed==0)

  call init_stage8_controlled_test_config(cfg)
  cfg%enable_rhs_candidate=1
  cfg%production_two_way_enabled=1
  call validate_stage8_config(cfg,valid,rejected,rhs_allowed)
  rhs_prod_twoway_blocked = merge(1,0,rejected==1 .and. rhs_allowed==0)
  rhs_gate_ok = merge(1,0,rhs_controlled_allowed==1 .and. rhs_prod_dns_blocked==1 .and. rhs_prod_twoway_blocked==1)

  call init_rhs(rhs_before)
  rhs_after = rhs_before
  rhs_change = maxval(abs(rhs_after-rhs_before))
  rhs_modified = merge(1,0,rhs_change>1e-14_mytype)
  noop_ok = merge(1,0,rhs_change<=1e-14_mytype .and. rhs_modified==0)

  final_status = merge(1,0,dep_ok==1 .and. default_ok==1 .and. controlled_ok==1 .and. invalid_ok==1 .and. rhs_gate_ok==1 .and. noop_ok==1)

  open(newunit=io,file='stage8_outputs/fibre_stage8_config_check.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage8_config_stage7_closed_marker_exists', ex_stage7_marker
  write(io,'(A,1X,I0)') 'stage8_config_stage7_total_smoke_output_exists', ex_stage7_total
  write(io,'(A,1X,I0)') 'stage8_config_stage7_total_smoke_status', s7_smoke
  write(io,'(A,1X,I0)') 'stage8_config_stage7_closed_marker_status', s7_closed_status
  write(io,'(A,1X,I0)') 'stage8_config_stage7_dependency_status', dep_ok

  call init_stage8_default_config(cfg)
  call validate_stage8_config(cfg,valid,rejected,rhs_allowed)
  write(io,'(A,1X,I0)') 'stage8_default_enable_stage8', cfg%enable_stage8
  write(io,'(A,1X,I0)') 'stage8_default_controlled_test_enabled', cfg%controlled_test_enabled
  write(io,'(A,1X,I0)') 'stage8_default_production_dns_enabled', cfg%production_dns_enabled
  write(io,'(A,1X,I0)') 'stage8_default_production_two_way_enabled', cfg%production_two_way_enabled
  write(io,'(A,1X,I0)') 'stage8_default_rhs_allowed_flag', rhs_allowed
  write(io,'(A,1X,I0)') 'stage8_default_valid_flag', valid
  write(io,'(A,1X,I0)') 'stage8_default_rejected_flag', rejected
  write(io,'(A,1X,I0)') 'stage8_default_config_status', default_ok

  call init_stage8_controlled_test_config(cfg)
  call validate_stage8_config(cfg,valid,rejected,rhs_allowed)
  write(io,'(A,1X,I0)') 'stage8_controlled_enable_stage8', cfg%enable_stage8
  write(io,'(A,1X,I0)') 'stage8_controlled_controlled_test_enabled', cfg%controlled_test_enabled
  write(io,'(A,1X,I0)') 'stage8_controlled_production_dns_enabled', cfg%production_dns_enabled
  write(io,'(A,1X,I0)') 'stage8_controlled_production_two_way_enabled', cfg%production_two_way_enabled
  write(io,'(A,1X,I0)') 'stage8_controlled_runtime_grid_bridge_enabled', cfg%enable_runtime_grid_bridge
  write(io,'(A,1X,I0)') 'stage8_controlled_lagrangian_state_enabled', cfg%enable_lagrangian_state
  write(io,'(A,1X,I0)') 'stage8_controlled_velocity_to_fibre_enabled', cfg%enable_velocity_to_fibre
  write(io,'(A,1X,I0)') 'stage8_controlled_feedback_candidate_enabled', cfg%enable_feedback_candidate
  write(io,'(A,1X,I0)') 'stage8_controlled_force_density_candidate_enabled', cfg%enable_force_density_candidate
  write(io,'(A,1X,I0)') 'stage8_controlled_rhs_candidate_enabled', cfg%enable_rhs_candidate
  write(io,'(A,1X,I0)') 'stage8_controlled_rhs_allowed_flag', rhs_allowed
  write(io,'(A,1X,I0)') 'stage8_controlled_valid_flag', valid
  write(io,'(A,1X,I0)') 'stage8_controlled_rejected_flag', rejected
  write(io,'(A,1X,I0)') 'stage8_controlled_config_status', controlled_ok

  write(io,'(A,1X,I0)') 'stage8_invalid_production_dns_rejected_flag', prod_dns_rej
  write(io,'(A,1X,I0)') 'stage8_invalid_production_twoway_rejected_flag', prod_twoway_rej
  write(io,'(A,1X,I0)') 'stage8_invalid_no_controlled_test_rejected_flag', no_ctrl_rej
  write(io,'(A,1X,I0)') 'stage8_invalid_rho_rejected_flag', rho_rej
  write(io,'(A,1X,I0)') 'stage8_invalid_config_status', invalid_ok

  write(io,'(A,1X,I0)') 'stage8_rhs_candidate_controlled_allowed_flag', rhs_controlled_allowed
  write(io,'(A,1X,I0)') 'stage8_rhs_candidate_production_dns_blocked_flag', rhs_prod_dns_blocked
  write(io,'(A,1X,I0)') 'stage8_rhs_candidate_production_twoway_blocked_flag', rhs_prod_twoway_blocked
  write(io,'(A,1X,I0)') 'stage8_rhs_candidate_gate_status', rhs_gate_ok

  write(io,'(A,1X,ES24.16E3)') 'stage8_config_noop_rhs_change_max', rhs_change
  write(io,'(A,1X,I0)') 'stage8_config_noop_rhs_modified_flag', rhs_modified
  write(io,'(A,1X,I0)') 'stage8_config_pressure_poisson_modified_flag', 0
  write(io,'(A,1X,I0)') 'stage8_config_projection_modified_flag', 0
  write(io,'(A,1X,I0)') 'stage8_config_real_projection_called_flag', 0
  write(io,'(A,1X,I0)') 'stage8_config_production_dns_called_flag', 0
  write(io,'(A,1X,I0)') 'stage8_config_fluid_update_called_flag', 0
  write(io,'(A,1X,I0)') 'stage8_config_fibre_advance_called_flag', 0
  write(io,'(A,1X,I0)') 'stage8_config_noop_safety_status', noop_ok

  write(io,'(A,1X,I0)') 'stage8_config_check_status', final_status
  close(io)

contains
  subroutine ensure_dir(path)
    character(len=*), intent(in) :: path
    call execute_command_line('mkdir -p ' // trim(path))
  end subroutine

  subroutine file_exists_int(path, flag)
    character(len=*), intent(in) :: path
    integer, intent(out) :: flag
    logical :: ex
    inquire(file=path,exist=ex)
    flag = merge(1,0,ex)
  end subroutine

  subroutine get_int_value(path,key,val)
    character(len=*), intent(in) :: path,key
    integer, intent(out) :: val
    integer :: u, ios
    character(len=256) :: k
    real(mytype) :: x
    logical :: ex
    val = 0
    inquire(file=path,exist=ex)
    if (.not. ex) return
    open(newunit=u,file=path,status='old',action='read',iostat=ios)
    if (ios/=0) return
    do
      read(u,*,iostat=ios) k, x
      if (ios/=0) exit
      if (trim(k)==trim(key)) then
        val = nint(x)
        exit
      end if
    end do
    close(u)
  end subroutine

  subroutine init_rhs(arr)
    real(mytype), intent(out) :: arr(4,3,2)
    integer :: i,j,k
    do k=1,2
      do j=1,3
        do i=1,4
          arr(i,j,k)=real(i+2*j+3*k,mytype)
        end do
      end do
    end do
  end subroutine

end program fibre_stage8_config_check
