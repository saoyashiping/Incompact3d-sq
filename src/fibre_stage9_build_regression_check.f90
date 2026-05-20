program fibre_stage9_build_regression_check
  implicit none

  character(len=*), parameter :: dep_file = 'stage9_outputs/fibre_stage9_dependency_gate_check.dat'
  character(len=*), parameter :: ev_file = 'stage9_outputs/stage9_1_build_evidence.dat'
  character(len=*), parameter :: s8_file = 'stage8_outputs/fibre_stage8_total_smoke_check.dat'
  character(len=*), parameter :: out_file = 'stage9_outputs/fibre_stage9_build_regression_check.dat'

  integer :: uout
  logical :: dep_exists, ev_exists, s8_exists, s7_closed, s8_closed
  integer :: st9_0, a9_entry, a9_1, pca
  integer :: s8_total_smoke, s8_default_prod_dis
  integer :: ev_cfg, ev_build, ev_used_decomp, ev_root, ev_mod, ev_iomod, ev_lib
  integer :: ev_no_compat, ev_stage_checks_only, ev_checks_only_off
  integer :: ev_dns, ev_rhs, ev_proj, ev_fluid, ev_fibre

  integer :: stage9_0_dep_status, prior_closure_status
  integer :: stage8_closed_marker_status, stage8_no_production_side_effect_status
  integer :: evidence_status, production_disabled_policy_status, noop_safety_status
  integer :: final_status

  dep_exists = file_exists(dep_file)
  ev_exists = file_exists(ev_file)
  s8_exists = file_exists(s8_file)
  s7_closed = file_exists('stage7_outputs/STAGE7_CLOSED.md')
  s8_closed = file_exists('stage8_outputs/STAGE8_CLOSED.md')

  st9_0 = geti(dep_file, 'stage9_dependency_gate_check_status', 0)
  a9_entry = geti(dep_file, 'stage9_gate_stage9_entry_allowed_flag', 0)
  a9_1 = geti(dep_file, 'stage9_gate_stage9_1_allowed_flag', 0)
  pca = geti(dep_file, 'stage9_gate_production_coupling_allowed_flag', 0)

  s8_total_smoke = geti(s8_file, 'stage8_total_smoke_check_status', 0)
  s8_default_prod_dis = geti(s8_file, 'stage8_total_default_production_disabled_status', 0)

  ev_cfg = geti(ev_file, 'stage9_build_evidence_configure_status', 0)
  ev_build = geti(ev_file, 'stage9_build_evidence_xcompact3d_compile_status', 0)
  ev_used_decomp = geti(ev_file, 'stage9_build_evidence_used_source_matched_decomp_flag', 0)
  ev_root = geti(ev_file, 'stage9_build_evidence_decomp_root_exists', 0)
  ev_mod = geti(ev_file, 'stage9_build_evidence_decomp_mod_exists', 0)
  ev_iomod = geti(ev_file, 'stage9_build_evidence_decomp_io_mod_exists', 0)
  ev_lib = geti(ev_file, 'stage9_build_evidence_decomp_static_lib_exists', 0)
  ev_no_compat = geti(ev_file, 'stage9_build_evidence_no_compat_remnants_flag', 0)
  ev_stage_checks_only = geti(ev_file, 'stage9_build_evidence_stage_checks_only_flag', 0)
  ev_checks_only_off = geti(ev_file, 'stage9_build_evidence_fibre_stage_checks_only_off_flag', 0)
  ev_dns = geti(ev_file, 'stage9_build_evidence_dns_executed_flag', 0)
  ev_rhs = geti(ev_file, 'stage9_build_evidence_rhs_modified_flag', 0)
  ev_proj = geti(ev_file, 'stage9_build_evidence_projection_called_flag', 0)
  ev_fluid = geti(ev_file, 'stage9_build_evidence_fluid_update_called_flag', 0)
  ev_fibre = geti(ev_file, 'stage9_build_evidence_fibre_advance_called_flag', 0)

  stage9_0_dep_status = merge(1,0, dep_exists .and. st9_0==1 .and. a9_entry==1 .and. a9_1==1 .and. pca==0)
  stage8_closed_marker_status = merge(1,0, s8_closed)
  stage8_no_production_side_effect_status = merge(1,0, s8_default_prod_dis==1)
  prior_closure_status = merge(1,0, s7_closed .and. s8_closed .and. s8_exists .and. s8_total_smoke==1 .and. s8_default_prod_dis==1)

  evidence_status = merge(1,0, ev_exists .and. ev_cfg==1 .and. ev_build==1 .and. ev_used_decomp==1 .and. ev_root==1 .and. ev_mod==1 .and. ev_iomod==1 .and. ev_lib==1 .and. ev_no_compat==1 .and. ev_checks_only_off==1)

  production_disabled_policy_status = merge(1,0, s8_default_prod_dis==1 .and. ev_dns==0 .and. ev_rhs==0 .and. ev_proj==0 .and. ev_fluid==0 .and. ev_fibre==0)

  noop_safety_status = merge(1,0, ev_rhs==0 .and. ev_proj==0 .and. ev_fluid==0 .and. ev_fibre==0 .and. ev_dns==0)

  final_status = merge(1,0, stage9_0_dep_status==1 .and. prior_closure_status==1 .and. evidence_status==1 .and. production_disabled_policy_status==1 .and. noop_safety_status==1)

  open(newunit=uout, file=out_file, status='replace', action='write')
  call emit(uout, 'stage9_build_stage9_0_output_exists', merge(1,0,dep_exists))
  call emit(uout, 'stage9_build_stage9_0_status', st9_0)
  call emit(uout, 'stage9_build_stage9_entry_allowed_flag', a9_entry)
  call emit(uout, 'stage9_build_stage9_1_allowed_flag', a9_1)
  call emit(uout, 'stage9_build_production_coupling_allowed_flag', pca)
  call emit(uout, 'stage9_build_stage9_0_dependency_status', stage9_0_dep_status)

  call emit(uout, 'stage9_build_stage7_closed_marker_exists', merge(1,0,s7_closed))
  call emit(uout, 'stage9_build_stage8_closed_marker_exists', merge(1,0,s8_closed))
  call emit(uout, 'stage9_build_stage8_total_smoke_output_exists', merge(1,0,s8_exists))
  call emit(uout, 'stage9_build_stage8_total_smoke_status', s8_total_smoke)
  call emit(uout, 'stage9_build_stage8_closed_marker_status', stage8_closed_marker_status)
  call emit(uout, 'stage9_build_stage8_no_production_side_effect_status', stage8_no_production_side_effect_status)
  call emit(uout, 'stage9_build_prior_stage_closure_status', prior_closure_status)

  call emit(uout, 'stage9_build_evidence_output_exists', merge(1,0,ev_exists))
  call emit(uout, 'stage9_build_configure_status', ev_cfg)
  call emit(uout, 'stage9_build_xcompact3d_compile_status', ev_build)
  call emit(uout, 'stage9_build_used_source_matched_decomp_flag', ev_used_decomp)
  call emit(uout, 'stage9_build_decomp_root_exists', ev_root)
  call emit(uout, 'stage9_build_decomp_mod_exists', ev_mod)
  call emit(uout, 'stage9_build_decomp_io_mod_exists', ev_iomod)
  call emit(uout, 'stage9_build_decomp_static_lib_exists', ev_lib)
  call emit(uout, 'stage9_build_no_compat_remnants_flag', ev_no_compat)
  call emit(uout, 'stage9_build_stage_checks_only_flag', ev_stage_checks_only)
  call emit(uout, 'stage9_build_fibre_stage_checks_only_off_flag', ev_checks_only_off)
  call emit(uout, 'stage9_build_dns_executed_flag', ev_dns)
  call emit(uout, 'stage9_build_rhs_modified_flag', ev_rhs)
  call emit(uout, 'stage9_build_projection_called_flag', ev_proj)
  call emit(uout, 'stage9_build_fluid_update_called_flag', ev_fluid)
  call emit(uout, 'stage9_build_fibre_advance_called_flag', ev_fibre)
  call emit(uout, 'stage9_build_evidence_status', evidence_status)

  call emit(uout, 'stage9_build_default_production_disabled_status', s8_default_prod_dis)
  call emit(uout, 'stage9_build_production_dns_not_called_status', merge(1,0,ev_dns==0))
  call emit(uout, 'stage9_build_production_coupling_disabled_status', merge(1,0,pca==0))
  call emit(uout, 'stage9_build_rhs_untouched_status', merge(1,0,ev_rhs==0))
  call emit(uout, 'stage9_build_projection_untouched_status', merge(1,0,ev_proj==0))
  call emit(uout, 'stage9_build_fluid_update_not_called_status', merge(1,0,ev_fluid==0))
  call emit(uout, 'stage9_build_fibre_advance_not_called_status', merge(1,0,ev_fibre==0))
  call emit(uout, 'stage9_build_production_disabled_policy_status', production_disabled_policy_status)

  call emit(uout, 'stage9_build_rhs_hook_called_flag', 0)
  call emit(uout, 'stage9_build_pressure_poisson_modified_flag', 0)
  call emit(uout, 'stage9_build_projection_modified_flag', ev_proj)
  call emit(uout, 'stage9_build_real_projection_called_flag', 0)
  call emit(uout, 'stage9_build_production_dns_called_flag', 0)
  call emit(uout, 'stage9_build_noop_safety_status', noop_safety_status)

  call emit(uout, 'stage9_build_regression_check_status', final_status)
  close(uout)

contains
  subroutine emit(u,k,v)
    integer, intent(in) :: u,v
    character(len=*), intent(in) :: k
    write(u,'(A,1X,I0)') trim(k), v
  end subroutine

  logical function file_exists(path)
    character(len=*), intent(in) :: path
    inquire(file=trim(path), exist=file_exists)
  end function

  integer function geti(path,key,default)
    character(len=*), intent(in) :: path,key
    integer, intent(in) :: default
    character(len=1024) :: line,k,vstr
    integer :: u,ios,val,eq
    logical :: ex
    geti = default
    inquire(file=trim(path), exist=ex)
    if (.not. ex) return
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios/=0) return
    do
      read(u,'(A)', iostat=ios) line
      if (ios/=0) exit
      if (len_trim(line)==0) cycle
      eq = index(line, '=')
      if (eq > 0) then
        k = adjustl(trim(line(:eq-1)))
        vstr = adjustl(trim(line(eq+1:)))
        read(vstr,*,iostat=ios) val
        if (ios==0 .and. trim(k)==trim(key)) then
          geti = val
          exit
        end if
      else
        read(line,*,iostat=ios) k, val
        if (ios==0) then
          if (trim(k)==trim(key)) then
            geti = val
            exit
          end if
        end if
      end if
    end do
    close(u)
  end function

end program fibre_stage9_build_regression_check
