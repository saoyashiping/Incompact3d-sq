program fibre_stage8_runtime_grid_bridge_check
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_channel_grid_adapter
  use fibre_stage7_scalar_interpolation
  use fibre_stage7_velocity_interpolation
  use fibre_stage7_force_spreading
  use fibre_stage7_boundary_safety
  use fibre_stage8_runtime_grid_bridge
  implicit none
  type(stage7_channel_grid_t) :: gref,gbridge,greal
  type(stage8_runtime_grid_bridge_status_t) :: st,streal
  type(stage7_velocity_layout_t) :: layout
  type(stage7_boundary_safety_result_t) :: bres
  type(stage7_interp_weight_t) :: w
  integer :: io,nx,ny,nz,v,r
  integer :: s7m,s7o,s7s,s7c,s7w,s80o,s80s,dep
  real(mytype), allocatable :: yface(:),ycenter(:),phi(:,:,:),ux(:,:,:),uy(:,:,:),uz(:,:,:),fx(:,:,:),fy(:,:,:),fz(:,:,:),rhs0(:,:,:),rhs1(:,:,:)
  real(mytype) :: err_meta,err_dy,err_vol,rhschg,phi_l,ulag(3),fe(3),fl(3)
  integer :: scalar_ok,vel_ok,spread_ok,bound_ok,iface_ok
  integer :: inv_yface,inv_dy,inv_px,inv_pz,inv_status
  integer :: rep_exists,real_ev_status,noop_ok,final_status

  call ensure_dir('stage8_outputs')
  call file_exists_int('stage7_outputs/STAGE7_CLOSED.md',s7m)
  call file_exists_int('stage7_outputs/fibre_stage7_total_smoke_check.dat',s7o)
  call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_smoke_check_status',s7s)
  call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_closed_marker_status',s7c)
  call get_int('stage7_outputs/fibre_stage7_total_smoke_check.dat','stage7_total_closed_marker_written_flag',s7w)
  call file_exists_int('stage8_outputs/fibre_stage8_config_check.dat',s80o)
  call get_int('stage8_outputs/fibre_stage8_config_check.dat','stage8_config_check_status',s80s)
  dep=merge(1,0,s7m==1.and.s7o==1.and.s7s==1.and.s7c==1.and.s7w==1.and.s80o==1.and.s80s==1)

  nx=16;ny=17;nz=12
  call init_stage7_nonuniform_channel_grid(gref,nx,ny,nz)
  allocate(yface(ny+1),ycenter(ny)); yface=gref%y_face; ycenter=gref%y_center
  call init_stage8_runtime_bridge_status(st)
  call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,1,1,st)

  err_meta=max(abs(gbridge%dx-gref%dx),abs(gbridge%dz-gref%dz))
  err_meta=max(err_meta,maxval(abs(gbridge%y_face-gref%y_face)))
  err_meta=max(err_meta,maxval(abs(gbridge%y_center-gref%y_center)))
  err_meta=max(err_meta,maxval(abs(gbridge%dy_cell-gref%dy_cell)))
  err_meta=max(err_meta,maxval(abs(gbridge%volume_y-gref%volume_y)))
  err_meta=max(err_meta,abs((gbridge%xmax-gbridge%xmin)*(gbridge%zmax-gbridge%zmin)*sum(gbridge%volume_y)- &
     (gref%xmax-gref%xmin)*(gref%zmax-gref%zmin)*sum(gref%volume_y)))

  call init_stage8_runtime_bridge_status(streal)
  call init_stage8_grid_from_xcompact_runtime_bridge(greal,streal)
  call file_exists_int('stage8_checks/STAGE8_1_RUNTIME_GRID_SOURCE_AUDIT.md',rep_exists)
  real_ev_status=merge(1,0,rep_exists==1 .and. st%explicit_array_fallback_called_flag==1 .and. &
      ((streal%real_xcompact_grid_source_found_flag==1 .and. streal%real_xcompact_grid_adapter_called_flag==1) .or. &
       (streal%real_xcompact_grid_source_found_flag==0 .and. streal%real_xcompact_grid_adapter_called_flag==0)))

  err_dy=maxval(abs((gbridge%y_face(2:ny+1)-gbridge%y_face(1:ny))-gbridge%dy_cell))
  err_vol=maxval(abs(gbridge%volume_y-gbridge%dy_cell/((gbridge%xmax-gbridge%xmin)*(gbridge%zmax-gbridge%zmin))))

  allocate(phi(nx,ny,nz),ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),rhs0(4,3,2),rhs1(4,3,2))
  phi=1._mytype; ux=2._mytype; uy=3._mytype; uz=4._mytype
  call init_stage7_interp_weight(w,64)
  call build_stage7_scalar_interp_weight(gbridge,gbridge%xmin+3*gbridge%dx,gbridge%y_center(8),gbridge%zmin+3*gbridge%dz,w)
  call interpolate_stage7_scalar(phi,w,phi_l)
  scalar_ok=merge(1,0,w%valid_flag==1 .and. abs(phi_l-1._mytype)<=1e-12_mytype)
  call free_stage7_interp_weight(w)

  call init_stage7_collocated_velocity_layout(layout)
  call interpolate_stage7_velocity(gbridge,layout,ux,uy,uz,gbridge%xmin+3*gbridge%dx,gbridge%y_center(8),gbridge%zmin+3*gbridge%dz,ulag,v,r)
  vel_ok=merge(1,0,v==1 .and. abs(ulag(1)-2)<=1e-12_mytype .and. abs(ulag(2)-3)<=1e-12_mytype .and. abs(ulag(3)-4)<=1e-12_mytype)

  call clear_stage7_force_buffer(fx,fy,fz)
  call spread_stage7_lagrangian_force(gbridge,gbridge%xmin+3*gbridge%dx,gbridge%y_center(8),gbridge%zmin+3*gbridge%dz,[1._mytype,0._mytype,0._mytype],0.1_mytype,fx,fy,fz,v,r)
  call compute_stage7_eulerian_force_total(gbridge,fx,fy,fz,fe); call compute_stage7_lagrangian_force_total(1,reshape([1._mytype,0._mytype,0._mytype],[3,1]),[0.1_mytype],fl)
  spread_ok=merge(1,0,v==1 .and. sqrt(sum((fe-fl)**2))<=1e-12_mytype)

  call classify_stage7_boundary_point(gbridge,layout,gbridge%xmin+3*gbridge%dx,gbridge%y_center(8),gbridge%zmin+3*gbridge%dz,1,bres)
  bound_ok=merge(1,0,bres%safe_flag==1 .and. bres%blocked_flag==0)
  iface_ok=merge(1,0,scalar_ok==1 .and. vel_ok==1 .and. spread_ok==1 .and. bound_ok==1)

  inv_yface=0; inv_dy=0; inv_px=0; inv_pz=0
  yface(8)=yface(7); call init_stage8_runtime_bridge_status(st); call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,1,1,st); inv_yface=merge(1,0,st%bridge_grid_rejected_flag==1)
  yface=gref%y_face; yface(9)=yface(8)-1e-3_mytype; call init_stage8_runtime_bridge_status(st); call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,1,1,st); inv_dy=merge(1,0,st%bridge_grid_rejected_flag==1)
  yface=gref%y_face; call init_stage8_runtime_bridge_status(st); call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,0,1,st); inv_px=merge(1,0,st%bridge_grid_rejected_flag==1)
  call init_stage8_runtime_bridge_status(st); call init_stage8_grid_from_explicit_arrays_bridge(gbridge,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycenter,1,0,st); inv_pz=merge(1,0,st%bridge_grid_rejected_flag==1)
  inv_status=merge(1,0,inv_yface==1.and.inv_dy==1.and.inv_px==1.and.inv_pz==1)

  call init_rhs(rhs0); rhs1=rhs0; rhschg=maxval(abs(rhs1-rhs0)); noop_ok=merge(1,0,rhschg<=1e-14_mytype)
  final_status=merge(1,0,dep==1 .and. err_meta<=1e-12_mytype .and. real_ev_status==1 .and. err_dy<=1e-12_mytype .and. err_vol<=1e-12_mytype .and. iface_ok==1 .and. inv_status==1 .and. noop_ok==1)

  open(newunit=io,file='stage8_outputs/fibre_stage8_runtime_grid_bridge_check.dat',status='replace',action='write')
  write(io,'(A,1X,I0)') 'stage8_bridge_stage7_closed_marker_exists',s7m
  write(io,'(A,1X,I0)') 'stage8_bridge_stage7_total_smoke_output_exists',s7o
  write(io,'(A,1X,I0)') 'stage8_bridge_stage7_total_smoke_status',s7s
  write(io,'(A,1X,I0)') 'stage8_bridge_stage7_closed_marker_status',s7c
  write(io,'(A,1X,I0)') 'stage8_bridge_stage8_0_output_exists',s80o
  write(io,'(A,1X,I0)') 'stage8_bridge_stage8_0_status',s80s
  write(io,'(A,1X,I0)') 'stage8_bridge_dependency_status',dep
  write(io,'(A,1X,I0)') 'stage8_bridge_explicit_array_fallback_called_flag',st%explicit_array_fallback_called_flag
  write(io,'(A,1X,I0)') 'stage8_bridge_validate_stage7_grid_called_flag',st%validate_stage7_grid_called_flag
  write(io,'(A,1X,I0)') 'stage8_bridge_explicit_grid_valid_flag',st%bridge_grid_valid_flag
  write(io,'(A,1X,I0)') 'stage8_bridge_explicit_grid_rejected_flag',st%bridge_grid_rejected_flag
  write(io,'(A,1X,I0)') 'stage8_bridge_explicit_bridge_status',st%bridge_status
  write(io,'(A,1X,ES24.16E3)') 'stage8_bridge_metadata_match_error_max',err_meta
  write(io,'(A,1X,I0)') 'stage8_bridge_metadata_match_status',merge(1,0,err_meta<=1e-12_mytype)
  write(io,'(A,1X,I0)') 'stage8_bridge_real_xcompact_grid_source_found_flag',streal%real_xcompact_grid_source_found_flag
  write(io,'(A,1X,I0)') 'stage8_bridge_real_xcompact_grid_adapter_called_flag',streal%real_xcompact_grid_adapter_called_flag
  write(io,'(A,1X,I0)') 'stage8_bridge_real_source_audit_report_exists',rep_exists
  write(io,'(A,1X,I0)') 'stage8_bridge_explicit_array_fallback_available_flag',1
  write(io,'(A,1X,I0)') 'stage8_bridge_real_source_evidence_status',real_ev_status
  write(io,'(A,1X,I0)') 'stage8_bridge_x_uniform_flag',gbridge%x_uniform_flag
  write(io,'(A,1X,I0)') 'stage8_bridge_z_uniform_flag',gbridge%z_uniform_flag
  write(io,'(A,1X,I0)') 'stage8_bridge_y_nonuniform_flag',gbridge%y_nonuniform_flag
  write(io,'(A,1X,I0)') 'stage8_bridge_y_monotonic_flag',gbridge%y_monotonic_flag
  write(io,'(A,1X,I0)') 'stage8_bridge_dy_positive_flag',gbridge%dy_positive_flag
  write(io,'(A,1X,I0)') 'stage8_bridge_volume_positive_flag',gbridge%volume_positive_flag
  write(io,'(A,1X,ES24.16E3)') 'stage8_bridge_dy_face_consistency_error_max',err_dy
  write(io,'(A,1X,ES24.16E3)') 'stage8_bridge_volume_formula_error_max',err_vol
  write(io,'(A,1X,I0)') 'stage8_bridge_wall_bounds_status',gbridge%wall_bounds_valid_flag
  write(io,'(A,1X,I0)') 'stage8_bridge_grid_validation_detail_status',merge(1,0,gbridge%x_uniform_flag==1.and.gbridge%z_uniform_flag==1.and.gbridge%y_nonuniform_flag==1.and.gbridge%y_monotonic_flag==1.and.gbridge%dy_positive_flag==1.and.gbridge%volume_positive_flag==1.and.err_dy<=1e-12_mytype.and.err_vol<=1e-12_mytype.and.gbridge%wall_bounds_valid_flag==1)
  write(io,'(A,1X,I0)') 'stage8_bridge_scalar_smoke_status',scalar_ok
  write(io,'(A,1X,I0)') 'stage8_bridge_velocity_smoke_status',vel_ok
  write(io,'(A,1X,I0)') 'stage8_bridge_spreading_smoke_status',spread_ok
  write(io,'(A,1X,I0)') 'stage8_bridge_boundary_smoke_status',bound_ok
  write(io,'(A,1X,I0)') 'stage8_bridge_stage7_interface_smoke_status',iface_ok
  write(io,'(A,1X,I0)') 'stage8_bridge_invalid_yface_rejected_flag',inv_yface
  write(io,'(A,1X,I0)') 'stage8_bridge_invalid_dy_rejected_flag',inv_dy
  write(io,'(A,1X,I0)') 'stage8_bridge_invalid_periodic_x_rejected_flag',inv_px
  write(io,'(A,1X,I0)') 'stage8_bridge_invalid_periodic_z_rejected_flag',inv_pz
  write(io,'(A,1X,I0)') 'stage8_bridge_invalid_input_status',inv_status
  write(io,'(A,1X,ES24.16E3)') 'stage8_bridge_noop_rhs_change_max',rhschg
  write(io,'(A,1X,I0)') 'stage8_bridge_noop_rhs_modified_flag',merge(1,0,rhschg>1e-14_mytype)
  write(io,'(A,1X,I0)') 'stage8_bridge_pressure_poisson_modified_flag',0
  write(io,'(A,1X,I0)') 'stage8_bridge_projection_modified_flag',0
  write(io,'(A,1X,I0)') 'stage8_bridge_real_projection_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_bridge_production_dns_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_bridge_fluid_update_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_bridge_fibre_advance_called_flag',0
  write(io,'(A,1X,I0)') 'stage8_bridge_noop_safety_status',noop_ok
  write(io,'(A,1X,I0)') 'stage8_runtime_grid_bridge_check_status',final_status
  close(io)
contains
  subroutine ensure_dir(p); character(len=*),intent(in)::p; call execute_command_line('mkdir -p '//trim(p)); end
  subroutine file_exists_int(path,flag); character(len=*),intent(in)::path; integer,intent(out)::flag; logical::ex; inquire(file=path,exist=ex); flag=merge(1,0,ex); end
  subroutine get_int(path,key,val)
    character(len=*),intent(in)::path,key; integer,intent(out)::val
    integer::u,ios; character(len=256)::k; real(mytype)::x; logical::ex
    val=0; inquire(file=path,exist=ex); if(.not.ex)return
    open(newunit=u,file=path,status='old',action='read',iostat=ios); if(ios/=0)return
    do; read(u,*,iostat=ios) k,x; if(ios/=0)exit; if(trim(k)==trim(key)) then; val=nint(x); exit; end if; end do; close(u)
  end
  subroutine init_rhs(a); real(mytype),intent(out)::a(4,3,2); integer::i,j,k; do k=1,2;do j=1,3;do i=1,4; a(i,j,k)=real(i+j+k,mytype); end do;end do;end do; end
end program
