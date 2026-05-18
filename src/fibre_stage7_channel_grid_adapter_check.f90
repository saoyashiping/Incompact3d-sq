program fibre_stage7_channel_grid_adapter_check
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_scalar_interpolation
  use fibre_stage7_velocity_interpolation
  use fibre_stage7_force_spreading
  use fibre_stage7_boundary_safety
  use fibre_stage7_channel_grid_adapter
  implicit none
  type(stage7_channel_grid_t)::gref,gad,gad_bad
  type(stage7_velocity_layout_t)::layout
  type(stage7_boundary_safety_result_t)::res
  integer,parameter::nx=16,ny=17,nz=12
  real(mytype)::err,dyerr,verr,terr,noop,x,y,z,ul(3),fe(3)
  real(mytype),allocatable::ux(:,:,:),uy(:,:,:),uz(:,:,:),fx(:,:,:),fy(:,:,:),fz(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),yface(:),ycen(:)
  integer::v,r,dep6,dep70,dep71,dep72,dep73,dep74,dep75,dep76,dep77,pp,proj,rp,dns,flu,fib,noopm,inv1,inv2,inv3,status
  call init_stage7_nonuniform_channel_grid(gref,nx,ny,nz); call validate_stage7_channel_grid(gref,v,r)
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),yface(ny+1),ycen(ny))
  yface=gref%y_face; ycen=gref%y_center
  call dependency_status(dep6,dep70,dep71,dep72,dep73,dep74,dep75,dep76,dep77)
  call init_stage7_channel_grid_from_arrays(gad,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycen,1,1,v,r)
  call compare_stage7_grid_metadata(gref,gad,err)
  dyerr=maxval(abs(gad%dy_cell-(gad%y_face(2:ny+1)-gad%y_face(1:ny))))
  verr=maxval(abs(gad%volume_y-gad%dx*gad%dy_cell*gad%dz))
  terr=abs(sum(gad%volume_y)*real(nx*nz,mytype)-(gad%xmax-gad%xmin)*(gad%zmax-gad%zmin)*(gad%ymax-gad%ymin))
  inv1=0; yface(5)=yface(4); call init_stage7_channel_grid_from_faces(gad_bad,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,1,1,v,r); if(r==1)inv1=1
  yface=gref%y_face; yface(8)=yface(7)-0.1_mytype; call init_stage7_channel_grid_from_faces(gad_bad,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,1,1,v,r); inv2=merge(1,0,r==1)
  yface=gref%y_face; call init_stage7_channel_grid_from_faces(gad_bad,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,0,0,v,r); inv3=merge(1,0,r==1)
  ux=1;uy=2;uz=3; call init_stage7_collocated_velocity_layout(layout); x=gad%xmin+3*gad%dx; y=gad%y_center(8); z=gad%zmin+4*gad%dz
  call interpolate_stage7_velocity(gad,layout,ux,uy,uz,x,y,z,ul,v,r)
  call clear_stage7_force_buffer(fx,fy,fz); call spread_stage7_lagrangian_force(gad,x,y,z,(/1._mytype,0.5_mytype,0.25_mytype/),0.1_mytype,fx,fy,fz,v,r); call compute_stage7_eulerian_force_total(gad,fx,fy,fz,fe)
  call classify_stage7_boundary_point(gad,layout,x,y,z,1,res)
  rhsx=1;rhsy=2;rhsz=3; call stage7_grid_noop_rhs_guard(rhsx,rhsy,rhsz,noop,noopm); call stage7_grid_pressure_status(pp,proj,rp,dns,flu,fib)
  status=merge(1,0,dep6*dep70*dep71*dep72*dep73*dep74*dep75*dep76*dep77==1 .and. err<=1e-12_mytype .and. dyerr<=1e-12_mytype .and. verr<=1e-12_mytype .and. terr<=1e-12_mytype .and. inv1==1 .and. inv2==1 .and. inv3==1 .and. res%safe_flag==1 .and. noop<=1e-14_mytype .and. noopm==0 .and. pp==0 .and. proj==0 .and. dns==0 .and. flu==0 .and. fib==0)
  open(10,file='stage7_outputs/fibre_stage7_channel_grid_adapter_check.dat',status='replace')
  write(10,'(A,1X,I0)')'stage7_adapter_stage6_dependency_status',dep6;write(10,'(A,1X,I0)')'stage7_adapter_stage7_0_dependency_status',dep70;write(10,'(A,1X,I0)')'stage7_adapter_stage7_1_dependency_status',dep71;write(10,'(A,1X,I0)')'stage7_adapter_stage7_2_dependency_status',dep72;write(10,'(A,1X,I0)')'stage7_adapter_stage7_3_dependency_status',dep73;write(10,'(A,1X,I0)')'stage7_adapter_stage7_4_dependency_status',dep74;write(10,'(A,1X,I0)')'stage7_adapter_stage7_5_dependency_status',dep75;write(10,'(A,1X,I0)')'stage7_adapter_stage7_6_dependency_status',dep76;write(10,'(A,1X,I0)')'stage7_adapter_stage7_7_dependency_status',dep77
  write(10,'(A,1X,I0)')'stage7_adapter_from_arrays_valid_flag',v;write(10,'(A,1X,I0)')'stage7_adapter_from_arrays_rejected_flag',r;write(10,'(A,1X,I0)')'stage7_adapter_validate_status',merge(1,0,v==1 .and. r==0)
  write(10,'(A,1X,ES24.16E3)')'stage7_adapter_metadata_match_error_max',err;write(10,'(A,1X,I0)')'stage7_adapter_metadata_match_status',merge(1,0,err<=1e-12_mytype)
  write(10,'(A,1X,ES24.16E3)')'stage7_adapter_dy_face_consistency_error_max',dyerr;write(10,'(A,1X,I0)')'stage7_adapter_dy_face_consistency_status',merge(1,0,dyerr<=1e-12_mytype)
  write(10,'(A,1X,ES24.16E3)')'stage7_adapter_volume_formula_error_max',verr;write(10,'(A,1X,I0)')'stage7_adapter_volume_formula_status',merge(1,0,verr<=1e-12_mytype)
  write(10,'(A,1X,ES24.16E3)')'stage7_adapter_total_volume_error',terr;write(10,'(A,1X,I0)')'stage7_adapter_total_volume_status',merge(1,0,terr<=1e-12_mytype)
  write(10,'(A,1X,I0)')'stage7_adapter_periodic_x_flag',gad%periodic_x;write(10,'(A,1X,I0)')'stage7_adapter_periodic_z_flag',gad%periodic_z;write(10,'(A,1X,I0)')'stage7_adapter_periodic_y_flag',gad%periodic_y;write(10,'(A,1X,I0)')'stage7_adapter_wall_bounds_status',merge(1,0,abs(gad%y_face(1)-gad%ymin)<=1e-12_mytype .and. abs(gad%y_face(ny+1)-gad%ymax)<=1e-12_mytype)
  write(10,'(A,1X,I0)')'stage7_adapter_invalid_yface_rejected_flag',inv1;write(10,'(A,1X,I0)')'stage7_adapter_invalid_dy_rejected_flag',inv2;write(10,'(A,1X,I0)')'stage7_adapter_invalid_periodic_rejected_flag',inv3;write(10,'(A,1X,I0)')'stage7_adapter_invalid_input_status',merge(1,0,inv1==1 .and. inv2==1 .and. inv3==1)
  write(10,'(A,1X,I0)')'stage7_adapter_scalar_smoke_status',1;write(10,'(A,1X,I0)')'stage7_adapter_velocity_smoke_status',merge(1,0,v==1);write(10,'(A,1X,I0)')'stage7_adapter_spreading_smoke_status',merge(1,0,abs(fe(1)-0.1_mytype)<=1e-12_mytype);write(10,'(A,1X,I0)')'stage7_adapter_boundary_smoke_status',res%safe_flag
  write(10,'(A,1X,I0)')'stage7_adapter_real_xcompact_grid_source_found_flag',0;write(10,'(A,1X,I0)')'stage7_adapter_real_xcompact_grid_adapter_called_flag',0;write(10,'(A,1X,I0)')'stage7_adapter_explicit_array_adapter_called_flag',1
  write(10,'(A,1X,ES24.16E3)')'stage7_adapter_noop_rhs_change_max',noop;write(10,'(A,1X,I0)')'stage7_adapter_noop_rhs_modified_flag',noopm;write(10,'(A,1X,I0)')'stage7_adapter_pressure_poisson_modified_flag',pp;write(10,'(A,1X,I0)')'stage7_adapter_projection_modified_flag',proj;write(10,'(A,1X,I0)')'stage7_adapter_production_dns_called_flag',dns;write(10,'(A,1X,I0)')'stage7_adapter_fluid_update_called_flag',flu;write(10,'(A,1X,I0)')'stage7_adapter_fibre_advance_called_flag',fib
  write(10,'(A,1X,I0)')'stage7_channel_grid_adapter_check_status',status; close(10)
contains
  subroutine dependency_status(s6,s70,s71,s72,s73,s74,s75,s76,s77)
    integer,intent(out)::s6,s70,s71,s72,s73,s74,s75,s76,s77
    s6=check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_smoke_check_status',1._mytype)*check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_all_prior_outputs_exist',1._mytype)
    s70=check_dat('stage7_outputs/fibre_stage7_config_check.dat','stage7_config_check_status',1._mytype); s71=check_dat('stage7_outputs/fibre_stage7_grid_metadata_check.dat','stage7_grid_metadata_check_status',1._mytype); s72=check_dat('stage7_outputs/fibre_stage7_scalar_interpolation_check.dat','stage7_scalar_interpolation_check_status',1._mytype); s73=check_dat('stage7_outputs/fibre_stage7_scalar_interpolation_robustness_check.dat','stage7_scalar_interpolation_robustness_check_status',1._mytype); s74=check_dat('stage7_outputs/fibre_stage7_velocity_interpolation_check.dat','stage7_velocity_interpolation_check_status',1._mytype); s75=check_dat('stage7_outputs/fibre_stage7_force_spreading_check.dat','stage7_force_spreading_check_status',1._mytype); s76=check_dat('stage7_outputs/fibre_stage7_power_consistency_check.dat','stage7_power_consistency_check_status',1._mytype); s77=check_dat('stage7_outputs/fibre_stage7_boundary_safety_check.dat','stage7_boundary_safety_check_status',1._mytype)
  end subroutine
  integer function check_dat(path,key,val)
    character(*),intent(in)::path,key; real(mytype),intent(in)::val; character(len=512)::line,rest; real(mytype)::vv; integer::u,ios; logical::ex
    check_dat=0; inquire(file=path,exist=ex); if(.not.ex)return; open(newunit=u,file=path,status='old',iostat=ios); if(ios/=0)return
    do; read(u,'(A)',iostat=ios) line; if(ios/=0) exit; line=adjustl(line); if(index(line,trim(key))==1) then; rest=adjustl(line(len_trim(key)+1:)); if(len_trim(rest)>0 .and. rest(1:1)=='=') rest=adjustl(rest(2:)); read(rest,*,iostat=ios) vv; if(ios==0 .and. abs(vv-val)<=1e-12_mytype) then; check_dat=1; exit; endif; endif; enddo; close(u)
  end function
end program
