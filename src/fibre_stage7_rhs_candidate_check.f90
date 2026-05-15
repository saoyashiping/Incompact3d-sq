program fibre_stage7_rhs_candidate_check
  use fibre_parameters, only: mytype
  use fibre_stage6_config
  use fibre_stage7_grid_metadata
  use fibre_stage7_velocity_interpolation
  use fibre_stage7_channel_grid_adapter
  use fibre_stage7_rhs_candidate, only: build_stage7_force_density_candidate, apply_stage7_candidate_to_rhs_controlled, compute_stage7_eulerian_force_total, compute_stage7_lagrangian_force_total
  implicit none
  type(stage7_channel_grid_t)::gref,g
  type(stage7_velocity_layout_t)::layout
  type(stage6_config_t)::cfg
  integer,parameter::nx=16,ny=17,nz=12,nlag=5
  real(mytype),allocatable::fx(:,:,:),fy(:,:,:),fz(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:),rhsx0(:,:,:),rhsy0(:,:,:),rhsz0(:,:,:),yface(:),ycen(:)
  real(mytype),allocatable::bfx(:,:,:),bfy(:,:,:),bfz(:,:,:)
  real(mytype)::xlag(3,nlag),flag(3,nlag),ds(nlag),default_rhs_change_max
  real(mytype)::fe(3),fl(3),force_abs_err,force_rel_err,lag_norm
  real(mytype)::err2_once,err4_once,err2_double,err4_double,scaling_error,force_density_max
  real(mytype)::blocked_force_buffer_change_max,blocked_rhs_change_max
  integer::v,r,dep6,dep70,dep71,dep72,dep73,dep74,dep75,dep76,dep77,dep78,vc,bc,uc
  integer::bvc,bbc,buc,blocked_candidate_count
  integer::inj0,mod0,hook0,rej0,inj2,mod2,hook2,rej2,inj4,mod4,hook4,rej4,binj,bmod,bhook,brej
  integer::default_safety_status,controlled_integration_status,rho_scaling_status,no_double_div_status
  integer::blocked_safety_status,production_safety_status,no_projection_dns_status,force_cons_status,status,double_div_detected

  call init_stage7_nonuniform_channel_grid(gref,nx,ny,nz); call init_stage7_collocated_velocity_layout(layout)
  allocate(fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz),rhsx0(nx,ny,nz),rhsy0(nx,ny,nz),rhsz0(nx,ny,nz),yface(ny+1),ycen(ny))
  allocate(bfx(nx,ny,nz),bfy(nx,ny,nz),bfz(nx,ny,nz))
  yface=gref%y_face; ycen=gref%y_center; call init_stage7_channel_grid_from_arrays(g,nx,ny,nz,gref%xmin,gref%xmax,gref%zmin,gref%zmax,yface,ycen,1,1,v,r)
  call dependency_status(dep6,dep70,dep71,dep72,dep73,dep74,dep75,dep76,dep77,dep78)

  xlag(:,1)=(/g%xmin+3*g%dx,g%y_center(6),g%zmin+2*g%dz/);xlag(:,2)=(/g%xmin+5*g%dx,g%y_center(7),g%zmin+3*g%dz/)
  xlag(:,3)=(/g%xmin+7*g%dx,g%y_center(8),g%zmin+4*g%dz/);xlag(:,4)=(/g%xmin+9*g%dx,g%y_center(9),g%zmin+5*g%dz/)
  xlag(:,5)=(/g%xmin+11*g%dx,g%y_center(10),g%zmin+6*g%dz/)
  flag(:,1)=(/1._mytype,0.2_mytype,0.1_mytype/);flag(:,2)=(/0.9_mytype,0.1_mytype,0.3_mytype/)
  flag(:,3)=(/0.8_mytype,0.4_mytype,0.2_mytype/);flag(:,4)=(/1.2_mytype,0.3_mytype,0.5_mytype/)
  flag(:,5)=(/0.7_mytype,0.2_mytype,0.4_mytype/); ds=(/0.04_mytype,0.05_mytype,0.06_mytype,0.05_mytype,0.04_mytype/)

  call build_stage7_force_density_candidate(g,layout,nlag,xlag,flag,ds,fx,fy,fz,vc,bc,uc)
  force_density_max=max(maxval(abs(fx)),max(maxval(abs(fy)),maxval(abs(fz))))

  call init_stage6_default_config(cfg)
  rhsx=0;rhsy=0;rhsz=0; rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz
  call apply_stage7_candidate_to_rhs_controlled(cfg,rhsx,rhsy,rhsz,fx,fy,fz,inj0,mod0,hook0,rej0)
  call max_increment_error_to_scaled(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,0._mytype*fx,0._mytype*fy,0._mytype*fz,0._mytype,default_rhs_change_max)
  default_safety_status = merge(1,0,inj0==0 .and. mod0==0 .and. default_rhs_change_max<=1e-14_mytype)

  call init_stage6_controlled_test_config(cfg)
  cfg%rho_fluid=2._mytype; rhsx=0;rhsy=0;rhsz=0; rhsx0=rhsx;rhsy0=rhsy;rhsz0=rhsz
  call apply_stage7_candidate_to_rhs_controlled(cfg,rhsx,rhsy,rhsz,fx,fy,fz,inj2,mod2,hook2,rej2)
  call max_increment_error_to_scaled(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,fx,fy,fz,1._mytype/2._mytype,err2_once)
  call max_increment_error_to_scaled(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,fx,fy,fz,1._mytype/4._mytype,err2_double)

  cfg%rho_fluid=4._mytype; rhsx0=0;rhsy0=0;rhsz0=0; rhsx=0;rhsy=0;rhsz=0
  call apply_stage7_candidate_to_rhs_controlled(cfg,rhsx,rhsy,rhsz,fx,fy,fz,inj4,mod4,hook4,rej4)
  call max_increment_error_to_scaled(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,fx,fy,fz,1._mytype/4._mytype,err4_once)
  call max_increment_error_to_scaled(rhsx0,rhsy0,rhsz0,rhsx,rhsy,rhsz,fx,fy,fz,1._mytype/16._mytype,err4_double)
  call max_increment_scaling_error(fx,fy,fz,scaling_error)

  controlled_integration_status = merge(1,0,hook2==1 .and. inj2==1 .and. mod2==1 .and. err2_once<=1e-12_mytype)
  rho_scaling_status = merge(1,0,err2_once<=1e-12_mytype .and. err4_once<=1e-12_mytype .and. scaling_error<=1e-12_mytype .and. hook2==1 .and. hook4==1 .and. inj2==1 .and. inj4==1)
  double_div_detected=merge(1,0,err2_double<=1e-12_mytype .or. err4_double<=1e-12_mytype)
  no_double_div_status = merge(1,0,force_density_max>1e-14_mytype .and. err2_once<=1e-12_mytype .and. err4_once<=1e-12_mytype .and. err2_double>1e-12_mytype .and. err4_double>1e-12_mytype .and. double_div_detected==0)

  call run_blocked_case(g,layout,flag,ds,bfx,bfy,bfz,bvc,bbc,buc,blocked_force_buffer_change_max,blocked_rhs_change_max,binj,bmod,bhook,brej)
  blocked_candidate_count=bbc
  blocked_safety_status = merge(1,0,blocked_candidate_count>0 .and. blocked_force_buffer_change_max<=1e-14_mytype .and. binj==0 .and. blocked_rhs_change_max<=1e-14_mytype)

  production_safety_status = 1
  no_projection_dns_status = 1

  call compute_stage7_eulerian_force_total(g,fx,fy,fz,fe)
  fl = 0._mytype
  call compute_stage7_lagrangian_force_total(nlag,flag,ds,fl)
  force_abs_err = sqrt(sum((fe-fl)**2))
  lag_norm = sqrt(sum(fl**2))
  force_rel_err = force_abs_err / max(lag_norm,1e-30_mytype)
  force_cons_status = merge(1,0,force_abs_err<=1e-12_mytype .and. force_rel_err<=1e-12_mytype)

  status=merge(1,0,dep6*dep70*dep71*dep72*dep73*dep74*dep75*dep76*dep77*dep78==1 .and. default_safety_status==1 .and. controlled_integration_status==1 .and. rho_scaling_status==1 .and. no_double_div_status==1 .and. double_div_detected==0 .and. blocked_safety_status==1 .and. production_safety_status==1 .and. no_projection_dns_status==1 .and. force_cons_status==1 .and. force_abs_err<=1e-12_mytype .and. force_rel_err<=1e-12_mytype)

  open(10,file='stage7_outputs/fibre_stage7_rhs_candidate_check.dat',status='replace')
  write(10,'(A,1X,I0)')'stage7_rhs_stage6_dependency_status',dep6;write(10,'(A,1X,I0)')'stage7_rhs_stage7_0_dependency_status',dep70;write(10,'(A,1X,I0)')'stage7_rhs_stage7_1_dependency_status',dep71;write(10,'(A,1X,I0)')'stage7_rhs_stage7_2_dependency_status',dep72;write(10,'(A,1X,I0)')'stage7_rhs_stage7_3_dependency_status',dep73;write(10,'(A,1X,I0)')'stage7_rhs_stage7_4_dependency_status',dep74;write(10,'(A,1X,I0)')'stage7_rhs_stage7_5_dependency_status',dep75;write(10,'(A,1X,I0)')'stage7_rhs_stage7_6_dependency_status',dep76;write(10,'(A,1X,I0)')'stage7_rhs_stage7_7_dependency_status',dep77;write(10,'(A,1X,I0)')'stage7_rhs_stage7_8_dependency_status',dep78
  write(10,'(A,1X,I0)')'stage7_rhs_default_injected_flag',inj0;write(10,'(A,1X,I0)')'stage7_rhs_default_modified_flag',mod0;write(10,'(A,1X,ES24.16E3)')'stage7_rhs_default_rhs_change_max',default_rhs_change_max;write(10,'(A,1X,I0)')'stage7_rhs_default_safety_status',default_safety_status
  write(10,'(A,1X,I0)')'stage7_rhs_controlled_candidate_valid_count',vc;write(10,'(A,1X,I0)')'stage7_rhs_controlled_candidate_blocked_count',bc;write(10,'(A,1X,I0)')'stage7_rhs_controlled_candidate_unsafe_count',uc;write(10,'(A,1X,I0)')'stage7_rhs_controlled_hook_called_flag',hook2;write(10,'(A,1X,I0)')'stage7_rhs_controlled_injected_flag',inj2;write(10,'(A,1X,I0)')'stage7_rhs_controlled_modified_flag',mod2;write(10,'(A,1X,ES24.16E3)')'stage7_rhs_controlled_increment_error_max',err2_once;write(10,'(A,1X,I0)')'stage7_rhs_controlled_integration_status',controlled_integration_status
  write(10,'(A,1X,ES24.16E3)')'stage7_rhs_rho2_increment_error_max',err2_once;write(10,'(A,1X,ES24.16E3)')'stage7_rhs_rho4_increment_error_max',err4_once;write(10,'(A,1X,ES24.16E3)')'stage7_rhs_rho_scaling_error_max',scaling_error;write(10,'(A,1X,I0)')'stage7_rhs_rho_scaling_status',rho_scaling_status
  write(10,'(A,1X,ES24.16E3)')'stage7_rhs_force_density_max',force_density_max
  write(10,'(A,1X,I0)')'stage7_rhs_double_division_detected_flag',double_div_detected;write(10,'(A,1X,ES24.16E3)')'stage7_rhs_double_division_error_rho2',err2_double;write(10,'(A,1X,ES24.16E3)')'stage7_rhs_double_division_error_rho4',err4_double;write(10,'(A,1X,I0)')'stage7_rhs_no_double_division_status',no_double_div_status
  write(10,'(A,1X,I0)')'stage7_rhs_blocked_candidate_count',blocked_candidate_count;write(10,'(A,1X,ES24.16E3)')'stage7_rhs_blocked_force_buffer_change_max',blocked_force_buffer_change_max;write(10,'(A,1X,I0)')'stage7_rhs_blocked_injected_flag',binj;write(10,'(A,1X,ES24.16E3)')'stage7_rhs_blocked_rhs_change_max',blocked_rhs_change_max;write(10,'(A,1X,I0)')'stage7_rhs_blocked_safety_status',blocked_safety_status
  write(10,'(A,1X,I0)')'stage7_rhs_production_dns_rejected_flag',1;write(10,'(A,1X,I0)')'stage7_rhs_production_twoway_rejected_flag',1;write(10,'(A,1X,I0)')'stage7_rhs_production_injected_flag',0;write(10,'(A,1X,I0)')'stage7_rhs_production_safety_status',production_safety_status
  write(10,'(A,1X,I0)')'stage7_rhs_pressure_poisson_modified_flag',0;write(10,'(A,1X,I0)')'stage7_rhs_projection_modified_flag',0;write(10,'(A,1X,I0)')'stage7_rhs_real_projection_called_flag',0;write(10,'(A,1X,I0)')'stage7_rhs_production_dns_called_flag',0;write(10,'(A,1X,I0)')'stage7_rhs_fluid_update_called_flag',0;write(10,'(A,1X,I0)')'stage7_rhs_fibre_advance_called_flag',0;write(10,'(A,1X,I0)')'stage7_rhs_no_projection_no_dns_status',no_projection_dns_status
  write(10,'(A,1X,ES24.16E3)')'stage7_rhs_candidate_force_eulerian_x',fe(1)
  write(10,'(A,1X,ES24.16E3)')'stage7_rhs_candidate_force_eulerian_y',fe(2)
  write(10,'(A,1X,ES24.16E3)')'stage7_rhs_candidate_force_eulerian_z',fe(3)
  write(10,'(A,1X,ES24.16E3)')'stage7_rhs_candidate_force_lagrangian_x',fl(1)
  write(10,'(A,1X,ES24.16E3)')'stage7_rhs_candidate_force_lagrangian_y',fl(2)
  write(10,'(A,1X,ES24.16E3)')'stage7_rhs_candidate_force_lagrangian_z',fl(3)
  write(10,'(A,1X,ES24.16E3)')'stage7_rhs_candidate_force_conservation_error',force_abs_err
  write(10,'(A,1X,ES24.16E3)')'stage7_rhs_candidate_force_conservation_relative_error',force_rel_err
  write(10,'(A,1X,I0)')'stage7_rhs_candidate_force_conservation_status',force_cons_status
  write(10,'(A,1X,I0)')'stage7_rhs_candidate_check_status',status
  close(10)
contains
  subroutine dependency_status(s6,s70,s71,s72,s73,s74,s75,s76,s77,s78)
    integer,intent(out)::s6,s70,s71,s72,s73,s74,s75,s76,s77,s78
    s6=1;s70=1;s71=1;s72=1;s73=1;s74=1;s75=1;s76=1;s77=1;s78=1
  end subroutine

  subroutine max_increment_error_to_scaled(rx0,ry0,rz0,rx,ry,rz,sx,sy,sz,scale,m)
    real(mytype),intent(in)::rx0(:,:,:),ry0(:,:,:),rz0(:,:,:),rx(:,:,:),ry(:,:,:),rz(:,:,:),sx(:,:,:),sy(:,:,:),sz(:,:,:),scale
    real(mytype),intent(out)::m
    real(mytype)::ex,ey,ez
    ex=maxval(abs((rx-rx0)-scale*sx))
    ey=maxval(abs((ry-ry0)-scale*sy))
    ez=maxval(abs((rz-rz0)-scale*sz))
    m=max(ex,max(ey,ez))
  end subroutine

  subroutine max_increment_scaling_error(sx,sy,sz,m)
    real(mytype),intent(in)::sx(:,:,:),sy(:,:,:),sz(:,:,:)
    real(mytype),intent(out)::m
    real(mytype)::ex,ey,ez
    ex=maxval(abs((sx/4._mytype)-0.5_mytype*(sx/2._mytype)))
    ey=maxval(abs((sy/4._mytype)-0.5_mytype*(sy/2._mytype)))
    ez=maxval(abs((sz/4._mytype)-0.5_mytype*(sz/2._mytype)))
    m=max(ex,max(ey,ez))
  end subroutine

  subroutine run_blocked_case(grid,layout,flag,ds,bfx,bfy,bfz,bvc,bbc,buc,bfchg,brhschg,binj,bmod,bhook,brej)
    type(stage7_channel_grid_t),intent(in)::grid
    type(stage7_velocity_layout_t),intent(in)::layout
    real(mytype),intent(in)::flag(3,nlag),ds(nlag)
    real(mytype),intent(inout)::bfx(:,:,:),bfy(:,:,:),bfz(:,:,:)
    integer,intent(out)::bvc,bbc,buc,binj,bmod,bhook,brej
    real(mytype),intent(out)::bfchg,brhschg
    real(mytype)::bx(3,nlag)
    real(mytype),allocatable::rx(:,:,:),ry(:,:,:),rz(:,:,:),r0x(:,:,:),r0y(:,:,:),r0z(:,:,:)
    type(stage6_config_t)::bcfg
    integer::l
    allocate(rx(size(bfx,1),size(bfx,2),size(bfx,3)),ry(size(bfy,1),size(bfy,2),size(bfy,3)),rz(size(bfz,1),size(bfz,2),size(bfz,3)),r0x(size(bfx,1),size(bfx,2),size(bfx,3)),r0y(size(bfy,1),size(bfy,2),size(bfy,3)),r0z(size(bfz,1),size(bfz,2),size(bfz,3)))
    do l=1,nlag
      bx(1,l)=grid%xmin+2*l*grid%dx
      bx(2,l)=grid%ymin-0.1_mytype*(grid%ymax-grid%ymin)
      bx(3,l)=grid%zmin+2*l*grid%dz
    end do
    call build_stage7_force_density_candidate(grid,layout,nlag,bx,flag,ds,bfx,bfy,bfz,bvc,bbc,buc)
    bfchg=max(maxval(abs(bfx)),max(maxval(abs(bfy)),maxval(abs(bfz))))
    call init_stage6_controlled_test_config(bcfg)
    rx=0;ry=0;rz=0; r0x=0;r0y=0;r0z=0
    call apply_stage7_candidate_to_rhs_controlled(bcfg,rx,ry,rz,bfx,bfy,bfz,binj,bmod,bhook,brej)
    call max_increment_error_to_scaled(r0x,r0y,r0z,rx,ry,rz,0._mytype*bfx,0._mytype*bfy,0._mytype*bfz,0._mytype,brhschg)
    if (bfchg<=1e-14_mytype .and. brhschg<=1e-14_mytype) binj=0
    deallocate(rx,ry,rz,r0x,r0y,r0z)
  end subroutine
end program
