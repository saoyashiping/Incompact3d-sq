program fibre_stage7_scalar_interpolation_robustness_check
  use fibre_parameters, only: mytype
  use fibre_stage7_grid_metadata
  use fibre_stage7_scalar_interpolation
  implicit none
  type(stage7_channel_grid_t)::grid
  type(stage7_interp_weight_t)::w1,w2
  integer,parameter::nx=16,ny=17,nz=12,np=6
  real(mytype),allocatable::phi(:,:,:),rhsx(:,:,:),rhsy(:,:,:),rhsz(:,:,:)
  real(mytype)::xp(np),yp(np),zp(np),v1,v2,expv,err
  real(mytype)::rep_w_err,rep_sum_err,wsum_err,qy_err,cy_err,mix_err,pshift_err,noop_change
  integer::i,j,k,p,q,valid,rej
  integer::dep6,dep70,dep71,dep72,rep_mismatch,rep_ok,wsum_ok,qy_ok,cy_ok,mix_ok,px_ok,pz_ok,ylo_ok,yhi_ok,yout_ok,nw_count,nw_ok
  integer::noop_mod,pp,proj,rp,dns,flu,fib,status,valid_count

  call init_stage7_nonuniform_channel_grid(grid,nx,ny,nz); call validate_stage7_channel_grid(grid,valid,rej)
  allocate(phi(nx,ny,nz),rhsx(nx,ny,nz),rhsy(nx,ny,nz),rhsz(nx,ny,nz)); call init_stage7_interp_weight(w1,64); call init_stage7_interp_weight(w2,64)
  call dependency_status(dep6,dep70,dep71,dep72)

  xp=(/grid%xmin+3.2_mytype*grid%dx,grid%xmin+6.4_mytype*grid%dx,grid%xmin+8.3_mytype*grid%dx,grid%xmin+10.2_mytype*grid%dx,grid%xmin+12.1_mytype*grid%dx,grid%xmin+14.0_mytype*grid%dx/)
  yp=(/grid%y_center(6),grid%y_center(7),grid%y_center(8),grid%y_center(9),grid%y_center(10),grid%y_center(11)/)
  zp=(/grid%zmin+2.1_mytype*grid%dz,grid%zmin+3.3_mytype*grid%dz,grid%zmin+4.7_mytype*grid%dz,grid%zmin+6.2_mytype*grid%dz,grid%zmin+7.9_mytype*grid%dz,grid%zmin+9.4_mytype*grid%dz/)

  rep_mismatch=0; rep_w_err=0; rep_sum_err=0; wsum_err=0; valid_count=0
  do p=1,np
    call build_stage7_scalar_interp_weight(grid,xp(p),yp(p),zp(p),w1)
    call build_stage7_scalar_interp_weight(grid,xp(p),yp(p),zp(p),w2)
    if (w1%n/=w2%n) rep_mismatch=rep_mismatch+1
    if (w1%n>0 .and. w2%n>0) then
      do q=1,w1%n
        if (w1%i(q)/=w2%i(q) .or. w1%j(q)/=w2%j(q) .or. w1%k(q)/=w2%k(q)) rep_mismatch=rep_mismatch+1
      end do
      rep_w_err=max(rep_w_err,maxval(abs(w1%w(1:w1%n)-w2%w(1:w2%n))))
      rep_sum_err=max(rep_sum_err,abs(w1%weight_sum-w2%weight_sum))
    end if
    if (w1%valid_flag==1) then
      valid_count=valid_count+1
      wsum_err=max(wsum_err,abs(w1%weight_sum-1._mytype))
    end if
  end do
  rep_ok=merge(1,0,rep_mismatch==0 .and. rep_w_err<=1e-14_mytype .and. rep_sum_err<=1e-14_mytype)
  wsum_ok=merge(1,0,wsum_err<=1e-12_mytype .and. valid_count>=5)

  do k=1,nz; do j=1,ny; do i=1,nx; phi(i,j,k)=grid%y_center(j)**2; end do; end do; end do
  qy_err=0; do p=1,np; call build_stage7_scalar_interp_weight(grid,xp(p),yp(p),zp(p),w1); call interpolate_stage7_scalar(phi,w1,v1); qy_err=max(qy_err,abs(v1-yp(p)**2)); end do
  qy_ok=merge(1,0,qy_err<=1e-12_mytype)

  do k=1,nz; do j=1,ny; do i=1,nx; phi(i,j,k)=grid%y_center(j)**3-0.2_mytype*grid%y_center(j)**2+0.1_mytype*grid%y_center(j)+0.3_mytype; end do; end do; end do
  cy_err=0; do p=1,np; call build_stage7_scalar_interp_weight(grid,xp(p),yp(p),zp(p),w1); call interpolate_stage7_scalar(phi,w1,v1); expv=yp(p)**3-0.2_mytype*yp(p)**2+0.1_mytype*yp(p)+0.3_mytype; cy_err=max(cy_err,abs(v1-expv)); end do
  cy_ok=merge(1,0,cy_err<=1e-11_mytype)

  do k=1,nz; do j=1,ny; do i=1,nx
    phi(i,j,k)=1._mytype+0.2_mytype*(grid%xmin+real(i-1,mytype)*grid%dx)-0.3_mytype*grid%y_center(j)+0.1_mytype*(grid%zmin+real(k-1,mytype)*grid%dz)+0.05_mytype*grid%y_center(j)**2-0.02_mytype*grid%y_center(j)**3
  end do; end do; end do
  mix_err=0; do p=1,4; call build_stage7_scalar_interp_weight(grid,xp(p),yp(p),zp(p),w1); call interpolate_stage7_scalar(phi,w1,v1); expv=1._mytype+0.2_mytype*xp(p)-0.3_mytype*yp(p)+0.1_mytype*zp(p)+0.05_mytype*yp(p)**2-0.02_mytype*yp(p)**3; mix_err=max(mix_err,abs(v1-expv)); end do
  mix_ok=merge(1,0,mix_err<=1e-11_mytype)

  do k=1,nz; do j=1,ny; do i=1,nx
    phi(i,j,k)=sin(2._mytype*stage7_pi*(grid%xmin+real(i-1,mytype)*grid%dx-grid%xmin)/(grid%xmax-grid%xmin))+cos(2._mytype*stage7_pi*(grid%zmin+real(k-1,mytype)*grid%dz-grid%zmin)/(grid%zmax-grid%zmin))+0.25_mytype*grid%y_center(j)
  end do; end do; end do
  pshift_err=0
  call build_stage7_scalar_interp_weight(grid,grid%xmin+0.1_mytype*grid%dx,yp(3),zp(3),w1); call interpolate_stage7_scalar(phi,w1,v1)
  call build_stage7_scalar_interp_weight(grid,grid%xmin+0.1_mytype*grid%dx+(grid%xmax-grid%xmin),yp(3),zp(3),w2); call interpolate_stage7_scalar(phi,w2,v2); pshift_err=max(pshift_err,abs(v1-v2))
  call build_stage7_scalar_interp_weight(grid,grid%xmin+0.1_mytype*grid%dx-(grid%xmax-grid%xmin),yp(3),zp(3),w2); call interpolate_stage7_scalar(phi,w2,v2); pshift_err=max(pshift_err,abs(v1-v2))
  call build_stage7_scalar_interp_weight(grid,xp(3),yp(3),grid%zmin+0.1_mytype*grid%dz,w1); call interpolate_stage7_scalar(phi,w1,v1)
  call build_stage7_scalar_interp_weight(grid,xp(3),yp(3),grid%zmin+0.1_mytype*grid%dz+(grid%zmax-grid%zmin),w2); call interpolate_stage7_scalar(phi,w2,v2); pshift_err=max(pshift_err,abs(v1-v2))
  call build_stage7_scalar_interp_weight(grid,xp(3),yp(3),grid%zmin+0.1_mytype*grid%dz-(grid%zmax-grid%zmin),w2); call interpolate_stage7_scalar(phi,w2,v2); pshift_err=max(pshift_err,abs(v1-v2))
  px_ok=merge(1,0,pshift_err<=1e-12_mytype); pz_ok=px_ok

  call build_stage7_scalar_interp_weight(grid,xp(1),grid%ymin-0.1_mytype*abs(grid%ymax-grid%ymin),zp(1),w1); ylo_ok=merge(1,0,w1%valid_flag==0 .and. w1%unsafe_flag==1)
  call build_stage7_scalar_interp_weight(grid,xp(1),grid%ymax+0.1_mytype*abs(grid%ymax-grid%ymin),zp(1),w1); yhi_ok=merge(1,0,w1%valid_flag==0 .and. w1%unsafe_flag==1)
  yout_ok=merge(1,0,ylo_ok==1 .and. yhi_ok==1)

  nw_count=0
  call build_stage7_scalar_interp_weight(grid,xp(1),grid%y_center(1),zp(1),w1); if (w1%valid_flag==0 .and. w1%unsafe_flag==1) nw_count=nw_count+1
  call build_stage7_scalar_interp_weight(grid,xp(1),grid%y_center(ny),zp(1),w1); if (w1%valid_flag==0 .and. w1%unsafe_flag==1) nw_count=nw_count+1
  nw_ok=merge(1,0,nw_count>0)

  rhsx=1; rhsy=2; rhsz=3
  call stage7_grid_noop_rhs_guard(rhsx,rhsy,rhsz,noop_change,noop_mod); call stage7_grid_pressure_status(pp,proj,rp,dns,flu,fib)

  status=merge(1,0,dep6==1 .and. dep70==1 .and. dep71==1 .and. dep72==1 .and. rep_ok==1 .and. wsum_ok==1 .and. qy_ok==1 .and. cy_ok==1 .and. mix_ok==1 .and. px_ok==1 .and. pz_ok==1 .and. yout_ok==1 .and. nw_ok==1 .and. noop_change<=1e-14_mytype .and. noop_mod==0 .and. pp==0 .and. proj==0 .and. dns==0 .and. flu==0 .and. fib==0)

  open(10,file='stage7_outputs/fibre_stage7_scalar_interpolation_robustness_check.dat',status='replace')
  write(10,'(A,1X,I0)') 'stage7_interp_robust_stage6_dependency_status',dep6
  write(10,'(A,1X,I0)') 'stage7_interp_robust_stage7_0_dependency_status',dep70
  write(10,'(A,1X,I0)') 'stage7_interp_robust_stage7_1_dependency_status',dep71
  write(10,'(A,1X,I0)') 'stage7_interp_robust_stage7_2_dependency_status',dep72
  write(10,'(A,1X,I0)') 'stage7_interp_robust_repeat_index_mismatch_count',rep_mismatch
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_robust_repeat_weight_error_max',rep_w_err
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_robust_repeat_weight_sum_error_max',rep_sum_err
  write(10,'(A,1X,I0)') 'stage7_interp_robust_repeatability_status',rep_ok
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_robust_weight_sum_error_max',wsum_err
  write(10,'(A,1X,I0)') 'stage7_interp_robust_weight_sum_status',wsum_ok
  write(10,'(A,1X,I0)') 'stage7_interp_robust_valid_weight_count',valid_count
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_robust_quadratic_y_error_max',qy_err
  write(10,'(A,1X,I0)') 'stage7_interp_robust_quadratic_y_status',qy_ok
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_robust_cubic_y_error_max',cy_err
  write(10,'(A,1X,I0)') 'stage7_interp_robust_cubic_y_status',cy_ok
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_robust_mixed_field_error_max',mix_err
  write(10,'(A,1X,I0)') 'stage7_interp_robust_mixed_field_status',mix_ok
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_robust_periodic_shift_error_max',pshift_err
  write(10,'(A,1X,I0)') 'stage7_interp_robust_periodic_x_shift_status',px_ok
  write(10,'(A,1X,I0)') 'stage7_interp_robust_periodic_z_shift_status',pz_ok
  write(10,'(A,1X,I0)') 'stage7_interp_robust_y_outside_low_blocked_flag',ylo_ok
  write(10,'(A,1X,I0)') 'stage7_interp_robust_y_outside_high_blocked_flag',yhi_ok
  write(10,'(A,1X,I0)') 'stage7_interp_robust_y_outside_status',yout_ok
  write(10,'(A,1X,I0)') 'stage7_interp_robust_nearwall_unsafe_count',nw_count
  write(10,'(A,1X,I0)') 'stage7_interp_robust_nearwall_blocked_flag',nw_ok
  write(10,'(A,1X,ES24.16E3)') 'stage7_interp_robust_noop_rhs_change_max',noop_change
  write(10,'(A,1X,I0)') 'stage7_interp_robust_noop_rhs_modified_flag',noop_mod
  write(10,'(A,1X,I0)') 'stage7_interp_robust_pressure_poisson_modified_flag',pp
  write(10,'(A,1X,I0)') 'stage7_interp_robust_projection_modified_flag',proj
  write(10,'(A,1X,I0)') 'stage7_interp_robust_production_dns_called_flag',dns
  write(10,'(A,1X,I0)') 'stage7_interp_robust_fluid_update_called_flag',flu
  write(10,'(A,1X,I0)') 'stage7_interp_robust_fibre_advance_called_flag',fib
  write(10,'(A,1X,I0)') 'stage7_scalar_interpolation_robustness_check_status',status
  close(10)
contains
  subroutine dependency_status(s6,s70,s71,s72)
    integer,intent(out)::s6,s70,s71,s72
    s6=check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_smoke_check_status',1._mytype)*check_dat('stage6_outputs/fibre_stage6_total_smoke_check.dat','stage6_total_all_prior_outputs_exist',1._mytype)
    s70=check_dat('stage7_outputs/fibre_stage7_config_check.dat','stage7_config_check_status',1._mytype)
    s71=check_dat('stage7_outputs/fibre_stage7_grid_metadata_check.dat','stage7_grid_metadata_check_status',1._mytype)
    s72=check_dat('stage7_outputs/fibre_stage7_scalar_interpolation_check.dat','stage7_scalar_interpolation_check_status',1._mytype)
  end subroutine
  integer function check_dat(path,key,val)
    character(*),intent(in)::path,key
    real(mytype),intent(in)::val
    character(len=512)::line,rest
    real(mytype)::v
    integer::u,ios
    check_dat=0; if (.not.file_exists(path)) return
    open(newunit=u,file=path,status='old',iostat=ios); if (ios/=0) return
    do
      read(u,'(A)',iostat=ios) line; if (ios/=0) exit
      line=adjustl(line)
      if (index(line,trim(key))==1) then
        rest=adjustl(line(len_trim(key)+1:)); if (len_trim(rest)>0 .and. rest(1:1)=='=') rest=adjustl(rest(2:))
        read(rest,*,iostat=ios) v
        if (ios==0 .and. abs(v-val)<=1e-12_mytype) then; check_dat=1; exit; end if
      end if
    end do
    close(u)
  end function
  logical function file_exists(path)
    character(*),intent(in)::path
    inquire(file=path,exist=file_exists)
  end function
end program
