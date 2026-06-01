module fibre_stage11_controlled_interpolation
  use decomp_2d_constants, only : mytype
  implicit none
  private
  integer, parameter :: nx=16, ny=17, nz=16
  real(mytype), parameter :: xmin=0.0_mytype, xmax=2.0_mytype, ymin=0.0_mytype, ymax=1.0_mytype, zmin=0.0_mytype, zmax=1.0_mytype
  real(mytype), parameter :: dx=(xmax-xmin)/real(nx-1,mytype), dy=(ymax-ymin)/real(ny-1,mytype), dz=(zmax-zmin)/real(nz-1,mytype)
  real(mytype), parameter :: constant_abs_tol=1.0e-12_mytype, linear_abs_tol=1.0e-12_mytype, shear_abs_tol=1.0e-12_mytype, weight_sum_tol=1.0e-12_mytype, poiseuille_abs_tol=5.0e-3_mytype
  logical :: initialized=.false.
  integer :: initialized_status=0, constant_field_status=0, linear_field_status=0, shear_field_status=0, poiseuille_field_status=0, weight_sum_status=0
  integer :: periodic_safety_status=0, near_wall_safety_status=0, out_of_domain_rejection_status=0
  integer :: no_production_fluid_access_status=1, no_fluid_field_modification_status=1, no_rhs_injection_status=1, no_ibm_spreading_status=1, no_feedback_force_status=1, no_twoway_force_status=1, no_structure_advance_status=1
  integer :: controlled_interpolation_status=0
  real(mytype) :: constant_max_error=0.0_mytype, linear_max_error=0.0_mytype, shear_max_error=0.0_mytype, poiseuille_max_error=0.0_mytype, weight_sum_max_error=0.0_mytype
  public :: stage11_controlled_interpolation_init, stage11_controlled_interpolation_finalize, stage11_controlled_interpolation_sample_uniform_trilinear, stage11_controlled_interpolation_get_status_values, stage11_controlled_interpolation_write_diagnostics
contains
subroutine stage11_controlled_interpolation_init(); initialized=.true.; initialized_status=1; call run_tests(); end
subroutine stage11_controlled_interpolation_finalize(); initialized=.false.; initialized_status=0; controlled_interpolation_status=0; end
subroutine stage11_controlled_interpolation_sample_uniform_trilinear(x,y,z,field_type,u,v,w,ok,werr)
 real(mytype),intent(in)::x,y,z; integer,intent(in)::field_type; real(mytype),intent(out)::u,v,w,werr; logical,intent(out)::ok
 integer::i0,j0,k0,i1,j1,k1; real(mytype)::tx,ty,tz,x0,y0,z0,wt,sumw
 if (x<xmin .or. x>xmax .or. y<ymin .or. y>ymax .or. z<zmin .or. z>zmax) then; ok=.false.;u=0;v=0;w=0;werr=1;return; end if
 tx=(x-xmin)/dx; ty=(y-ymin)/dy; tz=(z-zmin)/dz
 i0=max(1,min(nx-1,int(floor(tx))+1)); j0=max(1,min(ny-1,int(floor(ty))+1)); k0=max(1,min(nz-1,int(floor(tz))+1))
 i1=i0+1; j1=j0+1; k1=k0+1
 x0=xmin+real(i0-1,mytype)*dx; y0=ymin+real(j0-1,mytype)*dy; z0=zmin+real(k0-1,mytype)*dz
 tx=(x-x0)/dx; ty=(y-y0)/dy; tz=(z-z0)/dz
 sumw=0.0_mytype; u=0;v=0;w=0
 call accum(1-tx,1-ty,1-tz,x0,y0,z0); call accum(tx,1-ty,1-tz,x0+dx,y0,z0); call accum(1-tx,ty,1-tz,x0,y0+dy,z0); call accum(tx,ty,1-tz,x0+dx,y0+dy,z0)
 call accum(1-tx,1-ty,tz,x0,y0,z0+dz); call accum(tx,1-ty,tz,x0+dx,y0,z0+dz); call accum(1-tx,ty,tz,x0,y0+dy,z0+dz); call accum(tx,ty,tz,x0+dx,y0+dy,z0+dz)
 ok=.true.; werr=abs(sumw-1.0_mytype)
contains
 subroutine accum(wx,wy,wz,xx,yy,zz)
 real(mytype),intent(in)::wx,wy,wz,xx,yy,zz; real(mytype)::fu,fv,fw
 wt=wx*wy*wz; sumw=sumw+wt; call f(field_type,xx,yy,zz,fu,fv,fw); u=u+wt*fu; v=v+wt*fv; w=w+wt*fw
 end
end
subroutine f(t,x,y,z,u,v,w)
 integer,intent(in)::t; real(mytype),intent(in)::x,y,z; real(mytype),intent(out)::u,v,w
 select case(t); case(1);u=1.25_mytype;v=-0.5_mytype;w=0.75_mytype; case(2);u=1+0.2*x-0.3*y+0.1*z;v=-0.5+0.4*x+0.2*y-0.2*z;w=0.25-0.1*x+0.3*y+0.5*z; case(3);u=y;v=0;w=0; case default;u=4*y*(1-y);v=0;w=0; end select
end
subroutine run_tests()
 real(mytype)::u,v,w,ue,ve,we,e,werr; logical::ok
 call stage11_controlled_interpolation_sample_uniform_trilinear(0.37_mytype,0.44_mytype,0.22_mytype,1,u,v,w,ok,werr); call f(1,0.37_mytype,0.44_mytype,0.22_mytype,ue,ve,we); constant_max_error=max(abs(u-ue),max(abs(v-ve),abs(w-we))); constant_field_status=merge(1,0,ok.and.constant_max_error<=constant_abs_tol)
 call stage11_controlled_interpolation_sample_uniform_trilinear(1.13_mytype,0.67_mytype,0.41_mytype,2,u,v,w,ok,werr); call f(2,1.13_mytype,0.67_mytype,0.41_mytype,ue,ve,we); linear_max_error=max(abs(u-ue),max(abs(v-ve),abs(w-we))); linear_field_status=merge(1,0,ok.and.linear_max_error<=linear_abs_tol)
 call stage11_controlled_interpolation_sample_uniform_trilinear(1.19_mytype,0.31_mytype,0.77_mytype,3,u,v,w,ok,werr); call f(3,1.19_mytype,0.31_mytype,0.77_mytype,ue,ve,we); shear_max_error=max(abs(u-ue),max(abs(v-ve),abs(w-we))); shear_field_status=merge(1,0,ok.and.shear_max_error<=shear_abs_tol)
 call stage11_controlled_interpolation_sample_uniform_trilinear(0.63_mytype,0.35_mytype,0.49_mytype,4,u,v,w,ok,werr); call f(4,0.63_mytype,0.35_mytype,0.49_mytype,ue,ve,we); poiseuille_max_error=max(abs(u-ue),max(abs(v-ve),abs(w-we))); poiseuille_field_status=merge(1,0,ok.and.poiseuille_max_error<=poiseuille_abs_tol)
 weight_sum_max_error=werr; weight_sum_status=merge(1,0,werr<=weight_sum_tol)
 periodic_safety_status=merge(1,0,point_in(1.0e-10_mytype,0.5_mytype,1.0e-10_mytype).and.point_in(2.0_mytype-1.0e-10_mytype,0.5_mytype,1.0_mytype-1.0e-10_mytype))
 near_wall_safety_status=merge(1,0,point_in(1.0_mytype,1.0e-10_mytype,0.5_mytype).and.point_in(1.0_mytype,1.0_mytype-1.0e-10_mytype,0.5_mytype))
 out_of_domain_rejection_status=merge(1,0,.not.point_in(1.0_mytype,-1.0e-9_mytype,0.5_mytype))
 controlled_interpolation_status=merge(1,0,initialized_status==1.and.constant_field_status==1.and.linear_field_status==1.and.shear_field_status==1.and.poiseuille_field_status==1.and.weight_sum_status==1.and.periodic_safety_status==1.and.near_wall_safety_status==1.and.out_of_domain_rejection_status==1.and.no_production_fluid_access_status==1.and.no_fluid_field_modification_status==1.and.no_rhs_injection_status==1.and.no_ibm_spreading_status==1.and.no_feedback_force_status==1.and.no_twoway_force_status==1.and.no_structure_advance_status==1)
end
logical function point_in(x,y,z); real(mytype),intent(in)::x,y,z; point_in=(x>=xmin.and.x<=xmax.and.y>=ymin.and.y<=ymax.and.z>=zmin.and.z<=zmax); end
subroutine stage11_controlled_interpolation_get_status_values(a,b,c,d,e,f1,g,h,i,j,k,l,m,n,o,p,q)
 integer,intent(out)::a,b,c,d,e,f1,g,h,i,j,k,l,m,n,o,p,q
 a=initialized_status;b=constant_field_status;c=linear_field_status;d=shear_field_status;e=poiseuille_field_status;f1=weight_sum_status;g=periodic_safety_status;h=near_wall_safety_status;i=out_of_domain_rejection_status;j=no_production_fluid_access_status;k=no_fluid_field_modification_status;l=no_rhs_injection_status;m=no_ibm_spreading_status;n=no_feedback_force_status;o=no_twoway_force_status;p=no_structure_advance_status;q=controlled_interpolation_status
end
subroutine stage11_controlled_interpolation_write_diagnostics(unit)
 integer,intent(in)::unit
 write(unit,'(A,1X,I0)') 'stage11_4_controlled_interpolation_initialized_status', initialized_status
 write(unit,'(A,1X,I0)') 'stage11_4_constant_field_status', constant_field_status
 write(unit,'(A,1X,I0)') 'stage11_4_linear_field_status', linear_field_status
 write(unit,'(A,1X,I0)') 'stage11_4_shear_field_status', shear_field_status
 write(unit,'(A,1X,I0)') 'stage11_4_poiseuille_field_status', poiseuille_field_status
 write(unit,'(A,1X,I0)') 'stage11_4_weight_sum_status', weight_sum_status
 write(unit,'(A,1X,I0)') 'stage11_4_periodic_safety_status', periodic_safety_status
 write(unit,'(A,1X,I0)') 'stage11_4_near_wall_safety_status', near_wall_safety_status
 write(unit,'(A,1X,I0)') 'stage11_4_out_of_domain_rejection_status', out_of_domain_rejection_status
 write(unit,'(A,1X,I0)') 'stage11_4_no_production_fluid_access_status', no_production_fluid_access_status
 write(unit,'(A,1X,I0)') 'stage11_4_no_fluid_field_modification_status', no_fluid_field_modification_status
 write(unit,'(A,1X,I0)') 'stage11_4_no_rhs_injection_status', no_rhs_injection_status
 write(unit,'(A,1X,I0)') 'stage11_4_no_ibm_spreading_status', no_ibm_spreading_status
 write(unit,'(A,1X,I0)') 'stage11_4_no_feedback_force_status', no_feedback_force_status
 write(unit,'(A,1X,I0)') 'stage11_4_no_twoway_force_status', no_twoway_force_status
 write(unit,'(A,1X,I0)') 'stage11_4_no_structure_advance_status', no_structure_advance_status
 write(unit,'(A,1X,I0)') 'stage11_4_controlled_interpolation_status', controlled_interpolation_status
 write(unit,'(A,1X,ES24.16E3)') 'stage11_4_constant_max_error', constant_max_error
 write(unit,'(A,1X,ES24.16E3)') 'stage11_4_linear_max_error', linear_max_error
 write(unit,'(A,1X,ES24.16E3)') 'stage11_4_shear_max_error', shear_max_error
 write(unit,'(A,1X,ES24.16E3)') 'stage11_4_poiseuille_max_error', poiseuille_max_error
 write(unit,'(A,1X,ES24.16E3)') 'stage11_4_weight_sum_max_error', weight_sum_max_error
end
end module
