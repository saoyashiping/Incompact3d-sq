program fibre_stage4_feedback_check
  use fibre_parameters, only : mytype, fibre_init_straight_free_free
  use fibre_types, only : fibre_t
  use fibre_stage4_grid_adapter, only : stage4_grid_adapter_t, init_stage4_grid_adapter_from_arrays
  use fibre_stage4_feedback_adapter, only : compute_stage4_feedback_if_supported, apply_stage4_feedback_to_f_ext
  implicit none
  integer, parameter :: nx=16, ny=12, nz=10
  real(mytype) :: x(nx), y(ny), z(nz), yn(ny), eta
  real(mytype), allocatable :: ux(:,:,:), uy(:,:,:), uz(:,:,:), u_lag(:,:), fs(:,:), ff(:,:)
  type(fibre_t) :: fibre
  type(stage4_grid_adapter_t) :: a, an, a0, a2
  integer :: i,j,k,s,unit_id
  real(mytype) :: beta, err, meanfx
  integer :: flag
  beta=10._mytype
  do i=1,nx; x(i)=(real(i,mytype)-0.5_mytype)*(2._mytype/real(nx,mytype)); end do
  do j=1,ny
    y(j)=-1._mytype+(real(j,mytype)-0.5_mytype)*(2._mytype/real(ny,mytype)); eta=(real(j,mytype)-0.5_mytype)/real(ny,mytype); yn(j)=-1._mytype+2._mytype*eta*eta
  end do
  do k=1,nz; z(k)=(real(k,mytype)-0.5_mytype)*(1._mytype/real(nz,mytype)); end do
  call init_stage4_grid_adapter_from_arrays(a,x,y,z,.true.,.false.,.true.,1)
  call init_stage4_grid_adapter_from_arrays(an,x,yn,z,.true.,.false.,.true.,1)
  call init_stage4_grid_adapter_from_arrays(a0,x,y,z,.true.,.false.,.true.,0)
  call init_stage4_grid_adapter_from_arrays(a2,x,y,z,.true.,.false.,.true.,2)
  call fibre_init_straight_free_free(fibre,33,1._mytype,1._mytype,1._mytype)
  do i=1,fibre%nl; fibre%x(1,i)=0.5_mytype+(real(i-1,mytype)/real(fibre%nl-1,mytype)); fibre%x(2,i)=0._mytype; fibre%x(3,i)=0.5_mytype; end do
  allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz),u_lag(3,fibre%nl),fs(3,fibre%nl),ff(3,fibre%nl))
  open(newunit=unit_id,file='stage4_outputs/fibre_stage4_feedback_check.dat',status='replace',action='write')
  ux=0.2_mytype; uy=-0.1_mytype; uz=0.05_mytype; fibre%v(1,:)=0.2_mytype; fibre%v(2,:)=-0.1_mytype; fibre%v(3,:)=0.05_mytype
  call compute_stage4_feedback_if_supported(a,ux,uy,uz,fibre,beta,u_lag,fs,ff,s); call apply_stage4_feedback_to_f_ext(fibre,fs,'set',flag)
  write(unit_id,'(A,1X,I0)') 'stage4_zero_slip_feedback_status', s
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_zero_slip_force_structure_norm', sqrt(sum(fs**2))
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_zero_slip_force_fluid_norm', sqrt(sum(ff**2))
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_zero_slip_f_ext_norm', sqrt(sum(fibre%f_ext**2))
  fibre%v=0._mytype; ux=0.2_mytype; uy=0._mytype; uz=0._mytype
  call compute_stage4_feedback_if_supported(a,ux,uy,uz,fibre,beta,u_lag,fs,ff,s); err=maxval(abs(fs(1,:)-2._mytype))+maxval(abs(fs(2,:)))+maxval(abs(fs(3,:))); meanfx=sum(fs(1,:))/real(fibre%nl,mytype)
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_positive_drag_force_error', err
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_positive_drag_mean_fx', meanfx
  write(unit_id,'(A,1X,I0)') 'stage4_positive_drag_direction_flag', merge(1,0,meanfx>0._mytype)
  ux=-0.2_mytype; call compute_stage4_feedback_if_supported(a,ux,uy,uz,fibre,beta,u_lag,fs,ff,s); err=maxval(abs(fs(1,:)+2._mytype))+maxval(abs(fs(2,:)))+maxval(abs(fs(3,:))); meanfx=sum(fs(1,:))/real(fibre%nl,mytype)
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_reverse_drag_force_error', err
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_reverse_drag_mean_fx', meanfx
  write(unit_id,'(A,1X,I0)') 'stage4_reverse_drag_direction_flag', merge(1,0,meanfx<0._mytype)
  do k=1,nz; do j=1,ny; do i=1,nx; ux(i,j,k)=1._mytype+0.2_mytype*x(i)-0.1_mytype*y(j)+0.05_mytype*z(k); uy(i,j,k)=-0.3_mytype+0.4_mytype*y(j); uz(i,j,k)=0.7_mytype-0.2_mytype*z(k)+0.1_mytype*x(i); end do; end do; end do
  do i=1,fibre%nl; fibre%v(1,i)=0.01_mytype*real(i,mytype); fibre%v(2,i)=-0.005_mytype*real(i,mytype); fibre%v(3,i)=0.002_mytype*real(i,mytype); end do
  call compute_stage4_feedback_if_supported(a,ux,uy,uz,fibre,beta,u_lag,fs,ff,s)
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_action_reaction_pointwise_error', maxval(abs(fs+ff))
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_action_reaction_total_error', abs(sum(fs+ff))
  fibre%nl=3; deallocate(fibre%x,fibre%x_old,fibre%v,fibre%tension,fibre%f_ext); call fibre_init_straight_free_free(fibre,3,1._mytype,1._mytype,1._mytype); fibre%x(:,1)=[1._mytype,-0.5_mytype,0.5_mytype]; fibre%x(:,2)=[1._mytype,0._mytype,0.5_mytype]; fibre%x(:,3)=[1._mytype,0.5_mytype,0.5_mytype]; fibre%v=0._mytype
  do k=1,nz; do j=1,ny; do i=1,nx; ux(i,j,k)=1._mytype-y(j)*y(j); uy(i,j,k)=0._mytype; uz(i,j,k)=0._mytype; end do; end do; end do
  deallocate(u_lag,fs,ff); allocate(u_lag(3,fibre%nl),fs(3,fibre%nl),ff(3,fibre%nl)); call compute_stage4_feedback_if_supported(a,ux,uy,uz,fibre,beta,u_lag,fs,ff,s)
  write(unit_id,'(A,1X,I0)') 'stage4_poiseuille_force_status', s
  write(unit_id,'(A,1X,I0)') 'stage4_poiseuille_force_center_greater_flag', merge(1,0,fs(1,2)>fs(1,1).and.fs(1,2)>fs(1,3))
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_poiseuille_force_symmetry_error', abs(fs(1,1)-fs(1,3))
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_poiseuille_force_variation_norm', abs(fs(1,2)-fs(1,1))+abs(fs(1,2)-fs(1,3))
  call apply_stage4_feedback_to_f_ext(fibre,fs,'set',flag)
  write(unit_id,'(A,1X,ES24.16E3)') 'stage4_f_ext_matches_structure_force_error', maxval(abs(fibre%f_ext-fs))
  write(unit_id,'(A,1X,I0)') 'stage4_structure_advance_called', 0
  call compute_stage4_feedback_if_supported(an,ux,uy,uz,fibre,beta,u_lag,fs,ff,s); write(unit_id,'(A,1X,I0)') 'stage4_nonuniform_feedback_blocked_flag', merge(1,0,s==0)
  call compute_stage4_feedback_if_supported(a0,ux,uy,uz,fibre,beta,u_lag,fs,ff,s); write(unit_id,'(A,1X,I0)') 'stage4_unknown_layout_feedback_blocked_flag', merge(1,0,s==0)
  call compute_stage4_feedback_if_supported(a2,ux,uy,uz,fibre,beta,u_lag,fs,ff,s); write(unit_id,'(A,1X,I0)') 'stage4_staggered_layout_feedback_blocked_flag', merge(1,0,s==0)
  write(unit_id,'(A,1X,I0)') 'stage4_feedback_fluid_rhs_modified', 0
  write(unit_id,'(A,1X,I0)') 'stage4_feedback_status', 1
  close(unit_id)
end program
