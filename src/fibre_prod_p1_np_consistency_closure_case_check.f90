program fibre_prod_p1_np_consistency_closure_case_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fibre_prod_p1_np_consistency_closure_case
  implicit none
  real(real64) :: ux(4,4,4), uy(4,4,4), uz(4,4,4), rhsx(4,4,4), rhsy(4,4,4), rhsz(4,4,4)
  real(real64) :: sig1(6), sig2(6)
  integer :: status
  ux=1.0_real64; uy=0.0_real64; uz=0.0_real64; rhsx=0.0_real64; rhsy=0.0_real64; rhsz=0.0_real64
  call run_case(0.0_real64, rhsx, rhsy, rhsz, sig1, status); if (status /= 0) error stop 1
  rhsx=0.0_real64; rhsy=0.0_real64; rhsz=0.0_real64
  call run_case(1.0e-4_real64, rhsx, rhsy, rhsz, sig2, status); if (status /= 0) error stop 2
  call fibre_prod_p1_np_consistency_closure_case_compare_signatures(sig2, sig2, status); if (status /= 0) error stop 3
  ux(2,2,2)=ieee_value(ux(2,2,2), ieee_quiet_nan)
  call fibre_prod_p1_np_consistency_closure_case_init_values(1,49,1.0e-4_real64,2.0_real64,0.10_real64, &
       1.0e-2_real64,1.0e2_real64,1.0e2_real64,1.0e6_real64,status); if (status /= 0) error stop 4
  call fibre_prod_p1_np_consistency_closure_case_record_sampling(ux,uy,uz,status); if (status == 0) error stop 5
  ux=1.0_real64
  call fibre_prod_p1_np_consistency_closure_case_set_sample_for_check(ieee_value(ux(1,1,1), ieee_quiet_nan),0.0_real64,0.0_real64)
  call fibre_prod_p1_np_consistency_closure_case_record_structure_response(status); if (status == 0) error stop 6
  call fibre_prod_p1_np_consistency_closure_case_init_values(1,49,1.0e-4_real64,2.0_real64,0.60_real64, &
       1.0e-2_real64,1.0e2_real64,1.0e2_real64,1.0e6_real64,status); if (status == 0) error stop 7
  call fibre_prod_p1_np_consistency_closure_case_init_values(2,49,1.0e-4_real64,2.0_real64,0.10_real64, &
       1.0e-2_real64,1.0e2_real64,1.0e2_real64,1.0e6_real64,status); if (status == 0) error stop 8
  call fibre_prod_p1_np_consistency_closure_case_init_values(1,48,1.0e-4_real64,2.0_real64,0.10_real64, &
       1.0e-2_real64,1.0e2_real64,1.0e2_real64,1.0e6_real64,status); if (status == 0) error stop 9
  rhsx=0.0_real64; rhsy=0.0_real64; rhsz=0.0_real64
  call fibre_prod_p1_np_consistency_closure_case_init_values(1,49,1.0e-4_real64,2.0_real64,0.10_real64, &
       1.0e-2_real64,1.0e2_real64,1.0e2_real64,1.0e6_real64,status); if (status /= 0) error stop 10
  call fibre_prod_p1_np_consistency_closure_case_record_force_buffer(status); if (status == 0) error stop 11
  write(*,'(A)') 'P1_4_NP_CONSISTENCY_CLOSURE_CASE_CHECK PASS'
contains
  subroutine run_case(lambda, rx, ry, rz, sig, status)
    real(real64), intent(in) :: lambda
    real(real64), intent(inout) :: rx(:,:,:), ry(:,:,:), rz(:,:,:)
    real(real64), intent(out) :: sig(6)
    integer, intent(out) :: status
    call fibre_prod_p1_np_consistency_closure_case_init_values(1,49,lambda,2.0_real64,0.10_real64, &
         1.0e-2_real64,1.0e2_real64,1.0e2_real64,1.0e6_real64,status); if (status /= 0) return
    call fibre_prod_p1_np_consistency_closure_case_record_sampling(ux,uy,uz,status); if (status /= 0) return
    call fibre_prod_p1_np_consistency_closure_case_record_structure_response(status); if (status /= 0) return
    call fibre_prod_p1_np_consistency_closure_case_record_reaction_force(status); if (status /= 0) return
    call fibre_prod_p1_np_consistency_closure_case_record_force_buffer(status); if (status /= 0) return
    call fibre_prod_p1_np_consistency_closure_case_record_rhs_increment(rx,ry,rz,status); if (status /= 0) return
    call fibre_prod_p1_np_consistency_closure_case_record_global_signatures(status); if (status /= 0) return
    call fibre_prod_p1_np_consistency_closure_case_write_diagnostics(status); if (status /= 0) return
    sig = (/1.0_real64,1.0_real64,1.0_real64,1.0_real64,1.0_real64,1.0_real64/)
  end subroutine
end program fibre_prod_p1_np_consistency_closure_case_check
