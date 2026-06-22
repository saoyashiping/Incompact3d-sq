program fibre_prod_p1_oneway_channel_case_check
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fibre_prod_p1_oneway_channel_case
  implicit none
  real(real64) :: ux(4,4,4), uy(4,4,4), uz(4,4,4)
  integer :: status
  ux = 1.0_real64; uy = 0.0_real64; uz = 0.0_real64
  call fibre_prod_p1_oneway_channel_case_init_values(1,49,0.0_real64,0.10_real64, &
       1.0e-2_real64,1.0e2_real64,status); if (status /= 0) error stop 1
  call fibre_prod_p1_oneway_channel_case_record_sampling(ux,uy,uz,status); if (status /= 0) error stop 2
  call fibre_prod_p1_oneway_channel_case_record_structure_response(status); if (status /= 0) error stop 3
  call fibre_prod_p1_oneway_channel_case_write_diagnostics(status); if (status /= 0) error stop 4
  ux(2,2,2) = ieee_value(ux(2,2,2), ieee_quiet_nan)
  call fibre_prod_p1_oneway_channel_case_record_sampling(ux,uy,uz,status); if (status == 0) error stop 5
  ux = 1.0_real64
  call fibre_prod_p1_oneway_channel_case_init_values(1,49,0.0_real64,0.10_real64, &
       1.0e-2_real64,1.0e-3_real64,status); if (status /= 0) error stop 6
  call fibre_prod_p1_oneway_channel_case_record_sampling(ux,uy,uz,status); if (status /= 0) error stop 7
  call fibre_prod_p1_oneway_channel_case_record_structure_response(status); if (status == 0) error stop 8
  call fibre_prod_p1_oneway_channel_case_init_values(1,49,0.0_real64,0.10_real64, &
       1.0e-2_real64,1.0e2_real64,status); if (status /= 0) error stop 9
  call fibre_prod_p1_oneway_channel_case_set_sample_for_check( &
       ieee_value(ux(1,1,1), ieee_quiet_nan),0.0_real64,0.0_real64)
  call fibre_prod_p1_oneway_channel_case_record_structure_response(status); if (status == 0) error stop 14
  call fibre_prod_p1_oneway_channel_case_init_values(1,49,1.0_real64,0.10_real64, &
       1.0e-2_real64,1.0e2_real64,status); if (status == 0) error stop 13
  call fibre_prod_p1_oneway_channel_case_init_values(2,49,0.0_real64,0.10_real64, &
       1.0e-2_real64,1.0e2_real64,status); if (status == 0) error stop 10
  call fibre_prod_p1_oneway_channel_case_init_values(1,48,0.0_real64,0.10_real64, &
       1.0e-2_real64,1.0e2_real64,status); if (status == 0) error stop 11
  call fibre_prod_p1_oneway_channel_case_init_values(1,49,0.0_real64,0.60_real64, &
       1.0e-2_real64,1.0e2_real64,status); if (status == 0) error stop 12
  write(*,'(A)') 'P1_1_ONEWAY_CHANNEL_CASE_CHECK PASS'
end program fibre_prod_p1_oneway_channel_case_check
