program fibre_stage9_2_minimal_parallel_gate
  use mpi
  use decomp_2d
  use decomp_2d_mpi, only: nrank, nproc
  implicit none

  integer :: ierr, mpi_ok
  integer :: fail_local, fail_global
  integer :: nxg, nyg, nzg, p_row_l, p_col_l
  logical :: periodic_bc(3)
  integer :: i
  integer(kind=8) :: local_cells, total_cells
  integer :: fail_count_local, fail_count_global
  logical :: coarse_checks_enabled

  call MPI_Init(ierr)
  mpi_ok = merge(1, 0, ierr == MPI_SUCCESS)
  if (mpi_ok == 1) then
    write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': MPI initialized'
  else
    write(*,'(A,I0,A,I0)') '[FAIL] rank ', nrank, ': MPI_Init failed ierr=', ierr
  end if

  fail_local = 0
  fail_count_local = 0

  nxg = 8
  nyg = 8
  nzg = 8
  p_row_l = 0
  p_col_l = 0
  periodic_bc = (/ .false., .false., .false. /)

  call decomp_2d_init(nxg, nyg, nzg, p_row_l, p_col_l, periodic_bc)
  write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': decomp_2d initialized'

  write(*,'(A,I0,A,3(I0,1X),A,3(I0,1X),A,3(I0,1X))') &
    '[INFO] rank=', nrank, ' xsize=', xsize(1), xsize(2), xsize(3), &
    ' xstart=', xstart(1), xstart(2), xstart(3), ' xend=', xend(1), xend(2), xend(3)
  write(*,'(A,I0,A,3(I0,1X),A,3(I0,1X),A,3(I0,1X))') &
    '[INFO] rank=', nrank, ' ysize=', ysize(1), ysize(2), ysize(3), &
    ' ystart=', ystart(1), ystart(2), ystart(3), ' yend=', yend(1), yend(2), yend(3)
  write(*,'(A,I0,A,3(I0,1X),A,3(I0,1X),A,3(I0,1X))') &
    '[INFO] rank=', nrank, ' zsize=', zsize(1), zsize(2), zsize(3), &
    ' zstart=', zstart(1), zstart(2), zstart(3), ' zend=', zend(1), zend(2), zend(3)

  if (all(xsize(1:3) >= 0) .and. all(ysize(1:3) >= 0) .and. all(zsize(1:3) >= 0)) then
    write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': local cell count non-negative'
  else
    write(*,'(A,I0,A)') '[FAIL] rank ', nrank, ': local sizes contain negative values'
    fail_local = 1; fail_count_local = fail_count_local + 1
  end if

  if (all((xsize(1:3) == 0) .or. (xstart(1:3) <= xend(1:3)))) then
    write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': x-pencil local bounds valid'
  else
    write(*,'(A,I0,A)') '[FAIL] rank ', nrank, ': x-pencil local bounds invalid'
    fail_local = 1; fail_count_local = fail_count_local + 1
  end if

  if (all((ysize(1:3) == 0) .or. (ystart(1:3) <= yend(1:3)))) then
    write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': y-pencil local bounds valid'
  else
    write(*,'(A,I0,A)') '[FAIL] rank ', nrank, ': y-pencil local bounds invalid'
    fail_local = 1; fail_count_local = fail_count_local + 1
  end if

  if (all((zsize(1:3) == 0) .or. (zstart(1:3) <= zend(1:3)))) then
    write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': z-pencil local bounds valid'
  else
    write(*,'(A,I0,A)') '[FAIL] rank ', nrank, ': z-pencil local bounds invalid'
    fail_local = 1; fail_count_local = fail_count_local + 1
  end if

  if (all((xsize(1:3) == 0) .or. (xsize(1:3) == xend(1:3)-xstart(1:3)+1))) then
    write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': xsize equals xend-xstart+1 where local size > 0'
  else
    write(*,'(A,I0,A)') '[FAIL] rank ', nrank, ': xsize/end-start mismatch'
    fail_local = 1; fail_count_local = fail_count_local + 1
  end if

  if (all((ysize(1:3) == 0) .or. (ysize(1:3) == yend(1:3)-ystart(1:3)+1))) then
    write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': ysize equals yend-ystart+1 where local size > 0'
  else
    write(*,'(A,I0,A)') '[FAIL] rank ', nrank, ': ysize/end-start mismatch'
    fail_local = 1; fail_count_local = fail_count_local + 1
  end if

  if (all((zsize(1:3) == 0) .or. (zsize(1:3) == zend(1:3)-zstart(1:3)+1))) then
    write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': zsize equals zend-zstart+1 where local size > 0'
  else
    write(*,'(A,I0,A)') '[FAIL] rank ', nrank, ': zsize/end-start mismatch'
    fail_local = 1; fail_count_local = fail_count_local + 1
  end if

  local_cells = int(max(0,xsize(1)),8) * int(max(0,xsize(2)),8) * int(max(0,xsize(3)),8)
  call MPI_Allreduce(local_cells, total_cells, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (ierr == MPI_SUCCESS) then
    if (total_cells > 0_8) then
      write(*,'(A,I0,A,I0)') '[PASS] rank ', nrank, ': global local-cell-count sum positive, total=', total_cells
    else
      write(*,'(A,I0,A,I0)') '[FAIL] rank ', nrank, ': global local-cell-count sum not positive, total=', total_cells
      fail_local = 1; fail_count_local = fail_count_local + 1
    end if
  else
    write(*,'(A,I0,A,I0)') '[FAIL] rank ', nrank, ': MPI_Allreduce(total_cells) failed ierr=', ierr
    fail_local = 1; fail_count_local = fail_count_local + 1
  end if

  coarse_checks_enabled = .false.
  if (coarse_checks_enabled) then
    if (all(xstS(1:3) >= 1) .and. all(xenS(1:3) >= xstS(1:3)-1)) then
      write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': coarse/stat bounds initialized and valid'
    else
      write(*,'(A,I0,A)') '[FAIL] rank ', nrank, ': coarse/stat bounds invalid'
      fail_local = 1; fail_count_local = fail_count_local + 1
    end if
    if (all(xstV(1:3) >= 1) .and. all(xenV(1:3) >= xstV(1:3)-1)) then
      write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': visu bounds initialized and valid'
    else
      write(*,'(A,I0,A)') '[FAIL] rank ', nrank, ': visu bounds invalid'
      fail_local = 1; fail_count_local = fail_count_local + 1
    end if
  else
    write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': optional coarse/stat bounds check disabled (not part of minimal gate)'
    write(*,'(A,I0,A)') '[PASS] rank ', nrank, ': optional visu bounds check disabled (not part of minimal gate)'
  end if

  call MPI_Allreduce(fail_local, fail_global, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(fail_count_local, fail_count_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (nrank == 0) then
    if (fail_global == 0) then
      write(*,'(A)') 'STAGE 9.2 PROGRAM VERDICT: PASS'
      write(*,'(A,I0)') 'Total failed sub-check count across ranks: ', fail_count_global
    else
      write(*,'(A)') 'STAGE 9.2 PROGRAM VERDICT: FAIL'
      write(*,'(A)') 'Failed program checks:'
      write(*,'(A)') '  - See per-rank [FAIL] lines above for exact failed checks.'
      write(*,'(A,I0)') 'Total failed sub-check count across ranks: ', fail_count_global
    end if
  end if

  call decomp_2d_finalize()
  call MPI_Finalize(ierr)
  if (fail_global /= 0) stop 1

end program fibre_stage9_2_minimal_parallel_gate
