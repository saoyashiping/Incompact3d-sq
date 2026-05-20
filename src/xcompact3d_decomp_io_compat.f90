module xcompact3d_decomp_io_compat
  use decomp_2d_constants, only : mytype
  use decomp_2d, only : DECOMP_INFO
  use decomp_2d_io, only : &
       decomp_2d_write_mode, decomp_2d_read_mode, decomp_2d_append_mode, &
       real_decomp_2d_open_io   => decomp_2d_open_io, &
       real_decomp_2d_close_io  => decomp_2d_close_io, &
       real_decomp_2d_start_io  => decomp_2d_start_io, &
       real_decomp_2d_end_io    => decomp_2d_end_io, &
       real_decomp_2d_write_one => decomp_2d_write_one, &
       real_decomp_2d_read_one  => decomp_2d_read_one, &
       real_decomp_2d_write_plane => decomp_2d_write_plane, &
       fine_to_coarseS, fine_to_coarseV, gen_iodir_name

  implicit none
  private

  public :: decomp_2d_write_mode, decomp_2d_read_mode, decomp_2d_append_mode
  public :: decomp_2d_init_io, decomp_2d_register_variable
  public :: decomp_2d_open_io, decomp_2d_close_io
  public :: decomp_2d_start_io, decomp_2d_end_io
  public :: decomp_2d_write_one, decomp_2d_read_one
  public :: decomp_2d_write_plane
  public :: fine_to_coarseS, fine_to_coarseV
  public :: gen_iodir_name

  interface decomp_2d_write_one
    module procedure x3d_write_one_r3_simple
    module procedure x3d_write_one_r3_legacy
  end interface

  interface decomp_2d_read_one
    module procedure x3d_read_one_r3_simple
    module procedure x3d_read_one_r3_legacy
  end interface

  interface decomp_2d_write_plane
    module procedure x3d_write_plane_r3_simple
    module procedure x3d_write_plane_r3_legacy
  end interface

contains
  subroutine decomp_2d_init_io(io_name)
    character(len=*), intent(in) :: io_name
  end subroutine

  subroutine decomp_2d_register_variable(io_name, varname, ipencil, iscalar, output2D, dtype, opt_decomp, opt_nplanes)
    character(len=*), intent(in) :: io_name, varname
    integer, intent(in) :: ipencil, iscalar, output2D
    integer, intent(in) :: dtype
    type(DECOMP_INFO), intent(in), optional :: opt_decomp
    integer, intent(in), optional :: opt_nplanes
  end subroutine

  subroutine decomp_2d_open_io(io_name, dirname, mode)
    character(len=*), intent(in) :: io_name, dirname
    integer, intent(in) :: mode
    call real_decomp_2d_open_io(io_name, dirname, mode)
  end subroutine
  subroutine decomp_2d_close_io(io_name, dirname)
    character(len=*), intent(in) :: io_name, dirname
    call real_decomp_2d_close_io(io_name, dirname)
  end subroutine
  subroutine decomp_2d_start_io(io_name, dirname)
    character(len=*), intent(in) :: io_name, dirname
    call real_decomp_2d_start_io(io_name, dirname)
  end subroutine
  subroutine decomp_2d_end_io(io_name, dirname)
    character(len=*), intent(in) :: io_name, dirname
    call real_decomp_2d_end_io(io_name, dirname)
  end subroutine

  subroutine x3d_write_one_r3_simple(ipencil, array, filename, mode, opt_decomp, reduce_prec, opt_deferred_writes)
    integer, intent(in) :: ipencil
    real(mytype), intent(in) :: array(:,:,:)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mode
    type(DECOMP_INFO), intent(in), optional :: opt_decomp
    logical, intent(in), optional :: reduce_prec
    logical, intent(in), optional :: opt_deferred_writes
    call real_decomp_2d_write_one(ipencil, array, filename, mode)
  end subroutine

  subroutine x3d_write_one_r3_legacy(ipencil, array, dirname, filename, mode, io_name, opt_decomp, reduce_prec, opt_deferred_writes)
    integer, intent(in) :: ipencil
    real(mytype), intent(in) :: array(:,:,:)
    character(len=*), intent(in) :: dirname, filename
    integer, intent(in) :: mode
    character(len=*), intent(in), optional :: io_name
    type(DECOMP_INFO), intent(in), optional :: opt_decomp
    logical, intent(in), optional :: reduce_prec
    logical, intent(in), optional :: opt_deferred_writes
    character(len=:), allocatable :: fullpath
    fullpath = trim(dirname)//'/'//trim(filename)
    call x3d_write_one_r3_simple(ipencil, array, fullpath, mode, opt_decomp, reduce_prec, opt_deferred_writes)
  end subroutine

  subroutine x3d_read_one_r3_simple(ipencil, array, filename, opt_decomp, reduce_prec)
    integer, intent(in) :: ipencil
    real(mytype), intent(inout) :: array(:,:,:)
    character(len=*), intent(in) :: filename
    type(DECOMP_INFO), intent(in), optional :: opt_decomp
    logical, intent(in), optional :: reduce_prec
    call real_decomp_2d_read_one(ipencil, array, filename)
  end subroutine

  subroutine x3d_read_one_r3_legacy(ipencil, array, dirname, filename, io_name, opt_decomp, reduce_prec)
    integer, intent(in) :: ipencil
    real(mytype), intent(inout) :: array(:,:,:)
    character(len=*), intent(in) :: dirname, filename
    character(len=*), intent(in), optional :: io_name
    type(DECOMP_INFO), intent(in), optional :: opt_decomp
    logical, intent(in), optional :: reduce_prec
    character(len=:), allocatable :: fullpath
    fullpath = trim(dirname)//'/'//trim(filename)
    call x3d_read_one_r3_simple(ipencil, array, fullpath, opt_decomp, reduce_prec)
  end subroutine


  subroutine x3d_write_plane_r3_simple(ipencil, array, output2D, plane_index, filename)
    integer, intent(in) :: ipencil
    real(mytype), intent(in) :: array(:,:,:)
    integer, intent(in) :: output2D
    integer, intent(in) :: plane_index
    character(len=*), intent(in) :: filename

    call real_decomp_2d_write_plane(ipencil, array, output2D, plane_index, filename)
  end subroutine

  subroutine x3d_write_plane_r3_legacy(ipencil, array, output2D, plane_index, dirname, filename, io_name)
    integer, intent(in) :: ipencil
    real(mytype), intent(in) :: array(:,:,:)
    integer, intent(in) :: output2D
    integer, intent(in) :: plane_index
    character(len=*), intent(in) :: dirname, filename
    character(len=*), intent(in), optional :: io_name
    character(len=:), allocatable :: fullpath

    fullpath = trim(dirname)//'/'//trim(filename)
    call x3d_write_plane_r3_simple(ipencil, array, output2D, plane_index, fullpath)
  end subroutine

end module
