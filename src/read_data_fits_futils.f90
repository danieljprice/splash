module fitsutils
 use, intrinsic :: iso_c_binding, only:c_int,c_float,c_char
 !interface to the c versions
 interface
  subroutine read_fits_header(filename,npart,ncol,ierr) bind(c)
   import
   character(kind=c_char,len=1), intent(in) :: filename
   integer(kind=c_int), intent(out) :: npart,ncol,ierr
  end subroutine read_fits_header
 end interface

 interface
  subroutine read_fits_data(ierr) bind(c)
   import
   integer(kind=c_int), intent(out) :: ierr
  end subroutine read_fits_data
 end interface

end module fitsutils
