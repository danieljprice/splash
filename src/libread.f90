!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2020 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! Module providing library version of splash read_data routines
! specifies c interfaces to corresponding Fortran subroutines
!-------------------------------------------------------------------------
module libreaddata

 use readdata,         only:select_data_format
 use getdata,          only:get_data, get_labels
 use initialise,       only:defaults_set_initial
 use iso_c_binding,    only:c_float, c_int, c_bool, c_char, c_double
 use asciiutils,       only:fstring
 use filenames,        only:rootname, tagline, nfiles
 use particle_data,    only:dat, maxpart, maxcol, maxstep, npartoftype
 use settings_data,    only:ncolumns, ivegotdata
 use libutils,         only:ctypes_to_fstring, check_argcv

 implicit none

 public

contains
  subroutine check_argcv_c() bind(c, name='check_argcv')
   call check_argcv()
 end subroutine check_argcv_c

subroutine read_data_c(filename,fileformat,f_length, ff_length,&
                       sph_dat,npart,ncol,ierr) bind(c, name='read_data')
 integer(c_int),         intent(in)     :: f_length, ff_length
 character(kind=c_char), intent(in)     :: filename(f_length), fileformat(ff_length)
 integer(c_int),         intent(inout)  :: ncol, npart
 real(c_double),         intent(out)    :: sph_dat(ncol,npart)
 integer(c_int),         intent(out)    :: ierr

 character(len=f_length) :: filename_f
 character(len=ff_length)   :: format_f

 integer   :: i,j

 print*, tagline

 ierr = 0

 call defaults_set_initial

 nfiles = 1
 rootname(1) = ctypes_to_fstring(filename)
 format_f = ctypes_to_fstring(fileformat)

 print*, "format_f is ", format_f

call select_data_format(format_f,ierr)
print*, "ierr is:", ierr

 if (ierr == 0) then
   print*, "calling get_data"
   call get_data(1,.true.,.true.,1)
   print*, "called get_data"
   if (ivegotdata .and. maxpart>0) then
     print*, "setting ncol and npart"
     print*, "sum(nartpartoftype(:,1))", sum(npartoftype(:,1))
     print*, "size(sph_dat(:,1))", size(sph_dat(:,1))
     print*, "ncolumns", ncolumns
    npart = min(sum(npartoftype(:,1)), size(sph_dat(:,1)) )
    ncol = min(ncolumns, size(sph_dat(1,:)))
    print*, "npart and ncol in fortran are", npart, ncol
     if (ncol > 0) then
       do i=1,ncol
         do j=1,npart
           print*, "trying to write to sph_dat"
           print*, "sph_dat(i,j) is", sph_dat(i,j)
           print*, "dat(i,j,1) is", dat(i,j,1)
           sph_dat(i,j) = dat(i,j,1)
           print*,""
         end do
      end do
       ! print*, "attempting to write to sph_dat"
       ! print*, "dat(1:1,1:1,1)", dat(1:1,1:1,1)
       ! print*, "sph_dat(1:1,1:1)", sph_dat(1:1,1:1)
       ! sph_dat(1:npart,1:ncol) = dat(1:npart,1:ncol,1)
     endif
   endif
 else
   print*, "Error in selecting the data format"
 endif

end subroutine read_data_c

! subroutine select_data_format_c(string, ierr) bind(c, name='select_data_format')
!
!  integer(c_int), intent(out)   :: ierr
!
!  call select_data_format(string, ierr)
!
! end subroutine select_data_format_c

end module libreaddata
