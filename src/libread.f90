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
 use settings_data,    only:ncolumns, ivegotdata, required
 use libutils,         only:ctypes_to_fstring, check_argcv

 implicit none

 public

contains
  subroutine check_argcv_c() bind(c, name='check_argcv')
   call check_argcv()
 end subroutine check_argcv_c

subroutine read_data_c(filename,fileformat,f_length, ff_length,&
                       sph_dat,npart,ncol,read_header,verbose,ierr) bind(c, name='read_data')
 integer(c_int),         intent(in)     :: f_length, ff_length
 character(kind=c_char), intent(in)     :: filename(f_length), fileformat(ff_length)
 integer(c_int),         intent(inout)  :: ncol, npart
 real(c_double),         intent(out)    :: sph_dat(npart,ncol)
 integer(c_int),         intent(in)     :: read_header, verbose
 integer(c_int),         intent(out)    :: ierr

 character(len=f_length)    :: filename_f
 character(len=ff_length)   :: format_f

 integer   :: i,j

 if (verbose==1) print*, tagline

 ierr = 0

 call defaults_set_initial

 if (read_header==1) required = .false.

 nfiles = 1
 rootname(1) = ctypes_to_fstring(filename)
 format_f = ctypes_to_fstring(fileformat)

 if (verbose==1) then
   print*, "Received file format f ", format_f
   print*, "Size of sph_dat is ", size(sph_dat)
 endif

call select_data_format(format_f,ierr)

 if (ierr == 0) then
   if (verbose==1) print*, "Calling get_data"
   call get_data(1,.true.,.true.,1)
   if (ivegotdata .and. maxpart>0) then
     npart = min(sum(npartoftype(:,1)), size(sph_dat(:,1)) )
     ncol = min(ncolumns, size(sph_dat(1,:)))

   if (ncol > 0) then
       sph_dat(1:npart,1:ncol) = dat(1:npart,1:ncol,1)
    else
        if (verbose==1) print*, "Updating values for npart and ncol."
        npart = sum(npartoftype(:,1))
        ncol = ncolumns
     endif
   endif
 else
   print*, "Error in selecting the data format"
 endif

end subroutine read_data_c

end module libreaddata
