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
 use iso_c_binding,    only:c_float, c_int, c_bool, c_char
 use asciiutils,       only:fstring
 use filenames,        only:rootname, tagline, nfiles
 use particle_data,    only:dat, maxpart, maxcol, maxstep, npartoftype
 use settings_data,    only:ncolumns, ivegotdata

 implicit none

 public

contains
subroutine check_argcv_f() bind(c)
 include 'libinclude.f90'
end subroutine check_argcv_f

subroutine read_data_c(filename,fileformat,ierr) bind(c, name='read_data')
 character(kind=c_char), intent(in)   :: filename(*), fileformat(*)
 !real(c_double),         intent(out)  :: sph_dat(:,:)
 integer(c_int),         intent(out)  :: ierr

 character(len=120) :: filename_f
 character(len=20)   :: format_f
 integer :: np, nc

 print*, tagline

 ierr = 0

 call defaults_set_initial

 nfiles = 1
 rootname(1) = fstring(filename)
 format_f = fstring(fileformat)

call select_data_format(format_f,ierr)

if (ierr /= 0) then
  call get_data(1,.true.,.true.,1)
  if (ivegotdata .and. maxpart>0) then
    ! np = min(sum(npartoftype(:,1)), size(sph_dat(:,1)) )
    ! nc = min(ncolumns, size(sph_dat(1,:)))
    if (nc > 0) then
       !sph_dat(1:np,1:nc) = dat(1:np,1:nc,1)
    endif
  endif
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
