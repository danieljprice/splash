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
 use getdata,          only:get_data,get_labels
 use labels,           only:label,unitslabel,headertags,shortstring,lenlabel,lenunitslabel
 use params,           only:ltag
 use initialise,       only:defaults_set_initial
 use iso_c_binding,    only:c_float,c_int,c_bool,c_char,c_double
 use asciiutils,       only:fstring,cstring
 use filenames,        only:rootname, tagline, nfiles
 use particle_data,    only:dat,maxpart,maxcol,maxstep,npartoftype,iamtype,headervals
 use settings_data,    only:ncolumns,ivegotdata,required
 use libutils,         only:ctypes_to_fstring,check_argcv

 implicit none

 public

contains
  subroutine check_argcv_c() bind(c, name='check_argcv')
   call check_argcv()
 end subroutine check_argcv_c

subroutine get_labels_c(labels_out, ncol) bind(c, name='get_labels')
  integer(c_int),          intent(in)  :: ncol
  character(kind=c_char),  intent(out) :: labels_out(lenlabel, ncol)
  character(len=lenlabel)    :: temp_string
  character(kind=c_char)     :: temp_cstring
  integer :: i,j

  call get_labels

do i = 1, ncol
  temp_string = shortstring(label(i), unitslabel(i))
  do j = 1, lenlabel
    if (j .le. len(temp_string)) then
      labels_out(j, i) = temp_string(j:j)
    else
      labels_out(j, i) = ' '
    endif
  enddo
end do

  print*, labels_out

end subroutine get_labels_c

subroutine get_header_vals_size(taglength, vallength) bind(c)
  ! Need to get the correct size of the header array to allocate memory in Python
  integer(c_int), intent(out)  :: taglength, vallength
  ! taglength is the size of the tag array
  ! vallength is the length of the headerval array

  taglength = size(headertags)
  vallength = size(headervals)

end subroutine get_header_vals_size

subroutine get_headers(headertags_out, headervals_out, taglength, vallength) bind(c)
  integer(c_int),                  intent(in)  :: taglength, vallength
  character(kind=c_char),          intent(out) :: headertags_out(ltag, taglength)
  real(c_double),                  intent(out) :: headervals_out(vallength)
  integer :: i,j

print*, taglength
  do i = 1, taglength
    do j = 1, ltag
      if (j .le. len(headertags(i))) then
        headertags_out(j, i) = headertags(i)(j:j)
      else
        headertags_out(j, i) = ' '
      endif
    enddo
  enddo

  headervals_out(1:vallength) = headervals(1:vallength, 1)

end subroutine get_headers

subroutine read_data_c(filename,fileformat,f_length, ff_length,&
                       sph_dat,npart,ncol,read_header,verbose,ierr) bind(c, name='read_data')
 integer(c_int),         intent(in)     :: f_length, ff_length
 character(kind=c_char), intent(in)     :: filename(f_length), fileformat(ff_length)
 integer(c_int),         intent(inout)  :: ncol, npart
 real(c_double),         intent(out)    :: sph_dat(npart,ncol)
 integer(c_int),         intent(in)     :: read_header, verbose
 integer(c_int),         intent(out)    :: ierr

 character(len=ff_length)   :: format_f

 integer   :: ncolr,npartr

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
      npartr = sum(npartoftype(:,1))
      ncolr = ncolumns
      if (ncol/=0 .or. npart/=0) then
        if (ncol/=ncolr+1 .or. npart/=npartr) then
          print*, "WARNING: Array size given in libread is not equal to read array."
        endif
      else if (ncol==0 .or. npart==0) then
        ncol = ncolr
        npart = npartr
    endif

   if (ncol > 0 .and. read_header/=1) then
       sph_dat(1:npart,1:ncol-1) = dat(1:npart,1:ncol-1,1)
       sph_dat(1:npart,ncol) = iamtype(1:npart,1)
    else
        if (verbose==1) print*, "Updating values for npart and ncol."
        npart = sum(npartoftype(:,1))
        ncol = ncolumns + 1
     endif
   endif
 else
   print*, "Error in selecting the data format"
 endif

end subroutine read_data_c

end module libreaddata
