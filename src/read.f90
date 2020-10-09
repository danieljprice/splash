module read_test
 use readdata,         only:select_data_format
 use getdata,          only:get_data, get_labels
 use initialise,       only:defaults_set_initial
 use iso_c_binding,    only:c_float, c_int, c_bool, c_char, c_double
 use asciiutils,       only:fstring
 use filenames,        only:rootname, tagline, nfiles
 use particle_data,    only:dat, maxpart, maxcol, maxstep, npartoftype
 use settings_data,    only:ncolumns, ivegotdata
 use libutils,         only:ctypes_to_fstring

 implicit none

 public

contains

  subroutine read_data_wrap(filename,fileformat,&
                           sph_dat,npart,ncol,ierr)
   character(len=*),       intent(in)     :: filename, fileformat
   integer,         intent(inout)  :: ncol, npart
   real(8),         intent(out)    :: sph_dat(:,:)
   integer,         intent(out)    :: ierr

   integer   :: i,j

   print*, tagline

   ierr = 0

   call defaults_set_initial

   nfiles = 1
   rootname(1) = filename

   print*, "format is ", fileformat

  call select_data_format(fileformat,ierr)

   if (ierr == 0) then
     call get_data(1,.true.,.true.,1)
     if (ivegotdata .and. maxpart>0) then
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

end subroutine read_data_wrap

end module read_test
