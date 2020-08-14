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
!
! module containing routines to read/write images in .pfm format
! Info below taken from: http://www.pauldebevec.com/Research/HDR/PFM/
!
! The Portable FloatMap format was designed as a floating-point image format.
! The format begins with three lines of text specifying the image size and type,
!  and then continues with raw binary image data for the rest of the file.
! The text header of a .pfm file takes the following form:
! [type]
! [xres] [yres]
! [byte_order]
!
! Each of the three lines of text ends with a 1-byte Unix-style carriage return.
! "[type]" is one of "PF" for a 3-channel RGB color image, or "Pf" for a monochrome single-channel image. 
! "[xres] [yres]" indicates the x and y resolutions of the image. 
! "[byte_order]" is a number used to indicate the byte order within the file. 
! A positive number (e.g. "1.0") indicates big-endian, with the most significant byte of each 4-byte float first.
! If the number is negative (e.g. "-1.0") this indicates little-endian, with the least significant byte first. There are no comments in these files. For example:
!
! PF
! 768 512
! -1.0
! Indicates an RGB color image, 768 by 512 pixels, with little-endian byte order.
! After the final carriage return the file proceeds with a series of three 4-byte IEEE 754 single
! precision floating point numbers for each pixel, specified in left to right, bottom to top order.
!
!-----------------------------------------------------------------
module write_pfm
 implicit none

contains

subroutine write_pixmap_pfm(filename,nx,ny,datpix,ierr)
 use iso_c_binding, only:c_float
 character(len=*), intent(in)  :: filename
 integer,          intent(in)  :: nx,ny
 real(c_float),    intent(in)  :: datpix(nx,ny)
 integer,          intent(out) :: ierr
 integer :: iunit
 logical :: bigendian

 bigendian = IACHAR(TRANSFER(1,"a")) == 0

 ! write formatted header
 open(newunit=iunit,file=filename,action='write',status='replace',iostat=ierr)
 if (ierr /= 0) return
 write(iunit,"(a)",iostat=ierr) 'Pf'
 write(iunit,*,iostat=ierr) nx,ny
 if (bigendian) then
    write(iunit,*,iostat=ierr) 1.0
 else
    write(iunit,*,iostat=ierr) -1.0
 endif
 close(iunit)

 ! write unformatted data
 open(newunit=iunit,file=filename,action='write',status='old',position='append',form='unformatted',iostat=ierr)
 write(iunit) datpix(1:nx,1:ny)
 close(iunit)

end subroutine write_pixmap_pfm

end module write_pfm