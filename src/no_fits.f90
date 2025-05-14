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
!  Copyright (C) 2020- Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
!----------------------------------------------------------------------
!
!  Module handling read and write of FITS files
!  With thanks to Christophe Pinte
!
!----------------------------------------------------------------------
module readwrite_fits
 use iso_fortran_env, only:real32,real64
 implicit none
 public :: read_fits_image,write_fits_image,fits_error
 public :: read_fits_cube,write_fits_cube
 public :: read_fits_header
 public :: get_floats_from_fits_header,get_from_header,get_from_header_s
 public :: append_to_fits_cube
 public :: get_velocity_from_fits_header
 public :: flatten_header

 interface write_fits_image
  module procedure write_fits_image,write_fits_image64
 end interface write_fits_image

 interface write_fits_cube
  module procedure write_fits_cube,write_fits_cube64
 end interface write_fits_cube

 interface read_fits_image
  module procedure read_fits_image,read_fits_image64
 end interface read_fits_image

 interface read_fits_cube
  module procedure read_fits_cube,read_fits_cube64
 end interface read_fits_cube

 interface append_to_fits_cube
  module procedure append_to_fits_cube
 end interface append_to_fits_cube

 private

contains

!---------------------------------------------------
! subroutine to read image from FITS file
! using cfitsio library
!---------------------------------------------------
subroutine read_fits_image(filename,image,naxes,ierr,hdr)
 character(len=*), intent(in)   :: filename
 real(kind=real32), intent(out), allocatable :: image(:,:)
 character(len=80), intent(inout), allocatable, optional :: hdr(:)
 integer, intent(out) :: naxes(2),ierr

 ierr = 1
 naxes = 0
end subroutine read_fits_image

!---------------------------------------------------
! read FITS header from file
!---------------------------------------------------
subroutine read_fits_header(filename,hdr,ierr)
 character, intent(in)  :: filename
 character(len=80), allocatable, intent(out) :: hdr(:)
 integer, intent(out) :: ierr

 ierr = 1

end subroutine read_fits_header

!---------------------------------------------------
! internal subroutine to read FITS header information
!---------------------------------------------------
subroutine read_fits_head(iunit,hdr,ierr)
 integer, intent(in)  :: iunit
 integer, intent(out) :: ierr
 character(len=80), allocatable, intent(inout) :: hdr(:)

 ierr = 1

end subroutine read_fits_head

!---------------------------------------------------
! internal subroutine to write FITS header information
! excluding things we have changed
!---------------------------------------------------
subroutine write_fits_head(iunit,hdr,ierr)
 integer, intent(in) :: iunit
 character(len=80), intent(in) :: hdr(:)
 integer, intent(out) :: ierr

 ierr = 1

end subroutine write_fits_head

!---------------------------------------------------
! subroutine to read spectral cube from FITS file
! using cfitsio library
!---------------------------------------------------
subroutine read_fits_cube(filename,image,naxes,ierr,hdr,hdu,velocity)
 character(len=*), intent(in)   :: filename
 real(kind=real32), intent(out), allocatable :: image(:,:,:)
 character(len=80), intent(inout), allocatable, optional :: hdr(:)
 real(kind=real32), intent(out), allocatable, optional :: velocity(:)
 integer, intent(out) :: naxes(4),ierr
 integer, intent(in), optional :: hdu ! specify which hdu to read

 ierr = 1
 naxes = 0

end subroutine read_fits_cube

!---------------------------------------------------
! error code handling
!---------------------------------------------------
 character(len=50) function fits_error(ierr)
  integer, intent(in) :: ierr

  select case(ierr)
  case(1)
     fits_error = 'not compiled with fits lib'
  case default
     fits_error = 'unknown error'
  end select

 end function fits_error

!------------------------------------------------
! Writing new fits file
!------------------------------------------------
 subroutine write_fits_image(filename,image,naxes,ierr,hdr)
  character(len=*), intent(in) :: filename
  integer, intent(in)  :: naxes(2)
  real(kind=real32), intent(in) :: image(naxes(1),naxes(2))
  integer, intent(out) :: ierr
  character(len=80), intent(in), optional :: hdr(:)

  ierr = 1

 end subroutine write_fits_image

!-------------------------------------------------------------
! Writing new fits file (convert from double precision input)
!-------------------------------------------------------------
 subroutine write_fits_image64(filename,image,naxes,ierr,hdr)
  character(len=*), intent(in) :: filename
  integer,          intent(in) :: naxes(2)
  real(kind=real64),intent(in) :: image(naxes(1),naxes(2))
  real(kind=real32), allocatable :: img32(:,:)
  integer, intent(out) :: ierr
  character(len=80), intent(in), optional :: hdr(:)

  ierr = 1

 end subroutine write_fits_image64

!------------------------------------------------
! Writing new fits file
!------------------------------------------------
 subroutine write_fits_cube(filename,image,naxes,ierr,hdr)
   character(len=*), intent(in) :: filename
   integer, intent(in)  :: naxes(3)
   real(kind=real32), intent(in) :: image(naxes(1),naxes(2),naxes(3))
   integer, intent(out) :: ierr
   character(len=80), intent(in), optional :: hdr(:)

   ierr = 1

 end subroutine write_fits_cube


!------------------------------------------------
! Writing new fits file
!------------------------------------------------
 subroutine append_to_fits_cube(filename,image,naxes,ierr,hdr)
   character(len=*), intent(in) :: filename
   integer, intent(in)  :: naxes(3)
   real(kind=real32), intent(in) :: image(naxes(1),naxes(2),naxes(3))
   integer, intent(out) :: ierr
   character(len=80), intent(in), optional :: hdr(:)

   ierr = 1

 end subroutine append_to_fits_cube

!-------------------------------------------------------------
! Writing new fits file (convert from double precision input)
!-------------------------------------------------------------
subroutine write_fits_cube64(filename,image,naxes,ierr,hdr)
 character(len=*), intent(in) :: filename
 integer,          intent(in) :: naxes(3)
 real(kind=real64),intent(in) :: image(naxes(1),naxes(2),naxes(3))
 character(len=80), intent(in), optional :: hdr(:)
 integer,           intent(out) :: ierr

 ierr = 1

end subroutine write_fits_cube64

!-------------------------------------------------------------
! read fits file and convert to double precision
!-------------------------------------------------------------
subroutine read_fits_image64(filename,image,naxes,ierr,hdr)
 character(len=*), intent(in)   :: filename
 real(kind=real64), intent(out), allocatable :: image(:,:)
 character(len=80), intent(inout), allocatable, optional :: hdr(:)
 integer, intent(out) :: naxes(2),ierr

 ierr = 1
 naxes = 0

end subroutine read_fits_image64

!-------------------------------------------------------------
! read fits cube and convert to double precision
!-------------------------------------------------------------
subroutine read_fits_cube64(filename,image,naxes,ierr,hdr,velocity)
 character(len=*), intent(in)   :: filename
 real(kind=real64), intent(out), allocatable :: image(:,:,:)
 character(len=80), intent(inout), allocatable, optional :: hdr(:)
 real(kind=real64), intent(inout), allocatable, optional :: velocity(:)
 integer, intent(out) :: naxes(4),ierr

 ierr = 1
 naxes = 0

end subroutine read_fits_cube64

!--------------------------------------------------
! read all floating point variables from fits header
!--------------------------------------------------
subroutine get_floats_from_fits_header(hdr,tags,vals)
 character(len=80), intent(in) :: hdr(:)
 character(len=*),  intent(out) :: tags(:)
 real,              intent(out) :: vals(:)

 tags = ''
 vals = 0.

end subroutine get_floats_from_fits_header

!------------------------------------------------
! get tag:val pairs from fits header record
! will extract anything readable as a floating
! point number
!------------------------------------------------
subroutine get_fits_header_entry(record,key,rval,ierr)
 character(len=80), intent(in) :: record
 character(len=*),  intent(out) :: key
 real, intent(out) :: rval
 integer, intent(out) :: ierr

 ierr = 1
 rval = 0.

end subroutine get_fits_header_entry

!------------------------------------------------
! get tag:val pairs from fits header record
! returns the string value
!------------------------------------------------
subroutine get_fits_header_key_val(record,key,val,ierr)
 character(len=80), intent(in) :: record
 character(len=*),  intent(out) :: key
 character(len=*),  intent(out) :: val
 integer, intent(out) :: ierr

 ierr = 1
 key = ''
 val = ''

end subroutine get_fits_header_key_val

!------------------------------------------------
! search fits header to find a particular variable
! e.g. bmaj = get_from_header('BMAJ',hdr,ierr)
!------------------------------------------------
function get_from_header(key,hdr,ierr) result(val)
 character(len=*),  intent(in) :: key
 character(len=80), intent(in) :: hdr(:)
 integer,           intent(out) :: ierr
 real :: val

 ierr = 1
 val = 0.

end function get_from_header

!------------------------------------------------
! search fits header to find a particular string
! e.g. bmaj = get_from_header('BMAJ',hdr,ierr)
!------------------------------------------------
function get_from_header_s(key,hdr,ierr) result(val)
 character(len=*),  intent(in) :: key
 character(len=80), intent(in) :: hdr(:)
 integer,           intent(out) :: ierr
 character(len=1) :: val

 ierr = 1
 val = ''

end function get_from_header_s

!------------------------------------------------
! delete third dimension in the fits header
!------------------------------------------------
subroutine flatten_header(hdr)
 character(len=80), intent(inout) :: hdr(:)

end subroutine flatten_header

!------------------------------------------------
! get velocity grid from the fits file
!------------------------------------------------
subroutine get_velocity_from_fits_header(nv,vel,hdr,ierr)
 integer, intent(in)  :: nv
 real(kind=real32), intent(out) :: vel(nv)
 character(len=80), intent(in)  :: hdr(:)
 integer, intent(out) :: ierr

 ierr = 1
 vel = 0.

end subroutine get_velocity_from_fits_header

end module readwrite_fits
