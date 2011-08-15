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
!  Copyright (C) 2005-2011 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!  This module contributed by Andrew McLeod
!-----------------------------------------------------------------

module contours_module
 implicit none
 integer, parameter, private :: maxcontours = 50
 integer, parameter, public  :: lencontourtitles = 60
 real, dimension(maxcontours), public :: contours_list
 character(len=lencontourtitles), dimension(maxcontours), public :: contourtitles
 logical, public :: fixed_contours
 public :: read_contours

 private

contains
!
!--reads a list of contours (one per line), to be used on contour plots
!
subroutine read_contours(ncontours,ierr)
 use asciiutils, only:read_asciifile
 use filenames,  only:fileprefix
 implicit none
 integer, intent(out) :: ncontours, ierr
 character(len=50)    :: contourfile
 logical :: iexist

 contourfile = trim(fileprefix)//'.contours'
 ncontours = 0
 ierr = 0
 inquire(file=contourfile,exist=iexist)
 if (iexist) then
    call read_asciifile(contourfile,ncontours,contours_list,contourtitles)
 else
    contours_list(:) = 0.
    contourtitles(:) = ''
    ierr = 1
 endif
 if (ncontours.gt.0) then
    print "(1x,a)",'read contours and titles from file '//trim(contourfile)
 else
    ierr = -1
 endif

 return
end subroutine read_contours

end module contours_module
