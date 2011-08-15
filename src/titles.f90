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
!  Copyright (C) 2005-2009 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

module titles
 implicit none
 integer, parameter, private :: maxtitles = 50
 integer, parameter, private :: maxsteplegend = 100
 integer, parameter, public :: lensteplegend = 60
 integer, parameter, public :: lenpagetitles = 60
 character(len=lenpagetitles), dimension(maxtitles), public :: pagetitles
 character(len=lensteplegend), dimension(maxsteplegend), public :: steplegend
 public :: read_titles, read_steplegend

 private

contains
!
!--reads a list of titles (one per line), to be used to label each plot on page
!
subroutine read_titles(ntitles)
 use asciiutils, only:read_asciifile
 use filenames, only:fileprefix
 implicit none
 integer, intent(out) :: ntitles
 character(len=50) :: titlefile
 logical :: iexist

 titlefile = trim(fileprefix)//'.titles'
 !--also allow obsolete title filename
 inquire(file=titlefile,exist=iexist)
 if (.not.iexist) then
    inquire(file='titlelist',exist=iexist)
    if (iexist) titlefile='titlelist'
 endif
 ntitles = 0

 call read_asciifile(titlefile,ntitles,pagetitles)
 if (ntitles.gt.0) print "(a)",'read plot titles from file '//trim(titlefile)

 return
end subroutine read_titles
!
!--reads a list of labels (one per line) to be used in the timestep legend
! (ie. for multiple timesteps on same page)
!
subroutine read_steplegend(nsteptitles)
 use asciiutils, only:read_asciifile
 use filenames, only:fileprefix
 implicit none
 integer, intent(out) :: nsteptitles
 character(len=50) :: legendfile
 logical :: iexist

 legendfile = trim(fileprefix)//'.legend'
 !--also allow obsolete legend filename
 inquire(file=legendfile,exist=iexist)
 if (.not.iexist) then
    inquire(file='legend',exist=iexist)
    if (iexist) legendfile='legend'
 endif
 nsteptitles = 0

 call read_asciifile(legendfile,nsteptitles,steplegend)
 if (nsteptitles.gt.0) print "(a)"//'read legend text from file '''//trim(legendfile)//''''

 return
end subroutine read_steplegend

end module titles
