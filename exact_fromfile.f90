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

!------------------------------------------------
! reads an exact solution from a file
!
! file should contain two columns containing
! x axis data and y axis data
!
! this is plotted as a line on the chosen graph
!------------------------------------------------
module exactfromfile
 implicit none

contains

subroutine exact_fromfile(filename,xexact,yexact,iexactpts,ierr)
  implicit none
  character(len=*), intent(in) :: filename
  real, intent(out), dimension(:) :: xexact, yexact
  integer, intent(out) :: iexactpts, ierr
  integer :: i

  ierr = 0
  open(unit=33,file=filename,iostat=ierr,status='old',form='formatted')
  if (ierr /= 0) then
     ierr = 1
     print*,'error opening ',filename
     return
  endif
  do i=1,size(xexact)
     read(33,*,end=10,err=20) xexact(i),yexact(i)
  enddo
  print*,'WARNING: reached array limits in ',trim(filename),': partial solution read'
  ierr = -1
  close(33)
  return
10 continue
  iexactpts = i-1
  print*,'finished reading ',trim(filename),' : ',iexactpts,' read'
  close(33)
  return
20 print*,'error reading ',trim(filename),': partial solution read'
  iexactpts = i - 1
  ierr = -2
  close(33)
  return

end subroutine exact_fromfile

end module exactfromfile
