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
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!------------------------------------------------
! reads an exact solution from a file
!
! file should contain two or more columns containing
! x axis data and y axis data
!
! this is plotted as a line on the chosen graph
!------------------------------------------------
module exactfromfile
 implicit none

contains

subroutine exact_fromfile(filename,xexact,yexact,ixcolfile,iycolfile,iexactpts,ierr)
  use asciiutils, only:get_ncolumns
  implicit none
  character(len=*), intent(in) :: filename
  real, intent(out), dimension(:) :: xexact, yexact
  integer, intent(in)  :: ixcolfile,iycolfile
  integer, intent(out) :: iexactpts, ierr
  integer :: i,j,ncolumns,nheaderlines
  integer, parameter :: lu = 33
  character(len=10) :: str
  real :: dum

  ierr = 0
  open(unit=lu,file=filename,iostat=ierr,status='old',form='formatted')
  if (ierr /= 0) then
     ierr = 1
     print*,'error opening ',filename
     return
  endif

  !--query number of header lines
  call get_ncolumns(lu,ncolumns,nheaderlines)

  !--skip header lines
  do i=1,nheaderlines
     read(lu,*)
  enddo

  !--read data from file
  do i=1,size(xexact)
     if (ixcolfile.gt.iycolfile) then
        read(lu,*,end=10,err=20) (dum,j=1,iycolfile-1),yexact(i),(dum,j=iycolfile+1,ixcolfile-1),xexact(i)     
     elseif (ixcolfile.eq.iycolfile) then
        read(lu,*,end=10,err=20) (dum,j=1,ixcolfile-1),xexact(i)
        yexact(i) = xexact(i)
     else
        read(lu,*,end=10,err=20) (dum,j=1,ixcolfile-1),xexact(i),(dum,j=ixcolfile+1,iycolfile-1),yexact(i)
     endif
  enddo
  print*,'WARNING: reached array limits in ',trim(filename),': partial solution read'
  ierr = -1
  close(lu)
  return
10 continue
  iexactpts = i-1
  write(str,"(i10)") iexactpts
  print "(a)",' finished reading '//trim(filename)//': '//trim(adjustl(str))//' read'
  close(lu)
  return
20 print*,'error reading ',trim(filename),': partial solution read'
  iexactpts = i - 1
  ierr = -2
  close(lu)
  return

end subroutine exact_fromfile

end module exactfromfile
