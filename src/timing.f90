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
!  Copyright (C) 2005-2013 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!----------------------------------------------------------------
!
! This module contains utilities for code timings
!
!----------------------------------------------------------------
module timing
 implicit none
 integer, private :: istarttime(6)
 real,    private :: starttime

 data starttime/-1./

 public :: wall_time,print_time

 private

contains
!--------------------------------------------------------------------
!+
!  sets initial time
!+
!--------------------------------------------------------------------
 subroutine initialise_timing
  implicit none
  integer :: iday,imonth,iyear,ihour,imin,isec,imsec,ivalues(8)
  character(len=8)  :: date
  character(len=5)  :: zone
  character(len=10) :: time

  call date_and_time(date,time,zone,ivalues)
  iyear  = ivalues(1)
  imonth = ivalues(2)
  iday   = ivalues(3)
  ihour  = ivalues(5)
  imin   = ivalues(6)
  isec   = ivalues(7)
  imsec  = ivalues(8)
  istarttime(1) = iyear
  istarttime(2) = imonth
  istarttime(3) = iday
  istarttime(4) = ihour
  istarttime(5) = imin
  istarttime(6) = isec
  !istarttime(7) = imsec
  starttime = iday*86400. + ihour*3600. + imin*60. + isec + imsec*0.001

  return
 end subroutine initialise_timing

!--------------------------------------------------------------------
!+
!  Get time used since begining
!+
!--------------------------------------------------------------------
 subroutine wall_time(t)
  implicit none
  real, intent(out) :: t
  integer :: i,iday,imonth,ihour,imin,isec,imsec,ivalues(8)
  character(len=8)  :: date
  character(len=5)  :: zone
  character(len=10) :: time

  !--do self-initialisation the first time it is called
  if (starttime.lt.0.) call initialise_timing

  call date_and_time(date,time,zone,ivalues)
  iday   = ivalues(3)
  ihour  = ivalues(5)
  imin   = ivalues(6)
  isec   = ivalues(7)
  imsec  = ivalues(8)

  if (ivalues(2).lt.istarttime(2)) then
     ivalues(2) = ivalues(2) + 12
  endif
  do i = istarttime(2), ivalues(2) - 1
     imonth = mod(i,12)
     if (imonth.eq.4 .or. imonth.eq.6 .or. imonth.eq.9 .or. imonth.eq.11) then
        iday = iday + 30
     elseif (imonth.eq.2) then
        iday = iday + 28
     else
        iday = iday + 31
     endif
  end do
  t = iday*86400. + ihour*3600. + imin*60. + isec + imsec*0.001 - starttime

  return
 end subroutine wall_time

!--------------------------------------------------------------------
!
!  print a time, nicely formatted into hours, mins, seconds
!
!--------------------------------------------------------------------
 subroutine print_time(time,string,iunit)
  implicit none
  real, intent(in) :: time
  character(len=*), intent(in), optional :: string
  integer, intent(in), optional :: iunit
  character(len=64) :: newstring
  integer :: nhr,nmin,lunit
  real :: trem

  trem = time
  nhr = int(trem/3600.)
  if (nhr.gt.0) trem = trem - nhr*3600.

  nmin = int(trem/60.)
  if (nmin.gt.0) trem = trem - nmin*60.

  if (present(string)) then
     newstring = trim(string(1:min(len(newstring),len_trim(string))))
  else
     newstring = 'completed in'
  endif

  if (present(iunit)) then
     lunit = iunit
  else
     lunit = 6
  endif

  if (nhr.gt.0) then
     write(lunit,"(1x,a,1x,i3,a,i2,a,f6.2,a)") &
          trim(newstring),nhr,' hr, ',nmin,' min, ',trem,' s'
  elseif (nmin.gt.0) then
     write(lunit,"(1x,a,1x,i2,a,f6.2,a)") &
          trim(newstring),nmin,' min, ',trem,' s'
  else
     write(lunit,"(1x,a,1x,f6.2,a)") trim(newstring),trem,' s'
  endif

  return
 end subroutine print_time

end module timing
