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

! ----------------------------------------------------------------------
! Exact solution for Bondi flow
! ----------------------------------------------------------------------
module bondi
 implicit none
 public :: exact_bondi

contains

subroutine exact_bondi(iplot,time,gamma,r0,m,relativistic,xpts,ypts,ierr)
 integer, intent(in)  :: iplot
 integer, intent(out) :: ierr
 real, intent(in)     :: time,gamma,r0,m
 logical, intent(in)  :: relativistic
 real, dimension(:), intent(in) :: xpts
 real, dimension(size(xpts)), intent(out) :: ypts
 real, parameter :: pi = 3.1415926536
 integer :: i,npts
 real :: r,vr

 npts = size(xpts)

 print "(4(a,g8.2))",' Plotting exact Bondi solution at t = ',time

 !
 !--determine which parameter to plot
 !
 do i=1,npts
    r = xpts(i)
    vr = -(1. - 2.*m/r)/sqrt(1.-2.*m/r0)*sqrt(2.*m*(1./r - 1./r0))

    select case(iplot)
    case(1)
       ypts(i) = vr
    case default
       ypts(i) = 0.
    end select
 enddo
 ierr = 0

 return
end subroutine exact_bondi

end module bondi
