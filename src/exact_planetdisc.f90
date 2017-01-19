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
!  Copyright (C) 2005-2015 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!---------------------------------------------------------------------------
! compute exact solution for planet-disc interaction in phi-r plane
! taken from Ogilvie & Lubow (2002), MNRAS 330, 950
!---------------------------------------------------------------------------
module planetdisc
  implicit none
  public :: exact_planetdisc

contains

subroutine exact_planetdisc(iplot,time,HonR,rplanet,rplot,yplot,ierr)
  use plotlib, only:plot_line
  implicit none
  integer, intent(in)  :: iplot
  integer, intent(out) :: ierr
  real, intent(in) :: time, HonR, rplanet
  real, dimension(:), intent(inout) :: rplot
  real, dimension(size(rplot)), intent(out) :: yplot
  integer :: npts,iend,istart
  integer :: i,norbits
  real :: r,phase,dr,phi,rmin,rmax
  real, parameter :: pi = 4.*atan(1.)

  ierr = 0
  npts = size(rplot)
  norbits = int(time/(2.*pi))
  phase   = time - (2.*pi*norbits)
  print "(a,f6.2,a,f8.1,a)",' Lubow-Ogilvie planet-disc interaction: H/R=',HonR,' at ',time/(2.*pi),' orbits)'

  select case(iplot)
  case(2)
     ! in phi-r plane
     istart = 1
     do i=1,npts
        r = rplot(i)
        if (r > rplanet) then
           yplot(i) = phase - 2./(3.*HonR)*(sqrt(r**3) - 1.5*log(r) - 1.)
        else
           yplot(i) = phase + 2./(3.*HonR)*(sqrt(r**3) - 1.5*log(r) - 1.)
        endif
        if (yplot(i) > pi) then
           phase = phase - 2.*pi
           if (i > 1) then
              iend = i
              call plot_line(iend-istart+1,rplot(istart:iend),yplot(istart:iend))
              istart = i+1
           endif
        endif
        if (yplot(i) <= -pi) then
           phase = phase + 2.*pi
           ! plot separate line segments every time we cross the phase boundary
           if (i > 1) then
              iend = i
              call plot_line(iend-istart+1,rplot(istart:iend),yplot(istart:iend))
              istart = i+1
           endif
        endif
     enddo
     ierr = 1 ! do not plot outside this routine
     return
  case default ! in x-y plane
     ! define npts outside planet orbit
     rmin = 1.e-3
     rmax = max(maxval(rplot),abs(minval(rplot)))
     dr = (rmax - rmin)/npts
     ! outside planet
     do i=1,npts
        r = rmin + (i-1)*dr
        if (r > rplanet) then
           phi = phase - 2./(3.*HonR)*(sqrt(r**3) - 1.5*log(r) - 1.)
        else
           phi = phase + 2./(3.*HonR)*(sqrt(r**3) - 1.5*log(r) - 1.)        
        endif
        rplot(i) = r*cos(phi)
        yplot(i) = r*sin(phi)
     enddo
     call plot_line(npts,rplot,yplot)
  end select

end subroutine exact_planetdisc

end module planetdisc
