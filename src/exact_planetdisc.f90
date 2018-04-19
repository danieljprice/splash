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
!  Copyright (C) 2005-2017 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!---------------------------------------------------------------------------
! compute exact solution for various spiral structures:
!
! 1) planet-disc interaction in phi-r plane
!    (Ogilvie & Lubow (2002), MNRAS 330, 950)
! 2) planet-disc interaction from Rafikov (2002)
! 3) general parameterised spiral arms
!---------------------------------------------------------------------------
module planetdisc
  implicit none
  public :: exact_planetdisc
  integer, parameter :: maxspirals = 2
  integer, parameter :: maxcoeff = 5
  character(len=*), dimension(maxspirals), parameter, public :: labelspiral = &
    (/'Ogilvie-Lubow (2002) planet-disc interaction            ',&
      'Spiral arm fitting formula r(phi) = sum(a_i*phi^i,i=1,4)'/)

contains

subroutine exact_planetdisc(iplot,ispiral,time,HonR,rplanet,narms,params,rplot,yplot,ierr)
  use plotlib, only:plot_line
  implicit none
  integer, intent(in)  :: iplot,ispiral,narms
  integer, intent(out) :: ierr
  real,    intent(in)  :: time, HonR, rplanet, params(:,:)
  real, dimension(:),           intent(inout) :: rplot
  real, dimension(size(rplot)), intent(out)   :: yplot
  integer :: npts,iend,istart
  integer :: i,j,norbits,iarm
  real :: r,phase,dr,phi,rmin,rmax,phimin,phimax,dphi,coeff(maxcoeff)
  real, parameter :: pi = 4.*atan(1.)

  ierr = 0
  npts = size(rplot)
  norbits = int(time/(2.*pi))
  phase   = time - (2.*pi*norbits)
  select case(ispiral)
  case(2)
     print "(a,i2)",' Spiral arm fitting formula r = sum(a_i*phi^i,i=1,4) narms =',narms
  case default
     print "(a,f6.2,a,f8.1,a)",' Ogilvie-Lubow planet-disc interaction: H/R=',HonR,' at ',time/(2.*pi),' orbits)'
  end select

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
     rmax = min(max(maxval(rplot),abs(minval(rplot))),5.*rplanet)
     dr = (rmax - rmin)/npts
     do iarm=1,narms
        if (ispiral.eq.2) then
           phimin = params(1,iarm) + 90.  ! add 90 deg for East of North convention
           phimax = params(2,iarm) + 90.
           coeff(:) = params(3:,iarm)
           dphi = (phimax - phimin)/npts
           !print*,' GOT RMIN = ',rmin,phimin,phimax, 'COEFFS=',coeff
        endif
        ! outside planet
        do i=1,npts
           select case(ispiral)
           case(2)
              !
              ! Spiral arm fitting formula
              !
              phi = phimin + (i-1)*dphi
              r = 0. ! rmin
              do j=1,maxcoeff
                 r = r + coeff(j)*((phi-phimin)*pi/180.)**(j-1)
              enddo
           case default
              !
              ! Ogilvie & Lubow (2002)
              !
              r = rmin + (i-1)*dr
              if (r > rplanet) then
                 phi = phase - 2./(3.*HonR)*(sqrt(r**3) - 1.5*log(r) - 1.)
              else
                 phi = phase + 2./(3.*HonR)*(sqrt(r**3) - 1.5*log(r) - 1.)        
              endif
           end select
           rplot(i) = r*cos(phi*pi/180.)
           yplot(i) = r*sin(phi*pi/180.)
           !print*,'r, phi = ',r,phi,' : x, y = ',rplot(i),yplot(i)
        enddo
        call plot_line(npts,rplot,yplot)
     enddo
     ierr = 1 ! do not plot outside this routine
  end select

end subroutine exact_planetdisc

end module planetdisc
